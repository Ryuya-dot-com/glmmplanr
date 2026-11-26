# Practical pilot: fixed short forms, recalibration, and stimulus sampling.
# All items are synthetic; selected IDs are fixed before simulated outcomes.
library(glmmplanr)
args <- commandArgs(TRUE)
out <- if (length(args)) args[1L] else "validation/item_pool_planning"
n_rep <- if (length(args) > 1L) as.integer(args[2L]) else 128L
section <- if (length(args) > 2L) match.arg(args[3L], c("all", "calibration")) else "all"
stopifnot(n_rep > 0L)
if (dir.exists(out)) stop("Use a new output directory; existing evidence is not overwritten")
dir.create(out, recursive = TRUE)
capture.output(sessionInfo(), file = file.path(out, "environment.txt"))
files <- c("src/glmm_laplace.cpp", "R/rasch_precision.R", "R/engines.R",
           "validation/item_pool_planning.R", list.files(system.file("libs", package = "glmmplanr"),
             pattern = "\\.(so|dll|dylib)$", recursive = TRUE, full.names = TRUE))
write.csv(data.frame(file = files, md5 = unname(tools::md5sum(files))),
          file.path(out, "source.csv"), row.names = FALSE)
seed <- 92611001L
pool <- data.frame(item_id = sprintf("word_%03d", 1:150),
  domain = rep(c("A", "B", "C"), each = 50L),
  difficulty = rep(seq(-2.5, 2.5, length.out = 50L), 3L))
select_items <- function(k, rule) {
  unlist(lapply(split(seq_len(nrow(pool)), pool$domain), function(id) {
    if (rule == "central") id[order(abs(pool$difficulty[id]), id)[seq_len(k / 3L)]]
    else id[round(seq(1, length(id), length.out = k / 3L))]
  }), use.names = FALSE)
}
forms <- list(central30 = select_items(30L, "central"),
              spread30 = select_items(30L, "spread"),
              spread60 = select_items(60L, "spread"), full150 = 1:150)
stopifnot(all(vapply(forms, function(id) !anyDuplicated(id) &&
  length(unique(table(pool$domain[id]))) == 1L, logical(1L))))
write.csv(pool, file.path(out, "pool.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(names(forms), function(nm)
  data.frame(form = nm, pool[forms[[nm]], ]))), file.path(out, "forms.csv"), row.names = FALSE)
writeLines(c("Synthetic Rasch: a=1, no guessing, independent items given theta.",
  "Known-bank EAP scoring prior N(0,1); no item estimation uncertainty.",
  "Recalibration: theta mean zero, SD 0.4/0.8, N=120/240/480; MML.",
  "Joint error limits: bank RMSE <= .35 and max abs <= .60 logits.",
  "Assurance .80, pointwise Wilson 95% intervals, invalid = miss.",
  "GLMM: mean accuracy .50/.58; N=60/120, J=30/60/150, 30 per person.",
  "Person/item intercept SD .8/.6; slope SD 0/0 or .25/.35, rho=0.",
  "Marginal 8pp target; 95% CI width <= 20pp; native joint_delta.",
  sprintf("repetitions=%d master_seed=%d", n_rep, seed)), file.path(out, "protocol.txt"))

# Numerical expectation over the complete score distribution, not simulated
# people. The Poisson-binomial recursion retains all scores, including 0/J.
known_bank <- function(id, step = .025) {
  theta <- seq(-8, 8, by = step)
  b <- pool$difficulty[id]
  p <- plogis(outer(theta, b, "-"))
  score_prob <- matrix(1, length(theta), 1L)
  for (j in seq_along(id)) {
    score_prob <- cbind(score_prob * (1 - p[, j]), 0) + cbind(0, score_prob * p[, j])
  }
  prior <- dnorm(theta) * step
  mass <- colSums(score_prob * prior)
  eap <- drop(crossprod(theta * prior, score_prob)) / mass
  second <- drop(crossprod(theta^2 * prior, score_prob)) / mass
  bias <- drop(score_prob %*% eap) - theta
  mse <- rowSums(score_prob * sweep(matrix(eap, length(theta), length(eap), byrow = TRUE),
                                  1L, theta, "-")^2)
  stopifnot(max(abs(rowSums(score_prob) - 1)) < 1e-12, abs(sum(prior) - 1) < 1e-10)
  data.frame(theta, bias, rmse = sqrt(mse),
    mean_posterior_sd = drop(score_prob %*% sqrt(pmax(0, second - eap^2))),
    information = rowSums(p * (1 - p)))
}
curves <- do.call(rbind, lapply(names(forms), function(nm)
  data.frame(form = nm, n_item = length(forms[[nm]]), known_bank(forms[[nm]]))))
write.csv(curves, file.path(out, "known_bank_curves.csv"), row.names = FALSE)
fine <- known_bank(forms$full150, step = .0125)
coarse <- subset(curves, form == "full150")
matched <- fine[seq(1, nrow(fine), by = 2L), ]
quad_error <- max(abs(as.matrix(coarse[c("bias", "rmse")]) - as.matrix(matched[c("bias", "rmse")])))
stopifnot(quad_error < 1e-6)
writeLines(sprintf("Known-bank grid refinement maximum error: %.12g", quad_error),
           file.path(out, "integration_check.txt"))

layouts <- rbind(data.frame(form = names(forms), n_item = lengths(forms),
                            items_per_subject = lengths(forms)),
                 data.frame(form = "rotating150", n_item = 150L, items_per_subject = 30L))
cases <- merge(layouts, expand.grid(person_sd = c(.4, .8), n_subj = c(120L, 240L, 480L)))
cases$case <- seq_len(nrow(cases))
cases$responses <- with(cases, n_subj * items_per_subject)
cases$average_item_exposure <- with(cases, responses / n_item)
cases$n_rep <- n_rep
write.csv(cases, file.path(out, "calibration_scenarios.csv"), row.names = FALSE)
calibration <- lapply(seq_len(nrow(cases)), function(k) {
  x <- cases[k, ]
  id <- if (x$form == "rotating150") 1:150 else forms[[x$form]]
  b <- setNames(pool$difficulty[id], pool$item_id[id])
  elapsed <- system.time(plan <- plan_rasch_fixed_bank(b, x$n_subj,
    max_error = c(difficulty = .35, difficulty_max_abs_error = .60),
    person_sd = x$person_sd, items_per_subject = x$items_per_subject,
    method = "fixed_item", engine = "mml", nsim = n_rep,
    seed = seed, item_seed = 37L, workers = 1L))[["elapsed"]]
  saveRDS(plan, file.path(out, sprintf("calibration-%02d.rds", k)))
  cat(sprintf("Calibration %d/%d: %.2fs\n", k, nrow(cases), elapsed))
  data.frame(x[c("case", "form", "person_sd", "responses", "average_item_exposure")],
             elapsed_s = elapsed, as.data.frame(plan), row.names = NULL)
})
calibration <- do.call(rbind, calibration)
write.csv(calibration, file.path(out, "calibration_summary.csv"), row.names = FALSE)
if (section == "calibration") quit(save = "no")

# New items are sampled each replication in this experiment. These are
# population design scenarios, not estimates for the particular word IDs above.
experiment <- expand.grid(n_subj = c(60L, 120L), n_item = c(30L, 60L, 150L),
                          slopes = c("zero", "heterogeneous"), stringsAsFactors = FALSE)
experiment$items_per_subject <- 30L
experiment$responses <- experiment$n_subj * 30L
experiment$item_per_condition <- experiment$responses / (2L * experiment$n_item)
write.csv(experiment, file.path(out, "experiment_scenarios.csv"), row.names = FALSE)
checks <- lapply(seq_len(nrow(experiment)), function(k) {
  x <- experiment[k, ]
  des <- design_crossed(list(condition = c("control", "treatment")),
    allocation = "latin_square", items_per_subject = 30L, item_seed = 37L)
  varying <- x$slopes == "heterogeneous"
  dgp <- calibrate_dgp(des, c(control = .50, treatment = .58),
    subj_sd = c(.8, if (varying) .25 else 0),
    item_sd = c(.6, if (varying) .35 else 0), corr_structure = "zero")
  target <- specify_target(dgp, "condition", scale = "pp", criterion = "both",
    ci_width = 20, engine = "bayes_laplace", analysis_re = "match_dgp",
    prob_se_method = "joint_delta")
  elapsed <- system.time(check <- validate_design_point(target, x$n_subj, x$n_item,
    nsim = n_rep, seed = seed + 10000L))[["elapsed"]]
  saveRDS(check, file.path(out, sprintf("experiment-%02d.rds", k)))
  cat(sprintf("Experiment %d/%d: %.2fs\n", k, nrow(experiment), elapsed))
  data.frame(x, elapsed_s = elapsed, as.data.frame(check), row.names = NULL)
})
write.csv(do.call(rbind, checks), file.path(out, "experiment_summary.csv"), row.names = FALSE)
