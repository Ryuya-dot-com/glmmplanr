# Run with: sh validation/run-isolated.sh package_calibration
library(glmmplanr)
library(lme4)
stopifnot(requireNamespace("simr"), requireNamespace("PowRPriori"))
RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
set.seed(8210000)
n <- 192L
workers <- if (.Platform$OS.type == "windows") 1L else 6L
layouts <- c("latin_square", "between_subj")
effects <- c("null", "planned")
inputs <- c("direct", "calibrated")
levels_group <- c("A", "B")
model <- y ~ group + (1 | subj) + (1 | item)
test <- simr::fixed("groupB", "z")
jobs <- expand.grid(layout = layouts, effect = effects, input = inputs, stringsAsFactors = FALSE)
marginal <- function(eta) {
  vapply(eta, function(b) integrate(function(z) plogis(b + z) * dnorm(z),
                                     -Inf, Inf, rel.tol = 1e-11)$value, numeric(1))
}
prepare <- function(j) {
  job <- jobs[j, ]
  des <- design_crossed(list(group = levels_group), allocation = job$layout, contrast = "treatment")
  target <- c(.72, if (job$effect == "null") .72 else .80)
  make_dgp <- function(p) calibrate_dgp(des, cell_probs = p, chance = 0,
    subj_sd = .8, item_sd = .6, slope_sd_ratio = 0, corr_structure = 0)
  dgp <- make_dgp(target)
  diag <- check_dgp(dgp, plot = FALSE)$cell_diagnostics
  stopifnot(identical(diag$cell, levels_group))
  p <- if (job$input == "direct") target else diag$conditional_prob
  cell <- data.frame(group = factor(levels_group, levels = levels_group), y = p)
  b <- PowRPriori::fixed_effects_from_average_outcome(y ~ group, cell, center = FALSE, family = "binomial")
  X <- model.matrix(~ group, cell)
  stopifnot(identical(names(b), colnames(X)))
  implied <- unname(marginal(drop(X %*% unlist(b))))
  if (job$input == "direct") dgp <- make_dgp(implied)
  if (job$input == "calibrated") stopifnot(max(abs(implied - target)) < 1e-5)
  stopifnot(max(abs(dgp$mu_cells - drop(X %*% unlist(b)))) < 1e-5)
  sp <- specify_target(dgp, "group", engine = "bayes_laplace", scale = "logit",
                       analysis_re = "intercept_only", on_singular = "flag_only")
  schedule <- generate_design_schedule(des, 64, 24)
  tr <- schedule$trials
  ext <- tr[c("subj", "item", "group")]
  ext$group <- factor(as.character(tr$group), levels = levels_group)
  a <- xtabs(~ subj + group, ext)
  stopifnot(all(table(ext$subj, ext$item) == 1), all(colSums(a) == 768),
            if (job$layout == "latin_square") all(a == 12) else all(a %in% c(0,24)),
            all(xtabs(~ item + group, ext) == 32))
  template <- simr::makeGlmer(y ~ group + (1 | subj) + (1 | item), family = binomial(),
    fixef = unlist(b), VarCorr = list(subj = .8^2, item = .6^2), data = ext)
  stopifnot(max(abs(fixef(template) - unlist(b))) < 1e-10,
            abs(VarCorr(template)$subj[1,1] - .8^2) < 1e-10,
            abs(VarCorr(template)$item[1,1] - .6^2) < 1e-10)
  list(sp = sp, template = template, schedule = schedule, b = b,
       calibration = data.frame(job, group = levels_group, target, supplied = p, implied,
                                intercept = b[[1]], slope = b[[2]], row.names = NULL))
}
settings <- lapply(seq_len(nrow(jobs)), prepare)
calibration <- do.call(rbind, lapply(settings, "[[", "calibration"))
write.csv(calibration, "validation/package_calibration_cells.csv", row.names = FALSE)

fit_reference <- function(y, template) {
  tryCatch({
    f <- suppressMessages(simr::doFit(y, template, subset = seq_along(y),
                                      control = glmerControl(optimizer = "bobyqa")))
    conv <- f@optinfo$conv
    msg <- conv$lme4$messages
    valid <- all(conv$opt == 0L) && all(grepl("boundary (singular) fit", msg, fixed = TRUE))
    p <- as.numeric(simr::doTest(f, test))
    data.frame(simr_p = p, simr_valid = valid && is.finite(p), simr_singular = isSingular(f),
               simr_diagnostic = paste(msg, collapse = " | "))
  }, error = function(e) data.frame(simr_p = NA_real_, simr_valid = FALSE, simr_singular = NA,
                                     simr_diagnostic = conditionMessage(e)))
}
one_trial <- function(index) {
  j <- 1L + (index - 1L) %% nrow(jobs)
  trial <- 1L + (index - 1L) %/% nrow(jobs)
  job <- jobs[j, ]
  s <- settings[[j]]
  seed <- 8210000L + 10000L * match(job$layout, layouts) + 1000L * match(job$effect, effects) + trial
  set.seed(seed)
  tr <- s$schedule$trials
  tr$y <- simr::doSim(s$template)
  stopifnot(is.numeric(tr$y), length(tr$y) == nrow(tr), all(tr$y %in% 0:1))
  native <- glmmplanr:::fit_single_dataset(tr, s$sp)
  valid <- isTRUE(native$converged) && is.finite(native$est) && is.finite(native$se) && native$se > 0
  p <- if (valid) 2 * pnorm(-abs(native$est/native$se)) else NA_real_
  ref <- fit_reference(tr$y, s$template)
  set.seed(seed + 100000L)
  independent <- simulate_dgp_data(s$sp$dgp, s$schedule)
  cbind(data.frame(job, trial, seed, simr_A = mean(tr$y[tr$group == "A"]),
    simr_B = mean(tr$y[tr$group == "B"]), glmmplanr_A = mean(independent$y[independent$group == "A"]),
    glmmplanr_B = mean(independent$y[independent$group == "B"]),
    native_p = p, native_valid = valid, native_singular = native$singular,
    native_diagnostic = native$status), ref)
}
collect <- function(ids, fun) {
  x <- parallel::mclapply(ids, fun, mc.cores = workers, mc.set.seed = FALSE)
  failed <- which(!vapply(x, is.data.frame, logical(1)))
  if (length(failed)) stop("Incomplete job ", ids[failed[1]], ": ", as.character(x[[failed[1]]]))
  do.call(rbind, x)
}
raw <- collect(seq_len(nrow(jobs) * n), one_trial)
stopifnot(!anyDuplicated(raw[c(names(jobs), "trial")]))
write.csv(raw, "validation/package_calibration_fits.csv", row.names = FALSE)
cat("Completed paired fits:", nrow(raw), "data sets\n")

# PowRPriori's within option adds both conditions to each person-item pair.
# Its between-person public workflow has the same assignment up to person ID.
pow_jobs <- expand.grid(effect = effects, input = inputs, batch = 1:12, stringsAsFactors = FALSE)
pow_design <- PowRPriori::define_design(sample_size = list(subj = 64L, item = 24L),
  between = list(subj = list(group = factor(levels_group))))
pow_batch <- function(i) {
  job <- pow_jobs[i, ]
  s <- settings[[which(jobs$layout == "between_subj" & jobs$effect == job$effect & jobs$input == job$input)]]
  seed <- 8450000L + 1000L * match(job$effect, effects) + 100L * match(job$input, inputs) + job$batch
  set.seed(seed)
  fit <- tryCatch(suppressMessages(PowRPriori::power_sim(model, pow_design, test_parameter = "groupB",
    family = "binomial", fixed_effects = s$b, center = FALSE, adjust_p_value = FALSE,
    random_effects = list(subj = list(`(Intercept)` = .8), item = list(`(Intercept)` = .6)),
    along = "subj", n_start = 64L, n_increment = 2L, max_simulation_steps = 1L,
    n_sims = 16L, parallel_plan = "sequential")), error = identity)
  if (inherits(fit, "error")) return(data.frame(job, seed, attempted = 16L, reported_hits = NA_integer_,
    singular = NA_integer_, convergence = NA_integer_, identifiability = NA_integer_,
    other_warnings = NA_integer_, error = conditionMessage(fit)))
  t <- fit$power_table
  d <- fit$sample_data
  # sample_data is an extra illustration, possibly at the next size; power_table records the fitted N.
  stopifnot(nrow(t) == 1, t$subj == 64L, all(xtabs(~ subj + item, d) == 1),
            length(unique(d$item)) == 24, length(unique(table(d$group))) == 1L,
            all(xtabs(~ subj + group, d) %in% c(0,24)),
            is.finite(t$power_groupB), abs(16*t$power_groupB - round(16*t$power_groupB)) < 1e-8)
  data.frame(job, seed, attempted = 16L, reported_hits = round(16*t$power_groupB),
    singular = t$n_singular, convergence = t$n_convergence, identifiability = t$n_identifiable,
    other_warnings = t$n_other_warnings, error = "")
}
batches <- collect(seq_len(nrow(pow_jobs)), pow_batch)
write.csv(batches, "validation/package_calibration_pow.csv", row.names = FALSE)
summarize <- function(key, engine, attempted, valid, hits) {
  ci <- if (is.na(hits)) c(NA, NA) else binom.test(hits, attempted)$conf.int
  data.frame(key, engine, attempted, valid, detected = hits, rate = hits/attempted, lo = ci[1], hi = ci[2])
}
rates <- list()
for (j in seq_len(nrow(jobs))) {
  job <- jobs[j, ]
  r <- raw[raw$layout == job$layout & raw$effect == job$effect & raw$input == job$input, ]
  stopifnot(nrow(r) == n)
  for (engine in c("native","simr")) {
    ok <- r[[paste0(engine, "_valid")]]
    hits <- sum(ok & r[[paste0(engine, "_p")]] < .05)
    rates[[length(rates)+1L]] <- summarize(job, engine, n, sum(ok), hits)
  }
  if (job$layout == "between_subj") {
    b <- batches[batches$effect == job$effect & batches$input == job$input, ]
    clean <- all(!is.na(b$reported_hits) & b$convergence == 0 & b$identifiability == 0 & b$other_warnings == 0)
    rates[[length(rates)+1L]] <- summarize(job, "PowRPriori", sum(b$attempted),
      if (clean) sum(b$attempted) else NA_integer_, if (clean) sum(b$reported_hits) else NA_integer_)
  }
}
summary <- do.call(rbind, rates)
write.csv(summary, "validation/package_calibration_rates.csv", row.names = FALSE)
print(summary, row.names = FALSE)
