# Public irtsim pipeline and native MML on exactly the same response matrices.
source("validation/planning_comparison_helpers.R")
args <- commandArgs(TRUE)
stage <- if (length(args)) match.arg(args[1L], c("smoke", "screen", "confirm")) else "screen"
out <- if (length(args) > 1L) args[2L] else file.path("validation", paste0("item_pool_irtsim_", stage))
workers <- if (length(args) > 2L) as.integer(args[3L]) else 2L
if (.Platform$OS.type == "windows") workers <- 1L
n_rep <- switch(stage, smoke = 3L, screen = 128L, confirm = 512L)
cases <- rasch_cases
if (stage == "smoke") cases <- cases[c(1L, 5L, 9L, 50L), ]
if (stage == "confirm") {
  selection <- read.csv("validation/item_pool_irtsim_screen/selected.csv")
  search <- read.csv("validation/item_pool_search/search.csv")
  selected_cases <- unique(c(selection$case, search$case[search$task == "rasch"]))
  cases <- cases[cases$case %in% selected_cases[!is.na(selected_cases)], ]
}
checkpoint_dir(out, "validation/compare_item_pool_irtsim.R",
  list(stage = stage, n_rep = n_rep, cases = cases, target = c(rmse = .35, maximum = .60, assurance = .80)),
  c("glmmplanr", "irtsim", "mirt"))
write.csv(cases, file.path(out, "scenarios.csv"), row.names = FALSE)
run_case <- function(k) {
  case <- cases[k, ]; id <- case$case
  path <- file.path(out, sprintf("case-%02d.rds", id))
  if (file.exists(path)) return(readRDS(path)$summary)
  b <- bank_for(case)
  schedule <- generate_design_schedule(design_crossed(list(), items_per_subject = case$m, item_seed = 37L),
                                       case$n_subj, length(b))
  mask <- matrix(0L, case$n_subj, length(b))
  mask[cbind(schedule$trials$subj_id, schedule$trials$item_id)] <- 1L
  stopifnot(all(rowSums(mask) == case$m), all(colSums(mask) > 0L))
  # Independent ability streams avoid reusing the response RNG when irtsim's
  # sequential generator resets its cell seed inside generate_data().
  ability_seed <- 93710001L + if (stage == "confirm") 10000L else 0L
  ability_streams <- glmmplanr:::make_lecuyer_streams(ability_seed, n_rep)
  rep_id <- 0L
  theta_dist <- function(n) {
    rep_id <<- rep_id + 1L
    glmmplanr:::with_rng_stream(ability_streams[[rep_id]], rnorm(n, sd = case$person_sd))
  }
  design <- irtsim::irt_design("1PL", length(b), list(b = unname(b)), theta_dist = theta_dist)
  study <- if (case$m < length(b)) irtsim::irt_study(design, case$n_subj,
    missing = "booklet", test_design = list(booklet_matrix = mask)) else irtsim::irt_study(design, case$n_subj)
  seed <- 9370L + id + if (stage == "confirm") 1000L else 0L
  irt_time <- system.time(irt <- irtsim::irt_simulate(study, iterations = n_rep, seed = seed,
    progress = FALSE, parallel = FALSE, se = FALSE, compute_theta = FALSE))[["elapsed"]]
  stopifnot(rep_id == n_rep, nrow(irt$item_results) == n_rep * length(b))
  item_est <- matrix(NA_real_, n_rep, length(b))
  rows <- vector("list", n_rep)
  native_start <- proc.time()[["elapsed"]]
  for (r in seq_len(n_rep)) {
    # This is the installed irtsim 0.2.0 cell-seed expression, with one N per study.
    cell_seed <- as.integer(((seed * 31L + r) * 131L + 1L) %% .Machine$integer.max)
    theta <- glmmplanr:::with_rng_stream(ability_streams[[r]], rnorm(case$n_subj, sd = case$person_sd))
    y <- irtsim:::generate_data(design, case$n_subj, seed = cell_seed, theta = theta)
    y <- irtsim:::apply_missing(y, study, seed = cell_seed + .Machine$integer.max %/% 2L, theta = theta)
    stopifnot(identical(unname(!is.na(y)), unname(mask == 1L)))
    ir <- subset(irt$item_results, iteration == r)
    stopifnot(identical(ir$item, seq_along(b)), all(ir$true_value == unname(b)))
    # Refitting two reconstructed matrices checks exact response replay against
    # the public pipeline, beyond merely giving two generators the same seed.
    if (r <= 2L) {
      ref <- irtsim:::fit_model(y, "1PL", se = FALSE)
      stopifnot(identical(isTRUE(ref$converged), all(ir$converged)))
      if (isTRUE(ref$converged)) {
        reproduced <- -mirt::coef(ref$model, simplify = TRUE)$items[, "d"]
        stopifnot(max(abs(reproduced - ir$estimate)) < 1e-10)
      }
    }
    tr <- schedule$trials
    tr$y <- y[cbind(tr$subj_id, tr$item_id)]
    sim <- list(trials = tr, theta = theta, n_subj = case$n_subj, n_item = length(b))
    fit <- glmmplanr:::eval_fixed_bank_replication(sim, b, case$person_sd, "fixed_item", "mml")
    item_est[r, ] <- fit$b_est
    bv <- is.na(fit$metrics$fit_error) && all(is.finite(fit$b_est))
    iv <- all(ir$converged) && all(is.finite(ir$estimate))
    err <- list(native = fit$b_est - b, irtsim = ir$estimate - b)
    for (engine in names(err)) {
      valid <- if (engine == "native") bv else iv
      rmse <- if (valid) sqrt(mean(err[[engine]]^2)) else NA_real_
      maximum <- if (valid) max(abs(err[[engine]])) else NA_real_
      row <- data.frame(case = id, replicate = r, engine, valid, rmse, maximum,
        success = valid && rmse <= .35 && maximum <= .60,
        data_hash = digest::digest(y, algo = "xxhash64"),
        b_max_difference = if (bv && iv) max(abs(fit$b_est - ir$estimate)) else NA_real_)
      rows[[r]] <- rbind(rows[[r]], row)
    }
  }
  native_time <- proc.time()[["elapsed"]] - native_start
  rows <- do.call(rbind, rows)
  summary <- do.call(rbind, lapply(split(rows, rows$engine), function(z)
    data.frame(case, engine = z$engine[1], mc_summary(z$success, z$valid),
      elapsed_s = if (z$engine[1] == "native") native_time else irt_time)))
  # Native replay time includes two extra mirt fits used to audit pairing;
  # do not interpret it as a pure engine-speed benchmark.
  atomic_rds(list(summary = summary, trials = rows, irtsim_items = irt$item_results,
    native_estimates = item_est, seed = seed, ability_seed = ability_seed), path)
  cat(sprintf("%s case %d: irtsim %.2fs; native replay/audit %.2fs\n", stage, id, irt_time, native_time))
  summary
}
result <- parallel::mclapply(seq_len(nrow(cases)), run_case, mc.cores = workers, mc.set.seed = FALSE)
stopifnot(all(vapply(result, is.data.frame, logical(1))))
result <- do.call(rbind, result)
write.csv(result, file.path(out, "summary.csv"), row.names = FALSE)
rows <- do.call(rbind, lapply(cases$case, function(id) readRDS(file.path(out, sprintf("case-%02d.rds", id)))$trials))
con <- gzfile(file.path(out, "fits.csv.gz"), "wt"); write.csv(rows, con, row.names = FALSE); close(con)
selected <- do.call(rbind, lapply(split(result, interaction(result$family, result$engine)), function(z) {
  # Point-target proposals use the same .80 rule as the matched search study;
  # only the independent confirmation can supply a Wilson-qualified decision.
  pass <- z[z$attainment >= .80, ]
  if (nrow(pass)) pass[which.min(pass$n_subj), ] else {
    z <- z[which.max(z$n_subj), ]; z$case <- NA_integer_; z$n_subj <- NA_integer_; z
  }
}))
write.csv(selected, file.path(out, "selected.csv"), row.names = FALSE)
