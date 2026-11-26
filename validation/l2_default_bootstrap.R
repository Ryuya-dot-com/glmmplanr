# Confirm the actual default interval path after the covariance correction.
source("validation/planning_comparison_helpers.R")
args <- commandArgs(TRUE)
mode <- if (length(args)) match.arg(args[1], c("smoke", "confirm")) else "smoke"
out <- if (length(args) > 1L) args[2] else paste0("validation/l2_default_bootstrap_", mode)
workers <- if (length(args) > 2L) as.integer(args[3]) else 6L
if (.Platform$OS.type == "windows") workers <- 1L
stopifnot(workers >= 1L)
for (e in parse("validation/l2_replication_large/source_script.R")) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.symbol(e[[2]]) &&
      as.character(e[[2]]) %in% c("profiles", "make_spec")) eval(e)
}
chosen <- subset(read.csv("validation/l2_uncertainty_oracle/cases.csv"), case %in% c(31L, 50L, 25L))
jobs <- merge(chosen, data.frame(condition = c("null", "planned")))
jobs$job <- seq_len(nrow(jobs))
nsim <- if (mode == "confirm") 512L else 4L
n_sensitivity <- if (mode == "confirm") 32L else 2L
seed_root <- if (mode == "confirm") 941260001L else 941250001L
boot_root <- if (mode == "confirm") 942260001L else 942250001L
checkpoint_dir(out, c("validation/l2_default_bootstrap.R", "R/marginal_scale.R",
  "validation/l2_replication_large/source_script.R"),
  list(mode = mode, jobs = jobs, nsim = nsim, n_sensitivity = n_sensitivity,
    seed_root = seed_root, boot_root = boot_root, B = 99L, B_sensitivity = 399L,
    coverage_floor = .925, null_ceiling = .075, min_valid_rate = .99,
    point_tolerance_pp = 1e-7, se_tolerance_pp = 1e-7), c("glmmplanr", "digest"))
write.csv(jobs, file.path(out, "cases.csv"), row.names = FALSE)
cat("Loaded package:", find.package("glmmplanr"), "\n")

# Copy only the R call environment to retain the returned bootstrap vector.
# The installed namespace and all R/C++ function bodies remain unchanged.
ns <- asNamespace("glmmplanr")
capture <- new.env(parent = ns)
capture$last_boot <- numeric()
capture$cpp_prob_bootstrap_batch_ws <- function(...) {
  b <- get("cpp_prob_bootstrap_batch_ws", ns)(...)
  capture$last_boot <- b
  b
}
for (name in c("fit_via_cpp_cached", "summarize_probability_fit")) {
  f <- get(name, ns)
  environment(f) <- capture
  assign(name, f, capture)
  stopifnot(identical(body(f), body(get(name, ns))))
}

chunks <- merge(data.frame(job = jobs$job), data.frame(chunk = seq_len(ceiling(nsim / 32L))))
run <- function(k) {
  part <- chunks[k, ]; x <- jobs[jobs$job == part$job, ]
  path <- file.path(out, sprintf("job-%02d-chunk-%02d.rds", part$job, part$chunk))
  if (file.exists(path)) return(readRDS(path)$trials)
  spec <- make_spec(profiles[profiles$profile == x$profile, ])
  if (x$condition == "null") spec <- glmmplanr:::null_target_spec(spec)
  spec_boot <- spec
  spec_boot$prob_se_method <- "bootstrap"
  stopifnot(identical(spec_boot$prob_boot, 99L))
  schedule <- generate_design_schedule(spec$dgp$design, x$n_subj, x$n_item)
  cache <- glmmplanr:::prepare_cpp_schedule_cache(spec, schedule)
  bases <- glmmplanr:::allocate_stage_point_streams(seed_root + 10000L * x$case, 1L, 1:2, 2L)
  streams <- glmmplanr:::derive_rep_substreams(bases[[if (x$condition == "null") 1L else 2L]], nsim)
  boot_streams <- glmmplanr:::make_lecuyer_streams(boot_root + 10000L * x$case + x$job, nsim)
  rows <- draws <- list()
  for (i in seq.int((part$chunk - 1L) * 32L + 1L, min(part$chunk * 32L, nsim))) {
    sim <- glmmplanr:::with_rng_stream(streams[[i]], glmmplanr:::simulate_fast_y_keep(cache))
    jd <- glmmplanr:::fit_via_cpp_cached(sim, cache, spec, TRUE, spec$analysis_re, FALSE)
    point_valid <- glmmplanr:::valid_target_fit(jd$est, jd$se, jd$converged, jd$status,
      jd$prob_est, jd$prob_se, "pp")
    row <- function(method, se, valid, n_boot, n_inner_valid, elapsed, error) {
      data.frame(x, replicate = i, method = method, truth = spec$true_effect_target,
        estimate = 100 * jd$prob_est, se = se, valid = valid,
        point_status = jd$status, error = error, n_boot = n_boot,
        n_inner_valid = n_inner_valid, elapsed = elapsed,
        response_hash = digest::digest(sim$y, algo = "xxhash64"))
    }
    rows[[length(rows) + 1L]] <- row("joint_delta", 100 * jd$prob_se, point_valid,
      0L, 0L, NA_real_, jd$error_reason)
    B <- if (i <= n_sensitivity) 399L else 99L
    spec_boot$prob_boot <- B
    capture$last_boot <- numeric()
    timing <- system.time(bt <- glmmplanr:::with_rng_stream(boot_streams[[i]],
      capture$fit_via_cpp_cached(sim, cache, spec_boot, TRUE, spec$analysis_re, FALSE)))["elapsed"]
    boot <- capture$last_boot
    stopifnot(length(boot) %in% c(0L, B), !point_valid || length(boot) == B)
    if (point_valid) stopifnot(abs(100 * (bt$prob_est - jd$prob_est)) < 1e-7)
    draws[[as.character(i)]] <- boot
    for (nb in if (B == 399L) c(99L, 399L) else 99L) {
      b <- head(boot, nb); good <- sum(is.finite(b))
      valid <- point_valid && length(b) == nb && good >= ceiling(.9 * nb)
      se <- if (valid) 100 * sd(b[is.finite(b)]) else NA_real_
      valid <- valid && is.finite(se) && se > 0
      error <- if (valid) "" else if (!point_valid) jd$error_reason else
        sprintf("Marginal bootstrap: %d/%d valid refits", good, nb)
      rows[[length(rows) + 1L]] <- row(paste0("bootstrap_", nb), se, valid, nb,
        good, unname(timing), error)
      if (nb == B) {
        direct_valid <- glmmplanr:::valid_target_fit(bt$est, bt$se, bt$converged,
          bt$status, bt$prob_est, bt$prob_se, "pp")
        stopifnot(identical(valid, direct_valid))
        if (valid) stopifnot(abs(se - 100 * bt$prob_se) < 1e-7)
      }
    }
    # Independently call the unmodified entry point on fixed probe trials.
    if (i <= 2L) {
      spec99 <- spec_boot; spec99$prob_boot <- 99L
      original <- glmmplanr:::with_rng_stream(boot_streams[[i]],
        glmmplanr:::fit_via_cpp_cached(sim, cache, spec99, TRUE, spec$analysis_re, FALSE))
      probe <- rows[[length(rows) - as.integer(B == 399L)]]
      original_valid <- glmmplanr:::valid_target_fit(original$est, original$se,
        original$converged, original$status, original$prob_est, original$prob_se, "pp")
      stopifnot(identical(probe$valid, original_valid))
      if (original_valid) stopifnot(abs(probe$estimate - 100 * original$prob_est) < 1e-7,
        abs(probe$se - 100 * original$prob_se) < 1e-7)
    }
  }
  trials <- do.call(rbind, rows)
  atomic_rds(list(trials = trials, bootstrap_draws = draws), path)
  cat(sprintf("Default bootstrap %s: job %d chunk %d complete\n", mode, part$job, part$chunk))
  trials
}
result <- parallel::mclapply(seq_len(nrow(chunks)), run, mc.cores = workers,
  mc.set.seed = FALSE, mc.preschedule = FALSE)
stopifnot(all(vapply(result, is.data.frame, logical(1))))
trials <- do.call(rbind, result)
stopifnot(nrow(trials) == nrow(jobs) * (2L * nsim + n_sensitivity),
  !anyDuplicated(trials[c("case", "condition", "replicate", "method")]))
write.csv(trials, gzfile(file.path(out, "fits.csv.gz")), row.names = FALSE)
for (file in c("validation/l2_default_bootstrap.R", "R/engines.R", "src/glmm_laplace.cpp"))
  file.copy(file, file.path(out, paste0("source_", basename(file))), overwrite = TRUE)
cat("Completed", nrow(jobs) * nsim, "outer datasets.\n")
