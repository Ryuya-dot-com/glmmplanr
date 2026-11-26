# Compare fixed and random passage analyses on the same reader-randomized trials.
library(glmmplanr)
args <- commandArgs(TRUE)
out <- if (length(args)) args[1L] else "validation"
stage <- if (length(args) > 1L) match.arg(args[2L], c("pilot", "medium")) else "pilot"
workers <- if (length(args) > 2L) as.integer(args[3L]) else 4L
stopifnot(dir.exists(out), workers > 0L)
if (.Platform$OS.type == "windows") workers <- 1L
n_rep <- if (stage == "pilot") 8L else 256L
seed <- if (stage == "pilot") 91641001L else 91642001L
streams <- glmmplanr:::make_lecuyer_streams(seed, n_rep)
scenarios <- function(n, passages, slope_sd, local_sd) {
  expand.grid(n = n, passages = passages, slope_sd = slope_sd, local_sd = local_sd,
              effect_pp = c(0, 10), KEEP.OUT.ATTRS = FALSE)
}
cases <- rbind(scenarios(48L, 2:6, c(0, .4), c(0, .4)), scenarios(96L, c(2L, 5L), .4, .4))
cases$case <- seq_len(nrow(cases))
cases$master_seed <- seed
cases$n_rep <- n_rep
cases$n_local_check <- with(cases, ifelse(n == 48L & passages %in% c(2L, 5L) &
                                          slope_sd == .4 & local_sd == .4, min(n_rep, 64L), 0L))
stopifnot(nrow(cases) == 44L, sum(cases$n_local_check) == 4L * min(n_rep, 64L))
source("validation/few_shared_passages_helpers.R")
cases <- calibrate_shared_passage_cases(cases)

run_case <- function(k) {
  cfg <- cases[k, ]
  sched <- schedule(cfg$n, cfg$passages)
  random_cache <- glmmplanr:::prepare_cpp_schedule_cache(random_spec, sched)
  models <- list(fixed_common = fixed_setup(sched, common = TRUE),
    fixed_interaction = fixed_setup(sched, common = FALSE),
    random_intercept = list(spec = random_spec, cache = random_cache, structure = "intercept_only"),
    random_slope = list(spec = random_spec, cache = random_cache, structure = "uncorrelated"))
  if (cfg$n_local_check) models$fixed_interaction_local <- fixed_setup(sched, common = FALSE, local = TRUE)
  do.call(rbind, lapply(seq_len(n_rep), function(r) {
    sim <- draw_shared_passage_response(cfg, sched, streams[[r]], max(cases$n), max(cases$passages))
    finite <- sim$truth_sampled_pp
    use <- if (r <= cfg$n_local_check) models else models[seq_len(4L)]
    do.call(rbind, lapply(names(use), function(name) {
      fit <- fit_response(sim$y, use[[name]])
      if (is.null(fit$status)) stop("Missing native fit status")
      hw <- zcrit * fit$se_pp
      data.frame(case = k, replication = r, model = name, truth_sampled_pp = finite,
                 fit, covered_sampled = if (fit$valid) abs(fit$estimate_pp - finite) <= hw else NA,
                 covered_population = if (fit$valid) abs(fit$estimate_pp - cfg$effect_pp) <= hw else NA,
                 detected = fit$valid && abs(fit$estimate_pp) > hw)
    }))
  }))
}

prefix <- file.path(out, paste0("few_shared_passages_", stage))
checkpoint <- file.path(out, "checkpoints")
dir.create(checkpoint, showWarnings = FALSE)
manifest <- file.path(checkpoint, "scenarios.rds")
if (file.exists(manifest)) stopifnot(identical(readRDS(manifest), cases)) else saveRDS(cases, manifest)
write.csv(cases, paste0(prefix, "_scenarios.csv"), row.names = FALSE)
cat("Package:", normalizePath(find.package("glmmplanr")), "\n")
started <- proc.time()[["elapsed"]]
paths <- parallel::mclapply(seq_len(nrow(cases)), function(k) {
  path <- file.path(checkpoint, sprintf("case_%03d.rds", k))
  if (file.exists(path)) return(path)
  rows <- run_case(k)
  stopifnot(nrow(rows) == 4L * n_rep + cases$n_local_check[k])
  saveRDS(rows, paste0(path, ".new"))
  stopifnot(file.rename(paste0(path, ".new"), path))
  cat(sprintf("Case %d/%d, %d valid of %d\n", k, nrow(cases), sum(rows$valid), nrow(rows)))
  flush.console()
  path
}, mc.cores = workers, mc.set.seed = FALSE, mc.preschedule = FALSE)
completed <- vapply(paths, function(path) is.character(path) && length(path) == 1L && file.exists(path), logical(1L))
if (!all(completed)) {
  for (k in which(!completed)) message("Case ", k, ": ", paths[[k]])
  stop("Some cases failed. Completed cases are saved in ", checkpoint)
}
fits <- do.call(rbind, lapply(paths, readRDS))
stopifnot(!anyDuplicated(fits[c("case", "replication", "model")]),
          nrow(fits) == 4L * nrow(cases) * n_rep + sum(cases$n_local_check))
exact <- function(x, n) if (n) unname(binom.test(x, n)$conf.int) else c(NA_real_, NA_real_)
report <- function(rows) {
  cfg <- cases[rows$case[1L], ]
  m <- rows$model[1L]
  valid <- rows$valid
  primary <- if (startsWith(m, "fixed")) "sampled" else "population"
  sampled <- sum(rows$covered_sampled, na.rm = TRUE)
  population <- sum(rows$covered_population, na.rm = TRUE)
  detected <- sum(rows$detected)
  nv <- sum(valid)
  n <- nrow(rows)
  cis <- exact(sampled, nv)
  cip <- exact(population, nv)
  cid <- exact(detected, n)
  ci <- if (primary == "sampled") cis else cip
  coverage_gate <- if (nv / n < .90 || ci[2L] < .90) "fail" else if (ci[1L] >= .90) "pass" else "indeterminate"
  detect_gate <- if (cfg$effect_pp == 0 && primary == "sampled" && cfg$slope_sd > 0) "not_applicable" else
    if (cfg$effect_pp == 0) {
      if (cid[2L] <= .075) "pass" else if (cid[1L] > .075) "fail" else "indeterminate"
    } else if (primary == "population") {
      if (cid[1L] >= .80) "pass" else if (cid[2L] < .80) "fail" else "indeterminate"
    } else "descriptive"
  cbind(cfg, model = m, primary = primary, attempted = n, valid = nv, covered_sampled = sampled,
    covered_population = population, detected = detected, detection_lo = cid[1L], detection_hi = cid[2L],
    sampled_lo = cis[1L], sampled_hi = cis[2L], population_lo = cip[1L], population_hi = cip[2L],
    coverage_gate, detection_gate = detect_gate, mean_sampled_truth = mean(rows$truth_sampled_pp),
    sd_sampled_truth = sd(rows$truth_sampled_pp), bias_sampled = mean(rows$estimate_pp[valid] - rows$truth_sampled_pp[valid]),
    bias_population = mean(rows$estimate_pp[valid] - cfg$effect_pp), empirical_sd = sd(rows$estimate_pp[valid]),
    sampled_error_sd = sd(rows$estimate_pp[valid] - rows$truth_sampled_pp[valid]),
    mean_se = mean(rows$se_pp[valid]), mean_width = 2 * zcrit * mean(rows$se_pp[valid]),
    median_slope_sd_hat = median(rows$slope_sd_hat[valid], na.rm = TRUE))
}
summ <- do.call(rbind, lapply(split(fits, fits[c("case", "model")], drop = TRUE), report))
write.csv(summ, paste0(prefix, "_summary.csv"), row.names = FALSE)
con <- gzfile(paste0(prefix, "_fits.csv.gz"), "wt")
write.csv(fits, con, row.names = FALSE)
close(con)
cat(sprintf("Completed %d datasets, %d fits (%d valid), in %.1f seconds.\n", nrow(cases) * n_rep,
            nrow(fits), sum(fits$valid), proc.time()[["elapsed"]] - started))
