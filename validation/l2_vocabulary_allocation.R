# Binary meaning recall, with one response per participant-word pair.
# Full runs use 2,048 replications per DGP. See validation/README.md.
library(glmmplanr)

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1L] else "validation"
n_rep <- if (length(args) > 1L) as.integer(args[2L]) else 2048L
workers <- if (length(args) > 2L) as.integer(args[3L]) else
  if (.Platform$OS.type == "windows") 1L else 4L
stopifnot(dir.exists(out), n_rep > 0L, n_rep <= 2048L, workers > 0L)
master_seed <- if (n_rep == 2048L) 85321001L else 85322001L
z_crit <- stats::qnorm(0.975)

grid <- data.frame(n_subj = c(32L, 48L, 64L, 96L, 48L, 64L, 96L),
                   n_item = c(48L, 32L, 24L, 16L, 64L, 48L, 32L))
edges <- subset(grid, n_subj == 32L | n_item == 16L)
rotate <- data.frame(allocation = "latin_square",
                     profile = c("no_slopes", "balanced", "learner_dominant", "word_dominant"),
                     slope_subj = c(0, .25, .40, .15), slope_item = c(0, .25, .15, .40))
groups <- data.frame(allocation = "between_subj",
                     profile = c("no_slopes", "balanced", "word_dominant"),
                     slope_subj = 0, slope_item = c(0, .25, .40))
baseline <- rbind(rotate, groups)
core <- merge(merge(grid, baseline), data.frame(effect_pp = c(0, 8)))
core <- transform(core, phase = "core", rho_subj = 0, rho_item = 0)

balanced <- subset(baseline, profile == "balanced")
rho_subj <- data.frame(rho_subj = c(-.5, .5, 0, 0), rho_item = c(0, 0, -.5, .5))
rho_item <- subset(rho_subj, rho_subj == 0)
correlations <- rbind(merge(subset(balanced, allocation == "latin_square"), rho_subj),
                      merge(subset(balanced, allocation == "between_subj"), rho_item))
correlations <- merge(merge(edges, correlations), data.frame(effect_pp = c(0, 8)))
correlations$phase <- "correlation"
small <- transform(merge(grid, balanced), effect_pp = 5, phase = "small_effect",
                   rho_subj = 0, rho_item = 0)
word_fixed <- transform(merge(edges, data.frame(effect_pp = c(0, 8))),
                        allocation = "between_item", profile = "learner_only",
                        slope_subj = .25, slope_item = 0, phase = "word_assignment",
                        rho_subj = 0, rho_item = 0)
cases <- rbind(core, correlations, small, word_fixed)
cases$case <- seq_len(nrow(cases))
cases$subj_sd <- .8
cases$item_sd <- .6
cases$p_restudy <- .42
cases$p_retrieval <- .42 + cases$effect_pp / 100
cases$responses <- with(cases, n_subj * n_item)
cases$master_seed <- master_seed
cases$n_rep <- n_rep
edge_baseline <- with(cases, phase == "core" & profile == "balanced" &
                        (n_subj == 32L | n_item == 16L))
cases$fit_uncorrelated <- cases$phase == "correlation" | edge_baseline
cases$fit_intercept_only <- edge_baseline
stopifnot(nrow(cases) == 140L, sum(cases$fit_uncorrelated) == 32L,
          sum(cases$fit_intercept_only) == 8L, !anyDuplicated(cases[seq_len(10L)]))

# Same latent draws and uniforms for matching IDs at every size and allocation.
max_subj <- max(grid$n_subj)
max_item <- max(grid$n_item)
streams <- glmmplanr:::make_lecuyer_streams(master_seed, n_rep)
case_dir <- file.path(out, sprintf("l2_vocabulary_checkpoints_%d", master_seed))
dir.create(case_dir, showWarnings = FALSE)
manifest <- file.path(case_dir, "scenarios.rds")
if (file.exists(manifest)) {
  stopifnot(identical(readRDS(manifest), cases))
} else {
  saveRDS(cases, manifest)
}
utils::write.csv(cases, file.path(out, "l2_vocabulary_scenarios.csv"), row.names = FALSE)

run_case <- function(k) {
  cfg <- cases[k, ]
  path <- file.path(case_dir, sprintf("case_%03d.rds", k))
  if (file.exists(path)) return(path)
  design <- design_crossed(list(method = c("restudy", "retrieval")),
                           allocation = cfg$allocation, contrast = "sum")
  dgp <- calibrate_dgp(
    design, c(restudy = cfg$p_restudy, retrieval = cfg$p_retrieval),
    subj_sd = c(cfg$subj_sd, cfg$slope_subj)[seq_len(ncol(design$Z_subj_cells))],
    item_sd = c(cfg$item_sd, cfg$slope_item)[seq_len(ncol(design$Z_item_cells))],
    corr_structure = list(subj = cfg$rho_subj, item = cfg$rho_item)
  )
  spec <- specify_target(
    dgp, "method", scale = "pp", engine = "bayes_laplace", analysis_re = "maximal",
    prob_se_method = "joint_delta", control = target_control(on_singular = "flag_only")
  )
  schedule <- generate_design_schedule(design, cfg$n_subj, cfg$n_item)
  cache <- glmmplanr:::prepare_cpp_schedule_cache(spec, schedule)
  pairs <- cache$sid + (cache$iid - 1L) * max_subj
  is_restudy <- schedule$trials$method == "restudy"
  stopifnot(abs(spec$true_effect_target - cfg$effect_pp) < 1e-7,
            cache$n_tr == cfg$responses, sum(is_restudy) == cache$n_tr / 2L,
            !anyDuplicated(pairs))
  # Verify both calibrated probabilities by independent normal integration.
  cell_variance <- rowSums((design$Z_subj_cells %*% dgp$Sigma_subj) * design$Z_subj_cells) +
    rowSums((design$Z_item_cells %*% dgp$Sigma_item) * design$Z_item_cells)
  p_check <- vapply(seq_len(2L), function(j) {
    stats::integrate(function(z) stats::plogis(dgp$mu_cells[j] + sqrt(cell_variance[j]) * z) *
                       stats::dnorm(z), -Inf, Inf, rel.tol = 1e-10)$value
  }, numeric(1))
  stopifnot(max(abs(p_check - c(cfg$p_restudy, cfg$p_retrieval))) < 1e-7)

  models <- c("maximal", if (cfg$fit_uncorrelated) "uncorrelated",
              if (cfg$fit_intercept_only) "intercept_only")
  draw_group <- function(z, sd0, sd1, rho, ids, x) {
    sd0 * z[ids, 1L] + sd1 * (rho * z[ids, 1L] + sqrt(1 - rho^2) * z[ids, 2L]) * x
  }
  x <- ifelse(is_restudy, -.5, .5)
  fit_rows <- lapply(seq_len(n_rep), function(i) {
    y <- glmmplanr:::with_rng_stream(streams[[i]], {
      subj_z <- matrix(stats::rnorm(2L * max_subj), max_subj, 2L)
      item_z <- matrix(stats::rnorm(2L * max_item), max_item, 2L)
      eta <- cache$mu_trial +
        draw_group(subj_z, cfg$subj_sd, cfg$slope_subj, cfg$rho_subj, cache$sid, x) +
        draw_group(item_z, cfg$item_sd, cfg$slope_item, cfg$rho_item, cache$iid, x)
      as.integer(stats::runif(max_subj * max_item)[pairs] < stats::plogis(eta))
    })
    fitted <- lapply(models, function(model) {
      f <- glmmplanr:::fit_via_cpp_cached(list(y = y, keep = NULL), cache, spec,
                                         regularized = TRUE, re_structure = model,
                                         fast_stage = FALSE)
      valid <- glmmplanr:::valid_target_fit(f$est, f$se, f$converged, f$status,
                                            f$prob_est, f$prob_se, "pp")
      est <- if (valid) 100 * f$prob_est else NA_real_
      se <- if (valid) 100 * f$prob_se else NA_real_
      data.frame(model = model, valid = valid, estimate_pp = est, se_pp = se,
                 covered = if (valid) abs(est - cfg$effect_pp) <= z_crit * se else NA,
                 detected = valid && abs(est) > z_crit * se, status = f$status,
                 diagnostic = if (valid) "" else f$error_reason)
    })
    cbind(case = k, replication = i, p_restudy_observed = mean(y[is_restudy]),
          p_retrieval_observed = mean(y[!is_restudy]), do.call(rbind, fitted))
  })
  rows <- do.call(rbind, fit_rows)
  stopifnot(nrow(rows) == n_rep * length(models), !anyDuplicated(rows[c("replication", "model")]))
  tmp <- paste0(path, ".new")
  saveRDS(rows, tmp)
  stopifnot(file.rename(tmp, path))
  cat(sprintf("Case %3d/%d %-15s %-14s %2dx%2d %+d pp: %d/%d valid fits\n",
              k, nrow(cases), cfg$phase, cfg$allocation, cfg$n_subj, cfg$n_item,
              cfg$effect_pp, sum(rows$valid), nrow(rows)))
  flush.console()
  path
}

exact_ci <- function(k, n) {
  if (!n) return(c(NA_real_, NA_real_))
  stats::binom.test(k, n)$conf.int[1:2]
}

summarize_case <- function(rows) {
  cfg <- cases[rows$case[1L], ]
  n <- nrow(rows)
  n_valid <- sum(rows$valid)
  covered <- sum(rows$covered, na.rm = TRUE)
  detected <- sum(rows$detected)
  cov_ci <- exact_ci(covered, n_valid)
  det_ci <- exact_ci(detected, n)
  coverage_status <- if (n_valid / n < .90 || cov_ci[2L] < .90) "fail" else
    if (cov_ci[1L] >= .90) "pass" else "indeterminate"
  detection_status <- if (cfg$effect_pp == 0) {
    if (det_ci[2L] <= .075) "pass" else if (det_ci[1L] > .075) "fail" else "indeterminate"
  } else if (det_ci[1L] >= .80) "pass" else if (det_ci[2L] < .80) "fail" else "indeterminate"
  width <- 2 * z_crit * rows$se_pp[rows$valid]
  err <- rows$estimate_pp[rows$valid] - cfg$effect_pp
  cbind(cfg, model = rows$model[1L], attempted = n, valid = n_valid,
        covered = covered, detected = detected, coverage = if (n_valid) covered / n_valid else NA,
        coverage_lo = cov_ci[1L], coverage_hi = cov_ci[2L], coverage_status = coverage_status,
        detection = detected / n, detection_lo = det_ci[1L], detection_hi = det_ci[2L],
        detection_status = detection_status, bias_pp = mean(err), sd_pp = stats::sd(err),
        mean_se_pp = mean(rows$se_pp[rows$valid]), mean_width_pp = mean(width),
        width_q80_pp = unname(stats::quantile(width, .80)),
        p_restudy_observed = mean(rows$p_restudy_observed),
        p_retrieval_observed = mean(rows$p_retrieval_observed))
}

cat(sprintf("Master seed %d, %d datasets and %d fits per replication, %d replications, %d workers\n",
            master_seed, nrow(cases), nrow(cases) + 40L, n_rep, workers))
cat("Package:", normalizePath(find.package("glmmplanr")), "\n")
paths <- parallel::mclapply(seq_len(nrow(cases)), run_case, mc.cores = workers,
                            mc.preschedule = FALSE, mc.set.seed = FALSE)
stopifnot(all(vapply(paths, function(x) is.character(x) && length(x) == 1L &&
                      file.exists(x), logical(1))))
results <- lapply(paths, readRDS)
fits <- do.call(rbind, results)
stopifnot(!anyDuplicated(fits[c("case", "replication", "model")]),
          nrow(fits) == 180L * n_rep)
summary <- do.call(rbind, lapply(results, function(rows) {
  do.call(rbind, lapply(split(rows, rows$model), summarize_case))
}))
utils::write.csv(summary, file.path(out, "l2_vocabulary_summary.csv"), row.names = FALSE)
connection <- gzfile(file.path(out, "l2_vocabulary_fits.csv.gz"), "wt")
utils::write.csv(fits, connection, row.names = FALSE)
close(connection)
cat(sprintf("Completed %d datasets and %d fits.\n", nrow(cases) * n_rep, nrow(fits)))
