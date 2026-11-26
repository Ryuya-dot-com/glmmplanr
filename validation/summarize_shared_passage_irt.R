# Audit rejection and coverage from saved estimates and form paired contrasts.
args <- commandArgs(TRUE)
out <- if (length(args)) args[1L] else "validation"
stage <- if (length(args) > 1L) match.arg(args[2L], c("pilot", "large")) else "large"
part <- if (length(args) > 2L) match.arg(args[3L], c("population", "bank", "mirt")) else "population"
prefix <- file.path(out, paste("shared_passage_irt", stage, part, sep = "_"))
cases <- read.csv(paste0(prefix, "_scenarios.csv"))
fits <- read.csv(paste0(prefix, "_fits.csv.gz"))
z <- qnorm(.975)
attempts <- switch(part, population = cases$n_rep, bank = cases$n_bank, mirt = pmin(cases$n_bank, 128L))
models <- if (part == "population") c("native_intercept", "native_slope", "tmb_intercept", "tmb_slope") else
  c("native_bank", "tam_bank", "tmb_bank")
if (part == "mirt") {
  bank <- read.csv(file.path(out, paste0("shared_passage_irt_", stage, "_bank_fits.csv.gz")))
  fits <- rbind(fits, bank[bank$replication <= attempts[bank$case], ])
  models <- c(models, "mirt_bank")
}
stopifnot(!anyDuplicated(fits[c("case", "replication", "model")]),
          nrow(fits) == length(models) * sum(attempts), setequal(unique(fits$model), models),
          all(fits$target == if (part == "population") "population" else "bank"))
ci <- function(x, n) if (n) unname(binom.test(x, n)$conf.int) else rep(NA_real_, 2L)
mc_mean <- function(x) {
  half_width <- z * sd(x) / sqrt(length(x))
  c(difference = mean(x), lo = mean(x) - half_width, hi = mean(x) + half_width)
}
mc_rate <- function(x) {
  if (!length(x)) return(c(difference = NA_real_, lo = NA_real_, hi = NA_real_))
  positive <- binom.test(sum(x > 0), length(x), conf.level = .975)$conf.int
  negative <- binom.test(sum(x < 0), length(x), conf.level = .975)$conf.int
  c(difference = 100 * mean(x), lo = 100 * (positive[1L] - negative[2L]),
    hi = 100 * (positive[2L] - negative[1L]))
}
groups <- split(fits, paste(fits$case, fits$model))
fit_key <- paste(fits$case, fits$replication)
valid_by_dataset <- table(fit_key[fits$valid])
common_point_keys <- names(valid_by_dataset)[valid_by_dataset == length(models)]
summary <- lapply(groups, function(f) {
  cfg <- cases[f$case[1], ]
  m <- f$model[1L]
  stopifnot(identical(f$replication, seq_len(attempts[cfg$case])),
            all(!f$interval | f$valid), all(is.finite(f$estimate_pp[f$valid])),
            all(f$se_pp[f$interval] > 0),
            !any(f$detected[!f$interval]), all(is.na(f$covered[!f$interval])))
  if (part == "population") stopifnot(all(f$truth_pp == cfg$effect_pp))
  valid_int <- f[f$interval, ]
  valid_point <- f[f$valid, ]
  ni <- nrow(valid_int)
  nv <- nrow(valid_point)
  n <- nrow(f)
  err <- valid_point$estimate_pp - valid_point$truth_pp
  common <- paste(valid_point$case, valid_point$replication) %in% common_point_keys
  ct <- abs(valid_int$estimate_pp - valid_int$truth_pp) <= z * valid_int$se_pp
  dt <- abs(valid_int$estimate_pp) > z * valid_int$se_pp
  stopifnot(all(ct == valid_int$covered), all(dt == valid_int$detected))
  cc <- ci(sum(ct), ni)
  cd <- if (m == "tam_bank") c(NA_real_, NA_real_) else ci(sum(dt), n)
  cov_gate <- if (!ni) "no_interval" else if (ni / n < .9 || cc[2L] < .9) "fail" else
    if (cc[1L] >= .9) "pass" else "indeterminate"
  det_gate <- if (part != "population") "not_applicable" else if (cfg$effect_pp == 0) {
    if (cd[2L] <= .075) "pass" else if (cd[1L] > .075) "fail" else "indeterminate"
  } else if (cd[1L] >= .8) "pass" else if (cd[2L] < .8) "fail" else "indeterminate"
  cbind(cfg, model = m, scheduled = n, valid_points = nv, intervals = ni,
    boundary = if (all(is.na(f$boundary))) NA_integer_ else sum(f$boundary, na.rm = TRUE),
    valid_boundary = if (all(is.na(f$boundary))) NA_integer_ else sum(f$boundary & f$valid, na.rm = TRUE),
    mean_truth_pp = mean(f$truth_pp), bias_pp = mean(err), rmse_pp = sqrt(mean(err^2)),
    common_points = sum(common), common_rmse_pp = sqrt(mean(err[common]^2)),
    empirical_sd_pp = sd(err), mean_interval_width_pp = 2 * z * mean(valid_int$se_pp),
    covered = if (ni) sum(ct) else NA_integer_, coverage_lo = cc[1L], coverage_hi = cc[2L], coverage_gate = cov_gate,
    detected = if (m == "tam_bank") NA_integer_ else sum(dt), detection_lo = cd[1L], detection_hi = cd[2L],
    detection_gate = det_gate, median_slope_sd = median(valid_point$slope_sd_hat, na.rm = TRUE),
    median_local_sd = median(valid_point$local_sd_hat, na.rm = TRUE),
    median_cpu_s = median(f$cpu_s), total_cpu_s = sum(f$cpu_s), row.names = NULL)
})
paired <- lapply(which(attempts > 0), function(k) {
  compare <- function(left, right) {
    a <- groups[[paste(k, left)]]
    b <- groups[[paste(k, right)]]
    stopifnot(identical(a$replication, b$replication), identical(a$truth_pp, b$truth_pp))
    vp <- a$valid & b$valid
    vi <- a$interval & b$interval
    d <- rbind(estimate_pp = mc_mean(a$estimate_pp[vp] - b$estimate_pp[vp]),
      mse_pp_squared = mc_mean((a$estimate_pp[vp] - a$truth_pp[vp])^2 - (b$estimate_pp[vp] - b$truth_pp[vp])^2))
    if (any(vi)) {
      d <- rbind(d, width_pp = mc_mean(2 * z * (a$se_pp[vi] - b$se_pp[vi])),
        coverage_pp = mc_rate(a$covered[vi] - b$covered[vi]),
        detection_pp = mc_rate(a$detected - b$detected))
    }
    cbind(cases[k, ], left, right, paired_points = sum(vp), paired_intervals = sum(vi),
          measure = rownames(d), d, row.names = NULL)
  }
  if (part == "population") {
    rbind(compare("native_slope", "native_intercept"), compare("tmb_slope", "tmb_intercept"),
          compare("native_intercept", "tmb_intercept"), compare("native_slope", "tmb_slope"))
  } else if (part == "bank") {
    rbind(compare("native_bank", "tam_bank"), compare("tmb_bank", "tam_bank"), compare("native_bank", "tmb_bank"))
  } else {
    rbind(compare("mirt_bank", "native_bank"), compare("mirt_bank", "tmb_bank"), compare("mirt_bank", "tam_bank"))
  }
})
write.csv(do.call(rbind, summary), paste0(prefix, "_summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, paired), paste0(prefix, "_paired.csv"), row.names = FALSE)
cat(sprintf("Verified %d scenarios, %d scenario-model summaries, %d rows.\n",
            sum(attempts > 0), length(summary), nrow(fits)))
