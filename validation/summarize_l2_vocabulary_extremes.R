# Recompute fit counts and quantify paired differences between layouts and models.
args <- commandArgs(TRUE)
out <- if (length(args)) args[1L] else "validation"
stage <- if (length(args) > 1L) match.arg(args[2L], c("pilot", "large")) else "pilot"
prefix <- file.path(out, paste0("l2_vocabulary_extremes_", stage))
cases <- read.csv(paste0(prefix, "_scenarios.csv"))
summary <- read.csv(paste0(prefix, "_summary.csv"))
fits <- read.csv(paste0(prefix, "_fits.csv.gz"))
keys <- paste(fits$case, fits$model)
rows <- split(seq_len(nrow(fits)), keys)
z <- qnorm(.975)
exact_ci <- function(k, n) if (n > 0L) binom.test(k, n)$conf.int[1:2] else rep(NA_real_, 2L)
stopifnot(nrow(cases) == 304L, length(rows) == nrow(summary),
          !anyDuplicated(fits[c("case", "replication", "model")]))
for (k in seq_len(nrow(summary))) {
  s <- summary[k, ]
  f <- fits[rows[[paste(s$case, s$model)]], ]
  cov_ci <- exact_ci(sum(f$covered, na.rm = TRUE), sum(f$valid))
  det_ci <- exact_ci(sum(f$detected), nrow(f))
  stopifnot(identical(f$replication, seq_len(s$n_rep)), s$n_rep == s$attempted,
            sum(f$valid) == s$valid, sum(f$detected) == s$detected,
            sum(f$covered, na.rm = TRUE) == s$covered,
            sum(f$method_separation) == s$method_separated,
            isTRUE(all.equal(c(cov_ci, det_ci),
                             c(s$coverage_lo, s$coverage_hi, s$detection_lo, s$detection_hi))),
            all(f$detected[f$valid] == (abs(f$estimate_pp[f$valid]) > z * f$se_pp[f$valid])),
            all(f$covered[f$valid] == (abs(f$estimate_pp[f$valid] - s$true_pp) <= z * f$se_pp[f$valid])),
            abs(mean(f$estimate_pp[f$valid]) - s$true_pp - s$bias_pp) < 1e-10,
            abs(2 * z * mean(f$se_pp[f$valid]) - s$mean_width_pp) < 1e-10)
}

# Check the simulator against the independently integrated arm probabilities.
standardized_error <- function(observed, truth) {
  abs(mean(observed) - truth) / (sd(observed) / sqrt(length(observed)))
}
truth_errors <- vapply(seq_len(nrow(cases)), function(k) {
  f <- fits[rows[[paste(k, "unadjusted")]], ]
  c(standardized_error(f$p_restudy_obs, cases$p_restudy[k]),
    standardized_error(f$p_retrieval_obs, cases$p_retrieval[k]),
    standardized_error(f$p_retrieval_obs - f$p_restudy_obs, cases$true_pp[k] / 100))
}, numeric(3L))
stopifnot(max(truth_errors) < 6)

comparisons <- list()
pair <- function(left, right, name, model_left = "unadjusted", model_right = "unadjusted") {
  a <- fits[rows[[paste(left$case, model_left)]], ]
  b <- fits[rows[[paste(right$case, model_right)]], ]
  stopifnot(nrow(a) == nrow(b), identical(a$replication, b$replication),
            abs(left$true_pp - right$true_pp) < 1e-7)
  valid <- a$valid & b$valid
  # Report Monte Carlo uncertainty in the paired difference, in percentage points.
  mc <- function(dx) {
    diff <- mean(dx)
    c(diff, diff + c(-1, 1) * z * sd(dx) / sqrt(length(dx)))
  }
  d <- rbind(detection = mc(100 * (a$detected - b$detected)),
             coverage = mc(100 * (a$covered[valid] - b$covered[valid])),
             mean_width_pp = mc(2 * z * (a$se_pp[valid] - b$se_pp[valid])),
             bias_pp = mc(a$estimate_pp[valid] - b$estimate_pp[valid]))
  comparisons[[length(comparisons) + 1L]] <<- data.frame(
    comparison = name, case_left = left$case, case_right = right$case, model_left, model_right,
    p_restudy = left$p_restudy, effect = left$effect, tail_left = left$tail, tail_right = right$tail,
    slopes = left$slopes, covariates = left$covariates, left_id = left$id, right_id = right$id,
    n = nrow(a), both_valid = sum(valid), measure = rownames(d),
    difference = d[, 1L], lo = d[, 2L], hi = d[, 3L], row.names = NULL
  )
}

layouts <- data.frame(
  comparison = c("larger_bank", "all_words_vs_partial", "more_words_per_learner", "more_learners"),
  left = c("96_48_16", "32_48_48", "96_48_24", "144_48_16"),
  right = c("96_16_16", "96_48_16", "96_48_8", "48_48_16")
)
for (k in seq_len(nrow(layouts))) {
  l <- layouts[k, ]
  a <- subset(cases, id == l$left & covariates == "none")
  for (j in seq_len(nrow(a))) {
    a0 <- a[j, ]
    b <- subset(cases, id == l$right & p_restudy == a0$p_restudy & effect == a0$effect &
                  tail == a0$tail & slopes == a0$slopes & covariates == "none")
    if (nrow(b) == 1L) pair(a0, b, l$comparison)
  }
}
for (k in which(cases$phase == "local")) {
  a <- cases[k, ]
  group <- sub("_(high|low)$", "", a$tail)
  b <- subset(cases, id == a$id & effect == a$effect & tail == paste0(group, "_gaussian"))
  if (nrow(b) == 1L) pair(a, b, "mixture_vs_matched_gaussian")
}
for (k in which(cases$covariates != "none" | cases$fit_int_only)) {
  a <- cases[k, ]
  if (a$fit_int_only) pair(a, a, "omission_of_slopes", "intercept_only", "unadjusted")
  if (a$covariates != "none") pair(a, a, "covariate_adjustment", "adjusted", "unadjusted")
  if (a$method_wm != 0 | a$method_pos != 0) pair(a, a, "add_effect_modification", "moderation", "adjusted")
}
write.csv(do.call(rbind, comparisons), paste0(prefix, "_paired.csv"), row.names = FALSE)
cat(sprintf("Checked %d scenarios, %d scenario-model summaries, and %d individual fits.\n",
            nrow(cases), nrow(summary), nrow(fits)))
cat(sprintf("Largest arm or contrast deviation from the integrated truth: %.2f Monte Carlo SE.\n",
            max(truth_errors)))
