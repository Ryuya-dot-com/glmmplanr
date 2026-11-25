# Pre-specified null false-positive & allocation check for marginal probability intervals
# Matches coverage_recheck.R population (64 subj x 20 items, subj_sd = 0.5, item_sd = 0.3)
# under H0 (60% vs 60%) across between_subj and latin_square allocations, comparing
# C++ parametric bootstrap (B = 99) and analytical joint_delta on the same responses.
library(glmmplanr)
RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")

n_rep <- 96L
layouts <- c("between_subj", "latin_square")

run_layout <- function(layout_nm, base_id) {
  des <- design_crossed(list(condition = c("A", "B")), allocation = layout_nm, contrast = "treatment")
  dgp <- calibrate_dgp(
    des, cell_probs = c(A = 0.60, B = 0.60), chance = 0,
    subj_sd = 0.5, item_sd = 0.3, slope_sd_ratio = 0, corr_structure = 0,
    control = calibrate_control(seed = 417L)
  )
  sched <- generate_design_schedule(des, n_subj = 64L, n_item = 20L)

  sp_boot <- specify_target(
    dgp, effect = "condition", scale = "pp",
    engine = "bayes_laplace", analysis_re = "intercept_only",
    prob_se_method = "bootstrap",
    control = target_control(prob_estimand = "marginal", prob_boot = 99L, on_singular = "flag_only")
  )
  sp_jd <- specify_target(
    dgp, effect = "condition", scale = "pp",
    engine = "bayes_laplace", analysis_re = "intercept_only",
    prob_se_method = "joint_delta",
    control = target_control(prob_estimand = "marginal", on_singular = "flag_only")
  )

  rows <- lapply(seq_len(n_rep), function(r) {
    rep_seed <- base_id + r
    set.seed(rep_seed, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion", sample.kind = "Rejection")
    tr <- simulate_dgp_data(dgp, schedule = sched)

    set.seed(rep_seed + 500000L, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion", sample.kind = "Rejection")
    f_bt <- glmmplanr:::fit_single_dataset(tr, sp_boot, engine_use = "bayes_laplace", re_structure = "intercept_only")
    f_jd <- glmmplanr:::fit_single_dataset(tr, sp_jd, engine_use = "bayes_laplace", re_structure = "intercept_only")

    make_row <- function(method_nm, f) {
      valid <- isTRUE(f$converged) && is.finite(f$prob_est) && is.finite(f$prob_se) && f$prob_se > 0
      est_pp <- if (valid) f$prob_est * 100 else NA_real_
      se_pp  <- if (valid) f$prob_se * 100 else NA_real_
      z_crit <- stats::qnorm(0.975)
      false_pos <- valid && (abs(est_pp) > z_crit * se_pp)
      covered   <- if (valid) (abs(est_pp) <= z_crit * se_pp) else NA
      data.frame(
        allocation = layout_nm,
        se_method = method_nm,
        replication = r,
        seed = rep_seed,
        valid = valid,
        est_pp = est_pp,
        se_pp = se_pp,
        width_pp = if (valid) 2 * z_crit * se_pp else NA_real_,
        false_pos = false_pos,
        covered = covered,
        stringsAsFactors = FALSE
      )
    }
    rbind(make_row("bootstrap_99", f_bt), make_row("joint_delta", f_jd))
  })
  do.call(rbind, rows)
}

fits <- rbind(
  run_layout("between_subj", 20000L),
  run_layout("latin_square", 21000L)
)

by_grp <- split(fits, list(fits$allocation, fits$se_method), drop = TRUE)
summary_df <- do.call(rbind, lapply(by_grp, function(d) {
  n_att <- nrow(d)
  n_val <- sum(d$valid)
  k_fp  <- sum(d$false_pos)
  k_cov <- sum(d$covered, na.rm = TRUE)
  fp_ci <- as.numeric(stats::binom.test(k_fp, n_att)$conf.int)
  cov_ci <- if (n_val > 0L) as.numeric(stats::binom.test(k_cov, n_val)$conf.int) else c(NA_real_, NA_real_)
  data.frame(
    allocation = d$allocation[1L],
    se_method = d$se_method[1L],
    attempted = n_att,
    valid = n_val,
    false_positives = k_fp,
    type1_rate = k_fp / n_att,
    type1_lo = fp_ci[1L],
    type1_hi = fp_ci[2L],
    covered = k_cov,
    null_coverage = k_cov / n_val,
    coverage_lo = cov_ci[1L],
    coverage_hi = cov_ci[2L],
    empirical_sd_pp = stats::sd(d$est_pp[d$valid]),
    mean_se_pp = mean(d$se_pp[d$valid]),
    mean_width_pp = mean(d$width_pp[d$valid]),
    stringsAsFactors = FALSE
  )
}))
summary_df <- summary_df[order(summary_df$allocation, summary_df$se_method), ]
rownames(summary_df) <- NULL

write.csv(fits, "validation/marginal_null_recheck_fits.csv", row.names = FALSE)
write.csv(summary_df, "validation/marginal_null_recheck_summary.csv", row.names = FALSE)
print(summary_df, row.names = FALSE, digits = 4)
