# Pre-specified independent validation for sparse partial presentation and testlet designs
# (Addresses exploratory-scope gaps identified in the 0.1.0 validation matrix).
#
# Protocol fixed prior to execution (no data-dependent stopping or seed filtering):
# - Outer replications per cell: n_rep = 96 (nominal 95% coverage MCSE ~ 2.2 pp).
# - Validity gate: valid >= ceil(0.90 * attempted) (at least 87/96 valid fits).
# - Coverage gate (valid fits in denominator): exact binomial 95% CI lower bound >= 0.90 -> "pass",
#   upper bound < 0.90 -> "below_floor", otherwise -> "indeterminate".
# - Null Type I error (all 96 scheduled datasets in denominator): two-sided 5% test reported with
#   exact binomial 95% CI.
#
# Regimes evaluated:
# 1. partial_intercept: 64 subj (32 vs 32 between_subj) x 20-item bank, 5 items/subj (320 trials),
#    chance = 0, subj_sd = 0.5, item_sd = 0.3, planned 60% vs 55% (-5 pp) and null 60% vs 60% (0 pp).
#    Paired comparison of bayes_laplace (bootstrap B = 99), bayes_laplace (joint_delta), and
#    external lme4::glmer parametric bootstrap (B = 99). Seeds: 30001-30096 (null), 31001-31096 (planned).
# 2. partial_correlated_slopes: 64 subj x 24-item bank, 8 items/subj (512 trials), within-subj/item
#    Latin square, chance = 0, subj_sd = 0.8, item_sd = 0.6, slope_sd_ratio = 0.30, corr = 0.20,
#    planned 60% vs 70% (+10 pp) and null 60% vs 60% (0 pp). Evaluates maximal bayes_laplace on both
#    conditional logit and marginal pp (bootstrap B = 99 & joint_delta). Seeds: 40001-40096, 41001-41096.
# 3. testlet_nesting: 48 subj x 24 items (6 passages x 4 items = 1,152 trials), between-passage
#    condition assignment, chance = 0, subj_sd = 0.7, item_sd = 0.4, testlet_sd = 0.5, under
#    (a) matched generating model (participant_testlet_sd = 0.0; seeds 50001-50096, 51001-51096) and
#    (b) unmodeled participant-by-passage local dependence (participant_testlet_sd = 0.4; seeds
#    52001-52096, 53001-53096), comparing bayes_laplace (logit & joint_delta pp) and lme4::glmer (logit).

library(glmmplanr)
library(lme4)

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1L) args[1L] else "validation"
stopifnot(dir.exists(out_dir))

n_rep <- 96L
n_boot <- 99L
workers <- if (.Platform$OS.type == "windows") 1L else 4L
z_crit <- stats::qnorm(0.975)

glmer_marginal_boot <- function(d, refits = n_boot) {
  singular <- NA
  diagnostic <- ""
  ci <- tryCatch({
    model <- suppressMessages(glmer(
      y ~ group + (1 | subj) + (1 | item),
      data = d, family = binomial(),
      control = glmerControl(optimizer = "bobyqa")
    ))
    singular <- isSingular(model)
    effect_fn <- function(f) {
      conv <- f@optinfo$conv
      msg <- conv$lme4$messages
      if (any(conv$opt != 0L) || any(!grepl("boundary (singular) fit", msg, fixed = TRUE))) {
        return(NA_real_)
      }
      sd_re <- sqrt(sum(vapply(VarCorr(f), function(v) v[1, 1], numeric(1))))
      mu <- drop(getME(f, "X")[match(c("control", "treatment"), d$group), , drop = FALSE] %*% fixef(f))
      stats::integrate(
        function(z) (stats::plogis(mu[2] + sd_re * z) - stats::plogis(mu[1] + sd_re * z)) * stats::dnorm(z),
        -9, 9, abs.tol = 1e-9
      )$value
    }
    est <- effect_fn(model)
    if (!is.finite(est)) stop("Base glmer fit did not converge")
    boot <- bootMer(model, effect_fn, nsim = refits, type = "parametric", use.u = FALSE)
    usable <- is.finite(boot$t)
    if (sum(usable) < ceiling(0.9 * refits)) stop("Too few valid glmer bootstrap refits")
    c(estimate = est, se = stats::sd(boot$t[usable]))
  }, error = function(e) {
    diagnostic <<- conditionMessage(e)
    c(estimate = NA_real_, se = NA_real_)
  })
  list(estimate = unname(ci[1L]), se = unname(ci[2L]), singular = singular, diagnostic = diagnostic)
}

make_fit_row <- function(regime, condition, estimator, scale_nm, rep_id, seed, trials, true_val,
                         est, se, singular, diagnostic = "") {
  valid <- is.finite(est) && is.finite(se) && (se > 0)
  half <- if (valid) z_crit * se else NA_real_
  data.frame(
    regime = regime,
    condition = condition,
    estimator = estimator,
    scale = scale_nm,
    replication = rep_id,
    seed = seed,
    trials = trials,
    true_target = true_val,
    valid = valid,
    estimate = if (valid) est else NA_real_,
    se = if (valid) se else NA_real_,
    width = if (valid) 2 * half else NA_real_,
    detected = valid && (abs(est) > half),
    covered = if (valid) (abs(est - true_val) <= half) else NA,
    singular = singular,
    diagnostic = as.character(diagnostic),
    stringsAsFactors = FALSE
  )
}

# --- Regime 1: Partial presentation (5 of 20 items/person), random intercepts ---
run_regime1 <- function() {
  des <- design_crossed(
    list(group = c("control", "treatment")),
    allocation = "between_subj", items_per_subject = 5L, item_seed = 417L
  )
  sched <- generate_design_schedule(des, n_subj = 64L, n_item = 20L)
  conds <- list(
    null    = list(p2 = 0.60, base_seed = 30000L),
    planned = list(p2 = 0.55, base_seed = 31000L)
  )
  out <- lapply(names(conds), function(cnm) {
    cfg <- conds[[cnm]]
    dgp <- calibrate_dgp(
      des, c(control = 0.60, treatment = cfg$p2), chance = 0,
      subj_sd = 0.5, item_sd = 0.3, slope_sd_ratio = 0, corr_structure = 0
    )
    sp_bt <- specify_target(
      dgp, "group", scale = "pp", engine = "bayes_laplace",
      analysis_re = "intercept_only", prob_se_method = "bootstrap",
      control = target_control(prob_estimand = "marginal", prob_boot = n_boot, on_singular = "flag_only")
    )
    sp_jd <- specify_target(
      dgp, "group", scale = "pp", engine = "bayes_laplace",
      analysis_re = "intercept_only", prob_se_method = "joint_delta",
      control = target_control(prob_estimand = "marginal", on_singular = "flag_only")
    )
    target_pp <- sp_bt$true_effect_target
    reps <- parallel::mclapply(seq_len(n_rep), function(r) {
      s_gen <- cfg$base_seed + r
      set.seed(s_gen, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion", sample.kind = "Rejection")
      d <- simulate_dgp_data(dgp, schedule = sched)

      set.seed(s_gen + 500000L, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion", sample.kind = "Rejection")
      f_bt <- glmmplanr:::fit_single_dataset(d, sp_bt, engine_use = "bayes_laplace", re_structure = "intercept_only")
      f_jd <- glmmplanr:::fit_single_dataset(d, sp_jd, engine_use = "bayes_laplace", re_structure = "intercept_only")

      set.seed(s_gen + 600000L, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion", sample.kind = "Rejection")
      f_gl <- glmer_marginal_boot(d, refits = n_boot)

      rbind(
        make_fit_row("partial_intercept_64x20_k5", cnm, "bayes_laplace_boot99", "pp", r, s_gen, nrow(d),
                     target_pp, if (isTRUE(f_bt$converged)) 100 * f_bt$prob_est else NA_real_,
                     if (isTRUE(f_bt$converged)) 100 * f_bt$prob_se else NA_real_, f_bt$singular, f_bt$error_reason),
        make_fit_row("partial_intercept_64x20_k5", cnm, "bayes_laplace_joint_delta", "pp", r, s_gen, nrow(d),
                     target_pp, if (isTRUE(f_jd$converged)) 100 * f_jd$prob_est else NA_real_,
                     if (isTRUE(f_jd$converged)) 100 * f_jd$prob_se else NA_real_, f_jd$singular, f_jd$error_reason),
        make_fit_row("partial_intercept_64x20_k5", cnm, "glmer_boot99", "pp", r, s_gen, nrow(d),
                     target_pp, 100 * f_gl$estimate, 100 * f_gl$se, f_gl$singular, f_gl$diagnostic)
      )
    }, mc.cores = workers, mc.set.seed = FALSE)
    do.call(rbind, reps)
  })
  do.call(rbind, out)
}

# --- Regime 2: Partial presentation (8 of 24 items/person) with correlated random slopes ---
run_regime2 <- function() {
  des <- design_crossed(
    list(condition = c("A", "B")),
    allocation = "latin_square", items_per_subject = 8L, item_seed = 417L
  )
  sched <- generate_design_schedule(des, n_subj = 64L, n_item = 24L)
  conds <- list(
    null    = list(p2 = 0.60, base_seed = 40000L),
    planned = list(p2 = 0.70, base_seed = 41000L)
  )
  out <- lapply(names(conds), function(cnm) {
    cfg <- conds[[cnm]]
    dgp <- calibrate_dgp(
      des, c(A = 0.60, B = cfg$p2), chance = 0,
      subj_sd = 0.8, item_sd = 0.6, slope_sd_ratio = 0.30, corr_structure = 0.20
    )
    sp_bt <- specify_target(
      dgp, "condition", scale = "pp", engine = "bayes_laplace",
      analysis_re = "maximal", prob_se_method = "bootstrap",
      control = target_control(prob_estimand = "marginal", prob_boot = n_boot, on_singular = "flag_only")
    )
    sp_jd <- specify_target(
      dgp, "condition", scale = "pp", engine = "bayes_laplace",
      analysis_re = "maximal", prob_se_method = "joint_delta",
      control = target_control(prob_estimand = "marginal", on_singular = "flag_only")
    )
    reps <- parallel::mclapply(seq_len(n_rep), function(r) {
      s_gen <- cfg$base_seed + r
      set.seed(s_gen, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion", sample.kind = "Rejection")
      d <- simulate_dgp_data(dgp, schedule = sched)

      set.seed(s_gen + 500000L, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion", sample.kind = "Rejection")
      f_bt <- glmmplanr:::fit_single_dataset(d, sp_bt, engine_use = "bayes_laplace", re_structure = "maximal")
      f_jd <- glmmplanr:::fit_single_dataset(d, sp_jd, engine_use = "bayes_laplace", re_structure = "maximal")

      rbind(
        make_fit_row("partial_slopes_64x24_k8", cnm, "bayes_laplace_logit", "logit", r, s_gen, nrow(d),
                     sp_jd$true_beta, if (isTRUE(f_jd$converged)) f_jd$est else NA_real_,
                     if (isTRUE(f_jd$converged)) f_jd$se else NA_real_, f_jd$singular, f_jd$error_reason),
        make_fit_row("partial_slopes_64x24_k8", cnm, "bayes_laplace_joint_delta", "pp", r, s_gen, nrow(d),
                     sp_jd$true_effect_target, if (isTRUE(f_jd$converged)) 100 * f_jd$prob_est else NA_real_,
                     if (isTRUE(f_jd$converged)) 100 * f_jd$prob_se else NA_real_, f_jd$singular, f_jd$error_reason),
        make_fit_row("partial_slopes_64x24_k8", cnm, "bayes_laplace_boot99", "pp", r, s_gen, nrow(d),
                     sp_bt$true_effect_target, if (isTRUE(f_bt$converged)) 100 * f_bt$prob_est else NA_real_,
                     if (isTRUE(f_bt$converged)) 100 * f_bt$prob_se else NA_real_, f_bt$singular, f_bt$error_reason)
      )
    }, mc.cores = workers, mc.set.seed = FALSE)
    do.call(rbind, reps)
  })
  do.call(rbind, out)
}

# --- Regime 3: Testlet nesting (48 subj x 6 passages x 4 items) with/without local dependence ---
run_regime3 <- function() {
  des <- design_testlet(
    list(condition = c("A", "B")),
    testlet_size = 4L, between_item = "condition", contrast = "treatment"
  )
  sched <- generate_design_schedule(des, n_subj = 48L, n_item = 24L)
  configs <- list(
    list(reg = "testlet_matched_48x24_pt0", pt_sd = 0.0, cnm = "null",    p2 = 0.60, base_seed = 50000L),
    list(reg = "testlet_matched_48x24_pt0", pt_sd = 0.0, cnm = "planned", p2 = 0.70, base_seed = 51000L),
    list(reg = "testlet_localdep_48x24_pt04", pt_sd = 0.4, cnm = "null",    p2 = 0.60, base_seed = 52000L),
    list(reg = "testlet_localdep_48x24_pt04", pt_sd = 0.4, cnm = "planned", p2 = 0.70, base_seed = 53000L)
  )
  out <- lapply(configs, function(cfg) {
    dgp <- calibrate_dgp(
      des, c(A = 0.60, B = cfg$p2), chance = 0,
      subj_sd = 0.7, item_sd = 0.4, slope_sd_ratio = 0, corr_structure = 0,
      control = calibrate_control(testlet_sd = 0.5, participant_testlet_sd = cfg$pt_sd)
    )
    sp_jd <- specify_target(
      dgp, "condition", scale = "pp", engine = "bayes_laplace",
      analysis_re = "intercept_only", prob_se_method = "joint_delta",
      control = target_control(prob_estimand = "marginal", on_singular = "flag_only")
    )
    reps <- parallel::mclapply(seq_len(n_rep), function(r) {
      s_gen <- cfg$base_seed + r
      set.seed(s_gen, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion", sample.kind = "Rejection")
      d <- simulate_dgp_data(dgp, schedule = sched)

      f_jd <- glmmplanr:::fit_single_dataset(d, sp_jd, engine_use = "bayes_laplace", re_structure = "intercept_only")

      gl_sing <- NA
      gl_diag <- ""
      gl_coef <- tryCatch({
        m <- suppressMessages(glmer(
          y ~ condition + (1 | subj) + (1 | item) + (1 | testlet),
          data = d, family = binomial(),
          control = glmerControl(optimizer = "bobyqa")
        ))
        gl_sing <- isSingular(m)
        conv <- m@optinfo$conv
        msg <- conv$lme4$messages
        if (any(conv$opt != 0L) || any(!grepl("boundary (singular) fit", msg, fixed = TRUE))) {
          stop("glmer convergence warning")
        }
        cf <- summary(m)$coefficients
        idx <- grep("^condition", rownames(cf))[1L]
        c(est = unname(cf[idx, "Estimate"]), se = unname(cf[idx, "Std. Error"]))
      }, error = function(e) {
        gl_diag <<- conditionMessage(e)
        c(est = NA_real_, se = NA_real_)
      })

      rbind(
        make_fit_row(cfg$reg, cfg$cnm, "bayes_laplace_joint_delta", "pp", r, s_gen, nrow(d),
                     sp_jd$true_effect_target, if (isTRUE(f_jd$converged)) 100 * f_jd$prob_est else NA_real_,
                     if (isTRUE(f_jd$converged)) 100 * f_jd$prob_se else NA_real_, f_jd$singular, f_jd$error_reason),
        make_fit_row(cfg$reg, cfg$cnm, "bayes_laplace_logit", "logit", r, s_gen, nrow(d),
                     sp_jd$true_beta, if (isTRUE(f_jd$converged)) f_jd$est else NA_real_,
                     if (isTRUE(f_jd$converged)) f_jd$se else NA_real_, f_jd$singular, f_jd$error_reason),
        make_fit_row(cfg$reg, cfg$cnm, "glmer_logit", "logit", r, s_gen, nrow(d),
                     sp_jd$true_beta, gl_coef["est"], gl_coef["se"], gl_sing, gl_diag)
      )
    }, mc.cores = workers, mc.set.seed = FALSE)
    do.call(rbind, reps)
  })
  do.call(rbind, out)
}

fits <- rbind(run_regime1(), run_regime2(), run_regime3())
stopifnot(!anyDuplicated(fits[c("regime", "condition", "estimator", "scale", "replication")]))

keys <- c("regime", "condition", "estimator", "scale")
groups <- split(fits, fits[keys], drop = TRUE)
summary_df <- do.call(rbind, lapply(groups, function(d) {
  ok <- d$valid
  n_att <- nrow(d)
  n_val <- sum(ok)
  k_det <- sum(d$detected)
  k_cov <- sum(d$covered[ok])
  det_ci <- as.numeric(stats::binom.test(k_det, n_att)$conf.int)
  cov_ci <- if (n_val > 0L) as.numeric(stats::binom.test(k_cov, n_val)$conf.int) else c(NA_real_, NA_real_)
  status <- if (n_val < ceiling(0.9 * n_att)) {
    "insufficient_valid"
  } else if (cov_ci[1L] >= 0.90) {
    "pass"
  } else if (cov_ci[2L] < 0.90) {
    "below_floor"
  } else {
    "indeterminate"
  }
  sing_cnt <- if (all(is.na(d$singular[ok]))) NA_integer_ else sum(d$singular[ok], na.rm = TRUE)
  data.frame(
    regime = d$regime[1L],
    condition = d$condition[1L],
    estimator = d$estimator[1L],
    scale = d$scale[1L],
    true_target = d$true_target[1L],
    trials = d$trials[1L],
    attempted = n_att,
    valid = n_val,
    detected = k_det,
    detection = k_det / n_att,
    detection_lo = det_ci[1L],
    detection_hi = det_ci[2L],
    covered = k_cov,
    coverage = if (n_val > 0L) k_cov / n_val else NA_real_,
    coverage_lo = cov_ci[1L],
    coverage_hi = cov_ci[2L],
    emp_sd = if (n_val > 1L) stats::sd(d$estimate[ok]) else NA_real_,
    mean_se = if (n_val > 0L) mean(d$se[ok]) else NA_real_,
    mean_width = if (n_val > 0L) mean(d$width[ok]) else NA_real_,
    singular_valid = sing_cnt,
    floor_status = status,
    R_version = as.character(getRversion()),
    glmmplanr_version = as.character(packageVersion("glmmplanr")),
    lme4_version = as.character(packageVersion("lme4")),
    stringsAsFactors = FALSE
  )
}))
summary_df <- summary_df[order(summary_df$regime, summary_df$condition, summary_df$estimator), ]
rownames(summary_df) <- NULL

write.csv(fits, file.path(out_dir, "sparse_and_testlet_recheck_fits.csv"), row.names = FALSE)
write.csv(summary_df, file.path(out_dir, "sparse_and_testlet_recheck_summary.csv"), row.names = FALSE)
print(
  summary_df[c("regime", "condition", "estimator", "scale", "valid", "detected", "covered", "coverage_lo", "floor_status")],
  row.names = FALSE, digits = 3
)
