test_that("evaluate_stan_hmc_diagnostics refuses to pass without diagnostic backend and checks treedepth & all monitored draws", {
  set.seed(123)
  good_mat <- matrix(rnorm(800), nrow = 400, ncol = 2)
  bad_rhat <- cbind(rnorm(400, -3, 0.2), rnorm(400, 3, 0.2))

  # Missing both posterior and rstan must fail explicitly rather than faking Rhat=1, ESS=400
  no_pkg <- glmmplanr:::evaluate_stan_hmc_diagnostics(
    list("beta[2]" = good_mat), n_div = 0L, n_max_td = 0L,
    has_posterior = FALSE, has_rstan = FALSE
  )
  expect_false(no_pkg$ok)
  expect_identical(no_pkg$status, "engine_error")
  expect_match(no_pkg$reason, "posterior")

  if (requireNamespace("posterior", quietly = TRUE) || requireNamespace("rstan", quietly = TRUE)) {
    ok_res <- glmmplanr:::evaluate_stan_hmc_diagnostics(
      list("beta[2]" = good_mat, prob_contrast = good_mat, "tau_subj[1]" = abs(good_mat) + 0.5),
      n_div = 0L, n_max_td = 0L
    )
    expect_true(ok_res$ok)
    expect_identical(ok_res$status, "converged")

    # Exceeding max_treedepth must fail convergence gate
    td_fail <- glmmplanr:::evaluate_stan_hmc_diagnostics(
      list("beta[2]" = good_mat), n_div = 0L, n_max_td = 3L
    )
    expect_false(td_fail$ok)
    expect_identical(td_fail$status, "mcmc_not_converged")
    expect_match(td_fail$reason, "max_td=3")

    # Non-convergence in prob_contrast or variance component must fail gate
    contrast_fail <- glmmplanr:::evaluate_stan_hmc_diagnostics(
      list("beta[2]" = good_mat, prob_contrast = bad_rhat), n_div = 0L, n_max_td = 0L
    )
    expect_false(contrast_fail$ok)
    expect_match(contrast_fail$reason, "prob_contrast")
  }
})

test_that("stan_laplace emits deprecation warning and Bonferroni MC adjustment widens Stage-2 Wilson bounds", {
  des <- design_crossed(list(A = c("a", "b")), allocation = "latin_square")
  dgp <- calibrate_dgp(des, c(a = 0.45, b = 0.70), subj_sd = 0.4, item_sd = 0.3)

  if (glmmplanr:::has_cmdstan_backend() || requireNamespace("rstan", quietly = TRUE)) {
    expect_warning(
      sp_sl <- specify_target(dgp, "A", engine = "stan_laplace"),
      "deprecated alias"
    )
  } else {
    expect_error(
      suppressWarnings(specify_target(dgp, "A", engine = "stan_laplace")),
      "Stan backend"
    )
  }

  sp <- specify_target(dgp, "A", criterion = "power", power_min = 0.80,
                       engine = "bayes_laplace", analysis_re = "intercept_only",
                       control = target_control(on_singular = "flag_only"))
  res_none <- optimize_sample_size(
    sp, n_subj_range = c(12L, 20L), n_item_range = c(8L, 12L),
    method = "grid", nsim_final = 16L, workers = 1L, seed = 777L, verbose = FALSE,
    control = optimize_control(nsim_sensitivity = 0L, mc_adjustment = "none")
  )
  res_bonf <- optimize_sample_size(
    sp, n_subj_range = c(12L, 20L), n_item_range = c(8L, 12L),
    method = "grid", nsim_final = 16L, workers = 1L, seed = 777L, verbose = FALSE,
    control = optimize_control(nsim_sensitivity = 0L, mc_adjustment = "bonferroni")
  )

  expect_identical(res_bonf$search_meta$mc_adjustment, "bonferroni")
  expect_gt(res_bonf$search_meta$mc_conf_level_effective, sp$mc_conf_level)
  expect_true(all(res_bonf$verified$power_mc_lo <= res_none$verified$power_mc_lo + 1e-12))
  expect_true(all(res_bonf$verified$power_mc_hi >= res_none$verified$power_mc_hi - 1e-12))
  expect_match(report_plan(res_bonf, print_output = FALSE), "Bonferroni")
  expect_match(report_plan(res_none, print_output = FALSE), "pointwise")
})

test_that("target attainment and inferential validity are separated and low coverage triggers caution", {
  des <- design_crossed(list(A = c("a", "b")), allocation = "latin_square")
  dgp <- calibrate_dgp(des, c(a = 0.45, b = 0.70), subj_sd = 0.4, item_sd = 0.3)
  sp <- specify_target(dgp, "A", criterion = "power", power_min = 0.80,
                       engine = "bayes_laplace", analysis_re = "intercept_only",
                       control = target_control(on_singular = "flag_only"))
  res <- optimize_sample_size(
    sp, n_subj_range = 16L, n_item_range = 8L,
    method = "grid", nsim_final = 10L, workers = 1L, seed = 778L, verbose = FALSE,
    control = optimize_control(nsim_sensitivity = 0L)
  )

  rep_normal <- report_plan(res, print_output = FALSE)
  expect_match(rep_normal, "does not verify inferential validity")
  expect_match(rep_normal, "validate_design_point")

  # Simulate anti-conservative SEs yielding low empirical coverage (60% over 100 valid fits)
  res_low <- res
  res_low$optimal$n_valid <- 100L
  res_low$optimal$coverage <- 0.60
  res_low$optimal$prob_coverage <- 0.60
  res_low$verified$n_valid <- 100L
  res_low$verified$coverage <- 0.60
  res_low$verified$prob_coverage <- 0.60

  val_low <- glmmplanr:::assess_inferential_validity(0.60, 100L, 0.05)
  expect_identical(val_low$status, "caution_low_coverage")
  expect_identical(val_low$short_label, "Caution: Low Cov.")

  rep_low <- report_plan(res_low, print_output = FALSE)
  expect_match(rep_low, "Caution:")
  expect_match(rep_low, "underestimated standard errors")

  printed <- paste(testthat::capture_messages(print(res_low)), collapse = "\n")
  expect_match(printed, "Target attainment:")
  expect_match(printed, "Inferential validity:")
  expect_match(printed, "CAUTION")

  gt_tbl <- as_gt(res_low)
  expect_true("validity_status" %in% names(gt_tbl$`_data`))
  expect_true(all(gt_tbl$`_data`$validity_status == "Caution: Low Cov."))
})

test_that("independent top-up applies sequential multiplicity control across rounds", {
  des <- design_crossed(list(A = c("a", "b")), allocation = "latin_square")
  dgp <- calibrate_dgp(des, c(a = 0.45, b = 0.70), subj_sd = 0.4, item_sd = 0.3)
  sp <- specify_target(dgp, "A", criterion = "power", power_min = 0.80,
                       engine = "bayes_laplace", analysis_re = "intercept_only",
                       control = target_control(on_singular = "flag_only"))

  # Automatic Stage-3 independent top-up uses two-stage split (97.5% Wilson for 95% nominal)
  res_2stage <- optimize_sample_size(
    sp, n_subj_range = 12L, n_item_range = 8L,
    method = "grid", nsim_final = 8L, workers = 1L, seed = 779L, verbose = FALSE,
    control = optimize_control(nsim_topup = 8L, topup_policy = "independent", nsim_sensitivity = 0L)
  )
  expect_true(isTRUE(res_2stage$search_meta$sequential_topup_adjusted))
  expect_equal(res_2stage$search_meta$mc_conf_level_effective, 0.975)

  # Interactive topup_verification applies geometric alpha-spending and warns from round 2
  res_base <- optimize_sample_size(
    sp, n_subj_range = 12L, n_item_range = 8L,
    method = "grid", nsim_final = 8L, workers = 1L, seed = 780L, verbose = FALSE,
    control = optimize_control(nsim_topup = 0L, nsim_sensitivity = 0L)
  )
  res_base$verified$decision_status <- "indeterminate"
  r1 <- topup_verification(res_base, nsim_extra = 8L, verbose = FALSE)
  expect_equal(r1$search_meta$topup_rounds, 1L)
  expect_equal(r1$search_meta$topup_mc_conf_level, 1 - 0.05 / 4)

  r1$verified$decision_status <- "indeterminate"
  expect_warning(
    r2 <- topup_verification(r1, nsim_extra = 8L, verbose = FALSE),
    "round 2"
  )
  expect_equal(r2$search_meta$topup_rounds, 2L)
  expect_equal(r2$search_meta$topup_mc_conf_level, 1 - 0.05 / 8)
  expect_match(report_plan(r2, print_output = FALSE), "geometric sequential alpha-spending")

  # Unadjusted repeated top-up explicitly warns
  expect_warning(
    topup_verification(res_base, nsim_extra = 8L, sequential_adjustment = "none", verbose = FALSE),
    "inflates overall false-pass"
  )

  # Exact binomial check: geometric alpha-spending bounds 10-round cumulative false-pass risk
  p_true <- 0.799
  n_batch <- 1000L
  target <- 0.80
  k_vals <- 0:n_batch
  probs <- stats::dbinom(k_vals, n_batch, p_true)
  round_rates <- function(conf_vec) {
    surv <- 1.0
    cum_pass <- 0.0
    for (cf in conf_vec) {
      lo <- vapply(k_vals, function(k) glmmplanr:::wilson_mc_ci(k, n_batch, cf)[["lo"]], numeric(1L))
      hi <- vapply(k_vals, function(k) glmmplanr:::wilson_mc_ci(k, n_batch, cf)[["hi"]], numeric(1L))
      p_pass <- sum(probs[lo >= target])
      p_ind <- sum(probs[lo < target & hi >= target])
      cum_pass <- cum_pass + surv * p_pass
      surv <- surv * p_ind
    }
    cum_pass
  }
  unadj_10 <- round_rates(rep(0.95, 10L))
  geom_confs <- c(0.95, vapply(1:9, function(r) glmmplanr:::sequential_topup_conf_level(0.95, r, 1L, "geometric"), numeric(1L)))
  geom_10 <- round_rates(geom_confs)
  expect_gt(unadj_10, 0.14)
  expect_lt(geom_10, 0.035)
})

test_that("evaluate_design_point and validate_design_point accept custom data_fn for misspecification tests", {
  des <- design_crossed(list(A = c("a", "b")), allocation = "latin_square")
  dgp <- calibrate_dgp(des, c(a = 0.50, b = 0.72), subj_sd = 0.6, item_sd = 0.4, slope_sd_ratio = 0)
  sp <- specify_target(dgp, "A", scale = "logit", criterion = "power",
                       engine = "bayes_laplace", analysis_re = "intercept_only",
                       control = target_control(on_singular = "flag_only"))

  # Custom misspecified DGP: drop the lowest-performing 25% of subjects (MNAR-like dropout)
  mnar_fn <- function(schedule, spec) {
    df <- simulate_dgp_data(spec$dgp, schedule = schedule)
    means <- tapply(df$y, df$subj_id, mean)
    keep_s <- as.integer(names(means)[means >= stats::quantile(means, 0.25)])
    df[df$subj_id %in% keep_s, , drop = FALSE]
  }

  ev_custom <- evaluate_design_point(sp, 16L, 8L, nsim = 6L, stream_seeds = 901L, data_fn = mnar_fn)
  expect_equal(ev_custom$nsim, 6L)
  expect_lt(ev_custom$n_trials_obs_mean, 16L * 8L)
  expect_true(is.finite(ev_custom$power))

  val_custom <- validate_design_point(sp, 16L, 8L, nsim = 4L, seed = 902L, data_fn = mnar_fn)
  expect_equal(nrow(val_custom), 2L)
  expect_equal(val_custom$condition, c("null", "planned"))

  # Invalid data_fn output produces clear errors
  expect_error(
    evaluate_design_point(sp, 16L, 8L, nsim = 1L, data_fn = function(schedule, spec) list(a = 1)),
    "data.frame"
  )
  expect_error(
    evaluate_design_point(sp, 16L, 8L, nsim = 1L, data_fn = function(schedule, spec) schedule$trials),
    "binary"
  )
})

test_that("regularized singularity NA, single-cell validation, uneven-grid plot, and boundary/exploratory notes work", {
  # 1. Single-condition absolute probability target validates planned precision without aborting
  des_one <- design_crossed(list())
  dgp_one <- calibrate_dgp(des_one, cell_probs = 0.65, chance = 0, subj_sd = 0.6, item_sd = 0.5)
  sp_one <- specify_target(
    dgp_one, effect = "(Intercept)", scale = "pp", criterion = "ci_width",
    ci_width = 25, prob_se_method = "joint_delta", engine = "bayes_laplace"
  )
  expect_message(
    val_one <- validate_design_point(sp_one, n_subj = 24L, n_item = 12L, nsim = 5L, seed = 301L),
    "absolute probability"
  )
  expect_equal(nrow(val_one), 1L)
  expect_identical(val_one$condition, "planned")
  expect_true(is.finite(val_one$coverage))
  expect_true(is.finite(val_one$width_attainment))

  # 2. Uneven grid contour plot emits zero ggplot2 raster warnings; regularized singular_rate is NA
  des_lat <- design_crossed(list(learning = c("incidental", "instruction")),
                            allocation = "latin_square", items_per_subject = 8L, item_seed = 11L)
  dgp_lat <- calibrate_dgp(des_lat, c(incidental = 0.55, instruction = 0.72),
                           chance = 1 / 3, subj_sd = 0.7, item_sd = 0.5, slope_sd_ratio = 0)
  dgp_out <- paste(cli::ansi_strip(testthat::capture_messages(print(dgp_lat))), collapse = "\n")
  expect_match(dgp_out, "Implied 95% accuracy spans")
  expect_match(dgp_out, "Random-slope assumption: 0")

  sp_lat <- specify_target(
    dgp_lat, effect = "learning", scale = "pp", criterion = "power", power_min = 0.80,
    prob_se_method = "joint_delta", engine = "bayes_laplace", analysis_re = "intercept_only"
  )
  plan_lat <- optimize_sample_size(
    sp_lat, n_subj_range = c(24, 32, 48), n_item_range = c(12, 16, 24),
    method = "grid", nsim_final = 18L, workers = 1L, seed = 417L, verbose = FALSE,
    control = optimize_control(nsim_sensitivity = 0L)
  )
  expect_true(all(is.na(plan_lat$verified$singular_rate)))

  prt_msg <- paste(cli::ansi_strip(capture.output(print(plan_lat))), collapse = "\n")
  if (!nzchar(prt_msg)) {
    prt_msg <- paste(cli::ansi_strip(testthat::capture_messages(print(plan_lat))), collapse = "\n")
  }
  expect_match(prt_msg, "not applicable \\(regularized\\)")
  expect_match(prt_msg, "EXPLORATORY REGIME")

  rep_txt <- report_plan(plan_lat, print_output = FALSE)
  expect_match(rep_txt, "not applicable under regularized estimation")
  expect_match(rep_txt, "exploratory design regime")

  # Force a near-boundary pass (+1.0 pp margin) with a higher-margin alternative to verify boundary note
  plan_bnd <- plan_lat
  plan_bnd$feasible <- TRUE
  plan_bnd$verified$decision_status <- "fail"
  plan_bnd$verified$decision_status[1L] <- "pass"
  plan_bnd$verified$power[1L] <- 0.88
  plan_bnd$verified$power_mc_lo[1L] <- 0.81
  plan_bnd$verified$power_mc_hi[1L] <- 0.93
  plan_bnd$verified$decision_status[nrow(plan_bnd$verified)] <- "pass"
  plan_bnd$verified$power[nrow(plan_bnd$verified)] <- 0.98
  plan_bnd$verified$power_mc_lo[nrow(plan_bnd$verified)] <- 0.92
  plan_bnd$verified$power_mc_hi[nrow(plan_bnd$verified)] <- 0.995
  plan_bnd$optimal <- plan_bnd$verified[1L, , drop = FALSE]
  plan_bnd$decision_status <- "pass"
  rep_bnd <- report_plan(plan_bnd, print_output = FALSE)
  expect_match(rep_bnd, "Boundary sensitivity:")
  expect_match(rep_bnd, "higher-margin verified alternative")

  # Build ggplot on uneven grid and ensure zero warnings
  plot_warns <- character(0)
  withCallingHandlers({
    gp <- plot(plan_lat, type = "contour")
    built <- ggplot2::ggplot_build(gp)
  }, warning = function(w) {
    plot_warns <<- c(plot_warns, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_length(plot_warns, 0L)
  expect_s3_class(built, "ggplot_built")
})
