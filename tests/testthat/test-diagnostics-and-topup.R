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

  expect_warning(
    sp_sl <- specify_target(dgp, "A", engine = "stan_laplace"),
    "deprecated alias"
  )

  sp <- specify_target(dgp, "A", criterion = "power", power_min = 0.80,
                       engine = "bayes_laplace", analysis_re = "intercept_only",
                       on_singular = "flag_only")
  res_none <- optimize_sample_size(
    sp, n_subj_range = c(12L, 20L), n_item_range = c(8L, 12L),
    method = "grid", nsim_final = 16L, nsim_sensitivity = 0L,
    mc_adjustment = "none", workers = 1L, seed = 777L, verbose = FALSE
  )
  res_bonf <- optimize_sample_size(
    sp, n_subj_range = c(12L, 20L), n_item_range = c(8L, 12L),
    method = "grid", nsim_final = 16L, nsim_sensitivity = 0L,
    mc_adjustment = "bonferroni", workers = 1L, seed = 777L, verbose = FALSE
  )

  expect_identical(res_bonf$search_meta$mc_adjustment, "bonferroni")
  expect_gt(res_bonf$search_meta$mc_conf_level_effective, sp$mc_conf_level)
  expect_true(all(res_bonf$verified$power_mc_lo <= res_none$verified$power_mc_lo + 1e-12))
  expect_true(all(res_bonf$verified$power_mc_hi >= res_none$verified$power_mc_hi - 1e-12))
  expect_match(report_plan(res_bonf, print_output = FALSE), "Bonferroni")
  expect_match(report_plan(res_none, print_output = FALSE), "pointwise")
})
