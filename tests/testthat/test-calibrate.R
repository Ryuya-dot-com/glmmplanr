test_that("2D Gauss-Hermite marginal inversion recovers m-AFC cell probabilities accurately", {
  # Verify substantive SD conversion round-trip
  sd_s <- prob_range_to_sd(c(0.40, 0.85), chance = 1 / 3)
  expect_gt(sd_s, 0.8)
  expect_lt(sd_s, 0.95)
  rng_back <- sd_to_prob_range(sd_s, intercept = mean(mafc_qlogis(c(0.40, 0.85), 1 / 3)), chance = 1 / 3)
  expect_equal(rng_back, c(0.40, 0.85), tolerance = 1e-10)

  des <- design_crossed(
    factors = list(
      exposure = c("incidental", "intentional"),
      spacing  = c("massed", "distributed")
    ),
    allocation = "latin_square",
    contrast   = "sum",
    attrition  = 0.05
  )

  target_p <- c(
    incidental_massed       = 0.50,
    incidental_distributed  = 0.60,
    intentional_massed      = 0.65,
    intentional_distributed = 0.80
  )

  dgp <- calibrate_dgp(
    des,
    cell_probs      = target_p,
    chance          = 1 / 3,
    subj_prob_range = c(0.40, 0.85),
    item_prob_range = c(0.45, 0.88),
    slope_sd_ratio  = 0.35,
    corr_structure  = 0.20,
    check_nsim      = 8000L,
    seed            = 123L
  )

  expect_s3_class(dgp, "glmmplanr_dgp")
  expect_lt(dgp$solver_info$max_abs_residual, 1e-11)
  expect_equal(unname(dgp$achieved_cell_probs), unname(target_p), tolerance = 1e-10)

  # Non-collapsibility check: calibrated main effects must exceed naive unattenuated m-AFC logit
  expect_gt(abs(dgp$beta[["exposure1"]]), abs(dgp$beta_naive_mafc[["exposure1"]]))
  expect_gt(abs(dgp$beta[["spacing1"]]), abs(dgp$beta_naive_mafc[["spacing1"]]))

  # Empirical Monte Carlo marginal mean check matches target probabilities closely
  expect_equal(dgp$prior_check$mc_mean_prob, unname(target_p), tolerance = 0.012)

  # Sensitivity DGPs preserve exact target marginal probabilities
  sens <- sensitivity_dgps(dgp)
  expect_named(sens, c("optimistic", "standard", "conservative"))
  expect_equal(unname(sens$conservative$achieved_cell_probs), unname(target_p), tolerance = 1e-10)

  chk <- check_dgp(dgp, plot = TRUE, n_draws = 500L)
  expect_s3_class(chk$plot, "ggplot")
  conditional <- 1/3 + (2/3) * plogis(drop(des$X_cells %*% dgp$beta))
  expect_equal(chk$cell_diagnostics$conditional_prob, unname(conditional))
  expect_equal(chk$cell_diagnostics$quad_prob, unname(target_p), tolerance = 1e-10)
})
