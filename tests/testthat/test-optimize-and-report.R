test_that("2D GP surrogate optimizer is worker-invariant, monotone, and generates tables/plots/reports", {
  des <- design_crossed(
    factors = list(
      exposure = c("incidental", "intentional"),
      spacing  = c("massed", "distributed")
    ),
    allocation = "latin_square",
    contrast   = "sum",
    attrition  = 0.05
  )

  dgp <- calibrate_dgp(
    des,
    cell_probs = c(
      incidental_massed       = 0.50,
      incidental_distributed  = 0.58,
      intentional_massed      = 0.62,
      intentional_distributed = 0.85
    ),
    chance          = 1 / 3,
    subj_prob_range = c(0.42, 0.84),
    item_prob_range = c(0.45, 0.86),
    slope_sd_ratio  = 0.30,
    corr_structure  = 0.20
  )

  on_cran <- !identical(Sys.getenv("NOT_CRAN"), "true")
  spec <- specify_target(
    dgp,
    effect      = "exposure:spacing",
    criterion   = "both",
    power_min   = 0.75,
    ci_width    = 0.95,
    engine      = "glmmTMB",
    on_singular = if (on_cran) "flag_only" else "flag_and_compare"
  )

  subj_rng <- if (on_cran) c(16, 32) else c(24, 64)
  item_rng <- if (on_cran) c(12, 24) else c(16, 40)
  n_init   <- if (on_cran) 4L else 8L
  n_adapt  <- if (on_cran) 2L else 4L
  ns_init  <- if (on_cran) 4L else 15L
  ns_final <- if (on_cran) 6L else 25L
  n_ver    <- if (on_cran) 2L else 3L
  ns_sens  <- if (on_cran) 4L else 15L

  res_w1 <- optimize_sample_size(
    spec,
    n_subj_range     = subj_rng,
    n_item_range     = item_rng,
    prefer           = "balanced",
    method           = "surrogate_gp",
    n_eval_init      = n_init,
    n_eval_adapt     = n_adapt,
    nsim_init        = ns_init,
    nsim_final       = ns_final,
    n_verify         = n_ver,
    nsim_sensitivity = ns_sens,
    workers          = 1L,
    seed             = 20261119L,
    verbose          = FALSE
  )

  caller_rng <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
  res_w2 <- optimize_sample_size(
    spec,
    n_subj_range     = subj_rng,
    n_item_range     = item_rng,
    prefer           = "balanced",
    method           = "surrogate_gp",
    n_eval_init      = n_init,
    n_eval_adapt     = n_adapt,
    nsim_init        = ns_init,
    nsim_final       = ns_final,
    n_verify         = n_ver,
    nsim_sensitivity = ns_sens,
    workers          = 2L,
    seed             = 20261119L,
    verbose          = FALSE
  )

  # The RNG state and design decision are reproducible across worker counts.
  expect_identical(get0(".Random.seed", envir = globalenv(), inherits = FALSE), caller_rng)
  expect_equal(res_w1$optimal$n_subj, res_w2$optimal$n_subj)
  expect_equal(res_w1$optimal$n_item, res_w2$optimal$n_item)
  expect_identical(res_w1$optimal$decision_status, res_w2$optimal$decision_status)
  expect_lt(abs(res_w1$optimal$power - res_w2$optimal$power), 1e-8)
  expect_lt(max(abs(res_w1$surface$pred_power - res_w2$surface$pred_power)), 1e-8)

  # Multiples of 4 strictly enforced
  expect_equal(res_w1$optimal$n_subj %% 4L, 0L)
  expect_equal(res_w1$optimal$n_item %% 4L, 0L)

  # Verify 2D monotonicity of predicted GP surface
  surf <- res_w1$surface
  for (ss in unique(surf$n_subj)) {
    sub_s <- surf[surf$n_subj == ss, ]
    sub_s <- sub_s[order(sub_s$n_item), ]
    expect_true(all(diff(sub_s$pred_power) >= -1e-9))
    expect_true(all(diff(sub_s$pred_ci_width) <= 1e-9))
  }
  for (ii in unique(surf$n_item)) {
    sub_i <- surf[surf$n_item == ii, ]
    sub_i <- sub_i[order(sub_i$n_subj), ]
    expect_true(all(diff(sub_i$pred_power) >= -1e-9))
    expect_true(all(diff(sub_i$pred_ci_width) <= 1e-9))
  }

  # Verify S3 outputs: gt table, 2D contour, Pareto curve, 3D plotly, and the sample size report
  gt_tbl <- as_gt(res_w1)
  expect_s3_class(gt_tbl, "gt_tbl")

  p_cont <- plot(res_w1, type = "contour")
  expect_s3_class(p_cont, "ggplot")

  p_par <- plot(res_w1, type = "pareto")
  expect_s3_class(p_par, "ggplot")

  if (requireNamespace("plotly", quietly = TRUE)) {
    p_3d <- plot(res_w1, type = "surface3d")
    expect_s3_class(p_3d, "plotly")
  }

  report <- report_plan(res_w1, print_output = FALSE)
  expect_match(report, "Sample Size Justification")
  expect_output(printed <- withVisible(report_plan(res_w1)), "Sample Size Justification")
  expect_false(printed$visible)
  expect_identical(printed$value, report)
  expect_error(report_plan(res_w1, "en"), "print_output")
})
