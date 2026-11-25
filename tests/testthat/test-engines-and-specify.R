test_that("C++ Laplace, Regularized Bayes Laplace, and lme4 engines run and produce valid 3-layer metrics", {
  des <- design_crossed(
    factors = list(
      exposure = c("incidental", "intentional"),
      spacing  = c("massed", "distributed")
    ),
    allocation = "latin_square",
    contrast   = "sum"
  )

  dgp <- calibrate_dgp(
    des,
    cell_probs = c(
      incidental_massed       = 0.40,
      incidental_distributed  = 0.55,
      intentional_massed      = 0.55,
      intentional_distributed = 0.85
    ),
    chance          = 0,
    subj_prob_range = c(0.30, 0.80),
    item_prob_range = c(0.35, 0.80),
    slope_sd_ratio  = 0.30,
    corr_structure  = 0.15
  )

  # Permutation-invariant interaction lookup
  sp1 <- specify_target(dgp, effect = "exposure:spacing", criterion = "both", power_min = 0.80, ci_width = 0.80, engine = "cpp_laplace")
  sp2 <- specify_target(dgp, effect = "spacing:exposure", criterion = "both", power_min = 0.80, ci_width = 0.80, engine = "bayes_laplace")
  expect_equal(sp1$effect_col, sp2$effect_col)

  on_cran <- !identical(Sys.getenv("NOT_CRAN"), "true")
  ns_ev <- if (on_cran) 32L else 48L
  ni_ev <- if (on_cran) 20L else 32L
  n_rep <- if (on_cran) 5L else 12L

  ev_ml <- evaluate_design_point(sp1, n_subj = ns_ev, n_item = ni_ev, nsim = n_rep, stream_seeds = 99L)
  ev_by <- evaluate_design_point(sp2, n_subj = ns_ev, n_item = ni_ev, nsim = n_rep, stream_seeds = 99L)

  expect_equal(ev_ml$error_rate, 0)
  expect_equal(ev_ml$n_valid + as.integer(round(ev_ml$max_iter_rate * n_rep)), n_rep)
  expect_gt(ev_ml$n_valid, 0L)
  expect_equal(ev_by$n_valid, n_rep)
  # Regularized Bayesian engine marks boundary singularity check as not applicable (NA)
  expect_true(is.na(ev_by$singular_rate))
  expect_true(is.finite(ev_ml$power) && ev_ml$power >= 0 && ev_ml$power <= 1)
  expect_true(is.finite(ev_ml$ci_width_mean) && ev_ml$ci_width_mean > 0)
  expect_true("reduced_power" %in% names(ev_ml))

  # Verify lme4_fast backend integration
  sp_lme4 <- specify_target(dgp, effect = "exposure", criterion = "power", engine = "lme4_fast",
                            analysis_re = "uncorrelated",
                            control = target_control(on_singular = "flag_only"))
  ev_lme4 <- evaluate_design_point(sp_lme4, n_subj = 32, n_item = 20, nsim = 3L, stream_seeds = 7L)
  expect_equal(ev_lme4$n_valid, 3L)

  # Stress test models under maximal DGP
  st <- stress_test_models(
    dgp, effect = "exposure:spacing",
    n_subj = if (on_cran) 24L else 32L,
    n_item = if (on_cran) 16L else 20L,
    nsim = if (on_cran) 4L else 8L,
    include_null = TRUE, seed = 42L
  )
  expect_s3_class(st, "glmmplanr_stress_test")
  expect_equal(nrow(st), 4L)
  expect_s3_class(as_gt(st), "gt_tbl")
})

test_that("the regularized engine distinguishes small prior scales", {
  des <- design_crossed(list(arm = c("a", "b")), allocation = "between_subj",
                        contrast = "treatment")
  dgp <- calibrate_dgp(des, c(.72, .80), subj_sd = .8, item_sd = .6, slope_sd_ratio = 0)
  width <- function(s) {
    spec <- specify_target(dgp, "arm", engine = "bayes_laplace",
                           control = target_control(prior_tau_sd = s, on_singular = "flag_only"))
    evaluate_design_point(spec, 32, 12, nsim = 3, stream_seeds = 881)$logit_ci_width_mean
  }
  expect_lt(width(.1), width(.5))
  expect_error(width(0), "prior_tau_sd")
})
