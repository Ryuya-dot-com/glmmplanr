test_that("3AFC estimates agree with the independent Stan fit", {
  skip_if_not(identical(Sys.getenv("GLMMPLANR_TEST_STAN"), "true"),
              "Set GLMMPLANR_TEST_STAN=true to run the independent Stan fit")
  skip_if_not(has_cmdstan_backend() || requireNamespace("rstan", quietly = TRUE),
              "No Stan backend is installed")
  dgp <- three_choice_reference_dgp()
  spec <- specify_target(dgp, effect = "exposure", engine = "bayes_laplace",
                         control = target_control(on_singular = "flag_only"))
  schedule <- generate_design_schedule(dgp$design, n_subj = 48, n_item = 32)
  set.seed(303, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion", sample.kind = "Rejection")
  data <- simulate_dgp_data(dgp, schedule = schedule)
  cpp <- fit_single_dataset(data, spec, engine_use = "bayes_laplace")
  stan <- fit_single_dataset(data, spec, engine_use = "stan_laplace")
  expect_true(cpp$converged, info = cpp$error_reason)
  expect_true(stan$converged, info = stan$error_reason)
  expect_equal(cpp$est, stan$est, tolerance = 0.15)
  expect_equal(cpp$se, stan$se, tolerance = 0.06)
})
