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

test_that("calibration preserves a fresh and a Box-Muller caller, including the next normal", {
  # Do not alter the stream of the test runner itself.
  pkg_path <- getNamespaceInfo(asNamespace("glmmplanr"), "path")
  checked <- callr::r(function(path) {
    if (!file.exists(file.path(path, "Meta", "package.rds"))) pkgload::load_all(path, quiet = TRUE)
    else library(glmmplanr)
    d <- design_crossed(list())
    if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    model <- calibrate_dgp(d, .6, subj_sd = .4, item_sd = .3, check_nsim = 20L)
    empty <- !exists(".Random.seed", .GlobalEnv, inherits = FALSE)
    set.seed(18, normal.kind = "Box-Muller")
    expected <- rnorm(2)
    set.seed(18, normal.kind = "Box-Muller")
    first <- rnorm(1)
    before <- .Random.seed
    copy <- calibrate_dgp(d, .6, subj_sd = .4, item_sd = .3, check_nsim = 20L)
    preserved <- identical(.Random.seed, before) && identical(c(first, rnorm(1)), expected)
    set.seed(18, normal.kind = "Inversion")
    sample.int(.Machine$integer.max, 1L)
    after_one <- .Random.seed
    set.seed(18)
    invisible(check_dgp(model, n_draws = 20L, seed = NULL))
    c(empty, preserved, isTRUE(all.equal(copy$prior_check, model$prior_check, tolerance = 1e-8)),
      identical(after_one, .Random.seed))
  }, list(path = pkg_path))
  expect_true(all(checked))
})

test_that("wide normal effects recover the independently integrated probabilities", {
  des <- design_crossed(list(arm = c("a", "b")), allocation = "between_subj",
                        contrast = "treatment")
  g <- calibrate_dgp(des, c(.72, .80), chance = 1/3, subj_sd = 5, item_sd = 1,
                     slope_sd_ratio = 0)
  ref <- vapply(g$mu_cells, function(eta) {
    integrate(function(z) (1/3 + 2/3 * plogis(eta + sqrt(26) * z)) * dnorm(z),
              -Inf, Inf, rel.tol = 1e-12)$value
  }, numeric(1))
  expect_lt(max(abs(ref - c(.72, .80))), 1e-7)
  expect_lt(abs(true_prob_effect(g, c(-1, 1)) - diff(ref)), 2e-7)
  expect_gt(g$n_quad, g$solver_info$requested_nodes)
})
