test_that("roundoff does not stall a fixed-covariance 3AFC profile", {
  # Frozen synthetic case 12, replicate 1 from the 128-replication review,
  # with log-SD perturbations (+.01 on participant slope, -.01 on item slope).
  a <- readRDS(test_path("fixtures", "mafc-mode-roundoff.rds"))
  f <- do.call(cpp_fit_crossed_glmm, c(a, list(chance = 1/3, corr_mode = 0L,
    regularized = TRUE, fix_cov = TRUE, max_em_iter = 400L, max_pirls_iter = 44L,
    compute_joint_vcov = FALSE, return_ranef = TRUE)))
  expect_true(f$converged)
  expect_lt(f$inner_step, 2e-7)
  expect_true(is.na(f$covariance_score_norm))

  # Independent observed Hessian and conditional-mode score in standardized effects.
  Ls <- t(chol(f$Sigma_subj)); Li <- t(chol(f$Sigma_item))
  Z <- matrix(0, length(a$y), 2L * (a$n_subj + a$n_item))
  for (j in 1:2) {
    Z[cbind(seq_along(a$y), 2L * a$subj_idx + j)] <- (a$Z_subj %*% Ls)[, j]
    Z[cbind(seq_along(a$y), 2L * (a$n_subj + a$item_idx) + j)] <- (a$Z_item %*% Li)[, j]
  }
  v <- c(as.vector(solve(Ls, t(f$u_subj))), as.vector(solve(Li, t(f$u_item))))
  p <- plogis(drop(a$X %*% f$beta + Z %*% v)); mu <- 1/3 + 2/3 * p
  dp <- 2/3 * p * (1-p); d2p <- dp * (1-2*p)
  score <- (a$y-mu) * dp / (mu * (1-mu))
  w <- (dp^2-(a$y-mu)*d2p)/(mu*(1-mu)) +
    (a$y-mu)*dp^2*(1-2*mu)/(mu^2*(1-mu)^2)
  H <- diag(length(v)) + crossprod(Z, w * Z)
  ll <- sum(dbinom(a$y, 1, mu, log = TRUE)) - .5 * sum(v^2) - sum(log(diag(chol(H))))
  expect_equal(f$logLik, ll, tolerance = 1e-8)
  expect_lt(max(abs(solve(H, drop(crossprod(Z, score))-v))), 2e-7)
})

test_that("unavailable native diagnostics are missing rather than zero or working likelihoods", {
  a <- readRDS(test_path("fixtures", "mafc-mode-roundoff.rds"))
  for (chance in c(0, 1/3)) for (corr in 0:1) {
    f <- do.call(cpp_fit_crossed_glmm, c(a, list(chance = chance, corr_mode = corr,
      regularized = TRUE, max_em_iter = 400L, max_pirls_iter = 44L)))
    expect_true(f$converged)
    expect_true(is.finite(f$logLik))
    if (chance == 0 && corr == 0) expect_lt(f$covariance_score_norm, .005)
    else expect_true(is.na(f$covariance_score_norm))
  }
  failed <- do.call(cpp_fit_crossed_glmm, c(a, list(chance = 1/3, corr_mode = 0L,
    regularized = TRUE, max_em_iter = 1L, max_pirls_iter = 1L)))
  expect_false(failed$converged)
  expect_true(is.na(failed$logLik))
  expect_true(is.na(failed$covariance_score_norm))
})

test_that("analysis descriptions distinguish moment convergence from covariance scores", {
  des <- design_crossed(list(condition = c("a", "b")))
  for (g in c(0, 1/3)) for (rho in c(0, .2)) {
    dgp <- calibrate_dgp(des, c(a = .55, b = .65), chance = g,
      subj_sd = c(.8, .35), item_sd = c(.6, .4), corr_structure = rho)
    sp <- specify_target(dgp, "condition")
    description <- analysis_description(sp)
    expect_identical(grepl("approximate moment updates", description), g > 0 || rho != 0)
    if (g > 0 || rho != 0) expect_match(description, "does not include a marginal covariance-score check")
  }
})
