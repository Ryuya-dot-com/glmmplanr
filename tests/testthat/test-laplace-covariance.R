test_that("uncorrelated covariance updates solve the Laplace marginal objective", {
  des <- design_crossed(list(condition = c("control", "treatment")),
    allocation = "latin_square", items_per_subject = 24L)
  dgp <- calibrate_dgp(des, c(control = .42, treatment = .50),
    subj_sd = c(1.1, .4), item_sd = c(.9, .5), corr_structure = "zero")
  spec <- specify_target(dgp, "condition", engine = "bayes_laplace",
    analysis_re = "match_dgp", prob_se_method = "joint_delta")
  for (n in c(24L, 72L)) {
    j <- if (n == 24L) 72L else 24L
    schedule <- generate_design_schedule(des, n, j)
    tr <- with_rng_stream(621700L + n, simulate_dgp_data(dgp, schedule))
    X <- des$X_cells[tr$cell_id, , drop = FALSE]
    Zs <- des$Z_subj_cells[tr$cell_id, , drop = FALSE]
    Zi <- des$Z_item_cells[tr$cell_id, , drop = FALSE]
    ws <- cpp_create_glmm_workspace(X, tr$subj_id - 1L, Zs, n,
      tr$item_id - 1L, Zi, j, integer(), matrix(0, nrow(tr), 0L), 0L)
    fit <- function(theta = NULL) cpp_fit_crossed_glmm_ws(ws, tr$y,
      corr_mode = 0L, regularized = TRUE, prior_beta_sd = spec$prior_beta_sd,
      prior_tau_sd = spec$prior_tau_sd, max_em_iter = 160L, max_pirls_iter = 22L,
      tol = 1.5e-4, fix_cov = !is.null(theta), return_ranef = TRUE,
      init_Sigma_subj = if (!is.null(theta)) diag(exp(2 * theta[1:2])) else NULL,
      init_Sigma_item = if (!is.null(theta)) diag(exp(2 * theta[3:4])) else NULL)
    f <- fit()
    expect_true(f$converged)
    expect_lt(f$covariance_score_norm, .005)
    tau <- sqrt(c(diag(f$Sigma_subj), diag(f$Sigma_item)))
    theta <- log(tau)

    # Independent likelihood: dense conditional Hessian in standardized effects.
    Z <- matrix(0, nrow(tr), 2L * (n + j))
    for (k in 1:2) {
      Z[cbind(seq_len(nrow(tr)), 2L * (tr$subj_id - 1L) + k)] <- Zs[, k] * tau[k]
      Z[cbind(seq_len(nrow(tr)), 2L * n + 2L * (tr$item_id - 1L) + k)] <- Zi[, k] * tau[k + 2L]
    }
    v <- c(as.vector(t(sweep(f$u_subj, 2L, tau[1:2], "/"))),
      as.vector(t(sweep(f$u_item, 2L, tau[3:4], "/"))))
    eta <- drop(X %*% f$beta + Z %*% v)
    p <- plogis(eta)
    H <- diag(ncol(Z)) + crossprod(Z, (p * (1 - p)) * Z)
    ll <- sum(tr$y * eta - pmax(eta, 0) - log1p(exp(-abs(eta)))) -
      .5 * sum(v^2) - sum(log(diag(chol(H))))
    expect_equal(f$logLik, ll, tolerance = 1e-8)

    objective <- function(t) {
      z <- fit(t)
      expect_true(z$converged)
      -z$logLik - sum(dnorm(z$beta, 0, c(1.5, 1) * spec$prior_beta_sd, log = TRUE)) -
        sum(dt(exp(t) / spec$prior_tau_sd, 3, log = TRUE)) - sum(t)
    }
    h <- 1e-3
    gradient <- vapply(1:4, function(k) {
      step <- replace(numeric(4), k, h)
      (objective(theta + step) - objective(theta - step)) / (2 * h)
    }, numeric(1))
    info_sd <- sqrt(c(f$var_log_tau_subj, f$var_log_tau_item))
    # Small numerical differences arise from profiling beta and finite differences.
    expect_lt(max(abs(gradient * info_sd)), .01)
  }
})

test_that("testlet covariance scores agree with penalized and unpenalized profiles", {
  des <- design_testlet(list(condition = c("control", "treatment")), testlet_size = 4L)
  dgp <- calibrate_dgp(des, c(control = .42, treatment = .50),
    subj_sd = c(.8, 0), item_sd = c(.6, 0), corr_structure = "zero",
    control = calibrate_control(testlet_sd = c(.7, 0)))
  schedule <- generate_design_schedule(des, 48L, 48L)
  tr <- with_rng_stream(621799L, simulate_dgp_data(dgp, schedule))
  X <- des$X_cells[tr$cell_id, , drop = FALSE]
  Z <- matrix(1, nrow(tr), 1L)
  ws <- cpp_create_glmm_workspace(X, tr$subj_id - 1L, Z, 48L,
    tr$item_id - 1L, Z, 48L, tr$testlet_id - 1L, Z, schedule$n_testlet)
  for (regularized in c(TRUE, FALSE)) {
    fit <- function(theta = NULL) cpp_fit_crossed_glmm_ws(ws, tr$y,
      corr_mode = 0L, regularized = regularized, fix_cov = !is.null(theta),
      init_Sigma_subj = if (!is.null(theta)) matrix(exp(2 * theta[1]), 1L),
      init_Sigma_item = if (!is.null(theta)) matrix(exp(2 * theta[2]), 1L),
      init_Sigma_testlet = if (!is.null(theta)) matrix(exp(2 * theta[3]), 1L))
    f <- fit()
    expect_true(f$converged)
    expect_lt(f$covariance_score_norm, .005)
    theta <- log(sqrt(c(f$Sigma_subj[1, 1], f$Sigma_item[1, 1], f$Sigma_testlet[1, 1])))
    objective <- function(t) {
      z <- fit(t)
      expect_true(z$converged)
      prior <- if (regularized) {
        sum(dnorm(z$beta, 0, c(1.5, 1) * 2.5, log = TRUE)) +
          sum(dt(exp(t) / 1.5, 3, log = TRUE)) + sum(t)
      } else 0
      -z$logLik - prior
    }
    h <- 1e-3
    gradient <- vapply(1:3, function(k) {
      step <- replace(numeric(3), k, h)
      (objective(theta + step) - objective(theta - step)) / (2 * h)
    }, numeric(1))
    info_sd <- sqrt(c(f$var_log_tau_subj, f$var_log_tau_item, f$var_log_tau_testlet))
    expect_lt(max(abs(gradient * info_sd)), .01)
  }
})

test_that("early screening does not reject a valid corrected covariance fit", {
  des <- design_crossed(list(study_method = c("restudy", "retrieval")),
    allocation = "latin_square", items_per_subject = 30L, item_seed = 37L)
  dgp <- calibrate_dgp(des, c(restudy = .42, retrieval = .47),
    subj_sd = c(.8, .25), item_sd = c(.6, .35), corr_structure = "zero")
  spec <- specify_target(dgp, "study_method", scale = "pp", criterion = "joint",
    ci_width = 10, engine = "bayes_laplace", analysis_re = "match_dgp",
    prob_se_method = "joint_delta")
  schedule <- generate_design_schedule(des, 162L, 60L)
  cache <- prepare_cpp_schedule_cache(spec, schedule)
  bases <- allocate_stage_point_streams(912210001L, 1L, 1:2, 2L)
  stream <- derive_rep_substreams(bases[[2L]], 1L)[[1L]]
  sim <- with_rng_stream(stream, simulate_fast_y_keep(cache))
  final <- fit_via_cpp_cached(sim, cache, spec, TRUE, spec$analysis_re, FALSE)
  early <- fit_via_cpp_cached(sim, cache, spec, TRUE, spec$analysis_re, TRUE)
  expect_true(final$converged)
  expect_true(early$converged)
  fields <- c("est", "se", "prob_est", "prob_se", "Sigma_subj", "Sigma_item", "status")
  expect_equal(early[fields], final[fields], tolerance = 1e-10)
})

test_that("large contrasts move beta and satisfy the marginal beta score", {
  des <- design_crossed(list(condition = c("a", "b")))
  dgp <- calibrate_dgp(des, c(a = .1, b = .9), subj_sd = .4,
    item_sd = .3, corr_structure = "zero")
  schedule <- generate_design_schedule(des, 48L, 24L)
  tr <- with_rng_stream(37111L, simulate_dgp_data(dgp, schedule))
  X <- des$X_cells[tr$cell_id, , drop = FALSE]
  one <- matrix(1, nrow(tr), 1L)
  ws <- cpp_create_glmm_workspace(X, tr$subj_id - 1L, one, 48L,
    tr$item_id - 1L, one, 24L, integer(), matrix(0, nrow(tr), 0L), 0L)
  responses <- list(mixed = tr$y, separated = as.integer(tr$cell_id == 2L))
  for (kind in names(responses)) {
    y <- responses[[kind]]
    f <- cpp_fit_crossed_glmm_ws(ws, y, corr_mode = 0L,
      regularized = TRUE, return_ranef = TRUE)
    if (kind == "mixed") {
      expect_true(f$converged)
    } else {
      # A prior permits finite estimates; it does not guarantee convergence
      # within the iteration budget for fully separated data.
      expect_true(f$status %in% c("converged", "max_iter_finite"))
    }
    expect_gt(sum((des$X_cells[2L, ] - des$X_cells[1L, ]) * f$beta), 3)

    # Independent dense conditional-mode solve at perturbed beta values.
    tau <- sqrt(c(f$Sigma_subj[1, 1], f$Sigma_item[1, 1]))
    Z <- matrix(0, nrow(tr), 72L)
    Z[cbind(seq_len(nrow(tr)), tr$subj_id)] <- tau[1L]
    Z[cbind(seq_len(nrow(tr)), 48L + tr$item_id)] <- tau[2L]
    v0 <- c(f$u_subj / tau[1L], f$u_item / tau[2L])
    log_posterior <- function(beta) {
      v <- v0
      for (k in 1:20) {
        eta <- drop(X %*% beta + Z %*% v)
        p <- plogis(eta)
        H <- diag(72L) + crossprod(Z, (p * (1 - p)) * Z)
        step <- solve(H, drop(crossprod(Z, y - p)) - v)
        v <- v + step
        if (max(abs(step)) < 1e-11) break
      }
      eta <- drop(X %*% beta + Z %*% v)
      p <- plogis(eta)
      H <- diag(72L) + crossprod(Z, (p * (1 - p)) * Z)
      sum(y * eta - pmax(eta, 0) - log1p(exp(-abs(eta)))) -
        .5 * sum(v^2) - sum(log(diag(chol(H)))) +
        sum(dnorm(beta, 0, c(3.75, 2.5), log = TRUE))
    }
    h <- 1e-4
    gradient <- vapply(1:2, function(k) {
      step <- replace(numeric(2), k, h)
      (log_posterior(f$beta + step) - log_posterior(f$beta - step)) / (2 * h)
    }, numeric(1))
    expect_lt(max(abs(gradient * sqrt(diag(f$vcov)))), .01)
  }
})
