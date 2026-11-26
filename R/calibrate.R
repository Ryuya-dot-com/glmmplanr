#' Calibrate Conditional GLMM Parameters from Marginal Probabilities and m-AFC Guessing
#'
#' Numerically inverts researcher-specified marginal cell response probabilities
#' \eqn{P(Y = 1 \mid X)} into data-generating conditional fixed-effect
#' coefficients \eqn{\beta} on the guessing-adjusted logit scale using
#' exact 1D Gauss-Hermite quadrature over the combined normal random-effect
#' variance per cell and `nleqslv`, correcting both non-collapsibility
#' (random-effect attenuation toward the midpoint) and lower guessing
#' asymptotes (\eqn{g = 1/m} in forced-choice tasks).
#'
#' @param design A `"glmmplanr_design"` object created by [design_crossed()] or
#'   [design_testlet()].
#' @param cell_probs Numeric vector of target marginal accuracy probabilities
#'   (one per design cell), strictly inside `(chance, 1)`. Strongly recommended
#'   to pass a named vector whose names match factor level combinations joined
#'   by `"_"` or `":"` (e.g., `c(incidental_massed = 0.50, ...)`).
#' @param chance Guessing probability (lower asymptote \eqn{g \in [0, 1)},
#'   default `0`). Use `0` for an ordinary binary logit model -- including
#'   recall, production, yes/no or old/new judgments (even though there are two
#'   response options), and multiple-choice tasks analyzed with standard
#'   `lme4::glmer(..., family = binomial())`. Set `chance = 1/m` (e.g., `1/3`
#'   for 3AFC, `0.25` for 4AFC) only when you intend to model and analyze a
#'   fixed lower guessing asymptote \eqn{g + (1 - g)\,\mathrm{logit}^{-1}(\eta)}.
#' @param subj_sd Optional numeric scalar (intercept SD) or vector (intercept +
#'   slope SDs) on the conditional logit scale for subjects.
#' @param item_sd Optional numeric scalar (intercept SD) or vector (intercept +
#'   slope SDs) on the conditional logit scale for items.
#' @param subj_prob_range Optional numeric vector `c(low, high)` giving the 95%
#'   subject accuracy range used to specify the intercept SD via its width on
#'   the adjusted logit scale. Calibration sets the location from `cell_probs`,
#'   so the fitted endpoints can shift. Ignored when `subj_sd` is supplied.
#' @param item_prob_range Optional numeric vector `c(low, high)` giving the 95%
#'   item accuracy range used in the same way. Ignored if `item_sd` is
#'   supplied directly.
#' @param slope_sd_ratio Numeric scalar or vector specifying the ratio of random
#'   slope SDs to the corresponding random intercept SD (default `0`, i.e.,
#'   random intercepts only unless you explicitly specify a positive ratio such
#'   as `0.35` or pass full SD vectors in `subj_sd` / `item_sd`).
#' @param corr_structure Correlation specification for random effects: a single
#'   equicorrelation scalar \eqn{\rho} (default `0.20`), `"zero"` (uncorrelated
#'   random effects), a correlation matrix, or a named list
#'   `list(subj = ..., item = ...)`.
#' @param control Control list created by [calibrate_control()] holding optional
#'   testlet variance arguments (`testlet_sd`, `testlet_prob_range`,
#'   `participant_testlet_sd`) and numerical settings (`n_quad`, `check_nsim`,
#'   `seed`).
#'
#' @return An S3 object of class `"glmmplanr_dgp"` containing the calibrated
#'   conditional parameter vector `beta`, naive uncalibrated benchmarks,
#'   random-effect covariance matrices (`Sigma_subj`, `Sigma_item`,
#'   `Sigma_testlet`) and Cholesky factors, solver diagnostics, and prior
#'   predictive check results.
#'
#' @examples
#' des <- design_crossed(
#'   factors = list(
#'     exposure = c("incidental", "intentional"),
#'     spacing  = c("massed", "distributed")
#'   ),
#'   allocation = "latin_square"
#' )
#' dgp <- calibrate_dgp(
#'   des,
#'   cell_probs = c(
#'     incidental_massed       = 0.50,
#'     incidental_distributed  = 0.60,
#'     intentional_massed      = 0.65,
#'     intentional_distributed = 0.80
#'   ),
#'   chance          = 1 / 3,
#'   subj_prob_range = c(0.40, 0.85),
#'   item_prob_range = c(0.45, 0.88),
#'   slope_sd_ratio  = 0.35,
#'   corr_structure  = 0.20
#' )
#' print(dgp)
#' @export
calibrate_dgp <- function(
  design,
  cell_probs,
  chance = 0,
  subj_sd = NULL,
  item_sd = NULL,
  subj_prob_range = NULL,
  item_prob_range = NULL,
  slope_sd_ratio = 0,
  corr_structure = 0.20,
  control = calibrate_control()
) {
  if (!inherits(control, "glmmplanr_calibrate_control")) {
    cli::cli_abort("{.arg control} must be created by {.fn calibrate_control}.")
  }
  if (uses_box_muller()) {
    if (is.null(control$seed)) control$seed <- sample.int(.Machine$integer.max, 1L)
    return(isolated_seeded_call("calibrate_dgp", environment(), match.call()))
  }
  testlet_sd <- control$testlet_sd
  testlet_prob_range <- control$testlet_prob_range
  participant_testlet_sd <- control$participant_testlet_sd
  n_quad <- control$n_quad
  check_nsim <- control$check_nsim
  seed <- control$seed

  if (!inherits(design, "glmmplanr_design")) {
    cli::cli_abort("{.arg design} must be a {.cls glmmplanr_design} object.")
  }
  if (participant_testlet_sd > 0 && !design$has_testlet) {
    cli::cli_abort("{.arg participant_testlet_sd} must be nonnegative and requires a testlet design when positive.")
  }
  validate_chance(chance)

  target_probs <- align_cell_probs(design, cell_probs, chance)

  # Resolve subject and item random-effect SD vectors
  p_subj <- ncol(design$Z_subj_cells)
  p_item <- ncol(design$Z_item_cells)

  tau_subj <- resolve_re_sds(
    prob_range = subj_prob_range,
    direct_sd = subj_sd,
    n_re = p_subj,
    slope_sd_ratio = slope_sd_ratio,
    chance = chance,
    re_names = colnames(design$Z_subj_cells),
    label = "subj"
  )
  tau_item <- resolve_re_sds(
    prob_range = item_prob_range,
    direct_sd = item_sd,
    n_re = p_item,
    slope_sd_ratio = slope_sd_ratio,
    chance = chance,
    re_names = colnames(design$Z_item_cells),
    label = "item"
  )

  corr_subj_spec <- if (is.list(corr_structure) && !is.null(corr_structure$subj)) {
    corr_structure$subj
  } else {
    corr_structure
  }
  corr_item_spec <- if (is.list(corr_structure) && !is.null(corr_structure$item)) {
    corr_structure$item
  } else {
    corr_structure
  }

  R_subj <- build_corr_matrix(p_subj, corr_subj_spec, colnames(design$Z_subj_cells), "subj")
  R_item <- build_corr_matrix(p_item, corr_item_spec, colnames(design$Z_item_cells), "item")

  Sigma_subj <- diag(tau_subj, nrow = p_subj) %*% R_subj %*% diag(tau_subj, nrow = p_subj)
  Sigma_item <- diag(tau_item, nrow = p_item) %*% R_item %*% diag(tau_item, nrow = p_item)
  dimnames(Sigma_subj) <- list(colnames(design$Z_subj_cells), colnames(design$Z_subj_cells))
  dimnames(Sigma_item) <- list(colnames(design$Z_item_cells), colnames(design$Z_item_cells))

  L_subj <- safe_chol_lower(Sigma_subj)
  L_item <- safe_chol_lower(Sigma_item)

  if (design$has_testlet) {
    p_testlet <- ncol(design$Z_testlet_cells)
    if (is.null(testlet_prob_range) && is.null(testlet_sd)) {
      testlet_sd <- 0.5 * tau_item[1L]
    }
    tau_testlet <- resolve_re_sds(
      prob_range = testlet_prob_range,
      direct_sd = testlet_sd,
      n_re = p_testlet,
      slope_sd_ratio = slope_sd_ratio,
      chance = chance,
      re_names = colnames(design$Z_testlet_cells),
      label = "testlet"
    )
    corr_testlet_spec <- if (is.list(corr_structure) && !is.null(corr_structure$testlet)) {
      corr_structure$testlet
    } else {
      corr_structure
    }
    R_testlet <- build_corr_matrix(
      p_testlet, corr_testlet_spec, colnames(design$Z_testlet_cells), "testlet"
    )
    Sigma_testlet <- diag(tau_testlet, nrow = p_testlet) %*% R_testlet %*%
      diag(tau_testlet, nrow = p_testlet)
    dimnames(Sigma_testlet) <- list(
      colnames(design$Z_testlet_cells), colnames(design$Z_testlet_cells)
    )
    L_testlet <- safe_chol_lower(Sigma_testlet)
  } else {
    tau_testlet <- numeric(0)
    R_testlet <- NULL
    Sigma_testlet <- NULL
    L_testlet <- NULL
  }

  # Cell-specific random-effect total SDs along subject and item(+testlet) axes
  var_subj_cell <- rowSums((design$Z_subj_cells %*% Sigma_subj) * design$Z_subj_cells)
  var_item_cell <- rowSums((design$Z_item_cells %*% Sigma_item) * design$Z_item_cells)
  if (design$has_testlet) {
    var_testlet_cell <- rowSums(
      (design$Z_testlet_cells %*% Sigma_testlet) * design$Z_testlet_cells
    )
  } else {
    var_testlet_cell <- rep(0, design$n_cells)
  }
  var_testlet_cell <- var_testlet_cell + participant_testlet_sd^2
  sd_subj_cell <- sqrt(pmax(var_subj_cell, 0))
  sd_item_tot_cell <- sqrt(pmax(var_item_cell + var_testlet_cell, 0))

  # Solve for conditional beta via exact 1D Gauss-Hermite quadrature + nleqslv.
  # Independent normal participant, item, and testlet effects sum to N(0, sd_total^2) in each cell.
  if (!is.numeric(n_quad) || length(n_quad) != 1L || !is.finite(n_quad) ||
      n_quad < 1 || n_quad != floor(n_quad)) {
    cli::cli_abort("{.arg n_quad} must be a positive integer.")
  }
  sd_total <- sqrt(sd_subj_cell^2 + sd_item_tot_cell^2)
  X_cells <- design$X_cells

  # Analytical starting point using logistic-normal attenuation approximation
  c_att <- (16 * sqrt(3)) / (15 * pi)
  eta_naive_mafc <- mafc_qlogis(target_probs, chance = chance)
  eta_init <- eta_naive_mafc * sqrt(1 + (c_att^2) * (var_subj_cell + var_item_cell + var_testlet_cell))
  beta_init <- as.numeric(qr.solve(X_cells, eta_init))

  eval_marginal <- function(beta_vec) {
    mu_cell <- as.numeric(X_cells %*% beta_vec)
    # Vectorized across cells (M x Q)
    linpred <- outer(mu_cell, rep.int(1, length(quad1d$weights))) +
      outer(sd_total, quad1d$nodes)
    p_inner <- stats::plogis(linpred)
    probs <- chance + (1 - chance) * as.numeric(p_inner %*% quad1d$weights)
    probs - target_probs
  }

  eval_jacobian <- function(beta_vec) {
    mu_cell <- as.numeric(X_cells %*% beta_vec)
    linpred <- outer(mu_cell, rep.int(1, length(quad1d$weights))) +
      outer(sd_total, quad1d$nodes)
    p_inner <- stats::plogis(linpred)
    deriv_inner <- (1 - chance) * p_inner * (1 - p_inner)
    d_mu <- as.numeric(deriv_inner %*% quad1d$weights)
    X_cells * d_mu
  }

  requested_nodes <- n_quad
  repeat {
    quad1d <- gh_quad_1d(n_quad)
    sol <- nleqslv::nleqslv(beta_init, eval_marginal, jac = eval_jacobian,
      method = "Newton", control = list(ftol = 1e-13, xtol = 1e-13, maxit = 200L))
    max_err <- max(abs(sol$fvec))
    if (!is.finite(max_err) || max_err > 1e-7) {
      cli::cli_abort("Marginal inversion failed (maximum solver residual {signif(max_err, 4)}). Check the target probabilities and SDs.")
    }
    mu <- drop(X_cells %*% sol$x)
    reference <- vapply(seq_along(mu), function(i) {
      stats::integrate(function(z) mafc_plogis(mu[i] + sd_total[i] * z, chance) *
        stats::dnorm(z), -Inf, Inf, rel.tol = 1e-12)$value
    }, numeric(1))
    reported <- drop(mafc_plogis(mu + outer(sd_total, quad1d$nodes), chance) %*% quad1d$weights)
    integration_error <- max(abs(reference - target_probs), abs(reported - reference))
    if (integration_error <= 1e-10) break
    if (n_quad >= 256L) {
      if (integration_error <= 1e-7) break
      cli::cli_abort("Quadrature differs from adaptive integration by {signif(integration_error, 4)} at {n_quad} nodes. Increase {.arg n_quad} or review the random-effect scales.")
    }
    beta_init <- sol$x
    n_quad <- min(256L, 2L * n_quad)
  }

  beta <- stats::setNames(sol$x, colnames(X_cells))
  beta_naive_mafc <- stats::setNames(
    as.numeric(qr.solve(X_cells, eta_naive_mafc)),
    colnames(X_cells)
  )
  beta_naive_logit <- stats::setNames(
    as.numeric(qr.solve(X_cells, stats::qlogis(target_probs))),
    colnames(X_cells)
  )
  mu_cells <- stats::setNames(as.numeric(X_cells %*% beta), design$cells$cell_name)
  achieved_probs <- stats::setNames(reference, design$cells$cell_name)

  dgp_obj <- structure(
    list(
      design = design,
      target_cell_probs = target_probs,
      achieved_cell_probs = achieved_probs,
      chance = chance,
      beta = beta,
      beta_naive_mafc = beta_naive_mafc,
      beta_naive_logit = beta_naive_logit,
      mu_cells = mu_cells,
      tau_subj = tau_subj,
      tau_item = tau_item,
      tau_testlet = tau_testlet,
      participant_testlet_sd = participant_testlet_sd,
      R_subj = R_subj,
      R_item = R_item,
      R_testlet = R_testlet,
      Sigma_subj = Sigma_subj,
      Sigma_item = Sigma_item,
      Sigma_testlet = Sigma_testlet,
      L_subj = L_subj,
      L_item = L_item,
      L_testlet = L_testlet,
      slope_sd_ratio = slope_sd_ratio,
      corr_structure = corr_structure,
      subj_prob_range = subj_prob_range,
      item_prob_range = item_prob_range,
      n_quad = as.integer(n_quad),
      solver_info = list(
        requested_nodes = requested_nodes,
        integration_error = integration_error,
        iterations = sol$iter,
        max_abs_residual = max_err,
        termcd = sol$termcd,
        message = sol$message
      )
    ),
    class = "glmmplanr_dgp"
  )

  dgp_obj$prior_check <- compute_prior_predictive_check(
    dgp_obj, n_draws = check_nsim, seed = seed
  )
  dgp_obj
}

align_cell_probs <- function(design, cell_probs, chance) {
  M <- design$n_cells
  if (!is.numeric(cell_probs) || length(cell_probs) != M || anyNA(cell_probs)) {
    cli::cli_abort(
      "{.arg cell_probs} must be a numeric vector of length {M} matching the design cells: {.val {design$cells$cell_name}}."
    )
  }
  if (any(cell_probs <= chance) || any(cell_probs >= 1)) {
    hints <- c(
      "All {.arg cell_probs} must strictly exceed {.arg chance} ({round(chance, 4)}) and be < 1."
    )
    if (any(cell_probs > 1)) {
      hints <- c(hints, "i" = "Enter probabilities on the 0-1 scale (e.g., {.code 0.60} for 60%), not percentages.")
    }
    if (chance > 0 && any(cell_probs <= chance)) {
      hints <- c(
        hints,
        "i" = "With {.arg chance} = {round(chance, 4)}, every expected cell probability must exceed the guessing floor. For yes/no judgments or multiple-choice tasks analyzed with standard {.code lme4::glmer(..., family = binomial())}, set {.code chance = 0}."
      )
    }
    cli::cli_abort(hints)
  }

  canonical <- design$cells$cell_name
  nms <- names(cell_probs)
  if (!is.null(nms) && all(nzchar(nms))) {
    norm_nms <- gsub("[:./\\s]+", "_", trimws(nms), perl = TRUE)
    idx <- match(canonical, norm_nms)
    if (anyNA(idx)) {
      cli::cli_abort(c(
        "Names of {.arg cell_probs} do not match expected design cells.",
        "i" = "Expected: {.val {canonical}}",
        "x" = "Provided: {.val {nms}}"
      ))
    }
    out <- stats::setNames(as.numeric(cell_probs[idx]), canonical)
  } else {
    out <- stats::setNames(as.numeric(cell_probs), canonical)
  }
  out
}

resolve_re_sds <- function(
  prob_range,
  direct_sd,
  n_re,
  slope_sd_ratio,
  chance,
  re_names,
  label
) {
  if (!is.null(direct_sd)) {
    if (!is.numeric(direct_sd) || any(!is.finite(direct_sd)) || any(direct_sd < 0)) {
      cli::cli_abort("{.arg {label}_sd} must contain non-negative numeric SD values.")
    }
    if (length(direct_sd) == 1L) {
      tau0 <- direct_sd
      tau <- build_tau_from_ratio(tau0, n_re, slope_sd_ratio)
    } else if (length(direct_sd) == n_re) {
      tau <- as.numeric(direct_sd)
    } else {
      cli::cli_abort(
        "{.arg {label}_sd} must have length 1 (intercept SD) or {n_re} (all RE terms: {.val {re_names}})."
      )
    }
  } else if (!is.null(prob_range)) {
    tau0 <- prob_range_to_sd(prob_range, chance = chance)
    tau <- build_tau_from_ratio(tau0, n_re, slope_sd_ratio)
  } else {
    cli::cli_abort(
      "Supply either {.arg {label}_prob_range} or {.arg {label}_sd} to define {label} random-effect variability."
    )
  }
  stats::setNames(tau, re_names)
}

build_tau_from_ratio <- function(tau0, n_re, slope_sd_ratio) {
  if (n_re == 1L) return(as.numeric(tau0))
  if (!is.numeric(slope_sd_ratio) || anyNA(slope_sd_ratio) || any(slope_sd_ratio < 0)) {
    cli::cli_abort("{.arg slope_sd_ratio} must be non-negative numeric.")
  }
  if (length(slope_sd_ratio) == 1L) {
    c(tau0, rep(tau0 * slope_sd_ratio, n_re - 1L))
  } else if (length(slope_sd_ratio) == n_re - 1L) {
    c(tau0, tau0 * slope_sd_ratio)
  } else {
    cli::cli_abort(
      "{.arg slope_sd_ratio} must have length 1 or {n_re - 1L} (matching the number of random slopes)."
    )
  }
}

build_corr_matrix <- function(n_re, corr_spec, re_names, label) {
  if (n_re == 1L) {
    R <- matrix(1, nrow = 1L, ncol = 1L, dimnames = list(re_names, re_names))
    return(R)
  }
  if (identical(corr_spec, "zero") || identical(corr_spec, 0)) {
    R <- diag(1, nrow = n_re)
  } else if (is.numeric(corr_spec) && length(corr_spec) == 1L && !is.matrix(corr_spec)) {
    rho <- as.numeric(corr_spec)
    min_rho <- -1 / (n_re - 1L) + 1e-6
    if (is.na(rho) || rho <= min_rho || rho >= 1) {
      cli::cli_abort(
        "Equicorrelation {.arg corr_structure} ({rho}) for {label} ({n_re} RE terms) must lie strictly in ({round(min_rho, 3)}, 1) for positive definiteness."
      )
    }
    R <- matrix(rho, nrow = n_re, ncol = n_re)
    diag(R) <- 1
  } else if (is.matrix(corr_spec)) {
    if (!all(dim(corr_spec) == c(n_re, n_re))) {
      cli::cli_abort(
        "Correlation matrix for {label} must have dimensions {n_re} x {n_re}."
      )
    }
    R <- corr_spec
    if (max(abs(R - t(R))) > 1e-8 || max(abs(diag(R) - 1)) > 1e-8) {
      cli::cli_abort("Correlation matrix for {label} must be symmetric with unit diagonal.")
    }
    ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) <= 1e-8) {
      cli::cli_abort("Correlation matrix for {label} must be strictly positive definite.")
    }
  } else {
    cli::cli_abort(
      "Unsupported {.arg corr_structure} for {label}; use a scalar rho, {.val zero}, or a correlation matrix."
    )
  }
  dimnames(R) <- list(re_names, re_names)
  R
}

safe_chol_lower <- function(Sigma) {
  p <- nrow(Sigma)
  if (all(Sigma == 0)) {
    return(matrix(0, nrow = p, ncol = p, dimnames = dimnames(Sigma)))
  }
  # Handle semi-definite zero-variance slopes cleanly via pivoted Cholesky
  ch <- suppressWarnings(chol(Sigma, pivot = TRUE))
  piv <- attr(ch, "pivot")
  rnk <- attr(ch, "rank")
  if (rnk < p) {
    ch[(rnk + 1L):p, ] <- 0
  }
  inv_piv <- order(piv)
  U <- ch[, inv_piv, drop = FALSE]
  L_out <- t(U)
  dimnames(L_out) <- dimnames(Sigma)
  L_out
}

#' Prior-Predictive Sanity Check of Calibrated DGP
#' @noRd
compute_prior_predictive_check <- function(dgp, n_draws = 4000L, seed = 1L) {
  streams <- make_lecuyer_streams(seed, 1L)
  with_rng_stream(streams[[1L]], {
    des <- dgp$design
    p_s <- ncol(des$Z_subj_cells)
    p_i <- ncol(des$Z_item_cells)
    u_s <- matrix(stats::rnorm(n_draws * p_s), nrow = n_draws, ncol = p_s) %*% t(dgp$L_subj)
    u_i <- matrix(stats::rnorm(n_draws * p_i), nrow = n_draws, ncol = p_i) %*% t(dgp$L_item)
    if (des$has_testlet) {
      p_t <- ncol(des$Z_testlet_cells)
      u_t <- matrix(stats::rnorm(n_draws * p_t), nrow = n_draws, ncol = p_t) %*% t(dgp$L_testlet)
    }

    local_effect <- draw_testlet_intercepts(dgp, n_draws)
    cell_summary <- lapply(seq_len(des$n_cells), function(c_idx) {
      zs <- des$Z_subj_cells[c_idx, ]
      zi <- des$Z_item_cells[c_idx, ]
      eta_draw <- dgp$mu_cells[c_idx] + as.numeric(u_s %*% zs) + as.numeric(u_i %*% zi)
      if (des$has_testlet) {
        zt <- des$Z_testlet_cells[c_idx, ]
        eta_draw <- eta_draw + as.numeric(u_t %*% zt)
      }
      prob_draw <- mafc_plogis(eta_draw + local_effect, chance = dgp$chance)
      data.frame(
        cell = des$cells$cell_name[c_idx],
        target_prob = unname(dgp$target_cell_probs[c_idx]),
        quad_prob = unname(dgp$achieved_cell_probs[c_idx]),
        mc_mean_prob = mean(prob_draw),
        q025_prob = unname(stats::quantile(prob_draw, 0.025)),
        q500_prob = unname(stats::quantile(prob_draw, 0.500)),
        q975_prob = unname(stats::quantile(prob_draw, 0.975)),
        near_floor_pct = mean(prob_draw < (dgp$chance + 0.03)) * 100,
        near_ceil_pct = mean(prob_draw > 0.97) * 100,
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, cell_summary)
  })
}

#' Inspect and Plot Prior-Predictive Diagnostics of a Calibrated DGP
#'
#' Displays how marginal cell probabilities map onto conditional fixed effects,
#' quantifies non-collapsibility attenuation relative to naive logit conversion,
#' and optionally plots the prior-predictive trial probability distributions
#' across experimental cells.
#'
#' @param dgp A `"glmmplanr_dgp"` object (or `"glmmplanr_spec"` / `"glmmplanr_result"`).
#' @param plot Logical; if `TRUE` (default), includes a plot showing
#'   cell-level probability distributions and guessing floor.
#' @param n_draws Number of Monte Carlo draws per cell for density plotting
#'   (default `3000L`).
#' @param seed Integer seed for the plot draws. Use `NULL` to draw a master seed
#'   from the caller's random stream. Ignored when `plot = FALSE`.
#' @inheritParams plot.glmmplanr_result
#' @return A visible S3 list of class `glmmplanr_dgp_check` with the data frames
#'   `coefficients`, `cell_diagnostics`, `marginal_effects`, and `plot` (a ggplot
#'   or `NULL`). `marginal_effects` gives each non-intercept model-matrix
#'   contrast of the integrated cell probabilities in percentage points, using
#'   the design coding. With two binary factors coded -0.5/+0.5 the interaction
#'   is the full difference in differences.
#'   Its print method displays the tables and, interactively, the plot.
#'   In `cell_diagnostics`, `quad_prob` averages over the random effects and
#'   `conditional_prob` gives the response probability with all random effects
#'   set to zero. Both include the specified guessing floor.
#'   The plot uses [ggplot2::theme_bw()]. Add [ggplot2::labs()] to supply
#'   your own title or caption.
#' @examples
#' des <- design_crossed(list(condition = c("a", "b")))
#' dgp <- calibrate_dgp(des, c(a = 0.4, b = 0.6), subj_sd = 0.3, item_sd = 0.4)
#' check_dgp(dgp, plot = FALSE)$cell_diagnostics[c("cell", "quad_prob", "conditional_prob")]
#' @export
check_dgp <- function(dgp, plot = TRUE, n_draws = 3000L, seed = 42L,
                      show_annotations = FALSE) {
  if (isTRUE(plot) && uses_box_muller()) {
    if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1L)
    return(isolated_seeded_call("check_dgp", environment(), match.call()))
  }
  if (inherits(dgp, "glmmplanr_result")) dgp <- dgp$spec$dgp
  if (inherits(dgp, "glmmplanr_spec")) dgp <- dgp$dgp
  if (!inherits(dgp, "glmmplanr_dgp")) {
    cli::cli_abort("{.arg dgp} must be a {.cls glmmplanr_dgp}, {.cls glmmplanr_spec}, or {.cls glmmplanr_result} object.")
  }

  coef_df <- data.frame(
    term = names(dgp$beta),
    calibrated_beta = unname(dgp$beta),
    naive_mafc_beta = unname(dgp$beta_naive_mafc),
    naive_std_logit = unname(dgp$beta_naive_logit),
    attenuation_ratio = ifelse(
      abs(dgp$beta_naive_mafc) > 1e-8,
      unname(dgp$beta / dgp$beta_naive_mafc),
      NA_real_
    ),
    stringsAsFactors = FALSE
  )

  p_out <- NULL
  if (isTRUE(plot)) {
    streams <- make_lecuyer_streams(seed, 1L)
    draw_df <- with_rng_stream(streams[[1L]], {
      des <- dgp$design
      p_s <- ncol(des$Z_subj_cells)
      p_i <- ncol(des$Z_item_cells)
      u_s <- matrix(stats::rnorm(n_draws * p_s), nrow = n_draws, ncol = p_s) %*% t(dgp$L_subj)
      u_i <- matrix(stats::rnorm(n_draws * p_i), nrow = n_draws, ncol = p_i) %*% t(dgp$L_item)
      u_t <- if (des$has_testlet) {
        p_t <- ncol(des$Z_testlet_cells)
        matrix(stats::rnorm(n_draws * p_t), nrow = n_draws, ncol = p_t) %*% t(dgp$L_testlet)
      } else {
        NULL
      }
      local_effect <- draw_testlet_intercepts(dgp, n_draws)
      dfs <- lapply(seq_len(des$n_cells), function(c_idx) {
        eta <- dgp$mu_cells[c_idx] +
          as.numeric(u_s %*% des$Z_subj_cells[c_idx, ]) +
          as.numeric(u_i %*% des$Z_item_cells[c_idx, ])
        if (des$has_testlet) {
          eta <- eta + as.numeric(u_t %*% des$Z_testlet_cells[c_idx, ])
        }
        data.frame(
          cell = factor(des$cells$cell_name[c_idx], levels = des$cells$cell_name),
          prob = mafc_plogis(eta + local_effect, chance = dgp$chance),
          target = unname(dgp$target_cell_probs[c_idx])
        )
      })
      do.call(rbind, dfs)
    })

    targets_df <- data.frame(
      cell = factor(dgp$design$cells$cell_name, levels = dgp$design$cells$cell_name),
      target = unname(dgp$target_cell_probs)
    )

    p_out <- ggplot2::ggplot(draw_df, ggplot2::aes(x = .data$prob, fill = .data$cell)) +
      ggplot2::geom_density(alpha = 0.55, color = "#1f2937", linewidth = 0.4) +
      ggplot2::geom_vline(
        data = targets_df,
        ggplot2::aes(xintercept = .data$target),
        linetype = "dashed",
        linewidth = 0.7,
        color = "#111827"
      ) +
      ggplot2::facet_wrap(~ cell, ncol = 2L) +
      ggplot2::scale_x_continuous(limits = c(0, 1), labels = function(x) paste0(round(x * 100), "%")) +
      ggplot2::labs(
        title = "Prior-Predictive Trial Accuracy Distribution by Condition Cell",
        subtitle = sprintf(
          "Dashed line = target marginal accuracy | Lower guessing asymptote g = %.3f",
          dgp$chance
        ),
        x = "Response Probability P(Y = 1 | subj, item)",
        y = "Prior-Predictive Density"
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(legend.position = "none")

    if (dgp$chance > 0) {
      p_out <- p_out + ggplot2::geom_vline(
        xintercept = dgp$chance,
        linetype = "dotted",
        color = "#dc2626",
        linewidth = 0.7
      )
    }
    p_out <- style_plot_annotations(p_out, show_annotations)
  }

  out <- list(
    coefficients = coef_df,
    marginal_effects = data.frame(term = coef_df$term[-1L],
      effects_pp = 100 * unname(qr.solve(dgp$design$X_cells,
        dgp_cell_probabilities(dgp, "marginal"))[-1L])),
    cell_diagnostics = data.frame(dgp$prior_check,
      conditional_prob = dgp_cell_probabilities(dgp, "conditional")),
    plot = p_out
  )
  class(out) <- "glmmplanr_dgp_check"
  out
}

#' @rdname print.glmmplanr_design
#' @export
print.glmmplanr_dgp_check <- function(x, ...) {
  cli::cli_h2("glmmplanr DGP Calibration & Prior-Predictive Check")
  cli::cli_h3("Conditional Fixed-Effect Coefficients (Logit Scale)")
  cf <- x$coefficients
  cf[, -1L] <- lapply(cf[, -1L], round, 4)
  print(cf, row.names = FALSE)
  cli::cli_h3("Marginal and Conditional Cell Probabilities")
  cd <- x$cell_diagnostics
  cd[, -1L] <- lapply(cd[, -1L], round, 4)
  print(cd, row.names = FALSE)
  if (NROW(x$marginal_effects)) {
    cli::cli_h3("Design contrasts in population accuracy (percentage points)")
    print(x$marginal_effects, digits = 4, row.names = FALSE)
  }
  if (!is.null(x$plot) && interactive()) {
    print(x$plot)
  }
  invisible(x)
}

testlet_interaction_sd <- function(dgp) {
  if (is.null(dgp$participant_testlet_sd)) 0 else dgp$participant_testlet_sd
}

draw_testlet_intercepts <- function(dgp, n) {
  sd <- testlet_interaction_sd(dgp)
  if (sd == 0) 0 else stats::rnorm(n, sd = sd)
}

trial_testlet_effect <- function(dgp, ns, nt, sid, tid) {
  if (testlet_interaction_sd(dgp) == 0) return(0)
  draw_testlet_intercepts(dgp, ns * nt)[sid + (tid - 1L) * ns]
}

#' Simulate One Balanced Dataset from a Calibrated DGP
#'
#' Generates binary trial-level responses from the maximal crossed (and optional
#' testlet-nested) data-generating process with optional random trial attrition.
#'
#' @param dgp A `"glmmplanr_dgp"` object.
#' @param schedule Precomputed schedule from [generate_design_schedule()], or
#'   `NULL` if `n_subj` and `n_item` are provided.
#' @param n_subj Integer number of subjects (used when `schedule` is `NULL`).
#' @param n_item Integer number of items (used when `schedule` is `NULL`).
#' @return A `data.frame` of observed trials with binary response `y` (`0L`/`1L`),
#'   cluster factors `subj` and `item` (plus `testlet` if present), experimental
#'   factors with the selected contrasts attached, and cell index `cell_id`.
#' @examples
#' des <- design_crossed(list(condition = c("a", "b")))
#' dgp <- calibrate_dgp(des, c(a = 0.4, b = 0.6), subj_sd = 0.3, item_sd = 0.4)
#' trials <- simulate_dgp_data(dgp, n_subj = 8, n_item = 4)
#' head(trials[c("subj", "item", "condition", "y")])
#' @export
simulate_dgp_data <- function(dgp, schedule = NULL, n_subj = NULL, n_item = NULL) {
  if (!inherits(dgp, "glmmplanr_dgp")) {
    cli::cli_abort("{.arg dgp} must be a {.cls glmmplanr_dgp} object.")
  }
  des <- dgp$design
  if (is.null(schedule)) {
    if (is.null(n_subj) || is.null(n_item)) {
      cli::cli_abort("Supply either {.arg schedule} or both {.arg n_subj} and {.arg n_item}.")
    }
    schedule <- generate_design_schedule(des, n_subj = n_subj, n_item = n_item)
  }
  draw_dgp_trials(dgp, schedule)
}

draw_dgp_trials <- function(dgp, schedule, keep_person_item_effects = FALSE) {
  des <- dgp$design
  ns <- schedule$n_subj
  ni <- schedule$n_item
  p_s <- ncol(des$Z_subj_cells)
  p_i <- ncol(des$Z_item_cells)

  u_s <- matrix(stats::rnorm(ns * p_s), nrow = ns, ncol = p_s) %*% t(dgp$L_subj)
  u_i <- matrix(stats::rnorm(ni * p_i), nrow = ni, ncol = p_i) %*% t(dgp$L_item)

  tr <- schedule$trials
  cid <- tr$cell_id
  sid <- tr$subj_id
  iid <- tr$item_id

  eta <- if (p_s == 1L && p_i == 1L &&
             all(des$Z_subj_cells == 1.0) && all(des$Z_item_cells == 1.0)) {
    dgp$mu_cells[cid] + u_s[sid, 1L] + u_i[iid, 1L]
  } else {
    zs_trial <- des$Z_subj_cells[cid, , drop = FALSE]
    zi_trial <- des$Z_item_cells[cid, , drop = FALSE]
    dgp$mu_cells[cid] +
      rowSums(zs_trial * u_s[sid, , drop = FALSE]) +
      rowSums(zi_trial * u_i[iid, , drop = FALSE])
  }

  if (des$has_testlet) {
    nt <- schedule$n_testlet
    p_t <- ncol(des$Z_testlet_cells)
    u_t <- matrix(stats::rnorm(nt * p_t), nrow = nt, ncol = p_t) %*% t(dgp$L_testlet)
    zt_trial <- des$Z_testlet_cells[cid, , drop = FALSE]
    eta <- eta + rowSums(zt_trial * u_t[tr$testlet_id, , drop = FALSE]) +
      trial_testlet_effect(dgp, ns, nt, sid, tr$testlet_id)
  }

  prob <- mafc_plogis(eta, chance = dgp$chance)
  y <- stats::rbinom(length(prob), size = 1L, prob = prob)
  tr$y <- y

  keep <- observed_trials(des, sid, ns)
  if (!is.null(keep)) {
    tr <- tr[keep, , drop = FALSE]
    rownames(tr) <- NULL
  }

  # Re-attach contrasts to factor columns explicitly
  for (fn in des$factor_names) {
    stats::contrasts(tr[[fn]]) <- des$contrasts[[fn]]
  }
  if (keep_person_item_effects) return(list(trials = tr, subj = u_s, item = u_i))
  tr
}

#' Create Sensitivity Variants of a Calibrated DGP (Optimistic / Standard / Conservative)
#'
#' Scales random slope SDs by multipliers (e.g., `0.7`, `1.0`, `1.4`) and
#' recalibrates conditional fixed effects so marginal cell probabilities remain
#' invariant across sensitivity levels.
#'
#' @param dgp A `"glmmplanr_dgp"` object.
#' @param slope_multipliers Named numeric vector of multipliers applied to
#'   `slope_sd_ratio` (default `c(optimistic = 0.70, standard = 1.00, conservative = 1.35)`).
#' @return A named list of `"glmmplanr_dgp"` objects.
#' @examples
#' des <- design_crossed(list(condition = c("a", "b")))
#' dgp <- calibrate_dgp(des, c(a = 0.4, b = 0.6), subj_sd = 0.3, item_sd = 0.4,
#'                      slope_sd_ratio = 0.35)
#' models <- sensitivity_dgps(dgp)
#' vapply(models, function(x) unname(x$tau_subj[2L]), numeric(1))
#' @export
sensitivity_dgps <- function(
  dgp,
  slope_multipliers = c(optimistic = 0.70, standard = 1.00, conservative = 1.35)
) {
  if (!inherits(dgp, "glmmplanr_dgp")) {
    cli::cli_abort("{.arg dgp} must be a {.cls glmmplanr_dgp} object.")
  }
  lapply(slope_multipliers, function(mult) {
    s_sd <- dgp$tau_subj
    i_sd <- dgp$tau_item
    if (length(s_sd) > 1L) s_sd[-1L] <- s_sd[-1L] * mult
    if (length(i_sd) > 1L) i_sd[-1L] <- i_sd[-1L] * mult
    t_sd <- if (dgp$design$has_testlet) {
      ts <- dgp$tau_testlet
      if (length(ts) > 1L) ts[-1L] <- ts[-1L] * mult
      ts
    } else {
      NULL
    }
    calibrate_dgp(
      design = dgp$design,
      cell_probs = dgp$target_cell_probs,
      chance = dgp$chance,
      subj_sd = s_sd,
      item_sd = i_sd,
      corr_structure = dgp$corr_structure,
      control = calibrate_control(
        testlet_sd = t_sd,
        participant_testlet_sd = testlet_interaction_sd(dgp),
        n_quad = dgp$n_quad,
        check_nsim = 1000L
      )
    )
  })
}

#' @rdname print.glmmplanr_design
#' @export
print.glmmplanr_dgp <- function(x, ...) {
  cli::cli_h2("glmmplanr Calibrated Data-Generating Process (DGP)")
  link_line <- if (x$chance == 0) {
    "Response link: ordinary binary logit (chance = 0)"
  } else {
    sprintf("Response link: guessing-adjusted logit (chance = %s)", round(x$chance, 4))
  }
  mu_mid <- mean(x$mu_cells)
  s_span <- sd_to_prob_range(x$tau_subj[1L], intercept = mu_mid, chance = x$chance)
  i_span <- sd_to_prob_range(x$tau_item[1L], intercept = mu_mid, chance = x$chance)
  has_within_factors <- length(x$design$factor_names) > 0L &&
    (ncol(x$design$Z_subj_cells) > 1L || ncol(x$design$Z_item_cells) > 1L)
  no_slopes <- (length(x$tau_subj) <= 1L || all(abs(x$tau_subj[-1L]) < 1e-12)) &&
    (length(x$tau_item) <= 1L || all(abs(x$tau_item[-1L]) < 1e-12))
  cli::cli_ul(c(
    link_line,
    "Subject RE SDs: {paste(sprintf('%s=%.3f', names(x$tau_subj), x$tau_subj), collapse = ', ')}",
    "Item RE SDs: {paste(sprintf('%s=%.3f', names(x$tau_item), x$tau_item), collapse = ', ')}",
    sprintf(
      "Implied 95%% accuracy spans (at mean cell): participants = [%.0f%%, %.0f%%], items = [%.0f%%, %.0f%%]",
      100 * s_span[1L], 100 * s_span[2L], 100 * i_span[1L], 100 * i_span[2L]
    ),
    if (has_within_factors && no_slopes) {
      "Random-slope assumption: 0 (assumes identical condition effects across all participants and items; set slope_sd_ratio = 0.25-0.35 to allow heterogeneity)"
    },
    if (x$design$has_testlet) {
      "Testlet RE SDs: {paste(sprintf('%s=%.3f', names(x$tau_testlet), x$tau_testlet), collapse = ', ')}"
    },
    "Gauss-Hermite inversion ({x$n_quad} nodes, exact 1D variance-sum reduction): max |residual| = { format(x$solver_info$max_abs_residual, scientific = TRUE, digits = 2) } ({x$solver_info$iterations} Newton steps)"
  ))
  coef_tbl <- data.frame(
    Term = names(x$beta),
    Calibrated_Beta = round(unname(x$beta), 4),
    Naive_mAFC_Beta = round(unname(x$beta_naive_mafc), 4),
    Naive_Std_Logit = round(unname(x$beta_naive_logit), 4)
  )
  print(coef_tbl, row.names = FALSE)
  invisible(x)
}
