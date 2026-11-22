#' Extract Cell-Probability Contrast Weights for a Fixed Effect
#'
#' Maps a fixed-effect model-matrix column `effect_col` onto its
#' exact linear contrast weights \eqn{w \in \mathbb{R}^M} across the `M`
#' experimental condition cells via the Moore-Penrose / least-squares inverse of
#' `design$X_cells`.
#'
#' For example, in a \eqn{2 \times 2} centered sum-coded design:
#' * Main effect of factor A: \eqn{w = (-0.5, -0.5, +0.5, +0.5)} (marginal mean
#'   difference \eqn{\bar{P}_{A2} - \bar{P}_{A1}}),
#' * \eqn{2 \times 2} interaction `A:B`: \eqn{w = (+1, -1, -1, +1)}
#'   (difference-in-differences \eqn{(P_{22} - P_{21}) - (P_{12} - P_{11})}).
#'
#' @param design A `"glmmplanr_design"` object.
#' @param effect_col Model-matrix column name (resolved via `resolve_effect_column`).
#' @param custom_weights Optional numeric vector of length `design$n_cells` (named
#'   or unnamed) overriding automatic contrast extraction.
#' @return Named numeric vector of cell contrast weights of length `design$n_cells`.
#' @export
cell_contrast_weights <- function(design, effect_col, custom_weights = NULL) {
  if (!inherits(design, "glmmplanr_design")) {
    cli::cli_abort("{.arg design} must be a {.cls glmmplanr_design} object.")
  }
  cell_nms <- design$cells$cell_name
  if (!is.null(custom_weights)) {
    if (!is.numeric(custom_weights) || length(custom_weights) != length(cell_nms) ||
        anyNA(custom_weights)) {
      cli::cli_abort("{.arg custom_weights} must be a numeric vector of length {length(cell_nms)}.")
    }
    if (!is.null(names(custom_weights)) && all(nzchar(names(custom_weights)))) {
      norm_nms <- gsub("[:./\\s]+", "_", trimws(names(custom_weights)), perl = TRUE)
      idx <- match(cell_nms, norm_nms)
      if (anyNA(idx)) {
        cli::cli_abort("Names of {.arg custom_weights} must match {.val {cell_nms}}.")
      }
      return(stats::setNames(as.numeric(custom_weights[idx]), cell_nms))
    }
    return(stats::setNames(as.numeric(custom_weights), cell_nms))
  }

  X <- design$X_cells
  if (!(effect_col %in% colnames(X))) {
    effect_col <- resolve_effect_column(design, effect_col)
  }
  proj <- qr.solve(crossprod(X), t(X))
  w <- as.numeric(proj[effect_col, ])
  stats::setNames(w, cell_nms)
}

#' Compute True Population-Marginal or Conditional Probability Contrast from DGP
#'
#' @param dgp A `"glmmplanr_dgp"` object.
#' @param weights Numeric vector of cell contrast weights from [cell_contrast_weights()].
#' @param estimand Either `"marginal"` (default; integrates out subject, item,
#'   and optional testlet random effects via Gauss-Hermite quadrature) or
#'   `"conditional"` (evaluates at average cluster random effects \eqn{u = 0}).
#' @return Numeric scalar giving the true probability contrast \eqn{\Delta P \in (-1, 1)}
#'   (multiply by `100` for percentage points).
#' @export
true_prob_effect <- function(
  dgp,
  weights,
  estimand = c("marginal", "conditional")
) {
  sum(weights * dgp_cell_probabilities(dgp, match.arg(estimand)))
}

dgp_cell_probabilities <- function(dgp, estimand) {
  des <- dgp$design
  mu_cell <- as.numeric(des$X_cells %*% dgp$beta)

  if (estimand == "conditional") {
    p_cell <- mafc_plogis(mu_cell, chance = dgp$chance)
    return(p_cell)
  }

  var_s <- rowSums((des$Z_subj_cells %*% dgp$Sigma_subj) * des$Z_subj_cells)
  var_i <- rowSums((des$Z_item_cells %*% dgp$Sigma_item) * des$Z_item_cells)
  var_t <- if (des$has_testlet && !is.null(dgp$Sigma_testlet)) {
    rowSums((des$Z_testlet_cells %*% dgp$Sigma_testlet) * des$Z_testlet_cells)
  } else {
    rep(0.0, des$n_cells)
  }
  sd_tot <- sqrt(pmax(var_s + var_i + var_t + testlet_interaction_sd(dgp)^2, 0.0))

  q1 <- gh_quad_1d(dgp$n_quad)
  linpred <- outer(mu_cell, rep.int(1.0, length(q1$weights))) + outer(sd_tot, q1$nodes)
  p_inner <- stats::plogis(linpred)
  p_marg <- dgp$chance + (1.0 - dgp$chance) * as.numeric(p_inner %*% q1$weights)
  p_marg
}

#' Estimate Marginal or Conditional Probability Contrast & Joint Delta-Method SE
#'
#' Performs population standardization (G-computation via Gauss-Hermite
#' quadrature over estimated random-effect variances) or conditional evaluation
#' (\eqn{u = 0}) and computes the multivariate delta-method standard error.
#' When variance-component log-SD variances (`var_log_tau_*`) are supplied,
#' propagates joint uncertainty in both fixed effects \eqn{\beta} and
#' random-effect variance components \eqn{\theta = \log \tau} analytically,
#' accounting for non-collapsibility cross-covariance via `J_beta_theta`.
#'
#' @keywords internal
estimate_prob_effect_delta <- function(
  beta,
  vcov_beta,
  Sigma_subj,
  Sigma_item,
  Sigma_testlet,
  design,
  chance,
  weights,
  estimand = c("marginal", "conditional"),
  re_structure = "maximal",
  n_quad = 21L,
  var_log_tau_subj = NULL,
  var_log_tau_item = NULL,
  var_log_tau_testlet = NULL,
  J_beta_theta = NULL
) {
  estimand <- match.arg(estimand)
  X_c <- design$X_cells
  mu_c <- as.numeric(X_c %*% beta)

  if (estimand == "conditional") {
    p_in <- stats::plogis(mu_c)
    p_c <- chance + (1.0 - chance) * p_in
    d_c <- (1.0 - chance) * p_in * (1.0 - p_in)
    est_diff <- sum(weights * p_c)
    grad_beta <- as.numeric(crossprod(X_c, weights * d_c))
    var_diff <- as.numeric(crossprod(grad_beta, vcov_beta %*% grad_beta))
    return(list(est = est_diff, se = sqrt(max(var_diff, 1e-16))))
  }

  Z_s_c <- if (re_structure == "intercept_only") {
    matrix(1.0, nrow = design$n_cells, ncol = 1L)
  } else {
    design$Z_subj_cells
  }
  Z_i_c <- if (re_structure == "intercept_only") {
    matrix(1.0, nrow = design$n_cells, ncol = 1L)
  } else {
    design$Z_item_cells
  }
  Z_t_c <- if (design$has_testlet && !is.null(Sigma_testlet)) {
    if (re_structure == "intercept_only") {
      matrix(1.0, nrow = design$n_cells, ncol = 1L)
    } else {
      design$Z_testlet_cells
    }
  } else {
    NULL
  }

  ZS_s <- if (!is.null(Sigma_subj)) Z_s_c %*% Sigma_subj else NULL
  ZS_i <- if (!is.null(Sigma_item)) Z_i_c %*% Sigma_item else NULL
  ZS_t <- if (!is.null(Z_t_c)) Z_t_c %*% Sigma_testlet else NULL

  var_s <- if (!is.null(ZS_s)) rowSums(ZS_s * Z_s_c) else rep(0.0, design$n_cells)
  var_i <- if (!is.null(ZS_i)) rowSums(ZS_i * Z_i_c) else rep(0.0, design$n_cells)
  var_t <- if (!is.null(ZS_t)) rowSums(ZS_t * Z_t_c) else rep(0.0, design$n_cells)

  sd_tot <- sqrt(pmax(var_s + var_i + var_t, 0.0))
  q1 <- gh_quad_1d(n_quad)
  linpred <- outer(mu_c, rep.int(1.0, length(q1$weights))) + outer(sd_tot, q1$nodes)
  p_in <- stats::plogis(linpred)
  p_c <- chance + (1.0 - chance) * as.numeric(p_in %*% q1$weights)
  deriv_in <- (1.0 - chance) * p_in * (1.0 - p_in)
  d_c <- as.numeric(deriv_in %*% q1$weights)

  est_diff <- sum(weights * p_c)
  grad_beta <- as.numeric(crossprod(X_c, weights * d_c))
  var_diff <- as.numeric(crossprod(grad_beta, vcov_beta %*% grad_beta))

  # Joint delta propagation wrt random-effect log-SDs theta = log(tau)
  if (!is.null(var_log_tau_subj) && !is.null(var_log_tau_item)) {
    dP_dsd <- as.numeric(deriv_in %*% (q1$weights * q1$nodes))
    inv_sd <- ifelse(sd_tot > 1e-8, 1.0 / sd_tot, 0.0)
    w_dsd <- weights * dP_dsd * inv_sd

    g_theta_s <- if (!is.null(ZS_s)) as.numeric(crossprod(Z_s_c * ZS_s, w_dsd)) else numeric(0)
    g_theta_i <- if (!is.null(ZS_i)) as.numeric(crossprod(Z_i_c * ZS_i, w_dsd)) else numeric(0)
    g_theta_t <- if (!is.null(ZS_t) && !is.null(var_log_tau_testlet)) {
      as.numeric(crossprod(Z_t_c * ZS_t, w_dsd))
    } else {
      numeric(0)
    }

    g_theta <- c(g_theta_s, g_theta_i, g_theta_t)
    v_theta <- c(
      as.numeric(var_log_tau_subj),
      as.numeric(var_log_tau_item),
      if (!is.null(var_log_tau_testlet)) as.numeric(var_log_tau_testlet) else numeric(0)
    )

    if (length(g_theta) == length(v_theta) && length(v_theta) > 0L) {
      if (!is.null(J_beta_theta) && is.matrix(J_beta_theta) &&
          nrow(J_beta_theta) == length(grad_beta) && ncol(J_beta_theta) == length(v_theta)) {
        # vcov_beta = vcov_fixed + J * diag(v_theta) * J'
        # Total variance = g_b' vcov_fixed g_b + (J' g_b + g_theta)' diag(v_theta) (J' g_b + g_theta)
        a_vec <- as.numeric(crossprod(J_beta_theta, grad_beta))
        cross_extra <- sum(v_theta * ((a_vec + g_theta)^2 - (a_vec^2)))
        # Ensure total variance stays >= fixed-theta delta floor
        fixed_floor <- max(var_diff - sum(v_theta * (a_vec^2)), 0.0)
        var_diff <- max(var_diff + cross_extra, fixed_floor)
      } else {
        var_diff <- var_diff + sum(v_theta * (g_theta^2))
      }
    }
  }

  se_diff <- sqrt(max(var_diff, 1e-16))
  list(est = est_diff, se = se_diff)
}

#' Construct Cell Probabilities Directly from Percentage-Point Effects
#'
#' Builds a named vector of marginal `cell_probs` from a grand-mean (or
#' reference-cell) accuracy and substantive percentage-point effect sizes (e.g.,
#' `+8` points for exposure, `+5` points for spacing, `+5` points for their
#' interaction).
#'
#' @param design A `"glmmplanr_design"` object.
#' @param baseline_prob Numeric scalar in `(0, 1)`: grand-mean accuracy across
#'   cells when `design$contrast != "treatment"`, or reference-cell accuracy when
#'   `design$contrast == "treatment"`.
#' @param effects_pp Named numeric vector of marginal accuracy differences in
#'   **percentage points** (e.g., `c(exposure = 8, spacing = 5, "exposure:spacing" = 5)`).
#'   Omitted terms default to `0`.
#' @return Named numeric vector of cell probabilities ready for [calibrate_dgp()].
#' @examples
#' des <- design_crossed(list(
#'   method = c("control", "treatment")
#' ), allocation = "between_subj")
#' cell_probs_from_diffs(des, baseline_prob = 0.65, effects_pp = c(method = 8))
#' @export
cell_probs_from_diffs <- function(design, baseline_prob, effects_pp = numeric(0)) {
  if (!inherits(design, "glmmplanr_design")) {
    cli::cli_abort("{.arg design} must be a {.cls glmmplanr_design} object.")
  }
  if (!is.numeric(baseline_prob) || length(baseline_prob) != 1L ||
      is.na(baseline_prob) || baseline_prob <= 0 || baseline_prob >= 1) {
    cli::cli_abort("{.arg baseline_prob} must be a probability strictly in (0, 1).")
  }
  X <- design$X_cells
  b_prob <- stats::setNames(rep(0.0, ncol(X)), colnames(X))
  b_prob[1L] <- baseline_prob

  if (length(effects_pp) > 0L) {
    if (is.null(names(effects_pp)) || any(!nzchar(names(effects_pp)))) {
      cli::cli_abort("{.arg effects_pp} must be a named numeric vector of percentage-point differences.")
    }
    for (nm in names(effects_pp)) {
      col_nm <- resolve_effect_column(design, nm)
      b_prob[[col_nm]] <- as.numeric(effects_pp[[nm]]) / 100.0
    }
  }

  probs <- stats::setNames(as.numeric(X %*% b_prob), design$cells$cell_name)
  if (any(probs <= 0) || any(probs >= 1)) {
    cli::cli_abort(c(
      "Implied marginal cell probabilities fall outside (0, 1): {paste(sprintf('%s=%.3f', names(probs), probs), collapse=', ')}.",
      "i" = "Adjust {.arg baseline_prob} or reduce {.arg effects_pp}."
    ))
  }
  probs
}

#' Build Multi-Scenario DGPs Varying Percentage-Point Effect Sizes and/or Random Slopes
#'
#' Creates a named list of calibrated `"glmmplanr_dgp"` scenarios (for example,
#' conservative 5-point, standard 8-point, and optimistic 10-point improvements)
#' that can be passed directly to `optimize_sample_size(..., scenarios = ...)`
#' to find both the standard-assumption optimum and the worst-case robust
#' optimum along with the additional sample size required for robustness.
#'
#' @param base_dgp A calibrated `"glmmplanr_dgp"` object serving as the baseline
#'   reference template.
#' @param effect Focal effect name whose marginal percentage-point difference is
#'   varied across scenarios.
#' @param effect_diffs_pp Named numeric vector of marginal effect sizes in
#'   **percentage points** (e.g., `c(diff_5pp = 5, diff_8pp = 8, diff_10pp = 10)`).
#' @param slope_multipliers Optional numeric vector (length 1 or matching
#'   `effect_diffs_pp`) scaling random-slope SDs in each scenario (default `1.0`).
#' @return A named list of `"glmmplanr_dgp"` objects.
#' @export
scenario_dgps <- function(
  base_dgp,
  effect,
  effect_diffs_pp,
  slope_multipliers = 1.0
) {
  if (!inherits(base_dgp, "glmmplanr_dgp")) {
    cli::cli_abort("{.arg base_dgp} must be a {.cls glmmplanr_dgp} object.")
  }
  if (!is.numeric(effect_diffs_pp) || length(effect_diffs_pp) < 1L || anyNA(effect_diffs_pp)) {
    cli::cli_abort("{.arg effect_diffs_pp} must be a non-empty numeric vector of percentage-point differences.")
  }
  if (is.null(names(effect_diffs_pp)) || any(!nzchar(names(effect_diffs_pp)))) {
    names(effect_diffs_pp) <- paste0("diff_", effect_diffs_pp, "pp")
  }
  if (length(slope_multipliers) == 1L) {
    slope_multipliers <- rep(as.numeric(slope_multipliers), length(effect_diffs_pp))
  } else if (length(slope_multipliers) != length(effect_diffs_pp)) {
    cli::cli_abort("{.arg slope_multipliers} must have length 1 or match {.arg effect_diffs_pp}.")
  }

  des <- base_dgp$design
  eff_col <- resolve_effect_column(des, effect)
  X <- des$X_cells
  prob_coefs <- qr.solve(X, base_dgp$target_cell_probs)

  out <- vector("list", length(effect_diffs_pp))
  names(out) <- names(effect_diffs_pp)

  for (k in seq_along(effect_diffs_pp)) {
    pc_k <- prob_coefs
    pc_k[eff_col] <- as.numeric(effect_diffs_pp[k]) / 100.0
    new_probs <- stats::setNames(as.numeric(X %*% pc_k), des$cells$cell_name)
    if (any(new_probs <= base_dgp$chance) || any(new_probs >= 1)) {
      cli::cli_abort(
        "Scenario {.val {names(effect_diffs_pp)[k]}} ({effect_diffs_pp[k]} pp) produces cell probabilities outside (chance = {round(base_dgp$chance, 3)}, 1): {paste(round(new_probs, 3), collapse = ', ')}."
      )
    }
    mult <- slope_multipliers[k]
    s_sd <- base_dgp$tau_subj
    i_sd <- base_dgp$tau_item
    if (length(s_sd) > 1L) s_sd[-1L] <- s_sd[-1L] * mult
    if (length(i_sd) > 1L) i_sd[-1L] <- i_sd[-1L] * mult
    t_sd <- if (des$has_testlet) {
      ts <- base_dgp$tau_testlet
      if (length(ts) > 1L) ts[-1L] <- ts[-1L] * mult
      ts
    } else {
      NULL
    }

    out[[k]] <- calibrate_dgp(
      design = des,
      cell_probs = new_probs,
      chance = base_dgp$chance,
      subj_sd = s_sd,
      item_sd = i_sd,
      testlet_sd = t_sd,
      participant_testlet_sd = testlet_interaction_sd(base_dgp),
      corr_structure = base_dgp$corr_structure,
      n_quad = base_dgp$n_quad,
      check_nsim = 1000L
    )
  }
  out
}
