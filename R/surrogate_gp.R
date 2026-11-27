#' Space-Filling Balanced Initial Design (Maximin LHS + Corner Anchors)
#'
#' Selects `n_init` distinct balanced grid points `(n_subj, n_item)` combining
#' corner anchors with maximin Latin Hypercube Sampling so the Gaussian process
#' surrogate covers the full 2D search domain without extrapolation.
#'
#' @noRd
select_initial_lhs_points <- function(grid_df, n_init = 16L, seed = 1L) {
  M <- nrow(grid_df)
  n_init <- min(as.integer(n_init), M)
  if (n_init >= M) {
    return(seq_len(M))
  }

  streams <- make_lecuyer_streams(seed, 1L)
  with_rng_stream(streams[[1L]], {
    u1 <- normalize_unit(grid_df$n_subj)
    u2 <- normalize_unit(grid_df$n_item)
    U <- cbind(u1, u2)

    # Always include the 4 corners when n_init >= 6 to stabilize boundaries
    s_min <- min(grid_df$n_subj); s_max <- max(grid_df$n_subj)
    i_min <- min(grid_df$n_item); i_max <- max(grid_df$n_item)
    corners <- unique(c(
      which(grid_df$n_subj == s_min & grid_df$n_item == i_min),
      which(grid_df$n_subj == s_min & grid_df$n_item == i_max),
      which(grid_df$n_subj == s_max & grid_df$n_item == i_min),
      which(grid_df$n_subj == s_max & grid_df$n_item == i_max)
    ))

    if (n_init <= length(corners)) {
      return(corners[seq_len(n_init)])
    }

    n_rem <- n_init - length(corners)
    best_idx <- NULL
    best_min_dist <- -Inf

    for (cand in seq_len(80L)) {
      pts1 <- (sample.int(n_rem) - stats::runif(n_rem)) / n_rem
      pts2 <- (sample.int(n_rem) - stats::runif(n_rem)) / n_rem
      chosen <- corners
      avail <- setdiff(seq_len(M), chosen)
      for (k in seq_len(n_rem)) {
        if (length(avail) == 0L) break
        d2 <- (U[avail, 1L] - pts1[k])^2 + (U[avail, 2L] - pts2[k])^2
        pick <- avail[which.min(d2)]
        chosen <- c(chosen, pick)
        avail <- setdiff(avail, pick)
      }
      if (length(chosen) < n_init) next
      dmat <- as.matrix(stats::dist(U[chosen, , drop = FALSE]))
      diag(dmat) <- Inf
      md <- min(dmat)
      if (md > best_min_dist) {
        best_min_dist <- md
        best_idx <- chosen
      }
    }
    sort(best_idx)
  })
}

normalize_unit <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng)) || rng[2L] <= rng[1L]) return(rep(0.5, length(x)))
  (x - rng[1L]) / (rng[2L] - rng[1L])
}

#' Fit a Monotone Anisotropic Matern-5/2 Gaussian Process Surrogate on `[0, 1]^2`
#'
#' Combines an asymptotic GLMM inverse-sample-size prior mean trend with an
#' anisotropic Matern-5/2 Universal Kriging Gaussian process accounting for both
#' Monte Carlo binomial noise and trend-parameter estimation uncertainty,
#' followed by Dykstra 2D isotonic projection over `(n_subj, n_item)`.
#'
#' @noRd
fit_monotone_gp_prob <- function(
  eval_df,
  grid_df,
  prob_col = "power",
  nsim_col = "n_valid"
) {
  if (!(nsim_col %in% names(eval_df)) && ("nsim" %in% names(eval_df))) {
    nsim_col <- "nsim"
  }
  valid <- is.finite(eval_df[[prob_col]]) & !is.null(eval_df[[nsim_col]]) & (eval_df[[nsim_col]] > 0)
  ev <- eval_df[valid, , drop = FALSE]
  n_obs <- nrow(ev)

  s_vals <- sort(unique(grid_df$n_subj))
  i_vals <- sort(unique(grid_df$n_item))

  if (n_obs == 0L) {
    return(list(
      mean_prob = rep(NA_real_, nrow(grid_df)),
      sd_probit = rep(1.0, nrow(grid_df)),
      mean_probit = rep(0.0, nrow(grid_df)),
      n_ref = NA_real_
    ))
  }
  if (n_obs < 3L) {
    p_bar <- mean(ev[[prob_col]])
    n_med <- stats::median(pmax(ev[[nsim_col]], 5))
    return(list(
      mean_prob = rep(p_bar, nrow(grid_df)),
      sd_probit = rep(1.0, nrow(grid_df)),
      mean_probit = rep(stats::qnorm(pmin(pmax(p_bar, 0.01), 0.99)), nrow(grid_df)),
      n_ref = n_med
    ))
  }

  # Laplace-smoothed probability -> probit space + delta-method noise variance
  n_i <- pmax(ev[[nsim_col]], 5)
  k_i <- ev[[prob_col]] * n_i
  p_tilde <- (k_i + 0.75) / (n_i + 1.5)
  p_tilde <- pmin(pmax(p_tilde, 0.005), 0.995)
  y_probit <- stats::qnorm(p_tilde)
  phi_y <- pmax(stats::dnorm(y_probit), 0.02)
  noise_var <- (p_tilde * (1 - p_tilde)) / (n_i * (phi_y^2))
  noise_var <- pmin(pmax(noise_var, 1e-4), 2.5)

  # Monotone asymptotic GLMM basis: [1, -1/N_s, -1/N_i, -1/(N_s*N_i)]
  H_obs <- cbind(1, -1 / ev$n_subj, -1 / ev$n_item, -1 / (ev$n_subj * ev$n_item))
  H_grid <- cbind(
    1,
    -1 / grid_df$n_subj,
    -1 / grid_df$n_item,
    -1 / (grid_df$n_subj * grid_df$n_item)
  )
  trend_coef <- fit_nonneg_slope_wls(H_obs, y_probit, w = 1 / noise_var)
  m_obs <- as.numeric(H_obs %*% trend_coef)
  m_grid <- as.numeric(H_grid %*% trend_coef)
  resid_obs <- y_probit - m_obs

  # Fit Matern-5/2 GP on residuals with Universal Kriging trend-uncertainty propagation
  U_all <- cbind(normalize_unit(grid_df$n_subj), normalize_unit(grid_df$n_item))
  obs_match <- match(
    paste(ev$n_subj, ev$n_item, sep = "_"),
    paste(grid_df$n_subj, grid_df$n_item, sep = "_")
  )
  U_obs <- U_all[obs_match, , drop = FALSE]

  gp_pred <- fit_predict_matern52(
    U_obs, resid_obs, noise_var, U_all, H_obs = H_obs, H_pred = H_grid
  )
  raw_probit <- m_grid + gp_pred$mean

  # Project predicted surface onto 2D monotonically non-decreasing cone
  mono_probit <- project_2d_monotone_grid(
    vals = raw_probit,
    n_subj = grid_df$n_subj,
    n_item = grid_df$n_item,
    s_levels = s_vals,
    i_levels = i_vals,
    increasing = TRUE
  )

  list(
    mean_probit = mono_probit,
    sd_probit = pmax(gp_pred$sd, 0.015),
    mean_prob = stats::pnorm(mono_probit),
    n_ref = stats::median(n_i)
  )
}

smooth_prob_threshold <- function(p, n_ref = NA_real_) {
  if (!is.finite(n_ref) || n_ref <= 0) {
    return(min(max(p, 0.005), 0.995))
  }
  p_sm <- (p * n_ref + 0.75) / (n_ref + 1.5)
  min(max(p_sm, 0.005), 0.995)
}

#' Fit a Monotone GP Surrogate for Expected 95% CI Width
#' @noRd
fit_monotone_gp_width <- function(eval_df, grid_df, width_col = "ci_width_mean") {
  valid <- is.finite(eval_df[[width_col]]) & (eval_df[[width_col]] > 0)
  ev <- eval_df[valid, , drop = FALSE]
  s_vals <- sort(unique(grid_df$n_subj))
  i_vals <- sort(unique(grid_df$n_item))

  if (nrow(ev) == 0L) {
    return(list(mean_width = rep(NA_real_, nrow(grid_df)), sd_logw = rep(0.2, nrow(grid_df))))
  }
  if (nrow(ev) < 3L) {
    w_bar <- mean(ev[[width_col]])
    return(list(mean_width = rep(w_bar, nrow(grid_df)), sd_logw = rep(0.2, nrow(grid_df))))
  }

  y_logw <- log(ev[[width_col]])
  n_w <- if (!is.null(ev$n_valid)) ev$n_valid else ev$nsim
  noise_var <- pmax(1 / (2 * pmax(n_w, 10)), 1e-4)

  # Decreasing in (N_s, N_i) -> positive weights on (1/N_s, 1/N_i, 1/(N_s*N_i))
  H_obs <- cbind(1, 1 / ev$n_subj, 1 / ev$n_item, 1 / (ev$n_subj * ev$n_item))
  H_grid <- cbind(
    1,
    1 / grid_df$n_subj,
    1 / grid_df$n_item,
    1 / (grid_df$n_subj * grid_df$n_item)
  )
  trend_coef <- fit_nonneg_slope_wls(H_obs, y_logw, w = 1 / noise_var)
  m_obs <- as.numeric(H_obs %*% trend_coef)
  m_grid <- as.numeric(H_grid %*% trend_coef)
  resid_obs <- y_logw - m_obs

  U_all <- cbind(normalize_unit(grid_df$n_subj), normalize_unit(grid_df$n_item))
  obs_match <- match(
    paste(ev$n_subj, ev$n_item, sep = "_"),
    paste(grid_df$n_subj, grid_df$n_item, sep = "_")
  )
  U_obs <- U_all[obs_match, , drop = FALSE]

  gp_pred <- fit_predict_matern52(
    U_obs, resid_obs, noise_var, U_all, H_obs = H_obs, H_pred = H_grid
  )
  raw_logw <- m_grid + gp_pred$mean

  mono_logw <- project_2d_monotone_grid(
    vals = raw_logw,
    n_subj = grid_df$n_subj,
    n_item = grid_df$n_item,
    s_levels = s_vals,
    i_levels = i_vals,
    increasing = FALSE
  )

  list(
    mean_width = exp(mono_logw),
    sd_logw = pmax(gp_pred$sd, 0.01)
  )
}

fit_nonneg_slope_wls <- function(H, y, w) {
  p <- ncol(H)
  w_sum <- sum(w)
  if (!is.finite(w_sum) || w_sum <= 0) {
    return(c(mean(y), rep(0.0, max(0L, p - 1L))))
  }
  y_bar <- sum(w * y) / w_sum
  if (p <= 1L) {
    return(y_bar)
  }
  q <- p - 1L
  H_sl <- H[, -1L, drop = FALSE]
  h_bar <- as.numeric(crossprod(H_sl, w)) / w_sum
  H_c <- sweep(H_sl, 2L, h_bar, "-")
  y_c <- y - y_bar
  A <- crossprod(H_c, H_c * w) + diag(2e-4, nrow = q, ncol = q)
  rhs <- as.numeric(crossprod(H_c, w * y_c))

  best_b <- rep(0.0, q)
  best_val <- 0.0
  if (q <= 12L) {
    n_masks <- bitwShiftL(1L, q) - 1L
    for (mask in seq_len(n_masks)) {
      S <- which(as.logical(intToBits(mask)[seq_len(q)]))
      ch <- tryCatch(chol(A[S, S, drop = FALSE]), error = function(e) NULL)
      if (is.null(ch)) next
      b_S <- backsolve(ch, forwardsolve(t(ch), rhs[S]))
      if (any(!is.finite(b_S)) || any(b_S < -1e-10)) next
      b_S <- pmax(b_S, 0.0)
      b_cand <- rep(0.0, q)
      b_cand[S] <- b_S
      val <- 0.5 * sum(b_cand * as.numeric(A %*% b_cand)) - sum(rhs * b_cand)
      if (is.finite(val) && val < best_val) {
        best_val <- val
        best_b <- b_cand
      }
    }
  } else {
    # Projected coordinate descent fallback for high-dimensional bases
    d_A <- pmax(diag(A), 2e-4)
    for (iter in seq_len(200L)) {
      b_prev <- best_b
      for (j in seq_len(q)) {
        r_j <- rhs[j] - sum(A[j, -j] * best_b[-j])
        best_b[j] <- max(0.0, r_j / d_A[j])
      }
      if (max(abs(best_b - b_prev)) < 1e-12) break
    }
  }
  b0 <- y_bar - sum(h_bar * best_b)
  c(b0, best_b)
}

matern52_cov <- function(U1, U2, sigma_f, l1, l2) {
  d1 <- outer(U1[, 1L], U2[, 1L], "-") / l1
  d2 <- outer(U1[, 2L], U2[, 2L], "-") / l2
  r <- sqrt(d1^2 + d2^2)
  sr5 <- sqrt(5) * r
  (sigma_f^2) * (1 + sr5 + (5 / 3) * (r^2)) * exp(-sr5)
}

fit_predict_matern52 <- function(U_obs, y_resid, noise_var, U_pred, H_obs = NULL, H_pred = NULL) {
  sd_r <- stats::sd(y_resid)
  sd_emp <- max(if (is.finite(sd_r)) sd_r else 0.0, 0.08)

  neg_log_marg_lik <- function(par) {
    sf <- exp(par[1L])
    l1 <- exp(par[2L])
    l2 <- exp(par[3L])
    K <- matern52_cov(U_obs, U_obs, sf, l1, l2)
    diag(K) <- diag(K) + noise_var + 1e-6
    ch <- tryCatch(chol(K), error = function(e) NULL)
    if (is.null(ch)) return(1e10)
    alpha <- backsolve(ch, forwardsolve(t(ch), y_resid))
    logdet <- 2 * sum(log(diag(ch)))
    prior_pen <- 0.5 * ((par[2L] - log(0.45))^2 + (par[3L] - log(0.45))^2)
    val <- 0.5 * sum(y_resid * alpha) + 0.5 * logdet + prior_pen
    if (!is.finite(val)) 1e10 else val
  }

  init_par <- c(log(sd_emp), log(0.45), log(0.45))
  opt <- tryCatch(
    stats::optim(
      par = init_par,
      fn = neg_log_marg_lik,
      method = "L-BFGS-B",
      lower = c(log(0.01), log(0.12), log(0.12)),
      upper = c(log(3.0), log(3.0), log(3.0))
    ),
    error = function(e) list(par = init_par)
  )

  sf <- exp(opt$par[1L])
  l1 <- exp(opt$par[2L])
  l2 <- exp(opt$par[3L])

  K_oo <- matern52_cov(U_obs, U_obs, sf, l1, l2)
  diag(K_oo) <- diag(K_oo) + noise_var + 1e-6
  ch <- tryCatch(chol(K_oo), error = function(e) {
    diag(K_oo) <- diag(K_oo) + 1e-4
    chol(K_oo)
  })
  alpha <- backsolve(ch, forwardsolve(t(ch), y_resid))

  K_po <- matern52_cov(U_pred, U_obs, sf, l1, l2)
  pred_mean <- as.numeric(K_po %*% alpha)

  V <- forwardsolve(t(ch), t(K_po))
  pred_var <- pmax((sf^2) - colSums(V^2), 1e-6)

  # Universal Kriging trend-parameter uncertainty propagation
  if (!is.null(H_obs) && !is.null(H_pred) && ncol(H_obs) == ncol(H_pred)) {
    p_tr <- ncol(H_obs)
    F_mat <- forwardsolve(t(ch), H_obs)
    ridge_pen <- diag(c(1e-4, rep(1e-2, p_tr - 1L)), nrow = p_tr)
    A_tr <- crossprod(F_mat) + ridge_pen
    ch_tr <- tryCatch(chol(A_tr), error = function(e) NULL)
    if (!is.null(ch_tr)) {
      D_tr <- t(H_pred) - crossprod(F_mat, V)
      W_tr <- forwardsolve(t(ch_tr), D_tr)
      pred_var <- pred_var + colSums(W_tr^2)
    }
  }

  list(mean = pred_mean, sd = sqrt(pred_var), hyper = c(sigma_f = sf, l_subj = l1, l_item = l2))
}

#' 2D Isotonic Projection via Dykstra Alternating Projections + PAVA
#' @noRd
project_2d_monotone_grid <- function(
  vals,
  n_subj,
  n_item,
  s_levels,
  i_levels,
  increasing = TRUE
) {
  ns <- length(s_levels)
  ni <- length(i_levels)
  if (ns <= 1L && ni <= 1L) return(vals)

  sgn <- if (increasing) 1.0 else -1.0
  M <- matrix(NA_real_, nrow = ns, ncol = ni)
  r_idx <- match(n_subj, s_levels)
  c_idx <- match(n_item, i_levels)
  M[cbind(r_idx, c_idx)] <- sgn * vals

  P <- matrix(0, nrow = ns, ncol = ni)
  Q <- matrix(0, nrow = ns, ncol = ni)
  X <- M

  for (sweep in seq_len(12L)) {
    Y <- X + P
    if (ni > 1L) {
      for (r in seq_len(ns)) {
        Y[r, ] <- stats::isoreg(Y[r, ])$yf
      }
    }
    P <- X + P - Y

    Z <- Y + Q
    if (ns > 1L) {
      for (c_i in seq_len(ni)) {
        Z[, c_i] <- stats::isoreg(Z[, c_i])$yf
      }
    }
    Q <- Y + Q - Z

    if (max(abs(Z - X)) < 1e-7) {
      X <- Z
      break
    }
    X <- Z
  }

  sgn * X[cbind(r_idx, c_idx)]
}

min_positive_step <- function(u) {
  vals <- sort(unique(u[is.finite(u)]))
  if (length(vals) <= 1L) return(1.0)
  diffs <- diff(vals)
  pos <- diffs[diffs > 1e-12]
  if (!length(pos)) 1.0 else min(pos)
}

#' Compute Sample-Size Effort Score for `prefer` Ranking
#' @noRd
compute_preference_score <- function(
  grid_df,
  prefer = c("balanced", "fewer_subj", "fewer_items", "fewer_trials")
) {
  prefer <- match.arg(prefer)
  u_s <- normalize_unit(grid_df$n_subj)
  u_i <- normalize_unit(grid_df$n_item)
  trials_col <- if ("n_trials_scheduled_max" %in% names(grid_df)) {
    "n_trials_scheduled_max"
  } else {
    "n_trials_scheduled"
  }
  u_t <- normalize_unit(grid_df[[trials_col]])

  eps_s <- min(0.02, 0.25 * min_positive_step(u_s))
  eps_i <- min(0.02, 0.25 * min_positive_step(u_i))
  eps_t <- min(0.02, 0.25 * min_positive_step(u_t))

  switch(
    prefer,
    balanced = sqrt(0.5 * (u_s^2 + u_i^2)) + 0.01 * u_t,
    fewer_subj = u_s + eps_s * u_i,
    fewer_items = u_i + eps_i * u_s,
    fewer_trials = u_t + eps_t * sqrt(0.5 * (u_s^2 + u_i^2))
  )
}

preference_order <- function(
  df,
  prefer = NULL,
  tie1 = NULL,
  tie2 = NULL
) {
  n <- nrow(df)
  if (n <= 1L) return(seq_len(n))
  t1 <- if (is.null(tie1)) rep(0.0, n) else ifelse(is.na(tie1), Inf, tie1)
  t2 <- if (is.null(tie2)) rep(0.0, n) else ifelse(is.na(tie2), Inf, tie2)
  pref_sc <- if ("preference_score" %in% names(df)) {
    ifelse(is.na(df$preference_score), Inf, df$preference_score)
  } else {
    rep(0.0, n)
  }
  tr_max <- if ("n_trials_scheduled_max" %in% names(df)) {
    df$n_trials_scheduled_max
  } else if ("n_trials_scheduled" %in% names(df)) {
    df$n_trials_scheduled
  } else {
    as.numeric(df$n_subj) * as.numeric(df$n_item)
  }
  tr_mean <- if ("n_trials_scheduled_mean" %in% names(df)) df$n_trials_scheduled_mean else tr_max
  if (is.null(prefer)) {
    return(order(pref_sc, t1, t2, df$n_subj, df$n_item))
  }
  prefer <- match.arg(prefer, c("balanced", "fewer_subj", "fewer_items", "fewer_trials"))
  switch(
    prefer,
    fewer_subj = order(df$n_subj, df$n_item, tr_max, t1, t2),
    fewer_items = order(df$n_item, df$n_subj, tr_max, t1, t2),
    fewer_trials = order(tr_max, tr_mean, pref_sc, df$n_subj, df$n_item, t1, t2),
    balanced = order(pref_sc, tr_max, df$n_subj, df$n_item, t1, t2)
  )
}

#' Compute Contour-Straddle + Preference-Aware Acquisition Scores Over Grid
#' @noRd
compute_surrogate_acquisition <- function(
  grid_df,
  spec,
  gp_power,
  gp_aipe,
  gp_rope,
  evaluated_mask,
  prefer = "balanced",
  gp_joint = NULL,
  gp_valid = NULL
) {
  if (spec$criterion %in% c("joint", "power")) {
    if (spec$criterion == "joint") gp_power <- gp_joint
    thr <- stats::qnorm(spec$power_min)
    p_feas <- stats::pnorm((gp_power$mean_probit - thr) / gp_power$sd_probit)
    unc <- gp_power$sd_probit
  } else if (spec$criterion == "ci_width") {
    thr <- stats::qnorm(spec$ci_assurance)
    p_feas <- stats::pnorm((gp_aipe$mean_probit - thr) / gp_aipe$sd_probit)
    unc <- gp_aipe$sd_probit
  } else if (spec$criterion == "both") {
    thr_p <- stats::qnorm(spec$power_min)
    thr_w <- stats::qnorm(spec$ci_assurance)
    pf_p <- stats::pnorm((gp_power$mean_probit - thr_p) / gp_power$sd_probit)
    pf_w <- stats::pnorm((gp_aipe$mean_probit - thr_w) / gp_aipe$sd_probit)
    p_feas <- pf_p * pf_w
    unc <- 0.5 * (gp_power$sd_probit + gp_aipe$sd_probit)
  } else {
    thr <- stats::qnorm(spec$power_min)
    p_feas <- stats::pnorm((gp_rope$mean_probit - thr) / gp_rope$sd_probit)
    unc <- gp_rope$sd_probit
  }

  if (identical(spec$failure_policy, "conditional") &&
      !is.null(gp_valid) && any(is.finite(gp_valid$mean_prob))) {
    n_ref_v <- if (!is.null(gp_valid$n_ref)) gp_valid$n_ref else NA_real_
    thr_v <- stats::qnorm(smooth_prob_threshold(spec$min_valid_rate, n_ref_v))
    pf_v <- stats::pnorm((gp_valid$mean_probit - thr_v) / gp_valid$sd_probit)
    p_feas <- p_feas * pf_v
    unc <- 0.5 * (unc + gp_valid$sd_probit)
  }

  eff_score <- compute_preference_score(grid_df, prefer = prefer)
  plausible <- which(p_feas >= 0.35)
  e_ref <- if (length(plausible) > 0L) min(eff_score[plausible]) else max(eff_score)

  # Contour straddle peaks right along the target contour
  contour_straddle <- 4.0 * p_feas * (1.0 - p_feas) * (0.25 + unc)

  # Prioritize the lower-effort Pareto frontier
  excess_eff <- pmax(eff_score - e_ref, 0.0)
  pref_weight <- ifelse(eff_score <= e_ref, 1.20, 1.0) * exp(-1.8 * excess_eff)

  acq <- contour_straddle * pref_weight + 0.03 * unc
  acq[evaluated_mask] <- -Inf

  list(acq = acq, prob_feasible = p_feas, preference_score = eff_score)
}
