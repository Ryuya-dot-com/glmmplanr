# External Estimation Backends (lme4, glmmTMB, CmdStan/RStan)

fit_via_lme4 <- function(df, X, Z_s, Z_i, Z_t, chance, corr_mode, nagq, eff_idx) {
  fit_env <- build_numeric_glmm_df(df, X, Z_s, Z_i, Z_t, corr_mode, pkg = "lme4")
  fam <- mafc_binomial(chance)
  ctrl <- lme4::glmerControl(
    optimizer = "bobyqa",
    calc.derivs = (nagq > 0L),
    optCtrl = list(maxfun = if (nagq == 0L) 2500L else 8000L)
  )

  warn_msgs <- character(0)
  err_msg <- ""
  fit <- tryCatch(
    withCallingHandlers(
      lme4::glmer(
        formula = fit_env$formula,
        data = fit_env$data,
        family = fam,
        nAGQ = nagq,
        control = ctrl
      ),
      warning = function(w) {
        warn_msgs <<- c(warn_msgs, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      err_msg <<- conditionMessage(e)
      NULL
    }
  )

  if (is.null(fit)) {
    return(make_failed_fit("engine_error", err_msg))
  }

  b <- lme4::fixef(fit)
  V <- as.matrix(stats::vcov(fit))
  col_nm <- fit_env$x_names[eff_idx]
  est <- unname(b[col_nm])
  se <- sqrt(V[col_nm, col_nm])
  sing <- lme4::isSingular(fit, tol = 1e-4)

  bad_warn <- warn_msgs[grepl(
    "failed to converge|degenerate|unable to evaluate|maximum number of function evaluations",
    warn_msgs,
    ignore.case = TRUE
  )]
  nonconv <- length(bad_warn) > 0L
  conv <- !nonconv && is.finite(est) && is.finite(se) && (se > 0)

  vc <- lme4::VarCorr(fit)
  Sig_s <- extract_vc_matrix(vc, "subj", ncol(Z_s))
  Sig_i <- extract_vc_matrix(vc, "item", ncol(Z_i))
  Sig_t <- if (ncol(Z_t) > 0L) extract_vc_matrix(vc, "testlet", ncol(Z_t)) else NULL

  b_ord <- as.numeric(b[fit_env$x_names])
  V_ord <- V[fit_env$x_names, fit_env$x_names, drop = FALSE]

  list(
    est = est,
    se = se,
    beta_full = b_ord,
    vcov_full = V_ord,
    Sigma_subj = Sig_s,
    Sigma_item = Sig_i,
    Sigma_testlet = Sig_t,
    singular = isTRUE(sing),
    converged = conv,
    status = if (conv) "converged" else "optimizer_warning",
    error_reason = if (conv) "" else paste(unique(bad_warn), collapse = "; ")
  )
}

extract_vc_matrix <- function(vc_list, grp_name, q) {
  # Handles both correlated `(1 + ... | subj)` single entry and uncorrelated
  # `(1 | subj) + (0 + zs_1 | subj)` split entries (`subj`, `subj.1`, ...)
  nms <- names(vc_list)
  hits <- which(nms == grp_name | startsWith(nms, paste0(grp_name, ".")))
  if (length(hits) == 0L) {
    return(diag(0.0, nrow = q))
  }
  if (length(hits) == 1L && nrow(vc_list[[hits[1L]]]) == q) {
    return(as.matrix(vc_list[[hits[1L]]]))
  }
  vars <- unlist(lapply(hits, function(h) diag(as.matrix(vc_list[[h]]))))
  vars <- vars[seq_len(min(length(vars), q))]
  if (length(vars) < q) vars <- c(vars, rep(0.0, q - length(vars)))
  diag(as.numeric(vars), nrow = q)
}

fit_via_glmmtmb <- function(df, X, Z_s, Z_i, Z_t, chance, corr_mode, eff_idx) {
  fit_env <- build_numeric_glmm_df(df, X, Z_s, Z_i, Z_t, corr_mode, pkg = "glmmTMB")
  fam <- mafc_binomial(chance)

  warn_msgs <- character(0)
  err_msg <- ""
  fit <- tryCatch(
    withCallingHandlers(
      glmmTMB::glmmTMB(
        formula = fit_env$formula,
        data = fit_env$data,
        family = fam
      ),
      warning = function(w) {
        warn_msgs <<- c(warn_msgs, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      err_msg <<- conditionMessage(e)
      NULL
    }
  )
  if (is.null(fit)) {
    return(make_failed_fit("engine_error", err_msg))
  }

  b <- glmmTMB::fixef(fit)$cond
  V <- as.matrix(stats::vcov(fit)$cond)
  col_nm <- fit_env$x_names[eff_idx]
  est <- unname(b[col_nm])
  se <- sqrt(V[col_nm, col_nm])
  conv_code <- fit$fit$convergence
  pd_hess <- isTRUE(fit$sdr$pdHess)
  conv_ok <- identical(conv_code, 0L) && pd_hess && is.finite(est) && is.finite(se) && (se > 0)

  vc <- glmmTMB::VarCorr(fit)$cond
  sds <- unlist(lapply(vc, function(m) attr(m, "stddev")))
  cors <- unlist(lapply(vc, function(m) {
    cm <- attr(m, "correlation")
    if (is.null(cm) || nrow(cm) < 2L) numeric(0) else cm[lower.tri(cm)]
  }))
  sing <- any(sds < 1e-3, na.rm = TRUE) || any(abs(cors) > 0.995, na.rm = TRUE)

  Sig_s <- extract_vc_matrix(vc, "subj", ncol(Z_s))
  Sig_i <- extract_vc_matrix(vc, "item", ncol(Z_i))
  Sig_t <- if (ncol(Z_t) > 0L) extract_vc_matrix(vc, "testlet", ncol(Z_t)) else NULL

  st <- if (conv_ok) {
    "converged"
  } else if (!pd_hess) {
    "non_pd_hessian"
  } else {
    "optimizer_warning"
  }
  reason_txt <- if (conv_ok) {
    ""
  } else {
    msg_parts <- c(
      if (!identical(conv_code, 0L)) paste0("nlminb convergence=", conv_code, ": ", fit$fit$message),
      if (!pd_hess) "non-positive-definite Hessian",
      warn_msgs
    )
    paste(unique(msg_parts[nzchar(msg_parts)]), collapse = "; ")
  }

  list(
    est = est,
    se = se,
    beta_full = as.numeric(b[fit_env$x_names]),
    vcov_full = V[fit_env$x_names, fit_env$x_names, drop = FALSE],
    Sigma_subj = Sig_s,
    Sigma_item = Sig_i,
    Sigma_testlet = Sig_t,
    singular = isTRUE(sing),
    converged = conv_ok,
    status = st,
    error_reason = reason_txt
  )
}

build_stan_data <- function(df, X, Z_s, Z_i, Z_t, chance, corr_mode, spec) {
  s_fac <- as.integer(factor(df$subj_id))
  i_fac <- as.integer(factor(df$item_id))
  has_t <- ncol(Z_t) > 0L
  t_fac <- if (has_t) as.integer(factor(df$testlet_id)) else integer(0)

  des <- spec$dgp$design
  cell_z <- function(z, q) {
    if (q == 1L) matrix(1, des$n_cells, 1L) else z
  }
  quad <- gh_quad_1d(spec$dgp$n_quad)
  list(
    M = des$n_cells, X_cell = des$X_cells,
    Zs_cell = cell_z(des$Z_subj_cells, ncol(Z_s)),
    Zi_cell = cell_z(des$Z_item_cells, ncol(Z_i)),
    Zt_cell = if (has_t) cell_z(des$Z_testlet_cells, ncol(Z_t)) else matrix(0, des$n_cells, 0L),
    contrast = unname(spec$contrast_weights), K = length(quad$nodes),
    nodes = quad$nodes, weights = quad$weights,
    marginal = as.integer(spec$prob_estimand == "marginal"),
    N = nrow(df),
    y = as.integer(df$y),
    P = ncol(X),
    X = X,
    N_subj = max(s_fac),
    Q_subj = ncol(Z_s),
    subj_id = s_fac,
    Z_subj = Z_s,
    N_item = max(i_fac),
    Q_item = ncol(Z_i),
    item_id = i_fac,
    Z_item = Z_i,
    N_testlet = if (has_t) max(t_fac) else 0L,
    Q_testlet = if (has_t) ncol(Z_t) else 0L,
    testlet_id = t_fac,
    Z_testlet = if (has_t) Z_t else matrix(0.0, nrow = 0L, ncol = 0L),
    chance = chance,
    corr_mode = corr_mode,
    lkj_eta = spec$lkj_eta,
    prior_beta_sd = spec$prior_beta_sd,
    prior_tau_sd = spec$prior_tau_sd
  )
}

fit_via_stan_laplace <- function(
  df, X, Z_s, Z_i, Z_t, chance, corr_mode, spec, eff_idx
) {
  # Joint MAP/Laplace over hierarchical (beta, tau, z) suffers from Neal's funnel
  # degeneracy in crossed GLMMs; integrated Laplace is provided natively by
  # `bayes_laplace`, while Stan backends (`stan_hmc` and legacy alias `stan_laplace`)
  # execute full HMC/NUTS with strict convergence gates across CmdStan and RStan.
  fit_via_stan_hmc(df, X, Z_s, Z_i, Z_t, chance, corr_mode, spec, eff_idx)
}

evaluate_stan_hmc_diagnostics <- function(
  diag_draws,
  n_div,
  n_max_td,
  rhat_max = 1.01,
  aux_rhat_max = 1.02,
  aux_vars = character(0),
  ess_min = 100,
  has_posterior = requireNamespace("posterior", quietly = TRUE),
  has_rstan = requireNamespace("rstan", quietly = TRUE)
) {
  if (!isTRUE(has_posterior) && !isTRUE(has_rstan)) {
    return(list(
      ok = FALSE,
      status = "engine_error",
      reason = "Package 'posterior' (or 'rstan') is required to compute Stan MCMC convergence diagnostics."
    ))
  }
  if (!is.finite(n_div) || !is.finite(n_max_td)) {
    return(list(
      ok = FALSE,
      status = "mcmc_not_converged",
      reason = "HMC sampler diagnostics (divergences or treedepth) unavailable."
    ))
  }

  rhat_fn <- if (isTRUE(has_posterior)) posterior::rhat else rstan::Rhat
  ess_bulk_fn <- if (isTRUE(has_posterior)) posterior::ess_bulk else rstan::ess_bulk
  ess_tail_fn <- if (isTRUE(has_posterior)) posterior::ess_tail else rstan::ess_tail

  worst_rhat <- -Inf
  min_ess_b <- Inf
  min_ess_t <- Inf
  worst_par <- ""
  rhat_failed <- FALSE
  rhat_fail_par <- ""
  rhat_fail_val <- NA_real_

  for (nm in names(diag_draws)) {
    arr <- diag_draws[[nm]]
    if (is.null(arr) || length(arr) == 0L) {
      return(list(
        ok = FALSE,
        status = "mcmc_not_converged",
        reason = sprintf("Missing HMC draws for monitored parameter '%s'.", nm)
      ))
    }
    rh <- tryCatch(as.numeric(rhat_fn(arr))[1L], error = function(e) NA_real_)
    eb <- tryCatch(as.numeric(ess_bulk_fn(arr))[1L], error = function(e) NA_real_)
    et <- tryCatch(as.numeric(ess_tail_fn(arr))[1L], error = function(e) NA_real_)
    if (!is.finite(rh) || !is.finite(eb) || !is.finite(et)) {
      return(list(
        ok = FALSE,
        status = "mcmc_not_converged",
        reason = sprintf("Non-finite HMC diagnostics for '%s' (Rhat=%s, ESS_bulk=%s, ESS_tail=%s).", nm, rh, eb, et)
      ))
    }
    lim <- if (nm %in% aux_vars) aux_rhat_max else rhat_max
    if (rh > lim && (!rhat_failed || rh > rhat_fail_val)) {
      rhat_failed <- TRUE
      rhat_fail_par <- nm
      rhat_fail_val <- rh
    }
    if (rh > worst_rhat) {
      worst_rhat <- rh
      worst_par <- nm
    }
    if (eb < min_ess_b) min_ess_b <- eb
    if (et < min_ess_t) min_ess_t <- et
  }

  ok <- (!rhat_failed) && (min_ess_b >= ess_min) && (min_ess_t >= ess_min) &&
    (n_div == 0L) && (n_max_td == 0L)
  rep_par <- if (rhat_failed) rhat_fail_par else worst_par
  rep_rhat <- if (rhat_failed) rhat_fail_val else worst_rhat
  reason <- if (ok) {
    ""
  } else {
    sprintf(
      "HMC diagnostics failed (max Rhat=%.3f [%s], min ESS_bulk=%.0f, min ESS_tail=%.0f, div=%d, max_td=%d)",
      rep_rhat, rep_par, min_ess_b, min_ess_t, as.integer(n_div), as.integer(n_max_td)
    )
  }
  list(
    ok = ok,
    status = if (ok) "converged" else "mcmc_not_converged",
    reason = reason,
    max_rhat = worst_rhat,
    min_ess_bulk = min_ess_b,
    min_ess_tail = min_ess_t,
    n_div = as.integer(n_div),
    n_max_td = as.integer(n_max_td)
  )
}

fit_via_stan_hmc <- function(
  df, X, Z_s, Z_i, Z_t, chance, corr_mode, spec, eff_idx
) {
  stan_data <- build_stan_data(df, X, Z_s, Z_i, Z_t, chance, corr_mode, spec)
  mod_info <- get_compiled_stan_model(prefer = "any")
  if (is.null(mod_info)) return(make_failed_fit("engine_error", "Stan model unavailable"))

  focal_nm <- sprintf("beta[%d]", eff_idx)
  b_vars <- sprintf("beta[%d]", seq_len(ncol(X)))
  tau_vars <- c(
    sprintf("tau_subj[%d]", seq_len(ncol(Z_s))),
    sprintf("tau_item[%d]", seq_len(ncol(Z_i))),
    if (ncol(Z_t) > 0L) sprintf("tau_testlet[%d]", seq_len(ncol(Z_t)))
  )
  monitored_vars <- unique(c(
    focal_nm,
    if (spec$scale %in% c("pp", "prob")) "prob_contrast",
    tau_vars
  ))

  if (identical(mod_info$backend, "cmdstan") &&
      !requireNamespace("posterior", quietly = TRUE) &&
      !requireNamespace("rstan", quietly = TRUE)) {
    return(make_failed_fit(
      "engine_error",
      "Package 'posterior' is required for CmdStan MCMC convergence diagnostics."
    ))
  }

  run_attempt <- function(chains, iter_warmup, iter_sampling, adapt_delta, max_td) {
    if (identical(mod_info$backend, "cmdstan")) {
      err_msg <- ""
      fit <- tryCatch(
        mod_info$model$sample(
          data = stan_data, chains = chains, parallel_chains = 1L,
          iter_warmup = iter_warmup, iter_sampling = iter_sampling,
          adapt_delta = adapt_delta, max_treedepth = max_td,
          refresh = 0L, show_messages = FALSE, show_exceptions = FALSE
        ),
        error = function(e) {
          err_msg <<- conditionMessage(e)
          NULL
        }
      )
      if (is.null(fit)) return(make_failed_fit("engine_error", err_msg))
      rc <- tryCatch(fit$return_codes(), error = function(e) 0L)
      if (any(rc != 0L)) {
        return(make_failed_fit("mcmc_not_converged", sprintf("CmdStan HMC return_code=%s", paste(rc, collapse = ","))))
      }
      draws_b <- fit$draws(variables = b_vars, format = "matrix")
      contrast_draws <- as.numeric(fit$draws(variables = "prob_contrast", format = "matrix"))
      mon_arr <- fit$draws(variables = monitored_vars, format = "draws_array")
      diag_draws <- stats::setNames(
        lapply(monitored_vars, function(v) mon_arr[, , v, drop = TRUE]),
        monitored_vars
      )
      diag_df <- tryCatch(fit$diagnostic_summary(quiet = TRUE), error = function(e) NULL)
      n_div <- if (!is.null(diag_df) && !is.null(diag_df$num_divergent)) sum(diag_df$num_divergent) else NA_integer_
      n_max_td <- if (!is.null(diag_df) && !is.null(diag_df$num_max_treedepth)) sum(diag_df$num_max_treedepth) else NA_integer_
    } else {
      err_msg <- ""
      sfit <- tryCatch(
        suppressWarnings(rstan::sampling(
          mod_info$model,
          data = stan_data,
          chains = chains,
          iter = iter_warmup + iter_sampling,
          warmup = iter_warmup,
          cores = 1L,
          refresh = 0L,
          show_messages = FALSE,
          control = list(adapt_delta = adapt_delta, max_treedepth = max_td)
        )),
        error = function(e) {
          err_msg <<- conditionMessage(e)
          NULL
        }
      )
      if (is.null(sfit)) return(make_failed_fit("engine_error", err_msg))
      ext <- rstan::extract(sfit, pars = c("beta", "prob_contrast"))
      draws_b <- ext$beta
      contrast_draws <- as.numeric(ext$prob_contrast)
      arr_all <- rstan::extract(sfit, permuted = FALSE, inc_warmup = FALSE, pars = monitored_vars)
      diag_draws <- stats::setNames(
        lapply(monitored_vars, function(v) arr_all[, , v]),
        monitored_vars
      )
      sam <- do.call(rbind, rstan::get_sampler_params(sfit, inc_warmup = FALSE))
      n_div <- sum(sam[, "divergent__"])
      n_max_td <- sum(sam[, "treedepth__"] >= max_td)
    }

    b_est <- colMeans(draws_b)
    V_est <- stats::cov(draws_b)
    est <- unname(b_est[eff_idx])
    se <- sqrt(V_est[eff_idx, eff_idx])
    prob_est <- mean(contrast_draws)
    prob_se <- stats::sd(contrast_draws)

    diag_res <- evaluate_stan_hmc_diagnostics(
      diag_draws, n_div = n_div, n_max_td = n_max_td, aux_vars = tau_vars
    )
    num_ok <- is.finite(est) && is.finite(se) && (se > 0) &&
      is.finite(prob_est) && is.finite(prob_se) && (prob_se > 0)
    conv <- isTRUE(diag_res$ok) && num_ok

    reason_txt <- if (conv) {
      ""
    } else if (!isTRUE(diag_res$ok)) {
      diag_res$reason
    } else {
      "Non-finite HMC posterior summary"
    }

    list(
      est = est,
      se = se,
      beta_full = as.numeric(b_est),
      vcov_full = as.matrix(V_est),
      prob_est = prob_est,
      prob_se = prob_se,
      singular = NA,
      converged = conv,
      status = if (conv) "converged" else diag_res$status,
      error_reason = reason_txt
    )
  }

  # Initial fast HMC run; if convergence diagnostics fail, retry once with 4 chains
  # and tighter NUTS adaptation matching the validation protocol.
  res <- run_attempt(chains = 2L, iter_warmup = 400L, iter_sampling = 400L, adapt_delta = 0.90, max_td = 10L)
  if (!isTRUE(res$converged) && identical(res$status, "mcmc_not_converged")) {
    res <- run_attempt(chains = 4L, iter_warmup = 600L, iter_sampling = 800L, adapt_delta = 0.95, max_td = 12L)
  }
  res
}

get_compiled_stan_model <- function(prefer = c("any", "cmdstan", "rstan")) {
  prefer <- match.arg(prefer)
  stan_file <- system.file("stan", "crossed_binomial_mafc.stan", package = "glmmplanr")
  if (!nzchar(stan_file) || !file.exists(stan_file)) {
    stan_file <- file.path("inst", "stan", "crossed_binomial_mafc.stan")
  }

  if (prefer %in% c("any", "cmdstan")) {
    if (!is.null(.glmmplanr_env$stan_model_cmdstan)) {
      return(.glmmplanr_env$stan_model_cmdstan)
    }
    mod <- NULL
    if (file.exists(stan_file) && requireNamespace("cmdstanr", quietly = TRUE) &&
        !is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL))) {
      mod <- tryCatch(cmdstanr::cmdstan_model(stan_file, dir = tempdir()), error = function(e) NULL)
    } else if (requireNamespace("instantiate", quietly = TRUE) &&
               isTRUE(tryCatch(instantiate::stan_cmdstan_exists(), error = function(e) FALSE))) {
      mod <- tryCatch(instantiate::stan_package_model(
        name = "crossed_binomial_mafc", package = "glmmplanr", dir = tempdir()
      ), error = function(e) NULL)
    }
    if (!is.null(mod)) {
      res <- list(backend = "cmdstan", model = mod)
      .glmmplanr_env$stan_model_cmdstan <- res
      return(res)
    }
    if (prefer == "cmdstan") return(NULL)
  }

  if (prefer %in% c("any", "rstan")) {
    if (!is.null(.glmmplanr_env$stan_model_rstan)) {
      return(.glmmplanr_env$stan_model_rstan)
    }
    if (file.exists(stan_file) && requireNamespace("rstan", quietly = TRUE)) {
      mod <- tryCatch(rstan::stan_model(stan_file), error = function(e) NULL)
      if (!is.null(mod)) {
        res <- list(backend = "rstan", model = mod)
        .glmmplanr_env$stan_model_rstan <- res
        return(res)
      }
    }
  }
  NULL
}

build_numeric_glmm_df <- function(df, X, Z_s, Z_i, Z_t, corr_mode, pkg = "lme4") {
  out <- data.frame(
    y = as.integer(df$y),
    subj = factor(df$subj_id),
    item = factor(df$item_id)
  )
  if (ncol(Z_t) > 0L) {
    out$testlet <- factor(df$testlet_id)
  }
  x_names <- character(ncol(X))
  x_names[1L] <- "(Intercept)"
  if (ncol(X) > 1L) {
    for (j in 2L:ncol(X)) {
      nm <- paste0("fx_", j - 1L)
      out[[nm]] <- X[, j]
      x_names[j] <- nm
    }
  }
  fx_rhs <- if (ncol(X) > 1L) paste(x_names[-1L], collapse = " + ") else "1"

  re_block <- function(Z_mat, grp_nm, prefix) {
    if (ncol(Z_mat) == 1L) {
      return(sprintf("(1 | %s)", grp_nm))
    }
    z_nms <- character(ncol(Z_mat) - 1L)
    for (k in 2L:ncol(Z_mat)) {
      zn <- sprintf("%s_%d", prefix, k - 1L)
      out[[zn]] <<- Z_mat[, k]
      z_nms[k - 1L] <- zn
    }
    if (corr_mode == 1L) {
      sprintf("(1 + %s | %s)", paste(z_nms, collapse = " + "), grp_nm)
    } else if (pkg == "glmmTMB") {
      sprintf("diag(1 + %s | %s)", paste(z_nms, collapse = " + "), grp_nm)
    } else {
      parts <- c(sprintf("(1 | %s)", grp_nm), sprintf("(0 + %s | %s)", z_nms, grp_nm))
      paste(parts, collapse = " + ")
    }
  }

  re_s <- re_block(Z_s, "subj", "zs")
  re_i <- re_block(Z_i, "item", "zi")
  re_terms <- c(re_s, re_i)
  if (ncol(Z_t) > 0L) {
    re_terms <- c(re_terms, re_block(Z_t, "testlet", "zt"))
  }

  fml <- stats::as.formula(sprintf("y ~ %s + %s", fx_rhs, paste(re_terms, collapse = " + ")))
  list(formula = fml, data = out, x_names = x_names)
}
