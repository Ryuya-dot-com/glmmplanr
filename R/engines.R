compact_idx0 <- function(ids) {
  mx <- max(ids)
  present <- logical(mx)
  present[ids] <- TRUE
  map <- cumsum(present)
  list(idx0 = map[ids] - 1L, n = map[mx])
}

#' Fit a Single Simulated Dataset with the Configured GLMM Engine
#'
#' Returns point estimates and standard errors on both the conditional logit
#' scale (`est`, `se`) and the response-probability difference scale
#' (`prob_est`, `prob_se`), along with strict convergence status and captured
#' backend error/warning messages.
#'
#' @noRd
fit_single_dataset <- function(
  df,
  spec,
  engine_use = spec$engine,
  re_structure = spec$analysis_re,
  fast_stage = FALSE
) {
  des <- spec$dgp$design
  chance_fit <- if (spec$analysis_link == "match_dgp") spec$dgp$chance else 0.0

  cid <- df$cell_id
  X <- des$X_cells[cid, , drop = FALSE]

  if (!nrow(X) || (length(unique(cid)) < ncol(X))) {
    return(make_failed_fit("insufficient_data", "Observed trials omit required design cells."))
  }
  if (re_structure == "intercept_only") {
    Z_s <- matrix(1.0, nrow = nrow(df), ncol = 1L, dimnames = list(NULL, "(Intercept)"))
    Z_i <- matrix(1.0, nrow = nrow(df), ncol = 1L, dimnames = list(NULL, "(Intercept)"))
    Z_t <- if (des$has_testlet) {
      matrix(1.0, nrow = nrow(df), ncol = 1L, dimnames = list(NULL, "(Intercept)"))
    } else {
      matrix(0.0, nrow = nrow(df), ncol = 0L)
    }
    corr_mode <- 0L
  } else {
    Z_s <- des$Z_subj_cells[cid, , drop = FALSE]
    Z_i <- des$Z_item_cells[cid, , drop = FALSE]
    Z_t <- if (des$has_testlet) {
      des$Z_testlet_cells[cid, , drop = FALSE]
    } else {
      matrix(0.0, nrow = nrow(df), ncol = 0L)
    }
    corr_mode <- if (re_structure == "maximal") 1L else 0L
  }

  eff_idx <- match(spec$effect_col, colnames(X))

  warm_cpp <- NULL
  fit_response <- function(y) {
    df$y <- y
    if (engine_use %in% c("cpp_laplace", "bayes_laplace")) {
      res <- fit_via_cpp(
        df = df, X = X, Z_s = Z_s, Z_i = Z_i, Z_t = Z_t,
        chance = chance_fit, corr_mode = corr_mode,
        regularized = (engine_use == "bayes_laplace"),
        spec = spec, fast_stage = fast_stage, eff_idx = eff_idx,
        warm = warm_cpp
      )
      if (is.null(warm_cpp) && isTRUE(res$converged)) warm_cpp <<- res
      res
    } else if (engine_use %in% c("lme4_fast", "lme4")) {
      fit_via_lme4(
        df = df, X = X, Z_s = Z_s, Z_i = Z_i, Z_t = Z_t,
        chance = chance_fit, corr_mode = corr_mode,
        nagq = if (engine_use == "lme4_fast") 0L else 1L,
        eff_idx = eff_idx
      )
    } else if (engine_use == "glmmTMB") {
      fit_via_glmmtmb(
        df = df, X = X, Z_s = Z_s, Z_i = Z_i, Z_t = Z_t,
        chance = chance_fit, corr_mode = corr_mode, eff_idx = eff_idx
      )
    } else if (engine_use == "stan_laplace") {
      fit_via_stan_laplace(
        df = df, X = X, Z_s = Z_s, Z_i = Z_i, Z_t = Z_t,
        chance = chance_fit, corr_mode = corr_mode,
        spec = spec, eff_idx = eff_idx
      )
    } else if (engine_use == "stan_hmc") {
      fit_via_stan_hmc(
        df = df, X = X, Z_s = Z_s, Z_i = Z_i, Z_t = Z_t,
        chance = chance_fit, corr_mode = corr_mode,
        spec = spec, eff_idx = eff_idx
      )
    } else {
      cli::cli_abort("Unknown engine {.val {engine_use}}.")
    }
  }

  cpp_boot_ctx <- if (engine_use %in% c("cpp_laplace", "bayes_laplace")) {
    cs <- compact_idx0(df$subj_id)
    ci <- compact_idx0(df$item_id)
    ct <- if (ncol(Z_t) > 0L) compact_idx0(df$testlet_id) else list(idx0 = integer(0), n = 0L)
    list(
      X = X, s_idx0 = cs$idx0, Z_s = Z_s, n_s = cs$n,
      i_idx0 = ci$idx0, Z_i = Z_i, n_i = ci$n,
      t_idx0 = ct$idx0, Z_t = Z_t, n_t = ct$n,
      corr_mode = corr_mode, regularized = (engine_use == "bayes_laplace"),
      eff_idx = eff_idx
    )
  } else {
    NULL
  }

  summarize_probability_fit(
    df$y, fit_response, X, Z_s, Z_i, Z_t, df$subj_id, df$item_id, df$testlet_id,
    spec, re_structure, chance_fit, fast_stage, cpp_boot_ctx = cpp_boot_ctx,
    unpenalized_cells = if (engine_use %in% c("cpp_laplace", "lme4_fast", "lme4", "glmmTMB")) cid
  )
}

# Propagate joint uncertainty in beta and random-effect covariance matrices into
# a marginal contrast via either analytical joint delta or pre-allocated C++ / R bootstrap.
summarize_probability_fit <- function(
  y, fit_response, X, Z_s, Z_i, Z_t, subj, item, testlet, spec, re_structure, chance, fast_stage,
  cpp_boot_ctx = NULL, unpenalized_cells = NULL
) {
  refit <- function(response) {
    if (!is.null(unpenalized_cells)) {
      design <- spec$dgp$design
      cell <- cpp_constant_response_cell(response, unpenalized_cells, design$n_cells)
      if (cell > 0L) {
        return(make_failed_fit("separation", sprintf(
          "Separation in fixed-effect cell '%s': only one response category was observed. Its unpenalized cell logit has no finite estimate.",
          design$cells$cell_name[cell])))
      }
    }
    fit_response(response)
  }
  fit <- refit(y)
  if (!is.null(fit$prob_est)) return(fit)
  n_s_obs <- length(unique(subj))
  n_i_obs <- length(unique(item))
  n_t_obs <- if (ncol(Z_t) > 0L && length(testlet) > 0L) length(unique(testlet)) else 0L

  delta <- function(f, include_joint = TRUE) {
    v_ts <- if (include_joint) f$var_log_tau_subj else NULL
    v_ti <- if (include_joint) f$var_log_tau_item else NULL
    v_tt <- if (include_joint) f$var_log_tau_testlet else NULL
    J_bt <- if (include_joint) f$J_beta_theta else NULL
    if (include_joint && is.null(v_ts) && !is.null(f$Sigma_subj)) {
      v_ts <- rep(1.0 / max(1.0, 0.8 * n_s_obs), nrow(f$Sigma_subj))
      v_ti <- rep(1.0 / max(1.0, 0.8 * n_i_obs), nrow(f$Sigma_item))
      if (!is.null(f$Sigma_testlet) && n_t_obs > 0L) {
        v_tt <- rep(1.0 / max(1.0, 0.8 * n_t_obs), nrow(f$Sigma_testlet))
      }
    }
    estimate_prob_effect_delta(
      f$beta_full, f$vcov_full, f$Sigma_subj, f$Sigma_item, f$Sigma_testlet,
      spec$dgp$design, chance, spec$contrast_weights, spec$prob_estimand,
      re_structure, spec$dgp$n_quad,
      var_log_tau_subj = v_ts, var_log_tau_item = v_ti,
      var_log_tau_testlet = v_tt, J_beta_theta = J_bt
    )
  }

  fit$prob_est <- fit$prob_se <- NA_real_
  if (!isTRUE(fit$converged)) return(fit)
  tryCatch({
    point <- delta(fit, include_joint = TRUE)
    fit$prob_est <- point$est
    fit$prob_se <- point$se
    se_method <- if (is.null(spec$prob_se_method)) "bootstrap" else spec$prob_se_method

    if (spec$prob_estimand == "marginal" && spec$scale != "logit" &&
        !fast_stage && identical(se_method, "bootstrap")) {
      B <- if (is.null(spec$prob_boot)) 99L else as.integer(spec$prob_boot)
      if (!is.null(cpp_boot_ctx)) {
        des <- spec$dgp$design
        Zs_c <- if (re_structure == "intercept_only") matrix(1.0, des$n_cells, 1L) else des$Z_subj_cells
        Zi_c <- if (re_structure == "intercept_only") matrix(1.0, des$n_cells, 1L) else des$Z_item_cells
        Zt_c <- if (des$has_testlet) {
          if (re_structure == "intercept_only") matrix(1.0, des$n_cells, 1L) else des$Z_testlet_cells
        } else {
          matrix(0.0, des$n_cells, 0L)
        }
        q1 <- gh_quad_1d(spec$dgp$n_quad)
        max_em <- if (chance > 0) 400L else 160L
        max_inner <- if (chance > 0) 44L else 22L
        ws_b <- cpp_boot_ctx$ws_ptr
        boot <- if (!is.null(ws_b) && cpp_glmm_workspace_valid(ws_b)) {
          cpp_prob_bootstrap_batch_ws(
            ws_ptr = ws_b,
            n_boot = B,
            cell_id = as.integer(unpenalized_cells),
            beta_hat = fit$beta_full,
            Sigma_subj_hat = fit$Sigma_subj,
            Sigma_item_hat = fit$Sigma_item,
            Sigma_testlet_hat = fit$Sigma_testlet,
            X_cells = des$X_cells,
            Zs_cells = Zs_c,
            Zi_cells = Zi_c,
            Zt_cells = Zt_c,
            contrast_weights = unname(spec$contrast_weights),
            gh_nodes = q1$nodes,
            gh_weights = q1$weights,
            chance = chance,
            corr_mode = cpp_boot_ctx$corr_mode,
            regularized = cpp_boot_ctx$regularized,
            lkj_eta = spec$lkj_eta,
            prior_beta_sd = spec$prior_beta_sd,
            prior_tau_sd = spec$prior_tau_sd,
            max_em_iter = max_em,
            max_pirls_iter = max_inner,
            tol = 1.5e-4,
            eff_idx_1based = cpp_boot_ctx$eff_idx,
            marginal = TRUE
          )
        } else {
          cpp_prob_bootstrap_batch(
            n_boot = B,
            cell_id = as.integer(unpenalized_cells),
            X = cpp_boot_ctx$X,
            subj_idx = cpp_boot_ctx$s_idx0,
            Z_subj = cpp_boot_ctx$Z_s,
            n_subj = cpp_boot_ctx$n_s,
            item_idx = cpp_boot_ctx$i_idx0,
            Z_item = cpp_boot_ctx$Z_i,
            n_item = cpp_boot_ctx$n_i,
            testlet_idx = cpp_boot_ctx$t_idx0,
            Z_testlet = cpp_boot_ctx$Z_t,
            n_testlet = cpp_boot_ctx$n_t,
            beta_hat = fit$beta_full,
            Sigma_subj_hat = fit$Sigma_subj,
            Sigma_item_hat = fit$Sigma_item,
            Sigma_testlet_hat = fit$Sigma_testlet,
            X_cells = des$X_cells,
            Zs_cells = Zs_c,
            Zi_cells = Zi_c,
            Zt_cells = Zt_c,
            contrast_weights = unname(spec$contrast_weights),
            gh_nodes = q1$nodes,
            gh_weights = q1$weights,
            chance = chance,
            corr_mode = cpp_boot_ctx$corr_mode,
            regularized = cpp_boot_ctx$regularized,
            lkj_eta = spec$lkj_eta,
            prior_beta_sd = spec$prior_beta_sd,
            prior_tau_sd = spec$prior_tau_sd,
            max_em_iter = max_em,
            max_pirls_iter = max_inner,
            tol = 1.5e-4,
            eff_idx_1based = cpp_boot_ctx$eff_idx,
            marginal = TRUE
          )
        }
      } else {
        group <- function(Z, id, Sigma) {
          index <- match(id, unique(id))
          list(design = Z %*% safe_chol_lower(Sigma), index = index, n = max(index))
        }
        groups <- list(group(Z_s, subj, fit$Sigma_subj), group(Z_i, item, fit$Sigma_item))
        if (ncol(Z_t)) groups[[3L]] <- group(Z_t, testlet, fit$Sigma_testlet)
        fixed <- drop(X %*% fit$beta_full)
        boot <- replicate(B, {
          eta <- fixed
          for (g in groups) {
            u <- matrix(stats::rnorm(g$n * ncol(g$design)), nrow = g$n)
            eta <- eta + rowSums(g$design * u[g$index, , drop = FALSE])
          }
          y_star <- stats::rbinom(length(y), 1L, mafc_plogis(eta, chance))
          tryCatch({
            f <- refit(y_star)
            if (isTRUE(f$converged)) delta(f, include_joint = FALSE)$est else NA_real_
          }, error = function(e) NA_real_)
        })
      }
      ok <- is.finite(boot)
      if (sum(ok) < max(2L, ceiling(0.9 * B))) {
        stop(sprintf("Marginal bootstrap: %d/%d valid refits", sum(ok), B))
      }
      fit$prob_se <- stats::sd(boot[ok])
      if (!is.finite(fit$prob_se) || fit$prob_se <= 0) stop("Invalid marginal bootstrap SE")
    }
    fit
  }, error = function(e) {
    fit$prob_se <- NA_real_
    fit$error_reason <- conditionMessage(e)
    fit
  })
}

engine_supports_singularity <- function(engine) {
  !(engine %in% c("bayes_laplace", "stan_hmc", "stan_laplace"))
}

make_failed_fit <- function(status = "engine_error", reason = "", singular = FALSE) {
  list(
    est = NA_real_,
    se = NA_real_,
    prob_est = NA_real_,
    prob_se = NA_real_,
    beta_full = NULL,
    vcov_full = NULL,
    Sigma_subj = NULL,
    Sigma_item = NULL,
    Sigma_testlet = NULL,
    singular = if (is.na(singular)) NA else isTRUE(singular),
    converged = FALSE,
    status = status,
    error_reason = trimws(as.character(reason))
  )
}

fit_via_cpp <- function(
  df, X, Z_s, Z_i, Z_t, chance, corr_mode, regularized, spec, fast_stage, eff_idx,
  warm = NULL
) {
  cs <- compact_idx0(df$subj_id)
  ci <- compact_idx0(df$item_id)
  ct <- if (ncol(Z_t) > 0L) compact_idx0(df$testlet_id) else list(idx0 = integer(0), n = 0L)

  run_cpp_glmm_core(
    y = as.integer(df$y),
    X = X,
    s_idx0 = cs$idx0,
    Z_s = Z_s,
    n_s = cs$n,
    i_idx0 = ci$idx0,
    Z_i = Z_i,
    n_i = ci$n,
    t_idx0 = ct$idx0,
    Z_t = Z_t,
    n_t = ct$n,
    chance = chance,
    corr_mode = corr_mode,
    regularized = regularized,
    spec = spec,
    fast_stage = fast_stage,
    eff_idx = eff_idx,
    warm = warm
  )
}

run_cpp_glmm_core <- function(
  y, X, s_idx0, Z_s, n_s, i_idx0, Z_i, n_i, t_idx0, Z_t, n_t,
  chance, corr_mode, regularized, spec, fast_stage, eff_idx,
  warm = NULL, ws_ptr = NULL
) {
  # The corrected diagonal-covariance path uses the same numerical accuracy
  # when screening and confirming: an 85-iteration cap can discard otherwise
  # valid designs. Screening still uses fewer simulations and skips bootstrap.
  coarse_fit <- fast_stage && !(chance == 0 && corr_mode == 0L)
  max_em <- if (coarse_fit) 85L else if (chance > 0) 400L else 160L
  max_inner <- if (coarse_fit) 15L else if (chance > 0) 44L else 22L
  tol_em <- if (coarse_fit) 6.0e-4 else 1.5e-4

  err_msg <- ""
  fit <- tryCatch(
    if (!is.null(ws_ptr) && cpp_glmm_workspace_valid(ws_ptr)) {
      cpp_fit_crossed_glmm_ws(
        ws_ptr = ws_ptr,
        y = y,
        chance = chance,
        corr_mode = corr_mode,
        regularized = regularized,
        lkj_eta = spec$lkj_eta,
        prior_beta_sd = spec$prior_beta_sd,
        prior_tau_sd = spec$prior_tau_sd,
        max_em_iter = max_em,
        max_pirls_iter = max_inner,
        tol = tol_em,
        init_beta = if (!is.null(warm)) warm$beta_full else NULL,
        init_Sigma_subj = if (!is.null(warm)) warm$Sigma_subj else NULL,
        init_Sigma_item = if (!is.null(warm)) warm$Sigma_item else NULL,
        init_Sigma_testlet = if (!is.null(warm)) warm$Sigma_testlet else NULL
      )
    } else {
      cpp_fit_crossed_glmm(
        y = y,
        X = X,
        subj_idx = s_idx0,
        Z_subj = Z_s,
        n_subj = n_s,
        item_idx = i_idx0,
        Z_item = Z_i,
        n_item = n_i,
        testlet_idx = t_idx0,
        Z_testlet = Z_t,
        n_testlet = n_t,
        chance = chance,
        corr_mode = corr_mode,
        regularized = regularized,
        lkj_eta = spec$lkj_eta,
        prior_beta_sd = spec$prior_beta_sd,
        prior_tau_sd = spec$prior_tau_sd,
        max_em_iter = max_em,
        max_pirls_iter = max_inner,
        tol = tol_em,
        init_beta = if (!is.null(warm)) warm$beta_full else NULL,
        init_Sigma_subj = if (!is.null(warm)) warm$Sigma_subj else NULL,
        init_Sigma_item = if (!is.null(warm)) warm$Sigma_item else NULL,
        init_Sigma_testlet = if (!is.null(warm)) warm$Sigma_testlet else NULL
      )
    },
    error = function(e) {
      err_msg <<- conditionMessage(e)
      NULL
    }
  )

  if (is.null(fit)) {
    status <- if (startsWith(err_msg, "non_pd_hessian")) "non_pd_hessian" else "engine_error"
    return(make_failed_fit(status, err_msg, singular = if (regularized) NA else FALSE))
  }

  st <- as.character(fit$status)
  conv <- isTRUE(fit$converged) && identical(st, "converged") &&
    is.finite(fit$beta[eff_idx]) && is.finite(fit$se[eff_idx]) && (fit$se[eff_idx] > 0)

  reason_txt <- if (conv) {
    ""
  } else if (identical(st, "max_iter_finite")) {
    score <- fit$covariance_score_norm
    detail <- if (length(score) == 1L && is.finite(score) && score > 0) {
      sprintf(", scaled covariance score=%.2e", score)
    } else ""
    sprintf("max_em_iter (%d) reached before convergence (diff=%.2e%s)",
      max_em, fit$outer_diff, detail)
  } else {
    st
  }

  list(
    est = fit$beta[eff_idx],
    se = fit$se[eff_idx],
    beta_full = as.numeric(fit$beta),
    vcov_full = as.matrix(fit$vcov),
    Sigma_subj = as.matrix(fit$Sigma_subj),
    Sigma_item = as.matrix(fit$Sigma_item),
    Sigma_testlet = if (!is.null(fit$Sigma_testlet)) as.matrix(fit$Sigma_testlet) else NULL,
    var_log_tau_subj = as.numeric(fit$var_log_tau_subj),
    var_log_tau_item = as.numeric(fit$var_log_tau_item),
    var_log_tau_testlet = if (!is.null(fit$var_log_tau_testlet)) as.numeric(fit$var_log_tau_testlet) else NULL,
    J_beta_theta = as.matrix(fit$J_beta_theta),
    singular = if (regularized || is.na(fit$singular)) NA else isTRUE(fit$singular),
    converged = conv,
    status = st,
    error_reason = reason_txt
  )
}

#' Pre-Index Trial Matrices and Cluster Indices Once per (n_subj, n_item) Schedule
#' @noRd
prepare_cpp_schedule_cache <- function(spec, schedule) {
  dgp <- spec$dgp
  des <- dgp$design
  tr <- schedule$trials
  n_tr <- nrow(tr)
  ns <- schedule$n_subj
  ni <- schedule$n_item
  has_t <- isTRUE(des$has_testlet)
  nt <- if (has_t) schedule$n_testlet else 0L

  cid <- tr$cell_id
  sid <- tr$subj_id
  iid <- tr$item_id
  tid <- if (has_t) tr$testlet_id else integer(0)

  s_fac <- as.integer(factor(sid))
  i_fac <- as.integer(factor(iid))
  t_fac <- if (has_t) as.integer(factor(tid)) else integer(0)

  X_full <- des$X_cells[cid, , drop = FALSE]
  p_s <- ncol(des$Z_subj_cells)
  p_i <- ncol(des$Z_item_cells)

  list(
    n_tr = n_tr,
    ns = ns,
    ni = ni,
    has_t = has_t,
    nt = nt,
    p_s = p_s,
    p_i = p_i,
    p_t = if (has_t) ncol(des$Z_testlet_cells) else 0L,
    scalar_dgp = (p_s == 1L && p_i == 1L &&
                    all(des$Z_subj_cells == 1.0) && all(des$Z_item_cells == 1.0)),
    cell_id = cid,
    sid = sid,
    iid = iid,
    tid = tid,
    mu_trial = dgp$mu_cells[cid],
    zs_dgp = des$Z_subj_cells[cid, , drop = FALSE],
    zi_dgp = des$Z_item_cells[cid, , drop = FALSE],
    zt_dgp = if (has_t) des$Z_testlet_cells[cid, , drop = FALSE] else NULL,
    t_Ls = t(dgp$L_subj),
    t_Li = t(dgp$L_item),
    t_Lt = if (has_t) t(dgp$L_testlet) else NULL,
    chance_dgp = dgp$chance,
    chance_fit = if (spec$analysis_link == "match_dgp") dgp$chance else 0.0,
    design = des,
    participant_testlet_sd = testlet_interaction_sd(dgp),
    X_full = X_full,
    Z_int = matrix(1.0, nrow = n_tr, ncol = 1L),
    Z_empty = matrix(0.0, nrow = n_tr, ncol = 0L),
    eff_idx = match(spec$effect_col, colnames(X_full)),
    s_idx0_full = s_fac - 1L,
    n_s_full = max(s_fac),
    i_idx0_full = i_fac - 1L,
    n_i_full = max(i_fac),
    t_idx0_full = if (has_t) (t_fac - 1L) else integer(0),
    n_t_full = if (has_t) max(t_fac) else 0L,
    ws_env = new.env(parent = emptyenv())
  )
}

get_cached_cpp_ws <- function(cache, re_structure) {
  env <- cache$ws_env
  if (is.null(env)) return(NULL)
  key <- if (re_structure == "intercept_only") "int" else "slope"
  ws <- env[[key]]
  if (!is.null(ws) && cpp_glmm_workspace_valid(ws)) return(ws)
  Z_s <- if (re_structure == "intercept_only") cache$Z_int else cache$zs_dgp
  Z_i <- if (re_structure == "intercept_only") cache$Z_int else cache$zi_dgp
  Z_t <- if (cache$has_t) {
    if (re_structure == "intercept_only") cache$Z_int else cache$zt_dgp
  } else {
    cache$Z_empty
  }
  ws <- cpp_create_glmm_workspace(
    X = cache$X_full,
    subj_idx = cache$s_idx0_full,
    Z_subj = Z_s,
    n_subj = cache$n_s_full,
    item_idx = cache$i_idx0_full,
    Z_item = Z_i,
    n_item = cache$n_i_full,
    testlet_idx = cache$t_idx0_full,
    Z_testlet = Z_t,
    n_testlet = cache$n_t_full
  )
  env[[key]] <- ws
  ws
}

#' Zero-Allocation Response Simulation Matching simulate_dgp_data() Bit-for-Bit
#' @noRd
simulate_fast_y_keep <- function(cache) {
  u_s <- matrix(stats::rnorm(cache$ns * cache$p_s), nrow = cache$ns, ncol = cache$p_s) %*% cache$t_Ls
  u_i <- matrix(stats::rnorm(cache$ni * cache$p_i), nrow = cache$ni, ncol = cache$p_i) %*% cache$t_Li

  eta <- if (isTRUE(cache$scalar_dgp)) {
    cache$mu_trial + u_s[cache$sid, 1L] + u_i[cache$iid, 1L]
  } else {
    cache$mu_trial +
      rowSums(cache$zs_dgp * u_s[cache$sid, , drop = FALSE]) +
      rowSums(cache$zi_dgp * u_i[cache$iid, , drop = FALSE])
  }

  if (cache$has_t) {
    u_t <- matrix(stats::rnorm(cache$nt * cache$p_t), nrow = cache$nt, ncol = cache$p_t) %*% cache$t_Lt
    eta <- eta + rowSums(cache$zt_dgp * u_t[cache$tid, , drop = FALSE]) +
      trial_testlet_effect(cache, cache$ns, cache$nt, cache$sid, cache$tid)
  }

  prob <- cache$chance_dgp + (1.0 - cache$chance_dgp) * stats::plogis(eta)
  y <- stats::rbinom(length(prob), size = 1L, prob = prob)

  keep <- observed_trials(cache$design, cache$sid, cache$ns)
  if (!is.null(keep)) {
    list(y = y[keep], keep = keep, nobs = sum(keep))
  } else {
    list(y = y, keep = NULL, nobs = cache$n_tr)
  }
}

#' Fit a Cached Simulated Response Vector Directly via C++ Laplace
#' @noRd
fit_via_cpp_cached <- function(
  sim_dat, cache, spec, regularized, re_structure, fast_stage
) {
  if (re_structure == "intercept_only") {
    Z_s_full <- cache$Z_int
    Z_i_full <- cache$Z_int
    Z_t_full <- if (cache$has_t) cache$Z_int else cache$Z_empty
    corr_mode <- 0L
  } else {
    Z_s_full <- cache$zs_dgp
    Z_i_full <- cache$zi_dgp
    Z_t_full <- if (cache$has_t) cache$zt_dgp else cache$Z_empty
    corr_mode <- if (re_structure == "maximal") 1L else 0L
  }

  ws_ptr <- NULL
  if (is.null(sim_dat$keep)) {
    X <- cache$X_full
    Z_s <- Z_s_full
    Z_i <- Z_i_full
    Z_t <- Z_t_full
    s_idx0 <- cache$s_idx0_full
    n_s <- cache$n_s_full
    i_idx0 <- cache$i_idx0_full
    n_i <- cache$n_i_full
    t_idx0 <- cache$t_idx0_full
    n_t <- cache$n_t_full
    ws_ptr <- get_cached_cpp_ws(cache, re_structure)
  } else {
    kp <- sim_dat$keep
    if (length(unique(cache$cell_id[kp])) < ncol(cache$X_full)) {
      return(make_failed_fit("insufficient_data", "Observed trials omit required design cells."))
    }
    X <- cache$X_full[kp, , drop = FALSE]
    Z_s <- Z_s_full[kp, , drop = FALSE]
    Z_i <- Z_i_full[kp, , drop = FALSE]
    Z_t <- Z_t_full[kp, , drop = FALSE]
    cs <- compact_idx0(cache$sid[kp])
    s_idx0 <- cs$idx0
    n_s <- cs$n
    ci <- compact_idx0(cache$iid[kp])
    i_idx0 <- ci$idx0
    n_i <- ci$n
    if (ncol(Z_t) > 0L) {
      ct <- compact_idx0(cache$tid[kp])
      n_t <- ct$n
      t_idx0 <- ct$idx0
    } else {
      n_t <- 0L
      t_idx0 <- integer(0)
    }
  }

  warm_cpp <- NULL
  fit_response <- function(y) {
    res <- run_cpp_glmm_core(
      y = y,
      X = X,
      s_idx0 = s_idx0,
      Z_s = Z_s,
      n_s = n_s,
      i_idx0 = i_idx0,
      Z_i = Z_i,
      n_i = n_i,
      t_idx0 = t_idx0,
      Z_t = Z_t,
      n_t = n_t,
      chance = cache$chance_fit,
      corr_mode = corr_mode,
      regularized = regularized,
      spec = spec,
      fast_stage = fast_stage,
      eff_idx = cache$eff_idx,
      warm = warm_cpp,
      ws_ptr = ws_ptr
    )
    if (is.null(warm_cpp) && isTRUE(res$converged)) warm_cpp <<- res
    res
  }
  cpp_boot_ctx <- list(
    X = X, s_idx0 = s_idx0, Z_s = Z_s, n_s = n_s,
    i_idx0 = i_idx0, Z_i = Z_i, n_i = n_i,
    t_idx0 = t_idx0, Z_t = Z_t, n_t = n_t,
    corr_mode = corr_mode, regularized = regularized,
    eff_idx = cache$eff_idx, ws_ptr = ws_ptr
  )
  cid <- if (!regularized) {
    if (is.null(sim_dat$keep)) cache$cell_id else cache$cell_id[sim_dat$keep]
  }
  summarize_probability_fit(
    sim_dat$y, fit_response, X, Z_s, Z_i, Z_t, s_idx0 + 1L, i_idx0 + 1L, t_idx0 + 1L,
    spec, re_structure, cache$chance_fit, fast_stage, cpp_boot_ctx = cpp_boot_ctx,
    unpenalized_cells = cid
  )
}

validate_custom_data <- function(df_b, schedule) {
  if (!is.data.frame(df_b)) {
    cli::cli_abort("{.arg data_fn} must return a {.cls data.frame}.")
  }
  if (!("y" %in% names(df_b))) {
    cli::cli_abort("{.arg data_fn} must return a {.cls data.frame} containing a binary {.field y} column.")
  }
  if (nrow(df_b) > 0L) {
    y_num <- as.numeric(df_b$y)
    if (any(!is.finite(y_num)) || any(!(y_num %in% c(0, 1)))) {
      cli::cli_abort("{.field y} returned by {.arg data_fn} must contain only 0 and 1 values.")
    }
    df_b$y <- as.integer(y_num)
  }
  if (!("subj_id" %in% names(df_b)) && ("subj" %in% names(df_b))) {
    df_b$subj_id <- as.integer(as.factor(df_b$subj))
  }
  if (!("item_id" %in% names(df_b)) && ("item" %in% names(df_b))) {
    df_b$item_id <- as.integer(as.factor(df_b$item))
  }
  if (!("subj" %in% names(df_b)) && ("subj_id" %in% names(df_b))) {
    df_b$subj <- factor(sprintf("S%04d", as.integer(df_b$subj_id)))
  }
  if (!("item" %in% names(df_b)) && ("item_id" %in% names(df_b))) {
    df_b$item <- factor(sprintf("I%04d", as.integer(df_b$item_id)))
  }
  req_cols <- c("subj_id", "item_id", "cell_id")
  missing_cols <- setdiff(req_cols, names(df_b))
  if (length(missing_cols) > 0L) {
    cli::cli_abort(
      "{.arg data_fn} output is missing required column(s): {.val {missing_cols}}. Start from {.code schedule$trials} or {.fn simulate_dgp_data}."
    )
  }
  df_b
}

#' Execute One Monte Carlo Replication (Zero-Alloc C++ Path or External Backend)
#' @noRd
execute_single_replication <- function(
  spec,
  schedule,
  cpp_cache,
  engine_use,
  re_primary,
  fast_stage,
  do_compare,
  data_fn = NULL
) {
  if (!is.null(cpp_cache) && is.null(data_fn)) {
    sim_dat <- simulate_fast_y_keep(cpp_cache)
    reg_prim <- (engine_use == "bayes_laplace")
    primary <- fit_via_cpp_cached(
      sim_dat, cpp_cache, spec = spec,
      regularized = reg_prim, re_structure = re_primary, fast_stage = fast_stage
    )
    is_sing <- isTRUE(primary$singular)
    sing_init <- if (engine_supports_singularity(engine_use)) is_sing else NA

    if (is_sing && spec$on_singular == "reduce_corr" && re_primary == "maximal") {
      primary <- fit_via_cpp_cached(
        sim_dat, cpp_cache, spec = spec,
        regularized = reg_prim, re_structure = "uncorrelated", fast_stage = fast_stage
      )
    } else if (is_sing && spec$on_singular == "switch_bayes") {
      primary <- fit_via_cpp_cached(
        sim_dat, cpp_cache, spec = spec,
        regularized = TRUE, re_structure = re_primary, fast_stage = fast_stage
      )
    }

    reduced <- if (do_compare) {
      if (is_sing || !isTRUE(primary$converged)) {
        fit_via_cpp_cached(
          sim_dat, cpp_cache, spec = spec,
          regularized = reg_prim, re_structure = "uncorrelated", fast_stage = fast_stage
        )
      } else {
        primary
      }
    } else {
      NULL
    }

    subj <- cpp_cache$sid
    item <- cpp_cache$iid
    if (!is.null(sim_dat$keep)) {
      subj <- subj[sim_dat$keep]
      item <- item[sim_dat$keep]
    }
    patterns <- replication_patterns(sim_dat$y, subj, item, schedule$n_subj, schedule$n_item,
                                     primary, spec$scale)
    return(list(primary = primary, sing_init = sing_init, reduced = reduced, nobs = sim_dat$nobs,
                patterns = patterns))
  }

  df_b <- if (is.null(data_fn)) {
    simulate_dgp_data(spec$dgp, schedule = schedule)
  } else {
    validate_custom_data(data_fn(schedule, spec), schedule)
  }
  primary <- fit_single_dataset(
    df_b, spec = spec, engine_use = engine_use,
    re_structure = re_primary, fast_stage = fast_stage
  )
  is_sing <- isTRUE(primary$singular)
  sing_init <- if (engine_supports_singularity(engine_use)) is_sing else NA

  if (is_sing && spec$on_singular == "reduce_corr" && re_primary == "maximal") {
    primary <- fit_single_dataset(
      df_b, spec = spec, engine_use = engine_use,
      re_structure = "uncorrelated", fast_stage = fast_stage
    )
  } else if (is_sing && spec$on_singular == "switch_bayes") {
    primary <- fit_single_dataset(
      df_b, spec = spec, engine_use = "bayes_laplace",
      re_structure = re_primary, fast_stage = fast_stage
    )
  }

  reduced <- if (do_compare) {
    if (is_sing || !isTRUE(primary$converged)) {
      fit_single_dataset(
        df_b, spec = spec, engine_use = engine_use,
        re_structure = "uncorrelated", fast_stage = fast_stage
      )
    } else {
      primary
    }
  } else {
    NULL
  }

  patterns <- replication_patterns(df_b$y, df_b$subj_id, df_b$item_id, schedule$n_subj,
                                   schedule$n_item, primary, spec$scale)
  list(primary = primary, sing_init = sing_init, reduced = reduced, nobs = nrow(df_b),
       patterns = patterns)
}

#' Evaluate Monte Carlo Replications at a Single (n_subj, n_item) Design Point
#'
#' Runs `nsim` independent Monte Carlo simulations with reproducible
#' L'Ecuyer-CMRG substreams and computes both unconditional (`effective_*`) and
#' valid-fit conditional (`conditional_*`) inferential/precision metrics,
#' percentage-point accuracy difference intervals, convergence/error breakdowns,
#' and Morris et al. (2019) Monte Carlo uncertainty intervals.
#'
#' @param spec A `"glmmplanr_spec"` object.
#' @param n_subj Number of subjects.
#' @param n_item Number of items.
#' @param nsim Number of Monte Carlo replications.
#' @param engine_use Estimation backend to execute (`spec$stage1_engine` or
#'   `spec$engine`).
#' @param fast_stage Logical; whether to use Stage-1 fast tolerances.
#' @param stream_seeds List of `nsim` L'Ecuyer-CMRG substream seeds (or integer
#'   master seed).
#' @param data_fn Optional custom data-generating function `function(schedule, spec)`
#'   returning a trial `data.frame` (starting from `schedule$trials` or
#'   [simulate_dgp_data()] with binary `y`, `subj_id`, `item_id`, and `cell_id`)
#'   for evaluating structural misspecification such as ability-dependent
#'   dropout, item-specific guessing, or non-normal random effects. Default
#'   `NULL` uses the built-in parametric DGP in `spec$dgp`.
#' @param verbose Logical; display a live `cli` progress bar with ETA across
#'   replications (default `FALSE`).
#' @return A single-row `data.frame` summarizing all evaluation and diagnostic
#'   layers. `power` concerns the target-scale interval excluding zero.
#'   `effective_logit_power` counts significant conditional coefficients among all
#'   replications. Error counts are in `attr(res, "error_summary")`.
#'   `pattern_nsim` and `pattern_n_invalid` count diagnosed replications and
#'   invalid target fits. For each subject/item pattern (all correct, all
#'   incorrect, unobserved), `n_rep_*` counts replications with any such unit
#'   and `n_invalid_*` counts those with an invalid target fit. `max_obs_*`
#'   gives the largest observed response count for an all-correct or
#'   all-incorrect unit (0 if absent). See [response_patterns()] for unit counts.
#' @examples
#' des <- design_crossed(list(condition = c("a", "b")))
#' dgp <- calibrate_dgp(des, c(a = 0.4, b = 0.6), subj_sd = 0.3, item_sd = 0.4)
#' spec <- specify_target(dgp, "condition", criterion = "power", engine = "bayes_laplace")
#' # Small run to inspect the output.
#' r <- evaluate_design_point(spec, 16, 8, nsim = 5L)
#' r[c("nsim", "n_valid", "power")]
#' @export
evaluate_design_point <- function(
  spec,
  n_subj,
  n_item,
  nsim = 100L,
  engine_use = spec$engine,
  fast_stage = FALSE,
  stream_seeds = 1L,
  data_fn = NULL,
  verbose = FALSE
) {
  if (uses_box_muller()) {
    return(isolated_seeded_call("evaluate_design_point", environment(), match.call()))
  }
  if (!inherits(spec, "glmmplanr_spec")) {
    cli::cli_abort("{.arg spec} must be a {.cls glmmplanr_spec} object.")
  }
  if (!is.null(data_fn) && !is.function(data_fn)) {
    cli::cli_abort("{.arg data_fn} must be NULL or a function of {.code (schedule, spec)}.")
  }
  nsim <- as.integer(nsim)
  if (!is.list(stream_seeds)) {
    base_st <- allocate_stage_point_streams(stream_seeds, stage_slot = 1L, grid_ids = 1L, stride = 1L)
    stream_seeds <- derive_rep_substreams(base_st[[1L]], nsim = nsim, offset = 0L)
  }

  schedule <- generate_design_schedule(spec$dgp$design, n_subj = n_subj, n_item = n_item)
  re_primary <- if (fast_stage && !is.null(spec$stage1_re)) spec$stage1_re else spec$analysis_re
  cpp_cache <- if (is.null(data_fn) && engine_use %in% c("cpp_laplace", "bayes_laplace")) {
    prepare_cpp_schedule_cache(spec, schedule)
  } else {
    NULL
  }

  est_vec <- rep(NA_real_, nsim)
  se_vec <- rep(NA_real_, nsim)
  pest_vec <- rep(NA_real_, nsim)
  pse_vec <- rep(NA_real_, nsim)
  sing_vec <- rep(NA, nsim)
  conv_vec <- rep(FALSE, nsim)
  status_vec <- rep("engine_error", nsim)
  err_vec <- rep("", nsim)
  nobs_vec <- rep(0L, nsim)
  patterns <- empty_pattern_totals

  do_compare <- (spec$on_singular == "flag_and_compare") && (re_primary == "maximal")
  red_est_vec <- rep(NA_real_, nsim)
  red_se_vec <- rep(NA_real_, nsim)
  red_pest_vec <- rep(NA_real_, nsim)
  red_pse_vec <- rep(NA_real_, nsim)
  red_conv_vec <- rep(FALSE, nsim)
  red_status_vec <- rep("engine_error", nsim)

  pb_id <- if (isTRUE(verbose) && nsim >= 2L) {
    cli::cli_progress_bar(
      format = "Simulating ({n_subj}x{n_item}) {cli::pb_bar} {cli::pb_current}/{cli::pb_total} [{cli::pb_percent} | ETA: {cli::pb_eta}]",
      total = nsim,
      clear = FALSE
    )
  } else {
    NULL
  }
  if (!is.null(pb_id)) on.exit(cli::cli_progress_done(id = pb_id), add = TRUE)

  for (b in seq_len(nsim)) {
    rep_out <- with_rng_stream(
      stream_seeds[[b]],
      execute_single_replication(
        spec = spec,
        schedule = schedule,
        cpp_cache = cpp_cache,
        engine_use = engine_use,
        re_primary = re_primary,
        fast_stage = fast_stage,
        do_compare = do_compare,
        data_fn = data_fn
      )
    )

    est_vec[b] <- rep_out$primary$est
    se_vec[b] <- rep_out$primary$se
    pest_vec[b] <- rep_out$primary$prob_est
    pse_vec[b] <- rep_out$primary$prob_se
    sing_vec[b] <- rep_out$sing_init
    conv_vec[b] <- isTRUE(rep_out$primary$converged)
    status_vec[b] <- rep_out$primary$status
    err_vec[b] <- rep_out$primary$error_reason
    nobs_vec[b] <- rep_out$nobs
    patterns <- combine_pattern_totals(patterns, rep_out$patterns)

    if (do_compare && !is.null(rep_out$reduced)) {
      red_est_vec[b] <- rep_out$reduced$est
      red_se_vec[b] <- rep_out$reduced$se
      red_pest_vec[b] <- rep_out$reduced$prob_est
      red_pse_vec[b] <- rep_out$reduced$prob_se
      red_conv_vec[b] <- isTRUE(rep_out$reduced$converged)
      red_status_vec[b] <- rep_out$reduced$status
    }
    if (!is.null(pb_id)) cli::cli_progress_update(id = pb_id)
  }

  out <- summarize_mc_metrics(
    est = est_vec,
    se = se_vec,
    prob_est = pest_vec,
    prob_se = pse_vec,
    singular = sing_vec,
    converged = conv_vec,
    status = status_vec,
    error_reasons = err_vec,
    nobs = nobs_vec,
    spec = spec,
    n_subj = schedule$n_subj,
    n_item = schedule$n_item,
    nsim = nsim,
    engine_use = engine_use
  )

  if (do_compare) {
    red_sum <- summarize_mc_metrics(
      est = red_est_vec,
      se = red_se_vec,
      prob_est = red_pest_vec,
      prob_se = red_pse_vec,
      singular = rep(FALSE, nsim),
      converged = red_conv_vec,
      status = red_status_vec,
      error_reasons = rep("", nsim),
      nobs = nobs_vec,
      spec = spec,
      n_subj = schedule$n_subj,
      n_item = schedule$n_item,
      nsim = nsim,
      engine_use = engine_use
    )
    out$reduced_power <- red_sum$power
    out$reduced_ci_width_mean <- red_sum$ci_width_mean
    out$reduced_coverage <- red_sum$coverage
  }

  out[pattern_columns] <- as.list(patterns)
  out
}

valid_target_fit <- function(est, se, converged, status, prob_est, prob_se, scale) {
  valid <- converged & status == "converged" & is.finite(est) & is.finite(se) & se > 0
  if (scale %in% c("pp", "prob")) {
    if (is.null(prob_est) || is.null(prob_se)) return(rep(FALSE, length(valid)))
    valid <- valid & is.finite(prob_est) & is.finite(prob_se) & prob_se > 0
  }
  valid
}

summarize_mc_metrics <- function(
  est,
  se,
  prob_est = NULL,
  prob_se = NULL,
  singular,
  converged,
  status = NULL,
  error_reasons = NULL,
  nobs = NULL,
  spec,
  n_subj,
  n_item,
  nsim,
  engine_use
) {
  des <- spec$dgp$design
  z_crit <- stats::qnorm(1 - spec$alpha / 2)
  true_b <- spec$true_beta
  true_p <- spec$true_prob_diff

  if (is.null(status)) {
    status <- ifelse(converged, "converged", "non_converged")
  }
  if (is.null(error_reasons)) {
    error_reasons <- rep("", length(est))
  }

  valid_logit <- valid_target_fit(est, se, converged, status, prob_est, prob_se, "logit")
  logit_detection <- sum(abs(est[valid_logit] / se[valid_logit]) > z_crit) / nsim
  valid <- valid_target_fit(est, se, converged, status, prob_est, prob_se, spec$scale)

  n_valid <- sum(valid)
  valid_rate <- n_valid / nsim
  failure_rate <- 1.0 - valid_rate
  max_iter_rate <- mean(status == "max_iter_finite")
  error_rate <- mean(status == "engine_error")
  singular_rate <- if (all(is.na(singular))) NA_real_ else mean(singular, na.rm = TRUE)

  n_sched <- compute_scheduled_trials(des, n_subj, n_item)
  n_subj_expected <- n_subj * (1 - dropout_rate(des)) * (1 - des$attrition^(n_sched / n_subj))
  n_exp <- compute_expected_trials(des, n_subj, n_item)
  n_obs_mean <- if (!is.null(nobs) && length(nobs) > 0L) mean(nobs) else as.numeric(n_exp)

  # Summarize backend error / non-convergence reasons
  non_empty_err <- error_reasons[!valid & nzchar(error_reasons)]
  if (length(non_empty_err) > 0L) {
    err_tab <- sort(table(non_empty_err), decreasing = TRUE)
    top_reason <- sprintf("[%d/%d] %s", err_tab[1L], nsim, names(err_tab)[1L])
  } else if (n_valid < nsim) {
    st_tab <- sort(table(status[!valid]), decreasing = TRUE)
    top_reason <- sprintf("[%d/%d] %s", st_tab[1L], nsim, names(st_tab)[1L])
    err_tab <- st_tab
  } else {
    top_reason <- ""
    err_tab <- table(character(0))
  }

  if (n_valid == 0L) {
    mc0 <- wilson_mc_ci(0L, nsim, conf_level = spec$mc_conf_level)
    out0 <- data.frame(
      n_subj = as.integer(n_subj),
      n_item = as.integer(n_item),
      n_trials_scheduled = n_sched,
      n_trials_expected = n_exp,
      n_subj_expected = n_subj_expected,
      n_trials = n_exp,
      n_trials_obs_mean = n_obs_mean,
      nsim = as.integer(nsim),
      n_valid = 0L,
      valid_rate = 0.0,
      engine = engine_use,
      scale = spec$scale,
      power = if (spec$failure_policy == "effective") 0.0 else NA_real_,
      effective_power = 0.0,
      effective_logit_power = logit_detection,
      conditional_power = NA_real_,
      power_mcse = 0.0,
      power_mc_lo = unname(mc0["lo"]),
      power_mc_hi = unname(mc0["hi"]),
      ci_width_mean = NA_real_,
      ci_width_q80 = NA_real_,
      logit_ci_width_mean = NA_real_,
      prob_diff_true_pp = spec$true_prob_diff_pp,
      prob_diff_est_pp = NA_real_,
      prob_ci_width_mean_pp = NA_real_,
      prob_coverage = NA_real_,
      aipe_prob = if (!is.null(spec$ci_width) && spec$failure_policy == "effective") 0.0 else NA_real_,
      effective_aipe_prob = if (!is.null(spec$ci_width)) 0.0 else NA_real_,
      conditional_aipe_prob = NA_real_,
      aipe_mc_lo = if (!is.null(spec$ci_width)) unname(mc0["lo"]) else NA_real_,
      aipe_mc_hi = if (!is.null(spec$ci_width)) unname(mc0["hi"]) else NA_real_,
      k_joint = if (is.null(spec$ci_width)) NA_integer_ else 0L,
      joint_prob = if (!is.null(spec$ci_width) && spec$failure_policy == "effective") 0.0 else NA_real_,
      joint_mc_lo = if (!is.null(spec$ci_width) && spec$failure_policy == "effective") unname(mc0["lo"]) else NA_real_,
      joint_mc_hi = if (!is.null(spec$ci_width) && spec$failure_policy == "effective") unname(mc0["hi"]) else NA_real_,
      rope_excl_prob = if (spec$failure_policy == "effective") 0.0 else NA_real_,
      effective_rope_prob = 0.0,
      conditional_rope_prob = NA_real_,
      rope_mc_lo = unname(mc0["lo"]),
      rope_mc_hi = unname(mc0["hi"]),
      criterion_met = FALSE,
      decision_status = compute_decision_status(spec, spec$failure_policy == "effective",
                                                mc0, mc0, mc0, mc0),
      curtailed = FALSE,
      singular_rate = singular_rate,
      nonconv_rate = failure_rate,
      max_iter_rate = max_iter_rate,
      error_rate = error_rate,
      top_error_reason = top_reason,
      coverage = NA_real_,
      power_nonsingular = NA_real_,
      coverage_nonsingular = NA_real_,
      bias = NA_real_,
      rmse = NA_real_,
      type_m_error = NA_real_,
      type_s_error = NA_real_,
      stringsAsFactors = FALSE
    )
    attr(out0, "error_summary") <- err_tab
    return(out0)
  }

  e_v <- est[valid]
  s_v <- se[valid]
  logit_lo <- e_v - z_crit * s_v
  logit_hi <- e_v + z_crit * s_v
  logit_width_v <- 2 * z_crit * s_v

  has_prob <- !is.null(prob_est) && !is.null(prob_se) && all(is.finite(prob_est[valid])) &&
    all(is.finite(prob_se[valid]))
  if (has_prob) {
    pe_v <- prob_est[valid]
    ps_v <- prob_se[valid]
    prob_lo <- pe_v - z_crit * ps_v
    prob_hi <- pe_v + z_crit * ps_v
    prob_width_v <- 2 * z_crit * ps_v
    prob_diff_est_pp <- mean(pe_v) * 100.0
    prob_ci_width_mean_pp <- mean(prob_width_v) * 100.0
    prob_cov <- mean((prob_lo <= true_p) & (true_p <= prob_hi))
  } else {
    pe_v <- rep(NA_real_, n_valid)
    prob_lo <- rep(NA_real_, n_valid)
    prob_hi <- rep(NA_real_, n_valid)
    prob_width_v <- rep(NA_real_, n_valid)
    prob_diff_est_pp <- NA_real_
    prob_ci_width_mean_pp <- NA_real_
    prob_cov <- NA_real_
  }

  # Scale-matched interval bounds & widths for AIPE and ROPE evaluation
  if (spec$scale == "logit") {
    target_lo <- logit_lo
    target_hi <- logit_hi
    target_width_v <- logit_width_v
    true_target <- true_b
    est_target_v <- e_v
  } else if (spec$scale == "prob") {
    target_lo <- prob_lo
    target_hi <- prob_hi
    target_width_v <- prob_width_v
    true_target <- true_p
    est_target_v <- pe_v
  } else {
    target_lo <- prob_lo * 100.0
    target_hi <- prob_hi * 100.0
    target_width_v <- prob_width_v * 100.0
    true_target <- true_p * 100.0
    est_target_v <- pe_v * 100.0
  }

  sig_v <- target_lo > 0 | target_hi < 0
  k_sig <- sum(sig_v)
  effective_power <- k_sig / nsim
  conditional_power <- k_sig / n_valid

  if (!is.null(spec$ci_width)) {
    k_aipe <- sum(target_width_v <= spec$ci_width)
    effective_aipe <- k_aipe / nsim
    conditional_aipe <- k_aipe / n_valid
  } else {
    k_aipe <- NA_integer_
    effective_aipe <- NA_real_
    conditional_aipe <- NA_real_
  }

  rope_excl_v <- (target_lo > spec$rope[2L]) | (target_hi < spec$rope[1L])
  k_rope <- sum(rope_excl_v)
  effective_rope <- k_rope / nsim
  conditional_rope <- k_rope / n_valid

  if (spec$failure_policy == "effective") {
    denom_mc <- nsim
    power_use <- effective_power
    aipe_use <- effective_aipe
    rope_use <- effective_rope
    valid_gate_ok <- TRUE
  } else {
    denom_mc <- n_valid
    power_use <- conditional_power
    aipe_use <- conditional_aipe
    rope_use <- conditional_rope
    valid_gate_ok <- (valid_rate >= spec$min_valid_rate)
  }

  power_mcse <- sqrt(power_use * (1.0 - power_use) / denom_mc)
  mc_pow <- wilson_mc_ci(k_sig, denom_mc, conf_level = spec$mc_conf_level)
  mc_aipe <- if (!is.null(spec$ci_width)) {
    wilson_mc_ci(k_aipe, denom_mc, conf_level = spec$mc_conf_level)
  } else {
    c(lo = NA_real_, hi = NA_real_)
  }
  k_joint <- if (is.null(spec$ci_width)) NA_integer_ else sum(sig_v & target_width_v <= spec$ci_width)
  joint_prob <- k_joint / denom_mc
  mc_joint <- wilson_mc_ci(k_joint, denom_mc, conf_level = spec$mc_conf_level)
  mc_rope <- wilson_mc_ci(k_rope, denom_mc, conf_level = spec$mc_conf_level)

  ci_width_mean <- mean(target_width_v)
  ci_width_q80 <- unname(stats::quantile(target_width_v, probs = spec$ci_assurance))
  logit_ci_width_mean <- mean(logit_width_v)

  coverage <- mean((target_lo <= true_target) & (true_target <= target_hi))
  bias <- mean(est_target_v - true_target)
  rmse <- sqrt(mean((est_target_v - true_target)^2))

  if (all(is.na(singular))) {
    power_nonsing <- NA_real_
    cov_ns <- NA_real_
  } else {
    not_sing <- is.na(singular) | (!singular)
    nonsing_valid <- valid & not_sing
    if (any(nonsing_valid)) {
      ns <- not_sing[valid]
      power_nonsing <- mean(sig_v[ns])
      cov_ns <- mean(target_lo[ns] <= true_target & true_target <= target_hi[ns])
    } else {
      power_nonsing <- NA_real_
      cov_ns <- NA_real_
    }
  }

  if (any(sig_v) && abs(true_target) > 1e-8) {
    e_sig <- est_target_v[sig_v]
    type_m <- mean(abs(e_sig) / abs(true_target))
    type_s <- mean(sign(e_sig) != sign(true_target))
  } else {
    type_m <- NA_real_
    type_s <- NA_real_
  }

  criterion_met <- valid_gate_ok && check_criterion_met(
    spec = spec,
    power = power_use,
    aipe_prob = aipe_use,
    rope_excl_prob = rope_use, joint_prob = joint_prob
  )

  decision_status <- compute_decision_status(
    spec = spec,
    valid_gate_ok = valid_gate_ok,
    mc_pow = mc_pow,
    mc_aipe = mc_aipe,
    mc_rope = mc_rope, mc_joint = mc_joint
  )

  out <- data.frame(
    n_subj = as.integer(n_subj),
    n_item = as.integer(n_item),
    n_trials_scheduled = n_sched,
    n_trials_expected = n_exp,
    n_subj_expected = n_subj_expected,
    n_trials = n_exp,
    n_trials_obs_mean = n_obs_mean,
    nsim = as.integer(nsim),
    n_valid = as.integer(n_valid),
    valid_rate = valid_rate,
    engine = engine_use,
    scale = spec$scale,
    power = power_use,
    effective_power = effective_power,
    effective_logit_power = logit_detection,
    conditional_power = conditional_power,
    power_mcse = power_mcse,
    power_mc_lo = unname(mc_pow["lo"]),
    power_mc_hi = unname(mc_pow["hi"]),
    ci_width_mean = ci_width_mean,
    ci_width_q80 = ci_width_q80,
    logit_ci_width_mean = logit_ci_width_mean,
    prob_diff_true_pp = spec$true_prob_diff_pp,
    prob_diff_est_pp = prob_diff_est_pp,
    prob_ci_width_mean_pp = prob_ci_width_mean_pp,
    prob_coverage = prob_cov,
    aipe_prob = aipe_use,
    effective_aipe_prob = effective_aipe,
    conditional_aipe_prob = conditional_aipe,
    aipe_mc_lo = unname(mc_aipe["lo"]),
    aipe_mc_hi = unname(mc_aipe["hi"]),
    k_joint = k_joint, joint_prob = joint_prob,
    joint_mc_lo = unname(mc_joint["lo"]), joint_mc_hi = unname(mc_joint["hi"]),
    rope_excl_prob = rope_use,
    effective_rope_prob = effective_rope,
    conditional_rope_prob = conditional_rope,
    rope_mc_lo = unname(mc_rope["lo"]),
    rope_mc_hi = unname(mc_rope["hi"]),
    criterion_met = criterion_met,
    decision_status = decision_status,
    curtailed = FALSE,
    singular_rate = singular_rate,
    nonconv_rate = failure_rate,
    max_iter_rate = max_iter_rate,
    error_rate = error_rate,
    top_error_reason = top_reason,
    coverage = coverage,
    power_nonsingular = power_nonsing,
    coverage_nonsingular = cov_ns,
    bias = bias,
    rmse = rmse,
    type_m_error = type_m,
    type_s_error = type_s,
    stringsAsFactors = FALSE
  )
  attr(out, "error_summary") <- err_tab
  out
}

check_criterion_met <- function(spec, power, aipe_prob, rope_excl_prob, joint_prob = NA_real_) {
  if (is.na(power)) return(FALSE)
  switch(
    spec$criterion,
    power = (power >= spec$power_min),
    ci_width = (!is.na(aipe_prob) && aipe_prob >= spec$ci_assurance),
    both = (power >= spec$power_min) && (!is.na(aipe_prob) && aipe_prob >= spec$ci_assurance),
    joint = (!is.na(joint_prob) && joint_prob >= spec$power_min),
    rope = (!is.na(rope_excl_prob) && rope_excl_prob >= spec$power_min)
  )
}

compute_decision_status <- function(spec, valid_gate_ok, mc_pow, mc_aipe, mc_rope, mc_joint = NULL) {
  if (!isTRUE(valid_gate_ok)) return("fail")
  switch(
    spec$criterion,
    power = classify_mc_threshold(mc_pow["lo"], mc_pow["hi"], spec$power_min),
    ci_width = classify_mc_threshold(mc_aipe["lo"], mc_aipe["hi"], spec$ci_assurance),
    both = {
      d1 <- classify_mc_threshold(mc_pow["lo"], mc_pow["hi"], spec$power_min)
      d2 <- classify_mc_threshold(mc_aipe["lo"], mc_aipe["hi"], spec$ci_assurance)
      if (d1 == "fail" || d2 == "fail") {
        "fail"
      } else if (d1 == "pass" && d2 == "pass") {
        "pass"
      } else {
        "indeterminate"
      }
    },
    joint = classify_mc_threshold(mc_joint["lo"], mc_joint["hi"], spec$power_min),
    rope = classify_mc_threshold(mc_rope["lo"], mc_rope["hi"], spec$power_min)
  )
}

can_define_null_target <- function(spec) {
  if (identical(spec$scale, "logit")) return(TRUE)
  w <- spec$contrast_weights
  all(is.finite(w)) && any(w > 0) && any(w < 0)
}

null_target_spec <- function(spec) {
  dgp <- spec$dgp
  X <- dgp$design$X_cells
  if (spec$scale == "logit") {
    dgp$beta[[spec$effect_col]] <- 0
  } else {
    w <- spec$contrast_weights
    if (!all(is.finite(w)) || !any(w > 0) || !any(w < 0)) {
      cli::cli_abort(c(
        "A probability-scale null requires finite contrast weights of both signs (a difference between conditions).",
        "i" = "For a single-condition or absolute probability precision target, validate interval coverage and width attainment under the planned condition via {.fn validate_design_point} or {.fn evaluate_design_point}."
      ))
    }
    lower <- sum(w[w < 0]) + dgp$chance * sum(w[w > 0])
    upper <- sum(w[w > 0]) + dgp$chance * sum(w[w < 0])
    if (lower >= 0 || upper <= 0) {
      cli::cli_abort("Zero is outside the attainable probability contrast for this guessing floor and weights.")
    }
    # Tilt cell logits along the contrast, preserving orthogonal logit contrasts.
    direction <- drop(qr.solve(X, w / max(abs(w))))
    contrast_at <- function(shift) {
      trial <- dgp
      trial$beta <- dgp$beta + shift * direction
      true_prob_effect(trial, w, spec$prob_estimand)
    }
    limit <- 1
    while (limit < 2^30 && (contrast_at(-limit) > 0 || contrast_at(limit) < 0)) limit <- 2 * limit
    if (contrast_at(-limit) > 0 || contrast_at(limit) < 0) {
      cli::cli_abort("Zero is outside the attainable probability contrast for this guessing floor and weights.")
    }
    shift <- stats::uniroot(contrast_at, c(-limit, limit), tol = 1e-12)$root
    dgp$beta <- dgp$beta + shift * direction
  }
  cells <- dgp$design$cells$cell_name
  dgp$mu_cells <- stats::setNames(drop(X %*% dgp$beta), cells)
  dgp$target_cell_probs <- dgp$achieved_cell_probs <-
    stats::setNames(dgp_cell_probabilities(dgp, "marginal"), cells)
  spec$dgp <- dgp
  spec$true_beta <- unname(dgp$beta[[spec$effect_col]])
  spec$true_prob_diff <- true_prob_effect(dgp, spec$contrast_weights, spec$prob_estimand)
  spec$true_prob_diff_pp <- 100 * spec$true_prob_diff
  spec$true_effect_target <- switch(spec$scale, logit = spec$true_beta,
    prob = spec$true_prob_diff, pp = spec$true_prob_diff_pp)
  if (abs(spec$true_effect_target) > 1e-8) cli::cli_abort("Failed to calibrate a zero target contrast.")
  spec
}

#' Validate the Planned Analysis at a Fixed Design
#'
#' Evaluates the complete analysis specified by `spec` under the planned
#' effect and a null on the same target scale, using disjoint random streams.
#' On the probability scale the null shifts cell logits along the supplied
#' contrast until its probability difference is zero, holding the covariance
#' structure fixed. On the logit scale it sets the focal coefficient to zero.
#' @inheritParams evaluate_design_point
#' @param seed Integer seed for this independent fixed-size validation.
#' @return A data frame of class `glmmplanr_validation` with two rows
#'   (`"null"` and `"planned"`) for contrasts
#'   or logit targets, or one `"planned"` row when `spec` targets an absolute
#'   probability (e.g., a single-condition design on `scale = "pp"` or `"prob"`
#'   where no zero-difference null exists), including the true target, detection and
#'   its Wilson interval under the chosen failure policy, valid-fit coverage
#'   and its Wilson interval, and the valid-fit count. Null detection is the
#'   Type I error. Response-pattern counts follow [evaluate_design_point()].
#'   Use [plot.glmmplanr_validation()] to display these results. Prespecify
#'   validation tolerances and replication counts.
#' @examples
#' des <- design_crossed(list(condition = c("a", "b")))
#' dgp <- calibrate_dgp(des, c(a = 0.4, b = 0.6), subj_sd = 0.3, item_sd = 0.4)
#' spec <- specify_target(dgp, "condition", criterion = "power", engine = "bayes_laplace")
#' r <- validate_design_point(spec, 16, 8, nsim = 5L)
#' r[c("condition", "nsim", "n_valid", "detection", "coverage")]
#' @export
validate_design_point <- function(spec, n_subj, n_item, nsim = 200L, seed = 2026L, data_fn = NULL, verbose = FALSE) {
  if (uses_box_muller()) {
    return(isolated_seeded_call("validate_design_point", environment(), match.call()))
  }
  if (!inherits(spec, "glmmplanr_spec")) cli::cli_abort("Supply a target specification.")
  if (!is.null(data_fn) && !is.function(data_fn)) {
    cli::cli_abort("{.arg data_fn} must be NULL or a function of {.code (schedule, spec)}.")
  }
  if (!is.numeric(nsim) || length(nsim) != 1L || !is.finite(nsim) ||
      nsim < 1 || nsim != floor(nsim)) cli::cli_abort("{.arg nsim} must be a positive integer.")
  bases_pair <- allocate_stage_point_streams(seed, stage_slot = 1L, grid_ids = 1:2, stride = 2L)
  has_null <- can_define_null_target(spec)
  if (has_null) {
    specs <- list(null = null_target_spec(spec), planned = spec)
    bases <- bases_pair
  } else {
    cli::cli_inform(c(
      "i" = "Target {.val {spec$effect_input}} on the {.val {spec$scale}} scale estimates an absolute probability rather than a between-condition difference.",
      "i" = "No zero-difference null applies; validating interval coverage and width attainment under the {.val planned} condition only."
    ))
    specs <- list(planned = spec)
    bases <- bases_pair[2L]
  }
  rows <- lapply(seq_along(specs), function(i) {
    sp <- specs[[i]]
    if (isTRUE(verbose)) {
      cli::cli_inform("Validating {.val {names(specs)[i]}} condition at ({n_subj} x {n_item})...")
    }
    ev <- evaluate_design_point(sp, n_subj, n_item, nsim,
      stream_seeds = derive_rep_substreams(bases[[i]], nsim), data_fn = data_fn, verbose = verbose)
    cover <- wilson_mc_ci(round(ev$coverage * ev$n_valid), ev$n_valid, sp$mc_conf_level)
    data.frame(condition = names(specs)[i], scale = sp$scale, true_target = sp$true_effect_target,
      nsim = nsim, n_valid = ev$n_valid, valid_rate = ev$valid_rate,
      detection = ev$power,
      detection_lo = if (is.na(ev$power)) NA_real_ else ev$power_mc_lo,
      detection_hi = if (is.na(ev$power)) NA_real_ else ev$power_mc_hi,
      coverage = ev$coverage, coverage_lo = unname(cover["lo"]), coverage_hi = unname(cover["hi"]),
      width_attainment = ev$aipe_prob, joint_attainment = ev$joint_prob, ev[pattern_columns])
  })
  design <- spec$dgp$design
  structure(do.call(rbind, rows), class = c("glmmplanr_validation", "data.frame"),
    validation_settings = list(n_subj = n_subj, n_item = n_item,
      items_per_subject = if (is.null(design$items_per_subject)) n_item else design$items_per_subject,
      has_null = has_null, alpha = spec$alpha, mc_conf_level = spec$mc_conf_level,
      failure_policy = spec$failure_policy, criterion = spec$criterion, power_min = spec$power_min,
      ci_width = spec$ci_width, ci_assurance = spec$ci_assurance))
}

#' Stress-Test Analysis Models Under Maximal DGP Mismatch
#'
#' Generates paired datasets from the maximal correlated data-generating process
#' (optionally under a null focal effect to inspect Type I error inflation) and
#' compares competing analysis strategies side by side on the logit scale.
#' Use [validate_design_point()] to validate the actual target and analysis.
#'
#' @inheritParams evaluate_design_point
#' @param dgp A `"glmmplanr_dgp"` object.
#' @param effect Focal fixed effect name (e.g., `"exposure:spacing"`).
#' @param include_null Logical; if `TRUE` (default), also evaluates empirical
#'   Type I error when the focal conditional effect is zeroed out while keeping
#'   all maximal random slopes active.
#' @param alpha Significance level (default `0.05`).
#' @param seed Integer RNG seed (default `2026L`).
#' @return A `data.frame` of class `"glmmplanr_stress_test"`.
#' @examples
#' des <- design_crossed(list(condition = c("a", "b")))
#' dgp <- calibrate_dgp(des, c(a = 0.4, b = 0.6), subj_sd = 0.3, item_sd = 0.4)
#' comparison <- stress_test_models(dgp, "condition", 16, 8, nsim = 3L)
#' head(comparison)
#' @export
stress_test_models <- function(
  dgp,
  effect,
  n_subj,
  n_item,
  nsim = 200L,
  include_null = TRUE,
  alpha = 0.05,
  seed = 2026L
) {
  if (uses_box_muller()) {
    return(isolated_seeded_call("stress_test_models", environment(), match.call()))
  }
  if (!inherits(dgp, "glmmplanr_dgp")) {
    cli::cli_abort("{.arg dgp} must be a {.cls glmmplanr_dgp} object.")
  }
  candidates <- list(
    list(label = "Maximal ML", engine = "cpp_laplace", re = "maximal"),
    list(label = "Uncorrelated Slopes ML", engine = "cpp_laplace", re = "uncorrelated"),
    list(label = "Intercept-Only ML", engine = "cpp_laplace", re = "intercept_only"),
    list(label = "Regularized Laplace", engine = "bayes_laplace", re = "maximal")
  )

  base_st <- allocate_stage_point_streams(seed, stage_slot = 1L, grid_ids = 1L:2L, stride = 2L)
  streams_alt <- derive_rep_substreams(base_st[[1L]], nsim = nsim, offset = 0L)
  streams_null <- derive_rep_substreams(base_st[[2L]], nsim = nsim, offset = 0L)

  rows <- lapply(candidates, function(cand) {
    sp <- specify_target(
      dgp,
      effect = effect,
      criterion = "power",
      alpha = alpha,
      engine = cand$engine,
      analysis_re = cand$re,
      control = target_control(on_singular = "flag_only")
    )
    alt_res <- evaluate_design_point(
      sp, n_subj = n_subj, n_item = n_item, nsim = nsim,
      engine_use = sp$engine, fast_stage = FALSE, stream_seeds = streams_alt
    )
    out_row <- data.frame(
      model = cand$label,
      n_subj = alt_res$n_subj,
      n_item = alt_res$n_item,
      power = alt_res$power,
      conditional_power = alt_res$conditional_power,
      valid_rate = alt_res$valid_rate,
      type1_error = NA_real_,
      ci_width_mean = alt_res$ci_width_mean,
      prob_ci_width_mean_pp = alt_res$prob_ci_width_mean_pp,
      coverage = alt_res$coverage,
      singular_rate = alt_res$singular_rate,
      rmse = alt_res$rmse,
      type_m_error = alt_res$type_m_error,
      stringsAsFactors = FALSE
    )
    if (isTRUE(include_null)) {
      sp_null <- null_target_spec(sp)
      null_res <- evaluate_design_point(
        sp_null, n_subj = n_subj, n_item = n_item, nsim = nsim,
        engine_use = sp_null$engine, fast_stage = FALSE, stream_seeds = streams_null
      )
      out_row$type1_error <- null_res$power
    }
    out_row
  })

  res_df <- do.call(rbind, rows)
  class(res_df) <- c("glmmplanr_stress_test", "data.frame")
  res_df
}
