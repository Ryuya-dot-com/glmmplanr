# Multi-Scenario Robustness and Random-Slope Sensitivity Helpers

dgp_is_equivalent <- function(dgp1, dgp2, tol = 1e-9) {
  if (!inherits(dgp1, "glmmplanr_dgp") || !inherits(dgp2, "glmmplanr_dgp")) return(FALSE)
  # Formula environments can differ between otherwise identical designs.
  fields <- setdiff(names(dgp1$design), c("fixed_rhs", "subj_rhs", "item_rhs", "testlet_rhs",
                                        "items_per_subject", "item_seed"))
  pars <- c("chance", "beta", "mu_cells", "Sigma_subj", "Sigma_item", "Sigma_testlet",
            "L_subj", "L_item", "L_testlet", "n_quad")
  identical(testlet_interaction_sd(dgp1), testlet_interaction_sd(dgp2)) &&
    identical(dgp1$design$items_per_subject, dgp2$design$items_per_subject) &&
    identical(dgp1$design$item_seed, dgp2$design$item_seed) &&
    identical(dropout_rate(dgp1$design), dropout_rate(dgp2$design)) &&
    identical(dgp1$design[fields], dgp2$design[fields]) &&
    isTRUE(all.equal(dgp1[pars], dgp2[pars], tolerance = tol))
}

align_and_bind_rows <- function(df_list) {
  df_list <- Filter(function(x) !is.null(x) && nrow(x) > 0L, df_list)
  if (length(df_list) == 0L) return(data.frame())
  if (length(df_list) == 1L) return(df_list[[1L]])
  all_cols <- Reduce(union, lapply(df_list, names))
  list_cols <- all_cols[vapply(all_cols, function(cn) {
    any(vapply(df_list, function(x) cn %in% names(x) && is.list(x[[cn]]), logical(1L)))
  }, logical(1L))]
  aligned <- lapply(df_list, function(d) {
    miss <- setdiff(all_cols, names(d))
    for (cn in miss) {
      if (cn %in% list_cols) {
        d[[cn]] <- I(vector("list", nrow(d)))
      } else {
        d[[cn]] <- NA
      }
    }
    d[, all_cols, drop = FALSE]
  })
  out <- do.call(rbind, aligned)
  rownames(out) <- NULL
  out
}

format_scenario_counts <- function(vals) {
  paste(
    sprintf("%s = %s", names(vals), format(as.integer(vals), big.mark = ",", trim = TRUE)),
    collapse = "; "
  )
}

rebuild_multi_scenario_surface <- function(grid_df, all_scen_df, scen_specs, exploratory_evals = NULL, prefer = "balanced") {
  scen_names <- names(scen_specs)
  surfaces <- lapply(scen_names, function(snm) {
    sp <- scen_specs[[snm]]
    ver_s <- all_scen_df[all_scen_df$scenario == snm, , drop = FALSE]
    exp_s <- if (!is.null(exploratory_evals) && (snm %in% names(exploratory_evals))) {
      exploratory_evals[[snm]]
    } else {
      NULL
    }
    if (!is.null(exp_s) && nrow(exp_s) > 0L) {
      keep <- !(paste(exp_s$n_subj, exp_s$n_item) %in% paste(ver_s$n_subj, ver_s$n_item))
      eval_s <- align_and_bind_rows(list(exp_s[keep, , drop = FALSE], ver_s))
    } else {
      eval_s <- ver_s
    }
    build_surrogate_surface(grid_df, eval_s, sp, prefer = prefer)
  })
  combined <- robust_surface(grid_df, surfaces)
  cand_keys <- unique(paste(all_scen_df$n_subj, all_scen_df$n_item))
  combined$evaluated <- paste(combined$n_subj, combined$n_item) %in% cand_keys
  combined
}

#' Multi-Scenario Robust Sample-Size Comparison (Pure Sample-Size Increments)
#' @noRd
run_multi_scenario_analysis <- function(
  base_spec,
  scenarios,
  standard_scenario,
  grid_df,
  n_grid_stride,
  prefer,
  method,
  n_eval_init,
  n_eval_adapt,
  nsim_init,
  nsim_final,
  nsim_topup,
  topup_policy = "independent",
  mc_adjustment = "none",
  n_verify,
  early_stop = FALSE,
  core_res = NULL,
  workers,
  parallel_backend,
  seed,
  verbose
) {
  scen_names <- names(scenarios)
  if (is.null(scen_names) || any(!nzchar(scen_names))) {
    scen_names <- paste0("scenario_", seq_along(scenarios))
    names(scenarios) <- scen_names
  }

  if (anyDuplicated(scen_names)) cli::cli_abort("Scenario names must be unique.")
  if (is.null(standard_scenario)) {
    if ("standard" %in% scen_names) {
      standard_scenario <- "standard"
    } else {
      effs <- vapply(scenarios, function(d) {
        abs(true_prob_effect(d, base_spec$contrast_weights, base_spec$prob_estimand))
      }, numeric(1L))
      standard_scenario <- scen_names[order(effs)[ceiling(length(effs) / 2)]]
    }
  } else if (!(standard_scenario %in% scen_names)) {
    cli::cli_abort("{.arg standard_scenario} ({.val {standard_scenario}}) not found in {.val {scen_names}}.")
  }

  make_scen_spec <- function(dgp_s) {
    sp_s <- base_spec
    sp_s$dgp <- dgp_s
    sp_s$effect_col <- resolve_effect_column(dgp_s$design, base_spec$effect_input)
    sp_s$contrast_weights <- align_cell_weights(
      dgp_s$design, base_spec$contrast_weights, arg = "contrast_weights"
    )
    chance_s <- if (base_spec$analysis_link == "match_dgp") dgp_s$chance else 0.0
    sp_s$engine <- resolve_available_engine(base_spec$requested_engine, chance = chance_s, notify = FALSE)
    if (!is.null(base_spec$requested_analysis_re)) {
      sp_s$analysis_re <- resolve_analysis_re(dgp_s, base_spec$requested_analysis_re)
    }
    req_s1_re <- if (!is.null(base_spec$requested_stage1_re)) base_spec$requested_stage1_re else "match"
    sp_s$stage1_re <- if (identical(req_s1_re, "uncorrelated") && identical(sp_s$analysis_re, "maximal")) {
      "uncorrelated"
    } else {
      sp_s$analysis_re
    }
    if (is.null(base_spec$requested_stage1_engine)) {
      sp_s$stage1_engine <- if (sp_s$engine %in% c("stan_laplace", "stan_hmc", "bayes_laplace")) {
        "bayes_laplace"
      } else {
        "cpp_laplace"
      }
    } else {
      sp_s$stage1_engine <- resolve_available_engine(
        base_spec$requested_stage1_engine, chance = chance_s, notify = FALSE
      )
    }
    null_st <- probability_contrast_null_status(dgp_s, sp_s$contrast_weights, sp_s$scale, sp_s$prob_estimand)
    sp_s$has_null_target <- null_st$ok
    if (!null_st$ok && sp_s$criterion %in% c("power", "both", "joint")) {
      cli::cli_abort(null_st$reason)
    }
    sp_s$true_beta <- unname(dgp_s$beta[[sp_s$effect_col]])
    sp_s$true_prob_diff <- true_prob_effect(
      dgp_s, weights = sp_s$contrast_weights, estimand = base_spec$prob_estimand
    )
    sp_s$true_prob_diff_pp <- sp_s$true_prob_diff * 100.0
    sp_s$true_effect_target <- switch(
      base_spec$scale,
      logit = sp_s$true_beta,
      prob = sp_s$true_prob_diff,
      pp = sp_s$true_prob_diff_pp
    )
    sp_s
  }

  scen_specs <- lapply(scenarios, make_scen_spec)

  scen_sched_mat <- vapply(scen_specs, function(sp) {
    compute_scheduled_trials(sp$dgp$design, grid_df$n_subj, grid_df$n_item)
  }, integer(nrow(grid_df)))
  if (is.null(dim(scen_sched_mat))) {
    scen_sched_mat <- matrix(scen_sched_mat, nrow = nrow(grid_df))
  }
  grid_df$n_trials_scheduled_max <- as.integer(apply(scen_sched_mat, 1L, max))
  grid_df$n_trials_scheduled_mean <- rowMeans(scen_sched_mat)
  grid_df$preference_score <- compute_preference_score(grid_df, prefer = prefer)

  full_grid_confirm <- identical(method, "grid") || (nrow(grid_df) <= (n_eval_init + n_eval_adapt))

  if (full_grid_confirm) {
    # When every lattice point is confirmed directly, skip redundant exploratory
    # evaluations and proceed straight to independent per-scenario verification.
    cand_grid <- grid_df
    exploratory_evals <- NULL
  } else {
    # Explore each scenario, then verify designs on the intersection of their
    # predicted feasible regions, including compromises absent from individual fronts.
    searches <- lapply(seq_along(scen_specs), function(i) {
      sp <- scen_specs[[i]]
      if (!is.null(core_res) && dgp_is_equivalent(base_spec$dgp, sp$dgp)) return(core_res)
      if (isTRUE(verbose)) {
        cli::cli_inform("Exploring scenario {i}/{length(scen_specs)} ({.val {scen_names[i]}})...")
      }
      run_single_spec_search(
        sp, grid_df, n_grid_stride, prefer, method, n_eval_init, n_eval_adapt,
        nsim_init, nsim_final, 0L, topup_policy, "none", n_verify, early_stop, workers,
        parallel_backend, seed, stage_offset = 7L + 3L * (i - 1L), verbose = verbose)
    })
    names(searches) <- scen_names
    surfaces <- lapply(searches, `[[`, "surface")
    combined_s1 <- robust_surface(grid_df, surfaces)
    candidates <- lapply(searches, function(res) res$verified[c("n_subj", "n_item")])
    robust_idx <- select_stage2_candidates(combined_s1, base_spec, prefer, n_verify)
    cand_coords <- unique(rbind(do.call(rbind, candidates),
                                grid_df[robust_idx, c("n_subj", "n_item"), drop = FALSE]))
    cand_match <- match(
      paste(cand_coords$n_subj, cand_coords$n_item),
      paste(grid_df$n_subj, grid_df$n_item)
    )
    cand_grid <- grid_df[cand_match, , drop = FALSE]
    exploratory_evals <- lapply(searches, function(res) {
      ev <- res$evaluated
      if (is.null(ev) || !nrow(ev)) return(NULL)
      ev[!vapply(ev, is.list, logical(1L))]
    })
  }

  scen_evals <- vector("list", length(scenarios))
  names(scen_evals) <- scen_names

  seq_two_stage <- isTRUE(nsim_topup > 0L && identical(topup_policy, "independent"))

  verify_scenario <- function(sub_grid, snm, s_i) {
    sp <- scen_specs[[snm]]
    sp$mc_conf_level <- sequential_topup_conf_level(
      sp$mc_conf_level, round_num = 1L,
      n_decisions = if (identical(mc_adjustment, "bonferroni")) nrow(cand_grid) * length(scenarios) else 1L,
      mode = if (seq_two_stage) "two_stage" else "none")
    streams <- function(slot) allocate_stage_point_streams(
      seed, stage_slot = slot, grid_ids = sub_grid$orig_grid_id, stride = n_grid_stride)
    verification_offset <- 7L + 3L * length(scenarios)
    batch <- evaluate_grid_batch(sp, sub_grid, streams(verification_offset + s_i), nsim_final,
      sp$engine, FALSE, workers, parallel_backend, keep_raw = nsim_topup > 0L,
      allow_curtail = early_stop,
      stage_label = sprintf("Scenario %s (%d/%d)", snm, s_i, length(scenarios)),
      verbose = verbose)
    df <- attach_grid_meta(batch$summary, sub_grid)
    df$stage <- "scenario_final"
    if (nsim_topup > 0L) {
      df <- apply_indeterminate_topup(df, batch$raw, sub_grid,
        streams(verification_offset + length(scenarios) + s_i), sp, nsim_topup,
        topup_policy = topup_policy, workers = workers, parallel_backend = parallel_backend,
        verbose = verbose)
    }
    df
  }

  for (s_i in seq_along(scenarios)) {
    snm <- scen_names[s_i]
    # Confirm the selected family on fresh streams, with its multiplicity
    # fixed before decisions or optional top-up. Search outcomes are exploratory.
    ev_s <- verify_scenario(cand_grid, snm, s_i)
    ev_s$scenario <- snm
    rownames(ev_s) <- NULL
    scen_evals[[snm]] <- ev_s
  }

  all_scen_df <- align_and_bind_rows(scen_evals)
  rownames(all_scen_df) <- NULL

  if (identical(mc_adjustment, "bonferroni")) {
    conf_eff <- sequential_topup_conf_level(
      base_spec$mc_conf_level,
      round_num = 1L,
      n_decisions = nrow(all_scen_df),
      mode = if (seq_two_stage) "two_stage" else "none"
    )
    all_scen_df <- reclassify_verified_mc(all_scen_df, base_spec, conf_eff)
  }

  plans <- rank_scenarios(
    all_scen_df, standard_scenario, standard_scenario, base_spec,
    prefer = prefer, scen_order = scen_names
  )
  plans$surface <- rebuild_multi_scenario_surface(
    grid_df, all_scen_df, scen_specs, exploratory_evals, prefer = prefer
  )
  plans$specs <- scen_specs
  plans$exploratory_evaluations <- exploratory_evals
  plans
}

robust_surface <- function(grid, surfaces) {
  combine <- function(column, fun) Reduce(fun, lapply(surfaces, `[[`, column))
  for (nm in c("pred_power", "pred_aipe_prob", "pred_joint_prob", "pred_rope_prob",
               "pred_valid_rate", "obs_valid_rate", "obs_attainment", "prob_feasible")) {
    if (all(vapply(surfaces, function(s) nm %in% names(s), logical(1L)))) {
      grid[[nm]] <- combine(nm, pmin)
    }
  }
  for (nm in c("pred_ci_width", "pred_prob_ci_width_pp")) grid[[nm]] <- combine(nm, pmax)
  for (nm in c("pred_criterion_met", "obs_criterion_met")) {
    if (all(vapply(surfaces, function(s) nm %in% names(s), logical(1L)))) {
      grid[[nm]] <- combine(nm, `&`)
    }
  }
  grid
}

rank_scenarios <- function(all_scen_df, standard_scenario, toughest_name, spec,
                           prefer = NULL, scen_order = NULL) {
  ordered_names <- if (!is.null(scen_order)) {
    intersect(scen_order, unique(as.character(all_scen_df$scenario)))
  } else {
    unique(as.character(all_scen_df$scenario))
  }
  scen_evals_raw <- split(all_scen_df, all_scen_df$scenario)
  eval_order <- c(toughest_name, setdiff(ordered_names, toughest_name))
  scen_evals <- scen_evals_raw[eval_order]
  std_df <- scen_evals[[standard_scenario]]
  keys <- paste(std_df$n_subj, std_df$n_item)
  scen_evals <- lapply(scen_evals, function(d) {
    idx <- match(keys, paste(d$n_subj, d$n_item))
    if (anyNA(idx)) cli::cli_abort("Missing scenario evaluation at a candidate design.")
    d[idx, , drop = FALSE]
  })
  ordered_evals <- scen_evals[ordered_names]
  std_df$n_trials_scheduled_min <- as.integer(std_df$n_trials_scheduled)
  std_df$n_trials_scheduled_max <- as.integer(std_df$n_trials_scheduled)
  std_df$n_trials_scheduled_mean <- as.numeric(std_df$n_trials_scheduled)
  if ("n_trials_expected" %in% names(std_df)) {
    std_df$n_trials_expected_min <- as.integer(std_df$n_trials_expected)
    std_df$n_trials_expected_max <- as.integer(std_df$n_trials_expected)
  }
  standard <- select_recommended_from_verified(std_df, spec, prefer = prefer)
  feasible_std <- standard$feasible
  opt_std <- standard$optimal
  statuses <- c("fail", "indeterminate", "pass")
  limiting <- lapply(seq_along(keys), function(i) {
    rows <- do.call(rbind, lapply(scen_evals, function(d) d[i, , drop = FALSE]))
    ord <- order(match(rows$decision_status, statuses), attainment_score(rows, spec))
    best_r <- rows[ord[1L], , drop = FALSE]
    sched_i <- vapply(ordered_evals, function(d) as.integer(d$n_trials_scheduled[i]), integer(1L))
    best_r$n_trials_scheduled_min <- min(sched_i)
    best_r$n_trials_scheduled_max <- max(sched_i)
    best_r$n_trials_scheduled_mean <- mean(sched_i)
    if ("n_trials_expected" %in% names(best_r)) {
      exp_i <- vapply(ordered_evals, function(d) as.integer(d$n_trials_expected[i]), integer(1L))
      best_r$n_trials_expected_min <- min(exp_i)
      best_r$n_trials_expected_max <- max(exp_i)
    }
    best_r
  })
  verified <- do.call(rbind, limiting)
  verified$criterion_met <- Reduce(`&`, lapply(scen_evals, `[[`, "criterion_met"))
  robust <- select_recommended_from_verified(verified, spec, prefer = prefer)
  feasible_rob <- robust$feasible
  opt_rob <- robust$optimal
  pick_r <- match(paste(opt_rob$n_subj, opt_rob$n_item), keys)
  pick_s <- match(paste(opt_std$n_subj, opt_std$n_item), keys)

  rob_sched_by_scen <- stats::setNames(
    vapply(ordered_evals, function(d) as.integer(d$n_trials_scheduled[pick_r]), integer(1L)),
    ordered_names
  )
  rob_exp_by_scen <- if ("n_trials_expected" %in% names(std_df)) {
    stats::setNames(
      vapply(ordered_evals, function(d) as.integer(d$n_trials_expected[pick_r]), integer(1L)),
      ordered_names
    )
  } else {
    NULL
  }
  schedules_differ <- any(verified$n_trials_scheduled_min != verified$n_trials_scheduled_max)
  expected_differ <- !is.null(rob_exp_by_scen) &&
    any(verified$n_trials_expected_min != verified$n_trials_expected_max)

  delta_subj <- as.integer(opt_rob$n_subj - opt_std$n_subj)
  delta_item <- as.integer(opt_rob$n_item - opt_std$n_item)
  delta_trials <- if (schedules_differ) {
    NA_integer_
  } else {
    as.integer(opt_rob$n_trials_scheduled - opt_std$n_trials_scheduled)
  }

  if (!feasible_std || !feasible_rob) {
    delta_subj <- delta_item <- delta_trials <- NA_integer_
  }

  comparison_df <- data.frame(
    plan = c(
      sprintf("Standard Plan (%s)", standard_scenario),
      sprintf("Robust All-Scenario Plan (%d scenarios)", length(scen_evals))
    ),
    feasible = c(feasible_std, feasible_rob),
    n_subj = c(opt_std$n_subj, opt_rob$n_subj),
    n_item = c(opt_std$n_item, opt_rob$n_item),
    n_trials_scheduled = c(
      opt_std$n_trials_scheduled,
      if (schedules_differ) NA_integer_ else opt_rob$n_trials_scheduled
    ),
    delta_subj = c(0L, delta_subj),
    delta_item = c(0L, delta_item),
    delta_trials = c(0L, delta_trials),
    worst_case_power = c(
      min(vapply(scen_evals, function(d) d$power[pick_s], numeric(1L))),
      min(vapply(scen_evals, function(d) d$power[pick_r], numeric(1L)))
    ),
    worst_case_ci_width = c(
      max(vapply(scen_evals, function(d) d$ci_width_mean[pick_s], numeric(1L))),
      max(vapply(scen_evals, function(d) d$ci_width_mean[pick_r], numeric(1L)))
    ),
    stringsAsFactors = FALSE
  )
  if (schedules_differ) {
    comparison_df$scheduled_trials_by_scenario <- c(
      sprintf("%s = %s", standard_scenario, format(as.integer(opt_std$n_trials_scheduled), big.mark = ",", trim = TRUE)),
      format_scenario_counts(rob_sched_by_scen)
    )
  }

  list(
    standard_scenario = standard_scenario,
    toughest_scenario = as.character(opt_rob$scenario[1L]),
    optimal_standard = opt_std,
    optimal_robust = opt_rob,
    feasible_standard = feasible_std,
    feasible_robust = feasible_rob,
    delta_subj = delta_subj,
    delta_item = delta_item,
    delta_trials = delta_trials,
    schedules_differ = schedules_differ,
    expected_differ = expected_differ,
    robust_trials_by_scenario = rob_sched_by_scen,
    robust_expected_by_scenario = rob_exp_by_scen,
    comparison = comparison_df,
    verified = verified,
    scenario_evaluations = all_scen_df
  )
}

has_random_slopes <- function(dgp) {
  any(dgp$tau_subj[-1L] > 0) || any(dgp$tau_item[-1L] > 0) ||
    (isTRUE(dgp$design$has_testlet) && any(dgp$tau_testlet[-1L] > 0))
}

compute_sensitivity_table <- function(
  spec,
  n_subj,
  n_item,
  orig_grid_id,
  n_grid_stride,
  nsim,
  seed,
  workers,
  parallel_backend,
  verbose = FALSE,
  source_scenario = NA_character_
) {
  mults <- c(optimistic = 0.70, standard = 1.00, conservative = 1.35)
  if (!has_random_slopes(spec$dgp)) return(NULL)
  dgps <- sensitivity_dgps(spec$dgp, slope_multipliers = mults)
  grid_one <- data.frame(n_subj = n_subj, n_item = n_item, grid_id = 1L)

  rows <- lapply(seq_along(dgps), function(idx) {
    nm <- names(dgps)[idx]
    sp_lvl <- spec
    sp_lvl$dgp <- dgps[[idx]]
    sp_lvl$true_beta <- unname(dgps[[idx]]$beta[[spec$effect_col]])
    sp_lvl$true_prob_diff <- true_prob_effect(
      dgps[[idx]], weights = spec$contrast_weights, estimand = spec$prob_estimand
    )
    sp_lvl$true_prob_diff_pp <- sp_lvl$true_prob_diff * 100.0
    sp_lvl$true_effect_target <- switch(
      spec$scale,
      logit = sp_lvl$true_beta,
      prob = sp_lvl$true_prob_diff,
      pp = sp_lvl$true_prob_diff_pp
    )

    lvl_base <- allocate_stage_point_streams(
      seed, stage_slot = 3L + idx, grid_ids = orig_grid_id, stride = n_grid_stride
    )

    ev <- evaluate_grid_batch(
      spec = sp_lvl,
      grid_sub = grid_one,
      point_bases = lvl_base,
      nsim = nsim,
      engine_use = spec$engine,
      fast_stage = FALSE,
      workers = workers,
      parallel_backend = parallel_backend,
      stage_label = sprintf("Sensitivity (%s)", nm),
      verbose = verbose
    )$summary
    ev$scenario <- nm
    ev$slope_multiplier <- unname(mults[idx])
    ev$true_beta <- sp_lvl$true_beta
    ev
  })

  out <- do.call(rbind, rows)
  src_chr <- if (is.null(source_scenario) || !length(source_scenario)) {
    NA_character_
  } else {
    as.character(source_scenario[1L])
  }
  if (!is.na(src_chr) && nzchar(src_chr)) {
    out$source_scenario <- rep(src_chr, nrow(out))
  }
  cols_front <- c("scenario", "slope_multiplier", "n_subj", "n_item", "true_beta")
  out <- out[, c(cols_front, setdiff(names(out), cols_front)), drop = FALSE]
  attr(out, "source_scenario") <- src_chr
  out
}
