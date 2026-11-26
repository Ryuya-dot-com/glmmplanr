#' Optimize 2D Sample Size (Subjects x Items) via Monotone GP Surrogate & Pareto Frontier
#'
#' Searches the balanced `(N_subj, N_item)` lattice using a monotone anisotropic
#' Gaussian process surrogate (`method = "surrogate_gp"`) or exhaustive grid
#' search (`method = "grid"`), identifying the **non-dominated Pareto frontier**
#' of `(N_subj, N_item)` combinations that satisfy the target power/precision
#' criteria while accounting for estimation failures, Morris et al. (2019) Monte
#' Carlo uncertainty (`"pass"` / `"indeterminate"` / `"fail"`), optional
#' pre-specified independent top-up replications, and multi-scenario robustness.
#' For person abilities, item difficulties, or population variance components
#' in a random-item Rasch model, use [plan_rasch_precision()].
#'
#' @param spec A `"glmmplanr_spec"` object produced by [specify_target()].
#' @param n_subj_range Numeric length-2 range `c(min, max)` (or explicit vector)
#'   of candidate subject sample sizes, automatically snapped to multiples of
#'   `spec$dgp$design$k_subj`. Single fixed value (e.g., `n_subj_range = 48`) is
#'   also supported.
#' @param n_item_range Numeric length-2 range `c(min, max)` (or explicit vector)
#'   of candidate item sample sizes, automatically snapped to multiples of
#'   `spec$dgp$design$k_item`. Single fixed value (e.g., `n_item_range = 32`) is
#'   also supported.
#' @param prefer Selection rule for choosing the primary recommended design from
#'   among non-dominated Pareto-optimal `(N_subj, N_item)` pairs:
#'   * `"balanced"` (default): Knee/elbow of the feasible contour minimizing
#'     normalized distance from `(min_subj, min_item)`.
#'   * `"fewer_subj"`: Minimizes participant count `N_subj` first, then `N_item`.
#'   * `"fewer_items"`: Minimizes stimulus/item count `N_item` first, then `N_subj`.
#'   * `"fewer_trials"`: Minimizes scheduled total trials (`N_trials_scheduled`).
#' @param scenarios Optional named list of `"glmmplanr_dgp"` objects (e.g., from
#'   [scenario_dgps()] or [sensitivity_dgps()]) for **robust multi-scenario
#'   planning**. When supplied, identifies both the standard-scenario design and
#'   the all-scenario robust design along with the sample-size increments
#'   (`+delta_subj`, `+delta_item`, `+delta_trials`).
#' @param standard_scenario Optional character name selecting the reference
#'   `"standard"` scenario inside `scenarios`.
#' @param method Search strategy: `"surrogate_gp"` (default) or `"grid"`.
#' @param nsim_final Number of Monte Carlo simulations per candidate Pareto
#'   design in Stage 2 final verification (default `1000L`).
#' @param workers Integer number of parallel workers (default `1L`).
#' @param seed Integer master RNG seed driving strictly disjoint `(stage, point, rep)`
#'   L'Ecuyer-CMRG streams/substreams (default `20261119L`).
#' @param verbose Logical; print concise stage progress updates (default `TRUE`).
#' @param control Control list created by [optimize_control()] holding Stage-1
#'   screening, top-up, multiplicity, early-stopping, sensitivity, and parallel
#'   backend settings (`n_eval_init`, `n_eval_adapt`, `nsim_init`, `nsim_topup`,
#'   `topup_policy`, `mc_adjustment`, `n_verify`, `early_stop`,
#'   `nsim_sensitivity`, `parallel_backend`).
#'
#' @return An S3 object of class `"glmmplanr_result"`. A confirmed recommendation
#'   requires `decision_status == "pass"` in every requested scenario.
#'   If `feasible` is false, `optimal` identifies a diagnostic candidate and
#'   `optimal_high_confidence` is null. `verified_pareto` contains confirmed designs.
#'   `verified` includes the response-pattern diagnostics of [evaluate_design_point()].
#'   When pooling with a saved plan that omits response-pattern counts,
#'   `pattern_nsim` counts only diagnosed replications.
#' @examples
#' des <- design_crossed(list(condition = c("a", "b")))
#' dgp <- calibrate_dgp(des, c(a = 0.4, b = 0.6), subj_sd = 0.3, item_sd = 0.4)
#' spec <- specify_target(dgp, "condition", criterion = "power")
#' # Five replications demonstrate the output. Use more to make a design decision.
#' plan <- optimize_sample_size(
#'   spec, c(16, 18), 8, method = "grid", nsim_final = 5L,
#'   workers = 1L, verbose = FALSE,
#'   control = optimize_control(nsim_init = 5L, n_verify = 2L, nsim_sensitivity = 0L)
#' )
#' plan$verified[c("n_subj", "n_item", "nsim", "decision_status")]
#' @export
optimize_sample_size <- function(
  spec,
  n_subj_range,
  n_item_range,
  prefer = c("balanced", "fewer_subj", "fewer_items", "fewer_trials"),
  scenarios = NULL,
  standard_scenario = NULL,
  method = c("surrogate_gp", "grid"),
  nsim_final = 1000L,
  workers = 1L,
  seed = 20261119L,
  verbose = TRUE,
  control = optimize_control()
) {
  if (!inherits(control, "glmmplanr_optimize_control")) {
    cli::cli_abort("{.arg control} must be created by {.fn optimize_control}.")
  }
  if (uses_box_muller()) {
    return(isolated_seeded_call("optimize_sample_size", environment(), match.call()))
  }
  if (!inherits(spec, "glmmplanr_spec")) {
    cli::cli_abort("{.arg spec} must be a {.cls glmmplanr_spec} object.")
  }
  prefer <- match.arg(prefer)
  method <- match.arg(method)
  n_eval_init <- control$n_eval_init
  n_eval_adapt <- control$n_eval_adapt
  nsim_init <- control$nsim_init
  nsim_topup <- control$nsim_topup
  topup_policy <- control$topup_policy
  mc_adjustment <- control$mc_adjustment
  n_verify <- control$n_verify
  early_stop <- control$early_stop
  nsim_sensitivity <- control$nsim_sensitivity
  parallel_backend <- control$parallel_backend
  workers <- max(1L, as.integer(workers))
  seed <- as.integer(seed)
  if (nsim_topup > 0L && identical(topup_policy, "pooled")) {
    cli::cli_warn(c(
      "Using {.code topup_policy = 'pooled'} combines Stage-2 and adaptive top-up replications.",
      "i" = "Data-adaptive pooling invalidates nominal fixed-sample Wilson coverage; use {.code topup_policy = 'independent'} or confirm the final design with {.fn validate_design_point}."
    ))
  }

  # Start persistent mirai daemon pool once across all stages when applicable
  # (including macOS where fork/mclapply can deadlock after CoreText/GCD calls)
  is_macos <- identical(Sys.info()[["sysname"]], "Darwin")
  if (workers > 1L &&
      (parallel_backend == "mirai" ||
       (parallel_backend == "auto" && (.Platform$OS.type != "unix" || is_macos))) &&
      !isTRUE(mirai::daemons_set())) {
    chk <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
    w_pool <- if (nzchar(chk) && chk != "false") min(workers, 2L) else workers
    started_pool <- start_mirai_pool(w_pool)
    if (started_pool) {
      on.exit(try(mirai::daemons(0L), silent = TRUE), add = TRUE)
    }
  }

  des <- spec$dgp$design
  s_grid <- valid_multiples(n_subj_range, des$k_subj, "n_subj_range")
  i_grid <- valid_multiples(n_item_range, des$k_item, "n_item_range")
  if (!is.null(des$items_per_subject)) {
    i_grid <- i_grid[i_grid >= des$items_per_subject]
    if (!length(i_grid)) cli::cli_abort("No candidate item bank can supply {.arg items_per_subject}.")
  }

  if (!is.null(scenarios)) {
    if (!is.list(scenarios) || length(scenarios) < 2L ||
        !all(vapply(scenarios, inherits, logical(1L), "glmmplanr_dgp"))) {
      cli::cli_abort("{.arg scenarios} must contain at least two calibrated DGPs.")
    }
    for (scenario in scenarios) {
      d <- scenario$design
      if (!setequal(names(spec$contrast_weights), d$cells$cell_name)) {
        cli::cli_abort("All scenarios must define the same cells for the requested contrast.")
      }
      if (spec$scale == "logit" && !identical(des$X_cells, d$X_cells)) {
        cli::cli_abort("A logit-scale comparison requires the same coefficient coding in all scenarios.")
      }
      s_grid <- s_grid[s_grid %% d$k_subj == 0L]
      i_grid <- i_grid[i_grid %% d$k_item == 0L]
      if (!is.null(d$items_per_subject)) i_grid <- i_grid[i_grid >= d$items_per_subject]
    }
    if (!length(s_grid) || !length(i_grid)) cli::cli_abort("No candidate satisfies the allocation requirements in every scenario.")
  }

  grid_df <- expand.grid(n_subj = s_grid, n_item = i_grid, KEEP.OUT.ATTRS = FALSE)
  n_grid_stride <- nrow(grid_df)
  grid_df$orig_grid_id <- seq_len(n_grid_stride)
  grid_df$grid_id <- seq_len(n_grid_stride)
  grid_df$n_trials_scheduled <- compute_scheduled_trials(des, grid_df$n_subj, grid_df$n_item)
  grid_df$n_trials_expected <- compute_expected_trials(des, grid_df$n_subj, grid_df$n_item)
  grid_df$preference_score <- compute_preference_score(grid_df, prefer = prefer)

  t_start <- proc.time()
  if (isTRUE(verbose)) {
    cli::cli_inform(
      "Starting 2D sample size planning ({n_grid_stride} lattice points, prefer = {.val {prefer}}, workers = {workers})..."
    )
  }

  run_base <- is.null(scenarios) ||
    any(vapply(scenarios, function(sc) dgp_is_equivalent(spec$dgp, sc), logical(1L)))

  core_res <- if (run_base) {
    run_single_spec_search(
      spec = spec,
      grid_df = grid_df,
      n_grid_stride = n_grid_stride,
      prefer = prefer,
      method = method,
      n_eval_init = as.integer(n_eval_init),
      n_eval_adapt = as.integer(n_eval_adapt),
      nsim_init = as.integer(nsim_init),
      nsim_final = as.integer(nsim_final),
      nsim_topup = nsim_topup,
      topup_policy = topup_policy,
      mc_adjustment = if (is.null(scenarios)) mc_adjustment else "none",
      n_verify = as.integer(n_verify),
      early_stop = early_stop,
      workers = workers,
      parallel_backend = parallel_backend,
      seed = seed,
      stage_offset = 0L,
      verbose = verbose
    )
  } else {
    NULL
  }

  robust_info <- NULL
  if (!is.null(scenarios)) {
    robust_info <- run_multi_scenario_analysis(
      base_spec = spec,
      scenarios = scenarios,
      standard_scenario = standard_scenario,
      grid_df = grid_df,
      n_grid_stride = n_grid_stride,
      prefer = prefer,
      method = method,
      n_eval_init = as.integer(n_eval_init),
      n_eval_adapt = as.integer(n_eval_adapt),
      nsim_init = as.integer(nsim_init),
      nsim_final = as.integer(nsim_final),
      nsim_topup = nsim_topup,
      topup_policy = topup_policy,
      mc_adjustment = mc_adjustment,
      n_verify = as.integer(n_verify),
      early_stop = early_stop,
      core_res = core_res,
      workers = workers,
      parallel_backend = parallel_backend,
      seed = seed,
      verbose = verbose
    )
    optimal <- robust_info$optimal_robust
    feasible <- robust_info$feasible_robust
    optimal_hc <- if (feasible) optimal
    surf <- robust_info$surface
    robust_info$surface <- NULL
    eval_all <- robust_info$scenario_evaluations
    verified_df <- robust_info$verified
    verified_pareto <- extract_2d_pareto(verified_df[verified_df$decision_status == "pass", , drop = FALSE])
  } else {
    optimal <- core_res$optimal
    optimal_hc <- core_res$optimal_high_confidence
    feasible <- core_res$feasible
    surf <- core_res$surface
    eval_all <- core_res$evaluated
    verified_df <- core_res$verified
    verified_pareto <- extract_2d_pareto(verified_df[verified_df$decision_status == "pass", , drop = FALSE])
  }

  sens_spec <- if (!is.null(robust_info) && !is.null(robust_info$specs) &&
                   nzchar(robust_info$toughest_scenario) &&
                   (robust_info$toughest_scenario %in% names(robust_info$specs))) {
    robust_info$specs[[robust_info$toughest_scenario]]
  } else {
    spec
  }

  sens_df <- if (as.integer(nsim_sensitivity) > 0L && length(sens_spec$dgp$design$factor_names) > 0L) {
    if (verbose) {
      scen_lbl <- if (!is.null(robust_info)) paste0(" (", robust_info$toughest_scenario, ")") else ""
      cli::cli_inform(
        "Computing 3-level random-slope sensitivity analysis at N_subj = {optimal$n_subj}, N_item = {optimal$n_item}{scen_lbl}..."
      )
    }
    opt_orig_id <- grid_df$orig_grid_id[
      grid_df$n_subj == optimal$n_subj & grid_df$n_item == optimal$n_item
    ][1L]
    if (is.na(opt_orig_id)) opt_orig_id <- 1L

    compute_sensitivity_table(
      spec = sens_spec,
      n_subj = optimal$n_subj,
      n_item = optimal$n_item,
      orig_grid_id = opt_orig_id,
      n_grid_stride = n_grid_stride,
      nsim = as.integer(nsim_sensitivity),
      seed = seed,
      workers = workers,
      parallel_backend = parallel_backend,
      verbose = verbose
    )
  } else {
    NULL
  }

  pareto_df <- verified_pareto
  n_decisions <- if (!is.null(robust_info)) nrow(robust_info$scenario_evaluations) else nrow(verified_df)
  seq_two_stage <- isTRUE(nsim_topup > 0L && identical(topup_policy, "independent"))
  n_mult <- if (identical(mc_adjustment, "bonferroni")) n_decisions else 1L
  mc_conf_eff <- sequential_topup_conf_level(
    spec$mc_conf_level,
    round_num = 1L,
    n_decisions = n_mult,
    mode = if (seq_two_stage) "two_stage" else "none"
  )
  used_topup_rows <- if (!is.null(robust_info)) {
    ev_chk <- robust_info$scenario_evaluations
    "nsim_topup_used" %in% names(ev_chk) && any(ev_chk$nsim_topup_used > 0L, na.rm = TRUE)
  } else {
    "nsim_topup_used" %in% names(verified_df) && any(verified_df$nsim_topup_used > 0L, na.rm = TRUE)
  }
  adaptive_pooled <- isTRUE(used_topup_rows && identical(topup_policy, "pooled"))
  elapsed_sec <- unname((proc.time() - t_start)[3L])
  if (isTRUE(verbose)) {
    cli::cli_inform(
      "Optimization finished in {round(elapsed_sec, 2)}s -> {toupper(optimal$decision_status[1L])} at N_subj = {optimal$n_subj}, N_item = {optimal$n_item}."
    )
  }

  structure(
    list(
      spec = spec,
      prefer = prefer,
      optimal = optimal,
      optimal_high_confidence = optimal_hc,
      decision_status = optimal$decision_status[1L],
      feasible = feasible,
      surface = surf,
      evaluated = eval_all,
      verified = verified_df,
      verified_pareto = verified_pareto,
      pareto_frontier = pareto_df,
      sensitivity = sens_df,
      robust = robust_info,
      search_meta = list(
        method = method,
        prefer = prefer,
        n_subj_grid = s_grid,
        n_item_grid = i_grid,
        n_grid_total = nrow(grid_df),
        n_grid_verified = nrow(verified_df),
        n_stage2_decisions = n_decisions,
        simulation_version = 3L,
        n_grid_stride = n_grid_stride,
        n_eval_init = as.integer(n_eval_init),
        n_eval_adapt = as.integer(n_eval_adapt),
        nsim_init = as.integer(nsim_init),
        nsim_final = as.integer(nsim_final),
        nsim_topup = nsim_topup,
        topup_policy = topup_policy,
        topup_rounds = if (isTRUE(used_topup_rows)) 1L else 0L,
        sequential_topup_adjusted = seq_two_stage,
        topup_mc_conf_level = if (isTRUE(used_topup_rows)) mc_conf_eff else NA_real_,
        adaptive_pooled = adaptive_pooled,
        mc_adjustment = mc_adjustment,
        mc_conf_level_effective = mc_conf_eff,
        early_stop = early_stop,
        workers = workers,
        seed = seed,
        elapsed_sec = elapsed_sec
      )
    ),
    class = "glmmplanr_result"
  )
}

#' Switch the Preferred Trade-Off Rule Without Re-Running Simulations
#'
#' Re-ranks already verified Pareto-optimal `(N_subj, N_item)` candidates under
#' a different `prefer` rule (`"balanced"`, `"fewer_subj"`, `"fewer_items"`, or
#' `"fewer_trials"`).
#'
#' @param result A `"glmmplanr_result"` object.
#' @param prefer New preference rule.
#' @param verbose Logical; notify if stored random-slope sensitivity table is
#'   cleared because the selected candidate changed (default `FALSE`).
#' @return Updated `"glmmplanr_result"` object.
#' @examples
#' # Create the five-replication demonstration plan from the optimizer help.
#' example("optimize_sample_size", package = "glmmplanr", local = environment(), echo = FALSE)
#' reordered <- switch_preference(plan, "fewer_items", verbose = FALSE)
#' reordered$optimal[c("n_subj", "n_item", "decision_status")]
#' @export
switch_preference <- function(
  result,
  prefer = c("balanced", "fewer_subj", "fewer_items", "fewer_trials"),
  verbose = FALSE
) {
  if (!inherits(result, "glmmplanr_result")) {
    cli::cli_abort("{.arg result} must be a {.cls glmmplanr_result} object.")
  }
  if (is.null(result$spec$prob_boot) && result$spec$scale %in% c("pp", "prob")) {
    cli::cli_abort("Rerun this probability-scale plan to verify the target contrast under the current interval method.")
  }
  prefer <- match.arg(prefer)

  result$prefer <- prefer
  result$search_meta$prefer <- prefer
  result$surface$preference_score <- compute_preference_score(result$surface, prefer = prefer)

  if (!is.null(result$robust)) {
    rc <- result$robust
    ev <- attach_grid_meta(rc$scenario_evaluations, result$surface)
    ranked <- rank_scenarios(ev, rc$standard_scenario, rc$toughest_scenario, result$spec)
    ranked$specs <- rc$specs
    result$robust <- ranked
    result$verified <- ranked$verified
  } else {
    result$verified <- attach_grid_meta(result$verified, result$surface)
  }
  if (!is.null(result$evaluated) && nrow(result$evaluated) > 0L) {
    result$evaluated <- attach_grid_meta(result$evaluated, result$surface)
  }

  sel <- select_recommended_from_verified(result$verified, result$spec)
  result$optimal <- sel$optimal
  result$optimal_high_confidence <- sel$optimal_high_confidence
  result$decision_status <- sel$optimal$decision_status[1L]
  result$feasible <- sel$feasible
  confirmed <- (result$verified$decision_status == "pass")
  result$verified_pareto <- extract_2d_pareto(result$verified[confirmed, , drop = FALSE])
  result$pareto_frontier <- result$verified_pareto
  if (!is.null(result$sensitivity) &&
      any(result$sensitivity$n_subj != sel$optimal$n_subj | result$sensitivity$n_item != sel$optimal$n_item)) {
    if (isTRUE(verbose)) {
      cli::cli_inform(
        "Cleared stored random-slope sensitivity table because selected design moved to ({sel$optimal$n_subj}, {sel$optimal$n_item})."
      )
    }
    result$sensitivity <- NULL
  }
  result
}

lookup_pref_score <- function(df_sub, grid_df) {
  m <- match(paste(df_sub$n_subj, df_sub$n_item), paste(grid_df$n_subj, grid_df$n_item))
  grid_df$preference_score[m]
}

attach_grid_meta <- function(df_sub, grid_df) {
  df_sub$preference_score <- lookup_pref_score(df_sub, grid_df)
  df_sub
}

sequential_topup_conf_level <- function(
  base_conf,
  round_num = 1L,
  n_decisions = 1L,
  mode = c("two_stage", "geometric", "none")
) {
  mode <- match.arg(mode)
  alpha <- 1.0 - base_conf
  m <- max(1L, as.integer(n_decisions))
  denom <- switch(
    mode,
    none = as.numeric(m),
    two_stage = 2.0 * as.numeric(m),
    geometric = as.numeric(m) * (2.0^(max(1L, as.integer(round_num)) + 1L))
  )
  conf <- 1.0 - alpha / denom
  min(1.0 - 1e-12, max(base_conf, conf))
}

reclassify_verified_mc <- function(df, spec, conf_level) {
  if (is.null(df) || nrow(df) == 0L || !is.finite(conf_level)) return(df)
  for (i in seq_len(nrow(df))) {
    n_tot <- df$nsim[i]
    n_val <- df$n_valid[i]
    denom <- if (spec$failure_policy == "effective") n_tot else n_val
    valid_gate_ok <- if (spec$failure_policy == "effective") TRUE else ((n_val / n_tot) >= spec$min_valid_rate)

    k_sig <- if (is.finite(df$effective_power[i])) round(df$effective_power[i] * n_tot) else NA_integer_
    mc_pow <- wilson_mc_ci(k_sig, denom, conf_level = conf_level)
    df$power_mc_lo[i] <- unname(mc_pow["lo"])
    df$power_mc_hi[i] <- unname(mc_pow["hi"])

    mc_aipe <- c(lo = NA_real_, hi = NA_real_)
    if (!is.null(spec$ci_width) && "effective_aipe_prob" %in% names(df) && is.finite(df$effective_aipe_prob[i])) {
      k_aipe <- round(df$effective_aipe_prob[i] * n_tot)
      mc_aipe <- wilson_mc_ci(k_aipe, denom, conf_level = conf_level)
      df$aipe_mc_lo[i] <- unname(mc_aipe["lo"])
      df$aipe_mc_hi[i] <- unname(mc_aipe["hi"])
    }

    mc_joint <- c(lo = NA_real_, hi = NA_real_)
    if (!is.null(spec$ci_width) && "k_joint" %in% names(df) && is.finite(df$k_joint[i])) {
      mc_joint <- wilson_mc_ci(df$k_joint[i], denom, conf_level = conf_level)
      df$joint_mc_lo[i] <- unname(mc_joint["lo"])
      df$joint_mc_hi[i] <- unname(mc_joint["hi"])
    }

    k_rope <- if ("effective_rope_prob" %in% names(df) && is.finite(df$effective_rope_prob[i])) {
      round(df$effective_rope_prob[i] * n_tot)
    } else {
      NA_integer_
    }
    mc_rope <- wilson_mc_ci(k_rope, denom, conf_level = conf_level)
    if ("rope_mc_lo" %in% names(df)) {
      df$rope_mc_lo[i] <- unname(mc_rope["lo"])
      df$rope_mc_hi[i] <- unname(mc_rope["hi"])
    }

    df$decision_status[i] <- compute_decision_status(
      spec = spec,
      valid_gate_ok = valid_gate_ok,
      mc_pow = mc_pow,
      mc_aipe = mc_aipe,
      mc_rope = mc_rope,
      mc_joint = mc_joint
    )
  }
  df
}

#' Core Single-Specification Search + Pareto Verification + Top-Up Helper
#' @noRd
run_single_spec_search <- function(
  spec,
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
  workers,
  parallel_backend,
  seed,
  stage_offset = 0L,
  verbose = TRUE
) {
  des <- spec$dgp$design
  n_grid <- nrow(grid_df)

  slot_s1 <- stage_offset + 1L
  slot_s2 <- stage_offset + 2L
  slot_topup <- stage_offset + 3L

  point_base_s1 <- allocate_stage_point_streams(
    seed, stage_slot = slot_s1, grid_ids = grid_df$orig_grid_id, stride = n_grid_stride
  )
  point_base_s2 <- allocate_stage_point_streams(
    seed, stage_slot = slot_s2, grid_ids = grid_df$orig_grid_id, stride = n_grid_stride
  )
  point_base_topup <- if (nsim_topup > 0L) {
    allocate_stage_point_streams(
      seed, stage_slot = slot_topup, grid_ids = grid_df$orig_grid_id, stride = n_grid_stride
    )
  } else {
    NULL
  }

  attach_pref <- function(df_in) {
    attach_grid_meta(df_in, grid_df)
  }

  seq_two_stage <- isTRUE(nsim_topup > 0L && identical(topup_policy, "independent"))
  spec_look <- spec
  if (seq_two_stage) {
    spec_look$mc_conf_level <- sequential_topup_conf_level(
      spec$mc_conf_level, round_num = 1L, n_decisions = 1L, mode = "two_stage"
    )
  }

  if (method == "grid" || n_grid <= (n_eval_init + n_eval_adapt)) {
    if (verbose) {
      cli::cli_inform(
        "Evaluating full balanced grid ({n_grid} points, multiples of {des$k_subj} x {des$k_item}) with nsim = {nsim_final}..."
      )
    }
    s2_out <- evaluate_grid_batch(
      spec = spec_look,
      grid_sub = grid_df,
      point_bases = point_base_s2,
      nsim = nsim_final,
      engine_use = spec$engine,
      fast_stage = FALSE,
      workers = workers,
      parallel_backend = parallel_backend,
      keep_raw = (nsim_topup > 0L),
      allow_curtail = early_stop,
      stage_label = "Grid Verification",
      verbose = verbose
    )
    verified_df <- attach_pref(s2_out$summary)
    verified_df$stage <- "grid_final"
    verified_raw <- s2_out$raw

    if (nsim_topup > 0L) {
      verified_df <- apply_indeterminate_topup(
        verified_df = verified_df,
        verified_raw = verified_raw,
        grid_df = grid_df,
        point_base_topup = point_base_topup,
        spec = spec_look,
        nsim_topup = nsim_topup,
        topup_policy = topup_policy,
        workers = workers,
        parallel_backend = parallel_backend,
        verbose = verbose
      )
    }
    if (identical(mc_adjustment, "bonferroni")) {
      conf_eff <- sequential_topup_conf_level(
        spec$mc_conf_level,
        round_num = 1L,
        n_decisions = nrow(verified_df),
        mode = if (seq_two_stage) "two_stage" else "none"
      )
      verified_df <- reclassify_verified_mc(verified_df, spec, conf_eff)
    }
    eval_all <- verified_df
    surf <- build_surrogate_surface(grid_df, eval_all, spec, prefer = prefer)
  } else {
    init_idx <- select_initial_lhs_points(grid_df, n_init = n_eval_init, seed = seed)
    if (verbose) {
      cli::cli_inform(
        "Stage 1a: Evaluating {length(init_idx)}/{n_grid} space-filling LHS points ({nsim_init} sims/pt, engine = {.val {spec$stage1_engine}})..."
      )
    }
    eval_s1 <- evaluate_grid_batch(
      spec = spec,
      grid_sub = grid_df[init_idx, , drop = FALSE],
      point_bases = point_base_s1[init_idx],
      nsim = nsim_init,
      engine_use = spec$stage1_engine,
      fast_stage = TRUE,
      workers = workers,
      parallel_backend = parallel_backend,
      stage_label = "Stage 1a (LHS)",
      verbose = verbose
    )$summary
    eval_s1 <- attach_pref(eval_s1)
    eval_s1$stage <- "surrogate_init"
    evaluated_mask <- rep(FALSE, n_grid)
    evaluated_mask[init_idx] <- TRUE

    n_adapt_rem <- min(n_eval_adapt, sum(!evaluated_mask))
    batch_sz <- max(1L, min(2L, n_adapt_rem))
    if (n_adapt_rem > 0L && verbose) {
      cli::cli_inform(
        "Stage 1b: Running {n_adapt_rem} adaptive GP contour-focused evaluations..."
      )
    }
    adapt_done <- 0L
    adapt_total <- n_adapt_rem
    while (n_adapt_rem > 0L && any(!evaluated_mask)) {
      cur_k <- min(batch_sz, n_adapt_rem, sum(!evaluated_mask))
      surf_tmp <- build_surrogate_surface(grid_df, eval_s1, spec, evaluated_mask, prefer = prefer)
      next_idx <- select_diverse_top_acq(surf_tmp, cur_k, evaluated_mask)
      if (length(next_idx) == 0L) break

      batch_res <- evaluate_grid_batch(
        spec = spec,
        grid_sub = grid_df[next_idx, , drop = FALSE],
        point_bases = point_base_s1[next_idx],
        nsim = nsim_init,
        engine_use = spec$stage1_engine,
        fast_stage = TRUE,
        workers = workers,
        parallel_backend = parallel_backend,
        stage_label = sprintf("Stage 1b (%d/%d pts)", adapt_done + length(next_idx), adapt_total),
        verbose = verbose
      )$summary
      batch_res <- attach_pref(batch_res)
      batch_res$stage <- "surrogate_adapt"
      eval_s1 <- rbind(eval_s1, batch_res)
      evaluated_mask[next_idx] <- TRUE
      adapt_done <- adapt_done + length(next_idx)
      n_adapt_rem <- n_adapt_rem - length(next_idx)
    }

    surf_s1 <- build_surrogate_surface(grid_df, eval_s1, spec, evaluated_mask, prefer = prefer)
    verify_idx <- select_stage2_candidates(surf_s1, spec, prefer = prefer, n_verify = n_verify)

    if (verbose) {
      cli::cli_inform(
        "Stage 2: Verifying {length(verify_idx)} Pareto-frontier candidate design(s) with nsim_final = {nsim_final} (engine = {.val {spec$engine}})..."
      )
    }
    s2_out <- evaluate_grid_batch(
      spec = spec_look,
      grid_sub = grid_df[verify_idx, , drop = FALSE],
      point_bases = point_base_s2[verify_idx],
      nsim = nsim_final,
      engine_use = spec$engine,
      fast_stage = FALSE,
      workers = workers,
      parallel_backend = parallel_backend,
      keep_raw = (nsim_topup > 0L),
      allow_curtail = early_stop,
      stage_label = "Stage 2 Verification",
      verbose = verbose
    )
    verified_df <- attach_pref(s2_out$summary)
    verified_df$stage <- "final_verified"
    verified_raw <- s2_out$raw

    if (!any(verified_df$decision_status == "pass")) {
      # Refit GP surface incorporating Stage-2 evaluations so systematic
      # Stage-1 vs Stage-2 estimator shifts are corrected before step-up selection.
      s1_k0 <- !(paste(eval_s1$n_subj, eval_s1$n_item) %in%
                   paste(verified_df$n_subj, verified_df$n_item))
      eval_refit <- align_and_bind_rows(list(eval_s1[s1_k0, , drop = FALSE], verified_df))
      surf_refit <- build_surrogate_surface(grid_df, eval_refit, spec, prefer = prefer)
      gp_feas_idx <- which(surf_refit$pred_criterion_met &
        !(paste(surf_refit$n_subj, surf_refit$n_item) %in% paste(verified_df$n_subj, verified_df$n_item)))
      if (length(gp_feas_idx) > 0L) {
        best_backup <- gp_feas_idx[order(surf_refit$preference_score[gp_feas_idx])[seq_len(min(2L, length(gp_feas_idx)))]]
        ex_out <- evaluate_grid_batch(
          spec = spec_look,
          grid_sub = grid_df[best_backup, , drop = FALSE],
          point_bases = point_base_s2[best_backup],
          nsim = nsim_final,
          engine_use = spec$engine,
          fast_stage = FALSE,
          workers = workers,
          parallel_backend = parallel_backend,
          keep_raw = (nsim_topup > 0L),
          allow_curtail = early_stop,
          stage_label = "Stage 2 Step-Up",
          verbose = verbose
        )
        ex_df <- attach_pref(ex_out$summary)
        ex_df$stage <- "final_verified_stepup"
        verified_df <- rbind(verified_df, ex_df)
        if (nsim_topup > 0L) verified_raw <- c(verified_raw, ex_out$raw)
      }
    }

    if (!any(verified_df$decision_status == "pass")) {
      max_idx <- which.max(grid_df$n_subj * 1e5 + grid_df$n_item)
      already_max <- any(verified_df$n_subj == grid_df$n_subj[max_idx] &
                           verified_df$n_item == grid_df$n_item[max_idx])
      if (!already_max) {
        cn_out <- evaluate_grid_batch(
          spec = spec_look,
          grid_sub = grid_df[max_idx, , drop = FALSE],
          point_bases = point_base_s2[max_idx],
          nsim = nsim_final,
          engine_use = spec$engine,
          fast_stage = FALSE,
          workers = workers,
          parallel_backend = parallel_backend,
          keep_raw = (nsim_topup > 0L),
          stage_label = "Stage 2 Ceiling",
          verbose = verbose
        )
        cn_df <- attach_pref(cn_out$summary)
        cn_df$stage <- "final_verified_ceiling"
        verified_df <- rbind(verified_df, cn_df)
        if (nsim_topup > 0L) verified_raw <- c(verified_raw, cn_out$raw)
      }
    }

    if (nsim_topup > 0L) {
      verified_df <- apply_indeterminate_topup(
        verified_df = verified_df,
        verified_raw = verified_raw,
        grid_df = grid_df,
        point_base_topup = point_base_topup,
        spec = spec_look,
        nsim_topup = nsim_topup,
        topup_policy = topup_policy,
        workers = workers,
        parallel_backend = parallel_backend,
        verbose = verbose
      )
    }

    if (identical(mc_adjustment, "bonferroni")) {
      conf_eff <- sequential_topup_conf_level(
        spec$mc_conf_level,
        round_num = 1L,
        n_decisions = nrow(verified_df),
        mode = if (seq_two_stage) "two_stage" else "none"
      )
      verified_df <- reclassify_verified_mc(verified_df, spec, conf_eff)
    }

    s1_keep <- !(paste(eval_s1$n_subj, eval_s1$n_item) %in%
                   paste(verified_df$n_subj, verified_df$n_item))
    for (cn in setdiff(names(verified_df), names(eval_s1))) eval_s1[[cn]] <- NA
    for (cn in setdiff(names(eval_s1), names(verified_df))) verified_df[[cn]] <- NA
    eval_s1 <- eval_s1[, names(verified_df), drop = FALSE]

    eval_for_gp <- rbind(eval_s1[s1_keep, , drop = FALSE], verified_df)
    eval_all <- rbind(eval_s1, verified_df)
    surf <- build_surrogate_surface(grid_df, eval_for_gp, spec, prefer = prefer)
  }

  sel <- select_recommended_from_verified(verified_df, spec)

  list(
    optimal = sel$optimal,
    optimal_high_confidence = sel$optimal_high_confidence,
    feasible = sel$feasible,
    surface = surf,
    evaluated = eval_all,
    verified = verified_df
  )
}

select_recommended_from_verified <- function(verified_df, spec) {
  passed <- which(verified_df$decision_status == "pass")
  pending <- which(verified_df$decision_status == "indeterminate")

  if (length(passed) || length(pending)) {
    rows <- if (length(passed)) passed else pending
    best <- rows[order(verified_df$preference_score[rows], -verified_df$power[rows],
                       verified_df$ci_width_mean[rows])[1L]]
  } else {
    att <- attainment_score(verified_df, spec)
    best <- order(-att, verified_df$preference_score)[1L]
  }
  candidate <- verified_df[best, , drop = FALSE]
  rownames(candidate) <- NULL

  list(
    optimal = candidate,
    optimal_high_confidence = if (length(passed)) candidate,
    feasible = length(passed) > 0L
  )
}

#' Execute Pre-Specified Independent Top-Up Verification on Indeterminate Candidates
#' @noRd
apply_indeterminate_topup <- function(
  verified_df,
  verified_raw,
  grid_df,
  point_base_topup,
  spec,
  nsim_topup,
  topup_policy = "independent",
  workers,
  parallel_backend,
  verbose = TRUE
) {
  verified_df$nsim_initial <- verified_df$nsim
  verified_df$nsim_topup_used <- 0L
  verified_df$decision_initial <- verified_df$decision_status
  verified_df$topup_only_power <- NA_real_

  ind_idx <- which(verified_df$decision_status == "indeterminate")
  if (length(ind_idx) == 0L) {
    return(verified_df)
  }

  if (verbose) {
    cli::cli_inform(
      "Stage 3 (Top-Up, policy = {.val {topup_policy}}): Running {nsim_topup} independent replications on {length(ind_idx)} borderline ({.val indeterminate}) candidate(s)..."
    )
  }

  g_match <- match(
    paste(verified_df$n_subj[ind_idx], verified_df$n_item[ind_idx]),
    paste(grid_df$n_subj, grid_df$n_item)
  )
  sub_grid <- grid_df[g_match, , drop = FALSE]
  top_bases <- point_base_topup[g_match]

  top_out <- evaluate_grid_batch(
    spec = spec,
    grid_sub = sub_grid,
    point_bases = top_bases,
    nsim = nsim_topup,
    engine_use = spec$engine,
    fast_stage = FALSE,
    workers = workers,
    parallel_backend = parallel_backend,
    keep_raw = identical(topup_policy, "pooled"),
    stage_label = "Stage 3 Top-Up",
    verbose = verbose
  )

  for (j in seq_along(ind_idx)) {
    r_i <- ind_idx[j]
    if (identical(topup_policy, "independent")) {
      fresh <- attach_grid_meta(top_out$summary[j, , drop = FALSE], grid_df)
      fresh$stage <- paste0(verified_df$stage[r_i], "+independent")
      fresh$nsim_initial <- verified_df$nsim_initial[r_i]
      fresh$nsim_topup_used <- fresh$nsim
      fresh$decision_initial <- verified_df$decision_initial[r_i]
      fresh$topup_only_power <- fresh$power
      verified_df[r_i, names(fresh)] <- fresh
    } else {
      raw_init <- verified_raw[[r_i]]
      raw_top <- top_out$raw[[j]]

      pooled <- pool_replications(raw_init, raw_top, spec, verified_df[r_i, , drop = FALSE],
                                  top_out$summary[j, , drop = FALSE])

      pooled <- attach_grid_meta(pooled, grid_df)
      pooled$stage <- paste0(verified_df$stage[r_i], "+topup")
      pooled$nsim_initial <- raw_init$nsim
      pooled$nsim_topup_used <- raw_top$nsim
      pooled$decision_initial <- verified_df$decision_initial[r_i]
      pooled$topup_only_power <- top_out$summary$power[j]

      verified_df[r_i, names(pooled)] <- pooled
    }
  }

  verified_df
}

pool_replications <- function(initial, extra, spec, row, extra_summary) {
  r <- Map(c, initial, extra)
  r$nsim <- initial$nsim + extra$nsim
  summarize <- function(est, se, pest, pse, sing, conv, status, err) {
    summarize_mc_metrics(est, se, pest, pse, sing, conv, status, err, r$nobs,
                         spec, row$n_subj, row$n_item, r$nsim, spec$engine)
  }
  out <- summarize(r$est, r$se, r$prob_est, r$prob_se, r$sing, r$conv, r$status, r$err)
  if (length(r$red_est)) {
    reduced <- summarize(r$red_est, r$red_se, r$red_prob_est, r$red_prob_se,
                         rep(FALSE, r$nsim), r$red_conv, r$red_status, rep("", r$nsim))
    out$reduced_power <- reduced$power
    out$reduced_ci_width_mean <- reduced$ci_width_mean
    out$reduced_coverage <- reduced$coverage
  }
  out[pattern_columns] <- as.list(pooled_pattern_totals(row, extra_summary))
  out$replications <- I(list(if (out$decision_status == "indeterminate") r))
  out
}

#' Run Additional Independent Top-Up Replications on Borderline Candidates
#'
#' Executes `nsim_extra` fresh, independent Monte Carlo replications (using a
#' new L'Ecuyer-CMRG top-level stage slot) on any `"indeterminate"` candidates in
#' an existing `"glmmplanr_result"` object to narrow simulation uncertainty near
#' the target threshold.
#'
#' Even when each top-up batch is independent of earlier replications, repeatedly
#' re-evaluating only `"indeterminate"` candidates across multiple rounds inflates
#' the overall false-pass probability if unadjusted pointwise intervals are used.
#' By default (`sequential_adjustment = "geometric"`), each independent top-up
#' round \eqn{r = 1, 2, \dots} tightens its Wilson confidence level to
#' \eqn{1 - (1 - \text{mc\_conf\_level}) / 2^{r + 1}} (further divided by the
#' number of verified Stage-2 decisions when `mc_adjustment = "bonferroni"`), so
#' the cumulative added false-pass error across all subsequent rounds remains
#' bounded. For a single unpenalized final confirmation of a chosen design, use
#' [validate_design_point()].
#'
#' @param result A `"glmmplanr_result"` object.
#' @param nsim_extra Positive integer number of independent additional
#'   replications per borderline candidate (default `1000L`).
#' @param topup_policy How top-up replications are used to update borderline
#'   candidates:
#'   * `"independent"` (default): Evaluates each `"indeterminate"` candidate
#'     solely on the fresh `nsim_extra` batch without pooling prior Stage-2
#'     draws.
#'   * `"pooled"`: Pools prior and new replications (emits a warning and flags
#'     `search_meta$adaptive_pooled = TRUE` because data-adaptive pooling
#'     invalidates nominal fixed-sample Wilson coverage).
#' @param sequential_adjustment Multiplicity adjustment applied across repeated
#'   independent top-up rounds on `"indeterminate"` candidates:
#'   * `"geometric"` (default): Applies geometric alpha-spending
#'     (\eqn{\alpha_r = \alpha / 2^{r+1}}) to round \eqn{r} Wilson intervals so
#'     repeated re-evaluation does not accumulate unbounded false-pass risk.
#'   * `"none"`: Uses unadjusted round-wise Wilson intervals (emits a warning
#'     because repeated independent testing of `"indeterminate"` candidates
#'     inflates the overall false-pass rate).
#' @param workers Integer number of parallel workers.
#' @param parallel_backend Parallel backend (`"auto"`, `"mirai"`, `"multicore"`,
#'   or `"sequential"`).
#' @param verbose Logical; print progress message (default `TRUE`).
#' @return Updated `"glmmplanr_result"` object with refined candidate diagnostics.
#' @examples
#' # Create the five-replication demonstration plan from the optimizer help.
#' example("optimize_sample_size", package = "glmmplanr", local = environment(), echo = FALSE)
#' more <- topup_verification(plan, nsim_extra = 5L, verbose = FALSE)
#' more$verified[c("n_subj", "n_item", "nsim", "decision_status")]
#' @export
topup_verification <- function(
  result,
  nsim_extra = 1000L,
  topup_policy = c("independent", "pooled"),
  sequential_adjustment = c("geometric", "none"),
  workers = result$search_meta$workers,
  parallel_backend = "auto",
  verbose = TRUE
) {
  if (uses_box_muller()) {
    return(isolated_seeded_call("topup_verification", environment(), match.call()))
  }
  if (!inherits(result, "glmmplanr_result")) {
    cli::cli_abort("{.arg result} must be a {.cls glmmplanr_result} object.")
  }
  topup_policy <- match.arg(topup_policy)
  sequential_adjustment <- match.arg(sequential_adjustment)
  nsim_extra <- as.integer(nsim_extra)
  if (nsim_extra < 1L) cli::cli_abort("{.arg nsim_extra} must be a positive integer.")

  if (is.null(result$spec$prob_boot) && result$spec$scale %in% c("pp", "prob")) {
    cli::cli_abort("Rerun this probability-scale plan before pooling replications with the current interval method.")
  }
  meta <- result$search_meta
  old_specs <- c(list(result$spec), result$robust$specs)
  if (identical(topup_policy, "pooled") && is.null(meta$simulation_version) &&
      any(vapply(old_specs, function(sp) isTRUE(sp$dgp$design$attrition > 0), logical(1)))) {
    cli::cli_abort("Rerun this plan before pooling results from different trial missingness mechanisms.")
  }
  round_num <- if (is.null(meta$topup_rounds)) 1L else as.integer(meta$topup_rounds) + 1L
  rc <- result$robust
  decisions <- if (is.null(rc)) result$verified else rc$scenario_evaluations
  if (!any(decisions$decision_status == "indeterminate")) return(result)
  if (identical(topup_policy, "pooled") &&
      (is.null(meta$simulation_version) || meta$simulation_version < 3L)) {
    used_old_scale <- vapply(old_specs, function(sp) {
      isTRUE(sp$prior_tau_sd < 0.5) &&
        (identical(sp$engine, "bayes_laplace") || identical(sp$on_singular, "switch_bayes"))
    }, logical(1))
    if (any(used_old_scale)) {
      cli::cli_abort("Rerun this plan before pooling replications. The previous regularized estimator used 0.5 for {.arg prior_tau_sd} values below 0.5.")
    }
  }
  if (identical(topup_policy, "pooled")) {
    cli::cli_warn(c(
      "Pooling prior Stage-2 and adaptive top-up replications ({.code topup_policy = 'pooled'}).",
      "i" = "Data-adaptive pooling invalidates nominal fixed-sample Wilson coverage; use {.code topup_policy = 'independent'} or confirm the final design with {.fn validate_design_point}."
    ))
  } else {
    if (identical(sequential_adjustment, "none")) {
      cli::cli_warn(c(
        "Unadjusted repeated top-up verification ({.code sequential_adjustment = 'none'}) inflates overall false-pass probability across rounds.",
        "i" = "Use {.code sequential_adjustment = 'geometric'} or run a single independent confirmation via {.fn validate_design_point}."
      ))
    }
    if (round_num >= 2L) {
      cli::cli_warn(c(
        "Running sequential top-up round {round_num} on {.val indeterminate} candidates.",
        "i" = "Repeated re-evaluation triggers progressively stricter Wilson bounds; for a final confirmation on a chosen design, run {.fn validate_design_point} once."
      ))
    }
  }
  mc_adj <- if (!is.null(meta$mc_adjustment)) meta$mc_adjustment else "none"
  n_mult <- if (identical(mc_adj, "bonferroni")) nrow(decisions) else 1L
  conf_round <- if (identical(topup_policy, "independent")) {
    sequential_topup_conf_level(
      result$spec$mc_conf_level,
      round_num = round_num,
      n_decisions = n_mult,
      mode = sequential_adjustment
    )
  } else {
    sequential_topup_conf_level(
      result$spec$mc_conf_level,
      round_num = 1L,
      n_decisions = n_mult,
      mode = "none"
    )
  }

  if (is.null(rc)) {
    result$verified <- topup_rows(
      result$verified, result$spec, result, nsim_extra,
      50L + round_num, topup_policy, workers, parallel_backend, verbose,
      conf_level_override = conf_round
    )
  } else {
    if (is.null(rc$specs)) cli::cli_abort("Rerun this saved plan to recover the scenario models for top-up.")
    ev <- rc$scenario_evaluations
    names_s <- names(rc$specs)
    slots <- max(50L, 30L + 2L * length(names_s), 7L + 5L * length(names_s)) + (round_num - 1L) * length(names_s) + seq_along(names_s)
    updated <- lapply(seq_along(names_s), function(i) {
      sn <- names_s[i]
      topup_rows(
        ev[ev$scenario == sn, , drop = FALSE], rc$specs[[sn]], result,
        nsim_extra, slots[i], topup_policy, workers, parallel_backend, verbose,
        conf_level_override = conf_round
      )
    })
    ev_new <- align_and_bind_rows(updated)
    result$robust$scenario_evaluations <- ev_new
  }
  result$search_meta$topup_rounds <- round_num
  result$search_meta$topup_policy <- topup_policy
  result$search_meta$sequential_adjustment <- sequential_adjustment
  result$search_meta$topup_mc_conf_level <- conf_round
  if (identical(topup_policy, "independent") && identical(sequential_adjustment, "geometric")) {
    result$search_meta$sequential_topup_adjusted <- TRUE
  }
  if (identical(topup_policy, "pooled")) {
    result$search_meta$adaptive_pooled <- TRUE
  }
  switch_preference(result, result$prefer)
}

topup_rows <- function(
  ver, sp, result, nsim_extra, slot_extra, topup_policy, workers, parallel_backend, verbose,
  conf_level_override = sp$mc_conf_level
) {
  ind_idx <- which(ver$decision_status == "indeterminate")
  if (!length(ind_idx)) return(ver)
  sp_eval <- sp
  sp_eval$mc_conf_level <- conf_level_override
  stride <- result$search_meta$n_grid_stride

  surf <- result$surface
  g_match <- match(
    paste(ver$n_subj[ind_idx], ver$n_item[ind_idx]),
    paste(surf$n_subj, surf$n_item)
  )
  orig_ids <- surf$orig_grid_id[g_match]
  bases <- allocate_stage_point_streams(
    result$search_meta$seed, stage_slot = slot_extra, grid_ids = orig_ids, stride = stride
  )

  if (verbose) {
    cli::cli_inform(
      "Running {nsim_extra} independent top-up replications (policy = {.val {topup_policy}}, Wilson level = {sprintf('%.2f%%', 100 * conf_level_override)}) on {length(ind_idx)} indeterminate candidate(s)..."
    )
  }

  extra_out <- evaluate_grid_batch(
    spec = sp_eval,
    grid_sub = surf[g_match, , drop = FALSE],
    point_bases = bases,
    nsim = nsim_extra,
    engine_use = sp_eval$engine,
    fast_stage = FALSE,
    workers = workers,
    parallel_backend = parallel_backend,
    keep_raw = identical(topup_policy, "pooled"),
    stage_label = "Top-Up Verification",
    verbose = verbose
  )

  for (j in seq_along(ind_idx)) {
    r_i <- ind_idx[j]
    old_r <- ver[r_i, , drop = FALSE]
    new_r <- extra_out$summary[j, , drop = FALSE]
    if (identical(topup_policy, "independent")) {
      fresh <- attach_grid_meta(new_r, surf)
      if ("stage" %in% names(ver)) {
        fresh$stage <- paste0(old_r$stage, "+independent")
      }
      if ("scenario" %in% names(ver)) {
        fresh$scenario <- old_r$scenario
      }
      ver[r_i, names(fresh)] <- fresh
      next
    }
    if (is.list(old_r$replications) && is.list(old_r$replications[[1L]])) {
      summary <- attach_grid_meta(
        pool_replications(old_r$replications[[1L]], extra_out$raw[[j]], sp_eval, old_r, new_r),
        surf
      )
      ver[r_i, names(summary)] <- summary
      next
    }
    ver[r_i, pattern_columns] <- as.list(pooled_pattern_totals(old_r, new_r))
    # Fallback when a saved result retains only aggregate counts.
    n1 <- old_r$nsim; n2 <- new_r$nsim; ntot <- n1 + n2
    v1 <- old_r$n_valid; v2 <- new_r$n_valid; vtot <- v1 + v2

    k_sig <- round(old_r$effective_power * n1) + round(new_r$effective_power * n2)
    eff_pow <- k_sig / ntot
    cond_pow <- if (vtot > 0L) k_sig / vtot else NA_real_

    ver$nsim[r_i] <- ntot
    ver$n_valid[r_i] <- vtot
    ver$valid_rate[r_i] <- vtot / ntot
    ver$nonconv_rate[r_i] <- 1.0 - (vtot / ntot)
    if ("effective_logit_power" %in% names(ver)) {
      ver$effective_logit_power[r_i] <- (old_r$effective_logit_power * n1 + new_r$effective_logit_power * n2) / ntot
    }
    ver$effective_power[r_i] <- eff_pow
    ver$conditional_power[r_i] <- cond_pow

    denom_mc <- if (sp_eval$failure_policy == "effective") ntot else vtot
    pow_use <- if (sp_eval$failure_policy == "effective") eff_pow else cond_pow
    ver$power[r_i] <- pow_use
    ver$power_mcse[r_i] <- if (denom_mc > 0L) sqrt(pow_use * (1 - pow_use) / denom_mc) else NA_real_
    mc_pow <- wilson_mc_ci(k_sig, denom_mc, conf_level = sp_eval$mc_conf_level)
    ver$power_mc_lo[r_i] <- unname(mc_pow["lo"])
    ver$power_mc_hi[r_i] <- unname(mc_pow["hi"])

    if (!is.null(sp_eval$ci_width)) {
      k_aipe <- round(old_r$effective_aipe_prob * n1) + round(new_r$effective_aipe_prob * n2)
      eff_aipe <- k_aipe / ntot
      cond_aipe <- if (vtot > 0L) k_aipe / vtot else NA_real_
      ver$effective_aipe_prob[r_i] <- eff_aipe
      ver$conditional_aipe_prob[r_i] <- cond_aipe
      ver$aipe_prob[r_i] <- if (sp_eval$failure_policy == "effective") eff_aipe else cond_aipe
      mc_aipe <- wilson_mc_ci(k_aipe, denom_mc, conf_level = sp_eval$mc_conf_level)
      ver$aipe_mc_lo[r_i] <- unname(mc_aipe["lo"])
      ver$aipe_mc_hi[r_i] <- unname(mc_aipe["hi"])
    } else {
      mc_aipe <- c(lo = NA_real_, hi = NA_real_)
    }

    mc_joint <- c(lo = NA_real_, hi = NA_real_)
    if ("k_joint" %in% names(ver) && !is.null(sp_eval$ci_width)) {
      k_joint <- old_r$k_joint + new_r$k_joint
      ver$k_joint[r_i] <- k_joint
      ver$joint_prob[r_i] <- if (denom_mc > 0) k_joint / denom_mc else NA_real_
      mc_joint <- wilson_mc_ci(k_joint, denom_mc, sp_eval$mc_conf_level)
      ver$joint_mc_lo[r_i] <- mc_joint["lo"]
      ver$joint_mc_hi[r_i] <- mc_joint["hi"]
    }
    k_rope <- round(old_r$effective_rope_prob * n1) + round(new_r$effective_rope_prob * n2)
    ver$effective_rope_prob[r_i] <- k_rope / ntot
    ver$conditional_rope_prob[r_i] <- if (vtot > 0L) k_rope / vtot else NA_real_
    ver$rope_excl_prob[r_i] <- if (sp_eval$failure_policy == "effective") (k_rope / ntot) else ver$conditional_rope_prob[r_i]
    mc_rope <- wilson_mc_ci(k_rope, denom_mc, conf_level = sp_eval$mc_conf_level)
    ver$rope_mc_lo[r_i] <- unname(mc_rope["lo"])
    ver$rope_mc_hi[r_i] <- unname(mc_rope["hi"])

    pooled <- function(x, y, a, b) {
      if (a == 0L) return(y)
      if (b == 0L) return(x)
      (x * a + y * b) / (a + b)
    }
    for (nm in c("ci_width_mean", "prob_ci_width_mean_pp", "logit_ci_width_mean", "coverage",
                 "prob_coverage", "prob_diff_est_pp", "bias")) {
      ver[[nm]][r_i] <- pooled(old_r[[nm]], new_r[[nm]], v1, v2)
    }
    ver$rmse[r_i] <- sqrt(pooled(old_r$rmse^2, new_r$rmse^2, v1, v2))
    for (nm in c("singular_rate", "max_iter_rate", "error_rate", "n_trials_obs_mean")) {
      ver[[nm]][r_i] <- pooled(old_r[[nm]], new_r[[nm]], n1, n2)
    }
    # Marginal summaries do not retain the order statistics or subgroup counts.
    for (nm in intersect(c("ci_width_q80", "power_nonsingular", "coverage_nonsingular",
                           "type_m_error", "type_s_error", "reduced_power", "reduced_coverage",
                           "reduced_ci_width_mean"), names(ver))) ver[[nm]][r_i] <- NA_real_

    valid_gate_ok <- if (sp_eval$failure_policy == "effective") TRUE else ((vtot / ntot) >= sp_eval$min_valid_rate)
    ver$criterion_met[r_i] <- valid_gate_ok && check_criterion_met(
      sp_eval, ver$power[r_i], ver$aipe_prob[r_i], ver$rope_excl_prob[r_i], ver$joint_prob[r_i]
    )
    ver$decision_status[r_i] <- compute_decision_status(
      sp_eval, valid_gate_ok, mc_pow, mc_aipe, mc_rope, mc_joint
    )
  }

  ver
}

#' Chunked Reproducible Parallel Evaluator Over Candidate Lattice Points
#' @noRd
evaluate_grid_batch <- function(
  spec,
  grid_sub,
  point_bases,
  nsim,
  engine_use,
  fast_stage,
  workers,
  parallel_backend,
  keep_raw = FALSE,
  allow_curtail = FALSE,
  stage_label = NULL,
  verbose = FALSE
) {
  n_pts <- nrow(grid_sub)
  re_primary <- if (fast_stage && !is.null(spec$stage1_re)) spec$stage1_re else spec$analysis_re
  do_cmp <- (spec$on_singular == "flag_and_compare") && (re_primary == "maximal")

  all_streams <- lapply(seq_len(n_pts), function(p_i) {
    derive_rep_substreams(point_bases[[p_i]], nsim = nsim, offset = 0L)
  })

  build_tasks <- function(pt_indices, rep_range) {
    n_rep_sub <- length(rep_range)
    if (length(pt_indices) == 0L || n_rep_sub == 0L) return(list())
    target_tasks <- max(length(pt_indices), workers * 4L, if (isTRUE(verbose)) 16L else 1L)
    chunks_per_pt <- max(1L, min(n_rep_sub, ceiling(target_tasks / length(pt_indices))))
    local_bounds <- split_indices(n_rep_sub, chunks_per_pt)

    out_tasks <- list()
    tid <- 0L
    for (p_i in pt_indices) {
      st_p <- all_streams[[p_i]]
      for (ch in seq_along(local_bounds)) {
        idx_vec <- rep_range[local_bounds[[ch]]]
        tid <- tid + 1L
        out_tasks[[tid]] <- list(
          pt_idx = p_i,
          rep_idx = idx_vec,
          n_subj = grid_sub$n_subj[p_i],
          n_item = grid_sub$n_item[p_i],
          streams = st_p[idx_vec]
        )
      }
    }
    out_tasks
  }

  worker_fn <- function(tsk) {
    schedule <- generate_design_schedule(
      spec$dgp$design, n_subj = tsk$n_subj, n_item = tsk$n_item
    )
    cpp_cache <- if (engine_use %in% c("cpp_laplace", "bayes_laplace")) {
      prepare_cpp_schedule_cache(spec, schedule)
    } else {
      NULL
    }

    nb <- length(tsk$rep_idx)
    est <- rep(NA_real_, nb)
    se <- rep(NA_real_, nb)
    pest <- rep(NA_real_, nb)
    pse <- rep(NA_real_, nb)
    sing <- rep(NA, nb)
    conv <- rep(FALSE, nb)
    status <- rep("engine_error", nb)
    err <- rep("", nb)
    nobs <- rep(0L, nb)
    patterns <- empty_pattern_totals

    red_est <- rep(NA_real_, nb)
    red_se <- rep(NA_real_, nb)
    red_pest <- rep(NA_real_, nb)
    red_pse <- rep(NA_real_, nb)
    red_conv <- rep(FALSE, nb)
    red_status <- rep("engine_error", nb)

    for (k in seq_len(nb)) {
      ro <- with_rng_stream(
        tsk$streams[[k]],
        execute_single_replication(
          spec = spec,
          schedule = schedule,
          cpp_cache = cpp_cache,
          engine_use = engine_use,
          re_primary = re_primary,
          fast_stage = fast_stage,
          do_compare = do_cmp
        )
      )
      est[k] <- ro$primary$est
      se[k] <- ro$primary$se
      pest[k] <- ro$primary$prob_est
      pse[k] <- ro$primary$prob_se
      sing[k] <- ro$sing_init
      conv[k] <- isTRUE(ro$primary$converged)
      status[k] <- ro$primary$status
      err[k] <- ro$primary$error_reason
      nobs[k] <- ro$nobs
      patterns <- combine_pattern_totals(patterns, ro$patterns)
      if (do_cmp && !is.null(ro$reduced)) {
        red_est[k] <- ro$reduced$est
        red_se[k] <- ro$reduced$se
        red_pest[k] <- ro$reduced$prob_est
        red_pse[k] <- ro$reduced$prob_se
        red_conv[k] <- isTRUE(ro$reduced$converged)
        red_status[k] <- ro$reduced$status
      }
    }
    list(
      pt_idx = tsk$pt_idx,
      rep_idx = tsk$rep_idx,
      est = est, se = se, pest = pest, pse = pse,
      sing = sing, conv = conv, status = status, err = err, nobs = nobs, patterns = patterns,
      red_est = red_est, red_se = red_se, red_pest = red_pest, red_pse = red_pse,
      red_conv = red_conv, red_status = red_status
    )
  }

  est_mat <- matrix(NA_real_, nrow = n_pts, ncol = nsim)
  se_mat <- matrix(NA_real_, nrow = n_pts, ncol = nsim)
  pest_mat <- matrix(NA_real_, nrow = n_pts, ncol = nsim)
  pse_mat <- matrix(NA_real_, nrow = n_pts, ncol = nsim)
  sing_mat <- matrix(NA, nrow = n_pts, ncol = nsim)
  conv_mat <- matrix(FALSE, nrow = n_pts, ncol = nsim)
  status_mat <- matrix("engine_error", nrow = n_pts, ncol = nsim)
  err_mat <- matrix("", nrow = n_pts, ncol = nsim)
  nobs_mat <- matrix(0L, nrow = n_pts, ncol = nsim)
  point_patterns <- rep(list(empty_pattern_totals), n_pts)

  red_est_mat <- if (do_cmp) matrix(NA_real_, nrow = n_pts, ncol = nsim) else NULL
  red_se_mat <- if (do_cmp) matrix(NA_real_, nrow = n_pts, ncol = nsim) else NULL
  red_pest_mat <- if (do_cmp) matrix(NA_real_, nrow = n_pts, ncol = nsim) else NULL
  red_pse_mat <- if (do_cmp) matrix(NA_real_, nrow = n_pts, ncol = nsim) else NULL
  red_conv_mat <- if (do_cmp) matrix(FALSE, nrow = n_pts, ncol = nsim) else NULL
  red_status_mat <- if (do_cmp) matrix("engine_error", nrow = n_pts, ncol = nsim) else NULL

  store_chunks <- function(chunk_res_list) {
    for (cr in chunk_res_list) {
      p_i <- cr$pt_idx
      ix <- cr$rep_idx
      est_mat[p_i, ix] <<- cr$est
      se_mat[p_i, ix] <<- cr$se
      pest_mat[p_i, ix] <<- cr$pest
      pse_mat[p_i, ix] <<- cr$pse
      sing_mat[p_i, ix] <<- cr$sing
      conv_mat[p_i, ix] <<- cr$conv
      status_mat[p_i, ix] <<- cr$status
      err_mat[p_i, ix] <<- cr$err
      nobs_mat[p_i, ix] <<- cr$nobs
      point_patterns[[p_i]] <<- combine_pattern_totals(point_patterns[[p_i]], cr$patterns)
      if (do_cmp) {
        red_est_mat[p_i, ix] <<- cr$red_est
        red_se_mat[p_i, ix] <<- cr$red_se
        red_pest_mat[p_i, ix] <<- cr$red_pest
        red_pse_mat[p_i, ix] <<- cr$red_pse
        red_conv_mat[p_i, ix] <<- cr$red_conv
        red_status_mat[p_i, ix] <<- cr$red_status
      }
    }
  }

  use_curtail <- isTRUE(allow_curtail) && !fast_stage && (nsim >= 200L) && (n_pts >= 2L)
  curtailed_mask <- rep(FALSE, n_pts)
  nsim_used <- rep(nsim, n_pts)

  if (use_curtail) {
    n_interim <- max(100L, min(250L, as.integer(nsim %/% 4L)))
    t_phase1 <- build_tasks(seq_len(n_pts), seq_len(n_interim))
    lbl1 <- if (!is.null(stage_label)) paste0(stage_label, " (Interim)") else "Interim Screening"
    store_chunks(run_parallel_tasks(
      t_phase1, worker_fn, workers, parallel_backend,
      stage_label = lbl1, verbose = verbose
    ))

    spec_999 <- spec
    spec_999$mc_conf_level <- 0.999
    idx_int <- seq_len(n_interim)
    for (p_i in seq_len(n_pts)) {
      r999 <- summarize_mc_metrics(
        est = est_mat[p_i, idx_int],
        se = se_mat[p_i, idx_int],
        prob_est = pest_mat[p_i, idx_int],
        prob_se = pse_mat[p_i, idx_int],
        singular = sing_mat[p_i, idx_int],
        converged = conv_mat[p_i, idx_int],
        status = status_mat[p_i, idx_int],
        error_reasons = err_mat[p_i, idx_int],
        nobs = nobs_mat[p_i, idx_int],
        spec = spec_999,
        n_subj = grid_sub$n_subj[p_i],
        n_item = grid_sub$n_item[p_i],
        nsim = n_interim,
        engine_use = engine_use
      )
      if (identical(r999$decision_status, "fail")) {
        curtailed_mask[p_i] <- TRUE
      }
    }
    # Never curtail every candidate: always run full nsim on at least the ceiling point
    if (all(curtailed_mask)) {
      ceil_i <- which.max(grid_sub$n_subj * 1e5 + grid_sub$n_item)
      curtailed_mask[ceil_i] <- FALSE
    }
    nsim_used[curtailed_mask] <- n_interim

    cont_pts <- which(!curtailed_mask)
    if (length(cont_pts) > 0L && n_interim < nsim) {
      t_phase2 <- build_tasks(cont_pts, seq.int(n_interim + 1L, nsim))
      lbl2 <- if (!is.null(stage_label)) paste0(stage_label, " (Full)") else "Full Verification"
      store_chunks(run_parallel_tasks(
        t_phase2, worker_fn, workers, parallel_backend,
        stage_label = lbl2, verbose = verbose
      ))
    }
  } else {
    t_all <- build_tasks(seq_len(n_pts), seq_len(nsim))
    store_chunks(run_parallel_tasks(
      t_all, worker_fn, workers, parallel_backend,
      stage_label = stage_label, verbose = verbose
    ))
  }

  pt_rows <- vector("list", n_pts)
  raw_list <- if (keep_raw) vector("list", n_pts) else NULL

  for (p_i in seq_len(n_pts)) {
    n_u <- nsim_used[p_i]
    ix_u <- seq_len(n_u)

    est_p <- est_mat[p_i, ix_u]
    se_p <- se_mat[p_i, ix_u]
    pest_p <- pest_mat[p_i, ix_u]
    pse_p <- pse_mat[p_i, ix_u]
    sing_p <- sing_mat[p_i, ix_u]
    conv_p <- conv_mat[p_i, ix_u]
    status_p <- status_mat[p_i, ix_u]
    err_p <- err_mat[p_i, ix_u]
    nobs_p <- nobs_mat[p_i, ix_u]

    row_p <- summarize_mc_metrics(
      est = est_p, se = se_p, prob_est = pest_p, prob_se = pse_p,
      singular = sing_p, converged = conv_p,
      status = status_p, error_reasons = err_p, nobs = nobs_p,
      spec = spec, n_subj = grid_sub$n_subj[p_i], n_item = grid_sub$n_item[p_i],
      nsim = n_u, engine_use = engine_use
    )
    row_p[pattern_columns] <- as.list(point_patterns[[p_i]])
    row_p$curtailed <- curtailed_mask[p_i]

    if (do_cmp) {
      red_est_p <- red_est_mat[p_i, ix_u]
      red_se_p <- red_se_mat[p_i, ix_u]
      red_pest_p <- red_pest_mat[p_i, ix_u]
      red_pse_p <- red_pse_mat[p_i, ix_u]
      red_conv_p <- red_conv_mat[p_i, ix_u]
      red_status_p <- red_status_mat[p_i, ix_u]

      red_sum <- summarize_mc_metrics(
        est = red_est_p, se = red_se_p, prob_est = red_pest_p, prob_se = red_pse_p,
        singular = rep(FALSE, n_u), converged = red_conv_p,
        status = red_status_p, error_reasons = rep("", n_u), nobs = nobs_p,
        spec = spec, n_subj = grid_sub$n_subj[p_i], n_item = grid_sub$n_item[p_i],
        nsim = n_u, engine_use = engine_use
      )
      row_p$reduced_power <- red_sum$power
      row_p$reduced_ci_width_mean <- red_sum$ci_width_mean
      row_p$reduced_coverage <- red_sum$coverage
    } else {
      red_est_p <- NULL; red_se_p <- NULL
      red_pest_p <- NULL; red_pse_p <- NULL
      red_conv_p <- NULL; red_status_p <- NULL
    }
    pending <- !fast_stage && row_p$decision_status == "indeterminate"
    if (keep_raw || pending) {
      record <- list(
        nsim = n_u,
        est = est_p, se = se_p, prob_est = pest_p, prob_se = pse_p,
        sing = sing_p, conv = conv_p, status = status_p, err = err_p, nobs = nobs_p,
        red_est = red_est_p, red_se = red_se_p,
        red_prob_est = red_pest_p, red_prob_se = red_pse_p,
        red_conv = red_conv_p, red_status = red_status_p
      )
      if (keep_raw) raw_list[[p_i]] <- record
    }
    if (!fast_stage) row_p$replications <- I(list(if (pending) record))
    pt_rows[[p_i]] <- row_p
  }

  list(summary = do.call(rbind, pt_rows), raw = raw_list)
}

build_surrogate_surface <- function(grid_df, eval_df, spec, evaluated_mask = NULL, prefer = "balanced") {
  surf <- grid_df
  denom_col <- if (spec$failure_policy == "effective") "nsim" else "n_valid"

  gp_power <- fit_monotone_gp_prob(eval_df, grid_df, prob_col = "power", nsim_col = denom_col)
  gp_width <- fit_monotone_gp_width(eval_df, grid_df, width_col = "ci_width_mean")
  gp_width_pp <- fit_monotone_gp_width(eval_df, grid_df, width_col = "prob_ci_width_mean_pp")

  if (!is.null(spec$ci_width)) {
    gp_aipe <- fit_monotone_gp_prob(eval_df, grid_df, prob_col = "aipe_prob", nsim_col = denom_col)
  } else {
    gp_aipe <- list(
      mean_prob = rep(NA_real_, nrow(grid_df)),
      mean_probit = rep(0, nrow(grid_df)),
      sd_probit = rep(0.5, nrow(grid_df))
    )
  }
  gp_rope <- fit_monotone_gp_prob(eval_df, grid_df, prob_col = "rope_excl_prob", nsim_col = denom_col)

  if (is.null(evaluated_mask)) {
    evaluated_mask <- paste(grid_df$n_subj, grid_df$n_item) %in%
      paste(eval_df$n_subj, eval_df$n_item)
  }

  gp_joint <- if (spec$criterion == "joint") {
    fit_monotone_gp_prob(eval_df, grid_df, prob_col = "joint_prob", nsim_col = denom_col)
  }
  acq_info <- compute_surrogate_acquisition(
    grid_df = grid_df,
    spec = spec,
    gp_power = gp_power,
    gp_aipe = gp_aipe,
    gp_rope = gp_rope, gp_joint = gp_joint,
    evaluated_mask = evaluated_mask,
    prefer = prefer
  )

  surf$pred_power <- gp_power$mean_prob
  surf$pred_power_sd <- gp_power$sd_probit * stats::dnorm(gp_power$mean_probit)
  surf$pred_ci_width <- gp_width$mean_width
  surf$pred_prob_ci_width_pp <- gp_width_pp$mean_width
  surf$pred_joint_prob <- if (is.null(gp_joint)) NA_real_ else gp_joint$mean_prob
  surf$pred_aipe_prob <- gp_aipe$mean_prob
  surf$pred_rope_prob <- gp_rope$mean_prob
  surf$prob_feasible <- acq_info$prob_feasible
  surf$preference_score <- acq_info$preference_score
  surf$acquisition <- acq_info$acq
  surf$evaluated <- evaluated_mask

  surf$pred_criterion_met <- vapply(seq_len(nrow(surf)), function(i) {
    check_criterion_met(
      spec = spec,
      power = surf$pred_power[i],
      aipe_prob = surf$pred_aipe_prob[i],
      rope_excl_prob = surf$pred_rope_prob[i], joint_prob = surf$pred_joint_prob[i]
    )
  }, logical(1L))

  surf
}

select_diverse_top_acq <- function(surf, k, evaluated_mask) {
  avail <- which(!evaluated_mask & is.finite(surf$acquisition))
  if (length(avail) <= k) return(avail)
  u1 <- normalize_unit(surf$n_subj)
  u2 <- normalize_unit(surf$n_item)

  chosen <- integer(0)
  scores <- surf$acquisition
  for (step in seq_len(k)) {
    best <- which.max(scores)
    if (!is.finite(scores[best])) break
    chosen <- c(chosen, best)
    scores[best] <- -Inf
    d2 <- (u1 - u1[best])^2 + (u2 - u2[best])^2
    scores <- scores * (1 - 0.65 * exp(-d2 / 0.02))
  }
  chosen
}

select_stage2_candidates <- function(surf, spec, prefer = "balanced", n_verify = 4L) {
  n_verify <- max(1L, min(n_verify, nrow(surf)))

  # Extract predicted feasible Pareto frontier + high-confidence neighbors
  pf <- extract_pareto_frontier(surf, spec)
  if (nrow(pf) > 0L) {
    pf_idx <- match(paste(pf$n_subj, pf$n_item), paste(surf$n_subj, surf$n_item))
    # 1. Best under requested `prefer`
    best_pref <- pf_idx[order(surf$preference_score[pf_idx], -surf$prob_feasible[pf_idx])]
    # 2. High-confidence neighbor (`prob_feasible >= 0.80`) along `prefer`
    hc <- which(surf$prob_feasible >= 0.80)
    ord_hc <- hc[order(surf$preference_score[hc], -surf$prob_feasible[hc])]
    # 3. Spread points along Pareto frontier (fewer_items end, balanced knee, fewer_subj end)
    spread_idx <- unique(pf_idx[round(seq(1, length(pf_idx), length.out = min(length(pf_idx), n_verify)))])

    pool <- unique(c(
      best_pref[1L],
      if (length(ord_hc) > 0L) ord_hc[1L] else integer(0),
      spread_idx,
      best_pref,
      ord_hc
    ))
    return(pool[seq_len(min(n_verify, length(pool)))])
  }

  plausible <- which(surf$prob_feasible >= 0.35)
  if (length(plausible) > 0L) {
    ord_pl <- plausible[order(surf$preference_score[plausible], -surf$prob_feasible[plausible])]
    return(ord_pl[seq_len(min(n_verify, length(ord_pl)))])
  }

  ord <- order(-surf$prob_feasible, -surf$pred_power, surf$pred_ci_width, surf$preference_score)
  ord[seq_len(n_verify)]
}

attainment_score <- function(df, spec) {
  p_val <- ifelse(is.na(df$power), 0.0, df$power)
  a_val <- ifelse(is.na(df$aipe_prob), 0.0, df$aipe_prob)
  r_val <- ifelse(is.na(df$rope_excl_prob), 0.0, df$rope_excl_prob)
  switch(
    spec$criterion,
    power = p_val,
    ci_width = a_val,
    joint = ifelse(is.na(df$joint_prob), 0.0, df$joint_prob),
    both = pmin(p_val / spec$power_min, a_val / spec$ci_assurance),
    rope = r_val
  )
}

#' Extract Exact 2D Non-Dominated Sample-Size Pareto Frontier `(N_subj, N_item)`
#'
#' Keeps feasible `(N_subj, N_item)` pairs where neither `N_subj` nor `N_item`
#' can be reduced without increasing the other while still meeting the target.
#'
#' @noRd
extract_2d_pareto <- function(feas_df) {
  if (is.null(feas_df) || nrow(feas_df) == 0L) {
    return(feas_df)
  }
  ord <- order(feas_df$n_item, feas_df$n_subj)
  d_ord <- feas_df[ord, , drop = FALSE]

  keep <- logical(nrow(d_ord))
  min_subj_so_far <- Inf
  for (i in seq_len(nrow(d_ord))) {
    if (d_ord$n_subj[i] < min_subj_so_far) {
      keep[i] <- TRUE
      min_subj_so_far <- d_ord$n_subj[i]
    }
  }
  out <- d_ord[keep, , drop = FALSE]
  rownames(out) <- NULL
  out
}

extract_pareto_frontier <- function(surf, spec) {
  feas <- surf[surf$pred_criterion_met, , drop = FALSE]
  if (nrow(feas) > 0L) {
    return(extract_2d_pareto(feas))
  }
  max_s <- max(surf$n_subj)
  max_i <- max(surf$n_item)
  out <- surf[surf$n_subj == max_s | surf$n_item == max_i, , drop = FALSE]
  out <- out[order(out$n_item, out$n_subj), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' @rdname print.glmmplanr_design
#' @export
print.glmmplanr_result <- function(x, ...) {
  cli::cli_h2("glmmplanr 2D Sample Size Optimization Result")
  cli::cli_text("{analysis_description(x$spec)}")
  opt <- x$optimal
  unit_lbl <- switch(x$spec$scale, logit = "logit", pp = "pp", prob = "prob")

  status_txt <- switch(
    opt$decision_status[1L],
    pass = "PASS (Target attainment verified: MC lower bound >= target)",
    indeterminate = "INDETERMINATE / BORDERLINE (MC confidence interval straddles target; consider top-up)",
    fail = if (isTRUE(x$feasible)) {
      "FEASIBLE POINT ESTIMATE (Check MC interval)"
    } else {
      "No confirmed design among evaluated candidates"
    }
  )

  role <- if (isTRUE(x$feasible)) "Recommended design" else "Unconfirmed candidate"
  cov_show <- if (x$spec$scale == "logit") opt$coverage else opt$prob_coverage
  has_null <- can_define_null_target(x$spec)
  val_info <- assess_inferential_validity(cov_show, opt$n_valid, x$spec$alpha, has_null = has_null)
  one_cell <- isTRUE(x$spec$dgp$design$n_cells == 1L)
  acc_line <- if (one_cell) {
    sprintf(
      "Target accuracy (%s): true = %.1f%% | Mean CI width = {.strong %.2f pp} (%.3f logit)",
      x$spec$prob_estimand, opt$prob_diff_true_pp, opt$prob_ci_width_mean_pp, opt$logit_ci_width_mean
    )
  } else {
    sprintf(
      "Accuracy diff (%s): true = %+.1f pp | Mean CI width = {.strong %.2f pp} (%.3f logit)",
      x$spec$prob_estimand, opt$prob_diff_true_pp, opt$prob_ci_width_mean_pp, opt$logit_ci_width_mean
    )
  }

  bullets <- c(
    "Target attainment: {.strong {status_txt}}",
    "Inferential validity: {val_info$detail}",
    "{role} ({.code prefer = '{x$prefer}'}): N_subj = {.strong {opt$n_subj}}, N_item = {.strong {opt$n_item}} (Scheduled: {opt$n_trials_scheduled} trials | Expected post-attrition: {opt$n_trials_expected})",
    "Power ({x$spec$failure_policy}) = {.strong {sprintf('%.1f%%', opt$power * 100)}} [{sprintf('%.1f%%', opt$power_mc_lo * 100)}, {sprintf('%.1f%%', opt$power_mc_hi * 100)}] | Conditional Power = {ifelse(is.na(opt$conditional_power), 'NA', sprintf('%.1f%%', opt$conditional_power * 100))} (Valid fits: {sprintf('%.1f%%', opt$valid_rate * 100)})",
    acc_line
  )
  exp_regimes <- detect_exploratory_regimes(x$spec)
  if (length(exp_regimes) > 0L) {
    exp_str <- paste(exp_regimes, collapse = "; ")
    bullets <- c(
      bullets,
      "Validation scope: {.strong EXPLORATORY REGIME} ({exp_str}). Coverage and null Type I error are not pre-benchmarked across schedules; confirm via {.fn validate_design_point}."
    )
  }
  bnd_info <- assess_boundary_sensitivity(x)
  if (!is.null(bnd_info)) {
    buf_str <- if (!is.null(bnd_info$buffer)) {
      sprintf(
        " Higher-margin verified PASS: (%d subj x %d items, margin %+.1f pp).",
        bnd_info$buffer$n_subj, bnd_info$buffer$n_item, 100 * bnd_info$buffer$mc_margin
      )
    } else {
      ""
    }
    bullets <- c(
      bullets,
      sprintf(
        "Boundary sensitivity: Selected PASS is near the threshold (MC lower-bound margin = %+.1f pp); minimum passing (N_subj, N_item) can shift across search seeds.%s Report candidate range or confirm via {.fn validate_design_point}.",
        100 * bnd_info$opt_margin, buf_str
      )
    )
  }
  if (!is.null(opt$joint_prob) && !is.na(opt$joint_prob)) {
    cli::cli_inform(sprintf("Detection and width attainment in the same replication = %.1f%% [%.1f%%, %.1f%%]",
      100 * opt$joint_prob, 100 * opt$joint_mc_lo, 100 * opt$joint_mc_hi))
  }
  if (!is.na(opt$aipe_prob)) {
    bullets <- c(
      bullets,
      sprintf(
        "AIPE Assurance P(CI width <= %.2f %s) = %.1f%% [%.1f%%, %.1f%%]",
        x$spec$ci_width, unit_lbl, opt$aipe_prob * 100, opt$aipe_mc_lo * 100, opt$aipe_mc_hi * 100
      )
    )
  }
  if (x$spec$criterion == "rope" && !is.na(opt$rope_excl_prob)) {
    bullets <- c(
      bullets,
      sprintf(
        "ROPE Exclusion P(CI excludes [%.2f, %.2f] %s) = %.1f%% [%.1f%%, %.1f%%]",
        x$spec$rope[1L], x$spec$rope[2L], unit_lbl,
        opt$rope_excl_prob * 100, opt$rope_mc_lo * 100, opt$rope_mc_hi * 100
      )
    )
  }
  if (!is.null(x$pareto_frontier) && nrow(x$pareto_frontier) > 1L && isTRUE(x$feasible)) {
    pf_str <- paste(
      sprintf("(%d subj x %d items)", x$pareto_frontier$n_subj, x$pareto_frontier$n_item),
      collapse = ", "
    )
    bullets <- c(bullets, "Pareto trade-off frontier: {pf_str}")
  }
  sing_str <- if (is.na(opt$singular_rate)) {
    "not applicable (regularized)"
  } else {
    sprintf("%.1f%%", opt$singular_rate * 100)
  }
  bullets <- c(
    bullets,
    "Diagnostics: Singular = {sing_str} | Non-conv/Fail = {sprintf('%.1f%%', opt$nonconv_rate * 100)} (MaxIter = {sprintf('%.1f%%', opt$max_iter_rate * 100)}) | Stage-2 H1 Coverage = {ifelse(is.na(cov_show), 'NA', sprintf('%.1f%%', cov_show * 100))}"
  )
  if (nzchar(opt$top_error_reason)) {
    bullets <- c(bullets, "Top non-conv/error reason: {.val {opt$top_error_reason}}")
  }
  if (!is.null(x$optimal_high_confidence) && opt$decision_status == "indeterminate") {
    hc <- x$optimal_high_confidence
    bullets <- c(
      bullets,
      "High-confidence PASS alternative: N_subj = {hc$n_subj}, N_item = {hc$n_item} (Power = {sprintf('%.1f%%', hc$power * 100)} [{sprintf('%.1f%%', hc$power_mc_lo * 100)}, {sprintf('%.1f%%', hc$power_mc_hi * 100)}])"
    )
  }
  if (!is.null(x$robust)) {
    rc <- x$robust
    if (isTRUE(rc$feasible_robust)) {
      bullets <- c(
        bullets,
        "All-scenario robustness increment (vs {.val {rc$standard_scenario}}): {sprintf('%+d', rc$delta_subj)} subjects, {sprintf('%+d', rc$delta_item)} items ({sprintf('%+d', rc$delta_trials)} scheduled trials)"
      )
    }
  }
  if (isTRUE(x$search_meta$adaptive_pooled)) {
    bullets <- c(
      bullets,
      "Note: Adaptive top-up pooling was used; nominal fixed-sample Wilson coverage is not guaranteed. Confirm the selected design with {.fn validate_design_point}."
    )
  }
  if (isTRUE(x$search_meta$topup_rounds > 0L) &&
      identical(x$search_meta$topup_policy, "independent") &&
      !is.null(x$search_meta$topup_mc_conf_level) &&
      is.finite(x$search_meta$topup_mc_conf_level)) {
    seq_mode <- if (identical(x$search_meta$sequential_adjustment, "none")) {
      "unadjusted"
    } else if (identical(x$search_meta$sequential_adjustment, "geometric")) {
      "geometric alpha-spending"
    } else {
      "two-stage split"
    }
    bullets <- c(
      bullets,
      sprintf(
        "Sequential top-up: %d independent round(s) (%s, latest Wilson level = %.2f%%)",
        as.integer(x$search_meta$topup_rounds),
        seq_mode,
        100 * x$search_meta$topup_mc_conf_level
      )
    )
  }
  if (identical(x$search_meta$mc_adjustment, "bonferroni") &&
      !is.null(x$search_meta$mc_conf_level_effective)) {
    bullets <- c(
      bullets,
      sprintf(
        "MC multiplicity adjustment: Bonferroni across %d Stage-2 decision(s) (per-decision Wilson level = %.2f%%)",
        as.integer(x$search_meta$n_stage2_decisions),
        100 * x$search_meta$mc_conf_level_effective
      )
    )
  }
  bullets <- c(
    bullets,
    "Search time: {round(x$search_meta$elapsed_sec, 2)} sec ({nrow(x$evaluated)} point evaluations, method = {.val {x$search_meta$method}})"
  )
  cli::cli_ul(bullets)
  invisible(x)
}

#' @rdname print.glmmplanr_design
#' @export
summary.glmmplanr_result <- function(object, ...) {
  print.glmmplanr_result(object, ...)
  selected <- object$optimal
  if (isTRUE(selected$pattern_nsim > 0L)) {
    scenario <- if (is.null(selected$scenario)) "" else paste0(", ", selected$scenario)
    cli::cli_h3("Response patterns at the selected design{scenario} ({selected$pattern_nsim} replications, {selected$pattern_n_invalid} invalid target fits)")
    maxima <- rep(NA_integer_, length(unit_pattern_names))
    maxima[c(1L, 2L, 4L, 5L)] <- unlist(selected[observed_max_columns], use.names = FALSE)
    print(data.frame(pattern = unit_pattern_names,
      replications = unlist(selected[paste0("n_rep_", unit_pattern_names)], use.names = FALSE),
      invalid_fits = unlist(selected[paste0("n_invalid_", unit_pattern_names)], use.names = FALSE),
      max_observed = maxima), row.names = FALSE)
  }
  if (!is.null(object$verified_pareto) && nrow(object$verified_pareto) > 0L) {
    cli::cli_h3("Verified Non-Dominated Pareto Frontier (N_subj x N_item)")
    cols_vp <- c(
      "n_subj", "n_item", "n_trials_scheduled",
      "effective_power", "prob_ci_width_mean_pp", "decision_status"
    )
    vp <- object$verified_pareto[, cols_vp, drop = FALSE]
    vp$effective_power <- round(vp$effective_power, 3)
    vp$prob_ci_width_mean_pp <- round(vp$prob_ci_width_mean_pp, 2)
    print(vp, row.names = FALSE)
  }
  if (!is.null(object$robust)) {
    cli::cli_h3("Multi-Scenario Robust vs Standard Plan Comparison")
    print(object$robust$comparison, row.names = FALSE)
  }
  if (!is.null(object$sensitivity)) {
    cli::cli_h3("Random-Slope Sensitivity Analysis at Selected (N_subj, N_item)")
    s_tbl <- object$sensitivity[, c(
      "scenario", "slope_multiplier", "effective_power", "conditional_power",
      "valid_rate", "prob_ci_width_mean_pp", "logit_ci_width_mean",
      "singular_rate", "coverage", "decision_status"
    )]
    num_cols <- c(
      "effective_power", "conditional_power", "valid_rate",
      "prob_ci_width_mean_pp", "logit_ci_width_mean", "singular_rate", "coverage"
    )
    s_tbl[num_cols] <- lapply(s_tbl[num_cols], round, 3)
    print(s_tbl, row.names = FALSE)
  }
  invisible(object)
}
