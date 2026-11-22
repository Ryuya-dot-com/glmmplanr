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
#'   * `"min_cost"`: Minimizes total study cost computed from `cost` (selected
#'     automatically when `cost` is supplied and `prefer` is omitted).
#' @param cost Optional explicit study cost specification: either a named
#'   numeric vector/list `c(subj = ..., item = ..., trial = ..., fixed = 0)`
#'   giving unit recruitment, item-development, per-trial, and fixed overhead
#'   costs, or a custom `function(n_subj, n_item, n_trials_scheduled)`.
#' @param budget Optional positive numeric maximum total cost. When supplied
#'   together with `cost`, confirmed recommendations must satisfy
#'   `total_cost <= budget` in addition to `decision_status == "pass"`.
#' @param scenarios Optional named list of `"glmmplanr_dgp"` objects (e.g., from
#'   [scenario_dgps()] or [sensitivity_dgps()]) for **robust multi-scenario
#'   planning**. When supplied, identifies both the standard-scenario design and
#'   the all-scenario robust design along with the sample-size increments
#'   (`+delta_subj`, `+delta_item`, `+delta_trials`).
#' @param standard_scenario Optional character name selecting the reference
#'   `"standard"` scenario inside `scenarios`.
#' @param method Search strategy: `"surrogate_gp"` (default) or `"grid"`.
#' @param n_eval_init Number of space-filling initial lattice points evaluated
#'   in Stage 1 when `method = "surrogate_gp"` (default `16L`).
#' @param n_eval_adapt Number of adaptive contour-focused GP evaluations added
#'   after `n_eval_init` in Stage 1 (default `8L`).
#' @param nsim_init Number of Monte Carlo simulations per exploratory point in
#'   Stage 1 (default `80L`).
#' @param nsim_final Number of Monte Carlo simulations per candidate Pareto
#'   design in Stage 2 final verification (default `1000L`).
#' @param nsim_topup Pre-specified number of **independent** additional Monte
#'   Carlo replications automatically run on Stage-2 candidates whose Monte
#'   Carlo confidence interval straddles the target (`decision_status == "indeterminate"`),
#'   default `0L`.
#' @param n_verify Number of Pareto-frontier candidate designs verified at
#'   `nsim_final` replications in Stage 2 (default `4L`).
#' @param early_stop Logical; when `TRUE` (default) and `nsim_final >= 200L`,
#'   applies group-sequential futility curtailment during Stage-2 and
#'   multi-scenario verification: after an interim batch of
#'   `max(100L, min(250L, nsim_final %/% 4L))` replications, any candidate whose
#'   `99.9%` Wilson Monte Carlo confidence interval upper bound already lies
#'   strictly below the target can stop early. Pointwise Monte Carlo intervals
#'   do not provide simultaneous or sequential coverage guarantees.
#' @param nsim_sensitivity Number of Monte Carlo simulations per level used to
#'   build the 3-level random-slope sensitivity table (`optimistic`, `standard`,
#'   `conservative`) at the recommended `(N_subj, N_item)` (default `200L`; set
#'   to `0L` to skip).
#' @param workers Integer number of parallel workers (default `1L`).
#' @param parallel_backend Parallel execution engine: `"auto"` (default),
#'   `"mirai"`, `"multicore"`, or `"sequential"`.
#' @param seed Integer master RNG seed driving strictly disjoint `(stage, point, rep)`
#'   L'Ecuyer-CMRG streams/substreams (default `20261119L`).
#' @param verbose Logical; print concise stage progress updates (default `TRUE`).
#'
#' @return An S3 object of class `"glmmplanr_result"`. A confirmed recommendation
#'   requires `decision_status == "pass"` in every requested scenario.
#'   If `feasible` is false, `optimal` identifies a diagnostic candidate and
#'   `optimal_high_confidence` is null. `verified_pareto` contains confirmed designs.
#'   `verified` includes the response-pattern diagnostics of [evaluate_design_point()].
#'   After a top-up of an older plan, `pattern_nsim` counts only diagnosed replications.
#' @export
optimize_sample_size <- function(
  spec,
  n_subj_range,
  n_item_range,
  prefer = c("balanced", "fewer_subj", "fewer_items", "fewer_trials", "min_cost"),
  cost = NULL,
  budget = NULL,
  scenarios = NULL,
  standard_scenario = NULL,
  method = c("surrogate_gp", "grid"),
  n_eval_init = 16L,
  n_eval_adapt = 8L,
  nsim_init = 80L,
  nsim_final = 1000L,
  nsim_topup = 0L,
  n_verify = 4L,
  early_stop = TRUE,
  nsim_sensitivity = 200L,
  workers = 1L,
  parallel_backend = c("auto", "mirai", "multicore", "sequential"),
  seed = 20261119L,
  verbose = TRUE
) {
  if (!inherits(spec, "glmmplanr_spec")) {
    cli::cli_abort("{.arg spec} must be a {.cls glmmplanr_spec} object.")
  }
  prefer_missing <- missing(prefer)
  prefer <- match.arg(prefer)
  if (!is.null(cost) && prefer_missing) {
    prefer <- "min_cost"
  }
  if (prefer == "min_cost" && is.null(cost)) {
    cli::cli_abort("Supply {.arg cost} when using {.code prefer = 'min_cost'}.")
  }
  if (!is.null(budget)) {
    if (is.null(cost)) {
      cli::cli_abort("Supply {.arg cost} when specifying {.arg budget}.")
    }
    if (!is.numeric(budget) || length(budget) != 1L || !is.finite(budget) || budget <= 0) {
      cli::cli_abort("{.arg budget} must be a single positive finite number.")
    }
    budget <- as.numeric(budget)
  }
  method <- match.arg(method)
  parallel_backend <- match.arg(parallel_backend)
  workers <- max(1L, as.integer(workers))
  seed <- as.integer(seed)
  nsim_topup <- max(0L, as.integer(nsim_topup))
  early_stop <- isTRUE(early_stop)

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
  if (!is.null(cost)) {
    grid_df$total_cost <- evaluate_design_cost(grid_df, cost = cost)
    grid_df$within_budget <- if (!is.null(budget)) (grid_df$total_cost <= budget + 1e-9) else rep(TRUE, n_grid_stride)
    if (!any(grid_df$within_budget)) {
      cli::cli_abort(
        "No candidate in the search range fits within {.arg budget} = {signif(budget, 6)} (minimum candidate cost = {signif(min(grid_df$total_cost), 6)})."
      )
    }
  } else {
    grid_df$within_budget <- rep(TRUE, n_grid_stride)
  }
  grid_df$preference_score <- compute_preference_score(grid_df, prefer = prefer)

  t_start <- proc.time()
  if (isTRUE(verbose)) {
    cost_note <- if (!is.null(cost)) {
      if (!is.null(budget)) {
        sprintf(" | cost range [%.2f, %.2f], budget <= %.2f", min(grid_df$total_cost), max(grid_df$total_cost), budget)
      } else {
        sprintf(" | cost range [%.2f, %.2f]", min(grid_df$total_cost), max(grid_df$total_cost))
      }
    } else {
      ""
    }
    cli::cli_inform(
      "Starting 2D sample size planning ({n_grid_stride} lattice points, prefer = {.val {prefer}}, workers = {workers}{cost_note})..."
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
  min_cost_over_budget <- NULL
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
    min_cost_over_budget <- robust_info$min_cost_over_budget
    surf <- robust_info$surface
    robust_info$surface <- NULL
    eval_all <- robust_info$scenario_evaluations
    verified_df <- robust_info$verified
    in_b_ver <- if ("within_budget" %in% names(verified_df)) verified_df$within_budget else rep(TRUE, nrow(verified_df))
    verified_pareto <- extract_2d_pareto(verified_df[verified_df$decision_status == "pass" & in_b_ver, , drop = FALSE])
  } else {
    optimal <- core_res$optimal
    optimal_hc <- core_res$optimal_high_confidence
    feasible <- core_res$feasible
    min_cost_over_budget <- core_res$min_cost_over_budget
    surf <- core_res$surface
    eval_all <- core_res$evaluated
    verified_df <- core_res$verified
    in_b_ver <- if ("within_budget" %in% names(verified_df)) verified_df$within_budget else rep(TRUE, nrow(verified_df))
    verified_pareto <- extract_2d_pareto(verified_df[verified_df$decision_status == "pass" & in_b_ver, , drop = FALSE])
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
      cost = cost,
      budget = budget,
      optimal = optimal,
      optimal_high_confidence = optimal_hc,
      min_cost_over_budget = min_cost_over_budget,
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
        cost = cost,
        budget = budget,
        n_subj_grid = s_grid,
        n_item_grid = i_grid,
        n_grid_total = nrow(grid_df),
        n_grid_verified = nrow(verified_df),
        simulation_version = 2L,
        n_grid_stride = n_grid_stride,
        n_eval_init = as.integer(n_eval_init),
        n_eval_adapt = as.integer(n_eval_adapt),
        nsim_init = as.integer(nsim_init),
        nsim_final = as.integer(nsim_final),
        nsim_topup = nsim_topup,
        early_stop = early_stop,
        workers = workers,
        seed = seed,
        elapsed_sec = elapsed_sec
      )
    ),
    class = "glmmplanr_result"
  )
}

#' Switch the Preferred Trade-Off or Cost/Budget Rule Without Re-Running Simulations
#'
#' Re-ranks already verified Pareto-optimal `(N_subj, N_item)` candidates under
#' a different `prefer` rule (`"balanced"`, `"fewer_subj"`, `"fewer_items"`,
#' `"fewer_trials"`, or `"min_cost"`) or updated `cost` / `budget` assumptions.
#'
#' @param result A `"glmmplanr_result"` object.
#' @param prefer New preference rule.
#' @param cost Optional updated cost specification (defaults to `result$cost`).
#' @param budget Optional updated maximum total budget (defaults to `result$budget`;
#'   pass `Inf` to remove an existing budget cap).
#' @param verbose Logical; notify if stored random-slope sensitivity table is
#'   cleared because the selected candidate changed (default `FALSE`).
#' @return Updated `"glmmplanr_result"` object.
#' @export
switch_preference <- function(
  result,
  prefer = c("balanced", "fewer_subj", "fewer_items", "fewer_trials", "min_cost"),
  cost = result$cost,
  budget = result$budget,
  verbose = FALSE
) {
  if (!inherits(result, "glmmplanr_result")) {
    cli::cli_abort("{.arg result} must be a {.cls glmmplanr_result} object.")
  }
  if (is.null(result$spec$prob_boot) && result$spec$scale %in% c("pp", "prob")) {
    cli::cli_abort("Rerun this probability-scale plan to verify the target contrast under the current interval method.")
  }
  prefer_missing <- missing(prefer)
  prefer <- match.arg(prefer)
  if (!missing(cost) && !is.null(cost) && prefer_missing) {
    prefer <- "min_cost"
  }
  if (prefer == "min_cost" && is.null(cost)) {
    cli::cli_abort("Supply {.arg cost} when using {.code prefer = 'min_cost'}.")
  }
  if (!is.null(budget) && is.infinite(budget)) budget <- NULL
  if (!is.null(budget)) {
    if (is.null(cost)) cli::cli_abort("Supply {.arg cost} when specifying {.arg budget}.")
    if (!is.numeric(budget) || length(budget) != 1L || !is.finite(budget) || budget <= 0) {
      cli::cli_abort("{.arg budget} must be a single positive finite number.")
    }
  }

  result$prefer <- prefer
  result$cost <- cost
  result$budget <- budget
  result$search_meta$prefer <- prefer
  result$search_meta$cost <- cost
  result$search_meta$budget <- budget

  if (!is.null(cost)) {
    result$surface$total_cost <- evaluate_design_cost(result$surface, cost = cost)
    result$surface$within_budget <- if (!is.null(budget)) {
      result$surface$total_cost <= budget + 1e-9
    } else {
      rep(TRUE, nrow(result$surface))
    }
  } else {
    result$surface$total_cost <- NULL
    result$surface$within_budget <- rep(TRUE, nrow(result$surface))
  }
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
  result$min_cost_over_budget <- sel$min_cost_over_budget
  result$decision_status <- sel$optimal$decision_status[1L]
  result$feasible <- sel$feasible
  in_b <- if ("within_budget" %in% names(result$verified)) result$verified$within_budget else rep(TRUE, nrow(result$verified))
  confirmed <- (result$verified$decision_status == "pass") & in_b
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
  m <- match(paste(df_sub$n_subj, df_sub$n_item), paste(grid_df$n_subj, grid_df$n_item))
  df_sub$preference_score <- grid_df$preference_score[m]
  if ("total_cost" %in% names(grid_df)) {
    df_sub$total_cost <- grid_df$total_cost[m]
  } else {
    df_sub$total_cost <- NULL
  }
  if ("within_budget" %in% names(grid_df)) {
    df_sub$within_budget <- grid_df$within_budget[m]
  }
  df_sub
}

#' Core Single-Specification Search + Pareto Verification + Top-Up Helper
#' @keywords internal
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

  in_b_grid <- if ("within_budget" %in% names(grid_df)) grid_df$within_budget else rep(TRUE, n_grid)

  if (method == "grid" || n_grid <= (n_eval_init + n_eval_adapt)) {
    if (verbose) {
      cli::cli_inform(
        "Evaluating full balanced grid ({n_grid} points, multiples of {des$k_subj} x {des$k_item}) with nsim = {nsim_final}..."
      )
    }
    s2_out <- evaluate_grid_batch(
      spec = spec,
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
        spec = spec,
        nsim_topup = nsim_topup,
        workers = workers,
        parallel_backend = parallel_backend,
        verbose = verbose
      )
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
      spec = spec,
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

    in_b_ver <- if ("within_budget" %in% names(verified_df)) verified_df$within_budget else rep(TRUE, nrow(verified_df))
    if (!any(verified_df$decision_status == "pass" & in_b_ver)) {
      gp_feas_idx <- which(surf_s1$pred_criterion_met & in_b_grid &
        !(paste(surf_s1$n_subj, surf_s1$n_item) %in% paste(verified_df$n_subj, verified_df$n_item)))
      if (length(gp_feas_idx) > 0L) {
        best_backup <- gp_feas_idx[order(surf_s1$preference_score[gp_feas_idx])[seq_len(min(2L, length(gp_feas_idx)))]]
        ex_out <- evaluate_grid_batch(
          spec = spec,
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

    in_b_ver <- if ("within_budget" %in% names(verified_df)) verified_df$within_budget else rep(TRUE, nrow(verified_df))
    if (!any(verified_df$decision_status == "pass" & in_b_ver)) {
      elig <- which(in_b_grid)
      if (!length(elig)) elig <- seq_len(nrow(grid_df))
      max_idx <- elig[which.max(grid_df$n_subj[elig] * 1e5 + grid_df$n_item[elig])]
      already_max <- any(verified_df$n_subj == grid_df$n_subj[max_idx] &
                           verified_df$n_item == grid_df$n_item[max_idx])
      if (!already_max) {
        cn_out <- evaluate_grid_batch(
          spec = spec,
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
        spec = spec,
        nsim_topup = nsim_topup,
        workers = workers,
        parallel_backend = parallel_backend,
        verbose = verbose
      )
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
    min_cost_over_budget = sel$min_cost_over_budget,
    feasible = sel$feasible,
    surface = surf,
    evaluated = eval_all,
    verified = verified_df
  )
}

select_recommended_from_verified <- function(verified_df, spec) {
  in_b <- if ("within_budget" %in% names(verified_df)) {
    !is.na(verified_df$within_budget) & verified_df$within_budget
  } else {
    rep(TRUE, nrow(verified_df))
  }
  passed <- which(verified_df$decision_status == "pass" & in_b)
  pending <- which(verified_df$decision_status == "indeterminate" & in_b)

  if (length(passed) || length(pending)) {
    rows <- if (length(passed)) passed else pending
    best <- rows[order(verified_df$preference_score[rows], -verified_df$power[rows],
                       verified_df$ci_width_mean[rows])[1L]]
  } else {
    pool <- if (any(in_b)) which(in_b) else seq_len(nrow(verified_df))
    att <- attainment_score(verified_df[pool, , drop = FALSE], spec)
    best <- pool[order(-att, verified_df$preference_score[pool])[1L]]
  }
  candidate <- verified_df[best, , drop = FALSE]
  rownames(candidate) <- NULL

  over_pass <- which(verified_df$decision_status == "pass" & !in_b)
  min_over <- if (!length(passed) && length(over_pass) && ("total_cost" %in% names(verified_df))) {
    ov <- verified_df[over_pass[order(verified_df$total_cost[over_pass], verified_df$preference_score[over_pass])[1L]], , drop = FALSE]
    rownames(ov) <- NULL
    ov
  } else {
    NULL
  }

  list(
    optimal = candidate,
    optimal_high_confidence = if (length(passed)) candidate,
    min_cost_over_budget = min_over,
    feasible = length(passed) > 0L
  )
}

#' Execute Pre-Specified Independent Top-Up Verification on Indeterminate Candidates
#' @keywords internal
apply_indeterminate_topup <- function(
  verified_df,
  verified_raw,
  grid_df,
  point_base_topup,
  spec,
  nsim_topup,
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
      "Stage 3 (Top-Up): Running {nsim_topup} independent replications on {length(ind_idx)} borderline ({.val indeterminate}) candidate(s)..."
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
    keep_raw = TRUE,
    stage_label = "Stage 3 Top-Up",
    verbose = verbose
  )

  for (j in seq_along(ind_idx)) {
    r_i <- ind_idx[j]
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
  aligned <- lapply(df_list, function(d) {
    miss <- setdiff(all_cols, names(d))
    for (cn in miss) d[[cn]] <- NA
    d[, all_cols, drop = FALSE]
  })
  out <- do.call(rbind, aligned)
  rownames(out) <- NULL
  out
}

#' Multi-Scenario Robust Sample-Size Comparison (Pure Sample-Size Increments)
#' @keywords internal
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
        abs(true_prob_effect(d, base_spec$contrast_weights[d$design$cells$cell_name], base_spec$prob_estimand))
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
    sp_s$contrast_weights <- base_spec$contrast_weights[dgp_s$design$cells$cell_name]
    sp_s$engine <- resolve_available_engine(base_spec$requested_engine,
      chance = if (base_spec$analysis_link == "match_dgp") dgp_s$chance else 0, notify = FALSE)
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
      nsim_init, nsim_final, nsim_topup, n_verify, early_stop, workers,
      parallel_backend, seed, stage_offset = 7L + 3L * (i - 1L), verbose = verbose)
  })
  names(searches) <- scen_names
  surfaces <- lapply(searches, `[[`, "surface")
  combined <- robust_surface(grid_df, surfaces)
  candidates <- lapply(searches, function(res) res$verified[c("n_subj", "n_item")])
  robust_idx <- select_stage2_candidates(combined, base_spec, prefer, n_verify)
  cand_coords <- unique(rbind(do.call(rbind, candidates),
                              grid_df[robust_idx, c("n_subj", "n_item"), drop = FALSE]))
  cand_match <- match(
    paste(cand_coords$n_subj, cand_coords$n_item),
    paste(grid_df$n_subj, grid_df$n_item)
  )
  cand_grid <- grid_df[cand_match, , drop = FALSE]
  cand_keys <- paste(cand_grid$n_subj, cand_grid$n_item)

  scen_evals <- vector("list", length(scenarios))
  names(scen_evals) <- scen_names

  verify_scenario <- function(sub_grid, snm, s_i) {
    sp <- scen_specs[[snm]]
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
        streams(verification_offset + length(scenarios) + s_i), sp, nsim_topup, workers, parallel_backend, verbose)
    }
    df
  }

  for (s_i in seq_along(scenarios)) {
    snm <- scen_names[s_i]
    pre_df <- searches[[snm]]$verified
    pre_keys <- paste(pre_df$n_subj, pre_df$n_item)
    hit_idx <- match(cand_keys, pre_keys)
    miss_mask <- is.na(hit_idx)
    reused_part <- if (any(!miss_mask)) {
      pre_df[hit_idx[!miss_mask], , drop = FALSE]
    } else {
      NULL
    }
    new_part <- if (any(miss_mask)) {
      verify_scenario(cand_grid[miss_mask, , drop = FALSE], snm, s_i)
    } else {
      NULL
    }
    comb <- align_and_bind_rows(list(reused_part, new_part))
    ord_back <- match(cand_keys, paste(comb$n_subj, comb$n_item))
    ev_s <- attach_grid_meta(comb[ord_back, , drop = FALSE], cand_grid)
    ev_s$scenario <- snm
    rownames(ev_s) <- NULL
    scen_evals[[snm]] <- ev_s
  }

  all_scen_df <- align_and_bind_rows(scen_evals)
  rownames(all_scen_df) <- NULL

  plans <- rank_scenarios(all_scen_df, standard_scenario, standard_scenario, base_spec)
  combined$evaluated <- paste(combined$n_subj, combined$n_item) %in% cand_keys
  plans$surface <- combined
  plans$specs <- scen_specs
  plans
}

robust_surface <- function(grid, surfaces) {
  combine <- function(column, fun) Reduce(fun, lapply(surfaces, `[[`, column))
  for (nm in c("pred_power", "pred_aipe_prob", "pred_joint_prob", "pred_rope_prob", "prob_feasible")) {
    grid[[nm]] <- combine(nm, pmin)
  }
  for (nm in c("pred_ci_width", "pred_prob_ci_width_pp")) grid[[nm]] <- combine(nm, pmax)
  grid$pred_criterion_met <- combine("pred_criterion_met", `&`)
  grid
}

rank_scenarios <- function(all_scen_df, standard_scenario, toughest_name, spec) {
  scen_evals <- split(all_scen_df, all_scen_df$scenario)
  scen_evals <- scen_evals[c(toughest_name, setdiff(names(scen_evals), toughest_name))]
  std_df <- scen_evals[[standard_scenario]]
  keys <- paste(std_df$n_subj, std_df$n_item)
  scen_evals <- lapply(scen_evals, function(d) {
    idx <- match(keys, paste(d$n_subj, d$n_item))
    if (anyNA(idx)) cli::cli_abort("Missing scenario evaluation at a candidate design.")
    d[idx, , drop = FALSE]
  })
  standard <- select_recommended_from_verified(std_df, spec)
  feasible_std <- standard$feasible
  opt_std <- standard$optimal
  statuses <- c("fail", "indeterminate", "pass")
  limiting <- lapply(seq_along(keys), function(i) {
    rows <- do.call(rbind, lapply(scen_evals, function(d) d[i, , drop = FALSE]))
    ord <- order(match(rows$decision_status, statuses), attainment_score(rows, spec))
    rows[ord[1L], , drop = FALSE]
  })
  verified <- do.call(rbind, limiting)
  verified$criterion_met <- Reduce(`&`, lapply(scen_evals, `[[`, "criterion_met"))
  robust <- select_recommended_from_verified(verified, spec)
  feasible_rob <- robust$feasible
  opt_rob <- robust$optimal
  pick_r <- match(paste(opt_rob$n_subj, opt_rob$n_item), keys)

  delta_subj <- as.integer(opt_rob$n_subj - opt_std$n_subj)
  delta_item <- as.integer(opt_rob$n_item - opt_std$n_item)
  delta_trials <- as.integer(opt_rob$n_trials_scheduled - opt_std$n_trials_scheduled)
  has_cost <- ("total_cost" %in% names(opt_rob)) && is.finite(opt_rob$total_cost) && is.finite(opt_std$total_cost)
  delta_cost <- if (has_cost) as.numeric(opt_rob$total_cost - opt_std$total_cost) else NA_real_

  if (!feasible_std || !feasible_rob) {
    delta_subj <- delta_item <- delta_trials <- NA_integer_
    delta_cost <- NA_real_
  }

  comparison_df <- data.frame(
    plan = c(
      sprintf("Standard Plan (%s)", standard_scenario),
      sprintf("Robust All-Scenario Plan (%d scenarios)", length(scen_evals))
    ),
    feasible = c(feasible_std, feasible_rob),
    n_subj = c(opt_std$n_subj, opt_rob$n_subj),
    n_item = c(opt_std$n_item, opt_rob$n_item),
    n_trials_scheduled = c(opt_std$n_trials_scheduled, opt_rob$n_trials_scheduled),
    delta_subj = c(0L, delta_subj),
    delta_item = c(0L, delta_item),
    delta_trials = c(0L, delta_trials),
    worst_case_power = c(
      min(vapply(scen_evals, function(d) d$power[d$n_subj == opt_std$n_subj & d$n_item == opt_std$n_item][1L], numeric(1L))),
      min(vapply(scen_evals, function(d) d$power[pick_r], numeric(1L)))
    ),
    worst_case_ci_width = c(
      max(vapply(scen_evals, function(d) d$ci_width_mean[d$n_subj == opt_std$n_subj & d$n_item == opt_std$n_item][1L], numeric(1L))),
      max(vapply(scen_evals, function(d) d$ci_width_mean[pick_r], numeric(1L)))
    ),
    stringsAsFactors = FALSE
  )
  if (has_cost) {
    comparison_df$total_cost <- c(opt_std$total_cost, opt_rob$total_cost)
    comparison_df$delta_cost <- c(0.0, delta_cost)
  }

  list(
    standard_scenario = standard_scenario,
    toughest_scenario = as.character(opt_rob$scenario[1L]),
    optimal_standard = opt_std,
    optimal_robust = opt_rob,
    min_cost_over_budget = robust$min_cost_over_budget,
    feasible_standard = feasible_std,
    feasible_robust = feasible_rob,
    delta_subj = delta_subj,
    delta_item = delta_item,
    delta_trials = delta_trials,
    delta_cost = delta_cost,
    comparison = comparison_df,
    verified = verified,
    scenario_evaluations = all_scen_df
  )
}

#' Run Additional Independent Top-Up Replications on Borderline Candidates
#'
#' Executes `nsim_extra` fresh, independent Monte Carlo replications (using a
#' new L'Ecuyer-CMRG top-level stage slot) on any `"indeterminate"` candidates in
#' an existing `"glmmplanr_result"` object to narrow simulation uncertainty near
#' the target threshold.
#'
#' @param result A `"glmmplanr_result"` object.
#' @param nsim_extra Positive integer number of independent additional
#'   replications per borderline candidate (default `1000L`).
#' @param workers Integer number of parallel workers.
#' @param parallel_backend Parallel backend (`"auto"`, `"mirai"`, `"multicore"`,
#'   or `"sequential"`).
#' @param verbose Logical; print progress message (default `TRUE`).
#' @return Updated `"glmmplanr_result"` object with refined candidate diagnostics.
#' @export
topup_verification <- function(
  result,
  nsim_extra = 1000L,
  workers = result$search_meta$workers,
  parallel_backend = "auto",
  verbose = TRUE
) {
  if (!inherits(result, "glmmplanr_result")) {
    cli::cli_abort("{.arg result} must be a {.cls glmmplanr_result} object.")
  }
  nsim_extra <- as.integer(nsim_extra)
  if (nsim_extra < 1L) cli::cli_abort("{.arg nsim_extra} must be a positive integer.")

  if (is.null(result$spec$prob_boot) && result$spec$scale %in% c("pp", "prob")) {
    cli::cli_abort("Rerun this probability-scale plan before pooling replications with the current interval method.")
  }
  meta <- result$search_meta
  old_specs <- c(list(result$spec), result$robust$specs)
  if (is.null(meta$simulation_version) &&
      any(vapply(old_specs, function(sp) isTRUE(sp$dgp$design$attrition > 0), logical(1)))) {
    cli::cli_abort("Rerun this plan before pooling results from different trial missingness mechanisms.")
  }
  round_num <- if (is.null(meta$topup_rounds)) 1L else meta$topup_rounds + 1L
  rc <- result$robust
  decisions <- if (is.null(rc)) result$verified else rc$scenario_evaluations
  if (!any(decisions$decision_status == "indeterminate")) return(result)
  if (is.null(rc)) {
    result$verified <- topup_rows(result$verified, result$spec, result, nsim_extra,
                                 50L + round_num, workers, parallel_backend, verbose)
  } else {
    if (is.null(rc$specs)) cli::cli_abort("Rerun this legacy plan to recover the scenario models for top-up.")
    ev <- rc$scenario_evaluations
    names_s <- names(rc$specs)
    slots <- max(50L, 30L + 2L * length(names_s), 7L + 5L * length(names_s)) + (round_num - 1L) * length(names_s) + seq_along(names_s)
    updated <- lapply(seq_along(names_s), function(i) {
      sn <- names_s[i]
      topup_rows(ev[ev$scenario == sn, , drop = FALSE], rc$specs[[sn]], result,
                 nsim_extra, slots[i], workers, parallel_backend, verbose)
    })
    result$robust$scenario_evaluations <- align_and_bind_rows(updated)
  }
  result$search_meta$topup_rounds <- round_num
  switch_preference(result, result$prefer)
}

topup_rows <- function(ver, sp, result, nsim_extra, slot_extra, workers, parallel_backend, verbose) {
  ind_idx <- which(ver$decision_status == "indeterminate")
  if (!length(ind_idx)) return(ver)
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
      "Running {nsim_extra} independent top-up replications on {length(ind_idx)} indeterminate candidate(s)..."
    )
  }

  extra_out <- evaluate_grid_batch(
    spec = sp,
    grid_sub = surf[g_match, , drop = FALSE],
    point_bases = bases,
    nsim = nsim_extra,
    engine_use = sp$engine,
    fast_stage = FALSE,
    workers = workers,
    parallel_backend = parallel_backend,
    keep_raw = TRUE,
    stage_label = "Top-Up Verification",
    verbose = verbose
  )

  for (j in seq_along(ind_idx)) {
    r_i <- ind_idx[j]
    old_r <- ver[r_i, , drop = FALSE]
    new_r <- extra_out$summary[j, , drop = FALSE]
    if (is.list(old_r$replications) && is.list(old_r$replications[[1L]])) {
      summary <- attach_grid_meta(
        pool_replications(old_r$replications[[1L]], extra_out$raw[[j]], sp, old_r, new_r),
        surf
      )
      ver[r_i, names(summary)] <- summary
      next
    }
    ver[r_i, pattern_columns] <- as.list(pooled_pattern_totals(old_r, new_r))
    # Older logit-scale results only retain aggregate counts.
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

    denom_mc <- if (sp$failure_policy == "effective") ntot else vtot
    pow_use <- if (sp$failure_policy == "effective") eff_pow else cond_pow
    ver$power[r_i] <- pow_use
    ver$power_mcse[r_i] <- if (denom_mc > 0L) sqrt(pow_use * (1 - pow_use) / denom_mc) else NA_real_
    mc_pow <- wilson_mc_ci(k_sig, denom_mc, conf_level = sp$mc_conf_level)
    ver$power_mc_lo[r_i] <- unname(mc_pow["lo"])
    ver$power_mc_hi[r_i] <- unname(mc_pow["hi"])

    if (!is.null(sp$ci_width)) {
      k_aipe <- round(old_r$effective_aipe_prob * n1) + round(new_r$effective_aipe_prob * n2)
      eff_aipe <- k_aipe / ntot
      cond_aipe <- if (vtot > 0L) k_aipe / vtot else NA_real_
      ver$effective_aipe_prob[r_i] <- eff_aipe
      ver$conditional_aipe_prob[r_i] <- cond_aipe
      ver$aipe_prob[r_i] <- if (sp$failure_policy == "effective") eff_aipe else cond_aipe
      mc_aipe <- wilson_mc_ci(k_aipe, denom_mc, conf_level = sp$mc_conf_level)
      ver$aipe_mc_lo[r_i] <- unname(mc_aipe["lo"])
      ver$aipe_mc_hi[r_i] <- unname(mc_aipe["hi"])
    } else {
      mc_aipe <- c(lo = NA_real_, hi = NA_real_)
    }

    mc_joint <- c(lo = NA_real_, hi = NA_real_)
    if ("k_joint" %in% names(ver) && !is.null(sp$ci_width)) {
      k_joint <- old_r$k_joint + new_r$k_joint
      ver$k_joint[r_i] <- k_joint
      ver$joint_prob[r_i] <- if (denom_mc > 0) k_joint / denom_mc else NA_real_
      mc_joint <- wilson_mc_ci(k_joint, denom_mc, sp$mc_conf_level)
      ver$joint_mc_lo[r_i] <- mc_joint["lo"]
      ver$joint_mc_hi[r_i] <- mc_joint["hi"]
    }
    k_rope <- round(old_r$effective_rope_prob * n1) + round(new_r$effective_rope_prob * n2)
    ver$effective_rope_prob[r_i] <- k_rope / ntot
    ver$conditional_rope_prob[r_i] <- if (vtot > 0L) k_rope / vtot else NA_real_
    ver$rope_excl_prob[r_i] <- if (sp$failure_policy == "effective") (k_rope / ntot) else ver$conditional_rope_prob[r_i]
    mc_rope <- wilson_mc_ci(k_rope, denom_mc, conf_level = sp$mc_conf_level)
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

    valid_gate_ok <- if (sp$failure_policy == "effective") TRUE else ((vtot / ntot) >= sp$min_valid_rate)
    ver$criterion_met[r_i] <- valid_gate_ok && check_criterion_met(
      sp, ver$power[r_i], ver$aipe_prob[r_i], ver$rope_excl_prob[r_i], ver$joint_prob[r_i]
    )
    ver$decision_status[r_i] <- compute_decision_status(
      sp, valid_gate_ok, mc_pow, mc_aipe, mc_rope, mc_joint
    )
  }

  ver
}

#' Chunked Reproducible Parallel Evaluator Over Candidate Lattice Points
#' @keywords internal
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
    sing <- rep(FALSE, nb)
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
  sing_mat <- matrix(FALSE, nrow = n_pts, ncol = nsim)
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

split_indices <- function(n, k) {
  if (k <= 1L) return(list(seq_len(n)))
  grp <- cut(seq_len(n), breaks = k, labels = FALSE)
  unname(split(seq_len(n), grp))
}

start_mirai_pool <- function(workers) {
  tryCatch({
    mirai::daemons(workers)
    is_dev <- isNamespaceLoaded("pkgload") &&
      isTRUE(tryCatch(pkgload::is_dev_package("glmmplanr"), error = function(e) FALSE))
    if (is_dev) {
      pkg_root <- getNamespaceInfo(asNamespace("glmmplanr"), "path")
      init_res <- mirai::everywhere(
        pkgload::load_all(p, compile = FALSE, quiet = TRUE),
        p = pkg_root
      )[]
      if (any(vapply(init_res, mirai::is_error_value, logical(1L)))) {
        try(mirai::daemons(0L), silent = TRUE)
        return(FALSE)
      }
    }
    TRUE
  }, error = function(e) {
    try(mirai::daemons(0L), silent = TRUE)
    FALSE
  })
}

run_parallel_tasks <- function(
  tasks,
  worker_fn,
  workers,
  parallel_backend,
  stage_label = NULL,
  verbose = FALSE
) {
  n_tasks <- length(tasks)
  if (n_tasks == 0L) return(list())

  task_weights <- vapply(tasks, function(t) {
    if (is.list(t) && !is.null(t$rep_idx)) length(t$rep_idx) else 1L
  }, integer(1L))
  total_reps <- sum(task_weights)
  lbl <- if (!is.null(stage_label) && nzchar(stage_label)) stage_label else "Simulating"

  show_pb <- isTRUE(verbose) && total_reps > 1L && n_tasks > 1L
  pb_id <- NULL
  if (show_pb) {
    pb_id <- cli::cli_progress_bar(
      name = lbl,
      total = total_reps,
      format = paste0(
        "{cli::pb_spin} {cli::pb_name} [{cli::pb_bar}] ",
        "{cli::pb_current}/{cli::pb_total} sims ({cli::pb_percent}) | ",
        "Elapsed: {cli::pb_elapsed} | ETA: {cli::pb_eta}"
      ),
      clear = FALSE
    )
    on.exit(cli::cli_progress_done(id = pb_id), add = TRUE)
  }

  run_seq <- function(idx_vec) {
    out <- vector("list", length(idx_vec))
    for (k in seq_along(idx_vec)) {
      ii <- idx_vec[k]
      out[[k]] <- worker_fn(tasks[[ii]])
      if (!is.null(pb_id)) cli::cli_progress_update(inc = task_weights[ii], id = pb_id)
    }
    out
  }

  if (workers <= 1L || parallel_backend == "sequential" || n_tasks <= 1L) {
    return(run_seq(seq_len(n_tasks)))
  }

  chk <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  if (nzchar(chk) && chk != "false") {
    workers <- min(workers, 2L)
  }

  n_waves <- if (show_pb) min(n_tasks, max(1L, ceiling(n_tasks / max(workers * 2L, 2L)))) else 1L
  waves <- split_indices(n_tasks, n_waves)

  is_macos <- identical(Sys.info()[["sysname"]], "Darwin")
  prefer_mirai <- parallel_backend == "mirai" ||
    (parallel_backend == "auto" &&
       (isTRUE(mirai::daemons_set()) || .Platform$OS.type != "unix" || is_macos))

  if (prefer_mirai) {
    res <- tryCatch({
      started_daemons <- FALSE
      if (!isTRUE(mirai::daemons_set())) {
        if (!start_mirai_pool(workers)) stop("mirai pool init failed", call. = FALSE)
        started_daemons <- TRUE
      }
      on.exit(if (started_daemons) try(mirai::daemons(0L), silent = TRUE), add = TRUE)
      full_out <- vector("list", n_tasks)
      for (wv in waves) {
        sub_out <- mirai::mirai_map(
          tasks[wv], function(tsk, fn) fn(tsk), .args = list(fn = worker_fn)
        )[]
        if (any(vapply(sub_out, mirai::is_error_value, logical(1L)))) {
          stop("mirai worker error", call. = FALSE)
        }
        full_out[wv] <- sub_out
        if (!is.null(pb_id)) cli::cli_progress_update(inc = sum(task_weights[wv]), id = pb_id)
      }
      full_out
    }, error = function(e) NULL)
    if (!is.null(res)) return(res)
  }

  if (.Platform$OS.type == "unix") {
    old_omp <- Sys.getenv("OMP_NUM_THREADS", unset = NA_character_)
    old_vec <- Sys.getenv("VECLIB_MAXIMUM_THREADS", unset = NA_character_)
    Sys.setenv(OMP_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
    on.exit({
      if (is.na(old_omp)) Sys.unsetenv("OMP_NUM_THREADS") else Sys.setenv(OMP_NUM_THREADS = old_omp)
      if (is.na(old_vec)) Sys.unsetenv("VECLIB_MAXIMUM_THREADS") else Sys.setenv(VECLIB_MAXIMUM_THREADS = old_vec)
    }, add = TRUE)

    full_out <- vector("list", n_tasks)
    ok_mc <- TRUE
    for (wv in waves) {
      use_presched <- length(wv) > (workers * 4L)
      sub_out <- parallel::mclapply(
        tasks[wv], worker_fn, mc.cores = workers, mc.preschedule = use_presched
      )
      bad <- vapply(sub_out, function(x) inherits(x, "try-error") || is.null(x), logical(1L))
      if (any(bad)) {
        ok_mc <- FALSE
        break
      }
      full_out[wv] <- sub_out
      if (!is.null(pb_id)) cli::cli_progress_update(inc = sum(task_weights[wv]), id = pb_id)
    }
    if (ok_mc) return(full_out)
  }

  run_seq(seq_len(n_tasks))
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
  in_b <- if ("within_budget" %in% names(surf)) {
    !is.na(surf$within_budget) & surf$within_budget
  } else {
    rep(TRUE, nrow(surf))
  }

  over_feas <- which(!in_b & surf$pred_criterion_met)
  over_best <- if (!any(in_b & surf$pred_criterion_met) && length(over_feas) > 0L && ("total_cost" %in% names(surf))) {
    over_feas[order(surf$total_cost[over_feas], surf$preference_score[over_feas])[1L]]
  } else {
    integer(0)
  }

  # Extract predicted feasible Pareto frontier + high-confidence neighbors
  pf <- extract_pareto_frontier(surf, spec)
  if (nrow(pf) > 0L) {
    pf_idx <- match(paste(pf$n_subj, pf$n_item), paste(surf$n_subj, surf$n_item))
    # 1. Best under requested `prefer`
    best_pref <- pf_idx[order(surf$preference_score[pf_idx], -surf$prob_feasible[pf_idx])]
    # 2. High-confidence neighbor (`prob_feasible >= 0.80`) along `prefer` within budget
    hc <- which(surf$prob_feasible >= 0.80 & in_b)
    ord_hc <- hc[order(surf$preference_score[hc], -surf$prob_feasible[hc])]
    # 3. Spread points along Pareto frontier (fewer_items end, balanced knee, fewer_subj end)
    spread_idx <- unique(pf_idx[round(seq(1, length(pf_idx), length.out = min(length(pf_idx), n_verify)))])

    pool <- unique(c(
      best_pref[1L],
      over_best,
      if (length(ord_hc) > 0L) ord_hc[1L] else integer(0),
      spread_idx,
      best_pref,
      ord_hc
    ))
    return(pool[seq_len(min(n_verify, length(pool)))])
  }

  plausible <- which(surf$prob_feasible >= 0.35 & in_b)
  if (length(plausible) > 0L) {
    ord_pl <- unique(c(
      plausible[order(surf$preference_score[plausible], -surf$prob_feasible[plausible])],
      over_best
    ))
    return(ord_pl[seq_len(min(n_verify, length(ord_pl)))])
  }

  ord <- unique(c(
    over_best,
    order(-in_b, -surf$prob_feasible, -surf$pred_power, surf$pred_ci_width, surf$preference_score)
  ))
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
#' @keywords internal
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
  in_b <- if ("within_budget" %in% names(surf)) {
    !is.na(surf$within_budget) & surf$within_budget
  } else {
    rep(TRUE, nrow(surf))
  }
  feas <- surf[surf$pred_criterion_met & in_b, , drop = FALSE]
  if (nrow(feas) > 0L) {
    return(extract_2d_pareto(feas))
  }
  # Fallback when no lattice point reaches threshold within budget: return outer eligible edge
  elig <- if (any(in_b)) surf[in_b, , drop = FALSE] else surf
  if ("within_budget" %in% names(surf) && any(in_b) && !all(in_b)) {
    # Budget boundary is itself the outer feasible frontier
    ord_c <- order(elig$n_item, -elig$n_subj)
    e_ord <- elig[ord_c, , drop = FALSE]
    keep <- !duplicated(e_ord$n_item)
    out <- e_ord[keep, , drop = FALSE]
    out <- out[order(out$n_item, out$n_subj), , drop = FALSE]
    rownames(out) <- NULL
    return(out)
  }
  max_s <- max(elig$n_subj)
  max_i <- max(elig$n_item)
  out <- elig[elig$n_subj == max_s | elig$n_item == max_i, , drop = FALSE]
  out <- out[order(out$n_item, out$n_subj), , drop = FALSE]
  rownames(out) <- NULL
  out
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
  verbose = FALSE
) {
  mults <- c(optimistic = 0.70, standard = 1.00, conservative = 1.35)
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
  cols_front <- c("scenario", "slope_multiplier", "n_subj", "n_item", "true_beta")
  out[, c(cols_front, setdiff(names(out), cols_front)), drop = FALSE]
}

#' @export
print.glmmplanr_result <- function(x, ...) {
  cli::cli_h2("glmmplanr 2D Sample Size Optimization Result")
  opt <- x$optimal
  unit_lbl <- switch(x$spec$scale, logit = "logit", pp = "pp", prob = "prob")

  status_txt <- switch(
    opt$decision_status[1L],
    pass = if (isTRUE(x$feasible)) {
      "PASS (High-confidence verified attainment; MC lower bound >= target)"
    } else {
      "OVER BUDGET (Target met, but exceeds specified study budget)"
    },
    indeterminate = "INDETERMINATE / BORDERLINE (MC confidence interval straddles target; consider top-up)",
    fail = if (isTRUE(x$feasible)) {
      "FEASIBLE POINT ESTIMATE (Check MC interval)"
    } else {
      "No confirmed design among evaluated candidates"
    }
  )

  role <- if (isTRUE(x$feasible)) "Recommended design" else "Unconfirmed candidate"
  cost_str <- if (!is.null(x$cost) && ("total_cost" %in% names(opt)) && is.finite(opt$total_cost)) {
    if (!is.null(x$budget)) {
      sprintf(" | Total Cost = %s (Budget = %s)",
              format(round(opt$total_cost, 2), big.mark = ","),
              format(round(x$budget, 2), big.mark = ","))
    } else {
      sprintf(" | Total Cost = %s", format(round(opt$total_cost, 2), big.mark = ","))
    }
  } else {
    ""
  }

  bullets <- c(
    "Decision status: {.strong {status_txt}}",
    "{role} ({.code prefer = '{x$prefer}'}): N_subj = {.strong {opt$n_subj}}, N_item = {.strong {opt$n_item}} (Scheduled: {opt$n_trials_scheduled} trials | Expected post-attrition: {opt$n_trials_expected}{cost_str})",
    "Power ({x$spec$failure_policy}) = {.strong {sprintf('%.1f%%', opt$power * 100)}} [{sprintf('%.1f%%', opt$power_mc_lo * 100)}, {sprintf('%.1f%%', opt$power_mc_hi * 100)}] | Conditional Power = {ifelse(is.na(opt$conditional_power), 'NA', sprintf('%.1f%%', opt$conditional_power * 100))} (Valid fits: {sprintf('%.1f%%', opt$valid_rate * 100)})",
    "Accuracy diff ({x$spec$prob_estimand}): true = {sprintf('%+.1f pp', opt$prob_diff_true_pp)} | Mean CI width = {.strong {sprintf('%.2f pp', opt$prob_ci_width_mean_pp)}} ({sprintf('%.3f logit', opt$logit_ci_width_mean)})"
  )
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
  if (!is.null(x$min_cost_over_budget) && !isTRUE(x$feasible) && !is.null(x$budget)) {
    ov <- x$min_cost_over_budget
    bullets <- c(
      bullets,
      sprintf(
        "Minimum-cost verified PASS over budget: N_subj = %d, N_item = %d (Cost = %s vs Budget = %s)",
        ov$n_subj, ov$n_item,
        format(round(ov$total_cost, 2), big.mark = ","),
        format(round(x$budget, 2), big.mark = ",")
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
  cov_show <- if (x$spec$scale == "logit") opt$coverage else opt$prob_coverage
  bullets <- c(
    bullets,
    "Diagnostics: Singular = {sprintf('%.1f%%', opt$singular_rate * 100)} | Non-conv/Fail = {sprintf('%.1f%%', opt$nonconv_rate * 100)} (MaxIter = {sprintf('%.1f%%', opt$max_iter_rate * 100)}) | Coverage = {ifelse(is.na(cov_show), 'NA', sprintf('%.1f%%', cov_show * 100))}"
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
      dc_txt <- if (!is.null(rc$delta_cost) && is.finite(rc$delta_cost)) {
        sprintf(", %+.2f cost", rc$delta_cost)
      } else {
        ""
      }
      bullets <- c(
        bullets,
        "All-scenario robustness increment (vs {.val {rc$standard_scenario}}): {sprintf('%+d', rc$delta_subj)} subjects, {sprintf('%+d', rc$delta_item)} items ({sprintf('%+d', rc$delta_trials)} scheduled trials{dc_txt})"
      )
    }
  }
  bullets <- c(
    bullets,
    "Search time: {round(x$search_meta$elapsed_sec, 2)} sec ({nrow(x$evaluated)} point evaluations, method = {.val {x$search_meta$method}})"
  )
  cli::cli_ul(bullets)
  invisible(x)
}

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
      if ("total_cost" %in% names(object$verified_pareto)) "total_cost",
      "effective_power", "prob_ci_width_mean_pp", "decision_status"
    )
    vp <- object$verified_pareto[, cols_vp, drop = FALSE]
    vp$effective_power <- round(vp$effective_power, 3)
    vp$prob_ci_width_mean_pp <- round(vp$prob_ci_width_mean_pp, 2)
    if ("total_cost" %in% names(vp)) vp$total_cost <- round(vp$total_cost, 2)
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
