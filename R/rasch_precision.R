#' Plan Sample Sizes for Random-Item Rasch Parameter Recovery
#'
#' Evaluates each proposed allocation of people and items with the same set of
#' error tolerances. Fits a crossed random-intercept logistic model to each
#' simulated study using either `lme4::glmer()` (default) or the fast native
#' C++ Laplace engine (`"cpp_laplace"`). Person ability has population mean
#' zero. Item difficulty is the negative sum of the fixed and item intercepts.
#'
#' @param dgp A [calibrate_dgp()] model with `chance = 0`, no experimental
#'   factors, and participant and item intercepts. Each replication draws a
#'   new sample of both populations. Partial presentation uses the design's
#'   fixed allocation seed. Missingness is drawn anew.
#' @param designs Data frame with integer `n_subj` and `n_item` columns, each
#'   at least two. Use `expand.grid()` for a rectangular grid or supply pairs.
#'   Distinct pairs are evaluated in ascending participant and item order.
#' @param max_error Named vector of nonnegative tolerances. For `ability` and
#'   `difficulty` the error in a replication is the root mean squared error of
#'   the conditional-mode predictions across all planned people or items.
#'   `difficulty_max_abs_error` is the largest absolute difficulty error
#'   across all planned items in that replication.
#'   For `person_sd`, `item_sd`, `person_variance` and `item_variance` it is
#'   the absolute error of the fitted population SD or variance. Variance
#'   tolerances are in squared logits. The other tolerances are in logits.
#' @param assurance Required probability that a study meets its tolerance(s).
#' @param nsim Number of independent studies per design.
#' @param seed Integer master simulation seed.
#' @param mc_conf_level Pointwise confidence level of Wilson Monte Carlo intervals.
#' @param engine Estimation backend: `"lme4"` (default; `lme4::glmer()` with
#'   `"bobyqa"`) or `"cpp_laplace"` (fast native C++ Laplace / Empirical Bayes
#'   solver).
#' @param workers Integer number of parallel workers (default `1L`). Results are
#'   strictly worker-invariant via per-replication L'Ecuyer-CMRG substreams.
#' @param parallel_backend Parallel backend (`"auto"`, `"mirai"`, `"multicore"`,
#'   or `"sequential"`).
#' @param verbose Logical; display live progress bar with ETA (default `FALSE`).
#' @return A `glmmplanr_rasch` list (see [rasch_plan_methods]) with `summary` (per-target and joint attainment, counts, Wilson
#'   bounds, and decision), `pareto` (non-dominated designs passing the joint
#'   criterion), `replications` (errors and fit diagnostics), and the settings.
#'   The `rmse` column is the root mean square of the evaluable replication
#'   errors. A missing prediction for any planned unit makes that target
#'   unevaluable. Failed fits and unevaluable errors count as misses in the
#'   assurance denominator. Converged singular fits are included. Joint success
#'   requires every selected tolerance to be met in the same replication.
#'   A pass requires the lower Monte Carlo bound to reach `assurance`.
#' @examples
#' design <- design_crossed(list(), allocation = "crossed")
#' dgp <- calibrate_dgp(design, 0.65, subj_sd = 0.8, item_sd = 0.6)
#' plan <- plan_rasch_precision(dgp, data.frame(n_subj = 48, n_item = 12),
#'                              c(difficulty = 0.30), nsim = 4)
#' as.data.frame(plan)
#' @export
plan_rasch_precision <- function(
  dgp,
  designs,
  max_error,
  assurance = 0.80,
  nsim = 200L,
  seed = 2026L,
  mc_conf_level = 0.95,
  engine = c("lme4", "cpp_laplace"),
  workers = 1L,
  parallel_backend = c("auto", "mirai", "multicore", "sequential"),
  verbose = FALSE
) {
  if (uses_box_muller()) {
    return(isolated_seeded_call("plan_rasch_precision", environment(), match.call()))
  }
  if (!inherits(dgp, "glmmplanr_dgp") || dgp$chance != 0 ||
      length(dgp$design$factor_names) || dgp$design$has_testlet) {
    cli::cli_abort("Supply a Rasch DGP with {.code factors = list()}, {.code chance = 0}, and person and item intercepts.")
  }
  engine <- match.arg(engine)
  parallel_backend <- match.arg(parallel_backend)
  workers <- max(1L, as.integer(workers))

  targets <- c("ability", "difficulty", "difficulty_max_abs_error", "person_sd",
               "item_sd", "person_variance", "item_variance")
  if (!is.numeric(max_error) || !is.null(dim(max_error)) || !length(max_error) ||
      is.null(names(max_error)) || anyDuplicated(names(max_error)) ||
      !all(names(max_error) %in% targets) || any(!is.finite(max_error) | max_error < 0)) {
    cli::cli_abort("{.arg max_error} needs uniquely named nonnegative finite tolerances from {.val {targets}}.")
  }
  counts <- c("n_subj", "n_item")
  if (!is.data.frame(designs) || !nrow(designs) || !all(counts %in% names(designs))) {
    cli::cli_abort("{.arg designs} must include at least one row with {.field n_subj} and {.field n_item}.")
  }
  valid_count <- function(x, minimum) is.numeric(x) && is.null(dim(x)) &&
    all(is.finite(x) & x == trunc(x) & x >= minimum & x <= .Machine$integer.max)
  designs <- unique(as.data.frame(designs)[, counts, drop = FALSE])
  if (!all(vapply(designs, valid_count, logical(1), minimum = 2)) ||
      any(designs$n_item < dgp$design$items_per_subject)) {
    cli::cli_abort("Candidate counts must be integers of at least two and banks must cover {.arg items_per_subject}.")
  }
  if (length(nsim) != 1L || !valid_count(nsim, 1) ||
      length(seed) != 1L || !valid_count(seed, 0)) {
    cli::cli_abort("Supply a positive integer {.arg nsim} and a nonnegative integer {.arg seed}.")
  }
  nsim <- as.integer(nsim)
  seed <- as.integer(seed)
  for (p in list(assurance, mc_conf_level)) {
    if (!is.numeric(p) || length(p) != 1L || !is.finite(p) || p <= 0 || p >= 1) {
      cli::cli_abort("{.arg assurance} and {.arg mc_conf_level} must be in (0, 1).")
    }
  }
  designs <- designs[order(designs$n_subj, designs$n_item), , drop = FALSE]
  n_pts <- nrow(designs)
  points <- allocate_stage_point_streams(seed, 1L, seq_len(n_pts), n_pts)
  control <- if (engine == "lme4") lme4::glmerControl(optimizer = "bobyqa") else NULL

  all_streams <- lapply(seq_len(n_pts), function(i) derive_rep_substreams(points[[i]], nsim))
  target_tasks <- max(n_pts, workers * 4L, if (isTRUE(verbose)) 16L else 1L)
  chunks_per_pt <- max(1L, min(nsim, ceiling(target_tasks / n_pts)))
  local_bounds <- split_indices(nsim, chunks_per_pt)

  tasks <- list()
  tid <- 0L
  for (i in seq_len(n_pts)) {
    for (ch in seq_along(local_bounds)) {
      ix <- local_bounds[[ch]]
      tid <- tid + 1L
      tasks[[tid]] <- list(
        pt_idx = i,
        rep_idx = ix,
        n_subj = as.integer(designs$n_subj[i]),
        n_item = as.integer(designs$n_item[i]),
        streams = all_streams[[i]][ix]
      )
    }
  }

  worker_fn <- function(tsk) {
    schedule <- generate_design_schedule(dgp$design, tsk$n_subj, tsk$n_item)
    rasch_ctx <- if (identical(engine, "cpp_laplace")) {
      build_rasch_cpp_ctx(schedule, tsk$n_subj, tsk$n_item, dgp$design)
    } else {
      NULL
    }
    rows <- lapply(tsk$streams, function(stream) {
      with_rng_stream(stream, {
        draw <- if (identical(engine, "cpp_laplace")) {
          draw_rasch_trials_light(dgp, schedule)
        } else {
          draw_dgp_trials(dgp, schedule, keep_person_item_effects = TRUE)
        }
        rasch_recovery_errors(
          draw,
          dgp,
          control,
          engine = engine,
          rasch_ctx = rasch_ctx
        )
      })
    })
    list(pt_idx = tsk$pt_idx, rep_idx = tsk$rep_idx, df = do.call(rbind, rows))
  }

  chunk_out <- run_parallel_tasks(
    tasks, worker_fn, workers = workers, parallel_backend = parallel_backend,
    stage_label = "Rasch Precision", verbose = verbose
  )

  by_pt <- vector("list", n_pts)
  for (cr in chunk_out) {
    p_i <- cr$pt_idx
    if (is.null(by_pt[[p_i]])) by_pt[[p_i]] <- vector("list", nsim)
    for (k in seq_along(cr$rep_idx)) {
      by_pt[[p_i]][[cr$rep_idx[k]]] <- cr$df[k, , drop = FALSE]
    }
  }

  summaries <- samples <- vector("list", n_pts)
  for (i in seq_len(n_pts)) {
    size <- designs[i, , drop = FALSE]
    runs <- do.call(rbind, by_pt[[i]])
    rownames(runs) <- NULL
    losses <- runs[names(max_error)]
    meets <- sweep(as.matrix(losses), 2L, max_error, FUN = "<=")
    evaluable <- is.finite(as.matrix(losses))
    meets[!evaluable] <- FALSE
    meets <- cbind(meets, joint = rowSums(meets) == ncol(losses))
    hits <- colSums(meets)
    intervals <- vapply(hits, wilson_mc_ci, numeric(2), n = nsim, conf_level = mc_conf_level)
    rmse <- vapply(losses, function(x) {
      if (any(is.finite(x))) sqrt(mean(x[is.finite(x)]^2)) else NA_real_
    }, numeric(1))
    summaries[[i]] <- data.frame(
      n_subj = size$n_subj, n_item = size$n_item,
      target = names(hits), max_error = unname(c(max_error, NA_real_)), nsim = nsim,
      n_valid_fit = sum(is.na(runs$fit_error)), n_singular = sum(runs$singular, na.rm = TRUE),
      n_evaluable = c(colSums(evaluable), sum(rowSums(evaluable) == ncol(losses))),
      rmse = c(rmse, NA_real_), n_met = hits, attainment = hits / nsim,
      mc_lo = intervals[1L, ], mc_hi = intervals[2L, ],
      decision = mapply(classify_mc_threshold, intervals[1L, ], intervals[2L, ], assurance),
      row.names = NULL
    )
    samples[[i]] <- data.frame(
      n_subj = size$n_subj, n_item = size$n_item,
      replication = seq_len(nsim), losses, singular = runs$singular, fit_error = runs$fit_error
    )
  }
  summary <- do.call(rbind, summaries)
  passing <- summary$target == "joint" & summary$decision == "pass"
  structure(list(
    summary = summary,
    engine = engine,
    pareto = extract_2d_pareto(summary[passing, , drop = FALSE]),
    replications = do.call(rbind, samples),
    max_error = max_error,
    assurance = assurance,
    seed = seed,
    mc_conf_level = mc_conf_level
  ), class = "glmmplanr_rasch")
}

build_rasch_cpp_ctx <- function(schedule, n_subj, n_item, design) {
  no_missing <- (design$attrition == 0 && design$participant_dropout == 0)
  if (!no_missing) return(NULL)
  tr <- schedule$trials
  n_obs <- nrow(tr)
  if (n_obs == 0L) return(NULL)
  ones <- matrix(1.0, nrow = n_obs, ncol = 1L)
  ws <- cpp_create_glmm_workspace(
    X = ones,
    subj_idx = as.integer(tr$subj_id) - 1L,
    Z_subj = ones,
    n_subj = as.integer(n_subj),
    item_idx = as.integer(tr$item_id) - 1L,
    Z_item = ones,
    n_item = as.integer(n_item),
    testlet_idx = integer(n_obs),
    Z_testlet = matrix(0.0, nrow = n_obs, ncol = 0L),
    n_testlet = 0L
  )
  list(
    ws = ws,
    unobs_subj = tabulate(tr$subj_id, nbins = n_subj) == 0L,
    unobs_item = tabulate(tr$item_id, nbins = n_item) == 0L
  )
}

draw_rasch_trials_light <- function(dgp, schedule) {
  des <- dgp$design
  ns <- schedule$n_subj
  ni <- schedule$n_item
  u_s <- matrix(stats::rnorm(ns) * dgp$L_subj[1L, 1L], nrow = ns, ncol = 1L)
  u_i <- matrix(stats::rnorm(ni) * dgp$L_item[1L, 1L], nrow = ni, ncol = 1L)
  tr <- schedule$trials
  sid <- tr$subj_id
  iid <- tr$item_id
  eta <- dgp$mu_cells[tr$cell_id] + u_s[sid, 1L] + u_i[iid, 1L]
  y <- stats::rbinom(length(eta), size = 1L, prob = stats::plogis(eta))
  keep <- observed_trials(des, sid, ns)
  trials_out <- if (!is.null(keep)) {
    list(subj_id = sid[keep], item_id = iid[keep], y = y[keep], has_keep = TRUE)
  } else {
    list(subj_id = sid, item_id = iid, y = y, has_keep = FALSE)
  }
  list(trials = trials_out, subj = u_s, item = u_i)
}

fit_rasch_random_intercept_model <- function(trials, n_subj, n_item, engine = "lme4", control = NULL, rasch_ctx = NULL) {
  if (identical(engine, "cpp_laplace")) {
    n_obs <- length(trials$y)
    use_ctx <- !is.null(rasch_ctx) && !isTRUE(trials$has_keep) && cpp_glmm_workspace_valid(rasch_ctx$ws)
    fit <- tryCatch(
      if (use_ctx) {
        cpp_fit_crossed_glmm_ws(
          ws_ptr = rasch_ctx$ws,
          y = as.integer(trials$y),
          chance = 0.0,
          corr_mode = 0L,
          regularized = FALSE,
          lkj_eta = 2.0,
          prior_tau_sd = 0.0,
          prior_beta_sd = 0.0,
          max_em_iter = 160L,
          max_pirls_iter = 25L,
          tol = 1e-4,
          fix_cov = FALSE,
          compute_joint_vcov = FALSE,
          return_ranef = TRUE
        )
      } else {
        ones <- matrix(1.0, nrow = n_obs, ncol = 1L)
        cpp_fit_crossed_glmm(
          y = as.integer(trials$y),
          X = ones,
          Z_subj = ones,
          subj_idx = as.integer(trials$subj_id) - 1L,
          n_subj = n_subj,
          Z_item = ones,
          item_idx = as.integer(trials$item_id) - 1L,
          n_item = n_item,
          Z_testlet = matrix(0.0, nrow = n_obs, ncol = 0L),
          testlet_idx = integer(n_obs),
          n_testlet = 0L,
          chance = 0.0,
          corr_mode = 0L,
          regularized = FALSE,
          lkj_eta = 2.0,
          prior_tau_sd = 0.0,
          prior_beta_sd = 0.0,
          max_em_iter = 160L,
          max_pirls_iter = 25L,
          tol = 1e-4,
          fix_cov = FALSE,
          compute_joint_vcov = FALSE,
          return_ranef = TRUE
        )
      },
      error = identity
    )
    if (inherits(fit, "error")) {
      return(list(error = conditionMessage(fit)))
    }
    if (!isTRUE(fit$converged)) {
      return(list(error = paste0("C++ Laplace status: ", fit$status)))
    }
    ability <- as.numeric(fit$u_subj[, 1L])
    item <- as.numeric(fit$u_item[, 1L])
    if (use_ctx) {
      ability[rasch_ctx$unobs_subj] <- NA_real_
      item[rasch_ctx$unobs_item] <- NA_real_
    } else {
      ability[tabulate(trials$subj_id, nbins = n_subj) == 0L] <- NA_real_
      item[tabulate(trials$item_id, nbins = n_item) == 0L] <- NA_real_
    }
    v_s <- as.numeric(fit$Sigma_subj[1L, 1L])
    v_i <- as.numeric(fit$Sigma_item[1L, 1L])
    if (v_s <= 1e-8) v_s <- 0.0
    if (v_i <= 1e-8) v_i <- 0.0
    return(list(
      error = NA_character_,
      singular = isTRUE(fit$singular),
      intercept = as.numeric(fit$beta[1L]),
      ability = ability,
      item = item,
      person_var = v_s,
      item_var = v_i
    ))
  }

  if (is.null(control)) control <- lme4::glmerControl(optimizer = "bobyqa")
  warnings <- character()
  fit <- tryCatch(withCallingHandlers(
    suppressMessages(lme4::glmer(y ~ 1 + (1 | subj) + (1 | item), data = trials,
                                 family = stats::binomial(), control = control)),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = identity)
  if (inherits(fit, "error")) {
    return(list(error = conditionMessage(fit)))
  }
  conv <- fit@optinfo$conv
  diagnostics <- c(warnings, conv$lme4$messages)
  diagnostics <- diagnostics[!grepl("boundary \\(singular\\) fit", diagnostics)]
  if (any(c(conv$opt, conv$lme4$code) != 0) || length(diagnostics)) {
    return(list(error = paste(c("Fit diagnostic", unique(diagnostics)), collapse = " / ")))
  }
  modes <- lme4::ranef(fit, condVar = FALSE)
  # Preserve planned integer IDs when glmer drops an unobserved level.
  ability <- modes$subj[match(seq_len(n_subj), rownames(modes$subj)), 1L]
  item <- modes$item[match(seq_len(n_item), rownames(modes$item)), 1L]
  variance <- lme4::VarCorr(fit)
  list(
    error = NA_character_,
    singular = lme4::isSingular(fit),
    intercept = unname(lme4::fixef(fit)[1L]),
    ability = ability,
    item = item,
    person_var = as.numeric(variance$subj[1L, 1L]),
    item_var = as.numeric(variance$item[1L, 1L])
  )
}

rasch_recovery_errors <- function(draw, dgp, control = NULL, engine = "lme4", rasch_ctx = NULL) {
  data <- draw$trials
  result <- data.frame(
    ability = NA_real_, difficulty = NA_real_, difficulty_max_abs_error = NA_real_,
    person_sd = NA_real_, item_sd = NA_real_, person_variance = NA_real_, item_variance = NA_real_,
    singular = NA, fit_error = NA_character_
  )
  y_vec <- data$y
  if (length(y_vec) < 2L || !any(y_vec == 0L) || !any(y_vec == 1L)) {
    result$fit_error <- "Fewer than two observed response values"
    return(result)
  }

  n_subj <- nrow(draw$subj)
  n_item <- nrow(draw$item)
  fit <- fit_rasch_random_intercept_model(
    data, n_subj, n_item, engine = engine, control = control, rasch_ctx = rasch_ctx
  )
  if (!is.na(fit$error)) {
    result$fit_error <- fit$error
    return(result)
  }

  result$singular <- fit$singular
  result$ability <- sqrt(mean((fit$ability - draw$subj[, 1L])^2))
  difficulty <- -fit$intercept - fit$item
  truth <- -unname(dgp$beta[1L]) - draw$item[, 1L]
  difficulty_error <- difficulty - truth
  result$difficulty <- sqrt(mean(difficulty_error^2))
  result$difficulty_max_abs_error <- max(abs(difficulty_error))
  result$person_variance <- abs(fit$person_var - dgp$tau_subj[1L]^2)
  result$item_variance <- abs(fit$item_var - dgp$tau_item[1L]^2)
  result$person_sd <- abs(sqrt(fit$person_var) - dgp$tau_subj[1L])
  result$item_sd <- abs(sqrt(fit$item_var) - dgp$tau_item[1L])
  result
}

#' Plan Examinee Sample Sizes for Calibrating a Fixed Rasch Item Bank
#'
#' Evaluates examinee sample sizes (`n_subj`) for calibrating a fixed,
#' pre-specified vector of item difficulties `item_difficulties` (held constant
#' across replications, whereas examinee abilities \eqn{\theta_i \sim N(0, \sigma_p^2)}
#' are drawn anew in every replication). Supports both fixed-item Marginal
#' Maximum Likelihood (`method = "fixed_item"`) and random-item Empirical Bayes
#' shrinkage (`method = "random_item"`).
#'
#' In fixed-item calibration, an examinee with all responses correct or all
#' incorrect contributes to the likelihood if they answered at least one item
#' with both response categories. Their ability is scored by the posterior mean
#' (EAP) under the fitted normal ability distribution. A person with no such
#' responses has missing ability. A constant item has missing difficulty, and
#' an item with no responses has missing difficulty under either method.
#' A missing estimate counts as a miss for each target that requires it.
#' If the entire dataset contains fewer than two observed response categories,
#' the replication fails under both methods. All planned replications remain
#' in the assurance denominator.
#'
#' The zero mean of the examinee population identifies the origin. Comparisons
#' across populations with unknown mean shifts require a linking design, for
#' example common anchor items.
#'
#' @param item_difficulties Numeric vector of at least two finite true item
#'   difficulties (in logits) defining the fixed item bank.
#' @param n_subj Integer vector (or data frame containing `n_subj`) of candidate
#'   examinee sample sizes (\eqn{\ge 2}).
#' @param max_error Named vector of nonnegative tolerances from
#'   `c("ability", "difficulty", "difficulty_max_abs_error", "person_sd", "person_variance")`.
#'   `difficulty` is the bank-wide RMSE across all `length(item_difficulties)`
#'   items; `difficulty_max_abs_error` is the worst-item absolute error across
#'   the bank.
#' @param person_sd Nonnegative numeric population SD of examinee ability
#'   \eqn{\theta_i \sim N(0, \sigma_p^2)} (default `0.8`).
#' @param items_per_subject Optional integer number of distinct items assigned
#'   per examinee for incomplete-block / matrix-sampling designs (`NULL` = full
#'   crossing).
#' @param item_seed Positive integer seed governing deterministic incomplete
#'   item-block assignment (default `1L`).
#' @param attrition Trial-level random missingness rate in `[0, 1)` (default `0`).
#' @param participant_dropout Participant-level dropout rate in `[0, 1)` (default `0`).
#' @param method Calibration model:
#'   * `"fixed_item"` (default): Marginal Maximum Likelihood with fixed item
#'     parameters identified by centering person abilities at mean zero
#'     (`constraint = "cases"`). Items with 0% or 100% observed correct
#'     responses have no finite fixed-item MLE (`NA_real_`), rendering bank-wide
#'     difficulty targets unevaluable (counted as misses in the `nsim`
#'     assurance denominator).
#'   * `"random_item"`: Crossed random-intercept GLMM (`y ~ 1 + (1 | subj) + (1 | item)`),
#'     estimating person and item effects by empirical Bayes conditional modes.
#'     When the fit succeeds, an administered item with zero correct or zero
#'     incorrect responses receives a finite difficulty. Shrinkage toward the
#'     bank mean can be substantial for extreme items.
#' @param engine Estimation engine: `"auto"` (default; selects `"mml"` for
#'   `method = "fixed_item"` and `"cpp_laplace"` for `method = "random_item"`),
#'   `"mml"` (built-in 21-point Gauss-Hermite MML-EM solver), `"tam"`
#'   (`TAM::tam.mml(..., constraint = "cases")`), `"cpp_laplace"`, or `"lme4"`.
#' @param assurance Required probability that a calibration study meets its
#'   tolerance(s) (default `0.80`).
#' @param nsim Number of independent calibration studies per candidate `n_subj`
#'   (default `200L`).
#' @param seed Integer master simulation seed (default `2026L`).
#' @param mc_conf_level Confidence level for Wilson Monte Carlo intervals
#'   (default `0.95`).
#' @param workers Integer number of parallel workers (default `1L`).
#' @param parallel_backend Parallel backend (`"auto"`, `"mirai"`, `"multicore"`,
#'   or `"sequential"`).
#' @param verbose Logical; display live progress bar with ETA (default `FALSE`).
#' @return A `glmmplanr_rasch` list (see [rasch_plan_methods]) with:
#'   * `summary`: Per-target and `joint` attainment, `n_extreme_item_reps`,
#'     `n_evaluable`, `rmse`, Wilson bounds (`mc_lo`, `mc_hi`), and `decision`.
#'   * `recommended_n_subj`: Smallest candidate `n_subj` achieving `PASS` on
#'     `target == "joint"` (`NA_integer_` if none pass).
#'   * `item_summary`: Per-item bias, RMSE, MAE, finite-estimate count, and
#'     all-zero / all-one counts across replications at each `n_subj`.
#'     `n_scheduled` and `component` describe the pre-missingness schedule
#'     (see [generate_design_schedule()]).
#'   * `replications`: Replication-level error and diagnostic table.
#' @examples
#' bank <- seq(-2.0, 2.0, length.out = 12)
#' fb <- plan_rasch_fixed_bank(
#'   item_difficulties = bank,
#'   n_subj = c(40L, 100L),
#'   max_error = c(difficulty = 0.35, difficulty_max_abs_error = 0.60),
#'   method = "fixed_item",
#'   nsim = 10L
#' )
#' summary(fb)
#' @export
plan_rasch_fixed_bank <- function(
  item_difficulties,
  n_subj,
  max_error,
  person_sd = 0.8,
  items_per_subject = NULL,
  item_seed = 1L,
  attrition = 0,
  participant_dropout = 0,
  method = c("fixed_item", "random_item"),
  engine = c("auto", "mml", "tam", "cpp_laplace", "lme4"),
  assurance = 0.80,
  nsim = 200L,
  seed = 2026L,
  mc_conf_level = 0.95,
  workers = 1L,
  parallel_backend = c("auto", "mirai", "multicore", "sequential"),
  verbose = FALSE
) {
  if (uses_box_muller()) {
    return(isolated_seeded_call("plan_rasch_fixed_bank", environment(), match.call()))
  }
  method <- match.arg(method)
  engine <- match.arg(engine)
  parallel_backend <- match.arg(parallel_backend)
  workers <- max(1L, as.integer(workers))

  if (!is.numeric(item_difficulties) || !is.null(dim(item_difficulties)) ||
      length(item_difficulties) < 2L || any(!is.finite(item_difficulties))) {
    cli::cli_abort("{.arg item_difficulties} must be a numeric vector of at least two finite item difficulties.")
  }
  b_true <- as.numeric(item_difficulties)
  n_item <- length(b_true)

  if (!is.numeric(person_sd) || length(person_sd) != 1L || !is.finite(person_sd) || person_sd < 0) {
    cli::cli_abort("{.arg person_sd} must be a single nonnegative finite number.")
  }

  if (engine == "auto") {
    engine <- if (method == "fixed_item") "mml" else "cpp_laplace"
  }
  if (method == "fixed_item" && !(engine %in% c("mml", "tam"))) {
    cli::cli_abort("{.code method = 'fixed_item'} requires {.code engine} in {.val {c('auto', 'mml', 'tam')}}.")
  }
  if (method == "random_item" && !(engine %in% c("cpp_laplace", "lme4"))) {
    cli::cli_abort("{.code method = 'random_item'} requires {.code engine} in {.val {c('auto', 'cpp_laplace', 'lme4')}}.")
  }
  if (engine == "tam" && !requireNamespace("TAM", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg TAM} is required for {.code engine = 'tam'}.")
  }

  targets <- c("ability", "difficulty", "difficulty_max_abs_error", "person_sd", "person_variance")
  if (!is.numeric(max_error) || !is.null(dim(max_error)) || !length(max_error) ||
      is.null(names(max_error)) || anyDuplicated(names(max_error)) ||
      !all(names(max_error) %in% targets) || any(!is.finite(max_error) | max_error < 0)) {
    cli::cli_abort("{.arg max_error} needs uniquely named nonnegative finite tolerances from {.val {targets}}.")
  }

  if (is.data.frame(n_subj)) {
    if (!("n_subj" %in% names(n_subj))) {
      cli::cli_abort("Data frame {.arg n_subj} must contain column {.field n_subj}.")
    }
    n_subj <- n_subj$n_subj
  }
  valid_count <- function(x, minimum) is.numeric(x) && is.null(dim(x)) &&
    all(is.finite(x) & x == trunc(x) & x >= minimum & x <= .Machine$integer.max)
  if (!length(n_subj) || !valid_count(n_subj, 2L)) {
    cli::cli_abort("{.arg n_subj} must contain integers of at least two.")
  }
  n_subj_vec <- sort(unique(as.integer(n_subj)))

  if (!is.null(items_per_subject) && n_item < items_per_subject) {
    cli::cli_abort("Item bank length ({n_item}) must be at least {.arg items_per_subject} ({items_per_subject}).")
  }
  if (length(nsim) != 1L || !valid_count(nsim, 1L) ||
      length(seed) != 1L || !valid_count(seed, 0L)) {
    cli::cli_abort("Supply a positive integer {.arg nsim} and a nonnegative integer {.arg seed}.")
  }
  nsim <- as.integer(nsim)
  seed <- as.integer(seed)
  for (p in list(assurance, mc_conf_level)) {
    if (!is.numeric(p) || length(p) != 1L || !is.finite(p) || p <= 0 || p >= 1) {
      cli::cli_abort("{.arg assurance} and {.arg mc_conf_level} must be in (0, 1).")
    }
  }

  des <- design_crossed(
    factors = list(),
    allocation = "crossed",
    attrition = attrition,
    items_per_subject = items_per_subject,
    item_seed = item_seed,
    participant_dropout = participant_dropout
  )

  n_pts <- length(n_subj_vec)
  points <- allocate_stage_point_streams(seed, 1L, seq_len(n_pts), n_pts)
  all_streams <- lapply(seq_len(n_pts), function(i) derive_rep_substreams(points[[i]], nsim))
  control <- if (engine == "lme4") lme4::glmerControl(optimizer = "bobyqa") else NULL

  target_tasks <- max(n_pts, workers * 4L, if (isTRUE(verbose)) 16L else 1L)
  chunks_per_pt <- max(1L, min(nsim, ceiling(target_tasks / n_pts)))
  local_bounds <- split_indices(nsim, chunks_per_pt)

  tasks <- list()
  tid <- 0L
  for (i in seq_len(n_pts)) {
    for (ch in seq_along(local_bounds)) {
      ix <- local_bounds[[ch]]
      tid <- tid + 1L
      tasks[[tid]] <- list(
        pt_idx = i,
        rep_idx = ix,
        n_subj = n_subj_vec[i],
        streams = all_streams[[i]][ix]
      )
    }
  }

  worker_fn <- function(tsk) {
    ns <- tsk$n_subj
    schedule <- generate_design_schedule(des, ns, n_item, diagnostics = tsk$rep_idx[1L] == 1L)
    rasch_ctx <- if (method == "random_item" && identical(engine, "cpp_laplace")) {
      build_rasch_cpp_ctx(schedule, ns, n_item, des)
    } else {
      NULL
    }
    use_light <- !identical(engine, "lme4")
    nb <- length(tsk$rep_idx)
    rows <- vector("list", nb)
    b_est_mat <- matrix(NA_real_, nrow = nb, ncol = n_item)
    zero_mat <- matrix(FALSE, nrow = nb, ncol = n_item)
    all1_mat <- matrix(FALSE, nrow = nb, ncol = n_item)

    for (k in seq_len(nb)) {
      out_k <- with_rng_stream(tsk$streams[[k]], {
        sim <- simulate_fixed_bank_trials(
          des, schedule, b_true = b_true, person_sd = person_sd, lightweight = use_light
        )
        eval_fixed_bank_replication(
          sim = sim, b_true = b_true, person_sd = person_sd,
          method = method, engine = engine, control = control, rasch_ctx = rasch_ctx
        )
      })
      rows[[k]] <- out_k$metrics
      b_est_mat[k, ] <- out_k$b_est
      zero_mat[k, ] <- out_k$item_zero
      all1_mat[k, ] <- out_k$item_all1
    }
    list(
      pt_idx = tsk$pt_idx,
      rep_idx = tsk$rep_idx,
      df = do.call(rbind, rows),
      item_exposure = schedule$item_exposure,
      b_est = b_est_mat,
      item_zero = zero_mat,
      item_all1 = all1_mat
    )
  }

  chunk_out <- run_parallel_tasks(
    tasks, worker_fn, workers = workers, parallel_backend = parallel_backend,
    stage_label = sprintf("Fixed-Bank Rasch (%s)", method), verbose = verbose
  )

  exposure <- by_pt_rows <- vector("list", n_pts)
  by_pt_b <- lapply(seq_len(n_pts), function(i) matrix(NA_real_, nrow = nsim, ncol = n_item))
  by_pt_z <- lapply(seq_len(n_pts), function(i) matrix(FALSE, nrow = nsim, ncol = n_item))
  by_pt_a <- lapply(seq_len(n_pts), function(i) matrix(FALSE, nrow = nsim, ncol = n_item))

  for (cr in chunk_out) {
    p_i <- cr$pt_idx
    if (is.null(by_pt_rows[[p_i]])) by_pt_rows[[p_i]] <- vector("list", nsim)
    if (!is.null(cr$item_exposure)) exposure[[p_i]] <- cr$item_exposure
    ix <- cr$rep_idx
    for (k in seq_along(ix)) {
      by_pt_rows[[p_i]][[ix[k]]] <- cr$df[k, , drop = FALSE]
    }
    by_pt_b[[p_i]][ix, ] <- cr$b_est
    by_pt_z[[p_i]][ix, ] <- cr$item_zero
    by_pt_a[[p_i]][ix, ] <- cr$item_all1
  }

  summaries <- samples <- item_summaries <- vector("list", n_pts)
  for (i in seq_len(n_pts)) {
    ns <- n_subj_vec[i]
    runs <- do.call(rbind, by_pt_rows[[i]])
    rownames(runs) <- NULL
    losses <- runs[names(max_error)]
    meets <- sweep(as.matrix(losses), 2L, max_error, FUN = "<=")
    evaluable <- is.finite(as.matrix(losses))
    meets[!evaluable] <- FALSE
    meets <- cbind(meets, joint = rowSums(meets) == ncol(losses))
    hits <- colSums(meets)
    intervals <- vapply(hits, wilson_mc_ci, numeric(2), n = nsim, conf_level = mc_conf_level)
    rmse <- vapply(losses, function(x) {
      if (any(is.finite(x))) sqrt(mean(x[is.finite(x)]^2)) else NA_real_
    }, numeric(1))

    summaries[[i]] <- data.frame(
      n_subj = ns,
      n_item = n_item,
      method = method,
      engine = engine,
      target = names(hits),
      max_error = unname(c(max_error, NA_real_)),
      nsim = nsim,
      n_valid_fit = sum(is.na(runs$fit_error)),
      n_extreme_item_reps = sum(runs$extreme_item),
      n_singular = if (all(is.na(runs$singular))) NA_integer_ else sum(runs$singular, na.rm = TRUE),
      n_evaluable = c(colSums(evaluable), sum(rowSums(evaluable) == ncol(losses))),
      rmse = c(rmse, NA_real_),
      n_met = hits,
      attainment = hits / nsim,
      mc_lo = intervals[1L, ],
      mc_hi = intervals[2L, ],
      decision = mapply(classify_mc_threshold, intervals[1L, ], intervals[2L, ], assurance),
      row.names = NULL,
      stringsAsFactors = FALSE
    )

    samples[[i]] <- data.frame(
      n_subj = ns,
      n_item = n_item,
      method = method,
      replication = seq_len(nsim),
      losses,
      extreme_item = runs$extreme_item,
      singular = runs$singular,
      fit_error = runs$fit_error,
      stringsAsFactors = FALSE
    )

    b_mat <- by_pt_b[[i]]
    err_mat <- sweep(b_mat, 2L, b_true, "-")
    item_summaries[[i]] <- data.frame(
      n_subj = ns,
      item_id = seq_len(n_item),
      n_scheduled = exposure[[i]]$n_scheduled,
      component = exposure[[i]]$component,
      true_difficulty = b_true,
      n_finite = colSums(is.finite(err_mat)),
      n_zero_correct = colSums(by_pt_z[[i]]),
      n_all_correct = colSums(by_pt_a[[i]]),
      bias = vapply(seq_len(n_item), function(j) {
        ej <- err_mat[, j]; if (any(is.finite(ej))) mean(ej[is.finite(ej)]) else NA_real_
      }, numeric(1L)),
      rmse = vapply(seq_len(n_item), function(j) {
        ej <- err_mat[, j]; if (any(is.finite(ej))) sqrt(mean(ej[is.finite(ej)]^2)) else NA_real_
      }, numeric(1L)),
      mae = vapply(seq_len(n_item), function(j) {
        ej <- err_mat[, j]; if (any(is.finite(ej))) mean(abs(ej[is.finite(ej)])) else NA_real_
      }, numeric(1L)),
      stringsAsFactors = FALSE
    )
  }

  summary_df <- do.call(rbind, summaries)
  joint_pass <- summary_df$n_subj[summary_df$target == "joint" & summary_df$decision == "pass"]
  rec_n <- if (length(joint_pass)) min(joint_pass) else NA_integer_

  structure(list(
    summary = summary_df,
    recommended_n_subj = rec_n,
    item_summary = do.call(rbind, item_summaries),
    replications = do.call(rbind, samples),
    item_difficulties = b_true,
    person_sd = person_sd,
    method = method,
    engine = engine,
    max_error = max_error,
    assurance = assurance,
    seed = seed,
    mc_conf_level = mc_conf_level
  ), class = "glmmplanr_rasch")
}

simulate_fixed_bank_trials <- function(des, schedule, b_true, person_sd, lightweight = FALSE) {
  ns <- schedule$n_subj
  theta <- stats::rnorm(ns, mean = 0.0, sd = person_sd)
  tr <- schedule$trials
  sid <- tr$subj_id
  iid <- tr$item_id
  eta <- theta[sid] - b_true[iid]
  y <- stats::rbinom(length(eta), size = 1L, prob = stats::plogis(eta))
  keep <- observed_trials(des, sid, ns)
  if (isTRUE(lightweight)) {
    trials_out <- if (!is.null(keep)) {
      list(subj_id = sid[keep], item_id = iid[keep], y = y[keep], has_keep = TRUE)
    } else {
      list(subj_id = sid, item_id = iid, y = y, has_keep = FALSE)
    }
    return(list(trials = trials_out, theta = theta, n_subj = ns, n_item = length(b_true)))
  }
  tr$y <- y
  if (!is.null(keep)) {
    tr <- tr[keep, , drop = FALSE]
    rownames(tr) <- NULL
  }
  list(trials = tr, theta = theta, n_subj = ns, n_item = length(b_true))
}

eval_fixed_bank_replication <- function(sim, b_true, person_sd, method, engine, control = NULL, rasch_ctx = NULL) {
  tr <- sim$trials
  ns <- sim$n_subj
  ni <- sim$n_item
  metrics <- data.frame(
    ability = NA_real_, difficulty = NA_real_, difficulty_max_abs_error = NA_real_,
    person_sd = NA_real_, person_variance = NA_real_,
    extreme_item = FALSE, singular = NA, fit_error = NA_character_,
    stringsAsFactors = FALSE
  )

  y_vec <- tr$y
  counts <- unit_response_counts(tr$item_id, y_vec == 1L, ni)
  n_obs_item <- counts$observed
  item_zero <- counts$all_incorrect
  item_all1 <- counts$all_correct
  metrics$extreme_item <- any(counts$unobserved | item_zero | item_all1)
  if (length(y_vec) < 2L || !any(y_vec == 0L) || !any(y_vec == 1L)) {
    metrics$fit_error <- "Fewer than two observed response values"
    return(list(
      metrics = metrics,
      b_est = rep(NA_real_, ni),
      item_zero = item_zero,
      item_all1 = item_all1
    ))
  }

  if (method == "random_item") {
    fit <- fit_rasch_random_intercept_model(
      tr, n_subj = ns, n_item = ni, engine = engine, control = control, rasch_ctx = rasch_ctx
    )
    if (!is.na(fit$error)) {
      metrics$fit_error <- fit$error
      return(list(metrics = metrics, b_est = rep(NA_real_, ni), item_zero = item_zero, item_all1 = item_all1))
    }
    metrics$singular <- fit$singular
    ab_est <- fit$ability
    b_est <- -fit$intercept - fit$item
    v_s <- fit$person_var

    b_err <- b_est - b_true
    metrics$ability <- sqrt(mean((ab_est - sim$theta)^2))
    metrics$difficulty <- sqrt(mean(b_err^2))
    metrics$difficulty_max_abs_error <- max(abs(b_err))
    metrics$person_variance <- abs(v_s - person_sd^2)
    metrics$person_sd <- abs(sqrt(v_s) - person_sd)
    return(list(metrics = metrics, b_est = b_est, item_zero = item_zero, item_all1 = item_all1))
  }

  # Fixed-item Marginal Maximum Likelihood (`method == "fixed_item"`)
  estimable <- (n_obs_item > 0L) & !item_zero & !item_all1
  if (!any(estimable)) {
    metrics$fit_error <- "No item has both response categories"
    return(list(metrics = metrics, b_est = rep(NA_real_, ni), item_zero = item_zero, item_all1 = item_all1))
  }

  b_est <- rep(NA_real_, ni)
  if (engine == "tam") {
    item_names <- sprintf("i%04d", seq_len(ni))
    resp <- matrix(NA_integer_, nrow = ns, ncol = ni, dimnames = list(NULL, item_names))
    resp[cbind(tr$subj_id, tr$item_id)] <- as.integer(tr$y)
    obs_s <- which(rowSums(!is.na(resp[, estimable, drop = FALSE])) > 0L)
    if (length(obs_s) < 2L) {
      metrics$fit_error <- "Fewer than two observed participants"
      return(list(metrics = metrics, b_est = b_est, item_zero = item_zero, item_all1 = item_all1))
    }
    tam_fit <- tryCatch({
      fit <- TAM::tam.mml(
        resp[obs_s, estimable, drop = FALSE],
        constraint = "cases",
        verbose = FALSE
      )
      if (fit$iter >= fit$control$maxiter) stop("TAM reached its iteration limit")
      fit
    }, warning = identity, error = identity)
    if (inherits(tam_fit, "condition")) {
      metrics$fit_error <- conditionMessage(tam_fit)
      return(list(metrics = metrics, b_est = b_est, item_zero = item_zero, item_all1 = item_all1))
    }
    b_est[match(rownames(tam_fit$xsi), item_names)] <- tam_fit$xsi$xsi
    ab_est <- rep(NA_real_, ns)
    ab_est[obs_s] <- tam_fit$person$EAP
    v_s <- as.numeric(tam_fit$variance[1L, 1L])
    metrics$singular <- v_s < 1e-8
    if (!all(estimable)) {
      metrics$fit_error <- NA_character_ # Fit succeeded on estimable items; constant items yield NA in b_est
    }
  } else {
    mml_fit <- fit_rasch_fixed_mml(tr, n_subj = ns, n_item = ni, estimable = estimable)
    if (!is.na(mml_fit$error)) {
      metrics$fit_error <- mml_fit$error
      return(list(metrics = metrics, b_est = b_est, item_zero = item_zero, item_all1 = item_all1))
    }
    b_est <- mml_fit$b_est
    ab_est <- mml_fit$ability
    v_s <- mml_fit$person_var
    metrics$singular <- mml_fit$singular
  }

  b_err <- b_est - b_true
  metrics$ability <- sqrt(mean((ab_est - sim$theta)^2))
  metrics$difficulty <- sqrt(mean(b_err^2))
  metrics$difficulty_max_abs_error <- max(abs(b_err))
  metrics$person_variance <- abs(v_s - person_sd^2)
  metrics$person_sd <- abs(sqrt(v_s) - person_sd)
  list(metrics = metrics, b_est = b_est, item_zero = item_zero, item_all1 = item_all1)
}

fit_rasch_fixed_mml <- function(trials, n_subj, n_item, estimable, max_iter = 200L, tol = 1e-6) {
  q1 <- gh_quad_1d(21L)
  cpp_fit_rasch_fixed_mml(
    subj_id = as.integer(trials$subj_id),
    item_id = as.integer(trials$item_id),
    y = as.integer(trials$y),
    n_subj = as.integer(n_subj),
    n_item = as.integer(n_item),
    estimable = as.logical(estimable),
    gh_nodes = q1$nodes,
    gh_weights = q1$weights,
    max_iter = as.integer(max_iter),
    tol = as.numeric(tol)
  )
}
