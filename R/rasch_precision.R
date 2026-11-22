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
#' @return A list with `summary` (per-target and joint attainment, counts, Wilson
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
#' plan$summary
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
    rows <- lapply(tsk$streams, function(stream) {
      with_rng_stream(stream, {
        rasch_recovery_errors(
          draw_dgp_trials(dgp, schedule, keep_person_item_effects = TRUE),
          dgp,
          control,
          engine = engine
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
  list(
    summary = summary,
    pareto = extract_2d_pareto(summary[passing, , drop = FALSE]),
    replications = do.call(rbind, samples),
    max_error = max_error,
    assurance = assurance,
    seed = seed,
    mc_conf_level = mc_conf_level
  )
}

rasch_recovery_errors <- function(draw, dgp, control = NULL, engine = "lme4") {
  data <- draw$trials
  result <- data.frame(
    ability = NA_real_, difficulty = NA_real_, difficulty_max_abs_error = NA_real_,
    person_sd = NA_real_, item_sd = NA_real_, person_variance = NA_real_, item_variance = NA_real_,
    singular = NA, fit_error = NA_character_
  )
  if (length(unique(data$y)) < 2L) {
    result$fit_error <- "Fewer than two observed response values"
    return(result)
  }

  n_subj <- nrow(draw$subj)
  n_item <- nrow(draw$item)

  if (identical(engine, "cpp_laplace")) {
    n_obs <- nrow(data)
    ones <- matrix(1.0, nrow = n_obs, ncol = 1L)
    fit <- tryCatch(
      cpp_fit_crossed_glmm(
        y = as.integer(data$y),
        X = ones,
        Z_subj = ones,
        subj_idx = as.integer(data$subj_id) - 1L,
        n_subj = n_subj,
        Z_item = ones,
        item_idx = as.integer(data$item_id) - 1L,
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
      ),
      error = identity
    )
    if (inherits(fit, "error")) {
      result$fit_error <- conditionMessage(fit)
      return(result)
    }
    if (!isTRUE(fit$converged)) {
      result$fit_error <- paste0("C++ Laplace status: ", fit$status)
      return(result)
    }
    result$singular <- isTRUE(fit$singular)
    ability <- as.numeric(fit$u_subj[, 1L])
    item <- as.numeric(fit$u_item[, 1L])
    ability[tabulate(data$subj_id, nbins = n_subj) == 0L] <- NA_real_
    item[tabulate(data$item_id, nbins = n_item) == 0L] <- NA_real_

    result$ability <- sqrt(mean((ability - draw$subj[, 1L])^2))
    difficulty <- -as.numeric(fit$beta[1L]) - item
    truth <- -unname(dgp$beta[1L]) - draw$item[, 1L]
    difficulty_error <- difficulty - truth
    result$difficulty <- sqrt(mean(difficulty_error^2))
    result$difficulty_max_abs_error <- max(abs(difficulty_error))

    v_s <- as.numeric(fit$Sigma_subj[1L, 1L])
    v_i <- as.numeric(fit$Sigma_item[1L, 1L])
    if (v_s <= 1e-8) v_s <- 0.0
    if (v_i <= 1e-8) v_i <- 0.0
    result$person_variance <- abs(v_s - dgp$tau_subj[1L]^2)
    result$item_variance <- abs(v_i - dgp$tau_item[1L]^2)
    result$person_sd <- abs(sqrt(v_s) - dgp$tau_subj[1L])
    result$item_sd <- abs(sqrt(v_i) - dgp$tau_item[1L])
    return(result)
  }

  if (is.null(control)) control <- lme4::glmerControl(optimizer = "bobyqa")
  warnings <- character()
  fit <- tryCatch(withCallingHandlers(
    suppressMessages(lme4::glmer(y ~ 1 + (1 | subj) + (1 | item), data = data,
                                 family = stats::binomial(), control = control)),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = identity)
  if (inherits(fit, "error")) {
    result$fit_error <- conditionMessage(fit)
    return(result)
  }
  conv <- fit@optinfo$conv
  diagnostics <- c(warnings, conv$lme4$messages)
  diagnostics <- diagnostics[!grepl("boundary \\(singular\\) fit", diagnostics)]
  if (any(c(conv$opt, conv$lme4$code) != 0) || length(diagnostics)) {
    result$fit_error <- paste(c("Fit diagnostic", unique(diagnostics)), collapse = " / ")
    return(result)
  }
  result$singular <- lme4::isSingular(fit)
  modes <- lme4::ranef(fit, condVar = FALSE)
  # Preserve planned integer IDs when glmer drops an unobserved level.
  ability <- modes$subj[match(seq_len(n_subj), rownames(modes$subj)), 1L]
  item <- modes$item[match(seq_len(n_item), rownames(modes$item)), 1L]
  result$ability <- sqrt(mean((ability - draw$subj[, 1L])^2))
  difficulty <- -unname(lme4::fixef(fit)[1L]) - item
  truth <- -unname(dgp$beta[1L]) - draw$item[, 1L]
  difficulty_error <- difficulty - truth
  result$difficulty <- sqrt(mean(difficulty_error^2))
  result$difficulty_max_abs_error <- max(abs(difficulty_error))
  variance <- lme4::VarCorr(fit)
  result$person_variance <- abs(variance$subj[1L, 1L] - dgp$tau_subj[1L]^2)
  result$item_variance <- abs(variance$item[1L, 1L] - dgp$tau_item[1L]^2)
  result$person_sd <- abs(sqrt(variance$subj[1L, 1L]) - dgp$tau_subj[1L])
  result$item_sd <- abs(sqrt(variance$item[1L, 1L]) - dgp$tau_item[1L])
  result
}

#' Plan Examinee Sample Sizes for Calibrating a Fixed Rasch Item Bank
#'
#' Evaluates examinee sample sizes (`n_subj`) for calibrating a fixed,
#' pre-specified vector of item difficulties `item_difficulties` (held constant
#' across replications, whereas examinee abilities \eqn{\theta_i \sim N(0, \sigma_p^2)}
#' are drawn anew in every replication). Supports both fixed-item Marginal
#' Maximum Likelihood (`method = "fixed_item"`) and random-item Empirical Bayes
#' shrinkage (`method = "random_item"`), accounting explicitly for extreme
#' all-correct / all-incorrect item failures in the assurance denominator.
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
#'     yielding finite Empirical Bayes shrunken item difficulties even when an
#'     administered extreme item has zero correct or zero incorrect responses.
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
#' @return A list with:
#'   * `summary`: Per-target and `joint` attainment, `n_extreme_item_reps`,
#'     `n_evaluable`, `rmse`, Wilson bounds (`mc_lo`, `mc_hi`), and `decision`.
#'   * `recommended_n_subj`: Smallest candidate `n_subj` achieving `PASS` on
#'     `target == "joint"` (`NA_integer_` if none pass).
#'   * `item_summary`: Per-item bias, RMSE, MAE, finite-estimate count, and
#'     all-zero / all-one counts across replications at each `n_subj`.
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
#' fb$summary
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
    schedule <- generate_design_schedule(des, ns, n_item)
    nb <- length(tsk$rep_idx)
    rows <- vector("list", nb)
    b_est_mat <- matrix(NA_real_, nrow = nb, ncol = n_item)
    zero_mat <- matrix(FALSE, nrow = nb, ncol = n_item)
    all1_mat <- matrix(FALSE, nrow = nb, ncol = n_item)

    for (k in seq_len(nb)) {
      out_k <- with_rng_stream(tsk$streams[[k]], {
        sim <- simulate_fixed_bank_trials(des, schedule, b_true = b_true, person_sd = person_sd)
        eval_fixed_bank_replication(
          sim = sim, b_true = b_true, person_sd = person_sd,
          method = method, engine = engine, control = control
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
      b_est = b_est_mat,
      item_zero = zero_mat,
      item_all1 = all1_mat
    )
  }

  chunk_out <- run_parallel_tasks(
    tasks, worker_fn, workers = workers, parallel_backend = parallel_backend,
    stage_label = sprintf("Fixed-Bank Rasch (%s)", method), verbose = verbose
  )

  by_pt_rows <- vector("list", n_pts)
  by_pt_b <- lapply(seq_len(n_pts), function(i) matrix(NA_real_, nrow = nsim, ncol = n_item))
  by_pt_z <- lapply(seq_len(n_pts), function(i) matrix(FALSE, nrow = nsim, ncol = n_item))
  by_pt_a <- lapply(seq_len(n_pts), function(i) matrix(FALSE, nrow = nsim, ncol = n_item))

  for (cr in chunk_out) {
    p_i <- cr$pt_idx
    if (is.null(by_pt_rows[[p_i]])) by_pt_rows[[p_i]] <- vector("list", nsim)
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

  list(
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
  )
}

simulate_fixed_bank_trials <- function(des, schedule, b_true, person_sd) {
  ns <- schedule$n_subj
  theta <- stats::rnorm(ns, mean = 0.0, sd = person_sd)
  tr <- schedule$trials
  eta <- theta[tr$subj_id] - b_true[tr$item_id]
  tr$y <- stats::rbinom(nrow(tr), size = 1L, prob = stats::plogis(eta))
  keep <- observed_trials(des, tr$subj_id, ns)
  if (!is.null(keep)) {
    tr <- tr[keep, , drop = FALSE]
    rownames(tr) <- NULL
  }
  list(trials = tr, theta = theta, n_subj = ns, n_item = length(b_true))
}

eval_fixed_bank_replication <- function(sim, b_true, person_sd, method, engine, control = NULL) {
  tr <- sim$trials
  ns <- sim$n_subj
  ni <- sim$n_item
  metrics <- data.frame(
    ability = NA_real_, difficulty = NA_real_, difficulty_max_abs_error = NA_real_,
    person_sd = NA_real_, person_variance = NA_real_,
    extreme_item = FALSE, singular = NA, fit_error = NA_character_,
    stringsAsFactors = FALSE
  )

  if (nrow(tr) == 0L || length(unique(tr$y)) < 2L) {
    metrics$fit_error <- "Fewer than two observed response values"
    return(list(
      metrics = metrics,
      b_est = rep(NA_real_, ni),
      item_zero = rep(FALSE, ni),
      item_all1 = rep(FALSE, ni)
    ))
  }

  n_obs_item <- tabulate(tr$item_id, nbins = ni)
  s_obs_item <- as.numeric(rowsum(tr$y, factor(tr$item_id, levels = seq_len(ni)), reorder = FALSE))
  s_obs_item[is.na(s_obs_item)] <- 0.0
  item_zero <- (n_obs_item > 0L) & (s_obs_item == 0.0)
  item_all1 <- (n_obs_item > 0L) & (s_obs_item == n_obs_item)
  metrics$extreme_item <- any(n_obs_item == 0L | item_zero | item_all1)

  if (method == "random_item") {
    # Construct synthetic draw compatible with rasch_recovery_errors
    fake_dgp <- list(
      beta = c(`(Intercept)` = -mean(b_true)),
      tau_subj = person_sd,
      tau_item = stats::sd(b_true)
    )
    fake_draw <- list(
      trials = tr,
      subj = matrix(sim$theta, ncol = 1L),
      item = matrix(mean(b_true) - b_true, ncol = 1L)
    )
    # Extract difficulty estimates directly so item_summary is populated too
    if (engine == "cpp_laplace") {
      n_obs <- nrow(tr)
      ones <- matrix(1.0, nrow = n_obs, ncol = 1L)
      fit <- tryCatch(
        cpp_fit_crossed_glmm(
          y = as.integer(tr$y), X = ones,
          Z_subj = ones, subj_idx = as.integer(tr$subj_id) - 1L, n_subj = ns,
          Z_item = ones, item_idx = as.integer(tr$item_id) - 1L, n_item = ni,
          Z_testlet = matrix(0.0, nrow = n_obs, ncol = 0L),
          testlet_idx = integer(n_obs), n_testlet = 0L,
          chance = 0.0, corr_mode = 0L, regularized = FALSE,
          lkj_eta = 2.0, prior_tau_sd = 0.0, prior_beta_sd = 0.0,
          max_em_iter = 160L, max_pirls_iter = 25L, tol = 1e-4,
          fix_cov = FALSE, compute_joint_vcov = FALSE, return_ranef = TRUE
        ),
        error = identity
      )
      if (inherits(fit, "error") || !isTRUE(fit$converged)) {
        metrics$fit_error <- if (inherits(fit, "error")) conditionMessage(fit) else paste0("C++ Laplace status: ", fit$status)
        return(list(metrics = metrics, b_est = rep(NA_real_, ni), item_zero = item_zero, item_all1 = item_all1))
      }
      metrics$singular <- isTRUE(fit$singular)
      ab_est <- as.numeric(fit$u_subj[, 1L])
      it_est <- as.numeric(fit$u_item[, 1L])
      ab_est[tabulate(tr$subj_id, nbins = ns) == 0L] <- NA_real_
      it_est[n_obs_item == 0L] <- NA_real_
      b_est <- -as.numeric(fit$beta[1L]) - it_est
      v_s <- as.numeric(fit$Sigma_subj[1L, 1L])
      if (v_s <= 1e-8) v_s <- 0.0
    } else {
      if (is.null(control)) control <- lme4::glmerControl(optimizer = "bobyqa")
      warnings <- character()
      fit <- tryCatch(withCallingHandlers(
        suppressMessages(lme4::glmer(y ~ 1 + (1 | subj) + (1 | item), data = tr,
                                     family = stats::binomial(), control = control)),
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }), error = identity)
      if (inherits(fit, "error")) {
        metrics$fit_error <- conditionMessage(fit)
        return(list(metrics = metrics, b_est = rep(NA_real_, ni), item_zero = item_zero, item_all1 = item_all1))
      }
      conv <- fit@optinfo$conv
      diagnostics <- c(warnings, conv$lme4$messages)
      diagnostics <- diagnostics[!grepl("boundary \\(singular\\) fit", diagnostics)]
      if (any(c(conv$opt, conv$lme4$code) != 0) || length(diagnostics)) {
        metrics$fit_error <- paste(c("Fit diagnostic", unique(diagnostics)), collapse = " / ")
        return(list(metrics = metrics, b_est = rep(NA_real_, ni), item_zero = item_zero, item_all1 = item_all1))
      }
      metrics$singular <- lme4::isSingular(fit)
      modes <- lme4::ranef(fit, condVar = FALSE)
      ab_est <- modes$subj[match(seq_len(ns), rownames(modes$subj)), 1L]
      it_est <- modes$item[match(seq_len(ni), rownames(modes$item)), 1L]
      b_est <- -unname(lme4::fixef(fit)[1L]) - it_est
      v_s <- as.numeric(lme4::VarCorr(fit)$subj[1L, 1L])
    }

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
    obs_s <- which(rowSums(!is.na(resp)) > 0L)
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
  M <- matrix(0.0, nrow = n_subj, ncol = n_item)
  Y <- matrix(0.0, nrow = n_subj, ncol = n_item)
  idx <- cbind(trials$subj_id, trials$item_id)
  M[idx] <- 1.0
  Y[idx] <- as.numeric(trials$y)

  obs_s <- which(rowSums(M) > 0)
  if (length(obs_s) < 2L) {
    return(list(error = "Fewer than two observed participants"))
  }

  est_j <- which(estimable)
  M_sub <- M[obs_s, est_j, drop = FALSE]
  Y_sub <- Y[obs_s, est_j, drop = FALSE]
  n_sub <- colSums(M_sub)
  s_sub <- colSums(Y_sub)

  q1 <- gh_quad_1d(21L)
  z <- q1$nodes
  log_w <- log(q1$weights)

  p0 <- pmin(pmax(s_sub / n_sub, 0.02), 0.98)
  b <- -stats::qlogis(p0)
  sigma <- 0.80
  converged <- FALSE

  n_obs_s <- nrow(M_sub)
  all_obs <- all(M_sub == 1.0)
  log_w_mat <- matrix(log_w, nrow = n_obs_s, ncol = length(log_w), byrow = TRUE)

  for (it in seq_len(max_iter)) {
    b_old <- b
    sig_old <- sigma

    theta_q <- sigma * z
    eta_jq <- outer(-b, theta_q, "+")
    sp_jq <- ifelse(eta_jq > 0, eta_jq + log1p(exp(-eta_jq)), log1p(exp(eta_jq)))
    LL <- if (all_obs) {
      sweep(Y_sub %*% eta_jq, 2L, colSums(sp_jq), "-")
    } else {
      Y_sub %*% eta_jq - M_sub %*% sp_jq
    }
    LP <- LL + log_w_mat
    mx <- LP[cbind(seq_len(n_obs_s), max.col(LP, ties.method = "first"))]
    W_post <- exp(LP - mx)
    R_iq <- W_post / rowSums(W_post)

    # M-step for item difficulties b_j given posterior weights R_iq
    N_jq <- if (all_obs) {
      matrix(colSums(R_iq), nrow = length(b), ncol = length(z), byrow = TRUE)
    } else {
      crossprod(M_sub, R_iq)
    }
    for (nw in 1:4) {
      P_jq <- 1.0 / (1.0 + exp(-outer(-b, theta_q, "+")))
      grad_b <- rowSums(N_jq * P_jq) - s_sub
      hess_b <- pmax(rowSums(N_jq * P_jq * (1.0 - P_jq)), 1e-8)
      step_b <- pmin(pmax(grad_b / hess_b, -2.0), 2.0)
      b <- b + step_b
      if (max(abs(step_b)) < 1e-8) break
    }

    # Case-centering identification (`constraint = "cases"`): E_emp[EAP] = 0
    eap_i <- as.numeric(R_iq %*% theta_q)
    mu_case <- mean(eap_i)
    b <- b - mu_case
    var_new <- sum(colSums(R_iq) * (theta_q - mu_case)^2) / n_obs_s
    sigma <- sqrt(max(var_new, 1e-10))

    diff_max <- max(max(abs(b - b_old)), abs(sigma - sig_old))
    if (diff_max < tol && it >= 3L) {
      converged <- TRUE
      break
    }
  }

  if (!converged) {
    return(list(error = "MML-EM reached iteration limit"))
  }

  theta_q <- sigma * z
  eta_jq <- outer(-b, theta_q, "+")
  sp_jq <- ifelse(eta_jq > 0, eta_jq + log1p(exp(-eta_jq)), log1p(exp(eta_jq)))
  LL <- if (all_obs) {
    sweep(Y_sub %*% eta_jq, 2L, colSums(sp_jq), "-")
  } else {
    Y_sub %*% eta_jq - M_sub %*% sp_jq
  }
  LP <- LL + log_w_mat
  mx <- LP[cbind(seq_len(n_obs_s), max.col(LP, ties.method = "first"))]
  W_post <- exp(LP - mx)
  R_iq <- W_post / rowSums(W_post)
  eap_obs <- as.numeric(R_iq %*% theta_q)
  mu_case <- mean(eap_obs)
  b <- b - mu_case
  eap_obs <- eap_obs - mu_case

  b_full <- rep(NA_real_, n_item)
  b_full[est_j] <- b
  ab_full <- rep(NA_real_, n_subj)
  ab_full[obs_s] <- eap_obs
  v_s <- if (sigma <= 1e-4) 0.0 else sigma^2

  list(
    error = NA_character_,
    b_est = b_full,
    ability = ab_full,
    person_var = v_s,
    singular = sigma <= 1e-4
  )
}
