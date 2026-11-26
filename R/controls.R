#' Control Settings for Calibration, Target Specification, and Optimization
#'
#' Constructor functions bundling advanced numerical, testlet, prior, Stage-1
#' screening, and Monte Carlo control parameters for [calibrate_dgp()],
#' [specify_target()], and [optimize_sample_size()].
#'
#' @param testlet_sd Optional numeric scalar or vector of conditional logit SDs
#'   for passages/testlets when `design$has_testlet` is `TRUE` (defaults to
#'   `0.5 * item_intercept_sd` if neither `testlet_sd` nor `testlet_prob_range`
#'   is supplied).
#' @param testlet_prob_range Optional numeric vector `c(low, high)` for passage
#'   / testlet 95% accuracy range when `design$has_testlet` is `TRUE`.
#' @param participant_testlet_sd Conditional logit SD of an independent
#'   participant-by-passage intercept in the generating model. Use positive
#'   values to stress-test a passage-only analysis against residual local
#'   dependence (default `0`).
#' @param n_quad Integer number of 1D Gauss-Hermite quadrature nodes (initially
#'   `21L`, refined automatically up to `256L` if needed).
#' @param check_nsim Integer number of Monte Carlo draws per cell used for the
#'   prior-predictive sanity check (default `4000L`).
#' @param seed Integer seed for the prior-predictive check (default `1L`). Use
#'   `NULL` to draw a master seed from the caller's random stream.
#'
#' @return A validated control list of class `"glmmplanr_calibrate_control"`,
#'   `"glmmplanr_target_control"`, or `"glmmplanr_optimize_control"`.
#' @examples
#' calibrate_control(testlet_sd = 0.45, participant_testlet_sd = 0.30)
#' target_control(failure_policy = "effective", prob_boot = 99L)
#' optimize_control(nsim_init = 40L, nsim_sensitivity = 0L)
#' @export
calibrate_control <- function(
  testlet_sd = NULL,
  testlet_prob_range = NULL,
  participant_testlet_sd = 0,
  n_quad = 21L,
  check_nsim = 4000L,
  seed = 1L
) {
  if (!is.numeric(participant_testlet_sd) || length(participant_testlet_sd) != 1L ||
      !is.finite(participant_testlet_sd) || participant_testlet_sd < 0) {
    cli::cli_abort("{.arg participant_testlet_sd} must be a non-negative numeric scalar.")
  }
  if (!is.numeric(n_quad) || length(n_quad) != 1L || !is.finite(n_quad) ||
      n_quad < 1 || n_quad != floor(n_quad)) {
    cli::cli_abort("{.arg n_quad} must be a positive integer.")
  }
  if (!is.numeric(check_nsim) || length(check_nsim) != 1L || !is.finite(check_nsim) ||
      check_nsim < 1 || check_nsim != floor(check_nsim)) {
    cli::cli_abort("{.arg check_nsim} must be a positive integer.")
  }
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      abs(seed) > .Machine$integer.max || seed != trunc(seed))) {
    cli::cli_abort("{.arg seed} must be NULL or a single finite integer.")
  }
  structure(
    list(
      testlet_sd = testlet_sd,
      testlet_prob_range = testlet_prob_range,
      participant_testlet_sd = as.numeric(participant_testlet_sd),
      n_quad = as.integer(n_quad),
      check_nsim = as.integer(check_nsim),
      seed = if (is.null(seed)) NULL else as.integer(seed)
    ),
    class = "glmmplanr_calibrate_control"
  )
}

#' @rdname calibrate_control
#' @param prob_estimand Estimand used for probability-difference (`"pp"` /
#'   `"prob"`) calculations: `"marginal"` (default, population-standardized
#'   G-computation over random effects) or `"conditional"` (at zero random
#'   effects).
#' @param contrast_weights Optional custom numeric vector of cell contrast
#'   weights of length `dgp$design$n_cells` (defaults to
#'   [cell_contrast_weights()]).
#' @param ci_assurance Target probability \eqn{\Pr(W \le w_{\max})} for AIPE
#'   interval width assurance (default `0.80`).
#' @param rope Numeric vector `c(low, high)` defining the Region of Practical
#'   Equivalence in units matching `scale` (defaults to `c(-0.18, 0.18)` on
#'   `"logit"`, `c(-3, 3)` on `"pp"`, or `c(-0.03, 0.03)` on `"prob"`).
#' @param failure_policy How non-converged or failed replications enter sample
#'   size pass/fail decisions: `"effective"` (default; failed fits count as `0`)
#'   or `"conditional"` (evaluates power among valid fits while requiring
#'   `valid_rate >= min_valid_rate`).
#' @param min_valid_rate Minimum acceptable proportion of valid converged fits
#'   when `failure_policy = "conditional"` (default `0.90`).
#' @param mc_conf_level Confidence level for Wilson Monte Carlo uncertainty
#'   intervals classifying borderline candidates as `"pass"`, `"indeterminate"`,
#'   or `"fail"` (default `0.95`; Morris et al., 2019).
#' @param stage1_engine Fast estimation backend used during Stage-1 surrogate
#'   exploration before Stage-2 verification (`NULL` selects `"bayes_laplace"`
#'   for regularized/Bayesian primary engines and `"cpp_laplace"` otherwise).
#' @param stage1_re Random-effect covariance structure used during Stage-1
#'   surrogate exploration: `"match"` (default) or `"uncorrelated"`.
#' @param on_singular Policy when maximum-likelihood fits hit a singular
#'   random-effects boundary: `"flag_and_compare"` (default), `"flag_only"`,
#'   `"reduce_corr"`, or `"switch_bayes"`.
#' @param analysis_link Link function used in the analysis model: `"match_dgp"`
#'   (default) or `"logit"`.
#' @param lkj_eta Shape parameter \eqn{\eta} for the LKJ correlation prior in
#'   Stan and correlation shrinkage in `"bayes_laplace"` (default `2.0`).
#' @param prior_beta_sd Weakly informative Normal prior SD on non-intercept
#'   fixed effects in Bayesian engines (default `2.5`).
#' @param prior_tau_sd Scale parameter for the Half-Student-t(\eqn{\nu = 3})
#'   prior on random-effect SDs in `"stan_hmc"` and `"bayes_laplace"` (default
#'   `1.5`).
#' @param prob_boot Number of parametric bootstrap refits per marginal
#'   probability-scale replication when `prob_se_method = "bootstrap"` (default
#'   `99L`).
#' @export
target_control <- function(
  prob_estimand = c("marginal", "conditional"),
  contrast_weights = NULL,
  ci_assurance = 0.80,
  rope = NULL,
  failure_policy = c("effective", "conditional"),
  min_valid_rate = 0.90,
  mc_conf_level = 0.95,
  stage1_engine = NULL,
  stage1_re = c("match", "uncorrelated"),
  on_singular = c("flag_and_compare", "flag_only", "reduce_corr", "switch_bayes"),
  analysis_link = c("match_dgp", "logit"),
  lkj_eta = 2.0,
  prior_beta_sd = 2.5,
  prior_tau_sd = 1.5,
  prob_boot = 99L
) {
  prob_estimand <- match.arg(prob_estimand)
  failure_policy <- match.arg(failure_policy)
  stage1_re <- match.arg(stage1_re)
  on_singular <- match.arg(on_singular)
  analysis_link <- match.arg(analysis_link)

  if (!is.numeric(prob_boot) || length(prob_boot) != 1L || !is.finite(prob_boot) ||
      prob_boot < 2 || prob_boot > .Machine$integer.max || prob_boot != floor(prob_boot)) {
    cli::cli_abort("{.arg prob_boot} must be an integer of at least 2.")
  }
  for (scale_name in c("lkj_eta", "prior_beta_sd", "prior_tau_sd")) {
    value <- get(scale_name, inherits = FALSE)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value <= 0) {
      cli::cli_abort("{.arg {scale_name}} must be a finite positive scalar.")
    }
  }
  if (!is.numeric(ci_assurance) || length(ci_assurance) != 1L || is.na(ci_assurance) ||
      ci_assurance <= 0 || ci_assurance >= 1) {
    cli::cli_abort("{.arg ci_assurance} must lie strictly in (0, 1).")
  }
  if (!is.null(rope) && (!is.numeric(rope) || length(rope) != 2L || anyNA(rope) || rope[1L] >= rope[2L])) {
    cli::cli_abort("{.arg rope} must be a numeric vector c(low, high) with low < high.")
  }
  if (!is.numeric(min_valid_rate) || length(min_valid_rate) != 1L || is.na(min_valid_rate) ||
      min_valid_rate <= 0 || min_valid_rate > 1) {
    cli::cli_abort("{.arg min_valid_rate} must lie in (0, 1].")
  }
  if (!is.numeric(mc_conf_level) || length(mc_conf_level) != 1L || is.na(mc_conf_level) ||
      mc_conf_level <= 0.5 || mc_conf_level >= 1) {
    cli::cli_abort("{.arg mc_conf_level} must lie strictly in (0.5, 1).")
  }
  if (!is.null(stage1_engine)) {
    stage1_engine <- match.arg(
      stage1_engine,
      c("bayes_laplace", "cpp_laplace", "lme4", "lme4_fast", "glmmTMB", "stan_hmc", "stan_laplace")
    )
  }

  structure(
    list(
      prob_estimand = prob_estimand,
      contrast_weights = contrast_weights,
      ci_assurance = as.numeric(ci_assurance),
      rope = rope,
      failure_policy = failure_policy,
      min_valid_rate = as.numeric(min_valid_rate),
      mc_conf_level = as.numeric(mc_conf_level),
      stage1_engine = stage1_engine,
      stage1_re = stage1_re,
      on_singular = on_singular,
      analysis_link = analysis_link,
      lkj_eta = as.numeric(lkj_eta),
      prior_beta_sd = as.numeric(prior_beta_sd),
      prior_tau_sd = as.numeric(prior_tau_sd),
      prob_boot = as.integer(prob_boot)
    ),
    class = "glmmplanr_target_control"
  )
}

#' @rdname calibrate_control
#' @param n_eval_init Number of space-filling initial lattice points evaluated
#'   in Stage 1 when `method = "surrogate_gp"` (default `16L`).
#' @param n_eval_adapt Number of adaptive contour-focused GP evaluations added
#'   after `n_eval_init` in Stage 1 (default `8L`).
#' @param nsim_init Number of Monte Carlo simulations per exploratory point in
#'   Stage 1 (default `80L`).
#' @param nsim_topup Pre-specified number of independent additional Monte Carlo
#'   replications run on Stage-2 `"indeterminate"` candidates (default `0L`).
#' @param topup_policy How top-up replications are evaluated when
#'   `nsim_topup > 0L`: `"independent"` (default, fresh batch with two-stage
#'   Bonferroni split) or `"pooled"` (exploratory, emits warning).
#' @param mc_adjustment Multiplicity adjustment applied to Stage-2 Wilson Monte
#'   Carlo intervals across verified candidates: `"none"` (default) or
#'   `"bonferroni"`. Adjusted decisions are used before step-up, ceiling and
#'   top-up checks; adding candidates also updates earlier decisions. In
#'   multi-scenario planning the family comprises candidate-by-scenario pairs.
#' @param n_verify Number of Pareto-frontier candidate designs verified at
#'   `nsim_final` replications in Stage 2 (default `4L`).
#' @param early_stop Logical; when `TRUE` (default) and `nsim_final >= 200L`,
#'   applies group-sequential futility curtailment during Stage-2 verification.
#' @param nsim_sensitivity Number of Monte Carlo simulations per level used to
#'   build the 3-level random-slope sensitivity table at the recommended design
#'   (default `200L`; set to `0L` to skip).
#'   Skipped automatically when all participant, item and testlet slope SDs
#'   are zero or no random-slope columns exist.
#' @param parallel_backend Parallel execution engine: `"auto"` (default),
#'   `"mirai"`, `"multicore"`, or `"sequential"`.
#' @export
optimize_control <- function(
  n_eval_init = 16L,
  n_eval_adapt = 8L,
  nsim_init = 80L,
  nsim_topup = 0L,
  topup_policy = c("independent", "pooled"),
  mc_adjustment = c("none", "bonferroni"),
  n_verify = 4L,
  early_stop = TRUE,
  nsim_sensitivity = 200L,
  parallel_backend = c("auto", "mirai", "multicore", "sequential")
) {
  topup_policy <- match.arg(topup_policy)
  mc_adjustment <- match.arg(mc_adjustment)
  parallel_backend <- match.arg(parallel_backend)
  for (nm in c("n_eval_init", "n_eval_adapt", "nsim_init", "nsim_topup", "n_verify", "nsim_sensitivity")) {
    val <- get(nm, inherits = FALSE)
    if (!is.numeric(val) || length(val) != 1L || !is.finite(val) || val < 0 || val != floor(val)) {
      cli::cli_abort("{.arg {nm}} must be a non-negative integer.")
    }
  }
  if (n_verify < 1L) {
    cli::cli_abort("{.arg n_verify} must be an integer of at least 1.")
  }
  if (nsim_init < 1L) {
    cli::cli_abort("{.arg nsim_init} must be an integer of at least 1.")
  }
  if (!rlang::is_bool(early_stop)) {
    cli::cli_abort("{.arg early_stop} must be TRUE or FALSE.")
  }

  structure(
    list(
      n_eval_init = as.integer(n_eval_init),
      n_eval_adapt = as.integer(n_eval_adapt),
      nsim_init = as.integer(nsim_init),
      nsim_topup = as.integer(nsim_topup),
      topup_policy = topup_policy,
      mc_adjustment = mc_adjustment,
      n_verify = as.integer(n_verify),
      early_stop = isTRUE(early_stop),
      nsim_sensitivity = as.integer(nsim_sensitivity),
      parallel_backend = parallel_backend
    ),
    class = "glmmplanr_optimize_control"
  )
}
