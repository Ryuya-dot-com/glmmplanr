#' Specify Target Effect, Scale, Evaluation Criteria, and Estimation Strategy
#'
#' Defines the inferential or precision objective (statistical power, AIPE
#' confidence/credible interval width, joint power + precision, or ROPE
#' exclusion) on either the conditional log-odds (`"logit"`) scale or the
#' substantive response-probability difference (`"pp"` percentage points /
#' `"prob"`) scale, along with the estimation backend, estimation-failure
#' accounting policy, and Monte Carlo uncertainty confidence level.
#'
#' For a population difference such as 60% versus 70%, choose `scale = "pp"`
#' (10 percentage points). The default targets a conditional logit coefficient.
#' For a quick marginal pp grid, set `prob_se_method = "joint_delta"`.
#' A surrogate search already screens cheaply in Stage 1. Final grid or
#' surrogate verification with the default bootstrap can require 99 refits
#' per simulated study. To check one selected pair with bootstrap intervals,
#' specify the final method and call [validate_design_point()] with a fresh seed.
#' [formula()] and [family()] give the planned model structure and analysis link.
#'
#' @param dgp A `"glmmplanr_dgp"` object produced by [calibrate_dgp()].
#' @param effect Character string naming the focal fixed effect (e.g.,
#'   `"exposure"`, `"spacing"`, or `"exposure:spacing"`). Both factor-level term
#'   labels and model-matrix column names (e.g., `"exposure1:spacing1"`) are
#'   accepted. For a design with no experimental factors, use `"(Intercept)"`.
#'   On the logit scale its null is a zero conditional intercept (probability
#'   0.5 for an ordinary logit at zero random effects). On the probability
#'   scale the one cell gives an absolute probability: use `criterion = "ci_width"`
#'   for its precision, and confirm planned-condition coverage and width
#'   attainment via [validate_design_point()] or [evaluate_design_point()]
#'   (no zero-difference null applies). Variance components and individual
#'   scores require a separate analysis.
#' @param scale Scale on which detection, precision (`ci_width`) and equivalence (`rope`)
#'   targets are expressed and evaluated:
#'   * `"logit"` (default): Conditional log-odds scale.
#'   * `"pp"` (or `"percentage_point"`): Response accuracy difference in
#'     **percentage points** (e.g., `ci_width = 10` means a full 95% CI width of
#'     at most 10 percentage points, such as `+5.0 pp [0.0, 10.0]`).
#'   * `"prob"`: Response accuracy difference on `[0, 1]` probability scale
#'     (e.g., `ci_width = 0.10` for 10 percentage points).
#' @param criterion Evaluation objective:
#'   * `"power"`: Statistical power \eqn{\ge} `power_min`.
#'   * `"ci_width"`: Precision / AIPE criterion requiring the `(1 - alpha)`
#'     interval width to be \eqn{\le} `ci_width` with probability \eqn{\ge}
#'     `ci_assurance` (set in [target_control()]).
#'   * `"both"`: Requires both statistical power \eqn{\ge} `power_min` and
#'     interval width criterion (`ci_width`), assessed separately.
#'   * `"joint"`: Probability of detecting the effect and meeting the width bound
#'     in the same replication must reach `power_min`.
#'   * `"rope"`: Probability that the `(1 - alpha)` interval excludes `rope`
#'     is \eqn{\ge} `power_min`.
#' @param power_min Minimum target power (or ROPE exclusion probability),
#'   default `0.80`.
#' @param ci_width Maximum target `(1 - alpha)` confidence/credible interval
#'   width \eqn{w_{\max}} in units matching `scale` (`logit` when
#'   `scale = "logit"`, percentage points `0..100` when `scale = "pp"`, or
#'   probability `0..1` when `scale = "prob"`).
#' @param alpha Two-sided significance level / `1 - conf_level` (default `0.05`).
#' @param engine Primary estimation backend: `"bayes_laplace"` (default,
#'   built-in regularized C++ Laplace engine supporting both ordinary binary
#'   logit and `chance > 0` links), `"cpp_laplace"`, `"lme4"`, `"lme4_fast"`,
#'   `"glmmTMB"`, `"stan_hmc"`, or legacy alias `"stan_laplace"` (deprecated;
#'   executes full HMC/NUTS via `"stan_hmc"`). Optional external backends
#'   (`"glmmTMB"`, `"stan_hmc"`) abort rather than silently falling back if
#'   unavailable or incompatible with `chance > 0`, ensuring numerical
#'   consistency. Unpenalized engines flag separation if all observed responses
#'   in any fixed-effect design cell are 0 or all are 1. For `"cpp_laplace"`
#'   and `"bayes_laplace"`, analytical fixed-effect SEs combine the
#'   Tierney-Kadane Laplace Hessian with a diagonal log-SD variance propagation
#'   term; correlated-slope designs (especially under $m$-AFC or sparsity)
#'   should be cross-checked against `"stan_hmc"` or parametric bootstrap SEs.
#' @param analysis_re Random-effect structure of the fitted analysis model:
#'   `"match_dgp"` (default; automatically fits `"intercept_only"` if `dgp` has
#'   no non-zero random slopes, `"uncorrelated"` if slopes have zero
#'   correlations, or `"maximal"` otherwise), `"maximal"`, `"uncorrelated"`, or
#'   `"intercept_only"`.
#' @param prob_se_method Standard-error method for marginal probability-scale
#'   contrasts in Stage-2 final verification:
#'   * `"bootstrap"` (default): Parametric bootstrap refits (`prob_boot` draws
#'     in [target_control()]) executed via a pre-allocated C++ batch kernel for
#'     built-in engines.
#'   * `"joint_delta"`: Analytical Gauss-Hermite delta method propagating joint
#'     uncertainty in fixed effects \eqn{\beta} and random-effect log-SDs
#'     \eqn{\theta = \log \tau} (also used automatically in Stage-1 screening).
#' @param control Control list created by [target_control()] holding advanced
#'   settings (`prob_estimand`, `contrast_weights`, `ci_assurance`, `rope`,
#'   `failure_policy`, `min_valid_rate`, `mc_conf_level`, `stage1_engine`,
#'   `stage1_re`, `on_singular`, `analysis_link`, `lkj_eta`, `prior_beta_sd`,
#'   `prior_tau_sd`, `prob_boot`).
#' @return An S3 object of class `"glmmplanr_spec"`.
#' @examples
#' des <- design_crossed(list(condition = c("a", "b")))
#' dgp <- calibrate_dgp(des, c(a = 0.4, b = 0.6), subj_sd = 0.3, item_sd = 0.4)
#' spec <- specify_target(dgp, "condition", scale = "pp", criterion = "ci_width", ci_width = 12)
#' spec$true_effect_target
#' @export
specify_target <- function(
  dgp,
  effect,
  scale = c("logit", "pp", "percentage_point", "prob"),
  criterion = c("power", "ci_width", "both", "rope", "joint"),
  power_min = 0.80,
  ci_width = NULL,
  alpha = 0.05,
  engine = c("bayes_laplace", "cpp_laplace", "lme4", "lme4_fast", "glmmTMB", "stan_hmc", "stan_laplace"),
  analysis_re = c("match_dgp", "maximal", "uncorrelated", "intercept_only"),
  prob_se_method = c("bootstrap", "joint_delta"),
  control = target_control()
) {
  if (!inherits(dgp, "glmmplanr_dgp")) {
    cli::cli_abort("{.arg dgp} must be a {.cls glmmplanr_dgp} object.")
  }
  if (!inherits(control, "glmmplanr_target_control")) {
    cli::cli_abort("{.arg control} must be created by {.fn target_control}.")
  }

  engine_was_missing <- missing(engine)
  scale <- match.arg(scale)
  if (scale == "percentage_point") scale <- "pp"
  criterion <- match.arg(criterion)
  engine <- match.arg(engine)
  requested_analysis_re <- match.arg(analysis_re)
  analysis_re <- resolve_analysis_re(dgp, requested_analysis_re)
  prob_se_method <- match.arg(prob_se_method)

  prob_estimand <- control$prob_estimand
  contrast_weights <- control$contrast_weights
  ci_assurance <- control$ci_assurance
  rope <- control$rope
  failure_policy <- control$failure_policy
  min_valid_rate <- control$min_valid_rate
  mc_conf_level <- control$mc_conf_level
  stage1_engine <- control$stage1_engine
  stage1_re <- control$stage1_re
  on_singular <- control$on_singular
  analysis_link <- control$analysis_link
  lkj_eta <- control$lkj_eta
  prior_beta_sd <- control$prior_beta_sd
  prior_tau_sd <- control$prior_tau_sd
  prob_boot <- control$prob_boot

  chance_use <- if (analysis_link == "match_dgp") dgp$chance else 0.0

  resolved_stage1_re <- if (stage1_re == "uncorrelated" && analysis_re == "maximal") {
    "uncorrelated"
  } else {
    analysis_re
  }

  effect_col <- resolve_effect_column(dgp$design, effect)
  if (!length(dgp$design$factor_names) && scale != "logit" &&
      criterion %in% c("power", "both", "joint")) {
    cli::cli_abort("The single cell gives an absolute probability. Use {.code criterion = 'ci_width'} for precision or {.code scale = 'logit'} to test a zero conditional intercept.")
  }
  true_beta <- unname(dgp$beta[[effect_col]])
  c_weights <- cell_contrast_weights(dgp$design, effect_col, custom_weights = contrast_weights)
  true_prob_diff <- true_prob_effect(dgp, weights = c_weights, estimand = prob_estimand)
  true_effect_target <- switch(
    scale,
    logit = true_beta,
    prob = true_prob_diff,
    pp = true_prob_diff * 100.0
  )

  if (is.null(rope)) {
    rope <- switch(
      scale,
      logit = c(-0.18, 0.18),
      prob = c(-0.03, 0.03),
      pp = c(-3.0, 3.0)
    )
  }

  if (!is.numeric(power_min) || length(power_min) != 1L || is.na(power_min) ||
      power_min <= 0 || power_min >= 1) {
    cli::cli_abort("{.arg power_min} must be a probability strictly in (0, 1).")
  }
  if (criterion %in% c("ci_width", "both", "joint")) {
    if (is.null(ci_width) || !is.numeric(ci_width) || length(ci_width) != 1L ||
        is.na(ci_width) || ci_width <= 0) {
      cli::cli_abort(
        "{.arg ci_width} must be a positive scalar (on {.val {scale}} scale) when {.arg criterion} is {.val {criterion}}."
      )
    }
    if (scale == "prob" && ci_width >= 2.0) {
      cli::cli_warn(
        "{.arg ci_width} = {ci_width} with {.code scale = 'prob'} exceeds 1.0; did you mean {.code scale = 'pp'} (percentage points)?"
      )
    }
  }
  if (!is.numeric(ci_assurance) || length(ci_assurance) != 1L || is.na(ci_assurance) ||
      ci_assurance <= 0 || ci_assurance >= 1) {
    cli::cli_abort("{.arg ci_assurance} must lie strictly in (0, 1).")
  }
  if (!is.numeric(rope) || length(rope) != 2L || anyNA(rope) || rope[1L] >= rope[2L]) {
    cli::cli_abort("{.arg rope} must be a numeric vector c(low, high) with low < high.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || alpha <= 0 || alpha >= 0.5) {
    cli::cli_abort("{.arg alpha} must be a significance level strictly in (0, 0.5).")
  }
  if (!is.numeric(min_valid_rate) || length(min_valid_rate) != 1L || is.na(min_valid_rate) ||
      min_valid_rate <= 0 || min_valid_rate > 1) {
    cli::cli_abort("{.arg min_valid_rate} must lie in (0, 1].")
  }
  if (!is.numeric(mc_conf_level) || length(mc_conf_level) != 1L || is.na(mc_conf_level) ||
      mc_conf_level <= 0.5 || mc_conf_level >= 1) {
    cli::cli_abort("{.arg mc_conf_level} must lie strictly in (0.5, 1).")
  }

  if (identical(engine, "stan_laplace") && !engine_was_missing) {
    cli::cli_warn(c(
      "{.code engine = 'stan_laplace'} is a deprecated alias for {.code 'stan_hmc'}.",
      "i" = "It runs full Stan HMC/NUTS rather than Laplace approximation. Specify {.code engine = 'stan_hmc'} for Stan HMC or {.code engine = 'bayes_laplace'} for C++ regularized Laplace."
    ))
  }

  resolved_engine <- resolve_available_engine(
    engine, chance = chance_use, notify = !engine_was_missing
  )
  engine_fallback <- !identical(resolved_engine, engine)
  if (is.null(stage1_engine)) {
    stage1_engine <- if (resolved_engine %in% c("stan_laplace", "stan_hmc", "bayes_laplace")) {
      "bayes_laplace"
    } else {
      "cpp_laplace"
    }
  } else {
    stage1_engine <- resolve_available_engine(
      stage1_engine, chance = chance_use, notify = FALSE
    )
  }

  structure(
    list(
      dgp = dgp,
      effect_input = effect,
      effect_col = effect_col,
      scale = scale,
      prob_estimand = prob_estimand,
      contrast_weights = c_weights,
      true_beta = true_beta,
      true_prob_diff = true_prob_diff,
      true_prob_diff_pp = true_prob_diff * 100.0,
      true_effect_target = true_effect_target,
      criterion = criterion,
      power_min = power_min,
      ci_width = ci_width,
      ci_assurance = ci_assurance,
      rope = rope,
      alpha = alpha,
      failure_policy = failure_policy,
      min_valid_rate = min_valid_rate,
      mc_conf_level = mc_conf_level,
      requested_engine = engine,
      engine = resolved_engine,
      engine_fallback = engine_fallback,
      stage1_engine = stage1_engine,
      stage1_re = resolved_stage1_re,
      on_singular = on_singular,
      requested_analysis_re = requested_analysis_re,
      analysis_re = analysis_re,
      analysis_link = analysis_link,
      lkj_eta = lkj_eta,
      prior_beta_sd = prior_beta_sd,
      prob_boot = as.integer(prob_boot),
      prob_se_method = prob_se_method,
      prior_tau_sd = prior_tau_sd
    ),
    class = "glmmplanr_spec"
  )
}

resolve_effect_column <- function(design, effect) {
  if (!is.character(effect) || length(effect) != 1L || !nzchar(effect)) {
    cli::cli_abort("{.arg effect} must be a single non-empty character string.")
  }
  col_names <- colnames(design$X_cells)
  if (effect %in% col_names) {
    return(effect)
  }
  parts <- trimws(strsplit(effect, ":", fixed = TRUE)[[1L]])
  for (term_lbl in names(design$effect_map)) {
    term_parts <- trimws(strsplit(term_lbl, ":", fixed = TRUE)[[1L]])
    if (length(parts) == length(term_parts) && setequal(parts, term_parts)) {
      cols <- design$effect_map[[term_lbl]]
      if (length(cols) > 1L) {
        cli::cli_inform(c(
          "i" = "Term {.val {effect}} spans {length(cols)} contrast columns ({.val {cols}}); targeting first column {.val {cols[1L]}}.",
          "i" = "Pass a specific column name in {.arg effect} to target another contrast."
        ))
      }
      return(cols[1L])
    }
  }
  cli::cli_abort(c(
    "Effect {.val {effect}} not found in design.",
    "i" = "Available fixed-effect terms: {.val {names(design$effect_map)}}",
    "i" = "Available model-matrix columns: {.val {col_names}}"
  ))
}

has_cmdstan_backend <- function() {
  isTRUE(
    (requireNamespace("instantiate", quietly = TRUE) &&
       tryCatch(instantiate::stan_cmdstan_exists(), error = function(e) FALSE)) ||
      (requireNamespace("cmdstanr", quietly = TRUE) &&
         !is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)))
  )
}

resolve_analysis_re <- function(dgp, requested_re) {
  if (!identical(requested_re, "match_dgp")) {
    return(requested_re)
  }
  has_slopes <- (length(dgp$tau_subj) > 1L && any(dgp$tau_subj[-1L] > 0)) ||
    (length(dgp$tau_item) > 1L && any(dgp$tau_item[-1L] > 0)) ||
    (isTRUE(dgp$design$has_testlet) && length(dgp$tau_testlet) > 1L && any(dgp$tau_testlet[-1L] > 0))
  if (!has_slopes) {
    return("intercept_only")
  }
  has_corr <- function(R, tau) {
    if (is.null(R) || length(tau) <= 1L) return(FALSE)
    pos <- which(tau > 0)
    if (length(pos) <= 1L) return(FALSE)
    sub <- R[pos, pos, drop = FALSE]
    any(abs(sub[row(sub) != col(sub)]) > 1e-10)
  }
  if (has_corr(dgp$R_subj, dgp$tau_subj) ||
      has_corr(dgp$R_item, dgp$tau_item) ||
      (isTRUE(dgp$design$has_testlet) && has_corr(dgp$R_testlet, dgp$tau_testlet))) {
    "maximal"
  } else {
    "uncorrelated"
  }
}

resolve_available_engine <- function(engine, chance = 0.0, notify = TRUE) {
  if (engine == "glmmTMB") {
    if (chance > 0.0) {
      cli::cli_abort(c(
        "{.pkg glmmTMB} does not support custom guessing-adjusted links ({.code chance = {round(chance, 4)} > 0}).",
        "i" = "Choose {.code engine = 'bayes_laplace'}, {.code 'cpp_laplace'}, {.code 'lme4'}, or {.code 'stan_hmc'} for guessing links, or set {.code control = target_control(analysis_link = 'logit')} to test a canonical logit analysis."
      ))
    }
    if (!requireNamespace("glmmTMB", quietly = TRUE)) {
      cli::cli_abort(c(
        "{.code engine = 'glmmTMB'} requires the optional {.pkg glmmTMB} package (no silent fallback is performed so results stay numerically consistent).",
        "i" = "Install it with {.code install.packages('glmmTMB')} or select a built-in engine ({.code 'bayes_laplace'}, {.code 'cpp_laplace'}, or {.code 'lme4'})."
      ))
    }
    return("glmmTMB")
  }
  if (engine %in% c("stan_laplace", "stan_hmc")) {
    if (has_cmdstan_backend() || requireNamespace("rstan", quietly = TRUE)) {
      return(engine)
    }
    cli::cli_abort(c(
      "{.code engine = '{engine}'} requires a Stan backend ({.pkg cmdstanr} or {.pkg rstan}); no silent fallback to Laplace approximation is performed.",
      "i" = "Install {.pkg rstan} or {.pkg cmdstanr}, or explicitly select {.code engine = 'bayes_laplace'}."
    ))
  }
  engine
}

#' @rdname print.glmmplanr_design
#' @export
print.glmmplanr_spec <- function(x, ...) {
  cli::cli_h2("glmmplanr Target & Estimation Specification")
  unit_lbl <- switch(x$scale, logit = "logit", pp = "pp", prob = "prob")
  crit_desc <- switch(
    x$criterion,
    power = sprintf("Power >= %.2f (alpha = %.3f)", x$power_min, x$alpha),
    ci_width = sprintf("P(CI width <= %.3f %s) >= %.2f", x$ci_width, unit_lbl, x$ci_assurance),
    both = sprintf(
      "Power >= %.2f AND P(CI width <= %.3f %s) >= %.2f",
      x$power_min, x$ci_width, unit_lbl, x$ci_assurance
    ),
    joint = sprintf("P(detect and width <= %.3f %s) >= %.2f", x$ci_width, unit_lbl, x$power_min),
    rope = sprintf(
      "P(CI excludes ROPE [%.2f, %.2f] %s) >= %.2f",
      x$rope[1L], x$rope[2L], unit_lbl, x$power_min
    )
  )
  re_lbl <- if (identical(x$requested_analysis_re, "match_dgp")) {
    sprintf("%s (matched to DGP)", x$analysis_re)
  } else {
    x$analysis_re
  }
  cli::cli_ul(c(
    "Focal effect: {.val {x$effect_input}} (true beta = {round(x$true_beta, 4)} logit | {x$prob_estimand} diff = {sprintf('%+.2f pp', x$true_prob_diff_pp)})",
    "Target scale: {.val {x$scale}} | Criterion ({.val {x$criterion}}): {crit_desc}",
    "Failure accounting: {.val {x$failure_policy}} | MC confidence level: {round(x$mc_conf_level * 100, 1)}%",
    "Engine: Stage 1 = {.val {x$stage1_engine}} -> Final = {.val {x$engine}} | Analysis RE: {.val {re_lbl}}"
  ))
  cli::cli_text("{analysis_description(x)}")
  cli::cli_text("Model: {.code {deparse1(stats::formula(x))}}")
  invisible(x)
}
