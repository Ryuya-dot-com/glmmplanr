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
#'   scale the one cell gives an absolute probability. Use `ci_width` for its
#'   precision. Variance components and individual scores require a separate analysis.
#' @param scale Scale on which detection, precision (`ci_width`) and equivalence (`rope`)
#'   targets are expressed and evaluated:
#'   * `"logit"` (default): Conditional log-odds scale.
#'   * `"pp"` (or `"percentage_point"`): Response accuracy difference in
#'     **percentage points** (e.g., `ci_width = 10` means a full 95% CI width of
#'     at most 10 percentage points, such as `+5.0 pp [0.0, 10.0]`).
#'   * `"prob"`: Response accuracy difference on `[0, 1]` probability scale
#'     (e.g., `ci_width = 0.10` for 10 percentage points).
#' @param prob_estimand Estimand used for probability-difference (`"pp"` /
#'   `"prob"`) calculations (and always reported diagnostically even when
#'   `scale = "logit"`):
#'   * `"marginal"` (default): Population-standardized (G-computation) accuracy
#'     difference averaging over participant and item random-effect variances via
#'     Gauss-Hermite quadrature.
#'   * `"conditional"`: Accuracy difference for an average participant and
#'     average item (\eqn{u_S = 0, u_I = 0}).
#' @param contrast_weights Optional custom numeric vector of cell contrast
#'   weights of length `dgp$design$n_cells` (defaults to
#'   [cell_contrast_weights()]).
#' @param criterion Evaluation objective:
#'   * `"power"`: Statistical power \eqn{\ge} `power_min`.
#'   * `"ci_width"`: Precision / AIPE criterion requiring the `(1 - alpha)`
#'     interval width to be \eqn{\le} `ci_width` with probability \eqn{\ge}
#'     `ci_assurance`.
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
#' @param ci_assurance Target probability \eqn{\Pr(W \le w_{\max})} for AIPE
#'   interval width assurance (default `0.80`).
#' @param rope Numeric vector `c(low, high)` defining the Region of Practical
#'   Equivalence in units matching `scale` (defaults to `c(-0.18, 0.18)` on
#'   `"logit"`, `c(-3, 3)` on `"pp"`, or `c(-0.03, 0.03)` on `"prob"`).
#' @param alpha Two-sided significance level / `1 - conf_level` (default `0.05`).
#' @param failure_policy How non-converged or failed replications enter sample
#'   size pass/fail decisions:
#'   * `"effective"` (default): Conservative unconditional probability of
#'     obtaining a valid converged fit that also meets the criterion (failed
#'     replications count as `0` / non-attainment), while still reporting
#'     `conditional_power` alongside.
#'   * `"conditional"`: Evaluates `conditional_power` among valid fits while
#'     requiring `valid_rate >= min_valid_rate`.
#' @param min_valid_rate Minimum acceptable proportion of valid converged fits
#'   when `failure_policy = "conditional"` (default `0.90`).
#' @param mc_conf_level Confidence level for Wilson Monte Carlo uncertainty
#'   intervals used to classify borderline candidates as `"pass"`,
#'   `"indeterminate"`, or `"fail"` (default `0.95`; Morris et al., 2019).
#' @param engine Primary estimation backend: `"glmmTMB"`, `"lme4_fast"`,
#'   `"lme4"`, `"stan_hmc"`, `"bayes_laplace"`, `"cpp_laplace"`, or legacy alias
#'   `"stan_laplace"` (deprecated; executes full HMC/NUTS via `"stan_hmc"`
#'   rather than Laplace approximation). Unpenalized engines flag separation if
#'   all observed responses in any fixed-effect design cell are 0 or all are 1.
#'   Each cell has its own linear predictor in the factorial model. The primary
#'   fit is then invalid, and a bootstrap refit with this pattern is also
#'   invalid. Regularized engines can return finite estimates in these cases.
#'   Check sensitivity to the fixed-effect prior. An individual participant or
#'   item with a constant response pattern can contribute to either fit.
#' @param stage1_engine Fast estimation backend used during Stage-1 surrogate
#'   exploration before Stage-2 verification.
#' @param stage1_re Random-effect covariance structure used during fast Stage-1
#'   surrogate exploration before Stage-2 verification: `"match"` (default; uses
#'   `analysis_re`) or `"uncorrelated"` (uses fast diagonal random slopes in
#'   Stage 1 while still verifying the full `analysis_re` model in Stage 2).
#' @param on_singular Policy when maximum-likelihood fits hit a singular
#'   random-effects boundary: `"flag_and_compare"` (default), `"flag_only"`,
#'   `"reduce_corr"`, or `"switch_bayes"`.
#' @param analysis_re Random-effect structure of the fitted analysis model:
#'   `"maximal"` (default), `"uncorrelated"`, or `"intercept_only"`.
#' @param analysis_link Link function used in the analysis model: `"match_dgp"`
#'   (default) or `"logit"`.
#' @param lkj_eta Shape parameter \eqn{\eta} for the LKJ correlation prior in
#'   Stan (default `2.0`); in `"bayes_laplace"`, controls unit-diagonal
#'   correlation shrinkage \eqn{n / (n + 2\max(\eta - 1, 0))} toward the
#'   identity matrix.
#' @param prior_beta_sd Weakly informative Normal prior SD on non-intercept
#'   fixed effects in Bayesian engines (default `2.5`).
#' @param prior_tau_sd Scale parameter for the Half-Student-t(\eqn{\nu = 3})
#'   prior on random-effect SDs in both Stan (`"stan_hmc"`) and the built-in
#'   regularized Laplace engine (`"bayes_laplace"`), default `1.5`. Whereas
#'   Stan performs full joint HMC/NUTS over \eqn{(\beta, \tau, L, z)},
#'   `"bayes_laplace"` uses an EM-Laplace approximation with exact
#'   \eqn{\log \tau}-space Half-Student-t(3) MAP variance updates and
#'   LKJ-scaled correlation shrinkage, so finite-sample intervals can still
#'   differ between the two engines.
#' @param prob_boot Number of parametric bootstrap refits per final marginal
#'   probability-scale replication when `prob_se_method = "bootstrap"`
#'   (default `99`). At least 90% of refits must succeed.
#' @param prob_se_method Standard-error method for marginal probability-scale
#'   contrasts in Stage-2 final verification:
#'   * `"bootstrap"` (default): Parametric bootstrap refits (`prob_boot` draws)
#'     executed via a pre-allocated C++ batch kernel for built-in engines.
#'   * `"joint_delta"`: Analytical Gauss-Hermite delta method propagating joint
#'     uncertainty in fixed effects \eqn{\beta} and random-effect log-SDs
#'     \eqn{\theta = \log \tau} (also used automatically in Stage-1 screening).
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
  prob_estimand = c("marginal", "conditional"),
  contrast_weights = NULL,
  criterion = c("power", "ci_width", "both", "rope", "joint"),
  power_min = 0.80,
  ci_width = NULL,
  ci_assurance = 0.80,
  rope = NULL,
  alpha = 0.05,
  failure_policy = c("effective", "conditional"),
  min_valid_rate = 0.90,
  mc_conf_level = 0.95,
  engine = c("glmmTMB", "lme4_fast", "lme4", "stan_laplace", "stan_hmc", "bayes_laplace", "cpp_laplace"),
  stage1_engine = NULL,
  stage1_re = c("match", "uncorrelated"),
  on_singular = c("flag_and_compare", "flag_only", "reduce_corr", "switch_bayes"),
  analysis_re = c("maximal", "uncorrelated", "intercept_only"),
  analysis_link = c("match_dgp", "logit"),
  lkj_eta = 2.0,
  prior_beta_sd = 2.5,
  prior_tau_sd = 1.5,
  prob_boot = 99L,
  prob_se_method = c("bootstrap", "joint_delta")
) {
  if (!inherits(dgp, "glmmplanr_dgp")) {
    cli::cli_abort("{.arg dgp} must be a {.cls glmmplanr_dgp} object.")
  }
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
  engine_was_missing <- missing(engine)
  scale <- match.arg(scale)
  if (scale == "percentage_point") scale <- "pp"
  prob_estimand <- match.arg(prob_estimand)
  criterion <- match.arg(criterion)
  failure_policy <- match.arg(failure_policy)
  engine <- match.arg(engine)
  stage1_re <- match.arg(stage1_re)
  on_singular <- match.arg(on_singular)
  analysis_re <- match.arg(analysis_re)
  analysis_link <- match.arg(analysis_link)
  prob_se_method <- match.arg(prob_se_method)

  chance_use <- if (analysis_link == "match_dgp") dgp$chance else 0.0
  if (engine_was_missing && chance_use > 0.0) {
    engine <- "cpp_laplace"
  }

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

resolve_available_engine <- function(engine, chance = 0.0, notify = TRUE) {
  if (engine == "glmmTMB") {
    if (chance > 0.0) {
      if (notify) {
        cli::cli_inform(c(
          "i" = "{.pkg glmmTMB} only supports built-in C++ binomial links and rejects custom m-AFC lower-asymptote links ({.code chance = {round(chance, 4)} > 0}); delegating to {.val cpp_laplace}.",
          "i" = "Use {.code engine = 'cpp_laplace'}, {.code 'bayes_laplace'}, {.code 'lme4'}, or {.code 'stan_hmc'} directly, or set {.code analysis_link = 'logit'} to test canonical logit mismatch."
        ))
      }
      return("cpp_laplace")
    }
    if (requireNamespace("glmmTMB", quietly = TRUE)) {
      return("glmmTMB")
    }
    if (notify) {
      cli::cli_inform(c(
        "i" = "{.pkg glmmTMB} is not installed; delegating {.val glmmTMB} to glmmplanr's compiled C++ Laplace engine ({.val cpp_laplace})."
      ))
    }
    return("cpp_laplace")
  }
  if (engine %in% c("stan_laplace", "stan_hmc")) {
    if (has_cmdstan_backend() || requireNamespace("rstan", quietly = TRUE)) {
      return(engine)
    }
    if (notify) {
      cli::cli_warn(c(
        "No Stan backend ({.pkg cmdstanr}/{.pkg rstan}) is installed; falling back from {.val {engine}} to {.val bayes_laplace}.",
        "i" = "{.val bayes_laplace} uses EM-Laplace approximation rather than HMC/NUTS; finite-sample interval widths can differ from Stan."
      ))
    }
    return("bayes_laplace")
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
  cli::cli_ul(c(
    "Focal effect: {.val {x$effect_input}} (true beta = {round(x$true_beta, 4)} logit | {x$prob_estimand} diff = {sprintf('%+.2f pp', x$true_prob_diff_pp)})",
    "Target scale: {.val {x$scale}} | Criterion ({.val {x$criterion}}): {crit_desc}",
    "Failure accounting: {.val {x$failure_policy}} | MC confidence level: {round(x$mc_conf_level * 100, 1)}%",
    "Engine: Stage 1 = {.val {x$stage1_engine}} -> Final = {.val {x$engine}} (requested: {.val {x$requested_engine}})"
  ))
  cli::cli_text("{analysis_description(x)}")
  cli::cli_text("Model: {.code {deparse1(stats::formula(x))}}")
  invisible(x)
}
