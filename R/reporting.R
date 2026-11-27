#' Render Sample Size Optimization or Stress-Test Results as a Publication-Ready `gt` Table
#'
#' Formats verified Pareto-frontier sample sizes, scheduled vs. post-attrition
#' trial counts, effective and conditional power with Monte Carlo uncertainty
#' intervals, percentage-point and logit 95% interval widths, convergence /
#' singularity rates, multi-scenario robustness increments, and random-slope
#' sensitivity analyses into a clean `gt` table.
#'
#' @param x A GLMM optimization, Rasch planning, or stress-test result.
#' @param which Which table section(s) to display for `"glmmplanr_result"`:
#'   `"combined"` (default), `"verified"`, `"pareto"`, `"sensitivity"`, or
#'   `"robust"`.
#' @param ... Reserved for method extensions.
#' @return A `gt_tbl` object created by [gt::gt()].
#' @examples
#' # Create the five-replication demonstration plan from the optimizer help.
#' example("optimize_sample_size", package = "glmmplanr", local = environment(), echo = FALSE)
#' table <- as_gt(plan)
#' class(table)
#' @export
as_gt <- function(x, ...) {
  UseMethod("as_gt")
}

assess_inferential_validity <- function(cov, n_valid, alpha, has_null = TRUE) {
  nominal <- 1.0 - alpha
  if (is.null(cov) || is.null(n_valid) || length(cov) == 0L || length(n_valid) == 0L ||
      !is.finite(cov[1L]) || !is.finite(n_valid[1L]) || n_valid[1L] < 1L) {
    return(list(
      status = "unverified",
      short_label = "Unverified",
      detail = "UNVERIFIED (no valid fits to inspect coverage; run validate_design_point())",
      cov = NA_real_, cov_lo = NA_real_, cov_hi = NA_real_, nominal = nominal
    ))
  }
  cov1 <- cov[1L]
  nv1 <- as.integer(n_valid[1L])
  k_cov <- max(0L, min(nv1, as.integer(round(cov1 * nv1))))
  ci <- wilson_mc_ci(k_cov, nv1, conf_level = 0.95)
  cov_lo <- unname(ci["lo"])
  cov_hi <- unname(ci["hi"])
  if (nv1 >= 20L && is.finite(cov_hi) && cov_hi < (nominal - 0.02)) {
    return(list(
      status = "caution_low_coverage",
      short_label = "Caution: Low Cov.",
      detail = sprintf(
        "CAUTION (Stage-2 coverage %.1f%% [%.1f%%, %.1f%%] is below nominal %.0f%%; SEs may be anti-conservative; verify with validate_design_point())",
        100 * cov1, 100 * cov_lo, 100 * cov_hi, 100 * nominal
      ),
      cov = cov1, cov_lo = cov_lo, cov_hi = cov_hi, nominal = nominal
    ))
  }
  null_note <- if (isTRUE(has_null)) {
    "Stage-2 H1 coverage = %.1f%% [%.1f%%, %.1f%%]; null Type I error not checked; confirm with validate_design_point()"
  } else {
    "Stage-2 coverage = %.1f%% [%.1f%%, %.1f%%]; confirm planned-condition coverage and width attainment on a fresh seed with validate_design_point()"
  }
  list(
    status = "unverified",
    short_label = "Unverified",
    detail = sprintf(paste0("UNVERIFIED (", null_note, ")"), 100 * cov1, 100 * cov_lo, 100 * cov_hi),
    cov = cov1, cov_lo = cov_lo, cov_hi = cov_hi, nominal = nominal
  )
}

inferential_validity_labels <- function(cov_vec, n_valid_vec, alpha) {
  vapply(seq_along(cov_vec), function(i) {
    assess_inferential_validity(cov_vec[i], n_valid_vec[i], alpha)$short_label
  }, character(1L))
}

detect_exploratory_regimes <- function(spec) {
  dgp <- spec$dgp
  des <- dgp$design
  reasons <- character()
  has_slopes <- (!identical(spec$analysis_re, "intercept_only")) ||
    (length(dgp$tau_subj) > 1L && any(abs(dgp$tau_subj[-1L]) > 1e-12)) ||
    (length(dgp$tau_item) > 1L && any(abs(dgp$tau_item[-1L]) > 1e-12))
  if (!is.null(des$items_per_subject)) {
    if (has_slopes) {
      reasons <- c(
        reasons,
        sprintf("partial item presentation (items_per_subject = %d) with random slopes", des$items_per_subject)
      )
    } else {
      reasons <- c(
        reasons,
        sprintf("partial item presentation (items_per_subject = %d)", des$items_per_subject)
      )
    }
  }
  if (isTRUE(des$has_testlet)) {
    pt_sd <- testlet_interaction_sd(dgp)
    if (pt_sd > 0) {
      reasons <- c(
        reasons,
        sprintf("testlet nesting with participant-by-passage local dependence (SD = %.2f; can inflate Type I error)", pt_sd)
      )
    } else {
      reasons <- c(
        reasons,
        "passage/testlet nesting (design_testlet; small passage counts can inflate asymptotic Wald Type I error)"
      )
    }
  }
  n_colons <- lengths(regmatches(spec$effect_col, gregexpr(":", spec$effect_col, fixed = TRUE)))
  if (isTRUE(n_colons >= 2L)) {
    reasons <- c(reasons, sprintf("higher-order interaction (%s)", spec$effect_input))
  }
  reasons
}

assess_boundary_sensitivity <- function(x) {
  if (!isTRUE(x$feasible) || is.null(x$optimal) || nrow(x$optimal) == 0L) {
    return(NULL)
  }
  opt <- x$optimal[1L, , drop = FALSE]
  if (!identical(opt$decision_status, "pass")) return(NULL)
  sp <- x$spec
  if (!is.null(x$robust$specs) && !is.null(opt$scenario) &&
      (opt$scenario %in% names(x$robust$specs))) {
    sp <- x$robust$specs[[opt$scenario]]
  }
  margin_of_row <- function(r) {
    switch(
      sp$criterion,
      power = r$power_mc_lo - sp$power_min,
      ci_width = r$aipe_mc_lo - sp$ci_assurance,
      both = pmin(r$power_mc_lo - sp$power_min, r$aipe_mc_lo - sp$ci_assurance),
      joint = r$joint_mc_lo - sp$power_min,
      rope = r$rope_mc_lo - sp$power_min
    )
  }
  opt_margin <- margin_of_row(opt)
  if (!is.finite(opt_margin)) return(NULL)

  ver <- x$verified
  indet_comp <- if (!is.null(ver) && nrow(ver) > 0L) {
    ver[
      ver$decision_status == "indeterminate" &
        ver$n_trials_scheduled >= 0.75 * opt$n_trials_scheduled &
        ver$n_trials_scheduled <= 1.50 * opt$n_trials_scheduled,
      ,
      drop = FALSE
    ]
  } else {
    NULL
  }
  n_indet <- if (is.null(indet_comp)) 0L else nrow(indet_comp)

  near_boundary <- (opt_margin < 0.05) ||
    (!is.null(ver) && any(ver$decision_status == "indeterminate" &
                          ver$n_trials_scheduled >= opt$n_trials_scheduled))
  if (!near_boundary) return(NULL)

  buffer_row <- NULL
  if (!is.null(ver) && nrow(ver) > 0L) {
    pass_more <- ver[
      ver$decision_status == "pass" &
        ver$n_trials_scheduled > opt$n_trials_scheduled,
      ,
      drop = FALSE
    ]
    if (nrow(pass_more) > 0L) {
      m_more <- margin_of_row(pass_more)
      better <- which(is.finite(m_more) & m_more > opt_margin + 0.01)
      if (length(better) > 0L) {
        sub_b <- pass_more[better, , drop = FALSE]
        pref_rk <- integer(nrow(sub_b))
        pref_rk[preference_order(sub_b, prefer = x$prefer)] <- seq_len(nrow(sub_b))
        tr_col <- if ("n_trials_scheduled_max" %in% names(sub_b)) sub_b$n_trials_scheduled_max else sub_b$n_trials_scheduled
        ord <- order(tr_col, pref_rk)
        buffer_row <- sub_b[ord[1L], , drop = FALSE]
        buffer_row$mc_margin <- margin_of_row(buffer_row)
      }
    }
  }

  list(
    opt_margin = opt_margin,
    n_indeterminate_nearby = n_indet,
    buffer = buffer_row
  )
}

#' @rdname as_gt
#' @export
as_gt.glmmplanr_result <- function(
  x,
  which = c("combined", "verified", "pareto", "sensitivity", "robust"),
  ...
) {
  which <- match.arg(which)
  ver <- x$verified[order(x$verified$n_item, x$verified$n_subj), , drop = FALSE]

  is_opt <- (ver$n_subj == x$optimal$n_subj) & (ver$n_item == x$optimal$n_item)
  opt_idx <- which(is_opt)
  if (length(opt_idx) > 1L) is_opt[-opt_idx[1L]] <- FALSE

  pareto_keys <- if (!is.null(x$verified_pareto) && nrow(x$verified_pareto) > 0L) {
    paste(x$verified_pareto$n_subj, x$verified_pareto$n_item)
  } else {
    character(0)
  }
  is_par <- paste(ver$n_subj, ver$n_item) %in% pareto_keys

  role_lbl <- ifelse(
    is_opt,
    if (isTRUE(x$feasible)) sprintf("Recommended (%s) *", x$prefer) else "Unconfirmed candidate",
    ifelse(is_par, "Pareto Frontier", "Candidate")
  )

  fmt_mc_pow <- function(p, lo, hi) {
    ifelse(
      is.na(p),
      "-",
      sprintf("%.1f%% [%.1f, %.1f]", p * 100, lo * 100, hi * 100)
    )
  }

  ver_cov <- if (x$spec$scale == "logit") ver$coverage else ver$prob_coverage
  ver_clean <- data.frame(
    section = "Stage-2 Verified Candidate & Pareto Designs",
    label = role_lbl,
    n_subj = ver$n_subj,
    n_item = ver$n_item,
    n_trials_scheduled = ver$n_trials_scheduled,
    n_trials_expected = ver$n_trials_expected,
    decision = toupper(ver$decision_status),
    eff_power_ci = fmt_mc_pow(ver$power, ver$power_mc_lo, ver$power_mc_hi),
    rope_ci = if ("rope_excl_prob" %in% names(ver)) fmt_mc_pow(ver$rope_excl_prob, ver$rope_mc_lo, ver$rope_mc_hi) else rep("-", nrow(ver)),
    conditional_power = ver$conditional_power,
    prob_ci_width_pp = ver$prob_ci_width_mean_pp,
    logit_ci_width = ver$logit_ci_width_mean,
    aipe_prob = ver$aipe_prob,
    joint_prob = if (is.null(ver$joint_prob)) NA_real_ else ver$joint_prob,
    validity_status = inferential_validity_labels(ver_cov, ver$n_valid, x$spec$alpha),
    valid_rate = ver$valid_rate,
    singular_rate = ver$singular_rate,
    coverage = ver_cov,
    stringsAsFactors = FALSE
  )

  pareto_only <- ver_clean[is_par, , drop = FALSE]

  sens_clean <- NULL
  if (!is.null(x$sensitivity)) {
    s <- x$sensitivity
    cap <- paste0(toupper(substring(s$scenario, 1L, 1L)), substring(s$scenario, 2L))
    s_cov <- if (x$spec$scale == "logit") s$coverage else s$prob_coverage
    sens_src_lbl <- attr(s, "source_scenario")
    if ((is.null(sens_src_lbl) || !length(sens_src_lbl)) && "source_scenario" %in% names(s)) {
      sens_src_lbl <- as.character(s$source_scenario[1L])
    }
    sens_scen_tag <- if (!is.null(sens_src_lbl) && length(sens_src_lbl) > 0L &&
                         !is.na(sens_src_lbl[1L]) && nzchar(sens_src_lbl[1L])) {
      sprintf(" [scenario: %s]", sens_src_lbl[1L])
    } else {
      ""
    }
    sens_clean <- data.frame(
      section = sprintf(
        "Random-Slope Sensitivity at (N_subj = %d, N_item = %d)%s",
        x$optimal$n_subj, x$optimal$n_item, sens_scen_tag
      ),
      label = sprintf("%s (%.2fx slope SD)", cap, s$slope_multiplier),
      n_subj = s$n_subj,
      n_item = s$n_item,
      n_trials_scheduled = s$n_trials_scheduled,
      n_trials_expected = s$n_trials_expected,
      decision = toupper(s$decision_status),
      eff_power_ci = fmt_mc_pow(s$power, s$power_mc_lo, s$power_mc_hi),
      rope_ci = if ("rope_excl_prob" %in% names(s)) fmt_mc_pow(s$rope_excl_prob, s$rope_mc_lo, s$rope_mc_hi) else rep("-", nrow(s)),
      conditional_power = s$conditional_power,
      prob_ci_width_pp = s$prob_ci_width_mean_pp,
      logit_ci_width = s$logit_ci_width_mean,
      aipe_prob = s$aipe_prob,
      joint_prob = if (is.null(s$joint_prob)) NA_real_ else s$joint_prob,
      validity_status = inferential_validity_labels(s_cov, s$n_valid, x$spec$alpha),
      valid_rate = s$valid_rate,
      singular_rate = s$singular_rate,
      coverage = s_cov,
      stringsAsFactors = FALSE
    )
  }

  rob_clean <- NULL
  if (!is.null(x$robust)) {
    se <- x$robust$scenario_evaluations
    rc_verified <- isTRUE(x$robust$feasible_standard) && isTRUE(x$robust$feasible_robust)
    opt_s <- x$robust$optimal_standard
    opt_r <- x$robust$optimal_robust
    pick_rows <- se[
      (se$n_subj == opt_s$n_subj & se$n_item == opt_s$n_item) |
        (se$n_subj == opt_r$n_subj & se$n_item == opt_r$n_item),
      ,
      drop = FALSE
    ]
    if (nrow(pick_rows) > 0L) {
      p_cov <- if (x$spec$scale == "logit") pick_rows$coverage else pick_rows$prob_coverage
      rob_clean <- data.frame(
        section = if (!rc_verified) "Multi-Scenario Candidates (increment undetermined)" else sprintf(
          "Multi-Scenario Robustness (All-Scenario Increment: %+d Subj, %+d Items)",
          x$robust$delta_subj, x$robust$delta_item
        ),
        label = sprintf("%s @ (%d, %d)", pick_rows$scenario, pick_rows$n_subj, pick_rows$n_item),
        n_subj = pick_rows$n_subj,
        n_item = pick_rows$n_item,
        n_trials_scheduled = pick_rows$n_trials_scheduled,
        n_trials_expected = pick_rows$n_trials_expected,
        decision = toupper(pick_rows$decision_status),
        eff_power_ci = fmt_mc_pow(pick_rows$power, pick_rows$power_mc_lo, pick_rows$power_mc_hi),
        rope_ci = if ("rope_excl_prob" %in% names(pick_rows)) fmt_mc_pow(pick_rows$rope_excl_prob, pick_rows$rope_mc_lo, pick_rows$rope_mc_hi) else rep("-", nrow(pick_rows)),
        conditional_power = pick_rows$conditional_power,
        prob_ci_width_pp = pick_rows$prob_ci_width_mean_pp,
        logit_ci_width = pick_rows$logit_ci_width_mean,
        aipe_prob = pick_rows$aipe_prob,
        joint_prob = if (is.null(pick_rows$joint_prob)) NA_real_ else pick_rows$joint_prob,
        validity_status = inferential_validity_labels(p_cov, pick_rows$n_valid, x$spec$alpha),
        valid_rate = pick_rows$valid_rate,
        singular_rate = pick_rows$singular_rate,
        coverage = p_cov,
        stringsAsFactors = FALSE
      )
    }
  }

  if (identical(which, "sensitivity") && is.null(sens_clean)) {
    cli::cli_abort(c(
      "Random-slope sensitivity analysis was not run for {.arg x} (stored {.code sensitivity} is NULL).",
      "i" = "Run {.fn optimize_sample_size} with {.code nsim_sensitivity > 0} on a DGP with non-zero random slopes, or use {.code which = 'verified'} / {.code 'combined'}."
    ))
  }
  if (identical(which, "robust") && is.null(rob_clean)) {
    cli::cli_abort(c(
      "Multi-scenario robustness analysis was not run for {.arg x} (stored {.code robust} is NULL).",
      "i" = "Supply {.arg scenarios} to {.fn optimize_sample_size}, or use {.code which = 'verified'} / {.code 'combined'}."
    ))
  }

  tbl_df <- switch(
    which,
    verified = ver_clean,
    pareto = pareto_only,
    sensitivity = sens_clean,
    robust = rob_clean,
    combined = do.call(rbind, Filter(Negate(is.null), list(ver_clean, rob_clean, sens_clean)))
  )

  if (all(is.na(tbl_df$aipe_prob))) {
    tbl_df$aipe_prob <- tbl_df$joint_prob <- NULL
    has_aipe <- FALSE
  } else {
    has_aipe <- TRUE
  }

  size_cols <- c("n_subj", "n_item", "n_trials_scheduled", "n_trials_expected")

  sp_hdr <- x$spec
  scen_hdr_suffix <- ""
  if (!is.null(x$robust$specs) && !is.null(x$optimal$scenario) &&
      length(x$optimal$scenario) > 0L && (x$optimal$scenario[1L] %in% names(x$robust$specs))) {
    opt_sc <- as.character(x$optimal$scenario[1L])
    sp_hdr <- x$robust$specs[[opt_sc]]
    scen_hdr_suffix <- sprintf(
      " | scenario = %s [evaluated: %s]",
      opt_sc, paste(names(x$robust$specs), collapse = ", ")
    )
  }

  has_rope <- identical(sp_hdr$criterion, "rope")
  if (!has_rope) {
    tbl_df$rope_ci <- NULL
  }

  crit_hdr <- switch(
    sp_hdr$criterion,
    power = sprintf("criterion = power >= %.0f%%", sp_hdr$power_min * 100),
    ci_width = sprintf("criterion = ci_width <= %g %s (>= %.0f%%)", sp_hdr$ci_width, sp_hdr$scale, sp_hdr$ci_assurance * 100),
    both = sprintf("criterion = power >= %.0f%% & width <= %g %s (>= %.0f%%)", sp_hdr$power_min * 100, sp_hdr$ci_width, sp_hdr$scale, sp_hdr$ci_assurance * 100),
    joint = sprintf("criterion = joint detect & width <= %g %s >= %.0f%%", sp_hdr$ci_width, sp_hdr$scale, sp_hdr$power_min * 100),
    rope = sprintf("criterion = ROPE [%g, %g %s] excl. >= %.0f%%", sp_hdr$rope[1L], sp_hdr$rope[2L], sp_hdr$scale, sp_hdr$power_min * 100)
  )
  hdr_tail <- paste0(" | ", crit_hdr, scen_hdr_suffix)

  p_kind_gt <- prob_target_kind(sp_hdr)
  gt_sub_fmt <- if (identical(p_kind_gt, "absolute")) {
    "Focal target: %s (%s absolute probability = %.1f%% | cond. beta = %.3f logit | scale = %s | prefer = %s%s)"
  } else if (identical(p_kind_gt, "difference")) {
    "Focal effect: %s (%s diff = %+.1f pp | cond. beta = %.3f logit | scale = %s | prefer = %s%s)"
  } else if (can_define_null_target(sp_hdr)) {
    "Focal effect: %s (%s weighted contrast = %+.1f pp | cond. beta = %.3f logit | scale = %s | prefer = %s%s)"
  } else {
    "Focal target: %s (%s weighted contrast = %+.1f pp | cond. beta = %.3f logit | scale = %s | prefer = %s%s)"
  }
  gt_obj <- gt::gt(tbl_df, groupname_col = "section", rowname_col = "label") |>
    gt::tab_header(
      title = "Two-Dimensional GLMM Sample Size Plan & Trustworthiness Summary",
      subtitle = sprintf(
        gt_sub_fmt,
        sp_hdr$effect_input, sp_hdr$prob_estimand, sp_hdr$true_prob_diff_pp,
        sp_hdr$true_beta, sp_hdr$scale, x$prefer, hdr_tail
      )
    ) |>
    gt::cols_label(
      n_subj = "N_subj",
      n_item = "N_item",
      n_trials_scheduled = "Sched. Trials",
      n_trials_expected = "Exp. Obs",
      decision = "Target Status",
      eff_power_ci = if (sp_hdr$failure_policy == "effective") "Effective Power [MC CI]" else "Conditional Power [MC CI]",
      conditional_power = "Cond. Power",
      prob_ci_width_pp = "CI (pp)",
      logit_ci_width = "CI (logit)",
      validity_status = "Inferential Validity",
      valid_rate = "Valid",
      singular_rate = "Singular",
      coverage = if (can_define_null_target(sp_hdr)) "H1 Coverage" else "Coverage"
    ) |>
    gt::tab_spanner(
      label = "Balanced Sample Size",
      columns = size_cols
    ) |>
    gt::tab_spanner(
      label = "Inferential Validity & Estimation Health",
      columns = c("validity_status", "coverage", "valid_rate", "singular_rate")
    ) |>
    gt::fmt_number(
      columns = c("n_trials_scheduled", "n_trials_expected"),
      decimals = 0,
      use_seps = TRUE
    ) |>
    gt::fmt_percent(
      columns = c("conditional_power", "valid_rate", "singular_rate", "coverage"),
      decimals = 1
    ) |>
    gt::fmt_number(columns = "prob_ci_width_pp", decimals = 2) |>
    gt::fmt_number(columns = "logit_ci_width", decimals = 3) |>
    gt::sub_missing(missing_text = "-")

  if (has_rope) {
    rope_lbl <- sprintf(
      "%s ROPE Excl. [%g, %g %s] [MC CI]",
      if (sp_hdr$failure_policy == "effective") "Effective" else "Conditional",
      sp_hdr$rope[1L], sp_hdr$rope[2L], sp_hdr$scale
    )
    gt_obj <- gt::cols_label(gt_obj, rope_ci = rope_lbl)
  }

  dec_cols <- c(
    "decision", "eff_power_ci",
    if (has_rope) "rope_ci",
    "conditional_power", "prob_ci_width_pp", "logit_ci_width",
    if (has_aipe) c("aipe_prob", "joint_prob")
  )
  if (has_aipe) {
    gt_obj <- gt_obj |>
      gt::cols_label(aipe_prob = "Width attainment", joint_prob = "Detect and meet width") |>
      gt::fmt_percent(columns = c("aipe_prob", "joint_prob"), decimals = 1)
  }
  gt_obj <- gt_obj |>
    gt::tab_spanner(
      label = if (can_define_null_target(sp_hdr)) "Power, Precision & MC Decision" else "Precision & MC Decision",
      columns = intersect(dec_cols, names(tbl_df))
    )

  if (!can_define_null_target(sp_hdr)) {
    gt_obj <- gt::cols_hide(gt_obj, columns = intersect(c("eff_power_ci", "conditional_power", "joint_prob"), names(tbl_df)))
  }
  opt_sections_note <- if (identical(which, "combined")) {
    omitted <- c(
      if (is.null(rob_clean)) "multi-scenario robustness",
      if (is.null(sens_clean)) "random-slope sensitivity"
    )
    if (length(omitted) > 0L) {
      sprintf(" Not run in this plan: %s.", paste(omitted, collapse = ", "))
    } else {
      ""
    }
  } else {
    ""
  }
  sing_note <- if (any(is.na(tbl_df$singular_rate))) {
    " Boundary singularity ('Singular') is not applicable ('-') under regularized/Bayesian estimation."
  } else {
    ""
  }
  exp_regimes <- detect_exploratory_regimes(sp_hdr)
  exp_note <- if (length(exp_regimes) > 0L) {
    sprintf(
      " Exploratory regime (%s): nominal coverage and null Type I error are not pre-benchmarked across schedules; confirm at the chosen design via validate_design_point().",
      paste(exp_regimes, collapse = "; ")
    )
  } else {
    ""
  }
  sched_note <- if (!is.null(x$robust) && isTRUE(x$robust$schedules_differ)) {
    sprintf(
      " Scheduled trials differ across scenarios (at selected design: %s).",
      format_scenario_counts(x$robust$robust_trials_by_scenario)
    )
  } else {
    ""
  }
  gt_obj |>
    gt::tab_source_note(
      source_note = sprintf(
        "Target Status uses %.0f%% Wilson intervals (Morris et al., 2019) for %s: PASS (MC lower >= target), INDETERMINATE (straddles target), FAIL (MC upper < target). Target PASS does not confirm interval coverage or null Type I error; verify inferential validity with validate_design_point().%s%s%s%s Seed = %d.",
        sp_hdr$mc_conf_level * 100,
        crit_hdr,
        opt_sections_note,
        sing_note,
        exp_note,
        sched_note,
        x$search_meta$seed
      )
    )
}

#' @rdname as_gt
#' @export
as_gt.glmmplanr_stress_test <- function(x, ...) {
  df <- as.data.frame(x)
  gt::gt(df, rowname_col = "model") |>
    gt::tab_header(
      title = "GLMM Specification Mismatch Stress-Test",
      subtitle = sprintf("Paired Monte Carlo evaluation under Maximal DGP (N_subj = %d, N_item = %d)", df$n_subj[1L], df$n_item[1L])
    ) |>
    gt::cols_hide(columns = c("n_subj", "n_item")) |>
    gt::cols_label(
      power = "Eff. Power (H1)",
      conditional_power = "Cond. Power",
      valid_rate = "Valid",
      type1_error = "Type I Error (H0)",
      ci_width_mean = "CI (logit)",
      prob_ci_width_mean_pp = "CI (pp)",
      coverage = "Coverage",
      singular_rate = "Singular",
      rmse = "RMSE",
      type_m_error = "Type M"
    ) |>
    gt::fmt_percent(
      columns = c("power", "conditional_power", "valid_rate", "type1_error", "coverage", "singular_rate"),
      decimals = 1
    ) |>
    gt::fmt_number(columns = c("ci_width_mean", "rmse"), decimals = 3) |>
    gt::fmt_number(columns = c("prob_ci_width_mean_pp", "type_m_error"), decimals = 2) |>
    gt::sub_missing(missing_text = "-")
}

style_plot_annotations <- function(p, show_annotations) {
  if (!is.logical(show_annotations) || length(show_annotations) != 1L || is.na(show_annotations)) {
    cli::cli_abort("{.arg show_annotations} must be TRUE or FALSE.")
  }
  if (show_annotations) return(p)
  if (inherits(p, "plotly")) return(plotly::layout(p, title = list(text = "")))
  p + ggplot2::labs(title = NULL, subtitle = NULL, caption = NULL)
}

#' Plot Sample Size Trade-offs and Scenario Comparisons
#'
#' @param x A `"glmmplanr_result"` object.
#' @param type Visualization type:
#'   * `"contour"` (default): 2D `(N_subj x N_item)` response surface with
#'     iso-power/precision contours, non-dominated Pareto step frontier, and
#'     recommended design.
#'   * `"pareto"`: Non-dominated `(N_subj vs N_item)` Pareto trade-off curve
#'     showing the exact minimum participants required at each item count.
#'   * `"surface3d"`: Interactive `plotly` 3D response surface.
#'   * `"scenarios"`: Compare verified designs across the scenarios supplied to
#'     [optimize_sample_size()], showing detection with its stored Monte Carlo
#'     interval and mean interval width in separate panels. Absolute probability
#'     targets show width attainment in place of detection. Shapes give the
#'     stored decision under the plan's criterion.
#' @param metric Surface metric for `"contour"` and `"surface3d"`: `"power"`
#'   (default for contrasts and logit targets), `"ci_width"` (default for
#'   absolute probability targets; on `x$spec$scale`), `"prob_ci_width_pp"`
#'   (95% CI width in percentage points), or `"aipe_prob"` (width-attainment
#'   probability when `ci_width` is set).
#' @param show_pareto_frontier Logical; overlay the non-dominated `(N_subj, N_item)`
#'   Pareto frontier on `"contour"` plots (default `TRUE`).
#' @param ... Additional arguments (unused).
#' @param candidates For `"scenarios"`, an optional data frame with `n_subj`
#'   and `n_item` columns giving the verified design pairs to display. By
#'   default, show the designs selected under the standard scenario and
#'   across all scenarios, with one row of panels if the selections coincide.
#'   Increase the saved figure height when showing many designs.
#' @param show_annotations Logical; show automatic titles, subtitles and
#'   captions (default `FALSE`). Axis labels, legends, facets, and data labels
#'   remain visible.
#' @details Static figures use [ggplot2::theme_bw()]; add [ggplot2::labs()]
#'   to supply your own title or caption. The 3D plot uses a white background
#'   and can be edited with [plotly::layout()].
#' @return A `ggplot` object or `plotly` htmlwidget.
#' @examples
#' # Create the five-replication demonstration plan from the optimizer help.
#' example("optimize_sample_size", package = "glmmplanr", local = environment(), echo = FALSE)
#' plot(plan, type = "pareto")
#' @export
plot.glmmplanr_result <- function(
  x,
  type = c("contour", "surface3d", "pareto", "scenarios"),
  metric = c("power", "ci_width", "prob_ci_width_pp", "aipe_prob"),
  show_pareto_frontier = TRUE,
  ...,
  candidates = NULL,
  show_annotations = FALSE
) {
  type <- match.arg(type)
  metric_missing <- missing(metric)
  if (type == "scenarios") {
    p <- plot_scenario_comparison(x, candidates)
  } else {
    has_null <- can_define_null_target(x$spec)
    if (metric_missing) {
      metric <- if (has_null) "power" else "ci_width"
    } else {
      metric <- match.arg(metric)
      if (metric == "power" && !has_null) {
        cli::cli_abort("Detection power is undefined when the target has no attainable zero null; use {.code metric = 'ci_width'}, {.code 'prob_ci_width_pp'}, or {.code 'aipe_prob'}.")
      }
      if (metric == "aipe_prob" && is.null(x$spec$ci_width)) {
        cli::cli_abort("{.code metric = 'aipe_prob'} requires a specification with a finite {.arg ci_width} bound.")
      }
    }
    p <- switch(type,
      contour = plot_contour_2d(x, metric = metric, show_pareto_frontier = show_pareto_frontier),
      surface3d = plot_surface_3d(x, metric = metric),
      pareto = plot_pareto_curve(x))
  }
  style_plot_annotations(p, show_annotations)
}

plot_scenario_comparison <- function(x, candidates) {
  robust <- x$robust
  ev <- robust$scenario_evaluations
  if (is.null(ev) || !nrow(ev)) {
    cli::cli_abort("Scenario plots require a plan computed with multiple {.arg scenarios}.")
  }
  coord <- c("n_subj", "n_item")
  if (is.null(candidates)) {
    candidates <- rbind(robust$optimal_standard[coord], robust$optimal_robust[coord])
  }
  if (!is.data.frame(candidates) || !all(coord %in% names(candidates)) || !nrow(candidates)) {
    cli::cli_abort("Supply verified designs in a data frame with {.field n_subj} and {.field n_item} columns.")
  }
  candidates <- unique(candidates[coord])
  keys <- with(ev, paste(n_subj, n_item))
  chosen <- with(candidates, paste(n_subj, n_item))
  if (!all(chosen %in% keys)) cli::cli_abort("Every candidate must have verified scenario results.")
  ev <- ev[keys %in% chosen, , drop = FALSE]
  label <- function(d) sprintf("%s participants\n%s items in bank", d$n_subj, d$n_item)
  designs <- factor(label(ev), levels = label(candidates))
  presentation <- vapply(seq_len(nrow(ev)), function(i) {
    des <- robust$specs[[as.character(ev$scenario[i])]]$dgp$design
    if (is.null(des)) {
      return(sprintf("%g scheduled trials per participant", ev$n_trials_scheduled[i] / ev$n_subj[i]))
    }
    items <- if (is.null(des$items_per_subject)) ev$n_item[i] else des$items_per_subject
    sprintf("%s items per participant", items)
  }, character(1L))
  scenario <- paste(ev$scenario, presentation, sep = "\n")
  ordered <- order(ev$scenario != robust$standard_scenario, match(ev$scenario, unique(ev$scenario)))
  logit <- identical(x$spec$scale, "logit")
  has_null <- can_define_null_target(x$spec)
  rate <- if (has_null) "power" else "aipe_prob"
  mc_prefix <- if (has_null) "power" else "aipe"
  panels <- c(if (has_null) "Detection (%)" else "Width attainment (%)",
              if (logit) "Mean interval width (logits)" else "Mean interval width (pp)")
  width <- if (logit) ev$ci_width_mean else ev$prob_ci_width_mean_pp
  points <- data.frame(design = designs, scenario = factor(scenario, levels = rev(unique(scenario[ordered]))),
    decision = ev$decision_status, panel = factor(rep(panels, each = nrow(ev)), levels = panels),
    value = c(100 * ev[[rate]], width),
    lo = c(100 * ev[[paste0(mc_prefix, "_mc_lo")]], rep(NA_real_, nrow(ev))),
    hi = c(100 * ev[[paste0(mc_prefix, "_mc_hi")]], rep(NA_real_, nrow(ev))))
  ranges <- points[is.finite(points$value) & is.finite(points$lo) & is.finite(points$hi), ]
  p <- ggplot2::ggplot(points, ggplot2::aes(x = .data$value, y = .data$scenario))
  reference <- if (!has_null) x$spec$ci_assurance else
    if (x$spec$criterion %in% c("power", "both")) x$spec$power_min else NA_real_
  if (is.finite(reference)) {
    p <- p + ggplot2::geom_vline(data = data.frame(panel = panels[1L], target = 100 * reference),
      ggplot2::aes(xintercept = .data$target), linetype = "dashed", color = "grey55")
  }
  rate_sample <- if (x$spec$failure_policy == "effective") "all simulated studies" else "valid fits"
  p + ggplot2::geom_segment(data = ranges,
      ggplot2::aes(x = .data$lo, xend = .data$hi, yend = .data$scenario), color = "grey40") +
    ggplot2::geom_point(ggplot2::aes(shape = .data$decision), size = 2.8, na.rm = TRUE) +
    ggplot2::geom_text(data = points[!is.finite(points$value), ], x = 0,
      label = "Unavailable", hjust = 0, size = 3) +
    ggplot2::facet_grid(design ~ panel, scales = "free", space = "free_y", switch = "y") +
    ggplot2::expand_limits(x = 0) +
    ggplot2::scale_shape_manual(values = c(pass = 16, indeterminate = 17, fail = 4)) +
    ggplot2::labs(title = "Verified designs across scenarios", x = NULL, y = NULL, shape = "Decision",
      caption = sprintf(paste("Bars show stored Monte Carlo intervals for rates among %s.",
                              "Mean widths use valid fits. Shapes show the plan's criterion decision.", sep = "\n"),
                        rate_sample)) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(strip.placement = "outside", strip.text.y.left = ggplot2::element_text(angle = 0),
      legend.position = "bottom", plot.title.position = "plot",
      plot.caption.position = "plot", plot.caption = ggplot2::element_text(hjust = 0))
}

#' Plot Fixed-design Validation Results
#'
#' @param x A result from [validate_design_point()].
#' @param type `"rates"` (default) shows planned power, null Type I error, and
#'   interval coverage with Wilson Monte Carlo intervals. For absolute
#'   probability targets it displays coverage and any planned width attainment.
#'   `"patterns"` shows the share of diagnosed studies with at least one
#'   all-correct, all-incorrect, or unobserved participant or item. For the
#'   two extreme scores it also reports the largest observation count seen in
#'   a unit with that pattern. Use [response_patterns()] for counts by unit
#'   in a single dataset.
#' @param ... Additional arguments (unused).
#' @inheritParams plot.glmmplanr_result
#' @return An editable `ggplot` object. Dashed lines in the rate plot show
#'   nominal or planning targets. These lines are reference values.
#'   The default theme is [ggplot2::theme_bw()]. Add [ggplot2::labs()] to
#'   supply your own title or caption.
#' @examples
#' des <- design_crossed(list(condition = c("a", "b")))
#' dgp <- calibrate_dgp(des, c(a = 0.4, b = 0.6), subj_sd = 0.3, item_sd = 0.4)
#' spec <- specify_target(dgp, "condition", criterion = "power", engine = "bayes_laplace")
#' val <- validate_design_point(spec, 16, 8, nsim = 5L)
#' plot(val)
#' plot(val, type = "patterns")
#' @export
plot.glmmplanr_validation <- function(x, type = c("rates", "patterns"), ...,
                                      show_annotations = FALSE) {
  type <- match.arg(type)
  settings <- attr(x, "validation_settings")
  cols <- c("condition", "nsim", "n_valid", "detection", "detection_lo", "detection_hi",
            "coverage", "coverage_lo", "coverage_hi", "width_attainment", pattern_columns)
  if (is.null(settings) || !nrow(x) || !all(cols %in% names(x)) || anyDuplicated(x$condition)) {
    cli::cli_abort("Plot the complete result from {.fn validate_design_point}.")
  }
  subtitle <- sprintf("%s participants / %s items in bank / %s items per participant",
    settings$n_subj, settings$n_item, settings$items_per_subject)
  if (type == "patterns") return(style_plot_annotations(plot_validation_patterns(x, subtitle), show_annotations))

  condition <- ifelse(x$condition == "null", "Null", "Planned")
  studies <- sprintf("%s\n%s/%s valid fits", condition, x$n_valid, x$nsim)
  detection_denom <- if (settings$failure_policy == "effective") x$nsim else x$n_valid
  make_panel <- function(panel, rate, lo, hi, target) {
    data.frame(panel, studies, rate, lo, hi, target)
  }
  coverage <- make_panel("Interval coverage", x$coverage, x$coverage_lo, x$coverage_hi, 1 - settings$alpha)
  width_bounds <- function(df) {
    if (all(c("width_attainment_lo", "width_attainment_hi") %in% names(df))) {
      list(lo = df$width_attainment_lo, hi = df$width_attainment_hi)
    } else {
      denom <- if (settings$failure_policy == "effective") df$nsim else df$n_valid
      ci <- vapply(seq_len(nrow(df)), function(i) {
        k_i <- if (is.finite(df$width_attainment[i]) && denom[i] > 0L) {
          round(df$width_attainment[i] * denom[i])
        } else {
          NA_integer_
        }
        wilson_mc_ci(k_i, denom[i], settings$mc_conf_level)
      }, numeric(2L))
      list(lo = unname(ci["lo", ]), hi = unname(ci["hi", ]))
    }
  }
  if (settings$has_null) {
    power_target <- if (settings$criterion %in% c("power", "both")) settings$power_min else NA_real_
    row_zero <- vapply(seq_len(nrow(x)), function(i) {
      cnm <- as.character(x$condition[i])
      if (!is.null(settings$zero_truth_verified) && (cnm %in% names(settings$zero_truth_verified))) {
        return(isTRUE(settings$zero_truth_verified[[cnm]]))
      }
      if (identical(cnm, "null") && !is.null(settings$null_verified)) {
        return(isTRUE(settings$null_verified))
      }
      is.finite(x$true_target[i]) && (x$true_target[i] == 0)
    }, logical(1L))
    det_panel <- ifelse(
      row_zero,
      "Type I error",
      ifelse(x$condition == "null", "Null-input detection", "Power")
    )
    det_target <- ifelse(
      row_zero,
      settings$alpha,
      ifelse(
        x$condition == "null" | is.na(x$true_target),
        NA_real_,
        power_target
      )
    )
    detection <- make_panel(det_panel, x$detection, x$detection_lo, x$detection_hi, det_target)
    points <- rbind(detection, coverage)
    p_idx <- which(x$condition == "planned")
    if (length(p_idx) > 0L) {
      xp <- x[p_idx, , drop = FALSE]
      sp_p <- studies[p_idx]
      if (!is.null(settings$ci_width) && settings$criterion %in% c("ci_width", "both")) {
        wb <- width_bounds(xp)
        points <- rbind(points, data.frame(
          panel = "Width attainment", studies = sp_p,
          rate = xp$width_attainment, lo = wb$lo, hi = wb$hi,
          target = settings$ci_assurance
        ))
      } else if (!is.null(settings$ci_width) && settings$criterion == "joint" &&
                 all(c("joint_attainment", "joint_attainment_lo", "joint_attainment_hi") %in% names(xp))) {
        points <- rbind(points, data.frame(
          panel = "Joint attainment", studies = sp_p,
          rate = xp$joint_attainment, lo = xp$joint_attainment_lo, hi = xp$joint_attainment_hi,
          target = settings$power_min
        ))
      } else if (settings$criterion == "rope" &&
                 all(c("rope_attainment", "rope_attainment_lo", "rope_attainment_hi") %in% names(xp))) {
        points <- rbind(points, data.frame(
          panel = "ROPE exclusion", studies = sp_p,
          rate = xp$rope_attainment, lo = xp$rope_attainment_lo, hi = xp$rope_attainment_hi,
          target = settings$power_min
        ))
      }
    }
  } else {
    points <- coverage
    if (!is.null(settings$ci_width)) {
      wb <- width_bounds(x)
      points <- rbind(points, make_panel("Width attainment", x$width_attainment,
        wb$lo, wb$hi, settings$ci_assurance))
    }
    if (settings$criterion == "rope" &&
        all(c("rope_attainment", "rope_attainment_lo", "rope_attainment_hi") %in% names(x))) {
      points <- rbind(points, make_panel("ROPE exclusion", x$rope_attainment,
        x$rope_attainment_lo, x$rope_attainment_hi, settings$power_min))
    }
  }
  points$panel <- factor(
    points$panel,
    levels = c("Power", "Type I error", "Null-input detection", "Interval coverage",
               "Width attainment", "Joint attainment", "ROPE exclusion")
  )
  points$studies <- factor(points$studies, levels = rev(unique(studies[order(x$condition != "null")])))
  ranges <- points[is.finite(points$rate) & is.finite(points$lo) & is.finite(points$hi), ]
  targets <- unique(points[is.finite(points$target), c("panel", "target")])
  rate_denom <- if (settings$failure_policy == "effective") "all scheduled studies" else "valid fits"
  rate_label <- if (settings$has_null) {
    "Detection"
  } else if (settings$criterion == "rope") {
    "ROPE exclusion"
  } else {
    "Width attainment"
  }
  caption <- sprintf("%g%% Wilson Monte Carlo intervals. Coverage uses valid fits.", 100 * settings$mc_conf_level)
  if (settings$has_null || !is.null(settings$ci_width) || settings$criterion == "rope") {
    caption <- paste(caption, sprintf("%s uses %s.", rate_label, rate_denom), sep = "\n")
  }
  p <- ggplot2::ggplot(points, ggplot2::aes(x = .data$rate, y = .data$studies)) +
    ggplot2::geom_vline(data = targets, ggplot2::aes(xintercept = .data$target),
      linetype = "dashed", color = "grey55") +
    ggplot2::geom_segment(data = ranges,
      ggplot2::aes(x = .data$lo, xend = .data$hi, yend = .data$studies), color = "#2563eb") +
    ggplot2::geom_point(size = 2.6, color = "#2563eb", na.rm = TRUE) +
    ggplot2::geom_text(data = points[!is.finite(points$rate), ], x = 0.5, label = "Unavailable", size = 3) +
    ggplot2::facet_wrap(~panel, ncol = 1, scales = "free_y") +
    ggplot2::scale_x_continuous(limits = c(0, 1), labels = function(v) paste0(100 * v, "%")) +
    ggplot2::labs(title = "Fixed-design validation", subtitle = subtitle, x = "Share of studies", y = NULL,
      caption = paste(caption, "Dashed lines show nominal or planning targets.")) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(plot.title.position = "plot", plot.caption.position = "plot",
      plot.caption = ggplot2::element_text(hjust = 0))
  style_plot_annotations(p, show_annotations)
}

plot_validation_patterns <- function(x, subtitle) {
  labels <- rep(c("All correct", "All incorrect", "Unobserved"), 2L)
  units <- rep(c("Participants", "Items"), each = 3L)
  points <- do.call(rbind, lapply(seq_len(nrow(x)), function(i) {
    diagnosed <- x$pattern_nsim[i]
    counts <- unlist(x[i, paste0("n_rep_", unit_pattern_names)], use.names = FALSE)
    max_obs <- rep(NA_integer_, length(labels))
    max_obs[labels != "Unobserved"] <- unlist(x[i, observed_max_columns], use.names = FALSE)
    details <- paste0(counts, "/", diagnosed, " studies",
      ifelse(!is.na(max_obs) & counts > 0L, paste0(" / max ", max_obs, " responses"), ""))
    data.frame(condition = if (x$condition[i] == "null") "Null" else "Planned",
      pattern = paste0(units, " - ", labels, "\n", details),
      rate = if (diagnosed > 0L) counts / diagnosed else NA_real_)
  }))
  key <- paste(points$condition, points$pattern)
  points$row <- factor(key, levels = rev(key))
  ggplot2::ggplot(points, ggplot2::aes(x = .data$rate, y = .data$row)) +
    ggplot2::geom_point(size = 2.6, color = "#2563eb", na.rm = TRUE) +
    ggplot2::facet_wrap(~condition, ncol = 1, scales = "free_y") +
    ggplot2::scale_x_continuous(limits = c(0, 1), labels = function(v) paste0(100 * v, "%")) +
    ggplot2::scale_y_discrete(labels = stats::setNames(points$pattern, key)) +
    ggplot2::labs(title = "Response patterns in validation studies", subtitle = subtitle,
      x = "Share of diagnosed studies with at least one such unit", y = NULL,
      caption = "Max responses is the largest observation count in an extreme-scoring unit.") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(plot.title.position = "plot", plot.caption.position = "plot",
      plot.caption = ggplot2::element_text(hjust = 0))
}

axis_tile_bounds <- function(vals) {
  u <- sort(unique(vals))
  n <- length(u)
  if (n == 1L) {
    lo <- u - 0.5
    hi <- u + 0.5
  } else {
    mids <- 0.5 * (u[-n] + u[-1L])
    lo <- c(u[1L] - 0.5 * (u[2L] - u[1L]), mids)
    hi <- c(mids, u[n] + 0.5 * (u[n] - u[n - 1L]))
  }
  idx <- match(vals, u)
  list(lo = lo[idx], hi = hi[idx], breaks = u, range = c(lo[1L], hi[n]))
}

is_regular_grid_axis <- function(vals) {
  u <- sort(unique(vals))
  if (length(u) <= 2L) return(TRUE)
  d <- diff(u)
  max(abs(d - d[1L])) <= 1e-8 * max(1.0, abs(d[1L]))
}

plot_contour_2d <- function(x, metric = "power", show_pareto_frontier = TRUE) {
  surf <- x$surface
  opt <- x$optimal
  ev <- x$evaluated
  pf <- x$pareto_frontier

  fill_col <- switch(
    metric,
    power = "pred_power",
    ci_width = "pred_ci_width",
    prob_ci_width_pp = "pred_prob_ci_width_pp",
    aipe_prob = "pred_aipe_prob"
  )
  fill_lbl <- switch(
    metric,
    power = "Predicted\nPower",
    ci_width = sprintf("CI\n(%s)", x$spec$scale),
    prob_ci_width_pp = "CI\n(pp)",
    aipe_prob = "Predicted\nWidth Attainment"
  )

  ux <- sort(unique(surf$n_subj))
  uy <- sort(unique(surf$n_item))
  regular_xy <- is_regular_grid_axis(surf$n_subj) && is_regular_grid_axis(surf$n_item)
  use_tiles <- (!regular_xy) || (length(ux) * length(uy) <= 36L)

  if (use_tiles) {
    bx <- axis_tile_bounds(surf$n_subj)
    by <- axis_tile_bounds(surf$n_item)
    surf$xmin <- bx$lo
    surf$xmax <- bx$hi
    surf$ymin <- by$lo
    surf$ymax <- by$hi
    p <- ggplot2::ggplot(
      surf,
      ggplot2::aes(x = .data$n_subj, y = .data$n_item)
    ) +
      ggplot2::geom_rect(
        ggplot2::aes(
          xmin = .data$xmin,
          xmax = .data$xmax,
          ymin = .data$ymin,
          ymax = .data$ymax,
          fill = .data[[fill_col]]
        ),
        color = "white",
        linewidth = 0.35,
        alpha = 0.88
      )
    x_lim <- bx$range
    y_lim <- by$range
    coord_expand <- FALSE
  } else {
    p <- ggplot2::ggplot(
      surf,
      ggplot2::aes(x = .data$n_subj, y = .data$n_item)
    ) +
      ggplot2::geom_raster(
        ggplot2::aes(fill = .data[[fill_col]]),
        interpolate = TRUE,
        alpha = 0.88
      )
    x_lim <- range(surf$n_subj)
    y_lim <- range(surf$n_item)
    coord_expand <- TRUE
  }

  if (metric %in% c("power", "aipe_prob")) {
    p <- p + ggplot2::scale_fill_viridis_c(
      name = fill_lbl,
      limits = c(0, 1),
      labels = function(v) paste0(round(v * 100), "%"),
      option = "viridis"
    )
  } else {
    p <- p + ggplot2::scale_fill_viridis_c(
      name = fill_lbl,
      direction = -1,
      option = "magma"
    )
  }

  if (length(unique(surf$n_subj)) > 1L && length(unique(surf$n_item)) > 1L) {
    if (any(is.finite(surf$pred_power))) {
      p_rng <- range(surf$pred_power, na.rm = TRUE)
      gen_brks <- setdiff(seq(0.20, 0.95, by = 0.10), x$spec$power_min)
      gen_brks <- gen_brks[gen_brks > p_rng[1L] & gen_brks < p_rng[2L]]
      if (length(gen_brks) > 0L) {
        p <- p + ggplot2::geom_contour(
          ggplot2::aes(z = .data$pred_power),
          breaks = gen_brks,
          color = "white",
          alpha = 0.55,
          linewidth = 0.45
        )
      }

      if (x$spec$criterion %in% c("power", "both") &&
          x$spec$power_min > p_rng[1L] && x$spec$power_min < p_rng[2L]) {
        p <- p + ggplot2::geom_contour(
          ggplot2::aes(z = .data$pred_power),
          breaks = x$spec$power_min,
          color = "#ffffff",
          linewidth = 1.25
        )
      }
    }

    if (x$spec$criterion %in% c("joint", "rope")) {
      column <- if (x$spec$criterion == "joint") "pred_joint_prob" else "pred_rope_prob"
      values <- surf[[column]]
      values <- values[is.finite(values)]
      if (length(values) && min(values) < x$spec$power_min && max(values) > x$spec$power_min) {
        p <- p + ggplot2::geom_contour(ggplot2::aes(z = .data[[column]]),
          breaks = x$spec$power_min, color = "#fde047", linetype = "dashed")
      }
    }
    if (x$spec$criterion %in% c("ci_width", "both") && !is.null(x$spec$ci_width) &&
        any(is.finite(surf$pred_aipe_prob))) {
      a_rng <- range(surf$pred_aipe_prob, na.rm = TRUE)
      if (x$spec$ci_assurance > a_rng[1L] && x$spec$ci_assurance < a_rng[2L]) {
        p <- p + ggplot2::geom_contour(
          ggplot2::aes(z = .data$pred_aipe_prob),
          breaks = x$spec$ci_assurance,
          color = "#fde047",
          linetype = "dashed",
          linewidth = 1.15
        )
      }
    }
  }

  # Overlay non-dominated Pareto-efficient `(N_subj, N_item)` frontier
  if (isTRUE(show_pareto_frontier) && !is.null(pf) && nrow(pf) >= 1L) {
    pf_ord <- pf[order(pf$n_subj, -pf$n_item), , drop = FALSE]
    if (nrow(pf_ord) >= 2L) {
      p <- p + ggplot2::geom_line(
        data = pf_ord,
        ggplot2::aes(x = .data$n_subj, y = .data$n_item),
        inherit.aes = FALSE,
        color = "#38bdf8",
        linewidth = 1.0,
        linetype = "twodash"
      )
    }
    p <- p + ggplot2::geom_point(
      data = pf_ord,
      ggplot2::aes(x = .data$n_subj, y = .data$n_item),
      inherit.aes = FALSE,
      shape = 21,
      size = 3.0,
      fill = "#38bdf8",
      color = "#0f172a"
    )
  }

  has_null <- can_define_null_target(x$spec)
  att_sub <- if (has_null && is.finite(opt$power)) {
    sprintf("Power = %.1f%% [%.1f, %.1f]", opt$power * 100, opt$power_mc_lo * 100, opt$power_mc_hi * 100)
  } else if (!is.na(opt$aipe_prob)) {
    sprintf("Width attainment = %.1f%% [%.1f, %.1f]", opt$aipe_prob * 100, opt$aipe_mc_lo * 100, opt$aipe_mc_hi * 100)
  } else if (!is.na(opt$rope_excl_prob)) {
    sprintf("ROPE excl. = %.1f%% [%.1f, %.1f]", opt$rope_excl_prob * 100, opt$rope_mc_lo * 100, opt$rope_mc_hi * 100)
  } else {
    sprintf("Valid = %.1f%%", opt$valid_rate * 100)
  }

  p <- p +
    ggplot2::geom_point(
      data = ev,
      ggplot2::aes(x = .data$n_subj, y = .data$n_item),
      inherit.aes = FALSE,
      shape = 21,
      size = 2.0,
      fill = "#111827",
      color = "white",
      alpha = 0.70
    ) +
    ggplot2::geom_point(
      data = opt,
      ggplot2::aes(x = .data$n_subj, y = .data$n_item),
      inherit.aes = FALSE,
      shape = 23,
      size = 4.6,
      stroke = 1.2,
      fill = "#ef4444",
      color = "white"
    ) +
    ggplot2::coord_cartesian(
      xlim = x_lim,
      ylim = y_lim,
      expand = coord_expand,
      clip = "off"
    ) +
    ggplot2::labs(
      title = sprintf(
        "2D Sample Size Trade-off Surface: %s (%s)",
        x$spec$effect_input, toupper(opt$decision_status)
      ),
      subtitle = sprintf(
        "Selected candidate (%s): N_subj = %d, N_item = %d | %s | CI = %.2f pp",
        x$prefer, opt$n_subj, opt$n_item,
        att_sub, opt$prob_ci_width_mean_pp
      ),
      x = sprintf("Number of Subjects (N_subj, multiples of %d)", x$spec$dgp$design$k_subj),
      y = sprintf("Number of Items (N_item, multiples of %d)", x$spec$dgp$design$k_item),
      caption = paste0(
        switch(x$spec$criterion, power = "Solid white = target power",
          both = "Solid white = target power | Dashed yellow = target width assurance",
          ci_width = "Dashed yellow = target width assurance",
          joint = "Dashed yellow = detection and width in the same study",
          rope = "Dashed yellow = target ROPE exclusion"), "\n",
        "Cyan nodes = non-dominated Pareto designs | Red diamond = selected candidate"
      )
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.caption = ggplot2::element_text(hjust = 0, size = 8.5, color = "#4b5563")
    )

  if (length(ux) <= 10L) {
    p <- p + ggplot2::scale_x_continuous(breaks = ux)
  }
  if (length(uy) <= 10L) {
    p <- p + ggplot2::scale_y_continuous(breaks = uy)
  }

  p
}

plot_surface_3d <- function(x, metric = "power") {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    cli::cli_abort(
      "Package {.pkg plotly} is required for {.code plot(..., type = 'surface3d')}. Please install {.pkg plotly}."
    )
  }
  surf <- x$surface
  opt <- x$optimal
  s_vals <- sort(unique(surf$n_subj))
  i_vals <- sort(unique(surf$n_item))

  z_col <- switch(
    metric,
    power = "pred_power",
    ci_width = "pred_ci_width",
    prob_ci_width_pp = "pred_prob_ci_width_pp",
    aipe_prob = "pred_aipe_prob"
  )
  z_title <- switch(
    metric,
    power = "Power",
    ci_width = sprintf("CI Width (%s)", x$spec$scale),
    prob_ci_width_pp = "CI Width (pp)",
    aipe_prob = "Width Attainment"
  )

  Z_mat <- matrix(NA_real_, nrow = length(i_vals), ncol = length(s_vals))
  Tr_mat <- matrix(NA_integer_, nrow = length(i_vals), ncol = length(s_vals))
  W_pp_mat <- matrix(NA_real_, nrow = length(i_vals), ncol = length(s_vals))
  P_mat <- matrix(NA_real_, nrow = length(i_vals), ncol = length(s_vals))
  A_mat <- matrix(NA_real_, nrow = length(i_vals), ncol = length(s_vals))

  r_ix <- match(surf$n_item, i_vals)
  c_ix <- match(surf$n_subj, s_vals)
  Z_mat[cbind(r_ix, c_ix)] <- surf[[z_col]]
  Tr_mat[cbind(r_ix, c_ix)] <- surf$n_trials_scheduled
  W_pp_mat[cbind(r_ix, c_ix)] <- surf$pred_prob_ci_width_pp
  P_mat[cbind(r_ix, c_ix)] <- surf$pred_power
  A_mat[cbind(r_ix, c_ix)] <- surf$pred_aipe_prob

  key_mat <- outer(i_vals, s_vals, function(ii, ss) paste(ss, ii))
  idx_mat <- match(key_mat, paste(surf$n_subj, surf$n_item))
  has_null <- can_define_null_target(x$spec)
  rate_line <- if (has_null) {
    sprintf("Power: %.1f%%<br>", as.numeric(P_mat) * 100)
  } else if (any(is.finite(A_mat))) {
    sprintf("Width Attainment: %.1f%%<br>", as.numeric(A_mat) * 100)
  } else {
    rep("", length(idx_mat))
  }

  hover_txt <- matrix(
    sprintf(
      "N_subj: %d<br>N_item: %d<br>Sched. Trials: %s<br>%sCI Width: %.2f pp",
      surf$n_subj[idx_mat],
      surf$n_item[idx_mat],
      format(as.numeric(Tr_mat), big.mark = ","),
      rate_line,
      as.numeric(W_pp_mat)
    ),
    nrow = length(i_vals),
    ncol = length(s_vals)
  )

  z_opt <- switch(
    metric,
    power = opt$power,
    ci_width = opt$ci_width_mean,
    prob_ci_width_pp = opt$prob_ci_width_mean_pp,
    aipe_prob = opt$aipe_prob
  )
  thr_val <- if (metric == "power") {
    x$spec$power_min
  } else if (metric == "aipe_prob") {
    x$spec$ci_assurance
  } else if (metric == "ci_width") {
    x$spec$ci_width
  } else if (x$spec$scale == "pp") {
    x$spec$ci_width
  } else {
    NULL
  }

  fig <- plotly::plot_ly() |>
    plotly::add_surface(
      x = s_vals,
      y = i_vals,
      z = Z_mat,
      text = hover_txt,
      hoverinfo = "text",
      colorscale = if (metric %in% c("power", "aipe_prob")) "Viridis" else "Magma",
      opacity = 0.90,
      name = z_title
    )
  if (length(z_opt) == 1L && is.finite(z_opt)) {
    fig <- fig |>
      plotly::add_markers(
        x = opt$n_subj,
        y = opt$n_item,
        z = z_opt,
        marker = list(size = 7, color = "#ef4444", symbol = "diamond"),
        name = sprintf("Selected candidate (%d, %d)", opt$n_subj, opt$n_item)
      )
  }

  if (!is.null(thr_val) && is.finite(thr_val)) {
    Plane_mat <- matrix(thr_val, nrow = length(i_vals), ncol = length(s_vals))
    fig <- fig |>
      plotly::add_surface(
        x = s_vals,
        y = i_vals,
        z = Plane_mat,
        showscale = FALSE,
        opacity = 0.30,
        colorscale = list(c(0, "#fbbf24"), c(1, "#fbbf24")),
        hoverinfo = "skip",
        name = "Target Threshold"
      )
  }

  axis_style <- list(showbackground = TRUE, backgroundcolor = "white",
    gridcolor = "#EBEBEB", zerolinecolor = "#B3B3B3", showline = TRUE, linecolor = "black")
  fig |>
    plotly::layout(
      title = sprintf("Interactive 3D GLMM Response Surface (%s)", x$spec$effect_input),
      paper_bgcolor = "white", plot_bgcolor = "white", font = list(color = "black"),
      scene = list(
        bgcolor = "white",
        xaxis = c(list(title = "N_subj (Subjects)"), axis_style),
        yaxis = c(list(title = "N_item (Items)"), axis_style),
        zaxis = c(list(title = z_title), axis_style)
      )
    )
}

plot_pareto_curve <- function(x) {
  pf <- x$pareto_frontier[order(x$pareto_frontier$n_subj, -x$pareto_frontier$n_item), , drop = FALSE]
  pf$label <- sprintf("(%d, %d)", pf$n_subj, pf$n_item)
  opt <- x$optimal

  ggplot2::ggplot(pf, ggplot2::aes(x = .data$n_subj, y = .data$n_item)) +
    ggplot2::geom_line(color = "#2563eb", linewidth = 0.95) +
    ggplot2::geom_point(
      size = 3.2,
      shape = 21,
      fill = "#38bdf8",
      color = "#1e3a8a"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$label),
      vjust = -0.9,
      size = 3.1,
      color = "#1f2937"
    ) +
    ggplot2::geom_point(
      data = opt,
      ggplot2::aes(x = .data$n_subj, y = .data$n_item),
      inherit.aes = FALSE,
      shape = 23,
      size = 4.6,
      fill = "#ef4444",
      color = "white"
    ) +
    ggplot2::labs(
      title = if (nrow(pf)) "Non-dominated sample sizes" else "No confirmed Pareto designs",
      subtitle = sprintf(
        "Selected candidate (%s): N_subj = %d, N_item = %d (%s scheduled trials, MC Status = %s)",
        x$prefer, opt$n_subj, opt$n_item,
        format(opt$n_trials_scheduled, big.mark = ","), toupper(opt$decision_status)
      ),
      x = "Participants",
      y = "Items"
    ) +
    ggplot2::theme_bw(base_size = 11)
}

#' Write a Sample Size Justification
#'
#' Writes an English sample size justification describing the experimental design,
#' scheduled and expected post-attrition trials, the Pareto trade-off frontier,
#' percentage-point and conditional logit effects, and estimation failure accounting.
#' It reports Monte Carlo intervals (Morris et al., 2019) and available sensitivity
#' results for alternative scenarios and random-slope assumptions.
#'
#' For Rasch planning, the report gives the error limits, joint attainment,
#' number of evaluable studies, and confirmed designs.
#' @param x A result from [optimize_sample_size()], [plan_rasch_fixed_bank()],
#'   or [plan_rasch_precision()].
#' @param print_output Logical; if `TRUE` (default), prints the formatted text
#'   to the console and returns it invisibly.
#' @return Character string containing the methods report (invisibly when
#'   `print_output = TRUE`).
#' @examples
#' # Create the five-replication demonstration plan from the optimizer help.
#' example("optimize_sample_size", package = "glmmplanr", local = environment(), echo = FALSE)
#' report_plan(plan, print_output = FALSE)
#' @export
report_plan <- function(x, print_output = TRUE) {
  if (!rlang::is_bool(print_output)) {
    cli::cli_abort("{.arg print_output} must be TRUE or FALSE.")
  }

  if (inherits(x, "glmmplanr_rasch")) {
    txt <- rasch_report(x)
    if (!print_output) return(txt)
    cat(txt, "\n")
    return(invisible(txt))
  }
  if (!inherits(x, "glmmplanr_result")) {
    cli::cli_abort("{.arg x} must be a {.cls glmmplanr_result} or {.cls glmmplanr_rasch} object.")
  }
  opt <- x$optimal
  sp <- x$spec
  if (!is.null(x$robust$specs) && !is.null(opt$scenario)) sp <- x$robust$specs[[opt$scenario]]
  dgp <- sp$dgp
  des <- dgp$design
  unit_lbl <- switch(sp$scale, logit = "logit", pp = "percentage points", prob = "probability units")

  one_cell <- des$n_cells == 1L
  no_null_target <- !can_define_null_target(sp)
  p_kind <- prob_target_kind(sp)
  is_abs_prob <- identical(p_kind, "absolute")
  prob_label <- switch(
    p_kind,
    absolute = "accuracy",
    difference = "accuracy difference",
    weighted_contrast = "weighted accuracy contrast"
  )
  details <- character()
  if (!is.null(x$search_meta$n_grid_verified)) {
    total <- x$search_meta$n_grid_total
    verified <- x$search_meta$n_grid_verified
    details <- sprintf("Final verification covered %d of %d lattice points.", verified, total)
  }
  if (sp$criterion %in% c("both", "joint") && !is.null(opt$joint_prob)) {
    rule <- if (sp$criterion == "both") {
      sprintf("Separate thresholds are %.0f%% for detection and %.0f%% for width attainment.", 100 * sp$power_min, 100 * sp$ci_assurance)
    } else {
      sprintf("The same-replication target is %.0f%%.", 100 * sp$power_min)
    }
    details <- c(details, rule, sprintf(
      "Width limit %.2f %s, width attainment %.1f%%, same-replication detection and width attainment %.1f%% (MC interval %.1f to %.1f%%).",
      sp$ci_width, unit_lbl, 100 * opt$aipe_prob, 100 * opt$joint_prob,
      100 * opt$joint_mc_lo, 100 * opt$joint_mc_hi))
  }
  details <- c(details, analysis_description(sp), sprintf(
    "The target interval level was %g%%.", 100 * (1 - sp$alpha)))
  if (identical(x$search_meta$mc_adjustment, "bonferroni") &&
      !is.null(x$search_meta$mc_conf_level_effective)) {
    details <- c(details, sprintf(
      "Stage-2 Wilson intervals used Bonferroni multiplicity adjustment across %d verified decisions (initial Stage-2 level %.2f%%).",
      as.integer(x$search_meta$n_stage2_decisions),
      100 * x$search_meta$mc_conf_level_effective
    ))
  } else if (isTRUE(x$search_meta$n_stage2_decisions > 1L)) {
    details <- c(details, sprintf(
      "Stage-2 Wilson intervals are pointwise across %d verified decisions without simultaneous multiplicity adjustment.",
      as.integer(x$search_meta$n_stage2_decisions)
    ))
  }
  if (isTRUE(x$search_meta$adaptive_pooled)) {
    details <- c(details,
      "Exploratory adaptive top-up pooling was applied; confirm the selected design on a fresh independent seed with validate_design_point().")
  } else if (isTRUE(x$search_meta$topup_rounds > 0L) &&
             identical(x$search_meta$topup_policy, "independent") &&
             !is.null(x$search_meta$topup_mc_conf_level) &&
             is.finite(x$search_meta$topup_mc_conf_level)) {
    seq_mode <- if (identical(x$search_meta$sequential_adjustment, "none")) {
      "without sequential multiplicity adjustment"
    } else if (identical(x$search_meta$sequential_adjustment, "geometric")) {
      "with geometric sequential alpha-spending"
    } else {
      "with a two-stage Bonferroni split"
    }
    details <- c(details, sprintf(
      "Borderline candidates underwent %d independent top-up round(s) %s (latest Wilson level %.2f%%).",
      as.integer(x$search_meta$topup_rounds),
      seq_mode,
      100 * x$search_meta$topup_mc_conf_level
    ))
  }
  has_null <- can_define_null_target(sp)
  cov_val <- if (sp$scale == "logit") opt$coverage else opt$prob_coverage
  val_info <- assess_inferential_validity(cov_val, opt$n_valid, sp$alpha, has_null = has_null)
  if (identical(val_info$status, "caution_low_coverage")) {
    details <- c(details, sprintf(
      "Target attainment status reflects only the specified power/precision targets and does not verify inferential validity. Caution: Stage-2 empirical coverage was %.1f%% (95%% Wilson CI %.1f%% to %.1f%%), below the nominal %.0f%% level; narrow intervals or high detection may reflect underestimated standard errors. Confirm interval coverage%s in independent simulations with validate_design_point() or an external estimator.",
      100 * val_info$cov, 100 * val_info$cov_lo, 100 * val_info$cov_hi, 100 * val_info$nominal,
      if (has_null) " and null Type I error" else ""
    ))
  } else {
    cov_summary <- if (is.finite(val_info$cov)) {
      sprintf("%.1f%% (95%% Wilson CI %.1f%% to %.1f%%)", 100 * val_info$cov, 100 * val_info$cov_lo, 100 * val_info$cov_hi)
    } else {
      "NA"
    }
    if (has_null) {
      details <- c(details, sprintf(
        "Target attainment status reflects only the specified power/precision targets and does not verify inferential validity (Stage-2 H1 coverage was %s; null Type I error was not evaluated). Confirm interval coverage and Type I error in independent simulations at the selected design with validate_design_point().",
        cov_summary
      ))
    } else {
      null_st_rep <- probability_contrast_null_status(sp$dgp, sp$contrast_weights, sp$scale, sp$prob_estimand)
      no_null_reason <- if (is_abs_prob) {
        "no zero-difference null applies to an absolute probability target"
      } else if (identical(null_st_rep$kind, "numerical_boundary")) {
        "a zero-contrast null lies beyond reliable double-precision calibration for these extreme relative weights"
      } else {
        "zero lies outside the attainable range of this weighted probability contrast"
      }
      details <- c(details, sprintf(
        "Target attainment status reflects only the specified precision target and does not verify inferential validity (Stage-2 coverage was %s; %s). Confirm interval coverage and width attainment in independent simulations at the selected design with validate_design_point().",
        cov_summary, no_null_reason
      ))
    }
  }
  if (!is.null(des$items_per_subject)) {
    details <- c(details, sprintf(
      "Each participant is assigned %2$d distinct items from a bank of %1$d (allocation seed %3$d); reported operating characteristics and validate_design_point() condition on this fixed schedule (vary item_seed across scenarios if item subsets will be re-randomized at execution).",
      opt$n_item, des$items_per_subject, des$item_seed))
  }
  if (dropout_rate(des) > 0 && !is.null(opt$n_subj_expected)) {
    details <- c(details, sprintf(
      "Recruit %d participants, with %.1f%% participant dropout and %.1f expected to supply observations.",
      opt$n_subj, 100 * dropout_rate(des), opt$n_subj_expected))
  }
  if (des$has_testlet) {
    details <- c(details, sprintf(
      "%d passages with %d items each. The generating participant-by-passage SD is %.2f logit; the analysis uses a shared passage effect.",
      opt$n_item %/% des$testlet_size, des$testlet_size, testlet_interaction_sd(dgp)))
  }
  exp_regimes <- detect_exploratory_regimes(sp)
  if (length(exp_regimes) > 0L) {
    details <- c(details, sprintf(
      "Scope note: This plan uses an exploratory design regime (%s) outside the package's pre-benchmarked complete-presentation regimes; confirm empirical interval coverage%s at the planned sample size via validate_design_point() (and an external estimator when feasible) before pre-registration.",
      paste(exp_regimes, collapse = "; "),
      if (has_null) " and null Type I error" else ""
    ))
  }
  bnd_info <- assess_boundary_sensitivity(x)
  if (!is.null(bnd_info)) {
    ind_clause <- if (bnd_info$n_indeterminate_nearby > 0L) {
      sprintf(" and %d nearby candidate(s) remained indeterminate", bnd_info$n_indeterminate_nearby)
    } else {
      ""
    }
    buf_clause <- if (!is.null(bnd_info$buffer)) {
      sprintf(
        "; a higher-margin verified alternative was %d participants x %d items (%s scheduled trials, MC lower margin %+.1f percentage points)",
        bnd_info$buffer$n_subj, bnd_info$buffer$n_item,
        format(bnd_info$buffer$n_trials_scheduled, big.mark = ","),
        100 * bnd_info$buffer$mc_margin
      )
    } else {
      ""
    }
    details <- c(details, sprintf(
      "Boundary sensitivity: Because the selected candidate's Monte Carlo lower bound cleared the target threshold by %+.1f percentage points%s, the smallest passing (N_subj, N_item) pair can vary across search seeds%s. Report the candidate range alongside the selected design or confirm it on a fresh high-replication seed via validate_design_point().",
      100 * bnd_info$opt_margin, ind_clause, buf_clause
    ))
  }

  fmt_re_sds <- function(tau_vec) {
    if (length(tau_vec) == 0L) {
      return("0.000 (none)")
    }
    int_sd <- tau_vec[1L]
    slopes <- tau_vec[-1L]
    if (length(slopes) == 0L || all(abs(slopes) < 1e-12)) {
      return(sprintf("%.3f (intercept only)", int_sd))
    }
    if (length(unique(round(slopes, 6))) == 1L) {
      return(sprintf("%.3f (intercept) / %.3f (slopes)", int_sd, slopes[1L]))
    }
    sl_nms <- names(slopes)
    sl_str <- if (!is.null(sl_nms) && all(nzchar(sl_nms))) {
      paste(sprintf("%s = %.3f", sl_nms, slopes), collapse = ", ")
    } else {
      paste(sprintf("%.3f", slopes), collapse = ", ")
    }
    sprintf("%.3f (intercept) / slopes [%s]", int_sd, sl_str)
  }

  fmt_pct_or_na <- function(val) {
    if (is.null(val) || length(val) == 0L || is.na(val)) "NA" else sprintf("%.1f%%", val * 100)
  }

  cov_val <- if (sp$scale == "logit") opt$coverage else opt$prob_coverage
  cov_lbl <- if (sp$scale == "logit") "logit coverage" else sprintf("%s probability coverage", sp$prob_estimand)

  if (!isTRUE(x$feasible)) {
    infeas_att <- switch(
      sp$criterion,
      power = sprintf("Target-scale power was %.1f%% [%.1f%%, %.1f%%]", 100 * opt$power, 100 * opt$power_mc_lo, 100 * opt$power_mc_hi),
      ci_width = sprintf("Width-attainment probability P(CI width <= %.2f %s) was %.1f%% [%.1f%%, %.1f%%]", sp$ci_width, unit_lbl, 100 * opt$aipe_prob, 100 * opt$aipe_mc_lo, 100 * opt$aipe_mc_hi),
      rope = sprintf("ROPE exclusion probability was %.1f%% [%.1f%%, %.1f%%]", 100 * opt$rope_excl_prob, 100 * opt$rope_mc_lo, 100 * opt$rope_mc_hi),
      joint = sprintf("Joint detection-and-width attainment was %.1f%% [%.1f%%, %.1f%%]", 100 * opt$joint_prob, 100 * opt$joint_mc_lo, 100 * opt$joint_mc_hi),
      both = sprintf("Power was %.1f%% [%.1f%%, %.1f%%] and width attainment was %.1f%% [%.1f%%, %.1f%%]", 100 * opt$power, 100 * opt$power_mc_lo, 100 * opt$power_mc_hi, 100 * opt$aipe_prob, 100 * opt$aipe_mc_lo, 100 * opt$aipe_mc_hi)
    )
    if (sp$criterion %in% c("ci_width", "power", "rope")) {
      threshold <- if (sp$criterion == "ci_width") sp$ci_assurance else sp$power_min
      details <- c(details, sprintf("Required attainment was %g%%.", 100 * threshold))
    }
    txt <- sprintf(paste0("### Sample Size Justification\n\nNo evaluated design has confirmed attainment within constraints. ",
      "The candidate is %d participants x %d items (%s, %d replications). ",
      "%s, and mean %s interval width was %.2f percentage points (%.3f logit)."),
      opt$n_subj, opt$n_item, opt$decision_status, opt$nsim, infeas_att, prob_label, opt$prob_ci_width_mean_pp, opt$logit_ci_width_mean)
    if (!is.null(x$robust)) {
      if (isTRUE(x$robust$schedules_differ)) {
        txt <- paste(txt, sprintf(
          "Scheduled trials at this candidate by scenario: %s. The all-scenario sample-size increment is undetermined.",
          format_scenario_counts(x$robust$robust_trials_by_scenario)
        ))
      } else {
        txt <- paste(txt, "The all-scenario sample-size increment is undetermined.")
      }
    }
    txt <- paste(c(txt, details), collapse = " ")
    if (print_output) {
      cat(txt, "\n")
      return(invisible(txt))
    }
    return(txt)
  }

  fac_str <- paste(
    vapply(des$factor_names, function(fn) {
      sprintf("%s (%s)", fn, paste(des$factors[[fn]], collapse = " vs. "))
    }, character(1L)),
    collapse = " x "
  )
  if (one_cell) fac_str <- "intercept only"
  cell_str <- paste(
    sprintf("%s = %.1f%%", names(dgp$target_cell_probs), dgp$target_cell_probs * 100),
    collapse = ", "
  )

  pareto_txt <- ""
  if (!is.null(x$pareto_frontier) && nrow(x$pareto_frontier) > 1L && isTRUE(x$feasible)) {
    pf <- x$pareto_frontier
    pf_pairs <- paste(sprintf("%d participants x %d items", pf$n_subj, pf$n_item), collapse = "; ")
    pareto_txt <- sprintf(" Non-dominated Pareto-efficient alternatives along the target contour included: %s.", pf_pairs)
  }

  sens_txt <- ""
  if (!is.null(x$sensitivity) && nrow(x$sensitivity) >= 2L) {
    s <- x$sensitivity
    w_all <- s$prob_ci_width_mean_pp[is.finite(s$prob_ci_width_mean_pp)]
    lvl_word <- if (nrow(s) == 3L) "three" else as.character(nrow(s))
    if (length(w_all) > 0L) {
      if (no_null_target) {
        sens_txt <- sprintf(
          " Across %s random-slope sensitivity levels, mean CI width ranged from %.2f to %.2f percentage points.",
          lvl_word, min(w_all), max(w_all)
        )
      } else {
        p_all <- s$effective_power[is.finite(s$effective_power)]
        if (length(p_all) > 0L) {
          sens_txt <- sprintf(
            " Across %s random-slope sensitivity levels (0.70x optimistic to 1.35x conservative), effective power ranged from %.1f%% to %.1f%% (mean CI width: %.2f to %.2f percentage points).",
            lvl_word, min(p_all) * 100, max(p_all) * 100,
            min(w_all), max(w_all)
          )
        }
      }
    }
  }

  rob_txt <- ""
  if (!is.null(x$robust)) {
    rc <- x$robust
    rob_trials_clause <- if (isTRUE(rc$schedules_differ)) {
      sprintf("scheduled trials at the robust design by scenario: %s", format_scenario_counts(rc$robust_trials_by_scenario))
    } else {
      sprintf("%+d scheduled trials", rc$delta_trials)
    }
    rob_txt <- sprintf(
      " Compared against the standard-assumption plan (`%s`: N_subj = %d, N_item = %d, %s scheduled trials), ensuring robustness across all evaluated scenarios required %+d participants and %+d items (%s).",
      rc$standard_scenario, rc$optimal_standard$n_subj, rc$optimal_standard$n_item,
      format(rc$optimal_standard$n_trials_scheduled, big.mark = ","),
      rc$delta_subj, rc$delta_item, rob_trials_clause
    )
  }

  chance <- if (dgp$chance > 0) {
    sprintf(" accounting for a lower guessing asymptote of g = %.3f", dgp$chance)
  } else {
    " under a canonical logistic link"
  }

  crit <- switch(
    sp$criterion,
    power = sprintf("statistical power >= %.0f%% (alpha = %.2f)", sp$power_min * 100, sp$alpha),
    ci_width = sprintf("interval width <= %.2f %s with >= %.0f%% assurance (AIPE)", sp$ci_width, unit_lbl, sp$ci_assurance * 100),
    both = sprintf(
      "separate thresholds: statistical power >= %.0f%% (alpha = %.2f) and interval width <= %.2f %s (>= %.0f%% assurance)",
      sp$power_min * 100, sp$alpha, sp$ci_width, unit_lbl, sp$ci_assurance * 100
    ),
    joint = sprintf("detection and width <= %.2f %s in the same replication with probability >= %.0f%%",
                    sp$ci_width, unit_lbl, 100 * sp$power_min),
    rope = sprintf("ROPE [%.2f, %.2f] %s exclusion probability >= %.0f%%", sp$rope[1L], sp$rope[2L], unit_lbl, sp$power_min * 100)
  )

  # MC bounds refer to the denominator selected by the failure policy and primary criterion.
  crit_att <- switch(
    sp$criterion,
    power = sprintf("target-scale power = %.1f%% [%.1f%%, %.1f%%]", opt$power * 100, opt$power_mc_lo * 100, opt$power_mc_hi * 100),
    ci_width = sprintf(
      "AIPE width-attainment probability P(CI width <= %.2f %s) = %.1f%% [%.1f%%, %.1f%%] (power = %.1f%%)",
      sp$ci_width, unit_lbl, opt$aipe_prob * 100, opt$aipe_mc_lo * 100, opt$aipe_mc_hi * 100, opt$power * 100
    ),
    rope = sprintf(
      "ROPE [%.2f, %.2f] %s exclusion probability = %.1f%% [%.1f%%, %.1f%%] (power = %.1f%%)",
      sp$rope[1L], sp$rope[2L], unit_lbl, opt$rope_excl_prob * 100, opt$rope_mc_lo * 100, opt$rope_mc_hi * 100, opt$power * 100
    ),
    joint = sprintf(
      "same-replication detection-and-width probability = %.1f%% [%.1f%%, %.1f%%] (power = %.1f%%, width attainment = %.1f%%)",
      opt$joint_prob * 100, opt$joint_mc_lo * 100, opt$joint_mc_hi * 100, opt$power * 100, opt$aipe_prob * 100
    ),
    both = sprintf(
      "target-scale power = %.1f%% [%.1f%%, %.1f%%] and width-attainment probability = %.1f%% [%.1f%%, %.1f%%]",
      opt$power * 100, opt$power_mc_lo * 100, opt$power_mc_hi * 100,
      opt$aipe_prob * 100, opt$aipe_mc_lo * 100, opt$aipe_mc_hi * 100
    )
  )

  sing_clause <- if (is.na(opt$singular_rate)) {
    "boundary singularity check not applicable under regularized estimation"
  } else {
    sprintf("singular-fit rate = %.1f%%", opt$singular_rate * 100)
  }
  if (no_null_target) crit_att <- sub(" [(]power = [^)]*[)]", "", crit_att)
  opt_mc_lvl <- if ("mc_conf_level_used" %in% names(opt) &&
                    length(opt$mc_conf_level_used) == 1L &&
                    is.finite(opt$mc_conf_level_used)) {
    opt$mc_conf_level_used
  } else if (!is.null(x$search_meta$mc_conf_level_effective) &&
             is.finite(x$search_meta$mc_conf_level_effective)) {
    x$search_meta$mc_conf_level_effective
  } else {
    sp$mc_conf_level
  }
  opt_mc_pct <- 100 * opt_mc_lvl
  opt_mc_str <- if (abs(opt_mc_pct - round(opt_mc_pct)) < 1e-6) {
    sprintf("%.0f%%", opt_mc_pct)
  } else {
    sprintf("%.2f%%", opt_mc_pct)
  }
  mc_clause <- sprintf(
    "Monte Carlo decision status (Morris et al., 2019; %s Wilson CI) was `%s`: %s (valid-fit conditional power = %s across %.1f%% strictly converged replications; %s).",
    opt_mc_str, toupper(opt$decision_status),
    crit_att,
    fmt_pct_or_na(opt$conditional_power),
    opt$valid_rate * 100, sing_clause
  )

  if (no_null_target) {
    mc_clause <- sprintf("Monte Carlo decision status (Morris et al., 2019; %s Wilson CI) was `%s`: %s (%.1f%% strictly converged replications; %s).",
      opt_mc_str, toupper(opt$decision_status), crit_att, opt$valid_rate * 100, sing_clause)
  }

  sched_clause <- if (!is.null(x$robust) && isTRUE(x$robust$schedules_differ)) {
    sprintf("scheduled trials by scenario: %s", format_scenario_counts(x$robust$robust_trials_by_scenario))
  } else {
    sprintf("%s scheduled trials", format(opt$n_trials_scheduled, big.mark = ","))
  }
  exp_clause <- if (!is.null(x$robust) && isTRUE(x$robust$expected_differ) && !is.null(x$robust$robust_expected_by_scenario)) {
    sprintf("expected post-attrition observations by scenario: %s", format_scenario_counts(x$robust$robust_expected_by_scenario))
  } else {
    sprintf("%s expected post-attrition observations", format(opt$n_trials_expected, big.mark = ","))
  }

  txt <- paste0(
    "### Sample Size Justification (glmmplanr)\n\n",
    sprintf(
      "Sample size was planned a priori via two-dimensional simulation-based optimization for crossed binary GLMMs using `glmmplanr` (seed = %d, disjoint L'Ecuyer-CMRG stage/point streams). We specified a balanced `%s` (%s) design with `%s` contrasts (multiples of %d participants x %d items) and %.1f%% random trial attrition.\n\n",
      x$search_meta$seed, des$allocation, fac_str, des$contrast,
      des$k_subj, des$k_item, des$attrition * 100
    ),
    sprintf(
      "Expected marginal cell accuracies (%s) corresponded to a focal `%s` %s %s of %s (calibrated conditional beta = %.3f logit via %d-node Gauss-Hermite quadrature%s). Random-effect SDs were %s for participants and %s for items.\n\n",
      cell_str, sp$effect_input, sp$prob_estimand, prob_label,
      if (is_abs_prob) sprintf("%.2f%%", sp$true_prob_diff_pp) else sprintf("%+.2f percentage points", sp$true_prob_diff_pp),
      sp$true_beta, dgp$n_quad, chance,
      fmt_re_sds(dgp$tau_subj),
      fmt_re_sds(dgp$tau_item)
    ),
    sprintf(
      "Targeting %s under `%s` failure accounting (`%s` engine, `%s` search, `prefer = \"%s\"`, %d Stage-2 replications), the recommended balanced design was N_subj = %d participants and N_item = %d items (%s; %s). %s Mean interval width for the %s %s was %.2f percentage points (%.3f logit; %s = %s).%s%s%s\n",
      crit, sp$failure_policy, sp$engine, x$search_meta$method, x$prefer, opt$nsim,
      opt$n_subj, opt$n_item,
      sched_clause,
      exp_clause,
      mc_clause, sp$prob_estimand, prob_label, opt$prob_ci_width_mean_pp,
      opt$logit_ci_width_mean, cov_lbl, fmt_pct_or_na(cov_val),
      pareto_txt, rob_txt, sens_txt
    )
  )

  if (!is.null(x$robust) && !is.null(opt$scenario)) {
    txt <- paste(txt, sprintf("The reported power and interval width refer to the %s scenario.", opt$scenario))
  }
  txt <- paste(c(txt, details), collapse = " ")
  if (print_output) {
    cat(txt)
    return(invisible(txt))
  }
  txt
}
