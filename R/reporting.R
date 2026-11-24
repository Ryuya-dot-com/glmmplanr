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

assess_inferential_validity <- function(cov, n_valid, alpha) {
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
        "CAUTION (Stage-2 H1 coverage %.1f%% [%.1f%%, %.1f%%] is below nominal %.0f%%; SEs may be anti-conservative; verify with validate_design_point())",
        100 * cov1, 100 * cov_lo, 100 * cov_hi, 100 * nominal
      ),
      cov = cov1, cov_lo = cov_lo, cov_hi = cov_hi, nominal = nominal
    ))
  }
  list(
    status = "unverified",
    short_label = "Unverified",
    detail = sprintf(
      "UNVERIFIED (Stage-2 H1 coverage = %.1f%% [%.1f%%, %.1f%%]; null Type I error not checked; confirm with validate_design_point())",
      100 * cov1, 100 * cov_lo, 100 * cov_hi
    ),
    cov = cov1, cov_lo = cov_lo, cov_hi = cov_hi, nominal = nominal
  )
}

inferential_validity_labels <- function(cov_vec, n_valid_vec, alpha) {
  vapply(seq_along(cov_vec), function(i) {
    assess_inferential_validity(cov_vec[i], n_valid_vec[i], alpha)$short_label
  }, character(1L))
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
    sens_clean <- data.frame(
      section = sprintf(
        "Random-Slope Sensitivity at (N_subj = %d, N_item = %d)",
        x$optimal$n_subj, x$optimal$n_item
      ),
      label = sprintf("%s (%.2fx slope SD)", cap, s$slope_multiplier),
      n_subj = s$n_subj,
      n_item = s$n_item,
      n_trials_scheduled = s$n_trials_scheduled,
      n_trials_expected = s$n_trials_expected,
      decision = toupper(s$decision_status),
      eff_power_ci = fmt_mc_pow(s$power, s$power_mc_lo, s$power_mc_hi),
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

  tbl_df <- switch(
    which,
    verified = ver_clean,
    pareto = pareto_only,
    sensitivity = if (!is.null(sens_clean)) sens_clean else ver_clean,
    robust = if (!is.null(rob_clean)) rob_clean else ver_clean,
    combined = do.call(rbind, Filter(Negate(is.null), list(ver_clean, rob_clean, sens_clean)))
  )

  if (all(is.na(tbl_df$aipe_prob))) {
    tbl_df$aipe_prob <- tbl_df$joint_prob <- NULL
    has_aipe <- FALSE
  } else {
    has_aipe <- TRUE
  }

  size_cols <- c("n_subj", "n_item", "n_trials_scheduled", "n_trials_expected")

  gt_obj <- gt::gt(tbl_df, groupname_col = "section", rowname_col = "label") |>
    gt::tab_header(
      title = "Two-Dimensional GLMM Sample Size Plan & Trustworthiness Summary",
      subtitle = sprintf(
        "Focal effect: %s (%s diff = %+.1f pp | cond. beta = %.3f logit | scale = %s | prefer = %s)",
        x$spec$effect_input, x$spec$prob_estimand, x$spec$true_prob_diff_pp,
        x$spec$true_beta, x$spec$scale, x$prefer
      )
    ) |>
    gt::cols_label(
      n_subj = "N_subj",
      n_item = "N_item",
      n_trials_scheduled = "Sched. Trials",
      n_trials_expected = "Exp. Obs",
      decision = "Target Status",
      eff_power_ci = if (x$spec$failure_policy == "effective") "Effective Power [MC CI]" else "Conditional Power [MC CI]",
      conditional_power = "Cond. Power",
      prob_ci_width_pp = "CI (pp)",
      logit_ci_width = "CI (logit)",
      validity_status = "Inferential Validity",
      valid_rate = "Valid",
      singular_rate = "Singular",
      coverage = "H1 Coverage"
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

  if (has_aipe) {
    gt_obj <- gt_obj |>
      gt::cols_label(aipe_prob = "Width attainment", joint_prob = "Detect and meet width") |>
      gt::fmt_percent(columns = c("aipe_prob", "joint_prob"), decimals = 1) |>
      gt::tab_spanner(
        label = "Power, Precision & MC Decision",
        columns = c("decision", "eff_power_ci", "conditional_power", "prob_ci_width_pp", "logit_ci_width", "aipe_prob", "joint_prob")
      )
  } else {
    gt_obj <- gt_obj |>
      gt::tab_spanner(
        label = "Power, Precision & MC Decision",
        columns = c("decision", "eff_power_ci", "conditional_power", "prob_ci_width_pp", "logit_ci_width")
      )
  }

  gt_obj |>
    gt::tab_source_note(
      source_note = sprintf(
        "Target Status uses %.0f%% Wilson intervals (Morris et al., 2019) for power/width/ROPE targets only: PASS (MC lower >= target), INDETERMINATE (straddles target), FAIL (MC upper < target). Target PASS does not confirm interval coverage or null Type I error; verify inferential validity with validate_design_point(). Seed = %d.",
        x$spec$mc_conf_level * 100,
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

#' Plot 2D Iso-Power / Iso-Precision Contours, 3D Surfaces, and Pareto Trade-offs
#'
#' @param x A `"glmmplanr_result"` object.
#' @param type Visualization type:
#'   * `"contour"` (default): 2D `(N_subj x N_item)` response surface with
#'     iso-power/precision contours, non-dominated Pareto step frontier, and
#'     recommended design.
#'   * `"pareto"`: Non-dominated `(N_subj vs N_item)` Pareto trade-off curve
#'     showing the exact minimum participants required at each item count.
#'   * `"surface3d"`: Interactive `plotly` 3D response surface.
#' @param metric Surface metric to render: `"power"` (default), `"ci_width"`
#'   (on `x$spec$scale`), or `"prob_ci_width_pp"` (95% CI width in percentage
#'   points).
#' @param show_pareto_frontier Logical; overlay the non-dominated `(N_subj, N_item)`
#'   Pareto frontier on `"contour"` plots (default `TRUE`).
#' @param ... Additional arguments (unused).
#' @return A `ggplot` object or `plotly` htmlwidget.
#' @examples
#' # Create the five-replication demonstration plan from the optimizer help.
#' example("optimize_sample_size", package = "glmmplanr", local = environment(), echo = FALSE)
#' plot(plan, type = "pareto")
#' @export
plot.glmmplanr_result <- function(
  x,
  type = c("contour", "surface3d", "pareto"),
  metric = c("power", "ci_width", "prob_ci_width_pp"),
  show_pareto_frontier = TRUE,
  ...
) {
  type <- match.arg(type)
  metric <- match.arg(metric)

  if (type == "contour") {
    return(plot_contour_2d(x, metric = metric, show_pareto_frontier = show_pareto_frontier))
  }
  if (type == "surface3d") {
    return(plot_surface_3d(x, metric = metric))
  }
  plot_pareto_curve(x)
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
    prob_ci_width_pp = "pred_prob_ci_width_pp"
  )
  fill_lbl <- switch(
    metric,
    power = "Predicted\nPower",
    ci_width = sprintf("CI\n(%s)", x$spec$scale),
    prob_ci_width_pp = "CI\n(pp)"
  )

  p <- ggplot2::ggplot(
    surf,
    ggplot2::aes(x = .data$n_subj, y = .data$n_item)
  ) +
    ggplot2::geom_raster(ggplot2::aes(fill = .data[[fill_col]]), interpolate = TRUE, alpha = 0.88)

  if (metric == "power") {
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
      xlim = range(surf$n_subj),
      ylim = range(surf$n_item),
      expand = TRUE,
      clip = "off"
    ) +
    ggplot2::labs(
      title = sprintf(
        "2D Sample Size Trade-off Surface: %s (%s)",
        x$spec$effect_input, toupper(opt$decision_status)
      ),
      subtitle = sprintf(
        "Selected candidate (%s): N_subj = %d, N_item = %d | Power = %.1f%% [%.1f, %.1f] | CI = %.2f pp",
        x$prefer, opt$n_subj, opt$n_item,
        opt$power * 100, opt$power_mc_lo * 100, opt$power_mc_hi * 100,
        opt$prob_ci_width_mean_pp
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
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.caption = ggplot2::element_text(hjust = 0, size = 8.5, color = "#4b5563")
    )

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
    prob_ci_width_pp = "pred_prob_ci_width_pp"
  )
  z_title <- switch(
    metric,
    power = "Power",
    ci_width = sprintf("CI Width (%s)", x$spec$scale),
    prob_ci_width_pp = "CI Width (pp)"
  )

  Z_mat <- matrix(NA_real_, nrow = length(i_vals), ncol = length(s_vals))
  Tr_mat <- matrix(NA_integer_, nrow = length(i_vals), ncol = length(s_vals))
  W_pp_mat <- matrix(NA_real_, nrow = length(i_vals), ncol = length(s_vals))
  P_mat <- matrix(NA_real_, nrow = length(i_vals), ncol = length(s_vals))

  r_ix <- match(surf$n_item, i_vals)
  c_ix <- match(surf$n_subj, s_vals)
  Z_mat[cbind(r_ix, c_ix)] <- surf[[z_col]]
  Tr_mat[cbind(r_ix, c_ix)] <- surf$n_trials_scheduled
  W_pp_mat[cbind(r_ix, c_ix)] <- surf$pred_prob_ci_width_pp
  P_mat[cbind(r_ix, c_ix)] <- surf$pred_power

  key_mat <- outer(i_vals, s_vals, function(ii, ss) paste(ss, ii))
  idx_mat <- match(key_mat, paste(surf$n_subj, surf$n_item))

  hover_txt <- matrix(
    sprintf(
      "N_subj: %d<br>N_item: %d<br>Sched. Trials: %s<br>Power: %.1f%%<br>CI Width: %.2f pp",
      surf$n_subj[idx_mat],
      surf$n_item[idx_mat],
      format(as.numeric(Tr_mat), big.mark = ","),
      as.numeric(P_mat) * 100,
      as.numeric(W_pp_mat)
    ),
    nrow = length(i_vals),
    ncol = length(s_vals)
  )

  z_opt <- switch(
    metric,
    power = opt$power,
    ci_width = opt$ci_width_mean,
    prob_ci_width_pp = opt$prob_ci_width_mean_pp
  )
  thr_val <- if (metric == "power") {
    x$spec$power_min
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
      colorscale = if (metric == "power") "Viridis" else "Magma",
      opacity = 0.90,
      name = z_title
    ) |>
    plotly::add_markers(
      x = opt$n_subj,
      y = opt$n_item,
      z = z_opt,
      marker = list(size = 7, color = "#ef4444", symbol = "diamond"),
      name = sprintf("Selected candidate (%d, %d)", opt$n_subj, opt$n_item)
    )

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

  fig |>
    plotly::layout(
      title = sprintf("Interactive 3D GLMM Response Surface (%s)", x$spec$effect_input),
      scene = list(
        xaxis = list(title = "N_subj (Subjects)"),
        yaxis = list(title = "N_item (Items)"),
        zaxis = list(title = z_title)
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
    ggplot2::theme_minimal(base_size = 11)
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
  prob_label <- if (one_cell) "accuracy" else "accuracy difference"
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
      "Stage-2 Wilson intervals used Bonferroni multiplicity adjustment across %d verified decisions (per-decision level %.2f%%).",
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
  cov_val <- if (sp$scale == "logit") opt$coverage else opt$prob_coverage
  val_info <- assess_inferential_validity(cov_val, opt$n_valid, sp$alpha)
  if (identical(val_info$status, "caution_low_coverage")) {
    details <- c(details, sprintf(
      "Target attainment status reflects only the specified power/precision targets and does not verify inferential validity. Caution: Stage-2 empirical H1 coverage was %.1f%% (95%% Wilson CI %.1f%% to %.1f%%), below the nominal %.0f%% level; narrow intervals or high detection may reflect underestimated standard errors. Confirm interval coverage and null Type I error in independent simulations with validate_design_point() or an external estimator.",
      100 * val_info$cov, 100 * val_info$cov_lo, 100 * val_info$cov_hi, 100 * val_info$nominal
    ))
  } else {
    cov_summary <- if (is.finite(val_info$cov)) {
      sprintf("%.1f%% (95%% Wilson CI %.1f%% to %.1f%%)", 100 * val_info$cov, 100 * val_info$cov_lo, 100 * val_info$cov_hi)
    } else {
      "NA"
    }
    details <- c(details, sprintf(
      "Target attainment status reflects only the specified power/precision targets and does not verify inferential validity (Stage-2 H1 coverage was %s; null Type I error was not evaluated). Confirm interval coverage and Type I error in independent simulations at the selected design with validate_design_point().",
      cov_summary
    ))
  }
  if (!is.null(des$items_per_subject)) {
    details <- c(details, sprintf(
      "Each participant is assigned %2$d distinct items from a bank of %1$d (allocation seed %3$d).",
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
      txt <- paste(txt,
        "The all-scenario sample-size increment is undetermined.")
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
  if (!is.null(x$sensitivity) && nrow(x$sensitivity) >= 3L) {
    s <- x$sensitivity
    p_opt <- s$effective_power[s$scenario == "optimistic"]
    p_con <- s$effective_power[s$scenario == "conservative"]
    w_opt <- s$prob_ci_width_mean_pp[s$scenario == "optimistic"]
    w_con <- s$prob_ci_width_mean_pp[s$scenario == "conservative"]
    sens_txt <- sprintf(
      " Across three random-slope sensitivity levels (0.70x optimistic to 1.35x conservative), effective power ranged from %.1f%% to %.1f%% (mean CI width: %.2f to %.2f percentage points).",
      min(p_con, p_opt) * 100, max(p_con, p_opt) * 100,
      min(w_opt, w_con), max(w_opt, w_con)
    )
  }

  rob_txt <- ""
  if (!is.null(x$robust)) {
    rc <- x$robust
    rob_txt <- sprintf(
      " Compared against the standard-assumption plan (`%s`: N_subj = %d, N_item = %d, %s scheduled trials), ensuring robustness across all evaluated scenarios required %+d participants and %+d items (%+d scheduled trials).",
      rc$standard_scenario, rc$optimal_standard$n_subj, rc$optimal_standard$n_item,
      format(rc$optimal_standard$n_trials_scheduled, big.mark = ","),
      rc$delta_subj, rc$delta_item, rc$delta_trials
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

  mc_clause <- sprintf(
    "Monte Carlo decision status (Morris et al., 2019; %.0f%% Wilson CI) was `%s`: %s (valid-fit conditional power = %s across %.1f%% strictly converged replications; singular-fit rate = %.1f%%).",
    sp$mc_conf_level * 100, toupper(opt$decision_status),
    crit_att,
    fmt_pct_or_na(opt$conditional_power),
    opt$valid_rate * 100, opt$singular_rate * 100
  )

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
      if (one_cell) sprintf("%.2f%%", sp$true_prob_diff_pp) else sprintf("%+.2f percentage points", sp$true_prob_diff_pp),
      sp$true_beta, dgp$n_quad, chance,
      fmt_re_sds(dgp$tau_subj),
      fmt_re_sds(dgp$tau_item)
    ),
    sprintf(
      "Targeting %s under `%s` failure accounting (`%s` engine, `%s` search, `prefer = \"%s\"`, %d Stage-2 replications), the recommended balanced design was N_subj = %d participants and N_item = %d items (%s scheduled trials; %s expected post-attrition observations). %s Mean interval width for the %s %s was %.2f percentage points (%.3f logit; %s = %s).%s%s%s\n",
      crit, sp$failure_policy, sp$engine, x$search_meta$method, x$prefer, opt$nsim,
      opt$n_subj, opt$n_item,
      format(opt$n_trials_scheduled, big.mark = ","),
      format(opt$n_trials_expected, big.mark = ","),
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
