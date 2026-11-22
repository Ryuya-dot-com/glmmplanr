#' Render Sample Size Optimization or Stress-Test Results as a Publication-Ready `gt` Table
#'
#' Formats verified Pareto-frontier sample sizes, scheduled vs. post-attrition
#' trial counts, effective and conditional power with Monte Carlo uncertainty
#' intervals, percentage-point and logit 95% interval widths, convergence /
#' singularity rates, multi-scenario robustness increments, and random-slope
#' sensitivity analyses into a clean `gt` table.
#'
#' @param x A `"glmmplanr_result"` or `"glmmplanr_stress_test"` object.
#' @param which Which table section(s) to display for `"glmmplanr_result"`:
#'   `"combined"` (default), `"verified"`, `"pareto"`, `"sensitivity"`, or
#'   `"robust"`.
#' @param ... Reserved for method extensions.
#' @return A `gt_tbl` object created by [gt::gt()].
#' @export
as_gt <- function(x, ...) {
  UseMethod("as_gt")
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
    valid_rate = ver$valid_rate,
    singular_rate = ver$singular_rate,
    coverage = if (x$spec$scale == "logit") ver$coverage else ver$prob_coverage,
    stringsAsFactors = FALSE
  )

  pareto_only <- ver_clean[is_par, , drop = FALSE]

  sens_clean <- NULL
  if (!is.null(x$sensitivity)) {
    s <- x$sensitivity
    cap <- paste0(toupper(substring(s$scenario, 1L, 1L)), substring(s$scenario, 2L))
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
      valid_rate = s$valid_rate,
      singular_rate = s$singular_rate,
      coverage = if (x$spec$scale == "logit") s$coverage else s$prob_coverage,
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
        valid_rate = pick_rows$valid_rate,
        singular_rate = pick_rows$singular_rate,
        coverage = if (x$spec$scale == "logit") pick_rows$coverage else pick_rows$prob_coverage,
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
      decision = "MC Status",
      eff_power_ci = if (x$spec$failure_policy == "effective") "Effective Power [MC CI]" else "Conditional Power [MC CI]",
      conditional_power = "Cond. Power",
      prob_ci_width_pp = "CI (pp)",
      logit_ci_width = "CI (logit)",
      valid_rate = "Valid",
      singular_rate = "Singular",
      coverage = "Coverage"
    ) |>
    gt::tab_spanner(
      label = "Balanced Sample Size",
      columns = size_cols
    ) |>
    gt::tab_spanner(
      label = "Estimation Health & Calibration",
      columns = c("valid_rate", "singular_rate", "coverage")
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
        "MC Status uses %.0f%% Wilson intervals (Morris et al., 2019): PASS (MC lower >= target), INDETERMINATE (MC interval straddles target), FAIL (MC upper < target). Seed = %d.",
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

    if (x$spec$criterion == "joint" && any(is.finite(surf$pred_joint_prob)) &&
        min(surf$pred_joint_prob) < x$spec$power_min && max(surf$pred_joint_prob) > x$spec$power_min) {
      p <- p + ggplot2::geom_contour(ggplot2::aes(z = .data$pred_joint_prob),
        breaks = x$spec$power_min, color = "#fde047", linetype = "dashed")
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
        "Solid white = target power contour | Dashed yellow = target CI-width assurance\n",
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
      title = "Non-Dominated Sample-Size Pareto Frontier (Subjects vs. Items)",
      subtitle = sprintf(
        "Selected candidate (%s): N_subj = %d, N_item = %d (%s scheduled trials, MC Status = %s)",
        x$prefer, opt$n_subj, opt$n_item,
        format(opt$n_trials_scheduled, big.mark = ","), toupper(opt$decision_status)
      ),
      x = "Minimum Required Subjects (N_subj)",
      y = "Minimum Required Items (N_item)"
    ) +
    ggplot2::theme_minimal(base_size = 11)
}

#' Generate Copy-Paste Methods Justification for Registered Reports & Grant Proposals
#'
#' Produces a reproducible, transparent sample size justification narrative (in
#' English or Japanese) detailing the experimental design, scheduled vs.
#' post-attrition trials, Pareto trade-off frontier, percentage-point and
#' conditional logit effects, estimation failure accounting, Morris et al.
#' (2019) Monte Carlo uncertainty intervals, and multi-scenario / random-slope
#' sensitivity bounds.
#'
#' @param x A `"glmmplanr_result"` object.
#' @param lang Output language: `"en"` (default; English for Registered Reports /
#'   journal manuscripts) or `"ja"` (Japanese for JSPS Kakenhi grant proposals).
#' @param print_output Logical; if `TRUE` (default), prints the formatted text
#'   to the console and returns it invisibly.
#' @return Character string containing the methods report (invisibly when
#'   `print_output = TRUE`).
#' @export
report_plan <- function(x, lang = c("en", "ja"), print_output = TRUE) {
  if (!inherits(x, "glmmplanr_result")) {
    cli::cli_abort("{.arg x} must be a {.cls glmmplanr_result} object.")
  }
  lang <- match.arg(lang)

  opt <- x$optimal
  sp <- x$spec
  if (!is.null(x$robust$specs) && !is.null(opt$scenario)) sp <- x$robust$specs[[opt$scenario]]
  dgp <- sp$dgp
  des <- dgp$design
  unit_lbl <- switch(sp$scale, logit = "logit", pp = "percentage points", prob = "probability units")

  one_cell <- des$n_cells == 1L
  prob_label_en <- if (one_cell) "accuracy" else "accuracy difference"
  prob_label_ja <- if (one_cell) "\u6b63\u7b54\u7387" else "\u6b63\u7b54\u7387\u5dee"
  details <- character()
  if (!is.null(x$search_meta$n_grid_verified)) {
    total <- x$search_meta$n_grid_total
    verified <- x$search_meta$n_grid_verified
    details <- if (lang == "ja") {
      sprintf("\u5168%d\u70b9\u4e2d%d\u70b9\u3092\u6700\u7d42\u691c\u8a3c\u3057\u307e\u3057\u305f\u3002", total, verified)
    } else {
      sprintf("Final verification covered %d of %d lattice points.", verified, total)
    }
  }
  if (sp$criterion %in% c("both", "joint") && !is.null(opt$joint_prob)) {
    rule <- if (sp$criterion == "both") {
      if (lang == "ja") sprintf("\u691c\u51fa\u7387%.0f%%\u3068\u5e45\u9054\u6210\u7387%.0f%%\u3092\u500b\u5225\u306b\u5224\u5b9a\u3002", 100 * sp$power_min, 100 * sp$ci_assurance)
      else sprintf("Separate thresholds are %.0f%% for detection and %.0f%% for width attainment.", 100 * sp$power_min, 100 * sp$ci_assurance)
    } else {
      if (lang == "ja") sprintf("\u540c\u6642\u9054\u6210\u7387\u306e\u76ee\u6a19\u306f%.0f%%\u3002", 100 * sp$power_min)
      else sprintf("The same-replication target is %.0f%%.", 100 * sp$power_min)
    }
    details <- c(details, rule, sprintf(if (lang == "ja")
      "\u5e45\u306e\u4e0a\u9650\u306f%.2f %s\u3001\u5e45\u9054\u6210\u7387\u306f%.1f%%\u3001\u540c\u3058\u53cd\u5fa9\u3067\u306e\u691c\u51fa\u3068\u5e45\u306e\u540c\u6642\u9054\u6210\u7387\u306f%.1f%%\uff08MC\u533a\u9593 %.1f\u2013%.1f%%\uff09\u3002" else
      "Width limit %.2f %s, width attainment %.1f%%, same-replication detection and width attainment %.1f%% (MC interval %.1f to %.1f%%).",
      sp$ci_width, unit_lbl, 100 * opt$aipe_prob, 100 * opt$joint_prob,
      100 * opt$joint_mc_lo, 100 * opt$joint_mc_hi))
  }
  details <- c(details, sprintf(if (lang == "ja") "\u5bfe\u8c61\u52b9\u679c\u306e\u533a\u9593\u6c34\u6e96\u306f%g%%\u3067\u3059\u3002" else
    "The target interval level was %g%%.", 100 * (1 - sp$alpha)))
  details <- c(details, if (lang == "ja")
    "\u5408\u5426\u306f\u6307\u5b9a\u3057\u305f\u8a2d\u8a08\u76ee\u6a19\u306b\u3064\u3044\u3066\u306e\u5224\u5b9a\u3067\u3059\u3002\u533a\u9593\u306e\u88ab\u8986\u7387\u306f\u5019\u88dc\u3054\u3068\u306b\u72ec\u7acb\u3057\u305f\u30b7\u30df\u30e5\u30ec\u30fc\u30b7\u30e7\u30f3\u3067\u8a55\u4fa1\u3067\u304d\u307e\u3059\u3002" else
    "The decision concerns the specified design targets. Interval coverage can be assessed in independent simulations at the selected design.")
  if (!is.null(des$items_per_subject)) {
    details <- c(details, sprintf(if (lang == "ja")
      "\u5404\u53c2\u52a0\u8005\u306b\u5168%d\u9805\u76ee\u304b\u3089%d\u9805\u76ee\u3092\u63d0\u793a\u3057\u307e\u3059\uff08\u5272\u4ed8\u30b7\u30fc\u30c9%d\uff09\u3002" else
      "Each participant is assigned %2$d distinct items from a bank of %1$d (allocation seed %3$d).",
      opt$n_item, des$items_per_subject, des$item_seed))
  }
  if (dropout_rate(des) > 0 && !is.null(opt$n_subj_expected)) {
    details <- c(details, sprintf(if (lang == "ja")
      "\u52df\u96c6%d\u4eba\u3001\u53c2\u52a0\u8005\u8131\u843d\u7387%.1f%%\u3001\u89b3\u6e2c\u304c1\u4ef6\u4ee5\u4e0a\u6b8b\u308b\u671f\u5f85\u4eba\u6570%.1f\u4eba\u3002" else
      "Recruit %d participants, with %.1f%% participant dropout and %.1f expected to supply observations.",
      opt$n_subj, 100 * dropout_rate(des), opt$n_subj_expected))
  }
  if (des$has_testlet) {
    details <- c(details, sprintf(if (lang == "ja")
      "\u6587\u7ae0%d\u672c\u3001\u5404%d\u554f\u3002\u751f\u6210\u30e2\u30c7\u30eb\u306e\u53c2\u52a0\u8005\u00d7\u6587\u7ae0\u306e\u6a19\u6e96\u504f\u5dee\u306f%.2f logit\uff08\u89e3\u6790\u3067\u306f\u6587\u7ae0\u306e\u5171\u6709\u52b9\u679c\u3092\u4f7f\u7528\uff09\u3002" else
      "%d passages with %d items each. The generating participant-by-passage SD is %.2f logit; the analysis uses a shared passage effect.",
      opt$n_item %/% des$testlet_size, des$testlet_size, testlet_interaction_sd(dgp)))
  }

  fmt_re_sds <- function(tau_vec, lang_use = "en") {
    if (length(tau_vec) == 0L) {
      return(if (lang_use == "ja") "0.000\uff08\u5909\u91cf\u52b9\u679c\u306a\u3057\uff09" else "0.000 (none)")
    }
    int_sd <- tau_vec[1L]
    slopes <- tau_vec[-1L]
    if (length(slopes) == 0L || all(abs(slopes) < 1e-12)) {
      return(if (lang_use == "ja") {
        sprintf("\u5207\u7247 %.3f\uff08\u5909\u91cf\u50be\u304d\u306a\u3057\uff09", int_sd)
      } else {
        sprintf("%.3f (intercept only)", int_sd)
      })
    }
    if (length(unique(round(slopes, 6))) == 1L) {
      return(if (lang_use == "ja") {
        sprintf("\u5207\u7247 %.3f / \u50be\u304d %.3f", int_sd, slopes[1L])
      } else {
        sprintf("%.3f (intercept) / %.3f (slopes)", int_sd, slopes[1L])
      })
    }
    sl_nms <- names(slopes)
    sl_str <- if (!is.null(sl_nms) && all(nzchar(sl_nms))) {
      paste(sprintf("%s = %.3f", sl_nms, slopes), collapse = ", ")
    } else {
      paste(sprintf("%.3f", slopes), collapse = ", ")
    }
    if (lang_use == "ja") {
      sprintf("\u5207\u7247 %.3f / \u50be\u304d [%s]", int_sd, sl_str)
    } else {
      sprintf("%.3f (intercept) / slopes [%s]", int_sd, sl_str)
    }
  }

  fmt_pct_or_na <- function(val) {
    if (is.null(val) || length(val) == 0L || is.na(val)) "NA" else sprintf("%.1f%%", val * 100)
  }

  cov_val <- if (sp$scale == "logit") opt$coverage else opt$prob_coverage
  cov_lbl_en <- if (sp$scale == "logit") "logit coverage" else sprintf("%s probability coverage", sp$prob_estimand)
  cov_lbl_ja <- if (sp$scale == "logit") "logit\u88ab\u8986\u7387" else sprintf("%s\u78ba\u7387\u88ab\u8986\u7387", ifelse(sp$prob_estimand == "marginal", "\u5468\u8fba", "\u6761\u4ef6\u4ed8\u304d"))

  if (!isTRUE(x$feasible)) {
    infeas_att_en <- switch(
      sp$criterion,
      power = sprintf("Target-scale power was %.1f%% [%.1f%%, %.1f%%]", 100 * opt$power, 100 * opt$power_mc_lo, 100 * opt$power_mc_hi),
      ci_width = sprintf("Width-attainment probability P(CI width <= %.2f %s) was %.1f%% [%.1f%%, %.1f%%]", sp$ci_width, unit_lbl, 100 * opt$aipe_prob, 100 * opt$aipe_mc_lo, 100 * opt$aipe_mc_hi),
      rope = sprintf("ROPE exclusion probability was %.1f%% [%.1f%%, %.1f%%]", 100 * opt$rope_excl_prob, 100 * opt$rope_mc_lo, 100 * opt$rope_mc_hi),
      joint = sprintf("Joint detection-and-width attainment was %.1f%% [%.1f%%, %.1f%%]", 100 * opt$joint_prob, 100 * opt$joint_mc_lo, 100 * opt$joint_mc_hi),
      both = sprintf("Power was %.1f%% [%.1f%%, %.1f%%] and width attainment was %.1f%% [%.1f%%, %.1f%%]", 100 * opt$power, 100 * opt$power_mc_lo, 100 * opt$power_mc_hi, 100 * opt$aipe_prob, 100 * opt$aipe_mc_lo, 100 * opt$aipe_mc_hi)
    )
    infeas_att_ja <- switch(
      sp$criterion,
      power = sprintf("\u5bfe\u8c61\u5c3a\u5ea6\u306e\u691c\u51fa\u7387\u306f %.1f%%\uff3b%.1f%%, %.1f%%\uff3d", 100 * opt$power, 100 * opt$power_mc_lo, 100 * opt$power_mc_hi),
      ci_width = sprintf("\u533a\u9593\u5e45\u9054\u6210\u7387 P(\u5e45 <= %.2f %s) \u306f %.1f%%\uff3b%.1f%%, %.1f%%\uff3d", sp$ci_width, unit_lbl, 100 * opt$aipe_prob, 100 * opt$aipe_mc_lo, 100 * opt$aipe_mc_hi),
      rope = sprintf("ROPE\u68c4\u5374\u78ba\u7387\u306f %.1f%%\uff3b%.1f%%, %.1f%%\uff3d", 100 * opt$rope_excl_prob, 100 * opt$rope_mc_lo, 100 * opt$rope_mc_hi),
      joint = sprintf("\u691c\u51fa\u30fb\u533a\u9593\u5e45\u306e\u540c\u6642\u9054\u6210\u7387\u306f %.1f%%\uff3b%.1f%%, %.1f%%\uff3d", 100 * opt$joint_prob, 100 * opt$joint_mc_lo, 100 * opt$joint_mc_hi),
      both = sprintf("\u691c\u51fa\u7387 %.1f%%\uff3b%.1f%%, %.1f%%\uff3d\u30fb\u533a\u9593\u5e45\u9054\u6210\u7387 %.1f%%\uff3b%.1f%%, %.1f%%\uff3d", 100 * opt$power, 100 * opt$power_mc_lo, 100 * opt$power_mc_hi, 100 * opt$aipe_prob, 100 * opt$aipe_mc_lo, 100 * opt$aipe_mc_hi)
    )
    txt <- if (lang == "ja") {
      sprintf(paste0("### \u6a19\u672c\u6570\u8a2d\u8a08\u306e\u6839\u62e0\n\n\u76ee\u6a19\u9054\u6210\u3092\u78ba\u8a8d\u3067\u304d\u305f\u8a2d\u8a08\u306f\u3042\u308a\u307e\u305b\u3093\u3002",
        "\u691c\u8a3c\u5019\u88dc\u306f%d\u4eba\u00d7%d\u9805\u76ee\u3001\u5224\u5b9a\u306f %s \u3067\u3059\uff08%d\u53cd\u5fa9\uff09\u3002",
        "%s\u3001%s\u306e\u5e73\u5747\u533a\u9593\u5e45\u306f %.2f \u30d1\u30fc\u30bb\u30f3\u30c8\u30dd\u30a4\u30f3\u30c8\uff08%.3f logit\uff09\u3067\u3059\u3002"),
        opt$n_subj, opt$n_item, opt$decision_status, opt$nsim, infeas_att_ja, prob_label_ja, opt$prob_ci_width_mean_pp, opt$logit_ci_width_mean)
    } else {
      sprintf(paste0("### Sample Size Justification\n\nNo evaluated design has confirmed attainment within constraints. ",
        "The candidate is %d participants x %d items (%s, %d replications). ",
        "%s, and mean %s interval width was %.2f percentage points (%.3f logit)."),
        opt$n_subj, opt$n_item, opt$decision_status, opt$nsim, infeas_att_en, prob_label_en, opt$prob_ci_width_mean_pp, opt$logit_ci_width_mean)
    }
    if (!is.null(x$robust)) {
      txt <- paste(txt, if (lang == "ja") "\u5168\u30b7\u30ca\u30ea\u30aa\u3092\u6e80\u305f\u3059\u8a2d\u8a08\u306e\u7a4d\u307f\u5897\u3057\u91cf\u306f\u672a\u78ba\u5b9a\u3067\u3059\u3002" else
        "The all-scenario sample-size increment is undetermined.")
    }
    txt <- paste(c(txt, details), collapse = " ")
    if (print_output) { cat(txt, "\n"); return(invisible(txt)) }
    return(txt)
  }

  fac_str_en <- paste(
    vapply(des$factor_names, function(fn) {
      sprintf("%s (%s)", fn, paste(des$factors[[fn]], collapse = " vs. "))
    }, character(1L)),
    collapse = " x "
  )
  if (one_cell) fac_str_en <- "intercept only"
  cell_str <- paste(
    sprintf("%s = %.1f%%", names(dgp$target_cell_probs), dgp$target_cell_probs * 100),
    collapse = ", "
  )

  pareto_txt_en <- ""
  pareto_txt_ja <- ""
  if (!is.null(x$pareto_frontier) && nrow(x$pareto_frontier) > 1L && isTRUE(x$feasible)) {
    pf <- x$pareto_frontier
    pf_pairs_en <- paste(sprintf("%d participants x %d items", pf$n_subj, pf$n_item), collapse = "; ")
    pf_pairs_ja <- paste(sprintf("%d\u540d\u00d7%d\u9805\u76ee", pf$n_subj, pf$n_item), collapse = "\u3001")
    pareto_txt_en <- sprintf(" Non-dominated Pareto-efficient alternatives along the target contour included: %s.", pf_pairs_en)
    pareto_txt_ja <- sprintf("\u307e\u305f\u3001\u76ee\u6a19\u7b49\u9ad8\u7dda\u4e0a\u306e\u975e\u52a3\u89e3\uff08\u30d1\u30ec\u30fc\u30c8\u6700\u9069\u89e3\uff09\u3068\u3057\u3066 %s \u306e\u30c8\u30ec\u30fc\u30c9\u30aa\u30d5\u95a2\u4fc2\u304c\u5f97\u3089\u308c\u305f\u3002", pf_pairs_ja)
  }

  sens_txt_en <- ""
  sens_txt_ja <- ""
  if (!is.null(x$sensitivity) && nrow(x$sensitivity) >= 3L) {
    s <- x$sensitivity
    p_opt <- s$effective_power[s$scenario == "optimistic"]
    p_con <- s$effective_power[s$scenario == "conservative"]
    w_opt <- s$prob_ci_width_mean_pp[s$scenario == "optimistic"]
    w_con <- s$prob_ci_width_mean_pp[s$scenario == "conservative"]
    sens_txt_en <- sprintf(
      " Across three random-slope sensitivity levels (0.70x optimistic to 1.35x conservative), effective power ranged from %.1f%% to %.1f%% (mean CI width: %.2f to %.2f percentage points).",
      min(p_con, p_opt) * 100, max(p_con, p_opt) * 100,
      min(w_opt, w_con), max(w_opt, w_con)
    )
    sens_txt_ja <- sprintf(
      "\u3055\u3089\u306b\u5909\u91cf\u50be\u304d\u306e\u6a19\u6e96\u504f\u5dee\u3092\u697d\u89b3\u6761\u4ef6\uff080.70\u500d\uff09\u304b\u3089\u4fdd\u5b88\u6761\u4ef6\uff081.35\u500d\uff09\u306e3\u6c34\u6e96\u306b\u632f\u3063\u305f\u611f\u5ea6\u5206\u6790\u3067\u306f\u3001\u5b9f\u52b9\u691c\u51fa\u529b\u306f%.1f%%\u301c%.1f%%\uff08\u6b63\u7b54\u7387\u5dee\u306e\u5e73\u5747\u533a\u9593\u5e45: %.2f\u301c%.2f\u30dd\u30a4\u30f3\u30c8\uff09\u306e\u7bc4\u56f2\u306b\u4fdd\u305f\u308c\u305f\u3002",
      min(p_con, p_opt) * 100, max(p_con, p_opt) * 100,
      min(w_opt, w_con), max(w_opt, w_con)
    )
  }

  rob_txt_en <- ""
  rob_txt_ja <- ""
  if (!is.null(x$robust)) {
    rc <- x$robust
    rob_txt_en <- sprintf(
      " Compared against the standard-assumption plan (`%s`: N_subj = %d, N_item = %d, %s scheduled trials), ensuring robustness across all evaluated scenarios required %+d participants and %+d items (%+d scheduled trials).",
      rc$standard_scenario, rc$optimal_standard$n_subj, rc$optimal_standard$n_item,
      format(rc$optimal_standard$n_trials_scheduled, big.mark = ","),
      rc$delta_subj, rc$delta_item, rc$delta_trials
    )
    rob_txt_ja <- sprintf(
      "\u306a\u304a\u3001\u6a19\u6e96\u60f3\u5b9a\uff08`%s`\uff1aN_subj = %d\u540d\u30fbN_item = %d\u9805\u76ee\u3001\u4e88\u5b9a %s\u8a66\u884c\uff09\u306e\u307f\u3092\u6e80\u305f\u3059\u8a2d\u8a08\u3068\u6bd4\u8f03\u3059\u308b\u3068\u3001\u5168\u60f3\u5b9a\u30b7\u30ca\u30ea\u30aa\u3092\u540c\u6642\u306b\u6e80\u305f\u3059\u30ed\u30d0\u30b9\u30c8\u8a2d\u8a08\u306b\u306f\u88ab\u9a13\u8005 %+d\u540d\u30fb\u9805\u76ee %+d\u9805\u76ee\uff08\u4e88\u5b9a %+d\u8a66\u884c\uff09\u306e\u7a4d\u307f\u5897\u3057\u304c\u5fc5\u8981\u3067\u3042\u3063\u305f\u3002",
      rc$standard_scenario, rc$optimal_standard$n_subj, rc$optimal_standard$n_item,
      format(rc$optimal_standard$n_trials_scheduled, big.mark = ","),
      rc$delta_subj, rc$delta_item, rc$delta_trials
    )
  }

  chance_en <- if (dgp$chance > 0) {
    sprintf(" accounting for a lower guessing asymptote of g = %.3f", dgp$chance)
  } else {
    " under a canonical logistic link"
  }

  crit_en <- switch(
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
  crit_att_en <- switch(
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

  crit_att_ja <- switch(
    sp$criterion,
    power = sprintf("\u6307\u5b9a\u57fa\u6e96\u306e\u691c\u51fa\u529b\u306f %.1f%%\uff3b%.1f%%, %.1f%%\uff3d", opt$power * 100, opt$power_mc_lo * 100, opt$power_mc_hi * 100),
    ci_width = sprintf(
      "\u76ee\u6a19\u533a\u9593\u5e45\u9054\u6210\u78ba\u7387 P(\u533a\u9593\u5e45 <= %.2f %s) \u306f %.1f%%\uff3b%.1f%%, %.1f%%\uff3d\uff08\u691c\u51fa\u529b %.1f%%\uff09",
      sp$ci_width, unit_lbl, opt$aipe_prob * 100, opt$aipe_mc_lo * 100, opt$aipe_mc_hi * 100, opt$power * 100
    ),
    rope = sprintf(
      "ROPE [%.2f, %.2f] %s \u68c4\u5374\u78ba\u7387\u306f %.1f%%\uff3b%.1f%%, %.1f%%\uff3d\uff08\u691c\u51fa\u529b %.1f%%\uff09",
      sp$rope[1L], sp$rope[2L], unit_lbl, opt$rope_excl_prob * 100, opt$rope_mc_lo * 100, opt$rope_mc_hi * 100, opt$power * 100
    ),
    joint = sprintf(
      "\u691c\u51fa\u30fb\u533a\u9593\u5e45\u306e\u540c\u6642\u9054\u6210\u7387\u306f %.1f%%\uff3b%.1f%%, %.1f%%\uff3d\uff08\u691c\u51fa\u529b %.1f%%\u3001\u5e45\u9054\u6210\u7387 %.1f%%\uff09",
      opt$joint_prob * 100, opt$joint_mc_lo * 100, opt$joint_mc_hi * 100, opt$power * 100, opt$aipe_prob * 100
    ),
    both = sprintf(
      "\u691c\u51fa\u529b %.1f%%\uff3b%.1f%%, %.1f%%\uff3d\u304a\u3088\u3073\u533a\u9593\u5e45\u9054\u6210\u78ba\u7387 %.1f%%\uff3b%.1f%%, %.1f%%\uff3d",
      opt$power * 100, opt$power_mc_lo * 100, opt$power_mc_hi * 100,
      opt$aipe_prob * 100, opt$aipe_mc_lo * 100, opt$aipe_mc_hi * 100
    )
  )

  if (lang == "en") {
    mc_clause <- sprintf(
      "Monte Carlo decision status (Morris et al., 2019; %.0f%% Wilson CI) was `%s`: %s (valid-fit conditional power = %s across %.1f%% strictly converged replications; singular-fit rate = %.1f%%).",
      sp$mc_conf_level * 100, toupper(opt$decision_status),
      crit_att_en,
      fmt_pct_or_na(opt$conditional_power),
      opt$valid_rate * 100, opt$singular_rate * 100
    )

    txt <- paste0(
      "### Sample Size Justification (glmmplanr)\n\n",
      sprintf(
        "Sample size was planned a priori via two-dimensional simulation-based optimization for crossed binary GLMMs using `glmmplanr` (seed = %d, disjoint L'Ecuyer-CMRG stage/point streams). We specified a balanced `%s` (%s) design with `%s` contrasts (multiples of %d participants x %d items) and %.1f%% random trial attrition.\n\n",
        x$search_meta$seed, des$allocation, fac_str_en, des$contrast,
        des$k_subj, des$k_item, des$attrition * 100
      ),
      sprintf(
        "Expected marginal cell accuracies (%s) corresponded to a focal `%s` %s %s of %s (calibrated conditional beta = %.3f logit via %dx%d Gauss-Hermite quadrature%s). Random-effect SDs were %s for participants and %s for items.\n\n",
        cell_str, sp$effect_input, sp$prob_estimand, prob_label_en,
        if (one_cell) sprintf("%.2f%%", sp$true_prob_diff_pp) else sprintf("%+.2f percentage points", sp$true_prob_diff_pp),
        sp$true_beta, dgp$n_quad, dgp$n_quad, chance_en,
        fmt_re_sds(dgp$tau_subj, "en"),
        fmt_re_sds(dgp$tau_item, "en")
      ),
      sprintf(
        "Targeting %s under `%s` failure accounting (`%s` engine, `%s` search, `prefer = \"%s\"`, %d Stage-2 replications), the recommended balanced design was N_subj = %d participants and N_item = %d items (%s scheduled trials; %s expected post-attrition observations). %s Mean interval width for the %s %s was %.2f percentage points (%.3f logit; %s = %s).%s%s%s\n",
        crit_en, sp$failure_policy, sp$engine, x$search_meta$method, x$prefer, opt$nsim,
        opt$n_subj, opt$n_item,
        format(opt$n_trials_scheduled, big.mark = ","),
        format(opt$n_trials_expected, big.mark = ","),
        mc_clause, sp$prob_estimand, prob_label_en, opt$prob_ci_width_mean_pp,
        opt$logit_ci_width_mean, cov_lbl_en, fmt_pct_or_na(cov_val),
        pareto_txt_en, rob_txt_en, sens_txt_en
      )
    )
  } else {
    dec_ja <- switch(
      opt$decision_status,
      pass = "\u5408\u683c\uff08PASS\uff1aMC\u4fe1\u983c\u533a\u9593\u4e0b\u9650\u304c\u76ee\u6a19\u4ee5\u4e0a\uff09",
      indeterminate = "\u5224\u5b9a\u4fdd\u7559\uff08INDETERMINATE\uff1aMC\u4fe1\u983c\u533a\u9593\u304c\u76ee\u6a19\u95be\u5024\u3092\u307e\u305f\u3050\u5883\u754c\u9818\u57df\uff09",
      fail = "\u672a\u9054\uff08FAIL\uff1aMC\u4fe1\u983c\u533a\u9593\u4e0a\u9650\u304c\u76ee\u6a19\u672a\u6e80\uff09"
    )
    txt <- paste0(
      "### \u6a19\u672c\u6570\u8a2d\u8a08\u306e\u6839\u62e0\uff08glmmplanr \u306b\u3088\u308b\u4e8b\u524d\u691c\u51fa\u529b\u30fb\u7cbe\u5ea6\u30fb\u4fe1\u983c\u6027\u8a55\u4fa1\uff09\n\n",
      sprintf(
        "\u672c\u7814\u7a76\u306e\u6a19\u672c\u6570\uff08\u88ab\u9a13\u8005\u6570\u00d7\u9805\u76ee\u6570\uff09\u306f\u3001R\u30d1\u30c3\u30b1\u30fc\u30b8 `glmmplanr`\uff08\u4e71\u6570\u30b7\u30fc\u30c9 = %d\uff1b\u6bb5\u968e\u30fb\u8a2d\u8a08\u70b9\u3054\u3068\u306b\u72ec\u7acb\u306a L'Ecuyer-CMRG \u30b9\u30c8\u30ea\u30fc\u30e0\u5272\u5f53\uff09\u3092\u7528\u3044\u3066\u4e8b\u524d\u6c7a\u5b9a\u3057\u305f\u3002\u5b9f\u9a13\u8a08\u753b\u306f`%s` \u30b3\u30f3\u30c8\u30e9\u30b9\u30c8\u4ed8\u304d\u5747\u8861 `%s` \u30c7\u30b6\u30a4\u30f3\uff08\u88ab\u9a13\u8005%d\u306e\u500d\u6570\u00d7\u9805\u76ee%d\u306e\u500d\u6570\u3001\u8a66\u884c\u6b20\u6e2c\u7387%.1f%%\uff09\u3068\u3057\u305f\u3002\n\n",
        x$search_meta$seed, des$contrast, des$allocation,
        des$k_subj, des$k_item, des$attrition * 100
      ),
      sprintf(
        "\u5404\u6761\u4ef6\u306e\u671f\u5f85\u5468\u8fba\u6b63\u7b54\u7387\uff08%s\uff09\u306f\u3001\u5bfe\u8c61\u52b9\u679c `%s` \u306b\u304a\u3044\u3066 %s\u306e%s%s\u306b\u76f8\u5f53\u3057\u3001%dx%d \u30ac\u30a6\u30b9\u30fb\u30a8\u30eb\u30df\u30fc\u30c8\u6c42\u7a4d\u306b\u3088\u308a\u9006\u7b97\u3057\u305f\u6761\u4ef6\u4ed8\u304d\u6bcd\u6570\u306f beta = %.3f logit\uff08\u30c1\u30e3\u30f3\u30b9\u30ec\u30d9\u30eb g = %.3f\uff09\u3001\u5909\u91cf\u52b9\u679cSD\u306f\u88ab\u9a13\u8005\uff08%s\uff09\u30fb\u9805\u76ee\uff08%s\uff09\u3067\u3042\u3063\u305f\u3002\n\n",
        cell_str, sp$effect_input,
        if (one_cell) sprintf("%.2f%%", sp$true_prob_diff_pp) else sprintf("%+.2f \u30d1\u30fc\u30bb\u30f3\u30c8\u30dd\u30a4\u30f3\u30c8", sp$true_prob_diff_pp),
        ifelse(sp$prob_estimand == "marginal", "\u96c6\u56e3\u5468\u8fba", "\u6761\u4ef6\u4ed8\u304d"),
        prob_label_ja, dgp$n_quad, dgp$n_quad, sp$true_beta, dgp$chance,
        fmt_re_sds(dgp$tau_subj, "ja"),
        fmt_re_sds(dgp$tau_item, "ja")
      ),
      sprintf(
        "\u63a8\u5b9a\u5931\u6557\u306e\u6271\u3044\u306f `%s` \u3068\u3057\u3001 `%s` \u63a2\u7d22\u304a\u3088\u3073\u6700\u7d42\u691c\u8a3c\uff08\u53cd\u5fa9\u6570 %d\u56de\u3001\u30a8\u30f3\u30b8\u30f3 `%s`\u3001\u4ee3\u8868\u89e3\u9078\u629e\u898f\u5247 `prefer = \"%s\"`\uff09\u3092\u5b9f\u65bd\u3057\u305f\u7d50\u679c\u3001\u63a8\u5968\u30c7\u30b6\u30a4\u30f3\u306f\u88ab\u9a13\u8005 N_subj = %d\u540d\u30fb\u9805\u76ee N_item = %d\u9805\u76ee\uff08\u5272\u4ed8\u4e88\u5b9a\u8a66\u884c\u6570 %s\u8a66\u884c\u3001\u6b20\u6e2c\u5f8c\u671f\u5f85\u89b3\u6e2c\u6570 %s\u8a66\u884c\uff09\u3068\u306a\u3063\u305f\u3002\u30b7\u30df\u30e5\u30ec\u30fc\u30b7\u30e7\u30f3\u8aa4\u5dee\u3092\u8003\u616e\u3057\u305f\u5224\u5b9a\uff08Morris et al., 2019\uff1b%.0f%% Wilson\u533a\u9593\uff09\u306f\u300c%s\u300d\u3067\u3042\u308a\u3001%s\uff08\u6709\u52b9\u53ce\u675f\u7387 %.1f%% \u4e0b\u306e\u6761\u4ef6\u4ed8\u304d\u691c\u51fa\u7387 %s\u3001\u7279\u7570\u9069\u5408\u7387 %.1f%%\uff09\u3001%s\u306e\u5e73\u5747\u4fe1\u983c\u533a\u9593\u5168\u5e45\u306f %.2f \u30dd\u30a4\u30f3\u30c8\uff08%.3f logit\u3001%s %s\uff09\u3067\u3042\u3063\u305f\u3002%s%s%s\n",
        sp$failure_policy, x$search_meta$method, opt$nsim, sp$engine, x$prefer,
        opt$n_subj, opt$n_item,
        format(opt$n_trials_scheduled, big.mark = ","),
        format(opt$n_trials_expected, big.mark = ","),
        sp$mc_conf_level * 100, dec_ja,
        crit_att_ja,
        opt$valid_rate * 100,
        fmt_pct_or_na(opt$conditional_power),
        opt$singular_rate * 100, prob_label_ja,
        opt$prob_ci_width_mean_pp, opt$logit_ci_width_mean,
        cov_lbl_ja, fmt_pct_or_na(cov_val),
        pareto_txt_ja, rob_txt_ja, sens_txt_ja
      )
    )
  }

  if (!is.null(x$robust) && !is.null(opt$scenario)) {
    txt <- paste(txt, if (lang == "ja") sprintf("\u8a18\u8f09\u3057\u305f\u691c\u51fa\u7387\u3068\u533a\u9593\u5e45\u306f %s \u30b7\u30ca\u30ea\u30aa\u306e\u5024\u3067\u3059\u3002", opt$scenario)
                 else sprintf("The reported power and interval width refer to the %s scenario.", opt$scenario))
  }
  txt <- paste(c(txt, details), collapse = " ")
  if (isTRUE(print_output)) {
    cat(txt)
    return(invisible(txt))
  }
  txt
}
