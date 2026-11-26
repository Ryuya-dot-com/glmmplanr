#' Extract verified GLMM candidates as a data frame
#'
#' @param x A result from [optimize_sample_size()].
#' @param row.names,optional,... Passed to the data frame method.
#' @return One row per verified candidate, including failed and indeterminate
#'   candidates. Numeric values retain full precision. Nested replication
#'   records are available in the original object.
#' @export
as.data.frame.glmmplanr_result <- function(x, row.names = NULL, optional = FALSE, ...) {
  df <- x$verified
  as.data.frame(df[!vapply(df, is.list, logical(1))], row.names = row.names,
                optional = optional, ...)
}

#' Inspect and plot Rasch planning results
#'
#' @name rasch_plan_methods
#' @param x,object A result from [plan_rasch_fixed_bank()] or [plan_rasch_precision()].
#' @param row.names,optional Arguments to the data frame method.
#' @param type The default shows attainment for a fixed bank and decisions at
#'   the tested participant/item pairs for random items. Use `"items"` for
#'   item-level errors and finite-estimate counts in a fixed bank.
#' @param target Name of a recovery target, or `"joint"` for all tolerances in
#'   the same replication. Applies to the default plot.
#' @param metric Error measure for the item plot, in logits.
#' @inheritParams plot.glmmplanr_result
#' @param ... Additional arguments passed to data frame coercion, otherwise unused.
#' @return `summary()` and `as.data.frame()` return the full unrounded table,
#'   one row per candidate and target. `print()` returns the plan invisibly.
#'   `plot()` returns an editable ggplot. Item errors use finite estimates
#'   only. The count panel retains every item, including those with no estimate.
#'   Figures use [ggplot2::theme_bw()]. Add [ggplot2::labs()] to supply your
#'   own title or caption.
#' @examples
#' bank <- plan_rasch_fixed_bank(c(-1, 0, 1), 20L, c(difficulty = .5), nsim = 4L)
#' subset(as.data.frame(bank), target == "joint")
#' plot(bank, type = "items", metric = "bias")
NULL

#' @rdname rasch_plan_methods
#' @export
as.data.frame.glmmplanr_rasch <- function(x, row.names = NULL, optional = FALSE, ...) {
  as.data.frame(x$summary, row.names = row.names, optional = optional, ...)
}

#' @rdname rasch_plan_methods
#' @export
summary.glmmplanr_rasch <- function(object, ...) {
  as.data.frame(object)
}

rasch_decision_text <- function(x) {
  joint <- x$summary[x$summary$target == "joint", , drop = FALSE]
  if (!any(joint$decision == "pass")) {
    return(sprintf("No tested design has a confirmed pass (%d indeterminate, %d failed).",
      sum(joint$decision == "indeterminate"), sum(joint$decision == "fail")))
  }
  if (!is.null(x$item_summary)) {
    return(sprintf("Smallest tested sample with a confirmed pass: %d examinees.",
                   x$recommended_n_subj))
  }
  pairs <- sprintf("%d x %d", x$pareto$n_subj, x$pareto$n_item)
  paste("Non-dominated passing designs (examinees x items):", paste(pairs, collapse = ", "))
}

#' @rdname rasch_plan_methods
#' @export
print.glmmplanr_rasch <- function(x, ...) {
  joint <- x$summary[x$summary$target == "joint", , drop = FALSE]
  cat(if (is.null(x$item_summary)) "Random-item" else "Fixed-bank",
      "Rasch plan | Joint assurance target", sprintf("%g%%\n", 100 * x$assurance))
  cat(rasch_decision_text(x), "\n")
  cols <- c("n_subj", "n_item", "n_valid_fit", "n_evaluable", "n_met", "nsim",
            "mc_lo", "mc_hi", "decision")
  print(utils::head(joint[cols], 10L), digits = 3, row.names = FALSE)
  if (nrow(joint) > 10L) cat("First 10 candidates shown. See summary() for all candidates.\n")
  invisible(x)
}

rasch_report <- function(x) {
  joint <- x$summary[x$summary$target == "joint", , drop = FALSE]
  fixed <- !is.null(x$item_summary)
  setup <- if (fixed) {
    sprintf("The same %d item difficulties were used in every replication with newly sampled examinees (%s, %s).",
            length(x$item_difficulties), x$method, x$engine)
  } else {
    sprintf("Examinees and items were sampled anew in every replication (%s).", x$engine)
  }
  # Report the smallest passing bank design or the non-dominated random-item pairs.
  # Without a pass, identify the tested candidate with the largest lower bound.
  selected <- if (!any(joint$decision == "pass")) {
    joint[which.max(joint$mc_lo), , drop = FALSE]
  } else if (fixed) {
    joint[joint$n_subj == x$recommended_n_subj, , drop = FALSE]
  } else x$pareto
  results <- sprintf(paste0(
    "At %d examinees and %d items, %d/%d studies met all tolerances ",
    "(%g%% Wilson Monte Carlo interval %.1f%% to %.1f%%, %s). ",
    "Fit validity was %d/%d and joint evaluability was %d/%d."),
    selected$n_subj, selected$n_item, selected$n_met, selected$nsim, 100 * x$mc_conf_level,
    100 * selected$mc_lo, 100 * selected$mc_hi, selected$decision,
    selected$n_valid_fit, selected$nsim, selected$n_evaluable, selected$nsim)
  paste(setup, sprintf("The error limits were %s, with joint assurance of %g%%.",
    paste(sprintf("%s <= %g", names(x$max_error), x$max_error), collapse = ", "), 100 * x$assurance),
    "All planned studies contributed to the denominator.", rasch_decision_text(x),
    paste(results, collapse = " "))
}

#' @rdname rasch_plan_methods
#' @export
plot.glmmplanr_rasch <- function(x, type = c("design", "items"), target = "joint",
                                 metric = c("rmse", "bias"), ..., show_annotations = FALSE) {
  type <- match.arg(type)
  if (type == "items") {
    if (is.null(x$item_summary)) cli::cli_abort("Item-level errors require a fixed bank.")
    metric <- match.arg(metric)
    items <- x$item_summary
    items$sample <- factor(items$n_subj, levels = unique(items$n_subj),
                           labels = paste(unique(items$n_subj), "examinees"))
    errors <- items
    errors$value <- items[[metric]]
    errors$measure <- if (metric == "rmse") "RMSE (logits)" else "Bias (logits)"
    items$value <- items$n_finite
    items$measure <- "Finite estimates (studies)"
    df <- rbind(errors, items)
    p <- ggplot2::ggplot(df, ggplot2::aes(.data$true_difficulty, .data$value)) +
      ggplot2::geom_point(na.rm = TRUE) +
      ggplot2::facet_grid(measure ~ sample, scales = "free_y") +
      ggplot2::labs(x = "True item difficulty (logits)", y = NULL,
        title = "Item recovery across the fixed bank",
        subtitle = sprintf("Errors use finite estimates. Counts include every item (%d studies per design).",
                           x$summary$nsim[1L])) + ggplot2::theme_bw()
    return(style_plot_annotations(p, show_annotations))
  }
  target <- match.arg(target, unique(x$summary$target))
  df <- x$summary[x$summary$target == target, , drop = FALSE]
  if (!is.null(x$item_summary)) {
    p <- ggplot2::ggplot(df, ggplot2::aes(.data$n_subj, .data$attainment)) +
      ggplot2::geom_hline(yintercept = x$assurance, linetype = "dashed") +
      ggplot2::geom_linerange(ggplot2::aes(ymin = .data$mc_lo, ymax = .data$mc_hi)) +
      ggplot2::scale_y_continuous(limits = c(0, 1), labels = function(v) paste0(v * 100, "%")) +
      ggplot2::labs(y = "Studies meeting the target / all planned studies",
        subtitle = sprintf("%g%% Wilson Monte Carlo intervals. Dashed line at %g%% required assurance.",
                           100 * x$mc_conf_level, 100 * x$assurance))
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(.data$n_subj, .data$n_item)) +
      ggplot2::labs(y = "Items")
    if (target == "joint" && nrow(x$pareto)) {
      p <- p + ggplot2::geom_point(data = x$pareto, shape = 1, size = 6) +
        ggplot2::labs(subtitle = "Rings mark non-dominated passing designs.")
    }
  }
  p <- p + ggplot2::geom_point(ggplot2::aes(shape = .data$decision), size = 3) +
    ggplot2::scale_shape_manual(values = c(pass = 16, indeterminate = 17, fail = 4)) +
    ggplot2::labs(title = paste("Rasch recovery for", target), x = "Examinees", shape = "Decision") +
    ggplot2::theme_bw()
  style_plot_annotations(p, show_annotations)
}

#' @rdname as_gt
#' @export
as_gt.glmmplanr_rasch <- function(x, ...) {
  gt::gt(as.data.frame(x)) |>
    gt::fmt_percent(columns = c("attainment", "mc_lo", "mc_hi"), decimals = 1)
}

#' Use the planned model structure in an analysis of observed responses
#'
#' @param x A specification from [specify_target()].
#' @param object A specification for `family()`, or a calibrated DGP for `coef()`.
#' @param ... Reserved for method extensions.
#' @return `formula()` gives the primary mixed-model formula for binary response
#'   `y` grouped by `subj`, `item` and optionally `testlet`. Give the observed
#'   factors the levels and contrast matrices from the design. For independent
#'   slopes, the formula uses `cell_id` (the row index in `design$cells`) to
#'   look up each numeric random-slope column separately. `family()` supplies
#'   the planned analysis link and guessing floor. These describe the model
#'   structure; the planning engine, priors and interval method also affect
#'   inference. `coef()` gives the named conditional DGP coefficients in logits.
#' @examples
#' d <- design_crossed(list(learning = c("incidental", "instruction")))
#' g <- calibrate_dgp(d, c(.60, .70), chance = 1/3, subj_sd = .8, item_sd = .6)
#' s <- specify_target(g, "learning", analysis_re = "intercept_only")
#' formula(s)
#' family(s)
#' coef(g)
#' @name analysis_model
NULL

#' @rdname analysis_model
#' @importFrom stats formula
#' @export
formula.glmmplanr_spec <- function(x, ...) {
  d <- x$dgp$design
  # Keep only the required lookup matrices in the returned formula environment.
  env <- new.env(parent = asNamespace("stats"))
  groups <- c("subj", "item", if (d$has_testlet) "testlet")
  block <- function(group) {
    rhs <- d[[paste0(group, "_rhs")]][[2L]]
    if (x$analysis_re == "intercept_only" || identical(rhs, 1)) {
      return(sprintf("(1 | %s)", group))
    }
    if (x$analysis_re == "maximal") return(sprintf("(%s | %s)", deparse1(rhs), group))
    nm <- paste0("Z_", group, "_cells")
    z <- d[[nm]]
    env[[nm]] <- z
    c(sprintf("(1 | %s)", group),
      sprintf("(0 + I(%s[cell_id, %dL]) | %s)", nm, seq.int(2L, ncol(z)), group))
  }
  terms <- c(deparse1(d$fixed_rhs[[2L]]), unlist(lapply(groups, block)))
  stats::as.formula(paste("y ~", paste(terms, collapse = " + ")), env = env)
}

#' @rdname analysis_model
#' @importFrom stats family
#' @export
family.glmmplanr_spec <- function(object, ...) {
  mafc_binomial(if (object$analysis_link == "match_dgp") object$dgp$chance else 0)
}

#' @rdname analysis_model
#' @importFrom stats coef
#' @export
coef.glmmplanr_dgp <- function(object, ...) {
  object$beta
}

analysis_description <- function(sp) {
  stan <- sp$engine %in% c("stan_laplace", "stan_hmc")
  engine_label <- if (identical(sp$engine, "stan_laplace")) {
    "stan_hmc (via deprecated 'stan_laplace' alias)"
  } else if (!is.null(sp$requested_engine) && !identical(sp$requested_engine, sp$engine)) {
    sprintf("%s (fallback from requested '%s')", sp$engine, sp$requested_engine)
  } else {
    sp$engine
  }
  method <- if (stan) {
    "posterior draws"
  } else if (sp$scale == "logit") {
    "coefficient standard error"
  } else if (sp$prob_estimand == "conditional") {
    "conditional delta method"
  } else if (identical(sp$prob_se_method, "joint_delta")) {
    "joint delta method"
  } else {
    sprintf("parametric bootstrap (%d refits per successful replication)", if (is.null(sp$prob_boot)) 99L else sp$prob_boot)
  }
  sprintf("Primary analysis: %s, %s random effects, %s link. Target uncertainty: %s.",
    engine_label, sp$analysis_re, stats::family(sp)$link, method)
}
