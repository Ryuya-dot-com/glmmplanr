library(ggplot2)
args <- commandArgs(TRUE)
stage <- if (length(args)) match.arg(args[1L], c("all", "targeting", "replication")) else "all"
if (stage %in% c("all", "targeting")) {
  x <- read.csv("validation/l2_targeting_large/summary.csv")
  x <- subset(x, model == "Rasch" & engine == "native" & layout == "full50" &
    profile %in% c("reference", "lower", "lower_narrow", "severe_lower"))
  x$profile <- factor(x$profile, c("reference", "lower", "lower_narrow", "severe_lower"),
    c("Reference", "Mean -3", "Mean -3; half SD", "Mean -5"))
  p <- ggplot(x, aes(n_subj, attainment, color = profile, shape = profile)) +
    geom_hline(yintercept = .8, color = "grey55", linewidth = .4) +
    geom_errorbar(aes(ymin = mc_lo, ymax = mc_hi), width = .02, linewidth = .35) +
    geom_line(linewidth = .7) + geom_point(size = 2) + facet_wrap(~ format) +
    scale_color_manual(values = c("#0072B2", "#D55E00", "#CC79A7", "#333333")) +
    scale_x_continuous(trans = "log10", breaks = c(100, 200, 500, 1000, 2000)) +
    scale_y_continuous(limits = c(0, 1), labels = function(y) paste0(100 * y, "%")) +
    labs(title = "Calibration gains depend on population/item targeting",
      subtitle = "50-item Rasch banks refitted to public vocabulary responses; 1,024 fresh replications per condition",
      x = "Participants (log scale)", y = "Whole-bank difficulty attainment", color = NULL, shape = NULL,
      caption = "Pointwise Wilson 95% Monte Carlo intervals. Failures remain in the denominator.\nTarget: RMSE <= 0.35 and maximum item error <= 0.60 logits. Each format has its own reference scale.\nSource parameters are fixed; ability shifts are sensitivity assumptions, not age or CEFR labels.") +
    theme_minimal(base_size = 12) + theme(legend.position = "bottom", panel.grid.minor = element_blank(),
      plot.caption = element_text(hjust = 0))
  ggsave("validation/l2_targeting_large.png", p, width = 10.5, height = 5.6, dpi = 160, bg = "white")
}
if (stage %in% c("all", "replication")) {
  x <- read.csv("validation/l2_replication_large/design_decisions.csv")
  labels <- c(published_8pp = "8pp; intercepts only", slopes_8pp = "8pp; moderate slopes",
    conservative_5pp = "5pp; moderate slopes", high_variance_8pp = "8pp; larger variances")
  x$profile <- factor(x$profile, names(labels), labels)
  x$pool <- factor(x$n_item)
  p <- ggplot(x, aes(n_subj, joint_attainment, color = pool, shape = pool)) +
    geom_hline(yintercept = .8, color = "grey55", linewidth = .4) +
    geom_errorbar(aes(ymin = joint_lo, ymax = joint_hi), width = 3, linewidth = .35) +
    geom_line(linewidth = .7) + geom_point(size = 2) + facet_wrap(~ profile, ncol = 2) +
    scale_color_manual(values = c("#0072B2", "#D55E00", "#009E73")) +
    scale_x_continuous(breaks = c(54, 72, 108, 162, 216)) +
    scale_y_continuous(limits = c(0, 1), labels = function(y) paste0(100 * y, "%")) +
    labs(title = "Replication planning depends on effect size and random slopes",
      subtitle = "30 responses per person; detection AND 95% interval width <= 10 percentage points",
      x = "Participants", y = "Joint attainment", color = "Item pool", shape = "Item pool",
      caption = "2,048 independent planned replications per design; pointwise Wilson 95% Monte Carlo intervals.\nThe 80% planning target does not establish nominal coverage or Type I error control.\nPublished means inform this focused study; unpublished random-effect variances are assumptions.") +
    theme_minimal(base_size = 12) + theme(legend.position = "bottom", panel.grid.minor = element_blank(),
      plot.caption = element_text(hjust = 0))
  ggsave("validation/l2_replication_large.png", p, width = 10.5, height = 7.5, dpi = 160, bg = "white")
  diagnosis <- rbind(
    data.frame(x, metric = "Null rejection", value = x$null_detection, lo = x$null_detection_lo,
      hi = x$null_detection_hi, nominal = .05),
    data.frame(x, metric = "Planned coverage", value = x$coverage, lo = x$coverage_lo,
      hi = x$coverage_hi, nominal = .95))
  d <- ggplot(diagnosis, aes(n_subj, value, color = pool, shape = pool)) +
    geom_hline(aes(yintercept = nominal), color = "grey55", linewidth = .4) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 3, linewidth = .35) +
    geom_line(linewidth = .7) + geom_point(size = 1.6) +
    facet_grid(metric ~ profile, scales = "free_y") +
    scale_color_manual(values = c("#0072B2", "#D55E00", "#009E73")) +
    scale_x_continuous(breaks = c(54, 108, 162, 216)) +
    scale_y_continuous(labels = function(y) paste0(round(100 * y, 1), "%")) +
    labs(title = "Check nominal inference separately from the planning target",
      subtitle = "The same native estimator and stated generating models at every design",
      x = "Participants", y = NULL, color = "Item pool", shape = "Item pool",
      caption = "2,048 planned and 2,048 null replications per design. Pointwise 95% Monte Carlo intervals; no familywise adjustment.\nDetection uses all scheduled runs; coverage uses valid fits. Horizontal lines mark nominal 5% and 95%.") +
    theme_minimal(base_size = 11) + theme(legend.position = "bottom", panel.grid.minor = element_blank(),
      plot.caption = element_text(hjust = 0))
  ggsave("validation/l2_replication_large_diagnostics.png", d, width = 13, height = 6.4, dpi = 160, bg = "white")
}
