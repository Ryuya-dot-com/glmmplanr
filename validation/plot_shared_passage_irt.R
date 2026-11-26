# Plot the shared-response population comparison with exact Monte Carlo intervals.
args <- commandArgs(TRUE)
out <- if (length(args)) args[1L] else "validation"
s <- read.csv(file.path(out, "shared_passage_irt_large_population_summary.csv"))
s <- subset(s, n == 48 & local_sd == .4 & slope_sd == .4)
null <- subset(s, effect_pp == 0)
planned <- subset(s, effect_pp == 10)
row <- function(s, count, denominator, lo, hi, label) {
  data.frame(passages = s$passages, model = s$model, value = 100 * s[[count]] / s[[denominator]],
             lo = 100 * s[[lo]], hi = 100 * s[[hi]], measure = label)
}
d <- rbind(
  row(null, "detected", "scheduled", "detection_lo", "detection_hi", "Null detection (%)"),
  row(planned, "covered", "intervals", "coverage_lo", "coverage_hi", "Coverage at +10 pp (%)"),
  row(planned, "detected", "scheduled", "detection_lo", "detection_hi", "Power at +10 pp (%)")
)
d$measure <- factor(d$measure, levels = unique(d$measure))
d$engine <- ifelse(grepl("^native", d$model), "glmmplanr", "glmmTMB")
d$effect <- ifelse(grepl("intercept$", d$model), "Passage intercept", "Passage intercept and slope")
ref <- data.frame(measure = levels(d$measure), value = c(5, 95, 80))
library(ggplot2)
g <- ggplot(d, aes(passages, value, colour = engine, shape = effect, linetype = effect)) +
  geom_hline(data = ref, aes(yintercept = value), inherit.aes = FALSE, colour = "#8C969D") +
  geom_line(linewidth = .65) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = .23, linewidth = .5) +
  geom_point(size = 2.4) +
  facet_wrap(vars(measure), scales = "free_y") +
  scale_x_continuous(breaks = c(2, 5, 12)) +
  scale_colour_manual(values = c(glmmplanr = "#006E75", glmmTMB = "#AD4D18"), name = NULL) +
  labs(x = "Passages shared across both arms", y = NULL, shape = NULL, linetype = NULL,
       title = "More passages improve population inference at the same response count",
       subtitle = "48 readers, 24 questions, 1,024 datasets per scenario, heterogeneous passage effects",
       caption = "Bars show exact 95% Monte Carlo intervals. Gray lines mark 5% null detection, 95% coverage, and 80% power.") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", plot.title.position = "plot", panel.grid.minor = element_blank(),
        panel.spacing.x = grid::unit(1.2, "lines"))
ggsave(file.path(out, "shared_passage_irt_population.png"), g, width = 10, height = 4.5, dpi = 170, bg = "white")
