library(ggplot2)
x <- read.csv("validation/l2_targeting_pilot/summary.csv")
x <- subset(x, model == "Rasch" & engine == "native" & layout == "full50" &
  profile %in% c("reference", "lower", "lower_narrow", "severe_lower"))
x$profile <- factor(x$profile, c("reference", "lower", "lower_narrow", "severe_lower"),
  c("Reference", "Mean -3", "Mean -3; half SD", "Mean -5"))
p <- ggplot(x, aes(n_subj, attainment, color = profile, shape = profile)) +
  geom_hline(yintercept = .8, color = "grey55", linewidth = .4) +
  geom_errorbar(aes(ymin = mc_lo, ymax = mc_hi), width = .025, linewidth = .35) +
  geom_line(linewidth = .7) + geom_point(size = 2) + facet_wrap(~ format) +
  scale_color_manual(values = c("#0072B2", "#D55E00", "#CC79A7", "#333333")) +
  scale_x_continuous(trans = "log10", breaks = c(100, 200, 500, 1000)) +
  scale_y_continuous(limits = c(0, 1), labels = function(y) paste0(100 * y, "%")) +
  labs(title = "Participant count alone does not determine calibration accuracy",
    subtitle = "50 vocabulary items per format; Rasch refits to public responses define the simulated banks",
    x = "Participants (log scale)", y = "Whole-bank difficulty attainment", color = NULL, shape = NULL,
    caption = "64 synthetic replications per cell; Wilson 95% Monte Carlo intervals. Invalid banks count as misses.\nTarget: RMSE <= 0.35 and maximum item error <= 0.60 logits. Each format has its own reference scale.\nAbility shifts are sensitivity assumptions, not age/CEFR labels. The source study fitted 2PL, not Rasch.") +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom",
    panel.grid.minor = element_blank(), plot.caption = element_text(hjust = 0))
ggsave("validation/l2_targeting.png", p, width = 10, height = 5.6, dpi = 160, bg = "white")
