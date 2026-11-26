# Paired screening results; independent confirmations are reported separately.
library(ggplot2)
x <- read.csv("validation/item_pool_irtsim_screen/summary.csv")
x$form <- factor(x$form,
  c("central30", "spread30", "spread60", "full150", "rotating150"),
  c("Central 30", "Spread 30", "Spread 60", "Full 150", "150 pool; 30/person"))
x$engine <- factor(x$engine, c("native", "irtsim"), c("Native MML", "irtsim / mirt"))
x$ability <- factor(x$person_sd, c(.4, .8), c("Ability SD 0.4", "Ability SD 0.8"))
p <- ggplot(x, aes(n_subj, attainment, color = engine, linetype = engine, shape = engine)) +
  geom_hline(yintercept = .8, color = "grey50", linewidth = .4) +
  geom_errorbar(aes(ymin = mc_lo, ymax = mc_hi), width = .07, linewidth = .35) +
  geom_line(linewidth = .6) + geom_point(size = 2) +
  facet_grid(ability ~ form) +
  scale_x_continuous(trans = "log2", breaks = c(120, 480, 1920)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, .4, .8, 1),
    labels = function(z) paste0(100 * z, "%")) +
  scale_color_manual(values = c("#222222", "#D55E00")) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  scale_shape_manual(values = c(16, 1)) +
  labs(title = "Recalibration attainment depends on the retained form and item exposure",
    subtitle = "Same response matrices and Rasch scale for both estimators; 128 replications per condition",
    x = "Participants (log scale)", y = "Whole-bank joint attainment", color = NULL,
    linetype = NULL, shape = NULL,
    caption = "RMSE <= 0.35 and maximum item error <= 0.60 logits; invalid banks count as misses.\nBars: pointwise Wilson 95% Monte Carlo intervals. Independent 512-run checks determine confirmation.\nNear-overlapping curves show close agreement for these conditions, not universal estimator equivalence.") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom",
    panel.spacing.x = grid::unit(1.3, "lines"),
    panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0))
ggsave("validation/item_pool_comparison.png", p, width = 11, height = 6, dpi = 160, bg = "white")
