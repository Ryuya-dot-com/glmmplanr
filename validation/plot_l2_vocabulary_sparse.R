library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1L] else "validation"
results <- utils::read.csv(file.path(out, "l2_vocabulary_sparse_summary.csv"))
show <- subset(results, effect_pp == 0 &
                 profile %in% c("match", "allocation_bias", "match_mar", "outcome_mnar"))
show$profile <- factor(show$profile, c("outcome_mnar", "match_mar", "allocation_bias", "match"),
                       c("Outcome-related missingness", "Compatibility-related missingness",
                         "Imbalanced allocation", "Balanced allocation"))
show$layout <- factor(show$layout, c("full", "linked", "disjoint"),
                      c("Full, 16 words", "Linked, 48 words", "Disjoint, 48 words"))
show$model <- factor(show$model, c("unadjusted", "standardized_match"),
                     c("Unadjusted", "Compatibility-adjusted"))
show$lower <- show$bias_pp - stats::qnorm(.975) * show$bias_mcse_pp
show$upper <- show$bias_pp + stats::qnorm(.975) * show$bias_mcse_pp
offset <- position_dodge(width = .53)
figure <- ggplot(show, aes(bias_pp, profile, colour = model, shape = model)) +
  geom_vline(xintercept = 0, colour = "#545c65", linewidth = .5) +
  geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = .27,
                position = offset, linewidth = .6) +
  geom_point(position = offset, size = 1.9) +
  facet_wrap(~layout, nrow = 1) +
  scale_colour_manual(values = c("#a2432f", "#1f6285")) +
  scale_x_continuous(breaks = seq(0, 15, by = 5)) +
  labs(x = "Bias in estimated recall difference (percentage points)", y = NULL,
       colour = NULL, shape = NULL,
       caption = paste("True average effect of zero. Each estimate uses 3,072 datasets.",
                       "Bars show 95% Monte Carlo intervals for mean bias.")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(colour = "#ebedef"), strip.text = element_text(face = "bold"),
        panel.spacing = grid::unit(.9, "lines"), plot.caption = element_text(hjust = 0))
ggsave(file.path(out, "l2_vocabulary_sparse_bias.png"), figure,
       width = 11, height = 4.4, dpi = 160, bg = "white")
