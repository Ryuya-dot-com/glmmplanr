library(ggplot2)
args <- commandArgs(TRUE)
directory <- if (length(args)) args[1L] else "validation/item_pool_planning_verified"
destination <- if (length(args) > 1L) args[2L] else "validation/item_pool_planning.png"
x <- read.csv(file.path(directory, "known_bank_curves.csv"))
x <- subset(x, abs(theta) <= 3.5)
curves <- rbind(transform(x, metric = "Conditional person RMSE (logits)", value = rmse),
                transform(x, metric = "Test information", value = information))
curves$form <- factor(curves$form, c("central30", "spread30", "spread60", "full150"),
                      c("Central 30", "Spread 30", "Spread 60", "Full 150"))
p <- ggplot(curves, aes(theta, value, color = form)) + geom_line(linewidth = .85) +
  facet_wrap(~ metric, ncol = 1L, scales = "free_y") +
  scale_color_manual(values = c("#D55E00", "#0072B2", "#009E73", "#333333")) +
  labs(title = "A shorter test changes which abilities it measures well",
       subtitle = "Synthetic 150-item Rasch bank; equal content quotas in every candidate",
       x = "True ability (logits)", y = NULL, color = NULL,
       caption = "Known item parameters. RMSE includes EAP bias under a fixed N(0,1) scoring prior.\nNumerical expectation over all scores; this is not item-recalibration accuracy.") +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom",
    strip.text = element_text(hjust = 0, face = "bold"),
    plot.caption = element_text(hjust = 0), panel.grid.minor = element_blank())
ggsave(destination, p, width = 8.4, height = 7.6, dpi = 160, bg = "white")
