# Summarize the item/person trade-off and anchored-DIF comparisons.
args <- commandArgs(TRUE)
out <- if (length(args)) args[1L] else "validation/rasch_fit_benchmark"
library(ggplot2)

read_panel <- function(task, metrics, profiles) {
  a <- read.csv(file.path(out, paste0(task, "_summary.csv")))
  a <- a[a$metric %in% metrics & a$profile %in% profiles & a$engine == "TAM", ]
  a$layout <- factor(paste0(a$n, " × ", a$items_per_person),
                     levels = c("120 × 16", "240 × 16", "480 × 16", "480 × 8"))
  a
}
cal <- read_panel("calibration", c("b_rmse", "person_eap_rel"), c("targeted", "narrow", "off_target"))
cal$measure <- factor(ifelse(cal$metric == "b_rmse", "Item difficulty RMSE (logits)",
                            "Person EAP reliability"),
                     levels = c("Item difficulty RMSE (logits)", "Person EAP reliability"))
cal$profile <- factor(cal$profile, levels = c("targeted", "narrow", "off_target"),
                      labels = c("Targeted", "Narrow ability range", "Harder bank"))
d <- read_panel("dif", "rejected", c("group_impact", "group_impact_and_dif"))
d$profile <- factor(d$profile, levels = c("group_impact", "group_impact_and_dif"),
                    labels = c("No DIF", "+0.70 logit DIF"))
base <- theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(), strip.text = element_text(face = "bold"),
        plot.title.position = "plot", plot.caption.position = "plot")
palette <- c("#1D6890", "#BC572C", "#557047")
position <- position_dodge(width = .34)
series <- function(data) {
  ggplot(data, aes(layout, estimate, group = profile, color = profile)) +
    geom_line(position = position, linewidth = .65, alpha = .75) +
    geom_pointrange(aes(ymin = mc_lower, ymax = mc_upper), position = position, linewidth = .55) +
    scale_color_manual(values = palette) +
    labs(x = "Persons × items per person", color = NULL, y = NULL) + base
}
calibration <- series(cal) + facet_wrap(~measure, scales = "free_y") +
  labs(title = "Item calibration and person measurement respond to different design changes",
       subtitle = sprintf("Fixed 16-item bank · %d simulations per condition", max(cal$scheduled)),
       caption = "TAM estimates with common diagnostics. Bars are 95% Monte Carlo intervals.")
dif <- series(d) + geom_hline(yintercept = .05, color = "grey50", linetype = "dashed") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2),
                     labels = function(x) paste0(round(100 * x), "%")) +
  labs(title = "Uniform DIF detection across sample sizes and layouts",
       subtitle = sprintf("Group mean difference 0.60 logits · 15 invariant anchors · %d simulations per condition",
                          max(d$scheduled)),
       y = "Rejection rate", caption = "TAM likelihood-ratio test. Bars are exact 95% Monte Carlo intervals.")
ggsave(file.path(out, "calibration.png"), calibration, width = 8.6, height = 4.4, dpi = 160, bg = "white")
ggsave(file.path(out, "dif.png"), dif, width = 7.8, height = 4.2, dpi = 160, bg = "white")
