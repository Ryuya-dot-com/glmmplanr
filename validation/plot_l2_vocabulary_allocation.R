library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
dir <- if (length(args)) args[1L] else "validation"
data <- utils::read.csv(file.path(dir, "l2_vocabulary_summary.csv"))
data <- subset(data, phase == "core" & model == "maximal" & effect_pp == 8)
data$allocation <- factor(data$allocation, c("latin_square", "between_subj"),
                          c("Words rotate\nacross conditions", "Participants assigned\nto one condition"))
data$budget <- factor(data$responses, c(1536, 3072), c("1,536 responses", "3,072 responses"))
layouts <- with(data, unique(data.frame(n_subj, n_item)))
layouts <- layouts[order(layouts$n_subj, layouts$n_item), ]
data$layout <- factor(with(data, paste(n_subj, n_item, sep = " × ")),
                      with(layouts, paste(n_subj, n_item, sep = " × ")))
data$profile <- factor(data$profile,
                       c("no_slopes", "balanced", "learner_dominant", "word_dominant"),
                       c("No slopes", "Balanced", "Learner-dominant", "Word-dominant"))

figure <- ggplot(data, aes(layout, mean_width_pp, colour = profile, shape = profile,
                           group = profile)) +
  geom_line(linewidth = .8) +
  geom_point(size = 2.7) +
  facet_grid(allocation ~ budget, scales = "free", switch = "y") +
  scale_colour_manual(values = c("#626977", "#316592", "#a64028", "#137661")) +
  labs(x = "Participants × words", y = "Mean full 95% interval width (percentage points)",
       colour = "Generating slope SDs", shape = "Generating slope SDs",
       caption = paste("Planned effect of +8 percentage points. Each point uses 2,048 datasets. Vertical scales differ by row.",
                       "Balanced means SD 0.25 for each slope allowed by the allocation.", sep = "\n")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        strip.placement = "outside", strip.text.y.left = element_text(angle = 0),
        strip.text = element_text(face = "bold"), panel.spacing = grid::unit(1, "lines"),
        plot.caption = element_text(hjust = 0))
ggsave(file.path(dir, "l2_vocabulary_allocation.png"), figure,
       width = 11, height = 6.1, dpi = 160, bg = "white")
