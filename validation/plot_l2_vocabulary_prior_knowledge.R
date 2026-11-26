library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1L] else "validation"
stem <- "l2_vocabulary_prior_knowledge"
scenarios <- utils::read.csv(file.path(out, paste0(stem, "_scenarios.csv")))
exposure <- utils::read.csv(file.path(out, paste0(stem, "_exposure.csv")))
data <- merge(scenarios, exposure, by = "case")
data <- subset(data, effect_pp == 8 & profile %in% c("low", "strong", "ceiling") &
                 screen != "weaker" & policy %in% c("fixed_screened", "topup_screened"))
data$profile <- factor(data$profile, c("low", "strong", "ceiling"),
                       c("10% and 20% known", "10% and 50% known", "30% and 70% known"))
data$screen <- factor(data$screen, c("perfect", "imperfect"),
                      c("Perfect screening", "90% sensitivity, 95% specificity"))
data$screened <- data$upper_screened_n / 24
data$unknown <- data$upper_unknown_n / 24
data$policy <- factor(data$policy, c("fixed_screened", "topup_screened"),
                      c("Exclude from fixed list", "Top up from bank"))

figure <- ggplot(data, aes(screened, unknown, colour = screen, shape = policy)) +
  geom_hline(yintercept = 16, linewidth = .4, linetype = "dashed", colour = "#696e73") +
  geom_line(aes(group = screen), linewidth = .7, alpha = .78) +
  geom_point(size = 2.8) +
  facet_wrap(~profile, nrow = 1) +
  scale_colour_manual(values = c("#16627c", "#a9472d")) +
  scale_x_continuous(breaks = c(16, 24, 32, 40, 48), limits = c(15, 48)) +
  scale_y_continuous(breaks = seq(4, 16, by = 4), limits = c(3, 17)) +
  labs(x = "Mean words screened per participant", y = "Mean genuinely unknown words taught",
       colour = NULL, shape = NULL,
       caption = paste("Participants in the top vocabulary-size quartile, 2,048 datasets per scenario.",
                       "Facet headings give known-word rates at −1 and +1 vocabulary SD.",
                       sep = "\n")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", legend.box = "vertical", panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"), plot.caption = element_text(hjust = 0))
ggsave(file.path(out, paste0(stem, "_screening.png")), figure,
       width = 10.8, height = 5.2, dpi = 160, bg = "white")
