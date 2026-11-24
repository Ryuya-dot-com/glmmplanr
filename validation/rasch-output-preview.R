library(ggplot2)
library(patchwork)
library(glmmplanr)
# Run with the glmmplanr 0.1.0 working version installed.
# Rscript validation/rasch-output-preview.R /path/to/figure.png
args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args)) args[1L] else "rasch-output-preview.png"
limits <- c(difficulty = .35, difficulty_max_abs_error = .60)
fixed <- plan_rasch_fixed_bank(
  item_difficulties = c(-2, -1, -.6, -.2, .2, .6, 1, 2),
  n_subj = c(40L, 100L), person_sd = .8, max_error = limits,
  method = "fixed_item", engine = "mml", assurance = .80, nsim = 16L, seed = 417L
)
sparse <- plan_rasch_fixed_bank(
  item_difficulties = seq(-2, 2, length.out = 12L),
  n_subj = 8L, items_per_subject = 2L, item_seed = 1L,
  max_error = c(difficulty = .35), nsim = 16L, seed = 417L
)
figure <- (plot(fixed) / plot(sparse, type = "items")) +
  plot_layout(heights = c(1, 1.2))
ggsave(output, figure, width = 9, height = 8, dpi = 140, bg = "white")
