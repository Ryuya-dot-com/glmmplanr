# Joint design assurance under broad and hard anchor placement.
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
x <- read.csv(args[1L])
d <- subset(x, phase == "fixed_workload" & focal_mean == -.75 & drift == "none")
stopifnot(nrow(d) == 60L)
d$anchor_count <- factor(d$anchors)
d$placement <- factor(d$band, levels = c("broad", "hard"), labels = c("Broad anchors", "Hard anchors"))
d$reference_bank <- factor(d$reference, levels = c("oracle", "estimated"),
                            labels = c("Known anchor difficulties", "500-person reference calibration"))
d$attainment <- d$achieved / d$attempted
p <- ggplot2::ggplot(d, ggplot2::aes(n, attainment, colour = anchor_count, shape = anchor_count)) +
  ggplot2::geom_hline(yintercept = .8, linetype = "dashed", colour = "grey50") +
  ggplot2::geom_linerange(ggplot2::aes(ymin = mc_lo, ymax = mc_hi), linewidth = .65) +
  ggplot2::geom_line(linewidth = .7) + ggplot2::geom_point(size = 2.2) +
  ggplot2::facet_grid(reference_bank ~ placement) +
  ggplot2::scale_x_continuous(trans = "log2", breaks = sort(unique(d$n))) +
  ggplot2::scale_y_continuous(limits = c(0, 1), labels = function(v) paste0(round(v * 100), "%")) +
  ggplot2::scale_colour_manual(values = c("#236B9A", "#AC5223", "#357259")) +
  ggplot2::labs(title = "Joint accuracy of Rasch bank linking",
                subtitle = "16 items per examinee. New cohort mean is 0.75 logits below the reference.",
                x = "New examinees", y = "Calibration studies meeting all three tolerances",
                colour = "Anchors per examinee", shape = "Anchors per examinee",
                caption = "Bars are exact 95% Monte Carlo intervals. The dashed line is the 80% planning threshold.") +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(legend.position = "bottom", plot.caption = ggplot2::element_text(hjust = 0),
                 plot.caption.position = "plot")
ggplot2::ggsave(args[2L], p, width = 9, height = 6, dpi = 160, bg = "white")
