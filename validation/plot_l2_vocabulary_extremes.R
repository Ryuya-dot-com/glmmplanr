# Plot the effect-scale comparison and the slope sensitivity of a paired layout contrast.
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
summary <- read.csv(paste0(args[1L], "_summary.csv"))
paired <- read.csv(paste0(args[1L], "_paired.csv"))
fixed <- subset(summary, id == "96_48_16" & model == "unadjusted" &
                  phase == "global" & effect == "pp5")
constant_logit <- rbind(subset(summary, id == "96_48_16" & phase == "global_logit"),
                        subset(fixed, p_restudy == .42))
slope_labels <- c(balanced = "Balanced", learner_dominant = "Learner-dominant",
                  word_dominant = "Word-dominant")
width <- subset(paired, comparison == "larger_bank" & effect == "pp5" &
                  measure == "mean_width_pp" & p_restudy %in% c(.10, .90) &
                  slopes %in% names(slope_labels))
stopifnot(nrow(fixed) == 3L, nrow(constant_logit) == 3L, nrow(width) == 6L)

png(args[2L], width = 1800, height = 900, res = 150, bg = "white")
par(family = "sans", mfrow = c(1L, 2L), mar = c(5.1, 5.1, 4.4, 1.0), oma = c(.2, 0, .2, .5),
    las = 1, bty = "l", col.axis = "#243647", col.lab = "#243647", col.main = "#243647")
blue <- "#236B9A"
orange <- "#AC5223"
plot(c(.10, .90), c(0, 100), type = "n", xlim = c(.05, .96), xaxt = "n", yaxt = "n",
     xlab = "True restudy accuracy", ylab = "Detection rate (%)")
axis(1, at = c(.10, .42, .90), labels = c("10%", "42%", "90%"))
axis(2, at = seq(0, 100, 20))
abline(h = seq(20, 80, 20), col = "#DDE4E9", lty = 3)
for (k in 1:2) {
  d <- list(fixed, constant_logit)[[k]]
  d <- d[order(d$p_restudy), ]
  x <- d$p_restudy + if (k == 1L) -.005 else .005
  y <- 100 * d$detected / d$attempted
  col <- c(blue, orange)[k]
  segments(x, 100 * d$detection_lo, x, 100 * d$detection_hi, col = col, lwd = 2)
  lines(x, y, col = col, lty = k, lwd = 2)
  points(x, y, pch = c(16L, 17L)[k], col = col, cex = 1.2)
}
legend("top", legend = c("Fixed 5 percentage points", "Fixed log odds ratio"),
       col = c(blue, orange), lty = 1:2, pch = c(16, 17), bty = "n", cex = .84)
title(main = "A  Effect scale changes detection", line = 2.5, adj = 0)
mtext("96 learners, 48-word bank, 16 words per learner", side = 3, line = 1.1, adj = 0, cex = .83)

par(mar = c(5.1, 7.6, 4.4, 1.0))
plot(range(c(width$lo, width$hi, 0)) + c(-.04, .04), c(.65, 3.4), type = "n", yaxt = "n",
     xlab = "Difference in mean 95% interval width (pp)", ylab = "")
axis(2, at = 3:1, labels = unname(slope_labels), tick = FALSE, line = -.45, cex.axis = .85)
abline(h = 1:3, col = "#DDE4E9", lty = 3)
abline(v = 0, col = "#7A8790", lty = 2)
for (k in 1:2) {
  d <- subset(width, p_restudy == c(.10, .90)[k])
  y <- 4L - match(d$slopes, names(slope_labels)) + c(.12, -.12)[k]
  col <- c(blue, orange)[k]
  segments(d$lo, y, d$hi, y, col = col, lwd = 2)
  points(d$difference, y, pch = c(16, 17)[k], col = col, cex = 1.2)
}
legend("bottomright", legend = c("10% restudy", "90% restudy"),
       col = c(blue, orange), pch = c(16, 17), bty = "n", cex = .84)
title(main = "B  Precision gained by a larger word bank", line = 2.5, adj = 0)
mtext("48-word bank minus 16-word bank, fixed 5-point effect", side = 3, line = 1.1, adj = 0, cex = .83)
dev.off()
