# Independently recompute decisions and retain every attempted outer dataset.
library(glmmplanr)
args <- commandArgs(TRUE)
out <- if (length(args)) args[1] else "validation/l2_default_bootstrap_confirm"
m <- readRDS(file.path(out, "manifest.rds"))$settings
z <- read.csv(file.path(out, "fits.csv.gz"))
keys <- c("case", "condition", "replicate")
stopifnot(nrow(z) == 6L * (2L * m$nsim + m$n_sensitivity),
  !anyDuplicated(z[c(keys, "method")]))
z$covered <- z$valid & abs(z$estimate - z$truth) <= qnorm(.975) * z$se
z$detected <- z$valid & abs(z$estimate / z$se) > qnorm(.975)
z$width <- 2 * qnorm(.975) * z$se
z$width_met <- z$valid & z$width <= 10
z$joint <- z$detected & z$width_met
stopifnot(!anyNA(z[c("covered", "detected", "joint")]))
ci <- function(k, n) if (n) unname(glmmplanr:::wilson_mc_ci(k, n)) else c(NA_real_, NA_real_)
gate <- function(interval, value, lower = TRUE) {
  if (anyNA(interval)) return("indeterminate")
  if (lower) {
    if (interval[1] >= value) "pass" else if (interval[2] < value) "fail" else "indeterminate"
  } else {
    if (interval[2] <= value) "pass" else if (interval[1] > value) "fail" else "indeterminate"
  }
}
summary <- do.call(rbind, lapply(split(z, interaction(z$case, z$condition, z$method)), function(v) {
  ok <- v$valid; n <- nrow(v); nv <- sum(ok)
  vc <- ci(nv, n); cc <- ci(sum(v$covered), nv); dc <- ci(sum(v$detected), n)
  jc <- ci(sum(v$joint), n); wc <- ci(sum(v$width_met), n)
  data.frame(v[1, c("profile", "n_subj", "n_item", "case", "condition", "method")],
    nsim = n, n_valid = nv, valid_rate = nv / n, valid_lo = vc[1], valid_hi = vc[2],
    n_covered = sum(v$covered), coverage = sum(v$covered) / nv, coverage_lo = cc[1], coverage_hi = cc[2],
    n_detected = sum(v$detected), detection = mean(v$detected), detection_lo = dc[1], detection_hi = dc[2],
    detection_given_valid = sum(v$detected) / nv, n_joint = sum(v$joint),
    joint = mean(v$joint), joint_lo = jc[1], joint_hi = jc[2],
    width_attainment = mean(v$width_met), width_lo = wc[1], width_hi = wc[2],
    bias_pp = mean(v$estimate[ok] - v$truth[ok]), empirical_sd_pp = sd(v$estimate[ok]),
    rms_se_pp = sqrt(mean(v$se[ok]^2)), mean_width_pp = mean(v$width[ok]),
    inner_requested = sum(v$n_boot), inner_valid = sum(v$n_inner_valid),
    inner_valid_min = if (all(v$n_boot == 0)) NA_real_ else min(v$n_inner_valid / v$n_boot),
    valid_threshold_met = nv / n >= .99, coverage_floor = gate(cc, .925),
    null_ceiling = if (v$condition[1] == "null") gate(dc, .075, FALSE) else NA_character_,
    joint_decision = gate(jc, .8),
    nominal_adverse_flag = (cc[2] < .95) || (v$condition[1] == "null" && dc[1] > .05))
}))
write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)
write.csv(subset(z, !valid), file.path(out, "failures.csv"), row.names = FALSE)
paired_summary <- function(a, b) {
  p <- merge(subset(z, method == a), subset(z, method == b), by = keys, suffixes = c("_a", "_b"))
  stopifnot(all(p$response_hash_a == p$response_hash_b), all(p$truth_a == p$truth_b))
  do.call(rbind, lapply(split(p, interaction(p$case, p$condition)), function(v) {
    ok <- v$valid_a & v$valid_b
    stopifnot(max(abs(v$estimate_a[ok] - v$estimate_b[ok])) < 1e-7)
    dc <- as.numeric(v$covered_b) - as.numeric(v$covered_a)
    dd <- as.numeric(v$detected_b) - as.numeric(v$detected_a)
    dj <- as.numeric(v$joint_b) - as.numeric(v$joint_a)
    data.frame(case = v$case[1], condition = v$condition[1], method_a = a, method_b = b,
      n = nrow(v), n_both_valid = sum(ok), validity_changes = sum(v$valid_a != v$valid_b),
      rms_se_ratio = sqrt(mean(v$se_b[ok]^2) / mean(v$se_a[ok]^2)),
      se_ratio_q025 = unname(quantile(v$se_b[ok] / v$se_a[ok], .025)),
      se_ratio_q975 = unname(quantile(v$se_b[ok] / v$se_a[ok], .975)),
      coverage_change_both_valid = mean(dc[ok]), coverage_mcse = sd(dc[ok]) / sqrt(sum(ok)),
      detection_change = mean(dd), detection_mcse = sd(dd) / sqrt(nrow(v)),
      joint_change = mean(dj), joint_mcse = sd(dj) / sqrt(nrow(v)),
      n_joint_changes = sum(v$joint_a != v$joint_b))
  }))
}
paired <- rbind(paired_summary("joint_delta", "bootstrap_99"),
  paired_summary("bootstrap_99", "bootstrap_399"))
write.csv(paired, file.path(out, "paired.csv"), row.names = FALSE)

# Audit the retained vector lengths, failed inner fits, and recomputed SEs.
files <- list.files(out, pattern = "^job-.*[.]rds$", full.names = TRUE)
n_executed <- n_inner_valid <- 0L
for (file in files) {
  r <- readRDS(file)
  for (id in names(r$bootstrap_draws)) {
    b <- r$bootstrap_draws[[id]]
    n_executed <- n_executed + length(b); n_inner_valid <- n_inner_valid + sum(is.finite(b))
    v <- subset(r$trials, replicate == as.integer(id) & method != "joint_delta")
    for (i in seq_len(nrow(v))) {
      w <- head(b, v$n_boot[i]); good <- sum(is.finite(w))
      stopifnot(good == v$n_inner_valid[i])
      if (v$valid[i]) stopifnot(abs(100 * sd(w[is.finite(w)]) - v$se[i]) < 1e-7)
    }
  }
}
inner <- data.frame(n_outer = 6L * m$nsim, n_executed = n_executed, n_inner_valid = n_inner_valid,
  n_inner_failed = n_executed - n_inner_valid, duplicated_probe_attempts = 6L * 2L * 99L)
write.csv(inner, file.path(out, "inner_execution.csv"), row.names = FALSE)
library(ggplot2)
main <- subset(summary, method != "bootstrap_399")
planned <- subset(main, condition == "planned"); null <- subset(main, condition == "null")
panel <- function(v, metric, rate, lo, hi) {
  data.frame(case = v$case, method = v$method, metric = metric,
    rate = v[[rate]], lo = v[[lo]], hi = v[[hi]])
}
pd <- rbind(panel(planned, "Planned 95% interval coverage", "coverage", "coverage_lo", "coverage_hi"),
  panel(null, "Null rejection", "detection", "detection_lo", "detection_hi"),
  panel(planned, "Detection AND width <=10pp", "joint", "joint_lo", "joint_hi"))
labels <- c(`31` = "Intercepts only\nN=54, J=30", `50` = "Moderate slopes\nN=72, J=60",
  `25` = "Larger variances\nN=162, J=30")
pd$design <- factor(labels[as.character(pd$case)], levels = rev(labels))
pd$method <- factor(pd$method, c("joint_delta", "bootstrap_99"), c("Joint delta", "Default bootstrap (99)"))
reference <- data.frame(metric = c("Planned 95% interval coverage", "Null rejection",
  "Detection AND width <=10pp"), value = c(.95, .05, .8))
p <- ggplot(pd, aes(rate, design, color = method)) +
  geom_vline(data = reference, aes(xintercept = value), color = "grey55", linetype = "dashed") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = .15,
    position = position_dodge(width = .5)) + geom_point(position = position_dodge(width = .5), size = 2.5) +
  facet_wrap(~metric, scales = "free_x") +
  scale_x_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  scale_color_manual(values = c("Joint delta" = "#6D7784", "Default bootstrap (99)" = "#B74627")) +
  labs(x = NULL, y = NULL, color = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(out, "comparison.png"), p, width = 12, height = 5.5, dpi = 150)
print(summary[c("case", "condition", "method", "n_valid", "coverage", "detection", "joint",
  "coverage_floor", "null_ceiling", "joint_decision")], row.names = FALSE)
print(paired, row.names = FALSE)
print(inner, row.names = FALSE)
cat("All outer keys, paired hashes, denominators, and inner SEs audited.\n")
