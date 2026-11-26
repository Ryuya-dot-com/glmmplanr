# Verify the archived counts and compare models on matched response data.
args <- commandArgs(TRUE)
out <- if (length(args)) args[1L] else "validation"
stage <- if (length(args) > 1L) match.arg(args[2L], c("pilot", "medium")) else "medium"
prefix <- file.path(out, paste0("few_shared_passages_", stage))
cases <- read.csv(paste0(prefix, "_scenarios.csv"))
fits <- read.csv(paste0(prefix, "_fits.csv.gz"))
summary <- read.csv(paste0(prefix, "_summary.csv"))
groups <- split(fits, paste(fits$case, fits$model))
z <- qnorm(.975)
exact <- function(x, n) if (n) unname(binom.test(x, n)$conf.int) else c(NA_real_, NA_real_)
same_number <- function(x, y) isTRUE(all.equal(x, y, tolerance = 1e-10, check.attributes = FALSE))

stopifnot(nrow(cases) == 44L, identical(cases$case, seq_len(nrow(cases))),
          !anyDuplicated(fits[c("case", "replication", "model")]),
          !anyDuplicated(summary[c("case", "model")]), length(groups) == nrow(summary),
          nrow(fits) == sum(4L * cases$n_rep + cases$n_local_check))
for (k in seq_len(nrow(summary))) {
  s <- summary[k, ]
  f <- groups[[paste(s$case, s$model)]]
  v <- f[f$valid, ]
  nv <- nrow(v)
  expected <- if (s$model == "fixed_interaction_local") s$n_local_check else s$n_rep
  cs <- sum(abs(v$estimate_pp - v$truth_sampled_pp) <= z * v$se_pp)
  cp <- sum(abs(v$estimate_pp - s$effect_pp) <= z * v$se_pp)
  det <- sum(abs(v$estimate_pp) > z * v$se_pp)
  stopifnot(identical(f$replication, seq_len(expected)), nrow(f) == s$attempted, nv == s$valid,
            cs == s$covered_sampled, cp == s$covered_population, det == s$detected,
            cs == sum(f$covered_sampled, na.rm = TRUE), cp == sum(f$covered_population, na.rm = TRUE),
            det == sum(f$detected), !any(f$detected[!f$valid]),
            same_number(exact(cs, nv), c(s$sampled_lo, s$sampled_hi)),
            same_number(exact(cp, nv), c(s$population_lo, s$population_hi)),
            same_number(exact(det, nrow(f)), c(s$detection_lo, s$detection_hi)),
            same_number(mean(v$estimate_pp - v$truth_sampled_pp), s$bias_sampled),
            same_number(mean(v$estimate_pp - s$effect_pp), s$bias_population),
            same_number(2 * z * mean(v$se_pp), s$mean_width))
  if (s$slope_sd == 0 && s$effect_pp == 0) stopifnot(all(f$truth_sampled_pp == 0))
}

shared <- expand.grid(condition = factor(c("A", "B")), passage = factor(1:2))
one_per_condition <- subset(shared, as.integer(condition) == as.integer(passage))
rank <- c(shared_additive = qr(model.matrix(~ condition + passage, shared))$rank,
          shared_interaction = qr(model.matrix(~ condition * passage, shared))$rank,
          one_per_condition_additive = qr(model.matrix(~ condition + passage, one_per_condition))$rank)
stopifnot(identical(unname(rank), c(3L, 4L, 2L)))

paired_rate_ci <- function(dx) {
  if (!length(dx)) return(c(difference = NA_real_, lower = NA_real_, upper = NA_real_))
  # The difference is the positive discordance rate minus the negative discordance rate.
  positive <- binom.test(sum(dx > 0), length(dx), conf.level = .975)$conf.int
  negative <- binom.test(sum(dx < 0), length(dx), conf.level = .975)$conf.int
  c(difference = mean(dx), lower = 100 * (positive[1L] - negative[2L]),
    upper = 100 * (positive[2L] - negative[1L]))
}

paired_mean_ci <- function(dx) {
  delta <- mean(dx)
  error <- z * sd(dx) / sqrt(length(dx))
  c(difference = delta, lower = delta - error, upper = delta + error)
}

compare <- function(k, left, right, coverage) {
  a <- groups[[paste(k, left)]]
  b <- groups[[paste(k, right)]]
  b <- b[match(a$replication, b$replication), ]
  stopifnot(identical(a$replication, b$replication), identical(a$truth_sampled_pp, b$truth_sampled_pp))
  v <- a$valid & b$valid
  stats <- rbind(detection_pp = paired_rate_ci(100 * (a$detected - b$detected)),
                 coverage_pp = paired_rate_ci(100 * (a[[coverage]][v] - b[[coverage]][v])),
                 width_pp = paired_mean_ci(2 * z * (a$se_pp[v] - b$se_pp[v])),
                 estimate_pp = paired_mean_ci(a$estimate_pp[v] - b$estimate_pp[v]))
  cbind(cases[k, ], left, right, coverage, pairs = nrow(a), valid_pairs = sum(v),
        measure = rownames(stats), stats, row.names = NULL)
}
paired <- lapply(cases$case, function(k) {
  contrasts <- rbind(compare(k, "fixed_interaction", "fixed_common", "covered_sampled"),
                     compare(k, "random_slope", "random_intercept", "covered_population"))
  if (cases$n_local_check[k] > 0L) {
    contrasts <- rbind(contrasts, compare(k, "fixed_interaction_local", "fixed_interaction", "covered_sampled"))
  }
  contrasts
})
write.csv(do.call(rbind, paired), paste0(prefix, "_paired.csv"), row.names = FALSE)
cat(sprintf("Verified %d scenarios, %d scenario-model summaries, and %d individual fits.\n",
            nrow(cases), nrow(summary), nrow(fits)))
print(rank)
