# Run from the repository root against the installed version of glmmplanr.
library(glmmplanr)

design <- design_crossed(list(A = c("a0", "a1"), B = c("b0", "b1"), C = c("c0", "c1")),
                         between_subj = "A", between_item = "B", attrition = .05)
dgp <- calibrate_dgp(design,
  cell_probs_from_diffs(design, .65, c(A = 8, B = 6, C = 4, "A:B" = 8, "A:B:C" = 8)),
  chance = 1/3, subj_prob_range = c(.42, .85), item_prob_range = c(.45, .86),
  slope_sd_ratio = .30, corr_structure = .20)
planned <- specify_target(dgp, "A:B:C", scale = "pp", criterion = "joint", power_min = .80,
                          ci_width = 12, engine = "bayes_laplace", prob_boot = 19, on_singular = "flag_only")
specs <- list(null = glmmplanr:::null_target_spec(planned), planned = planned)
scenarios <- data.frame(scenario = c("natural", "perfect_person", "zero_item", "both_item_priority"),
                        persons = c(FALSE, TRUE, FALSE, TRUE), items = c(FALSE, FALSE, TRUE, TRUE))
# Sample exactly 1/16 of the relevant units, independently of allocation and responses.
fraction <- 1/16
n_fit <- 12L                 # worst-case binomial MCSE 0.145 for downstream rates
n_frequency <- 400L          # worst-case MCSE 0.025 for natural pattern prevalence
cases <- expand.grid(n_subj = c(32L, 64L), condition = names(specs), stringsAsFactors = FALSE)
cases$n_item <- 1024L %/% cases$n_subj

count_extremes <- function(id, y, n) {
  attempts <- tabulate(id, nbins = n)
  correct <- tabulate(id[y == 1L], nbins = n)
  c(perfect = sum(attempts > 0L & correct == attempts),
    zero = sum(attempts > 0L & correct == 0L), unobserved = sum(attempts == 0L))
}

trial_counts <- function(d, ns, ni) {
  s <- count_extremes(d$subj_id, d$y, ns)
  i <- count_extremes(d$item_id, d$y, ni)
  as.data.frame(as.list(c(setNames(s, paste0("person_", names(s))), setNames(i, paste0("item_", names(i))))))
}

stopifnot(identical(unname(count_extremes(c(1L, 1L, 2L), c(1L, 1L, 0L), 3L)), c(1L, 1L, 1L)))

run_case <- function(j) {
  case <- cases[j, ]
  spec <- specs[[case$condition]]
  ns <- case$n_subj
  ni <- case$n_item
  schedule <- generate_design_schedule(design, ns, ni)
  stopifnot(nrow(schedule$trials) == 1024L, abs(sum(spec$contrast_weights)) < 1e-9)
  frequencies <- vector("list", n_frequency)
  fits <- vector("list", n_fit * nrow(scenarios))
  for (r in seq_len(n_frequency)) {
    seed <- 731000L + 1000L * j + r
    set.seed(seed, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion", sample.kind = "Rejection")
    base <- simulate_dgp_data(spec$dgp, schedule)
    counts <- trial_counts(base, ns, ni)
    frequencies[[r]] <- data.frame(case, replication = r, seed, n_obs = nrow(base), counts)
    if (r > n_fit) next
    set.seed(seed + 10000L)
    perfect <- sample.int(ns, ns * fraction)
    zero <- sample.int(ni, ni * fraction)
    for (k in seq_len(nrow(scenarios))) {
      sc <- scenarios[k, ]
      d <- base
      if (sc$persons) d$y[d$subj_id %in% perfect] <- 1L
      if (sc$items) d$y[d$item_id %in% zero] <- 0L
      # An always-wrong item takes precedence at row/column intersections.
      # The resulting target is the population contrast averaged over class assignments.
      p <- spec$dgp$achieved_cell_probs
      if (sc$persons) p <- fraction + (1 - fraction) * p
      if (sc$items) p <- (1 - fraction) * p
      target <- 100 * sum(spec$contrast_weights * p)
      stopifnot(abs(target - spec$true_effect_target * (1 - fraction)^(sc$persons + sc$items)) < 1e-6)
      set.seed(seed + 20000L)   # paired bootstrap streams within each generated dataset
      fit <- glmmplanr:::fit_single_dataset(d, spec)
      valid <- isTRUE(fit$converged) && identical(fit$status, "converged") &&
        is.finite(fit$se) && fit$se > 0 && is.finite(fit$est) &&
        is.finite(fit$prob_est) && is.finite(fit$prob_se) && fit$prob_se > 0
      estimate <- if (valid) 100 * fit$prob_est else NA_real_
      half_width <- if (valid) qnorm(1 - spec$alpha / 2) * 100 * fit$prob_se else NA_real_
      detected <- valid && abs(estimate) > half_width
      counts_fit <- trial_counts(d, ns, ni)
      stopifnot(nrow(d) == nrow(base), identical(counts_fit$person_unobserved, counts$person_unobserved),
                identical(counts_fit$item_unobserved, counts$item_unobserved))
      fits[[(r - 1L) * nrow(scenarios) + k]] <- data.frame(case, replication = r, seed,
        scenario = sc$scenario, forced_persons = length(perfect) * sc$persons, forced_items = length(zero) * sc$items,
        n_obs = nrow(d), changes = sum(d$y != base$y), counts_fit, true_pp = target,
        estimated_pp = estimate, width_pp = 2 * half_width, valid, detected,
        covered = if (valid) abs(estimate - target) <= half_width else NA,
        joint = detected && 2 * half_width <= spec$ci_width, singular = isTRUE(fit$singular),
        status = fit$status, diagnostic = fit$error_reason)
    }
  }
  list(frequency = do.call(rbind, frequencies), fits = do.call(rbind, fits))
}

workers <- if (.Platform$OS.type == "windows") 1L else 4L
runs <- parallel::mclapply(seq_len(nrow(cases)), run_case, mc.cores = workers, mc.set.seed = FALSE)
stopifnot(all(vapply(runs, is.list, logical(1))))
raw <- do.call(rbind, lapply(runs, "[[", "fits"))
frequencies <- do.call(rbind, lapply(runs, "[[", "frequency"))
stopifnot(nrow(raw) == nrow(cases) * n_fit * nrow(scenarios),
          nrow(frequencies) == nrow(cases) * n_frequency)
grouped <- function(d, keys) split(d, do.call(interaction, c(d[keys], list(drop = TRUE))))
keys <- c("n_subj", "condition", "n_item")
freq_summary <- do.call(rbind, lapply(grouped(frequencies, keys), function(d) {
  do.call(rbind, lapply(setdiff(names(d), c(keys, "replication", "seed", "n_obs")), function(name) {
    n <- sum(d[[name]] > 0)
    ci <- binom.test(n, nrow(d))$conf.int
    data.frame(d[1, keys], pattern = name, generated = nrow(d), datasets_with_pattern = n,
               prevalence = n/nrow(d), lo = ci[1], hi = ci[2], total_units_with_pattern = sum(d[[name]]))
  }))
}))
keys <- c(keys, "scenario")
summary_rows <- do.call(rbind, lapply(grouped(raw, keys), function(d) {
  n <- nrow(d)
  usable <- subset(d, valid)
  interval <- function(x, total) if (total > 0) binom.test(x, total)$conf.int else c(NA_real_, NA_real_)
  detection_ci <- interval(sum(d$detected), n)
  coverage_ci <- interval(sum(usable$covered), nrow(usable))
  joint_ci <- interval(sum(d$joint), n)
  err <- usable$estimated_pp - usable$true_pp
  data.frame(d[1, keys], true_pp = d$true_pp[1], attempted = n, valid = nrow(usable),
    detected = sum(d$detected), detected_lo = detection_ci[1], detected_hi = detection_ci[2],
    covered = sum(usable$covered), coverage_lo = coverage_ci[1], coverage_hi = coverage_ci[2],
    joint = sum(d$joint), joint_lo = joint_ci[1], joint_hi = joint_ci[2],
    mean_width_pp = if (nrow(usable)) mean(usable$width_pp) else NA_real_,
    bias_pp = if (length(err)) mean(err) else NA_real_, rmse_pp = if (length(err)) sqrt(mean(err^2)) else NA_real_,
    singular = sum(d$singular), R_version = as.character(getRversion()),
    package_version = as.character(packageVersion("glmmplanr")))
}))
write.csv(raw, "validation/extreme_fits.csv", row.names = FALSE)
write.csv(freq_summary, "validation/extreme_frequency.csv", row.names = FALSE)
write.csv(summary_rows, "validation/extreme_summary.csv", row.names = FALSE)
print(summary_rows[c(keys, "attempted", "valid", "detected", "covered", "joint", "mean_width_pp")], row.names = FALSE)

# Refit the original responses for two planned datasets per layout to diagnose
# inner iteration failures. This checks convergence of a single fit, not interval calibration.
iteration_cases <- subset(raw, condition == "planned" & scenario == "zero_item" & replication <= 2L)
iteration_rows <- lapply(seq_len(nrow(iteration_cases)), function(j) {
  case <- iteration_cases[j, ]
  ns <- case$n_subj
  ni <- case$n_item
  set.seed(case$seed, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion", sample.kind = "Rejection")
  base <- simulate_dgp_data(planned$dgp, n_subj = ns, n_item = ni)
  set.seed(case$seed + 10000L)
  sampled_persons <- sample.int(ns, ns * fraction) # reproduce draw preceding the item assignment
  zero <- sample.int(ni, ni * fraction)
  y <- base$y
  y[base$item_id %in% zero] <- 0L
  stopifnot(nrow(base) == case$n_obs, sum(y != base$y) == case$changes)
  cid <- base$cell_id
  controls <- data.frame(chance = c(1/3, 1/3, 0), max_inner = c(44L, 400L, 44L))
  do.call(rbind, lapply(seq_len(nrow(controls)), function(h) {
    ctl <- controls[h, ]
    f <- glmmplanr:::cpp_fit_crossed_glmm(y, design$X_cells[cid, , drop = FALSE],
      as.integer(factor(base$subj_id)) - 1L, design$Z_subj_cells[cid, , drop = FALSE], ns,
      as.integer(factor(base$item_id)) - 1L, design$Z_item_cells[cid, , drop = FALSE], ni,
      integer(), matrix(0, nrow(base), 0L), 0L, ctl$chance, 1L, TRUE,
      planned$lkj_eta, planned$prior_beta_sd, planned$prior_tau_sd, 400L, ctl$max_inner, 1.5e-4)
    if (h == 1L) stopifnot(identical(f$status, case$status))
    data.frame(n_subj = ns, n_item = ni, seed = case$seed, replication = case$replication,
               ctl, status = f$status, outer_iterations = f$iterations, final_inner_step = f$inner_step)
  }))
})
write.csv(do.call(rbind, iteration_rows), "validation/extreme_iterations.csv", row.names = FALSE)
