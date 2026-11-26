# Fixed before the full run. See validation/README.md for the estimand and gates.
library(glmmplanr)

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1L] else "validation"
n_rep <- if (length(args) > 1L) as.integer(args[2L]) else 1024L
stopifnot(dir.exists(out), !is.na(n_rep), n_rep > 0L, n_rep <= 1024L)
workers <- if (.Platform$OS.type == "windows") 1L else 4L
n_subj <- 48L
n_item <- 24L
material_var <- 0.41
z_crit <- stats::qnorm(0.975)
master_seed <- if (n_rep == 1024L) 83251001L else 83252001L

layouts <- data.frame(
  passages = c(24L, 12L, 8L, 6L, 4L, 2L, 3L),
  assignment = c(rep("passage", 6L), "participant")
)
layouts$questions <- n_item %/% layouts$passages
base <- data.frame(profile = "both_constant", subj_sd = 0.7, passage_sd = 0.5,
                   interaction_sd = 0.4, baseline = 0.60, effect = 0.10)
core <- rbind(
  transform(base, profile = "independent", passage_sd = 0, interaction_sd = 0),
  transform(base, profile = "passage_only", interaction_sd = 0),
  transform(base, profile = "interaction_only", passage_sd = 0),
  base,
  transform(base, profile = "both_size_dependent")
)
sensitivity <- rbind(
  transform(base, profile = "passage_low", passage_sd = 0.25),
  transform(base, profile = "passage_high", passage_sd = 0.60),
  transform(base, profile = "interaction_low", interaction_sd = 0.20),
  transform(base, profile = "interaction_high", interaction_sd = 0.60),
  transform(base, profile = "subj_low", subj_sd = 0.40),
  transform(base, profile = "subj_high", subj_sd = 1.00),
  transform(base, profile = "ceiling", baseline = 0.80),
  transform(base, profile = "smaller_effect", effect = 0.05),
  transform(base, profile = "both_high", passage_sd = 0.60, interaction_sd = 0.60)
)
scenarios <- rbind(
  merge(layouts, core),
  merge(subset(layouts, passages %in% c(6L, 12L)), sensitivity)
)
scenarios <- merge(scenarios, data.frame(condition = c("null", "planned")))
scenarios$p_a <- scenarios$baseline
scenarios$p_b <- scenarios$baseline + ifelse(scenarios$condition == "null", 0, scenarios$effect)
by_size <- scenarios$profile == "both_size_dependent"
scenarios$interaction_sd[by_size] <- 0.4 * sqrt((scenarios$questions[by_size] - 1) / 3)
# A one-item passage has the same global difficulty term as the item itself.
scenarios$passage_sd[scenarios$questions == 1L] <- 0
scenarios$item_sd <- sqrt(material_var - scenarios$passage_sd^2)
keys <- with(scenarios, paste(assignment, passages, condition, p_a, p_b,
                             subj_sd, passage_sd, interaction_sd))
scenarios$cell <- match(keys, unique(keys))
cases <- scenarios[!duplicated(scenarios$cell), ]
stopifnot(nrow(scenarios) == 106L, nrow(cases) == 96L)

calibrate_logit <- function(p, sd_total) {
  mean_response <- function(mu) {
    stats::integrate(function(z) stats::plogis(mu + sd_total * z) * stats::dnorm(z),
                     -Inf, Inf, rel.tol = 1e-10)$value
  }
  mu <- stats::uniroot(function(x) mean_response(x) - p, c(-20, 20), tol = 1e-11)$root
  stopifnot(abs(mean_response(mu) - p) < 1e-9)
  mu
}

streams <- glmmplanr:::allocate_stage_point_streams(master_seed, 1L, 1:2, 2L)
streams <- setNames(lapply(streams, glmmplanr:::derive_rep_substreams, nsim = n_rep),
                    c("null", "planned"))

# Positions within each arm associate the same item and passage draws across layouts.
arm_slots <- function(arm) {
  slots <- integer(length(arm))
  slots[arm == "A"] <- seq_len(sum(arm == "A"))
  slots[arm == "B"] <- n_item / 2L + seq_len(sum(arm == "B"))
  slots
}

run_case <- function(k) {
  cfg <- cases[k, ]
  factors <- list(condition = c("A", "B"))
  design <- if (cfg$questions == 1L) {
    design_crossed(factors, allocation = "between_item", contrast = "treatment")
  } else if (cfg$assignment == "passage") {
    design_testlet(factors, testlet_size = cfg$questions, between_item = "condition",
                   contrast = "treatment")
  } else {
    design_testlet(factors, testlet_size = cfg$questions, between_subj = "condition",
                   contrast = "treatment")
  }
  probs <- c(A = cfg$p_a, B = cfg$p_b)
  dgp <- calibrate_dgp(
    design, probs, chance = 0, subj_sd = cfg$subj_sd, item_sd = cfg$item_sd,
    slope_sd_ratio = 0, control = calibrate_control(
      testlet_sd = cfg$passage_sd,
      participant_testlet_sd = if (design$has_testlet) cfg$interaction_sd else 0
    )
  )
  spec <- specify_target(
    dgp, "condition", scale = "pp", engine = "bayes_laplace", analysis_re = "intercept_only",
    prob_se_method = "joint_delta", control = target_control(on_singular = "flag_only")
  )
  sd_total <- sqrt(cfg$subj_sd^2 + material_var + cfg$interaction_sd^2)
  mu <- vapply(probs, calibrate_logit, numeric(1), sd_total = sd_total)
  truth <- 100 * (cfg$p_b - cfg$p_a)
  stopifnot(abs(spec$true_effect_target - truth) < 1e-7)
  if (design$has_testlet || cfg$interaction_sd == 0) {
    stopifnot(max(abs(mu - dgp$mu_cells)) < 1e-7)
  }

  trials <- generate_design_schedule(design, n_subj, n_item)$trials
  passage_id <- if (design$has_testlet) trials$testlet_id else trials$item_id
  is_a <- trials$condition == "A"
  stopifnot(nrow(trials) == n_subj * n_item, sum(is_a) == nrow(trials) / 2L)
  if (cfg$assignment == "passage") {
    item_slots <- arm_slots(trials$condition[match(seq_len(n_item), trials$item_id)])
    passage_slots <- arm_slots(trials$condition[match(seq_len(cfg$passages), passage_id)])
  } else {
    item_slots <- seq_len(n_item)
    passage_slots <- seq_len(cfg$passages)
  }
  item <- item_slots[trials$item_id]
  passage <- passage_slots[passage_id]
  reader <- trials$subj_id
  reader_item <- reader + (item - 1L) * n_subj
  reader_passage <- reader + (passage - 1L) * n_subj
  mu_trial <- mu[as.character(trials$condition)]
  n_trial <- nrow(trials)
  stopifnot(!anyNA(mu_trial), length(unique(reader_item)) == n_trial)

  rows <- lapply(seq_len(n_rep), function(i) {
    sim <- glmmplanr:::with_rng_stream(streams[[cfg$condition]][[i]], {
      a <- stats::rnorm(n_subj, sd = cfg$subj_sd)
      d <- cfg$passage_sd * stats::rnorm(n_item)
      b <- stats::rnorm(n_item, sd = cfg$item_sd)
      u <- cfg$interaction_sd * stats::rnorm(n_trial)
      eps <- stats::runif(n_trial)
      # With one item per passage the local effect is independent for every answer.
      eta <- mu_trial + a[reader] + d[passage] + b[item] + u[reader_passage]
      as.integer(eps[reader_item] < stats::plogis(eta))
    })
    trials$y <- sim
    fit <- glmmplanr:::fit_single_dataset(trials, spec)
    valid <- glmmplanr:::valid_target_fit(fit$est, fit$se, isTRUE(fit$converged), fit$status,
                                          fit$prob_est, fit$prob_se, spec$scale)
    est <- if (valid) 100 * fit$prob_est else NA_real_
    se <- if (valid) 100 * fit$prob_se else NA_real_
    lower <- est - z_crit * se
    upper <- est + z_crit * se
    data.frame(replication = i, valid = valid, estimate_pp = est, se_pp = se,
               observed_a = mean(sim[is_a]), observed_b = mean(sim[!is_a]),
               covered = if (valid) lower <= truth && truth <= upper else NA,
               detected = valid && (lower > 0 || upper < 0), status = fit$status,
               diagnostic = fit$error_reason)
  })
  result <- do.call(rbind, rows)
  cat(sprintf("Cell %2d %2d x %2d %-19s %-7s valid %d/%d\n", cfg$cell, cfg$passages,
              cfg$questions, cfg$profile, cfg$condition, sum(result$valid), n_rep))
  flush.console()
  list(fits = result, logit_a = unname(mu['A']), logit_b = unname(mu['B']))
}

exact_ci <- function(k, n) {
  if (n == 0L) return(c(NA_real_, NA_real_))
  stats::binom.test(k, n)$conf.int[1:2]
}

summarize_case <- function(fits, cfg) {
  n <- nrow(fits)
  valid <- fits$valid
  n_valid <- sum(valid)
  covered <- sum(fits$covered, na.rm = TRUE)
  detected <- sum(fits$detected)
  cov_ci <- exact_ci(covered, n_valid)
  det_ci <- exact_ci(detected, n)
  coverage_status <- if (n_valid / n < 0.90 || (!is.na(cov_ci[2L]) && cov_ci[2L] < 0.90)) {
    "fail"
  } else if (cov_ci[1L] >= 0.90) "pass" else "indeterminate"
  detection_status <- if (cfg$condition == "null") {
    if (det_ci[2L] <= 0.075) "pass" else if (det_ci[1L] > 0.075) "fail" else "indeterminate"
  } else {
    if (det_ci[1L] >= 0.80) "pass" else if (det_ci[2L] < 0.80) "fail" else "indeterminate"
  }
  data.frame(attempted = n, valid = n_valid, covered = covered, detected = detected,
             coverage = if (n_valid) covered / n_valid else NA_real_,
             coverage_lo = cov_ci[1L], coverage_hi = cov_ci[2L], coverage_status,
             detection = detected / n, detection_lo = det_ci[1L], detection_hi = det_ci[2L],
             detection_status, bias_pp = mean(fits$estimate_pp[valid]) - 100 * (cfg$p_b - cfg$p_a),
             sd_pp = stats::sd(fits$estimate_pp[valid]), mean_se_pp = mean(fits$se_pp[valid]),
             mean_width_pp = mean(2 * z_crit * fits$se_pp[valid]),
             p_a_observed = mean(fits$observed_a), p_b_observed = mean(fits$observed_b))
}

cat(sprintf("Master seed %d, %d distinct cells x %d replications, %d workers\n",
            master_seed, nrow(cases), n_rep, workers))
cat("Package:", normalizePath(find.package("glmmplanr")), "\n")
results <- parallel::mclapply(seq_len(nrow(cases)), run_case, mc.cores = workers,
                              mc.set.seed = FALSE)
stopifnot(all(vapply(results, function(r) is.list(r) && nrow(r$fits) == n_rep, logical(1))))
fits <- do.call(rbind, lapply(seq_along(results), function(k) {
  data.frame(cell = k, results[[k]]$fits)
}))
stopifnot(!anyDuplicated(fits[c("cell", "replication")]), nrow(fits) == nrow(cases) * n_rep)
metrics <- do.call(rbind, lapply(seq_along(results), function(k) {
  cbind(logit_a = results[[k]]$logit_a, logit_b = results[[k]]$logit_b,
        summarize_case(results[[k]]$fits, cases[k, ]))
}))
summary <- cbind(cases, metrics)
scenarios$master_seed <- master_seed
summary$master_seed <- master_seed
utils::write.csv(scenarios, file.path(out, "testlet_allocation_scenarios.csv"), row.names = FALSE)
utils::write.csv(summary, file.path(out, "testlet_allocation_summary.csv"), row.names = FALSE)
connection <- gzfile(file.path(out, "testlet_allocation_fits.csv.gz"), "wt")
utils::write.csv(fits, connection, row.names = FALSE)
close(connection)
cat(sprintf("Completed %d simulated datasets.\n", nrow(fits)))
