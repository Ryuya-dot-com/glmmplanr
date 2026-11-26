# Independent count/scale audits; never pool pilot or search draws into confirmation.
source("validation/planning_comparison_helpers.R")
args <- commandArgs(TRUE)
stage <- if (length(args)) match.arg(args[1L], c("all", "targeting", "replication")) else "all"
close_enough <- function(a, b, tolerance = 1e-12) {
  stopifnot(length(a) == length(b), identical(is.na(a), is.na(b)),
    all(abs(a[!is.na(a)] - b[!is.na(b)]) <= tolerance))
}
if (stage %in% c("all", "targeting")) {
  out <- "validation/l2_targeting_large"
  s <- read.csv(file.path(out, "summary.csv"))
  z <- read.csv(file.path(out, "fits.csv.gz"))
  stopifnot(nrow(s) == 440L, nrow(z) == 343040L,
    !anyDuplicated(z[c("case", "replicate", "model", "engine")]),
    sum(z$engine == "native") == 307200L,
    sum(z$model == "Rasch" & z$engine == "mirt") == 17920L,
    sum(z$model == "2PL") == 17920L)
  groups <- split(z, interaction(z$case, z$model, z$engine, drop = TRUE))
  recovery <- do.call(rbind, lapply(groups, function(v) {
    joint <- v$valid & v$rmse <= .35 & v$max_error <= .60
    parameter <- joint & (v$model == "Rasch" | (v$a_rmse <= .25 & v$a_max_error <= .50))
    curve <- v$valid & v$probability_rmse <= .05
    stopifnot(identical(joint, v$success), identical(parameter, v$parameter_success),
      identical(curve, v$probability_success))
    key <- s$case == v$case[1L] & s$model == v$model[1L] & s$engine == v$engine[1L]
    a <- s[key, ]; b <- mc_summary(joint, v$valid)
    stopifnot(nrow(a) == 1L, a$nsim == nrow(v), a$n_valid == sum(v$valid), a$n_met == sum(joint),
      a$decision == b$decision)
    close_enough(c(a$attainment, a$mc_lo, a$mc_hi), c(b$attainment, b$mc_lo, b$mc_hi))
    metrics <- c("difficulty", "parameters", "probability", if (v$model[1L] == "2PL")
      c("difficulty_rmse_only_diagnostic", "discrimination_rmse_only_diagnostic"))
    do.call(rbind, lapply(metrics, function(metric) {
      met <- switch(metric, difficulty = joint, parameters = parameter, probability = curve,
        difficulty_rmse_only_diagnostic = v$valid & v$rmse <= .35,
        discrimination_rmse_only_diagnostic = v$valid & v$a_rmse <= .25)
      data.frame(a[c("case", "format", "layout", "profile", "n_subj", "model", "engine")],
        metric, mc_summary(met, v$valid))
    }))
  }))
  write.csv(recovery, file.path(out, "recovery_metrics.csv"), row.names = FALSE)
  native <- subset(s, engine == "native")
  minimum <- do.call(rbind, lapply(split(native,
    interaction(native$format, native$layout, native$profile, drop = TRUE)), function(v) {
      v <- v[order(v$n_subj), ]
      passing <- which(v$decision == "pass")
      chosen <- if (length(passing)) passing[1L] else NA_integer_
      data.frame(v[1L, c("format", "layout", "profile")],
        smallest_tested_passing_n = if (!is.na(chosen)) v$n_subj[chosen] else NA_integer_,
        attainment = if (!is.na(chosen)) v$attainment[chosen] else NA_real_,
        lower = if (!is.na(chosen)) v$mc_lo[chosen] else NA_real_,
        upper = if (!is.na(chosen)) v$mc_hi[chosen] else NA_real_,
        smaller_tested_indeterminate = if (!is.na(chosen))
          any(v$n_subj < v$n_subj[chosen] & v$decision == "indeterminate") else NA,
        largest_tested_n = max(v$n_subj),
        largest_tested_attainment = tail(v$attainment, 1L))
    }))
  write.csv(minimum, file.path(out, "smallest_confirmed.csv"), row.names = FALSE)
  pairs <- read.csv(file.path(out, "paired_rasch.csv"))
  stopifnot(nrow(pairs) == 17920L, all(pairs$data_hash_native == pairs$data_hash_mirt))
  paired_summary <- do.call(rbind, lapply(split(pairs, pairs$case), function(v) {
    agree <- v$success_native == v$success_mirt
    native_check <- mc_summary(v$success_native, v$valid_native)
    mirt_check <- mc_summary(v$success_mirt, v$valid_mirt)
    data.frame(case = v$case[1L], nsim = nrow(v),
      valid_disagreements = sum(v$valid_native != v$valid_mirt),
      native_only = sum(v$success_native & !v$success_mirt),
      mirt_only = sum(!v$success_native & v$success_mirt),
      agreement = mean(agree), agreement_lo = glmmplanr:::wilson_mc_ci(sum(agree), nrow(v))[1L],
      agreement_hi = glmmplanr:::wilson_mc_ci(sum(agree), nrow(v))[2L],
      native_met = native_check$n_met, mirt_met = mirt_check$n_met,
      native_decision = native_check$decision, mirt_decision = mirt_check$decision)
  }))
  write.csv(paired_summary, file.path(out, "paired_summary.csv"), row.names = FALSE)
  bank <- read.csv("validation/l2_empirical_bank/items.csv")
  full <- subset(s, model == "Rasch" & engine == "native" & layout == "full50")
  q <- glmmplanr:::gh_quad_1d(61L)
  stopifnot(abs(sum(q$weights) - 1) < 1e-12, abs(sum(q$weights * q$nodes^2) - 1) < 1e-12)
  analytic <- do.call(rbind, lapply(seq_len(nrow(full)), function(i) {
    a <- full[i, ]; b <- bank$b[bank$model == "Rasch" & bank$format == a$format]
    p <- plogis(outer(a$person_mean + a$person_sd * q$nodes, b, "-"))
    pbar <- drop(crossprod(q$weights, p))
    v <- subset(z, case == a$case & model == "Rasch" & engine == "native")
    data.frame(a, expected_response = mean(pbar), response_mcse = sd(v$mean_response)/sqrt(nrow(v)),
      expected_extreme_items = sum(pbar^a$n_subj + exp(a$n_subj * log1p(-pbar))),
      extreme_mcse = sd(v$n_extreme_items)/sqrt(nrow(v)))
  }))
  analytic$response_z <- with(analytic, (mean_response - expected_response)/response_mcse)
  write.csv(analytic, file.path(out, "full_form_targeting.csv"), row.names = FALSE)
  cat("Targeting: all counts, hashes and decisions verified. Largest analytic response |z|:",
    max(abs(analytic$response_z)), "\n")
  stopifnot(all(is.finite(analytic$response_z)), max(abs(analytic$response_z)) < 4)
  old <- read.csv("validation/l2_targeting_pilot/summary.csv")
  pilot <- merge(s, old, by = c("format", "layout", "n_subj", "profile", "model", "engine"),
    suffixes = c("_large", "_pilot"))
  pilot$difference <- pilot$attainment_large - pilot$attainment_pilot
  write.csv(pilot, file.path(out, "pilot_comparison.csv"), row.names = FALSE)
  bad <- subset(z, !valid)
  bad$error[is.na(bad$error) | !nzchar(bad$error)] <- "Incomplete bank (nonfinite item estimate)"
  errors <- aggregate(list(n = rep(1L, nrow(bad))), bad[c("model", "engine", "error")], sum)
  stopifnot(sum(errors$n) == nrow(bad))
  write.csv(errors, file.path(out, "failure_counts.csv"), row.names = FALSE)
  failures <- subset(bad, engine == "native" & error != "Incomplete bank (nonfinite item estimate)")
  write.csv(merge(failures, s[c("case", "model", "engine", "format", "layout", "profile", "n_subj")],
    by = c("case", "model", "engine")), file.path(out, "native_numerical_failures.csv"), row.names = FALSE)
}
if (stage %in% c("all", "replication")) {
  out <- "validation/l2_replication_large"
  summary <- read.csv(file.path(out, "confirmation.csv"))
  search <- read.csv(file.path(out, "search.csv"))
  stopifnot(nrow(summary) == 120L, all(summary$nsim == 2048L),
    sum(summary$nsim) == 245760L, !anyDuplicated(summary[c("case", "condition")]), nrow(search) == 80L)
  counts <- split(search$evaluations, search$job)
  stopifnot(all(vapply(counts, function(x) length(x) == 2L && x[1L] == x[2L], logical(1))))
  trials <- do.call(rbind, lapply(sort(unique(summary$case)), function(k) {
    v <- readRDS(file.path(out, sprintf("confirm-%02d.rds", k)))$trials
    v$case <- k; v
  }))
  stopifnot(nrow(trials) == 245760L, !anyDuplicated(trials[c("case", "condition", "replicate")]))
  summary$rms_se_pp <- summary$standardized_error_sd <- NA_real_
  for (i in seq_len(nrow(summary))) {
    a <- summary[i, ]; v <- subset(trials, case == a$case & condition == a$condition)
    d <- v$valid & abs(v$estimate / v$se) > qnorm(.975)
    w <- v$valid & 2 * qnorm(.975) * v$se <= 10
    covered <- v$valid & abs(v$estimate - v$truth) <= qnorm(.975) * v$se
    close_enough(c(mean(d), mean(w), mean(d & w), sum(covered)/sum(v$valid)),
      c(a$detection, a$width_attainment, a$joint_attainment, a$coverage))
    b <- mc_summary(d & w, v$valid)
    stopifnot(a$n_valid == sum(v$valid), b$decision == a$joint_decision)
    close_enough(c(a$joint_lo, a$joint_hi), c(b$mc_lo, b$mc_hi))
    close_enough(a$bias_pp, mean(v$estimate[v$valid] - v$truth[v$valid]))
    close_enough(a$mean_se_pp, mean(v$se[v$valid]))
    summary$rms_se_pp[i] <- sqrt(mean(v$se[v$valid]^2))
    summary$standardized_error_sd[i] <- sd((v$estimate[v$valid] - v$truth[v$valid]) / v$se[v$valid])
  }
  write.csv(summary, file.path(out, "diagnostic_summary.csv"), row.names = FALSE)
  con <- gzfile(file.path(out, "fits.csv.gz"), "wt")
  write.csv(trials, con, row.names = FALSE); close(con)
  planned <- subset(summary, condition == "planned")
  null <- subset(summary, condition == "null", c("profile", "n_subj", "n_item", "detection", "detection_lo", "detection_hi"))
  names(null)[4:6] <- paste0("null_", names(null)[4:6])
  decisions <- merge(planned, null, by = c("profile", "n_subj", "n_item"))
  decisions$null_interval_above_nominal <- decisions$null_detection_lo > .05
  decisions$coverage_interval_below_nominal <- decisions$coverage_hi < .95
  decisions$se_to_empirical_sd <- decisions$mean_se_pp / decisions$empirical_sd_pp
  decisions$rms_se_to_empirical_sd <- decisions$rms_se_pp / decisions$empirical_sd_pp
  write.csv(decisions, file.path(out, "design_decisions.csv"), row.names = FALSE)
  robustness <- do.call(rbind, lapply(split(decisions,
    interaction(decisions$n_subj, decisions$n_item, drop = TRUE)), function(v) {
      stopifnot(nrow(v) == 4L, length(unique(v$profile)) == 4L)
      data.frame(n_subj = v$n_subj[1L], n_item = v$n_item[1L],
        profiles_passing = sum(v$joint_decision == "pass"),
        minimum_joint_attainment = min(v$joint_attainment), minimum_joint_lower = min(v$joint_lo),
        largest_null_rejection = max(v$null_detection), lowest_planned_coverage = min(v$coverage),
        profiles_with_inference_flags = sum(v$null_interval_above_nominal | v$coverage_interval_below_nominal))
    }))
  write.csv(robustness, file.path(out, "scenario_robustness.csv"), row.names = FALSE)
  checked <- merge(search, decisions, by = c("profile", "n_subj", "n_item"), all.x = TRUE)
  stopifnot(nrow(checked) == 80L, all(is.na(checked$n_subj) | !is.na(checked$joint_decision)))
  write.csv(checked, file.path(out, "confirmed_search.csv"), row.names = FALSE)
  search_summary <- do.call(rbind, lapply(split(checked,
    interaction(checked$profile, checked$preference, checked$method, drop = TRUE)), function(v) {
      data.frame(v[1L, c("profile", "preference", "method")], searches = nrow(v),
        evaluations = sum(v$evaluations),
        confirmed = sum(v$joint_decision == "pass", na.rm = TRUE),
        indeterminate = sum(v$joint_decision == "indeterminate", na.rm = TRUE),
        failed = sum(v$joint_decision == "fail", na.rm = TRUE),
        no_candidate = sum(is.na(v$n_subj)),
        search_errors = sum(!is.na(v$error) & nzchar(v$error)),
        boundary_warnings = sum(!is.na(v$warning) & grepl("edge of the parameter space", v$warning)),
        proposals_with_inference_flags = sum(v$null_interval_above_nominal |
          v$coverage_interval_below_nominal, na.rm = TRUE))
    }))
  write.csv(search_summary, file.path(out, "search_summary.csv"), row.names = FALSE)
  print(with(checked, table(method, joint_decision, useNA = "ifany")))
  cat("Replication: all 245760 trial outcomes, counts, diagnostic summaries and matched budgets verified\n")
  old <- read.csv("validation/l2_replication_full/confirmation.csv")
  pilot <- merge(summary, old, by = c("profile", "n_subj", "n_item", "condition"), suffixes = c("_large", "_pilot"))
  write.csv(pilot, file.path(out, "pilot_comparison.csv"), row.names = FALSE)
}
