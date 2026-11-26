# Numerical audits and decision tables for the empirical-bank/replication examples.
source("validation/planning_comparison_helpers.R")
x <- read.csv("validation/l2_targeting_pilot/summary.csv")
trials <- read.csv("validation/l2_targeting_pilot/fits.csv.gz")
bank <- read.csv("validation/l2_empirical_bank/items.csv")
stopifnot(nrow(subset(x, model == "Rasch" & engine == "native")) == 240L,
  sum(subset(x, model == "Rasch" & engine == "native")$nsim) == 15360L,
  sum(subset(x, model == "Rasch" & engine == "mirt")$nsim) == 768L,
  sum(subset(x, model == "2PL")$nsim) == 768L,
  !anyDuplicated(trials[c("case", "replicate", "model", "engine")]))
paired <- read.csv("validation/l2_targeting_pilot/paired_rasch.csv")
stopifnot(nrow(paired) == 768L, all(paired$data_hash_native == paired$data_hash_mirt))
full <- subset(x, model == "Rasch" & engine == "native" & layout == "full50")
q <- glmmplanr:::gh_quad_1d(61L)
stopifnot(abs(sum(q$weights) - 1) < 1e-12, abs(sum(q$weights * q$nodes^2) - 1) < 1e-12)
analytic <- do.call(rbind, lapply(seq_len(nrow(full)), function(k) {
  z <- full[k, ]; b <- subset(bank, model == "Rasch" & format == z$format)$b
  p <- plogis(outer(z$person_mean + z$person_sd * q$nodes, b, "-"))
  pbar <- drop(crossprod(q$weights, p))
  information <- z$n_subj * drop(crossprod(q$weights, p * (1 - p)))
  runs <- subset(trials, case == z$case & model == "Rasch" & engine == "native")
  response_mcse <- sd(runs$mean_response) / sqrt(nrow(runs))
  data.frame(z, expected_response = mean(pbar), response_mcse,
    response_z = (z$mean_response - mean(pbar)) / response_mcse,
    expected_extreme_items = sum(exp(z$n_subj * log1p(-pbar)) + pbar^z$n_subj),
    extreme_items_mcse = sd(runs$n_extreme_items) / sqrt(nrow(runs)),
    min_known_ability_information = min(information),
    max_known_ability_se_bound = max(1 / sqrt(information)))
}))
write.csv(analytic, "validation/l2_targeting_pilot/full_form_targeting.csv", row.names = FALSE)
cat("Largest absolute response-mean discrepancy in MCSE units:", max(abs(analytic$response_z)), "\n")
stopifnot(all(is.finite(analytic$response_z)), max(abs(analytic$response_z)) < 4)
search <- read.csv("validation/l2_replication_full/search.csv")
confirmation <- read.csv("validation/l2_replication_full/confirmation.csv")
counts <- split(search$evaluations, search$job)
stopifnot(nrow(search) == 16L, all(vapply(counts, function(z) length(z) == 2L && z[1] == z[2], logical(1))),
  !anyDuplicated(confirmation[c("case", "condition")]))
planned <- subset(confirmation, condition == "planned")
checked <- merge(search, planned, by = c("profile", "n_subj", "n_item"), all.x = TRUE)
null <- subset(confirmation, condition == "null",
  c("profile", "n_subj", "n_item", "detection", "detection_lo", "detection_hi"))
names(null)[4:6] <- paste0("null_", names(null)[4:6])
checked <- merge(checked, null, by = c("profile", "n_subj", "n_item"), all.x = TRUE)
checked$null_interval_above_nominal <- checked$null_detection_lo > .05
checked$coverage_interval_below_nominal <- checked$coverage_hi < .95
stopifnot(nrow(checked) == nrow(search), all(is.na(checked$n_subj) | !is.na(checked$joint_decision)))
write.csv(checked, "validation/l2_replication_full/confirmed_search.csv", row.names = FALSE)
print(checked[, c("profile", "preference", "method", "n_subj", "n_item", "joint_attainment",
  "joint_lo", "joint_hi", "joint_decision")], row.names = FALSE)
