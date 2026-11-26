# Audit pairing, actual budgets, and independent confirmation before reporting.
source("validation/planning_comparison_helpers.R")
paired_summary <- function(directory) {
  trials <- read.csv(file.path(directory, "fits.csv.gz"))
  cases <- read.csv(file.path(directory, "scenarios.csv"))
  a <- subset(trials, engine == "native")
  b <- subset(trials, engine == "irtsim")
  stopifnot(!anyDuplicated(a[c("case", "replicate")]), !anyDuplicated(b[c("case", "replicate")]))
  paired <- merge(a, b, by = c("case", "replicate"), suffixes = c("_native", "_irtsim"))
  stopifnot(nrow(paired) == nrow(a), nrow(paired) == nrow(b),
            all(paired$data_hash_native == paired$data_hash_irtsim))
  summary <- do.call(rbind, lapply(split(paired, paired$case), function(z) {
    n <- nrow(z)
    native_only <- sum(z$success_native & !z$success_irtsim)
    irtsim_only <- sum(!z$success_native & z$success_irtsim)
    # Conservative pointwise 95% interval for the paired difference: two
    # exact 97.5% binomial intervals for the discordant probabilities.
    a <- binom.test(native_only, n, conf.level = .975)$conf.int
    b <- binom.test(irtsim_only, n, conf.level = .975)$conf.int
    common <- z$valid_native & z$valid_irtsim
    delta <- z$rmse_native[common] - z$rmse_irtsim[common]
    data.frame(cases[cases$case == z$case[1], ], nsim = n,
      n_native = sum(z$valid_native), n_irtsim = sum(z$valid_irtsim), n_common = sum(common),
      native_only, irtsim_only, difference = (native_only - irtsim_only) / n,
      difference_lo = a[1] - b[2], difference_hi = a[2] - b[1],
      mean_rmse_difference = if (length(delta)) mean(delta) else NA_real_,
      rmse_difference_mcse = if (length(delta) > 1L) sd(delta) / sqrt(length(delta)) else NA_real_,
      max_item_difference = if (any(common)) max(z$b_max_difference_native[common]) else NA_real_)
  }))
  write.csv(summary, file.path(directory, "paired.csv"), row.names = FALSE)
  summary
}
screen_pairs <- paired_summary("validation/item_pool_irtsim_screen")
confirm_pairs <- paired_summary("validation/item_pool_irtsim_confirm")
search <- read.csv("validation/item_pool_search/search.csv")
search$predicted_attainment <- vapply(seq_len(nrow(search)), function(i) {
  if (search$method[i] != "mlpwr") return(NA_real_)
  x <- readRDS(sprintf("validation/item_pool_search/search-%03d.rds", search$job[i]))
  if (is.null(x$mlpwr_final)) NA_real_ else as.numeric(x$mlpwr_final$power)
}, numeric(1))
search$proposal_status <- ifelse(nzchar(search$error), "search_error",
  ifelse(is.na(search$n_subj), "no_candidate",
    ifelse(search$method == "mlpwr" &
      (grepl("No good design found", search$warning) |
       is.na(search$predicted_attainment) | search$predicted_attainment < .80),
      "unqualified_return", "candidate")))
counts <- split(search$actual_evaluations, search$job)
stopifnot(all(lengths(counts) == 2L), all(vapply(counts, function(x) x[1] == x[2], logical(1))),
          !anyDuplicated(search[c("job", "method")]))
rc <- read.csv("validation/item_pool_irtsim_confirm/summary.csv")
gc <- subset(read.csv("validation/item_pool_search_confirmation/summary.csv"), condition == "planned")
sensitivity <- read.csv("validation/item_pool_search/surrogate_sensitivity.csv")
sensitivity <- merge(sensitivity, subset(rc, engine == "native",
  c("case", "nsim", "n_valid", "n_met", "attainment", "mc_lo", "mc_hi", "decision")),
  by = "case", all.x = TRUE)
stopifnot(nrow(sensitivity) == 10L, !anyNA(sensitivity$nsim),
  all(sensitivity$prior_evaluations == 600L), all(sensitivity$additional_evaluations == 0L))
write.csv(sensitivity, "validation/item_pool_search/confirmed_surrogate_sensitivity.csv", row.names = FALSE)
checked <- lapply(seq_len(nrow(search)), function(i) {
  s <- search[i, ]
  if (is.na(s$n_subj)) return(data.frame(s, confirm_n = 0L, confirm_valid = 0L,
    confirm_attainment = NA_real_, confirm_lo = NA_real_, confirm_hi = NA_real_, confirm_decision = "no_proposal"))
  if (s$task == "rasch") {
    x <- subset(rc, case == s$case & engine == "native")
    stopifnot(nrow(x) == 1L)
    data.frame(s, confirm_n = x$nsim, confirm_valid = x$n_valid,
      confirm_attainment = x$attainment, confirm_lo = x$mc_lo, confirm_hi = x$mc_hi, confirm_decision = x$decision)
  } else {
    x <- subset(gc, family == s$family & n_subj == s$n_subj & n_item == s$n_item)
    stopifnot(nrow(x) == 1L)
    data.frame(s, confirm_n = x$nsim, confirm_valid = x$n_valid,
      confirm_attainment = x$joint_attainment, confirm_lo = x$joint_lo,
      confirm_hi = x$joint_hi, confirm_decision = x$joint_decision)
  }
})
checked <- do.call(rbind, checked)
# Preserve checks of fallback coordinates but do not label them recommendations.
checked$returned_point_decision <- checked$confirm_decision
checked$confirm_decision[checked$proposal_status != "candidate"] <- "no_proposal"
write.csv(checked, "validation/item_pool_search/confirmed_search.csv", row.names = FALSE)
aggregate <- do.call(rbind, lapply(split(checked, interaction(checked$task, checked$family,
                                                             checked$preference, checked$method, drop = TRUE)), function(x)
  data.frame(x[1, c("task", "family", "preference", "method")], runs = nrow(x),
    proposals = sum(x$proposal_status == "candidate"),
    search_errors = sum(x$proposal_status == "search_error"),
    unqualified_returns = sum(x$proposal_status == "unqualified_return"),
    confirmed = sum(x$confirm_decision == "pass"),
    indeterminate = sum(x$confirm_decision == "indeterminate"), failed = sum(x$confirm_decision == "fail"),
    evaluations = sum(x$actual_evaluations), elapsed_s = sum(x$elapsed_s),
    returned_designs = paste(paste(x$n_subj, x$n_item, sep = "x"), collapse = ";"),
    qualified_designs = paste(paste(x$n_subj[x$proposal_status == "candidate"],
      x$n_item[x$proposal_status == "candidate"], sep = "x"), collapse = ";"),
    mean_cost = if (any(x$proposal_status == "candidate"))
      mean(x$cost[x$proposal_status == "candidate"]) else NA_real_)))
write.csv(aggregate, "validation/item_pool_search/summary.csv", row.names = FALSE)
proposals <- read.csv("validation/item_pool_irtsim_screen/selected.csv")
selected <- do.call(rbind, lapply(seq_len(nrow(proposals)), function(i) {
  s <- proposals[i, ]
  if (is.na(s$case)) return(data.frame(family = s$family, engine = s$engine,
    proposed_n = NA_integer_, n_met = NA_integer_, nsim = 0L, attainment = NA_real_,
    mc_lo = NA_real_, mc_hi = NA_real_, decision = "no_proposal"))
  x <- subset(rc, case == s$case & engine == s$engine)
  stopifnot(nrow(x) == 1L)
  data.frame(family = s$family, engine = s$engine, proposed_n = s$n_subj,
    x[c("n_met", "nsim", "attainment", "mc_lo", "mc_hi", "decision")])
}))
write.csv(selected, "validation/item_pool_irtsim_confirm/confirmed_proposals.csv", row.names = FALSE)
cat("Paired IRT responses and equal search budgets verified.\n")
print(aggregate[, c("task", "family", "preference", "method", "runs", "confirmed", "indeterminate", "failed")], row.names = FALSE)
