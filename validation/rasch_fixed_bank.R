# Run with glmmplanr, lme4 and TAM installed.
library(glmmplanr)
library(lme4)
if (!requireNamespace("TAM", quietly = TRUE)) stop("Install TAM to run the fixed-bank comparison")

central <- seq(-1, 1, length.out = 14)
banks <- list(moderate = c(-1.5, central, 1.5), wide = c(-3.5, central, 3.5))
n_item <- length(banks$wide)
item_names <- sprintf("i%02d", seq_len(n_item))
n_rep <- 64L
limit <- 0.50

compare_bank <- function(trials, responses, difficulty) {
  trials$y <- responses[cbind(trials$subj_id, trials$item_id)]
  scores <- colSums(responses)
  observed_both <- scores > 0 & scores < nrow(responses)
  stopifnot(identical(colnames(responses), item_names), length(difficulty) == n_item)
  estimates <- matrix(NA_real_, 2L, n_item, dimnames = list(c("random_item", "fixed_item"), item_names))
  diagnostic <- rep(NA_character_, 2L)
  singular <- NA
  mixed <- tryCatch({
    fit <- glmer(y ~ 1 + (1 | subj) + (1 | item), trials, binomial(),
                 control = glmerControl(optimizer = "bobyqa"))
    conv <- fit@optinfo$conv
    messages <- conv$lme4$messages
    if (any(c(conv$opt, conv$lme4$code) != 0) ||
        any(!grepl("boundary (singular) fit", messages, fixed = TRUE))) {
      stop("glmer convergence diagnostic: ", paste(c(conv$opt, messages), collapse = " "))
    }
    fit
  }, warning = identity, error = identity)
  if (inherits(mixed, "condition")) {
    diagnostic[1] <- conditionMessage(mixed)
  } else {
    singular <- isSingular(mixed)
    modes <- ranef(mixed, condVar = FALSE)$item
    estimates[1, ] <- -unname(fixef(mixed)[1]) - modes[match(seq_len(n_item), rownames(modes)), 1L]
  }
  fixed <- tryCatch({
    if (!any(observed_both)) stop("No item has both response categories")
    fit <- TAM::tam.mml(responses[, observed_both, drop = FALSE],
                        constraint = "cases", verbose = FALSE)
    if (fit$iter >= fit$control$maxiter) stop("TAM reached its iteration limit")
    fit
  }, warning = identity, error = identity)
  if (inherits(fixed, "condition")) {
    diagnostic[2] <- conditionMessage(fixed)
  } else {
    stopifnot(setequal(rownames(fixed$xsi), item_names[observed_both]))
    estimates[2, ] <- fixed$xsi$xsi[match(item_names, rownames(fixed$xsi))]
    if (!all(observed_both)) diagnostic[2] <- "Constant item has no finite fixed difficulty"
  }
  errors <- sweep(estimates, 2L, difficulty, "-")
  complete <- rowSums(is.finite(errors)) == n_item
  diagnostic[!complete & is.na(diagnostic)] <- "Missing or non-finite difficulty"
  data.frame(method = rownames(errors), complete = complete,
             both_complete = all(complete), hard_zero = scores[n_item] == 0,
             hard_error = errors[, n_item], bank_rmse = sqrt(rowMeans(errors^2)),
             worst_error = apply(abs(errors), 1, max),
             singular = c(singular, NA), diagnostic = diagnostic, row.names = NULL)
}

design <- design_crossed(list())
constant_trials <- generate_design_schedule(design, 30L, n_item)$trials
for (response in 0:1) {
  constant <- matrix(response, 30L, n_item, dimnames = list(NULL, item_names))
  check <- compare_bank(constant_trials, constant, banks$wide)
  stopifnot(nrow(check) == 2L, !any(check$complete), all(!is.na(check$diagnostic)))
}
all_runs <- lapply(c(30L, 100L), function(n_subj) {
  trials <- generate_design_schedule(design, n_subj, n_item)$trials
  do.call(rbind, lapply(seq_len(n_rep), function(replication) {
    set.seed(811000L + 100L * n_subj + replication)
    ability <- rnorm(n_subj, sd = 0.8)
    uniform <- matrix(runif(n_subj * n_item), n_subj, dimnames = list(NULL, item_names))
    # Common abilities and uniform draws keep the central 14 response columns identical.
    do.call(rbind, lapply(names(banks), function(bank) {
      probability <- plogis(outer(ability, banks[[bank]], "-"))
      responses <- 1L * (uniform < probability)
      metrics <- compare_bank(trials, responses, banks[[bank]])
      data.frame(n_subj, bank, replication, metrics)
    }))
  }))
})
runs <- do.call(rbind, all_runs)
by_design <- split(runs, list(runs$n_subj, runs$bank, runs$method), drop = TRUE)
results <- do.call(rbind, lapply(by_design, function(x) {
  stopifnot(nrow(x) == n_rep, !anyDuplicated(x$replication))
  bank_pass <- x$complete & x$bank_rmse <= limit
  every_pass <- x$complete & x$worst_error <= limit
  stopifnot(all(!every_pass | bank_pass))
  paired <- x$hard_error[x$both_complete]
  hard_zero_errors <- x$hard_error[x$hard_zero & is.finite(x$hard_error)]
  interval <- glmmplanr:::wilson_mc_ci(sum(every_pass), n_rep)
  data.frame(x[1, c("n_subj", "bank", "method")], finite_bank = sum(x$complete),
             hard_zero = sum(x$hard_zero), bank_pass = sum(bank_pass),
             every_item_pass = sum(every_pass), mc_lo = interval[1], mc_hi = interval[2],
             paired_n = length(paired),
             paired_hard_bias = mean(paired), paired_hard_mcse = sd(paired) / sqrt(length(paired)),
             hard_zero_evaluable = length(hard_zero_errors),
             hard_bias_when_zero = if (length(hard_zero_errors)) mean(hard_zero_errors) else NA_real_,
             singular = if (all(is.na(x$singular))) NA_integer_ else sum(x$singular))
}))
print(results[order(results$bank, results$n_subj, results$method), ], row.names = FALSE, digits = 4)
if (any(!is.na(runs$diagnostic))) {
  print(subset(runs, !is.na(diagnostic), select = c(n_subj, bank, replication, method, diagnostic)),
        row.names = FALSE)
}
wide_100 <- subset(runs, bank == "wide" & n_subj == 100L)
paired_pass <- with(wide_100, tapply(complete & worst_error <= limit, list(replication, method), identity))
print(table(random_pass = paired_pass[, "random_item"], fixed_pass = paired_pass[, "fixed_item"]))

# Exercise the public planner with items drawn anew in each replication.
random_dgp <- calibrate_dgp(design, 0.65, subj_sd = 0.8, item_sd = 0.6)
plan <- plan_rasch_precision(random_dgp, data.frame(n_subj = 100L, n_item = 16L),
                             c(difficulty = limit, difficulty_max_abs_error = limit),
                             nsim = 32L, seed = 890123L, engine = "lme4")
results <- as.data.frame(plan)
stopifnot(all(plan$replications$difficulty <= plan$replications$difficulty_max_abs_error),
          with(results, n_met[target == "difficulty_max_abs_error"] == n_met[target == "joint"]))
print(results[c("target", "n_evaluable", "n_met", "mc_lo", "mc_hi", "decision")],
      row.names = FALSE)
print(vapply(c("glmmplanr", "lme4", "TAM"), function(p) as.character(packageVersion(p)), ""))
