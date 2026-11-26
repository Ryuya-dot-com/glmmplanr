# Exploratory sensitivity after failed default logistic surrogates.
# Refit the already collected 600 outcomes: no new simulation budget is used.
source("validation/planning_comparison_helpers.R")
versions <- readRDS("validation/item_pool_search/manifest.rds")$packages
stopifnot(identical(versions, vapply(names(versions),
  function(p) as.character(packageVersion(p)), "")))
search <- read.csv("validation/item_pool_search/search.csv")
jobs <- subset(search, task == "rasch" & method == "mlpwr" &
  grepl("No good design found", warning))
stopifnot(nrow(jobs) == 10L)
results <- lapply(seq_len(nrow(jobs)), function(k) {
  job <- jobs[k, ]
  saved <- readRDS(sprintf("validation/item_pool_search/search-%03d.rds", job$job))
  obs <- saved$mlpwr_observations; grid <- saved$grid
  dat <- lapply(split(obs, obs$pos), function(z)
    list(x = c(i = grid$i[z$pos[1L]]), y = z$success))
  calls <- 0L
  simfun <- function(i) { calls <<- calls + 1L; stop("No remaining simulation budget") }
  costfun <- function(i) grid$cost[as.integer(i)]
  warnings <- character()
  set.seed(9470000L + 1000L * job$job + 202L)
  elapsed <- system.time(result <- withCallingHandlers(mlpwr::find.design(
    simfun, boundaries = list(i = c(1, 5)), power = .8, costfun = costfun,
    dat = dat, evaluations = 0L, surrogate = "gpr", silent = TRUE),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    }))[["elapsed"]]
  stopifnot(calls == 0L, result$evaluations_used == nrow(obs))
  pos <- as.integer(result$final$design$i)
  good <- is.finite(result$final$power) && result$final$power >= .8 &&
    !any(grepl("No good design found", warnings))
  data.frame(job = job$job, family = job$family, run = job$run,
    case = grid$case[pos], n_subj = grid$n_subj[pos], n_item = grid$n_item[pos],
    prior_evaluations = nrow(obs), additional_evaluations = calls,
    predicted_attainment = result$final$power, qualified = good,
    elapsed_s = elapsed, warning = paste(unique(warnings), collapse = " | "))
})
result <- do.call(rbind, results)
write.csv(result, "validation/item_pool_search/surrogate_sensitivity.csv", row.names = FALSE)
print(result[, c("job", "family", "n_subj", "predicted_attainment", "qualified")], row.names = FALSE)
