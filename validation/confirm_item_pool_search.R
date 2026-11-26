# Independent planned/null checks for the union of GLMM search proposals.
source("validation/planning_comparison_helpers.R")
args <- commandArgs(TRUE)
out <- if (length(args)) args[1L] else "validation/item_pool_search_confirmation"
workers <- if (length(args) > 1L) as.integer(args[2L]) else 2L
if (.Platform$OS.type == "windows") workers <- 1L
search <- read.csv("validation/item_pool_search/search.csv")
selected <- subset(search, task == "glmm" & !is.na(n_subj))
cases <- unique(selected[c("family", "n_subj", "n_item")])
cases <- cases[order(cases$family, cases$n_subj, cases$n_item), ]
cases$case <- seq_len(nrow(cases))
checkpoint_dir(out, c("validation/confirm_item_pool_search.R", "validation/item_pool_search/search.csv"),
  list(cases = cases, nsim = 512L, seed = 95710001L, target = .8), "glmmplanr")
write.csv(cases, file.path(out, "scenarios.csv"), row.names = FALSE)
run <- function(k) {
  x <- cases[k, ]; path <- file.path(out, sprintf("case-%02d.rds", k))
  if (file.exists(path)) return(readRDS(path)$summary)
  spec <- experiment_spec(x$family)
  timing <- system.time(check <- validate_design_point(spec, x$n_subj, x$n_item,
    nsim = 512L, seed = 95710001L + 1000L * k))[["elapsed"]]
  summary <- data.frame(x, elapsed_s = timing, as.data.frame(check))
  summary$joint_n_met <- round(summary$joint_attainment * summary$nsim)
  ci <- t(vapply(summary$joint_n_met, glmmplanr:::wilson_mc_ci, numeric(2L), n = 512L))
  summary$joint_lo <- ci[, 1]; summary$joint_hi <- ci[, 2]
  summary$joint_decision <- ifelse(summary$joint_lo >= .8, "pass",
                                   ifelse(summary$joint_hi < .8, "fail", "indeterminate"))
  atomic_rds(list(summary = summary, check = check), path)
  cat(sprintf("GLMM confirmation %d/%d: %.1fs\n", k, nrow(cases), timing))
  summary
}
result <- parallel::mclapply(seq_len(nrow(cases)), run, mc.cores = workers, mc.set.seed = FALSE)
stopifnot(all(vapply(result, is.data.frame, logical(1))))
write.csv(do.call(rbind, result), file.path(out, "summary.csv"), row.names = FALSE)
