# Run from the project root. An optional, existing output directory receives the two CSVs.
source("validation/partial_presentation.R")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) <= 1L, !length(args) || dir.exists(args[1L]))

# Fixed before the rerun: 96 new outer datasets, 99 refits, no interim stopping.
# At 95% coverage and 96 valid intervals, the nominal Monte Carlo SE is 0.022.
ids <- 10000L + seq_len(96L)
workers <- if (.Platform$OS.type == "windows") 1L else 4L
reps <- parallel::mclapply(ids, function(i) run_rep(i, "planned", "full", 99L),
                           mc.cores = workers, mc.set.seed = FALSE)
stopifnot(all(vapply(reps, is.data.frame, logical(1))))
raw <- do.call(rbind, reps)
stopifnot(!anyDuplicated(raw[c("replication", "engine")]),
          all(table(factor(raw$replication, levels = ids),
                    factor(raw$engine, levels = c("bayes_laplace", "glmer"))) == 1L))
result <- summarize_presentation(raw)
# Exploratory validation tolerance, distinct from the nominal 95% confidence level.
result$status <- with(result, ifelse(valid < ceiling(0.9 * attempted), "insufficient_valid",
                            ifelse(coverage_lo >= 0.90, "pass",
                            ifelse(coverage_hi < 0.90, "below_floor", "indeterminate"))))
result$glmmplanr_version <- as.character(packageVersion("glmmplanr"))
if (length(args)) {
  write.csv(raw, file.path(args[1], "coverage_recheck_fits.csv"), row.names = FALSE)
  write.csv(result, file.path(args[1], "coverage_recheck_summary.csv"), row.names = FALSE)
}
print(result[c("engine", "attempted", "valid", "covered", "coverage_lo", "coverage_hi", "status")],
      row.names = FALSE)
