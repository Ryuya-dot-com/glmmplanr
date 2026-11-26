# Check the 0.3.2 DFF API on one fixed replicate per archived DIF condition.
args <- commandArgs(TRUE)
ref <- if (length(args)) args[1L] else "validation/rasch_fit_benchmark"
out <- if (length(args) > 1L) args[2L] else ref
workers <- if (length(args) > 2L) as.integer(args[3L]) else 2L
stopifnot(workers > 0L, packageVersion("mfrmr") == "0.3.2")
if (.Platform$OS.type == "windows") workers <- 1L
source("validation/rasch_fit_benchmark_helpers.R")
cases <- read.csv(file.path(ref, "dif_scenarios.csv"))
forms <- as.matrix(read.csv("validation/two_pl_fixed_bank_booklets_forms.csv"))
replicate <- 4L
stream <- glmmplanr:::make_lecuyer_streams(unique(cases$seed), replicate)[[replicate]]

check <- function(i) {
  case <- cases[i, ]
  data <- rasch_data(case, stream, forms)
  long <- long_rasch(data$y)
  long$g <- data$g[long$id]
  f <- fit_rasch("mfrmr", data$y, long, data$g)
  start <- proc.time()
  ans <- mfrmr::analyze_dff(f$fit, facet = "item", group = "g", focal = "1",
                            method = "mml_lrt", p_adjust = "holm", threads = 1L)
  elapsed <- proc.time() - start
  tab <- ans$table
  row <- tab[tab$Level == items[focal_dif], ]
  stopifnot(nrow(row) == 1L, nrow(tab) == j_bank, row$ReferenceGroup == "0", row$FocalGroup == "1",
            row$ObsRef == case$responses_per_item / 2L, row$ObsFocal == row$ObsRef)
  data.frame(case = case$case, replicate = replicate, focal_n = row$ObsFocal,
             estimate = -row$Contrast, lr = row$Statistic, p = row$pValue,
             p_holm = row$pAdj, status = row$InferenceStatus, reason = row$InferenceReason,
             tests_valid = sum(is.finite(tab$pValue)),
             cpu_s = unname(sum(elapsed[c("user.self", "sys.self")])),
             elapsed_s = unname(elapsed[["elapsed"]]))
}

rows <- parallel::mclapply(seq_len(nrow(cases)), check, mc.cores = workers, mc.set.seed = FALSE)
stopifnot(all(vapply(rows, is.data.frame, logical(1L))))
write.csv(do.call(rbind, rows), file.path(out, "dif_public_api.csv"), row.names = FALSE)
