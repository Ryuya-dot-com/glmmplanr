# Large, paired fixed-bank Rasch diagnostics and anchored DIF benchmark.
library(glmmplanr)
args <- commandArgs(TRUE)
mode <- if (length(args)) match.arg(args[1L], c("pilot", "calibration", "dif")) else "pilot"
out <- if (length(args) > 1L) args[2L] else "validation"
workers <- if (length(args) > 2L) as.integer(args[3L]) else 4L
stopifnot(workers > 0L)
if (.Platform$OS.type == "windows") workers <- 1L
dir.create(out, showWarnings = FALSE, recursive = TRUE)
source("validation/rasch_fit_benchmark_helpers.R")
engines <- if (length(args) > 3L) {
  match.arg(args[4L], "mfrmr")
} else {
  c("glmmplanr", "TAM", "mfrmr")
}
stopifnot(all(vapply(engines, requireNamespace, logical(1L), quietly = TRUE)),
          as.character(packageVersion("glmmplanr")) == "0.1.0")
forms <- as.matrix(read.csv("validation/two_pl_fixed_bank_booklets_forms.csv"))
stopifnot(nrow(forms) == 30L, all(rowSums(forms) == 8L), all(colSums(forms) == 15L))

layouts <- data.frame(n = c(120L, 240L, 480L, 480L), items_per_person = c(16L, 16L, 16L, 8L))
make_case <- function(layout, profile, sd = .8, difficulty_shift = 0, focal_slope = 1,
                      group_shift = 0, dif = 0) {
  cbind(layouts[layout, ], profile, sd, difficulty_shift, focal_slope, group_shift, dif,
        responses = layouts$n[layout] * layouts$items_per_person[layout],
        responses_per_item = layouts$n[layout] * layouts$items_per_person[layout] / j_bank)
}
calibration <- rbind(make_case(1:4, "targeted"), make_case(1:4, "narrow", sd = .4),
                     make_case(1:4, "off_target", difficulty_shift = 1.4),
                     make_case(c(1L, 3L), "low_discrimination", focal_slope = .35))
dif_cases <- rbind(make_case(c(1L, 3L, 4L), "group_impact", group_shift = .6),
                   make_case(c(1L, 3L, 4L), "group_impact_and_dif", group_shift = .6, dif = .7),
                   make_case(3L, "no_group_impact"),
                   make_case(3L, "dif_only", dif = .7))
seeds <- c(pilot = 91660001L, calibration = 91661001L, dif = 91662001L)
n_rep <- switch(mode, pilot = 3L, calibration = 512L, dif = 384L)
n_native <- if (mode == "pilot") 3L else 32L
cases <- switch(mode, calibration = calibration, dif = dif_cases,
                pilot = rbind(calibration[c(2L, 4L, 6L, 10L, 13L), ], dif_cases[c(2L, 3L, 6L), ]))
cases$task <- ifelse(cases$profile %in% dif_cases$profile, "dif", "calibration")
cases$case <- seq_len(nrow(cases))
cases$n_rep <- n_rep
cases$seed <- seeds[[mode]]
cases$native_rep <- ifelse(cases$task == "calibration", n_native, 0L)
# Check the executable and design before touching saved scenarios or batches.
# Old batches without a manifest remain evidence, but cannot be resumed safely.
manifest_file <- file.path(out, paste0(mode, "_manifest.rds"))
inputs <- c("validation/rasch_fit_benchmark.R", "validation/rasch_fit_benchmark_helpers.R",
            "validation/two_pl_fixed_bank_booklets_forms.csv")
binary <- list.files(system.file("libs", package = "glmmplanr"),
                     pattern = "\\.(so|dll|dylib)$", recursive = TRUE, full.names = TRUE)
manifest <- list(cases = cases, engines = engines, R = R.version.string,
  versions = vapply(engines, function(x) as.character(packageVersion(x)), ""),
  script_hashes = tools::md5sum(inputs), binary_hash = unname(tools::md5sum(binary)),
  native_adapter = paste(deparse(body(glmmplanr:::fit_rasch_fixed_mml)), collapse = "\n"))
if (file.exists(manifest_file)) {
  if (!identical(readRDS(manifest_file), manifest)) stop("Resume configuration differs; use a new output directory")
} else {
  existing <- list.files(out, pattern = paste0("^", mode, "_.*(batches|fits\\.csv\\.gz)$"))
  if (length(existing)) stop("Existing results have no resume manifest; use a new output directory")
  saveRDS(manifest, manifest_file)
}
write.csv(cases, file.path(out, paste0(mode, "_scenarios.csv")), row.names = FALSE)
streams <- glmmplanr:::make_lecuyer_streams(seeds[[mode]], n_rep)
writeLines(c(capture.output(sessionInfo()), "",
             paste(engines, vapply(engines, function(x) as.character(packageVersion(x)), character(1L))),
             sprintf("mode=%s workers=%d reps=%d seed=%d", mode, workers, n_rep, seeds[[mode]])),
           file.path(out, paste0(mode, "_environment.txt")))

measure <- function(expr) {
  start <- proc.time()
  warnings <- character()
  value <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) e)
  elapsed <- proc.time() - start
  list(value = value, cpu = unname(sum(elapsed[c("user.self", "sys.self")])),
       warning = paste(unique(warnings), collapse = " | "),
       error = if (inherits(value, "error")) conditionMessage(value) else "")
}

calibration_trial <- function(case, rep, truth, long) {
  y <- truth$y
  pat <- rasch_patterns(y)
  do.call(rbind, lapply(engines, function(engine) {
    fit <- measure(fit_rasch(engine, y, long))
    valid <- !nzchar(fit$error)
    metrics <- if (valid) common_diagnostics(y, pat, fit$value, truth)
    if (!valid) {
      # The row is retained even when an engine cannot fit the scheduled data.
      metrics <- as.data.frame(as.list(setNames(rep(NA_real_, 18L),
        c("sigma", "b_bias", "b_rmse", "b_max", "bank_ok", "theta_rmse", "person_eap_rel",
          "person_eap_ratio", "item_rel", "item_separation", "mean_item_infit", "mean_item_outfit",
          "mean_person_infit", "mean_person_outfit", "focal_infit", "focal_outfit",
          "person_overfit_pct", "person_underfit_pct"))))
    }
    native_names <- paste0("native_", c("item_infit", "item_outfit", "focal_infit", "focal_outfit",
      "person_infit", "person_outfit", "person_rel", "person_separation", "item_rel", "item_separation"))
    native_values <- as.data.frame(as.list(setNames(rep(NA_real_, length(native_names)), native_names)))
    scheduled <- engine != "glmmplanr" && rep <= case$native_rep
    native_error <- native_warning <- ""
    native_cpu <- NA_real_
    if (scheduled && valid) {
      native <- measure(native_diagnostics(engine, fit$value$fit))
      native_error <- native$error
      native_warning <- native$warning
      native_cpu <- native$cpu
      if (!nzchar(native$error)) native_values <- native$value
    }
    ll <- if (valid) rasch_ll(y, pat, fit$value$b, fit$value$sigma) else NA_real_
    data.frame(engine, valid, fit_cpu_s = fit$cpu, error = fit$error, warning = fit$warning,
               native_scheduled = scheduled, native_cpu_s = native_cpu,
               native_error, native_warning, ll_common = ll,
               correct = mean(y, na.rm = TRUE), metrics, native_values)
  }))
}

dif_trial <- function(case, truth, long) {
  y <- truth$y
  g <- truth$g
  y_alt <- cbind(y, i17 = ifelse(g == 1L, y[, focal_dif], NA_integer_))
  y_alt[g == 1L, focal_dif] <- NA_integer_
  long_alt <- long_rasch(y_alt)
  pat_null <- rasch_patterns(y, g)
  pat_alt <- rasch_patterns(y_alt, g)
  do.call(rbind, lapply(engines, function(engine) {
    null <- measure(fit_rasch(engine, y, long, g))
    alt <- measure(fit_rasch(engine, y_alt, long_alt, g))
    null_ok <- !nzchar(null$error)
    alt_ok <- !nzchar(alt$error)
    valid <- null_ok && alt_ok
    diff <- ll_delta <- ll_null_gap <- ll_alt_gap <- estimate <- group <- sigma <- NA_real_
    if (valid) {
      f <- null$value
      a <- alt$value
      diff <- 2 * (a$ll - f$ll)
      common_null <- rasch_ll(y, pat_null, f$b, f$sigma, f$group)
      common_alt <- rasch_ll(y_alt, pat_alt, a$b, a$sigma, a$group)
      ll_delta <- 2 * (common_alt - common_null)
      ll_null_gap <- f$ll - common_null
      ll_alt_gap <- a$ll - common_alt
      estimate <- a$b[focal_dif] - a$b[j_bank + 1L]
      group <- a$group
      sigma <- a$sigma
    }
    # A negative likelihood difference of more than roundoff makes this test
    # unavailable. Small negative discrepancies are set to zero explicitly.
    test_valid <- valid && is.finite(diff) && diff >= -1e-4
    p <- if (test_valid) pchisq(max(0, diff), df = 1L, lower.tail = FALSE) else NA_real_
    data.frame(engine, null_ok, alt_ok, valid, test_valid, fit_cpu_s = null$cpu + alt$cpu,
               null_error = null$error, alt_error = alt$error, null_warning = null$warning,
               alt_warning = alt$warning, estimate, group, sigma, lr = diff, lr_common = ll_delta,
               ll_null_gap, ll_alt_gap, p, rejected = p < .05,
               group_n = min(table(g)), focal_n = min(table(g[!is.na(y[, focal_dif])])))
  }))
}

run_trial <- function(job) {
  case <- cases[job[1L], ]
  rep <- job[2L]
  truth <- rasch_data(case, streams[[rep]], forms)
  long <- long_rasch(truth$y)
  result <- if (case$task == "calibration") {
    calibration_trial(case, rep, truth, long)
  } else {
    dif_trial(case, truth, long)
  }
  data.frame(case = case$case, replicate = rep, result)
}

start <- Sys.time()
for (task in unique(cases$task)) {
  selected <- cases$case[cases$task == task]
  batches <- split(seq_len(n_rep), ceiling(seq_len(n_rep) / 32L))
  file <- file.path(out, paste0(mode, "_", task, "_fits.csv.gz"))
  block_dir <- file.path(out, paste0(mode, "_", task, "_batches"))
  dir.create(block_dir, showWarnings = FALSE)
  block_files <- file.path(block_dir, sprintf("batch-%02d.rds", seq_along(batches)))
  for (batch in seq_along(batches)) {
    if (file.exists(block_files[batch])) next
    jobs <- expand.grid(case = selected, rep = batches[[batch]])
    index <- split(jobs, rep_len(seq_len(min(workers, nrow(jobs))), nrow(jobs)))
    rows <- parallel::mclapply(index, function(part) {
      lapply(seq_len(nrow(part)), function(i) run_trial(as.integer(part[i, ])))
    }, mc.cores = workers, mc.set.seed = FALSE)
    children_ok <- vapply(rows, function(x) !inherits(x, "try-error") && !is.null(x), logical(1L))
    if (!all(children_ok)) stop(paste("Failed worker(s):", paste(rows[!children_ok], collapse = " ")))
    block <- do.call(rbind, unlist(rows, recursive = FALSE))
    block <- block[order(block$case, block$replicate, match(block$engine, engines)), ]
    stopifnot(nrow(block) == nrow(jobs) * length(engines),
              !anyDuplicated(block[c("case", "replicate", "engine")]))
    tmp <- paste0(block_files[batch], ".tmp")
    saveRDS(block, tmp)
    stopifnot(file.rename(tmp, block_files[batch]))
    cat(sprintf("%s reps %d/%d | %.1fs elapsed\n", task, max(batches[[batch]]), n_rep,
                as.numeric(difftime(Sys.time(), start, units = "secs"))))
  }
  con <- gzfile(file, "wt")
  for (batch in seq_along(batches)) {
    write.table(readRDS(block_files[batch]), con, sep = ",", row.names = FALSE,
                col.names = batch == 1L, na = "")
  }
  close(con)
}
