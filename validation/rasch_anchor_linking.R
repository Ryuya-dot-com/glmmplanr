# Rasch bank expansion at a fixed per-examinee response budget.
library(glmmplanr)
args <- commandArgs(TRUE)
out <- if (length(args)) args[1L] else "validation"
stage <- if (length(args) > 1L) match.arg(args[2L], c("pilot", "large")) else "pilot"
workers <- if (length(args) > 2L) as.integer(args[3L]) else 4L
stopifnot(dir.exists(out), workers > 0L)
if (.Platform$OS.type == "windows") workers <- 1L
n_rep <- if (stage == "pilot") 32L else 2048L
seed <- if (stage == "pilot") 90541001L else 90542001L
if (stage == "large") stopifnot(requireNamespace("TAM", quietly = TRUE))
streams <- glmmplanr:::make_lecuyer_streams(seed, n_rep)

n_grid <- c(120L, 240L, 480L, 960L, 1920L)
k_grid <- c(4L, 8L, 12L)
core <- expand.grid(n = n_grid, anchors = k_grid, band = c("broad", "hard"),
                    focal_mean = c(0, -.75), drift = "none", stringsAsFactors = FALSE)
core$new_per_person <- 16L - core$anchors
core$phase <- "fixed_workload"
drift <- expand.grid(n = c(480L, 960L, 1920L), anchors = k_grid, band = "broad", focal_mean = -.75,
                     drift = c("quarter_0.6", "half_0.8", "all_0.4"), stringsAsFactors = FALSE)
drift$new_per_person <- 16L - drift$anchors
drift$phase <- "drift"
extra <- expand.grid(n = c(480L, 960L), anchors = c(4L, 12L), band = c("broad", "hard"),
                     focal_mean = -.75, drift = "none", stringsAsFactors = FALSE)
extra$new_per_person <- 8L
extra$phase <- "fixed_new_exposure"
cases <- rbind(core, drift, extra)
cases$case <- seq_len(nrow(cases))
cases$per_person <- cases$anchors + cases$new_per_person
cases$per_new_item <- cases$n * cases$new_per_person / 12L
cases$master_seed <- seed
cases$n_rep <- n_rep
cases$tam <- stage == "large" & with(cases, phase == "fixed_workload" & focal_mean == -.75 &
  n %in% c(240L, 960L) & ((band == "broad" & anchors < 12L) | (band == "hard" & anchors == 12L)))
stopifnot(nrow(cases) == 95L, sum(cases$tam) == if (stage == "large") 6L else 0L)
ref <- c(seq(-2.2, 2.2, length.out = 12L), seq(1.3, 3.5, length.out = 12L))
new <- c(seq(-2, 2, length.out = 10L), -3, 3)
item_columns <- paste0("new_centered_", seq_along(new))
theta_grid <- seq(-2, 2, by = .25)
expected_score <- function(b) rowSums(plogis(outer(theta_grid, b, "-")))
true_score <- expected_score(new)

schedule <- function(cfg) {
  k <- cfg$anchors
  n <- cfg$n
  selected <- as.integer(round(seq(1, 12, length.out = k)))
  if (cfg$band == "hard") selected <- selected + 12L
  new_mask <- outer(seq_len(n) - 1L, seq_along(new) - 1L,
                    function(s, j) (j - s) %% length(new) < cfg$new_per_person)
  tr <- which(cbind(matrix(TRUE, n, k), new_mask), arr.ind = TRUE)
  stopifnot(all(tabulate(tr[, 1L], n) == cfg$per_person),
            identical(as.numeric(tabulate(tr[, 2L], k + length(new))),
                      c(rep(n, k), rep(cfg$per_new_item, length(new)))))
  shifts <- numeric(k)
  if (cfg$drift != "none") {
    sizes <- c("quarter_0.6" = .25, "half_0.8" = .50, "all_0.4" = 1)
    values <- c("quarter_0.6" = .60, "half_0.8" = .80, "all_0.4" = .40)
    changed <- as.integer(round(seq(1, k, length.out = k * sizes[[cfg$drift]])))
    shifts[changed] <- values[[cfg$drift]]
  }
  bank <- c(ref[selected] + shifts, new)
  global <- c(selected, length(ref) + seq_along(new))
  list(n = n, j = length(bank), k = k, sid = tr[, 1L], iid = tr[, 2L],
       index = tr[, 1L] + (global[tr[, 2L]] - 1L) * max(n_grid), selected = selected,
       b = bank, shifts = shifts, exposure = tabulate(tr[, 2L], length(bank)))
}
schedules <- lapply(seq_len(nrow(cases)), function(i) schedule(cases[i, ]))
cases$drift_mean <- vapply(schedules, function(s) mean(s$shifts), numeric(1L))
cases$anchor_min <- vapply(schedules, function(s) min(ref[s$selected]), numeric(1L))
cases$anchor_max <- vapply(schedules, function(s) max(ref[s$selected]), numeric(1L))
fit_native <- function(y, sid, iid, n, j, finite) {
  fit <- glmmplanr:::fit_rasch_fixed_mml(list(subj_id = sid, item_id = iid, y = y), n, j, finite)
  if (!is.na(fit$error)) return(rep(NA_real_, j))
  fit$b_est
}
fit_tam <- function(y, s, finite) {
  nm <- sprintf("i%02d", seq_len(s$j))
  resp <- matrix(NA_integer_, s$n, s$j, dimnames = list(NULL, nm))
  resp[cbind(s$sid, s$iid)] <- y
  usable <- resp[, finite, drop = FALSE]
  keep <- rowSums(!is.na(usable)) > 0L
  out <- rep(NA_real_, s$j)
  if (sum(keep) < 2L) return(out)
  fit <- tryCatch(TAM::tam.mml(usable[keep, , drop = FALSE], constraint = "cases", verbose = FALSE),
                  error = identity, warning = identity)
  if (inherits(fit, "condition") || fit$iter >= fit$control$maxiter) return(out)
  stopifnot(setequal(rownames(fit$xsi), nm[finite]))
  fit$xsi$xsi[match(nm, rownames(fit$xsi))]
}
link <- function(b, anchor_b, mu) {
  k <- length(anchor_b)
  anchors <- seq_len(k)
  delta <- anchor_b - b[anchors]
  n_anchor <- sum(is.finite(delta))
  shift <- if (n_anchor >= 2L) mean(delta[is.finite(delta)]) else NA_real_
  difficulty <- b[-anchors] + shift
  error <- difficulty - new
  complete <- all(is.finite(error))
  rmse <- sqrt(mean(error^2))
  worst <- max(abs(error))
  score_error <- if (complete) max(abs(expected_score(difficulty) - true_score)) else NA_real_
  hit <- complete && abs(shift - mu) <= .15 && rmse <= .35 && worst <= .60
  c(anchors_used = n_anchor, shift = shift, rmse = rmse, worst = worst,
    tcc_error = score_error, complete = complete, achieved = hit)
}
stopifnot(abs(link(c(ref[1:4], new) + .75, ref[1:4], -.75)[["shift"]] + .75) < 1e-12,
          link(c(ref[1:4], new) + .75, ref[1:4], -.75)[["achieved"]] == 1,
          link(c(ref[1:4], NA, new[-1]), ref[1:4], 0)[["achieved"]] == 0)

reference_n <- 500L
ref_sid <- rep(seq_len(reference_n), length(ref))
ref_iid <- rep(seq_along(ref), each = reference_n)
clock <- function() sum(proc.time()[1:2])
run_replication <- function(r) {
  sim <- glmmplanr:::with_rng_stream(streams[[r]], {
    ref_theta <- rnorm(reference_n)
    ref_y <- as.integer(runif(length(ref_sid)) < plogis(ref_theta[ref_sid] - ref[ref_iid]))
    list(ref_y = ref_y, z = rnorm(max(n_grid)), u = runif(max(n_grid) * (length(ref) + length(new))))
  })
  ref_success <- tabulate(ref_iid[sim$ref_y == 1L], length(ref))
  ref_est <- fit_native(sim$ref_y, ref_sid, ref_iid, reference_n, length(ref),
                        ref_success > 0L & ref_success < reference_n)
  rows <- vector("list", nrow(cases))
  ext <- list()
  for (i in seq_len(nrow(cases))) {
    cfg <- cases[i, ]
    s <- schedules[[i]]
    y <- as.integer(sim$u[s$index] < plogis(sim$z[s$sid] + cfg$focal_mean - s$b[s$iid]))
    success <- tabulate(s$iid[y == 1L], s$j)
    finite <- success > 0L & success < s$exposure
    timed <- cfg$tam && r <= 64L
    if (timed) start <- clock()
    b <- fit_native(y, s$sid, s$iid, s$n, s$j, finite)
    if (timed) native_time <- clock() - start
    oracle <- link(b, ref[s$selected], cfg$focal_mean)
    estimated <- link(b, ref_est[s$selected], cfg$focal_mean)
    rows[[i]] <- c(case = i, replication = r, native_fit = any(is.finite(b)),
                   anchors_constant = sum(!finite[seq_len(s$k)]), new_constant = sum(!finite[-seq_len(s$k)]),
                   unlinked_bias = mean(b[-seq_len(s$k)] - new),
                   setNames(oracle, paste0("oracle_", names(oracle))),
                   setNames(estimated, paste0("estimated_", names(estimated))),
                   setNames(b[-seq_len(s$k)], item_columns))
    if (timed) {
      start <- clock()
      tt <- fit_tam(y, s, finite)
      tam_time <- clock() - start
      tam <- link(tt, ref[s$selected], cfg$focal_mean)
      both <- oracle[["complete"]] == 1 && tam[["complete"]] == 1
      shift_difference <- oracle[["shift"]] - tam[["shift"]]
      max_item_difference <- if (both) max(abs((b - tt)[-seq_len(s$k)] + shift_difference)) else NA_real_
      ext[[length(ext) + 1L]] <- c(case = i, replication = r, native_cpu = native_time, tam_cpu = tam_time,
        native_complete = oracle[["complete"]], tam_complete = tam[["complete"]],
        native_hit = oracle[["achieved"]], tam_hit = tam[["achieved"]], shift_difference = shift_difference,
        max_item_difference = max_item_difference)
    }
  }
  list(native = do.call(rbind, rows), external = do.call(rbind, ext))
}

prefix <- file.path(out, paste0("rasch_anchor_", stage))
checkpoint <- file.path(out, "checkpoints")
dir.create(checkpoint, showWarnings = FALSE)
manifest <- file.path(checkpoint, "scenarios.rds")
if (file.exists(manifest)) stopifnot(identical(readRDS(manifest), cases)) else saveRDS(cases, manifest)
write.csv(cases, paste0(prefix, "_scenarios.csv"), row.names = FALSE)
write.csv(data.frame(item = seq_along(new), difficulty = new), paste0(prefix, "_new_bank.csv"), row.names = FALSE)
blocks <- split(seq_len(n_rep), ceiling(seq_len(n_rep) / if (stage == "large") 16L else 4L))
run_block <- function(i) {
  path <- file.path(checkpoint, sprintf("block_%03d.rds", i))
  if (file.exists(path)) return(path)
  rows <- lapply(blocks[[i]], run_replication)
  res <- list(native = do.call(rbind, lapply(rows, `[[`, "native")),
              external = do.call(rbind, lapply(rows, `[[`, "external")))
  tmp <- paste0(path, ".new")
  saveRDS(res, tmp)
  stopifnot(file.rename(tmp, path))
  cat(sprintf("%s block %d/%d, through replication %d\n", stage, i, length(blocks), max(blocks[[i]])))
  flush.console()
  path
}
started <- proc.time()[["elapsed"]]
cat("Package:", normalizePath(find.package("glmmplanr")), "\n")
paths <- parallel::mclapply(seq_along(blocks), run_block, mc.cores = workers,
                            mc.preschedule = FALSE, mc.set.seed = FALSE)
stopifnot(all(vapply(paths, function(p) is.character(p) && length(p) == 1L && file.exists(p), logical(1))))
results <- lapply(paths, readRDS)
fits <- as.data.frame(do.call(rbind, lapply(results, `[[`, "native")))
stopifnot(nrow(fits) == nrow(cases) * n_rep, !anyDuplicated(fits[c("case", "replication")]),
          all(table(fits$case) == n_rep))

summarize <- function(d, ref_type) {
  cfg <- cases[d$case[1L], ]
  x <- function(t) d[[paste0(ref_type, "_", t)]]
  hit <- sum(x("achieved"))
  ci <- unname(binom.test(hit, nrow(d))$conf.int)
  shift_error <- x("shift") - cfg$focal_mean
  cbind(cfg, reference = ref_type, attempted = nrow(d), native_fit = sum(d$native_fit),
        complete = sum(x("complete")), achieved = hit, mc_lo = ci[1], mc_hi = ci[2],
        decision = if (ci[1] >= .8) "pass" else if (ci[2] < .8) "fail" else "indeterminate",
        n_linked = sum(is.finite(shift_error)), mean_anchors_used = mean(x("anchors_used")),
        link_bias = mean(shift_error, na.rm = TRUE), link_rmse = sqrt(mean(shift_error^2, na.rm = TRUE)),
        origin_met = sum(abs(shift_error) <= .15, na.rm = TRUE),
        mean_rmse = mean(x("rmse"), na.rm = TRUE), mean_worst = mean(x("worst"), na.rm = TRUE),
        mean_max_tcc_error = mean(x("tcc_error"), na.rm = TRUE),
        unlinked_bias = mean(d$unlinked_bias, na.rm = TRUE),
        any_constant_anchor = sum(d$anchors_constant > 0), any_constant_new = sum(d$new_constant > 0))
}
summary <- do.call(rbind, lapply(split(fits, fits$case), function(d) {
  rbind(summarize(d, "oracle"), summarize(d, "estimated"))
}))
write.csv(summary, paste0(prefix, "_summary.csv"), row.names = FALSE)
con <- gzfile(paste0(prefix, "_fits.csv.gz"), "wt")
write.csv(fits, con, row.names = FALSE)
close(con)
if (stage == "large") {
  external <- as.data.frame(do.call(rbind, lapply(results, `[[`, "external")))
  stopifnot(nrow(external) == 384L, !anyDuplicated(external[c("case", "replication")]))
  write.csv(external, paste0(prefix, "_tam.csv"), row.names = FALSE)
}
cat(sprintf("Completed %d calibration datasets in %.1f seconds.\n", nrow(fits),
            proc.time()[["elapsed"]] - started))
