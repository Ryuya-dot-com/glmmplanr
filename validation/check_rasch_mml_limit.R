# Refit only the binary-bank datasets that hit the release EM iteration limit.
args <- commandArgs(TRUE)
out <- if (length(args)) args[1L] else "validation/rasch_fit_benchmark"
workers <- if (length(args) > 1L) as.integer(args[2L]) else 2L
stopifnot(workers > 0L, as.character(packageVersion("glmmplanr")) == "0.1.0")
if (.Platform$OS.type == "windows") workers <- 1L
source("validation/rasch_fit_benchmark_helpers.R")
cases <- read.csv(file.path(out, "calibration_scenarios.csv"))
fits <- read.csv(file.path(out, "calibration_calibration_fits.csv.gz"))
need <- fits[fits$engine == "glmmplanr" & fits$error == "MML-EM reached iteration limit", ]
if (!nrow(need)) {
  message("No fits reached the release iteration limit")
  quit(save = "no")
}
mfr <- fits[fits$engine == "mfrmr", c("case", "replicate", "sigma", "ll_common", "bank_ok")]
need <- merge(need[c("case", "replicate")], mfr, by = c("case", "replicate"), all.x = TRUE)
forms <- as.matrix(read.csv("validation/two_pl_fixed_bank_booklets_forms.csv"))
streams <- glmmplanr:::make_lecuyer_streams(unique(cases$seed), unique(cases$n_rep))
refit <- function(k) {
  x <- need[k, ]
  case <- cases[cases$case == x$case, ]
  truth <- rasch_data(case, streams[[x$replicate]], forms)
  long <- long_rasch(truth$y)
  trials <- data.frame(subj_id = long$id, item_id = long$item_id, y = long$y)
  start <- proc.time()
  fit <- glmmplanr:::fit_rasch_fixed_mml(trials, nrow(truth$y), ncol(truth$y),
                                        rep(TRUE, ncol(truth$y)), max_iter = 5000L)
  t <- proc.time() - start
  ok <- is.na(fit$error)
  b0 <- -qlogis(colMeans(truth$y, na.rm = TRUE))
  p0 <- plogis(-b0)
  observed <- !is.na(truth$y)
  # Profile likelihood right derivative with respect to the variance
  # at zero, in the model with independently estimated item intercepts.
  centered_score <- rowSums(truth$y, na.rm = TRUE) - drop(observed %*% p0)
  boundary_score <- .5 * sum(centered_score^2 - drop(observed %*% (p0 * (1 - p0))))
  sigma <- rmse <- bank_ok <- ll <- NA_real_
  if (ok) {
    sigma <- sqrt(fit$person_var)
    rmse <- sqrt(mean((fit$b_est - truth$b)^2))
    bank_ok <- as.integer(rmse <= .35 && max(abs(fit$b_est - truth$b)) <= .6)
    ll <- rasch_ll(truth$y, rasch_patterns(truth$y), fit$b_est, sigma)
  }
  data.frame(case = x$case, replicate = x$replicate, extended_valid = ok, error = fit$error,
             boundary_score, extended_sigma = sigma, mfrmr_sigma = x$sigma,
             extended_rmse = rmse, extended_bank_ok = bank_ok, mfrmr_bank_ok = x$bank_ok,
             extended_ll_minus_mfrmr = ll - x$ll_common,
             cpu_s = unname(sum(t[c("user.self", "sys.self")])))
}
rows <- parallel::mclapply(seq_len(nrow(need)), refit, mc.cores = workers, mc.set.seed = FALSE)
stopifnot(length(rows) > 0L, all(vapply(rows, is.data.frame, logical(1L))))
write.csv(do.call(rbind, rows), file.path(out, "calibration_iteration_limit.csv"), row.names = FALSE)
