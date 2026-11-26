# Replay saved failures and successful controls; no replacement of old evidence.
library(glmmplanr)
args <- commandArgs(TRUE)
out <- if (length(args)) args[1L] else "validation/rasch_mml_acceleration"
references_only <- length(args) > 1L && args[2L] == "references"
if (dir.exists(out) && !references_only) stop("Use a new output directory; existing evidence is not overwritten")
if (!dir.exists(out)) dir.create(out, recursive = TRUE)
source("validation/rasch_fit_benchmark_helpers.R")
cases <- read.csv("validation/rasch_fit_benchmark/calibration_scenarios.csv")
old <- read.csv("validation/rasch_fit_benchmark/calibration_calibration_fits.csv.gz")
native <- subset(old, engine == "glmmplanr")
failed <- subset(native, !valid & error == "MML-EM reached iteration limit")
controls <- do.call(rbind, lapply(split(subset(native, valid), native$case[native$valid]), head, 8L))
jobs <- rbind(failed, controls)
forms <- as.matrix(read.csv("validation/two_pl_fixed_bank_booklets_forms.csv"))
streams <- glmmplanr:::make_lecuyer_streams(unique(cases$seed), unique(cases$n_rep))
# Predeclared comparison tolerances: 0.002 logits for b and sigma, 0.002
# log-likelihood units (same 41-point integral). Decisions use unrounded errors.
refit <- function(k) {
  x <- jobs[k, ]
  case <- cases[cases$case == x$case, ]
  truth <- rasch_data(case, streams[[x$replicate]], forms)
  long <- long_rasch(truth$y)
  timing <- system.time(f <- tryCatch(fit_native(truth$y, long), error = identity))
  ok <- !inherits(f, "error")
  reference <- subset(old, engine == "TAM" & case == x$case & replicate == x$replicate)
  data.frame(case = x$case, replicate = x$replicate, old_valid = x$valid,
    valid = ok, error = if (ok) NA_character_ else conditionMessage(f),
    sigma = if (ok) f$sigma else NA_real_,
    b_rmse = if (ok) sqrt(mean((f$b - truth$b)^2)) else NA_real_,
    b_max = if (ok) max(abs(f$b - truth$b)) else NA_real_,
    ll_common = if (ok) rasch_ll(truth$y, rasch_patterns(truth$y), f$b, f$sigma) else NA_real_,
    old_b_rmse = x$b_rmse, old_sigma = x$sigma, old_ll = x$ll_common,
    tam_valid = reference$valid, tam_rmse = reference$b_rmse,
    tam_sigma = reference$sigma, tam_ll = reference$ll_common,
    tam_bank_ok = reference$bank_ok, cpu_s = sum(timing[1:2]))
}
if (!references_only) {
rows <- parallel::mclapply(seq_len(nrow(jobs)), refit,
  mc.cores = if (.Platform$OS.type == "windows") 1L else 2L, mc.set.seed = FALSE)
stopifnot(all(vapply(rows, is.data.frame, logical(1L))))
fits <- do.call(rbind, rows)
fits$bank_ok <- with(fits, valid & b_rmse <= .35 & b_max <= .6)
write.csv(fits, file.path(out, "replay.csv"), row.names = FALSE)
print(aggregate(cbind(valid, bank_ok) ~ old_valid, fits, sum))
print(aggregate(cpu_s ~ old_valid, fits, sum))
} else fits <- read.csv(file.path(out, "replay.csv"))

# Independent likelihood optimization for the first recovered dataset per
# case and the smallest positive variance and boundary cases. High-resolution
# integration and a different optimizer supply a reference, not a new engine.
pick <- unique(c(unlist(lapply(split(which(!fits$old_valid), fits$case[!fits$old_valid]), head, 1L)),
                 which(fits$case == 5L & fits$replicate == 329L),
                 head(which(fits$valid & fits$sigma == 0), 1L),
                 head(order(ifelse(fits$sigma > 0, fits$sigma, Inf)), 2L)))
references <- lapply(pick, function(k) {
  x <- fits[k, ]
  truth <- rasch_data(cases[cases$case == x$case, ], streams[[x$replicate]], forms)
  f <- tryCatch(fit_native(truth$y, long_rasch(truth$y)), error = function(e)
    list(b = rep(NA_real_, 16L), sigma = NA_real_))
  pat <- rasch_patterns(truth$y)
  objective <- function(par) -rasch_ll(truth$y, pat, par[1:16], sqrt(par[17]))
  opt <- optim(c(-qlogis(colMeans(truth$y, na.rm = TRUE)), .16), objective, method = "L-BFGS-B",
               lower = c(rep(-Inf, 16L), 0), control = list(factr = 1e3, pgtol = 1e-7, maxit = 2000L))
  r <- data.frame(case = x$case, replicate = x$replicate, engine = "direct_41",
    valid = opt$convergence == 0L, reference_sigma = sqrt(opt$par[17]),
    b_max_difference = max(abs(f$b - opt$par[1:16])),
    sigma_difference = f$sigma - sqrt(opt$par[17]),
    ll_difference = rasch_ll(truth$y, pat, f$b, f$sigma) + opt$value)
  if (requireNamespace("mirt", quietly = TRUE)) {
    m <- tryCatch(mirt::mirt(truth$y, 1L, itemtype = "Rasch", SE = FALSE,
      quadpts = 61L, verbose = FALSE, TOL = 1e-6, technical = list(NCYCLES = 2000L)), error = identity)
    if (!inherits(m, "error")) {
      cf <- mirt::coef(m, simplify = TRUE)
      mb <- -cf$items[, "d"]
      ms <- sqrt(cf$cov[1, 1])
      r <- rbind(r, data.frame(case = x$case, replicate = x$replicate, engine = "mirt",
        valid = isTRUE(m@OptimInfo$converged), reference_sigma = ms,
        b_max_difference = max(abs(f$b - mb)),
        sigma_difference = f$sigma - ms,
        ll_difference = rasch_ll(truth$y, pat, f$b, f$sigma) - rasch_ll(truth$y, pat, mb, ms)))
    } else r <- rbind(r, data.frame(case = x$case, replicate = x$replicate,
      engine = "mirt", valid = FALSE, reference_sigma = NA,
      b_max_difference = NA, sigma_difference = NA, ll_difference = NA))
  }
  r
})
write.csv(do.call(rbind, references), file.path(out, "references.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out, "environment.txt"))
write.csv(data.frame(file = c("src/glmm_laplace.cpp", "R/rasch_precision.R"),
  md5 = unname(tools::md5sum(c("src/glmm_laplace.cpp", "R/rasch_precision.R")))),
  file.path(out, "source.csv"), row.names = FALSE)
