# Empirical-bank targeting pilot: native Rasch grid, paired mirt sentinels,
# and a separate correctly specified 2PL arm (not a native 2PL capability).
source("validation/planning_comparison_helpers.R")
args <- commandArgs(TRUE)
mode <- if (length(args)) match.arg(args[1L], c("smoke", "pilot")) else "pilot"
out <- if (length(args) > 1L) args[2L] else paste0("validation/l2_targeting_", mode)
workers <- if (length(args) > 2L) as.integer(args[3L]) else 4L
if (.Platform$OS.type == "windows") workers <- 1L
bank <- read.csv("validation/l2_empirical_bank/items.csv")
empirical <- read.csv("validation/l2_empirical_bank/summary.csv")
booklets <- read.csv("validation/l2_empirical_bank/booklets.csv")
profiles <- data.frame(profile = c("reference", "narrow", "lower", "lower_narrow", "upper", "severe_lower"),
  mean = c(0, 0, -3, -3, 3, -5), sd_multiplier = c(1, .5, 1, .5, 1, 1))
cases <- merge(expand.grid(format = c("Recall", "Recognition"),
  layout = c("full50", "balanced25", "difficult25", "booklets_equal", "booklets_observed"),
  n_subj = c(100L, 200L, 500L, 1000L), stringsAsFactors = FALSE), profiles)
cases$case <- seq_len(nrow(cases))
cases$external <- with(cases, n_subj %in% c(200L, 1000L) &
  profile %in% c("reference", "lower_narrow", "severe_lower") &
  layout %in% c("full50", "booklets_observed"))
if (mode == "smoke") cases <- cases[c(1L, which(cases$external)[1L],
  tail(which(cases$external), 1L)), ]
n_rep <- if (mode == "smoke") 2L else 64L
n_external <- if (mode == "smoke") 2L else 32L
checkpoint_dir(out, c("validation/l2_targeting_simulation.R", "validation/l2_empirical_bank/items.csv",
  "validation/l2_empirical_bank/summary.csv", "validation/l2_empirical_bank/booklets.csv"),
  list(mode = mode, cases = cases, n_rep = n_rep, n_external = n_external,
    seed = 96720001L, rmse = .35, maximum = .60, probability_rmse = .05), c("glmmplanr", "mirt"))
write.csv(cases, file.path(out, "scenarios.csv"), row.names = FALSE)
run <- function(k) {
  x <- cases[k, ]; path <- file.path(out, sprintf("case-%03d.rds", x$case))
  if (file.exists(path)) return(readRDS(path)$summary)
  ip <- subset(bank, format == x$format & model == "Rasch")
  if (x$layout != "full50" && !grepl("booklets", x$layout)) {
    keep <- unlist(lapply(split(seq_len(nrow(ip)), ip$frequency_band), function(i)
      if (x$layout == "balanced25") i[c(1L, 3L, 5L, 7L, 9L)] else i[order(ip$b[i], decreasing = TRUE)[1:5]]))
    ip <- ip[sort(keep), ]
  }
  j <- nrow(ip); n <- x$n_subj
  mask <- matrix(TRUE, n, j)
  if (grepl("booklets", x$layout)) {
    bs <- subset(booklets, format == x$format)
    bm <- matrix(bs$present, nrow = 50L, ncol = 5L)
    weights <- if (x$layout == "booklets_equal") rep(1, 5L) else unique(bs[c("booklet", "observed_booklet_n")])$observed_booklet_n
    quota <- n * weights / sum(weights); counts <- floor(quota)
    remainder <- n - sum(counts)
    if (remainder) counts[order(quota - counts, decreasing = TRUE)[seq_len(remainder)]] <-
      counts[order(quota - counts, decreasing = TRUE)[seq_len(remainder)]] + 1L
    mask <- t(bm[, rep(seq_len(5L), counts), drop = FALSE]) == 1L
  }
  idx <- which(mask, arr.ind = TRUE)
  sigma <- sqrt(subset(empirical, format == x$format & model == "Rasch")$person_var) * x$sd_multiplier
  streams <- glmmplanr:::make_lecuyer_streams(96720001L + 1000L * x$case, n_rep)
  rows <- list(); row_index <- 0L
  assess <- function(y, model, engine, a_true, b_true, fit_b, fit_a, valid, error) {
    err_b <- fit_b - b_true
    theta_grid <- x$mean + (if (model == "Rasch") sigma else x$sd_multiplier) * seq(-2, 2, length.out = 9L)
    truth <- plogis(outer(theta_grid, b_true, "-") * rep(a_true, each = length(theta_grid)))
    estimate <- plogis(outer(theta_grid, fit_b, "-") * rep(fit_a, each = length(theta_grid)))
    n_obs <- colSums(!is.na(y)); correct <- colSums(y, na.rm = TRUE)
    by_person <- rowSums(y, na.rm = TRUE); person_obs <- rowSums(!is.na(y))
    data.frame(case = x$case, replicate = r, model, engine, valid,
      rmse = if (valid) sqrt(mean(err_b^2)) else NA_real_,
      max_error = if (valid) max(abs(err_b)) else NA_real_,
      a_rmse = if (valid) sqrt(mean((fit_a - a_true)^2)) else NA_real_,
      a_max_error = if (valid) max(abs(fit_a - a_true)) else NA_real_,
      probability_rmse = if (valid) sqrt(mean((estimate - truth)^2)) else NA_real_,
      success = valid && sqrt(mean(err_b^2)) <= .35 && max(abs(err_b)) <= .60,
      parameter_success = valid && sqrt(mean(err_b^2)) <= .35 && max(abs(err_b)) <= .60 &&
        (model == "Rasch" || (sqrt(mean((fit_a - a_true)^2)) <= .25 && max(abs(fit_a - a_true)) <= .50)),
      probability_success = valid && sqrt(mean((estimate - truth)^2)) <= .05,
      n_extreme_items = sum(correct == 0 | correct == n_obs),
      n_extreme_people = sum(by_person == 0 | by_person == person_obs),
      mean_response = mean(y, na.rm = TRUE), exposure_min = min(n_obs), exposure_max = max(n_obs),
      error, data_hash = digest::digest(y, algo = "xxhash64"))
  }
  for (r in seq_len(n_rep)) {
    dat <- glmmplanr:::with_rng_stream(streams[[r]], {
      u <- rnorm(n); theta <- x$mean + sigma * u
      y <- matrix(rbinom(n * j, 1, plogis(outer(theta, ip$b, "-"))), n, j)
      y[!mask] <- NA
      list(u = u, y = y)
    })
    y <- dat$y
    tr <- data.frame(subj_id = idx[, 1L], item_id = idx[, 2L], y = y[idx])
    eligible <- colSums(y, na.rm = TRUE) > 0 & colSums(y, na.rm = TRUE) < colSums(mask)
    native <- glmmplanr:::fit_rasch_fixed_mml(tr, n, j, eligible)
    valid <- is.na(native$error) && all(is.finite(native$b_est))
    row_index <- row_index + 1L
    rows[[row_index]] <- assess(y, "Rasch", "native", rep(1, j), ip$b,
      native$b_est + x$mean, rep(1, j), valid, native$error)
    if (isTRUE(x$external) && r <= n_external) {
      for (model in c("Rasch", "2PL")) {
        p <- if (model == "Rasch") ip else subset(bank, format == x$format & model == "2PL")
        if (model == "2PL") {
          # Known generating location/scale is used only for recovery scoring.
          y <- glmmplanr:::with_rng_stream(streams[[r]], {
            u <- rnorm(n)
            eta <- outer(x$mean + x$sd_multiplier * u, p$b, "-") * rep(p$a, each = n)
            ans <- matrix(rbinom(n * j, 1, plogis(eta)), n, j); ans[!mask] <- NA; ans
          })
        }
        colnames(y) <- p$item_id
        set.seed(96750001L + 10000L * x$case + 100L * r + match(model, c("Rasch", "2PL")))
        constant <- any(colSums(y, na.rm = TRUE) %in% 0L | colSums(y, na.rm = TRUE) == colSums(!is.na(y)))
        fit <- if (constant) simpleError("Constant item: complete bank not estimable") else
          tryCatch(suppressWarnings(mirt::mirt(y, 1, itemtype = model, verbose = FALSE)), error = identity)
        good <- !inherits(fit, "error") && isTRUE(mirt::extract.mirt(fit, "converged"))
        est <- if (good) mirt::coef(fit, IRTpars = TRUE, simplify = TRUE)$items else matrix(NA_real_, j, 4,
          dimnames = list(NULL, c("a", "b", "g", "u")))
        good <- good && all(is.finite(est[, c("a", "b")]))
        fit_b <- est[, "b"] * (if (model == "2PL") x$sd_multiplier else 1) + x$mean
        fit_a <- est[, "a"] / (if (model == "2PL") x$sd_multiplier else 1)
        # Hash only response values and missingness, independent of column labels.
        dimnames(y) <- NULL
        row_index <- row_index + 1L
        rows[[row_index]] <- assess(y, model, "mirt", p$a, p$b, fit_b, fit_a, good,
          if (inherits(fit, "error")) conditionMessage(fit) else if (!good) "Nonconvergence" else NA_character_)
      }
    }
  }
  rows <- do.call(rbind, rows)
  summary <- do.call(rbind, lapply(split(rows, interaction(rows$model, rows$engine, drop = TRUE)), function(z)
    data.frame(x, model = z$model[1L], engine = z$engine[1L], items = j,
      person_mean = x$mean, person_sd = if (z$model[1L] == "Rasch") sigma else x$sd_multiplier,
      mc_summary(z$success, z$valid), mean_rmse = mean(z$rmse, na.rm = TRUE),
      mean_a_rmse = mean(z$a_rmse, na.rm = TRUE),
      parameter_attainment = mean(z$parameter_success),
      parameter_mc_lo = glmmplanr:::wilson_mc_ci(sum(z$parameter_success), nrow(z))[1L],
      parameter_mc_hi = glmmplanr:::wilson_mc_ci(sum(z$parameter_success), nrow(z))[2L],
      probability_attainment = mean(z$probability_success),
      mean_probability_rmse = mean(z$probability_rmse, na.rm = TRUE),
      mean_response = mean(z$mean_response), mean_extreme_items = mean(z$n_extreme_items),
      mean_extreme_people = mean(z$n_extreme_people), exposure_min = min(z$exposure_min), exposure_max = max(z$exposure_max))))
  atomic_rds(list(summary = summary, trials = rows), path)
  cat(sprintf("Targeting %d/%d (case %d) complete\n", k, nrow(cases), x$case))
  summary
}
result <- parallel::mclapply(seq_len(nrow(cases)), run, mc.cores = workers, mc.set.seed = FALSE)
stopifnot(all(vapply(result, is.data.frame, logical(1))))
write.csv(do.call(rbind, result), file.path(out, "summary.csv"), row.names = FALSE)
trials <- do.call(rbind, lapply(cases$case, function(id) readRDS(file.path(out, sprintf("case-%03d.rds", id)))$trials))
con <- gzfile(file.path(out, "fits.csv.gz"), "wt"); write.csv(trials, con, row.names = FALSE); close(con)
native <- subset(trials, model == "Rasch" & engine == "native")
external <- subset(trials, model == "Rasch" & engine == "mirt")
paired <- merge(native, external, by = c("case", "replicate"), suffixes = c("_native", "_mirt"))
stopifnot(nrow(paired) == nrow(external), all(paired$data_hash_native == paired$data_hash_mirt))
write.csv(paired, file.path(out, "paired_rasch.csv"), row.names = FALSE)
