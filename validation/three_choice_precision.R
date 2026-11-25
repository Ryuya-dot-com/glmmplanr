# Conditional exposure coefficient in a balanced 3AFC study. Run from the repo root.
library(glmmplanr)
source("tests/testthat/helper-stan.R")
stopifnot(requireNamespace("posterior", quietly = TRUE), requireNamespace("rstan", quietly = TRUE))
dgp <- three_choice_reference_dgp()
spec <- specify_target(dgp, "exposure", engine = "bayes_laplace",
                       control = target_control(on_singular = "flag_only"))
designs <- data.frame(persons = c(32L, 48L, 72L, 48L), items = c(32L, 32L, 32L, 48L))
schedules <- Map(generate_design_schedule, design = list(dgp$design),
                 n_subj = designs$persons, n_item = designs$items)
jobs <- expand.grid(k = 1:4, trial = 1:192)
jobs$seed <- with(jobs, 7310000L + 1000L * k + trial)
idx <- vapply(list(dgp$design$X_cells, dgp$design$Z_subj_cells, dgp$design$Z_item_cells),
              function(x) match(spec$effect_col, colnames(x)), integer(1))
stopifnot(!anyNA(idx), !anyDuplicated(jobs$seed))
pars <- sprintf("%s[%d]", c("beta", "tau_subj", "tau_item"), idx)
trial_data <- function(j) {
  set.seed(jobs$seed[j], kind = "L'Ecuyer-CMRG", normal.kind = "Inversion", sample.kind = "Rejection")
  simulate_dgp_data(dgp, schedules[[jobs$k[j]]])
}
workers <- if (.Platform$OS.type == "windows") 1L else 8L
map_jobs <- function(rows, fun, preschedule = TRUE) {
  ans <- parallel::mclapply(rows, fun, mc.cores = workers, mc.set.seed = FALSE,
                            mc.preschedule = preschedule)
  stopifnot(all(vapply(ans, is.data.frame, logical(1))))
  do.call(rbind, ans)
}
# Record a checksum of each observed vector to verify pairing when regenerating it.
cpp_one <- function(j) {
  tr <- trial_data(j)
  f <- tryCatch(glmmplanr:::fit_single_dataset(tr, spec), error = function(e) {
    list(converged = FALSE, status = conditionMessage(e))
  })
  valid <- isTRUE(f$converged)
  data.frame(jobs[j, ], designs[jobs$k[j], ], truth = spec$true_beta,
             response_hash = digest::digest(tr$y, algo = "xxhash64"), cpp_valid = valid,
             cpp_status = f$status, cpp_singular = if (is.null(f$singular)) NA else f$singular,
             cpp_est = if (valid) f$est else NA_real_, cpp_se = if (valid) f$se else NA_real_,
             cpp_subj_sd = if (valid) sqrt(f$Sigma_subj[idx[2], idx[2]]) else NA_real_,
             cpp_item_sd = if (valid) sqrt(f$Sigma_item[idx[3], idx[3]]) else NA_real_)
}
native <- map_jobs(seq_len(nrow(jobs)), cpp_one)
stopifnot(identical(native$seed, jobs$seed), all(is.finite(native$truth)))
write.csv(native, "validation/three_choice_precision_native.csv", row.names = FALSE)
cat("Native fits complete:", nrow(native), "\n")
model <- rstan::stan_model(system.file("stan", "crossed_binomial_mafc.stan", package = "glmmplanr"))
run_hmc <- function(dat, seed, retry) {
  warnings <- character()
  result <- tryCatch(withCallingHandlers({
    chains <- if (retry) 4L else 2L
    iter <- if (retry) 1200L else 800L
    fit <- rstan::sampling(model, data = dat, seed = seed, chains = chains, iter = iter,
      warmup = iter %/% 2L, cores = 1L, refresh = 0L, pars = c("beta", "tau_subj", "tau_item"),
      control = list(adapt_delta = if (retry) 0.95 else 0.8, max_treedepth = 10L))
    draws <- rstan::extract(fit, permuted = FALSE, inc_warmup = FALSE, pars = pars)
    beta <- draws[, , pars[1]]
    stopifnot(identical(dim(beta), c(iter %/% 2L, chains)))
    sam <- do.call(rbind, rstan::get_sampler_params(fit, inc_warmup = FALSE))
    stopifnot(nrow(sam) == length(beta))
    rhat <- posterior::rhat(beta)
    bulk <- posterior::ess_bulk(beta)
    tail <- posterior::ess_tail(beta)
    div <- sum(sam[, "divergent__"])
    depth <- sum(sam[, "treedepth__"] >= 10L)
    checks <- c(rhat = rhat <= 1.01, bulk = bulk >= 100, tail = tail >= 100,
                divergent = div == 0, treedepth = depth == 0, sd = is.finite(sd(beta)) && sd(beta) > 0)
    valid <- all(!is.na(checks) & checks)
    data.frame(valid = valid, status = if (valid) "ok" else paste(names(checks)[is.na(checks) | !checks], collapse = "+"),
      est = mean(beta), se = sd(beta), subj_sd = mean(draws[, , pars[2]]),
      item_sd = mean(draws[, , pars[3]]), rhat = rhat, ess_bulk = bulk, ess_tail = tail,
      divergent = div, depth = depth)
  }, warning = function(w) {
    warnings <<- unique(c(warnings, conditionMessage(w)))
    invokeRestart("muffleWarning")
  }), error = function(e) data.frame(valid = FALSE, status = paste("error", conditionMessage(e)),
    est = NA_real_, se = NA_real_, subj_sd = NA_real_, item_sd = NA_real_, rhat = NA_real_,
    ess_bulk = NA_real_, ess_tail = NA_real_, divergent = NA_integer_, depth = NA_integer_))
  result$warnings <- paste(gsub("[\r\n]+", " ", warnings), collapse = " | ")
  result
}
stan_one <- function(j) {
  tr <- trial_data(j)
  stopifnot(identical(digest::digest(tr$y, algo = "xxhash64"), native$response_hash[j]))
  dat <- three_choice_stan_data(tr, spec)
  seed <- jobs$seed[j] + 100000L
  first <- run_hmc(dat, seed, FALSE)
  last <- if (first$valid) first else run_hmc(dat, seed + 100000L, TRUE)
  ans <- cbind(seed = jobs$seed[j], retried = !first$valid,
               setNames(first, paste0("initial_", names(first))),
               setNames(last, paste0("stan_", names(last))))
  saveRDS(ans, file.path(checkpoints, paste0(jobs$seed[j], ".rds")))
  ans
}
# Individual completed fits survive an interruption in the isolated run directory.
checkpoints <- "validation/three_choice_hmc_checkpoints"
dir.create(checkpoints, showWarnings = FALSE)
ref <- map_jobs(which(jobs$trial <= 40L), stan_one, preschedule = FALSE)
stopifnot(identical(ref$seed, jobs$seed[jobs$trial <= 40L]))
write.csv(ref, "validation/three_choice_precision_stan.csv", row.names = FALSE)
cat("Stan fits complete:", nrow(ref), "\n")
