# Published-mean planning example, not reanalysis of unavailable Huang data.
source("validation/planning_comparison_helpers.R")
args <- commandArgs(TRUE)
mode <- if (length(args)) match.arg(args[1L], c("smoke", "full", "large_smoke", "large")) else "full"
large <- mode %in% c("large_smoke", "large")
smoke <- mode %in% c("smoke", "large_smoke")
out <- if (length(args) > 1L) args[2L] else paste0("validation/l2_replication_", mode)
workers <- if (length(args) > 2L) as.integer(args[3L]) else 4L
if (.Platform$OS.type == "windows") workers <- 1L
profiles <- data.frame(profile = c("published_8pp", "slopes_8pp", "conservative_5pp", "high_variance_8pp"),
  effect_pp = c(8, 8, 5, 8), subj_intercept = c(.8, .8, .8, 1.1),
  item_intercept = c(.6, .6, .6, .9), subj_slope = c(0, .25, .25, .4),
  item_slope = c(0, .35, .35, .5))
jobs <- merge(profiles, data.frame(preference = c("fewer_subj", "fewer_items")))
if (large) jobs <- merge(jobs, data.frame(search_seed = seq_len(5L)))
jobs$job <- seq_len(nrow(jobs))
if (smoke) jobs <- jobs[c(1L, nrow(jobs)), ]
make_spec <- function(x) {
  des <- design_crossed(list(study_method = c("restudy", "retrieval")),
    allocation = "latin_square", items_per_subject = 30L, item_seed = 37L)
  dgp <- calibrate_dgp(des, c(restudy = .42, retrieval = .42 + x$effect_pp / 100),
    subj_sd = c(x$subj_intercept, x$subj_slope), item_sd = c(x$item_intercept, x$item_slope),
    corr_structure = "zero")
  specify_target(dgp, "study_method", scale = "pp", criterion = "joint", ci_width = 10,
    engine = "bayes_laplace", analysis_re = "match_dgp", prob_se_method = "joint_delta")
}
grid <- expand.grid(n_subj = c(54L, 72L, 108L, 162L, 216L), n_item = c(30L, 60L, 90L))
grid$i <- match(grid$n_subj, unique(grid$n_subj)); grid$j <- match(grid$n_item, unique(grid$n_item))
checkpoint_dir(out, "validation/l2_replication_planning.R", list(mode = mode, jobs = jobs,
  grid = grid, nsim_confirm = if (large) 2048L else 512L,
  seed = if (large) 610100001L else 96810001L,
  confirmation_seed = if (large) 810900001L else 96890001L, budget = 1080L, ci_width = 10),
  c("glmmplanr", "mlpwr", "DiceKriging", "randtoolbox", "rngWELL", "rgenoud"))
write.csv(jobs, file.path(out, "scenarios.csv"), row.names = FALSE)
run <- function(k) {
  x <- jobs[k, ]; path <- file.path(out, sprintf("search-%02d.rds", x$job))
  if (file.exists(path)) return(readRDS(path)$selection)
  spec <- make_spec(x)
  costs <- if (x$preference == "fewer_subj") (grid$n_subj - 54) / 162 + .02 * (grid$n_item - 30) / 60 else
    (grid$n_item - 30) / 60 + .02 * (grid$n_subj - 54) / 162
  costfun <- function(i, j) costs[match(paste(i, j), paste(grid$i, grid$j))]
  runners <- lapply(seq_len(nrow(grid)), function(g) {
    schedule <- generate_design_schedule(spec$dgp$design, grid$n_subj[g], grid$n_item[g])
    cache <- glmmplanr:::prepare_cpp_schedule_cache(spec, schedule)
    function() {
      fit <- glmmplanr:::execute_single_replication(spec, schedule, cache,
        "bayes_laplace", spec$analysis_re, FALSE, FALSE)$primary
      valid <- isTRUE(glmmplanr:::valid_target_fit(fit$est, fit$se, fit$converged, fit$status,
        fit$prob_est, fit$prob_se, "pp"))
      data.frame(valid, success = valid && abs(fit$prob_est / fit$prob_se) > qnorm(.975) &&
        2 * qnorm(.975) * fit$prob_se * 100 <= 10)
    }
  })
  seed <- if (large) 610100001L + 100000L * x$job else 96810001L + 10000L * x$job
  audit_stream <- glmmplanr:::make_lecuyer_streams(seed + 101L, 1L)[[1L]]
  audit <- glmmplanr:::with_rng_stream(audit_stream, runners[[1L]]())
  public <- evaluate_design_point(spec, grid$n_subj[1L], grid$n_item[1L], nsim = 1L,
    stream_seeds = list(audit_stream))
  stopifnot(audit$valid == (public$n_valid == 1L), audit$success == (public$k_joint == 1L))
  budget <- if (smoke) 120L else 1080L
  streams <- glmmplanr:::make_lecuyer_streams(seed, budget + 60L)
  observed <- list(); count <- 0L
  evaluate <- function(i, j) {
    g <- match(paste(as.integer(i), as.integer(j)), paste(grid$i, grid$j))
    stopifnot(!is.na(g)); count <<- count + 1L
    z <- glmmplanr:::with_rng_stream(streams[[count]], runners[[g]]())
    observed[[count]] <<- data.frame(pos = g, z); isTRUE(z$success)
  }
  # Explicit unique starts through the public API avoid duplicate rounded
  # initial points on a small grid. Their cost is included in the SAME budget.
  initial <- lapply(seq_len(nrow(grid)), function(g)
    list(x = c(i = grid$i[g], j = grid$j[g]), y = replicate(if (smoke) 2L else 10L,
      evaluate(grid$i[g], grid$j[g]))))
  set.seed(seed + 201L); warnings <- character()
  ml <- tryCatch(withCallingHandlers(mlpwr::find.design(evaluate,
    boundaries = list(i = c(1, 5), j = c(1, 3)), dat = initial,
    power = .8, costfun = costfun, evaluations = budget - count,
    setsize = if (smoke) 10L else 30L, surrogate = "gpr", silent = TRUE),
    warning = function(w) {warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")}),
    error = identity)
  stopifnot(count >= 30L)
  counts <- rep(count %/% nrow(grid), nrow(grid))
  if (count %% nrow(grid)) counts[seq_len(count %% nrow(grid))] <- counts[seq_len(count %% nrow(grid))] + 1L
  gs <- do.call(rbind, lapply(seq_len(nrow(grid)), function(g) {
    z <- evaluate_design_point(spec, grid$n_subj[g], grid$n_item[g], nsim = counts[g],
      stream_seeds = seed + 300L + g)
    data.frame(grid[g, ], cost = costs[g], nsim = counts[g], n_valid = z$n_valid,
      n_met = z$k_joint, attainment = z$k_joint / counts[g])
  }))
  stopifnot(sum(gs$nsim) == count)
  eligible <- which(gs$attainment >= .8)
  gp <- if (length(eligible)) eligible[which.min(costs[eligible])] else NA_integer_
  mp <- NA_integer_
  if (!inherits(ml, "error") && is.finite(ml$final$power) && ml$final$power >= .8 &&
      !any(grepl("No good design found", warnings)))
    mp <- match(paste(ml$final$design$i, ml$final$design$j), paste(grid$i, grid$j))
  selected <- do.call(rbind, lapply(c("glmmplanr_grid", "mlpwr"), function(method) {
    pos <- if (method == "mlpwr") mp else gp
    data.frame(x, method, evaluations = count, n_subj = if (!is.na(pos)) grid$n_subj[pos] else NA_integer_,
      n_item = if (!is.na(pos)) grid$n_item[pos] else NA_integer_,
      error = if (method == "mlpwr" && inherits(ml, "error")) conditionMessage(ml) else "",
      warning = if (method == "mlpwr") paste(unique(warnings), collapse = " | ") else "")
  }))
  atomic_rds(list(selection = selected, grid = gs, mlpwr_observations = do.call(rbind, observed),
    mlpwr_final = if (!inherits(ml, "error")) ml$final else NULL), path)
  cat(sprintf("Replication search %d: %d evaluations per method\n", x$job, count))
  selected
}
search <- parallel::mclapply(seq_len(nrow(jobs)), run, mc.cores = workers,
  mc.set.seed = FALSE, mc.preschedule = !large)
stopifnot(all(vapply(search, is.data.frame, logical(1))))
search <- do.call(rbind, search)
write.csv(search, file.path(out, "search.csv"), row.names = FALSE)
# A single illustrative search per preference/profile, followed by genuinely
# independent planned/null confirmation; no claim about optimizer seed stability.
selected <- unique(subset(search, !is.na(n_subj), c("profile", "n_subj", "n_item")))
if (large && !smoke) selected <- merge(profiles["profile"], grid[c("n_subj", "n_item")])
selected <- selected[order(selected$profile, selected$n_subj, selected$n_item), ]
selected$case <- seq_len(nrow(selected))
write.csv(selected, file.path(out, "confirmation_scenarios.csv"), row.names = FALSE)
validate_with_trials <- function(spec, n, j, nsim, seed) {
  # Local copies retain the public calculation, while recording its inputs
  # to the summary. The installed namespace and estimator remain untouched.
  capture <- new.env(parent = asNamespace("glmmplanr"))
  capture$records <- list()
  capture$summarize_mc_metrics <- function(...) {
    z <- list(...)
    capture$records[[length(capture$records) + 1L]] <- data.frame(
      replicate = seq_along(z$est), est_logit = z$est, se_logit = z$se,
      estimate = 100 * z$prob_est, se = 100 * z$prob_se,
      valid = glmmplanr:::valid_target_fit(z$est, z$se, z$converged,
        z$status, z$prob_est, z$prob_se, "pp"), status = z$status,
      error = z$error_reasons)
    glmmplanr:::summarize_mc_metrics(...)
  }
  evaluator <- glmmplanr::evaluate_design_point
  environment(evaluator) <- capture
  capture$evaluate_design_point <- evaluator
  validator <- glmmplanr::validate_design_point
  environment(validator) <- capture
  check <- validator(spec, n, j, nsim = nsim, seed = seed)
  stopifnot(length(capture$records) == nrow(check))
  trials <- do.call(rbind, lapply(seq_len(nrow(check)), function(i) {
    z <- capture$records[[i]]; z$condition <- check$condition[i]
    z$truth <- check$true_target[i]
    z$detected <- z$valid & abs(z$estimate / z$se) > qnorm(.975)
    z$width_ok <- z$valid & 2 * qnorm(.975) * z$se <= 10
    z$covered <- z$valid & abs(z$estimate - z$truth) <= qnorm(.975) * z$se
    stopifnot(sum(z$valid) == check$n_valid[i],
      abs(mean(z$detected) - check$detection[i]) < 1e-12,
      abs(mean(z$detected & z$width_ok) - check$joint_attainment[i]) < 1e-12,
      abs(sum(z$covered) / sum(z$valid) - check$coverage[i]) < 1e-12)
    z
  }))
  if (smoke) stopifnot(identical(check,
    glmmplanr::validate_design_point(spec, n, j, nsim = nsim, seed = seed)))
  list(check = check, trials = trials)
}
confirm <- function(k) {
  x <- selected[k, ]; path <- file.path(out, sprintf("confirm-%02d.rds", k))
  if (file.exists(path)) return(readRDS(path)$summary)
  spec <- make_spec(profiles[profiles$profile == x$profile, ])
  n_rep <- if (smoke) 4L else if (large) 2048L else 512L
  confirmation_seed <- if (large) 810900001L + 10000L * k else 96890001L + 1000L * k
  detail <- if (large) validate_with_trials(spec, x$n_subj, x$n_item, n_rep, confirmation_seed) else NULL
  check <- if (large) detail$check else
    validate_design_point(spec, x$n_subj, x$n_item, nsim = n_rep, seed = confirmation_seed)
  summary <- data.frame(x, as.data.frame(check))
  if (large) {
    diagnostics <- do.call(rbind, lapply(summary$condition, function(condition) {
      z <- detail$trials[detail$trials$condition == condition & detail$trials$valid, ]
      data.frame(bias_pp = mean(z$estimate - z$truth),
        rmse_pp = sqrt(mean((z$estimate - z$truth)^2)),
        empirical_sd_pp = sd(z$estimate), mean_se_pp = mean(z$se),
        mean_width_pp = mean(2 * qnorm(.975) * z$se))
    }))
    summary <- cbind(summary, diagnostics)
  }
  n_met <- round(summary$joint_attainment * n_rep)
  ci <- t(vapply(n_met, glmmplanr:::wilson_mc_ci, numeric(2), n = n_rep))
  summary$joint_lo <- ci[, 1L]; summary$joint_hi <- ci[, 2L]
  summary$joint_decision <- ifelse(ci[, 1L] >= .8, "pass", ifelse(ci[, 2L] < .8, "fail", "indeterminate"))
  atomic_rds(list(summary = summary, check = check,
    trials = if (large) detail$trials else NULL), path)
  cat(sprintf("Replication confirmation %d/%d complete\n", k, nrow(selected)))
  summary
}
checks <- parallel::mclapply(seq_len(nrow(selected)), confirm, mc.cores = workers,
  mc.set.seed = FALSE, mc.preschedule = !large)
stopifnot(all(vapply(checks, is.data.frame, logical(1))))
write.csv(do.call(rbind, checks), file.path(out, "confirmation.csv"), row.names = FALSE)
