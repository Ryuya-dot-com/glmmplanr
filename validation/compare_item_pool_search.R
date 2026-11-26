# Equal-budget grid evaluation versus mlpwr proposals, using native estimators.
source("validation/planning_comparison_helpers.R")
args <- commandArgs(TRUE)
mode <- if (length(args)) match.arg(args[1L], c("smoke", "full")) else "full"
out <- if (length(args) > 1L) args[2L] else "validation/item_pool_search"
workers <- if (length(args) > 2L) as.integer(args[3L]) else 2L
if (.Platform$OS.type == "windows") workers <- 1L
stopifnot(requireNamespace("mlpwr", quietly = TRUE))
families <- unique(rasch_cases$family)
jobs <- rbind(data.frame(task = "rasch", family = families, preference = "n_subj"),
  expand.grid(task = "glmm", family = c("zero", "heterogeneous"),
              preference = c("fewer_subj", "fewer_items"), stringsAsFactors = FALSE))
if (mode == "smoke") jobs <- jobs[c(1L, 5L, 11L, 12L), ]
jobs <- merge(jobs, data.frame(run = seq_len(if (mode == "full") 5L else 1L)))
jobs$job <- seq_len(nrow(jobs))
checkpoint_dir(out, "validation/compare_item_pool_search.R", list(mode = mode, jobs = jobs,
  target = .80, rasch_budget = 600L, glmm_budget = 1080L, confirmation = 512L),
  c("glmmplanr", "mlpwr", "DiceKriging", "randtoolbox", "rngWELL", "rgenoud"))
write.csv(jobs, file.path(out, "jobs.csv"), row.names = FALSE)
run_job <- function(k) {
  job <- jobs[k, ]; path <- file.path(out, sprintf("search-%03d.rds", job$job))
  if (file.exists(path)) return(readRDS(path)$selection)
  is_rasch <- job$task == "rasch"
  if (is_rasch) {
    grid <- rasch_cases[rasch_cases$family == job$family, ]
    grid <- grid[order(grid$n_subj), ]; rownames(grid) <- NULL
    grid$i <- seq_len(nrow(grid)); grid$j <- 1L
    grid$cost <- grid$n_subj
    runners <- lapply(seq_len(nrow(grid)), function(i) rasch_runner(grid[i, ]))
    costfun <- function(i) grid$cost[as.integer(i)]
    boundaries <- list(i = c(1, nrow(grid)))
  } else {
    spec <- experiment_spec(job$family)
    grid <- expand.grid(n_subj = c(30L, 60L, 120L), n_item = c(30L, 60L, 150L))
    grid$i <- match(grid$n_subj, c(30L, 60L, 120L))
    grid$j <- match(grid$n_item, c(30L, 60L, 150L))
    ns_unit <- (grid$n_subj - 30) / 90; ni_unit <- (grid$n_item - 30) / 120
    grid$cost <- if (job$preference == "fewer_subj") ns_unit + .02 * ni_unit else ni_unit + .02 * ns_unit
    runners <- lapply(seq_len(nrow(grid)), function(i) experiment_runner(spec, grid$n_subj[i], grid$n_item[i]))
    costfun <- function(i, j) grid$cost[match(paste(i, j), paste(grid$i, grid$j))]
    boundaries <- list(i = c(1, 3), j = c(1, 3))
  }
  target_budget <- if (mode == "smoke") if (is_rasch) 100L else 180L else if (is_rasch) 600L else 1080L
  seed <- 9470000L + 1000L * job$job
  streams <- glmmplanr:::make_lecuyer_streams(seed, target_budget + 200L)
  observations <- list(); counter <- 0L
  evaluate <- function(i, j = 1L) {
    i <- as.integer(i); j <- as.integer(j)
    pos <- match(paste(i, j), paste(grid$i, grid$j))
    stopifnot(!is.na(pos))
    counter <<- counter + 1L
    stopifnot(counter <= length(streams))
    x <- glmmplanr:::with_rng_stream(streams[[counter]], runners[[pos]]())
    observations[[counter]] <<- data.frame(pos = pos, x)
    isTRUE(x$success)
  }
  simfun <- if (is_rasch) function(i) evaluate(i) else function(i, j) evaluate(i, j)
  # Verify the shared success indicator against a public planning call before
  # either search starts. Its data stream is separate from search/confirmation.
  audit_seed <- glmmplanr:::make_lecuyer_streams(seed + 101L, 1L)[[1L]]
  if (!is_rasch) {
    direct <- glmmplanr:::with_rng_stream(audit_seed, runners[[1L]]())
    public <- evaluate_design_point(spec, grid$n_subj[1], grid$n_item[1], nsim = 1L,
                                    stream_seeds = list(audit_seed))
    stopifnot(direct$valid == (public$n_valid == 1L), direct$success == (public$k_joint == 1L))
  }
  warnings <- character()
  set.seed(seed + 201L)
  ml_time <- system.time(ml <- tryCatch(withCallingHandlers(mlpwr::find.design(
    simfun = simfun, boundaries = boundaries, power = .8, costfun = costfun,
    evaluations = target_budget, setsize = if (mode == "smoke") 20L else 60L,
    n.startsets = if (is_rasch) 4L else 3L, silent = TRUE), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = identity))[["elapsed"]]
  obs <- do.call(rbind, observations)
  stopifnot(counter == nrow(obs), counter > 0L)
  # mlpwr rounds some batches upward. Charge its actual number of simulations
  # to BOTH methods; distribute the grid remainder in a fixed, outcome-blind order.
  counts <- rep(counter %/% nrow(grid), nrow(grid))
  counts[seq_len(counter %% nrow(grid))] <- counts[seq_len(counter %% nrow(grid))] + 1L
  stopifnot(sum(counts) == counter, all(counts > 0L))
  grid_time <- system.time(gs <- lapply(seq_len(nrow(grid)), function(g) {
    x <- grid[g, ]
    if (is_rasch) {
      plan <- plan_rasch_fixed_bank(bank_for(x), x$n_subj,
        c(difficulty = .35, difficulty_max_abs_error = .60), person_sd = x$person_sd,
        items_per_subject = x$m, item_seed = 37L, nsim = counts[g], seed = seed + 300L + g, engine = "mml")
      z <- subset(as.data.frame(plan), target == "joint")
      return(data.frame(pos = g, nsim = counts[g], n_met = z$n_met, n_valid = z$n_evaluable))
    }
    z <- evaluate_design_point(spec, x$n_subj, x$n_item, nsim = counts[g], stream_seeds = seed + 300L + g)
    data.frame(pos = g, nsim = counts[g], n_met = z$k_joint, n_valid = z$n_valid)
  }))[["elapsed"]]
  gs <- do.call(rbind, gs); gs$attainment <- gs$n_met / gs$nsim
  stopifnot(sum(gs$nsim) == counter)
  # Both searches propose for the SAME point target .80. Only the independent
  # 512-run Wilson check can turn a proposal into a confirmed recommendation.
  eligible <- which(gs$attainment >= .80)
  grid_pos <- if (length(eligible)) eligible[order(grid$cost[eligible], grid$n_subj[eligible],
                                                 grid$n_item[eligible])[1L]] else NA_integer_
  ml_pos <- NA_integer_
  if (!inherits(ml, "error") && length(ml$final$design)) {
    idx <- as.numeric(ml$final$design)
    if (all(is.finite(idx))) ml_pos <- match(paste(as.integer(idx[1]), if (is_rasch) 1L else as.integer(idx[2])),
                                             paste(grid$i, grid$j))
  }
  selected <- do.call(rbind, lapply(c("glmmplanr_grid", "mlpwr"), function(method) {
    pos <- if (method == "mlpwr") ml_pos else grid_pos
    data.frame(job, method, target_budget, actual_evaluations = counter,
      case = if (is_rasch && !is.na(pos)) grid$case[pos] else NA_integer_,
      n_subj = if (!is.na(pos)) grid$n_subj[pos] else NA_integer_,
      n_item = if (!is.na(pos)) grid$n_item[pos] else NA_integer_,
      cost = if (!is.na(pos)) grid$cost[pos] else NA_real_,
      elapsed_s = if (method == "mlpwr") ml_time else grid_time,
      error = if (method == "mlpwr" && inherits(ml, "error")) conditionMessage(ml) else "",
      warning = if (method == "mlpwr") paste(unique(warnings), collapse = " | ") else "")
  }))
  atomic_rds(list(selection = selected, grid = grid, grid_results = gs, mlpwr_observations = obs,
    mlpwr_final = if (!inherits(ml, "error")) ml$final else NULL), path)
  cat(sprintf("%s %s %s run %d: %d evaluations EACH, %.1fs grid / %.1fs mlpwr\n",
    job$task, job$family, job$preference, job$run, counter, grid_time, ml_time))
  selected
}
results <- parallel::mclapply(seq_len(nrow(jobs)), run_job, mc.cores = workers, mc.set.seed = FALSE)
stopifnot(all(vapply(results, is.data.frame, logical(1))))
write.csv(do.call(rbind, results), file.path(out, "search.csv"), row.names = FALSE)
