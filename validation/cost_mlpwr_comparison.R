# Pre-specified cost-constrained 2D design search comparison: glmmplanr vs. mlpwr
# Aligns both solvers on the EXACT same 25-point discrete candidate lattice,
# identical 1,200-evaluation total budget, and 5 independent search seeds per regime,
# separating three distinct decision rules:
#   1) point_80: matched 80% point-estimate target (glmmplanr Stage-2 power >= 0.80 vs mlpwr power = 0.80)
#   2) high_target_86: mlpwr raised point target power = 0.860 (exact minimum out of 200 whose
#      95% Wilson lower bound is >= 0.80: 172/200 -> lo = 0.8050; 171/200 -> lo = 0.7995 < 0.80)
#   3) wilson_lo_80: explicit 200-rep 95% Wilson lower bound >= 0.80 confirmation
#      (glmmplanr Stage-2 Wilson gate vs mlpwr 1,000-eval search + 200-rep independent Wilson check)
# Run with: sh validation/run-isolated.sh cost_mlpwr_comparison
library(glmmplanr)
stopifnot(requireNamespace("mlpwr", quietly = TRUE))
RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")

select_min_cost <- function(pool, cost_fn) {
  if (is.null(pool) || nrow(pool) == 0L) return(NULL)
  costs <- cost_fn(pool$n_subj, pool$n_item)
  pool[which.min(costs), , drop = FALSE]
}

# Fixed population: 2-list Latin square, 68% vs 80% marginal accuracy, chance = 0
des <- design_crossed(list(group = c("A", "B")), allocation = "latin_square", contrast = "treatment")
dgp <- calibrate_dgp(
  des, cell_probs = c(A = 0.68, B = 0.80), chance = 0,
  subj_sd = 0.8, item_sd = 0.6, slope_sd_ratio = 0, corr_structure = 0,
  control = calibrate_control(seed = 8601L)
)
sp <- specify_target(
  dgp, effect = "group", scale = "logit", criterion = "power", power_min = 0.80,
  engine = "bayes_laplace", analysis_re = "intercept_only",
  control = target_control(on_singular = "flag_only")
)

# Shared 5 x 5 = 25-point discrete candidate lattice for ALL solvers
subj_grid <- c(16L, 24L, 32L, 40L, 48L)
item_grid <- c(8L, 12L, 16L, 20L, 24L)
snap_subj <- function(idx_s) subj_grid[min(length(subj_grid), max(1L, as.integer(round(idx_s))))]
snap_item <- function(idx_i) item_grid[min(length(item_grid), max(1L, as.integer(round(idx_i))))]

regimes <- list(
  expensive_items = list(vec = c(subj = 10, item = 80), fn = function(n_subj, n_item) 10 * n_subj + 80 * n_item),
  expensive_subj  = list(vec = c(subj = 80, item = 10), fn = function(n_subj, n_item) 80 * n_subj + 10 * n_item)
)

# 1. Exhaustive high-precision reference grid (nsim = 400 per point, 10,000 sims total)
ref_plan <- optimize_sample_size(
  sp, n_subj_range = subj_grid, n_item_range = item_grid,
  method = "grid", nsim_final = 400L,
  workers = if (.Platform$OS.type == "windows") 1L else 4L, seed = 861000L, verbose = FALSE,
  control = optimize_control(nsim_topup = 0L, nsim_sensitivity = 0L)
)

# 2. Separate independent out-of-sample evaluation grid (nsim = 400 per point, fresh seed 869000)
eval_plan <- optimize_sample_size(
  sp, n_subj_range = subj_grid, n_item_range = item_grid,
  method = "grid", nsim_final = 400L,
  workers = if (.Platform$OS.type == "windows") 1L else 4L, seed = 869000L, verbose = FALSE,
  control = optimize_control(nsim_topup = 0L, nsim_sensitivity = 0L)
)
eval_tbl <- eval_plan$verified

lookup_eval <- function(ns, ni) {
  hit <- eval_tbl[eval_tbl$n_subj == ns & eval_tbl$n_item == ni, , drop = FALSE]
  stopifnot(nrow(hit) == 1L)
  hit
}

# Helper for mlpwr binary success simulation mapped to the exact same 25-point lattice
simfun_mlpwr <- function(idx_s, idx_i) {
  ns <- snap_subj(idx_s)
  ni <- snap_item(idx_i)
  sched <- generate_design_schedule(des, ns, ni)
  tr <- simulate_dgp_data(dgp, sched)
  f <- glmmplanr:::fit_single_dataset(tr, sp, engine_use = "bayes_laplace", re_structure = "intercept_only")
  if (!isTRUE(f$converged) || !is.finite(f$est) || !is.finite(f$se) || f$se <= 0) return(FALSE)
  (2 * stats::pnorm(-abs(f$est / f$se))) < 0.05
}

n_seeds <- 5L
run_rows <- list()
row_i <- 0L

for (r_idx in seq_along(regimes)) {
  r_nm <- names(regimes)[r_idx]
  rg <- regimes[[r_nm]]
  pref_rule <- if (rg$vec[["item"]] > rg$vec[["subj"]]) "fewer_items" else "fewer_subj"

  ref_point  <- select_min_cost(ref_plan$verified[ref_plan$verified$power >= 0.80, , drop = FALSE], rg$fn)
  ref_wilson <- select_min_cost(ref_plan$verified[ref_plan$verified$decision_status == "pass", , drop = FALSE], rg$fn)

  add_record <- function(solver, criterion, seed_id, ns, ni, internal_wilson_pass = NA) {
    ev <- lookup_eval(ns, ni)
    row_i <<- row_i + 1L
    run_rows[[row_i]] <<- data.frame(
      regime = r_nm,
      criterion = criterion,
      solver = solver,
      seed_id = seed_id,
      n_subj = ns,
      n_item = ni,
      cost = rg$fn(ns, ni),
      internal_wilson_pass = internal_wilson_pass,
      eval_nsim = ev$nsim,
      eval_valid = ev$n_valid,
      eval_power = ev$power,
      eval_mc_lo = ev$power_mc_lo,
      eval_mc_hi = ev$power_mc_hi,
      eval_point_ok = ev$power >= 0.80,
      eval_wilson_pass = ev$decision_status == "pass",
      ref_point_cost = ref_point$n_subj * rg$vec[["subj"]] + ref_point$n_item * rg$vec[["item"]],
      ref_wilson_cost = ref_wilson$n_subj * rg$vec[["subj"]] + ref_wilson$n_item * rg$vec[["item"]],
      stringsAsFactors = FALSE
    )
  }

  add_record("reference_grid_400", "point_80", 0L, ref_point$n_subj, ref_point$n_item, NA)
  add_record("reference_grid_400", "wilson_lo_80", 0L, ref_wilson$n_subj, ref_wilson$n_item, TRUE)

  for (s_id in seq_len(n_seeds)) {
    # A. glmmplanr surrogate_gp: 10*40 + 4*200 = 1,200 total evaluations on 25-point lattice
    gp_res <- optimize_sample_size(
      sp, n_subj_range = subj_grid, n_item_range = item_grid,
      prefer = pref_rule, method = "surrogate_gp",
      nsim_final = 200L,
      workers = 1L, seed = 862000L + 100L * r_idx + s_id, verbose = FALSE,
      control = optimize_control(
        n_eval_init = 6L, n_eval_adapt = 4L, nsim_init = 40L,
        n_verify = 4L, nsim_topup = 0L, nsim_sensitivity = 0L
      )
    )
    gp_pt_pool <- gp_res$verified[gp_res$verified$power >= 0.80, , drop = FALSE]
    gp_pt_opt <- if (nrow(gp_pt_pool) > 0L) {
      select_min_cost(gp_pt_pool, rg$fn)
    } else {
      gp_res$verified[which.max(gp_res$verified$power), , drop = FALSE]
    }
    gp_wl_pool <- gp_res$verified[gp_res$verified$decision_status == "pass", , drop = FALSE]
    gp_wl_opt <- if (nrow(gp_wl_pool) > 0L) {
      select_min_cost(gp_wl_pool, rg$fn)
    } else {
      gp_res$verified[which.max(gp_res$verified$power_mc_lo), , drop = FALSE]
    }
    add_record("glmmplanr_1200", "point_80", s_id, gp_pt_opt$n_subj, gp_pt_opt$n_item,
               gp_pt_opt$decision_status == "pass")
    add_record("glmmplanr_1200", "wilson_lo_80", s_id, gp_wl_opt$n_subj, gp_wl_opt$n_item,
               gp_wl_opt$decision_status == "pass")

    # B. mlpwr::find.design on EXACT SAME 25-point lattice, 1,200 evaluations, target power = 0.80
    set.seed(863000L + 100L * r_idx + s_id)
    ml_pt <- suppressWarnings(mlpwr::find.design(
      simfun = simfun_mlpwr,
      boundaries = list(idx_s = c(1, 5), idx_i = c(1, 5)),
      power = 0.80,
      evaluations = 1200,
      costfun = function(idx_s, idx_i) rg$fn(snap_subj(idx_s), snap_item(idx_i)),
      silent = TRUE
    ))
    ml_pt_d <- as.numeric(ml_pt$final$design)
    add_record("mlpwr_1200", "point_80", s_id, snap_subj(ml_pt_d[1L]), snap_item(ml_pt_d[2L]), NA)

    # C. mlpwr::find.design with raised point target power = 0.860 (1,200 evals; 172/200 -> Wilson lo = 0.8050)
    #    Note: raising the surrogate point target is NOT equivalent to applying a Wilson lower-bound gate.
    set.seed(864000L + 100L * r_idx + s_id)
    ml_86 <- suppressWarnings(mlpwr::find.design(
      simfun = simfun_mlpwr,
      boundaries = list(idx_s = c(1, 5), idx_i = c(1, 5)),
      power = 0.860,
      evaluations = 1200,
      costfun = function(idx_s, idx_i) rg$fn(snap_subj(idx_s), snap_item(idx_i)),
      silent = TRUE
    ))
    ml_86_d <- as.numeric(ml_86$final$design)
    add_record("mlpwr_target86_1200", "high_target_86", s_id, snap_subj(ml_86_d[1L]), snap_item(ml_86_d[2L]), NA)

    # D. mlpwr 1,000-eval search (power = 0.80) + 200-rep independent Wilson lower-bound check = 1,200 budget
    set.seed(865000L + 100L * r_idx + s_id)
    ml_1k <- suppressWarnings(mlpwr::find.design(
      simfun = simfun_mlpwr,
      boundaries = list(idx_s = c(1, 5), idx_i = c(1, 5)),
      power = 0.80,
      evaluations = 1000,
      costfun = function(idx_s, idx_i) rg$fn(snap_subj(idx_s), snap_item(idx_i)),
      silent = TRUE
    ))
    ml_1k_d <- as.numeric(ml_1k$final$design)
    ml_ns <- snap_subj(ml_1k_d[1L])
    ml_ni <- snap_item(ml_1k_d[2L])
    ml_chk <- evaluate_design_point(
      sp, n_subj = ml_ns, n_item = ml_ni, nsim = 200L,
      stream_seeds = 866000L + 100L * r_idx + s_id
    )
    add_record("mlpwr_1000_wilson200", "wilson_lo_80", s_id, ml_ns, ml_ni,
               ml_chk$decision_status == "pass")
  }
}

runs_df <- do.call(rbind, run_rows)
rownames(runs_df) <- NULL
write.csv(runs_df, "validation/cost_mlpwr_runs.csv", row.names = FALSE)

agg_key <- paste(runs_df$regime, runs_df$criterion, runs_df$solver, sep = "|")
agg_list <- lapply(split(runs_df, agg_key), function(d) {
  n_r <- nrow(d)
  designs_str <- paste(sprintf("%dx%d", d$n_subj, d$n_item), collapse = ", ")
  int_pass <- if (all(is.na(d$internal_wilson_pass))) NA_integer_ else sum(d$internal_wilson_pass, na.rm = TRUE)
  data.frame(
    regime = d$regime[1L],
    criterion = d$criterion[1L],
    solver = d$solver[1L],
    n_seeds = n_r,
    n_internal_wilson_pass = int_pass,
    mean_cost = mean(d$cost),
    min_cost = min(d$cost),
    max_cost = max(d$cost),
    mean_eval_power = mean(d$eval_power),
    min_eval_power = min(d$eval_power),
    n_point_ok = sum(d$eval_point_ok),
    n_wilson_pass = sum(d$eval_wilson_pass),
    designs = designs_str,
    stringsAsFactors = FALSE
  )
})
sum_df <- do.call(rbind, agg_list)
sum_df <- sum_df[order(sum_df$regime, sum_df$criterion, sum_df$solver), ]
rownames(sum_df) <- NULL
write.csv(sum_df, "validation/cost_mlpwr_summary.csv", row.names = FALSE)
print(sum_df, row.names = FALSE, digits = 4)
