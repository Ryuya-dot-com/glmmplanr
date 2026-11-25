# Pre-specified 3AFC boundary/necessity validation and high-N strict conditional-null calibration
# Part A: Boundary check on 5 smaller adjacent candidates (< 2,304 trials) under planned effect (192 reps each)
# Part B: Multi-seed optimize_sample_size() boundary hit-rate check across 5 independent search seeds
# Part C: High-N (512 reps each) strict conditional-null (beta_exposure = 0) Type I error & coverage confirmation
#         with pre-specified calibration gates (Type I error exact 95% CI upper <= 7.5%, coverage lower >= 92.5%)
# Run with: sh validation/run-isolated.sh three_choice_boundary_and_null
library(glmmplanr)

mc_cores <- if (.Platform$OS.type == "windows") 1L else min(4L, max(1L, parallel::detectCores() - 1L))
z <- stats::qnorm(0.975)

des <- design_crossed(list(
  exposure = c("incidental", "intentional"),
  spacing = c("massed", "distributed")
))
dgp_planned <- calibrate_dgp(
  des,
  cell_probs = c(
    incidental_massed       = 0.50,
    incidental_distributed  = 0.62,
    intentional_massed      = 0.68,
    intentional_distributed = 0.84
  ),
  chance = 1 / 3,
  subj_prob_range = c(0.42, 0.84),
  item_prob_range = c(0.45, 0.86),
  slope_sd_ratio = 0.30,
  corr_structure = 0.15
)
sp_planned <- specify_target(
  dgp_planned, effect = "exposure", scale = "logit",
  criterion = "ci_width", ci_width = 1.00,
  engine = "bayes_laplace",
  control = target_control(ci_assurance = 0.80, on_singular = "flag_only")
)
sp_strict_null <- glmmplanr:::null_target_spec(sp_planned)
stopifnot(abs(sp_strict_null$true_beta) < 1e-12)

run_trials <- function(sp, effect_label, candidates, n_rep, seed_base) {
  out_list <- vector("list", nrow(candidates))
  for (c_i in seq_len(nrow(candidates))) {
    ns <- candidates$persons[c_i]
    ni <- candidates$items[c_i]
    sched <- generate_design_schedule(des, n_subj = ns, n_item = ni)
    worker_fn <- function(tr_i) {
      sd_val <- seed_base + 1000L * c_i + tr_i
      set.seed(sd_val, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion", sample.kind = "Rejection")
      dat <- simulate_dgp_data(sp$dgp, schedule = sched)
      f <- glmmplanr:::fit_single_dataset(dat, sp, engine_use = "bayes_laplace")
      valid <- isTRUE(f$converged) && is.finite(f$est) && is.finite(f$se) && f$se > 0
      data.frame(
        effect = effect_label,
        persons = ns,
        items = ni,
        trials = ns * ni,
        trial = tr_i,
        seed = sd_val,
        truth = sp$true_beta,
        valid = valid,
        est = if (valid) f$est else NA_real_,
        se = if (valid) f$se else NA_real_,
        stringsAsFactors = FALSE
      )
    }
    reps <- if (mc_cores > 1L) {
      parallel::mclapply(seq_len(n_rep), worker_fn, mc.cores = mc_cores)
    } else {
      lapply(seq_len(n_rep), worker_fn)
    }
    out_list[[c_i]] <- do.call(rbind, reps)
  }
  do.call(rbind, out_list)
}

judge_lower <- function(ci, floor_val) {
  if (anyNA(ci)) return("indeterminate")
  if (ci[1L] >= floor_val) return("pass")
  if (ci[2L] < floor_val) return("below_floor")
  "indeterminate"
}

judge_upper <- function(ci, ceil_val) {
  if (anyNA(ci)) return("indeterminate")
  if (ci[2L] <= ceil_val) return("pass")
  if (ci[1L] > ceil_val) return("above_ceiling")
  "indeterminate"
}

summarize_trials <- function(fits, cov_floor = 0.90, wid_floor = 0.80, type1_ceiling = NA_real_) {
  groups <- split(fits, list(fits$effect, fits$persons, fits$items), drop = TRUE)
  res <- do.call(rbind, lapply(groups, function(d) {
    n <- nrow(d)
    nv <- sum(d$valid)
    err <- d$est[d$valid] - d$truth[d$valid]
    se <- d$se[d$valid]
    width <- 2 * z * se
    covered <- sum(abs(err) <= z * se)
    narrow <- sum(width <= 1.0)
    detected <- sum(abs(d$est[d$valid]) > z * se)

    cov_ci <- as.numeric(stats::binom.test(covered, nv)$conf.int)
    wid_ci <- as.numeric(stats::binom.test(narrow, n)$conf.int)
    det_ci <- as.numeric(stats::binom.test(detected, n)$conf.int)

    cov_st <- judge_lower(cov_ci, cov_floor)
    wid_st <- judge_lower(wid_ci, wid_floor)
    t1_st <- if (is.finite(type1_ceiling)) judge_upper(det_ci, type1_ceiling) else NA_character_

    data.frame(
      effect = d$effect[1L],
      persons = d$persons[1L],
      items = d$items[1L],
      trials = d$trials[1L],
      attempted = n,
      valid = nv,
      valid_rate = nv / n,
      covered = covered,
      coverage = covered / nv,
      coverage_lo = cov_ci[1L],
      coverage_hi = cov_ci[2L],
      coverage_status = cov_st,
      narrow = narrow,
      assurance = narrow / n,
      width_lo = wid_ci[1L],
      width_hi = wid_ci[2L],
      width_status = wid_st,
      detected = detected,
      detect_rate = detected / n,
      detect_lo = det_ci[1L],
      detect_hi = det_ci[2L],
      type1_status = t1_st,
      bias = mean(err),
      empirical_sd = stats::sd(err),
      mean_se = mean(se),
      mean_width = mean(width),
      stringsAsFactors = FALSE
    )
  }))
  res <- res[order(res$trials, res$persons), ]
  rownames(res) <- NULL
  res
}

# --- Part A: Boundary check on 6 smaller adjacent candidates (planned effect, n_rep = 192) ---
# Includes 72x24 (1,728 trials) so every design selected by the Part B optimizer has
# an independent 192-rep out-of-sample confirmation.
cand_boundary <- data.frame(
  cand_id = 1:6,
  persons = c(32L, 48L, 48L, 64L, 48L, 72L),
  items   = c(32L, 24L, 32L, 24L, 40L, 24L)
)
fits_boundary <- run_trials(sp_planned, "planned_boundary", cand_boundary, n_rep = 192L, seed_base = 7610000L)
sum_boundary <- summarize_trials(fits_boundary, cov_floor = 0.90, wid_floor = 0.80)
write.csv(fits_boundary, "validation/three_choice_boundary_fits.csv", row.names = FALSE)
write.csv(sum_boundary, "validation/three_choice_boundary_summary.csv", row.names = FALSE)

cat("\n=== Part A: 3AFC Boundary Check (planned effect, n_rep = 192) ===\n")
print(sum_boundary[, c("persons", "items", "trials", "valid", "covered", "coverage_lo", "narrow", "assurance", "width_lo", "width_status", "mean_width")],
      row.names = FALSE, digits = 4)

# --- Part B: Optimizer boundary search across 5 independent seeds ---
# Lattice spans both failing (32x24, 32x32, 48x24) and passing (48x32..72x48) designs.
subj_lat <- c(32L, 48L, 64L, 72L)
item_lat <- c(24L, 32L, 40L, 48L)
opt_rows <- lapply(1:5, function(s_id) {
  plan <- optimize_sample_size(
    sp_planned,
    n_subj_range = subj_lat,
    n_item_range = item_lat,
    prefer = "fewer_trials",
    method = "surrogate_gp",
    nsim_final = 96L,
    workers = mc_cores,
    seed = 7650000L + s_id,
    verbose = FALSE,
    control = optimize_control(
      n_eval_init = 8L,
      n_eval_adapt = 4L,
      nsim_init = 24L,
      n_verify = 4L,
      nsim_topup = 0L,
      nsim_sensitivity = 0L
    )
  )
  opt <- plan$optimal
  data.frame(
    seed_id = s_id,
    seed = 7650000L + s_id,
    feasible = isTRUE(plan$feasible),
    n_subj = opt$n_subj,
    n_item = opt$n_item,
    trials = opt$n_trials_scheduled,
    stage2_aipe = opt$aipe_prob,
    stage2_aipe_lo = opt$aipe_mc_lo,
    stage2_mean_width = opt$logit_ci_width_mean,
    decision = opt$decision_status,
    n_pareto = nrow(plan$verified_pareto),
    pareto_designs = if (nrow(plan$verified_pareto) > 0L) {
      paste(sprintf("%dx%d", plan$verified_pareto$n_subj, plan$verified_pareto$n_item), collapse = ", ")
    } else "",
    stringsAsFactors = FALSE
  )
})
opt_df <- do.call(rbind, opt_rows)
write.csv(opt_df, "validation/three_choice_optimizer_boundary.csv", row.names = FALSE)

cat("\n=== Part B: Multi-Seed Optimizer Boundary Search (5 seeds, prefer = 'fewer_trials') ===\n")
print(opt_df[, c("seed_id", "feasible", "n_subj", "n_item", "trials", "stage2_aipe", "stage2_aipe_lo", "decision", "pareto_designs")],
      row.names = FALSE, digits = 4)

# --- Part C: High-N Strict Conditional-Null Confirmation (beta = 0, n_rep = 512) ---
# Pre-specified calibration gates:
#   - Type I error exact 95% CI upper bound <= 7.5% (all 512 attempted trials in denominator)
#   - 95% CI coverage exact 95% CI lower bound >= 92.5% (valid fits in denominator)
#   - Width <= 1.00 logit exact 95% CI lower bound >= 80%
cand_null <- data.frame(cand_id = 1:2, persons = c(72L, 48L), items = c(32L, 48L))
fits_null512 <- run_trials(sp_strict_null, "conditional_null_512", cand_null, n_rep = 512L, seed_base = 7710000L)
sum_null512 <- summarize_trials(fits_null512, cov_floor = 0.925, wid_floor = 0.80, type1_ceiling = 0.075)
write.csv(fits_null512, "validation/three_choice_strict_null_512_fits.csv", row.names = FALSE)
write.csv(sum_null512, "validation/three_choice_strict_null_512_summary.csv", row.names = FALSE)

cat("\n=== Part C: High-N Strict Conditional-Null Confirmation (beta = 0, n_rep = 512) ===\n")
print(sum_null512[, c("persons", "items", "valid", "covered", "coverage", "coverage_lo", "coverage_status", "detected", "detect_rate", "detect_hi", "type1_status")],
      row.names = FALSE, digits = 4)
