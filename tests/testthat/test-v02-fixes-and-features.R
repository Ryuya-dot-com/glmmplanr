test_that("Fix 1: Non-converged/failed fits never inflate power to 100% and log error reasons", {
  des <- design_crossed(list(
    exposure = c("incidental", "intentional"),
    spacing  = c("massed", "distributed")
  ))
  dgp <- calibrate_dgp(
    des,
    cell_probs = c(
      incidental_massed       = 0.50,
      incidental_distributed  = 0.60,
      intentional_massed      = 0.65,
      intentional_distributed = 0.80
    ),
    chance = 1 / 3,
    subj_prob_range = c(0.40, 0.85),
    item_prob_range = c(0.45, 0.88)
  )

  sp_eff <- specify_target(
    dgp, effect = "exposure", criterion = "power", power_min = 0.80,
    failure_policy = "effective"
  )

  # Exact review scenario: 10 reps, 8 failed, 2 converged & significant
  est <- c(rep(NA_real_, 8), 1.50, 1.60)
  se  <- c(rep(NA_real_, 8), 0.25, 0.25)
  pest <- c(rep(NA_real_, 8), 0.18, 0.19)
  pse  <- c(rep(NA_real_, 8), 0.04, 0.04)
  conv <- c(rep(FALSE, 8), TRUE, TRUE)
  stat <- c(rep("engine_error", 8), "converged", "converged")
  errs <- c(rep("custom mafc link not supported", 8), "", "")

  row_eff <- glmmplanr:::summarize_mc_metrics(
    est = est, se = se, prob_est = pest, prob_se = pse,
    singular = rep(FALSE, 10), converged = conv,
    status = stat, error_reasons = errs, nobs = rep(384L, 10),
    spec = sp_eff, n_subj = 24L, n_item = 16L, nsim = 10L, engine_use = "glmmTMB"
  )

  expect_equal(row_eff$valid_rate, 0.20)
  expect_equal(row_eff$conditional_power, 1.00)
  expect_equal(row_eff$effective_power, 0.20)
  expect_equal(row_eff$power, 0.20)
  expect_false(row_eff$criterion_met)
  expect_equal(row_eff$decision_status, "fail")
  expect_match(row_eff$top_error_reason, "\\[8/10\\] custom mafc link not supported")

  # Even under conditional policy, min_valid_rate (0.90) blocks 2/10 valid fits from passing
  sp_cond <- specify_target(
    dgp, effect = "exposure", criterion = "power", power_min = 0.80,
    failure_policy = "conditional", min_valid_rate = 0.90
  )
  row_cond <- glmmplanr:::summarize_mc_metrics(
    est = est, se = se, prob_est = pest, prob_se = pse,
    singular = rep(FALSE, 10), converged = conv,
    status = stat, error_reasons = errs, nobs = rep(384L, 10),
    spec = sp_cond, n_subj = 24L, n_item = 16L, nsim = 10L, engine_use = "glmmTMB"
  )
  expect_equal(row_cond$power, 1.00)
  expect_false(row_cond$criterion_met)
  expect_equal(row_cond$decision_status, "fail")
})

test_that("Fix 2: Stage x Point x Replication L'Ecuyer-CMRG streams have zero collisions", {
  stride <- 20L
  all_keys <- character(0)
  for (stage_slot in 1L:5L) {
    bases <- glmmplanr:::allocate_stage_point_streams(
      seed = 20261119L, stage_slot = stage_slot, grid_ids = seq_len(stride), stride = stride
    )
    for (g in seq_len(stride)) {
      reps <- glmmplanr:::derive_rep_substreams(bases[[g]], nsim = 15L, offset = 0L)
      keys <- vapply(reps, function(s) paste(s, collapse = ","), character(1L))
      all_keys <- c(all_keys, keys)
    }
  }
  # All 5 * 20 * 15 = 1,500 RNG states must be 100% distinct (including Stage 1 rep 2 vs Stage 2 rep 1!)
  expect_equal(anyDuplicated(all_keys), 0L)
})

test_that("Fix 3: Scheduled trials and post-attrition expected trials match exact allocation", {
  # 2x2 Full Crossed: 24 subj x 16 item x 4 cells = 1,536 scheduled trials
  des_cr <- design_crossed(
    factors = list(
      A = c("a1", "a2"),
      B = c("b1", "b2")
    ),
    allocation = "crossed",
    attrition = 0.10
  )
  sched_cr <- generate_design_schedule(des_cr, n_subj = 24, n_item = 16)
  expect_equal(nrow(sched_cr$trials), 1536L)
  expect_equal(compute_scheduled_trials(des_cr, 24, 16), 1536L)
  expect_equal(compute_expected_trials(des_cr, 24, 16), as.integer(round(1536 * 0.90)))

  # 2x2 Latin Square: 24 subj x 16 item = 384 scheduled trials
  des_ls <- design_crossed(
    factors = list(
      A = c("a1", "a2"),
      B = c("b1", "b2")
    ),
    allocation = "latin_square",
    attrition = 0.10
  )
  sched_ls <- generate_design_schedule(des_ls, n_subj = 24, n_item = 16)
  expect_equal(nrow(sched_ls$trials), 384L)
  expect_equal(compute_scheduled_trials(des_ls, 24, 16), 384L)
  expect_equal(compute_expected_trials(des_ls, 24, 16), as.integer(round(384 * 0.90)))
})

test_that("Fix 4: C++ Laplace separates max_iter_finite and agrees with lme4", {
  des <- design_crossed(list(
    exposure = c("incidental", "intentional"),
    spacing  = c("massed", "distributed")
  ))
  dgp0 <- calibrate_dgp(
    des,
    cell_probs = c(
      incidental_massed       = 0.35,
      incidental_distributed  = 0.50,
      intentional_massed      = 0.55,
      intentional_distributed = 0.78
    ),
    chance = 0,
    subj_sd = 0.70,
    item_sd = 0.60,
    slope_sd_ratio = 0.30,
    corr_structure = 0.15
  )

  sched <- generate_design_schedule(des, n_subj = 48, n_item = 32)
  set.seed(2026)
  df0 <- simulate_dgp_data(dgp0, schedule = sched)
  X <- des$X_cells[df0$cell_id, , drop = FALSE]

  # 4a. Strict separation of max_iter_finite vs converged
  f_cap <- glmmplanr:::cpp_fit_crossed_glmm(
    y = as.integer(df0$y), X = X,
    subj_idx = df0$subj_id - 1L, Z_subj = X, n_subj = 48L,
    item_idx = df0$item_id - 1L, Z_item = X, n_item = 32L,
    testlet_idx = integer(0), Z_testlet = matrix(0, nrow(df0), 0), n_testlet = 0L,
    chance = 0.0, corr_mode = 1L, regularized = FALSE,
    max_em_iter = 2L, tol = 1e-7
  )
  expect_false(f_cap$converged)
  expect_equal(f_cap$status, "max_iter_finite")

  # 4b. Canonical logit agreement with lme4::glmer(..., nAGQ = 1)
  sp0_cpp <- specify_target(dgp0, effect = "exposure", engine = "cpp_laplace", on_singular = "flag_only")
  sp0_lme <- specify_target(dgp0, effect = "exposure", engine = "lme4", on_singular = "flag_only")
  fit_cpp <- glmmplanr:::fit_single_dataset(df0, sp0_cpp, engine_use = "cpp_laplace")
  fit_lme <- glmmplanr:::fit_single_dataset(df0, sp0_lme, engine_use = "lme4")

  expect_true(fit_cpp$converged)
  expect_equal(fit_cpp$status, "converged")
  expect_equal(fit_cpp$beta_full, fit_lme$beta_full, tolerance = 0.06)
  expect_equal(fit_cpp$se, fit_lme$se, tolerance = 0.025)
})

test_that("3AFC calibration and coverage meet the small regression fixture", {
  dgp3 <- three_choice_reference_dgp()
  sp3 <- specify_target(dgp3, effect = "exposure", engine = "bayes_laplace", on_singular = "flag_only")
  ev3 <- evaluate_design_point(sp3, n_subj = 56, n_item = 32, nsim = 35L, stream_seeds = 303L)
  expect_gt(ev3$valid_rate, 0.95)
  expect_gte(ev3$coverage, 0.82)
  expect_gte(ev3$prob_coverage, 0.82)
})

test_that("Features 1-3: Percentage-point precision, indeterminate top-up, and multi-scenario robustness increment", {
  # Morris et al. (2019) borderline decision check: 801/1000 and 799/1000 are BOTH indeterminate!
  ci_801 <- glmmplanr:::wilson_mc_ci(801L, 1000L, conf_level = 0.95)
  ci_799 <- glmmplanr:::wilson_mc_ci(799L, 1000L, conf_level = 0.95)
  ci_840 <- glmmplanr:::wilson_mc_ci(840L, 1000L, conf_level = 0.95)
  ci_750 <- glmmplanr:::wilson_mc_ci(750L, 1000L, conf_level = 0.95)
  expect_equal(glmmplanr:::classify_mc_threshold(ci_801["lo"], ci_801["hi"], 0.80), "indeterminate")
  expect_equal(glmmplanr:::classify_mc_threshold(ci_799["lo"], ci_799["hi"], 0.80), "indeterminate")
  expect_equal(glmmplanr:::classify_mc_threshold(ci_840["lo"], ci_840["hi"], 0.80), "pass")
  expect_equal(glmmplanr:::classify_mc_threshold(ci_750["lo"], ci_750["hi"], 0.80), "fail")

  # Feature 1: Percentage-point (pp) accuracy difference design & AIPE precision
  des <- design_crossed(
    factors = list(method = c("control", "treatment")),
    allocation = "latin_square",
    contrast = "sum"
  )
  c_probs <- cell_probs_from_diffs(des, baseline_prob = 0.68, effects_pp = c(method = 10))
  expect_equal(unname(c_probs), c(0.63, 0.73), tolerance = 1e-10)

  dgp_base <- calibrate_dgp(
    des,
    cell_probs = c_probs,
    chance = 0.25,
    subj_prob_range = c(0.45, 0.88),
    item_prob_range = c(0.48, 0.88),
    slope_sd_ratio = 0.30
  )

  spec_pp <- specify_target(
    dgp_base,
    effect = "method",
    scale = "pp",
    prob_estimand = "marginal",
    criterion = "both",
    power_min = 0.80,
    ci_width = 14.0,
    prob_boot = 9L,
    engine = "bayes_laplace",
    on_singular = "flag_only"
  )
  expect_equal(spec_pp$true_prob_diff_pp, 10.0, tolerance = 1e-6)

  skip_on_cran()

  # Feature 2 + 3: Multi-scenario robust search (6pp, 10pp, 12pp) + pre-specified top-up
  scens <- scenario_dgps(
    dgp_base,
    effect = "method",
    effect_diffs_pp = c(diff_6pp = 6, diff_10pp = 10, diff_12pp = 12)
  )

  res_rob <- optimize_sample_size(
    spec_pp,
    n_subj_range = c(24, 64),
    n_item_range = c(16, 40),
    prefer = "balanced",
    scenarios = scens,
    standard_scenario = "diff_10pp",
    method = "surrogate_gp",
    n_eval_init = 6L,
    n_eval_adapt = 2L,
    nsim_init = 12L,
    nsim_final = 20L,
    nsim_topup = 15L,
    n_verify = 3L,
    nsim_sensitivity = 12L,
    workers = 2L,
    seed = 20261119L,
    verbose = FALSE
  )

  expect_s3_class(res_rob, "glmmplanr_result")
  expect_false(is.null(res_rob$robust))
  if (res_rob$robust$feasible_robust) {
    expect_gte(res_rob$robust$delta_trials, 0L)
  } else {
    expect_true(is.na(res_rob$robust$delta_trials))
  }
  expect_true(res_rob$optimal$decision_status %in% c("pass", "indeterminate", "fail"))
  expect_true(is.finite(res_rob$optimal$prob_ci_width_mean_pp))

  # Verify preference switching ("fewer_items", "fewer_subj", "fewer_trials")
  res_fi <- switch_preference(res_rob, prefer = "fewer_items")
  res_fs <- switch_preference(res_rob, prefer = "fewer_subj")
  res_ft <- switch_preference(res_rob, prefer = "fewer_trials")
  expect_equal(res_fi$prefer, "fewer_items")
  expect_equal(res_fs$prefer, "fewer_subj")
  expect_equal(res_ft$prefer, "fewer_trials")
  if (res_rob$feasible) {
    expect_lte(res_fi$optimal$n_item, res_fs$optimal$n_item)
    expect_lte(res_fs$optimal$n_subj, res_fi$optimal$n_subj)
  }

  # Check plots on pp scale, gt table with robust section, and EN/JA reports
  p_pp <- plot(res_rob, type = "contour", metric = "prob_ci_width_pp")
  expect_s3_class(p_pp, "ggplot")
  gt_rob <- as_gt(res_rob)
  expect_s3_class(gt_rob, "gt_tbl")
  rep_ja <- report_plan(res_rob, lang = "ja", print_output = FALSE)
  expect_match(rep_ja, "パーセントポイント")
  expect_match(rep_ja, "積み増し")
})

test_that("Speedup architecture preserves bit-for-bit C++ accuracy, worker invariance, and unbiased futility curtailment", {
  des <- design_crossed(
    factors = list(
      exposure = c("incidental", "intentional"),
      spacing  = c("massed", "distributed")
    ),
    allocation = "latin_square",
    attrition = 0.08
  )
  dgp <- calibrate_dgp(
    des,
    cell_probs = c(
      incidental_massed       = 0.48,
      incidental_distributed  = 0.58,
      intentional_massed      = 0.64,
      intentional_distributed = 0.82
    ),
    chance = 1 / 3,
    subj_prob_range = c(0.42, 0.84),
    item_prob_range = c(0.45, 0.86),
    slope_sd_ratio = 0.30,
    corr_structure = 0.15
  )
  spec <- specify_target(
    dgp,
    effect = "exposure",
    criterion = "power",
    power_min = 0.85,
    engine = "bayes_laplace",
    stage1_re = "uncorrelated",
    on_singular = "flag_and_compare"
  )

  # 1. Bit-for-bit numerical identity between zero-alloc cached path and full data.frame path
  sched <- generate_design_schedule(des, n_subj = 32, n_item = 24)
  cache <- glmmplanr:::prepare_cpp_schedule_cache(spec, sched)
  st_list <- glmmplanr:::derive_rep_substreams(
    glmmplanr:::allocate_stage_point_streams(999L, 1L, 1L, 1L)[[1L]], nsim = 5L
  )
  for (b in seq_len(5L)) {
    fast_r <- glmmplanr:::with_rng_stream(
      st_list[[b]],
      glmmplanr:::execute_single_replication(
        spec, sched, cpp_cache = cache, engine_use = "bayes_laplace",
        re_primary = "maximal", fast_stage = FALSE, do_compare = TRUE
      )
    )
    slow_r <- glmmplanr:::with_rng_stream(
      st_list[[b]],
      glmmplanr:::execute_single_replication(
        spec, sched, cpp_cache = NULL, engine_use = "bayes_laplace",
        re_primary = "maximal", fast_stage = FALSE, do_compare = TRUE
      )
    )
    expect_equal(fast_r$primary$est, slow_r$primary$est, tolerance = 1e-14)
    expect_equal(fast_r$primary$se, slow_r$primary$se, tolerance = 1e-14)
    expect_equal(fast_r$primary$prob_est, slow_r$primary$prob_est, tolerance = 1e-14)
    expect_equal(fast_r$nobs, slow_r$nobs)
  }

  # 2. Group-sequential futility curtailment stops tiny underpowered grid points early
  #    while matching early_stop = FALSE bit-for-bit on viable points
  skip_on_cran()
  res_es_on <- optimize_sample_size(
    spec,
    n_subj_range = c(12, 20, 72),
    n_item_range = c(12, 20, 44),
    method = "grid",
    nsim_final = 200L,
    early_stop = TRUE,
    nsim_sensitivity = 0L,
    workers = 2L,
    seed = 20261119L,
    verbose = FALSE
  )
  res_es_off <- optimize_sample_size(
    spec,
    n_subj_range = c(12, 20, 72),
    n_item_range = c(12, 20, 44),
    method = "grid",
    nsim_final = 200L,
    early_stop = FALSE,
    nsim_sensitivity = 0L,
    workers = 1L,
    seed = 20261119L,
    verbose = FALSE
  )
  expect_true(any(res_es_on$verified$curtailed))
  expect_false(res_es_on$optimal$curtailed)
  expect_equal(res_es_on$optimal$n_subj, res_es_off$optimal$n_subj)
  expect_equal(res_es_on$optimal$n_item, res_es_off$optimal$n_item)
  expect_equal(res_es_on$optimal$power, res_es_off$optimal$power, tolerance = 1e-14)
})
