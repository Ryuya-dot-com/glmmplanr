test_that("Universal Kriging propagates trend-parameter covariance in GP surrogates", {
  grid_df <- expand.grid(n_subj = c(24L, 48L, 72L, 96L), n_item = c(16L, 32L, 48L))
  eval_df <- data.frame(
    n_subj = c(24L, 48L, 96L, 48L),
    n_item = c(16L, 32L, 48L, 16L),
    power = c(0.35, 0.72, 0.96, 0.50),
    ci_width_mean = c(1.10, 0.72, 0.45, 0.88),
    nsim = rep(40L, 4L),
    n_valid = rep(40L, 4L)
  )

  gp_p <- fit_monotone_gp_prob(eval_df, grid_df)
  gp_w <- fit_monotone_gp_width(eval_df, grid_df)

  expect_true(all(is.finite(gp_p$sd_probit) & gp_p$sd_probit > 0))
  expect_true(all(is.finite(gp_w$sd_logw) & gp_w$sd_logw > 0))

  # Far extrapolated corner (24, 48) or (96, 16) has higher uncertainty than evaluated interior (48, 32)
  idx_eval <- which(grid_df$n_subj == 48L & grid_df$n_item == 32L)
  idx_unobs <- which(grid_df$n_subj == 96L & grid_df$n_item == 16L)
  expect_gt(gp_p$sd_probit[idx_unobs], gp_p$sd_probit[idx_eval])
  expect_gt(gp_w$sd_logw[idx_unobs], gp_w$sd_logw[idx_eval])
})

test_that("Cost and budget optimization shifts Pareto selection and identifies over-budget minimum-cost PASS", {
  des <- design_crossed(
    factors = list(exposure = c("incidental", "intentional")),
    allocation = "latin_square"
  )
  dgp <- calibrate_dgp(
    des,
    cell_probs = c(incidental = 0.52, intentional = 0.72),
    subj_prob_range = c(0.35, 0.80),
    item_prob_range = c(0.35, 0.80),
    slope_sd_ratio = 0.30
  )
  spec <- specify_target(
    dgp, effect = "exposure", criterion = "power", power_min = 0.75,
    engine = "cpp_laplace", on_singular = "flag_only"
  )

  # Cheap subjects / expensive items -> favors fewer items, more subjects
  res_exp_items <- optimize_sample_size(
    spec,
    n_subj_range = c(24L, 44L, 64L),
    n_item_range = c(12L, 24L, 36L),
    cost = c(subj = 10, item = 250, trial = 0.1, fixed = 100),
    method = "grid",
    nsim_final = 24L,
    nsim_topup = 0L,
    nsim_sensitivity = 0L,
    seed = 404L,
    verbose = FALSE
  )
  expect_identical(res_exp_items$prefer, "min_cost")
  expect_true("total_cost" %in% names(res_exp_items$optimal))
  expect_true(is.finite(res_exp_items$optimal$total_cost))

  # Re-rank same verified surface with expensive subjects / cheap items
  res_exp_subj <- switch_preference(
    res_exp_items,
    cost = c(subj = 250, item = 10, trial = 0.1, fixed = 100)
  )
  expect_lte(res_exp_items$optimal$n_item, res_exp_subj$optimal$n_item)
  expect_gte(res_exp_items$optimal$n_subj, res_exp_subj$optimal$n_subj)

  # Custom cost function equivalence
  fn_cost <- function(n_subj, n_item, n_trials_scheduled) {
    100 + 250 * n_subj + 10 * n_item + 0.1 * n_trials_scheduled
  }
  res_fn <- switch_preference(res_exp_items, cost = fn_cost)
  expect_equal(res_fn$surface$total_cost, res_exp_subj$surface$total_cost, tolerance = 1e-10)

  # Tight budget below the minimum PASS cost flags infeasibility and reports min_cost_over_budget
  min_pass_cost <- min(res_exp_subj$verified$total_cost[res_exp_subj$verified$decision_status == "pass"])
  tight_budget <- min_pass_cost - 50
  res_tight <- switch_preference(res_exp_subj, budget = tight_budget)
  expect_false(res_tight$feasible)
  expect_false(is.null(res_tight$min_cost_over_budget))
  expect_equal(res_tight$min_cost_over_budget$total_cost, min_pass_cost, tolerance = 1e-10)

  expect_s3_class(plot(res_tight, type = "contour", metric = "cost"), "ggplot")
  expect_s3_class(as_gt(res_tight), "gt_tbl")
  rep_en <- report_plan(res_tight, lang = "en", print_output = FALSE)
  rep_ja <- report_plan(res_tight, lang = "ja", print_output = FALSE)
  expect_match(rep_en, "minimum-cost verified PASS")
  expect_match(rep_ja, "最小費用で合格")
})

test_that("Marginal joint_delta SE and C++ batch bootstrap SE agree and preserve worker-invariance", {
  des <- design_crossed(list(exposure = c("incidental", "intentional")))
  dgp <- calibrate_dgp(
    des,
    cell_probs = c(incidental = 0.50, intentional = 0.68),
    subj_prob_range = c(0.32, 0.80),
    item_prob_range = c(0.35, 0.82),
    slope_sd_ratio = 0.35
  )

  sp_jd <- specify_target(
    dgp, effect = "exposure", scale = "pp", engine = "bayes_laplace",
    prob_se_method = "joint_delta", on_singular = "flag_only"
  )
  sp_bt <- specify_target(
    dgp, effect = "exposure", scale = "pp", engine = "bayes_laplace",
    prob_se_method = "bootstrap", prob_boot = 49L, on_singular = "flag_only"
  )

  ev_jd <- evaluate_design_point(sp_jd, n_subj = 40L, n_item = 24L, nsim = 8L, stream_seeds = 901L)
  ev_bt1 <- evaluate_design_point(sp_bt, n_subj = 40L, n_item = 24L, nsim = 8L, stream_seeds = 901L)
  expect_equal(ev_jd$n_valid, 8L)
  expect_equal(ev_bt1$n_valid, 8L)
  expect_equal(ev_jd$prob_diff_est_pp, ev_bt1$prob_diff_est_pp, tolerance = 1e-10)
  # Mean 95% CI widths in percentage points agree closely between analytical joint delta and parametric bootstrap
  expect_lt(abs(ev_jd$ci_width_mean - ev_bt1$ci_width_mean) / ev_bt1$ci_width_mean, 0.20)
})

test_that("3AFC joint Laplace VCOV widens fixed-Sigma Schur SE and converges on extreme zero_item data", {
  dgp <- three_choice_reference_dgp()
  # Default engine omits warning and resolves cleanly to cpp_laplace when chance > 0
  expect_no_warning(sp_def <- specify_target(dgp, effect = "exposure"))
  expect_identical(sp_def$engine, "cpp_laplace")
  expect_message(sp_tmb <- specify_target(dgp, effect = "exposure", engine = "glmmTMB"), "glmmTMB")
  expect_identical(sp_tmb$engine, "cpp_laplace")

  sp_by <- specify_target(dgp, effect = "exposure", engine = "bayes_laplace", on_singular = "flag_only")
  sched <- generate_design_schedule(dgp$design, n_subj = 48L, n_item = 32L)
  set.seed(7312001L, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion", sample.kind = "Rejection")
  dat <- simulate_dgp_data(dgp, schedule = sched)

  X <- dgp$design$X_cells[dat$cell_id, , drop = FALSE]
  Zs <- dgp$design$Z_subj_cells[dat$cell_id, , drop = FALSE]
  Zi <- dgp$design$Z_item_cells[dat$cell_id, , drop = FALSE]
  raw_fit <- cpp_fit_crossed_glmm(
    y = as.integer(dat$y), X = X,
    subj_idx = as.integer(dat$subj_id) - 1L, Z_subj = Zs, n_subj = 48L,
    item_idx = as.integer(dat$item_id) - 1L, Z_item = Zi, n_item = 32L,
    testlet_idx = integer(nrow(dat)), Z_testlet = matrix(0, nrow(dat), 0L), n_testlet = 0L,
    chance = 1 / 3, corr_mode = 1L, regularized = TRUE,
    compute_joint_vcov = TRUE
  )
  expect_true(raw_fit$converged)
  eff_i <- match(sp_by$effect_col, colnames(X))
  expect_gt(sqrt(raw_fit$vcov[eff_i, eff_i]), sqrt(raw_fit$vcov_cond[eff_i, eff_i]))
})

test_that("report_plan accurately reports ci_width, rope, logit coverage, and heterogeneous slope SDs", {
  des <- design_crossed(list(A = c("a1", "a2"), B = c("b1", "b2")))
  dgp <- calibrate_dgp(
    des,
    cell_probs = c(a1_b1 = 0.50, a1_b2 = 0.55, a2_b1 = 0.65, a2_b2 = 0.80),
    subj_sd = c(0.70, 0.15, 0.25, 0.40),
    item_sd = c(0.60, 0.10, 0.20, 0.30)
  )

  sp_w <- specify_target(
    dgp, effect = "A:B", scale = "logit", criterion = "ci_width",
    ci_width = 2.00, ci_assurance = 0.70, engine = "cpp_laplace", on_singular = "flag_only"
  )
  res_w <- optimize_sample_size(
    sp_w, n_subj_range = c(32L, 40L, 48L), n_item_range = c(24L, 28L, 32L),
    method = "grid", nsim_final = 16L, nsim_topup = 0L, nsim_sensitivity = 0L,
    seed = 515L, verbose = FALSE
  )
  rep_w <- report_plan(res_w, lang = "en", print_output = FALSE)
  expect_match(rep_w, "AIPE width-attainment probability")
  expect_match(rep_w, "logit coverage")
  expect_match(rep_w, "0.400")

  sp_r <- specify_target(
    dgp, effect = "A:B", scale = "logit", criterion = "rope",
    rope = c(-0.10, 0.10), power_min = 0.60, engine = "cpp_laplace", on_singular = "flag_only"
  )
  res_r <- optimize_sample_size(
    sp_r, n_subj_range = c(32L, 40L, 48L), n_item_range = c(24L, 28L, 32L),
    method = "grid", nsim_final = 16L, nsim_topup = 0L, nsim_sensitivity = 0L,
    seed = 516L, verbose = FALSE
  )
  rep_r <- report_plan(res_r, lang = "en", print_output = FALSE)
  expect_match(rep_r, "ROPE")
})

test_that("plan_rasch_precision(cpp_laplace) and plan_rasch_fixed_bank handle extreme items accurately", {
  des <- design_crossed(list())
  dgp <- calibrate_dgp(des, 0.65, subj_sd = 0.8, item_sd = 0.6)
  sizes <- data.frame(n_subj = c(40L, 80L), n_item = c(12L, 16L))
  lim <- c(ability = 0.55, difficulty = 0.45, difficulty_max_abs_error = 0.80)

  p_cpp1 <- plan_rasch_precision(dgp, sizes, lim, nsim = 12L, seed = 88L, engine = "cpp_laplace", workers = 1L)
  p_cpp2 <- plan_rasch_precision(dgp, sizes, lim, nsim = 12L, seed = 88L, engine = "cpp_laplace", workers = 2L)
  expect_equal(p_cpp1$replications, p_cpp2$replications, tolerance = 1e-12)

  p_lme4 <- plan_rasch_precision(dgp, sizes[1L, , drop = FALSE], lim, nsim = 8L, seed = 88L, engine = "lme4")
  expect_lt(abs(p_cpp1$summary$rmse[1L] - p_lme4$summary$rmse[1L]), 0.03)

  # Fixed bank with extreme item (+4.5 logit -> frequent zero-correct at N_subj = 24)
  wide_bank <- c(seq(-1.2, 1.2, length.out = 11L), 4.5)
  fb_fix <- plan_rasch_fixed_bank(
    wide_bank, n_subj = c(24L, 80L),
    max_error = c(difficulty = 0.50, difficulty_max_abs_error = 1.00),
    method = "fixed_item", engine = "mml", nsim = 12L, seed = 314L
  )
  fb_rnd <- plan_rasch_fixed_bank(
    wide_bank, n_subj = c(24L, 80L),
    max_error = c(difficulty = 0.50, difficulty_max_abs_error = 1.00),
    method = "random_item", engine = "cpp_laplace", nsim = 12L, seed = 314L
  )

  # At N_subj = 24, extreme item has zero-correct replications: fixed_item marks bank unevaluable, random_item shrinks finitely
  s_fix_24 <- subset(fb_fix$summary, n_subj == 24L & target == "difficulty")
  s_rnd_24 <- subset(fb_rnd$summary, n_subj == 24L & target == "difficulty")
  expect_gt(s_fix_24$n_extreme_item_reps, 0L)
  expect_lt(s_fix_24$n_evaluable, 12L)
  expect_equal(s_rnd_24$n_evaluable, 12L)

  if (requireNamespace("TAM", quietly = TRUE)) {
    mod_bank <- seq(-1.5, 1.5, length.out = 10L)
    fb_mml <- plan_rasch_fixed_bank(
      mod_bank, n_subj = 80L, max_error = c(difficulty = 0.40),
      method = "fixed_item", engine = "mml", nsim = 5L, seed = 777L
    )
    fb_tam <- plan_rasch_fixed_bank(
      mod_bank, n_subj = 80L, max_error = c(difficulty = 0.40),
      method = "fixed_item", engine = "tam", nsim = 5L, seed = 777L
    )
    expect_equal(fb_mml$replications$difficulty, fb_tam$replications$difficulty, tolerance = 0.02)
  }
})
