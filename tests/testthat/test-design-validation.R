test_that("joint attainment counts the intersection and pools replication counts", {
  g <- calibrate_dgp(design_crossed(list(A = c("a", "b"))), c(.4, .6),
                     subj_sd = .3, item_sd = .4)
  sp <- specify_target(g, "A", criterion = "joint", power_min = .8, ci_width = 1,
                       engine = "bayes_laplace", analysis_re = "intercept_only")
  ev <- function(s) summarize_mc_metrics(c(rep(2, 85), rep(0, 15)),
    c(rep(.1, 70), rep(.5, 15), rep(.1, 15)), singular = rep(FALSE, 100),
    converged = rep(TRUE, 100), spec = s, n_subj = 8, n_item = 4, nsim = 100,
    engine_use = s$engine)
  r <- ev(sp)
  expect_equal(c(r$power, r$aipe_prob, r$joint_prob), c(.85, .85, .70))
  expect_equal(r$k_joint, 70L)
  expect_equal(r$decision_status, "fail")
  sp$criterion <- "both"
  expect_true(ev(sp)$criterion_met)
  sp$criterion <- "joint"
  sp$ci_width <- 100
  sp$power_min <- .5
  res <- optimize_sample_size(sp, 16, 8, method = "grid", nsim_final = 4L,
    nsim_sensitivity = 0L, workers = 1L, verbose = FALSE)
  res$spec$power_min <- mean(c(res$verified$joint_mc_lo, res$verified$joint_mc_hi))
  res$verified$decision_status <- "indeterminate"
  updated <- topup_verification(res, 3L, verbose = FALSE)$verified
  expect_equal(updated$nsim, 7L)
  expect_equal(updated$pattern_nsim, 7L)
  expect_equal(updated$pattern_n_invalid, updated$nsim - updated$n_valid)
  legacy <- res
  legacy$verified <- legacy$verified[setdiff(names(legacy$verified), c(pattern_columns, "replications"))]
  topped <- topup_verification(legacy, 3L, verbose = FALSE)$verified
  expect_equal(topped$nsim, 7L)
  expect_equal(topped$pattern_nsim, 3L)
  expect_equal(updated$joint_prob, updated$k_joint / updated$nsim)
  expect_equal(updated$decision_status, classify_mc_threshold(updated$joint_mc_lo,
    updated$joint_mc_hi, res$spec$power_min))
})

test_that("the null matches the actual estimand and preserves variance assumptions", {
  d <- design_crossed(list(A = c("a", "b"), B = c("c", "d")))
  g <- calibrate_dgp(d, c(.36, .46, .57, .83), chance = 1/3, subj_sd = .45, item_sd = .3)
  for (scale in c("pp", "prob", "logit")) {
    for (estimand in c("marginal", "conditional")) {
      sp <- specify_target(g, "A:B", scale = scale, prob_estimand = estimand,
                           engine = "bayes_laplace")
      null <- null_target_spec(sp)
      expect_lt(abs(null$true_effect_target), 1e-8)
      expect_identical(null$dgp$Sigma_subj, g$Sigma_subj)
      if (scale == "pp") expect_gt(abs(null$true_beta), .01)
    }
  }
  sp <- specify_target(g, "A", scale = "pp", contrast_weights = c(-1, 0, 0, 1),
                       engine = "bayes_laplace", analysis_re = "intercept_only", prob_boot = 3)
  out <- validate_design_point(sp, 16, 8, nsim = 3, seed = 48)
  expect_equal(out$true_target, c(0, 47), tolerance = 1e-6)
  expect_equal(out$condition, c("null", "planned"))
  expect_equal(out$detection[1], (1 - out$coverage[1]) * out$valid_rate[1])
})

test_that("missing clusters, empty datasets and passage interactions agree across generators", {
  expect_error(design_crossed(list(A = c("a", "b")), between_subj = "A",
                              between_item = "A"), "both between")
  d <- design_testlet(list(A = c("a", "b")), testlet_size = 2,
                      attrition = .4, participant_dropout = .5)
  g <- calibrate_dgp(d, c(.5,.7), subj_sd = .35, item_sd = .4, participant_testlet_sd = .8)
  sp <- specify_target(g, "A", scale = "pp", engine = "bayes_laplace")
  sched <- generate_design_schedule(d, 4, 4)
  expect_equal(sched$n_trials_expected, round(16 * .6 * .5))
  expect_equal(sp$true_prob_diff_pp, 20, tolerance = 1e-6)
  cache <- prepare_cpp_schedule_cache(sp, sched)
  for (seed in 1:4) {
    set.seed(seed)
    df <- simulate_dgp_data(g, sched)
    set.seed(seed)
    fast <- simulate_fast_y_keep(cache)
    expect_equal(df$y, fast$y)
    expect_equal(df$subj_id, sched$trials$subj_id[fast$keep])
  }
  local <- trial_testlet_effect(g, 2, 2, c(1,1,2,2,1), c(1,1,1,1,2))
  expect_equal(local[c(1,3)], local[c(2,4)])
  expect_equal(length(unique(local)), 3L)
  g$design$participant_dropout <- 1 - 1e-12
  sp$dgp <- g
  for (engine in c("bayes_laplace", "lme4")) {
    r <- evaluate_design_point(sp, 4, 4, nsim = 2, engine_use = engine)
    expect_equal(r$n_valid, 0)
    expect_equal(r$n_trials_obs_mean, 0)
    expect_equal(r$n_invalid_subj_unobserved, 2L)
    expect_equal(r$n_invalid_item_unobserved, 2L)
    expect_equal(r$n_rep_item_all_incorrect, 0L)
    expect_equal(r$effective_power, 0)
  }
})

test_that("the combined frontier includes a compromise absent from individual frontiers", {
  grid <- expand.grid(n_subj = c(8,16,24), n_item = c(8,16,24))
  grid$preference_score <- grid$n_subj + grid$n_item
  surface <- function(ok) {
    x <- grid
    for (nm in c("pred_power","pred_aipe_prob","pred_joint_prob","pred_rope_prob","prob_feasible"))
      x[[nm]] <- ifelse(ok, .9, .4)
    x$pred_ci_width <- x$pred_prob_ci_width_pp <- 1
    x$pred_criterion_met <- ok
    x
  }
  surfaces <- list(surface(grid$n_subj >= 16), surface(grid$n_item >= 16))
  r <- robust_surface(grid, surfaces)
  found <- select_stage2_candidates(r, list(), n_verify = 1)
  expect_equal(c(r$n_subj[found], r$n_item[found]), c(16,16))
})
