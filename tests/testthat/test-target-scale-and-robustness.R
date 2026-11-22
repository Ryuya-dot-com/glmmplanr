test_that("detection, coverage and sign errors use the requested contrast", {
  des <- design_crossed(list(A = c("a", "b"), B = c("c", "d")))
  dgp <- calibrate_dgp(des, c(.02, .26, .26, .50), subj_sd = .3, item_sd = .3,
                       slope_sd_ratio = 0, corr_structure = 0)
  for (scale in c("prob", "pp", "logit")) {
    spec <- specify_target(dgp, "A:B", scale = scale, engine = "bayes_laplace")
    expect_equal(spec$true_prob_diff_pp, 0, tolerance = 1e-8)
    out <- summarize_mc_metrics(c(-2, -2), c(.1, .1), c(.01, -.1), c(.02, .02),
                                c(FALSE, FALSE), c(TRUE, TRUE), spec = spec,
                                n_subj = 8, n_item = 8, nsim = 2, engine_use = spec$engine)
    expect_equal(out$power, if (scale == "logit") 1 else .5)
    expect_equal(out$power_nonsingular, out$power)
    expect_equal(out$coverage_nonsingular, out$coverage)
    expect_equal(out$effective_logit_power, 1)
  }
  sp <- specify_target(dgp, "A:B", scale = "pp", engine = "bayes_laplace",
                       analysis_re = "intercept_only", on_singular = "flag_only", prob_boot = 9L)
  ev <- evaluate_design_point(sp, 64, 32, nsim = 4, stream_seeds = 1357)
  expect_equal(ev$n_valid, 4L)
  expect_equal(ev$effective_power, 1 - ev$prob_coverage)
  expect_gt(ev$effective_logit_power, ev$effective_power)

})

test_that("joint refits vary the marginal contrast even when beta stays fixed", {
  des <- design_crossed(list(A = c("a", "b")))
  dgp <- calibrate_dgp(des, c(.4, .7), subj_sd = .4, item_sd = .3)
  spec <- specify_target(dgp, "A", scale = "pp", engine = "bayes_laplace", prob_boot = 6)
  n <- 12
  x <- des$X_cells[rep(1:2, 6), ]
  counter <- 0L
  variances <- c(1, .25, .5, .75, 1, 1.25, 1.5)
  refit <- function(y) {
    counter <<- counter + 1L
    list(converged = TRUE, beta_full = c(0, 1), vcov_full = diag(0, 2),
         Sigma_subj = matrix(variances[counter]), Sigma_item = matrix(.5))
  }
  calc <- function(fun) summarize_probability_fit(rep(0L, n), fun, x, matrix(1, n),
    matrix(1, n), matrix(0, n, 0), rep(1:4, each = 3), rep(1:3, 4), NULL,
    spec, "intercept_only", 0, FALSE)
  set.seed(12)
  out <- calc(refit)
  expected <- vapply(variances[-1], function(v) estimate_prob_effect_delta(
    c(0, 1), diag(0, 2), matrix(v), matrix(.5), NULL, des, 0,
    spec$contrast_weights, re_structure = "intercept_only")$est, numeric(1))
  expect_equal(counter, 7L)
  expect_equal(out$prob_se, sd(expected))
  counter <- 0L
  failed <- calc(function(y) {
    fit <- refit(y)
    fit$converged <- counter == 1L
    fit
  })
  expect_true(is.na(failed$prob_se))
  expect_match(failed$error_reason, "0/6 valid refits")
})

test_that("only confirmed points can be recommended and preference preserves all scenarios", {
  d <- data.frame(scenario = "base", n_subj = c(8, 16, 24), n_item = c(20, 12, 8),
                  n_trials_scheduled = c(160, 192, 192), power = c(.801, .86, .88),
                  ci_width_mean = 12, preference_score = 1:3, criterion_met = TRUE,
                  decision_status = c("indeterminate", "pass", "pass"))
  sp <- list(criterion = "power", scale = "logit")
  expect_equal(select_recommended_from_verified(d, sp)$optimal$n_subj, 16)
  weak <- d
  weak$scenario <- "weak"
  weak$decision_status <- c("fail", "pass", "fail")
  weak$power <- c(.6, .85, .55)
  weak$criterion_met <- weak$power >= .8
  rc <- rank_scenarios(rbind(d, weak), "base", "weak", sp)
  res <- structure(list(spec = sp, robust = rc, verified = rc$verified,
                        surface = d, search_meta = list(), prefer = "balanced"),
                   class = "glmmplanr_result")
  for (prefer in c("fewer_subj", "fewer_items")) {
    result <- switch_preference(res, prefer)
    expect_true(result$feasible)
    expect_equal(result$optimal$n_subj, 16)
    expect_equal(result$optimal, result$robust$optimal_robust)
    expect_equal(result$optimal$n_subj, result$robust$comparison$n_subj[2])
    expect_equal(result$verified_pareto$n_subj, 16)
  }
  weak$decision_status[2] <- "indeterminate"
  rc <- rank_scenarios(rbind(d, weak), "base", "weak", sp)
  expect_false(rc$feasible_robust)
  expect_true(is.na(rc$delta_trials))
})

test_that("otherwise identical dropout scenarios are simulated separately", {
  make_dgp <- function(attrition) calibrate_dgp(
    design_crossed(list(A = c("a", "b")), attrition = attrition),
    c(.4, .7), subj_sd = .4, item_sd = .3)
  full <- make_dgp(0)
  half <- make_dgp(.5)
  expect_true(dgp_is_equivalent(full, make_dgp(0)))
  expect_false(dgp_is_equivalent(full, half))
  sp <- specify_target(full, "A", engine = "bayes_laplace", analysis_re = "intercept_only")
  result <- optimize_sample_size(sp, 12, 8, scenarios = list(base = full, half = half),
    standard_scenario = "base", method = "grid", nsim_final = 3, nsim_sensitivity = 0,
    workers = 1, verbose = FALSE)
  ev <- result$robust$scenario_evaluations
  expect_equal(ev$n_trials_expected[ev$scenario == "half"], 48L)
  expect_lt(ev$n_trials_obs_mean[ev$scenario == "half"], 96)
  expect_false(result$feasible)
  expect_match(report_plan(result, "en", FALSE), "No evaluated design")
  expect_match(report_plan(result, "ja", FALSE), "目標達成を確認できた設計はありません")
  expect_match(report_plan(result, "en", FALSE), "independent simulations")
  expect_match(report_plan(result, "ja", FALSE), "区間の被覆率は候補ごとに")
  topped <- topup_verification(result, nsim_extra = 3L, verbose = FALSE)
  expect_equal(topped$optimal, topped$robust$optimal_robust)
})

test_that("3AFC determinant matches an independently differentiated joint likelihood", {
  des <- design_crossed(list(A = c("a", "b")))
  dgp <- calibrate_dgp(des, c(.51, .75), chance = 1/3, subj_sd = .65,
                       item_sd = .5, slope_sd_ratio = 0)
  set.seed(920)
  df <- simulate_dgp_data(dgp, n_subj = 16, n_item = 10)
  X <- des$X_cells[df$cell_id, ]
  n <- nrow(X)
  f <- cpp_fit_crossed_glmm(df$y, X, df$subj_id - 1L, matrix(1, n), 16L,
    df$item_id - 1L, matrix(1, n), 10L, integer(0), matrix(0, n, 0), 0L,
    chance = 1/3, regularized = TRUE)
  expect_true(f$converged)
  Z <- cbind(model.matrix(~ 0 + subj, df) * sqrt(f$Sigma_subj[1]),
             model.matrix(~ 0 + item, df) * sqrt(f$Sigma_item[1]))
  nll <- deriv(~ -y * log(1/3 + 2/3 / (1 + exp(-eta))) -
                 (1-y) * log(2/3 / (1 + exp(eta))), "eta",
               function.arg = c("eta", "y"), hessian = TRUE)
  expect_equal(as.numeric(attr(nll(c(0, 0), 0:1), "hessian")), c(.25, .0625))
  likelihood <- function(u) nll(drop(X %*% f$beta + Z %*% u), df$y)
  objective <- function(u) sum(likelihood(u)) + sum(u^2)/2
  gradient <- function(u) drop(crossprod(Z, attr(likelihood(u), "gradient"))) + u
  mode <- optim(numeric(ncol(Z)), objective, gradient, method = "BFGS",
                control = list(reltol = 1e-13))$par
  weights <- as.numeric(attr(likelihood(mode), "hessian"))
  H <- diag(ncol(Z)) + crossprod(Z, weights * Z)
  expect_equal(f$logLik, -objective(mode) - sum(log(diag(chol(H)))), tolerance = 1e-7)
})
