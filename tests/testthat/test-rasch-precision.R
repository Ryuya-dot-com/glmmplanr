test_that("Rasch losses agree with direct IRT simulation and an independent fitted model", {
  design <- design_crossed(list())
  g <- calibrate_dgp(design, 0.65, subj_sd = 0.8, item_sd = 0.6)
  sizes <- data.frame(n_subj = 48, n_item = 16)
  limits <- c(ability = 0.6, difficulty = 0.4, difficulty_max_abs_error = 0.4,
              person_sd = 0.2, item_sd = 0.2,
              person_variance = 0.25, item_variance = 0.25)
  set.seed(231)
  old_seed <- .Random.seed
  plan <- plan_rasch_precision(g, sizes, limits, nsim = 3, seed = 64)
  expect_identical(.Random.seed, old_seed)
  schedule <- generate_design_schedule(design, 48, 16)
  stream <- derive_rep_substreams(allocate_stage_point_streams(64, 1L, 1L, 1L)[[1L]], 1L)[[1L]]
  reference <- with_rng_stream(stream, {
    theta <- rnorm(48, sd = 0.8)
    b <- -unname(g$beta[1]) - rnorm(16, sd = 0.6)
    tr <- schedule$trials
    tr$y <- rbinom(nrow(tr), 1L, plogis(theta[tr$subj_id] - b[tr$item_id]))
    model <- lme4::glmer(y ~ 1 + (1 | item) + (1 | subj), data = tr,
                         family = binomial(), control = lme4::glmerControl(optimizer = "bobyqa"))
    modes <- lme4::ranef(model, condVar = FALSE)
    v <- lme4::VarCorr(model)
    item_error <- b + lme4::fixef(model)[1] + modes$item[, 1]
    c(sqrt(mean((theta - modes$subj[, 1])^2)), sqrt(mean(item_error^2)), max(abs(item_error)),
      abs(sqrt(v$subj[1, 1]) - 0.8), abs(sqrt(v$item[1, 1]) - 0.6),
      abs(v$subj[1, 1] - 0.8^2), abs(v$item[1, 1] - 0.6^2))
  })
  expect_equal(as.numeric(plan$replications[1, names(limits)]), unname(reference), tolerance = 1e-5)
  losses <- as.matrix(plan$replications[names(limits)])
  met <- sweep(losses, 2L, limits, "<=")
  expect_equal(plan$summary$n_met, unname(c(colSums(met), sum(rowSums(met) == length(limits)))))
  expect_equal(head(plan$summary$rmse, -1), unname(sqrt(colMeans(losses^2))))
  expect_equal(tail(plan$summary$mc_lo, 1), unname(wilson_mc_ci(tail(plan$summary$n_met, 1), 3)[1]))
})

test_that("Rasch assurance counts unavailable predictions and failed fits as misses", {
  d <- design_crossed(list(), items_per_subject = 3L, item_seed = 67)
  g <- calibrate_dgp(d, 0.65, subj_sd = 0.8, item_sd = 0.6)
  sizes <- data.frame(n_subj = 30, n_item = 96)
  # Ninety planned responses cannot cover 96 different items.
  plan <- plan_rasch_precision(g, sizes, c(ability = 10, difficulty = 10,
                                          difficulty_max_abs_error = 10, item_sd = 10),
                               nsim = 4, seed = 218)
  s <- plan$summary
  expect_gt(s$n_valid_fit[1], 0)
  expect_equal(s$n_evaluable, c(s$n_valid_fit[1], 0, 0, s$n_valid_fit[1], 0))
  expect_equal(s$n_met, s$n_evaluable)
  expect_equal(s$n_singular[1], sum(plan$replications$singular, na.rm = TRUE))
  expect_identical(s$decision[s$target == "joint"], "fail")
  expect_equal(nrow(plan$pareto), 0L)

  dropout <- design_crossed(list(), participant_dropout = 1 - 1e-12)
  g <- calibrate_dgp(dropout, 0.65, subj_sd = 0.8, item_sd = 0.6)
  none <- plan_rasch_precision(g, sizes, c(ability = 10), nsim = 3, seed = 7)
  expect_true(all(!is.na(none$replications$fit_error)))
  expect_true(all(none$summary$n_met == 0 & none$summary$n_valid_fit == 0))
  expect_true(all(is.na(none$summary$rmse)))
})

test_that("boundary estimates and unobserved participants have distinct recovery outcomes", {
  design <- design_crossed(list())
  g <- calibrate_dgp(design, 0.55, subj_sd = 0, item_sd = 0)
  schedule <- generate_design_schedule(design, 32, 12)
  draw <- with_rng_stream(77, draw_dgp_trials(g, schedule, keep_person_item_effects = TRUE))
  draw$trials$y <- with(draw$trials, (subj_id + item_id) %% 2L)
  control <- lme4::glmerControl(optimizer = "bobyqa")
  errors <- rasch_recovery_errors(draw, g, control)
  expect_true(errors$singular)
  expect_true(is.na(errors$fit_error))
  expect_equal(c(errors$person_sd, errors$item_sd), c(0, 0), tolerance = 1e-5)
  draw$trials <- draw$trials[draw$trials$subj_id != 4L, ]
  errors <- rasch_recovery_errors(draw, g, control)
  expect_true(is.na(errors$ability))
  expect_true(is.finite(errors$difficulty))
})

test_that("Rasch search validates its inputs and retains only non-dominated passes", {
  g <- calibrate_dgp(design_crossed(list()), 0.6, subj_sd = 0.6, item_sd = 0.6)
  sizes <- data.frame(n_subj = c(40, 32, 40), n_item = c(12, 12, 12))
  plan <- plan_rasch_precision(g, sizes, c(ability = 10), nsim = 3, assurance = 0.2)
  expect_equal(plan$pareto$n_subj, 32)
  expect_equal(nrow(plan$replications), 6L)
  expect_equal(plan, plan_rasch_precision(g, sizes[3:1, ], c(ability = 10), nsim = 3, assurance = 0.2))
  expect_error(plan_rasch_precision(g, sizes, c(x = 0.1)), "max_error")
  expect_error(plan_rasch_precision(g, sizes, c(ability = -1)), "max_error")
  expect_error(plan_rasch_precision(g, sizes / 7, c(item_sd = 0.1)), "integer")
  expect_error(plan_rasch_precision(g, sizes, c(ability = 1), nsim = 2.5), "nsim")
  expect_error(plan_rasch_precision(g, sizes, c(ability = 1), assurance = 1), "assurance")
  for (invalid in list(calibrate_dgp(design_crossed(list(A = c("a", "b"))),
                                    c(0.5, 0.6), subj_sd = 0.8, item_sd = 0.6),
                       calibrate_dgp(design_testlet(list()), 0.6, subj_sd = 0.8,
                                     item_sd = 0.6, testlet_sd = 0.3),
                       calibrate_dgp(design_crossed(list()), 0.65, chance = 1/3,
                                     subj_sd = 0.8, item_sd = 0.6))) {
    expect_error(plan_rasch_precision(invalid, sizes, c(ability = 1)), "Rasch")
  }
})
