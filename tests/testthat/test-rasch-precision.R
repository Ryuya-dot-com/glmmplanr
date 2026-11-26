test_that("Rasch losses agree with direct IRT simulation and an independent fitted model", {
  design <- design_crossed(list())
  g <- calibrate_dgp(design, 0.65, subj_sd = 0.8, item_sd = 0.6)
  sizes <- data.frame(n_subj = 48, n_item = 16)
  limits <- c(ability = 0.6, difficulty = 0.4, difficulty_max_abs_error = 0.4,
              person_sd = 0.2, item_sd = 0.2,
              person_variance = 0.25, item_variance = 0.25)
  set.seed(231)
  old_seed <- .Random.seed
  plan <- plan_rasch_precision(g, sizes, limits, nsim = 3, seed = 64, engine = "lme4")
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
  errors <- rasch_recovery_errors(draw, g, control, engine = "lme4")
  expect_true(errors$singular)
  expect_true(is.na(errors$fit_error))
  expect_equal(c(errors$person_sd, errors$item_sd), c(0, 0), tolerance = 1e-5)
  draw$trials <- draw$trials[draw$trials$subj_id != 4L, ]
  errors <- rasch_recovery_errors(draw, g, control, engine = "lme4")
  expect_true(is.na(errors$ability))
  expect_true(is.finite(errors$difficulty))
})

test_that("Rasch search validates its inputs and retains only non-dominated passes", {
  g <- calibrate_dgp(design_crossed(list()), 0.6, subj_sd = 0.6, item_sd = 0.6)
  sizes <- data.frame(n_subj = c(40, 32, 40), n_item = c(12, 12, 12))
  plan <- plan_rasch_precision(g, sizes, c(ability = 10), nsim = 3, assurance = 0.2)
  expect_match(report_plan(plan, FALSE), "Non-dominated passing designs.*32 x 12")
  expect_s3_class(plot(plan), "ggplot")
  expect_error(plot(plan, type = "items"), "fixed bank")
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
                                     item_sd = 0.6, control = calibrate_control(testlet_sd = 0.3)),
                       calibrate_dgp(design_crossed(list()), 0.65, chance = 1/3,
                                     subj_sd = 0.8, item_sd = 0.6))) {
    expect_error(plan_rasch_precision(invalid, sizes, c(ability = 1)), "Rasch")
  }
})

test_that("fixed-bank diagnostics count each item's observed rows, including failed fits", {
  b <- seq(-1.5, 1.5, length.out = 6L)
  des <- design_crossed(list(), allocation = "crossed", items_per_subject = 2L)
  schedule <- generate_design_schedule(des, 12L, length(b))
  stream <- derive_rep_substreams(allocate_stage_point_streams(1L, 1L, 1L, 1L)[[1L]], 1L)[[1L]]
  sim <- with_rng_stream(stream, simulate_fixed_bank_trials(des, schedule, b, .8))
  trials <- sim$trials
  n <- tabulate(trials$item_id, length(b))
  yes <- tabulate(trials$item_id[trials$y == 1], length(b))
  out <- plan_rasch_fixed_bank(b, 12L, c(difficulty = 1), items_per_subject = 2L,
                               nsim = 1L, seed = 1L)
  expect_equal(out$item_summary$n_zero_correct, as.integer(n > 0 & yes == 0))
  expect_equal(out$item_summary$n_all_correct, as.integer(n > 0 & n == yes))
  for (value in 0:1) {
    sim$trials <- trials[trials$item_id != 2L, ]
    sim$trials$y[] <- value
    obs <- tabulate(sim$trials$item_id, length(b)) > 0L
    fit <- eval_fixed_bank_replication(sim, b, .8, "fixed_item", "mml")
    expect_true(fit$metrics$extreme_item)
    expect_identical(fit$item_all1, obs & value == 1)
    expect_identical(fit$item_zero, obs & value == 0)
  }
})

test_that("fixed-bank EAP at both extreme scores agrees with adaptive integration", {
  set.seed(400)
  n <- 48L
  b <- seq(-1.5, 1.5, length.out = 12L)
  theta <- rnorm(n)
  tr <- expand.grid(subj_id = seq_len(n), item_id = seq_along(b))
  tr$y <- rbinom(nrow(tr), 1L, plogis(theta[tr$subj_id] - b[tr$item_id]))
  tr$y[tr$subj_id == 1L] <- 0L
  tr$y[tr$subj_id == 2L] <- 1L
  fit <- fit_rasch_fixed_mml(tr, n, length(b), rep(TRUE, length(b)))
  expect_true(is.na(fit$error))
  for (value in 0:1) {
    posterior <- function(x) {
      likelihood <- vapply(x, function(t) prod(plogis(t - fit$b_est, lower.tail = value == 1L)),
                           numeric(1))
      likelihood * dnorm(x, sd = sqrt(fit$person_var))
    }
    mass <- integrate(posterior, -Inf, Inf, abs.tol = 1e-10)$value
    mean_theta <- integrate(function(x) x * posterior(x), -Inf, Inf, abs.tol = 1e-10)$value / mass
    expect_lt(abs(fit$ability[value + 1L] - mean_theta), 5e-5)
  }
})

test_that("a participant with only constant items has no fixed-bank ability estimate", {
  b <- seq(-1.4, 1.4, length.out = 8L)
  design <- design_crossed(list(), items_per_subject = 3L)
  schedule <- generate_design_schedule(design, 40L, 8L)
  stream <- derive_rep_substreams(allocate_stage_point_streams(1L, 1L, 1L, 1L)[[1L]], 1L)[[1L]]
  sim <- with_rng_stream(stream, simulate_fixed_bank_trials(design, schedule, b, .8))
  sim$trials <- rbind(sim$trials[c("subj_id", "item_id", "y")],
                      data.frame(subj_id = 41L, item_id = 9L, y = 1L))
  sim$n_subj <- 41L
  sim$n_item <- 9L
  sim$theta <- c(sim$theta, 0)
  for (engine in c("mml", if (requireNamespace("TAM", quietly = TRUE)) "tam")) {
    fit <- eval_fixed_bank_replication(sim, c(b, 0), .8, "fixed_item", engine)
    expect_true(is.na(fit$metrics$ability))
    expect_true(is.na(fit$b_est[9L]))
    expect_true(all(is.finite(fit$b_est[1:8])))
  }
})

test_that("slow fixed-bank EM agrees with an independent marginal likelihood optimizer", {
  # Saved 200-iteration failure: narrow population, case 5, replication 7.
  stream <- make_lecuyer_streams(91661001L, 7L)[[7L]]
  y <- with_rng_stream(stream, {
    theta <- .4 * rnorm(480L)
    uniform <- matrix(runif(480L * 16L), 480L, 16L)
    1L * (uniform[1:120, ] < plogis(outer(theta[1:120], seq(-1.8, 1.8, length.out = 16), "-")))
  })
  tr <- expand.grid(subj_id = 1:120, item_id = 1:16)
  tr$y <- as.integer(y)
  fit <- fit_rasch_fixed_mml(tr, 120L, 16L, rep(TRUE, 16L))
  expect_true(is.na(fit$error))
  expect_gt(fit$person_var, .01)
  # Higher-order quadrature, direct observed likelihood, no EM map.
  gh <- gh_quad_1d(61L)
  nll <- function(par) {
    eta <- outer(sqrt(par[17]) * gh$nodes, -par[1:16], "+")
    ll <- y %*% t(eta) + matrix(rowSums(log1p(-plogis(eta))), 120, 61, byrow = TRUE)
    ll <- sweep(ll, 2L, log(gh$weights), "+")
    peak <- apply(ll, 1L, max)
    -sum(peak + log(rowSums(exp(ll - peak))))
  }
  ref <- optim(c(-qlogis(colMeans(y)), .16), nll, method = "L-BFGS-B",
               lower = c(rep(-Inf, 16L), 0), control = list(factr = 1e3))
  expect_equal(ref$convergence, 0L)
  expect_equal(fit$b_est, ref$par[1:16], tolerance = 2e-4)
  expect_equal(sqrt(fit$person_var), sqrt(ref$par[17]), tolerance = 2e-4)
  expect_lt(abs(nll(c(fit$b_est, fit$person_var)) - ref$value), 1e-5)
})

test_that("fixed-bank variance boundary is finite for complete and partial forms", {
  for (partial in c(FALSE, TRUE)) {
    tr <- expand.grid(subj_id = 1:40, item_id = 1:12)
    tr$y <- with(tr, (subj_id + item_id) %% 2L)
    if (partial) tr <- subset(tr, subj_id <= 39 & (subj_id + item_id) %% 3 != 0)
    fit <- fit_rasch_fixed_mml(tr, 40L, 12L, rep(TRUE, 12L))
    expect_true(is.na(fit$error))
    expect_true(fit$singular)
    expect_identical(fit$person_var, 0)
    exact <- -qlogis(tapply(tr$y, tr$item_id, mean))
    expect_equal(fit$b_est, as.numeric(exact), tolerance = 1e-8)
    expect_true(all(fit$ability[unique(tr$subj_id)] == 0))
    if (partial) expect_true(is.na(fit$ability[40]))
  }
})
