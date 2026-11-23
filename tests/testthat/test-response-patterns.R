test_that("response patterns retain observation counts and scheduled but unseen units", {
  id <- factor(c("p1", "p1", "p2", "p2", "p3"), levels = paste0("p", 1:4))
  item <- c("i1", "i2", "i2", "i3", "i1")
  d <- response_patterns(c(TRUE, TRUE, FALSE, FALSE, NA), id, item)
  expect_equal(d$id, c(levels(id), "i1", "i2", "i3"))
  expect_equal(d$n_observed, c(2, 2, 0, 0, 1, 2, 1))
  expect_equal(d$n_correct, c(2, 0, 0, 0, 1, 1, 0))
  expect_equal(d$pattern, c("all_correct", "all_incorrect", "unobserved", "unobserved",
                            "all_correct", "mixed", "all_incorrect"))
  expect_equal(response_patterns(numeric(), id[FALSE], character())$pattern, rep("unobserved", 4))
  expect_error(response_patterns(2, "p", "i"), "0, 1")
  expect_error(response_patterns(1, NA_character_, "i"), "participant")
  expect_error(response_patterns(1, "p", 1:2), "item")
})

test_that("unpenalized factorial fits distinguish separated cells from constant units", {
  d <- design_crossed(list(group = c("A", "B")), allocation = "latin_square")
  g <- calibrate_dgp(d, c(.55, .75), chance = 1/3, subj_sd = .6, item_sd = .45,
                     slope_sd_ratio = 0)
  s <- generate_design_schedule(d, 24, 12)
  sp <- specify_target(g, "group", engine = "cpp_laplace", analysis_re = "intercept_only")
  set.seed(131)
  tr <- simulate_dgp_data(g, s)
  tr$y[tr$subj_id == 1L | tr$item_id == 1L] <- 0L
  expect_true(fit_single_dataset(tr, sp)$converged)

  cache <- prepare_cpp_schedule_cache(sp, s)
  keep <- tr$subj_id != 1L
  for (value in 0:1) {
    tr$y[tr$cell_id == 2L] <- value
    for (engine in c("cpp_laplace", "lme4", "lme4_fast", "glmmTMB")) {
      fit <- fit_single_dataset(tr, sp, engine_use = engine)
      expect_identical(fit$status, "separation")
      expect_match(fit$error_reason, "cell 'B'")
    }
    cached <- fit_via_cpp_cached(list(y = tr$y[keep], keep = keep), cache, sp,
                                 FALSE, "intercept_only", FALSE)
    expect_identical(cached$status, "separation")
  }
  expect_true(fit_single_dataset(tr, sp, engine_use = "bayes_laplace")$converged)
})

test_that("design diagnostics agree with observed trials across execution paths and top-ups", {
  des <- design_crossed(list(A = c("a", "b")), attrition = .2, participant_dropout = .3)
  g <- calibrate_dgp(des, c(.6, .8), subj_sd = .4, item_sd = .3)
  sp <- specify_target(g, "A", engine = "bayes_laplace", analysis_re = "intercept_only",
                       on_singular = "flag_only")
  sched <- generate_design_schedule(des, 8, 4)
  base <- allocate_stage_point_streams(6042L, 1L, 1L, 1L)
  streams <- derive_rep_substreams(base[[1]], 5L)
  totals <- empty_pattern_totals
  for (seed in streams) {
    d <- with_rng_stream(seed, simulate_dgp_data(g, sched))
    by <- with(d, response_patterns(y, subj, item))
    rep <- with_rng_stream(seed, execute_single_replication(sp, sched, NULL, sp$engine,
                                                            sp$analysis_re, FALSE, FALSE))
    present <- lapply(c("participant", "item"), function(unit)
      c("all_correct", "all_incorrect", "unobserved") %in% by$pattern[by$unit == unit])
    expect_equal(unname(rep$patterns[paste0("n_rep_", unit_pattern_names)]),
                 as.integer(unlist(present)))
    expect_equal(unname(rep$patterns["max_obs_item_all_incorrect"]),
                 max(0L, by$n_observed[by$unit == "item" & by$pattern == "all_incorrect"]))
    totals <- combine_pattern_totals(totals, rep$patterns)
  }
  point <- evaluate_design_point(sp, 8, 4, nsim = 5, stream_seeds = streams)
  expect_equal(unlist(point[pattern_columns]), totals)
  expect_equal(point$pattern_n_invalid, point$nsim - point$n_valid)
  expect_equal(validate_design_point(sp, 8, 4, nsim = 2)$pattern_nsim, c(2L, 2L))

  run_batch <- function(start, n) evaluate_grid_batch(sp, data.frame(n_subj = 8, n_item = 4),
    start, n, sp$engine, FALSE, 1L, "sequential", keep_raw = TRUE)
  first <- run_batch(base, 3L)
  second <- run_batch(list(streams[[3]]), 2L)
  combined <- pool_replications(first$raw[[1]], second$raw[[1]], sp, first$summary, second$summary)
  expect_equal(unlist(combined[pattern_columns]), totals)
  legacy <- first$summary[setdiff(names(first$summary), pattern_columns)]
  continued <- pool_replications(first$raw[[1]], second$raw[[1]], sp, legacy, second$summary)
  expect_equal(continued$nsim, 5L)
  expect_equal(continued[pattern_columns], second$summary[pattern_columns])
})

test_that("failed probability intervals count against validity, including wholly unobserved data", {
  fit <- list(est = 0, se = 1, converged = TRUE, status = "converged", prob_est = .1, prob_se = NA_real_)
  patterns <- replication_patterns(c(0L, 0L, 1L), c(1L, 2L, 1L), c(1L, 1L, 2L), 3, 3, fit, "pp")
  expect_equal(unname(patterns[c("pattern_n_invalid", "n_invalid_item_all_incorrect",
    "n_invalid_subj_unobserved", "max_obs_item_all_incorrect")]), c(1, 1, 1, 2))
  expect_equal(replication_patterns(0L, 1L, 1L, 1, 1, fit, "logit")[["pattern_n_invalid"]], 0)
})
