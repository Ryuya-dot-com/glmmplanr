test_that("scenario plots use verified rates and identify full and partial presentation", {
  full <- design_crossed(list(method = c("restudy", "retrieval")))
  partial <- design_crossed(list(method = c("restudy", "retrieval")), items_per_subject = 4L)
  calibrate <- function(design) calibrate_dgp(design, c(restudy = .42, retrieval = .55),
    subj_sd = .3, item_sd = .4)
  dgps <- list(full = calibrate(full), partial = calibrate(partial))
  spec <- specify_target(dgps$full, "method", scale = "pp", engine = "bayes_laplace",
                         prob_se_method = "joint_delta")
  plan <- optimize_sample_size(spec, 16L, c(8L, 10L), scenarios = dgps,
    standard_scenario = "full", method = "grid", nsim_final = 4L, workers = 1L, verbose = FALSE,
    control = optimize_control(nsim_init = 4L, n_verify = 2L, nsim_sensitivity = 0L))

  selected <- unique(rbind(plan$robust$optimal_standard, plan$robust$optimal_robust)[c("n_subj", "n_item")])
  expect_equal(nlevels(plot(plan, type = "scenarios")$data$design), nrow(selected))
  candidates <- as.data.frame(plan)[c("n_subj", "n_item")]
  figure <- plot(plan, type = "scenarios", candidates = candidates)
  expect_s3_class(figure, "ggplot")
  expect_equal(nlevels(figure$data$design), nrow(candidates))
  expect_true(all(c("full\n8 items per participant", "partial\n4 items per participant") %in%
                  as.character(figure$data$scenario)))
  evaluated <- plan$robust$scenario_evaluations
  detection <- subset(figure$data, panel == "Detection (%)")
  expect_equal(detection$value, 100 * evaluated$power)
  expect_equal(detection$lo, 100 * evaluated$power_mc_lo)
  expect_equal(detection$hi, 100 * evaluated$power_mc_hi)
  expect_equal(detection$decision, evaluated$decision_status)
  width <- subset(figure$data, panel == "Mean interval width (pp)")
  expect_equal(width$value, evaluated$prob_ci_width_mean_pp)
  expect_true(all(is.na(width$lo)))
  expect_warning(ggplot2::ggplot_build(figure), NA)
  expect_error(plot(plan, type = "scenarios", candidates = data.frame(n_subj = 18L, n_item = 8L)),
               "verified scenario results")

  plan$robust$specs <- NULL
  expect_true("partial\n4 scheduled trials per participant" %in%
    as.character(plot(plan, type = "scenarios", candidates = candidates)$data$scenario))
  plan$robust <- NULL
  expect_error(plot(plan, type = "scenarios"), "multiple `scenarios`")
})

test_that("validation plots carry stored intervals, denominators and extreme scores", {
  design <- design_crossed(list(method = c("a", "b")), items_per_subject = 4L)
  dgp <- calibrate_dgp(design, c(a = .4, b = .75), subj_sd = .3, item_sd = .4)
  spec <- specify_target(dgp, "method", engine = "bayes_laplace", alpha = .1,
                         control = target_control(mc_conf_level = .9))
  check <- validate_design_point(spec, 16L, 8L, nsim = 3L, seed = 43L)
  expect_s3_class(check, c("glmmplanr_validation", "data.frame"))
  expect_identical(check$condition, c("null", "planned"))
  rates <- plot(check)
  expect_match(rates$labels$subtitle, "8 items in bank / 4 items per participant")
  expect_match(rates$labels$caption, "90% Wilson")
  expect_match(rates$labels$caption, "Detection uses all scheduled studies")
  expect_equal(rates$data$rate, c(check$detection, check$coverage))
  expect_equal(rates$data$lo, c(check$detection_lo, check$coverage_lo))
  expect_equal(rates$data$hi, c(check$detection_hi, check$coverage_hi))
  expect_equal(rates$data$target, c(.1, spec$power_min, .9, .9))
  expect_true(all(grepl("3/3 valid fits", rates$data$studies, fixed = TRUE)))
  expect_warning(ggplot2::ggplot_build(rates), NA)

  patterns <- plot(check, type = "patterns")
  expect_equal(nrow(patterns$data), 12L)
  expect_equal(patterns$data$rate[1L], check$n_rep_subj_all_correct[1L] / check$pattern_nsim[1L])
  expect_true(any(grepl("max 4 responses", patterns$data$pattern, fixed = TRUE)))
  expect_true(all(vapply(split(as.integer(patterns$data$row), patterns$data$condition),
    function(rows) all(diff(rows) < 0L), logical(1L))))
  expect_warning(ggplot2::ggplot_build(patterns), NA)

  spec$failure_policy <- "conditional"
  spec$dgp$design$participant_dropout <- 1 - 1e-12
  empty <- validate_design_point(spec, 16L, 8L, nsim = 2L, seed = 43L)
  missing <- plot(empty)
  expect_true(all(grepl("0/2 valid fits", missing$data$studies, fixed = TRUE)))
  expect_true(all(is.na(missing$data$rate)))
  expect_match(missing$labels$caption, "Detection uses valid fits")
  expect_warning(ggplot2::ggplot_build(missing), NA)
})

test_that("absolute probability validation plots coverage and width attainment", {
  design <- design_crossed(list())
  dgp <- calibrate_dgp(design, .65, subj_sd = .3, item_sd = .4)
  spec <- specify_target(dgp, "(Intercept)", scale = "prob", criterion = "ci_width",
                         ci_width = .20, engine = "bayes_laplace", prob_se_method = "joint_delta")
  check <- suppressMessages(validate_design_point(spec, 16L, 8L, nsim = 3L))
  precision <- plot(check)
  expect_identical(as.character(precision$data$panel), c("Interval coverage", "Width attainment"))
  expect_equal(precision$data$rate, c(check$coverage, check$width_attainment))
  expect_equal(precision$data$target, c(1 - spec$alpha, spec$ci_assurance))
  expect_warning(ggplot2::ggplot_build(precision), NA)

  scenarios <- list(base = dgp, lower = calibrate_dgp(design, .45, subj_sd = .3, item_sd = .4))
  plan <- optimize_sample_size(spec, 16L, 8L, scenarios = scenarios, method = "grid",
    nsim_final = 3L, workers = 1L, verbose = FALSE, control = optimize_control(nsim_sensitivity = 0L))
  panels <- plot(plan, type = "scenarios")$data
  tested <- plan$robust$scenario_evaluations
  attainment <- subset(panels, panel == "Width attainment (%)")
  expect_equal(attainment$value, 100 * tested$aipe_prob)
  expect_equal(attainment$lo, 100 * tested$aipe_mc_lo)
})
