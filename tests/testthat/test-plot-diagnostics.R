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
  annotated <- plot(check, show_annotations = TRUE)
  expect_match(annotated$labels$subtitle, "8 items in bank / 4 items per participant")
  expect_match(annotated$labels$caption, "90% Wilson")
  expect_match(annotated$labels$caption, "Detection uses all scheduled studies")
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
  missing <- plot(empty, show_annotations = TRUE)
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

test_that("contour halo, single-point Pareto curve, and Rasch attainment labels build cleanly", {
  des <- design_crossed(list(condition = c("a", "b")))
  dgp <- calibrate_dgp(des, c(a = .40, b = .65), subj_sd = .3, item_sd = .3)
  sp_w <- specify_target(dgp, "condition", scale = "pp", criterion = "ci_width",
                         ci_width = 25, prob_se_method = "joint_delta")
  plan_w <- optimize_sample_size(
    sp_w, n_subj_range = c(12L, 36L), n_item_range = c(8L, 20L),
    method = "grid", nsim_final = 6L, workers = 1L, seed = 91L, verbose = FALSE,
    control = optimize_control(nsim_sensitivity = 0L)
  )
  p_cont <- plot(plan_w, type = "contour", show_annotations = TRUE)
  expect_match(p_cont$labels$caption, "Dashed black = target width assurance", fixed = TRUE)
  expect_warning(ggplot2::ggplot_build(p_cont), NA)

  # Single-point Pareto frontier omits geom_line and uses integer breaks without warning
  plan_1 <- plan_w
  plan_1$pareto_frontier <- plan_w$optimal
  p_par1 <- plot(plan_1, type = "pareto")
  expect_warning(b_par1 <- ggplot2::ggplot_build(p_par1), NA)
  expect_equal(glmmplanr:::integer_breaks(c(12.21, 12.39)), 12)
  expect_true(all(abs(glmmplanr:::integer_breaks(c(16, 24)) %% 1) < 1e-8))

  # Random-item and fixed-bank Rasch plots use exact candidate breaks and guard NA attainment
  dgp_r <- calibrate_dgp(design_crossed(list()), 0.55, subj_sd = 0.8, item_sd = 0.6)
  rp <- plan_rasch_precision(
    dgp_r, expand.grid(n_subj = c(16L, 24L), n_item = c(8L, 12L)),
    max_error = c(ability = 0.9), nsim = 4L, seed = 12L, verbose = FALSE
  )
  p_rp <- plot(rp, show_annotations = TRUE)
  expect_identical(p_rp$labels$title, "Rasch recovery: joint target")
  expect_warning(ggplot2::ggplot_build(p_rp), NA)

  rp_na <- rp
  rp_na$summary$attainment[1L] <- NA_real_
  p_rp_na <- plot(rp_na)
  expect_warning(b_rp_na <- ggplot2::ggplot_build(p_rp_na), NA)
  expect_false(any(b_rp_na$data[[1L]]$label == "NA%", na.rm = TRUE))

  fb <- plan_rasch_fixed_bank(c(-1, 0, 1), c(18L, 26L), c(difficulty = 0.8), nsim = 4L, seed = 13L)
  p_fb <- plot(fb)
  expect_warning(b_fb <- ggplot2::ggplot_build(p_fb), NA)
  expect_equal(b_fb$layout$panel_params[[1L]]$x$breaks, c(18, 26))
})
