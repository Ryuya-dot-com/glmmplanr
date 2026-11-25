test_that("2x2 Latin square enforces 4-list balance and orthogonal sum contrasts", {
  des <- design_crossed(
    factors = list(
      exposure = c("incidental", "intentional"),
      spacing  = c("massed", "distributed")
    ),
    allocation = "latin_square",
    contrast   = "sum",
    attrition  = 0.05
  )

  expect_s3_class(des, "glmmplanr_design")
  expect_equal(des$k_subj, 4L)
  expect_equal(des$k_item, 4L)
  expect_equal(des$n_cells, 4L)

  # Sum contrast is (-0.5, +0.5)
  expect_equal(as.numeric(des$contrasts$exposure[, 1L]), c(-0.5, 0.5))
  expect_equal(as.numeric(des$contrasts$spacing[, 1L]), c(-0.5, 0.5))

  # Non-multiples rejected when strict = TRUE
  expect_error(generate_design_schedule(des, n_subj = 25, n_item = 16))
  expect_error(generate_design_schedule(des, n_subj = 24, n_item = 18))

  sched <- generate_design_schedule(des, n_subj = 24, n_item = 16)
  tr <- sched$trials
  expect_equal(nrow(tr), 24L * 16L)

  # Every subject sees each of the 4 cells exactly 16/4 = 4 times
  subj_cell_tab <- table(tr$subj_id, tr$cell_id)
  expect_true(all(subj_cell_tab == 4L))

  # Every item appears in each of the 4 cells exactly 24/4 = 6 times
  item_cell_tab <- table(tr$item_id, tr$cell_id)
  expect_true(all(item_cell_tab == 6L))

  # No subject sees the same item twice
  expect_true(all(table(tr$subj_id, tr$item_id) == 1L))
})

test_that("Between-subject x within-item and testlet presets enforce proper nesting and REs", {
  # Between-subject group x within-subject & within-item condition
  des_bw <- design_crossed(
    factors = list(
      group = c("control", "treatment"),
      time  = c("pretest", "posttest")
    ),
    between_subj = "group",
    contrast = "helmert"
  )
  # 2 between-subject groups x 2 Latin-square lists for `time` = 4 subject lists, 2 item sets
  expect_equal(des_bw$k_subj, 4L)
  expect_equal(des_bw$k_item, 2L)
  # Subject RE must NOT include between-subject `group`, only `time`
  expect_equal(ncol(des_bw$Z_subj_cells), 2L)
  expect_equal(ncol(des_bw$Z_item_cells), 4L)

  # Passage / testlet nested design
  des_tl <- design_testlet(
    factors = list(
      modality = c("reading", "listening"),
      gloss    = c("absent", "present")
    ),
    testlet_size = 3L,
    contrast = "sum"
  )
  expect_true(des_tl$has_testlet)
  expect_equal(des_tl$k_subj, 4L)
  expect_equal(des_tl$k_item, 12L) # 4 Latin square lists x 3 questions per passage

  sched_tl <- generate_design_schedule(des_tl, n_subj = 16, n_item = 24)
  expect_equal(sched_tl$n_testlet, 8L)
  # Within a subject and passage, all 3 nested questions share the same condition cell
  by_sp <- split(sched_tl$trials$cell_id, list(sched_tl$trials$subj_id, sched_tl$trials$testlet_id))
  expect_true(all(vapply(by_sp, function(v) length(unique(v)) == 1L, logical(1L))))
})

test_that("participant-specific item samples use a reproducible allocation before attrition", {
  f <- list(group = c("a", "b", "c"))
  des <- design_crossed(f, "between_subj", items_per_subject = 5, item_seed = 417,
                        attrition = .2, participant_dropout = .1)
  set.seed(77)
  before <- .Random.seed
  tr <- generate_design_schedule(des, 642, 20)$trials
  expect_identical(.Random.seed, before)
  expect_true(all(table(tr$subj) == 5L))
  expect_equal(unname(as.integer(table(tr$group))), rep(1070L, 3))
  expect_identical(anyDuplicated(tr[c("subj_id", "item_id")]), 0L)
  expect_identical(tr, generate_design_schedule(des, 642, 20)$trials)
  des$item_seed <- 418L
  expect_false(identical(tr$item_id, generate_design_schedule(des, 642, 20)$trials$item_id))
  expect_equal(compute_expected_trials(des, 642, c(20, 40)), rep(round(3210*.8*.9), 2))
  expect_error(generate_design_schedule(des, 642, 4), "items_per_subject")
  expect_error(design_crossed(f, items_per_subject = 1.5), "items_per_subject")
  expect_error(design_crossed(f, items_per_subject = 2, item_seed = NA), "item_seed")
  # If the whole bank is drawn, preserve the original schedule and row order.
  full <- design_crossed(f, "crossed")
  quota <- design_crossed(f, "crossed", items_per_subject = 6)
  expect_identical(generate_design_schedule(full, 6, 6), generate_design_schedule(quota, 6, 6))
})

test_that("partial presentation preserves repetitions, whole testlets and unselected bank levels", {
  f <- list(A = c("a", "b"))
  des <- design_crossed(f, "crossed", items_per_subject = 1)
  tr <- generate_design_schedule(des, 2, 40)$trials
  expect_identical(nlevels(tr$item), 40L)
  expect_true(all(table(tr$subj, tr$A) == 1))
  expect_equal(length(unique(tr$item_id[tr$subj_id == 1])), 1)
  expect_equal(compute_scheduled_trials(des, 2, 40), nrow(tr))
  des <- design_testlet(f, testlet_size = 2, items_per_subject = 4)
  tr <- generate_design_schedule(des, 8, 20)$trials
  expect_true(all(table(tr$subj, tr$testlet) %in% c(0, 2)))
  expect_true(all(table(tr$subj) == 4))
  expect_error(design_testlet(f, testlet_size = 2, items_per_subject = 3), "items_per_subject")
})

test_that("search and validation use admissible bank sizes and per-person workload", {
  des <- design_crossed(list(group = c("a","b")), "between_subj", items_per_subject = 5)
  dgp <- calibrate_dgp(des, c(.6,.55), subj_sd = .5, item_sd = .3, slope_sd_ratio = 0)
  sp <- specify_target(dgp, "group", engine = "bayes_laplace", analysis_re = "intercept_only")
  plan <- optimize_sample_size(sp, 16, c(4, 6), method = "grid", nsim_final = 2,
                               verbose = FALSE, control = optimize_control(nsim_sensitivity = 0))
  expect_equal(sort(unique(plan$verified$n_item)), 5:6)
  expect_true(all(plan$verified$n_trials_scheduled == 80))
  expect_match(report_plan(plan, print_output = FALSE), "5 distinct items")
  expect_error(optimize_sample_size(sp, 16, 4, verbose = FALSE), "bank")
  sparse <- generate_design_schedule(des, 16, 100)
  cache <- glmmplanr:::prepare_cpp_schedule_cache(sp, sparse)
  set.seed(744)
  cached <- glmmplanr:::simulate_fast_y_keep(cache)
  set.seed(744)
  expect_identical(cached$y, simulate_dgp_data(dgp, sparse)$y)
})
test_that("an empty factor list generates a one-cell crossed random-item Rasch model", {
  design <- design_crossed(list(), allocation = "crossed")
  expect_identical(design$cells$cell_name, "overall")
  for (matrix in list(design$X_cells, design$Z_subj_cells, design$Z_item_cells)) {
    expect_identical(colnames(matrix), "(Intercept)")
    expect_equal(as.numeric(matrix), 1)
  }
  expect_error(design_crossed(list(c("a", "b"))), "named list")
  expect_error(design_crossed(NULL), "named list")
  dgp <- calibrate_dgp(design, c(overall = 0.65), subj_sd = 0.8, item_sd = 0.6)
  beta <- unname(dgp$beta["(Intercept)"])
  marginal <- integrate(function(z) plogis(beta + z) * dnorm(z), -Inf, Inf)$value
  expect_equal(marginal, 0.65, tolerance = 1e-7)
  set.seed(415)
  data <- simulate_dgp_data(dgp, n_subj = 36, n_item = 12)
  expect_equal(nrow(data), 36 * 12)
  expect_true(all(data$cell_id == 1L))
  set.seed(415)
  ability <- rnorm(36, sd = 0.8)
  difficulty <- -beta - rnorm(12, sd = 0.6)
  prob <- plogis(ability[data$subj_id] - difficulty[data$item_id])
  expect_identical(data$y, rbinom(nrow(data), 1, prob))
  partial <- design_testlet(list(), testlet_size = 2L, items_per_subject = 4L)
  trials <- generate_design_schedule(partial, 7, 8)$trials
  expect_equal(nrow(trials), 28L)
  expect_true(all(trials$cell_id == 1L))
  expect_true(all(as.vector(table(trials$subj, trials$testlet)) %in% c(0L, 2L)))
  expect_error(specify_target(dgp, "(Intercept)", scale = "pp"), "absolute probability")
  precision <- specify_target(dgp, "(Intercept)", scale = "prob", criterion = "ci_width",
                              ci_width = 0.15, engine = "bayes_laplace")
  expect_equal(precision$true_effect_target, 0.65, tolerance = 1e-7)
  spec <- specify_target(dgp, "(Intercept)", engine = "bayes_laplace")
  validation <- validate_design_point(spec, 36, 12, nsim = 2L, seed = 432)
  expect_equal(validation$true_target, c(0, beta))
  expect_equal(validation$n_valid, c(2L, 2L))
})

test_that("design_crossed print and error messages guide allocation and balanced multiples", {
  des <- design_crossed(
    factors = list(study = c("restudy", "retrieval"), rounds = c("three", "four")),
    between_subj = "rounds"
  )
  combined <- paste(testthat::capture_messages(print(des)), collapse = "\n")
  expect_match(combined, "cell_probs", fixed = TRUE)
  expect_match(combined, "restudy_three", fixed = TRUE)
  expect_match(combined, "n_subj = 4k, n_item = 2m", fixed = TRUE)
  expect_match(combined, "each participant sees each item once", fixed = TRUE)

  expect_error(
    generate_design_schedule(des, n_subj = 10L, n_item = 10L),
    "Nearest balanced values: 8 or 12"
  )
  expect_error(
    generate_design_schedule(des, n_subj = 12L, n_item = 9L),
    "Nearest balanced values: 8 or 10"
  )
  expect_error(
    design_crossed(list(item = c("a", "b"))),
    "reserved column names"
  )
})


