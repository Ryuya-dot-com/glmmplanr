test_that("GLMM tabular export retains decisions and omits nested replication records", {
  verified <- data.frame(n_subj = c(8L, 16L, 24L), decision_status = c("fail", "indeterminate", "pass"),
                         scenario = c("low", "high", "low"), power = c(0.25, 0.5, 1))
  verified$replications <- I(list(NULL, data.frame(est = 1:5), NULL))
  plan <- structure(list(verified = verified), class = "glmmplanr_result")
  tab <- as.data.frame(plan)
  expect_identical(tab, verified[setdiff(names(verified), "replications")])
  expect_identical(rownames(as.data.frame(plan, row.names = letters[1:3])), letters[1:3])
  path <- tempfile()
  utils::write.csv(plan, path, row.names = FALSE)
  expect_equal(utils::read.csv(path), tab)
})

test_that("Rasch presentation preserves unevaluable items and all planned studies", {
  p <- plan_rasch_fixed_bank(seq(-2, 2, length.out = 12), 8L, c(difficulty = .35),
                            items_per_subject = 2L, item_seed = 1L, nsim = 4L)
  expect_identical(as.data.frame(p), p$summary)
  expect_identical(summary(p), p$summary)
  expect_identical(sum(p$item_summary$n_scheduled), 16L)
  unused <- p$item_summary$n_scheduled == 0
  expect_true(any(unused))
  expect_true(all(p$item_summary$n_finite[unused] == 0))
  expect_identical(as.data.frame(p)$n_met, c(0, 0))
  expect_false(any(grepl("replications", capture.output(print(p)))))
  expect_match(report_plan(p, FALSE), "0/4 studies met")
  expect_match(report_plan(p, FALSE), "No tested design has a confirmed pass")
  expect_error(report_plan(p, lang = "en"), "unused argument")
  path <- tempfile()
  saveRDS(p, path)
  expect_identical(summary(readRDS(path)), as.data.frame(p))
  expect_s3_class(as_gt(p), "gt_tbl")
  expect_s3_class(ggplot2::ggplot_build(plot(p)), "ggplot_built")
  fig <- plot(p, type = "items", metric = "bias")
  counts <- subset(fig$data, measure == "Finite estimates (studies)")
  expect_identical(counts$value, p$item_summary$n_finite)
  expect_identical(sum(is.finite(fig$data$value[fig$data$measure == "Bias (logits)"])),
                   sum(p$item_summary$n_finite > 0))
})

test_that("schedule components join indirect overlaps and retain isolated items", {
  tr <- data.frame(subj_id = c(1, 1, 2, 2, 3, 3, 4), item_id = c(1L,3L,2L,4L,3L,4L,5L))
  expect_identical(item_components(tr, 6L), c(1L,1L,1L,1L,5L,6L))
  des <- design_crossed(list(a = c("x", "y")), allocation = "crossed")
  out <- generate_design_schedule(des, 4L, 6L, diagnostics = TRUE)
  expect_equal(out$item_exposure$n_scheduled, rep(8, 6))
  expect_identical(out$item_exposure$component, rep(1L, 6))
})

test_that("marginal 2x2 contrasts and analysis formulas reproduce the design coding", {
  d <- design_testlet(list(A = c("a", "b"), B = c("c", "d")),
                      testlet_size = 2, between_subj = "B")
  cells <- cell_probs_from_diffs(d, .65, c(A = 8, B = 4, "A:B" = 4))
  g <- calibrate_dgp(d, cells, chance = 1/3, subj_sd = .8, item_sd = .6, testlet_sd = .2)
  expect_equal(check_dgp(g, FALSE)$marginal_effects$effects_pp, c(8, 4, 4), tolerance = 1e-6)
  expect_identical(coef(g), g$beta)
  tr <- with_rng_stream(15L, simulate_dgp_data(g, n_subj = 24L, n_item = 16L))
  for (re in c("maximal", "uncorrelated", "intercept_only")) {
    sp <- specify_target(g, "A", analysis_re = re)
    parts <- lme4::glFormula(formula(sp), tr, family = family(sp))
    expected_X <- d$X_cells[tr$cell_id, , drop = FALSE]
    expect_identical(colnames(parts$X), colnames(expected_X))
    expect_equal(as.numeric(parts$X), as.numeric(expected_X))
    re_cols <- function(z) {
      cols <- if (re == "intercept_only") 1L else seq_len(ncol(z))
      z[tr$cell_id, cols, drop = FALSE]
    }
    numeric_fit <- build_numeric_glmm_df(tr, expected_X, re_cols(d$Z_subj_cells),
      re_cols(d$Z_item_cells), re_cols(d$Z_testlet_cells), as.integer(re == "maximal"))
    reference <- lme4::glFormula(numeric_fit$formula, numeric_fit$data, family = family(sp))
    expect_equal(as.numeric(parts$reTrms$Zt), as.numeric(reference$reTrms$Zt))
    expect_equal(family(sp)$linkinv(0), 2/3)
    terms <- parts$reTrms$cnms
    expected <- switch(re, maximal = c(2L, 4L, 4L), intercept_only = rep(1L, 3L),
                       uncorrelated = rep(1L, 10L))
    expect_equal(sort(unname(lengths(terms))), sort(expected))
  }
  expect_identical(family(specify_target(g, "A", analysis_link = "logit"))$link, "logit")
})
