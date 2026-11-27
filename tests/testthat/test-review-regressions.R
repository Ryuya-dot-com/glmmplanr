test_that("slope sensitivity requires a positive slope SD", {
  d <- design_crossed(list(condition=c('a','b')))
  g <- calibrate_dgp(d, c(a=.4,b=.6), subj_sd=.4, item_sd=.3)
  sp <- specify_target(g, 'condition', prob_se_method='joint_delta')
  plan <- optimize_sample_size(sp, 16L, 8L, method='grid', nsim_final=3L,
    verbose=FALSE, control=optimize_control(nsim_sensitivity=2L))
  expect_null(plan$sensitivity)
  expect_false(has_random_slopes(g))
  g2 <- calibrate_dgp(d, c(a=.4,b=.6), subj_sd=c(.4,.2), item_sd=c(.3,0))
  expect_true(has_random_slopes(g2))
  td <- design_testlet(list(condition=c('a','b')), testlet_size=2L)
  tg <- calibrate_dgp(td, c(a=.4,b=.6), subj_sd=.4, item_sd=.3,
    control=calibrate_control(testlet_sd=c(.3,.2)))
  expect_true(has_random_slopes(tg))
})

test_that("custom data validates IDs and recovers the scheduled passages", {
  d <- design_testlet(list(condition=c('a','b')), testlet_size=2L)
  g <- calibrate_dgp(d, c(a=.4,b=.6), subj_sd=.4, item_sd=.3)
  sp <- specify_target(g,'condition',prob_se_method='joint_delta')
  schedule <- generate_design_schedule(d,8L,8L)
  tr <- schedule$trials[c('subj_id','item_id','cell_id')]
  tr$y <- rep(0:1, length.out=nrow(tr))
  fixed <- validate_custom_data(tr, schedule)
  expect_identical(fixed$testlet_id, schedule$trials$testlet_id)
  expect_identical(as.integer(fixed$testlet), fixed$testlet_id)
  ev <- evaluate_design_point(sp,8L,8L,nsim=1L,stream_seeds=15L,
    data_fn=function(schedule,spec) tr)
  expect_equal(ev$n_trials_obs_mean,nrow(tr))
  for (col in c('subj_id','item_id','cell_id')) for(value in c(0,-1,1.5,NA_real_,Inf,100)) {
    bad <- tr; bad[[col]][1] <- value
    expect_error(validate_custom_data(bad,schedule), 'integer IDs')
  }
  bad <- tr; bad$testlet_id <- rep(1L,nrow(tr))
  expect_error(validate_custom_data(bad,schedule), 'item-to-passage')
  subset <- fixed[fixed$testlet_id==3L, ]; subset$testlet_id <- NULL
  expect_true(all(validate_custom_data(subset,schedule)$testlet_id==3L))
  subset$testlet <- factor('1')
  expect_error(validate_custom_data(subset,schedule), 'item-to-passage')
})

test_that("absolute probability targets have precision but no zero-difference power", {
  d <- design_crossed(list(condition=c('a','b')))
  g <- calibrate_dgp(d,c(a=.4,b=.6),subj_sd=.3,item_sd=.3)
  for (criterion in c('power','both','joint')) {
    expect_error(specify_target(g,'condition',scale='pp',criterion=criterion,
      ci_width=80,control=target_control(contrast_weights=c(.5,.5))),'absolute probability')
  }
  sp <- specify_target(g,'condition',scale='pp',criterion='ci_width',
    ci_width=80,prob_se_method='joint_delta',
    control=target_control(ci_assurance=.3,contrast_weights=c(.5,.5)))
  for (ok in c(TRUE,FALSE)) {
    ev <- summarize_mc_metrics(rep(0,10),rep(.1,10),rep(.5,10),rep(.01,10),
      rep(NA,10),rep(ok,10),spec=sp,n_subj=16L,n_item=8L,nsim=10L,engine_use=sp$engine)
    expect_true(all(is.na(ev[c('power','effective_power','conditional_power','power_mc_lo','joint_prob','k_joint')])))
    expect_equal(ev$aipe_prob,as.numeric(ok))
  }
  plan <- optimize_sample_size(sp,16L,8L,method='grid',nsim_final=4L,
    verbose=FALSE,control=optimize_control(nsim_sensitivity=0L))
  expect_true(plan$feasible)
  printed <- paste(capture_messages(capture.output(print(plan))),collapse='\n')
  expect_match(printed,'Target accuracy')
  expect_false(grepl('Power [(]|H1 Coverage|Detection and width',printed))
  report <- report_plan(plan,FALSE)
  expect_match(report,'accuracy')
  expect_false(grepl('power =|accuracy difference',report))
  old_opt <- options(sass.cache = FALSE)
  on.exit(options(old_opt), add = TRUE)
  tab <- gt::as_raw_html(as_gt(plan))
  expect_match(tab,'absolute probability')
  expect_false(grepl('Effective Power|H1 Coverage',tab))
  searched <- optimize_sample_size(sp,c(8L,14L),8L,method='surrogate_gp',
    nsim_final=4L,verbose=FALSE,control=optimize_control(nsim_init=3L,
      n_eval_init=2L,n_eval_adapt=0L,n_verify=2L,nsim_sensitivity=0L))
  expect_true(all(is.na(searched$verified$power)))
  expect_true(any(is.finite(searched$verified$aipe_prob)))
})

# Deterministic rates straddle the adjusted Wilson boundary. These tests inspect
# the optimizer's control flow without relying on noisy model fitting.
test_that("multiplicity is applied before top-up and follow-up metadata is refreshed", {
  d <- design_crossed(list(condition=c('a','b')))
  g <- calibrate_dgp(d,c(a=.4,b=.6),subj_sd=.3,item_sd=.3)
  sp <- specify_target(g,'condition',power_min=.8)
  calls <- list()
  search_probe <- FALSE
  stage1_rate <- .9
  fake_batch <- function(spec,grid_sub,point_bases,nsim,engine_use,fast_stage,
      workers,parallel_backend,keep_raw=FALSE,allow_curtail=FALSE,stage_label='',verbose=FALSE,...) {
    calls[[length(calls)+1L]] <<- list(stage=stage_label,conf=spec$mc_conf_level,n=nsim,subjects=grid_sub$n_subj)
    k <- if(nsim==100L) 89L else floor(.99*nsim)
    rows <- lapply(seq_len(nrow(grid_sub)),function(i) {
      if (search_probe) k <- if (fast_stage) floor(stage1_rate*nsim) else if (grid_sub$n_subj[i]==14L) nsim else 88L
      summarize_mc_metrics(
      c(rep(1,k),rep(0,nsim-k)),rep(.1,nsim),rep(.1,nsim),rep(.01,nsim),
      rep(NA,nsim),rep(TRUE,nsim),spec=spec,n_subj=grid_sub$n_subj[i],n_item=grid_sub$n_item[i],nsim=nsim,engine_use=engine_use)
    })
    list(summary=do.call(rbind,rows),raw=rep(list(NULL),nrow(grid_sub)))
  }
  local_mocked_bindings(evaluate_grid_batch=fake_batch,.package='glmmplanr')
  plan <- optimize_sample_size(sp,c(8L,10L),8L,method='grid',nsim_final=100L,
    verbose=FALSE,control=optimize_control(mc_adjustment='bonferroni',nsim_topup=200L,nsim_sensitivity=0L))
  expect_true(any(vapply(calls,function(x) x$stage=='Stage 3 Top-Up',logical(1))))
  expect_equal(plan$verified$nsim_topup_used,rep(200L,2L))
  expect_true(all(plan$verified$decision_initial=='indeterminate'))
  expect_true(all(plan$verified$decision_status=='pass'))
  expect_equal(calls[[1]]$conf,.9875)
  plan$verified$decision_status <- 'indeterminate'
  expect_warning(more <- topup_verification(plan,nsim_extra=150L,topup_policy='independent',verbose=FALSE), 'sequential top-up round 2')
  expect_equal(more$verified$nsim,rep(150L,2L))
  expect_equal(more$verified$nsim_topup_used,rep(150L,2L))
  expect_equal(more$verified$topup_only_power,more$verified$power)
  calls <- list()
  sc <- optimize_sample_size(sp,c(8L,10L),8L,method='grid',nsim_final=100L,
    scenarios=list(base=g,other=g),standard_scenario='base',verbose=FALSE,
    control=optimize_control(mc_adjustment='bonferroni',nsim_topup=200L,nsim_sensitivity=0L))
  expect_true(all(sc$robust$scenario_evaluations$nsim_topup_used==200L))
  conf <- vapply(Filter(function(x) startsWith(x$stage,'Scenario '),calls),`[[`,numeric(1),'conf')
  expect_equal(conf,rep(.99375,2L))
  search_probe <- TRUE
  local_mocked_bindings(select_stage2_candidates=function(...) 1:2,.package='glmmplanr')
  for (rate in c(.9,.2)) {
    calls <- list(); stage1_rate <- rate
    searched <- optimize_sample_size(sp,c(8L,14L),8L,method='surrogate_gp',
      nsim_final=100L,verbose=FALSE,control=optimize_control(nsim_init=100L,
        n_eval_init=2L,n_eval_adapt=0L,n_verify=2L,mc_adjustment='bonferroni',nsim_sensitivity=0L))
    expect_true(searched$feasible)
    expect_equal(searched$optimal$n_subj,14L)
    expected <- if (rate==.9) 'Stage 2 Step-Up' else 'Stage 2 Ceiling'
    expect_true(any(vapply(calls,function(x) x$stage==expected,logical(1))),
      info=paste(expected,paste(vapply(calls,`[[`,character(1),'stage'),collapse=', ')))
  }
})

test_that("fixed-bank MML rejects duplicates and mislabeled constant items", {
  tr <- expand.grid(subj_id=1:4,item_id=1:2)
  tr$y <- c(0L,1L,1L,0L,1L,0L,0L,1L)
  fit <- function(x) fit_rasch_fixed_mml(x,4L,2L,c(TRUE,TRUE))
  expect_true(is.na(fit(tr)$error))
  for (row in c(1L,2L)) expect_match(fit(rbind(tr,tr[row,]))$error,'Duplicate')
  for (y in 0:1) {
    bad <- tr; bad$y[bad$item_id==1L] <- y
    expect_match(fit(bad)$error,'both response categories')
  }
})

test_that("missing-condition partial schedules fail identically on cached and data-frame paths", {
  des <- design_crossed(
    list(condition = c("a", "b")),
    allocation = "between_item",
    items_per_subject = 1L,
    item_seed = 8L
  )
  sched <- generate_design_schedule(des, n_subj = 4L, n_item = 6L, diagnostics = TRUE)
  expect_false(sched$identifiable)
  expect_true(length(sched$missing_cells) >= 1L)
  expect_true(any(sched$cell_exposure$n_scheduled == 0L))

  dgp <- calibrate_dgp(des, c(a = 0.4, b = 0.6), subj_sd = 0.3, item_sd = 0.4)
  spec <- specify_target(
    dgp, "condition", scale = "pp", criterion = "ci_width",
    ci_width = 20, control = target_control(prob_boot = 9L)
  )

  ev_fast <- evaluate_design_point(spec, 4L, 6L, nsim = 10L, stream_seeds = 41L)
  ev_slow <- evaluate_design_point(
    spec, 4L, 6L, nsim = 10L, stream_seeds = 41L,
    data_fn = function(schedule, spec) simulate_dgp_data(spec$dgp, schedule = schedule)
  )
  expect_equal(ev_fast$n_valid, 0L)
  expect_equal(ev_slow$n_valid, 0L)
  expect_identical(ev_fast$decision_status, "fail")
  expect_identical(ev_slow$decision_status, "fail")
  expect_match(ev_fast$top_error_reason, "Observed trials omit required design cells")
  expect_identical(ev_fast$top_error_reason, ev_slow$top_error_reason)
})

test_that("custom data_fn truth and validate_design_point MC intervals stay aligned", {
  des <- design_crossed(list(condition = c("a", "b")))
  dgp15 <- calibrate_dgp(des, c(a = 0.45, b = 0.60), subj_sd = 0.3, item_sd = 0.3)
  dgp70 <- calibrate_dgp(des, c(a = 0.15, b = 0.85), subj_sd = 0.3, item_sd = 0.3)
  sp <- specify_target(
    dgp15, "condition", scale = "pp", criterion = "joint",
    ci_width = 35, prob_se_method = "joint_delta"
  )

  fn70 <- function(schedule, spec) simulate_dgp_data(dgp70, schedule = schedule)

  ev_unverified <- evaluate_design_point(sp, 24L, 12L, nsim = 5L, stream_seeds = 71L, data_fn = fn70)
  expect_true(all(is.na(ev_unverified[c("coverage", "prob_coverage", "bias", "rmse", "prob_diff_true_pp")])))

  ev_verified <- evaluate_design_point(
    sp, 24L, 12L, nsim = 5L, stream_seeds = 71L, data_fn = fn70, true_target = 70
  )
  expect_equal(ev_verified$prob_diff_true_pp, 70)
  expect_true(is.finite(ev_verified$coverage))

  expect_message(
    val_unverified <- validate_design_point(sp, 24L, 12L, nsim = 5L, seed = 72L, data_fn = fn70),
    "unverified"
  )
  expect_true(all(is.na(val_unverified$true_target)))
  expect_true(all(is.na(val_unverified$coverage)))
  expect_false("Type I error" %in% as.character(plot(val_unverified)$data$panel))

  val_known <- validate_design_point(
    sp, 24L, 12L, nsim = 5L, seed = 72L, data_fn = fn70,
    true_target = c(null = 70, planned = 70)
  )
  expect_equal(val_known$true_target, c(70, 70))
  expect_true(all(is.finite(val_known$coverage)))
  expect_false("Type I error" %in% as.character(plot(val_known)$data$panel))

  val_builtin <- validate_design_point(sp, 24L, 12L, nsim = 6L, seed = 73L)
  expect_true(all(c(
    "decision_status", "width_attainment_lo", "width_attainment_hi",
    "joint_attainment_lo", "joint_attainment_hi"
  ) %in% names(val_builtin)))
  expect_true(is.na(val_builtin$decision_status[1L]))
  expect_true(val_builtin$decision_status[2L] %in% c("pass", "indeterminate", "fail"))
  p_j <- plot(val_builtin)
  expect_true("Joint attainment" %in% as.character(p_j$data$panel))
  expect_warning(ggplot2::ggplot_build(p_j), NA)
})

test_that("interim curtailment only drops valid_rate when remaining reps cannot reach floor", {
  des <- design_crossed(list(condition = c("a", "b")))
  dgp <- calibrate_dgp(des, c(a = 0.35, b = 0.70), subj_sd = 0.3, item_sd = 0.3)
  sp <- specify_target(
    dgp, "condition", criterion = "power", power_min = 0.80,
    control = target_control(failure_policy = "conditional", min_valid_rate = 0.90)
  )

  counters <- new.env(parent = emptyenv())
  fake_rep <- function(spec, schedule, ...) {
    key <- as.character(schedule$n_subj)
    idx <- (if (is.null(counters[[key]])) 0L else counters[[key]]) + 1L
    counters[[key]] <- idx
    ok <- switch(
      key,
      "6" = idx <= 70L,
      "8" = !(idx >= 90L && idx <= 100L),
      "10" = TRUE
    )
    prim <- if (ok) {
      list(est = 2.0, se = 0.1, prob_est = 0.35, prob_se = 0.02,
           singular = NA, converged = TRUE, status = "converged", error_reason = "")
    } else {
      make_failed_fit("non_pd_hessian", "non_pd_hessian", singular = NA)
    }
    list(primary = prim, sing_init = NA, reduced = NULL, nobs = nrow(schedule$trials),
         patterns = empty_pattern_totals)
  }

  local_mocked_bindings(execute_single_replication = fake_rep, .package = "glmmplanr")
  rm(list = ls(counters), envir = counters)
  res_on <- optimize_sample_size(
    sp, c(6L, 8L, 10L), 4L, prefer = "fewer_subj", method = "grid",
    nsim_final = 200L, workers = 1L, verbose = FALSE,
    control = optimize_control(early_stop = TRUE, nsim_sensitivity = 0L)
  )
  rm(list = ls(counters), envir = counters)
  res_off <- optimize_sample_size(
    sp, c(6L, 8L, 10L), 4L, prefer = "fewer_subj", method = "grid",
    nsim_final = 200L, workers = 1L, verbose = FALSE,
    control = optimize_control(early_stop = FALSE, nsim_sensitivity = 0L)
  )

  expect_true(res_on$verified$curtailed[res_on$verified$n_subj == 6L])
  expect_false(res_on$verified$curtailed[res_on$verified$n_subj == 8L])
  expect_equal(res_on$verified$n_valid[res_on$verified$n_subj == 8L], 189L)
  expect_equal(res_on$optimal$n_subj, 8L)
  expect_equal(res_on$optimal$n_subj, res_off$optimal$n_subj)
})

test_that("preference rules are strictly lexicographic and scenario trial counts stay distinct", {
  des <- design_crossed(list(condition = c("a", "b")))
  dgp <- calibrate_dgp(des, c(a = 0.4, b = 0.7), subj_sd = 0.3, item_sd = 0.3)
  sp <- specify_target(dgp, "condition", prob_se_method = "joint_delta")

  grid_wide <- expand.grid(n_subj = c(100L, 102L, 500L), n_item = c(8L, 1000L), KEEP.OUT.ATTRS = FALSE)
  grid_wide$n_trials_scheduled <- grid_wide$n_subj * grid_wide$n_item
  sc_fs <- compute_preference_score(grid_wide, "fewer_subj")
  expect_lt(
    sc_fs[grid_wide$n_subj == 100L & grid_wide$n_item == 1000L],
    sc_fs[grid_wide$n_subj == 102L & grid_wide$n_item == 8L]
  )

  ver <- data.frame(
    n_subj = c(100L, 102L),
    n_item = c(1000L, 8L),
    n_trials_scheduled = c(100000L, 816L),
    power = c(0.90, 0.92),
    ci_width_mean = c(0.3, 0.3),
    preference_score = c(0.9, 0.1),
    criterion_met = c(TRUE, TRUE),
    decision_status = c("pass", "pass")
  )
  expect_equal(select_recommended_from_verified(ver, sp, prefer = "fewer_subj")$optimal$n_subj, 100L)
  expect_equal(select_recommended_from_verified(ver, sp, prefer = "fewer_items")$optimal$n_item, 8L)
  expect_equal(select_recommended_from_verified(ver, sp, prefer = "fewer_trials")$optimal$n_trials_scheduled, 816L)

  des_part <- design_crossed(list(condition = c("a", "b")), items_per_subject = 4L)
  dgp_part <- calibrate_dgp(des_part, c(a = 0.20, b = 0.90), subj_sd = 0.2, item_sd = 0.2)
  dgp_full <- calibrate_dgp(des, c(a = 0.20, b = 0.90), subj_sd = 0.2, item_sd = 0.2)
  sp_part <- specify_target(dgp_part, "condition", power_min = 0.50, prob_se_method = "joint_delta")

  batch_calls <- 0L
  orig_batch <- evaluate_grid_batch
  counting_batch <- function(...) {
    batch_calls <<- batch_calls + 1L
    orig_batch(...)
  }
  local_mocked_bindings(evaluate_grid_batch = counting_batch, .package = "glmmplanr")

  plan_sc <- optimize_sample_size(
    sp_part, c(8L, 10L), 8L,
    scenarios = list(partial = dgp_part, full = dgp_full),
    standard_scenario = "partial",
    method = "grid", nsim_final = 10L, workers = 1L, verbose = FALSE,
    control = optimize_control(nsim_sensitivity = 0L)
  )
  # Full-grid multi-scenario skips duplicate exploratory pass: exactly 1 batch per scenario
  expect_equal(batch_calls, 2L)
  expect_true(isTRUE(plan_sc$feasible))
  expect_true(isTRUE(plan_sc$robust$schedules_differ))
  expect_identical(names(plan_sc$robust$optimal_standard), names(plan_sc$robust$optimal_robust))
  expect_true("scheduled_trials_by_scenario" %in% names(plan_sc$robust$comparison))
  expect_match(plan_sc$robust$comparison$scheduled_trials_by_scenario[2L], "partial = ")
  expect_match(plan_sc$robust$comparison$scheduled_trials_by_scenario[2L], "full = ")
  rep_sc <- report_plan(plan_sc, print_output = FALSE)
  expect_match(rep_sc, "partial = ")
  expect_match(rep_sc, "full = ")
})

test_that("topup_verification synchronizes evaluated and surface and absolute plots avoid 50% power", {
  des <- design_crossed(list(condition = c("a", "b")))
  dgp <- calibrate_dgp(des, c(a = 0.4, b = 0.7), subj_sd = 0.3, item_sd = 0.3)
  sp <- specify_target(dgp, "condition", scale = "pp", criterion = "ci_width",
                       ci_width = 30, prob_se_method = "joint_delta")

  round_flag <- 1L
  fake_batch <- function(spec, grid_sub, point_bases, nsim, engine_use, fast_stage, ...) {
    w_se <- if (round_flag == 1L) 0.075 else 0.040
    k_ok <- if (round_flag == 1L) 7L else nsim
    rows <- lapply(seq_len(nrow(grid_sub)), function(i) {
      se_v <- c(rep(w_se, k_ok), rep(0.20, nsim - k_ok))
      summarize_mc_metrics(
        rep(0.8, nsim), rep(0.1, nsim), rep(0.30, nsim), se_v,
        rep(NA, nsim), rep(TRUE, nsim),
        spec = spec, n_subj = grid_sub$n_subj[i], n_item = grid_sub$n_item[i],
        nsim = nsim, engine_use = engine_use
      )
    })
    list(summary = do.call(rbind, rows), raw = rep(list(NULL), nrow(grid_sub)))
  }
  local_mocked_bindings(evaluate_grid_batch = fake_batch, .package = "glmmplanr")

  plan8 <- optimize_sample_size(
    sp, c(8L, 12L), c(8L, 12L), method = "grid", nsim_final = 8L,
    workers = 1L, verbose = FALSE, control = optimize_control(nsim_sensitivity = 0L)
  )
  expect_true(all(plan8$verified$decision_status == "indeterminate"))
  old_width <- plan8$surface$pred_prob_ci_width_pp[1L]

  round_flag <- 2L
  plan30 <- topup_verification(plan8, nsim_extra = 30L, topup_policy = "independent", verbose = FALSE)
  tf <- tempfile(fileext = ".rds")
  on.exit(unlink(tf), add = TRUE)
  saveRDS(plan30, tf)
  reloaded <- readRDS(tf)

  expect_true(all(reloaded$verified$decision_status == "pass"))
  expect_identical(reloaded$evaluated$decision_status, reloaded$verified$decision_status)
  expect_equal(reloaded$evaluated$nsim, rep(30L, nrow(reloaded$verified)))
  expect_lt(reloaded$surface$pred_prob_ci_width_pp[1L], old_width - 5)

  # Absolute probability target: pred_power is NA and default plot shows CI width
  des_one <- design_crossed(list())
  dgp_one <- calibrate_dgp(des_one, 0.65, subj_sd = 0.3, item_sd = 0.4)
  sp_one <- specify_target(
    dgp_one, "(Intercept)", scale = "pp", criterion = "ci_width",
    ci_width = 25, prob_se_method = "joint_delta"
  )
  rm_mock <- local_mocked_bindings(evaluate_grid_batch = evaluate_grid_batch, .package = "glmmplanr")
  plan_one <- optimize_sample_size(
    sp_one, c(12L, 16L), c(8L, 12L), method = "grid", nsim_final = 4L,
    workers = 1L, verbose = FALSE, control = optimize_control(nsim_sensitivity = 0L)
  )
  expect_true(all(is.na(plan_one$surface$pred_power)))
  p_def <- plot(plan_one, show_annotations = TRUE)
  expect_warning(ggplot2::ggplot_build(p_def), NA)
  expect_false(grepl("Power = NA", p_def$labels$subtitle, fixed = TRUE))
  expect_error(plot(plan_one, metric = "power"), "undefined")
  p_aipe <- plot(plan_one, metric = "aipe_prob", show_annotations = TRUE)
  expect_warning(ggplot2::ggplot_build(p_aipe), NA)
  if (requireNamespace("plotly", quietly = TRUE)) {
    expect_warning(plotly::plotly_build(plot(plan_one, type = "surface3d")), NA)
  }
})

test_that("factor responses and non-integer or non-positive replication counts are rejected", {
  des <- design_crossed(list(condition = c("a", "b")))
  dgp <- calibrate_dgp(des, c(a = 0.4, b = 0.6), subj_sd = 0.3, item_sd = 0.3)
  sp <- specify_target(dgp, "condition", prob_se_method = "joint_delta")
  sched <- generate_design_schedule(des, 8L, 4L)

  tr_fac <- sched$trials
  tr_fac$y <- factor(rep("0", nrow(tr_fac)), levels = c("0", "1"))
  expect_error(validate_custom_data(tr_fac, sched), "numeric/integer 0/1")
  expect_error(
    evaluate_design_point(sp, 8L, 4L, nsim = 2L, data_fn = function(schedule, spec) tr_fac),
    "numeric/integer 0/1"
  )

  for (bad_n in list(0, -2, 1.8, NA_real_, Inf)) {
    expect_error(optimize_sample_size(sp, 8L, 4L, nsim_final = bad_n, verbose = FALSE), "positive integer")
    expect_error(evaluate_design_point(sp, 8L, 4L, nsim = bad_n), "positive integer")
    expect_error(validate_design_point(sp, 8L, 4L, nsim = bad_n), "positive integer")
  }
  plan <- optimize_sample_size(
    sp, 8L, 4L, method = "grid", nsim_final = 3L, verbose = FALSE,
    control = optimize_control(nsim_sensitivity = 0L)
  )
  for (bad_n in list(0, 1.8, NA_real_)) {
    expect_error(topup_verification(plan, nsim_extra = bad_n, verbose = FALSE), "positive integer")
  }
  expect_error(optimize_sample_size(sp, 8L, 4L, nsim_final = 3L, workers = 0, verbose = FALSE), "positive integer")
})

test_that("fixed-bank Rasch MML invalidates replications with unidentifiable person variance", {
  bank4 <- stats::setNames(seq(-0.6, 0.6, length.out = 4L), paste0("item", 1:4))
  res1 <- plan_rasch_fixed_bank(
    item_difficulties = bank4, n_subj = 400L,
    max_error = c(person_sd = 0.01), person_sd = 0.8,
    items_per_subject = 1L, item_seed = 11L,
    method = "fixed_item", engine = "mml",
    nsim = 32L, seed = 901L, verbose = FALSE
  )
  psd_row <- res1$summary[res1$summary$target == "person_sd", , drop = FALSE]
  jnt_row <- res1$summary[res1$summary$target == "joint", , drop = FALSE]
  expect_equal(psd_row$n_evaluable, 0L)
  expect_equal(psd_row$attainment, 0.0)
  expect_identical(psd_row$decision, "fail")
  expect_identical(jnt_row$decision, "fail")
  expect_true(is.na(res1$recommended_n_subj))

  # Even under complete presentation of 2 items, if one item is constant (all 0),
  # only 1 estimable item remains so no participant has >= 2 estimable responses.
  des2 <- design_crossed(list())
  bank2 <- c(item1 = 0.0, item2 = 3.0)
  sched2 <- generate_design_schedule(des2, n_subj = 12L, n_item = 2L)
  tr2 <- sched2$trials
  tr2$y <- ifelse(tr2$item_id == 1L, rep(c(0L, 1L), length.out = nrow(tr2)), 0L)
  sim2 <- list(trials = tr2, theta = rep(0, 12L), n_subj = 12L, n_item = 2L)
  fit_const <- eval_fixed_bank_replication(sim2, bank2, person_sd = 0.8, method = "fixed_item", engine = "mml")
  expect_true(grepl("two or more responses on estimable items", fit_const$metrics$fit_error, fixed = TRUE))
  expect_true(is.na(fit_const$metrics$person_sd))
})

test_that("unattainable zero probability contrasts and custom-truth overrides without data_fn are rejected", {
  des <- design_crossed(list(condition = c("a", "b")))
  dgp_3afc <- calibrate_dgp(des, c(a = 0.7, b = 0.8), subj_sd = 0.5, item_sd = 0.4, chance = 1 / 3)
  w_bad <- c(a = -1, b = 0.2)

  for (crit in c("power", "both", "joint")) {
    expect_error(
      specify_target(
        dgp_3afc, "condition",
        scale = "pp", criterion = crit, ci_width = 20, prob_se_method = "joint_delta",
        control = target_control(contrast_weights = w_bad)
      ),
      "outside the attainable probability contrast"
    )
  }

  sp_width <- specify_target(
    dgp_3afc, "condition",
    scale = "pp", criterion = "ci_width", ci_width = 35, prob_se_method = "joint_delta",
    control = target_control(contrast_weights = w_bad)
  )
  expect_false(sp_width$has_null_target)
  set.seed(41L)
  ev_w <- evaluate_design_point(sp_width, 16L, 12L, nsim = 3L)
  expect_true(is.na(ev_w$power))
  val_w <- validate_design_point(sp_width, 16L, 12L, nsim = 3L, seed = 42L)
  expect_equal(nrow(val_w), 1L)
  expect_identical(val_w$condition, "planned")

  # Scenario override where base DGP allows zero (-1 * 0.1 + 0.2 * 0.5 = 0) but scenario chance = 1/3 does not
  dgp_ch0 <- calibrate_dgp(des, c(a = 0.7, b = 0.8), subj_sd = 0.5, item_sd = 0.4, chance = 0)
  sp_ch0 <- specify_target(
    dgp_ch0, "condition",
    scale = "pp", criterion = "power", prob_se_method = "joint_delta",
    control = target_control(contrast_weights = w_bad)
  )
  expect_true(sp_ch0$has_null_target)
  expect_error(
    optimize_sample_size(
      sp_ch0, n_subj_range = c(12L, 16L), n_item_range = 8L,
      scenarios = list(base = dgp_ch0, three_choice = dgp_3afc),
      method = "grid", nsim_final = 2L, verbose = FALSE,
      control = optimize_control(nsim_sensitivity = 0L)
    ),
    "outside the attainable probability contrast"
  )

  # Disallow true_target when data_fn is NULL
  sp_ok <- specify_target(dgp_ch0, "condition", scale = "pp", prob_se_method = "joint_delta")
  expect_error(evaluate_design_point(sp_ok, 12L, 8L, nsim = 2L, true_target = 0), "data_fn")
  expect_error(validate_design_point(sp_ok, 12L, 8L, nsim = 2L, true_target = c(null = 0, planned = 0)), "data_fn")
})

test_that("surrogate search enforces min_valid_rate under conditional failure policy and resolves stage1_re per scenario", {
  des <- design_crossed(list(condition = c("a", "b")))
  dgp_int <- calibrate_dgp(des, c(a = 0.45, b = 0.65), subj_sd = 0.5, item_sd = 0.4)
  dgp_slp <- calibrate_dgp(des, c(a = 0.45, b = 0.65), subj_sd = c(0.5, 0.3), item_sd = c(0.4, 0.25))

  sp_cond <- specify_target(
    dgp_int, "condition", scale = "pp", criterion = "power", power_min = 0.80,
    prob_se_method = "joint_delta",
    control = target_control(failure_policy = "conditional", min_valid_rate = 0.90)
  )
  expect_identical(sp_cond$stage1_re, "intercept_only")

  real_eval <- evaluate_grid_batch
  seen_stage1_re <- list()
  local_mocked_bindings(
    evaluate_grid_batch = function(spec, grid_sub, point_bases, nsim, engine_use = spec$engine, fast_stage = FALSE, ...) {
      if (isTRUE(fast_stage)) {
        seen_stage1_re[[length(seen_stage1_re) + 1L]] <<- list(
          has_slopes = any(spec$dgp$tau_subj[-1L] > 0),
          stage1_re = spec$stage1_re
        )
      }
      base <- real_eval(spec, grid_sub[1L, , drop = FALSE], point_bases[1L], nsim = 2L,
                        engine_use = spec$engine, fast_stage = TRUE, workers = 1L,
                        parallel_backend = "serial", verbose = FALSE)$summary
      rows <- lapply(seq_len(nrow(grid_sub)), function(i) {
        r <- base
        ns <- grid_sub$n_subj[i]
        ni <- grid_sub$n_item[i]
        high_valid <- (ns >= 12L && ni >= 10L)
        vrate <- if (high_valid) 0.96 else 0.50
        nv <- max(1L, as.integer(round(nsim * vrate)))
        ksig <- if (high_valid) as.integer(round(nv * 0.96)) else nv
        r$n_subj <- ns
        r$n_item <- ni
        r$nsim <- as.integer(nsim)
        r$n_valid <- nv
        r$valid_rate <- nv / nsim
        r$nonconv_rate <- 1 - r$valid_rate
        r$conditional_power <- ksig / nv
        r$effective_power <- ksig / nsim
        r$effective_rope_prob <- ksig / nsim
        r$power <- r$conditional_power
        mc_p <- wilson_mc_ci(ksig, nv, conf_level = spec$mc_conf_level)
        r$power_mc_lo <- unname(mc_p["lo"])
        r$power_mc_hi <- unname(mc_p["hi"])
        valid_ok <- r$valid_rate >= spec$min_valid_rate
        r$criterion_met <- valid_ok && (r$power >= spec$power_min)
        r$decision_status <- compute_decision_status(
          spec, valid_ok, mc_p, c(lo = NA_real_, hi = NA_real_),
          c(lo = NA_real_, hi = NA_real_), c(lo = NA_real_, hi = NA_real_)
        )
        r$mc_conf_level_used <- spec$mc_conf_level
        r
      })
      list(summary = do.call(rbind, rows), raw = vector("list", nrow(grid_sub)))
    },
    .package = "glmmplanr"
  )

  plan <- optimize_sample_size(
    sp_cond,
    n_subj_range = c(8L, 12L),
    n_item_range = c(6L, 10L),
    prefer = "fewer_trials",
    method = "surrogate_gp",
    nsim_final = 50L,
    workers = 1L, verbose = FALSE,
    control = optimize_control(
      n_eval_init = 7L, n_eval_adapt = 1L, n_verify = 4L,
      nsim_init = 20L, nsim_sensitivity = 0L
    )
  )
  low_rows <- plan$surface$n_subj < 12L | plan$surface$n_item < 10L
  expect_false(any(plan$surface$pred_criterion_met[low_rows]))
  expect_true(plan$feasible)
  expect_equal(plan$optimal$n_subj, 12L)
  expect_equal(plan$optimal$n_item, 10L)
  expect_identical(plan$optimal$decision_status, "pass")

  seen_stage1_re <- list()
  plan_ms <- optimize_sample_size(
    sp_cond,
    n_subj_range = c(8L, 16L),
    n_item_range = c(8L, 12L),
    scenarios = list(base = dgp_int, slopes = dgp_slp),
    method = "surrogate_gp",
    nsim_final = 50L,
    workers = 1L, verbose = FALSE,
    control = optimize_control(
      n_eval_init = 3L, n_eval_adapt = 1L, n_verify = 2L,
      nsim_init = 20L, nsim_sensitivity = 0L
    )
  )
  expect_identical(plan_ms$robust$specs$base$stage1_re, "intercept_only")
  expect_identical(plan_ms$robust$specs$slopes$stage1_re, "maximal")
  slope_calls <- Filter(function(x) isTRUE(x$has_slopes), seen_stage1_re)
  expect_gt(length(slope_calls), 0L)
  expect_true(all(vapply(slope_calls, function(x) identical(x$stage1_re, "maximal"), logical(1L))))
})

test_that("report_plan preserves per-candidate Wilson levels after top-up and validate_design_point reports ROPE intervals", {
  des <- design_crossed(list(condition = c("a", "b")), items_per_subject = 6L, item_seed = 19L)
  dgp <- calibrate_dgp(des, c(a = 0.45, b = 0.65), subj_sd = 0.4, item_sd = 0.3)
  sp <- specify_target(dgp, "condition", scale = "pp", criterion = "power", prob_se_method = "joint_delta")

  real_eval <- evaluate_grid_batch
  call_no <- 0L
  local_mocked_bindings(
    evaluate_grid_batch = function(spec, grid_sub, point_bases, nsim, ...) {
      call_no <<- call_no + 1L
      base <- real_eval(spec, grid_sub[1L, , drop = FALSE], point_bases[1L], nsim = 2L,
                        engine_use = spec$engine, fast_stage = FALSE, workers = 1L,
                        parallel_backend = "serial", verbose = FALSE)$summary
      rows <- lapply(seq_len(nrow(grid_sub)), function(i) {
        r <- base
        ns <- grid_sub$n_subj[i]
        ni <- grid_sub$n_item[i]
        ksig <- if (call_no == 1L) {
          if (ns == 14L) as.integer(nsim) else as.integer(round(0.84 * nsim))
        } else {
          as.integer(nsim)
        }
        r$n_subj <- ns
        r$n_item <- ni
        r$nsim <- as.integer(nsim)
        r$n_valid <- as.integer(nsim)
        r$valid_rate <- 1.0
        r$nonconv_rate <- 0.0
        r$conditional_power <- ksig / nsim
        r$effective_power <- ksig / nsim
        r$power <- ksig / nsim
        mc_p <- wilson_mc_ci(ksig, nsim, conf_level = spec$mc_conf_level)
        r$power_mc_lo <- unname(mc_p["lo"])
        r$power_mc_hi <- unname(mc_p["hi"])
        r$criterion_met <- r$power >= spec$power_min
        r$decision_status <- compute_decision_status(
          spec, TRUE, mc_p, c(lo = NA_real_, hi = NA_real_),
          c(lo = NA_real_, hi = NA_real_), c(lo = NA_real_, hi = NA_real_)
        )
        r$mc_conf_level_used <- spec$mc_conf_level
        r
      })
      list(summary = do.call(rbind, rows), raw = vector("list", nrow(grid_sub)))
    },
    .package = "glmmplanr"
  )

  plan0 <- optimize_sample_size(
    sp, n_subj_range = c(12L, 14L), n_item_range = 12L,
    method = "grid", nsim_final = 40L,
    workers = 1L, verbose = FALSE,
    control = optimize_control(mc_adjustment = "bonferroni", nsim_sensitivity = 0L)
  )
  expect_equal(plan0$search_meta$mc_conf_level_effective, 0.975, tolerance = 1e-8)
  expect_identical(plan0$verified$decision_status, c("indeterminate", "pass"))

  plan_tu <- topup_verification(plan0, nsim_extra = 60L, topup_policy = "independent", verbose = FALSE)
  expect_equal(plan_tu$search_meta$mc_conf_level_effective, 0.975, tolerance = 1e-8)
  expect_equal(plan_tu$search_meta$topup_mc_conf_level, 0.99375, tolerance = 1e-8)
  # Now (12, 12) passed via top-up (99.38% Wilson CI) while (14, 12) kept its initial 97.50% Wilson CI
  expect_equal(plan_tu$optimal$n_subj, 12L)
  expect_equal(plan_tu$optimal$mc_conf_level_used, 0.99375, tolerance = 1e-8)
  rpt_tu <- report_plan(plan_tu, print_output = FALSE)
  expect_match(rpt_tu, "initial Stage-2 level 97.50%", fixed = TRUE)
  expect_match(rpt_tu, "latest Wilson level 99.38%", fixed = TRUE)
  expect_match(rpt_tu, "99.38% Wilson CI", fixed = TRUE)
  expect_match(rpt_tu, "condition on this fixed schedule", fixed = TRUE)

  # If the initially passing candidate (14, 12) remains optimal, report_plan cites 97.50% Wilson CI
  plan_keep <- plan_tu
  plan_keep$optimal <- plan_tu$verified[plan_tu$verified$n_subj == 14L, , drop = FALSE]
  rpt_keep <- report_plan(plan_keep, print_output = FALSE)
  expect_match(rpt_keep, "97.50% Wilson CI", fixed = TRUE)

  # ROPE validation intervals and plot
  rm_mock <- local_mocked_bindings(evaluate_grid_batch = real_eval, .package = "glmmplanr")
  sp_rope <- specify_target(
    dgp, "condition", scale = "pp", criterion = "rope",
    prob_se_method = "joint_delta",
    control = target_control(rope = c(-3, 3))
  )
  val_rope <- validate_design_point(sp_rope, 16L, 12L, nsim = 6L, seed = 88L)
  expect_true(all(c("rope_attainment", "rope_attainment_lo", "rope_attainment_hi") %in% names(val_rope)))
  expect_true(all(is.finite(val_rope$rope_attainment)))
  expect_true(all(val_rope$rope_attainment_lo <= val_rope$rope_attainment + 1e-12))
  expect_true(all(val_rope$rope_attainment_hi >= val_rope$rope_attainment - 1e-12))
  p_rope <- plot(val_rope)
  expect_warning(gb <- ggplot2::ggplot_build(p_rope), NA)
  expect_true("ROPE exclusion" %in% as.character(p_rope$data$panel))
})

test_that("observed-data cell separation check rejects separated bootstrap refits accepted by finite-only glmer", {
  skip_if_not_installed("lme4")
  des <- design_crossed(list(condition = c("a", "b")))
  sched <- generate_design_schedule(des, 8L, 4L)
  tr_sep <- sched$trials
  # Force complete separation in cell "a" (all 0) while cell "b" has mixed 0/1
  tr_sep$y <- ifelse(tr_sep$condition == "a", 0L, rep(c(0L, 1L), length.out = nrow(tr_sep)))
  cell_ok <- all(tapply(tr_sep$y, tr_sep$condition, function(v) any(v == 0L) && any(v == 1L)))
  expect_false(cell_ok)
  dgp <- calibrate_dgp(des, c(a = 0.3, b = 0.7), subj_sd = 0.4, item_sd = 0.3)
  sp_lme4 <- specify_target(dgp, "condition", scale = "pp", engine = "lme4", prob_se_method = "joint_delta")
  fit_res <- fit_single_dataset(tr_sep, sp_lme4, engine_use = "lme4")
  expect_false(isTRUE(fit_res$converged))
  expect_identical(fit_res$status, "separation")
})

test_that("observed-analysis vignette chunk executes on valid 54x30 data and catches separation with NA rows", {
  skip_if_not_installed("lme4")
  rmd_candidates <- c(
    "vignettes/planning-workflows.Rmd",
    file.path("..", "..", "vignettes", "planning-workflows.Rmd"),
    system.file("doc", "planning-workflows.Rmd", package = "glmmplanr")
  )
  rmd_candidates <- rmd_candidates[nzchar(rmd_candidates) & file.exists(rmd_candidates)]
  skip_if(length(rmd_candidates) == 0L, "planning-workflows.Rmd not found")

  lines <- readLines(rmd_candidates[1L], warn = FALSE)
  start_idx <- grep("^```\\{r observed-analysis", lines)
  expect_length(start_idx, 1L)
  end_candidates <- grep("^```\\s*$", lines)
  end_idx <- end_candidates[end_candidates > start_idx][1L]
  chunk_src <- paste(lines[(start_idx + 1L):(end_idx - 1L)], collapse = "\n")

  # 1. Execute the literal vignette chunk on a 54 x 30 synthetic dataset
  design <- design_crossed(list(study_method = c("restudy", "retrieval")))
  dgp <- calibrate_dgp(design, c(restudy = 0.42, retrieval = 0.50), subj_sd = 0.8, item_sd = 0.6)
  target <- specify_target(
    dgp, "study_method", scale = "pp", engine = "lme4",
    control = target_control(prob_boot = 12L)
  )
  set.seed(8201L)
  sim_df <- simulate_dgp_data(dgp, n_subj = 54L, n_item = 30L)
  names(sim_df)[names(sim_df) == "subj_id"] <- "subj"
  names(sim_df)[names(sim_df) == "item_id"] <- "item"

  env54 <- new.env(parent = globalenv())
  env54$design <- design
  env54$target <- target
  env54$responses <- sim_df
  out54 <- suppressMessages(eval(parse(text = chunk_src), envir = env54))
  expect_identical(typeof(env54$fit@optinfo$conv$opt), "double")
  expect_equal(env54$fit@optinfo$conv$opt, 0)
  expect_true(env54$valid_glmer_fit(env54$fit))
  expect_named(out54, c("estimate_pp", "se_pp", "ci_lo_pp", "ci_hi_pp", "width_pp"))
  expect_true(all(is.finite(out54)))
  expect_gt(out54[["width_pp"]], 0)

  # 2. Verify row-dropping alignment when responses contain NA rows on an 8 x 4 dataset
  dgp_ext <- calibrate_dgp(design, c(restudy = 0.12, retrieval = 0.88), subj_sd = 0.8, item_sd = 0.6)
  set.seed(104L)
  sim8 <- simulate_dgp_data(dgp_ext, n_subj = 8L, n_item = 4L)
  names(sim8)[names(sim8) == "subj_id"] <- "subj"
  names(sim8)[names(sim8) == "item_id"] <- "item"
  # Ensure initial dataset has both 0 and 1 in each cell, then drop row 3 via NA
  sim8$y[sim8$study_method == "restudy"][1:2] <- c(0L, 1L)
  sim8$y[sim8$study_method == "retrieval"][1:2] <- c(0L, 1L)
  sim8$y[3L] <- NA_integer_
  sim8$study_method <- factor(sim8$study_method, levels = design$factors$study_method)
  stats::contrasts(sim8$study_method) <- design$contrasts$study_method

  env8 <- new.env(parent = globalenv())
  env8$design <- design
  env8$target <- specify_target(
    dgp_ext, "study_method", scale = "pp", engine = "lme4",
    control = target_control(prob_boot = 40L)
  )
  env8$responses <- sim8
  # Evaluate up through marginal_diff_pp definition, then inspect bootMer draws directly
  fit8 <- suppressMessages(lme4::glmer(stats::formula(env8$target), data = sim8, family = stats::family(env8$target)))
  expect_equal(nrow(stats::model.frame(fit8)), 31L)
  env8$valid_glmer_fit <- env54$valid_glmer_fit
  environment(env8$valid_glmer_fit) <- env8

  # Force separated response vector on the 31 retained rows (all 0 in "restudy")
  mf8 <- stats::model.frame(fit8)
  y_sep <- ifelse(as.character(mf8$study_method) == "restudy", 0L, rep(c(0L, 1L), length.out = nrow(mf8)))
  fit8_sep <- suppressWarnings(suppressMessages(lme4::refit(fit8, y_sep)))
  expect_false(env8$valid_glmer_fit(fit8_sep))

  # Across 40 parametric draws on 8x4 with 1 dropped NA row, every truly separated draw is caught
  set.seed(312L)
  sim_draws <- stats::simulate(fit8, nsim = 40L, use.u = FALSE)
  n_true_sep <- 0L
  n_missed_sep <- 0L
  for (j in seq_len(ncol(sim_draws))) {
    yj <- as.integer(sim_draws[[j]])
    is_sep <- !all(tapply(yj, as.character(mf8$study_method), function(v) any(v == 0L) && any(v == 1L)))
    if (is_sep) {
      n_true_sep <- n_true_sep + 1L
      mj <- suppressWarnings(suppressMessages(tryCatch(lme4::refit(fit8, yj), error = function(e) NULL)))
      if (!is.null(mj) && isTRUE(env8$valid_glmer_fit(mj))) {
        n_missed_sep <- n_missed_sep + 1L
      }
    }
  }
  expect_gt(n_true_sep, 0L)
  expect_equal(n_missed_sep, 0L)
})

test_that("surrogate search handles min_valid_rate >= 0.99 and empty predicted feasibility without skipping 8x8", {
  des <- design_crossed(list(condition = c("a", "b")))
  dgp <- calibrate_dgp(des, c(a = 0.35, b = 0.70), subj_sd = 0.4, item_sd = 0.3)

  real_eval <- evaluate_grid_batch
  local_mocked_bindings(
    evaluate_grid_batch = function(spec, grid_sub, point_bases, nsim, ...) {
      base <- real_eval(spec, grid_sub[1L, , drop = FALSE], point_bases[1L], nsim = 2L,
                        engine_use = spec$engine, fast_stage = TRUE, workers = 1L,
                        parallel_backend = "serial", verbose = FALSE)$summary
      rows <- lapply(seq_len(nrow(grid_sub)), function(i) {
        r <- base
        ns <- grid_sub$n_subj[i]
        ni <- grid_sub$n_item[i]
        r$n_subj <- ns
        r$n_item <- ni
        r$nsim <- as.integer(nsim)
        r$n_valid <- as.integer(nsim)
        r$valid_rate <- 1.0
        r$nonconv_rate <- 0.0
        r$conditional_power <- 1.0
        r$effective_power <- 1.0
        r$effective_rope_prob <- 1.0
        r$power <- 1.0
        mc_p <- wilson_mc_ci(nsim, nsim, conf_level = spec$mc_conf_level)
        r$power_mc_lo <- unname(mc_p["lo"])
        r$power_mc_hi <- unname(mc_p["hi"])
        r$criterion_met <- TRUE
        r$decision_status <- compute_decision_status(
          spec, TRUE, mc_p, c(lo = NA_real_, hi = NA_real_),
          c(lo = NA_real_, hi = NA_real_), c(lo = NA_real_, hi = NA_real_)
        )
        r$mc_conf_level_used <- spec$mc_conf_level
        r
      })
      list(summary = do.call(rbind, rows), raw = vector("list", nrow(grid_sub)))
    },
    .package = "glmmplanr"
  )

  for (mvr in c(0.99, 1.00)) {
    sp_hi <- specify_target(
      dgp, "condition", scale = "pp", criterion = "power", power_min = 0.80,
      prob_se_method = "joint_delta",
      control = target_control(failure_policy = "conditional", min_valid_rate = mvr)
    )
    plan_hi <- optimize_sample_size(
      sp_hi,
      n_subj_range = c(8L, 14L),
      n_item_range = c(8L, 12L),
      prefer = "fewer_trials",
      method = "surrogate_gp",
      nsim_final = 40L,
      workers = 1L, verbose = FALSE,
      control = optimize_control(
        n_eval_init = 6L, n_eval_adapt = 2L, n_verify = 4L,
        nsim_init = 40L, nsim_sensitivity = 0L
      )
    )
    expect_true(plan_hi$feasible)
    expect_equal(plan_hi$optimal$n_subj, 8L)
    expect_equal(plan_hi$optimal$n_item, 8L)
    expect_identical(plan_hi$optimal$decision_status, "pass")
  }

  # Also test select_stage2_candidates directly when pred_criterion_met is all FALSE
  # while an observed candidate (8, 8) met its criterion in Stage 1
  surf_zero <- plan_hi$surface
  surf_zero$pred_criterion_met <- rep(FALSE, nrow(surf_zero))
  surf_zero$obs_criterion_met <- surf_zero$n_subj == 8L & surf_zero$n_item == 8L
  cand_idx <- select_stage2_candidates(surf_zero, sp_hi, prefer = "fewer_trials", n_verify = 2L)
  idx_8x8 <- which(surf_zero$n_subj == 8L & surf_zero$n_item == 8L)
  expect_equal(cand_idx[1L], idx_8x8)
})

test_that("unattainable-zero weighted contrasts display as weighted contrasts in pp rather than negative percentages", {
  des <- design_crossed(list(condition = c("a", "b")))
  dgp <- calibrate_dgp(des, c(a = 0.53, b = 0.80), subj_sd = 0.5, item_sd = 0.4, chance = 1 / 3)
  sp_w <- specify_target(
    dgp, "condition",
    scale = "pp", criterion = "ci_width", ci_width = 35, prob_se_method = "joint_delta",
    control = target_control(contrast_weights = c(a = -1, b = 0.2))
  )
  expect_equal(sp_w$true_prob_diff_pp, -37.0, tolerance = 1e-6)
  expect_identical(prob_target_kind(sp_w), "weighted_contrast")

  sp_txt <- gsub("\\s+", " ", paste(cli::ansi_strip(cli::cli_fmt(print(sp_w))), collapse = " "))
  expect_match(sp_txt, "weighted contrast = -37.00 pp", fixed = TRUE)
  expect_false(grepl("absolute probability", sp_txt, fixed = TRUE))
  expect_false(grepl("-37.00%", sp_txt, fixed = TRUE))

  plan_w <- optimize_sample_size(
    sp_w, n_subj_range = c(12L, 16L), n_item_range = 8L,
    method = "grid", nsim_final = 20L, workers = 1L, verbose = FALSE,
    control = optimize_control(nsim_sensitivity = 0L)
  )
  pl_txt <- gsub("\\s+", " ", paste(cli::ansi_strip(cli::cli_fmt(print(plan_w))), collapse = " "))
  expect_match(pl_txt, "Weighted contrast (marginal): true = -37.0 pp", fixed = TRUE)
  expect_false(grepl("Target accuracy", pl_txt, fixed = TRUE))
  expect_false(grepl("-37.0%", pl_txt, fixed = TRUE))

  rpt_w <- report_plan(plan_w, print_output = FALSE)
  expect_match(rpt_w, "marginal weighted accuracy contrast of -37.00 percentage points", fixed = TRUE)
  expect_match(rpt_w, "zero lies outside the attainable range of this weighted probability contrast", fixed = TRUE)
  expect_false(grepl("accuracy of -37.00%", rpt_w, fixed = TRUE))
  expect_false(grepl("absolute probability target", rpt_w, fixed = TRUE))

  if (requireNamespace("gt", quietly = TRUE)) {
    old_opt <- options(sass.cache = FALSE)
    on.exit(options(old_opt), add = TRUE)
    gt_html <- gt::as_raw_html(as_gt(plan_w))
    expect_match(gt_html, "weighted contrast = -37.0 pp", fixed = TRUE)
    expect_false(grepl("absolute probability", gt_html, fixed = TRUE))
  }
})

test_that("package build-time source MD5 stamp is embedded in the installed namespace", {
  expect_true(is.character(.glmmplanr_source_md5))
  expect_gt(length(.glmmplanr_source_md5), 10L)
  expect_true("R/design.R" %in% names(.glmmplanr_source_md5))
  expect_true("src/glmm_laplace.cpp" %in% names(.glmmplanr_source_md5))
  expect_true(all(nzchar(.glmmplanr_source_md5)))
})

test_that("single-cell and tiny contrast weights preserve display semantics and exact scale invariance", {
  des1 <- design_crossed(list())
  dgp1 <- calibrate_dgp(des1, 0.60, subj_sd = 0.4, item_sd = 0.3)

  # Reject all-zero or non-finite weights at entry points
  expect_error(target_control(contrast_weights = 0), "non-zero")
  expect_error(target_control(contrast_weights = c(0, 0)), "non-zero")
  expect_error(target_control(contrast_weights = Inf), "finite")
  expect_error(cell_contrast_weights(des1, "(Intercept)", custom_weights = 0), "non-zero")
  expect_error(cell_contrast_weights(des1, "(Intercept)", custom_weights = NaN), "finite")

  # Default weight 1 is an absolute probability; any weight != 1 is a weighted contrast in pp
  sp1 <- specify_target(dgp1, "(Intercept)", scale = "pp", criterion = "ci_width",
                        ci_width = 25, prob_se_method = "joint_delta")
  expect_identical(prob_target_kind(sp1), "absolute")
  txt1 <- gsub("\\s+", " ", paste(cli::ansi_strip(cli::cli_fmt(print(sp1))), collapse = " "))
  expect_match(txt1, "absolute probability = 60.00%", fixed = TRUE)

  for (w_val in c(0.5, 1.3, 2.0, -2.0, 1e-7)) {
    sp_k <- specify_target(
      dgp1, "(Intercept)", scale = "pp", criterion = "ci_width",
      ci_width = 25 * abs(w_val), prob_se_method = "joint_delta",
      control = target_control(contrast_weights = w_val)
    )
    expect_identical(prob_target_kind(sp_k), "weighted_contrast")
    txt_k <- gsub("\\s+", " ", paste(cli::ansi_strip(cli::cli_fmt(print(sp_k))), collapse = " "))
    expect_match(txt_k, "weighted contrast =", fixed = TRUE)
    expect_false(grepl("absolute probability", txt_k, fixed = TRUE))
  }

  # Scale invariance across w = 1, 0.5, 2, -2, 1e-7, -1e-7 under joint_delta and bootstrap
  weights_test <- c(1.0, 0.5, 2.0, -2.0, 1e-7, -1e-7)
  for (se_m in c("joint_delta", "bootstrap")) {
    nsim_use <- if (se_m == "joint_delta") 16L else 8L
    res_list <- lapply(weights_test, function(w_val) {
      sp_k <- specify_target(
        dgp1, "(Intercept)", scale = "pp", criterion = "ci_width",
        ci_width = 25 * abs(w_val), prob_se_method = se_m,
        control = target_control(contrast_weights = w_val, prob_boot = 19L)
      )
      set.seed(7711L)
      ev_k <- evaluate_design_point(sp_k, n_subj = 20L, n_item = 12L, nsim = nsim_use)
      list(
        w = w_val,
        n_valid = ev_k$n_valid,
        est_unit = ev_k$prob_diff_est_pp / w_val,
        width_unit = ev_k$prob_ci_width_mean_pp / abs(w_val),
        aipe = ev_k$aipe_prob,
        cov = ev_k$prob_coverage,
        decision = ev_k$decision_status
      )
    })
    ref <- res_list[[1L]]
    for (k in 2:length(res_list)) {
      cur <- res_list[[k]]
      expect_equal(cur$n_valid, ref$n_valid)
      expect_equal(cur$est_unit, ref$est_unit, tolerance = 1e-11)
      expect_equal(cur$width_unit, ref$width_unit, tolerance = 1e-11)
      expect_equal(cur$aipe, ref$aipe, tolerance = 1e-12)
      expect_equal(cur$cov, ref$cov, tolerance = 1e-12)
      expect_identical(cur$decision, ref$decision)
    }
  }
})

test_that("n_eval_init = 0 is rejected only when surrogate_gp explores a strict subset, and n_eval_init = 1 succeeds", {
  des <- design_crossed(list(condition = c("a", "b")))
  dgp <- calibrate_dgp(des, c(a = 0.48, b = 0.52), subj_sd = 0.4, item_sd = 0.3)
  sp <- specify_target(dgp, "condition", scale = "pp", prob_se_method = "joint_delta")

  # 1. Explicit method = "grid" with n_eval_init = 0, n_eval_adapt = 0 succeeds
  p_grid <- optimize_sample_size(
    sp, n_subj_range = c(8L, 12L), n_item_range = c(8L, 12L),
    method = "grid", nsim_final = 4L, workers = 1L, verbose = FALSE,
    control = optimize_control(n_eval_init = 0L, n_eval_adapt = 0L, nsim_sensitivity = 0L)
  )
  expect_equal(nrow(p_grid$verified), 9L)

  # 2. surrogate_gp with n_eval_init = 0 succeeds when post-constraint grid <= n_eval_adapt (4 <= 4)
  p_auto <- optimize_sample_size(
    sp, n_subj_range = c(8L, 10L), n_item_range = c(8L, 10L),
    method = "surrogate_gp", nsim_final = 4L, workers = 1L, verbose = FALSE,
    control = optimize_control(n_eval_init = 0L, n_eval_adapt = 4L, nsim_sensitivity = 0L)
  )
  expect_equal(nrow(p_auto$verified), 4L)

  # 3. surrogate_gp with n_eval_init = 0 aborts when post-constraint grid > n_eval_adapt (4 > 2)
  expect_error(
    optimize_sample_size(
      sp, n_subj_range = c(8L, 10L), n_item_range = c(8L, 10L),
      method = "surrogate_gp", nsim_final = 4L, workers = 1L, verbose = FALSE,
      control = optimize_control(n_eval_init = 0L, n_eval_adapt = 2L, nsim_sensitivity = 0L)
    ),
    "n_eval_init"
  )

  # Multi-scenario also aborts cleanly when n_grid > n_eval_init + n_eval_adapt and n_eval_init == 0
  dgp2 <- calibrate_dgp(des, c(a = 0.45, b = 0.55), subj_sd = 0.4, item_sd = 0.3)
  expect_error(
    optimize_sample_size(
      sp, n_subj_range = c(8L, 10L), n_item_range = c(8L, 10L),
      scenarios = list(s1 = dgp, s2 = dgp2),
      method = "surrogate_gp", nsim_final = 4L, workers = 1L, verbose = FALSE,
      control = optimize_control(n_eval_init = 0L, n_eval_adapt = 2L, nsim_sensitivity = 0L)
    ),
    "n_eval_init"
  )

  # 4. Exact reproducer for n_eval_init = 1, n_eval_adapt = 2 on 8..10 x 8..10 (no non-finite optim failure)
  for (nv in c(1L, 2L)) {
    p_1_2 <- optimize_sample_size(
      sp, n_subj_range = c(8L, 10L), n_item_range = c(8L, 10L),
      method = "surrogate_gp", nsim_final = 8L, seed = 1L,
      workers = 1L, verbose = FALSE,
      control = optimize_control(
        n_eval_init = 1L, n_eval_adapt = 2L, n_verify = nv,
        nsim_init = 8L, nsim_sensitivity = 0L
      )
    )
    expect_true(all(is.finite(p_1_2$surface$pred_power)))
    expect_true(all(is.finite(p_1_2$surface$pred_ci_width)))
  }

  # Direct rank-deficient 3-point constant-response check on fit_nonneg_slope_wls
  H3 <- cbind(1, -1 / c(8, 8, 10), -1 / c(8, 10, 8), -1 / c(64, 80, 80))
  b_const <- fit_nonneg_slope_wls(H3, rep(-1.405, 3L), rep(12.0, 3L))
  expect_true(all(is.finite(b_const)))
  expect_equal(b_const, c(-1.405, 0, 0, 0), tolerance = 1e-12)
})

test_that("validation plots respect small nonzero custom true_target and record calibrated nulls", {
  des <- design_crossed(list(condition = c("a", "b")))
  dgp_plan <- calibrate_dgp(des, c(a = 0.45, b = 0.65), subj_sd = 0.3, item_sd = 0.3)

  # Custom generator produces +20 pp on BOTH null and planned rows; scale weights & truth by 1e-8
  gen_20pp <- function(schedule, spec) {
    simulate_dgp_data(dgp_plan, schedule = schedule)
  }
  k_tiny <- 1e-8
  sp_prob <- specify_target(
    dgp_plan, "condition", scale = "prob", criterion = "power",
    prob_se_method = "joint_delta",
    control = target_control(contrast_weights = k_tiny * c(a = -1, b = 1))
  )
  val_nonnull <- validate_design_point(
    sp_prob, n_subj = 24L, n_item = 16L, nsim = 8L, seed = 401L,
    data_fn = gen_20pp,
    true_target = c(null = 0.20 * k_tiny, planned = 0.20 * k_tiny)
  )
  expect_false(isTRUE(attr(val_nonnull, "validation_settings")$null_verified))
  p_nonnull <- plot(val_nonnull)
  expect_warning(ggplot2::ggplot_build(p_nonnull), NA)
  panels_nn <- as.character(p_nonnull$data$panel)
  expect_true("Null-input detection" %in% panels_nn)
  expect_false("Type I error" %in% panels_nn)

  # When custom data_fn explicitly declares true_target["null"] == 0, Type I error IS shown
  val_zero <- validate_design_point(
    sp_prob, n_subj = 24L, n_item = 16L, nsim = 4L, seed = 402L,
    data_fn = gen_20pp,
    true_target = c(null = 0, planned = 0.20 * k_tiny)
  )
  expect_true(isTRUE(attr(val_zero, "validation_settings")$null_verified))
  p_zero <- plot(val_zero)
  expect_true("Type I error" %in% as.character(p_zero$data$panel))
})

test_that("named contrast weights align consistently across public helpers and reject partial names", {
  des_ab <- design_crossed(list(condition = c("a", "b")))
  dgp_ab <- calibrate_dgp(des_ab, c(a = 0.45, b = 0.65), subj_sd = 0.4, item_sd = 0.3)

  w_rev <- c(b = 1, a = -1)
  sp_rev <- specify_target(
    dgp_ab, "condition", scale = "pp", prob_se_method = "joint_delta",
    control = target_control(contrast_weights = w_rev)
  )
  expect_equal(sp_rev$true_prob_diff_pp, 20.0, tolerance = 1e-6)
  expect_equal(100 * true_prob_effect(dgp_ab, w_rev), 20.0, tolerance = 1e-6)
  expect_equal(cell_contrast_weights(des_ab, "condition", custom_weights = w_rev), c(a = -1, b = 1))

  # Passing named weights extracted from a reversed-order design preserves exact cell direction
  des_ba <- design_crossed(list(condition = c("b", "a")))
  dgp_ba <- calibrate_dgp(des_ba, c(b = 0.65, a = 0.45), subj_sd = 0.4, item_sd = 0.3)
  w_from_ba <- cell_contrast_weights(des_ba, "condition")
  sp_from_ba <- specify_target(
    dgp_ab, "condition", scale = "pp", prob_se_method = "joint_delta",
    control = target_control(contrast_weights = w_from_ba)
  )
  expect_equal(sp_from_ba$true_prob_diff_pp, 100 * true_prob_effect(dgp_ab, w_from_ba), tolerance = 1e-8)

  # Multi-scenario optimization with reversed cell ordering aligns contrast weights per scenario
  plan_ms_ord <- optimize_sample_size(
    sp_rev, n_subj_range = c(12L, 16L), n_item_range = 8L,
    scenarios = list(ab = dgp_ab, ba = dgp_ba),
    method = "grid", nsim_final = 4L, workers = 1L, verbose = FALSE,
    control = optimize_control(nsim_sensitivity = 0L)
  )
  expect_equal(plan_ms_ord$robust$specs$ab$true_prob_diff_pp, 20.0, tolerance = 1e-6)
  expect_equal(plan_ms_ord$robust$specs$ba$true_prob_diff_pp, 20.0, tolerance = 1e-6)

  # Partial or duplicate names are rejected consistently
  expect_error(target_control(contrast_weights = c(b = 1, -1)), "name every element")
  expect_error(cell_contrast_weights(des_ab, "condition", custom_weights = c(b = 1, -1)), "name every element")
  expect_error(true_prob_effect(dgp_ab, c(b = 1, -1)), "name every element")
  expect_error(true_prob_effect(dgp_ab, c(a = 1, a = -1)), "without duplicates")
})

test_that("null calibration and Type M / Type S diagnostics are invariant across weight scales", {
  des <- design_crossed(list(condition = c("a", "b")))
  dgp <- calibrate_dgp(des, c(a = 0.35, b = 0.70), subj_sd = 0.35, item_sd = 0.3)

  for (k_mult in c(1.0, 1e-12, 1e5)) {
    sp_pow <- specify_target(
      dgp, "condition", scale = "pp", criterion = "power",
      prob_se_method = "joint_delta",
      control = target_control(contrast_weights = k_mult * c(a = -1, b = 1))
    )
    expect_true(sp_pow$has_null_target)

    sp_ci <- specify_target(
      dgp, "condition", scale = "pp", criterion = "ci_width",
      ci_width = 30 * k_mult, prob_se_method = "joint_delta",
      control = target_control(contrast_weights = k_mult * c(a = -1, b = 1))
    )
    val_k <- validate_design_point(sp_ci, n_subj = 16L, n_item = 10L, nsim = 4L, seed = 611L)
    expect_equal(nrow(val_k), 2L)
    expect_identical(val_k$condition, c("null", "planned"))
    expect_true(isTRUE(attr(val_k, "validation_settings")$null_verified))
    p_k <- plot(val_k)
    expect_true("Type I error" %in% as.character(p_k$data$panel))
  }

  # At k = 1e-9, Type M and Type S diagnostics are finite and identical between prob and pp
  sp_tiny_prob <- specify_target(
    dgp, "condition", scale = "prob", criterion = "power",
    prob_se_method = "joint_delta",
    control = target_control(contrast_weights = 1e-9 * c(a = -1, b = 1))
  )
  sp_tiny_pp <- specify_target(
    dgp, "condition", scale = "pp", criterion = "power",
    prob_se_method = "joint_delta",
    control = target_control(contrast_weights = 1e-9 * c(a = -1, b = 1))
  )
  set.seed(712L)
  ev_prob <- evaluate_design_point(sp_tiny_prob, n_subj = 24L, n_item = 16L, nsim = 12L)
  set.seed(712L)
  ev_pp <- evaluate_design_point(sp_tiny_pp, n_subj = 24L, n_item = 16L, nsim = 12L)
  expect_gt(ev_prob$power, 0)
  expect_true(is.finite(ev_prob$type_m_error))
  expect_true(is.finite(ev_prob$type_s_error))
  expect_equal(ev_prob$type_m_error, ev_pp$type_m_error, tolerance = 1e-12)
  expect_equal(ev_prob$type_s_error, ev_pp$type_s_error, tolerance = 1e-12)
})

test_that("pooled top-up rejects pre-v4 plans affected by unnormalized weight SEs or partial fast-path", {
  des1 <- design_crossed(list())
  dgp1 <- calibrate_dgp(des1, 0.60, subj_sd = 0.4, item_sd = 0.3)
  sp_tiny <- specify_target(
    dgp1, "(Intercept)", scale = "pp", criterion = "ci_width",
    ci_width = 30 * 1e-7, prob_se_method = "joint_delta",
    control = target_control(contrast_weights = 1e-7, ci_assurance = 0.50)
  )
  plan_v4 <- optimize_sample_size(
    sp_tiny, n_subj_range = 16L, n_item_range = 8L,
    method = "grid", nsim_final = 3L, seed = 919L,
    workers = 1L, verbose = FALSE, control = optimize_control(nsim_sensitivity = 0L)
  )
  expect_equal(plan_v4$search_meta$simulation_version, 5L)
  expect_identical(plan_v4$verified$decision_status, "indeterminate")

  plan_legacy <- plan_v4
  plan_legacy$search_meta$simulation_version <- 3L
  expect_error(
    topup_verification(plan_legacy, nsim_extra = 9L, topup_policy = "pooled", verbose = FALSE),
    "unnormalized contrast weights"
  )

  # Both pooled and independent top-up reject incompatible pre-v4 plans so stale
  # pass/fail candidate rows from unnormalized SEs can never mix with new rows
  expect_error(
    topup_verification(plan_legacy, nsim_extra = 9L, topup_policy = "independent", verbose = FALSE),
    "unnormalized contrast weights"
  )

  # Multi-candidate pre-v4 plan with both a stale fail row and an indeterminate row is also rejected
  des2 <- design_crossed(list(condition = c("a", "b")))
  dgp2 <- calibrate_dgp(des2, c(a = 0.45, b = 0.60), subj_sd = 0.35, item_sd = 0.30)
  sp_w <- specify_target(
    dgp2, "condition", scale = "pp", criterion = "power",
    prob_se_method = "joint_delta",
    control = target_control(contrast_weights = 1e-7 * c(a = -1, b = 1))
  )
  plan_multi <- optimize_sample_size(
    sp_w, n_subj_range = c(12L, 24L), n_item_range = 16L,
    method = "grid", nsim_final = 16L, seed = 404L,
    workers = 1L, verbose = FALSE, control = optimize_control(nsim_sensitivity = 0L)
  )
  legacy_multi <- plan_multi
  legacy_multi$search_meta$simulation_version <- 3L
  legacy_multi$verified$decision_status[legacy_multi$verified$n_subj == 24L] <- "fail"
  legacy_multi$verified$decision_status[legacy_multi$verified$n_subj == 12L] <- "indeterminate"
  expect_error(
    topup_verification(legacy_multi, nsim_extra = 24L, topup_policy = "independent", verbose = FALSE),
    "unnormalized contrast weights"
  )

  # Current v4 plan pools cleanly (with standard adaptive-pooling warning)
  expect_warning(
    plan_pooled <- topup_verification(plan_v4, nsim_extra = 9L, topup_policy = "pooled", verbose = FALSE),
    "adaptive"
  )
  expect_equal(plan_pooled$verified$nsim, 12L)
  expect_identical(plan_pooled$verified$decision_status, "pass")
})

test_that("factor names cannot collide with reserved schedule or formula identifiers", {
  reserved <- c(
    "subj", "item", "testlet",
    "subj_id", "item_id", "testlet_id",
    "list_id", "cell_id", "cell_name", "trial_id", "y",
    "Z_subj_cells", "Z_item_cells", "Z_testlet_cells"
  )
  for (bad_nm in reserved) {
    fac <- stats::setNames(list(c("a", "b")), bad_nm)
    expect_error(design_crossed(fac), "reserved")
    expect_error(design_testlet(fac, testlet_size = 2L), "reserved")
    fac_ws <- stats::setNames(list(c("a", "b")), paste0(" ", bad_nm, " "))
    expect_error(design_crossed(fac_ws), "reserved")
  }
})

test_that("random-item Rasch rejects variance-unidentifiable one-to-one and one-item-per-person schedules", {
  # Find an item_seed where 8x8 with items_per_subject = 1 presents each item once (1-to-1 bijection)
  seed_bij <- NULL
  for (s_try in 1:5000) {
    d_try <- design_crossed(list(), allocation = "crossed", items_per_subject = 1L, item_seed = s_try)
    sc_try <- generate_design_schedule(d_try, 8L, 8L)
    if (all(tabulate(sc_try$trials$item_id, 8L) == 1L)) {
      seed_bij <- s_try
      break
    }
  }
  expect_false(is.null(seed_bij))

  des_bij <- design_crossed(list(), allocation = "crossed", items_per_subject = 1L, item_seed = seed_bij)
  dgp_bij <- calibrate_dgp(des_bij, 0.50, subj_sd = 0.27, item_sd = 0.26)
  for (eng in c("cpp_laplace", "lme4")) {
    res_bij <- plan_rasch_precision(
      dgp_bij,
      data.frame(n_subj = 8L, n_item = 8L),
      max_error = c(person_sd = 0.01, item_sd = 0.01),
      assurance = 0.20,
      nsim = 16L,
      seed = 101L,
      engine = eng
    )
    expect_true(all(res_bij$summary$n_valid_fit == 0L))
    expect_true(all(res_bij$summary$n_evaluable == 0L))
    expect_true(all(res_bij$summary$n_met == 0L))
    expect_true(all(res_bij$summary$decision == "fail"))
    expect_true(all(grepl("unidentifiable|Fewer than two", res_bij$replications$fit_error)))
  }

  # 1 item per person on a 4-item bank: person variance is unidentifiable under random_item too
  des_1 <- design_crossed(list(), allocation = "crossed", items_per_subject = 1L, item_seed = 1L)
  dgp_1 <- calibrate_dgp(des_1, 0.55, subj_sd = 0.80, item_sd = 0.60)
  res_p1 <- plan_rasch_precision(
    dgp_1, data.frame(n_subj = 40L, n_item = 4L),
    max_error = c(person_sd = 0.50, item_sd = 0.50),
    assurance = 0.20, nsim = 16L, seed = 202L
  )
  expect_true(all(res_p1$summary$n_valid_fit == 0L))
  expect_true(all(res_p1$summary$decision == "fail"))

  fb_rnd_1 <- plan_rasch_fixed_bank(
    item_difficulties = c(-1.0, -0.3, 0.3, 1.0),
    n_subj = 40L,
    max_error = c(person_sd = 0.50, difficulty = 0.50),
    items_per_subject = 1L,
    method = "random_item",
    assurance = 0.20,
    nsim = 16L,
    seed = 203L
  )
  expect_true(all(fb_rnd_1$summary$n_valid_fit == 0L))
  expect_true(all(fb_rnd_1$summary$decision == "fail"))
})

test_that("simulate_dgp_data rejects schedules from incompatible or level-reversed designs", {
  des_ab <- design_crossed(list(condition = c("a", "b")))
  des_ba <- design_crossed(list(condition = c("b", "a")))
  dgp_ab <- calibrate_dgp(des_ab, c(a = 0.20, b = 0.80), subj_sd = 0.2, item_sd = 0.2)

  sched_ba <- generate_design_schedule(des_ba, 64L, 24L)
  expect_error(
    simulate_dgp_data(dgp_ab, schedule = sched_ba),
    "match|levels|Cell ordering"
  )

  # Even if schedule$cells metadata is stripped (simulating a legacy schedule), factor levels catch the reversal
  sched_ba_legacy <- sched_ba
  sched_ba_legacy$cells <- NULL
  expect_error(
    simulate_dgp_data(dgp_ab, schedule = sched_ba_legacy),
    "levels|Cell-to-condition"
  )

  # Equivalent design with identical factor & level ordering succeeds and preserves true direction
  des_ab2 <- design_crossed(list(condition = c("a", "b")))
  sched_ab2 <- generate_design_schedule(des_ab2, 64L, 24L)
  set.seed(501L)
  df_ok <- simulate_dgp_data(dgp_ab, schedule = sched_ab2)
  means <- tapply(df_ok$y, df_ok$condition, mean)
  expect_lt(means[["a"]], 0.40)
  expect_gt(means[["b"]], 0.60)
})

test_that("cell_probs and contrast_weights share strict partial-name rejection and exact cell-name matching", {
  des_ab <- design_crossed(list(condition = c("a", "b")))
  expect_error(
    calibrate_dgp(des_ab, c(b = 0.65, 0.45), subj_sd = 0.3, item_sd = 0.3),
    "name every element"
  )
  expect_error(
    calibrate_dgp(des_ab, c(a = 0.45, a = 0.65), subj_sd = 0.3, item_sd = 0.3),
    "without duplicates"
  )
  dgp_rev <- calibrate_dgp(des_ab, c(b = 0.65, a = 0.45), subj_sd = 0.3, item_sd = 0.3)
  expect_equal(unname(dgp_rev$target_cell_probs[c("a", "b")]), c(0.45, 0.65), tolerance = 1e-12)
  sp_rev <- specify_target(dgp_rev, "condition", scale = "pp", prob_se_method = "joint_delta")
  expect_equal(sp_rev$true_prob_diff_pp, 20.0, tolerance = 1e-6)

  # Factor levels containing spaces or dots work end-to-end on both logit and pp scales
  for (lvls in list(c("pre test", "post test"), c("pre.test", "post.test"))) {
    des_sp <- design_crossed(list(time = lvls))
    dgp_unnamed <- calibrate_dgp(des_sp, c(0.40, 0.65), subj_sd = 0.3, item_sd = 0.3)
    named_p <- stats::setNames(c(0.65, 0.40), rev(lvls))
    dgp_named <- calibrate_dgp(des_sp, named_p, subj_sd = 0.3, item_sd = 0.3)
    expect_equal(dgp_unnamed$beta, dgp_named$beta, tolerance = 1e-10)

    for (sc in c("logit", "pp")) {
      sp <- specify_target(dgp_named, "time", scale = sc, prob_se_method = "joint_delta")
      ev <- evaluate_design_point(sp, n_subj = 12L, n_item = 8L, nsim = 3L, stream_seeds = 77L)
      expect_equal(ev$n_valid, 3L)
    }
  }
})

test_that("multi-scenario sensitivity updates when limiting scenario changes at the same (N_subj, N_item)", {
  des <- design_crossed(list(condition = c("a", "b")))
  dgp_s10 <- calibrate_dgp(des, c(a = 0.50, b = 0.60), subj_sd = c(0.5, 0.25), item_sd = c(0.4, 0.20))
  dgp_s03 <- calibrate_dgp(des, c(a = 0.50, b = 0.53), subj_sd = c(0.5, 0.25), item_sd = c(0.4, 0.20))
  sp_base <- specify_target(dgp_s10, "condition", scale = "pp", prob_se_method = "joint_delta")

  plan <- optimize_sample_size(
    sp_base,
    n_subj_range = 16L,
    n_item_range = 8L,
    scenarios = list(s10 = dgp_s10, s03 = dgp_s03),
    method = "grid",
    nsim_final = 4L,
    seed = 808L,
    workers = 1L,
    verbose = FALSE,
    control = optimize_control(nsim_sensitivity = 3L)
  )
  expect_false(is.null(plan$sensitivity))
  expect_identical(attr(plan$sensitivity, "source_scenario"), as.character(plan$optimal$scenario[1L]))

  # Simulate scenario flip at the same (16, 8) candidate before top-up vs switch_preference
  flipped <- plan
  ev <- flipped$robust$scenario_evaluations
  # Make s03 pass and s10 indeterminate so s10 becomes the limiting scenario at (16, 8)
  i10 <- which(ev$scenario == "s10")
  i03 <- which(ev$scenario == "s03")
  ev$decision_status[i03] <- "pass"
  ev$power[i03] <- 0.95
  ev$power_mc_lo[i03] <- 0.85
  ev$power_mc_hi[i03] <- 0.99
  ev$criterion_met[i03] <- TRUE
  ev$decision_status[i10] <- "indeterminate"
  ev$power[i10] <- 0.75
  ev$power_mc_lo[i10] <- 0.40
  ev$power_mc_hi[i10] <- 0.92
  ev$criterion_met[i10] <- FALSE
  flipped$robust$scenario_evaluations <- ev
  # Force stored sensitivity to belong to s03 before re-ranking / top-up
  attr(flipped$sensitivity, "source_scenario") <- "s03"
  flipped$sensitivity$source_scenario <- rep("s03", nrow(flipped$sensitivity))
  flipped$sensitivity$prob_diff_true_pp <- rep(3.0, nrow(flipped$sensitivity))

  # switch_preference clears stale sensitivity when limiting scenario changed from s03 to s10
  reranked <- switch_preference(flipped, "balanced", verbose = FALSE)
  expect_identical(as.character(reranked$optimal$scenario[1L]), "s10")
  expect_null(reranked$sensitivity)

  # topup_verification recomputes sensitivity for the newly limiting scenario at (16, 8)
  topped <- topup_verification(flipped, nsim_extra = 4L, topup_policy = "independent", verbose = FALSE)
  expect_false(is.null(topped$sensitivity))
  expect_identical(attr(topped$sensitivity, "source_scenario"), as.character(topped$optimal$scenario[1L]))
  expect_equal(
    topped$sensitivity$prob_diff_true_pp[1L],
    topped$robust$specs[[as.character(topped$optimal$scenario[1L])]]$true_prob_diff_pp,
    tolerance = 1e-6
  )
})

test_that("as_gt header reflects evaluated scenario and report_plan sensitivity range includes standard row", {
  des <- design_crossed(list(condition = c("a", "b")))
  dgp_unused60 <- calibrate_dgp(des, c(a = 0.20, b = 0.80), subj_sd = c(0.5, 0.25), item_sd = c(0.4, 0.20))
  dgp_10 <- calibrate_dgp(des, c(a = 0.50, b = 0.60), subj_sd = c(0.5, 0.25), item_sd = c(0.4, 0.20))
  dgp_03 <- calibrate_dgp(des, c(a = 0.50, b = 0.53), subj_sd = c(0.5, 0.25), item_sd = c(0.4, 0.20))
  sp_unused <- specify_target(dgp_unused60, "condition", scale = "pp", prob_se_method = "joint_delta")

  plan <- optimize_sample_size(
    sp_unused,
    n_subj_range = 16L,
    n_item_range = 8L,
    scenarios = list(ten = dgp_10, three = dgp_03),
    method = "grid",
    nsim_final = 3L,
    seed = 909L,
    workers = 1L,
    verbose = FALSE,
    control = optimize_control(nsim_sensitivity = 3L)
  )
  gt_tbl <- as_gt(plan)
  sub_txt <- as.character(gt_tbl$`_heading`$subtitle)
  expect_false(grepl("60\\.0 pp", sub_txt))
  lim_sc <- as.character(plan$optimal$scenario[1L])
  lim_pp <- sprintf("%+.1f pp", plan$robust$specs[[lim_sc]]$true_prob_diff_pp)
  expect_match(sub_txt, lim_pp, fixed = TRUE)
  expect_match(sub_txt, paste0("scenario = ", lim_sc), fixed = TRUE)

  # Force a feasible plan whose standard sensitivity row has the minimum CI width (39.19 < 39.71 < 40.20)
  plan_feas <- plan
  plan_feas$feasible <- TRUE
  plan_feas$optimal$decision_status <- "pass"
  plan_feas$sensitivity$prob_ci_width_mean_pp <- c(39.71, 39.19, 40.20)
  plan_feas$sensitivity$effective_power <- c(0.84, 0.88, 0.79)
  rep_txt <- report_plan(plan_feas, print_output = FALSE)
  expect_match(rep_txt, "mean CI width: 39.19 to 40.20 percentage points", fixed = TRUE)
  expect_match(rep_txt, "effective power ranged from 79.0% to 88.0%", fixed = TRUE)
})

test_that("validate_custom_data maps factor item and subj labels to schedule IDs without re-indexing", {
  des_t <- design_testlet(list(condition = c("a", "b")), testlet_size = 2L)
  dgp_t <- calibrate_dgp(
    des_t, c(a = 0.40, b = 0.65), subj_sd = 0.3, item_sd = 0.3,
    control = calibrate_control(testlet_sd = 0.25)
  )
  sp_t <- specify_target(dgp_t, "condition", scale = "logit")
  sched <- generate_design_schedule(des_t, n_subj = 8L, n_item = 8L)

  # Drop item 1, omit item_id, and drop unused factor levels
  tr_sub <- droplevels(sched$trials[sched$trials$item_id != 1L, ])
  tr_sub$item_id <- NULL
  tr_sub$y <- rep_len(c(0L, 1L), nrow(tr_sub))

  # Keeping original testlet column succeeds without spurious conflict
  val_keep_t <- validate_custom_data(tr_sub, sched)
  expect_identical(sort(unique(val_keep_t$item_id)), 2L:8L)
  expect_true(all(val_keep_t$testlet_id[val_keep_t$item_id == 2L] == 1L))
  expect_true(all(val_keep_t$testlet_id[val_keep_t$item_id == 3L] == 2L))

  # Omitting testlet and testlet_id still assigns item 2 -> passage 1 and item 3 -> passage 2
  tr_no_t <- tr_sub
  tr_no_t$testlet <- NULL
  tr_no_t$testlet_id <- NULL
  val_no_t <- validate_custom_data(tr_no_t, sched)
  expect_true(all(val_no_t$testlet_id[val_no_t$item_id == 2L] == 1L))
  expect_true(all(val_no_t$testlet_id[val_no_t$item_id == 3L] == 2L))

  # Unknown item label aborts cleanly
  tr_bad <- tr_no_t
  tr_bad$item <- factor("unknown_item")
  expect_error(validate_custom_data(tr_bad, sched), "not present in")

  ev_t <- evaluate_design_point(
    sp_t, n_subj = 8L, n_item = 8L, nsim = 2L, stream_seeds = 55L,
    data_fn = function(schedule, spec) {
      d <- droplevels(schedule$trials[schedule$trials$item_id != 1L, ])
      d$item_id <- NULL
      d$y <- stats::rbinom(nrow(d), 1L, 0.55)
      d
    }
  )
  expect_equal(ev_t$nsim, 2L)
})

test_that("isolated_seeded_call under Box-Muller preserves glmmplanr library path and .GlobalEnv data_fn bindings", {
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  had_g <- exists(".glmmplanr_bm_test_dgp", envir = .GlobalEnv, inherits = FALSE)
  if (had_g) old_g <- get(".glmmplanr_bm_test_dgp", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    if (had_g) assign(".glmmplanr_bm_test_dgp", old_g, envir = .GlobalEnv)
    else if (exists(".glmmplanr_bm_test_dgp", envir = .GlobalEnv, inherits = FALSE)) rm(".glmmplanr_bm_test_dgp", envir = .GlobalEnv)
  }, add = TRUE)

  des <- design_crossed(list(condition = c("a", "b")))
  dgp <- calibrate_dgp(des, c(a = 0.40, b = 0.65), subj_sd = 0.3, item_sd = 0.3)
  sp <- specify_target(dgp, "condition", scale = "pp", prob_se_method = "joint_delta")
  assign(".glmmplanr_bm_test_dgp", dgp, envir = .GlobalEnv)

  # Build data_fn whose enclosing environment is .GlobalEnv referencing `.glmmplanr_bm_test_dgp`
  fn_global <- eval(
    quote(function(schedule, spec) glmmplanr::simulate_dgp_data(.glmmplanr_bm_test_dgp, schedule = schedule)),
    envir = .GlobalEnv
  )

  suppressWarnings(set.seed(1234L, kind = "Marsaglia-Multicarry", normal.kind = "Box-Muller"))
  val_bm <- validate_design_point(
    sp, n_subj = 12L, n_item = 8L, nsim = 3L, seed = 42L,
    data_fn = fn_global, true_target = c(null = 25, planned = 25)
  )
  expect_equal(val_bm$n_valid, c(3L, 3L))

  plan_bm <- optimize_sample_size(
    sp, n_subj_range = 12L, n_item_range = 8L,
    method = "grid", nsim_final = 3L, seed = 42L,
    workers = 1L, verbose = FALSE, control = optimize_control(nsim_sensitivity = 0L)
  )
  expect_equal(plan_bm$search_meta$simulation_version, 5L)
})

test_that("extreme asymmetric contrast weights calibrate valid nulls or distinguish numerical precision boundaries", {
  des <- design_crossed(list(condition = c("a", "b")))
  dgp0 <- calibrate_dgp(des, c(a = 0.40, b = 0.60), subj_sd = 0.3, item_sd = 0.3, chance = 0)

  # (-1, 1e-13) with chance = 0 has an exact representable interior null (p_a ~ 5e-14, p_b ~ 0.5)
  sp_asym <- specify_target(
    dgp0, "condition", scale = "pp", criterion = "power",
    prob_se_method = "joint_delta",
    control = target_control(contrast_weights = c(a = -1, b = 1e-13))
  )
  expect_true(sp_asym$has_null_target)
  val_asym <- validate_design_point(sp_asym, n_subj = 12L, n_item = 8L, nsim = 3L, seed = 313L)
  expect_equal(nrow(val_asym), 2L)
  expect_identical(val_asym$condition, c("null", "planned"))
  expect_equal(val_asym$true_target[1L], 0.0, tolerance = 1e-20)

  # Sub-epsilon weight ratio against nonzero chance = 0.25 hits IEEE-754 double-precision boundary
  dgp_g <- calibrate_dgp(des, c(a = 0.45, b = 0.70), subj_sd = 0.3, item_sd = 0.3, chance = 0.25)
  st_num <- probability_contrast_null_status(dgp_g, c(a = -0.25, b = 1 - 1e-16), "pp", "marginal")
  expect_false(st_num$ok)
  expect_identical(st_num$kind, "numerical_boundary")
  expect_error(
    specify_target(
      dgp_g, "condition", scale = "pp", criterion = "power",
      control = target_control(contrast_weights = c(a = -0.25, b = 1 - 1e-16))
    ),
    "double precision"
  )

  # Mathematically unattainable zero still reports "outside the attainable probability contrast"
  dgp_3 <- calibrate_dgp(des, c(a = 0.45, b = 0.70), subj_sd = 0.3, item_sd = 0.3, chance = 1 / 3)
  expect_error(
    specify_target(
      dgp_3, "condition", scale = "pp", criterion = "power",
      control = target_control(contrast_weights = c(a = -1, b = 3))
    ),
    "outside the attainable probability contrast"
  )
})

test_that("plan_rasch_fixed_bank preserves item names in item_difficulties and item_summary", {
  bank_named <- c(walk = -1.0, read = 0.0, write = 1.0)
  fb <- plan_rasch_fixed_bank(
    item_difficulties = bank_named,
    n_subj = c(20L, 40L),
    max_error = c(difficulty = 0.60),
    nsim = 4L,
    seed = 515L
  )
  expect_identical(names(fb$item_difficulties), c("walk", "read", "write"))
  expect_equal(unname(fb$item_difficulties), c(-1.0, 0.0, 1.0))
  expect_true("item_name" %in% names(fb$item_summary))
  expect_identical(fb$item_summary$item_id, c(1L:3L, 1L:3L))
  expect_identical(fb$item_summary$item_name, rep(c("walk", "read", "write"), 2L))
  expect_s3_class(plot(fb, type = "items"), "ggplot")

  # Unnamed bank retains unnamed vector and character 1..J item_name
  fb_un <- plan_rasch_fixed_bank(c(-1.0, 0.0, 1.0), n_subj = 20L, max_error = c(difficulty = 0.60), nsim = 3L, seed = 515L)
  expect_null(names(fb_un$item_difficulties))
  expect_identical(fb_un$item_summary$item_name, c("1", "2", "3"))

  # Partial or duplicate item names are rejected
  expect_error(
    plan_rasch_fixed_bank(c(walk = -1.0, 0.0, write = 1.0), n_subj = 20L, max_error = c(difficulty = 0.60), nsim = 2L),
    "unique[[:space:]]+and[[:space:]]+non-empty"
  )
  expect_error(
    plan_rasch_fixed_bank(c(walk = -1.0, walk = 1.0), n_subj = 20L, max_error = c(difficulty = 0.60), nsim = 2L),
    "unique[[:space:]]+and[[:space:]]+non-empty"
  )

  # Export contract check: cell_contrast_weights is exported and align_cell_weights is internal
  exports <- getNamespaceExports("glmmplanr")
  expect_true("cell_contrast_weights" %in% exports)
  expect_false("align_cell_weights" %in% exports)
  expect_false("align_named_cell_vector" %in% exports)
})

test_that("mirai pool detects workers loaded with a mismatched glmmplanr build and matches sequential when aligned", {
  skip_on_cran()
  on.exit(try(mirai::daemons(0L), silent = TRUE), add = TRUE)

  des <- design_crossed(list(condition = c("a", "b")))
  dgp <- calibrate_dgp(des, c(a = 0.45, b = 0.65), subj_sd = 0.35, item_sd = 0.30)
  sp <- specify_target(
    dgp, "condition", scale = "pp", criterion = "power",
    prob_se_method = "joint_delta",
    control = target_control(contrast_weights = 1e-7 * c(a = -1, b = 1))
  )

  # Aligned fresh daemons produce identical results to sequential execution
  mirai::daemons(2L)
  plan_seq <- optimize_sample_size(
    sp, n_subj_range = c(12L, 16L), n_item_range = 8L,
    method = "grid", nsim_final = 4L, seed = 606L,
    workers = 1L, verbose = FALSE,
    control = optimize_control(nsim_sensitivity = 0L, parallel_backend = "sequential")
  )
  plan_par <- optimize_sample_size(
    sp, n_subj_range = c(12L, 16L), n_item_range = 8L,
    method = "grid", nsim_final = 4L, seed = 606L,
    workers = 2L, verbose = FALSE,
    control = optimize_control(nsim_sensitivity = 0L, parallel_backend = "mirai")
  )
  expect_identical(plan_seq$verified$decision_status, plan_par$verified$decision_status)
  expect_equal(plan_seq$verified$power, plan_par$verified$power, tolerance = 1e-12)

  # Simulate pre-existing workers that loaded an older/mismatched glmmplanr build
  mirai::everywhere({
    ns <- asNamespace("glmmplanr")
    unlockBinding(".glmmplanr_source_md5", ns)
    assign(".glmmplanr_source_md5", c(stale = "0000"), envir = ns)
    lockBinding(".glmmplanr_source_md5", ns)
  })[]
  expect_error(
    optimize_sample_size(
      sp, n_subj_range = c(12L, 16L), n_item_range = 8L,
      method = "grid", nsim_final = 4L, seed = 606L,
      workers = 2L, verbose = FALSE,
      control = optimize_control(nsim_sensitivity = 0L, parallel_backend = "mirai")
    ),
    "different.*glmmplanr.*mirai::daemons\\(0\\)"
  )
  mirai::daemons(0L)
})

test_that("named random-effect SD vectors and correlation matrices align by name rather than position", {
  des1 <- design_crossed(list(condition = c("a", "b")))
  dgp_fwd <- calibrate_dgp(
    des1, c(a = 0.45, b = 0.65),
    subj_sd = c("(Intercept)" = 0.05, condition = 0.95),
    item_sd = c("(Intercept)" = 0.30, condition1 = 0.20)
  )
  dgp_rev <- calibrate_dgp(
    des1, c(a = 0.45, b = 0.65),
    subj_sd = c(condition = 0.95, "(Intercept)" = 0.05),
    item_sd = c(condition1 = 0.20, "(Intercept)" = 0.30)
  )
  expect_equal(dgp_fwd$tau_subj, dgp_rev$tau_subj, tolerance = 1e-12)
  expect_equal(dgp_fwd$tau_item, dgp_rev$tau_item, tolerance = 1e-12)
  expect_equal(dgp_fwd$Sigma_subj, dgp_rev$Sigma_subj, tolerance = 1e-12)
  expect_equal(dgp_fwd$beta, dgp_rev$beta, tolerance = 1e-12)

  # Identical evaluation results regardless of input order of named subj_sd
  sp_fwd <- specify_target(dgp_fwd, "condition", scale = "pp", criterion = "ci_width", ci_width = 21, prob_se_method = "joint_delta")
  sp_rev <- specify_target(dgp_rev, "condition", scale = "pp", criterion = "ci_width", ci_width = 21, prob_se_method = "joint_delta")
  ev_fwd <- evaluate_design_point(sp_fwd, n_subj = 16L, n_item = 8L, nsim = 8L, stream_seeds = 707L)
  ev_rev <- evaluate_design_point(sp_rev, n_subj = 16L, n_item = 8L, nsim = 8L, stream_seeds = 707L)
  expect_equal(ev_fwd$aipe_prob, ev_rev$aipe_prob, tolerance = 1e-12)
  expect_identical(ev_fwd$decision_status, ev_rev$decision_status)

  # Named correlation matrix aligns row/col names to design RE order
  des2 <- design_crossed(list(A = c("a1", "a2"), B = c("b1", "b2")))
  # Canonical RE order: (Intercept), A1, B1, A1:B1
  # Supply matrix in order: (Intercept), B, A, B:A with corr((Intercept), A) = 0.8 and corr((Intercept), B) = -0.2
  nms_perm <- c("(Intercept)", "B", "A", "B:A")
  R_perm <- diag(1.0, 4L)
  dimnames(R_perm) <- list(nms_perm, nms_perm)
  R_perm["(Intercept)", "A"] <- R_perm["A", "(Intercept)"] <- 0.80
  R_perm["(Intercept)", "B"] <- R_perm["B", "(Intercept)"] <- -0.20

  dgp_corr <- calibrate_dgp(
    des2,
    c(a1_b1 = 0.40, a1_b2 = 0.50, a2_b1 = 0.55, a2_b2 = 0.65),
    subj_sd = c("(Intercept)" = 0.5, B = 0.2, A = 0.4, "A:B" = 0.1),
    item_sd = 0.3,
    corr_structure = list(subj = R_perm, item = "zero")
  )
  expect_equal(dgp_corr$R_subj["(Intercept)", "A1"], 0.80, tolerance = 1e-12)
  expect_equal(dgp_corr$R_subj["(Intercept)", "B1"], -0.20, tolerance = 1e-12)
  expect_equal(unname(dgp_corr$tau_subj[c("(Intercept)", "A1", "B1", "A1:B1")]), c(0.5, 0.4, 0.2, 0.1))

  # Partial, duplicate, or mismatched names are rejected
  expect_error(
    calibrate_dgp(des1, c(0.45, 0.65), subj_sd = c(condition = 0.5), item_sd = 0.3),
    "Intercept"
  )
  expect_error(
    calibrate_dgp(des1, c(0.45, 0.65), subj_sd = c(condition = 0.5, 0.2), item_sd = 0.3),
    "name every element"
  )
  R_bad <- R_perm
  colnames(R_bad) <- c("(Intercept)", "A", "B", "A:B")
  expect_error(
    calibrate_dgp(des2, c(0.40, 0.50, 0.55, 0.65), subj_sd = 0.4, item_sd = 0.3, corr_structure = list(subj = R_bad, item = "zero")),
    "same order"
  )
})

test_that("simulate_dgp_data rejects schedules from different allocations while allowing re-seeded partial schedules", {
  des_ls <- design_crossed(list(condition = c("a", "b")), allocation = "latin_square")
  des_bs <- design_crossed(list(condition = c("a", "b")), allocation = "between_subj")
  des_cr <- design_crossed(list(condition = c("a", "b")), allocation = "crossed")
  dgp_ls <- calibrate_dgp(des_ls, c(a = 0.40, b = 0.65), subj_sd = 0.3, item_sd = 0.3)

  sched_bs <- generate_design_schedule(des_bs, 16L, 8L)
  sched_cr <- generate_design_schedule(des_cr, 16L, 8L)
  expect_error(simulate_dgp_data(dgp_ls, schedule = sched_bs), "Allocation|within-subject")
  expect_error(simulate_dgp_data(dgp_ls, schedule = sched_cr), "Allocation|scheduled trial")

  # Even if allocation metadata is stripped, trial-level structure catches both mismatches
  sched_bs_strip <- sched_bs
  sched_bs_strip$allocation <- NULL
  sched_bs_strip$between_subj <- NULL
  expect_error(simulate_dgp_data(dgp_ls, schedule = sched_bs_strip), "between-subject")

  # Same partial-presentation design with a different item_seed remains compatible
  des_p1 <- design_crossed(list(condition = c("a", "b")), items_per_subject = 4L, item_seed = 1L)
  des_p2 <- design_crossed(list(condition = c("a", "b")), items_per_subject = 4L, item_seed = 99L)
  dgp_p1 <- calibrate_dgp(des_p1, c(a = 0.40, b = 0.65), subj_sd = 0.3, item_sd = 0.3)
  sched_p2 <- generate_design_schedule(des_p2, 16L, 8L)
  df_p <- simulate_dgp_data(dgp_p1, schedule = sched_p2)
  expect_equal(nrow(df_p), 64L)
})

test_that("high guessing floors chance >= 0.95 initialize and converge cleanly in native C++ Laplace", {
  des <- design_crossed(list(condition = c("a", "b")))
  for (g_hi in c(0.95, 0.96, 0.995)) {
    p_a <- g_hi + 0.35 * (1.0 - g_hi)
    p_b <- g_hi + 0.70 * (1.0 - g_hi)
    dgp_hi <- calibrate_dgp(des, c(a = p_a, b = p_b), chance = g_hi, subj_sd = 0.3, item_sd = 0.3)
    sp_hi <- specify_target(dgp_hi, "condition", scale = "logit", engine = "bayes_laplace")
    ev_hi <- evaluate_design_point(sp_hi, n_subj = 64L, n_item = 24L, nsim = 6L, stream_seeds = 818L)
    expect_equal(ev_hi$n_valid, 6L)
  }

  stan_lines <- readLines(system.file("stan", "crossed_binomial_mafc.stan", package = "glmmplanr"))
  expect_true(any(grepl("real<lower=0,\\s*upper=1>\\s*chance;", stan_lines)))
})

test_that("cell_probs_from_diffs rejects duplicate, conflicting, or (Intercept) effect specifications", {
  des <- design_crossed(list(C = c("c0", "c1")))
  expect_error(cell_probs_from_diffs(des, 0.50, c(C = 30, C1 = 8)), "Conflicting")
  expect_error(cell_probs_from_diffs(des, 0.50, c(C1 = 8, C = 30)), "Conflicting")
  expect_error(cell_probs_from_diffs(des, 0.50, c(C = 8, C = 30)), "Duplicate")
  expect_error(cell_probs_from_diffs(des, 0.50, c("(Intercept)" = 70)), "baseline_prob")
  expect_equal(unname(cell_probs_from_diffs(des, 0.50, c(C = 10))), c(0.45, 0.55), tolerance = 1e-12)
})

test_that("conditional failure accounting returns NA intervals and MCSE when n_valid == 0", {
  des <- design_crossed(list(condition = c("a", "b")), participant_dropout = 1 - 1e-12)
  dgp <- calibrate_dgp(des, c(a = 0.45, b = 0.65), subj_sd = 0.3, item_sd = 0.3)
  sp_cond <- specify_target(
    dgp, "condition", scale = "pp", criterion = "both", ci_width = 15,
    prob_se_method = "joint_delta",
    control = target_control(failure_policy = "conditional")
  )
  ev0 <- evaluate_design_point(sp_cond, n_subj = 12L, n_item = 8L, nsim = 10L, stream_seeds = 11L)
  expect_equal(ev0$n_valid, 0L)
  expect_true(is.na(ev0$power))
  expect_true(is.na(ev0$power_mcse))
  expect_true(is.na(ev0$power_mc_lo))
  expect_true(is.na(ev0$power_mc_hi))
  expect_true(is.na(ev0$aipe_prob))
  expect_true(is.na(ev0$aipe_mc_lo))
  expect_true(is.na(ev0$aipe_mc_hi))
  expect_true(is.na(ev0$joint_prob))
  expect_true(is.na(ev0$joint_mc_lo))
  expect_true(is.na(ev0$joint_mc_hi))
  expect_true(is.na(ev0$rope_excl_prob))
  expect_true(is.na(ev0$rope_mc_lo))
  expect_true(is.na(ev0$rope_mc_hi))
  expect_equal(ev0$effective_power, 0.0)
  expect_equal(ev0$effective_aipe_prob, 0.0)
  expect_equal(ev0$effective_rope_prob, 0.0)
  expect_identical(ev0$decision_status, "fail")
})

test_that("validation plot labels zero-truth planned condition as Type I error rather than 80% Power", {
  des <- design_crossed(list(condition = c("a", "b")))
  dgp_zero <- calibrate_dgp(des, c(a = 0.60, b = 0.60), subj_sd = 0.3, item_sd = 0.3)
  sp_zero <- specify_target(dgp_zero, "condition", scale = "pp", criterion = "power", prob_se_method = "joint_delta")
  val_zero <- validate_design_point(sp_zero, n_subj = 16L, n_item = 8L, nsim = 6L, seed = 44L)
  p_zero <- plot(val_zero)
  panels <- unique(as.character(p_zero$data$panel))
  expect_false("Power" %in% panels)
  expect_true("Type I error" %in% panels)
  t_rows <- p_zero$data[p_zero$data$panel == "Type I error", , drop = FALSE]
  expect_equal(nrow(t_rows), 2L)
  expect_true(all(t_rows$target == sp_zero$alpha))
})

test_that("as_gt displays ROPE exclusion column for ROPE plans and rejects uncomputed sensitivity/robust views", {
  des <- design_crossed(list(condition = c("a", "b")))
  dgp <- calibrate_dgp(des, c(a = 0.30, b = 0.70), subj_sd = 0.3, item_sd = 0.3)
  sp_rope <- specify_target(
    dgp, "condition", scale = "pp", criterion = "rope",
    prob_se_method = "joint_delta",
    control = target_control(rope = c(-60, 60))
  )
  plan_rope <- optimize_sample_size(
    sp_rope, n_subj_range = 16L, n_item_range = 8L,
    method = "grid", nsim_final = 8L, seed = 321L,
    workers = 1L, verbose = FALSE, control = optimize_control(nsim_sensitivity = 0L)
  )
  gt_rope <- as_gt(plan_rope)
  expect_true("rope_ci" %in% names(gt_rope$`_data`))
  html_rope <- gt::as_raw_html(gt_rope)
  expect_match(html_rope, "ROPE Excl.", fixed = TRUE)

  expect_error(as_gt(plan_rope, which = "sensitivity"), "not run")
  expect_error(as_gt(plan_rope, which = "robust"), "not run")
})

test_that("ci_width_q80 stays fixed at the 80th percentile across ci_assurance settings", {
  des <- design_crossed(list(condition = c("a", "b")))
  dgp <- calibrate_dgp(des, c(a = 0.45, b = 0.65), subj_sd = 0.35, item_sd = 0.30)
  q80s <- numeric(3L)
  qass <- numeric(3L)
  ass_vals <- c(0.50, 0.80, 0.95)
  for (k in seq_along(ass_vals)) {
    sp_k <- specify_target(
      dgp, "condition", scale = "pp", criterion = "ci_width", ci_width = 30,
      prob_se_method = "joint_delta",
      control = target_control(ci_assurance = ass_vals[k])
    )
    ev_k <- evaluate_design_point(sp_k, n_subj = 16L, n_item = 8L, nsim = 16L, stream_seeds = 999L)
    q80s[k] <- ev_k$ci_width_q80
    qass[k] <- ev_k$ci_width_q_assurance
  }
  expect_equal(q80s[1L], q80s[2L], tolerance = 1e-12)
  expect_equal(q80s[2L], q80s[3L], tolerance = 1e-12)
  expect_lt(qass[1L], qass[2L])
  expect_lt(qass[2L], qass[3L])
  expect_equal(q80s[2L], qass[2L], tolerance = 1e-12)
})

test_that("switch_preference warns on legacy pre-v4 saved plans, Rasch plans retain DGP/allocation settings, Box-Muller captures formal defaults, and invalid dimensions abort", {
  des <- design_crossed(list(condition = c("a", "b")))
  dgp <- calibrate_dgp(des, c(a = 0.40, b = 0.65), subj_sd = 0.3, item_sd = 0.3)
  sp_w <- specify_target(
    dgp, "condition", scale = "pp", criterion = "ci_width",
    ci_width = 40 * 1e-7, prob_se_method = "joint_delta",
    control = target_control(contrast_weights = 1e-7 * c(a = -1, b = 1))
  )
  plan_cur <- optimize_sample_size(
    sp_w, n_subj_range = c(12L, 16L), n_item_range = 8L,
    method = "grid", nsim_final = 4L, seed = 111L,
    workers = 1L, verbose = FALSE, control = optimize_control(nsim_sensitivity = 0L)
  )
  plan_old <- plan_cur
  plan_old$search_meta$simulation_version <- 3L
  expect_warning(
    switch_preference(plan_old, "fewer_subj"),
    "earlier.*simulation_version = 3|unnormalized contrast weights"
  )

  # Rasch plans retain generating/allocation/missingness parameters
  des_r <- design_crossed(list(), allocation = "crossed", items_per_subject = 4L, item_seed = 77L,
                          attrition = 0.10, participant_dropout = 0.25)
  dgp_r <- calibrate_dgp(des_r, 0.55, subj_sd = 0.7, item_sd = 0.5)
  rp <- plan_rasch_precision(dgp_r, data.frame(n_subj = 24L, n_item = 8L), c(ability = 0.8), nsim = 3L, seed = 12L)
  expect_identical(rp$items_per_subject, 4L)
  expect_identical(rp$item_seed, 77L)
  expect_equal(rp$attrition, 0.10)
  expect_equal(rp$participant_dropout, 0.25)
  expect_s3_class(rp$dgp, "glmmplanr_dgp")

  fb <- plan_rasch_fixed_bank(
    seq(-1, 1, length.out = 8L), n_subj = 24L, max_error = c(difficulty = 0.8),
    items_per_subject = 4L, item_seed = 88L, attrition = 0.05, participant_dropout = 0.20,
    nsim = 3L, seed = 12L
  )
  expect_identical(fb$items_per_subject, 4L)
  expect_identical(fb$item_seed, 88L)
  expect_equal(fb$attrition, 0.05)
  expect_equal(fb$participant_dropout, 0.20)
  expect_match(report_plan(fb, print_output = FALSE), "participant dropout = 20.0%", fixed = TRUE)

  # Box-Muller isolated_seeded_call captures globals referenced only in default argument expressions
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    if (exists(".glmmplanr_default_arg_p", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".glmmplanr_default_arg_p", envir = .GlobalEnv)
    }
  }, add = TRUE)

  assign(".glmmplanr_default_arg_p", 0.62, envir = .GlobalEnv)
  fn_def <- eval(
    quote(function(schedule, spec, p = .glmmplanr_default_arg_p) {
      d <- schedule$trials
      d$y <- stats::rbinom(nrow(d), 1L, p)
      d
    }),
    envir = .GlobalEnv
  )
  suppressWarnings(set.seed(314L, kind = "Marsaglia-Multicarry", normal.kind = "Box-Muller"))
  val_def <- validate_design_point(
    sp_w, n_subj = 12L, n_item = 8L, nsim = 2L, seed = 19L,
    data_fn = fn_def, true_target = c(null = 0, planned = 0)
  )
  expect_equal(val_def$n_valid, c(2L, 2L))

  # Negative and non-integer design dimensions abort before rounding/coercion
  expect_error(design_testlet(list(condition = c("a", "b")), testlet_size = 2.9), "whole integer")
  expect_error(design_crossed(list(condition = c("a", "b")), items_per_subject = 2.5), "integer")
  expect_error(valid_multiples(-4, 2), "positive whole integers")
  expect_error(valid_multiples(c(-4, -6), 2), "positive whole integers")
  expect_error(valid_multiples(2.5, 2), "positive whole integers")
  expect_error(
    optimize_sample_size(sp_w, n_subj_range = -4, n_item_range = -6, method = "grid", nsim_final = 2L, verbose = FALSE),
    "positive whole integers"
  )
})

test_that("v4 positive-chance plans warn in switch_preference, crossed duplicate conditions abort, and mafc_binomial round-trips without truncation", {
  # 1. simulation_version = 5L distinguishes pre-fix v4 positive-chance plans
  des_hi <- design_crossed(list(condition = c("A", "B")))
  dgp_hi <- calibrate_dgp(des_hi, c(A = 0.9675, B = 0.985), chance = 0.95, subj_sd = 0.3, item_sd = 0.3)
  sp_hi <- specify_target(
    dgp_hi, "condition", scale = "pp", criterion = "ci_width", ci_width = 5,
    prob_se_method = "joint_delta", engine = "bayes_laplace"
  )
  plan_v5 <- optimize_sample_size(
    sp_hi, n_subj_range = 64L, n_item_range = 24L,
    method = "grid", nsim_final = 24L, seed = 2026L,
    workers = 1L, verbose = FALSE, control = optimize_control(nsim_sensitivity = 0L)
  )
  expect_equal(plan_v5$search_meta$simulation_version, 5L)
  expect_equal(plan_v5$verified$n_valid, 24L)
  expect_identical(plan_v5$verified$decision_status, "pass")

  saved_v4_hi <- plan_v5
  saved_v4_hi$search_meta$simulation_version <- 4L
  saved_v4_hi$verified$n_valid <- 0L
  saved_v4_hi$verified$decision_status <- "fail"
  saved_v4_hi$decision_status <- "fail"
  saved_v4_hi$feasible <- FALSE
  expect_warning(
    switch_preference(saved_v4_hi, "fewer_subj"),
    "simulation_version = 4.*positive-chance|simulation_version < 5"
  )
  expect_error(
    topup_verification(saved_v4_hi, nsim_extra = 8L, verbose = FALSE),
    "simulation_version < 5"
  )

  # 2. Crossed 8x6x2 schedule where one (subj_id, item_id) pair has A,B replaced by A,A aborts
  des_cr <- design_crossed(list(condition = c("A", "B")), allocation = "crossed")
  dgp_cr <- calibrate_dgp(des_cr, c(A = 0.45, B = 0.65), subj_sd = 0.3, item_sd = 0.3)
  sched_cr <- generate_design_schedule(des_cr, n_subj = 8L, n_item = 6L)
  expect_equal(nrow(sched_cr$trials), 96L)

  bad_pair <- sched_cr
  idx_p11 <- which(bad_pair$trials$subj_id == 1L & bad_pair$trials$item_id == 1L)
  expect_equal(length(idx_p11), 2L)
  bad_pair$trials$cell_id[idx_p11] <- 1L
  bad_pair$trials$condition[idx_p11] <- factor(c("A", "A"), levels = c("A", "B"))
  expect_error(
    simulate_dgp_data(dgp_cr, schedule = bad_pair),
    "repeats the same experimental condition"
  )

  # Latin-square schedule with an unbalanced condition swap also aborts
  des_ls <- design_crossed(list(condition = c("A", "B")), allocation = "latin_square")
  dgp_ls <- calibrate_dgp(des_ls, c(A = 0.45, B = 0.65), subj_sd = 0.3, item_sd = 0.3)
  bad_ls <- generate_design_schedule(des_ls, n_subj = 8L, n_item = 6L)
  i_b <- which(bad_ls$trials$subj_id == 1L & bad_ls$trials$condition == "B")[1L]
  bad_ls$trials$cell_id[i_b] <- 1L
  bad_ls$trials$condition[i_b] <- factor("A", levels = c("A", "B"))
  expect_error(
    simulate_dgp_data(dgp_ls, schedule = bad_ls),
    "unbalanced"
  )

  # 3. mafc_binomial preserves round-trip precision at chance = 1 - 1e-8 and eta = -18
  fam_hi <- mafc_binomial(1 - 1e-8)
  rt0 <- fam_hi$linkfun(fam_hi$linkinv(0))
  expect_true(is.finite(rt0))
  expect_equal(rt0, 0, tolerance = 1e-7)
  expect_equal(fam_hi$linkfun(fam_hi$linkinv(c(-6, 0, 6))), c(-6, 0, 6), tolerance = 1e-4)
  expect_true(all(fam_hi$linkinv(c(-100, 0, 100)) > (1 - 1e-8)))
  expect_true(all(fam_hi$linkinv(c(-100, 0, 100)) < 1))

  fam_3 <- mafc_binomial(1 / 3)
  expect_equal(fam_3$linkfun(fam_3$linkinv(-18)), -18, tolerance = 1e-7)
  expect_equal(fam_3$linkfun(fam_3$linkinv(c(-25, -18, 0, 18, 25))), c(-25, -18, 0, 18, 25), tolerance = 1e-4)

  # External glm fit succeeds cleanly with high guessing floor
  df_glm <- simulate_dgp_data(dgp_hi, n_subj = 32L, n_item = 12L)
  fit_glm <- stats::glm(y ~ condition, family = mafc_binomial(0.95), data = df_glm)
  expect_true(all(is.finite(stats::coef(fit_glm))))
})
