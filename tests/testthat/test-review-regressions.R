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
  withr::local_options(sass.cache=FALSE)
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
