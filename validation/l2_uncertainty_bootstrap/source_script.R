# Mechanism diagnostics on preserved trials; oracle covariance is NOT a usable estimator.
source("validation/planning_comparison_helpers.R")
args <- commandArgs(TRUE)
mode <- if (length(args)) match.arg(args[1L], c("audit", "smoke", "oracle", "bootstrap")) else "audit"
out <- if (length(args) > 1L) args[2L] else paste0("validation/l2_uncertainty_", mode)
workers <- if (length(args) > 2L) as.integer(args[3L]) else 6L
if (.Platform$OS.type == "windows") workers <- 1L
source_run <- "validation/l2_replication_large"
# Reuse only the frozen DGP definitions, not the original runner's execution.
for (e in parse(file.path(source_run, "source_script.R"))) {
  if (is.call(e) && identical(e[[1L]], as.name("<-")) && is.symbol(e[[2L]]) &&
      as.character(e[[2L]]) %in% c("profiles", "make_spec")) eval(e)
}
cases <- read.csv(file.path(source_run, "confirmation_scenarios.csv"))
chosen <- data.frame(profile = c("published_8pp", "slopes_8pp", "conservative_5pp",
  rep("high_variance_8pp", 3)), n_subj = c(54L, 72L, 162L, 72L, 162L, 216L),
  n_item = c(30L, 60L, 60L, 90L, 30L, 90L))
chosen <- merge(chosen, cases, by = c("profile", "n_subj", "n_item"))
stopifnot(nrow(chosen) == 6L)
nsim <- switch(mode, audit = 2048L, smoke = 3L, oracle = 512L, bootstrap = 128L)
B <- if (mode == "bootstrap") 399L else 0L
checkpoint_dir(out, c("validation/l2_uncertainty_diagnostics.R",
  file.path(source_run, c("source_script.R", "fits.csv.gz", "manifest.rds")),
  "R/marginal_scale.R", "R/calibrate.R", "R/utils_math.R"),
  list(mode = mode, chosen = chosen, nsim = nsim, bootstrap = B,
    bootstrap_seed = 910100001L, replay_tolerance_pp = 1e-7), c("glmmplanr", "digest"))
write.csv(chosen, file.path(out, "cases.csv"), row.names = FALSE)
saved <- read.csv(file.path(source_run, "fits.csv.gz"))
stopifnot(nrow(saved) == 245760L, all(saved$valid),
  !anyDuplicated(saved[c("case", "condition", "replicate")]))
summarize <- function(v) {
  ok <- v$valid & is.finite(v$estimate) & is.finite(v$se) & v$se > 0
  z <- v[ok, ]; err <- z$estimate - z$truth; t <- err/z$se
  ci <- function(k, n) unname(glmmplanr:::wilson_mc_ci(k, n))
  cover <- abs(t) <= qnorm(.975)
  detect <- abs(z$estimate/z$se) > qnorm(.975)
  joint <- detect & 2*qnorm(.975)*z$se <= 10
  cc <- ci(sum(cover), sum(ok)); dc <- ci(sum(detect), nrow(v)); jc <- ci(sum(joint), nrow(v))
  data.frame(nsim = nrow(v), n_valid = sum(ok), bias_pp = mean(err),
    empirical_sd_pp = sd(err), rms_se_pp = sqrt(mean(z$se^2)),
    coverage = mean(cover), coverage_lo = cc[1], coverage_hi = cc[2],
    detection = sum(detect)/nrow(v), detection_lo = dc[1], detection_hi = dc[2],
    joint = sum(joint)/nrow(v), joint_lo = jc[1], joint_hi = jc[2],
    joint_decision = if (jc[1] >= .8) "pass" else if (jc[2] < .8) "fail" else "indeterminate",
    mean_width_pp = mean(2*qnorm(.975)*z$se),
    z_mean = mean(t), z_sd = sd(t), z_q025 = unname(quantile(t,.025)),
    z_q975 = unname(quantile(t,.975)),
    lower_miss = mean(t < -qnorm(.975)), upper_miss = mean(t > qnorm(.975)),
    z_excess_kurtosis = mean((t-mean(t))^4)/mean((t-mean(t))^2)^2-3,
    error_squared_se_squared_correlation = cor(err^2, z$se^2),
    bias_removed_coverage_diagnostic = mean(abs(err-mean(err)) <= qnorm(.975)*z$se),
    fixed_empirical_se_coverage_diagnostic = mean(abs(err) <= qnorm(.975)*sd(err)))
}
if (mode == "audit") {
  summary <- do.call(rbind, lapply(split(saved, interaction(saved$case,saved$condition)), function(v)
    cbind(cases[cases$case == v$case[1], ], condition = v$condition[1], summarize(v))))
  quintiles <- do.call(rbind, lapply(split(saved, interaction(saved$case,saved$condition)), function(v) {
    v$quintile <- pmin(5L, ceiling(rank(v$se, ties.method = "first")/nrow(v)*5))
    do.call(rbind, lapply(split(v,v$quintile), function(w)
      data.frame(case = w$case[1], condition = w$condition[1], quintile = w$quintile[1],
        nsim = nrow(w), se_min = min(w$se), se_max = max(w$se),
        rms_error = sqrt(mean((w$estimate-w$truth)^2)), rms_se = sqrt(mean(w$se^2)),
        coverage = mean(w$covered))))
  }))
  previous <- read.csv(file.path(source_run,"diagnostic_summary.csv"))
  matched <- merge(summary,previous,by=c("case","condition"),suffixes=c("_new","_old"))
  for (nm in c("bias_pp","empirical_sd_pp","rms_se_pp","coverage","detection"))
    stopifnot(max(abs(matched[[paste0(nm,"_new")]]-matched[[paste0(nm,"_old")]])) < 1e-10)
  write.csv(summary,file.path(out,"summary.csv"),row.names=FALSE)
  write.csv(quintiles,file.path(out,"se_quintiles.csv"),row.names=FALSE)
  print(subset(summary,case %in% chosen$case,
    c(profile,n_subj,n_item,condition,bias_pp,rms_se_pp,coverage,z_sd,z_q025,z_q975,
      bias_removed_coverage_diagnostic,fixed_empirical_se_coverage_diagnostic)),row.names=FALSE)
} else {
  jobs <- merge(chosen,data.frame(condition=c("null","planned")))
  jobs$job <- seq_len(nrow(jobs))
  run <- function(k) {
    x <- jobs[k, ]; path <- file.path(out,sprintf("case-%02d.rds",k))
    if (file.exists(path)) return(readRDS(path)$trials)
    spec <- make_spec(profiles[profiles$profile == x$profile, ])
    if (x$condition == "null") spec <- glmmplanr:::null_target_spec(spec)
    schedule <- generate_design_schedule(spec$dgp$design,x$n_subj,x$n_item)
    cache <- glmmplanr:::prepare_cpp_schedule_cache(spec,schedule)
    ws <- glmmplanr:::get_cached_cpp_ws(cache,spec$analysis_re)
    bases <- glmmplanr:::allocate_stage_point_streams(810900001L+10000L*x$case,
      stage_slot=1L,grid_ids=1:2,stride=2L)
    streams <- glmmplanr:::derive_rep_substreams(bases[[if(x$condition=="null") 1 else 2]],nsim)
    old <- subset(saved,case==x$case & condition==x$condition & replicate<=nsim)
    old <- old[order(old$replicate), ]; stopifnot(nrow(old)==nsim)
    boot_streams <- glmmplanr:::make_lecuyer_streams(910100001L+10000L*k,nsim)
    diagnostics <- list(); rows <- list()
    for (i in seq_len(nsim)) {
      sim <- glmmplanr:::with_rng_stream(streams[[i]],glmmplanr:::simulate_fast_y_keep(cache))
      current <- glmmplanr:::fit_via_cpp_cached(sim,cache,spec,TRUE,spec$analysis_re,FALSE)
      stopifnot(current$converged,abs(100*current$prob_est-old$estimate[i])<1e-7,
        abs(100*current$prob_se-old$se[i])<1e-7)
      Ss <- spec$dgp$Sigma_subj; Si <- spec$dgp$Sigma_item
      if (spec$analysis_re=="intercept_only") { Ss <- Ss[1,1,drop=FALSE]; Si <- Si[1,1,drop=FALSE] }
      fixed <- glmmplanr:::cpp_fit_crossed_glmm_ws(ws,sim$y,chance=0,corr_mode=0L,
        regularized=TRUE,lkj_eta=spec$lkj_eta,prior_beta_sd=spec$prior_beta_sd,
        prior_tau_sd=spec$prior_tau_sd,max_em_iter=160L,max_pirls_iter=22L,tol=1.5e-4,
        init_Sigma_subj=Ss,init_Sigma_item=Si,fix_cov=TRUE,compute_joint_vcov=TRUE)
      delta <- glmmplanr:::estimate_prob_effect_delta(fixed$beta,fixed$vcov,Ss,Si,NULL,
        spec$dgp$design,0,spec$contrast_weights,"marginal",spec$analysis_re,spec$dgp$n_quad)
      rows[[length(rows)+1L]] <- data.frame(x,replicate=i,method="current",estimate=100*current$prob_est,
        se=100*current$prob_se,valid=current$converged,truth=old$truth[i],
        response_hash=digest::digest(sim$y,algo="xxhash64"))
      rows[[length(rows)+1L]] <- data.frame(x,replicate=i,method="fixed_true_covariance",
        estimate=100*delta$est,se=100*delta$se,valid=fixed$converged,truth=old$truth[i],
        response_hash=digest::digest(sim$y,algo="xxhash64"))
      if (B>0L) {
        boot_spec <- spec; boot_spec$prob_se_method <- "bootstrap"; boot_spec$prob_boot <- B
        boot <- glmmplanr:::with_rng_stream(boot_streams[[i]],
          glmmplanr:::fit_via_cpp_cached(sim,cache,boot_spec,TRUE,spec$analysis_re,FALSE))
        stopifnot(abs(boot$prob_est-current$prob_est)<1e-9)
        rows[[length(rows)+1L]] <- data.frame(x,replicate=i,method="bootstrap",
          estimate=100*boot$prob_est,se=100*boot$prob_se,valid=boot$converged && is.finite(boot$prob_se),
          truth=old$truth[i],response_hash=digest::digest(sim$y,algo="xxhash64"))
      }
      diagnostics[[i]] <- list(replicate=i,current=current,fixed=fixed)
    }
    trials <- do.call(rbind,rows)
    atomic_rds(list(trials=trials,diagnostics=diagnostics),path)
    cat(sprintf("Uncertainty %s: case %d %s, %d paired trials complete\n",mode,x$case,x$condition,nsim))
    trials
  }
  rows <- parallel::mclapply(seq_len(nrow(jobs)),run,mc.cores=workers,mc.set.seed=FALSE,
    mc.preschedule=FALSE)
  stopifnot(all(vapply(rows,is.data.frame,logical(1))))
  trials <- do.call(rbind,rows)
  stopifnot(!anyDuplicated(trials[c("case","condition","replicate","method")]))
  summary <- do.call(rbind,lapply(split(trials,interaction(trials$case,trials$condition,trials$method)),
    function(v) cbind(v[1,c("profile","n_subj","n_item","case","condition","method")],summarize(v))))
  write.csv(trials,gzfile(file.path(out,"fits.csv.gz")),row.names=FALSE)
  write.csv(summary,file.path(out,"summary.csv"),row.names=FALSE)
  print(summary[c("profile","n_subj","n_item","condition","method","n_valid","coverage","detection","rms_se_pp")],row.names=FALSE)
}
file.copy("validation/l2_uncertainty_diagnostics.R",file.path(out,"source_script.R"),overwrite=TRUE)
