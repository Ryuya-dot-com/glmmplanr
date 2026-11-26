# Same-data replay and independent confirmation of a native covariance-update candidate.
source("validation/planning_comparison_helpers.R")
args<-commandArgs(TRUE)
mode<-if(length(args)) match.arg(args[1],c("smoke","paired","confirm")) else "smoke"
out<-if(length(args)>1) args[2] else paste0("validation/l2_uncertainty_confirmation_",mode)
variant<-if(length(args)>2) args[3] else "candidate"
workers<-if(length(args)>3) as.integer(args[4]) else 6L
if(.Platform$OS.type=="windows") workers<-1L
for(e in parse("validation/l2_replication_large/source_script.R")) {
  if(is.call(e) && identical(e[[1]],as.name("<-")) && is.symbol(e[[2]]) &&
     as.character(e[[2]]) %in% c("profiles","make_spec")) eval(e)
}
chosen<-read.csv("validation/l2_uncertainty_oracle/cases.csv")
jobs<-merge(chosen,data.frame(condition=c("null","planned"))); jobs$job<-seq_len(nrow(jobs))
nsim<-switch(mode,smoke=4L,paired=512L,confirm=2048L)
seed_root<-if(mode=="paired") 810900001L else 912100001L
manifest<-checkpoint_dir(out,c("validation/l2_uncertainty_confirmation.R","R/marginal_scale.R",
  "validation/l2_replication_large/source_script.R"),
  list(mode=mode,variant=variant,jobs=jobs,nsim=nsim,seed_root=seed_root,
    acceptable_coverage=.925,acceptable_false_positive=.075),c("glmmplanr","digest"))
write.csv(jobs,file.path(out,"cases.csv"),row.names=FALSE)
run<-function(k) {
  x<-jobs[k,]; path<-file.path(out,sprintf("case-%02d.rds",k))
  if(file.exists(path)) return(readRDS(path))
  progress<-sub("rds$","progress.rds",path)
  rows<-if(file.exists(progress)) readRDS(progress) else list()
  spec<-make_spec(profiles[profiles$profile==x$profile,])
  if(x$condition=="null") spec<-glmmplanr:::null_target_spec(spec)
  schedule<-generate_design_schedule(spec$dgp$design,x$n_subj,x$n_item)
  cache<-glmmplanr:::prepare_cpp_schedule_cache(spec,schedule)
  bases<-glmmplanr:::allocate_stage_point_streams(seed_root+10000L*x$case,1L,1:2,2L)
  streams<-glmmplanr:::derive_rep_substreams(bases[[if(x$condition=="null") 1 else 2]],nsim)
  if(length(rows)<nsim) for(i in seq.int(length(rows)+1L,nsim)) {
    sim<-glmmplanr:::with_rng_stream(streams[[i]],glmmplanr:::simulate_fast_y_keep(cache))
    timing<-system.time(f<-glmmplanr:::fit_via_cpp_cached(sim,cache,spec,TRUE,spec$analysis_re,FALSE))["elapsed"]
    valid<-glmmplanr:::valid_target_fit(f$est,f$se,f$converged,f$status,f$prob_est,f$prob_se,"pp")
    rows[[i]]<-data.frame(x,variant=variant,replicate=i,truth=spec$true_effect_target,
      estimate=100*f$prob_est,se=100*f$prob_se,est_logit=f$est,se_logit=f$se,
      valid=valid,status=f$status,error=f$error_reason,
      sd_subj_intercept=sqrt(f$Sigma_subj[1,1]),
      sd_subj_slope=if(nrow(f$Sigma_subj)>1) sqrt(f$Sigma_subj[2,2]) else 0,
      sd_item_intercept=sqrt(f$Sigma_item[1,1]),
      sd_item_slope=if(nrow(f$Sigma_item)>1) sqrt(f$Sigma_item[2,2]) else 0,
      elapsed=unname(timing),response_hash=digest::digest(sim$y,algo="xxhash64"))
    if(i%%128L==0L) atomic_rds(rows,progress)
  }
  trials<-do.call(rbind,rows); atomic_rds(trials,path)
  if(file.exists(progress)) unlink(progress)
  cat(sprintf("Confirmation %s %s: case %d %s, %d trials complete\n",variant,mode,x$case,x$condition,nsim))
  trials
}
res<-parallel::mclapply(seq_len(nrow(jobs)),run,mc.cores=workers,mc.set.seed=FALSE,mc.preschedule=FALSE)
if(any(vapply(res,inherits,logical(1),"try-error"))) print(res[vapply(res,inherits,logical(1),"try-error")])
stopifnot(all(vapply(res,is.data.frame,logical(1))))
z<-do.call(rbind,res); stopifnot(nrow(z)==12L*nsim,!anyDuplicated(z[c("case","condition","replicate")]))
summary<-do.call(rbind,lapply(split(z,interaction(z$case,z$condition)),function(v) {
  ok<-v$valid & is.finite(v$estimate) & is.finite(v$se) & v$se>0
  cover<-ok & abs(v$estimate-v$truth)<=qnorm(.975)*v$se
  detect<-ok & abs(v$estimate/v$se)>qnorm(.975)
  joint<-detect & 2*qnorm(.975)*v$se<=10
  ci<-function(k,n) unname(glmmplanr:::wilson_mc_ci(k,n))
  cc<-ci(sum(cover),sum(ok)); dc<-ci(sum(detect),nrow(v)); jc<-ci(sum(joint),nrow(v))
  data.frame(v[1,c("profile","n_subj","n_item","case","condition","variant")],
    nsim=nrow(v),n_valid=sum(ok),bias_pp=mean(v$estimate[ok]-v$truth[ok]),
    empirical_sd_pp=sd(v$estimate[ok]),rms_se_pp=sqrt(mean(v$se[ok]^2)),
    coverage=mean(cover[ok]),coverage_lo=cc[1],coverage_hi=cc[2],
    detection=mean(detect),detection_lo=dc[1],detection_hi=dc[2],
    joint=mean(joint),joint_lo=jc[1],joint_hi=jc[2],
    joint_decision=if(jc[1]>=.8) "pass" else if(jc[2]<.8) "fail" else "indeterminate",
    coverage_floor_decision=if(cc[1]>=.925) "pass" else if(cc[2]<.925) "fail" else "indeterminate",
    null_ceiling_decision=if(v$condition[1]!="null") NA_character_ else
      if(dc[2]<=.075) "pass" else if(dc[1]>.075) "fail" else "indeterminate",
    nominal_adverse_flag=cc[2]<.95 || (v$condition[1]=="null" && dc[1]>.05),
    mean_width_pp=mean(2*qnorm(.975)*v$se[ok]),median_elapsed=median(v$elapsed))
}))
write.csv(z,gzfile(file.path(out,"fits.csv.gz")),row.names=FALSE)
write.csv(summary,file.path(out,"summary.csv"),row.names=FALSE)
file.copy("validation/l2_uncertainty_confirmation.R",file.path(out,"source_script.R"),overwrite=TRUE)
file.copy("src/glmm_laplace.cpp",file.path(out,"source_glmm_laplace.cpp"),overwrite=TRUE)
writeLines(normalizePath("."),file.path(out,"source_directory.txt"))
print(summary,row.names=FALSE)
