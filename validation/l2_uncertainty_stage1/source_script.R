# Check the changed kernel at the existing early-stage limits; no new DGP or grid search.
source("validation/planning_comparison_helpers.R")
for(e in parse("validation/l2_replication_large/source_script.R"))
  if(is.call(e) && identical(e[[1]],as.name("<-")) && is.symbol(e[[2]]) &&
     as.character(e[[2]]) %in% c("profiles","make_spec")) eval(e)
out<-"validation/l2_uncertainty_stage1"
full<-read.csv("validation/l2_uncertainty_confirmation_efficient/fits.csv.gz")
jobs<-read.csv("validation/l2_uncertainty_confirmation_efficient/cases.csv")
checkpoint_dir(out,c("validation/l2_uncertainty_stage1.R",
  "validation/l2_uncertainty_confirmation_efficient/fits.csv.gz",
  "validation/l2_replication_large/source_script.R"),
  list(nsim=64L,seed_root=912100001L,point_tolerance_pp=.01,se_tolerance_pp=.01),
  c("glmmplanr","digest"))
# Use the working R error diagnostics with the numerically identical tested binary.
current<-new.env(parent=asNamespace("glmmplanr"));sys.source("R/engines.R",current)
run<-function(k) {
  path<-file.path(out,sprintf("case-%02d.rds",k));if(file.exists(path))return(readRDS(path))
  x<-jobs[k,];spec<-make_spec(profiles[profiles$profile==x$profile,])
  if(x$condition=="null")spec<-glmmplanr:::null_target_spec(spec)
  sch<-generate_design_schedule(spec$dgp$design,x$n_subj,x$n_item)
  cache<-glmmplanr:::prepare_cpp_schedule_cache(spec,sch)
  b<-glmmplanr:::allocate_stage_point_streams(912100001L+10000L*x$case,1L,1:2,2L)
  streams<-glmmplanr:::derive_rep_substreams(b[[if(x$condition=="null")1 else 2]],64L)
  rows<-lapply(1:64,function(i){
    sim<-glmmplanr:::with_rng_stream(streams[[i]],glmmplanr:::simulate_fast_y_keep(cache))
    f<-current$fit_via_cpp_cached(sim,cache,spec,TRUE,spec$analysis_re,TRUE)
    data.frame(x,replicate=i,valid=f$converged && is.finite(f$prob_est) && is.finite(f$prob_se) && f$prob_se>0,
      estimate=100*f$prob_est,se=100*f$prob_se,status=f$status,error=f$error_reason,
      response_hash=digest::digest(sim$y,algo="xxhash64"))
  })
  z<-do.call(rbind,rows);atomic_rds(z,path);z
}
z<-do.call(rbind,parallel::mclapply(1:12,run,mc.cores=if(.Platform$OS.type=="windows")1L else 3L,mc.set.seed=FALSE))
p<-merge(z,subset(full,replicate<=64L),by=c("case","condition","replicate"),suffixes=c("_early","_final"))
stopifnot(nrow(p)==768L,all(p$response_hash_early==p$response_hash_final))
p$both_valid<-p$valid_early & p$valid_final
p$point_difference_pp<-p$estimate_early-p$estimate_final
p$se_difference_pp<-p$se_early-p$se_final
p$joint_early<-p$valid_early & abs(p$estimate_early/p$se_early)>qnorm(.975) & 2*qnorm(.975)*p$se_early<=10
p$joint_final<-p$valid_final & abs(p$estimate_final/p$se_final)>qnorm(.975) & 2*qnorm(.975)*p$se_final<=10
write.csv(p,gzfile(file.path(out,"fits.csv.gz")),row.names=FALSE)
s<-do.call(rbind,lapply(split(p,interaction(p$case,p$condition)),function(v)data.frame(
  case=v$case[1],condition=v$condition[1],n=nrow(v),n_early_valid=sum(v$valid_early),n_final_valid=sum(v$valid_final),
  both_valid=sum(v$both_valid),max_point_difference_pp=max(abs(v$point_difference_pp[v$both_valid])),
  max_se_difference_pp=max(abs(v$se_difference_pp[v$both_valid])),
  differences_over_tolerance=sum(v$both_valid & (abs(v$point_difference_pp)>.01|abs(v$se_difference_pp)>.01)),
  joint_disagreements=sum(v$joint_early!=v$joint_final))))
write.csv(s,file.path(out,"summary.csv"),row.names=FALSE)
file.copy("validation/l2_uncertainty_stage1.R",file.path(out,"source_script.R"),overwrite=TRUE)
print(s,row.names=FALSE)
