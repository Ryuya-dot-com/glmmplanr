# Recalculate trial-level decisions; keep exploratory and independent evidence separate.
source("validation/planning_comparison_helpers.R")
include_bootstrap<-!("--independent-only" %in% commandArgs(TRUE))
out<-"validation/l2_uncertainty_comparison"
dir.create(out,showWarnings=FALSE)
old<-read.csv("validation/l2_uncertainty_confirmation_baseline/fits.csv.gz")
new<-read.csv("validation/l2_uncertainty_confirmation_guarded/fits.csv.gz")
keys<-c("case","condition","replicate")
stopifnot(nrow(old)==24576L,nrow(new)==24576L,!anyDuplicated(old[keys]),!anyDuplicated(new[keys]))
both<-merge(old,new,by=keys,suffixes=c("_old","_new"))
stopifnot(nrow(both)==24576L,all(both$response_hash_old==both$response_hash_new),
  max(abs(both$truth_old-both$truth_new))<1e-12)
outcomes<-function(z) {
  z$valid<-z$valid & is.finite(z$estimate) & is.finite(z$se) & z$se>0
  z$covered<-z$valid & abs(z$estimate-z$truth)<=qnorm(.975)*z$se
  z$detected<-z$valid & abs(z$estimate/z$se)>qnorm(.975)
  z$joint<-z$detected & 2*qnorm(.975)*z$se<=10
  z
}
old<-outcomes(old);new<-outcomes(new)
intermediate<-outcomes(read.csv("validation/l2_uncertainty_confirmation_efficient/fits.csv.gz"))
replay<-merge(intermediate,new,by=keys,suffixes=c("_before","_guarded"))
stopifnot(nrow(replay)==24576L,all(replay$response_hash_before==replay$response_hash_guarded))
# Declared before inspection: 1e-8 pp for point/SE equivalence; audit every decision separately.
replay_audit<-do.call(rbind,lapply(split(replay,interaction(replay$case,replay$condition)),function(v){
  ok<-v$valid_before & v$valid_guarded
  point<-abs(v$estimate_guarded[ok]-v$estimate_before[ok]);se<-abs(v$se_guarded[ok]-v$se_before[ok])
  data.frame(case=v$case[1],condition=v$condition[1],n=nrow(v),n_both_valid=sum(ok),
    n_validity_changed=sum(v$valid_before!=v$valid_guarded),
    n_status_changed=sum(v$status_before!=v$status_guarded),
    max_point_difference_pp=max(point),max_se_difference_pp=max(se),
    n_beyond_tolerance=sum(point>1e-8|se>1e-8),
    n_coverage_changed=sum(v$covered_before!=v$covered_guarded),
    n_detection_changed=sum(v$detected_before!=v$detected_guarded),
    n_joint_changed=sum(v$joint_before!=v$joint_guarded))
}))
write.csv(replay_audit,file.path(out,"beta_guard_replay.csv"),row.names=FALSE)
for(pair in list(list(old,"baseline"),list(new,"guarded"))) {
  z<-pair[[1]]
  s<-read.csv(paste0("validation/l2_uncertainty_confirmation_",pair[[2]],"/summary.csv"))
  for(k in seq_len(nrow(s))) {
    v<-subset(z,case==s$case[k]&condition==s$condition[k])
    stopifnot(nrow(v)==2048L,sum(v$valid)==s$n_valid[k],
      abs(mean(v$detected)-s$detection[k])<1e-12,
      abs(mean(v$covered[v$valid])-s$coverage[k])<1e-12,
      abs(mean(v$joint)-s$joint[k])<1e-12)
    ci<-glmmplanr:::wilson_mc_ci(sum(v$joint),nrow(v))
    decision<-if(ci[1]>=.8) "pass" else if(ci[2]<.8) "fail" else "indeterminate"
    stopifnot(decision==s$joint_decision[k])
  }
}
paired<-merge(old,new,by=keys,suffixes=c("_old","_new"))
differences<-do.call(rbind,lapply(split(paired,interaction(paired$case,paired$condition)),function(v) {
  ok<-v$valid_old & v$valid_new
  dcover<-as.numeric(v$covered_new)-as.numeric(v$covered_old)
  ddetect<-as.numeric(v$detected_new)-as.numeric(v$detected_old)
  djoint<-as.numeric(v$joint_new)-as.numeric(v$joint_old)
  data.frame(case=v$case[1],condition=v$condition[1],n=nrow(v),n_both_valid=sum(ok),
    mean_point_change_pp=mean(v$estimate_new[ok]-v$estimate_old[ok]),
    se_ratio_geometric=exp(mean(log(v$se_new[ok]/v$se_old[ok]))),
    coverage_change_all_attempted=mean(dcover),
    coverage_change_all_attempted_mcse=sd(dcover)/sqrt(nrow(v)),
    coverage_change_both_valid=mean(dcover[ok]),
    coverage_change_both_valid_mcse=sd(dcover[ok])/sqrt(sum(ok)),
    detection_change=mean(ddetect),detection_change_mcse=sd(ddetect)/sqrt(nrow(v)),
    joint_change=mean(djoint),joint_change_mcse=sd(djoint)/sqrt(nrow(v)))
}))
write.csv(differences,file.path(out,"paired_differences.csv"),row.names=FALSE)
summary<-rbind(read.csv("validation/l2_uncertainty_confirmation_baseline/summary.csv"),
  read.csv("validation/l2_uncertainty_confirmation_guarded/summary.csv"))
all_trials<-rbind(old,new)
summary$valid_rate<-summary$n_valid/summary$nsim
for(k in seq_len(nrow(summary))) {
  v<-subset(all_trials,case==summary$case[k]&condition==summary$condition[k]&variant==summary$variant[k])
  vc<-glmmplanr:::wilson_mc_ci(sum(v$valid),nrow(v))
  dc<-glmmplanr:::wilson_mc_ci(sum(v$detected),sum(v$valid))
  summary$valid_lo[k]<-vc[1];summary$valid_hi[k]<-vc[2]
  summary$detection_given_valid[k]<-sum(v$detected)/sum(v$valid)
  summary$detection_given_valid_lo[k]<-dc[1];summary$detection_given_valid_hi[k]<-dc[2]
}
write.csv(summary,file.path(out,"independent_summary.csv"),row.names=FALSE)
write.csv(subset(rbind(old,new),!valid),file.path(out,"failures.csv"),row.names=FALSE)
bench0<-read.csv("validation/l2_uncertainty_benchmark_baseline/fits.csv")
bench1<-read.csv("validation/l2_uncertainty_benchmark_guarded/fits.csv")
b<-merge(bench0,bench1,by=c("case","replicate"),suffixes=c("_old","_new"))
stopifnot(nrow(b)==128L,all(b$hash_old==b$hash_new))
untouched<-subset(b,kind_old!="uncorrelated")
stopifnot(all(untouched$valid_old==untouched$valid_new),
  max(abs(untouched$estimate_old-untouched$estimate_new),na.rm=TRUE)<1e-10,
  max(abs(untouched$se_old-untouched$se_new),na.rm=TRUE)<1e-10)
timing<-do.call(rbind,lapply(split(b,b$case),function(v)
  data.frame(case=v$case[1],profile=v$profile_old[1],n_subj=v$n_subj_old[1],n_item=v$n_item_old[1],
    n=nrow(v),old_valid=sum(v$valid_old),new_valid=sum(v$valid_new),
    old_median_cpu=median(v$cpu_old),new_median_cpu=median(v$cpu_new),
    aggregate_cpu_ratio=sum(v$cpu_new)/sum(v$cpu_old))))
write.csv(timing,file.path(out,"cpu_comparison.csv"),row.names=FALSE)

profile<-read.csv("validation/l2_uncertainty_profile_optimized/summary.csv")
candidate<-read.csv("validation/l2_uncertainty_efficient_profile/summary.csv")
numerical<-merge(profile,candidate,by=keys,suffixes=c("_reference","_candidate"))
stopifnot(nrow(numerical)==96L,
  all(numerical$response_hash_reference==numerical$response_hash_candidate),
  all(numerical$optimization_code_reference==0L))
numerical$point_difference_pp<-numerical$current_estimate_candidate-numerical$profile_estimate_reference
stopifnot(max(abs(numerical$point_difference_pp))<.001,
  sum(is.finite(numerical$likelihood_difference_reference))==12L,
  max(abs(numerical$likelihood_difference_reference),na.rm=TRUE)<1e-8,
  max(abs(numerical$fine_se_difference_reference),na.rm=TRUE)<1e-5)
write.csv(numerical,file.path(out,"independent_profile_comparison.csv"),row.names=FALSE)
profile_files<-list.files("validation/l2_uncertainty_efficient_profile",
  pattern="^case-[0-9]+\\.rds$",full.names=TRUE)
stopifnot(length(profile_files)==12L,
  all(vapply(profile_files,function(path)
    all(vapply(readRDS(path)$details,function(x)isTRUE(x$current$converged),logical(1))),logical(1))))

oracle<-read.csv("validation/l2_uncertainty_oracle/summary.csv")
write.csv(oracle,file.path(out,"exploratory_oracle.csv"),row.names=FALSE)
if(include_bootstrap) {
bootstrap<-read.csv("validation/l2_uncertainty_bootstrap/summary.csv")
write.csv(bootstrap,file.path(out,"exploratory_bootstrap.csv"),row.names=FALSE)
boot_trials<-read.csv("validation/l2_uncertainty_bootstrap/fits.csv.gz")
oracle_trials<-read.csv("validation/l2_uncertainty_oracle/fits.csv.gz")
stopifnot(nrow(boot_trials)==4608L,
  !anyDuplicated(boot_trials[c(keys,"method")]))
replayed<-merge(subset(boot_trials,method!="bootstrap"),oracle_trials,
  by=c(keys,"method"),suffixes=c("_bootstrap_run","_oracle_run"))
stopifnot(nrow(replayed)==3072L,
  all(replayed$response_hash_bootstrap_run==replayed$response_hash_oracle_run),
  all(replayed$valid_bootstrap_run==replayed$valid_oracle_run),
  max(abs(replayed$estimate_bootstrap_run-replayed$estimate_oracle_run),na.rm=TRUE)<1e-7,
  max(abs(replayed$se_bootstrap_run-replayed$se_oracle_run),na.rm=TRUE)<1e-7)
}

library(ggplot2)
labels<-c(published_8pp="8pp, intercepts only",slopes_8pp="8pp, moderate slopes",
  conservative_5pp="5pp, moderate slopes",high_variance_8pp="8pp, larger variances")
summary$design<-paste0(labels[summary$profile],"\nN=",summary$n_subj,", J=",summary$n_item)
summary$design<-factor(summary$design,
  levels=rev(summary$design[match(c(31L,50L,11L,21L,25L,30L),summary$case)]))
summary$method<-factor(summary$variant,levels=c("baseline","guarded"),labels=c("Previous","Corrected"))
cov<-subset(summary,condition=="planned")
cov$metric<-"Coverage of planned 95% interval"; cov$rate<-cov$coverage; cov$lo<-cov$coverage_lo;cov$hi<-cov$coverage_hi
fp<-subset(summary,condition=="null")
fp$metric<-"Null rejection at 5%";fp$rate<-fp$detection;fp$lo<-fp$detection_lo;fp$hi<-fp$detection_hi
plotdata<-rbind(cov,fp)
reference<-data.frame(metric=c("Coverage of planned 95% interval","Null rejection at 5%"),reference=c(.95,.05))
p<-ggplot(plotdata,aes(rate,design,color=method))+geom_vline(data=reference,aes(xintercept=reference),
  color="grey50",linetype="dashed")+
  geom_errorbar(aes(xmin=lo,xmax=hi),width=.15,position=position_dodge(width=.5),orientation="y")+
  geom_point(position=position_dodge(width=.5),size=2)+facet_wrap(~metric,scales="free_x")+
  scale_x_continuous(labels=function(x)paste0(round(100*x),"%"))+
  scale_color_manual(values=c("Previous"="#6D7784","Corrected"="#B74627"))+
  labs(x=NULL,y=NULL,color=NULL,title="Independent confirmation of native interval calibration",
    subtitle="Same responses, 2,048 replications per design and condition; every person answers 30 words",
    caption="Pointwise 95% Wilson intervals: coverage uses valid fits; null rejection uses all attempted fits.\nDashed lines are nominal values, not equivalence limits. Selected designs, not a repeat of the full grid.")+
  theme_minimal(base_size=11)+theme(legend.position="bottom",panel.grid.minor=element_blank(),
    plot.caption=element_text(hjust=0),strip.text=element_text(face="bold"))
ggsave("validation/l2_uncertainty_confirmation.png",p,width=11,height=7,dpi=160)
print(summary[c("case","condition","variant","n_valid","coverage","detection","joint_decision",
  "coverage_floor_decision","null_ceiling_decision","nominal_adverse_flag")],row.names=FALSE)
print(timing,row.names=FALSE)
print(replay_audit,row.names=FALSE)
cat("All 24,576 response pairs, decisions and unchanged-path values audited.\n")
