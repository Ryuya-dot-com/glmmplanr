# Inspect all tolerance exceedances; do not raise the prespecified tolerance.
library(glmmplanr)
args <- commandArgs(TRUE)
out <- if(length(args)) args[1] else "validation/covariance_mode_replay"
old <- "validation/covariance_review_medium"
x <- read.csv(file.path(out,"comparison.csv"))
r <- read.csv(file.path(out,"reference_summary.csv"))
stopifnot(nrow(x)==2048L,!anyDuplicated(x[c("case","replicate")]),all(x$valid),all(x$early_valid),
  !any(x$coverage_changed|x$detection_changed|x$width_changed|x$joint_changed))
for(e in parse(file.path(old,"source_script.R"))) {
  if(is.call(e) && identical(e[[1]],as.name("<-")) && identical(e[[2]],as.name("dense_check"))) eval(e)
}
exceed <- subset(x,abs(estimate_change_pp)>1e-4 | abs(se_change_pp)>1e-4)
rows <- list()
for(k in seq_len(nrow(exceed))) {
  q<-exceed[k,]; block<-(q$replicate-1L)%/%16L+1L
  a<-readRDS(file.path(old,sprintf("case-%02d-block-%02d.rds",q$case,block)))
  z<-a$details[[which(vapply(a$details,function(v) v$rows$replicate[1],0)==q$replicate)]]
  cache<-glmmplanr:::prepare_cpp_schedule_cache(a$spec,a$schedule)
  ws<-glmmplanr:::get_cached_cpp_ws(cache,a$spec$analysis_re)
  f<-glmmplanr:::cpp_fit_crossed_glmm_ws(ws,as.integer(z$y),chance=a$case$chance,
    corr_mode=as.integer(a$case$correlation!=0),regularized=TRUE,return_ranef=TRUE,
    prior_beta_sd=a$spec$prior_beta_sd,prior_tau_sd=a$spec$prior_tau_sd,lkj_eta=a$spec$lkj_eta,
    max_em_iter=if(a$case$chance>0)400L else 160L,max_pirls_iter=if(a$case$chance>0)44L else 22L)
  chk<-dense_check(f,cache,as.integer(z$y),a$case$chance)
  stopifnot(f$converged,abs(chk['ll_difference'])<1e-8,chk['mode_score']<2e-5)
  rows[[k]]<-cbind(q,data.frame(before_iterations=z$native$iterations,after_iterations=f$iterations,
    after_covariance_score=f$covariance_score_norm,after_inner_step=f$inner_step,
    dense_ll_difference=chk[1],dense_mode_score=chk[2]))
}
if(length(rows)) {
  inspected<-do.call(rbind,rows)
  write.csv(inspected,file.path(out,"tolerance_exceedances.csv"),row.names=FALSE)
  print(inspected[c("case","replicate","estimate_change_pp","se_change_pp",
    "before_iterations","after_iterations","after_covariance_score","dense_ll_difference")],row.names=FALSE)
}
paths<-list.files(out,pattern="^case-.*-block-.*rds$",full.names=TRUE)
first<-list()
for(p in paths) for(z in readRDS(p)$references) if(z$rows$replicate[1]==1L) first[[length(first)+1L]]<-z
dense<-do.call(rbind,lapply(first,function(z)do.call(rbind,z$dense)))
stopifnot(length(first)==16L,nrow(dense)==32L,all(is.finite(dense)),
  max(abs(dense[,"ll_difference"]))<1e-6,max(abs(dense[,"mode_score"]))<2e-5,
  all(r$profile_valid),all(r$failed_evaluations==0),max(r$profile_score)<.01,
  max(abs(r$fine_se_difference_pp),na.rm=TRUE)<.05,
  max(abs(r$alternate_objective_difference),na.rm=TRUE)<1e-4)
cat('Three types of evidence are separate: 2,048 native replays; tolerance exceptions inspected;',
  nrow(r),'selected profile replays.\n')
cat('Max first-case dense likelihood difference:',max(abs(dense[,"ll_difference"])),
  '; mode score:',max(abs(dense[,"mode_score"])), '\n')
cat('Nominal inference was not independently revalidated. The original comparison tolerance was not met in',
  nrow(exceed),'cases; no saved decision changed.\n')
