# Replay saved responses after the mode-roundoff/diagnostic-contract repair.
source("validation/planning_comparison_helpers.R")
old <- "validation/covariance_review_medium"
out <- "validation/covariance_mode_replay"
cases <- read.csv(file.path(old,"cases.csv"))
diagnostics <- read.csv(file.path(old,"diagnostics.csv"))
scores <- read.csv(file.path(old,"score_checks.csv"))
selected <- unique(rbind(diagnostics[!diagnostics$profile_valid,c("case","replicate")],
  scores[!scores$score_evaluable,c("case","replicate")],
  data.frame(case=cases$case,replicate=1L)))
paths <- list.files(old,pattern="^case-.*-block-.*rds$",full.names=TRUE)
checkpoint_dir(out,c("validation/covariance_mode_replay.R",paths,
  file.path(old,c("source_script.R","manifest.rds","score_checks.csv")),"R/plan_methods.R"),
  list(nsim=2048L,selected=selected,point_tolerance_pp=1e-4,se_tolerance_pp=1e-4,
    reference="unchanged archived optimizer, priors and acceptance rules"),c("glmmplanr","digest"))
for(e in parse(file.path(old,"source_script.R"))) {
  if(is.call(e) && identical(e[[1]],as.name("<-")) && is.symbol(e[[2]]) &&
     as.character(e[[2]]) %in% c("make_case","unpack","dense_check","run_trial")) eval(e)
}
atomic_rds(list(selected=selected,cases=cases),file.path(out,"cases.rds"))
run_block <- function(path) {
  target <- file.path(out,basename(path)); if(file.exists(target)) return(target)
  a <- readRDS(path); obj <- list(spec=a$spec,schedule=a$schedule,
    cache=glmmplanr:::prepare_cpp_schedule_cache(a$spec,a$schedule))
  x <- a$case; rows <- list(); references <- list()
  bases <- glmmplanr:::allocate_stage_point_streams(952260001L,1L,cases$case,nrow(cases))
  streams <- glmmplanr:::derive_rep_substreams(bases[[x$case]],128L)
  for(z in a$details) {
    before <- z$rows[z$rows$method=="native",]; i <- before$replicate
    sim <- list(y=as.integer(z$y))
    f <- glmmplanr:::fit_via_cpp_cached(sim,obj$cache,a$spec,TRUE,a$spec$analysis_re,FALSE)
    early <- glmmplanr:::fit_via_cpp_cached(sim,obj$cache,a$spec,TRUE,a$spec$analysis_re,TRUE)
    est<-100*f$prob_est; se<-100*f$prob_se; valid<-isTRUE(f$converged) && is.finite(se) && se>0
    covered<-valid && abs(est-before$truth)<=qnorm(.975)*se
    detect<-valid && abs(est)>qnorm(.975)*se
    width<-valid && 2*qnorm(.975)*se<=20
    rows[[length(rows)+1L]] <- data.frame(case=x$case,replicate=i,valid=valid,
      early_valid=isTRUE(early$converged),before_valid=before$valid,
      estimate=est,se=se,estimate_change_pp=est-before$estimate,se_change_pp=se-before$se,
      covariance_max_change=max(abs(f$Sigma_subj-z$native$Sigma_subj),
        abs(f$Sigma_item-z$native$Sigma_item)),covered=covered,detected=detect,width_ok=width,
      joint=detect && width,coverage_changed=covered!=before$covered,
      detection_changed=detect!=before$detected,width_changed=width!=before$width_ok,
      joint_changed=(detect && width)!=before$joint,response_hash=before$response_hash)
    if(any(selected$case==x$case & selected$replicate==i)) {
      ref <- run_trial(x,i,obj,streams[[i]])
      stopifnot(identical(ref$y,z$y),all(ref$rows$response_hash==before$response_hash))
      ref$before_profile_valid <- z$diagnostic$profile_valid
      ref$before_score_evaluable <- scores$score_evaluable[scores$case==x$case & scores$replicate==i]
      references[[length(references)+1L]] <- ref
    }
  }
  atomic_rds(list(comparison=do.call(rbind,rows),references=references),target)
  cat(basename(path),"saved\n"); target
}
results <- parallel::mclapply(paths,run_block,mc.cores=if(.Platform$OS.type=="windows") 1L else 6L,
  mc.set.seed=FALSE,mc.preschedule=FALSE)
stopifnot(all(vapply(results,function(p) is.character(p) && length(p)==1L && file.exists(p),logical(1))))
batches<-lapply(results,readRDS)
comparison<-do.call(rbind,lapply(batches,`[[`,"comparison"))
ref<-unlist(lapply(batches,`[[`,"references"),recursive=FALSE)
ref_summary<-do.call(rbind,lapply(ref,function(z) cbind(z$diagnostic,
  before_profile_valid=z$before_profile_valid,before_score_evaluable=z$before_score_evaluable)))
stopifnot(nrow(comparison)==2048L,nrow(ref_summary)==nrow(selected),
  !anyDuplicated(comparison[c("case","replicate")]),!anyDuplicated(ref_summary[c("case","replicate")]))
write.csv(comparison,file.path(out,"comparison.csv"),row.names=FALSE)
write.csv(ref_summary,file.path(out,"reference_summary.csv"),row.names=FALSE)
write.csv(do.call(rbind,lapply(ref,`[[`,"rows")),file.path(out,"reference_fits.csv"),row.names=FALSE)
file.copy("validation/covariance_mode_replay.R",file.path(out,"source_script.R"),overwrite=TRUE)
file.copy("src/glmm_laplace.cpp",file.path(out,"source_glmm_laplace.cpp"),overwrite=TRUE)
cat("Native valid:",sum(comparison$valid),"/2048; early valid:",sum(comparison$early_valid),"\n")
print(sapply(comparison[c("estimate_change_pp","se_change_pp","covariance_max_change")],function(z) max(abs(z))))
print(colSums(comparison[c("coverage_changed","detection_changed","width_changed","joint_changed")]))
print(with(ref_summary,table(before_profile_valid,profile_valid)))
