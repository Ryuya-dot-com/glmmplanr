# Re-evaluate ONLY initial score perturbations, retaining failures as missing scores.
# The original optimizer's finite penalty is not a likelihood derivative.
library(glmmplanr)
args <- commandArgs(TRUE)
out <- if(length(args)) args[1] else "validation/covariance_review_medium"
paths <- list.files(out,pattern="^case-.*-block-.*rds$",full.names=TRUE)
run <- function(path) {
  batch <- readRDS(path); spec <- batch$spec; x <- batch$case
  cache <- glmmplanr:::prepare_cpp_schedule_cache(spec,batch$schedule)
  ws <- glmmplanr:::get_cached_cpp_ws(cache,spec$analysis_re)
  do.call(rbind,lapply(batch$details,function(z) {
    n <- z$native; rho <- c(cov2cor(n$Sigma_subj)[1,2],cov2cor(n$Sigma_item)[1,2])
    theta <- log(sqrt(c(diag(n$Sigma_subj),diag(n$Sigma_item))))
    if(x$correlation!=0) theta <- c(theta,atanh(rho))
    evaluate <- function(t) {
      rr <- if(length(t)==6) tanh(t[5:6]) else c(0,0)
      cv <- function(sd,r) outer(sd,sd)*matrix(c(1,r,r,1),2)
      f <- tryCatch(glmmplanr:::cpp_fit_crossed_glmm_ws(ws,as.integer(z$y),chance=x$chance,
        corr_mode=if(x$correlation==0) 0L else 1L,regularized=TRUE,fix_cov=TRUE,
        compute_joint_vcov=FALSE,init_beta=n$beta,
        init_Sigma_subj=cv(exp(t[1:2]),rr[1]),init_Sigma_item=cv(exp(t[3:4]),rr[2]),
        prior_beta_sd=spec$prior_beta_sd,prior_tau_sd=spec$prior_tau_sd,lkj_eta=spec$lkj_eta,
        max_em_iter=if(x$chance>0) 400L else 160L,
        max_pirls_iter=if(x$chance>0) 44L else 22L,tol=1.5e-4),error=function(e) NULL)
      if(is.null(f) || !isTRUE(f$converged)) return(list(value=NA_real_,
        status=if(is.null(f)) "error" else f$status,step=if(is.null(f)) NA_real_ else f$inner_step))
      prior <- sum(dnorm(f$beta,0,c(1.5,1)*spec$prior_beta_sd,log=TRUE)) +
        sum(dt(exp(t[1:4])/spec$prior_tau_sd,3,log=TRUE)) + sum(t[1:4]) +
        (spec$lkj_eta-1)*sum(log1p(-rr^2))
      list(value=-f$logLik-prior,status=f$status,step=f$inner_step)
    }
    evaluations <- list(); gradient <- numeric(length(theta))
    for(j in seq_along(theta)) {
      step <- replace(numeric(length(theta)),j,.001)
      up<-evaluate(theta+step);dn<-evaluate(theta-step)
      gradient[j]<-(up$value-dn$value)/.002
      evaluations<-c(evaluations,list(up,dn))
    }
    good <- all(is.finite(gradient))
    failures<-Filter(function(v) !is.finite(v$value),evaluations)
    data.frame(case=x$case,replicate=z$rows$replicate[1],score_evaluable=good,
      scaled_logsd_score=if(good) max(abs(gradient[1:4])*
        sqrt(c(n$var_log_tau_subj,n$var_log_tau_item))) else NA_real_,
      rho_score=if(good && length(theta)==6) max(abs(gradient[5:6])) else NA_real_,
      failed_perturbations=length(failures),
      failure_status=paste(unique(vapply(failures,`[[`,"","status")),collapse=";"),
      max_failed_inner_step=if(length(failures)) max(vapply(failures,`[[`,0,"step"),na.rm=TRUE) else NA_real_)
  }))
}
rows<-parallel::mclapply(paths,run,mc.cores=if(.Platform$OS.type=="windows") 1L else 6L,
  mc.preschedule=FALSE,mc.set.seed=FALSE)
stopifnot(all(vapply(rows,is.data.frame,logical(1))))
scores<-do.call(rbind,rows)
old<-read.csv(file.path(out,"diagnostics.csv"))
p<-merge(old,scores,by=c("case","replicate"))
stopifnot(nrow(p)==nrow(old),!anyDuplicated(scores[c("case","replicate")]),
  all(is.na(scores$scaled_logsd_score)==!scores$score_evaluable),
  max(abs(p$native_logsd_score[p$score_evaluable]-p$scaled_logsd_score[p$score_evaluable])) < 1e-6)
write.csv(scores,file.path(out,"score_checks.csv"),row.names=FALSE)
print(aggregate(cbind(score_evaluable,failed_perturbations)~case,scores,sum),row.names=FALSE)
cat('Valid scores reproduce within 1e-6; unevaluable scores retained as NA.\n')
