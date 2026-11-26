# Recompute summaries from retained trials; no model refitting.
args <- commandArgs(TRUE)
out <- if(length(args)) args[1] else "validation/covariance_review_medium"
library(glmmplanr)
x <- read.csv(file.path(out,"fits.csv.gz"))
d <- read.csv(file.path(out,"diagnostics.csv"))
if(file.exists(file.path(out,"score_checks.csv"))) {
  checked <- read.csv(file.path(out,"score_checks.csv"))
  d <- merge(d,checked,by=c("case","replicate"))
  d$native_logsd_score <- d$scaled_logsd_score
}
manifest <- readRDS(file.path(out,"manifest.rds"))
nsim <- manifest$settings$nsim
stopifnot(nrow(x)==16L*nsim*2L,nrow(d)==16L*nsim,
  !anyDuplicated(x[c("case","replicate","method")]),
  !anyDuplicated(d[c("case","replicate")]),!anyNA(x$valid))
z <- qnorm(.975)
stopifnot(identical(x$covered,x$valid & abs(x$estimate-x$truth)<=z*x$se),
  identical(x$detected,x$valid & abs(x$estimate)>z*x$se),
  identical(x$width_ok,x$valid & 2*z*x$se<=20),
  identical(x$joint,x$detected & x$width_ok))
ci <- function(k,n) unname(glmmplanr:::wilson_mc_ci(k,n))
rate <- function(y,denom=length(y)) {
  k <- sum(y); c(n=denom,k=k,rate=if(denom>0) k/denom else NA_real_,ci(k,denom))
}
decision <- function(k,n) {
  v <- ci(k,n)
  if(v[1]>=.8) "pass" else if(v[2]<.8) "fail" else "indeterminate"
}
summary <- do.call(rbind,lapply(split(x,interaction(x$case,x$method)),function(v) {
  ok <- v$valid; a <- v[ok,]; err <- a$estimate-a$truth
  r <- data.frame(v[1,c("case","layout","n_subj","n_item","per_person","chance",
    "correlation","effect_pp","method")],nsim=nrow(v),n_valid=sum(ok),
    bias_pp=mean(err),bias_mcse_pp=sd(err)/sqrt(length(err)),rmse_pp=sqrt(mean(err^2)),
    empirical_sd_pp=sd(a$estimate),rms_se_pp=sqrt(mean(a$se^2)),
    mean_width_pp=mean(2*z*a$se),logit_bias=mean(a$logit-a$logit_truth))
  for(nm in c("valid","covered","detected","width_ok","joint","logit_covered")) {
    rr <- rate(v[[nm]])
    r[[paste0(nm,"_count")]] <- rr[2]
    r[[paste0(nm,"_rate")]] <- rr[3]
    r[[paste0(nm,"_lo")]] <- rr[4]
    r[[paste0(nm,"_hi")]] <- rr[5]
  }
  cc <- rate(a$covered)
  r$coverage_valid <- cc[3]; r$coverage_valid_lo <- cc[4]; r$coverage_valid_hi <- cc[5]
  r$joint_decision <- decision(sum(v$joint),nrow(v))
  r
}))
a <- subset(x,method=="native"); b <- subset(x,method=="profile")
p <- merge(a,b,by=c("case","replicate"),suffixes=c("_native","_profile"))
stopifnot(nrow(p)==16L*nsim,all(p$response_hash_native==p$response_hash_profile))
paired <- do.call(rbind,lapply(split(p,p$case),function(v) {
  ok <- v$valid_native & v$valid_profile; a <- v[ok,]
  diff <- a$estimate_profile-a$estimate_native
  r <- data.frame(case=v$case[1],nsim=nrow(v),n_common=sum(ok),
    mean_estimate_change_pp=mean(diff),change_mcse_pp=sd(diff)/sqrt(length(diff)),
    rms_estimate_change_pp=sqrt(mean(diff^2)),max_abs_estimate_change_pp=max(abs(diff),0),
    rms_se_ratio=sqrt(mean(a$se_profile^2)/mean(a$se_native^2)))
  for(nm in c("valid","covered","detected","width_ok","joint")) {
    dd <- as.numeric(v[[paste0(nm,"_profile")]])-as.numeric(v[[paste0(nm,"_native")]])
    r[[paste0(nm,"_change")]] <- mean(dd)
    r[[paste0(nm,"_mcse")]] <- sd(dd)/sqrt(length(dd))
    r[[paste0(nm,"_changed")]] <- sum(dd!=0)
  }
  r$coverage_common_change <- mean(a$covered_profile-a$covered_native)
  r$coverage_common_mcse <- sd(a$covered_profile-a$covered_native)/sqrt(nrow(a))
  r
}))
diagnostic <- do.call(rbind,lapply(split(d,d$case),function(v) {
  finite <- is.finite(v$native_logsd_score)
  data.frame(case=v$case[1],nsim=nrow(v),native_valid=sum(v$native_valid),
    early_valid=sum(v$early_valid),early_failed_final_valid=sum(!v$early_valid & v$native_valid),
    profile_valid=sum(v$profile_valid),score_available=sum(finite),
    native_score_median=median(v$native_logsd_score[finite]),
    native_score_max=max(v$native_logsd_score[finite]),
    native_score_above_01=sum(v$native_logsd_score[finite]>.01),
    profile_score_max=max(v$profile_score[v$profile_valid]),
    improvement_median=median(v$objective_improvement,na.rm=TRUE),
    boundary_count=sum(v$boundary,na.rm=TRUE),optimizer_failures=sum(v$optimizer_code!=0,na.rm=TRUE),
    failed_evaluations=sum(v$failed_evaluations),elapsed_sum=sum(v$elapsed),
    fine_se_difference_pp=v$fine_se_difference_pp[1],
    alternate_objective_difference=v$alternate_objective_difference[1])
}))
for (nm in c("summary","paired","diagnostic"))
  write.csv(get(nm),file.path(out,paste0(nm,"_summary.csv")),row.names=FALSE)
parameters <- list()
for(path in list.files(out,pattern="^case-.*-block-.*rds$",full.names=TRUE)) {
  batch <- readRDS(path)
  truth <- c(.8,.35,.6,.4,batch$case$correlation,batch$case$correlation)
  for (trial in batch$details) for (method in c("native","profile")) {
    fit <- if(method=="native") trial$native else trial$profile$fit
    if(is.null(fit)) next
    val <- c(sqrt(c(diag(fit$Sigma_subj),diag(fit$Sigma_item))),
      cov2cor(fit$Sigma_subj)[1,2],cov2cor(fit$Sigma_item)[1,2])
    parameters[[length(parameters)+1L]] <- data.frame(case=batch$case$case,
      replicate=trial$rows$replicate[1],method=method,
      valid=trial$rows$valid[trial$rows$method==method],
      parameter=c("participant_intercept_sd","participant_slope_sd","item_intercept_sd",
        "item_slope_sd","participant_correlation","item_correlation"),estimate=val,truth=truth)
  }
}
parameters <- do.call(rbind,parameters)
parameter_summary <- do.call(rbind,lapply(split(parameters,
  interaction(parameters$case,parameters$method,parameters$parameter)),function(v) {
    a <- v[v$valid,]; error <- a$estimate-a$truth
    data.frame(v[1,c("case","method","parameter","truth")],n_planned=nsim,n_valid=nrow(a),
      mean=mean(a$estimate),bias=mean(error),bias_mcse=sd(error)/sqrt(length(error)),
      rmse=sqrt(mean(error^2)))
  }))
write.csv(parameter_summary,file.path(out,"parameter_summary.csv"),row.names=FALSE)
print(summary[c("case","method","n_valid","bias_pp","rms_se_pp","coverage_valid",
  "detected_rate","joint_rate","joint_decision")],row.names=FALSE)
cat("Integrity checks passed: keys, counts, response pairs, and all saved decisions.\n")

if (identical(manifest$settings$mode,"medium")) {
  library(ggplot2)
  s <- summary
  s$method <- factor(s$method,c("native","profile"),c("Native joint delta","Profile diagnostic"))
  s$guessing <- factor(s$chance,c(0,1/3),c("0","1/3"))
  s$effect <- factor(s$effect_pp,c(0,10),c("Null: 0 pp","Planned: 10 pp"))
  s$layout <- factor(s$layout,c("complete","partial"),
    c("32 x 24, complete","72 x 60, 30/person"))
  s$rho <- factor(s$correlation,c(0,.2),c("Correlation 0","Correlation 0.2"))
  pos <- position_dodge(width=.3)
  fig <- ggplot(s,aes(guessing,coverage_valid,colour=method,shape=method)) +
    geom_hline(yintercept=.95,linetype=2,colour="grey45") +
    geom_errorbar(aes(ymin=coverage_valid_lo,ymax=coverage_valid_hi),width=.12,position=pos) +
    geom_point(position=pos,size=2) + facet_grid(effect ~ layout + rho) +
    scale_colour_manual(values=c("#0072B2","#D55E00")) +
    scale_y_continuous(labels=function(v) paste0(round(100*v),"%")) +
    labs(x="Guessing floor",y="Coverage among valid 95% intervals",colour=NULL,shape=NULL) +
    theme_bw() + theme(legend.position="bottom")
  ggsave(file.path(out,"coverage.png"),fig,width=11,height=5.5,dpi=160)
  ggsave(file.path(out,"coverage.pdf"),fig,width=11,height=5.5)
}
