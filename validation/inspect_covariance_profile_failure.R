# Follow-up of a declared failure; does not replace any main-run fit.
library(glmmplanr)
args <- commandArgs(TRUE)
out <- if(length(args)) args[1] else "validation/covariance_review_medium"
a <- readRDS(file.path(out,"case-12-block-01.rds")); z <- a$details[[1]]
cache <- glmmplanr:::prepare_cpp_schedule_cache(a$spec,a$schedule)
ws <- glmmplanr:::get_cached_cpp_ws(cache,a$spec$analysis_re)
t <- z$theta; t[4] <- t[4]-.01; t[2] <- t[2]+.01
beta <- z$native$beta
for(e in parse(file.path(out,"source_script.R"))) {
  if(is.call(e) && identical(e[[1]],as.name("<-")) &&
     identical(e[[2]],as.name("dense_check"))) eval(e)
}
rows <- list()
for(k in 1:2) {
  f <- glmmplanr:::cpp_fit_crossed_glmm_ws(ws,as.integer(z$y),chance=a$case$chance,
    corr_mode=0L,regularized=TRUE,fix_cov=TRUE,compute_joint_vcov=FALSE,return_ranef=TRUE,
    init_beta=beta,init_Sigma_subj=diag(exp(2*t[1:2])),init_Sigma_item=diag(exp(2*t[3:4])),
    max_em_iter=400L,max_pirls_iter=44L,tol=1.5e-4)
  check <- dense_check(f,cache,as.integer(z$y),a$case$chance)
  rows[[k]] <- data.frame(case=12L,replicate=1L,attempt=k,status=f$status,
    inner_step=f$inner_step,logLik=f$logLik,beta_intercept=f$beta[1],beta_effect=f$beta[2],
    dense_ll_difference=check[1],dense_mode_score=check[2])
  beta <- f$beta
}
result <- do.call(rbind,rows)
stopifnot(all(is.finite(result$logLik)),all(is.finite(result$dense_ll_difference)))
write.csv(result,file.path(out,"failure_inspection.csv"),row.names=FALSE)
print(result,row.names=FALSE)
cat("Unchanged main-run failures retained; a small mode step alone does not validate a profiled optimum.\n")
