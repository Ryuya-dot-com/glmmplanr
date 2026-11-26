# Exploratory full profile curvature; compares the SAME native likelihood and priors.
source("validation/planning_comparison_helpers.R")
args <- commandArgs(TRUE)
nsim <- if (length(args)) as.integer(args[1]) else 4L
out <- if (length(args)>1) args[2] else "validation/l2_uncertainty_profile_pilot"
workers <- if (length(args)>2) as.integer(args[3]) else 4L
stopifnot(nsim>0L,nsim<=512L)
if (.Platform$OS.type=="windows") workers <- 1L
for (e in parse("validation/l2_replication_large/source_script.R")) {
  if (is.call(e) && identical(e[[1]],as.name("<-")) && is.symbol(e[[2]]) &&
      as.character(e[[2]]) %in% c("profiles","make_spec")) eval(e)
}
chosen <- read.csv("validation/l2_uncertainty_oracle/cases.csv")
jobs <- merge(chosen,data.frame(condition=c("null","planned")))
jobs$job <- seq_len(nrow(jobs))
checkpoint_dir(out,c("validation/l2_uncertainty_profile.R","R/marginal_scale.R",
  "validation/l2_replication_large/source_script.R"),list(nsim=nsim,jobs=jobs,h=.01),c("glmmplanr","digest"))
run <- function(k) {
  path <- file.path(out,sprintf("case-%02d.rds",k)); if(file.exists(path)) return(readRDS(path)$summary)
  x<-jobs[k,]; spec<-make_spec(profiles[profiles$profile==x$profile,])
  if(x$condition=="null") spec<-glmmplanr:::null_target_spec(spec)
  schedule<-generate_design_schedule(spec$dgp$design,x$n_subj,x$n_item)
  cache<-glmmplanr:::prepare_cpp_schedule_cache(spec,schedule)
  ws<-glmmplanr:::get_cached_cpp_ws(cache,spec$analysis_re)
  bases<-glmmplanr:::allocate_stage_point_streams(810900001L+10000L*x$case,1L,1:2,2L)
  streams<-glmmplanr:::derive_rep_substreams(bases[[if(x$condition=="null") 1 else 2]],nsim)
  details<-vector("list",nsim); rows<-vector("list",nsim)
  for(i in seq_len(nsim)) {
    sim<-glmmplanr:::with_rng_stream(streams[[i]],glmmplanr:::simulate_fast_y_keep(cache))
    current<-glmmplanr:::fit_via_cpp_cached(sim,cache,spec,TRUE,spec$analysis_re,FALSE)
    qs<-nrow(current$Sigma_subj); qi<-nrow(current$Sigma_item)
    theta<-log(sqrt(c(diag(current$Sigma_subj),diag(current$Sigma_item)))); d<-length(theta)
    fit<-function(t,joint=FALSE) glmmplanr:::cpp_fit_crossed_glmm_ws(ws,sim$y,chance=0,
      corr_mode=0L,regularized=TRUE,lkj_eta=spec$lkj_eta,prior_beta_sd=spec$prior_beta_sd,
      prior_tau_sd=spec$prior_tau_sd,max_em_iter=160L,max_pirls_iter=22L,tol=1.5e-4,
      init_beta=current$beta_full,init_Sigma_subj=diag(exp(2*t[seq_len(qs)]),qs),
      init_Sigma_item=diag(exp(2*t[qs+seq_len(qi)]),qi),fix_cov=TRUE,compute_joint_vcov=joint)
    objective<-function(t,f) {
      # Native covariance MAP is in log-SD coordinates: retain the +sum(theta) Jacobian.
      -f$logLik-sum(dnorm(f$beta,0,c(1.5*spec$prior_beta_sd,rep(spec$prior_beta_sd,length(f$beta)-1)),log=TRUE))-
        sum(dt(exp(t)/spec$prior_tau_sd,3,log=TRUE))-sum(t)
    }
    effect<-function(f) glmmplanr:::estimate_prob_effect_delta(f$beta,f$vcov,f$Sigma_subj,
      f$Sigma_item,NULL,spec$dgp$design,0,spec$contrast_weights,"marginal",spec$analysis_re,spec$dgp$n_quad)
    f0<-fit(theta,TRUE); stopifnot(f0$converged); p0<-effect(f0); value0<-objective(theta,f0)
    curvature<-function(h) {
      H<-matrix(0,d,d); gradient<-numeric(d); gp<-numeric(d); J<-matrix(0,length(f0$beta),d)
      eval_at<-function(t) { f<-fit(t); stopifnot(f$converged); list(v=objective(t,f),p=effect(f)$est,b=f$beta) }
      for(a in seq_len(d)) {
        step<-replace(numeric(d),a,h); up<-eval_at(theta+step); dn<-eval_at(theta-step)
        H[a,a]<-(up$v-2*value0+dn$v)/h^2; gradient[a]<-(up$v-dn$v)/(2*h)
        gp[a]<-(up$p-dn$p)/(2*h); J[,a]<-(up$b-dn$b)/(2*h)
        if(a>1L) for(b in seq_len(a-1L)) {
          other<-replace(numeric(d),b,h)
          H[a,b]<-H[b,a]<-(eval_at(theta+step+other)$v-eval_at(theta+step-other)$v-
            eval_at(theta-step+other)$v+eval_at(theta-step-other)$v)/(4*h^2)
        }
      }
      eig<-eigen(H,symmetric=TRUE,only.values=TRUE)$values
      V<-if(min(eig)>1e-6) solve(H) else matrix(NA_real_,d,d)
      se<-if(all(is.finite(V))) sqrt(p0$se^2+drop(crossprod(gp,V%*%gp))) else NA_real_
      list(H=H,V=V,gradient=gradient,g_prob_profile=gp,J=J,se=se,min_eigen=min(eig))
    }
    timing<-system.time(curve<-curvature(.01))["elapsed"]
    fine<-if(i==1L) curvature(.005) else NULL
    rows[[i]]<-data.frame(x,replicate=i,truth=spec$true_effect_target,
      current_estimate=100*current$prob_est,current_se=100*current$prob_se,
      profile_estimate=100*p0$est,fixed_estimated_se=100*p0$se,profile_se=100*curve$se,
      score_max=max(abs(curve$gradient)),min_eigen=curve$min_eigen,
      J_max_difference=max(abs(curve$J-current$J_beta_theta)),
      fine_se_difference=if(!is.null(fine)) 100*(curve$se-fine$se) else NA_real_,
      elapsed=unname(timing),response_hash=digest::digest(sim$y,algo="xxhash64"))
    details[[i]]<-list(current=current,profile=f0,theta=theta,curve=curve,fine=fine)
  }
  summary<-do.call(rbind,rows); atomic_rds(list(summary=summary,details=details),path)
  cat(sprintf("Profile curvature case %d %s: %d trials complete\n",x$case,x$condition,nsim)); summary
}
results<-parallel::mclapply(seq_len(nrow(jobs)),run,mc.cores=workers,mc.set.seed=FALSE,mc.preschedule=FALSE)
stopifnot(all(vapply(results,is.data.frame,logical(1))))
summary<-do.call(rbind,results); write.csv(summary,file.path(out,"summary.csv"),row.names=FALSE)
file.copy("validation/l2_uncertainty_profile.R",file.path(out,"source_script.R"),overwrite=TRUE)
print(summary,row.names=FALSE)
