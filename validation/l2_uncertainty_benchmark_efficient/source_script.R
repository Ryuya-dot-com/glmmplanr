# Matched response/CPU check, including paths deliberately left unchanged.
source("validation/planning_comparison_helpers.R")
args<-commandArgs(TRUE); out<-args[1]; variant<-args[2]
for(e in parse("validation/l2_replication_large/source_script.R")) {
  if(is.call(e)&&identical(e[[1]],as.name("<-"))&&is.symbol(e[[2]])&&
    as.character(e[[2]])%in%c("profiles","make_spec")) eval(e)
}
chosen<-read.csv("validation/l2_uncertainty_oracle/cases.csv"); chosen$kind<-"uncorrelated"
controls<-data.frame(profile=c("correlated","guessing"),n_subj=48L,n_item=32L,
  case=c(101L,102L),kind=c("correlated","guessing"))
jobs<-rbind(chosen,controls); jobs$job<-seq_len(nrow(jobs))
checkpoint_dir(out,c("validation/l2_uncertainty_benchmark.R","R/marginal_scale.R"),
  list(jobs=jobs,nsim=16L,seed=913100001L,variant=variant),c("glmmplanr","digest"))
z<-do.call(rbind,lapply(seq_len(nrow(jobs)),function(k) {
  x<-jobs[k,]
  if(x$kind=="uncorrelated") spec<-make_spec(profiles[profiles$profile==x$profile,]) else {
    des<-design_crossed(list(condition=c("control","treatment")))
    dgp<-calibrate_dgp(des,c(control=.55,treatment=.65),subj_sd=c(.8,.3),item_sd=c(.6,.4),
      corr_structure=if(x$kind=="correlated") .2 else "zero",chance=if(x$kind=="guessing") 1/3 else 0)
    spec<-specify_target(dgp,"condition",scale="pp",engine="bayes_laplace",analysis_re="match_dgp",
      prob_se_method="joint_delta")
  }
  schedule<-generate_design_schedule(spec$dgp$design,x$n_subj,x$n_item)
  cache<-glmmplanr:::prepare_cpp_schedule_cache(spec,schedule)
  streams<-glmmplanr:::make_lecuyer_streams(913100001L+1000L*k,16L)
  do.call(rbind,lapply(seq_len(16L),function(i) {
    sim<-glmmplanr:::with_rng_stream(streams[[i]],glmmplanr:::simulate_fast_y_keep(cache))
    t<-system.time(f<-glmmplanr:::fit_via_cpp_cached(sim,cache,spec,TRUE,spec$analysis_re,FALSE))
    data.frame(x,variant=variant,replicate=i,valid=f$converged,
      estimate=100*f$prob_est,se=100*f$prob_se,
      cpu=unname(t["user.self"]+t["sys.self"]),elapsed=unname(t["elapsed"]),
      hash=digest::digest(sim$y,algo="xxhash64"))
  }))
}))
write.csv(z,file.path(out,"fits.csv"),row.names=FALSE)
file.copy("validation/l2_uncertainty_benchmark.R",file.path(out,"source_script.R"),overwrite=TRUE)
file.copy("src/glmm_laplace.cpp",file.path(out,"source_glmm_laplace.cpp"),overwrite=TRUE)
print(aggregate(cbind(cpu,elapsed,valid)~profile,z,mean),row.names=FALSE)
