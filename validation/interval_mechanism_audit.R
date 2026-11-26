# Exploratory reanalysis of frozen estimates/draws; no GLMM is fitted here.
source("validation/planning_comparison_helpers.R")
input <- "validation/l2_default_bootstrap_confirm"
args <- commandArgs(TRUE)
out <- if (length(args)) args[1] else "validation/interval_mechanism_audit"
files <- list.files(input, "^job-.*[.]rds$", full.names = TRUE)
checkpoint_dir(out, c("validation/interval_mechanism_audit.R", files,
  file.path(input, c("manifest.rds", "fits.csv.gz", "summary.csv", "source_glmm_laplace.cpp"))),
  list(analysis_only = TRUE, outer_per_cell = 512L, long_draw_outer_per_cell = 32L,
    empirical_resamples = 2000L, resample_seed = 953260001L,
    quantiles = "boot.ci normal-quantile interpolation", width_pp = 10,
    match_tolerance_pp = 1e-7), c("glmmplanr", "boot"))
old_manifest <- readRDS(file.path(input, "manifest.rds"))
stopifnot(unname(tools::md5sum(file.path(input, "source_glmm_laplace.cpp"))) ==
  unname(old_manifest$files["src/glmm_laplace.cpp"]))
original <- read.csv(file.path(input, "fits.csv.gz"), colClasses=c(error="character"))
zcrit <- qnorm(.975)
keys <- c("case", "condition", "replicate")
metadata <- c("case", "profile", "n_subj", "n_item", "condition", "replicate", "response_hash")
ci <- function(k, n) if (n) unname(glmmplanr:::wilson_mc_ci(k, n)) else c(NA_real_, NA_real_)
gate <- function(x, p, lower = TRUE) {
  if (anyNA(x)) return("indeterminate")
  if (lower) { if (x[1] >= p) "pass" else if (x[2] < p) "fail" else "indeterminate" }
  else { if (x[2] <= p) "pass" else if (x[1] > p) "fail" else "indeterminate" }
}
interval_row <- function(v, method, lo, hi, valid, se = NA_real_) {
  valid <- valid && all(is.finite(c(lo, hi))) && lo < hi
  data.frame(v[1, metadata], method, truth = v$truth[1], estimate = v$estimate[1], se,
    valid, lower = lo, upper = hi, width = hi-lo,
    covered = valid && lo <= v$truth[1] && hi >= v$truth[1],
    detected = valid && (lo > 0 || hi < 0),
    width_ok = valid && hi-lo <= 10,
    joint = valid && (lo > 0 || hi < 0) && hi-lo <= 10)
}
saved <- intervals <- pools <- list()
for (file in files) {
  r <- readRDS(file)
  saved[[length(saved)+1L]] <- r$trials
  for (id in names(r$bootstrap_draws)) {
    a <- subset(r$trials, replicate == as.integer(id))
    jd <- subset(a, method == "joint_delta")
    stopifnot(nrow(jd) == 1L)
    e <- jd$estimate; b_all <- 100*r$bootstrap_draws[[id]]
    intervals[[length(intervals)+1L]] <- interval_row(jd, "joint_delta",
      e-zcrit*jd$se, e+zcrit*jd$se, jd$valid, jd$se)
    for (B in if (length(b_all) == 399L) c(99L, 399L) else 99L) {
      b <- head(b_all, B); good <- is.finite(b); b <- b[good]
      old <- subset(a, method == paste0("bootstrap_", B))
      valid <- jd$valid && sum(good) >= ceiling(.9*B)
      se <- if (valid) sd(b) else NA_real_
      stopifnot(nrow(old) == 1L, sum(good) == old$n_inner_valid,
        identical(valid, old$valid), !valid || abs(se-old$se) < 1e-7)
      intervals[[length(intervals)+1L]] <- interval_row(jd, paste0("wald_", B),
        e-zcrit*se, e+zcrit*se, valid, se)
      # Assemble a boot object from saved parametric draws, without boot::boot().
      if (valid) {
        bc <- boot::boot.ci(structure(list(t0=e, t=matrix(b,ncol=1L), R=length(b)),
          class="boot"), type=c("perc", "basic"), index=1L)
        perc <- unname(bc$percent[4:5]); basic <- unname(bc$basic[4:5])
        stopifnot(max(abs(basic - (2*e-rev(perc)))) < 1e-10)
      } else perc <- basic <- c(NA_real_, NA_real_)
      intervals[[length(intervals)+1L]] <- interval_row(jd, paste0("percentile_", B),
        perc[1], perc[2], valid)
      intervals[[length(intervals)+1L]] <- interval_row(jd, paste0("basic_", B),
        basic[1], basic[2], valid)
    }
    if (length(b_all) == 399L) pools[[length(pools)+1L]] <- list(row=jd, draws=b_all)
  }
}
saved <- do.call(rbind, saved)
ord <- function(d) do.call(order, unname(d[c(keys, "method")]))
stopifnot(isTRUE(all.equal(saved[ord(saved), ], original[ord(original), ],
  check.attributes=FALSE, tolerance=1e-12)), length(pools) == 192L)
iv <- do.call(rbind, intervals)
stopifnot(nrow(iv) == 3072L*4L + 192L*3L, !anyDuplicated(iv[c(keys,"method")]))
write.csv(iv, gzfile(file.path(out,"intervals.csv.gz")), row.names=FALSE)
summarize_intervals <- function(v) {
  ok <- v$valid; n <- nrow(v); nv <- sum(ok)
  cc <- ci(sum(v$covered), nv); dc <- ci(sum(v$detected), n); jc <- ci(sum(v$joint), n)
  data.frame(v[1,c("case","profile","n_subj","n_item","condition","method")],
    n, n_valid=nv, coverage=sum(v$covered)/nv, coverage_lo=cc[1], coverage_hi=cc[2],
    n_covered=sum(v$covered), detection=mean(v$detected), detection_lo=dc[1], detection_hi=dc[2],
    n_detected=sum(v$detected), width_attainment=mean(v$width_ok), n_width=sum(v$width_ok),
    joint=mean(v$joint), joint_lo=jc[1], joint_hi=jc[2], n_joint=sum(v$joint),
    mean_width=mean(v$width[ok]), coverage_floor=gate(cc,.925),
    null_ceiling=if(v$condition[1]=="null") gate(dc,.075,FALSE) else NA_character_,
    joint_decision=gate(jc,.8))
}
summary <- do.call(rbind,lapply(split(iv,interaction(iv$case,iv$condition,iv$method)),summarize_intervals))
write.csv(summary,file.path(out,"interval_summary.csv"),row.names=FALSE)
old <- read.csv(file.path(input,"summary.csv"))
old$method <- sub("bootstrap_","wald_",old$method)
match <- merge(summary,old,by=c("case","condition","method"),suffixes=c("_new","_old"))
stopifnot(nrow(match)==18L,all(match$n==match$nsim),all(match$n_valid_new==match$n_valid_old),
  all(match$n_covered_new==match$n_covered_old),all(match$n_detected_new==match$n_detected_old),
  all(match$n_joint_new==match$n_joint_old),all(match$joint_decision_new==match$joint_decision_old))

paired <- function(a,b) {
  p <- merge(subset(iv,method==a),subset(iv,method==b),by=keys,suffixes=c("_a","_b"))
  stopifnot(all(p$response_hash_a==p$response_hash_b))
  do.call(rbind,lapply(split(p,interaction(p$case,p$condition)),function(v) {
    delta <- function(name) as.numeric(v[[paste0(name,"_b")]])-as.numeric(v[[paste0(name,"_a")]])
    d <- delta("joint"); c <- delta("covered"); r <- delta("detected")
    data.frame(v[1,c("case","condition")],a,b,n=nrow(v),n_both_valid=sum(v$valid_a&v$valid_b),
      joint_change=mean(d),joint_mcse=sd(d)/sqrt(length(d)),joint_lost=sum(d<0),joint_gained=sum(d>0),
      lost_width=sum(d<0 & !v$width_ok_b & v$detected_b),
      lost_detection=sum(d<0 & v$width_ok_b & !v$detected_b),
      lost_both=sum(d<0 & !v$width_ok_b & !v$detected_b),
      coverage_change=mean(c),coverage_mcse=sd(c)/sqrt(length(c)),
      detection_change=mean(r),detection_mcse=sd(r)/sqrt(length(r)))
  }))
}
pairs <- rbind(paired("joint_delta","wald_99"), paired("wald_99","percentile_99"),
  paired("wald_99","basic_99"),paired("wald_99","wald_399"),
  paired("wald_399","percentile_399"),paired("wald_399","basic_399"))
write.csv(pairs,file.path(out,"paired.csv"),row.names=FALSE)

w <- subset(iv,method %in% c("joint_delta","wald_99"))
diagnostics <- strata <- list()
for (v in split(w,interaction(w$case,w$condition,w$method))) {
  v <- subset(v,valid); err <- v$estimate-v$truth; bias <- mean(err)
  empirical_sd <- sd(err); rms <- sqrt(mean(v$se^2)); z <- err/v$se
  diagnostics[[length(diagnostics)+1L]] <- data.frame(v[1,c("case","condition","method")],
    n=nrow(v),bias,bias_mcse=empirical_sd/sqrt(nrow(v)),empirical_sd,rms_se=rms,
    rms_se_over_sd=rms/empirical_sd,z_mean=mean(z),z_sd=sd(z),
    z_q025=unname(quantile(z,.025)),z_q975=unname(quantile(z,.975)),
    z_left=sum(z < -zcrit),z_right=sum(z > zcrit),corr_abs_error_se=cor(abs(err),v$se),
    oracle_debiased_native_miss=mean(abs(err-bias)>zcrit*v$se),
    oracle_constant_sd_miss=mean(abs(err)>zcrit*empirical_sd),
    oracle_debiased_constant_sd_miss=mean(abs(err-bias)>zcrit*empirical_sd))
  groups <- ceiling(rank(v$se,ties.method="first")/nrow(v)*5L)
  for (g in 1:5) {
    take <- groups==g
    strata[[length(strata)+1L]] <- data.frame(v[1,c("case","condition","method")],
      se_quintile=g,n=sum(take),bias=mean(err[take]),rmse=sqrt(mean(err[take]^2)),
      rms_se=sqrt(mean(v$se[take]^2)),noncoverage=mean(abs(z[take])>zcrit))
  }
}
diagnostics <- do.call(rbind,diagnostics)
write.csv(diagnostics,file.path(out,"error_diagnostics.csv"),row.names=FALSE)
write.csv(do.call(rbind,strata),file.path(out,"se_strata.csv"),row.names=FALSE)

# Conditional finite-pool experiments: these are NOT extra outer replications.
RNGkind("L'Ecuyer-CMRG","Inversion","Rejection"); set.seed(953260001L)
M <- 2000L; block_rows <- pool_rows <- list()
for (p in pools) {
  v <- p$row; e <- v$estimate; b <- p$draws
  stopifnot(all(is.finite(b)))
  for (k in 1:4) {
    se <- sd(b[seq.int((k-1L)*99L+1L,k*99L)])
    block_rows[[length(block_rows)+1L]] <- cbind(interval_row(v,"wald_block_99",
      e-zcrit*se,e+zcrit*se,v$valid,se),block=k)
  }
  draws <- matrix(sample(b,399L*M,replace=TRUE),nrow=399L)
  se99 <- apply(draws[1:99,,drop=FALSE],2,sd); se399 <- apply(draws,2,sd)
  det99 <- abs(e)>zcrit*se99; det399 <- abs(e)>zcrit*se399
  width99 <- 2*zcrit*se99<=10; width399 <- 2*zcrit*se399<=10
  joint99 <- det99&width99; joint399 <- det399&width399
  d <- as.numeric(joint399)-as.numeric(joint99)
  pool_rows[[length(pool_rows)+1L]] <- data.frame(v[1,metadata],estimate=e,
    observed_se399=sd(b),pool_mean=mean(b),pool_bias=mean(b)-e,
    observed_joint399=abs(e)>zcrit*sd(b) && 2*zcrit*sd(b)<=10,
    p_joint99=mean(joint99),p_joint399=mean(joint399),
    p_width99=mean(width99),p_width399=mean(width399),
    p_detect99=mean(det99),p_detect399=mean(det399),
    conditional_change_mcse=sd(d)/sqrt(M),
    se99_cv=sd(se99)/mean(se99),se399_cv=sd(se399)/mean(se399),
    conditional_disagreement=mean(joint99!=joint399))
}
blocks <- do.call(rbind,block_rows); pool <- do.call(rbind,pool_rows)
write.csv(blocks,file.path(out,"disjoint_blocks.csv"),row.names=FALSE)
write.csv(pool,file.path(out,"finite_pool.csv"),row.names=FALSE)
pool_summary <- do.call(rbind,lapply(split(pool,interaction(pool$case,pool$condition)),function(v) {
  d <- v$p_joint399-v$p_joint99
  data.frame(v[1,c("case","condition")],n_outer=nrow(v),resamples_per_outer=M,
    observed_joint399=mean(v$observed_joint399),mean_p_joint99=mean(v$p_joint99),
    mean_p_joint399=mean(v$p_joint399),mean_change=mean(d),
    conditional_change_mcse=sqrt(sum(v$conditional_change_mcse^2))/nrow(v),
    outer_se_of_change=sd(d)/sqrt(nrow(v)),mean_se99_cv=mean(v$se99_cv),
    mean_se399_cv=mean(v$se399_cv),mean_disagreement=mean(v$conditional_disagreement),
    n_ambiguous99=sum(v$p_joint99>.1 & v$p_joint99<.9),
    n_ambiguous399=sum(v$p_joint399>.1 & v$p_joint399<.9))
}))
write.csv(pool_summary,file.path(out,"finite_pool_summary.csv"),row.names=FALSE)
block_summary <- do.call(rbind,lapply(split(blocks,interaction(blocks$case,blocks$condition,blocks$block)),
  function(v)cbind(summarize_intervals(v),block=v$block[1])))
write.csv(block_summary,file.path(out,"block_summary.csv"),row.names=FALSE)
first <- subset(blocks,block==1L)
ref <- subset(iv,method=="wald_99" & replicate<=32L)
p <- merge(first,ref,by=keys,suffixes=c("_a","_b"))
stopifnot(nrow(p)==192L,all(p$joint_a==p$joint_b),max(abs(p$se_a-p$se_b))<1e-10)

library(ggplot2)
main <- subset(summary,method %in% c("joint_delta","wald_99","percentile_99","basic_99"))
panel <- function(v,label,metric) data.frame(v[c("case","method")],metric=label,
  rate=v[[metric]],lo=v[[paste0(metric,"_lo")]],hi=v[[paste0(metric,"_hi")]])
pd <- rbind(panel(subset(main,condition=="null"),"Null rejection","detection"),
  panel(subset(main,condition=="planned"),"Planned coverage","coverage"),
  panel(subset(main,condition=="planned"),"Detection AND width <=10pp","joint"))
pd$design <- factor(pd$case,c(25,50,31),c("Larger variances: 162 x 30","Moderate slopes: 72 x 60","Intercepts only: 54 x 30"))
pd$method <- factor(pd$method,c("joint_delta","wald_99","percentile_99","basic_99"),
  c("Joint delta","Bootstrap Wald (99)","Percentile (99)","Basic (99)"))
refs <- data.frame(metric=c("Null rejection","Planned coverage","Detection AND width <=10pp"),value=c(.05,.95,.8))
fig <- ggplot(pd,aes(rate,design,color=method))+
  geom_vline(data=refs,aes(xintercept=value),linetype="dashed",color="grey55")+
  geom_errorbar(aes(xmin=lo,xmax=hi),orientation="y",position=position_dodge(.65),width=.12)+
  geom_point(position=position_dodge(.65),size=2)+facet_wrap(~metric,scales="free_x")+
  scale_x_continuous(labels=function(x)paste0(round(100*x),"%"))+
  scale_color_manual(values=c("#555555","#0072B2","#D55E00","#009E73"))+
  labs(x=NULL,y=NULL,color=NULL)+theme_bw(base_size=11)+theme(legend.position="bottom")
ggsave(file.path(out,"interval_comparison.png"),fig,width=13,height=5,dpi=150)
ggsave(file.path(out,"interval_comparison.pdf"),fig,width=13,height=5)
print(diagnostics,row.names=FALSE)
print(subset(pairs,condition=="planned" & case==50),row.names=FALSE)
print(pool_summary,row.names=FALSE)
cat("All saved counts, draw SEs, original summary decisions and disjoint first blocks match.\n")
cat("No GLMM fits or new outer response datasets; alternative intervals and pool experiments are exploratory.\n")
