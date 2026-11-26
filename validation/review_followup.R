# Focused checks of review claims; not a Monte Carlo performance study.
pkgload::load_all(quiet = TRUE)
ns <- asNamespace('glmmplanr')
cat('Source MD5:', unname(tools::md5sum('src/glmm_laplace.cpp')), '\n')
cat('t25 tail at 1.96:', 2 * pt(-1.96, 25), '\n')
cat('MCSE for p=.061, n=2048:', sqrt(.061 * .939 / 2048), '\n')
toy <- expand.grid(B = c(99L, 399L), bound_over_sigma = c(.9, 1.1))
toy$pass <- with(toy, pchisq((B-1)*bound_over_sigma^2, df=B-1))
print(toy)
cat('Normal-bootstrap SD relative MCSE B99:', 1/sqrt(196), '\n')
cat('Four Stage-2 points plus three sensitivity levels, no curtailment:', (4*1000+3*200)*100, '\n')

sparse <- design_crossed(list(condition=c('a','b')), items_per_subject=2L)
schedule <- generate_design_schedule(sparse, 48L, 4L)
tr <- schedule$trials
patterns <- vapply(split(tr, tr$subj_id), function(x) paste(x$cell_id[order(x$item_id)], x$item_id[order(x$item_id)], collapse=';'), '')
cat('Partial schedule: identical item/cell groups=', length(unique(patterns)), '; rank-1 group limit=', max(8L,48L/2L), '\n')
cat('Per-participant condition counts:\n'); print(table(rowSums(table(tr$subj_id, tr$cell_id)[,1L,drop=FALSE])))
des <- design_crossed(list(condition=c('a','b')))
one_sided <- calibrate_dgp(des, c(a=.42,b=.5), subj_sd=c(.8,.35), item_sd=c(.6,0))
sp <- specify_target(one_sided, 'condition')
cat('One-sided DGP resolves to:', sp$analysis_re, '\n'); print(formula(sp))

# A necessary stationarity check: profile beta and conditional modes while
# perturbing log SDs, holding correlations fixed. The LKJ term is constant
# along these perturbations. The objective uses log-SD coordinates, matching
# the existing diagonal-path test, not a claimed exact HMC posterior fit.
rows <- list()
for (g in c(0,1/3)) for (rho in c(0,.2)) {
  dgp <- calibrate_dgp(des, c(a=.55,b=.65), chance=g,
    subj_sd=c(.8,.35), item_sd=c(.6,.4), corr_structure=rho)
  spec <- specify_target(dgp,'condition',prob_se_method='joint_delta')
  sch <- generate_design_schedule(des,32L,24L)
  data <- get('with_rng_stream',ns)(921126L, get('simulate_dgp_data',ns)(dgp,sch))
  X <- des$X_cells[data$cell_id,,drop=FALSE]
  Zs <- des$Z_subj_cells[data$cell_id,,drop=FALSE]
  Zi <- des$Z_item_cells[data$cell_id,,drop=FALSE]
  ws <- get('cpp_create_glmm_workspace',ns)(X,data$subj_id-1L,Zs,32L,
    data$item_id-1L,Zi,24L,integer(),matrix(0,nrow(data),0L),0L)
  mode <- if (rho==0) 0L else 1L
  fit <- function(Ss=NULL,Si=NULL) get('cpp_fit_crossed_glmm_ws',ns)(ws,data$y,
    chance=g,corr_mode=mode,regularized=TRUE,lkj_eta=spec$lkj_eta,
    prior_beta_sd=spec$prior_beta_sd,prior_tau_sd=spec$prior_tau_sd,
    max_em_iter=if(g>0) 400L else 160L,max_pirls_iter=if(g>0) 44L else 22L,
    tol=1.5e-4,fix_cov=!is.null(Ss),init_Sigma_subj=Ss,init_Sigma_item=Si)
  f <- fit()
  Rs <- cov2cor(f$Sigma_subj); Ri <- cov2cor(f$Sigma_item)
  theta <- log(sqrt(c(diag(f$Sigma_subj),diag(f$Sigma_item))))
  valid <- TRUE
  objective <- function(t) {
    Ds <- diag(exp(t[1:2])); Di <- diag(exp(t[3:4]))
    z <- fit(Ds %*% Rs %*% Ds, Di %*% Ri %*% Di)
    valid <<- valid && isTRUE(z$converged)
    -z$logLik-sum(dnorm(z$beta,0,c(1.5,1)*spec$prior_beta_sd,log=TRUE))-
      sum(dt(exp(t)/spec$prior_tau_sd,3,log=TRUE))-sum(t)
  }
  gradient <- function(h) vapply(seq_along(theta),function(k){
    step <- replace(numeric(4),k,h)
    (objective(theta+step)-objective(theta-step))/(2*h)
  },numeric(1))
  a <- gradient(.001); b <- gradient(.0005)
  info_sd <- sqrt(c(f$var_log_tau_subj,f$var_log_tau_item))
  row <- data.frame(chance=g,corr=rho,converged=f$converged,reported_score=f$covariance_score_norm,
    max_scaled_logsd_gradient=max(abs(b*info_sd)),gradient_step_change=max(abs(a-b)),
    perturbations_converged=valid)
  rows[[length(rows)+1L]] <- row
}
cat('Four single-dataset stationarity probes (no accuracy-rate claim):\n')
print(do.call(rbind,rows), digits=7,row.names=FALSE)
