# Paired numerical-path study; the profile fit is a diagnostic, not a public engine.
source("validation/planning_comparison_helpers.R")
mode <- match.arg(Sys.getenv("GLMMPLANR_REVIEW_MODE", "medium"), c("smoke", "medium"))
out <- paste0("validation/covariance_review_", mode)
nsim <- if (mode == "smoke") 1L else 128L
workers <- if (.Platform$OS.type == "windows") 1L else 6L
seed <- if (mode == "smoke") 951260001L else 952260001L
cases <- expand.grid(layout = c("complete", "partial"), chance = c(0, 1/3),
  correlation = c(0, .2), effect_pp = c(0, 10), stringsAsFactors = FALSE)
cases$case <- seq_len(nrow(cases))
cases$n_subj <- ifelse(cases$layout == "complete", 32L, 72L)
cases$n_item <- ifelse(cases$layout == "complete", 24L, 60L)
cases$per_person <- ifelse(cases$layout == "complete", 24L, 30L)
settings <- list(mode = mode, nsim = nsim, cases = cases, seed = seed, item_seed = 37L,
  h = .01, h_fine = .005, score_tolerance = .01, dense_ll_tolerance = 1e-6,
  width_pp = 20, assurance = .8, alpha = .05,
  optimizer = list(method = "L-BFGS-B", maxit = 200L, factr = 1e6, pgtol = 1e-4),
  coordinate = "log SD, Fisher rho; density in log SD and rho (no Fisher Jacobian)")
checkpoint_dir(out, c("validation/covariance_review_medium.R", "R/calibrate.R",
  "R/marginal_scale.R", "R/utils_math.R", "R/design.R"), settings,
  c("glmmplanr", "digest", "Rcpp", "RcppEigen"))
write.csv(cases, file.path(out, "cases.csv"), row.names = FALSE)
file.copy("validation/covariance_review_medium.R", file.path(out, "source_script.R"), overwrite = TRUE)
file.copy("src/glmm_laplace.cpp", file.path(out, "source_glmm_laplace.cpp"), overwrite = TRUE)

make_case <- function(x) {
  des <- design_crossed(list(condition = c("a", "b")), allocation = "latin_square",
    items_per_subject = x$per_person, item_seed = 37L)
  dgp <- calibrate_dgp(des, c(a = .55, b = .55 + x$effect_pp/100), chance = x$chance,
    subj_sd = c(.8, .35), item_sd = c(.6, .4), corr_structure = x$correlation)
  spec <- specify_target(dgp, "condition", scale = "pp", criterion = "joint",
    ci_width = 20, prob_se_method = "joint_delta", engine = "bayes_laplace")
  schedule <- generate_design_schedule(des, x$n_subj, x$n_item)
  list(spec = spec, schedule = schedule,
    cache = glmmplanr:::prepare_cpp_schedule_cache(spec, schedule))
}

unpack <- function(t) {
  rho <- if (length(t) == 6L) tanh(t[5:6]) else c(0, 0)
  covar <- function(sd, r) outer(sd, sd) * matrix(c(1, r, r, 1), 2)
  list(Ss = covar(exp(t[1:2]), rho[1]), Si = covar(exp(t[3:4]), rho[2]), rho = rho)
}

dense_check <- function(f, cache, y, chance) {
  ns <- cache$n_s_full; ni <- cache$n_i_full
  Ls <- t(chol(f$Sigma_subj)); Li <- t(chol(f$Sigma_item))
  Z <- matrix(0, length(y), 2L*(ns + ni))
  zs <- cache$zs_dgp %*% Ls; zi <- cache$zi_dgp %*% Li
  for (j in 1:2) {
    Z[cbind(seq_along(y), 2L*cache$s_idx0_full + j)] <- zs[,j]
    Z[cbind(seq_along(y), 2L*ns + 2L*cache$i_idx0_full + j)] <- zi[,j]
  }
  v <- c(as.vector(solve(Ls, t(f$u_subj))), as.vector(solve(Li, t(f$u_item))))
  eta <- drop(cache$X_full %*% f$beta + Z %*% v)
  p <- plogis(eta); mu <- chance + (1-chance)*p
  dp <- (1-chance)*p*(1-p); d2p <- dp*(1-2*p)
  score <- (y-mu)*dp/(mu*(1-mu))
  w <- (dp^2-(y-mu)*d2p)/(mu*(1-mu)) +
    (y-mu)*dp^2*(1-2*mu)/(mu^2*(1-mu)^2)
  H <- diag(length(v)) + crossprod(Z, w*Z)
  ll <- sum(dbinom(y, 1, mu, log = TRUE)) - .5*sum(v^2) - sum(log(diag(chol(H))))
  c(ll_difference = f$logLik-ll, mode_score = max(abs(drop(crossprod(Z, score))-v)))
}

run_trial <- function(x, i, obj, stream) {
  started <- proc.time()[[3]]
  spec <- obj$spec; cache <- obj$cache; des <- spec$dgp$design
  sim <- glmmplanr:::with_rng_stream(stream, glmmplanr:::simulate_fast_y_keep(cache))
  ws <- glmmplanr:::get_cached_cpp_ws(cache, spec$analysis_re)
  corr <- if (x$correlation == 0) 0L else 1L
  core <- function(t = NULL, joint = TRUE, ranef = FALSE, beta = NULL, early = FALSE) {
    cv <- if (!is.null(t)) unpack(t) else NULL
    coarse <- early && !(x$chance == 0 && corr == 0L)
    glmmplanr:::cpp_fit_crossed_glmm_ws(ws, sim$y, chance = x$chance, corr_mode = corr,
      regularized = TRUE, lkj_eta = spec$lkj_eta, prior_beta_sd = spec$prior_beta_sd,
      prior_tau_sd = spec$prior_tau_sd,
      max_em_iter = if (coarse) 85L else if (x$chance > 0) 400L else 160L,
      max_pirls_iter = if (coarse) 15L else if (x$chance > 0) 44L else 22L,
      tol = if (coarse) 6e-4 else 1.5e-4, init_beta = beta,
      init_Sigma_subj = cv$Ss, init_Sigma_item = cv$Si,
      fix_cov = !is.null(t), compute_joint_vcov = joint, return_ranef = ranef)
  }
  point <- function(f, joint = FALSE) glmmplanr:::estimate_prob_effect_delta(
    f$beta, f$vcov, f$Sigma_subj, f$Sigma_item, NULL, des, x$chance,
    spec$contrast_weights, "marginal", spec$analysis_re, spec$dgp$n_quad,
    var_log_tau_subj = if (joint) f$var_log_tau_subj else NULL,
    var_log_tau_item = if (joint) f$var_log_tau_item else NULL,
    J_beta_theta = if (joint) f$J_beta_theta else NULL)
  error_text <- character()
  attempt <- function(expr) tryCatch(expr, error = function(e) {
    error_text <<- c(error_text, conditionMessage(e)); NULL
  })
  native <- attempt(core(ranef = i == 1L))
  early <- attempt(core(early = TRUE))
  p_native <- if (!is.null(native)) attempt(point(native, TRUE)) else NULL
  start <- if (!is.null(native) && all(is.finite(native$beta))) native else
    list(beta = c(0,0), Sigma_subj = diag(c(.8,.35)^2), Sigma_item = diag(c(.6,.4)^2))
  t0 <- log(sqrt(c(diag(start$Sigma_subj), diag(start$Sigma_item))))
  if (corr == 1L) t0 <- c(t0, atanh(c(cov2cor(start$Sigma_subj)[1,2],
    cov2cor(start$Sigma_item)[1,2])))
  lower <- c(rep(-5,4), if (corr == 1L) rep(-atanh(.98),2))
  upper <- c(rep(2,4), if (corr == 1L) rep(atanh(.98),2))
  t0 <- pmax(lower + 1e-5, pmin(upper - 1e-5, t0))
  evaluations <- 0L; failed_evaluations <- 0L
  evaluate <- function(t, joint = FALSE, ranef = FALSE) {
    evaluations <<- evaluations + 1L
    f <- attempt(core(t, joint, ranef, beta = start$beta))
    if (is.null(f) || !isTRUE(f$converged) || !is.finite(f$logLik)) {
      failed_evaluations <<- failed_evaluations + 1L
      return(list(value = 1e12, fit = NULL))
    }
    rho <- unpack(t)$rho
    prior <- sum(dnorm(f$beta, 0, c(1.5,1)*spec$prior_beta_sd, log = TRUE)) +
      sum(dt(exp(t[1:4])/spec$prior_tau_sd, 3, log = TRUE)) + sum(t[1:4]) +
      (spec$lkj_eta-1)*sum(log1p(-rho^2))
    list(value = -f$logLik-prior, fit = f)
  }
  objective <- function(t) evaluate(t)$value
  gradient <- function(t, h = .001) vapply(seq_along(t), function(j) {
    step <- replace(numeric(length(t)), j, h)
    (objective(t+step)-objective(t-step))/(2*h)
  }, numeric(1))
  value0 <- objective(t0); score0 <- gradient(t0)
  fitopt <- attempt(optim(t0, objective, method = "L-BFGS-B", lower = lower, upper = upper,
    control = list(maxit = 200L, factr = 1e6, pgtol = 1e-4, ndeps = rep(.001,length(t0)))))
  profile <- NULL; curve <- NULL; fine <- NULL; alt <- NULL
  if (!is.null(fitopt)) {
    t <- fitopt$par; profile <- evaluate(t, TRUE, i == 1L)
    if (!is.null(profile$fit)) {
      curve_at <- function(h) {
        d <- length(t); H <- matrix(0,d,d); gp <- numeric(d)
        J <- matrix(0,2,d); grad <- numeric(d); all_ok <- TRUE
        at <- function(v) {
          z <- evaluate(v); all_ok <<- all_ok && !is.null(z$fit)
          list(v = z$value, p = if (!is.null(z$fit)) point(z$fit)$est else NA_real_,
            b = if (!is.null(z$fit)) z$fit$beta else c(NA_real_, NA_real_))
        }
        for (a in seq_len(d)) {
          step <- replace(numeric(d),a,h); up <- at(t+step); dn <- at(t-step)
          H[a,a] <- (up$v-2*profile$value+dn$v)/h^2
          grad[a] <- (up$v-dn$v)/(2*h); gp[a] <- (up$p-dn$p)/(2*h)
          J[,a] <- (up$b-dn$b)/(2*h)
          if (a>1L) for (b in seq_len(a-1L)) {
            other <- replace(numeric(d),b,h)
            H[a,b] <- H[b,a] <- (at(t+step+other)$v-at(t+step-other)$v-
              at(t-step+other)$v+at(t-step-other)$v)/(4*h^2)
          }
        }
        eig <- eigen(H,symmetric=TRUE,only.values=TRUE)$values
        V <- if (all_ok && min(eig)>1e-6) solve(H) else matrix(NA_real_,d,d)
        fixed <- point(profile$fit)
        se <- sqrt(fixed$se^2 + drop(crossprod(gp,V%*%gp)))
        beta_se <- sqrt(diag(profile$fit$vcov + J%*%V%*%t(J)))
        list(H=H, V=V, gradient=grad, se=se, beta_se=beta_se,
          point=fixed$est, fixed_se=fixed$se, J=J, gp=gp, all_ok=all_ok,
          min_eigen=min(eig), scaled_score=if(all(is.finite(V)))
            max(abs(grad)*sqrt(diag(V))) else Inf)
      }
      curve <- attempt(curve_at(.01))
      if (i == 1L) {
        fine <- attempt(curve_at(.005))
        alt_start <- t0; if(corr == 1L) alt_start[5:6] <- 0
        alt_start[1:4] <- alt_start[1:4] + c(.15,-.15,.15,-.15)
        alt <- attempt(optim(pmax(lower+1e-5,pmin(upper-1e-5,alt_start)), objective,
          method="L-BFGS-B",lower=lower,upper=upper,
          control=list(maxit=200L,factr=1e6,pgtol=1e-4,ndeps=rep(.001,length(t)))))
      }
    }
  }
  boundary <- if (!is.null(fitopt)) any(fitopt$par-lower < 1e-4 | upper-fitopt$par < 1e-4) else NA
  native_ok <- !is.null(native) && isTRUE(native$converged) && !is.null(p_native) &&
    is.finite(p_native$se) && p_native$se>0
  profile_ok <- !is.null(curve) && isTRUE(curve$all_ok) && is.finite(curve$se) &&
    curve$se>0 && !isTRUE(boundary) && curve$scaled_score < .01 &&
    profile$value <= value0+1e-5 && fitopt$convergence == 0L
  make_row <- function(method, valid, est, se, beta, beta_se) {
    z <- qnorm(.975); covered <- valid && abs(est-x$effect_pp)<=z*se
    detected <- valid && abs(est)>z*se; width_ok <- valid && 2*z*se<=20
    data.frame(x, replicate=i, method=method, valid=valid, estimate=est, se=se,
      truth=x$effect_pp, covered=covered, detected=detected, width_ok=width_ok,
      joint=detected && width_ok, logit=beta, logit_se=beta_se,
      logit_truth=unname(spec$dgp$beta[2]),
      logit_covered=valid && abs(beta-spec$dgp$beta[2])<=z*beta_se,
      response_hash=digest::digest(sim$y,algo="xxhash64"))
  }
  rows <- rbind(make_row("native", native_ok,
    if(!is.null(p_native)) 100*p_native$est else NA_real_,
    if(!is.null(p_native)) 100*p_native$se else NA_real_,
    if(!is.null(native)) native$beta[2] else NA_real_,
    if(!is.null(native)) native$se[2] else NA_real_),
    make_row("profile",profile_ok,if(!is.null(curve)) 100*curve$point else NA_real_,
      if(!is.null(curve)) 100*curve$se else NA_real_,
      if(!is.null(profile$fit)) profile$fit$beta[2] else NA_real_,
      if(!is.null(curve)) curve$beta_se[2] else NA_real_))
  dense <- if (i == 1L) list(native=if(!is.null(native))
    attempt(dense_check(native,cache,sim$y,x$chance)), profile=if(!is.null(profile$fit))
    attempt(dense_check(profile$fit,cache,sim$y,x$chance))) else NULL
  diagnostic <- data.frame(case=x$case,replicate=i,
    native_valid=native_ok, profile_valid=profile_ok, early_valid=isTRUE(early$converged),
    native_reported_score=if(!is.null(native)) native$covariance_score_norm else NA_real_,
    native_logsd_score=if(!is.null(native)) max(abs(score0[1:4])*
      sqrt(c(native$var_log_tau_subj,native$var_log_tau_item))) else NA_real_,
    native_rho_score=if(corr==1L) max(abs(score0[5:6])) else NA_real_,
    profile_score=if(!is.null(curve)) curve$scaled_score else NA_real_,
    objective_improvement=if(!is.null(profile)) value0-profile$value else NA_real_,
    min_eigen=if(!is.null(curve)) curve$min_eigen else NA_real_, boundary=boundary,
    optimizer_code=if(!is.null(fitopt)) fitopt$convergence else NA_integer_,
    fine_se_difference_pp=if(!is.null(fine) && !is.null(curve)) 100*(fine$se-curve$se) else NA_real_,
    alternate_objective_difference=if(!is.null(alt)) alt$value-profile$value else NA_real_,
    evaluations=evaluations, failed_evaluations=failed_evaluations,
    elapsed=proc.time()[[3]]-started, error=paste(unique(error_text),collapse="; "))
  list(rows=rows,diagnostic=diagnostic,y=as.raw(sim$y),native=native,profile=profile,
    theta=if(!is.null(fitopt)) fitopt$par else NULL,curve=curve,fine=fine,
    dense=dense,alternate=alt,early_status=early$status)
}

jobs <- expand.grid(case=cases$case, block=seq_len(ceiling(nsim/16)))
run <- function(k) {
  x <- cases[jobs$case[k],]; block <- jobs$block[k]
  path <- file.path(out,sprintf("case-%02d-block-%02d.rds",x$case,block))
  if(file.exists(path)) return(path)
  obj <- make_case(x)
  bases <- glmmplanr:::allocate_stage_point_streams(seed,1L,cases$case,nrow(cases))
  streams <- glmmplanr:::derive_rep_substreams(bases[[x$case]],nsim)
  ids <- seq.int((block-1L)*16L+1L,min(block*16L,nsim))
  details <- lapply(ids,function(i) run_trial(x,i,obj,streams[[i]]))
  atomic_rds(list(case=x,spec=obj$spec,schedule=obj$schedule,details=details),path)
  cat(sprintf("%s case %02d block %02d: %d paired datasets saved\n",mode,x$case,block,length(ids)))
  path
}
paths <- parallel::mclapply(seq_len(nrow(jobs)),run,mc.cores=workers,
  mc.set.seed=FALSE,mc.preschedule=FALSE)
stopifnot(all(vapply(paths,function(p) is.character(p) && length(p)==1L && file.exists(p),logical(1))))
details <- unlist(lapply(paths,function(p) readRDS(p)$details),recursive=FALSE)
trials <- do.call(rbind,lapply(details,`[[`,"rows"))
diagnostics <- do.call(rbind,lapply(details,`[[`,"diagnostic"))
stopifnot(nrow(trials)==nrow(cases)*nsim*2L,!anyDuplicated(trials[c("case","replicate","method")]))
write.csv(trials,gzfile(file.path(out,"fits.csv.gz")),row.names=FALSE)
write.csv(diagnostics,file.path(out,"diagnostics.csv"),row.names=FALSE)
cat("Completed",nrow(cases)*nsim,"paired datasets; native valid",sum(trials$valid[trials$method=="native"]),
  "; profile valid",sum(trials$valid[trials$method=="profile"]),"\n")
