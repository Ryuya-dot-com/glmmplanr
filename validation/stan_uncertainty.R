# Longer reference fit for the same single dataset as test-stan-reference.R.
library(glmmplanr)
source("tests/testthat/helper-stan.R")
dgp <- three_choice_reference_dgp()
spec <- specify_target(dgp, "exposure", engine = "bayes_laplace", on_singular = "flag_only")
des <- dgp$design
set.seed(303, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion", sample.kind = "Rejection")
tr <- simulate_dgp_data(dgp, generate_design_schedule(des, 48L, 32L))
cpp <- glmmplanr:::fit_single_dataset(tr, spec, engine_use = "bayes_laplace")
stopifnot(cpp$converged, spec$analysis_re == "maximal", !des$has_testlet)
dat <- three_choice_stan_data(tr, spec)
model <- rstan::stan_model(system.file("stan", "crossed_binomial_mafc.stan", package = "glmmplanr"))
fit <- rstan::sampling(model, data = dat, chains = 4L, cores = 2L, iter = 1500L, warmup = 750L,
                       seed = 913003L, refresh = 0L, pars = "beta")
idx <- match(spec$effect_col, colnames(des$X_cells))
summary_beta <- rstan::summary(fit, pars = sprintf("beta[%d]", idx))$summary
divergences <- sum(vapply(rstan::get_sampler_params(fit, inc_warmup = FALSE),
                          function(chain) sum(chain[, "divergent__"]), numeric(1)))
print(data.frame(cpp_est = cpp$est, cpp_se = cpp$se, stan_mean = summary_beta[1, "mean"],
                 stan_sd = summary_beta[1, "sd"], n_eff = summary_beta[1, "n_eff"],
                 rhat = summary_beta[1, "Rhat"], divergences = divergences))
