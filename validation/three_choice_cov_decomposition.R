# Exploratory variance-decomposition check on the 40 paired 48x32 3AFC datasets
# (seeds 7312001-7312040): compares fixed true covariance vs. estimated covariance
# across C++ (conditional Schur SE vs. joint Laplace SE) and RStan HMC.
library(glmmplanr)
source("tests/testthat/helper-stan.R")
stopifnot(requireNamespace("posterior", quietly = TRUE),
          requireNamespace("rstan", quietly = TRUE),
          requireNamespace("digest", quietly = TRUE))

dgp <- three_choice_reference_dgp()
spec <- specify_target(dgp, "exposure", engine = "bayes_laplace",
                       control = target_control(on_singular = "flag_only"))
schedule <- generate_design_schedule(dgp$design, n_subj = 48L, n_item = 32L)
eff_idx <- match(spec$effect_col, colnames(dgp$design$X_cells))

native_prev <- subset(read.csv("validation/three_choice_precision_native.csv"),
                      persons == 48L & items == 32L & trial <= 40L)
stan_prev <- subset(read.csv("validation/three_choice_precision_stan.csv"),
                    seed %in% native_prev$seed)
stopifnot(nrow(native_prev) == 40L, nrow(stan_prev) == 40L,
          identical(native_prev$seed, stan_prev$seed))

fixed_stan_code <- "
data {
  int<lower=1> N;
  array[N] int<lower=0, upper=1> y;
  int<lower=1> P;
  matrix[N, P] X;
  int<lower=1> N_subj;
  int<lower=1> Q_subj;
  array[N] int<lower=1, upper=N_subj> subj_id;
  matrix[N, Q_subj] Z_subj;
  matrix<lower=-10, upper=10>[Q_subj, Q_subj] L_subj_true;
  int<lower=1> N_item;
  int<lower=1> Q_item;
  array[N] int<lower=1, upper=N_item> item_id;
  matrix[N, Q_item] Z_item;
  matrix<lower=-10, upper=10>[Q_item, Q_item] L_item_true;
  real<lower=0, upper=0.99> chance;
  real<lower=0.1> prior_beta_sd;
}
parameters {
  vector[P] beta;
  matrix[Q_subj, N_subj] z_subj;
  matrix[Q_item, N_item] z_item;
}
model {
  matrix[N_subj, Q_subj] u_subj = (L_subj_true * z_subj)';
  matrix[N_item, Q_item] u_item = (L_item_true * z_item)';
  vector[N] eta = X * beta +
    rows_dot_product(Z_subj, u_subj[subj_id]) +
    rows_dot_product(Z_item, u_item[item_id]);
  beta[1] ~ normal(0, 1.5 * prior_beta_sd);
  if (P > 1) beta[2:P] ~ normal(0, prior_beta_sd);
  to_vector(z_subj) ~ std_normal();
  to_vector(z_item) ~ std_normal();
  y ~ bernoulli(chance + (1.0 - chance) * inv_logit(eta));
}
"
fixed_mod <- rstan::stan_model(model_code = fixed_stan_code)

run_pair <- function(idx_row) {
  sd_val <- native_prev$seed[idx_row]
  set.seed(sd_val, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion", sample.kind = "Rejection")
  tr <- simulate_dgp_data(dgp, schedule)
  stopifnot(identical(digest::digest(tr$y, algo = "xxhash64"), native_prev$response_hash[idx_row]))

  X <- dgp$design$X_cells[tr$cell_id, , drop = FALSE]
  Z_s <- dgp$design$Z_subj_cells[tr$cell_id, , drop = FALSE]
  Z_i <- dgp$design$Z_item_cells[tr$cell_id, , drop = FALSE]

  f_fix <- glmmplanr:::cpp_fit_crossed_glmm(
    y = as.integer(tr$y), X = X,
    subj_idx = tr$subj_id - 1L, Z_subj = Z_s, n_subj = 48L,
    item_idx = tr$item_id - 1L, Z_item = Z_i, n_item = 32L,
    testlet_idx = integer(0), Z_testlet = matrix(0, nrow(tr), 0L), n_testlet = 0L,
    chance = dgp$chance, corr_mode = 1L, regularized = TRUE,
    lkj_eta = spec$lkj_eta, prior_beta_sd = spec$prior_beta_sd, prior_tau_sd = spec$prior_tau_sd,
    init_Sigma_subj = dgp$Sigma_subj, init_Sigma_item = dgp$Sigma_item,
    fix_cov = TRUE, compute_joint_vcov = TRUE
  )

  f_est <- glmmplanr:::cpp_fit_crossed_glmm(
    y = as.integer(tr$y), X = X,
    subj_idx = tr$subj_id - 1L, Z_subj = Z_s, n_subj = 48L,
    item_idx = tr$item_id - 1L, Z_item = Z_i, n_item = 32L,
    testlet_idx = integer(0), Z_testlet = matrix(0, nrow(tr), 0L), n_testlet = 0L,
    chance = dgp$chance, corr_mode = 1L, regularized = TRUE,
    lkj_eta = spec$lkj_eta, prior_beta_sd = spec$prior_beta_sd, prior_tau_sd = spec$prior_tau_sd,
    max_em_iter = 400L, max_pirls_iter = 44L,
    fix_cov = FALSE, compute_joint_vcov = TRUE
  )

  sdat <- list(
    N = nrow(tr), y = as.integer(tr$y), P = ncol(X), X = X,
    N_subj = 48L, Q_subj = ncol(Z_s), subj_id = tr$subj_id, Z_subj = Z_s,
    L_subj_true = dgp$L_subj,
    N_item = 32L, Q_item = ncol(Z_i), item_id = tr$item_id, Z_item = Z_i,
    L_item_true = dgp$L_item,
    chance = dgp$chance, prior_beta_sd = spec$prior_beta_sd
  )
  sfit <- suppressWarnings(rstan::sampling(
    fixed_mod, data = sdat, seed = sd_val + 300000L,
    chains = 2L, iter = 800L, warmup = 400L, cores = 1L, refresh = 0L, pars = "beta"
  ))
  par_nm <- sprintf("beta[%d]", eff_idx)
  draws <- rstan::extract(sfit, permuted = FALSE, inc_warmup = FALSE, pars = par_nm)[, , par_nm]
  sam <- do.call(rbind, rstan::get_sampler_params(sfit, inc_warmup = FALSE))
  s_valid <- posterior::rhat(draws) <= 1.01 &&
    posterior::ess_bulk(draws) >= 100 &&
    posterior::ess_tail(draws) >= 100 &&
    sum(sam[, "divergent__"]) == 0L

  data.frame(
    seed = sd_val,
    truth = spec$true_beta,
    cpp_fix_valid = isTRUE(f_fix$converged),
    cpp_fix_est = f_fix$beta[eff_idx],
    cpp_fix_schur_se = sqrt(f_fix$vcov_cond[eff_idx, eff_idx]),
    cpp_fix_lap_se = f_fix$se[eff_idx],
    stan_fix_valid = s_valid,
    stan_fix_est = mean(draws),
    stan_fix_se = sd(draws),
    cpp_est_valid = isTRUE(f_est$converged),
    cpp_est = f_est$beta[eff_idx],
    cpp_est_schur_se = sqrt(f_est$vcov_cond[eff_idx, eff_idx]),
    cpp_est_joint_se = f_est$se[eff_idx],
    stan_est_valid = stan_prev$stan_valid[idx_row],
    stan_est = stan_prev$stan_est[idx_row],
    stan_est_se = stan_prev$stan_se[idx_row]
  )
}

workers <- if (.Platform$OS.type == "windows") 1L else 8L
rows <- parallel::mclapply(seq_len(40L), run_pair, mc.cores = workers, mc.set.seed = FALSE)
stopifnot(all(vapply(rows, is.data.frame, logical(1))))
paired <- do.call(rbind, rows)

truth <- spec$true_beta
sum_method <- function(label, valid, est, se) {
  ok <- is.finite(est) & is.finite(se) & valid
  e <- est[ok]; s <- se[ok]
  cov_k <- sum(abs(e - truth) <= 1.96 * s)
  nar_k <- sum(2 * 1.96 * s <= 1.0)
  data.frame(
    method = label,
    attempted = length(valid),
    valid = sum(ok),
    mean_est = mean(e),
    bias = mean(e - truth),
    empirical_sd = sd(e),
    mean_se = mean(s),
    mean_width = mean(2 * 1.96 * s),
    covered = cov_k,
    narrow = nar_k
  )
}

decomp_summary <- rbind(
  sum_method("cpp_fixed_true_schur", paired$cpp_fix_valid, paired$cpp_fix_est, paired$cpp_fix_schur_se),
  sum_method("cpp_fixed_true_laplace", paired$cpp_fix_valid, paired$cpp_fix_est, paired$cpp_fix_lap_se),
  sum_method("stan_fixed_true_hmc", paired$stan_fix_valid, paired$stan_fix_est, paired$stan_fix_se),
  sum_method("cpp_estimated_schur", paired$cpp_est_valid, paired$cpp_est, paired$cpp_est_schur_se),
  sum_method("cpp_estimated_joint", paired$cpp_est_valid, paired$cpp_est, paired$cpp_est_joint_se),
  sum_method("stan_estimated_hmc", paired$stan_est_valid, paired$stan_est, paired$stan_est_se)
)

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[1L] else "validation"
write.csv(paired, file.path(out_dir, "three_choice_cov_decomp_paired.csv"), row.names = FALSE)
write.csv(decomp_summary, file.path(out_dir, "three_choice_cov_decomp_summary.csv"), row.names = FALSE)
print(decomp_summary, row.names = FALSE, digits = 4)
