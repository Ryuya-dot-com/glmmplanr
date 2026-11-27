// Bundled Stan model for glmmplanr: Crossed Binary GLMM with m-AFC Guessing
// Compatible with cmdstanr and rstan.
data {
  int<lower=1> N;
  array[N] int<lower=0, upper=1> y;
  int<lower=1> P;
  matrix[N, P] X;

  int<lower=1> N_subj;
  int<lower=1> Q_subj;
  array[N] int<lower=1, upper=N_subj> subj_id;
  matrix[N, Q_subj] Z_subj;

  int<lower=1> N_item;
  int<lower=1> Q_item;
  array[N] int<lower=1, upper=N_item> item_id;
  matrix[N, Q_item] Z_item;

  int<lower=0> N_testlet;
  int<lower=0> Q_testlet;
  array[N_testlet > 0 ? N : 0] int<lower=1> testlet_id;
  matrix[N_testlet > 0 ? N : 0, Q_testlet] Z_testlet;

  int<lower=1> M;
  matrix[M, P] X_cell;
  matrix[M, Q_subj] Zs_cell;
  matrix[M, Q_item] Zi_cell;
  matrix[M, Q_testlet] Zt_cell;
  vector[M] contrast;
  int<lower=1> K;
  vector[K] nodes;
  vector[K] weights;
  int<lower=0, upper=1> marginal;
  real<lower=0, upper=1> chance;
  int<lower=0, upper=1> corr_mode;
  real<lower=0> lkj_eta;
  real<lower=0> prior_beta_sd;
  real<lower=0> prior_tau_sd;
}
parameters {
  vector[P] beta;

  vector<lower=0>[Q_subj] tau_subj;
  matrix[Q_subj, N_subj] z_subj;
  cholesky_factor_corr[corr_mode == 1 ? Q_subj : 1] L_corr_subj;

  vector<lower=0>[Q_item] tau_item;
  matrix[Q_item, N_item] z_item;
  cholesky_factor_corr[corr_mode == 1 ? Q_item : 1] L_corr_item;

  vector<lower=0>[N_testlet > 0 ? Q_testlet : 0] tau_testlet;
  matrix[N_testlet > 0 ? Q_testlet : 0, N_testlet] z_testlet;
  cholesky_factor_corr[(N_testlet > 0 && corr_mode == 1) ? Q_testlet : 1] L_corr_testlet;
}
transformed parameters {
  matrix[N_subj, Q_subj] u_subj;
  matrix[N_item, Q_item] u_item;
  matrix[N_testlet, N_testlet > 0 ? Q_testlet : 0] u_testlet;

  if (corr_mode == 1 && Q_subj > 1) {
    u_subj = (diag_pre_multiply(tau_subj, L_corr_subj) * z_subj)';
  } else {
    u_subj = (diag_matrix(tau_subj) * z_subj)';
  }

  if (corr_mode == 1 && Q_item > 1) {
    u_item = (diag_pre_multiply(tau_item, L_corr_item) * z_item)';
  } else {
    u_item = (diag_matrix(tau_item) * z_item)';
  }

  if (N_testlet > 0) {
    if (corr_mode == 1 && Q_testlet > 1) {
      u_testlet = (diag_pre_multiply(tau_testlet, L_corr_testlet) * z_testlet)';
    } else {
      u_testlet = (diag_matrix(tau_testlet) * z_testlet)';
    }
  }
}
model {
  // Weakly informative priors preventing boundary collapse
  beta[1] ~ normal(0, 1.5 * prior_beta_sd);
  if (P > 1) {
    beta[2:P] ~ normal(0, prior_beta_sd);
  }

  tau_subj ~ student_t(3, 0, prior_tau_sd);
  to_vector(z_subj) ~ std_normal();
  if (corr_mode == 1 && Q_subj > 1) {
    L_corr_subj ~ lkj_corr_cholesky(lkj_eta);
  }

  tau_item ~ student_t(3, 0, prior_tau_sd);
  to_vector(z_item) ~ std_normal();
  if (corr_mode == 1 && Q_item > 1) {
    L_corr_item ~ lkj_corr_cholesky(lkj_eta);
  }

  if (N_testlet > 0) {
    tau_testlet ~ student_t(3, 0, prior_tau_sd);
    to_vector(z_testlet) ~ std_normal();
    if (corr_mode == 1 && Q_testlet > 1) {
      L_corr_testlet ~ lkj_corr_cholesky(lkj_eta);
    }
  }

  vector[N] eta = X * beta +
    rows_dot_product(Z_subj, u_subj[subj_id]) +
    rows_dot_product(Z_item, u_item[item_id]);

  if (N_testlet > 0) {
    eta += rows_dot_product(Z_testlet, u_testlet[testlet_id]);
  }

  if (chance == 0) {
    y ~ bernoulli_logit(eta);
  } else {
    vector[N] p_cond = inv_logit(eta);
    vector[N] mu = chance + (1.0 - chance) * p_cond;
    y ~ bernoulli(mu);
  }
}

generated quantities {
  real prob_contrast = 0;
  for (m in 1:M) {
    real mu = X_cell[m] * beta;
    real v = 0;
    if (marginal == 1) {
      row_vector[Q_subj] zs = Zs_cell[m] .* tau_subj';
      row_vector[Q_item] zi = Zi_cell[m] .* tau_item';
      if (corr_mode == 1) {
        zs = zs * L_corr_subj;
        zi = zi * L_corr_item;
      }
      v = dot_self(zs) + dot_self(zi);
      if (N_testlet > 0) {
        row_vector[Q_testlet] zt = Zt_cell[m] .* tau_testlet';
        if (corr_mode == 1) zt = zt * L_corr_testlet;
        v += dot_self(zt);
      }
    }
    prob_contrast += contrast[m] * (chance + (1 - chance) *
      dot_product(weights, inv_logit(mu + sqrt(v) * nodes)));
  }
}
