#include <RcppEigen.h>
#include <Rmath.h>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

// [[Rcpp::depends(RcppEigen)]]

namespace {

using RowMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

inline double inv_logit(double x) {
  if (x >= 0.0) {
    double e = std::exp(-x);
    return 1.0 / (1.0 + e);
  } else {
    double e = std::exp(x);
    return e / (1.0 + e);
  }
}

Eigen::MatrixXd safe_chol(const Eigen::MatrixXd& S, double min_diag = 1e-10) {
  int q = S.rows();
  if (q == 0) return Eigen::MatrixXd(0, 0);
  if (q == 1) {
    Eigen::MatrixXd L(1, 1);
    L(0, 0) = std::sqrt(std::max(S(0, 0), min_diag));
    return L;
  }
  Eigen::MatrixXd M = 0.5 * (S + S.transpose());
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(M);
  Eigen::VectorXd ev = es.eigenvalues().cwiseMax(min_diag);
  Eigen::MatrixXd S_pos = es.eigenvectors() * ev.asDiagonal() * es.eigenvectors().transpose();
  Eigen::LLT<Eigen::MatrixXd> llt(S_pos);
  if (llt.info() == Eigen::Success) {
    return llt.matrixL();
  }
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(q, q);
  return (S_pos + 1e-7 * I).llt().matrixL();
}

Eigen::MatrixXd scale_cov_single_sd(const Eigen::MatrixXd& Sigma, int j, double mult) {
  Eigen::MatrixXd Out = Sigma;
  Out.row(j) *= mult;
  Out.col(j) *= mult;
  return Out;
}

Eigen::MatrixXd update_cov_matrix(
    const Eigen::MatrixXd& sum_outer,
    const Eigen::MatrixXd& Sigma_prev,
    int n_clusters,
    int corr_mode,
    bool regularized,
    double lkj_eta,
    double prior_tau_sd
) {
  int q = sum_outer.rows();
  Eigen::MatrixXd S = 0.5 * (sum_outer + sum_outer.transpose());

  if (!regularized) {
    Eigen::MatrixXd Sigma = S / static_cast<double>(n_clusters);
    for (int j = 0; j < q; ++j) {
      if (Sigma(j, j) < 2.5e-5 && Sigma(j, j) <= Sigma_prev(j, j)) {
        Sigma.row(j).setZero();
        Sigma.col(j).setZero();
        Sigma(j, j) = 1e-10;
      }
    }
    if (corr_mode == 0) {
      Eigen::VectorXd d = Sigma.diagonal().cwiseMax(1e-10);
      return d.asDiagonal();
    }
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Sigma);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_prev(0.5 * (Sigma_prev + Sigma_prev.transpose()));
    Eigen::VectorXd ev = es.eigenvalues();
    Eigen::VectorXd ev_prev = es_prev.eigenvalues();
    for (int j = 0; j < q; ++j) {
      if (ev(j) < 5e-5 && ev(j) <= ev_prev(j)) {
        ev(j) = 1e-10;
      } else if (ev(j) < 1e-10) {
        ev(j) = 1e-10;
      }
    }
    return es.eigenvectors() * ev.asDiagonal() * es.eigenvectors().transpose();
  }

  double s0_sq = std::max(prior_tau_sd * prior_tau_sd, 0.25);
  Eigen::VectorXd new_sd(q);
  for (int j = 0; j < q; ++j) {
    double a_exp = (j == 0) ? 1.0 : 0.75;
    double c_lin = static_cast<double>(n_clusters) - a_exp;
    double s_jj = std::max(S(j, j), 1e-8);
    double inv_s0 = 1.0 / s0_sq;
    double disc = c_lin * c_lin + 4.0 * inv_s0 * s_jj;
    double v_hat = (-c_lin + std::sqrt(std::max(disc, 0.0))) / (2.0 * inv_s0);
    new_sd(j) = std::sqrt(std::max(v_hat, 1e-6));
  }

  if (corr_mode == 0 || q == 1) {
    return (new_sd.array().square()).matrix().asDiagonal();
  }

  Eigen::VectorXd emp_sd = S.diagonal().cwiseMax(1e-10).cwiseSqrt();
  Eigen::MatrixXd R = Eigen::MatrixXd::Identity(q, q);
  double shrink = static_cast<double>(n_clusters) /
    (static_cast<double>(n_clusters) + 2.0 * std::max(lkj_eta - 1.0, 0.0) + static_cast<double>(q));
  shrink = std::min(std::max(shrink, 0.0), 0.995);

  for (int r = 0; r < q; ++r) {
    for (int c = 0; c < r; ++c) {
      double cor_rc = shrink * (S(r, c) / (emp_sd(r) * emp_sd(c)));
      cor_rc = std::min(std::max(cor_rc, -0.98), 0.98);
      R(r, c) = cor_rc;
      R(c, r) = cor_rc;
    }
  }
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_r(R);
  if (es_r.eigenvalues().minCoeff() < 1e-4) {
    Eigen::VectorXd ev = es_r.eigenvalues().cwiseMax(1e-3);
    Eigen::MatrixXd R_fix = es_r.eigenvectors() * ev.asDiagonal() * es_r.eigenvectors().transpose();
    Eigen::VectorXd d_inv = R_fix.diagonal().cwiseSqrt().cwiseInverse();
    R = d_inv.asDiagonal() * R_fix * d_inv.asDiagonal();
  }

  return new_sd.asDiagonal() * R * new_sd.asDiagonal();
}

bool is_cov_singular(const Eigen::MatrixXd& Sigma, int corr_mode, double tol_sd = 1e-2) {
  int q = Sigma.rows();
  Eigen::VectorXd sd = Sigma.diagonal().cwiseMax(0.0).cwiseSqrt();
  for (int j = 0; j < q; ++j) {
    if (sd(j) < tol_sd) return true;
  }
  if (corr_mode == 1 && q > 1) {
    Eigen::MatrixXd R = Eigen::MatrixXd::Identity(q, q);
    for (int r = 0; r < q; ++r) {
      for (int c = 0; c < r; ++c) {
        double cor = Sigma(r, c) / std::max(sd(r) * sd(c), 1e-12);
        R(r, c) = cor;
        R(c, r) = cor;
        if (std::abs(cor) > 0.99) return true;
      }
    }
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(R);
    if (es.eigenvalues().minCoeff() < 1e-3) return true;
  }
  return false;
}

Eigen::MatrixXd squarem_extrapolate_cov(
    const Eigen::MatrixXd& S0,
    const Eigen::MatrixXd& S1,
    const Eigen::MatrixXd& S2,
    double alpha_sq,
    int corr_mode
) {
  int q = S0.rows();
  Eigen::VectorXd sd0 = S0.diagonal().cwiseMax(0.0).cwiseSqrt();
  Eigen::VectorXd sd1 = S1.diagonal().cwiseMax(0.0).cwiseSqrt();
  Eigen::VectorXd sd2 = S2.diagonal().cwiseMax(0.0).cwiseSqrt();

  Eigen::VectorXd r_sd = sd1 - sd0;
  Eigen::VectorXd v_sd = (sd2 - sd1) - r_sd;
  Eigen::VectorXd sd_prop = sd0 - 2.0 * alpha_sq * r_sd + (alpha_sq * alpha_sq) * v_sd;
  for (int j = 0; j < q; ++j) {
    sd_prop(j) = std::min(std::max(sd_prop(j), 1e-5), 4.5);
  }

  if (corr_mode == 0 || q == 1) {
    return (sd_prop.array().square()).matrix().asDiagonal();
  }

  auto extract_corr_chol = [&](const Eigen::MatrixXd& S, const Eigen::VectorXd& sd) {
    Eigen::MatrixXd R = Eigen::MatrixXd::Identity(q, q);
    for (int r = 0; r < q; ++r) {
      if (sd(r) <= 1e-4) continue;
      for (int c = 0; c < r; ++c) {
        if (sd(c) <= 1e-4) continue;
        double cor = S(r, c) / std::max(sd(r) * sd(c), 1e-12);
        cor = std::min(std::max(cor, -0.995), 0.995);
        R(r, c) = cor;
        R(c, r) = cor;
      }
    }
    return safe_chol(R);
  };

  Eigen::MatrixXd C0 = extract_corr_chol(S0, sd0);
  Eigen::MatrixXd C1 = extract_corr_chol(S1, sd1);
  Eigen::MatrixXd C2 = extract_corr_chol(S2, sd2);
  Eigen::MatrixXd C_prop = C0 - 2.0 * alpha_sq * (C1 - C0) +
    (alpha_sq * alpha_sq) * ((C2 - C1) - (C1 - C0));

  for (int r = 0; r < q; ++r) {
    if (sd_prop(r) <= 2e-5) {
      C_prop.row(r).setZero();
      C_prop(r, r) = 1.0;
      continue;
    }
    for (int c = 0; c < q; ++c) {
      if (c > r) {
        C_prop(r, c) = 0.0;
      } else if (c == r) {
        C_prop(r, c) = std::max(C_prop(r, c), 2e-2);
      } else if (sd_prop(c) <= 2e-5) {
        C_prop(r, c) = 0.0;
      }
    }
    double rn = C_prop.row(r).norm();
    if (rn > 1e-12) C_prop.row(r) /= rn;
  }

  Eigen::MatrixXd R_prop = C_prop * C_prop.transpose();
  return sd_prop.asDiagonal() * R_prop * sd_prop.asDiagonal();
}

struct FitResult {
  Eigen::VectorXd beta;
  Eigen::VectorXd se_beta;
  Eigen::MatrixXd vcov_beta;
  Eigen::MatrixXd vcov_cond;
  Eigen::MatrixXd Sigma_s;
  Eigen::MatrixXd Sigma_i;
  Eigen::MatrixXd Sigma_t;
  Eigen::VectorXd var_log_tau_s;
  Eigen::VectorXd var_log_tau_i;
  Eigen::VectorXd var_log_tau_t;
  Eigen::MatrixXd J_beta_theta;
  Eigen::MatrixXd u_subj;
  Eigen::MatrixXd u_item;
  Eigen::MatrixXd u_testlet;
  double log_lik_laplace = NA_REAL;
  double final_outer_diff = 1.0;
  double final_inner_step = 1.0;
  int outer_used = 0;
  bool singular = false;
  bool strictly_converged = false;
  std::string status = "engine_error";
};

struct GLMMWorkspace {
  int N, P, q_s, q_i, q_t;
  int n_subj, n_item, n_testlet;
  bool has_testlet;
  int dim_s, dim_i, dim_t, dim_2;

  Eigen::VectorXi subj_idx, item_idx, testlet_idx;
  std::vector<std::vector<int>> subj_trials;

  RowMatrixXd X_rm, Zs_rm, Zi_rm, Zt_rm;
  RowMatrixXd ZsL, ZiL, ZtL;

  Eigen::VectorXd w_vec, wp_vec, g_vec;
  std::vector<Eigen::MatrixXd> A_ss_inv;
  Eigen::MatrixXd L_22_saved, S_22, Bt, Bt_Ass_inv, C_bs, C_bs_Ass_inv, S_2b_saved;
  Eigen::MatrixXd H_bb_schur, L22_inv, W_mat;
  Eigen::VectorXd c_lap, q_tmp, cs_vec, c2_vec;

  GLMMWorkspace(
      const Eigen::MatrixXd& X,
      const Eigen::VectorXi& s_idx,
      const Eigen::MatrixXd& Z_subj,
      int ns,
      const Eigen::VectorXi& i_idx,
      const Eigen::MatrixXd& Z_item,
      int ni,
      const Eigen::VectorXi& t_idx,
      const Eigen::MatrixXd& Z_testlet,
      int nt
  ) : N(X.rows()), P(X.cols()),
      q_s(Z_subj.cols()), q_i(Z_item.cols()),
      n_subj(ns), n_item(ni), n_testlet(nt),
      has_testlet(nt > 0 && Z_testlet.cols() > 0),
      subj_idx(s_idx), item_idx(i_idx), testlet_idx(t_idx),
      subj_trials(ns),
      X_rm(X), Zs_rm(Z_subj), Zi_rm(Z_item),
      Zt_rm(has_testlet ? RowMatrixXd(Z_testlet) : RowMatrixXd(0, 0))
  {
    q_t = has_testlet ? Z_testlet.cols() : 0;
    dim_s = n_subj * q_s;
    dim_i = n_item * q_i;
    dim_t = has_testlet ? (n_testlet * q_t) : 0;
    dim_2 = dim_i + dim_t;

    for (int n = 0; n < N; ++n) {
      subj_trials[subj_idx(n)].push_back(n);
    }

    ZsL.resize(N, q_s);
    ZiL.resize(N, q_i);
    ZtL.resize(has_testlet ? N : 0, has_testlet ? q_t : 0);

    w_vec.resize(N);
    wp_vec.resize(N);
    g_vec.resize(N);

    A_ss_inv.assign(n_subj, Eigen::MatrixXd::Identity(q_s, q_s));
    L_22_saved = Eigen::MatrixXd::Identity(dim_2, dim_2);
    S_22.resize(dim_2, dim_2);
    Bt.resize(dim_2, dim_s);
    Bt_Ass_inv.resize(dim_2, dim_s);
    C_bs.resize(P, dim_s);
    C_bs_Ass_inv.resize(P, dim_s);
    S_2b_saved.resize(dim_2, P);
    H_bb_schur = Eigen::MatrixXd::Identity(P, P);
    L22_inv.resize(dim_2, dim_2);
    W_mat.resize(dim_2, dim_s);
    c_lap.resize(N);
    q_tmp.resize(dim_2);
    cs_vec.resize(dim_s);
    c2_vec.resize(dim_2);
  }

  double eval_data_ll_only(
      const Eigen::VectorXi& y,
      const Eigen::VectorXd& beta,
      const Eigen::VectorXd& v_s,
      const Eigen::VectorXd& v_2,
      double chance
  ) const {
    double ll_data = 0.0;
    for (int n = 0; n < N; ++n) {
      int s = subj_idx(n);
      int i = item_idx(n);
      double eta = X_rm.row(n).dot(beta) +
        ZsL.row(n).dot(v_s.segment(s * q_s, q_s)) +
        ZiL.row(n).dot(v_2.segment(i * q_i, q_i));
      if (has_testlet) {
        int t = testlet_idx(n);
        eta += ZtL.row(n).dot(v_2.segment(dim_i + t * q_t, q_t));
      }
      eta = std::min(std::max(eta, -25.0), 25.0);
      double p = inv_logit(eta);
      double mu = chance + (1.0 - chance) * p;
      mu = std::min(std::max(mu, chance + 1e-12), 1.0 - 1e-12);
      ll_data += (y(n) == 1) ? std::log(mu) : std::log1p(-mu);
    }
    return ll_data;
  }

  double eval_trials(
      const Eigen::VectorXi& y,
      const Eigen::VectorXd& beta,
      const Eigen::VectorXd& v_s,
      const Eigen::VectorXd& v_2,
      double chance,
      bool observed
  ) {
    double ll_data = 0.0;
    for (int n = 0; n < N; ++n) {
      int s = subj_idx(n);
      int i = item_idx(n);
      double eta = X_rm.row(n).dot(beta) +
        ZsL.row(n).dot(v_s.segment(s * q_s, q_s)) +
        ZiL.row(n).dot(v_2.segment(i * q_i, q_i));
      if (has_testlet) {
        int t = testlet_idx(n);
        eta += ZtL.row(n).dot(v_2.segment(dim_i + t * q_t, q_t));
      }
      eta = std::min(std::max(eta, -25.0), 25.0);

      double p = inv_logit(eta);
      double omp = 1.0 - p;
      double mu = chance + (1.0 - chance) * p;
      mu = std::min(std::max(mu, chance + 1e-12), 1.0 - 1e-12);

      const bool correct = (y(n) == 1);
      ll_data += correct ? std::log(mu) : std::log1p(-mu);
      g_vec(n) = correct ? ((1.0 - chance) * p / mu) * omp : -p;
      if (observed) {
        // Exact observed curvature and third derivative of log P(Y = y_n | eta_n)
        const double q = correct ? p / mu : 0.0;
        const double v = q * (1.0 - q);
        w_vec(n) = p * omp - v;
        wp_vec(n) = p * omp * (1.0 - 2.0 * p) - v * (1.0 - 2.0 * q);
      } else {
        // Safeguarded concave curvature:
        // For y = 0, -d^2/d_eta^2 log(1 - mu) = p*(1 - p) > 0 everywhere even when chance > 0!
        // For y = 1, use max(observed curvature, Fisher scoring curvature) so search direction
        // stays strictly positive-definite while preserving exact curvature in concave regions.
        if (chance <= 0.0 || !correct) {
          w_vec(n) = std::max(p * omp, 1e-9);
        } else {
          const double q = p / mu;
          const double w_obs = p * omp - q * (1.0 - q);
          const double w_fs = ((1.0 - chance) * p * p * omp) / mu;
          w_vec(n) = std::max({w_obs, w_fs, 1e-9});
        }
      }
    }
    return ll_data;
  }

  bool assemble_and_factor_core(
      const Eigen::VectorXd& beta,
      const Eigen::VectorXd& v_s,
      const Eigen::VectorXd& v_2,
      const Eigen::MatrixXd& P_beta,
      Eigen::VectorXd& rhs_s,
      Eigen::VectorXd& rhs_2,
      Eigen::VectorXd& rhs_b,
      Eigen::LLT<Eigen::MatrixXd>& llt_22,
      double& logdet_out
  ) {
    S_22.setIdentity();
    Bt.setZero();
    C_bs.setZero();
    Eigen::MatrixXd C_2b = Eigen::MatrixXd::Zero(dim_2, P);
    Eigen::MatrixXd H_bb = P_beta;

    rhs_s = -v_s;
    rhs_2 = -v_2;
    rhs_b.noalias() = X_rm.transpose() * g_vec - P_beta * beta;

    double logdet_Ass = 0.0;

    for (int s = 0; s < n_subj; ++s) {
      Eigen::MatrixXd Ass = Eigen::MatrixXd::Identity(q_s, q_s);
      const int off_s = s * q_s;
      for (int n : subj_trials[s]) {
        const double wn = w_vec[n];
        const double gn = g_vec[n];
        const double* zs = ZsL.data() + static_cast<std::ptrdiff_t>(n) * q_s;
        const double* zi = ZiL.data() + static_cast<std::ptrdiff_t>(n) * q_i;
        const double* xn = X_rm.data() + static_cast<std::ptrdiff_t>(n) * P;
        const int off_i = item_idx[n] * q_i;

        for (int a = 0; a < q_s; ++a) {
          const double w_zsa = wn * zs[a];
          rhs_s[off_s + a] += gn * zs[a];
          for (int b = 0; b < q_s; ++b) {
            Ass(b, a) += w_zsa * zs[b];
          }
          for (int b = 0; b < q_i; ++b) {
            Bt(off_i + b, off_s + a) += w_zsa * zi[b];
          }
          for (int p = 0; p < P; ++p) {
            C_bs(p, off_s + a) += w_zsa * xn[p];
          }
        }

        for (int a = 0; a < q_i; ++a) {
          const double w_zia = wn * zi[a];
          rhs_2[off_i + a] += gn * zi[a];
          for (int b = 0; b < q_i; ++b) {
            S_22(off_i + b, off_i + a) += w_zia * zi[b];
          }
          for (int p = 0; p < P; ++p) {
            C_2b(off_i + a, p) += w_zia * xn[p];
          }
        }

        for (int p1 = 0; p1 < P; ++p1) {
          const double w_xp1 = wn * xn[p1];
          for (int p2 = 0; p2 < P; ++p2) {
            H_bb(p2, p1) += w_xp1 * xn[p2];
          }
        }

        if (has_testlet) {
          const double* zt = ZtL.data() + static_cast<std::ptrdiff_t>(n) * q_t;
          const int off_t = dim_i + testlet_idx[n] * q_t;
          for (int a = 0; a < q_t; ++a) {
            const double w_zta = wn * zt[a];
            rhs_2[off_t + a] += gn * zt[a];
            for (int b = 0; b < q_t; ++b) {
              S_22(off_t + b, off_t + a) += w_zta * zt[b];
            }
            for (int b = 0; b < q_i; ++b) {
              const double v_it = w_zta * zi[b];
              S_22(off_i + b, off_t + a) += v_it;
              S_22(off_t + a, off_i + b) += v_it;
            }
            for (int b = 0; b < q_s; ++b) {
              Bt(off_t + a, off_s + b) += w_zta * zs[b];
            }
            for (int p = 0; p < P; ++p) {
              C_2b(off_t + a, p) += w_zta * xn[p];
            }
          }
        }
      }

      Eigen::LLT<Eigen::MatrixXd> llt_s(Ass);
      if (llt_s.info() != Eigen::Success) return false;
      const Eigen::MatrixXd& L_ass = llt_s.matrixL();
      for (int j = 0; j < q_s; ++j) logdet_Ass += 2.0 * std::log(L_ass(j, j));
      A_ss_inv[s] = llt_s.solve(Eigen::MatrixXd::Identity(q_s, q_s));

      Bt_Ass_inv.middleCols(off_s, q_s).noalias() = Bt.middleCols(off_s, q_s) * A_ss_inv[s];
      C_bs_Ass_inv.middleCols(off_s, q_s).noalias() = C_bs.middleCols(off_s, q_s) * A_ss_inv[s];
    }

    S_22.noalias() -= Bt_Ass_inv * Bt.transpose();
    S_2b_saved.noalias() = C_2b - Bt_Ass_inv * C_bs.transpose();

    llt_22.compute(S_22);
    if (llt_22.info() != Eigen::Success) return false;
    L_22_saved = llt_22.matrixL();
    double logdet_S22 = 0.0;
    for (int j = 0; j < dim_2; ++j) logdet_S22 += 2.0 * std::log(L_22_saved(j, j));

    Eigen::MatrixXd S22_inv_S2b = llt_22.solve(S_2b_saved);
    H_bb_schur.noalias() = H_bb - C_bs_Ass_inv * C_bs.transpose() - S_2b_saved.transpose() * S22_inv_S2b;

    logdet_out = logdet_Ass + logdet_S22;
    return true;
  }

  double assemble_and_factor(
      const Eigen::VectorXd& beta,
      const Eigen::VectorXd& v_s,
      const Eigen::VectorXd& v_2,
      const Eigen::MatrixXd& P_beta,
      Eigen::VectorXd& rhs_s,
      Eigen::VectorXd& rhs_2,
      Eigen::VectorXd& rhs_b,
      Eigen::LLT<Eigen::MatrixXd>& llt_22,
      bool allow_clamp = false
  ) {
    double logdet = 0.0;
    if (assemble_and_factor_core(beta, v_s, v_2, P_beta, rhs_s, rhs_2, rhs_b, llt_22, logdet)) {
      return logdet;
    }
    if (allow_clamp) {
      for (int n = 0; n < N; ++n) {
        if (w_vec(n) < 1e-6) w_vec(n) = 1e-6;
      }
      if (assemble_and_factor_core(beta, v_s, v_2, P_beta, rhs_s, rhs_2, rhs_b, llt_22, logdet)) {
        return logdet;
      }
    }
    Rcpp::stop("non_pd_hessian: random-effects curvature block");
  }

  Eigen::VectorXd compute_laplace_grad_beta(
      const Eigen::VectorXi& y,
      const Eigen::VectorXd& beta_eval,
      const Eigen::MatrixXd& L_s_eval,
      const Eigen::MatrixXd& L_i_eval,
      const Eigen::MatrixXd& L_t_eval,
      const Eigen::MatrixXd& P_beta,
      double chance,
      Eigen::VectorXd v_s_cur,
      Eigen::VectorXd v_2_cur
  ) {
    ZsL.noalias() = Zs_rm * L_s_eval;
    ZiL.noalias() = Zi_rm * L_i_eval;
    if (has_testlet) ZtL.noalias() = Zt_rm * L_t_eval;

    Eigen::VectorXd rhs_s(dim_s), rhs_2(dim_2), rhs_b(P);
    Eigen::LLT<Eigen::MatrixXd> llt_22;

    // Solve conditional mode (v_s, v_2) at fixed (beta_eval, L_*_eval)
    for (int it = 0; it < 6; ++it) {
      double ll_cur = eval_trials(y, beta_eval, v_s_cur, v_2_cur, chance, false);
      double phi_cur = ll_cur - 0.5 * (v_s_cur.squaredNorm() + v_2_cur.squaredNorm());
      assemble_and_factor(beta_eval, v_s_cur, v_2_cur, P_beta, rhs_s, rhs_2, rhs_b, llt_22, true);
      Eigen::VectorXd a(dim_s);
      for (int s = 0; s < n_subj; ++s) {
        a.segment(s * q_s, q_s).noalias() = A_ss_inv[s] * rhs_s.segment(s * q_s, q_s);
      }
      Eigen::VectorXd step_2 = llt_22.solve(rhs_2 - Bt * a);
      Eigen::VectorXd step_s = a - Bt_Ass_inv.transpose() * step_2;
      double st_max = std::max(step_s.lpNorm<Eigen::Infinity>(), step_2.lpNorm<Eigen::Infinity>());
      if (st_max < 1e-8) break;
      double sc = std::min(1.0, 4.0 / st_max);
      double slope = rhs_s.dot(step_s) + rhs_2.dot(step_2);
      for (int bt = 0; bt < 8; ++bt) {
        Eigen::VectorXd vs_try = v_s_cur + sc * step_s;
        Eigen::VectorXd v2_try = v_2_cur + sc * step_2;
        double phi_try = eval_data_ll_only(y, beta_eval, vs_try, v2_try, chance) -
          0.5 * (vs_try.squaredNorm() + v2_try.squaredNorm());
        if (phi_try >= phi_cur + 1e-4 * sc * std::max(slope, 0.0)) {
          v_s_cur = vs_try;
          v_2_cur = v2_try;
          break;
        }
        sc *= 0.5;
        if (bt == 7) {
          v_s_cur += sc * step_s;
          v_2_cur += sc * step_2;
        }
      }
    }

    eval_trials(y, beta_eval, v_s_cur, v_2_cur, chance, true);
    assemble_and_factor(beta_eval, v_s_cur, v_2_cur, P_beta, rhs_s, rhs_2, rhs_b, llt_22, true);

    L22_inv.setIdentity();
    L_22_saved.triangularView<Eigen::Lower>().solveInPlace(L22_inv);
    W_mat.noalias() = L22_inv * Bt_Ass_inv;

    for (int n = 0; n < N; ++n) {
      int s = subj_idx(n), i = item_idx(n);
      Eigen::VectorXd zs_t = ZsL.row(n).transpose();
      Eigen::VectorXd zi_t = ZiL.row(n).transpose();
      double h_s = zs_t.dot(A_ss_inv[s] * zs_t);
      q_tmp.noalias() = L22_inv.middleCols(i * q_i, q_i) * zi_t -
        W_mat.middleCols(s * q_s, q_s) * zs_t;
      if (has_testlet) {
        int t = testlet_idx(n);
        q_tmp.noalias() += L22_inv.middleCols(dim_i + t * q_t, q_t) * ZtL.row(n).transpose();
      }
      c_lap(n) = -0.5 * (h_s + q_tmp.squaredNorm()) * wp_vec(n);
    }

    cs_vec.setZero();
    c2_vec.setZero();
    for (int n = 0; n < N; ++n) {
      double cn = c_lap(n);
      cs_vec.segment(subj_idx(n) * q_s, q_s).noalias() += cn * ZsL.row(n).transpose();
      c2_vec.segment(item_idx(n) * q_i, q_i).noalias() += cn * ZiL.row(n).transpose();
      if (has_testlet) {
        c2_vec.segment(dim_i + testlet_idx(n) * q_t, q_t).noalias() += cn * ZtL.row(n).transpose();
      }
    }

    Eigen::VectorXd Ass_inv_cs(dim_s), Ass_inv_rs(dim_s);
    for (int s = 0; s < n_subj; ++s) {
      Ass_inv_cs.segment(s * q_s, q_s).noalias() = A_ss_inv[s] * cs_vec.segment(s * q_s, q_s);
      Ass_inv_rs.segment(s * q_s, q_s).noalias() = A_ss_inv[s] * rhs_s.segment(s * q_s, q_s);
    }
    Eigen::VectorXd S22_inv_c2 = llt_22.solve(c2_vec - Bt * Ass_inv_cs);
    Eigen::VectorXd S22_inv_r2 = llt_22.solve(rhs_2 - Bt * Ass_inv_rs);

    return (rhs_b + X_rm.transpose() * c_lap) -
      C_bs * (Ass_inv_rs + Ass_inv_cs) - S_2b_saved.transpose() * (S22_inv_r2 + S22_inv_c2);
  }

  FitResult fit(
      const Eigen::VectorXi& y,
      double chance,
      int corr_mode,
      bool regularized,
      double lkj_eta,
      double prior_beta_sd,
      double prior_tau_sd,
      int max_em_iter,
      int max_pirls_iter,
      double tol,
      const Eigen::VectorXd* init_beta,
      const Eigen::MatrixXd* init_Sigma_s,
      const Eigen::MatrixXd* init_Sigma_i,
      const Eigen::MatrixXd* init_Sigma_t,
      bool fix_cov = false,
      bool compute_joint_vcov = true,
      bool return_ranef = false
  ) {
    FitResult res;

    Eigen::VectorXd init_s = Eigen::VectorXd::Constant(q_s, 0.25);
    init_s(0) = 0.65;
    Eigen::MatrixXd Sigma_s = (init_s.array().square()).matrix().asDiagonal();
    if (init_Sigma_s != nullptr && init_Sigma_s->rows() == q_s && init_Sigma_s->cols() == q_s && init_Sigma_s->allFinite()) {
      Sigma_s = fix_cov ? (*init_Sigma_s) : (0.92 * (*init_Sigma_s) + 0.08 * Sigma_s);
    }

    Eigen::VectorXd init_i = Eigen::VectorXd::Constant(q_i, 0.25);
    init_i(0) = 0.55;
    Eigen::MatrixXd Sigma_i = (init_i.array().square()).matrix().asDiagonal();
    if (init_Sigma_i != nullptr && init_Sigma_i->rows() == q_i && init_Sigma_i->cols() == q_i && init_Sigma_i->allFinite()) {
      Sigma_i = fix_cov ? (*init_Sigma_i) : (0.92 * (*init_Sigma_i) + 0.08 * Sigma_i);
    }

    Eigen::MatrixXd Sigma_t;
    if (has_testlet) {
      Eigen::VectorXd init_t = Eigen::VectorXd::Constant(q_t, 0.20);
      init_t(0) = 0.40;
      Sigma_t = (init_t.array().square()).matrix().asDiagonal();
      if (init_Sigma_t != nullptr && init_Sigma_t->rows() == q_t && init_Sigma_t->cols() == q_t && init_Sigma_t->allFinite()) {
        Sigma_t = fix_cov ? (*init_Sigma_t) : (0.92 * (*init_Sigma_t) + 0.08 * Sigma_t);
      }
    }

    Eigen::MatrixXd P_beta = Eigen::MatrixXd::Zero(P, P);
    if (regularized && prior_beta_sd > 0.0) {
      double prec = 1.0 / (prior_beta_sd * prior_beta_sd);
      P_beta = prec * Eigen::MatrixXd::Identity(P, P);
      P_beta(0, 0) = 1.0 / ((1.5 * prior_beta_sd) * (1.5 * prior_beta_sd));
    }

    double mean_y = std::min(std::max(y.cast<double>().mean(), chance + 0.05), 0.95);
    Eigen::VectorXd beta = Eigen::VectorXd::Zero(P);
    beta(0) = std::log((mean_y - chance) / (1.0 - mean_y));
    if (init_beta != nullptr && init_beta->size() == P && init_beta->allFinite()) {
      beta = *init_beta;
    }

    Eigen::VectorXd v_s = Eigen::VectorXd::Zero(dim_s);
    Eigen::VectorXd v_2 = Eigen::VectorXd::Zero(dim_2);

    double log_lik_laplace = NA_REAL;
    double final_inner_step = 1.0;
    double final_outer_diff = 1.0;
    bool inner_converged = false;
    bool outer_converged = false;
    int outer_used = 0;

    Eigen::VectorXd rel_sum_s = Eigen::VectorXd::Zero(q_s);
    Eigen::VectorXd rel_sum_i = Eigen::VectorXd::Zero(q_i);
    Eigen::VectorXd rel_sum_t = Eigen::VectorXd::Zero(q_t);

    auto run_em_map = [&](
        const Eigen::MatrixXd& Sig_s_in,
        const Eigen::MatrixXd& Sig_i_in,
        const Eigen::MatrixXd& Sig_t_in,
        Eigen::MatrixXd& Sig_s_out,
        Eigen::MatrixXd& Sig_i_out,
        Eigen::MatrixXd& Sig_t_out,
        bool update_sigma
    ) {
      Eigen::MatrixXd L_s = safe_chol(Sig_s_in);
      Eigen::MatrixXd L_i = safe_chol(Sig_i_in);
      Eigen::MatrixXd L_t = has_testlet ? safe_chol(Sig_t_in) : Eigen::MatrixXd();

      ZsL.noalias() = Zs_rm * L_s;
      ZiL.noalias() = Zi_rm * L_i;
      if (has_testlet) ZtL.noalias() = Zt_rm * L_t;

      inner_converged = false;
      Eigen::VectorXd rhs_s(dim_s), rhs_2(dim_2), rhs_b(P);
      Eigen::LLT<Eigen::MatrixXd> llt_22;

      for (int it = 0; it < max_pirls_iter; ++it) {
        double ll_data = eval_trials(y, beta, v_s, v_2, chance, false);
        double logdet_H = assemble_and_factor(beta, v_s, v_2, P_beta, rhs_s, rhs_2, rhs_b, llt_22, false);
        double phi_cur = ll_data - 0.5 * (v_s.squaredNorm() + v_2.squaredNorm() + beta.dot(P_beta * beta));

        Eigen::VectorXd Ass_inv_rs(dim_s);
        for (int s = 0; s < n_subj; ++s) {
          Ass_inv_rs.segment(s * q_s, q_s).noalias() = A_ss_inv[s] * rhs_s.segment(s * q_s, q_s);
        }
        Eigen::VectorXd r_2_tilde = rhs_2 - Bt * Ass_inv_rs;
        Eigen::VectorXd S22_inv_r2 = llt_22.solve(r_2_tilde);
        Eigen::MatrixXd S22_inv_S2b = llt_22.solve(S_2b_saved);
        Eigen::VectorXd r_b_tilde = rhs_b - C_bs * Ass_inv_rs - S_2b_saved.transpose() * S22_inv_r2;

        Eigen::LDLT<Eigen::MatrixXd> ldlt_b(H_bb_schur);
        Eigen::VectorXd d_beta = ldlt_b.solve(r_b_tilde);
        Eigen::VectorXd d_v2 = S22_inv_r2 - S22_inv_S2b * d_beta;
        Eigen::VectorXd d_vs = Ass_inv_rs - Bt_Ass_inv.transpose() * d_v2 - C_bs_Ass_inv.transpose() * d_beta;

        double max_step = std::max({
          d_beta.lpNorm<Eigen::Infinity>(),
          d_v2.lpNorm<Eigen::Infinity>(),
          d_vs.lpNorm<Eigen::Infinity>()
        });
        final_inner_step = max_step;

        if (max_step < 1.5e-5 && (it >= 1 || outer_used >= 3)) {
          log_lik_laplace = ll_data - 0.5 * (v_s.squaredNorm() + v_2.squaredNorm() + logdet_H);
          inner_converged = true;
          break;
        }

        double step_scale = (max_step > 4.0) ? (4.0 / max_step) : 1.0;
        double slope = rhs_b.dot(d_beta) + rhs_s.dot(d_vs) + rhs_2.dot(d_v2);
        double accepted_ll = ll_data;

        // Armijo backtracking line search on penalized conditional log-posterior
        for (int bt = 0; bt < 10; ++bt) {
          Eigen::VectorXd b_try = beta + step_scale * d_beta;
          Eigen::VectorXd v2_try = v_2 + step_scale * d_v2;
          Eigen::VectorXd vs_try = v_s + step_scale * d_vs;
          double ll_try = eval_data_ll_only(y, b_try, vs_try, v2_try, chance);
          double phi_try = ll_try - 0.5 * (vs_try.squaredNorm() + v2_try.squaredNorm() + b_try.dot(P_beta * b_try));
          if (phi_try >= phi_cur + 1e-4 * step_scale * std::max(slope, 0.0) || (step_scale * max_step) < 1e-7) {
            beta = b_try;
            v_2 = v2_try;
            v_s = vs_try;
            accepted_ll = ll_try;
            break;
          }
          step_scale *= 0.5;
          if (bt == 9) {
            beta = b_try;
            v_2 = v2_try;
            v_s = vs_try;
            accepted_ll = ll_try;
          }
        }

        log_lik_laplace = accepted_ll - 0.5 * (v_s.squaredNorm() + v_2.squaredNorm() + logdet_H);

        if ((step_scale * max_step) < 1.5e-5 && max_step < 5.0e-4 && (it >= 1 || outer_used >= 3)) {
          inner_converged = true;
          break;
        }
      }

      if (chance > 0.0 && inner_converged) {
        eval_trials(y, beta, v_s, v_2, chance, true);
        assemble_and_factor(beta, v_s, v_2, P_beta, rhs_s, rhs_2, rhs_b, llt_22, true);
      }

      L22_inv.setIdentity();
      L_22_saved.triangularView<Eigen::Lower>().solveInPlace(L22_inv);
      W_mat.noalias() = L22_inv * Bt_Ass_inv;

      rel_sum_s.setZero();
      Eigen::MatrixXd sum_uu_s = Eigen::MatrixXd::Zero(q_s, q_s);
      for (int s = 0; s < n_subj; ++s) {
        const int off_s = s * q_s;
        Eigen::VectorXd us = L_s * v_s.segment(off_s, q_s);
        const auto Ws = W_mat.middleCols(off_s, q_s);
        Eigen::MatrixXd V_vs = A_ss_inv[s] + Ws.transpose() * Ws;
        for (int j = 0; j < q_s; ++j) {
          double rel = std::min(std::max(1.0 - V_vs(j, j), 0.0), 1.0);
          rel_sum_s(j) += 2.0 * rel * rel;
        }
        sum_uu_s.noalias() += us * us.transpose() + L_s * V_vs * L_s.transpose();
      }

      rel_sum_i.setZero();
      Eigen::MatrixXd sum_uu_i = Eigen::MatrixXd::Zero(q_i, q_i);
      for (int i = 0; i < n_item; ++i) {
        const int off_i = i * q_i;
        Eigen::VectorXd ui = L_i * v_2.segment(off_i, q_i);
        const auto Li_cols = L22_inv.middleCols(off_i, q_i);
        Eigen::MatrixXd V_vi = Li_cols.transpose() * Li_cols;
        for (int j = 0; j < q_i; ++j) {
          double rel = std::min(std::max(1.0 - V_vi(j, j), 0.0), 1.0);
          rel_sum_i(j) += 2.0 * rel * rel;
        }
        sum_uu_i.noalias() += ui * ui.transpose() + L_i * V_vi * L_i.transpose();
      }

      rel_sum_t.setZero();
      Eigen::MatrixXd sum_uu_t;
      if (has_testlet) {
        sum_uu_t = Eigen::MatrixXd::Zero(q_t, q_t);
        for (int t = 0; t < n_testlet; ++t) {
          const int off_t = dim_i + t * q_t;
          Eigen::VectorXd ut = L_t * v_2.segment(off_t, q_t);
          const auto Lt_cols = L22_inv.middleCols(off_t, q_t);
          Eigen::MatrixXd V_vt = Lt_cols.transpose() * Lt_cols;
          for (int j = 0; j < q_t; ++j) {
            double rel = std::min(std::max(1.0 - V_vt(j, j), 0.0), 1.0);
            rel_sum_t(j) += 2.0 * rel * rel;
          }
          sum_uu_t.noalias() += ut * ut.transpose() + L_t * V_vt * L_t.transpose();
        }
      }

      if (!update_sigma) {
        Sig_s_out = Sig_s_in;
        Sig_i_out = Sig_i_in;
        if (has_testlet) Sig_t_out = Sig_t_in;
        return;
      }

      Sig_s_out = update_cov_matrix(sum_uu_s, Sig_s_in, n_subj, corr_mode, regularized, lkj_eta, prior_tau_sd);
      Sig_i_out = update_cov_matrix(sum_uu_i, Sig_i_in, n_item, corr_mode, regularized, lkj_eta, prior_tau_sd);
      if (has_testlet) {
        Sig_t_out = update_cov_matrix(sum_uu_t, Sig_t_in, n_testlet, corr_mode, regularized, lkj_eta, prior_tau_sd);
      }

      Eigen::MatrixXd map_s = safe_chol(Sig_s_out).triangularView<Eigen::Lower>().solve(L_s);
      Eigen::MatrixXd map_i = safe_chol(Sig_i_out).triangularView<Eigen::Lower>().solve(L_i);
      for (int s = 0; s < n_subj; ++s) v_s.segment(s * q_s, q_s) = map_s * v_s.segment(s * q_s, q_s);
      for (int i = 0; i < n_item; ++i) v_2.segment(i * q_i, q_i) = map_i * v_2.segment(i * q_i, q_i);
      if (has_testlet) {
        Eigen::MatrixXd map_t = safe_chol(Sig_t_out).triangularView<Eigen::Lower>().solve(L_t);
        for (int t = 0; t < n_testlet; ++t) {
          v_2.segment(dim_i + t * q_t, q_t) = map_t * v_2.segment(dim_i + t * q_t, q_t);
        }
      }
    };

    if (fix_cov) {
      Eigen::MatrixXd S_s1, S_i1, S_t1;
      outer_used = 3;
      run_em_map(Sigma_s, Sigma_i, Sigma_t, S_s1, S_i1, S_t1, false);
      outer_used = 1;
      final_outer_diff = 0.0;
      outer_converged = inner_converged;
    } else {
      while (outer_used < max_em_iter) {
        if ((outer_used & 7) == 0) Rcpp::checkUserInterrupt();
        Eigen::VectorXd beta_prev = beta;
        Eigen::MatrixXd S_s1, S_i1, S_t1;
        outer_used++;
        run_em_map(Sigma_s, Sigma_i, Sigma_t, S_s1, S_i1, S_t1, true);

        double d1 = std::max(
          (S_s1 - Sigma_s).lpNorm<Eigen::Infinity>(),
          (S_i1 - Sigma_i).lpNorm<Eigen::Infinity>()
        );
        if (has_testlet) d1 = std::max(d1, (S_t1 - Sigma_t).lpNorm<Eigen::Infinity>());
        final_outer_diff = std::max(d1, (beta - beta_prev).lpNorm<Eigen::Infinity>());

        if (final_outer_diff < tol && outer_used >= 3 && inner_converged) {
          Sigma_s = S_s1; Sigma_i = S_i1; if (has_testlet) Sigma_t = S_t1;
          outer_converged = true; break;
        }
        if (outer_used >= max_em_iter) {
          Sigma_s = S_s1; Sigma_i = S_i1; if (has_testlet) Sigma_t = S_t1; break;
        }

        beta_prev = beta;
        Eigen::MatrixXd S_s2, S_i2, S_t2;
        outer_used++;
        run_em_map(S_s1, S_i1, S_t1, S_s2, S_i2, S_t2, true);

        double d2 = std::max(
          (S_s2 - S_s1).lpNorm<Eigen::Infinity>(),
          (S_i2 - S_i1).lpNorm<Eigen::Infinity>()
        );
        if (has_testlet) d2 = std::max(d2, (S_t2 - S_t1).lpNorm<Eigen::Infinity>());
        final_outer_diff = std::max(d2, (beta - beta_prev).lpNorm<Eigen::Infinity>());

        if (final_outer_diff < tol && outer_used >= 3 && inner_converged) {
          Sigma_s = S_s2; Sigma_i = S_i2; if (has_testlet) Sigma_t = S_t2;
          outer_converged = true; break;
        }
        if (outer_used >= max_em_iter) {
          Sigma_s = S_s2; Sigma_i = S_i2; if (has_testlet) Sigma_t = S_t2; break;
        }

        Eigen::MatrixXd L_s0 = safe_chol(Sigma_s), L_s1 = safe_chol(S_s1), L_s2 = safe_chol(S_s2);
        Eigen::MatrixXd L_i0 = safe_chol(Sigma_i), L_i1 = safe_chol(S_i1), L_i2 = safe_chol(S_i2);
        double r_norm2 = (L_s1 - L_s0).squaredNorm() + (L_i1 - L_i0).squaredNorm();
        double v_norm2 = ((L_s2 - L_s1) - (L_s1 - L_s0)).squaredNorm() +
          ((L_i2 - L_i1) - (L_i1 - L_i0)).squaredNorm();
        if (has_testlet) {
          Eigen::MatrixXd L_t0 = safe_chol(Sigma_t), L_t1 = safe_chol(S_t1), L_t2 = safe_chol(S_t2);
          r_norm2 += (L_t1 - L_t0).squaredNorm();
          v_norm2 += ((L_t2 - L_t1) - (L_t1 - L_t0)).squaredNorm();
        }

        if (v_norm2 > 1e-16 && r_norm2 > 1e-16 && outer_used < max_em_iter) {
          double alpha_sq = std::max(std::min(-std::sqrt(r_norm2 / v_norm2), -1.0), -20.0);
          if (alpha_sq < -1.05) {
            Eigen::VectorXd v_s_bak = v_s, v_2_bak = v_2, beta_bak = beta;
            Eigen::MatrixXd Sp_s = squarem_extrapolate_cov(Sigma_s, S_s1, S_s2, alpha_sq, corr_mode);
            Eigen::MatrixXd Sp_i = squarem_extrapolate_cov(Sigma_i, S_i1, S_i2, alpha_sq, corr_mode);
            Eigen::MatrixXd Sp_t = has_testlet ?
              squarem_extrapolate_cov(Sigma_t, S_t1, S_t2, alpha_sq, corr_mode) : Eigen::MatrixXd();

            Eigen::MatrixXd map_s = safe_chol(Sp_s).triangularView<Eigen::Lower>().solve(L_s2);
            Eigen::MatrixXd map_i = safe_chol(Sp_i).triangularView<Eigen::Lower>().solve(L_i2);
            for (int s = 0; s < n_subj; ++s) v_s.segment(s * q_s, q_s) = map_s * v_s.segment(s * q_s, q_s);
            for (int i = 0; i < n_item; ++i) v_2.segment(i * q_i, q_i) = map_i * v_2.segment(i * q_i, q_i);
            if (has_testlet) {
              Eigen::MatrixXd L_t2 = safe_chol(S_t2);
              Eigen::MatrixXd map_t = safe_chol(Sp_t).triangularView<Eigen::Lower>().solve(L_t2);
              for (int t = 0; t < n_testlet; ++t) {
                v_2.segment(dim_i + t * q_t, q_t) = map_t * v_2.segment(dim_i + t * q_t, q_t);
              }
            }

            Eigen::MatrixXd S_s3, S_i3, S_t3;
            outer_used++;
            run_em_map(Sp_s, Sp_i, Sp_t, S_s3, S_i3, S_t3, true);

            double d3 = std::max(
              (S_s3 - Sp_s).lpNorm<Eigen::Infinity>(),
              (S_i3 - Sp_i).lpNorm<Eigen::Infinity>()
            );
            if (has_testlet) d3 = std::max(d3, (S_t3 - Sp_t).lpNorm<Eigen::Infinity>());

            if (d3 <= 4.0 * std::max(d1, d2) || d3 < 5e-3) {
              final_outer_diff = d3;
              Sigma_s = S_s3; Sigma_i = S_i3; if (has_testlet) Sigma_t = S_t3;
              if (d3 < 0.5 * tol && inner_converged) {
                outer_converged = true; break;
              }
              continue;
            } else {
              v_s = v_s_bak; v_2 = v_2_bak; beta = beta_bak;
            }
          }
        }

        Sigma_s = S_s2; Sigma_i = S_i2; if (has_testlet) Sigma_t = S_t2;
      }
    }

    Eigen::MatrixXd L_s_fin = safe_chol(Sigma_s);
    Eigen::MatrixXd L_i_fin = safe_chol(Sigma_i);
    Eigen::MatrixXd L_t_fin = has_testlet ? safe_chol(Sigma_t) : Eigen::MatrixXd();

    // Post-convergence Observed-Hessian Laplace Polish (up to 3 iterations)
    if (inner_converged && outer_converged) {
      ZsL.noalias() = Zs_rm * L_s_fin;
      ZiL.noalias() = Zi_rm * L_i_fin;
      if (has_testlet) ZtL.noalias() = Zt_rm * L_t_fin;

      for (int polish = 0; polish < 3; ++polish) {
        eval_trials(y, beta, v_s, v_2, chance, true);
        Eigen::VectorXd rhs_s(dim_s), rhs_2(dim_2), rhs_b(P);
        Eigen::LLT<Eigen::MatrixXd> llt_22;
        assemble_and_factor(beta, v_s, v_2, P_beta, rhs_s, rhs_2, rhs_b, llt_22, true);

        L22_inv.setIdentity();
        L_22_saved.triangularView<Eigen::Lower>().solveInPlace(L22_inv);
        W_mat.noalias() = L22_inv * Bt_Ass_inv;

        for (int n = 0; n < N; ++n) {
          int s = subj_idx(n), i = item_idx(n);
          Eigen::VectorXd zs_t = ZsL.row(n).transpose();
          Eigen::VectorXd zi_t = ZiL.row(n).transpose();
          double h_s = zs_t.dot(A_ss_inv[s] * zs_t);
          q_tmp.noalias() = L22_inv.middleCols(i * q_i, q_i) * zi_t -
            W_mat.middleCols(s * q_s, q_s) * zs_t;
          if (has_testlet) {
            int t = testlet_idx(n);
            q_tmp.noalias() += L22_inv.middleCols(dim_i + t * q_t, q_t) * ZtL.row(n).transpose();
          }
          c_lap(n) = -0.5 * (h_s + q_tmp.squaredNorm()) * wp_vec(n);
        }

        cs_vec.setZero();
        c2_vec.setZero();
        for (int n = 0; n < N; ++n) {
          double cn = c_lap(n);
          cs_vec.segment(subj_idx(n) * q_s, q_s).noalias() += cn * ZsL.row(n).transpose();
          c2_vec.segment(item_idx(n) * q_i, q_i).noalias() += cn * ZiL.row(n).transpose();
          if (has_testlet) {
            c2_vec.segment(dim_i + testlet_idx(n) * q_t, q_t).noalias() += cn * ZtL.row(n).transpose();
          }
        }

        Eigen::VectorXd Ass_inv_cs(dim_s), Ass_inv_rs(dim_s);
        for (int s = 0; s < n_subj; ++s) {
          Ass_inv_cs.segment(s * q_s, q_s).noalias() = A_ss_inv[s] * cs_vec.segment(s * q_s, q_s);
          Ass_inv_rs.segment(s * q_s, q_s).noalias() = A_ss_inv[s] * rhs_s.segment(s * q_s, q_s);
        }
        Eigen::VectorXd S22_inv_c2 = llt_22.solve(c2_vec - Bt * Ass_inv_cs);
        Eigen::VectorXd S22_inv_r2 = llt_22.solve(rhs_2 - Bt * Ass_inv_rs);
        Eigen::MatrixXd S22_inv_S2b = llt_22.solve(S_2b_saved);

        Eigen::VectorXd grad_lap_b = (rhs_b + X_rm.transpose() * c_lap) -
          C_bs * (Ass_inv_rs + Ass_inv_cs) - S_2b_saved.transpose() * (S22_inv_r2 + S22_inv_c2);

        Eigen::VectorXd d_beta_lap = H_bb_schur.ldlt().solve(grad_lap_b);
        double lap_norm = d_beta_lap.lpNorm<Eigen::Infinity>();
        if (lap_norm < 1.5) {
          Eigen::VectorXd d_v2 = S22_inv_r2 - S22_inv_S2b * d_beta_lap;
          Eigen::VectorXd d_vs = Ass_inv_rs - Bt_Ass_inv.transpose() * d_v2 - C_bs_Ass_inv.transpose() * d_beta_lap;
          beta += d_beta_lap;
          v_2 += d_v2;
          v_s += d_vs;
        }
        if (lap_norm < 1e-5) break;
      }

      // Evaluate exact conditional mode, determinant, and Laplace log-likelihood
      // at the final corrected beta and covariance matrices (with Armijo safeguard).
      Eigen::VectorXd rhs_s(dim_s), rhs_2(dim_2), rhs_b(P);
      Eigen::LLT<Eigen::MatrixXd> llt_22;
      inner_converged = false;
      for (int it = 0; it < max_pirls_iter; ++it) {
        double ll_cur = eval_trials(y, beta, v_s, v_2, chance, false);
        double phi_cur = ll_cur - 0.5 * (v_s.squaredNorm() + v_2.squaredNorm());
        assemble_and_factor(beta, v_s, v_2, P_beta, rhs_s, rhs_2, rhs_b, llt_22, true);
        Eigen::VectorXd a(dim_s);
        for (int s = 0; s < n_subj; ++s) {
          a.segment(s * q_s, q_s).noalias() = A_ss_inv[s] * rhs_s.segment(s * q_s, q_s);
        }
        Eigen::VectorXd step_2 = llt_22.solve(rhs_2 - Bt * a);
        Eigen::VectorXd step_s = a - Bt_Ass_inv.transpose() * step_2;
        final_inner_step = std::max(step_s.lpNorm<Eigen::Infinity>(), step_2.lpNorm<Eigen::Infinity>());
        if (final_inner_step < 1.5e-5) {
          const double ll_obs = eval_trials(y, beta, v_s, v_2, chance, true);
          const double logdet_obs = assemble_and_factor(beta, v_s, v_2, P_beta, rhs_s, rhs_2, rhs_b, llt_22, true);
          log_lik_laplace = ll_obs - 0.5 * (v_s.squaredNorm() + v_2.squaredNorm() + logdet_obs);
          inner_converged = true;
          break;
        }
        double sc = std::min(1.0, 4.0 / final_inner_step);
        double slope = rhs_s.dot(step_s) + rhs_2.dot(step_2);
        for (int bt = 0; bt < 8; ++bt) {
          Eigen::VectorXd vs_try = v_s + sc * step_s;
          Eigen::VectorXd v2_try = v_2 + sc * step_2;
          double phi_try = eval_data_ll_only(y, beta, vs_try, v2_try, chance) -
            0.5 * (vs_try.squaredNorm() + v2_try.squaredNorm());
          if (phi_try >= phi_cur + 1e-4 * sc * std::max(slope, 0.0) || (sc * final_inner_step) < 1e-7) {
            v_s = vs_try;
            v_2 = v2_try;
            break;
          }
          sc *= 0.5;
          if (bt == 7) {
            v_s += sc * step_s;
            v_2 += sc * step_2;
          }
        }
      }
    }

    Eigen::MatrixXd H_bb_saved = H_bb_schur;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_cond(0.5 * (H_bb_saved + H_bb_saved.transpose()));
    bool pos_def_b = (es_cond.eigenvalues().minCoeff() > 1e-8) && std::isfinite(log_lik_laplace);
    Eigen::VectorXd ev_cond = es_cond.eigenvalues().cwiseMax(1e-8);
    Eigen::MatrixXd vcov_cond = es_cond.eigenvectors() * ev_cond.cwiseInverse().asDiagonal() *
      es_cond.eigenvectors().transpose();

    // Asymptotic variances of log(tau) from per-cluster reliability & prior curvature
    double s0_sq = std::max(prior_tau_sd * prior_tau_sd, 0.25);
    auto make_var_log_tau = [&](const Eigen::MatrixXd& Sig, const Eigen::VectorXd& rel_sum) {
      int q = Sig.rows();
      Eigen::VectorXd out(q);
      for (int j = 0; j < q; ++j) {
        double tau2 = std::max(Sig(j, j), 1e-10);
        double info_prior = 0.0;
        if (regularized) {
          // Half-Student-t(df=3, scale=s0) prior curvature on log(tau) + boundary stabilization
          double r = tau2 / (3.0 * s0_sq);
          info_prior = (8.0 * r) / ((1.0 + r) * (1.0 + r)) + 2.0 * (tau2 / s0_sq);
        }
        double info_tot = std::max(rel_sum(j) + info_prior, 0.75);
        out(j) = 1.0 / info_tot;
      }
      return out;
    };

    res.var_log_tau_s = make_var_log_tau(Sigma_s, rel_sum_s);
    res.var_log_tau_i = make_var_log_tau(Sigma_i, rel_sum_i);
    res.var_log_tau_t = has_testlet ? make_var_log_tau(Sigma_t, rel_sum_t) : Eigen::VectorXd(0);

    int n_theta = q_s + q_i + q_t;
    res.J_beta_theta = Eigen::MatrixXd::Zero(P, n_theta);
    Eigen::MatrixXd vcov_beta = vcov_cond;

    // Exact Tierney-Kadane Laplace Hessian wrt beta + profile covariance sensitivity J_beta_theta
    if (compute_joint_vcov && inner_converged && outer_converged && pos_def_b) {
      const Eigen::VectorXd v_s_opt = v_s;
      const Eigen::VectorXd v_2_opt = v_2;
      const double eps_b = 5e-4;

      const Eigen::VectorXd g_base = compute_laplace_grad_beta(
        y, beta, L_s_fin, L_i_fin, L_t_fin, P_beta, chance, v_s_opt, v_2_opt
      );

      Eigen::MatrixXd H_lap_bb(P, P);
      for (int p = 0; p < P; ++p) {
        Eigen::VectorXd b_hi = beta;
        b_hi(p) += eps_b;
        Eigen::VectorXd g_hi = compute_laplace_grad_beta(
          y, b_hi, L_s_fin, L_i_fin, L_t_fin, P_beta, chance, v_s_opt, v_2_opt
        );
        H_lap_bb.col(p) = -(g_hi - g_base) / eps_b;
      }
      H_lap_bb = 0.5 * (H_lap_bb + H_lap_bb.transpose());

      // Safeguard Laplace Hessian eigenvalues against pathological finite-sample flatness
      double min_floor = std::max(0.25 * ev_cond.minCoeff(), 1e-6);
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_lap(H_lap_bb);
      Eigen::VectorXd ev_lap = es_lap.eigenvalues().cwiseMax(min_floor);
      Eigen::MatrixXd vcov_lap_fixed = es_lap.eigenvectors() * ev_lap.cwiseInverse().asDiagonal() *
        es_lap.eigenvectors().transpose();

      if (fix_cov) {
        vcov_beta = vcov_lap_fixed;
      } else {
        const double eps_t = 1e-3;
        const double mult_hi = std::exp(eps_t);
        Eigen::VectorXd var_theta_all(n_theta);
        int col_k = 0;

        for (int j = 0; j < q_s; ++j, ++col_k) {
          var_theta_all(col_k) = res.var_log_tau_s(j);
          Eigen::MatrixXd Ls_hi = safe_chol(scale_cov_single_sd(Sigma_s, j, mult_hi));
          Eigen::VectorXd g_hi = compute_laplace_grad_beta(
            y, beta, Ls_hi, L_i_fin, L_t_fin, P_beta, chance, v_s_opt, v_2_opt
          );
          res.J_beta_theta.col(col_k) = vcov_lap_fixed * ((g_hi - g_base) / eps_t);
        }

        for (int j = 0; j < q_i; ++j, ++col_k) {
          var_theta_all(col_k) = res.var_log_tau_i(j);
          Eigen::MatrixXd Li_hi = safe_chol(scale_cov_single_sd(Sigma_i, j, mult_hi));
          Eigen::VectorXd g_hi = compute_laplace_grad_beta(
            y, beta, L_s_fin, Li_hi, L_t_fin, P_beta, chance, v_s_opt, v_2_opt
          );
          res.J_beta_theta.col(col_k) = vcov_lap_fixed * ((g_hi - g_base) / eps_t);
        }

        for (int j = 0; j < q_t; ++j, ++col_k) {
          var_theta_all(col_k) = res.var_log_tau_t(j);
          Eigen::MatrixXd Lt_hi = safe_chol(scale_cov_single_sd(Sigma_t, j, mult_hi));
          Eigen::VectorXd g_hi = compute_laplace_grad_beta(
            y, beta, L_s_fin, L_i_fin, Lt_hi, P_beta, chance, v_s_opt, v_2_opt
          );
          res.J_beta_theta.col(col_k) = vcov_lap_fixed * ((g_hi - g_base) / eps_t);
        }

        vcov_beta = vcov_lap_fixed +
          res.J_beta_theta * var_theta_all.asDiagonal() * res.J_beta_theta.transpose();
      }
      v_s = v_s_opt;
      v_2 = v_2_opt;
    }

    Eigen::VectorXd se_beta = vcov_beta.diagonal().cwiseMax(1e-12).cwiseSqrt();

    bool singular = (!regularized) && (!fix_cov) && (
      is_cov_singular(Sigma_s, corr_mode) ||
      is_cov_singular(Sigma_i, corr_mode) ||
      (has_testlet && is_cov_singular(Sigma_t, corr_mode))
    );

    if (!pos_def_b || !beta.allFinite() || !se_beta.allFinite()) {
      res.status = "non_pd_hessian";
      res.strictly_converged = false;
    } else if (!inner_converged) {
      res.status = "inner_not_converged";
      res.strictly_converged = false;
    } else if (!outer_converged) {
      res.status = "max_iter_finite";
      res.strictly_converged = false;
    } else {
      res.status = "converged";
      res.strictly_converged = true;
    }

    if (return_ranef) {
      res.u_subj.resize(n_subj, q_s);
      for (int s = 0; s < n_subj; ++s) {
        res.u_subj.row(s) = (L_s_fin * v_s.segment(s * q_s, q_s)).transpose();
      }
      res.u_item.resize(n_item, q_i);
      for (int i = 0; i < n_item; ++i) {
        res.u_item.row(i) = (L_i_fin * v_2.segment(i * q_i, q_i)).transpose();
      }
      if (has_testlet) {
        res.u_testlet.resize(n_testlet, q_t);
        for (int t = 0; t < n_testlet; ++t) {
          res.u_testlet.row(t) = (L_t_fin * v_2.segment(dim_i + t * q_t, q_t)).transpose();
        }
      }
    }

    res.beta = beta;
    res.se_beta = se_beta;
    res.vcov_beta = vcov_beta;
    res.vcov_cond = vcov_cond;
    res.Sigma_s = Sigma_s;
    res.Sigma_i = Sigma_i;
    res.Sigma_t = Sigma_t;
    res.log_lik_laplace = log_lik_laplace;
    res.final_outer_diff = final_outer_diff;
    res.final_inner_step = final_inner_step;
    res.outer_used = outer_used;
    res.singular = singular;
    return res;
  }
};

double eval_marginal_contrast_cpp(
    const Eigen::VectorXd& beta,
    const Eigen::MatrixXd& Sigma_s,
    const Eigen::MatrixXd& Sigma_i,
    const Eigen::MatrixXd& Sigma_t,
    bool has_testlet,
    const Eigen::MatrixXd& X_cells,
    const Eigen::MatrixXd& Zs_cells,
    const Eigen::MatrixXd& Zi_cells,
    const Eigen::MatrixXd& Zt_cells,
    const Eigen::VectorXd& weights_contrast,
    const Eigen::VectorXd& gh_nodes,
    const Eigen::VectorXd& gh_weights,
    double chance,
    bool marginal
) {
  int M = X_cells.rows();
  int K = gh_nodes.size();
  Eigen::VectorXd mu_c = X_cells * beta;

  if (!marginal) {
    double diff = 0.0;
    for (int m = 0; m < M; ++m) {
      double p_m = chance + (1.0 - chance) * inv_logit(mu_c(m));
      diff += weights_contrast(m) * p_m;
    }
    return diff;
  }

  double diff = 0.0;
  for (int m = 0; m < M; ++m) {
    double var_tot = Zs_cells.row(m).dot(Sigma_s * Zs_cells.row(m).transpose()) +
      Zi_cells.row(m).dot(Sigma_i * Zi_cells.row(m).transpose());
    if (has_testlet && Zt_cells.cols() > 0 && Sigma_t.rows() > 0) {
      var_tot += Zt_cells.row(m).dot(Sigma_t * Zt_cells.row(m).transpose());
    }
    double sd_tot = std::sqrt(std::max(var_tot, 0.0));
    double p_avg = 0.0;
    for (int k = 0; k < K; ++k) {
      p_avg += gh_weights(k) * inv_logit(mu_c(m) + sd_tot * gh_nodes(k));
    }
    diff += weights_contrast(m) * (chance + (1.0 - chance) * p_avg);
  }
  return diff;
}

// Match R's pivoted Cholesky lower factor `safe_chol_lower` for positive-definite matrices
Eigen::MatrixXd chol_lower_pd(const Eigen::MatrixXd& Sigma) {
  int p = Sigma.rows();
  if (p == 0) return Eigen::MatrixXd(0, 0);
  if (Sigma.cwiseAbs().maxCoeff() == 0.0) return Eigen::MatrixXd::Zero(p, p);
  Eigen::LLT<Eigen::MatrixXd> llt(Sigma);
  if (llt.info() == Eigen::Success) {
    return llt.matrixL();
  }
  return safe_chol(Sigma, 1e-12);
}

} // anonymous namespace

// Cache-Aligned Safeguarded SQUAREM + Observed-Hessian Laplace (Tierney-Kadane)
// Estimator for Crossed Binary & m-AFC GLMMs with Profile Covariance Propagation.
//
// [[Rcpp::export]]
Rcpp::List cpp_fit_crossed_glmm(
    const Eigen::VectorXi& y,
    const Eigen::MatrixXd& X,
    const Eigen::VectorXi& subj_idx,
    const Eigen::MatrixXd& Z_subj,
    int n_subj,
    const Eigen::VectorXi& item_idx,
    const Eigen::MatrixXd& Z_item,
    int n_item,
    const Eigen::VectorXi& testlet_idx,
    const Eigen::MatrixXd& Z_testlet,
    int n_testlet,
    double chance = 0.0,
    int corr_mode = 1,
    bool regularized = false,
    double lkj_eta = 2.0,
    double prior_beta_sd = 2.5,
    double prior_tau_sd = 1.5,
    int max_em_iter = 160,
    int max_pirls_iter = 22,
    double tol = 1.5e-4,
    Rcpp::Nullable<Eigen::VectorXd> init_beta = R_NilValue,
    Rcpp::Nullable<Eigen::MatrixXd> init_Sigma_subj = R_NilValue,
    Rcpp::Nullable<Eigen::MatrixXd> init_Sigma_item = R_NilValue,
    Rcpp::Nullable<Eigen::MatrixXd> init_Sigma_testlet = R_NilValue,
    bool fix_cov = false,
    bool compute_joint_vcov = true,
    bool return_ranef = false
) {
  GLMMWorkspace ws(X, subj_idx, Z_subj, n_subj, item_idx, Z_item, n_item, testlet_idx, Z_testlet, n_testlet);

  Eigen::VectorXd bw;
  Eigen::MatrixXd Ss_w, Si_w, St_w;
  const Eigen::VectorXd* p_bw = nullptr;
  const Eigen::MatrixXd* p_Ss = nullptr;
  const Eigen::MatrixXd* p_Si = nullptr;
  const Eigen::MatrixXd* p_St = nullptr;

  if (init_beta.isNotNull()) { bw = Rcpp::as<Eigen::VectorXd>(init_beta); p_bw = &bw; }
  if (init_Sigma_subj.isNotNull()) { Ss_w = Rcpp::as<Eigen::MatrixXd>(init_Sigma_subj); p_Ss = &Ss_w; }
  if (init_Sigma_item.isNotNull()) { Si_w = Rcpp::as<Eigen::MatrixXd>(init_Sigma_item); p_Si = &Si_w; }
  if (init_Sigma_testlet.isNotNull()) { St_w = Rcpp::as<Eigen::MatrixXd>(init_Sigma_testlet); p_St = &St_w; }

  FitResult res = ws.fit(
    y, chance, corr_mode, regularized, lkj_eta, prior_beta_sd, prior_tau_sd,
    max_em_iter, max_pirls_iter, tol, p_bw, p_Ss, p_Si, p_St,
    fix_cov, compute_joint_vcov, return_ranef
  );

  return Rcpp::List::create(
    Rcpp::Named("beta") = res.beta,
    Rcpp::Named("se") = res.se_beta,
    Rcpp::Named("vcov") = res.vcov_beta,
    Rcpp::Named("vcov_cond") = res.vcov_cond,
    Rcpp::Named("Sigma_subj") = res.Sigma_s,
    Rcpp::Named("Sigma_item") = res.Sigma_i,
    Rcpp::Named("Sigma_testlet") = ws.has_testlet ? Rcpp::wrap(res.Sigma_t) : R_NilValue,
    Rcpp::Named("var_log_tau_subj") = res.var_log_tau_s,
    Rcpp::Named("var_log_tau_item") = res.var_log_tau_i,
    Rcpp::Named("var_log_tau_testlet") = ws.has_testlet ? Rcpp::wrap(res.var_log_tau_t) : R_NilValue,
    Rcpp::Named("J_beta_theta") = res.J_beta_theta,
    Rcpp::Named("u_subj") = return_ranef ? Rcpp::wrap(res.u_subj) : R_NilValue,
    Rcpp::Named("u_item") = return_ranef ? Rcpp::wrap(res.u_item) : R_NilValue,
    Rcpp::Named("u_testlet") = (return_ranef && ws.has_testlet) ? Rcpp::wrap(res.u_testlet) : R_NilValue,
    Rcpp::Named("logLik") = res.log_lik_laplace,
    Rcpp::Named("singular") = res.singular,
    Rcpp::Named("converged") = res.strictly_converged,
    Rcpp::Named("status") = res.status,
    Rcpp::Named("iterations") = res.outer_used,
    Rcpp::Named("outer_diff") = res.final_outer_diff,
    Rcpp::Named("inner_step") = res.final_inner_step
  );
}

// Pre-Allocated Zero-Copy C++ Parametric Bootstrap Batch for Marginal Probability Contrasts
//
// [[Rcpp::export]]
Rcpp::NumericVector cpp_prob_bootstrap_batch(
    int n_boot,
    const Eigen::MatrixXd& X,
    const Eigen::VectorXi& subj_idx,
    const Eigen::MatrixXd& Z_subj,
    int n_subj,
    const Eigen::VectorXi& item_idx,
    const Eigen::MatrixXd& Z_item,
    int n_item,
    const Eigen::VectorXi& testlet_idx,
    const Eigen::MatrixXd& Z_testlet,
    int n_testlet,
    const Eigen::VectorXd& beta_hat,
    const Eigen::MatrixXd& Sigma_subj_hat,
    const Eigen::MatrixXd& Sigma_item_hat,
    Rcpp::Nullable<Eigen::MatrixXd> Sigma_testlet_hat,
    const Eigen::MatrixXd& X_cells,
    const Eigen::MatrixXd& Zs_cells,
    const Eigen::MatrixXd& Zi_cells,
    const Eigen::MatrixXd& Zt_cells,
    const Eigen::VectorXd& contrast_weights,
    const Eigen::VectorXd& gh_nodes,
    const Eigen::VectorXd& gh_weights,
    double chance = 0.0,
    int corr_mode = 1,
    bool regularized = false,
    double lkj_eta = 2.0,
    double prior_beta_sd = 2.5,
    double prior_tau_sd = 1.5,
    int max_em_iter = 160,
    int max_pirls_iter = 22,
    double tol = 1.5e-4,
    int eff_idx_1based = 1,
    bool marginal = true
) {
  Rcpp::RNGScope rng_scope;

  GLMMWorkspace ws(X, subj_idx, Z_subj, n_subj, item_idx, Z_item, n_item, testlet_idx, Z_testlet, n_testlet);
  const int N = X.rows();
  const int q_s = Z_subj.cols();
  const int q_i = Z_item.cols();
  const bool has_t = ws.has_testlet;
  const int q_t = ws.q_t;
  const int eff_idx0 = std::max(0, eff_idx_1based - 1);

  Eigen::MatrixXd St_hat = (has_t && Sigma_testlet_hat.isNotNull()) ?
    Rcpp::as<Eigen::MatrixXd>(Sigma_testlet_hat) : Eigen::MatrixXd(0, 0);

  // Precompute Z %*% chol_lower(Sigma_hat) matching R's `Z %*% safe_chol_lower(Sigma)`
  RowMatrixXd Gs = ws.Zs_rm * chol_lower_pd(Sigma_subj_hat);
  RowMatrixXd Gi = ws.Zi_rm * chol_lower_pd(Sigma_item_hat);
  RowMatrixXd Gt = has_t ? RowMatrixXd(ws.Zt_rm * chol_lower_pd(St_hat)) : RowMatrixXd(0, 0);
  Eigen::VectorXd fixed_eta = X * beta_hat;

  Eigen::MatrixXd u_s(n_subj, q_s);
  Eigen::MatrixXd u_i(n_item, q_i);
  Eigen::MatrixXd u_t(has_t ? n_testlet : 0, has_t ? q_t : 0);
  Eigen::VectorXi y_star(N);

  const Eigen::MatrixXd* p_St = has_t ? &St_hat : nullptr;
  Rcpp::NumericVector boot_est(n_boot, NA_REAL);

  for (int b = 0; b < n_boot; ++b) {
    if ((b & 3) == 0) Rcpp::checkUserInterrupt();
    // Draw standard normal cluster effects in exact column-major order matching R's matrix(rnorm(n*q), nrow=n)
    for (int c = 0; c < q_s; ++c) {
      for (int r = 0; r < n_subj; ++r) u_s(r, c) = R::rnorm(0.0, 1.0);
    }
    for (int c = 0; c < q_i; ++c) {
      for (int r = 0; r < n_item; ++r) u_i(r, c) = R::rnorm(0.0, 1.0);
    }
    if (has_t) {
      for (int c = 0; c < q_t; ++c) {
        for (int r = 0; r < n_testlet; ++r) u_t(r, c) = R::rnorm(0.0, 1.0);
      }
    }

    for (int n = 0; n < N; ++n) {
      double eta = fixed_eta(n) +
        Gs.row(n).dot(u_s.row(subj_idx(n))) +
        Gi.row(n).dot(u_i.row(item_idx(n)));
      if (has_t) {
        eta += Gt.row(n).dot(u_t.row(testlet_idx(n)));
      }
      double prob = chance + (1.0 - chance) * inv_logit(eta);
      prob = std::min(std::max(prob, 0.0), 1.0);
      y_star(n) = static_cast<int>(R::rbinom(1.0, prob));
    }

    try {
      FitResult f = ws.fit(
        y_star, chance, corr_mode, regularized, lkj_eta, prior_beta_sd, prior_tau_sd,
        max_em_iter, max_pirls_iter, tol, &beta_hat, &Sigma_subj_hat, &Sigma_item_hat, p_St,
        false, false, false
      );
      bool conv = f.strictly_converged && (f.status == "converged") &&
        std::isfinite(f.beta(eff_idx0)) && std::isfinite(f.se_beta(eff_idx0)) && (f.se_beta(eff_idx0) > 0.0);
      if (conv) {
        boot_est[b] = eval_marginal_contrast_cpp(
          f.beta, f.Sigma_s, f.Sigma_i, f.Sigma_t, has_t,
          X_cells, Zs_cells, Zi_cells, Zt_cells,
          contrast_weights, gh_nodes, gh_weights, chance, marginal
        );
      }
    } catch (const Rcpp::internal::InterruptedException&) {
      throw;
    } catch (...) {
      boot_est[b] = NA_REAL;
    }
  }

  return boot_est;
}
