#include <RcppEigen.h>
#include <Rmath.h>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>

// [[Rcpp::depends(RcppEigen)]]

namespace {

// The Laplace determinant is first-order sensitive to error in the conditional mode.
constexpr double conditional_mode_tol = 2e-7;
constexpr int max_cell_archetypes = 32;

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

inline double softplus(double x) {
  return (x > 0.0) ? (x + std::log1p(std::exp(-x))) : std::log1p(std::exp(x));
}

// Recursive GEMM-based lower-triangular inversion using a preallocated scratch buffer (zero heap allocs)
void invert_lower_rec_buf(
    const Eigen::Ref<const Eigen::MatrixXd>& L,
    Eigen::Ref<Eigen::MatrixXd> Linv,
    Eigen::MatrixXd& scratch
) {
  const int n = L.rows();
  if (n <= 32) {
    Linv.setIdentity();
    L.triangularView<Eigen::Lower>().solveInPlace(Linv);
    return;
  }
  const int n1 = (n / 2 + 3) & ~3;
  const int n2 = n - n1;
  invert_lower_rec_buf(L.topLeftCorner(n1, n1), Linv.topLeftCorner(n1, n1), scratch);
  invert_lower_rec_buf(L.bottomRightCorner(n2, n2), Linv.bottomRightCorner(n2, n2), scratch);
  Linv.topRightCorner(n1, n2).setZero();
  auto tmp = scratch.topLeftCorner(n2, n1);
  tmp.noalias() = -L.bottomLeftCorner(n2, n1) * Linv.topLeftCorner(n1, n1).triangularView<Eigen::Lower>();
  Linv.bottomLeftCorner(n2, n1).noalias() =
    Linv.bottomRightCorner(n2, n2).triangularView<Eigen::Lower>() * tmp;
}

// Fast unrolled Cholesky + lower-triangular inverse R = L^{-1} for small q x q blocks (q = 1..4)
inline bool factor_invert_small_spd(
    const Eigen::MatrixXd& A,
    Eigen::Ref<Eigen::MatrixXd> R_out,
    Eigen::LLT<Eigen::MatrixXd>& llt_fallback,
    double& logdet_acc
) {
  const int q = A.rows();
  if (q == 1) {
    const double a00 = A(0, 0);
    if (!(a00 > 0.0)) return false;
    logdet_acc += std::log(a00);
    R_out(0, 0) = 1.0 / std::sqrt(a00);
    return true;
  }
  if (q == 2) {
    const double a00 = A(0, 0), a10 = A(1, 0), a11 = A(1, 1);
    if (!(a00 > 0.0)) return false;
    const double l00 = std::sqrt(a00);
    const double r00 = 1.0 / l00;
    const double l10 = a10 * r00;
    const double d11 = a11 - l10 * l10;
    if (!(d11 > 0.0)) return false;
    const double l11 = std::sqrt(d11);
    const double r11 = 1.0 / l11;
    logdet_acc += 2.0 * (std::log(l00) + std::log(l11));
    R_out(0, 0) = r00; R_out(0, 1) = 0.0;
    R_out(1, 0) = -l10 * r00 * r11; R_out(1, 1) = r11;
    return true;
  }
  if (q <= 4) {
    double L[4][4] = {{0.0}};
    double R[4][4] = {{0.0}};
    double ld = 0.0;
    for (int j = 0; j < q; ++j) {
      double d = A(j, j);
      for (int k = 0; k < j; ++k) d -= L[j][k] * L[j][k];
      if (!(d > 0.0)) return false;
      const double ljj = std::sqrt(d);
      L[j][j] = ljj;
      ld += 2.0 * std::log(ljj);
      const double inv_ljj = 1.0 / ljj;
      R[j][j] = inv_ljj;
      for (int i = j + 1; i < q; ++i) {
        double s = A(i, j);
        for (int k = 0; k < j; ++k) s -= L[i][k] * L[j][k];
        L[i][j] = s * inv_ljj;
      }
    }
    for (int j = 0; j < q; ++j) {
      for (int i = j + 1; i < q; ++i) {
        double s = 0.0;
        for (int k = j; k < i; ++k) s -= L[i][k] * R[k][j];
        R[i][j] = s * R[i][i];
      }
    }
    for (int c = 0; c < q; ++c) {
      for (int r = 0; r < q; ++r) {
        R_out(r, c) = (r >= c) ? R[r][c] : 0.0;
      }
    }
    logdet_acc += ld;
    return true;
  }
  llt_fallback.compute(A);
  if (llt_fallback.info() != Eigen::Success) return false;
  const auto L_ass = llt_fallback.matrixL();
  for (int j = 0; j < q; ++j) logdet_acc += 2.0 * std::log(L_ass(j, j));
  R_out.setIdentity();
  L_ass.solveInPlace(R_out);
  return true;
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

  // Exact MAP update in theta_j = log(tau_j) coordinates under Half-Student-t(df = 3, scale = prior_tau_sd)
  const double df_t = 3.0;
  const double df_s0_sq = df_t * prior_tau_sd * prior_tau_sd;
  const double a_quad = static_cast<double>(n_clusters) + df_t;
  Eigen::VectorXd new_sd(q);
  for (int j = 0; j < q; ++j) {
    const double s_jj = std::max(S(j, j), 1e-8);
    const double b_quad = static_cast<double>(n_clusters - 1) * df_s0_sq - s_jj;
    const double c_pos = df_s0_sq * s_jj;
    const double disc = std::sqrt(b_quad * b_quad + 4.0 * a_quad * c_pos);
    const double v_hat = (b_quad >= 0.0)
      ? (2.0 * c_pos) / (b_quad + disc)
      : (-b_quad + disc) / (2.0 * a_quad);
    new_sd(j) = std::sqrt(std::max(v_hat, 1e-6));
  }

  if (corr_mode == 0 || q == 1) {
    return (new_sd.array().square()).matrix().asDiagonal();
  }

  Eigen::VectorXd emp_sd = S.diagonal().cwiseMax(1e-10).cwiseSqrt();
  Eigen::MatrixXd R = Eigen::MatrixXd::Identity(q, q);
  double shrink = static_cast<double>(n_clusters) /
    (static_cast<double>(n_clusters) + 2.0 * std::max(lkj_eta - 1.0, 0.0));
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
  double covariance_score_norm = NA_REAL;
  double log_lik_laplace = NA_REAL;
  double final_outer_diff = 1.0;
  double final_inner_step = 1.0;
  int outer_used = 0;
  bool singular = false;
  bool strictly_converged = false;
  std::string status = "engine_error";
};

struct ScheduleGroup {
  std::vector<int> subjs;
  std::vector<int> items;
  std::vector<int> cells;
  Eigen::MatrixXd H;
  Eigen::MatrixXd Q;
  Eigen::MatrixXd B;
  Eigen::MatrixXd W_grp;
  Eigen::MatrixXd C_bs_grp;
  Eigen::MatrixXd T_2b;
  Eigen::VectorXd x_grp;
  Eigen::VectorXd v_grp;
};

struct GLMMWorkspace {
  int N = 0, P = 0, q_s = 0, q_i = 0, q_t = 0;
  int n_subj = 0, n_item = 0, n_testlet = 0;
  bool has_testlet = false;
  bool scalar_two_way = false;
  bool unique_subj_item = true;
  bool use_archetypes = false;
  bool use_group_rank1 = false;
  bool use_sparse_items = false;
  int n_arch = 0;
  int dim_s = 0, dim_i = 0, dim_t = 0, dim_2 = 0;

  // Contiguous subject-major CSR permutation and indices
  std::vector<int> perm;
  std::vector<int> subj_ptr;
  Eigen::VectorXi y_csr;
  Eigen::VectorXi subj_idx, item_idx, testlet_idx;
  Eigen::VectorXi cell_arch_idx;

  // Unique cell-archetype tables (M_u <= max_cell_archetypes)
  RowMatrixXd X_u, Zs_u, Zi_u, Zt_u;
  RowMatrixXd ZsL_u, ZiL_u, ZtL_u;
  Eigen::VectorXd Xb_u, w_arch, g_arch;
  Eigen::MatrixXd w_sc_buf, g_sc_buf, w_ic_buf, g_ic_buf;
  std::vector<Eigen::MatrixXd> ZZs_u, ZZi_u;

  // Schedule groups for Rank-1 Schur reduction & sparse subject-item lists
  std::vector<ScheduleGroup> groups;
  std::vector<std::vector<int>> subj_u_items;

  RowMatrixXd X_rm, Zs_rm, Zi_rm, Zt_rm;
  RowMatrixXd ZsL, ZiL, ZtL;

  Eigen::VectorXd w_vec, wp_vec, g_vec;
  Eigen::VectorXd Xb_vec, eta0_vec, deta_vec;

  // Flat contiguous block slabs replacing std::vector<Eigen::MatrixXd>
  Eigen::MatrixXd R_ass_slab, V_ss_slab, V_ii_slab, V_tt_slab;
  Eigen::MatrixXd Ass_tmp, V_tmp_s, V_tmp_i, V_tmp_t;
  Eigen::MatrixXd L_22_saved, S_22, G_mat, C_bs, C_2b_buf, S_2b_saved, S22_inv_S2b_saved, H_bb_buf;
  Eigen::MatrixXd H_bb_schur, L22_inv, inv_scratch_buf, W_mat, M_cross, M_it;
  Eigen::MatrixXd sum_uu_s_buf, sum_uu_i_buf, sum_uu_t_buf;
  Eigen::VectorXd c_lap, cs_vec, c2_vec;
  Eigen::VectorXd rhs_s_buf, rhs_2_buf, rhs_b_buf, rs_tilde_buf, cs_tilde_buf;
  Eigen::VectorXd r_2_tilde_buf, S22_inv_r2_buf, S22_inv_c2_buf, r_b_tilde_buf;
  Eigen::VectorXd d_beta_buf, d_v2_buf, d_vs_buf, b_try_buf, step_2_buf, step_s_buf;
  Eigen::VectorXd us_buf, ui_buf, ut_buf;
  Eigen::LLT<Eigen::MatrixXd> llt_s_ws, llt_22_ws, llt_22_fin_ws, llt_22_grad_ws;
  Eigen::LDLT<Eigen::MatrixXd> ldlt_b_ws;

  GLMMWorkspace(
      const Eigen::MatrixXd& X,
      const Eigen::VectorXi& s_idx_in,
      const Eigen::MatrixXd& Z_subj,
      int ns,
      const Eigen::VectorXi& i_idx_in,
      const Eigen::MatrixXd& Z_item,
      int ni,
      const Eigen::VectorXi& t_idx_in,
      const Eigen::MatrixXd& Z_testlet,
      int nt
  ) : N(X.rows()), P(X.cols()),
      q_s(Z_subj.cols()), q_i(Z_item.cols()),
      n_subj(ns), n_item(ni), n_testlet(nt),
      has_testlet(nt > 0 && Z_testlet.cols() > 0),
      llt_s_ws(Z_subj.cols()),
      ldlt_b_ws(X.cols())
  {
    q_t = has_testlet ? Z_testlet.cols() : 0;
    scalar_two_way = (!has_testlet && q_s == 1 && q_i == 1);
    dim_s = n_subj * q_s;
    dim_i = n_item * q_i;
    dim_t = has_testlet ? (n_testlet * q_t) : 0;
    dim_2 = dim_i + dim_t;

    llt_22_ws = Eigen::LLT<Eigen::MatrixXd>(dim_2);
    llt_22_fin_ws = Eigen::LLT<Eigen::MatrixXd>(dim_2);
    llt_22_grad_ws = Eigen::LLT<Eigen::MatrixXd>(dim_2);

    // Build subject-major CSR permutation sorted by (subj_idx, item_idx, orig_idx)
    perm.resize(N);
    std::iota(perm.begin(), perm.end(), 0);
    std::sort(perm.begin(), perm.end(), [&](int a, int b) {
      if (s_idx_in[a] != s_idx_in[b]) return s_idx_in[a] < s_idx_in[b];
      if (i_idx_in[a] != i_idx_in[b]) return i_idx_in[a] < i_idx_in[b];
      return a < b;
    });

    subj_ptr.assign(n_subj + 1, 0);
    y_csr.resize(N);
    subj_idx.resize(N);
    item_idx.resize(N);
    if (has_testlet) testlet_idx.resize(N);
    X_rm.resize(N, P);
    Zs_rm.resize(N, q_s);
    Zi_rm.resize(N, q_i);
    if (has_testlet) Zt_rm.resize(N, q_t);

    for (int k = 0; k < N; ++k) {
      const int orig = perm[k];
      const int s = s_idx_in[orig];
      subj_ptr[s + 1]++;
      subj_idx[k] = s;
      item_idx[k] = i_idx_in[orig];
      if (has_testlet) testlet_idx[k] = t_idx_in[orig];
      X_rm.row(k) = X.row(orig);
      Zs_rm.row(k) = Z_subj.row(orig);
      Zi_rm.row(k) = Z_item.row(orig);
      if (has_testlet) Zt_rm.row(k) = Z_testlet.row(orig);
    }
    for (int s = 0; s < n_subj; ++s) {
      subj_ptr[s + 1] += subj_ptr[s];
    }

    // Check unique (subj, item) pairs and collect per-subject unique sorted items
    subj_u_items.resize(n_subj);
    int max_u_items = 0;
    unique_subj_item = true;
    for (int s = 0; s < n_subj; ++s) {
      auto& ui = subj_u_items[s];
      ui.reserve(subj_ptr[s + 1] - subj_ptr[s]);
      int prev_i = -1;
      for (int k = subj_ptr[s]; k < subj_ptr[s + 1]; ++k) {
        const int ii = item_idx[k];
        if (ii == prev_i) {
          unique_subj_item = false;
        } else {
          ui.push_back(ii);
          prev_i = ii;
        }
      }
      if (static_cast<int>(ui.size()) > max_u_items) {
        max_u_items = static_cast<int>(ui.size());
      }
    }

    // Deduplicate (X, Z_s, Z_i, Z_t) rows into compact cell archetypes when <= max_cell_archetypes
    cell_arch_idx.resize(N);
    std::vector<int> rep_rows;
    rep_rows.reserve(max_cell_archetypes + 1);
    bool arch_ok = true;
    for (int k = 0; k < N; ++k) {
      int found = -1;
      for (int c = 0; c < static_cast<int>(rep_rows.size()); ++c) {
        const int r0 = rep_rows[c];
        if (X_rm.row(k) == X_rm.row(r0) &&
            Zs_rm.row(k) == Zs_rm.row(r0) &&
            Zi_rm.row(k) == Zi_rm.row(r0) &&
            (!has_testlet || Zt_rm.row(k) == Zt_rm.row(r0))) {
          found = c;
          break;
        }
      }
      if (found < 0) {
        if (static_cast<int>(rep_rows.size()) >= max_cell_archetypes) {
          arch_ok = false;
          break;
        }
        found = static_cast<int>(rep_rows.size());
        rep_rows.push_back(k);
      }
      cell_arch_idx[k] = found;
    }

    if (arch_ok && !rep_rows.empty()) {
      use_archetypes = true;
      n_arch = static_cast<int>(rep_rows.size());
      X_u.resize(n_arch, P);
      Zs_u.resize(n_arch, q_s);
      Zi_u.resize(n_arch, q_i);
      if (has_testlet) Zt_u.resize(n_arch, q_t);
      for (int c = 0; c < n_arch; ++c) {
        const int r0 = rep_rows[c];
        X_u.row(c) = X_rm.row(r0);
        Zs_u.row(c) = Zs_rm.row(r0);
        Zi_u.row(c) = Zi_rm.row(r0);
        if (has_testlet) Zt_u.row(c) = Zt_rm.row(r0);
      }
      ZsL_u.resize(n_arch, q_s);
      ZiL_u.resize(n_arch, q_i);
      if (has_testlet) ZtL_u.resize(n_arch, q_t);
      Xb_u.resize(n_arch);
      w_arch.resize(n_arch);
      g_arch.resize(n_arch);
      w_sc_buf.resize(n_arch, n_subj);
      g_sc_buf.resize(n_arch, n_subj);
      w_ic_buf.resize(n_arch, n_item);
      g_ic_buf.resize(n_arch, n_item);
      ZZs_u.assign(n_arch, Eigen::MatrixXd::Zero(q_s, q_s));
      ZZi_u.assign(n_arch, Eigen::MatrixXd::Zero(q_i, q_i));
    }

    // Detect identical item-cell schedule groups for Rank-1 Schur reduction when q_i >= 2
    if (use_archetypes && !has_testlet && unique_subj_item && q_i >= 2 && n_subj >= 4) {
      std::vector<ScheduleGroup> cand_groups;
      for (int s = 0; s < n_subj; ++s) {
        const int beg = subj_ptr[s];
        const int len = subj_ptr[s + 1] - beg;
        if (len == 0) continue;
        int match_g = -1;
        for (int g = 0; g < static_cast<int>(cand_groups.size()); ++g) {
          if (static_cast<int>(cand_groups[g].items.size()) != len) continue;
          bool same = true;
          for (int r = 0; r < len; ++r) {
            if (cand_groups[g].items[r] != item_idx[beg + r] ||
                cand_groups[g].cells[r] != cell_arch_idx[beg + r]) {
              same = false;
              break;
            }
          }
          if (same) { match_g = g; break; }
        }
        if (match_g < 0) {
          ScheduleGroup ng;
          ng.subjs.push_back(s);
          ng.items.resize(len);
          ng.cells.resize(len);
          for (int r = 0; r < len; ++r) {
            ng.items[r] = item_idx[beg + r];
            ng.cells[r] = cell_arch_idx[beg + r];
          }
          cand_groups.push_back(std::move(ng));
          if (static_cast<int>(cand_groups.size()) > std::max(8, n_subj / 2)) break;
        } else {
          cand_groups[match_g].subjs.push_back(s);
        }
      }
      int covered = 0;
      for (const auto& g : cand_groups) covered += static_cast<int>(g.subjs.size());
      int non_empty_subj = 0;
      for (int s = 0; s < n_subj; ++s) if (subj_ptr[s + 1] > subj_ptr[s]) non_empty_subj++;
      if (covered == non_empty_subj && !cand_groups.empty() &&
          static_cast<int>(cand_groups.size()) <= std::max(8, n_subj / 2)) {
        use_group_rank1 = true;
        groups = std::move(cand_groups);
        for (auto& g : groups) {
          const int mg = static_cast<int>(g.items.size());
          const int cols_g = static_cast<int>(g.subjs.size()) * q_s;
          g.H.resize(mg, cols_g);
          g.Q.resize(mg, mg);
          g.B.resize(dim_2, mg);
          g.W_grp.resize(dim_2, cols_g);
          g.C_bs_grp.resize(P, cols_g);
          g.T_2b.resize(mg, P);
          g.x_grp.resize(cols_g);
          g.v_grp.resize(mg);
        }
      }
    }

    if (!use_group_rank1 && !has_testlet && n_item >= 8 && max_u_items <= (n_item * 3) / 5) {
      use_sparse_items = true;
    }

    if (!use_archetypes) {
      ZsL.resize(N, q_s);
      ZiL.resize(N, q_i);
      ZtL.resize(has_testlet ? N : 0, has_testlet ? q_t : 0);
    }

    w_vec.resize(N);
    wp_vec.resize(N);
    g_vec.resize(N);
    if (!use_archetypes) Xb_vec.resize(N);
    eta0_vec.resize(N);
    deta_vec.resize(N);

    R_ass_slab = Eigen::MatrixXd::Zero(q_s, std::max(dim_s, q_s));
    V_ss_slab = Eigen::MatrixXd::Zero(q_s, std::max(dim_s, q_s));
    V_ii_slab = Eigen::MatrixXd::Zero(q_i, std::max(dim_i, q_i));
    if (has_testlet) V_tt_slab = Eigen::MatrixXd::Zero(q_t, std::max(dim_t, q_t));
    Ass_tmp.resize(q_s, q_s);
    V_tmp_s.resize(q_s, q_s);
    V_tmp_i.resize(q_i, q_i);
    if (has_testlet) V_tmp_t.resize(q_t, q_t);

    L_22_saved = Eigen::MatrixXd::Identity(dim_2, dim_2);
    S_22 = Eigen::MatrixXd::Zero(dim_2, dim_2);
    if (!use_group_rank1) G_mat = Eigen::MatrixXd::Zero(dim_2, dim_s);
    C_bs = Eigen::MatrixXd::Zero(P, dim_s);
    C_2b_buf = Eigen::MatrixXd::Zero(dim_2, P);
    S_2b_saved.resize(dim_2, P);
    S22_inv_S2b_saved.resize(dim_2, P);
    H_bb_buf.resize(P, P);
    H_bb_schur = Eigen::MatrixXd::Identity(P, P);
    L22_inv = Eigen::MatrixXd::Zero(dim_2, dim_2);
    inv_scratch_buf.resize(dim_2, (dim_2 / 2) + 8);
    W_mat = Eigen::MatrixXd::Zero(dim_2, dim_s);
    if (!use_group_rank1) M_cross.resize(dim_2, dim_s);
    if (has_testlet) M_it.resize(dim_i, dim_t);
    sum_uu_s_buf.resize(q_s, q_s);
    sum_uu_i_buf.resize(q_i, q_i);
    if (has_testlet) sum_uu_t_buf.resize(q_t, q_t);

    c_lap.resize(N);
    cs_vec.resize(dim_s);
    c2_vec.resize(dim_2);
    rhs_s_buf.resize(dim_s);
    rhs_2_buf.resize(dim_2);
    rhs_b_buf.resize(P);
    rs_tilde_buf.resize(dim_s);
    cs_tilde_buf.resize(dim_s);
    r_2_tilde_buf.resize(dim_2);
    S22_inv_r2_buf.resize(dim_2);
    S22_inv_c2_buf.resize(dim_2);
    r_b_tilde_buf.resize(P);
    d_beta_buf.resize(P);
    d_v2_buf.resize(dim_2);
    d_vs_buf.resize(dim_s);
    b_try_buf.resize(P);
    step_2_buf.resize(dim_2);
    step_s_buf.resize(dim_s);
    us_buf.resize(q_s);
    ui_buf.resize(q_i);
    if (has_testlet) ut_buf.resize(q_t);
  }

  void update_ZL(const Eigen::MatrixXd& L_s, const Eigen::MatrixXd& L_i, const Eigen::MatrixXd& L_t, int mask = 7) {
    if (use_archetypes) {
      if (mask & 1) {
        ZsL_u.noalias() = Zs_u * L_s;
        for (int c = 0; c < n_arch; ++c) {
          ZZs_u[c].noalias() = ZsL_u.row(c).transpose() * ZsL_u.row(c);
        }
      }
      if (mask & 2) {
        ZiL_u.noalias() = Zi_u * L_i;
        for (int c = 0; c < n_arch; ++c) {
          ZZi_u[c].noalias() = ZiL_u.row(c).transpose() * ZiL_u.row(c);
        }
      }
      if (has_testlet && (mask & 4)) ZtL_u.noalias() = Zt_u * L_t;
    } else {
      if (mask & 1) ZsL.noalias() = Zs_rm * L_s;
      if (mask & 2) ZiL.noalias() = Zi_rm * L_i;
      if (has_testlet && (mask & 4)) ZtL.noalias() = Zt_rm * L_t;
    }
  }

  void apply_R_ass(const Eigen::VectorXd& x, Eigen::VectorXd& out) const {
    if (q_s == 1) {
      const double* rp = R_ass_slab.data();
      for (int s = 0; s < n_subj; ++s) out[s] = rp[s] * x[s];
      return;
    }
    for (int s = 0; s < n_subj; ++s) {
      const double* rp = R_ass_slab.data() + static_cast<std::ptrdiff_t>(s * q_s) * q_s;
      const double* xp = x.data() + s * q_s;
      double* op = out.data() + s * q_s;
      for (int a = 0; a < q_s; ++a) {
        double acc = 0.0;
        for (int b = 0; b <= a; ++b) acc += rp[a + b * q_s] * xp[b];
        op[a] = acc;
      }
    }
  }

  void apply_Rt_ass_in_place(Eigen::VectorXd& x) const {
    if (q_s == 1) {
      const double* rp = R_ass_slab.data();
      for (int s = 0; s < n_subj; ++s) x[s] *= rp[s];
      return;
    }
    for (int s = 0; s < n_subj; ++s) {
      const double* rp = R_ass_slab.data() + static_cast<std::ptrdiff_t>(s * q_s) * q_s;
      double* xp = x.data() + s * q_s;
      for (int a = 0; a < q_s; ++a) {
        double acc = rp[a + a * q_s] * xp[a];
        for (int b = a + 1; b < q_s; ++b) acc += rp[b + a * q_s] * xp[b];
        xp[a] = acc;
      }
    }
  }

  // Compute out_2 = rhs_2 - G * x_s
  void sub_G_mul_vec(const Eigen::VectorXd& rhs_2, const Eigen::VectorXd& x_s, Eigen::VectorXd& out_2) {
    if (!use_group_rank1) {
      out_2.noalias() = rhs_2 - G_mat * x_s;
      return;
    }
    out_2 = rhs_2;
    for (auto& g : groups) {
      const int ns_g = static_cast<int>(g.subjs.size());
      for (int k = 0; k < ns_g; ++k) {
        g.x_grp.segment(k * q_s, q_s) = x_s.segment(g.subjs[k] * q_s, q_s);
      }
      g.v_grp.noalias() = g.H * g.x_grp;
      const int mg = static_cast<int>(g.items.size());
      for (int r = 0; r < mg; ++r) {
        const double vr = g.v_grp[r];
        const double* zi = ZiL_u.data() + static_cast<std::ptrdiff_t>(g.cells[r]) * q_i;
        double* dst = out_2.data() + g.items[r] * q_i;
        for (int b = 0; b < q_i; ++b) dst[b] -= vr * zi[b];
      }
    }
  }

  // Compute out_s = rhs_s - G^T * y_2 (and optionally - C_bs^T * d_beta)
  void sub_Gt_mul_vec(
      const Eigen::VectorXd& rhs_s,
      const Eigen::VectorXd& y_2,
      const Eigen::VectorXd* d_beta,
      Eigen::VectorXd& out_s
  ) {
    if (!use_group_rank1) {
      if (d_beta != nullptr) {
        out_s.noalias() = rhs_s - G_mat.transpose() * y_2 - C_bs.transpose() * (*d_beta);
      } else {
        out_s.noalias() = rhs_s - G_mat.transpose() * y_2;
      }
      return;
    }
    if (d_beta != nullptr) {
      out_s.noalias() = rhs_s - C_bs.transpose() * (*d_beta);
    } else {
      out_s = rhs_s;
    }
    for (auto& g : groups) {
      const int mg = static_cast<int>(g.items.size());
      for (int r = 0; r < mg; ++r) {
        const double* zi = ZiL_u.data() + static_cast<std::ptrdiff_t>(g.cells[r]) * q_i;
        const double* y2p = y_2.data() + g.items[r] * q_i;
        double acc = 0.0;
        for (int b = 0; b < q_i; ++b) acc += zi[b] * y2p[b];
        g.v_grp[r] = acc;
      }
      g.x_grp.noalias() = g.H.transpose() * g.v_grp;
      const int ns_g = static_cast<int>(g.subjs.size());
      for (int k = 0; k < ns_g; ++k) {
        out_s.segment(g.subjs[k] * q_s, q_s) -= g.x_grp.segment(k * q_s, q_s);
      }
    }
  }

  void compute_deta_vec(
      const Eigen::VectorXd* d_beta,
      const Eigen::VectorXd& d_vs,
      const Eigen::VectorXd& d_v2
  ) {
    if (use_archetypes) {
      if (d_beta != nullptr) {
        Xb_u.noalias() = X_u * (*d_beta);
      } else {
        Xb_u.setZero();
      }
      const int* c_ptr = cell_arch_idx.data();
      const int* i_ptr = item_idx.data();
      double* de_ptr = deta_vec.data();

      if (scalar_two_way) {
        const double ls0 = ZsL_u(0, 0);
        const double li0 = ZiL_u(0, 0);
        const double* xb_p = Xb_u.data();
        const double* dv2_p = d_v2.data();
        for (int s = 0; s < n_subj; ++s) {
          const double vs_term = ls0 * d_vs[s];
          for (int k = subj_ptr[s]; k < subj_ptr[s + 1]; ++k) {
            de_ptr[k] = xb_p[c_ptr[k]] + vs_term + li0 * dv2_p[i_ptr[k]];
          }
        }
        return;
      }

      for (int s = 0; s < n_subj; ++s) {
        const double* dvs = d_vs.data() + s * q_s;
        for (int k = subj_ptr[s]; k < subj_ptr[s + 1]; ++k) {
          const int c = c_ptr[k];
          const double* zs = ZsL_u.data() + static_cast<std::ptrdiff_t>(c) * q_s;
          const double* zi = ZiL_u.data() + static_cast<std::ptrdiff_t>(c) * q_i;
          const double* dv2_i = d_v2.data() + i_ptr[k] * q_i;
          double acc = Xb_u[c];
          for (int a = 0; a < q_s; ++a) acc += zs[a] * dvs[a];
          for (int a = 0; a < q_i; ++a) acc += zi[a] * dv2_i[a];
          if (has_testlet) {
            const double* zt = ZtL_u.data() + static_cast<std::ptrdiff_t>(c) * q_t;
            const double* dv2_t = d_v2.data() + dim_i + testlet_idx[k] * q_t;
            for (int a = 0; a < q_t; ++a) acc += zt[a] * dv2_t[a];
          }
          de_ptr[k] = acc;
        }
      }
      return;
    }

    if (d_beta != nullptr) {
      deta_vec.noalias() = X_rm * (*d_beta);
    } else {
      deta_vec.setZero();
    }
    const double* zs_ptr = ZsL.data();
    const double* zi_ptr = ZiL.data();
    const double* zt_ptr = has_testlet ? ZtL.data() : nullptr;
    for (int s = 0; s < n_subj; ++s) {
      const double* dvs = d_vs.data() + s * q_s;
      for (int k = subj_ptr[s]; k < subj_ptr[s + 1]; ++k) {
        const double* dv2_i = d_v2.data() + item_idx[k] * q_i;
        double acc = 0.0;
        for (int a = 0; a < q_s; ++a) acc += zs_ptr[a] * dvs[a];
        for (int a = 0; a < q_i; ++a) acc += zi_ptr[a] * dv2_i[a];
        zs_ptr += q_s;
        zi_ptr += q_i;
        if (has_testlet) {
          const double* dv2_t = d_v2.data() + dim_i + testlet_idx[k] * q_t;
          for (int a = 0; a < q_t; ++a) acc += zt_ptr[a] * dv2_t[a];
          zt_ptr += q_t;
        }
        deta_vec[k] += acc;
      }
    }
  }

  double eval_data_ll_step(double scale, double chance) const {
    double ll_data = 0.0;
    const int* y_ptr = y_csr.data();
    const double* e0_ptr = eta0_vec.data();
    const double* de_ptr = deta_vec.data();
    if (chance <= 0.0) {
      for (int k = 0; k < N; ++k) {
        double eta = e0_ptr[k] + scale * de_ptr[k];
        eta = std::min(std::max(eta, -25.0), 25.0);
        ll_data += (y_ptr[k] == 1) ? -softplus(-eta) : -softplus(eta);
      }
      return ll_data;
    }
    for (int k = 0; k < N; ++k) {
      double eta = e0_ptr[k] + scale * de_ptr[k];
      eta = std::min(std::max(eta, -25.0), 25.0);
      double p = inv_logit(eta);
      double mu = chance + (1.0 - chance) * p;
      mu = std::min(std::max(mu, chance + 1e-12), 1.0 - 1e-12);
      ll_data += (y_ptr[k] == 1) ? std::log(mu) : std::log1p(-mu);
    }
    return ll_data;
  }

  void invert_L22_lower() {
    invert_lower_rec_buf(L_22_saved, L22_inv, inv_scratch_buf);
  }

  void compute_W_mat_from_L22_inv() {
    if (use_group_rank1) {
      for (auto& g : groups) {
        const int mg = static_cast<int>(g.items.size());
        for (int r = 0; r < mg; ++r) {
          const int j = g.items[r];
          const int off_j = j * q_i;
          if (off_j > 0) g.B.col(r).head(off_j).setZero();
          const Eigen::RowVectorXd zi = ZiL_u.row(g.cells[r]);
          g.B.col(r).tail(dim_2 - off_j).noalias() =
            L22_inv.block(off_j, off_j, dim_2 - off_j, q_i) * zi.transpose();
        }
        g.W_grp.noalias() = g.B * g.H;
        const int ns_g = static_cast<int>(g.subjs.size());
        for (int k = 0; k < ns_g; ++k) {
          W_mat.middleCols(g.subjs[k] * q_s, q_s) = g.W_grp.middleCols(k * q_s, q_s);
        }
      }
      return;
    }

    if (use_sparse_items) {
      W_mat.setZero();
      for (int s = 0; s < n_subj; ++s) {
        const int off_s = s * q_s;
        auto Ws = W_mat.middleCols(off_s, q_s);
        for (int ii : subj_u_items[s]) {
          const int off_i = ii * q_i;
          Ws.bottomRows(dim_2 - off_i).noalias() +=
            L22_inv.block(off_i, off_i, dim_2 - off_i, q_i) *
            G_mat.block(off_i, off_s, q_i, q_s);
        }
      }
      return;
    }

    W_mat.noalias() = L22_inv.triangularView<Eigen::Lower>() * G_mat;
  }

  void compute_c_lap_fast() {
    invert_L22_lower();
    compute_W_mat_from_L22_inv();

    for (int s = 0; s < n_subj; ++s) {
      const int off_s = s * q_s;
      const auto Ws = W_mat.middleCols(off_s, q_s);
      const auto Rs = R_ass_slab.middleCols(off_s, q_s);
      V_tmp_s.noalias() = Ws.transpose() * Ws;
      for (int j = 0; j < q_s; ++j) V_tmp_s(j, j) += 1.0;
      V_ss_slab.middleCols(off_s, q_s).noalias() = Rs.transpose() * V_tmp_s * Rs;
    }
    for (int i = 0; i < n_item; ++i) {
      const int off_i = i * q_i;
      const auto Li_sub = L22_inv.block(off_i, off_i, dim_2 - off_i, q_i);
      V_ii_slab.middleCols(off_i, q_i).noalias() = Li_sub.transpose() * Li_sub;
    }
    if (has_testlet) {
      for (int t = 0; t < n_testlet; ++t) {
        const int off_t_rel = t * q_t;
        const int off_t = dim_i + off_t_rel;
        const auto Lt_sub = L22_inv.block(off_t, off_t, dim_2 - off_t, q_t);
        V_tt_slab.middleCols(off_t_rel, q_t).noalias() = Lt_sub.transpose() * Lt_sub;
      }
      M_it.noalias() = L22_inv.block(dim_i, 0, dim_t, dim_i).transpose() *
        L22_inv.block(dim_i, dim_i, dim_t, dim_t);
    }

    if (use_group_rank1) {
      // Compute cross_is directly via B_g^T * W_grp_g (m_g x n_g*q_s) without forming dense M_cross!
      for (auto& g : groups) {
        const int mg = static_cast<int>(g.items.size());
        const int ns_g = static_cast<int>(g.subjs.size());
        Eigen::MatrixXd BtW(mg, ns_g * q_s);
        BtW.noalias() = g.B.transpose() * g.W_grp;

        for (int k_s = 0; k_s < ns_g; ++k_s) {
          const int s = g.subjs[k_s];
          const int off_s = s * q_s;
          const auto Vss = V_ss_slab.middleCols(off_s, q_s);
          const auto Rs = R_ass_slab.middleCols(off_s, q_s);
          const int beg = subj_ptr[s];

          for (int r = 0; r < mg; ++r) {
            const int k = beg + r;
            const int i = g.items[r];
            const int c = g.cells[r];
            const auto Vii = V_ii_slab.middleCols(i * q_i, q_i);
            const double* zs = ZsL_u.data() + static_cast<std::ptrdiff_t>(c) * q_s;
            const double* zi = ZiL_u.data() + static_cast<std::ptrdiff_t>(c) * q_i;

            double quad = 0.0;
            for (int a = 0; a < q_s; ++a) {
              double tmp = 0.0;
              for (int b = 0; b < q_s; ++b) tmp += Vss(b, a) * zs[b];
              quad += zs[a] * tmp;
            }
            for (int a = 0; a < q_i; ++a) {
              double tmp = 0.0;
              for (int b = 0; b < q_i; ++b) tmp += Vii(b, a) * zi[b];
              quad += zi[a] * tmp;
            }
            // R_s * zs (R_s lower triangular)
            double cross_is = 0.0;
            for (int b = 0; b < q_s; ++b) {
              double r_zs_b = 0.0;
              for (int a = 0; a <= b; ++a) r_zs_b += Rs(b, a) * zs[a];
              cross_is += BtW(r, k_s * q_s + b) * r_zs_b;
            }
            quad -= 2.0 * cross_is;
            c_lap[k] = -0.5 * std::max(quad, 0.0) * wp_vec[k];
          }
        }
      }
    } else if (use_sparse_items && use_archetypes && !has_testlet) {
      // Sparse cross-covariance evaluation without forming dense M_cross = L22_inv^T * W_mat
      Eigen::MatrixXd B_ic(dim_2, n_item * n_arch);
      for (int i = 0; i < n_item; ++i) {
        const int off_i = i * q_i;
        const auto Li_sub = L22_inv.block(off_i, off_i, dim_2 - off_i, q_i);
        for (int c = 0; c < n_arch; ++c) {
          Eigen::VectorXd zi_col = ZiL_u.row(c).transpose();
          B_ic.col(i * n_arch + c).tail(dim_2 - off_i).noalias() = Li_sub * zi_col;
        }
      }
      Eigen::MatrixXd W_sc(dim_2, n_arch);
      std::vector<unsigned char> seen_c(n_arch, 0);

      for (int s = 0; s < n_subj; ++s) {
        const int off_s = s * q_s;
        const auto Vss = V_ss_slab.middleCols(off_s, q_s);
        const auto Ws = W_mat.middleCols(off_s, q_s);
        const auto Rs = R_ass_slab.middleCols(off_s, q_s);
        std::fill(seen_c.begin(), seen_c.end(), 0);

        for (int k = subj_ptr[s]; k < subj_ptr[s + 1]; ++k) {
          const int i = item_idx[k];
          const int off_i = i * q_i;
          const int c = cell_arch_idx[k];
          const auto Vii = V_ii_slab.middleCols(off_i, q_i);
          const double* zs = ZsL_u.data() + static_cast<std::ptrdiff_t>(c) * q_s;
          const double* zi = ZiL_u.data() + static_cast<std::ptrdiff_t>(c) * q_i;

          if (!seen_c[c]) {
            for (int b = 0; b < q_s; ++b) {
              double acc = 0.0;
              for (int a = 0; a <= b; ++a) acc += Rs(b, a) * zs[a];
              us_buf[b] = acc;
            }
            W_sc.col(c).noalias() = Ws * us_buf;
            seen_c[c] = 1;
          }

          double quad = 0.0;
          for (int a = 0; a < q_s; ++a) {
            double tmp = 0.0;
            for (int b = 0; b < q_s; ++b) tmp += Vss(b, a) * zs[b];
            quad += zs[a] * tmp;
          }
          for (int a = 0; a < q_i; ++a) {
            double tmp = 0.0;
            for (int b = 0; b < q_i; ++b) tmp += Vii(b, a) * zi[b];
            quad += zi[a] * tmp;
          }
          const int rem = dim_2 - off_i;
          const double cross_is = B_ic.col(i * n_arch + c).tail(rem).dot(W_sc.col(c).tail(rem));
          quad -= 2.0 * cross_is;
          c_lap[k] = -0.5 * std::max(quad, 0.0) * wp_vec[k];
        }
      }
    } else {
      M_cross.noalias() = L22_inv.triangularView<Eigen::Lower>().transpose() * W_mat;
      for (int s = 0; s < n_subj; ++s) {
        const int off_s = s * q_s;
        const auto Rs = R_ass_slab.middleCols(off_s, q_s);
        for (int j = 0; j < q_s; ++j) {
          M_cross.col(off_s + j) *= Rs(j, j);
          for (int k = j + 1; k < q_s; ++k) {
            M_cross.col(off_s + j).noalias() += Rs(k, j) * M_cross.col(off_s + k);
          }
        }
      }

      const double* zs_ptr = use_archetypes ? nullptr : ZsL.data();
      const double* zi_ptr = use_archetypes ? nullptr : ZiL.data();
      const double* zt_ptr = (!use_archetypes && has_testlet) ? ZtL.data() : nullptr;

      for (int s = 0; s < n_subj; ++s) {
        const int off_s = s * q_s;
        const auto Vss = V_ss_slab.middleCols(off_s, q_s);
        for (int k = subj_ptr[s]; k < subj_ptr[s + 1]; ++k) {
          const int i = item_idx[k];
          const int off_i = i * q_i;
          const auto Vii = V_ii_slab.middleCols(off_i, q_i);
          const int c = use_archetypes ? cell_arch_idx[k] : 0;
          const double* zs = use_archetypes ? (ZsL_u.data() + static_cast<std::ptrdiff_t>(c) * q_s) : zs_ptr;
          const double* zi = use_archetypes ? (ZiL_u.data() + static_cast<std::ptrdiff_t>(c) * q_i) : zi_ptr;

          double quad = 0.0;
          for (int a = 0; a < q_s; ++a) {
            double tmp = 0.0;
            for (int b = 0; b < q_s; ++b) tmp += Vss(b, a) * zs[b];
            quad += zs[a] * tmp;
          }
          for (int a = 0; a < q_i; ++a) {
            double tmp = 0.0;
            for (int b = 0; b < q_i; ++b) tmp += Vii(b, a) * zi[b];
            quad += zi[a] * tmp;
          }
          double cross_is = 0.0;
          for (int a = 0; a < q_s; ++a) {
            double tmp = 0.0;
            for (int b = 0; b < q_i; ++b) tmp += M_cross(off_i + b, off_s + a) * zi[b];
            cross_is += tmp * zs[a];
          }
          quad -= 2.0 * cross_is;

          if (has_testlet) {
            const int t = testlet_idx[k];
            const int off_t_rel = t * q_t;
            const int off_t = dim_i + off_t_rel;
            const auto Vtt = V_tt_slab.middleCols(off_t_rel, q_t);
            const double* zt = use_archetypes ? (ZtL_u.data() + static_cast<std::ptrdiff_t>(c) * q_t) : zt_ptr;
            for (int a = 0; a < q_t; ++a) {
              double tmp = 0.0;
              for (int b = 0; b < q_t; ++b) tmp += Vtt(b, a) * zt[b];
              quad += zt[a] * tmp;
            }
            double cross_ts = 0.0;
            for (int a = 0; a < q_s; ++a) {
              double tmp = 0.0;
              for (int b = 0; b < q_t; ++b) tmp += M_cross(off_t + b, off_s + a) * zt[b];
              cross_ts += tmp * zs[a];
            }
            double cross_it = 0.0;
            for (int a = 0; a < q_t; ++a) {
              double tmp = 0.0;
              for (int b = 0; b < q_i; ++b) tmp += M_it(off_i + b, off_t_rel + a) * zi[b];
              cross_it += tmp * zt[a];
            }
            quad += 2.0 * (cross_it - cross_ts);
            if (!use_archetypes) zt_ptr += q_t;
          }

          c_lap[k] = -0.5 * std::max(quad, 0.0) * wp_vec[k];
          if (!use_archetypes) {
            zs_ptr += q_s;
            zi_ptr += q_i;
          }
        }
      }
    }

    cs_vec.setZero();
    c2_vec.setZero();
    if (use_archetypes && !has_testlet) {
      w_sc_buf.setZero();
      w_ic_buf.setZero();
      for (int s = 0; s < n_subj; ++s) {
        double* cs_c = w_sc_buf.data() + static_cast<std::ptrdiff_t>(s) * n_arch;
        for (int k = subj_ptr[s]; k < subj_ptr[s + 1]; ++k) {
          const double cn = c_lap[k];
          const int c = cell_arch_idx[k];
          cs_c[c] += cn;
          w_ic_buf(c, item_idx[k]) += cn;
        }
      }
      for (int s = 0; s < n_subj; ++s) {
        cs_vec.segment(s * q_s, q_s).noalias() = ZsL_u.transpose() * w_sc_buf.col(s);
      }
      for (int i = 0; i < n_item; ++i) {
        c2_vec.segment(i * q_i, q_i).noalias() = ZiL_u.transpose() * w_ic_buf.col(i);
      }
      return;
    }

    const double* zs_ptr = use_archetypes ? nullptr : ZsL.data();
    const double* zi_ptr = use_archetypes ? nullptr : ZiL.data();
    const double* zt_ptr = (!use_archetypes && has_testlet) ? ZtL.data() : nullptr;
    for (int s = 0; s < n_subj; ++s) {
      double* cs_p = cs_vec.data() + s * q_s;
      for (int k = subj_ptr[s]; k < subj_ptr[s + 1]; ++k) {
        const double cn = c_lap[k];
        const int c = use_archetypes ? cell_arch_idx[k] : 0;
        const double* zs = use_archetypes ? (ZsL_u.data() + static_cast<std::ptrdiff_t>(c) * q_s) : zs_ptr;
        const double* zi = use_archetypes ? (ZiL_u.data() + static_cast<std::ptrdiff_t>(c) * q_i) : zi_ptr;
        double* c2_p = c2_vec.data() + item_idx[k] * q_i;
        for (int a = 0; a < q_s; ++a) cs_p[a] += cn * zs[a];
        for (int a = 0; a < q_i; ++a) c2_p[a] += cn * zi[a];
        if (has_testlet) {
          const double* zt = use_archetypes ? (ZtL_u.data() + static_cast<std::ptrdiff_t>(c) * q_t) : zt_ptr;
          double* ct_p = c2_vec.data() + dim_i + testlet_idx[k] * q_t;
          for (int a = 0; a < q_t; ++a) ct_p[a] += cn * zt[a];
          if (!use_archetypes) zt_ptr += q_t;
        }
        if (!use_archetypes) {
          zs_ptr += q_s;
          zi_ptr += q_i;
        }
      }
    }
  }

  double eval_trials(
      const Eigen::VectorXd& beta,
      const Eigen::VectorXd& v_s,
      const Eigen::VectorXd& v_2,
      double chance,
      bool observed
  ) {
    // Step 1: Compute linear predictors eta0_vec in contiguous CSR order
    if (use_archetypes) {
      Xb_u.noalias() = X_u * beta;
      const int* c_ptr = cell_arch_idx.data();
      const int* i_ptr = item_idx.data();
      double* e0_ptr = eta0_vec.data();

      if (scalar_two_way) {
        const double ls0 = ZsL_u(0, 0);
        const double li0 = ZiL_u(0, 0);
        const double* xb_p = Xb_u.data();
        const double* v2_p = v_2.data();
        for (int s = 0; s < n_subj; ++s) {
          const double vs_term = ls0 * v_s[s];
          for (int k = subj_ptr[s]; k < subj_ptr[s + 1]; ++k) {
            e0_ptr[k] = xb_p[c_ptr[k]] + vs_term + li0 * v2_p[i_ptr[k]];
          }
        }
      } else {
        for (int s = 0; s < n_subj; ++s) {
          const double* vs_p = v_s.data() + s * q_s;
          for (int k = subj_ptr[s]; k < subj_ptr[s + 1]; ++k) {
            const int c = c_ptr[k];
            const double* zs = ZsL_u.data() + static_cast<std::ptrdiff_t>(c) * q_s;
            const double* zi = ZiL_u.data() + static_cast<std::ptrdiff_t>(c) * q_i;
            const double* v2_i = v_2.data() + i_ptr[k] * q_i;
            double eta_raw = Xb_u[c];
            for (int a = 0; a < q_s; ++a) eta_raw += zs[a] * vs_p[a];
            for (int a = 0; a < q_i; ++a) eta_raw += zi[a] * v2_i[a];
            if (has_testlet) {
              const double* zt = ZtL_u.data() + static_cast<std::ptrdiff_t>(c) * q_t;
              const double* v2_t = v_2.data() + dim_i + testlet_idx[k] * q_t;
              for (int a = 0; a < q_t; ++a) eta_raw += zt[a] * v2_t[a];
            }
            e0_ptr[k] = eta_raw;
          }
        }
      }
    } else {
      Xb_vec.noalias() = X_rm * beta;
      const double* zs_ptr = ZsL.data();
      const double* zi_ptr = ZiL.data();
      const double* zt_ptr = has_testlet ? ZtL.data() : nullptr;
      for (int s = 0; s < n_subj; ++s) {
        const double* vs_p = v_s.data() + s * q_s;
        for (int k = subj_ptr[s]; k < subj_ptr[s + 1]; ++k) {
          const double* v2_i = v_2.data() + item_idx[k] * q_i;
          double eta_raw = Xb_vec[k];
          for (int a = 0; a < q_s; ++a) eta_raw += zs_ptr[a] * vs_p[a];
          for (int a = 0; a < q_i; ++a) eta_raw += zi_ptr[a] * v2_i[a];
          zs_ptr += q_s;
          zi_ptr += q_i;
          if (has_testlet) {
            const double* v2_t = v_2.data() + dim_i + testlet_idx[k] * q_t;
            for (int a = 0; a < q_t; ++a) eta_raw += zt_ptr[a] * v2_t[a];
            zt_ptr += q_t;
          }
          eta0_vec[k] = eta_raw;
        }
      }
    }

    // Step 2: Evaluate likelihood, gradient, and curvature weights over contiguous arrays
    double ll_data = 0.0;
    const int* y_ptr = y_csr.data();
    const double* e0_ptr = eta0_vec.data();
    double* g_ptr = g_vec.data();
    double* w_ptr = w_vec.data();
    double* wp_ptr = wp_vec.data();

    if (chance <= 0.0) {
      for (int k = 0; k < N; ++k) {
        const double eta = std::min(std::max(e0_ptr[k], -25.0), 25.0);
        const double abs_eta = std::abs(eta);
        const double e = std::exp(-abs_eta);
        const double inv1pe = 1.0 / (1.0 + e);
        const double p = (eta >= 0.0) ? inv1pe : (e * inv1pe);
        const double omp = (eta >= 0.0) ? (e * inv1pe) : inv1pe;
        const double sp = std::log1p(e);
        const bool correct = (y_ptr[k] == 1);
        ll_data += correct ? ((eta >= 0.0) ? -sp : (eta - sp)) : ((eta >= 0.0) ? (-eta - sp) : -sp);
        g_ptr[k] = correct ? omp : -p;
        if (observed) {
          const double w = p * omp;
          w_ptr[k] = w;
          wp_ptr[k] = w * (1.0 - 2.0 * p);
        } else {
          w_ptr[k] = std::max(p * omp, 1e-9);
        }
      }
      return ll_data;
    }

    for (int k = 0; k < N; ++k) {
      const double eta = std::min(std::max(e0_ptr[k], -25.0), 25.0);
      const double p = inv_logit(eta);
      const double omp = 1.0 - p;
      const double mu = std::min(std::max(chance + (1.0 - chance) * p, chance + 1e-12), 1.0 - 1e-12);
      const bool correct = (y_ptr[k] == 1);
      ll_data += correct ? std::log(mu) : std::log1p(-mu);
      g_ptr[k] = correct ? ((1.0 - chance) * p / mu) * omp : -p;
      if (observed) {
        const double q = correct ? p / mu : 0.0;
        const double v = q * (1.0 - q);
        w_ptr[k] = p * omp - v;
        wp_ptr[k] = p * omp * (1.0 - 2.0 * p) - v * (1.0 - 2.0 * q);
      } else {
        if (!correct) {
          w_ptr[k] = std::max(p * omp, 1e-9);
        } else {
          const double q = p / mu;
          const double w_obs = p * omp - q * (1.0 - q);
          const double w_fs = ((1.0 - chance) * p * p * omp) / mu;
          w_ptr[k] = std::max({w_obs, w_fs, 1e-9});
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
    S_22.triangularView<Eigen::Lower>().setZero();
    S_22.diagonal().setOnes();

    double logdet_Ass = 0.0;

    if (use_archetypes && !has_testlet) {
      // Fast sufficient-weight aggregation across compact cell archetypes
      w_arch.setZero();
      g_arch.setZero();
      w_sc_buf.setZero();
      g_sc_buf.setZero();
      w_ic_buf.setZero();
      g_ic_buf.setZero();

      if (!use_group_rank1) {
        if (use_sparse_items) {
          for (int s = 0; s < n_subj; ++s) {
            const int off_s = s * q_s;
            for (int ii : subj_u_items[s]) {
              G_mat.block(ii * q_i, off_s, q_i, q_s).setZero();
            }
          }
        } else if (!(unique_subj_item && N == n_subj * n_item)) {
          G_mat.setZero();
        }
      }

      const int* c_ptr = cell_arch_idx.data();
      const int* i_ptr = item_idx.data();
      const double* w_ptr = w_vec.data();
      const double* g_ptr = g_vec.data();

      for (int s = 0; s < n_subj; ++s) {
        double* wsc = w_sc_buf.data() + static_cast<std::ptrdiff_t>(s) * n_arch;
        double* gsc = g_sc_buf.data() + static_cast<std::ptrdiff_t>(s) * n_arch;
        for (int k = subj_ptr[s]; k < subj_ptr[s + 1]; ++k) {
          const int c = c_ptr[k];
          const int ii = i_ptr[k];
          const double wn = w_ptr[k];
          const double gn = g_ptr[k];
          wsc[c] += wn;
          gsc[c] += gn;
          w_ic_buf(c, ii) += wn;
          g_ic_buf(c, ii) += gn;
        }
        for (int c = 0; c < n_arch; ++c) {
          w_arch[c] += wsc[c];
          g_arch[c] += gsc[c];
        }
      }

      H_bb_buf = P_beta;
      for (int c = 0; c < n_arch; ++c) {
        H_bb_buf.noalias() += w_arch[c] * (X_u.row(c).transpose() * X_u.row(c));
      }
      rhs_b.noalias() = X_u.transpose() * g_arch - P_beta * beta;

      for (int i = 0; i < n_item; ++i) {
        const int off_i = i * q_i;
        rhs_2.segment(off_i, q_i).noalias() =
          ZiL_u.transpose() * g_ic_buf.col(i) - v_2.segment(off_i, q_i);
        C_2b_buf.middleRows(off_i, q_i).noalias() =
          ZiL_u.transpose() * w_ic_buf.col(i).asDiagonal() * X_u;
        auto S_ii = S_22.block(off_i, off_i, q_i, q_i);
        for (int c = 0; c < n_arch; ++c) {
          const double wic = w_ic_buf(c, i);
          if (wic != 0.0) S_ii.noalias() += wic * ZZi_u[c];
        }
      }

      for (int s = 0; s < n_subj; ++s) {
        const int off_s = s * q_s;
        rhs_s.segment(off_s, q_s).noalias() =
          ZsL_u.transpose() * g_sc_buf.col(s) - v_s.segment(off_s, q_s);
        Ass_tmp.setIdentity();
        for (int c = 0; c < n_arch; ++c) {
          const double wsc = w_sc_buf(c, s);
          if (wsc != 0.0) Ass_tmp.noalias() += wsc * ZZs_u[c];
        }
        auto Rs = R_ass_slab.middleCols(off_s, q_s);
        if (!factor_invert_small_spd(Ass_tmp, Rs, llt_s_ws, logdet_Ass)) return false;

        // C_bs(:, off_s..off_s+q_s-1) = (X_u^T * diag(w_sc) * ZsL_u) * Rs^T
        V_tmp_s.resize(P, q_s); // reuse small stack-like matrix if P == q_s else compute directly
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> C_raw =
          X_u.transpose() * w_sc_buf.col(s).asDiagonal() * ZsL_u;
        C_bs.middleCols(off_s, q_s).noalias() = C_raw * Rs.transpose();
      }
      V_tmp_s.resize(q_s, q_s);

      if (use_group_rank1) {
        S_2b_saved = C_2b_buf;
        double u_sc[max_cell_archetypes][8]; // q_s <= 8 fast path; dynamic fallback if q_s > 8
        Eigen::MatrixXd U_dyn;
        if (q_s > 8) U_dyn.resize(q_s, n_arch);

        for (auto& g : groups) {
          const int mg = static_cast<int>(g.items.size());
          const int ns_g = static_cast<int>(g.subjs.size());
          for (int k_s = 0; k_s < ns_g; ++k_s) {
            const int s = g.subjs[k_s];
            const int off_s = s * q_s;
            const auto Rs = R_ass_slab.middleCols(off_s, q_s);
            g.C_bs_grp.middleCols(k_s * q_s, q_s) = C_bs.middleCols(off_s, q_s);

            if (q_s <= 8) {
              for (int c = 0; c < n_arch; ++c) {
                const double* zs = ZsL_u.data() + static_cast<std::ptrdiff_t>(c) * q_s;
                for (int a = 0; a < q_s; ++a) {
                  double acc = 0.0;
                  for (int b = 0; b <= a; ++b) acc += Rs(a, b) * zs[b];
                  u_sc[c][a] = acc;
                }
              }
              const int beg = subj_ptr[s];
              for (int r = 0; r < mg; ++r) {
                const double wn = w_ptr[beg + r];
                const int c = g.cells[r];
                for (int a = 0; a < q_s; ++a) {
                  g.H(r, k_s * q_s + a) = wn * u_sc[c][a];
                }
              }
            } else {
              U_dyn.noalias() = Rs * ZsL_u.transpose();
              const int beg = subj_ptr[s];
              for (int r = 0; r < mg; ++r) {
                const double wn = w_ptr[beg + r];
                g.H.block(r, k_s * q_s, 1, q_s) = wn * U_dyn.col(g.cells[r]).transpose();
              }
            }
          }

          g.Q.triangularView<Eigen::Lower>().setZero();
          g.Q.selfadjointView<Eigen::Lower>().rankUpdate(g.H, 1.0);

          for (int r1 = 0; r1 < mg; ++r1) {
            const int i1 = g.items[r1];
            const int off_1 = i1 * q_i;
            const double* z1 = ZiL_u.data() + static_cast<std::ptrdiff_t>(g.cells[r1]) * q_i;
            for (int r2 = 0; r2 <= r1; ++r2) {
              const double q_val = g.Q(r1, r2);
              if (q_val == 0.0) continue;
              const int i2 = g.items[r2];
              const int off_2 = i2 * q_i;
              const double* z2 = ZiL_u.data() + static_cast<std::ptrdiff_t>(g.cells[r2]) * q_i;
              for (int b = 0; b < q_i; ++b) {
                const double qb = q_val * z2[b];
                double* col_p = S_22.data() + static_cast<std::ptrdiff_t>(off_2 + b) * dim_2 + off_1;
                for (int a = 0; a < q_i; ++a) {
                  col_p[a] -= qb * z1[a];
                }
              }
            }
          }

          g.T_2b.noalias() = g.H * g.C_bs_grp.transpose();
          for (int r = 0; r < mg; ++r) {
            const int off_i = g.items[r] * q_i;
            S_2b_saved.middleRows(off_i, q_i).noalias() -=
              ZiL_u.row(g.cells[r]).transpose() * g.T_2b.row(r);
          }
        }
      } else {
        // Populate G_mat = C_is * Rs^T directly in one pass
        const bool overwrite_all = (unique_subj_item && N == n_subj * n_item);
        if (scalar_two_way && unique_subj_item && (overwrite_all || use_sparse_items)) {
          const double ls_li = ZsL_u(0, 0) * ZiL_u(0, 0);
          const double* rp = R_ass_slab.data();
          for (int s = 0; s < n_subj; ++s) {
            const double mult_s = ls_li * rp[s];
            double* g_col = G_mat.data() + static_cast<std::ptrdiff_t>(s) * dim_2;
            for (int k = subj_ptr[s]; k < subj_ptr[s + 1]; ++k) {
              g_col[i_ptr[k]] = mult_s * w_ptr[k];
            }
          }
        } else {
          for (int s = 0; s < n_subj; ++s) {
            const int off_s = s * q_s;
            const auto Rs = R_ass_slab.middleCols(off_s, q_s);
            for (int k = subj_ptr[s]; k < subj_ptr[s + 1]; ++k) {
              const double wn = w_ptr[k];
              const int c = c_ptr[k];
              const int off_i = i_ptr[k] * q_i;
              const double* zs = ZsL_u.data() + static_cast<std::ptrdiff_t>(c) * q_s;
              const double* zi = ZiL_u.data() + static_cast<std::ptrdiff_t>(c) * q_i;
              for (int a = 0; a < q_s; ++a) {
                double r_zs = 0.0;
                for (int m = 0; m <= a; ++m) r_zs += Rs(a, m) * zs[m];
                const double scale = wn * r_zs;
                double* g_col = G_mat.data() + static_cast<std::ptrdiff_t>(off_s + a) * dim_2 + off_i;
                if ((overwrite_all || use_sparse_items) && unique_subj_item) {
                  for (int b = 0; b < q_i; ++b) g_col[b] = scale * zi[b];
                } else {
                  for (int b = 0; b < q_i; ++b) g_col[b] += scale * zi[b];
                }
              }
            }
          }
        }

        if (use_sparse_items) {
          for (int s = 0; s < n_subj; ++s) {
            const int off_s = s * q_s;
            const auto& ui = subj_u_items[s];
            const int m_s = static_cast<int>(ui.size());
            for (int r1 = 0; r1 < m_s; ++r1) {
              const int off_1 = ui[r1] * q_i;
              const auto G1 = G_mat.block(off_1, off_s, q_i, q_s);
              for (int r2 = 0; r2 <= r1; ++r2) {
                const int off_2 = ui[r2] * q_i;
                S_22.block(off_1, off_2, q_i, q_i).noalias() -=
                  G1 * G_mat.block(off_2, off_s, q_i, q_s).transpose();
              }
            }
          }
        } else {
          S_22.selfadjointView<Eigen::Lower>().rankUpdate(G_mat, -1.0);
        }
        S_2b_saved.noalias() = C_2b_buf - G_mat * C_bs.transpose();
      }
    } else {
      // General path (including testlets or >32 continuous covariate rows)
      G_mat.setZero();
      C_bs.setZero();
      C_2b_buf.setZero();
      H_bb_buf = P_beta;

      rhs_s = -v_s;
      rhs_2 = -v_2;
      rhs_b.noalias() = X_rm.transpose() * g_vec - P_beta * beta;

      if (use_archetypes) {
        for (int s = 0; s < n_subj; ++s) {
          Ass_tmp.setIdentity();
          const int off_s = s * q_s;
          for (int k = subj_ptr[s]; k < subj_ptr[s + 1]; ++k) {
            const double wn = w_vec[k];
            const double gn = g_vec[k];
            const int c = cell_arch_idx[k];
            const double* zs = ZsL_u.data() + static_cast<std::ptrdiff_t>(c) * q_s;
            const double* zi = ZiL_u.data() + static_cast<std::ptrdiff_t>(c) * q_i;
            const double* xn = X_u.data() + static_cast<std::ptrdiff_t>(c) * P;
            const int off_i = item_idx[k] * q_i;

            for (int a = 0; a < q_s; ++a) {
              const double w_zsa = wn * zs[a];
              rhs_s[off_s + a] += gn * zs[a];
              for (int b = 0; b < q_s; ++b) Ass_tmp(b, a) += w_zsa * zs[b];
              for (int b = 0; b < q_i; ++b) G_mat(off_i + b, off_s + a) += w_zsa * zi[b];
              for (int p = 0; p < P; ++p) C_bs(p, off_s + a) += w_zsa * xn[p];
            }
            for (int a = 0; a < q_i; ++a) {
              const double w_zia = wn * zi[a];
              rhs_2[off_i + a] += gn * zi[a];
              for (int b = 0; b < q_i; ++b) S_22(off_i + b, off_i + a) += w_zia * zi[b];
              for (int p = 0; p < P; ++p) C_2b_buf(off_i + a, p) += w_zia * xn[p];
            }
            for (int p1 = 0; p1 < P; ++p1) {
              const double w_xp1 = wn * xn[p1];
              for (int p2 = 0; p2 < P; ++p2) H_bb_buf(p2, p1) += w_xp1 * xn[p2];
            }
            if (has_testlet) {
              const double* zt = ZtL_u.data() + static_cast<std::ptrdiff_t>(c) * q_t;
              const int off_t = dim_i + testlet_idx[k] * q_t;
              for (int a = 0; a < q_t; ++a) {
                const double w_zta = wn * zt[a];
                rhs_2[off_t + a] += gn * zt[a];
                for (int b = 0; b < q_t; ++b) S_22(off_t + b, off_t + a) += w_zta * zt[b];
                for (int b = 0; b < q_i; ++b) {
                  const double v_it = w_zta * zi[b];
                  S_22(off_i + b, off_t + a) += v_it;
                  S_22(off_t + a, off_i + b) += v_it;
                }
                for (int b = 0; b < q_s; ++b) G_mat(off_t + a, off_s + b) += w_zta * zs[b];
                for (int p = 0; p < P; ++p) C_2b_buf(off_t + a, p) += w_zta * xn[p];
              }
            }
          }
          auto Rs = R_ass_slab.middleCols(off_s, q_s);
          if (!factor_invert_small_spd(Ass_tmp, Rs, llt_s_ws, logdet_Ass)) return false;
          for (int j = q_s - 1; j >= 0; --j) {
            G_mat.col(off_s + j) *= Rs(j, j);
            C_bs.col(off_s + j) *= Rs(j, j);
            for (int m = 0; m < j; ++m) {
              G_mat.col(off_s + j).noalias() += Rs(j, m) * G_mat.col(off_s + m);
              C_bs.col(off_s + j).noalias() += Rs(j, m) * C_bs.col(off_s + m);
            }
          }
        }
      } else {
        for (int s = 0; s < n_subj; ++s) {
          Ass_tmp.setIdentity();
          const int off_s = s * q_s;
          for (int k = subj_ptr[s]; k < subj_ptr[s + 1]; ++k) {
            const double wn = w_vec[k];
            const double gn = g_vec[k];
            const double* zs = ZsL.data() + static_cast<std::ptrdiff_t>(k) * q_s;
            const double* zi = ZiL.data() + static_cast<std::ptrdiff_t>(k) * q_i;
            const double* xn = X_rm.data() + static_cast<std::ptrdiff_t>(k) * P;
            const int off_i = item_idx[k] * q_i;

            for (int a = 0; a < q_s; ++a) {
              const double w_zsa = wn * zs[a];
              rhs_s[off_s + a] += gn * zs[a];
              for (int b = 0; b < q_s; ++b) Ass_tmp(b, a) += w_zsa * zs[b];
              for (int b = 0; b < q_i; ++b) G_mat(off_i + b, off_s + a) += w_zsa * zi[b];
              for (int p = 0; p < P; ++p) C_bs(p, off_s + a) += w_zsa * xn[p];
            }
            for (int a = 0; a < q_i; ++a) {
              const double w_zia = wn * zi[a];
              rhs_2[off_i + a] += gn * zi[a];
              for (int b = 0; b < q_i; ++b) S_22(off_i + b, off_i + a) += w_zia * zi[b];
              for (int p = 0; p < P; ++p) C_2b_buf(off_i + a, p) += w_zia * xn[p];
            }
            for (int p1 = 0; p1 < P; ++p1) {
              const double w_xp1 = wn * xn[p1];
              for (int p2 = 0; p2 < P; ++p2) H_bb_buf(p2, p1) += w_xp1 * xn[p2];
            }
            if (has_testlet) {
              const double* zt = ZtL.data() + static_cast<std::ptrdiff_t>(k) * q_t;
              const int off_t = dim_i + testlet_idx[k] * q_t;
              for (int a = 0; a < q_t; ++a) {
                const double w_zta = wn * zt[a];
                rhs_2[off_t + a] += gn * zt[a];
                for (int b = 0; b < q_t; ++b) S_22(off_t + b, off_t + a) += w_zta * zt[b];
                for (int b = 0; b < q_i; ++b) {
                  const double v_it = w_zta * zi[b];
                  S_22(off_i + b, off_t + a) += v_it;
                  S_22(off_t + a, off_i + b) += v_it;
                }
                for (int b = 0; b < q_s; ++b) G_mat(off_t + a, off_s + b) += w_zta * zs[b];
                for (int p = 0; p < P; ++p) C_2b_buf(off_t + a, p) += w_zta * xn[p];
              }
            }
          }
          auto Rs = R_ass_slab.middleCols(off_s, q_s);
          if (!factor_invert_small_spd(Ass_tmp, Rs, llt_s_ws, logdet_Ass)) return false;
          for (int j = q_s - 1; j >= 0; --j) {
            G_mat.col(off_s + j) *= Rs(j, j);
            C_bs.col(off_s + j) *= Rs(j, j);
            for (int m = 0; m < j; ++m) {
              G_mat.col(off_s + j).noalias() += Rs(j, m) * G_mat.col(off_s + m);
              C_bs.col(off_s + j).noalias() += Rs(j, m) * C_bs.col(off_s + m);
            }
          }
        }
      }

      S_22.selfadjointView<Eigen::Lower>().rankUpdate(G_mat, -1.0);
      S_2b_saved.noalias() = C_2b_buf - G_mat * C_bs.transpose();
    }

    llt_22.compute(S_22);
    if (llt_22.info() != Eigen::Success) return false;
    L_22_saved = llt_22.matrixL();
    double logdet_S22 = 0.0;
    for (int j = 0; j < dim_2; ++j) logdet_S22 += 2.0 * std::log(L_22_saved(j, j));

    S22_inv_S2b_saved.noalias() = llt_22.solve(S_2b_saved);
    H_bb_schur.noalias() = H_bb_buf - C_bs * C_bs.transpose() - S_2b_saved.transpose() * S22_inv_S2b_saved;

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
      for (int k = 0; k < N; ++k) {
        if (w_vec(k) < 1e-6) w_vec(k) = 1e-6;
      }
      if (assemble_and_factor_core(beta, v_s, v_2, P_beta, rhs_s, rhs_2, rhs_b, llt_22, logdet)) {
        return logdet;
      }
    }
    Rcpp::stop("non_pd_hessian: random-effects curvature block");
  }

  Eigen::VectorXd grad_beta_from_factored(
      const Eigen::VectorXd& rhs_s,
      const Eigen::VectorXd& rhs_2,
      const Eigen::VectorXd& rhs_b,
      const Eigen::LLT<Eigen::MatrixXd>& llt_22,
      Eigen::VectorXd* out_S22_inv_r2 = nullptr
  ) {
    compute_c_lap_fast();
    apply_R_ass(cs_vec, cs_tilde_buf);
    apply_R_ass(rhs_s, rs_tilde_buf);
    sub_G_mul_vec(c2_vec, cs_tilde_buf, r_2_tilde_buf);
    S22_inv_c2_buf.noalias() = llt_22.solve(r_2_tilde_buf);
    sub_G_mul_vec(rhs_2, rs_tilde_buf, r_2_tilde_buf);
    S22_inv_r2_buf.noalias() = llt_22.solve(r_2_tilde_buf);
    if (out_S22_inv_r2 != nullptr) *out_S22_inv_r2 = S22_inv_r2_buf;

    Eigen::VectorXd X_clap(P);
    if (use_archetypes) {
      w_arch.setZero();
      for (int k = 0; k < N; ++k) w_arch[cell_arch_idx[k]] += c_lap[k];
      X_clap.noalias() = X_u.transpose() * w_arch;
    } else {
      X_clap.noalias() = X_rm.transpose() * c_lap;
    }

    return (rhs_b + X_clap) -
      C_bs * (rs_tilde_buf + cs_tilde_buf) - S_2b_saved.transpose() * (S22_inv_r2_buf + S22_inv_c2_buf);
  }

  Eigen::VectorXd compute_laplace_grad_beta(
      const Eigen::VectorXd& beta_eval,
      const Eigen::MatrixXd& L_s_eval,
      const Eigen::MatrixXd& L_i_eval,
      const Eigen::MatrixXd& L_t_eval,
      const Eigen::MatrixXd& P_beta,
      double chance,
      Eigen::VectorXd v_s_cur,
      Eigen::VectorXd v_2_cur,
      int update_zl_mask = 7
  ) {
    if (update_zl_mask != 0) {
      update_ZL(L_s_eval, L_i_eval, L_t_eval, update_zl_mask);
    }

    bool obs_factored = false;
    for (int it = 0; it < 8; ++it) {
      double ll_cur = eval_trials(beta_eval, v_s_cur, v_2_cur, chance, true);
      double phi_cur = ll_cur - 0.5 * (v_s_cur.squaredNorm() + v_2_cur.squaredNorm());
      double logdet_tmp = 0.0;
      bool ok_obs = assemble_and_factor_core(
        beta_eval, v_s_cur, v_2_cur, P_beta, rhs_s_buf, rhs_2_buf, rhs_b_buf, llt_22_grad_ws, logdet_tmp
      );
      double st_max = 0.0;
      double slope = 0.0;
      if (ok_obs) {
        apply_R_ass(rhs_s_buf, rs_tilde_buf);
        sub_G_mul_vec(rhs_2_buf, rs_tilde_buf, r_2_tilde_buf);
        step_2_buf.noalias() = llt_22_grad_ws.solve(r_2_tilde_buf);
        sub_Gt_mul_vec(rs_tilde_buf, step_2_buf, nullptr, step_s_buf);
        apply_Rt_ass_in_place(step_s_buf);
        st_max = std::max(step_s_buf.lpNorm<Eigen::Infinity>(), step_2_buf.lpNorm<Eigen::Infinity>());
        if (st_max < conditional_mode_tol) {
          obs_factored = true;
          break;
        }
        slope = rhs_s_buf.dot(step_s_buf) + rhs_2_buf.dot(step_2_buf);
      }
      if (!ok_obs || !(slope > 0.0)) {
        ll_cur = eval_trials(beta_eval, v_s_cur, v_2_cur, chance, false);
        assemble_and_factor(beta_eval, v_s_cur, v_2_cur, P_beta, rhs_s_buf, rhs_2_buf, rhs_b_buf, llt_22_grad_ws, true);
        apply_R_ass(rhs_s_buf, rs_tilde_buf);
        sub_G_mul_vec(rhs_2_buf, rs_tilde_buf, r_2_tilde_buf);
        step_2_buf.noalias() = llt_22_grad_ws.solve(r_2_tilde_buf);
        sub_Gt_mul_vec(rs_tilde_buf, step_2_buf, nullptr, step_s_buf);
        apply_Rt_ass_in_place(step_s_buf);
        st_max = std::max(step_s_buf.lpNorm<Eigen::Infinity>(), step_2_buf.lpNorm<Eigen::Infinity>());
        if (st_max < conditional_mode_tol) break;
        slope = rhs_s_buf.dot(step_s_buf) + rhs_2_buf.dot(step_2_buf);
      }
      double sc = std::min(1.0, 4.0 / st_max);
      compute_deta_vec(nullptr, step_s_buf, step_2_buf);
      const double vs_norm2 = v_s_cur.squaredNorm() + v_2_cur.squaredNorm();
      const double vs_dot = v_s_cur.dot(step_s_buf) + v_2_cur.dot(step_2_buf);
      const double step_norm2 = step_s_buf.squaredNorm() + step_2_buf.squaredNorm();
      bool accepted_bt = false;
      double best_sc = 0.0;
      double best_phi = phi_cur;
      for (int bt = 0; bt < 8; ++bt) {
        double pen_try = vs_norm2 + 2.0 * sc * vs_dot + (sc * sc) * step_norm2;
        double phi_try = eval_data_ll_step(sc, chance) - 0.5 * pen_try;
        if (std::isfinite(phi_try) && phi_try > best_phi) {
          best_phi = phi_try;
          best_sc = sc;
        }
        if (std::isfinite(phi_try) && phi_try >= phi_cur + 1e-4 * sc * std::max(slope, 0.0)) {
          v_s_cur.noalias() += sc * step_s_buf;
          v_2_cur.noalias() += sc * step_2_buf;
          accepted_bt = true;
          break;
        }
        sc *= 0.5;
      }
      if (!accepted_bt) {
        if (best_sc > 0.0) {
          v_s_cur.noalias() += best_sc * step_s_buf;
          v_2_cur.noalias() += best_sc * step_2_buf;
        } else {
          break;
        }
      }
    }

    if (!obs_factored) {
      eval_trials(beta_eval, v_s_cur, v_2_cur, chance, true);
      assemble_and_factor(beta_eval, v_s_cur, v_2_cur, P_beta, rhs_s_buf, rhs_2_buf, rhs_b_buf, llt_22_grad_ws, true);
    }
    return grad_beta_from_factored(rhs_s_buf, rhs_2_buf, rhs_b_buf, llt_22_grad_ws);
  }

  FitResult fit(
      const Eigen::VectorXi& y_orig,
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
      bool return_ranef = false,
      double default_sd0_s = 0.65,
      double default_sd0_i = 0.55
  ) {
    FitResult res;

    for (int k = 0; k < N; ++k) {
      y_csr[k] = y_orig[perm[k]];
    }

    Eigen::VectorXd init_s = Eigen::VectorXd::Constant(q_s, 0.25);
    init_s(0) = default_sd0_s;
    Eigen::MatrixXd Sigma_s = (init_s.array().square()).matrix().asDiagonal();
    if (init_Sigma_s != nullptr && init_Sigma_s->rows() == q_s && init_Sigma_s->cols() == q_s && init_Sigma_s->allFinite()) {
      Sigma_s = fix_cov ? (*init_Sigma_s) : (0.92 * (*init_Sigma_s) + 0.08 * Sigma_s);
    }

    Eigen::VectorXd init_i = Eigen::VectorXd::Constant(q_i, 0.25);
    init_i(0) = default_sd0_i;
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

    double sum_y = 0.0;
    for (int k = 0; k < N; ++k) sum_y += static_cast<double>(y_csr[k]);
    double mean_y = std::min(std::max(sum_y / std::max(N, 1), chance + 0.05), 0.95);
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

    bool fin_obs_ready = false;
    const bool profile_covariance = (chance == 0.0 && corr_mode == 0);
    double covariance_score_norm = profile_covariance && !fix_cov ? INFINITY : 0.0;
    auto covariance_stationary = [&]() {
      // A score below 0.5% of the local information scale is numerically
      // negligible relative to sampling uncertainty; verify against direct
      // profile optimization rather than tightening EM iterations indefinitely.
      return !profile_covariance || covariance_score_norm < 5e-3;
    };
    // Reuse the Laplace beta/mode refinement inside covariance updates so
    // covariance and beta solve the same approximate marginal objective.
    auto polish_laplace = [&](const Eigen::MatrixXd& L_s_fin,
                              const Eigen::MatrixXd& L_i_fin,
                              const Eigen::MatrixXd& L_t_fin, int beta_steps = 3) {
      fin_obs_ready = false;
      update_ZL(L_s_fin, L_i_fin, L_t_fin, 7);

      Eigen::VectorXd S22_inv_r2(dim_2);
      for (int polish = 0; polish < beta_steps; ++polish) {
        eval_trials(beta, v_s, v_2, chance, true);
        assemble_and_factor(beta, v_s, v_2, P_beta, rhs_s_buf, rhs_2_buf, rhs_b_buf, llt_22_fin_ws, true);
        Eigen::VectorXd grad_lap_b = grad_beta_from_factored(rhs_s_buf, rhs_2_buf, rhs_b_buf, llt_22_fin_ws, &S22_inv_r2);

        ldlt_b_ws.compute(H_bb_schur);
        d_beta_buf.noalias() = ldlt_b_ws.solve(grad_lap_b);
        double lap_norm = d_beta_buf.lpNorm<Eigen::Infinity>();
        if (lap_norm < 1.5) {
          d_v2_buf.noalias() = S22_inv_r2 - S22_inv_S2b_saved * d_beta_buf;
          sub_Gt_mul_vec(rs_tilde_buf, d_v2_buf, &d_beta_buf, d_vs_buf);
          apply_Rt_ass_in_place(d_vs_buf);
          beta += d_beta_buf;
          v_2 += d_v2_buf;
          v_s += d_vs_buf;
        }
        if (lap_norm < 1e-5) break;
      }

      inner_converged = false;
      for (int it = 0; it < max_pirls_iter; ++it) {
        double ll_cur = eval_trials(beta, v_s, v_2, chance, true);
        double phi_cur = ll_cur - 0.5 * (v_s.squaredNorm() + v_2.squaredNorm());
        double logdet_obs = 0.0;
        bool ok_obs = assemble_and_factor_core(
          beta, v_s, v_2, P_beta, rhs_s_buf, rhs_2_buf, rhs_b_buf, llt_22_fin_ws, logdet_obs
        );
        double slope = 0.0;
        if (ok_obs) {
          apply_R_ass(rhs_s_buf, rs_tilde_buf);
          sub_G_mul_vec(rhs_2_buf, rs_tilde_buf, r_2_tilde_buf);
          step_2_buf.noalias() = llt_22_fin_ws.solve(r_2_tilde_buf);
          sub_Gt_mul_vec(rs_tilde_buf, step_2_buf, nullptr, step_s_buf);
          apply_Rt_ass_in_place(step_s_buf);
          final_inner_step = std::max(step_s_buf.lpNorm<Eigen::Infinity>(), step_2_buf.lpNorm<Eigen::Infinity>());
          if (final_inner_step < conditional_mode_tol) {
            log_lik_laplace = ll_cur - 0.5 * (v_s.squaredNorm() + v_2.squaredNorm() + logdet_obs);
            inner_converged = true;
            fin_obs_ready = true;
            break;
          }
          slope = rhs_s_buf.dot(step_s_buf) + rhs_2_buf.dot(step_2_buf);
        }
        if (!ok_obs || !(slope > 0.0)) {
          ll_cur = eval_trials(beta, v_s, v_2, chance, false);
          assemble_and_factor(beta, v_s, v_2, P_beta, rhs_s_buf, rhs_2_buf, rhs_b_buf, llt_22_fin_ws, true);
          apply_R_ass(rhs_s_buf, rs_tilde_buf);
          sub_G_mul_vec(rhs_2_buf, rs_tilde_buf, r_2_tilde_buf);
          step_2_buf.noalias() = llt_22_fin_ws.solve(r_2_tilde_buf);
          sub_Gt_mul_vec(rs_tilde_buf, step_2_buf, nullptr, step_s_buf);
          apply_Rt_ass_in_place(step_s_buf);
          final_inner_step = std::max(step_s_buf.lpNorm<Eigen::Infinity>(), step_2_buf.lpNorm<Eigen::Infinity>());
          if (final_inner_step < conditional_mode_tol) {
            const double ll_obs = eval_trials(beta, v_s, v_2, chance, true);
            logdet_obs = assemble_and_factor(beta, v_s, v_2, P_beta, rhs_s_buf, rhs_2_buf, rhs_b_buf, llt_22_fin_ws, true);
            log_lik_laplace = ll_obs - 0.5 * (v_s.squaredNorm() + v_2.squaredNorm() + logdet_obs);
            inner_converged = true;
            fin_obs_ready = true;
            break;
          }
          slope = rhs_s_buf.dot(step_s_buf) + rhs_2_buf.dot(step_2_buf);
        }
        double sc = std::min(1.0, 4.0 / final_inner_step);
        compute_deta_vec(nullptr, step_s_buf, step_2_buf);
        const double vs_norm2 = v_s.squaredNorm() + v_2.squaredNorm();
        const double vs_dot = v_s.dot(step_s_buf) + v_2.dot(step_2_buf);
        const double step_norm2 = step_s_buf.squaredNorm() + step_2_buf.squaredNorm();
        bool accepted_bt = false;
        double best_sc = 0.0;
        double best_phi = phi_cur;
        for (int bt = 0; bt < 8; ++bt) {
          double pen_try = vs_norm2 + 2.0 * sc * vs_dot + (sc * sc) * step_norm2;
          double phi_try = eval_data_ll_step(sc, chance) - 0.5 * pen_try;
          if (std::isfinite(phi_try) && phi_try > best_phi) {
            best_phi = phi_try;
            best_sc = sc;
          }
          if (std::isfinite(phi_try) &&
              (phi_try >= phi_cur + 1e-4 * sc * std::max(slope, 0.0) ||
               ((sc * final_inner_step) < 1e-7 && phi_try >= phi_cur))) {
            v_s.noalias() += sc * step_s_buf;
            v_2.noalias() += sc * step_2_buf;
            accepted_bt = true;
            break;
          }
          sc *= 0.5;
        }
        if (!accepted_bt) {
          if (best_sc > 0.0) {
            v_s.noalias() += best_sc * step_s_buf;
            v_2.noalias() += best_sc * step_2_buf;
          } else {
            break;
          }
        }
      }
    };

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

      update_ZL(L_s, L_i, L_t, 7);

      inner_converged = false;
      bool obs_factored = false;

      for (int it = 0; it < max_pirls_iter; ++it) {
        double ll_data = 0.0;
        double logdet_H = 0.0;
        double max_step = 0.0;
        double slope = 0.0;
        bool used_obs = false;

        if (chance > 0.0 && (outer_used >= 2 || it >= 2)) {
          ll_data = eval_trials(beta, v_s, v_2, chance, true);
          bool ok_obs = assemble_and_factor_core(
            beta, v_s, v_2, P_beta, rhs_s_buf, rhs_2_buf, rhs_b_buf, llt_22_ws, logdet_H
          );
          if (ok_obs) {
            ldlt_b_ws.compute(H_bb_schur);
            ok_obs = (ldlt_b_ws.info() == Eigen::Success) && (ldlt_b_ws.vectorD().minCoeff() > 1e-8);
          }
          if (ok_obs) {
            apply_R_ass(rhs_s_buf, rs_tilde_buf);
            sub_G_mul_vec(rhs_2_buf, rs_tilde_buf, r_2_tilde_buf);
            S22_inv_r2_buf.noalias() = llt_22_ws.solve(r_2_tilde_buf);
            r_b_tilde_buf.noalias() = rhs_b_buf - C_bs * rs_tilde_buf - S_2b_saved.transpose() * S22_inv_r2_buf;
            d_beta_buf.noalias() = ldlt_b_ws.solve(r_b_tilde_buf);
          if (profile_covariance && update_sigma) d_beta_buf.setZero();
            d_v2_buf.noalias() = S22_inv_r2_buf - S22_inv_S2b_saved * d_beta_buf;
            sub_Gt_mul_vec(rs_tilde_buf, d_v2_buf, &d_beta_buf, d_vs_buf);
            apply_Rt_ass_in_place(d_vs_buf);

            max_step = std::max({
              d_beta_buf.lpNorm<Eigen::Infinity>(),
              d_v2_buf.lpNorm<Eigen::Infinity>(),
              d_vs_buf.lpNorm<Eigen::Infinity>()
            });
            slope = rhs_b_buf.dot(d_beta_buf) + rhs_s_buf.dot(d_vs_buf) + rhs_2_buf.dot(d_v2_buf);
            if (std::isfinite(max_step) && (max_step < 1.5e-5 || slope > 0.0)) {
              used_obs = true;
            }
          }
        }

        if (!used_obs) {
          ll_data = eval_trials(beta, v_s, v_2, chance, false);
          logdet_H = assemble_and_factor(beta, v_s, v_2, P_beta, rhs_s_buf, rhs_2_buf, rhs_b_buf, llt_22_ws, false);

          apply_R_ass(rhs_s_buf, rs_tilde_buf);
          sub_G_mul_vec(rhs_2_buf, rs_tilde_buf, r_2_tilde_buf);
          S22_inv_r2_buf.noalias() = llt_22_ws.solve(r_2_tilde_buf);
          r_b_tilde_buf.noalias() = rhs_b_buf - C_bs * rs_tilde_buf - S_2b_saved.transpose() * S22_inv_r2_buf;

          ldlt_b_ws.compute(H_bb_schur);
          d_beta_buf.noalias() = ldlt_b_ws.solve(r_b_tilde_buf);
          if (profile_covariance && update_sigma) d_beta_buf.setZero();
          d_v2_buf.noalias() = S22_inv_r2_buf - S22_inv_S2b_saved * d_beta_buf;
          sub_Gt_mul_vec(rs_tilde_buf, d_v2_buf, &d_beta_buf, d_vs_buf);
          apply_Rt_ass_in_place(d_vs_buf);

          max_step = std::max({
            d_beta_buf.lpNorm<Eigen::Infinity>(),
            d_v2_buf.lpNorm<Eigen::Infinity>(),
            d_vs_buf.lpNorm<Eigen::Infinity>()
          });
          slope = rhs_b_buf.dot(d_beta_buf) + rhs_s_buf.dot(d_vs_buf) + rhs_2_buf.dot(d_v2_buf);
        }

        double phi_cur = ll_data - 0.5 * (v_s.squaredNorm() + v_2.squaredNorm() + beta.dot(P_beta * beta));
        final_inner_step = max_step;

        if (!std::isfinite(max_step)) {
          inner_converged = false;
          break;
        }
        if (max_step < 1.5e-5 && (it >= 1 || outer_used >= 3)) {
          log_lik_laplace = ll_data - 0.5 * (v_s.squaredNorm() + v_2.squaredNorm() + logdet_H);
          inner_converged = true;
          obs_factored = used_obs;
          break;
        }

        double step_scale = (max_step > 4.0) ? (4.0 / max_step) : 1.0;
        double accepted_ll = ll_data;

        compute_deta_vec(&d_beta_buf, d_vs_buf, d_v2_buf);
        const double v_norm2 = v_s.squaredNorm() + v_2.squaredNorm();
        const double v_dot = v_s.dot(d_vs_buf) + v_2.dot(d_v2_buf);
        const double dv_norm2 = d_vs_buf.squaredNorm() + d_v2_buf.squaredNorm();

        bool accepted_bt = false;
        double best_scale = 0.0;
        double best_phi = phi_cur;
        double best_ll = ll_data;
        for (int bt = 0; bt < 10; ++bt) {
          b_try_buf.noalias() = beta + step_scale * d_beta_buf;
          double ll_try = eval_data_ll_step(step_scale, chance);
          double v_pen_try = v_norm2 + 2.0 * step_scale * v_dot + (step_scale * step_scale) * dv_norm2;
          double phi_try = ll_try - 0.5 * (v_pen_try + b_try_buf.dot(P_beta * b_try_buf));
          if (std::isfinite(phi_try) && phi_try > best_phi) {
            best_phi = phi_try;
            best_scale = step_scale;
            best_ll = ll_try;
          }
          if (std::isfinite(phi_try) &&
              (phi_try >= phi_cur + 1e-4 * step_scale * std::max(slope, 0.0) ||
               ((step_scale * max_step) < 1e-7 && phi_try >= phi_cur))) {
            beta = b_try_buf;
            v_2.noalias() += step_scale * d_v2_buf;
            v_s.noalias() += step_scale * d_vs_buf;
            accepted_ll = ll_try;
            accepted_bt = true;
            break;
          }
          step_scale *= 0.5;
        }
        if (!accepted_bt) {
          if (best_scale > 0.0) {
            step_scale = best_scale;
            beta.noalias() += step_scale * d_beta_buf;
            v_2.noalias() += step_scale * d_v2_buf;
            v_s.noalias() += step_scale * d_vs_buf;
            accepted_ll = best_ll;
          } else {
            step_scale = 0.0;
            accepted_ll = ll_data;
            log_lik_laplace = accepted_ll - 0.5 * (v_s.squaredNorm() + v_2.squaredNorm() + logdet_H);
            break;
          }
        }

        log_lik_laplace = accepted_ll - 0.5 * (v_s.squaredNorm() + v_2.squaredNorm() + logdet_H);

        if ((step_scale * max_step) < 1.5e-5 && max_step < 5.0e-4 && (it >= 1 || outer_used >= 3)) {
          inner_converged = true;
          break;
        }
      }

      if (profile_covariance && update_sigma && inner_converged) {
        polish_laplace(L_s, L_i, L_t, 1);
        obs_factored = fin_obs_ready;
      }

      if (chance > 0.0 && inner_converged && !obs_factored) {
        eval_trials(beta, v_s, v_2, chance, true);
        assemble_and_factor(beta, v_s, v_2, P_beta, rhs_s_buf, rhs_2_buf, rhs_b_buf, llt_22_ws, true);
      }

      Eigen::VectorXd mean_shift_s = Eigen::VectorXd::Zero(dim_s);
      Eigen::VectorXd mean_shift_2 = Eigen::VectorXd::Zero(dim_2);
      if (profile_covariance && update_sigma && inner_converged && fin_obs_ready) {
        // Include movement of the conditional mode in the Laplace determinant
        // derivative: u * mean_shift' + mean_shift * u' in E(uu').
        compute_c_lap_fast();
        apply_R_ass(cs_vec, cs_tilde_buf);
        sub_G_mul_vec(c2_vec, cs_tilde_buf, r_2_tilde_buf);
        mean_shift_2.noalias() = llt_22_fin_ws.solve(r_2_tilde_buf);
        sub_Gt_mul_vec(cs_tilde_buf, mean_shift_2, nullptr, mean_shift_s);
        apply_Rt_ass_in_place(mean_shift_s);
      } else {
        invert_L22_lower();
        compute_W_mat_from_L22_inv();
      }

      rel_sum_s.setZero();
      sum_uu_s_buf.setZero();
      for (int s = 0; s < n_subj; ++s) {
        const int off_s = s * q_s;
        us_buf.noalias() = L_s * v_s.segment(off_s, q_s);
        const auto Ws = W_mat.middleCols(off_s, q_s);
        const auto Rs = R_ass_slab.middleCols(off_s, q_s);
        V_tmp_s.noalias() = Ws.transpose() * Ws;
        for (int j = 0; j < q_s; ++j) V_tmp_s(j, j) += 1.0;
        auto V_vs = V_ss_slab.middleCols(off_s, q_s);
        V_vs.noalias() = Rs.transpose() * V_tmp_s * Rs;
        for (int j = 0; j < q_s; ++j) {
          double rel = std::min(std::max(1.0 - V_vs(j, j), 0.0), 1.0);
          rel_sum_s(j) += 2.0 * rel * rel;
        }
        sum_uu_s_buf.noalias() += us_buf * us_buf.transpose() + L_s * V_vs * L_s.transpose();
        if (profile_covariance && update_sigma) {
          Eigen::VectorXd shift = L_s * mean_shift_s.segment(off_s, q_s);
          sum_uu_s_buf.noalias() += us_buf * shift.transpose() + shift * us_buf.transpose();
        }
      }

      rel_sum_i.setZero();
      sum_uu_i_buf.setZero();
      for (int i = 0; i < n_item; ++i) {
        const int off_i = i * q_i;
        ui_buf.noalias() = L_i * v_2.segment(off_i, q_i);
        const auto Li_sub = L22_inv.block(off_i, off_i, dim_2 - off_i, q_i);
        V_tmp_i.noalias() = Li_sub.transpose() * Li_sub;
        for (int j = 0; j < q_i; ++j) {
          double rel = std::min(std::max(1.0 - V_tmp_i(j, j), 0.0), 1.0);
          rel_sum_i(j) += 2.0 * rel * rel;
        }
        sum_uu_i_buf.noalias() += ui_buf * ui_buf.transpose() + L_i * V_tmp_i * L_i.transpose();
        if (profile_covariance && update_sigma) {
          Eigen::VectorXd shift = L_i * mean_shift_2.segment(off_i, q_i);
          sum_uu_i_buf.noalias() += ui_buf * shift.transpose() + shift * ui_buf.transpose();
        }
      }

      rel_sum_t.setZero();
      if (has_testlet) {
        sum_uu_t_buf.setZero();
        for (int t = 0; t < n_testlet; ++t) {
          const int off_t = dim_i + t * q_t;
          ut_buf.noalias() = L_t * v_2.segment(off_t, q_t);
          const auto Lt_sub = L22_inv.block(off_t, off_t, dim_2 - off_t, q_t);
          V_tmp_t.noalias() = Lt_sub.transpose() * Lt_sub;
          for (int j = 0; j < q_t; ++j) {
            double rel = std::min(std::max(1.0 - V_tmp_t(j, j), 0.0), 1.0);
            rel_sum_t(j) += 2.0 * rel * rel;
          }
          sum_uu_t_buf.noalias() += ut_buf * ut_buf.transpose() + L_t * V_tmp_t * L_t.transpose();
          if (profile_covariance && update_sigma) {
            Eigen::VectorXd shift = L_t * mean_shift_2.segment(off_t, q_t);
            sum_uu_t_buf.noalias() += ut_buf * shift.transpose() + shift * ut_buf.transpose();
          }
        }
      }

      if (profile_covariance && update_sigma) {
        covariance_score_norm = 0.0;
        auto check_score = [&](const Eigen::MatrixXd& moment, const Eigen::MatrixXd& Sigma,
                               const Eigen::VectorXd& rel_sum, int n) {
          for (int j = 0; j < Sigma.rows(); ++j) {
            double tau2 = std::max(Sigma(j, j), 1e-10);
            double score = moment(j, j) / tau2 - n;
            double info_prior = 0.0;
            if (regularized) {
              double r = tau2 / (3.0 * prior_tau_sd * prior_tau_sd);
              score += 1.0 - 4.0 * r / (1.0 + r);
              info_prior = 8.0 * r / ((1.0 + r) * (1.0 + r));
            }
            // Unpenalized zero-variance solutions use a one-sided boundary score.
            if (!regularized && tau2 <= 1e-10 && score < 0.0) score = 0.0;
            covariance_score_norm = std::max(covariance_score_norm,
              std::abs(score) / std::sqrt(std::max(rel_sum(j) + info_prior, 0.75)));
          }
        };
        check_score(sum_uu_s_buf, Sig_s_in, rel_sum_s, n_subj);
        check_score(sum_uu_i_buf, Sig_i_in, rel_sum_i, n_item);
        if (has_testlet) check_score(sum_uu_t_buf, Sig_t_in, rel_sum_t, n_testlet);
      }

      if (!update_sigma) {
        Sig_s_out = Sig_s_in;
        Sig_i_out = Sig_i_in;
        if (has_testlet) Sig_t_out = Sig_t_in;
        return;
      }

      Sig_s_out = update_cov_matrix(sum_uu_s_buf, Sig_s_in, n_subj, corr_mode, regularized, lkj_eta, prior_tau_sd);
      Sig_i_out = update_cov_matrix(sum_uu_i_buf, Sig_i_in, n_item, corr_mode, regularized, lkj_eta, prior_tau_sd);
      if (has_testlet) {
        Sig_t_out = update_cov_matrix(sum_uu_t_buf, Sig_t_in, n_testlet, corr_mode, regularized, lkj_eta, prior_tau_sd);
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

        if (final_outer_diff < tol && outer_used >= 3 && inner_converged && covariance_stationary()) {
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

        if (final_outer_diff < tol && outer_used >= 3 && inner_converged && covariance_stationary()) {
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
          if (profile_covariance) d3 = std::max(d3, (beta - beta_bak).lpNorm<Eigen::Infinity>());

            if (d3 <= 4.0 * std::max(d1, d2) || d3 < 5e-3) {
              final_outer_diff = d3;
              Sigma_s = S_s3; Sigma_i = S_i3; if (has_testlet) Sigma_t = S_t3;
              if (d3 < 0.5 * tol && inner_converged && covariance_stationary()) {
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

    if (inner_converged && outer_converged) polish_laplace(L_s_fin, L_i_fin, L_t_fin);

    Eigen::MatrixXd H_bb_saved = H_bb_schur;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_cond(0.5 * (H_bb_saved + H_bb_saved.transpose()));
    bool pos_def_b = (es_cond.eigenvalues().minCoeff() > 1e-8) && std::isfinite(log_lik_laplace);
    Eigen::VectorXd ev_cond = es_cond.eigenvalues().cwiseMax(1e-8);
    Eigen::MatrixXd vcov_cond = es_cond.eigenvectors() * ev_cond.cwiseInverse().asDiagonal() *
      es_cond.eigenvectors().transpose();

    double s0_sq = prior_tau_sd * prior_tau_sd;
    auto make_var_log_tau = [&](const Eigen::MatrixXd& Sig, const Eigen::VectorXd& rel_sum) {
      int q = Sig.rows();
      Eigen::VectorXd out(q);
      for (int j = 0; j < q; ++j) {
        double tau2 = std::max(Sig(j, j), 1e-10);
        double info_prior = 0.0;
        if (regularized && s0_sq > 0.0) {
          double r = tau2 / (3.0 * s0_sq);
          info_prior = (8.0 * r) / ((1.0 + r) * (1.0 + r));
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

      // Precompute exact Newton mode tangents wrt beta from the already-factored base Schur system
      Eigen::MatrixXd dv2_dbeta(dim_2, P);
      Eigen::MatrixXd dvs_dbeta(dim_s, P);
      Eigen::VectorXd g_base;
      if (fin_obs_ready) {
        dv2_dbeta = -S22_inv_S2b_saved;
        for (int p = 0; p < P; ++p) {
          d_beta_buf.setZero();
          d_beta_buf[p] = 1.0;
          sub_Gt_mul_vec(Eigen::VectorXd::Zero(dim_s), dv2_dbeta.col(p), &d_beta_buf, step_s_buf);
          apply_Rt_ass_in_place(step_s_buf);
          dvs_dbeta.col(p) = step_s_buf;
        }
        g_base = grad_beta_from_factored(rhs_s_buf, rhs_2_buf, rhs_b_buf, llt_22_fin_ws);
      } else {
        dv2_dbeta.setZero();
        dvs_dbeta.setZero();
        g_base = compute_laplace_grad_beta(
          beta, L_s_fin, L_i_fin, L_t_fin, P_beta, chance, v_s_opt, v_2_opt, 7
        );
      }

      Eigen::MatrixXd H_lap_bb(P, P);
      for (int p = 0; p < P; ++p) {
        Eigen::VectorXd b_hi = beta;
        b_hi(p) += eps_b;
        Eigen::VectorXd vs_init = v_s_opt + eps_b * dvs_dbeta.col(p);
        Eigen::VectorXd v2_init = v_2_opt + eps_b * dv2_dbeta.col(p);
        Eigen::VectorXd g_hi = compute_laplace_grad_beta(
          b_hi, L_s_fin, L_i_fin, L_t_fin, P_beta, chance, vs_init, v2_init, 0
        );
        H_lap_bb.col(p) = -(g_hi - g_base) / eps_b;
      }
      H_lap_bb = 0.5 * (H_lap_bb + H_lap_bb.transpose());

      double min_floor = std::max(0.25 * ev_cond.minCoeff(), 1e-6);
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_lap(H_lap_bb);
      Eigen::VectorXd ev_lap = es_lap.eigenvalues().cwiseMax(min_floor);
      Eigen::MatrixXd vcov_lap_fixed = es_lap.eigenvectors() * ev_lap.cwiseInverse().asDiagonal() *
        es_lap.eigenvectors().transpose();

      if (fix_cov) {
        vcov_beta = vcov_lap_fixed;
      } else {
        // Restore base Schur system once to warm-start covariance perturbations via 1-step Schur tangent
        if (fin_obs_ready) {
          eval_trials(beta, v_s_opt, v_2_opt, chance, true);
          double ld_tmp = 0.0;
          fin_obs_ready = assemble_and_factor_core(
            beta, v_s_opt, v_2_opt, P_beta, rhs_s_buf, rhs_2_buf, rhs_b_buf, llt_22_fin_ws, ld_tmp
          );
        }
        const Eigen::MatrixXd R_ass_base = R_ass_slab;
        const std::vector<ScheduleGroup> groups_base = groups;
        const Eigen::MatrixXd G_mat_base = G_mat;

        auto schur_solve_base = [&](const Eigen::VectorXd& r_s_in, const Eigen::VectorXd& r_2_in,
                                    Eigen::VectorXd& dv_s_out, Eigen::VectorXd& dv_2_out) {
          R_ass_slab = R_ass_base;
          if (use_group_rank1) groups = groups_base; else G_mat = G_mat_base;
          apply_R_ass(r_s_in, rs_tilde_buf);
          sub_G_mul_vec(r_2_in, rs_tilde_buf, r_2_tilde_buf);
          dv_2_out.noalias() = llt_22_fin_ws.solve(r_2_tilde_buf);
          sub_Gt_mul_vec(rs_tilde_buf, dv_2_out, nullptr, dv_s_out);
          apply_Rt_ass_in_place(dv_s_out);
        };

        const double eps_t = 1e-3;
        const double mult_hi = std::exp(eps_t);
        Eigen::VectorXd var_theta_all(n_theta);
        int col_k = 0;

        for (int j = 0; j < q_s; ++j, ++col_k) {
          var_theta_all(col_k) = res.var_log_tau_s(j);
          Eigen::MatrixXd Ls_hi = safe_chol(scale_cov_single_sd(Sigma_s, j, mult_hi));
          Eigen::VectorXd vs_init = v_s_opt;
          Eigen::VectorXd v2_init = v_2_opt;
          if (fin_obs_ready) {
            Eigen::MatrixXd map_inv = Ls_hi.triangularView<Eigen::Lower>().solve(L_s_fin);
            Eigen::MatrixXd map_tr = Ls_hi.transpose() * L_s_fin.triangularView<Eigen::Lower>().transpose().solve(Eigen::MatrixXd::Identity(q_s, q_s));
            Eigen::VectorXd r_s_pert(dim_s);
            for (int s = 0; s < n_subj; ++s) {
              const auto vs_seg = v_s_opt.segment(s * q_s, q_s);
              vs_init.segment(s * q_s, q_s).noalias() = map_inv * vs_seg;
              r_s_pert.segment(s * q_s, q_s).noalias() = map_tr * vs_seg - vs_init.segment(s * q_s, q_s);
            }
            Eigen::VectorXd dvs_t(dim_s), dv2_t(dim_2);
            schur_solve_base(r_s_pert, Eigen::VectorXd::Zero(dim_2), dvs_t, dv2_t);
            vs_init += dvs_t;
            v2_init += dv2_t;
          }
          Eigen::VectorXd g_hi = compute_laplace_grad_beta(
            beta, Ls_hi, L_i_fin, L_t_fin, P_beta, chance, vs_init, v2_init, 1
          );
          res.J_beta_theta.col(col_k) = vcov_lap_fixed * ((g_hi - g_base) / eps_t);
        }
        if (q_s > 0 && (q_i > 0 || q_t > 0)) update_ZL(L_s_fin, L_i_fin, L_t_fin, 1);

        for (int j = 0; j < q_i; ++j, ++col_k) {
          var_theta_all(col_k) = res.var_log_tau_i(j);
          Eigen::MatrixXd Li_hi = safe_chol(scale_cov_single_sd(Sigma_i, j, mult_hi));
          Eigen::VectorXd vs_init = v_s_opt;
          Eigen::VectorXd v2_init = v_2_opt;
          if (fin_obs_ready) {
            Eigen::MatrixXd map_inv = Li_hi.triangularView<Eigen::Lower>().solve(L_i_fin);
            Eigen::MatrixXd map_tr = Li_hi.transpose() * L_i_fin.triangularView<Eigen::Lower>().transpose().solve(Eigen::MatrixXd::Identity(q_i, q_i));
            Eigen::VectorXd r_2_pert = Eigen::VectorXd::Zero(dim_2);
            for (int i = 0; i < n_item; ++i) {
              const auto vi_seg = v_2_opt.segment(i * q_i, q_i);
              v2_init.segment(i * q_i, q_i).noalias() = map_inv * vi_seg;
              r_2_pert.segment(i * q_i, q_i).noalias() = map_tr * vi_seg - v2_init.segment(i * q_i, q_i);
            }
            Eigen::VectorXd dvs_t(dim_s), dv2_t(dim_2);
            schur_solve_base(Eigen::VectorXd::Zero(dim_s), r_2_pert, dvs_t, dv2_t);
            vs_init += dvs_t;
            v2_init += dv2_t;
          }
          Eigen::VectorXd g_hi = compute_laplace_grad_beta(
            beta, L_s_fin, Li_hi, L_t_fin, P_beta, chance, vs_init, v2_init, 2
          );
          res.J_beta_theta.col(col_k) = vcov_lap_fixed * ((g_hi - g_base) / eps_t);
        }
        if (q_i > 0 && q_t > 0) update_ZL(L_s_fin, L_i_fin, L_t_fin, 2);

        for (int j = 0; j < q_t; ++j, ++col_k) {
          var_theta_all(col_k) = res.var_log_tau_t(j);
          Eigen::MatrixXd Lt_hi = safe_chol(scale_cov_single_sd(Sigma_t, j, mult_hi));
          Eigen::VectorXd vs_init = v_s_opt;
          Eigen::VectorXd v2_init = v_2_opt;
          if (fin_obs_ready) {
            Eigen::MatrixXd map_inv = Lt_hi.triangularView<Eigen::Lower>().solve(L_t_fin);
            Eigen::MatrixXd map_tr = Lt_hi.transpose() * L_t_fin.triangularView<Eigen::Lower>().transpose().solve(Eigen::MatrixXd::Identity(q_t, q_t));
            Eigen::VectorXd r_2_pert = Eigen::VectorXd::Zero(dim_2);
            for (int t = 0; t < n_testlet; ++t) {
              const int off_t = dim_i + t * q_t;
              const auto vt_seg = v_2_opt.segment(off_t, q_t);
              v2_init.segment(off_t, q_t).noalias() = map_inv * vt_seg;
              r_2_pert.segment(off_t, q_t).noalias() = map_tr * vt_seg - v2_init.segment(off_t, q_t);
            }
            Eigen::VectorXd dvs_t(dim_s), dv2_t(dim_2);
            schur_solve_base(Eigen::VectorXd::Zero(dim_s), r_2_pert, dvs_t, dv2_t);
            vs_init += dvs_t;
            v2_init += dv2_t;
          }
          Eigen::VectorXd g_hi = compute_laplace_grad_beta(
            beta, L_s_fin, L_i_fin, Lt_hi, P_beta, chance, vs_init, v2_init, 4
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
    res.covariance_score_norm = covariance_score_norm;
    res.log_lik_laplace = log_lik_laplace;
    res.final_outer_diff = final_outer_diff;
    res.final_inner_step = final_inner_step;
    res.outer_used = outer_used;
    res.singular = singular;
    return res;
  }
};

// Wrapper holding GLMMWorkspace + automatic (subj <-> item) orientation so S_22 always has min(dim_s, dim_i)
struct OrientedGLMMWorkspace {
  bool swapped = false;
  int orig_n_subj = 0, orig_n_item = 0, orig_n_testlet = 0;
  int orig_q_s = 0, orig_q_i = 0, orig_q_t = 0;
  Eigen::VectorXi orig_subj_idx, orig_item_idx, orig_testlet_idx;
  RowMatrixXd orig_X, orig_Zs, orig_Zi, orig_Zt;
  GLMMWorkspace ws;

  OrientedGLMMWorkspace(
      const Eigen::MatrixXd& X,
      const Eigen::VectorXi& subj_idx,
      const Eigen::MatrixXd& Z_subj,
      int n_subj,
      const Eigen::VectorXi& item_idx,
      const Eigen::MatrixXd& Z_item,
      int n_item,
      const Eigen::VectorXi& testlet_idx,
      const Eigen::MatrixXd& Z_testlet,
      int n_testlet
  ) : swapped((n_testlet == 0 || Z_testlet.cols() == 0) &&
              (static_cast<long long>(n_item) * Z_item.cols() >
               static_cast<long long>(n_subj) * Z_subj.cols())),
      orig_n_subj(n_subj), orig_n_item(n_item), orig_n_testlet(n_testlet),
      orig_q_s(Z_subj.cols()), orig_q_i(Z_item.cols()),
      orig_q_t((n_testlet > 0 && Z_testlet.cols() > 0) ? Z_testlet.cols() : 0),
      orig_subj_idx(subj_idx), orig_item_idx(item_idx), orig_testlet_idx(testlet_idx),
      orig_X(X), orig_Zs(Z_subj), orig_Zi(Z_item),
      orig_Zt((n_testlet > 0 && Z_testlet.cols() > 0) ? RowMatrixXd(Z_testlet) : RowMatrixXd(0, 0)),
      ws(X,
         swapped ? item_idx : subj_idx,
         swapped ? Z_item : Z_subj,
         swapped ? n_item : n_subj,
         swapped ? subj_idx : item_idx,
         swapped ? Z_subj : Z_item,
         swapped ? n_subj : n_item,
         testlet_idx, Z_testlet, n_testlet)
  {}

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
      const Eigen::MatrixXd* init_Sigma_subj,
      const Eigen::MatrixXd* init_Sigma_item,
      const Eigen::MatrixXd* init_Sigma_testlet,
      bool fix_cov = false,
      bool compute_joint_vcov = true,
      bool return_ranef = false
  ) {
    if (!swapped) {
      return ws.fit(
        y, chance, corr_mode, regularized, lkj_eta, prior_beta_sd, prior_tau_sd,
        max_em_iter, max_pirls_iter, tol, init_beta, init_Sigma_subj, init_Sigma_item, init_Sigma_testlet,
        fix_cov, compute_joint_vcov, return_ranef, 0.65, 0.55
      );
    }
    FitResult r = ws.fit(
      y, chance, corr_mode, regularized, lkj_eta, prior_beta_sd, prior_tau_sd,
      max_em_iter, max_pirls_iter, tol, init_beta, init_Sigma_item, init_Sigma_subj, init_Sigma_testlet,
      fix_cov, compute_joint_vcov, return_ranef, 0.55, 0.65
    );
    std::swap(r.Sigma_s, r.Sigma_i);
    std::swap(r.var_log_tau_s, r.var_log_tau_i);
    if (return_ranef) std::swap(r.u_subj, r.u_item);
    if (r.J_beta_theta.cols() == orig_q_s + orig_q_i) {
      Eigen::MatrixXd J_orig(r.J_beta_theta.rows(), orig_q_s + orig_q_i);
      J_orig.leftCols(orig_q_s) = r.J_beta_theta.rightCols(orig_q_s);
      J_orig.rightCols(orig_q_i) = r.J_beta_theta.leftCols(orig_q_i);
      r.J_beta_theta = std::move(J_orig);
    }
    return r;
  }
};

Rcpp::List wrap_fit_result(const FitResult& res, bool has_testlet, bool regularized, bool fix_cov, bool return_ranef) {
  return Rcpp::List::create(
    Rcpp::Named("beta") = res.beta,
    Rcpp::Named("se") = res.se_beta,
    Rcpp::Named("vcov") = res.vcov_beta,
    Rcpp::Named("vcov_cond") = res.vcov_cond,
    Rcpp::Named("Sigma_subj") = res.Sigma_s,
    Rcpp::Named("Sigma_item") = res.Sigma_i,
    Rcpp::Named("Sigma_testlet") = has_testlet ? Rcpp::wrap(res.Sigma_t) : R_NilValue,
    Rcpp::Named("var_log_tau_subj") = res.var_log_tau_s,
    Rcpp::Named("var_log_tau_item") = res.var_log_tau_i,
    Rcpp::Named("var_log_tau_testlet") = has_testlet ? Rcpp::wrap(res.var_log_tau_t) : R_NilValue,
    Rcpp::Named("J_beta_theta") = res.J_beta_theta,
    Rcpp::Named("u_subj") = return_ranef ? Rcpp::wrap(res.u_subj) : R_NilValue,
    Rcpp::Named("u_item") = return_ranef ? Rcpp::wrap(res.u_item) : R_NilValue,
    Rcpp::Named("u_testlet") = (return_ranef && has_testlet) ? Rcpp::wrap(res.u_testlet) : R_NilValue,
    Rcpp::Named("logLik") = res.log_lik_laplace,
    Rcpp::Named("covariance_score_norm") = res.covariance_score_norm,
    Rcpp::Named("singular") = (regularized || fix_cov)
      ? Rcpp::LogicalVector::create(NA_LOGICAL)
      : Rcpp::LogicalVector::create(res.singular),
    Rcpp::Named("converged") = res.strictly_converged,
    Rcpp::Named("status") = res.status,
    Rcpp::Named("iterations") = res.outer_used,
    Rcpp::Named("outer_diff") = res.final_outer_diff,
    Rcpp::Named("inner_step") = res.final_inner_step
  );
}

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

int constant_response_cell_impl(const Eigen::VectorXi& y, const Eigen::VectorXi& cell_id, int n_cells) {
  if (y.size() != cell_id.size() || n_cells < 1) Rcpp::stop("Invalid response cell dimensions");
  std::vector<unsigned char> mask(n_cells, 0);
  const int* yp = y.data();
  const int* cp = cell_id.data();
  const int n = y.size();
  for (int i = 0; i < n; ++i) {
    const int cell = cp[i] - 1;
    if (cell < 0 || cell >= n_cells) Rcpp::stop("Invalid response cell index");
    mask[cell] |= (yp[i] == 1) ? 2u : 1u;
  }
  for (int cell = 0; cell < n_cells; ++cell) {
    if (mask[cell] == 1u || mask[cell] == 2u) return cell + 1;
  }
  return 0;
}

Rcpp::NumericVector run_prob_bootstrap_on_ws(
    OrientedGLMMWorkspace& ows,
    int n_boot,
    const Eigen::VectorXi& cell_id,
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
    double chance,
    int corr_mode,
    bool regularized,
    double lkj_eta,
    double prior_beta_sd,
    double prior_tau_sd,
    int max_em_iter,
    int max_pirls_iter,
    double tol,
    int eff_idx_1based,
    bool marginal
) {
  Rcpp::RNGScope rng_scope;

  const int N = ows.orig_X.rows();
  const int n_subj = ows.orig_n_subj;
  const int n_item = ows.orig_n_item;
  const int n_testlet = ows.orig_n_testlet;
  const int q_s = ows.orig_q_s;
  const int q_i = ows.orig_q_i;
  const bool has_t = (ows.orig_q_t > 0);
  const int q_t = ows.orig_q_t;
  const int eff_idx0 = std::max(0, eff_idx_1based - 1);

  Eigen::MatrixXd St_hat = (has_t && Sigma_testlet_hat.isNotNull()) ?
    Rcpp::as<Eigen::MatrixXd>(Sigma_testlet_hat) : Eigen::MatrixXd(0, 0);

  RowMatrixXd Gs = ows.orig_Zs * chol_lower_pd(Sigma_subj_hat);
  RowMatrixXd Gi = ows.orig_Zi * chol_lower_pd(Sigma_item_hat);
  RowMatrixXd Gt = has_t ? RowMatrixXd(ows.orig_Zt * chol_lower_pd(St_hat)) : RowMatrixXd(0, 0);
  Eigen::VectorXd fixed_eta = ows.orig_X * beta_hat;

  Eigen::MatrixXd u_s(n_subj, q_s);
  Eigen::MatrixXd u_i(n_item, q_i);
  Eigen::MatrixXd u_t(has_t ? n_testlet : 0, has_t ? q_t : 0);
  Eigen::VectorXi y_star(N);

  const Eigen::MatrixXd* p_St = has_t ? &St_hat : nullptr;
  Rcpp::NumericVector boot_est(n_boot, NA_REAL);

  for (int b = 0; b < n_boot; ++b) {
    if ((b & 3) == 0) Rcpp::checkUserInterrupt();
    // Preserve exact column-major normal draw order matching R's matrix(rnorm(n*q), nrow=n)
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
        Gs.row(n).dot(u_s.row(ows.orig_subj_idx(n))) +
        Gi.row(n).dot(u_i.row(ows.orig_item_idx(n)));
      if (has_t) {
        eta += Gt.row(n).dot(u_t.row(ows.orig_testlet_idx(n)));
      }
      double prob = chance + (1.0 - chance) * inv_logit(eta);
      prob = std::min(std::max(prob, 0.0), 1.0);
      y_star(n) = static_cast<int>(R::rbinom(1.0, prob));
    }

    if (!regularized && constant_response_cell_impl(y_star, cell_id, X_cells.rows())) continue;
    try {
      FitResult f = ows.fit(
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

} // anonymous namespace

// Pre-Allocate Reusable Schedule Workspace (XPtr) Across Replications
// [[Rcpp::export(rng = false)]]
SEXP cpp_create_glmm_workspace(
    const Eigen::MatrixXd& X,
    const Eigen::VectorXi& subj_idx,
    const Eigen::MatrixXd& Z_subj,
    int n_subj,
    const Eigen::VectorXi& item_idx,
    const Eigen::MatrixXd& Z_item,
    int n_item,
    const Eigen::VectorXi& testlet_idx,
    const Eigen::MatrixXd& Z_testlet,
    int n_testlet
) {
  Rcpp::XPtr<OrientedGLMMWorkspace> ptr(
    new OrientedGLMMWorkspace(X, subj_idx, Z_subj, n_subj, item_idx, Z_item, n_item, testlet_idx, Z_testlet, n_testlet),
    true
  );
  return ptr;
}

// [[Rcpp::export(rng = false)]]
bool cpp_glmm_workspace_valid(SEXP ws_ptr) {
  return TYPEOF(ws_ptr) == EXTPTRSXP && R_ExternalPtrAddr(ws_ptr) != nullptr;
}

// [[Rcpp::export]]
Rcpp::List cpp_fit_crossed_glmm_ws(
    SEXP ws_ptr,
    const Eigen::VectorXi& y,
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
  Rcpp::XPtr<OrientedGLMMWorkspace> ows(ws_ptr);
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

  FitResult res = ows->fit(
    y, chance, corr_mode, regularized, lkj_eta, prior_beta_sd, prior_tau_sd,
    max_em_iter, max_pirls_iter, tol, p_bw, p_Ss, p_Si, p_St,
    fix_cov, compute_joint_vcov, return_ranef
  );
  return wrap_fit_result(res, ows->orig_q_t > 0, regularized, fix_cov, return_ranef);
}

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
  OrientedGLMMWorkspace ows(X, subj_idx, Z_subj, n_subj, item_idx, Z_item, n_item, testlet_idx, Z_testlet, n_testlet);

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

  FitResult res = ows.fit(
    y, chance, corr_mode, regularized, lkj_eta, prior_beta_sd, prior_tau_sd,
    max_em_iter, max_pirls_iter, tol, p_bw, p_Ss, p_Si, p_St,
    fix_cov, compute_joint_vcov, return_ranef
  );
  return wrap_fit_result(res, ows.orig_q_t > 0, regularized, fix_cov, return_ranef);
}

// [[Rcpp::export(rng = false)]]
int cpp_constant_response_cell(const Eigen::VectorXi& y, const Eigen::VectorXi& cell_id,
                               int n_cells) {
  return constant_response_cell_impl(y, cell_id, n_cells);
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_prob_bootstrap_batch_ws(
    SEXP ws_ptr,
    int n_boot,
    const Eigen::VectorXi& cell_id,
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
  Rcpp::XPtr<OrientedGLMMWorkspace> ows(ws_ptr);
  return run_prob_bootstrap_on_ws(
    *ows, n_boot, cell_id, beta_hat, Sigma_subj_hat, Sigma_item_hat, Sigma_testlet_hat,
    X_cells, Zs_cells, Zi_cells, Zt_cells, contrast_weights, gh_nodes, gh_weights,
    chance, corr_mode, regularized, lkj_eta, prior_beta_sd, prior_tau_sd,
    max_em_iter, max_pirls_iter, tol, eff_idx_1based, marginal
  );
}

// Pre-Allocated Zero-Copy C++ Parametric Bootstrap Batch for Marginal Probability Contrasts
//
// [[Rcpp::export]]
Rcpp::NumericVector cpp_prob_bootstrap_batch(
    int n_boot,
    const Eigen::VectorXi& cell_id,
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
  OrientedGLMMWorkspace ows(X, subj_idx, Z_subj, n_subj, item_idx, Z_item, n_item, testlet_idx, Z_testlet, n_testlet);
  return run_prob_bootstrap_on_ws(
    ows, n_boot, cell_id, beta_hat, Sigma_subj_hat, Sigma_item_hat, Sigma_testlet_hat,
    X_cells, Zs_cells, Zi_cells, Zt_cells, contrast_weights, gh_nodes, gh_weights,
    chance, corr_mode, regularized, lkj_eta, prior_beta_sd, prior_tau_sd,
    max_em_iter, max_pirls_iter, tol, eff_idx_1based, marginal
  );
}

// Single-pass zero-allocation unit response tally (1-based IDs)
// [[Rcpp::export(rng = false)]]
Rcpp::List cpp_unit_response_counts(
    const Rcpp::IntegerVector& id,
    const Rcpp::LogicalVector& success,
    int n
) {
  const int len = id.size();
  Rcpp::IntegerVector observed(n, 0);
  Rcpp::IntegerVector correct(n, 0);
  const int* id_p = id.begin();
  const int* suc_p = success.begin();
  int* obs_p = observed.begin();
  int* cor_p = correct.begin();

  for (int k = 0; k < len; ++k) {
    const int idx = id_p[k] - 1;
    if (idx >= 0 && idx < n) {
      obs_p[idx]++;
      if (suc_p[k] == TRUE) cor_p[idx]++;
    }
  }

  Rcpp::LogicalVector all_correct(n);
  Rcpp::LogicalVector all_incorrect(n);
  Rcpp::LogicalVector unobserved(n);
  int* ac_p = all_correct.begin();
  int* ai_p = all_incorrect.begin();
  int* un_p = unobserved.begin();

  for (int j = 0; j < n; ++j) {
    const int o = obs_p[j];
    const int c = cor_p[j];
    un_p[j] = (o == 0);
    ac_p[j] = (o > 0 && c == o);
    ai_p[j] = (o > 0 && c == 0);
  }

  return Rcpp::List::create(
    Rcpp::Named("observed") = observed,
    Rcpp::Named("correct") = correct,
    Rcpp::Named("all_correct") = all_correct,
    Rcpp::Named("all_incorrect") = all_incorrect,
    Rcpp::Named("unobserved") = unobserved
  );
}

// Fast C++ Gauss-Hermite MML-EM for Fixed-Bank Rasch Calibration (Sufficient-Score Accelerated)
// [[Rcpp::export(rng = false)]]
Rcpp::List cpp_fit_rasch_fixed_mml(
    const Eigen::VectorXi& subj_id,
    const Eigen::VectorXi& item_id,
    const Eigen::VectorXi& y,
    int n_subj,
    int n_item,
    const Rcpp::LogicalVector& estimable,
    const Eigen::VectorXd& gh_nodes,
    const Eigen::VectorXd& gh_weights,
    int max_iter = 200,
    double tol = 1e-6
) {
  std::vector<int> est_j;
  est_j.reserve(n_item);
  std::vector<int> item_to_col(n_item, -1);
  for (int j = 0; j < n_item; ++j) {
    if (estimable[j] == TRUE) {
      item_to_col[j] = static_cast<int>(est_j.size());
      est_j.push_back(j);
    }
  }
  const int J = static_cast<int>(est_j.size());
  if (J == 0) {
    return Rcpp::List::create(Rcpp::Named("error") = "No item has both response categories");
  }

  // Count observed estimable items per subject first
  std::vector<int> subj_obs_cnt(n_subj, 0);
  std::vector<int> subj_raw_score(n_subj, 0);
  const int n_tr = y.size();
  for (int k = 0; k < n_tr; ++k) {
    const int c = item_to_col[item_id[k] - 1];
    if (c >= 0) {
      const int s = subj_id[k] - 1;
      subj_obs_cnt[s]++;
      subj_raw_score[s] += y[k];
    }
  }

  std::vector<int> obs_s;
  obs_s.reserve(n_subj);
  std::vector<int> subj_to_row(n_subj, -1);
  bool all_obs = true;
  for (int s = 0; s < n_subj; ++s) {
    if (subj_obs_cnt[s] > 0) {
      subj_to_row[s] = static_cast<int>(obs_s.size());
      obs_s.push_back(s);
      if (subj_obs_cnt[s] != J) all_obs = false;
    }
  }
  const int S_obs = static_cast<int>(obs_s.size());
  if (S_obs < 2) {
    return Rcpp::List::create(Rcpp::Named("error") = "Fewer than two observed participants");
  }

  Eigen::VectorXd n_sub = Eigen::VectorXd::Zero(J);
  Eigen::VectorXd s_sub = Eigen::VectorXd::Zero(J);
  Eigen::MatrixXd M_sub, Y_sub;
  if (!all_obs) {
    M_sub = Eigen::MatrixXd::Zero(S_obs, J);
    Y_sub = Eigen::MatrixXd::Zero(S_obs, J);
  }
  for (int k = 0; k < n_tr; ++k) {
    const int c = item_to_col[item_id[k] - 1];
    if (c < 0) continue;
    const int r = subj_to_row[subj_id[k] - 1];
    if (r < 0) continue;
    const double yv = static_cast<double>(y[k]);
    if (!all_obs) {
      if (M_sub(r, c) == 0.0) {
        M_sub(r, c) = 1.0;
        n_sub[c] += 1.0;
        s_sub[c] += yv;
      }
      Y_sub(r, c) = yv;
    } else {
      s_sub[c] += yv;
    }
  }
  if (all_obs) {
    n_sub.setConstant(static_cast<double>(S_obs));
  }

  // double check all_obs in case of duplicate rows in trials
  if (all_obs) {
    for (int j = 0; j < J; ++j) {
      if (s_sub[j] < 0.0 || s_sub[j] > static_cast<double>(S_obs)) all_obs = false;
    }
  }

  const int Q = gh_nodes.size();
  Eigen::VectorXd log_w = gh_weights.array().log();
  Eigen::VectorXd b(J);
  for (int j = 0; j < J; ++j) {
    double p0 = std::min(std::max(s_sub[j] / n_sub[j], 0.02), 0.98);
    b[j] = -std::log(p0 / (1.0 - p0));
  }
  double sigma = 0.80;
  bool converged = false;

  Eigen::VectorXd theta_q(Q), col_sp(Q), col_R(Q);
  Eigen::MatrixXd eta_jq(J, Q), sp_jq(J, Q), N_jq(J, Q);

  // Sufficient-score grouping when all observed participants responded to every estimable item
  std::vector<int> active_scores;
  std::vector<double> score_counts;
  Eigen::MatrixXd R_sq;
  Eigen::MatrixXd LL_iq, R_iq;

  if (all_obs) {
    std::vector<int> cnt(J + 1, 0);
    for (int s : obs_s) cnt[subj_raw_score[s]]++;
    for (int sc = 0; sc <= J; ++sc) {
      if (cnt[sc] > 0) {
        active_scores.push_back(sc);
        score_counts.push_back(static_cast<double>(cnt[sc]));
      }
    }
    R_sq.resize(static_cast<int>(active_scores.size()), Q);
  } else {
    LL_iq.resize(S_obs, Q);
    R_iq.resize(S_obs, Q);
  }

  auto e_step = [&]() {
    double log_likelihood = all_obs ? -s_sub.dot(b) : 0.0;
    theta_q = sigma * gh_nodes;
    for (int q = 0; q < Q; ++q) {
      double sp_sum = 0.0;
      const double tq = theta_q[q];
      for (int j = 0; j < J; ++j) {
        const double eta = -b[j] + tq;
        eta_jq(j, q) = eta;
        const double sp = softplus(eta);
        sp_jq(j, q) = sp;
        sp_sum += sp;
      }
      col_sp[q] = sp_sum;
    }

    if (all_obs) {
      col_R.setZero();
      const int n_sc = static_cast<int>(active_scores.size());
      for (int idx_s = 0; idx_s < n_sc; ++idx_s) {
        const double sc = static_cast<double>(active_scores[idx_s]);
        double mx = -1e300;
        for (int q = 0; q < Q; ++q) {
          const double lp = sc * theta_q[q] - col_sp[q] + log_w[q];
          R_sq(idx_s, q) = lp;
          if (lp > mx) mx = lp;
        }
        double rsum = 0.0;
        for (int q = 0; q < Q; ++q) {
          const double w = std::exp(R_sq(idx_s, q) - mx);
          R_sq(idx_s, q) = w;
          rsum += w;
        }
        const double inv_rsum = 1.0 / rsum;
        const double mult = score_counts[idx_s];
        log_likelihood += mult * (mx + std::log(rsum));
        for (int q = 0; q < Q; ++q) {
          R_sq(idx_s, q) *= inv_rsum;
          col_R[q] += mult * R_sq(idx_s, q);
        }
      }
    } else {
      LL_iq.noalias() = Y_sub * eta_jq - M_sub * sp_jq;
      col_R.setZero();
      for (int i = 0; i < S_obs; ++i) {
        double mx = -1e300;
        for (int q = 0; q < Q; ++q) {
          const double lp = LL_iq(i, q) + log_w[q];
          R_iq(i, q) = lp;
          if (lp > mx) mx = lp;
        }
        double rsum = 0.0;
        for (int q = 0; q < Q; ++q) {
          const double w = std::exp(R_iq(i, q) - mx);
          R_iq(i, q) = w;
          rsum += w;
        }
        const double inv_rsum = 1.0 / rsum;
        log_likelihood += mx + std::log(rsum);
        for (int q = 0; q < Q; ++q) {
          R_iq(i, q) *= inv_rsum;
          col_R[q] += R_iq(i, q);
        }
      }
    }
    return log_likelihood;
  };

  auto em_step = [&]() {
    e_step();

    if (!all_obs) {
      N_jq.noalias() = M_sub.transpose() * R_iq;
    }

    for (int nw = 0; nw < 4; ++nw) {
      double max_st = 0.0;
      for (int j = 0; j < J; ++j) {
        const double bj = b[j];
        double grad = -s_sub[j];
        double hess = 0.0;
        if (all_obs) {
          for (int q = 0; q < Q; ++q) {
            const double p = inv_logit(-bj + theta_q[q]);
            const double np = col_R[q] * p;
            grad += np;
            hess += np * (1.0 - p);
          }
        } else {
          for (int q = 0; q < Q; ++q) {
            const double p = inv_logit(-bj + theta_q[q]);
            const double np = N_jq(j, q) * p;
            grad += np;
            hess += np * (1.0 - p);
          }
        }
        hess = std::max(hess, 1e-8);
        const double step = std::min(std::max(grad / hess, -2.0), 2.0);
        b[j] += step;
        max_st = std::max(max_st, std::abs(step));
      }
      if (max_st < 1e-8) break;
    }

    const double inv_S = 1.0 / static_cast<double>(S_obs);
    const double mu_case = col_R.dot(theta_q) * inv_S;
    b.array() -= mu_case;
    double var_new = 0.0;
    for (int q = 0; q < Q; ++q) {
      const double d = theta_q[q] - mu_case;
      var_new += col_R[q] * d * d;
    }
    sigma = std::sqrt(std::max(var_new * inv_S, 1e-10));
  };

  for (int it = 1; it <= max_iter; ++it) {
    const Eigen::VectorXd b_old = b;
    const double sig_old = sigma;
    em_step();
    const double diff_max = std::max((b - b_old).lpNorm<Eigen::Infinity>(), std::abs(sigma - sig_old));
    if (diff_max < tol && it >= 3) {
      converged = true;
      break;
    }
  }

  if (!converged) {
    // At variance zero the item MLEs are available exactly. A nonpositive
    // profile right derivative and no loss of likelihood are both required;
    // a small positive variance alone is not evidence for a boundary solution.
    Eigen::VectorXd b_zero(J), p_zero(J), bern_var(J);
    for (int j = 0; j < J; ++j) {
      p_zero[j] = s_sub[j] / n_sub[j];
      b_zero[j] = -std::log(p_zero[j] / (1.0 - p_zero[j]));
      bern_var[j] = p_zero[j] * (1.0 - p_zero[j]);
    }
    double boundary_score = 0.0;
    for (int i = 0; i < S_obs; ++i) {
      const double mean = all_obs ? p_zero.sum() : M_sub.row(i).dot(p_zero);
      const double variance = all_obs ? bern_var.sum() : M_sub.row(i).dot(bern_var);
      const double residual = subj_raw_score[obs_s[i]] - mean;
      boundary_score += 0.5 * (residual * residual - variance);
    }
    double boundary_ll = -s_sub.dot(b_zero);
    for (int j = 0; j < J; ++j) boundary_ll -= n_sub[j] * softplus(-b_zero[j]);

    // Safeguarded squared extrapolation of the EM map (Varadhan & Roland,
    // 2008). Work in variance coordinates, then check the ordinary EM
    // residual. Backtracking cannot itself declare convergence.
    for (int it = 0; it < max_iter && !converged; ++it) {
      const double current_ll = e_step();
      if (boundary_score <= 0.0 && boundary_ll >= current_ll - 1e-9) {
        b = b_zero;
        sigma = 0.0;
        converged = true;
        break;
      }
      Eigen::VectorXd x0(J + 1), x1(J + 1), x2(J + 1);
      x0.head(J) = b;
      x0[J] = sigma * sigma;
      const double sig_old = sigma;
      em_step();
      x1.head(J) = b;
      x1[J] = sigma * sigma;
      // Near zero an EM SD step can be tiny far from the likelihood optimum.
      // The posterior second moment gives the variance score per participant.
      e_step();
      const double variance = sigma * sigma;
      const double variance_score = (col_R.dot(theta_q.array().square().matrix()) /
        static_cast<double>(S_obs) - variance) / (2.0 * variance * variance);
      if ((b - x0.head(J)).lpNorm<Eigen::Infinity>() < tol &&
          std::abs(sigma - sig_old) < tol && std::abs(variance_score) < tol) {
        converged = true;
        break;
      }
      em_step();
      x2.head(J) = b;
      x2[J] = sigma * sigma;
      const double baseline_ll = e_step();
      const Eigen::VectorXd r = x1 - x0;
      const Eigen::VectorXd v = x2 - 2.0 * x1 + x0;
      double alpha = v.squaredNorm() > 0.0 ?
        -std::sqrt(r.squaredNorm() / v.squaredNorm()) : -1.0;
      alpha = std::isfinite(alpha) ? std::min(-1.0, alpha) : -1.0;
      bool accepted = false;
      for (int bt = 0; bt < 12; ++bt) {
        const Eigen::VectorXd candidate = x0 - 2.0 * alpha * r + alpha * alpha * v;
        if (candidate.allFinite() && candidate[J] > 0.0) {
          b = candidate.head(J);
          sigma = std::sqrt(candidate[J]);
          em_step();
          const double candidate_ll = e_step();
          if (b.allFinite() && std::isfinite(sigma) && std::isfinite(candidate_ll) &&
              candidate_ll >= baseline_ll - 1e-9) {
            accepted = true;
            break;
          }
        }
        alpha = 0.5 * (alpha - 1.0);
      }
      if (!accepted) {
        b = x2.head(J);
        sigma = std::sqrt(x2[J]);
      }
    }
  }

  if (!converged) {
    return Rcpp::List::create(Rcpp::Named("error") = "MML-EM reached iteration limit");
  }

  e_step();
  const double inv_S = 1.0 / static_cast<double>(S_obs);
  const double mu_case = col_R.dot(theta_q) * inv_S;
  b.array() -= mu_case;

  Rcpp::NumericVector b_full(n_item, NA_REAL);
  for (int c = 0; c < J; ++c) b_full[est_j[c]] = b[c];

  Rcpp::NumericVector ab_full(n_subj, NA_REAL);
  if (all_obs) {
    std::vector<double> score_eap(J + 1, 0.0);
    const int n_sc = static_cast<int>(active_scores.size());
    for (int idx_s = 0; idx_s < n_sc; ++idx_s) {
      score_eap[active_scores[idx_s]] = R_sq.row(idx_s).dot(theta_q) - mu_case;
    }
    for (int s : obs_s) {
      ab_full[s] = score_eap[subj_raw_score[s]];
    }
  } else {
    Eigen::VectorXd eap_obs = R_iq * theta_q;
    for (int r = 0; r < S_obs; ++r) {
      ab_full[obs_s[r]] = eap_obs[r] - mu_case;
    }
  }

  const bool sing = (sigma <= 1e-4);
  const double v_s = sing ? 0.0 : (sigma * sigma);

  return Rcpp::List::create(
    Rcpp::Named("error") = Rcpp::CharacterVector::create(NA_STRING),
    Rcpp::Named("b_est") = b_full,
    Rcpp::Named("ability") = ab_full,
    Rcpp::Named("person_var") = v_s,
    Rcpp::Named("singular") = sing
  );
}
