three_choice_reference_dgp <- function() {
  des <- design_crossed(list(
    exposure = c("incidental", "intentional"),
    spacing = c("massed", "distributed")
  ))
  calibrate_dgp(
    des,
    cell_probs = c(
      incidental_massed       = 0.50,
      incidental_distributed  = 0.62,
      intentional_massed      = 0.68,
      intentional_distributed = 0.84
    ),
    chance = 1 / 3,
    subj_prob_range = c(0.42, 0.84),
    item_prob_range = c(0.45, 0.86),
    slope_sd_ratio = 0.30,
    corr_structure = 0.15
  )
}

three_choice_stan_data <- function(tr, spec) {
  dgp <- spec$dgp
  des <- dgp$design
  stopifnot(spec$analysis_re == "maximal", !des$has_testlet)
  q <- glmmplanr:::gh_quad_1d(dgp$n_quad)
  list(N = nrow(tr), y = tr$y, P = ncol(des$X_cells), X = des$X_cells[tr$cell_id, ],
    N_subj = nlevels(tr$subj), Q_subj = ncol(des$Z_subj_cells), subj_id = tr$subj_id,
    Z_subj = des$Z_subj_cells[tr$cell_id, ], N_item = nlevels(tr$item),
    Q_item = ncol(des$Z_item_cells), item_id = tr$item_id, Z_item = des$Z_item_cells[tr$cell_id, ],
    N_testlet = 0L, Q_testlet = 0L, testlet_id = integer(), Z_testlet = matrix(numeric(), 0, 0),
    M = des$n_cells, X_cell = des$X_cells, Zs_cell = des$Z_subj_cells, Zi_cell = des$Z_item_cells,
    Zt_cell = matrix(numeric(), des$n_cells, 0), contrast = unname(spec$contrast_weights),
    K = length(q$nodes), nodes = q$nodes, weights = q$weights, marginal = 1L,
    chance = dgp$chance, corr_mode = 1L, lkj_eta = spec$lkj_eta,
    prior_beta_sd = spec$prior_beta_sd, prior_tau_sd = spec$prior_tau_sd)
}
