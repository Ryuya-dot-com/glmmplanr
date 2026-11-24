# Fixed 32 x 16 allocation, two conditions, ordinary logit, random intercepts only.
# Compare the variance of the unadjusted difference in condition means to the
# variance with mutually independent Bernoulli trials of the same marginal risks.
library(glmmplanr)

compare_allocation <- function(allocation, n_quad = 31L, subj_sd = .7, item_sd = .4) {
  design <- design_crossed(list(A = c("control", "treatment")), allocation = allocation)
  p <- c(.60, .68)
  g <- calibrate_dgp(design, p, chance = 0, subj_sd = subj_sd, item_sd = item_sd,
                     slope_sd_ratio = 0, corr_structure = 0, n_quad = n_quad)
  d <- generate_design_schedule(design, 32, 16)$trials
  stopifnot(all(table(d$subj_id, d$item_id) == 1L), all(table(d$cell_id) == 256L))
  levels_per_subj <- rowSums(table(d$subj_id, d$A) > 0L)
  stopifnot(all(levels_per_subj == if (allocation == "between_subj") 1L else 2L))

  # Integrate each trial over its unshared intercept, then take the covariance
  # over the shared intercept. The resulting 2 x 2 matrices are on the 0/1 scale.
  q <- glmmplanr:::gh_quad_1d(n_quad)
  pair_cov <- function(shared_sd, other_sd) {
    offsets <- outer(shared_sd * q$nodes, other_sd * q$nodes, "+")
    means <- vapply(g$mu_cells, function(mu) drop(plogis(mu + offsets) %*% q$weights),
                    numeric(n_quad))
    marginal <- drop(crossprod(q$weights, means))
    stopifnot(max(abs(marginal - p)) < 1e-7)
    centered <- sweep(means, 2L, marginal)
    crossprod(centered, q$weights * centered)
  }
  cs <- pair_cov(subj_sd, item_sd)
  ci <- pair_cov(item_sd, subj_sd)
  cell <- d$cell_id
  V <- cs[cell, cell] * outer(d$subj_id, d$subj_id, "==") +
       ci[cell, cell] * outer(d$item_id, d$item_id, "==")
  diag(V) <- p[cell] * (1 - p[cell])
  w <- c(-1, 1)[cell] / as.numeric(table(cell))[cell]
  var_independent <- sum(w^2 * diag(V))
  var_design <- drop(crossprod(w, V %*% w))
  latent_var <- subj_sd^2 + item_sd^2 + pi^2 / 3
  data.frame(allocation, levels_per_participant = unique(levels_per_subj),
             latent_icc_participant = subj_sd^2 / latent_var, latent_icc_item = item_sd^2 / latent_var,
             response_icc_control_participant = cs[1,1] / (p[1]*(1-p[1])),
             response_icc_control_item = ci[1,1] / (p[1]*(1-p[1])),
             se_independent_pp = 100 * sqrt(var_independent), se_design_pp = 100 * sqrt(var_design),
             deff = var_design / var_independent)
}

layouts <- c("latin_square", "between_subj")
results <- do.call(rbind, lapply(layouts, compare_allocation))
stopifnot(isTRUE(all.equal(results, do.call(rbind, lapply(layouts, compare_allocation, n_quad = 41L)),
                           tolerance = 1e-7)))
for (layout in layouts) stopifnot(abs(compare_allocation(layout, subj_sd = 0, item_sd = 0)$deff - 1) < 1e-9)
print(results, row.names = FALSE)
