# Research-only use of the native fixed-effect matrix and population standardization.
library(glmmplanr)

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1L] else "validation"
n_rep <- if (length(args) > 1L) as.integer(args[2L]) else 1024L
workers <- if (length(args) > 2L) as.integer(args[3L]) else
  if (.Platform$OS.type == "windows") 1L else 2L
stopifnot(dir.exists(out), n_rep > 0L, n_rep <= 1024L, workers > 0L)
master_seed <- if (n_rep == 1024L) 85323001L else 85324001L
streams <- glmmplanr:::make_lecuyer_streams(master_seed, n_rep)
grid <- data.frame(n_subj = c(32L, 96L), n_item = c(48L, 16L))
profiles <- data.frame(
  allocation = c("latin_square", "between_subj", "latin_square"),
  vocab_by_method = c(0, 0, .35)
)
cases <- merge(merge(grid, profiles), data.frame(effect_pp = c(0, 8)))
cases$case <- seq_len(nrow(cases))
cases$vocab_logit <- .50
cases$length_logit <- -.35
cases$master_seed <- master_seed
cases$n_rep <- n_rep
z_crit <- stats::qnorm(.975)
q <- glmmplanr:::gh_quad_1d(11L)
nodes <- expand.grid(vocab = seq_along(q$nodes), length = seq_along(q$nodes))
node_weights <- q$weights[nodes$vocab] * q$weights[nodes$length]
stopifnot(abs(sum(node_weights) - 1) < 1e-12)

covariate_cols <- function(vocab, length, x, interaction) {
  cols <- cbind(vocabulary = vocab, word_length = length)
  if (interaction) cols <- cbind(cols, method_vocabulary = x * vocab)
  cols
}

population_spec <- function(spec, interaction) {
  design <- spec$dgp$design
  cells <- rep(1:2, each = nrow(nodes))
  X <- design$X_cells[cells, , drop = FALSE]
  design$X_cells <- cbind(X, covariate_cols(rep(q$nodes[nodes$vocab], 2L),
                                          rep(q$nodes[nodes$length], 2L),
                                          X[, 2L], interaction))
  design$Z_subj_cells <- design$Z_subj_cells[cells, , drop = FALSE]
  design$Z_item_cells <- design$Z_item_cells[cells, , drop = FALSE]
  design$n_cells <- length(cells)
  spec$dgp$design <- design
  spec$contrast_weights <- rep(spec$contrast_weights, each = nrow(nodes)) * rep(node_weights, 2L)
  spec
}

run_case <- function(k) {
  cfg <- cases[k, ]
  ns <- cfg$n_subj
  ni <- cfg$n_item
  design <- design_crossed(list(method = c("restudy", "retrieval")), allocation = cfg$allocation)
  slope_subj <- if (cfg$allocation == "latin_square") .25 else 0
  slope_item <- .25
  subj_residual <- sqrt(.8^2 - cfg$vocab_logit^2)
  item_residual <- sqrt(.6^2 - cfg$length_logit^2)
  aggregate_slope <- sqrt(slope_subj^2 + cfg$vocab_by_method^2)
  aggregate_corr <- if (aggregate_slope > 0) {
    cfg$vocab_logit * cfg$vocab_by_method / (.8 * aggregate_slope)
  } else 0
  probs <- c(restudy = .42, retrieval = .42 + cfg$effect_pp / 100)
  dgp <- calibrate_dgp(
    design, probs, subj_sd = if (slope_subj > 0) c(.8, aggregate_slope) else .8,
    item_sd = c(.6, slope_item), corr_structure = list(subj = aggregate_corr, item = 0)
  )
  spec <- specify_target(dgp, "method", scale = "pp", engine = "bayes_laplace",
                         analysis_re = "maximal", prob_se_method = "joint_delta",
                         control = target_control(on_singular = "flag_only"))
  cache <- glmmplanr:::prepare_cpp_schedule_cache(spec, generate_design_schedule(design, ns, ni))
  x <- cache$X_full[, 2L]
  is_restudy <- x < 0
  adjusted <- population_spec(spec, interaction = FALSE)
  moderator <- population_spec(spec, interaction = TRUE)

  # Check equal truth by the adjusted-model integration and the aggregate normal DGP.
  ps <- ncol(design$Z_subj_cells)
  sig_subj <- diag(c(subj_residual, slope_subj)[seq_len(ps)]^2, nrow = ps)
  sig_item <- diag(c(item_residual, slope_item)^2)
  population <- if (cfg$vocab_by_method != 0) moderator else adjusted
  beta <- c(dgp$beta, cfg$vocab_logit, cfg$length_logit,
            if (cfg$vocab_by_method != 0) cfg$vocab_by_method)
  zero <- matrix(0, length(beta), length(beta))
  pop_check <- glmmplanr:::estimate_prob_effect_delta(
    beta, zero, sig_subj, sig_item, NULL, population$dgp$design, 0,
    population$contrast_weights, n_quad = dgp$n_quad
  )$est
  stopifnot(abs(pop_check * 100 - cfg$effect_pp) < 1e-6, sum(is_restudy) == ns * ni / 2L)

  rows <- lapply(seq_len(n_rep), function(i) {
    draw <- glmmplanr:::with_rng_stream(streams[[i]], {
      subjects <- matrix(stats::rnorm(3L * max(grid$n_subj)), max(grid$n_subj), 3L)
      items <- matrix(stats::rnorm(3L * max(grid$n_item)), max(grid$n_item), 3L)
      vocab <- subjects[cache$sid, 1L]
      length <- items[cache$iid, 1L]
      eta <- cache$mu_trial + (cfg$vocab_logit + cfg$vocab_by_method * x) * vocab +
        cfg$length_logit * length +
        subj_residual * subjects[cache$sid, 2L] + slope_subj * subjects[cache$sid, 3L] * x +
        item_residual * items[cache$iid, 2L] + slope_item * items[cache$iid, 3L] * x
      pairs <- cache$sid + (cache$iid - 1L) * max(grid$n_subj)
      u <- stats::runif(max(grid$n_subj) * max(grid$n_item))[pairs]
      list(y = as.integer(u < stats::plogis(eta)), vocab = vocab, length = length)
    })
    response <- list(y = draw$y, keep = NULL)
    fit <- function(fit_cache, fit_spec) {
      glmmplanr:::fit_via_cpp_cached(response, fit_cache, fit_spec, regularized = TRUE,
                                    re_structure = "maximal", fast_stage = FALSE)
    }
    fits <- list(unadjusted = fit(cache, spec))
    fit_cache <- cache
    fit_cache$ws_env <- NULL
    fit_cache$X_full <- cbind(cache$X_full, covariate_cols(draw$vocab, draw$length, x, FALSE))
    fits$adjusted <- fit(fit_cache, adjusted)
    if (cfg$vocab_by_method != 0) {
      fit_cache$X_full <- cbind(fit_cache$X_full, method_vocabulary = x * draw$vocab)
      fits$moderation <- fit(fit_cache, moderator)
    }
    do.call(rbind, lapply(names(fits), function(model) {
      f <- fits[[model]]
      valid <- glmmplanr:::valid_target_fit(f$est, f$se, f$converged, f$status,
                                           f$prob_est, f$prob_se, "pp")
      est <- if (valid) 100 * f$prob_est else NA_real_
      se <- if (valid) 100 * f$prob_se else NA_real_
      interaction <- valid && model == "moderation"
      data.frame(case = k, replication = i, model = model, valid = valid,
                 estimate_pp = est, se_pp = se,
                 covered = if (valid) abs(est - cfg$effect_pp) <= z_crit * se else NA,
                 detected = valid && abs(est) > z_crit * se,
                 interaction_est = if (interaction) f$beta_full[5L] else NA_real_,
                 interaction_se = if (interaction) sqrt(f$vcov_full[5L, 5L]) else NA_real_,
                 status = f$status, diagnostic = if (valid) "" else f$error_reason)
    }))
  })
  fits <- do.call(rbind, rows)
  cat(sprintf("Covariate case %2d/%d %-14s %2dx%2d %+d pp slope interaction %.2f: %d/%d valid\n",
              k, nrow(cases), cfg$allocation, ns, ni, cfg$effect_pp, cfg$vocab_by_method,
              sum(fits$valid), nrow(fits)))
  flush.console()
  fits
}

summarize <- function(rows) {
  cfg <- cases[rows$case[1L], ]
  valid <- rows$valid
  n <- nrow(rows)
  n_valid <- sum(valid)
  covered <- sum(rows$covered, na.rm = TRUE)
  detected <- sum(rows$detected)
  cov_ci <- if (n_valid) stats::binom.test(covered, n_valid)$conf.int else c(NA, NA)
  det_ci <- stats::binom.test(detected, n)$conf.int
  interaction_covered <- abs(rows$interaction_est - cfg$vocab_by_method) <= z_crit * rows$interaction_se
  cbind(cfg, model = rows$model[1L], attempted = n, valid = n_valid,
        covered = covered, detected = detected, coverage = covered / n_valid,
        coverage_lo = cov_ci[1L], coverage_hi = cov_ci[2L], detection = detected / n,
        detection_lo = det_ci[1L], detection_hi = det_ci[2L],
        bias_pp = mean(rows$estimate_pp[valid]) - cfg$effect_pp,
        sd_pp = stats::sd(rows$estimate_pp[valid]),
        mean_width_pp = 2 * z_crit * mean(rows$se_pp[valid]),
        interaction_mean = mean(rows$interaction_est, na.rm = TRUE),
        interaction_covered = if (any(!is.na(interaction_covered))) {
          sum(interaction_covered, na.rm = TRUE)
        } else NA_integer_)
}

cat(sprintf("Covariates: seed %d, %d DGPs and 28 fits per replication, %d replications\n",
            master_seed, nrow(cases), n_rep))
cat("Package:", normalizePath(find.package("glmmplanr")), "\n")
res <- parallel::mclapply(seq_len(nrow(cases)), run_case, mc.cores = workers,
                          mc.preschedule = FALSE, mc.set.seed = FALSE)
stopifnot(all(vapply(res, is.data.frame, logical(1))))
fits <- do.call(rbind, res)
stopifnot(nrow(cases) == 12L, nrow(fits) == 28L * n_rep,
          !anyDuplicated(fits[c("case", "replication", "model")]))
summary <- do.call(rbind, lapply(split(fits, fits[c("case", "model")], drop = TRUE), summarize))
utils::write.csv(cases, file.path(out, "l2_vocabulary_covariate_scenarios.csv"), row.names = FALSE)
utils::write.csv(summary, file.path(out, "l2_vocabulary_covariate_summary.csv"), row.names = FALSE)
connection <- gzfile(file.path(out, "l2_vocabulary_covariate_fits.csv.gz"), "wt")
utils::write.csv(fits, connection, row.names = FALSE)
close(connection)
cat(sprintf("Completed %d datasets and %d fits.\n", nrow(cases) * n_rep, nrow(fits)))
