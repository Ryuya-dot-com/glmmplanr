# Same scheduled response budget with full, linked sparse, and disjoint sparse layouts.
library(glmmplanr)

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1L] else "validation"
n_rep <- if (length(args) > 1L) as.integer(args[2L]) else 3072L
workers <- if (length(args) > 2L) as.integer(args[3L]) else
  if (.Platform$OS.type == "windows") 1L else 4L
stopifnot(dir.exists(out), n_rep > 0, n_rep <= 3072L, workers > 0)
n_subj <- 96L
max_words <- 48L
z_crit <- stats::qnorm(.975)
seed <- if (n_rep == 3072L) 86614001L else 86615001L
streams <- glmmplanr:::make_lecuyer_streams(seed, n_rep)

forms <- list(
  full = list(c(1L, 2L)),
  linked = lapply(seq_len(6L), function(i) c(i, i %% 6L + 1L)),
  disjoint = list(c(1L, 2L), c(3L, 4L), c(5L, 6L))
)
profiles <- data.frame(
  profile = c("baseline", "learner_dominant", "word_dominant", "match", "moderation",
              "allocation_bias", "mcar", "match_mar", "outcome_mnar"),
  slope_subj = c(.25, .40, .15, rep(.25, 6L)),
  slope_item = c(.25, .15, .40, rep(.25, 6L)),
  match_logit = c(0, 0, 0, rep(.65, 6L)),
  moderation_logit = c(0, 0, 0, 0, .60, 0, 0, 0, 0),
  matched_version_share = c(rep(.5, 5L), .75, rep(.5, 3L)),
  observation_rule = c(rep("complete", 6L), "mcar", "match_mar", "outcome_mnar")
)
cases <- merge(merge(data.frame(layout = names(forms)), profiles), data.frame(effect_pp = c(0, 8)))
cases <- cases[order(match(cases$profile, profiles$profile), cases$effect_pp,
                     match(cases$layout, names(forms))), ]
rownames(cases) <- NULL
cases$case <- seq_len(nrow(cases))
cases$n_subj <- n_subj
cases$n_item <- ifelse(cases$layout == "full", 16L, 48L)
cases$n_scheduled <- n_subj * 16L
cases$p_restudy <- .42
cases$p_retrieval <- .42 + cases$effect_pp / 100
cases$master_seed <- seed
cases$n_rep <- n_rep
cases$fit_adjusted <- cases$profile %in% c("match", "moderation", "allocation_bias",
                                           "match_mar", "outcome_mnar")
stopifnot(nrow(cases) == 54L, sum(cases$fit_adjusted) == 30L)

make_schedule <- function(cfg) {
  edges <- forms[[cfg$layout]]
  n_form <- n_subj %/% length(edges)
  half <- n_form %/% 2L
  rows <- lapply(seq_along(edges), function(k) {
    blocks <- edges[[k]]
    items <- c((blocks[1L] - 1L) * 8L + seq_len(8L), (blocks[2L] - 1L) * 8L + seq_len(8L))
    g <- expand.grid(form_subj = seq_len(n_form), item_id = items)
    g$subj_id <- (k - 1L) * n_form + g$form_subj
    subgroup <- ifelse(g$form_subj <= half, 1, -1)
    word_type <- ifelse(((g$item_id - 1L) %/% 8L) %% 2L == 0L, 1, -1)
    subgroup_position <- (g$form_subj - 1L) %% half + 1L
    version <- ifelse(subgroup_position <= half * cfg$matched_version_share, subgroup, -subgroup)
    g$cell_id <- ifelse(version * word_type == 1, 2L, 1L)
    g$match_score <- subgroup * word_type
    g$form <- k
    g
  })
  trials <- do.call(rbind, rows)
  both <- table(trials$cell_id, trials$match_score)
  exposure <- nrow(trials) / (2L * cfg$n_item)
  stopifnot(nrow(trials) == cfg$n_scheduled, all(table(trials$subj_id, trials$cell_id) == 8L),
            all(table(trials$item_id, trials$cell_id) == exposure),
            !anyDuplicated(trials[c("subj_id", "item_id")]),
            both["2", "1"] / sum(both["2", ]) == cfg$matched_version_share,
            both["1", "1"] / sum(both["1", ]) == 1 - cfg$matched_version_share)
  list(trials = trials, n_subj = cfg$n_subj, n_item = cfg$n_item)
}

design <- design_crossed(list(method = c("restudy", "retrieval")))
match_design <- design
cells <- rep(1:2, each = 2L)
match_scores <- rep(c(-1, 1), 2L)
X <- design$X_cells[cells, , drop = FALSE]
match_design$X_cells <- cbind(X, match_score = match_scores, method_match = X[, 2L] * match_scores)
match_design$Z_subj_cells <- design$Z_subj_cells[cells, , drop = FALSE]
match_design$Z_item_cells <- design$Z_item_cells[cells, , drop = FALSE]
match_design$n_cells <- length(cells)

calibrate_truth <- function(cfg) {
  x <- c(-.5, .5)
  targets <- c(cfg$p_restudy, cfg$p_retrieval)
  sd <- sqrt(.8^2 + .6^2 + x^2 * (cfg$slope_subj^2 + cfg$slope_item^2))
  match_term <- cfg$match_logit + cfg$moderation_logit * x
  integrate_arm <- function(mu, k, m) {
    stats::integrate(function(z) stats::plogis(mu + m * match_term[k] + sd[k] * z) *
                       stats::dnorm(z), -Inf, Inf, rel.tol = 1e-10)$value
  }
  mu <- vapply(seq_len(2L), function(k) {
    stats::uniroot(function(b) {
      mean(vapply(c(-1, 1), function(m) integrate_arm(b, k, m), numeric(1))) - targets[k]
    },
                   c(-12, 12), tol = 1e-11)$root
  }, numeric(1))
  conditional <- vapply(seq_along(cells), function(i) {
    integrate_arm(mu[cells[i]], cells[i], match_scores[i])
  }, numeric(1))
  stopifnot(max(abs(rowMeans(matrix(conditional, nrow = 2L, byrow = TRUE)) - targets)) < 1e-8)
  list(mu = mu, conditional = conditional)
}

case_dir <- file.path(out, sprintf("l2_vocabulary_sparse_checkpoints_%d", seed))
dir.create(case_dir, showWarnings = FALSE)
manifest <- file.path(case_dir, "scenarios.rds")
if (file.exists(manifest)) {
  stopifnot(identical(readRDS(manifest), cases))
} else {
  saveRDS(cases, manifest)
}
utils::write.csv(cases, file.path(out, "l2_vocabulary_sparse_scenarios.csv"), row.names = FALSE)

run_case <- function(k) {
  cfg <- cases[k, ]
  path <- file.path(case_dir, sprintf("case_%03d.rds", k))
  if (file.exists(path)) return(path)
  schedule <- make_schedule(cfg)
  tr <- schedule$trials
  n_components <- length(unique(glmmplanr:::item_components(tr, cfg$n_item)))
  stopifnot(n_components == if (cfg$layout == "disjoint") 3L else 1L)
  truth <- calibrate_truth(cfg)
  dgp <- calibrate_dgp(design, c(restudy = cfg$p_restudy, retrieval = cfg$p_retrieval),
                       subj_sd = c(.8, cfg$slope_subj), item_sd = c(.6, cfg$slope_item),
                       corr_structure = 0)
  spec <- specify_target(dgp, "method", scale = "pp", engine = "bayes_laplace",
                         analysis_re = "maximal", prob_se_method = "joint_delta",
                         control = target_control(on_singular = "flag_only"))
  cache <- glmmplanr:::prepare_cpp_schedule_cache(spec, schedule)
  adjusted_spec <- spec
  adjusted_spec$dgp$design <- match_design
  adjusted_spec$contrast_weights <- rep(spec$contrast_weights / 2, each = 2L)
  x <- cache$X_full[, 2L]
  m <- tr$match_score
  is_b <- tr$cell_id == 2L
  pairs <- tr$subj_id + n_subj * (tr$item_id - 1L)
  targets <- c(cfg$p_restudy, cfg$p_retrieval)[tr$cell_id]
  fixed <- truth$mu[tr$cell_id] + (cfg$match_logit + cfg$moderation_logit * x) * m
  adjusted_cache <- cache
  adjusted_cache$X_full <- cbind(cache$X_full, match_score = m, method_match = x * m)
  adjusted_cache$cell_id <- 2L * (tr$cell_id - 1L) + ifelse(m == 1, 2L, 1L)
  adjusted_cache$ws_env <- new.env(parent = emptyenv())
  stopifnot(identical(unname(adjusted_cache$X_full),
                      unname(match_design$X_cells[adjusted_cache$cell_id, ])))
  beta <- c(mean(truth$mu), diff(truth$mu), cfg$match_logit, cfg$moderation_logit)
  adjusted_truth <- glmmplanr:::estimate_prob_effect_delta(
    beta, matrix(0, 4L, 4L), dgp$Sigma_subj, dgp$Sigma_item, NULL,
    match_design, 0, adjusted_spec$contrast_weights, n_quad = dgp$n_quad
  )$est
  stopifnot(abs(100 * adjusted_truth - cfg$effect_pp) < 1e-6)

  fit_row <- function(response, model, model_cache, model_spec) {
    f <- glmmplanr:::fit_via_cpp_cached(response, model_cache, model_spec,
                                        regularized = TRUE, re_structure = "maximal",
                                        fast_stage = FALSE)
    valid <- glmmplanr:::valid_target_fit(f$est, f$se, f$converged, f$status,
                                          f$prob_est, f$prob_se, "pp")
    estimate <- if (valid) 100 * f$prob_est else NA_real_
    se <- if (valid) 100 * f$prob_se else NA_real_
    data.frame(model = model, valid = valid, estimate_pp = estimate, se_pp = se,
               covered = if (valid) abs(estimate - cfg$effect_pp) <= z_crit * se else NA,
               detected = valid && abs(estimate) > z_crit * se, status = f$status,
               diagnostic = if (valid) "" else f$error_reason)
  }
  rows <- lapply(seq_len(n_rep), function(i) {
    sim <- glmmplanr:::with_rng_stream(streams[[i]], {
      zs <- matrix(stats::rnorm(2L * n_subj), n_subj, 2L)[tr$subj_id, , drop = FALSE]
      zi <- matrix(stats::rnorm(2L * max_words), max_words, 2L)[tr$item_id, , drop = FALSE]
      y_u <- stats::runif(n_subj * max_words)[pairs]
      keep_u <- stats::runif(n_subj * max_words)[pairs]
      eta <- fixed + .8 * zs[, 1L] + .6 * zi[, 1L] +
        (cfg$slope_subj * zs[, 2L] + cfg$slope_item * zi[, 2L]) * x
      y <- as.integer(y_u < stats::plogis(eta))
      keep_prob <- switch(cfg$observation_rule, complete = 1, mcar = .85,
                          match_mar = .85 + .12 * m * (2 * x),
                          outcome_mnar = .85 + .12 * (y - targets) * (2 * x))
      list(y = y, keep = if (cfg$observation_rule == "complete") seq_along(y) else
        which(keep_u < keep_prob))
    })
    keep <- sim$keep
    y <- sim$y[keep]
    b <- is_b[keep]
    dat <- list(y = y, keep = if (length(keep) == length(is_b)) NULL else keep)
    fitted <- fit_row(dat, "unadjusted", cache, spec)
    if (cfg$fit_adjusted) {
      fitted <- rbind(fitted, fit_row(dat, "standardized_match", adjusted_cache, adjusted_spec))
    }
    cbind(case = k, replication = i, observed_a = sum(!b), observed_b = sum(b),
          recalled_a = sum(y[!b]), recalled_b = sum(y[b]),
          matched_a = sum(m[keep][!b] == 1), matched_b = sum(m[keep][b] == 1), fitted)
  })
  fits <- do.call(rbind, rows)
  stopifnot(nrow(fits) == n_rep * (1L + cfg$fit_adjusted),
            !anyDuplicated(fits[c("replication", "model")]))
  tmp <- paste0(path, ".new")
  saveRDS(fits, tmp)
  stopifnot(file.rename(tmp, path))
  cat(sprintf("Case %02d/%d %-17s %-8s %+d pp, %d/%d valid\n", k, nrow(cases),
              cfg$profile, cfg$layout, cfg$effect_pp, sum(fits$valid), nrow(fits)))
  flush.console()
  path
}

exact_ci <- function(k, n) {
  if (!n) return(c(NA_real_, NA_real_))
  stats::binom.test(k, n)$conf.int[1:2]
}

summarize <- function(fits) {
  cfg <- cases[fits$case[1L], ]
  n <- nrow(fits)
  valid <- fits$valid
  n_valid <- sum(valid)
  covered <- sum(fits$covered, na.rm = TRUE)
  detected <- sum(fits$detected)
  cov_ci <- exact_ci(covered, n_valid)
  det_ci <- exact_ci(detected, n)
  coverage_status <- if (n_valid / n < .90 || cov_ci[2L] < .90) "fail" else
    if (cov_ci[1L] >= .90) "pass" else "indeterminate"
  detection_status <- if (cfg$effect_pp == 0) {
    if (det_ci[2L] <= .075) "pass" else if (det_ci[1L] > .075) "fail" else "indeterminate"
  } else if (det_ci[1L] >= .80) "pass" else if (det_ci[2L] < .80) "fail" else "indeterminate"
  err <- fits$estimate_pp[valid] - cfg$effect_pp
  cbind(cfg, model = fits$model[1L], attempted = n, valid = n_valid, covered = covered,
        detected = detected, coverage = covered / n_valid, coverage_lo = cov_ci[1L],
        coverage_hi = cov_ci[2L], coverage_status = coverage_status, detection = detected / n,
        detection_lo = det_ci[1L], detection_hi = det_ci[2L], detection_status = detection_status,
        bias_pp = mean(err), bias_mcse_pp = stats::sd(err) / sqrt(n_valid), sd_pp = stats::sd(err),
        mean_se_pp = mean(fits$se_pp[valid]), mean_width_pp = 2 * z_crit * mean(fits$se_pp[valid]),
        mean_observed = mean(fits$observed_a + fits$observed_b),
        matched_a = mean(fits$matched_a / fits$observed_a),
        matched_b = mean(fits$matched_b / fits$observed_b),
        raw_diff_pp = 100 * mean(fits$recalled_b / fits$observed_b -
                                  fits$recalled_a / fits$observed_a))
}

cat(sprintf("Seed %d, %d scenarios x %d replications, %d workers\n",
            seed, nrow(cases), n_rep, workers))
cat("Package:", normalizePath(find.package("glmmplanr")), "\n")
paths <- parallel::mclapply(seq_len(nrow(cases)), run_case, mc.cores = workers,
                            mc.preschedule = FALSE, mc.set.seed = FALSE)
stopifnot(all(vapply(paths, function(x) is.character(x) && length(x) == 1L &&
                      file.exists(x), logical(1))))
result <- lapply(paths, readRDS)
fits <- do.call(rbind, result)
stopifnot(nrow(fits) == (nrow(cases) + sum(cases$fit_adjusted)) * n_rep,
          !anyDuplicated(fits[c("case", "replication", "model")]))
summary <- do.call(rbind, lapply(result, function(rows) {
  do.call(rbind, lapply(split(rows, rows$model), summarize))
}))
utils::write.csv(summary, file.path(out, "l2_vocabulary_sparse_summary.csv"), row.names = FALSE)
connection <- gzfile(file.path(out, "l2_vocabulary_sparse_fits.csv.gz"), "wt")
utils::write.csv(fits, connection, row.names = FALSE)
close(connection)
cat(sprintf("Completed %d datasets and %d fits.\n", nrow(cases) * n_rep, nrow(fits)))
