item_n <- 24L
zcrit <- qnorm(.975)
q <- glmmplanr:::gh_quad_1d(41L)
mean_prob <- function(mu, sd) as.numeric(plogis(outer(mu, sd * q$nodes, "+")) %*% q$weights)
calibrate <- function(p, sd) {
  mu <- uniroot(function(mu) mean_prob(mu, sd) - p, c(-15, 15), tol = 1e-11)$root
  check <- integrate(function(x) plogis(mu + sd * x) * dnorm(x), -Inf, Inf, rel.tol = 1e-10)$value
  stopifnot(abs(check - p) < 1e-7)
  mu
}
factorial <- list(condition = c("A", "B"))
base <- design_crossed(factorial, allocation = "between_subj")
testlet <- design_testlet(factorial, testlet_size = 4L, between_subj = "condition")
base$Z_item_cells <- base$Z_item_cells[, 1L, drop = FALSE]
testlet$Z_item_cells <- testlet$Z_item_cells[, 1L, drop = FALSE]
make_spec <- function(design) {
  dgp <- calibrate_dgp(design, c(A = .6, B = .7), subj_sd = .7, item_sd = .4, corr_structure = 0,
                       control = calibrate_control(testlet_sd = if (design$has_testlet) c(.5, .4) else NULL))
  specify_target(dgp, "condition", scale = "pp", engine = "bayes_laplace", analysis_re = "uncorrelated",
                 prob_se_method = "joint_delta", control = target_control(on_singular = "flag_only"))
}
base_spec <- make_spec(base)
random_spec <- make_spec(testlet)
stopifnot(identical(as.numeric(base$X_cells[, 2L]), c(-.5, .5)))

schedule <- function(n, P) {
  size <- rep(item_n %/% P, P)
  size[seq_len(item_n %% P)] <- size[seq_len(item_n %% P)] + 1L
  tr <- expand.grid(subj_id = seq_len(n), item_id = seq_len(item_n), KEEP.OUT.ATTRS = FALSE)
  tr$cell_id <- ifelse(tr$subj_id <= n / 2L, 1L, 2L)
  tr$testlet_id <- rep(seq_len(P), size)[tr$item_id]
  stopifnot(all(table(tr$item_id, tr$cell_id) == n/2L), sum(size) == item_n)
  list(trials = tr, n_subj = n, n_item = item_n, n_testlet = P)
}

fixed_setup <- function(sched, common, local = FALSE) {
  P <- sched$n_testlet
  tr <- sched$trials
  spec <- if (local) random_spec else base_spec
  if (local) {
    # The third native grouping factor indexes reader-by-passage here.
    sched$trials$testlet_id <- tr$subj_id + (tr$testlet_id - 1L) * sched$n_subj
    sched$n_testlet <- sched$n_subj * P
    spec$dgp$design$Z_testlet_cells <- spec$dgp$design$Z_testlet_cells[, 1L, drop = FALSE]
  }
  cache <- glmmplanr:::prepare_cpp_schedule_cache(spec, sched)
  cells <- rep(1:2, each = P)
  X <- base$X_cells[cells, , drop = FALSE]
  passages <- contr.sum(P)[rep(seq_len(P), 2L), , drop = FALSE]
  X <- cbind(X, passages)
  if (!common) X <- cbind(X, passages * X[, 2L])
  rownames(X) <- NULL
  idx <- (tr$cell_id - 1L) * P + tr$testlet_id
  cache$X_full <- X[idx, , drop = FALSE]
  spec$dgp$design$X_cells <- X
  for (name in c("Z_subj_cells", "Z_item_cells", if (local) "Z_testlet_cells")) {
    spec$dgp$design[[name]] <- spec$dgp$design[[name]][cells, , drop = FALSE]
  }
  spec$dgp$design$n_cells <- length(cells)
  spec$contrast_weights <- rep(base_spec$contrast_weights, each = P) / P
  stopifnot(qr(X)$rank == ncol(X), ncol(X) == if (common) P + 1L else 2L * P)
  list(spec = spec, cache = cache, structure = "intercept_only")
}

fit_response <- function(y, model) {
  f <- glmmplanr:::fit_via_cpp_cached(list(y = y, keep = NULL), model$cache, model$spec,
                                      regularized = TRUE, re_structure = model$structure, fast_stage = FALSE)
  valid <- glmmplanr:::valid_target_fit(f$est, f$se, f$converged, f$status, f$prob_est, f$prob_se, "pp")
  fitted_passage_slope <- valid && !is.null(f$Sigma_testlet) && nrow(f$Sigma_testlet) == 2L
  list(valid = valid, estimate_pp = if (valid) 100 * f$prob_est else NA_real_,
       se_pp = if (valid) 100 * f$prob_se else NA_real_,
       slope_sd_hat = if (fitted_passage_slope) sqrt(f$Sigma_testlet[2L, 2L]) else NA_real_,
       status = f$status)
}

calibrate_shared_passage_cases <- function(cases) {
  calibration <- t(vapply(seq_len(nrow(cases)), function(k) {
    cfg <- cases[k, ]
    total_sd <- sqrt(.7^2 + .4^2 + .5^2 + cfg$local_sd^2 + (.5 * cfg$slope_sd)^2)
    vapply(c(.60, .60 + cfg$effect_pp / 100), calibrate, numeric(1L), sd = total_sd)
  }, numeric(2L)))
  cases$mu_a <- calibration[, 1L]
  cases$mu_b <- calibration[, 2L]
  cases
}

draw_shared_passage_response <- function(cfg, sched, stream, max_n, max_p) {
  P <- sched$n_testlet
  tr <- sched$trials
  sid <- tr$subj_id
  iid <- tr$item_id
  cid <- tr$cell_id
  x <- c(-.5, .5)
  mu <- c(cfg$mu_a, cfg$mu_b)
  sd_finite <- sqrt(.7^2 + .4^2 + cfg$local_sd^2)
  item_pos <- sid + (iid - 1L) * max_n
  draw <- glmmplanr:::with_rng_stream(stream, {
    person <- .7 * rnorm(max_n)
    item <- .4 * rnorm(item_n)
    text <- .5 * rnorm(max_p)[seq_len(P)]
    slope <- cfg$slope_sd * rnorm(max_p)[seq_len(P)]
    local <- cfg$local_sd * rnorm(max_n * max_p)
    u <- runif(max_n * item_n)
    short_slot <- as.integer(runif(1L) * P) + 1L
    list(person = person, item = item, text = text, slope = slope, local = local, u = u,
         permutation = if (P == 5L) c(seq_len(P)[-short_slot], short_slot) else seq_len(P))
  })
  # Exchangeable passage draws let the shorter fifth block occupy any passage slot.
  pid <- draw$permutation[tr$testlet_id]
  finite <- mean(mean_prob(mu[2L] + draw$text + x[2L] * draw$slope, sd_finite) -
                 mean_prob(mu[1L] + draw$text + x[1L] * draw$slope, sd_finite)) * 100
  eta <- mu[cid] + draw$person[sid] + draw$item[iid] + draw$text[pid] + x[cid] * draw$slope[pid] +
    draw$local[sid + (pid - 1L) * max_n]
  y <- as.integer(draw$u[item_pos] < plogis(eta))
  list(y = y, truth_sampled_pp = finite, items = draw$item, text = draw$text[draw$permutation],
       slope = draw$slope[draw$permutation])
}
