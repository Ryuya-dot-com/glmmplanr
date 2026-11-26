# Screening and allocation precede any learning or posttest outcome.
library(glmmplanr)

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1L] else "validation"
n_rep <- if (length(args) > 1L) as.integer(args[2L]) else 2048L
workers <- if (length(args) > 2L) as.integer(args[3L]) else
  if (.Platform$OS.type == "windows") 1L else 4L
stopifnot(dir.exists(out), n_rep > 0L, n_rep <= 2048L, workers > 0L)
master_seed <- if (n_rep == 2048L) 87761001L else 87762001L
streams <- glmmplanr:::make_lecuyer_streams(master_seed, n_rep)
ns <- 96L
ni <- 48L
quota <- 16L
z_crit <- stats::qnorm(.975)
profiles <- data.frame(
  profile = c("low", "strong", "ceiling", "learning", "moderation"),
  known_low = c(.10, .10, .30, .10, .10), known_high = c(.20, .50, .70, .50, .50),
  vocab_new = c(0, 0, 0, .45, .45), vocab_by_method = c(0, 0, 0, 0, .35)
)
screen <- data.frame(screen = c("perfect", "imperfect", "weaker"),
                     sensitivity = c(1, .90, .75), specificity = c(1, .95, .85))
conditions <- merge(profiles, screen)
conditions <- subset(conditions, screen != "weaker" | profile %in% c("ceiling", "moderation"))
cases <- merge(conditions, data.frame(effect_pp = c(0, 8)))
cases <- cases[order(match(cases$screen, screen$screen), match(cases$profile, profiles$profile),
                     cases$effect_pp), ]
rownames(cases) <- NULL
cases$case <- seq_len(nrow(cases))
cases$n_rep <- n_rep
cases$master_seed <- master_seed
stopifnot(nrow(cases) == 24L, sum(cases$screen != "perfect") == 14L)

design <- design_crossed(list(method = c("restudy", "retrieval")))
variance <- diag(c(.8, .25)^2)
word_variance <- diag(c(.6, .25)^2)
latent_sd <- sqrt(.8^2 + .6^2 + .5^2 * 2 * .25^2)
fit_dgp <- calibrate_dgp(design, c(restudy = .42, retrieval = .50),
                         subj_sd = c(.8, .25), item_sd = c(.6, .25), corr_structure = 0)
fit_spec <- specify_target(fit_dgp, "method", scale = "pp", engine = "bayes_laplace",
                           analysis_re = "maximal", prob_se_method = "joint_delta",
                           control = target_control(on_singular = "flag_only"))

population <- function(n) {
  q <- glmmplanr:::gh_quad_1d(n)
  i <- rep(seq_len(n), times = n)
  j <- rep(seq_len(n), each = n)
  list(vocab = q$nodes[i], ease = q$nodes[j], weight = q$weights[i] * q$weights[j],
       node = q$nodes, gh_weight = q$weights)
}
fine <- population(31L)
standard <- population(13L)
norm_logit_mean <- function(eta) {
  drop(stats::plogis(outer(eta, latent_sd * fine$node, "+")) %*% fine$gh_weight)
}
root <- function(expectation, truth) {
  stats::uniroot(function(mu) expectation(mu) - truth, c(-15, 15), tol = 1e-11)$root
}
known_mu <- root(norm_logit_mean, .90)

calibrate <- function(cfg) {
  k_at <- vapply(c(cfg$known_low, cfg$known_high), function(target) {
    root(function(mu) sum(fine$gh_weight * stats::plogis(mu + .70 * fine$node)), target)
  }, numeric(1))
  k_coeff <- c(mean(k_at), diff(k_at) / 2)
  k_prob <- function(g) stats::plogis(k_coeff[1L] + k_coeff[2L] * g$vocab + .70 * g$ease)
  population_known <- sum(fine$weight * k_prob(fine))
  w <- fine$weight * (1 - k_prob(fine)) / (1 - population_known)
  mu <- vapply(c(-.5, .5), function(x) {
    offset <- (cfg$vocab_new + x * cfg$vocab_by_method) * fine$vocab + .35 * fine$ease
    target <- .42 + (x + .5) * cfg$effect_pp / 100
    root(function(b) sum(w * norm_logit_mean(b + offset)), target)
  }, numeric(1))
  standard_unknown <- standard$weight * (1 - k_prob(standard))
  standard_unknown <- standard_unknown / sum(standard_unknown)
  arm_diff <- function(v) {
    wv <- fine$gh_weight * (1 - stats::plogis(k_coeff[1L] + k_coeff[2L] * v + .70 * fine$node))
    eta <- .35 * fine$node + cfg$vocab_new * v
    sum(wv * (norm_logit_mean(mu[2L] + eta + .5 * cfg$vocab_by_method * v) -
                norm_logit_mean(mu[1L] + eta - .5 * cfg$vocab_by_method * v))) / sum(wv)
  }
  list(k_coeff = k_coeff, population_known = population_known, mu = mu,
       standard_unknown = standard_unknown,
       subgroup_pp = 100 * vapply(c(-1, 1), arm_diff, numeric(1)))
}
truths <- lapply(seq_len(nrow(cases)), function(k) calibrate(cases[k, ]))
cases$population_known <- vapply(truths, function(t) t$population_known, numeric(1))
cases$truth_all_pp <- (1 - cases$population_known) * cases$effect_pp
cases$truth_new_low_pp <- vapply(truths, function(t) t$subgroup_pp[1L], numeric(1))
cases$truth_new_high_pp <- vapply(truths, function(t) t$subgroup_pp[2L], numeric(1))

# The linked AB, BC, CD, DE, EF, FA lists give every word 16 exposures per arm.
initial_words <- lapply(seq_len(ns), function(i) {
  form <- (i - 1L) %/% 16L + 1L
  as.integer(c((form - 1L) * 8L, (form %% 6L) * 8L) |> rep(each = 8L)) + rep(seq_len(8L), 2L)
})
initial_items <- unlist(initial_words, use.names = FALSE)
initial_subj <- rep(seq_len(ns), each = quota)
initial_pair <- initial_subj + ns * (initial_items - 1L)
initial_x <- rep(rep(c(-.5, .5), each = quota / 2L), ns) * rep(rep(c(1, -1), each = 8L), each = quota)
stopifnot(length(initial_x) == length(initial_pair), !anyDuplicated(initial_pair),
          all(table(initial_subj, initial_x) == 8L), all(table(initial_items, initial_x) == 16L))
other_words <- lapply(initial_words, function(w) setdiff(seq_len(ni), w))
pair_subj <- rep(seq_len(ns), times = ni)
pair_item <- rep(seq_len(ni), each = ns)

case_dir <- file.path(out, sprintf("l2_vocabulary_prior_knowledge_checkpoints_%d", master_seed))
dir.create(case_dir, showWarnings = FALSE)
manifest <- file.path(case_dir, "scenarios.rds")
if (file.exists(manifest)) {
  stopifnot(identical(readRDS(manifest), cases))
} else {
  saveRDS(cases, manifest)
}
utils::write.csv(cases, file.path(out, "l2_vocabulary_prior_knowledge_scenarios.csv"), row.names = FALSE)

new_matrix <- function(x, vocab, ease) {
  cbind(1, x, vocab, ease, x * vocab)
}

run_case <- function(k) {
  cfg <- cases[k, ]
  path <- file.path(case_dir, sprintf("case_%03d.rds", k))
  if (file.exists(path)) return(path)
  truth <- truths[[k]]
  reference <- NULL
  if (cfg$screen != "perfect") {
    ref_case <- cases$case[cases$profile == cfg$profile & cases$effect_pp == cfg$effect_pp &
                           cases$screen == "perfect"]
    ref_fits <- readRDS(file.path(case_dir, sprintf("case_%03d.rds", ref_case)))$fits
    reference <- lapply(split(ref_fits, ref_fits$model), function(rows) {
      stopifnot(identical(rows$replication, seq_len(n_rep)))
      rows[setdiff(names(rows), c("case", "replication"))]
    })
  }
  n_grid <- length(standard$weight)
  marginal <- design
  marginal$X_cells <- new_matrix(rep(c(-.5, .5), each = n_grid),
                                  rep(standard$vocab, 2L), rep(standard$ease, 2L))
  marginal$n_cells <- nrow(marginal$X_cells)
  marginal$Z_subj_cells <- marginal$Z_item_cells <- marginal$X_cells[, 1:2, drop = FALSE]
  new_weights <- c(-truth$standard_unknown, truth$standard_unknown)
  true_beta <- c(mean(truth$mu), diff(truth$mu), cfg$vocab_new, .35, cfg$vocab_by_method)
  check <- glmmplanr:::estimate_prob_effect_delta(
    true_beta, matrix(0, 5L, 5L), variance, word_variance, NULL, marginal, 0, new_weights
  )$est
  stopifnot(abs(100 * check - cfg$effect_pp) < 1e-6)

  native_fit <- function(pair, x, y, vocab, ease, model) {
    is_all <- model == "fixed_all"
    z <- cbind(1, x)
    X <- if (is_all) z else new_matrix(x, vocab[pair_subj[pair]], ease[pair_item[pair]])
    valid <- FALSE
    status <- "unsupported_design"
    diagnostic <- "The retained observations do not identify the fixed effects."
    estimate <- se <- NA_real_
    if (nrow(X) > ncol(X) && qr(X)$rank == ncol(X)) {
      df <- data.frame(y = y, subj_id = pair_subj[pair], item_id = pair_item[pair])
      fit <- glmmplanr:::fit_via_cpp(df, X, z, z, matrix(0, length(y), 0L),
                                      0, 1L, TRUE, fit_spec, FALSE, 2L)
      status <- fit$status
      diagnostic <- fit$error_reason
      if (isTRUE(fit$converged)) {
        point <- glmmplanr:::estimate_prob_effect_delta(
          fit$beta_full, fit$vcov_full, fit$Sigma_subj, fit$Sigma_item, NULL,
          if (is_all) design else marginal, 0,
          if (is_all) fit_spec$contrast_weights else new_weights,
          var_log_tau_subj = fit$var_log_tau_subj, var_log_tau_item = fit$var_log_tau_item,
          J_beta_theta = fit$J_beta_theta
        )
        valid <- glmmplanr:::valid_target_fit(fit$est, fit$se, fit$converged, fit$status,
                                             point$est, point$se, "pp")
        if (valid) {
          estimate <- 100 * point$est
          se <- 100 * point$se
        }
      }
    }
    target <- if (is_all) cfg$truth_all_pp else cfg$effect_pp
    data.frame(model = model, truth_pp = target, valid = valid, estimate_pp = estimate, se_pp = se,
               covered = if (valid) abs(estimate - target) <= z_crit * se else NA,
               detected = valid && abs(estimate) > z_crit * se, status = status,
               diagnostic = if (valid) "" else diagnostic)
  }

  rows <- lapply(seq_len(n_rep), function(i) {
    draw <- glmmplanr:::with_rng_stream(streams[[i]], {
      subj <- matrix(stats::rnorm(3L * ns), ns, 3L)
      item <- matrix(stats::rnorm(3L * ni), ni, 3L)
      k_eta <- truth$k_coeff[1L] + truth$k_coeff[2L] * subj[pair_subj, 1L] + .70 * item[pair_item, 1L]
      known <- stats::runif(ns * ni) < stats::plogis(k_eta)
      pre_u <- stats::runif(ns * ni)
      list(subj = subj, item = item, known = known,
           screen_known = pre_u < ifelse(known, cfg$sensitivity, 1 - cfg$specificity),
           post_u = stats::runif(ns * ni), screen_order = matrix(stats::runif(ns * ni), ns, ni),
           assignment_order = stats::runif(ns * ni))
    })
    vocab <- draw$subj[, 1L]
    ease <- draw$item[, 1L]
    high <- rank(vocab) > 3L * ns / 4L
    post <- function(pair, x) {
      s <- pair_subj[pair]
      j <- pair_item[pair]
      mu <- truth$mu[1L + (x > 0)] +
        (cfg$vocab_new + cfg$vocab_by_method * x) * vocab[s] + .35 * ease[j]
      mu[draw$known[pair]] <- known_mu
      eta <- mu + .8 * draw$subj[s, 2L] + .6 * draw$item[j, 2L] +
        .25 * x * (draw$subj[s, 3L] + draw$item[j, 3L])
      as.integer(draw$post_u[pair] < stats::plogis(eta))
    }

    screened <- !draw$screen_known[initial_pair]
    initial_y <- post(initial_pair, initial_x)
    topup <- lapply(seq_len(ns), function(s) {
      others <- other_words[[s]]
      others <- others[order(draw$screen_order[s, others])]
      candidates <- s + ns * (c(initial_words[[s]], others) - 1L)
      eligible <- which(!draw$screen_known[candidates])
      inspected <- if (length(eligible) >= quota) eligible[quota] else ni
      selected <- candidates[head(eligible, quota)]
      random_order <- order(draw$assignment_order[selected])
      x <- numeric(length(selected))
      x[random_order] <- rep(c(-.5, .5), length.out = length(selected)) *
        if (draw$screen_order[s, initial_words[[s]][1L]] < .5) -1 else 1
      list(pair = selected, x = x, inspected = inspected)
    })
    top_pair <- unlist(lapply(topup, `[[`, "pair"), use.names = FALSE)
    top_x <- unlist(lapply(topup, `[[`, "x"), use.names = FALSE)
    top_screen <- vapply(topup, `[[`, integer(1), "inspected")
    top_y <- post(top_pair, top_x)
    stopifnot(!anyDuplicated(top_pair), all(!draw$screen_known[top_pair]),
              length(top_pair) <= quota * ns, all(top_screen >= quota))
    exposure <- function(policy, pair, x, inspected) {
      u <- !draw$known[pair]
      s <- pair_subj[pair]
      arm_counts <- table(factor(s[u], levels = seq_len(ns)), factor(x[u], levels = c(-.5, .5)))
      item_counts <- table(factor(pair_item[pair[u]], levels = seq_len(ni)),
                           factor(x[u], levels = c(-.5, .5)))
      taught_by_s <- tabulate(s, nbins = ns)
      data.frame(policy = policy, screened_n = sum(inspected), taught_n = length(pair),
                 unknown_n = sum(u), known_n = sum(!u),
                 unknown_a = sum(u & x < 0), unknown_b = sum(u & x > 0),
                 upper_screened_n = sum(inspected[high]), upper_taught_n = sum(taught_by_s[high]),
                 upper_unknown_n = sum(u & high[s]),
                 shortfall_n = sum(taught_by_s < quota),
                 upper_shortfall_n = sum(taught_by_s[high] < quota),
                 no_unknown_one_arm = sum(rowSums(arm_counts > 0) < 2),
                 upper_no_unknown_one_arm = sum(rowSums(arm_counts[high, ] > 0) < 2),
                 shared_unknown_words = sum(rowSums(item_counts > 0) == 2),
                 min_unknown_word_arm = min(item_counts))
    }
    counts <- rbind(exposure("fixed_all", initial_pair, initial_x, rep(quota, ns)),
                    exposure("fixed_screened", initial_pair[screened], initial_x[screened],
                             rep(quota, ns)),
                    exposure("topup_screened", top_pair, top_x, top_screen))

    all_fit <- if (is.null(reference)) {
      native_fit(initial_pair, initial_x, initial_y, vocab, ease, "fixed_all")
    } else reference$fixed_all[i, ]
    rows <- rbind(all_fit, native_fit(initial_pair[screened], initial_x[screened], initial_y[screened],
                             vocab, ease, "fixed_screened"),
                  native_fit(top_pair, top_x, top_y, vocab, ease, "topup_screened"))
    if (!is.null(reference)) {
      fixed_oracle <- reference$fixed_screened[i, ]
      fixed_oracle$model <- "fixed_oracle"
      top_unknown <- !draw$known[top_pair]
      rows <- rbind(rows, fixed_oracle,
                    native_fit(top_pair[top_unknown], top_x[top_unknown], top_y[top_unknown],
                               vocab, ease, "topup_oracle"))
    }
    list(fits = cbind(case = k, replication = i, rows),
         exposure = cbind(case = k, replication = i, counts))
  })
  result <- list(fits = do.call(rbind, lapply(rows, `[[`, "fits")),
                 exposure = do.call(rbind, lapply(rows, `[[`, "exposure")))
  tmp <- paste0(path, ".new")
  saveRDS(result, tmp)
  stopifnot(file.rename(tmp, path))
  cat(sprintf("Case %02d/%d %-10s %-10s %+d pp, %d/%d valid\n",
              k, nrow(cases), cfg$profile, cfg$screen, cfg$effect_pp,
              sum(result$fits$valid), nrow(result$fits)))
  flush.console()
  path
}

exact_ci <- function(k, n) {
  if (!n) return(c(NA_real_, NA_real_))
  stats::binom.test(k, n)$conf.int[1:2]
}

summarize <- function(rows) {
  cfg <- cases[rows$case[1L], ]
  n <- nrow(rows)
  valid <- rows$valid
  n_valid <- sum(valid)
  covered <- sum(rows$covered, na.rm = TRUE)
  detected <- sum(rows$detected)
  cc <- exact_ci(covered, n_valid)
  dc <- exact_ci(detected, n)
  coverage_status <- if (n_valid / n < .90 || cc[2L] < .90) "fail" else
    if (cc[1L] >= .90) "pass" else "indeterminate"
  detection_status <- if (cfg$effect_pp == 0) {
    if (dc[2L] <= .075) "pass" else if (dc[1L] > .075) "fail" else "indeterminate"
  } else if (dc[1L] >= .80) "pass" else if (dc[2L] < .80) "fail" else "indeterminate"
  err <- rows$estimate_pp[valid] - rows$truth_pp[1L]
  cbind(cfg, model = rows$model[1L], truth_pp = rows$truth_pp[1L], attempted = n, valid = n_valid,
        covered = covered, detected = detected, coverage = covered / n_valid,
        coverage_lo = cc[1L], coverage_hi = cc[2L], coverage_status = coverage_status,
        detection = detected / n, detection_lo = dc[1L], detection_hi = dc[2L],
        detection_status = detection_status, bias_pp = mean(err),
        bias_mcse_pp = stats::sd(err) / sqrt(n_valid), mean_estimate_pp = mean(rows$estimate_pp[valid]),
        empirical_sd_pp = stats::sd(err), mean_width_pp = 2 * z_crit * mean(rows$se_pp[valid]))
}

cat(sprintf("Master seed %d, %d scenarios × %d replications, %d workers\n",
            master_seed, nrow(cases), n_rep, workers))
cat("Package:", normalizePath(find.package("glmmplanr")), "\n")
run_batch <- function(ids) {
  paths <- parallel::mclapply(ids, run_case, mc.cores = workers,
                              mc.preschedule = FALSE, mc.set.seed = FALSE)
  stopifnot(all(vapply(paths, function(x) is.character(x) && length(x) == 1L && file.exists(x),
                        logical(1))))
  paths
}
paths <- c(run_batch(which(cases$screen == "perfect")), run_batch(which(cases$screen != "perfect")))
results <- lapply(paths, readRDS)
fits <- do.call(rbind, lapply(results, `[[`, "fits"))
exposures <- do.call(rbind, lapply(results, `[[`, "exposure"))
stopifnot(nrow(fits) == 100L * n_rep, nrow(exposures) == 72L * n_rep,
          !anyDuplicated(fits[c("case", "replication", "model")]),
          !anyDuplicated(exposures[c("case", "replication", "policy")]))
summary <- do.call(rbind, lapply(results, function(r) {
  do.call(rbind, lapply(split(r$fits, r$fits$model), summarize))
}))
counts <- do.call(rbind, lapply(split(exposures, exposures[c("case", "policy")], drop = TRUE),
  function(rows) {
    keys <- rows[1L, c("case", "policy")]
    cbind(keys, as.list(colMeans(rows[setdiff(names(rows), c("case", "policy", "replication"))])))
  }))
utils::write.csv(summary, file.path(out, "l2_vocabulary_prior_knowledge_summary.csv"), row.names = FALSE)
utils::write.csv(counts, file.path(out, "l2_vocabulary_prior_knowledge_exposure.csv"), row.names = FALSE)
for (name in c("fits", "exposures")) {
  connection <- gzfile(file.path(out, paste0("l2_vocabulary_prior_knowledge_", name, ".csv.gz")), "wt")
  utils::write.csv(get(name), connection, row.names = FALSE)
  close(connection)
}
cat(sprintf("Completed %d datasets and %d fit records.\n", nrow(cases) * n_rep, nrow(fits)))
