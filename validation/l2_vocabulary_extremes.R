# Binary meaning recall under global and concentrated floor/ceiling effects.
library(glmmplanr)

args <- commandArgs(TRUE)
out <- if (length(args)) args[1L] else "validation"
stage <- if (length(args) > 1L) match.arg(args[2L], c("pilot", "large")) else "pilot"
workers <- if (length(args) > 2L) as.integer(args[3L]) else
  if (.Platform$OS.type == "windows") 1L else 4L
stopifnot(dir.exists(out), workers > 0L)
n_rep <- if (stage == "pilot") 48L else 2048L
seed <- if (stage == "pilot") 88741001L else 88742001L
streams <- glmmplanr:::make_lecuyer_streams(seed, n_rep)
z <- qnorm(.975)

layouts <- rbind(transform(expand.grid(n_subj = c(48L, 96L, 144L),
                                       n_item = c(16L, 32L, 48L)), m = 16L),
                 data.frame(n_subj = c(96L, 96L, 32L, 48L, 64L),
                            n_item = c(48L, 48L, 48L, 32L, 24L),
                            m = c(8L, 24L, 48L, 32L, 24L)))
layouts$id <- with(layouts, paste(n_subj, n_item, m, sep = "_"))
core <- c("48_48_16", "96_48_16", "96_16_16")
paired <- core[-1L]
local <- c(core, "144_48_16", "32_48_48", "96_48_8")

scenarios <- function(phase, ids, p = .42, effect = c("null", "pp5"),
                      tail = "normal", slopes = "balanced", covariates = "none") {
  transform(expand.grid(id = ids, p_restudy = p, effect = effect, tail = tail,
                        slopes = slopes, covariates = covariates, stringsAsFactors = FALSE),
            phase = phase)
}
tails <- c("learner_high", "learner_low", "word_high", "word_low")
dominant <- c("learner_dominant", "word_dominant")
cov_profiles <- c("prognostic", "audio", "audio_strong", "method_wm", "method_pos", "joint")
cases <- rbind(
  scenarios("global", layouts$id, c(.10, .42, .90)),
  scenarios("global_logit", layouts$id, c(.10, .90), "fixed_logit"),
  scenarios("local", local, tail = tails),
  scenarios("variance_control", core, tail = c("learner_gaussian", "word_gaussian")),
  scenarios("mixture_strength", paired, tail = c("learner_high_mild", "word_low_mild")),
  scenarios("slope_sd", core, c(.10, .90), slopes = dominant),
  scenarios("slope_sd_local", paired, tail = c("learner_high", "word_low"), slopes = dominant),
  scenarios("zero_slopes", paired, c(.10, .90), slopes = "zero"),
  scenarios("correlation", paired, c(.10, .90),
            slopes = c("learner_rho_neg", "learner_rho_pos", "word_rho_neg", "word_rho_pos")),
  scenarios("covariates", core, covariates = cov_profiles),
  scenarios("covariates_local", paired, tail = "learner_high", covariates = "audio"),
  scenarios("covariates_local", paired, tail = "word_low", covariates = "method_wm")
)
cases <- cbind(case = seq_len(nrow(cases)), cases,
               layouts[match(cases$id, layouts$id), c("n_subj", "n_item", "m")])
cases$n_responses <- with(cases, n_subj * m)
cases$word_per_arm <- with(cases, n_responses / (2L * n_item))
cases$master_seed <- seed
cases$n_rep <- n_rep
cases$fit_int_only <- with(cases, id %in% paired & p_restudy != .42 &
                            phase %in% c("global", "zero_slopes"))
stopifnot(nrow(layouts) == 14L, nrow(cases) == 304L,
          !anyDuplicated(cases[c("id", "p_restudy", "effect", "tail", "slopes", "covariates")]))

slope_parameters <- rbind(
  balanced = c(.25, .25, 0, 0), learner_dominant = c(.40, .15, 0, 0),
  word_dominant = c(.15, .40, 0, 0), zero = c(0, 0, 0, 0),
  learner_rho_neg = c(.25, .25, -.5, 0), learner_rho_pos = c(.25, .25, .5, 0),
  word_rho_neg = c(.25, .25, 0, -.5), word_rho_pos = c(.25, .25, 0, .5)
)
colnames(slope_parameters) <- c("slope_subj", "slope_item", "rho_subj", "rho_item")
cov_parameters <- rbind(
  none = c(0, 0, 0, 0, 0, 0, 0, 0),
  prognostic = c(.45, .35, .20, 0, 0, 0, 0, 0),
  audio = c(.45, .35, .20, .15, -.30, .35, 0, 0),
  audio_strong = c(.45, .35, .20, .15, -.30, .70, 0, 0),
  method_wm = c(.45, .35, .20, .15, -.30, .35, .35, 0),
  method_pos = c(.45, .35, .20, .15, -.30, .35, 0, .35),
  joint = c(.45, .35, .20, .15, -.30, .70, .70, -.35)
)
colnames(cov_parameters) <- c("wm", "concrete", "pos", "discrimination", "contrast_word",
                              "phonological_match", "method_wm", "method_pos")
cases <- cbind(cases, slope_parameters[cases$slopes, ], cov_parameters[cases$covariates, ])
rownames(cases) <- NULL
cases$tail_share <- ifelse(grepl("high|low", cases$tail), .10, 0)
cases$tail_shift <- ifelse(grepl("high|low", cases$tail),
                           ifelse(grepl("high", cases$tail), 1, -1) *
                             ifelse(grepl("mild", cases$tail), 1.5, 3), 0)
cases$sd_subj <- sqrt(.8^2 - cases$wm^2 + .81 * (cases$tail == "learner_gaussian"))
cases$sd_item <- sqrt(.6^2 - cases$concrete^2 + .81 * (cases$tail == "word_gaussian"))

# The reference schedule and fitted random-effects structure are shared by the DGPs.
design <- design_crossed(list(method = c("restudy", "retrieval")), allocation = "latin_square")
dgp <- calibrate_dgp(design, c(restudy = .42, retrieval = .47), subj_sd = c(.8, .25),
                     item_sd = c(.6, .25), corr_structure = 0)
spec <- specify_target(dgp, "method", scale = "pp", engine = "bayes_laplace",
                       analysis_re = "maximal", prob_se_method = "joint_delta",
                       control = target_control(on_singular = "flag_only"))
make_schedule <- function(d) {
  if (d$m == d$n_item) {
    schedule <- generate_design_schedule(design, d$n_subj, d$n_item)
  } else {
    n_forms <- d$n_item %/% 4L
    stopifnot(d$n_item %% 4L == 0L, d$m %% 4L == 0L, d$n_subj %% (2L * n_forms) == 0L)
    trials <- do.call(rbind, lapply(seq_len(n_forms), function(f) {
      blocks <- (f - 1L + seq_len(d$m %/% 4L) - 1L) %% n_forms
      words <- as.vector(outer(seq_len(4L), 4L * blocks, "+"))
      size <- d$n_subj %/% n_forms
      g <- expand.grid(subj_id = seq_len(size) + (f - 1L) * size, item_id = words)
      g$cell_id <- as.integer((g$subj_id + g$item_id) %% 2L + 1L)
      g
    }))
    schedule <- list(trials = trials, n_subj = d$n_subj, n_item = d$n_item)
  }
  tr <- schedule$trials
  stopifnot(!anyDuplicated(tr[c("subj_id", "item_id")]), nrow(tr) == d$n_subj * d$m,
            all(table(tr$subj_id, tr$cell_id) == d$m / 2L),
            all(table(tr$item_id, tr$cell_id) == d$n_subj * d$m / (2L * d$n_item)),
            length(unique(glmmplanr:::item_components(tr, d$n_item))) == 1L)
  schedule
}
schedules <- setNames(lapply(seq_len(nrow(layouts)), function(i) make_schedule(layouts[i, ])),
                      layouts$id)

# The independent normal covariates can be integrated as one normal variable per arm.
discrete <- expand.grid(pos = c(-1, 1), discrimination = c(-1, 1), contrast_word = 0:1)
normal_part <- function(cfg, x) {
  s0 <- cfg$sd_subj
  i0 <- cfg$sd_item
  c(s0^2 + (x * cfg$slope_subj)^2 + 2 * x * cfg$rho_subj * s0 * cfg$slope_subj,
    i0^2 + (x * cfg$slope_item)^2 + 2 * x * cfg$rho_item * i0 * cfg$slope_item)
}
calibrate <- function(cfg) {
  arm_probability <- function(mu, x, method = "integrate") {
    sd <- sqrt(sum(normal_part(cfg, x)) + (cfg$wm + x * cfg$method_wm)^2 + cfg$concrete^2)
    shifts <- with(discrete, (cfg$pos + x * cfg$method_pos) * pos +
                     cfg$discrimination * discrimination + cfg$contrast_word * contrast_word +
                     cfg$phonological_match * discrimination * contrast_word)
    shifts <- outer(shifts, cfg$tail_shift * (0:1 - cfg$tail_share), "+")
    weights <- rep(c(1 - cfg$tail_share, cfg$tail_share) / nrow(discrete), each = nrow(discrete))
    use <- weights > 0
    shifts <- as.vector(shifts)[use]
    weights <- weights[use]
    if (method == "integrate") {
      return(sum(weights * vapply(shifts, function(s) {
        integrate(function(z) plogis(mu + s + sd * z) * dnorm(z), -Inf, Inf,
                  rel.tol = 1e-10)$value
      }, numeric(1))))
    }
    q <- glmmplanr:::gh_quad_1d(41L)
    sum(weights * as.numeric(plogis(outer(mu + shifts, sd * q$nodes, "+")) %*% q$weights))
  }
  from_p <- function(p, x) {
    uniroot(function(mu) arm_probability(mu, x) - p, c(-15, 15), tol = 1e-10)$root
  }
  mu0 <- from_p(cfg$p_restudy, -.5)
  mu1 <- switch(cfg$effect, null = from_p(cfg$p_restudy, .5),
                pp5 = from_p(cfg$p_restudy + .05, .5), fixed_logit = mu0 + dgp$beta[2L])
  p0 <- arm_probability(mu0, -.5, "quadrature")
  p1 <- arm_probability(mu1, .5, "quadrature")
  stopifnot(abs(p0 - cfg$p_restudy) < 1e-7,
            abs(p1 - arm_probability(mu1, .5)) < 1e-7)
  if (cfg$effect != "fixed_logit") {
    stopifnot(abs(p1 - p0 - ifelse(cfg$effect == "pp5", .05, 0)) < 1e-7)
  }
  c(mu_restudy = mu0, mu_retrieval = mu1, p_retrieval = p1, true_pp = 100 * (p1 - p0),
    true_logit = mu1 - mu0)
}
cases <- cbind(cases, t(vapply(seq_len(nrow(cases)), function(i) calibrate(cases[i, ]), numeric(5))))

features <- function(wm, concrete, pos, discrimination, contrast_word, x) {
  cbind(wm = wm, concrete = concrete, pos = pos, discrimination = discrimination,
        contrast_word = contrast_word, phonological_match = discrimination * contrast_word,
        method_wm = x * wm, method_pos = x * pos)
}

# Fixed-effect columns used for fitting and target standardization must coincide.
standardize <- function(columns, order = 9L) {
  q <- glmmplanr:::gh_quad_1d(order)
  grid <- expand.grid(wm_node = seq_along(q$nodes), concrete_node = seq_along(q$nodes),
                      discrete_row = if ("discrimination" %in% columns) seq_len(nrow(discrete))
                        else which(discrete$discrimination == -1 & discrete$contrast_word == 0))
  weights <- q$weights[grid$wm_node] * q$weights[grid$concrete_node] /
    length(unique(grid$discrete_row))
  id <- rep(1:2, each = nrow(grid))
  X <- design$X_cells[id, , drop = FALSE]
  groups <- discrete[rep(grid$discrete_row, 2L), ]
  F <- features(rep(q$nodes[grid$wm_node], 2L), rep(q$nodes[grid$concrete_node], 2L),
                groups$pos, groups$discrimination, groups$contrast_word, X[, 2L])
  s <- spec
  s$dgp$design$X_cells <- cbind(X, F[, columns, drop = FALSE])
  s$dgp$design$Z_subj_cells <- design$Z_subj_cells[id, , drop = FALSE]
  s$dgp$design$Z_item_cells <- design$Z_item_cells[id, , drop = FALSE]
  s$dgp$design$n_cells <- length(id)
  s$contrast_weights <- spec$contrast_weights[id] * rep(weights, 2L)
  stopifnot(abs(sum(weights) - 1) < 1e-12)
  s
}

max_subj <- max(layouts$n_subj)
max_item <- max(layouts$n_item)
exact_ci <- function(k, n) if (n > 0L) binom.test(k, n)$conf.int[1:2] else rep(NA_real_, 2L)
checkpoint <- file.path(out, "checkpoints")
dir.create(checkpoint, showWarnings = FALSE)
manifest <- file.path(checkpoint, "scenarios.rds")
if (file.exists(manifest)) stopifnot(identical(readRDS(manifest), cases)) else saveRDS(cases, manifest)
prefix <- file.path(out, paste0("l2_vocabulary_extremes_", stage))
write.csv(cases, paste0(prefix, "_scenarios.csv"), row.names = FALSE)

run_case <- function(k) {
  cfg <- cases[k, ]
  path <- file.path(checkpoint, sprintf("case_%03d.rds", k))
  if (file.exists(path)) return(path)
  cache <- glmmplanr:::prepare_cpp_schedule_cache(spec, schedules[[cfg$id]])
  s <- cache$sid
  i <- cache$iid
  cell <- cache$cell_id
  x <- cache$X_full[, 2L]
  pairs <- s + max_subj * (i - 1L)
  active <- names(cov_parameters[cfg$covariates, ])[cov_parameters[cfg$covariates, ] != 0]
  main_cols <- setdiff(active, c("method_wm", "method_pos"))
  moderator <- length(main_cols) < length(active)
  fit_specs <- list(unadjusted = spec)
  if (length(main_cols)) fit_specs$adjusted <- standardize(main_cols)
  if (moderator) fit_specs$moderation <- standardize(active)
  if (cfg$fit_int_only) fit_specs$intercept_only <- spec
  if (length(active) && cfg$tail_share == 0) {
    full <- fit_specs[[if (moderator) "moderation" else "adjusted"]]
    beta <- c(mean(c(cfg$mu_restudy, cfg$mu_retrieval)), cfg$true_logit,
              cov_parameters[cfg$covariates, active])
    p <- glmmplanr:::estimate_prob_effect_delta(
      beta, matrix(0, length(beta), length(beta)), diag(c(cfg$sd_subj, cfg$slope_subj)^2),
      diag(c(cfg$sd_item, cfg$slope_item)^2), NULL, full$dgp$design, 0,
      full$contrast_weights, n_quad = dgp$n_quad
    )$est
    stopifnot(abs(100 * p - cfg$true_pp) < 1e-5)
  }
  mu <- c(cfg$mu_restudy, cfg$mu_retrieval)[cell]
  if (length(active)) coeff <- cov_parameters[cfg$covariates, active]
  ns <- cfg$n_subj
  ni <- cfg$n_item
  item_exp <- 2L * cfg$word_per_arm
  tallies <- function(y, ids, units, exposures) {
    successes <- tabulate(ids[y == 1L], units)
    c(all_wrong = sum(successes == 0L), all_correct = sum(successes == exposures))
  }
  result <- lapply(seq_len(n_rep), function(r) {
    draw <- glmmplanr:::with_rng_stream(streams[[r]], {
      p <- matrix(rnorm(4L * max_subj), max_subj, 4L)
      w <- matrix(rnorm(5L * max_item), max_item, 5L)
      up <- runif(max_subj)
      uw <- runif(max_item)
      u <- runif(max_subj * max_item)
      eta <- mu + cfg$sd_subj * p[s, 1L] + cfg$sd_item * w[i, 1L] + x * (
        cfg$slope_subj * (cfg$rho_subj * p[s, 1L] + sqrt(1 - cfg$rho_subj^2) * p[s, 2L]) +
          cfg$slope_item * (cfg$rho_item * w[i, 1L] + sqrt(1 - cfg$rho_item^2) * w[i, 2L]))
      if (cfg$tail_share > 0) {
        v <- if (startsWith(cfg$tail, "learner")) up[s] else uw[i]
        eta <- eta + cfg$tail_shift * ((v < cfg$tail_share) - cfg$tail_share)
      }
      F <- if (length(active)) {
        features(p[s, 3L], w[i, 3L], ifelse(w[i, 4L] >= 0, 1, -1),
                 ifelse(p[s, 4L] >= 0, 1, -1), as.integer(w[i, 5L] >= 0), x)
      } else NULL
      if (length(active)) eta <- eta + as.numeric(F[, active, drop = FALSE] %*% coeff)
      list(y = as.integer(u[pairs] < plogis(eta)), features = F)
    })
    y <- draw$y
    s_count <- tallies(y, s, ns, cfg$m)
    i_count <- tallies(y, i, ni, item_exp)
    arm_count <- tabulate(cell[y == 1L], 2L)
    separation <- any(arm_count == 0L | arm_count == length(y) / 2L)
    empty_strata <- constant_strata <- NA_integer_
    if (length(active)) {
      f <- draw$features
      has_audio <- "discrimination" %in% active
      strata <- if (has_audio) {
        cell + 2L * (f[, "pos"] > 0) + 4L * (f[, "discrimination"] > 0) +
          8L * f[, "contrast_word"]
      } else cell + 2L * (f[, "pos"] > 0)
      n_strata <- if (has_audio) 16L else 4L
      n <- tabulate(strata, n_strata)
      correct <- tabulate(strata[y == 1L], n_strata)
      empty_strata <- sum(n == 0L)
      constant_strata <- sum(n > 0L & (correct == 0L | correct == n))
    }
    fits <- do.call(rbind, lapply(names(fit_specs), function(model) {
      fit_cache <- cache
      if (model %in% c("adjusted", "moderation")) {
        cols <- if (model == "adjusted") main_cols else active
        fit_cache$X_full <- cbind(cache$X_full, draw$features[, cols, drop = FALSE])
        fit_cache$ws_env <- NULL
      }
      f <- glmmplanr:::fit_via_cpp_cached(
        list(y = y, keep = NULL), fit_cache, fit_specs[[model]], regularized = TRUE,
        re_structure = if (model == "intercept_only") "intercept_only" else "maximal",
        fast_stage = FALSE
      )
      valid <- glmmplanr:::valid_target_fit(f$est, f$se, f$converged, f$status,
                                            f$prob_est, f$prob_se, "pp")
      est <- if (valid) 100 * f$prob_est else NA_real_
      se <- if (valid) 100 * f$prob_se else NA_real_
      contrast_coef <- function(term) {
        if (valid && model == "moderation" && term %in% active) {
          h <- 2L + match(term, active)
          return(c(f$beta_full[h], sqrt(f$vcov_full[h, h])))
        }
        c(NA_real_, NA_real_)
      }
      wm <- contrast_coef("method_wm")
      pos <- contrast_coef("method_pos")
      data.frame(model, valid, estimate_pp = est, se_pp = se,
                 covered = if (valid) abs(est - cfg$true_pp) <= z * se else NA,
                 detected = valid && abs(est) > z * se,
                 wm_est = wm[1L], wm_se = wm[2L], pos_est = pos[1L], pos_se = pos[2L],
                 status = f$status, diagnostic = if (valid) "" else f$error_reason)
    }))
    cbind(case = k, replication = r, p_restudy_obs = 2 * arm_count[1L] / length(y),
          p_retrieval_obs = 2 * arm_count[2L] / length(y),
          learner_wrong = unname(s_count[1L]), learner_correct = unname(s_count[2L]),
          word_wrong = unname(i_count[1L]), word_correct = unname(i_count[2L]),
          method_separation = separation, empty_strata, constant_strata, fits)
  })
  rows <- do.call(rbind, result)
  stopifnot(nrow(rows) == n_rep * length(fit_specs),
            !anyDuplicated(rows[c("replication", "model")]))
  tmp <- paste0(path, ".new")
  saveRDS(rows, tmp)
  stopifnot(file.rename(tmp, path))
  cat(sprintf("%s %3d/%d %s %s %s: %d/%d valid\n", stage, k, nrow(cases),
              cfg$id, cfg$phase, cfg$effect, sum(rows$valid), nrow(rows)))
  flush.console()
  path
}

summarize <- function(rows) {
  cfg <- cases[rows$case[1L], ]
  v <- rows$valid
  nv <- sum(v)
  n <- nrow(rows)
  cov <- sum(rows$covered, na.rm = TRUE)
  det <- sum(rows$detected)
  cov_ci <- exact_ci(cov, nv)
  det_ci <- exact_ci(det, n)
  coverage_gate <- if (nv / n < .90 || cov_ci[2L] < .90) "fail" else
    if (cov_ci[1L] >= .90) "pass" else "indeterminate"
  rate_gate <- if (cfg$effect == "null") {
    if (det_ci[2L] <= .075) "pass" else if (det_ci[1L] > .075) "fail" else "indeterminate"
  } else if (det_ci[1L] >= .80) "pass" else if (det_ci[2L] < .80) "fail" else "indeterminate"
  interactions <- function(est, se, truth) {
    active <- is.finite(est) & is.finite(se)
    if (!any(active)) return(rep(NA_real_, 4L))
    k <- sum(abs(est[active] - truth) <= z * se[active])
    c(mean(est[active]), k, exact_ci(k, sum(active)))
  }
  wm <- interactions(rows$wm_est, rows$wm_se, cfg$method_wm)
  pos <- interactions(rows$pos_est, rows$pos_se, cfg$method_pos)
  error <- rows$estimate_pp[v] - cfg$true_pp
  width <- 2 * z * rows$se_pp[v]
  cbind(cfg, model = rows$model[1L], attempted = n, valid = nv, covered = cov, detected = det,
        coverage_lo = cov_ci[1L], coverage_hi = cov_ci[2L], detection_lo = det_ci[1L],
        detection_hi = det_ci[2L], coverage_gate, rate_gate,
        bias_pp = mean(error), bias_mcse_pp = sd(error) / sqrt(nv), sd_pp = sd(error),
        mean_width_pp = mean(width), p80_width_pp = unname(quantile(width, .8)),
        mean_learner_wrong = mean(rows$learner_wrong), mean_learner_correct = mean(rows$learner_correct),
        mean_word_wrong = mean(rows$word_wrong), mean_word_correct = mean(rows$word_correct),
        method_separated = sum(rows$method_separation),
        trials_empty_strata = sum(rows$empty_strata > 0L, na.rm = TRUE),
        trials_constant_strata = sum(rows$constant_strata > 0L, na.rm = TRUE),
        mean_wm_est = wm[1L], wm_covered = wm[2L], wm_coverage_lo = wm[3L], wm_coverage_hi = wm[4L],
        mean_pos_est = pos[1L], pos_covered = pos[2L], pos_coverage_lo = pos[3L], pos_coverage_hi = pos[4L])
}

cat("Package:", normalizePath(find.package("glmmplanr")), "\n")
cat(sprintf("Stage %s: seed %d, %d scenarios, %d replications per scenario\n",
            stage, seed, nrow(cases), n_rep))
started <- proc.time()[["elapsed"]]
paths <- parallel::mclapply(seq_len(nrow(cases)), run_case, mc.cores = workers,
                            mc.preschedule = FALSE, mc.set.seed = FALSE)
stopifnot(all(vapply(paths, function(p) is.character(p) && length(p) == 1L && file.exists(p), logical(1))))
summary <- vector("list", length(paths))
connection <- gzfile(paste0(prefix, "_fits.csv.gz"), "wt")
for (k in seq_along(paths)) {
  rows <- readRDS(paths[[k]])
  stopifnot(identical(unique(rows$case), k), identical(unique(rows$replication), seq_len(n_rep)))
  summary[[k]] <- do.call(rbind, lapply(split(rows, rows$model), summarize))
  write.table(rows, connection, sep = ",", row.names = FALSE, col.names = k == 1L)
}
close(connection)
summary <- do.call(rbind, summary)
write.csv(summary, paste0(prefix, "_summary.csv"), row.names = FALSE)
cat(sprintf("Completed %d datasets and %d fits (%d valid) in %.1f seconds.\n",
            nrow(cases) * n_rep, sum(summary$attempted), sum(summary$valid),
            proc.time()[["elapsed"]] - started))
