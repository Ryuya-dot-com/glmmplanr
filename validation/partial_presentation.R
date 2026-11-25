# Run from the repository root against the updated installed package.
library(glmmplanr)
library(lme4)

n_subj <- 64L
n_item <- 20L
n_rep <- 32L
n_boot <- 29L
designs <- list(full = design_crossed(list(group = c("control", "treatment")),
                                       allocation = "between_subj"),
                partial = design_crossed(list(group = c("control", "treatment")),
                                          allocation = "between_subj", items_per_subject = 5,
                                          item_seed = 417L))
schedules <- lapply(designs, generate_design_schedule, n_subj = n_subj, n_item = n_item)
stopifnot(all(table(schedules$partial$trials$subj) == 5L),
          all(table(schedules$partial$trials$group) == 160L),
          !anyDuplicated(schedules$partial$trials[c("subj_id", "item_id")]))
full <- schedules$full$trials
partial <- schedules$partial$trials
# Match the identical participant-item rows in both schedules.
subset_rows <- match((partial$item_id - 1L)*n_subj + partial$subj_id,
                     (full$item_id - 1L)*n_subj + full$subj_id)
stopifnot(!anyNA(subset_rows), identical(full$group[subset_rows], partial$group))

# External reference: marginalize over both fitted random intercepts by numeric
# integration. Bootstrap fresh random intercepts and binary responses on the
# observed design; each successful refit supplies one marginal contrast.
reference_fit <- function(d, refits = n_boot) {
  singular <- NA
  diagnostic <- ""
  ci <- tryCatch({
    model <- suppressMessages(glmer(y ~ group + (1 | subj) + (1 | item), data = d,
      family = binomial(), control = glmerControl(optimizer = "bobyqa")))
    singular <- isSingular(model)
    effect <- function(f) {
      conv <- f@optinfo$conv
      msg <- conv$lme4$messages
      if (any(conv$opt != 0L) || any(!grepl("boundary (singular) fit", msg, fixed = TRUE))) return(NA_real_)
      sd_re <- sqrt(sum(vapply(VarCorr(f), function(v) v[1,1], numeric(1))))
      mu <- drop(getME(f, "X")[match(c("control", "treatment"), d$group), , drop = FALSE] %*% fixef(f))
      integrate(function(z) (plogis(mu[2] + sd_re*z) - plogis(mu[1] + sd_re*z))*dnorm(z),
                -9, 9, abs.tol = 1e-9)$value
    }
    est <- effect(model)
    if (!is.finite(est)) stop("Base fit did not converge")
    boot <- bootMer(model, effect, nsim = refits, type = "parametric", use.u = FALSE)
    usable <- is.finite(boot$t)
    if (sum(usable) < ceiling(.9*refits)) stop("Too few valid bootstrap refits")
    c(estimate = est, se = sd(boot$t[usable]))
  }, error = function(e) {
    diagnostic <<- conditionMessage(e)
    c(estimate = NA_real_, se = NA_real_)
  })
  list(estimate = unname(ci[1]), se = unname(ci[2]), singular = singular, diagnostic = diagnostic)
}

specs <- lapply(c(null = .60, planned = .55), function(p) lapply(designs, function(des) {
  dgp <- calibrate_dgp(des, c(.60, p), chance = 0, subj_sd = .5, item_sd = .3,
                       slope_sd_ratio = 0, corr_structure = 0)
  specify_target(dgp, "group", scale = "pp", engine = "bayes_laplace",
                 analysis_re = "intercept_only",
                 control = target_control(prob_estimand = "marginal", on_singular = "flag_only", prob_boot = n_boot))
}))

run_rep <- function(rep_id, conditions = names(specs), layouts = names(designs), refits = n_boot) {
  rows <- list()
  for (condition in conditions) {
    spec <- specs[[condition]]
    target <- if (condition == "null") 0 else -5
    stopifnot(abs(spec$full$true_effect_target - target) < 1e-6)
    seed <- 871500L + rep_id
    set.seed(seed, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion", sample.kind = "Rejection")
    data <- simulate_dgp_data(spec$full$dgp, schedules$full)
    for (layout in layouts) {
      d <- if (layout == "partial") data[subset_rows, ] else data
      set.seed(872500L + rep_id)
      sp <- spec[[layout]]
      sp$prob_boot <- refits
      fit <- glmmplanr:::fit_single_dataset(d, sp)
      native <- list(estimate = fit$prob_est, se = fit$prob_se, singular = fit$singular,
                     diagnostic = fit$error_reason)
      if (!isTRUE(fit$converged)) native$se <- NA_real_
      set.seed(873500L + rep_id)
      fits <- list(bayes_laplace = native, glmer = reference_fit(d, refits))
      for (engine in names(fits)) {
        f <- fits[[engine]]
        est <- 100*f$estimate
        se <- 100*f$se
        valid <- is.finite(est) && is.finite(se) && se > 0
        half <- if (valid) qnorm(.975)*se else NA_real_
        rows[[length(rows)+1L]] <- data.frame(replication = rep_id, seed, condition, target_pp = target,
          layout, engine, bootstrap_refits = refits, trials = nrow(d), items_min = min(table(d$item)), items_max = max(table(d$item)),
          valid, estimate_pp = est, width_pp = 2*half, detected = valid && abs(est) > half,
          covered = if (valid) abs(est - target) <= half else NA, singular = f$singular,
          diagnostic = f$diagnostic)
      }
    }
  }
  do.call(rbind, rows)
}

summarize_presentation <- function(raw) {
  keys <- c("condition","layout","engine","bootstrap_refits")
  groups <- split(raw, raw[keys], drop = TRUE)
  do.call(rbind, lapply(groups, function(d) {
    ok <- d$valid
    n <- nrow(d)
    valid <- sum(ok)
    hits <- sum(d$detected)
    cover <- sum(d$covered[ok])
    di <- binom.test(hits, n)$conf.int
    ci <- if (valid) binom.test(cover, valid)$conf.int else c(NA,NA)
    data.frame(d[1, c(keys, "target_pp","trials","items_min","items_max")], attempted = n,
      valid, detected = hits, detection = hits/n, detection_lo = di[1], detection_hi = di[2],
      covered = cover, coverage = cover/valid, coverage_lo = ci[1], coverage_hi = ci[2],
      estimate_pp = mean(d$estimate_pp[ok]), width_pp = mean(d$width_pp[ok]),
      singular_valid = if (all(is.na(d$singular[ok]))) NA_integer_ else sum(d$singular[ok], na.rm = TRUE),
      R_version = as.character(getRversion()), lme4_version = as.character(packageVersion("lme4")))
  }))
}

if (sys.nframe() == 0L) {
  run_both <- function(i) rbind(run_rep(i), run_rep(i, "planned", "full", 99L))
  runs <- parallel::mclapply(seq_len(n_rep), run_both,
    mc.cores = if (.Platform$OS.type == "windows") 1L else 2L, mc.set.seed = FALSE)
  stopifnot(all(vapply(runs, is.data.frame, logical(1))))
  raw <- do.call(rbind, runs)
  summary <- summarize_presentation(raw)
  stopifnot(nrow(raw) == 10L*n_rep, nrow(summary) == 10L,
            !anyDuplicated(raw[c("replication","condition","layout","engine","bootstrap_refits")]))
  write.csv(raw, "validation/partial_presentation_fits.csv", row.names = FALSE)
  write.csv(summary, "validation/partial_presentation_summary.csv", row.names = FALSE)
  print(summary[c("condition", "layout", "engine", "bootstrap_refits", "valid", "detected",
                  "covered", "estimate_pp", "width_pp")], row.names = FALSE)
}
