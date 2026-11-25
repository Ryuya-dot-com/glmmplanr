# Pre-specified confirmation of conditional exposure precision and zero exclusion in 3AFC.
# Candidates: (72 x 32) and (48 x 48) [both 2,304 scheduled trials], 192 fresh datasets each,
# Disjoint seed series 7510000 + 10000 * effect_id + 1000 * cand_id + trial.
library(glmmplanr)
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
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !dir.exists(args[[1L]])) {
  stop("Supply the path to an existing output directory.")
}
out_dir <- args[[1L]]

dgp_planned <- three_choice_reference_dgp()
des <- dgp_planned$design
# Historical CSV label "null" denotes equal marginal accuracies within spacing.
# The calibrated conditional exposure coefficient is about +0.004097 logit.
dgp_null <- calibrate_dgp(
  des,
  cell_probs = c(
    incidental_massed       = 0.59,
    incidental_distributed  = 0.73,
    intentional_massed      = 0.59,
    intentional_distributed = 0.73
  ),
  chance = 1 / 3,
  subj_prob_range = c(0.42, 0.84),
  item_prob_range = c(0.45, 0.86),
  slope_sd_ratio = 0.30,
  corr_structure = 0.15,
  control = calibrate_control(seed = 1L)
)

dgps <- list(null = dgp_null, planned = dgp_planned)
specs <- lapply(dgps, function(g) {
  specify_target(
    g, effect = "exposure", scale = "logit", engine = "bayes_laplace",
    control = target_control(on_singular = "flag_only")
  )
})
# Strict conditional null (true conditional exposure coefficient = 0 logit)
specs$conditional_null <- glmmplanr:::null_target_spec(specs$planned)

candidates <- data.frame(cand_id = 1:2, persons = c(72L, 48L), items = c(32L, 48L))
n_rep <- 192L
z <- stats::qnorm(0.975)

all_rows <- list()
idx <- 0L
for (e_id in seq_along(specs)) {
  e_nm <- names(specs)[e_id]
  sp <- specs[[e_nm]]
  for (c_i in seq_len(nrow(candidates))) {
    ns <- candidates$persons[c_i]
    ni <- candidates$items[c_i]
    sched <- generate_design_schedule(des, n_subj = ns, n_item = ni)

    reps <- lapply(seq_len(n_rep), function(tr_i) {
      sd_val <- 7510000L + 10000L * e_id + 1000L * c_i + tr_i
      set.seed(sd_val, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion", sample.kind = "Rejection")
      dat <- simulate_dgp_data(sp$dgp, schedule = sched)
      f <- glmmplanr:::fit_single_dataset(dat, sp, engine_use = "bayes_laplace")
      valid <- isTRUE(f$converged) && is.finite(f$est) && is.finite(f$se) && f$se > 0
      data.frame(
        effect = e_nm,
        persons = ns,
        items = ni,
        trial = tr_i,
        seed = sd_val,
        truth = sp$true_beta,
        valid = valid,
        est = if (valid) f$est else NA_real_,
        se = if (valid) f$se else NA_real_,
        stringsAsFactors = FALSE
      )
    })
    idx <- idx + 1L
    all_rows[[idx]] <- do.call(rbind, reps)
  }
}

fits <- do.call(rbind, all_rows)
by_cell <- split(fits, list(fits$effect, fits$persons, fits$items), drop = TRUE)

judge <- function(ci, threshold) {
  if (anyNA(ci)) return("indeterminate")
  if (ci[1L] >= threshold) return("pass")
  if (ci[2L] < threshold) return("below_floor")
  "indeterminate"
}

summary_df <- do.call(rbind, lapply(by_cell, function(d) {
  n <- nrow(d)
  nv <- sum(d$valid)
  err <- d$est[d$valid] - d$truth[d$valid]
  se  <- d$se[d$valid]
  width <- 2 * z * se
  covered <- sum(abs(err) <= z * se)
  narrow  <- sum(width <= 1.0)
  detected <- sum(abs(d$est[d$valid]) > z * se)

  cov_ci <- as.numeric(stats::binom.test(covered, nv)$conf.int)
  wid_ci <- as.numeric(stats::binom.test(narrow, n)$conf.int)
  det_ci <- as.numeric(stats::binom.test(detected, n)$conf.int)

  cov_st <- judge(cov_ci, 0.90)
  wid_st <- judge(wid_ci, 0.80)
  joint_st <- if (nv / n < 0.90) {
    "insufficient_valid"
  } else if ("below_floor" %in% c(cov_st, wid_st)) {
    "below_floor"
  } else if (cov_st == "pass" && wid_st == "pass") {
    "pass"
  } else {
    "indeterminate"
  }

  data.frame(
    effect = d$effect[1L],
    persons = d$persons[1L],
    items = d$items[1L],
    attempted = n,
    valid = nv,
    covered = covered,
    coverage = covered / nv,
    coverage_lo = cov_ci[1L],
    coverage_hi = cov_ci[2L],
    coverage_status = cov_st,
    narrow = narrow,
    assurance = narrow / n,
    width_lo = wid_ci[1L],
    width_hi = wid_ci[2L],
    width_status = wid_st,
    detected = detected,
    detect_rate = detected / n,
    detect_lo = det_ci[1L],
    detect_hi = det_ci[2L],
    joint_status = joint_st,
    bias = mean(err),
    empirical_sd = stats::sd(err),
    mean_se = mean(se),
    mean_width = mean(width),
    stringsAsFactors = FALSE
  )
}))
summary_df <- summary_df[order(summary_df$effect, summary_df$persons), ]
rownames(summary_df) <- NULL

write.csv(fits, file.path(out_dir, "three_choice_confirmation_fits.csv"), row.names = FALSE)
write.csv(summary_df, file.path(out_dir, "three_choice_confirmation_summary.csv"), row.names = FALSE)
print(summary_df[, c("effect", "persons", "items", "valid", "covered", "coverage_lo", "narrow", "width_lo", "detected", "detect_rate", "joint_status")],
      row.names = FALSE, digits = 4)
