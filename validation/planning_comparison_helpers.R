# Common targets and bookkeeping for the matched package comparisons.
library(glmmplanr)
RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
pool <- read.csv("validation/item_pool_planning_verified/pool.csv")
forms <- read.csv("validation/item_pool_planning_verified/forms.csv")
layouts <- data.frame(form = c("central30", "spread30", "spread60", "full150", "rotating150"),
                      n_item = c(30L, 30L, 60L, 150L, 150L), m = c(30L, 30L, 60L, 150L, 30L))
rasch_cases <- merge(layouts, expand.grid(person_sd = c(.4, .8), n_subj = c(120L, 240L, 480L, 960L, 1920L)))
rasch_cases$case <- seq_len(nrow(rasch_cases))
rasch_cases$family <- paste(rasch_cases$form, rasch_cases$person_sd, sep = "_")
bank_for <- function(case) {
  if (case$form == "rotating150") return(setNames(pool$difficulty, pool$item_id))
  x <- subset(forms, form == case$form)
  setNames(x$difficulty, x$item_id)
}
mc_summary <- function(success, valid = rep(TRUE, length(success))) {
  stopifnot(!anyNA(success), !anyNA(valid), length(success) == length(valid), all(!success | valid))
  n <- length(success); k <- sum(success)
  ci <- glmmplanr:::wilson_mc_ci(k, n)
  data.frame(nsim = n, n_valid = sum(valid), n_met = k, attainment = k / n,
    mc_lo = unname(ci[1]), mc_hi = unname(ci[2]),
    decision = if (ci[1] >= .8) "pass" else if (ci[2] < .8) "fail" else "indeterminate")
}
checkpoint_dir <- function(out, scripts, settings, packages) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  files <- c(scripts, "validation/planning_comparison_helpers.R",
    "validation/item_pool_planning_verified/pool.csv", "validation/item_pool_planning_verified/forms.csv",
    "src/glmm_laplace.cpp", "R/rasch_precision.R", "R/engines.R", "R/optimize.R")
  binary <- list.files(system.file("libs", package = "glmmplanr"), pattern = "\\.(so|dll)$",
                       recursive = TRUE, full.names = TRUE)
  manifest <- list(settings = settings, files = tools::md5sum(files),
    binary = unname(tools::md5sum(binary)), R = R.version.string,
    packages = vapply(packages, function(p) as.character(packageVersion(p)), ""))
  path <- file.path(out, "manifest.rds")
  if (file.exists(path)) {
    if (!identical(readRDS(path), manifest)) stop("Resume manifest differs: ", out)
  } else {
    if (length(list.files(out))) stop("Existing directory has no matching manifest: ", out)
    saveRDS(manifest, path)
    capture.output(sessionInfo(), file = file.path(out, "environment.txt"))
  }
  invisible(manifest)
}
atomic_rds <- function(x, path) {
  tmp <- paste0(path, ".tmp")
  saveRDS(x, tmp)
  stopifnot(file.rename(tmp, path))
}
rasch_runner <- function(case) {
  b <- bank_for(case)
  design <- design_crossed(list(), items_per_subject = case$m, item_seed = 37L)
  schedule <- generate_design_schedule(design, case$n_subj, length(b))
  function() {
    sim <- glmmplanr:::simulate_fixed_bank_trials(design, schedule, b, case$person_sd, lightweight = TRUE)
    fit <- glmmplanr:::eval_fixed_bank_replication(sim, b, case$person_sd, "fixed_item", "mml")
    m <- fit$metrics
    valid <- is.na(m$fit_error) && all(is.finite(fit$b_est))
    data.frame(valid = valid, success = valid && m$difficulty <= .35 && m$difficulty_max_abs_error <= .60)
  }
}
experiment_spec <- function(slopes) {
  des <- design_crossed(list(condition = c("control", "treatment")),
    allocation = "latin_square", items_per_subject = 30L, item_seed = 37L)
  varying <- slopes == "heterogeneous"
  dgp <- calibrate_dgp(des, c(control = .50, treatment = .58),
    subj_sd = c(.8, if (varying) .25 else 0), item_sd = c(.6, if (varying) .35 else 0),
    corr_structure = "zero")
  specify_target(dgp, "condition", scale = "pp", criterion = "joint", ci_width = 20,
    engine = "bayes_laplace", analysis_re = "match_dgp", prob_se_method = "joint_delta")
}
experiment_runner <- function(spec, n, j) {
  schedule <- generate_design_schedule(spec$dgp$design, n, j)
  cache <- glmmplanr:::prepare_cpp_schedule_cache(spec, schedule)
  function() {
    fit <- glmmplanr:::execute_single_replication(spec, schedule, cache,
      "bayes_laplace", spec$analysis_re, FALSE, FALSE)$primary
    valid <- isTRUE(glmmplanr:::valid_target_fit(fit$est, fit$se, fit$converged,
      fit$status, fit$prob_est, fit$prob_se, "pp"))
    detected <- valid && abs(fit$prob_est / fit$prob_se) > qnorm(.975)
    width_ok <- valid && 2 * qnorm(.975) * fit$prob_se * 100 <= 20
    data.frame(valid, success = detected && width_ok)
  }
}
