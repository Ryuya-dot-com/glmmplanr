# Run from the repository root: Rscript --vanilla validation/native_interactions.R
# Requires the installed version of glmmplanr that is being evaluated.
library(glmmplanr)

baseline <- design_crossed(
  list(A = c("incidental", "intentional"), B = c("massed", "distributed")),
  attrition = .05
)
baseline_dgp <- calibrate_dgp(baseline,
  cell_probs_from_diffs(baseline, .65, c(A = 8, B = 6, "A:B" = 4)), chance = 1/3,
  subj_prob_range = c(.42, .85), item_prob_range = c(.45, .86),
  slope_sd_ratio = .30, corr_structure = .20)

# One between-participant factor, one between-item factor, and one within both.
# Both interaction contrasts are evaluated in the same full analysis model.
mixed <- design_crossed(list(A = c("a0", "a1"), B = c("b0", "b1"), C = c("c0", "c1")),
                        between_subj = "A", between_item = "B", attrition = .05)
mixed_dgp <- calibrate_dgp(mixed,
  cell_probs_from_diffs(mixed, .65, c(A = 8, B = 6, C = 4, "A:B" = 8, "A:B:C" = 8)),
  chance = 1/3, subj_prob_range = c(.42, .85), item_prob_range = c(.45, .86),
  slope_sd_ratio = .30, corr_structure = .20)

cases <- data.frame(design = c("baseline", "baseline", "mixed", "mixed"),
                    effect = c("A", "A", "A:B", "A:B:C"), n_item = c(16, 32, 16, 16),
                    seed = c(81231L, 81232L, 81233L, 81233L))
results <- lapply(seq_len(nrow(cases)), function(i) {
  case <- cases[i, ]
  dgp <- if (case$design == "baseline") baseline_dgp else mixed_dgp
  spec <- specify_target(dgp, case$effect, scale = "pp", criterion = "joint",
                         power_min = .80, ci_width = 12, engine = "bayes_laplace", prob_boot = 19)
  stopifnot(abs(spec$true_effect_target - 8) < 1e-5)
  row <- validate_design_point(spec, n_subj = 32, n_item = case$n_item, nsim = 16, seed = case$seed)
  stopifnot(abs(row$true_target[row$condition == "null"]) < 1e-6)
  message(sprintf("Finished %s, %s, 32 x %d", case$design, case$effect, case$n_item))
  data.frame(design = case$design, effect = case$effect, n_subj = 32L, n_item = case$n_item,
             seed = case$seed, bootstrap_refits = spec$prob_boot, row,
             n_detected = round(row$nsim * row$detection),
             n_covered = round(row$n_valid * row$coverage),
             n_joint = round(row$nsim * row$joint_attainment),
             package_version = as.character(packageVersion("glmmplanr")), R_version = as.character(getRversion()))
})
results <- do.call(rbind, results)
write.csv(results, "validation/native_results.csv", row.names = FALSE)
print(results[c("design", "effect", "n_item", "condition", "n_valid", "n_covered", "n_detected", "n_joint")],
      row.names = FALSE)
