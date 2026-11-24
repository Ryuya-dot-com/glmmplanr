# Run with Rscript --vanilla validation/rasch_precision.R after installation.
library(glmmplanr)
dgp <- calibrate_dgp(design_crossed(list()), 0.65, subj_sd = 0.8, item_sd = 0.6)
designs <- data.frame(n_subj = c(144, 72, 144), n_item = c(12, 24, 24))
tolerance <- c(ability = 0.50, difficulty = 0.30, item_sd = 0.20)
plan <- plan_rasch_precision(dgp, designs, tolerance, nsim = 64, seed = 518307)
cols <- c("n_subj", "n_item", "target", "n_valid_fit", "n_singular", "n_met",
          "rmse", "mc_lo", "mc_hi", "decision")
print(as.data.frame(plan)[cols], row.names = FALSE)
print(plan$pareto[c("n_subj", "n_item")], row.names = FALSE)

# Use fresh simulated studies for a prespecified repeat of the largest design.
repeat_plan <- plan_rasch_precision(dgp, data.frame(n_subj = 144, n_item = 24),
                                    tolerance, nsim = 64, seed = 924031)
print(as.data.frame(repeat_plan)[cols], row.names = FALSE)

# Sensitivity check added after the baseline comparison: keep the marginal
# accuracy, person SD and absolute tolerances, and widen the item distribution.
wide_bank <- calibrate_dgp(dgp$design, 0.65, subj_sd = 0.8, item_sd = 1.0)
sensitivity <- plan_rasch_precision(wide_bank, data.frame(n_subj = 144, n_item = 24),
                                    tolerance, nsim = 64, seed = 714823)
print(as.data.frame(sensitivity)[cols], row.names = FALSE)
