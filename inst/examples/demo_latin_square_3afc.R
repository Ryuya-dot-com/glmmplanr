# End-to-End Demo: 2x2 Latin-Square 3AFC GLMM Sample Size Plan (Pareto & pp Scale)
library(glmmplanr)

# 1. 2x2 Latin-square design (within-subject x within-item, 4-list balance, 5% trial attrition)
des <- design_crossed(
  factors = list(
    exposure = c("incidental", "intentional"),
    spacing  = c("massed", "distributed")
  ),
  allocation = "latin_square",
  contrast   = "sum",
  attrition  = 0.05
)

# 2. Construct cell probabilities from percentage-point effects & invert 3AFC (1/3) DGP
c_probs <- cell_probs_from_diffs(
  des,
  baseline_prob = 0.65,
  effects_pp    = c(exposure = 8, spacing = 6, "exposure:spacing" = 4)
)

dgp <- des |>
  calibrate_dgp(
    cell_probs      = c_probs,
    chance          = 1 / 3,
    subj_prob_range = c(0.42, 0.85),
    item_prob_range = c(0.45, 0.86),
    slope_sd_ratio  = 0.30,
    corr_structure  = 0.20
  )

check_dgp(dgp, plot = FALSE)$marginal_effects

# 3. Design directly in percentage points: Power >= 80% & 95% CI width <= 12 pp
spec <- dgp |>
  specify_target(
    effect         = "exposure",
    scale          = "pp",
    prob_estimand  = "marginal",
    criterion      = "both",
    power_min      = 0.80,
    ci_width       = 12,
    failure_policy = "effective",
    engine         = "bayes_laplace"
  )

# 4. Multi-scenario robust planning (5 / 8 / 10 pp) + Pareto frontier exploration
scens <- scenario_dgps(
  dgp,
  effect          = "exposure",
  effect_diffs_pp = c(diff_5pp = 5, diff_8pp = 8, diff_10pp = 10)
)

plan_res <- spec |>
  optimize_sample_size(
    n_subj_range      = c(24, 120),
    n_item_range      = c(16, 64),
    prefer            = "balanced",
    scenarios         = scens,
    standard_scenario = "diff_8pp",
    method            = "surrogate_gp",
    n_eval_init       = 16,
    nsim_final        = 1000,
    nsim_topup        = 500,
    workers           = 8,
    seed              = 20261119
  )

# 5. Summary, Pareto switching, gt table, 2D/3D plots, and the sample size report
summary(plan_res)
head(as.data.frame(plan_res)[c("n_subj", "n_item", "power", "decision_status")])
switch_preference(plan_res, prefer = "fewer_items")
as_gt(plan_res)
plot(plan_res, type = "contour", metric = "prob_ci_width_pp")
plot(plan_res, type = "pareto")
plot(plan_res, type = "surface3d", metric = "power")
report_plan(plan_res)
