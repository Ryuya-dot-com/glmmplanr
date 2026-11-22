#!/usr/bin/env Rscript
# ==============================================================================
# Build Visual & Interactive Gallery Assets for `glmmplanr`
# Generates:
#   - man/figures/gallery_01_latin_square_3afc.png
#   - man/figures/gallery_02_incomplete_block_dgp.png
#   - man/figures/gallery_03_testlet_robustness.png
#   - man/figures/gallery_04_cost_reversal.png
#   - man/figures/gallery_05_rasch_calibration.png
#   - gallery_payload.json in tempdir() (or --json-out=PATH)
# ==============================================================================

suppressPackageStartupMessages({
  library(glmmplanr)
  library(ggplot2)
  library(grid)
})

dir.create("man/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("inst/gallery", recursive = TRUE, showWarnings = FALSE)

save_two_panel <- function(path, p_left, p_right, width_px = 1860, height_px = 700, res = 132) {
  pad_theme <- ggplot2::theme(
    plot.margin = ggplot2::margin(12, 22, 12, 18),
    plot.title = ggplot2::element_text(size = 11.5, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 9.5, color = "#334155")
  )
  grDevices::png(filename = path, width = width_px, height = height_px, res = res)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    layout = grid::grid.layout(nrow = 1L, ncol = 2L, widths = grid::unit(c(1.04, 0.96), "null"))
  ))
  print(p_left + pad_theme, vp = grid::viewport(layout.pos.row = 1L, layout.pos.col = 1L))
  print(p_right + pad_theme, vp = grid::viewport(layout.pos.row = 1L, layout.pos.col = 2L))
}

# Minimal base-R JSON serializer (avoids extra package requirements)
to_json_val <- function(x) {
  if (is.null(x)) return("null")
  if (is.data.frame(x)) {
    rows <- lapply(seq_len(nrow(x)), function(i) {
      vals <- vapply(names(x), function(nm) {
        paste0('"', nm, '":', to_json_val(x[[nm]][i]))
      }, character(1L))
      paste0("{", paste(vals, collapse = ","), "}")
    })
    return(paste0("[", paste(rows, collapse = ","), "]"))
  }
  if (is.list(x)) {
    nms <- names(x)
    if (is.null(nms)) {
      return(paste0("[", paste(vapply(x, to_json_val, character(1L)), collapse = ","), "]"))
    }
    kv <- vapply(seq_along(x), function(i) {
      paste0('"', nms[i], '":', to_json_val(x[[i]]))
    }, character(1L))
    return(paste0("{", paste(kv, collapse = ","), "}"))
  }
  if (is.logical(x)) {
    s <- ifelse(is.na(x), "null", ifelse(x, "true", "false"))
    return(if (length(x) == 1L) s else paste0("[", paste(s, collapse = ","), "]"))
  }
  if (is.numeric(x)) {
    s <- ifelse(!is.finite(x), "null", format(round(x, 5L), trim = TRUE, scientific = FALSE))
    return(if (length(x) == 1L) s else paste0("[", paste(s, collapse = ","), "]"))
  }
  esc <- function(v) {
    v <- gsub("\\\\", "\\\\\\\\", as.character(v))
    v <- gsub('"', '\\\\"', v)
    v <- gsub("\n", "\\\\n", v)
    v <- gsub("\r", "", v)
    paste0('"', v, '"')
  }
  s <- ifelse(is.na(x), "null", esc(x))
  if (length(x) == 1L) s else paste0("[", paste(s, collapse = ","), "]")
}

workers_n <- min(4L, max(1L, parallel::detectCores(logical = FALSE)))
cat(sprintf("=== Building glmmplanr Gallery (workers = %d) ===\n", workers_n))

# ------------------------------------------------------------------------------
# RECIPE 1: Psycholinguistics — 2x2 Latin-Square 3AFC Vocabulary Experiment
# ------------------------------------------------------------------------------
cat("[1/5] Recipe 1: Balanced 2x2 Latin-Square 3AFC Experiment...\n")
des1 <- design_crossed(
  factors = list(
    exposure = c("incidental", "intentional"),
    spacing  = c("massed", "distributed")
  ),
  allocation = "latin_square",
  contrast   = "sum",
  attrition  = 0.05
)
probs1 <- cell_probs_from_diffs(
  des1,
  baseline_prob = 0.64,
  effects_pp    = c(exposure = 9, spacing = 5, "exposure:spacing" = 3)
)
dgp1 <- calibrate_dgp(
  des1,
  cell_probs      = probs1,
  chance          = 1 / 3,
  subj_prob_range = c(0.45, 0.84),
  item_prob_range = c(0.46, 0.85),
  slope_sd_ratio  = 0.25,
  corr_structure  = 0.15
)
spec1 <- specify_target(
  dgp1,
  effect         = "exposure",
  scale          = "pp",
  prob_estimand  = "marginal",
  criterion      = "both",
  power_min      = 0.80,
  ci_width       = 12.0,
  prob_se_method = "joint_delta",
  engine         = "bayes_laplace"
)
res1 <- optimize_sample_size(
  spec1,
  n_subj_range = c(24, 80),
  n_item_range = c(16, 48),
  prefer       = "balanced",
  method       = "grid",
  nsim_init    = 60,
  nsim_final   = 180,
  nsim_topup   = 80,
  workers      = workers_n,
  seed         = 810101,
  verbose      = FALSE
)

p1_left <- plot(res1, type = "contour", metric = "power") +
  labs(
    title = "A. 2D Power & Precision Surface (Latin-Square 3AFC)",
    subtitle = sprintf(
      "Balanced winner: (%d subj, %d items) | Power = %.1f%% | Mean CI = %.1f pp",
      res1$optimal$n_subj, res1$optimal$n_item,
      res1$optimal$power * 100, res1$optimal$prob_ci_width_mean_pp
    )
  )
p1_right <- plot(res1, type = "pareto") +
  scale_x_continuous(expand = expansion(mult = c(0.10, 0.12))) +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.16))) +
  labs(
    title = "B. Non-Dominated Pareto Frontier",
    subtitle = "Cyan nodes satisfy both Power >= 80% and 95% CI Width <= 12 pp"
  )
save_two_panel("man/figures/gallery_01_latin_square_3afc.png", p1_left, p1_right)

# ------------------------------------------------------------------------------
# RECIPE 2: Behavioral Experiment — Between-Subject Incomplete-Block & Calibration
# ------------------------------------------------------------------------------
cat("[2/5] Recipe 2: Between-Subject Incomplete-Block & Non-Collapsibility...\n")
des2 <- design_crossed(
  factors             = list(arm = c("control", "intervention")),
  allocation          = "between_subj",
  between_subj        = "arm",
  contrast            = "treatment",
  items_per_subject   = 16L,
  attrition           = 0.05,
  participant_dropout = 0.04,
  item_seed           = 402L
)
dgp2 <- calibrate_dgp(
  des2,
  cell_probs     = c(control = 0.60, intervention = 0.72),
  chance         = 0,
  subj_sd        = 1.05,
  item_sd        = 0.80,
  slope_sd_ratio = 0,
  corr_structure = 0
)
chk2 <- check_dgp(dgp2, plot = FALSE)

# Compute what happens if a user naively feeds qlogis(c(0.60, 0.72)) into a simulator
dgp2_naive <- dgp2
dgp2_naive$beta <- dgp2$beta_naive_logit
dgp2_naive$mu_cells <- stats::setNames(drop(des2$X_cells %*% dgp2_naive$beta), des2$cells$cell_name)
naive_marg <- vapply(seq_len(des2$n_cells), function(c_i) {
  w <- numeric(des2$n_cells); w[c_i] <- 1.0
  true_prob_effect(dgp2_naive, w, "marginal")
}, numeric(1L))

scale_df2 <- data.frame(
  cell = rep(c("Control (60%)", "Intervention (72%)"), each = 3L),
  quantity = factor(
    rep(c("Target Marginal (quad_prob)", "Conditional at u = 0", "Naive qlogis() Marginal"), 2L),
    levels = c("Target Marginal (quad_prob)", "Conditional at u = 0", "Naive qlogis() Marginal")
  ),
  prob = c(
    chk2$cell_diagnostics$quad_prob[1], chk2$cell_diagnostics$conditional_prob[1], naive_marg[1],
    chk2$cell_diagnostics$quad_prob[2], chk2$cell_diagnostics$conditional_prob[2], naive_marg[2]
  )
)

spec2 <- specify_target(
  dgp2,
  effect         = "arm",
  scale          = "pp",
  prob_estimand  = "marginal",
  criterion      = "power",
  power_min      = 0.80,
  prob_se_method = "joint_delta",
  engine         = "bayes_laplace"
)
n_vec2 <- c(48L, 72L, 96L, 120L)
val_rows2 <- lapply(seq_along(n_vec2), function(k) {
  ns <- n_vec2[k]
  vd <- validate_design_point(spec2, n_subj = ns, n_item = 32L, nsim = 140L, seed = 820000L + k * 10L)
  vd$n_subj <- ns
  vd$n_item <- 32L
  vd
})
val_df2 <- do.call(rbind, val_rows2)

p2_left <- ggplot(scale_df2, aes(x = cell, y = prob * 100, fill = quantity)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64, color = "#0f172a", linewidth = 0.3) +
  geom_text(
    aes(label = sprintf("%.1f%%", prob * 100)),
    position = position_dodge(width = 0.72), vjust = -0.45, size = 3.2, fontface = "bold"
  ) +
  scale_fill_manual(
    name = NULL,
    values = c("#0284c7", "#38bdf8", "#f97316")
  ) +
  guides(fill = guide_legend(nrow = 2L)) +
  coord_cartesian(ylim = c(0, 88)) +
  labs(
    title = "A. Scale Non-Collapsibility (SD_subj = 1.05, SD_item = 0.80)",
    subtitle = sprintf(
      "Calibrated: +12.0 pp (beta = %.3f) vs. Naive qlogis: +%.1f pp (beta = %.3f)",
      dgp2$beta[2], (naive_marg[2] - naive_marg[1]) * 100, dgp2$beta_naive_logit[2]
    ),
    x = NULL, y = "Response Probability (%)"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

val_plot_df2 <- rbind(
  data.frame(
    n_subj = val_df2$n_subj[val_df2$condition == "planned"],
    metric = "Planned Power (12 pp)",
    val = val_df2$detection[val_df2$condition == "planned"],
    lo = val_df2$detection_lo[val_df2$condition == "planned"],
    hi = val_df2$detection_hi[val_df2$condition == "planned"]
  ),
  data.frame(
    n_subj = val_df2$n_subj[val_df2$condition == "planned"],
    metric = "95% CI Coverage",
    val = val_df2$coverage[val_df2$condition == "planned"],
    lo = val_df2$coverage_lo[val_df2$condition == "planned"],
    hi = val_df2$coverage_hi[val_df2$condition == "planned"]
  ),
  data.frame(
    n_subj = val_df2$n_subj[val_df2$condition == "null"],
    metric = "Null Type I Error (0 pp)",
    val = val_df2$detection[val_df2$condition == "null"],
    lo = val_df2$detection_lo[val_df2$condition == "null"],
    hi = val_df2$detection_hi[val_df2$condition == "null"]
  )
)

p2_right <- ggplot(val_plot_df2, aes(x = n_subj, y = val * 100, color = metric, shape = metric)) +
  geom_hline(yintercept = c(5, 80, 95), linetype = "dashed", color = "#94a3b8", linewidth = 0.45) +
  geom_errorbar(aes(ymin = lo * 100, ymax = hi * 100), width = 3.0, linewidth = 0.65, alpha = 0.8) +
  geom_line(linewidth = 0.95) +
  geom_point(size = 3.0) +
  scale_color_manual(
    name = NULL,
    values = c("Planned Power (12 pp)" = "#2563eb", "95% CI Coverage" = "#059669", "Null Type I Error (0 pp)" = "#dc2626")
  ) +
  scale_shape_manual(name = NULL, values = c(16, 17, 15)) +
  guides(color = guide_legend(nrow = 2L), shape = guide_legend(nrow = 2L)) +
  scale_x_continuous(breaks = n_vec2) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(
    title = "B. Incomplete-Block Validation (16 of 32 Items/Subj)",
    subtitle = "Paired null & planned operating characteristics (validate_design_point)",
    x = "Recruited Subjects (N_subj; 4% dropout + 5% trial missingness)",
    y = "Empirical Rate (%)"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

save_two_panel("man/figures/gallery_02_incomplete_block_dgp.png", p2_left, p2_right)

# ------------------------------------------------------------------------------
# RECIPE 3: Reading Comprehension & Explanatory IRT — Testlets & Multi-Scenario
# ------------------------------------------------------------------------------
cat("[3/5] Recipe 3: Testlet-Nested Passages & Multi-Scenario Robustness...\n")
des3 <- design_testlet(
  factors      = list(glossing = c("none", "marginal_gloss")),
  testlet_size = 4L,
  contrast     = "sum",
  attrition    = 0.03
)
probs3 <- cell_probs_from_diffs(des3, baseline_prob = 0.66, effects_pp = c(glossing = 8))
dgp3 <- calibrate_dgp(
  des3,
  cell_probs             = probs3,
  chance                 = 0.25,
  subj_sd                = 0.75,
  item_sd                = 0.50,
  testlet_sd             = 0.45,
  participant_testlet_sd = 0.30,
  slope_sd_ratio         = 0.20,
  corr_structure         = 0.10
)
spec3 <- specify_target(
  dgp3,
  effect         = "glossing",
  scale          = "pp",
  prob_estimand  = "marginal",
  criterion      = "power",
  power_min      = 0.80,
  prob_se_method = "joint_delta",
  engine         = "bayes_laplace"
)
scens3 <- scenario_dgps(
  dgp3,
  effect          = "glossing",
  effect_diffs_pp = c(conservative_6pp = 6, standard_8pp = 8, optimistic_10pp = 10)
)
res3 <- optimize_sample_size(
  spec3,
  n_subj_range      = c(24, 72),
  n_item_range      = c(16, 48),
  prefer            = "balanced",
  scenarios         = scens3,
  standard_scenario = "standard_8pp",
  method            = "grid",
  nsim_init         = 50,
  nsim_final        = 150,
  nsim_topup        = 60,
  workers           = workers_n,
  seed              = 830303,
  verbose           = FALSE
)

se3 <- res3$robust$scenario_evaluations
# Slice at fixed N_subj = optimal_robust$n_subj (or median N_subj) to show item/passage curve
slice_subj3 <- res3$robust$optimal_standard$n_subj[1L]
se3_slice <- subset(se3, n_subj == slice_subj3)
if (nrow(se3_slice) < 3L) {
  top_s <- as.integer(names(sort(table(se3$n_subj), decreasing = TRUE))[1L])
  slice_subj3 <- top_s
  se3_slice <- subset(se3, n_subj == slice_subj3)
}

p3_left <- ggplot(se3, aes(x = n_trials_scheduled, y = power * 100, color = scenario, shape = scenario)) +
  geom_hline(yintercept = 80, linetype = "dashed", color = "#64748b", linewidth = 0.6) +
  geom_point(size = 2.6, alpha = 0.88) +
  geom_errorbar(aes(ymin = power_mc_lo * 100, ymax = power_mc_hi * 100), width = 35, alpha = 0.45) +
  scale_color_manual(
    name = "Scenario",
    values = c("conservative_6pp" = "#dc2626", "standard_8pp" = "#2563eb", "optimistic_10pp" = "#059669")
  ) +
  scale_shape_manual(
    name = "Scenario",
    values = c("conservative_6pp" = 16, "standard_8pp" = 15, "optimistic_10pp" = 17)
  ) +
  coord_cartesian(ylim = c(20, 100)) +
  labs(
    title = "A. Multi-Scenario Power (Testlets: 4 Items/Passage)",
    subtitle = sprintf(
      "Standard 8pp: (%d, %d) vs. Robust: (%d, %d) [+%d subj, +%d items]",
      res3$robust$optimal_standard$n_subj, res3$robust$optimal_standard$n_item,
      res3$robust$optimal_robust$n_subj, res3$robust$optimal_robust$n_item,
      res3$robust$delta_subj, res3$robust$delta_item
    ),
    x = "Scheduled Trials (N_subj x N_item; 4AFC Reading Passages)",
    y = "Verified Power (%) [95% Wilson CI]"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

sens3 <- res3$sensitivity
p3_right <- ggplot(sens3, aes(x = factor(sprintf("%.2fx", slope_multiplier)), y = power * 100, fill = decision_status)) +
  geom_hline(yintercept = 80, linetype = "dashed", color = "#64748b", linewidth = 0.6) +
  geom_col(width = 0.52, color = "#0f172a", linewidth = 0.35) +
  geom_errorbar(aes(ymin = power_mc_lo * 100, ymax = power_mc_hi * 100), width = 0.16, linewidth = 0.65) +
  geom_text(
    aes(label = sprintf("%.1f%%\n(CI %.1f pp)", power * 100, prob_ci_width_mean_pp)),
    vjust = 1.45, color = "white", fontface = "bold", size = 3.2
  ) +
  scale_fill_manual(
    name = "Decision",
    values = c("pass" = "#0284c7", "indeterminate" = "#f59e0b", "fail" = "#ef4444")
  ) +
  coord_cartesian(ylim = c(0, 105)) +
  labs(
    title = sprintf(
      "B. Random-Slope Sensitivity at (%d Subj x %d Items)",
      res3$optimal$n_subj, res3$optimal$n_item
    ),
    subtitle = "Evaluated under conservative 6 pp scenario across slope SD multipliers",
    x = "Random-Slope SD Multiplier",
    y = "Effective Power (%)"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

save_two_panel("man/figures/gallery_03_testlet_robustness.png", p3_left, p3_right)

# ------------------------------------------------------------------------------
# RECIPE 4: Two-Dimensional Cost & Budget Optimization — Allocation Reversal
# ------------------------------------------------------------------------------
cat("[4/5] Recipe 4: Cost & Budget Allocation Reversal (Lab vs. Field IRT)...\n")
cost_lab <- c(subj = 60, item = 8, trial = 0.02, fixed = 200)
cost_irt <- c(subj = 6, item = 85, trial = 0.01, fixed = 200)

res4_lab <- switch_preference(res1, prefer = "min_cost", cost = cost_lab, budget = 3200)
res4_irt <- switch_preference(res1, prefer = "min_cost", cost = cost_irt, budget = 3200)

p4_left <- plot(res4_lab, type = "contour", metric = "cost") +
  labs(
    title = "A. In-Lab Psycholinguistics ($60/subj, $8/item)",
    subtitle = sprintf(
      "Budget $3,200 | Winner: (%d, %d) | Cost = $%.0f | Power = %.1f%%",
      res4_lab$optimal$n_subj, res4_lab$optimal$n_item,
      res4_lab$optimal$total_cost, res4_lab$optimal$power * 100
    )
  )
p4_right <- plot(res4_irt, type = "contour", metric = "cost") +
  labs(
    title = "B. Field IRT Item Piloting ($6/subj, $85/item)",
    subtitle = sprintf(
      "Budget $3,200 | Reverses to: (%d, %d) | Cost = $%.0f | Power = %.1f%%",
      res4_irt$optimal$n_subj, res4_irt$optimal$n_item,
      res4_irt$optimal$total_cost, res4_irt$optimal$power * 100
    )
  )
save_two_panel("man/figures/gallery_04_cost_reversal.png", p4_left, p4_right)

# ------------------------------------------------------------------------------
# RECIPE 5: Psychometrics — Random-Item vs. Fixed-Item-Bank Rasch Calibration
# ------------------------------------------------------------------------------
cat("[5/5] Recipe 5: Random-Item vs. Fixed-Item-Bank Rasch Calibration...\n")
des5 <- design_crossed(list(), allocation = "crossed")
dgp5 <- calibrate_dgp(des5, cell_probs = 0.65, chance = 0, subj_sd = 0.85, item_sd = 0.65)
rasch_grid5 <- expand.grid(n_subj = c(60L, 90L, 120L, 160L), n_item = c(16L, 24L, 36L))
rasch_rand5 <- plan_rasch_precision(
  dgp5,
  designs   = rasch_grid5,
  max_error = c(ability = 0.50, difficulty = 0.30, item_sd = 0.22),
  assurance = 0.80,
  nsim      = 120L,
  engine    = "cpp_laplace",
  workers   = workers_n,
  seed      = 850501L
)

b_bank5 <- c(-3.0, seq(-1.1, 1.1, length.out = 14L), 3.0)
n_bank_subj5 <- c(40L, 80L, 120L, 180L, 240L)
bank_mml5 <- plan_rasch_fixed_bank(
  item_difficulties = b_bank5,
  n_subj            = n_bank_subj5,
  person_sd         = 0.85,
  max_error         = c(difficulty = 0.45, difficulty_max_abs_error = 0.60),
  method            = "fixed_item",
  engine            = "mml",
  assurance         = 0.80,
  nsim              = 120L,
  workers           = workers_n,
  seed              = 850502L
)
bank_eb5 <- plan_rasch_fixed_bank(
  item_difficulties = b_bank5,
  n_subj            = n_bank_subj5,
  person_sd         = 0.85,
  max_error         = c(difficulty = 0.45, difficulty_max_abs_error = 0.60),
  method            = "random_item",
  engine            = "cpp_laplace",
  assurance         = 0.80,
  nsim              = 120L,
  workers           = workers_n,
  seed              = 850502L
)

rj5 <- subset(rasch_rand5$summary, target == "joint")
p5_left <- ggplot(rj5, aes(x = n_subj, y = attainment * 100, color = factor(n_item), group = factor(n_item))) +
  geom_hline(yintercept = 80, linetype = "dashed", color = "#64748b", linewidth = 0.6) +
  geom_errorbar(aes(ymin = mc_lo * 100, ymax = mc_hi * 100), width = 4.0, alpha = 0.65) +
  geom_line(linewidth = 1.0) +
  geom_point(aes(shape = decision), size = 3.2) +
  scale_color_manual(name = "Items (N_item)", values = c("16" = "#dc2626", "24" = "#f59e0b", "36" = "#0284c7")) +
  scale_shape_discrete(name = "Decision") +
  scale_x_continuous(breaks = c(60, 90, 120, 160)) +
  coord_cartesian(ylim = c(0, 102)) +
  labs(
    title = "A. Random-Item Rasch Precision",
    subtitle = "Joint assurance: Ability <= 0.50, Diff <= 0.30, Item SD <= 0.22",
    x = "Number of Examinees (N_subj)",
    y = "Joint Attainment Rate (%)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.spacing.y = grid::unit(1, "pt"),
    plot.title = element_text(face = "bold")
  )

b_comb5 <- rbind(
  transform(subset(bank_mml5$summary, target %in% c("difficulty", "difficulty_max_abs_error")),
            estimator = "Fixed-Item MML"),
  transform(subset(bank_eb5$summary, target %in% c("difficulty", "difficulty_max_abs_error")),
            estimator = "Random-Item EB")
)
b_comb5$series <- paste0(
  b_comb5$estimator, " — ",
  ifelse(b_comb5$target == "difficulty", "Bank RMSE <= 0.45", "All 16 Items |err| <= 0.60")
)

p5_right <- ggplot(b_comb5, aes(x = n_subj, y = attainment * 100, color = series, linetype = series)) +
  geom_hline(yintercept = 80, linetype = "dashed", color = "#64748b", linewidth = 0.6) +
  geom_line(linewidth = 0.95) +
  geom_point(size = 2.8) +
  scale_color_manual(
    name = NULL,
    values = c(
      "Fixed-Item MML — Bank RMSE <= 0.45" = "#0284c7",
      "Random-Item EB — Bank RMSE <= 0.45" = "#059669",
      "Fixed-Item MML — All 16 Items |err| <= 0.60" = "#dc2626",
      "Random-Item EB — All 16 Items |err| <= 0.60" = "#7c3aed"
    )
  ) +
  scale_linetype_manual(
    name = NULL,
    values = c(
      "Fixed-Item MML — Bank RMSE <= 0.45" = "solid",
      "Random-Item EB — Bank RMSE <= 0.45" = "solid",
      "Fixed-Item MML — All 16 Items |err| <= 0.60" = "twodash",
      "Random-Item EB — All 16 Items |err| <= 0.60" = "twodash"
    )
  ) +
  scale_x_continuous(breaks = n_bank_subj5) +
  coord_cartesian(ylim = c(0, 102)) +
  guides(color = guide_legend(ncol = 1L), linetype = guide_legend(ncol = 1L)) +
  labs(
    title = "B. Fixed 16-Item Bank (+/-3.0 Logit Anchors)",
    subtitle = "Overall RMSE converges fast; worst-item max error needs larger N_subj",
    x = "Number of Examinees (N_subj)",
    y = "Attainment Rate (%)"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", legend.text = element_text(size = 8.8), plot.title = element_text(face = "bold"))

save_two_panel("man/figures/gallery_05_rasch_calibration.png", p5_left, p5_right)

# ------------------------------------------------------------------------------
# Export JSON payload for Interactive HTML Gallery
# ------------------------------------------------------------------------------
ver1_cols <- c(
  "n_subj", "n_item", "n_trials_scheduled", "n_trials_expected",
  "power", "power_mc_lo", "power_mc_hi",
  "prob_ci_width_mean_pp", "aipe_prob", "joint_prob",
  "valid_rate", "singular_rate", "decision_status"
)
ver1_df <- res1$verified[order(res1$verified$n_item, res1$verified$n_subj), ver1_cols, drop = FALSE]
par1_keys <- paste(res1$verified_pareto$n_subj, res1$verified_pareto$n_item)
ver1_df$is_pareto <- paste(ver1_df$n_subj, ver1_df$n_item) %in% par1_keys

payload <- list(
  recipe1 = list(
    optimal = as.list(res1$optimal[1L, ver1_cols, drop = FALSE]),
    grid = ver1_df,
    report_en = report_plan(res1, lang = "en")
  ),
  recipe2 = list(
    scale_shift = scale_df2,
    beta_calibrated = unname(dgp2$beta[2]),
    beta_naive = unname(dgp2$beta_naive_logit[2]),
    naive_diff_pp = (naive_marg[2] - naive_marg[1]) * 100,
    validation = val_df2[, c("condition", "n_subj", "n_item", "nsim", "n_valid", "detection", "detection_lo", "detection_hi", "coverage", "coverage_lo", "coverage_hi")]
  ),
  recipe3 = list(
    optimal_standard = as.list(res3$robust$optimal_standard[1L, c("n_subj", "n_item", "n_trials_scheduled", "power", "prob_ci_width_mean_pp")]),
    optimal_robust = as.list(res3$robust$optimal_robust[1L, c("n_subj", "n_item", "n_trials_scheduled", "power", "prob_ci_width_mean_pp")]),
    delta_subj = res3$robust$delta_subj,
    delta_item = res3$robust$delta_item,
    scenarios = se3[, c("scenario", "n_subj", "n_item", "n_trials_scheduled", "power", "power_mc_lo", "power_mc_hi", "prob_ci_width_mean_pp", "decision_status")],
    sensitivity = sens3[, c("scenario", "slope_multiplier", "power", "power_mc_lo", "power_mc_hi", "prob_ci_width_mean_pp", "decision_status")],
    report_en = report_plan(res3, lang = "en")
  ),
  recipe4 = list(
    lab_optimal = as.list(res4_lab$optimal[1L, c("n_subj", "n_item", "total_cost", "power", "prob_ci_width_mean_pp")]),
    irt_optimal = as.list(res4_irt$optimal[1L, c("n_subj", "n_item", "total_cost", "power", "prob_ci_width_mean_pp")]),
    report_lab_en = report_plan(res4_lab, lang = "en"),
    report_irt_en = report_plan(res4_irt, lang = "en")
  ),
  recipe5 = list(
    random_rasch = rj5[, c("n_subj", "n_item", "nsim", "attainment", "mc_lo", "mc_hi", "decision")],
    fixed_bank = b_comb5[, c("estimator", "target", "n_subj", "n_item", "n_extreme_item_reps", "attainment", "mc_lo", "mc_hi", "decision")]
  )
)

args <- commandArgs(trailingOnly = TRUE)
hit_json <- grep("^--json-out=", args, value = TRUE)
json_out <- if (length(hit_json)) {
  sub("^--json-out=", "", hit_json[1L])
} else {
  file.path(tempdir(), "gallery_payload.json")
}
writeLines(to_json_val(payload), json_out)
cat(sprintf("=== Gallery build completed! JSON payload written to: %s ===\n", json_out))
print(list.files("man/figures", full.names = TRUE))
