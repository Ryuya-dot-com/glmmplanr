# glmmplanr — Two-Dimensional Sample Size and Precision Planning for Crossed Binary Designs and Explanatory Item Response Models

[![R-CMD-check](https://github.com/Ryuya-dot-com/glmmplanr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Ryuya-dot-com/glmmplanr/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-0284c7.svg)](https://github.com/Ryuya-dot-com/glmmplanr/blob/main/LICENSE.md)
[![R >= 4.1.0](https://img.shields.io/badge/R-%3E%3D%204.1.0-276DC3.svg)](https://www.r-project.org/)
[![Backend: RcppEigen](https://img.shields.io/badge/Backend-C%2B%2B%20%2F%20RcppEigen-059669.svg)](src/glmm_laplace.cpp)

`glmmplanr` helps answer two questions. How many learners and words would you need to detect an improvement from 60% to 70% correct? How many examinees would you need to estimate the difficulty of every word on a fixed list? The input is a planned study with a **binary response on each trial** (correct/incorrect or yes/no). A three-choice task uses correct or incorrect scoring with a guessing probability of `1/3` (`chance = 0` for binary endorsement/judgment tasks). Ordinal Likert ratings, reaction times, counts, and distractor-specific choices are outside the implemented response models; see [Scope Across Research Designs](inst/validation-guide.md#cross-field-designs).

---

## Which Function Should I Use?

Most workflows follow four steps: **1. Declare design** → **2. Calibrate probabilities** → **3. Specify target** → **4. Search or verify sample sizes**.

Planning uses built-in estimators by default: `bayes_laplace` for crossed
GLMM effects, `cpp_laplace` for random-item Rasch recovery, and `mml` for
fixed-item Rasch banks. External engines can be selected explicitly for
comparison or to reproduce earlier analyses. Each recommendation applies
to its specified model, estimator, and success criterion.

| Your planning goal | Start here |
|---|---|
| Find combinations of **participants × items** that achieve target power or CI width | `design_crossed()` → `calibrate_dgp()` → `specify_target()` → `optimize_sample_size()` |
| Check Type I error, power, and CI coverage at **one specific design** `(N_subj, N_item)` | `validate_design_point()` |
| Plan for **passages / testlets** (multiple questions nested per reading/listening text) | `design_testlet()` (replaces `design_crossed()`) |
| Plan participants and items for **random-item Rasch** recovery (new people & new items) | `plan_rasch_precision()` |
| Plan examinees to calibrate a **fixed list of Rasch items** (e.g., an operational word list) | `plan_rasch_fixed_bank()` |
| Inspect all-correct, all-incorrect, or unobserved participants/items in trial data | `response_patterns()` |

For a large item bank or stimulus pool, first define admissible candidate
sets, then compare their measurement or experimental precision and required
sample sizes. The [150-item worked example](inst/validation-guide.md#item-pool-planning)
distinguishes a fixed short form from rotating partial presentation of the
whole bank, and known-bank scoring precision from item recalibration.

---

## Installation

`glmmplanr` compiles C++ code via `RcppEigen`, so you need a standard C++ toolchain (Rtools on Windows, Xcode Command Line Tools on macOS, or `g++`/`clang++` on Linux).

```r
# install.packages("remotes")
remotes::install_github("Ryuya-dot-com/glmmplanr")

# Optional: also build local vignettes (requires knitr, rmarkdown, and pandoc)
remotes::install_github("Ryuya-dot-com/glmmplanr", build_vignettes = TRUE)
```

Start with the [Planning Workflows Tutorial](vignettes/planning-workflows.Rmd) to map a paper's mean accuracies, odds ratios or logit coefficients, random-effect SDs, and model formula to the package arguments. It includes separate routes for condition means and GLMM tables, with checks on the resulting inputs. You can also run `vignette("planning-workflows", package = "glmmplanr")` if you built vignettes. For statistical benchmarks and a comparison with 19 related R packages, see [Scope, Validation & Related Packages](inst/validation-guide.md#related-packages).

---

## Quickstart (Runs in Seconds)

Suppose learners complete a 3-choice vocabulary task (`chance = 1/3`) under two counterbalanced conditions (`incidental` vs. `instruction`). You expect **60% vs. 70% population-averaged accuracy** (+10 percentage points) and want **at least 80% power**.

```r
library(glmmplanr)

# 1. Declare a counterbalanced within-subject / within-item Latin-square design
des <- design_crossed(
  factors    = list(learning = c("incidental", "instruction")),
  allocation = "latin_square"
)

# 2. Convert expected population accuracies (60% vs. 70%) into GLMM parameters
dgp <- calibrate_dgp(
  des,
  cell_probs     = c(incidental = 0.60, instruction = 0.70),
  chance         = 1 / 3,
  subj_sd        = 0.8,
  item_sd        = 0.6,
  slope_sd_ratio = 0
)
check_dgp(dgp, plot = FALSE)$cell_diagnostics[c("cell", "target_prob", "quad_prob")]

# 3. Define the target: >= 80% power on the marginal percentage-point scale
spec <- specify_target(
  dgp,
  effect         = "learning",
  scale          = "pp",
  criterion      = "power",
  power_min      = 0.80,
  prob_se_method = "joint_delta",  # Fast analytical SE for quick exploration
  engine         = "bayes_laplace",
  analysis_re    = "intercept_only"
)

# 4. Evaluate candidate (N_subj x N_item) designs (small nsim for a fast demo)
plan <- optimize_sample_size(
  spec,
  n_subj_range = c(24, 32, 48),
  n_item_range = c(12, 16, 24),
  method       = "grid",
  nsim_final   = 40L,
  workers      = 1L,
  seed         = 417L,
  verbose      = FALSE,
  control      = optimize_control(nsim_sensitivity = 0L)
)

plan                                    # Recommended design and decision status
plan$verified_pareto[c("n_subj", "n_item", "power", "power_mc_lo", "decision_status")]
plot(plan, type = "contour")            # 2D power surface + Pareto frontier
```

### How to Read the Result

- **`decision_status` (Target attainment only)**: Uses a Wilson Monte Carlo confidence interval (default 95%) for attainment of your power, width, or ROPE target.
  - **`pass`**: Even the lower Monte Carlo bound (`power_mc_lo`) reaches your target (e.g., $\ge 80\%$). Only `pass` designs are confirmed recommendations.
  - **`indeterminate`**: The Monte Carlo interval straddles the target; increase `nsim_final`, use `topup_verification()` (which applies geometric sequential alpha-spending across independent rounds by default), or confirm a chosen design once via `validate_design_point()`.
  - **`fail`**: The upper Monte Carlo bound falls below the target.
- **Target attainment vs. inferential validity**: A `pass` status confirms only your specified power/precision target—it does **not** guarantee nominal interval coverage or null Type I error, because an estimator that underestimates standard errors can produce narrow intervals and high detection. Console summaries, `as_gt()`, and `report_plan()` report **Inferential validity** separately (flagging low Stage-2 coverage and pointing to `validate_design_point()`).
- **`power` vs. `conditional_power`**: By default (`failure_policy = "effective"`), non-converged or separated fits count as misses in `power`, while `conditional_power` reports power among valid fits only.
- **`plan$verified_pareto`**: Lists non-dominated passing `(N_subj, N_item)` trade-offs (reducing participants requires more items, and vice versa). Switch priorities without re-simulating via `switch_preference(plan, prefer = "fewer_subj")` or `"fewer_items"`. Near the pass/fail threshold, the smallest passing pair can shift across search seeds—console summaries and `report_plan()` flag narrow Monte Carlo margins and suggest reporting the candidate range alongside a higher-margin alternative.
- **Final confirmation**: Before pre-registration, confirm your chosen design on a fresh seed with a larger `nsim` and your planned final interval method (e.g., `prob_se_method = "bootstrap"`) via `validate_design_point()`, which checks both Type I error under the null and power/coverage under the planned effect (or planned-condition precision alone when estimating a single condition's accuracy).

### Adapting the Quickstart Assumptions to Your Own Study

- **Translating `subj_sd` and `item_sd` without prior GLMM tables**: Printing `dgp` shows the 95% participant and item accuracy spans implied by your SDs. If a prior paper does not report `VarCorr()`, pass expected 95% accuracy ranges directly via `subj_prob_range = c(0.42, 0.85)` and `item_prob_range = c(0.45, 0.86)` (or convert with `prob_range_to_sd()`), and compare low/moderate/high variance brackets via `scenarios`—see [Choosing starting inputs when the source paper omits `VarCorr`](vignettes/planning-workflows.Rmd).
- **Do not keep `slope_sd_ratio = 0` uncritically**: Setting `slope_sd_ratio = 0` assumes the condition effect is *identical* for every learner and every word, which can yield optimistic sample sizes in within-participant or within-item designs. Unless you intentionally plan an intercept-only approximation and accept mismatch risk, test `slope_sd_ratio = 0.25`–`0.35` with `analysis_re = "match_dgp"` (or run `stress_test_models()`).

---

## Visual Recipe Gallery & Interactive Explorer

Explore four worked examples with figures, R scripts, and sample size justifications covering **(1) Psycholinguistic $2 \times 2$ Latin-Square 3AFC Experiments**, **(2) Between-Subject Incomplete-Block Designs & Logistic Non-Collapsibility**, **(3) Testlet-Nested Reading Passages & Multi-Scenario Robustness**, and **(4) Random-Item vs. Fixed-Item-Bank Rasch Calibration**.

- **GitHub Visual Recipe Gallery**: [`inst/gallery/GALLERY.md`](inst/gallery/GALLERY.md)
- **Interactive HTML Explorer**: [**Launch Online Demo**](https://ryuya-dot-com.github.io/glmmplanr/) or open locally from R via:
  ```r
  browseURL(system.file("gallery", "index.html", package = "glmmplanr"))
  ```

![2D Power & Precision Surface and Non-Dominated Pareto Frontier](man/figures/gallery_01_latin_square_3afc.png)

---

## Planning Rasch / 1PL Models

Choose the Rasch function based on whether you generalize to **new items** or calibrate **one fixed list**:

| Question | Function | Items across replications |
|---|---|---|
| How many examinees & items to recover person ability, item difficulty, and population item SD? | `plan_rasch_precision()` | Drawn anew on every replication |
| How many examinees to calibrate a specific list of operational items (including extreme items)? | `plan_rasch_fixed_bank()` | Held fixed across replications |

### Example A: Fixed Item Bank (`plan_rasch_fixed_bank()`)

```r
library(glmmplanr)

# Fixed 8-item list with true difficulties from -2 to +2 logits
bank_diffs <- c(-2.0, -1.0, -0.6, -0.2, 0.2, 0.6, 1.0, 2.0)

bank <- plan_rasch_fixed_bank(
  item_difficulties = bank_diffs,
  n_subj            = c(40L, 100L),
  person_sd         = 0.8,
  max_error         = c(difficulty = 0.35, difficulty_max_abs_error = 0.60),
  method            = "fixed_item",  # Marginal Maximum Likelihood ("mml" or "tam")
  assurance         = 0.80,
  nsim              = 32L,
  seed              = 417L
)

subset(as.data.frame(bank), target == "joint")
plot(bank, type = "items", metric = "bias")  # Per-item bias and finite-estimate counts
report_plan(bank)
```

### Example B: Random-Item Generalization (`plan_rasch_precision()`)

```r
des_r <- design_crossed(list(), allocation = "crossed")
dgp_r <- calibrate_dgp(des_r, cell_probs = 0.65, chance = 0, subj_sd = 0.8, item_sd = 0.6)

rasch_plan <- plan_rasch_precision(
  dgp_r,
  designs   = data.frame(n_subj = c(144, 72, 144), n_item = c(12, 24, 24)),
  max_error = c(ability = 0.50, difficulty = 0.30, item_sd = 0.20),
  assurance = 0.80,
  nsim      = 64L,
  engine    = "cpp_laplace",
  seed      = 518307L
)

subset(as.data.frame(rasch_plan), target == "joint")
plot(rasch_plan)
rasch_plan$pareto
```

Unobserved items and constant (0% or 100% correct) items in fixed-item MML remain in the assurance denominator as misses; see `?plan_rasch_fixed_bank` and [Fixed-Bank Comparison with TAM](inst/validation-guide.md#rasch-fixed-bank).

---

## Full Workflow: 2×2 3AFC Precision, Multi-Scenario Robustness & Pareto Switching

For a realistic study plan—combining a $2 \times 2$ Latin square, 3AFC guessing, simultaneous power + CI-width targets, and robustness across conservative/standard/optimistic effect sizes—adjust `workers` to your CPU cores. Note that `prob_se_method = "bootstrap"` (default for `scale = "pp"`) runs 99 refits per simulated study and takes substantially longer than `prob_se_method = "joint_delta"`; use `"joint_delta"` for rapid screening and `"bootstrap"` for final confirmation.

```r
library(glmmplanr)

# 1. Balanced 2x2 Latin-square design (4 lists, 5% random trial attrition)
des <- design_crossed(
  factors = list(
    exposure = c("incidental", "intentional"),
    spacing  = c("massed", "distributed")
  ),
  allocation = "latin_square",
  contrast   = "sum",
  attrition  = 0.05
)

# 2. Specify cell accuracies from percentage-point effects & calibrate 3AFC DGP
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

# 3. Target both >= 80% power AND 95% CI width <= 12 pp
spec <- dgp |>
  specify_target(
    effect         = "exposure",
    scale          = "pp",
    criterion      = "both",
    power_min      = 0.80,
    ci_width       = 12,
    prob_se_method = "joint_delta",  # Switch to "bootstrap" for final confirmation
    engine         = "bayes_laplace"
  )

# 4. Optimize across +5 / +8 / +10 pp scenarios
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
    nsim_final        = 400,
    workers           = 4,            # Match to your available CPU cores
    seed              = 20261119,
    control           = optimize_control(n_eval_init = 16L)
  )

# Instantaneously switch trade-off preferences without re-simulating
switch_preference(plan_res, prefer = "fewer_subj")
switch_preference(plan_res, prefer = "fewer_items")

# Tables, plots, CSV export, and English Methods text
summary(plan_res)
head(as.data.frame(plan_res))
as_gt(plan_res)
plot(plan_res, type = "contour", metric = "prob_ci_width_pp")
plot(plan_res, type = "pareto")
plot(plan_res, type = "scenarios")
# plot(plan_res, type = "surface3d", metric = "power")  # Requires optional 'plotly'
report_plan(plan_res)

# Independent null and planned confirmation on a fresh seed
check <- validate_design_point(spec, n_subj = 48, n_item = 32, nsim = 200, seed = 74011)
plot(check)
plot(check, type = "patterns")
```

The scenario plot compares the designs selected under the standard scenario
and across all scenarios. Pass a data frame with `n_subj` and `n_item` as
`candidates` to compare specific verified pairs. The validation plot shows
Monte Carlo intervals for power, Type I error, and coverage with the valid-fit
counts. The pattern plot reports the frequency of extreme scores and
unobserved units by participant and item.

Static figures use `ggplot2::theme_bw()` and hide automatic titles,
subtitles and captions by default. Axis labels, legends, scenario labels,
Monte Carlo intervals and valid-fit counts remain visible. Use
`plot(check, show_annotations = TRUE)` to include the generated explanations,
or `plot(check) + ggplot2::labs(title = "My study", caption = "My caption")`
to supply your own. The same option applies to GLMM and Rasch plans and
`check_dgp()`; interactive 3D plots use a white background.

The [visualization rationale](inst/validation-guide.md#visualization)
connects the plots to research on presenting simulation results.

---

## Key Capabilities

- **Flexible Crossed, Incomplete-Block & Testlet Designs (`design_crossed()`, `design_testlet()`)**: Supports Latin-square, between-subject (`between_subj`), between-item (`between_item`), full crossed, partial item-bank presentation (`items_per_subject`), and passage/testlet nesting with optional local-dependence stress testing (`participant_testlet_sd`). Tracks trial-level missingness (`attrition`) and whole-participant dropout (`participant_dropout`) separately.
- **Exact Marginal Probability & $m$-AFC Calibration (`calibrate_dgp()`, `check_dgp()`)**: Inverts population-averaged cell probabilities into conditional GLMM fixed effects via exact 1D Gauss-Hermite quadrature over the combined cell normal variance, correcting logistic non-collapsibility and lower guessing asymptotes (`chance = 1/m`).
- **Percentage-Point Power, Precision (AIPE) & ROPE (`specify_target()`)**: Plan on `"pp"`, `"prob"`, or `"logit"` scales using `"power"`, `"ci_width"`, `"both"`, `"joint"` (same-replication detection + width), or `"rope"` criteria, with either fast analytical `joint_delta` SEs or C++ batch parametric `bootstrap` SEs.
- **Multiple Backends & Conservative Monte Carlo Decisions**: Built-in compiled C++ engines (`cpp_laplace`, regularized `bayes_laplace` using Half-Student-$t_3$ log-SD MAP updates) plus `lme4`, `glmmTMB`, and Stan HMC (`stan_hmc`). Offers optional Bonferroni multiplicity adjustment (`mc_adjustment = "bonferroni"`), independent fresh-batch top-ups (`topup_policy = "independent"`), and `mirai` parallel workers.
- **Extreme Response Pattern Diagnostics (`response_patterns()`)**: Summarizes `all_correct`, `all_incorrect`, and `unobserved` participants/items alongside observation counts (`n_observed`), while unpenalized GLMM engines automatically flag separation when an entire fixed-effect cell has only one outcome.

For detailed coverage checks, 3AFC estimator comparisons (Laplace vs. Stan HMC), and comparisons with `lme4`, `TAM`, `PowRPriori`, `simr`, and `mlpwr`, see [Scope, Statistical Validation, and Related R Packages](inst/validation-guide.md).

---

## Validation Status & Modeling Assumptions

- **Benchmarked vs. exploratory scope**: External cross-validation and pre-specified confirmation benchmarks currently cover complete-presentation canonical logit GLMMs (`lme4`, `simr`, `PowRPriori`, floor-checked), complete-presentation $2 \times 2$ 3AFC conditional-logit models (`rstan`, `bayes_laplace`, including an $N = 512$ strict conditional-null calibration check and boundary analysis), and complete-presentation random-item and fixed-bank Rasch models (`TAM`, `lme4`, floor-checked). **Partial presentation (`items_per_subject`), sparse designs with correlated random slopes, testlet local dependence (`design_testlet()`), and higher-order interactions are currently exploratory**—verify coverage and null Type I error at your planned design with `validate_design_point()` (and an external estimator when feasible). Full reproduction scripts and CSVs live in the repository [`validation/` directory](https://github.com/Ryuya-dot-com/glmmplanr/tree/main/validation); see the [Validation Status Matrix](inst/validation-guide.md#validation-status).
- **Generating-model assumptions & custom DGPs**: Built-in data generation assumes multivariate normal participant/item random effects, a shared lower guessing asymptote (`chance`), and random trial attrition / participant dropout independent of latent ability and responses (MCAR). Multi-scenario planning (`scenario_dgps()`) assesses parameter sensitivity within this model family, not structural misspecification; pass a custom `data_fn(schedule, spec)` to `evaluate_design_point()` or `validate_design_point()` to test ability-dependent dropout, item-specific guessing, or non-normal random effects. See [Generating-Model Assumptions and Custom Misspecification Testing](inst/validation-guide.md#modeling-assumptions).

---

## License & Disclaimer

`glmmplanr` is open-source software released under the **[MIT License](https://github.com/Ryuya-dot-com/glmmplanr/blob/main/LICENSE.md)** (Copyright © 2026 Ryuya Komuro).

You are completely free to use, copy, modify, merge, publish, and distribute this software for academic, educational, or commercial purposes. **However, the software is provided "AS IS", without warranty of any kind, express or implied.** Sample size and precision outputs depend on user-supplied data-generating assumptions, model specifications, and Monte Carlo tolerances; **users assume sole responsibility** for validating assumptions, interpreting simulation diagnostics, and making design, ethical, or financial decisions in their own research.
