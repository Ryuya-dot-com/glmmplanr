# Statistical Validation & Reproducibility Suite (`glmmplanr` 0.1.0)

This directory archives the simulation scripts, trial-level CSVs, and summary tables supporting the benchmarks summarized in [`inst/validation-guide.md`](../inst/validation-guide.md). It is tracked in Git for full reproducibility alongside the package source, and excluded from CRAN source tarballs via `.Rbuildignore` to keep package downloads compact (the core 3AFC confirmation script and CSVs are also shipped inside `inst/validation/`).

## Reference Environment & Isolated Reproduction

All recent `glmmplanr 0.1.0` confirmation runs were executed under:

- **R**: `4.6.1` (`aarch64-apple-darwin25.4.0`)
- **Packages**: `glmmplanr 0.1.0`, `lme4 2.1.0`, `TAM 4.3.25`, `rstan 2.32.7`, `cmdstanr 0.9.0`, `simr`, `PowRPriori`, `mlpwr`

Use the isolated runner from the repository root to install the current source into a clean temporary library, verify package resolution, record `sessionInfo()` and package versions, and run any script without mutating your global R library:

```sh
sh validation/run-isolated.sh coverage_recheck
sh validation/run-isolated.sh marginal_null_recheck
sh validation/run-isolated.sh three_choice_confirmation
sh validation/run-isolated.sh three_choice_boundary_and_null
sh validation/run-isolated.sh package_calibration
sh validation/run-isolated.sh cost_mlpwr_comparison
sh validation/run-isolated.sh rasch_fixed_bank
```

When comparing two runs on identical seeds, match each trial by `(condition/effect, design, seed)`, verify convergence/validity per pair, record maximum and median absolute differences in point estimates and SEs on valid pairs, and recompute binary pass/fail outcomes from unrounded values. For fresh seeds, interpret rate differences using exact binomial or Wilson 95% Monte Carlo intervals ([Morris et al., 2019](https://doi.org/10.1002/sim.8086)).

---

## 1. Marginal 95% Interval Coverage & Paired Null Checks

- **Scripts & CSVs**: [`coverage_recheck.R`](coverage_recheck.R) ([`coverage_recheck_summary.csv`](coverage_recheck_summary.csv), [`coverage_recheck_fits.csv`](coverage_recheck_fits.csv)), [`marginal_null_recheck.R`](marginal_null_recheck.R) ([`marginal_null_recheck_summary.csv`](marginal_null_recheck_summary.csv), [`marginal_null_recheck_fits.csv`](marginal_null_recheck_fits.csv))
- **Setup**: Canonical binary logit (`chance = 0`), complete presentation `64 x 20`, random intercept SDs `0.5` (participants) and `0.3` (items), 99 parametric bootstrap refits per dataset.
- **Planned contrast** (`60%` vs. `55%`, `-5 pp`, between-subject, 96 independent confirmation datasets, seeds `10001–10096`):
  - Native `bayes_laplace` bootstrap: `96/96` valid, coverage **`95/96` (`99.0%`, exact 95% CI `[94.3%, 100.0%]`)**.
  - External `lme4::glmer` on identical datasets: `96/96` valid, coverage **`94/96` (`97.9%`, `[92.7%, 99.7%]`)**.
  - Both exceeded the pre-specified `90%` lower-bound floor.
- **Paired null checks** (`60%` vs. `60%`, 96 datasets per allocation, seeds `20001–20096` and `21001–21096`, all 96 datasets in denominator; floor-checked rather than equivalence-tested at 5%):
  - Between-subject: `5/96` (`5.2%`, `[1.7%, 11.7%]`) with C++ batch bootstrap (`prob_boot = 99`); `3/96` (`3.1%`, `[0.6%, 8.9%]`) with analytical `joint_delta`.
  - Within-subject Latin square: `7/96` (`7.3%`, `[3.0%, 14.4%]`) with bootstrap; `8/96` (`8.3%`, `[3.7%, 15.8%]`) with `joint_delta`.

---

## 2. Conditional-Logit Uncertainty, Boundary, and Strict Null in 3AFC Tasks

- **Scripts & CSVs**:
  - Paired Laplace vs. RStan HMC & covariance decomposition: [`stan_uncertainty.R`](stan_uncertainty.R), [`three_choice_precision.R`](three_choice_precision.R), [`three_choice_cov_decomposition.R`](three_choice_cov_decomposition.R) ([`three_choice_precision_summary.csv`](three_choice_precision_summary.csv), [`three_choice_cov_decomp_summary.csv`](three_choice_cov_decomp_summary.csv))
  - Sufficiency & 3-condition confirmation (`1,152` trials): [`three_choice_confirmation.R`](three_choice_confirmation.R) ([`three_choice_confirmation_summary.csv`](three_choice_confirmation_summary.csv), [`three_choice_confirmation_fits.csv`](three_choice_confirmation_fits.csv))
  - Boundary / necessity, multi-seed optimizer hit-rate, and high-$N$ strict conditional null (`N = 512`): [`three_choice_boundary_and_null.R`](three_choice_boundary_and_null.R) ([`three_choice_boundary_summary.csv`](three_choice_boundary_summary.csv), [`three_choice_optimizer_boundary.csv`](three_choice_optimizer_boundary.csv), [`three_choice_strict_null_512_summary.csv`](three_choice_strict_null_512_summary.csv))
- **Setup**: Within-subject $2 \times 2$ Latin square (`exposure x spacing`), `chance = 1/3`, correlated random intercepts and slopes (`slope_sd_ratio = 0.30`, `corr_structure = 0.15`), target estimand = conditional `exposure` logit coefficient (`+1.616374` planned; `0` strict conditional null via `null_target_spec()`).
- **SE approximation caveat**: In `src/glmm_laplace.cpp`, the analytical variance-uncertainty correction propagates marginal log-SD curvature diagonally (`J_beta_theta * diag(var_log_tau) * J_beta_theta'`), omitting off-diagonal SD-SD covariances and random-slope correlation uncertainty. On 40 paired `48 x 32` datasets, fixing covariance at true values gave mean Laplace SE `0.201` vs. HMC posterior SD `0.204` (`38/40` valid HMC), whereas estimating covariance gave mean Laplace SE `0.227` vs. HMC posterior SD `0.270`.

### 2A. Sufficiency & 3-Condition Confirmation (`N = 192` per cell, `bayes_laplace`)

| Condition | Design | Valid | 95% Coverage (Exact 95% CI) | Width $\le 1.00$ | Zero Excluded (Exact 95% CI) | Emp. SD / Mean SE | Floor Status |
|---|---|---:|---:|---:|---:|---:|---|
| Planned (`+1.616`) | `72 x 32` | 192/192 | 183/192 (`95.3% [91.3, 97.8]`) | 192/192 | 192/192 (`100%`) | 0.183 / 0.191 | pass |
| Planned (`+1.616`) | `48 x 48` | 192/192 | 187/192 (`97.4% [94.0, 99.1]`) | 192/192 | 192/192 (`100%`) | 0.169 / 0.190 | pass |
| Strict conditional null (`0.000`) | `72 x 32` | 190/192 | 184/190 (`96.8% [93.3, 98.8]`) | 190/192 | 6/192 (`3.1% [1.2, 6.7]`) | 0.154 / 0.168 | pass |
| Strict conditional null (`0.000`) | `48 x 48` | 191/192 | 183/191 (`95.8% [91.9, 98.2]`) | 191/192 | 8/192 (`4.2% [1.8, 8.0]`) | 0.163 / 0.165 | pass |
| Equal-marginal null (`+0.0041`) | `72 x 32` | 190/192 | 185/190 (`97.4% [94.0, 99.1]`) | 190/192 | 6/192 (`3.1% [1.2, 6.7]`) | 0.151 / 0.165 | pass |
| Equal-marginal null (`+0.0041`) | `48 x 48` | 192/192 | 182/192 (`94.8% [90.6, 97.5]`) | 192/192 | 10/192 (`5.2% [2.5, 9.4]`) | 0.163 / 0.165 | pass |

### 2B. High-$N$ Strict Conditional-Null Calibration (`N = 512` per design, $\beta = 0$)

Pre-specified calibration gates: exact 95% Type I error upper bound $\le 7.5\%$ (denominator = all `512` attempted datasets) and exact 95% coverage lower bound $\ge 92.5\%$ (denominator = valid fits).

| Design | Valid | 95% Coverage (Exact 95% CI) | Coverage Gate ($\ge 92.5\%$) | Type I Error (Exact 95% CI) | Type I Gate ($\le 7.5\%$) | Emp. SD / Mean SE |
|---|---:|---:|---|---:|---|---:|
| `48 x 48` | 512/512 | 495/512 (`96.7% [94.7, 98.1]`) | **pass** | 17/512 (`3.32% [1.95, 5.26]`) | **pass** | 0.158 / 0.166 |
| `72 x 32` | 511/512 | 493/511 (`96.5% [94.5, 97.9]`) | **pass** | 18/512 (`3.52% [2.10, 5.50]`) | **pass** | 0.153 / 0.167 |

### 2C. Necessity / Pass-Fail Boundary & Multi-Seed Optimizer Search (`bayes_laplace`)

Adjacent smaller candidates evaluated on `192` fresh datasets each (planned effect `+1.616` logit, width target $\le 1.00$ logit with 80% assurance floor, seeds `7611001–7615192`):

| Design | Trials | Valid | 95% Coverage (Exact 95% CI) | Width $\le 1.00$ (Exact 95% CI) | Mean Width | Width Status |
|---|---:|---:|---:|---:|---:|---|
| `32 x 32` | 1,024 | 191/192 | 185/191 (`96.9% [93.3, 98.8]`) | 10/192 (`5.2% [2.5, 9.4]`) | 1.150 | **below_floor** |
| `48 x 24` | 1,152 | 192/192 | 186/192 (`96.9% [93.3, 98.8]`) | 25/192 (`13.0% [8.6, 18.6]`) | 1.090 | **below_floor** |
| `48 x 32` | 1,536 | 192/192 | 188/192 (`97.9% [94.8, 99.4]`) | 163/192 (`84.9% [79.0, 89.6]`) | 0.928 | **indeterminate** |
| `64 x 24` | 1,536 | 192/192 | 187/192 (`97.4% [94.0, 99.1]`) | 163/192 (`84.9% [79.0, 89.6]`) | 0.930 | **indeterminate** |
| `48 x 40` | 1,920 | 192/192 | 185/192 (`96.4% [92.6, 98.5]`) | 191/192 (`99.5% [97.1, 100.0]`) | 0.823 | **pass** |

Across 5 independent search seeds (`optimize_sample_size()`, `method = "surrogate_gp"`, `prefer = "fewer_trials"`, `nsim_final = 96`, lattice `{32,48,64,72} x {24,32,40,48}`, seeds `7650001–7650005`), the optimizer selected `72 x 24` (`1,728` trials, Stage-2 width attainment `92.7%–96.9%`) on **4/5 seeds** and `48 x 32` (`1,536` trials, Stage-2 width attainment `88.5%`) on **1/5 seed**, with **0/5** failing designs selected.

---

## 3. External Package Calibration (`PowRPriori`, `simr`) & Cost-Constrained Search (`mlpwr`)

### 3A. Marginal Probability Calibration (`PowRPriori` & `simr`)

- **Script & CSVs**: [`package_calibration.R`](package_calibration.R) ([`package_calibration_cells.csv`](package_calibration_cells.csv), [`package_calibration_rates.csv`](package_calibration_rates.csv), [`package_calibration_fits.csv`](package_calibration_fits.csv), [`package_calibration_pow.csv`](package_calibration_pow.csv))
- **Setup**: Canonical logit (`chance = 0`), `64 x 24`, random intercept SDs `0.8` (participants) and `0.6` (items), target marginal rates `72%` vs. `80%` (planned) and `72%` vs. `72%` (null), `192` datasets per cell, two-sided 5% Wald test.
- Passing uncalibrated `72%`/`80%` directly into logit link shrunk true population-averaged rates to `68.68%`/`76.13%` (`7.45 pp`), whereas passing `glmmplanr::calibrate_dgp()` conditional probabilities (`75.65%`/`83.89%`) matched `72.00%`/`80.00%` within `1e-5`.

| Allocation | Input Mode | Planned Detections (`/192`: Native / `simr` / `PowRPriori`) | Null Detections (`/192`, same order) |
|---|---|---:|---:|
| Within-subject (Latin square) | Direct logit | 183 / 182 / — | 9 / 9 / — |
| Within-subject (Latin square) | Calibrated | 188 / 188 / — | 6 / 4 / — |
| Between-subject | Direct logit | 95 / 92 / 96 | 10 / 10 / 8 |
| Between-subject | Calibrated | 115 / 111 / 112 | 9 / 9 / 13 |

### 3B. Aligned 25-Point Cost Search Comparison (`glmmplanr` vs. `mlpwr`)

- **Script & CSVs**: [`cost_mlpwr_comparison.R`](cost_mlpwr_comparison.R) ([`cost_mlpwr_summary.csv`](cost_mlpwr_summary.csv), [`cost_mlpwr_runs.csv`](cost_mlpwr_runs.csv))
- **Setup**: Same Latin-square DGP (`68%` vs. `80%`, `chance = 0`), **identical 25-point discrete lattice** (`{16,24,32,40,48} x {8,12,16,20,24}`), **identical 1,200 total simulation budget**, **5 independent search seeds** per cost regime (`expensive_items`: $10 N_s + 80 N_i$; `expensive_subj`: $80 N_s + 10 N_i$), evaluated on a shared independent `400`-rep test lattice (seed `869000`). Three distinct decision rules are separated explicitly:
  1. `point_80`: Matched 80% point-estimate target (`glmmplanr` Stage-2 `power >= 0.80` vs. `mlpwr` `power = 0.80`).
  2. `high_target_86`: `mlpwr` raised point target `power = 0.860` (exact minimum count `172/200` whose 95% Wilson lower bound reaches `0.80`, namely `0.8050`; note that `171/200 = 0.855` gives `0.7995 < 0.80`, and raising a surrogate point target is *not* equivalent to applying a finite-sample Wilson gate).
  3. `wilson_lo_80`: Explicit 200-rep 95% Wilson lower bound $\ge 0.80$ gate (`glmmplanr` Stage-2 vs. `mlpwr` 1,000-eval search at `power = 0.80` followed by an independent 200-rep Wilson check on the selected point, totaling 1,200 evaluations).

| Regime | Rule | Solver (Budget) | Internal Wilson Pass | Mean Cost [Min–Max] | Out-of-Sample Power (Min) | OOS $\ge 80\%$ / OOS Wilson Pass | Selected Designs (5 Seeds) |
|---|---|---|---:|---:|---:|---:|---|
| Expensive items | `point_80` | Reference grid (10,000) | — | 1,280 | 84.0% (84.0%) | 1/1 / 1/1 | `32x12` |
| Expensive items | `point_80` | `glmmplanr` (1,200) | 1/5 | 1,264 [1,120–1,360] | 85.7% (84.0%) | 5/5 / 5/5 | `32x12`(×3), `48x8`, `40x12` |
| Expensive items | `point_80` | `mlpwr` (1,200) | — | 1,184 [1,120–1,280] | 85.4% (84.0%) | 5/5 / 5/5 | `48x8`(×3), `32x12`(×2) |
| Expensive items | `high_target_86` | `mlpwr` `power=0.86` (1,200) | — | 1,360 [1,360–1,360] | 90.3% (90.3%) | 5/5 / 5/5 | `40x12`(×5) |
| Expensive items | `wilson_lo_80` | `glmmplanr` (1,200) | **5/5** | 1,360 [1,360–1,360] | 90.3% (90.3%) | 5/5 / 5/5 | `40x12`(×5) |
| Expensive items | `wilson_lo_80` | `mlpwr` 1k + Wilson 200 (1,200) | **0/5** | 1,152 [1,120–1,280] | 85.8% (84.0%) | 5/5 / 5/5 | `48x8`(×4), `32x12` |
| Expensive subj | `point_80` | Reference grid (10,000) | — | 1,520 | 82.5% (82.5%) | 1/1 / 0/1 | `16x24` |
| Expensive subj | `point_80` | `glmmplanr` (1,200) | 1/5 | 2,088 [2,080–2,120] | 84.2% (82.3%) | 5/5 / 1/5 | `24x16`(×4), `24x20` |
| Expensive subj | `point_80` | `mlpwr` (1,200) | — | 1,864 [1,520–2,120] | 84.3% (82.3%) | 5/5 / 1/5 | `16x24`(×2), `24x16`(×2), `24x20` |
| Expensive subj | `high_target_86` | `mlpwr` `power=0.86` (1,200) | — | 2,120 [2,120–2,120] | 92.0% (92.0%) | 5/5 / 5/5 | `24x20`(×5) |
| Expensive subj | `wilson_lo_80` | `glmmplanr` (1,200) | **5/5** | 2,120 [2,120–2,120] | 92.0% (92.0%) | 5/5 / 5/5 | `24x20`(×5) |
| Expensive subj | `wilson_lo_80` | `mlpwr` 1k + Wilson 200 (1,200) | **1/5** | 1,632 [1,520–2,080] | 82.5% (82.3%) | 5/5 / 0/5 | `16x24`(×4), `24x16` |

Under the matched 80% point target (`point_80`), `mlpwr`—which incorporates an arbitrary linear `costfun` directly during search—achieved lower mean cost by locating corner designs (`48x8`, `16x24`). However, when an explicit 200-rep 95% Wilson lower-bound confirmation is enforced (`wilson_lo_80`), searching at `power = 0.80` left `mlpwr`'s selected design unconfirmed (`0/5` and `1/5` passing its own 200-rep Wilson gate), whereas `glmmplanr` confirmed passing designs on `5/5` seeds (matching `mlpwr` when its point target was raised to `0.860`).

---

## 4. Rasch Recovery, Fixed-Bank Calibration, and Exploratory Benchmarks

| Script | Summary / Fits CSV | Scope & Key Outcome |
|---|---|---|
| [`rasch.R`](rasch.R), [`rasch_precision.R`](rasch_precision.R) | — | Random-item 1PL recovery (`144x12` vs. `72x24` ability/difficulty RMSE trade-off) and `144x24` confirmation (`62/64` joint attainment). |
| [`rasch_fixed_bank.R`](rasch_fixed_bank.R) | — | 16-item fixed bank with $\pm 3.5$ logit anchors cross-validated against `TAM::tam.mml()` (`64` replications; bank RMSE `64/64` EB vs. `62/64` TAM; worst-item `19/64` vs. `15/64`). |
| [`partial_presentation.R`](partial_presentation.R) | [`partial_presentation_summary.csv`](partial_presentation_summary.csv), [`partial_presentation_fits.csv`](partial_presentation_fits.csv) | Exploratory pilot: `64 x 20` (5 items/person), random intercepts, `30/32` coverage (exact 95% CI `79.2–99.2%`). |
| [`links_covariates.R`](links_covariates.R), [`probit_null.R`](probit_null.R) | [`link_summary.csv`](link_summary.csv), [`probit_null_summary.csv`](probit_null_summary.csv) | External `lme4` logit vs. probit link and correlated covariate sensitivity checks (`959/960` valid; probit null recheck `5/64`–`7/64`). |
| [`extreme_patterns.R`](extreme_patterns.R) | [`extreme_summary.csv`](extreme_summary.csv), [`extreme_frequency.csv`](extreme_frequency.csv), [`extreme_iterations.csv`](extreme_iterations.csv) | All-correct / all-incorrect unit frequency and Armijo-backtracking convergence checks under `chance = 1/3`. |
| [`native_interactions.R`](native_interactions.R), [`trial_covariates.R`](trial_covariates.R), [`covariate_allocation.R`](covariate_allocation.R), [`within_between.R`](within_between.R) | [`native_results.csv`](native_results.csv), [`trial_covariate_summary.csv`](trial_covariate_summary.csv), [`covariate_allocation_summary.csv`](covariate_allocation_summary.csv) | Exploratory interaction, trial-level covariate, and within/between allocation benchmarks. |
