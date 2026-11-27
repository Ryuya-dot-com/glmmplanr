# Statistical Validation & Reproducibility Suite (`glmmplanr` 0.1.0)

This directory archives the simulation scripts, trial-level CSVs, and summary tables supporting the benchmarks summarized in [`inst/validation-guide.md`](../inst/validation-guide.md). It is tracked in Git for full reproducibility alongside the package source, and excluded from CRAN source tarballs via `.Rbuildignore` to keep package downloads compact (the core 3AFC confirmation script and CSVs are also shipped inside `inst/validation/`).

Current development centers on built-in estimation for participant-by-item
sample-size and precision planning. TAM and mirt are the primary IRT
comparators for fixed-bank calibration; comparisons with mfrmr in Sections
4L and 4M are retained as checks of a binary special case of its many-facet
framework. A fixed-bank MML change must be evaluated on its own numerical
path and does not validate crossed-GLMM inference. Historical benchmark
scripts explicitly retain the engines used for their saved results, including
`engine = "lme4"` in `rasch_precision.R` and the random-item check in
`rasch_fixed_bank.R`, following the change to the public random-item default.

The [native covariance correction](#l2-uncertainty-correction) supersedes the
zero-correlation, zero-guessing crossed-GLMM implementation used in earlier
benchmarks, including the 60-design vocabulary study. Those results retain
their original executable and remain historical evidence, not confirmation
of the corrected implementation over the entire old grid. Fixed-bank MML
is byte-identical across this correction; its existing evidence is retained.

<a id="identifiability-and-ranking-repair"></a>

## Identifiability, preference ranking, validation intervals and replay separation (2026-11-27)

The follow-up to `daa4d83` fixes planning, identifiability, worker-version,
reporting and reproducibility issues and separates read-only archive checks from
fresh isolated recomputations. In `src/glmm_laplace.cpp`, `cpp_fit_rasch_fixed_mml`
rejects unidentifiable person-variance data and `fit_glmm_core_ws` scales
above-chance starting-mean bounds (`[chance + 0.05(1 - chance), chance + 0.95(1 - chance)]`)
and numerical boundary epsilons by `1 - chance` so valid high guessing floors
(`chance >= 0.95`) initialize cleanly without changing `chance = 0` or interior
`chance = 1/3` fits; priors, covariance updates and default interval methods are
unchanged.

| Issue | Corrected behavior and check |
| --- | --- |
| Worker build verification & `simulation_version = 5L` | `run_parallel_tasks()` prepends the parent `glmmplanr` library path when starting `mirai` pools and verifies `.glmmplanr_source_md5` across all active daemons (aborting before stamping `simulation_version = 5L` if pre-existing daemons hold a mismatched build). `simulation_version = 5L` ensures `topup_verification()` blocks and `switch_preference()` warns on incompatible pre-v5 saved candidate tables (including `simulation_version = 4` plans affected by positive-`chance` clamping, non-0.80 `ci_assurance` quantiles in `ci_width_q80`, or `n_valid == 0` conditional intervals, alongside pre-v4 unnormalized contrast weights), and `isolated_seeded_call()` extracts global symbols from default argument expressions in `formals(fn)` as well as `body(fn)` under Box-Muller RNG. |
| Named RE alignment, effect conflict guards & schedule condition/balance checks | `calibrate_dgp()` aligns named `subj_sd`, `item_sd`, `testlet_sd`, `slope_sd_ratio`, and named `corr_structure` matrices by RE column or single-df term name (rejecting partial, duplicate, unknown, or mismatched row/col names); `cell_probs_from_diffs()` rejects duplicate, conflicting-alias (`c(C = 30, C1 = 8)`), or `"(Intercept)"` specifications; `simulate_dgp_data()` verifies schedule `allocation`, `between_subj`, `between_item`, `items_per_subject`, trial counts, unit-ID ranges, per-subject distinct-item quotas, distinct within-pair conditions (`(subj_id, item_id, cell_id)`, catching e.g. `A,B -> A,A` pair mutations), whole-testlet selection, and within-unit Latin-square balance against `dgp$design`; and integer dimension/range validators reject negative or fractional inputs before coercion. |
| High-chance initialization, `mafc_binomial()` link precision, `n_valid == 0` conditional intervals, width quantiles & reporting | Above-chance C++ initialization, IEEE-754 machine-epsilon bounds on `(chance, 1)` in `mafc_binomial()` (preserving round-trip `linkfun(linkinv(eta))` at both `chance = 1 - 1e-8` and `eta = -18`), and Stan `real<lower=0, upper=1> chance;` support `chance` across `[0, 1)`. Under `failure_policy = "conditional"`, `n_valid == 0` sets `power_mcse` and conditional Wilson bounds to `NA_real_` (keeping `effective_* = 0` and `fail`); `ci_width_q80` stays fixed at the 80th percentile alongside `ci_width_q_assurance`; `plot(validate_design_point(...))` labels zero-truth `"planned"` rows as `"Type I error"` (`alpha`); `as_gt()` adds `rope_ci` for `criterion = "rope"` plans and aborts on uncomputed `"sensitivity"`/`"robust"` requests; and Rasch plan objects retain `dgp`, `items_per_subject`, `item_seed`, `attrition`, and `participant_dropout`. |
| Missing-condition fast path | Unify fixed-effect cell and rank identifiability checks (`assess_trial_identifiability()`) across `generate_design_schedule()`, cached C++ fits (even when `attrition = 0` and `participant_dropout = 0`) and data-frame fits. Partial-presentation schedules omitting required cells return `insufficient_data` (`0/nsim` valid, `fail`) on both paths. |
| Fixed-bank & random-item Rasch variance identifiability + item names | Fixed-item Rasch MML (`engine = "mml"` or `"tam"`) marks replications invalid whenever fewer than two items are estimable or no examinee has `>= 2` observed responses on estimable (non-constant) items, and random-item Rasch (`plan_rasch_precision()` and `plan_rasch_fixed_bank(method = "random_item")`) requires `>= 2` observed participants, `>= 2` observed items, and `>= 2` observed responses on at least one participant and at least one item, preventing unidentifiable 1-to-1 or 1-item-per-person schedules from passing SD tolerances. `plan_rasch_fixed_bank()` retains input names on `item_difficulties` and adds `item_name` to `item_summary`. |
| Schedule/design compatibility, reserved factor names & custom testlet IDs | `design_crossed()` / `design_testlet()` reserve all internal schedule and formula columns (`subj_id`, `item_id`, `testlet_id`, `cell_name`, `trial_id`, `Z_*_cells`, etc.); `simulate_dgp_data()` verifies schedule cell/level/testlet compatibility against `dgp$design`; and `validate_custom_data()` maps omitted `subj_id`/`item_id` by matching `schedule$trials` labels rather than `as.factor()` re-indexing. |
| Attainable-zero probability contrasts, shared cell-name alignment, saved-plan version checks & sensitivity/reporting sync | `specify_target()`, `null_target_spec()`, Type M/S diagnostics and `prob_target_kind()` evaluate bounds and residuals on unit-normalized weights `w / max(|w|)` (calibrating representable asymmetric nulls such as `c(-1, 1e-13)` at `chance = 0` while distinguishing double-precision underflow boundaries from mathematical unattainability). `align_named_cell_vector()` unifies exact-first cell-name matching and partial/duplicate-name rejection across `calibrate_dgp(cell_probs = ...)` and `align_cell_weights()` (supporting spaces and dots in factor levels, and keeping `cell_contrast_weights` cleanly exported under `roxygen2::roxygenise()`). `simulation_version = 4L` blocks `topup_verification()` under both `"independent"` and `"pooled"` from mixing or pooling incompatible pre-v4 candidate tables; `sensitivity` tracks `source_scenario` across `switch_preference()` and `topup_verification()`; `as_gt()` headers show the limiting scenario's effect; `report_plan()` computes sensitivity ranges across all slope levels including `"standard"`; and `isolated_seeded_call()` propagates `lib.loc`/MD5 checks and `.GlobalEnv` callback bindings under Box-Muller RNG. |
| Surrogate `min_valid_rate` gating, `n_eval_init` guards & exact WLS trend | Under `failure_policy = "conditional"`, `build_surrogate_surface()` fits a monotone GP for `valid_rate` (`pred_valid_rate`), uses observed `valid_rate` at evaluated grid points and Laplace-smoothed thresholds on unevaluated points (so `min_valid_rate >= 0.99` does not reject `100%`-valid surfaces), and `select_stage2_candidates()` retains observed feasible candidates when `pred_criterion_met` is empty rather than defaulting to the grid ceiling. `optimize_sample_size()` requires `n_eval_init >= 1` once post-constraint candidate counts enter the `surrogate_gp` subset branch (while allowing `n_eval_init = 0` for `grid` or auto full-grid confirmation), `fit_nonneg_slope_wls()` solves centered active-set Cholesky analytically rather than calling `L-BFGS-B` (preventing `non-finite value supplied by optim` on 3-point underdetermined designs), and multi-scenario `optimize_sample_size()` resolves `stage1_re` and `stage1_engine` independently per scenario DGP. |
| Validation intervals & post-top-up Wilson reporting | `validate_design_point()` returns `decision_status` and Wilson intervals for `width_attainment`, `joint_attainment` and `rope_attainment` (plotted as `"ROPE exclusion"`). `topup_verification()` preserves `search_meta$mc_conf_level_effective` as the initial Stage-2 level while storing per-candidate `mc_conf_level_used`, so `report_plan()` cites the exact Wilson level of the selected candidate alongside fixed-`item_seed` partial-presentation conditioning notes. |
| Conditional-validity interim curtailment | Under `failure_policy = "conditional"`, `early_stop` curtails for `valid_rate` only when `(n_valid_interim + n_remaining) / nsim < min_valid_rate` (deterministic impossibility), preserving candidates that reach the validity floor over the full run. |
| Preference ordering and scenario trial costs | Make `fewer_subj`, `fewer_items` and `fewer_trials` strictly lexicographic across arbitrary candidate ranges. When multi-scenario designs differ in scheduled or expected trial counts at the same `(n_subj, n_item)`, report scenario-specific counts explicitly rather than attributing one scenario's count to all scenarios. |
| Full-grid multi-scenario deduplication | Skip redundant exploratory searches when `method = "grid"` or `nrow(grid_df) <= n_eval_init + n_eval_adapt`, confirming the fixed candidate-by-scenario grid once on fresh streams. |
| Post-top-up state and absolute-target plots | `topup_verification()` refreshes `evaluated` and `surface` alongside `verified` (single- and multi-scenario, preserved across `saveRDS()`/`readRDS()`). Absolute probability targets keep `pred_power = NA` rather than `0.5` and default `plot()` to `ci_width` (also supporting `metric = "aipe_prob"`). |
| Pre-coercion input checks | Reject factor/character `y` in `validate_custom_data()` before integer conversion, and validate positive whole-integer replication and worker counts before coercion. |
| Observed-data interval workflow & CI | Enforce `model.frame()`-aligned cell-separation checks (robust to dropped `NA` rows), numeric optimizer-convergence checks (`opt == 0`) and finite-vcov checks on both `glmer` and `bootMer` refits in `vignettes/planning-workflows.Rmd`, verified by executing the vignette chunk in `test-review-regressions.R`; add `release/**`, `workflow_dispatch`, `R 4.1` (declared minimum) and `R-devel` to `.github/workflows/R-CMD-check.yaml`. |

### Separating archived-result verification from fresh recomputation

1. **Read-only audit of saved replay CSVs (no re-simulation or file modification):**
   ```sh
   python3 validation/verify_replay_against_archive.py
   ```
   `verify_replay_against_archive.py` checks header schemas, internal marginal
   and logit interval/detection consistency, and all marginal, conditional
   (`logit`, `logit_se`, `logit_truth`, `logit_covered`) and profile-diagnostic
   columns without Python `assert` statements, so `python3 -O` fails non-zero
   identically on any `NaN`/`Inf`, tolerance, or decision mismatch.
2. **Fresh recomputation into a new empty directory (no block reuse):**
   ```sh
   # Isolated install + fresh replay into a temporary directory (~2 MB staged inputs, not 213 MB)
   sh validation/run-isolated.sh covariance_mode_replay
   # Verify the newly computed directory against the archived baseline:
   python3 validation/verify_replay_against_archive.py /path/to/glmmplanr-validation.XXXXXXXX/replay_out
   ```
   Passing `--fresh <new_output_dir>` to `validation/covariance_mode_replay.R`
   refuses non-empty directories so cached `.rds` blocks are never silently
   reused. When resuming an interrupted run without `--fresh`,
   `checkpoint_dir()` hashes all `R/*.R` and `src/*.{cpp,h}` package source
   files alongside script and input MD5s (including `cases.csv` and
   `diagnostics.csv`), `settings` (including `cases`), R version and package
   versions, verifies per-block `case` metadata on resume, and verifies that
   the loaded `glmmplanr` namespace's build-time `.glmmplanr_source_md5`
   matches the current working-tree sources (rejecting stale or mismatched
   `R_LIBS` installations while still permitting separate temporary installs
   built from the identical source).

3. **Input dependencies staged by `validation/run-isolated.sh`:**
   By default `run-isolated.sh` copies package sources plus `validation/*.R`
   and parses `validation/$name.R` together with transitively `source()`d
   helper scripts to stage every existing `validation/...` input file or
   directory automatically (set `VALIDATION_COPY_ALL=1` only if a full copy of
   all ~213 MB archived outputs is needed):
   - `rasch_fixed_bank`, `coverage_recheck`, `stan_uncertainty`,
     `three_choice_precision`, `package_calibration`: package source +
     `validation/*.R` only (no archived subdirectories required).
   - Scripts sourcing `planning_comparison_helpers.R`: stages
     `validation/item_pool_planning_verified/pool.csv` and `forms.csv`, plus
     script-specific inputs (e.g., `validation/item_pool_search/search.csv` for
     `confirm_item_pool_search` and `compare_item_pool_irtsim`).
   - `covariance_mode_replay` / `inspect_covariance_profile_failure` /
     `check_covariance_scores`: stages `validation/item_pool_planning_verified/`
     and `validation/covariance_review_medium/`.
   - `interval_mechanism_audit` / `summarize_l2_default_bootstrap`: stages
     `validation/l2_default_bootstrap_confirm/`.
   - `summarize_l2_vocabulary_extremes` / `summarize_few_shared_passages` /
     `check_two_pl_*`: stages matching stem-prefixed scenario/summary/fit
     files (e.g., `validation/l2_vocabulary_extremes_*`) constructed via
     `file.path()` / `paste0()`.
   - `three_choice_cov_decomposition`: stages
     `validation/three_choice_precision_native.csv` and
     `validation/three_choice_precision_stan.csv`.

<a id="planning-contract-repair"></a>

## Planning decisions, custom data and numerical failure contracts (2026-11-26)

The follow-up to `e5149b8` fixes six bounded implementation/reporting issues;
it does not change the response-model roadmap, priors, interval methods or
nominal-inference recommendations.

| Issue | Corrected behavior and check |
| --- | --- |
| Zero random slopes | Skip the three sensitivity DGPs when every participant/item/testlet slope SD is zero. Positive slopes in any block remain eligible. |
| Late MC adjustment | Apply Bonferroni-adjusted decisions before step-up, ceiling and top-up; reclassify earlier candidates when the family grows. Deterministic rates straddling the Wilson boundary exercise all three branches. Scenario candidate selection precedes fresh confirmation with the candidate-by-scenario multiplicity fixed; exploratory searches do not run top-ups. |
| Stale top-up metadata | Refresh `nsim_topup_used` and `topup_only_power` on subsequent independent or pooled top-up calls. An independent batch replaces the inference batch rather than accumulating its denominator. |
| Custom testlet responses | Validate finite positive integer IDs and bounds; derive omitted passage IDs from scheduled items and reject conflicting passage labels/IDs. Exercise the documented four-column callback and retained passage IDs after subsetting. |
| Absolute probabilities and printing | Reject probability-scale power/both/joint targets without weights of both signs. Undefined power/joint fields are NA; precision, coverage and reports remain available. Joint-attainment console text is emitted as a separate bullet. |
| Numerical failure guards | Reset working likelihood on polishing; track nonpositive observed curvature separately from iteration failure; use observed curvature for polishing gradients. Derivative failures return a failed fit rather than throwing. Fixed-bank MML validates bounds, duplicate person-item rows (including zeros) and nonconstant items marked estimable before boundary calculations. |

The new regression cases are in
[`test-review-regressions.R`](../tests/testthat/test-review-regressions.R),
alongside the existing curvature/saddle and Rasch precision tests. Mocked
Monte Carlo counts test the optimizer's control flow, not statistical power
or exact familywise coverage of adaptively selected Wilson intervals. The
fresh scenario confirmation can require more fits than earlier result reuse;
removing redundant zero-slope sensitivity eliminates a separate unnecessary
cost without changing the DGP.

For the numerical kernel, reuse the preceding 2,048 stored responses and 84
selected profiles with unchanged estimands, starting values, fit settings,
random streams, failure denominators and unrounded decision rules:

```sh
# Use R_LIBS pointing to an isolated installation of this source.
Rscript --vanilla validation/covariance_mode_replay.R validation/review_contract_replay
python3 validation/review_contract_replay/verify.py
```

The replay script retains comparisons against the original medium run;
`verify.py` additionally compares directly against `e5149b8`'s
`curvature_contract_replay/`. All 2,048 native final/coarse fits and 84 profiles
remain valid, with zero failed profile objective evaluations. Point estimates,
SEs, response hashes and all four interval/decision indicators are identical
to the preceding repair (maximum point/SE difference 0pp; the prior 1e-4pp
tolerance is not widened). See
[`verification.json`](review_contract_replay/verification.json). The three
older pre-roundoff tolerance exceptions remain historical evidence.

An additional bounded four-person/two-item probe enumerated 254 nonconstant
response patterns at three fixed slope variances (762 cases): both versions
returned 666 converged fits and 96 curvature failures. A 2,001-point variance
scan for the existing saddle response compared joint and conditional
uncertainty in the old version: both returned 1,840 converged and 161 failed
fits, without a derivative exception. Thus these probes do not reproduce the
specific perturbation-only exception or establish its frequency. Its removal
is checked at the implementation level: every derivative caller propagates
the failure flag into `non_pd_hessian`/missing likelihood. The saved-response
replay and existing negative-curvature controls assess regression behavior;
they are not a new nominal-coverage study or a global-mode guarantee.

Win-builder, GitHub CI and release packaging remain deferred until the feature
set is finalized. No new API or independent large simulation is introduced.
The isolated optimized macOS/R 4.6.1 installation passes one full invocation
with `NOT_CRAN=true`: 824 assertions, zero failures/warnings and one opt-in
Stan skip. The native source and compiler options are unchanged across the
replay build and the subsequent R-documentation installation; their binary
hashes are recorded separately rather than claiming binary identity.

<a id="conditional-curvature-contract"></a>

## Rejecting conditional saddles in native fits (2026-11-26)

A review found a failure of the final curvature contract, separate from the
roundoff repair: a small Fisher step could set both convergence and observed
curvature readiness after negative observed weights had been clipped. The
four-person, two-item, 3AFC counterexample has one all-correct person and three
all-incorrect people, condition coding -1/+1, fixed participant intercept/slope
SDs .1/4 and item intercept SD .1. Its standardized conditional negative-log
density Hessian has minimum eigenvalue -1.081659; moving .01 along that
direction increases the log density by 5.408e-5. The old fit reported
convergence and logLik -6.020682 despite this saddle. This deterministic
counterexample does not estimate a failure rate under a planning DGP.

The final check now requires a small observed-Newton step and unmodified
positive-definite observed curvature. A small Fisher step alone is rejected,
and derivative evaluations no longer substitute clipped observed curvature.
The counterexample returns `non_pd_hessian`, `converged=FALSE` and `logLik=NA`,
with and without joint covariance propagation. Regression tests also exchange
the participant/item blocks to exercise automatic orientation, and retain a
positive-curvature control. The old code fails the 12 rejection assertions;
the repaired code passes them. The unchanged 2e-7 tolerance, priors, covariance
updates and iteration budgets remain in use. This is a validity guard, not a
new saddle-escape algorithm or a global-mode/coverage guarantee.

Reproduction: install the source in an isolated library, set `R_LIBS` to that
library, then run
`Rscript --vanilla validation/covariance_mode_replay.R validation/curvature_contract_replay`.
The script retains the original saved-response comparison; its optional output
directory prevents overwriting the preceding repair. The additional
[`verify.py`](curvature_contract_replay/verify.py) compares the new CSVs directly
with `covariance_mode_replay/`, not with the older, pre-roundoff executable.
The same predeclared 1e-4pp point/SE tolerance and unrounded decisions apply.

All 2,048 native final/coarse fits remain valid, with identical stored point
estimates, SEs, response hashes and coverage/detection/width/joint decisions
relative to the immediately preceding repair. All 84 selected profiles remain
valid, with zero failed objective evaluations and identical point/SE and
interval decisions. Maximum point and SE differences are zero. No new DGP
responses, bootstrap draws or independent nominal-inference confirmation were
generated. See [`verification.json`](curvature_contract_replay/verification.json)
and the [dense-curvature reproduction](curvature_contract_replay/counterexample.R).
Historical raw results and the three earlier tolerance exceptions are retained.

Local validation used an isolated optimized installation on macOS/R 4.6.1.
The curvature/status groups passed 94 assertions; the other test groups passed
626, and the two CRAN-skipped cases passed 54 in a separate run. Their 34
previously completed assertions are not counted twice: 740 distinct assertions
pass, with zero failures/warnings and only the opt-in Stan check unexecuted.
This is combined test evidence, not one uninterrupted invocation. An initial
test-selection command failed before running tests because of string escaping;
the selector was corrected without changing tests or restarting completed groups.
The updated validation vignette renders successfully and contains the current
guard and statistical distinctions. Package installation, the tests and the
render are recorded in [`local-validation.json`](curvature_contract_replay/local-validation.json).
A new source archive or full `R CMD check` is not claimed.

Printed optimization results now include the same analysis description as
specifications and plan reports. The maintainer has deferred Win-builder and
GitHub CI until the feature set is finalized; no push or external submission
is made for this correction. This repository-only validation directory is not
bundled in CRAN archives. The current sections are not yet on published `main`;
their publication and the guide's external links must be verified at release.

<a id="interval-mechanism-audit"></a>

## Saved-draw audit of interval calibration and planning thresholds (2026-11-26)

This exploratory follow-up asks why the 72×60 planned joint decision changes
between joint-delta and default bootstrap intervals, and why the 162×30 null
rejection remains high. It uses the completed default-bootstrap study's
three frozen designs × null/planned × 512 datasets, with 99 saved parametric
refits for every dataset and 399 for the first 32 per cell. The estimand is
the same marginal effect in percentage points (null 0, planned 8); schedules,
data, fitted estimates, priors and simulation seeds are unchanged. The
54×30 intercept-only design is a control. These are the covariance-corrected
fits before the later mode-roundoff repair; source identities are retained,
and no new current-kernel fitting or independent confirmation is implied.

The data are simulated binary responses under the frozen vocabulary-study
planning assumptions, not raw participant data. Each person sees 30 items
in a Latin-square schedule (`item_seed=37`); guessing and random-effect
correlations are zero. Planned marginal probabilities are .42/.50. The
participant/item intercept SDs are .8/.6 for 54×30 and 72×60, and 1.1/.9
for 162×30; corresponding slope SDs are 0/0, .25/.35 and .4/.5. The source
study's native regularized estimator, priors and equal-marginal null are
retained exactly.

Before reanalysis, the following comparisons were fixed. First, inspect
bias, empirical SD versus RMS SE, standardized-error tails, and the joint
detection/width decision's two components. Mean-bias removal and replacing
SE by the same cell's empirical SD are descriptive oracle diagnostics, not
available interval corrections. Second, compare the existing Wald interval
with percentile and basic intervals from exactly the same parametric draws,
using boot's normal-quantile interpolation. No bootstrap-t or BCa calculation
is attempted because per-refit variance/influence information is not saved.
Both alternative intervals are exploratory, not new public methods. The
bootstrap Wald interval remains centered at the original estimate; it is
not boot.ci's bias-adjusted normal interval.

Third, on only the 32 datasets per cell with 399 draws, inspect four disjoint
99-draw blocks (three unused draws), and 2,000 nested empirical resamples
of size 99/399 from each saved 399-draw pool (with replacement; root
953260001). This isolates variability of the interval computation conditional
on the finite pool. The 399-draw pool is not an infinite-bootstrap truth.
The 2,000 resamples are not additional independent outer datasets; sampling
uncertainty across the 32 outer datasets remains distinct from conditional
resampling error. No selection or extension of outer replications is made.

Valid intervals require valid original fits, >=90% successful inner refits,
and finite ordered endpoints; failures remain in detection, width and joint
denominators. Coverage uses valid intervals. Report pointwise Wilson 95%
MC intervals, paired differences/MCSEs, and the existing width <=10pp,
joint lower bound >=.80, coverage lower bound >=.925 and null upper bound
<=.075 rules. These practical gates are distinct from nominal .95/.05.
The saved original Wald SEs/counts must match within 1e-7pp and exactly,
respectively. Published claims and defaults will not change on exploratory
improvements alone. See the [boot interval documentation](https://stat.ethz.ch/R-manual/R-devel/library/boot/html/boot.ci.html)
for interval definitions and required inputs; the present generating scheme
remains parametric and uses the native fitted GLMM.

**Completed results and interpretation.** No new response datasets or GLMM
fits were generated. All 3,072 original datasets and 361,728 saved refits
were retained. The 12,864 computed intervals include four methods for all
512 datasets per cell and three 399-draw methods for only 32 per cell.
Original Wald SEs, valid counts, decisions and paired response hashes match.

At moderate-slope 72×60, joint attainment falls from 434/512 with joint delta
to 395/512 with the default Wald bootstrap: 58 losses and 19 gains. Of the
58 losses, 54 retain detection but fail width <=10pp, and four lose detection
while retaining the width criterion. Detection actually rises from 484 to
486; width attainment falls from 457 to 415. Thus lower joint attainment
in this batch primarily concerns the width threshold, not lower power.
Similar RMS SEs do not imply similar per-dataset threshold decisions.

The long-draw subset contains 32 datasets per cell, not 512. At 72×60,
the four disjoint 99-draw batches give planned joint counts 22, 25, 25 and
23 out of the same 32 datasets, versus 26/32 using all 399. In the conditional
finite-pool experiment, mean SE CV falls from 7.15% to 3.54% as B rises from
99 to 399. Mean conditional joint probabilities are 74.86% and 78.74%, a
3.88pp difference (conditional resampling MCSE .134pp; across-outer SE of
the difference 1.414pp). The latter is not the uncertainty of a fresh
512-dataset confirmation. The empirical pools are themselves noisy and
shared between the comparisons. These results establish sensitivity to
inner simulation; they neither attribute the entire full-sample -7.62pp
difference to B nor establish that B=399 passes the 80% planning gate.

For the high-variance 162×30 null, the mean error is .01099pp (MCSE .09284),
the empirical SD is 2.10072pp, and bootstrap RMS SE is 2.10517pp (ratio 1.00212).
The standardized-error SD is nevertheless 1.04976, with 18 left-tail and 22
right-tail exclusions. The lowest-SE fifth has RMS SE 1.68590pp, error RMSE
2.25845pp and 17/102 exclusions; the highest fifth has RMS SE 2.51965pp,
error RMSE 2.20089pp and 6/103 exclusions. These are post hoc descriptive
strata. A near-correct aggregate SE does not demonstrate that each fitted
SE adequately reflects that dataset's error. Removing the observed mean
bias gives 41/512 exclusions, while substituting the cell's empirical SD
gives 31/512. Neither oracle operation is available for an actual new study,
and neither identifies a unique cause or a universally correct multiplier
or degrees of freedom. In particular, this is not evidence that finite-B
noise alone explains the 7.81% null rejection.

The unmodified joint-delta null rejection is 35/512 (6.84%, pointwise 95%
Wilson interval 4.96–9.36%). The oracle 31/512 (6.05%) above must not be
reported as joint-delta performance. Likewise, 434/512 versus 395/512 compares
joint delta with B=99, not B=399 with B=99; the width limit is 10pp and the
planned effect is 8pp. The experimental model has one two-level factor.

The alternative intervals use the same parametric refits and original point
estimate; both remain exploratory. On the full 512-dataset cells:

| Metric | Default Wald (99) | Percentile (99) | Basic (99) |
| --- | --- | --- | --- |
| 162×30 null rejection | 40/512 (7.81%) | 32/512 (6.25%) | 36/512 (7.03%) |
| 162×30 planned coverage | 480/512 (93.75%) | 485/512 (94.73%) | 488/512 (95.31%) |
| 72×60 planned coverage | 489/512 (95.51%) | 489/512 (95.51%) | 494/512 (96.48%) |
| 72×60 detection AND width <=10pp | 395/512 (77.15%) | 318/512 (62.11%) | 317/512 (61.91%) |

The percentile null interval at 162×30 is 4.46–8.69%, still indeterminate
against the practical 7.5% upper gate; its paired improvement is not a
general calibration guarantee. At 72×60 its mean width is 9.60pp versus
9.20pp for Wald, and its per-dataset widths also vary more at B=99. The
399-draw subset is insufficient to rank nominal coverage. Selecting the
method with the largest planning pass rate would ignore inferential validity;
selecting percentile/basic solely for a smaller observed null rate would
ignore both uncertainty and precision cost. No new interval option, default
change, t correction, or general participant/item recommendation is adopted.
Independent confirmation follows only a sufficiently justified candidate;
these results do not meet that adoption criterion.

[Analysis](interval_mechanism_audit.R) ·
[independent arithmetic verification](verify_interval_mechanism_audit.py) ·
[all interval summaries and Wilson bounds](interval_mechanism_audit/interval_summary.csv) ·
[paired changes and MCSEs](interval_mechanism_audit/paired.csv) ·
[error/SE diagnostics](interval_mechanism_audit/error_diagnostics.csv) ·
[SE strata](interval_mechanism_audit/se_strata.csv) ·
[finite-pool sensitivity](interval_mechanism_audit/finite_pool_summary.csv) ·
[figure](interval_mechanism_audit/interval_comparison.png).

To reproduce with an installed glmmplanr and boot, use a fresh output directory
so the retained evidence is not overwritten:

```sh
Rscript --vanilla validation/interval_mechanism_audit.R /tmp/glmmplanr-interval-audit
python3 validation/verify_interval_mechanism_audit.py /tmp/glmmplanr-interval-audit
```

A fresh-directory execution of these entry points reproduced all nine numeric
tables and the PNG byte-for-byte (comparing decompressed CSV content for gzip
files), and passed the independent verification again. The original response
and refit files were not modified; see [reproduction record](interval_mechanism_audit/reproduction.json).

R assertions and an independent Python implementation verify all 42 summary
rows, 36 paired comparisons, draw-based interval endpoints, counts, Wilson
bounds and error moments. The largest endpoint discrepancy is 1.51e-13pp
(tolerance 1e-10pp). The figure was visually inspected: theme_bw, no automatic
title/caption. Two initial analysis-only attempts stopped at key ordering and
CSV empty-string type checks; those bookkeeping issues were corrected before
the successful analysis, without altering input data or the model. Logs and
[verification details](interval_mechanism_audit/verification.json) are retained.
The package kernel/API and tests are unchanged; successful 740-assertion
evidence is reused rather than rerun. This study does not extend the old
candidate's platform-check coverage or replace its existing limitations.

The guide update is included in a documentation-only rebuild (1,258,394 bytes;
SHA-256 `add111e217880877426f11e971d65389b60800b377f5b981c6165a0da05127d8`).
R CMD build and vignette reconstruction passed. Complete archive comparison
against `ca751f6a...` finds changes only in the guide, its HTML, the build
timestamp and three printed environment addresses in another vignette.
All 40 R/C++/Stan/test files are unchanged, and all 13 README relative links
resolve. This combines the previous numerical checks with the new document
build and identity audit, rather than claiming a new full check. See
[archive comparison](interval_mechanism_audit/archive-identity.json) and
[build log](interval_mechanism_audit/build.log).

<a id="covariance-mode-repair"></a>

## Conditional-mode roundoff and diagnostic contract repair (2026-11-26)

The medium study's representative failed perturbation is frozen as a small
synthetic regression fixture. The old installed source reproduces
`inner_not_converged`, step 2.886e-7 and a finite, unusable likelihood.
The repair permits a conditional-mode Newton step only within summation
roundoff when its predicted improvement is itself below that scale and the
unscaled step is small. The next iteration must meet the unchanged 2e-7
mode criterion. Iteration budgets, priors and covariance update formulas
are unchanged. Uncomputed covariance scores now return NA; the likelihood
is NA unless strict convergence and final observed curvature are available.
Specifications and reports explicitly describe the approximate moment path.

Before replay, the point-estimate and SE comparison tolerances were fixed
at 1e-4 percentage points; every coverage, detection and width decision is
compared before rounding, even within those tolerances.
[`covariance_mode_replay.R`](covariance_mode_replay.R) reuses all 2,048 saved
response vectors for native final/coarse fits. It reruns the unchanged archived
profile optimizer only for its 69 failures, the 15 initially unevaluable
scores and the first dataset of each cell (overlaps removed). The remainder
is not recomputed as another full profile study. New and old valid fits,
failures, response identity, likelihood/mode checks and changes in design
decisions remain separate. This replay assesses a numerical repair on saved
data; it is not an independent Monte Carlo confirmation of inference.

The modified conditional-mode routine also participates in fixed-effect
curvature calculations, so the full regression suite is covered in addition
to the targeted regression and same-response replay. Fixed-bank MML code is
unchanged; no new Rasch performance simulation is required. Historical
platform checks remain evidence for their original source.

**Replay results.** Both final and coarse native settings retained 2,048/2,048
valid fits, with zero changes to coverage, detection, width or joint decisions.
Maximum absolute point/SE changes were .000221762/.000380149pp. The
prespecified 1e-4pp comparison tolerance was exceeded in three cases: cell 2,
replications 32, 85 and 118. It was not increased after seeing those results.
All three have zero guessing/correlation. Their covariance iteration counts
changed 33→36, 22→26 and 28→23, respectively; each meets the .005 score rule.
Independent dense likelihood differences were zero at printed precision,
and mode checks passed. These small, inspected numerical differences are
not strict numerical equivalence. The largest covariance-element change
was .000381501. All results and exceptions remain in the denominator.

The selected profile replay comprises 84 unique datasets: all 69 failures,
all 15 initially unavailable scores (already among those failures), and 15
additional first-case controls. All 84 profile fits passed the unchanged
acceptance rules; no failed objective evaluations occurred. The largest
scaled profile score was .0008411. The original 69 failures are preserved in
their old run, not overwritten. In the 16 first cases, dense likelihood
differences were at most 7.4e-13 and mode scores at most 1.83e-7; halving
curvature steps changed SE by at most 2.97e-6pp, and alternate-start objective
differences were at most 1.27e-7. All 15 previously missing initial scores
are now evaluable and above .01, so recovering convergence does not remove
the separate moment-update nonstationarity finding.

The frozen failure and diagnostic semantics passed the initial 71 targeted
assertions. [`verify_covariance_mode_replay.R`](verify_covariance_mode_replay.R)
checks the saved replay, its numerical exceptions and selected reference
fits without starting a new simulation. Results:
[native comparisons](covariance_mode_replay/comparison.csv),
[reference replay](covariance_mode_replay/reference_summary.csv),
[tolerance exceptions](covariance_mode_replay/tolerance_exceedances.csv),
[verification log](covariance_mode_replay/verification.log) and
[execution log](covariance_mode_replay/run.log).

**Regression and packaging.** All test groups are covered by 740 passing
assertions, zero failures/warnings and one opt-in Stan skip. The source
`test_local()` invocation used pkgload's default debug build (-O0): after
34 calibrate assertions passed, it was interrupted during core-workflows
and the remaining groups ran against the matching optimized installation
(706 assertions, 294.8 seconds). A single frozen-response timing probe gave
.030/.026 seconds for old/new optimized builds and 1.602 seconds for the
debug build, each with 9 iterations; this is not a general speed benchmark.
The incomplete invocation is retained, and no claim of a single uninterrupted
full run is made. See [execution settings](covariance_mode_replay/regression-environment.txt),
[interrupted log](covariance_mode_replay/regression-debug-interrupted.log)
and [completed remaining groups](covariance_mode_replay/regression-installed.log).

The new 1,257,141-byte archive has SHA-256
`ca751f6a0405f8f053eee6f57976dbd7b37647fd1412efb254f566da7d528355`.
On macOS/R 4.6.1, `R CMD check --as-cran --no-manual --no-tests` passed with
0 errors, 0 warnings and 1 NOTE (new submission; optional cmdstanr repository,
availability confirmed). Tests use the separate evidence above rather than
duplicating the full suite. Installation, compiled code, examples, help,
URLs and vignette reconstruction passed. Optional instantiate was unavailable.
The initial sandboxed check stopped at DNS resolution before installation;
the network-enabled check completed. The archive's 43 code/test/selected-doc
files match the source, the optimized installation's 21 R/C++/Stan source
files match, and all 13 relative README links resolve. See
[archive identity](covariance_mode_replay/archive-identity.json) and
[package check](covariance_mode_replay/00check.log).
No new Win-builder upload, three-platform CI or CRAN submission was performed;
the prior archive's platform results are not assigned to this changed kernel.

<a id="covariance-review-medium"></a>

## Paired covariance-path study: prespecified medium run (2026-11-26)

The authorized follow-up asks whether the review's nonstationarity finding
affects point estimates, intervals, and planning decisions. It does not
change the package estimator or establish a new public inference method.
[`covariance_review_medium.R`](covariance_review_medium.R) evaluates 16 cells:
complete 32 participants by 24 items, or partial 72 by 60 with 30 per person;
guessing floor 0 or 1/3; random-effect correlation 0 or .2; and marginal
effect 0 or 10pp with baseline .55. Participant intercept/slope SDs are
.8/.35 and item SDs .6/.4. Latin-square schedules use item seed 37 and are
fixed across replications. Layout and sample size change together, so this
study does not identify a pure effect of partial presentation.

Each cell has 128 independent outer replications (2,048 datasets total),
using seed root 952260001 and separate case/replication streams. A separate
one-per-cell technical run uses root 951260001 and is excluded from rates.
Each pair uses identical responses. The current regularized native fit uses
final-stage budgets and joint-delta SEs; bootstrap intervals are outside this
comparison. Coarse-stage convergence is recorded on the same responses.
Failure counts remain in the scheduled denominator for detection and joint
attainment; coverage conditional on valid intervals and all-trial coverage
are both reported. Joint attainment means zero excluded and a 95% interval
width at most 20pp. Pointwise Wilson 95% bounds versus .80 give
pass/indeterminate/fail, separately from coverage and null rejection.
At 128 replications, MCSE is 1.93pp near 5% and 3.54pp near 80%; small
departures cannot be resolved reliably.

The diagnostic comparator profiles beta and conditional modes with the same
native likelihood evaluator, then numerically optimizes four log SDs and,
for correlated fits, two Fisher-transformed correlations. It retains the
native beta/Half-t priors, log-SD Jacobian, and LKJ density on the correlation
matrices. Fisher coordinates are an optimization device, not an additional
prior Jacobian. This comparator is not independent of the native likelihood,
not HMC, and not a demonstrated global optimizer. An independent dense
likelihood/mode-score calculation checks the first dataset in each cell.
An alternate optimizer start and halved curvature step are also checked
on those first datasets. Bounds are SD exp(-5) to exp(2), correlation
[-.98,.98]; boundary optima within 1e-4 are not accepted as interior reference
fits. No eigenvalue repair is used for the full profile Hessian: nonpositive
curvature, failed perturbations, optimizer failure, scaled score >= .01,
or an objective increase > 1e-5 produce an invalid diagnostic interval.

Reference uncertainty combines conditional beta uncertainty with the full
profile covariance of log SDs and correlations. It changes both the optimum
and uncertainty propagation, so any coverage change is not attributable to
one component alone. Paired point changes, SE ratios, individual decision
changes, bias, empirical SD, coverage/null rates and MC uncertainty are
retained. Likelihood check tolerance is 1e-6; first-case curvature SE changes
above .05pp or alternate-start objective differences above 1e-4 require
investigation before accepting the technical setup. Responses, estimates,
curvatures, failures, settings, source/binary hashes and environment are
checkpointed every 16 datasets. Resume requires the same manifest.

Run `GLMMPLANR_REVIEW_MODE=smoke sh validation/run-isolated.sh covariance_review_medium`
first. After checking the technical output, reuse that isolated installation
and copied source with `GLMMPLANR_REVIEW_MODE=medium Rscript --vanilla
validation/covariance_review_medium.R`. The main-run settings and source
must remain fixed after launch. This is a numerical validation study, not
a comparison of newly selected optimal sample sizes or a change to the
manuscript's scope.

### Completed results and interpretation

The main run completed all 2,048 datasets in approximately 290 seconds with
six local workers (macOS, R 4.6.1), excluding installation, the separate smoke
run and subsequent checks. It retained 389,759 profile-evaluator calls,
including 930 unsuccessful evaluations; this is not a benchmark of only
4,096 individual GLMM calls. Native final-stage and coarse-stage fits were
both valid for 2,048/2,048 datasets. The diagnostic profile procedure yielded
1,979/2,048 valid intervals (96.63%). All 69 invalid references occurred with
nonzero guessing, with cell-specific validity from 104/128 to 128/128.
All scheduled datasets remain in detection/attainment denominators.

The central numerical finding replicates beyond the four initial examples.
For zero guessing and zero correlation, all 512 numerical log-SD scores were
below .005 (maximum .004915). In the other paths, all 1,521 evaluable scores
exceeded .01; cell medians were .353 to .693. Fifteen score evaluations were
unavailable because a perturbation did not converge. They are **not** zero
scores or arbitrarily large gradients. The original optimization runner used
a finite failure penalty; that penalty contaminated some recorded finite
differences. [`check_covariance_scores.R`](check_covariance_scores.R) replays
only these initial perturbations, preserves missing scores as NA, and
reproduces the valid scores within 1e-6. Use
[`score_checks.csv`](covariance_review_medium/score_checks.csv) and the
corrected diagnostic summary for score interpretation; the original raw
diagnostics are retained for provenance.

Among the 1,979 common valid pairs, the change in the marginal condition
contrast was small in this DGP. For noncanonical paths, cell-specific RMS
changes were .032 to .120pp, the largest individual change was .538pp, and
the RMS-SE ratio was .985 to 1.004. For the corrected zero/zero path, the
largest point change was .000645pp. Eight common valid pairs changed coverage
classification, six changed detection, and eight changed joint attainment.
These counts do not include invalid diagnostic fits as disagreements about
an available interval. The all-trial differences are separately retained.
Changing the optimum and full covariance propagation did not consistently
improve empirical coverage. Correlation/SD recovery is in
[`parameter_summary.csv`](covariance_review_medium/parameter_summary.csv);
the profile optimizer does not remove small-sample regularization bias.

For the native method, marginal 95% coverage ranged from 93.75% to 98.44%
over 16 cells. Null rejection ranged from 4/128 to 8/128 (3.13% to 6.25%);
all eight pointwise Wilson intervals include .05. This resolution does not
establish exact nominal calibration. Bias ranged from -.405 to +.724pp;
the largest positive value is the zero-correlation, 3AFC, complete-layout
null cell, with MCSE .337pp. Planned-effect bias ranged from -.363 to +.317pp.

| Layout | Guessing | Correlation | Native joint count / 128 | 95% MC interval | Planning decision |
| --- | ---: | ---: | ---: | --- | --- |
| Complete 32 x 24 | 0 | 0 | 86 | 58.66–74.72% | Fail |
| Complete 32 x 24 | 1/3 | 0 | 90 | 61.90–77.54% | Fail |
| Complete 32 x 24 | 0 | .2 | 90 | 61.90–77.54% | Fail |
| Complete 32 x 24 | 1/3 | .2 | 100 | 70.20–84.41% | Indeterminate |
| Partial 72 x 60 | 0 | 0 | 126 | 94.48–99.57% | Pass |
| Partial 72 x 60 | 1/3 | 0 | 127 | 95.71–99.86% | Pass |
| Partial 72 x 60 | 0 | .2 | 126 | 94.48–99.57% | Pass |
| Partial 72 x 60 | 1/3 | .2 | 127 | 95.71–99.86% | Pass |

This table concerns a true 10pp contrast and width <=20pp. The profile
procedure, counting its failures as unmet, gives the same eight categories.
These are two evaluated layouts, not an optimized required sample size;
different total observations and cluster counts prevent attributing their
difference to partial presentation alone. The generous 20pp width criterion
also must not be transferred to the earlier 10pp-width L2 study.

**Reference limitations and verification.** Of the 69 rejected profile
intervals, 34 had optimizer code 52; 35 had code 0 but failed the additional
curvature/score checks. Every one encountered an unsuccessful evaluation.
Huge negative raw Hessian eigenvalues can reflect the failure penalty, not
negative curvature of the statistical likelihood. The first dataset in
case 12 failed the reference checks, including alternate-start agreement;
its alternate objective was .348 lower. Its rejected perturbation reports
`inner_not_converged`, step 2.886e-7 versus the native 2e-7 threshold. Repeating
with its returned beta did not resolve it, and its returned likelihood
differed by -.780 from the independent observed-curvature calculation.
Thus accepting that failed return merely because its step is small is not
justified. This representative investigation does not assign the same cause
to every failed evaluation. No failed result was replaced or silently retried
into the main performance table.

All 16 first-dataset native fits matched the current planning core to 1e-8pp.
Dense likelihood checks at the saved native/profile points differed by at
most 7.4e-13, with maximum conditional-mode score 1.83e-7. The 15 valid
first-dataset references passed the halved-step/alternate-start checks; the
invalid case remains excluded. An initial checking assertion requiring all
16 references to succeed was inappropriate for a study that retains failures;
both the initial and corrected check logs are saved. The failed-perturbation
inspection likewise retains its initial log rather than presenting that
unconverged likelihood as a successful equality check. Independent R and
Python aggregation agreed on all 4,096 method rows, response pairs, counts,
width/coverage/detection decisions and Wilson bounds. No full regression,
Rasch simulation, package rebuild or platform recheck was needed because
the numerical implementation is unchanged. The updated validation vignette
was rendered separately and its new section/results checked in HTML; this
does not update or recheck the existing release archive.

Results and reproducibility: [all trials](covariance_review_medium/fits.csv.gz),
[performance summary](covariance_review_medium/summary_summary.csv),
[paired changes](covariance_review_medium/paired_summary.csv),
[diagnostics](covariance_review_medium/diagnostic_summary.csv),
[verification](covariance_review_medium/verification.log),
[failure inspection](covariance_review_medium/failure-inspection.log), and
[run log](covariance_review_medium/run.log). The RDS batches retain responses,
specifications, schedules, fits and curvature matrices. Summary reproduction
uses [`summarize_covariance_review.R`](summarize_covariance_review.R);
[`check_covariance_review.R`](check_covariance_review.R) and
[`inspect_covariance_profile_failure.R`](inspect_covariance_profile_failure.R)
check the first cases and the declared failure, respectively. Invoke each
script with the results directory, using the recorded isolated library.

The [coverage figure](covariance_review_medium/coverage.png)
([PDF](covariance_review_medium/coverage.pdf)) uses theme_bw without a title
or caption. Bars are pointwise MC intervals for coverage among valid
intervals; the profile denominator can differ. The conclusion is limited:
nonstationarity is reproducible, but its impact on this marginal contrast
was modest among comparable fits. The diagnostic optimizer is not yet a
production replacement or an independent gold standard. Resolve the native
objective/diagnostic contract and strengthen the reference evaluator before
changing defaults or treating numerical stationarity as an accuracy guarantee.

<a id="review-followup"></a>

## Technical review follow-up (2026-11-26)

This follow-up checks the supplied review against the implementation and saved
results. It is not a new exhaustive audit or a Monte Carlo performance study.
The numerical source is commit `50f2c4d`; `src/glmm_laplace.cpp` has MD5
`91c9cf10d5730fb57ad8f81230fd358e`. No estimator, default, public API, or
release archive changed. Existing platform checks remain valid for their
recorded sources, but do not establish the numerical property examined here.
The local environment is R 4.6.1, Rcpp 1.1.2, RcppEigen 0.3.4.0.2 and
pkgload 1.5.3. The archived script/log match the executed temporary files;
syntax, new links and recorded perturbation convergence/step stability were
checked without rerunning fits. Full tests, builds and platform checks were
not repeated for these validation-only records.

**Covariance stationarity is an unresolved release assessment.** The
mean-shift correction and covariance-score stopping rule apply only when
`chance == 0 && corr_mode == 0`. Elsewhere, the reported zero score is an
initial value for a check that is bypassed, not an evaluated zero derivative.
The default nonzero correlation selects the other path when active random
slopes are matched. Correlation shrinkage/eigenvalue repair is not an exact
joint optimizer of the Half-t/LKJ Laplace objective. The joint-delta
calculation propagates individual log-SD uncertainties and their beta
sensitivities, not a full covariance matrix including correlation parameters.
These approximation limits were already partly disclosed in the shipped
guide; the new finding below supplies a concrete numerical diagnostic.

[`review_followup.R`](review_followup.R) and its
[`execution log`](review_followup.log) retain four single-dataset probes.
Run from the repository root with
`Rscript --vanilla validation/review_followup.R`; this loads the development
source with `pkgload`, not an independently installed reference estimator.
Each design has 32 participants, 24 items, two conditions with marginal
probabilities .55/.65, participant SDs .8/.35, item SDs .6/.4, and response
seed 921126. The four cells cross guessing floors 0 and 1/3 with correlations
0 and .2. Regularization and final-stage numerical budgets match the current
native settings. All four point fits and all fixed-covariance perturbation
fits converged according to their implementation-specific stopping rules.

| Guessing floor | Correlation | Reported covariance score | Maximum scaled numerical log-SD gradient |
| ---: | ---: | ---: | ---: |
| 0 | 0 | .004927 | .004718 |
| 0 | .2 | 0 (not evaluated) | .393855 |
| 1/3 | 0 | 0 (not evaluated) | .311132 |
| 1/3 | .2 | 0 (not evaluated) | .312848 |

The negative profile objective is minus the native Laplace log likelihood
and fixed-effect/SD log priors, including the log-SD coordinate Jacobian.
Each perturbation holds fitted correlation matrices fixed and refits beta
and conditional modes. The LKJ term is constant along these directions.
Central differences at steps .001 and .0005 differ by at most 3.1e-6 in
unscaled gradients. Scaling uses the native approximate log-SD standard
errors, as in the diagonal-path regression test; it is not full observed
information. These are necessary stationarity checks using the native
profile evaluator, not independent likelihood optimization, correlation
gradient checks, HMC comparisons, or estimates of bias/coverage. They show
why convergence must not be described as common-objective stationarity
across all paths. Existing 3AFC performance evidence is retained within its
tested conditions. Resolve the objective/diagnostic contract and check its
effect on design decisions before treating the release as statistically
ready; do not silently replace correlated DGPs with zero correlations.

**Two causal interpretations in the review exceed the evidence.**

- `2 * pt(-1.96, 25)` is .06123257, but matching the saved 6.10% null rate
  does not identify 25 effective degrees of freedom. At 2,048 replications
  the MCSE of a 6.1% rate is .529 percentage points. Neither that match nor
  unbiased SE alone establishes a t pivot. The claim that every design with
  at most 30 items or 12 passages necessarily misses nominal inference is
  unsupported. Finite clusters are a plausible contributor alongside
  covariance estimation, regularization and nonlinear interval construction.
  A GLMM Satterthwaite approximation needs derivation and validation, not
  substitution of the smallest cluster count minus one. The
  [lme4 inference documentation](https://lme4.github.io/lme4/reference/pvalues.html)
  distinguishes available finite-sample procedures from GLMM Wald inference.
- Under independent normal bootstrap replicates, the relative MCSE of the
  estimated SD at B=99 is approximately 7.1%. Its effect on width attainment
  can have either sign: for a threshold of .9 times the true SD, attainment
  is 8.43% at B=99 versus .22% at B=399; for 1.1 times the true SD it is
  92.29% versus 99.74%. These are analytic chi-square examples, not observed
  GLMM rates. The saved 72-by-60 paired drop of 7.62pp (MCSE 1.68pp) is real;
  the 399-refit comparison uses only 32 outer datasets, changing 22/32 to
  26/32. It does not establish recovery over all 512 datasets, superior
  coverage, or a sufficient new bootstrap default. The implemented interval
  uses bootstrap SE with a normal critical value, not bootstrap quantiles.

**API and computation findings.** Partial presentation does not guarantee
condition/exposure balance; a fixed item seed conditions the simulation on
one schedule. However, testlets are sampled as whole units, and repeated
partial schedules can qualify for the rank-1 path. The retained 48-person,
4-item, 2-per-person example has 12 item/cell patterns, below its 24-group
limit. It satisfies the grouping condition; no native timing/flag benchmark
was performed. Its counts in the first condition are 0/1/2 for 9/30/9 people,
illustrating the separate balance issue. A one-sided zero slope SD in the
DGP does not remove that fitted slope under `match_dgp` in a within-both
design. Explicit balanced schedules and separate participant/item analysis
structures are useful enhancement candidates, not grounds to promote every
research-only procedure into the 0.1.0 API.

The external-engine joint-delta fallback uses heuristic log-SD variances
and lacks the native beta sensitivity; the .75 information floor caps
log-SD variance; and the 85-iteration coarse path remains for noncanonical
fits. These are verified code properties, not measured effect sizes on
selection or coverage. Expose the approximation/floor status and examine
affected saved cases before changing numerical rules. The 460,000-fit
calculation is correct for four fully evaluated 1,000-replication final
points plus three 200-replication sensitivity levels, each with 99 inner
fits. It excludes screening, top-ups and additional scenarios, and is not a
global maximum; curtailment can reduce it. Logit remains a legitimate
estimand, so its default is not itself a statistical defect. The surrogate
uses a monotonized mean with unconstrained kriging uncertainty; describe it
as a search approximation, not an exact constrained-GP posterior.

**Publication claims must reflect the matched comparisons.** In the saved
item-pool comparison, native-grid and mlpwr GLMM searches each produced 19
confirmed proposals from 20 runs (respectively one indeterminate result and
one search error). For Rasch, native grid had 46 confirmed/4 indeterminate;
mlpwr had 34 confirmed/5 indeterminate/1 failed and 10 explicitly unqualified
returns. GP refitting of those ten saved searches recovered confirmed
candidates without new simulation. Moreover, in the separate large L2
comparison, mlpwr confirmed 34/40 versus grid 30/40, with equal evaluation
budgets. These historical benchmarks do not support general superiority or
the claim that mlpwr proposals frequently fail a Wilson gate. Their original
kernel/scope labels remain in force. See the
[item-pool results](item_pool_search/summary.csv) and
[large L2 comparison](l2_replication_large/search_summary.csv).
The [official mlpwr IRT tutorial](https://cran.r-project.org/web/packages/mlpwr/vignettes/IRT_Vignette.html)
also explicitly demonstrates Rasch-versus-2PL and DIF power analysis.
The manuscript's useful central question remains how estimand, allocation
and success criterion change a researcher's design choice. Calibration,
failure denominators and MC uncertainty support that question; they are not
automatically new mathematical methods. Wilson bounds with alpha spending
are not an exact finite-sample simultaneous guarantee merely by construction.
Keep current user guidance separate from historical development evidence,
and distinguish public workflows from research-script examples.

Other wording qualifications: workspace reuse reduces allocations but is
not allocation-free; RNG stream/substream lengths are 2^127/2^76, whereas
the generator period is approximately 2^191. Conditional logit-range inputs
are not the distribution of observed total scores or participant averages
integrated over other random effects. Fixed-threshold cumulative-link
concavity does not prove global concavity of the full marginal objective
with estimated thresholds/covariances, nor does it make a free-discrimination
GRM identical to a common-slope ordinal GLMM. Cubic category scaling describes
a possible dense algebra cost, not every multinomial implementation's total
runtime. The ordinal/multinomial roadmap priorities are unchanged.

<a id="item-pool-planning"></a>

## Practical item-pool pilot and MML repair (2026-11-26)

- [`item_pool_planning.R`](item_pool_planning.R) fixes synthetic item IDs and
  equal content quotas before evaluating central-30, spread-30, spread-60,
  full-150, and rotating-150 (30 per person) designs. Known-bank EAP bias,
  information and conditional RMSE use numerical integration of the complete
  score distribution, with N(0,1) scoring prior and a grid refinement check.
  Recalibration uses N=120/240/480, person SD=0.4/0.8, 128 replications,
  native fixed-item MML, RMSE <= .35 and maximum error <= .60 logits, joint
  assurance .80 with pointwise Wilson 95% lower bounds. All fits and
  unevaluable banks remain in the attainment denominator. The master seed
  is 92611001 and the fixed partial-allocation seed is 37.
- Current recalibration results, per-replication RDS objects, IDs, input and
  executable hashes, and environment are in
  [`item_pool_planning_verified/`](item_pool_planning_verified/).
  There are 30 scenarios / 3,840 fits, all valid; some banks contain constant
  items and are not wholly evaluable. See
  [`calibration_summary.csv`](item_pool_planning_verified/calibration_summary.csv).
- Experimental checks in
  [`item_pool_planning/experiment_summary.csv`](item_pool_planning/experiment_summary.csv)
  used the then-current crossed kernel: N=60/120, J=30/60/150, 30 per person,
  planned marginal .50/.58 and a separate zero-difference null. Intercept
  SDs are .8/.6, slope SDs zero/zero or .25/.35, correlations zero,
  `bayes_laplace`, matching analysis and `joint_delta`. Each cell has 128
  runs, seed 92621001: 3,072/3,072 valid. The 20pp width limit is illustrative.
  Counts, Wilson intervals, coverage and response-pattern diagnostics are
  retained. Existing experimental results were reused after the fixed-bank
  convergence-only update; the original recalibration files in that directory
  are superseded by `item_pool_planning_verified`.
- This is a practical pilot. It does not optimize actual word identities,
  establish content validity, propagate calibration uncertainty, or prove
  an optimal stimulus count. The matched `mlpwr` and `irtsim` extension is
  reported below; the original pilot results remain separate.

The [conditional-precision figure](item_pool_planning.png) is produced by
[`plot_item_pool_planning.R`](plot_item_pool_planning.R). It shows why the
same item count can imply different accuracy across the ability range.

Run a new directory with
`Rscript --vanilla validation/item_pool_planning.R OUTPUT 128`; append
`calibration` to repeat only the affected calibration section. The script
refuses to overwrite an existing output directory. The released native
21-point quadrature is unchanged; conditional accuracy for a known bank uses
separate fine-grid integration, not that calibration approximation.

[`rasch_mml_acceleration.R`](rasch_mml_acceleration.R) replays the 457 saved
200-iteration failures and 112 successful controls from the existing
16-item calibration study. Final evidence is in
[`rasch_mml_verified/`](rasch_mml_verified/). The guarded fallback recovers
456/457; case 5, replication 329 remains a declared iteration-limit failure.
Previously successful controls have exactly unchanged RMSE and SD. The 564
common valid pairs with archived TAM fits have no bank-decision differences;
the largest RMSE difference is .00448 logits, SD difference .15377, and
native-minus-TAM common log-likelihood ranges from approximately -0.0000001
to +0.28024. These are not strict parameter-equivalence results. Whole-bank
accuracy and numerical validity are distinct outcomes.

The independent direct likelihood reference uses 41-point quadrature and
L-BFGS-B initialized from empirical proportions and variance .16; mirt uses
Rasch constraints and 61-point quadrature. Comparisons use the same generated
responses, mean-zero ability scale, discrimination 1 and free variance.
The declared numerical tolerances are .002 logits for item/SD differences
and .002 common log-likelihood units; discrepancies and nonconvergence must
be retained, including mirt's nonzero variance floor. Reference-only reruns
can use `Rscript --vanilla validation/rasch_mml_acceleration.R OUTPUT references`.
The earlier update-only fallback returned 457 fits but was rejected because
near-boundary estimates did not satisfy the variance score. This failure
investigation is distinct from the final convergence rule.

The acceleration is an independently implemented safeguarded squared
extrapolation of the EM map; see Varadhan and Roland (2008) in the
[SQUAREM reference manual](https://cran.r-project.org/web/packages/SQUAREM/SQUAREM.pdf).
No additional runtime package dependency is introduced. The archived
`rasch_fit_benchmark.R` now checks a manifest before reusing batches, including
the design, seeds, engines/versions, script/form hashes and installed native
binary hash. Unidentified historical batches cannot be silently resumed.

The missing-manifest and mismatched-manifest guards were exercised without
starting a simulation. Targeted tests passed 143 assertions; the macOS source
archive passed 624 assertions with three expected skips (two CRAN-only skips
and the opt-in Stan reference). Examples and vignettes built successfully.
After correcting only the NEWS heading, archive contents differed only in
NEWS and the DESCRIPTION packaging timestamp; those successful numerical
checks were reused. The final targeted package check had zero errors/warnings
and one environment NOTE for Apple's `xcrun_db` temporary cache. Remote time
and incoming checks were disabled. See ROADMAP.md for artifact hashes and
both check logs. Win-builder and other OS checks were not run in this step.

<a id="item-pool-package-comparison"></a>

## Matched item-pool comparisons with irtsim and mlpwr (2026-11-26)

The practical questions are whether an independent Rasch planning workflow
supports the same recalibration sample sizes, and whether an adaptive search
finds acceptable participant/stimulus combinations under the same simulation
budget. Estimator agreement and search performance are separate comparisons.
These comparisons do not assess person-score reliability or select actual
word identities.

**Shared Rasch conditions.** The pilot's fixed item IDs and content quotas
are retained. Five layouts crossed with person SD .4/.8 and N=120, 240, 480,
960, 1920 give 50 conditions. Responses follow the same unidimensional
Rasch model, discrimination 1, ability mean 0 and estimated variance.
`irtsim` 0.2.0 runs its public `irt_design()` / `irt_study()` /
`irt_simulate()` workflow, with the mirt Rasch backend. Its per-iteration
estimates are used to calculate the same whole-bank event as native MML:
RMSE <= .35 AND maximum absolute error <= .60 logits. Every scheduled
replication remains in the denominator; incomplete or nonconverged banks
are misses. "Valid" here means a complete finite bank with an accepted fit,
not simply that an estimator returned an object.

Both estimators receive identical responses and missingness masks. In the
installed serial irtsim generator, ability and response generation reset the
cell seed. A public custom `theta_dist` therefore draws ability from a
separate L'Ecuyer stream and restores the response RNG. The response matrices
are replayed using the installed cell-seed expression. The first two matrices
in every condition/stage are refitted through the same mirt wrapper and must
reproduce the public results within 1e-10. Per-trial data hashes, item IDs,
mask equality and mean/variance/discrimination constraints are checked.
The two packages retain their own numerical estimation defaults; aligned
models do not imply bit-identical estimates.

Screening uses 128 replications per condition. The least N whose point
attainment reaches .80 is a proposal. The union of proposals from this screen
and the two search methods is evaluated on 512 fresh response matrices per
condition, with both IRT estimators. A confirmed design requires a pointwise
Wilson 95% lower bound >= .80; an upper bound < .80 is a failure, otherwise
the result is indeterminate. This validates tested grid candidates, not the
exact minimum N between grid values or the global optimum. Paired differences
use all scheduled trials and conservative pointwise 95% intervals formed from
two exact 97.5% intervals for the discordant probabilities. RMSE differences
use common valid fits and have a separate MCSE.

**Shared search conditions.** `mlpwr` 1.1.1 receives the same native simulation,
estimator and binary success function as glmmplanr's public grid evaluation.
This explicitly includes an IRT use of mlpwr; its interface accepts a custom
simulation function and its [official IRT tutorial](https://CRAN.R-project.org/package=mlpwr/vignettes/IRT_Vignette.html)
already illustrates IRT power analysis. Here the Bernoulli event is joint
recovery/precision attainment, rather than only rejection of a null.
For Rasch, each of the ten form/SD families uses the five N values above and
cost=N. For GLMM, the grid is N=30/60/120 by J=30/60/150, with 30 stimuli per
person and a two-condition Latin square. The marginal effect is 8pp
(.50 versus .58); success requires both two-sided 5% detection and 95%
interval width <=20pp. Person/item intercept SDs are .8/.6; slope SDs are
zero/zero or .25/.35, with zero correlations and matching native
`bayes_laplace` analysis using `joint_delta`.

Each family/preference has five independent search seeds. Rasch has a nominal
600 evaluations per search and GLMM 1,080. The two GLMM preference costs are
normalized N + .02 normalized J and the reverse, identically supplied to both
methods. The discrete indices refer to exactly the same candidate designs.
The baseline is uniform grid evaluation, **not** glmmplanr's GP optimizer.
mlpwr uses its default 1D logistic or 2D GP surrogate, with setsize 60 and
four 1D or three-per-dimension 2D starting sets. Its actual callback count,
including any batch rounding, is also the grid budget; remainders are
allocated in a fixed, outcome-independent order. If a search errors early,
its paired grid receives only that many evaluations and the incomplete
budget is reported separately.

Both searches target a point attainment of .80. mlpwr warnings, errors,
predicted attainment and fallback coordinates are retained: a returned
coordinate accompanied by "No good design found" is not a qualified
recommendation. A shared independent 512-replication confirmation evaluates
the union of returned designs. GLMM confirmation also uses 512 null runs
per design and reports coverage and false-positive rates. Repeated searches
that return the same design share its confirmation; their five selection
outcomes are not five independent estimates of that design's precision.
Warnings/errors and independent confirmation results must be reported
separately. No tuning based on confirmation results is applied.

**Paired screening results.** All 6,400 response matrices were evaluated by
both IRT engines. Each had 6,095 complete valid banks. Only 2/6,400
per-replication joint decisions differed (both at N=240, full-150, SD=.8,
close to the .60 maximum-error boundary, one in each direction). Consequently,
all 50 condition-level attainment counts and all ten proposed sample sizes
agree. The largest item-estimate difference among common valid fits is
.02620 logits: these are close planning outcomes, not exact numerical
equivalence. The finite sample sizes and parameter ranges still bound any
agreement claim.

The [screening curves](item_pool_comparison.png), generated by
[`plot_item_pool_comparison.R`](plot_item_pool_comparison.R), show the
layout effect and the near-overlapping estimates. Machine-readable results
are in [screening summaries](item_pool_irtsim_screen/summary.csv),
[paired differences](item_pool_irtsim_screen/paired.csv) and
[trial-level results](item_pool_irtsim_screen/fits.csv.gz).

**Independent Rasch confirmation.** Fifteen distinct returned points were
checked on 512 fresh response matrices each (7,680 pairs). Each estimator
had 6,696 complete valid banks; the other 984 were the sparse N=120 fallback
checks. There were eight discordant joint events, confined to central-30,
SD=.8, N=120 and full-150, SD=.8, N=240. All fifteen pass/fail/indeterminate
classifications agree. Among the ten minimum-N screening proposals, each
engine confirms nine and leaves spread-60, SD=.4, N=240 indeterminate.

| Form | Proposed N | Native attainment, SD=.4 (Wilson 95%) | Native attainment, SD=.8 (Wilson 95%) |
| --- | ---: | --- | --- |
| Central 30 | 120 | 470/512, 91.8% (89.1–93.9), pass | 440/512, 85.9% (82.7–88.7), pass |
| Spread 30 | 240 | 428/512, 83.6% (80.1–86.5), pass | 451/512, 88.1% (85.0–90.6), pass |
| Spread 60 | 240 | 413/512, 80.7% (77.0–83.9), indeterminate | 433/512, 84.6% (81.2–87.4), pass |
| Full 150 | 480 | 507/512, 99.0% (97.7–99.6), pass | 502/512, 98.0% (96.4–98.9), pass |
| Rotating 150; 30/person | 1920 | 480/512, 93.8% (91.3–95.5), pass | 487/512, 95.1% (92.9–96.7), pass |

irtsim differs in this table only for central-30, SD=.8: 441/512,
86.1% (82.9–88.9), also a pass. The separately proposed spread-60,
SD=.4, N=480 passes with 506/512 (98.8%, 97.5–99.5%). This is a confirmed
tested alternative, not a search for the exact intermediate minimum.
These recalibration targets concern all retained items and do not establish
person-score accuracy or content validity.
[Confirmation results](item_pool_irtsim_confirm/summary.csv) ·
[paired differences](item_pool_irtsim_confirm/paired.csv) ·
[screening proposals and confirmation](item_pool_irtsim_confirm/confirmed_proposals.csv).

**Search diagnostics and sensitivity.** All 70 paired searches were completed
or retained as declared failures. Each method used 30,000 Rasch and 20,880
GLMM simulation evaluations. Nineteen of the twenty GLMM searches used the
full 1,080 budget. One mlpwr run stopped after 360 initial evaluations when
two rounded starting coordinates were identical; `addval()` then used a
multi-element match as recursive list indexing. The same failure is reproduced
using saved outcomes only by
[`replay_item_pool_search_error.R`](replay_item_pool_search_error.R), with the
[recorded diagnostic](item_pool_search/search-053-error-replay.log).
This is a diagnosed failure on this small discrete grid, not evidence about
all mlpwr applications. Its paired grid budget was also 360.

In all ten rotating-pool searches, the default 1D logistic surrogate issued
"No good design found" despite observed successes at N=1920. Its returned
N=120 is retained as a fallback coordinate, not a sample-size recommendation.
An exploratory [GP refit](refit_item_pool_search_surrogate.R) uses exactly the
same 600 saved outcomes per run through mlpwr's public `dat` and
`surrogate="gpr"` options, with zero additional simulation calls. All ten
refits return N=1920 with predicted attainment >=.80. This sensitivity was
chosen after inspecting default-search failures; it is not a prespecified
full GP search comparison. It shows why the default logistic failure cannot
be generalized to mlpwr's ability to plan IRT studies. The independently
confirmed precision at the returned point is evaluated separately.
[Refit results](item_pool_search/surrogate_sensitivity.csv).

Across the 50 Rasch search seeds/families, uniform grid evaluation returned
46 confirmed and four indeterminate candidates. Default mlpwr returned
34 confirmed, five indeterminate and one failed candidate, plus ten
explicitly unqualified returns. The failure was N=240 for full-150, SD=.8
(325/512, 63.5%, 59.2–67.5% on native confirmation). The exploratory GP refits
for the ten unqualified returns all select the independently confirmed
N=1920. Thus the apparent default-setting gap depends substantially on the
surrogate; it is not evidence that mlpwr cannot solve the IRT planning task.
The five search seeds and limited grid do not support a general ranking.
[Search summary](item_pool_search/summary.csv) ·
[all search/confirmation rows](item_pool_search/confirmed_search.csv) ·
[GP refit confirmations](item_pool_search/confirmed_surrogate_sensitivity.csv).

**Independent GLMM confirmation.** The union of returned experimental designs
contains three distinct conditions; all 3,072 planned/null fits were valid.
Counts below are per 512 replications, with pointwise Wilson 95% intervals.

| Slope profile; N × J | Joint detection/width attainment | 95% interval | Confirmation | Planned coverage | Null rejection |
| --- | ---: | --- | --- | ---: | ---: |
| Zero; 30 × 30 | 405/512 (79.1%) | 75.4–82.4% | Indeterminate | 485/512 (94.7%) | 16/512 (3.1%; 1.9–5.0%) |
| Zero; 60 × 30 | 494/512 (96.5%) | 94.5–97.8% | Pass | 482/512 (94.1%) | 24/512 (4.7%; 3.2–6.9%) |
| Heterogeneous; 60 × 30 | 428/512 (83.6%) | 80.1–86.5% | Pass | 490/512 (95.7%) | 29/512 (5.7%; 4.0–8.0%) |

The heterogeneous candidate passes narrowly. All planned/null confirmation
intervals met the illustrative 20pp width limit, so the observed joint-rate
differences here are driven by detection. This is not a comparison in which
the width requirement actively constrains the selected designs. The null
interval still includes rates above 7.5%; these checks do not guarantee Type I
error control at all parameter values. More stimuli, other covariances, or different effect sizes
are outside this particular search comparison. The full results, including
validity and extreme-response patterns, are in
[`item_pool_search_confirmation/summary.csv`](item_pool_search_confirmation/summary.csv).

Scripts and reproduction, from the repository root with glmmplanr, irtsim,
mirt and mlpwr installed:

```sh
Rscript --vanilla validation/compare_item_pool_irtsim.R screen
Rscript --vanilla validation/compare_item_pool_search.R full
Rscript --vanilla validation/refit_item_pool_search_surrogate.R
Rscript --vanilla validation/confirm_item_pool_search.R
Rscript --vanilla validation/compare_item_pool_irtsim.R confirm
Rscript --vanilla validation/summarize_item_pool_comparison.R
```

The first two commands are independent; confirmation requires both screening
and search files. Each computational script checkpoints completed cases,
validates source/binary/settings manifests on resumption, and records its
environment. Shared definitions are in
[`planning_comparison_helpers.R`](planning_comparison_helpers.R).
Response seeds are 9370 + case for screening and 10370 + case for confirmation;
ability-stream master seeds are 93710001 and 93720001. Search master seeds
are 9470000 + 1000 * job, with separately recorded stream offsets; independent
GLMM checks use 95710001 + 1000 * confirmation case.

This run used source builds of mlpwr 1.1.1, DiceKriging 1.6.1,
rgenoud 5.9-0.11, randtoolbox 2.0.5 and rngWELL 0.10-10 under R 4.6.1.
On this Mac, rngWELL's optional message translation was disabled at compilation
because libintl headers were unavailable; no RNG/statistical source was
changed. No R 4.5 binary was used in the completed calculations. Native MML
replay timings include two external refits per case for the pairing audit;
they are not a pure estimator-speed benchmark. These validation scripts do
not change the package's public API or estimation kernel.

<a id="l2-empirical-planning"></a>

## Empirical vocabulary banks, targeting, and a replication-planning example

The extension requested during the matched comparison asks which examinee
population and item set make a proposed sample informative. Participant count,
item-parameter precision, response-curve recovery and hypothesis-test power
are distinct quantities. The often-cited Rasch sample-size guidance itself
conditions on targeting and allowable calibration error
([Linacre's note](https://www.rasch.org/rmt/rmt74m.htm)). A per-item precision
rule is not a guarantee that every item in a long bank meets a maximum-error
criterion. For illustration, an item with marginal success probability .001
has only one expected correct response from 1,000 independent examinees and
probability .999^1000 = .3677 of no correct response.

**Public vocabulary data and scale checks.** Hamada and colleagues' study
([paper](https://doi.org/10.1093/applin/amae087),
[authors' data/code](https://osf.io/t8ax6/)) supplies a concrete item bank.
The workbook has 1,220 participants, 50 word identities and five booklet
patterns: 20 recall and 30 recognition responses per person. There are 394
missing CEFR labels; these participants are retained. One full-width blank
in the recall responses is treated as an additional missing response, not
as an incorrect answer or a new booklet.

[`prepare_l2_vocabulary_bank.R`](prepare_l2_vocabulary_bank.R) follows the
authors' separate 2PL fits for the two formats and also fits separate Rasch
working models for the native planning examples. These Rasch parameters are
new refits, not parameters reported as Rasch in the source study. mirt's
probability traces agree with `plogis(a * (theta - b))` to 2.3e-16; adding a
1.702 multiplier would change the fitted model's scale. Native and mirt Rasch
refits have maximum item differences .000382 (recognition) and .000865
(recall). The fitted Rasch variances are .7328 and 1.5339, respectively.
The two formats retain separate reference scales and are not equated here.

The empirical 2PL discriminations range from .307 to 2.727 for recognition
and .366 to 3.634 for recall; difficulties range from -10.799 to 1.629 and
-4.589 to 2.456. These are estimates from this dataset, not known population
truth. They become fixed generating values for the conditional simulation.
Original raw responses and test texts are kept outside the repository;
[derived parameters, booklets, provenance and fit audits](l2_empirical_bank/)
record item IDs and numerical summaries. Source downloads are
[workbook](https://osf.io/download/rvcfs/) and
[author code](https://osf.io/download/uekcw/).

**Expanded targeting pilot.**
[`l2_targeting_simulation.R`](l2_targeting_simulation.R) crosses two formats,
five layouts, six ability profiles and four sample sizes (N=100/200/500/1000):
240 native Rasch conditions ×64 replications = 15,360 evaluations. Layouts
are all 50 items, 25 items retaining five per frequency band, the five
hardest per band (25 total), and the original five booklet masks allocated
equally or in their observed proportions. Candidate IDs are fixed before
synthetic responses are generated. Ability profiles use the refitted Rasch
SD, half that SD, means -3/+3/-5, and mean -3 with half SD. They are sensitivity
assumptions, not elementary-school/high-school or CEFR conversions.

A subset of 24 conditions (reference, lower/narrow and severe lower profiles;
N=200/1000; full or observed-proportion booklets) has 32 paired native/mirt
Rasch runs and 32 separate correctly specified mirt 2PL runs each. All 768
Rasch pairs have identical response hashes and validity decisions; one joint
difficulty decision differs. The 768 external 2PL evaluations do not add 2PL
estimation to glmmplanr. The mean shifts are in each model's reference scale;
these separate DGP arms are not a ranking of Rasch against 2PL. The two arms
reuse the underlying normal/uniform random streams within a replication;
they generate different response matrices and are not independent of each other.

For Rasch, the difficulty target is RMSE <=.35 AND maximum error <=.60.
For 2PL, difficulty recovery is reported separately from joint recovery that
also requires discrimination RMSE <=.25 and maximum error <=.50. Response
curve RMSE <=.05 is a separate illustrative target, evaluated over nine
ability points within two SDs of the generating mean. This is not individual
ability-score precision. Known generating location/scale transformations
are used for recovery scoring; no unknown population is linked to an external
scale without anchors. All scheduled runs remain denominators. A constant
item makes a complete bank unevaluable, even if other items can be fitted.

At N=1000 with all 50 items, native Rasch difficulty attainment was:

| Format; ability profile | Complete valid banks | Successes / 64 | Attainment (Wilson 95%) | Mean constant items |
| --- | ---: | ---: | --- | ---: |
| Recall; reference | 64 | 60 | 93.8% (85.0–97.5%) | 0.00 |
| Recall; mean -3 | 56 | 4 | 6.3% (2.5–15.0%) | 0.13 |
| Recall; mean -3, half SD | 24 | 1 | 1.6% (0.3–8.3%) | 0.84 |
| Recall; mean -5 | 0 | 0 | 0.0% (0.0–5.7%) | 4.73 |
| Recognition; reference | 57 | 24 | 37.5% (26.7–49.7%) | 0.11 |
| Recognition; mean -3 | 64 | 61 | 95.3% (87.1–98.4%) | 0.00 |

Thus lowering the ability mean can harm calibration for one bank and improve
it for an easy bank that otherwise produces ceiling responses. This is a
within-format targeting comparison, not an absolute comparison of learner
ability between recognition and recall. At N=1000, reference 2PL recall
runs recovered the response curves within .05 in 32/32 runs, but met the
whole-bank difficulty/discrimination criterion in only 1/32 (3.1%, 0.6–15.7%).
Recognition gave 29/32 and 0/32 (0%, 0–10.7%), respectively. In particular,
very easy items with shallow slopes can have poorly recovered difficulty
thresholds far outside the sampled ability range while their response
probabilities within that range are recovered well. "Sufficient N" requires
the intended use and tolerated error to be specified.

Exposure also changes the answer. For reference recall at N=1000, the same
50-item bank attains the criterion in 60/64 complete forms (93.8%, 85.0–97.5%),
45/64 equal-booklet runs (70.3%, 58.2–80.1%), and 36/64 observed-proportion
booklet runs (56.3%, 44.1–67.7%). Each person answers 20 recall items in the
booklets, giving 400 responses per item under equal allocation and 316–518
under the observed proportions, versus 1,000 under full presentation.
The frequency-balanced 25-item form attains 62/64 (96.9%, 89.3–99.1%), but
calibrating fewer parameters is a different requirement: this does not
establish equivalent vocabulary-size measurement or content coverage.

This is a **64/32-replication pilot**, not a field-wide sample-size guarantee.
The source item parameters are treated as fixed; their sampling uncertainty,
non-normal population distributions, a guessing asymptote, DIF and actual
learner-age mappings are not propagated. A larger confirmatory study should select decision-relevant
boundaries from this pilot, use fresh response streams, and, when bank
uncertainty is part of the question, resample calibration participants within
booklets. This scope keeps sampling uncertainty separate from post hoc form
selection. [Complete summaries](l2_targeting_pilot/summary.csv) ·
[trial results](l2_targeting_pilot/fits.csv.gz) ·
[paired Rasch checks](l2_targeting_pilot/paired_rasch.csv) ·
[targeting figure](l2_targeting.png).

**Published-mean replication planning.** Huang et al.'s
[2025 study](https://doi.org/10.3390/bs15111540) supplies the existing
vignette's three-round computer-test means (.42 restudy, .50 retrieval).
The source study had 108 participants across two learning-round groups and
60 Swahili–Chinese pairs. Its data are available on author request and its
random-effect variances are not printed. This example therefore plans a
focused two-condition study from published means, not an exact replication
of the full interaction test or its four-sublists schedule. Per-person
accuracy SDs are not substituted for latent random-effect SDs.

[`l2_replication_planning.R`](l2_replication_planning.R) compares
N=54/72/108/162/216 and J=30/60/90, with 30 items per person in a Latin square.
Four explicit planning profiles are: published 8pp difference with intercept
SDs .8/.6 and no slopes; the same means with slope SDs .25/.35; a smaller
5pp difference with those slopes; and 8pp with larger intercept SDs 1.1/.9
and slopes .4/.5. All correlations are zero and analysis matches generation.
The joint target is two-sided 5% detection AND a 95% interval no wider than
10pp. Both normalized participant-priority and stimulus-priority costs are
used; each method receives 1,080 simulations for each profile/preference.

Following the [mlpwr tutorial](https://CRAN.R-project.org/package=mlpwr/vignettes/IRT_Vignette.html),
a function generates data, fits the analysis and returns a logical success
indicator to `find.design()`. The tutorial's IRT model-comparison/DIF testing
goals differ from this recovery/precision goal. Here mlpwr uses a GP and
explicit unique initial designs through its public `dat` argument, including
all initial draws in its budget. This avoids the duplicate integer initial
points observed earlier. There is one illustrative search seed per profile
and preference; this is not another search-stability benchmark. The union
of proposals receives 512 fresh planned and 512 null replications per point.

All eight searches per method used exactly 1,080 evaluations (8,640 each).
There were no search errors. The union contained eight distinct candidates;
all 8,192 planned/null confirmation fits were valid. The independently checked
joint attainment is shown below; confidence intervals are pointwise Wilson
95% Monte Carlo intervals, based on 512 planned runs per point.

| Planning profile | N × item pool J (30/person) | Joint attainment | 95% interval | Planning confirmation |
| --- | --- | ---: | --- | --- |
| Published 8pp; intercepts only | 54 × 30 | 485/512, 94.7% | 92.4–96.4% | Pass |
| 8pp; slopes .25/.35 | 72 × 60 | 434/512, 84.8% | 81.4–87.6% | Pass |
| 8pp; slopes .25/.35 | 108 × 30 | 466/512, 91.0% | 88.2–93.2% | Pass |
| 5pp; slopes .25/.35 | 162 × 30 | 403/512, 78.7% | 75.0–82.0% | Indeterminate |
| 5pp; slopes .25/.35 | 162 × 60 | 433/512, 84.6% | 81.2–87.4% | Pass |
| 8pp; larger variances | 72 × 60 | 409/512, 79.9% | 76.2–83.1% | Indeterminate |
| 8pp; larger variances | 72 × 90 | 464/512, 90.6% | 87.8–92.9% | Pass; inference concerns below |
| 8pp; larger variances | 162 × 30 | 486/512, 94.9% | 92.7–96.5% | Pass |

In this single illustrative search per preference/profile, mlpwr's eight
proposals pass the planning confirmation; the grid's five pass and three
are indeterminate. Two preferences can return the same point, so these
counts do not represent eight independent holdout studies. Under the 5pp
profile, the grid's 162×30 candidate is indeterminate, whereas mlpwr's
162×60 passes. Under larger variances, its 72×90 candidate passes the
planning target, but **is not an unqualified recommendation**: null rejection
is 37/512 (7.2%, 5.3–9.8%) and planned coverage is 471/512
(92.0%, 89.3–94.0%). Both intervals exclude their nominal values in the
unfavorable direction. The 162×30 alternative has null rejection 25/512
(4.9%, 3.3–7.1%) and planned coverage 482/512 (94.1%, 91.8–95.9%).
These diagnostics concern the common native estimator, not the search
algorithm. No estimator changes were made to improve this benchmark.

**Update from the larger independent run:** the 72×90 warning was not
reproduced in 2,048 new trials per condition (null 5.0%, coverage 94.7%).
Conversely, 162×30 gave null rejection 7.1% (6.0–8.3%). The historical
512-run results above remain intact; use the
[larger confirmation](#l2-large-confirmation) for the updated assessment.

Unlike the initial 20pp-width pilot, the 10pp width restriction is active:
for the moderate-slope 72×60 design, detection is 92.4%, width attainment
91.8%, and joint attainment 84.8%. A power-only recommendation would miss
that distinction. These are finite-grid candidates under assumed variances;
they neither recover the unpublished variances nor determine an exact
minimum replication size. The full original study's learning-method ×
learning-round interaction is a different target and requires its own plan.
[Search and confirmation decisions](l2_replication_full/confirmed_search.csv) ·
[all planned/null diagnostics](l2_replication_full/confirmation.csv).

The [analytic targeting audit](l2_targeting_pilot/full_form_targeting.csv)
compares the 48 complete-form response means with 61-node integration;
the largest discrepancy is 2.09 empirical Monte Carlo standard errors.
Expected constant-item counts also follow the independent-examinee
probabilities p_bar^N + (1-p_bar)^N. The reported known-ability Fisher
information is an optimistic diagnostic; its inverse square root is not
presented as an estimated MML standard error. These checks, pairing and
budget assertions are implemented in
[`summarize_l2_planning.R`](summarize_l2_planning.R).

Reproduction, after downloading the original workbook outside the repository:

```sh
Rscript --vanilla validation/prepare_l2_vocabulary_bank.R PATH_TO_WORKBOOK
Rscript --vanilla validation/l2_targeting_simulation.R pilot
Rscript --vanilla validation/l2_replication_planning.R full
Rscript --vanilla validation/summarize_l2_planning.R
Rscript --vanilla validation/plot_l2_targeting.R
```

The prepared-bank output records the workbook hash and calibration code;
subsequent runs check source, parameter-file and executable manifests before
resuming completed cases. The authors' complete analysis script is inspected
for model/data conventions, not executed wholesale. No public API or
estimation kernel is changed by these examples.

<a id="l2-large-confirmation"></a>

### Larger Independent Confirmation of the Vocabulary Examples

The authorized large run preserves the pilot results and uses fresh streams.
It asks whether the targeting, exposure and replication-planning conclusions
persist at higher Monte Carlo precision, and whether N=2000 changes the bank
calibration decision. Source item parameters remain fixed; this is not a
bootstrap over the original calibration sample or an unconditional guarantee
for future vocabulary tests. No estimator or public API is changed.

The native Rasch grid has 300 conditions: two formats × five layouts × six
ability profiles × N=100/200/500/1000/2000, with 1,024 replications per cell.
The 70 external conditions have 256 paired mirt Rasch and 256 correctly
specified mirt 2PL evaluations each. They comprise both formats, N=200/500/1000/2000,
reference/lower/lower-narrow/severe-lower profiles, and full/observed-booklet
layouts (64 cells), plus reference balanced-25 forms at N=500/1000/2000
(six cells). All item IDs are fixed before simulation. Within a replication,
the two generating models share underlying random streams but produce
different response matrices; only the native/mirt Rasch comparison shares
the exact response matrix. Models and formats retain their own scales.

Targets and failure denominators are unchanged: difficulty RMSE <=.35 AND
maximum error <=.60; joint 2PL recovery additionally requires discrimination
RMSE <=.25 and maximum error <=.50. Probability-curve RMSE <=.05 over nine
points within two generating SDs, averaged over all items and those points,
is reported separately; it is not a maximum-item probability-error guarantee. Every scheduled
evaluation remains in the denominator, including incomplete constant-item
banks and nonconvergence. Mirt uses the same default fitting controls as the
pilot; the study does not rank alternative optimizer settings. The known
generating origin/scale is restored for recovery scoring, without claiming
operational anchor-free equating.

For vocabulary learning, all 60 profile/design combinations (four variance/
effect profiles × five participant counts × three item pools) receive 2,048
planned and 2,048 null replications. Every person answers 30 items. The same
native `bayes_laplace` / `joint_delta` calculation is used, with no inner
bootstrap refits. Confirmation passes when the Wilson lower bound for
detection AND 95% interval width <=10pp is at least .80; an upper bound below
.80 fails, otherwise the result is indeterminate. False-positive rate,
coverage, bias, empirical SD and mean SE are separate diagnostics. Detection
and joint attainment use all scheduled runs; coverage uses valid runs with
the valid count reported. Pointwise 95% intervals are not familywise-adjusted
tests over all designs.

The search comparison repeats each of four profiles × two cost preferences
under five fresh search seeds. Each method receives 1,080 evaluations per
search, including the 150 explicitly initialized evaluations. Early errors
retain matching actual budgets and no-candidate outcomes. Both methods'
proposals use the same independent 60-design confirmation grid. Repeated
proposals do not multiply its independent evidence. This is uniform grid
search versus mlpwr GP, not glmmplanr's GP versus mlpwr's GP. Integer indices
map to the five participant counts and three item-pool counts; the GP
operates on these indices, not on every integer N between the endpoints.

Random-stream roots are 110200001 (bank responses), 410500001 (external fits),
610100001 (search), and 810900001 (learning confirmation), with nonoverlapping
seed ranges and condition/replication offsets recorded in the scripts.
The nominal MCSE near 80% is 1.25pp for 1,024 replications and 2.5pp for 256;
near a 5% null rate it is 0.48pp for 2,048. The run keeps the checked native
executable and package versions used by the pilot. Bank cases checkpoint
every 128 replications; resuming a partial checkpoint reproduced every keyed
trial value, response hash and summary in the smoke check.

Confirmation records the inputs to the public function's summary through
local function copies; it does not modify the installed namespace. The smoke
check reproduced the public validation output exactly. Before confirmation
began, this logging was added after preserving 14 finished searches. All
search expressions and other source hashes were unchanged; the original
manifest and runner are retained with the results. Completed searches were
not repeated.

```sh
Rscript --vanilla validation/l2_targeting_simulation.R large validation/l2_targeting_large 6
Rscript --vanilla validation/l2_replication_planning.R large validation/l2_replication_large 6
Rscript --vanilla validation/summarize_l2_large.R
Rscript --vanilla validation/plot_l2_large.R
```

**Completed bank results.** All 307,200 native, 17,920 paired mirt Rasch,
and 17,920 mirt 2PL evaluations are saved. Complete valid banks number
148,296, 10,569 and 7,559, respectively; these aggregates combine deliberately
extreme and reference populations and are not general convergence rates.
Native incomplete banks account for 158,902 misses. Two further native runs
reached the iteration limit, but already had 14 and 13 constant items, so
even successful fitting of their remaining items could not satisfy the
complete-bank target. Mirt Rasch has 7,351 constant-item banks and no other
fit failures; mirt 2PL has 9,623 constant-item banks and 738 nonconvergences
under its default fitting controls. Every failure remains a miss.

All 17,920 native/mirt Rasch response hashes and validity decisions match.
Five trial-level attainment decisions differ (four native-only and one
mirt-only), but all 70 condition-level confirmation classifications agree
when both methods use the same 256 paired trials. The separate native grid
uses 1,024 trials per condition; its narrower intervals are not used as the
denominator for that paired method comparison.

For all 50 items and N=2000, the native difficulty target gives:

| Format; ability profile | Complete valid / 1024 | Attainment | Wilson 95% MC interval |
| --- | ---: | ---: | --- |
| Recall; reference | 1024 | 1023/1024, 99.9% | 99.4–100.0% |
| Recall; mean -3 | 1010 | 504/1024, 49.2% | 46.2–52.3% |
| Recall; mean -5 | 89 | 0/1024, 0% | 0–0.4% |
| Recognition; reference | 1007 | 666/1024, 65.0% | 62.1–67.9% |
| Recognition; mean -3 | 1024 | 1018/1024, 99.4% | 98.7–99.7% |
| Recognition; mean -5 | 1005 | 363/1024, 35.4% | 32.6–38.4% |

This confirms that the direction of a targeting effect depends on the bank.
At N=1000, reference recall attainment is 984/1024 (96.1%, 94.7–97.1%) under
full presentation, 660/1024 (64.5%, 61.5–67.3%) under equal booklets, and
582/1024 (56.8%, 53.8–59.8%) under observed-proportion booklets. At N=2000,
the two booklet allocations rise to 955/1024 (93.3%, 91.6–94.6%) and
915/1024 (89.4%, 87.3–91.1%). Under the reference recognition distribution,
neither the full form nor either booklet allocation passes at any tested N
through 2000. The frequency-balanced 25-item form first passes at N=2000;
the harder-25 form passes at N=200. These are the smallest tested passing
points, not exact minima or evidence of equivalent vocabulary-size measurement.
For recall balanced-25, N=500 is indeterminate (834/1024, 81.4%, 78.9–83.7%),
so passing N=1000 does not prove that every smaller sample would be inadequate.

The distinction between parameters and prediction persists for 2PL. At
N=2000, reference recall/full50 has 256/256 valid banks and probability-curve
successes, but only 54/256 joint difficulty/discrimination successes
(21.1%, 16.5–26.5%). Reference recognition/full50 has 253/256 valid banks
and probability successes, but 0/256 joint parameter successes (0–1.5%).
At N=1000, the corresponding recall joint parameter count is 7/256
(2.7%, 1.3–5.5%), despite 256/256 probability successes. These are the stated
whole-bank tolerances, not the itemwise precision criteria behind generic
sample-size rules. They do not demonstrate that 2PL fitting or person-score
measurement is generally impossible at those sample sizes.

A descriptive breakdown of the saved trials makes the aggregation issue
explicit. At N=2000, recall difficulty RMSE alone meets .35 in 249/256 runs,
and discrimination RMSE alone meets .25 in 256/256; adding the maximum-error
limits reduces joint success to 54/256. For recognition, difficulty RMSE
alone succeeds in only 2/256, while discrimination RMSE alone succeeds in
246/256. Its most frequent worst-recovered difficulties correspond to source
thresholds around -7 to -11, outside the evaluated population range. Thus
the parameter/prediction contrast reflects both the chosen aggregate versus
maximum-error criteria and extrapolation of difficulty thresholds. This
post-run diagnostic uses existing fits and does not replace the prespecified
planning criterion. Probability RMSE is averaged across items and ability
points; success does not certify every individual item curve.

The numerical audit recomputed all successes, denominators and decisions
from saved trials. Complete-form response means agree with a separate
61-node integral to within 3.31 empirical Monte Carlo standard errors.
[Full summary](l2_targeting_large/summary.csv) ·
[separate recovery metrics](l2_targeting_large/recovery_metrics.csv) ·
[smallest tested passing points](l2_targeting_large/smallest_confirmed.csv) ·
[paired estimator comparison](l2_targeting_large/paired_summary.csv) ·
[failure counts](l2_targeting_large/failure_counts.csv) ·
[native numerical failures](l2_targeting_large/native_numerical_failures.csv) ·
[trial results](l2_targeting_large/fits.csv.gz) ·
[targeting figure](l2_targeting_large.png).

**Completed learning-design results, before the covariance correction.** All 245,760 confirmation fits
(60 designs × two conditions × 2,048) are valid. Independent recalculation
from saved estimates/SEs reproduces detection, coverage, width attainment,
joint counts and decisions. The table shows selected, profile-specific
candidates, with pointwise Wilson 95% Monte Carlo intervals throughout.

| Assumed profile | N × pool J | Joint successes / 2048 | Joint rate (95%) | Null rejection (95%) | Planned coverage (95%) |
| --- | --- | ---: | --- | --- | --- |
| 8pp; intercepts only | 54 × 30 | 1935 | 94.5% (93.4–95.4%) | 4.7% (3.9–5.7%) | 95.4% (94.4–96.2%) |
| 8pp; moderate slopes | 72 × 60 | 1730 | 84.5% (82.8–86.0%) | 4.0% (3.2–4.9%) | 95.8% (94.8–96.5%) |
| 8pp; moderate slopes | 108 × 30 | 1906 | 93.1% (91.9–94.1%) | 5.0% (4.2–6.1%) | 94.9% (93.8–95.7%) |
| 5pp; moderate slopes | 162 × 60 | 1756 | 85.7% (84.2–87.2%) | 5.0% (4.2–6.1%) | 93.9% (92.8–94.9%) |
| 5pp; moderate slopes | 216 × 60 | 1882 | 91.9% (90.6–93.0%) | 5.6% (4.7–6.6%) | 95.2% (94.2–96.1%) |
| 8pp; larger variances | 72 × 90 | 1867 | 91.2% (89.9–92.3%) | 5.0% (4.1–6.0%) | 94.7% (93.7–95.6%) |
| 8pp; larger variances | 162 × 30 | 1928 | 94.1% (93.0–95.1%) | 7.1% (6.0–8.3%) | 94.3% (93.2–95.2%) |

The pilot's high-variance 72×90 warning was not reproduced, whereas the
162×30 alternative now has an elevated null rate. The conservative 162×60
candidate passes the planning criterion but its coverage interval is below
95%; its 216×60 alternative has no adverse nominal-departure flag in that
profile. These are updates from fresh streams, not changes to the estimator
or generating assumptions. The historical pilot estimates are retained.
Under moderate slopes at 72×60, detection is 93.7%, width attainment 89.8%,
and joint attainment 84.5%, again showing why power alone is insufficient.

Across all 60 designs, 39 pass the planning criterion. Eighteen have a null
interval wholly above .05 or a planned-coverage interval wholly below .95;
14 of those 18 also pass planning. The flags comprise four conservative,
nine high-variance and five moderate-slope designs; none occurs in the
intercept-only profile. These are pointwise diagnostics, not 60-design
familywise tests or automatic declarations that a method is invalid.
They concern nominal calibration, separately from historical floor-based
coverage checks and the power/width planning target.

A descriptive intersection across the four assumptions finds four layouts
that pass all four planning checks: 162×60, 162×90, 216×60 and 216×90.
Every one has at least one adverse nominal-departure flag among the four
profiles. Thus this grid does not supply a single design that both passes
all four planning checks and has no such flag. Absence of a flag would not
itself prove exact nominal calibration, and this post-run intersection is
not a simultaneous confidence guarantee across scenarios. Increasing N or
J alone is not a demonstrated remedy for all observed interval limitations.

Saved estimates help characterize the discrepancy without assigning an
unverified cause. At high-variance 216×90, planned bias is +0.176pp, empirical
SD is 1.548pp, RMS estimated SE is 1.480pp, and coverage is 93.6%. The
standardized-error SD is 1.059. These observations warrant focused checking
of the native uncertainty approximation; they do not isolate regularization,
Laplace approximation or delta-method contributions. The two native
fixed-bank MML failures above cannot explain this separate GLMM result.

All 40 search pairs completed without errors and each used 1,080 evaluations
per method (43,200 per method). The grid returned 30 independently confirmed,
six indeterminate and four failed proposals; mlpwr returned 34 confirmed and
six indeterminate proposals. All returned a candidate. Mlpwr recorded
boundary warnings in 29 searches; these indicate the specified grid limits,
not failed evaluations. Five seeds per profile/preference and shared
confirmation points do not establish general optimizer superiority.
Both search methods use the same estimator, so passing more planning
checks does not resolve its interval-calibration limitations.

[All design decisions](l2_replication_large/design_decisions.csv) ·
[diagnostic estimates and SEs](l2_replication_large/diagnostic_summary.csv) ·
[scenario intersection](l2_replication_large/scenario_robustness.csv) ·
[search summary](l2_replication_large/search_summary.csv) ·
[candidate confirmation](l2_replication_large/confirmed_search.csv) ·
[trial results](l2_replication_large/fits.csv.gz) ·
[planning figure](l2_replication_large.png) ·
[inference diagnostics figure](l2_replication_large_diagnostics.png).

The earlier pilot scripts are preserved inside their result directories as
`source_script.R`; changed runner hashes intentionally prevent silently
resuming a pilot directory with the large-run code.

<a id="l2-uncertainty-correction"></a>

### Diagnosing Native GLMM Interval Calibration

This follow-up asks whether the covariance update and uncertainty approximation
support the participant/item decisions in the vocabulary example. It retains
the earlier executable and results. All generating profiles, probability
contrasts, item assignments and priors are unchanged; the focused scenarios
have zero random-effect correlations and no guessing lower asymptote.

The saved 245,760 estimates are first inspected for bias, standardized-error
tails, and coverage stratified by estimated SE. Removing the observed mean
bias or substituting the empirical across-replication SD is descriptive only:
neither is an available correction for a new dataset. The six diagnostic
designs are intercept-only 54×30, moderate-slope 72×60, conservative 162×60,
and high-variance 72×90, 162×30 and 216×90. Their selection uses previous
results, so it is not evidence for an arbitrary design population.

The first 512 saved null and planned trials at each design are reconstructed
from their original RNG substreams. The current point estimate and SE must
match the stored results within 1e-7 percentage points. On exactly the same
responses, a covariance-fixed fit uses the generating covariance matrices
and propagates only fixed-effect uncertainty. This oracle comparison is a
mechanism check, not a proposed practical estimator. An additional comparison
uses the first 128 trials with 399 existing native bootstrap refits per trial;
its smaller outer replication count is not a nominal-coverage guarantee.
This uses the former estimator, including its former bootstrap refits. Each
bootstrap SE requires at least 90% of the 399 inner refits to be valid. The
interval is centered on the original point estimate and uses normal 95%
limits with the SD of valid parametric refits as its SE. Conclusions concern
this interval construction; percentile and bootstrap-t intervals were not tested. The wrapper does not retain the exact number of successful inner refits here;
outer validity therefore establishes that threshold, not zero inner failures.
It does not evaluate bootstrap intervals under the corrected estimator.

Full profile curvature is first evaluated on four trials per design/condition.
The same native Laplace likelihood is then independently optimized over
log-SDs on eight trials per design/condition. The objective retains the normal
fixed-effect prior, Half-Student-t(3) SD prior, and log-SD Jacobian used by the
native update. A separate sparse conditional Hessian reconstructs the marginal
Laplace likelihood, and central differences at two step sizes check numerical
stability. Failed line-search proposals receive infinite objective values;
they are not accepted optimizers or omitted simulation outcomes.

The candidate carries movement of the conditional random-effect mode into
the covariance score. In standardized coordinates, the correction to each
second moment includes `v * delta' + delta * v'`, where `delta` solves the
conditional Hessian against the log-determinant derivative. Fixed effects
and covariance components are updated toward the same approximate marginal
objective. A scaled covariance-score check supplements the parameter-change
stopping rule. Its numerical tolerance is checked against independent profile
optimization; it is separate from all statistical acceptance criteria.
The candidate is restricted to zero-correlation, zero-guessing paths; the
correlated and guessing-adjusted paths have paired preservation checks.

For a diagonal variance component, the check can be derived without the EM
implementation. In centered coordinates let `H = Z' W Z + G^-1` be the
conditional precision at the mode `u`, and define
`c = -0.5 * diag(Z H^-1 Z') * dW/deta` and `delta = H^-1 Z' c`.
Differentiating the Laplace likelihood with respect to `log(tau_j)` gives
the sum over units of
`(u_j^2 + (H^-1)_jj + 2*u_j*delta_j)/tau_j^2 - 1`.
The last term in the numerator carries the mode's movement through the
log-determinant and was absent from the former second-moment update.
For the native log-SD MAP objective, the additional prior/Jacobian score is
`1 - 4*r/(1+r)`, with `r = tau_j^2/(3*prior_tau_sd^2)`.
The code computes the equivalent correction in standardized coordinates
using the existing Schur factors. This identity concerns the Laplace
objective, not an exact posterior second moment.

The numerical score limit is 0.005 after scaling by the existing approximate
log-SD information. A preliminary 0.001 limit caused iteration-limit failures
even near independently profiled optima; the adopted numerical tolerance was
chosen before fresh confirmation and does not relax any coverage, rejection,
or planning threshold. The efficient implementation holds beta fixed during
conditional-mode updates, then takes one marginal beta step per covariance
update and retains the existing final polish. It avoids repeatedly resetting
beta to the conditional joint mode.

Across the 96 diagnostic trials, all development-candidate fits are valid and
their probability-effect estimates differ from independently optimized
profiles by less than 0.00023 percentage points. The separate Hessian
calculation on the first trial in each of the 12 design/condition cells
agrees with native Laplace log likelihood within 1.55e-11.
At a fixed point estimate, adding all log-SD cross-curvature changed SE
little; this does not support a uniform SE multiplier. Independent finite
differences also verify the corrected covariance score for testlet fits,
both with and without regularization. These numerical comparisons establish
consistency of the approximate objective, not exact marginal likelihood or
frequentist 95% coverage.

Independent confirmation uses 2,048 null and 2,048 planned trials at each of
the six designs, with a new seed root 912100001. Previous and candidate
executables use the same responses. Every attempted trial remains in the
detection and joint-attainment denominators, and coverage reports its valid
denominator. The prespecified practical checks are valid rate >=99%, Wilson
95% lower coverage bound >=92.5%, and Wilson 95% upper null-rejection bound
<=7.5%. These are tolerance checks, not proof of exact 95%/5% calibration;
pointwise nominal-departure flags remain separate. Joint planning still uses
detection AND interval width <=10pp, with Wilson lower bound >=80%. The
nominal null-rate MCSE is about 0.48pp. This is confirmation of selected
designs, not a repeat of the entire 60-cell search grid.

**Independent confirmation completed.** The comparison contains 24,576
matched datasets (six designs × null/planned × 2,048), each analyzed by both
versions. All response hashes and unrounded planning decisions were audited.
The previous version has 24,576 valid fits; the correction has 24,556.
All 20 iteration-limit failures occur at high-variance 216×90: 13 null and
seven planned. They remain in all-attempt detection and joint denominators;
coverage uses 2,035 and 2,041 valid fits in those cells, and 2,048 elsewhere.

| Profile; N×J | Null rejection %, old → corrected (corrected 95% MC interval) | Planned coverage %, old → corrected (corrected 95% MC interval) | Corrected valid fits, null / planned |
| --- | --- | --- | --- |
| 8pp, intercepts only; 54×30 | 4.35 → 4.35% (3.54–5.32) | 95.46 → 95.41% (94.42–96.23) | 2048 / 2048 |
| 8pp, moderate slopes; 72×60 | 4.10 → 3.96% (3.19–4.89) | 95.80 → 95.95% (95.00–96.72) | 2048 / 2048 |
| 5pp, moderate slopes; 162×60 | 5.27 → 4.98% (4.12–6.01) | 95.02 → 95.26% (94.26–96.10) | 2048 / 2048 |
| 8pp, larger variances; 72×90 | 6.25 → 5.86% (4.92–6.96) | 95.31 → 95.46% (94.47–96.28) | 2048 / 2048 |
| 8pp, larger variances; 162×30 | 6.84 → 6.10% (5.15–7.22) | 93.55 → 94.48% (93.41–95.39) | 2048 / 2048 |
| 8pp, larger variances; 216×90 | 6.25 → 5.52% (4.61–6.59) | 93.75 → 94.17% (93.07–95.11) | 2035 / 2041 |

| Profile; N×J | Joint attainment %, old → corrected | Corrected 95% MC lower bound | Cached-fit median CPU seconds, old → corrected |
| --- | --- | ---: | --- |
| 8pp, intercepts only; 54×30 | 94.68 → 94.68 | 93.62% | 0.0020 → 0.0030 |
| 8pp, moderate slopes; 72×60 | 84.57 → 81.79 | 80.06% | 0.0530 → 0.1150 |
| 5pp, moderate slopes; 162×60 | 85.84 → 85.25 | 83.65% | 0.0985 → 0.2925 |
| 8pp, larger variances; 72×90 | 91.60 → 89.55 | 88.15% | 0.0675 → 0.1305 |
| 8pp, larger variances; 162×30 | 94.63 → 92.77 | 91.57% | 0.0190 → 0.0455 |
| 8pp, larger variances; 216×90 | 99.85 → 99.51 | 99.10% | 0.2050 → 0.4275 |

Every corrected cell satisfies the prespecified **observed** valid-rate
threshold of 99%, and all selected designs satisfy the practical coverage
floor and null ceiling using their pointwise Wilson intervals. These are
not exact nominal-calibration claims. High-variance 162×30 still has null
rejection 6.10% (5.15–7.22%), above nominal 5%; high-variance 216×90 planned
coverage is 94.17% (93.07–95.11%). Absence of an adverse flag does not prove
95% coverage. Valid-rate intervals and rejection conditional on a valid fit
are also retained; the 216×90 null conditional rejection is 113/2035 (5.55%),
versus 113/2048 (5.52%) with all attempts in the denominator.

All six planned candidates still pass detection AND 95%-interval width
<=10pp with a Wilson lower bound >=80%. The moderate-slope 72×60 design now
has only 80.06% at that lower bound, so its planning margin is small. Its
joint attainment falls from 84.57% to 81.79%. Larger-variance designs also
lose some joint attainment as estimates and intervals change. This is a
meaningful planning consequence even where the categorical decision stays
PASS; it is not a reason to claim a new minimum N or J for the full grid.

Paired changes isolate the implementation effect from Monte Carlo variation.
For high-variance 162×30, null rejection changes by -0.732 percentage points
(paired MCSE 0.188pp) and planned coverage by +0.928pp (MCSE 0.243pp).
The earlier 5pp/162×60 coverage warning does not recur even in the previous
executable's fresh data (95.02%); that change cannot be credited to the fix.
The complete paired table distinguishes coverage across all attempts from
coverage changes restricted to pairs where both fits are valid.

Timing uses 16 separate, matched datasets per design, with simulation and
workspace construction outside the timed cached fit. CPU means user plus
system time, not end-to-end planning wall time. Median fit CPU rises
about 2–3-fold in the slope scenarios shown above, with no material increase
in the small intercept-only scenario. An additional 16 correlated and 16
guessing-adjusted fits have exactly unchanged estimates and SEs. These are
preservation checks, not new calibration evidence for those paths. Timing
varies with concurrent load and processor frequency; it is not a precise
end-to-end runtime promise. The
fixed-bank MML function and following code are byte-identical to the old
source; its large Rasch comparisons were not rerun.

[Independent counts and intervals](l2_uncertainty_comparison/independent_summary.csv) ·
[paired changes and Monte Carlo SEs](l2_uncertainty_comparison/paired_differences.csv) ·
[retained failures](l2_uncertainty_comparison/failures.csv) ·
[independent profile comparison](l2_uncertainty_comparison/independent_profile_comparison.csv) ·
[CPU comparison](l2_uncertainty_comparison/cpu_comparison.csv) ·
[previous trial results](l2_uncertainty_confirmation_baseline/fits.csv.gz) ·
[corrected trial results](l2_uncertainty_confirmation_guarded/fits.csv.gz) ·
[confirmation figure](l2_uncertainty_confirmation.png).

**Strong-contrast counterexample and safeguard.** Reusing a small final
beta polish as the primary beta update introduced a false-convergence case
in the development candidate. With 48 people, 24 items, generating means
.10/.90 and intercept SDs .4/.3, seed 37111 gives observed means .100694/.913194.
The candidate discarded a proposed coefficient update larger than 1.5 and
returned a zero effect with `converged = TRUE`. The final solver bounds that
primary step to size 1 and continues updating; nonfinite primary steps cannot
establish convergence. The example now converges at a conditional logit
contrast 4.8552 and marginal difference 80.922pp. An independent dense
conditional-mode solve and finite differences verify its marginal beta score.

In a second stress input setting every control response to zero and every
treatment response to one, the old candidate also falsely accepted a zero
effect. The guarded solver returns `max_iter_finite`, with finite coefficient
contrast 11.7939 but no valid marginal estimate/interval. Its covariance
score and parameter changes do not meet the stopping rule at 160 updates.
This failure is retained; a proper prior is not treated as a guarantee of
algorithmic convergence. These two inputs are numerical regression cases,
not a coverage study or evidence for unpenalized MLE existence under separation.
The independent beta/covariance and screening checks now contain 51 assertions.
Replaying all 24,576 confirmation datasets with the guarded solver gives
identical validity, status, coverage, detection, and joint decisions to the
pre-safeguard candidate. Common-valid point estimates and SEs are also
identical in the saved values (predeclared tolerance 1e-8 pp). This establishes
preservation on the ordinary confirmation data; the replay does not add
independent Monte Carlo replications.

[Complete safeguard replay audit](l2_uncertainty_comparison/beta_guard_replay.csv).

[Reproducible stress cases](l2_uncertainty_beta_guard.R) ·
[before safeguarding](l2_uncertainty_beta_before/results.csv) ·
[guarded results](l2_uncertainty_beta_guarded/results.csv).

**Exploratory bootstrap comparison completed.** All 1,536 outer datasets
have valid original, true-covariance-fixed, and bootstrap-SE results. The
bootstrap requested 399 refits per dataset, or 612,864 inner attempts. Across
the 12 cells, its RMS SE is between 0.998 and 1.013 times the original joint
SE. Coverage changes are not uniformly positive: for example, high-variance
162×30 null coverage is 92.19% with the original SE and 91.41% with the
bootstrap SE (bootstrap Wilson interval 85.27–95.13%). Conservative 162×60
planned coverage is 94.53% versus 96.09% (91.18–98.32%). With 128 outer
trials the nominal-coverage MCSE is about 1.9pp, so these are mechanism
checks, not a precision ranking. No universal SE multiplier, automatic
bootstrap remedy, or claim about the guarded solver's bootstrap intervals
is inferred. This exploratory series uses the former covariance estimator.

[Bootstrap summaries](l2_uncertainty_bootstrap/summary.csv) ·
[outer trial records](l2_uncertainty_bootstrap/fits.csv.gz).

**Screening must preserve the same numerical fit.** A focused check used the
first 64 fresh trials per design/condition, paired with the already completed
final-stage fits. The former early-stage settings (85 outer updates, 15 inner
steps, parameter tolerance 0.0006) produced only 732/768 valid fits, versus
766/768 in final confirmation. At 5pp/162×60, planned validity was 49/64 versus
64/64, and 18 joint-attainment decisions differed across the full probe.
The common valid estimates were close; the extra iteration-limit failures
were the material issue for candidate screening.

On the corrected zero-correlation, zero-guessing path, the R wrapper now
uses the final-stage numerical settings in screening too: 160 outer updates,
22 inner steps, tolerance 0.00015, with the same 0.005 scaled-score check.
It retains the existing smaller simulation budgets and skips bootstrap in
screening. No new iteration extension or relaxed score criterion was added
to final confirmation, whose 24,576 datasets therefore remain applicable.
The repeated probe has identical validity (766/768) and all 768 joint
decisions; common-valid estimate and SE differences are below 5e-14 pp
(the stored CSVs introduce decimal rounding). The predeclared numerical
comparison tolerances were 0.01pp for both quantities. This verifies the
entry-point correction on those data, not a new Monte Carlo guarantee or
an optimizer comparison over the whole grid.

[Screening runner](l2_uncertainty_stage1.R) ·
[before alignment](l2_uncertainty_stage1/summary.csv) ·
[after alignment](l2_uncertainty_stage1_fixed/summary.csv) ·
[paired trial records](l2_uncertainty_stage1_fixed/fits.csv.gz).

Reproducible runners are [saved-trial diagnostics](l2_uncertainty_diagnostics.R),
[profile calculations](l2_uncertainty_profile.R),
[independent confirmation](l2_uncertainty_confirmation.R),
[matched CPU and unchanged-path checks](l2_uncertainty_benchmark.R), and
[trial-level audit and plot](summarize_l2_uncertainty.R). Run each executable
in its own library and source directory; result manifests identify the
binary, source hashes, settings and seeds. Both errors and completed
checkpoints are retained. Runners refuse incompatible resume conditions.
The `candidate_*` directories retain earlier development variants. The
`confirmation_efficient` run precedes the large-step safeguard; the final
`confirmation_guarded` run reuses its same 24,576 response datasets and is
not pooled with it as additional independent replications. The 96-trial
profile comparison used the development candidate before that safeguard.
The final source is checked with separate beta/covariance score tests and
same-response confirmation. Preserved C++ and R-engine snapshots identify
the tested source; manifests identify the actual simulation binaries.
Reinstalled/archive binaries are separate builds, not asserted byte-identical.

The final source passed **713 regression assertions**, including 51 checks
of the covariance objective, screening consistency, and strong contrasts;
one optional Stan test was skipped. The archive's 37 R/C++/test files match
that tested source and the working tree byte-for-byte. Package build,
examples and vignette reconstruction passed on macOS/R 4.6.1. The archive
check had zero errors/warnings and one `xcrun_db` temporary-cache NOTE.
Its `--no-tests` check reuses the separately completed 713 assertions;
this combines explicitly matched stage evidence, not a new single-command
all-environment check. The optional `instantiate` package was unavailable.
Remote incoming/time checks were disabled; no new Stan, Windows/Linux,
or Win-builder run was performed. The pre-safeguard checks remain labeled
as intermediate evidence.

[Final regression log](l2_uncertainty_comparison/full-regression.log) ·
[final archive check](l2_uncertainty_comparison/package-check.log) ·
[archive identity and source matching](l2_uncertainty_comparison/archive-identity.json).

<a id="default-bootstrap-confirmation"></a>

## Default-bootstrap confirmation before submission

The final covariance-corrected solver was checked on its default
probability-scale interval path. The preceding 24,576-dataset study used
`joint_delta`; the 399-refit exploratory series used the previous solver.
Neither is relabeled as validation of the corrected default bootstrap.

[Runner](l2_default_bootstrap.R) · [independent aggregation](summarize_l2_default_bootstrap.R).
Three frozen profiles use 54×30 (intercepts only), 72×60 (moderate slopes),
and 162×30 (larger variances), each with 30 words per person, no guessing,
zero assumed correlations, and a .42/.50 planned mean or equal-marginal null.
Each of the six cells has 512 new outer datasets. Response and inner RNG
roots are 941260001 and 942260001. The default uses 99 inner refits; the
first 32 datasets per cell use 399, with their first 99 draws supplying the
paired default comparison. Two fixed probes per cell additionally use the
unmodified 99-refit entry point. The function bodies and installed namespace
are unchanged; a local R environment retains the returned C++ vector.

Detection and detection/width attainment use every attempted outer dataset;
coverage uses valid fits. Inner failures and the existing >=90% inner-success
rule are retained explicitly. Prespecified practical checks remain observed
valid rate >=99%, Wilson lower coverage bound >=92.5%, and upper null bound
<=7.5%; joint planning uses detection AND width <=10pp, lower bound >=80%.
At 512 replications the nominal-null MCSE is about 0.96pp. Indeterminate
decisions are retained without result-dependent extensions. The 32-pair
99/399 check measures inner simulation sensitivity, not comparative coverage
with precise outer Monte Carlo uncertainty.

**Completed confirmation.** All 3,072 outer datasets returned valid fits
under both interval methods. All 361,728 retained inner refits were valid;
the unmodified-call probes requested another 1,188 refits. Those probes
matched the captured 99-refit point estimates, SEs and validity decisions.
All paired response hashes, trial keys, denominators and retained-vector
SEs passed independent aggregation checks. The nearest coverage, detection
and width decision boundaries were respectively 0.000961, 0.000163 and
0.002101 pp away, exceeding the 1e-7 pp matching tolerance.

| Profile; N×J | Interval method | Null rejection (95% Wilson MC interval) | Planned coverage (95% MC interval) | Detection AND width <=10pp (95% MC interval) |
| --- | --- | --- | --- | --- |
| Intercepts only; 54×30 | Joint delta | 7.23% (5.29–9.80) | 96.09% (94.04–97.46) | 95.51% (93.35–96.95) |
| Intercepts only; 54×30 | Default bootstrap, 99 | 6.25% (4.46–8.69) | 94.92% (92.66–96.51) | 91.99% (89.32–94.04) |
| Moderate slopes; 72×60 | Joint delta | 3.91% (2.54–5.96) | 95.70% (93.58–97.15) | 84.77% (81.39–87.62) |
| Moderate slopes; 72×60 | Default bootstrap, 99 | 3.71% (2.39–5.72) | 95.51% (93.35–96.99) | 77.15% (73.32–80.58) |
| Larger variances; 162×30 | Joint delta | 6.84% (4.96–9.36) | 94.14% (91.76–95.87) | 92.19% (89.54–94.21) |
| Larger variances; 162×30 | Default bootstrap, 99 | 7.81% (5.79–10.46) | 93.75% (91.31–95.54) | 90.23% (87.36–92.51) |

Every row uses 512 datasets per condition. Bootstrap joint attainment passes
the 80% lower-bound gate at 54×30 and 162×30; 72×60 is **indeterminate**, not
a confirmed pass or confirmed failure. Its paired change from joint delta
is -7.62pp (MCSE 1.68pp), with 77/512 individual joint decisions changing.
Its RMS SE ratio is 0.998, and its mean width is 9.20 versus 9.24pp. Thus
similar average SEs or widths do not establish interchangeable threshold
decisions. The per-dataset bootstrap/joint-delta SE ratio spans 0.863–1.130
at the empirical 2.5th–97.5th percentiles in this cell.

The default's practical coverage/null checks both pass at 72×60. At 54×30,
planned coverage passes but the null check is indeterminate. At 162×30,
planned coverage and the null check are indeterminate; its null interval
lies above nominal 5%, an adverse calibration flag. Joint delta also has a
nominal-null flag at 54×30 in this independent batch. These results do not
supersede the earlier 2,048-per-cell results or justify claiming exact
95%/5% calibration from all successful fits.

On the 32-dataset planned subsets, increasing 99 to 399 refits changes
joint attainment from 31/32 to 31/32, 22/32 to 26/32, and 27/32 to 30/32
for 54×30, 72×60 and 162×30 respectively. Individual joint decisions change
in 2, 6 and 3 pairs. This is sensitivity evidence, with insufficient outer
precision to establish 399-refit coverage or a new default. The public
defaults remain unchanged. A final planning confirmation must use the
interval method and inner-refit count intended for the analysis; a
joint-delta pass is not transferred to the default bootstrap. The three
selected layouts neither establish new minimum sample sizes nor provide
uniform inferential calibration.

[Trial records](l2_default_bootstrap_confirm/fits.csv.gz) ·
[summary](l2_default_bootstrap_confirm/summary.csv) ·
[paired differences and MCSEs](l2_default_bootstrap_confirm/paired.csv) ·
[actual inner counts](l2_default_bootstrap_confirm/inner_execution.csv) ·
[comparison figure](l2_default_bootstrap_confirm/comparison.png) ·
[aggregation audit](l2_default_bootstrap_confirm/summary.log).

The separate technical smoke run used roots 941250001/942250001 and four
outer datasets per cell. All 24 outer datasets and 5,976 retained inner
refits were valid; the unmodified-call probes agreed within 1e-7 pp.
These small counts establish operation and matching, not calibration.
[Smoke records](l2_default_bootstrap_smoke/fits.csv.gz) ·
[inner counts](l2_default_bootstrap_smoke/inner_execution.csv).

[Public-API workflow checks](release_workflows.R) additionally exercise
default-bootstrap full/partial GLMM planning and fresh-design confirmation,
random-item Rasch recovery and fixed-bank MML. Tables, plots and reports
survive save/read round trips, and the partial-presented GLMM numeric
results agree with one and two configured workers. The checks use six
replications per design for integration only; their passes, failures and
indeterminate decisions are not sample-size recommendations.
[Workflow check log](release_workflows/check.log).

The resulting source archive is 1,254,503 bytes, SHA-256
`b9b025ca54039867b32795661b63c20f02f9f50631b72c98e19002eb07b31087`.
Its 37 R/C++/test files match the source that passed 713 regression
assertions (one opt-in Stan skip). The local final difference check has
zero errors/warnings and one incoming NOTE for the new submission and
`cmdstanr` additional repository. Remote URL checking is enabled and
passes. Examples and vignette reconstruction passed on the immediately
preceding archive; the final build also reconstructed vignettes. The final
archive changes only packaging metadata, four validation-guide URLs and
their HTML, and three printed environment addresses in the otherwise
identical workflow HTML. This is combined stage evidence, not one full
check command. Win-builder results for this source and the preceding
pre-link-fix archive have now been retrieved; see the results below.
[Archive identity](release_candidate/archive-identity.json) ·
[final difference check](release_candidate/package-check.log) ·
[preceding example/vignette check](release_candidate/pre-link-fix-check.log).

## Plot defaults and presentation checks

GLMM/Rasch plotting and `check_dgp()` now use `theme_bw()` and omit automatic
titles, subtitles and captions by default. `show_annotations = TRUE` restores
them; axis/legend/facet labels and data remain present. The 3D plot uses a
white background. The [presentation rationale and references](../inst/validation-guide.md#visualization)
draw on Morris et al. (2019), Rücker and Schwarzer (2014), and rsimsum's
documentation; no estimator or simulation criterion changes.

[Display checks](plot_defaults.R) reuse the saved release-workflow objects.
A four-replication scenario fixture is for display checks only. Ten static
plot variants render with and without generated annotations, with identical
built layer data; custom titles can be added with `labs()`. The Plotly layout
and trace data are also checked. Saved input objects remain unchanged.
The current default-bootstrap figure is regenerated from the same saved
trials and summaries; its numerical results are unchanged. Earlier archived
figures retain their original presentation.
[Targeted test log](plot_defaults/tests.log) · [render checks](plot_defaults/render.log) ·
[validation example](plot_defaults/rates.png) · [scenario example](plot_defaults/scenarios.png) ·
[Rasch example](plot_defaults/rasch_rmse.png).

The targeted plot/plan-method tests pass 79 assertions without failures,
warnings or skips. The local `--as-cran --no-manual --no-tests` check passes
installation, examples, help and vignette reconstruction with 0 errors,
0 warnings and the known incoming NOTE (new submission / optional `cmdstanr`).
After that check, one section anchor was moved to its correct heading and
the archive rebuilt; code, tests and help are identical to the checked archive.
All 13 README-relative links resolve in the final archive. These are combined
stage results; the preceding three-platform CI and completed Win-builder checks
do not include the new plot defaults. See the [package check](plot_defaults/package-check.log),
[build log](plot_defaults/build.log), [archive identity](plot_defaults/archive-identity.json)
and [source comparison](plot_defaults/source-comparison.txt).

## Win-builder results received on 2026-11-26

These checks predate the plot-default changes. Both use Windows Server 2022,
x86_64-w64-mingw32, and GCC 14.3.0.

| Environment | Submitted archive SHA-256 prefix | Result | Tests |
| --- | --- | --- | --- |
| R-devel, 2026-11-17 r90926 ucrt | `b9b025ca...` | 0 errors / 0 warnings / 1 NOTE | 675 pass / 0 fail / 0 warn / 3 skip |
| R-release, 4.6.2 (2026-10-29 ucrt) | `219ebc6c...` | 0 errors / 0 warnings / 1 NOTE | 675 pass / 0 fail / 0 warn / 3 skip |

The skips are two CRAN-only exclusions and the opt-in Stan fit. Examples,
vignette reconstruction, and PDF/HTML manuals pass on both. The common incoming
NOTE concerns the new submission, spelling candidates (author names and
statistical terms), and the optional `cmdstanr` additional repository. That
repository is available, but `cmdstanr` was not installed for checking (INFO).
R-release also flags four URLs: two noncanonical CRAN links and two unpublished
figure links. These were already corrected in the R-devel and current
archives; R-devel and the current local check have no URL findings.

The logs were retrieved before expiry. Submitted archive hashes, binary
DESCRIPTION fields and the bundled guide were matched, allowing newline and
installation metadata differences. Raw installation logs are losslessly
compressed; their warnings come from RcppEigen headers, with no compiler errors.
This is package-operation evidence, not additional simulation calibration.
No tests or simulations were rerun solely to record these results.

[R-devel check](release_candidate/win-builder/devel/00check.log) ·
[R-devel tests](release_candidate/win-builder/devel/testthat.Rout) ·
[R-release check](release_candidate/win-builder/release/00check.log) ·
[R-release tests](release_candidate/win-builder/release/testthat.Rout) ·
[URLs, archive identities and file hashes](release_candidate/win-builder/results.json).

## Reference Environment & Isolated Reproduction

All recent `glmmplanr 0.1.0` confirmation runs were executed under:

- **R**: `4.6.1` (`aarch64-apple-darwin23`, as recorded in the saved session information)
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
sh validation/run-isolated.sh sparse_and_testlet_recheck
sh validation/run-isolated.sh testlet_allocation
sh validation/run-isolated.sh l2_vocabulary_allocation
sh validation/run-isolated.sh l2_vocabulary_covariates
sh validation/run-isolated.sh l2_vocabulary_sparse
sh validation/run-isolated.sh l2_vocabulary_prior_knowledge
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
| Strict conditional null (`0.000`) | `72 x 32` | 192/192 | 186/192 (`96.9% [93.3, 98.8]`) | 192/192 | 6/192 (`3.1% [1.2, 6.7]`) | 0.155 / 0.168 | pass |
| Strict conditional null (`0.000`) | `48 x 48` | 192/192 | 184/192 (`95.8% [92.0, 98.2]`) | 192/192 | 8/192 (`4.2% [1.8, 8.0]`) | 0.164 / 0.165 | pass |
| Equal-marginal null (`+0.0041`) | `72 x 32` | 192/192 | 187/192 (`97.4% [94.0, 99.1]`) | 192/192 | 6/192 (`3.1% [1.2, 6.7]`) | 0.151 / 0.165 | pass |
| Equal-marginal null (`+0.0041`) | `48 x 48` | 192/192 | 182/192 (`94.8% [90.6, 97.5]`) | 192/192 | 10/192 (`5.2% [2.5, 9.4]`) | 0.163 / 0.165 | pass |

### 2B. High-$N$ Strict Conditional-Null Calibration (`N = 512` per design, $\beta = 0$)

Pre-specified calibration gates: exact 95% Type I error upper bound $\le 7.5\%$ (denominator = all `512` attempted datasets) and exact 95% coverage lower bound $\ge 92.5\%$ (denominator = valid fits).

| Design | Valid | 95% Coverage (Exact 95% CI) | Coverage Gate ($\ge 92.5\%$) | Type I Error (Exact 95% CI) | Type I Gate ($\le 7.5\%$) | Emp. SD / Mean SE |
|---|---:|---:|---|---:|---|---:|
| `48 x 48` | 512/512 | 495/512 (`96.7% [94.7, 98.1]`) | **pass** | 17/512 (`3.32% [1.95, 5.26]`) | **pass** | 0.158 / 0.166 |
| `72 x 32` | 511/512 | 493/511 (`96.5% [94.5, 97.9]`) | **pass** | 18/512 (`3.52% [2.10, 5.50]`) | **pass** | 0.153 / 0.167 |

### 2C. Necessity / Pass-Fail Boundary & Multi-Seed Optimizer Search (`bayes_laplace`)

Six smaller adjacent candidates evaluated on `192` fresh datasets each (planned effect `+1.616` logit, width target $\le 1.00$ logit with 80% assurance floor, seeds `7611001–7616192`):

| Design | Trials | Valid | 95% Coverage (Exact 95% CI) | Width $\le 1.00$ (Exact 95% CI) | Mean Width | Width Status |
|---|---:|---:|---:|---:|---:|---|
| `32 x 32` | 1,024 | 191/192 | 185/191 (`96.9% [93.3, 98.8]`) | 10/192 (`5.2% [2.5, 9.4]`) | 1.150 | **below_floor** |
| `48 x 24` | 1,152 | 192/192 | 186/192 (`96.9% [93.3, 98.8]`) | 25/192 (`13.0% [8.6, 18.6]`) | 1.090 | **below_floor** |
| `48 x 32` | 1,536 | 192/192 | 188/192 (`97.9% [94.8, 99.4]`) | 163/192 (`84.9% [79.0, 89.6]`) | 0.928 | **indeterminate** |
| `64 x 24` | 1,536 | 192/192 | 187/192 (`97.4% [94.0, 99.1]`) | 163/192 (`84.9% [79.0, 89.6]`) | 0.930 | **indeterminate** |
| `72 x 24` | 1,728 | 192/192 | 192/192 (`100.0% [98.1, 100.0]`) | 187/192 (`97.4% [94.0, 99.1]`) | 0.882 | **pass** |
| `48 x 40` | 1,920 | 192/192 | 185/192 (`96.4% [92.6, 98.5]`) | 191/192 (`99.5% [97.1, 100.0]`) | 0.823 | **pass** |

Across 5 independent search seeds (`optimize_sample_size()`, `method = "surrogate_gp"`, `prefer = "fewer_trials"`, `nsim_final = 96`, lattice `{32,48,64,72} x {24,32,40,48}`, seeds `7650001–7650005`), **5/5 seeds** reached an internal Stage-2 `pass`, selecting `72 x 24` (`1,728` trials, Stage-2 width attainment `92.7%–96.9%`) on **4/5 seeds** and `48 x 32` (`1,536` trials, Stage-2 width attainment `88.5%`) on **1/5 seed**. Cross-checked against the independent 192-rep boundary table above, **4/5** selections (`72 x 24`) confirmed as independent `pass`, **1/5** (`48 x 32`) was `indeterminate` on the boundary, and **0/5** fell `below_floor`.

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

This historical study retains its original results and differing confirmation
procedures. For the current comparison using equal actual simulation budgets
and a shared independent confirmation rule, see the
[matched item-pool comparison](#item-pool-package-comparison) above. Do not
combine the two studies as a single ranking of search algorithms.

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
| [`partial_presentation.R`](partial_presentation.R) | [`partial_presentation_summary.csv`](partial_presentation_summary.csv), [`partial_presentation_fits.csv`](partial_presentation_fits.csv) | Initial 32-rep pilot: `64 x 20` (5 items/person), random intercepts, `30/32` coverage (exact 95% CI `79.2–99.2%`). |
| [`sparse_and_testlet_recheck.R`](sparse_and_testlet_recheck.R) | [`sparse_and_testlet_recheck_summary.csv`](sparse_and_testlet_recheck_summary.csv), [`sparse_and_testlet_recheck_fits.csv`](sparse_and_testlet_recheck_fits.csv) | Pre-specified 96-rep confirmation across sparse partial presentation (intercepts & correlated slopes) and testlet nesting (see Section 4A below). |
| [`links_covariates.R`](links_covariates.R), [`probit_null.R`](probit_null.R) | [`link_summary.csv`](link_summary.csv), [`probit_null_summary.csv`](probit_null_summary.csv) | External `lme4` logit vs. probit link and correlated covariate sensitivity checks (`959/960` valid; probit null recheck `5/64`–`7/64`). |
| [`extreme_patterns.R`](extreme_patterns.R) | [`extreme_summary.csv`](extreme_summary.csv), [`extreme_frequency.csv`](extreme_frequency.csv), [`extreme_iterations.csv`](extreme_iterations.csv) | All-correct / all-incorrect unit frequency and Armijo-backtracking convergence checks under `chance = 1/3`. |
| [`native_interactions.R`](native_interactions.R), [`trial_covariates.R`](trial_covariates.R), [`covariate_allocation.R`](covariate_allocation.R), [`within_between.R`](within_between.R) | [`native_results.csv`](native_results.csv), [`trial_covariate_summary.csv`](trial_covariate_summary.csv), [`covariate_allocation_summary.csv`](covariate_allocation_summary.csv) | Exploratory interaction, trial-level covariate, and within/between allocation benchmarks. |

### 4A. Sparse Partial Presentation & Testlet Nesting Recheck (`N = 96` per cell)

- **Script & CSVs**: [`sparse_and_testlet_recheck.R`](sparse_and_testlet_recheck.R) ([`sparse_and_testlet_recheck_summary.csv`](sparse_and_testlet_recheck_summary.csv), [`sparse_and_testlet_recheck_fits.csv`](sparse_and_testlet_recheck_fits.csv))
- **Pre-specified gates**: `96` independent outer datasets per cell (`99` parametric bootstrap refits where applicable), validity $\ge 90\%$ (`87/96`), exact binomial 95% CI lower bound on valid-fit 95% coverage $\ge 90\%$ (`pass`), upper bound $< 90\%$ (`below_floor`), otherwise `indeterminate`; two-sided 5% null detection uses all `96` scheduled datasets in the denominator. Regularized `bayes_laplace` records boundary singularity as `NA` (not applicable).

| Regime (Trials; Seeds) | Condition | Estimator (`scale`) | Valid (Sing.) | 95% Coverage (Exact 95% CI) | Detection (Exact 95% CI) | Emp. SD / Mean SE | Floor Status |
|---|---|---|---:|---:|---:|---:|---|
| `partial_intercept_64x20_k5` (`320`; `30001–30096`) | Null (`0 pp`) | `bayes_laplace` boot99 (`pp`) | 96/96 (`NA`) | 90/96 (`93.8% [86.9, 97.7]`) | 6/96 (`6.2% [2.3, 13.1]`) | 5.32 / 5.91 | indeterminate |
| `partial_intercept_64x20_k5` (`320`; `30001–30096`) | Null (`0 pp`) | `bayes_laplace` joint\_delta (`pp`) | 96/96 (`NA`) | 91/96 (`94.8% [88.3, 98.3]`) | 5/96 (`5.2% [1.7, 11.7]`) | 5.32 / 5.90 | indeterminate |
| `partial_intercept_64x20_k5` (`320`; `30001–30096`) | Null (`0 pp`) | `lme4::glmer` boot99 (`pp`) | 96/96 (`43`) | 91/96 (`94.8% [88.3, 98.3]`) | 5/96 (`5.2% [1.7, 11.7]`) | 5.37 / 5.85 | indeterminate |
| `partial_intercept_64x20_k5` (`320`; `31001–31096`) | Planned (`-5 pp`) | `bayes_laplace` boot99 (`pp`) | 96/96 (`NA`) | 89/96 (`92.7% [85.6, 97.0]`) | 10/96 (`10.4% [5.1, 18.3]`) | 5.76 / 5.85 | indeterminate |
| `partial_intercept_64x20_k5` (`320`; `31001–31096`) | Planned (`-5 pp`) | `bayes_laplace` joint\_delta (`pp`) | 96/96 (`NA`) | 89/96 (`92.7% [85.6, 97.0]`) | 10/96 (`10.4% [5.1, 18.3]`) | 5.76 / 5.95 | indeterminate |
| `partial_intercept_64x20_k5` (`320`; `31001–31096`) | Planned (`-5 pp`) | `lme4::glmer` boot99 (`pp`) | 96/96 (`48`) | 88/96 (`91.7% [84.2, 96.3]`) | 10/96 (`10.4% [5.1, 18.3]`) | 5.81 / 5.86 | indeterminate |
| `partial_slopes_64x24_k8` (`512`; `40001–40096`) | Null (`0 pp`) | `bayes_laplace` boot99 (`pp`) | 96/96 (`NA`) | 94/96 (`97.9% [92.7, 99.7]`) | 2/96 (`2.1% [0.3, 7.3]`) | 4.10 / 4.74 | **pass** |
| `partial_slopes_64x24_k8` (`512`; `40001–40096`) | Null (`0 pp`) | `bayes_laplace` logit (`logit`) | 96/96 (`NA`) | 95/96 (`99.0% [94.3, 100.0]`) | 2/96 (`2.1% [0.3, 7.3]`) | 0.204 / 0.239 | **pass** |
| `partial_slopes_64x24_k8` (`512`; `41001–41096`) | Planned (`+10 pp`) | `bayes_laplace` boot99 (`pp`) | 96/96 (`NA`) | 95/96 (`99.0% [94.3, 100.0]`) | 54/96 (`56.2% [45.7, 66.4]`) | 4.02 / 4.63 | **pass** |
| `partial_slopes_64x24_k8` (`512`; `41001–41096`) | Planned (`+10 pp`) | `bayes_laplace` logit (`logit`) | 96/96 (`NA`) | 94/96 (`97.9% [92.7, 99.7]`) | 58/96 (`60.4% [49.9, 70.3]`) | 0.216 / 0.245 | **pass** |
| `testlet_matched_48x24_pt0` (`1,152`, 6 passages; `50001–50096`) | Null (`0 pp`) | `bayes_laplace` logit (`logit`) | 96/96 (`NA`) | 82/96 (`85.4% [76.7, 91.8]`) | 14/96 (`14.6% [8.2, 23.3]`) | 0.539 / 0.397 | indeterminate |
| `testlet_matched_48x24_pt0` (`1,152`, 6 passages; `50001–50096`) | Null (`0 pp`) | `lme4::glmer` logit (`logit`) | 96/96 (`11`) | 76/96 (`79.2% [69.7, 86.8]`) | 20/96 (`20.8% [13.2, 30.3]`) | 0.560 / 0.372 | **below_floor** |
| `testlet_matched_48x24_pt0` (`1,152`, 6 passages; `51001–51096`) | Planned (`+10 pp`) | `bayes_laplace` joint\_delta (`pp`) | 96/96 (`NA`) | 91/96 (`94.8% [88.3, 98.3]`) | 25/96 (`26.0% [17.6, 36.0]`) | 6.98 / 7.77 | indeterminate |
| `testlet_localdep_48x24_pt04` (`1,152`, 6 passages; `52001–52096`) | Null (`0 pp`) | `bayes_laplace` logit (`logit`) | 96/96 (`NA`) | 81/96 (`84.4% [75.5, 91.0]`) | 15/96 (`15.6% [9.0, 24.5]`) | 0.476 / 0.383 | indeterminate |
| `testlet_localdep_48x24_pt04` (`1,152`, 6 passages; `52001–52096`) | Null (`0 pp`) | `lme4::glmer` logit (`logit`) | 96/96 (`13`) | 79/96 (`82.3% [73.2, 89.3]`) | 17/96 (`17.7% [10.7, 26.8]`) | 0.489 / 0.356 | **below_floor** |
| `testlet_localdep_48x24_pt04` (`1,152`, 6 passages; `53001–53096`) | Planned (`+10 pp`) | `bayes_laplace` joint\_delta (`pp`) | 96/96 (`NA`) | 77/96 (`80.2% [70.8, 87.6]`) | 41/96 (`42.7% [32.7, 53.2]`) | 9.41 / 7.50 | **below_floor** |
| `testlet_localdep_48x24_pt04` (`1,152`, 6 passages; `53001–53096`) | Planned (`+10 pp`) | `lme4::glmer` logit (`logit`) | 96/96 (`7`) | 73/96 (`76.0% [66.3, 84.2]`) | 46/96 (`47.9% [37.6, 58.4]`) | 0.499 / 0.365 | **below_floor** |

Key takeaways: (1) Under `64 x 20` between-subject partial presentation (`5` items/person), `bayes_laplace` (`89/96 = 92.7%` planned; `90–91/96` null) closely matched `lme4::glmer` (`88/96 = 91.7%`, with `48/96` singular ML fits), narrowing the 32-rep pilot CI (`79.2–99.2%`) to `85.6–97.0%` (`indeterminate` at $N = 96$ relative to the 90% floor). (2) Within-unit `64 x 24` Latin-square partial presentation (`8` items/person) with correlated random slopes passed the 90% floor under `bayes_laplace` (`94–95/96`, slightly conservative SEs). (3) In between-passage testlet assignment with only `6` passages (`3` per condition), asymptotic normal Wald intervals undercovered for both `lme4::glmer` (`20.8%` null Type I error, `below_floor`) and `bayes_laplace` (`14.6%`), and unmodeled participant-by-passage local dependence (`SD = 0.4`) pushed planned-effect coverage `below_floor` (`76.0–82.3%`) despite `1,152` total trials.

### 4B. Passage allocation at fixed response count

The planned comparison holds 48 readers, 24 questions per reader, and 1,152 responses fixed. Between-passage designs have 24×1, 12×2, 8×3, 6×4, 4×6, or 2×12 passages × questions, with an equal number of passages per arm. The 2×12 design is a boundary case. A separate 3×8 design assigns 24 readers to each arm and gives both arms the same three passages. All conditions fully present the assigned task, use binary logit (`chance = 0`), and estimate the population-averaged difference in response probabilities.

The generating model sums independent, centered normal reader, passage, item-within-passage, and reader-by-passage effects on the logit scale. Reader SD is 0.7 in the main comparison. The raw marginal item-difficulty variance is 0.41 (`item_sd = sqrt(0.41 - passage_sd^2)`). Cross passage SD 0 or 0.5 with reader-by-passage SD 0 or 0.4. The fifth profile uses passage SD 0.5 and reader-by-passage SD `0.4 * sqrt((questions - 1) / 3)`. At 24×1 the passage and item effects reduce to one item intercept, and the reader-by-passage effect is an independent contribution for each answer. The generator retains its variance under the constant-SD profiles. For 3×8, the treatment coefficient is constant across passages.

Both arms have true marginal probability 0.60 under the null and the planned probabilities are 0.60 and 0.70. An independent ordinary-R integral calibrates the conditional intercepts and checks the target probabilities before any fit. The 6×4 and 12×2 sensitivity profiles vary passage SD (0.25 or 0.60), reader-by-passage SD (0.20 or 0.60), reader SD (0.40 or 1.00), baseline probability (0.80 with planned probability 0.90), or effect (0.05). The joint high profile sets passage and reader-by-passage SDs to 0.60. Each sensitivity profile uses the remaining values from the constant-dependence reference. The smaller effect reuses the reference null because their generating distributions coincide.

The 106 design/profile/outcome combinations reduce to 96 distinct cells once equivalent DGPs are identified. Each cell generates 1,024 independent datasets with master seed `83251001`. Null and planned outcomes use separate L'Ecuyer-CMRG streams. Configurations within an outcome use matched draws for readers, item and passage positions within each arm, reader-by-passage effects, and response uniforms wherever the layouts permit. The analysis uses the built-in `bayes_laplace` engine with intercepts and `joint_delta` 95% intervals (zero bootstrap refits). Its passage intercept captures shared difficulty. The reader-by-passage effect enters the generator only.

Use exact two-sided 95% binomial intervals pointwise for the simulation rates. Coverage uses valid datasets and passes when its lower bound is at least 90%; an upper bound below 90% or a valid rate below 90% gives a failure. Return indeterminate in the remaining cases. Detection uses all scheduled datasets. The null gates for pass and failure are an upper bound at most 7.5% and a lower bound above 7.5%, respectively. Compare planned detection with the package's default 80% power target using its exact Monte Carlo bounds. Report bias, empirical SD, mean SE, and mean width on the percentage-point scale alongside the rate diagnostics.

The isolated run under R 4.6.1 used the estimator source at `bcd5ee0`. All **98,304** simulated datasets yielded valid estimates and intervals. The [scenario table](testlet_allocation_scenarios.csv) maps equivalent profiles to the same cell. The [summary](testlet_allocation_summary.csv) gives exact Monte Carlo bounds for all distinct cells. The [compressed trial data](testlet_allocation_fits.csv.gz) and [reproduction script](testlet_allocation.R) provide replication-level results and diagnostics.

| Allocation and layout | Null detection (95% MC CI) | Planned 95% coverage (95% MC CI) | Planned power | Mean full width (pp) |
|---|---:|---:|---:|---:|
| Between passages, 24×1 | 67/1,024 (6.5% [5.1, 8.2]) | 949/1,024 (92.7% [90.9, 94.2]) | 48.3% | 20.34 |
| Between passages, 12×2 | 82/1,024 (8.0% [6.4, 9.8]) | 947/1,024 (92.5% [90.7, 94.0]) | 37.8% | 24.67 |
| Between passages, 8×3 | 116/1,024 (11.3% [9.5, 13.4]) | 933/1,024 (91.1% [89.2, 92.8]) | 32.0% | 27.49 |
| Between passages, 6×4 | 130/1,024 (12.7% [10.7, 14.9]) | 897/1,024 (87.6% [85.4, 89.6]) | 30.9% | 29.39 |
| Between passages, 4×6 | 181/1,024 (17.7% [15.4, 20.2]) | 858/1,024 (83.8% [81.4, 86.0]) | 31.9% | 31.31 |
| Between passages, 2×12 | 420/1,024 (41.0% [38.0, 44.1]) | 622/1,024 (60.7% [57.7, 63.7]) | 49.3% | 23.53 |
| Between readers, 3×8 | 60/1,024 (5.9% [4.5, 7.478]) | 972/1,024 (94.9% [93.4, 96.2]) | 59.9% | 17.97 |

The table uses the constant-dependence profile. In the paired 12×2 minus 6×4 comparison, null detection fell by 4.69 pp (approximate 95% paired MC CI −6.97 to −2.41), planned coverage increased by 4.88 pp (2.75 to 7.01), and mean interval width fell by 4.72 pp (−5.17 to −4.27). The 12×2 design met the 90% coverage floor under both effects and had an indeterminate null detection gate. The 2×12 boundary had a smaller fitted SE than the 4×6 layout and markedly greater empirical dispersion. Its planned empirical SD was 13.62 pp and mean SE was 6.00 pp.

Shared passage difficulty mattered more for error rates in this complete-presentation comparison. At 6×4 with reader-by-passage SD 0.4, passage SDs of 0, 0.25, 0.50, and 0.60 gave null detection of 2.7%, 5.8%, 12.7%, and 14.9%, respectively. Holding passage SD at 0.50 and varying the reader-by-passage SD from 0 to 0.60 gave 12.4–13.3% null detection. Both 6×4 and 12×2 met the coverage and null gates at passage SD 0.25 with interaction SD 0.4. With both SDs at 0.60, null detection was 14.7% at 6×4 and 9.8% at 12×2, with both intervals entirely above the 7.5% upper tolerance. Size-dependent reader-by-passage variation gave similar decisions to constant dependence for the between-passage layouts tested. In the independence control, the matched data are identical across the between-passage layouts, and the specified passage random effects account for differences across their fits.

The 3×8 reader-randomized design met the calibration gates under constant and size-dependent local variation with a common treatment effect across the three shared passages. A passage-randomized design has a different unit of assignment. Across the 49 distinct planned cells, the highest power was 710/1,024 (69.3%, 95% MC CI 66.4–72.1%) for 12×2 at 80% versus 90%. All planned cells fell below the 80% power target. These results apply to the specified allocation, full presentation, effect scale, and generating models.

### 4C. L2 vocabulary recall and the allocation of participants and words

The planning target is the population-averaged difference in correct meaning recall between retrieval practice and restudy. The generator uses one binary response per participant–word pair after a fixed delay, full presentation, and a standard logit link. The null uses marginal recall of 42% in both arms. The planned probabilities are 42% and 50%, based on the three-round computer-test means in [Huang et al. (2025), Table 3](https://www.mdpi.com/2076-328X/15/11/1540). The 42% versus 47% condition probes a smaller effect. Calibrating each generating model to these probabilities holds the estimand fixed across allocations and variance profiles. The benchmark checks calibration against a separate normal integral.

Compare 32×48, 48×32, 64×24, and 96×16 at 1,536 responses and 48×64, 64×48, and 96×32 at 3,072 responses. A two-list Latin square places each word and each participant in both conditions across the full schedule. The participant-randomized design assigns each participant to one arm and shares every word across arms. A sensitivity check assigns words to arms with exchangeable word sampling and presents all words to every participant (32×48 and 96×16).

Use centered unit-step condition coding (−0.5, +0.5), Gaussian intercept SDs 0.8 for participants and 0.6 for words, and the following participant/word condition-slope SDs. The core Latin square uses 0/0, 0.25/0.25, 0.40/0.15, and 0.15/0.40. Participant assignment uses 0/0, 0/0.25, and 0/0.40. Word assignment uses 0.25/0. With slope SD 0.25 in all eligible groups, sensitivity at 32×48 and 96×16 changes one intercept–slope correlation at a time from zero to ±0.5. The full allocation grid examines the five-point effect with slope SD 0.25 in each eligible group. The SDs and correlations are planning assumptions. For an example of reported estimates in real vocabulary data, see [Wang & Yum (2025), Tables 3–4](https://doi.org/10.1371/journal.pone.0334038).

The allocation run uses 140 DGPs × 2,048 datasets (master seed 85321001). Independent L'Ecuyer-CMRG replication streams provide matching underlying normal draws and uniform draws for the same IDs across scenarios. The primary analysis uses the native regularized `bayes_laplace` estimator, allowing the slopes and correlations supported by the allocation, including when true slope variances are zero. It uses joint-delta 95% intervals on the percentage-point scale, with zero bootstrap refits. At 32×48 and 96×16 the balanced-slope profiles for the null and eight-point effect also use zero-correlation analyses on the same responses. The matching zero-correlation DGPs at those two layouts also use random-intercept analyses. These extra fits bring the total to 368,640. The script preserves completed scenarios in a checkpoint directory under its output directory. [Allocation script](l2_vocabulary_allocation.R)

Report exact pointwise 95% Monte Carlo intervals for coverage among valid fits and for detection among all scheduled fits. Use the allocation benchmark gates in §4B for validity, coverage, null detection, and 80% power. Report mean and 80th percentile of full interval width, bias, empirical SD, and mean SE in percentage points. The comparison concerns candidate schedules and selected DGPs, with no absolute width requirement.

The independent covariate check uses 1,024 datasets per DGP (master seed 85323001) at 32×48 and 96×16. Each participant has a continuous standard-normal baseline vocabulary score and each word an independent standard-normal proxy for centered, standardized length. Their fixed logit effects are +0.50 and −0.35. The residual intercept variances are reduced by 0.50² and 0.35², respectively, so the corresponding total intercept variances match the core study. For the Latin square, a second DGP adds a +0.35 condition × vocabulary logit interaction per vocabulary SD. The same simulated responses are fitted with the unadjusted and adjusted models, and with the interaction when present. All analyses target the same population-average contrast, including standardization to the known covariate distribution in the adjusted models. The script passes continuous columns to the existing native kernel and checks the standardized generating contrast against the aggregate-normal calibration. This research interface uses internal functions. The public planning interface accepts fixed-effect factors. [Covariate script](l2_vocabulary_covariates.R)

#### Participants versus words

The isolated macOS run under R 4.6.1 used the estimator source at `bcd5ee0`. All 368,640 fits across 286,720 datasets were valid. Recomputing coverage, detection, width, and the exact Monte Carlo gates from the saved estimates reproduced the summary. The main run took about 9 minutes 33 seconds with four workers, including result serialization. The covariate calculation used the same installed package and took about 1 minute 45 seconds with two workers, overlapping the main run. [Scenarios](l2_vocabulary_scenarios.csv) · [Summary](l2_vocabulary_summary.csv) · [Replication-level results](l2_vocabulary_fits.csv.gz) · [Interval-width figure](l2_vocabulary_allocation.png) · [Plot script](plot_l2_vocabulary_allocation.R)

At 1,536 responses, the two-list designs reversed their relative performance when the location of slope variation changed. The table shows results for an eight-point effect, with interval widths on the percentage-point scale.

| Generating slope SDs, participant/word | 32×48 width | 96×16 width | 32×48 power | 96×16 power |
|---|---:|---:|---:|---:|
| Participant-dominant, 0.40/0.15 | 11.66 | 11.36 | 77.6% | 81.9% |
| Word-dominant, 0.15/0.40 | 11.38 | 12.60 | 81.0% | 71.8% |

At 1,536 responses, 48×32 was the only layout that passed all gates for the balanced slopes (power 82.4%, 95% MC CI 80.7–84.0%, coverage 96.4%, null detection 3.2%). No layout at this budget passed the gates for all four slope profiles. All three 3,072-response layouts passed across the four core slope profiles for word rotation. With participant assignment, every core comparison failed the 80% power gate. Under balanced slopes, 32×48 to 96×16 reduced mean width from 23.23 to 16.63 pp and raised power from 29.7% to 48.0%.

Switching from 32×48 to 96×16 with word rotation increased power by 4.35 points (approximate paired 95% MC CI 2.24 to 6.45) under participant-dominant slopes and decreased it by 9.23 points (−11.51 to −6.95) under word-dominant slopes. At 3,072 responses, power for the eight-point effect ranged from 94.8% to 99.4% across the four profiles. For the five-point effect under balanced slopes, the highest power was 71.7% at 64×48 (95% MC CI 69.7–73.6%). In the word-assigned sensitivity check, moving from 32×48 to 96×16 raised mean interval width from 16.71 to 23.54 pp and lowered power from 47.4% to 31.2%. At 96×16, null detection was 7.0% (95% MC CI 6.0–8.2%), which left the gate indeterminate.

#### Correlations and omission of slopes

With slope SD 0.25, changing one true intercept–slope correlation from zero to ±0.5 changed eight-point power by at most 1.13 percentage points under the primary fit. Those generating changes increased mean interval width by 0.01–0.11 pp. With nonzero true correlation, fixing the fitted correlation at zero on the same responses shortened mean width by 0.05–0.35 pp and raised estimated power by at most 1.71 points. The comparisons cover two layouts with 1,536 responses.

In the balanced, zero-correlation DGP under word rotation, the intercept-only fit gave null detection rates of 8.9% at 32×48 (95% MC CI 7.74–10.25%) and 8.7% at 96×16 (7.55–10.05%). The lower limits exceeded 7.5% in both cases. The slope-and-correlation fit gave 3.3% and 3.6%. Under the planned effect, coverage at 32×48 fell from 96.1% to 90.6% when both slopes were omitted. The higher detection rate under the intercept-only fit came with narrower, less well calibrated intervals.

#### Continuous covariates

All 28,672 fits across 12,288 datasets were valid. In the additive covariate DGP, adjustment reduced mean full width from 23.27 to 19.33 pp at 32×48 and from 16.64 to 14.82 pp at 96×16 under participant assignment. Detection of the eight-point effect increased by 7.81 pp (approximate paired Monte Carlo 95% CI 5.26 to 10.37) and 9.86 pp (7.36 to 12.37), respectively. Under word rotation, mean width fell by less than 0.03 pp at both layouts and the paired detection changes were within one point of zero.

With a vocabulary interaction of +0.35 logits per SD, the fitted interaction averaged 0.347–0.355 across the two layouts and both marginal effect sizes. Its interval coverage ranged from 979/1,024 to 982/1,024. Modeling the interaction also shortened the marginal effect interval from 11.98 to 11.27 pp at 32×48 and from 11.90 to 11.68 pp at 96×16 under the planned effect. These comparisons use the [covariate scenarios](l2_vocabulary_covariate_scenarios.csv), [summary](l2_vocabulary_covariate_summary.csv), and [replication-level estimates](l2_vocabulary_covariate_fits.csv.gz).

### 4D. Partial presentation and participant–word compatibility

All layouts schedule 1,536 binary recall responses from 96 participants, with 16 words per participant and one response per participant–word pair. Full presentation uses the same 16 words. Two partial layouts use a bank of 48 words in six eight-word blocks. The linked layout uses AB, BC, CD, DE, EF, and FA. The disjoint layout uses AB, CD, and EF. Each list has two versions with reversed learning methods. Each participant gets eight words per method. In full presentation, each word is scheduled for 48 participants per method. In the partial layouts, each word is scheduled for 16 per method. The scheduled person–word graph has one component in each of the full and linked layouts and three in the disjoint layout.

The target is the marginal retrieval-minus-restudy recall difference in an equally weighted population of compatible and incompatible pairs. A pair has compatibility score `m=+1` when its participant and word type match, and `m=-1` otherwise. Both participant types appear equally on each list, and each list contains eight words of each type. Treatment is coded `x=-0.5/+0.5`. The logit includes independent normal participant/word intercepts (SD 0.8/0.6) and slopes (SD 0.25/0.25 in the reference). Each DGP calibrates both marginal rates by normal integration to 42% vs. 42% or 42% vs. 50%. The moderated marginal null permits nonzero conditional differences within compatibility strata.

| Generating profile | Change from the balanced-slope reference |
|---|---|
| Baseline | No compatibility effect, equal versions, complete observation |
| Learner- or word-dominant | Slope SDs 0.40/0.15 or 0.15/0.40, respectively |
| Match | Adds `0.65*m` to the logit, with versions equally assigned within participant type |
| Moderation | Adds `0.65*m + 0.60*x*m` to the logit |
| Allocation bias | Uses the Match outcome model. Compatible pairs are 75% in retrieval and 25% in restudy. Overall participant×method and word×method counts remain balanced. All four method×compatibility cells are exposed. |
| MCAR | Uses the Match outcome model and retains each response with probability 0.85 |
| Compatibility-related missingness (MAR) | Uses the Match outcome model and retains with probability `0.85 + 0.12*m*(2*x)`. Compatibility is known for scheduled rows. |
| Outcome-related missingness (MNAR) | Uses the Match outcome model and retains with probability `0.85 + 0.12*(y-p_x)*(2*x)`, with `p_x` the true marginal outcome probability for the method. |

The allocation run uses 54 DGPs × 3,072 datasets (seed 86614001). All nine generating profiles are crossed with the three layouts and both effect sizes. Independent L'Ecuyer-CMRG streams use shared participant and word draws, outcome uniforms, and retention uniforms at matching IDs. Analyze all cases with native `bayes_laplace`, correlated intercepts/slopes as supported by the schedule, and joint-delta 95% intervals. For Match, Moderation, Allocation bias, MAR, and MNAR, add a fit to the identical observed responses with `m` and `x*m` as fixed effects. Standardize this fit to a fixed 50:50 mix of the two compatibility strata. The test estimates a population contrast and uses zero bootstrap refits. A research-only interface accepts a precomputed schedule and numeric design matrices through internal functions. [Reproduction script](l2_vocabulary_sparse.R)

Coverage is conditional on valid fits. Detection uses all planned datasets. Exact pointwise 95% Monte Carlo bounds and the gates in §4B apply. Report empirical bias, its MCSE, width, exposure balance, and observed response counts. For biased allocation and selective observation, coverage and error refer to the same equal-stratum target as in the randomized complete conditions. Those scenarios test deviations from the DGP assumptions used in routine planning.

All 258,048 native fits across 165,888 datasets were valid. The run reused the isolated macOS R 4.6.1 build of estimator source `bcd5ee0` and took 14 minutes 40 seconds with four workers, including output serialization. Independent recomputation from the saved fits matched all coverage, detection, bias, width, and Monte Carlo gate values. A separate integration of the outcome and retention probabilities matched the mean raw contrasts to within 0.064 percentage points across the 54 scenarios. [Scenarios](l2_vocabulary_sparse_scenarios.csv) · [Summary](l2_vocabulary_sparse_summary.csv) · [Replication-level results](l2_vocabulary_sparse_fits.csv.gz) · [Bias figure](l2_vocabulary_sparse_bias.png) · [Plot script](plot_l2_vocabulary_sparse.R)

#### Full versus partial presentation

The table gives mean full interval width in percentage points and detection of the eight-point effect when pair compatibility has no effect on recall.

| Slope SDs, participant/word | Full, 16 words | Linked, 48 words | Disjoint, 48 words |
|---|---:|---:|---:|
| Balanced, 0.25/0.25 | 11.68 / 79.2% | 10.55 / 86.7% | 10.53 / 86.3% |
| Participant-dominant, 0.40/0.15 | 11.40 / 82.0% | 10.58 / 86.5% | 10.56 / 86.8% |
| Word-dominant, 0.15/0.40 | 12.67 / 71.4% | 10.81 / 83.9% | 10.80 / 83.7% |

For word-dominant slopes, linked partial presentation reduced width by 1.86 pp (paired approximate 95% Monte Carlo CI 1.80–1.92) and increased power by 12.53 points (10.58–14.48) relative to full presentation. Linked versus disjoint width differed by less than 0.02 pp in all three profiles. The corresponding paired power differences were smaller than 0.4 points in magnitude, and their 95% intervals included zero. Both partial designs cross the learning methods within every form and assume common random-effect distributions across forms. Form-specific difficulty is fixed at zero, and the reported intervals concern the population learning-method contrast.

Both partial layouts passed the coverage, null detection, and 80% power gates for all three slope profiles in this comparison. Full presentation passed all three gates in the participant-dominant case, was indeterminate for 80% power in the balanced case, and fell below the power gate in the word-dominant case.

#### Compatibility and effect modification

The balanced assignment averages over both directions of the compatibility advantage. Because the two method versions switch whole blocks, compatibility can still account for part of the difference between methods within each participant. The linked layout with a compatibility main effect gave planned intervals of width 13.64 pp and power 70.0% without adjustment. Modeling compatibility and standardizing to the same target shortened the interval to 10.22 pp and increased power to 88.4%. Coverage was 99.4% and 96.2%, respectively.

In the moderation DGP, the conditional recall differences after integrating the participant and word effects were −11.21 and +11.21 pp for incompatible and compatible pairs under the equal-stratum null. Under the eight-point marginal effect, they were −3.71 and +19.71 pp. Testing the average contrast and estimating how compatibility changes the learning effect are different questions.

#### Imbalanced allocation and missingness

At a true mean effect of zero, the 75:25 allocation imbalance gave unadjusted biases of 12.87–12.88 pp across the three layouts. For the linked layout, unadjusted null detection was 3,058/3,072 (99.54%, exact 95% MC CI 99.24–99.75%). Compatibility adjustment and standardization gave mean bias −0.03 pp and null detection 125/3,072 (4.07%, 3.40–4.83%). The margins by participant and by word were exactly balanced in every dataset.

MCAR retained about 1,306 of 1,536 scheduled responses and did not introduce material bias. With the linked layout, planned power under the unadjusted fit went from 70.0% with complete observation to 66.2% under MCAR. With compatibility-related missingness the same average response count left 42.9% compatible pairs in observed restudy and 57.1% in observed retrieval. This shifted the unadjusted estimate by +2.17 pp at the null, with 96.7% coverage and 3.3% null detection. Compatibility adjustment reduced bias to +0.01 pp. The bias therefore adds information beyond coverage or the null detection gate.

Outcome-related missingness left unadjusted null bias of 6.07–6.19 pp across layouts. Compatibility adjustment left 5.75–5.86 pp. In the linked layout, null detection rose from 38.64% (95% MC CI 36.91–40.39%) to 56.64% (54.87–58.40%) after compatibility adjustment. Adjustment narrowed the intervals around an estimate that retained substantial bias. The two missingness mechanisms demonstrate the importance of whether the variable driving selection is available at analysis time.

### 4E. Prior word knowledge and vocabulary screening

The primary estimand is the marginal retrieval-minus-restudy difference among genuinely unknown participant–word pairs in a common reference population. The secondary estimand averages over all candidate pairs, including retention of previously known words. The linked starting lists from §4D give 96 participants 16 initial words each from a 48-word bank. All policies use the same screening accuracy and the same per-pair outcome model.

| Policy | Screening and teaching |
|---|---|
| Fixed all | Screen the initial 16 words for measurement and teach all 16 with the original two counterbalanced versions. |
| Fixed screened | Screen those 16, teach only words classified as unknown, and keep the original method assignments. |
| Top-up screened | Start with those 16 and screen the remaining bank words in a random order until 16 appear unknown or the 48-word bank is exhausted. Randomize methods across the selected words within each participant, with counts differing by at most one. |

Let standardized vocabulary size `V` and word familiarity `F` be independent standard normal variables. The logit of prior knowledge is `alpha + gamma*V + 0.70*F`. Calibrate `alpha` and `gamma` so that true known rates at vocabulary size −1 and +1 SD, averaging over words, are 10%/20%, 10%/50%, or 30%/70%. Generate the screen separately conditional on true knowledge. Sensitivity/specificity are 1/1 or 0.90/0.95 in all profiles. The 30%/70% profile and the effect-modification profile also use 0.75/0.85.

The binary recall logit uses participant/word intercept SDs 0.8/0.6, slope SDs 0.25/0.25, and independent centered normal random effects. For unknown pairs, the fixed logit includes `0.35*F`. Two profiles with 10%/50% known rates add `0.45*V`, and one of those also adds `0.35*x*V` for centered unit-step method coding `x=-0.5/+0.5`. Recalibrate the two unknown-pair averages to 42% versus 42% or 42% versus 50% at every profile. For known pairs, both marginal recall probabilities are 90%. The modular outcome DGP assumes the same per-pair learning opportunity, posttest delay, and screening response procedure in every policy. The whole-bank effect is the unknown-pair effect times the population share of unknown pairs.

Use the native regularized Laplace estimator with crossed intercepts/slopes and a joint-delta 95% interval. The fixed-all outcome analysis uses method as its fixed predictor and evaluates the all-pair contrast. For unknown-word analyses, use method, `V`, `F`, and method×`V`, and standardize to the same generating distribution of `V,F` among truly unknown pairs. Thus operational screening is noisy, while the reference distribution and the two continuous covariates are available by construction. The noisy-screen scenarios also fit the genuinely unknown pairs on the starting lists and among the top-up words as oracle comparisons. Within each data-generating scenario the policies share the people, words, knowledge states, screening draws, and posttest uniforms. Compare coverage, bias, width, null detection, and power on each corresponding scale. The fit count omits duplicates of the oracle analyses under perfect classification.

Run 24 DGPs × 2,048 independent replications with master seed 87761001 and zero bootstrap refits. Use the exact pointwise Monte Carlo intervals and gates in §4B. At positive effects, a higher detection rate can accompany bias against the true unknown-pair target. Record screening and teaching counts, counts of truly unknown words by arm, coverage of words in both arms among the genuinely unknown responses, and shortages for the top quartile of each sample's vocabulary size distribution. The re-execution script reuses fixed-all and fixed-oracle records across screening accuracies when the generating profile and effect are identical. The computational comparison uses precomputed lists, continuous matrices, and the existing native estimator through internal functions. [Reproduction script](l2_vocabulary_prior_knowledge.R)

The run produced 204,800 valid fit records across 49,152 datasets. After preserving the first 12 checkpoints, the remaining cases reused the identical fixed controls, for 155,648 native fits in the completed run. It took about 15 minutes with four workers, including the checkpoint restart and output writing. The engine was the isolated R 4.6.1 build of source `bcd5ee0`. Independent aggregation verified the replication keys, fixed-control equality, screening counts, interval metrics, and decision gates. A separate normal integral confirmed the calibrated population known fractions and the all-pair effects. [Scenario truths](l2_vocabulary_prior_knowledge_scenarios.csv) · [Fit summary](l2_vocabulary_prior_knowledge_summary.csv) · [Exposure summary](l2_vocabulary_prior_knowledge_exposure.csv) · [Replication-level fits](l2_vocabulary_prior_knowledge_fits.csv.gz) · [Replication-level exposures](l2_vocabulary_prior_knowledge_exposures.csv.gz) · [Screening-effort figure](l2_vocabulary_prior_knowledge_screening.png) · [Plot script](plot_l2_vocabulary_prior_knowledge.R)

#### Screening effort and precision with perfect classification

The table shows the number of genuinely unknown words taught per participant, within each replication's top vocabulary quartile. The two learning policies target the same 8 pp unknown-pair effect. The fixed list always screens 16 words. Top-up effort is the total number screened per upper-quartile participant.

| Known rates at −1/+1 vocabulary SD | Fixed unknown words | Top-up unknown words | Top-up screening | Top-up cannot reach 16 | Fixed / top-up power for the unknown effect |
|---|---:|---:|---:|---:|---:|
| 10%/20% | 12.50 | 16.00 | 20.54 | <0.01% | 82.0% / 87.1% |
| 10%/50% | 6.98 | 14.98 | 36.22 | 21.1% | 74.1% / 84.5% |
| 30%/70% | 4.20 | 11.98 | 46.12 | 72.6% | 56.7% / 82.1% |

For the 30%/70% profile, 20.6% of upper-quartile participants had zero truly unknown words in at least one method on the fixed lists. Top-up reduced this proportion to 0.35%. Power increased by 25.39 points, with approximate paired 95% Monte Carlo limits of 22.75–28.04. The power targets the population contrast standardized among unknown pairs, with all 96 participants sampled across vocabulary levels.

The all-word contrast in the 30%/70% profile is 4 pp because the population known fraction is 50% and retention has a common truth in both methods. The fixed-all fit averaged 4.02 pp. In the 10%/50% profile, the corresponding truth is 5.66 pp, and the fixed-all fit averaged 5.72 pp. These are separate conclusions about the posttest score across the candidate bank.

#### Imperfect screening

With 90% sensitivity and 95% specificity, the 10%/50% profile gave unknown-effect detection of 69.9% with fixed exclusions and 81.3% with top-up. The exact 95% Monte Carlo interval for top-up was 79.6–83.0%, leaving the 80% gate indeterminate. In the 30%/70% profile, upper-quartile top-up needed an average of 43.82 screens and taught 13.99 apparently unknown words, of which only 10.73 were genuinely unknown. Shortages below 16 fell to 49.5%, while the number of truly unknown words taught also fell relative to perfect screening. The estimated eight-point acquisition effect averaged 7.42 pp and detection fell to 73.9% (95% MC interval 71.9–75.8%). The oracle fit within the same top-up selection estimated 8.17 pp and detected 78.3%.

At 75% sensitivity and 85% specificity in the same high-knowledge profile, upper-quartile participants learned an apparent 15.60 words on average after screening just 38.57, yet only 8.41 of these words were truly unknown. The corresponding estimated acquisition effect averaged 6.48 pp, with 61.3% detection and 93.9% coverage. The oracle fit on genuinely unknown top-up words averaged 8.14 pp, with 73.4% detection. Across the benchmark, nominal-null detection ranged from 2.7% to 5.0%. Good null calibration can coexist with attenuation when the true effect is positive.

With effect modification, the calibrated eight-point acquisition contrast averages over truly unknown pairs. Its true value at vocabulary size −1 SD is +2.78 pp and at +1 SD is +16.68 pp. The same model family supports a marginal null with opposing subgroup effects (−5.12 and +8.62 pp). This separates effect modification in acquisition from the changing contribution of retention to posttest scores.

### 4F. Vocabulary recall at the floor and ceiling

The target is the average retrieval-practice minus restudy difference in delayed binary meaning recall for new learners and words. All learner–word pairs are initially unknown. Give each learner the same number of opportunities per word and hold the feedback and test delay fixed. A layout `N_J_m` uses `N` learners, a `J`-word bank, and `m` taught words per learner. Cross 48, 96, and 144 learners with 16-, 32-, and 48-word banks at `m = 16`. Add `96_48_8`, `96_48_24`, `32_48_48`, `48_32_32`, and `64_24_24`. The 14 layouts let us vary the number of learners, distinct words, and observations per word. Partial presentation uses overlapping lists of four-word blocks, with two equally frequent method versions per list. All schedules have one connected component and exactly `m/2` responses per learner per method and `N*m/(2*J)` per word per method. Each observed learner–word pair contributes one response.

**Global floor and ceiling.** Set true population restudy accuracy to 10%, 42%, or 90%. Use an equal-marginal null and a five-percentage-point retrieval advantage at each baseline. A separate positive-effect comparison holds the conditional method coefficient at the value that gives 42% vs. 47% under the reference DGP. Reuse its 42% cell. This comparison lets the percentage-point effect shrink at the floor and ceiling. The generator is standard logit with no guessing asymptote.

**Concentrated floor and ceiling.** Hold restudy at 42% and the true marginal contrast at zero or five points. Use six representative layouts. In the main mixtures, independently assign each member of one group (learners or words) a 10% chance of a centered +3 or −3 logit shift. Generate the extreme response totals by Bernoulli draws. At three layouts, the normal controls add 0.81 to the intercept variance of that group, matching the mixture's full intercept variance. At the two 96-learner layouts with `m = 16`, use an additional ±1.5 logit check for high learners and low words. Group membership, shift, and variances are planning assumptions.

**Random slopes.** Code method as −0.5/+0.5 so each slope is a retrieval-minus-restudy log-odds difference. Baseline learner/word random intercept SDs are 0.8/0.6 and slope SDs are 0.25/0.25, with both intercept–slope correlations zero. Switch slope SDs to 0.40/0.15 or 0.15/0.40 at global floor and ceiling and at representative concentrated extremes. At the two 96-learner layouts with `m = 16`, also use 0/0 slopes and vary one correlation at a time to −0.5 or +0.5 at the global floor and ceiling. A correlation changes the marginal arm probabilities even if the conditional coefficient is zero. Recalibrate both arms to their target means after changing the covariance, including for the equal-marginal null. Fits retain both random slopes and their correlations. A paired sensitivity fit omits the slopes in selected cells with both balanced and zero generating slopes.

**Learner and word covariates.** In three layouts at 42%, draw independent standard normal working-memory and concreteness scores with logit coefficients 0.45 and 0.35. Center the two word classes at −1/+1 with coefficient 0.20. Reduce the residual intercept variances by the corresponding continuous main-effect variances. In the audio profiles, the learner discrimination group `D` is −1/+1 and the word feature `F` indicates a confusable contrast (0/1). Use fixed terms `0.15*D - 0.30*F + gamma*D*F`, with `gamma` at 0.35 or 0.70. This interaction changes baseline recall for specific learner–word pairs. Method interactions are separate generating terms: working memory × method at 0.35 logits per SD, word class × method at 0.35 logits per unit (0.70 logits from one class to the other), or a joint setting with 0.70 and −0.35, respectively, and stronger phonological matching. Two joint checks combine a measured predictor profile with an unobserved concentrated tail. Coefficients and independent covariate distributions are sensitivity assumptions. For each simulated dataset, compare method-only, adjustment for measured predictors (with `D*F` when present), and addition of the relevant method interactions. Standardize each adjusted fit to the same generating distribution of learner and word covariates. Generate equal-marginal nulls even when treatment effects vary by learner or word class.

**Validation and reporting.** Before the pilot, fix the 304 generating scenarios, the fits, and the 48- and 2,048-dataset stages. Use independent master seeds `88741001` and `88742001` for the two stages, and pair the underlying draws by replication across scenarios within a stage. Use an isolated installation of `glmmplanr` with the native `bayes_laplace` estimator and joint-delta 95% intervals (zero bootstrap refits). Verify the arm truths by ordinary R integration and separate Gaussian quadrature. Use the exact pointwise Monte Carlo intervals and gates in §4B. Coverage is among valid fits and detection is among all scheduled fits. Mean interval width, its 80th percentile, empirical SD, and bias are in percentage points. In the method-only model, an arm with exclusively successes or failures implies separation for the unpenalized fixed-effects design. Also count empty and constant-response strata formed by the categorical predictors and method. This latter count measures data support in the strata, and separation for a model with continuous predictors depends on its full design matrix. Record all-wrong/all-correct learner and word counts alongside the unit exposures (`m` and `N*m/J`). [Reproduction script](l2_vocabulary_extremes.R) · [Independent summary and paired comparisons](summarize_l2_vocabulary_extremes.R)

The isolated pilot on macOS and R 4.6.1 produced 14,592 datasets and all 18,528 fits were valid. Independent aggregation reproduced the fit counts and interval decisions. The pilot was used to check generation, schedules, and numerics. The large run uses the same settings and a new seed. [Pilot scenarios](l2_vocabulary_extremes_pilot_scenarios.csv) · [Pilot summary](l2_vocabulary_extremes_pilot_summary.csv) · [Pilot fit-level records](l2_vocabulary_extremes_pilot_fits.csv.gz) · [Pilot paired contrasts](l2_vocabulary_extremes_pilot_paired.csv)

The isolated large run used the estimator code at `bcd5ee0` and took 22 minutes 28 seconds with eight workers. It produced 622,592 datasets and 790,528 valid fits. Recomputing the summaries from the saved trials reproduced the counts and interval decisions. Across 386 scenario–model combinations, the observed coverage ranged from 92.77% to 98.39%, and each exact pointwise lower bound exceeded 90%. The lower bound was at least 80% for 4 of 14 layouts at 10% versus 15% and for 10 of 14 at 90% versus 95%. The highest detection at 42% versus 47% under the reference slopes was 1,313/2,048 (64.1%, 95% MC interval 62.0–66.2%), for `144_48_16`. Every equal-marginal null passed the 7.5% upper-bound gate using the models with random slopes. The intercept-only fit was indeterminate in one null configuration. [Generating scenarios](l2_vocabulary_extremes_large_scenarios.csv) · [Fit summary](l2_vocabulary_extremes_large_summary.csv) · [Paired contrasts](l2_vocabulary_extremes_large_paired.csv) · [Replication-level fits](l2_vocabulary_extremes_large_fits.csv.gz) · [Figure with 95% MC intervals](l2_vocabulary_extremes.png) · [Plot script](plot_l2_vocabulary_extremes.R)

Each method arm contained both correct and incorrect answers in every generated dataset. Among the 90,112 covariate-profile datasets, 1,106 had at least one empty categorical stratum, and 160 had at least one observed stratum with only one response category. These stratum checks describe the data support for adjustment.

The default isolated command runs the pilot. To reproduce both stages using the same isolated installation, use a directory for each stage and pass its path, the stage (`pilot` or `large`), and the worker count to `Rscript --vanilla validation/l2_vocabulary_extremes.R`. Set `R_LIBS` to the isolated run's `library` directory before calling `Rscript`. The script saves completed scenarios into its output directory for resuming with the same code and settings. Run `Rscript --vanilla validation/summarize_l2_vocabulary_extremes.R OUT STAGE` to reproduce the paired contrasts and audit the fit summaries.

#### Effect scale at the global floor and ceiling

For 96 learners with 16 words each from a 48-word bank, the conditional log odds ratio required for a five-point effect was 0.518 at a 10% restudy rate, 0.246 at 42%, and 0.809 at 90%. Holding the coefficient at 0.246 gave true effects of 2.17, 5.00, and 1.84 pp. In the five-point comparison the intervals on the probability scale were narrower toward the endpoints. With the log odds ratio held fixed, the absolute effect also shrank and detection was lower toward both endpoints.

| True restudy rate | Detection at fixed +5 pp (95% MC CI) | Detection with fixed log odds ratio (95% MC CI) | Mean full width at fixed +5 pp |
|---|---:|---:|---:|
| 10% | 1,613/2,048 (78.8% [76.9%, 80.5%]) | 403/2,048 (19.7% [18.0%, 21.5%]) | 7.18 pp |
| 42% | 995/2,048 (48.6% [46.4%, 50.8%]) | 995/2,048 (48.6% [46.4%, 50.8%]) | 10.51 pp |
| 90% | 1,919/2,048 (93.7% [92.6%, 94.7%]) | 363/2,048 (17.7% [16.1%, 19.4%]) | 5.84 pp |

In this layout at 90% versus 95%, a mean of 36.04 out of 96 learners answered all 16 words correctly, and 7.07 out of 48 words were correct on all 32 presentations. The planned 95% interval covered the true effect in 1,976/2,048 datasets (96.5%). Under the equal-marginal 90% null, the error rate was 54/2,048 (2.6%, 95% MC interval 2.0–3.4%). At 42% versus 47%, the four independently simulated complete-presentation layouts with 1,536 responses agreed with §4C within 0.04 pp for mean interval width, and each detection-rate difference had an approximate independent 95% MC interval containing zero.

#### Concentrated extremes and variance-matched controls

With the same 96-learner sparse layout and a 42% versus 47% marginal contrast, high learners under the 10% mixture produced on average 2.58 all-correct learners versus 0.92 under the learner Gaussian control. Low words produced 0.88 all-wrong words versus 0.17 under the word Gaussian control. The differences in the mean full interval width between the mixture and its variance-matched control were +0.05 to +0.13 pp across the four concentrated tails in this layout. Thus the number of extremes can change visibly while 95% interval width for the average treatment effect changes modestly in the same design.

For the low-word mixture, hold 96 learners and 16 words per learner fixed. Increasing the bank from 16 to 48 reduced mean interval width from 10.84 to 10.13 pp with learner-dominant slopes and from 11.99 to 10.36 pp with word-dominant slopes. The corresponding power changes were 41.5% to 48.7% and 36.2% to 48.3%. Which unit varies more in its response to treatment affects the gain from increasing the number of distinct words, even at the same test average and exposure per learner.

#### Sensitivity to omitting the random slopes

In the 96-by-16 complete design at a 10% marginal null, retaining both random slopes gave a 2.34% false positive rate and 97.66% coverage when both generating slope SDs were 0.25. Omitting them in the same responses increased the false positive rate to 7.23%, reduced coverage to 92.77%, and narrowed the mean interval from 7.07 to 5.75 pp. When both generating slope SDs were zero, the intercept-only fit gave 4.30% null detection and 95.70% coverage. The maximal fit gave 1.61% and 98.39% in that zero-slope DGP, with mean width 6.89 pp compared with 5.74 pp for the intercept-only fit.

Changing one intercept–slope correlation at a time to ±0.5 in the two 96-learner layouts gave 1.9–3.6% null detection and 95.7–98.1% coverage across the equal-marginal null and five-point cases. The equal-marginal null required conditional method coefficients between −0.067 and +0.067 logits when the arm-specific random-effect variances differed.

#### Measured predictors and treatment effect modification

For 96 learners with 16 words from a 48-word bank at 42% versus 47%, adjusting for the prognostic learner and word predictors shifted mean width from 10.49 to 10.45 pp. With the stronger phonological matching term the change was from 10.34 to 10.16 pp. The largest change in this layout was under the joint moderation profile. The method-only, predictor-adjusted, and effect-modification fits gave mean widths of 11.45, 11.18, and 10.17 pp and detection rates of 40.8%, 41.0%, and 47.8%. Adding moderation to the predictor-adjusted fit narrowed the interval by 1.01 pp (paired 95% MC interval 0.99–1.03 pp) and increased detection by 6.79 points (4.93–8.65). The effect-modification fit recovered mean interaction coefficients of 0.706 and −0.347 against generating values of 0.70 and −0.35. All three estimated the same five-point marginal target, with absolute mean bias below 0.08 pp in this balanced layout.

### 4G. Common-item linking for Rasch bank expansion

This benchmark plans a calibration study to place 12 new Rasch items on an existing bank's scale. The operational target is recovery of the new item difficulties and the reference-scale origin. The design compares the number and location of anchors at a fixed examinee workload. `equateIRT` supports test score equating after the forms have been placed on a common scale. [Battauz (2015)](https://www.jstatsoft.org/article/view/v068i07)

The reference population has normal ability with mean 0 and SD 1. Its 24 bank difficulties consist of 12 evenly spaced values between −2.2 and 2.2 logits and 12 between 1.3 and 3.5 logits. The new bank has 10 items evenly spaced between −2 and 2 plus two at −3 and 3. The new population has normal ability with SD 1 and mean either 0 or −0.75. Cross 120, 240, 480, 960, and 1,920 new examinees with 4, 8, or 12 anchors spread evenly across the broad or hard reference difficulty bands. Each examinee answers all assigned anchors and 16 minus the number of anchors in new items. A cyclic assignment gives every new item exactly the same exposure. There is no dropout, guessing floor, or local dependence. Additional cells hold new-item exposure at eight per examinee, or drift one quarter of the anchors by +0.60 logits, one half by +0.80, or all by +0.40. These changes affect the new cohort only and are uniform shifts of item difficulty.

Fit the existing C++ Rasch fixed-item MML model in each population with its origin at the population's estimated center. On the new data, obtain the additive Rasch linking shift as the mean reference-minus-new difference among finite anchors. Require at least two usable anchors. Do this both with true reference difficulties and with bank estimates from a separate calibration of 500 examinees. The same new-data fit serves both comparisons. Record the centered estimates and shifts, allowing re-use of the retrieved bank in a downstream assessment simulation. Constant-response new items have no finite fixed-item MLE and remain in the study denominator. The linked performance target per study is absolute scale-origin error ≤0.15 logits, new-item RMSE ≤0.35 logits, and all 12 absolute difficulty errors ≤0.60 logits. A candidate passes when a two-sided pointwise exact 95% Monte Carlo lower bound for the joint attainment rate is at least 80%, fails when the upper bound is below 80%, and is otherwise indeterminate. The tolerances are benchmark inputs. Also report the maximum absolute difference in the 12-item expected total score across ability −2 to 2 in 0.25-logit increments. This outcome is a test characteristic curve diagnostic.

The pilot uses 32 independent replications per scenario (master seed `90541001`). The large run uses 2,048 per scenario (seed `90542001`). Scenarios within a stage share calibration draws where applicable. In six selected scenarios, compare native and TAM fixed-item MML on the same first 64 large-run response datasets. Select both 240 and 960 examinees, each at 4 broad anchors, 8 broad anchors, and 12 hard anchors, with no drift and true new-population mean −0.75. Compare the same Rasch mean-anchor shift, full linked-bank differences, and fit CPU time under both engines. This is a paired comparison of two calibrators followed by the same transformation. [TAM estimation reference](https://search.r-project.org/CRAN/refmans/TAM/help/tam.mml.html)

#### Results

The isolated R 4.6.1 run used the previously checked `glmmplanr` 0.1.0 installation. The pilot completed 3,040 datasets. The large run completed 194,560 datasets in 93 seconds using eight workers. The native fit returned finite item estimates in 194,559 datasets, and 192,389 had a complete linked new bank under both anchor references. The other 2,171 datasets all contained at least one constant-response new item. In one such dataset (case 56, replication 265), the native EM solver also reached its iteration limit. The denominator retains these datasets. The independent audit recalculated item errors and attainment from the saved centered bank estimates and shifts, and verified the 95% Monte Carlo intervals. [Pilot summary](rasch_anchor_pilot_summary.csv) · [Large scenario table](rasch_anchor_large_scenarios.csv) · [Large summary](rasch_anchor_large_summary.csv) · [Paired comparisons](rasch_anchor_large_paired.csv) · [Trial records and centered bank estimates](rasch_anchor_large_fits.csv.gz) · [True new bank](rasch_anchor_large_new_bank.csv)

A cohort whose mean ability was 0.75 below the reference yielded a mean bias of +0.757 logits before linking in the 480-examinee, eight-broad-anchor cell. With stable anchors, mean bias in the origin was −0.003 logits with known reference difficulties and −0.004 when they came from 500 reference examinees.

| Anchors per examinee at new N = 480 | Observations per new item | Joint attainment with known anchors | Joint attainment with 500-person reference estimates |
|---:|---:|---:|---:|
| 4 broad | 480 | 90.1% | 78.0% |
| 8 broad | 320 | 89.3% | 82.5% |
| 12 broad | 160 | 67.7% | 62.8% |
| 8 hard | 320 | 77.7% | 68.8% |

The eight broad anchors achieved 1,689/2,048 (82.5%, 95% MC interval 80.8–84.1%) with the previously calibrated reference. The paired gain over four anchors at equal workload was 4.44 points (approximate 95% paired MC interval 2.42–6.47). With twelve anchors, the gain in origin precision was outweighed by the reduction in new-item exposure. Holding new-item exposure at 320 instead gave 89.3% for eight broad and 90.3% for twelve broad anchors with true difficulties. The latter assigned 20 items rather than 16 per person. With the 500-person reference, the minimum tested sample sizes under broad, stable anchors were 960, 480, and 960 for four, eight, and twelve anchors, respectively. For eight hard anchors, the minimum was 960. [Joint assurance figure](rasch_anchor_linking.png) · [Plot script](plot_rasch_anchor_linking.R)

At `N = 1920` with eight broad anchors and known reference difficulties, stable anchors met all criteria in 2,047/2,048 studies. When a quarter shifted by +0.60 logits, attainment was 987/2,048 and mean origin error was −0.151 logits. When half shifted by +0.80 or all shifted by +0.40, the mean error was about −0.40 and neither scenario met the criteria in any replication. A common anchor difficulty shift can be offset exactly by translating the new cohort ability and new items together. Separation of these changes requires information about invariance beyond this response design.

On the 384 matched datasets, native and TAM each returned 379 complete new banks. The maximum corresponding linked-item difference was 0.00074 logits and both engines gave the same attainment decision for every study. Summed process CPU time for the fit calls was 1.03 s native and 5.76 s TAM, including matrix assembly for TAM. This comparison covers the specified Rasch fits on the same responses. [TAM paired results](rasch_anchor_large_tam.csv)

#### Connection to downstream assessment

An offline `mirtCAT` 1.14 smoke check passed one linked bank (24 existing and 12 new items with `a1 = 1`, `d = -difficulty`) to `generate.mirt_object()`. It used separate simulated response patterns from the true bank for 24 examinees and administered eight items per examinee with maximum information selection and EAP scoring. The saved linked estimates can be reconstructed by adding `oracle_shift` or `estimated_shift` to the twelve `new_centered_*` columns. This is a connectivity and parameterization check. A CAT design study would propagate errors from multiple calibrated banks to ability estimates, stopping decisions, exposure and classification. A shortened fixed form can be evaluated on the same held-out population with its content requirements carried through item selection. [Offline CAT check](rasch_anchor_cat_smoke.R) · [mirtCAT methods](https://www.jstatsoft.org/article/view/v071i05) · [TestDesign fixed-form assembly](https://choi-phd.github.io/TestDesign/)

To reproduce, use `sh validation/run-isolated.sh rasch_anchor_linking` for the pilot. For the large run, use the same isolated installation, set `R_LIBS` to its `library` directory and run `Rscript --vanilla validation/rasch_anchor_linking.R OUT large WORKERS` with a separate existing output directory and TAM installed. Blocks resume from the same code and settings. For each completed stage, run `Rscript --vanilla validation/summarize_rasch_anchor_linking.R OUT STAGE` to verify the summaries and compute paired contrasts. [Simulation script](rasch_anchor_linking.R) · [Independent audit](summarize_rasch_anchor_linking.R)

### 4H. Fixed and random passage effects with two to six shared passages

The benchmark asks how to interpret planning output with few passages when every passage is used in both conditions. Assign 48 readers equally to A/B. Everyone answers 24 binary questions. Use 2, 3, 4, 5, or 6 passages, with 12, 8, 6, 5/5/5/5/4, or 4 questions per passage. Randomly allocate the shorter passage to a position in the five-passage bank within each replication. Weight the passages equally for the fixed-bank contrast. A sensitivity check doubles the number of readers for two and five passages in the heterogeneous, locally dependent profile. The one-passage-per-condition boundary is reported from §4B, with a rank check for the design including fixed passage and condition effects.

In the binary logit DGP, independent reader, question, and passage intercept SDs are 0.7, 0.4, and 0.5. Add a reader-by-passage random intercept (SD 0 or 0.4) and a passage-specific condition slope (SD 0 or 0.4). Code A/B as −0.5/+0.5 with equal variance under both conditions. Adjust the two population logits to give marginal probabilities 60% versus 60% or 60% versus 70% after integrating over all random effects. Resample passages in each replication. The **sampled-passage contrast** conditions on the realized passage intercepts and slopes while averaging over new readers, questions within those passages, and reader-by-passage deviations. The **population-passage contrast** also integrates over fresh passages. Under the population null, the sampled-passage contrast can be nonzero when the generating slopes vary. Report detection of zero on those datasets as such. Under the common-effect null both contrasts are exactly zero.

Use four paired native regularized Laplace fits per response dataset, with random reader/question intercepts. (1) Fixed passage intercepts with a common condition effect. (2) Fixed passage intercepts and condition-by-passage effects. Both contrasts standardize across the sampled passages. (3) Random passage intercepts with a common condition effect, as in the built-in intercept-only plan. (4) Random passage intercepts and independent random condition slopes, integrating over the passage distribution. All four primary fits leave the reader-by-passage effect out. In a sensitivity check on the first 64 replications of two and five passages with 48 readers, heterogeneous slopes and local dependence, add a reader-by-passage random intercept to model 2 using the native kernel directly. This internal fit is a research comparison. Use joint-delta 95% intervals with zero bootstrap refits throughout. The native random-slope model estimates the passage slope SD with the existing regularization.

Use 44 conditions and 256 replications per condition in the medium run (master seed `91642001`), after an eight-replication technical pilot (seed `91641001`). Within a stage use independent streams across replication IDs and shared draws across scenarios where their design permits. Report valid estimates, bias, width, empirical SD, coverage against both truths, and exact pointwise 95% Monte Carlo intervals. The primary coverage target for the fixed fits is the sampled-passage contrast, and for the random fits it is the population-passage contrast. Use the §4B gates for valid-fit coverage (90% floor), null detection (7.5% upper tolerance for the applicable null), and power (80%). Detection uses all scheduled fits. Estimates that exclude zero under a population null with heterogeneous slopes can reflect a nonzero sampled-passage contrast.

The technical pilot completed 352 datasets and 1,440 valid fits. The medium run used `glmmplanr` 0.1.0 under R 4.6.1 on arm64 macOS. It completed 11,264 datasets and 45,312 valid fits in 22.7 seconds with six workers. A separate aggregation of the saved fits reproduced coverage, detection, and the exact binomial intervals for all 180 scenario-model combinations. The design-matrix rank check for two shared passages gave 3/3 for the additive fixed model and 4/4 with the interaction. Assigning a different passage to each arm gave rank 2/3 for the additive model. The corresponding earlier one-passage-per-arm simulation gave 420/1,024 null detections (41.0%) and 622/1,024 planned intervals covering the population contrast (60.7%). [Allocation benchmark](#4b-passage-allocation-at-fixed-response-count)

With 48 readers and both the passage-slope and reader-by-passage SDs set to 0.4, the results below show which truth each interval addresses. Percentages in brackets are exact pointwise 95% Monte Carlo intervals.

| Shared passages | Sampled-passage coverage, fixed interaction, +10 pp population effect | Population null detection, random passage intercept only | Population null detection, with random passage slopes | Population coverage, with random passage slopes, +10 pp |
|---:|---:|---:|---:|---:|
| 2 | 244/256 (95.3% [92.0–97.6]) | 63/256 (24.6% [19.5–30.4]) | 12/256 (4.7% [2.4–8.0]) | 245/256 (95.7% [92.4–97.8]) |
| 3 | 238/256 (93.0% [89.1–95.8]) | 42/256 (16.4% [12.1–21.5]) | 11/256 (4.3% [2.2–7.6]) | 244/256 (95.3% [92.0–97.6]) |
| 4 | 238/256 (93.0% [89.1–95.8]) | 34/256 (13.3% [9.4–18.1]) | 15/256 (5.9% [3.3–9.5]) | 244/256 (95.3% [92.0–97.6]) |
| 5 | 241/256 (94.1% [90.5–96.4]) | 31/256 (12.1% [8.4–16.7]) | 10/256 (3.9% [1.9–7.1]) | 242/256 (94.5% [91.0–97.0]) |
| 6 | 244/256 (95.3% [92.0–97.6]) | 30/256 (11.7% [8.0–16.3]) | 14/256 (5.5% [3.0–9.0]) | 238/256 (93.0% [89.1–95.8]) |

- **The same interval can perform very differently against the two truths.** For two shared passages at +10 pp, the fixed interaction covered the sampled-passage average in 244/256 datasets and the population average in 197/256 (77.0%, 95% MC interval 71.3–82.0%). At 96 readers, the two rates were 242/256 (94.5%) and 175/256 (68.4%). The sampled-passage truth varied across datasets with SD 5.59 pp under this design, regardless of reader count. For equal exposure to the same questions, adding the fixed interaction changed sampled-passage coverage by at most 1.2 points across the medium-run design grid. This comparison concerns the average across passages.
- **Accounting for heterogeneous passage slopes changed population inference.** At two shared passages and +10 pp, adding the random slope raised population coverage by 16.8 points (paired 95% conservative MC interval 10.2–22.7) and increased mean interval width from 18.26 to 28.65 pp. Across 2–6 passages in the table, detection of +10 pp with random slopes ranged from 77/256 to 118/256 (30.1–46.1%). Doubling readers for two passages raised detection to 102/256 (39.8%, 95% MC interval 33.8–46.1%) and gave 24/256 population null detections (9.4%, 6.1–13.6%). The null gate at 48 readers passed only for five passages in the table. That result pertains to the specified DGP and gate. Coverage at six passages under +10 pp was indeterminate against the 90% floor.
- **The common-effect scenario explains when the simpler model behaved differently.** With slope SD 0 and reader-by-passage SD 0.4, random passage intercepts gave 3.9–5.1% null detection and 94.9–96.1% planned population coverage across 2–6 passages. Fitting random slopes under true slope SD 0 yielded median regularized SD estimates of 0.23–0.28 in this profile. In the reader-by-passage sensitivity check, adding the random intercept for each reader–passage pair changed mean width by 0.09–0.21 pp and gave the same sampled-passage coverage decisions in the four sets of 64 pairs. These checks cover balanced complete presentation at the specified variances.

For reporting, state whether the interval refers to the presented passages or averages over passages from a population. For the latter, state how the planning model represents condition effects that vary by passage. Here valid native regularized fits alone gave no indication of the large population coverage losses in the random-intercept model. The benchmark provides an interpretation example for the existing 0.1.0 output.

To reproduce, run `sh validation/run-isolated.sh few_shared_passages` for the pilot. With the same isolated `R_LIBS`, run `Rscript --vanilla validation/few_shared_passages.R OUT medium WORKERS` in the repository root, with a separate existing output directory. Interrupted runs with the same code and settings resume at the next incomplete scenario. Verify the archive and generate paired comparisons with `Rscript --vanilla validation/summarize_few_shared_passages.R OUT medium`. Paired rate-difference intervals use exact binomial limits on the two discordance rates with a Bonferroni adjustment. [Simulation](few_shared_passages.R) · [Audit](summarize_few_shared_passages.R) · [Scenarios](few_shared_passages_medium_scenarios.csv) · [Summary and exact intervals](few_shared_passages_medium_summary.csv) · [Paired comparisons](few_shared_passages_medium_paired.csv) · [Trial-level fits](few_shared_passages_medium_fits.csv.gz)

### 4I. Shared-passage comparisons with external IRT and GLMM estimators

This independent-seed confirmation asks whether the loss of population coverage without passage-specific condition slopes also occurs with a second fitting engine. It also checks whether an external Rasch testlet model recovers the condition contrast for the actual test form. Keep the shared-passage binary response generator from §4H. At 48 readers, compare 2, 5 and 12 passages with 24 items total, slope and reader-by-passage SDs both 0.4. Use two and five passages for the homogeneous-slope and 96-reader checks, and two passages for the zero local-dependence check. The null and planned marginal population contrasts are zero and ten percentage points in all eight profiles. Use 1,024 datasets per condition (seed 91643001) with new readers, passages and questions in each; use independent seed 91643501 for the eight-dataset technical pilot. Match response vectors across fitters. Retain all scheduled fits in detection denominators and use valid intervals in coverage denominators.

For the population contrast, fit both random-intercept and independent passage-random-slope models under the native regularized Laplace estimator and under unpenalized `glmmTMB`. Both engines integrate over new readers, passages and items, with the same omission of reader-by-passage variation in these fits. For `glmmTMB`, propagate uncertainty from its joint fixed-effect and log-SD covariance to the same marginal probability difference. Record optimizer and Hessian failures, and flag any estimated SD below 0.001. A boundary fit with a valid interval contributes to the primary `glmmTMB` summaries. Native fitting and its regularization are as in §4H. All reported intervals are Wald on the percentage-point scale with zero bootstrap refits. Keep the §4H coverage, null detection and power gates and report exact pointwise 95% Monte Carlo intervals.

For the fixed test form, use the first 512 paired datasets from each of the four 48-reader, two-passage conditions with heterogeneous passage slopes, crossing the local SDs and the two population contrasts. Condition on the realized item difficulties and passage effects when computing truth. Average with equal weight per passage and equal weight per item within each passage, over new readers and their reader-by-passage variation. Fit fixed item intercepts, fixed condition slopes by passage, and normally distributed reader and independent reader-by-passage effects. The native and `glmmTMB` fits use common reader-by-passage variance. In `TAM`, represent each reader's assigned arm with linked A/B item columns and a shared item intercept, and use a Q-matrix for the general ability and the two orthogonal passage factors. Constrain the passage variances to be equal. Use TAM as a point-estimate reference and report bias and RMSE against the fixed-form truth. Its standard error routine omits covariances needed for the joint marginal-contrast interval. Under the population null with heterogeneous slopes, the true form contrast can be nonzero, so call zero exclusion by the fixed-form intervals detection. [TAM fixed-design MML](https://alexanderrobitzsch.github.io/TAM/reference/tam.mml.html) · [TAM standard errors](https://alexanderrobitzsch.github.io/TAM/reference/tam.se.html) · [glmmTMB reference](https://glmmtmb.github.io/glmmTMB/reference/glmmTMB.html)

The Rasch testlet fits use 1,024 adaptive quasi-Monte Carlo integration points in TAM. The technical pilot checks the transformation of the full joint fixed-effect and log-SD covariance from glmmTMB by independent integration and numerical derivatives, and checks the sensitivity of eight TAM fits to using 4,096 points. A second IRT fit uses `mirt::mixedmirt()` with the same fixed item intercepts, fixed condition effects by passage, and common reader-by-passage variance. Use the first 128 datasets per fixed-form case and compare the other models on this same subset. Use the MH-RM algorithm with a maximum of 1,200 cycles and separate L'Ecuyer streams (master seed 91646300 for the main run). Transform the two estimated variance parameters to log SD before computing the joint delta interval. Report point validity and interval validity separately. [mixedmirt model and estimation](https://philchalmers.github.io/mirt/docs/reference/mixedmirt.html)

**Population results.** All 16,384 datasets and all 65,536 scheduled fits finished. The native engine returned all 32,768 intervals. `glmmTMB` returned 16,376/16,384 intercept-only intervals and 16,307/16,384 slope intervals. The remaining 85 fits had a non-positive-definite Hessian or an invalid delta variance. The estimates below come from the 48-reader, heterogeneous-slope and locally dependent profile. Coverage uses the +10 pp condition. Values in brackets are exact pointwise 95% Monte Carlo intervals.

| Shared passages | Native null detection, intercept → slope | glmmTMB null detection, intercept → slope | Native planned coverage, slope | glmmTMB planned coverage, slope |
|---:|---:|---:|---:|---:|
| 2 | 193/1,024 (18.8%) → 59/1,024 (5.8% [4.4–7.4]) | 189/1,024 (18.5%) → 128/1,024 (12.5% [10.5–14.7]) | 966/1,024 (94.3% [92.7–95.7]) | 885/1,010 (87.6% [85.4–89.6]) |
| 5 | 114/1,024 (11.1%) → 51/1,024 (5.0% [3.7–6.5]) | 109/1,024 (10.6%) → 70/1,024 (6.8% [5.4–8.6]) | 973/1,024 (95.0% [93.5–96.3]) | 946/1,021 (92.7% [90.9–94.2]) |
| 12 | 84/1,024 (8.2%) → 55/1,024 (5.4% [4.1–6.9]) | 79/1,024 (7.7%) → 58/1,024 (5.7% [4.3–7.3]) | 968/1,024 (94.5% [93.0–95.8]) | 965/1,024 (94.2% [92.6–95.6]) |

With two passages, the slope estimate under `glmmTMB` had a median SD of 0.16 logit under the null versus 0.36 for the native regularized fit (true SD 0.4). The latter used a mean null interval 5.12 points wider on valid pairs (paired 95% MC interval 5.06–5.18). At 12 passages both engines passed the coverage floor and null tolerance with slopes. The planned interval widths were 19.70 and 19.46 points, respectively. The corresponding power was 533/1,024 (52.1%, 48.9–55.1%) and 560/1,024 (54.7%, 51.6–57.8%). Both slope analyses fell short of the 80% power target in all eight +10 pp cases.

At two passages and 96 readers, the slope-model null detection rose to 120/1,024 (11.7%, 9.8–13.8%) for the native estimate and 208/1,024 (20.3%, 17.9–22.9%) for `glmmTMB`. The homogeneous-slope 48-reader profiles gave null detection of 5.3–6.3% (native) and 5.2–6.2% (`glmmTMB`) using intercept-only fits. Removing local dependence in the heterogeneous two-passage profile still gave 20.8% and 20.3% null detection with only intercepts. These checks put the excess population null detection in context: which variance is omitted, how many passages are sampled, and how the fitting and interval methods handle few levels all matter.

The population fits on arm64 macOS took roughly 7.3 minutes of elapsed time with eight workers. Summed process CPU for the 16,384 intercept fits was 48.7 s native and 1,420.0 s `glmmTMB`; the slope fits took 59.2 s and 1,751.2 s, respectively. These are times for the full recorded fit calls, including the external formula setup and joint interval, on this design grid. [Population figure](shared_passage_irt_population.png) · [Population summary](shared_passage_irt_large_population_summary.csv) · [Population paired contrasts](shared_passage_irt_large_population_paired.csv)

**Fixed-form results.** The fixed-form comparison completed 2,048 paired datasets and 6,144 native, TAM and glmmTMB fits in about 8.3 minutes with six workers. The native and `glmmTMB` fits returned intervals on all 2,048 datasets. `TAM` returned 1,978/2,048 valid points (96.6%), with 38 iteration-limit and 32 fitting errors. The subset for `mirt` used 512 of those same datasets. `mirt` returned 465/512 valid points (90.8%) and 379/512 valid intervals (74.0%) in 5.5 minutes with eight workers. The 133 interval failures comprise 86 unusable joint covariance estimates, 43 nonconverged fits and 4 fitting errors. These elapsed times overlap.

| Local SD | Population effect in DGP (pp) | All four valid points on the mirt subset | Common-set RMSE in pp, native / TAM / glmmTMB / mirt | mirt valid intervals | mirt coverage among valid intervals |
|---:|---:|---:|---:|---:|---:|
| 0.4 | 0 | 108 | 4.52 / 4.55 / 4.55 / 4.51 | 85/128 | 78/85 (91.8% [83.8–96.6]) |
| 0.4 | 10 | 114 | 4.46 / 4.50 / 4.50 / 4.53 | 93/128 | 87/93 (93.5% [86.5–97.6]) |
| 0 | 0 | 111 | 4.63 / 4.65 / 4.66 / 4.66 | 100/128 | 91/100 (91.0% [83.6–95.8]) |
| 0 | 10 | 117 | 4.35 / 4.38 / 4.38 / 4.32 | 101/128 | 96/101 (95.0% [88.8–98.4]) |

The four point estimates had closely matched RMSE on their common valid datasets. At local SD 0.4 and a +10 pp population effect, the native fixed-form intervals covered 476/512 (93.0%, 95% MC interval 90.4–95.0%) and glmmTMB covered 478/512 (93.4%, 90.8–95.4%). Across all four cases, the native coverage was 92.4–93.0%, and glmmTMB coverage was 92.6–93.4%. Coverage conditional on valid intervals in the mirt subset has wider Monte Carlo uncertainty. With true local SD zero, the estimated median was about 0.22 under TAM's QMC fit and near zero under glmmTMB. Conclusions about the fixed-form average are limited to that estimand. The two IRT fits assess reader-by-passage dependence and condition effects for the presented items. Population-interval evidence in this benchmark comes from the matched native and glmmTMB fits.

The 1,024- versus 4,096-node sensitivity check for the eight pilot TAM datasets changed the estimated condition contrast by at most 0.33 pp. An initial coarse fixed integration grid was replaced before the main bank comparison because it treated small passage variances poorly. The original §4H pilot scenarios and fit rows reproduced exactly after extracting the shared generator. For case 1, replication 1 of the main stage, separate repeat fits from each fixed-form branch reproduced the saved results to a relative tolerance of 1e-10, excluding CPU time.

The stage scripts use the recorded versions of glmmplanr 0.1.0, glmmTMB 1.1.15.2, TAM 4.3.25 and mirt 1.48 with R 4.6.1 on arm64 macOS. Run the commands from the repository root with an existing output directory and the isolated installation on `R_LIBS`.

```sh
Rscript --vanilla validation/shared_passage_irt_comparison.R OUT large 8 population
Rscript --vanilla validation/shared_passage_irt_comparison.R OUT large 6 bank
Rscript --vanilla validation/shared_passage_irt_comparison.R OUT large 8 mirt
Rscript --vanilla validation/summarize_shared_passage_irt.R OUT large population
Rscript --vanilla validation/summarize_shared_passage_irt.R OUT large bank
Rscript --vanilla validation/summarize_shared_passage_irt.R OUT large mirt
Rscript --vanilla validation/plot_shared_passage_irt.R OUT
```

For a technical pilot, use `pilot` in place of `large`. Each branch keeps completed blocks in an output directory specific to its stage and branch. The mirt summary uses matching rows from the bank archive. [Comparison script](shared_passage_irt_comparison.R) · [Audit and paired summaries](summarize_shared_passage_irt.R) · [Population trial results](shared_passage_irt_large_population_fits.csv.gz) · [Fixed-form summary](shared_passage_irt_large_bank_summary.csv) · [Fixed-form trial results](shared_passage_irt_large_bank_fits.csv.gz) · [mirt matched-subset summary](shared_passage_irt_large_mirt_summary.csv) · [mirt trial results](shared_passage_irt_large_mirt_fits.csv.gz)

### 4J. Fixed-bank 2PL calibration with GLMM and IRT estimators

This benchmark asks whether the reduced-rank GLMM and specialist 2PL fits recover the same item response curves, and how a common-discrimination Rasch calibration behaves when true discrimination varies. Generate complete binary responses on four fixed-bank layouts (120 × 16, 240 × 8, 240 × 16, and 480 × 8). The eight-item bank is a common subset of the 16-item bank. Set latent ability to standard normal and difficulty between −1.8 and 1.8. Use three discrimination profiles: all 1, four values log-spaced from 0.75 to 1.35, and four values log-spaced from 0.55 to 1.80. Balance the discrimination values across positive and negative difficulty. Each replicate draws new readers and responses. All models use identical responses, and the cases share latent draws and response uniforms where their design permits. Use an eight-replication technical pilot (seed `91651001`) and 128 independent replications per case in the medium run (seed `91652001`). Use a separate L’Ecuyer stream per replication for the random residual initialization of `glmmTMB`, with master seed 100 above the data seed.

Fit the rank-one logit GLMM (`glmmTMB::glmmTMB()` with `rr(0 + item | person, d = 1)`), unidimensional `TAM::tam.mml.2pl()` and `mirt::mirt()` 2PL. Also fit a GLMM with discriminations supplied at their generating values (`lme4::glmer()` with Laplace integration), estimating its ability SD. This known-discrimination fit assesses the information from known relative loadings. The two common-slope models are the internal `mml` engine used in `plan_rasch_fixed_bank()` and the corresponding `mirt` Rasch fit. Both estimate ability variance under a common raw discrimination of 1. Convert every estimate to the mean-zero, unit-SD ability scale and derive difficulty as `-d/a`. Reflect the common latent direction if needed. A fit with a mixed-sign or effectively zero loading cannot supply a complete positively oriented bank in this comparison. Record optimizer failures, nonpositive loadings, constant-response items, and the glmmTMB Hessian diagnostic. Point validity uses convergence and a complete positive-loading bank; the Hessian diagnostic is tabulated separately.

Record discrimination and difficulty RMSE, worst-item absolute difficulty error, response-curve RMSE integrated over standard normal ability, and a common 61-point Gauss–Hermite marginal log likelihood. The 61-point likelihood evaluates the fitted banks under the same ability distribution. Report difficulty criteria on the eight common items as well as the entire presented bank, so comparisons of layouts do not depend only on how many items enter the maximum. The full-bank descriptive attainment criterion is difficulty RMSE ≤ 0.35 and maximum absolute difficulty error ≤ 0.60 on the standardized ability scale. The joint criterion for the three free 2PL fits adds discrimination RMSE ≤ 0.30. Use exact pointwise 95% Monte Carlo intervals over all scheduled fits for attainment and validity. The design comparison threshold is a lower limit of 80% for attainment. Report paired differences over the same response datasets. No interval refits are scheduled.

As a separate prediction check, compute the change in expected mean item score when the ability mean moves from 0 to 0.4 while its SD stays 1. Average both predicted scores over the same items. This is a projection from the calibrated test curves on a specified scale. The data generator contains one calibration group. [TAM 2PL](https://alexanderrobitzsch.github.io/TAM/reference/tam.mml.html) · [mirt](https://philchalmers.github.io/mirt/docs/reference/mirt.html) · [glmmTMB reduced rank](https://glmmtmb.github.io/glmmTMB/articles/covstruct.html#reduced-rank)

**Results.** The medium run completed 1,536 datasets, 9,216 fits and 110,592 item rows under R 4.6.1 on arm64 macOS. The native Rasch, mirt Rasch and known-discrimination fits each gave 1,536 complete banks. The free fits gave 1,476 (glmmTMB), 1,488 (TAM) and 1,486 (mirt) complete positively oriented banks. There were 60 nonpositive-loading banks in glmmTMB, and 17 each in TAM and mirt. The latter two had 31 and 33 optimizer failures, respectively. The first medium pass took 306.7 elapsed seconds with eight workers. Only the 1,536 glmmTMB fits were repeated with the explicitly seeded start, in 173.4 seconds with six workers, using the original responses. The archived results use these seeded glmmTMB fits and the other five original fits.

- **Matched Rasch fits.** On the unit-SD ability scale, the maximum difference in difficulty between native and mirt Rasch fits across the whole grid was 0.00082. With all true discriminations equal to one at 240 × 8, both met the difficulty criteria in 122/128 runs (95.3%, exact 95% MC interval 90.1–98.3%). The free TAM and mirt 2PL fits each met the same difficulty criteria in 66/128 runs (51.6%, 42.6–60.5%). At 240 × 16 under moderate heterogeneity, median curve RMSE was 4.38 pp under the native Rasch fit and 4.44 pp under mirt 2PL. This setting illustrates the sampling cost of estimating a slope for every item.
- **Specialist 2PL fits.** TAM and mirt had closely matched fitted banks when both were valid. Across matched items, their median absolute differences were about 0.000003 for discrimination and 0.000007 for difficulty. With wide heterogeneity at 480 × 8, both met the difficulty criteria in 90/128 (70.3%, 61.6–78.1%) and the joint difficulty plus discrimination criterion in 68/128 (53.1%, 44.1–62.0%). The known-discrimination GLMM met the difficulty criterion in 127/128 (99.2%, 95.7–100.0%). That comparison assumes the relative loadings are known at their generating values. Across the grid, no free 2PL fit passed the 80% attainment criterion with its pointwise Monte Carlo lower limit.
- **Prediction with wide discrimination heterogeneity.** Each layout below has 3,840 responses. Values are medians of the root mean squared error in item-level correct-response probabilities, integrated over standard normal ability and reported in percentage points. They use complete fits for each method. The native Rasch fit met the full-bank difficulty criterion in 0/128 trials for both layouts. In the 480 × 8 layout, its RMSE for the projected whole-test score gap following a +0.4 ability shift was 0.35 pp, compared with 0.34 pp under mirt 2PL. This distinction reflects the different targets of the curve and the aggregate score projections.

| Layout | Native Rasch | glmmTMB rank 1 | TAM 2PL | mirt 2PL | GLMM with known relative discrimination |
|---|---:|---:|---:|---:|---:|
| 240 × 16 | 6.68 | 4.46 | 4.47 | 4.47 | 3.28 |
| 480 × 8 | 5.75 | 29.75 | 3.61 | 3.61 | 2.33 |

**Likelihood diagnosis.** For 240 × 8 and 480 × 8 under wide heterogeneity, glmmTMB returned a loading above 5 in 78 and 105 complete fits out of 128 scheduled trials, respectively. The largest true loading was 1.80. For the first such fit in each case with a matched valid mirt fit, refit the same RR model using its saved seed and compare the likelihoods at the two banks. In the 240 × 8 example, Laplace integration favored the glmmTMB solution by 8.73 log-likelihood units, while direct one-dimensional integration favored mirt by 9.87. The largest estimated loading was 23.37 versus 1.68. In the 480 × 8 example the two margins were 74.17 in favor of glmmTMB under Laplace and 425.08 in favor of mirt under direct integration. This check shows an approximation discrepancy in two selected datasets. Use these cases to interpret the extreme estimates. The paired 61-point likelihood summaries use common valid pairs with all fitted loadings at most 5 because direct integration detected insufficient quadrature resolution at extreme loadings. The attained-rate criteria always use all scheduled fits. [Likelihood check script](check_two_pl_likelihood.R) · [Likelihood check results](two_pl_fixed_bank_medium_likelihood_check.csv)

The evidence supports using the native Rasch estimator for the benchmarked common-discrimination calibration and treating pilots with free discrimination as a distinct design problem. A rank-one covariance gives the 2PL response structure. The achieved error under that structure also depends on the integration method and amount of within-person item information. Whole-test expected-score accuracy under a known ability shift addresses a different goal from calibrating each item or verifying condition-effect power.

To reproduce, run `sh validation/run-isolated.sh two_pl_fixed_bank` for the pilot with the external packages installed. Start a medium run with `Rscript --vanilla validation/two_pl_fixed_bank.R OUT medium WORKERS` in the repository root with the isolated package library on `R_LIBS`. Use `Rscript --vanilla validation/summarize_two_pl_fixed_bank.R OUT medium` to verify the item-based error calculations and aggregate the saved fits. Run `Rscript --vanilla validation/check_two_pl_likelihood.R OUT` under the same isolated package to reproduce the numerical diagnosis. [Benchmark script](two_pl_fixed_bank.R) · [Truth](two_pl_fixed_bank_medium_truth.csv) · [Summary and exact intervals](two_pl_fixed_bank_medium_summary.csv) · [Paired contrasts](two_pl_fixed_bank_medium_paired.csv) · [Fit outcomes](two_pl_fixed_bank_medium_fits.csv.gz) · [Item parameters](two_pl_fixed_bank_medium_items.csv.gz) · [Environment](two_pl_fixed_bank_medium_environment.txt)

### 4K. Complete and incomplete presentation at equal item exposure

This follow-up holds the 16-item bank and every item's exposure fixed to distinguish the number of responses from the number of item pairs observed per reader. Compare 240 readers answering all 16 items with 480 readers each answering 8. Each design provides 3,840 responses and 240 observations per item. For incomplete presentation use 30 booklets repeated 16 times each. Define the 16 items as the points of a four-dimensional binary vector space. For each of the 15 nonzero binary directions, take the two eight-item sets with inner product 0 or 1, modulo 2. Each item is in 15 booklets and each item pair in seven. The resulting mask gives exactly 112 co-observations per pair, compared with 240 under full presentation. Assignment cycles through this fixed design independently of the generated abilities and responses.

Use the same bank and three discrimination profiles as §4J (1 throughout, 0.75–1.35, 0.55–1.80). The complete arm reuses its 240 × 16 archives, including all 128 replications per profile. Generate each incomplete replicate from the same L'Ecuyer stream (master seed 91652001), drawing the same abilities and response uniforms so that the first 240 readers share potential responses before taking the booklet mask. The complete and incomplete arms both average any error over the same 16 item identities. Fit the native Rasch and mirt Rasch models and the mirt and TAM 2PL models as in §4J. Apply the same complete-positive-bank criterion, difficulty RMSE and maximum absolute error criteria, and joint 2PL criterion. No interval refits are scheduled. Report per-arm exact 95% Monte Carlo intervals and paired differences in attainment with conservative pointwise 95% intervals based on the two discordance rates. For interpretable measures of accuracy, retain the bank and item-curve RMSE medians and the projected mean-score contrast error. Compare bank validity and numerical failures separately.

Use the public irtsim 0.2.0 workflow in a separate 128-replicate-per-layout, wide-discrimination check (master seed 91654, within the integer range of its serial seed calculation). Supply the same 16 true item parameters and the same 30-row booklet matrix. Check the exact masking schedule and derive the same full-bank criteria from its per-iteration item results. It generates a separate response sample and uses its own mirt defaults. This independently checks the design and the usability of the public irtsim workflow. [irtsim reference](https://CRAN.R-project.org/package=irtsim/refman/irtsim.html)

The new arm completed 384 datasets and 1,536 fits in 40.1 seconds with six workers under R 4.6.1 on arm64 macOS. All 768 Rasch fits returned a complete bank. The incomplete designs returned complete positively oriented 2PL banks in 371/384 fits for mirt and 376/384 for TAM. A separate aggregation checked all 24,576 fitted item rows, the identity of the bank across the two layouts, and the paired replication keys.

| Generating discrimination | Joint 2PL attainment, 240 × 16 complete | Joint 2PL attainment, 480 × 8 over 16 items | Paired change in attainment (percentage points, 95% MC interval) |
|---|---:|---:|---:|
| Equal | 70/128 (54.7%) | 37/128 (28.9%) | −25.8 [−42.0, −8.0] |
| Moderate heterogeneity | 53/128 (41.4%) | 28/128 (21.9%) | −19.5 [−34.9, −3.0] |
| Wide heterogeneity | 27/128 (21.1%) | 5/128 (3.9%) | −17.2 [−27.6, −5.4] |

The joint counts and paired intervals coincide for mirt and TAM. For wide heterogeneity, the exact 95% Monte Carlo intervals for the two attainment rates are 14.4–29.2% and 1.3–8.9%. For mirt, median discrimination RMSE increased from 0.279 to 0.370, and median item-curve RMSE increased from 4.47 to 5.03 percentage points on valid banks. Across the 114 pairs valid under both layouts, the median change in discrimination RMSE was +0.077 (order-statistic 95% MC interval +0.054 to +0.118). The difficulty-only criterion was met in 36/128 complete runs and in 23/128 mirt and 24/128 TAM booklet runs. The paired intervals for difficulty-only attainment included zero in this profile. The projected whole-bank score contrast for a +0.4 ability shift had mirt RMSE of 0.36 and 0.35 percentage points, respectively.

When the true discrimination was common and the analysis used a Rasch model, difficulty attainment was 124/128 with complete presentation and 125/128 with booklets (paired difference +0.8 percentage points, 95% MC interval −5.8 to +7.3) under both Rasch estimators. The data support a change in the precision of free discrimination for the 2PL designs here. The Rasch difficulty attainment showed little observed change for the common-discrimination bank under the tested presentation schedules.

The public irtsim 0.2.0 workflow honored the same booklet matrix and the exact 240-per-item allocation. Its separate wide-profile run gave 32/128 joint successes for the complete form (25.0%, 17.8–33.4%) and 10/128 for the booklets (7.8%, 3.8–13.9%). There were 126/128 and 113/128 valid complete banks, respectively. The package provides per-item summaries and per-iteration item estimates. The full-bank success rate was calculated by evaluating the joint error criterion on the item estimates from each iteration and retaining invalid iterations as misses. This comparison uses irtsim's public design and simulation interface and mirt backend. Differences from the paired experiment include its response sample and optimizer defaults. [irtsim run summary](two_pl_fixed_bank_booklets_irtsim_summary.csv) · [irtsim per-iteration results](two_pl_fixed_bank_booklets_irtsim_fits.csv)

The planning result concerns how the same item exposure is arranged across people. In the fixed-bank 2PL example, candidate designs should specify both item exposure and how many items each person answers, with the actual overlap between booklets. The public irtsim workflow already supports the booklet layout and exposes the relevant per-iteration estimates. This provides an external planning route and a benchmark for any future design-search integration. The native Rasch estimate in these custom booklets uses the internal validation entry point. The benchmark reports accuracy and attainment at the tested designs.

To reproduce from the repository root, set `R_LIBS` to an isolated glmmplanr installation and pass an existing output directory. The paired aggregation reads §4J's complete trials from the second directory.

```sh
Rscript --vanilla validation/two_pl_fixed_bank.R OUT booklets WORKERS
Rscript --vanilla validation/summarize_two_pl_fixed_bank.R OUT booklets
Rscript --vanilla validation/compare_two_pl_booklets.R OUT validation
Rscript --vanilla validation/check_two_pl_booklet_irtsim.R OUT
```

[Booklet matrix](two_pl_fixed_bank_booklets_forms.csv) · [Parameters and layouts](two_pl_fixed_bank_booklets_truth.csv) · [Booklet summary](two_pl_fixed_bank_booklets_summary.csv) · [Paired layout differences](two_pl_fixed_bank_booklets_layout_differences.csv) · [Fit outcomes](two_pl_fixed_bank_booklets_fits.csv.gz) · [Item estimates](two_pl_fixed_bank_booklets_items.csv.gz)

### 4L. Rasch calibration, fit, separation, and anchored DIF

This benchmark asks how sample size and items per person affect three
decisions: whether 16 item difficulties can be recovered, how well the
same responses support item diagnostics and person measurement, and
whether a prespecified item shows uniform group DIF. It compares
`glmmplanr` 0.1.0, `TAM` 4.3.25, and `mfrmr` 0.2.4 on each simulated dataset.

- Calibration uses 120, 240, or 480 persons answering all 16 items, or
  480 answering eight items in each of 30 balanced booklets. At equal
  240 responses per item, the complete and booklet layouts give 240 and
  112 joint observations per item pair. The three Rasch profiles use
  evenly spaced difficulties from -1.8 to 1.8 and a mean-zero normal
  ability distribution. Ability SD is 0.8 in the targeted profile and
  0.4 in the narrow profile. The off-target profile has SD 0.8 and
  difficulties shifted up by 1.4 logits. Two further complete layouts
  at 120 and 480 persons set discrimination to 0.35 for item i08,
  leaving all other items at 1. All 14 conditions use 512 replications
  with seed 91661001.
- The bank criterion requires item RMSE at most 0.35 logits and maximum
  absolute item error at most 0.60 logits in the same replication.
  Common diagnostics use 41-point quadrature for EAP estimates and
  posterior variances, EAP plug-in fit residuals, and the observed
  information (Louis identity) for item marginal SEs. EAP reliability
  is the sample variance of EAP estimates divided by that variance plus
  the mean posterior variance. Item separation is the square root of
  positive-part error-adjusted sample variance divided by mean squared
  item SE. Item reliability is that adjusted variance divided by sample
  variance. Infit and outfit above 1.2 are descriptive screening events
  for the prespecified item i08. The first 32 replications of each
  condition also record TAM's native item and person fit and EAP
  reliability and mfrmr's legacy fit and model-SE separation and
  reliability. These native diagnostics use their own definitions.
- The DIF benchmark uses 120 and 480 persons answering all items and
  480 in the eight-item booklets. Groups A and B are balanced within each
  form. Both have residual ability SD 0.8. Group B has a mean ability
  advantage of 0.6 logits. Item i16 has a group-B conditional log-odds
  shift of 0 or +0.7 logits. Two additional 480-person complete
  conditions use equal group means. All eight DIF conditions use 384
  replications with seed 91662001. In both fitted models, a group
  coefficient accounts for mean ability and 15 common items identify
  the scale. The alternative permits separate focal difficulties by
  group. The one-degree-of-freedom likelihood-ratio test uses nominal
  alpha 0.05. TAM and mfrmr use direct marginal likelihood for this
  comparison. The glmmplanr validation script uses the private Laplace
  GLMM engine with fixed items and a person random intercept. It
  records common 41-point likelihoods as a check on the approximation.
- A L'Ecuyer stream gives new persons and responses in each replication.
  Conditions within a replication share standard-normal and uniform
  draws, and all estimators receive the same observed response matrix
  within a condition. The native Rasch engine keeps its 200-iteration
  release limit. MML versions of all three models estimate the latent
  variance, including an intercept-only population model in mfrmr for
  calibration. The DIF GLMM uses 500 EM iterations and tolerance
  `2e-5`. The preliminary technical check uses a separate seed 91660001.
- The script records failure reasons and unavailability of diagnostics.
  Attainment and detection rates have scheduled denominators, with
  valid-result rates and valid-only rates also reported. The simulation
  uses exact 95% Monte Carlo intervals for proportions. Fit differences,
  null rejection, power, and precision are specific to these bank,
  population, and allocation settings.

Run from the repository root with the three packages installed. Supply
a fresh output directory for a new run. A restart in the same directory
reuses completed 32-replication batches for the same source, settings,
and versions.

```sh
Rscript --vanilla validation/rasch_fit_benchmark.R calibration /path/to/results 8
Rscript --vanilla validation/rasch_fit_benchmark.R dif /path/to/results 8
Rscript --vanilla validation/summarize_rasch_fit_benchmark.R /path/to/results
Rscript --vanilla validation/check_rasch_mml_limit.R /path/to/results
Rscript --vanilla validation/plot_rasch_fit_benchmark.R /path/to/results
```

#### Fixed-bank calibration and person measurement

The targeted profile gave the following results. Mean values use TAM
fits and common diagnostics. All three engines agreed on each bank
attainment decision for these four layouts.

| Persons × items per person | Item difficulty RMSE | Bank attainment (95% MC interval) | Person EAP reliability | Item reliability |
| --- | ---: | ---: | ---: | ---: |
| 120 × 16 | 0.231 | 409/512 (79.9%, 76.1–83.3%) | 0.649 | 0.960 |
| 240 × 16 | 0.161 | 509/512 (99.4%, 98.3–99.9%) | 0.648 | 0.980 |
| 480 × 16 | 0.114 | 512/512 (100.0%, 99.3–100.0%) | 0.648 | 0.990 |
| 480 × 8 | 0.162 | 503/512 (98.2%, 96.7–99.2%) | 0.479 | 0.980 |

For 480 × 8 minus 240 × 16, the paired difference in item RMSE was
+0.0011 logits (95% MC interval −0.0022 to +0.0044). Person EAP
reliability was lower by 0.169 (0.166–0.172) and person estimation
RMSE rose from 0.477 to 0.579 logits. At the same equal-exposure
comparison, EAP reliability was 0.330 vs. 0.196 with ability SD 0.4
and 0.592 vs. 0.422 with the harder bank. For the harder bank at
120 × 16, only 235/512 datasets met the item criterion (45.9%, 41.5–50.3%).

Across the 7,168 datasets, glmmplanr, TAM and mfrmr yielded 6,710,
7,163 and 7,167 valid calibrations. The off-target condition had one
dataset with an all-incorrect item, counted as non-attaining for all
engines. Four TAM fits reached the 1,000-iteration limit. In 16
further narrow-range datasets, TAM reported an ability variance of
0.001. Their re-evaluated log likelihoods were 0.013–0.280 below
mfrmr's on common quadrature. For jointly valid glmmplanr–mfrmr pairs,
item RMSE differed by at most 0.000029 logits and the bank decisions
agreed across all such pairs.

All 457 glmmplanr iteration-limit failures occurred at true ability
SD 0.4. In 480 × 8 the release engine returned 143/512 valid fits,
with 142/512 meeting the bank criterion. TAM and mfrmr met that
criterion in 506/512 and 510/512. Raising the limit to 5,000 in the
follow-up recovered 454/457 fits, including 367/369 in 480 × 8.
The 454 fits agreed with mfrmr on every bank decision. Their common
log likelihoods were within 0.000032 of mfrmr. The three remaining
fits corresponded to mfrmr SD estimates of 0.005, 0.024 and 0.090.
The release benchmark keeps all original failures in its denominators.
[Simulation-level fits](rasch_fit_benchmark/calibration_calibration_fits.csv.gz) ·
[Monte Carlo intervals](rasch_fit_benchmark/calibration_summary.csv) ·
[Matched contrasts](rasch_fit_benchmark/calibration_paired.csv) ·
[Iteration-limit follow-up](rasch_fit_benchmark/calibration_iteration_limit.csv) ·
[Calibration figure](rasch_fit_benchmark/calibration.png).

#### Package-reported diagnostics

Both packages returned the requested fit indices in the first 32
replications per condition. At 480 persons answering all 16 items,
changing i08 discrimination from 1 to 0.35 moved its mean infit and
outfit from 1.002 and 1.002 to 1.086 and 1.106 in TAM's
posterior-integrated `msq.itemfit()`. In mfrmr's legacy diagnostics the
same change moved the means from 0.911 and 0.894 to 0.997 and 0.993.
The 1.2 threshold flagged 0/32 datasets with either native index in
this 480-person low-discrimination condition. These residual
diagnostics have different reference behavior. A threshold such as
1.2 answers a limited screening question under the distribution and
definition used here. In the targeted Rasch profile, the fraction of
persons above the same threshold with common EAP plug-in infit averaged
13.6% at 240 × 16 and 17.8% at 480 × 8. These are screening frequencies
when the response model generated all observations.

At 240 × 16 with ability SD 0.4, mean TAM EAP reliability was 0.337.
The mfrmr legacy model-SE person reliability and separation both
averaged zero because the estimated error variance exhausted the
observed person-estimate spread. In the targeted 240 × 16 scenario,
the mean item and person separation statistics from mfrmr were 7.65
and 0.96. EAP reliability uses `v / (v + s)`. The legacy
variance-subtraction reliability uses `max(0, v - s) / v` with its
own estimates and SEs.
[TAM item fit](https://alexanderrobitzsch.github.io/TAM/reference/msq.itemfit.html) ·
[TAM reliability](https://alexanderrobitzsch.github.io/TAM/reference/WLErel.html) ·
[mfrmr documentation](https://cran.r-project.org/web/packages/mfrmr/refman/mfrmr.html).

#### Anchored uniform DIF results

TAM and mfrmr returned valid tests and identical decisions on all 3,072
datasets. When group B had a 0.60-logit mean advantage, their results were:

| Persons × items per person | Focal responses per group | Null rejection | Detection at +0.70 logit DIF (95% MC interval) | DIF estimate RMSE (logits) |
| --- | ---: | ---: | ---: | ---: |
| 120 × 16 | 60 | 18/384 (4.7%) | 137/384 (35.7%, 30.9–40.7%) | 0.514 |
| 480 × 16 | 240 | 22/384 (5.7%) | 329/384 (85.7%, 81.8–89.0%) | 0.238 |
| 480 × 8 | 120 | 20/384 (5.2%) | 213/384 (55.5%, 50.3–60.5%) | 0.351 |

With equal group means and 480 complete responders, the two MML fits
rejected 23/384 nulls (6.0%, 95% MC interval 3.8–8.9%) and detected
320/384 alternatives (83.3%, 79.2–86.9%). DIF is identified as a group-by-item
fixed effect on the common scale. At equal numbers of people, the booklet
design changes both the number of responses to the focal item and the
number of anchor responses informing each person's latent location.

The Laplace GLMM had 3,064 valid tests. In eight datasets, both fits were
finite and the alternative's approximate log likelihood was lower than
the null beyond the numerical tolerance. It rejected seven fewer
datasets in total, with all seven also rejected by direct MML. The
largest differences from TAM were 0.049 logits in estimated DIF and
0.31 in the likelihood-ratio statistic. This is a comparison for a
prespecified single item with invariant anchors and equal group variances.
[Simulation-level fits](rasch_fit_benchmark/dif_dif_fits.csv.gz) ·
[95% Monte Carlo intervals](rasch_fit_benchmark/dif_summary.csv) ·
[Matched contrasts](rasch_fit_benchmark/dif_paired.csv) ·
[DIF figure](rasch_fit_benchmark/dif.png).

### 4M. Paired replay with the mfrmr 0.3.2 release candidate

The exact candidate archive was `mfrmr_0.3.2.tar.gz` with SHA-256
`f0ebe1a9880e6ffaf41f96654408a9f9305604f7a95b120917cfc147c59ec4ab`.
It was installed into a temporary library. The replay uses the 4L response
generator, forms, seeds, contrasts, and denominators. The glmmplanr 0.1.0,
TAM 4.3.25, and mfrmr 0.2.4 fits come from the saved 4L files. The
replay compares the release candidate's default C++ MML with the same
31-node integration order and population formulas. The common 41-node
calculation checks likelihoods, common EAP-based diagnostics, and the
point estimates under the same definitions. The first 32 draws per
calibration condition again supply the native diagnostics.

The primary DIF comparison refits the null and the prespecified i16 split
model under the same grouping and likelihood criteria as 4L. A separate
public-API check runs `analyze_dff(method = "mml_lrt")` on replication 4 of
each of the eight conditions and records the i16 unadjusted and Holm-adjusted
results from the native 16-item screen. Its eight datasets check the fit and
test output on the same responses. The main replay supplies the estimates
of test size and power. Runtime comparisons use process CPU time for the
public fitting calls on the same machine. The median, mean, and 95th
percentile characterize the typical time and the effect of slow fits.
The candidate computes some diagnostic values during fitting, so the
calibration files also report the combined fit and diagnostic CPU time
for the planned diagnostic draws.

With both glmmplanr 0.1.0 and the specified mfrmr 0.3.2 archive installed,
run the replay and its comparison from the repository root.

```sh
Rscript --vanilla validation/rasch_fit_benchmark.R calibration /path/to/replay 10 mfrmr
Rscript --vanilla validation/rasch_fit_benchmark.R dif /path/to/replay 10 mfrmr
Rscript --vanilla validation/summarize_rasch_fit_benchmark.R /path/to/replay /path/to/replay validation/rasch_fit_benchmark
Rscript --vanilla validation/check_rasch_mfrmr_dff.R validation/rasch_fit_benchmark /path/to/replay
```

#### Calibration and CPU time

On the 7,168 matched calibration datasets, both mfrmr versions returned
7,167 valid fits and the same 6,404 bank attainments. The remaining dataset
contained an all-incorrect item and was retained in the scheduled
denominator. On valid pairs, the largest changes in item RMSE and
re-evaluated log likelihood were 0.000011 logits and 0.000010,
respectively. The largest fitted ability-SD change occurred in narrow-range
booklet replication 399, from 0.00538 to 0.00000615. Both versions met the
bank criterion on that dataset. In the narrow-range 480 × 8 condition the
attainment counts were 142/512 for glmmplanr, 506/512 for TAM, and 510/512
for each mfrmr version.

For the targeted calibration profile, fit CPU time per dataset was as follows.

| Persons × items per person | 0.2.4 median / 95th percentile (s) | 0.3.2 median / 95th percentile (s) |
| --- | ---: | ---: |
| 120 × 16 | 0.890 / 1.033 | 0.140 / 0.288 |
| 240 × 16 | 1.472 / 1.667 | 0.297 / 1.485 |
| 480 × 16 | 2.673 / 3.607 | 0.986 / 28.416 |
| 480 × 8 | 2.099 / 2.429 | 0.352 / 2.380 |

Across all calibration fits the total recorded CPU time was 13,301 seconds
in 0.2.4 and 11,809 seconds in 0.3.2. In the targeted 480 × 16 layout
the mean increased from 2.733 to 4.027 seconds as the upper tail grew.
[Candidate fits](rasch_mfrmr_0_3_2/calibration_calibration_fits.csv.gz) ·
[Four-engine summaries](rasch_mfrmr_0_3_2/calibration_summary.csv) ·
[Matched contrasts](rasch_mfrmr_0_3_2/calibration_paired.csv).

#### Native reliability and fit diagnostics

In the 448 matched diagnostic datasets, the candidate's model-based person
reliability agreed with the common EAP calculation within 0.00000031.
The following means use the same 32 diagnostic datasets per row.

| Calibration condition | Person reliability in 0.2.4 | Person EAP reliability in 0.3.2 |
| --- | ---: | ---: |
| Targeted, 240 × 16 | 0.477 | 0.658 |
| Targeted, 480 × 8 | 0.014 | 0.478 |
| Ability SD 0.4, 240 × 16 | 0.000 | 0.338 |
| Ability SD 0.4, 480 × 8 | 0.000 | 0.187 |

The release candidate reports the EAP variance divided by the sum of that
variance and the mean posterior variance. In 0.2.4 the corresponding
reported measure used error-variance subtraction. The focal i08 native
infit and outfit each differed by at most 0.000002 between versions
across the 448 matched diagnostic datasets.

**Winsteps definitions vs. `mfrmr` 0.3.2 diagnostics.** Per the Winsteps
and Facets manuals (`reliability.htm`, `table3_1.htm`) and Wright & Masters
(1982, pp. 92, 106; 2002, *RMT* 16:3, p. 888), Rasch separation, strata, and
reliability are defined on direct (unshrunk) logit measures using the
**population** observed variance $\text{P.SD}^2 = \frac{1}{N}\sum_{i=1}^N (\hat{\beta}_i - \bar{\hat{\beta}})^2$
(divisor $N$; sample counterparts `S.SD`, `S.Sep`, `S.Strata`, and
$\text{S.Rel} = \{1 + (N - 1)\text{P.Rel}\}/N$ use divisor $N - 1$),
$\text{Model RMSE} = \sqrt{\frac{1}{N}\sum_{i=1}^N \text{Model SE}_i^2}$, and
$\text{Real RMSE} = \sqrt{\frac{1}{N}\sum_{i=1}^N \text{Real SE}_i^2}$ with
$\text{Real SE}_i = \text{Model SE}_i\sqrt{\max(1, \text{INFIT MnSq}_i)}$:
$\text{True SD} = \sqrt{\max(0, \text{P.SD}^2 - \text{RMSE}^2)}$,
$\text{Separation } G = \text{True SD} / \text{RMSE}$,
$\text{Strata } H = (4G + 1)/3$, and
$\text{Reliability } R = \max(0, \text{P.SD}^2 - \text{RMSE}^2) / \text{P.SD}^2 = G^2 / (1 + G^2)$
(reported separately for non-extreme and extreme+non-extreme elements).
By contrast, `mfrmr` 0.3.2's `diagnose_mfrm()` (`.mfrm_diagnose_fit()`) departs
from Winsteps in three ways: (1) under `score_method = "EAP"` (its MML default),
it sets `v_true = stats::var(EAP)` (shrunk posterior-mean variance) and plugs
misfit-inflated `RealSE` into the additive EAP ratio
`RealReliability = v_obs / (v_obs + mean(RealSE^2))` and
`RealSeparation = sqrt(v_obs / mean(RealSE^2))`, which breaks the total-variance
identity and prevents `RealReliability` and `RealSeparation` from ever reaching
`0` (or `Strata` from reaching `1/3`) even when `Real RMSE >= Observed SD`;
(2) across all facets and methods (`JML` and `MML`), it uses sample variance
`stats::var()` (divisor $N - 1$, `S.SD^2`) alongside `mse = mean(s^2)` (divisor
$N$) without distinguishing Winsteps/Facets's default population definitions
(`P.SD^2`, divisor $N$), inflating separation, strata, and reliability for small
facet counts; and (3) it reports only `RealSeparation`-based strata under an
unlabeled `Strata` column without `ModelStrata`. In our 4L helper script
(`validation/rasch_fit_benchmark_helpers.R`), `person_eap_rel` is marginal EAP
reliability (`var(eap) / (var(eap) + mean(post_var))`), `person_eap_ratio` is
not interpreted as Winsteps separation, and `item_rel` / `item_separation` use
sample variance (`var(b)`, divisor $J - 1$, i.e. `S.Rel` / `S.Sep`) with MML
observed-information SEs.

#### Anchored DIF replay

The 3,072 matched tests were valid and the 1,082 rejections agreed
dataset by dataset between mfrmr versions and with the saved TAM results.
The largest changes in focal DIF and the likelihood-ratio statistic were
0.000020 logits and 0.000000004. Under 0.60-logit population impact and
a true DIF shift of +0.70 logits, both mfrmr versions detected DIF in
137/384 (35.7%) datasets at 120 × 16, 329/384 (85.7%) at 480 × 16,
and 213/384 (55.5%) at 480 × 8.

The table shows the measured CPU time for both fits per replicate in this
DIF alternative. Across all DIF fits the total CPU time was 12,997 seconds
for 0.2.4 and 15,495 seconds for 0.3.2.

| Persons × items per person | 0.2.4 median / 95th percentile (s) | 0.3.2 median / 95th percentile (s) |
| --- | ---: | ---: |
| 120 × 16 | 1.759 / 2.073 | 0.307 / 0.695 |
| 480 × 16 | 5.098 / 7.084 | 2.510 / 33.070 |
| 480 × 8 | 4.371 / 5.029 | 1.028 / 15.081 |

An independently seeded 480 × 16 pilot split fit took 26 CPU seconds
and recorded 250 L-BFGS iterations. R profiling attributed over 99% of
sampled time to the candidate's C++ fitting call for that example.
[Candidate fits](rasch_mfrmr_0_3_2/dif_dif_fits.csv.gz) ·
[Four-engine summaries](rasch_mfrmr_0_3_2/dif_summary.csv) ·
[Matched contrasts](rasch_mfrmr_0_3_2/dif_paired.csv).

#### Direct check of the public DIF API

All 16 item tests were available in each of the eight API-check datasets.
For the prespecified i16 contrast, the largest difference from the archived
0.2.4 refits was 0.000006 logits for the effect and below 0.000000001 for
the likelihood-ratio statistic. The candidate manual refits gave the same
results to numerical precision. Sign was aligned to express higher
response log odds in the focal group. In the 480-person booklet example
with planted DIF, the raw p-value was 0.0223 and the 16-item Holm-adjusted
p-value was 0.334. This comparison quantifies the distinction between
testing a nominated item and screening all items in that example.
[Public-API trial results](rasch_mfrmr_0_3_2/dif_public_api.csv).
