# glmmplanr 0.1.0

* Initial public release of `glmmplanr`.
* Verify build-time `.glmmplanr_source_md5` across `mirai` worker pools (prepending
  the caller's `glmmplanr` library path on newly started daemons and aborting
  before stamping `simulation_version = 5L` when pre-existing daemons have a
  different `glmmplanr` build loaded), bump `simulation_version` to `5L` so
  `topup_verification()` blocks and `switch_preference()` warns on incompatible
  pre-v5 saved candidate tables (including `simulation_version = 4` plans
  affected by positive-`chance` clamping, non-0.80 `ci_assurance` quantiles in
  `ci_width_q80`, or `n_valid == 0` conditional intervals, alongside pre-v4
  unnormalized contrast weights), and extract global
  bindings from default argument expressions in `formals(fn)` alongside `body(fn)`
  when isolating Box-Muller RNG callbacks in child R sessions.
* Align named random-effect SD vectors (`subj_sd`, `item_sd`, `testlet_sd`,
  `slope_sd_ratio`) and named correlation matrices (`corr_structure`) by design
  random-effect column/term name rather than positional order in `calibrate_dgp()`
  (rejecting partial, duplicate, unknown, or mismatched row/column names), reject
  duplicate, conflicting-alias, or `"(Intercept)"` effect specifications in
  `cell_probs_from_diffs()`, enforce positive whole-integer design dimensions and
  range bounds before integer coercion across `design_crossed()`,
  `design_testlet()`, `valid_multiples()`, `evaluate_design_point()`, and
  `validate_design_point()`, and validate schedule `allocation`, `between_subj`,
  `between_item`, `items_per_subject`, trial counts, unit-ID ranges, per-subject
  item quotas, distinct within-pair conditions (`(subj_id, item_id, cell_id)`),
  and within-unit balance in `simulate_dgp_data()`.
* Initialize above-chance starting intercepts on `[chance + 0.05(1 - chance),
  chance + 0.95(1 - chance)]` with `(1 - chance)`-scaled boundary guards in
  native C++ Laplace estimation, replace fixed `1e-6` / `1e-10` clamping in
  `mafc_binomial()` with IEEE-754 machine-epsilon bounds on `(chance, 1)` (so
  `linkfun(linkinv(eta))` round-trips accurately across both high `chance` near
  `1` and extreme `eta` such as `-18`), and allow `chance` in `[0, 1)` in the
  bundled Stan model;
  return `NA` Monte Carlo intervals and `power_mcse` under
  `failure_policy = "conditional"` when `n_valid == 0L`; keep `ci_width_q80`
  fixed at the 80th percentile alongside `ci_width_q_assurance` at
  `ci_assurance`; plot `"planned"` validation rows with verified zero truth as
  `"Type I error"` against `alpha`; display `rope_ci` in `as_gt()` for
  `criterion = "rope"` plans while rejecting uncomputed `"sensitivity"` or
  `"robust"` table requests; and retain generating DGP, partial-presentation,
  and missingness fields (`dgp`, `items_per_subject`, `item_seed`, `attrition`,
  `participant_dropout`) on returned Rasch planning objects and reports.
* Unify fixed-effect cell and rank identifiability checks across schedule
  generation, the cached C++ simulation path (even without trial attrition),
  and the data-frame path so schedules omitting required conditions fail
  consistently rather than passing regularized width checks. Fixed-item Rasch
  MML (`plan_rasch_fixed_bank(method = "fixed_item")`) marks replications
  invalid when no examinee has two or more observed responses on estimable
  (non-constant) items, and random-item Rasch (`plan_rasch_precision()` and
  `plan_rasch_fixed_bank(method = "random_item")`) requires at least two
  observed participants, two observed items, at least one participant with two
  or more observed responses, and at least one item with two or more observed
  responses so unidentifiable variance components cannot pass tolerance checks.
  `plan_rasch_fixed_bank()` also retains input names on `item_difficulties` and
  includes `item_name` in `item_summary`.
* Reserve all internal schedule and formula identifiers (`subj`, `item`,
  `testlet`, `subj_id`, `item_id`, `testlet_id`, `list_id`, `cell_id`,
  `cell_name`, `trial_id`, `y`, `Z_subj_cells`, `Z_item_cells`,
  `Z_testlet_cells`) against experimental factor name collisions in
  `design_crossed()` and `design_testlet()`, and validate in
  `simulate_dgp_data()` that supplied schedules match `dgp$design` factor
  names, factor level order, cell mappings, and testlet nesting before
  generating responses.
* Check at `specify_target()` (and per scenario in `optimize_sample_size()`)
  that zero lies strictly inside the unit-normalized attainable probability
  contrast for `power`, `both`, and `joint` criteria given `chance` and
  `contrast_weights` (distinguishing mathematical unattainability from
  double-precision calibration limits on extreme weight ratios), and calibrate
  probability-scale nulls and Type M / Type S diagnostics on the same
  normalized scale so weight rescaling never alters null availability or
  diagnostic availability. Share strict cell-name alignment across
  `calibrate_dgp(cell_probs = ...)`, `specify_target()`,
  `cell_contrast_weights()`, `true_prob_effect()`, and multi-scenario specs:
  prefer exact matches against `design$cells$cell_name` (supporting spaces and
  dots in factor levels) before separator-normalized fallbacks, and reject
  partial or duplicate names. Reject all-zero or non-finite weights, compute
  delta, bootstrap, and HMC standard errors on unit-normalized weights before
  scaling back, and bump `simulation_version` to `4L` so `topup_verification()`
  (under both `"independent"` and `"pooled"` policies) refuses to update or
  merge incompatible pre-v4 candidate tables. Track the limiting scenario in
  random-slope `sensitivity` tables across `switch_preference()` and
  `topup_verification()`, display the selected scenario's target effect in
  `as_gt()` headers, compute sensitivity ranges in `report_plan()` across all
  evaluated slope levels (including `"standard"`), recover omitted custom
  `subj_id` / `item_id` values in `validate_custom_data()` by matching
  `schedule$trials` labels rather than re-numbering factor levels, and preserve
  both library paths / MD5 verification and `.GlobalEnv` callback bindings when
  isolating Box-Muller RNG calls in child R sessions. Validate custom `data_fn`
  responses against factor coercion, validate integer replication/worker counts
  before conversion, and allow `true_target` in `evaluate_design_point()` and
  `validate_design_point()` only when a custom `data_fn` is supplied (recording
  built-in null calibration explicitly so tiny non-zero custom `true_target`
  values are never mislabeled as nominal `"Type I error"`).
* Model `valid_rate` explicitly in monotone GP surrogate surfaces under
  `failure_policy = "conditional"` so `min_valid_rate` gates candidate
  acquisition and `pred_criterion_met` (preserving observed `valid_rate` at
  evaluated points and aligning smoothed GP thresholds at `min_valid_rate >= 0.99`,
  while retaining observed feasible candidates when predicted feasibility is
  empty), resolve `stage1_re` and `stage1_engine` per scenario DGP in
  multi-scenario `optimize_sample_size()`, require `n_eval_init >= 1` when
  `method = "surrogate_gp"` explores a strict subset of the post-constraint
  lattice, and solve non-negative monotone GP trend regressions via exact
  centered active-set Cholesky rather than `L-BFGS-B` so 3-point collinear or
  constant surfaces never fail with non-finite `optim` evaluations.
* Return Wilson Monte Carlo intervals for `width_attainment`,
  `joint_attainment`, and `rope_attainment` plus `decision_status` in
  `validate_design_point()`. Preserve per-candidate Wilson levels across
  `topup_verification()` in `report_plan()`, distinguish weighted probability
  contrasts without an attainable zero from absolute probability targets in
  specifications, tables, and reports, note fixed-`item_seed` schedule
  conditioning under partial item presentation, and extract aligned responses
  and cells from `model.frame()` with numeric optimizer-status checks in the
  observed-data `lme4::bootMer()` workflow vignette.
* Restrict interim `early_stop` validity-floor curtailment under
  `failure_policy = "conditional"` to cases where remaining replications cannot
  reach `min_valid_rate`; make `fewer_subj`, `fewer_items` and `fewer_trials`
  strictly lexicographic; report scenario-specific trial counts when schedules
  differ across scenarios; skip duplicate exploratory passes in full-grid
  multi-scenario confirmation; synchronize `evaluated` and `surface` after
  `topup_verification()`; and avoid `0.5` power placeholders on absolute
  probability plots.
* Skip random-slope sensitivity simulations when every slope SD is zero.
  Bonferroni-adjusted decisions now precede step-up, ceiling and top-up checks;
  multi-scenario confirmation uses fresh streams after selecting candidates.
  Subsequent top-ups refresh their replication-count and power metadata.
* Validate custom response IDs and recover omitted testlet IDs from the
  scheduled item-to-passage mapping. Fixed-bank MML rejects duplicate
  person-item responses and constant items incorrectly marked estimable.
* Absolute probability targets report precision and coverage without
  zero-difference power or joint detection claims. Probability-scale power
  criteria require contrast weights of both signs. Fix joint-attainment
  output being joined to the next console bullet.
* Failed observed-curvature checks during polishing or finite-difference
  derivatives now return an invalid-fit status and missing likelihood, rather
  than retaining a working likelihood or throwing a derivative exception.
* Native GLMM fits reject conditional saddles even when a Fisher step is small.
  Final mode checks and derivative evaluations require unclamped observed
  curvature; failed checks cannot return a converged, finite Laplace likelihood.
  Optimization-result printing now includes the same analysis description as
  specification printing and plan reports.
* Native conditional-mode line searches now allow only floating-point roundoff
  when the predicted Newton improvement is below that scale. The original
  convergence tolerance and iteration budgets remain unchanged. Uncomputed
  covariance scores and likelihoods from unconverged fits are returned as `NA`.
  Specification printing and plan reports identify covariance moment-update
  paths, whose convergence does not certify marginal covariance stationarity.
* All planning entry points now use built-in estimators by default. Random-item Rasch planning defaults to `engine = "cpp_laplace"`; use `engine = "lme4"` explicitly to reproduce earlier runs that used the previous default. Fixed-bank MML and the regularized crossed-GLMM defaults are unchanged.
* Declarative crossed, balanced Latin-square, incomplete-block (`items_per_subject`), and testlet-nested (`design_testlet()`) experimental and assessment schedules with separate trial-level missingness (`attrition`) and whole-participant dropout (`participant_dropout`).
* Exact 1D Gauss-Hermite quadrature inversion (`calibrate_dgp()`, `cell_probs_from_diffs()`) over the combined cell normal random-effect variance mapping substantive marginal cell probabilities and $m$-AFC lower guessing asymptotes (`chance = 1/m`) to conditional GLMM parameters, with adaptive quadrature verification against numerical integration.
* Two-dimensional Pareto-optimal sample size, interval-width (AIPE), and ROPE planning (`optimize_sample_size()`, `switch_preference()`) via Universal Kriging monotone Gaussian process surrogates, stabilized C++ Laplace / SQUAREM estimation (`RcppEigen`, Half-Student-t(3) regularization and approximate covariance updates), batch C++ parametric bootstrap, optional Bonferroni multiplicity adjustment (`mc_adjustment`), and independent fresh-batch top-up verification (`topup_policy = "independent"`).
* Scenario comparison plots for verified sample sizes (`plot(plan, type = "scenarios")`) and fixed-design validation plots for detection, coverage, valid-fit counts, and response patterns (`plot(check)` and `plot(check, type = "patterns")`).
* Static plots now use `theme_bw()` and hide automatic titles, subtitles and captions by default. Set `show_annotations = TRUE` to restore them, or add custom `ggplot2::labs()`. This applies to GLMM/Rasch plots and `check_dgp()`; interactive 3D plots use a white background and omit their title by default. Simulation values, intervals and decisions are unchanged.
* Random-item (`plan_rasch_precision()`) and fixed-item-bank (`plan_rasch_fixed_bank()`) Rasch parameter recovery planning with per-item exposure and connected-component diagnostics, extreme response pattern diagnostics (`response_patterns()`), publication-ready visualizations/tables (`plot()`, `as_gt()`), CSV export (`as.data.frame()`), and English sample size justifications (`report_plan()`).
* Analysis-model inspection methods (`formula(spec)`, `family(spec)`, `coef(dgp)`, `check_dgp()`), separation detection for constant fixed-effect cells in unpenalized GLMM fits, and caller-RNG preservation across seeded workflows (including uninitialized sessions and Box-Muller generators).
* Bundled Stan reference model (`inst/stan/crossed_binomial_mafc.stan`) with comprehensive HMC convergence diagnostics and trial-level confirmation outputs (`inst/validation/`).

## Practical item-pool planning

- Added a worked 150-item example separating short-form measurement precision,
  fixed-item recalibration, and participant-by-stimulus experimental planning.
  Candidate content quotas, retained IDs, partial exposure, random slopes,
  and the population to which results apply are explicit. No new public API.
- Added matched irtsim/native Rasch comparisons and equal-budget mlpwr/grid
  planning examples, with shared independent confirmation, retained failures,
  and a surrogate sensitivity check using saved simulation outcomes.
- Added public-data vocabulary-bank targeting examples and published-mean
  replication planning, separating parameter recovery, probability prediction,
  and detection/interval-width goals. Independent checks retain uncertainty
  and cases where planning targets pass but nominal inference is inadequate.
- Expanded these examples to 300 Rasch conditions, 70 paired external checks,
  and 60 learning designs with independent high-replication confirmation.
  Updated pilot interpretations and documented conditional limitations of
  native GLMM interval calibration. Those large runs used the previous
  estimator and are preserved separately from the correction below.
- Corrected zero-correlation, zero-guessing crossed-GLMM covariance updates
  to include movement of the conditional mode in the Laplace determinant.
  Fixed effects and variances now target the same approximate marginal
  objective, and a scaled covariance-score check supplements the stopping
  rule. Iteration-limit diagnostics report that score when applicable.
  Large primary fixed-effect steps are bounded so strong contrasts continue
  to update; numerical failures remain failures even with regularization.
  Screening uses the same numerical settings as confirmation on this path,
  preventing iteration-limit failures unique to screening.
  It retains fewer simulation replications and omits bootstrap refitting.
  This affects both regularized and unpenalized fits, including testlets;
  fixed-bank MML and the correlated/guessing-adjusted paths are unchanged.
  Independent likelihood/gradient checks and same-response confirmation
  assess numerical accuracy, interval calibration, and additional CPU cost.
  No public API was added; planning attainment is still distinct from
  adequate coverage and null rejection.
- Added a separate confirmation of the default probability-scale bootstrap
  (99 inner refits), paired with joint-delta intervals on the same responses.
  Inner failures, Monte Carlo uncertainty, and a 399-refit sensitivity subset
  are retained. Public-workflow checks also cover default full/partial
  planning, independent validation, serialization, and Rasch planning.
- Fixed-bank MML now uses safeguarded EM acceleration after ordinary EM reaches
  its limit and can identify an exact zero-variance solution. The fallback
  checks the variance score as well as parameter updates, avoiding premature
  convergence near zero. Remaining failures continue to count as misses.
- The archived Rasch benchmark refuses to resume batches with incompatible
  settings or without a matching executable/configuration manifest.
