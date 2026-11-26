# glmmplanr 0.1.0

* Initial public release of `glmmplanr`.
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
