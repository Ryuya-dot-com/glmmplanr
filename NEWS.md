# glmmplanr 0.1.0

* Initial CRAN release.
* Declarative two-dimensional experimental and assessment design specifications
  for crossed participant-by-item (`design_crossed()`) and passage/testlet-nested
  (`design_testlet()`) binary response studies, including Latin-square,
  between-within, and incomplete-block (`items_per_subject`) schedules with
  separate trial-level attrition and participant dropout.
* Exact 1D Gauss-Hermite quadrature calibration (`calibrate_dgp()`,
  `cell_probs_from_diffs()`, `check_dgp()`) mapping population-averaged
  condition probabilities and $m$-AFC guessing asymptotes (`chance`) to
  conditional GLMM parameters.
* Target and precision specification (`specify_target()`, `target_control()`)
  on percentage-point (`"pp"`), probability (`"prob"`), or conditional logit
  (`"logit"`) scales across power, confidence-interval width (AIPE),
  combined/joint detection-and-width, and ROPE criteria.
* Pareto-optimal `(N_subj x N_item)` sample size planning
  (`optimize_sample_size()`, `switch_preference()`, `topup_verification()`,
  `scenario_dgps()`) using Universal Kriging monotone Gaussian process
  surrogates, compiled `RcppEigen` Laplace estimation (`bayes_laplace`,
  `cpp_laplace`), C++ batch parametric bootstrap standard errors, Wilson Monte
  Carlo uncertainty intervals, and optional `mirai` parallel workers.
* Fixed-design verification and diagnostics (`evaluate_design_point()`,
  `validate_design_point()`, `response_patterns()`) plus optional external
  backend cross-checks (`lme4`, `glmmTMB`, and Stan HMC via `cmdstanr` or
  `rstan`).
* Random-item (`plan_rasch_precision()`) and fixed-item-bank
  (`plan_rasch_fixed_bank()`) Rasch parameter recovery planning, along with
  `ggplot2` visualizations (`plot()`), `gt` summary tables (`as_gt()`),
  data-frame export (`as.data.frame()`), and English methods summaries
  (`report_plan()`).
