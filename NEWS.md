# glmmplanr 0.1.0

* Initial public release of `glmmplanr`.
* Declarative crossed, balanced Latin-square, incomplete-block (`items_per_subject`), and testlet-nested (`design_testlet()`) experimental and assessment schedules with separate trial-level missingness (`attrition`) and whole-participant dropout (`participant_dropout`).
* Two-dimensional Gauss-Hermite quadrature inversion (`calibrate_dgp()`, `cell_probs_from_diffs()`) mapping substantive marginal cell probabilities and $m$-AFC lower guessing asymptotes (`chance = 1/m`) to conditional GLMM parameters.
* Two-dimensional Pareto-optimal sample size, interval-width (AIPE), and ROPE planning (`optimize_sample_size()`, `switch_preference()`) via Universal Kriging monotone Gaussian process surrogates, stabilized C++ Laplace / SQUAREM estimation (`RcppEigen`), batch C++ parametric bootstrap, and Wilson Monte Carlo confidence bounds.
* Random-item (`plan_rasch_precision()`) and fixed-item-bank (`plan_rasch_fixed_bank()`) Rasch parameter recovery planning, extreme response pattern diagnostics (`response_patterns()`), publication-ready visualizations/tables (`plot()`, `as_gt()`), and Registered Report text generation (`report_plan()`).
