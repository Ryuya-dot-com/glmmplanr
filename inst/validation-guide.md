# Scope, Statistical Validation, and Related R Packages (`glmmplanr`)

<a id="cross-field-designs"></a>

## Scope Across Research Designs and Units of Assignment

In experiments and assessments where multiple participants respond (`0` / `1`) to a shared pool of items, experimental conditions can be assigned within or between participants and within or between items, with either complete or partial item presentation per participant. Both accuracy in psycholinguistic / cognitive tasks and binary preference judgments in behavioral experiments fall within this crossed binary structure. By contrast, when interventions are assigned to intact schools or clinics, the number of clusters and cluster size form distinct hierarchical design axes. Cluster-randomized school trials, multi-rater scoring allocations, and longitudinal within-person time series with serial dependence are outside the scope of `glmmplanr`'s built-in sample-size optimizers; probit models are covered via external-estimator comparisons.

<a id="coverage-recheck"></a>

## Marginal Percentage-Point 95% Interval Coverage

We evaluated 64 participants (32 per group) responding to all 20 items under a true population-averaged group difference of $-5$ percentage points ($60\%$ vs. $55\%$), with random intercept SDs of `0.5` (participants) and `0.3` (items) on the conditional logit scale and 99 parametric bootstrap refits per dataset. Across an initial pilot set of 32 datasets, the built-in engine covered the true contrast in `25/32` datasets; in an independent pre-specified confirmation set of 96 datasets, coverage was `95/96` (`98.96%`, exact binomial 95% CI `[94.33%, 99.97%]`). On the exact same 96 datasets, `lme4::glmer` achieved `94/96` coverage (`[92.68%, 99.75%]`). All 96 intervals converged validly under both estimators and exceeded the pre-specified `90%` lower-bound benchmark.

Under the paired null condition ($60\%$ vs. $60\%$ on the same $64 \times 20$ schedule, 96 datasets per allocation), Type I error rates for the C++ batch parametric bootstrap (`prob_boot = 99`) were `5/96` (`5.2%`) for between-participant assignment and `7/96` (`7.3%`) for within-participant Latin-square assignment. Without refitting, the analytical joint delta method (`prob_se_method = "joint_delta"`, propagating both fixed-effect curvature and variance-component uncertainty) yielded `3/96` (`3.1%`) and `8/96` (`8.3%`) on the identical datasets, agreeing in mean interval width within `0.07` percentage points. Because candidate `pass` status in `optimize_sample_size()` reflects user-specified power or width targets, coverage under alternative sample sizes, missingness rates, or partial-presentation schedules should be verified for the planned design via `validate_design_point()`.

<a id="three-choice-precision"></a>

## Conditional-Logit Uncertainty in 3-Choice (3AFC) Tasks

On a single $48 \times 32$ 3AFC dataset, an early fixed-covariance Laplace approximation returned a conditional fixed-effect SE of `0.207` logit whereas 4-chain RStan HMC returned a posterior SD of `0.291` logit. Decomposing this discrepancy across 40 paired datasets showed that when participant and item covariance matrices were fixed at their data-generating values, C++ Tierney-Kadane Laplace curvature (mean SE `0.201`) and RStan HMC (mean posterior SD `0.204`) agreed closely; the discrepancy arose from (1) omitting covariance-estimation sensitivity and (2) small-sample right-tail integration over half-Student-$t_3$ random-slope priors in HMC (mean posterior slope stretch `+0.104` logit). In `glmmplanr` 0.3.0, `bayes_laplace` combines Tierney-Kadane Laplace curvature with variance-component sensitivity, yielding mean SE `0.227` (`40/40` coverage across the 40 paired datasets; SE `0.237` on the seed-303 benchmark dataset).

When four candidate $(N_{\text{subj}}, N_{\text{item}})$ designs were compared on the same within-participant $2 \times 2$ 3AFC task (true conditional coefficient `1.616` logit) across 40 paired datasets per candidate, the smallest scheduled trial count achieving $\ge 80\%$ attainment of full 95% interval width $\le 1.00$ logit under RStan HMC was `2,304` trials ($72 \times 32$ and $48 \times 48$).

Evaluating the joint-covariance `bayes_laplace` engine on 192 independent confirmation datasets at $72 \times 32$ and $48 \times 48$ yielded `192/192` valid fits and `192/192` width attainment ($\le 1.00$ logit) at both designs, with empirical 95% coverage of `182/192` (`94.8%`, exact 95% CI `[90.6%, 97.5%]`) at $72 \times 32$ and `185/192` (`96.4%`, `[92.6%, 98.5%]`) at $48 \times 48$; paired null Type I error counts were `7/192` and `11/192`, respectively. When RStan HMC is the planned final estimator or when small-sample candidates are considered, verify candidate width and coverage directly under the target estimator (e.g., `engine = "stan_hmc"`).

<a id="rasch-recovery"></a>

## Binary Rasch and Multilevel Parameter Recovery

Using a 1PL Rasch logit model with discrimination `1`, lower asymptote `0`, and population mean person ability `0`, we drew fresh persons and items on each of 32 replications at a fixed total of `1,728` responses. Mean ability RMSE was `0.518` at $144 \times 12$ vs. `0.416` at $72 \times 24$, whereas mean item difficulty RMSE reversed to `0.187` vs. `0.275`. A four-variance-component multilevel illustration (12 schools, 36 classrooms) was also verified against external `lme4` estimation; selecting school-level sample sizes requires dedicated hierarchical schedules ([De Boeck et al., 2011](https://www.jstatsoft.org/article/view/v039i12)).

<a id="rasch-precision"></a>

## Generalizing Across Newly Sampled Examinees and Items

Assuming complete presentation without missingness, `65%` marginal accuracy, person ability SD `0.8`, and item difficulty SD `0.6`, we required simultaneous attainment of ability RMSE $\le 0.50$, difficulty RMSE $\le 0.30$, and item population SD absolute error $\le 0.20$ logit. At $144 \times 24$, all three tolerances were simultaneously satisfied in `63/64` replications and in `62/64` independent confirmation replications (95% Wilson CI `[89.3%, 99.1%]`). Widening the true item difficulty SD to `1.0` reduced joint attainment at $144 \times 24$ to `47/64` (`[61.5%, 82.7%]`).

<a id="rasch-fixed-bank"></a>

## Calibrating a Fixed 16-Item Operational Bank

Fixing 14 core item difficulties evenly between $-1$ and $+1$ logit plus two extreme anchor items at $\pm 3.5$ logit (ability SD `0.8`), 64 replications at $N_{\text{subj}} = 100$ achieved bank-wide difficulty RMSE $\le 0.50$ in `64/64` replications under random-item Empirical Bayes shrinkage and `62/64` under `TAM` fixed-item Marginal Maximum Likelihood (MML), whereas simultaneous worst-item accuracy ($|\hat{b}_i - b_i| \le 0.50$ across all 16 items) dropped to `19/64` and `15/64`. At $N_{\text{subj}} = 30$, the hardest anchor item yielded zero correct responses in `17/64` replications. Under unpenalized fixed-difficulty estimation, all-correct and all-incorrect items are retained in the evaluation denominator as non-attaining replications. `plan_rasch_fixed_bank()` plans examinee sample sizes for a fixed difficulty vector via MML (built-in `"mml"` or `"tam"`) or Empirical Bayes shrinkage (`"cpp_laplace"` or `"lme4"`), whereas `plan_rasch_precision()` resamples items across replications ([TAM MML reference](https://search.r-project.org/CRAN/refmans/TAM/html/tam.mml.html)).

<a id="related-packages"></a>

## Related R Packages by Research Goal

To situate `glmmplanr` 0.3.0 within the R ecosystem, the tables below summarize 19 CRAN and R-universe packages that address statistical power, multi-axis design search, interval precision, or psychometric calibration for binary and multilevel data, based on each package's documented interface and primary methodological reference. Empirical benchmarks in this guide cover marginal probability intervals (`lme4::glmer`), 3AFC conditional-logit intervals (`rstan`), fixed-bank Rasch calibration (`TAM`), and [probability-to-logit input calibration (`PowRPriori`, `simr`, and `mlpwr`)](#package-calibration).

`glmmplanr` is specifically designed for studies that start from substantive population-averaged accuracy differences (percentage points), lower guessing asymptotes, and crossed presentation schedules, and then identify $(N_{\text{subj}}, N_{\text{item}})$ designs meeting power, interval-width, and cost criteria. For arbitrary numbers of hierarchical random factors, group-mean centering of continuous predictors, or free 2PL/3PL slope and guessing estimation, the packages below provide complementary capabilities.

### Power, Precision, and Design Search for Binary Mixed Models

| Package | Primary Inputs & Capabilities | Relationship & Complementarity with `glmmplanr` |
|---|---|---|
| **PowRPriori** | Defines crossed/nested structures via `define_design()`, converts expected binary probabilities into fixed effects, and increments a chosen design axis in `power_sim()`, including within/between predictor centering. | Closest counterpart for canonical binary GLMMs. Standard vignettes apply direct logit transformation, whereas `glmmplanr::calibrate_dgp()` inverts 2D Gauss-Hermite quadrature so marginal probabilities match under the specified random-effect variances ([PowRPriori workflow](https://cran.r-project.org/web/packages/PowRPriori/vignettes/Workflow-Example.html)). |
| **simr** | Evaluates power via `powerSim()` and 1D `powerCurve()` from either fitted `lme4` models or synthetic `makeGlmer()` objects given an external schedule, fixed effects, and variance components. | Standard reference for canonical binomial GLMMs without pilot data when supplied with calibrated conditional parameters and explicit schedule data frames ([Green & MacLeod, 2016](https://doi.org/10.1111/2041-210X.12504); [from-scratch vignette](https://cran.r-project.org/web/packages/simr/vignettes/fromscratch.html)). |
| **mixedpower** | Computes simulation-based power for fixed effects and interactions from a fitted `lme4` model and pilot dataset; `R2power()` scans a second random factor while holding the first at a specified size. | Direct precedent for participant-by-stimulus trade-offs when pilot data and fitted `lme4` models are available ([Kumle et al., 2021](https://doi.org/10.3758/s13428-021-01546-0); [manual](https://lkumle.github.io/power_notebooks/intro/Introduction_to_mixedpower.pdf)). |
| **mlpwr** | Accepts a user-written binary simulation function and cost function in `find.design()` to search multi-parameter design spaces via surrogate modeling, with IRT and multilevel examples. | Surrogate-search counterpart. `glmmplanr::optimize_sample_size(cost = ..., budget = ..., prefer = "min_cost")` adds declarative crossed schedules, $m$-AFC quadrature calibration, and Stage-2 Wilson lower-bound verification to prevent optimistic point-estimate selection ([Zimmer et al., 2024](https://doi.org/10.3758/s13428-023-02269-0); [manual](https://cran.r-project.org/web/packages/mlpwr/refman/mlpwr.html)). |
| **glmmrOptim** | Selects $c$-optimal observation subsets from candidate schedules by minimizing asymptotic linear-contrast variance via the GLMM information matrix. | Useful for optimal condition-level observation allocation; asymptotic variance criteria complement simulation-verified finite-sample CI-width assurance and coverage ([glmmrOptim manual](https://samuel-watson.r-universe.dev/glmmrOptim/doc/manual.html)). |
| **pilotr** | Generates crossed participant-by-item experimental datasets (including Bernoulli outcomes) and evaluates Type S/M errors, mean CI width, and ROPE metrics. | Precision-reporting counterpart; `power_mixed()` in v0.3.1 targets linear mixed models via `lmerTest`/Satterthwaite rather than binary GLMM refitting ([pilotr manual](https://cran.r-project.org/web/packages/pilotr/refman/pilotr.html); [precision & ROPE vignette](https://cran.r-project.org/web/packages/pilotr/vignettes/precision-rope.html)). |
| **simglm** | Simulates configurable regression/multilevel data-generating and fitting models, summarizing power, Type I error, and parameter precision via `compute_statistics()`. | Flexible engine for custom covariate distributions and misspecification stress tests when conditional coefficients are specified directly ([simglm manual](https://cran.r-project.org/web/packages/simglm/refman/simglm.html)). |

### Assignment Units, Guessing Asymptotes, and Design Priors

| Package | Target Design Question | Role Relative to `glmmplanr` |
|---|---|---|
| **designr** | Constructs general crossed, nested, and counterbalanced schedules across arbitrary numbers of fixed and random factors (participants, stimuli, raters). | Reference for multi-factor schedule construction beyond participant-item-testlet hierarchies ([designr reference](https://search.r-project.org/CRAN/refmans/designr/html/designr.html); [overview](https://cran.r-project.org/web/packages/designr/readme/README.html)). |
| **CRTSize** | Analytical sample size formulas (`n4props()`, etc.) for two-arm binary cluster-randomized trials given event rates, cluster size, and ICC. | Appropriate when intact schools/clinics are randomized; single-ICC design effects do not transfer directly to crossed participant-by-item variance structures ([n4props reference](https://search.r-project.org/CRAN/refmans/CRTSize/html/n4props.html)). |
| **quickpsy** | Fits psychometric threshold functions incorporating lower guessing and upper attentional lapse asymptotes. | Shares the lower-asymptote response formulation ($g = 1/3$, lapse $= 0$ under a logistic core) while targeting psychophysical threshold curves rather than crossed GLMM condition contrasts ([Linares & López-Moliner, 2016](https://journal.r-project.org/articles/RJ-2016-008/index.html); [reference](https://search.r-project.org/CRAN/refmans/quickpsy/html/quickpsy.html)). |
| **bayesassurance** | Computes Bayesian assurance (unconditional probability of study success) integrated over continuous design priors on unknown effect sizes. | `glmmplanr` evaluates conditional assurance at each scenario and worst-case multi-scenario robustness (`scenario_dgps()`); prior-weighted integration addresses unconditional assurance ([Pan & Banerjee, 2023](https://journal.r-project.org/articles/RJ-2023-066/)). |

### Rasch, IRT, and Educational Measurement

| Package | Target Measurement Question | Role Relative to `glmmplanr` |
|---|---|---|
| **irtpwr** | Analytical/asymptotic sample size calculation for IRT model-comparison and item-fit hypotheses (e.g., Rasch vs. 2PL). | Suited to hypothesis testing of IRT model structure, whereas `plan_rasch_precision()` bounds parameter RMSE and variance recovery tolerances ([irtpwr vignette](https://cran.r-project.org/web/packages/irtpwr/vignettes/demo.html)). |
| **pwrRasch** | Simulation-based power to detect Rasch model misfit via group-by-item interaction tests. | Addresses DIF / model-fit testing rather than individual item calibration error bounds ([pwrRasch manual](https://cran.r-project.org/web/packages/pwrRasch/refman/pwrRasch.html)). |
| **TAM** & **eRm** | Item calibration for operational item banks via Marginal Maximum Likelihood (`TAM`) or Conditional Maximum Likelihood (`eRm`). | `plan_rasch_fixed_bank()` cross-validates directly against `TAM::tam.mml()` on identical response matrices; `eRm` provides CML sensitivity checks ([TAM reference](https://search.r-project.org/CRAN/refmans/TAM/html/tam.mml.html); [Mair & Hatzinger, 2007](https://www.jstatsoft.org/article/view/v020i09)). |
| **mirt** & **eirm** | `mirt` simulates and estimates multidimensional 1PL–4PL IRT models; `eirm` fits explanatory IRT models linking item/person covariates to responses via `lme4`. | `mirt::simdata()` supports 2PL/3PL misspecification studies; `eirm` complements `glmmplanr`'s explanatory IRT sample-size planning at analysis time ([mirt::simdata](https://philchalmers.github.io/mirt/reference/simdata.html); [Bulut et al., 2021](https://doi.org/10.3390/psych3030023)). |
| **optical** | Optimal item-to-examinee block assignment using preliminary item parameters to maximize calibration information. | Useful for adaptive/targeted incomplete-block routing, whereas `items_per_subject` in `glmmplanr` plans uniform random incomplete-block sampling ([optical CRAN page](https://cran.r-project.org/web/packages/optical/index.html)). |
| **gtheory** | Generalizability (G) and Decision (D) studies quantifying variance components and composite score reliability across items, raters, and occasions. | Preferred when planning rater counts for Classical Test / G-theory dependability coefficients rather than logit Rasch recovery ([gstudy](https://search.r-project.org/CRAN/refmans/gtheory/html/gstudy.html); [dstudy](https://search.r-project.org/CRAN/refmans/gtheory/html/dstudy.html)). |

### Alignment Checklist When Comparing Sample-Size Tools

1. **Marginal vs. Conditional Input Probabilities.** Under a canonical logit GLMM with a $\mathcal{N}(0, 1^2)$ random intercept, setting $\beta_0 = \text{qlogis}(0.80) = 1.3863$ yields $80\%$ accuracy at zero random effect (`conditional_prob`) but shrinks the population-averaged accuracy (`quad_prob`) to $76.13\%$; achieving an $80\%$ population-averaged rate requires $\beta_0 = 1.6497$. Always verify simulated cell means when passing nominal percentages across packages.
2. **Success Criteria and Failure Denominators.** Hypothesis-test power, mean CI width, CI-width assurance ($\Pr(\text{width} \le w)$), and empirical interval coverage answer distinct questions ([Pargent et al., 2024](https://doi.org/10.1177/25152459241287132)). Align link functions, contrast coding, interval estimators, and whether non-convergent replications remain in the denominator ([see 3AFC comparison above](#three-choice-precision)).
3. **Target Item Population.** Calibrating a fixed 16-word list keeps item difficulties fixed across replications (`plan_rasch_fixed_bank()`), whereas generalizing to new lexical samples redraws item difficulties on every replication (`plan_rasch_precision()`; [Dunn, 2024](https://www.sciencedirect.com/science/article/pii/S2772766124000491)). Report both bank-average RMSE and all-item worst-case attainment.

<a id="package-calibration"></a>

## Passing Calibrated Marginal Probabilities to `PowRPriori` and `simr`

Suppose we plan for population-averaged cell accuracies of $72\%$ vs. $80\%$ with participant and item random-intercept SDs of `0.8` and `0.6` logit (`chance = 0`, `treatment` coding). In `check_dgp()`, `quad_prob` is the population-averaged accuracy across new participants and items, whereas `conditional_prob` is the accuracy when all random effects equal zero (`75.65%` vs. `83.89%` after quadrature calibration):

```r
library(glmmplanr)
design <- design_crossed(list(group = c("A", "B")), allocation = "between_subj",
                         contrast = "treatment")
dgp <- calibrate_dgp(design, cell_probs = c(A = 0.72, B = 0.80),
                     subj_sd = 0.8, item_sd = 0.6, slope_sd_ratio = 0, corr_structure = 0)
cells <- check_dgp(dgp, plot = FALSE)$cell_diagnostics
cells[c("cell", "target_prob", "quad_prob", "conditional_prob")]

b <- PowRPriori::fixed_effects_from_average_outcome(
  y ~ group, data.frame(group = factor(cells$cell, levels = c("A", "B")),
                        y = cells$conditional_prob),
  family = "binomial", center = FALSE
)
```

Passing `cells$conditional_prob` yields intercept `1.13345` and group contrast `0.51625`, which can be supplied to `simr::makeGlmer(fixef = unlist(b), ...)` alongside random-intercept variances `0.64` and `0.36`. By contrast, transforming `0.72` and `0.80` directly via `qlogis()` yields slope `0.44183` and shrinks the true population-averaged cell probabilities to `68.68%` vs. `76.13%` (a `7.45 pp` rather than `8.00 pp` marginal difference; see [PowRPriori workflow](https://cran.r-project.org/web/packages/PowRPriori/vignettes/Workflow-Example.html) and [simr from-scratch vignette](https://cran.r-project.org/web/packages/simr/vignettes/fromscratch.html)).

Across 192 datasets on a complete $64 \times 24$ schedule (`1,536` observations, two-sided 5% Wald test), calibrating the 8-percentage-point marginal difference produced `188/192` (`97.9%`) detections for both `glmmplanr` and `simr` under a within-participant Latin square, and `115/192` (`59.9%`, `glmmplanr`), `111/192` (`57.8%`, `simr` on identical responses), and `112/192` (`58.3%`, independent `PowRPriori::power_sim()`) under between-participant assignment (`32` vs. `32`). Without quadrature calibration, between-participant detections dropped to `95`, `92`, and `96` of `192`, illustrating how conflating marginal and conditional probabilities alters both the simulated population and resulting power estimates.

On a $5 \times 5$ within-participant integer lattice ($68\%$ vs. $80\%$, total search budget `1,200` replications) under asymmetric cost functions ($10 N_{\text{subj}} + 80 N_{\text{item}}$ vs. $80 N_{\text{subj}} + 10 N_{\text{item}}$), `glmmplanr::optimize_sample_size(method = "surrogate_gp", cost = ..., prefer = "min_cost")` selected the exact exhaustive-grid minimum-cost passing designs ($(40, 12)$ at cost `1,360` and $(16, 24)$ at cost `1,520`), achieving `90.3%` `[87.0%, 92.8%]` and `84.0%` `[80.1%, 87.3%]` power in 400-replication out-of-sample confirmation, whereas targeting the 80% point-estimate contour directly without a Stage-2 lower-confidence gate selected $(34, 10)$ with `73.8%` `[69.2%, 77.8%]` out-of-sample power.
