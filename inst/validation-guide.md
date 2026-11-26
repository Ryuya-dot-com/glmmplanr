# Scope, Statistical Validation, and Related R Packages (`glmmplanr`)

## Built-in Estimation and External Checks

The default planning engines are `bayes_laplace` for crossed binary GLMM
effects, `cpp_laplace` for random-item Rasch recovery, and `mml` for fixed-item
Rasch banks. Fixed-bank Rasch has fixed item intercepts and a random person
intercept; its one-dimensional quadrature MML solver is a separate numerical
path from crossed-GLMM Laplace estimation. A correction or validation of one
path does not establish the accuracy of the other.

TAM and mirt are the primary external references for fixed-bank IRT checks.
Match the response matrix, latent scale, discrimination and variance
constraints, and failure denominators. Crossed random-item effects require
a matching mixed model for external comparison. Previous binary special-case
comparisons with mfrmr remain supporting evidence; its main application is
many-facet Rasch measurement. Existing results keep their recorded engines.
To reproduce earlier random-item planning with the previous default, specify
`engine = "lme4"` explicitly.

The [native covariance-update correction](#l2-uncertainty-correction) changes
the zero-correlation, zero-guessing crossed-GLMM path. Earlier GLMM benchmarks
below retain their original implementation and should not be read as a fresh
validation of every old design under the correction. Fixed-bank MML is
unchanged by this update and retains its existing calibration evidence.

**Current numerical diagnostics.** Native covariance-score stopping is
implemented for zero-guessing, uncorrelated fits. Other paths retain fast
approximate moment updates; their convergence does not certify stationarity
of the same marginal covariance objective. Specifications and plan reports
identify those paths. An uncomputed covariance score is `NA`, including fits
with covariance held fixed. Native likelihoods are also `NA` unless the fit
converges and the final observed-curvature calculation is ready. Do not use
intermediate values from failed fits as objective evaluations.

The conditional-mode line search now accounts for floating-point roundoff
only when a small Newton step has predicted improvement below that scale.
It then rechecks the original 2e-7 mode tolerance; neither this tolerance nor
the iteration budget is relaxed. The study below used the preceding source;
its raw records retain the old zero diagnostic sentinel and failed returns.
The same-response repair replay is recorded separately in the validation
suite. All 69 previously failed diagnostic optimizations recovered under the
same acceptance rules, and 15 selected valid controls remained valid. Native
replays of all 2,048 responses changed no coverage, detection or width decision.
Three fits exceeded the predeclared 1e-4pp numerical comparison tolerance:
maximum point/SE changes were .000222/.000381pp. Independent likelihood and
mode checks passed for those exceptions; they are not reported as bitwise
identical fits. This change does not turn the moment updates into exact joint
MAP estimation or introduce a new public estimator. The selected repair
replay is not an independent confirmation of nominal inference.

<a id="covariance-review-medium"></a>

## Focused Check of Covariance Paths

A paired study used 16 conditions with 128 datasets each: complete 32-person
by 24-item presentation or partial presentation to 72 people over 60 items,
with 30 per person,
guessing floor 0 or 1/3, random-effect correlation 0 or .2, and marginal
contrasts 0 or 10pp from baseline .55. Participant intercept/slope SDs were
.8/.35; item SDs .6/.4. Fixed Latin-square schedules used seed 37. This is
a limited DGP comparison, not an estimate of an optimal required sample size.

All 2,048 native fits were valid. A diagnostic optimizer of the same Laplace
objective, with full profile covariance propagation, supplied 1,979 valid
intervals. It shares the native likelihood evaluator and is not an independent
reference or a public engine. Numerical checks found small log-SD scores
on the corrected zero-guessing/zero-correlation path; all 1,521 evaluable
scores on the other paths exceeded the diagnostic .01 threshold. Fifteen
additional scores could not be evaluated because of failed perturbations.
In that historical run, a zero returned covariance score on those native
paths denoted a bypassed check. The current implementation returns `NA`.

Among common valid pairs, noncanonical-path effect changes had cell-specific
RMS .032–.120pp, maximum .538pp, and RMS-SE ratios .985–1.004. Eight of 1,979
pairs changed joint detection/width classification. Native coverage ranged
93.75–98.44%; null rejection ranged 3.13–6.25%, with all eight pointwise MC
intervals including .05. At 128 datasets, MCSE near 5% is 1.93pp, so these
results cannot establish exact nominal calibration. Failed reference fits
remain unmet in all-trial detection/attainment rates; valid-only comparisons
do not remove their limitations.

For a 10pp effect and width <=20pp, the complete layout had three failed
planning checks and one indeterminate check; the larger partial layout
passed all four. The diagnostic procedure gave the same categories despite
its failures. Counts, sample sizes and total exposure change together, so
this does not demonstrate that partial presentation itself improves precision.
The study did not test the default bootstrap intervals or justify switching
defaults. Numerical stationarity, successful fitting and inferential accuracy
remain distinct requirements. Detailed scripts, responses, failures and
MC intervals are retained in the repository's validation suite.

<a id="visualization"></a>

## Reading and Customizing Simulation Figures

Static plots use `ggplot2::theme_bw()` with automatic titles, subtitles and
captions hidden by default. This styling choice follows the package user's
preference; the papers below do not prescribe this ggplot theme. Axis names,
legends, scenario panels, interval bars and diagnostic counts remain visible.
`show_annotations = TRUE` restores the generated explanations in GLMM/Rasch
`plot()` methods and `check_dgp()`. Add `ggplot2::labs()` for custom text.
The optional 3D display uses a white background and no default title;
customize it with `plotly::layout()`.

[Morris, White and Crowther (2019)](https://doi.org/10.1002/sim.8086)
organize simulation reporting by estimand, generating mechanism, method and
performance measure, and show performance estimates with Monte Carlo
uncertainty. Accordingly, our validation and scenario figures retain their
stored Monte Carlo intervals and separate rates from interval widths.
These bars quantify simulation uncertainty in a performance estimate;
they are not the analysis intervals whose coverage is being evaluated.

[Rücker and Schwarzer (2014)](https://doi.org/10.1186/1471-2288-14-129)
use nested-loop and Trellis displays to organize many simulation conditions.
Their work informs our emphasis on explicit scenario labels and panels.
The current API retains small-multiple comparisons of selected designs;
it does not implement their nested-loop method or imply that a few selected
designs represent an entire parameter grid.

[Gasparini's rsimsum (2018)](https://doi.org/10.21105/joss.00739) and its
[visualization documentation](https://ellessenne.github.io/rsimsum/articles/C-plotting.html)
provide related examples of performance plots grouped by method and
generating conditions. These are presentation references, not new estimator
benchmarks or dependencies. There is no new simulation performance claim.

When exporting an unannotated figure, describe the DGP, sample/exposure
schedule, replication count, interval level and failure denominator in the
surrounding text or manuscript caption. Coverage and mean widths use valid
fits; detection denominators follow the selected failure policy. Mean-width
panels do not gain an uncertainty interval when none was stored. Plotting
does not change simulated outcomes, denominators or pass/indeterminate/fail
decisions. The selected marker on a contour or Pareto plot can represent an
indeterminate or failed candidate: consult `summary()`/`report_plan()` or
enable annotations for its status before describing it as recommended.

<a id="item-pool-planning"></a>

## Practical Pilot: Short Forms and Stimulus Pools

A synthetic 150-item Rasch pool has three content strata of 50 items, each
spanning -2.5 to 2.5 logits. Candidate forms retain equal stratum counts:
30 central-difficulty items, 30 or 60 spread-difficulty items, or all 150.
Their identities are fixed before response generation. This is a worked
planning example, not validation of an actual short form.

**Measurement with known items.** Under a fixed N(0,1) EAP scoring prior,
integration over every possible score gives the following conditional person
RMSE. Item parameters are treated as known; item estimation uncertainty is
not included. Halving the numerical grid spacing changed bias/RMSE by less
than 0.000001 logits.

| Form | RMSE at ability 0 | RMSE at ability 3 |
|---|---:|---:|
| Central 30 | 0.341 | 0.907 |
| Spread 30 | 0.385 | 0.805 |
| Spread 60 | 0.292 | 0.534 |
| Full 150 | 0.193 | 0.316 |

Which form is acceptable depends on the ability range and allowable error.
EAP bias and information are also saved: neither information alone nor a
single population-average RMSE describes conditional scoring accuracy.

**Recalibration sample size.** Native fixed-item MML evaluated 30 scenarios:
five layouts (the four forms plus 150 items with 30 per person), N=120/240/480,
and person SD=0.4/0.8, each with 128 replications. The illustrative joint target
is bank RMSE <= 0.35 and every-item absolute error <= 0.60 logits; a pass
requires the pointwise Wilson 95% lower bound to reach assurance 0.80.
All 3,840 fits completed, but some banks had constant items and were
unevaluable for whole-bank accuracy; these remain misses.

At SD=0.8 and N=120, the central 30 passed in 113/128 runs (88.3%, Wilson
95% interval 81.6–92.8%); the spread 30 passed in 53/128 (41.4%, 33.2–50.1%).
The same central form had worse measurement at ability 3. Ease of calibration
therefore cannot substitute for short-form measurement evidence.
At N=480, all four fixed forms passed the assurance criterion for both SDs;
the rotating 150-item pool met the joint tolerance only 1/128 times at each
SD. Its mean exposure was 96 per item, compared with 480 for a fixed form.
Requiring all 150 items to meet the maximum-error criterion also makes the
joint target more demanding than for 30 items.

**Experimental stimuli.** A separate native crossed-GLMM pilot used N=60/120,
J=30/60/150 and 30 items per person in a two-condition Latin square. The
marginal planned difference was 8 percentage points (50% versus 58%), with
person/item intercept SDs 0.8/0.6. Slope SDs were either zero or 0.25/0.35,
with zero correlations and a matching analysis. Each of 12 layouts/profiles
had 128 planned and 128 null runs; all 3,072 fits were valid. For N=60 and
heterogeneous slopes, detection was 82.8%, 89.1%, and 85.2% at J=30,60,150,
with overlapping Wilson intervals. This does not establish an optimal J.
Across the pilot, null detection ranged from 1.6% to 8.6%; 128 replications
do not tightly establish Type I error control. New items were sampled each
replication: these results do not estimate the effect of selecting particular
word identities or changing the target population to high-frequency words.

This example uses existing APIs. Content-constrained assembly can supply
candidate forms, for example through [TestDesign](https://choi-phd.github.io/TestDesign/).
glmmplanr evaluates their design consequences. The matched extension below
now includes `mlpwr` and `irtsim`, with shared objectives, candidate designs,
failure denominators and independent confirmation. The two comparisons
address estimation and design search separately.

Scripts, parameters, item IDs, full curves, trial-level plan objects and
Monte Carlo intervals are recorded in the repository's
[item-pool validation](https://github.com/Ryuya-dot-com/glmmplanr/tree/main/validation#item-pool-planning).
The macOS pilot took about 11 seconds for recalibration and 80 seconds for
the experimental planned/null checks, excluding installation and compilation;
these are observed times for these settings, not runtime guarantees.

<a id="item-pool-package-comparison"></a>

### Matched Comparisons with irtsim and mlpwr

The extension asks whether another Rasch workflow supports the same
recalibration plans, and whether another search algorithm identifies useful
designs for the same native estimator. It uses the synthetic pool and content
quotas above, without adding a public API.

**irtsim 0.2.0.** Its public 1PL workflow with mirt and native MML evaluates
identical response matrices, item IDs, masks, discrimination 1 and a mean-zero,
free-variance ability scale. Five layouts, two ability SDs and five sample
sizes (120/240/480/960/1920) give 50 conditions ×128 replications.
Both estimators have 6,095/6,400 complete valid banks. Only two per-trial
joint decisions differ, one in each direction at the .60 maximum-error
boundary. All 50 attainment counts and ten proposed sample sizes agree.
Item estimates are not identical (maximum paired difference .02620 logits).
The [screening figure](https://github.com/Ryuya-dot-com/glmmplanr/blob/f35e579103825f9e8c4b1d49f3f0a6772f1e5834/validation/item_pool_comparison.png)
shows the layout effect; it is not evidence of equivalence outside these
Rasch conditions.

Independent confirmation used 15 distinct returned points ×512 fresh
replications with both IRT engines. All 15 decision classifications agree
(eight per-trial differences among 7,680 pairs). Of the ten minimum-N screen
proposals, nine pass and one remains indeterminate: spread-60, SD=.4,
N=240 reaches 80.7% (95% interval 77.0–83.9%). A separately proposed N=480
passes at 98.8% (97.5–99.5%). Rotating 150 items with 30 per person requires
the tested N=1920 candidate to pass: 93.8% (91.3–95.5%) at SD=.4 and
95.1% (92.9–96.7%) at SD=.8. These are tested candidates under specified
whole-bank tolerances, not universal minimum sample sizes.

**mlpwr 1.1.1.** It receives the same native success function, candidate grid
and cost as glmmplanr's uniform grid evaluation. This includes Rasch recovery
as well as GLMM detection/precision; mlpwr already supports custom IRT
simulations, as illustrated in its [IRT tutorial](https://CRAN.R-project.org/package=mlpwr/vignettes/IRT_Vignette.html).
This comparison does not evaluate glmmplanr's GP optimizer. There are 50
Rasch searches and 20 GLMM searches per method, using five search seeds per
family/preference. Actual simulation budgets match at 30,000 Rasch and
20,880 GLMM evaluations per method. Errors, warnings and unqualified
fallback coordinates remain in the report.

For the GLMM search, both methods selected independently confirmed designs
in 19/20 runs. The grid's other candidate (zero slopes, 30×30) is
indeterminate; mlpwr's other run encountered a reproducible duplicate-index
error after 360 evaluations. With heterogeneous slopes, 60×30 attained the
joint detection/20pp-width target in 428/512 fresh trials (83.6%, Wilson
95% interval 80.1–86.5%). Its null rejection was 29/512 (5.7%, 4.0–8.0%);
passing the planning target does not guarantee nominal Type I error control.
All confirmation intervals met the 20pp width limit, so detection drives the
joint attainment differences in this example.

In ten rotating-pool searches, mlpwr's default logistic surrogate warned
that no good design was found. Replacing only that surrogate with a GP,
using the same saved 600 outcomes and no additional simulations, returned
N=1920 in all ten. This is an exploratory refit after a default-setting
failure, not a prespecified GP search benchmark. It illustrates why a
particular surrogate's failure is not a general limitation of mlpwr. Both
returned N=1920 designs pass the independent confirmation.

The [full protocol, results and reproducible scripts](https://github.com/Ryuya-dot-com/glmmplanr/tree/main/validation#item-pool-package-comparison)
record seeds, common data checks, failures, actual budgets and Monte Carlo
intervals. Repeated searches returning the same candidate share its
independent confirmation and must not be counted as independent precision
replications.

<a id="l2-empirical-planning"></a>

### Empirical Vocabulary Banks and Replication Planning

An additional pilot uses [Hamada et al.'s public vocabulary data](https://osf.io/t8ax6/):
1,220 participants, 50 word identities, and five booklets with 20 recall and
30 recognition responses per person. Separate 2PL fits follow the authors'
model calls. Separate Rasch working-model refits supply the native simulation
banks; they are not the source paper's reported model. Their native/mirt item
differences are below .001 logits. The two formats retain their own scales.

The native pilot has 240 conditions ×64 replications, crossing two formats,
five 25/50-item layouts (including the original booklets), six ability
profiles and N=100/200/500/1000. A 24-condition subset has 32 paired mirt Rasch
and 32 separate mirt 2PL runs. Ability shifts are sensitivity assumptions,
not age or CEFR labels. All failures and constant-item banks remain misses.

For the recall Rasch bank with all 50 items and N=1000, difficulty attainment
(RMSE <=.35 and maximum error <=.60) falls from 60/64 (93.8%, Wilson
95% interval 85.0–97.5%) in the reference distribution to 4/64
(6.3%, 2.5–15.0%) when its mean is shifted -3 logits, and 0/64
(0–5.7%) at -5 logits. The latter has 4.73 constant items on average.
For the easier recognition bank, the -3 shift improves attainment from
24/64 (37.5%) to 61/64 (95.3%) by reducing ceiling responses. Thus the
relevant quantity is the population/item match, not simply low ability.
The [targeting figure](https://github.com/Ryuya-dot-com/glmmplanr/blob/f35e579103825f9e8c4b1d49f3f0a6772f1e5834/validation/l2_targeting.png)
shows these within-format comparisons.

In the separate correctly specified 2PL recall arm at N=1000, response-curve
RMSE <=.05 is attained in 32/32 trials, but the whole-bank difficulty and
discrimination tolerances jointly in only 1/32. These outcomes concern
different uses of the estimates. This pilot does not add native 2PL fitting,
validate person-score precision, propagate calibration-parameter uncertainty,
or establish a universal required sample size.

A second example makes the existing Huang et al. vocabulary-learning vignette
concrete. It plans a focused three-round computer-test comparison from
published means .42/.50. Missing random-effect variances are explicit
assumptions; the original full interaction design is not reproduced.
The target is simultaneous detection and a 95% interval width <=10pp,
with 30 responses per person. Equal-budget grid and mlpwr GP proposals are
independently checked with 512 planned and 512 null replications per point.

| Assumed profile | Confirmed candidate N × pool J | Joint attainment (Wilson 95%) |
| --- | --- | --- |
| 8pp, intercept SDs .8/.6, no slopes | 54 × 30 | 94.7% (92.4–96.4%) |
| 8pp, same intercepts, slope SDs .25/.35 | 72 × 60 | 84.8% (81.4–87.6%) |
| Same slope profile, fewer items prioritized | 108 × 30 | 91.0% (88.2–93.2%) |
| 5pp, slope SDs .25/.35 | 162 × 60 | 84.6% (81.2–87.4%) |
| 8pp, intercept SDs 1.1/.9, slopes .4/.5 | 162 × 30 | 94.9% (92.7–96.5%) |

Passing the planning target is not the same as validated nominal inference.
For example, another high-variance candidate, 72×90, passes at 90.6% but has
null rejection 7.2% (5.3–9.8%) and planned coverage 92.0% (89.3–94.0%). It
is therefore not an unqualified recommendation. The 162×30 alternative has
null rejection 4.9% (3.3–7.1%) and coverage 94.1% (91.8–95.9%). The 10pp
width condition also affects decisions: moderate slopes at 72×60 give
92.4% detection, 91.8% width attainment and only 84.8% joint attainment.

These are the initial 512-run results. In the subsequent independent
2,048-run confirmation, the 72×90 warning was not reproduced (null 5.0%,
coverage 94.7%), while 162×30 gave null rejection 7.1% (6.0–8.3%). The
[large-run evidence](#l2-large-confirmation) updates the candidate assessment;
the pilot's apparent ranking of those two designs should not be retained.

[Full methods, parameter provenance, results and scripts](https://github.com/Ryuya-dot-com/glmmplanr/tree/main/validation#l2-empirical-planning)
separate the 64/32-replication targeting pilot from the 512-replication design
confirmations. These are tested grid candidates under stated assumptions,
not universal minima or exact reconstructions of the source studies.

<a id="l2-large-confirmation"></a>

### Large Independent Vocabulary Confirmation

The public-data targeting study now has 300 native Rasch conditions with
1,024 fresh replications each, including N=2000. Seventy representative
conditions also have 256 paired mirt Rasch and 256 separate mirt 2PL trials.
All 17,920 native/mirt Rasch response hashes and validity decisions match;
five trial-level precision decisions differ, but all 70 condition-level
classifications agree at the same 256-trial denominator.

At N=2000 with all 50 items, native Rasch difficulty attainment is 99.9%
for reference recall, 49.2% after shifting the mean -3, and 0% after -5.
For reference recognition it is 65.0%, rising to 99.4% after the -3 shift
reduces ceiling responses. These are each format's own fitted scale, not
age or CEFR differences. At N=1000, reference recall gives 96.1% under full
presentation, 64.5% under equal booklets and 56.8% under the original booklet
proportions. Both booklet allocations pass at N=2000. Selection of fewer
or harder items changes the calibration problem and does not establish
equivalent vocabulary-size measurement.

For the correctly specified 2PL recall bank at N=2000, all 256 runs satisfy
the probability-curve criterion but only 54/256 (21.1%, 16.5–26.5%) satisfy
the joint difficulty/discrimination criterion. For recognition these counts
are 253/256 and 0/256. Whole-bank parameter tolerances, in-population
probability prediction and person-score measurement are distinct purposes.
For recall, difficulty RMSE alone succeeds in 249/256 trials: the maximum
item-error limits account for much of the stricter criterion's low success.
The probability RMSE is averaged across items and ability points and does
not guarantee a small error for every individual item curve.
The examples do not implement native 2PL or replace conditional planning
with a universal required N.

Before the covariance correction, the learning study confirmed all 60 designs with 2,048 planned
and 2,048 null runs each; all 245,760 fits are valid. Representative joint
attainment results are 94.5% (93.4–95.4%) for intercept-only 54×30, 84.5%
(82.8–86.0%) for moderate-slope 72×60, and 93.1% (91.9–94.1%) for moderate-slope
108×30. Every person still answers 30 items. Under the 5pp effect, 162×60
passes planning at 85.7% but has coverage 93.9% (92.8–94.9%); 216×60 has
91.9% joint attainment and 95.2% (94.2–96.1%) coverage in that profile.

For larger variances, the old 72×90 warning is not reproduced: joint
attainment is 91.2% (89.9–92.3%), null rejection 5.0% (4.1–6.0%) and
coverage 94.7% (93.7–95.6%). Conversely, 162×30 has null rejection 7.1%
(6.0–8.3%) despite 94.1% joint attainment. These updates preserve the old
512-run results while giving precedence to the fresh, larger confirmation.

Thirty-nine of 60 designs pass planning, but 14 of those also have an
adverse nominal-departure flag; 18 designs are flagged overall. These use
pointwise intervals, not familywise tests. Four layouts pass planning in
all four generating profiles, but each has at least one inference flag
across those profiles. This grid therefore does not establish a universally
adequate plan across the tested variance/effect assumptions. At high-variance
216×90, for example, RMS estimated SE is 1.480pp versus empirical SD 1.548pp,
with bias +0.176pp and 93.6% coverage. The subsequent
[mechanism checks](#l2-uncertainty-correction) identified an inconsistency in
the covariance update; residual interval approximation is assessed separately.

The equal-budget search comparison now uses five seeds per profile/preference:
40 searches and 43,200 evaluations per method. Uniform grid proposals give
30 confirmations, six indeterminate and four failed decisions; mlpwr GP gives
34 confirmations and six indeterminate decisions. Neither search has an
error; 29 mlpwr boundary warnings are retained. Both use the same native
estimator and shared independent confirmation, so this is a limited search
comparison and does not certify the uncertainty approximation.

[Protocol, full counts, intervals, failures and reproducible scripts](https://github.com/Ryuya-dot-com/glmmplanr/tree/main/validation#l2-large-confirmation)
preserve the smaller pilot separately. Parameters fitted to the source data
are fixed here; their calibration-sample uncertainty is not propagated.

<a id="l2-uncertainty-correction"></a>

### Native Covariance Updates and Interval Calibration

The learning example exposed a practical distinction: a design can achieve
the requested detection and width targets while its nominal interval
under-covers. Saved trials were used to investigate this before expanding
the design grid. Fixing variance components at their generating values was
a diagnostic comparison only, not an estimator available in practice.

For zero assumed random-effect correlations and no guessing lower asymptote,
the covariance update now includes movement of the conditional random-effect
mode in the Laplace determinant. Fixed effects and variances target the same
approximate marginal objective. A small parameter change alone is no longer
enough to declare convergence; the scaled covariance score must also be small.
An iteration-limit failure remains a failure in planning denominators, even
when finite coefficients can be returned. The correction applies to native
regularized and unpenalized crossed fits, including testlets. It does not
change fixed-bank MML or the correlated/guessing-adjusted calculation paths.

Screening uses the same numerical settings as confirmation on the corrected
path. Its former smaller iteration budget could reject a dataset that fitted
successfully in confirmation and thereby underestimate a design's attainment.
Screening remains cheaper through fewer simulation replications and omission
of bootstrap refitting; a smaller numerical budget is not used to obtain
that saving on this path.

Strong contrasts are also checked against a separate calculation of the
fixed-effect score. Large coefficient updates continue in bounded steps.
For example, a simulated sample with about 10% versus 91% correct gives an
estimated marginal difference of 80.9pp. A fully separated stress example
still reaches the iteration limit and remains invalid for planning. A prior
permits finite coefficient estimates without guaranteeing convergence within
the numerical budget; finite coefficients alone do not establish a valid fit.

Independent profile optimization, a separate conditional Hessian calculation,
and finite-difference scores check the numerical change. They establish
agreement with the same approximate objective, not exact integration or
nominal frequentist coverage. The follow-up compares six selected layouts
on fresh, matched responses; it does not identify new minimum sample sizes
over the full original grid. Planning PASS and inferential diagnostics should
continue to be read separately.

The fresh comparison used `joint_delta` intervals, with 2,048 null and
2,048 planned trials per design,
with the same responses for previous and corrected executables. Each person
answers 30 words. The corrected results are:

| Profile; N×J | Null rejection (95% MC interval) | Planned coverage (95% MC interval) | Detection AND width (95% MC interval) |
| --- | --- | --- | --- | --- |
| 8pp, intercepts only; 54×30 | 4.35% (3.54–5.32) | 95.41% (94.42–96.23) | 94.68% (93.62–95.57) |
| 8pp, moderate slopes; 72×60 | 3.96% (3.19–4.89) | 95.95% (95.00–96.72) | 81.79% (80.06–83.40) |
| 5pp, moderate slopes; 162×60 | 4.98% (4.12–6.01) | 95.26% (94.26–96.10) | 85.25% (83.65–86.72) |
| 8pp, larger variances; 72×90 | 5.86% (4.92–6.96) | 95.46% (94.47–96.28) | 89.55% (88.15–90.80) |
| 8pp, larger variances; 162×30 | 6.10% (5.15–7.22) | 94.48% (93.41–95.39) | 92.77% (91.57–93.82) |
| 8pp, larger variances; 216×90 | 5.52% (4.61–6.59) | 94.17% (93.07–95.11) | 99.51% (99.10–99.73) |

Coverage uses valid fits; rejection and joint attainment use all 2,048
attempts. There are 24,556/24,576 valid corrected fits, compared with all
24,576 previously. The 20 iteration-limit failures are retained: 13 null
and seven planned at the larger-variance 216×90 design. Every observed cell
valid rate exceeds 99%. The six layouts meet the prespecified practical
checks: coverage interval lower bound >=92.5% and null-rejection interval
upper bound <=7.5%. These bounds are more permissive than exact 95%/5%
calibration. The 162×30 null rate still has its interval above nominal 5%.

All six still pass the joint planning target (95% interval width <=10pp,
with detection/width joint-rate lower bound >=80%). However, moderate-slope
72×60 falls from 84.57% joint attainment to 81.79%, with a lower bound of
80.06%. It is a borderline planning pass, not a robust minimum. The old
60-design grid has not been rerun with this correction, and no new universal
participant/item recommendation follows from these six selected layouts.

The numerical correction has a cost: in separate matched timing checks,
median cached-fit CPU time rises from 0.053 to 0.115 seconds at 72×60 and
from 0.205 to 0.428 seconds at high-variance 216×90. These are per-fit CPU
times on the recorded machine, excluding generation and setup, not a
promise of total simulation wall time. The small intercept-only design
remains around 0.002–0.003 seconds per fit. CPU timings also vary with
concurrent load and processor frequency.

[Detailed protocol, executable provenance and paired results](https://github.com/Ryuya-dot-com/glmmplanr/tree/main/validation#l2-uncertainty-correction)
are retained with the earlier results.

<a id="default-bootstrap-confirmation"></a>

### Confirming the Default Probability-Scale Bootstrap

The default uses `prob_se_method = "bootstrap"` with 99 inner refits.
The corrected solver was separately checked on three of the preceding
profiles, with 512 new null and 512 new planned datasets per layout.
Both methods use the same responses and point estimates. All 3,072 outer
fits and 361,728 retained inner refits were valid. This confirms operation;
it does not establish nominal coverage or attainment of every planning goal.

| Profile; N×J | Default null rejection (95% MC interval) | Default planned coverage (95% MC interval) | Default detection AND width <=10pp (95% MC interval) |
| --- | --- | --- | --- |
| Intercepts only; 54×30 | 6.25% (4.46–8.69) | 94.92% (92.66–96.51) | 91.99% (89.32–94.04) |
| Moderate slopes; 72×60 | 3.71% (2.39–5.72) | 95.51% (93.35–96.99) | 77.15% (73.32–80.58) |
| Larger variances; 162×30 | 7.81% (5.79–10.46) | 93.75% (91.31–95.54) | 90.23% (87.36–92.51) |

Every rate uses 512 datasets per condition and pointwise Wilson intervals.
At 72×60, the same-response joint-delta attainment is 84.77%
(81.39–87.62), compared with 77.15% for the default. The paired difference
is -7.62pp (MCSE 1.68pp). The bootstrap result is indeterminate against the
80% planning target; a joint-delta pass cannot be transferred to it.
Average SEs are very similar (RMS ratio 0.998), but the per-dataset SE ratio
has empirical 2.5th–97.5th percentiles of 0.863–1.130. Similar average
interval widths therefore do not guarantee the same threshold decisions.

The default's coverage and null checks pass their practical bounds at
72×60. The null check is indeterminate at 54×30. Both checks are
indeterminate at 162×30, whose null interval is above nominal 5%.
These limitations remain visible alongside the two joint-planning passes;
all valid fits do not imply adequately calibrated inference in every profile.

The first 32 planned datasets per layout also use 399 refits. Joint
attainment changes from 31/32 to 31/32, 22/32 to 26/32, and 27/32 to
30/32 respectively. This small subset demonstrates sensitivity to inner
simulation, not superior coverage or justification for a new default.
The defaults remain unchanged. Specify the same interval method and
`control = target_control(prob_boot = ...)` in final planning confirmation
as intended for the analysis. No new minimum participant/item recommendation
follows from these selected layouts.

**What the saved draws explain.** Reanalysis of the same 3,072 datasets
separates interval computation from model fitting. At 72×60, 54 of the 58
joint-delta successes lost by the default bootstrap still detect the effect
but exceed the 10pp width limit. Four disjoint 99-draw batches give 22–25
joint successes among the same 32 long-draw datasets; all 399 draws give 26.
This shows inner-simulation sensitivity, not a confirmed benefit or an
80% guarantee for 399 refits. At the 162×30 null, RMS bootstrap SE matches
the empirical SD closely (ratio 1.002), but the lowest-SE fifth has 17/102
exclusions. Aggregate SE agreement therefore does not establish adequate
calibration of dataset-specific intervals.

Exploratory percentile/basic intervals reconstructed from the saved draws
are not public interval options. With 99 refits, percentile intervals reduce
162×30 null rejection from 40/512 to 32/512 (95% MC interval 4.46–8.69%),
while 72×60 joint attainment falls from 395/512 to 318/512; planned coverage
there remains 489/512. These tradeoffs do not justify an automatic method
switch. The results use the covariance-corrected executable before the
mode-roundoff repair and are not new independent confirmation of current
inference. Defaults and the need to match planning/analysis settings remain
unchanged. [Detailed saved-draw audit](https://github.com/Ryuya-dot-com/glmmplanr/tree/main/validation#interval-mechanism-audit).

[Protocol, paired results, retained inner counts and figure](https://github.com/Ryuya-dot-com/glmmplanr/tree/main/validation#default-bootstrap-confirmation).

<a id="rasch-mml-acceleration"></a>

### Slow Fixed-Bank MML and the Variance Boundary

The native solver now attempts safeguarded squared EM extrapolation after
the ordinary iteration limit. It accepts an exact zero-variance solution only
when the profile right derivative is nonpositive and its likelihood is not
lower. Fallback convergence also checks the variance score, since tiny EM
updates near zero can misleadingly suggest convergence.

Replaying 457 previously capped fits recovered 456 (99.8%); one near-zero
positive-variance case still reached the limit and remains a miss. An initial
update-only stopping rule returned all 457, but independent likelihood checks
revealed premature convergence; that rule was rejected. The 112 previously
successful controls retained exactly the same item RMSE and person SD.
On 564 common valid pairs with archived TAM results, all bank pass/fail
decisions agreed. Parameter and likelihood differences remain, especially
near the boundary; neither package is treated as the truth. Direct marginal
likelihood optimization and mirt checks are saved separately. This repair
does not validate crossed-GLMM intervals or all possible item banks.

<a id="cross-field-designs"></a>

## Scope Across Research Designs and Units of Assignment

In experiments and assessments where multiple participants respond (`0` / `1`) to a shared pool of items, experimental conditions can be assigned within or between participants and within or between items, with either complete or partial item presentation per participant. Both accuracy in psycholinguistic / cognitive tasks and binary preference judgments in behavioral experiments fall within this crossed binary structure. By contrast, when interventions are assigned to intact schools or clinics, the number of clusters and cluster size form distinct hierarchical design axes. Cluster-randomized school trials, multi-rater scoring allocations, and longitudinal within-person time series with serial dependence are outside the scope of `glmmplanr`'s built-in sample-size optimizers; probit models are covered via external-estimator comparisons.

<a id="validation-status"></a>

### Validation Status by Design and Estimator

Because `glmmplanr` supports more design structures than have undergone large-sample external confirmation, the table below separates **floor-checked / calibrated** benchmark regimes from **exploratory** regimes where users should explicitly verify interval coverage, null Type I error, and estimator agreement at their target design via `validate_design_point()`.

Distinguish two validation tiers below:

- **Floor-checked**: Pre-specified confirmation ruled out large calibration breakdown (e.g., exact 95% CI lower bound on empirical 95% coverage $\ge 90\%$, width attainment $\ge 80\%$); null Type I error is reported with its 95% CI rather than equivalence-tested tightly around 5% (e.g., at $N = 96\text{–}192$, exact 95% upper bounds on a 4–5% rate remain around 7–11%).
- **Calibrated (high-$N$ Type I error gate)**: Pre-specified confirmation at $N = 512$ datasets per design with explicit formal gates on both coverage (exact 95% CI lower bound $\ge 92.5\%$) and strict null Type I error (exact 95% CI upper bound $\le 7.5\%$).

All validation scripts and trial-level CSVs are archived in the public repository's [`validation/` directory](https://github.com/Ryuya-dot-com/glmmplanr/tree/main/validation) (excluded from CRAN source tarballs to keep installs compact; 3AFC confirmation CSVs are also shipped inside `inst/validation/`). Re-running any script via `sh validation/run-isolated.sh <script>` writes `environment.txt` and execution logs to an isolated temporary directory.

| Design / Model Regime | Estimators Benchmarked | Validation Status | Empirical Basis & Caveats |
|---|---|---|---|
| Complete-presentation canonical logit (`chance = 0`), random intercepts, marginal pp & logit contrasts | `bayes_laplace`, `lme4::glmer`, `simr`, `PowRPriori` | **Floor-checked** (benchmarked regimes) | 96-dataset paired coverage confirmation (`95/96` native, `94/96` `lme4`), 96-dataset paired null checks (between-subject & Latin square; Type I error 3.1–8.3%, wide 95% CIs), and 192-dataset paired calibration against `simr`/`PowRPriori`. |
| Complete-presentation $2 \times 2$ Latin-square 3AFC (`chance = 1/3`), correlated random intercepts & slopes, conditional logit | `bayes_laplace`, `rstan` HMC | **Calibrated (strict null, $N=512$) / Floor-checked (planned)** | 40-dataset paired HMC decomposition (`48 x 32`), 192-dataset sufficiency (`72 x 32`, `48 x 48`) and boundary (`32 x 32` to `48 x 40`) checks, plus 512-dataset strict conditional-null ($\beta = 0$) confirmation (Type I error `3.3% [1.9%, 5.3%]` and `3.5% [2.1%, 5.5%]`; slightly conservative). Laplace SEs omit SD-SD covariances and correlation uncertainty; cross-check with `stan_hmc` if HMC is the final analysis. |
| Random-item 1PL Rasch recovery (complete presentation, `chance = 0`) | `cpp_laplace`, `lme4` | **Floor-checked** (benchmarked regimes) | 32-dataset recovery checks (`144 x 12` vs. `72 x 24`) and 64-dataset independent confirmation (`62/64` joint attainment at `144 x 24`). |
| Fixed-bank 1PL Rasch calibration (complete presentation, including extreme anchor items) | `mml`, `TAM::tam.mml`, `cpp_laplace`, `lme4` | **Floor-checked** (benchmarked regimes) | 64-dataset paired cross-validation against `TAM` on a 16-item bank with extreme anchors ($\pm 3.5$ logit); constant-response items retained in denominator. |
| Partial item presentation (`items_per_subject`), random intercepts (`chance = 0`) | `bayes_laplace`, `lme4::glmer` | **Exploratory** (`indeterminate` at $N=96$ vs. 90% floor) | On `64 x 20` (5 items/person = `320` trials, between-subject), 96-dataset paired confirmation gave `89/96` (`92.7% [85.6%, 97.0%]`) planned and `90/96`–`91/96` (`93.8–94.8%`, Type I error `5.2–6.2%`) null coverage under `bayes_laplace`, matching `lme4::glmer` (`88/96` planned, `48/96` singular ML fits); verify coverage for your specific schedule via `validate_design_point()`. |
| Sparse / partial presentation with **correlated random slopes**, or higher-order interactions (`A:B:C`) | `bayes_laplace`, `cpp_laplace` | **Exploratory** (schedule-dependent) | On within-unit `64 x 24` Latin square (`8` items/person = `512` trials, `slope_sd_ratio = 0.30`, `corr = 0.20`), `bayes_laplace` achieved `94/96`–`95/96` (`97.9–99.0%`) coverage and `2.1%` null detection (`pass`, conservative Half-$t$ regularization), but sparser schedules or `chance > 0` remain unestablished; verify with `validate_design_point()`. |
| Passage / testlet nesting (`design_testlet()`) and participant-by-testlet local dependence | `bayes_laplace` with `joint_delta`; pilot `lme4::glmer` | **Exploratory** (allocation- and variance-dependent) | In a 1,024-dataset confirmation per scenario at `48 x 24`, `6 x 4` and `12 x 2` passages × questions gave 12.7% and 8.0% null detection with passage SD 0.5 and participant-by-passage SD 0.4. Planned coverage was 87.6% and 92.5%, respectively. The `12 x 2` design passed the coverage floor and had an indeterminate Type I error gate. |

At the same 1,152 responses, a single passage per arm (`2 x 12`) gave 41.0% null detection and 60.7% planned coverage in the same benchmark. With weaker shared passage variation (SD 0.25), `6 x 4` and `12 x 2` both passed coverage and null gates. A separate `3 x 8` benchmark shared all passages across arms and randomized readers, with 94.9% coverage and 5.9% null detection under a common effect across passages. All tested designs fell below the 80% power target. The [allocation benchmark](https://github.com/Ryuya-dot-com/glmmplanr/tree/main/validation#4b-passage-allocation-at-fixed-response-count) uses binary logit with full presentation and lists the exact Monte Carlo intervals for each scenario.

When both arms use the same passages, state whether the interval describes the average for those passages or for a population of passages. In a shared-passage benchmark with 48 readers, 24 questions and heterogeneous condition effects, the fixed passage-by-condition analysis with two passages covered the presented-passage average in 244/256 runs (95.3%) and the population average in 197/256 (77.0%). In the population-null scenario, adding a random passage slope reduced false detection from 63/256 (24.6%) to 12/256 (4.7%). At the planned ten-point effect, it also increased mean interval width from 18.26 to 28.65 percentage points. All fits were valid and the results varied with the effect heterogeneity and the intended scope of inference. The benchmark uses the native kernel and research code to specify the fixed passage effects. [Shared-passage benchmark and Monte Carlo intervals](https://github.com/Ryuya-dot-com/glmmplanr/tree/main/validation#4h-fixed-and-random-passage-effects-with-two-to-six-shared-passages)

In the independent external comparison, the two-passage, 48-reader population-null scenario with heterogeneous slopes gave 18.8% and 18.5% null detection using only random passage intercepts under the native engine and `glmmTMB`. With passage slopes, the rates were 5.8% and 12.5%. At 12 passages and the same 1,152 responses, the slope-model rates were 5.4% and 5.7%, planned coverage was 94.5% and 94.2%, and planned power was 52.1% and 54.7%, respectively. A separate analysis of the two presented passages and their actual items also used `TAM` and `mirt`. Fixed-form point-estimate RMSE was similar on the common valid datasets. Interval availability varied across methods, with `mirt` returning 379/512 intervals. [External validation and exact Monte Carlo intervals](https://github.com/Ryuya-dot-com/glmmplanr/tree/main/validation#4i-shared-passage-comparisons-with-external-irt-and-glmm-estimators)

In an exploratory fixed-bank comparison with 128 datasets per scenario, the native common-discrimination calibration agreed with `mirt` Rasch to within 0.00082 standardized difficulty units. With 240 readers, eight items and equal true discrimination, both met the difficulty criteria in 122/128 runs. With wide discrimination heterogeneity at 480 readers and eight items, `TAM` and `mirt` 2PL met the difficulty criteria in 90/128 runs, and the Rasch fits met them in 0/128. The likelihood diagnostics also found extreme positive loadings in the rank-one Laplace GLMM under short tests. This benchmark holds the response data and ability scale constant and compares the item curves with the known truth. [2PL comparison and reproduction](https://github.com/Ryuya-dot-com/glmmplanr/tree/main/validation#4j-fixed-bank-2pl-calibration-with-glmm-and-irt-estimators)

For the same 16-item bank and 240 responses per item, 240 readers answering all items gave 240 co-observations per item pair, and 480 readers answering eight items in balanced booklets gave 112. With wide discrimination heterogeneity, the joint difficulty and discrimination criterion under mirt and TAM 2PL was met in 27/128 complete runs and 5/128 booklet runs. Under equal true discrimination, the Rasch difficulty criterion was met in 124/128 and 125/128 runs. An independent 2PL workflow in `irtsim` using the same booklet matrix also gave lower joint attainment for the booklets. [Equal-exposure comparison and reproduction](https://github.com/Ryuya-dot-com/glmmplanr/tree/main/validation#4k-complete-and-incomplete-presentation-at-equal-item-exposure)

In a complete-presentation vocabulary-learning benchmark with 2,048 datasets per generating model, the relative precision of 32 participants × 48 words and 96 × 16 reversed when the dominant source of treatment-effect variation switched between participants and words. With words rotated across conditions, 48 × 32 passed coverage, null detection, and 80% power gates for an eight-point marginal recall effect under balanced slopes. The three 3,072-response layouts passed across all four core slope profiles. The highest power for a five-point effect was 71.7%. These checks use binary logit, `bayes_laplace`, and `joint_delta`. Numerical results for complete presentation also include participant and word assignment, correlation sensitivity, and a separate internal example with continuous covariates. [L2 vocabulary benchmark](https://github.com/Ryuya-dot-com/glmmplanr/tree/main/validation#4c-l2-vocabulary-recall-and-the-allocation-of-participants-and-words)

A complementary L2 vocabulary benchmark uses 3,072 datasets per generating model and 1,536 scheduled responses per dataset. At the word-dominant slope setting, spreading 96 participants' responses across 48 words in linked lists reduced mean interval width from 12.67 to 10.81 percentage points compared with giving everyone the same 16 words. Linked and disjoint lists gave similar precision for the marginal method effect under the common-form assumptions. An allocation with identical participant and word margins still produced about 12.9 points of bias when compatible participant–word pairs were overrepresented in one arm. Modeling observed compatibility corrected that bias in the benchmark. Outcome-dependent missingness left bias after adjustment. This benchmark uses precomputed balanced schedules and internal interfaces to the native estimator. [Partial-presentation and compatibility benchmark](https://github.com/Ryuya-dot-com/glmmplanr/tree/main/validation#4d-partial-presentation-and-participantword-compatibility)

In a related 2,048-dataset-per-scenario benchmark, a 48-word bank supplied initial 16-word lists to 96 participants. When known-word rates were 30% at −1 vocabulary SD and 70% at +1 SD, topping up to 16 words under perfect screening raised eight-point unknown-word power from 56.7% to 82.1%. The top vocabulary quartile screened an average of 46.1 words to learn 12.0 genuinely unknown words. With 90% sensitivity and 95% specificity, the same quartile learned 10.7 genuinely unknown words, and population unknown-word power fell to 73.9%. The analysis uses a common generating reference distribution for the unknown-word contrast and internal covariate matrices. [Prior knowledge and screening benchmark](https://github.com/Ryuya-dot-com/glmmplanr/tree/main/validation#4e-prior-word-knowledge-and-vocabulary-screening)

In a floor and ceiling vocabulary benchmark with 2,048 datasets per scenario, 96 learners studied 16 words each from a 48-word bank. At reference accuracies of 10%, 42%, and 90%, detection of a fixed five-point effect was 78.8%, 48.6%, and 93.7%. Holding the conditional log odds ratio constant gave smaller true effects at the endpoints and detection of 19.7%, 48.6%, and 17.7%. The benchmark compares 14 allocations and also examines concentrated learner and word extremes, random-slope and correlation sensitivity, and measured covariates. It reports counts of all-correct and all-wrong learners or words alongside the observations per unit. [Floor and ceiling benchmark](https://github.com/Ryuya-dot-com/glmmplanr/tree/main/validation#4f-vocabulary-recall-at-the-floor-and-ceiling)

<a id="modeling-assumptions"></a>

### Generating-Model Assumptions and Custom Misspecification Testing

Built-in simulation and calibration (`calibrate_dgp()`, `simulate_dgp_data()`, `optimize_sample_size()`) rely on the following structural assumptions:

1. **Multivariate normal random effects**: Participant and item random intercepts and slopes are drawn from zero-mean multivariate normal distributions (`MASS`/Cholesky draws). Heavy-tailed, skewed, or mixture ability/difficulty distributions are not generated by the built-in DGP.
2. **Shared guessing floor (`chance`) and zero upper lapse**: All items share one fixed lower asymptote $g = \text{chance}$ and an upper asymptote of $1.0$ (no item-specific guessing parameters, no 2PL/3PL discrimination heterogeneity, and no attentional lapses).
3. **Missingness independent of ability and responses (MCAR)**: Trial-level missingness (`attrition`), whole-participant dropout (`participant_dropout`), and partial item assignment (`items_per_subject`) are sampled independently of latent random effects and trial outcomes. Ability-dependent dropout (MNAR) or adaptive routing is outside the built-in DGP.
4. **Scope of `scenario_dgps()` and custom `data_fn` testing**: Multi-scenario planning varies effect sizes, variance components, or missingness rates *within this same parametric family*. Passing across multiple scenarios does not guarantee robustness when the structural assumptions above are violated. To stress-test structural misspecification (e.g., ability-dependent dropout or item-specific guessing) directly under your planned analysis model, pass a custom `data_fn(schedule, spec)` function to `evaluate_design_point()` or `validate_design_point()`:

```r
# Example: evaluate robustness to ability-dependent (MNAR) dropout via data_fn
mnar_dropout <- function(schedule, spec) {
  df <- simulate_dgp_data(spec$dgp, schedule = schedule)
  subj_acc <- tapply(df$y, df$subj_id, mean)
  keep <- as.integer(names(subj_acc)[subj_acc >= stats::quantile(subj_acc, 0.20)])
  df[df$subj_id %in% keep, , drop = FALSE]
}
# evaluate_design_point(spec, n_subj = 48, n_item = 24, nsim = 100, data_fn = mnar_dropout)
```

### Sparse responses and the limits of a converged fit

Check the number of observed responses for each learner and each item. A correct or incorrect response contributes even when it is the only one recorded for that unit. In a three-choice task at chance, five independent responses are all wrong with probability $(2/3)^5 = 0.132$. With 32 responses the probability is about 0.00000232. These are upper bounds under the fitted guessing model with conditionally independent trials. A mechanism that forces a particular item to be wrong on every presentation would fall outside that model.

The ceiling and the guessing floor are two different limits. Write the conditional correct probability as $p=g+(1-g)s$, with $g=1/m$ and $s=\operatorname{logit}^{-1}(\eta)$. The contribution of one binary response to expected information about $\eta$ is $[(1-g)s(1-s)]^2/[p(1-p)]$. This approaches zero near either limit. Few responses per learner, few per item, or little within-unit variation in condition can leave variances and slope correlations weakly determined. First-order Laplace approximation uses local curvature, and increasing the total trial count while keeping each unit sparse can leave appreciable approximation error. Its direction and size depend on the design. [Joe (2008)](https://doi.org/10.1016/j.csda.2008.05.002)

In the initial 32-dataset partial-presentation pilot (`64` people each seeing `5` of `20` items = `320` trials, between-subject, random intercepts, `chance = 0`, true marginal difference $-5$ percentage points, `29` bootstrap refits), `bayes_laplace` covered the truth in `30/32` datasets (`93.8%`, exact 95% CI `79.2–99.2%`). A subsequent pre-specified 96-dataset paired confirmation (`validation/sparse_and_testlet_recheck.R`, `99` bootstrap refits, seeds `30001–31096`) yielded `89/96` (`92.7%`, `[85.6%, 97.0%]`) planned coverage and `6/96` (`6.2%`, `[2.3%, 13.1%]`) null Type I error under `bayes_laplace` (`5/96` with `joint_delta`), closely matching `lme4::glmer` (`88/96 = 91.7%` planned coverage; `48/96` singular ML fits), while leaving the 90% floor `indeterminate` at $N = 96$. On a `64 x 24` within-unit Latin-square schedule (`8` items/person = `512` trials) with correlated random slopes (`slope_sd_ratio = 0.30`, `corr = 0.20`), `bayes_laplace` achieved `94/96`–`95/96` (`97.9–99.0%`) coverage (`2/96` null detection), whereas a between-passage testlet design with only `6` passages (`48 x 24 = 1,152` trials) suffered small-cluster Wald undercoverage under both `lme4::glmer` (`20/96 = 20.8%` null Type I error, `below_floor`) and `bayes_laplace` (`14/96 = 14.6%`), worsening further under unmodeled participant-by-passage local dependence (`SD = 0.4`, planned coverage `76.0–82.3%`, `below_floor`). Always inspect `n_valid`, coverage, and null detection from `validate_design_point()` for your specific sparse or clustered schedule, and note that `bayes_laplace` reports boundary singularity as not applicable (`NA`) because Half-$t$ priors prevent zero-variance boundary collapse even when data provide little variance information.

<a id="coverage-recheck"></a>

## Marginal Percentage-Point 95% Interval Coverage

We evaluated 64 participants (32 per group) responding to all 20 items under a true population-averaged group difference of $-5$ percentage points ($60\%$ vs. $55\%$), with random intercept SDs of `0.5` (participants) and `0.3` (items) on the conditional logit scale and 99 parametric bootstrap refits per dataset. Across an initial pilot set of 32 datasets, the built-in engine covered the true contrast in `25/32` datasets; in an independent pre-specified confirmation set of 96 datasets, coverage was `95/96` (`98.96%`, exact binomial 95% CI `[94.33%, 99.97%]`). On the exact same 96 datasets, `lme4::glmer` achieved `94/96` coverage (`[92.68%, 99.75%]`). All 96 intervals converged validly under both estimators and exceeded the pre-specified `90%` lower-bound benchmark.

Under the paired null condition ($60\%$ vs. $60\%$ on the same $64 \times 20$ schedule, 96 datasets per allocation), Type I error rates for the C++ batch parametric bootstrap (`prob_boot = 99`) were `5/96` (`5.2%`) for between-participant assignment and `7/96` (`7.3%`) for within-participant Latin-square assignment. Without refitting, the analytical joint delta method (`prob_se_method = "joint_delta"`, propagating both fixed-effect curvature and variance-component uncertainty) yielded `3/96` (`3.1%`) and `8/96` (`8.3%`) on the identical datasets, agreeing in mean interval width within `0.07` percentage points. Because candidate `pass` status in `optimize_sample_size()` reflects user-specified power or width targets, coverage under alternative sample sizes, missingness rates, or partial-presentation schedules should be verified for the planned design via `validate_design_point()`.

<a id="three-choice-precision"></a>

## Conditional-Logit Uncertainty in 3-Choice (3AFC) Tasks

In a paired benchmark across 40 simulated 48-by-32 datasets fitted using both built-in Laplace and RStan HMC estimators, fixing covariance at the generating values yielded a mean Laplace SE of 0.201 logit versus a mean posterior SD of 0.204 (38 valid HMC fits). With covariance estimated, the mean joint Laplace SE was 0.227 and the mean HMC posterior SD was 0.270 (40 paired fits). In `src/glmm_laplace.cpp`, the analytical variance-uncertainty correction propagates marginal log-SD curvature diagonally (`J_beta_theta * diag(var_log_tau) * J_beta_theta'`), omitting off-diagonal covariances among SD estimates and estimation uncertainty in random-slope correlation parameters. Because this diagonal approximation and full joint HMC differ in covariance propagation, numerical integration, and regularization, finite-sample interval widths can still differ between the two estimators.

When four candidate $(N_{\text{subj}}, N_{\text{item}})$ designs were compared on the same within-participant $2 \times 2$ 3AFC task (true conditional coefficient `1.616` logit) across 40 paired datasets per candidate, the smallest scheduled trial count achieving $\ge 80\%$ attainment of full 95% interval width $\le 1.00$ logit under RStan HMC was `2,304` trials ($72 \times 32$ and $48 \times 48$).

**Sufficiency confirmation (`72 x 32` and `48 x 48`, 192 datasets each).** Evaluating the joint-covariance `bayes_laplace` engine on 192 independent confirmation datasets at $72 \times 32$ and $48 \times 48$ yielded `192/192` valid fits and `192/192` width attainment ($\le 1.00$ logit) at both designs under the planned effect (`1.616` logit), with empirical 95% coverage of `183/192` (`95.3%`, exact 95% CI `[91.3%, 97.8%]`) at $72 \times 32$ and `187/192` (`97.4%`, `[94.0%, 99.1%]`) at $48 \times 48$. Under the strict conditional null generated via `null_target_spec()` (`conditional_null`, where the true conditional exposure coefficient is identically `0` logit), `bayes_laplace` produced `192/192` valid fits at both $72 \times 32$ and $48 \times 48$, with 95% coverage of `186/192` (`96.9%`, `[93.3%, 98.8%]`) and `184/192` (`95.8%`, `[92.0%, 98.2%]`), width attainment of `192/192` at both designs, and strict Type I error rates (zero excluded, all 192 datasets in the denominator) of `6/192` (`3.1%`, `[1.2%, 6.7%]`) at $72 \times 32$ and `8/192` (`4.2%`, `[1.8%, 8.0%]`) at $48 \times 48$. In the historical CSV condition labeled `null` (equal marginal accuracies of 59% and 73% within each spacing condition, corresponding to a small nonzero conditional exposure coefficient of `+0.004097` logit), zero was excluded in `6/192` (`3.1%`, `192/192` valid) and `10/192` (`5.2%`, `192/192` valid) trials, respectively.

**High-$N$ strict conditional-null calibration ($N = 512$ datasets per design).** Because 192 replications leave wide exact 95% upper bounds on a 5% Type I error rate (6.7–8.0%), a separate pre-specified confirmation (`validation/three_choice_boundary_and_null.R`, seeds `7711001–7712512`) evaluated `512` fresh strict conditional-null ($\beta = 0$) datasets at `72 x 32` and `48 x 48` against tighter gates (exact 95% Type I error upper bound $\le 7.5\%$ and coverage lower bound $\ge 92.5\%$). At $48 \times 48$, `512/512` fits were valid, coverage was `495/512` (`96.7%`, `[94.7%, 98.1%]`), and Type I error was `17/512` (`3.32%`, `[1.95%, 5.26%]`). At $72 \times 32$, `511/512` fits were valid, coverage was `493/511` (`96.5%`, `[94.5%, 97.9%]`), and Type I error was `18/512` (`3.52%`, `[2.10%, 5.50%]`). Both designs passed the tighter gates, showing slightly conservative calibration under the strict null. Note that the `<= 7.5%` upper-bound gate rules out gross Type I error inflation rather than proving exact equivalence to `5.0%`, and this confirmation applies to `72 x 32` and `48 x 48` under the stated variance/correlation settings (not automatically to `72 x 24` or other regimes).

**Necessity and pass/fail boundary check (`32 x 32` to `48 x 40`, 192 datasets each).** Because `192/192` width attainment at `2,304` trials confirms sufficiency rather than minimum required sample size, six smaller adjacent candidates were evaluated under the same planned DGP (`validation/three_choice_boundary_and_null.R`, seeds `7611001–7616192`). Width attainment ($\le 1.00$ logit) fell cleanly below the 80% floor at `32 x 32` (`1,024` trials: `10/192 = 5.2%`, `[2.5%, 9.4%]`) and `48 x 24` (`1,152` trials: `25/192 = 13.0%`, `[8.6%, 18.6%]`), straddled the boundary as `indeterminate` at `48 x 32` (`1,536` trials: `163/192 = 84.9%`, `[79.0%, 89.6%]`) and `64 x 24` (`1,536` trials: `163/192 = 84.9%`, `[79.0%, 89.6%]`), and passed at `72 x 24` (`1,728` trials: `187/192 = 97.4%`, `[94.0%, 99.1%]`) and `48 x 40` (`1,920` trials: `191/192 = 99.5%`, `[97.1%, 100.0%]`). Across five independent search seeds (`nsim_final = 96`, `prefer = "fewer_trials"`), `optimize_sample_size()` achieved an internal Stage-2 `pass` on `5/5` seeds, selecting `72 x 24` (`1,728` trials) on 4/5 seeds and `48 x 32` (`1,536` trials) on 1/5 seed. Checked against the independent 192-rep boundary evaluations above, 4/5 selected designs (`72 x 24`) confirmed as independent `pass`, 1/5 (`48 x 32`) was `indeterminate` on the boundary, and 0/5 were `below_floor`. When RStan HMC is the planned final estimator or when small-sample candidates are considered, verify candidate width and coverage directly under the target estimator (e.g., `engine = "stan_hmc"`).

The 0.1.0 confirmation covers all 1,152 trial seeds (`planned`, `null`, and `conditional_null`; 192 datasets per condition and candidate design) in R 4.6.1. These are operating characteristics for the stated generating model, analysis, and two candidate designs. The following code reads the trial-level results shipped with the package. Each row identifies the seed, true coefficient, estimate, SE, and fit validity.

```r
path <- system.file("validation", package = "glmmplanr")
fits <- read.csv(file.path(path, "three_choice_confirmation_fits.csv"))
counts <- read.csv(file.path(path, "three_choice_confirmation_summary.csv"))
counts[c("effect", "persons", "items", "attempted", "valid", "covered", "narrow", "detected")]
```

To repeat the fits, pass the installed `three_choice_confirmation.R` script and an existing output directory to `Rscript --vanilla`. The script writes fresh trial and summary CSV files there. Use an empty output directory to preserve earlier runs.

When comparing two runs on the same seeds, match each trial by effect, design and seed. Check validity on each matched trial and report the largest absolute differences in the coefficient and SE among pairs of valid fits, in logits. Recompute coverage, detection and width attainment from the unrounded values. Inspect any trial whose outcome changes near a threshold, even if the numerical difference is small. For a fresh set of seeds, use Monte Carlo intervals to interpret differences in rates.

<a id="rasch-recovery"></a>

## Binary Rasch and Multilevel Parameter Recovery

Using a 1PL Rasch logit model with discrimination `1`, lower asymptote `0`, and population mean person ability `0`, we drew fresh persons and items on each of 32 replications at a fixed total of `1,728` responses. Mean ability RMSE was `0.518` at $144 \times 12$ vs. `0.416` at $72 \times 24$, whereas mean item difficulty RMSE reversed to `0.187` vs. `0.275`. A four-variance-component multilevel illustration (12 schools, 36 classrooms) was also verified against external `lme4` estimation; selecting school-level sample sizes requires dedicated hierarchical schedules ([De Boeck et al., 2011](https://www.jstatsoft.org/article/view/v039i12)).

<a id="rasch-precision"></a>

## Generalizing Across Newly Sampled Examinees and Items

Assuming complete presentation without missingness, `65%` marginal accuracy, person ability SD `0.8`, and item difficulty SD `0.6`, we required simultaneous attainment of ability RMSE $\le 0.50$, difficulty RMSE $\le 0.30$, and item population SD absolute error $\le 0.20$ logit. At $144 \times 24$, all three tolerances were simultaneously satisfied in `63/64` replications and in `62/64` independent confirmation replications (95% Wilson CI `[89.3%, 99.1%]`). Widening the true item difficulty SD to `1.0` reduced joint attainment at $144 \times 24$ to `47/64` (`[61.5%, 82.7%]`).

<a id="rasch-fixed-bank"></a>

## Calibrating a Fixed 16-Item Operational Bank

Fixing 14 core item difficulties evenly between $-1$ and $+1$ logit plus two extreme anchor items at $\pm 3.5$ logit (ability SD `0.8`), 64 replications at $N_{\text{subj}} = 100$ achieved bank-wide difficulty RMSE $\le 0.50$ in `64/64` replications under random-item Empirical Bayes shrinkage and `62/64` under `TAM` fixed-item Marginal Maximum Likelihood (MML), whereas simultaneous worst-item accuracy ($|\hat{b}_i - b_i| \le 0.50$ across all 16 items) dropped to `19/64` and `15/64`. At $N_{\text{subj}} = 30$, the hardest anchor item yielded zero correct responses in `17/64` replications. Under unpenalized fixed-difficulty estimation, all-correct and all-incorrect items are retained in the evaluation denominator as non-attaining replications. `plan_rasch_fixed_bank()` plans examinee sample sizes for a fixed difficulty vector via MML (built-in `"mml"` or `"tam"`) or Empirical Bayes shrinkage (`"cpp_laplace"` or `"lme4"`), whereas `plan_rasch_precision()` resamples items across replications ([TAM MML reference](https://search.r-project.org/CRAN/refmans/TAM/html/tam.mml.html)).

A separate Rasch bank-expansion benchmark applies a common-item location shift after fixed-item MML. At 480 new examinees with 16 items per person and a 500-examinee reference calibration, four, eight, and twelve stable anchors spanning the difficulty scale yielded joint attainment of 78.0%, 82.5%, and 62.8%, respectively. Joint attainment includes limits on the origin error, new-item RMSE, and worst new-item error. The research script uses the built-in calibrator for the simulation and compares paired fits from `TAM`. The fixed-bank planning API identifies the origin by centering the ability distribution. [Common-item linking benchmark](https://github.com/Ryuya-dot-com/glmmplanr/tree/main/validation#4g-common-item-linking-for-rasch-bank-expansion)

A paired 16-item benchmark with glmmplanr 0.1.0, TAM and mfrmr compared
14 calibration and eight DIF conditions on 10,240 simulated datasets.
At 240 responses per item, the targeted 240-person complete and 480-person
eight-item booklet designs gave difficulty RMSEs of 0.161 and 0.162
logits and person EAP reliabilities of 0.648 and 0.479. With ability SD
0.4, the release MML iteration cap allowed 143/512 booklet fits. A
follow-up at a higher cap recovered 367/369 of those failed fits. For
one preselected item with uniform DIF of 0.70 logits and 15 invariant
anchors, TAM and mfrmr detected DIF in 329/384 complete and 213/384
booklet datasets at 480 persons. The benchmark computes EAP and
residual diagnostics from the fitted models and compares them with
package-reported values. The glmmplanr DIF comparison uses the
validation script's Laplace GLMM. [Rasch fit, separation, and anchored DIF benchmark](https://github.com/Ryuya-dot-com/glmmplanr/tree/main/validation#4l-rasch-calibration-fit-separation-and-anchored-dif)

<a id="related-packages"></a>

## Related R Packages by Research Goal

To situate `glmmplanr` 0.1.0 within the R ecosystem, the tables below summarize CRAN and R-universe packages that address statistical power, multi-axis design search, interval precision, or psychometric calibration for binary and multilevel data, based on each package's documented interface and primary methodological reference. Empirical benchmarks include the [matched irtsim/mlpwr item-pool comparison](#item-pool-package-comparison), marginal probability intervals (`lme4::glmer`), 3AFC conditional-logit intervals (`rstan`), fixed-bank Rasch calibration (`TAM`), and [probability-to-logit input calibration (`PowRPriori`, `simr`, and `mlpwr`)](#package-calibration).

`glmmplanr` is specifically designed for studies that start from substantive population-averaged accuracy differences (percentage points), lower guessing asymptotes, and crossed presentation schedules, and then identify non-dominated Pareto-optimal $(N_{\text{subj}}, N_{\text{item}})$ designs meeting power and interval-width criteria. For arbitrary numbers of hierarchical random factors, group-mean centering of continuous predictors, or free 2PL/3PL slope and guessing estimation, the packages below provide complementary capabilities.

### Power, Precision, and Design Search for Binary Mixed Models

| Package | Primary Inputs & Capabilities | Relationship & Complementarity with `glmmplanr` |
|---|---|---|
| **PowRPriori** | Defines crossed/nested structures via `define_design()`, converts expected binary probabilities into fixed effects, and increments a chosen design axis in `power_sim()`, including within/between predictor centering. | Closest counterpart for canonical binary GLMMs. Standard vignettes apply direct logit transformation, whereas `glmmplanr::calibrate_dgp()` inverts Gauss-Hermite quadrature over the combined cell normal variance so marginal probabilities match under the specified random-effect variances ([PowRPriori workflow](https://CRAN.R-project.org/package=PowRPriori/vignettes/Workflow-Example.html)). |
| **simr** | Evaluates power via `powerSim()` and 1D `powerCurve()` from either fitted `lme4` models or synthetic `makeGlmer()` objects given an external schedule, fixed effects, and variance components. | Standard reference for canonical binomial GLMMs without pilot data when supplied with calibrated conditional parameters and explicit schedule data frames ([Green & MacLeod, 2016](https://doi.org/10.1111/2041-210X.12504); [from-scratch vignette](https://CRAN.R-project.org/package=simr/vignettes/fromscratch.html)). |
| **mixedpower** | Computes simulation-based power for fixed effects and interactions from a fitted `lme4` model and pilot dataset; `R2power()` scans a second random factor while holding the first at a specified size. | Direct precedent for participant-by-stimulus trade-offs when pilot data and fitted `lme4` models are available ([Kumle et al., 2021](https://doi.org/10.3758/s13428-021-01546-0); [manual](https://lkumle.github.io/power_notebooks/intro/Introduction_to_mixedpower.pdf)). |
| **mlpwr** | Accepts a user-written binary simulation function in `find.design()` to search multi-parameter design spaces via surrogate modeling, with IRT and multilevel examples. | Surrogate-search counterpart that can use glmmplanr's native simulation as its callback. The [matched comparison](#item-pool-package-comparison) applies the same independent Wilson confirmation to both methods; glmmplanr supplies the declarative schedules and effect/precision targets ([Zimmer et al., 2024](https://doi.org/10.3758/s13428-023-02269-0); [manual](https://CRAN.R-project.org/package=mlpwr/refman/mlpwr.html)). |
| **glmmrOptim** | Selects $c$-optimal observation subsets from candidate schedules by minimizing asymptotic linear-contrast variance via the GLMM information matrix. | Useful for optimal condition-level observation allocation; asymptotic variance criteria complement simulation-verified finite-sample CI-width assurance and coverage ([glmmrOptim manual](https://samuel-watson.r-universe.dev/glmmrOptim/doc/manual.html)). |
| **pilotr** | Generates crossed participant-by-item experimental datasets (including Bernoulli outcomes) and evaluates Type S/M errors, mean CI width, and ROPE metrics. | Precision-reporting counterpart; `power_mixed()` in v0.3.1 targets linear mixed models via `lmerTest`/Satterthwaite rather than binary GLMM refitting ([pilotr manual](https://CRAN.R-project.org/package=pilotr/refman/pilotr.html); [precision & ROPE vignette](https://CRAN.R-project.org/package=pilotr/vignettes/precision-rope.html)). |
| **simglm** | Simulates configurable regression/multilevel data-generating and fitting models, summarizing power, Type I error, and parameter precision via `compute_statistics()`. | Flexible engine for custom covariate distributions and misspecification stress tests when conditional coefficients are specified directly ([simglm manual](https://CRAN.R-project.org/package=simglm/refman/simglm.html)). |

### Assignment Units, Guessing Asymptotes, and Design Priors

| Package | Target Design Question | Role Relative to `glmmplanr` |
|---|---|---|
| **designr** | Constructs general crossed, nested, and counterbalanced schedules across arbitrary numbers of fixed and random factors (participants, stimuli, raters). | Reference for multi-factor schedule construction beyond participant-item-testlet hierarchies ([designr reference](https://search.r-project.org/CRAN/refmans/designr/html/designr.html); [overview](https://CRAN.R-project.org/package=designr/readme/README.html)). |
| **CRTSize** | Analytical sample size formulas (`n4props()`, etc.) for two-arm binary cluster-randomized trials given event rates, cluster size, and ICC. | Appropriate when intact schools/clinics are randomized; single-ICC design effects do not transfer directly to crossed participant-by-item variance structures ([n4props reference](https://search.r-project.org/CRAN/refmans/CRTSize/html/n4props.html)). |
| **quickpsy** | Fits psychometric threshold functions incorporating lower guessing and upper attentional lapse asymptotes. | Shares the lower-asymptote response formulation ($g = 1/3$, lapse $= 0$ under a logistic core) while targeting psychophysical threshold curves rather than crossed GLMM condition contrasts ([Linares & López-Moliner, 2016](https://journal.r-project.org/articles/RJ-2016-008/index.html); [reference](https://search.r-project.org/CRAN/refmans/quickpsy/html/quickpsy.html)). |
| **bayesassurance** | Computes Bayesian assurance (unconditional probability of study success) integrated over continuous design priors on unknown effect sizes. | `glmmplanr` evaluates conditional assurance at each scenario and worst-case multi-scenario robustness (`scenario_dgps()`); prior-weighted integration addresses unconditional assurance ([Pan & Banerjee, 2023](https://journal.r-project.org/articles/RJ-2023-066/)). |

### Rasch, IRT, and Educational Measurement

| Package | Target Measurement Question | Role Relative to `glmmplanr` |
|---|---|---|
| **irtsim** | Simulation-based IRT study planning with item-parameter recovery and missingness/booklet designs. | Independent planning workflow with a mirt backend. The [matched Rasch comparison](#item-pool-package-comparison) uses the same responses, scale and whole-bank success event ([manual](https://CRAN.R-project.org/package=irtsim/refman/irtsim.html)). |
| **irtpwr** | Analytical/asymptotic sample size calculation for IRT model-comparison and item-fit hypotheses (e.g., Rasch vs. 2PL). | Suited to hypothesis testing of IRT model structure, whereas `plan_rasch_precision()` bounds parameter RMSE and variance recovery tolerances ([irtpwr vignette](https://CRAN.R-project.org/package=irtpwr/vignettes/demo.html)). |
| **pwrRasch** | Simulation-based power to detect Rasch model misfit via group-by-item interaction tests. | Addresses DIF / model-fit testing rather than individual item calibration error bounds ([pwrRasch manual](https://CRAN.R-project.org/package=pwrRasch/refman/pwrRasch.html)). |
| **TAM** & **eRm** | Item calibration for operational item banks via Marginal Maximum Likelihood (`TAM`) or Conditional Maximum Likelihood (`eRm`). | `plan_rasch_fixed_bank()` cross-validates directly against `TAM::tam.mml()` on identical response matrices; `eRm` provides CML sensitivity checks ([TAM reference](https://search.r-project.org/CRAN/refmans/TAM/html/tam.mml.html); [Mair & Hatzinger, 2007](https://www.jstatsoft.org/article/view/v020i09)). |
| **mirt** & **eirm** | `mirt` simulates and estimates multidimensional 1PL–4PL IRT models; `eirm` fits explanatory IRT models linking item/person covariates to responses via `lme4`. | `mirt::simdata()` supports 2PL/3PL misspecification studies; `eirm` complements `glmmplanr`'s explanatory IRT sample-size planning at analysis time ([mirt::simdata](https://philchalmers.github.io/mirt/reference/simdata.html); [Bulut et al., 2021](https://doi.org/10.3390/psych3030023)). |
| **optical** | Optimal item-to-examinee block assignment using preliminary item parameters to maximize calibration information. | Useful for adaptive/targeted incomplete-block routing, whereas `items_per_subject` in `glmmplanr` plans uniform random incomplete-block sampling ([optical CRAN page](https://CRAN.R-project.org/package=optical)). |
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

Passing `cells$conditional_prob` yields intercept `1.13345` and group contrast `0.51625`, which can be supplied to `simr::makeGlmer(fixef = unlist(b), ...)` alongside random-intercept variances `0.64` and `0.36`. By contrast, transforming `0.72` and `0.80` directly via `qlogis()` yields slope `0.44183` and shrinks the true population-averaged cell probabilities to `68.68%` vs. `76.13%` (a `7.45 pp` rather than `8.00 pp` marginal difference; see [PowRPriori workflow](https://CRAN.R-project.org/package=PowRPriori/vignettes/Workflow-Example.html) and [simr from-scratch vignette](https://CRAN.R-project.org/package=simr/vignettes/fromscratch.html)).

Across 192 datasets on a complete $64 \times 24$ schedule (`1,536` observations, two-sided 5% Wald test), calibrating the 8-percentage-point marginal difference produced `188/192` (`97.9%`) detections for both `glmmplanr` and `simr` under a within-participant Latin square, and `115/192` (`59.9%`, `glmmplanr`), `111/192` (`57.8%`, `simr` on identical responses), and `112/192` (`58.3%`, independent `PowRPriori::power_sim()`) under between-participant assignment (`32` vs. `32`). Without quadrature calibration, between-participant detections dropped to `95`, `92`, and `96` of `192`, illustrating how conflating marginal and conditional probabilities alters both the simulated population and resulting power estimates.
