pkgname <- "glmmplanr"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
base::assign(".ExTimings", "glmmplanr-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('glmmplanr')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("analysis_model")
### * analysis_model

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: analysis_model
### Title: Use the planned model structure in an analysis of observed
###   responses
### Aliases: analysis_model formula.glmmplanr_spec family.glmmplanr_spec
###   coef.glmmplanr_dgp

### ** Examples

d <- design_crossed(list(learning = c("incidental", "instruction")))
g <- calibrate_dgp(d, c(.60, .70), chance = 1/3, subj_sd = .8, item_sd = .6)
s <- specify_target(g, "learning", analysis_re = "intercept_only")
formula(s)
family(s)
coef(g)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("analysis_model", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("as_gt")
### * as_gt

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: as_gt.glmmplanr_rasch
### Title: Render Sample Size Optimization or Stress-Test Results as a
###   Publication-Ready 'gt' Table
### Aliases: as_gt.glmmplanr_rasch as_gt as_gt.glmmplanr_result
###   as_gt.glmmplanr_stress_test

### ** Examples

# Create the five-replication demonstration plan from the optimizer help.
example("optimize_sample_size", package = "glmmplanr", local = environment(), echo = FALSE)
table <- as_gt(plan)
class(table)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("as_gt", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("calibrate_control")
### * calibrate_control

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: calibrate_control
### Title: Control Settings for Calibration, Target Specification, and
###   Optimization
### Aliases: calibrate_control target_control optimize_control

### ** Examples

calibrate_control(testlet_sd = 0.45, participant_testlet_sd = 0.30)
target_control(failure_policy = "effective", prob_boot = 99L)
optimize_control(nsim_init = 40L, nsim_sensitivity = 0L)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("calibrate_control", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("calibrate_dgp")
### * calibrate_dgp

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: calibrate_dgp
### Title: Calibrate Conditional GLMM Parameters from Marginal
###   Probabilities and m-AFC Guessing
### Aliases: calibrate_dgp

### ** Examples

des <- design_crossed(
  factors = list(
    exposure = c("incidental", "intentional"),
    spacing  = c("massed", "distributed")
  ),
  allocation = "latin_square"
)
dgp <- calibrate_dgp(
  des,
  cell_probs = c(
    incidental_massed       = 0.50,
    incidental_distributed  = 0.60,
    intentional_massed      = 0.65,
    intentional_distributed = 0.80
  ),
  chance          = 1 / 3,
  subj_prob_range = c(0.40, 0.85),
  item_prob_range = c(0.45, 0.88),
  slope_sd_ratio  = 0.35,
  corr_structure  = 0.20
)
print(dgp)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("calibrate_dgp", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("cell_contrast_weights")
### * cell_contrast_weights

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: cell_contrast_weights
### Title: Extract Cell-Probability Contrast Weights for a Fixed Effect
### Aliases: cell_contrast_weights

### ** Examples

des <- design_crossed(list(condition = c("a", "b")))
cell_contrast_weights(des, "condition")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cell_contrast_weights", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("cell_probs_from_diffs")
### * cell_probs_from_diffs

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: cell_probs_from_diffs
### Title: Construct Cell Probabilities Directly from Percentage-Point
###   Effects
### Aliases: cell_probs_from_diffs

### ** Examples

des <- design_crossed(list(
  method = c("control", "treatment")
), allocation = "between_subj")
cell_probs_from_diffs(des, baseline_prob = 0.65, effects_pp = c(method = 8))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cell_probs_from_diffs", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("check_dgp")
### * check_dgp

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: check_dgp
### Title: Inspect and Plot Prior-Predictive Diagnostics of a Calibrated
###   DGP
### Aliases: check_dgp

### ** Examples

des <- design_crossed(list(condition = c("a", "b")))
dgp <- calibrate_dgp(des, c(a = 0.4, b = 0.6), subj_sd = 0.3, item_sd = 0.4)
check_dgp(dgp, plot = FALSE)$cell_diagnostics[c("cell", "quad_prob", "conditional_prob")]



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("check_dgp", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("compute_scheduled_trials")
### * compute_scheduled_trials

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: compute_scheduled_trials
### Title: Compute Exact Scheduled Trials and Post-Attrition Expected
###   Observations
### Aliases: compute_scheduled_trials compute_expected_trials

### ** Examples

des <- design_crossed(list(), attrition = 0.1)
compute_scheduled_trials(des, 20, 10)
compute_expected_trials(des, 20, 10)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("compute_scheduled_trials", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("design_crossed")
### * design_crossed

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: design_crossed
### Title: Create a Balanced Crossed, Between-Within, or Testlet
###   Experimental Design
### Aliases: design_crossed

### ** Examples

des <- design_crossed(
  factors = list(
    exposure = c("incidental", "intentional"),
    spacing  = c("massed", "distributed")
  ),
  allocation = "latin_square",
  contrast   = "sum",
  attrition  = 0.05
)
print(des)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("design_crossed", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("design_testlet")
### * design_testlet

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: design_testlet
### Title: Convenience Constructor for Passage / Testlet-Nested Crossed
###   Designs
### Aliases: design_testlet

### ** Examples

des <- design_testlet(list(condition = c("a", "b")), testlet_size = 4L)
des$k_item



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("design_testlet", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("evaluate_design_point")
### * evaluate_design_point

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: evaluate_design_point
### Title: Evaluate Monte Carlo Replications at a Single (n_subj, n_item)
###   Design Point
### Aliases: evaluate_design_point

### ** Examples

des <- design_crossed(list(condition = c("a", "b")))
dgp <- calibrate_dgp(des, c(a = 0.4, b = 0.6), subj_sd = 0.3, item_sd = 0.4)
spec <- specify_target(dgp, "condition", criterion = "power", engine = "bayes_laplace")
# Small run to inspect the output.
r <- evaluate_design_point(spec, 16, 8, nsim = 5L)
r[c("nsim", "n_valid", "power")]



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("evaluate_design_point", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("generate_design_schedule")
### * generate_design_schedule

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: generate_design_schedule
### Title: Generate a Balanced Experimental Trial Schedule for (n_subj,
###   n_item)
### Aliases: generate_design_schedule

### ** Examples

des <- design_crossed(list(condition = c("a", "b")))
schedule <- generate_design_schedule(des, n_subj = 4, n_item = 4)
head(schedule$trials)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("generate_design_schedule", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mafc_binomial")
### * mafc_binomial

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mafc_binomial
### Title: Custom GLM Binomial Family with Lower Guessing Asymptote (m-AFC)
### Aliases: mafc_binomial

### ** Examples

family <- mafc_binomial(chance = 1/3)
family$linkinv(0)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mafc_binomial", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mafc_qlogis")
### * mafc_qlogis

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mafc_qlogis
### Title: Modified Logit Link and Inverse Link for m-AFC Tasks
### Aliases: mafc_qlogis mafc_plogis

### ** Examples

eta <- mafc_qlogis(0.7, chance = 1/3)
mafc_plogis(eta, chance = 1/3)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mafc_qlogis", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("optimize_sample_size")
### * optimize_sample_size

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: optimize_sample_size
### Title: Optimize 2D Sample Size (Subjects x Items) via Monotone GP
###   Surrogate & Pareto Frontier
### Aliases: optimize_sample_size

### ** Examples

des <- design_crossed(list(condition = c("a", "b")))
dgp <- calibrate_dgp(des, c(a = 0.4, b = 0.6), subj_sd = 0.3, item_sd = 0.4)
spec <- specify_target(dgp, "condition", criterion = "power")
# Five replications demonstrate the output. Use more to make a design decision.
plan <- optimize_sample_size(
  spec, c(16, 18), 8, method = "grid", nsim_final = 5L,
  workers = 1L, verbose = FALSE,
  control = optimize_control(nsim_init = 5L, n_verify = 2L, nsim_sensitivity = 0L)
)
plan$verified[c("n_subj", "n_item", "nsim", "decision_status")]



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("optimize_sample_size", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plan_rasch_fixed_bank")
### * plan_rasch_fixed_bank

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plan_rasch_fixed_bank
### Title: Plan Examinee Sample Sizes for Calibrating a Fixed Rasch Item
###   Bank
### Aliases: plan_rasch_fixed_bank

### ** Examples

bank <- seq(-2.0, 2.0, length.out = 12)
fb <- plan_rasch_fixed_bank(
  item_difficulties = bank,
  n_subj = c(40L, 100L),
  max_error = c(difficulty = 0.35, difficulty_max_abs_error = 0.60),
  method = "fixed_item",
  nsim = 10L
)
summary(fb)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plan_rasch_fixed_bank", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plan_rasch_precision")
### * plan_rasch_precision

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plan_rasch_precision
### Title: Plan Sample Sizes for Random-Item Rasch Parameter Recovery
### Aliases: plan_rasch_precision

### ** Examples

design <- design_crossed(list(), allocation = "crossed")
dgp <- calibrate_dgp(design, 0.65, subj_sd = 0.8, item_sd = 0.6)
plan <- plan_rasch_precision(dgp, data.frame(n_subj = 48, n_item = 12),
                             c(difficulty = 0.30), nsim = 4)
as.data.frame(plan)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plan_rasch_precision", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot.glmmplanr_result")
### * plot.glmmplanr_result

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot.glmmplanr_result
### Title: Plot Sample Size Trade-offs and Scenario Comparisons
### Aliases: plot.glmmplanr_result

### ** Examples

# Create the five-replication demonstration plan from the optimizer help.
example("optimize_sample_size", package = "glmmplanr", local = environment(), echo = FALSE)
plot(plan, type = "pareto")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot.glmmplanr_result", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("print.glmmplanr_design")
### * print.glmmplanr_design

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: print.glmmplanr_dgp_check
### Title: Print or Summarize a Planning Object
### Aliases: print.glmmplanr_dgp_check print.glmmplanr_dgp
###   print.glmmplanr_design print.glmmplanr_result
###   summary.glmmplanr_result print.glmmplanr_spec

### ** Examples

des <- design_crossed(list(condition = c("a", "b")))
print(des)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("print.glmmplanr_design", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("prob_range_to_sd")
### * prob_range_to_sd

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: prob_range_to_sd
### Title: Convert Intuitive Accuracy Range to Random-Effect Standard
###   Deviation
### Aliases: prob_range_to_sd

### ** Examples

prob_range_to_sd(c(0.40, 0.85), chance = 1 / 3)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("prob_range_to_sd", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("rasch_plan_methods")
### * rasch_plan_methods

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: rasch_plan_methods
### Title: Inspect and plot Rasch planning results
### Aliases: rasch_plan_methods as.data.frame.glmmplanr_rasch
###   summary.glmmplanr_rasch print.glmmplanr_rasch plot.glmmplanr_rasch

### ** Examples

bank <- plan_rasch_fixed_bank(c(-1, 0, 1), 20L, c(difficulty = .5), nsim = 4L)
subset(as.data.frame(bank), target == "joint")
plot(bank, type = "items", metric = "bias")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("rasch_plan_methods", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("report_plan")
### * report_plan

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: report_plan
### Title: Write a Sample Size Justification
### Aliases: report_plan

### ** Examples

# Create the five-replication demonstration plan from the optimizer help.
example("optimize_sample_size", package = "glmmplanr", local = environment(), echo = FALSE)
report_plan(plan, print_output = FALSE)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("report_plan", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("response_patterns")
### * response_patterns

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: response_patterns
### Title: Count Observed and Correct Responses by Participant and Item
### Aliases: response_patterns

### ** Examples

response_patterns(c(1, 1, 0, NA), factor(c("A", "A", "B", "C")),
                  c("easy", "other", "hard", "hard"))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("response_patterns", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("scenario_dgps")
### * scenario_dgps

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: scenario_dgps
### Title: Build Multi-Scenario DGPs Varying Percentage-Point Effect Sizes
###   and/or Random Slopes
### Aliases: scenario_dgps

### ** Examples

des <- design_crossed(list(condition = c("a", "b")))
dgp <- calibrate_dgp(des, c(a = 0.4, b = 0.6), subj_sd = 0.3, item_sd = 0.4)
models <- scenario_dgps(dgp, "condition", c(smaller = 5, larger = 10))
lapply(models, function(x) x$target_cell_probs)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("scenario_dgps", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("sd_to_prob_range")
### * sd_to_prob_range

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: sd_to_prob_range
### Title: Convert Random-Effect SD Back to Conditional Probability
###   Interval
### Aliases: sd_to_prob_range

### ** Examples

sd_to_prob_range(sd = 0.5, intercept = qlogis(0.7))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("sd_to_prob_range", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("sensitivity_dgps")
### * sensitivity_dgps

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: sensitivity_dgps
### Title: Create Sensitivity Variants of a Calibrated DGP (Optimistic /
###   Standard / Conservative)
### Aliases: sensitivity_dgps

### ** Examples

des <- design_crossed(list(condition = c("a", "b")))
dgp <- calibrate_dgp(des, c(a = 0.4, b = 0.6), subj_sd = 0.3, item_sd = 0.4,
                     slope_sd_ratio = 0.35)
models <- sensitivity_dgps(dgp)
vapply(models, function(x) unname(x$tau_subj[2L]), numeric(1))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("sensitivity_dgps", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("simulate_dgp_data")
### * simulate_dgp_data

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: simulate_dgp_data
### Title: Simulate One Balanced Dataset from a Calibrated DGP
### Aliases: simulate_dgp_data

### ** Examples

des <- design_crossed(list(condition = c("a", "b")))
dgp <- calibrate_dgp(des, c(a = 0.4, b = 0.6), subj_sd = 0.3, item_sd = 0.4)
trials <- simulate_dgp_data(dgp, n_subj = 8, n_item = 4)
head(trials[c("subj", "item", "condition", "y")])



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("simulate_dgp_data", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("specify_target")
### * specify_target

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: specify_target
### Title: Specify Target Effect, Scale, Evaluation Criteria, and
###   Estimation Strategy
### Aliases: specify_target

### ** Examples

des <- design_crossed(list(condition = c("a", "b")))
dgp <- calibrate_dgp(des, c(a = 0.4, b = 0.6), subj_sd = 0.3, item_sd = 0.4)
spec <- specify_target(dgp, "condition", scale = "pp", criterion = "ci_width", ci_width = 12)
spec$true_effect_target



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("specify_target", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("stress_test_models")
### * stress_test_models

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: stress_test_models
### Title: Stress-Test Analysis Models Under Maximal DGP Mismatch
### Aliases: stress_test_models

### ** Examples

des <- design_crossed(list(condition = c("a", "b")))
dgp <- calibrate_dgp(des, c(a = 0.4, b = 0.6), subj_sd = 0.3, item_sd = 0.4)
comparison <- stress_test_models(dgp, "condition", 16, 8, nsim = 3L)
head(comparison)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("stress_test_models", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("switch_preference")
### * switch_preference

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: switch_preference
### Title: Switch the Preferred Trade-Off Rule Without Re-Running
###   Simulations
### Aliases: switch_preference

### ** Examples

# Create the five-replication demonstration plan from the optimizer help.
example("optimize_sample_size", package = "glmmplanr", local = environment(), echo = FALSE)
reordered <- switch_preference(plan, "fewer_items", verbose = FALSE)
reordered$optimal[c("n_subj", "n_item", "decision_status")]



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("switch_preference", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("topup_verification")
### * topup_verification

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: topup_verification
### Title: Run Additional Independent Top-Up Replications on Borderline
###   Candidates
### Aliases: topup_verification

### ** Examples

# Create the five-replication demonstration plan from the optimizer help.
example("optimize_sample_size", package = "glmmplanr", local = environment(), echo = FALSE)
more <- topup_verification(plan, nsim_extra = 5L, verbose = FALSE)
more$verified[c("n_subj", "n_item", "nsim", "decision_status")]



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("topup_verification", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("true_prob_effect")
### * true_prob_effect

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: true_prob_effect
### Title: Compute True Population-Marginal or Conditional Probability
###   Contrast from DGP
### Aliases: true_prob_effect

### ** Examples

des <- design_crossed(list(condition = c("a", "b")))
dgp <- calibrate_dgp(des, c(a = 0.4, b = 0.6), subj_sd = 0.3, item_sd = 0.4)
w <- cell_contrast_weights(des, "condition")
c(marginal = true_prob_effect(dgp, w), conditional = true_prob_effect(dgp, w, "conditional"))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("true_prob_effect", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("valid_multiples")
### * valid_multiples

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: valid_multiples
### Title: Snap Candidate Sample Sizes to Balanced Multiples of List Count
###   K
### Aliases: valid_multiples

### ** Examples

valid_multiples(c(10, 20), k = 4)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("valid_multiples", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("validate_design_point")
### * validate_design_point

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: validate_design_point
### Title: Validate the Planned Analysis at a Fixed Design
### Aliases: validate_design_point

### ** Examples

des <- design_crossed(list(condition = c("a", "b")))
dgp <- calibrate_dgp(des, c(a = 0.4, b = 0.6), subj_sd = 0.3, item_sd = 0.4)
spec <- specify_target(dgp, "condition", criterion = "power", engine = "bayes_laplace")
r <- validate_design_point(spec, 16, 8, nsim = 5L)
r[c("condition", "nsim", "n_valid", "detection", "coverage")]



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("validate_design_point", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
