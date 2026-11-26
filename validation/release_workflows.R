# Public-API integration checks, not Monte Carlo calibration evidence.
library(glmmplanr)
args <- commandArgs(TRUE)
out <- if (length(args)) args[1] else "validation/release_workflows"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
roundtrip <- function(x, name, type = NULL, report = TRUE) {
  path <- file.path(out, paste0(name, ".rds"))
  saveRDS(x, path)
  y <- readRDS(path)
  stopifnot(identical(as.data.frame(x), as.data.frame(y)))
  write.csv(as.data.frame(y), file.path(out, paste0(name, ".csv")), row.names = FALSE)
  if (report) {
    text <- report_plan(x, print_output = FALSE)
    stopifnot(identical(text, report_plan(y, print_output = FALSE)), inherits(as_gt(y), "gt_tbl"))
    writeLines(text, file.path(out, paste0(name, ".txt")))
  }
  fig <- if (is.null(type)) plot(y) else plot(y, type = type)
  stopifnot(inherits(ggplot2::ggplot_build(fig), "ggplot_built"))
  ggplot2::ggsave(file.path(out, paste0(name, ".png")), fig, width = 8, height = 5, dpi = 120)
  y
}
specs <- list()
for (layout in c("full", "partial")) {
  des <- design_crossed(list(method = c("control", "treatment")),
    items_per_subject = if (layout == "partial") 8L else NULL)
  dgp <- calibrate_dgp(des, c(control = .42, treatment = .55),
    subj_sd = c(.5, .15), item_sd = c(.4, .2), corr_structure = "zero")
  spec <- specify_target(dgp, "method", scale = "pp", criterion = "joint", ci_width = 30)
  stopifnot(spec$engine == "bayes_laplace", spec$prob_se_method == "bootstrap", spec$prob_boot == 99L)
  specs[[layout]] <- spec
  prior <- file.path(out, paste0("glmm_", layout, ".rds"))
  plan <- if (file.exists(prior)) readRDS(prior) else
    optimize_sample_size(spec, c(16L, 24L), c(8L, 12L), method = "grid",
      nsim_final = 6L, workers = 1L, seed = 943260001L, verbose = FALSE,
      control = optimize_control(nsim_sensitivity = 0L))
  stopifnot(isTRUE(all.equal(plan$spec, spec, tolerance = 0)))
  saved <- roundtrip(plan, paste0("glmm_", layout), "contour")
  stopifnot(grepl("does not verify inferential validity", report_plan(saved, FALSE), fixed = TRUE))
  check <- validate_design_point(spec, 24L, 12L, nsim = 6L, seed = 943260002L)
  roundtrip(check, paste0("confirmation_", layout), report = FALSE)
  stopifnot(identical(check$condition, c("null", "planned")), all(check$nsim == 6L))
}
# Reuse the partial-presented default path to check sequential/worker agreement.
sp <- specs$partial
serial <- as.data.frame(readRDS(file.path(out, "glmm_partial.rds")))
parallel <- as.data.frame(optimize_sample_size(sp, c(16L, 24L), c(8L, 12L), method = "grid",
  nsim_final = 6L, workers = 2L, seed = 943260001L, verbose = FALSE,
  control = optimize_control(nsim_sensitivity = 0L)))
metrics <- intersect(c("n_valid", "power", "coverage", "ci_width_mean", "prob_diff_est_pp"), names(serial))
stopifnot(length(metrics) >= 4L, isTRUE(all.equal(serial[metrics], parallel[metrics], tolerance = 1e-10)))

g <- calibrate_dgp(design_crossed(list()), .6, subj_sd = .6, item_sd = .5)
random <- plan_rasch_precision(g, data.frame(n_subj = c(32L, 48L), n_item = 12L),
  c(difficulty = .5), nsim = 6L, seed = 943260004L)
roundtrip(random, "rasch_random")
fixed <- plan_rasch_fixed_bank(setNames(seq(-2, 2, length.out = 12), paste0("item", 1:12)),
  c(32L, 48L), c(difficulty = .5), items_per_subject = 6L, nsim = 6L, seed = 943260005L)
roundtrip(fixed, "rasch_fixed", "items")
capture.output(sessionInfo(), file = file.path(out, "environment.txt"))
file.copy("validation/release_workflows.R", file.path(out, "source_script.R"), overwrite = TRUE)
cat("Public default GLMM, confirmation, two Rasch targets, serialization, plots, reports, and worker checks passed.\n")
