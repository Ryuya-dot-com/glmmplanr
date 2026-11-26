# Presentation checks reuse saved results; the small scenario plan is a display fixture.
pkgload::load_all(quiet = TRUE)
library(ggplot2)
out <- "validation/plot_defaults"
dir.create(out, showWarnings = FALSE)
read_saved <- function(name) readRDS(file.path("validation/release_workflows", paste0(name, ".rds")))
plan <- read_saved("glmm_partial")
validation <- read_saved("confirmation_partial")
fixed <- read_saved("rasch_fixed")
random <- read_saved("rasch_random")
before <- serialize(list(plan, validation, fixed, random), NULL)

full <- design_crossed(list(method = c("restudy", "retrieval")))
partial <- design_crossed(list(method = c("restudy", "retrieval")), items_per_subject = 4L)
calibrate <- function(design) calibrate_dgp(design, c(restudy = .42, retrieval = .55),
  subj_sd = .3, item_sd = .4)
dgps <- list(full = calibrate(full), partial = calibrate(partial))
spec <- specify_target(dgps$full, "method", scale = "pp", prob_se_method = "joint_delta")
scenarios <- optimize_sample_size(spec, 16L, c(8L, 10L), scenarios = dgps,
  standard_scenario = "full", method = "grid", nsim_final = 4L, workers = 1L,
  verbose = FALSE, control = optimize_control(nsim_init = 4L, n_verify = 2L, nsim_sensitivity = 0L))

plots <- list(
  contour = function(show) plot(plan, show_annotations = show),
  pareto = function(show) plot(plan, type = "pareto", show_annotations = show),
  scenarios = function(show) plot(scenarios, type = "scenarios", show_annotations = show),
  rates = function(show) plot(validation, show_annotations = show),
  patterns = function(show) plot(validation, type = "patterns", show_annotations = show),
  rasch_fixed = function(show) plot(fixed, show_annotations = show),
  rasch_rmse = function(show) plot(fixed, type = "items", show_annotations = show),
  rasch_bias = function(show) plot(fixed, type = "items", metric = "bias", show_annotations = show),
  rasch_random = function(show) plot(random, show_annotations = show),
  dgp = function(show) check_dgp(plan$spec$dgp, n_draws = 1000L, show_annotations = show)$plot
)
for (name in names(plots)) {
  plain <- plots[[name]](FALSE); annotated <- plots[[name]](TRUE)
  stopifnot(all(vapply(c("title", "subtitle", "caption"), function(k) is.null(plain$labels[[k]]), logical(1))),
    !is.null(annotated$labels$title),
    identical(plain$theme$panel.background$fill, "white"),
    !inherits(plain$theme$panel.border, "element_blank"),
    identical(ggplot_build(plain)$data, ggplot_build(annotated)$data),
    identical((plain + labs(title = "Custom title"))$labels$title, "Custom title"))
  ggsave(file.path(out, paste0(name, ".png")), plain, width = 9,
    height = if (name == "patterns") 8 else 6, dpi = 150)
  cat(name, ": white panel/border, annotations optional, plotted values unchanged\n")
}
if (requireNamespace("plotly", quietly = TRUE)) {
  plain <- plotly::plotly_build(plot(plan, type = "surface3d"))
  annotated <- plotly::plotly_build(plot(plan, type = "surface3d", show_annotations = TRUE))
  stopifnot(identical(plain$x$layout$title$text, ""), nzchar(annotated$x$layout$title),
    identical(plain$x$layout$paper_bgcolor, "white"),
    identical(plain$x$layout$scene$bgcolor, "white"),
    identical(plain$x$data, annotated$x$data))
  dput(plain$x$layout[c("title", "paper_bgcolor", "scene")], file.path(out, "surface3d-layout.R"))
  cat("surface3d: white background, optional title, trace values unchanged\n")
}
stopifnot(identical(before, serialize(list(plan, validation, fixed, random), NULL)))
capture.output(sessionInfo(), file = file.path(out, "environment.txt"))
cat("Saved plan and validation objects were not modified.\n")
