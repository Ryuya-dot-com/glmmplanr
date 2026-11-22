#' @keywords internal
#' @section Inspecting objects:
#' `print()` displays a design, calibrated generator, generator check, target
#' specification, or sample-size plan. For a plan, `summary()` also shows
#' response-pattern frequencies, invalid fits, and the largest observed
#' counts for extreme participant and item responses at the selected design.
#' @aliases print.glmmplanr_design print.glmmplanr_dgp print.glmmplanr_dgp_check print.glmmplanr_spec print.glmmplanr_result summary.glmmplanr_result
#' @exportS3Method print glmmplanr_design
#' @exportS3Method print glmmplanr_dgp
#' @exportS3Method print glmmplanr_dgp_check
#' @exportS3Method print glmmplanr_spec
#' @exportS3Method print glmmplanr_result
#' @exportS3Method summary glmmplanr_result
#' @useDynLib glmmplanr, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom rlang .data
"_PACKAGE"

utils::globalVariables("p")
