#' Count Observed and Correct Responses by Participant and Item
#'
#' @param y A numeric or logical vector of binary responses (1 is correct).
#'   Missing responses are omitted from the observed count.
#' @param participant,item Participant and item identifiers, one per response.
#'   Supply factors with all planned levels to include units with no trial rows.
#'   The `subj` and `item` factors from [simulate_dgp_data()] retain these levels.
#' @return One row per unit with `unit`, `id`, `n_observed`, `n_correct`, and
#'   `pattern` (all_correct, all_incorrect, mixed, or unobserved).
#' @examples
#' response_patterns(c(1, 1, 0, NA), factor(c("A", "A", "B", "C")),
#'                   c("easy", "other", "hard", "hard"))
#' @export
response_patterns <- function(y, participant, item) {
  if (!(is.numeric(y) || is.logical(y)) || !is.null(dim(y)) ||
      any(!is.na(y) & !(y %in% 0:1))) {
    cli::cli_abort("{.arg y} must be a vector of 0, 1 or missing responses.")
  }
  observed <- !is.na(y)
  by_unit <- function(id, unit) {
    if (!(is.character(id) || is.numeric(id) || is.factor(id)) || !is.null(dim(id)) ||
        length(id) != length(y) || anyNA(id) || anyNA(levels(id))) {
      cli::cli_abort("Supply one non-missing {unit} identifier per response.")
    }
    id <- as.factor(id)
    counts <- unit_response_counts(as.integer(id)[observed], y[observed] == 1L, nlevels(id))
    pattern <- rep("mixed", nlevels(id))
    pattern[counts$all_correct] <- "all_correct"
    pattern[counts$all_incorrect] <- "all_incorrect"
    pattern[counts$unobserved] <- "unobserved"
    data.frame(unit = rep(unit, nlevels(id)), id = levels(id), n_observed = counts$observed,
               n_correct = counts$correct, pattern = pattern)
  }
  rbind(by_unit(participant, "participant"), by_unit(item, "item"))
}

unit_response_counts <- function(id, success, n) {
  observed <- tabulate(id, nbins = n)
  correct <- tabulate(id[success], nbins = n)
  list(observed = observed, correct = correct, all_correct = observed > 0L & correct == observed,
       all_incorrect = observed > 0L & correct == 0L, unobserved = observed == 0L)
}

unit_pattern_names <- c("subj_all_correct", "subj_all_incorrect", "subj_unobserved",
                        "item_all_correct", "item_all_incorrect", "item_unobserved")
observed_max_columns <- paste0("max_obs_", unit_pattern_names[c(1L, 2L, 4L, 5L)])
pattern_columns <- c("pattern_nsim", "pattern_n_invalid", paste0("n_rep_", unit_pattern_names),
                     paste0("n_invalid_", unit_pattern_names), observed_max_columns)
empty_pattern_totals <- stats::setNames(integer(length(pattern_columns)), pattern_columns)

replication_patterns <- function(y, subj, item, n_subj, n_item, fit, scale) {
  success <- y == 1L
  s <- unit_response_counts(subj, success, n_subj)
  i <- unit_response_counts(item, success, n_item)
  present <- c(any(s$all_correct), any(s$all_incorrect), any(s$unobserved),
               any(i$all_correct), any(i$all_incorrect), any(i$unobserved))
  maxima <- c(max(0L, s$observed[s$all_correct]), max(0L, s$observed[s$all_incorrect]),
              max(0L, i$observed[i$all_correct]), max(0L, i$observed[i$all_incorrect]))
  invalid <- !valid_target_fit(fit$est, fit$se, isTRUE(fit$converged), fit$status,
                               fit$prob_est, fit$prob_se, scale)
  stats::setNames(c(1L, as.integer(invalid), present, present & invalid, maxima), pattern_columns)
}

combine_pattern_totals <- function(a, b) {
  totals <- a + b
  totals[observed_max_columns] <- pmax(a[observed_max_columns], b[observed_max_columns])
  totals
}

# Old saved plans have no pattern measurements for their original replications.
pooled_pattern_totals <- function(old, extra) {
  previous <- empty_pattern_totals
  if (isTRUE(old$pattern_nsim > 0L)) previous <- unlist(old[pattern_columns], use.names = TRUE)
  combine_pattern_totals(previous, unlist(extra[pattern_columns], use.names = TRUE))
}
