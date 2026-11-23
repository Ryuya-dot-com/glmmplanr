#' Create a Balanced Crossed, Between-Within, or Testlet Experimental Design
#'
#' Declares the factorial, counterbalancing (Latin square), nesting (testlet),
#' contrast coding, and trial-level attrition structure of a linguistic or
#' psycholinguistic experiment without requiring pilot data or pre-fitted models.
#'
#' @param factors Named list mapping factor names to character vectors of at
#'   least 2 level names (ordered from reference/control to focal/treatment).
#'   Use `list()` for a single cell named `overall` with intercepts only.
#' @param allocation Character string specifying the experimental allocation
#'   preset:
#'   * `"latin_square"` (default): Counterbalanced Latin square across crossed
#'     within-subject and within-item factors so each subject sees each item at
#'     most once and list multiplicity is enforced automatically.
#'   * `"between_subj"`: All factors (unless overridden by `between_item`) are
#'     treated as between-subject groups completing a common set of items.
#'   * `"between_item"`: All factors (unless overridden by `between_subj`)
#'     partition items while subjects see all items.
#'   * `"crossed"`: Full factorial within-subject and within-item repetition
#'     (every subject sees every item under every within-unit cell).
#'   * `"testlet"`: Same as `"latin_square"` with passages/testlets grouping
#'     multiple items (requires `testlet_size >= 2`).
#' @param between_subj Optional character vector of factor names varied between
#'   participants (e.g., `"group"` for L1 vs L2 learners). Classroom-level
#'   assignment would require an additional clustering level.
#' @param between_item Optional character vector of factor names varied between
#'   items (e.g., `"word_freq"` for low- vs high-frequency target words).
#' @param testlet_size Optional integer \eqn{\ge 2} specifying the number of
#'   questions nested within each reading/listening passage (testlet).
#' @param contrast Contrast coding scheme automatically attached to all design
#'   factors: `"sum"` (default, centered indicators with `(-0.5, +0.5)` for two
#'   levels),
#'   `"helmert"` (unit-scaled orthogonal Helmert coding), or `"treatment"`.
#' @param attrition Numeric scalar in `[0, 0.9)` giving the independent
#'   trial missingness rate among retained participants (default `0`).
#' @param participant_dropout Probability in `[0, 1)` that a recruited participant
#'   provides no trials, independently of outcomes and trial missingness.
#'
#' @param items_per_subject Number of distinct items drawn without replacement for
#'   each participant from the bank of `n_item` items. NULL uses the entire bank.
#'   Draws are independent across participants, so item and condition exposure
#'   can vary. Nested designs draw complete testlets. Crossed designs present
#'   all within-pair conditions for every selected item.
#' @param item_seed Integer seed for the item allocation (default 1). A design
#'   point uses the same schedule across response replications and engines.
#'
#' @return An S3 object of class `"glmmplanr_design"` containing cell tables,
#'   contrast matrices, design matrices (`X_cells`, `Z_subj_cells`,
#'   `Z_item_cells`), and list-balance multiplicity constants (`k_subj`,
#'   `k_item`).
#'
#' @examples
#' des <- design_crossed(
#'   factors = list(
#'     exposure = c("incidental", "intentional"),
#'     spacing  = c("massed", "distributed")
#'   ),
#'   allocation = "latin_square",
#'   contrast   = "sum",
#'   attrition  = 0.05
#' )
#' print(des)
#' @export
design_crossed <- function(
  factors,
  allocation = c("latin_square", "between_subj", "between_item", "crossed", "testlet"),
  between_subj = NULL,
  between_item = NULL,
  testlet_size = NULL,
  contrast = c("sum", "helmert", "treatment"),
  attrition = 0,
  participant_dropout = 0,
  items_per_subject = NULL,
  item_seed = 1L
) {
  allocation <- match.arg(allocation)
  contrast <- match.arg(contrast)

  if (!is.list(factors) ||
      (length(factors) && (is.null(names(factors)) || anyNA(names(factors)) ||
                          any(!nzchar(names(factors)))))) {
    cli::cli_abort("{.arg factors} must be a named list of factor levels, or an empty list for an intercept-only design.")
  }
  factor_names <- as.character(names(factors))
  if (anyDuplicated(factor_names)) {
    cli::cli_abort("Factor names in {.arg factors} must be unique.")
  }
  if (any(factor_names %in% c("subj", "item", "testlet", "list_id", "cell_id", "y"))) {
    cli::cli_abort(
      "Factor names cannot use reserved column names ({.val subj}, {.val item}, {.val testlet}, {.val list_id}, {.val cell_id}, {.val y})."
    )
  }
  for (fn in factor_names) {
    lv <- as.character(factors[[fn]])
    if (length(lv) < 2L || anyNA(lv) || any(!nzchar(lv)) || anyDuplicated(lv)) {
      cli::cli_abort("Factor {.field {fn}} must have at least 2 unique non-empty level names.")
    }
    factors[[fn]] <- lv
  }

  if (!is.numeric(attrition) || length(attrition) != 1L || is.na(attrition) ||
      attrition < 0 || attrition >= 0.9) {
    cli::cli_abort("{.arg attrition} must be a numeric scalar in [0, 0.9).")
  }

  if (length(participant_dropout) != 1L || !is.numeric(participant_dropout) ||
      !is.finite(participant_dropout) || participant_dropout < 0 || participant_dropout >= 1) {
    cli::cli_abort("{.arg participant_dropout} must be a probability in [0, 1).")
  }

  # Preset handling for allocation
  if (allocation == "between_subj" && is.null(between_subj)) {
    between_subj <- setdiff(factor_names, as.character(between_item))
  }
  if (allocation == "between_item" && is.null(between_item)) {
    between_item <- setdiff(factor_names, as.character(between_subj))
  }
  if (allocation == "testlet" && is.null(testlet_size)) {
    testlet_size <- 4L
  }

  between_subj <- if (is.null(between_subj)) character(0) else as.character(between_subj)
  between_item <- if (is.null(between_item)) character(0) else as.character(between_item)

  overlap <- intersect(between_subj, between_item)
  if (length(overlap)) {
    cli::cli_abort("A factor cannot be both between-subject and between-item: {.val {overlap}}.")
  }
  unknown_bs <- setdiff(between_subj, factor_names)
  unknown_bi <- setdiff(between_item, factor_names)
  if (length(unknown_bs) > 0L) {
    cli::cli_abort("Unknown factor(s) in {.arg between_subj}: {.val {unknown_bs}}.")
  }
  if (length(unknown_bi) > 0L) {
    cli::cli_abort("Unknown factor(s) in {.arg between_item}: {.val {unknown_bi}}.")
  }

  has_testlet <- !is.null(testlet_size)
  if (has_testlet) {
    testlet_size <- as.integer(testlet_size)
    if (length(testlet_size) != 1L || is.na(testlet_size) || testlet_size < 2L) {
      cli::cli_abort("{.arg testlet_size} must be an integer >= 2 when specified.")
    }
  } else {
    testlet_size <- 1L
  }

  if (!is.null(items_per_subject)) {
    if (length(items_per_subject) != 1L || !is.numeric(items_per_subject) ||
        !is.finite(items_per_subject) || items_per_subject < 1 ||
        items_per_subject > .Machine$integer.max || items_per_subject %% testlet_size != 0) {
      cli::cli_abort("{.arg items_per_subject} must be a positive integer multiple of the testlet size ({testlet_size}).")
    }
    if (length(item_seed) != 1L || !is.numeric(item_seed) || !is.finite(item_seed) ||
        abs(item_seed) > .Machine$integer.max || item_seed != trunc(item_seed)) {
      cli::cli_abort("{.arg item_seed} must be a single finite integer seed.")
    }
    items_per_subject <- as.integer(items_per_subject)
  }

  within_subj <- setdiff(factor_names, between_subj)
  within_item <- setdiff(factor_names, between_item)
  ls_factors <- intersect(within_subj, within_item)

  # Build canonical cell grid in outer-to-inner factor order so
  # cell 1..M follows natural reading order (e.g., incidental_massed,
  # incidental_distributed, intentional_massed, intentional_distributed).
  cells <- subfactor_grid(factors, factor_names)

  contr_list <- list()
  for (fn in factor_names) {
    lv <- factors[[fn]]
    f_col <- factor(cells[[fn]], levels = lv)
    cm <- build_contrast_matrix(lv, scheme = contrast)
    stats::contrasts(f_col) <- cm
    cells[[fn]] <- f_col
    contr_list[[fn]] <- cm
  }

  cells$cell_name <- if (length(factor_names)) {
    do.call(paste, c(unname(as.list(cells)), sep = "_"))
  } else {
    "overall"
  }
  cells$cell_id <- seq_len(nrow(cells))
  rownames(cells) <- cells$cell_name

  # Build fixed and random design matrices on the unique cells
  fixed_rhs <- if (length(factor_names)) {
    stats::reformulate(paste(factor_names, collapse = " * "))
  } else {
    ~ 1
  }
  X_cells <- stats::model.matrix(fixed_rhs, data = cells)
  effect_map <- build_effect_map(fixed_rhs, cells, X_cells)

  subj_rhs <- if (length(within_subj) > 0L) {
    stats::reformulate(paste(within_subj, collapse = " * "))
  } else {
    ~ 1
  }
  item_rhs <- if (length(within_item) > 0L) {
    stats::reformulate(paste(within_item, collapse = " * "))
  } else {
    ~ 1
  }
  Z_subj_cells <- stats::model.matrix(subj_rhs, data = cells)
  Z_item_cells <- stats::model.matrix(item_rhs, data = cells)

  if (has_testlet) {
    Z_testlet_cells <- Z_item_cells
    testlet_rhs <- item_rhs
  } else {
    Z_testlet_cells <- NULL
    testlet_rhs <- NULL
  }

  # Compute list-balance multiplicity constants
  n_bs_cells <- if (length(between_subj) > 0L) {
    prod(vapply(factors[between_subj], length, integer(1L)))
  } else {
    1L
  }
  n_bi_cells <- if (length(between_item) > 0L) {
    prod(vapply(factors[between_item], length, integer(1L)))
  } else {
    1L
  }
  n_ls_lists <- if (allocation != "crossed" && length(ls_factors) > 0L) {
    prod(vapply(factors[ls_factors], length, integer(1L)))
  } else {
    1L
  }
  trials_per_pair <- if (allocation == "crossed" && length(ls_factors) > 0L) {
    as.integer(prod(vapply(factors[ls_factors], length, integer(1L))))
  } else {
    1L
  }

  k_subj <- as.integer(n_bs_cells * n_ls_lists)
  k_item <- as.integer(n_bi_cells * n_ls_lists * testlet_size)

  structure(
    list(
      factors = factors,
      factor_names = factor_names,
      allocation = allocation,
      between_subj = between_subj,
      between_item = between_item,
      within_subj = within_subj,
      within_item = within_item,
      ls_factors = ls_factors,
      has_testlet = has_testlet,
      testlet_size = testlet_size,
      contrast = contrast,
      contrasts = contr_list,
      attrition = attrition,
      participant_dropout = participant_dropout,
      items_per_subject = items_per_subject,
      item_seed = if (!is.null(items_per_subject)) as.integer(item_seed),
      cells = cells,
      n_cells = nrow(cells),
      fixed_rhs = fixed_rhs,
      subj_rhs = subj_rhs,
      item_rhs = item_rhs,
      testlet_rhs = testlet_rhs,
      X_cells = X_cells,
      Z_subj_cells = Z_subj_cells,
      Z_item_cells = Z_item_cells,
      Z_testlet_cells = Z_testlet_cells,
      effect_map = effect_map,
      n_ls_lists = as.integer(n_ls_lists),
      n_bs_cells = as.integer(n_bs_cells),
      n_bi_cells = as.integer(n_bi_cells),
      trials_per_pair = trials_per_pair,
      k_subj = k_subj,
      k_item = k_item
    ),
    class = "glmmplanr_design"
  )
}

#' Convenience Constructor for Passage / Testlet-Nested Crossed Designs
#'
#' @inheritParams design_crossed
#' @param testlet_size Integer \eqn{\ge 2} giving the number of items nested per
#'   passage/testlet (default `4`).
#' @return An object of class `"glmmplanr_design"`.
#' @examples
#' des <- design_testlet(list(condition = c("a", "b")), testlet_size = 4L)
#' des$k_item
#' @export
design_testlet <- function(
  factors,
  testlet_size = 4L,
  between_subj = NULL,
  between_item = NULL,
  contrast = c("sum", "helmert", "treatment"),
  attrition = 0,
  participant_dropout = 0,
  items_per_subject = NULL,
  item_seed = 1L
) {
  design_crossed(
    factors = factors,
    allocation = "testlet",
    between_subj = between_subj,
    between_item = between_item,
    testlet_size = testlet_size,
    contrast = contrast,
    attrition = attrition,
    participant_dropout = participant_dropout,
    items_per_subject = items_per_subject,
    item_seed = item_seed
  )
}

#' Build a Factor Contrast Matrix
#' @noRd
build_contrast_matrix <- function(levels_chr, scheme = c("sum", "helmert", "treatment")) {
  scheme <- match.arg(scheme)
  L <- length(levels_chr)
  if (scheme == "sum") {
    # Centered deviation coding: columns sum to 0, each column j compares
    # level (j + 1) vs level 1 with unit step (-0.5, +0.5 for L = 2).
    cm <- matrix(-1 / L, nrow = L, ncol = L - 1L)
    for (j in seq_len(L - 1L)) {
      cm[j + 1L, j] <- (L - 1L) / L
    }
    colnames(cm) <- as.character(seq_len(L - 1L))
  } else if (scheme == "helmert") {
    # Unit-scaled orthogonal Helmert coding (Schad et al., 2020)
    cm <- stats::contr.helmert(L)
    for (j in seq_len(L - 1L)) {
      cm[, j] <- cm[, j] / (j + 1L)
    }
    colnames(cm) <- as.character(seq_len(L - 1L))
  } else {
    cm <- stats::contr.treatment(L)
  }
  rownames(cm) <- levels_chr
  cm
}

#' Map Term Labels (e.g. "exposure:spacing") to Model Matrix Columns
#' @noRd
build_effect_map <- function(fixed_rhs, cells, X_cells) {
  assign_idx <- attr(X_cells, "assign")
  term_labels <- attr(stats::terms(fixed_rhs), "term.labels")
  col_names <- colnames(X_cells)
  map <- list("(Intercept)" = col_names[assign_idx == 0L])
  for (i in seq_along(term_labels)) {
    map[[term_labels[i]]] <- col_names[assign_idx == i]
  }
  map
}

#' Snap Candidate Sample Sizes to Balanced Multiples of List Count K
#'
#' @param n_range Integer scalar or length-2 range `c(min, max)` or explicit vector.
#' @param k Positive integer multiplicity required for exact list balance.
#' @param arg_name Name for error reporting.
#' @return Sorted integer vector of valid multiples of `k` (or full min/max grid).
#' @examples
#' valid_multiples(c(10, 20), k = 4)
#' @export
valid_multiples <- function(n_range, k, arg_name = "n_range") {
  k <- as.integer(k)
  if (!is.numeric(n_range) || length(n_range) < 1L || anyNA(n_range)) {
    cli::cli_abort("{.arg {arg_name}} must be a numeric scalar or range c(min, max).")
  }
  if (length(n_range) == 1L) {
    val <- as.integer(round(n_range / k) * k)
    if (val < k) val <- k
    return(val)
  }
  if (length(n_range) == 2L) {
    lo <- as.integer(ceiling(min(n_range) / k) * k)
    hi <- as.integer(floor(max(n_range) / k) * k)
    if (lo < k) lo <- k
    if (hi < lo) {
      cli::cli_abort(
        "No multiple of list balance factor {.val {k}} exists inside {.arg {arg_name}} = c({min(n_range)}, {max(n_range)})."
      )
    }
    return(seq.int(lo, hi, by = k))
  }
  vals <- sort(unique(as.integer(round(n_range / k) * k)))
  vals[vals >= k]
}

#' Generate a Balanced Experimental Trial Schedule for (n_subj, n_item)
#'
#' Constructs the exact trial-by-trial layout (before optional random trial
#' attrition) enforcing Latin-square list balance, between-unit partitions, and
#' passage/testlet nesting.
#'
#' @param design A `"glmmplanr_design"` object.
#' @param n_subj Number of subjects (automatically checked/snapped to a multiple
#'   of `design$k_subj`).
#' @param n_item Number of items (automatically checked/snapped to a multiple
#'   of `design$k_item`).
#' @param strict Logical; if `TRUE` (default), abort if `n_subj` or `n_item` is
#'   not an exact multiple of `k_subj` / `k_item`.
#' @param diagnostics If TRUE, include `item_exposure` with one row per planned
#'   item. `n_scheduled` counts trials before missingness or dropout. Items share
#'   a `component` if a chain of participants was scheduled to see them in common.
#'   Unscheduled items each form an isolated component.
#' @return A list containing the trial schedule data frame `trials` (with
#'   `subj_id`, `item_id`, `testlet_id`, `list_id`, `cell_id`, and factor
#'   columns) plus row index mappings into `design$X_cells`.
#' @examples
#' des <- design_crossed(list(condition = c("a", "b")))
#' schedule <- generate_design_schedule(des, n_subj = 4, n_item = 4)
#' head(schedule$trials)
#' @export
generate_design_schedule <- function(design, n_subj, n_item, strict = TRUE, diagnostics = FALSE) {
  if (!rlang::is_bool(diagnostics)) cli::cli_abort("{.arg diagnostics} must be TRUE or FALSE.")
  if (!inherits(design, "glmmplanr_design")) {
    cli::cli_abort("{.arg design} must be a {.cls glmmplanr_design} object.")
  }
  n_subj <- as.integer(n_subj)
  n_item <- as.integer(n_item)
  if (n_subj < design$k_subj || n_item < design$k_item) {
    cli::cli_abort(
      "Sample sizes ({n_subj}, {n_item}) must be at least ({design$k_subj}, {design$k_item})."
    )
  }
  if (strict) {
    if (n_subj %% design$k_subj != 0L) {
      cli::cli_abort(
        "{.arg n_subj} ({n_subj}) must be a multiple of list count {.val {design$k_subj}}."
      )
    }
    if (n_item %% design$k_item != 0L) {
      cli::cli_abort(
        "{.arg n_item} ({n_item}) must be a multiple of item-list balance factor {.val {design$k_item}}."
      )
    }
  } else {
    n_subj <- as.integer(round(n_subj / design$k_subj) * design$k_subj)
    n_item <- as.integer(round(n_item / design$k_item) * design$k_item)
  }

  # Partition subjects into between-subject groups and Latin-square lists
  subj_ids <- seq_len(n_subj)
  bs_group_s <- if (design$n_bs_cells > 1L) {
    ((subj_ids - 1L) %% design$n_bs_cells) + 1L
  } else {
    rep.int(1L, n_subj)
  }
  ls_list_s <- if (design$n_ls_lists > 1L) {
    (((subj_ids - 1L) %/% design$n_bs_cells) %% design$n_ls_lists) + 1L
  } else {
    rep.int(1L, n_subj)
  }

  # Partition items (or testlets, with items nested inside testlets)
  n_passages <- n_item %/% design$testlet_size
  item_ids <- seq_len(n_item)
  testlet_of_item <- ((item_ids - 1L) %/% design$testlet_size) + 1L

  bi_group_p <- if (design$n_bi_cells > 1L) {
    ((seq_len(n_passages) - 1L) %% design$n_bi_cells) + 1L
  } else {
    rep.int(1L, n_passages)
  }
  ls_set_p <- if (design$n_ls_lists > 1L) {
    (((seq_len(n_passages) - 1L) %/% design$n_bi_cells) %% design$n_ls_lists) + 1L
  } else {
    rep.int(1L, n_passages)
  }

  bi_group_i <- bi_group_p[testlet_of_item]
  ls_set_i <- ls_set_p[testlet_of_item]

  # Build sub-grids for BS, BI, and LS factor levels
  bs_grid <- subfactor_grid(design$factors, design$between_subj)
  bi_grid <- subfactor_grid(design$factors, design$between_item)
  ls_grid <- subfactor_grid(design$factors, design$ls_factors)

  # Map (bs_idx, bi_idx, ls_idx) -> cell_id in design$cells
  cell_lookup <- build_cell_lookup(design, bs_grid, bi_grid, ls_grid)

  m <- design$items_per_subject
  if (!is.null(m) && m > n_item) {
    cli::cli_abort("{.arg n_item} must be at least {.arg items_per_subject} ({m}).")
  }
  if (!is.null(m) && m < n_item) {
    if (uses_box_muller()) {
      return(isolated_seeded_call("generate_design_schedule", environment(), match.call()))
    }
    # Draw whole passages, with ordinary items as passages of size one.
    size <- design$testlet_size
    passages <- with_rng_stream(design$item_seed,
      replicate(n_subj, sample.int(n_passages, m %/% size)))
    selected <- rep((as.integer(passages) - 1L) * size, each = size) + seq_len(size)
    grid <- data.frame(subj_id = rep(subj_ids, each = m), item_id = selected)
  } else {
    grid <- expand.grid(subj_id = subj_ids, item_id = item_ids, KEEP.OUT.ATTRS = FALSE)
  }
  if (design$allocation == "crossed" && length(design$ls_factors) > 0L) {
    n_pair <- nrow(grid)
    ls_idx <- rep(seq_len(nrow(ls_grid)), each = n_pair)
    grid <- data.frame(subj_id = rep.int(grid$subj_id, nrow(ls_grid)),
                       item_id = rep.int(grid$item_id, nrow(ls_grid)))
    list_id <- rep.int(1L, nrow(grid))
  } else {
    ls_s <- ls_list_s[grid$subj_id]
    ls_i <- ls_set_i[grid$item_id]
    ls_idx <- if (design$n_ls_lists > 1L) {
      ((ls_s - 1L + ls_i - 1L) %% design$n_ls_lists) + 1L
    } else {
      rep.int(1L, nrow(grid))
    }
    list_id <- (bs_group_s[grid$subj_id] - 1L) * design$n_ls_lists + ls_s
  }

  bs_idx <- bs_group_s[grid$subj_id]
  bi_idx <- bi_group_i[grid$item_id]
  cell_id <- cell_lookup[cbind(bs_idx, bi_idx, ls_idx)]
  trials <- data.frame(
    subj = factor(grid$subj_id),
    item = factor(grid$item_id, levels = item_ids),
    subj_id = grid$subj_id,
    item_id = grid$item_id,
    list_id = list_id,
    cell_id = cell_id
  )
  if (design$has_testlet) {
    trials$testlet_id <- testlet_of_item[grid$item_id]
    trials$testlet <- factor(trials$testlet_id, levels = seq_len(n_passages))
  }
  for (fn in design$factor_names) {
    trials[[fn]] <- design$cells[[fn]][cell_id]
  }

  schedule <- list(
    trials = trials,
    n_subj = n_subj,
    n_item = n_item,
    n_testlet = if (design$has_testlet) n_passages else 0L,
    n_trials_full = nrow(trials),
    n_trials_scheduled = nrow(trials),
    n_trials_expected = compute_expected_trials(design, n_subj, n_item)
  )
  if (diagnostics) {
    components <- if (is.null(m) || m == n_item) rep.int(1L, n_item) else {
      item_components(trials, n_item)
    }
    schedule$item_exposure <- data.frame(item_id = item_ids,
      n_scheduled = tabulate(trials$item_id, n_item), component = components)
  }
  schedule
}

item_components <- function(trials, n_item) {
  parent <- seq_len(n_item)
  root <- function(i) {
    while (parent[i] != i) i <- parent[i]
    i
  }
  for (items in split(trials$item_id, trials$subj_id)) {
    roots <- unique(vapply(unique(items), root, integer(1)))
    parent[roots] <- min(roots)
  }
  vapply(seq_len(n_item), root, integer(1))
}

#' Compute Exact Scheduled Trials and Post-Attrition Expected Observations
#'
#' Computes scheduled trial counts directly from the experimental allocation
#' (`N_subj * N_item` for Latin-square / between-unit designs; `N_subj * N_item * K_rep`
#' for full crossed designs). With partial presentation, `items_per_subject`
#' replaces the bank size in both formulas. The design tracks expected post-attrition observations
#' separately. Expected observations are rounded to the nearest whole trial.
#'
#' @param design A `"glmmplanr_design"` object.
#' @param n_subj Integer scalar or vector of subject counts.
#' @param n_item Integer scalar or vector of item counts.
#' @return Integer vector matching `n_subj` and `n_item`.
#' @examples
#' des <- design_crossed(list(), attrition = 0.1)
#' compute_scheduled_trials(des, 20, 10)
#' compute_expected_trials(des, 20, 10)
#' @export
compute_scheduled_trials <- function(design, n_subj, n_item) {
  if (!inherits(design, "glmmplanr_design")) {
    cli::cli_abort("{.arg design} must be a {.cls glmmplanr_design} object.")
  }
  mult <- if (!is.null(design$trials_per_pair)) design$trials_per_pair else 1L
  if (!is.null(design$items_per_subject)) {
    if (any(n_item < design$items_per_subject)) {
      cli::cli_abort("{.arg n_item} must be at least {.arg items_per_subject} ({design$items_per_subject}).")
    }
    n_item <- rep.int(design$items_per_subject, length(n_item))
  }
  as.integer(round(as.numeric(n_subj) * as.numeric(n_item) * as.numeric(mult)))
}

#' @rdname compute_scheduled_trials
#' @export
compute_expected_trials <- function(design, n_subj, n_item) {
  sched <- compute_scheduled_trials(design, n_subj, n_item)
  as.integer(round(as.numeric(sched) * (1 - design$attrition) * (1 - dropout_rate(design))))
}

dropout_rate <- function(design) {
  if (is.null(design$participant_dropout)) 0 else design$participant_dropout
}

observed_trials <- function(design, subj_id, n_subj) {
  dropout <- dropout_rate(design)
  if (dropout == 0 && design$attrition == 0) return(NULL)
  retained <- if (dropout > 0) (stats::runif(n_subj) >= dropout)[subj_id] else TRUE
  if (design$attrition > 0) return(retained & stats::runif(length(subj_id)) >= design$attrition)
  retained
}

subfactor_grid <- function(factors, sub_names) {
  if (length(sub_names) == 0L) {
    return(data.frame(row.names = 1L))
  }
  rev_sub <- rev(factors[sub_names])
  g <- expand.grid(rev_sub, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  g[, sub_names, drop = FALSE]
}

build_cell_lookup <- function(design, bs_grid, bi_grid, ls_grid) {
  if (!length(design$factor_names)) return(array(1L, dim = c(1L, 1L, 1L)))
  n_bs <- nrow(bs_grid)
  n_bi <- nrow(bi_grid)
  n_ls <- nrow(ls_grid)
  arr <- array(0L, dim = c(n_bs, n_bi, n_ls))
  cells_chr <- as.data.frame(
    lapply(design$cells[, design$factor_names, drop = FALSE], as.character),
    stringsAsFactors = FALSE
  )
  key_cells <- apply(cells_chr, 1L, paste, collapse = "\r")

  for (a in seq_len(n_bs)) {
    for (b in seq_len(n_bi)) {
      for (c_idx in seq_len(n_ls)) {
        vals <- character(length(design$factor_names))
        names(vals) <- design$factor_names
        for (fn in design$between_subj) vals[fn] <- bs_grid[a, fn]
        for (fn in design$between_item) vals[fn] <- bi_grid[b, fn]
        for (fn in design$ls_factors)   vals[fn] <- ls_grid[c_idx, fn]
        k <- paste(vals, collapse = "\r")
        arr[a, b, c_idx] <- match(k, key_cells)
      }
    }
  }
  arr
}

#' Print or Summarize a Planning Object
#'
#' The print methods display the design, calibration, target, or selected plan.
#' A plan summary also displays available response-pattern counts, verified
#' designs, and sensitivity results.
#' @name print.glmmplanr_design
#' @param x,object A design, calibrated generator, generator check, target
#'   specification, or planning result of the corresponding class.
#' @param ... Additional arguments passed to the method.
#' @return The input object with its original class, invisibly.
#' @examples
#' des <- design_crossed(list(condition = c("a", "b")))
#' print(des)
#' @export
print.glmmplanr_design <- function(x, ...) {
  cli::cli_h2("glmmplanr Experimental Design")
  fac_desc <- vapply(
    x$factor_names,
    function(fn) sprintf("%s (%s) [%s-participant, %s-item]", fn,
                         paste(x$factors[[fn]], collapse = ", "),
                         if (fn %in% x$between_subj) "between" else "within",
                         if (fn %in% x$between_item) "between" else "within"),
    character(1L)
  )
  cli::cli_ul(c(
    if (length(fac_desc)) "Factors ({length(x$factor_names)}): {paste(fac_desc, collapse = ' x ')}" else "Intercept-only design",
    "Allocation: {.val {x$allocation}} ({x$n_cells} condition cells)",
    "List-balance multiples: N_subj = {x$k_subj}k, N_item = {x$k_item}m",
    if (!is.null(x$items_per_subject)) "Items drawn per participant: {x$items_per_subject} (allocation seed {x$item_seed})",
    if (x$has_testlet) "Testlet nesting: {x$testlet_size} items per passage/testlet",
    "Factor coding: {.val {x$contrast}} | Trial missingness: {round(x$attrition * 100, 1)}% | Participant dropout: {round(dropout_rate(x) * 100, 1)}%",
    "Fixed formula: {.code {deparse1(x$fixed_rhs)}}",
    "Subject RE: {.code {deparse1(x$subj_rhs)}} | Item RE: {.code {deparse1(x$item_rhs)}}"
  ))
  invisible(x)
}
