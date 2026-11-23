# Internal cache for 1D Gauss-Hermite quadrature rules and compiled Stan models
.glmmplanr_env <- new.env(parent = emptyenv())
.glmmplanr_env$gh1d <- list()

#' Standard-Normal Gauss-Hermite Quadrature Rule (1D)
#'
#' Computes nodes and weights for integrating against the standard normal
#' density \eqn{\phi(z; 0, 1)} via the Golub-Welsch algorithm on the
#' probabilist Hermite Jacobi matrix.
#'
#' @param n_nodes Integer number of quadrature points (default 21).
#' @return A list with numeric vectors `nodes` and `weights` (summing to 1).
#' @noRd
gh_quad_1d <- function(n_nodes = 21L) {
  n_nodes <- as.integer(n_nodes)
  key <- as.character(n_nodes)
  if (!is.null(.glmmplanr_env$gh1d[[key]])) {
    return(.glmmplanr_env$gh1d[[key]])
  }
  if (n_nodes < 1L) {
    cli::cli_abort("{.arg n_nodes} must be a positive integer.")
  }
  if (n_nodes == 1L) {
    res <- list(nodes = 0, weights = 1)
    .glmmplanr_env$gh1d[[key]] <- res
    return(res)
  }
  off_diag <- sqrt(seq_len(n_nodes - 1L))
  jacobi <- matrix(0, nrow = n_nodes, ncol = n_nodes)
  idx <- cbind(seq_len(n_nodes - 1L), seq.int(2L, n_nodes))
  jacobi[idx] <- off_diag
  jacobi[idx[, c(2L, 1L), drop = FALSE]] <- off_diag
  eig <- eigen(jacobi, symmetric = TRUE)
  ord <- order(eig$values)
  nodes <- eig$values[ord]
  weights <- (eig$vectors[1L, ord])^2
  weights <- weights / sum(weights)
  res <- list(nodes = nodes, weights = weights)
  .glmmplanr_env$gh1d[[key]] <- res
  res
}

#' Modified Logit Link and Inverse Link for m-AFC Tasks
#'
#' Maps between linear predictor \eqn{\eta \in (-\infty, \infty)} and response
#' probability \eqn{P \in (g, 1)} with lower guessing asymptote \eqn{g \in [0, 1)}:
#' \eqn{P = g + (1 - g)\,\mathrm{logit}^{-1}(\eta)}.
#'
#' @param p Numeric vector of probabilities strictly in `(chance, 1)`.
#' @param eta Numeric vector of conditional linear predictors.
#' @param chance Lower asymptote (guessing probability \eqn{g \in [0, 1)}, e.g. `1/3` for 3AFC).
#' @return Numeric vector on the conditional logit or probability scale.
#' @examples
#' eta <- mafc_qlogis(0.7, chance = 1/3)
#' mafc_plogis(eta, chance = 1/3)
#' @export
mafc_qlogis <- function(p, chance = 0) {
  validate_chance(chance)
  if (any(is.na(p)) || any(p <= chance) || any(p >= 1)) {
    cli::cli_abort(
      "All probabilities {.arg p} must strictly lie in ({round(chance, 4)}, 1)."
    )
  }
  log((p - chance) / (1 - p))
}

#' @rdname mafc_qlogis
#' @export
mafc_plogis <- function(eta, chance = 0) {
  validate_chance(chance)
  chance + (1 - chance) * stats::plogis(eta)
}

#' Convert Intuitive Accuracy Range to Random-Effect Standard Deviation
#'
#' Converts a substantive empirical range of subject or item response
#' probabilities (e.g., 95% of subjects fall between 40% and 85% accuracy at a
#' baseline condition) into the corresponding random-intercept standard
#' deviation \eqn{\tau_0} on the conditional (guessing-adjusted) logit scale.
#'
#' @param prob_range Numeric vector `c(p_low, p_high)` with `chance < p_low < p_high < 1`.
#' @param chance Guessing probability \eqn{g \in [0, 1)} (default `0`).
#' @param ci_level Central interval probability mass (default `0.95`).
#' @return Positive numeric scalar giving the random-effect SD on the conditional logit scale.
#' @examples
#' prob_range_to_sd(c(0.40, 0.85), chance = 1 / 3)
#' @export
prob_range_to_sd <- function(prob_range, chance = 0, ci_level = 0.95) {
  validate_chance(chance)
  if (!is.numeric(prob_range) || length(prob_range) != 2L || anyNA(prob_range)) {
    cli::cli_abort("{.arg prob_range} must be a numeric vector of length 2: c(low, high).")
  }
  if (prob_range[1L] >= prob_range[2L]) {
    cli::cli_abort("{.arg prob_range[1]} must be strictly less than {.arg prob_range[2]}.")
  }
  if (prob_range[1L] <= chance || prob_range[2L] >= 1) {
    cli::cli_abort(
      "{.arg prob_range} ({round(prob_range[1], 3)}, {round(prob_range[2], 3)}) must lie strictly inside ({round(chance, 3)}, 1)."
    )
  }
  if (!is.numeric(ci_level) || length(ci_level) != 1L || ci_level <= 0 || ci_level >= 1) {
    cli::cli_abort("{.arg ci_level} must be a probability strictly between 0 and 1.")
  }
  z_crit <- stats::qnorm(1 - (1 - ci_level) / 2)
  eta_bounds <- mafc_qlogis(prob_range, chance = chance)
  unname((eta_bounds[2L] - eta_bounds[1L]) / (2 * z_crit))
}

#' Convert Random-Effect SD Back to Conditional Probability Interval
#'
#' @param sd Positive numeric standard deviation on the conditional logit scale.
#' @param intercept Conditional logit mean at reference condition (default `0`).
#' @param chance Guessing probability \eqn{g \in [0, 1)} (default `0`).
#' @param ci_level Central interval probability mass (default `0.95`).
#' @return Numeric vector `c(low, high)` of response probabilities.
#' @examples
#' sd_to_prob_range(sd = 0.5, intercept = qlogis(0.7))
#' @export
sd_to_prob_range <- function(sd, intercept = 0, chance = 0, ci_level = 0.95) {
  validate_chance(chance)
  if (!is.numeric(sd) || length(sd) != 1L || is.na(sd) || sd < 0) {
    cli::cli_abort("{.arg sd} must be a non-negative scalar.")
  }
  z_crit <- stats::qnorm(1 - (1 - ci_level) / 2)
  mafc_plogis(intercept + c(-z_crit, z_crit) * sd, chance = chance)
}

#' Custom GLM Binomial Family with Lower Guessing Asymptote (m-AFC)
#'
#' Constructs a `family` object compatible with [stats::glm()] and
#' [lme4::glmer()] for binary responses with known lower asymptote `chance`.
#'
#' @param chance Guessing probability \eqn{g \in [0, 1)}. When `chance == 0`,
#'   returns standard `stats::binomial("logit")` so `lme4` uses its compiled C++
#'   fast path.
#' @return An object of class `"family"`.
#' @examples
#' family <- mafc_binomial(chance = 1/3)
#' family$linkinv(0)
#' @export
mafc_binomial <- function(chance = 0) {
  validate_chance(chance)
  if (chance == 0) {
    return(stats::binomial(link = "logit"))
  }
  g <- as.numeric(chance)
  eps <- 1e-10
  link <- structure(
    list(
      linkfun = function(mu) {
        mu_c <- pmin(pmax(mu, g + 1e-6), 1 - 1e-6)
        log((mu_c - g) / (1 - mu_c))
      },
      linkinv = function(eta) {
        p <- stats::plogis(eta)
        pmin(pmax(g + (1 - g) * p, g + eps), 1 - eps)
      },
      mu.eta = function(eta) {
        p <- stats::plogis(eta)
        pmax((1 - g) * p * (1 - p), eps)
      },
      valideta = function(eta) TRUE,
      name = sprintf("mafc_logit(%.4f)", g)
    ),
    class = "link-glm"
  )
  fam <- stats::binomial(link = link)
  fam$initialize <- bquote({
    if (NCOL(y) == 1L) {
      if (is.factor(y)) y <- y != levels(y)[1L]
      n <- rep.int(1, nobs)
      y[weights == 0] <- 0
      mustart <- .(g) + (1 - .(g)) * (weights * y + 0.5) / (weights + 1)
      m <- weights * y
    } else {
      stop("mafc_binomial requires binary 0/1 response vector.", call. = FALSE)
    }
  })
  fam
}

validate_chance <- function(chance) {
  if (!is.numeric(chance) || length(chance) != 1L || is.na(chance) ||
      chance < 0 || chance >= 1) {
    cli::cli_abort("{.arg chance} must be a numeric scalar in [0, 1).")
  }
  invisible(chance)
}

# Box-Muller keeps the second normal draw outside .Random.seed. R exposes
# no way to save that value before set.seed(), even when the seed is restored.
uses_box_muller <- function() {
  exists(".Random.seed", .GlobalEnv, inherits = FALSE) && RNGkind()[2L] == "Box-Muller"
}

isolated_seeded_call <- function(fun, frame, call) {
  supplied <- names(as.list(call)[-1L])
  args <- mget(supplied, frame, inherits = FALSE)
  path <- getNamespaceInfo(asNamespace("glmmplanr"), "path")
  callr::r(function(fun, args, path) {
    if (!file.exists(file.path(path, "Meta", "package.rds"))) {
      pkgload::load_all(path, quiet = TRUE)
    }
    do.call(get(fun, asNamespace("glmmplanr"), inherits = FALSE), args)
  }, list(fun = fun, args = args, path = path), show = TRUE)
}

#' Generate Deterministic L'Ecuyer-CMRG Stream Seeds
#'
#' Creates independent reproducible L'Ecuyer-CMRG `.Random.seed` vectors while
#' preserving an existing caller stream when a master seed is supplied.
#'
#' @param seed Integer scalar master seed (or `NULL` to draw one seed).
#' @param n_streams Number of independent streams required.
#' @return A list of length `n_streams`, each an integer L'Ecuyer-CMRG `.Random.seed`.
#' @noRd
make_lecuyer_streams <- function(seed, n_streams) {
  n_streams <- as.integer(n_streams)
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1L)
  s <- with_rng_stream(seed, get(".Random.seed", .GlobalEnv), kind = "L'Ecuyer-CMRG")
  streams <- vector("list", n_streams)
  for (i in seq_len(n_streams)) {
    s <- parallel::nextRNGStream(s)
    streams[[i]] <- s
  }
  streams
}

#' Run an Expression with an Isolated Stream or Allocation Seed
#' @param stream_seed Integer seed or L'Ecuyer-CMRG state vector.
#' @param expr Expression to evaluate.
#' @param kind Uniform generator used for an integer seed.
#' @return The result of evaluating `expr`.
#' @noRd
with_rng_stream <- function(stream_seed, expr, kind = "Mersenne-Twister") {
  has_old <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_kind <- RNGkind()
  if (has_old) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (has_old) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else rm(".Random.seed", envir = .GlobalEnv)
  })
  if (length(stream_seed) == 1L) {
    set.seed(stream_seed, kind, "Inversion", "Rejection")
  } else {
    RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
    assign(".Random.seed", stream_seed, envir = .GlobalEnv)
  }
  force(expr)
}

#' Allocate Strictly Disjoint L'Ecuyer-CMRG Base Streams by (Stage, Point)
#'
#' Maps each `(stage_slot, grid_id)` pair to a unique top-level stream index
#' `K = (stage_slot - 1L) * stride + grid_id` via `parallel::nextRNGStream`
#' (step `2^127`), while individual replications within a `(stage_slot, grid_id)`
#' advance via `parallel::nextRNGSubStream` (step `2^76` in
#' `derive_rep_substreams()`). Because `1 <= rep <= 2^51`, two states coincide
#' if and only if both `stage_slot`, `grid_id`, and `rep` are identical.
#'
#' @param seed Integer scalar master seed.
#' @param stage_slot Positive integer stage/scenario identifier (`1L` = Stage 1
#'   exploration, `2L` = Stage 2 verification, `3L` = Stage 2 top-up, `4L+` =
#'   sensitivity / robust scenarios).
#' @param grid_ids Integer vector of `1`-based lattice point IDs (`<= stride`).
#' @param stride Positive integer >= maximum `grid_id` across the lattice.
#' @return A list of L'Ecuyer-CMRG base `.Random.seed` vectors aligned with
#'   `grid_ids`.
#' @noRd
allocate_stage_point_streams <- function(seed, stage_slot, grid_ids, stride) {
  stage_slot <- as.integer(stage_slot)
  grid_ids <- as.integer(grid_ids)
  stride <- as.integer(stride)
  if (stage_slot < 1L || stride < 1L || any(grid_ids < 1L) || any(grid_ids > stride)) {
    cli::cli_abort("Invalid (stage_slot, grid_ids, stride) in RNG stream allocation.")
  }

  target_k <- (stage_slot - 1L) * stride + grid_ids
  max_k <- max(target_k)

  s <- with_rng_stream(seed, get(".Random.seed", .GlobalEnv), kind = "L'Ecuyer-CMRG")

  lookup <- vector("list", max_k)
  want <- logical(max_k)
  want[target_k] <- TRUE
  for (k in seq_len(max_k)) {
    s <- parallel::nextRNGStream(s)
    if (want[k]) lookup[[k]] <- s
  }
  lookup[target_k]
}

#' Derive Independent Per-Replication Substreams Within a (Stage, Point) Stream
#'
#' Advances `parallel::nextRNGSubStream` (`2^76` period spacing) so replications
#' never collide with adjacent `(stage_slot, grid_id)` top-level streams (`2^127`
#' spacing).
#'
#' @param point_base L'Ecuyer-CMRG base seed for `(stage_slot, grid_id)`.
#' @param nsim Number of replication substreams to generate.
#' @param offset Non-negative integer number of initial substreams to skip
#'   (used when appending independent top-up replications).
#' @return List of length `nsim` of `.Random.seed` vectors.
#' @noRd
derive_rep_substreams <- function(point_base, nsim, offset = 0L) {
  nsim <- as.integer(nsim)
  offset <- as.integer(offset)
  s <- point_base
  if (offset > 0L) {
    for (j in seq_len(offset)) {
      s <- parallel::nextRNGSubStream(s)
    }
  }
  out <- vector("list", nsim)
  for (b in seq_len(nsim)) {
    s <- parallel::nextRNGSubStream(s)
    out[[b]] <- s
  }
  out
}

#' Wilson Score Monte Carlo Confidence Interval and Three-Way Decision Rule
#'
#' Computes Monte Carlo uncertainty intervals for simulated proportions
#' following Morris, White, and Crowther (2019, `doi:10.1002/sim.8086`) and
#' classifies attainment against `threshold` into `"pass"`, `"fail"`, or
#' `"indeterminate"` (when the MC interval straddles `threshold`).
#'
#' @noRd
wilson_mc_ci <- function(k, n, conf_level = 0.95) {
  if (is.na(k) || is.na(n) || n <= 0L) {
    return(c(lo = NA_real_, hi = NA_real_))
  }
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  z2 <- z * z
  p <- k / n
  denom <- 1 + z2 / n
  center <- (p + z2 / (2 * n)) / denom
  half <- (z * sqrt((p * (1 - p) + z2 / (4 * n)) / n)) / denom
  c(lo = max(0.0, center - half), hi = min(1.0, center + half))
}

bonferroni_mc_conf_level <- function(conf_level, n_tests = 1L) {
  n_tests <- max(1L, as.integer(n_tests))
  1.0 - (1.0 - conf_level) / n_tests
}

classify_mc_threshold <- function(ci_lo, ci_hi, threshold) {
  if (is.na(ci_lo) || is.na(ci_hi) || is.na(threshold)) {
    return("fail")
  }
  if (ci_lo >= threshold) {
    "pass"
  } else if (ci_hi < threshold) {
    "fail"
  } else {
    "indeterminate"
  }
}
