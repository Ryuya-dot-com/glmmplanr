# Exploratory case check of Laplace vs. direct numerical integration on extreme RR fits.
library(glmmplanr)
args <- commandArgs(TRUE)
out <- if (length(args)) args[1L] else "validation"
prefix <- file.path(out, "two_pl_fixed_bank_medium")
cases <- read.csv(paste0(prefix, "_scenarios.csv"))
truth <- read.csv(paste0(prefix, "_truth.csv"))
fits <- read.csv(paste0(prefix, "_fits.csv.gz"))
items <- read.csv(paste0(prefix, "_items.csv.gz"))
stopifnot(length(unique(cases$seed)) == 1L, length(unique(cases$n_rep)) == 1L)
streams <- glmmplanr:::make_lecuyer_streams(cases$seed[1L], cases$n_rep[1L])
rr_streams <- glmmplanr:::make_lecuyer_streams(cases$seed[1L] + 100L, cases$n_rep[1L])
max_n <- max(cases$n)
max_j <- max(cases$items)
logits <- function(theta, a, d) sweep(outer(theta, a), 2L, d, "+")

integrated_loglik <- function(y, a, d) {
  patterns <- apply(y, 1L, paste0, collapse = "")
  keys <- unique(patterns)
  rows <- match(keys, patterns)
  count <- tabulate(match(patterns, keys), length(keys))
  b <- -d / a
  # Splitting the integral at the response-curve midpoints resolves sharp transitions.
  bounds <- sort(unique(c(-10, b[b > -10 & b < 10], 10)))
  prob <- vapply(rows, function(i) {
    observed <- y[i, ]
    density <- function(theta) {
      eta <- logits(theta, a, d)
      lp <- plogis(eta, log.p = TRUE) %*% observed +
        plogis(eta, lower.tail = FALSE, log.p = TRUE) %*% (1 - observed)
      exp(as.numeric(lp) + dnorm(theta, log = TRUE))
    }
    sum(vapply(seq_len(length(bounds) - 1L), function(s) {
      integrate(density, bounds[s], bounds[s + 1L], rel.tol = 1e-9, abs.tol = 1e-14)$value
    }, numeric(1L)))
  }, numeric(1L))
  sum(count * log(prob))
}

selected <- subset(cases, profile == "wide" & items == 8L)
results <- lapply(selected$case, function(k) {
  n <- cases$n[k]
  j <- cases$items[k]
  bank <- subset(truth, case == k)
  matched <- subset(fits, case == k & valid & model %in% c("glmmTMB_rr", "mirt_2pl"))
  both <- as.integer(names(which(table(matched$replication) == 2L)))
  r <- min(subset(items, case == k & model == "glmmTMB_rr" & replication %in% both & a > 5)$replication)
  stopifnot(is.finite(r))
  draw <- glmmplanr:::with_rng_stream(streams[[r]], list(theta = rnorm(max_n),
    u = matrix(runif(max_n * max_j), max_n, max_j)))
  y <- 1L * (draw$u[seq_len(n), seq_len(j)] < plogis(logits(draw$theta[seq_len(n)], bank$a_true, bank$d_true)))
  labels <- paste0("j", seq_len(j))
  dat <- data.frame(person = factor(rep(seq_len(n), j)),
                    item = factor(rep(labels, each = n), levels = labels), y = as.vector(y))
  fitted <- glmmplanr:::with_rng_stream(rr_streams[[r]], glmmTMB::glmmTMB(
    y ~ 0 + item + rr(0 + item | person, d = 1), data = dat, family = binomial(),
    control = glmmTMB::glmmTMBControl(start_method = list(method = "res"), parallel = 1L)))
  par <- fitted$fit$par
  b_idx <- names(par) == "beta"
  a_idx <- names(par) == "theta"
  direction <- sign(sum(par[a_idx]))
  rr <- subset(items, case == k & replication == r & model == "glmmTMB_rr")
  mirt <- subset(items, case == k & replication == r & model == "mirt_2pl")
  stopifnot(length(par) == 2L * j, max(abs(direction * par[a_idx] - rr$a)) < 1e-7,
            max(abs(par[b_idx] - rr$d)) < 1e-7)
  alt <- par
  alt[a_idx] <- direction * mirt$a
  alt[b_idx] <- mirt$d
  data.frame(case = k, n, replication = r, max_a_rr = max(rr$a), max_a_mirt = max(mirt$a),
    laplace_ll_rr = -fitted$obj$fn(par), laplace_ll_mirt = -fitted$obj$fn(alt),
    integrated_ll_rr = integrated_loglik(y, rr$a, rr$d), integrated_ll_mirt = integrated_loglik(y, mirt$a, mirt$d))
})
write.csv(do.call(rbind, results), paste0(prefix, "_likelihood_check.csv"), row.names = FALSE)
