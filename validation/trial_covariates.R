# Run at the repository root with Rscript --vanilla validation/trial_covariates.R
library(lme4)

n_person <- 64L
n_item <- 16L
n_rep <- 64L
sd_person_x <- .5
sd_trial_x <- 1
beta_within <- .6
beta_between <- -.4
split_model <- y ~ within + between + (1 | person) + (1 | item)
pooled_model <- y ~ pooled + (1 | person) + (1 | item)
source("validation/covariate_fit.R")

simulate_one <- function(replication) {
  seed <- 613800L + replication
  set.seed(seed)
  d <- expand.grid(person = factor(seq_len(n_person)), item = factor(seq_len(n_item)))
  pid <- as.integer(d$person)
  iid <- as.integer(d$item)
  mu <- rnorm(n_person, sd = sd_person_x)
  trial <- rnorm(nrow(d), sd = sd_trial_x)
  d$x <- mu[pid] + trial
  eta <- .35 + beta_within * trial + beta_between * mu[pid] +
    rnorm(n_person, sd = .55)[pid] + rnorm(n_item, sd = .35)[iid]
  u <- runif(nrow(d))
  # Nested subsets: all 16 items remain represented and equally frequent.
  phase <- sample(rep(0:(n_item - 1L), each = n_person/n_item))
  item_order <- sample.int(n_item) - 1L
  rank <- (phase[pid] + item_order[iid]) %% n_item

  rows <- list()
  for (m in c(4L, 16L)) {
    keep <- rank < m
    z <- d[keep, ]
    stopifnot(all(table(z$person) == m), all(table(z$item) == n_person*m/n_item))
    person_means <- as.numeric(tapply(z$x, z$person, mean))
    sample_mean <- person_means[as.integer(z$person)]
    grand <- mean(z$x)
    z$pooled <- z$x - grand
    rel <- sd_person_x^2 / (sd_person_x^2 + sd_trial_x^2/m)
    # Balanced one-way ANOVA on the covariate. A zero between variance makes
    # the plug-in correction unidentified and is retained as a failed estimate.
    ms_within <- sum((z$x - sample_mean)^2) / (n_person * (m - 1))
    ms_between <- m * var(person_means)
    rel_hat <- max(0, 1 - ms_within/ms_between)
    observed_between <- beta_within + rel * (beta_between - beta_within)

    for (link in c("logit", "probit")) {
      z$y <- as.integer(u[keep] < binomial(link)$linkinv(eta[keep]))
      make_row <- function(method, effect, fit, weights, target, interval = TRUE) {
        valid <- fit$ok && all(is.finite(weights))
        estimate <- se <- NA_real_
        if (valid) {
          estimate <- sum(weights * fit$b)
          if (interval) se <- sqrt(drop(crossprod(weights, fit$V %*% weights)))
        }
        covered <- if (valid && interval && is.finite(target)) abs(estimate-target) <= qnorm(.975)*se else NA
        data.frame(replication, seed, link, items_per_person = m, method, effect,
                   target, estimate, se, covered, valid, singular = fit$singular,
                   reliability = rel, estimated_reliability = rel_hat,
                   diagnostic = if (!all(is.finite(weights))) "Zero estimated reliability" else fit$diagnostic)
      }

      z$within <- z$x - mu[as.integer(z$person)]
      z$between <- mu[as.integer(z$person)] - grand
      known <- fit_slopes(z, split_model, link, c("within", "between"))
      z$within <- z$x - sample_mean
      z$between <- sample_mean - grand
      stopifnot(max(abs(tapply(z$within, z$person, mean))) < 1e-10,
                max(abs(z$within + z$between - z$pooled)) < 1e-10)
      observed <- fit_slopes(z, split_model, link, c("within", "between"))
      pooled <- fit_slopes(z, pooled_model, link, "pooled")
      # Shrinking the person mean by rel just reparameterizes the split fit.
      # Its latent between slope uses the covariance of BOTH fitted slopes.
      correction <- c(1 - 1/rel, 1/rel)
      estimated_correction <- if (rel_hat > 0) c(1 - 1/rel_hat, 1/rel_hat) else c(NA_real_, NA_real_)
      rows[[length(rows) + 1L]] <- rbind(
        make_row("latent_mean_known", "within", known, c(1, 0), beta_within),
        make_row("latent_mean_known", "between", known, c(0, 1), beta_between),
        make_row("observed_mean", "within", observed, c(1, 0), beta_within),
        make_row("observed_mean", "between", observed, c(0, 1), observed_between),
        make_row("known_reliability", "between", observed, correction, beta_between),
        make_row("estimated_reliability", "between", observed, estimated_correction, beta_between, interval = FALSE),
        make_row("grand_mean_only", "pooled", pooled, 1, NA_real_))
      if (replication == 1L) {
        # Check the analytic reparameterization against an independent fit.
        eb <- rel * sample_mean
        z$within <- z$x - eb
        z$between <- eb - grand
        check <- fit_slopes(z, split_model, link, c("within", "between"))
        stopifnot(observed$ok, check$ok,
          abs(check$b[2] - sum(correction * observed$b)) < 1e-3,
          abs(check$V[2,2] - drop(crossprod(correction, observed$V %*% correction))) < 1e-3)
      }
    }
  }
  do.call(rbind, rows)
}

runs <- parallel::mclapply(seq_len(n_rep), simulate_one,
  mc.cores = if (.Platform$OS.type == "windows") 1L else 2L, mc.set.seed = FALSE)
stopifnot(all(vapply(runs, is.data.frame, logical(1))))
raw <- do.call(rbind, runs)
keys <- c("link", "items_per_person", "method", "effect")
groups <- split(raw, do.call(interaction, c(raw[keys], list(drop = TRUE))))
summary <- do.call(rbind, lapply(groups, function(d) {
  z <- d[d$valid, ]
  measured <- sum(!is.na(z$covered))
  n_covered <- if (measured) sum(z$covered, na.rm = TRUE) else NA_integer_
  ci <- if (measured) binom.test(n_covered, measured)$conf.int else c(NA_real_, NA_real_)
  error <- z$estimate - d$target[1]
  data.frame(d[1, keys], target = d$target[1], attempted = nrow(d), valid = nrow(z),
    singular_valid = sum(z$singular), mean_estimate = mean(z$estimate), mean_mcse = sd(z$estimate)/sqrt(nrow(z)),
    bias = mean(error), rmse = sqrt(mean(error^2)), mean_se = if (measured) mean(z$se) else NA_real_,
    n_covered, n_intervals = measured, coverage = if (measured) n_covered/measured else NA_real_,
    coverage_lo = ci[1], coverage_hi = ci[2], R_version = as.character(getRversion()),
    lme4_version = as.character(packageVersion("lme4")))
}))
stopifnot(nrow(raw) == 28L * n_rep, nrow(summary) == 28L,
          all(summary$attempted == n_rep), all(is.finite(raw$estimate[raw$valid])),
          !anyDuplicated(raw[c("replication", keys)]))
write.csv(raw, "validation/trial_covariate_fits.csv", row.names = FALSE)
write.csv(summary, "validation/trial_covariate_summary.csv", row.names = FALSE)
print(summary[c(keys, "target", "valid", "mean_estimate", "rmse", "n_covered")], row.names = FALSE)
