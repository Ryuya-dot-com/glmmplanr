# Run at the repository root with Rscript --vanilla validation/covariate_allocation.R
library(lme4)
source("validation/covariate_fit.R")

target <- -.4
split_model <- y ~ within + between + (1 | person) + (1 | item)

# Gaussian random-intercept likelihood for X only. Profile out the trial
# variance. The remaining parameters are the population mean and log ratio
# of trial variance to between-person variance.
covariate_ml <- function(count, xbar, within_ss) {
  n <- sum(count)
  nll <- function(par) {
    k <- exp(par[2])
    q <- within_ss + sum(count * (xbar - par[1])^2 / (1 + count/k))
    .5 * (n * log(q/n) + sum(log1p(count/k)))
  }
  mu <- function(logk) weighted.mean(xbar, count / (1 + count/exp(logk)))
  profile <- function(logk) nll(c(mu(logk), logk))
  bounds <- c(-15, 15)
  logk <- optimize(profile, bounds, tol = 1e-8)$minimum
  par <- c(mu = mu(logk), logk = logk)
  result <- list(par = par, V = NULL, diagnostic = "")
  if (min(abs(logk - bounds)) < .001) {
    result$par[] <- NA_real_
    result$diagnostic <- "Variance ratio at search boundary"
    return(result)
  }
  info <- optimHess(par, nll)
  if (all(is.finite(info)) && min(eigen(info, symmetric = TRUE, only.values = TRUE)$values) > 0) {
    result$V <- chol2inv(chol(info))
  } else result$diagnostic <- "Nonpositive covariate information"
  result
}

allocation_intervals <- function(z, count, xbar, nuisance, link) {
  idx <- as.integer(z$person)
  oracle <- c(mu = 0, logk = log(4))
  refit <- function(par) {
    rel <- count/(count + exp(par[2]))
    eb <- par[1] + rel * (xbar - par[1])
    z$between <- eb[idx]
    z$within <- z$x - z$between
    fit_slopes(z, split_model, link, c("within", "between"))
  }
  balanced <- all(count == count[1])
  if (balanced) observed <- refit(c(0, -Inf))
  result <- list()
  for (method in c("known", "estimated")) {
    par <- if (method == "known") oracle else nuisance$par
    fit <- list(ok = FALSE, singular = NA, diagnostic = nuisance$diagnostic)
    estimate <- se <- propagated <- NA_real_
    extra_diagnostic <- ""
    if (all(is.finite(par))) {
      fit <- if (balanced) observed else refit(par)
      if (fit$ok) {
        weights <- if (balanced) c(0, 1) + c(-1, 1)*exp(par[2])/count[1] else c(0, 1)
        estimate <- sum(weights * fit$b)
        se <- sqrt(drop(crossprod(weights, fit$V %*% weights)))
        if (method == "estimated" && !is.null(nuisance$V)) {
          if (balanced) {
            gradient <- c(0, exp(par[2])/count[1] * (fit$b[2] - fit$b[1]))
          } else {
            # Change each first-stage parameter and re-fit the outcome model.
            gradient <- numeric(2)
            for (j in 1:2) {
              step <- numeric(2)
              step[j] <- .02
              plus <- refit(par + step)
              minus <- refit(par - step)
              if (plus$ok && minus$ok) {
                gradient[j] <- (plus$b[2] - minus$b[2])/(2 * step[j])
              } else {
                gradient[j] <- NA_real_
                extra_diagnostic <- paste("Derivative fit failed", j, plus$diagnostic, minus$diagnostic)
              }
            }
          }
          propagated <- sqrt(se^2 + drop(crossprod(gradient, nuisance$V %*% gradient)))
        } else if (method == "estimated") extra_diagnostic <- nuisance$diagnostic
      }
    }
    ses <- if (method == "known") c(known = se) else c(plugin = se, delta = propagated)
    result[[method]] <- data.frame(method = names(ses), estimate, se = unname(ses),
      singular = fit$singular, diagnostic = ifelse(names(ses) == "delta" & nzchar(extra_diagnostic),
        extra_diagnostic, fit$diagnostic))
  }
  do.call(rbind, result)
}

simulate_allocation <- function(replication) {
  seed <- 772500L + replication
  set.seed(seed)
  d <- expand.grid(person = seq_len(128L), item = seq_len(16L))
  pid <- d$person
  latent <- rnorm(128L, sd = .5)
  noise <- rnorm(nrow(d))
  d$x <- latent[pid] + noise
  eta <- .35 + .6*noise + target*latent[pid] +
    rnorm(128L, sd = .55)[pid] + rnorm(16L, sd = .35)[d$item]
  u <- runif(nrow(d))
  d$logit <- as.integer(u < plogis(eta))
  d$probit <- as.integer(u < pnorm(eta))
  retain <- runif(nrow(d))
  phase <- unlist(replicate(8L, sample(0:15), simplify = FALSE))
  rank <- (phase[pid] + sample(0:15)[d$item]) %% 16L
  rows <- list()
  for (m in c(16L, 8L, 4L)) {
    n_person <- 512L %/% m
    scheduled <- pid <= n_person & rank < m
    stopifnot(sum(scheduled) == 512L, all(table(pid[scheduled]) == m),
              all(table(d$item[scheduled]) == 32L))
    for (missing in c(0, .25)) {
      z <- d[scheduled & retain >= missing, ]
      z$person <- factor(z$person)
      z$item <- factor(z$item)
      idx <- as.integer(z$person)
      count <- tabulate(idx)
      xbar <- as.numeric(tapply(z$x, z$person, mean))
      nuisance <- covariate_ml(count, xbar, sum((z$x - xbar[idx])^2))
      k_hat <- exp(nuisance$par[2])
      if (replication == 1L) {
        # Independent check of the profiled ML covariate stage.
        reference <- lmer(x ~ 1 + (1 | person), z, REML = FALSE)
        stopifnot(abs(k_hat * getME(reference, "theta")^2 - 1) < .001,
                  abs(nuisance$par[1] - fixef(reference)[1]) < .001)
      }
      for (link in c("logit", "probit")) {
        z$y <- z[[link]]
        res <- allocation_intervals(z, count, xbar, nuisance, link)
        rows[[length(rows) + 1L]] <- data.frame(replication, seed, n_person,
          items_per_person = m, missing, link, scheduled = 512L, observed = nrow(z),
          active_persons = length(count), min_observed = if (length(count) < n_person) 0L else min(count),
          max_observed = max(count), mean_reliability = mean(count/(count + 4)),
          estimated_ratio = unname(k_hat), estimated_mean = unname(nuisance$par[1]), res)
      }
    }
  }
  raw <- do.call(rbind, rows)
  raw$lower <- raw$estimate - qnorm(.975)*raw$se
  raw$upper <- raw$estimate + qnorm(.975)*raw$se
  raw$covered <- ifelse(is.finite(raw$lower) & is.finite(raw$upper),
                        raw$lower <= target & target <= raw$upper, NA)
  raw
}

summarize_allocation <- function(raw) {
  keys <- c("link", "n_person", "items_per_person", "missing", "method")
  groups <- split(raw, do.call(interaction, c(raw[keys], list(drop = TRUE))))
  ans <- do.call(rbind, lapply(groups, function(d) {
    point <- is.finite(d$estimate)
    valid <- !is.na(d$covered)
    n <- sum(valid)
    hits <- sum(d$covered[valid])
    ci <- if (n) binom.test(hits, n)$conf.int else c(NA_real_, NA_real_)
    error <- d$estimate[point] - target
    data.frame(d[1, keys], target, attempted = nrow(d), points = sum(point),
      intervals = n, singular_intervals = sum(d$singular[valid]), covered = hits,
      coverage = hits/n, coverage_lo = ci[1], coverage_hi = ci[2],
      bias = mean(error), rmse = sqrt(mean(error^2)), mean_width = mean(2*qnorm(.975)*d$se[valid]),
      mean_observed = mean(d$observed), mean_active_persons = mean(d$active_persons),
      mean_reliability = mean(d$mean_reliability), R_version = as.character(getRversion()),
      lme4_version = as.character(packageVersion("lme4")))
  }))
  ans[order(ans$link, ans$missing, ans$n_person, ans$method), ]
}

if (sys.nframe() == 0L) {
  n_rep <- 64L
  runs <- parallel::mclapply(seq_len(n_rep), simulate_allocation,
    mc.cores = if (.Platform$OS.type == "windows") 1L else 2L, mc.set.seed = FALSE)
  stopifnot(all(vapply(runs, is.data.frame, logical(1))))
  raw <- do.call(rbind, runs)
  summary <- summarize_allocation(raw)
  delta <- raw[raw$method == "delta", ]
  plugin <- raw[raw$method == "plugin", ]
  stopifnot(nrow(summary) == 36L, all(summary$attempted == n_rep),
    !anyDuplicated(raw[c("replication", "n_person", "missing", "link", "method")]),
    identical(delta$estimate, plugin$estimate), all(delta$se >= plugin$se, na.rm = TRUE),
    all(raw$observed[raw$missing == 0] == 512L))
  write.csv(raw, "validation/covariate_allocation_fits.csv", row.names = FALSE)
  write.csv(summary, "validation/covariate_allocation_summary.csv", row.names = FALSE)
  print(summary[summary$method == "delta", c("link", "n_person", "missing", "intervals",
        "covered", "rmse", "mean_width")], row.names = FALSE)
}
