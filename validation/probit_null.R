# Independent check of the 86% baseline, null treatment contrast in links_covariates.R.
library(lme4)

covariates <- paste0("x", 1:8)
formulas <- lapply(c(1L, 8L), function(k) reformulate(
  c("treatment", covariates[seq_len(k)], "(1 | subj)", "(1 | item)"), response = "y"))
intercept <- qnorm(.86) * sqrt(1 + .45^2 + .35^2 + .55^2)
n_rep <- 64L  # nominal 5% rejection has MCSE 2.7 percentage points
check_fit <- function(d, formula) {
  out <- data.frame(valid = FALSE, singular = NA, detected = FALSE, diagnostic = "")
  warnings <- character()
  tryCatch({
    fit <- withCallingHandlers(suppressMessages(glmer(formula, d, family = binomial("probit"),
      nAGQ = 1, control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 30000)))),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    msg <- fit@optinfo$conv$lme4$messages
    msg <- msg[!grepl("boundary (singular) fit", msg, fixed = TRUE)]
    out$singular <- isSingular(fit)
    out$diagnostic <- paste(unique(c(msg, warnings)), collapse = " | ")
    if (any(fit@optinfo$conv$opt != 0L) || length(msg)) return(out)
    b <- fixef(fit, add.dropped = TRUE)
    se <- sqrt(as.matrix(vcov(fit))["treatment", "treatment"])
    if (!all(is.finite(b)) || !is.finite(se) || se <= 0) stop("Invalid coefficient or treatment SE")
    out$valid <- TRUE
    out$detected <- abs(b["treatment"]) > qnorm(.975) * se
    out
  }, error = function(e) {
    out$diagnostic <- conditionMessage(e)
    out
  })
}
run_rep <- function(ns, r) {
  seed <- 953000L + 100L * ns + r
  set.seed(seed)
  ni <- 16L
  z <- matrix(rnorm(ns * 8L), nrow = ns)
  allocation <- sample(rep(0:1, each = ns / 2L))
  d <- expand.grid(subj = factor(seq_len(ns)), item = factor(seq_len(ni)))
  sid <- as.integer(d$subj)
  iid <- as.integer(d$item)
  d$treatment <- allocation[sid]
  lp <- intercept + .45*rnorm(ns)[sid] + .35*rnorm(ni)[iid] + .55*z[sid, 1]
  d$y <- rbinom(nrow(d), 1L, pnorm(lp))
  d$x1 <- z[sid, 1]
  one <- check_fit(d, formulas[[1L]])
  rows <- lapply(c(0, .85), function(rho) {
    sigma <- matrix(rho, 8, 8)
    diag(sigma) <- 1
    X <- z %*% chol(sigma)
    stopifnot(identical(X[sid, 1], d$x1))
    d[covariates] <- as.data.frame(X[sid, , drop = FALSE])
    data.frame(n_subj = ns, n_item = ni, replication = r, seed, rho, n_covariates = c(1L, 8L),
      n_perfect_person = sum(tabulate(sid[d$y == 1L], ns) == ni),
      n_zero_item = sum(tabulate(iid[d$y == 1L], ni) == 0L), rbind(one, check_fit(d, formulas[[2L]])))
  })
  do.call(rbind, rows)
}
raw <- do.call(rbind, lapply(c(32L, 64L), function(ns)
  do.call(rbind, lapply(seq_len(n_rep), function(r) run_rep(ns, r)))))
keys <- c("n_subj", "n_item", "rho", "n_covariates")
groups <- split(raw, do.call(interaction, c(raw[keys], list(drop = TRUE))))
summary <- do.call(rbind, lapply(groups, function(d) {
  ci <- binom.test(sum(d$detected), nrow(d))$conf.int
  data.frame(d[1, keys], attempted = nrow(d), valid = sum(d$valid), detected = sum(d$detected),
    rate = mean(d$detected), lo = ci[1], hi = ci[2], singular_valid = sum(d$valid & d$singular),
    datasets_with_perfect_person = sum(d$n_perfect_person > 0),
    datasets_with_zero_item = sum(d$n_zero_item > 0), R_version = as.character(getRversion()),
    lme4_version = as.character(packageVersion("lme4")))
}))
stopifnot(nrow(raw) == 8L * n_rep, nrow(summary) == 8L, all(summary$attempted == n_rep))
write.csv(raw, "validation/probit_null_fits.csv", row.names = FALSE)
write.csv(summary, "validation/probit_null_summary.csv", row.names = FALSE)
print(summary, row.names = FALSE)
