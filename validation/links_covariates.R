# Run from the repository root: Rscript --vanilla validation/links_covariates.R
# External benchmark for ordinary binary responses. Integrate over new
# participants, items and the known population distribution of the covariates.
library(lme4)

n_subj <- 32L
n_item <- 16L
n_rep <- 24L
covariates <- paste0("x", 1:8)
Sigma <- matrix(.85, 8, 8)
diag(Sigma) <- 1
dimnames(Sigma) <- list(covariates, covariates)
analysis_formulas <- function(target) {
  setNames(lapply(c(1L, 8L), function(k) {
    reformulate(c(target, covariates[seq_len(k)], "(1 | subj)", "(1 | item)"), response = "y")
  }), c("1", "8"))
}
main_formulas <- analysis_formulas("treatment")
interaction_formulas <- analysis_formulas("treatment * B * C")

marginal <- function(mu, sd, link) {
  if (link == "probit") return(pnorm(mu / sqrt(1 + sd^2)))
  integrate(function(z) plogis(mu + sd * z) * dnorm(z), -9, 9, abs.tol = 1e-10)$value
}

# Express heterogeneity relative to the latent residual SD for each link.
conditions <- expand.grid(link = c("logit", "probit"), baseline = c(.60, .86),
                          true_pp = c(0, 8), stringsAsFactors = FALSE)
conditions$latent_scale <- ifelse(conditions$link == "logit", pi / sqrt(3), 1)
eta <- function(p, link, scale) {
  sd <- scale * sqrt(.45^2 + .35^2 + .55^2)
  mu <- uniroot(function(b) marginal(b, sd, link) - p, c(-20, 20), tol = 1e-10)$root
  stopifnot(abs(marginal(mu, sd, link) - p) < 1e-8)
  mu
}
conditions$intercept <- mapply(eta, conditions$baseline, conditions$link, conditions$latent_scale)
conditions$effect <- mapply(eta, conditions$baseline + conditions$true_pp / 100,
                            conditions$link, conditions$latent_scale) - conditions$intercept

fit_once <- function(d, truth, fitted_link, k, cell_design = NULL, contrast = NULL) {
  formulas <- if (is.null(cell_design)) main_formulas else interaction_formulas
  warning_text <- character()
  row <- data.frame(valid = FALSE, singular = NA, detected = if (is.null(cell_design)) FALSE else NA,
                    estimated_pp = NA_real_, coef_covered = NA, diagnostic = "")
  tryCatch({
    fit <- withCallingHandlers(suppressMessages(glmer(formulas[[as.character(k)]], data = d,
      family = binomial(link = fitted_link), nAGQ = 1,
      control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 30000)))),
      warning = function(w) {
        warning_text <<- c(warning_text, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    row$singular <- isSingular(fit)
    messages <- fit@optinfo$conv$lme4$messages
    messages <- messages[!grepl("boundary \\(singular\\) fit", messages)]
    row$diagnostic <- paste(unique(c(messages, warning_text)), collapse = " | ")
    nonconverged <- any(fit@optinfo$conv$opt != 0L) || length(messages) > 0L
    if (nonconverged) return(row)
    b <- fixef(fit, add.dropped = TRUE)
    if (!all(is.finite(b))) stop("Invalid coefficient")
    xnames <- covariates[seq_len(k)]
    cov_var <- drop(crossprod(b[xnames], Sigma[xnames, xnames, drop = FALSE] %*% b[xnames]))
    re_var <- sum(vapply(VarCorr(fit), function(V) V[1, 1], numeric(1)))
    sd <- sqrt(re_var + cov_var)
    if (is.null(cell_design)) {
      cell_mu <- b[1] + c(0, b["treatment"])
      contrast <- c(-1, 1)
      se <- sqrt(as.matrix(vcov(fit))["treatment", "treatment"])
      if (!is.finite(se) || se <= 0) stop("Invalid treatment SE")
    } else {
      cell_mu <- drop(cell_design %*% b[colnames(cell_design)])
    }
    risk <- vapply(cell_mu, marginal, numeric(1), sd = sd, link = fitted_link)
    if (!all(is.finite(risk))) stop("Invalid standardized risk")
    row$estimated_pp <- 100 * sum(contrast * risk)
    row$valid <- TRUE
    if (is.null(cell_design)) {
      row$detected <- abs(b["treatment"]) > qnorm(.975) * se
      if (truth$link == fitted_link)
        row$coef_covered <- abs(b["treatment"] - truth$effect) <= qnorm(.975) * se
    }
    row
  }, error = function(e) {
    row$diagnostic <- conditionMessage(e)
    row
  })
}

records <- vector("list", n_rep)
for (rep_id in seq_len(n_rep)) {
  # Reuse the treatment assignment, predictors, latent draws and uniforms
  # across model conditions in a replication. Replications use separate seeds.
  set.seed(87100L + rep_id)
  X <- matrix(rnorm(n_subj * 8), n_subj, 8) %*% chol(Sigma)
  colnames(X) <- covariates
  d <- expand.grid(subj = factor(seq_len(n_subj)), item = factor(seq_len(n_item)))
  sid <- as.integer(d$subj)
  iid <- as.integer(d$item)
  treatment <- sample(rep(0:1, each = n_subj / 2))
  d$treatment <- treatment[sid]
  d <- cbind(d, X[sid, , drop = FALSE])
  nuisance <- .45 * rnorm(n_subj)[sid] + .35 * rnorm(n_item)[iid] + .55 * d$x1
  u <- runif(nrow(d))
  rows <- list()
  for (j in seq_len(nrow(conditions))) {
    truth <- conditions[j, ]
    lp <- truth$intercept + truth$effect * d$treatment + truth$latent_scale * nuisance
    d$y <- as.integer(u < if (truth$link == "logit") plogis(lp) else pnorm(lp))
    for (fitted_link in c("logit", "probit")) {
      for (k in c(1L, 8L)) {
        rows[[length(rows) + 1L]] <- data.frame(replication = rep_id, target = "main", dgp = truth$link,
          baseline = truth$baseline, true_pp = truth$true_pp, fitted_link = fitted_link,
          n_covariates = k, fit_once(d, truth, fitted_link, k))
      }
    }
  }
  records[[rep_id]] <- do.call(rbind, rows)
}
# A joint stress case: a full three-way interaction and up to eight covariates.
# Compare marginal point estimates; a latent interaction Wald test would test
# a different null from the response-scale difference-in-differences-in-differences.
cells <- expand.grid(treatment = 0:1, B = 0:1, C = 0:1)
cell_design <- model.matrix(~ treatment * B * C, cells)
contrast <- with(cells, (2*treatment-1) * (2*B-1) * (2*C-1))
interaction_conditions <- subset(conditions, baseline == .60)
cell_mu <- lapply(seq_len(nrow(interaction_conditions)), function(i) {
  cond <- interaction_conditions[i, ]
  p <- with(cells, .60 + .06*(treatment-.5) + .04*(B-.5) + .04*(C-.5) +
              .04*(treatment-.5)*(B-.5) + cond$true_pp/100*(treatment-.5)*(B-.5)*(C-.5))
  stopifnot(abs(sum(contrast*p)*100 - cond$true_pp) < 1e-8)
  vapply(p, eta, numeric(1), link = cond$link, scale = cond$latent_scale)
})
for (rep_id in seq_len(12L)) {
  set.seed(89100L + rep_id)
  X <- matrix(rnorm(n_subj * 8), n_subj, 8) %*% chol(Sigma)
  colnames(X) <- covariates
  d <- expand.grid(subj = factor(seq_len(n_subj)), item = factor(seq_len(n_item)))
  sid <- as.integer(d$subj)
  iid <- as.integer(d$item)
  subj_list <- sample(rep(0:3, each = n_subj / 4))[sid]
  item_list <- sample(rep(0:3, each = n_item / 4))[iid]
  d$treatment <- subj_list %/% 2L
  d$B <- item_list %/% 2L
  d$C <- (subj_list + item_list) %% 2L
  cell_id <- 1L + d$treatment + 2L*d$B + 4L*d$C
  stopifnot(all(tabulate(cell_id, 8L) == nrow(d) / 8))
  d <- cbind(d, X[sid, , drop = FALSE])
  nuisance <- .45*rnorm(n_subj)[sid] + .35*rnorm(n_item)[iid] + .55*d$x1
  u <- runif(nrow(d))
  rows <- list()
  for (j in seq_len(nrow(interaction_conditions))) {
    truth <- interaction_conditions[j, ]
    lp <- cell_mu[[j]][cell_id] + truth$latent_scale * nuisance
    d$y <- as.integer(u < if (truth$link == "logit") plogis(lp) else pnorm(lp))
    for (fitted_link in c("logit", "probit")) {
      for (k in c(1L, 8L)) {
        rows[[length(rows) + 1L]] <- data.frame(replication = rep_id, target = "three_way",
          dgp = truth$link, baseline = truth$baseline, true_pp = truth$true_pp, fitted_link = fitted_link,
          n_covariates = k, fit_once(d, truth, fitted_link, k, cell_design, contrast))
      }
    }
  }
  records[[n_rep + rep_id]] <- do.call(rbind, rows)
}
raw <- do.call(rbind, records)
keys <- c("target", "baseline", "true_pp", "dgp", "fitted_link", "n_covariates")
groups <- split(raw, do.call(interaction, c(raw[keys], list(drop = TRUE))))
summary_rows <- lapply(groups, function(d) {
  usable <- d[d$valid, ]
  err <- usable$estimated_pp - usable$true_pp
  detected <- sum(d$detected)
  ci <- if (is.na(detected)) c(NA_real_, NA_real_) else binom.test(detected, nrow(d))$conf.int
  covered <- usable$coef_covered[!is.na(usable$coef_covered)]
  data.frame(d[1, keys], attempted = nrow(d), valid = nrow(usable),
    failed = nrow(d) - nrow(usable), singular_valid = sum(usable$singular),
    n_detected = detected, detection = detected / nrow(d), detection_lo = ci[1], detection_hi = ci[2],
    bias_pp = if (length(err)) mean(err) else NA_real_, rmse_pp = if (length(err)) sqrt(mean(err^2)) else NA_real_,
    conditional_coef_coverage = if (length(covered)) mean(covered) else NA_real_,
    R_version = as.character(getRversion()), lme4_version = as.character(packageVersion("lme4")))
})
summary_table <- do.call(rbind, summary_rows)
summary_table <- summary_table[do.call(order, summary_table[keys]), ]
stopifnot(nrow(summary_table) == 48, nrow(raw) == 32 * n_rep + 16 * 12,
          all(summary_table$attempted == ifelse(summary_table$target == "main", n_rep, 12L)))
write.csv(raw, "validation/link_fits.csv", row.names = FALSE)
write.csv(summary_table, "validation/link_summary.csv", row.names = FALSE)
print(summary_table[c(keys, "valid", "singular_valid", "n_detected", "rmse_pp")], row.names = FALSE)
