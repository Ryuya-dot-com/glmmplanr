# Shared Rasch simulations and fit adapters for rasch_fit_benchmark.R.
j_bank <- 16L
items <- sprintf("i%02d", seq_len(j_bank + 1L))
b_bank <- seq(-1.8, 1.8, length.out = j_bank)
focal_fit <- 8L
focal_dif <- 16L
quad <- glmmplanr:::gh_quad_1d(41L)
tam_maxiter <- 1000L

rasch_data <- function(case, stream, forms) {
  assign(".Random.seed", stream, envir = .GlobalEnv)
  standard <- rnorm(480L)
  uniform <- matrix(runif(480L * j_bank), 480L, j_bank)
  ids <- seq_len(case$n)
  g <- rep(0:1, each = case$n / 2L)
  theta <- case$sd * standard[ids] + case$group_shift * g
  b <- b_bank + case$difficulty_shift
  slopes <- rep(1, j_bank)
  slopes[focal_fit] <- case$focal_slope
  eta <- sweep(outer(theta, slopes), 2L, slopes * b, "-")
  eta[, focal_dif] <- eta[, focal_dif] + case$dif * g
  y <- 1L * (uniform[ids, ] < plogis(eta))
  colnames(y) <- items[seq_len(j_bank)]
  if (case$items_per_person < j_bank) {
    present <- forms[rep(seq_len(nrow(forms)), length.out = case$n), , drop = FALSE]
    y[!present] <- NA_integer_
  }
  list(y = y, theta = theta, b = b, g = g)
}

long_rasch <- function(y) {
  index <- which(!is.na(y), arr.ind = TRUE)
  data.frame(id = index[, 1L], item = factor(colnames(y)[index[, 2L]], levels = colnames(y)),
             y = y[index], item_id = index[, 2L])
}

fit_tam <- function(y, g = NULL) {
  f <- TAM::tam.mml(y, Y = if (!is.null(g)) matrix(g, ncol = 1L, dimnames = list(NULL, "g")),
                   constraint = "cases", verbose = FALSE, item.elim = FALSE,
                   control = list(nodes = seq(-6, 6, by = .25), maxiter = tam_maxiter,
                                  conv = 1e-5, convD = 1e-7, convM = 1e-5))
  if (f$iter >= tam_maxiter) stop("TAM iteration limit")
  list(b = f$xsi$xsi, sigma = sqrt(f$variance[1L, 1L]),
       group = if (is.null(g)) 0 else as.numeric(tail(f$beta, 1L)),
       ll = -f$deviance / 2, fit = f)
}

fit_mfrmr <- function(long, n, g = NULL) {
  person <- data.frame(id = seq_len(n), g = if (is.null(g)) 0 else g)
  f <- mfrmr::fit_mfrm(long, person = "id", facets = "item", score = "y",
                       rating_min = 0, rating_max = 1, model = "RSM", method = "MML",
                       mml_engine = "direct", quad_points = 31L,
                       population_formula = if (is.null(g)) ~1 else ~g,
                       person_data = person, person_id = "id",
                       attach_diagnostics = FALSE, progress = FALSE)
  if (inherits(f, "mfrm_fast_fit")) {
    if (!isTRUE(f$converged)) stop("mfrmr C++ population optimizer")
    item <- f$facets[f$facets$Facet == "item", ]
    sigma <- f$population$sigma
  } else {
    if (!isTRUE(f$population$estimation_converged)) stop("mfrmr population optimizer")
    item <- f$facets$others
    sigma <- sqrt(f$population$sigma2)
  }
  co <- f$population$coefficients
  b <- item$Estimate[match(levels(long$item), item$Level)] + f$steps$Estimate[1L] - co[[1L]]
  list(b = b, sigma = sigma, group = if (is.null(g)) 0 else co[["g"]],
       ll = as.numeric(logLik(f)), fit = f)
}

fit_native <- function(y, long, g = NULL) {
  if (is.null(g)) {
    trials <- data.frame(subj_id = long$id, item_id = long$item_id, y = long$y)
    f <- glmmplanr:::fit_rasch_fixed_mml(trials, nrow(y), ncol(y), rep(TRUE, ncol(y)))
    if (!is.na(f$error)) stop(f$error)
    return(list(b = f$b_est, sigma = sqrt(f$person_var), group = 0, ll = NA_real_))
  }
  # Fixed item columns plus group; the installed crossed engine needs an item
  # index. The zero-valued item design contributes nothing to the predictor.
  x <- cbind(model.matrix(~0 + item, long), g = g[long$id])
  n_obs <- nrow(long)
  f <- glmmplanr:::cpp_fit_crossed_glmm(
    y = long$y, X = x, subj_idx = long$id - 1L, Z_subj = matrix(1, n_obs, 1L), n_subj = nrow(y),
    item_idx = rep(0L, n_obs), Z_item = matrix(0, n_obs, 1L), n_item = 1L,
    testlet_idx = integer(), Z_testlet = matrix(0, n_obs, 0L), n_testlet = 0L,
    corr_mode = 0L, regularized = FALSE, max_em_iter = 500L, max_pirls_iter = 28L,
    tol = 2e-5, compute_joint_vcov = FALSE)
  if (!isTRUE(f$converged)) stop(paste("glmmplanr Laplace", f$status))
  list(b = -f$beta[seq_len(ncol(y))], sigma = sqrt(f$Sigma_subj[1L, 1L]),
       group = tail(f$beta, 1L), ll = f$logLik)
}

fit_rasch <- function(engine, y, long, g = NULL) {
  if (any(colSums(y == 1L, na.rm = TRUE) == 0L | colSums(y == 0L, na.rm = TRUE) == 0L)) {
    stop("extreme item")
  }
  f <- switch(engine,
    glmmplanr = fit_native(y, long, g),
    TAM = fit_tam(y, g),
    mfrmr = fit_mfrmr(long, nrow(y), g))
  if (length(f$b) != ncol(y) || any(!is.finite(c(f$b, f$sigma, f$group))) || f$sigma < 0) {
    stop("invalid parameters")
  }
  f
}

# Unique booklet/score combinations supply sufficient statistics for posterior
# moments and Louis observed information under the fitted Rasch model.
rasch_patterns <- function(y, group = rep(0, nrow(y))) {
  mask <- !is.na(y)
  score <- rowSums(y, na.rm = TRUE)
  form <- drop(mask %*% 2^(seq_len(ncol(y)) - 1L)) + 2^ncol(y) * group
  members <- split(seq_len(nrow(y)), form)
  lapply(members, function(id) {
    counts <- table(score[id])
    list(id = id, j = which(mask[id[1L], ]), g = group[id[1L]],
         score = as.numeric(names(counts)), n = as.numeric(counts), raw = score[id])
  })
}

rasch_posterior <- function(patterns, b, sigma, group = 0, moments = FALSE, information = FALSE) {
  theta <- sigma * quad$nodes
  size <- length(b)
  if (moments) eap <- variance <- numeric(sum(vapply(patterns, function(z) length(z$id), integer(1L))))
  if (information) info <- matrix(0, size + 1L, size + 1L)
  ll <- 0
  for (pat in patterns) {
    ability <- theta + group * pat$g
    eta <- outer(ability, -b[pat$j], "+")
    log_mass <- outer(ability, pat$score) + log(quad$weights) +
      rowSums(plogis(eta, lower.tail = FALSE, log.p = TRUE))
    # Log weights have baseline all incorrect. The score contributes ability
    # for each correct response; difficulty constants are added by rasch_ll().
    peak <- apply(log_mass, 2L, max)
    weight <- exp(sweep(log_mass, 2L, peak, "-"))
    weight_sum <- colSums(weight)
    posterior <- sweep(weight, 2L, weight_sum, "/")
    ll <- ll + sum(pat$n * (peak + log(weight_sum)))
    if (moments) {
      posterior_mean <- drop(crossprod(posterior, ability))
      posterior_var <- pmax(0, drop(crossprod(posterior, ability^2)) - posterior_mean^2)
      pos <- match(pat$raw, pat$score)
      eap[pat$id] <- posterior_mean[pos]
      variance[pat$id] <- posterior_var[pos]
    }
    if (information) {
      p <- plogis(eta)
      v <- p * (1 - p)
      k <- length(pat$j)
      ix <- c(pat$j, size + 1L)
      for (s in seq_along(pat$score)) {
        w <- posterior[, s]
        slope_score <- theta * (pat$score[s] - rowSums(p))
        score <- cbind(p, slope_score)
        centered <- sweep(score, 2L, drop(crossprod(w, score)), "-")
        mean_cross <- -drop(crossprod(w * theta, v))
        expected <- diag(c(drop(crossprod(w, v)), sum(w * (theta^2 * rowSums(v) - slope_score))))
        expected[k + 1L, seq_len(k)] <- expected[seq_len(k), k + 1L] <- mean_cross
        info[ix, ix] <- info[ix, ix] + pat$n[s] * (expected - crossprod(centered, centered * w))
      }
    }
  }
  list(ll = ll, eap = if (moments) eap, variance = if (moments) variance,
       information = if (information) info)
}

rasch_ll <- function(y, patterns, b, sigma, group = 0) {
  rasch_posterior(patterns, b, sigma, group)$ll - sum(colSums(y, na.rm = TRUE) * b)
}

positive_ratio <- function(a, b) if (is.finite(b) && b > 0) a / b else NA_real_

common_diagnostics <- function(y, pat, fit, truth) {
  post <- rasch_posterior(pat, fit$b, fit$sigma, moments = TRUE, information = fit$sigma > 0)
  p <- plogis(outer(post$eap, -fit$b, "+"))
  res2 <- (y - p)^2
  v <- p * (1 - p)
  v[is.na(y)] <- NA_real_
  item_infit <- colSums(res2, na.rm = TRUE) / colSums(v, na.rm = TRUE)
  item_outfit <- colMeans(res2 / v, na.rm = TRUE)
  person_infit <- rowSums(res2, na.rm = TRUE) / rowSums(v, na.rm = TRUE)
  person_outfit <- rowMeans(res2 / v, na.rm = TRUE)
  v_person <- var(post$eap)
  s_person <- mean(post$variance)
  cov <- if (!is.null(post$information)) tryCatch(chol2inv(chol(post$information)), error = function(e) NULL)
  s_item <- if (is.null(cov)) NA_real_ else mean(diag(cov)[seq_along(fit$b)])
  v_item <- var(fit$b)
  true_item <- max(0, v_item - s_item)
  delta <- fit$b - truth$b
  data.frame(
    sigma = fit$sigma, b_bias = mean(delta), b_rmse = sqrt(mean(delta^2)), b_max = max(abs(delta)),
    bank_ok = as.integer(sqrt(mean(delta^2)) <= .35 && max(abs(delta)) <= .60),
    theta_rmse = sqrt(mean((post$eap - truth$theta)^2)),
    person_eap_rel = positive_ratio(v_person, v_person + s_person),
    person_eap_ratio = sqrt(positive_ratio(v_person, s_person)),
    item_rel = positive_ratio(true_item, v_item), item_separation = sqrt(positive_ratio(true_item, s_item)),
    mean_item_infit = mean(item_infit), mean_item_outfit = mean(item_outfit),
    mean_person_infit = mean(person_infit), mean_person_outfit = mean(person_outfit),
    focal_infit = item_infit[focal_fit], focal_outfit = item_outfit[focal_fit],
    person_overfit_pct = mean(person_infit < .8), person_underfit_pct = mean(person_infit > 1.2))
}

native_diagnostics <- function(engine, f) {
  if (engine == "TAM") {
    it <- TAM::msq.itemfit(f)$itemfit
    pe <- TAM::tam.personfit(f)
    i <- match(items[focal_fit], it$item)
    return(data.frame(native_item_infit = mean(it$Infit), native_item_outfit = mean(it$Outfit),
                      native_focal_infit = it$Infit[i], native_focal_outfit = it$Outfit[i],
                      native_person_infit = mean(pe$infitPerson), native_person_outfit = mean(pe$outfitPerson),
                      native_person_rel = f$EAP.rel, native_person_separation = NA_real_,
                      native_item_rel = NA_real_, native_item_separation = NA_real_))
  }
  d <- mfrmr::diagnose_mfrm(f, diagnostic_mode = "legacy")
  fast <- inherits(f, "mfrm_fast_fit")
  measures <- if (fast) d$measures else d$fit
  it <- measures[measures$Facet == "item", ]
  pe <- measures[measures$Facet == "Person", ]
  if (!nrow(it) || !nrow(pe)) stop("mfrmr fit tables absent")
  i <- match(items[focal_fit], it$Level)
  item <- d$reliability[d$reliability$Facet == "item", ]
  person <- d$reliability[d$reliability$Facet == "Person", ]
  item_infit <- if (fast) mean(it$Infit) else item$MeanInfit
  item_outfit <- if (fast) mean(it$Outfit) else item$MeanOutfit
  person_infit <- if (fast) mean(pe$Infit) else person$MeanInfit
  person_outfit <- if (fast) mean(pe$Outfit) else person$MeanOutfit
  data.frame(native_item_infit = item_infit, native_item_outfit = item_outfit,
             native_focal_infit = it$Infit[i], native_focal_outfit = it$Outfit[i],
             native_person_infit = person_infit, native_person_outfit = person_outfit,
             native_person_rel = person$ModelReliability, native_person_separation = person$ModelSeparation,
             native_item_rel = item$ModelReliability, native_item_separation = item$ModelSeparation)
}
