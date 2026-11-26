# Paired external comparisons for the shared-passage design.
library(glmmplanr)
source("validation/few_shared_passages_helpers.R")
args <- commandArgs(TRUE)
out <- if (length(args)) args[1L] else "validation"
stage <- if (length(args) > 1L) match.arg(args[2L], c("pilot", "large")) else "pilot"
workers <- if (length(args) > 2L) as.integer(args[3L]) else 4L
part <- if (length(args) > 3L) match.arg(args[4L], c("population", "bank", "mirt")) else "population"
packages <- switch(part, population = "glmmTMB", bank = c("glmmTMB", "TAM"), mirt = "mirt")
stopifnot(dir.exists(out), workers > 0L)
for (package in packages) {
  if (!requireNamespace(package, quietly = TRUE)) stop("Install ", package, " for the ", part, " comparison.")
}
if (.Platform$OS.type == "windows") workers <- 1L
n_rep <- if (stage == "pilot") 8L else 1024L
seed <- if (stage == "pilot") 91643501L else 91643001L
grid <- function(n, passages, slope_sd, local_sd) {
  expand.grid(n = n, passages = passages, slope_sd = slope_sd, local_sd = local_sd,
              effect_pp = c(0, 10), KEEP.OUT.ATTRS = FALSE)
}
cases <- rbind(grid(48L, c(2L, 5L, 12L), .4, .4),
               grid(48L, c(2L, 5L), 0, .4),
               grid(96L, c(2L, 5L), .4, .4),
               grid(48L, 2L, .4, 0))
cases$case <- seq_len(nrow(cases))
cases$n_rep <- n_rep
cases$seed <- seed
cases$n_bank <- with(cases, ifelse(n == 48L & passages == 2L & slope_sd == .4, min(n_rep, 512L), 0L))
cases <- calibrate_shared_passage_cases(cases)
streams <- glmmplanr:::make_lecuyer_streams(seed, n_rep)
if (part == "mirt") mirt_streams <- glmmplanr:::make_lecuyer_streams(seed + 3299L, min(n_rep, 128L))
max_n <- max(cases$n)
max_p <- max(cases$passages)
stopifnot(nrow(cases) == 16L, sum(cases$n_bank > 0) == 4L)

# Differentiate the same marginal probability functional for both TMB models.
delta <- function(par, X, W, weights, V = NULL) {
  p <- ncol(X)
  vars <- exp(2 * par[-seq_len(p)])
  v <- as.numeric(W %*% vars)
  eta <- as.numeric(X %*% par[seq_len(p)])
  pr <- plogis(sweep(outer(sqrt(v), q$nodes), 1L, eta, "+"))
  est <- sum(weights * as.numeric(pr %*% q$weights))
  if (is.null(V)) return(est)
  d <- pr * (1 - pr)
  g_beta <- crossprod(X, weights * as.numeric(d %*% q$weights))
  g_sd <- as.numeric(crossprod(W, weights * as.numeric(d %*% (q$nodes * q$weights)) / sqrt(v))) * vars
  g <- c(g_beta, g_sd)
  var <- sum(g * as.numeric(V %*% g))
  c(estimate_pp = 100 * est, se_pp = if (is.finite(var) && var > 0) 100 * sqrt(var) else NA_real_)
}

# For the fixed bank, average over readers and reader-by-passage variation.
bank_basis <- function(sched, data) {
  pid <- as.integer(data$passage[match(seq_len(item_n), data$item)])
  contrast <- rep(c(-1, 1), each = item_n) / (sched$n_testlet * rep(tabulate(pid)[pid], 2L))
  new <- data.frame(item = factor(rep(seq_len(item_n), 2L)),
                    passage = factor(rep(pid, 2L)), x = rep(c(-.5, .5), each = item_n))
  item_contrast <- contr.sum(item_n)[as.integer(new$item), , drop = FALSE]
  passage_contrast <- contr.sum(sched$n_testlet)[as.integer(new$passage), , drop = FALSE]
  base_model <- base$X_cells[rep(1:2, each = item_n), , drop = FALSE]
  X <- cbind(base_model, item_contrast, passage_contrast * new$x)
  cache_sched <- sched
  cache_sched$n_testlet <- sched$n_subj * sched$n_testlet
  cache_sched$trials$testlet_id <- as.integer(interaction(data$reader, data$passage, drop = TRUE))
  spec <- random_spec
  spec$dgp$design$Z_testlet_cells <- spec$dgp$design$Z_testlet_cells[, 1L, drop = FALSE]
  spec$dgp$design$Z_item_cells[,] <- 0
  cache <- glmmplanr:::prepare_cpp_schedule_cache(spec, cache_sched)
  idx <- (sched$trials$cell_id - 1L) * item_n + sched$trials$item_id
  cache$X_full <- X[idx, , drop = FALSE]
  spec$dgp$design$X_cells <- X
  for (name in c("Z_subj_cells", "Z_item_cells", "Z_testlet_cells")) {
    spec$dgp$design[[name]] <- spec$dgp$design[[name]][rep(1:2, each = item_n), , drop = FALSE]
  }
  spec$dgp$design$n_cells <- nrow(X)
  spec$contrast_weights <- contrast
  stopifnot(qr(X)$rank == item_n + sched$n_testlet, abs(sum(contrast[contrast > 0]) - 1) < 1e-12)
  list(native = list(spec = spec, cache = cache, structure = "uncorrelated"),
       pid = pid, new = new, contrast = contrast)
}

empty_fit <- function(status) {
  list(valid = FALSE, estimate_pp = NA_real_, se_pp = NA_real_, slope_sd_hat = NA_real_,
       local_sd_hat = NA_real_, boundary = NA, status = status)
}

tmb <- function(data, formula, X, W, weights, slope = FALSE, bank = FALSE) {
  f <- suppressWarnings(glmmTMB::glmmTMB(formula, data = data, family = binomial(),
    control = glmmTMB::glmmTMBControl(rank_check = "stop", parallel = 1L)))
  s <- f$sdr
  stopifnot(identical(names(s$par.fixed), c(rep("beta", ncol(X)), rep("theta", ncol(W)))))
  vc <- glmmTMB::VarCorr(f)$cond
  stopifnot(isTRUE(all.equal(unname(exp(s$par.fixed[-seq_len(ncol(X))])),
    unname(unlist(lapply(vc, attr, "stddev"))), tolerance = 1e-7)))
  point <- delta(s$par.fixed, X, W, weights, s$cov.fixed)
  if (!all(is.finite(point))) return(empty_fit("invalid_delta"))
  boundary <- any(exp(tail(s$par.fixed, ncol(W))) < .001)
  status <- if (f$fit$convergence != 0L) "optimizer" else if (!isTRUE(s$pdHess)) "hessian" else "ok"
  list(valid = status == "ok", estimate_pp = unname(point["estimate_pp"]),
       se_pp = unname(point["se_pp"]),
       slope_sd_hat = if (slope) exp(tail(s$par.fixed, 1L)) else NA_real_,
       local_sd_hat = if (bank) exp(tail(s$par.fixed, 1L)) else NA_real_,
       boundary = boundary, status = status)
}

tam_setup <- function(basis, n) {
  X <- model.matrix(~ 0 + item + x:passage, basis$new)
  A <- array(0, c(nrow(X), 2, ncol(X)),
    dimnames = list(paste0("j", seq_len(nrow(X))), c("0", "1"), paste0("b", seq_len(ncol(X)))))
  A[, 2, ] <- X
  P <- max(basis$pid)
  Q <- cbind(1, diag(P)[rep(basis$pid, 2L), , drop = FALSE])
  resp <- matrix(NA_integer_, n, nrow(X), dimnames = list(NULL, dimnames(A)[[1L]]))
  list(A = A, Q = Q, resp = resp, a = seq_len(n / 2L),
       b = seq.int(n / 2L + 1L, n), X = X, contrast = basis$contrast)
}
equal_local_variance <- function(variance) {
  diag(c(variance[1L, 1L], rep(mean(diag(variance)[-1L]), ncol(variance) - 1L)))
}
tam <- function(y, setup, integration_nodes = 1024L) {
  resp <- setup$resp
  rr <- matrix(y, nrow(resp), item_n)
  resp[setup$a, seq_len(item_n)] <- rr[setup$a, ]
  resp[setup$b, item_n + seq_len(item_n)] <- rr[setup$b, ]
  if (any(colSums(rr) %in% c(0, nrow(rr)))) {
    return(empty_fit("constant_item"))
  }
  max_iter <- 320L
  f <- suppressWarnings(TAM::tam.mml(resp, A = setup$A, Q = setup$Q,
    variance.inits = diag(c(.5, .13, .13)),
    userfct.variance = equal_local_variance, variance.Npars = 2L, item.elim = FALSE, verbose = FALSE,
    control = list(snodes = integration_nodes, QMC = TRUE, maxiter = max_iter,
                   convD = 1e-4, conv = 2e-4, convM = 1e-5)))
  logit <- as.numeric(setup$X %*% f$xsi$xsi)
  stopifnot(max(abs(logit + f$AXsi_[, 2L])) < 1e-7,
            max(abs(f$beta)) < 1e-8, max(abs(f$variance - equal_local_variance(f$variance))) < 1e-10)
  sig <- diag(f$variance)
  point <- 100 * sum(setup$contrast * mean_prob(logit, sqrt(sum(sig[1:2]))))
  status <- if (f$iter >= max_iter) "optimizer" else if (!is.finite(point)) "invalid_estimate" else "ok"
  list(valid = status == "ok", estimate_pp = point, se_pp = NA_real_, slope_sd_hat = NA_real_,
       local_sd_hat = sqrt(sig[2L]), boundary = any(sig < 1e-6), status = status)
}

mirt_fit <- function(y, setup, stream, X, W, weights) {
  resp <- setup$data
  resp[] <- y
  if (any(colSums(resp) %in% c(0, nrow(resp)))) return(empty_fit("constant_item"))
  f <- glmmplanr:::with_rng_stream(stream, suppressWarnings(mirt::mixedmirt(
    resp, covdata = setup$readers, model = 1, fixed = ~0 + items + x:passage,
    random = ~1 | reader:passage, itemdesign = setup$passages, SE = TRUE, verbose = FALSE,
    technical = list(NCYCLES = 1200L))))
  if (!isTRUE(mirt::extract.mirt(f, "converged"))) return(empty_fit("optimizer"))
  p <- mirt::mod2values(f)
  item <- match(paste0("j", seq_len(item_n)), p$item[p$name == "d"])
  all_items <- which(p$name == "d")
  slopes <- match(colnames(X)[-(seq_len(item_n))], p$name)
  idx <- c(all_items[item], slopes, which(p$name == "COV_11" & p$item == "GROUP"), which(p$item == "RANDOM"))
  stopifnot(length(idx) == ncol(X) + ncol(W), !anyNA(idx), all(p$est[idx]))
  theta <- p$value[idx]
  sd_rows <- ncol(X) + seq_len(ncol(W))
  v <- theta[sd_rows]
  if (any(v <= 0)) return(empty_fit("invalid_variance"))
  theta[sd_rows] <- .5 * log(v)
  point <- 100 * delta(theta, X, W, weights)
  V <- mirt::extract.mirt(f, "vcov")
  ids <- vapply(strsplit(rownames(V), ".", fixed = TRUE), function(x) as.integer(x[2L]), integer(1L))
  order <- match(p$parnum[idx], ids)
  status <- "hessian"
  se <- NA_real_
  if (!anyNA(order) && all(is.finite(V))) {
    d <- c(rep(1, ncol(X)), .5 / v)
    V <- V[order, order, drop = FALSE] * outer(d, d)
    eig <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
    if (all(eig > 0)) {
      se <- unname(delta(theta, X, W, weights, V)["se_pp"])
      if (is.finite(se)) status <- "ok"
    }
  }
  list(valid = is.finite(point), estimate_pp = point, se_pp = se, slope_sd_hat = NA_real_,
       local_sd_hat = sqrt(v[2L]), boundary = any(v < 1e-6), status = status)
}

run_block <- function(cfg, reps) {
  sched <- schedule(cfg$n, cfg$passages)
  tr <- sched$trials
  data <- data.frame(x = c(-.5, .5)[tr$cell_id], reader = factor(tr$subj_id),
                     item = factor(tr$item_id), passage = factor(tr$testlet_id))
  cache <- glmmplanr:::prepare_cpp_schedule_cache(random_spec, sched)
  native <- function(re) list(spec = random_spec, cache = cache, structure = re)
  X <- base$X_cells
  W <- matrix(1, 2L, 3L)
  W_slope <- cbind(W, X[, 2L]^2)
  if (part != "population") {
    basis <- bank_basis(sched, data)
    if (part == "bank") tam_design <- tam_setup(basis, cfg$n)
    X_bank <- model.matrix(~ 0 + item + x:passage, basis$new)
    W_bank <- matrix(1, nrow(X_bank), 2L)
    if (part == "mirt") {
      item_labels <- paste0("j", seq_len(item_n))
      mirt_design <- list(data = matrix(0L, cfg$n, item_n, dimnames = list(NULL, item_labels)),
        readers = data[seq_len(cfg$n), c("x", "reader")],
        passages = data.frame(passage = factor(basis$pid), row.names = item_labels))
    }
  }
  do.call(rbind, lapply(reps, function(r) {
    sim <- draw_shared_passage_response(cfg, sched, streams[[r]], max_n, max_p)
    data$y <- sim$y
    evaluate <- function(model, target, fn) {
      clock <- proc.time()
      result <- tryCatch(fn(), error = function(e) empty_fit(paste("error", conditionMessage(e))))
      elapsed <- proc.time() - clock
      row <- modifyList(empty_fit(result$status), result)
      truth <- if (target == "population") cfg$effect_pp else bank_truth
      half_width <- zcrit * row$se_pp
      has_interval <- row$valid && is.finite(half_width)
      data.frame(case = cfg$case, replication = r, model, target, truth_pp = truth, row,
        interval = has_interval,
        covered = if (has_interval) abs(row$estimate_pp - truth) <= half_width else NA,
        detected = has_interval && abs(row$estimate_pp) > half_width,
        cpu_s = sum(elapsed[1:2]))
    }
    if (part == "population") return(do.call(rbind, list(
      evaluate("native_intercept", "population", function() fit_response(sim$y, native("intercept_only"))),
      evaluate("native_slope", "population", function() fit_response(sim$y, native("uncorrelated"))),
      evaluate("tmb_intercept", "population", function() {
        tmb(data, y ~ x + (1 | reader) + (1 | item) + (1 | passage), X, W, c(-1, 1))
      }),
      evaluate("tmb_slope", "population", function() {
        tmb(data, y ~ x + (1 | reader) + (1 | item) + diag(1 + x | passage), X, W_slope, c(-1, 1), slope = TRUE)
      })
    )))
    if (part != "population") {
      pid <- as.integer(basis$new$passage)
      x <- basis$new$x
      logit <- c(cfg$mu_a, cfg$mu_b)[as.integer(x > 0) + 1L] + sim$text[pid] +
        rep(sim$items, 2L) + x * sim$slope[pid]
      bank_truth <- 100 * sum(basis$contrast * mean_prob(logit, sqrt(.7^2 + cfg$local_sd^2)))
      if (part == "mirt") return(evaluate("mirt_bank", "bank", function() {
        mirt_fit(sim$y, mirt_design, mirt_streams[[r]], X_bank, W_bank, basis$contrast)
      }))
      fits <- list(
        evaluate("native_bank", "bank", function() fit_response(sim$y, basis$native)),
        evaluate("tam_bank", "bank", function() tam(sim$y, tam_design)),
        evaluate("tmb_bank", "bank", function() {
          tmb(data, y ~ 0 + item + x:passage + (1 | reader) + (1 | reader:passage),
              X_bank, W_bank, basis$contrast, bank = TRUE)
        })
      )
    }
    do.call(rbind, fits)
  }))
}

prefix <- file.path(out, paste("shared_passage_irt", stage, part, sep = "_"))
checkpoint <- paste0(prefix, "_checkpoints")
dir.create(checkpoint, showWarnings = FALSE)
manifest <- file.path(checkpoint, "scenarios.rds")
if (file.exists(manifest)) stopifnot(identical(readRDS(manifest), cases)) else saveRDS(cases, manifest)
write.csv(cases, paste0(prefix, "_scenarios.csv"), row.names = FALSE)
capture.output(sessionInfo(), sapply(c("glmmplanr", packages), packageVersion),
               file = paste0(prefix, "_environment.txt"))
cat("Package:", normalizePath(find.package("glmmplanr")), "\n")
attempts <- switch(part, population = cases$n_rep, bank = cases$n_bank, mirt = pmin(cases$n_bank, 128L))
blocks <- split(seq_len(max(attempts)), ceiling(seq_len(max(attempts)) / 32L))
jobs <- expand.grid(case = which(attempts > 0), block = seq_along(blocks))
model_count <- switch(part, population = 4L, bank = 3L, mirt = 1L)
started <- proc.time()[["elapsed"]]
paths <- parallel::mclapply(seq_len(nrow(jobs)), function(j) {
  k <- jobs$case[j]; b <- jobs$block[j]
  path <- file.path(checkpoint, sprintf("case_%03d_block_%03d.rds", k, b))
  if (file.exists(path)) return(path)
  rr <- blocks[[b]]
  rows <- run_block(cases[k, ], rr)
  stopifnot(nrow(rows) == model_count * length(rr))
  saveRDS(rows, paste0(path, ".new"))
  stopifnot(file.rename(paste0(path, ".new"), path))
  if (stage == "pilot" || b %% 4L == 0L || b == 1L) {
    cat(sprintf("Case %d block %d/%d, %d/%d valid points\n", k, b, length(blocks), sum(rows$valid), nrow(rows)))
  }
  path
}, mc.cores = workers, mc.set.seed = FALSE, mc.preschedule = FALSE)
success <- vapply(paths, function(x) is.character(x) && length(x) == 1L && file.exists(x), logical(1L))
if (!all(success)) stop(paste(unlist(paths[!success]), collapse = "\n"))
fits <- do.call(rbind, lapply(paths, readRDS))
stopifnot(!anyDuplicated(fits[c("case", "replication", "model")]),
          nrow(fits) == model_count * sum(attempts))
con <- gzfile(paste0(prefix, "_fits.csv.gz"), "wt")
write.csv(fits, con, row.names = FALSE)
close(con)
cat(sprintf("Completed %d datasets, %d/%d valid points, %d intervals, in %.1f seconds.\n",
            sum(attempts), sum(fits$valid), nrow(fits), sum(fits$interval),
            proc.time()[["elapsed"]] - started))
print(subset(as.data.frame(table(fits$model, fits$status)), Freq > 0))
