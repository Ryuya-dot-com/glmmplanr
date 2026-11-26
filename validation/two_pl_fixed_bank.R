# Paired item-calibration benchmark for three free-loading and two common-loading fits.
library(glmmplanr)
args <- commandArgs(TRUE)
out <- if (length(args)) args[1L] else "validation"
stage <- if (length(args) > 1L) match.arg(args[2L], c("pilot", "medium", "booklets")) else "pilot"
workers <- if (length(args) > 2L) as.integer(args[3L]) else 4L
stopifnot(dir.exists(out), workers > 0L)
if (.Platform$OS.type == "windows") workers <- 1L
booklet_run <- stage == "booklets"
packages <- c("glmmplanr", "TAM", "mirt", if (!booklet_run) c("glmmTMB", "lme4"))
stopifnot(all(vapply(packages, requireNamespace, logical(1L), quietly = TRUE)))
n_rep <- if (stage == "pilot") 8L else 128L
seed <- if (stage == "pilot") 91651001L else 91652001L
streams <- glmmplanr:::make_lecuyer_streams(seed, n_rep)
rr_streams <- glmmplanr:::make_lecuyer_streams(seed + 100L, n_rep)
layouts <- data.frame(n = c(120L, 240L, 240L, 480L), items = c(16L, 8L, 16L, 8L))
profiles <- list(common = c(1, 1), moderate = c(.75, 1.35), wide = c(.55, 1.8))
cases <- do.call(rbind, lapply(names(profiles), function(profile) cbind(layouts, profile)))
cases$case <- seq_len(nrow(cases))
cases$n_rep <- n_rep
cases$seed <- seed
max_n <- max(cases$n)
max_j <- max(cases$items)
stopifnot(nrow(cases) == 12L, max_j == 16L)
if (booklet_run) {
  cases <- subset(cases, n == 240L & items == max_j)
  cases$complete_case <- cases$case
  cases$case <- seq_len(nrow(cases))
  cases$n <- max_n
  cases$items_per_subject <- 8L
  # The two affine hyperplanes for each nonzero normal form complementary booklets.
  bits <- vapply(0L:3L, function(b) bitwAnd(bitwShiftR(0L:15L, b), 1L), integer(max_j))
  side <- (bits[-1L, ] %*% t(bits)) %% 2L
  forms <- 1L * rbind(side == 0L, side == 1L)
  colnames(forms) <- paste0("j", seq_len(max_j))
  present <- forms[rep(seq_len(nrow(forms)), length.out = max_n), , drop = FALSE] == 1L
  joint <- crossprod(present)
  stopifnot(nrow(forms) == 30L, all(rowSums(present) == 8L), all(diag(joint) == 240L),
            all(joint[lower.tri(joint)] == 112L), sum(present) == 3840L)
}
q <- glmmplanr:::gh_quad_1d(61L)
b_bank <- c(-1.8, -1.2, -.6, -.2, .2, .6, 1.2, 1.8, -1.6, -1, -.5, -.3, .3, .5, 1, 1.6)
loading_rank <- c(4, 1, 3, 2, 2, 3, 1, 4, 1, 3, 2, 4, 4, 2, 3, 1)
banks <- do.call(rbind, lapply(seq_len(nrow(cases)), function(k) {
  cfg <- cases[k, ]
  j <- seq_len(cfg$items)
  ends <- profiles[[cfg$profile]]
  a <- exp(seq(log(ends[1L]), log(ends[2L]), length.out = 4L))[loading_rank[j]]
  data.frame(case = k, item = j, a_true = a, b_true = b_bank[j], d_true = -a * b_bank[j])
}))

logits <- function(theta, a, d) sweep(outer(theta, a), 2L, d, "+")
predict_bank <- function(theta, a, d) plogis(logits(theta, a, d))
score_gap <- function(a, d) {
  100 * sum(q$weights * rowMeans(predict_bank(q$nodes + .4, a, d) - predict_bank(q$nodes, a, d)))
}
reference_loglik <- function(y, a, d) {
  eta <- logits(q$nodes, a, d)
  observed_correct <- !is.na(y) & y == 1L
  observed_wrong <- !is.na(y) & y == 0L
  loglik <- sweep(observed_correct %*% t(plogis(eta, log.p = TRUE)) +
                    observed_wrong %*% t(plogis(eta, lower.tail = FALSE, log.p = TRUE)),
                  2L, log(q$weights), "+")
  m <- apply(loglik, 1L, max)
  sum(m + log(rowSums(exp(loglik - m))))
}

rasch_native <- function(long, y) {
  f <- glmmplanr:::fit_rasch_fixed_mml(long, nrow(y), ncol(y), rep(TRUE, ncol(y)))
  if (!is.na(f$error)) return(list(status = "optimizer"))
  list(a = rep(sqrt(f$person_var), ncol(y)), d = -f$b_est, status = "ok")
}
mirt_model <- function(y, type) {
  f <- mirt::mirt(y, 1, itemtype = type, method = "EM", quadpts = 49L, TOL = 1e-5,
    verbose = FALSE, technical = list(NCYCLES = 1500L))
  vals <- mirt::mod2values(f)
  a <- vals$value[vals$name == "a1"]
  d <- vals$value[vals$name == "d"]
  sig <- vals$value[vals$name == "COV_11"]
  stopifnot(length(a) == ncol(y), length(d) == ncol(y), length(sig) == 1L,
            abs(vals$value[vals$name == "MEAN_1"]) < 1e-8)
  list(a = a * sqrt(sig), d = d, status = if (isTRUE(mirt::extract.mirt(f, "converged"))) "ok" else "optimizer")
}
tam_model <- function(y) {
  maxiter <- 1000L
  f <- TAM::tam.mml.2pl(y, irtmodel = "2PL", est.variance = FALSE, item.elim = FALSE,
    verbose = FALSE, control = list(nodes = seq(-6, 6, by = .25), maxiter = maxiter,
      convD = 1e-5, conv = 1e-5, convM = 1e-5))
  a <- f$B[, 2L, 1L] - f$B[, 1L, 1L]
  d <- f$AXsi[, 2L] - f$AXsi[, 1L]
  stopifnot(max(abs(f$beta)) < 1e-8, abs(f$variance - 1) < 1e-8,
            isTRUE(all.equal(unname(a), f$item_irt$alpha)),
            max(abs(-a * f$item_irt$beta - d)) < 1e-8)
  list(a = a, d = d, status = if (f$iter >= maxiter) "optimizer" else "ok")
}
rr_model <- function(long, stream) {
  f <- glmmplanr:::with_rng_stream(stream, glmmTMB::glmmTMB(y ~ 0 + item + rr(0 + item | person, d = 1),
    data = long, family = binomial(),
    control = glmmTMB::glmmTMBControl(start_method = list(method = "res"), parallel = 1L)))
  a <- drop(f$obj$env$report(f$fit$parfull)$fact_load[[1L]])
  list(a = a, d = unname(glmmTMB::fixef(f)$cond),
       status = if (f$fit$convergence == 0L) "ok" else "optimizer", pd_hessian = isTRUE(f$sdr$pdHess))
}
known_model <- function(long, a) {
  f <- lme4::glmer(y ~ 0 + item + (0 + a_known | person), data = long,
                  family = binomial(), nAGQ = 1L,
                  control = lme4::glmerControl(optimizer = "bobyqa"))
  sig <- as.numeric(attr(lme4::VarCorr(f)$person, "stddev"))
  conv <- f@optinfo$conv
  list(a = a * sig, d = unname(lme4::fixef(f)),
       status = if (conv$opt == 0L && is.null(conv$lme4$messages)) "ok" else
         if (sig <= 1e-8) "zero_loading" else "optimizer")
}

model_names <- c("glmmplanr_rasch", "mirt_rasch", "mirt_2pl", "tam_2pl", "glmmTMB_rr", "glmer_known")
if (booklet_run) model_names <- model_names[seq_len(4L)]
free_models <- c("mirt_2pl", "tam_2pl", "glmmTMB_rr")
simulate <- function(k, r) {
  cfg <- cases[k, ]
  bank <- banks[banks$case == k, ]
  n <- cfg$n
  j <- cfg$items
  labels <- paste0("j", seq_len(j))
  draw <- glmmplanr:::with_rng_stream(streams[[r]], list(theta = rnorm(max_n),
    u = matrix(runif(max_n * max_j), max_n, max_j)))
  theta <- draw$theta[seq_len(n)]
  y <- 1L * (draw$u[seq_len(n), seq_len(j)] < predict_bank(theta, bank$a_true, bank$d_true))
  colnames(y) <- labels
  if (booklet_run) y[!present] <- NA_integer_
  long <- expand.grid(subj_id = seq_len(n), item_id = seq_len(j), KEEP.OUT.ATTRS = FALSE)
  long$person <- factor(long$subj_id)
  long$item <- factor(labels[long$item_id], levels = labels)
  long$a_known <- bank$a_true[long$item_id]
  long$y <- as.vector(y)
  if (booklet_run) long <- long[!is.na(long$y), , drop = FALSE]
  correct <- colSums(y, na.rm = TRUE)
  constant <- sum(correct == 0L | correct == colSums(!is.na(y)))
  base_curve <- predict_bank(q$nodes, bank$a_true, bank$d_true)
  true_gap <- score_gap(bank$a_true, bank$d_true)
  fn <- list(function() rasch_native(long, y), function() mirt_model(y, "Rasch"),
    function() mirt_model(y, "2PL"), function() tam_model(y), function() rr_model(long, rr_streams[[r]]),
    function() known_model(long, bank$a_true))
  results <- lapply(seq_along(model_names), function(m) {
    method <- model_names[m]
    clock <- proc.time()
    f <- if (constant) list(status = "constant_item") else tryCatch(suppressMessages(suppressWarnings(fn[[m]]())),
      error = function(e) list(status = paste("error:", conditionMessage(e))))
    elapsed <- proc.time() - clock
    a <- if (is.null(f$a)) rep(NA_real_, j) else unname(f$a)
    d <- if (is.null(f$d)) rep(NA_real_, j) else unname(f$d)
    stopifnot(length(a) == j, length(d) == j)
    # The one-factor model can be reflected without changing its response likelihood.
    if (all(is.finite(a)) && sum(a) < 0) a <- -a
    status <- f$status
    if (status == "ok" && !all(is.finite(c(a, d)))) status <- "nonfinite"
    if (status == "ok" && any(a <= 1e-8)) status <- "nonpositive_loading"
    valid <- status == "ok"
    b <- if (valid) -d / a else rep(NA_real_, j)
    ea <- a - bank$a_true
    eb <- b - bank$b_true
    ra <- if (valid) sqrt(mean(ea^2)) else NA_real_
    rb <- sqrt(mean(eb^2))
    max_b <- max(abs(eb))
    shared_err <- eb[seq_len(8L)]
    rb_shared <- sqrt(mean(shared_err^2))
    max_b_shared <- max(abs(shared_err))
    attained <- valid && rb <= .35 && max_b <= .60
    curve <- if (valid) predict_bank(q$nodes, a, d) else matrix(NA_real_, length(q$nodes), j)
    c_error_sq <- rowMeans((curve - base_curve)^2)
    est_gap <- if (valid) score_gap(a, d) else NA_real_
    items <- data.frame(case = k, replication = r, model = method, item = seq_len(j),
                        a = a, d = d, b = b)
    metrics <- data.frame(case = k, replication = r, model = method, valid = valid, status = status,
      pd_hessian = if (is.null(f$pd_hessian)) NA else f$pd_hessian, constant_items = constant,
      rmse_a = ra, rmse_b = rb, max_b = max_b, rmse_b_shared = rb_shared, max_b_shared = max_b_shared,
      curve_rmse_pp = 100 * sqrt(sum(q$weights * c_error_sq)),
      gap_true_pp = true_gap, gap_hat_pp = est_gap,
      ref_loglik = if (valid) reference_loglik(y, a, d) else NA_real_,
      attained_b = attained, attained_b_shared = valid && rb_shared <= .35 && max_b_shared <= .60,
      attained_ab = if (method %in% free_models) attained && ra <= .30 else NA,
      cpu_s = sum(elapsed[1:2]))
    list(items = items, metrics = metrics)
  })
  list(items = do.call(rbind, lapply(results, "[[", "items")),
       fits = do.call(rbind, lapply(results, "[[", "metrics")))
}

prefix <- file.path(out, paste("two_pl_fixed_bank", stage, sep = "_"))
checkpoint <- paste0(prefix, "_checkpoints")
dir.create(checkpoint, showWarnings = FALSE)
manifest <- file.path(checkpoint, "scenarios.rds")
if (file.exists(manifest)) stopifnot(identical(readRDS(manifest), list(cases, banks))) else
  saveRDS(list(cases, banks), manifest)
write.csv(cases, paste0(prefix, "_scenarios.csv"), row.names = FALSE)
write.csv(banks, paste0(prefix, "_truth.csv"), row.names = FALSE)
if (booklet_run) write.csv(forms, paste0(prefix, "_forms.csv"), row.names = FALSE)
capture.output(sessionInfo(), sapply(packages, packageVersion), file = paste0(prefix, "_environment.txt"))
cat("Package:", normalizePath(find.package("glmmplanr")), "\n")
blocks <- split(seq_len(n_rep), ceiling(seq_len(n_rep) / 16L))
jobs <- expand.grid(case = cases$case, block = seq_along(blocks))
started <- proc.time()[["elapsed"]]
paths <- parallel::mclapply(seq_len(nrow(jobs)), function(i) {
  k <- jobs$case[i]
  block <- jobs$block[i]
  path <- file.path(checkpoint, sprintf("case_%03d_block_%03d.rds", k, block))
  if (file.exists(path)) return(path)
  results <- lapply(blocks[[block]], function(r) simulate(k, r))
  fits <- do.call(rbind, lapply(results, "[[", "fits"))
  items <- do.call(rbind, lapply(results, "[[", "items"))
  stopifnot(nrow(fits) == length(model_names) * length(blocks[[block]]),
            nrow(items) == nrow(fits) * cases$items[k])
  saveRDS(list(fits = fits, items = items), paste0(path, ".part"))
  stopifnot(file.rename(paste0(path, ".part"), path))
  cat(sprintf("Case %d block %d/%d, %d/%d valid\n", k, block, length(blocks), sum(fits$valid), nrow(fits)))
  path
}, mc.cores = workers, mc.set.seed = FALSE, mc.preschedule = FALSE)
completed <- vapply(paths, function(p) is.character(p) && length(p) == 1L && file.exists(p), logical(1L))
if (!all(completed)) stop(paste(unlist(paths[!completed]), collapse = "\n"))
all_blocks <- lapply(paths, readRDS)
for (kind in c("fits", "items")) {
  tab <- do.call(rbind, lapply(all_blocks, "[[", kind))
  con <- gzfile(paste0(prefix, "_", kind, ".csv.gz"), "wt")
  write.csv(tab, con, row.names = FALSE)
  close(con)
  if (kind == "fits") print(subset(as.data.frame(table(tab$model, tab$status)), Freq > 0))
}
cat(sprintf("Finished %d datasets and %d fits in %.2f seconds.\n",
            nrow(cases) * n_rep, nrow(cases) * n_rep * length(model_names), proc.time()[["elapsed"]] - started))
