# Recompute the item-based criteria from fitted banks, then compare matched response datasets.
args <- commandArgs(TRUE)
out <- if (length(args)) args[1L] else "validation"
stage <- if (length(args) > 1L) match.arg(args[2L], c("pilot", "medium", "booklets")) else "medium"
prefix <- file.path(out, paste("two_pl_fixed_bank", stage, sep = "_"))
cases <- read.csv(paste0(prefix, "_scenarios.csv"))
truth <- read.csv(paste0(prefix, "_truth.csv"))
fits <- read.csv(paste0(prefix, "_fits.csv.gz"))
items <- read.csv(paste0(prefix, "_items.csv.gz"))
key_cols <- c("case", "replication", "model")
key <- function(x) paste(x$case, x$replication, x$model)
free <- c("glmmTMB_rr", "tam_2pl", "mirt_2pl")
models <- c("glmmplanr_rasch", "mirt_rasch", free, "glmer_known")
if (stage == "booklets") models <- setdiff(models, c("glmmTMB_rr", "glmer_known"))
n_bank <- sum(cases$n_rep * cases$items)
stopifnot(!anyDuplicated(fits[key_cols]), !anyDuplicated(items[c(key_cols, "item")]),
          identical(sort(unique(fits$model)), sort(models)),
          nrow(fits) == length(models) * sum(cases$n_rep), nrow(items) == length(models) * n_bank)
expected <- paste(rep(cases$case, cases$items), sequence(cases$items))
stopifnot(identical(paste(truth$case, truth$item), expected))
i_pos <- split(seq_len(nrow(items)), key(items))
selected <- seq_len(nrow(fits))
fit_key <- key(fits)
stopifnot(setequal(names(i_pos), fit_key))
truth_key <- paste(truth$case, truth$item)
j_pos <- match(paste(items$case, items$item), truth_key)
stopifnot(!anyNA(j_pos), all(abs(truth$d_true + truth$a_true * truth$b_true) < 1e-12))
check <- function(computed, stored) stopifnot(isTRUE(all.equal(computed, stored, tolerance = 1e-10)))
rmse <- function(x) sqrt(mean(x^2))
fits$max_a <- NA_real_
for (r in selected) {
  f <- fits[r, ]
  pos <- i_pos[[fit_key[r]]]
  b <- items[pos, ]
  t <- truth[j_pos[pos], ]
  stopifnot(identical(b$item, seq_len(cases$items[f$case])), nrow(t) == cases$items[f$case])
  if (f$valid) {
    stopifnot(f$status == "ok", all(b$a > 1e-8))
    fits$max_a[r] <- max(b$a)
    a_err <- b$a - t$a_true
    b_diff <- -b$d / b$a
    b_err <- b_diff - t$b_true
    shared <- b_err[seq_len(8L)]
    check(b_diff, b$b)
    check(c(rmse(a_err), rmse(b_err), max(abs(b_err)), rmse(shared), max(abs(shared))),
          c(f$rmse_a, f$rmse_b, f$max_b, f$rmse_b_shared, f$max_b_shared))
    stopifnot(f$attained_b == (rmse(b_err) <= .35 && max(abs(b_err)) <= .60),
              f$attained_b_shared == (rmse(shared) <= .35 && max(abs(shared)) <= .60))
    if (f$model %in% free) stopifnot(f$attained_ab == (f$attained_b && rmse(a_err) <= .30))
  } else {
    stopifnot(!f$attained_b, !f$attained_b_shared, all(is.na(b$b)))
    if (f$model %in% free) stopifnot(!f$attained_ab)
  }
  if (!f$model %in% free) stopifnot(is.na(f$attained_ab))
}

cp <- function(k, n) if (is.na(k)) c(NA_real_, NA_real_) else unname(binom.test(k, n)$conf.int)
gate <- function(ci) {
  if (is.na(ci[1L])) return("not_applicable")
  if (ci[1L] >= .80) return("pass")
  if (ci[2L] < .80) return("fail")
  "indeterminate"
}
groups <- split(fits, paste(fits$case, fits$model))
summaries <- lapply(groups, function(f) {
  cfg <- cases[f$case[1L], ]
  stopifnot(identical(f$replication, seq_len(cfg$n_rep)))
  n <- nrow(f)
  v <- subset(f, valid)
  vb <- sum(f$attained_b)
  va <- if (all(is.na(f$attained_ab))) NA_integer_ else sum(f$attained_ab)
  shared <- sum(f$attained_b_shared)
  i <- cp(nrow(v), n)
  ib <- cp(vb, n)
  ia <- cp(va, n)
  is <- cp(shared, n)
  cbind(cfg, model = f$model[1L], total = n, valid = nrow(v), valid_lo = i[1L], valid_hi = i[2L],
        pd_hessian = if (all(is.na(f$pd_hessian))) NA_integer_ else sum(f$pd_hessian, na.rm = TRUE),
        b_attained = vb, b_lo = ib[1L], b_hi = ib[2L], b_gate = gate(ib),
        ab_attained = va, ab_lo = ia[1L], ab_hi = ia[2L], ab_gate = gate(ia),
        shared_attained = shared, shared_lo = is[1L], shared_hi = is[2L],
        med_rmse_a = median(v$rmse_a), med_rmse_b = median(v$rmse_b), max_a_gt_5 = sum(v$max_a > 5),
        med_max_b = median(v$max_b), med_rmse_b_shared = median(v$rmse_b_shared),
        med_curve_rmse_pp = median(v$curve_rmse_pp),
        gap_true_pp = f$gap_true_pp[1L], gap_bias_pp = mean(v$gap_hat_pp - v$gap_true_pp),
        gap_rmse_pp = rmse(v$gap_hat_pp - v$gap_true_pp),
        median_cpu_s = median(f$cpu_s), total_cpu_s = sum(f$cpu_s), row.names = NULL)
})
paired_ci <- function(diff) {
  n <- length(diff)
  if (!n) return(c(difference = NA_real_, lo = NA_real_, hi = NA_real_))
  pos <- binom.test(sum(diff > 0), n, conf.level = .975)$conf.int
  neg <- binom.test(sum(diff < 0), n, conf.level = .975)$conf.int
  c(difference = 100 * mean(diff), lo = 100 * (pos[1L] - neg[2L]), hi = 100 * (pos[2L] - neg[1L]))
}
mean_ci <- function(diff) {
  m <- mean(diff)
  e <- qnorm(.975) * sd(diff) / sqrt(length(diff))
  c(difference = m, lo = m - e, hi = m + e)
}
paired <- lapply(cases$case, function(k) {
  cfg <- cases[k, ]
  compare <- function(left, right) {
    a <- groups[[paste(k, left)]]
    b <- groups[[paste(k, right)]]
    stopifnot(identical(a$replication, b$replication), identical(a$gap_true_pp, b$gap_true_pp))
    v <- a$valid & b$valid
    q_usable <- v & a$max_a <= 5 & b$max_a <= 5
    delta <- rbind(
      valid_pp = paired_ci(a$valid - b$valid),
      b_attained_pp = paired_ci(a$attained_b - b$attained_b),
      shared_attained_pp = paired_ci(a$attained_b_shared - b$attained_b_shared),
      curve_rmse_pp = mean_ci(a$curve_rmse_pp[v] - b$curve_rmse_pp[v]),
      gh_loglik_per_person = mean_ci((a$ref_loglik[q_usable] - b$ref_loglik[q_usable]) / cfg$n),
      score_gap_pp = mean_ci(a$gap_hat_pp[v] - b$gap_hat_pp[v]),
      rmse_a = mean_ci(a$rmse_a[v] - b$rmse_a[v]),
      rmse_b = mean_ci(a$rmse_b[v] - b$rmse_b[v]))
    if (left %in% free && right %in% free) {
      delta <- rbind(delta, ab_attained_pp = paired_ci(a$attained_ab - b$attained_ab))
    }
    cbind(cfg, left, right, paired_valid = sum(v), gh_loglik_pairs = sum(q_usable),
          measure = rownames(delta), delta, row.names = NULL)
  }
  common <- rbind(compare("glmmplanr_rasch", "mirt_rasch"), compare("tam_2pl", "mirt_2pl"))
  if (stage != "booklets") {
    common <- rbind(common, compare("glmmTMB_rr", "mirt_2pl"), compare("glmer_known", "mirt_2pl"))
  }
  rbind(common, compare("mirt_2pl", "glmmplanr_rasch"))
})
write.csv(do.call(rbind, summaries), paste0(prefix, "_summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, paired), paste0(prefix, "_paired.csv"), row.names = FALSE)
cat(sprintf("Verified %d scenario-model summaries, %d fits, and %d item estimates.\n",
            length(summaries), nrow(fits), nrow(items)))
