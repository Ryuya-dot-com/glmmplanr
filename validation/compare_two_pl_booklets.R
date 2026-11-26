# Pair the booklet trials to the archived complete-presentation trials.
args <- commandArgs(TRUE)
out <- if (length(args)) args[1L] else "validation"
base <- if (length(args) > 1L) args[2L] else "validation"
booklet <- file.path(out, "two_pl_fixed_bank_booklets")
complete <- file.path(base, "two_pl_fixed_bank_medium")
cases <- read.csv(paste0(booklet, "_scenarios.csv"))
archived <- read.csv(paste0(complete, "_scenarios.csv"))
truth <- read.csv(paste0(booklet, "_truth.csv"))
truth_archive <- read.csv(paste0(complete, "_truth.csv"))
forms <- as.matrix(read.csv(paste0(booklet, "_forms.csv")))
a <- read.csv(paste0(booklet, "_fits.csv.gz"))
b <- read.csv(paste0(complete, "_fits.csv.gz"))
keys <- c("case", "replication", "model")
stopifnot(!anyDuplicated(a[keys]), !anyDuplicated(b[keys]),
          nrow(a) == 4L * sum(cases$n_rep), identical(dim(forms), c(30L, 16L)),
          all(forms %in% 0:1), all(rowSums(forms) == 8L))
forms480 <- forms[rep(seq_len(nrow(forms)), length.out = 480L), ]
cross <- crossprod(forms480)
stopifnot(sum(cross[lower.tri(cross)]) == 13440L,
          all(diag(cross) == 240L), all(cross[lower.tri(cross)] == 112L))
paired_rate <- function(diff) {
  pos <- binom.test(sum(diff > 0), length(diff), conf.level = .975)$conf.int
  neg <- binom.test(sum(diff < 0), length(diff), conf.level = .975)$conf.int
  c(difference = 100 * mean(diff), lo = 100 * (pos[1L] - neg[2L]),
    hi = 100 * (pos[2L] - neg[1L]))
}
paired_median <- function(diff) {
  if (!length(diff)) return(c(difference = NA_real_, lo = NA_real_, hi = NA_real_))
  d <- sort(diff)
  k <- qbinom(.025, length(d), .5)
  c(difference = median(d), lo = if (k) d[k] else -Inf,
    hi = if (k) d[length(d) - k + 1L] else Inf)
}

results <- lapply(seq_len(nrow(cases)), function(k) {
  cfg <- cases[k, ]
  old <- archived[match(cfg$complete_case, archived$case), ]
  stopifnot(cfg$n == 480L, old$n == 240L, cfg$items == 16L, old$items == 16L,
            cfg$n_rep == old$n_rep, cfg$seed == old$seed, cfg$profile == old$profile)
  t_args <- c("item", "a_true", "b_true", "d_true")
  stopifnot(isTRUE(all.equal(truth[truth$case == cfg$case, t_args],
    truth_archive[truth_archive$case == cfg$complete_case, t_args], check.attributes = FALSE)))
  do.call(rbind, lapply(unique(a$model), function(m) {
    incomplete <- a[a$case == k & a$model == m, ]
    full <- b[b$case == cfg$complete_case & b$model == m, ]
    stopifnot(identical(incomplete$replication, seq_len(cfg$n_rep)),
              identical(incomplete$replication, full$replication),
              identical(incomplete$gap_true_pp, full$gap_true_pp))
    valid <- incomplete$valid & full$valid
    contrasts <- rbind(valid_pp = paired_rate(incomplete$valid - full$valid),
      difficulty_pass_pp = paired_rate(incomplete$attained_b - full$attained_b))
    if (!all(is.na(incomplete$attained_ab))) {
      contrasts <- rbind(contrasts, joint_pass_pp = paired_rate(incomplete$attained_ab - full$attained_ab))
    }
    for (field in c("rmse_a", "rmse_b", "curve_rmse_pp")) {
      z <- paired_median(incomplete[[field]][valid] - full[[field]][valid])
      contrasts <- rbind(contrasts, matrix(z, nrow = 1L, dimnames = list(field, names(z))))
    }
    cbind(cfg, model = m, paired_valid = sum(valid), measure = rownames(contrasts), contrasts, row.names = NULL)
  }))
})
write.csv(do.call(rbind, results), paste0(booklet, "_layout_differences.csv"), row.names = FALSE)
cat(sprintf("Paired %d booklet trials to the complete-presentation archive and verified the common banks and exposure.\n",
            nrow(a)))
