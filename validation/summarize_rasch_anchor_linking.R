# Audit the saved linked banks and compute paired Monte Carlo comparisons.
args <- commandArgs(TRUE)
out <- if (length(args)) args[1L] else "validation"
stage <- if (length(args) > 1L) match.arg(args[2L], c("pilot", "large")) else "pilot"
prefix <- file.path(out, paste0("rasch_anchor_", stage))
cases <- read.csv(paste0(prefix, "_scenarios.csv"))
summ <- read.csv(paste0(prefix, "_summary.csv"))
fits <- read.csv(paste0(prefix, "_fits.csv.gz"))
bank <- read.csv(paste0(prefix, "_new_bank.csv"))$difficulty
rows <- split(seq_len(nrow(fits)), fits$case)
stopifnot(length(rows) == nrow(cases), !anyDuplicated(fits[c("case", "replication")]))

for (k in seq_len(nrow(summ))) {
  s <- summ[k, ]
  f <- fits[rows[[as.character(s$case)]], ]
  centered <- as.matrix(f[paste0("new_centered_", seq_along(bank))])
  shift <- f[[paste0(s$reference, "_shift")]]
  err <- sweep(centered + shift, 2L, bank, "-")
  complete <- rowSums(!is.finite(err)) == 0L
  rmse <- sqrt(rowMeans(err^2))
  worst <- apply(abs(err), 1L, max)
  hit <- complete & abs(shift - s$focal_mean) <= .15 & rmse <= .35 & worst <= .60
  stopifnot(identical(f$replication, seq_len(s$n_rep)), nrow(f) == s$attempted,
            all(hit == f[[paste0(s$reference, "_achieved")]]), sum(hit) == s$achieved,
            sum(complete) == s$complete, sum(f$native_fit) == s$native_fit,
            max(abs(binom.test(sum(hit), nrow(f))$conf.int - c(s$mc_lo, s$mc_hi))) < 1e-12,
            isTRUE(all.equal(mean(rmse, na.rm = TRUE), s$mean_rmse)))
  good <- which(complete)
  if (length(good)) {
    index <- good[1L]
    grid <- seq(-2, 2, by = .25)
    tcc <- function(b) vapply(grid, function(t) sum(plogis(t - b)), numeric(1L))
    score_error <- max(abs(tcc(centered[index, ] + shift[index]) - tcc(bank)))
    stopifnot(abs(score_error - f[[paste0(s$reference, "_tcc_error")]][index]) < 1e-10)
  }
}

pair <- function(left, right, name, ref_left, ref_right = ref_left) {
  a <- fits[rows[[as.character(left$case)]], ]
  b <- fits[rows[[as.character(right$case)]], ]
  stopifnot(identical(a$replication, b$replication))
  dx <- 100 * (a[[paste0(ref_left, "_achieved")]] - b[[paste0(ref_right, "_achieved")]])
  limit <- qnorm(.975) * sd(dx) / sqrt(length(dx))
  data.frame(comparison = name, left = left$case, right = right$case,
    ref_left, ref_right, n = length(dx), difference_pp = mean(dx),
    lo_pp = mean(dx) - limit, hi_pp = mean(dx) + limit)
}
contrasts <- lapply(seq_len(nrow(cases)), function(k) {
  s <- cases[k, ]
  pair(s, s, "estimated_minus_oracle", "estimated", "oracle")
})
core <- subset(cases, phase == "fixed_workload")
for (k in which(core$anchors > 4L)) {
  a <- core[k, ]
  b <- subset(core, n == a$n & band == a$band & focal_mean == a$focal_mean & anchors == a$anchors - 4L)
  stopifnot(nrow(b) == 1L)
  for (ref in c("oracle", "estimated")) {
    contrasts[[length(contrasts) + 1L]] <- pair(a, b, "additional_anchors_fixed_workload", ref)
  }
}
write.csv(do.call(rbind, contrasts), paste0(prefix, "_paired.csv"), row.names = FALSE)
if (stage == "large") {
  tam <- read.csv(paste0(prefix, "_tam.csv"))
  stopifnot(nrow(tam) == 384L, !anyDuplicated(tam[c("case", "replication")]),
            all(table(tam$case) == 64L), all(tam$replication %in% seq_len(64L)))
  common <- tam$native_complete == 1 & tam$tam_complete == 1
  cat(sprintf(paste0("TAM paired checks: %d/%d complete in both engines; %d attainment changes; ",
                     "max linked-item difference = %.5f logits.\n"),
              sum(common), nrow(tam), sum(tam$native_hit != tam$tam_hit),
              max(tam$max_item_difference[common])))
}
cat(sprintf("Checked %d calibration datasets and %d summaries.\n", nrow(fits), nrow(summ)))
