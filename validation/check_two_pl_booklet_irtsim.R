# Independently reproduce the exposure design and point-precision comparison in irtsim.
args <- commandArgs(TRUE)
out <- if (length(args)) args[1L] else "validation"
prefix <- file.path(out, "two_pl_fixed_bank_booklets")
cases <- read.csv(paste0(prefix, "_scenarios.csv"))
wide <- subset(cases, profile == "wide")
truth <- subset(read.csv(paste0(prefix, "_truth.csv")), case == wide$case)
forms <- as.matrix(read.csv(paste0(prefix, "_forms.csv")))
bank <- irtsim::irt_design("2PL", 16L, list(a = truth$a_true, b = truth$b_true))
studies <- list(complete = irtsim::irt_study(bank, sample_sizes = 240L),
  booklets = irtsim::irt_study(bank, sample_sizes = 480L, missing = "booklet",
    test_design = list(booklet_matrix = forms)))
mask <- !is.na(irtsim:::apply_missing(matrix(1L, 480L, 16L), studies$booklets))
expected <- forms[rep(seq_len(nrow(forms)), length.out = 480L), , drop = FALSE] == 1L
stopifnot(identical(unname(mask), unname(expected)), all(colSums(mask) == 240L))
rng <- 91654L
fits <- lapply(names(studies), function(layout) {
  x <- irtsim::irt_simulate(studies[[layout]], iterations = wide$n_rep, seed = rng,
    progress = FALSE, se = FALSE, compute_theta = FALSE)
  s <- summary(x, criterion = "rmse")
  stopifnot(nrow(s$item_summary) == 32L, "rmse" %in% names(s$item_summary))
  tab <- x$item_results
  do.call(rbind, lapply(split(tab, tab$iteration), function(rep) {
    a <- subset(rep, param == "a")
    b <- subset(rep, param == "b")
    stopifnot(identical(a$item, seq_len(16L)), identical(a$item, b$item),
              isTRUE(all.equal(a$true_value, truth$a_true)), isTRUE(all.equal(b$true_value, truth$b_true)))
    direction <- if (all(is.finite(a$estimate)) && sum(a$estimate) < 0) -1 else 1
    slope <- a$estimate * direction
    diff <- b$estimate * direction - b$true_value
    valid <- all(rep$converged) && all(is.finite(c(slope, diff))) && all(slope > 1e-8)
    rmse_a <- if (valid) sqrt(mean((slope - a$true_value)^2)) else NA_real_
    rmse_b <- if (valid) sqrt(mean(diff^2)) else NA_real_
    worst_b <- if (valid) max(abs(diff)) else NA_real_
    b_pass <- valid && rmse_b <= .35 && worst_b <= .60
    data.frame(layout, n = rep$sample_size[1L], replication = rep$iteration[1L],
      converged = all(rep$converged), valid, rmse_a, rmse_b, worst_b,
      attained_b = b_pass, attained_ab = b_pass && rmse_a <= .30)
  }))
})
rows <- do.call(rbind, fits)
write.csv(rows, paste0(prefix, "_irtsim_fits.csv"), row.names = FALSE)
cp <- function(k, n) unname(binom.test(k, n)$conf.int)
s <- do.call(rbind, lapply(fits, function(f) {
  n <- nrow(f)
  ib <- cp(sum(f$attained_b), n)
  ij <- cp(sum(f$attained_ab), n)
  data.frame(layout = f$layout[1L], n = f$n[1L], seed = rng, irtsim = as.character(packageVersion("irtsim")),
    mirt = as.character(packageVersion("mirt")), total = n, converged = sum(f$converged), valid = sum(f$valid),
    med_rmse_b = median(f$rmse_b, na.rm = TRUE), med_rmse_a = median(f$rmse_a, na.rm = TRUE),
    difficulty_pass = sum(f$attained_b), b_lo = ib[1L], b_hi = ib[2L],
    joint_pass = sum(f$attained_ab), joint_lo = ij[1L], joint_hi = ij[2L])
}))
write.csv(s, paste0(prefix, "_irtsim_summary.csv"), row.names = FALSE)
print(s, row.names = FALSE)
