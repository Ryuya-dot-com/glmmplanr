# Check the handoff of one saved linked Rasch bank to mirtCAT's offline mode.
args <- commandArgs(TRUE)
out <- if (length(args)) args[1L] else "validation"
prefix <- file.path(out, "rasch_anchor_large")
cases <- read.csv(paste0(prefix, "_scenarios.csv"))
fits <- read.csv(paste0(prefix, "_fits.csv.gz"))
chosen <- subset(cases, n == 480L & anchors == 8L & band == "broad" &
                   focal_mean == -.75 & phase == "fixed_workload")
stopifnot(nrow(chosen) == 1L)
fit <- subset(fits, case == chosen$case & replication == 1L & oracle_complete == 1L)
stopifnot(nrow(fit) == 1L)
old <- c(seq(-2.2, 2.2, length.out = 12L), seq(1.3, 3.5, length.out = 12L))
truth <- c(old, read.csv(paste0(prefix, "_new_bank.csv"))$difficulty)
new <- unlist(fit[paste0("new_centered_", seq_len(12L))], use.names = FALSE) + fit$oracle_shift
bank <- c(old, new)
model <- mirtCAT::generate.mirt_object(data.frame(a1 = 1, d = -bank), "2PL")
check_theta <- matrix(seq(-2, 2, by = .5))
for (i in seq_along(bank)) {
  p <- mirt::probtrace(mirt::extract.item(model, i), check_theta)[, 2L]
  stopifnot(max(abs(p - plogis(check_theta[, 1L] - bank[i]))) < 1e-12)
}

set.seed(90545001L)
ability <- rnorm(24L)
u <- matrix(runif(length(ability) * length(truth)), nrow = length(ability))
responses <- 1L * (u < plogis(outer(ability, truth, "-")))
results <- mirtCAT::mirtCAT(mo = model, local_pattern = responses, method = "EAP",
  criteria = "MI", start_item = "MI", design = list(min_items = 8L, max_items = 8L))
stopifnot(length(results) == 24L, all(vapply(results, function(person) {
  length(unique(person$items_answered)) == 8L && all(is.finite(person$thetas)) &&
    all(is.finite(person$SE_thetas))
}, logical(1L))))
cat(sprintf("mirtCAT %s, %d items in the linked bank, %d examinees completed an eight-item CAT.\n",
            packageVersion("mirtCAT"), length(bank), length(results)))
