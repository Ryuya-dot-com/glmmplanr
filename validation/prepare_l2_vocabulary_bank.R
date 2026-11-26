# Derive reusable item parameters from the authors' public response workbook.
# Raw responses and test texts stay outside the repository; no scoring claims.
library(readxl)
library(mirt)
library(glmmplanr)
args <- commandArgs(TRUE)
input <- if (length(args)) args[1L] else "/private/tmp/glmmplanr-hamada-vst.xlsx"
out <- if (length(args) > 1L) args[2L] else "validation/l2_empirical_bank"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
stopifnot(!file.exists(file.path(out, "items.csv")))
set.seed(96710001L)
item_rows <- list(); summaries <- list(); masks <- list(); audits <- list()
recall_bands <- list(A = c(3L, 4L), B = c(1L, 3L), C = c(1L, 5L),
                     D = c(2L, 5L), E = c(2L, 4L))
for (format in c("Recognition", "Recall")) {
  raw <- read_excel(input, sheet = format, na = c("", "NA"))
  y <- vapply(raw[-c(1L, 2L)], function(z) {
    z <- trimws(as.character(z), whitespace = "[\\h\\v]")
    z[z == ""] <- NA_character_
    stopifnot(all(is.na(z) | z %in% c("0", "1")))
    as.numeric(z)
  }, numeric(nrow(raw)))
  stopifnot(nrow(y) == 1220L, ncol(y) == 50L)
  bands <- as.integer(sub("L([0-9]+)_.*", "\\1", colnames(y)))
  mask <- sapply(recall_bands, function(b) if (format == "Recall") bands %in% b else !bands %in% b)
  # A single unanswered item must not be mistaken for a different booklet.
  assignment <- apply(!is.na(y), 1L, function(present) {
    fits <- which(colSums(present & !mask) == 0L)
    stopifnot(length(fits) == 1L); fits
  })
  expected <- t(mask[, assignment, drop = FALSE])
  stopifnot(all(is.na(y) | expected))
  masks[[format]] <- data.frame(format, booklet = rep(names(recall_bands), each = 50L),
    item_id = rep(colnames(y), 5L), present = as.integer(mask),
    observed_booklet_n = rep(tabulate(assignment, nbins = 5L), each = 50L))
  for (model in c("2PL", "Rasch")) {
    fit <- mirt(y, 1, itemtype = model, verbose = FALSE,
                technical = list(NCYCLES = 2000L))
    ip <- coef(fit, IRTpars = TRUE, simplify = TRUE)
    converged <- isTRUE(extract.mirt(fit, "converged"))
    stopifnot(converged, all(is.finite(ip$items[, c("a", "b")])) )
    theta <- matrix(seq(-4, 4, length.out = 41L))
    error <- vapply(seq_len(ncol(y)), function(i)
      max(abs(probtrace(extract.item(fit, i), theta)[, 2L] -
        plogis(ip$items[i, "a"] * (theta[, 1L] - ip$items[i, "b"])))), 0.)
    stopifnot(max(error) < 1e-12)
    key <- paste(format, model, sep = "_")
    item_rows[[key]] <- data.frame(format, model, item_id = colnames(y), frequency_band = bands,
      a = ip$items[, "a"], b = ip$items[, "b"],
      observed = colSums(!is.na(y)), correct = colSums(y, na.rm = TRUE))
    summaries[[key]] <- data.frame(format, model, participants = nrow(y), items = ncol(y),
      responses_per_person = if (format == "Recall") 20L else 30L,
      incidental_missing = sum(expected & is.na(y)),
      exposure_min = min(colSums(!is.na(y))), exposure_max = max(colSums(!is.na(y))),
      converged, person_var = ip$cov[1L, 1L],
      a_min = min(ip$items[, "a"]), a_max = max(ip$items[, "a"]),
      b_min = min(ip$items[, "b"]), b_max = max(ip$items[, "b"]),
      trace_max_abs_error = max(error), logLik = as.numeric(logLik(fit)))
    if (model == "Rasch") {
      idx <- which(!is.na(y), arr.ind = TRUE)
      tr <- data.frame(subj_id = idx[, 1L], item_id = idx[, 2L], y = y[idx])
      native <- glmmplanr:::fit_rasch_fixed_mml(tr, nrow(y), ncol(y), rep(TRUE, ncol(y)))
      audits[[format]] <- data.frame(format, native_error = native$error,
        native_person_var = native$person_var, mirt_person_var = ip$cov[1L, 1L],
        max_b_difference = max(abs(native$b_est - ip$items[, "b"])))
    }
  }
}
write.csv(do.call(rbind, item_rows), file.path(out, "items.csv"), row.names = FALSE)
write.csv(do.call(rbind, summaries), file.path(out, "summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, masks), file.path(out, "booklets.csv"), row.names = FALSE)
write.csv(do.call(rbind, audits), file.path(out, "native_rasch_audit.csv"), row.names = FALSE)
saveRDS(list(source_url = "https://osf.io/download/rvcfs/", project = "https://osf.io/t8ax6/",
  paper = "https://doi.org/10.1093/applin/amae087", input_md5 = unname(tools::md5sum(input)),
  script_md5 = tools::md5sum("validation/prepare_l2_vocabulary_bank.R"), seed = 96710001L,
  packages = vapply(c("mirt", "readxl", "glmmplanr"), function(p) as.character(packageVersion(p)), "")),
  file.path(out, "provenance.rds"))
capture.output(sessionInfo(), file = file.path(out, "environment.txt"))
print(do.call(rbind, summaries), row.names = FALSE)
