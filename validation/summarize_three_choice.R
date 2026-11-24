# Usage: Rscript --vanilla validation/summarize_three_choice.R [directory containing the two fit CSVs]
args <- commandArgs(TRUE)
stopifnot(length(args) <= 1L)
directory <- if (length(args)) args[1] else "validation"
path <- function(table) file.path(directory, paste0("three_choice_precision_", table, ".csv"))
cpp <- read.csv(path("native"))
stan <- read.csv(path("stan"))
sched <- expand.grid(k = 1:4, trial = 1:192)
stopifnot(identical(cpp$seed, with(sched, 7310000L + 1000L * k + trial)),
          identical(cpp$k, sched$k), identical(cpp$trial, sched$trial),
          identical(stan$seed, cpp$seed[cpp$trial <= 40]), !anyNA(cpp$cpp_valid),
          !anyNA(stan$stan_valid), identical(stan$retried, !stan$initial_valid))
pairs <- merge(cpp, stan, by = "seed", sort = FALSE)
z <- qnorm(0.975)
cp <- function(hits, n) {
  if (!n) return(c(NA_real_, NA_real_))
  as.numeric(binom.test(hits, n)$conf.int)
}
judge <- function(ci, threshold) {
  if (anyNA(ci)) return("indeterminate")
  if (ci[1] >= threshold) return("pass")
  if (ci[2] < threshold) return("below_floor")
  "indeterminate"
}
summarize_fits <- function(df, engine, series) {
  valid <- df[[paste0(engine, "_valid")]]
  se <- df[[paste0(engine, "_se")]][valid]
  err <- df[[paste0(engine, "_est")]][valid] - df$truth[valid]
  stopifnot(all(is.finite(err)), all(is.finite(se) & se > 0))
  n <- nrow(df)
  nv <- sum(valid)
  covered <- sum(abs(err) <= z * se)
  width <- 2 * z * se
  narrow <- sum(width <= 1)
  coverage_ci <- cp(covered, nv)
  width_ci <- cp(narrow, n)
  coverage_status <- judge(coverage_ci, 0.90)
  width_status <- judge(width_ci, 0.80)
  joint <- if (nv / n < 0.9) "insufficient_valid" else if ("below_floor" %in% c(coverage_status, width_status)) {
    "below_floor"
  } else if (all(c(coverage_status, width_status) == "pass")) "pass" else "indeterminate"
  data.frame(persons = df$persons[1], items = df$items[1], series = series,
    attempted = n, valid = nv, covered = covered, coverage = covered / nv,
    coverage_lo = coverage_ci[1], coverage_hi = coverage_ci[2], coverage_status = coverage_status,
    narrow = narrow, assurance = narrow / n, width_lo = width_ci[1], width_hi = width_ci[2],
    width_status = width_status, joint_status = joint, bias = mean(err), rmse = sqrt(mean(err^2)),
    empirical_sd = sd(err), mean_uncertainty = mean(se), mean_width = mean(width))
}
summaries <- paired <- vector("list", 4L)
for (k in seq_along(summaries)) {
  all <- cpp[cpp$k == k, ]
  matched <- pairs[pairs$k == k, ]
  summaries[[k]] <- rbind(summarize_fits(all, "cpp", "native_all"),
                          summarize_fits(matched, "cpp", "native_matched"),
                          summarize_fits(matched, "stan", "stan_matched"))
  m <- matched[matched$cpp_valid & matched$stan_valid, ]
  log_ratio <- log(m$stan_se / m$cpp_se)
  avg <- mean(log_ratio)
  half <- if (nrow(m) > 1L) qt(.975, nrow(m) - 1L) * sd(log_ratio) / sqrt(nrow(m)) else NA_real_
  paired[[k]] <- data.frame(persons = all$persons[1], items = all$items[1], pairs = nrow(m),
    ratio = exp(avg), ratio_lo = exp(avg - half), ratio_hi = exp(avg + half),
    stan_wider = sum(log_ratio > 0), native_subj_sd = mean(m$cpp_subj_sd),
    stan_subj_sd = mean(m$stan_subj_sd), native_item_sd = mean(m$cpp_item_sd),
    stan_item_sd = mean(m$stan_item_sd), stan_retries = sum(matched$retried),
    stan_invalid = sum(!matched$stan_valid),
    stan_final_warned = sum(!is.na(matched$stan_warnings) & nzchar(matched$stan_warnings)))
}
out <- do.call(rbind, summaries)
p <- do.call(rbind, paired)
write.csv(out, path("summary"), row.names = FALSE)
write.csv(p, path("paired"), row.names = FALSE)
print(out[c("persons", "items", "series", "valid", "covered", "narrow", "width_status", "joint_status")], row.names = FALSE)
print(p, row.names = FALSE)
