library(glmmplanr)
args <- commandArgs(TRUE); out <- args[1]
p <- list.files(out, pattern='^case-.*-block-01.rds$', full.names=TRUE)
objects <- lapply(p, readRDS)
details <- lapply(objects, function(z) z$details[[1]])
d <- read.csv(file.path(out,'diagnostics.csv'))
d <- d[d$replicate == 1L, ]
print(d[,c('case','native_logsd_score','profile_score','objective_improvement','fine_se_difference_pp','alternate_objective_difference','elapsed','evaluations')],row.names=FALSE)
dense <- do.call(rbind,lapply(details, function(z) do.call(rbind,z$dense)))
print(apply(abs(dense),2,max))
stopifnot(nrow(dense)==32,all(is.finite(dense)),max(abs(dense[,'ll_difference'])) < 1e-6,
  max(abs(dense[,'mode_score'])) < 2e-5,
  all(d$native_valid),
  all(abs(d$fine_se_difference_pp[d$profile_valid]) < .05),
  all(abs(d$alternate_objective_difference[d$profile_valid]) < 1e-4))
if (readRDS(file.path(out,'manifest.rds'))$settings$mode == 'smoke')
  stopifnot(all(d$profile_valid))
for(j in seq_along(objects)) {
  o <- objects[[j]]; z <- o$details[[1]]
  cache <- glmmplanr:::prepare_cpp_schedule_cache(o$spec,o$schedule)
  f <- glmmplanr:::fit_via_cpp_cached(list(y=as.integer(z$y)),cache,o$spec,TRUE,o$spec$analysis_re,FALSE)
  row <- z$rows[z$rows$method=='native',]
  stopifnot(f$converged, abs(100*f$prob_est-row$estimate) < 1e-8,
    abs(100*f$prob_se-row$se)<1e-8)
}
cat('All 16 dense likelihood/mode and public-core consistency checks passed.\n')
cat('Curvature-step and alternate-start checks passed for',sum(d$profile_valid),
  'valid first-dataset profile fits; invalid cases remain excluded:',
  paste(d$case[!d$profile_valid],collapse=','),'\n')
