library(testthat)
x <- parse('tests/testthat/test-core-workflows.R', keep.source=TRUE)
selected <- Filter(function(e) any(grepl('skip_on_cran', deparse(e), fixed=TRUE)), as.list(x))
stopifnot(length(selected) == 2L)
dir.create('/private/tmp/glmmplanr-curvature-repair/remaining-tests', showWarnings=FALSE)
writeLines(unlist(lapply(selected, deparse, width.cutoff=120L)),
  '/private/tmp/glmmplanr-curvature-repair/remaining-tests/test-core-remaining.R')
Sys.setenv(NOT_CRAN='true')
r <- test_dir('/private/tmp/glmmplanr-curvature-repair/remaining-tests',
  package='glmmplanr', load_package='installed', reporter='summary')
saveRDS(r, '/private/tmp/glmmplanr-curvature-repair/remaining-results.rds')
