files <- list.files('tests/testthat', '^test-.*[.]R$')
groups <- sub('[.]R$', '', sub('^test-', '', files))
remaining <- setdiff(groups, c('laplace-status', 'laplace-covariance'))
stopifnot(length(remaining) == length(files)-2L)
cat('Installed package:', find.package('glmmplanr'), '\n')
print(sessionInfo())
r <- testthat::test_dir('tests/testthat', filter = paste(remaining, collapse = '|'),
  package = 'glmmplanr', load_package = 'installed', reporter = 'summary')
saveRDS(r, '/private/tmp/glmmplanr-curvature-repair/regression-results.rds')
