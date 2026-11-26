# Reproduce the recorded mlpwr failure without generating or fitting new data.
RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
x <- readRDS("validation/item_pool_search/search-053.rds")
obs <- x$mlpwr_observations; grid <- x$grid; count <- 0L
simfun <- function(i, j) {
  count <<- count + 1L
  stopifnot(count <= nrow(obs), obs$pos[count] == match(paste(i, j), paste(grid$i, grid$j)))
  obs$success[count]
}
costfun <- function(i, j) grid$cost[match(paste(i, j), paste(grid$i, grid$j))]
set.seed(9470000L + 1000L * 53L + 201L)
failure <- tryCatch(withCallingHandlers(mlpwr::find.design(simfun,
  list(i = c(1, 3), j = c(1, 3)), power = .8, costfun = costfun,
  evaluations = 1080L, setsize = 60L, n.startsets = 3L, silent = TRUE),
  error = function(e) {
    cat(conditionMessage(e), "\n")
    for (f in sys.frames()) if (exists("xvalues", f, inherits = FALSE)) {
      cat("Matching indices in mlpwr::addval:\n"); print(f$ind)
      cat("Requested integer points:\n"); print(f$points)
      cat("Stored design points:\n"); print(lapply(f$dat, function(z) z$x))
    }
  }), error = identity)
stopifnot(inherits(failure, "error"), count == 360L,
  identical(conditionMessage(failure), "$ operator is invalid for atomic vectors"))
