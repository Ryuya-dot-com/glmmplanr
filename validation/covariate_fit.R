# Shared fitter for the external trial-covariate benchmarks.
control <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 30000))

fit_slopes <- function(data, formula, link, terms) {
  messages <- character()
  result <- list(ok = FALSE, singular = NA, diagnostic = "", b = numeric(), V = NULL)
  tryCatch({
    fit <- withCallingHandlers(suppressMessages(glmer(formula, data, family = binomial(link),
      nAGQ = 1, control = control)), warning = function(w) {
        messages <<- c(messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    opt <- fit@optinfo$conv
    conv <- opt$lme4$messages
    conv <- conv[!grepl("boundary (singular) fit", conv, fixed = TRUE)]
    result$singular <- isSingular(fit)
    result$diagnostic <- paste(unique(c(conv, messages)), collapse = " | ")
    if (any(opt$opt != 0) || length(conv)) return(result)
    result$b <- fixef(fit, add.dropped = TRUE)[terms]
    result$V <- as.matrix(vcov(fit))[terms, terms, drop = FALSE]
    result$ok <- all(is.finite(result$b)) && all(is.finite(result$V)) && all(diag(result$V) > 0)
    if (!result$ok) result$diagnostic <- "Nonfinite slope or invalid covariance"
    result
  }, error = function(e) {
    result$diagnostic <- paste(c(messages, conditionMessage(e)), collapse = " | ")
    result
  })
}
