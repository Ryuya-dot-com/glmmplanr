script <- system.file("validation", "three_choice_confirmation.R", package = "glmmplanr", mustWork = TRUE)
status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script), "validation"))
if (status != 0L) stop("Three-choice confirmation failed.")
