# Summaries and matched Monte Carlo contrasts for rasch_fit_benchmark.R.
args <- commandArgs(TRUE)
src <- if (length(args)) args[1L] else "validation"
dest <- if (length(args) > 1L) args[2L] else src
reference <- if (length(args) > 2L) args[3L] else NULL
engines <- if (is.null(reference)) {
  c("glmmplanr", "TAM", "mfrmr")
} else {
  c("glmmplanr", "TAM", "mfrmr 0.2.4", "mfrmr 0.3.2")
}
read_result <- function(task) {
  cases <- read.csv(file.path(src, paste0(task, "_scenarios.csv")))
  file <- paste0(task, "_", task, "_fits.csv.gz")
  fits <- read.csv(file.path(src, file))
  if (!is.null(reference)) {
    original_cases <- read.csv(file.path(reference, paste0(task, "_scenarios.csv")))
    original <- read.csv(file.path(reference, file))
    stopifnot(identical(cases, original_cases), all(fits$engine == "mfrmr"),
              grepl("mfrmr 0.3.2", paste(readLines(file.path(src, paste0(task, "_environment.txt"))),
                                        collapse = "\n"), fixed = TRUE))
    data_check <- if (task == "calibration") "correct" else c("group_n", "focal_n")
    columns <- c("case", "replicate", data_check)
    saved <- original[original$engine == "mfrmr", columns]
    current <- fits[columns]
    checked <- merge(saved, current, by = c("case", "replicate"))
    stopifnot(nrow(checked) == nrow(saved), nrow(checked) == nrow(current))
    for (column in data_check) {
      stopifnot(identical(checked[[paste0(column, ".x")]], checked[[paste0(column, ".y")]]))
    }
    original$engine[original$engine == "mfrmr"] <- "mfrmr 0.2.4"
    fits$engine <- "mfrmr 0.3.2"
    fits <- rbind(original, fits)
  }
  key <- fits[c("case", "replicate", "engine")]
  want <- expand.grid(case = cases$case, replicate = seq_len(unique(cases$n_rep)), engine = engines,
                      stringsAsFactors = FALSE)
  colnames(want) <- names(key)
  stopifnot(!anyDuplicated(key), nrow(want) == nrow(key),
            nrow(merge(want, key)) == nrow(key))
  rownames(cases) <- cases$case
  list(cases = cases, fits = fits)
}

interval <- function(value, binary = FALSE, scheduled = length(value)) {
  good <- value[is.finite(value)]
  n <- length(good)
  if (binary) {
    successes <- sum(good)
    bounds <- if (n) binom.test(successes, n)$conf.int else c(NA_real_, NA_real_)
    return(data.frame(scheduled, valid = n, successes, estimate = if (n) mean(good) else NA_real_,
                      mc_lower = bounds[1L], mc_upper = bounds[2L]))
  }
  center <- if (n) mean(good) else NA_real_
  half <- if (n > 1L) qt(.975, n - 1L) * sd(good) / sqrt(n) else NA_real_
  data.frame(scheduled, valid = n, successes = NA_integer_, estimate = center,
             mc_lower = center - half, mc_upper = center + half)
}

rate <- function(x, valid_only = FALSE) {
  if (!valid_only) x[is.na(x)] <- FALSE
  interval(x, binary = TRUE)
}

summ_row <- function(z, case, task) {
  if (task == "calibration") {
    z$bank_attainment <- z$bank_ok == 1L
    z$infit_screen <- z$focal_infit > 1.2
    z$outfit_screen <- z$focal_outfit > 1.2
    metrics <- c("correct", "sigma", "b_bias", "b_rmse", "b_max", "theta_rmse",
      "person_eap_rel", "person_eap_ratio", "item_rel", "item_separation", "mean_item_infit",
      "mean_item_outfit", "mean_person_infit", "mean_person_outfit", "focal_infit", "focal_outfit",
      "person_overfit_pct", "person_underfit_pct")
    proportions <- c("valid", "bank_attainment", "infit_screen", "outfit_screen")
    if (any(z$native_scheduled)) {
      zz <- z[z$native_scheduled, ]
      native_names <- grep("^native_(item_|person_|focal_)", names(zz), value = TRUE)
      native <- c(lapply(zz[native_names], interval),
                  list(native_infit_screen = rate(zz$native_focal_infit > 1.2),
                       native_outfit_screen = rate(zz$native_focal_outfit > 1.2),
                       native_fit_available = rate(is.finite(zz$native_focal_infit))))
    } else {
      native <- list()
    }
  } else {
    z$dif_bias <- z$estimate - case$dif
    z$dif_absolute_error <- abs(z$dif_bias)
    z$dif_squared_error <- z$dif_bias^2
    z$group_bias <- z$group - case$group_shift
    z$lr_approx_difference <- z$lr - z$lr_common
    metrics <- c("estimate", "dif_bias", "dif_absolute_error", "group_bias", "sigma",
                 "lr", "lr_common", "lr_approx_difference")
    proportions <- c("null_ok", "alt_ok", "valid", "test_valid", "rejected")
    native <- list(dif_rmse = transform(interval(z$dif_squared_error),
      estimate = sqrt(estimate), mc_lower = sqrt(pmax(0, mc_lower)), mc_upper = sqrt(pmax(0, mc_upper))))
  }
  valid_names <- setdiff(proportions, c("valid", "null_ok", "alt_ok", "test_valid"))
  valid_only <- setNames(lapply(z[valid_names], rate, valid_only = TRUE), paste0(valid_names, "_valid"))
  data <- c(lapply(z[metrics], interval), lapply(z[proportions], rate), valid_only, native)
  data$fit_cpu_median_s <- interval(rep(NA_real_, nrow(z)))
  data$fit_cpu_median_s$estimate <- median(z$fit_cpu_s)
  data$fit_cpu_median_s$valid <- sum(is.finite(z$fit_cpu_s))
  if (!is.null(reference)) {
    data$fit_cpu_mean_s <- interval(z$fit_cpu_s)
    data$fit_cpu_p95_s <- data$fit_cpu_median_s
    data$fit_cpu_p95_s$estimate <- unname(quantile(z$fit_cpu_s, .95))
    if (task == "calibration" && any(z$native_scheduled)) {
      zz <- z[z$native_scheduled, ]
      data$fit_and_native_cpu_mean_s <- interval(zz$fit_cpu_s + zz$native_cpu_s)
    }
  }
  rows <- do.call(rbind, data)
  data.frame(case, engine = z$engine[1L], metric = rownames(rows), rows, row.names = NULL)
}

comparison <- function(a, b, task) {
  stopifnot(identical(sort(a$replicate), sort(b$replicate)), !anyDuplicated(a$replicate))
  z <- merge(a, b, by = "replicate", suffixes = c("_1", "_2"))
  if (task == "calibration") {
    numeric_metrics <- c("b_rmse", "theta_rmse", "person_eap_rel", "item_rel", "item_separation",
                         "focal_infit", "focal_outfit", "ll_common")
    binary_metrics <- c("valid", "bank_ok")
  } else {
    numeric_metrics <- c("estimate", "group", "lr")
    binary_metrics <- c("test_valid", "rejected")
  }
  if (!is.null(reference)) {
    numeric_metrics <- c(numeric_metrics, "sigma", "fit_cpu_s")
    if (task == "calibration" && any(z$native_scheduled_1 & z$native_scheduled_2)) {
      numeric_metrics <- c(numeric_metrics,
        "native_focal_infit", "native_focal_outfit", "native_person_infit", "native_person_outfit",
        "native_person_rel", "native_person_separation", "native_item_rel", "native_item_separation",
        "native_cpu_s")
    }
  }
  do.call(rbind, lapply(c(numeric_metrics, binary_metrics), function(metric) {
    x <- z[[paste0(metric, "_1")]]
    y <- z[[paste0(metric, "_2")]]
    binary <- metric %in% binary_metrics
    if (binary) {
      x[is.na(x)] <- FALSE
      y[is.na(y)] <- FALSE
    }
    d <- x - y
    data.frame(metric, interval(d), first_only = if (binary) sum(d == 1) else NA_integer_,
               second_only = if (binary) sum(d == -1) else NA_integer_,
               max_abs = if (any(is.finite(d))) max(abs(d), na.rm = TRUE) else NA_real_)
  }))
}

for (task in c("calibration", "dif")) {
  result <- read_result(task)
  cases <- result$cases
  fits <- result$fits
  blocks <- split(fits, list(fits$case, fits$engine), drop = TRUE)
  rows <- lapply(blocks, function(z) summ_row(z, cases[as.character(z$case[1L]), ], task))
  summary <- do.call(rbind, rows)
  write.csv(summary, file.path(dest, paste0(task, "_summary.csv")), row.names = FALSE)
  by_case <- split(fits, fits$case)
  paired <- lapply(by_case, function(z) {
    do.call(rbind, lapply(combn(engines, 2L, simplify = FALSE), function(pair) {
      data.frame(comparison = "engine", case1 = z$case[1L], case2 = z$case[1L],
                 engine1 = pair[1L], engine2 = pair[2L],
                 comparison(z[z$engine == pair[1L], ], z[z$engine == pair[2L], ], task))
    }))
  })
  # The equally exposed contrast holds the item bank, groups, responses per
  # item, and generated profile constant within each replication.
  at <- cases[cases$n == 480L & cases$items_per_person == 8L, ]
  same_exposure <- cases[cases$responses_per_item == 240L & cases$items_per_person == 16L, ]
  if (task == "dif") {
    # The 480-person full comparison holds total persons and subgroup means fixed.
    same_exposure <- cases[cases$n == 480L & cases$items_per_person == 16L, ]
  }
  contrasts <- merge(at[c("case", "profile")], same_exposure[c("case", "profile")], by = "profile")
  by_layout <- lapply(seq_len(nrow(contrasts)), function(k) {
    a <- contrasts$case.x[k]
    b <- contrasts$case.y[k]
    do.call(rbind, lapply(engines, function(engine) {
      data.frame(comparison = "booklets_minus_complete", case1 = a, case2 = b,
                 engine1 = engine, engine2 = engine,
                 comparison(fits[fits$case == a & fits$engine == engine, ],
                            fits[fits$case == b & fits$engine == engine, ], task))
    }))
  })
  paired <- do.call(rbind, c(paired, by_layout))
  write.csv(paired, file.path(dest, paste0(task, "_paired.csv")), row.names = FALSE)
}
