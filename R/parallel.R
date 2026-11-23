# Worker-Invariant Parallel Task Execution (mirai / multicore / sequential)

split_indices <- function(n, k) {
  if (k <= 1L) return(list(seq_len(n)))
  grp <- cut(seq_len(n), breaks = k, labels = FALSE)
  unname(split(seq_len(n), grp))
}

start_mirai_pool <- function(workers) {
  tryCatch({
    # Pool setup may draw from the caller even when all tasks have explicit seeds.
    with_rng_stream(1L, mirai::daemons(workers))
    is_dev <- isNamespaceLoaded("pkgload") &&
      isTRUE(tryCatch(pkgload::is_dev_package("glmmplanr"), error = function(e) FALSE))
    if (is_dev) {
      pkg_root <- getNamespaceInfo(asNamespace("glmmplanr"), "path")
      init_res <- mirai::everywhere(
        pkgload::load_all(p, compile = FALSE, quiet = TRUE),
        p = pkg_root
      )[]
      if (any(vapply(init_res, mirai::is_error_value, logical(1L)))) {
        try(mirai::daemons(0L), silent = TRUE)
        return(FALSE)
      }
    }
    TRUE
  }, error = function(e) {
    try(mirai::daemons(0L), silent = TRUE)
    FALSE
  })
}

run_parallel_tasks <- function(
  tasks,
  worker_fn,
  workers,
  parallel_backend,
  stage_label = NULL,
  verbose = FALSE
) {
  n_tasks <- length(tasks)
  if (n_tasks == 0L) return(list())

  task_weights <- vapply(tasks, function(t) {
    if (is.list(t) && !is.null(t$rep_idx)) length(t$rep_idx) else 1L
  }, integer(1L))
  total_reps <- sum(task_weights)
  lbl <- if (!is.null(stage_label) && nzchar(stage_label)) stage_label else "Simulating"

  show_pb <- isTRUE(verbose) && total_reps > 1L && n_tasks > 1L
  pb_id <- NULL
  if (show_pb) {
    pb_id <- cli::cli_progress_bar(
      name = lbl,
      total = total_reps,
      format = paste0(
        "{cli::pb_spin} {cli::pb_name} [{cli::pb_bar}] ",
        "{cli::pb_current}/{cli::pb_total} sims ({cli::pb_percent}) | ",
        "Elapsed: {cli::pb_elapsed} | ETA: {cli::pb_eta}"
      ),
      clear = FALSE
    )
    on.exit(cli::cli_progress_done(id = pb_id), add = TRUE)
  }

  run_seq <- function(idx_vec) {
    out <- vector("list", length(idx_vec))
    for (k in seq_along(idx_vec)) {
      ii <- idx_vec[k]
      out[[k]] <- worker_fn(tasks[[ii]])
      if (!is.null(pb_id)) cli::cli_progress_update(inc = task_weights[ii], id = pb_id)
    }
    out
  }

  if (workers <= 1L || parallel_backend == "sequential" || n_tasks <= 1L) {
    return(run_seq(seq_len(n_tasks)))
  }

  chk <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  if (nzchar(chk) && chk != "false") {
    workers <- min(workers, 2L)
  }

  n_waves <- if (show_pb) min(n_tasks, max(1L, ceiling(n_tasks / max(workers * 2L, 2L)))) else 1L
  waves <- split_indices(n_tasks, n_waves)

  is_macos <- identical(Sys.info()[["sysname"]], "Darwin")
  prefer_mirai <- parallel_backend == "mirai" ||
    (parallel_backend == "auto" &&
       (isTRUE(mirai::daemons_set()) || .Platform$OS.type != "unix" || is_macos))

  if (prefer_mirai) {
    res <- tryCatch({
      started_daemons <- FALSE
      if (!isTRUE(mirai::daemons_set())) {
        if (!start_mirai_pool(workers)) stop("mirai pool init failed", call. = FALSE)
        started_daemons <- TRUE
      }
      on.exit(if (started_daemons) try(mirai::daemons(0L), silent = TRUE), add = TRUE)
      full_out <- vector("list", n_tasks)
      # A pool supplied by the caller can be larger than the requested limit.
      if (!started_daemons && nzchar(chk) && chk != "false") {
        waves <- split(seq_len(n_tasks), ceiling(seq_len(n_tasks) / workers))
      }
      for (wv in waves) {
        sub_out <- mirai::mirai_map(
          tasks[wv], function(tsk, fn) fn(tsk), .args = list(fn = worker_fn)
        )[]
        if (any(vapply(sub_out, mirai::is_error_value, logical(1L)))) {
          stop("mirai worker error", call. = FALSE)
        }
        full_out[wv] <- sub_out
        if (!is.null(pb_id)) cli::cli_progress_update(inc = sum(task_weights[wv]), id = pb_id)
      }
      full_out
    }, error = function(e) NULL)
    if (!is.null(res)) return(res)
    # Platform or backend selection already favored non-forked workers.
    return(run_seq(seq_len(n_tasks)))
  }

  if (.Platform$OS.type == "unix") {
    old_omp <- Sys.getenv("OMP_NUM_THREADS", unset = NA_character_)
    old_vec <- Sys.getenv("VECLIB_MAXIMUM_THREADS", unset = NA_character_)
    Sys.setenv(OMP_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
    on.exit({
      if (is.na(old_omp)) Sys.unsetenv("OMP_NUM_THREADS") else Sys.setenv(OMP_NUM_THREADS = old_omp)
      if (is.na(old_vec)) Sys.unsetenv("VECLIB_MAXIMUM_THREADS") else Sys.setenv(VECLIB_MAXIMUM_THREADS = old_vec)
    }, add = TRUE)

    full_out <- vector("list", n_tasks)
    ok_mc <- TRUE
    for (wv in waves) {
      use_presched <- length(wv) > (workers * 4L)
      sub_out <- parallel::mclapply(
        tasks[wv], worker_fn, mc.cores = workers, mc.preschedule = use_presched,
        mc.set.seed = FALSE # Each replication installs its assigned substream.
      )
      bad <- vapply(sub_out, function(x) inherits(x, "try-error") || is.null(x), logical(1L))
      if (any(bad)) {
        ok_mc <- FALSE
        break
      }
      full_out[wv] <- sub_out
      if (!is.null(pb_id)) cli::cli_progress_update(inc = sum(task_weights[wv]), id = pb_id)
    }
    if (ok_mc) return(full_out)
  }

  run_seq(seq_len(n_tasks))
}
