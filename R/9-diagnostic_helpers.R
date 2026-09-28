# The vocabulary and the machinery the refit diagnostics share: the argument
# documentation retrospective(), jitter() and self_test() inherit, and the
# parallel dispatcher those three plus profile() and run_mse() refit through.

#' Shared refit-diagnostic arguments
#'
#' A common vocabulary, documented here once and inherited by the diagnostics
#' that refit the model: [retrospective()], [jitter()] and [self_test()]. Each
#' argument means the same thing wherever it appears. `phase`, `getsd` and
#' `timeout` are deliberately **not** here, their defaults and their
#' consequences for the diagnostic differ, so each function documents its own.
#' [profile()] refits too but takes only `cores` and a bare `getsd`.
#'
#' @section What `fit_control()` reaches:
#'
#' These diagnostics refit through an internal helper that forwards `phase` and
#' `getsd` and nothing else. Every other field either comes off the fitted
#' model's own `data_list` (`loopnum`, `comp_offset`, the bias-adjustment flags)
#' or takes `fit_control()`'s default. Setting one of those is an error rather
#' than a silent no-op, so a field you set is a field the refit used.
#'
#' Only the fields you set are applied, named in the `fit_control()` call, or
#' assigned to afterwards (`ctl$getsd <- TRUE`). One you never touch keeps the
#' diagnostic's own default, which is not always `fit_control()`'s,
#' [retrospective()] phases its peels where `fit_control()` does not, so
#' `fit_control(getsd = FALSE)` asks about standard errors and changes nothing
#' else. Setting a field to the value that is already `fit_control()`'s default
#' still counts as setting it.
#'
#' @param object an Rceattle model fit using [fit_mod()].
#' @param Rceattle deprecated name for `object`, still accepted so existing
#'   scripts keep working. Supplying both is an error.
#' @param cores Number of cores for the parallel refits. Default `NULL` picks
#'   `parallel::detectCores() - 6`, capped at 2 under `R CMD check` (which sets
#'   `_R_CHECK_LIMIT_CORES_`). Set to 1 to force sequential execution.
#' @param fit_control optional [fit_control()] bundle for the refits. Only
#'   `phase` and `getsd` are read; see **What `fit_control()` reaches**.
#' @name rceattle-refit-args
NULL


#' Parallel `lapply` over a cluster, using FORK where the platform allows.
#'
#' @description
#' On non-Windows platforms a FORK cluster inherits the parent process's memory
#' (the already-loaded `Rceattle` namespace and every local, including the large
#' fitted OM/EM objects) via copy-on-write, so it needs neither a per-worker
#' `library(Rceattle)` nor a `clusterExport()` of those objects. The PSOCK path
#' pays both on every worker: a package load and serialization/transfer of the
#' exported bindings. FORK therefore removes the dominant cluster-startup cost
#' for the retrospective / jitter / MSE dispatchers. PSOCK is the cross-platform
#' fallback (Windows), reproducing the previous behavior exactly.
#'
#' Each worker fits whole models, so several large ones at once can exhaust
#' memory and the operating system kills one. All the parent sees is
#' `Error in unserialize(node$con) : error reading from connection`, and every
#' item finished up to that point is lost with it. Fall back to running the items
#' sequentially instead, as `osa_residuals()` does when its parallel
#' one-step-ahead loop fails.
#'
#' Two things follow. The retry starts from the beginning, results from a
#' cluster that has lost a worker cannot be recovered, so a caller with side
#' effects repeats them; `run_mse()` re-writes the `.rds` of any simulation that
#' had already finished, which is safe only because it seeds each simulation
#' separately and so reproduces it exactly. And an error raised by `fun` itself
#' propagates out of the retry unchanged, so a real bug still surfaces as itself
#' rather than as a cluster failure.
#'
#' @param items Vector iterated over; each element is passed to `fun`.
#' @param fun Worker closure.
#' @param n_workers Number of cluster workers.
#' @param export_env For the PSOCK fallback only, the environment whose bindings
#'   are exported to the workers (the caller's `environment()`). Ignored for
#'   FORK, where the workers inherit it directly.
#' @return A list of `fun` applied to each element of `items`.
#' @noRd
.parallel_lapply <- function(items, fun, n_workers, export_env) {
  fork <- .Platform$OS.type != "windows"

  run_clustered <- function() {
    cl <- if (fork) {
      parallel::makeCluster(n_workers, type = "FORK")
    } else {
      parallel::makeCluster(n_workers)
    }
    on.exit(parallel::stopCluster(cl), add = TRUE)
    if (!fork) {
      parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(Rceattle)))
      parallel::clusterExport(cl, varlist = ls(envir = export_env), envir = export_env)
    }
    parallel::parLapply(cl, items, fun)
  }

  tryCatch(run_clustered(), error = function(e) {
    message("Parallel run failed (", conditionMessage(e), "); computing the ",
            length(items), " tasks sequentially instead. A worker most often ",
            "dies because it ran out of memory, since each one fits a full ",
            "model -- try a lower cores, or cores = 1 to skip the parallel ",
            "attempt.")
    lapply(items, fun)
  })
}
