#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Renamed-argument support
#
# The diagnostics used six different names for the same thing -- the fitted
# model they operate on (`Rceattle`, `fit`, `object`, ...). They now all take
# `object`. The old spelling is still accepted, so every existing assessment
# script keeps running unchanged.
#
# The old name is declared LAST in each signature, so a positional call such as
# `retrospective(ss_run, peels = 10)` still binds the model to `object`. The old
# name can therefore only be supplied by name, which is how the scripts in
# ../Rceattle-models and ../GOA-ATF-ESP already write it. An abbreviation of the
# old name (`osa_residuals(fi = m)`) is not supported -- it never was, for the
# functions that take `...`.
#
# BOTH formals must carry a default. `retrospective()`, `jitter()` and
# `self_test()` export the caller's frame to PSOCK workers with
# `parallel::clusterExport(cl, ls(envir = export_env), export_env)`, and `ls()`
# on a function frame lists formals that were never supplied. A formal with no
# default therefore makes `clusterExport()` throw "argument is missing", which
# `.parallel_lapply()` catches -- so every peel, jitter and self-test would
# silently drop to serial on Windows while reporting a memory problem. FORK
# (macOS/Linux) skips clusterExport entirely and would never show it.
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

#' Should a renamed argument warn?
#'
#' @description
#' One switch for the whole package. The old argument names are accepted
#' silently in 5.11.0-5.12.0 so that a release which renames nothing a user can see
#' does not start printing warnings inside long assessment and MSE loops. The
#' warning turns on in 5.13.0 by flipping this to `TRUE`, one edit, not ten.
#'
#' @return `TRUE` if deprecated argument names should warn, otherwise `FALSE`.
#' @keywords internal
#' @noRd
.rce_deprecate_warn <- function() {
  isTRUE(getOption("Rceattle.warn_deprecated_args", default = FALSE))
}

#' Resolve a renamed argument
#'
#' @description
#' Called only when the caller supplied the deprecated name. Returns the value
#' to use for the new argument.
#'
#' Supplying both names is an error rather than a silent preference: the two
#' would have to be the same fitted model for the call to mean anything, and
#' guessing which one the user meant is exactly the kind of assumption that
#' produces a plausible wrong answer.
#'
#' @param old value bound to the deprecated argument name.
#' @param new_supplied `TRUE` if the caller also gave the current name. Pass
#'   `!missing(<new>)`, not `!is.null(<new>)`, an explicit
#'   `object = NULL` alongside the old name is still two spellings of one
#'   argument, and should be rejected rather than quietly resolved.
#' @param old_name,new_name the two spellings, for the message.
#' @param fn calling function, for the message.
#'
#' @return `old`, the value supplied under the deprecated name.
#' @keywords internal
#' @noRd
.rce_deprecated_arg <- function(old, new_supplied, old_name, new_name, fn) {
  if (new_supplied) {
    stop(sprintf(
      "%s(): give the fitted model once. `%s` is the current name for it and `%s` is the deprecated one, but both were supplied.",
      fn, new_name, old_name), call. = FALSE)
  }
  if (.rce_deprecate_warn()) {
    warning(sprintf(
      "%s(): argument `%s` has been renamed to `%s`. The old name still works and is scheduled for removal in Rceattle 6.0.0.",
      fn, old_name, new_name), call. = FALSE)
  }
  old
}
