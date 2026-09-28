# =============================================================================
# model_config() -- a defaulted, validated model-configuration slot
# =============================================================================
#
# fit_mod() takes the model *structure* (predation mode, initialization, HCR,
# and the build_*() process specifications) as call-time arguments. That makes a
# configuration live only in the fit_mod() call, not on the data object, so two
# runs of the same model are only comparable by eye-diffing their fit_mod()
# calls. model_config() lets the configuration travel WITH the data list as a
# `$model_config` slot, so a data object records how it is meant to be fit.
#
# fit_mod()'s signature is unchanged. When a data list carries a model_config,
# fit_mod() uses each stored field ONLY for arguments the caller did not supply
# (detected with missing()); an explicitly-passed argument always wins. With no
# model_config present, fit_mod() reads every setting from its own arguments.

#' Build a model-configuration slot for a data list
#'
#' Bundles the model-structure arguments of [fit_mod()], predation mode,
#' initialization, harvest control rule, and the `build_*()` process
#' specifications, into a single validated object that can be stored on a data
#' list (`data_list$model_config`, e.g. via
#' `build_data(..., model_config = model_config(...))`). A configuration then
#' travels with the data instead of living only in the [fit_mod()] call.
#'
#' The defaults are exactly [fit_mod()]'s own argument defaults, so a data list
#' holding `model_config()` fits identically to one with no slot at all. When a
#' data list has a `model_config`, [fit_mod()] reads each field only for
#' arguments the caller did **not** pass explicitly, an argument passed to
#' [fit_mod()] always wins, even when passed at its default. Omit an argument to
#' let the stored configuration take effect for that field.
#'
#' @param msmMode Predation-mortality mode (see [fit_mod()]). Default `0`.
#' @param initMode Population initialization (see [fit_mod()]). Default
#'   `"NonEquilibrium"`.
#' @param avgnMode Average-N mode. Default `0`.
#' @param suitMode Suitability mode. Default `0`.
#' @param niter Number of predation iterations. Default `3`.
#' @param HCR A harvest-control-rule specification from [build_hcr()].
#' @param recFun A stock-recruit specification from [build_srr()].
#' @param M1Fun A natural-mortality specification from [build_M1()].
#' @param growthFun A growth specification from [build_growth()].
#' @param qFun A catchability specification from [build_catchability()].
#' @param selFun A selectivity specification from [build_selectivity()].
#' @param compFun A composition specification from [build_composition()].
#'
#' @return An object of class `"Rceattle_model_config"` (a named list) to be
#'   stored as `data_list$model_config`.
#'
#' @section Persistence:
#' The configuration is code-side model structure, not one of the workbook data
#' sheets, so it is **not** written by [write_data()] and does not survive an
#' xlsx round-trip, `build_data(base = x, model_config = cfg)` piped through
#' [write_data()] then [read_data()] returns without the slot. Re-attach it in
#' code, store it alongside the data, or persist it as a documented, git-diffable
#' YAML with [save_config()] / [load_config()] and apply it with
#' `fit_mod(data_list, config = load_config("run.yaml"))`.
#'
#' @seealso [fit_mod()], [build_data()], [save_config()], [load_config()].
#' @examples
#' cfg <- model_config(msmMode = 1, initMode = "NonEquilibrium")
#' dat <- build_data(base = BS2017MS, model_config = cfg)
#' # fit_mod(dat) would then fit as multispecies without passing msmMode.
#' @export
model_config <- function(msmMode = 0,
                         initMode = "NonEquilibrium",
                         avgnMode = 0,
                         suitMode = 0,
                         niter = 3,
                         HCR = build_hcr(),
                         recFun = build_srr(),
                         M1Fun = build_M1(),
                         growthFun = build_growth(),
                         qFun = build_catchability(),
                         selFun = build_selectivity(),
                         compFun = build_composition()) {

  # Light validation of the scalar switches; the build_*() objects validate
  # themselves at construction.
  # model_config() takes settings, not data -- but almost every other entry point
  # in the package takes a data_list first, so name that mistake explicitly.
  if (is.list(msmMode)) {
    stop("`msmMode` must be a single value; a list was supplied. model_config() ",
         "takes model settings, not a data_list -- pass the data separately, as ",
         "build_data(base = <data_list>, model_config = model_config(...)).",
         call. = FALSE)
  }
  if (length(msmMode) != 1 || !is.numeric(msmMode) && !is.character(msmMode)) {
    stop("`msmMode` must be a single value.", call. = FALSE)
  }
  if (length(niter) != 1 || !is.numeric(niter) || niter < 1) {
    stop("`niter` must be a single positive integer.", call. = FALSE)
  }
  # Catch a mistyped switch here, at the point it was written, rather than
  # several functions later inside fit_mod(). Values are checked against the
  # switch maps but not converted -- the config records what the caller wrote.
  # suitMode is legitimately per-species, so these accept vectors.
  .check_switch(msmMode,  msmMode_map,  "msmMode")
  .check_switch(initMode, initMode_map, "initMode")
  .check_switch(suitMode, suitMode_map, "suitMode")
  # avgnMode has no string aliases -- it is inert in the C++ (only 0 is
  # implemented), so there is nothing to name.
  if (!is.null(avgnMode) && length(avgnMode) &&
      !all(is.na(avgnMode) | avgnMode %in% 0:2)) {
    stop("`avgnMode` must be 0, 1, or 2 (only 0 is implemented).", call. = FALSE)
  }

  cfg <- list(
    msmMode   = msmMode,
    initMode  = initMode,
    avgnMode  = avgnMode,
    suitMode  = suitMode,
    niter     = niter,
    HCR       = HCR,
    recFun    = recFun,
    M1Fun     = M1Fun,
    growthFun = growthFun,
    qFun      = qFun,
    selFun    = selFun,
    compFun   = compFun
  )
  # Which fields the caller set, defaults included: fit_mod(config =) overlays
  # only these onto a data object's own model_config, and save_config() writes
  # them even at their default so a saved run reproduces.
  attr(cfg, "set") <- intersect(.RCE_MODEL_CONFIG_FIELDS, names(match.call())[-1])
  class(cfg) <- c("Rceattle_model_config", "list")
  cfg
}

# The model_config fields fit_mod() resolves, in fit_mod() argument order. Used
# by fit_mod() to overlay the slot onto unsupplied arguments, and available to
# the spec-tree printer.
.RCE_MODEL_CONFIG_FIELDS <- c(
  "msmMode", "initMode", "avgnMode", "suitMode", "niter",
  "HCR", "recFun", "M1Fun", "growthFun", "qFun", "selFun", "compFun"
)

#' @export
print.Rceattle_model_config <- function(x, ...) {
  cat("<Rceattle model_config>\n")
  cat("  msmMode  :", x$msmMode, "\n")
  cat("  initMode :", x$initMode, "\n")
  cat("  avgnMode :", x$avgnMode, "\n")
  cat("  suitMode :", paste(x$suitMode, collapse = ", "), "\n")
  cat("  niter    :", x$niter, "\n")
  cat("  HCR      :", if (!is.null(x$HCR$HCR)) x$HCR$HCR else "(build_hcr)", "\n")
  cat("  recFun   : build_srr(srr_fun =", x$recFun$srr_fun, ")\n")
  cat("  M1Fun    : build_M1(M1_model =",
      paste(x$M1Fun$M1_model, collapse = ", "), ")\n")
  cat("  growthFun: build_growth(fun =",
      if (!is.null(x$growthFun$fun)) x$growthFun$fun else x$growthFun$growth_model,
      ")\n")
  # Only the fields recorded here overlay a data object's own model_config in
  # fit_mod(config = ); the rest keep the data's. Print them so the overlay can
  # be checked before fitting rather than from the per-field warnings it raises.
  set <- attr(x, "set")
  if (!is.null(set)) {
    cat("  imposes :", if (length(set)) paste(set, collapse = ", ") else "(nothing)", "\n")
  }
  invisible(x)
}
