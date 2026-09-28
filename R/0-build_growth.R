#' Allowed growth functions for `fun` in [build_growth()]
#'
#' Currently:
#'   * `"empirical"`: use the empirical weight-at-age input directly,
#'     no parameters estimated.
#'   * `"vonBertalanffy"`: sex-specific von Bertalanffy (parameters
#'     `K`, `L1`, `Linf`).
#'   * `"Richards"`: sex-specific Richards (adds `m`).
#'
#' @keywords internal
GROWTH_FUNS <- c("empirical", "vonBertalanffy", "Richards")


#' Internal: integer codes consumed by the TMB template
#' @keywords internal
#' @noRd
.GROWTH_FUN_TO_INT <- c(empirical = 0L, vonBertalanffy = 1L, Richards = 2L)


#' Internal: plus-group SD-at-age treatment codes consumed by the TMB template
#'
#' `"WHAM"` (1) pins the oldest age class's SD-at-age to the upper anchor
#' `exp(sd_Linf)` (the WHAM SDAA convention). `"SS3"` (2) instead interpolates
#' it by length like any interior age. Only affects estimated growth
#' (`growth_model > 0`); the default is WHAM.
#' @keywords internal
#' @noRd
.GROWTH_SD_STYLE <- c(WHAM = 1L, SS3 = 2L)


#' Allowed growth-parameter names for `linkages` in [build_growth()]
#'
#' Natural-scale names of the underlying growth-function parameters.
#' Von Bertalanffy uses `K`, `L1`, `Linf`; Richards adds `m`.
#' `sd_L1` / `sd_Linf` are the standard deviations of length-at-age
#' anchored at `L1` and `Linf` (the SD-at-age interpolation endpoints).
#' Only intercept-only specs (`~ 1`) are honored
#' on the SD endpoints, they thread through `init` / `bounds` /
#' `priors` onto the growth SD-at-age but do not vary by year.
#' The empirical weight-at-age model admits no linkages.
#'
#' @keywords internal
GROWTH_LINKAGE_PARAMS <- c("K", "L1", "Linf", "m", "sd_L1", "sd_Linf")


#' Mean-growth subset of [GROWTH_LINKAGE_PARAMS]
#'
#' Names that index into `log_growth_pars` / `growth_parameters` /
#' `growth_linkage_offset` along their last dim. The SD endpoints live
#' on a separate parameter (`growth_log_sd`) and are intentionally
#' excluded here.
#'
#' @keywords internal
.GROWTH_MEAN_PARAMS <- c("K", "L1", "Linf", "m")


#' Map mean-growth linkage param names to slices along the third dim of
#' `log_growth_pars` (`[nspp, nsex, n_growth_pars]`).
#' @keywords internal
#' @noRd
.GROWTH_PARAM_TO_INDEX <- c(K = 1L, L1 = 2L, Linf = 3L, m = 4L)


#' Map growth SD linkage param names to slices along the third dim of
#' `growth_log_sd` (`[nspp, nsex, 2]`).
#' @keywords internal
#' @noRd
.GROWTH_SD_PARAM_TO_INDEX <- c(sd_L1 = 1L, sd_Linf = 2L)


#' Specify the growth model for Rceattle
#'
#' @param fun Growth function. Either a string ([GROWTH_FUNS]:
#'   `"empirical"` (default), `"vonBertalanffy"`, `"Richards"`) or the
#'   equivalent integer code (`0`, `1`, `2`). The canonical string form
#'   is stored on the returned object.
#' @param growth_age_L1 Von Bertalanffy / Richards anchor age (the age at
#'   which mean length equals `L1`). Matches SS3's `Growth_Age_for_L1`
#'   control input. Scalar (recycled across species) or a length-`nspp`
#'   vector for per-species values. Default `NA` inherits
#'   `data_list$growth_age_L1` if supplied (e.g. from the SS3 converter),
#'   otherwise falls back to `max(0.5, minage[sp])` so `minage >= 1`
#'   models stay backwards-compatible and `minage = 0` models pick up an
#'   SS3-consistent half-year anchor.
#' @param sd_plus_group How the oldest age class's SD-at-age is treated (only
#'   affects estimated growth, `fun != "empirical"`). `"WHAM"` pins the
#'   plus-group SD to the upper anchor `exp(sd_Linf)` (the WHAM SDAA
#'   convention); `"SS3"` instead interpolates it by length like any interior
#'   age. Accepts a string or the integer code (`1`/`2`), scalar or a
#'   length-`nspp` vector. Default `NA` inherits `data_list$growth_sd_style`
#'   if present (so a refit keeps the original choice), otherwise `"WHAM"`.
#' @param linkages Optional named list of [linkage_spec()] objects
#'   keyed by parameter name (must be one of [GROWTH_LINKAGE_PARAMS]).
#'   The mean-growth keys (`K`, `L1`, `Linf`, `m`)
#'   accept arbitrary one-sided formulas and make that growth parameter
#'   year-varying (a per-year offset around its mean). The SD-endpoint keys
#'   (`sd_L1`, `sd_Linf`) only honor intercept-bearing
#'   formulas (typically `~ 1`), they thread `init`, `bounds`, and
#'   `priors` onto the growth SD-at-age, giving
#'   the SDs the same prior/fix/initial-value contract as the mean
#'   parameters. Slope rows on SD specs raise a warning and have no
#'   effect; slope-only formulas (`~ 0 + temp`) error.
#'
#' @return A list of switches defining the growth model.
#' @export
#'
#' @examples
#' \donttest{
#' # Sex-specific von Bertalanffy with temperature on K, by species + sex
#' build_growth(
#'   fun = "vonBertalanffy",   # or fun = 1
#'   linkages = list(
#'     K = linkage_spec(
#'       formula = ~ temp,
#'       by      = ~ species + sex,
#'       priors  = list(temp = normal(0, 1))
#'     )
#'   )
#' )
#' }
build_growth <- function(fun = "empirical",
                         growth_age_L1 = NA,
                         sd_plus_group = NA,
                         linkages = NULL) {
  fun <- .coerce_growth_fun(fun)
  # NA (the default) means "inherit": fit_mod() resolves it from
  # data_list$growth_sd_style if present, else the WHAM fallback -- exactly like
  # growth_age_L1. This keeps a refit that rebuilds growth via build_growth(fun=)
  # from silently resetting an SS3 model to WHAM.
  sd_style_int <- if (all(is.na(sd_plus_group))) {
    NA_integer_
  } else {
    .coerce_switch_arg(sd_plus_group, .GROWTH_SD_STYLE, "sd_plus_group")
  }
  linkages <- .validate_growth_linkages(linkages, fun)
  list(
    fun            = fun,
    linkages       = linkages,
    # Internal integer code consumed by the TMB template until the
    # linkage-driven path replaces it. Vectorized so per-species
    # growth functions (e.g. fun = c("vonBertalanffy", "Richards")) each
    # propagate their own code downstream.
    growth_model   = unname(.GROWTH_FUN_TO_INT[fun]),
    # Plus-group SD-at-age treatment. Like `fun`/`growth_model`, the canonical
    # string is kept under the argument name (for save_config round-trip) and
    # the int code (1 = WHAM, 2 = SS3) feeds the TMB template. Per-species.
    # NA_integer_ = inherit (resolved in fit_mod); the string field is NA too and
    # is dropped from a saved config, so unspecified models re-derive on load.
    sd_plus_group   = names(.GROWTH_SD_STYLE)[match(sd_style_int, .GROWTH_SD_STYLE)],
    growth_sd_style = sd_style_int,
    # VB anchor age (= age at which `l1` is the length). Matches SS3's
    # `Growth_Age_for_L1` ctl input. Scalar or length-nspp; pass NA
    # (default) to inherit max(0.5, minage[sp]) downstream so old
    # configurations stay unchanged at minage >= 1 and minage = 0
    # models pick up an SS3-consistent half-year anchor.
    growth_age_L1  = growth_age_L1
  )
}


#' Coerce a `fun` argument (string or integer) to canonical strings
#'
#' Accepts either a [GROWTH_FUNS] string or its integer code from
#' `.GROWTH_FUN_TO_INT`. Length-N vectors are allowed for per-species
#' specifications. Returns a character vector of canonical names.
#'
#' @keywords internal
#' @noRd
.coerce_growth_fun <- function(fun) {
  if (length(fun) == 0L) {
    stop("`fun` must have length >= 1", call. = FALSE)
  }
  if (is.numeric(fun)) {
    int <- as.integer(fun)
    hit <- match(int, .GROWTH_FUN_TO_INT)
    if (anyNA(hit)) {
      stop("integer `fun` must be one of: ",
           paste(.GROWTH_FUN_TO_INT, collapse = ", "),
           " (= ",
           paste(names(.GROWTH_FUN_TO_INT), collapse = "/"),
           ")", call. = FALSE)
    }
    return(names(.GROWTH_FUN_TO_INT)[hit])
  }
  if (is.character(fun)) {
    bad <- setdiff(fun, GROWTH_FUNS)
    if (length(bad) > 0L) {
      stop("unknown growth fun(s): ",
           paste(unique(bad), collapse = ", "),
           "; allowed: ", paste(GROWTH_FUNS, collapse = ", "),
           call. = FALSE)
    }
    return(fun)
  }
  stop("`fun` must be a string or integer; got ", class(fun)[1],
       call. = FALSE)
}


#' Validate and canonicalize the `linkages` argument of [build_growth()]
#'
#' Returns either `NULL` (no linkages) or a named list of fully-typed
#' `Rceattle_linkage_spec` objects with `param` filled in from the list
#' keys. Errors loudly on invalid param names so the user catches typos
#' at build time rather than at materialization time.
#'
#' @keywords internal
#' @noRd
.validate_growth_linkages <- function(linkages, fun) {
  linkages <- .validate_process_linkages(
    linkages, GROWTH_LINKAGE_PARAMS, "growth"
  )
  if (is.null(linkages)) return(NULL)
  if (any(fun == "empirical")) {
    warning("at least one species uses fun = 'empirical', which does ",
            "not consume linkages; the supplied specs will be ",
            "retained on the object but ignored at fit time for those ",
            "species.", call. = FALSE)
  }
  if (!all(fun == "Richards") && "m" %in% names(linkages)) {
    stop("linkages$m is only valid when every species uses ",
         "fun = 'Richards'; von Bertalanffy has no `m` parameter.",
         call. = FALSE)
  }
  # SD endpoints (log_sd_L1, log_sd_Linf) plug into `growth_log_sd`,
  # which has no year dim. Only intercept-bearing formulas (~ 1 plus
  # any covariates that the user wants to anchor on the intercept) are
  # honored: the intercept init/prior/bounds threads onto
  # `growth_log_sd`, and slope rows are silently dropped at the TMB
  # encoder (no offset path exists for SD). Block slope-only formulas
  # (`~ 0 + temp`) outright -- they'd leave the SD pinned at its
  # `build_params()` default with no escape.
  sd_keys <- intersect(names(linkages), names(.GROWTH_SD_PARAM_TO_INDEX))
  for (nm in sd_keys) {
    specs <- linkages[[nm]]
    if (inherits(specs, "Rceattle_linkage_spec")) specs <- list(specs)
    for (sp_i in specs) {
      f <- sp_i$formula
      rhs <- attr(stats::terms(f), "intercept")
      has_slope <- length(attr(stats::terms(f), "term.labels")) > 0L
      if (!isTRUE(rhs == 1L)) {
        stop(sprintf(
          "linkages$%s must use an intercept-bearing formula (e.g. ~ 1); ",
          nm),
          "year-varying offsets on growth SD are not wired through ",
          "growth.hpp. Drop the `0 +` from the formula, or move the ",
          "env effect onto K / L1 / Linf / m.",
          call. = FALSE)
      }
      if (has_slope) {
        warning(sprintf(
          "linkages$%s carries slope terms but growth SD has no year ",
          "dim in growth.hpp; the intercept init/prior/bounds will be ",
          "honored on `growth_log_sd`, but slope coefficients have no ",
          "effect on the SD.", nm), call. = FALSE)
      }
    }
  }
  linkages
}
