#' Allowed M-parameter names for `linkages` in [build_M1()]
#'
#' Natural-scale names of the underlying natural-mortality
#' parameters that the linkage system can address. Currently just
#' `M1`, with the default log link the offset is added to M1 on the
#' log scale (applied across all ages unless the linkage row
#' pins a specific `age_bin`).
#'
#' Note: predation mortality `M2` is a derived quantity in CEATTLE
#' (a function of predator abundance, suitability, and ration), not
#' a parameter. There is no `M2` linkage target; environmental
#' effects on predation are mediated upstream via recruitment,
#' growth, suitability, or ration inputs.
#'
#' @keywords internal
M_LINKAGE_PARAMS <- c("M1")


#' String<->integer mapping for `M1_model` in [build_M1()]
#'
#' Either form is accepted by [build_M1()]; the canonical integer
#' code is what the TMB template ultimately consumes. The
#' env-driven integer codes 4 and 5 still work with a deprecation
#' warning, their structural part is identical to 1 and 2
#' respectively, and the env effect is expressed via the
#' \code{linkages} argument to [build_M1()] (see
#' `vignette("environmental-linkages-and-priors")`). No string alias is offered
#' for 4 or 5 to discourage their use in new code.
#'
#' @keywords internal
.M1_MODELS <- c(
  fixed             = 0L,
  sex_age_invariant = 1L,
  sex_specific      = 2L,
  sex_age_specific  = 3L
)


#' Deprecated env-driven `M1_model` integer codes.
#' @keywords internal
#' @noRd
.M1_DEPRECATED_MODELS <- c(4L, 5L)


#' String<->integer mapping for `M1_re` in [build_M1()]
#'
#' Either form is accepted by [build_M1()]; the canonical integer
#' code is what the TMB template ultimately consumes. The
#' env-driven options (4, 5) on M1_model can also be expressed
#' through the `linkages` argument to [build_M1()]; the integer
#' codes continue to work for backwards compatibility.
#'
#' @keywords internal
.M1_RES <- c(
  none         = 0L,
  iid_age      = 1L,
  iid_year     = 2L,
  iid_age_year = 3L,
  ar1_age      = 4L,
  ar1_year     = 5L,
  ar1_age_year = 6L
)


#' Coerce an `M1_model` or `M1_re` value to canonical integer code(s).
#'
#' Accepts either a string from the corresponding map (length-1 or
#' length-N for per-species values) or a length-1/length-N integer
#' vector that is already in the allowed set. Errors loudly on
#' anything else.
#'
#' For `M1_model` only, the env-driven integer codes 4
#' and 5 (controlled by `M1_indices`) are accepted for backwards
#' compatibility but emit a soft-deprecation warning pointing users
#' at the linkage table, see [build_M1()].
#'
#' @keywords internal
#' @noRd
.coerce_M1_arg <- function(x, map, what) {
  deprecated <- if (identical(what, "M1_model")) .M1_DEPRECATED_MODELS else integer(0)
  .coerce_switch_arg(
    x, map = map, what = what,
    deprecated = deprecated,
    warn_fn = function(int) .warn_M1_model_deprecation(int),
    length_exact_one = FALSE)
}


#' @keywords internal
#' @noRd
.warn_M1_model_deprecation <- function(int) {
  hits <- intersect(int, .M1_DEPRECATED_MODELS)
  warning(
    "M1_model %in% c(", paste(hits, collapse = ", "),
    ") is soft-deprecated: the structural part of mode 4 is ",
    "identical to mode 1, and mode 5 to mode 2 (or 1 in single-",
    "sex models). The environmental effect formerly driven by ",
    "`M1_indices` is now better expressed through the linkages ",
    "argument to build_M1():\n\n",
    "  build_M1(M1_model = c(1, 2, 1),\n",
    "           linkages = list(M1 = linkage_spec(\n",
    "             formula = ~ <env_col>, by = ~species,\n",
    "             species = <which species had the env effect>)))\n\n",
    "See vignette('environmental-linkages-and-priors').",
    call. = FALSE
  )
}


#' @keywords internal
#' @noRd
.warn_M1_indices_deprecation <- function() {
  warning(
    "`M1_indices` is soft-deprecated. The same environmental ",
    "effect is now better expressed through the linkages argument ",
    "to build_M1():\n\n",
    "  build_M1(M1_model = ...,\n",
    "           linkages = list(M1 = linkage_spec(\n",
    "             formula = ~ <env_col>, by = ~species)))\n\n",
    "Both paths add additively to log_M1 on the log scale, so do ",
    "NOT supply both for the same coefficient or you will ",
    "double-count. See vignette('environmental-linkages-and-priors').",
    call. = FALSE
  )
}


#' Specify the residual natural mortality (M1) model for Rceattle
#'
#' @param M1_model Vector or scalar specifying the M1 structural fixed-
#'   effects model. Either an integer code or the equivalent string
#'   alias (both forms are accepted; the integer code is canonical):
#'   * `0` / `"fixed"`, use the input `M1_base` (no estimation).
#'   * `1` / `"sex_age_invariant"`, estimate one `M1_{spp}`.
#'   * `2` / `"sex_specific"`, estimate `M1_{spp, sex}`.
#'   * `3` / `"sex_age_specific"`, estimate `M1_{spp, sex, age}`.
#'   * `4`, `5`, soft-deprecated env-driven codes; use the
#'     `linkages` argument instead. See `vignette("environmental-linkages-and-priors")`.
#' @param M1_re Vector or scalar specifying the M1 random-effects
#'   model. Either an integer code or the equivalent string alias:
#'   `0` / `"none"`, `1` / `"iid_age"`, `2` / `"iid_year"`,
#'   `3` / `"iid_age_year"`, `4` / `"ar1_age"`, `5` / `"ar1_year"`,
#'   `6` / `"ar1_age_year"`.
#' @param updateM1 If using initial parameters, use M1 fixed effects
#'   from data (`M1_base`) instead. Default `FALSE`.
#' @param M1_use_prior Vector or scalar; if `TRUE` and `M2_use_prior` is `FALSE`, apply the
#'   lognormal `M_prior` / `M_prior_sd` to `M1` directly, which is refused alongside a prior
#'   on that species' `M1` linkage intercept.
#' @param M2_use_prior Vector or scalar; if `TRUE`, apply the
#'   lognormal prior to `M1 + M2` in multi-species models.
#' @param M_prior Natural-scale mean of the lognormal prior on M (its median when `bias_adjust_proc = FALSE`).
#' @param M_prior_sd SD (log-scale) of the lognormal prior on M.
#' @param M1_indices Soft-deprecated. Vector of column indices into
#'   `env_data` (excluding `Year`) for environmentally linked M1 when
#'   `M1_model %in% c(4, 5)`. Use the `linkages` argument instead;
#'   see \code{vignette("environmental-linkages-and-priors")}.
#' @param linkages Optional named list of [linkage_spec()] objects
#'   keyed by M parameter name (currently the only valid key is
#'   `"M1"`). Each spec describes how `M1` depends on
#'   environmental covariates and on stratifying factors (species,
#'   sex, age). The offset enters additively (on the log scale)
#'   inside the `M1_at_age` compute. A row's `age_bin == NA`
#'   broadcasts the offset across ages; specific values pin it to
#'   that age slice.
#'
#' @return A list of switches for defining the M1 model.
#' @export
#'
#' @examples
#' \donttest{
#' # Sex/age-invariant M with a temperature linkage on M1
#' build_M1(
#'   M1_model = "sex_age_invariant",
#'   linkages = list(
#'     M1 = linkage_spec(
#'       formula = ~ temp,
#'       by      = ~ species,
#'       priors  = list(temp = normal(0, 0.5))
#'     )
#'   )
#' )
#' }
build_M1 <- function(M1_model = 0,
                     M1_re = 0,
                     updateM1 = FALSE,
                     M1_use_prior = FALSE,
                     M2_use_prior = FALSE,
                     M_prior = 0.40,
                     M_prior_sd = 0.35,
                     M1_indices = NA,
                     linkages = NULL){
  M1_model <- .coerce_M1_arg(M1_model, .M1_MODELS, "M1_model")
  M1_re    <- .coerce_M1_arg(M1_re,    .M1_RES,    "M1_re")
  linkages <- .validate_M_linkages(linkages)
  # `M1_indices` is soft-deprecated in favour of `linkages = list(log_M
  # = ...)`. NA / NA_real_ / NA_integer_ are all "not supplied". Any
  # other value triggers the deprecation note.
  if (!(length(M1_indices) == 1L && is.na(M1_indices))) {
    .warn_M1_indices_deprecation()
  }
  list(
    M1_model     = M1_model,
    M1_re        = M1_re,
    updateM1     = updateM1,
    M1_use_prior = M1_use_prior,
    M2_use_prior = M2_use_prior,
    M_prior      = M_prior,
    M_prior_sd   = M_prior_sd,
    M1_indices   = M1_indices,
    linkages     = linkages
  )
}


#' Validate and canonicalize the `linkages` argument of [build_M1()]
#'
#' Returns either `NULL` (no linkages) or a named list of
#' `Rceattle_linkage_spec` objects (or lists thereof) with `param`
#' filled in from the list keys. Errors loudly on invalid param
#' names so the user catches typos at build time.
#'
#' @keywords internal
#' @noRd
.validate_M_linkages <- function(linkages) {
  .validate_process_linkages(linkages, M_LINKAGE_PARAMS, "M")
}


# An M1 linkage intercept prior and M1_use_prior both penalize that species' log M1, so the
# prior would count twice. A row with no species targets species 1, as in the template.
.check_M_linkage_prior <- function(linkage_table, M1_use_prior, M2_use_prior, spnames) {
  if (is.null(linkage_table) || nrow(linkage_table) == 0L) return(invisible())
  pri <- linkage_table[linkage_table$process == "M" &
                         linkage_table$design_col == "(Intercept)" &
                         !is.na(linkage_table$prior_family) &
                         linkage_table$prior_family != "none", , drop = FALSE]
  if (nrow(pri) == 0L) return(invisible())
  sp <- unique(ifelse(is.na(pri$species), 1L, as.integer(pri$species)))
  # With M2_use_prior on, the M prior is on total M, not log M1.
  both <- sp[(M1_use_prior[sp] %in% 1) & (M2_use_prior[sp] %in% 0)]
  if (length(both) > 0L) {
    stop(sprintf(paste0(
      "species %s: M1_use_prior = TRUE and a prior on the M1 linkage's intercept are ",
      "both priors on the same M1, so it would be counted twice. Keep one: set ",
      "M1_use_prior = FALSE, or drop the intercept prior."),
      paste(spnames[both] %||% both, collapse = ", ")), call. = FALSE)
  }
  invisible()
}
