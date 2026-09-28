# =============================================================================
# build_data() -- code-first constructor for an Rceattle data list
# =============================================================================
#
# read_data() builds a data list from a fully-populated multi-sheet xlsx;
# build_data() lets a user assemble (or edit) one in R, supplying only the
# blocks a model actually uses and letting clean_data() default-fill the rest.
# It complements -- does not replace -- read_data()/write_data(): the returned
# object is the same bare list read_data() returns and round-trips through
# write_data() unchanged.

# Legacy top-level element names accepted by read_data(); build_data() maps them
# to the canonical names so scripts that used the old names port unchanged.
.RCE_DATA_ALIASES <- c(
  fsh_biom       = "catch_data",
  srv_biom       = "index_data",
  wt             = "weight",
  pmature        = "maturity",
  Pyrs           = "ration_data",
  UobsWtAge      = "diet_data",
  stom_prop_data = "diet_data",
  est_M1         = "M1_model"
)

#' Recognized top-level data-list element names
#'
#' Union of the requirement-table catalogue, the bioenergetics scalars, and the
#' auxiliary / configuration / legacy names that legitimately appear on a data
#' list. Used by build_data() to distinguish a genuine (rare) element from a
#' typo. Kept generous on purpose: a name NOT here is only flagged when it is a
#' near-miss of a known name.
#' @keywords internal
#' @noRd
.rce_known_data_names <- function() {
  unique(c(
    names(.rce_requirement_table()),
    .RCE_BIOENERGETICS_SCALARS,
    # auxiliary data / parameters
    "alpha_wt_len", "beta_wt_len", "aLW", "sigma_rec",
    "srr_prior_mean", "srr_prior_sd",
    "Diet_distribution", "Diet_comp_weights", "est_diet", "stom_tau",
    "UobsAge", "UobsWtAge", "minNByage", "MSSB0", "MSB0",
    # configuration / switches commonly stored on a data list
    "msmMode", "initMode", "avgnMode", "suitMode", "niter", "random_rec",
    "estDynamics", "est_M1", "M1_model", "M1_re", "debug", "proj_F",
    "HCR", "comp_offset", "model_config",
    # legacy top-level aliases (also normalised by read_data)
    names(.RCE_DATA_ALIASES)
  ))
}

#' Build an Rceattle data list in R
#'
#' Assemble (or edit) an Rceattle data list in R code rather than from an xlsx
#' workbook. Supply only the data
#' blocks a model uses, dimensions, biology (`weight`, `maturity`,
#' `sex_ratio`, `M1_base`), a `fleet_control`, and the observation tables
#' (`catch_data`, `index_data`, `comp_data`, ...), and the optional blocks a
#' single-species model does not need (`caal_data`, `emp_sel`, `diet_data`, ...)
#' are default-filled by [clean_data()]. The result is the same bare list
#' [read_data()] returns and round-trips through [write_data()] unchanged.
#'
#' Three entry points, which may be combined:
#' \itemize{
#'   \item **from blocks**: pass the elements as named arguments:
#'     `build_data(nspp = 1, styr = 1977, ..., fleet_control = fc, catch_data = catch)`.
#'   \item **from a file**: `build_data(file = "model.xlsx", projyr = 2060)`
#'     reads the workbook via [read_data()], then applies the overrides.
#'   \item **from an existing object**: `build_data(base = BS2017SS, projyr = 2060)`
#'     starts from a data list (e.g. a bundled dataset or a [combine_data()]
#'     result) and overrides the named blocks. This is the common
#'     copy-and-edit / combine-and-restamp workflow.
#' }
#'
#' Element names are checked against the recognized schema: an override whose
#' name is a near-miss of a known element (e.g. `maturty`) errors with a
#' suggestion, so a typo is caught here rather than surfacing much later in a
#' fit. Legacy top-level names (`fsh_biom`, `srv_biom`, `wt`, `pmature`, `Pyrs`)
#' are mapped to their canonical equivalents.
#'
#' Full validation runs at fit time, inside `fit_mod()`. `build_data()` runs only
#' a light presence pre-check (`.check = TRUE`) so a missing *required* block is
#' reported at construction with a clear message. The pre-check reads an attached
#' [model_config()], so a configuration held on the object is accounted for.
#' Requirements that depend on fit-time settings passed directly to [fit_mod()]
#' and stored nowhere on the data list are not knowable here and are left to
#' that fit-time validation; see [data_requirements()] to preview them.
#'
#' @param base Optional existing Rceattle data list to start from and override.
#' @param file Optional path to an Rceattle xlsx workbook; read via [read_data()]
#'   to form the starting object. Supply at most one of `base` / `file`.
#' @param ... Named top-level data-list elements to set or override (e.g.
#'   `nspp`, `styr`, `endyr`, `fleet_control`, `weight`, `maturity`,
#'   `catch_data`, `index_data`, `comp_data`). Every argument must be named.
#' @param .check Logical; run the presence pre-check (default `TRUE`). Set
#'   `FALSE` to assemble a deliberately partial object.
#'
#' @return An Rceattle data list (a bare `list`, as from [read_data()]).
#'
#' @seealso [read_data()], [clean_data()], [data_requirements()],
#'   [combine_data()], [fit_mod()].
#' @examples
#' # Copy-and-edit a bundled dataset.
#' dat <- build_data(base = BS2017SS, projyr = 2060)
#'
#' # Preview what a configuration needs before building it.
#' data_requirements(msmMode = 0)
#' @export
build_data <- function(base = NULL, file = NULL, ..., .check = TRUE) {

  # ---- Resolve the starting object ------------------------------------------
  if (!is.null(file)) {
    if (!is.null(base)) {
      stop("Supply at most one of `base` or `file`, not both.", call. = FALSE)
    }
    # `file` precedes `...`, so a *second* unnamed positional argument lands here
    # rather than in `...`. Reject a non-path to catch that (the first positional
    # is caught by the `base` type guard above) and to steer overrides toward the
    # required named form.
    if (!is.character(file) || length(file) != 1 || is.na(file)) {
      stop("`file` must be a single path to an .xlsx workbook. ",
           "Pass data blocks as named arguments, e.g. ",
           "build_data(fleet_control = fc, catch_data = catch).", call. = FALSE)
    }
    base <- Rceattle::read_data(file)
  }
  if (is.null(base)) {
    dl <- list()
  } else {
    # A data.frame is a list, so it would slip past a bare is.list() check; the
    # first unnamed positional argument lands in `base` (not `file`), so a stray
    # data frame passed positionally must be rejected here.
    if (is.data.frame(base) || !is.list(base)) {
      stop("`base` must be an Rceattle data list (a plain list), not a ",
           class(base)[1], ". Pass data blocks as named arguments, e.g. ",
           "build_data(fleet_control = fc, catch_data = catch).", call. = FALSE)
    }
    dl <- base
  }

  # ---- Collect, validate, and apply the named overrides ---------------------
  overrides <- list(...)
  if (length(overrides) > 0) {
    nms <- names(overrides)
    if (is.null(nms) || any(!nzchar(nms))) {
      stop("All data blocks passed to build_data() via `...` must be named.",
           call. = FALSE)
    }

    # Canonicalise each override name: map legacy aliases, catch typos. Matching
    # is case-sensitive (the model reads exact names like `maturity`), but a
    # case-variant of a known name (`Maturity`, `WT`) is a typo, not a novel
    # element -- flag it with the correctly-cased suggestion rather than letting
    # it through as junk or mis-suggesting a distant near-miss.
    known   <- .rce_known_data_names()
    aliases <- .RCE_DATA_ALIASES
    canon   <- character(length(nms))
    for (i in seq_along(nms)) {
      nm <- nms[i]
      if (nm %in% names(aliases)) {          # exact legacy alias
        canon[i] <- aliases[[nm]]
      } else if (nm %in% known) {            # exact canonical name
        canon[i] <- nm
      } else {                               # unknown: alias/known case-variant?
        ai <- match(tolower(nm), tolower(names(aliases)))
        ki <- match(tolower(nm), tolower(known))
        if (!is.na(ai)) {
          stop("Unrecognised data element '", nm, "'. Did you mean '",
               names(aliases)[ai], "' (names are case-sensitive)?", call. = FALSE)
        } else if (!is.na(ki)) {
          stop("Unrecognised data element '", nm, "'. Did you mean '",
               known[ki], "' (names are case-sensitive)?", call. = FALSE)
        }
        d <- utils::adist(nm, known, ignore.case = TRUE)[1, ]   # near-miss typo?
        close <- known[d > 0 & d <= 2]
        if (length(close) > 0) {
          stop("Unrecognised data element '", nm, "'. Did you mean: ",
               paste(sort(close), collapse = ", "), "?", call. = FALSE)
        }
        canon[i] <- nm                       # genuinely novel: pass through
      }
    }

    # Duplicate check on the CANONICAL names, so a collision introduced by alias
    # mapping (e.g. `wt` and `weight`, or two diet aliases) is caught rather than
    # silently resolved last-wins.
    if (any(duplicated(canon))) {
      dup <- unique(canon[duplicated(canon)])
      stop("Multiple arguments resolve to the same data element(s): ",
           paste(dup, collapse = ", "),
           ". Supply each block once (mind legacy aliases such as ",
           "wt/weight, fsh_biom/catch_data).", call. = FALSE)
    }

    for (i in seq_along(overrides)) dl[[canon[i]]] <- overrides[[i]]
  }

  # ---- Thin eager presence pre-check ----------------------------------------
  # Run BEFORE clean_data(): a missing required *core* block (weight,
  # fleet_control, ...) would otherwise crash clean_data()'s unguarded
  # dereferences with a cryptic error before this friendly message can fire.
  # Required core blocks are never default-filled, so checking pre-fill is
  # equivalent; the optional blocks clean_data fills are, by definition, not
  # required for this configuration.
  if (isTRUE(.check)) .rce_build_precheck(dl)

  # ---- Default-fill the optional blocks -------------------------------------
  # With .check = TRUE the required core is guaranteed present, so clean_data()
  # will not trip. With .check = FALSE the object may be deliberately partial;
  # keep whatever assembled if clean_data() cannot complete.
  dl <- tryCatch(Rceattle::clean_data(dl), error = function(e) {
    if (isTRUE(.check)) stop(e)
    dl
  })

  .rce_as_data(dl)
}

#' Presence pre-check for build_data()
#'
#' Reports data-list elements that are Required for the object's stored
#' configuration but absent, reusing the requirement table's classification so
#' the message can never diverge from data_check(). A subset of data_check()
#' (presence only), deliberately, the authoritative validation runs at fit
#' time.
#'
#' @param dl A (clean_data-filled) data list.
#' @return Invisibly `TRUE`; errors listing the missing required elements.
#' @keywords internal
#' @noRd
.rce_build_precheck <- function(dl) {
  # Classify against the configuration the fit will actually use: an attached
  # model_config() supplies the switches, and fleet_control may still hold
  # integer codes at this point (this runs before switch_check()).
  dl <- .rce_classify_view(dl)
  present <- function(x) {
    if (is.null(x)) return(FALSE)
    # A data.frame is present when it has rows; any other element (scalar,
    # vector) when it has length -- an empty scalar (numeric(0)) is NOT present
    # and would otherwise crash clean_data() with a cryptic error.
    if (is.data.frame(x)) nrow(x) > 0 else length(x) > 0
  }
  tbl <- .rce_requirement_table()
  missing_req <- character(0)
  conditions  <- character(0)
  for (nm in names(tbl)) {
    row <- tbl[[nm]]
    if (nm == "bioenergetics") next   # a group; checked at fit time
    if (identical(.rce_classify(row, dl), "Required") && !present(dl[[nm]])) {
      missing_req <- c(missing_req, nm)
      conditions  <- c(conditions,
                       if (isTRUE(row$always_required)) "always required"
                       else row$condition_label %||% "required")
    }
  }
  if (length(missing_req) > 0) {
    stop("build_data(): the following required data element(s) are missing:\n",
         paste0("  - ", missing_req, "  (", conditions, ")", collapse = "\n"),
         "\nSupply them, or see data_requirements() for the full picture. ",
         "(Pass .check = FALSE to build a partial object.)", call. = FALSE)
  }
  invisible(TRUE)
}
