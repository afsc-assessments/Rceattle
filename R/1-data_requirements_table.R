# =============================================================================
# Declarative data-requirement table
# =============================================================================
#
# The *conditional presence requirements* of the ~40-element Rceattle data
# list: e.g. "diet_data is required when msmMode > 0", "NByageFixed is required
# when estDynamics > 0". Those conditions are exactly what a user needs to know
# up front, and what `data_requirements()` reports and `build_data()`
# pre-checks. So each condition lives once, here, and `data_check()` and
# `data_requirements()` both *consume* it rather than restating it.
#
# Scope, deliberately narrow: this table drives only the **pure
# presence-requirement** gates in `data_check()`. Dimension / value /
# referential / structural checks and the
# two mirroring-dependent *adequacy* gates (index_cov-MVN, comp/caal-vs-estimated-
# selectivity) stay imperative in `data_check()` -- they depend on per-fleet row
# counts and cross-fleet lookups that do not reduce to a declarative row. Rows
# marked `driven = FALSE` are read by `data_requirements()` for classification
# only and are NOT consulted by `data_check()`; the authoritative check for those
# stays in `data_check()`.
#
# Each row is a list with fields:
#   element        chr  data_list element name (the row's identity / lookup key).
#   category       chr  grouping for display ("dimensions", "biology", "fishery",
#                       "composition", "predation", "environment").
#   always_required lgl TRUE for the core backbone the model always needs
#                       (dereferenced unconditionally); such rows carry no
#                       condition.
#   required_when  fn   function(dl) -> logical; TRUE when the element is
#                       required for this configuration. Evaluated with isTRUE()
#                       by the callers, so a length-0 / NA result reads as "not
#                       required" (robust to a not-yet-defaulted switch).
#   ignored_when   fn   function(dl) -> logical; TRUE when the feature that would
#                       consume the element is switched OFF, so the element is
#                       neither required nor used. Optional; defaults to never.
#   condition_label chr human-readable form of `required_when`, for reports.
#   optional_status chr when neither required nor ignored: "defaulted"
#                       (clean_data fills a safe default -> reported "Optional")
#                       or "none".
#   default_label  chr  human description of the clean_data default, for reports.
#   driven         lgl  TRUE if `data_check()` emits this row's requirement via
#                       the evaluators below; FALSE = classification-only.
#   severity       chr  ("error" | "message") how data_check surfaces an unmet
#                       requirement. Only meaningful when driven = TRUE.
#   adequate       fn   function(dl) -> logical; TRUE when the element is present
#                       and adequate enough to satisfy the *presence* gate (the
#                       separate dimension/value checks stay in data_check).
#                       Only required when driven = TRUE.
#   message        fn   function(dl) -> chr; the exact text data_check emits when
#                       required_when && !adequate. Only required when driven.
#
# NOTE: `required_when` reproduces each gate's original guard. Where the original
# guard would error on a NULL switch (e.g. the diet gate's bare `msmMode > 0`),
# the isTRUE() wrapping here reads NULL as "not required" instead -- which is
# unreachable in the real pipeline (switch_check defaults msmMode to 0 before
# data_check runs) and strictly more robust. Flagged for review.

# Canonical presence predicate: a data.frame element counts as "present" only
# when it is non-NULL and has at least one row. Mirrors the `has_data()` helper
# local to data_check(); defined here so the table and data_check share one
# definition.
.rce_has_data <- function(df) !is.null(df) && nrow(df) > 0

# The 12 bioenergetics scalars are required together (one gate, one message
# listing whichever are missing), so they share a single group row keyed on
# "bioenergetics".
.RCE_BIOENERGETICS_SCALARS <- c(
  "Ceq", "Cindex", "Pvalue", "fday", "CA", "CB",
  "Qc", "Tco", "Tcm", "Tcl", "CK1", "CK4"
)

# TRUE when any fleet requests an MVN survey covariance. Tolerates the integer
# code as well as the canonical string, for the same reason .rce_is_fixed_sel()
# does: this layer also classifies objects switch_check() has not yet seen.
.rce_any_mvn <- function(dl) {
  fc <- dl$fleet_control
  if (is.null(fc) || !"Index_distribution" %in% colnames(fc)) return(FALSE)
  v <- as.character(fc$Index_distribution)
  num <- suppressWarnings(as.numeric(v))
  v[!is.na(num)] <- as.character(num[!is.na(num)])
  any(!is.na(v) & v %in% c("MVN", "MVNORM", "1", "2"))
}

# TRUE for each Selectivity entry naming the "Fixed" (empirical) form, whether
# the column still holds the integer code or the canonical string. switch_check()
# converts 0 -> "Fixed", but this layer also classifies raw build_data() objects
# that have not been through it, and workbook reads can leave the code as a
# numeric-looking string. Mirrors the code/string tolerance .rce_any_mvn() has.
.rce_is_fixed_sel <- function(x) {
  num <- suppressWarnings(as.numeric(as.character(x)))
  !is.na(x) & (as.character(x) == "Fixed" | (!is.na(num) & num == 0))
}

# TRUE when any fleet fixes selectivity to an empirical curve.
.rce_any_fixed_sel <- function(dl) {
  fc <- dl$fleet_control
  !is.null(fc) && "Selectivity" %in% colnames(fc) && any(.rce_is_fixed_sel(fc$Selectivity))
}

# Normalize a data list for requirement classification.
#
# A model_config() slot carries the model switches, but the requirement
# predicates read the top-level slots. fit_mod() resolves the config and writes
# the result back onto data_list, so the requirement report has to see the same
# object the fit will, or print() and data_requirements() disagree about the
# same data list.
#
# Overlay, not fill: an attached config wins over a stored top-level value,
# matching fit_mod(). Scope is limited to what model_config() actually carries --
# it has no estDynamics or Ceq field, so NByageFixed and the bioenergetics
# inputs still classify from the top-level slots alone.
.rce_classify_view <- function(dl) {
  cfg <- dl$model_config
  # Absent, or a malformed non-list slot: classify from the top-level switches.
  # build_data() does not type-check a `model_config` override, so an atomic can
  # reach here and `cfg[[nm]]` would abort with "subscript out of bounds".
  if (!is.list(cfg)) return(dl)
  for (nm in c("msmMode", "initMode", "avgnMode", "suitMode", "niter")) {
    if (!is.null(cfg[[nm]])) dl[[nm]] <- cfg[[nm]]
  }
  # fit_mod() overwrites growth_model from growthFun unconditionally, so a
  # default model_config() (empirical growth) genuinely does make caal_data
  # optional. Match the fit rather than second-guessing it.
  if (!is.null(cfg$growthFun$growth_model)) dl$growth_model <- cfg$growthFun$growth_model
  dl
}

#' Full data-element catalogue: requirement conditions + classification metadata
#'
#' @return A named list of element rows (see file header for the schema), keyed
#'   by the row `element`. Consumed by `data_check()` (the `driven` rows) and by
#'   `data_requirements()` (all rows).
#' @keywords internal
#' @noRd
.rce_requirement_table <- function() {
  # ---- Conditional rows (the six data_check-driven gates + two adequacy) -----
  rows <- list(

    # ---- Predation / diet (msmMode > 0) -----------------------------------
    list(
      element        = "diet_data",
      category       = "predation",
      required_when  = function(dl) dl$msmMode > 0,
      ignored_when   = function(dl) !isTRUE(dl$msmMode > 0),
      condition_label = "msmMode > 0 (multispecies)",
      driven         = TRUE,
      severity       = "error",
      adequate       = function(dl) .rce_has_data(dl$diet_data),
      message        = function(dl) "No diet data included",
      optional_status = "defaulted",
      default_label  = "empty diet_data (clean_data)"
    ),

    list(
      element        = "ration_data",
      category       = "predation",
      required_when  = function(dl) dl$msmMode > 0,
      ignored_when   = function(dl) !isTRUE(dl$msmMode > 0),
      condition_label = "msmMode > 0 (multispecies)",
      driven         = TRUE,
      severity       = "message",
      adequate       = function(dl) .rce_has_data(dl$ration_data),
      message        = function(dl) "No ration data",
      optional_status = "defaulted",
      default_label  = "empty ration_data (clean_data)"
    ),

    list(
      element        = "bioenergetics",
      category       = "predation",
      required_when  = function(dl) !is.null(dl$msmMode) && dl$msmMode > 0,
      ignored_when   = function(dl) !isTRUE(dl$msmMode > 0),
      condition_label = "msmMode > 0 (multispecies)",
      driven         = TRUE,
      severity       = "error",
      adequate       = function(dl) {
        all(vapply(.RCE_BIOENERGETICS_SCALARS, function(nm) {
          !is.null(dl[[nm]]) && length(dl[[nm]]) == dl$nspp
        }, logical(1)))
      },
      message        = function(dl) {
        bad <- .RCE_BIOENERGETICS_SCALARS[vapply(.RCE_BIOENERGETICS_SCALARS,
          function(nm) is.null(dl[[nm]]) || length(dl[[nm]]) != dl$nspp,
          logical(1))]
        paste0(
          "msmMode > 0 requires bioenergetics scalars of length nspp = ",
          dl$nspp, "; missing or wrong length: ", paste(bad, collapse = ", ")
        )
      },
      optional_status = "defaulted",
      default_label  = "safe sentinels in single-species mode (switch_check)"
    ),

    # ---- Fixed selectivity requires empirical selectivity ------------------
    list(
      element        = "emp_sel",
      category       = "fishery",
      required_when  = function(dl) .rce_any_fixed_sel(dl),
      ignored_when   = function(dl) !.rce_any_fixed_sel(dl),
      condition_label = "any fleet Selectivity == 'Fixed'",
      driven         = TRUE,
      severity       = "error",
      adequate       = function(dl) .rce_has_data(dl$emp_sel),
      message        = function(dl) {
        fc <- dl$fleet_control
        fixed_flts <- fc[.rce_is_fixed_sel(fc$Selectivity), , drop = FALSE]
        paste0(
          "Fleet(s) with Selectivity = 'Fixed' (",
          paste(fixed_flts$Fleet_name, collapse = ", "),
          ") require emp_sel data"
        )
      },
      optional_status = "defaulted",
      default_label  = "empty emp_sel (clean_data)"
    ),

    # ---- Fixed numbers-at-age requires NByageFixed -------------------------
    list(
      element        = "NByageFixed",
      category       = "dimensions",
      required_when  = function(dl) any(dl$estDynamics > 0),
      ignored_when   = function(dl) !isTRUE(any(dl$estDynamics > 0)),
      condition_label = "any estDynamics > 0 (fixed numbers-at-age)",
      driven         = TRUE,
      severity       = "error",
      adequate       = function(dl) .rce_has_data(dl$NByageFixed),
      message        = function(dl)
        "estDynamics > 0 requires NByageFixed data to be provided",
      optional_status = "defaulted",
      default_label  = "empty NByageFixed (clean_data)"
    ),

    # ---- Growth estimation requires conditional age-at-length --------------
    # caal_data is ALSO required (as an *adequacy* gate, kept imperative in
    # data_check) for fleets with estimated selectivity and no comp data; that
    # cross-fleet, mirroring-aware check does not reduce to a declarative row.
    # When growth is fixed, caal_data is Optional (used for CAAL composition if
    # supplied), so it carries no ignored_when.
    list(
      element        = "caal_data",
      category       = "composition",
      required_when  = function(dl) any(dl$growth_model > 0),
      condition_label = "any growth_model > 0 (estimated growth)",
      driven         = TRUE,
      severity       = "error",
      adequate       = function(dl) .rce_has_data(dl$caal_data),
      message        = function(dl)
        "Growth estimation (growth_model > 0) requires caal_data to be provided",
      optional_status = "defaulted",
      default_label  = "empty caal_data (clean_data)"
    ),

    # ---- Environmental index required by temperature-dependent consumption -
    # Classification-only: the real check (per-species Cindex column count) is a
    # dimension gate and stays imperative in data_check(). Optional otherwise
    # (clean_data fills a Year-only frame).
    list(
      element        = "env_data",
      category       = "environment",
      required_when  = function(dl) !is.null(dl$Ceq) && any(dl$Ceq > 1),
      condition_label = "any Ceq > 1 (temperature-dependent consumption)",
      driven         = FALSE,
      optional_status = "defaulted",
      default_label  = "Year-only env_data (clean_data)"
    ),

    # ---- MVN survey covariance --------------------------------------------
    # Classification-only: presence + square + dimension + symmetry all live in
    # data_check() (they need the per-fleet fitted-observation count). Ignored
    # unless a fleet requests MVN.
    list(
      element        = "index_cov",
      category       = "fishery",
      required_when  = function(dl) .rce_any_mvn(dl),
      ignored_when   = function(dl) !.rce_any_mvn(dl),
      condition_label = "any fleet Index_distribution == 'MVN'",
      driven         = FALSE,
      optional_status = "defaulted",
      default_label  = "empty list (clean_data)"
    )
  )

  # ---- Core backbone: always required, no condition ------------------------
  core <- list(
    nspp        = "dimensions", styr        = "dimensions",
    endyr       = "dimensions", projyr      = "dimensions",
    spnames     = "dimensions", nsex        = "dimensions",
    spawn_month = "dimensions", nages       = "dimensions",
    minage      = "dimensions", nlengths    = "dimensions",
    other_food  = "dimensions",
    pop_wt_index = "biology",   ssb_wt_index = "biology",
    pop_age_transition_index = "biology",
    weight   = "biology", maturity = "biology", sex_ratio = "biology",
    M1_base  = "biology", age_trans_matrix = "biology", age_error = "biology",
    fleet_control = "fishery", catch_data = "fishery"
  )
  for (nm in names(core)) {
    rows[[length(rows) + 1L]] <- list(
      element = nm, category = core[[nm]], always_required = TRUE,
      condition_label = "always required", driven = FALSE,
      optional_status = "none"
    )
  }

  # ---- Optional-and-defaulted (used if supplied, no condition) -------------
  optional <- list(
    index_data = list(cat = "fishery",     def = "no survey index fitted"),
    comp_data  = list(cat = "composition", def = "empty comp_data (clean_data)")
  )
  for (nm in names(optional)) {
    rows[[length(rows) + 1L]] <- list(
      element = nm, category = optional[[nm]]$cat,
      condition_label = "optional (used if supplied)", driven = FALSE,
      optional_status = "defaulted", default_label = optional[[nm]]$def
    )
  }

  names(rows) <- vapply(rows, function(r) r$element, character(1))
  rows
}

#' Classify one element's status for a given (normalized) data list
#'
#' @param row One requirement-table row.
#' @param dl A data list whose switches are populated (ideally post switch_check).
#' @return "Required", "Optional", or "Ignored".
#' @keywords internal
#' @noRd
.rce_classify <- function(row, dl) {
  if (isTRUE(row$always_required)) return("Required")
  if (!is.null(row$required_when) && isTRUE(row$required_when(dl))) return("Required")
  if (!is.null(row$ignored_when)  && isTRUE(row$ignored_when(dl)))  return("Ignored")
  if (identical(row$optional_status, "defaulted")) return("Optional")
  "Ignored"
}

# -----------------------------------------------------------------------------
# Evaluators used by data_check() at each driven gate site. Placing the call at
# the original site preserves the accumulation order of `errors`; the table
# supplies the condition, adequacy test, and message text.
# -----------------------------------------------------------------------------

#' Presence-requirement check for a `severity = "error"` row
#'
#' @param dl Rceattle data list.
#' @param element Row key (data_list element name).
#' @return `character(0)` if the requirement is not active or is satisfied;
#'   otherwise the single error string for the caller to accumulate.
#' @keywords internal
#' @noRd
.rce_check_presence <- function(dl, element) {
  row <- .rce_requirement_table()[[element]]
  if (is.null(row) || !isTRUE(row$driven) || row$severity != "error") {
    stop("Internal: no driven error requirement '", element, "'")
  }
  if (isTRUE(row$required_when(dl)) && !isTRUE(row$adequate(dl))) {
    row$message(dl)
  } else {
    character(0)
  }
}

#' Presence-requirement notice for a `severity = "message"` row
#'
#' Emits the row message via [message()] when the requirement is active and
#' unmet (matching the historical `message("No ration data")` behavior).
#'
#' @param dl Rceattle data list.
#' @param element Row key.
#' @return Invisibly `NULL`.
#' @keywords internal
#' @noRd
.rce_notify_presence <- function(dl, element) {
  row <- .rce_requirement_table()[[element]]
  if (is.null(row) || !isTRUE(row$driven) || row$severity != "message") {
    stop("Internal: no driven message requirement '", element, "'")
  }
  if (isTRUE(row$required_when(dl)) && !isTRUE(row$adequate(dl))) {
    message(row$message(dl))
  }
  invisible(NULL)
}
