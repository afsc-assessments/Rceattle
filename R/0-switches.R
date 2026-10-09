# =============================================================================
# Configuration "switch" lifecycle for Rceattle
# =============================================================================
# A switch (selectivity form, catchability form, composition likelihood, HCR,
# init mode, ...) may be supplied as a string ("Logistic") or as the integer
# code the TMB template reads (1L); the maps below are the single source of
# truth for that correspondence. Inside fit_mod() the four functions run in
# this order:
#
#   switch_check()      fill missing switches, normalise to canonical strings
#   revert_switches()   integer code -> canonical string, for older data files
#   validate_switches() refuse an unknown code/string; from data_check()
#   convert_switches()  canonical string -> integer; from rearrange_data()
#
# revert_switches() and convert_switches() are inverses, not duplicates, and
# their per-column logic differs (Time_varying_q doubles as an environmental
# index column), so the two are kept in step by eye. The switch system is also
# written up in the developer guide.
# =============================================================================


#' Rceattle configuration switch maps
#'
#' Internal lookup tables that translate the human-readable switch strings a
#' user supplies (e.g. \code{"Logistic"}, \code{"NPFMC"}) into the integer
#' codes the TMB template consumes (and back again, via
#' \code{revert_switches()}). They are the single source of truth for the set
#' of valid switch values; \code{validate_switches()} errors on anything not
#' listed here. See the lifecycle note at the top of \code{R/0-switches.R} for
#' how the maps are applied.
#'
#' @details
#' \describe{
#'   \item{\code{sel_map}}{Selectivity form (\code{fleet_control$Selectivity}).}
#'   \item{\code{tv_sel_map}}{Time-varying selectivity structure
#'     (\code{fleet_control$Time_varying_sel}).}
#'   \item{\code{sel_norm_scope_map}}{Whether selectivity normalization pools its
#'     reference across sexes (\code{fleet_control$Sel_norm_scope}).}
#'   \item{\code{q_map}}{Catchability form (\code{fleet_control$Catchability}).}
#'   \item{\code{tv_q_map}}{Time-varying catchability structure
#'     (\code{fleet_control$Time_varying_q}).}
#'   \item{\code{comp_loglike_map}}{Composition likelihood
#'     (\code{fleet_control$Comp_distribution} and \code{CAAL_distribution}).}
#'   \item{\code{fleet_map}}{Fleet type (\code{fleet_control$Fleet_type}).}
#'   \item{\code{initMode_map}}{Initial age-structure mode
#'     (\code{data_list$initMode}; see \code{\link{fit_mod}}).}
#'   \item{\code{suitMode_map}}{Predator-prey suitability mode, per predator
#'     (\code{data_list$suitMode}).}
#'   \item{\code{hcr_map}}{Harvest control rule (see \code{\link{build_hcr}}).}
#' }
#'
#' @name rceattle_switch_maps
#' @keywords internal
NULL

# --- Forward Mappings (String -> Integer) ---
# - Maintains backwards compatibility
sel_map <- c(
  "Fixed" = 0,
  "Logistic" = 1,
  "NonParametric" = 2,
  "DoubleLogistic" = 3,
  "DescendingLogistic" = 4,
  "Hake" = 5,
  "2DAR1" = 6,
  "3DAR1" = 7,
  "DoubleNormal" = 8,
  "NonParametricPM" = 9,  # Ianelli non-parametric, ADMB AMAK ("pm") selectivity penalty
  "LogisticPM" = 11,      # ADMB AMAK ("pm") BTS: logistic (multiplicative inflection/slope devs) + free age-1 log-selectivity
  # Ianelli base curve whose deviations have a proper density, so random_sel can
  # integrate them. Time_varying_sel picks the structure: "Off" (none), "IID" or
  # "RandomWalk". Codes 10, 12 and 14 are unused.
  "NonParametricIntegrable" = 13,
  "DoubleNormalSS3"  = 15  # Stock Synthesis size pattern 24: six-parameter double normal, own parameter array (sel_dn6)
)

# Which selectivity forms read each Sel_curve_pen slot as a WEIGHT, and the sd
# column that sets the same weight safely. A weight multiplies a squared
# deviation, so a negative one rewards it. "2DAR1" (6) and "3DAR1" (7) are
# absent: they reuse these columns as AR1 correlations.
#
# `what` names the term for the error message, `forms` the forms the slot is
# charged on, `dev_forms` those charged on the time-varying DEVIATES -- under
# Time_varying_sel = "Off" build_map() maps those deviates away, the term is
# identically zero, and a negative weight there is inert rather than wrong, so
# it is allowed (measured in test-selectivity-penalty-sd.R).
#
# The rule throughout: refuse a negative weight only where THIS fleet's wiring
# reaches it. The template gates the penalty block on
# `flt_type(flt) > 0 && flt_sel_lead(flt) == 1`, so an "Off" fleet and one
# following another's Selectivity_index are skipped too; .rce_sel_pen_lead()
# resolves the lead. The four cases a negative weight is NOT refused in are
# tabulated in vignette("model-parameterizations").
.RCE_SEL_PEN_POSITIVE <- list(
  # "Non-parametric" is the legacy spelling of 2 that switch_check() still
  # accepts; it is not in sel_map, so it has to be listed alongside the canonical
  # name or a workbook using it would skip the check.
  Sel_curve_pen1 = list(sd_col = "Sel_shape_sd",
                        # One slot weights a different term per form, so the
                        # message is grouped by form and names that form's term.
                        what   = c(default    = "the shape/smoothness penalty",
                                   LogisticPM = "the random walk on realized log-selectivity"),
                        forms  = c("NonParametric", "Non-parametric", "NonParametricPM",
                                   "LogisticPM", "NonParametricIntegrable"),
                        dev_forms = "LogisticPM"),
  Sel_curve_pen2 = list(sd_col = "Sel_curvature_sd",
                        what   = c(default = "the 2nd-difference curvature penalty"),
                        forms  = c("NonParametric", "Non-parametric", "NonParametricPM",
                                   "NonParametricIntegrable"),
                        dev_forms = character(0)),
  Sel_curve_pen3 = list(sd_col = "Sel_devmag_sd",
                        what   = c(default    = "the deviation-magnitude penalty",
                                   LogisticPM = "the random walk on the age-1 deviates"),
                        forms  = c("NonParametricPM", "LogisticPM"),
                        dev_forms = c("NonParametricPM", "LogisticPM"))
)

# switch_check() canonicalizes only at its last statement (revert_switches()), so
# the checks below see whatever the workbook holds -- every bundled data set
# stores the integer code. .canon_switch() further down this file reads a value
# through its own map and is what resolves them; there is no second list of
# spellings.

#' Whether a fleet's Sel_curve_pen1 may legitimately be negative
#'
#' "NonParametricPM" (9) is the one form whose shape term branches on the sign
#' of slot 1, and only under Sel_shape_mode = "Directional". Its "Smooth" mode
#' applies the weight two-sided, so a negative one rewards the shape there as it
#' would anywhere else. An unset Sel_shape_mode is "Directional", as in
#' rearrange_data().
#'
#' @param fleet_control the fleet control table.
#' @param flt the fleet's row index.
#' @return TRUE when a negative Sel_curve_pen1 is meaningful on that fleet.
#' @keywords internal
#' @noRd
.rce_sel_pen1_sign_aware <- function(fleet_control, flt) {
  if (!identical(.canon_switch(fleet_control$Selectivity[flt], sel_map),
                 "NonParametricPM")) return(FALSE)
  mode <- if (is.null(fleet_control$Sel_shape_mode)) NA else
    as.character(fleet_control$Sel_shape_mode[flt])
  !(!is.na(mode) && mode %in% c("1", "Smooth", "smooth"))
}

#' Which fleets carry their selectivity group's penalty, as the template counts it
#'
#' `ceattle.cpp` gates the penalty block on `flt_sel_lead`, which
#' `rearrange_data()` groups by `Selectivity_index` AND the selectivity form: two
#' fleets sharing an index but not a form are two groups, and each is charged.
#' The map's donor rule (`adjust_map_shared_params()`, and `.shared_block_lead()`
#' with it) groups by the index alone, so it is the wrong rule to predict the
#' penalty with -- it would call the second fleet a follower and skip a weight
#' the template does read.
#'
#' @param fleet_control the fleet control table.
#' @return a logical vector, TRUE where the template charges that fleet's penalty.
#' @keywords internal
#' @noRd
.rce_sel_pen_lead <- function(fleet_control) {
  n <- nrow(fleet_control)
  # Without the index column there is no group to share, so every fleet leads.
  if (is.null(fleet_control$Selectivity_index)) return(rep(TRUE, n))
  off  <- vapply(seq_len(n), function(i)
    identical(.canon_switch(fleet_control$Fleet_type[i], fleet_map), "Off"), logical(1))
  # The canonical name stands in for the integer code rearrange_data() pastes,
  # which is not quite one-to-one: .canon_switch() trims and .pull_int() does
  # not, so " NonParametric" resolves here but reaches the template as NA, a
  # different group there. Anything that does not resolve the template's way
  # gets a key of its own, so it leads here too and its weight is checked.
  # data_check() does reach this with an invalid code, since it accumulates
  # errors before refusing one, but every out-of-range integer canonicalizes to
  # "<blank>", which is in no .RCE_SEL_PEN_POSITIVE form list, so the caller's
  # form test skips those fleets and the lead never decides anything for them.
  raw   <- as.character(fleet_control$Selectivity)
  clean <- raw %in% names(sel_map) | !is.na(suppressWarnings(as.integer(raw)))
  form  <- vapply(seq_len(n), function(i)
    .canon_switch(fleet_control$Selectivity[i], sel_map), character(1))
  form[!clean] <- paste0("<unresolved ", which(!clean), ">")
  .group_lead(paste(fleet_control$Selectivity_index, form), off) == 1L
}

#' Report Sel_curve_pen weights that are negative where the template reads them
#'
#' One rule, two callers: `data_check()` passes the `fleet_control` columns and
#' `fit_mod()` the `sel_curve_pen` parameter actually in use. The second matters
#' because `sel_curve_pen` is a PARAMETER, so `inits` from a stored fit override
#' the columns and would otherwise carry a negative weight past the column check.
#'
#' @param fleet_control the fleet control table.
#' @param pen optional `[nrow(fleet_control), 3]` matrix of the weights in use;
#'   the `Sel_curve_pen1/2/3` columns are read when it is NULL.
#' @param source one clause naming where the offending value came from.
#' @return a character vector of error messages, one per penalty slot, empty when
#'   every weight is fine.
#' @keywords internal
#' @noRd
.rce_sel_pen_sign_errors <- function(fleet_control, pen = NULL, source = NULL) {
  errs <- character(0)
  # A `pen` that is not the expected matrix cannot be checked, and returning
  # clean would pass a negative weight through. Refuse rather than fail open.
  # Checked before the lead is resolved, so a malformed `pen` still gets this
  # message rather than one about a fleet_control column.
  if (!is.null(pen) && !identical(dim(pen), c(nrow(fleet_control), 3L))) stop(
    "`inits$sel_curve_pen` must be a ", nrow(fleet_control), " x 3 matrix, one ",
    "row per fleet and one column per Sel_curve_pen slot.", call. = FALSE)
  pen_lead <- .rce_sel_pen_lead(fleet_control)
  for (col in names(.RCE_SEL_PEN_POSITIVE)) {
    # Slot read by name, so reordering the registry cannot silently shift it.
    slot <- as.integer(sub("^Sel_curve_pen", "", col))
    bad  <- integer(0)
    for (flt in seq_len(nrow(fleet_control))) {
      # The template charges no selectivity penalty for an "Off" fleet.
      if (identical(.canon_switch(fleet_control$Fleet_type[flt], fleet_map), "Off")) next
      form <- .canon_switch(fleet_control$Selectivity[flt], sel_map)
      val  <- if (is.null(pen)) suppressWarnings(as.numeric(as.character(fleet_control[[col]][flt]))) else
        pen[flt, slot]
      if (length(val) != 1L || is.na(val) || val >= 0) next
      if (!form %in% .RCE_SEL_PEN_POSITIVE[[col]]$forms) next
      if (col == "Sel_curve_pen1" && .rce_sel_pen1_sign_aware(fleet_control, flt)) next
      # Slots charged on the deviates are never read when there are none.
      if (form %in% .RCE_SEL_PEN_POSITIVE[[col]]$dev_forms &&
          identical(.canon_switch(fleet_control$Time_varying_sel[flt], tv_sel_map), "Off")) next
      # The group's penalty is charged once, on the lead, so a follower's weight
      # is never read (ceattle.cpp gates on flt_sel_lead == 1 as well as flt_type).
      if (!pen_lead[flt]) next
      bad <- c(bad, flt)
    }
    if (!length(bad)) next
    # Grouped by form as well as slot: the same slot weights a different term on
    # each form, so one message spanning forms could not name the term.
    form_bad <- vapply(bad, function(f) .canon_switch(fleet_control$Selectivity[f], sel_map),
                       character(1))
    for (frm in unique(form_bad)) {
      b    <- bad[form_bad == frm]
      what <- .RCE_SEL_PEN_POSITIVE[[col]]$what
      errs <- c(errs, sprintf(
        paste0(
          "%s is negative on fleet(s) %s%s. It weights %s on %s, and a weight ",
          "multiplies a squared deviation, so a negative one rewards the deviation ",
          "it is meant to penalize. Use a positive weight, or set '%s' instead, ",
          "which is a standard deviation and cannot go negative."),
        col,
        paste(sprintf("'%s' (%s)", fleet_control$Fleet_name[b],
                      format(if (is.null(pen)) suppressWarnings(as.numeric(as.character(fleet_control[[col]][b])))
                             else pen[b, slot], trim = TRUE)),
              collapse = ", "),
        if (is.null(source)) "" else paste0(" ", source),
        if (frm %in% names(what)) what[[frm]] else what[["default"]],
        frm,
        .RCE_SEL_PEN_POSITIVE[[col]]$sd_col))
    }
  }
  errs
}

# Whether selectivity normalization pools its reference across sexes. Orthogonal
# to WHERE the reference is taken (Sel_norm_bin: a named bin, or the max), so the
# two combine rather than multiply into more columns.
sel_norm_scope_map <- c(
  "WithinSex"   = 0,  # each sex divided by its own reference; both reach 1
  "AcrossSexes" = 1   # one pooled reference; relative sex selectivity retained
)

# "AR1" (2) was removed in 5.16.0 and is refused by data_check(). It stays in
# the map, as q_map keeps its removed "AR1" (6), so a workbook carrying the
# integer 2 still canonicalizes and the refusal can name the fleet and the
# replacement instead of a bare "invalid value". It was never an AR1 -- value 2
# is scored with value 1's independent normal penalty, and there is no
# correlation parameter -- so the replacement is a linkage, ar1(1 | Year).
tv_sel_map <-c(
  "Off" = 0,
  "IID" = 1,
  "AR1" = 2,
  "Block" = 3,
  "RandomWalk" = 4,
  "RandomWalkAscending" = 5
)

q_map <- c(
  "Fixed" = 0,
  "Estimated" = 1,
  "Estimated-with-prior" = 2,
  "Analytical" = 3,
  "PowerEquation" = 4,
  "Environmental" = 5,
  "AR1" = 6,
  "AnalyticalArith" = 7
)

# "AR1" (2) is REMOVED and refused by data_check(), for the same reason and on
# the same terms as tv_sel_map's -- the model gives value 2 the identical
# independent normal penalty as value 1 (`index_varying_q == 1 || == 2`), and
# the QAR1 catchability path that held a correlation parameter is gone (5.37.0).
# An AR1 on catchability is a q linkage -- ar1(1 | Year).
tv_q_map <- c(
  "Off" = 0,
  "IID" = 1,
  "AR1" = 2,
  "Block" = 3,
  "RandomWalk" = 4
)

# Whether a fleet's selectivity is indexed by age or by length bin. Read in
# rearrange_data(), build_params() and data_check(); a value outside this set
# silently became NA and entered the fit as a missing dimension.
sel_dimension_map <- c(
  "Age"    = 0,
  "Length" = 1
)

# Direction of the (directional-mode) NonParametric shape penalty. The integer
# side is the ADMB sign-of-the-weight convention that .is_incr reads below:
# -1 means increasing, because the penalty weight is negated. Writing 1 here
# would be read as "not increasing" and apply the penalty in the opposite
# direction, so the map must carry the convention the code implements.
sel_shape_dir_map <- c(
  "Decreasing" =  0,
  "Increasing" = -1
)

# Shape-penalty mode: one-sided directional (ADMB/AMAK) or two-sided smooth.
sel_shape_mode_map <- c(
  "Directional" = 0,
  "Smooth"      = 1
)


comp_loglike_map <- c(
  "MultinomialAFSC" = -1,
  "Multinomial" = 0,
  "DirichletMultinomial" = 1
)

# Diet (stomach-content) composition likelihood family, per predator species. The
# diet likelihood supports the multinomial and Dirichlet-multinomial only -- there
# is no AFSC variant -- so this is comp_loglike_map without the -1 case.
diet_loglike_map <- c(
  "Multinomial" = 0,
  "DirichletMultinomial" = 1
)

# Survey/index biomass observation likelihood family (per fleet). Lognormal is
# on log(obs); the other four are on the natural-scale residual (obs - q*pred),
# where Log_sd is read as an ABSOLUTE sd rather than a log-scale CV. MVN and
# MVNORM differ only by the normalizing constant, so they fit identically and
# MVN is the one that matches ADMB. The allowed values, their formulas and
# when to prefer each are in vignette("model-options-and-functionality").
index_distribution_map <- c(
  "Lognormal" = 0,
  "MVN" = 1,
  "MVNORM" = 2,
  "Normal" = 3,
  "TruncatedNormal" = 4
)

#' Read a switch column by canonical name, whichever spelling it arrived in
#'
#' `switch_check()` -> `revert_switches()` upgrades integer codes to names, and
#' every path through `fit_mod()` runs it. `read_data()` does not, and
#' `data_check()` is callable on a list straight from it, so its checks
#' canonicalize first rather than test one spelling and silently pass the other.
#'
#' @param x Column values, names or integer codes.
#' @param map Named integer vector for this switch (e.g. `q_map`).
#' @return Canonical names; `"<blank>"` for `NA` or an out-of-map integer, and
#'   the trimmed value itself for anything else it cannot resolve.
#' @noRd
.canon_switch <- function(x, map) {
  x <- trimws(as.character(x))
  num <- suppressWarnings(as.integer(x))
  out <- ifelse(!is.na(num), names(map)[match(num, map)], x)
  out[is.na(out)] <- "<blank>"
  out
}


#' Which index rows are fitted on the NATURAL scale?
#'
#' `Index_distribution` splits into two scales, and almost everything downstream
#' of the likelihood assumed the lognormal one. `Lognormal` (0) is a log-scale
#' family whose `Log_sd` is a CV / log-sd; `MVN` (1), `MVNORM` (2), `Normal` (3)
#' and `TruncatedNormal` (4) are natural-scale families whose sd is ABSOLUTE, in
#' the units of the index. Applying a log-scale formula to the second group does
#' not error, it silently returns nonsense, because `sigma^2 / 2` is then a
#' number the size of the index squared.
#'
#' A new natural-scale family has to be added to the vector below as well as to
#' `index_distribution_map`, or every fleet using it silently reverts to the
#' log-scale treatment this function exists to prevent.
#'
#' @param data_list A `data_list` holding `fleet_control` and `index_data`.
#' @return Logical, one per `index_data` row; `FALSE` where the fleet is
#'   lognormal or cannot be resolved.
#' @keywords internal
#' @noRd
.index_rows_natural_scale <- function(data_list) {
  idx <- data_list$index_data
  fc  <- data_list$fleet_control
  if (is.null(idx) || !nrow(idx) || is.null(fc)) return(logical(0))
  nat <- .index_fleets_natural_scale(fc)
  if (!length(nat)) return(rep(FALSE, nrow(idx)))
  out <- nat[match(idx[["Fleet_code"]], fc[["Fleet_code"]])]
  out[is.na(out)] <- FALSE
  out
}


# The per-fleet form of .index_rows_natural_scale(), which calls it; the one
# place the family codes are classified, so a new family is added here only.
.index_fleets_natural_scale <- function(fleet_control) {
  fam <- fleet_control[["Index_distribution"]]
  if (is.null(fam)) return(rep(FALSE, nrow(fleet_control)))
  chr <- trimws(as.character(fam))
  num <- suppressWarnings(as.integer(chr))
  code <- ifelse(!is.na(num), num, as.integer(index_distribution_map[chr]))
  code[is.na(code)] <- 0L
  code %in% c(1L, 2L, 3L, 4L)   # MVN, MVNORM, Normal, TruncatedNormal
}


#' Fleet codes that hold survey-index observations the model fits
#'
#' An index is a property of the data, not of the fleet type: the model scores
#' an `index_data` row for any fleet that is not `Off`, so a fishery with a CPUE
#' series is fitted like a survey, on that fleet's own selectivity. Selecting on
#' `Fleet_type == "Survey"` instead would leave such a fleet with its
#' catchability frozen and its index absent from `plot_index()`.
#'
#' @param data_list A `data_list` holding `fleet_control` and `index_data`.
#' @param fitted_only Keep only rows the likelihood uses (positive `Year`, at or
#'   before `endyr`). A prediction-only row is not an observation and must not,
#'   for instance, make catchability estimable.
#' @return Integer vector of `Fleet_code`s, possibly empty.
#' @keywords internal
#' @noRd
.fleets_with_index <- function(data_list, fitted_only = TRUE) {
  idx <- data_list$index_data
  fc  <- data_list$fleet_control
  # nrow() is NULL for anything without dimensions -- a data_list carrying a
  # scalar NA where a table is expected reaches here, and `!NULL` errors.
  if (is.null(idx) || is.null(nrow(idx)) || !nrow(idx) || is.null(fc)) {
    return(integer(0))
  }
  keep <- rep(TRUE, nrow(idx))
  if (fitted_only) {
    keep <- idx$Year > 0
    if (!is.null(data_list$endyr)) keep <- keep & idx$Year <= data_list$endyr
    if (!is.null(idx$Observation)) keep <- keep & !is.na(idx$Observation)
  }
  live <- fc$Fleet_code[!(fc$Fleet_type %in% c("Off", 0))]
  as.integer(sort(intersect(unique(idx$Fleet_code[keep]), live)))
}


fleet_map <- c(
  "Fishery" = 1,
  "Survey" = 2,
  "Off" = 0
)

# How the initial age structure is built; 2 is the default. Modes 1 and 5 both
# start from R_init, and 5 displaces it by the year-1 recruitment deviation,
# the only term separating them (init_log_scalar in ceattle.cpp; Cole Monnahan
# / AFSC GOA pollock convention). Modes 3, 4 and 6 estimate Finit, and only 6
# decays with sum(M1 + Finit * sel(a)), Stock Synthesis's InitF convention and
# the one equilibrium under a size-selective fishery, where Finit is the apical
# initial F. Every mode is tabulated in
# vignette("model-options-and-functionality").
initMode_map <- c(
  "FreeParams"                 = 0,
  "Equilibrium"                = 1,
  "NonEquilibrium"             = 2,
  "FishedNonEquilibrium"       = 3,
  "FishedNonEquilibriumScaled" = 4,
  "OffsetEquilibrium"          = 5,
  "FishedNonEquilibriumSelected" = 6
)

# Predator-prey suitability mode (per predator species)
# 1 and 3 are blocked in data_check() until growth-model validation is added
suitMode_map <- c(
  "Empirical"       = 0,
  "GammaLength"     = 1,   # NOT YET AVAILABLE
  "GammaWeight"     = 2,
  "LognormalLength" = 3,   # NOT YET AVAILABLE
  "LognormalWeight" = 4,
  "NormalLength"    = 5,
  "NormalWeight"    = 6
)

# Harvest control rule
hcr_map <- c(
  "NoFishing"    = 0,
  "CMSY"         = 1,
  "ConstantF"    = 2,
  "ConstantFSSB" = 3,
  "ConstantFSPR" = 4,
  "NPFMC"        = 5,
  "PFMC"         = 6,
  "SESSF"        = 7
)

# Numbers-at-age estimation mode, per species (data_list$estDynamics).
# 0 = estimate the population dynamics; 1 = use fixed input numbers-at-age from
# NByageFixed; 2 = scale the fixed input by one estimated coefficient, under
# predation only. Named per the Fixed/Estimated convention so "0 = not
# estimated" reads plainly. Code 3 (an age-specific scalar) was retired in
# 5.35.0; it always fitted as 2.
estDynamics_map <- c(
  "Estimated"   = 0,
  "Fixed"       = 1,
  "FixedScaled" = 2
)

# Validate a switch value against its map WITHOUT converting it, so a typo is
# reported where it was written rather than several functions later. It accepts
# every legal form -- canonical string, integer code, either as a per-species
# vector, NULL/NA for an unset switch -- because fit_mod() rebuilds a
# model_config() from resolved values once the fit is done, and anything
# stricter would throw there and discard the finished fit.
.check_switch <- function(x, map, name) {
  if (is.null(x) || !length(x)) return(invisible(x))
  vals <- as.character(x[!is.na(x)])
  if (!length(vals)) return(invisible(x))
  # Treat 2, 2L and "2" alike: normalise anything numeric-looking to its
  # canonical numeric string before matching.
  num <- suppressWarnings(as.numeric(vals))
  vals[!is.na(num)] <- as.character(num[!is.na(num)])
  bad <- unique(vals[!vals %in% c(names(map), as.character(as.numeric(map)))])
  if (length(bad)) {
    stop(sprintf(
      "Invalid `%s` value(s): %s.\nUse one of: %s (or the integer codes %s).",
      name, paste(bad, collapse = ", "),
      paste(names(map), collapse = ", "),
      paste(unname(map), collapse = ", ")), call. = FALSE)
  }
  invisible(x)
}

# Apply a string->integer switch map to a (possibly character) switch value,
# passing integers through unchanged and erroring clearly on an unknown string.
# Used for switches consumed as bare integers before convert_switches() runs
# (e.g. estDynamics is read numerically by build_map()), so the string must be
# resolved early, in switch_check().
.map_switch <- function(x, map, name) {
  if (is.null(x)) return(x)
  # A factor (read.csv with stringsAsFactors = TRUE) would otherwise pass through
  # and be read as its level index; a numeric-looking string is the code it names.
  if (is.factor(x)) x <- as.character(x)
  if (!is.character(x)) return(x)
  num    <- suppressWarnings(as.numeric(x))
  is_num <- !is.na(num)
  bad <- c(setdiff(x[!is.na(x) & !is_num], names(map)),
           x[is_num & !num %in% unname(map)])
  if (length(bad) > 0) {
    stop(sprintf("Invalid '%s' value(s): %s. Options: %s (or the integer codes %s).",
                 name, paste(unique(bad), collapse = ", "),
                 paste(names(map), collapse = ", "),
                 paste(unname(map), collapse = ", ")), call. = FALSE)
  }
  out <- unname(map[x])   # map[NA] is NA, so off-fleet NAs pass through
  out[is_num] <- num[is_num]
  out
}

# `est_M1` was renamed to `M1_model`. Fold the deprecated name into `M1_model`
# so it is honoured everywhere `M1_model` is (fit_mod's M1Fun reconciliation,
# switch_check's default, combine_data). Only fires when `M1_model` is unset and
# `est_M1` carries a real (non-all-NA) value -- an all-NA `est_M1` (e.g. the
# GOA2018SS placeholder) means "not specified" and must fall through to the
# `M1_model` default, not pin it to NA.
.alias_est_M1 <- function(data_list) {
  if (is.null(data_list$M1_model) && !is.null(data_list$est_M1) &&
      !all(is.na(data_list$est_M1))) {
    data_list$M1_model <- data_list$est_M1
    message("'est_M1' is deprecated; use 'M1_model'.")
  }
  data_list
}

# Observation-SD estimation mode for a survey index or catch series
# (fleet_control$Estimate_index_sd / Estimate_catch_sd). 0 = use the fixed SD
# implied by the data CV (not estimated); 1 = estimate; 2 = analytical
# (Walters-Ludwig 1994). "Fixed" reads plainly for the 0 = not-estimated case.
estimate_sd_map <- c(
  "Fixed"      = 0,
  "Estimated"  = 1,
  "Analytical" = 2
)

# Stock-recruit hyperparameter estimation mode (data_list$srr_est_mode; set via
# build_srr()). 0 = fix alpha/steepness at its prior mean (not estimated);
# 1 = freely estimate; 2 = lognormal prior; 3 = beta prior.
srr_est_mode_map <- c(
  "Fixed"          = 0,
  "Estimated"      = 1,
  "LognormalPrior" = 2,
  "BetaPrior"      = 3
)

# What fit_mod() does with the model (the `estimateMode` argument).
# 0 = fit the hindcast and the HCR projection; 1 = fit the hindcast only;
# 2 = project only, from the supplied inits; 3 = build MakeADFun without
# optimizing (returns the real objective, for pre-fit diagnostics); 4 = optimize
# with every hindcast parameter mapped out (a plumbing smoke test).
estimateMode_map <- c(
  "Estimate"      = 0,
  "Hindcast"      = 1,
  "Projection"    = 2,
  "DebugBuild"    = 3,
  "DebugOptimize" = 4
)

# Predation-mortality mode (the `msmMode` argument). 0 = single-species (no
# predation); 1 = Type II MSVPA (Holsman et al. 2015); 2 = Type III MSVPA. Higher
# integer codes (Kinzey-Punt, Holling forms) are not yet implemented.
msmMode_map <- c(
  "SingleSpecies" = 0,
  "MSVPA"         = 1,
  "TypeIIIMSVPA"  = 2
)


# The words Sel_norm_bin and Sel_norm_bin_upper accept, matched in any case.
.RCE_SEL_NORM_WORDS <- list(
  off = c("off", "none", "no", "na", "nan", ""),   # do not normalize
  max = c("max", "maximum"),                       # normalize by the largest value
  all = c("all")                                   # LogisticPM: the whole penalty range
)

# Refuse a Sel_norm_bin the model cannot use, naming the fleet. Shared so
# rearrange_data() fails the same way switch_check() does: a value past the last
# selected bin indexes off the end of the selectivity array, and the model builds
# safebounds = FALSE.
#' @noRd
.rce_stop_bad_sel_norm <- function(code, col, fleet_control, hi = NULL) {
  bad <- attr(code, "unrecognised")
  if (length(bad)) {
    stop(sprintf("Fleet '%s': '%s' is '%s', which is not a bin or one of: Max, Off, All.",
                 fleet_control$Fleet_name[bad[1]], col,
                 as.character(fleet_control[[col]][bad[1]])), call. = FALSE)
  }
  if (!is.null(hi)) {
    over <- which(!is.na(code) & !is.na(hi) & code > 0 & code > hi)
    if (length(over)) {
      stop(sprintf("Fleet '%s': '%s' is %s, past the last bin it is selected over (%s).",
                   fleet_control$Fleet_name[over[1]], col, code[over[1]],
                   hi[over[1]]), call. = FALSE)
    }
  }
  invisible()
}


#' First and last bin a fleet may normalize at, on the column's own scale
#'
#' Below `Bin_first_selected` the curve is zeroed, so a reference taken there
#' divides by nothing. `Bin_first_selected` is a 1-based bin ordinal while
#' `Sel_norm_bin` is an absolute age, so on an age-based fleet the first
#' selected age is `minage + Bin_first_selected - 1`; on a length-based fleet
#' both are bin ordinals and no conversion applies.
#'
#' @param data_list An Rceattle data list.
#' @param flt Row of `fleet_control`.
#' @param age_selex Whether the fleet is age-based.
#' @param max_bin The fleet's bin count (`nages` or `nlengths`).
#' @return A list of `lo` and `hi`; either may be `NA` if the dimension is unset.
#' @keywords internal
#' @noRd
.rce_sel_norm_bounds <- function(data_list, flt, age_selex, max_bin) {
  first <- suppressWarnings(as.numeric(
    data_list$fleet_control$Bin_first_selected[flt]))
  if (length(first) == 0 || is.na(first) || first < 1) first <- 1
  if (isTRUE(age_selex)) {
    minage <- suppressWarnings(as.numeric(
      data_list$minage[data_list$fleet_control$Species[flt]]))
    if (length(minage) == 0) minage <- NA_real_
    list(lo = minage + first - 1, hi = minage + max_bin - 1)
  } else {
    list(lo = first, hi = max_bin)
  }
}


#' Where a fleet normalizes its selectivity, as a number
#'
#' The one place the column is read, since it holds either a word or a bin.
#' Blank means do not normalize, and "Max" the largest value. A value below the
#' fleet's first bin also means the largest value: that is what a negative has
#' always meant, so `lo` must be the fleet's own first bin: the species' minage
#' on an age-based fleet, 1 on a length-based one. On a stock recruiting at age 0
#' a `Sel_norm_bin` of 0 is the first age, not a flag.
#'
#' @param x The raw column.
#' @param lo The fleet's first bin; anything below it means the maximum.
#' @param allow_all Whether "All" is accepted, i.e. a LogisticPM fleet.
#' @return Numeric: `NA` to not normalize, `-1` for the maximum, else the bin.
#'   An unrecognised word is `NA`, its position in the "unrecognised" attribute.
#' @keywords internal
#' @noRd
.rce_sel_norm_code <- function(x, lo = 1L, allow_all = FALSE) {
  if (is.null(x)) return(NULL)
  chr <- trimws(tolower(as.character(x)))
  out <- suppressWarnings(as.numeric(chr))

  words <- .RCE_SEL_NORM_WORDS
  if (!isTRUE(any(allow_all))) words$all <- character(0)
  is_off <- chr %in% words$off | is.na(chr)
  is_max <- chr %in% c(words$max, words$all)

  out[is_off] <- NA_real_
  out[is_max] <- -1
  lo <- suppressWarnings(as.numeric(lo))
  lo[is.na(lo)] <- 1                                  # unknown dimension: assume 1
  below <- !is.na(out) & !is_max & out < lo
  out[below] <- -1

  bad <- is.na(out) & !is_off
  attr(out, "unrecognised") <- which(bad)
  out
}


#' Resolve factor columns in a fleet_control to their labels
#'
#' `read.csv(stringsAsFactors = TRUE)` factors every column. A factor switch is
#' then read by its LEVEL INDEX downstream -- `convert_switches()` passes it
#' through and `as.integer()`s it -- so `Fleet_type = factor(c("Off","Survey"))`
#' reaches the template as 1, 2: the Off fleet fits as a fishery. The label is
#' what the user wrote, so resolve it before anything reads the column.
#'
#' @param fleet_control the fleet control table.
#' @return The table with every factor column as character.
#' @keywords internal
#' @noRd
.rce_defactor_fleet_control <- function(fleet_control) {
  if (is.null(fleet_control)) return(fleet_control)
  fac <- vapply(fleet_control, is.factor, logical(1))
  if (any(fac)) fleet_control[fac] <- lapply(fleet_control[fac], as.character)
  fleet_control
}

#' Name fleets for a message: Fleet_name, or "row N" where it is blank too
#'
#' `Fleet_name` has no schema default either, so a message that interpolates it
#' can name nothing at all.
#'
#' @param fleet_control the fleet control table.
#' @param i row indices to name.
#' @return A character vector, one entry per element of `i`.
#' @keywords internal
#' @noRd
.rce_fleet_who <- function(fleet_control, i) {
  who <- paste("row", i)
  if (!is.null(fleet_control$Fleet_name)) {
    nm <- as.character(fleet_control$Fleet_name)[i]
    who[!is.na(nm) & trimws(nm) != ""] <- nm[!is.na(nm) & trimws(nm) != ""]
  }
  who
}

#' Refuse a blank Fleet_type, naming the fleets
#'
#' `Fleet_type` has no schema default, so a blank one is not "unset, take the
#' default" -- it is a fleet whose role in the likelihood nobody stated, and the
#' package's two halves read it differently. Called from `switch_check()` and
#' again from `rearrange_data()`, which is exported and reachable without it.
#'
#' @param fleet_control the fleet control table.
#' @return `invisible(NULL)`; stops when any Fleet_type is blank.
#' @keywords internal
#' @noRd
.rce_stop_blank_fleet_type <- function(fleet_control) {
  if (is.null(fleet_control) || is.null(fleet_control$Fleet_type)) return(invisible(NULL))
  bad <- which(is.na(fleet_control$Fleet_type) |
                 trimws(as.character(fleet_control$Fleet_type)) == "")
  if (!length(bad)) return(invisible(NULL))
  who <- .rce_fleet_who(fleet_control, bad)
  stop("'Fleet_type' is blank for fleet(s) ", paste(who, collapse = ", "),
       ". It has no default: set it to one of ", paste(names(fleet_map), collapse = ", "),
       " (or the integer codes ", paste(fleet_map, collapse = ", "),
       "). Use \"Off\" for a fleet the model should carry but not fit.", call. = FALSE)
}

#' Refuse a Fleet_type the package cannot resolve, naming the fleet and the value
#'
#' `convert_switches()` passes an unrecognized value through and then
#' `as.integer()`s it, so `"Fisherie"` reaches the template as `NA` and `-1` or
#' `3` reach it as themselves -- a fleet whose role in the likelihood nobody
#' stated either way. The allowed set is the one `validate_switches()` applies on
#' the `data_check()` path; this is for `rearrange_data()`, which is exported and
#' has no such pass. A blank is the other guard's to report.
#'
#' @param fleet_control the fleet control table.
#' @return `invisible(NULL)`; stops when a Fleet_type does not resolve.
#' @keywords internal
#' @noRd
.rce_stop_unreadable_fleet_type <- function(fleet_control) {
  ft <- fleet_control$Fleet_type
  if (is.null(ft)) return(invisible(NULL))
  blank <- is.na(ft) | trimws(as.character(ft)) == ""
  bad   <- which(!blank & !ft %in% c(fleet_map, names(fleet_map)))
  if (!length(bad)) return(invisible(NULL))
  stop("'Fleet_type' could not be read for fleet(s) ",
       paste(sprintf("%s ('%s')", .rce_fleet_who(fleet_control, bad),
                     as.character(ft)[bad]), collapse = ", "),
       ". Use one of ", paste(names(fleet_map), collapse = ", "),
       " or the integer codes ", paste(fleet_map, collapse = ", "), ".",
       call. = FALSE)
}

#' Function to check for missing switches for map and parameter functions
#'
#' @param data_list Rceattle data list
#'
#' @export
#'

switch_check <- function(data_list){

  # Helper to set defaults and notify. Pass msg = NULL to fill the default
  # silently (for low-level switches that should not add console noise).
  set_default <- function(val, default, msg = NULL) {
    if(is.null(val)) {
      if(!is.null(msg)) message(msg)
      return(default)
    }
    return(val)
  }

  # Upgrade deprecated fleet_control column names to their canonical spellings,
  # from the schema's `aliases` field -- the one list of accepted old names.
  # First in switch_check() so the rename lands before build_params() reads the
  # columns and before the non-parametric penalty migration below.
  data_list$fleet_control <-
    .rce_upgrade_fleet_control_aliases(data_list$fleet_control)
  # Before anything reads or writes a switch: assigning "Off" into a factor whose
  # levels lack it writes NA, and every switch column is read by level index.
  data_list$fleet_control <-
    .rce_defactor_fleet_control(data_list$fleet_control)

  # Whether a fleet is fit at all is not a thing to infer, and the two halves of
  # the package disagree about a blank: `Fleet_type != "Off"` is NA, while
  # build_map_selectivity() reads NA as estimated. Refuse before either runs.
  .rce_stop_blank_fleet_type(data_list$fleet_control)
  # ...and the deprecated control / bioenergetics element names (e.g.
  # `sigma_rec_prior` -> `sigma_rec`, `Diet_loglike` -> `Diet_distribution`).
  data_list <- .rce_upgrade_data_list_aliases(data_list)

  # `Accumulation_age_lower` / `Accumulation_age_upper` are removed. They were
  # only ever range-validated and never applied (no composition-age grouping
  # exists anywhere in the R pipeline or the C++ template), so the columns did
  # nothing. Drop them -- once, with a soft-deprecation message -- if an old
  # workbook / data list still carries either spelling, so a stale column does
  # not silently propagate.
  drop_deprecated_col <- function(fc, old) {
    if (!is.null(fc[[old]])) {
      fc[[old]] <- NULL
      message(sprintf(
        "'%s' is deprecated and ignored; it was never applied to composition data and has been removed.",
        old))
    }
    fc
  }
  data_list$fleet_control <- drop_deprecated_col(
    data_list$fleet_control, "Accumulation_age_lower")
  data_list$fleet_control <- drop_deprecated_col(
    data_list$fleet_control, "Accumulation_age_upper")

  # `growth_re` and `growth_indices` are removed. `growth_re` was documented as
  # the way to put random effects on growth, but nothing consumed it: the
  # deviation array was mapped off in every configuration and the model gave
  # it no density, so setting it changed no fit. Drop them with a message rather
  # than silently, so a data list still carrying `growth_re = 1` says where the
  # feature went instead of appearing to work.
  for (.old in c("growth_re", "growth_indices")) {
    if (!is.null(data_list[[.old]])) {
      data_list[[.old]] <- NULL
      message(sprintf(
        paste0("'%s' is deprecated and ignored; it never had an effect on any ",
               "fit. Specify time-varying growth with build_growth(linkages = ) ",
               "-- see vignette('environmental-linkages-and-priors')."),
        .old))
    }
  }

  if(is.null(data_list$srr_fun)){
    data_list$srr_fun <- 0
    data_list$srr_pred_fun <- 0
    data_list$srr_est_mode <- 0
    message("'srr_fun' are not included in data, assuming 0")
  }

  # Gate the "assuming <default>" messages so each fires only where the model
  # reads that input, rather than nagging about inputs the configuration never
  # consumes. Read before the defaults below overwrite the raw fleet_control
  # values: weight-length (alpha_wt_len / beta_wt_len) and the selectivity
  # dimension are read only when growth is estimated (growth_model > 0), the
  # CAAL defaults only with CAAL data, and the normalization upper bin only on a
  # fleet normalizing at a named bin (Sel_norm_bin >= 0), not max or off.
  .dflt_when <- list(
    growth_estimated = isTRUE(any(data_list$growth_model > 0)),
    has_caal         = isTRUE(nrow(data_list$caal_data) > 0),
    sel_norm_upper   = isTRUE(any(.rce_sel_norm_code(
      data_list$fleet_control$Sel_norm_bin) > 0, na.rm = TRUE)),
    # sel_norm_scope_flip: the one configuration the "AcrossSexes" default
    # changes -- a two-sex fleet normalizing at a named bin. Restricted to the
    # fleets selectivity.hpp runs its normalization block on: an "Off" fleet is
    # skipped, Hake normalizes in its own year/sex block, and LogisticPM reads
    # Sel_norm_bin1/2 as a penalty age range, so without the gate the message
    # cries wolf on any AMAK-style model that sets a penalty range.
    sel_norm_scope_flip = isTRUE(any(
      data_list$nsex[data_list$fleet_control$Species] == 2 &
        .rce_sel_norm_code(data_list$fleet_control$Sel_norm_bin) > 0 &
        !data_list$fleet_control$Fleet_type %in% c(0, "Off") &
        !data_list$fleet_control$Selectivity %in%
          c(5, "Hake", 11, "LogisticPM"), na.rm = TRUE)),
    multispecies     = isTRUE(any(data_list$msmMode > 0))
  )

  # Model and multi-species switches
  data_list$estDynamics <- set_default(data_list$estDynamics, rep(0, data_list$nspp), "'estDynamics' are not included in data, assuming 0")
  # Code 3 was retired in 5.35.0: it never estimated an age-specific multiplier
  # and fitted as 2, so a stored fit is rebuilt with 2 and is unchanged.
  .ed3 <- which(data_list$estDynamics %in% c(3, "3", "FixedScaledByAge"))
  if (length(.ed3)) {
    # estDynamics is per species, so name the ones to edit rather than the value.
    # Only when it really is per species: a scalar applies to every species, and
    # naming the first would send the user to the wrong cell.
    .per_sp <- length(data_list$estDynamics) == data_list$nspp
    .who <- if (!.per_sp) {
      "every species"
    } else if (!is.null(data_list$spnames) &&
               length(data_list$spnames) >= max(.ed3)) {
      paste(data_list$spnames[.ed3], collapse = ", ")
    } else {
      paste0("species ", paste(.ed3, collapse = ", "))
    }
    stop("estDynamics = 3 ('FixedScaledByAge') for ", .who,
         " was retired in 5.35.0: it never estimated an age-specific multiplier ",
         "and fitted as 2 ('FixedScaled'). Set estDynamics = 2 for ", .who,
         "; the fit is unchanged.", call. = FALSE)
  }
  # Resolve readable strings ("Fixed"/"Estimated"/...) to integer codes now --
  # build_map()/build_params() read estDynamics numerically, before
  # convert_switches() runs.
  data_list$estDynamics <- .map_switch(data_list$estDynamics, estDynamics_map, "estDynamics")
  .bad_ed <- setdiff(stats::na.omit(suppressWarnings(as.numeric(data_list$estDynamics))),
                     estDynamics_map)
  if (length(.bad_ed)) {
    stop(sprintf("Invalid 'estDynamics' value(s): %s. Options: %s (or the integer codes %s).",
                 paste(.bad_ed, collapse = ", "),
                 paste(names(estDynamics_map), collapse = ", "),
                 paste(unname(estDynamics_map), collapse = ", ")), call. = FALSE)
  }
  data_list$suitMode <- .map_switch(data_list$suitMode, suitMode_map, "suitMode")
  if (!is.null(data_list$fleet_control$Estimate_index_sd)) {
    data_list$fleet_control$Estimate_index_sd <-
      .map_switch(data_list$fleet_control$Estimate_index_sd, estimate_sd_map, "Estimate_index_sd")
  }
  if (!is.null(data_list$fleet_control$Estimate_catch_sd)) {
    data_list$fleet_control$Estimate_catch_sd <-
      .map_switch(data_list$fleet_control$Estimate_catch_sd, estimate_sd_map, "Estimate_catch_sd")
  }
  # Diet inputs are consumed only under predation; fill silently otherwise.
  data_list$Diet_comp_weights <- set_default(data_list$Diet_comp_weights, rep(1, data_list$nspp),
    if (.dflt_when$multispecies) "'Diet_comp_weights' are not included in data, assuming 1")
  data_list$Diet_distribution <- set_default(data_list$Diet_distribution, rep(0, data_list$nspp),
    if (.dflt_when$multispecies) "'Diet_distribution' are not included in data, assuming 'Multinomial'")
  # Resolve readable strings ("Multinomial" / "DirichletMultinomial") to the integer
  # codes the C++ diet likelihood reads; integer input passes through unchanged.
  data_list$Diet_distribution <- .map_switch(data_list$Diet_distribution, diet_loglike_map, "Diet_distribution")
  # alpha_wt_len / beta_wt_len (the length -> weight conversion) are only used when
  # growth is estimated, so only announce the default then; otherwise fill silently.
  data_list$alpha_wt_len <- set_default(data_list$alpha_wt_len, 1e-6,
    if (.dflt_when$growth_estimated) "'alpha_wt_len' not specified in data, assuming 1e-6")
  data_list$beta_wt_len <- set_default(data_list$beta_wt_len, 3,
    if (.dflt_when$growth_estimated) "'beta_wt_len' not specified in data, assuming 3")
  data_list <- .alias_est_M1(data_list)
  data_list$M1_model <- set_default(data_list$M1_model, rep(0, data_list$nspp), "'M1_model' is not included in data, assuming 0")
  data_list$msmMode <- set_default(data_list$msmMode, 0, "'msmMode' is not included in data, assuming single-species (0)")
  data_list$M1_re <- set_default(data_list$M1_re, rep(0, data_list$nspp), "'M1_re' is not in data, assuming 0 for all species")
  # Canonicalise to the integer codes build_map_m1() compares against: a string
  # matches no arm (`"iid_year" %in% c(2, 5)` is FALSE), which would map every
  # deviation out and leave the random effect silently absent. build_M1()
  # canonicalises its own arguments, but a list from build_data(),
  # combine_data(), the deprecated est_M1 alias or built by hand never passes
  # through it. .map_switch() rather than build_M1()'s stricter
  # .coerce_M1_arg(), because it also takes a factor, a numeric-looking string
  # and an NA (test-switches-map-switch-factor.R).
  data_list$M1_model <- .map_switch(data_list$M1_model, .M1_MODELS, "M1_model")
  data_list$M1_re    <- .map_switch(data_list$M1_re,    .M1_RES,    "M1_re")
  data_list$initMode <- set_default(data_list$initMode, 2, "'initMode' is not in the data, setting to 2 (default)")
  data_list$comp_offset <- set_default(data_list$comp_offset, 1e-5, NULL) # Composition proportion offset (added to comp/caal before the multinomial. Filled silently; fit_control(comp_offset=) can override at fit time.

  # Bioenergetics scalars: TMB declares them as DATA_VECTOR length-nspp, so they
  # must exist even where nothing reads them. Under msmMode == 0 no consumption
  # quantity enters the likelihood, so a missing entry is filled silently; under
  # msmMode > 0 they are left for data_check() to report as missing or
  # wrong-length. The consumption code still runs in single-species mode --
  # ceattle.cpp calls calculate_temperature_function() for every species on
  # every fit -- so the sentinel is Ceq 4 (constant, fT = 1), the one equation
  # reading no environmental index; Ceq 1 would send the template to
  # env_index(year, Cindex) with no covariate supplied. Cindex is 1-based, so 1
  # is its lowest valid value even though Ceq 4 never reads it.
  if(data_list$msmMode == 0){
    bioenergetics_defaults <- list(
      Ceq    = rep(4L,  data_list$nspp),
      Cindex = rep(1L,  data_list$nspp),
      Pvalue = rep(1,   data_list$nspp),
      fday   = rep(365, data_list$nspp),
      CA     = rep(0,   data_list$nspp),
      CB     = rep(0,   data_list$nspp),
      Qc     = rep(0,   data_list$nspp),
      Tco    = rep(0,   data_list$nspp),
      Tcm    = rep(0,   data_list$nspp),
      Tcl    = rep(0,   data_list$nspp),
      CK1    = rep(0,   data_list$nspp),
      CK4    = rep(0,   data_list$nspp)
    )
    for(nm in names(bioenergetics_defaults)){
      if(is.null(data_list[[nm]])){
        data_list[[nm]] <- bioenergetics_defaults[[nm]]
      }
    }
  }

  # 1. Fleet Control defaults ----
  # Per-column defaults + messages come from the canonical schema
  # (R/0-column_schema.R) via .rce_apply_default(), so they live in one place
  # instead of being hand-copied here. The order of the calls (and thus the
  # message order) is unchanged.
  .sch <- .rce_column_schema()
  # `[[`, not `$`, for the three columns with exactly one longer sibling
  # (Sel_norm_bin, Time_varying_sel, Time_varying_q): `$` partial-matches the
  # sibling, so .rce_apply_default() is handed a non-NULL value, returns early,
  # and the sibling's values stand in for the default. Every other name here is
  # unique, or ambiguous and so safe -- `$` gives NULL on two candidates.
  data_list$fleet_control$Sel_norm_bin <- .rce_apply_default(data_list$fleet_control[["Sel_norm_bin"]], "Sel_norm_bin", .sch)
  data_list$fleet_control$Sel_norm_bin_upper <- .rce_apply_default(data_list$fleet_control[["Sel_norm_bin_upper"]], "Sel_norm_bin_upper", .sch, conditions = .dflt_when)
  # Absent means "no time variation" for either switch. This has to run BEFORE
  # the pre-4.4 non-parametric format migration below, which reads
  # Time_varying_sel to decide whether a workbook is legacy.
  data_list$fleet_control$Time_varying_sel <- .rce_apply_default(data_list$fleet_control[["Time_varying_sel"]], "Time_varying_sel", .sch)
  # Time_varying_q doubles as the env_data column INDEX for a fleet whose
  # Catchability is "Environmental" or "AR1", where "Off" is not a meaningful
  # value -- so default the column only when no fleet needs it as an index, and
  # say so plainly otherwise. Without this the fill produced "Off" for such a
  # fleet and the model died later on `missing value where TRUE/FALSE needed`.
  if (is.null(data_list$fleet_control[["Time_varying_q"]])) {
    .q_env <- data_list$fleet_control$Catchability %in% c("AR1", "Environmental")
    if (any(.q_env, na.rm = TRUE)) {
      stop("`fleet_control` has no `Time_varying_q` column, but fleet(s) ",
           paste(data_list$fleet_control$Fleet_name[which(.q_env)],
                 collapse = ", "),
           " set Catchability to \"Environmental\" or \"AR1\", which read ",
           "Time_varying_q as the env_data column index rather than as a ",
           "switch. Add the column: an index for those fleets, \"Off\" for ",
           "the rest.", call. = FALSE)
    }
  }
  data_list$fleet_control$Time_varying_q <- .rce_apply_default(data_list$fleet_control[["Time_varying_q"]], "Time_varying_q", .sch)
  # Sel_norm_scope defaults per-cell, not just per-column: .rce_apply_default()
  # returns early once the column exists, but a blank cell is the same silent
  # behaviour flip as a missing column -- it would reach the TMB integer vector
  # as NA, which the C++ reads as WithinSex. Announce whichever case applies, so
  # the "your two-sex named-bin fleet now pools" notice cannot be skipped by
  # supplying the column and leaving one row empty.
  data_list$fleet_control$Sel_norm_scope <- .rce_apply_default(
    data_list$fleet_control$Sel_norm_scope, "Sel_norm_scope", .sch,
    conditions = .dflt_when)
  if (anyNA(data_list$fleet_control$Sel_norm_scope)) {
    data_list$fleet_control$Sel_norm_scope[is.na(data_list$fleet_control$Sel_norm_scope)] <-
      .sch[["Sel_norm_scope"]]$default
    if (isTRUE(.dflt_when$sel_norm_scope_flip))
      message(.sch[["Sel_norm_scope"]]$default_msg)
  }
  # Sel_curve_pen1/2/3 only matter for non-parametric (Ianelli/PM), Hake, or
  # LogisticPM (random-walk weights) selectivity; only warn about missing columns
  # when such a fleet is present, otherwise default silently (avoids noise for
  # logistic-only models).
  .np_hake <- any(data_list$fleet_control$Selectivity %in%
                    c(2, "NonParametric", "Non-parametric", 9, "NonParametricPM", 13, "NonParametricIntegrable", 5, "Hake", 11, "LogisticPM"))
  # Alternative to the cryptic penalty WEIGHTS: a standard deviation. Every such
  # penalty is a Gaussian SSQ `weight * x^2 = x^2 / (2*sd^2)`, so
  # `weight = 1/(2*sd^2)`. Converted here only where the fleet left the legacy
  # weight unset, so legacy models stay bit-identical. An SD on a form whose
  # slot is not a weight (a logit-scale 2D/3D-AR1 correlation, or unused) is
  # refused rather than converted, as is a non-positive or non-finite SD (sd = 0
  # gives an Inf penalty; a negative SD squares to a spurious weight). See
  # vignette("model-parameterizations"), "Selectivity penalty weights".
  # Whether the workbook SHIPPED the penalty columns, recorded before the
  # defaults below create them. This is what tells a pre-4.4 non-parametric
  # workbook from a modern one (see the format upgrade further down); a blank
  # cell in a column that exists does not.
  .had_sel_curve_pen <- "Sel_curve_pen1" %in% names(data_list$fleet_control)
  .fc  <- data_list$fleet_control
  .col <- function(nm) if (is.null(.fc[[nm]])) rep(NA_real_, nrow(.fc)) else suppressWarnings(as.numeric(.fc[[nm]]))
  .np  <- c(2, "NonParametric", "Non-parametric", 9, "NonParametricPM", 13, "NonParametricIntegrable")
  .lpm <- c(11, "LogisticPM")
  # SD column -> (target Sel_curve_pen slot, forms that use it as a weight).
  # Slot 1 (shape) is read by every non-parametric form and by LogisticPM; slot 2
  # (curvature) by the non-parametric forms only. Slot 3 (deviation magnitude) is
  # read by "NonParametricPM" (9) and "LogisticPM" (11) alone: 2 and 13 also
  # estimate selectivity deviates, but score them with a Gaussian density on
  # Time_varying_sel_sd, so Sel_devmag_sd is inert on both.
  .np_pen3 <- c(9, "NonParametricPM")
  .sd_specs <- list(
    Sel_shape_sd     = list(pen = "Sel_curve_pen1", forms = c(.np, .lpm),
                            label = "NonParametric/NonParametricPM/NonParametricIntegrable/LogisticPM"),
    Sel_curvature_sd = list(pen = "Sel_curve_pen2", forms = .np,
                            label = "NonParametric/NonParametricPM/NonParametricIntegrable"),
    Sel_devmag_sd    = list(pen = "Sel_curve_pen3", forms = c(.np_pen3, .lpm),
                            label = "NonParametricPM/LogisticPM")
  )
  .pen <- list(Sel_curve_pen1 = .col("Sel_curve_pen1"),
               Sel_curve_pen2 = .col("Sel_curve_pen2"),
               Sel_curve_pen3 = .col("Sel_curve_pen3"))
  .dir <- if (is.null(.fc$Sel_shape_dir)) rep("Decreasing", nrow(.fc)) else as.character(.fc$Sel_shape_dir)
  .is_incr <- .dir %in% c("Increasing", "increasing", "-1")
  for (nm in names(.sd_specs)) {
    v <- .col(nm)
    if (all(is.na(v))) next
    spec <- .sd_specs[[nm]]
    bad_form <- which(!is.na(v) & !(.fc$Selectivity %in% spec$forms))
    if (length(bad_form) > 0) {
      stop(sprintf(paste0(
        "'%s' (a penalty-SD column) is set on fleet(s) %s whose Selectivity does ",
        "not use that penalty slot as a weight; it applies to %s only (use ",
        "Sel_curve_pen for other forms)."), nm, paste(bad_form, collapse = ", "),
        spec$label), call. = FALSE)
    }
    bad_val <- which(!is.na(v) & !(is.finite(v) & v > 0))
    if (length(bad_val) > 0) {
      stop(sprintf(paste0(
        "'%s' must be a positive standard deviation; fleet(s) %s are non-positive ",
        "or non-finite."), nm, paste(bad_val, collapse = ", ")), call. = FALSE)
    }
    w <- 1 / (2 * v^2)
    # The sign of pen1 sets penalize-decreasing vs -increasing on
    # "NonParametricPM" alone, which is the only form whose shape term branches
    # on it. On every other form the term is applied as written, so a negative
    # weight would reward the shape without bound rather than penalize it.
    if (nm == "Sel_shape_sd") {
      # Only "NonParametricPM" (9) under Sel_shape_mode = "Directional" branches
      # on the sign; its "Smooth" mode is two-sided, like every other form.
      # Checked on every fleet, including Fleet_type = "Off": the direction is a
      # statement about what the FORM can express, and write_data() persists the
      # weight, so exempting an off fleet would leave a sign of the opposite
      # meaning in the workbook for whenever it is switched on.
      dir_ok  <- vapply(seq_len(nrow(.fc)), function(i) .rce_sel_pen1_sign_aware(.fc, i),
                        logical(1))
      bad_dir <- which(!is.na(v) & !dir_ok & .is_incr)
      if (length(bad_dir) > 0) {
        stop(sprintf(paste0(
          "Sel_shape_dir = \"Increasing\" on fleet(s) %s with Selectivity = %s, ",
          "whose shape penalty does not read the sign: an increasing direction ",
          "would reward a decreasing curve rather than penalize it. If you did ",
          "not mean an increasing penalty, drop the column -- the default is ",
          "\"Decreasing\". Only \"NonParametricPM\" under Sel_shape_mode = ",
          "\"Directional\" implements it, and it is not a drop-in; see ",
          "vignette('model-parameterizations')."),
          paste(sprintf("'%s'", .fc$Fleet_name[bad_dir]), collapse = ", "),
          paste(unique(vapply(bad_dir, function(i) .canon_switch(.fc$Selectivity[i], sel_map),
                              character(1))), collapse = "/")), call. = FALSE)
      }
      w <- ifelse(dir_ok & .is_incr, -w, w)
    }
    u <- !is.na(v) & is.na(.pen[[spec$pen]])
    .pen[[spec$pen]][u] <- w[u]
  }
  data_list$fleet_control$Sel_curve_pen1 <- .pen$Sel_curve_pen1
  data_list$fleet_control$Sel_curve_pen2 <- .pen$Sel_curve_pen2
  data_list$fleet_control$Sel_curve_pen3 <- .pen$Sel_curve_pen3
  data_list$fleet_control$Sel_curve_pen1 <- .rce_apply_default(data_list$fleet_control$Sel_curve_pen1, "Sel_curve_pen1", .sch, .np_hake)
  data_list$fleet_control$Sel_curve_pen2 <- .rce_apply_default(data_list$fleet_control$Sel_curve_pen2, "Sel_curve_pen2", .sch, .np_hake)
  data_list$fleet_control$Sel_curve_pen3 <- .rce_apply_default(data_list$fleet_control$Sel_curve_pen3, "Sel_curve_pen3", .sch, .np_hake)
  data_list$fleet_control$Sel_start_year <- .rce_apply_default(data_list$fleet_control$Sel_start_year, "Sel_start_year", .sch)  # per-fleet selectivity penalty start year (NA -> styr); used by LogisticPM
  # The `Sel_pen_first_age` / `Sel_pen_last_age` / `Sel_cap_age` back-compat
  # renames (these bins were named *_age before they were generalised to work on
  # either dimension) are handled by .rce_upgrade_fleet_control_aliases() at the
  # top of switch_check(), from the schema `aliases` field.
  data_list$fleet_control$Sel_pen_first_bin <- .rce_apply_default(data_list$fleet_control$Sel_pen_first_bin, "Sel_pen_first_bin", .sch)  # first bin (age or length) for the non-parametric shape penalty (NA -> bin_first_selected)
  data_list$fleet_control$Sel_pen_last_bin <- .rce_apply_default(data_list$fleet_control$Sel_pen_last_bin, "Sel_pen_last_bin", .sch)  # last (left) bin of the shape-penalty pairs (NA -> nbins-2)
  data_list$fleet_control$Sel_shape_mode <- .rce_apply_default(data_list$fleet_control$Sel_shape_mode, "Sel_shape_mode", .sch)  # shape-penalty mode: "Directional" (default) or "Smooth" (two-sided d^2, RTMB)
  data_list$fleet_control$Sel_avgsel_pen <- .rce_apply_default(data_list$fleet_control$Sel_avgsel_pen, "Sel_avgsel_pen", .sch)  # weight on the AMAK avgsel base-level penalty (type 9): weight * (log(mean(exp(base coffs))))^2; 0 = off (default), 10 matches AMAK
  data_list$fleet_control$Sel_cap_bin <- .rce_apply_default(data_list$fleet_control$Sel_cap_bin, "Sel_cap_bin", .sch)  # NonParametricPM bin cap (NA -> no cap)
  data_list$fleet_control$Selectivity_dimension <- .rce_apply_default(data_list$fleet_control$Selectivity_dimension, "Selectivity_dimension", .sch, conditions = .dflt_when)
  data_list$fleet_control$Comp_distribution <- .rce_apply_default(data_list$fleet_control$Comp_distribution, "Comp_distribution", .sch)
  data_list$fleet_control$CAAL_distribution <- .rce_apply_default(data_list$fleet_control$CAAL_distribution, "CAAL_distribution", .sch, conditions = .dflt_when)
  data_list$fleet_control$Index_distribution <- .rce_apply_default(data_list$fleet_control$Index_distribution, "Index_distribution", .sch)  # survey index likelihood family; default preserves the historical lognormal fit
  # Also fill per-element NAs: setting Index_distribution for one fleet (e.g. only
  # the covariance survey) leaves the other rows NA, which should take the same
  # default. Read it from the schema rather than repeating the literal -- the
  # value belongs in one place, and this line sat three lines below the
  # schema-driven call that already knows it.
  data_list$fleet_control$Index_distribution[is.na(data_list$fleet_control$Index_distribution)] <-
    .sch[["Index_distribution"]]$default
  # Same for Selectivity_dimension: .rce_apply_default() fills a MISSING column,
  # never a blank cell, and a partially-assigned column
  # (fleet_control$Selectivity_dimension[i] <- "Length") is a live idiom in the
  # assessment scripts, so the blank rows are filled rather than refused.
  # Announced fleet by fleet, and under the growth-estimated gate: a blank there
  # is where the author meant "Length", and an age-based curve on a length-based
  # fleet changes the fit without saying so.
  .sel_dim_blank <- which(is.na(data_list$fleet_control$Selectivity_dimension))
  if (length(.sel_dim_blank) > 0) {
    if (isTRUE(.dflt_when$growth_estimated)) {
      .lbl <- data_list$fleet_control$Fleet_name
      .lbl <- if (is.null(.lbl)) .sel_dim_blank else .lbl[.sel_dim_blank]
      message(sprintf(
        "'Selectivity_dimension' is blank for fleet(s) %s in 'fleet_control', assuming '%s'",
        paste(.lbl, collapse = ", "), .sch[["Selectivity_dimension"]]$default))
    }
    data_list$fleet_control$Selectivity_dimension[.sel_dim_blank] <-
      .sch[["Selectivity_dimension"]]$default
  }
  data_list$fleet_control$CAAL_weights <- .rce_apply_default(data_list$fleet_control$CAAL_weights, "CAAL_weights", .sch, conditions = .dflt_when)
  data_list$fleet_control$Comp_accum_young <- .rce_apply_default(data_list$fleet_control$Comp_accum_young, "Comp_accum_young", .sch)  # young-tail composition accumulation bin (NA -> no accumulation)
  data_list$fleet_control$Comp_accum_old <- .rce_apply_default(data_list$fleet_control$Comp_accum_old, "Comp_accum_old", .sch)  # old-tail composition accumulation bin (NA -> no accumulation)
  data_list$fleet_control$Month <- .rce_apply_default(data_list$fleet_control$Month, "Month", .sch)

  # Format upgrade for a pre-4.4 non-parametric workbook, which carries the two
  # AMAK penalty WEIGHTS in Time_varying_sel and Time_varying_sel_sd instead of
  # the Sel_curve_pen columns; move them across. The trigger is the MISSING
  # Sel_curve_pen1 column, not the value: an AMAK shape weight of 4 and the
  # RandomWalk code 4 are the same number, so only the column's absence tells an
  # un-upgraded workbook from a modern one. 0 and 1 are modes, not weights -- a
  # fleet meaning "time-invariant" is not asking for a shape weight of 0.
  np_idx <- data_list$fleet_control$Selectivity %in% c(2, "NonParametric", "Non-parametric", 9, "NonParametricPM", 13, "NonParametricIntegrable")
  .tv_num <- suppressWarnings(as.numeric(data_list$fleet_control$Time_varying_sel))
  # Only the two original codes can appear in a pre-4.4 workbook.
  np_legacy_idx <- data_list$fleet_control$Selectivity %in% c(2, "NonParametric", "Non-parametric", 9, "NonParametricPM")
  legacy <- np_legacy_idx & !.had_sel_curve_pen &
    !is.na(.tv_num) & !(.tv_num %in% c(0, 1))
  if(any(legacy)){
    data_list$fleet_control <- data_list$fleet_control |>
      dplyr::mutate(
        Sel_curve_pen1      = ifelse(legacy, Time_varying_sel, Sel_curve_pen1),
        Sel_curve_pen2      = ifelse(legacy, Time_varying_sel_sd, Sel_curve_pen2),
        Time_varying_sel    = ifelse(legacy, 0, Time_varying_sel),
        Time_varying_sel_sd = ifelse(legacy, 0, Time_varying_sel_sd)
      )
    message("Updating format where 'Selectivity == Non-parametric'. Moving non-parametric penalties to 'Sel_curve_pen1' and 'Sel_curve_pen2'.")
  }

  if(any(np_idx & is.na(data_list$fleet_control$Sel_curve_pen1))) stop("'Sel_curve_pen1' is NA in 'fleet_control' for fleet with non-parametric selectivity")
  if(any(np_idx & is.na(data_list$fleet_control$Sel_curve_pen2))) stop("'Sel_curve_pen2' is NA in 'fleet_control' for fleet with non-parametric selectivity")


  # 2. Sel bins ----
  for(flt in 1:nrow(data_list$fleet_control)){

    sp_idx <- data_list$fleet_control$Species[flt]
    age_selex = data_list$fleet_control$Selectivity_dimension[flt] == "Age"
    selex_text <- ifelse(age_selex, "nages", "nlengths")
    max_bin <- ifelse(age_selex,
                      data_list$nages[sp_idx],
                      data_list$nlengths[sp_idx])

    # The bins the fleet is actually selected over. Below Bin_first_selected the
    # curve is zeroed, so normalizing there divides by ~0, and above the oldest
    # age there is nothing to divide by. Bin_first_selected is a 1-based bin
    # ordinal while Sel_norm_bin is an absolute age, so convert before comparing.
    lo <- .rce_sel_norm_bounds(data_list, flt, age_selex, max_bin)
    hi <- lo$hi; lo <- lo$lo
    is_pm <- data_list$fleet_control$Selectivity[flt] %in% c(11, "LogisticPM")
    unit <- if (isTRUE(age_selex)) "age" else "length bin"

    # A Fixed curve is read from emp_sel_obs and never normalized, so its value
    # may be stale. Every other form is normalized whether or not the fleet is
    # fitted -- selectivity.hpp skips on Selectivity, not Fleet_type -- so an
    # "Off" fleet's bin still indexes the array and must be in range.
    .reads_norm <-
      !data_list$fleet_control$Selectivity[flt] %in% c(0, "Fixed", "Empirical")

    for (col in c("Sel_norm_bin", "Sel_norm_bin_upper")) {
      if (!.reads_norm) next
      raw <- data_list$fleet_control[[col]][flt]
      val <- .rce_sel_norm_code(raw, lo = lo, allow_all = is_pm)
      if (length(attr(val, "unrecognised"))) {
        stop(sprintf(
          "Fleet '%s': '%s' is '%s', which is not a valid %s or one of: %s.",
          data_list$fleet_control$Fleet_name[flt], col, raw, unit,
          paste(c("Max", "Off", if (is_pm) "All"), collapse = ", ")),
          call. = FALSE)
      }
      if (!is.na(val) && !is.na(hi) && val > 0 && val > hi) {
        stop(sprintf(paste0(
          "Fleet '%s': '%s' is %s, past the last %s it is selected over (%s). ",
          "Use %s-%s, or %s."),
          data_list$fleet_control$Fleet_name[flt], col, val, unit, hi, lo, hi,
          paste(c("'Max'", "'Off'", if (is_pm) "'All'"), collapse = " / ")),
          call. = FALSE)
      }

      # Reinterpreting a value the user typed is worth saying out loud: below
      # the first selected bin the curve is zeroed, so it cannot be a reference.
      raw_n <- suppressWarnings(as.numeric(raw))
      if (!is.na(val) && val < 0 && !is.na(raw_n) && raw_n >= 0 && raw_n < lo) {
        message(sprintf(
          "Fleet '%s': '%s' is %s, before the first %s it is selected over (%s); reading it as the maximum.",
          data_list$fleet_control$Fleet_name[flt], col, raw_n, unit, lo))
      }

      # Write the instruction back as the word it means, so a saved workbook
      # says what it does rather than leaving the reader to decode a sentinel.
      # The upper bin has no maximum: a negative there means no range at all.
      data_list$fleet_control[[col]][flt] <- if (is.na(val)) {
        "Off"
      } else if (val < 0) {
        if (col == "Sel_norm_bin_upper") "Off" else if (is_pm) "All" else "Max"
      } else {
        as.character(val)
      }
    }

    # - N bins
    if(any(data_list$fleet_control$N_sel_bins[flt] > max_bin, na.rm = TRUE)){
      data_list$fleet_control$N_sel_bins[flt] <- max_bin
      message(paste0("'N_sel_bins' for fleet ", flt, " is greater than ", selex_text,", setting to ", selex_text))
    }
  }

  # Update switches to text-based if necessary for older files
  data_list <- revert_switches(data_list)

  # Auto-Off fleets with no observations contributing to the likelihood.
  # catch_data, index_data, comp_data,or caal_data should have
  # at least one row with Year > 0 referencing its
  # Fleet_code.
  active_codes <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(integer(0))
    if (!all(c("Fleet_code", "Year") %in% colnames(df))) return(integer(0))
    df <- df[!is.na(df$Year) & df$Year > 0 & !is.na(df$Fleet_code), , drop = FALSE]
    if (nrow(df) == 0) return(integer(0))
    unique(as.integer(df$Fleet_code))
  }
  active <- unique(c(
    active_codes(data_list$catch_data),
    active_codes(data_list$index_data),
    active_codes(data_list$comp_data),
    active_codes(data_list$caal_data)
  ))
  fc <- data_list$fleet_control
  inactive_idx <- which(!fc$Fleet_code %in% active & fc$Fleet_type != "Off")
  if (length(inactive_idx) > 0) {
    message(paste0(
      "Auto-Off fleet(s) with no observations in the likelihood (all Year < 0 or missing): ",
      paste(fc$Fleet_name[inactive_idx], collapse = ", "),
      ". Setting Fleet_type = 'Off'."
    ))
    data_list$fleet_control$Fleet_type[inactive_idx]   <- "Off"
  }

  # A fishery's index is predicted from the year-averaged numbers, which has no
  # instant to be taken at, so Month is inert on its index rows (a survey's is
  # read). Say so rather than let a workbook carry a month that does nothing.
  .fsh_idx <- .fleets_with_index(data_list, fitted_only = TRUE)
  if (length(.fsh_idx)) {
    .fc <- data_list$fleet_control
    .fsh_idx <- .fsh_idx[.fsh_idx %in% .fc$Fleet_code[.fc$Fleet_type %in% c(1, "1", "Fishery")]]
    .im <- data_list$index_data
    .dead_month <- .fsh_idx[vapply(.fsh_idx, function(fl) {
      m <- suppressWarnings(as.numeric(.im$Month[.im$Fleet_code == fl & .im$Year > 0]))
      any(!is.na(m) & m != 0)
    }, logical(1))]
    if (length(.dead_month)) {
      message(paste0(
        "'Month' is not read for the index of fishery fleet(s): ",
        paste(.fc$Fleet_name[match(.dead_month, .fc$Fleet_code)], collapse = ", "),
        ". A fishery index is predicted from the year-averaged numbers ",
        "N*(1-exp(-Z))/Z. Put the index on its own Survey fleet, mirroring the ",
        "fishery's Selectivity_index, if it should be timed to a month."))
    }
  }

  # Default Sel_start_year to the fleet's first year of data (not styr). Earlier
  # deviations have neither data nor a penalty (every cpp selectivity penalty is
  # anchored at start_yr), so they are unidentified. Only affects time-varying
  # selectivity; an explicit per-fleet value wins.
  fleet_obs_yrs <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    if (!all(c("Fleet_code", "Year") %in% colnames(df))) return(NULL)
    df <- df[!is.na(df$Year) & df$Year > 0 & !is.na(df$Fleet_code), , drop = FALSE]
    if (nrow(df) == 0) return(NULL)
    data.frame(Fleet_code = as.integer(df$Fleet_code), Year = as.integer(df$Year))
  }
  obs_yrs <- do.call(rbind, list(
    fleet_obs_yrs(data_list$catch_data),
    fleet_obs_yrs(data_list$index_data),
    fleet_obs_yrs(data_list$comp_data),
    fleet_obs_yrs(data_list$caal_data)
  ))
  if (!is.null(obs_yrs) && nrow(obs_yrs) > 0) {
    first_obs <- tapply(obs_yrs$Year, obs_yrs$Fleet_code, min)
    fc2 <- data_list$fleet_control
    # Fleets sharing a Selectivity_index share ONE selectivity curve, so the start
    # year must be the earliest first-observation across the whole group -- not the
    # fleet's own. E.g. AVO's data starts 2006 but it mirrors the ATS curve, whose
    # data starts 1994; keying off AVO alone would drop the 1994-2005 deviations
    # that the ATS data actually informs.
    sel_grp <- if (!is.null(fc2$Selectivity_index)) as.character(fc2$Selectivity_index) else
      as.character(fc2$Fleet_code)
    fleet_first <- first_obs[as.character(fc2$Fleet_code)]
    grp_first <- tapply(fleet_first, sel_grp, function(v) if (all(is.na(v))) NA_integer_ else min(v, na.rm = TRUE))
    unset_idx <- which(is.na(fc2$Sel_start_year))
    for (i in unset_idx) {
      fy <- grp_first[sel_grp[i]]
      if (!is.na(fy)) {
        data_list$fleet_control$Sel_start_year[i] <- max(as.integer(data_list$styr), as.integer(fy))
      }
    }
  }

  return(data_list)
}


#' Convert integer switches to intuitive text strings. Maintains backwards compatability.
#'
#' @param data_list Rceattle data list
#' @importFrom rlang .data
#' @keywords internal
#' @noRd
# Legacy spellings of a canonical switch VALUE (as distinct from the schema's
# column aliases, which rename a whole column). Upgraded on the way in, because
# switch_check() accepted these while validate_switches() rejected them -- so a
# model written with one loaded, ran switch_check() unchanged, and then failed
# its own data check with "Invalid 'Selectivity'".
.rce_switch_value_aliases <- c(
  "Non-parametric" = "NonParametric"
)


revert_switches <- function(data_list) {

  # Helper: convert integer codes to map names. Numeric-looking strings ("0",
  # " 0", "0.000000000") are canonicalized via as.numeric so they match the
  # integer keys in `map`. Strings that are already a map name (or anything
  # else) pass through unchanged.
  .conv <- function(x, map) {
    x_char <- gsub(" ", "", as.character(x))
    x_num <- suppressWarnings(as.numeric(x_char))
    is_num <- !is.na(x_num)
    x_char[is_num] <- as.character(x_num[is_num])
    # Upgrade a legacy spelling before matching, so it lands on the canonical
    # name rather than passing through to be rejected later.
    alias_hit <- match(x_char, names(.rce_switch_value_aliases))
    if (any(!is.na(alias_hit))) {
      x_char[!is.na(alias_hit)] <-
        unname(.rce_switch_value_aliases[alias_hit[!is.na(alias_hit)]])
    }
    idx <- match(x_char, as.character(map))
    matched <- !is.na(idx)
    if (any(matched)) {
      x_char[matched] <- names(map)[idx[matched]]
    }
    return(x_char)
  }

  # - Fleet switches
  # Guard: Index_distribution and Sel_norm_scope are newer columns; a hand-built
  # fleet_control (or a data_list passed straight to rearrange_data(), bypassing
  # switch_check) may lack them. Fill from the schema so the conversion below
  # works without them, and so the value cannot drift from switch_check()'s.
  .sch_defaults <- .rce_column_schema()
  for (.col in c("Index_distribution", "Sel_norm_scope")) {
    if (is.null(data_list$fleet_control[[.col]]))
      data_list$fleet_control[[.col]] <- .sch_defaults[[.col]]$default
  }

  data_list$fleet_control <- data_list$fleet_control |>
    dplyr::mutate(
      Fleet_type = .conv(.data$Fleet_type, fleet_map),
      Selectivity = .conv(.data$Selectivity, sel_map),
      Time_varying_sel = .conv(.data$Time_varying_sel, tv_sel_map),
      Sel_norm_scope = .conv(.data$Sel_norm_scope, sel_norm_scope_map),
      Catchability = .conv(.data$Catchability, q_map),
      Comp_distribution = .conv(.data$Comp_distribution, comp_loglike_map),
      CAAL_distribution = .conv(.data$CAAL_distribution, comp_loglike_map),
      Index_distribution = .conv(.data$Index_distribution, index_distribution_map)
    )

  # Time_varying_q doubles as an environmental-index column when Catchability
  # is "AR1" or "Environmental", so only convert the rows that hold a switch.
  # `[[` rather than `$` for the reason above.
  tv_q <- data_list$fleet_control[["Time_varying_q"]]
  non_env_idx <- !data_list$fleet_control$Catchability %in% c("AR1", "Environmental")
  if (any(non_env_idx)) {
    tv_q[non_env_idx] <- .conv(tv_q[non_env_idx], tv_q_map)
    data_list$fleet_control[["Time_varying_q"]] <- tv_q
  }

  # - Population dynamics switches
  data_list$initMode <- .conv(data_list$initMode, initMode_map)
  data_list$HCR <- .conv(data_list$HCR, hcr_map)

  return(data_list)
}


#' The value map a column is validated against, from the schema
#'
#' @description
#' `.rce_column_schema()` names the governing map on each `type = "switch"` row
#' (`allowed = "sel_map"`). Resolving it here is what makes the schema
#' authoritative for *which* values a column may take: adding a switch column
#' means adding a row, not another hardcoded map reference in this file.
#'
#' The per-column subset predicate and the wording of the error stay at the call
#' site. They are not uniform, `Time_varying_q` is exempt while `Catchability`
#' is `"Environmental"`, `Catchability` itself allows `NA`, and the newer
#' columns may be absent entirely on a list `switch_check()` has not yet seen,
#' and flattening that into one generic loop would lose real behaviour.
#'
#' @param col Canonical column name.
#' @return The named map, or `NULL` if the schema declares none.
#' @keywords internal
#' @noRd
.rce_allowed_map <- function(col) {
  row <- .rce_column_schema()[[col]]
  if (is.null(row)) {
    stop("Internal: '", col, "' is not in the column schema.", call. = FALSE)
  }
  nm <- row$allowed
  if (is.null(nm) || all(is.na(nm))) {
    # Returning NULL here would make every value invalid and print
    # "Please use one of:" with nothing after it -- loud, but pointing the user
    # at their data instead of at the schema row that is actually missing.
    stop("Internal: the schema declares no `allowed` map for '", col, "'.",
         call. = FALSE)
  }
  get(nm, envir = asNamespace("Rceattle"))
}


#' Validates switches are correct
#'
#' @param data_list Rceattle data list
#'
#' @keywords internal
#' @noRd
validate_switches <- function(data_list = NULL){
  errors <- character(0)

  # Governing maps, resolved from the schema's `allowed` field.
  fleet_map              <- .rce_allowed_map("Fleet_type")
  sel_map                <- .rce_allowed_map("Selectivity")
  tv_sel_map             <- .rce_allowed_map("Time_varying_sel")
  sel_norm_scope_map     <- .rce_allowed_map("Sel_norm_scope")
  q_map                  <- .rce_allowed_map("Catchability")
  tv_q_map               <- .rce_allowed_map("Time_varying_q")
  comp_loglike_map       <- .rce_allowed_map("Comp_distribution")
  index_distribution_map <- .rce_allowed_map("Index_distribution")
  sel_dimension_map      <- .rce_allowed_map("Selectivity_dimension")
  sel_shape_dir_map      <- .rce_allowed_map("Sel_shape_dir")
  sel_shape_mode_map     <- .rce_allowed_map("Sel_shape_mode")

  # Validate fleet_control inputs ----
  # The newer switch columns can be absent here: data_check() is callable on a
  # data list read straight from a workbook, before switch_check() has filled
  # its schema defaults. A column that is not there will be defaulted to a valid
  # value, so there is nothing to validate -- skip the check rather than fail on
  # a missing column.
  .fc_none <- data_list$fleet_control[0, , drop = FALSE]
  .fc_has <- function(col) !is.null(data_list$fleet_control[[col]])

  # Every predicate below exempts an Off fleet. data_check() is callable on a
  # list read straight from a workbook, where Fleet_type may still be the
  # integer 0 -- and `0 != "Off"` coerces to `"0" != "Off"`, which is TRUE, so a
  # bare string comparison would validate the fleets it means to skip.
  # Canonicalize once, and compare against that.
  .fleet_type_canon <- .canon_switch(data_list$fleet_control$Fleet_type, fleet_map)
  data_list$fleet_control$.rce_is_on <- .fleet_type_canon != "Off"

  invalid_flt_type <- data_list$fleet_control |>
    dplyr::filter(!.data$Fleet_type %in% c(fleet_map, names(fleet_map)))

  invalid_sel <- data_list$fleet_control |>
    dplyr::filter(.data$.rce_is_on & !.data$Selectivity %in% c(sel_map, names(sel_map)))

  invalid_tv_sel <- data_list$fleet_control |>
    dplyr::filter(.data$.rce_is_on & !.data$Time_varying_sel %in% c(tv_sel_map, names(tv_sel_map)))

  invalid_q <- data_list$fleet_control |>
    dplyr::filter(!.data$Catchability %in% c(NA, q_map, names(q_map)))

  invalid_tv_q <- data_list$fleet_control |>
    dplyr::filter(.data$.rce_is_on & .data$Catchability != "Environmental" & !.data$Time_varying_q %in% c(NA, tv_q_map, names(tv_q_map)))

  invalid_sel_norm_scope <- if (.fc_has("Sel_norm_scope")) {
    data_list$fleet_control |>
      dplyr::filter(.data$.rce_is_on &
                      !.data$Sel_norm_scope %in% c(sel_norm_scope_map, names(sel_norm_scope_map)))
  } else .fc_none

  # Selectivity_dimension and the two Sel_shape_* columns are validated here
  # because nothing downstream re-checks them: a typo resolves to NA rather than
  # erroring, leaving a missing selectivity dimension or penalty mode. Each is
  # checked against what its CONSUMER implements, not the map alone --
  # rearrange_data() matches Selectivity_dimension on the exact strings and
  # yields NA for anything else, an integer included, so the integer side must
  # not validate; the Sel_shape_* columns are read case-insensitively in one
  # spelling each ("Smooth"/"smooth", not "SMOOTH"), and
  # Sel_shape_dir also accepts "-1" (the ADMB sign convention).
  invalid_sel_dim <- if (.fc_has("Selectivity_dimension")) {
    data_list$fleet_control |>
      dplyr::filter(.data$.rce_is_on & !is.na(.data$Selectivity_dimension) &
                      !.data$Selectivity_dimension %in% names(sel_dimension_map))
  } else .fc_none

  .dir_ok <- c(names(sel_shape_dir_map), "increasing", "decreasing", "-1", "0")
  invalid_sel_shape_dir <- if (.fc_has("Sel_shape_dir")) {
    data_list$fleet_control |>
      dplyr::filter(.data$.rce_is_on & !is.na(.data$Sel_shape_dir) &
                      !.data$Sel_shape_dir %in% .dir_ok)
  } else .fc_none

  .mode_ok <- c(names(sel_shape_mode_map), "smooth", "directional",
                unname(sel_shape_mode_map))
  invalid_sel_shape_mode <- if (.fc_has("Sel_shape_mode")) {
    data_list$fleet_control |>
      dplyr::filter(.data$.rce_is_on & !is.na(.data$Sel_shape_mode) &
                      !.data$Sel_shape_mode %in% .mode_ok)
  } else .fc_none

  invalid_comp_ll <- data_list$fleet_control |>
    dplyr::filter(.data$.rce_is_on & !.data$Comp_distribution %in% c(comp_loglike_map, names(comp_loglike_map)))

  invalid_caal_ll <- data_list$fleet_control |>
    dplyr::filter(.data$.rce_is_on & !.data$CAAL_distribution %in% c(comp_loglike_map, names(comp_loglike_map)))

  invalid_index_ll <- if (.fc_has("Index_distribution")) {
    data_list$fleet_control |>
      dplyr::filter(.data$.rce_is_on &
                      !.data$Index_distribution %in% c(index_distribution_map, names(index_distribution_map)))
  } else .fc_none

  # Throw clear errors to guide the user
  if(nrow(invalid_flt_type) > 0) {
    errors <- c(errors, paste("Invalid 'Fleet_type' specified for fleets:",
                              paste(invalid_flt_type$Fleet_name, collapse = ", "),
                              ".\nPlease use an integer code ",paste(range(fleet_map), collapse = ":")," or one of:",
                              paste(names(fleet_map), collapse = ", ")))
  }

  if(nrow(invalid_sel) > 0) {
    errors <- c(errors, paste("Invalid 'Selectivity' specified for fleets:",
                              paste(invalid_sel$Fleet_name, collapse = ", "),
                              ".\nPlease use an integer code ",paste(range(sel_map), collapse = ":")," or one of:",
                              paste(names(sel_map), collapse = ", ")))
  }

  if(nrow(invalid_tv_sel) > 0) {
    errors <- c(errors, paste("Invalid 'Time_varying_sel' specified for fleets:",
                              paste(invalid_tv_sel$Fleet_name, collapse = ", "),
                              ".\nPlease use an integer code ",paste(range(tv_sel_map), collapse = ":")," or one of:",
                              paste(names(tv_sel_map), collapse = ", ")))
  }

  if(nrow(invalid_sel_norm_scope) > 0) {
    errors <- c(errors, paste("Invalid 'Sel_norm_scope' specified for fleets:",
                              paste(invalid_sel_norm_scope$Fleet_name, collapse = ", "),
                              ".\nPlease use one of:",
                              paste(names(sel_norm_scope_map), collapse = ", ")))
  }

  if(nrow(invalid_q) > 0) {
    errors <- c(errors, paste("Invalid 'Catchability' specified for fleets:",
                              paste(invalid_q$Fleet_name, collapse = ", "),
                              ".\nPlease use an integer code ", paste(range(q_map), collapse = ":")," or one of:",
                              paste(names(q_map), collapse = ", ")))
  }

  if(nrow(invalid_tv_q) > 0) {
    errors <- c(errors, paste("Invalid 'Time_varying_q' specified for fleets:",
                              paste(invalid_tv_q$Fleet_name, collapse = ", "),
                              ".\nPlease use an integer code ", paste(range(tv_q_map), collapse = ":")," or one of:",
                              paste(names(tv_q_map), collapse = ", ")))
  }

  if(nrow(invalid_sel_dim) > 0) {
    errors <- c(errors, paste("Invalid 'Selectivity_dimension' specified for fleets:",
                              paste(invalid_sel_dim$Fleet_name, collapse = ", "),
                              ".\nPlease use one of:",
                              paste(names(sel_dimension_map), collapse = ", ")))
  }

  if(nrow(invalid_sel_shape_dir) > 0) {
    errors <- c(errors, paste("Invalid 'Sel_shape_dir' specified for fleets:",
                              paste(invalid_sel_shape_dir$Fleet_name, collapse = ", "),
                              ".\nPlease use one of:",
                              paste(names(sel_shape_dir_map), collapse = ", ")))
  }

  if(nrow(invalid_sel_shape_mode) > 0) {
    errors <- c(errors, paste("Invalid 'Sel_shape_mode' specified for fleets:",
                              paste(invalid_sel_shape_mode$Fleet_name, collapse = ", "),
                              ".\nPlease use one of:",
                              paste(names(sel_shape_mode_map), collapse = ", ")))
  }

  if(nrow(invalid_comp_ll) > 0) {
    errors <- c(errors, paste("Invalid 'Comp_distribution' specified for fleets:",
                              paste(invalid_comp_ll$Fleet_name, collapse = ", "),
                              ".\nPlease use an integer code ", paste(range(comp_loglike_map), collapse = ":")," or one of:",
                              paste(names(comp_loglike_map), collapse = ", ")))
  }

  if(nrow(invalid_caal_ll) > 0) {
    errors <- c(errors, paste("Invalid 'CAAL_distribution' specified for fleets:",
                              paste(invalid_caal_ll$Fleet_name, collapse = ", "),
                              ".\nPlease use an integer code ", paste(range(comp_loglike_map), collapse = ":")," or one of:",
                              paste(names(comp_loglike_map), collapse = ", ")))
  }

  if(nrow(invalid_index_ll) > 0) {
    errors <- c(errors, paste("Invalid 'Index_distribution' specified for fleets:",
                              paste(invalid_index_ll$Fleet_name, collapse = ", "),
                              ".\nPlease use an integer code ", paste(range(index_distribution_map), collapse = ":")," or one of:",
                              paste(names(index_distribution_map), collapse = ", ")))
  }

  # Validate pop dy controls ----
  # * initMode ----
  invalid_initMode <- (!data_list$initMode %in% c(initMode_map, names(initMode_map)))

  if(sum(invalid_initMode) > 0) {
    errors <- c(errors, paste0("Invalid 'initMode' specified: ",
                              paste(unique(data_list$initMode[invalid_initMode]), collapse = ", "),
                              ".\nPlease use an integer code ", paste(range(initMode_map), collapse = ":"), " or one of: ",
                              paste(names(initMode_map), collapse = ", ")))
  }

  # * HCR ----
  # No HCR when the list came from a workbook -- `NULL %in% ...` is logical(0),
  # which propagates through `&` and leaves the `if` below with nothing to test.
  if (length(data_list$HCR)) {
    invalid_hcr <- (!data_list$HCR %in% c(hcr_map, names(hcr_map)))

    if(sum(invalid_hcr) > 0) {
      errors <- c(errors, paste0("Invalid 'HCR' specified: ",
                                paste(unique(data_list$HCR[invalid_hcr]), collapse = ", "),
                                ".\nPlease use an integer code ", paste(range(hcr_map), collapse = ":"), " or one of: ",
                                paste(names(hcr_map), collapse = ", ")))
    }

    if (isTRUE(data_list$msmMode > 0) & !data_list$HCR %in% c("NoFishing", "CMSY", "ConstantF", "ConstantFSSB", "PFMC")) {
      errors <- c(errors, 'Only HCRs "NoFishing" (0), "CMSY" (1), "ConstantF" (2), "ConstantFSSB" (3), or "PFMC" (6) work in multi-species mode currently')
    }
  }

  return(errors)
}


#' Convert intuitive text strings to integer switches for TMB
#'
#' @param data_list Rceattle data list
#'
#' @importFrom rlang .data
#' @keywords internal
#' @noRd
convert_switches <- function(data_list) {

  # Helper: convert a single string value using a map, pass integers through unchanged
  .conv_single <- function(x, map) {
    if (is.character(x) && x %in% names(map)) unname(map[[x]]) else x
  }
  .conv <- Vectorize(.conv_single, vectorize.args = "x", USE.NAMES = FALSE)

  # A factor reaches .conv() as its level index, not its label, and all nine
  # switch columns below would take it. rearrange_data() is exported, so this
  # entry point needs the rule too.
  data_list$fleet_control <-
    .rce_defactor_fleet_control(data_list$fleet_control)

  # Fleet controls ----
  # Default the newer switch columns when a caller supplies a fleet_control that
  # never went through switch_check(): a direct rearrange_data() call (exported,
  # and documented to work on a list read straight from a workbook), or a
  # data_check() before fitting. Every pre-5.8.0 data list, the bundled ones
  # included, lacks Sel_norm_scope, and the .data pronoun below fails on a
  # missing column. Defaults come from the schema so they cannot drift from
  # switch_check()'s.
  .sch_defaults <- .rce_column_schema()
  for (.col in c("Index_distribution", "Sel_norm_scope")) {
    if (is.null(data_list$fleet_control[[.col]]))
      data_list$fleet_control[[.col]] <- .sch_defaults[[.col]]$default
  }
  # A blank cell in a supplied Sel_norm_scope column would reach TMB as NA, which
  # the C++ reads as WithinSex -- the opposite of the documented default. Fill it
  # here as switch_check() does, so both entry points agree.
  data_list$fleet_control$Sel_norm_scope[is.na(data_list$fleet_control$Sel_norm_scope)] <-
    .sch_defaults[["Sel_norm_scope"]]$default
  # If vector is a string that exists in our map, replace it with the integer.
  data_list$fleet_control <- data_list$fleet_control %>%
    dplyr::mutate(
      Fleet_type = .conv(.data$Fleet_type, fleet_map),
      Selectivity = .conv(.data$Selectivity, sel_map),
      Time_varying_sel = .conv(.data$Time_varying_sel, tv_sel_map),
      Sel_norm_scope = .conv(.data$Sel_norm_scope, sel_norm_scope_map),
      Catchability = .conv(.data$Catchability, q_map),
      Time_varying_q = .conv(.data$Time_varying_q, tv_q_map),
      Comp_distribution = .conv(.data$Comp_distribution, comp_loglike_map),
      CAAL_distribution = .conv(.data$CAAL_distribution, comp_loglike_map),
      Index_distribution = .conv(.data$Index_distribution, index_distribution_map)
    ) %>%
    # CRITICAL: Force columns back to integers so TMB doesn't crash expecting ints but getting chars
    dplyr::mutate(
      Fleet_type = as.integer(.data$Fleet_type),
      Selectivity = as.integer(.data$Selectivity),
      Time_varying_sel = as.integer(.data$Time_varying_sel),
      Sel_norm_scope = as.integer(.data$Sel_norm_scope),
      Catchability = as.integer(.data$Catchability),
      Comp_distribution = as.integer(.data$Comp_distribution),
      CAAL_distribution = as.integer(.data$CAAL_distribution),
      Index_distribution = as.integer(.data$Index_distribution)
    )

  # `Time_varying_q` is excluded from the blanket coercion above: for "AR1" (6)
  # and "Environmental" (5) catchability it holds 1-based `env_data` column
  # indices rather than a tv_q_map mode, and "Environmental" may carry a
  # comma-separated list ("1,3"). Coerce only the rows holding a mode code.
  .q_env_rows <- data_list$fleet_control$Catchability %in% c(5L, 6L)
  if (!any(.q_env_rows)) {
    data_list$fleet_control$Time_varying_q <-
      as.integer(data_list$fleet_control$Time_varying_q)
  } else if (any(!.q_env_rows)) {
    data_list$fleet_control$Time_varying_q[!.q_env_rows] <-
      as.integer(data_list$fleet_control$Time_varying_q[!.q_env_rows])
  }

  # Pop dy controls ----
  data_list$initMode <- as.integer(.conv(data_list$initMode, initMode_map))
  data_list$HCR <- as.integer(.conv(data_list$HCR, hcr_map))

  return(data_list)
}


#' Report composition weights that are on the log scale
#'
#' @description
#' A composition weight is read on a scale set by its distribution. Under a
#' multinomial it multiplies the input sample size directly; under a
#' Dirichlet-multinomial the model uses `exp(weight)` (`ceattle.cpp`,
#' `DM_pars_comp = exp(comp_weights)`), so the column holds the LOG of the
#' starting weight. A value of 1 is therefore a starting weight of e, and a
#' weight of 1 is written as 0.
#'
#' `write_template()` seeds these columns with 1, so a model built from the
#' template and switched to a Dirichlet-multinomial starts at e without anything
#' saying so. This reports it once per fit, naming the fleets, so the value is a
#' choice rather than an inherited default.
#'
#' Called from `fit_mod()`, not from `switch_check()`: `switch_check()` is a
#' re-entrant normalizer that `build_params()` and `build_map()` also call, so
#' reporting there printed the same message three times per fit and roughly
#' twenty times per `retrospective()`.
#'
#' Only an exact 1 is reported, the model value. Any other number was typed
#' deliberately and needs no comment. Off fleets and fleets holding no data for
#' the composition in question are skipped: their weight is never read.
#'
#' @param data_list a data list, after `switch_check()` has resolved the
#'   distribution columns.
#' @return `data_list`, unchanged. This reports; it never edits.
#' @keywords internal
#' @noRd
.rce_flag_dm_weight_scale <- function(data_list) {
  is_dm <- function(x) !is.na(x) & as.character(x) %in% c("1", "DirichletMultinomial")

  # Fleet codes with at least one row actually fitted (Year > 0 marks a fitted
  # row; a negative year is carried but not fitted).
  fitted_codes <- function(df) {
    if (is.null(df) || !nrow(df)) return(integer(0))
    if (!all(c("Fleet_code", "Year") %in% names(df))) return(integer(0))
    unique(as.integer(df$Fleet_code[!is.na(df$Year) & df$Year > 0]))
  }

  fc <- data_list$fleet_control
  if (!is.null(fc)) {
    on <- if (is.null(fc$Fleet_type)) rep(TRUE, nrow(fc)) else
      as.character(fc$Fleet_type) != "Off"

    for (spec in list(list("Comp_distribution", "Comp_weights", data_list$comp_data),
                      list("CAAL_distribution", "CAAL_weights", data_list$caal_data))) {
      dist_col <- spec[[1]]; wt_col <- spec[[2]]
      if (is.null(fc[[dist_col]]) || is.null(fc[[wt_col]])) next
      has_data <- fc$Fleet_code %in% fitted_codes(spec[[3]])
      hit <- which(on & has_data & is_dm(fc[[dist_col]]) &
                     !is.na(fc[[wt_col]]) & fc[[wt_col]] == 1)
      if (length(hit)) {
        who <- if (!is.null(fc$Fleet_name)) fc$Fleet_name[hit] else fc$Fleet_code[hit]
        message("'", wt_col, "' is 1 on Dirichlet-multinomial fleet(s) ",
                paste(who, collapse = ", "),
                ". That likelihood reads the column as a log, so this is a starting ",
                "weight of e (", signif(exp(1), 4), "). Use 0 for a starting weight of 1.")
      }
    }
  }

  dw <- data_list$Diet_comp_weights
  dd <- data_list$Diet_distribution
  # Length-guarded: this runs before data_check(), and recycling a ragged pair
  # here would raise a warning ahead of the check that names the real problem.
  if (!is.null(dw) && !is.null(dd) && length(dw) == length(dd)) {
    hit <- which(is_dm(dd) & !is.na(dw) & dw == 1)
    if (length(hit)) {
      who <- if (!is.null(data_list$spnames)) data_list$spnames[hit] else hit
      message("'Diet_comp_weights' is 1 for Dirichlet-multinomial predator(s) ",
              paste(who, collapse = ", "),
              ". That likelihood reads the value as a log, so this is a starting ",
              "weight of e (", signif(exp(1), 4), "). Use 0 for a starting weight of 1.")
    }
  }

  data_list
}


#' Coerce a string/integer switch argument to its canonical integer code.
#'
#' Shared validation for the integer-returning `build_*()` switch arguments
#' (`srr_fun`, `M1_model`, `sd_plus_group`). Accepts a canonical string key from `map` or a
#' legacy integer code (a value of `map`, or a `deprecated` code); errors on
#' anything else and calls `warn_fn(int)` when a deprecated code is supplied.
#'
#' @param x user input (length-1+ character or numeric).
#' @param map named integer vector mapping canonical string -> integer code.
#' @param what argument name, used in messages.
#' @param deprecated integer codes that are accepted but soft-deprecated.
#' @param warn_fn `function(int)` emitting the deprecation warning, or `NULL`.
#' @param length_exact_one require length 1 (`srr_fun`) vs length >= 1 (`M1_model`).
#' @param legacy_note text appended inside the integer-range error message.
#' @keywords internal
#' @noRd
.coerce_switch_arg <- function(x, map, what,
                               deprecated = integer(0),
                               warn_fn = NULL,
                               length_exact_one = FALSE,
                               legacy_note = "") {
  if (length_exact_one) {
    if (length(x) != 1L) {
      stop(sprintf("`%s` must be length 1", what), call. = FALSE)
    }
  } else if (length(x) == 0L) {
    stop(sprintf("`%s` must have length >= 1", what), call. = FALSE)
  }
  if (is.numeric(x)) {
    int <- as.integer(x)
    allowed <- c(unname(map), deprecated)
    if (anyNA(int) || any(!int %in% allowed)) {
      stop(sprintf(
        "integer `%s` must be one of: %s (= %s%s)",
        what,
        paste(map, collapse = ", "),
        paste(names(map), collapse = "/"),
        legacy_note), call. = FALSE)
    }
    if (length(deprecated) > 0L && any(int %in% deprecated) && !is.null(warn_fn)) {
      warn_fn(int)
    }
    return(int)
  }
  if (is.character(x)) {
    bad <- setdiff(x, names(map))
    if (length(bad) > 0L) {
      stop(sprintf(
        "unknown `%s` value(s): %s; allowed: %s",
        what,
        paste(unique(bad), collapse = ", "),
        paste(names(map), collapse = ", ")), call. = FALSE)
    }
    return(unname(map[x]))
  }
  stop(sprintf("`%s` must be a string or integer; got %s",
               what, class(x)[1]), call. = FALSE)
}
