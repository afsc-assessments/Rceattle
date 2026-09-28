#' Selectivity parameters that accept a linkage
#'
#' Slot names shared across the parametric forms, plus DoubleNormal aliases.
#' @keywords internal
#' @noRd
SEL_LINKAGE_PARAMS <- c("slp_asc", "slp_desc", "inf_asc", "inf_desc", "coff",
                        "sigma_asc", "sigma_desc", "peak", "right_floor",
                        "apical")


#' @keywords internal
#' @noRd
.validate_sel_linkages <- function(linkages) {
  .validate_process_linkages(linkages, SEL_LINKAGE_PARAMS, "sel")
}


#' Selectivity specification
#'
#' @description
#' Holds environmental linkages on selectivity parameters. The effect on a
#' parameter is written as a formula and composes additively with any
#' `Time_varying_sel` process error on the same fleet (the two are separate
#' mechanisms: a covariate effect versus a deviation).
#'
#' The parameter names are the shape parameters of the parametric selectivity
#' forms:
#' \describe{
#'   \item{`slp_asc`, `slp_desc`}{ascending / descending logistic slope (log
#'     scale); for a double-normal the ascending / descending width, aliased
#'     `sigma_asc` / `sigma_desc`.}
#'   \item{`inf_asc`, `inf_desc`}{ascending / descending inflection age/length
#'     (natural scale); for a double-normal the peak and the logit right-floor,
#'     aliased `peak` / `right_floor`.}
#'   \item{`coff`}{non-parametric selectivity-at-bin coefficients.}
#'   \item{`apical`}{a multiplier on one sex's whole curve (log scale),
#'     applied after the form and before normalization, so every estimated
#'     form takes it. Name the fleet and the sex that holds it
#'     (`by = ~ fleet + sex`, `fleet = 3`, `sex = "male"`), as Stock
#'     Synthesis's male-offset option does; the other sex is the reference.
#'     See Details.}
#' }
#'
#' Every parameter but `apical` accepts `link = "log"` (multiplicative on the
#' natural parameter) or `link = "identity"` (additive), like the other
#' processes; `apical` is a multiplier already and takes `"log"` only.
#'
#' @details
#' **The `apical` offset.** Fishing mortality is one `log_F` per fleet and
#' year shared by the sexes, so a sex difference in F can only come from
#' selectivity, and no form has a height parameter: the logistic family and
#' DoubleNormal peak at 1 for every sex, the non-parametric forms re-centre
#' each sex, Hake normalizes each sex by its own maximum. `apical` multiplies
#' one sex's curve by `exp(log_sel_apical)`, bin by bin. That equals the ratio
#' of the sexes' peak heights only where their shapes peak equally (the
#' logistic family on an age axis); for a dome with sex-specific shape, read
#' it as the multiplier on that sex's curve and take the peak ratio from
#' `fit$quantities$sel_at_age`. Only the contrast between the sexes is
#' identified (the common level is `log_F`), so one sex holds it and the fit
#' is refused if both do, if no fleet or no sex is named, if the species has
#' one sex, on a `Fixed`, AR1 or `Fleet_type = "Off"` fleet, on a fleet that
#' shares another's `Selectivity_index` block, or under `link = "identity"`,
#' which could drive the multiplier negative; use the default `link = "log"`.
#' It is also refused where
#' `Sel_norm_scope = "WithinSex"` normalization would divide it straight back
#' out; use `"AcrossSexes"`, under which the more-selected sex peaks at 1, or
#' turn `Sel_norm_bin` off. The contrast is informed only by joint composition
#' (`comp_data$Sex = 3`); with single-sex compositions it rests on its prior,
#' and `fit_mod()` warns when a fleet has neither. Read the fitted multiplier
#' with `exp(fit$estimated_params$log_sel_apical[fleet, sex])`, and the
#' realized ratio of the sexes' maxima from `fit$quantities$sel_at_age`.
#' Naming `fleet` and `sex` is enough: `by` defaults to `~ fleet + sex` for
#' this parameter.
#' An intercept prior is on the multiplier's natural scale (`lognormal()`
#' centred on 1 means no offset). Like every selectivity linkage, a covariate
#' on it acts in the hindcast years; projection years hold the last hindcast
#' year's curve.
#'
#' **Priors on a selectivity parameter.** An intercept-only formula (`~ 1`) with
#' a `priors` entry places a prior on the selectivity parameter itself (no
#' year-to-year offset is added). Read the prior on the parameter's own scale:
#' the slopes (`slp_asc` / `slp_desc`) are on the log scale (use `lognormal()`),
#' the inflections (`inf_asc` / `inf_desc`) on the natural scale (use `normal()`).
#' See Examples for a normal prior on the ascending inflection. This mirrors the
#' prior-only [build_composition()] path.
#'
#' A selectivity prior targets one parameter, so in a two-sex model an
#' unstratified `~ 1` prior constrains sex 1 only, use `by = ~ sex` for a
#' per-sex prior. An `init` on a selectivity intercept has no effect (the
#' starting value comes from the data), and a prior on the double-normal
#' `right_floor` is not supported.
#' Fleets sharing a `Selectivity_index` estimate one parameter block, so place
#' the prior on the group's lead fleet (its first fleet that is not `Off`); a
#' prior on a follower would penalize the shared block once per sharing fleet.
#' A prior on an `Off` fleet is refused: its selectivity is not estimated.
#'
#' @param linkages Optional named list of [linkage_spec()] objects keyed by
#'   selectivity parameter. Coefficients are per fleet by default
#'   (`by = ~ fleet`); use the `fleet` argument of [linkage_spec()] to restrict
#'   a spec to particular fleets.
#'
#' @return A list of selectivity settings for [fit_mod()].
#'
#' @examples
#' \donttest{
#' # A cold-pool effect on the ascending inflection of a logistic fleet
#' build_selectivity(linkages = list(
#'   inf_asc = linkage_spec(~ cold_pool, by = ~ fleet)))
#'
#' # A normal prior on the ascending inflection (intercept-only formula)
#' build_selectivity(linkages = list(
#'   inf_asc = linkage_spec(~ 1, priors = list(`(Intercept)` = normal(0, 3)))))
#' }
#'
#' @export
build_selectivity <- function(linkages = NULL) {
  list(linkages = .validate_sel_linkages(linkages))
}


# Map a selectivity linkage param name to (array, 1-based slot). log_sel_slp
# and sel_inf are [2, fleet, sex]; sel_coff is [fleet, sex, bin].
.SEL_PARAM_TO_SLOT <- list(
  slp_asc     = list(arr = "log_sel_slp", slot = 1L),
  slp_desc    = list(arr = "log_sel_slp", slot = 2L),
  sigma_asc   = list(arr = "log_sel_slp", slot = 1L),
  sigma_desc  = list(arr = "log_sel_slp", slot = 2L),
  inf_asc     = list(arr = "sel_inf",     slot = 1L),
  inf_desc    = list(arr = "sel_inf",     slot = 2L),
  peak        = list(arr = "sel_inf",     slot = 1L),
  right_floor = list(arr = "sel_inf",     slot = 2L),
  coff        = list(arr = "sel_coff",    slot = NA_integer_),
  apical      = list(arr = "log_sel_apical", slot = NA_integer_)  # [fleet, sex]
)


# Selectivity forms whose consume site reads the linkage offset: every
# PARAMETRIC form. DoubleNormal reuses the slp/inf slots (peak/sigma/floor
# aliases) and LogisticPM's multiplicative deviates carry the offset inside
# their exp. The non-parametric forms are excluded on purpose -- see the
# `coff` note in .check_sel_linkage_support().
.SEL_LINKAGE_WIRED_FORMS <- c("Logistic", "DoubleLogistic", "DescendingLogistic",
                              "DoubleNormal", "LogisticPM")
.SEL_LINKAGE_WIRED_PARAMS <- c("slp_asc", "slp_desc", "inf_asc", "inf_desc",
                               "sigma_asc", "sigma_desc", "peak", "right_floor",
                               "apical")


#' Reject selectivity linkages the model does not yet consume
#'
#' @param linkage_table pooled linkage table (may be NULL / empty).
#' @param fleet_control the fleet control table.
#' @param nsex sexes per species (`data_list$nsex`); NULL skips the sex checks.
#' @param comp_data the composition data; NULL skips the joint-composition check.
#' @return invisibly NULL; errors on an unsupported sel linkage.
#' @keywords internal
#' @noRd
.check_sel_linkage_support <- function(linkage_table, fleet_control, nsex = NULL,
                                       comp_data = NULL) {
  if (is.null(linkage_table) || nrow(linkage_table) == 0L) return(invisible())
  sel <- linkage_table[linkage_table$process == "sel", , drop = FALSE]
  if (nrow(sel) == 0L) return(invisible())

  # `apical` multiplies the finished curve, so it needs no form-specific consume
  # site and is checked on its own terms below; the form check is for the rest.
  ap  <- sel[sel$param == "apical", , drop = FALSE]
  sel <- sel[sel$param != "apical", , drop = FALSE]
  if (nrow(ap) > 0L) .check_sel_apical_rows(ap, fleet_control, nsex, comp_data)
  if (nrow(sel) == 0L) return(invisible())

  bad_param <- setdiff(unique(sel$param), .SEL_LINKAGE_WIRED_PARAMS)
  if (length(bad_param) > 0) {
    extra <- if ("coff" %in% bad_param) paste0(
      "\n  `coff` (non-parametric selectivity) cannot carry a linkage: those ",
      "forms mean-centre their coefficients each year, so a per-year offset ",
      "applied across all bins cancels exactly. A meaningful effect would need ",
      "a per-bin covariate, which the formula grammar does not express.") else ""
    stop(sprintf(
      paste0("selectivity linkage parameter(s) not supported: %s.\n",
             "  Supported: %s.%s"),
      paste(bad_param, collapse = ", "),
      paste(.SEL_LINKAGE_WIRED_PARAMS, collapse = ", "), extra), call. = FALSE)
  }

  # Every fleet a sel row targets must use a wired selectivity form.
  flts <- unique(sel$fleet)
  flts <- flts[!is.na(flts)]
  if (length(flts) == 0L) flts <- seq_len(nrow(fleet_control))  # NA = all fleets
  forms <- as.character(fleet_control$Selectivity[flts])
  bad_flt <- flts[!forms %in% .SEL_LINKAGE_WIRED_FORMS]
  if (length(bad_flt) > 0) {
    stop(sprintf(
      paste0("selectivity linkage on fleet(s) %s whose form (%s) is not yet ",
             "wired for linkages.\n  Wired forms: %s."),
      paste(fleet_control$Fleet_name[bad_flt], collapse = ", "),
      paste(unique(as.character(fleet_control$Selectivity[bad_flt])),
            collapse = ", "),
      paste(.SEL_LINKAGE_WIRED_FORMS, collapse = ", ")), call. = FALSE)
  }

  # A PRIOR on a selectivity intercept re-targets the base parameter, whose scale
  # depends on the fleet's form and whose ownership depends on Selectivity_index.
  # Two cases the re-target cannot yet express correctly are rejected up front
  # (a covariate linkage on the same slot is fine -- only priors are affected).
  prior_rows <- sel[!is.na(sel$prior_family) & sel$prior_family != "none", ,
                    drop = FALSE]
  if (nrow(prior_rows) > 0L) {
    row_flt <- function(f) if (is.na(f)) 1L else as.integer(f)   # NA fleet = cell 1 (cpp default)

    # (a) DoubleNormal stores sel_inf(1) as logit(right_floor), so a natural-scale
    # prior on `inf_desc` / `right_floor` would be evaluated on the logit scale.
    # Reject until the logit transform is wired -- the ascending peak (`inf_asc`)
    # and the sigmas/slopes (log scale) are unaffected.
    dn <- prior_rows[prior_rows$param %in% c("inf_desc", "right_floor"), , drop = FALSE]
    dn_flt <- unique(vapply(dn$fleet, row_flt, integer(1)))
    dn_flt <- dn_flt[as.character(fleet_control$Selectivity[dn_flt]) == "DoubleNormal"]
    if (length(dn_flt) > 0L) {
      stop(sprintf(paste0(
        "prior on `inf_desc` / `right_floor` for DoubleNormal fleet(s) %s is not ",
        "supported: that slot holds logit(right_floor), so a natural-scale prior ",
        "would be applied on the logit scale. Prior the ascending peak / sigmas ",
        "instead."),
        paste(fleet_control$Fleet_name[dn_flt], collapse = ", ")), call. = FALSE)
    }

    prior_flt <- unique(vapply(prior_rows$fleet, row_flt, integer(1)))

    # (b) An "Off" fleet's selectivity parameters are mapped off, so a prior on
    # one is evaluated against a fixed value: it adds a constant to the objective
    # and informs nothing. Same reason the apical offset is refused there.
    off_flt <- prior_flt[vapply(prior_flt, function(f)
      identical(.canon_switch(fleet_control$Fleet_type[f], fleet_map), "Off"), logical(1))]
    if (length(off_flt) > 0L) {
      stop(sprintf(paste0(
        "selectivity prior on fleet(s) %s with Fleet_type = \"Off\": that fleet's ",
        "selectivity parameters are not estimated, so the prior would add a ",
        "constant to the objective and constrain nothing."),
        paste(fleet_control$Fleet_name[off_flt], collapse = ", ")), call. = FALSE)
    }

    # (c) Fleets sharing a Selectivity_index estimate one parameter block, so a
    # prior on a follower penalizes that block once per sharing fleet.
    mir_lead <- vapply(prior_flt, function(f)
      .shared_block_lead(list(fleet_control = fleet_control), f, "sel"), integer(1))
    if (any(!is.na(mir_lead))) {
      # Name each follower with its own lead: two followers of one lead would
      # otherwise print that lead twice, in a list nothing pairs to the first.
      .bad <- !is.na(mir_lead)
      stop(sprintf(paste0(
        "selectivity prior on %s, which share a Selectivity_index with the lead ",
        "fleet named and take its selectivity block: the shared block would be ",
        "penalized once per sharing fleet. Place the prior on the lead fleet ",
        "instead."),
        paste(sprintf("'%s' (lead: '%s')",
                      fleet_control$Fleet_name[prior_flt[.bad]],
                      fleet_control$Fleet_name[mir_lead[.bad]]),
              collapse = ", ")), call. = FALSE)
    }

    # (c) A prior on a limb the fleet's own curve never uses. Logistic reads only
    # the ascending slots, DescendingLogistic only the descending ones; the other
    # pair stays at its build default and never enters selectivity-at-age. The
    # prior would still be added to the objective -- a constant that shifts the
    # reported likelihood and moves with an unrelated default, while doing
    # nothing to the fit. Silently accepting it is how a reconciliation against
    # another model picks up an unexplained offset.
    used <- list(Logistic           = c("slp_asc", "inf_asc"),
                 DescendingLogistic = c("slp_desc", "inf_desc"))
    for (form in names(used)) {
      f_rows <- prior_rows[vapply(prior_rows$fleet, function(f)
        as.character(fleet_control$Selectivity[row_flt(f)]) == form,
        logical(1)), , drop = FALSE]
      unused <- f_rows[!f_rows$param %in% used[[form]], , drop = FALSE]
      if (nrow(unused) > 0L) {
        stop(sprintf(paste0(
          "selectivity prior on `%s` for %s fleet(s) %s: that %s curve does not ",
          "use those parameters, so the prior would add a constant to the ",
          "objective without affecting the fit. Prior %s instead, or drop the ",
          "fleet from this prior's fleet list."),
          paste(unique(unused$param), collapse = "`, `"), form,
          paste(unique(fleet_control$Fleet_name[
            vapply(unused$fleet, row_flt, integer(1))]), collapse = ", "),
          form, paste(used[[form]], collapse = " / ")), call. = FALSE)
      }
    }
  }
  invisible()
}


#' Refuse an `apical` selectivity linkage the model cannot identify
#'
#' The offset scales one sex's curve. `log_F` is shared by the sexes, so only
#' the ratio between them is identified, and normalization within a sex divides
#' the offset straight back out.
#'
#' @param ap the `apical` rows of the pooled linkage table.
#' @param fleet_control the fleet control table, canonical switch strings.
#' @param nsex sexes per species (`data_list$nsex`); NULL skips the sex checks.
#' @param comp_data the composition data; NULL skips the joint-composition check.
#' @return invisibly NULL; errors on an unidentified offset.
#' @keywords internal
#' @noRd
.check_sel_apical_rows <- function(ap, fleet_control, nsex = NULL, comp_data = NULL) {
  fc  <- fleet_control
  refuse <- function(fmt, flts) {
    stop(sprintf(fmt, paste(fc$Fleet_name[flts], collapse = ", ")), call. = FALSE)
  }

  # The offset is one parameter per fleet and sex, and its prior lands on the
  # fleet the row names; a row for every fleet would free one cell per fleet
  # under a single prior. So every row names its fleet.
  if (anyNA(ap$fleet)) stop(
    "apical selectivity linkage names no fleet: the offset is per fleet, so ",
    "write by = ~ fleet + sex with fleet = <Fleet_code>.", call. = FALSE)
  # A natural-scale offset is added to the multiplier, so below -1 it would
  # make selectivity, F and the predicted catch negative.
  if (any(ap$link == "identity")) stop(
    "apical selectivity linkage with link = \"identity\": the offset is a ",
    "multiplier on the curve, so use link = \"log\" (the default).", call. = FALSE)

  for (i in seq_len(nrow(ap))) {
    flts <- as.integer(ap$fleet[i])
    # Nothing is fit to an "Off" fleet, so the offset is a flat direction: a
    # singular Hessian and a failed getsd with nothing naming the cause. An NA
    # Fleet_type is not "Off" -- build_map_selectivity() treats it as estimated.
    off <- flts[which(vapply(flts, function(f)
      identical(.canon_switch(fc$Fleet_type[f], fleet_map), "Off"), logical(1)))]
    if (length(off)) refuse(paste0(
      "apical selectivity linkage on fleet(s) %s with Fleet_type = \"Off\": no ",
      "data are fit to that fleet, so nothing informs the offset."), off)
    fixed <- flts[as.character(fc$Selectivity[flts]) == "Fixed"]
    if (length(fixed)) refuse(paste0(
      "apical selectivity linkage on fleet(s) %s with Selectivity = \"Fixed\": an ",
      "input curve has no estimated height to offset."), fixed)
    # The AR1 forms already carry a free per-sex level in sel_coff and are held
    # in (0, 1); a multiplier on top is confounded with it.
    ar1 <- flts[as.character(fc$Selectivity[flts]) %in% c("2DAR1", "3DAR1")]
    if (length(ar1)) refuse(paste0(
      "apical selectivity linkage on fleet(s) %s with an AR1 selectivity form: ",
      "those forms estimate a per-sex level in sel_coff already."), ar1)

    # Fleets sharing a Selectivity_index estimate ONE block, this offset
    # included, so a linkage on a follower would free nothing.
    lead <- vapply(flts, function(f) .shared_block_lead(list(fleet_control = fc), f, "sel"),
                   integer(1))
    if (any(!is.na(lead))) {
      # `flts` is one fleet here (ap is read row by row), so name it with its
      # lead rather than printing two lists to be paired up by position.
      .bad <- !is.na(lead)
      stop(sprintf(paste0(
        "apical selectivity linkage on %s, which shares a Selectivity_index with ",
        "the lead fleet named and takes its selectivity block. Place the offset ",
        "on the lead fleet instead; the fleets sharing the index inherit it."),
        paste(sprintf("'%s' (lead: '%s')", fc$Fleet_name[flts[.bad]],
                      fc$Fleet_name[lead[.bad]]), collapse = ", ")), call. = FALSE)
    }

    if (!is.null(nsex)) {
      one_sex <- flts[nsex[fc$Species[flts]] == 1]
      if (length(one_sex)) refuse(paste0(
        "apical selectivity linkage on one-sex fleet(s) %s: with one sex the ",
        "offset is the common selectivity level, which log_F already carries, ",
        "so it is not identified."), one_sex)
      if (is.na(ap$sex[i])) refuse(paste0(
        "apical selectivity linkage on fleet(s) %s names no sex, so both sexes ",
        "would carry the offset and only their ratio is identified. Use ",
        "by = ~ fleet + sex with sex = \"male\" (or \"female\"); the other sex ",
        "is the reference."), flts)
    }

    # Normalization within a sex rescales each sex to its own reference, which
    # removes a whole-curve multiplier exactly. Hake and LogisticPM never reach
    # the shared normalizer.
    norm_on <- !is.na(.rce_sel_norm_code(fc$Sel_norm_bin[flts], allow_all = TRUE))
    within  <- !fc$Sel_norm_scope[flts] %in% c("AcrossSexes", sel_norm_scope_map[["AcrossSexes"]])
    shared  <- !as.character(fc$Selectivity[flts]) %in% c("Hake", "LogisticPM")
    cancel  <- flts[norm_on & within & shared]
    if (length(cancel)) refuse(paste0(
      "apical selectivity linkage on fleet(s) %s whose Sel_norm_scope is ",
      "\"WithinSex\": normalizing each sex to its own reference divides the offset ",
      "out. Set Sel_norm_scope = \"AcrossSexes\", or Sel_norm_bin = \"Off\"."), cancel)
  }

  # Both sexes named across rows on one fleet is the same ridge as naming none:
  # the union of named sexes per fleet must be one sex.
  if (!is.null(nsex)) {
    named  <- ap[!is.na(ap$sex), , drop = FALSE]
    by_flt <- split(as.integer(named$sex), as.integer(named$fleet))
    both   <- as.integer(names(by_flt)[vapply(by_flt, function(s) length(unique(s)) > 1L,
                                              logical(1))])
    if (length(both)) refuse(paste0(
      "apical selectivity linkage on fleet(s) %s names both sexes; only their ",
      "ratio is identified, so one sex carries the offset and the other is the ",
      "reference."), both)
  }

  # The sexes' ratio is informed only by joint-sex compositions. With neither
  # those nor a prior the offset is a free parameter in a flat direction: the
  # fit converges and reports a number the data never constrained.
  if (!is.null(comp_data) && !is.null(comp_data$Sex) && !is.null(comp_data$Fleet_code)) {
    # as.character() first: a factor Sex (read.csv with stringsAsFactors) would
    # otherwise compare as its level index rather than as the code it names.
    sex_code <- suppressWarnings(as.integer(as.character(comp_data$Sex)))
    joint <- unique(as.integer(comp_data$Fleet_code[!is.na(sex_code) & sex_code == 3L]))
    free  <- is.na(ap$est_phase) | as.integer(ap$est_phase) != 0L
    base  <- if (is.null(ap$design_col)) rep(TRUE, nrow(ap)) else
      ap$design_col == "(Intercept)"
    none  <- if (is.null(ap$prior_family)) TRUE else
      is.na(ap$prior_family) | ap$prior_family %in% c("none", "")
    blind <- unique(as.integer(ap$fleet[free & base & none &
                                          !as.integer(ap$fleet) %in% joint]))
    if (length(blind)) warning(sprintf(paste0(
      "apical selectivity linkage on fleet(s) %s: that fleet has no joint-sex ",
      "composition rows (comp_data Sex = 3) and the offset carries no prior, so ",
      "nothing informs the sexes' ratio and the estimate is whatever the ",
      "optimizer leaves. Add priors = list(intercept = lognormal(0, 0.5)), or fit ",
      "that fleet's compositions jointly."),
      paste(fc$Fleet_name[blind], collapse = ", ")), call. = FALSE)
  }
  invisible()
}
