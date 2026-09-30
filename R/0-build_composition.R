#' @keywords internal
#' @noRd
COMP_LINKAGE_PARAMS <- c("theta_comp", "theta_caal", "theta_diet")


#' @keywords internal
#' @noRd
.validate_comp_linkages <- function(linkages) {
  .validate_process_linkages(linkages, COMP_LINKAGE_PARAMS, "comp")
}


#' Composition-weighting specification
#'
#' @description
#' Holds **priors** on the Dirichlet-multinomial composition-weighting
#' overdispersion. The DM weight (the "theta" that scales the effective sample
#' size) is otherwise an unpenalized free parameter; a linkage lets you put a
#' prior on it through the same grammar as every other parameter. The three
#' parameters target the three DM likelihoods:
#' \describe{
#'   \item{`theta_comp`}{age / length composition weight (`comp_weights`), per
#'     fleet.}
#'   \item{`theta_caal`}{conditional-age-at-length weight (`caal_weights`), per
#'     fleet.}
#'   \item{`theta_diet`}{diet composition weight (`diet_comp_weights`), per
#'     predator.}
#' }
#'
#' This is a **prior-only** process: the DM weight is a scalar, not a
#' year-varying quantity, so only intercept formulas (`~ 1`) with a `priors`
#' entry are meaningful. The prior is placed on the natural-scale DM weight
#' `theta = exp(weight)`, so a `gamma()` prior reads naturally. A linkage on a
#' fleet whose `Comp_distribution` (or `CAAL_distribution`) is not
#' `"DirichletMultinomial"` errors, since the weight has no effect there.
#'
#' @param linkages Optional named list of [linkage_spec()] objects keyed by
#'   `theta_comp` / `theta_caal` (per fleet by default, `by = ~ fleet`) or
#'   `theta_diet` (per predator by default, `by = ~ species`).
#'
#' @return A list of composition-weighting settings for [fit_mod()].
#'
#' @examples
#' \donttest{
#' # A weak gamma prior on the DM overdispersion of every fleet's age comps
#' build_composition(linkages = list(
#'   theta_comp = linkage_spec(~ 1, by = ~ fleet,
#'                             priors = list(`(Intercept)` = gamma(2, 0.5)))))
#' }
#'
#' @export
build_composition <- function(linkages = NULL) {
  list(linkages = .validate_comp_linkages(linkages))
}


#' Reject composition-weighting (comp) linkages that cannot take effect
#'
#' @description
#' `comp` linkages are prior-only (the Dirichlet-multinomial weight is a scalar,
#' not a year-varying quantity), so only intercept formulas are allowed; a
#' covariate slope would be estimated to no effect. And the weight is only a
#' free parameter under a `"DirichletMultinomial"` likelihood, so a prior on a
#' non-DM fleet/predator would target a fixed parameter. Both error up front.
#'
#' @param linkage_table pooled linkage table (may be NULL / empty).
#' @param data_list the data list (for `fleet_control` and `Diet_distribution`).
#' @return invisibly NULL; errors on an unsupported comp linkage.
#' @keywords internal
#' @noRd
.check_comp_linkage_support <- function(linkage_table, data_list) {
  if (is.null(linkage_table) || nrow(linkage_table) == 0L) return(invisible())
  cmp <- linkage_table[linkage_table$process == "comp", , drop = FALSE]
  if (nrow(cmp) == 0L) return(invisible())
  fc <- data_list$fleet_control

  # (a) prior-only: only intercept rows (theta is a scalar; no accumulator).
  if (any(cmp$design_col != "(Intercept)")) {
    stop("composition-weighting (comp) linkages are prior-only: use an ",
         "intercept formula `~ 1` with a `priors` entry (the DM weight is a ",
         "scalar, not year-varying), not a covariate slope.", call. = FALSE)
  }

  # (a2) `est_phase` still does not apply: the intercept coefficient it would
  # control is mapped out at 0 for comp, so phasing it does nothing. `init` and
  # `bounds` DO apply -- they re-target the DM weight itself, the same contract
  # every other process gives the intercept (see build_params() "Push
  # (Intercept) inits to the base parameter").
  if (any(cmp$est_phase != 1L)) {
    stop("composition-weighting (comp) linkages are prior-only: `est_phase` on ",
         "the spec does not fix or phase the DM weight (the intercept ",
         "coefficient it controls is mapped out at 0); fix the weight via ",
         "`map` instead.", call. = FALSE)
  }

  # (b) each row must name the ONE weight it targets. The cpp prior loop runs
  # once per linkage row and addresses a single weight -- comp/caal weights are
  # fleet-indexed, diet weights predator(species)-indexed -- so a row must carry
  # a concrete fleet (theta_comp / theta_caal) or species (theta_diet). There is
  # NO "all fleets/species" broadcast: an unstratified sentinel (NA) collapses
  # to the first index in the cpp, silently attaching the prior to fleet 1 /
  # predator 1 (and dropping the rest). `theta_comp` / `theta_caal` therefore
  # need `by = ~ fleet` (comp weights are not species-indexed); `theta_diet`
  # needs `by = ~ species`. Reject the wrong / missing stratum up front.
  dm_str <- function(v) as.character(v) == "DirichletMultinomial"
  for (r in seq_len(nrow(cmp))) {
    prm <- cmp$param[r]
    fl  <- cmp$fleet[r]
    sp  <- cmp$species[r]
    if (prm %in% c("theta_comp", "theta_caal")) {
      if (is.na(fl)) {
        stop(sprintf(paste0(
          "%s is a fleet-indexed DM weight: stratify the linkage by fleet ",
          "(`by = ~ fleet`, naming fleets via `fleet = `). An unstratified ",
          "spec would attach the prior to a single fleet only."), prm),
          call. = FALSE)
      }
      col <- if (prm == "theta_comp") "Comp_distribution" else "CAAL_distribution"
      if (!dm_str(fc[[col]][fl])) {
        stop(sprintf(paste0(
          "%s linkage on fleet %s whose %s is not 'DirichletMultinomial'; ",
          "the DM weight is not estimated there, so the prior has no effect."),
          prm, fc$Fleet_name[fl], col), call. = FALSE)
      }
    } else if (prm == "theta_diet") {
      if (is.na(sp)) {
        stop(paste0(
          "theta_diet is a predator-indexed DM weight: stratify the linkage ",
          "by species (`by = ~ species`). An unstratified spec would attach ",
          "the prior to a single predator only."), call. = FALSE)
      }
      if (as.integer(data_list$Diet_distribution[sp]) != 1L) {
        stop(sprintf(paste0(
          "theta_diet linkage on predator %s whose Diet_distribution is not ",
          "DirichletMultinomial; the DM weight is not estimated there."), sp),
          call. = FALSE)
      }
    }
  }
  invisible()
}


#' Drop comp priors whose DM weight is fixed in this configuration
#'
#' @description
#' A `comp` linkage row is prior-only, so when the map fixes the DM weight it
#' targets the prior is a constant: it shifts the reported `jnll` without moving
#' an estimate, and makes likelihoods non-comparable across configurations.
#' Such rows are set to `prior_family = "none"`, which the model skips.
#'
#' `.check_comp_linkage_support()` rejects a prior that can never apply to the
#' data at hand (a non-DM `Comp_distribution` / `Diet_distribution`). Whether a
#' weight is estimated *in a given fit* also depends on `msmMode`, `suitMode`,
#' and the fleet setup, and one `compFun` is routinely shared across the
#' single-species and multispecies fits of a stock, so those are reported and
#' ignored rather than rejected.
#'
#' Inertness is read off the finished `map` so it stays in step with
#' [build_map()] and honors a user-supplied `map`. Rows are kept, not dropped:
#' `beta_linkage` is dimensioned by `nrow(linkage_table)`, so dropping them
#' would break `inits` reuse between fits sharing a `compFun`.
#'
#' @param linkage_table pooled linkage table (may be NULL / empty).
#' @param map the map object from [build_map()] (uses `$mapList`).
#' @param verbose integer; 0 silences the message.
#' @return the linkage table, with inert comp priors neutralized.
#' @keywords internal
#' @noRd
.neutralize_inert_comp_priors <- function(linkage_table, map, verbose = 1) {
  if (is.null(linkage_table) || nrow(linkage_table) == 0L) return(linkage_table)
  if (is.null(map) || is.null(map$mapList)) return(linkage_table)

  # comp param -> the map slot holding that DM weight, and how it is indexed.
  slots <- list(theta_comp = "comp_weights",
                theta_caal = "caal_weights",
                theta_diet = "diet_comp_weights")

  inert <- rep(FALSE, nrow(linkage_table))
  for (i in which(linkage_table$process == "comp" &
                  linkage_table$prior_family != "none")) {
    prm  <- linkage_table$param[i]
    slot <- slots[[prm]]
    if (is.null(slot)) next
    m <- map$mapList[[slot]]
    if (is.null(m)) next
    # theta_comp / theta_caal are fleet-indexed, theta_diet species-indexed.
    idx <- if (prm == "theta_diet") linkage_table$species[i] else linkage_table$fleet[i]
    if (is.na(idx) || idx < 1L || idx > length(m)) next
    inert[i] <- is.na(m[[idx]])
  }

  if (any(inert)) {
    if (verbose > 0) {
      message(sprintf(
        paste0("Ignoring %d composition-weighting prior(s) on a DM weight that ",
               "is not estimated in this configuration (%s). A prior on a fixed ",
               "parameter only adds a constant to the objective."),
        sum(inert),
        paste(sprintf("%s[%s]", linkage_table$param[inert],
                      ifelse(linkage_table$param[inert] == "theta_diet",
                             linkage_table$species[inert],
                             linkage_table$fleet[inert])),
              collapse = ", ")))
    }
    linkage_table$prior_family[inert] <- "none"
    linkage_table$prior_p1[inert]     <- NA_real_
    linkage_table$prior_p2[inert]     <- NA_real_
  }
  linkage_table
}
