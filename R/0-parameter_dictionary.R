#' Parameter dictionary: internal name -> natural-scale name -> plain meaning
#'
#' @description
#' CEATTLE's TMB parameter vector uses transformed, abbreviated names
#' (`log_sel_slp_dev`, `M1_beta`, `M1_rho`, ...). This table is the single
#' place mapping them to the quantity a user actually chose to estimate, so
#' error messages, diagnostics and output summaries share one vocabulary.
#'
#' Columns:
#' \describe{
#'   \item{internal}{Name in the TMB parameter list / `sdreport` output.}
#'   \item{natural}{The quantity on its natural (untransformed) scale, i.e.
#'     what the internal parameter becomes after `exp()` / `plogis()` /
#'     identity. Equal to `internal` where no transform is applied.}
#'   \item{process}{Which part of the model it belongs to. One of
#'     "recruitment", "mortality", "growth", "fishing", "catchability",
#'     "selectivity", "observation", "predation", "linkage", "internal".}
#'   \item{meaning}{One sentence a fisheries scientist can act on.}
#'   \item{dims}{Dimensions, in the model's own notation.}
#' }
#'
#' @keywords internal
#' @noRd
.PAR_INFO <- local({
  r <- function(internal, natural, process, meaning, dims)
    data.frame(internal = internal, natural = natural, process = process,
               meaning = meaning, dims = dims, stringsAsFactors = FALSE)
  rbind(
    # -- internal / bookkeeping -------------------------------------------
    r("dummy", "dummy", "internal",
      "Placeholder parameter; the only free parameter under estimateMode = 4.", "[1]"),
    r("log_pop_scalar", "pop_scalar", "internal",
      "Multiplier on user-supplied numbers-at-age; estimated for estDynamics = 2 under predation, 1 otherwise.", "[nspp]"),

    # -- recruitment -------------------------------------------------------
    r("rec_pars", "R0 / alpha / beta", "recruitment",
      "Stock-recruit parameters: mean recruitment (or R0), and the SRR alpha and beta.", "[nspp, 3]"),
    r("rec_dev", "recruitment deviations", "recruitment",
      "Annual log-scale recruitment deviations.", "[nspp, nyrs]"),
    r("R_log_sd", "sigma_R", "recruitment",
      "Standard deviation of the recruitment deviations; estimated only when random_rec = TRUE.", "[nspp]"),
    r("init_dev", "initial-age deviations", "recruitment",
      "Deviations defining the initial (first-year) age structure; one per age from minage + 1 upward, since age minage in the first year is recruitment.", "[nspp, nages-1]"),

    # -- mortality ---------------------------------------------------------
    r("log_M1", "M1", "mortality",
      "Natural mortality-at-age; residual M in multispecies mode, total M in single-species mode.", "[nspp, nsex, nages]"),
    r("log_M1_dev", "M1 deviations", "mortality",
      "Annual/age deviations on natural mortality; random effects when M1_re > 0.", "[nspp, nsex, nages, nyrs]"),
    r("M1_beta", "M1 covariate effect", "mortality",
      "LEGACY environmental regression coefficients on M1; superseded by beta_linkage.", "[nspp, nsex, n_env]"),
    r("M1_rho", "rho_M1", "mortality",
      "AR1 correlation for the M1 random effects, over age and over year.", "[nspp, nsex, 2]"),
    r("M1_dev_log_sd", "sigma_M1", "mortality",
      "Innovation standard deviation of the M1 random effects; the marginal sd is larger.", "[nspp, nsex]"),

    # -- growth ------------------------------------------------------------
    r("log_growth_pars", "K / L1 / Linf / m", "growth",
      "Mean growth-curve parameters (von Bertalanffy or Richards).", "[nspp, nsex, 4]"),
    r("growth_log_sd", "sigma_L1 / sigma_Linf", "growth",
      "Standard deviation of length-at-age at the minimum and maximum age.", "[nspp, nsex, 2]"),
    r("weight_length_pars", "alpha_wt / beta_wt", "growth",
      "Weight-length relationship coefficients.", "[nspp, 2]"),

    # -- environmental linkage --------------------------------------------
    r("beta_linkage", "linkage coefficients", "linkage",
      "Coefficients of the formula-driven environmental linkages; one per row of the linkage table.", "[n_linkage]"),
    r("beta_linkage_re", "linkage random-effect deviations", "linkage",
      "Random-effect deviation coefficients for `~ (1|group)` / rw() / ar1() linkage terms; entered into the Laplace approximation. Length 0 without a random linkage.", "[n_linkage_re]"),
    r("beta_linkage_re_pen", "penalized linkage deviations", "linkage",
      "Deviation coefficients for linkage terms declared `linkage_spec(integrate = FALSE)`: estimated as penalized fixed effects rather than integrated out, so they carry the same density as `beta_linkage_re` but as a plain penalty and are reported with standard errors. Requires a fixed SD. Length 0 without a penalized random linkage.", "[n_linkage_re_pen]"),
    r("log_sigma_linkage", "linkage RE log-SD", "linkage",
      "Log standard deviation of each random-effect linkage group; the variance of the time-varying deviations. Can be input (init) or given a prior. Length 0 without a random linkage.", "[n_re_group]"),
    r("trans_rho_linkage", "linkage RE correlation", "linkage",
      "Transformed autocorrelation for each ar1() linkage group. Length 0 unless a correlated (ar1) random linkage is used.", "[n_re_ar1_group]"),
    r("beta_linkage_obs", "linkage QAR1 effect size", "linkage",
      "Rogers et al. (2024) QAR1 effect size: one per observed ar1 linkage group (`observe = `), scaling the latent deviate into the linked parameter. Length 0 unless a state-space covariate is used.", "[n_re_obs_group]"),
    r("log_obs_sd_linkage", "linkage QAR1 observation log-SD", "linkage",
      "Rogers et al. (2024) QAR1 observation log-SD: one per observed ar1 linkage group (`observe = `), the SD of the covariate around the latent deviate. Fixed by default; estimated with `linkage_spec(obs_sd_est = TRUE)`. Length 0 unless a state-space covariate is used.", "[n_re_obs_group]"),

    # -- fishing mortality -------------------------------------------------
    r("log_F", "F", "fishing",
      "Annual fleet-specific fishing mortality.", "[n_flt, nyrs_hind]"),
    r("log_Flimit", "Flimit", "fishing",
      "Limit fishing mortality (FOFL proxy) used in projections.", "[nspp]"),
    r("log_Ftarget", "Ftarget", "fishing",
      "Target fishing mortality (max FABC proxy) used in projections.", "[nspp]"),
    r("log_Finit", "Finit", "fishing",
      "Fishing mortality generating a non-equilibrium initial age structure.", "[nspp]"),
    r("proj_F_prop", "proj_F_prop", "fishing",
      "Fixed proportion of a species' F allocated to each fishery in projections.", "[n_flt]"),

    # -- catchability ------------------------------------------------------
    r("index_log_q", "q", "catchability",
      "Survey/index catchability.", "[n_flt]"),
    r("index_q_beta", "q covariate effect", "catchability",
      "LEGACY environmental regression coefficients on q; superseded by beta_linkage.", "[n_flt, n_env]"),
    r("index_q_dev", "q deviations", "catchability",
      "Annual deviations on catchability when q is time-varying.", "[n_flt, nyrs_hind]"),
    r("index_q_log_sd", "sigma_q_prior", "catchability",
      "Standard deviation of the prior on catchability (used when Estimate_q = 2).", "[n_flt]"),
    r("index_q_dev_log_sd", "sigma_q_dev", "catchability",
      "Standard deviation of the time-varying catchability deviations.", "[n_flt]"),

    # -- selectivity -------------------------------------------------------
    r("log_sel_slp", "selectivity slope", "selectivity",
      "Logistic-family selectivity slope; row 1 ascending, row 2 descending.", "[2, n_sel, nsex]"),
    r("sel_inf", "selectivity inflection", "selectivity",
      "Logistic-family age/length at 50% selection; row 1 ascending, row 2 descending.", "[2, n_sel, nsex]"),
    r("log_sel_apical", "selectivity apical height", "selectivity",
      "Log multiplier on one sex's whole curve, applied after the form and before normalization; 0 = no offset. Estimated only through a selectivity linkage on `apical`.", "[n_sel, nsex]"),
    r("log_sel_slp_dev", "slope deviations", "selectivity",
      "Annual deviations on the selectivity slope.", "[2, n_sel, nsex, nyrs_hind]"),
    r("sel_inf_dev", "inflection deviations", "selectivity",
      "Annual deviations on the selectivity inflection point.", "[2, n_sel, nsex, nyrs_hind]"),
    r("sel_coff", "selectivity coefficients", "selectivity",
      "Non-parametric selectivity-at-bin coefficients.", "[n_sel, nsex, n_sel_bins]"),
    r("sel_coff_dev", "coefficient deviations", "selectivity",
      "Annual deviations on the non-parametric selectivity coefficients.", "[n_sel, nsex, n_sel_bins, nyrs_hind]"),
    r("sel_dev_log_sd", "sigma_sel", "selectivity",
      "Standard deviation of the time-varying selectivity deviations.", "[n_sel]"),
    r("sel_curve_pen", "selectivity penalty weights / AR1 correlations", "selectivity",
      "Meaning depends on the fleet's Selectivity. For the non-parametric forms these are weights on the shape and curvature penalties, supplied via fleet_control and not estimated. For '2DAR1' and '3DAR1' the same slots hold estimated correlations: slot 1 across selectivity BINS, slot 2 across YEARS, and for 3DAR1 slot 3 across COHORTS.", "[n_sel, 3]"),

    # -- observation error / data weighting -------------------------------
    r("index_log_sd", "sigma_index", "observation",
      "Observation error standard deviation of the survey/index series.", "[n_flt]"),
    r("catch_log_sd", "sigma_catch", "observation",
      "Observation error standard deviation of the fishery catch series.", "[n_flt]"),
    r("comp_weights", "comp weight / DM theta", "observation",
      "Composition data weight under a multinomial likelihood; log Dirichlet-multinomial theta under DirichletMultinomial.", "[n_flt]"),
    r("caal_weights", "CAAL weight / DM theta", "observation",
      "As comp_weights, for conditional age-at-length data.", "[n_flt]"),
    r("diet_comp_weights", "diet weight / DM theta", "observation",
      "As comp_weights, for diet composition data.", "[nspp]"),

    # -- predation ---------------------------------------------------------
    r("log_gam_a", "gamma a", "predation",
      "Shape parameter of the parametric predator size-preference function.", "[nspp]"),
    r("log_gam_b", "gamma b", "predation",
      "Scale parameter of the parametric predator size-preference function.", "[nspp]"),
    r("log_phi", "prey preference", "predation",
      "Predator-prey preference coefficients (multinomial-logit transformed to vulnerability).", "[nspp, nspp]")
  )
})


#' Describe a CEATTLE parameter in plain language
#'
#' Looks a parameter up in the dictionary. Used by error messages and
#' diagnostics so the same vocabulary is used everywhere.
#'
#' @param internal character vector of internal parameter names. `NULL`
#'   (default) returns the whole table.
#' @return A `data.frame` slice of the dictionary. Unknown names return a row
#'   with `NA` in the descriptive columns rather than erroring.
#' @keywords internal
#' @noRd
.par_info <- function(internal = NULL) {
  if (is.null(internal)) return(.PAR_INFO)
  idx <- match(internal, .PAR_INFO$internal)
  out <- .PAR_INFO[idx, , drop = FALSE]
  out$internal <- internal
  rownames(out) <- NULL
  out
}


#' One-line human description of a parameter, for use inside messages
#'
#' @param internal character scalar internal parameter name.
#' @return A character scalar, e.g.
#'   `"sel_inf (selectivity inflection): ..."`. Falls back to the bare name
#'   when the parameter is not in the dictionary.
#' @keywords internal
#' @noRd
.par_description <- function(internal) {
  info <- .par_info(internal)
  if (nrow(info) != 1L || is.na(info$meaning)) return(internal)
  sprintf("%s (%s): %s", info$internal, info$natural, info$meaning)
}


#' Look up what a CEATTLE parameter is
#'
#' The TMB parameter vector uses transformed, abbreviated names,
#' `log_M1`, `R_log_sd`, `index_log_q`. This returns the table mapping each one
#' to the quantity it represents on its natural scale, which process it belongs
#' to, what it means, and its dimensions.
#'
#' @param internal Character vector of internal parameter names, as they appear
#'   in `names(fitted$estimated_params)` or in an `sdreport()`. Default `NULL`
#'   returns every parameter.
#' @param process Character vector restricting the table to one or more parts of
#'   the model: `"recruitment"`, `"mortality"`, `"growth"`, `"fishing"`,
#'   `"catchability"`, `"selectivity"`, `"observation"`, `"predation"`,
#'   `"linkage"` or `"internal"`.
#'
#' @return A data frame with one row per parameter and columns `internal`,
#'   `natural`, `process`, `meaning` and `dims`.
#'
#' @details
#' A name beginning `log_`, or ending `_log_sd`, is estimated on the log scale,
#' so `exp()` turns it into the quantity named in the `natural` column. Most
#' fitted values are already back-transformed in `fitted$quantities`, so this is
#' mainly needed when reading `estimated_params` or an `sdreport()` directly.
#'
#' Dimensions use the model's own notation: `nspp` species, `nsex` sexes,
#' `nages` age bins, `nyrs` years including the projection, `nyrs_hind` hindcast
#' years only, `n_flt` fleets and `n_sel` selectivity blocks. `n_sel` equals
#' `n_flt`: fleets sharing a `Selectivity_index` share one block through the
#' map, not through a narrower array.
#'
#' @examples
#' # The whole table
#' head(parameter_dictionary())
#'
#' # What is R_log_sd?
#' parameter_dictionary("R_log_sd")
#'
#' # Everything governing selectivity
#' parameter_dictionary(process = "selectivity")
#'
#' # Search the meanings by keyword
#' dict <- parameter_dictionary()
#' dict[grep("catchability", dict$meaning, ignore.case = TRUE), ]
#' @seealso [fit_mod()] for the fitted object these names come from.
#' @export
parameter_dictionary <- function(internal = NULL, process = NULL) {
  out <- .PAR_INFO

  if (!is.null(process)) {
    known <- unique(.PAR_INFO$process)
    unknown <- setdiff(process, known)
    if (length(unknown)) {
      stop("Unknown process: ", paste(unknown, collapse = ", "),
           ". Use one of: ", paste(known, collapse = ", "))
    }
    out <- out[out$process %in% process, , drop = FALSE]
  }

  if (!is.null(internal)) {
    # Warn rather than stop, so passing every name from a fit still works when
    # that fit carries parameters from a branch this dictionary predates.
    unknown <- setdiff(internal, .PAR_INFO$internal)
    if (length(unknown)) {
      warning("Not in the parameter dictionary: ",
              paste(unknown, collapse = ", "), call. = FALSE)
    }
    out <- out[out$internal %in% internal, , drop = FALSE]
  }

  rownames(out) <- NULL
  out
}
