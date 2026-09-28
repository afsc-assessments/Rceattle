#' Quantity dictionary: report name -> meaning, units, dimensions
#'
#' @description
#' `fit$quantities` holds every derived quantity the TMB model reports, under
#' the model's own abbreviated names (`ssb`, `F_spp`, `NByageF`, ...). This
#' table is the single place saying what each one is, the units it is in, how it
#' is shaped, whether it holds a standard error, and what the same quantity is
#' called in the NOAA standardized assessment output.
#'
#' Columns:
#' \describe{
#'   \item{quantity}{Name in `fit$quantities`.}
#'   \item{process}{Which part of the model it belongs to. One of
#'     "population", "recruitment", "mortality", "fishing", "reference_points",
#'     "catchability", "selectivity", "growth", "composition", "predation",
#'     "likelihood", "linkage", "internal".}
#'   \item{meaning}{One sentence a fisheries scientist can act on.}
#'   \item{units}{The units the value is in, or "unitless" / "proportion".}
#'   \item{dims}{Dimensions, in the model's own notation.}
#'   \item{se}{Whether [TMB::sdreport()] gives a standard error for it, i.e.
#'     whether the template `ADREPORT`s it. `FALSE` means `fit$sdrep` holds
#'     nothing for this quantity and any interval must come from elsewhere.}
#'   \item{standard_label}{The `label` this quantity takes in the NOAA
#'     standardized assessment output consumed by `stockplotr` and `asar`, or
#'     `NA` where that standard has no equivalent.}
#' }
#'
#' Units follow the model throughout: numbers-at-age are **thousands of fish**
#' and weight-at-age is **kg**, so their product is **mt**; see the
#' `Observation_units` column in `R/0-column_schema.R`.
#'
#' @keywords internal
#' @noRd
.QUANT_INFO <- local({
  r <- function(quantity, process, meaning, units, dims, se = FALSE,
                standard_label = NA_character_)
    data.frame(quantity = quantity, process = process, meaning = meaning,
               units = units, dims = dims, se = se,
               standard_label = standard_label, stringsAsFactors = FALSE)
  rbind(
    # -- population --------------------------------------------------------
    r("biomass", "population",
      "Total stock biomass, summed over sexes and ages.",
      "mt", "[nspp, nyrs]", TRUE, "biomass"),
    r("log_biomass", "population",
      "Total stock biomass on the log scale; its standard error is the CV of biomass, which is the form an ABC buffer wants.",
      "log mt", "[nspp, nyrs]", TRUE, NA_character_),
    r("ssb", "population",
      "Female spawning-stock biomass at the spawning month.",
      "mt", "[nspp, nyrs]", TRUE, "spawning_biomass"),
    r("log_ssb", "population",
      "Female spawning-stock biomass on the log scale; its standard error is the CV of SSB.",
      "log mt", "[nspp, nyrs]", TRUE, NA_character_),
    r("exploitable_biomass", "population",
      "Biomass selected by the fisheries, summed over fishery fleets; exactly zero for a survey-only species or when proj_F_prop is zero.",
      "mt", "[nspp, nyrs]", TRUE, NA_character_),
    r("biomass_depletion", "population",
      "Total biomass relative to unfished biomass, biomass / B0.",
      "proportion", "[nspp, nyrs]", TRUE, NA_character_),
    r("ssb_depletion", "population",
      "Female spawning biomass relative to unfished, ssb / SB0; the quantity a Tier 3 harvest control rule compares against B40%. NA for a species with input numbers-at-age in single-species mode unless DynamicHCR = TRUE, where it is those numbers relative to themselves (1 at spawn_month = 0).",
      "proportion", "[nspp, nyrs]", TRUE, "relative_spawning_biomass"),
    r("N_at_age", "population",
      "Numbers at age at the start of the year.",
      "thousands of fish", "[nspp, nsex, nages, nyrs]", FALSE, "abundance"),
    r("avgN_at_age", "population",
      "Average numbers at age over the year, the abundance the catch and consumption equations use.",
      "thousands of fish", "[nspp, nsex, nages, nyrs]", FALSE, "abundance_midyear"),
    r("biomass_at_age", "population",
      "Biomass at age, numbers at age times weight at age.",
      "mt", "[nspp, nsex, nages, nyrs]", FALSE, "biomass_at_age"),

    # -- recruitment -------------------------------------------------------
    r("R", "recruitment",
      "Recruitment: numbers entering at the youngest age bin. For a species with input numbers-at-age (estDynamics > 0), the input recruits, with no standard error.",
      "thousands of fish", "[nspp, nyrs]", TRUE, "recruitment"),
    r("log_R", "recruitment",
      "Recruitment on the log scale; its standard error is the CV of recruitment.",
      "log thousands of fish", "[nspp, nyrs]", TRUE, NA_character_),
    r("R_hat", "recruitment",
      "Recruitment predicted by the stock-recruit curve, before the annual deviation is applied. Under msmMode > 0 the first year is R_init.",
      "thousands of fish", "[nspp, nyrs]", FALSE, NA_character_),
    r("R0", "recruitment",
      "Equilibrium recruitment at F = 0. Under msmMode > 0 with a hindcast curve, the estimated initial recruitment level.",
      "thousands of fish", "[nspp, nyrs]", FALSE, NA_character_),
    r("R_init", "recruitment",
      "Equilibrium recruitment at F = Finit, the non-equilibrium initial age structure. Under msmMode > 0 with a hindcast curve, the estimated initial recruitment level.",
      "thousands of fish", "[nspp]", FALSE, NA_character_),
    r("avg_R", "recruitment",
      "Mean recruitment over the hindcast, used for the reference-point age structures.",
      "thousands of fish", "[nspp]", FALSE, NA_character_),
    r("R_sd", "recruitment",
      "Standard deviation of the recruitment deviations, sigma_R, on the natural scale.",
      "log scale sd", "[nspp]", TRUE, NA_character_),
    r("steepness", "recruitment",
      "Expected fraction of R0 produced at 20% of unfished spawning biomass. Zero under msmMode > 0 with a stock-recruit curve, where spawning biomass per recruit is undefined.",
      "proportion", "[nspp, nyrs]", FALSE, NA_character_),

    # -- mortality ---------------------------------------------------------
    r("M_at_age", "mortality",
      "Total natural mortality at age; residual M1 plus predation M2 in multispecies mode.",
      "yr^-1", "[nspp, nsex, nages, nyrs]", FALSE, "natural_mortality"),
    r("M1_at_age", "mortality",
      "Residual natural mortality at age in multispecies mode; total natural mortality in single-species mode.",
      "yr^-1", "[nspp, nsex, nages, nyrs]", FALSE, NA_character_),
    r("M2_at_age", "mortality",
      "Predation mortality at age; zero under msmMode = 0.",
      "yr^-1", "[nspp, nsex, nages, nyrs]", FALSE, NA_character_),
    r("M2_prop", "mortality",
      "Share of a prey age's predation mortality attributable to each predator age.",
      "proportion", "[nspp*nsex, nspp*nsex, nages, nages, nyrs]", FALSE, NA_character_),
    r("Z_at_age", "mortality",
      "Total mortality at age, natural plus fishing.",
      "yr^-1", "[nspp, nsex, nages, nyrs]", FALSE, NA_character_),
    r("mort_sum", "mortality",
      "Cumulative mortality used to build the initial age structure.",
      "yr^-1", "[nspp, nages]", FALSE, NA_character_),

    # -- fishing -----------------------------------------------------------
    r("F_spp", "fishing",
      "Fully selected fishing mortality by species, summed over its fishery fleets.",
      "yr^-1", "[nspp, nyrs]", FALSE, "fishing_mortality"),
    r("F_flt", "fishing",
      "Fully selected fishing mortality by fleet.",
      "yr^-1", "[n_flt, nyrs]", FALSE, NA_character_),
    r("F_at_age", "fishing",
      "Fishing mortality at age for each species, summed over its fishery fleets.",
      "yr^-1", "[nspp, nsex, nages, nyrs]", FALSE, NA_character_),
    r("F_flt_age", "fishing",
      "Fishing mortality at age for each fishery fleet.",
      "yr^-1", "[n_flt, nsex, nages, nyrs]", FALSE, NA_character_),
    r("proj_F", "fishing",
      "Fishing mortality the harvest control rule set for each projection year.",
      "yr^-1", "[nspp, nyrs]", FALSE, NA_character_),
    r("catch_hat", "fishing",
      "Predicted fishery catch, one value per row of catch_data; weight or numbers per the fleet's Observation_units.",
      "mt or thousands of fish", "[nrow(catch_data)]", FALSE, "catch_predicted"),
    r("max_catch_hat", "fishing",
      "Predicted exploitable biomass or numbers available to the fleet, the ceiling the catch equation can take.",
      "mt or thousands of fish", "[nrow(catch_data)]", FALSE, NA_character_),
    r("catch_sd", "fishing",
      "Observation-error sd the catch likelihood used for each row. Catch is lognormal throughout, so this is a log-scale sd; it is not a log of anything.",
      "log scale sd", "[nrow(catch_data)]", FALSE, NA_character_),
    r("log_catch_sd", "fishing",
      "Deprecated spelling of catch_sd, kept one release; identical values. Use catch_sd.",
      "log scale sd", "[nrow(catch_data)]", FALSE, NA_character_),
    r("catch_analytical_sd", "fishing",
      "Concentrated (analytical) catch sd used when Estimate_catch_sd = 2; a log-scale sd, not a log.",
      "log scale sd", "[n_flt]", FALSE, NA_character_),

    # -- reference points --------------------------------------------------
    # The SPR block runs only under msmMode = 0, so every per-recruit quantity
    # below is exactly zero on a multispecies fit -- M there carries predation
    # M2, which scales with predator abundance, so spawning output per recruit
    # is not a property of the prey stock alone.
    r("Flimit", "reference_points",
      "Limit fishing mortality, the FOFL proxy (F35% under Tier 3) used in projections; NA for a species with input numbers-at-age.",
      "yr^-1", "[nspp]", FALSE, NA_character_),
    r("Ftarget", "reference_points",
      "Target fishing mortality, the maximum FABC proxy (F40% under Tier 3) used in projections; NA for a species with input numbers-at-age.",
      "yr^-1", "[nspp]", FALSE, NA_character_),
    r("B0", "reference_points",
      "Total biomass at F = 0, carrying the estimated stock-recruit relationship.",
      "mt", "[nspp, nyrs]", FALSE, NA_character_),
    r("SB0", "reference_points",
      "Female spawning biomass at F = 0; the B100% the Tier 3 proxies are taken from.",
      "mt", "[nspp, nyrs]", FALSE, "spawning_biomass_zero"),
    r("SBF", "reference_points",
      "Female spawning biomass at F = Ftarget; NA for a species with input numbers-at-age.",
      "mt", "[nspp, nyrs]", FALSE, NA_character_),
    r("DynamicB0", "reference_points",
      "Total biomass under the realized recruitment history with F set to zero. A species with input numbers-at-age keeps them.",
      "mt", "[nspp, nyrs]", FALSE, NA_character_),
    r("DynamicSB0", "reference_points",
      "Female spawning biomass under the realized recruitment history with F set to zero. A species with input numbers-at-age keeps them.",
      "mt", "[nspp, nyrs]", FALSE, NA_character_),
    r("DynamicSBF", "reference_points",
      "Female spawning biomass under the realized recruitment history with F set to Ftarget; NA for a species with input numbers-at-age.",
      "mt", "[nspp, nyrs]", FALSE, NA_character_),
    r("SPR0", "reference_points",
      "Spawning biomass per recruit at F = 0. Zero under msmMode > 0.",
      "kg per recruit", "[nspp]", FALSE, NA_character_),
    r("SPRlimit", "reference_points",
      "Spawning biomass per recruit at F = Flimit. Zero under msmMode > 0; NA for a species with input numbers-at-age.",
      "kg per recruit", "[nspp]", FALSE, NA_character_),
    r("SPRtarget", "reference_points",
      "Spawning biomass per recruit at F = Ftarget. Zero under msmMode > 0; NA for a species with input numbers-at-age.",
      "kg per recruit", "[nspp]", FALSE, NA_character_),
    r("SPRFinit", "reference_points",
      "Spawning biomass per recruit at F = Finit. Zero under msmMode > 0.",
      "kg per recruit", "[nspp]", FALSE, NA_character_),
    r("NbyageSPR", "reference_points",
      "Survivors per recruit behind the four SPR schedules; the first dimension is F = 0, Flimit, Ftarget, Finit. Zero under msmMode > 0.",
      "fish per recruit", "[4, nspp, nages]", FALSE, NA_character_),
    r("NByage0", "reference_points",
      "Numbers at age at mean recruitment and F = 0, the age structure behind B0 and SB0.",
      "thousands of fish", "[nspp, nsex, nages, nyrs]", FALSE, NA_character_),
    r("NByageF", "reference_points",
      "Numbers at age at mean recruitment and F = Ftarget, the age structure behind SBF.",
      "thousands of fish", "[nspp, nsex, nages, nyrs]", FALSE, NA_character_),

    # -- catchability / index ----------------------------------------------
    r("index_q", "catchability",
      "Survey catchability by fleet and hindcast year; time-varying where the fleet estimates q deviations.",
      "unitless", "[n_flt, nyrs_hind]", FALSE, NA_character_),
    r("index_hat", "catchability",
      "Predicted survey index, one value per row of index_data; weight or numbers per the fleet's Observation_units.",
      "mt or thousands of fish", "[nrow(index_data)]", FALSE, "index_predicted"),
    r("log_index_hat", "catchability",
      "Predicted survey index on the log scale.",
      "log mt or log thousands of fish", "[nrow(index_data)]", TRUE, NA_character_),
    r("index_sd", "catchability",
      "Observation-error sd the index likelihood used for each row. It is NOT a log: a log-scale sd for a Lognormal fleet, and an ABSOLUTE sd in index units for a natural-scale Index_distribution.",
      "log scale sd or index units", "[nrow(index_data)]", FALSE, NA_character_),
    r("log_index_sd", "catchability",
      "Deprecated spelling of index_sd, kept one release; identical values. Use index_sd.",
      "log scale sd or index units", "[nrow(index_data)]", FALSE, NA_character_),

    # -- selectivity -------------------------------------------------------
    r("sel_at_age", "selectivity",
      paste("Selectivity at age by fleet, normalized per Sel_norm_bin: a maximum",
            "of one, one at a nominated bin, or -- for the non-parametric forms,",
            "which are the default -- a mean of one, where values above one are",
            "ordinary. Zero below the fleet's first selected bin. Carries a",
            "standard error",
            "only via log_sel_at_age."),
      "unitless", "[n_flt, nsex, nages, nyrs]", FALSE, NA_character_),
    r("sel_at_length", "selectivity",
      paste("Selectivity at length bin by fleet, normalized as sel_at_age is;",
            "see there."),
      "unitless", "[n_flt, nsex, nlengths, nyrs]", FALSE, NA_character_),
    r("log_sel_at_age", "selectivity",
      paste("Log selectivity at age, reported only under",
            "fit_control(selectivity_se = TRUE), for a delta-method interval on",
            "the curve. One entry per estimated age-based lead fleet, sex, age",
            "at or above the fleet's first selected bin, and hindcast year;",
            "log_sel_at_age_index says which. Use exp(value +/- 1.96 * sd)."),
      "log scale", "[nrow(log_sel_at_age_index)]", TRUE, NA_character_),
    r("log_sel_at_age_index", "selectivity",
      paste("Fleet, sex, absolute age and calendar year of each log_sel_at_age",
            "entry, in that column order. Reported alongside it, and the only way",
            "to place a row, since which cells are reported depends on the model."),
      "unitless", "[nrow(log_sel_at_age_index), 4]", FALSE, NA_character_),

    # -- growth ------------------------------------------------------------
    r("weight_hat", "growth",
      "Weight at age. The first dimension runs two rows per species (biomass weight, spawning weight) then one row per fleet.",
      "kg", "[nspp*2 + n_flt, nsex, nages, nyrs]", FALSE, NA_character_),
    r("length_hat", "growth",
      "Length at age, indexed like weight_hat.",
      "cm", "[nspp*2 + n_flt, nsex, nages, nyrs]", FALSE, NA_character_),
    r("growth_matrix", "growth",
      "Age-length transition matrix, indexed like weight_hat.",
      "proportion", "[nspp*2 + n_flt, nsex, nages, nlengths, nyrs]", FALSE, NA_character_),
    r("growth_parameters", "growth",
      "Realized mean growth parameters by year, in the order log_K, log_L1, log_Linf, log_m.",
      "log scale", "[nspp, nsex, nyrs, 4]", FALSE, NA_character_),

    # -- composition -------------------------------------------------------
    r("comp_obs", "composition",
      "Observed age or length composition as proportions, one row per row of comp_data.",
      "proportion", "[nrow(comp_data), max bins]", FALSE, "composition_observed"),
    r("comp_hat", "composition",
      "Predicted age or length composition as proportions, aligned row-for-row with comp_obs.",
      "proportion", "[nrow(comp_data), max bins]", FALSE, "composition_predicted"),
    r("caal_obs", "composition",
      "Observed conditional age-at-length composition as proportions.",
      "proportion", "[nrow(caal_data), nages]", FALSE, NA_character_),
    r("caal_hat", "composition",
      "Predicted conditional age-at-length composition as proportions.",
      "proportion", "[nrow(caal_data), nages]", FALSE, NA_character_),
    r("age_hat", "composition",
      "Predicted catch at true age, before ageing error is applied.",
      "thousands of fish", "[nrow(comp_data), max age cols]", FALSE, NA_character_),
    r("age_obs_hat", "composition",
      "Predicted catch at observed age, after the ageing-error matrix is applied.",
      "thousands of fish", "[nrow(comp_data), max age cols]", FALSE, NA_character_),
    r("pred_CAAL", "composition",
      "Predicted conditional age-at-length by fleet, before it is matched to the observed rows.",
      "proportion", "[n_flt, nsex, nages, nlengths, nyrs]", FALSE, NA_character_),

    # -- predation ---------------------------------------------------------
    r("consumption_at_age", "predation",
      "Annual ration: biomass consumed per predator at age.",
      "kg per fish per yr", "[nspp, nsex, nages, nyrs]", FALSE, NA_character_),
    r("B_eaten", "predation",
      "Biomass of prey at age eaten by predator at age.",
      "mt", "[nspp*nsex, nspp*nsex, nages, nages, nyrs]", FALSE, NA_character_),
    r("B_eaten_as_prey", "predation",
      "Biomass of a species at age eaten by all predators in the model.",
      "mt", "[nspp, nsex, nages, nyrs]", FALSE, NA_character_),
    r("avail_food", "predation",
      "Biomass of prey available to a predator at age, including other food.",
      "mt", "[nspp, nsex, nages, nyrs]", FALSE, NA_character_),
    r("suitability", "predation",
      "Suitability of prey at age to predator at age.",
      "proportion", "[nspp*nsex, nspp*nsex, nages, nages, nyrs]", FALSE, NA_character_),
    r("suit_other", "predation",
      "Suitability of prey outside the model, the remainder of a predator's diet.",
      "proportion", "[nspp, nsex, nages, nyrs]", FALSE, NA_character_),
    r("stom_div_bio", "predation",
      "Observed stomach proportion divided by prey biomass, the empirical suitability numerator.",
      "mt^-1", "[nspp*nsex, nspp*nsex, nages, nages, nyrs]", FALSE, NA_character_),
    r("diet_obs", "predation",
      "Observed stomach-content proportion by weight, one row per row of diet_data.",
      "proportion", "[nrow(diet_data), 2]", FALSE, NA_character_),
    r("diet_hat", "predation",
      "Predicted stomach-content proportion by weight, aligned with diet_obs.",
      "proportion", "[nrow(diet_data), 2]", FALSE, NA_character_),
    r("vulnerability", "predation",
      "Predator-prey preference coefficients.",
      "proportion", "[nspp, nspp]", FALSE, NA_character_),
    r("vulnerability_other", "predation",
      "Predator preference for food outside the model.",
      "proportion", "[nspp]", FALSE, NA_character_),
    r("gam_a", "predation",
      "Predator size-preference shape parameter for gamma suitability; the mean for a normal of logs.",
      "unitless", "[nspp]", FALSE, NA_character_),
    r("gam_b", "predation",
      "Predator size-preference scale parameter for gamma suitability; the sd for a normal of logs.",
      "unitless", "[nspp]", FALSE, NA_character_),
    r("fT", "predation",
      "Temperature scalar on consumption from the bioenergetics submodel.",
      "unitless", "[nspp, nyrs]", FALSE, NA_character_),

    # -- likelihood --------------------------------------------------------
    r("jnll", "likelihood",
      "Joint negative log-likelihood, the objective the optimizer minimized.",
      "unitless", "[1]", FALSE, NA_character_),
    r("jnll_comp", "likelihood",
      "Weighted negative log-likelihood by component and fleet or species. Rows are named; columns count FLEETS on rows 1-8, SPECIES on rows 9-20 and neither on row 21, so rowSums() pools across different axes.",
      "unitless", "[21, max(n_flt, nspp)]", FALSE, "likelihood"),
    r("unweighted_jnll_comp", "likelihood",
      "As jnll_comp before data weights are applied. Written for only 5 of its 21 rows -- composition, CAAL, stomach and the two linkage rows; every other row is structurally zero here, not small.",
      "unitless", "[21, max(n_flt, nspp)]", FALSE, NA_character_),

    # -- linkage -----------------------------------------------------------
    r("beta_linkage", "linkage",
      "Fitted coefficients of the environmental linkage formulas, one per row of the linkage table.",
      "varies by linked parameter", "[n_linkage]", TRUE, NA_character_),
    r("beta_linkage_re", "linkage",
      "Random-effect deviation coefficients for random linkage terms; length 0 without one.",
      "varies by linked parameter", "[n_linkage_re]", FALSE, NA_character_),
    r("beta_linkage_re_pen", "linkage",
      "Penalized (not integrated) linkage deviation coefficients; length 0 without one.",
      "varies by linked parameter", "[n_linkage_re_pen]", FALSE, NA_character_),
    r("beta_linkage_re_all", "linkage",
      "Every linkage deviation coefficient, integrated and penalized together, in linkage-table order.",
      "varies by linked parameter", "[n_linkage_re + n_linkage_re_pen]", FALSE, NA_character_),
    r("beta_linkage_obs", "linkage",
      "QAR1 effect size scaling the latent deviate into the linked parameter; length 0 without a state-space covariate.",
      "varies by linked parameter", "[n_re_obs_group]", TRUE, NA_character_),
    r("M_linkage_offset", "linkage",
      "Linkage offset applied to natural mortality, on the log scale.",
      "log scale", "[nspp, nsex, nages, nyrs]", FALSE, NA_character_),
    r("M_linkage_offset_nat", "linkage",
      "Linkage offset applied to natural mortality, on the natural scale.",
      "multiplier", "[nspp, nsex, nages, nyrs]", FALSE, NA_character_),
    r("q_linkage_offset", "linkage",
      "Linkage offset applied to catchability, on the log scale.",
      "log scale", "[n_flt, nyrs_hind]", FALSE, NA_character_),
    r("q_linkage_offset_nat", "linkage",
      "Linkage offset applied to catchability, on the natural scale.",
      "multiplier", "[n_flt, nyrs_hind]", FALSE, NA_character_),
    r("recruitment_linkage_offset", "linkage",
      "Linkage offset applied to each recruitment parameter, on the log scale.",
      "log scale", "[nspp, 3, nyrs]", FALSE, NA_character_),
    r("recruitment_linkage_offset_nat", "linkage",
      "Linkage offset applied to each recruitment parameter, on the natural scale.",
      "multiplier", "[nspp, 3, nyrs]", FALSE, NA_character_),
    r("growth_linkage_offset", "linkage",
      "Linkage offset applied to each mean growth parameter, on the log scale.",
      "log scale", "[nspp, nsex, nyrs, 4]", FALSE, NA_character_),
    r("growth_linkage_offset_nat", "linkage",
      "Linkage offset applied to each mean growth parameter, on the natural scale.",
      "multiplier", "[nspp, nsex, nyrs, 4]", FALSE, NA_character_),

    # -- internal ----------------------------------------------------------
    r("pop_scalar", "internal",
      "Multiplier on user-supplied numbers-at-age; estimated for estDynamics = 2 under predation, 1 otherwise.",
      "multiplier", "[nspp]", TRUE, NA_character_),
    r("rec_srr_single_density", "internal",
      "Flag recording whether the stock-recruit prior was evaluated as a single density.",
      "unitless", "[1]", FALSE, NA_character_)
  )
})


#' Look up what a CEATTLE reported quantity is
#'
#' `fit$quantities` uses the model's own abbreviated names. This returns the
#' table mapping each one to what it means, the units it is in, how it is
#' shaped, whether it holds a standard error, and what the same quantity is
#' called in the NOAA standardized assessment output.
#'
#' @param quantity Report names as they appear in `names(fit$quantities)`,
#'   default `NULL` for every documented quantity.
#' @param process Character vector restricting the table to one or more parts of
#'   the model: `"population"`, `"recruitment"`, `"mortality"`, `"fishing"`,
#'   `"reference_points"`, `"catchability"`, `"selectivity"`, `"growth"`,
#'   `"composition"`, `"predation"`, `"likelihood"`, `"linkage"` or
#'   `"internal"`.
#'
#' @return A data frame with one row per quantity and columns `quantity`,
#'   `process`, `meaning`, `units`, `dims`, `se` and `standard_label`.
#'
#' @details
#' Units follow the model throughout: numbers-at-age are **thousands of fish**
#' and weight-at-age is **kg**, so their product is **mt**. A quantity whose
#' units read `"mt or thousands of fish"` follows its fleet's
#' `Observation_units` column.
#'
#' `se = TRUE` means the TMB template `ADREPORT`s the quantity, so `fit$sdrep`
#' holds a standard error for it and [as.data.frame.Rceattle()] can fill `se`,
#' `lwr` and `upr`. `se = FALSE` means no standard error exists anywhere on the
#' fit for that quantity. Nothing has a standard error when the fit was produced
#' with `fit_control(getsd = FALSE)`, which leaves `sdrep` NULL.
#'
#' `standard_label` gives the name the quantity takes in the NOAA standardized
#' assessment output that `stockplotr` and `asar` consume, so a table can be
#' relabelled for those tools; `NA` marks a quantity that standard has no
#' equivalent for. Rceattle's own names are kept as the primary vocabulary
#' because several quantities (predation, multispecies) have no standard name.
#'
#' Every per-recruit reference point (`SPR0`, `SPRlimit`, `SPRtarget`,
#' `SPRFinit`, `NbyageSPR`) is computed only under `msmMode = 0` and is exactly
#' **zero on a multispecies fit**: M there holds predation mortality, which
#' scales with predator abundance, so spawning output per recruit is not a
#' property of the prey stock alone.
#'
#' @examples
#' # The whole table
#' head(quantity_dictionary())
#'
#' # What is ssb_depletion, and what units is it in?
#' quantity_dictionary("ssb_depletion")
#'
#' # Everything that holds a standard error
#' dict <- quantity_dictionary()
#' dict[dict$se, c("quantity", "meaning")]
#'
#' # Every reference point
#' quantity_dictionary(process = "reference_points")
#'
#' # Search the meanings by keyword
#' dict[grep("depletion", dict$meaning, ignore.case = TRUE), ]
#' @seealso [parameter_dictionary()] for the estimated parameters,
#'   [as.data.frame.Rceattle()] for these quantities in tidy form, and
#'   [fit_mod()] for the fitted object they come from.
#' @export
quantity_dictionary <- function(quantity = NULL, process = NULL) {
  out <- .QUANT_INFO

  if (!is.null(process)) {
    known <- unique(.QUANT_INFO$process)
    unknown <- setdiff(process, known)
    if (length(unknown)) {
      stop("Unknown process: ", paste(unknown, collapse = ", "),
           ". Use one of: ", paste(known, collapse = ", "))
    }
    out <- out[out$process %in% process, , drop = FALSE]
  }

  if (!is.null(quantity)) {
    # Warn rather than stop, so passing every name from a fit still works when
    # that fit carries quantities from a branch this dictionary predates.
    unknown <- setdiff(quantity, .QUANT_INFO$quantity)
    if (length(unknown)) {
      warning("Not in the quantity dictionary: ",
              paste(unknown, collapse = ", "), call. = FALSE)
    }
    out <- out[out$quantity %in% quantity, , drop = FALSE]
  }

  rownames(out) <- NULL
  out
}
