## Time-varying selectivity options in Rceattle
##
## A runnable walkthrough of every Time_varying_sel option on GOApollock, kept
## as a developer example alongside two-sex-selectivity-example.R. Not run by
## the suite -- it fits real models.
##
## What it exercises:
##   1. NonParametric + IID, as a FIXED-effect fit. Under random_sel = TRUE this
##      pairing is refused (.rce_np_unintegrable_fleets(), R/3-build_map.R): the
##      always-on average-selectivity penalty tilts the integrated density. Use
##      Selectivity = "NonParametricIntegrable" (form 13, shipped in 5.40.0) for
##      deviations with a proper density that random_sel can integrate.
##   2. NonParametric + RandomWalk, whose "NA/NaN function evaluation" overflow
##      is handled by log_sum_exp.
##
## Every fit below passes random_sel = FALSE for that reason. Some models may
## not converge, so read the warnings.

library(Rceattle)
data("GOApollock")
d <- GOApollock

## 0. Model setup ----
d$env_data$QcovPol <- scale(d$env_data$BTempC)
SHELIKOF <- 1L; BOTTOM_TRAWL <- 2L; FISHERY <- 8L; ADFG <- 3L # Fleet indices

# * Control updates ----
# Updating the control from old setup
d$fleet_control$Catchability[SHELIKOF]   <- "Estimated"
d$fleet_control$Time_varying_q[SHELIKOF] <- "Off"
d$fleet_control$Selectivity[6] <- "DescendingLogistic"
d$fleet_control$Catchability[BOTTOM_TRAWL] <- "Estimated-with-prior"
d$fleet_control$Catchability_init[BOTTOM_TRAWL] <- 0.85
d$fleet_control$Catchability_prior_sd[BOTTOM_TRAWL] <- 0.1

# Composition young-age accumulation (fishery age-1 -> 2; Shelikof ages 1-2 -> 3).
d$fleet_control$Comp_accum_young <- 1L
d$fleet_control$Comp_accum_old   <- 0L
d$fleet_control$Comp_accum_young[FISHERY]  <- 2L
d$fleet_control$Comp_accum_young[SHELIKOF] <- 3L

# Time-varying bits off
d$fleet_control$Time_varying_q[ADFG] <- "Off"   # Using rw(1 | Year) on q, below
d$fleet_control$Time_varying_sel[FISHERY] <- "Off"  # Changing below


# * Priors ----
# Setting up Cole's priors to help with logistic selectivity identifiability
SHELIKOF_ACOUSTIC <- "Pollock_survey_1_shelikof_acoustic"
ADFG_SURVEY       <- "Pollock_survey_3_adfg"

ASC_LIMB_PRIOR_FLEETS <- c("Pollock_survey_2_bottom_trawl",
                           "Pollock_survey_3_adfg")
DESC_LIMB_PRIOR_FLEETS <- c("Pollock_survey_1_shelikof_acoustic",
                            "Pollock_survey_6_summer_acoustic")

# * Catchability: QAR1 (Shelikof) + random-walk q (adfg) ----
# Shelikof is the state-space QAR1; the adfg survey is a random walk on q. The
# random-walk SD is fixed at 0.05.
q_spec <- build_catchability(linkages = list(
  q = list(
    linkage_spec(~ ar1(1 | Year),
                 fleet = SHELIKOF_ACOUSTIC,
                 observe = "QcovPol", # Column in env_data you fit to
                 obs_sd = 0.1),
    linkage_spec(~ rw(1 | Year),
                 fleet = ADFG_SURVEY,
                 init = list(sigma = 0.05)))))

# * Selectivity priors + fishery ascending random walk ----
# Priors mirror goa_pk
sel_spec <- build_selectivity(linkages = list(
  slp_asc  = list(
    linkage_spec(~ 1,
                 fleet = ASC_LIMB_PRIOR_FLEETS,
                 priors = list(intercept = lognormal(-1, 1.5)))),
  inf_asc  = list(
    linkage_spec(~ 1,
                 fleet = ASC_LIMB_PRIOR_FLEETS,
                 priors = list(intercept = normal(0, 3)))),
  slp_desc = linkage_spec(~ 1,
                          fleet = DESC_LIMB_PRIOR_FLEETS,
                          priors = list(intercept = lognormal(-1, 1.5))),
  inf_desc = linkage_spec(~ 1,
                          fleet = DESC_LIMB_PRIOR_FLEETS,
                          priors = list(intercept = normal(10, 3)))))

## --- 1. Non-parametric, Ianelli style, IID deviations (AMAK) --------------
## The non-parametric penalties are Gaussian and can be input as a
## weight (Sel_curve_pen1/2) or as a standard deviation (Sel_shape_sd / Sel_curvature_sd). Sel_shape_dir is the direction of the one-sided shape penalty ("Decreasing"/"Increasing"). Don't set both weights and SDs.

d1 <- d
d1$fleet_control$Selectivity[FISHERY]         <- "NonParametric"
d1$fleet_control$Sel_shape_sd[FISHERY]        <- 1 / sqrt(2 * 1)     # Sel_curve_pen1 = 1
d1$fleet_control$Sel_shape_dir[FISHERY]       <- "Decreasing"
d1$fleet_control$Sel_curvature_sd[FISHERY]    <- 1 / sqrt(2 * 0.05)  # Sel_curve_pen2 = 0.05
d1$fleet_control$N_sel_bins[FISHERY]    <- 8
d1$fleet_control$Time_varying_sel[FISHERY]    <- "IID"
d1$fleet_control$Time_varying_sel_sd[FISHERY] <- 0.35
m1 <- fit_mod(data_list = d1, msmMode = 0, estimateMode = "Hindcast",
              qFun = q_spec,
              selFun = sel_spec,
              fit_control = fit_control(phase = TRUE))

## --- 2. Non-parametric, Ianelli style, random walk (AMAK) ---------
d2 <- d
d2$fleet_control$Selectivity[FISHERY]         <- "NonParametric"
d2$fleet_control$Sel_shape_sd[FISHERY]        <- 1 / sqrt(2 * 2)     # Sel_curve_pen1 = 2
d2$fleet_control$Sel_shape_dir[FISHERY]       <- "Decreasing"
d2$fleet_control$Sel_curvature_sd[FISHERY]    <- 1 / sqrt(2 * 0.05)  # Sel_curve_pen2 = 0.05
d2$fleet_control$N_sel_bins[FISHERY]    <- 8
d2$fleet_control$Time_varying_sel[FISHERY]    <- "RandomWalk"
d2$fleet_control$Time_varying_sel_sd[FISHERY] <- 0.1
m2 <- fit_mod(data_list = d2, msmMode = 0, estimateMode = "Hindcast",
              qFun = q_spec,
              selFun = sel_spec,
              fit_control = fit_control(phase = TRUE))

## 3. Non-parametric, ADMB "pm" form, random walk ----
## Three additional switches
##
##   Sel_avgsel_pen  AMAK's base-level penalty, weight * log(mean(exp(base)))^2.
##   Sel_cap_bin     Bin above which the realized curve is held flat.
##   Sel_start_year  Random walk start year. NA = earliest year with data.
d3 <- d
d3$fleet_control$Selectivity[FISHERY]         <- "NonParametricPM"
d3$fleet_control$Sel_shape_sd[FISHERY]        <- 1 / sqrt(2 * 2)     # Sel_curve_pen1 = 2
d3$fleet_control$Sel_shape_dir[FISHERY]       <- "Decreasing"
d3$fleet_control$Sel_curvature_sd[FISHERY]    <- 1 / sqrt(2 * 0.05)  # Sel_curve_pen2 = 0.05
d3$fleet_control$N_sel_bins[FISHERY]    <- 8
d3$fleet_control$Time_varying_sel[FISHERY]    <- "RandomWalk"
d3$fleet_control$Time_varying_sel_sd[FISHERY] <- 0.1
d3$fleet_control$Sel_avgsel_pen               <- rep(0, nrow(d3$fleet_control))
d3$fleet_control$Sel_avgsel_pen[FISHERY]      <- 10   # 0 = off; 10 matches AMAK
d3$fleet_control$Sel_cap_bin                  <- rep(NA, nrow(d3$fleet_control))
d3$fleet_control$Sel_start_year               <- rep(NA, nrow(d3$fleet_control))
d3$fleet_control$Sel_start_year[FISHERY]      <- 1990
m3 <- fit_mod(data_list = d3, msmMode = 0, estimateMode = "Hindcast",
              qFun = q_spec,
              selFun = sel_spec,
              fit_control = fit_control(phase = TRUE))

## 4. Hake selectivity (Taylor), IID deviations (SS3) ----
d4 <- d
d4$fleet_control$Selectivity[FISHERY]         <- "Hake"
d4$fleet_control$N_sel_bins[FISHERY]    <- 8
d4$fleet_control$Time_varying_sel[FISHERY]    <- "IID"
d4$fleet_control$Time_varying_sel_sd[FISHERY] <- 0.35
m4 <- fit_mod(data_list = d4, msmMode = 0, estimateMode = "Hindcast",
              qFun = q_spec,
              selFun = sel_spec,
              fit_control = fit_control(phase = TRUE))

## --- 5. Double normal, IID deviations ----
d5 <- d
d5$fleet_control$Selectivity[FISHERY]         <- "DoubleNormal"
d5$fleet_control$Time_varying_sel[FISHERY]    <- "IID"
d5$fleet_control$Time_varying_sel_sd[FISHERY] <- 0.05
m5 <- fit_mod(data_list = d5, msmMode = 0, estimateMode = "Hindcast",
              qFun = q_spec,
              selFun = sel_spec,
              fit_control = fit_control(phase = TRUE))

## --- 6. 2DAR1 over age x year ----
## Sel_curve_pen1/2 become estimated correlations on the logit scale.
## Time_varying_sel is ignored: the field always has deviations, at sd Time_varying_sel_sd.
d6 <- d
d6$fleet_control$Selectivity[FISHERY]      <- "2DAR1"
d6$fleet_control$N_sel_bins[FISHERY]       <- 8
d6$fleet_control$Time_varying_sel[FISHERY] <- "Off"
d6$fleet_control$Sel_norm_bin[FISHERY]     <- 0   # normalize by max (moving to "max" or 0 on next release)
d6$fleet_control$Sel_curve_pen1[FISHERY]   <- 0   # correlation across bins
d6$fleet_control$Sel_curve_pen2[FISHERY]   <- 0   # correlation across years
d6$fleet_control$Time_varying_sel_sd[FISHERY]   <- 0.2  # Fixed conditional sd
m6 <- fit_mod(data_list = d6, msmMode = 0, estimateMode = "Hindcast",
              qFun = q_spec,
              selFun = sel_spec,
              fit_control = fit_control(phase = TRUE),
              random_sel = FALSE)                # Set to true to estimate sigma

## --- 7. 3DAR1 over age x year x cohort (Cheng et al. 2024) ----
## As 2DAR1, plus Sel_curve_pen3 as the cohort correlation.
d7 <- d
d7$fleet_control$Selectivity[FISHERY]      <- "3DAR1"
d7$fleet_control$N_sel_bins[FISHERY]       <- 8
d7$fleet_control$Time_varying_sel[FISHERY] <- "Off"
d7$fleet_control$Sel_norm_bin[FISHERY]     <- 0    # normalize by max (moving to "max" or 0 on next release)
d7$fleet_control$Sel_curve_pen1[FISHERY]   <- 0   # correlation across bins
d7$fleet_control$Sel_curve_pen2[FISHERY]   <- 0   # correlation across years
d7$fleet_control$Sel_curve_pen3[FISHERY]   <- 0   # correlation across cohorts
d7$fleet_control$Time_varying_sel_sd[FISHERY]   <- 0.2  # Fixed conditional sd
m7 <- fit_mod(data_list = d7, msmMode = 0, estimateMode = "Hindcast",
              qFun = q_spec,
              selFun = sel_spec,
              fit_control = fit_control(phase = TRUE),
              random_sel = FALSE)

mod_list <- list(NP_IID      = m1,
                 NP_RW       = m2,
                 NPpm_RW     = m3,
                 Hake_IID    = m4,
                 DblNorm_IID = m5,
                 AR2D        = m6,
                 AR3D        = m7)

## Comparison ----
tabs <- report_tables(mod_list)
tabs$model[, c("model", "marginal_nll", "joint_nll",
               "n_random", "max_gradient", "pdHess")]

plot_biomass(mod_list)
plot_selectivity(mod_list)

## Fitted catchability by year, for every fleet with index_data.
plot_catchability(m2)

## Report ----
## Meaning, units, dimensions, and whether sdreport() gives a standard error.
quantity_dictionary(process = "selectivity")
quantity_dictionary()[quantity_dictionary()$se, "quantity"]

## Convergence check.
m1$convergence
