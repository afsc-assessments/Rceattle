## Time-varying selectivity options in Rceattle on two-sex model
## Also options for different scale between sexes
##
## Status note (re-checked 2026-09-21, 5.41.0). The header used to say
## "NonParametric + IID now works". It does not, and deliberately so:
## 5.35.0 REFUSES random_sel = TRUE on NonParametric and NonParametricPM,
## because the shape penalties are charged on each year's realized curve,
## so the reported deviation sd is not the sd of the deviations. 5.40.0
## added the form that does integrate, NonParametricIntegrable (13); use it
## when you want random_sel = TRUE.
## The sections below run because random_sel defaults to FALSE.
##
## Some models may not converge, so check warning!

library(Rceattle)
data("GOAatf")
d <- GOAatf

## 0. Model setup ----
FISHERY <- 3L

# Time-varying
d$fleet_control$Time_varying_q <- "Off"
d$fleet_control$Time_varying_sel <- "Off"
d$fleet_control$Sel_curve_pen2 <- NA
d$fleet_control$Sel_curve_pen1 <- NA

## --- 1. Non-parametric, Ianelli style, IID deviations (AMAK) --------------
## The non-parametric penalties are normal and can be input as a
## WEIGHT (Sel_curve_pen1/2) or as a standard deviation (Sel_shape_sd / Sel_curvature_sd). Sel_shape_dir is the direction of the one-sided shape penalty ("Decreasing"/"Increasing"). Don't set both weights and SDs.

d1 <- d
d1$fleet_control$Selectivity[FISHERY]         <- "NonParametric"
d1$fleet_control$Sel_shape_sd[FISHERY]        <- 1 / sqrt(2 * 100)     # Sel_curve_pen1 = 100
d1$fleet_control$Sel_shape_dir[FISHERY]       <- "Decreasing"
d1$fleet_control$Sel_curvature_sd[FISHERY]    <- 1 / sqrt(2 * 40)  # Sel_curve_pen2 = 40
d1$fleet_control$N_sel_bins[FISHERY]    <- 19
d1$fleet_control$Time_varying_sel[FISHERY]    <- "IID"
d1$fleet_control$Time_varying_sel_sd[FISHERY] <- 0.05
m1 <- fit_mod(data_list = d1, msmMode = 0, estimateMode = "Hindcast",
              fit_control = fit_control(phase = TRUE))

## --- 2. Non-parametric, Ianelli style, random walk (AMAK) ---------
d2 <- d
d2$fleet_control$Selectivity[FISHERY]         <- "NonParametric"
d2$fleet_control$Sel_shape_sd[FISHERY]        <- 1 / sqrt(2 * 100)     # Sel_curve_pen1 = 100
d2$fleet_control$Sel_shape_dir[FISHERY]       <- "Decreasing"
d2$fleet_control$Sel_curvature_sd[FISHERY]    <- 1 / sqrt(2 * 40)  # Sel_curve_pen2 = 40
d2$fleet_control$N_sel_bins[FISHERY]    <- 19
d2$fleet_control$Time_varying_sel[FISHERY]    <- "RandomWalk"
d2$fleet_control$Time_varying_sel_sd[FISHERY] <- 0.05
m2 <- fit_mod(data_list = d2, msmMode = 0, estimateMode = "Hindcast",
              fit_control = fit_control(phase = TRUE))

## 3. Non-parametric, ADMB "pm" form, random walk ----
## Three additional switches
##
##   Sel_avgsel_pen  AMAK's base-level penalty, weight * log(mean(exp(base)))^2.
##   Sel_cap_bin     Bin above which the realized curve is held flat.
##   Sel_start_year  Random walk start year. NA = earliest year with data.
d3 <- d
d3$fleet_control$Selectivity[FISHERY]         <- "NonParametricPM"
d3$fleet_control$Sel_shape_sd[FISHERY]        <- 1 / sqrt(2 * 100)     # Sel_curve_pen1 = 2
d3$fleet_control$Sel_shape_dir[FISHERY]       <- "Decreasing"
d3$fleet_control$Sel_curvature_sd[FISHERY]    <- 1 / sqrt(2 * 40)  # Sel_curve_pen2 = 0.05
d3$fleet_control$N_sel_bins[FISHERY]    <- 19
d3$fleet_control$Time_varying_sel[FISHERY]    <- "RandomWalk"
d3$fleet_control$Time_varying_sel_sd[FISHERY] <- 0.1
d3$fleet_control$Sel_avgsel_pen               <- rep(0, nrow(d3$fleet_control))
d3$fleet_control$Sel_avgsel_pen[FISHERY]      <- 10   # 0 = off; 10 matches AMAK
d3$fleet_control$Sel_cap_bin                  <- rep(NA, nrow(d3$fleet_control))
d3$fleet_control$Sel_start_year               <- rep(NA, nrow(d3$fleet_control))
d3$fleet_control$Sel_start_year[FISHERY]      <- 1990
m3 <- fit_mod(data_list = d3, msmMode = 0, estimateMode = "Hindcast",
              fit_control = fit_control(phase = TRUE))

## 4. Hake selectivity (Taylor), IID deviations (SS3) ----
d4 <- d
d4$fleet_control$Selectivity[FISHERY]         <- "Hake"
d4$fleet_control$N_sel_bins[FISHERY]    <- 19
d4$fleet_control$Time_varying_sel[FISHERY]    <- "IID"
d4$fleet_control$Time_varying_sel_sd[FISHERY] <- 0.35
m4 <- fit_mod(data_list = d4, msmMode = 0, estimateMode = "Hindcast",
              fit_control = fit_control(phase = TRUE))

## --- 5. Double normal, IID deviations ----
d5 <- d
d5$fleet_control$Selectivity[FISHERY]         <- "DoubleNormal"
d5$fleet_control$Time_varying_sel[FISHERY]    <- "IID"
d5$fleet_control$Time_varying_sel_sd[FISHERY] <- 0.05
m5 <- fit_mod(data_list = d5, msmMode = 0, estimateMode = "Hindcast",
              fit_control = fit_control(phase = TRUE))

## --- 6. 2DAR1 over age x year ----
## Sel_curve_pen1/2 become estimated correlations on the logit scale.
## Time_varying_sel is ignored: the field supplies its own deviations and sd.
d6 <- d
d6$fleet_control$Selectivity[FISHERY]      <- "2DAR1"
d6$fleet_control$N_sel_bins[FISHERY]       <- 19
d6$fleet_control$Time_varying_sel[FISHERY] <- "Off"
d6$fleet_control$Sel_norm_bin[FISHERY]     <- 0    # normalize by max ("Max" is accepted by Sel_norm_bin since 5.35.0)
d6$fleet_control$Sel_norm_scope[FISHERY]   <- "AcrossSexes" # compared to "WithinSex"
d6$fleet_control$Sel_curve_pen1[FISHERY]   <- 0   # correlation across bins

d6$fleet_control$Sel_curve_pen2[FISHERY]   <- 0   # correlation across years
d6$fleet_control$Time_varying_sel_sd[FISHERY]   <- 0.1  # Fixed conditional sd
m6 <- fit_mod(data_list = d6, msmMode = 0, estimateMode = "Hindcast",
              fit_control = fit_control(phase = TRUE),
              random_sel = FALSE)                # Set to true to estimate sigma

## --- 7. 3DAR1 over age x year x cohort (Cheng et al. 2024) ----
## As 2DAR1, plus Sel_curve_pen3 as the cohort correlation.
d7 <- d
d7$fleet_control$Selectivity[FISHERY]      <- "3DAR1"
d7$fleet_control$N_sel_bins[FISHERY]       <- 19
d7$fleet_control$Time_varying_sel[FISHERY] <- "Off"
d7$fleet_control$Sel_norm_bin[FISHERY]     <- 0    # normalize by max ("Max" is accepted by Sel_norm_bin since 5.35.0)
d7$fleet_control$Sel_norm_scope[FISHERY]   <- "AcrossSexes" # compared to "WithinSex"
d7$fleet_control$Sel_curve_pen1[FISHERY]   <- 0   # correlation across bins
d7$fleet_control$Sel_curve_pen2[FISHERY]   <- 0   # correlation across years
d7$fleet_control$Sel_curve_pen3[FISHERY]   <- 0   # correlation across cohorts
d7$fleet_control$Time_varying_sel_sd[FISHERY]   <- 0.1  # Fixed conditional sd
m7 <- fit_mod(data_list = d7, msmMode = 0, estimateMode = "Hindcast",
              fit_control = fit_control(phase = TRUE),
              random_sel = FALSE)

mod_list <- list(NP_IID      = m1,
                 NP_RW       = m2,
                 NPpm_RW     = m3,
                 Hake_IID    = m4,
                 DblNorm_IID = m5,
                 AR2D        = m6,
                 AR3D        = m7)

## 8. Different max selectivity for males and females --------------
## May won't work for the following:
##   Logistic / DescendingLogistic / LogisticPM   both sexes asymptote to 1, but could asymptote lower
##   NonParametric / NonParametricPM              each sex re-centered to mean 1
##   Hake                                         each sex scaled by its own max
##                                                (Sel_norm_scope is inert here --
##                                                 inst/dev/TODO-selectivity.md) (open: the capability gap, not the silence --
##                                                 data_check() announces it since 5.35.0)
##   DoubleNormal                                 both sexes peak at exactly 1
##                                                (only the old-age plateau differs)
##
## DoubleLogistic, 2DAR1 and 3DAR1 work when Sel_norm_scope = "AcrossSexes"

## Estimated maximum selectivity per sex.
sex_max <- function(fit, flt = FISHERY, yr = 1) {
  s <- fit$quantities$sel_at_age[flt, , , yr]
  c(female = max(s[1, ]), male = max(s[2, ]),
    ratio  = max(s[2, ]) / max(s[1, ]))
}


## 8a. Logistic ----
# - Both sexes max at 1
d8 <- d
d8$fleet_control$Selectivity[FISHERY]         <- "Logistic"
d8$fleet_control$Sel_norm_bin[FISHERY]     <- 0    # normalize by max ("Max" is accepted by Sel_norm_bin since 5.35.0)
d8$fleet_control$Sel_norm_scope[FISHERY]   <- "AcrossSexes" # compared to "WithinSex"
m8 <- fit_mod(data_list = d8, msmMode = 0, estimateMode = "Hindcast",
              fit_control = fit_control(phase = TRUE))

## 8b. Double logistic ----
d9 <- d
d9$fleet_control$Selectivity[FISHERY]         <- "DoubleLogistic"
d9$fleet_control$Sel_norm_bin[FISHERY]     <- 0    # normalize by max ("Max" is accepted by Sel_norm_bin since 5.35.0)
d9$fleet_control$Sel_norm_scope[FISHERY]   <- "AcrossSexes" # compared to "WithinSex"
m9 <- fit_mod(data_list = d9, msmMode = 0, estimateMode = "Hindcast",
              fit_control = fit_control(phase = TRUE))

## 8c. Double logistic, constraining the lesser-selected sex ----
# Nothing above needed a linkage: every form is sex-specific already, so 8b's
# two sexes each have their own four parameters and find their own peak height.
# A linkage is for saying which sex is less selected, rather than letting two
# free curves decide -- the male-as-offset-from-female setup in Stock Synthesis.
# Put the prior on that sex alone with by = ~ fleet + sex.
# Needs Rceattle >= 5.28.0: before that a per-sex selectivity linkage also fixed
# the OTHER sex's parameter, which is the reference the offset is measured from.
sel_off <- build_selectivity(linkages = list(
  inf_desc = linkage_spec(~ 1, by = ~ fleet + sex,
                          fleet = FISHERY, sex = "male",   # male/2 = the male curve
                          priors = list(intercept = normal(8, 0.1)))))

d10 <- d
d10$fleet_control$Selectivity[FISHERY]        <- "DoubleLogistic"
d10$fleet_control$Sel_norm_bin[FISHERY]       <- 0
d10$fleet_control$Sel_norm_scope[FISHERY]     <- "AcrossSexes"
d10$fleet_control$Selectivity_index[FISHERY]  <- FISHERY  # a prior needs the lead fleet
m10 <- fit_mod(data_list = d10, msmMode = 0, estimateMode = "Hindcast",
               selFun = sel_off,
               fit_control = fit_control(phase = TRUE))

# A loose prior is worse than none: N(8, 1) here reaches a ratio of 2.80 with no
# sdreport at all. Check the Hessian before believing a sex ratio.
m10$convergence$checks$hessian_conditioning

## 8d. Logistic with a per-sex apical height (5.38.0) ----
# 8a cannot produce a sex level difference at all: both logistic curves
# asymptote to 1. `apical` multiplies one sex's whole curve, after the form and
# before the across-sex normalization, so the sex contrast becomes a parameter
# with a standard error. Name the fleet and the sex that carries it; the other
# sex is the reference. For two saturating logistic curves the multiplier is
# also the ratio of maxima; for a dome with sex-specific shape it is not, and
# sex_max() is the ratio to quote.
sel_apical <- build_selectivity(linkages = list(
  apical = linkage_spec(~ 1, by = ~ fleet + sex, fleet = FISHERY, sex = "male",
                        priors = list(intercept = lognormal(0, 0.5)))))
d11 <- d8
d11$fleet_control$Selectivity_index[FISHERY] <- FISHERY   # the offset sits on the lead fleet
m11 <- fit_mod(data_list = d11, msmMode = 0, estimateMode = "Hindcast",
               selFun = sel_apical,
               fit_control = fit_control(phase = TRUE))
sex_max(m11)                                             # 5.38.0: ratio 2.15
exp(m11$estimated_params$log_sel_apical[FISHERY, 2])     # the multiplier; equal here, logistic
summary(m11$sdrep)["log_sel_apical", ]                   # log 0.77, SE 0.21

## Sex comparison ----
## The objective for 8c is not comparable to 8a/8b -- a prior adds its own term,
## so they are different objective functions rather than better and worse fits.
sex_mods <- list(Logistic        = m8,
                 DblLogis_across = m9,
                 DblLogis_offset = m10,
                 AR2D            = m6)
round(t(sapply(sex_mods, sex_max)), 4)
plot_selectivity(sex_mods)

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
