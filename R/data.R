#' Data inputs for single species CEATTLE of the Bering Sea from 1979 to 2017
#'
#' A data list of inputs for the single-species Eastern Bering Sea CEATTLE model
#' (walleye pollock, Pacific cod, and arrowtooth flounder), covering 1979-2017.
#'
#' @format
#' Control
#' \describe{
#' \item{nspp}{Number of species included in CEATTLE}
#' \item{styr}{Start year of the hindcast}
#' \item{endyr}{End year of the hindcast}
#' \item{projyr}{End year of the forecast}
#' \item{nsex}{Number of sexes to model in the population (1 = sexes combined, 2 = females and males modeled separately)}
#' \item{spawn_month}{Spawning month of the population to adjust the numbers spawning}
#' \item{R_sexr}{Percent of recruitment that is female (ignored if nsex = 1)}
#' \item{nages}{Number of ages of each species included in the hindcast}
#' \item{minage}{Minimum age for each population (i.e. does recruitment correspond to age 0, 1, 2?)}
#' \item{nlengths}{Number of lengths of each species included in the hindcast}
#' \item{pop_wt_index}{Weight-at-age (weight) index to use for calculation of each species population derived quantities (SSB, Consumption/Ration, Suitability, etc)}
#' \item{ssb_wt_index}{Weight-at-age index used to compute each species' female spawning-stock biomass (SSB)}
#' \item{pop_age_transition_index}{Age transition matrix (e.g. growth trajectory) index to use for derived quantities of the population to convert age to length (also used in length-based predation estimation)}
#' \item{sigma_rec}{Standard deviation of the log-scale recruitment deviations (fixed, or an initial value if estimated)}
#' \item{other_food}{Other food in the ecosystem for each species}
#' \item{estDynamics}{Estimate or fix numbers-at-age. Accepts integer codes or the equivalent readable strings: 0 / "Estimated" = estimate dynamics; 1 / "Fixed" = use input numbers-at-age in NByageFixed; 2 / "FixedScaled" = multiply input numbers-at-age (NByageFixed) by a single scaling coefficient, estimated under predation (msmMode > 0) and fixed at 1 otherwise.}
#' \item{M1_model}{Estimate residual (multi-species mode) or total natural mortality (single-species mode). 0 = use fixed natural mortality from M1_base, 1 = estimate sex- and age-invariant M1, 2 = sex-specific (two-sex model), age-invariant M1, 3 = estimate sex- and age-specific M1.}
#' \item{fleet_control}{Survey and fishery data specifications}
#' \item{index_data}{Survey index in weight (mt) or numbers (thousands of fish) data}
#' \item{index_cov}{Optional named list of survey-index variance-covariance matrices, keyed by Fleet_name (or Fleet_code), used only by fleets whose fleet_control \code{Index_distribution == "MVN"}. Each matrix must be square and symmetric with dimension equal to the number of fitted survey observations for that fleet (Year in \code{styr:endyr}, Observation > 0), ordered as in index_data. Inverted once internally to the precision matrix used in 0.5 * r' Sigma^-1 r. Leave unset for the default lognormal likelihood.}
#' \item{catch_data}{Total catch in weight (mt) or numbers (thousands of fish) data}
#' \item{comp_data}{Survey/fishery age or length composition data. Note if sex is 3, put female composition data then male composition data (similar to SS).}
#' \item{emp_sel}{Empirical/fixed selectivity for surveys and fisheries (leave empty if not used)}
#' \item{age_trans_matrix}{Age transition matrix (e.g. growth trajectory) used to convert age to length for length comp data. Can have multiple matrices for a species specified by Age_transition_index.}
#' \item{age_error}{Aging error matrices. Can have only one per species.}
#' \item{weight}{Weight-at-age (weight) to use for calculation of derived quantities (SSB, Consumption/Ration, Suitability, Total Catch, Survey Biomass, etc). Can have multiple weight-at-age data-sets for each species, but must be full for all years of the hindcast.}
#' \item{maturity}{Maturity-at-age for each species}
#' \item{sex_ratio}{Percent female at age for each species}
#' \item{M1_base}{Residual natural mortality for each species}
#' \item{aLW}{Parameters for weight-at-length power function for each species. Used when estimating time-variant length-based gamma suitability (suitMode = 2) or time-variant length-based lognormal suitability (suitMode = 5)}
#' \item{Ceq}{Which bioenergetics equation to use for each species for ft to scale max consumption: 1 = Exponential (Stewart et al 1983), 2 = Temperature-dependence for warm-water species (Kitchell et al 1977; sensu Holsman et al 2015), 3 = temperature dependence for cool and cold-water species (Thornton and Lessem 1979)}
#' \item{Cindex}{Which environmental index in env_data to use to drive bioenergetics}
#' \item{Pvalue}{This scales the maximum consumption used for ration for each species; Pvalue is in Cmax*fT*Pvalue*ration_data}
#' \item{fday}{Number of foraging days per year for each species}
#' \item{CA}{Intercept of allometric mass function for calculating maximum consumption: CA * Weight ^ CB}
#' \item{CB}{Slope of allometric mass function for calculating maximum consumption: CA * Weight ^ CB}
#' \item{Qc}{Parameter for temperature scaling function of maximum consumption specified by Ceq}
#' \item{Tco}{Parameter for temperature scaling function of maximum consumption specified by Ceq}
#' \item{Tcm}{Parameter for temperature scaling function of maximum consumption specified by Ceq}
#' \item{Tcl}{Parameter for temperature scaling function of maximum consumption specified by Ceq}
#' \item{CK1}{Parameter for temperature scaling function of maximum consumption specified by Ceq}
#' \item{CK4}{Parameter for temperature scaling function of maximum consumption specified by Ceq}
#' \item{env_data}{Environmental indices such as bottom temperature data to incorporate into ration equation specified by Ceq and Cindex. Also used to drive catchability if Catchability = 5 or 6. Will use the mean for missing years. Temperature should be in Celsius.}
#' \item{ration_data}{Annual relative foraging rate by age or input consumption at age. Multiplied by pvalue and fday to scale maximum consumption to the number of days in a year that foraging occurs.}
#' \item{UobsAge}{Stomach proportion by numbers for each predator, prey, predator age, prey age combination. Can also be year specific by including the column, "Year"}
#' \item{UobsWtAge}{Stomach proportion by weight for each predator, prey, predator age, prey age combination. Can also be year specific by including the column, "Year"}
#' }
#'
#' fleet_control: controls for survey and fisheries data
#'\describe{
#' \item{Fleet_name}{Name of survey or fishery}
#' \item{Fleet_code}{Index of survey/fishery ACROSS species}
#' \item{Fleet_type}{0 or 'Off' = Do not estimate; 1 = 'Fishery'; 2 = 'Survey'}
#' \item{Species}{Species number}
#' \item{Selectivity_index}{Index used to give fleets the same selectivity (otherwise, same as Fleet_code). Fleets sharing a value share one selectivity parameter block, with its penalties and priors accumulated once on the group's first non-Off fleet. Sharing the parameters is not enough on its own: the columns that shape the curve are read per fleet and must agree across the group, or the fleets end up with different selectivities. \code{data_check()} reports any that differ. To mirror a fleet, copy its fleet_control row and change only the identity and catchability columns. See \code{vignette("model-options-and-functionality")}.}
#' \item{Selectivity}{Selectivity form: 0 = "Fixed"; 1 = "Logistic"; 2 = "NonParametric" (Ianelli et al. 2018); 3 = "DoubleLogistic"; 4 = "DescendingLogistic"; 5 = "Hake" non-parametric (Taylor et al. 2014); 6 = "2DAR1"; 7 = "3DAR1" (Cheng et al. 2024); 8 = "DoubleNormal"; 9 = "NonParametricPM" (non-parametric random walk with penalized increments, AMAK "pm"; not integrable under random_sel); 11 = "LogisticPM" (logistic with a free age-1 selectivity, AMAK "pm"); 13 = "NonParametricIntegrable" (the Ianelli base curve whose deviations have a proper density, so random_sel integrates them; Time_varying_sel picks the structure: "Off", "IID" or "RandomWalk"). Whether a form is age- or length-based is set by Selectivity_dimension.}
#' \item{Selectivity_dimension}{"Age" or "Length".}
#' \item{N_sel_bins}{Number of bins (ages or length bins) to estimate for non-parametric selectivity.}
#' \item{Time_varying_sel}{Whether a time-varying selectivity should be estimated. 0 = "Off"; 1 = "IID" penalized deviates given \code{Time_varying_sel_sd}, or random effects with that sd estimated under \code{random_sel}; 2 = "AR1" REMOVED in 5.16.0, it was scored as "IID", so \code{data_check()} errors on it; use "IID" for the same fit, or \code{linkage_spec(~ ar1(1 | Year))} for a real one; 3 = "Block" time blocks with no penalty; 4 = "RandomWalk" following Dorn; 5 = "RandomWalkAscending" on the ascending limb of double logistic only.}
#' \item{Time_varying_sel_sd}{The fixed or initial sd to use for time varying selectivity.}
#' \item{Sel_start_year}{First year the fleet's time-varying selectivity deviations are estimated and penalized. Deviations before this year are fixed at 0: they have no data and no penalty (the selectivity penalties are anchored at this year), so estimating them would leave unidentified flat directions. Defaults to the earliest year of data across all fleets sharing the fleet's \code{Selectivity_index} (fleets sharing an index share one selectivity curve), bounded below by \code{styr}. Only used for time-varying selectivity; set explicitly to override.}
#' \item{Bin_first_selected}{Age/length bin at which selectivity is non-zero}
#' \item{Sel_shape_sd, Sel_curvature_sd, Sel_devmag_sd}{Selectivity penalties expressed as standard deviations (an intuitive alternative to the raw penalty weights in Sel_curve_pen1/2/3). Each penalty is a Gaussian SSQ, so the weight is \code{1 / (2 * sd^2)}: Sel_shape_sd is the shape/smoothness penalty (Sel_curve_pen1), Sel_curvature_sd the 2nd-difference curvature penalty (Sel_curve_pen2), Sel_devmag_sd the deviation-magnitude penalty (Sel_curve_pen3). Supplying an SD sets the corresponding Sel_curve_pen column; a legacy Sel_curve_pen value is never overwritten. Sel_shape_sd applies to NonParametric (2/9/13) and LogisticPM (11), Sel_devmag_sd to "NonParametricPM" (9) and "LogisticPM" (11) only, the slot being unread by "NonParametric" (2) and "NonParametricIntegrable" (13); Sel_curvature_sd is NonParametric-only (2/9/13; LogisticPM does not use Sel_curve_pen2). The 2DAR1/3DAR1 forms reuse Sel_curve_pen for logit-scale correlations and keep those columns.}
#' \item{Sel_shape_dir}{Direction of the shape penalty when Sel_shape_sd is used: "Decreasing" (default, penalize decreasing selectivity) or "Increasing" (sets the sign of Sel_curve_pen1). Only "NonParametricPM" (9) with Sel_shape_mode = "Directional" branches on the sign; every other form, that one under "Smooth" included, applies Sel_curve_pen1 as written, so "Increasing" is refused there rather than negating it.}
#' \item{Sel_curve_pen1}{Raw shape/smoothness penalty weight for non-parametric (type 2/9/13) and LogisticPM (11) selectivity. The intuitive alternative is Sel_shape_sd, which sets this to \code{1 / (2 * sd^2)}. The 2DAR1/3DAR1 forms reuse this column as a logit-scale correlation.}
#' \item{Sel_curve_pen2}{Raw 2nd-difference curvature penalty weight for non-parametric selectivity (the intuitive alternative is Sel_curvature_sd). Reused as a logit-scale correlation by the 2DAR1/3DAR1 forms.}
#' \item{Sel_curve_pen3}{Raw deviation-magnitude penalty weight, read only by "NonParametricPM" (9) and "LogisticPM" (11) (the intuitive alternative is Sel_devmag_sd).}
#' \item{Sel_pen_first_bin}{First bin (age or length, on the fleet's Selectivity_dimension) for the non-parametric shape penalty. NA defaults to Bin_first_selected.}
#' \item{Sel_pen_last_bin}{Last (left) bin of the shape-penalty pairs. NA defaults to nbins - 2.}
#' \item{Sel_shape_mode}{Shape-penalty mode: "Directional" (default) or "Smooth" (two-sided second-difference penalty, RTMB).}
#' \item{Sel_avgsel_pen}{Weight on the AMAK average-selectivity base-level penalty (type 9 only): 0 = off (default), 10 matches AMAK.}
#' \item{Sel_cap_bin}{NonParametricPM selectivity bin cap. NA (default) applies no cap.}
#' \item{Sel_norm_bin}{Where selectivity is normalized to 1. \code{"Max"}
#'   normalizes by the largest value, \code{"Off"} (or blank, the default) not at
#'   all, and \code{"All"} on a LogisticPM fleet covers its whole penalty range.
#'   A number is an absolute age for an age-based fleet (6 means age 6, not the
#'   6th bin) or a 1-based length-bin ordinal for a length-based one, and must
#'   lie in the range the fleet is selected over; below it the curve is zeroed,
#'   so a value there is read as \code{"Max"}. In a two-sex model this says only
#'   where the reference is taken, whether it is pooled across the sexes is
#'   \code{Sel_norm_scope}. See the sex-structure section of
#'   \code{vignette("model-options-and-functionality")}.}
#' \item{Sel_norm_bin_upper}{Optional upper age/length bin. \code{"Off"} (or blank, the default) normalizes at a single bin; a number normalizes by mean selectivity between Sel_norm_bin and this bin, on the same scale and within the same range.}
#' \item{Sel_norm_scope}{Whether selectivity normalization pools its reference
#'   across sexes, orthogonal to \code{Sel_norm_bin} (which says where the
#'   reference is taken). \code{"WithinSex"} divides each sex by its own
#'   reference, so both reach 1 and only the shape differs by sex;
#'   \code{"AcrossSexes"} (default) pools one reference over both sexes, so the
#'   less-selected sex stays below 1 and relative sex-specific selectivity is
#'   retained. No effect on a one-sex species or where \code{Sel_norm_bin} is
#'   NA.}
#' \item{Observation_units}{Units of the observation: 1 = weight (mt), 2 = numbers caught (thousands of fish). Drives both catch and index prediction.}
#' \item{Weight_index}{Weight-at-age (weight) index to use for calculation of derived quantities}
#' \item{Age_transition_index}{Age transition matrix (e.g. growth trajectory) index to use for derived quantities to convert age to length}
#' \item{Catchability_index}{Index used to give fleets the same catchability (otherwise, same as Fleet_code). Fleets sharing a value share one q parameter, so the group holds only one answer to whether q is estimated: the group's first non-Off fleet decides for all of them, regardless of fleet type. Give a fleet its own value when it should be estimated independently. "Analytical" and "AnalyticalArith" are the exception, they solve q from each fleet's own index observations, so a group containing one does \emph{not} share a catchability. \code{data_check()} reports that case, a Fixed lead that leaves fleets on different inits, and a \code{Catchability} or \code{Time_varying_q} that differs within a group. See \code{vignette("model-options-and-functionality")}.}
#' \item{Catchability}{Estimate catchability? 0 or "Fixed" = fixed at prior; 1 or "Estimated" = estimate a single parameter; 2 or "Estimated-with-prior" = estimate a single parameter with a prior; 3 or "Analytical" = geometric-mean analytical q (Ludwig and Walters 1994); 4 = power equation (not yet implemented); 5 or "Environmental" = linear equation log(q_y) = q_mu + beta * index_y; 6 or "AR1" = removed, and \code{data_check()} errors on it, the deviates were never estimated, so q came back constant; express the Rogers et al. (2024) form as a q linkage with ar1(1 | Year) and \code{observe}, passed to \code{fit_mod(qFun = )}; 7 or "AnalyticalArith" = arithmetic-mean analytical q = mean(obs)/mean(pred), the AMAK/ebswp form used with the MVN survey likelihood.}
#' \item{Catchability_init}{Starting value (or fixed value) for catchability. When Catchability = "Estimated-with-prior" this also centers the q prior.}
#' \item{Catchability_prior_sd}{Log-scale standard deviation of the q prior: dnorm(log_q, log(Catchability_init) - bias_adjust_proc * sd^2/2, sd). Used only when Catchability = "Estimated-with-prior".}
#' \item{Time_varying_q}{Whether a time-varying q should be estimated. 0 = "Off"; 1 = "IID" penalized deviate or random effect; 2 = "AR1" REMOVED in 5.16.0, it was scored as "IID", so \code{data_check()} errors on it; use "IID" for the same fit, or \code{linkage_spec(~ ar1(1 | Year))} for a real one; 3 = "Block" time blocks with no penalty; 4 = "RandomWalk" from the mean following Dorn 2018. Under \code{Catchability = "Environmental"} this column instead holds the \code{env_data} column index, not a mode.}
#' \item{Time_varying_q_sd}{The sd to use for the random walk of time varying q if set to 1}
#' \item{Estimate_index_sd}{Estimate survey/index observation SD. Accepts integer codes or readable strings: 0 / "Fixed" = use the fixed SD from the data CV; 1 / "Estimated" = estimate; 2 / "Analytical" = analytically estimate following Ludwig and Walters 1994.}
#' \item{Index_sd}{Starting/fixed value for the index observation SD, used when Estimate_index_sd = "Estimated".}
#' \item{Estimate_catch_sd}{Estimate fishery/catch observation SD. Accepts integer codes or readable strings: 0 / "Fixed" = use the fixed SD from the data CV; 1 / "Estimated" = estimate; 2 / "Analytical" = analytically estimate following Ludwig and Walters 1994.}
#' \item{Catch_sd}{Starting/fixed value for the catch observation SD, used when Estimate_catch_sd = "Estimated".}
#' \item{Index_distribution}{Survey/index biomass observation likelihood family. 0 or "Lognormal" (default) = independent lognormal on the log observation. 1 or "MVN" and 2 or "MVNORM" = multivariate normal on the natural-scale residual (obs - q*pred), with a user-supplied covariance in \code{index_cov} (e.g. a VAST-derived Sigma); "MVN" reports the bare quadratic form matching AMAK/ebswp and ADMB, "MVNORM" the properly normalized density, and the two give an identical fit. Pair either with Catchability = "AnalyticalArith". 3 or "Normal" = independent normal on the same residual, with an ABSOLUTE sd from Log_sd, not a log-scale CV; reproduces AMAK avo_like/cpue_like term for term. 4 or "TruncatedNormal" = that normal left-truncated at zero. Prefer it to "Normal" unless an exact ADMB comparison is needed: an index cannot be negative, so it is the only natural-scale family whose likelihood and simulator are the same distribution. See \code{vignette("model-options-and-functionality")}.}
#' \item{Comp_distribution}{Age/length composition data distribution. Accepts integer codes or readable strings: -1 or "MultinomialAFSC" (default) = AFSC multinomial; 0 or "Multinomial" = full multinomial; 1 or "DirichletMultinomial".}
#' \item{Comp_weights}{Composition weight. Its scale follows \code{Comp_distribution}: under a multinomial it is the natural-scale multiplier on the input sample size, while under a Dirichlet-multinomial the model uses \code{exp(Comp_weights)}, so the column holds the \emph{log} of the starting weight, a value of 1 is a starting weight of \eqn{e}; use 0 for a weight of 1. It is read when a model is built from scratch, so a refit keeps whatever weight it is given. A fit reports the weight its residuals imply, a harmonic mean of the ratio of effective to input sample size (McAllister & Ianelli 1997), in \code{Comp_weights_mcallister}; \code{\link{reweight_comps}} tunes towards it.}
#' \item{CAAL_distribution}{Conditional age-at-length composition distribution. Accepts integer codes or readable strings: -1 or "MultinomialAFSC" = AFSC multinomial; 0 or "Multinomial" (default) = full multinomial; 1 or "DirichletMultinomial".}
#' \item{CAAL_weights}{Conditional age-at-length (CAAL) composition weight. Its scale follows \code{CAAL_distribution} exactly as \code{Comp_weights} follows \code{Comp_distribution}: a natural-scale sample-size multiplier under a multinomial, the \emph{log} of the starting weight under a Dirichlet-multinomial.}
#' \item{Comp_accum_young}{Young-tail composition accumulation bin (AFSC ac_yng): a 1-based ordinal on the fleet's composition dimension; age/length bins below it are folded into it before the likelihood. NA (default) or 1 = no young accumulation.}
#' \item{Comp_accum_old}{Old-tail composition accumulation bin (AFSC ac_old): a 1-based ordinal on the fleet's composition dimension; age/length bins above it are folded into it before the likelihood. NA (default), 0, or a value at/above the number of bins = no old accumulation.}
#' \item{Month}{Observation month for the fleet (0 = not specified).}
#' \item{Proj_F_proportion}{The proportion of future fishing mortality assigned to this fleet.}
#'}
"BS2017SS"


#' Data inputs for multispecies CEATTLE of the Bering Sea from 1979 to 2017
#'
#' A data list containing inputs for the three-species (walleye pollock,
#' Pacific cod, arrowtooth flounder) multispecies CEATTLE model for the
#' Eastern Bering Sea. See \code{\link{BS2017SS}} for format details.
#'
#' @format A list with the same structure as \code{\link{BS2017SS}}.
"BS2017MS"


#' Fitted single-species CEATTLE model for the Eastern Bering Sea
#'
#' A fitted \code{Rceattle} model object for the single-species run of CEATTLE
#' for the Eastern Bering Sea (walleye pollock, Pacific cod, arrowtooth
#' flounder) without estimating natural mortality.
#'
#' @format An object of class \code{Rceattle}.
"EBS_ss_run"


#' Fitted single-species CEATTLE model with estimated M for the Eastern Bering Sea
#'
#' A fitted \code{Rceattle} model object for the single-species run of CEATTLE
#' for the Eastern Bering Sea with natural mortality estimated.
#'
#' @format An object of class \code{Rceattle}.
"EBS_ss_M_run"


#' Fitted multispecies CEATTLE model for the Eastern Bering Sea
#'
#' A fitted \code{Rceattle} model object for the multispecies run of CEATTLE
#' for the Eastern Bering Sea (walleye pollock, Pacific cod, arrowtooth
#' flounder).
#'
#' @format An object of class \code{Rceattle}.
"EBS_ms_run"


#' Data inputs for a single-species Gulf of Alaska CEATTLE model (2018)
#'
#' A data list containing inputs for a single-species CEATTLE model fit to
#' Gulf of Alaska groundfish data through 2018. See \code{\link{BS2017SS}}
#' for format details.
#'
#' @format A list with the same structure as \code{\link{BS2017SS}}.
"GOA2018SS"


#' Data inputs for Gulf of Alaska arrowtooth flounder CEATTLE model
#'
#' A data list containing inputs for a single-species CEATTLE model fit to
#' Gulf of Alaska arrowtooth flounder data. See \code{\link{BS2017SS}}
#' for format details.
#'
#' @format A list with the same structure as \code{\link{BS2017SS}}.
"GOAatf"


#' Data inputs for Gulf of Alaska arrowtooth flounder CEATTLE model (2023)
#'
#' A data list containing inputs for a single-species CEATTLE model fit to
#' Gulf of Alaska arrowtooth flounder data through 2023. See
#' \code{\link{BS2017SS}} for format details.
#'
#' @format A list with the same structure as \code{\link{BS2017SS}}.
"GOAatf2023"


#' Data inputs for Gulf of Alaska Pacific cod CEATTLE model
#'
#' A data list containing inputs for a single-species CEATTLE model fit to
#' Gulf of Alaska Pacific cod data. See \code{\link{BS2017SS}} for format
#' details.
#'
#' @format A list with the same structure as \code{\link{BS2017SS}}.
"GOAcod"


#' Data inputs for Gulf of Alaska walleye pollock CEATTLE model
#'
#' A data list containing inputs for a single-species CEATTLE model fit to
#' Gulf of Alaska walleye pollock data. See \code{\link{BS2017SS}} for
#' format details.
#'
#' @format A list with the same structure as \code{\link{BS2017SS}}.
"GOApollock"


#' Gulf of Alaska 2018 SAFE report reference values
#'
#' A list containing biomass, spawning stock biomass, and recruitment
#' reference values from the 2018 Gulf of Alaska Stock Assessment and
#' Fishery Evaluation (SAFE) report.
#'
#' @format A list with components:
#' \describe{
#'   \item{biomass}{Total biomass time series}
#'   \item{ssb}{Spawning stock biomass time series}
#'   \item{recruitment}{Recruitment time series}
#' }
"GOAsafe2018"


#' Data inputs for a three-species Georges Bank CEATTLE model
#'
#' A data list containing inputs for a multispecies CEATTLE model fit to
#' three groundfish species on Georges Bank. See \code{\link{BS2017SS}} for
#' format details.
#'
#' @format A list with the same structure as \code{\link{BS2017SS}}.
"GeorgesBank3spp"


#' Data inputs for Northern Rockfish CEATTLE model (2022)
#'
#' A data list containing inputs for a single-species CEATTLE model fit to
#' Gulf of Alaska northern rockfish data through 2022. See
#' \code{\link{BS2017SS}} for format details.
#'
#' @format A list with the same structure as \code{\link{BS2017SS}}.
"NorthernRockfish2022"


#' Data inputs for Atka mackerel CEATTLE model (2022)
#'
#' A data list containing inputs for a single-species CEATTLE model fit to
#' Aleutian Islands Atka mackerel data through 2022. See
#' \code{\link{BS2017SS}} for format details.
#'
#' @format A list with the same structure as \code{\link{BS2017SS}}.
"Atka2022"


#' Data inputs for CEATTLE model with WHAM-estimated growth
#'
#' A data list containing inputs for a CEATTLE model that uses growth
#' estimated from the Woods Hole Assessment Model (WHAM). See
#' \code{\link{BS2017SS}} for format details.
#'
#' @format A list with the same structure as \code{\link{BS2017SS}}.
"whamGrowthData"
