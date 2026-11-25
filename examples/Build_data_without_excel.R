# =============================================================================
# Build a Rceattle data object in R (without Excel)
# =============================================================================
#
# RECOMMENDED PATH: build_data()
#   For most uses you do not need to assemble the list by hand. build_data()
#   takes only the blocks a model uses and default-fills the rest, checks block
#   names against the schema (catching typos early), and prints a spec tree:
#
#     dat <- build_data(base = BS2017SS, projyr = 2060)          # copy-and-edit
#     dat <- build_data(nspp = 1, styr = 1977, endyr = 2023,     # from blocks
#                       fleet_control = fc, weight = wt, catch_data = catch)
#     data_requirements(msmMode = 0)   # which blocks does my model need?
#
#   See vignette("data-without-excel") for the full build_data()/data_requirements()
#   walkthrough. The script below shows the underlying by-hand construction --
#   useful for simulating data or understanding exactly what each block contains.
#
# PURPOSE (by-hand construction):
#   Rceattle can read input data from an Excel workbook using read_data().
#   This script shows how to build the same data object directly in R as a
#   named list.  This is useful when:
#     - you want to simulate data for testing or MSE operating models
#     - you are programmatically generating inputs (e.g., multi-scenario loops)
#     - you do not have an Excel workbook available
#
# WORKFLOW:
#   1. Define model dimensions (species, ages, years, etc.)
#   2. Build each component of the data list
#   3. Pass the completed list to fit_mod()
#
# RELATED FUNCTIONS:
#   build_data()   -- constructor that wraps this by-hand assembly (recommended)
#   data_requirements() -- reports which blocks a configuration needs
#   read_data()    -- reads the Excel template into this same list format
#   write_data()   -- writes a data list back to Excel (see end of script)
#   data_check()   -- validates a data list before fitting
#   fit_mod()      -- fits the CEATTLE model given a data list
# =============================================================================

library(Rceattle)

# =============================================================================
# 1. Model dimensions
# =============================================================================
# Set the core dimensions that drive the size of all downstream arrays.

nspp     <- 2   # Number of species (or stocks)
nages    <- 10  # Maximum number of age classes (same for all species here)
nlengths <- 20  # Number of length bins (used for age-length keys / CAAL data)
nyrs     <- 30  # Number of hindcast years
years    <- 1:nyrs

# =============================================================================
# 2. Initialise the data list
# =============================================================================
simData <- list()

# =============================================================================
# 3. Species / model controls
# =============================================================================
# These scalars and vectors describe the biological structure of each species.

simData$nspp    <- nspp          # Number of species
simData$styr    <- 1             # First year of the hindcast
simData$endyr   <- nyrs          # Last year of the hindcast
simData$projyr  <- nyrs + 10     # Last year of the projection period

simData$spnames <- paste0("Species", 1:nspp)  # Character labels used in output

# nsex: number of sexes modelled per species (1 = combined, 2 = sex-structured)
simData$nsex <- rep(1, nspp)

# spawn_month: calendar month of spawning (6 = mid-year approximation)
simData$spawn_month <- rep(6, nspp)

# nages / minage: number of age bins and min-age  per species
simData$nages   <- rep(nages, nspp) # number, not max age
simData$minage  <- rep(1, nspp)

# nlengths: number of length bins per species
simData$nlengths <- rep(nlengths, nspp)

# estDynamics: 0 = estimate population dynamics (standard)
#              1 = fix numbers-at-age exactly to NByageFixed
#              2 = scale NByageFixed by a single (age-invariant) estimated scalar per species
#              3 = scale NByageFixed by age-specific estimated scalars per species
simData$estDynamics <- rep(0, nspp)

# pop_wt_index / ssb_wt_index: Wt_index value from simData$weight to use for
#   total biomass (pop) and spawning biomass (ssb) calculations.
#   Multiple weight datasets can exist (e.g. fleet-specific), each with a unique
#   Wt_index; these vectors select which one to use per species.
simData$pop_wt_index <- 1:nspp
simData$ssb_wt_index <- 1:nspp

# Length-weight relationship coefficients: W = alpha * L ^ beta (metric tonnes)
simData$alpha_wt_len <- rep(0.00001, nspp)
simData$beta_wt_len  <- rep(3,       nspp)

# pop_age_transition_index: row index into simData$age_trans_matrix used for
#   the age-to-length transition (growth matrix) for biomass calculations
simData$pop_age_transition_index <- rep(1, nspp)

# sigma_rec: prior standard deviation for log-recruitment deviations
simData$sigma_rec <- rep(1, nspp)

# other_food: non-modelled prey biomass available to each predator (used in
#   multi-species mode to prevent unrealistically high predation mortality)
simData$other_food <- rep(1e5, nspp)


# =============================================================================
# 4. Fleet control table
# =============================================================================
# One row per fleet. Columns control how each fleet is modelled.
#
# Fleet_type:    "Off" (= 0), "Fishery" (= 1), or "Survey" (= 2)
# Selectivity:   "Logistic", "NonParametric", "DoubleLogistic",
#                "DescendingLogistic", "Hake", "2DAR1", "3DAR1", or int 0-7
# Selectivity_dimension: "Age" or "Length"
# Catchability:  "Fixed" (0) = fixed at Catchability_init,
#                "Estimated" (1) = freely estimated,
#                "Estimated-with-prior" (2) = free w/ normal prior (Catchability_init, Catchability_prior_sd),
#                "Analytical" (3) = Ludwig & Walters/Martell 1994,
#                "Environmental" (5) = mu_q + X * beta,
#                "AR1" (6) = sensu Rogers et al. (2024),
#                NA = not applicable (fisheries)
# Proj_F_proportion:   Proportion of projected F per fleet; must sum to 1 per
#                species across fishery fleets; NA for surveys
# Comp_distribution:  "Multinomial" (0),  "DirichletMultinomial" (1),
#                "MultinomialAFSC" (-1)
# CAAL_distribution:  "Multinomial" (0) or "DirichletMultinomial" (1)

simData$fleet_control <- data.frame(
  Fleet_name      = paste0(c("Survey", "Fishery"),
                           rep(paste(" Species", 1:nspp), each = 2)),
  Fleet_code      = 1:(nspp * 2),          # Unique integer ID for each fleet
  Fleet_type      = rep(c("Survey", "Fishery"), nspp),
  Species         = rep(1:nspp, each = 2), # Species this fleet targets

  Month           = 0,                     # Month of sampling (0 = mid-year)

  # Selectivity settings
  Selectivity_index         = 1:(nspp * 2), # Links to selectivity across fleets
  Selectivity               = "Logistic",
  Selectivity_dimension     = "Age",
  N_sel_bins                = NA,           # Used for non-parametric only
  Sel_curve_pen1            = NA,           # Smoothness penalty for Ianelli selectivity
  Sel_curve_pen2            = NA,
  Time_varying_sel          = 0,            # 0 = time-invariant
  # 1 = penalized deviates (fixed at "Time_varying_sel_sd_prior" if random_sel = FALSE in fit_mod or estimated if TRUE)
  # 3 = time blocks (R-side mapping only)
  # 4 = random walk
  # 5 = random walk on ascending limb of double logistic only
  Time_varying_sel_sd_prior = 1,
  Bin_first_selected        = 1,            # Youngest fully-selected age/length
  Sel_norm_bin             = NA,           # Bin to normalize from (optional)
  Sel_norm_bin_upper             = NA,           # Upper bin to normalize from a range (can be NA)

  # Composition data settings
  Comp_distribution    = "Multinomial",
  Comp_weights    = 1,                      # Input effective sample size weight
  CAAL_distribution    = "Multinomial",
  CAAL_weights    = 1,

  # Data weighting / index units
  Observation_units = 1,                    # 1 = weight (mt), 2 = numbers
  Weight_index     = rep(1:nspp, each = 2),# Index of weight for biomass calc
  Age_transition_index = 1,                # Index of age_trans_matrix

  # Catchability settings
  Catchability_index              = 1:(nspp * 2),         # Links to q across fleets
  Catchability         = rep(c("Estimated", NA), nspp), # NA = not applicable
  Catchability_init              = rep(c(1,  NA), nspp),
  Catchability_prior_sd           = rep(c(0.2,NA), nspp),
  Time_varying_q       = rep(c(0,  NA), nspp), # 0/"None" = none, 1/"IID" = penalized lognormal deviates with sd set at "Time_varying_q_sd_prior" if random_q = FALSE in fit_mod
  # 3 /"Block" = time blocks, 4/ "RandomWalk"
  Time_varying_q_sd_prior = rep(c(1, NA), nspp),

  # Survey index uncertainty
  Estimate_index_sd    = rep(c(0,  NA), nspp), # 0 = fix at "Log_sd" from index_data,
  # 1 = estimate sd,
  # 2 = analytical sigma (Ludwig & Walters 1994)
  Index_sd_prior       = rep(c(1,  NA), nspp), # Starting value if estimated

  # Catch uncertainty
  Estimate_catch_sd    = rep(c(NA, 0), nspp),  # 0 = fix at "Log_sd" in catch_data, 1 = estimate
  Catch_sd_prior       = rep(c(NA, 1), nspp),  # Starting value if estimated

  # Projection fishing mortality
  Proj_F_proportion  = rep(c(NA, 1), nspp) # NA for surveys; must sum to 1 per species
)


# =============================================================================
# 5. Survey index data
# =============================================================================
# One row per fleet-year combination for survey fleets.
# Observation: log-scale biomass index (or numbers if Observation_units = 2)
# Log_sd:      observation standard deviation for lognormal

simData$index_data <- data.frame(
  Fleet_name       = rep(paste("Survey Species", 1:nspp), each = nyrs),
  Fleet_code       = rep(seq(1, nspp * 2, by = 2), each = nyrs),
  Species          = rep(1:nspp, each = nyrs),
  Year             = rep(years, nspp), # Negative year will predict but not include in likelihood
  Month            = 0,
  Selectivity_block = 1,               # Selectivity/q parameter block for this year (used
  # when Time_varying_sel or Time_varying_q = 3)
  Observation      = rnorm(nyrs * nspp),   # Replace with real log-index values
  Log_sd           = 0.1
)


# =============================================================================
# 6. Catch data
# =============================================================================
# One row per fishery fleet-year combination.
# Catch:   total catch in metric tonnes (or numbers if Observation_units = 2)
# Log_sd:  observation sd for lognormal (often set small, e.g. 0.01-0.05,
#          when catch is treated as known)

simData$catch_data <- data.frame(
  Fleet_name       = rep(paste("Fishery Species", 1:nspp), each = nyrs),
  Fleet_code       = rep(seq(2, nspp * 2, by = 2), each = nyrs),
  Species          = rep(1:nspp, each = nyrs),
  Year             = rep(years, nspp),
  Month            = 0,
  Selectivity_block = 1,
  Catch            = abs(rnorm(nyrs * nspp)), # Replace with real catch values
  Log_sd           = 0.05
)


# =============================================================================
# 7. Composition data (age or length) [optional]
# =============================================================================
# One row per observation. Composition columns (Comp_1, Comp_2, ...) hold
# observed proportions OR raw counts -- Rceattle normalises internally.
#
# Age0_Length1: 0 = age composition,  1 = length composition
# Sex:          0 = combined,  1 = females,  2 = males, 3 = joint
# Sample_size:  input sample size

comp_matrix <- matrix(abs(rnorm(nyrs * nspp * nages * 2)),
                      nrow = nyrs * nspp * 2, ncol = nages)
colnames(comp_matrix) <- paste0("Comp_", 1:nages)

simData$comp_data <- cbind(
  data.frame(
    Fleet_name  = c(rep(paste("Survey Species",  1:nspp), each = nyrs),
                    rep(paste("Fishery Species", 1:nspp), each = nyrs)),
    Fleet_code  = rep(c(seq(1, nspp * 2, 2), seq(2, nspp * 2, 2)), each = nyrs),
    Species     = c(rep(1:nspp, each = nyrs), rep(1:nspp, each = nyrs)),
    Sex         = 0,
    Age0_Length1 = 0,            # 0 = age comp
    Year        = rep(years, 2 * nspp), # Negative year will predict but not include in likelihood
    Month       = 0,
    Sample_size = 200
  ),
  comp_matrix
)


# =============================================================================
# 8. Conditional age-at-length (CAAL) data  [optional]
# =============================================================================
# CAAL provides age observations within each length bin; this allows the model
# to jointly fit length and age compositions. Leave empty if not available.
#
# One row per fleet-sex-year-length-bin combination.
# CAAL_1, CAAL_2, ... hold age frequencies within that length bin.
# Sex:  0 = combined,  1 = females,  2 = males

caal_matrix <- matrix(abs(rnorm(nyrs * nspp * nages * 2 * nlengths)),
                      nrow = nyrs * nspp * 2 * nlengths, ncol = nages)
colnames(caal_matrix) <- paste0("CAAL_", 1:nages)

simData$caal_data <- cbind(
  data.frame(
    Fleet_name  = c(rep(paste("Survey Species",  1:nspp), each = nyrs * nlengths),
                    rep(paste("Fishery Species", 1:nspp), each = nyrs * nlengths)),
    Fleet_code  = rep(c(seq(1, nspp * 2, 2), seq(2, nspp * 2, 2)),
                      each = nyrs * nlengths),
    Species     = c(rep(1:nspp, each = nyrs * nlengths),
                    rep(1:nspp, each = nyrs * nlengths)),
    Sex         = 0,
    # Negative year will predict but not include in likelihood
    Year        = rep(years,      2 * nspp * nlengths),
    Length      = rep(1:nlengths, 2 * nspp * nyrs),
    Sample_size = 200
  ),
  caal_matrix
)


# =============================================================================
# 9. Empirical selectivity  [optional -- leave empty if not used]
# =============================================================================
# Pre-specified selectivity curves that are fixed (not estimated).
# If provided, fixed selectivity-at-age values.

simData$emp_sel <- data.frame(
  matrix(NA, nrow = 0, ncol = 5 + nages)
)
colnames(simData$emp_sel) <- c("Fleet_name", "Fleet_code", "Species",
                               "Sex", "Year", paste0("Comp_", 1:nages))


# =============================================================================
# 10. Fixed numbers-at-age  [optional -- leave empty unless estDynamics = 1, 2, or 3]
# =============================================================================
# Used when estDynamics > 0:
#   1 = numbers-at-age fixed exactly to these values
#   2 = numbers-at-age scaled by a single estimated scalar per species
#   3 = numbers-at-age scaled by age-specific estimated scalars per species

simData$NByageFixed <- data.frame(
  matrix(NA, nrow = 0, ncol = 4 + nages)
)
colnames(simData$NByageFixed) <- c("Species_name", "Species", "Sex", "Year",
                                   paste0("Age", 1:nages))


# =============================================================================
# 11. Age-to-length transition matrix (growth matrix)
# =============================================================================
# Probability of an individual of age a being observed in length bin l.
# Used when fitting length compositions and for biomass-at-length calculations.
# Each row is an age; Length_1 ... Length_nlengths columns sum to 1 across row.
#
# Here we use a diagonal identity (each age maps to a unique length bin) as a
# placeholder. Replace with empirically derived or model-based growth matrices.

age_trans_list <- lapply(1:nspp, function(sp) {
  tmp <- matrix(0, nrow = nages, ncol = nlengths)
  diag(tmp[1:min(nages, nlengths), 1:min(nages, nlengths)]) <- 1
  colnames(tmp) <- paste0("Length_", 1:nlengths)
  cbind(
    data.frame(
      Age_transition_name  = paste("Species", sp),
      Age_transition_index = sp,
      Species              = sp,
      Sex                  = 0,
      Age                  = 1:nages
    ),
    tmp
  )
})
simData$age_trans_matrix <- do.call(rbind, age_trans_list)


# =============================================================================
# 12. Age-reading error matrix
# =============================================================================
# Probability that a fish of true age t is read as observed age o.
# An identity matrix (no ageing error) is a reasonable starting point.
# Rows = true ages, columns (Obs_age1 ... Obs_ageN) = observed ages.

age_error_list <- lapply(1:nspp, function(sp) {
  tmp <- as.data.frame(diag(1, nages))
  colnames(tmp) <- paste0("Obs_age", 1:nages)
  cbind(
    data.frame(Species  = sp,
               True_age = 1:nages),
    tmp
  )
})
simData$age_error <- do.call(rbind, age_error_list)


# =============================================================================
# 13. Weight-at-age
# =============================================================================
# Used for biomass calculations. Linked to species or fleets.
# Set Year = 0 to apply the same weight-at-age to all years (time-invariant).
# Age_1 ... Age_nages columns hold mean weight (usually kg) at each age.

weight_matrix <- matrix(
  (1:nages / nages)^3 * 0.01,   # Simple power curve placeholder
  nrow = nspp, ncol = nages, byrow = TRUE
)
colnames(weight_matrix) <- paste0("Age", 1:nages)

simData$weight <- cbind(
  data.frame(
    Wt_name  = paste("Species", 1:nspp),
    Wt_index = 1:nspp,
    Species  = 1:nspp,
    Sex      = 0,
    Year     = 0    # Year = 0 applies to all years (time-invariant)
  ),
  weight_matrix
)


# =============================================================================
# 14. Maturity-at-age
# =============================================================================
# Proportion mature at each age, used to compute spawning stock biomass.
# Age_1 ... Age_nages columns; one row per species.

mat_matrix <- matrix(1, nrow = nspp, ncol = nages)  # All ages mature (placeholder)
colnames(mat_matrix) <- paste0("Age", 1:nages)

simData$maturity <- cbind(
  data.frame(Species = 1:nspp),
  mat_matrix
)


# =============================================================================
# 15. Sex ratio at age
# =============================================================================
# Proportion female at each age. Used when nsex = 1 (combined-sex model)
# Proportion female at recruitment. Used when nsex = 2 (two-sex model)
# to partition biomass for SSB calculations.

sexratio_matrix <- matrix(0.5, nrow = nspp, ncol = nages)
colnames(sexratio_matrix) <- paste0("Age", 1:nages)

simData$sex_ratio <- cbind(
  data.frame(Species = 1:nspp),
  sexratio_matrix
)


# =============================================================================
# 16. Natural mortality (M1 base)
# =============================================================================
# Baseline natural mortality-at-age before any predation component (M2).
# In single-species mode, total mortality M = M1.
# In multi-species mode, M = M1 + M2 (predation).

m_matrix <- matrix(0.2, nrow = nspp, ncol = nages)   # Constant M = 0.2 placeholder
colnames(m_matrix) <- paste0("Age", 1:nages)

simData$M1_base <- cbind(
  data.frame(Species = 1:nspp,
             Sex     = 0),
  m_matrix
)


# =============================================================================
# 17. Environmental covariate data
# =============================================================================
# Optional time series used in environment-recruitment,-M1, or -Q
# relationships.  Missing years will be filled in with the mean.
# EnvData columns can be named anything; column names are referenced by column order.

simData$env_data <- data.frame(
  Year    = 1:nyrs,
  EnvData = 1       # Constant = no environmental effect (placeholder)
)


# =============================================================================
# 18. Ration data  [optional -- leave empty if not used]
# =============================================================================
# Observed mean ration-at-age or relative foragine rate from bioenergetics studies.
# Used only in multi-species mode when ration is empirically constrained.

simData$ration_data <- data.frame(
  matrix(NA, nrow = 0, ncol = 3 + nages)
)
colnames(simData$ration_data) <- c("Species", "Sex", "Year",
                                   paste0("Age", 1:nages))


# =============================================================================
# 19. Bioenergetics parameters  [multi-species mode only]
# =============================================================================
# Parameters for the Wisconsin bioenergetics model used to scale predation
# mortality in multi-species mode.  In single-species mode these are unused.
# See ?build_M1 and the CEATTLE technical documentation for details.

simData$Ceq   <- rep(4,  nspp)  # Consumption equation form (Hanson et al. 1997)
simData$Cindex <- rep(1, nspp)  # Temperature index type
simData$Pvalue <- rep(1, nspp)  # Proportion of maximum consumption realised
simData$fday   <- rep(1, nspp)  # Foraging days per year

# Wisconsin model coefficients (CA, CB) and temperature dependence (Qc, Tco, etc.)
simData$CA  <- rep(1,  nspp)
simData$CB  <- rep(-1, nspp)
simData$Qc  <- rep(1,  nspp)
simData$Tco <- rep(1,  nspp)
simData$Tcm <- rep(1,  nspp)
simData$Tcl <- rep(1,  nspp)
simData$CK1 <- rep(1,  nspp)
simData$CK4 <- rep(1,  nspp)

simData$Diet_distribution      <- rep(0, nspp)  # 0 = Multinomial, 1 = Dirichlet-multinomial
# (only used when suitMode > 0)
simData$Diet_comp_weights <- rep(1, nspp)


# =============================================================================
# 20. Diet composition data  [multi-species mode only -- leave empty if unused]
# =============================================================================
# Observed stomach content data. One row per predator-prey-year-age combination.
# Stomach_proportion_by_weight: proportion of prey in predator's diet by weight.

simData$diet_data <- data.frame(
  matrix(NA, nrow = 0, ncol = 9)
)
colnames(simData$diet_data) <- c("Pred", "Prey", "Pred_sex", "Prey_sex",
                                 "Pred_age", "Prey_age", "Year",
                                 "Sample_size", "Stomach_proportion_by_weight")


# =============================================================================
# 21. (Optional) Validate the data list before fitting
# =============================================================================
# data_check() prints a summary of what it finds in the data list and flags
# common issues (mismatched dimensions, missing columns, etc.)

# data_check(simData)


# =============================================================================
# 22. Fit the model
# =============================================================================
# Key fit_mod() switches (all have defaults; shown here for documentation):
#
# estimateMode: 0 = full fit (hindcast + projection)
#               1 = hindcast only
#               2 = projection only
#               3 = debug (MakeADFun only, skips optimisation)
#               4 = debug + nlminb
#
# msmMode:      0 = single-species (no predation mortality)
#               1 = Type II MSVPA (sensu Holsman et al. 2016)
#               2 = Type III MSVPA
#
# suitMode:     Predator-prey suitability (one value per predator species):
#               0 = empirical (MSVPA-style from diet data)
#               2 = gamma selectivity on weight
#               4 = lognormal on weight
#               5 = normal on length
#               6 = normal on weight
#
# initMode:     0 = free parameters for initial N-at-age
#               1 = equilibrium, Finit = 0, no init devs (unfished)
#               2 = equilibrium, Finit = 0, with init devs (default)
#               3 = non-equilibrium: Finit estimated, init devs included
#               4 = non-equilibrium: Finit scales R0
#
# avgnMode:     Abundance used in predation mortality calculations:
#               0 = N/Z * (1 - exp(-Z))   (default)
#               1 = N * exp(-Z/2)
#               2 = N
#
# random_rec:   TRUE = treat log-recruitment deviations as random effects
#               (marginalized via Laplace approximation); FALSE = penalized fixed effects
# random_sel:   TRUE = treat selectivity deviations as random effects
# random_q:     TRUE = treat catchability deviations as random effects
#
# phase:        TRUE = use parameter phasing to improve convergence (recommended)

sim_run <- Rceattle::fit_mod(
  data_list    = simData,
  estimateMode = 3,       # Change to 0 for a real fit
  msmMode      = 0,       # Single-species
  initMode     = 2,       # Equilibrium with init devs (default)
  random_rec   = FALSE,
  fit_control = fit_control(
    phase        = TRUE,
    verbose      = 1))


# =============================================================================
# 23. (Optional) Export back to Excel
# =============================================================================
# write_data() writes the data list to an Excel workbook in the standard
# Rceattle template format, which can then be edited and re-read with read_data().

# write_data(simData, file = "mydata.xlsx")
