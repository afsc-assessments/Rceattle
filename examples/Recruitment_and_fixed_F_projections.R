#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Setup ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
library(Rceattle)
library(dplyr)
library(ggplot2)


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Data ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Example
# To run the 2017 single species assessment for the Bering Sea, a data file must first be loaded:
data(BS2017SS) # ?BS2017SS for more information on the data
BS2017SS$projyr <- 2060


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Estimation ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Then the model can be fit by setting `msmMode = 0` using the `Rceattle` function:
ss_run <- Rceattle::fit_mod(data_list = BS2017SS,
                            file = NULL,
                            inits = NULL, # Initial parameters = 0
                            estimateMode = 0, # Estimate
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            fit_control = fit_control(
                              phase   = TRUE, # Phase the model
                              verbose = 1))
# Type ?fit_mod for more details


# For the a multispecies model starting from the single species parameters, the following can be specified to load the data:
data("BS2017MS") # Note: the only difference is the residual mortality (M1_base) is lower
# Or we can use the previous data set
BS2017MS$projyr <- 2060

ms_run <- Rceattle::fit_mod(data_list = BS2017MS,
                            inits = ss_run$estimated_params, # Initial parameters from single species MLEs
                            file = NULL, # Don't save
                            estimateMode = 0, # Estimate
                            niter = 5, # 5 iterations around population and predation dynamics
                            random_rec = FALSE, # No random recruitment
                            msmMode = 1, # MSVPA based
                            suitMode = 0, # empirical suitability
                            fit_control = fit_control(
                              verbose = 1))


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Projections ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# * Changing future F ----
# Rceattle automatically projects the population forward when estimating under no fishing
BS2017MS$projyr # Year the population is projected forward
BS2017MS$fleet_control$Proj_F_proportion <- 1 # 1 fishing fleet per species

# Re-run, without estimating
ms_run_proj <- fit_mod(data_list = BS2017MS,
                       inits = ms_run$estimated_params, # Initial parameters from single species ests
                       file = NULL, # Don't save
                       estimateMode = 2, # Run in projection mode
                       HCR = build_hcr(HCR = 2,
                                       Ftarget = c(0.2342936, 0.513, 0.0774777)), # Set projection F mean historical F
                       niter = 5, # 5 iterations around population and predation dynamics
                       random_rec = FALSE, # No random recruitment
                       msmMode = 1, # MSVPA based
                       suitMode = 0, # empirical suitability
                       fit_control = fit_control(
                         verbose = 1))
plot_catch(ms_run_proj, incl_proj = T)




# * Future recruitment (by hand) ----
# Change recruitment deviations -

# Get years from projection
nyrs <- BS2017MS$endyr - BS2017MS$styr + 1
nyrs_proj <- BS2017MS$projyr - BS2017MS$styr + 1
yrs_proj <- (nyrs + 1):nyrs_proj

# Replace future rec_devs with numbers
ms_run$estimated_params$rec_dev[,yrs_proj] <- replace(
  ms_run$estimated_params$rec_dev[,yrs_proj],
  values = stats::rnorm( length(ms_run$estimated_params$rec_dev[,yrs_proj]),
                         mean = 0,
                         sd = 0.707) # Assumed value from penalized likelihood
)

ms_run_proj2 <- Rceattle::fit_mod(data_list = BS2017MS,
                                  inits = ms_run$estimated_params, # Initial parameters from single species ests
                                  file = NULL, # Don't save
                                  estimateMode = 3, # Do not estimate. Not changing parameters right now
                                  recFun = build_srr(proj_mean_rec = FALSE), # Project from rec_dev, not avg_R
                                  niter = 5, # 5 iterations around population and predation dynamics
                                  random_rec = FALSE, # No random recruitment
                                  msmMode = 1, # MSVPA based
                                  suitMode = 0, # empirical suitability
                                  fit_control = fit_control(
                                    verbose = 1))



# * Changing future recruitment (function) ----
# Function sample historical recruitment deviates and place in the projection
# "rec_trend" adds a linear trend in recruitment. In this case, a 50% decrease by 2026
ms_run_proj3 <- sample_rec(ms_run_proj2, sample_rec = TRUE, update_model = TRUE, rec_trend = -0.5) # fitted with proj_mean_rec = FALSE


# * Plot ----
mod_list <- list(ms_run, ms_run_proj, ms_run_proj2, ms_run_proj3)
mod_names <- c("MS no F", "MS mean historical F", "MS mean historical F w random rec", "MS mean historical F w trend")

# Plot biomass trajectory
plot_biomass(Rceattle = mod_list, model_names = mod_names, incl_proj = TRUE)
plot_recruitment(Rceattle = mod_list, model_names = mod_names, incl_proj = TRUE)
plot_catch(Rceattle = mod_list, model_names = mod_names, incl_proj = TRUE)



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Model variants ----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# For recruitment, the model can estimate recruitment deviates as random effects
ss_re <- Rceattle::fit_mod(
  data_list = BS2017SS,
  inits = NULL, # Initial parameters = 0
  file = NULL, # Don't save
  estimateMode = 0, # Estimate
  random_rec = TRUE, # Set recruitment deviations as random effects
  msmMode = 0, # Single species mode
  fit_control = fit_control(
    verbose = 1,
    phase = FALSE))

