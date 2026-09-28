# Comparison of Rceattle and WHAM fit to conditional age-at-length (CAAL) data
# Modified from https://giancarlomcorrea.netlify.app/labs/OFI_WK_2023/examples/case3.R
# Note that the plus group in year 1 is wrong on WHAM line 1720 of the cpp and needs adjustment

# remotes::install_github(repo = 'GiancarloMCorrea/wham', ref='growth', INSTALL_opts = c("--no-docs", "--no-multiarch", "--no-demo"))
library(readxl)
library(wham)
library(ggplot2)
library(dplyr)
source("tests/comparison/WHAM_growth_comparison/helper.R")
dir.create("tests/comparison/WHAM_growth_comparison")
runmodels = TRUE

# WHAM ------------------------------------------------------------------
catch_df = readxl::read_xlsx(path = 'tests/comparison/WHAM_growth_comparison/case3.xlsx', sheet = 1)
index_df = readxl::read_xlsx(path = 'tests/comparison/WHAM_growth_comparison/case3.xlsx', sheet = 2)
fsh_lcomp_df = readxl::read_xlsx(path = 'tests/comparison/WHAM_growth_comparison/case3.xlsx', sheet = 3)
index_lcomp_df = readxl::read_xlsx(path = 'tests/comparison/WHAM_growth_comparison/case3.xlsx', sheet = 4)
caal_df = readxl::read_xlsx(path = 'tests/comparison/WHAM_growth_comparison/case3.xlsx', sheet = 5)
maturity_df = readxl::read_xlsx(path = 'tests/comparison/WHAM_growth_comparison/case3.xlsx', sheet = 6)

# * Make input data: ----
input_data = list()
input_data$ages = 1:10 # ages
input_data$years = 1976:2020 # years
input_data$lengths = seq(from = 2, to = 130, by = 2) # lengths
input_data$n_fleets = 1 # number of fleets
input_data$n_indices = 1 # number of indices
n_years = length(input_data$years)
n_ages = length(input_data$ages)
n_lengths = length(input_data$lengths)

# Agg catch:
input_data$agg_catch = matrix(catch_df$catch, ncol = input_data$n_fleets, nrow = n_years) # Obs
input_data$catch_cv = matrix(catch_df$cv, ncol = input_data$n_fleets, nrow = n_years) # Obs error

# Length comps (fishery)
input_data$catch_pal = array(as.matrix(fsh_lcomp_df[,3:67]), dim = c(input_data$n_fleets, n_years, n_lengths)) # Obs
input_data$catch_NeffL = matrix(fsh_lcomp_df$Nsamp, ncol = input_data$n_fleets, nrow = n_years) # Obs error
input_data$use_catch_pal = matrix(1, nrow = n_years, ncol = input_data$n_fleets) # 1 = fit, 0 = don't fit

# Agg index:
input_data$agg_indices = matrix(index_df$index, ncol = input_data$n_indices, nrow = n_years) # Obs
input_data$index_cv = matrix(index_df$cv, ncol = input_data$n_indices, nrow = n_years) # Obs error

# Additional information:
input_data$units_indices = matrix(1L, nrow = n_years, ncol = input_data$n_indices) # 0 = numbers, 1 = biomass
input_data$fracyr_indices = matrix(0, ncol = input_data$n_indices, nrow = n_years) # fraction of the year when survey occurs

# Length comps (index):
input_data$index_pal = array(as.matrix(index_lcomp_df[,3:67]), dim = c(input_data$n_indices, n_years, n_lengths)) # Obs
input_data$index_NeffL = matrix(index_lcomp_df$Nsamp, ncol = input_data$n_indices, nrow = n_years) # Obs error
input_data$use_index_pal = matrix(1, nrow = n_years, ncol = input_data$n_indices) # 1 = fit, 0 = don't fit

# CAAL index:
input_data$index_caal = array(0, dim = c(input_data$n_indices, n_years, n_lengths, n_ages))
for(i in seq_along(input_data$years)){
  tmp = caal_df |> filter(year == input_data$years[i])
  input_data$index_caal[1,i,,] = as.matrix(tmp[,4:13])
}
# CAAL Neff index:
input_data$index_caal_Neff = array(0, dim = c(n_years, input_data$n_indices, n_lengths))
for(i in seq_along(input_data$years)){
  tmp = caal_df |> filter(year == input_data$years[i])
  input_data$index_caal_Neff[i,1,] = as.matrix(tmp[,3])
}
# CAAL index use/not use (use only when Nsamp > 0):
input_data$use_index_caal = array(0, dim = c(n_years, input_data$n_indices, n_lengths))
input_data$use_index_caal[input_data$index_caal_Neff > 0] = 1
# Dont forget to turn off the use of paa (default):
input_data$use_catch_paa = matrix(0, nrow = n_years, ncol = input_data$n_fleets) # 1 = fit, 0 = don't fit
input_data$use_index_paa = matrix(0, nrow = n_years, ncol = input_data$n_indices) # 1 = fit, 0 = don't fit
# Selex pointers:
input_data$selblock_pointer_fleets = matrix(1L, ncol = input_data$n_fleets, nrow = n_years)
input_data$selblock_pointer_indices = matrix(2L, ncol = input_data$n_indices, nrow = n_years)
# weight-at-age pointers information:
input_data$waa_pointer_fleets = 2
input_data$waa_pointer_indices = 1
input_data$waa_pointer_totcatch = 2
input_data$waa_pointer_ssb = 1
input_data$waa_pointer_jan1 = 1
# More information:
input_data$maturity = as.matrix(maturity_df[,2:11]) # maturity
input_data$fracyr_SSB = matrix(0, ncol = 1, nrow = n_years) # spawning fraction (0 = spawn at beginning of year)
input_data$Fbar_ages = 1:10 # es to include in mean F calculation
input_data$bias_correct_process = 1 # do process bias correction, 0 = no, 1 = yes
input_data$bias_correct_observation = 1 # do obs bias correction, 0 = no, 1 = yes

# * Create WHAM input ----
# - (LAA parametric approach)
input3 = wham::prepare_wham_input(model_name = "Case_3",
                                  basic_info = input_data,
                                  NAA_re = list(N1_model = 1, N1_pars = c(1e+05, 0),
                                                recruit_model = 2, recruit_pars = 1e+05,
                                                sigma = 'rec', cor = 'iid'), # Recruitment parameters
                                  M = list(model = 'constant', initial_means = 0.35), # M parameter
                                  selectivity = list(model = c('len-logistic', 'len-logistic'),
                                                     initial_pars = list(c(15, 3), c(15, 3)),
                                                     n_selblocks = 2, fix_pars = list(NULL, NULL) # GRANT: converted both to length-logistic
                                  ), # Selectivity parameter
                                  catchability = list(initial_q = 1), # Catchability parameter
                                  growth = list(model = 'vB_classic',
                                                init_vals = c(0.2, 90, 10), est_pars = 1:3,
                                                SD_vals = c(1, 9)),
                                  LW = list(init_vals = c(5.56e-06, 3.2))
)

# Fix some parameters:
input3$map$log_N1_pars = factor(c(1, NA))
input3$map$logit_q = factor(NA)
# Optional: fix sigmaR and turn off NAA (recruitment) as random variable (to make it faster):
# input3$par$log_NAA_sigma = log(0.6)
# input3$map$log_NAA_sigma = factor(NA)
# input3$random = NULL

# * Run model ----
if(isTRUE(runmodels)) {
  wham_model = wham::fit_wham(MakeADFun.silent = TRUE, input = input3, do.retro = FALSE, do.osa = FALSE)
  saveRDS(wham_model, file = "tests/comparison/WHAM_growth_comparison/wham_model.RDS")
} else {
  wham_model <- readRDS(file = "tests/comparison/WHAM_growth_comparison/wham_model.RDS")
}

wham_model$opt
wham_model$sdrep
wham_model$rep
names(wham_model)

# * Make plots ----
# plot_wham_output(wham_model, dir.main = "tests/comparison/WHAM_growth_comparison", out.type = 'pdf')






# Rceattle -------------------------------------------------------------
library(Rceattle)
data("GOAcod")
simData <- GOAcod


# * Make input data ----
# * Data controls
simData$spnames <- "WHAM data"
simData$nspp <- 1
years = 1976:2020 # years
nyrs <- length(years)
simData$styr <- 1976
simData$endyr <- 2020
simData$projyr <- 2030
simData$nsex <- 1
simData$nages <- 10
simData$minage <- 1
lengths = seq(from = 2, to = 130, by = 2) # lengths
simData$nlengths <- length(lengths)
simData$pop_wt_index <- 1
simData$ssb_wt_index <- 1
simData$spawn_month <- 0
simData$alpha_wt_len <- 5.56e-06
simData$beta_wt_len <- 3.2
simData$pop_age_transition_index <- 1
simData$sigma_rec <- 0.6
nages <- simData$nages
nlengths <- simData$nlengths

# * Fleet control
simData$fleet_control <- simData$fleet_control[c(1,3),] # BT and Trawl Fishery are both simple logistic
simData$fleet_control$Fleet_name <- c("Survey", "Fishery")
simData$fleet_control$Fleet_code <- 1:2
simData$fleet_control$Selectivity <- "Logistic"                  # Length-based logistic
simData$fleet_control$Selectivity_dimension <- "Length"          # Length-based logistic
simData$fleet_control$Selectivity_index <- 1:2
simData$fleet_control$Weight_index <- 1
simData$fleet_control$Observation_units <- 1
simData$fleet_control$Month <- c(0, 6)
simData$fleet_control$Bin_first_selected <- 1
simData$fleet_control$Catchability <- c(0, NA)
simData$fleet_control$Catchability_init <- c(1, NA)

# * Index data
simData$index_data <- data.frame(Fleet_name = "Survey",
                                 Fleet_code = 1,
                                 Species = 1,
                                 Year = years,
                                 Month = index_df$fr_yr,
                                 Selectivity_block = 1,
                                 Q_block = 1,
                                 Observation = index_df$index,
                                 Log_sd = sqrt(log(index_df$cv^2 + 1))) # WHAM converts CV -> sigma

# * Catch data
simData$catch_data <- data.frame(Fleet_name = "Fishery",
                                 Fleet_code = 2,
                                 Species = 1,
                                 Year = years,
                                 Month = 0,
                                 Selectivity_block = 1,
                                 Catch = catch_df$catch,
                                 Log_sd = sqrt(log(catch_df$cv^2 + 1))) # WHAM converts CV -> sigma

# * Comp data
# - Index length comp
tmp <- as.matrix(index_lcomp_df[,3:67])
colnames(tmp) <- paste0("Comp_",1:simData$nlengths)
index_comp <- cbind(data.frame(Fleet_name = "Survey",
                               Fleet_code = 1,
                               Species = 1,
                               Sex = 0,
                               Age0_Length1 = 1,
                               Year = years, # Negative turns data off
                               Month = 0,
                               Sample_size = index_lcomp_df$Nsamp),
                    tmp
)

# - Fishery length comp
tmp <- as.matrix(fsh_lcomp_df[,3:67])
colnames(tmp) <- paste0("Comp_",1:simData$nlengths)
fishery_comp <- cbind(data.frame(Fleet_name = "Fishery",
                                 Fleet_code = 2,
                                 Species = 1,
                                 Sex = 0,
                                 Age0_Length1 = 1,
                                 Year = years, # Negative turns data off
                                 Month = 0,
                                 Sample_size = fsh_lcomp_df$Nsamp),
                      tmp
)

simData$comp_data <- rbind(index_comp, fishery_comp)


# * CAAL data
# - Index
index_caal <- caal_df
colnames(index_caal) <- c("Year", "Length", "Sample_size", paste0("CAAL_",1:simData$nages))
simData$caal_data <- index_caal |>
  dplyr::mutate(Fleet_name = "Survey",
                Fleet_code = 1,
                Species = 1,
                Sex = 0) |>
  select(Fleet_name, Fleet_code, Species, Sex, Year, Length, Sample_size, paste0("CAAL_",1:simData$nages)) |>
  as.data.frame()
# - No fishery CAAL

# * Empirical selectivity
simData$emp_sel[] <- NA

# * Fixed numbers
simData$NByageFixed[] <- NA


# * Age transition matrix
tmp <- as.data.frame(diag(1,simData$nlengths))[1:simData$nages,]
colnames(tmp) <- paste0("Length_",1:simData$nlengths)
simData$age_trans_matrix <- cbind(data.frame(Age_transition_name = "Base",
                                             Age_transition_index = 1,
                                             Species = 1,
                                             Sex = 0,
                                             Age = 1:simData$nages),
                                  tmp
)


# * Age error
tmp <- as.data.frame(diag(1,simData$nages))
colnames(tmp) <- paste0("Obs_age",1:simData$nages)
simData$age_error <- cbind(data.frame(Species = 1,
                                      True_age = 1:simData$nages),
                           tmp
)

# * Empirical weight-at-age
WAA <- as.data.frame(matrix(1:simData$nages, ncol = simData$nages))
colnames(WAA) <- paste0("Age",1:simData$nages)
simData$weight <- cbind(data.frame(Wt_name = "Base",
                                   Wt_index = 1,
                                   Species = 1,
                                   Sex = 0,
                                   Year = 0), # 0 fills all years
                        WAA
)

# * Maturity (time-invariant in Rceattle; the WHAM input repeats one ogive over
# all 45 years, so row 1 is the whole schedule). as.matrix() first: a tibble row
# is a list, and matrix() would build a list-matrix that reads back as NA.
MatAA <- matrix(as.matrix(maturity_df[1,2:11]), ncol = simData$nages)
colnames(MatAA) <- paste0("Age",1:simData$nages)
simData$maturity <- cbind(data.frame(Species = 1),
                          MatAA
) |>
  as.data.frame()

# * Sex ratio
sexratio <- as.data.frame(matrix(0.5, ncol = simData$nages))
colnames(sexratio) <- paste0("Age",1:simData$nages)
simData$sex_ratio <- cbind(data.frame(Species = 1),
                           sexratio
)

# * Input mortality-at-age
mort <- as.data.frame(matrix(0.35, ncol = simData$nages))
colnames(mort) <- paste0("Age",1:simData$nages)
simData$M1_base <- cbind(data.frame(Species = 1,
                                    Sex = 0),
                         mort
)

# * Environmental data
simData$env_data <- data.frame(Year = years,
                               EnvData = 1)


# * Ration data
pyrs <- as.data.frame(matrix(1, nrow = 1, ncol = simData$nages))
colnames(pyrs) <- paste0("Age",1:simData$nages)
simData$ration_data <- cbind(data.frame(Species = 1,
                                        Sex = 0,
                                        Year = 0), # Zero fills out all years
                             pyrs
)


# Growth fn matched to WHAM: fix the growth SD (WHAM fixes SDgrowth_par) via the
# sd_L1 / sd_Linf linkage endpoints with est_phase = 0 (fix at init).  Init values
# are log-SD endpoints, threaded onto growth_log_sd (= WHAM SDgrowth_par = log(c(1,9))).
growthFun_vb <- build_growth(fun = "vonBertalanffy",
                             linkages = list(
                               sd_L1   = linkage_spec(~1, init = list(`(Intercept)` = 1), est_phase = 0L),
                               sd_Linf = linkage_spec(~1, init = list(`(Intercept)` = 9), est_phase = 0L)))

# * Null Rceattle ----
ss_fix <- Rceattle::fit_mod(data_list = simData,
                            inits = NULL, # Initial parameters = 0
                            estimateMode = 3, # Don't estimate
                            growthFun = growthFun_vb, # Von Bert, growth SD fixed to match WHAM
                            random_rec = FALSE, # No random recruitment
                            msmMode = 0, # Single species mode
                            phase = FALSE,
                            verbose = 1)


# * Fix-parameters and check derived quantities ----
inits <- ss_fix$estimated_params # get paramter object from NULL model
names(wham_model)
names(wham_model$parList)

# - Rec
inits$rec_pars[1,1] <- wham_model$parList$mean_rec_pars
inits$rec_dev[1,1] <- wham_model$parList$log_N1_pars[1] -  wham_model$parList$mean_rec_pars
inits$rec_dev[1,2:nyrs] <- wham_model$parList$log_NAA[,1] -  wham_model$parList$mean_rec_pars
inits$init_dev[1,] <- 0 # equilibrium initial age structure (ages 2+ at R0); year-1 recruitment enters via rec_dev[1,1]. Matches WHAM N1_model=1 (equilibrium-at-mean-rec) with penalized year-1 recruitment.
inits$R_log_sd <- wham_model$parList$log_NAA_sigma

# - F (random walk in WHAM)
inits$log_F[2,1]  <- wham_model$parList$log_F1
for(y in 2:nyrs){
  inits$log_F[2,y] <- inits$log_F[2,y-1] + wham_model$parList$F_devs[y-1,1]
}
inits$log_Finit[1] <- -Inf

# - Selectivity
selpars <-  wham_model$env$data$selpars_lower + ( wham_model$env$data$selpars_upper -  wham_model$env$data$selpars_lower) / (1.0 + exp(-(wham_model$parList$logit_selpars)))
inits$log_sel_slp[1,,1] <- rev(log(1/selpars[,24]))
inits$sel_inf[1,,1] <- rev(selpars[,23])

# - Q
inits$index_log_q[1] <- log(wham_model$env$data$q_lower + (wham_model$env$data$q_upper - wham_model$env$data$q_lower) / (1 + exp(-wham_model$parList$logit_q)))

# - Growth
inits$log_growth_pars[1,1,1:3] <- wham_model$parList$growth_a[c(1,3,2),1]
inits$growth_log_sd[1,1,] <- wham_model$parList$SDgrowth_par


ss_inits <- Rceattle::fit_mod(data_list = simData,
                              inits = inits, # Initial parameters = 0
                              estimateMode = 3, # Do not estimate
                              growthFun = growthFun_vb, # Von Bert, growth SD fixed to match WHAM
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              initMode = 1,
                              # comp_offset = 0 -> WHAM-style multinomial (no proportion offset)
                              fit_control = fit_control(phase = FALSE, verbose = 0, comp_offset = 0))


# SSB
plot(x = ss_inits$quantities$ssb[,1:nyrs] * 2, y = wham_model$rep$SSB); abline(0,1)
round(ss_inits$quantities$ssb[,1:nyrs] * 2 - wham_model$rep$SSB, 5) # Plus group in WHAM has a bug

# Rec
plot(x = ss_inits$quantities$N_at_age[1,1,1,1:nyrs], y = wham_model$rep$NAA[,1]); abline(0,1)
round(ss_inits$quantities$N_at_age[1,1,,1] - t(wham_model$rep$NAA)[,1], 6)

# N-at-age
round(ss_inits$quantities$N_at_age[1,1,,1:nyrs] - t(wham_model$rep$NAA), 6)

# Selex matches
plot(ss_inits$quantities$sel_at_length[1,1,,1], y = wham_model$rep$selLL[[2]][1,], xlab = "Rceattle survey sel"); abline(0,1)
plot(ss_inits$quantities$sel_at_length[2,1,,1], y = wham_model$rep$selLL[[1]][1,], xlab = "Rceattle fishery sel"); abline(0,1)
plot(ss_inits$quantities$sel_at_length[2,1,,10], y = wham_model$rep$selLL[[1]][10,], xlab = "Rceattle fishery sel"); abline(0,1)

# F
plot(x = ss_inits$quantities$F_spp[1,1:nyrs], y = wham_model$rep$F[,1]); abline(0,1)
ss_inits$quantities$F_spp[1,1:nyrs] - wham_model$rep$F[,1]

# F-at-age (good enough)
check <- outer(wham_model$rep$F[,1], t(wham_model$rep$catch_phi_mat[,,1]) %*% wham_model$rep$selLL[[1]][11,])
round(wham_model$rep$FAA_tot - check[,,1], 10)
round(t(wham_model$rep$FAA_tot) - ss_inits$quantities$F_spp_at_age[1,1,,1:nyrs], 10)

# ZAA (good enough)
round(t(wham_model$rep$ZAA) - ss_inits$quantities$Z_at_age[1,1,,1:nyrs], 10)

# M
plot(x = ss_inits$quantities$M1_at_age[1,1,1,1:nyrs], y = wham_model$rep$MAA[,1]); abline(0,1)

# Laa
plot(y = wham_model$rep$LAA[1,], x = t(ss_inits$quantities$length_hat[1,1,,1])); abline(0,1)

# Phi matrix
t(wham_model$rep$ssb_phi_mat[,,1]) - ss_inits$quantities$growth_matrix[1,1,,,1]
sum(t(wham_model$rep$catch_phi_mat[,,5]) - ss_inits$quantities$growth_matrix[4,1,,,5])
sum(t(wham_model$rep$fix_phi_mat[,,5]) - ss_inits$quantities$growth_matrix[4,1,,,5])

# Waa
plot(y = wham_model$rep$pred_waa[1,1,], x = t(ss_inits$quantities$weight_hat[3,1,,1])); abline(0,1) # Survey/Jan1 weight
plot(y = wham_model$rep$pred_waa[2,1,], x = t(ss_inits$quantities$weight_hat[4,1,,1])); abline(0,1) # Fishery weight

# Catch
plot(ss_inits$quantities$catch_hat[1:nyrs], wham_model$rep$pred_catch[,1]); abline(0,1)

# Index
plot(ss_inits$quantities$index_hat, wham_model$rep$pred_indices[,1]); abline(0,1)

# - IAAL
yr = 10
plot(x = t(wham_model$rep$pred_IAAL[yr,1,,]), y = ss_inits$quantities$pred_CAAL[1,1,,,yr]); abline(0, 1)

# - CAAL
yr = 10
plot(x = t(wham_model$rep$pred_CAAL[yr,1,,]), y = ss_inits$quantities$pred_CAAL[2,1,,,yr]); abline(0, 1)


# - COMP
par(mfrow = c(1, 2))
yr = 10
plot(x = wham_model$rep$pred_index_pal[,1,], y = ss_inits$quantities$comp_hat[1:nyrs,], xlab = "WHAM index comp", ylab = "Rceattle (fixed) index comp"); abline(0, 1) # Index

plot(x = wham_model$rep$pred_catch_pal[,1,], y = ss_inits$quantities$comp_hat[(nyrs+1):(nyrs*2),]); abline(0, 1) # Fishery


# - Likelihood
ss_inits$quantities$jnll_comp
sum(wham_model$rep$nll_agg_catch)
sum(wham_model$rep$nll_catch_lcomp)

sum(wham_model$rep$nll_agg_indices)
sum(wham_model$rep$nll_index_lcomp) # Comps use LL model 1
sum(wham_model$rep$nll_index_caal)

wham_model$rep$nll_NAA # Uses lognormal bias correction


# * Estimate Rceattle ----
ss_est <- Rceattle::fit_mod(data_list = simData,
                            inits = inits, # Initial parameters = 0
                            estimateMode = 0, # estimate
                            growthFun = growthFun_vb, # Von Bert, growth SD fixed to match WHAM
                            random_rec = TRUE, # No random recruitment
                            msmMode = 0, # Single species mode
                            initMode = 1, # Equilibrium (matches WHAM N1_model=1, init_dev off)
                            # comp_offset = 0 -> WHAM-style multinomial (no proportion offset)
                            fit_control = fit_control(phase = TRUE, verbose = 0, comp_offset = 0))


# Compare ----
wham_ests <- ss_est
wham_ests$quantities$ssb[1,1:nyrs] <- wham_model$rep$SSB/2

plot_ssb(list(ss_est, ss_inits, wham_ests), model_names = c("Rceattle", "Rceattle - fixed", "WHAM"))
#
#
#
# # * Make index by hand ----
# pred_IAAL <- matrix(0, nlengths, nages)
# index1 <- c()
# pred_IAAL2 <- matrix(0, nlengths, nages)
# index2 <- c()
# selLL <- wham_model$rep$selLL[[2]]
# out_phi_mat <- wham_model$rep$ssb_phi_mat
# NAA <- wham_model$rep$NAA
# pred_waa <- wham_model$rep$pred_waa[1,,]
# yr = 10
#
#
# selLL2 <- ss_inits$quantities$sel_at_length[1,1,,yr]
# out_phi_mat2 <- ss_inits$quantities$growth_matrix[3,1,,,yr]
# NAA2 <- ss_inits$quantities$N_at_age[1,1,,yr]
# pred_waa2 <- ss_inits$quantities$weight_hat[3,1,,yr]
#
# for(a in 1:nages) {
#   for(l in 1:nlengths) {
#     pred_IAAL[l,a] = selLL[yr,l]*out_phi_mat[l,a,yr]*NAA[yr,a];
#     pred_IAAL2[l,a] = selLL2[l]*out_phi_mat2[a,l]*NAA2[a];
#   }
#   index1[a] <- sum(pred_IAAL[,a]) * pred_waa[yr,a]
#   index2[a] <- sum(pred_IAAL2[,a]) # * pred_waa2[a]
# }
# t(pred_IAAL)
# t(pred_IAAL2)
#
# sum(wham_model$rep$pred_IAA[yr,1,] * wham_model$rep$pred_waa[1,yr,])
# wham_model$rep$pred_indices[yr,1]
# ss_inits$quantities$index_hat[yr]
# sum(index1)
# sum(index2)
#
# # - Catch
# pred_CAAL <- matrix(0, nlengths, nages)
# pred_CAAL2 <- matrix(0, nlengths, nages)
# selLL <- wham_model$rep$selLL[[1]]
# selAA <- wham_model$rep$selAA[[1]]
# out_phi_mat <- wham_model$rep$catch_phi_mat
# NAA <- wham_model$rep$NAA
# ZAA <- wham_model$rep$ZAA
# MAA <- wham_model$rep$MAA
# Frate <- wham_model$rep$F[,1]
# yr = 10
# for(l in 1:nlengths) {
#   for(a in 1:nages) {
#     pred_CAAL[l,a] = selAA[yr,a] * selLL[yr,l] * Frate[yr] * out_phi_mat[l,a,yr] * NAA[yr,a] / ZAA[yr,a] * (1-exp(-ZAA[yr,a]))
#   }
# }
# t(pred_CAAL)
# t(wham_model$rep$pred_CAAL[yr,1,,]) - ss_inits$quantities$pred_CAAL[2,1,,,yr]
#
