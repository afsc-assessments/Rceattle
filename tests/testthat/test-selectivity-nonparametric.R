
# Integration test (fits a CEATTLE model): skipped on CRAN to keep R CMD check
# fast; the coverage workflow runs it in full (NOT_CRAN=true). See README.md.
testthat::skip_on_cran()

testthat::test_that("Test age-based non-parametric selectivity not normalized", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("GOA2018SS") # Single-species data. ?BS2017SS for more information on the data

  # Adjust data
  rows_use <- which(GOA2018SS$fleet_control$Selectivity != 0)
  GOA2018SS$fleet_control$Selectivity[rows_use] <- "NonParametric"
  GOA2018SS$fleet_control$N_sel_bins[rows_use] <- 8
  GOA2018SS$fleet_control$Age_first_selected[rows_use] <- 1
  GOA2018SS$fleet_control$Sel_curve_pen1[rows_use] <- 5
  GOA2018SS$fleet_control$Sel_curve_pen2[rows_use] <- 10
  GOA2018SS$fleet_control$Sel_norm_bin[rows_use] <- NA
  GOA2018SS$fleet_control$Time_varying_sel[rows_use] <- 0
  GOA2018SS$fleet_control$Time_varying_sel_sd[rows_use] <- 1

  mod0 <- suppressMessages( fit_mod(data_list = GOA2018SS, inits = NULL, estimateMode = 3, random_rec = FALSE, msmMode = 0, fit_control = fit_control(verbose = 0)) )
  inits <- mod0$estimated_params


  # ADMB code
  # dvar_vector sel_coffs_tmp(1,n_sel_bins_fsh(k));
  # for (i=styr;i<=endyr;i++)
  # {
  #   if (i==yrs_sel_ch_fsh(k,isel_ch_tmp))
  #   {
  #     sel_coffs_tmp.initialize();
  #     sel_coffs_tmp = log_selcoffs_fsh(k,isel_ch_tmp);
  #     avgsel_fsh(k,isel_ch_tmp)              = log(mean(mfexp(sel_coffs_tmp)));
  #   }
  #   log_sel_fsh(k,i)(1,n_sel_bins_fsh(k))        = sel_coffs_tmp;
  #   log_sel_fsh(k,i)(n_sel_bins_fsh(k),nages)    = log_sel_fsh(k,i,n_sel_bins_fsh(k));
  #   log_sel_fsh(k,i)                                  -= log(mean(mfexp(log_sel_fsh(k,i) )));
  # }

  # Specify selectivity
  n_sel_bins <- 8
  log_selcoffs <- rnorm(n_sel_bins)
  log_selcoffs2 <- rnorm(n_sel_bins)
  avgsel_fsh <- log(mean(exp(log_selcoffs)))

  # - Pk
  log_sel_fsh <- c(log_selcoffs, rep(log_selcoffs[n_sel_bins], GOA2018SS$nages[1]-n_sel_bins))
  log_sel_fsh <- log_sel_fsh - log(mean(exp(log_sel_fsh)))

  # - ATF
  log_sel_fsh2 <- c(log_selcoffs, rep(log_selcoffs[n_sel_bins], GOA2018SS$nages[2]-n_sel_bins))
  log_sel_fsh2 <- log_sel_fsh2 - log(mean(exp(log_sel_fsh2)))

  log_sel_fsh2m <- c(log_selcoffs2, rep(log_selcoffs2[n_sel_bins], GOA2018SS$nages[2]-n_sel_bins))
  log_sel_fsh2m <- log_sel_fsh2m - log(mean(exp(log_sel_fsh2m)))

  # - Cod
  log_sel_fsh3 <- c(log_selcoffs, rep(log_selcoffs[n_sel_bins], GOA2018SS$nages[3]-n_sel_bins))
  log_sel_fsh3 <- log_sel_fsh3 - log(mean(exp(log_sel_fsh3)))

  # Set params
  inits$sel_coff[,1,1:8] <- rep(log_selcoffs, each = dim(inits$sel_coff)[1])
  inits$sel_coff[9:11,2,1:8] <- rep(log_selcoffs2, each = 3) # Males

  # Run
  ss_run <- Rceattle::fit_mod(data_list = GOA2018SS,
                              inits = inits, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              fit_control = fit_control(
                                verbose = 0))

  # Check selectivity
  # - Pollock
  apply(ss_run$quantities$sel_at_age[1,1,1:10,], 2, function(x) testthat::expect_equal(as.numeric(x), exp(log_sel_fsh), tolerance = 0.0001))

  # - ATF
  apply(ss_run$quantities$sel_at_age[9,1,1:21,], 2, function(x) testthat::expect_equal(as.numeric(x), exp(log_sel_fsh2), tolerance = 0.0001))
  apply(ss_run$quantities$sel_at_age[9,2,1:21,], 2, function(x) testthat::expect_equal(as.numeric(x), exp(log_sel_fsh2m), tolerance = 0.0001))

  # - Cod
  apply(ss_run$quantities$sel_at_age[12,1,1:12,], 2, function(x) testthat::expect_equal(as.numeric(x), exp(log_sel_fsh3), tolerance = 0.0001))

  # Check penalties
  # - Curvature penalty
  pen1 <- sum(10*(diff(diff(log_sel_fsh))^2))

  # - Descending penalty
  difftmp = -diff(log_sel_fsh) # Decreasing will be positive
  difftmp = (abs(difftmp) + difftmp)/2
  pen2 <- sum(difftmp^2) * 5

  # - Avgsel
  pen3 <- 2 * log(mean(exp(log_selcoffs)))^2

  testthat::expect_equal(pen1+pen2+pen3, as.numeric(ss_run$quantities$jnll_comp[5,8]), tolerance = 0.0001)
})

testthat::test_that("2DAR1 selectivity map and likelihood", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Data
  data("GOA2018SS")
  nyrs <- length(GOA2018SS$styr:GOA2018SS$endyr)
  rows_use <- which(GOA2018SS$fleet_control$Selectivity != 0 & GOA2018SS$fleet_control$Fleet_type != 0)
  GOA2018SS$fleet_control$Selectivity[rows_use] <- "2DAR1" # Age-based 2DAR1
  GOA2018SS$fleet_control$Selectivity_index <- 1:nrow(GOA2018SS$fleet_control)
  GOA2018SS$fleet_control$Time_varying_sel_sd[rows_use] <- 1
  GOA2018SS$fleet_control$Age_first_selected[rows_use] <- 1
  GOA2018SS$fleet_control$N_sel_bins[rows_use] <- 8
  GOA2018SS$fleet_control$Sel_curve_pen1[rows_use] <- 0
  GOA2018SS$fleet_control$Sel_curve_pen2[rows_use] <- 0
  GOA2018SS$fleet_control$Sel_norm_bin[rows_use] <- NA # Do not normalize
  GOA2018SS$catch_data$Catch <- 1e6 # If catch is zero, sel devs are turned off

  # Run
  ss_run <- Rceattle::fit_mod(data_list = GOA2018SS,
                              inits = NULL, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              random_sel = TRUE, # Turn on laplace for sel devs
                              msmMode = 0, # Single species mode
                              fit_control = fit_control(
                                phase = TRUE,
                                verbose = 0))

  # Hyper parameters (on except for not estimated fleet)
  na_out <- rep(NA, 16)
  na_out[rows_use] = rows_use
  testthat::expect_equal(as.numeric(ss_run$map$mapList$sel_dev_log_sd), na_out)
  testthat::expect_equal(as.numeric(ss_run$map$mapList$sel_curve_pen[,1]), na_out)
  testthat::expect_equal(as.numeric(ss_run$map$mapList$sel_curve_pen[,2]), na_out + 16)
  testthat::expect_equal(as.numeric(ss_run$map$mapList$sel_curve_pen[,3]), as.numeric(rep(NA, 16)))


  na_out <- rep(-Inf, 16)
  na_out[rows_use] = 0
  testthat::expect_equal(as.numeric(ss_run$estimated_params$sel_dev_log_sd), na_out)
  expect_all_true(as.numeric(ss_run$estimated_params$sel_curve_pen[,1]) == 0)

  # Fixed effects
  testthat::expect_equal(as.numeric(ss_run$estimated_params$sel_coff[,,1:8]), rep(0, 2 * 8 * 16))
  # - Females (sex combined)
  expect_all_true(!is.na(as.numeric(ss_run$map$mapList$sel_coff[rows_use,1,1:8])))
  expect_all_true(is.na(as.numeric(ss_run$map$mapList$sel_coff[7,1,1:8])))
  testthat::expect_equal(length(unique(ss_run$map$mapList$sel_coff[-7,1,1:8])), 8 * 13)

  # - Males
  flt2sex <- GOA2018SS$nsex[GOA2018SS$fleet_control$Species] == 2
  expect_all_true(!is.na(as.numeric(ss_run$map$mapList$sel_coff[flt2sex,2,1:8])))
  expect_all_true(is.na(as.numeric(ss_run$map$mapList$sel_coff[!flt2sex,2,1:8])))
  testthat::expect_equal(length(unique(ss_run$map$mapList$sel_coff[flt2sex,2,1:8])), 8 * sum(flt2sex))

  # Random effects
  testthat::expect_equal(as.numeric(ss_run$estimated_params$sel_coff[,1,1:8]), rep(0, 8 * 16))

  # - Females (sex combined)
  expect_all_true(!is.na(as.numeric(ss_run$map$mapList$sel_coff_dev[rows_use,1,1:8,1])))
  expect_all_true(is.na(as.numeric(ss_run$map$mapList$sel_coff_dev[-rows_use,1,1:8,])))
  testthat::expect_equal(length(unique(ss_run$map$mapList$sel_coff_dev[-7,1,1:8,])), 8 * 13 * nyrs)


  # - Males
  expect_all_true(!is.na(as.numeric(ss_run$map$mapList$sel_coff_dev[flt2sex,2,1:8,])))
  expect_all_true(is.na(as.numeric(ss_run$map$mapList$sel_coff_dev[!flt2sex,2,1:8,])))
  testthat::expect_equal(length(unique(ss_run$map$mapList$sel_coff_dev[flt2sex,2,1:8,])), 8 * sum(flt2sex) * nyrs)

  # TMB object
  testthat::expect_equal(length(unique(ss_run$obj$env$random)),  8 * 15 * nyrs)

  # JNLL
  expect_all_true(ss_run$quantities$jnll_comp[6,rows_use] != 0)
  expect_all_true(ss_run$quantities$jnll_comp[6,-rows_use] == 0)
  expect_all_true(!is.na(ss_run$quantities$jnll_comp[6,]))
})

testthat::test_that("3DAR1 selectivity map and likelihood", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Data
  data("GOA2018SS")
  nyrs <- length(GOA2018SS$styr:GOA2018SS$endyr)
  rows_use <- which(GOA2018SS$fleet_control$Selectivity != 0 & GOA2018SS$fleet_control$Fleet_type != 0)
  GOA2018SS$fleet_control$Selectivity <-0
  GOA2018SS$fleet_control$Selectivity[rows_use] <- "3DAR1" # Age-based 2DAR1
  GOA2018SS$fleet_control$Selectivity_index <- 1:nrow(GOA2018SS$fleet_control)
  GOA2018SS$fleet_control$Time_varying_sel_sd[rows_use] <- 1
  GOA2018SS$fleet_control$Age_first_selected[rows_use] <- 1
  GOA2018SS$fleet_control$N_sel_bins[rows_use] <- 8
  GOA2018SS$fleet_control$Sel_curve_pen1[rows_use] <- 0
  GOA2018SS$fleet_control$Sel_curve_pen2[rows_use] <- 0
  GOA2018SS$fleet_control$Sel_norm_bin[rows_use] <- NA # Do not normalize
  GOA2018SS$catch_data$Catch <- 1e6 # If catch is zero, sel devs are turned off

  # Run
  ss_run <- Rceattle::fit_mod(data_list = GOA2018SS,
                              inits = NULL, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              random_sel = TRUE, # Turn on laplace for sel devs
                              msmMode = 0, # Single species mode
                              fit_control = fit_control(
                                phase = TRUE,
                                verbose = 0))

  # Hyper parameters (on except for not estimated fleet)
  na_out <- rep(NA, 16)
  na_out[rows_use] = rows_use
  testthat::expect_equal(as.numeric(ss_run$map$mapList$sel_dev_log_sd), na_out)
  testthat::expect_equal(as.numeric(ss_run$map$mapList$sel_curve_pen[,1]), na_out)
  testthat::expect_equal(as.numeric(ss_run$map$mapList$sel_curve_pen[,2]), na_out + 16)
  testthat::expect_equal(as.numeric(ss_run$map$mapList$sel_curve_pen[,3]), na_out + 16*2)


  na_out <- rep(-Inf, 16)
  na_out[rows_use] = 0
  testthat::expect_equal(as.numeric(ss_run$estimated_params$sel_dev_log_sd), na_out)
  expect_all_true(as.numeric(ss_run$estimated_params$sel_curve_pen[,1]) == 0)

  # Fixed effects
  testthat::expect_equal(as.numeric(ss_run$estimated_params$sel_coff[,,1:8]), rep(0, 2 * 8 * 16))
  # - Females (sex combined)
  expect_all_true(!is.na(as.numeric(ss_run$map$mapList$sel_coff[rows_use,1,1:8])))
  expect_all_true(is.na(as.numeric(ss_run$map$mapList$sel_coff[7,1,1:8])))
  testthat::expect_equal(length(unique(ss_run$map$mapList$sel_coff[-7,1,1:8])), 8 * 13)

  # - Males
  flt2sex <- GOA2018SS$nsex[GOA2018SS$fleet_control$Species] == 2
  expect_all_true(!is.na(as.numeric(ss_run$map$mapList$sel_coff[flt2sex,2,1:8])))
  expect_all_true(is.na(as.numeric(ss_run$map$mapList$sel_coff[!flt2sex,2,1:8])))
  testthat::expect_equal(length(unique(ss_run$map$mapList$sel_coff[flt2sex,2,1:8])), 8 * sum(flt2sex))

  # Random effects
  testthat::expect_equal(as.numeric(ss_run$estimated_params$sel_coff[,1,1:8]), rep(0, 8 * 16))

  # - Females (sex combined)
  expect_all_true(!is.na(as.numeric(ss_run$map$mapList$sel_coff_dev[rows_use,1,1:8,1])))
  expect_all_true(is.na(as.numeric(ss_run$map$mapList$sel_coff_dev[-rows_use,1,1:8,])))
  testthat::expect_equal(length(unique(ss_run$map$mapList$sel_coff_dev[-7,1,1:8,])), 8 * 13 * nyrs)


  # - Males
  expect_all_true(!is.na(as.numeric(ss_run$map$mapList$sel_coff_dev[flt2sex,2,1:8,])))
  expect_all_true(is.na(as.numeric(ss_run$map$mapList$sel_coff_dev[!flt2sex,2,1:8,])))
  testthat::expect_equal(length(unique(ss_run$map$mapList$sel_coff_dev[flt2sex,2,1:8,])), 8 * sum(flt2sex) * nyrs)

  # TMB object
  testthat::expect_equal(length(unique(ss_run$obj$env$random)),  8 * 15 * nyrs)

  # JNLL
  expect_all_true(ss_run$quantities$jnll_comp[6,rows_use] != 0)
  expect_all_true(ss_run$quantities$jnll_comp[6,-rows_use] == 0)
  expect_all_true(!is.na(ss_run$quantities$jnll_comp[6,]))
})

testthat::test_that("2DAR1/3DAR1 reject saturating Sel_curve_pen rho values", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Data
  data("GOA2018SS")
  nyrs <- length(GOA2018SS$styr:GOA2018SS$endyr)
  rows_use <- which(GOA2018SS$fleet_control$Selectivity != 0 & GOA2018SS$fleet_control$Fleet_type != 0)
  GOA2018SS$fleet_control$Selectivity <-0
  GOA2018SS$fleet_control$Selectivity[rows_use] <- "2DAR1" # Age-based 2DAR1
  GOA2018SS$fleet_control$Selectivity_index <- 1:nrow(GOA2018SS$fleet_control)
  GOA2018SS$fleet_control$Time_varying_sel_sd[rows_use] <- 1
  GOA2018SS$fleet_control$Age_first_selected[rows_use] <- 1
  GOA2018SS$fleet_control$N_sel_bins[rows_use] <- 8
  GOA2018SS$fleet_control$Sel_curve_pen1[rows_use] <- 0
  GOA2018SS$fleet_control$Sel_curve_pen2[rows_use] <- 0
  GOA2018SS$fleet_control$Sel_norm_bin[rows_use] <- NA # Do not normalize
  GOA2018SS$catch_data$Catch <- 1e6 # If catch is zero, sel devs are turned off

  # Sel_curve_pen1 above the saturation threshold should error for 2DAR1
  bad_pen1 <- GOA2018SS
  bad_pen1$fleet_control$Sel_curve_pen1[rows_use] <- 20
  testthat::expect_error(
    suppressMessages(Rceattle::fit_mod(data_list = bad_pen1, inits = NULL, file = NULL,
                                       estimateMode = 3, random_rec = FALSE, random_sel = TRUE,
                                       msmMode = 0, fit_control = fit_control(phase = TRUE, verbose = 0))),
    regexp = "Sel_curve_pen1"
  )

  # Sel_curve_pen2 below the negative saturation threshold should error for 3DAR1
  bad_pen2 <- GOA2018SS
  bad_pen2$fleet_control$Selectivity[rows_use] <- "3DAR1"
  bad_pen2$fleet_control$Sel_curve_pen2[rows_use] <- -12.5
  testthat::expect_error(
    suppressMessages(Rceattle::fit_mod(data_list = bad_pen2, inits = NULL, file = NULL,
                                       estimateMode = 3, random_rec = FALSE, random_sel = TRUE,
                                       msmMode = 0, fit_control = fit_control(phase = TRUE, verbose = 0))),
    regexp = "Sel_curve_pen2"
  )

  # NonParametric selectivity uses Sel_curve_pen as a real penalty weight; large values are valid
  ok_np <- GOA2018SS
  ok_np$fleet_control$Selectivity[rows_use] <- "NonParametric"
  ok_np$fleet_control$Sel_curve_pen1[rows_use] <- 50
  # GOA2018SS's fishery is RandomWalkAscending, which the non-parametric forms do
  # not take. Set a mode they do, so this tests the penalty weight and not the
  # mode: until 5.25.0 the mode was silently read as a legacy AMAK weight, which
  # also overwrote the Sel_curve_pen1 this test sets.
  ok_np$fleet_control$Time_varying_sel[rows_use] <- "Off"
  testthat::expect_no_error(
    suppressMessages(Rceattle::fit_mod(data_list = ok_np, inits = NULL, file = NULL,
                                       estimateMode = 3, random_rec = FALSE,
                                       msmMode = 0, fit_control = fit_control(verbose = 0)))
  )
})
