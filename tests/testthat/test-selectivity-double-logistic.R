# Integration test (fits a CEATTLE model): skipped on CRAN to keep R CMD check
# fast; the coverage workflow runs it in full (NOT_CRAN=true). See README.md.
testthat::skip_on_cran()

testthat::test_that("Sex-specific age-based double logistic selectivity not normalized", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Data
  data("GOA2018SS")
  rows_use <- which(GOA2018SS$fleet_control$Selectivity != 0 & GOA2018SS$fleet_control$Fleet_type != 0)
  GOA2018SS$fleet_control$Selectivity <-0
  GOA2018SS$fleet_control$Selectivity[rows_use] <- "DoubleLogistic"
  GOA2018SS$fleet_control$Time_varying_sel <- 0
  GOA2018SS$fleet_control$Selectivity_index <- 1:nrow(GOA2018SS$fleet_control)
  GOA2018SS$fleet_control$Bin_first_selected <- 1
  GOA2018SS$fleet_control$Sel_norm_bin <- NA

  # Specify double logistic selectivity
  inf1 = 7; alpha = 0.2
  inf2 = 15
  ages <- 1:21
  sel <- 1/(1+exp(-alpha*(ages-inf1))) * (1/(1+exp(alpha*(ages-inf2))))
  sel2 <- 1/(1+exp(-alpha*(ages-inf1)))*(1/(1+exp(alpha*(ages-16))))

  # Set params to double logistic
  mod0 <- suppressMessages( fit_mod(data_list = GOA2018SS, inits = NULL, estimateMode = 3, random_rec = FALSE, msmMode = 0, fit_control = fit_control(verbose = 0)) )
  inits <- mod0$estimated_params
  inits$log_sel_slp[,,] <- log(alpha)
  inits$sel_inf[1,,1] <- inf1     # Females
  inits$sel_inf[2,,1] <- inf2     # Females
  inits$sel_inf[1,9:11,2] <- inf1 # Males
  inits$sel_inf[2,9:11,2] <- 16 # Males

  # Run
  ss_run <- suppressMessages(
    Rceattle::fit_mod(data_list = GOA2018SS,
                      inits = inits, # Initial parameters = 0
                      file = NULL, # Don't save
                      estimateMode = 3, # Don't estimate
                      random_rec = FALSE, # No random recruitment
                      msmMode = 0, # Single species mode
                      fit_control = fit_control(
                        verbose = 1))
  )

  # Map
  # - Slopes
  testthat::expect_equal(sum(!is.na(c(ss_run$map$mapList$log_sel_slp[,,]))), 15 * 2)
  testthat::expect_equal(length(unique(c(ss_run$map$mapList$log_sel_slp[,,]))), 15 * 2 + 1) # +1 for NAs

  # - Asympts
  testthat::expect_equal(sum(!is.na(c(ss_run$map$mapList$sel_inf[,,]))), 15 * 2)
  testthat::expect_equal(length(unique(c(ss_run$map$mapList$sel_inf[,,]))), 15 * 2 + 1) # +1 for NAs

  # - Devs
  testthat::expect_equal(c(ss_run$map$mapList$log_sel_slp_dev), as.numeric(rep(NA, 2688)))
  testthat::expect_equal(c(ss_run$map$mapList$sel_inf_dev), as.numeric(rep(NA, 2688)))
  testthat::expect_equal(as.numeric(ss_run$map$mapList$sel_dev_log_sd), as.numeric(rep(NA, 16))) # Dev sigma


  # Check selectivity
  # - Pollock
  apply(ss_run$quantities$sel_at_age[1,1,1:10,], 2, function(x) testthat::expect_equal(as.numeric(x), sel[1:10], tolerance = 0.0001))

  # - ATF
  apply(ss_run$quantities$sel_at_age[9,1,1:21,], 2, function(x) testthat::expect_equal(as.numeric(x), sel[1:21], tolerance = 0.0001))
  apply(ss_run$quantities$sel_at_age[9,2,1:21,], 2, function(x) testthat::expect_equal(as.numeric(x), sel2[1:21], tolerance = 0.0001))
})


testthat::test_that("Sex-specific age-based time-varying double logistic selectivity not normalized", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Data
  data("GOA2018SS")
  rows_use <- which(GOA2018SS$fleet_control$Selectivity != 0 & GOA2018SS$fleet_control$Fleet_type != 0)
  GOA2018SS$fleet_control$Selectivity <-0
  GOA2018SS$fleet_control$Selectivity[rows_use] <- "DoubleLogistic"
  GOA2018SS$fleet_control$Selectivity_index <- 1:nrow(GOA2018SS$fleet_control)
  GOA2018SS$fleet_control$Time_varying_sel <- 1
  GOA2018SS$fleet_control$Time_varying_sel_sd <- 1
  GOA2018SS$fleet_control$Bin_first_selected <- 1
  GOA2018SS$fleet_control$Sel_norm_bin <- NA # Do not normalize
  GOA2018SS$catch_data$Catch <- 1e6 # If catch is zero, sel devs are turned off

  # Specify double logistic selectivity
  nyrs <- length(GOA2018SS$styr:GOA2018SS$endyr)
  inf = 7
  inf2 = 15
  inf_dev <- rnorm(nyrs)
  inf_dev_desc <- rnorm(nyrs)
  log_slp_dev <- rnorm(nyrs)
  log_slp_dev_desc <- rnorm(nyrs)

  alpha = 0.5
  ages <- 1:21
  sel <- apply(cbind(log_slp_dev, inf_dev, log_slp_dev_desc, inf_dev_desc), 1,
               function(x) 1/(1+exp(-(alpha)*exp(x[1]) * (ages - inf - x[2]))) * (1/(1+exp((alpha)*exp(x[3]) * (ages - inf2 - x[4])))))
  sel2 <- apply(cbind(log_slp_dev, inf_dev, log_slp_dev_desc, inf_dev_desc), 1,
                function(x) 1/(1+exp(-(alpha+1)*exp(x[1]) * (ages - inf - x[2]))) * (1/(1+exp((alpha+1)*exp(x[3]) * (ages - inf2 - x[4])))))


  # Set params to double logistic
  mod0 <- suppressMessages( fit_mod(data_list = GOA2018SS, inits = NULL, estimateMode = 3, random_rec = FALSE, msmMode = 0, fit_control = fit_control(verbose = 0)) )
  inits <- mod0$estimated_params
  inits$log_sel_slp[,,1] <- log(alpha) # Females
  inits$log_sel_slp[,,2] <- log(alpha+1) # Males
  inits$sel_inf[1,,] <- inf
  inits$sel_inf[2,,] <- inf2
  for(i in 1:dim(inits$log_sel_slp_dev[1,,,])[1]){
    for(j in 1:dim(inits$log_sel_slp_dev[1,,,])[2]){
      inits$log_sel_slp_dev[1,i,j,] <- log_slp_dev
      inits$sel_inf_dev[1,i,j,] <- inf_dev
      inits$log_sel_slp_dev[2,i,j,] <- log_slp_dev_desc
      inits$sel_inf_dev[2,i,j,] <- inf_dev_desc
    }
  }

  # Run
  ss_run <- suppressMessages(
    Rceattle::fit_mod(data_list = GOA2018SS,
                      inits = inits, # Initial parameters = 0
                      file = NULL, # Don't save
                      estimateMode = 3, # Don't estimate
                      random_rec = FALSE, # No random recruitment
                      msmMode = 0, # Single species mode
                      fit_control = fit_control(
                        verbose = 0))
  )

  # Map
  # - Slopes
  testthat::expect_equal(sum(!is.na(c(ss_run$map$mapList$log_sel_slp[,,]))), 15 * 2)
  testthat::expect_equal(length(unique(c(ss_run$map$mapList$log_sel_slp[,,]))), 15 * 2 + 1) # +1 for NAs

  # - Asympts
  testthat::expect_equal(sum(!is.na(c(ss_run$map$mapList$sel_inf[,,]))), 15 * 2)
  testthat::expect_equal(length(unique(c(ss_run$map$mapList$sel_inf[,,]))), 15 * 2 + 1) # +1 for NAs

  # - Devs
  # -- double
  testthat::expect_equal(sum(!is.na(c(ss_run$map$mapList$log_sel_slp_dev[,,,]))), 15 * 2 * nyrs) # slope
  testthat::expect_equal(sum(!is.na(c(ss_run$map$mapList$sel_inf_dev[,,,]))), 15 * 2 * nyrs) # asymptote

  testthat::expect_equal(length(unique(c(ss_run$map$mapList$log_sel_slp_dev[,,,]))), 15 * 2 * nyrs + 1) # slope
  testthat::expect_equal(length(unique(c(ss_run$map$mapList$sel_inf_dev[,,,]))), 15 * 2 * nyrs + 1) # asymptote

  testthat::expect_equal(as.numeric(ss_run$map$mapList$sel_dev_log_sd), as.numeric(rep(NA, 16))) # Dev sigma turned off


  # Check selectivity
  # - Pollock
  testthat::expect_equal(as.numeric(ss_run$quantities$sel_at_age[1,1,1:10,1:42]), as.numeric(sel[1:10,]), tolerance = 0.0001)

  # - ATF
  testthat::expect_equal(as.numeric(ss_run$quantities$sel_at_age[9,1,1:21,1:42]), as.numeric(sel[1:21,]), tolerance = 0.0001) # Females
  testthat::expect_equal(as.numeric(ss_run$quantities$sel_at_age[9,2,1:21,1:42]), as.numeric(sel2[1:21,]), tolerance = 0.0001) # Males
})


testthat::test_that("Sex-specific age-based time-varying double logistic selectivity not normalized (devs as random effects)", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Data
  data("GOA2018SS")
  nyrs <- length(GOA2018SS$styr:GOA2018SS$endyr)
  rows_use <- which(GOA2018SS$fleet_control$Selectivity != 0 & GOA2018SS$fleet_control$Fleet_type != 0)
  GOA2018SS$fleet_control$Selectivity <-0
  GOA2018SS$fleet_control$Selectivity[rows_use] <- "DoubleLogistic"
  GOA2018SS$fleet_control$Selectivity_index <- 1:nrow(GOA2018SS$fleet_control)
  GOA2018SS$fleet_control$Time_varying_sel <- 1
  GOA2018SS$fleet_control$Time_varying_sel_sd <- 1
  GOA2018SS$fleet_control$Bin_first_selected <- 1
  GOA2018SS$fleet_control$Sel_norm_bin <- NA # Do not normalize
  GOA2018SS$catch_data$Catch <- 1e6 # If catch is zero, sel devs are turned off

  # Run
  ss_run <- suppressMessages(
    Rceattle::fit_mod(data_list = GOA2018SS,
                      file = NULL, # Don't save
                      estimateMode = 3, # Don't estimate
                      random_rec = FALSE, # No random recruitment
                      random_sel = TRUE, # Turn on laplace for sel devs
                      msmMode = 0, # Single species mode
                      fit_control = fit_control(
                        verbose = 0))
  )

  # Sigmas
  testthat::expect_equal(as.numeric(ss_run$map$mapList$sel_dev_log_sd), c(1:3, rep(NA,4), 8:16)) # Dev sigma turned on except for not estimated fleet

  # TMB object
  testthat::expect_equal(length(unique(ss_run$obj$env$random)),  2 * 15 * nyrs * 2)

})


testthat::test_that("Time-varying double logistic selectivity likelihood", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Data
  data("GOA2018SS")
  rows_use <- which(GOA2018SS$fleet_control$Selectivity != 0 & GOA2018SS$fleet_control$Fleet_type != 0)
  GOA2018SS$fleet_control$Selectivity <-0
  GOA2018SS$fleet_control$Selectivity[rows_use] <- "DoubleLogistic"
  GOA2018SS$fleet_control$Selectivity_index <- 1:nrow(GOA2018SS$fleet_control)
  GOA2018SS$fleet_control$Time_varying_sel <- 1
  GOA2018SS$fleet_control$Time_varying_sel_sd <- 1 # Note that inf does 4*sd prior, should adapt to scale-invariant sd
  #FIXME: move to scale invariant setup
  GOA2018SS$fleet_control$Bin_first_selected <- 1
  GOA2018SS$fleet_control$Sel_norm_bin <- NA # Do not normalize
  GOA2018SS$catch_data$Catch <- 1e6 # If catch is zero, sel devs are turned off

  # Specify double logistic selectivity
  nyrs <- length(GOA2018SS$styr:GOA2018SS$endyr)
  inf = 7
  inf2 = 15
  inf_dev <- rnorm(nyrs)
  inf_dev_desc <- rnorm(nyrs)
  log_slp_dev <- rnorm(nyrs)
  log_slp_dev_desc <- rnorm(nyrs)

  alpha = 0.5
  ages <- 1:21
  sel <- apply(cbind(log_slp_dev, inf_dev, log_slp_dev_desc, inf_dev_desc), 1,
               function(x) 1/(1+exp(-(alpha)*exp(x[1]) * (ages - inf - x[2]))) * (1/(1+exp((alpha)*exp(x[3]) * (ages - inf2 - x[4])))))
  sel2 <- apply(cbind(log_slp_dev, inf_dev, log_slp_dev_desc, inf_dev_desc), 1,
                function(x) 1/(1+exp(-(alpha+1)*exp(x[1]) * (ages - inf - x[2]))) * (1/(1+exp((alpha+1)*exp(x[3]) * (ages - inf2 - x[4])))))


  # Set params to double logistic
  mod0 <- suppressMessages( fit_mod(data_list = GOA2018SS, inits = NULL, estimateMode = 3, random_rec = FALSE, random_sel = TRUE, msmMode = 0, fit_control = fit_control(verbose = 0)) )
  inits <- mod0$estimated_params
  inits$log_sel_slp[,,1] <- log(alpha) # Females
  inits$log_sel_slp[,,2] <- log(alpha+1) # Males
  inits$sel_inf[1,,] <- inf
  inits$sel_inf[2,,] <- inf2
  for(i in 1:dim(inits$log_sel_slp_dev[1,,,])[1]){
    for(j in 1:dim(inits$log_sel_slp_dev[1,,,])[2]){
      inits$log_sel_slp_dev[1,i,j,] <- log_slp_dev
      inits$sel_inf_dev[1,i,j,] <- inf_dev
      inits$log_sel_slp_dev[2,i,j,] <- log_slp_dev_desc
      inits$sel_inf_dev[2,i,j,] <- inf_dev_desc
    }
  }

  # Run
  ss_run <- suppressMessages(
    Rceattle::fit_mod(data_list = GOA2018SS,
                      inits = inits, # Initial parameters = 0
                      file = NULL, # Don't save
                      estimateMode = 3, # Don't estimate
                      random_rec = FALSE, # No random recruitment
                      random_sel = TRUE, # Turn on laplace for sel devs
                      msmMode = 0, # Single species mode
                      fit_control = fit_control(
                        verbose = 0))
  )

  # Nll
  rcnll <- sum(ss_run$quantities$jnll_comp[6,])
  single_dv_nll <- sum(-dnorm(inf_dev, 0, 1, log = TRUE) - dnorm(inf_dev_desc, 0, 1, log = TRUE) - dnorm(log_slp_dev, 0, 4, log = TRUE) - dnorm(log_slp_dev_desc, 0, 4, log = TRUE)) * 15

  testthat::expect_equal(rcnll, single_dv_nll)
})
