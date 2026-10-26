# Integration test (fits a CEATTLE model): skipped on CRAN to keep R CMD check
# fast; the coverage workflow runs it in full (NOT_CRAN=true). See README.md.
testthat::skip_on_cran()

testthat::test_that("Test age-based non-parametric selectivity bin-first selected", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("GOA2018SS") # Single-species data. ?BS2017SS for more information on the data

  # Adjust data
  n_sel_bins <- 8
  rows_use <- which(GOA2018SS$fleet_control$Selectivity != 0 & GOA2018SS$fleet_control$Fleet_type != 0)
  GOA2018SS$fleet_control$Fleet_type[-rows_use] <- 0
  GOA2018SS$fleet_control$Selectivity <- "Logistic" # Age-based logistic
  GOA2018SS$fleet_control$N_sel_bins <- n_sel_bins
  GOA2018SS$fleet_control$Sel_curve_pen1 <- 5
  GOA2018SS$fleet_control$Sel_curve_pen2 <- 10
  GOA2018SS$fleet_control$Sel_norm_bin <- NA
  GOA2018SS$fleet_control$Bin_first_selected <- 3
  GOA2018SS$fleet_control$Time_varying_sel <- 0
  GOA2018SS$fleet_control$Time_varying_sel_sd <- 1


  # Set params
  mod0 <- suppressMessages( fit_mod(data_list = GOA2018SS, inits = NULL, estimateMode = 3, random_rec = FALSE, msmMode = 0, fit_control = fit_control(verbose = 0)) )
  inits <- mod0$estimated_params
  log_selcoffs <- rnorm(n_sel_bins)
  log_selcoffs2 <- rnorm(n_sel_bins)
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

  # Check selectivity is 0 for bins 1-2
  output <- c(ss_run$quantities$sel_at_age[,1,1:2,])
  testthat::expect_equal(output, rep(0, length(output)))

  # Map for sel coff is NA?
  outputmap <- c(ss_run$map$mapList$sel_coff[,,1:2])
  testthat::expect_equal(outputmap, as.numeric(rep(NA, length(outputmap))))
})


testthat::test_that("Time-varying double logistic selectivity bin first selected", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Data
  data("GOA2018SS")
  rows_use <- which(GOA2018SS$fleet_control$Selectivity != 0 & GOA2018SS$fleet_control$Fleet_type != 0)
  GOA2018SS$fleet_control$Fleet_type[-rows_use] <-0
  GOA2018SS$fleet_control$Selectivity <- "DoubleLogistic" # Age-based double logistic
  GOA2018SS$fleet_control$Selectivity_index <- 1:nrow(GOA2018SS$fleet_control)
  GOA2018SS$fleet_control$Time_varying_sel <- 1
  GOA2018SS$fleet_control$Time_varying_sel_sd <- 1 # Note that inf does 4*sd prior, should adapt to scale-invariant sd
  #FIXME: move to scale invariant setup
  GOA2018SS$fleet_control$Bin_first_selected <- 4
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

  # Check selectivity is 0 for bins 1-2
  output <- c(ss_run$quantities$sel_at_age[,1,1:3,])
  testthat::expect_equal(output, rep(0, length(output)))
})
