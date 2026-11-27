# Tests for various forms of logistic selectivity
# - Test correct value
# - Test correct mapping

# Integration test (fits a CEATTLE model): skipped on CRAN to keep R CMD check
# fast; the coverage workflow runs it in full (NOT_CRAN=true). See README.md.
testthat::skip_on_cran()

testthat::test_that("Sex-specific age-based logistic selectivity not normalized", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Data
  data("GOA2018SS")
  rows_use <- which(GOA2018SS$fleet_control$Selectivity != 0 & GOA2018SS$fleet_control$Fleet_type != 0)
  GOA2018SS$fleet_control$Selectivity <-0
  GOA2018SS$fleet_control$Selectivity[rows_use] <- "Logistic" # Age-based logistic
  GOA2018SS$fleet_control$Time_varying_sel <- 0
  GOA2018SS$fleet_control$Selectivity_index <- 1:nrow(GOA2018SS$fleet_control)
  GOA2018SS$fleet_control$Age_first_selected <- 1
  GOA2018SS$fleet_control$Sel_norm_bin <- NA

  # Specify logistic selectivity
  inf = 13; alpha = 0.2
  ages <- 1:21
  sel <- 1/(1+exp(-alpha*(ages-inf)))
  sel2 <- 1/(1+exp(-alpha*(ages-inf-1)))
  # curve(1/(1+exp(-alpha*(x-inf))), from = 0, to = 21)

  # Set params to logistic
  ss_run <- suppressMessages(
    Rceattle::fit_mod(data_list = GOA2018SS,
                      estimateMode = 3, # Don't estimate
                      msmMode = 0, # Single species mode
                      fit_control = fit_control(
                        verbose = 1))
  )
  inits <- ss_run$estimated_params
  inits$log_sel_slp[1,,] <- log(alpha)
  inits$sel_inf[1,,] <- inf     # Females
  inits$sel_inf[1,9:11,2] <- inf + 1 # Males

  # Run
  ss_run <- suppressMessages(
    Rceattle::fit_mod(data_list = GOA2018SS,
                      inits = inits, # Initial parameters = 0
                      estimateMode = 3, # Don't estimate
                      msmMode = 0, # Single species mode
                      fit_control = fit_control(
                        verbose = 1))
  )

  # Map
  # - Slope
  testthat::expect_equal(sum(!is.na(ss_run$map$mapList$log_sel_slp[1,,])), 15)
  testthat::expect_equal(length(unique(as.numeric(ss_run$map$mapList$log_sel_slp[1,,]))), 16)
  # - Desc slope
  testthat::expect_equal(c(ss_run$map$mapList$log_sel_slp[2,,]), as.numeric(rep(NA, 16*2)))

  # - Asympt
  testthat::expect_equal(sum(!is.na(ss_run$map$mapList$sel_inf[1,,])), 15)
  testthat::expect_equal(length(unique(as.numeric(ss_run$map$mapList$sel_inf[1,,]))), 16)
  # - Desc Asympt
  testthat::expect_equal(c(ss_run$map$mapList$sel_inf[2,,]), as.numeric(rep(NA, 16*2)))

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


testthat::test_that("Sex-specific age-based time-varying logistic selectivity not normalized", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Data
  data("GOA2018SS")
  rows_use <- which(GOA2018SS$fleet_control$Selectivity != 0 & GOA2018SS$fleet_control$Fleet_type != 0)
  GOA2018SS$fleet_control$Selectivity <-0
  GOA2018SS$fleet_control$Selectivity[rows_use] <- "Logistic" # Age-based logistic
  GOA2018SS$fleet_control$Selectivity_index <- 1:nrow(GOA2018SS$fleet_control)
  GOA2018SS$fleet_control$Time_varying_sel <- 1
  GOA2018SS$fleet_control$Time_varying_sel_sd <- 1
  GOA2018SS$fleet_control$Age_first_selected <- 1
  GOA2018SS$fleet_control$Sel_norm_bin <- NA # Do not normalize
  GOA2018SS$catch_data$Catch <- 1e6 # If catch is zero, sel devs are turned off

  # Specify logistic selectivity
  nyrs <- length(GOA2018SS$styr:GOA2018SS$endyr)
  inf = 10
  inf_dev <- rnorm(nyrs)
  log_slp_dev <- rnorm(nyrs)

  alpha = 0.5
  ages <- 1:21
  sel <- apply(cbind(log_slp_dev, inf_dev), 1, function(x) 1/(1+exp(-alpha*exp(x[1]) * (ages - inf - x[2]))))
  sel2 <- apply(cbind(log_slp_dev, inf_dev), 1, function(x) 1/(1+exp(-(alpha+1)*exp(x[1]) * (ages - inf - x[2]))))


  # Set params to logistic
  ss_run <- suppressMessages(
    Rceattle::fit_mod(data_list = GOA2018SS,
                      estimateMode = 3, # Don't estimate
                      msmMode = 0, # Single species mode
                      fit_control = fit_control(
                        verbose = 1))
  )
  inits <- ss_run$estimated_params
  inits$log_sel_slp[1,,1] <- log(alpha) # Females
  inits$log_sel_slp[1,,2] <- log(alpha+1) # Males
  inits$sel_inf[1,,] <- inf
  for(i in 1:dim(inits$log_sel_slp_dev[1,,,])[1]){
    for(j in 1:dim(inits$log_sel_slp_dev[1,,,])[2]){
      inits$log_sel_slp_dev[1,i,j,] <- log_slp_dev
      inits$sel_inf_dev[1,i,j,] <- inf_dev
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
  # - Slope
  testthat::expect_equal(sum(!is.na(ss_run$map$mapList$log_sel_slp[1,,])), 15)
  testthat::expect_equal(length(unique(as.numeric(ss_run$map$mapList$log_sel_slp[1,,]))), 16)
  # - Desc slope
  testthat::expect_equal(c(ss_run$map$mapList$log_sel_slp[2,,]), as.numeric(rep(NA, 16*2)))

  # - Asympt
  testthat::expect_equal(sum(!is.na(ss_run$map$mapList$sel_inf[1,,])), 15)
  testthat::expect_equal(length(unique(as.numeric(ss_run$map$mapList$sel_inf[1,,]))), 16)
  # - Desc Asympt
  testthat::expect_equal(c(ss_run$map$mapList$sel_inf[2,,]), as.numeric(rep(NA, 16*2)))

  # - Devs
  # -- Ascending
  testthat::expect_equal(sum(!is.na(c(ss_run$map$mapList$log_sel_slp_dev[1,,,]))), 15 * nyrs) # slope
  testthat::expect_equal(sum(!is.na(c(ss_run$map$mapList$sel_inf_dev[1,,,]))), 15 * nyrs) # asymptote

  testthat::expect_equal(length(unique(c(ss_run$map$mapList$log_sel_slp_dev[1,,,]))), 15 * nyrs + 1) # slope
  testthat::expect_equal(length(unique(c(ss_run$map$mapList$sel_inf_dev[1,,,]))), 15 * nyrs + 1) # asymptote

  testthat::expect_equal(c(ss_run$map$mapList$log_sel_slp_dev[2,,,]), as.numeric(rep(NA, 2688/2))) # Descending slope
  testthat::expect_equal(c(ss_run$map$mapList$sel_inf_dev[2,,,]), as.numeric(rep(NA, 2688/2))) # Descending asymptote

  testthat::expect_equal(as.numeric(ss_run$map$mapList$sel_dev_log_sd), as.numeric(rep(NA, 16))) # Dev sigma turned off


  # Check selectivity
  # - Pollock
  testthat::expect_equal(as.numeric(ss_run$quantities$sel_at_age[1,1,1:10,1:42]), as.numeric(sel[1:10,]), tolerance = 0.0001)

  # - ATF
  testthat::expect_equal(as.numeric(ss_run$quantities$sel_at_age[9,1,1:21,1:42]), as.numeric(sel[1:21,]), tolerance = 0.0001) # Females
  testthat::expect_equal(as.numeric(ss_run$quantities$sel_at_age[9,2,1:21,1:42]), as.numeric(sel2[1:21,]), tolerance = 0.0001) # Males
})


testthat::test_that("Sex-specific age-based time-varying logistic selectivity not normalized (devs as random effects)", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Data
  data("GOA2018SS")
  nyrs <- length(GOA2018SS$styr:GOA2018SS$endyr)
  rows_use <- which(GOA2018SS$fleet_control$Selectivity != 0 & GOA2018SS$fleet_control$Fleet_type != 0)
  GOA2018SS$fleet_control$Selectivity <-0
  GOA2018SS$fleet_control$Selectivity[rows_use] <- "Logistic" # Age-based logistic
  GOA2018SS$fleet_control$Selectivity_index <- 1:nrow(GOA2018SS$fleet_control)
  GOA2018SS$fleet_control$Time_varying_sel <- 1
  GOA2018SS$fleet_control$Time_varying_sel_sd <- 1
  GOA2018SS$fleet_control$Age_first_selected <- 1
  GOA2018SS$fleet_control$Sel_norm_bin <- NA # Do not normalize
  GOA2018SS$catch_data$Catch <- 1e6 # If catch is zero, sel devs are turned off

  # Run
  ss_run <- suppressMessages(
    Rceattle::fit_mod(data_list = GOA2018SS,
                      inits = NULL, # Initial parameters = 0
                      file = NULL, # Don't save
                      estimateMode = 3, # Don't estimate
                      random_rec = FALSE, # No random recruitment
                      random_sel = TRUE, # Turn on laplace for sel devs
                      msmMode = 0, # Single species mode
                      fit_control = fit_control(
                        verbose = 0))
  )

  testthat::expect_equal(as.numeric(ss_run$map$mapList$sel_dev_log_sd), c(1:3,NA,NA,NA,NA, 8:16)) # Dev sigma turned on except for not estimated fleet

  # TMB object
  testthat::expect_equal(length(unique(ss_run$obj$env$random)),  2 * 15 * nyrs)
})


testthat::test_that("Time-varying logistic selectivity likelihood", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Data
  data("GOA2018SS")
  rows_use <- which(GOA2018SS$fleet_control$Selectivity != 0 & GOA2018SS$fleet_control$Fleet_type != 0)
  GOA2018SS$fleet_control$Selectivity <-0
  GOA2018SS$fleet_control$Selectivity[rows_use] <- "Logistic" # Age-based logistic
  GOA2018SS$fleet_control$Selectivity_index <- 1:nrow(GOA2018SS$fleet_control)
  GOA2018SS$fleet_control$Time_varying_sel <- 1
  GOA2018SS$fleet_control$Time_varying_sel_sd <- 1 # Note that inf does 4*sd prior, should adapt to scale-invariant sd
  #FIXME: move to scale invariant setup
  GOA2018SS$fleet_control$Age_first_selected <- 1
  GOA2018SS$fleet_control$Sel_norm_bin <- NA # Do not normalize
  GOA2018SS$catch_data$Catch <- 1e6 # If catch is zero, sel devs are turned off

  # Specify logistic selectivity
  nyrs <- length(GOA2018SS$styr:GOA2018SS$endyr)
  inf = 10
  inf_dev <- rnorm(nyrs)
  log_slp_dev <- rnorm(nyrs)

  # Set params to logistic
  alpha = 0.5
  ss_run <- suppressMessages(
    Rceattle::fit_mod(data_list = GOA2018SS,
                      estimateMode = 3, # Don't estimate
                      msmMode = 0, # Single species mode
                      fit_control = fit_control(
                        verbose = 1))
  )
  inits <- ss_run$estimated_params
  inits$log_sel_slp[1,,1] <- log(alpha) # Females
  inits$log_sel_slp[1,,2] <- log(alpha+1) # Males
  inits$sel_inf[1,,] <- inf
  for(i in 1:dim(inits$log_sel_slp_dev[1,,,])[1]){
    for(j in 1:dim(inits$log_sel_slp_dev[1,,,])[2]){
      inits$log_sel_slp_dev[1,i,j,] <- log_slp_dev
      inits$sel_inf_dev[1,i,j,] <- inf_dev
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
  single_dv_nll <- sum(-dnorm(inf_dev, 0, 1, log = TRUE) - dnorm(log_slp_dev, 0, 4, log = TRUE))

  testthat::expect_equal(rcnll, single_dv_nll * 15) # 15 selectivities
})


testthat::test_that("Logistic selectivity blocks", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Data
  # 1) Set up simulation
  nyrs = 30
  nspp = 2
  Fmort <- c(seq(0.02, 0.3, length.out = nyrs/2), seq(0.3, 0.05, length.out = nyrs/2))
  Fmort2 <- seq(0.02, 0.3, length.out = nyrs)

  # First, simulate some data for the model
  set.seed(123)
  sim <- make_msm_test_data(
    years = 1:nyrs,
    Fmort = matrix(c(Fmort, Fmort2), nspp, nyrs, byrow = TRUE)
  )

  # Set up Rceattle data
  simData <- sim$data_list
  simData$fleet_control$Selectivity <- "Logistic" # Age-based logistic
  simData$fleet_control$Time_varying_sel <- 3 # Blocks
  simData$fleet_control$Age_first_selected <- 1
  simData$fleet_control$Sel_norm_bin <- NA

  # Specify blocks
  # - Define a function to generate 3 sequential random blocks
  assign_random_blocks <- function(years) {
    unique_years <- sort(unique(years))
    n_years <- length(unique_years)

    # Safety check: if a fleet has fewer than 3 years, assign everything to block 1
    if(n_years < 3) {
      return(rep(1, length(years)))
    }

    # Pick 2 random breakpoints from the 2nd year up to the last year.
    # These represent the starting years of Block 2 and Block 3.
    breakpoints <- sort(sample(unique_years[2:n_years], 2))

    # Assign the blocks based on where the year falls relative to the breakpoints
    blocks <- integer(length(years))
    blocks[years < breakpoints[1]] <- 1
    blocks[years >= breakpoints[1] & years < breakpoints[2]] <- 2
    blocks[years >= breakpoints[2]] <- 3

    return(blocks)
  }

  # - Apply the function by fleet
  simData$index_data <- simData$index_data |>
    dplyr::group_by(Fleet_code) |>
    dplyr::mutate(Selectivity_block = assign_random_blocks(Year),
                  Selectivity_block = as.integer(Selectivity_block)) |>
    dplyr::ungroup()

  simData$catch_data <- simData$catch_data |>
    dplyr::group_by(Fleet_code) |>
    dplyr::mutate(Selectivity_block = assign_random_blocks(Year),
                  Selectivity_block = as.integer(Selectivity_block)) |>
    dplyr::ungroup()


  # Run
  ss_run <- suppressMessages(
    Rceattle::fit_mod(data_list = simData,
                      estimateMode = 3, # Don't estimate
                      random_rec = FALSE, # No random recruitment
                      random_sel = FALSE, # Turn on laplace for sel devs
                      msmMode = 0, # Single species mode
                      fit_control = fit_control(
                        verbose = 0))
  )

  # Map
  # - Slope
  expect_all_true(is.na(as.numeric(ss_run$map$mapList$log_sel_slp)))
  # - Asympt
  expect_all_true(is.na(as.numeric(ss_run$map$mapList$sel_inf)))

  # - Devs
  # -- Ascending
  testthat::expect_equal(sum(!is.na(ss_run$map$mapList$log_sel_slp_dev[1,,,])), 4*nyrs) # slope
  testthat::expect_equal(sum(!is.na(c(ss_run$map$mapList$sel_inf_dev[1,,,]))), 4*nyrs) # asymptote

  testthat::expect_equal(length(unique(c(ss_run$map$mapList$log_sel_slp_dev[1,,,]))), 4 * 3) # slope
  testthat::expect_equal(length(unique(c(ss_run$map$mapList$sel_inf_dev[1,,,]))), 4 * 3) # asymptote

  expect_all_true(is.na(as.numeric(ss_run$map$mapList$log_sel_slp_dev[2,,,]))) # Descending slope
  expect_all_true(is.na(as.numeric(ss_run$map$mapList$sel_inf_dev[2,,,]))) # Descending asymptote

  expect_all_true(is.na(ss_run$map$mapList$sel_dev_log_sd)) # Dev sigma turned off
})

testthat::test_that("Invalid selectivity integer", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Data
  data("GOA2018SS")
  nflt <- nrow(GOA2018SS$fleet_control)
  GOA2018SS$fleet_control$Selectivity <- 9 # Not in scope

  # Run
  testthat::expect_error(
    Rceattle::fit_mod(data_list = GOA2018SS,
                      inits = inits, # Initial parameters = 0
                      file = NULL, # Don't save
                      estimateMode = 3, # Don't estimate
                      random_rec = FALSE, # No random recruitment
                      msmMode = 0, # Single species mode
                      fit_control = fit_control(
                        verbose = 0))
  )

})

testthat::test_that("Invalid selectivity string", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Data
  data("GOA2018SS")
  nflt <- nrow(GOA2018SS$fleet_control)
  GOA2018SS$fleet_control$Selectivity <- "logistic" # Not in scope

  # Run
  testthat::expect_error(
    Rceattle::fit_mod(data_list = GOA2018SS,
                      file = NULL, # Don't save
                      estimateMode = 3, # Don't estimate
                      random_rec = FALSE, # No random recruitment
                      msmMode = 0, # Single species mode
                      fit_control = fit_control(
                        verbose = 0))
  )

})



