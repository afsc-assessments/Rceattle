# Integration test (fits a CEATTLE model): skipped on CRAN to keep R CMD check
# fast; the coverage workflow runs it in full (NOT_CRAN=true). See README.md.
testthat::skip_on_cran()

testthat::test_that("Fixed catchability", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Data
  data("GOA2018SS")
  GOA2018SS$fleet_control$Catchability <- 0 # Fixed

  # Set params
  mod0 <- suppressMessages( fit_mod(data_list = GOA2018SS, inits = NULL, estimateMode = 3, random_rec = FALSE, msmMode = 0, fit_control = fit_control(verbose = 0)) )
  inits <- mod0$estimated_params
  inits$index_log_q[] <- 0

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
  nflt <- nrow(GOA2018SS$fleet_control)
  testthat::expect_equal(as.numeric(ss_run$map$mapList$index_log_q), as.numeric(rep(NA, nflt)))
  testthat::expect_equal(as.numeric(ss_run$map$mapList$index_log_q), as.numeric(rep(NA, nflt)))
  testthat::expect_equal(as.numeric(ss_run$map$mapList$index_q_beta), as.numeric(rep(NA, length(ss_run$map$mapList$index_q_beta))))
  testthat::expect_equal(as.numeric(ss_run$map$mapList$index_q_dev), as.numeric(rep(NA, length(ss_run$map$mapList$index_q_dev))))
  testthat::expect_equal(as.numeric(ss_run$map$mapList$index_q_dev_log_sd), as.numeric(rep(NA, length(ss_run$map$mapList$index_q_dev_log_sd))))
  testthat::expect_equal(as.numeric(ss_run$map$mapList$index_q_log_sd), as.numeric(rep(NA, length(ss_run$map$mapList$index_q_log_sd))))

  # Check q
  # - Pollock
  testthat::expect_equal(as.numeric(ss_run$quantities$index_q), as.numeric(rep(1, length(ss_run$quantities$index_q))))
})

testthat::test_that("Estimated catchability", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Data
  data("GOA2018SS")
  nflt <- nrow(GOA2018SS$fleet_control)
  GOA2018SS$fleet_control$Catchability <- 1 # Estimated
  GOA2018SS$fleet_control$Time_varying_q <- "Off" # Estimated
  GOA2018SS$fleet_control$Catchability_index <- 1:nflt

  # Set params
  mod0 <- suppressMessages( fit_mod(data_list = GOA2018SS, inits = NULL, estimateMode = 3, random_rec = FALSE, msmMode = 0, fit_control = fit_control(verbose = 0)) )
  inits <- mod0$estimated_params
  inits$index_log_q[] <- 0

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

  # Map. Catchability follows the DATA, not the fleet type: a fleet is given a q
  # only where it carries fitted index_data. GOA2018SS fleet 10
  # (ATF_bottom_trawl_length_comp) is a composition-only survey with no index
  # rows, so it has no q to estimate -- a catchability with no index to inform it
  # is a flat direction. Before 5.9.0 every Survey got one regardless.
  fleets <- 1:nflt
  fleets[ss_run$data_list$fleet_control$Fleet_type != "Survey"] <- NA
  fleets[!(seq_len(nflt) %in% Rceattle:::.fleets_with_index(ss_run$data_list))] <- NA
  testthat::expect_equal(as.numeric(ss_run$map$mapList$index_log_q), fleets)
  # ...and that really did drop a fleet, or the line above is doing nothing.
  testthat::expect_true(
    any(ss_run$data_list$fleet_control$Fleet_type == "Survey" &
          !(seq_len(nflt) %in% Rceattle:::.fleets_with_index(ss_run$data_list))))
  testthat::expect_equal(as.numeric(ss_run$map$mapList$index_q_beta), as.numeric(rep(NA, length(ss_run$map$mapList$index_q_beta))))

  testthat::expect_equal(as.numeric(ss_run$map$mapList$index_q_dev), as.numeric(rep(NA, length(ss_run$map$mapList$index_q_dev))))
  testthat::expect_equal(as.numeric(ss_run$map$mapList$index_q_dev_log_sd), as.numeric(rep(NA, length(ss_run$map$mapList$index_q_dev_log_sd))))
  testthat::expect_equal(as.numeric(ss_run$map$mapList$index_q_log_sd), as.numeric(rep(NA, length(ss_run$map$mapList$index_q_log_sd))))

  # Check q
  # - Pollock
  testthat::expect_equal(as.numeric(ss_run$quantities$index_q), as.numeric(rep(1, length(ss_run$quantities$index_q))))
})

testthat::test_that("Invalid catchability", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Data
  data("GOA2018SS")
  nflt <- nrow(GOA2018SS$fleet_control)
  GOA2018SS$fleet_control$Catchability <- 9 # Not in scope
  GOA2018SS$fleet_control$Catchability_index <- 1:nflt

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

