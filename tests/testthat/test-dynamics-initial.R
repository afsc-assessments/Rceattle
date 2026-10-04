# Integration test (fits a CEATTLE model): skipped on CRAN to keep R CMD check
# fast; the coverage workflow runs it in full (NOT_CRAN=true). See README.md.
testthat::skip_on_cran()

testthat::test_that("Invalid initMode", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Load helper and create small deterministic test data
  #source(file.path("tests", "testthat", "helpers.R"))
  dat <- make_test_data(nyrs = 20, nages = 5, seed = 123)

  # initMode = 5 (OffsetEquilibrium) and 6 (FishedNonEquilibriumSelected, added
  # in 5.46.0) are both valid; 7 is out of the 0:6 range and must error.
  testthat::expect_error( Rceattle::fit_mod(data_list = dat,
                                              initMode = 7,
                                              estimateMode = 3)
  )

  # Run
  testthat::expect_error( Rceattle::fit_mod(data_list = dat,
                                            initMode = "wrong",
                                            estimateMode = 3)
  )
})

testthat::test_that("Free parameter initMode", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Load helper and create small deterministic test data
  #source(file.path("tests", "testthat", "helpers.R"))
  nages = 5
  dat <- make_test_data(nyrs = 20, nages = nages, seed = 123)

  # Run integer
  mod <- Rceattle::fit_mod(data_list = dat,
                           initMode = 0,
                           estimateMode = 3)

  natage1 <- mod$quantities$N_at_age[1,1,,1]
  testthat::expect_equal(as.numeric(natage1), c(exp(9), rep(1, nages-1)))

  testthat::expect_equal(c(mod$map$mapList$init_dev), c(as.numeric(1:(nages-1)), NA))
  testthat::expect_equal(as.numeric(mod$map$mapList$log_Finit), as.numeric(rep(NA, length(mod$map$mapList$log_Finit))))

  # Run string
  mod <- Rceattle::fit_mod(data_list = dat,
                           initMode = "FreeParams",
                           estimateMode = 3)

  natage1 <- mod$quantities$N_at_age[1,1,,1]
  testthat::expect_equal(as.numeric(natage1), c(exp(9), rep(1, nages-1)))

  testthat::expect_equal(c(mod$map$mapList$init_dev), c(as.numeric(1:(nages-1)), NA))
  testthat::expect_equal(as.numeric(mod$map$mapList$log_Finit), as.numeric(rep(NA, length(mod$map$mapList$log_Finit))))
})

testthat::test_that("Equilibrium initMode", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Load helper and create small deterministic test data
  #source(file.path("tests", "testthat", "helpers.R"))
  dat <- make_test_data(nyrs = 20, nages = 5, seed = 123)

  # Run integer
  mod <- Rceattle::fit_mod(data_list = dat,
                           initMode = 1,
                           estimateMode = 3)

  natage1 <- mod$quantities$N_at_age[1,1,,1]
  testthat::expect_equal(as.numeric(natage1), c(exp(9-0.2 * 0:3), exp(9)/(1-exp(-0.2)) * exp(-0.2*4)))

  testthat::expect_equal(c(mod$map$mapList$init_dev), as.numeric(rep(NA, length(mod$map$mapList$init_dev))))
  testthat::expect_equal(as.numeric(mod$map$mapList$log_Finit), as.numeric(rep(NA, length(mod$map$mapList$log_Finit))))

  # Run string
  mod <- Rceattle::fit_mod(data_list = dat,
                           initMode = "Equilibrium",
                           estimateMode = 3)

  natage1 <- mod$quantities$N_at_age[1,1,,1]
  testthat::expect_equal(as.numeric(natage1), c(exp(9-0.2 * 0:3), exp(9)/(1-exp(-0.2)) * exp(-0.2*4)))

  testthat::expect_equal(c(mod$map$mapList$init_dev), as.numeric(rep(NA, length(mod$map$mapList$init_dev))))
  testthat::expect_equal(as.numeric(mod$map$mapList$log_Finit), as.numeric(rep(NA, length(mod$map$mapList$log_Finit))))
})


testthat::test_that("NonEquilibrium initMode", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Load helper and create small deterministic test data
  #source(file.path("tests", "testthat", "helpers.R"))
  nages = 5
  dat <- make_test_data(nyrs = 20, nages = nages, seed = 123)

  # Run integer
  mod <- Rceattle::fit_mod(data_list = dat,
                           initMode = 2,
                           estimateMode = 3)

  natage1 <- mod$quantities$N_at_age[1,1,,1]
  testthat::expect_equal(as.numeric(natage1), c(exp(9-0.2 * 0:3), exp(9)/(1-exp(-0.2)) * exp(-0.2*4)))

  testthat::expect_equal(c(mod$map$mapList$init_dev), c(as.numeric(1:(nages-1)), NA))
  testthat::expect_equal(as.numeric(mod$map$mapList$log_Finit), as.numeric(rep(NA, length(mod$map$mapList$log_Finit))))

  # Run string
  mod <- Rceattle::fit_mod(data_list = dat,
                           initMode = "NonEquilibrium",
                           estimateMode = 3)

  natage1 <- mod$quantities$N_at_age[1,1,,1]
  testthat::expect_equal(as.numeric(natage1), c(exp(9-0.2 * 0:3), exp(9)/(1-exp(-0.2)) * exp(-0.2*4)))

  testthat::expect_equal(c(mod$map$mapList$init_dev), c(as.numeric(1:(nages-1)), NA))
  testthat::expect_equal(as.numeric(mod$map$mapList$log_Finit), as.numeric(rep(NA, length(mod$map$mapList$log_Finit))))

  # - Likelihood
  testthat::expect_equal(as.numeric(mod$quantities$jnll_comp[10,1]), -dnorm(0, 1/2, sd = 1, log = TRUE) * (nages - 1)) # Lognormal bias correction with sd = 1
})


testthat::test_that("FishedNonEquilibrium initMode", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Load helper and create small deterministic test data
  #source(file.path("tests", "testthat", "helpers.R"))
  nages = 5
  dat <- make_test_data(nyrs = 20, nages = nages, seed = 123)

  # Run integer
  mod <- Rceattle::fit_mod(data_list = dat,
                           initMode = 3,
                           estimateMode = 3)

  natage1 <- mod$quantities$N_at_age[1,1,,1]
  Finit <- as.numeric(exp(mod$estimated_params$log_Finit))
  testthat::expect_equal(as.numeric(natage1), c(exp(9-(0.2+Finit) * 0:3), exp(9) * exp(-(0.2+Finit)*4)/(1-exp(-0.2 - Finit))))

  testthat::expect_equal(c(mod$map$mapList$init_dev), c(as.numeric(1:(nages-1)), NA))
  testthat::expect_equal(as.numeric(mod$map$mapList$log_Finit), 1)

  # - Likelihood
  testthat::expect_equal(as.numeric(mod$quantities$jnll_comp[10,1]), -dnorm(0, 1/2, sd = 1, log = TRUE) * (nages - 1)) # Lognormal bias correction with sd = 1

  # Run string
  mod <- Rceattle::fit_mod(data_list = dat,
                           initMode = "FishedNonEquilibrium",
                           estimateMode = 3)

  natage1 <- mod$quantities$N_at_age[1,1,,1]
  Finit <- as.numeric(exp(mod$estimated_params$log_Finit))
  testthat::expect_equal(as.numeric(natage1), c(exp(9-(0.2+Finit) * 0:3), exp(9) * exp(-(0.2+Finit)*4)/(1-exp(-0.2 - Finit))))

  testthat::expect_equal(c(mod$map$mapList$init_dev), c(as.numeric(1:(nages-1)), NA))
  testthat::expect_equal(as.numeric(mod$map$mapList$log_Finit), 1)

  # - Likelihood
  testthat::expect_equal(as.numeric(mod$quantities$jnll_comp[10,1]), -dnorm(0, 1/2, sd = 1, log = TRUE) * (nages - 1)) # Lognormal bias correction with sd = 1
})


testthat::test_that("FishedNonEquilibriumScaled initMode", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Load helper and create small deterministic test data
  #source(file.path("tests", "testthat", "helpers.R"))
  nages = 5
  dat <- make_test_data(nyrs = 20, nages = nages, seed = 123)

  # Run integer
  mod <- Rceattle::fit_mod(data_list = dat,
                           initMode = 4,
                           estimateMode = 3)

  natage1 <- mod$quantities$N_at_age[1,1,,1]
  Finit <- as.numeric(exp(mod$estimated_params$log_Finit))
  testthat::expect_equal(as.numeric(natage1), c(exp(9), exp(9-(0.2) * 1:3)* exp(-Finit), exp(9) * exp(-(0.2)*4)/(1-exp(-0.2 - Finit))* exp(-Finit)))

  testthat::expect_equal(c(mod$map$mapList$init_dev), c(as.numeric(1:(nages-1)), NA))
  testthat::expect_equal(as.numeric(mod$map$mapList$log_Finit), 1)

  # - Likelihood
  testthat::expect_equal(as.numeric(mod$quantities$jnll_comp[10,1]), -dnorm(0, 1/2, sd = 1, log = TRUE) * (nages - 1)) # Lognormal bias correction with sd = 1


  # Run string
  mod <- Rceattle::fit_mod(data_list = dat,
                           initMode = "FishedNonEquilibriumScaled",
                           estimateMode = 3)

  natage1 <- mod$quantities$N_at_age[1,1,,1]
  Finit <- as.numeric(exp(mod$estimated_params$log_Finit))
  testthat::expect_equal(as.numeric(natage1), c(exp(9), exp(9-(0.2) * 1:3)* exp(-Finit), exp(9) * exp(-(0.2)*4)/(1-exp(-0.2 - Finit))* exp(-Finit)))

  testthat::expect_equal(c(mod$map$mapList$init_dev), c(as.numeric(1:(nages-1)), NA))
  testthat::expect_equal(as.numeric(mod$map$mapList$log_Finit), 1)

  # - Likelihood
  testthat::expect_equal(as.numeric(mod$quantities$jnll_comp[10,1]), -dnorm(0, 1/2, sd = 1, log = TRUE) * (nages - 1)) # Lognormal bias correction with sd = 1
})
