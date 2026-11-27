# Tests for DoubleNormal selectivity (type 8)
#
# The double-normal blends a Gaussian ascending limb left of the peak with a
# Gaussian descending limb right of the peak using a steep logistic weight,
# and allows the right tail to approach a floor instead of zero:
#
#   right_floor = 1 / (1 + exp(-logit_floor))
#   val_desc(x) = right_floor + (1 - right_floor) * exp(-0.5 * ((x-peak)/sigma_desc)^2)
#   sel(x)      = (1-w) * exp(-0.5 * ((x-peak)/sigma_asc)^2)  +  w * val_desc(x)
#   w(x)        = 1 / (1 + exp(-20*(x - peak)))
#
# Parameters:
#   sel_inf[1]     = peak
#   sel_inf[2]     = logit(right_floor)  [analogous to SS3 P6 / end_logit]
#   log_sel_slp[1] = log(sigma_ascending)
#   log_sel_slp[2] = log(sigma_descending)
#
# Normalization: the model normalises sel to max = 1 when Sel_norm_bin = NA.

# Integration test (fits a CEATTLE model): skipped on CRAN to keep R CMD check
# fast; the coverage workflow runs it in full (NOT_CRAN=true). See README.md.
testthat::skip_on_cran()

# Helper: expected double-normal selectivity (R equivalent of TMB case 8)
double_normal_sel <- function(x, peak, sigma_asc, sigma_desc, logit_floor = 0) {
  right_floor <- 1 / (1 + exp(-logit_floor))
  w           <- 1 / (1 + exp(-20 * (x - peak)))
  asc_gauss   <- exp(-0.5 * ((x - peak) / sigma_asc)^2)
  desc_gauss  <- exp(-0.5 * ((x - peak) / sigma_desc)^2)
  val_desc    <- right_floor + (1 - right_floor) * desc_gauss
  val         <- (1 - w) * asc_gauss + w * val_desc
  val / max(val)   # max-normalise (matches Sel_norm_bin < 0 path)
}


testthat::test_that("Age-based double-normal: static dome shape recovers correctly", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("GOA2018SS")
  rows_use <- which(GOA2018SS$fleet_control$Selectivity != 0 & GOA2018SS$fleet_control$Fleet_type != 0)
  GOA2018SS$fleet_control$Fleet_type[-rows_use] <-0
  GOA2018SS$fleet_control$Selectivity <- "DoubleNormal"
  GOA2018SS$fleet_control$Time_varying_sel   <- 0
  GOA2018SS$fleet_control$Selectivity_index  <- seq_len(nrow(GOA2018SS$fleet_control))
  GOA2018SS$fleet_control$Age_first_selected <- 1
  GOA2018SS$fleet_control$Sel_norm_bin      <- NA  # max-normalise

  peak        <- 8
  sigma_asc   <- 3
  sigma_desc  <- 6
  logit_floor <- -10  # right_floor ~ 0 -> fully dome-shaped
  ages        <- 1:21
  expected    <- double_normal_sel(ages, peak, sigma_asc, sigma_desc, logit_floor)

  mod0 <- suppressMessages(
    Rceattle::fit_mod(data_list = GOA2018SS, inits = NULL, estimateMode = 3,
            random_rec = FALSE, msmMode = 0,
            fit_control = fit_control(verbose = 0))
  )
  inits <- mod0$estimated_params
  inits$sel_inf[1, , ]     <- peak
  inits$sel_inf[2, , ]     <- logit_floor
  inits$log_sel_slp[1, , ] <- log(sigma_asc)
  inits$log_sel_slp[2, , ] <- log(sigma_desc)

  ss_run <- suppressMessages(
    Rceattle::fit_mod(data_list = GOA2018SS, inits = inits, estimateMode = 3,
            random_rec = FALSE, msmMode = 0,
            fit_control = fit_control(verbose = 0))
  )

  # --- Map checks: all four base params estimated, no time-varying deviates ---
  # - Females/sex combined
  testthat::expect_true(all(!is.na(ss_run$map$mapList$sel_inf[1, rows_use, 1])))
  testthat::expect_true(all(!is.na(ss_run$map$mapList$sel_inf[2, rows_use, 1])))
  testthat::expect_true(all(!is.na(ss_run$map$mapList$log_sel_slp[1, rows_use, 1])))
  testthat::expect_true(all(!is.na(ss_run$map$mapList$log_sel_slp[2, rows_use, 1])))

  # - Males
  testthat::expect_true(all(!is.na(ss_run$map$mapList$sel_inf[1, 9:11, 2])))
  testthat::expect_true(all(!is.na(ss_run$map$mapList$sel_inf[2, 9:11, 2])))
  testthat::expect_true(all(!is.na(ss_run$map$mapList$log_sel_slp[1, 9:11, 2])))
  testthat::expect_true(all(!is.na(ss_run$map$mapList$log_sel_slp[2, 9:11, 2])))

  testthat::expect_true(all(is.na(c(ss_run$map$mapList$log_sel_slp_dev))))
  testthat::expect_true(all(is.na(c(ss_run$map$mapList$sel_inf_dev))))

  # --- Shape: Pollock fleet 1, ATF fleet 9 ---
  sel_out <- as.numeric(ss_run$quantities$sel_at_age[1, 1, 1:10, 1])
  testthat::expect_equal(sel_out, expected[1:10], tolerance = 1e-4)

  sel_atf <- as.numeric(ss_run$quantities$sel_at_age[9, 1, 1:21, 1])
  testthat::expect_equal(sel_atf, expected[1:21], tolerance = 1e-4)
})


testthat::test_that("Age-based double-normal: logit_floor=+10 produces near-logistic shape", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("GOA2018SS")
  rows_use <- which(GOA2018SS$fleet_control$Selectivity != 0 & GOA2018SS$fleet_control$Fleet_type != 0)
  GOA2018SS$fleet_control$Fleet_type[-rows_use] <-0
  GOA2018SS$fleet_control$Selectivity <- "DoubleNormal"
  GOA2018SS$fleet_control$Time_varying_sel   <- 0
  GOA2018SS$fleet_control$Selectivity_index  <- seq_len(nrow(GOA2018SS$fleet_control))
  GOA2018SS$fleet_control$Age_first_selected <- 1
  GOA2018SS$fleet_control$Sel_norm_bin      <- NA

  peak        <- 8
  sigma_asc   <- 3
  sigma_desc  <- 6
  logit_floor <- 10  # right_floor ~ 1 -> descending side flat at 1 (logistic)
  ages        <- 1:21
  expected    <- double_normal_sel(ages, peak, sigma_asc, sigma_desc, logit_floor)

  mod0 <- suppressMessages(
    fit_mod(data_list = GOA2018SS, inits = NULL, estimateMode = 3,
            random_rec = FALSE, msmMode = 0,
            fit_control = fit_control(verbose = 0))
  )
  inits <- mod0$estimated_params
  inits$sel_inf[1, , ]     <- peak
  inits$sel_inf[2, , ]     <- logit_floor
  inits$log_sel_slp[1, , ] <- log(sigma_asc)
  inits$log_sel_slp[2, , ] <- log(sigma_desc)

  ss_run <- suppressMessages(
    fit_mod(data_list = GOA2018SS, inits = inits, estimateMode = 3,
            random_rec = FALSE, msmMode = 0,
            fit_control = fit_control(verbose = 0))
  )

  # Use fleet 9 (ATF, species 2) which has 21 ages; fleet 1 (Pollock) only has 10.
  sel_out <- as.numeric(ss_run$quantities$sel_at_age[9, 1, 1:21, 1])
  testthat::expect_equal(sel_out, expected[1:21], tolerance = 1e-4)

  # With logit_floor=10 the selectivity should plateau near 1 at large ages
  testthat::expect_true(min(sel_out[ages >= peak]) > 0.98)
})


testthat::test_that("Length-based double-normal: static shape non-trivial", {
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
  simData$fleet_control$Selectivity           <- "DoubleNormal"
  simData$fleet_control$Selectivity_dimension <- "Length"
  simData$fleet_control$Time_varying_sel      <- 0
  simData$fleet_control$Selectivity_index     <- seq_len(nrow(simData$fleet_control))
  simData$fleet_control$Bin_first_selected    <- 1
  simData$fleet_control$Sel_norm_bin         <- NA

  peak        <- 60   # cm
  sigma_asc   <- 10
  sigma_desc  <- 20
  logit_floor <- -10  # dome-shaped
  true_sel <- sapply(1:10, function(x) double_normal_sel(x, peak, sigma_asc, sigma_desc, logit_floor))


  mod0 <- suppressMessages(
    fit_mod(data_list = simData, inits = NULL, estimateMode = 3,
            random_rec = FALSE, msmMode = 0,
            fit_control = fit_control(verbose = 0))
  )
  inits <- mod0$estimated_params
  inits$sel_inf[1, , ]     <- peak
  inits$sel_inf[2, , ]     <- logit_floor
  inits$log_sel_slp[1, , ] <- log(sigma_asc)
  inits$log_sel_slp[2, , ] <- log(sigma_desc)

  ss_run <- suppressMessages(
    fit_mod(data_list = simData, inits = inits, estimateMode = 3,
            random_rec = FALSE, msmMode = 0,
            fit_control = fit_control(verbose = 0))
  )

  # Both sel_inf[1] and sel_inf[2] estimated
  # - Females/sex combined
  testthat::expect_true(all(!is.na(ss_run$map$mapList$sel_inf[, , 1])))
  testthat::expect_true(all(!is.na(ss_run$map$mapList$log_sel_slp[, , 1])))

  # sel_at_age should be non-trivial and dome-shaped
  sel_out <- as.numeric(ss_run$quantities$sel_at_length[1, 1, , 1])
  testthat::expect_true(max(sel_out) > 0.9)   # max near 1 after normalisation
  testthat::expect_true(min(sel_out) < 0.5)   # dome-shaped: tails are low
})


testthat::test_that("IID time-varying double-normal: map has correct number of deviates", {
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
  simData$fleet_control$Selectivity              <- "DoubleNormal"
  simData$fleet_control$Time_varying_sel         <- "IID"
  simData$fleet_control$Time_varying_sel_sd <- 1
  simData$fleet_control$Selectivity_index        <- seq_len(nrow(simData$fleet_control))
  simData$fleet_control$Age_first_selected       <- 1
  simData$fleet_control$Sel_norm_bin            <- NA
  simData$catch_data$Catch <- 1e6  # non-zero catch keeps sel devs on

  n_active <- sum(simData$fleet_control$Fleet_type != "Off")
  nsex     <- 1  # simData is single-sex

  ss_run <- suppressMessages(
    fit_mod(data_list = simData, inits = NULL, estimateMode = 3,
            random_rec = FALSE, msmMode = 0,
            fit_control = fit_control(verbose = 0))
  )

  # Each active fleet contributes nyrs deviates for each of the 4 time-varying params:
  #   sel_inf_dev[1]: peak deviates
  #   sel_inf_dev[2]: right-floor logit deviates
  #   log_sel_slp_dev[1]: sigma_asc deviates
  #   log_sel_slp_dev[2]: sigma_desc deviates
  n_peak_devs    <- sum(!is.na(c(ss_run$map$mapList$sel_inf_dev[1, , , ])))
  n_floor_devs   <- sum(!is.na(c(ss_run$map$mapList$sel_inf_dev[2, , , ])))
  n_asc_sd_devs  <- sum(!is.na(c(ss_run$map$mapList$log_sel_slp_dev[1, , , ])))
  n_desc_sd_devs <- sum(!is.na(c(ss_run$map$mapList$log_sel_slp_dev[2, , , ])))

  testthat::expect_equal(n_peak_devs,    n_active * nsex * nyrs)
  testthat::expect_equal(n_floor_devs,   n_active * nsex * nyrs)
  testthat::expect_equal(n_asc_sd_devs,  n_active * nsex * nyrs)
  testthat::expect_equal(n_desc_sd_devs, n_active * nsex * nyrs)
})


testthat::test_that("IID time-varying double-normal: penalty likelihood matches R reference", {
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
  simData$fleet_control$Selectivity              <- "DoubleNormal"
  simData$fleet_control$Time_varying_sel         <- "IID"
  simData$fleet_control$Time_varying_sel_sd <- 1
  simData$fleet_control$Selectivity_index        <- seq_len(nrow(simData$fleet_control))
  simData$fleet_control$Age_first_selected       <- 1
  simData$fleet_control$Sel_norm_bin            <- NA
  simData$catch_data$Catch <- 1e6

  peak      <- 8; sigma_asc <- 3; sigma_desc <- 6; logit_floor <- 0
  peak_devs        <- rnorm(nyrs)
  floor_devs       <- rnorm(nyrs)
  asc_sd_devs      <- rnorm(nyrs)
  desc_sd_devs     <- rnorm(nyrs)

  sd_prior <- 1  # Time_varying_sel_sd

  mod0 <- suppressMessages(
    fit_mod(data_list = simData, inits = NULL, estimateMode = 3,
            random_rec = FALSE, msmMode = 0,
            fit_control = fit_control(verbose = 0))
  )
  inits <- mod0$estimated_params
  inits$sel_inf[1, , ]     <- peak
  inits$sel_inf[2, , ]     <- logit_floor
  inits$log_sel_slp[1, , ] <- log(sigma_asc)
  inits$log_sel_slp[2, , ] <- log(sigma_desc)

  n_flt <- nrow(simData$fleet_control)
  for (i in seq_len(n_flt)) {
    inits$sel_inf_dev[1, i, 1, ]     <- peak_devs
    inits$sel_inf_dev[2, i, 1, ]     <- floor_devs
    inits$log_sel_slp_dev[1, i, 1, ] <- asc_sd_devs
    inits$log_sel_slp_dev[2, i, 1, ] <- desc_sd_devs
  }

  ss_run <- suppressMessages(
    fit_mod(data_list = simData, inits = inits, estimateMode = 3,
            random_rec = FALSE, msmMode = 0,
            fit_control = fit_control(verbose = 0))
  )

  n_active <- sum(simData$fleet_control$Fleet_type != "Off")

  # NLL from TMB (jnll_comp row 6 = penalty for sel deviates, 1-indexed from R)
  rcnll <- sum(ss_run$quantities$jnll_comp[6, ])

  # R reference: for each active fleet/sex, penalise all four deviate streams.
  # Peak and floor deviates use sd=sel_dev_sd; sigma deviates use sd=4*sel_dev_sd.
  ref_nll <- n_active * (
    sum(-dnorm(peak_devs,    0, sd_prior,     log = TRUE)) +   # peak devs
    sum(-dnorm(floor_devs,   0, sd_prior,     log = TRUE)) +   # right-floor logit devs
    sum(-dnorm(asc_sd_devs,  0, 4 * sd_prior, log = TRUE)) +   # sigma_asc devs
    sum(-dnorm(desc_sd_devs, 0, 4 * sd_prior, log = TRUE))     # sigma_desc devs
  )

  testthat::expect_equal(rcnll, ref_nll, tolerance = 1e-4)
})


testthat::test_that("Block time-varying double-normal: map assigns correct block indices", {
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

  simData$fleet_control$Selectivity              <- "DoubleNormal"
  simData$fleet_control$Time_varying_sel         <- "Block"
  simData$fleet_control$Selectivity_index        <- seq_len(nrow(simData$fleet_control))
  simData$fleet_control$Age_first_selected       <- 1
  simData$fleet_control$Sel_norm_bin            <- NA

  # Assign 2 blocks: years 1..half in block 1, rest in block 2
  half <- floor(nyrs / 2)
  simData$catch_data <- simData$catch_data |>
    dplyr::mutate(Selectivity_block = dplyr::if_else(
      Year - simData$styr + 1 <= half, 1L, 2L))
  simData$index_data <- simData$index_data |>
    dplyr::mutate(Selectivity_block = dplyr::if_else(
      Year - simData$styr + 1 <= half, 1L, 2L))

  ss_run <- suppressMessages(
    fit_mod(data_list = simData, inits = NULL, estimateMode = 3,
            random_rec = FALSE, msmMode = 0,
            fit_control = fit_control(verbose = 0))
  )

  n_active <- sum(simData$fleet_control$Fleet_type != "Off")
  n_blocks <- 2

  # Base params (sel_inf[1,2], log_sel_slp[1,2]) must be fixed (NA) when blocks are active
  testthat::expect_true(all(is.na(ss_run$map$mapList$sel_inf[1, , ])))
  testthat::expect_true(all(is.na(ss_run$map$mapList$sel_inf[2, , ])))
  testthat::expect_true(all(is.na(ss_run$map$mapList$log_sel_slp[1, , ])))
  testthat::expect_true(all(is.na(ss_run$map$mapList$log_sel_slp[2, , ])))

  # Each parameter stream (peak, floor, asc-SD, desc-SD) should have n_active * n_blocks unique indices
  n_peak_devs  <- length(unique(na.omit(c(ss_run$map$mapList$sel_inf_dev[1, , , ]))))
  n_floor_devs <- length(unique(na.omit(c(ss_run$map$mapList$sel_inf_dev[2, , , ]))))
  n_asc_devs   <- length(unique(na.omit(c(ss_run$map$mapList$log_sel_slp_dev[1, , , ]))))
  n_desc_devs  <- length(unique(na.omit(c(ss_run$map$mapList$log_sel_slp_dev[2, , , ]))))

  testthat::expect_equal(n_peak_devs,  n_active * n_blocks)
  testthat::expect_equal(n_floor_devs, n_active * n_blocks)
  testthat::expect_equal(n_asc_devs,   n_active * n_blocks)
  testthat::expect_equal(n_desc_devs,  n_active * n_blocks)
})
