
# Tests for "LogisticPM" selectivity (type 11): the ADMB ("pm"/AMAK) bottom-trawl
# survey form. A 2-parameter logistic evaluated at the MID-AGE age_vector(j)=j+0.5
# with MULTIPLICATIVE time-varying deviates on slope and inflection, PLUS a free
# first-bin (age-1) log-selectivity independent of the logistic:
#   sel(age) = 1 / (1 + exp(-slp*exp(slp_dev) * ((age+0.5) - a50*exp(a50_dev))))
#   sel(1)   = exp(age_one * exp(age_one_dev))                       [overrides bin 0]
# Parameters reuse the logistic slots: log_sel_slp[1] = log(slope),
# sel_inf[1] = a50 (multiplicative dev in *_dev[1]); the free age-1 base/deviate
# live in the unused descending-limb slots sel_inf[2] / sel_inf_dev[2].
# Penalty (jnll_comp(5) -> R row 6) is a random walk (norm2 of first differences)
# on each deviate, weighted by Sel_curve_pen1/2/3 (slope / inflection / age-1).

# Integration test (fits a CEATTLE model): skipped on CRAN to keep R CMD check
# fast; the coverage workflow runs it in full (NOT_CRAN=true). See README.md.
testthat::skip_on_cran()

testthat::test_that("LogisticPM construction matches the AMAK mid-age logistic + free age-1", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("GOA2018SS")
  rows_use <- which(GOA2018SS$fleet_control$Selectivity != 0)

  GOA2018SS$fleet_control$Selectivity[rows_use]               <- "LogisticPM"
  GOA2018SS$fleet_control$Age_first_selected[rows_use]        <- 1   # age-1 selected (free)
  GOA2018SS$fleet_control$Sel_norm_bin[rows_use]             <- NA  # no normalization
  GOA2018SS$fleet_control$Time_varying_sel[rows_use]          <- 0   # time-invariant
  GOA2018SS$fleet_control$Sel_curve_pen1[rows_use]            <- 0
  GOA2018SS$fleet_control$Sel_curve_pen2[rows_use]            <- 0
  GOA2018SS$fleet_control$Sel_curve_pen3[rows_use]            <- 0

  mod0 <- suppressMessages(fit_mod(data_list = GOA2018SS, inits = NULL, estimateMode = 3,
                                   random_rec = FALSE, msmMode = 0,
                                   fit_control = fit_control(verbose = 0)))
  inits <- mod0$estimated_params

  slope_v <- 1.2; a50_v <- 5.0; age1_log <- -1.0
  inits$log_sel_slp[1, rows_use, ] <- log(slope_v)
  inits$sel_inf[1, rows_use, ]     <- a50_v
  inits$sel_inf[2, rows_use, ]     <- age1_log

  run <- suppressMessages(fit_mod(data_list = GOA2018SS, inits = inits, estimateMode = 3,  # quantities are reported at the input parameters
                                  random_rec = FALSE, msmMode = 0,
                                  fit_control = fit_control(verbose = 0)))

  flt   <- rows_use[1]
  nages <- GOA2018SS$nages[1]
  sel   <- run$quantities$sel_at_age[flt, 1, 1:nages, 1]

  expected <- numeric(nages)
  for (a in 1:nages) expected[a] <- 1 / (1 + exp(-slope_v * ((a + 0.5) - a50_v)))
  expected[1] <- exp(age1_log)   # free age-1 override

  testthat::expect_equal(as.numeric(sel), expected, tolerance = 1e-6)
})


testthat::test_that("LogisticPM penalty = realized-logsel RW (age-range, start-year) + age-1 dev RW", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("GOA2018SS")
  rows_use <- which(GOA2018SS$fleet_control$Selectivity != 0)
  w_sel <- 2; w_a1 <- 8           # realized-logsel RW weight (ctrl_flag26); age-1 dev RW weight
  start_year <- GOA2018SS$styr + 3 # exercise the start-year offset

  GOA2018SS$fleet_control$Selectivity[rows_use]               <- "LogisticPM"
  GOA2018SS$fleet_control$Age_first_selected[rows_use]        <- 1
  GOA2018SS$fleet_control$Sel_norm_bin[rows_use]             <- NA  # default age-range = whole selected range
  GOA2018SS$fleet_control$Sel_norm_bin_upper[rows_use]             <- NA
  GOA2018SS$fleet_control$Sel_start_year[rows_use]            <- start_year
  GOA2018SS$fleet_control$Time_varying_sel[rows_use]          <- "RandomWalk"
  GOA2018SS$fleet_control$Time_varying_sel_sd[rows_use] <- 1
  GOA2018SS$fleet_control$Sel_curve_pen1[rows_use]            <- w_sel
  GOA2018SS$fleet_control$Sel_curve_pen2[rows_use]            <- 0
  GOA2018SS$fleet_control$Sel_curve_pen3[rows_use]            <- w_a1

  mod0 <- suppressMessages(fit_mod(data_list = GOA2018SS, inits = NULL, estimateMode = 3,
                                   random_rec = FALSE, msmMode = 0,
                                   fit_control = fit_control(verbose = 0)))
  inits <- mod0$estimated_params

  flt   <- rows_use[1]
  nyrs  <- dim(inits$sel_inf_dev)[4]
  set.seed(42)
  slp_dev <- rnorm(nyrs, 0, 0.2); slp_dev[1] <- 0   # first dev fixed (RandomWalk)
  a50_dev <- rnorm(nyrs, 0, 0.2); a50_dev[1] <- 0
  a1_dev  <- rnorm(nyrs, 0, 0.2); a1_dev[1]  <- 0

  inits$log_sel_slp[1, rows_use, ] <- log(1.2)
  inits$sel_inf[1, rows_use, ]     <- 5.0
  inits$sel_inf[2, rows_use, ]     <- -1.0
  for (r in rows_use) {
    inits$log_sel_slp_dev[1, r, 1, ] <- slp_dev
    inits$sel_inf_dev[1, r, 1, ]     <- a50_dev
    inits$sel_inf_dev[2, r, 1, ]     <- a1_dev
  }

  run <- suppressMessages(fit_mod(data_list = GOA2018SS, inits = inits, estimateMode = 3,  # quantities are reported at the input parameters
                                  random_rec = FALSE, msmMode = 0,
                                  fit_control = fit_control(verbose = 0)))

  # Expected: penalty starts the year AFTER start_year (0-based index start_idx+1)
  start_idx <- start_year - GOA2018SS$styr            # 0-based
  nages     <- GOA2018SS$nages[1]
  ls        <- log(run$quantities$sel_at_age[flt, 1, 1:nages, ])  # [age, yr]
  yr_rng    <- (start_idx + 2):nyrs                    # R 1-based years contributing a first-difference
  pen_sel   <- w_sel * sum(sapply(yr_rng, function(y) sum((ls[, y] - ls[, y - 1])^2)))
  pen_a1    <- w_a1  * sum(sapply(yr_rng, function(y) (a1_dev[y] - a1_dev[y - 1])^2))

  # jnll_comp(5) (C++) -> R row 6
  testthat::expect_equal(as.numeric(run$quantities$jnll_comp[6, flt]), pen_sel + pen_a1, tolerance = 1e-4)
})


testthat::test_that("LogisticPM rejects non-random-walk time variation", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("GOA2018SS")
  rows_use <- which(GOA2018SS$fleet_control$Selectivity != 0)
  GOA2018SS$fleet_control$Selectivity[rows_use]      <- "LogisticPM"
  GOA2018SS$fleet_control$Time_varying_sel[rows_use] <- "IID"   # invalid for LogisticPM

  testthat::expect_error(
    suppressMessages(fit_mod(data_list = GOA2018SS, inits = NULL, file = NULL,
                             estimateMode = 3, random_rec = FALSE, msmMode = 0,
                             fit_control = fit_control(verbose = 0))),
    regexp = "RandomWalk"
  )
})

testthat::test_that("LogisticPM defaults the free age-1 log-selectivity slot to a log scale, not the inflection-age default", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Regression test. sel_inf[2] is normally the descending-limb inflection *age*
  # and build_params defaults it to 10 for every fleet. LogisticPM repurposes that
  # slot as the free first-bin (age-1) LOG-selectivity, so the shared default gave
  # age-1 selectivity = exp(10) = 22026 (ADMB's sel_age_one_bts = -3.19 -> 0.04).
  # With age-1 selected that swamps the survey-index prediction and the optimizer
  # diverges (gradient blow-up on sel_inf), so the default must be a log-scale 0.
  data("GOA2018SS")
  rows_use <- which(GOA2018SS$fleet_control$Selectivity != 0)
  GOA2018SS$fleet_control$Selectivity[rows_use]        <- "LogisticPM"
  GOA2018SS$fleet_control$Age_first_selected[rows_use] <- 1   # age-1 selected (free)
  GOA2018SS$fleet_control$Sel_norm_bin[rows_use]      <- NA  # no normalization
  GOA2018SS$fleet_control$Time_varying_sel[rows_use]   <- 0   # time-invariant

  mod0 <- suppressMessages(fit_mod(data_list = GOA2018SS, inits = NULL, estimateMode = 3,
                                   random_rec = FALSE, msmMode = 0,
                                   fit_control = fit_control(verbose = 0)))

  # Default init is on the log scale (0), NOT the inflection-age default (10)
  testthat::expect_equal(unname(mod0$estimated_params$sel_inf[2, rows_use[1], 1]), 0)

  # ...so realized age-1 selectivity is exp(0) = 1, and the curve has no
  # exp() underflow to exactly zero. sel_at_age is dimensioned to the max nages
  # across species, so only this fleet's own species' ages carry a curve --
  # bins beyond nages[sp] are structurally 0 and are not part of the check.
  sp <- GOA2018SS$fleet_control$Species[rows_use[1]]
  s  <- mod0$obj$report()$sel_at_age[rows_use[1], 1, seq_len(GOA2018SS$nages[sp]), 1]
  testthat::expect_equal(unname(s[1]), 1)
  testthat::expect_true(all(is.finite(s)))
  testthat::expect_true(all(s > 0))
})
