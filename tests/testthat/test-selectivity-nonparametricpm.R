
# Tests for "NonParametricPM" selectivity (type 9): identical construction to
# the Ianelli "NonParametric" (type 2) selectivity, but with the ADMB ("pm"/
# AMAK) Selectivity_Likelihood penalty form:
#   jnll_comp(4) [R row 5] = decreasing-only penalty  (ADMB sel_like(1),
#                            weight Sel_curve_pen1), differentiable (|d|+d)/2;
#   jnll_comp(5) [R row 6] = curvature (2nd-diff, weight Sel_curve_pen2)
#                          + bare-SSQ random walk over all ages (no dnorm const)
#                          + dev-magnitude on coefficient increments (Sel_curve_pen3).
# (The type-2 penalty instead keeps curvature in jnll_comp(4) and adds the
#  2*avg_sel^2 normalization term, and uses dnorm() for the random walk.)

# Integration test (fits a CEATTLE model): skipped on CRAN to keep R CMD check
# fast; the coverage workflow runs it in full (NOT_CRAN=true). See README.md.
testthat::skip_on_cran()

testthat::test_that("NonParametricPM construction matches NonParametric", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("GOA2018SS")
  rows_use <- which(GOA2018SS$fleet_control$Selectivity != 0)
  n_sel_bins <- 8
  log_selcoffs <- rnorm(n_sel_bins)

  setup <- function(dat, sel){
    dat$fleet_control$Selectivity[rows_use]                <- sel
    dat$fleet_control$N_sel_bins[rows_use]                 <- n_sel_bins
    dat$fleet_control$Age_first_selected[rows_use]         <- 1
    dat$fleet_control$Sel_curve_pen1[rows_use]             <- 5
    dat$fleet_control$Sel_curve_pen2[rows_use]             <- 10
    dat$fleet_control$Sel_curve_pen3[rows_use]             <- 0
    dat$fleet_control$Sel_norm_bin[rows_use]              <- NA
    dat$fleet_control$Time_varying_sel[rows_use]           <- 0   # time-invariant
    dat$fleet_control$Time_varying_sel_sd[rows_use]  <- 1
    dat
  }
  run_with <- function(dat){
    m0 <- suppressMessages(fit_mod(data_list = dat, inits = NULL, estimateMode = 3,
                                   random_rec = FALSE, msmMode = 0,
                                   fit_control = fit_control(verbose = 0)))
    inits <- m0$estimated_params
    inits$sel_coff[, 1, 1:n_sel_bins] <- rep(log_selcoffs, each = dim(inits$sel_coff)[1])
    suppressMessages(fit_mod(data_list = dat, inits = inits, estimateMode = 3,
                             random_rec = FALSE, msmMode = 0,
                             fit_control = fit_control(verbose = 0)))
  }

  np <- run_with(setup(GOA2018SS, "NonParametric"))
  pm <- run_with(setup(GOA2018SS, "NonParametricPM"))

  # Selectivity-at-age must be identical (shared construction; only the penalty differs)
  testthat::expect_equal(as.numeric(pm$quantities$sel_at_age),
                         as.numeric(np$quantities$sel_at_age), tolerance = 1e-10)
})


testthat::test_that("NonParametricPM penalty matches the ADMB form (time-invariant)", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("GOA2018SS")
  rows_use <- which(GOA2018SS$fleet_control$Selectivity != 0)
  n_sel_bins <- 8
  log_selcoffs <- rnorm(n_sel_bins)

  GOA2018SS$fleet_control$Selectivity[rows_use]               <- "NonParametricPM"
  GOA2018SS$fleet_control$N_sel_bins[rows_use]                <- n_sel_bins
  GOA2018SS$fleet_control$Age_first_selected[rows_use]        <- 1
  GOA2018SS$fleet_control$Sel_curve_pen1[rows_use]            <- 5    # decreasing weight
  GOA2018SS$fleet_control$Sel_curve_pen2[rows_use]            <- 10   # curvature weight
  GOA2018SS$fleet_control$Sel_curve_pen3[rows_use]            <- 0
  GOA2018SS$fleet_control$Sel_norm_bin[rows_use]             <- NA
  GOA2018SS$fleet_control$Time_varying_sel[rows_use]          <- 0
  GOA2018SS$fleet_control$Time_varying_sel_sd[rows_use] <- 1

  mod0 <- suppressMessages(fit_mod(data_list = GOA2018SS, inits = NULL, estimateMode = 3,
                                   random_rec = FALSE, msmMode = 0,
                                   fit_control = fit_control(verbose = 0)))
  inits <- mod0$estimated_params
  inits$sel_coff[, 1, 1:n_sel_bins] <- rep(log_selcoffs, each = dim(inits$sel_coff)[1])

  ss_run <- suppressMessages(fit_mod(data_list = GOA2018SS, inits = inits, estimateMode = 3,
                                     random_rec = FALSE, msmMode = 0,
                                     fit_control = fit_control(verbose = 0)))

  # Pollock fishery selectivity (sel_at_age index 1, jnll_comp column 8 - same as
  # the type-2 test test-selectivity-nonparametric.R)
  log_sel_fsh <- c(log_selcoffs, rep(log_selcoffs[n_sel_bins], GOA2018SS$nages[1] - n_sel_bins))
  log_sel_fsh <- log_sel_fsh - log(mean(exp(log_sel_fsh)))

  # Decreasing penalty -> jnll_comp(4) [R row 5], weight Sel_curve_pen1 = 5
  difftmp <- -diff(log_sel_fsh)                 # decreasing is positive
  difftmp <- (abs(difftmp) + difftmp) / 2       # differentiable max(., 0)
  pen_dec <- 5 * sum(difftmp^2)

  # Curvature penalty -> jnll_comp(5) [R row 6], weight Sel_curve_pen2 = 10
  pen_curv <- 10 * sum(diff(diff(log_sel_fsh))^2)

  flt <- 8
  # ADMB-form split: decreasing in row 5, curvature in row 6, NO 2*avg_sel^2 term
  testthat::expect_equal(as.numeric(ss_run$quantities$jnll_comp[5, flt]), pen_dec,  tolerance = 1e-4)
  testthat::expect_equal(as.numeric(ss_run$quantities$jnll_comp[6, flt]), pen_curv, tolerance = 1e-4)
})


testthat::test_that("Length-based NonParametricPM penalty spans the length bins, not the ages", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Length-based selectivity here has nlengths (20) bins against nages (15), so
  # the shape/curvature penalties must span all 20 length bins. One free
  # coefficient per bin keeps bins 16-20 informative (no plateau), making the
  # check sensitive to the range beyond nages.
  set.seed(11)
  sim <- make_msm_test_data(use_size_sel = TRUE)
  dat <- sim$data_list
  nlengths <- dat$nlengths[1]                      # 20
  n_sel_bins <- nlengths                           # one coefficient per length bin
  coffs <- rnorm(n_sel_bins)

  dat$fleet_control$Selectivity               <- "NonParametricPM"
  dat$fleet_control$Selectivity_dimension     <- "Length"
  dat$fleet_control$N_sel_bins                <- n_sel_bins
  dat$fleet_control$Bin_first_selected        <- 1
  dat$fleet_control$Sel_curve_pen1            <- 5     # decreasing weight
  dat$fleet_control$Sel_curve_pen2            <- 10    # curvature weight
  dat$fleet_control$Sel_curve_pen3            <- 0
  dat$fleet_control$Sel_norm_bin             <- NA
  dat$fleet_control$Time_varying_sel          <- 0
  dat$fleet_control$Time_varying_sel_sd <- 1

  flt <- 1                                          # Survey Species 1
  mod0 <- suppressMessages(fit_mod(data_list = dat, inits = NULL, estimateMode = 3,
                                   random_rec = FALSE, msmMode = 0,
                                   fit_control = fit_control(verbose = 0)))
  inits <- mod0$estimated_params
  inits$sel_coff[flt, 1, 1:n_sel_bins] <- coffs
  ss <- suppressMessages(fit_mod(data_list = dat, inits = inits, estimateMode = 3,
                                 random_rec = FALSE, msmMode = 0,
                                 fit_control = fit_control(verbose = 0)))

  # Realized log-selectivity over all length bins: coefficients plateau past
  # N_sel_bins, then centred so mean(exp(.)) = 1 (matches case 9 construction).
  log_sel <- c(coffs, rep(coffs[n_sel_bins], nlengths - n_sel_bins))
  log_sel <- log_sel - log(mean(exp(log_sel)))

  d       <- -diff(log_sel); d <- (abs(d) + d) / 2  # differentiable max(decreasing, 0)
  pen_dec <- 5  * sum(d^2)                           # jnll row 5, over 19 length pairs
  pen_curv <- 10 * sum(diff(diff(log_sel))^2)        # jnll row 6, over 20 length bins

  testthat::expect_equal(as.numeric(ss$quantities$jnll_comp[5, flt]), pen_dec,  tolerance = 1e-4)
  testthat::expect_equal(as.numeric(ss$quantities$jnll_comp[6, flt]), pen_curv, tolerance = 1e-4)
})


testthat::test_that("NonParametricPM rejects non-random-walk time variation", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("GOA2018SS")
  rows_use <- which(GOA2018SS$fleet_control$Selectivity != 0)
  GOA2018SS$fleet_control$Selectivity[rows_use]      <- "NonParametricPM"
  GOA2018SS$fleet_control$N_sel_bins[rows_use]       <- 8
  GOA2018SS$fleet_control$Sel_curve_pen1[rows_use]   <- 5
  GOA2018SS$fleet_control$Sel_curve_pen2[rows_use]   <- 10
  GOA2018SS$fleet_control$Time_varying_sel[rows_use] <- "IID"   # invalid for non-parametric

  testthat::expect_error(
    suppressMessages(fit_mod(data_list = GOA2018SS, inits = NULL, file = NULL,
                             estimateMode = 3, random_rec = FALSE, msmMode = 0,
                             fit_control = fit_control(verbose = 0))),
    regexp = "RandomWalk"
  )
})
