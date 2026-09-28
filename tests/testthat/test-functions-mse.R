testthat::test_that("Test MSE - Tier 3 w no uncertainty", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")
  testthat::skip_on_cran()


  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Data ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  data(BS2017SS) # ?BS2017SS for more information on the data
  BS2017SS$projyr <- 2040
  BS2017SS$fleet_control$Proj_F_proportion <-rep(1,7)


  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Operating Model ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Single-species with fixed M
  ss_run <- Rceattle::fit_mod(data_list = BS2017SS,
                              inits = NULL, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 1, # Estimate hindcast only
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              fit_control = fit_control(
                                phase = TRUE,
                                getsd = FALSE,
                                verbose = 0))

  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # NPFMC Tier 3 ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  ss_run_Tier3 <- suppressWarnings(
    Rceattle::fit_mod(data_list = BS2017SS,
                      inits = ss_run$estimated_params, # Initial parameters from ss_run
                      estimateMode = 2, # Run projection only
                      HCR = build_hcr(HCR = 5, # Tier3 HCR
                                      Ftarget = 0.4, # F40%
                                      Flimit = 0.35, # F35%
                                      Plimit = 0.2, # No fishing when SB<SB20
                                      Alpha = 0.05),
                      msmMode = 0, # Single species mode
                      fit_control = fit_control(
                        getsd = FALSE,
                        verbose = 0))
  )


  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # MSE ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  mse <- Rceattle::run_mse(om = ss_run, em = ss_run_Tier3, nsim = 1, assessment_period = 1, sampling_period = 1, simulate_data = FALSE, sample_rec = FALSE)

  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Tests ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  nyrs <- length(1979:2040)
  testthat::expect_equal(24, length(mse$Sim_1$EM))
  testthat::expect_equal(rep(0.4, 3), as.numeric(mse$Sim_1$OM$quantities$ssb_depletion[,nyrs]), tolerance = 0.005)

  # OM_no_F is the OM with no fishing after the original OM's terminal year: the
  # two share the hindcast, and the no-F projection is unfished. remove_F() on
  # the advanced OM with its default start would begin after 2040 and remove
  # nothing.
  om_yrs <- mse$Sim_1$OM$data_list$styr:mse$Sim_1$OM$data_list$projyr
  pre    <- which(om_yrs <= ss_run$data_list$endyr)
  post   <- which(om_yrs >  ss_run$data_list$endyr)
  testthat::expect_equal(unname(mse$Sim_1$OM_no_F$quantities$ssb[, pre]),
                         unname(mse$Sim_1$OM$quantities$ssb[, pre]), tolerance = 1e-8)
  testthat::expect_true(all(mse$Sim_1$OM_no_F$quantities$F_spp[, post] < 1e-8))
  testthat::expect_gt(max(mse$Sim_1$OM$quantities$F_spp[, post]), 1e-3)

  # mse_summary() derives the performance-metric table from the MSE run. It is
  # exported but was previously untested; guard that it runs end-to-end and
  # returns the expected per-species metrics (it does substantial index math and
  # references the package hcr_map constant).
  summ <- Rceattle::mse_summary(mse)
  testthat::expect_named(summ, c("species", "fleet", "total", "meta"))

  # Per-species conservation metrics: one row per species, no NA padding.
  testthat::expect_s3_class(summ$species, "data.frame")
  testthat::expect_equal(nrow(summ$species), mse[[1]]$OM$data_list$nspp)
  testthat::expect_true(all(c("om_avg_depletion", "om_sims_collapsed") %in%
                              names(summ$species)))
  dep <- summ$species[["om_avg_depletion"]]
  testthat::expect_true(all(is.finite(dep)))
  testthat::expect_true(all(dep > 0 & dep < 1.5))   # near the Tier-3 40% target

  # Per-fleet catch metrics: one row per fishery fleet, integer-keyed, finite.
  testthat::expect_s3_class(summ$fleet, "data.frame")
  testthat::expect_true("avg_catch" %in% names(summ$fleet))
  n_fishery <- sum(mse[[1]]$OM$data_list$fleet_control$Fleet_type == "Fishery")
  testthat::expect_equal(nrow(summ$fleet), n_fishery)
  testthat::expect_type(summ$fleet$Fleet_code, "integer")
  testthat::expect_true(all(is.finite(summ$fleet[["avg_catch"]])))

  # Across-fleet totals.
  testthat::expect_true(is.finite(summ$total[["avg_catch"]]))

  # -- A species with no fishery --------------------------------------------
  # A multispecies model can carry a predator for its consumption alone, with no
  # fishery and so no `catch_data` rows (arrowtooth and sablefish do this in the
  # hake model). Every per-species catch metric then works on a zero-length
  # vector: `Average Catch` warned "argument is not numeric or logical" and
  # returned NA, while `Catch IAV` and `P(Closed)` divided by a length of zero
  # and returned NaN with no warning at all.
  #
  # Doctored from the fitted `mse` rather than refit: species 3's fishery is
  # removed from `fleet_control` and its catch rows dropped, which is exactly
  # the shape an unfished predator arrives in.
  mse_unfished <- mse
  for (i in seq_along(mse_unfished)) {
    fc <- mse_unfished[[i]]$OM$data_list$fleet_control
    drop_flt <- fc$Fleet_code[fc$Species == 3 & fc$Fleet_type == "Fishery"]
    mse_unfished[[i]]$OM$data_list$fleet_control <-
      fc[!(fc$Fleet_code %in% drop_flt), , drop = FALSE]
    cd <- mse_unfished[[i]]$OM$data_list$catch_data
    mse_unfished[[i]]$OM$data_list$catch_data <-
      cd[cd$Species != 3, , drop = FALSE]
  }

  summ_u <- testthat::expect_no_warning(Rceattle::mse_summary(mse_unfished))

  # Zero catch is the real answer for a species nobody fishes.
  testthat::expect_equal(summ_u$species[["avg_catch"]][3], 0)
  # Catch variability and P(fishery closed) are undefined without a fishery --
  # NA, and specifically not the NaN the ratios used to produce.
  testthat::expect_true(is.na(summ_u$species[["catch_iav"]][3]))
  testthat::expect_true(is.na(summ_u$species[["p_closed"]][3]))
  testthat::expect_false(is.nan(summ_u$species[["catch_iav"]][3]))
  testthat::expect_false(is.nan(summ_u$species[["p_closed"]][3]))

  # The fished species must be untouched by that branch.
  for (cc in c("avg_catch", "catch_iav", "p_closed")) {
    testthat::expect_equal(summ_u$species[[cc]][1:2], summ$species[[cc]][1:2])
  }
  # ...and the unfished species drops out of the per-fleet frame entirely.
  testthat::expect_false(3 %in% summ_u$fleet$Fleet_code[
    summ_u$fleet$Fleet_name %in% "ATF"])
})

testthat::test_that("Test MSE - Tier 3 parallel", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")
  testthat::skip_on_cran()


  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Data ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  data(BS2017SS) # ?BS2017SS for more information on the data
  BS2017SS$projyr <- 2027
  BS2017SS$fleet_control$Proj_F_proportion <-rep(1,7)


  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Operating Model ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Single-species with fixed M
  ss_run <- Rceattle::fit_mod(data_list = BS2017SS,
                              inits = NULL, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 1, # Estimate hindcast only
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              fit_control = fit_control(
                                phase = TRUE,
                                getsd = FALSE,
                                verbose = 0))

  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # NPFMC Tier 3 ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  ss_run_Tier3 <- suppressWarnings(
    Rceattle::fit_mod(data_list = BS2017SS,
                      inits = ss_run$estimated_params, # Initial parameters from ss_run
                      estimateMode = 2, # Run projection only
                      HCR = build_hcr(HCR = 5, # Tier3 HCR
                                      Ftarget = 0.4, # F40%
                                      Flimit = 0.35, # F35%
                                      Plimit = 0.2, # No fishing when SB<SB20
                                      Alpha = 0.05),
                      msmMode = 0, # Single species mode
                      fit_control = fit_control(
                        getsd = FALSE,
                        verbose = 0))
  )


  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # MSE ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  mse <- Rceattle::run_mse(om = ss_run, em = ss_run_Tier3, nsim = 2, assessment_period = 1, sampling_period = 1, simulate_data = FALSE, sample_rec = FALSE)

  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Tests ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  testthat::expect_equal(2, length(mse))
  # length() alone holds whether the simulations ran or came back as failure
  # markers, and a failure no longer propagates out of the dispatch -- so assert
  # they actually completed. This is the only nsim > 1 case, and therefore the
  # only one that exercises the parallel path.
  testthat::expect_true(all(vapply(mse, function(x) isTRUE(x$use_sim), logical(1))))
  testthat::expect_true(all(vapply(mse, function(x) !is.null(x$OM), logical(1))))
})

# A survey fleet legitimately carries NA Proj_F_proportion -- it takes no catch,
# so there is no F to apportion, and read_data()/build_data() write NA there.
# The guard that some fleet takes projected F summed the column without na.rm,
# so a single NA made the sum NA and the `if` errored with "missing value where
# TRUE/FALSE needed" before the MSE started.
testthat::test_that("run_mse tolerates NA Proj_F_proportion on non-catch fleets", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_on_cran()

  data(BS2017SS)
  BS2017SS$projyr <- 2027
  fc <- BS2017SS$fleet_control
  BS2017SS$fleet_control$Proj_F_proportion <-
    ifelse(fc$Fleet_type %in% c(1, "Fishery"), 1, NA_real_)
  testthat::expect_true(anyNA(BS2017SS$fleet_control$Proj_F_proportion))

  ss <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = BS2017SS, inits = NULL, file = NULL, estimateMode = 1,
    random_rec = FALSE, msmMode = 0,
    fit_control = Rceattle::fit_control(phase = TRUE, getsd = FALSE, verbose = 0))))
  om <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = BS2017SS, inits = ss$estimated_params, file = NULL,
    estimateMode = 0, random_rec = FALSE, msmMode = 0,
    HCR = Rceattle::build_hcr(HCR = 5, Ftarget = 0.4, Flimit = 0.35,
                              Plimit = 0.2, Alpha = 0.05),
    fit_control = Rceattle::fit_control(phase = TRUE, getsd = FALSE, verbose = 0))))

  # The guard is reached before any projection work, so a 1-simulation run is
  # enough to show it no longer errors on the NA.
  testthat::expect_no_error(
    suppressMessages(suppressWarnings(
      Rceattle::run_mse(om = om, em = om, nsim = 1, assessment_period = 5,
                        sampling_period = 5, simulate_data = FALSE,
                        sample_rec = FALSE))))
})

# All-zero really is an error: nothing would take the projected catch.
testthat::test_that("run_mse still errors when no fleet takes projected F", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_on_cran()

  data(BS2017SS)
  BS2017SS$projyr <- 2027
  BS2017SS$fleet_control$Proj_F_proportion <- rep(0, 7)

  om <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = BS2017SS, inits = NULL, file = NULL, estimateMode = 1,
    random_rec = FALSE, msmMode = 0,
    fit_control = Rceattle::fit_control(phase = TRUE, getsd = FALSE, verbose = 0))))

  testthat::expect_error(
    suppressMessages(suppressWarnings(
      Rceattle::run_mse(om = om, em = om, nsim = 1, assessment_period = 5,
                        sampling_period = 5))),
    "Proj_F_proportion")
})
