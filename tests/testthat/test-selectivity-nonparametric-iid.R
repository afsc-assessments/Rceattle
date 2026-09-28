# Time-varying non-parametric (Ianelli, type 2) selectivity.
#
# Provenance. `Time_varying_sel = "IID"` on a non-parametric fleet was refused by
# build_map() and data_check() even though the bundled `Atka2022` -- and the ADMB
# bridge it comes from (`Rceattle-models/BSAI atka mackerel/2022 atka
# bridging.R`) -- is configured exactly that way, so that dataset could not be
# fitted as shipped. 5.25.0 scores the mode instead of refusing it:
# `dnorm(sel_coff_dev, 0, sel_dev_sd)` on each estimated coefficient.
#
# Why the two modes are not interchangeable, and why only one may be integrated:
# "RandomWalk" scores the year-to-year change in the REALIZED log-selectivity,
# which is renormalized to mean 1 within each year. The level of a year's
# coefficients does not enter that density, so `sel_dev_sd` is not identified
# from it and Laplace-integrating the deviates drives it to zero -- measured at
# 2.7e-8 on Atka2022, which is a time-invariant selectivity reported as a
# time-varying one. "IID" scores the deviations themselves and does identify the
# sd. But the AMAK shape penalties and the average-selectivity term are charged
# on each year's realized curve and do not scale with the sd, so the density
# integrated is tilted and the reported sd is biased low: fit_mod() refuses
# `random_sel = TRUE` for the walk because the density is improper, and for IID
# because the sd it would report is not the sd of the deviations.
testthat::skip_on_cran()

testthat::test_that("NonParametric accepts IID and the bundled Atka2022 fits", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("Atka2022")
  # As shipped: fishery is Selectivity 2 with Time_varying_sel 1 (IID), sd 0.35.
  testthat::expect_equal(as.numeric(Atka2022$fleet_control$Selectivity[2]), 2)
  testthat::expect_equal(as.numeric(Atka2022$fleet_control$Time_varying_sel[2]), 1)

  fit <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = Atka2022, inits = NULL, msmMode = 0, estimateMode = "Hindcast",
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0))))
  testthat::expect_true(is.finite(fit$opt$objective))

  # The deviates are estimated and the curve actually moves between years.
  sel <- fit$quantities$sel_at_age[2, 1, , ]
  testthat::expect_gt(max(apply(sel, 1, function(x) diff(range(x)))), 0)
})


testthat::test_that("random_sel = TRUE is refused for both non-parametric modes", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("Atka2022")
  d <- Atka2022
  d$fleet_control$Time_varying_sel <- as.character(d$fleet_control$Time_varying_sel)

  # The walk cannot be integrated: its density is blind to the level of each
  # year's coefficients, so the sd collapses. Refused by name, not by a bare TMB
  # "NA/NaN gradient evaluation" that mentions neither selectivity nor random_sel.
  d$fleet_control$Time_varying_sel[2] <- "RandomWalk"
  testthat::expect_error(
    suppressMessages(suppressWarnings(Rceattle::fit_mod(
      data_list = d, inits = NULL, msmMode = 0, estimateMode = "Hindcast",
      random_sel = TRUE,
      fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0)))),
    "random_sel"
  )

  # ... but it still fits as the penalized effects the AMAK formulation intends.
  testthat::expect_true(is.finite(suppressMessages(suppressWarnings(
    Rceattle::fit_mod(data_list = d, inits = NULL, msmMode = 0,
                      estimateMode = "Hindcast",
                      fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0))
  ))$opt$objective))

  # IID is refused too, whatever the shape penalties are set to: the
  # average-selectivity penalty is always charged on the realized curve, so the
  # reported deviation sd is the sd of a tilted density. Until 5.35.0 the
  # refusal was lifted at Sel_curve_pen1 = 0, which removed the kink (Atka2022
  # stopped at a maximum gradient of 6.8 with the penalty on) but not the bias.
  d$fleet_control$Time_varying_sel[2] <- "IID"
  for (pen1 in list(d$fleet_control$Sel_curve_pen1[2], 0)) {
    d$fleet_control$Sel_curve_pen1[2] <- pen1
    testthat::expect_error(
      suppressMessages(suppressWarnings(Rceattle::fit_mod(
        data_list = d, inits = NULL, msmMode = 0, estimateMode = "Hindcast",
        random_sel = TRUE,
        fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0)))),
      "average-selectivity"
    )
  }
})


testthat::test_that("NonParametricPM still refuses IID, naming the alternative", {
  testthat::skip_if_not_installed("Rceattle")

  data("Atka2022")
  d <- Atka2022
  d$fleet_control$Selectivity <- as.character(d$fleet_control$Selectivity)
  d$fleet_control$Selectivity[2] <- "NonParametricPM"
  # Its sel_coff_dev ARE walk increments (selectivity.hpp case 9 builds each
  # year from the previous one), so an independent-deviate reading of them would
  # describe a different curve than the model draws.
  testthat::expect_error(
    suppressMessages(suppressWarnings(Rceattle::fit_mod(
      data_list = d, inits = NULL, msmMode = 0, estimateMode = 3,
      fit_control = Rceattle::fit_control(verbose = 0)))),
    "NonParametricPM"
  )
})


# The random_sel guard puts NonParametricPM in its walk branch
# (.rce_np_unintegrable_fleets, R/3-build_map.R), but only the IID arm was
# exercised, so the form-9 walk arm could have stopped refusing unnoticed.
testthat::test_that("NonParametricPM under RandomWalk refuses random_sel, naming the integrable form", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("Rceattle")

  data("Atka2022")
  d <- Atka2022
  d$fleet_control$Selectivity <- as.character(d$fleet_control$Selectivity)
  d$fleet_control$Time_varying_sel <- as.character(d$fleet_control$Time_varying_sel)
  d$fleet_control$Selectivity[2]       <- "NonParametricPM"
  d$fleet_control$Time_varying_sel[2]  <- "RandomWalk"
  testthat::expect_error(
    suppressMessages(suppressWarnings(Rceattle::fit_mod(
      data_list = d, inits = NULL, msmMode = 0, estimateMode = 3,
      random_sel = TRUE,
      fit_control = Rceattle::fit_control(verbose = 0)))),
    "NonParametricIntegrable"
  )
})
