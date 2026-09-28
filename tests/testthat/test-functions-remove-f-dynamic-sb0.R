# remove_F() refits with F = 0 from `styr`, by default the year after endyr, so
# the hindcast is unchanged. Up to 5.33.0 it started the year after the latest
# suit_endyr, which removed fishing inside the hindcast whenever the suitability
# window ended early -- on the Pacific hake MSE, 2020-2023 of a 2023 hindcast.
# A start inside the suitability window is refused under predation: empirical
# suitability is computed from the fitted abundance, so it would change.
# mse_summary() reads dynamic depletion from the OM's DynamicSB0 for single- and
# multispecies models alike.

.msm_early_suit_fit <- function() {
  set.seed(123)
  d <- make_msm_test_data()$data_list
  suppressMessages(suppressWarnings(fit_mod(
    data_list = d, inits = NULL, estimateMode = 3, msmMode = 1, suitMode = 0,
    initMode = "NonEquilibrium", random_rec = FALSE,
    suit_styr = d$styr, suit_endyr = d$endyr - 5,
    fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE))))
}

test_that("remove_F() keeps the hindcast and removes fishing from styr", {
  testthat::skip_on_cran()
  fit <- .msm_early_suit_fit()
  dl  <- fit$data_list
  nh  <- length(dl$styr:dl$endyr)

  no_f <- suppressMessages(suppressWarnings(remove_F(fit)))
  expect_equal(unname(no_f$quantities$ssb[, 1:nh]), unname(fit$quantities$ssb[, 1:nh]),
               tolerance = 1e-10)
  expect_equal(unname(no_f$quantities$F_spp[, 1:nh]), unname(fit$quantities$F_spp[, 1:nh]),
               tolerance = 1e-10)

  # Starting in the hindcast, after the suitability window, removes that fishing.
  sy    <- max(dl$suit_endyr) + 1
  early <- suppressMessages(suppressWarnings(remove_F(fit, styr = sy)))
  cols  <- (sy - dl$styr + 1):nh
  expect_true(all(early$quantities$F_spp[, cols] < 1e-10))
  expect_gt(max(fit$quantities$F_spp[, cols]), 1e-3)

  expect_error(remove_F(fit, styr = max(dl$suit_endyr)), "suitability window")
  expect_error(remove_F(fit, styr = dl$styr), "styr")
  # The projection is always unfished, so a later start would be ignored.
  expect_error(remove_F(fit, styr = dl$endyr + 2), "year after endyr")
})

test_that("remove_F() leaves the projection unfished under a harvest control rule", {
  testthat::skip_on_cran()
  # Fitted (estimateMode 0) so the rule is evaluated and the projection is fished.
  d   <- make_test_data()
  fit <- suppressMessages(suppressWarnings(fit_mod(
    data_list = d, inits = NULL, estimateMode = 0, random_rec = FALSE,
    HCR = build_hcr(HCR = "ConstantF", Ftarget = 0.2),
    fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE))))
  proj <- (d$endyr + 1):d$projyr - d$styr + 1
  expect_gt(max(fit$quantities$F_spp[, proj]), 1e-3)   # the rule fishes the projection

  no_f <- suppressMessages(suppressWarnings(remove_F(fit)))
  expect_true(all(no_f$quantities$F_spp[, proj] < 1e-10))
})

test_that("parametric suitability does not restrict styr", {
  testthat::skip_on_cran()
  set.seed(123)
  d   <- make_msm_test_data()$data_list
  fit <- suppressMessages(suppressWarnings(fit_mod(
    data_list = d, inits = NULL, estimateMode = 3, msmMode = 1, suitMode = 4,
    initMode = "NonEquilibrium", random_rec = FALSE,
    suit_styr = d$styr, suit_endyr = d$endyr - 5,
    fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE))))
  # Suitability from predator-prey weight ratios does not read abundance, so a
  # start inside its window leaves it unchanged and is allowed.
  expect_no_error(suppressMessages(suppressWarnings(remove_F(fit, styr = d$endyr - 5))))
})

test_that("a prey-only species' suitability window does not restrict styr", {
  testthat::skip_on_cran()
  set.seed(123)
  d <- make_msm_test_data()$data_list
  d$diet_data <- d$diet_data[d$diet_data$Pred != 2, ]   # species 2 eats nothing
  fit <- suppressMessages(suppressWarnings(fit_mod(
    data_list = d, inits = NULL, estimateMode = 3, msmMode = 1, suitMode = 0,
    initMode = "NonEquilibrium", random_rec = FALSE,
    suit_styr = d$styr, suit_endyr = c(d$endyr - 5, d$endyr),
    fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE))))
  s <- fit$quantities$suitability                        # [pred, prey, pred age, prey age, yr]
  expect_gt(max(abs(s[1, , , , ])), 0)                   # species 1 eats, so has suitability
  expect_true(all(s[2, , , , ] == 0))                    # species 2 eats nothing
  # Only species 1's window (to endyr - 5) binds; suitability is unchanged after it.
  no_f <- suppressMessages(suppressWarnings(remove_F(fit, styr = d$endyr - 2)))
  expect_equal(no_f$quantities$suitability, fit$quantities$suitability, tolerance = 1e-12)
  expect_error(remove_F(fit, styr = d$endyr - 5), "suitability window")
})

test_that("mse_summary() takes multispecies dynamic depletion from DynamicSB0", {
  testthat::skip_on_cran()
  fit <- .msm_early_suit_fit()
  dl  <- fit$data_list
  mse <- list(Sim_1 = list(EM = list(EM = fit), use_sim = TRUE, failure = NA,
                           OM = fit, OM_no_F = suppressMessages(suppressWarnings(remove_F(fit)))))
  summ <- suppressWarnings(mse_summary(mse, om_only = TRUE))
  col  <- dl$projyr - dl$styr + 1
  expect_equal(unname(summ$species$om_terminal_dynamic_sb0),
               unname(fit$quantities$DynamicSB0[, col]), tolerance = 1e-10)
  expect_equal(unname(summ$species$om_terminal_depletion_dynamic),
               unname(fit$quantities$ssb[, col] / fit$quantities$DynamicSB0[, col]), tolerance = 1e-10)
})
