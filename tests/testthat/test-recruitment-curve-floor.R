# An identity-link recruitment offset can drive the curve to or below zero.
# Measured on make_test_data() before the floor (5.34.0), Ianelli form with an
# identity alpha offset of -100 per unit covariate: R_hat reached -89.2, and the
# objective, the stock-recruit penalty and dynamic B0 were all NaN. Recruitment,
# the curve and R_hat are now floored at one fish (1e-3 thousand) through
# posfun(), with the excursion charged to the zero-N penalty row.

testthat::skip_on_cran()

curve_floor_fixture <- function() {
  set.seed(7)
  d <- make_test_data()
  if (!is.null(d$data_list)) d <- d$data_list
  yrs <- d$styr:d$projyr
  d$env_data <- data.frame(Year = yrs, EnvData = round(stats::rnorm(length(yrs)), 3))
  d
}

curve_floor_build <- function(d, recFun) {
  fit_mod(data_list = d, inits = NULL, estimateMode = 3, msmMode = 0,
          random_rec = FALSE, recFun = recFun,
          fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE))
}

testthat::test_that("a non-positive curve is floored with a penalty, and the offset is warned about", {
  d <- curve_floor_fixture()
  rec <- build_srr(srr_fun = "mean", srr_pred_fun = "BevertonHolt",
                   linkages = list(alpha = linkage_spec(~ 1 + EnvData, link = "identity",
                                                        init = list(EnvData = -100),
                                                        data = d$env_data)))
  testthat::expect_warning(
    m <- suppressMessages(curve_floor_build(d, rec)),
    "identity-link recruitment linkage")
  testthat::expect_true(is.finite(m$obj$fn()))
  testthat::expect_true(all(is.finite(m$obj$gr())))
  testthat::expect_gt(m$quantities$jnll_comp["Zero n-at-age penalty", 1], 0)
  testthat::expect_true(all(is.finite(m$quantities$DynamicSB0)))
  # Known gap (inst/dev/TODO-srr-multispecies.md item 14): the recruitment behind
  # dynamic B0 is not floored, so it goes negative here. This records the gap and
  # flips when it is closed.
  testthat::expect_true(any(m$quantities$DynamicSB0 < 0))
})

testthat::test_that("the Ricker log arguments are kept positive too", {
  d <- curve_floor_fixture()
  rec <- build_srr(srr_fun = "mean", srr_pred_fun = "Ricker",
                   linkages = list(alpha = linkage_spec(~ 1 + EnvData, link = "identity",
                                                        init = list(EnvData = -100),
                                                        data = d$env_data)))
  m <- suppressWarnings(suppressMessages(curve_floor_build(d, rec)))
  testthat::expect_true(is.finite(m$obj$fn()))
  testthat::expect_true(all(is.finite(m$obj$gr())))
  testthat::expect_true(all(is.finite(m$quantities$R_hat)))
  testthat::expect_true(all(is.finite(m$quantities$steepness)))
  # A non-finite penalty row is a FAIL, not silence.
  fake <- m; fake$quantities$jnll_comp["Zero n-at-age penalty", 1] <- NaN
  testthat::expect_equal(Rceattle:::.check_zero_n_penalty(fake)$zero_n_penalty$severity, "FAIL")
})

testthat::test_that("a curve fitted in the hindcast is floored too, and check_convergence() reports it", {
  d <- curve_floor_fixture()
  rec <- build_srr(srr_fun = "BevertonHolt",
                   linkages = list(alpha = linkage_spec(~ 1 + EnvData, link = "identity",
                                                        init = list(EnvData = -100),
                                                        data = d$env_data)))
  m <- suppressWarnings(suppressMessages(curve_floor_build(d, rec)))
  testthat::expect_true(is.finite(m$obj$fn()))
  testthat::expect_true(all(m$quantities$R > 0))
  testthat::expect_true(all(m$quantities$R_hat > 0))
  testthat::expect_equal(Rceattle:::.check_zero_n_penalty(m)$zero_n_penalty$severity, "FAIL")
})

testthat::test_that("the floor is inert where the curve is positive", {
  d <- curve_floor_fixture()
  rec <- build_srr(srr_fun = "mean", srr_pred_fun = "BevertonHolt",
                   linkages = list(alpha = linkage_spec(~ 1 + EnvData, link = "identity",
                                                        init = list(EnvData = 0),
                                                        data = d$env_data)))
  m <- suppressWarnings(suppressMessages(curve_floor_build(d, rec)))
  testthat::expect_equal(unname(m$quantities$jnll_comp["Zero n-at-age penalty", 1]), 0)
})
