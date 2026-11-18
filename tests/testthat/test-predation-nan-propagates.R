# TMB evaluates an `if` on a parameter-dependent value once, when it records the
# tape, and every later objective evaluation keeps the branch taken there. The
# available-food and prey-biomass tests in predation.hpp read `> 0` until 5.67.11,
# so NaN failed them and that predator or prey age dropped out of the tape.
# BS2017MS (msmMode = 1, niter = 3) recorded at pollock init_dev = -720 (lower
# bound -1000) returned 7,372,253 there instead of NaN, and 1,564,196.9 instead of
# 1,784,463.7 back at the default start. Fixing only the available-food test still
# returned 7,128,177 and 1,394,720.3: the next iteration's prey-biomass test
# dropped the NaN prey instead.

nan_msm_build <- function(inits = NULL) {
  data("BS2017MS", package = "Rceattle", envir = environment())
  suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = BS2017MS, inits = inits, file = NULL, estimateMode = 3,
    niter = 3, msmMode = 1, suitMode = rep(0, 3), random_rec = FALSE,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0))))
}

testthat::test_that("NaN available food reaches the objective of a tape recorded there", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  start <- nan_msm_build()
  q <- start$estimated_params
  q$init_dev[1, seq_len(start$data_list$nages[1] - 1)] <- -720  # pollock, ages 2+
  bad <- nan_msm_build(q)
  xq <- bad$obj$par

  # Positive control: this input does make available food NaN.
  testthat::expect_gt(sum(is.nan(bad$obj$report(xq)$avail_food)), 0)
  testthat::expect_true(is.nan(bad$obj$fn(xq)))

  # Recorded there, the tape is still the full model back at the default start.
  x0 <- start$obj$par
  testthat::expect_identical(names(bad$obj$par), names(x0))
  testthat::expect_equal(bad$obj$fn(x0), start$obj$fn(x0), tolerance = 1e-12)
  testthat::expect_equal(bad$obj$gr(x0), start$obj$gr(x0), tolerance = 1e-12)
})
