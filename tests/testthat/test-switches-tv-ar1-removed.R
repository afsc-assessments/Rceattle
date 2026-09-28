# Time_varying_sel / Time_varying_q = "AR1" (value 2) are removed in 5.16.0.
#
# Provenance: the value was accepted by validate_switches() and then scored by
# the template with the SAME independent normal penalty as "IID"
# (`flt_varying_sel == 1 || == 2`, `index_varying_q == 1 || == 2`). No
# correlation parameter was ever read -- the one that existed belonged to the
# Catchability = "AR1" (QAR1) path, removed in 5.37.0, and selectivity has no equivalent. So a
# fleet set to "AR1" fitted independent deviations under a name that promised
# correlated ones, silently. These tests pin the refusal, the redirect, and the
# fact that "IID" is the exact fit the removed value gave.

testthat::skip_on_cran()

tv_dat <- function(col, value, fleet = 2L) {
  data("GOA2018SS", envir = environment())
  d <- suppressWarnings(suppressMessages(Rceattle::switch_check(GOA2018SS)))
  d$fleet_control[[col]][fleet] <- value
  d
}

run_check <- function(d) {
  suppressWarnings(suppressMessages(Rceattle:::data_check(d)))
}


testthat::test_that("Time_varying_q = 'AR1' is refused, by name and by code", {
  for (v in list("AR1", 2, "2")) {
    testthat::expect_error(
      run_check(tv_dat("Time_varying_q", v)),
      regexp = "Time_varying_q = 'AR1' is removed")
  }
})


testthat::test_that("Time_varying_sel = 'AR1' is refused, by name and by code", {
  for (v in list("AR1", 2, "2")) {
    testthat::expect_error(
      run_check(tv_dat("Time_varying_sel", v)),
      regexp = "Time_varying_sel = 'AR1' is removed")
  }
})


testthat::test_that("the refusal names the fleet and both replacements", {
  msg <- tryCatch(run_check(tv_dat("Time_varying_q", "AR1")),
                  error = function(e) conditionMessage(e))
  # The fleet, so an assessor knows which row to change...
  testthat::expect_match(msg, "Pollock_survey_2_bottom_trawl", fixed = TRUE)
  # ...the switch value that reproduces the fit they actually had...
  testthat::expect_match(msg, "Set Time_varying_q = 'IID'", fixed = TRUE)
  # ...and the grammar that fits the AR1 the name promised.
  testthat::expect_match(msg, "ar1(1 | Year)", fixed = TRUE)

  msg_sel <- tryCatch(run_check(tv_dat("Time_varying_sel", "AR1")),
                      error = function(e) conditionMessage(e))
  testthat::expect_match(msg_sel, "build_selectivity", fixed = TRUE)
  testthat::expect_match(msg_sel, "ar1(1 | Year)", fixed = TRUE)
})


testthat::test_that("an env-column index is not mistaken for the AR1 mode", {
  # Time_varying_q is OVERLOADED: under Catchability = "Environmental" it holds
  # a 1-based env_data column index, not a mode. A fleet naming env column 2
  # therefore carries a literal 2, which canonicalizes to "AR1". Refusing it
  # would reject a working environmental-catchability model for a setting the
  # assessor never made. validate_switches() and the soft-deprecation carry the
  # same exemption.
  for (idx in c(1, 2, 3)) {
    d <- tv_dat("Time_varying_q", idx)
    d$fleet_control$Catchability[2] <- "Environmental"
    testthat::expect_no_error(run_check(d))
  }

  # ...but an ESTIMATED-q fleet reads the same column as a mode, so 2 is refused.
  d <- tv_dat("Time_varying_q", 2)
  d$fleet_control$Catchability[2] <- "Estimated"
  testthat::expect_error(run_check(d), regexp = "Time_varying_q = 'AR1' is removed")
})


testthat::test_that("the env-column exemption fails towards refusing", {
  # The exemption is read off Catchability. `.canon_switch()` returns length 0
  # for an absent column, and `hit & !logical(0)` is logical(0) -- so an
  # unguarded exemption would make the refusal vanish whenever a DIFFERENT
  # column is missing. A removed value must not become acceptable that way.
  d <- tv_dat("Time_varying_q", "AR1")
  d$fleet_control$Catchability <- rep(NA, nrow(d$fleet_control))
  # Something must still complain; what matters is that nothing passes silently.
  testthat::expect_error(run_check(d))
})


testthat::test_that("only 'AR1' is refused -- the other modes still pass", {
  for (v in c("Off", "IID", "Block", "RandomWalk")) {
    testthat::expect_no_error(run_check(tv_dat("Time_varying_q", v)))
  }
  for (v in c("Off", "Block")) {
    testthat::expect_no_error(run_check(tv_dat("Time_varying_sel", v)))
  }
})


testthat::test_that("'IID' reproduces the removed value's fit, bit for bit", {
  testthat::skip_if_not_installed("TMB")

  # The error tells the assessor "Set Time_varying_* = 'IID' to keep the fit you
  # had". That is a promise about the TEMPLATE, so test it there rather than
  # inferring it from the `|| == 2` branch: build the model under "IID", flip
  # the two switch codes to 2 in the TMB data -- which is exactly what the
  # removed value produced -- and rebuild on the same parameters and map.
  d <- suppressWarnings(suppressMessages(Rceattle::switch_check(GOA2018SS)))
  d$fleet_control$Time_varying_sel[8]    <- "IID"   # DoubleLogistic fishery
  d$fleet_control$Time_varying_sel_sd[8] <- 0.05
  d$fleet_control$Time_varying_q[2]      <- "IID"   # estimated-q survey
  d$fleet_control$Time_varying_q_sd[2]   <- 0.1

  iid <- suppressWarnings(suppressMessages(
    Rceattle::fit_mod(data_list = d, inits = NULL, file = NULL,
                      estimateMode = 3, random_rec = FALSE, msmMode = 0,
                      fit_control = Rceattle::fit_control(verbose = 0))))

  dat <- iid$obj$env$data
  testthat::expect_equal(dat$flt_varying_sel[8], 1L)
  testthat::expect_equal(dat$index_varying_q[2], 1L)
  dat$flt_varying_sel[8] <- 2L
  dat$index_varying_q[2] <- 2L
  attr(dat, "check.passed") <- NULL

  as_ar1 <- TMB::MakeADFun(
    dat, parameters = iid$obj$env$parList(iid$obj$env$last.par.best),
    map = iid$obj$env$map, DLL = iid$TMBfilename, silent = TRUE)

  p <- iid$obj$env$last.par.best[iid$obj$env$lfixed()]
  testthat::expect_identical(iid$obj$fn(p), as_ar1$fn(p))
  testthat::expect_equal(max(abs(iid$obj$gr(p) - as_ar1$gr(p))), 0)
  testthat::expect_equal(
    max(abs(iid$obj$report(p)$jnll_comp - as_ar1$report(p)$jnll_comp)), 0)
})
