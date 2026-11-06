# =============================================================================
# Time_varying_sel compatibility warning for "Hake" (Taylor et al. 2014)
# selectivity, in build_map_selectivity().
#
# The Hake form supports exactly two time-varying modes -- "Off" (0), i.e. the
# bin coefficients with no annual deviates, and "IID" (1), penalized coefficient
# deviates -- and data_check() enforces that set. The map builder's guard used to
# be `tv_sel != "IID"`, which fired on "Off" as well and so warned about a
# configuration its own message (and the validator) call legal.
#
# switch_check() canonicalizes the integer 0 to the string "Off", so a fleet
# entered as Time_varying_sel = 0 reaches the guard as "Off" -- which is what
# made the spurious warning read "Current value: Off". fit_mod() wraps build_map()
# in suppressWarnings(), so it surfaced only in the callers that build the map
# directly (run_mse(), retrospective()).
# =============================================================================

.hake_data <- function(tv_sel) {
  dat <- make_test_data(nyrs = 6, nprojyrs = 2, nages = 5)
  # fit_mod() normally supplies growth_model from growthFun(); build_map() is
  # called directly here, so set the empirical-growth default explicitly.
  dat$growth_model <- rep(0, dat$nspp)
  dat$fleet_control$Selectivity <- "Hake"
  dat$fleet_control$N_sel_bins <- 4
  dat$fleet_control$Age_first_selected <- 1
  dat$fleet_control$Time_varying_sel <- tv_sel
  dat
}

# Warnings raised while building the map, with build_params() noise excluded.
.map_warnings <- function(dat) {
  params <- suppressMessages(suppressWarnings(Rceattle::build_params(dat)))
  testthat::capture_warnings(
    suppressMessages(Rceattle::build_map(dat, params))
  )
}

testthat::test_that("Hake selectivity does not warn for Time_varying_sel 'Off' or 'IID'", {
  # The user-facing case: one fleet with deviates, one without (0 -> "Off").
  w <- .map_warnings(.hake_data(c(1, 0)))
  testthat::expect_length(grep("Time_varying_sel", w), 0)

  # ... and each mode on its own, entered as both the integer and the string.
  for (tv in list(0, 1, "Off", "IID")) {
    w <- .map_warnings(.hake_data(tv))
    testthat::expect_length(grep("Time_varying_sel", w), 0)
  }
})

testthat::test_that("Hake selectivity still warns for a mode it cannot represent", {
  w <- .map_warnings(.hake_data(c("IID", "Block")))
  hits <- grep("Time_varying_sel", w, value = TRUE)
  testthat::expect_length(hits, 1)
  # Names the offending fleet, and the canonical "Off" -- not "None", which is
  # not a Time_varying_sel value.
  testthat::expect_match(hits, "fleet 2")
  testthat::expect_match(hits, "Current value: Block")
  testthat::expect_match(hits, "'Off' \\(0\\)")
})

testthat::test_that("the warned set is exactly the complement of data_check()'s", {
  # data_check() is the authority on which modes the Hake form accepts
  # (R/1-data_check.R: "must be 'Off' or 'IID'"). The map guard must not
  # contradict it in either direction, across every Time_varying_sel mode.
  allowed <- c("Off", "IID")
  for (tv in names(Rceattle:::tv_sel_map)) {
    n <- length(grep("Time_varying_sel", .map_warnings(.hake_data(tv))))
    testthat::expect_equal(n, if (tv %in% allowed) 0L else 2L, info = tv)
  }
})
