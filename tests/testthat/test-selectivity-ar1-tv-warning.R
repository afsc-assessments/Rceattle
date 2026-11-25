# =============================================================================
# Time_varying_sel compatibility warning for the 2DAR1 / 3DAR1 selectivity
# forms, in build_map_selectivity().
#
# The AR1 forms carry their own annual deviations -- the field supplies them,
# with the correlations estimated from Sel_curve_pen1/2/3 -- so a
# Time_varying_sel setting has nothing left to do, and the builder says which
# fleet it overrode. "Off" asks for no extra deviation layer, so it agrees with
# the AR1 field rather than conflicting with it, and is silent.
#
# Provenance: a live GOA pollock 2DAR1 run setting Time_varying_sel = "Off"
# warned that its choice was ignored. The same class on the Hake form is
# test-selectivity-hake-tv-warning.R, which is why both AR1 forms are pinned
# here rather than only the one reported.
#
# NA is not covered: build_map_selectivity() reaches `if (tv_sel == "Block")`
# before any AR1 branch and errors on a missing value, for every selectivity
# form. That is a separate gap.
# =============================================================================

.ar1_data <- function(tv_sel, sel = "2DAR1") {
  dat <- make_test_data(nyrs = 6, nprojyrs = 2, nages = 5)
  # fit_mod() normally supplies growth_model from growthFun(); build_map() is
  # called directly here, so set the empirical-growth default explicitly.
  dat$growth_model <- rep(0, dat$nspp)
  dat$fleet_control$Selectivity <- sel
  dat$fleet_control$N_sel_bins <- 4
  dat$fleet_control$Age_first_selected <- 1
  dat$fleet_control$Time_varying_sel <- tv_sel
  dat
}

# Warnings raised while building the map, with build_params() noise excluded.
.ar1_map_warnings <- function(dat, random_sel = TRUE) {
  params <- suppressMessages(suppressWarnings(Rceattle::build_params(dat)))
  testthat::capture_warnings(
    suppressMessages(Rceattle::build_map(dat, params, random_sel = random_sel))
  )
}

testthat::test_that("an AR1 fleet does not warn for Time_varying_sel 'Off'", {
  for (sel in c("2DAR1", "3DAR1")) {
    # Entered as both the integer and the canonical string: switch_check()
    # turns 0 into "Off", so a workbook using either reaches the same guard.
    for (tv in list(0, "Off")) {
      w <- .ar1_map_warnings(.ar1_data(tv, sel))
      testthat::expect_length(grep("Time_varying_sel", w), 0)
    }
  }
})

testthat::test_that("an AR1 fleet still warns for a setting it overrides", {
  for (sel in c("2DAR1", "3DAR1")) {
    w <- .ar1_map_warnings(.ar1_data("IID", sel))
    hits <- grep("Time_varying_sel", w, value = TRUE)
    testthat::expect_gt(length(hits), 0)
    testthat::expect_match(hits[1], paste("ignored for", sel), fixed = TRUE)
  }
})

testthat::test_that("only 'Off' and an absent column are silent", {
  # The AR1 forms override every other mode, so the warned set is every
  # Time_varying_sel value except "Off". Pinned across the whole switch map so a
  # new mode cannot quietly join the silent side.
  for (sel in c("2DAR1", "3DAR1")) {
    for (tv in names(Rceattle:::tv_sel_map)) {
      n <- length(grep("Time_varying_sel", .ar1_map_warnings(.ar1_data(tv, sel))))
      testthat::expect_equal(n > 0, tv != "Off", info = paste(sel, tv))
    }
  }
})
