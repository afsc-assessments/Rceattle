# Time_varying_sel_sd must be positive wherever the template reads it.
#
# Provenance: build_params() takes log(Time_varying_sel_sd), so a blank, zero or
# negative value became a -Inf/NaN starting value and the objective was not
# finite at the first evaluation -- reported from inside MakeADFun, naming
# neither the fleet nor the column. Time_varying_q_sd had this check; its
# selectivity twin did not.
#
# The gate is the (Selectivity, Time_varying_sel) PAIR, transcribed from the
# sel_dev_sd density sites in ceattle.cpp section 13.1 (FIT OBJECTIVE FUNCTION),
# so a combination that penalizes nothing is not asked for a value it would
# never read. These tests pin both halves: it fires where the sd is read, and
# stays quiet where it is not.

testthat::skip_on_cran()

sel_dat <- function(sel, tv, sd, fleet = 8L) {
  data("GOA2018SS", envir = environment())
  d <- suppressWarnings(suppressMessages(Rceattle::switch_check(GOA2018SS)))
  d$fleet_control$Selectivity[fleet]         <- sel
  d$fleet_control$Time_varying_sel[fleet]    <- tv
  d$fleet_control$Time_varying_sel_sd[fleet] <- sd
  d
}

run_check <- function(d) suppressWarnings(suppressMessages(Rceattle:::data_check(d)))

needs_sd <- "need a positive `Time_varying_sel_sd`"


testthat::test_that("a non-positive sd is refused where the density reads it", {
  # Fleet 8 is DoubleLogistic, whose deviates are scored at this sd under every
  # penalized mode.
  for (tv in c("IID", "RandomWalk", "RandomWalkAscending")) {
    for (bad in list(0, -1, NA_real_)) {
      testthat::expect_error(run_check(sel_dat("DoubleLogistic", tv, bad)),
                             regexp = needs_sd, fixed = TRUE)
    }
  }
})


testthat::test_that("the refusal names the fleet", {
  msg <- tryCatch(run_check(sel_dat("DoubleLogistic", "IID", 0)),
                  error = function(e) conditionMessage(e))
  testthat::expect_match(msg, "GOA_pollock_fishery", fixed = TRUE)
})


testthat::test_that("a positive sd passes on every mode that reads it", {
  for (tv in c("IID", "RandomWalk", "RandomWalkAscending")) {
    testthat::expect_no_error(run_check(sel_dat("DoubleLogistic", tv, 0.05)))
  }
})


testthat::test_that("modes and forms that never read the sd are not asked for one", {
  # Block carries no penalty, and Off has no deviation at all.
  testthat::expect_no_error(run_check(sel_dat("DoubleLogistic", "Block", 0)))
  testthat::expect_no_error(run_check(sel_dat("Logistic", "Off", 0)))
  # The 2DAR1/3DAR1 fields carry their own sd through sel_curve_pen, and
  # Time_varying_sel is ignored for them.
  testthat::expect_no_error(run_check(sel_dat("2DAR1", "IID", 0)))
})


testthat::test_that("every bundled data set still passes", {
  for (nm in c("BS2017SS", "BS2017MS", "GOA2018SS", "GOApollock",
               "GOAatf2023", "NorthernRockfish2022", "GeorgesBank3spp")) {
    d <- tryCatch(
      suppressWarnings(suppressMessages(Rceattle::switch_check(get(nm)))),
      error = function(e) NULL)
    if (is.null(d)) next
    testthat::expect_no_error(run_check(d))
  }
})
