# Time_varying_sel_sd must be positive wherever the template reads it.
#
# Provenance: build_params() takes log(Time_varying_sel_sd), so a blank, zero or
# negative value became a -Inf/NaN starting value and the objective was not
# finite at the first evaluation -- reported from inside MakeADFun, naming
# neither the fleet nor the column. Time_varying_q_sd had this check; its
# selectivity twin did not.
#
# The gate is the (Selectivity, Time_varying_sel) PAIR, transcribed from the
# sel_dev_sd density sites in ceattle.cpp section 13.7 (SELECTIVITY PENALTIES),
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
})

# 2DAR1 and 3DAR1 score their deviation field at sel_dev_sd under every
# Time_varying_sel. When data_check() exempted them, a lead fleet with a blank
# sd built a model whose objective was NA (sd = 0: NaN) on GOA2018SS, with no error.
testthat::test_that("2DAR1 and 3DAR1 need the sd under every mode", {
  for (form in c("2DAR1", "3DAR1")) {
    for (tv in c("Off", "IID")) {
      for (bad in list(0, -1, NA_real_)) {
        testthat::expect_error(run_check(sel_dat(form, tv, bad)),
                               regexp = needs_sd, fixed = TRUE)
      }
      testthat::expect_no_error(run_check(sel_dat(form, tv, 0.5)))
    }
  }
  # An Off fleet scores nothing, so it is not asked.
  d <- sel_dat("2DAR1", "Off", NA_real_)
  d$fleet_control$Fleet_type[8] <- "Off"
  testthat::expect_no_error(run_check(d))
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


# Fleets sharing a Selectivity_index share one 2DAR1 block, and every member is
# asked, so a group cannot pass on its lead's value with a blank on a mirror.
testthat::test_that("every fleet of a shared 2DAR1 block needs the sd", {
  # GOA2018SS fleets 9 and 10 share Selectivity_index 8.
  shared <- function(sd9, sd10) {
    d <- sel_dat("2DAR1", "Off", sd9, fleet = 9L)
    d$fleet_control$Selectivity[10]         <- "2DAR1"
    d$fleet_control$Time_varying_sel[10]    <- "Off"
    d$fleet_control$Time_varying_sel_sd[10] <- sd10
    d
  }
  msg <- function(d) tryCatch({ run_check(d); "" },
                              error = function(e) conditionMessage(e))

  testthat::expect_match(msg(shared(0.5, NA_real_)),
                         "Fleet(s) ATF_bottom_trawl_length_comp need", fixed = TRUE)
  testthat::expect_match(msg(shared(NA_real_, 0.5)),
                         "Fleet(s) ATF_bottom_trawl need", fixed = TRUE)
  testthat::expect_match(msg(shared(0, 0)),
                         "Fleet(s) ATF_bottom_trawl, ATF_bottom_trawl_length_comp need",
                         fixed = TRUE)
  testthat::expect_no_error(run_check(shared(0.5, 0.5)))
})


testthat::test_that("the integer codes for 2DAR1 and 3DAR1 are checked too", {
  for (form in c("6", "7")) {
    testthat::expect_error(run_check(sel_dat(form, "0", NA_real_)),
                           regexp = needs_sd, fixed = TRUE)
  }
})


testthat::test_that("an AR1 or invalid mode on a 2DAR1 fleet gets its own refusal", {
  msg <- function(tv) tryCatch(run_check(sel_dat("2DAR1", tv, NA_real_)),
                               error = function(e) conditionMessage(e))
  m <- msg("AR1")
  testthat::expect_match(m, "Time_varying_sel = 'AR1' is removed", fixed = TRUE)
  testthat::expect_false(grepl(needs_sd, m, fixed = TRUE))
  m <- msg(NA)
  testthat::expect_match(m, "Invalid 'Time_varying_sel'", fixed = TRUE)
  testthat::expect_false(grepl(needs_sd, m, fixed = TRUE))
})
