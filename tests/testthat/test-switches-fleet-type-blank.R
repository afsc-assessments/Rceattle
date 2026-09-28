# Fleet_type has no schema default, so a blank one is not "unset, take the
# default" -- it is a fleet whose role in the likelihood nobody stated. Before
# 5.43.0 nothing refused it by name and the two halves of the package read it
# differently: `Fleet_type != "Off"` is NA rather than FALSE, so data_check()'s
# estimated-selectivity subset kept an all-NA row and died on
# "missing value where TRUE/FALSE needed" naming no fleet, while
# build_map_selectivity() treated NA as estimated (`.on[is.na(.on)] <- TRUE`).
# Found reviewing the 5.34.0-5.42.1 release (#158).

testthat::test_that("switch_check() refuses a blank Fleet_type and names the fleet", {
  d <- Rceattle::Atka2022

  for (blank in list(NA, NA_character_, "", "  ")) {
    bad <- d
    bad$fleet_control$Fleet_type[1] <- blank
    testthat::expect_error(suppressMessages(Rceattle:::switch_check(bad)),
                           "'Fleet_type' is blank for fleet\\(s\\) Bottom_trawl")
  }

  # The message names what to set, including the Off case, so the fix does not
  # require reading the schema.
  bad <- d; bad$fleet_control$Fleet_type[1] <- NA
  msg <- tryCatch(suppressMessages(Rceattle:::switch_check(bad)),
                  error = function(e) conditionMessage(e))
  testthat::expect_match(msg, "Fishery")
  testthat::expect_match(msg, "Survey")
  testthat::expect_match(msg, "Off")
})

testthat::test_that("a stated Fleet_type is untouched, in every spelling", {
  d <- Rceattle::Atka2022

  # Canonical names, integer codes and character codes all still pass AND all
  # resolve to the same thing -- the refusal must not narrow what is accepted,
  # and the three spellings must not diverge. Compared, not just run.
  got <- lapply(list(c("Fishery", "Survey"), c(1, 2), c("1", "2")), function(v) {
    ok <- d
    ok$fleet_control$Fleet_type <- rep(v, length.out = nrow(ok$fleet_control))
    suppressMessages(Rceattle:::switch_check(ok))$fleet_control$Fleet_type
  })
  testthat::expect_identical(got[[1]], rep(c("Fishery", "Survey"),
                                           length.out = nrow(d$fleet_control)))
  testthat::expect_identical(got[[2]], got[[1]])
  testthat::expect_identical(got[[3]], got[[1]])

  # An "Off" fleet is a stated role, not a blank one.
  off <- d; off$fleet_control$Fleet_type[1] <- "Off"
  testthat::expect_no_error(suppressMessages(Rceattle:::switch_check(off)))

  # A fleet_control with no Fleet_type column at all is a different case, and
  # the guard passes it through: nothing tolerates a missing column, so
  # dplyr::mutate() raises its own error one step later. Pinned so the guard is
  # never blamed for it, and so the message is noticed if it moves.
  none <- d; none$fleet_control$Fleet_type <- NULL
  testthat::expect_error(suppressMessages(Rceattle:::switch_check(none)),
                         "Fleet_type")
})

testthat::test_that("rearrange_data() refuses a blank Fleet_type as well", {
  # It is exported and switch_check() does not run on every path in, so before
  # 5.43.0 a blank reached the template as flt_type = NA.
  d <- suppressMessages(Rceattle:::switch_check(Rceattle::Atka2022))
  ok <- suppressMessages(suppressWarnings(Rceattle::rearrange_data(d)))
  testthat::expect_false(anyNA(ok$flt_type))

  d$fleet_control$Fleet_type[1] <- NA
  testthat::expect_error(suppressMessages(suppressWarnings(Rceattle::rearrange_data(d))),
                         "'Fleet_type' is blank for fleet\\(s\\)")
})

testthat::test_that("a Fleet_type outside the allowed set is refused", {
  # convert_switches() passes an out-of-range code through to the template:
  # 3 fits as a survey with its catch dropped, 2.7 truncates to Survey, and -1
  # contributes nothing AND leads its Selectivity_index group, so the group's
  # shape and curvature penalty is never charged. data_check() refuses all three,
  # so this closes the exported rearrange_data() path to match.
  d <- suppressMessages(Rceattle:::switch_check(Rceattle::Atka2022))
  for (v in list(3, -1, 2.7)) {
    bad <- d; bad$fleet_control$Fleet_type[1] <- v
    testthat::expect_error(
      suppressMessages(suppressWarnings(Rceattle::rearrange_data(bad))),
      "could not be read for fleet\\(s\\) Bottom_trawl")
  }

  # A factor is read downstream by its LEVEL INDEX, so Atka2022's Off survey and
  # its fishery reached the template as 2, 1 -- the Off fleet fitted as a fishery.
  # Resolved by label now. Every switch column takes the same route, so check one
  # more: the fishery's Selectivity must not read as its level index either.
  fac <- d
  fac$fleet_control$Fleet_type  <- factor(c("Off", "Fishery"))
  fac$fleet_control$Selectivity <- factor(c("DoubleLogistic", "Logistic"))
  out <- suppressMessages(suppressWarnings(Rceattle::rearrange_data(fac)))
  testthat::expect_identical(out$flt_type, c(0L, 1L))
  testthat::expect_identical(out$flt_sel_type, c(3L, 1L))
})

testthat::test_that("a factor Time_varying_q is read by label on the fit path", {
  # revert_switches() resolves eight of the nine switch columns, so switch_check()
  # covers them; Time_varying_q is the one it does not, and a factor there reached
  # the template by LEVEL INDEX: factor("Off") became IID (1), estimating
  # time-varying catchability deviations nobody asked for.
  d <- Rceattle::Atka2022
  d$fleet_control$Time_varying_q <- factor(c("Off", "Off"))
  out <- suppressMessages(suppressWarnings(
    Rceattle::rearrange_data(suppressMessages(Rceattle:::switch_check(d)))))
  testthat::expect_identical(out$index_varying_q, c(0L, 0L))
})

testthat::test_that("the de-factor rule covers every column it is handed", {
  fc <- data.frame(Fleet_name = factor(c("a", "b")), Fleet_type = factor(c("Off", "Fishery")),
                   Selectivity = c(2, 1), stringsAsFactors = FALSE)
  out <- Rceattle:::.rce_defactor_fleet_control(fc)
  testthat::expect_false(any(vapply(out, is.factor, logical(1))))
  testthat::expect_identical(out$Fleet_type, c("Off", "Fishery"))
  testthat::expect_identical(out$Selectivity, c(2, 1))   # non-factors untouched
  testthat::expect_null(Rceattle:::.rce_defactor_fleet_control(NULL))
})

testthat::test_that("the refusal names the row when Fleet_name is blank as well", {
  # The blank guard falls back to "row N" because Fleet_name has no default
  # either; this one has to do the same or it names nothing.
  d <- suppressMessages(Rceattle:::switch_check(Rceattle::Atka2022))
  d$fleet_control$Fleet_type[1] <- "Fisherie"
  d$fleet_control$Fleet_name[1] <- NA
  testthat::expect_error(
    suppressMessages(suppressWarnings(Rceattle::rearrange_data(d))),
    "could not be read for fleet\\(s\\) row 1")
})

testthat::test_that("a Fleet_type that cannot be read is named, not called blank", {
  # convert_switches() maps an unrecognised type to NA, so a typo used to be
  # reported as a blank cell -- the wrong thing to go looking for. data_check()
  # has validate_switches() for this; rearrange_data() is a separate entry point.
  d <- suppressMessages(Rceattle:::switch_check(Rceattle::Atka2022))
  d$fleet_control$Fleet_type[1] <- "Fisherie"
  err <- tryCatch(suppressMessages(suppressWarnings(Rceattle::rearrange_data(d))),
                  error = function(e) conditionMessage(e))
  testthat::expect_match(err, "could not be read")
  testthat::expect_match(err, "Fisherie")          # names the value ...
  testthat::expect_match(err, "Bottom_trawl")      # ... and the fleet
  testthat::expect_false(grepl("is blank", err))
})

testthat::test_that("the blank guard names the row when Fleet_name is blank too", {
  d <- Rceattle::Atka2022
  d$fleet_control$Fleet_type[1] <- NA
  d$fleet_control$Fleet_name[1] <- NA
  testthat::expect_error(suppressMessages(Rceattle:::switch_check(d)),
                         "blank for fleet\\(s\\) row 1")
})

testthat::test_that("every bundled data set still passes switch_check()", {
  # Named explicitly rather than fetched with get() from the namespace: the
  # bundled data are lazy-loaded, so the namespace environment does not hold
  # them under every reporter, and the dynamic form passed alone and errored
  # inside the full suite.
  bundled <- list(BS2017SS = Rceattle::BS2017SS, BS2017MS = Rceattle::BS2017MS,
                  GOA2018SS = Rceattle::GOA2018SS, GOApollock = Rceattle::GOApollock,
                  GOAatf = Rceattle::GOAatf, GOAcod = Rceattle::GOAcod,
                  Atka2022 = Rceattle::Atka2022,
                  GeorgesBank3spp = Rceattle::GeorgesBank3spp,
                  NorthernRockfish2022 = Rceattle::NorthernRockfish2022)
  for (nm in names(bundled)) {
    testthat::expect_no_error(suppressMessages(Rceattle:::switch_check(bundled[[nm]])))
  }
})
