# `Sel_norm_bin` says where selectivity is normalized to 1. It took a number that
# encoded three different intents with nothing in the value to say which: blank
# meant "do not normalize", anything negative meant "by the maximum", and a
# positive number was an absolute age (or a length bin). It now takes the words
# too, and refuses a bin the fleet is not selected over.
#
# The values every workbook in the ecosystem already carries -- blank, 0, -1,
# -999 and a positive bin -- must keep reaching TMB as exactly the same code.
# That parity is the point of most of this file.

.snb <- function(minage = 1, nages = 5, bfs = NA, dim = "Age", sel = NULL) {
  d <- make_test_data(nyrs = 6, nprojyrs = 2, nages = nages)
  d$growth_model <- rep(0, d$nspp)
  d$minage <- rep(minage, d$nspp)
  d$nages  <- rep(nages, d$nspp)
  # `bfs` is a bin ordinal: an age-based fleet takes it as the age it names.
  d$fleet_control$Bin_first_selected    <- if (dim == "Length") bfs else NA
  d$fleet_control$Age_first_selected    <- if (dim == "Length") NA else minage + bfs - 1
  d$fleet_control$Selectivity_dimension <- dim
  if (!is.null(sel)) d$fleet_control$Selectivity <- sel
  d
}

.code <- function(d, v, col = "Sel_norm_bin") {
  d$fleet_control[[col]] <- v
  out <- suppressMessages(suppressWarnings(rearrange_data(switch_check(d))))
  out[[if (col == "Sel_norm_bin") "sel_norm_bin1" else "sel_norm_bin2"]][1]
}


testthat::test_that("every legacy value reaches TMB as the code it always did", {
  d <- .snb()
  # -999 = do not normalize, -99 = by the maximum, >= 0 = that 0-based bin.
  # Pinned from the released behaviour before the words were added.
  expected <- list(`NA` = -999, `0` = -99, `-1` = -99, `-99` = -99,
                   `-999` = -99, `1` = 0, `3` = 2, `5` = 4)
  for (nm in names(expected)) {
    v <- if (nm == "NA") NA else as.numeric(nm)
    testthat::expect_equal(.code(d, v), expected[[nm]],
                           info = paste("Sel_norm_bin =", nm))
  }
})


testthat::test_that("the words mean what the schema says", {
  d <- .snb()
  for (w in c("Max", "max", "MAX", "Maximum")) {
    testthat::expect_equal(.code(d, w), -99, info = w)
  }
  for (w in c("Off", "off", "None", "none", "No", "NA", "")) {
    testthat::expect_equal(.code(d, w), -999, info = w)
  }
})


testthat::test_that("a value below the first selected bin normalizes by the maximum", {
  # selectivity.hpp zeroes the curve below the first selected bin, so a
  # reference taken there would divide by nothing. Age 3 is the first age that
  # can be the reference when it is the first selected age.
  d <- .snb(bfs = 3)
  testthat::expect_equal(.code(d, 1), -99)
  testthat::expect_equal(.code(d, 2), -99)
  testthat::expect_equal(.code(d, 3), 2)     # first selected age -> its 0-based bin
  testthat::expect_equal(.code(d, 5), 4)
})


testthat::test_that("age 0 is a real age on a stock that recruits at age 0", {
  # The rule is "below the fleet's first bin", not "at or below zero": with
  # minage = 0 the first age IS 0, so it is a reference bin and not a flag.
  f <- Rceattle:::.rce_sel_norm_code
  testthat::expect_equal(as.numeric(f(0, lo = 0)), 0, ignore_attr = TRUE)     # a bin
  testthat::expect_equal(as.numeric(f(0, lo = 1)), -1, ignore_attr = TRUE)    # the maximum, as always
  testthat::expect_equal(as.numeric(f(c(1, 2, 3), lo = 3)), c(-1, -1, 3), ignore_attr = TRUE)
  # One bound per fleet, not one for the whole column.
  testthat::expect_equal(as.numeric(f(c(0, 0, 0), lo = c(0, 1, 3))), c(0, -1, -1), ignore_attr = TRUE)
})


testthat::test_that("a bin past the last selected one is refused, naming the fleet", {
  d <- .snb()
  testthat::expect_error(.code(d, 9), "Survey")
  testthat::expect_error(.code(d, 9), "past the last age")
  # An unreadable word is refused rather than read as "do not normalize".
  testthat::expect_error(.code(d, "Maxx"), "not a valid age")
  testthat::expect_error(.code(d, "Maxx"), "Max, Off")
})


testthat::test_that("rearrange_data() refuses an unreadable value on its own", {
  # It is exported and reached directly, and switch_check() does not run on
  # every path into it. Without its own guard a typo resolved to NA and turned
  # normalization off silently.
  d <- suppressMessages(suppressWarnings(switch_check(.snb())))
  d$fleet_control$Sel_norm_bin <- "Maxx"
  testthat::expect_error(suppressMessages(rearrange_data(d)), "Maxx")
})


testthat::test_that("a saved workbook keeps the word", {
  d <- .snb()
  canon <- function(v) {
    d$fleet_control$Sel_norm_bin <- v
    suppressMessages(suppressWarnings(switch_check(d)))$fleet_control$Sel_norm_bin[1]
  }
  testthat::expect_equal(canon(NA), "Off")
  testthat::expect_equal(canon(-1), "Max")
  testthat::expect_equal(canon(-999), "Max")
  testthat::expect_equal(canon(0), "Max")     # 0 is below minage = 1
  testthat::expect_equal(canon("max"), "Max")
  testthat::expect_equal(canon("none"), "Off")
  testthat::expect_equal(canon(3), "3")

  # Canonicalizing an already-canonical workbook changes nothing.
  d$fleet_control$Sel_norm_bin <- c(-1, NA)
  once  <- suppressMessages(suppressWarnings(switch_check(d)))
  twice <- suppressMessages(suppressWarnings(switch_check(once)))
  testthat::expect_identical(once$fleet_control$Sel_norm_bin,
                             twice$fleet_control$Sel_norm_bin)
})


testthat::test_that("only a Fixed curve is exempt, not an Off fleet", {
  # A Fixed curve is read from emp_sel_obs and never normalized, so a stale
  # value there is harmless. An Off fleet is NOT exempt: selectivity.hpp skips
  # the normalizer on Selectivity, not on Fleet_type, so its bin still indexes
  # the array.
  fix <- .snb(sel = "Fixed"); fix$fleet_control$Sel_norm_bin <- 99
  testthat::expect_no_error(suppressMessages(suppressWarnings(switch_check(fix))))

  off <- .snb(); off$fleet_control$Fleet_type <- "Off"
  off$fleet_control$Selectivity  <- "Logistic"
  off$fleet_control$Sel_norm_bin <- 99
  testthat::expect_error(suppressMessages(suppressWarnings(switch_check(off))),
                         "past the last age")
})


testthat::test_that("LogisticPM takes All for its whole penalty range", {
  # On LogisticPM these columns are a penalty age-range, not a normalization
  # reference, and the model reads any negative there as "the whole selected
  # range". "All" is the word for that; "Max"/"Off" would describe something the
  # model does not do for this form.
  f <- Rceattle:::.rce_sel_norm_code
  testthat::expect_equal(as.numeric(f("All", lo = 1, allow_all = TRUE)), -1, ignore_attr = TRUE)
  testthat::expect_true(is.na(f("All", lo = 1, allow_all = FALSE)))
  testthat::expect_length(attr(f("All", lo = 1, allow_all = FALSE), "unrecognised"), 1)

  pm <- .snb(sel = "LogisticPM")
  testthat::expect_no_error(
    suppressMessages(suppressWarnings(switch_check(
      `$<-`(pm, "fleet_control", `$<-`(pm$fleet_control, "Sel_norm_bin", "All"))))))
})


testthat::test_that("fleets sharing a block are compared on meaning, not spelling", {
  # data_check() warns when fleets sharing a Selectivity_index differ in a
  # shaping column. "Max" and -1 are the same instruction, so they must not be
  # reported as a difference.
  d <- .snb(); d$fleet_control$Selectivity_index <- 1
  differs <- function(v) {
    d$fleet_control$Sel_norm_bin <- v
    hit <- FALSE
    withCallingHandlers(
      suppressMessages(try(data_check(d), silent = TRUE)),
      warning = function(w) {
        if (grepl("differ in Sel_norm", conditionMessage(w))) hit <<- TRUE
        invokeRestart("muffleWarning")
      })
    hit
  }
  testthat::expect_false(differs(c("Max", "-1")))
  testthat::expect_false(differs(c("Off", NA)))
  testthat::expect_false(differs(c("Max", "max")))
  testthat::expect_true(differs(c("3", "5")))     # a real difference still warns
})


testthat::test_that("a length-based fleet counts bins, not ages", {
  # Sel_norm_bin is a 1-based length-bin ordinal there, so minage does not enter.
  d <- .snb(dim = "Length")
  testthat::expect_equal(.code(d, 1), 0)
  testthat::expect_equal(.code(d, 3), 2)
  testthat::expect_equal(.code(d, NA), -999)
  testthat::expect_equal(.code(d, -1), -99)
})


testthat::test_that("a workbook using the deprecated column names behaves the same", {
  # Age_max_selected / Sel_norm_bin1 (and the _upper pair) are upgraded to the
  # canonical name before anything reads the column, so an old workbook gets the
  # words, the canonicalization and the range check identically. Pinned because
  # the upgrade and the new reads are in two different files.
  canon_and_code <- function(col, v) {
    d <- .snb()
    d$fleet_control$Sel_norm_bin <- NULL
    d$fleet_control[[col]] <- v
    o <- suppressMessages(suppressWarnings(switch_check(d)))
    list(canon = o$fleet_control$Sel_norm_bin[1],
         code  = suppressMessages(suppressWarnings(rearrange_data(o)))$sel_norm_bin1[1])
  }
  for (alias in c("Age_max_selected", "Sel_norm_bin1")) {
    testthat::expect_equal(canon_and_code(alias, "Max"),
                           list(canon = "Max", code = -99), info = alias)
    testthat::expect_equal(canon_and_code(alias, 3),
                           list(canon = "3", code = 2), info = alias)
    testthat::expect_equal(canon_and_code(alias, NA),
                           list(canon = "Off", code = -999), info = alias)
    testthat::expect_equal(canon_and_code(alias, -1),
                           list(canon = "Max", code = -99), info = alias)
  }

  d <- .snb()
  d$fleet_control$Sel_norm_bin_upper <- NULL
  d$fleet_control$Sel_norm_bin <- 2
  d$fleet_control$Age_max_selected_upper <- 4
  o <- suppressMessages(suppressWarnings(switch_check(d)))
  testthat::expect_equal(o$fleet_control$Sel_norm_bin_upper[1], "4")
  testthat::expect_equal(
    suppressMessages(suppressWarnings(rearrange_data(o)))$sel_norm_bin2[1], 3)
})


testthat::test_that("an Off fleet's bin is still bounded", {
  # selectivity.hpp skips the normalizer on Selectivity == "Fixed", NOT on
  # Fleet_type, so an Off fleet with an estimated form still indexes the
  # selectivity array with this value -- in a safebounds = FALSE build.
  d <- .snb(); d$fleet_control$Fleet_type <- "Off"
  d$fleet_control$Selectivity <- "Logistic"
  testthat::expect_error(.code(d, 99), "past the last age")
  testthat::expect_equal(.code(d, 3), 2)      # an in-range value still works
})


testthat::test_that("both entry points skip the same fleets", {
  # switch_check() and rearrange_data() must agree, or the checker passes a data
  # list the reshape then refuses.
  d <- .snb(sel = "Fixed")
  testthat::expect_no_error(.code(d, 99))
})


testthat::test_that("rearrange_data() enforces the ceiling on its own", {
  # It is exported and does not always run behind switch_check(). An unreadable
  # value fails safe; an over-large one indexes off the end of the array, so it
  # is the one that has to be refused here.
  d <- suppressMessages(suppressWarnings(switch_check(.snb())))
  d$fleet_control$Sel_norm_bin <- 99
  testthat::expect_error(suppressMessages(rearrange_data(d)), "past the last bin")
})


testthat::test_that("the upper bin is written back as Off, never Max", {
  # The model reads a negative upper bin as "no range" -- there is no maximum
  # for an upper bound, so "Max" would name something the column cannot do.
  d <- .snb()
  d$fleet_control$Sel_norm_bin <- 2
  d$fleet_control$Sel_norm_bin_upper <- -1
  o <- suppressMessages(suppressWarnings(switch_check(d)))
  testthat::expect_equal(o$fleet_control$Sel_norm_bin_upper[1], "Off")
  testthat::expect_equal(
    suppressMessages(suppressWarnings(rearrange_data(o)))$sel_norm_bin2[1], -999)
})


testthat::test_that("reinterpreting a value below the first selected bin is announced", {
  # The clamp this replaced printed a message. Moving where the curve is
  # anchored changes the scale of q, so it must not happen quietly.
  d <- .snb(bfs = 4)
  d$fleet_control$Sel_norm_bin <- 2
  testthat::expect_message(suppressWarnings(switch_check(d)),
                           "before the first age it is selected over")
})
