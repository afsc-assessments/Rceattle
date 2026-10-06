# `$` partial-matches on a list AND on a data.frame, silently, unless
# `warnPartialMatchDollar` is on -- which it is not by default. Where a schema
# column is absent and the prefix match is UNIQUE, the read hands back the
# longer column instead of NULL, so the model runs on a different setting than
# the workbook asked for. CLAUDE.md records the mechanism and the example: with
# `Time_varying_sel` absent, `fleet_control$Time_varying_sel` returns
# `Time_varying_sel_sd` -- an sd read as a mode.
#
# The schema has six names that are a prefix of another, ten (short, long)
# pairs in all. Three of the six have exactly ONE longer sibling, and those are
# the dangerous ones: an ambiguous prefix returns NULL and so fails loudly.
#
#   Time_varying_sel -> Time_varying_sel_sd   unique: resolves silently
#   Time_varying_q   -> Time_varying_q_sd     unique: resolves silently
#   Sel_norm_bin     -> Sel_norm_bin_upper    unique: resolves silently
#   Selectivity      -> Selectivity_{index,dimension}        ambiguous
#   Catchability     -> Catchability_{index,init,prior_sd}    ambiguous
#   CA               -> CAAL_{distribution,weights}           ambiguous
#
# `Time_varying_sel` and `Time_varying_q` carry no schema default, so they can
# genuinely be absent from a hand-built `fleet_control`.
#
# These tests turn the option on LOCALLY, not for the whole suite: set
# globally it fires inside readxl, dplyr and TMB for names this package does
# not own, which would be red on arrival and unfixable here. Warnings are then
# filtered by the name matched TO -- if the schema or a data list owns that
# name, the read is ours. Narrower and more honest than walking the call
# stack, and it is exactly the documented hazard.
#
# The right pattern already exists in the package: `.pull_int()` in
# R/5-rearrange_data.R uses `fc[[col]]`.

# Every name this package owns on an object a `$` read might target.
.rce_owned_names <- function(data_list = NULL) {
  sc <- Rceattle:::.rce_column_schema()
  nm <- vapply(sc, function(r) r$name, character(1), USE.NAMES = FALSE)
  if (!is.null(data_list)) {
    nm <- c(nm, names(data_list), names(data_list$fleet_control))
  }
  unique(nm[nzchar(nm)])
}

# Run `expr` with the partial-match warning on, and return the matched-to
# names of any partial match that landed on a name this package owns.
.rce_partial_match_hits <- function(expr, owned) {
  old <- options(warnPartialMatchDollar = TRUE)
  on.exit(options(old), add = TRUE)

  hits <- character()
  withCallingHandlers(
    force(expr),
    warning = function(w) {
      msg <- conditionMessage(w)
      m <- regmatches(
        msg, regexec("^partial match of '([^']+)' to '([^']+)'", msg))[[1]]
      if (length(m) == 3L) {
        if (m[3] %in% owned) hits <<- c(hits, paste0(m[2], " -> ", m[3]))
        invokeRestart("muffleWarning")
      }
    }
  )
  unique(hits)
}


# The detector must be able to see a real partial match, or every test below
# passes vacuously. Positive control: a data frame carrying only the longer
# sibling, read by the shorter name -- exactly a fleet_control missing
# `Time_varying_sel`.
testthat::test_that("the partial-match detector catches a known-bad read", {
  fc <- data.frame(Time_varying_sel_sd = c(0.2, 0.2))
  owned <- "Time_varying_sel_sd"

  hits <- .rce_partial_match_hits(invisible(fc$Time_varying_sel), owned)

  testthat::expect_length(hits, 1L)
  testthat::expect_identical(hits, "Time_varying_sel -> Time_varying_sel_sd")

  # And it must stay quiet on a name the package does not own, which is what
  # keeps it from firing inside readxl / dplyr / TMB.
  foreign <- data.frame(some_other_column = 1)
  testthat::expect_length(
    .rce_partial_match_hits(invisible(foreign$some_other), owned), 0L)
})


# The option must actually be off by default, or the premise of this file is
# wrong and the hazard would already be visible.
testthat::test_that("warnPartialMatchDollar is off by default", {
  testthat::expect_false(isTRUE(getOption("warnPartialMatchDollar")))
})


# The schema's prefix structure is worth pinning: a new column that creates a
# NEW unique prefix pair adds a silent-read hazard, and whoever adds it should
# have to acknowledge that here.
testthat::test_that("the schema's silent-prefix set is the documented three", {
  nm <- vapply(Rceattle:::.rce_column_schema(),
               function(r) r$name, character(1), USE.NAMES = FALSE)
  testthat::expect_gt(length(nm), 50)  # not vacuous on an empty schema

  sibs <- lapply(stats::setNames(nm, nm),
                 function(a) setdiff(nm[startsWith(nm, a)], a))
  sibs <- sibs[lengths(sibs) > 0]

  # Exactly one longer sibling resolves silently; more than one returns NULL.
  silent <- sort(names(sibs)[lengths(sibs) == 1L])
  testthat::expect_identical(
    silent, sort(c("Sel_norm_bin", "Time_varying_q", "Time_varying_sel")),
    info = paste0(
      "A schema column has created or removed a silent `$` prefix pair. ",
      "Read it with [[ ]] everywhere, then update this list. Current: ",
      paste(silent, collapse = ", ")))
})


# The net: the pipeline must make no partial-match read on a name the package
# owns, on any bundled data set.
testthat::test_that("data_check() makes no partial-match read it owns", {
  datasets <- c("BS2017SS", "BS2017MS", "GOA2018SS", "GOAatf",
                "NorthernRockfish2022", "GOApollock", "Atka2022", "GOAcod")
  checked <- 0L
  for (nm in datasets) {
    d <- tryCatch(get(nm, envir = asNamespace("Rceattle")),
                  error = function(e) NULL)
    if (is.null(d)) next
    checked <- checked + 1L
    hits <- .rce_partial_match_hits(
      suppressWarnings(suppressMessages(
        tryCatch(Rceattle:::data_check(d), error = function(e) NULL))),
      owned = .rce_owned_names(d))
    testthat::expect_identical(
      hits, character(0),
      info = paste0(nm, ": data_check() partially matched ",
                    paste(hits, collapse = "; ")))
  }
  # A loop that checked nothing would pass; fail if the data sets move.
  testthat::expect_gt(checked, 5L)
})


# WHY the net above is green -- the part actually worth protecting.
#
# There are ~50 bare `$` reads of the three silent-prefix names across R/.
# None is reachable with the column absent, because `switch_check()` fills all
# three from the schema before anything downstream reads them:
# `Sel_norm_bin`, `Time_varying_sel` and `Time_varying_q` all default to
# "Off" -- an absent column means no normalisation bin and no time variation.
#
# The fill itself has to read with `[[`. Reading the short name with `$` when
# it is absent returns the ONE longer sibling, so `.rce_apply_default()` is
# handed a non-NULL value and returns early: the default is never applied and
# the sibling's values are used instead. Measured on `Sel_norm_bin` with
# `Sel_norm_bin_upper` set to 7 -- a bare `$` yielded `7, 7, 7, ...` where the
# default is "Off". `Sel_norm_bin` is an absolute age, so that is a different
# selectivity normalisation on every fleet, which moves q and hence the
# advice.
#
# So the accessor in the FILL is load-bearing in a way the ~50 downstream
# reads are not: those run after the column exists. Rewriting one of these
# `[[` reads as `$` would silently restore the defect, which is what the
# block below pins.

testthat::test_that("an absent silent-prefix column gets its schema default", {
  # Not the sibling's values. On each column in turn: give the `_sd` or
  # `_upper` sibling a distinctive value, drop the short name, and require the
  # result to be the schema default rather than the sibling.
  sch <- .rce_column_schema()
  cases <- list(
    list(col = "Sel_norm_bin",     sib = "Sel_norm_bin_upper",  mark = 7),
    list(col = "Time_varying_sel", sib = "Time_varying_sel_sd", mark = 7),
    list(col = "Time_varying_q",   sib = "Time_varying_q_sd",   mark = 7))

  for (k in cases) {
    d <- Rceattle::GOA2018SS
    d$fleet_control[[k$sib]] <- k$mark
    d$fleet_control[[k$col]] <- NULL
    sc <- suppressWarnings(suppressMessages(Rceattle::switch_check(d)))
    got <- sc$fleet_control[[k$col]]

    testthat::expect_false(
      is.null(got),
      info = paste0(k$col, " was neither supplied nor defaulted"))
    # The sibling's marker must not appear anywhere in the short column.
    testthat::expect_false(
      any(as.character(got) == as.character(k$mark)),
      info = paste0(k$col, " took ", k$sib, "'s values -- the `$` partial ",
                    "match defeated the schema default"))
    # And it must be the schema default, whatever that is.
    testthat::expect_true(
      all(as.character(got) == as.character(sch[[k$col]]$default)),
      info = paste0(k$col, " is not its schema default (",
                    sch[[k$col]]$default, ")"))
  }
})


testthat::test_that("a column that is present keeps its own values", {
  # The default must fire only on absence. GOA2018SS carries a non-Off value
  # in each: Time_varying_sel is RandomWalkAscending on fleet 8,
  # Time_varying_q is RandomWalk on fleets 1 and 3.
  sc <- suppressWarnings(suppressMessages(
    Rceattle::switch_check(Rceattle::GOA2018SS)))
  tvs <- sc$fleet_control$Time_varying_sel
  testthat::expect_true("RandomWalkAscending" %in% tvs)
  testthat::expect_true("RandomWalk" %in% sc$fleet_control$Time_varying_q)
})
