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
# There are ~50 bare `$` reads of the three silent-prefix names across R/,
# against 2 `[[ ]]` reads. None is reachable with the column absent, because
# something louder happens first:
#
#   * `switch_check()` and `data_check()` reach these columns through the dplyr
#     pronoun, `.data$Time_varying_q`, and rlang's pronoun ERRORS on a missing
#     column rather than partial-matching. Measured: dropping
#     `Time_varying_sel` makes `switch_check()` abort with "Column
#     `Time_varying_sel` not found in `.data`".
#   * where the pronoun is not the first reader, `switch_check()` has already
#     supplied the column, so the bare `$` reads downstream find it.
#
# The protection is therefore a property of WHICH ACCESSOR COMES FIRST, not of
# the accessors themselves. Rewriting a `.data$` read as a bare `$` -- an
# innocuous-looking tidy-up -- would remove the loud failure and expose the
# silent one. This test pins the loudness: dropping a silent-prefix column must
# error, not quietly substitute.
testthat::test_that("dropping a silent-prefix column fails loudly", {
  for (col in c("Time_varying_sel", "Time_varying_q", "Sel_norm_bin")) {
    d <- Rceattle::BS2017SS
    d$fleet_control[[col]] <- NULL

    sc <- tryCatch(
      suppressWarnings(suppressMessages(Rceattle::switch_check(d))),
      error = function(e) e)

    # Either switch_check() refuses the missing column, or it supplies it.
    # Both are loud enough; a silent wrong read downstream is not.
    if (!inherits(sc, "condition")) {
      testthat::expect_true(
        col %in% names(sc$fleet_control),
        info = paste0(
          col, ": switch_check() neither refused the missing column nor ",
          "supplied it, so downstream `$` reads of it are a silent partial ",
          "match waiting for an input."))
    } else {
      testthat::succeed()
    }
  }
})
