# Unit tests for the MASE aggregation behind `retrospective()$mase` and
# `hindcast_skill()$mase`. These run on a constructed `by_year` frame rather
# than a fit, because the arithmetic is what has been wrong twice: once by
# grouping on the wrong axis (per peel instead of per horizon, which is not
# Kell et al. 2021 eq. 5), and once by letting a single NA null a whole horizon
# while `n_peels` still advertised every peel. Neither needed an optimization to
# reproduce, and an end-to-end fit is covered by test-hindcast-skill.R.

by_year <- function(...) {
  d <- data.frame(..., stringsAsFactors = FALSE)
  stopifnot(all(c("peel", "years_ahead", "species", "quantity",
                  "forecast", "reference", "naive") %in% names(d)))
  d
}

testthat::test_that("MASE averages across peels at a fixed horizon, not within a peel", {
  # Two peels, two horizons. Eq. 5 sums over t -- the peels -- at fixed h, so
  # horizon 1 must pool peels 1 and 2, giving 2 rows of 2 peels rather than
  # 2 rows of 2 horizons.
  d <- by_year(
    peel        = c(1, 1, 2, 2),
    years_ahead = c(1, 2, 1, 2),
    species     = 1, quantity = "R",
    forecast    = c(11, 12, 13, 14),
    reference   = c(10, 10, 10, 10),
    naive       = c(12, 14, 16, 18))
  m <- Rceattle:::.rce_mase_aggregate(d)

  testthat::expect_identical(nrow(m), 2L)
  testthat::expect_identical(m$years_ahead, c(1, 2))
  testthat::expect_identical(m$n_peels, c(2L, 2L))
  # horizon 1: |11-10| and |13-10| -> 2; naive |12-10| and |16-10| -> 4
  testthat::expect_equal(m$mase[m$years_ahead == 1], 2 / 4)
  # horizon 2: |12-10| and |14-10| -> 3; naive |14-10| and |18-10| -> 6
  testthat::expect_equal(m$mase[m$years_ahead == 2], 3 / 6)
})

testthat::test_that("an unscoreable peel drops out and the counts say so", {
  # Peel 2 has no forecast. The horizon must still be scored on peel 1 alone,
  # and n_peels / peels_used must describe the peel that actually contributed
  # -- the failure was reporting n_peels = 2 beside a number from one peel, or
  # nulling the row entirely.
  d <- by_year(
    peel        = c(1, 2, 3),
    years_ahead = 1, species = 1, quantity = "SSB",
    forecast    = c(11, NA_real_, 13),
    reference   = 10,
    naive       = c(12, 16, 14))
  m <- Rceattle:::.rce_mase_aggregate(d)

  testthat::expect_identical(nrow(m), 1L)
  testthat::expect_identical(m$n_peels, 2L)
  testthat::expect_identical(m$peels_used, "1,3")
  # mean(|11-10|, |13-10|) / mean(|12-10|, |14-10|) = 2 / 3
  testthat::expect_equal(m$mase, 2 / 3)
})

testthat::test_that("a horizon with no scoreable peel is NA, not a mean of nothing", {
  d <- by_year(
    peel = 1, years_ahead = 4, species = 1, quantity = "SSB",
    forecast = NA_real_, reference = 10, naive = NA_real_)
  m <- Rceattle:::.rce_mase_aggregate(d)

  testthat::expect_identical(m$n_peels, 0L)
  testthat::expect_identical(m$peels_used, "")
  testthat::expect_true(is.na(m$mase))
})

testthat::test_that("a naive baseline of exactly 0 gives NA rather than Inf", {
  # Persistence was exactly right. An undefined ratio is not infinitely bad
  # skill, and an Inf would poison any mean taken over the column.
  d <- by_year(
    peel = 1, years_ahead = 1, species = 1, quantity = "R",
    forecast = 11, reference = 10, naive = 10)
  m <- Rceattle:::.rce_mase_aggregate(d)

  testthat::expect_equal(m$mae_naive, 0)
  testthat::expect_true(is.na(m$mase))
})

testthat::test_that("rows are kept apart by species and quantity as well as horizon", {
  d <- by_year(
    peel        = rep(1:2, each = 4),
    years_ahead = 1,
    species     = rep(c(1, 1, 2, 2), 2),
    quantity    = rep(c("SSB", "R"), 4),
    forecast    = 11, reference = 10, naive = 12)
  m <- Rceattle:::.rce_mase_aggregate(d)

  # 1 horizon x 2 species x 2 quantities, each resting on both peels.
  testthat::expect_identical(nrow(m), 4L)
  testthat::expect_identical(sort(unique(m$species)), c(1, 2))
  testthat::expect_identical(sort(unique(m$quantity)), c("R", "SSB"))
  testthat::expect_true(all(m$n_peels == 2L))
})

testthat::test_that("two rows from one peel at one horizon count as one peel", {
  # A fleet may legally carry two index rows in the same Year at different
  # Months, so both land in the same years_ahead from the same peel. Averaging
  # within the peel first stops that peel being weighted twice, and keeps
  # n_peels a count of PEELS rather than of rows.
  d <- by_year(
    peel = c(1, 1), years_ahead = 1, species = 1, quantity = "R",
    forecast = c(11, 13), reference = 10, naive = 12)
  m <- Rceattle:::.rce_mase_aggregate(d)

  testthat::expect_identical(m$n_peels, 1L)
  testthat::expect_identical(m$peels_used, "1")
  # mean(|11-10|, |13-10|) = 2, not a sum of 4
  testthat::expect_equal(m$mae_forecast, 2)
})

# A supplied `retro` already holds its forecast years, computed under the rule
# its own call was given, and hindcast_skill() recomputes none of them. So
# scoring one retro twice under two values of `forecast_rec` returned the SAME
# numbers both times, which reads as "the projection rule makes no difference"
# -- the conclusion the comparison exists to reach -- rather than as an argument
# that was ignored. This is the shape a user writes when told both rules can be
# scored from one fit: true of the FIT, false of the PEELS.
testthat::test_that("hindcast_skill() refuses a forecast_rec its retro contradicts", {
  fake <- structure(list(
    Rceattle_list = list(structure(list(), class = "Rceattle"),
                         structure(list(), class = "Rceattle")),
    mase = data.frame(forecast_rec = "mean", years_ahead = 1,
                      stringsAsFactors = FALSE)),
    class = "Rceattle_retro")
  obj <- structure(list(quantities = list(R = 1), data_list = list()),
                   class = "Rceattle")

  testthat::expect_error(
    Rceattle::hindcast_skill(obj, retro = fake, quantity = "R",
                             forecast_rec = "model"),
    "cannot be applied to a `retro` that was already fitted")
  # The error has to say what to do instead, or it just blocks the comparison.
  testthat::expect_error(
    Rceattle::hindcast_skill(obj, retro = fake, quantity = "R",
                             forecast_rec = "model"),
    "one retrospective under each")
  # Agreeing must get PAST this guard. The stub is not a real fit, so the call
  # still fails further in -- what matters is that it no longer fails HERE, so
  # assert on the message rather than on there being no error at all.
  agree <- tryCatch(
    suppressWarnings(Rceattle::hindcast_skill(obj, retro = fake,
                                              quantity = "R",
                                              forecast_rec = "mean")),
    error = conditionMessage)
  testthat::expect_false(grepl("cannot be applied to a `retro`", agree,
                               fixed = TRUE))
  # Same for leaving it unset: `missing(forecast_rec)` is what gates the guard.
  unset <- tryCatch(
    suppressWarnings(Rceattle::hindcast_skill(obj, retro = fake,
                                              quantity = "R")),
    error = conditionMessage)
  testthat::expect_false(grepl("cannot be applied to a `retro`", unset,
                               fixed = TRUE))
})


# `N` is the number of observations the statistic is computed from, in BOTH
# $mohns and $mase. Mohn's rho used to count a peel that contributed a
# NON-FINITE relative error: `(peel - base)/base` is NaN or Inf where the full
# model's base is 0 -- an unfished year for `F_spp`, a collapsed stock for `ssb`
# or `R` -- and summing that made rho NaN for every peel at that horizon, for
# that species, while N still advertised every peel. The bundled data cannot
# reach a zero base, so the base is injected here.
testthat::test_that("Mohn's rho counts only the peels it could score", {
  testthat::skip_on_cran()
  fit <- suppressMessages(suppressWarnings(fit_mod(
    data_list = Rceattle::BS2017SS, file = NULL, estimateMode = 1,
    fit_control = fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))

  # The full model is the BASE of every relative error, and rho at forecast year
  # 0 reads the column of the PEEL's terminal year -- a different column per
  # peel -- so the whole species row is zeroed rather than one cell. That makes
  # species 2's F_spp unscoreable at every horizon and leaves every other
  # species and quantity untouched, which is the half that used to break.
  bad <- fit
  bad$quantities$F_spp[2, ] <- 0

  r <- suppressMessages(suppressWarnings(
    retrospective(bad, peels = 2:3, nyrs_forecast = 1, getsd = FALSE,
                  cores = 1)))
  f0 <- r$mohns[r$mohns$Object == "F_spp" & r$mohns[["Forecast year"]] == 0, ]
  testthat::expect_identical(nrow(f0), 3L)

  sp <- fit$data_list$spnames
  hit  <- f0[f0$species == sp[2], ]
  rest <- f0[f0$species != sp[2], ]

  # The unscoreable species: no contributing observation, so no answer -- and N
  # says 0 rather than claiming the peels it could not use.
  testthat::expect_identical(hit$N, 0L)
  testthat::expect_true(is.na(hit$rho))

  # Its neighbours are unaffected. This is the half that used to break: one
  # species' zero base does not reach another's rho.
  testthat::expect_true(all(rest$N > 0L))
  testthat::expect_true(all(is.finite(rest$rho)))

  # And N never exceeds the peels that reached the horizon.
  testthat::expect_true(all(r$mohns$N <= length(r$Rceattle_list) - 1L))
})


# $mohns was WIDE before 5.23.0.9007 -- Object | Forecast year | N | one column
# per species -- and retrospectives are saved to disk and reloaded months later
# (GOA-multispecies-assessment/R/03_diagnostics.R writes retrospectives.RData
# and 07_figures_tables.R reloads it). print() read `rho`, which is NULL on that
# shape, so `rho` was numeric(0), no value could be outside the band, and the
# header said "status: OK" for a 55% retrospective bias before erroring further
# down. A false clean bill of health is worse than a failure.
testthat::test_that("print() reads a $mohns saved in the old wide shape", {
  wide <- data.frame(Object = c("biomass", "ssb"), "Forecast year" = c(0, 0),
                     N = c(3L, 3L), Pollock = c(0.55, 0.61),
                     check.names = FALSE)
  obj <- structure(list(Rceattle_list = list(1, 2, 3, 4), mohns = wide,
                        peels_requested = 3L, peel_depths = 1:3),
                   class = "Rceattle_retro")

  out <- paste(utils::capture.output(print(obj)), collapse = "\n")
  testthat::expect_match(out, "status: WARN")
  testthat::expect_match(out, "2 of 2")
  testthat::expect_match(out, "0.55")
  testthat::expect_match(out, "0.61")

  # And the reshape itself, so the rows are not merely counted but correct.
  long <- Rceattle:::.rce_mohns_as_long(wide)
  testthat::expect_identical(names(long),
    c("Object", "Forecast year", "N", "species", "rho"))
  testthat::expect_identical(long$rho, c(0.55, 0.61))
  testthat::expect_identical(unique(long$species), "Pollock")
  # A table already long is returned untouched.
  testthat::expect_identical(Rceattle:::.rce_mohns_as_long(long), long)
})

testthat::test_that("print() refuses a $mohns shape it cannot read", {
  # No `rho` AND no species column to derive one from. Printing a band verdict
  # over zero values would report OK.
  bad <- structure(list(
    Rceattle_list = list(1, 2),
    mohns = data.frame(Object = "ssb", "Forecast year" = 0, N = 3L,
                       check.names = FALSE)),
    class = "Rceattle_retro")
  testthat::expect_error(print(bad), "not a shape print\\(\\) recognizes")
})

# `$mase` is NULL whenever no peel converged, and absent entirely on a retro
# saved before it existed -- which is exactly where a silently ignored
# `forecast_rec` does the most damage. The rule is therefore recorded on the
# object, not only on the table.
testthat::test_that("the forecast_rec guard reads the object, not just $mase", {
  fake <- structure(list(
    Rceattle_list = list(structure(list(), class = "Rceattle"),
                         structure(list(), class = "Rceattle")),
    mase = NULL,                      # no peel scored
    forecast_rec = "mean"),
    class = "Rceattle_retro")
  obj <- structure(list(quantities = list(R = 1), data_list = list()),
                   class = "Rceattle")

  testthat::expect_error(
    suppressWarnings(Rceattle::hindcast_skill(obj, retro = fake, quantity = "R",
                                              forecast_rec = "model")),
    "cannot be applied to a `retro` that was already fitted")
})
