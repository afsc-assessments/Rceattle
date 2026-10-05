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
