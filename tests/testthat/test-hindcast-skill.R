# hindcast_skill(): forecast skill across retrospective peels.
#
# This answers a different question from Mohn's rho. Rho measures how the
# ESTIMATE of a year moves as data accumulate -- estimation consistency. This
# measures how well a peel PROJECTED years it could not see, given the catch
# that was actually taken, which is what a recruitment-projection assumption
# (proj_mean_rec TRUE/FALSE, or a DSEM) should be judged on.

# hindcast_skill() does not return the peel models, so recompute the one we need
# rather than asserting against numbers the function itself produced.
retro_peel_for <- function(hs, fit, peel) {
  r <- suppressWarnings(suppressMessages(
    Rceattle::retrospective(fit, peels = peel, cores = 1, getsd = FALSE,
                            forecast_rec = "model")))
  r$Rceattle_list[[1]]
}

testthat::test_that("hindcast_skill scores peels against the full model", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = Rceattle::BS2017SS, inits = NULL, file = NULL,
    estimateMode = 0, random_rec = TRUE, msmMode = 0,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0))))
  hs <- suppressWarnings(suppressMessages(
    Rceattle::hindcast_skill(fit, peels = 2, quantity = "ssb",
                             reference = "model", cores = 1, getsd = FALSE)))

  testthat::expect_true(all(c("mase", "by_year") %in% names(hs)))
  # Keyed by years_ahead, not peel: the scoring averages ACROSS peels at a fixed
  # horizon (Kell et al. 2021 eq. 5), so n_peels is the count that matters.
  testthat::expect_true(all(c("years_ahead", "species", "quantity", "n_peels",
                              "mae_forecast", "mae_naive", "mase") %in%
                              names(hs$mase)))
  testthat::expect_false("peel" %in% names(hs$mase))
  testthat::expect_gt(nrow(hs$mase), 0)

  # NOTE on what is worth asserting here. `years_ahead == year - (endyr - peel)`,
  # `year > endyr - peel` and `mase == mae_forecast/mae_naive` all restate the
  # lines that computed them and cannot fail; `all(is.finite(x) | is.na(x))`
  # passes for any numeric vector. They are omitted deliberately. What CAN fail
  # is whether the three series are the right quantities from the right objects,
  # which is what follows.
  endyr <- fit$data_list$endyr
  yr1 <- hs$by_year[hs$by_year$peel == 1, ]
  sp1 <- yr1[yr1$species == fit$data_list$spnames[1], ]
  col <- sp1$year - fit$data_list$styr + 1L

  # reference: the FULL model's estimate at those years.
  testthat::expect_equal(sp1$reference,
                         as.numeric(fit$quantities$ssb[1, col]),
                         tolerance = 1e-8)

  # forecast: the PEEL's own estimate at those years, and NOT the full model's.
  # Nothing tested this, so a forecast column silently taken from the wrong
  # object would have passed everything above.
  peel1 <- retro_peel_for(hs, fit, peel = 1)
  testthat::expect_equal(sp1$forecast,
                         as.numeric(peel1$quantities$ssb[1, col]),
                         tolerance = 1e-8)
  testthat::expect_false(isTRUE(all.equal(sp1$forecast, sp1$reference,
                                          tolerance = 1e-6)))

  # naive: one value held flat, and specifically the peel's TERMINAL estimate.
  testthat::expect_equal(length(unique(sp1$naive)), 1L)
  testthat::expect_equal(unique(sp1$naive),
                         as.numeric(peel1$quantities$ssb[1, (endyr - 1) -
                                                           fit$data_list$styr + 1L]),
                         tolerance = 1e-8)
})

testthat::test_that("hindcast_skill can score the held-out index instead", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = Rceattle::BS2017SS, inits = NULL, file = NULL,
    estimateMode = 0, random_rec = TRUE, msmMode = 0,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0))))
  hs <- suppressWarnings(suppressMessages(
    Rceattle::hindcast_skill(fit, peels = 2, reference = "observed",
                             cores = 1, getsd = FALSE)))

  idx <- hs$mase[grepl("^index_fleet", hs$mase$quantity), ]
  testthat::expect_gt(nrow(idx), 0)

  # The reference here must be the OBSERVED index, not a model quantity -- every
  # scored value has to appear in index_data at that fleet and year.
  by <- hs$by_year[grepl("^index_fleet", hs$by_year$quantity), ]
  obs <- fit$data_list$index_data
  # Match on FLEET and year, not year alone: with several fleets reporting in
  # the same year, a year-only match is satisfied by any of them and would pass
  # even if the reference were taken from the wrong fleet entirely.
  hit <- vapply(seq_len(min(nrow(by), 20L)), function(i) {
    flt <- as.integer(sub("^index_fleet_", "", by$quantity[i]))
    rows <- obs$Year == by$year[i] & obs$Fleet_code == flt
    any(rows) && any(abs(obs$Observation[rows] - by$reference[i]) < 1e-8)
  }, logical(1))
  testthat::expect_true(all(hit))

  # And the peel must be predicting years it did not fit.
  testthat::expect_true(all(by$years_ahead >= 1))
})


# MASE averages ACROSS PEELS at a fixed horizon -- Kell et al. (2021) eq. 5.
#
# Regression. It was grouped the other way round: by peel, averaging over
# horizons, which is the transpose of the published definition. Two consequences
# made it look like a property of the statistic rather than a defect. It reported
# one MASE per peel where Kell's Table 2 reports one per steps-ahead; and it made
# the h = 1 denominator a SINGLE |naive - reference| instead of a mean over
# peels, which is where the old "a one-year MASE is dominated by its own
# denominator, compare the deeper peels" advice came from. That advice was
# compensating for the grouping.
#
# Pinned by recomputing eq. 5 from $by_year rather than by restating the line
# that produced $mase: the point is the GROUPING, so the check has to do its own
# grouping.
testthat::test_that("MASE is computed per horizon, averaging over peels", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  d <- make_test_data()
  fit <- suppressMessages(suppressWarnings(fit_mod(
    data_list = d, file = NULL, estimateMode = 1,
    fit_control = fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
  hs <- suppressMessages(suppressWarnings(
    Rceattle::hindcast_skill(fit, peels = 2:4, quantity = "R",
                             reference = "model", cores = 1, getsd = FALSE)))

  by <- hs$by_year
  testthat::expect_gt(length(unique(by$peel)), 1L)   # or there is nothing to average

  # One row per (years_ahead, species, quantity) -- NOT per peel.
  key <- paste(by$years_ahead, by$species, by$quantity)
  testthat::expect_equal(nrow(hs$mase), length(unique(key)))

  # n_peels is the number of peels reaching that horizon, and it FALLS as the
  # horizon grows, because a horizon h needs a peel at least h deep.
  for (i in seq_len(nrow(hs$mase))) {
    z <- by[by$years_ahead == hs$mase$years_ahead[i] &
            by$species     == hs$mase$species[i] &
            by$quantity    == hs$mase$quantity[i], ]
    testthat::expect_equal(hs$mase$n_peels[i], nrow(z))
    testthat::expect_equal(hs$mase$n_peels[i], length(unique(z$peel)))
    # eq. 5, recomputed here
    testthat::expect_equal(hs$mase$mase[i],
                           mean(abs(z$forecast - z$reference)) /
                           mean(abs(z$naive    - z$reference)))
  }
  m1 <- hs$mase[hs$mase$years_ahead == min(hs$mase$years_ahead), ]
  mN <- hs$mase[hs$mase$years_ahead == max(hs$mase$years_ahead), ]
  testthat::expect_gt(m1$n_peels[1], mN$n_peels[1])
})
