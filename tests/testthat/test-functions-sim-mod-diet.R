# Stomach contents are drawn by the TMB model's SIMULATE block, beside the diet
# density (ceattle.cpp, section 13.6). Before that they were the one observation
# type never simulated, so a multispecies self_test() refit against the same
# stomachs every replicate and recovered suitability from data that never varied.
#
# The check that matters is the ROUND TRIP: take the simulated stomach
# proportions, put them through the transform a refit applies (append the "other
# prey" balance, add comp_offset, renormalize), and compare their mean to the
# proportions the density predicts. That is what the estimator actually sees.
#
# Comparing the RAW simulated proportions to the offset-normalized prediction
# instead -- the obvious-looking check -- cannot detect a draw that is centred on
# the wrong one of the two, because at the default comp_offset the two targets
# differ by only ~0.02%. Raise comp_offset to give the test power.

# The fixture and helpers live in helpers-diet-sim.R; the expensive round-trip,
# over-dispersion and weight tests live in the sibling files named for them.
testthat::skip_on_cran()

testthat::test_that("simulated stomach contents are a valid diet composition", {
  testthat::skip_if_not_installed("TMB")
  # 28 s here, but covr rebuilds the TMB model at -O0 and this multispecies
  # simulate-and-refit is the file's dominant cost there: shard 5 of the coverage
  # run spent ~100 of its 120 minutes inside this file and was cancelled before
  # writing a timing row for it. The two draw tests below are the same code path
  # at a tenth of the cost.
  testthat::skip_on_covr()

  mod <- .diet_sim_fit("Multinomial", 1)
  set.seed(2)
  sim <- suppressWarnings(Rceattle::sim_mod(mod, simulate = TRUE))
  x <- sim$diet_data$Stomach_proportion_by_weight

  testthat::expect_equal(nrow(sim$diet_data), nrow(mod$data_list$diet_data))
  testthat::expect_false(anyNA(x))
  testthat::expect_true(all(x >= 0))

  # Prey proportions are a share of the stomach; the balance is "other prey",
  # which is not stored and is recomputed on the next fit. Tested at exactly 1,
  # matching data_check(): a looser bound here would pass tables the pipeline
  # rejects, and a rejected refit surfaces as a failed self_test() replicate
  # rather than as bad data.
  sid <- Rceattle::rearrange_data(mod$data_list)$stomach_id
  testthat::expect_true(all(tapply(x, sid, sum) <= 1))

  # The draw actually varies -- diet used to pass through untouched.
  testthat::expect_gt(sum(abs(x - mod$data_list$diet_data$Stomach_proportion_by_weight) > 1e-12), 0)
})


testthat::test_that("a model without an estimated suitability leaves diet alone", {
  testthat::skip_if_not_installed("TMB")

  # suitMode = 0 is empirical suitability: the diet likelihood is not evaluated,
  # nothing is predicted, and so nothing should be redrawn.
  utils::data("BS2017SS", package = "Rceattle", envir = environment())
  ss <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    BS2017SS, inits = NULL, msmMode = 0, estimateMode = 1,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))

  set.seed(3)
  sim <- suppressWarnings(Rceattle::sim_mod(ss, simulate = TRUE))
  testthat::expect_identical(sim$diet_data, ss$data_list$diet_data)
})


testthat::test_that("a fractional Sample_size still yields a valid composition", {
  testthat::skip_if_not_installed("TMB")
  # The file's other 27-second fit, and the same -O0 story as the first block.
  testthat::skip_on_covr()

  # The draw places a whole number of observations, so dividing the counts by a
  # fractional Sample_size would give proportions summing to round(N)/N > 1 --
  # which data_check() rejects, and self_test() then reports as a replicate that
  # failed to converge. Divide by the number actually placed instead.
  utils::data("BS2017SS", package = "Rceattle", envir = environment())
  utils::data("BS2017MS", package = "Rceattle", envir = environment())
  ss <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    BS2017SS, inits = NULL, msmMode = 0, estimateMode = 1,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
  d <- BS2017MS
  d$diet_data$Sample_size <- 20.6
  mod <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    d, inits = ss$estimated_params, msmMode = 1, suitMode = 4,
    estimateMode = 1, niter = 3,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))

  sid <- Rceattle::rearrange_data(mod$data_list)$stomach_id
  set.seed(12)
  worst <- 0
  for (i in 1:15) {
    x <- suppressWarnings(Rceattle::sim_mod(mod, simulate = TRUE))$diet_data$Stomach_proportion_by_weight
    worst <- max(worst, max(tapply(x, sid, sum)))
  }
  testthat::expect_lte(worst, 1)
})


testthat::test_that("a diet row whose weighted sample size rounds to zero comes back empty", {
  testthat::skip_if_not_installed("TMB")

  # Regression: the write-back used to be gated on the draw having placed
  # something, so a stomach whose N_s * Diet_comp_weights fell below one
  # observation silently kept its OBSERVED proportions. A self_test() then
  # refitted against the real diet on those rows and read suitability recovery
  # as better than it was. The per-predator "frozen" warning could not see it --
  # a predator with some stomachs drawn and some rounded away is not frozen.
  # Compositions have always zeroed such a row; diet now does too.
  mod <- .diet_sim_fit("Multinomial", 1e-4)

  obs <- mod$data_list$diet_data$Stomach_proportion_by_weight
  sim <- suppressWarnings(
    Rceattle::sim_mod(mod, simulate = TRUE))$diet_data$Stomach_proportion_by_weight

  # Every row the model fits is emptied, not carried over.
  had <- obs > 0
  testthat::expect_gt(sum(had), 0)
  testthat::expect_true(all(sim[had] == 0))

  # And it is announced rather than silent.
  testthat::expect_warning(Rceattle::sim_mod(mod, simulate = TRUE),
                           "came back empty")
})
