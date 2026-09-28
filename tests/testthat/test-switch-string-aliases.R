# =============================================================================
# Integer-coded model switches accept readable string aliases (the CIE review's
# "rename options to improve interpretability, e.g. constant -> not estimated").
# switch_check() resolves the strings to the integer codes the model consumes,
# so a string spec must produce a BIT-IDENTICAL fit to the integer spec. New
# aliases use a consistent convention: the not-estimated value is "Fixed",
# estimated is "Estimated".
# =============================================================================

testthat::test_that("estDynamics resolves readable strings to integer codes", {
  testthat::skip_if_not_installed("Rceattle")
  d <- Rceattle::BS2017SS
  for (pair in list(c("Estimated", "0"), c("Fixed", "1"), c("FixedScaled", "2"))) {
    d$estDynamics <- pair[1]
    out <- suppressMessages(Rceattle::switch_check(d))
    testthat::expect_equal(out$estDynamics, as.integer(pair[2]), info = pair[1])
  }
})

testthat::test_that("an unknown estDynamics string errors clearly", {
  testthat::skip_if_not_installed("Rceattle")
  d <- Rceattle::BS2017SS
  d$estDynamics <- "Frozen"
  testthat::expect_error(suppressMessages(Rceattle::switch_check(d)),
                         "Invalid 'estDynamics'")
})

testthat::test_that("Estimate_index_sd / Estimate_catch_sd resolve readable strings", {
  testthat::skip_if_not_installed("Rceattle")
  d <- Rceattle::BS2017SS
  n <- nrow(d$fleet_control)
  d$fleet_control$Estimate_index_sd <- rep("Estimated", n)
  d$fleet_control$Estimate_catch_sd <- rep("Analytical", n)
  out <- suppressMessages(Rceattle::switch_check(d))
  testthat::expect_equal(unique(out$fleet_control$Estimate_index_sd), 1L)
  testthat::expect_equal(unique(out$fleet_control$Estimate_catch_sd), 2L)
})

testthat::test_that("Estimate_index_sd passes NA (off fleets) through unmapped", {
  testthat::skip_if_not_installed("Rceattle")
  d <- Rceattle::BS2017SS
  col <- rep(c("Fixed", NA), length.out = nrow(d$fleet_control))
  d$fleet_control$Estimate_index_sd <- col
  out <- suppressMessages(Rceattle::switch_check(d))
  testthat::expect_equal(out$fleet_control$Estimate_index_sd[!is.na(col)],
                         rep(0L, sum(!is.na(col))))
  testthat::expect_true(all(is.na(out$fleet_control$Estimate_index_sd[is.na(col)])))
})

testthat::test_that("suitMode wires its (previously dead) string map", {
  testthat::skip_if_not_installed("Rceattle")
  d <- Rceattle::BS2017SS
  d$suitMode <- "GammaWeight"
  out <- suppressMessages(Rceattle::switch_check(d))
  testthat::expect_equal(out$suitMode, 2L)
})

testthat::test_that("build_srr resolves srr_est_mode strings", {
  testthat::skip_if_not_installed("Rceattle")
  testthat::expect_equal(suppressWarnings(Rceattle::build_srr(srr_est_mode = "Fixed"))$srr_est_mode, 0L)
  testthat::expect_equal(Rceattle::build_srr(srr_est_mode = "Estimated")$srr_est_mode, 1L)
  testthat::expect_equal(Rceattle::build_srr(srr_est_mode = "LognormalPrior")$srr_est_mode, 2L)
  testthat::expect_equal(Rceattle::build_srr(srr_est_mode = "BetaPrior")$srr_est_mode, 3L)
})

testthat::test_that("estDynamics string spec fits bit-identically to the integer spec", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  ctl <- Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0)
  d_int <- make_test_data(seed = 3)
  d_int$estDynamics <- 0L
  d_str <- d_int
  d_str$estDynamics <- "Estimated"

  f_int <- suppressMessages(Rceattle::fit_mod(
    data_list = d_int, estimateMode = 3, msmMode = 0, random_rec = FALSE, fit_control = ctl))
  f_str <- suppressMessages(Rceattle::fit_mod(
    data_list = d_str, estimateMode = 3, msmMode = 0, random_rec = FALSE, fit_control = ctl))
  testthat::expect_equal(sum(f_str$quantities$jnll_comp),
                         sum(f_int$quantities$jnll_comp), tolerance = 0)
})
