# Under predation steepness does not exist, so .check_stock_recruit() used to
# return a NOTE without looking at the curve. A curve the data do not inform
# runs to a flat ridge (recruitment independent of SSB over the data) or a
# linear one (beta at zero), both with a positive-definite Hessian. The check
# now reads the curve over the hindcast SSB range and the +/-30 overflow bound
# (5.39.0). estimateMode = 3 evaluates the report at the starting values, so a
# curve is placed on a ridge or in the data through `inits`. build_srr()'s
# srr_alpha_init / srr_beta_init override supplied inits, so they are not used.

srr_dg_data <- function() {
  set.seed(123)
  make_msm_test_data()$data_list
}

srr_dg_build <- function(recFun, d = srr_dg_data(), inits = NULL) {
  suppressMessages(suppressWarnings(fit_mod(
    data_list = d, inits = inits, estimateMode = 3, msmMode = 1, suitMode = 0,
    initMode = "NonEquilibrium", random_rec = FALSE, recFun = recFun,
    fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE))))
}

# Rebuild with log alpha / log beta placed by hand.
srr_dg_at <- function(recFun, rec_pars, d = srr_dg_data()) {
  ini <- srr_dg_build(recFun, d)$initial_params
  ini$rec_pars[, 2:3] <- rec_pars
  srr_dg_build(recFun, d, inits = ini)
}

srr_check <- function(m) convergence_diagnostics(m)$checks$stock_recruit
hind_ssb  <- function(m) m$quantities$ssb[, seq_len(m$data_list$endyr - m$data_list$styr + 1)]

bh <- build_srr(srr_fun = "BevertonHolt")
rk <- build_srr(srr_fun = "Ricker")

testthat::test_that("a Beverton-Holt curve bending inside the SSB range passes", {
  testthat::skip_on_cran()
  # beta = 1 / median SSB puts predicted/asymptote at 0.5 in the middle of the
  # data; alpha keeps recruitment there near the fixture's mean.
  ssb <- hind_ssb(srr_dg_build(bh))
  m   <- srr_dg_at(bh, cbind(log(c(1, 10)), log(1 / apply(ssb, 1, stats::median))))
  chk <- srr_check(m)
  testthat::expect_identical(chk$severity, "OK")
  testthat::expect_match(chk$message, "bends within")
  testthat::expect_lt(chk$data$Species1$dd_at_smin, 0.9)
  testthat::expect_gt(chk$data$Species1$dd_at_smax, 0.1)
  # The record names its endpoints by SSB: the lower end is less saturated.
  testthat::expect_lt(chk$data$Species1$dd_at_smin, chk$data$Species1$dd_at_smax)
})

testthat::test_that("the flat and the linear ridge are both reported", {
  testthat::skip_on_cran()
  ssb  <- hind_ssb(srr_dg_build(bh))
  # Flat: beta so large that predicted/asymptote exceeds 0.9 at the lowest SSB.
  flat <- srr_check(srr_dg_at(bh, cbind(log(c(1, 10)), log(100 / apply(ssb, 1, min)))))
  testthat::expect_identical(flat$severity, "WARN")
  testthat::expect_match(flat$message, "flat over the observed SSB range")
  testthat::expect_match(flat$message, "Species1")
  testthat::expect_match(flat$message, "Species2")

  # Linear: beta at the overflow bound, no density dependence in the data.
  linear <- srr_check(srr_dg_at(bh, cbind(log(c(1, 10)), c(-30, -30))))
  testthat::expect_identical(linear$severity, "WARN")
  testthat::expect_match(linear$message, "linear over the observed SSB range")
  testthat::expect_match(linear$message, "overflow bound")
})

testthat::test_that("a Ricker curve peaking below the data, or with no bend, is reported", {
  testthat::skip_on_cran()
  ssb <- hind_ssb(srr_dg_build(rk))
  # Ricker beta is stored per 1e6 t: the peak sits at 1e6 / beta. Put it at a
  # tenth of the lowest observed SSB.
  chk <- srr_check(srr_dg_at(rk, cbind(log(c(1, 10)), log(1e6 / (0.1 * apply(ssb, 1, min))))))
  testthat::expect_identical(chk$severity, "WARN")
  testthat::expect_match(chk$message, "descending limb")
  lin <- srr_check(srr_dg_at(rk, cbind(log(c(1, 10)), c(-20, -20))))
  testthat::expect_match(lin$message, "linear over the observed SSB range")
})

testthat::test_that("a species with input numbers-at-age is skipped", {
  testthat::skip_on_cran()
  d   <- fixed_natage_data(c(1, 0))
  chk <- srr_check(fixed_natage_build(d, bh))
  testthat::expect_false("Species1" %in% names(chk$data))
  testthat::expect_true("Species2" %in% names(chk$data))
})
