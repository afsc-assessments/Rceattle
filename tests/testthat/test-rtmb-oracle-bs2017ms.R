# The RTMB model's predation (R/6-rtmb_predation.R) against the C++ template on
# BS2017MS: the configurations stage 2 of the port covers (msmMode 1 and 2, suitMode 0
# and 4). The full gate is tools/verify/verify-rtmb-oracle.R fixture=bs2017ms; this keeps
# three of its points from drifting. The helpers are in helpers-rtmb-oracle.R.
#
# Decision 7 (5.67.11): NaN available food gives a NaN objective in both models. Pollock
# year-1 numbers underflow to ~1e-309 at init_dev -720 (inside its -1000 bound), so prey
# biomass, suitability and available food are NaN there.

test_that("Type II predation with empirical suitability matches the C++ at the golden parameters", {
  testthat::skip_on_cran()
  ref <- readRDS(test_path("fixtures", "golden-reference.rds"))
  expect_rtmb_matches(rtmb_pair(BS2017MS, ref$ms$params, msmMode = 1, suitMode = 0, niter = 5))
})

test_that("estimated suitability and the stomach likelihood match the C++", {
  testthat::skip_on_cran()
  ref <- readRDS(test_path("fixtures", "golden-reference.rds"))
  o <- rtmb_pair(BS2017MS, ref$ss$params, msmMode = 2, suitMode = 4, niter = 3)
  expect_rtmb_matches(o)
  stomach <- o$rtmb$report(o$cpp$par)$jnll_comp[19, 1:3]
  expect_true(all(stomach > 0))   # positive control: the stomach row is scored
})

test_that("NaN available food gives a NaN objective on both sides", {
  testthat::skip_on_cran()
  start <- rtmb_pair(BS2017MS, msmMode = 1, suitMode = 0, niter = 3)
  q <- start$params
  q$init_dev[1, seq_len(BS2017MS$nages[1] - 1)] <- -720
  o <- rtmb_pair(BS2017MS, q, msmMode = 1, suitMode = 0, niter = 3)
  x <- o$cpp$par
  expect_true(any(is.na(o$rtmb$report(x)$avail_food)))   # positive control
  expect_true(is.nan(o$cpp$fn(x)))
  expect_true(is.nan(o$rtmb$fn(x)))
})

test_that("the RTMB model refuses a predation configuration it does not implement", {
  testthat::skip_on_cran()
  o <- rtmb_pair(BS2017MS, msmMode = 1, suitMode = 4, niter = 3)
  refuses <- function(change, pattern) expect_error(.rtmb_prepare_data(change(o$data)), pattern)
  refuses(function(d) { d$msmMode <- 3; d }, "msmMode other than")
  refuses(function(d) { d$suitMode[] <- 2; d }, "suitMode other than")
  refuses(function(d) { d$suitMode[] <- c(4, 0, 0); d }, "suitMode differing")
  refuses(function(d) { d$diet_ctl[1, 7] <- 1990; d }, "Year other than 0")
  refuses(function(d) { d$niter <- 0; d }, "niter other than")
  refuses(function(d) { d$other_food[1] <- -1; d }, "other food of 0 or less")
  refuses(function(d) { d$diet_ll_type[] <- 1; d }, "Dirichlet-multinomial")
  refuses(function(d) { d$weight_obs[d$pop_wt_index[1], 1, 2, 1] <- -1; d }, "negative weight-at-age")
})
