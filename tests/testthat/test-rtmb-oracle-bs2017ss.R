# The RTMB model (R/6-rtmb_*.R) against the C++ template, on BS2017SS at msmMode = 0:
# the configuration stage 1 of the port covers. The full gate, with a whole optimizer
# path, perturbed inputs, derivatives and planted defects, is
# tools/verify/verify-rtmb-oracle.R; this keeps two of its points from drifting. The
# helpers are in helpers-rtmb-oracle.R.

test_that("the RTMB objective matches the C++ at the default start", {
  testthat::skip_on_cran()
  expect_rtmb_matches(rtmb_pair())
})

test_that("the RTMB objective matches the C++ at the golden parameters", {
  testthat::skip_on_cran()
  ref <- readRDS(test_path("fixtures", "golden-reference.rds"))
  expect_rtmb_matches(rtmb_pair(BS2017SS, inits = ref$ss$params))
})

# Each refusal below is an input the guard once let through, where the C++ gave a
# different objective (by 0.26 to 8.8 for catchability on a fishery, 3.4% for a catch row
# read under the wrong species, 0 against 1.5e6 at estimateMode 4).
test_that("the RTMB model refuses a configuration it does not implement", {
  testthat::skip_on_cran()
  o <- rtmb_pair()
  refuses <- function(change, pattern) {
    d <- o$data
    d <- change(d)
    expect_error(.rtmb_prepare_data(d), pattern)
  }
  refuses(function(d) { d$msmMode <- 3; d }, "msmMode other than")
  refuses(function(d) { d$flt_sel_type[1] <- 4; d }, "selectivity forms")
  refuses(function(d) { d$estimateMode <- 4; d }, "estimateMode 4")
  refuses(function(d) { d$est_index_q[1] <- 2; d }, "catchability forms")
  refuses(function(d) { d$est_index_q[1] <- 1; d$index_varying_q[1] <- 1; d }, "time-varying catchability")
  refuses(function(d) { d$catch_ctl[d$catch_ctl[, 1] == 1, 2] <- 2; d }, "catch row whose Species")
  refuses(function(d) { d$emp_sel_ctl[, 3] <- 2; d }, "emp_sel Sex 2")
  par <- o$cpp$env$parList(o$cpp$env$last.par.best)[names(o$cpp$env$parameters)]
  expect_error(.rtmb_make_adfun(o$data, par, random = "rec_dev"), "random effects")
})
