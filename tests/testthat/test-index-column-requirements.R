# Each q column is required by the switch that reads it -- and only by that one.
#
# All three are taken as log() when the parameter list is built, so a missing or
# non-positive value gives NaN or -Inf and the objective is non-finite at the
# first evaluation. It fails loudly, but from inside MakeADFun, whose message
# names neither the fleet nor the column. Each case below was confirmed to
# produce a non-finite objective before the check existed.
testthat::skip_on_cran()

.q_fixture <- function(tweak) {
  d <- make_test_data(nyrs = 12, nages = 5, seed = 42)
  i <- which(d$fleet_control$Fleet_type == "Survey")[1]
  d$fleet_control$Catchability[i]      <- "Estimated"
  d$fleet_control$Catchability_init[i] <- 1
  tweak(d, i)
}
.fit <- function(d) suppressMessages(suppressWarnings(Rceattle::fit_mod(
  d, file = NULL, estimateMode = 3, msmMode = 0, random_rec = FALSE,
  fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0, phase = FALSE))))

testthat::test_that("a missing sd is caught for the switch that reads it", {
  testthat::skip_if_not_installed("TMB")

  # Estimate_index_sd = 1 starts the estimated observation sd at log(Index_sd).
  testthat::expect_error(.fit(.q_fixture(function(d, i) {
    d$fleet_control$Estimate_index_sd[i] <- 1
    d$fleet_control$Index_sd[i] <- NA
    d })), "Index_sd")

  # The catchability prior is scored at Catchability_prior_sd.
  testthat::expect_error(.fit(.q_fixture(function(d, i) {
    d$fleet_control$Catchability[i] <- "Estimated-with-prior"
    d$fleet_control$Catchability_prior_sd[i] <- NA
    d })), "Catchability_prior_sd")

  # AR1 catchability estimates the sd, starting it at log(Catchability_prior_sd).
  testthat::expect_error(.fit(.q_fixture(function(d, i) {
    d$fleet_control$Catchability[i] <- "AR1"
    d$fleet_control$Catchability_prior_sd[i] <- NA
    d })), "Catchability_prior_sd")

  # Time-varying q deviations are penalized at Time_varying_q_sd.
  for (tv in c("IID", "RandomWalk")) {
    testthat::expect_error(.fit(.q_fixture(function(d, i) {
      d$fleet_control$Time_varying_q[i] <- tv
      d$fleet_control$Time_varying_q_sd[i] <- NA
      d })), "Time_varying_q_sd", info = tv)
  }

  # log(Catchability_init) has to be finite, so zero is rejected as well as NA.
  testthat::expect_error(.fit(.q_fixture(function(d, i) {
    d$fleet_control$Catchability_init[i] <- 0
    d })), "Catchability_init")
})

testthat::test_that("a time BLOCK needs no sd, and valid settings still fit", {
  testthat::skip_if_not_installed("TMB")

  # Block is time-varying but carries no penalty, so it reads no sd. Requiring
  # one here would reject a working configuration -- the reason the check keys
  # on the switch value rather than on "is time-varying".
  fit <- .fit(.q_fixture(function(d, i) {
    d$fleet_control$Time_varying_q[i] <- "Block"
    d$fleet_control$Time_varying_q_sd[i] <- NA
    d }))
  testthat::expect_true(is.finite(fit$obj$fn()))

  # Analytical and AnalyticalArith solve q in closed form and overwrite index_q
  # (ceattle.cpp section 8.2), so index_log_q -- and therefore Catchability_init
  # -- is never read. Several GOA hake configurations leave it at 0 and fit; the
  # check has to let them through.
  for (form in c("Analytical", "AnalyticalArith")) {
    fit <- .fit(.q_fixture(function(d, i) {
      d$fleet_control$Catchability[i] <- form
      d$fleet_control$Catchability_init[i] <- 0
      d }))
    testthat::expect_true(is.finite(fit$obj$fn()), info = form)
  }

  # Estimate_index_sd = 2 derives the sd analytically (Walters-Ludwig 1994)
  # instead of estimating it, so it reads no starting value either.
  fit <- .fit(.q_fixture(function(d, i) {
    d$fleet_control$Estimate_index_sd[i] <- 2
    d$fleet_control$Index_sd[i] <- NA
    d }))
  testthat::expect_true(is.finite(fit$obj$fn()))

  # ...and each switch fits once its own column is supplied.
  ok <- list(
    function(d, i) { d$fleet_control$Estimate_index_sd[i] <- 1
                     d$fleet_control$Index_sd[i] <- 0.2; d },
    function(d, i) { d$fleet_control$Catchability[i] <- "Estimated-with-prior"
                     d$fleet_control$Catchability_prior_sd[i] <- 0.5; d },
    function(d, i) { d$fleet_control$Time_varying_q[i] <- "IID"
                     d$fleet_control$Time_varying_q_sd[i] <- 0.2; d })
  for (k in seq_along(ok)) {
    testthat::expect_true(is.finite(.fit(.q_fixture(ok[[k]]))$obj$fn()), info = k)
  }
})

testthat::test_that("a missing Ceq is named rather than crashing the check", {
  # A workbook with the bioenergetics columns left blank used to reach
  # `if (Ceq[sp] > 1)` with NA, and data_check() died on R's bare "missing
  # value where TRUE/FALSE needed" -- which names neither the column nor the
  # species. GOA hake's input_nodiet.xlsx is the case that found this.
  d <- make_test_data(nyrs = 12, nages = 5, seed = 42)
  d$Ceq[1] <- NA
  testthat::expect_error(Rceattle:::data_check(Rceattle:::switch_check(Rceattle::clean_data(d))),
                         "Ceq")
})
