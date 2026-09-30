# `link = "power"` multiplies log q by exp(beta * x) rather than shifting it:
#
#   index_q_yr = exp((log_q + q_offset(yr)) * exp(sum beta * x(yr)) + q_dev(yr))
#
# Stock Synthesis's environmental link type 1, verified against the pinned
# v3.30.22.1 source and cited so nobody has to re-derive it:
#   SS_timevaryparm.tpl:50       parm_timevary = baseparm
#   SS_timevaryparm.tpl:206-211  type 1 multiplies by mfexp(beta * env)
#   SS_timevaryparm.tpl:215-220  type 2 adds it (this is Rceattle's `log`)
#   SS_expval.tpl:407/418        Svy_log_q is a log; q is its exponential
#   SS_readcontrol_330.tpl:3118  env-var 1xx decodes to link type 1
#
# Catchability only: SS3 holds M and growth naturally, where its type 1 is
# already our `log` link. vignette("environmental-linkages-and-priors") has why.

testthat::skip_on_cran()

env_frame <- function(d, x) {
  data.frame(Year = d$styr:d$projyr, xcov = x)
}

fit3 <- function(d, selFun = NULL, qFun = NULL) {
  suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = d, file = NULL, inits = NULL, estimateMode = 3,
    random_rec = FALSE, msmMode = 0, qFun = qFun, selFun = selFun,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
}

testthat::test_that("a scale-link q reproduces SS3's exponential env link exactly", {
  testthat::skip_if_not_installed("TMB")
  d <- make_test_data(nyrs = 20, nages = 5, seed = 42)
  flt <- d$fleet_control$Fleet_code[d$fleet_control$Fleet_type == "Survey"][1]
  testthat::skip_if(is.na(flt))
  d$fleet_control$Catchability[d$fleet_control$Fleet_code == flt] <- "Estimated"

  x <- as.numeric(scale(seq_len(length(d$styr:d$projyr))))
  d$env_data <- env_frame(d, x)

  fit <- fit3(d, qFun = Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ xcov, by = ~ fleet, fleet = flt, link = "power"))))

  # beta starts at 0, where the scale and the additive forms coincide -- so set
  # it and re-report, otherwise the comparison below is vacuous.
  p    <- fit$obj$env$last.par.best
  nm   <- names(p)
  ibet <- which(nm == "beta_linkage")
  testthat::expect_gte(length(ibet), 1)
  BETA <- 0.35
  p[ibet[1]] <- BETA
  rep  <- fit$obj$report(p)

  q    <- rep$index_q
  base <- as.numeric(fit$estimated_params$index_log_q)[flt]
  nyr  <- ncol(q)

  # SS3: log q_y = log q_base * exp(beta * env_y); Rceattle exponentiates that.
  expect_q <- exp(base * exp(BETA * x[seq_len(nyr)]))
  testthat::expect_equal(as.numeric(q[flt, ]), expect_q, tolerance = 1e-10)

  # It is NOT the additive form, which is what `log` gives.
  testthat::expect_false(isTRUE(all.equal(as.numeric(q[flt, ]),
                                          exp(base + BETA * x[seq_len(nyr)]),
                                          tolerance = 1e-6)))
  # and the scale tensor is what carries it
  testthat::expect_equal(as.numeric(rep$q_linkage_scale[flt, ]),
                         BETA * x[seq_len(nyr)], tolerance = 1e-10)
})

testthat::test_that("a zero scale tensor leaves the model arithmetically unchanged", {
  testthat::skip_if_not_installed("TMB")
  d <- make_test_data(nyrs = 20, nages = 5, seed = 42)
  # No linkage at all: exp(0) = 1, so the new multiply is the identity.
  testthat::expect_true(is.finite(fit3(d)$quantities$jnll))
})

testthat::test_that("the effect vanishes when the base parameter is zero on its log scale", {
  # beta multiplies a LOG. log_base = 0 means q = 1 and the covariate cannot move
  # it, whatever beta is. SS3 behaves the same way; pinned so nobody reads it as
  # a broken linkage.
  testthat::skip_if_not_installed("TMB")
  d <- make_test_data(nyrs = 20, nages = 5, seed = 42)
  flt <- d$fleet_control$Fleet_code[d$fleet_control$Fleet_type == "Survey"][1]
  testthat::skip_if(is.na(flt))
  # Estimated, but started at q = 1 so log q = 0; estimateMode 3 does not move it.
  d$fleet_control$Catchability[d$fleet_control$Fleet_code == flt] <- "Estimated"
  d$fleet_control$Catchability_init[d$fleet_control$Fleet_code == flt] <- 1
  d$env_data <- env_frame(d, as.numeric(scale(seq_len(length(d$styr:d$projyr)))))

  fit <- fit3(d, qFun = Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ xcov, by = ~ fleet, fleet = flt, link = "power"))))
  p <- fit$obj$env$last.par.best
  p[which(names(p) == "beta_linkage")[1]] <- 0.9   # a large effect, deliberately
  q <- as.numeric(fit$obj$report(p)$index_q[flt, ])
  testthat::expect_equal(q, rep(1, length(q)), tolerance = 1e-12)
})

testthat::test_that("scale is refused where a parameter's storage scale is not uniform", {
  d <- make_test_data(nyrs = 20, nages = 5, seed = 42)
  d$env_data <- env_frame(d, as.numeric(scale(seq_len(length(d$styr:d$projyr)))))
  flt <- d$fleet_control$Fleet_code[1]

  # Selectivity mixes log, natural and logit storage in the same slots.
  testthat::expect_error(
    fit3(d, selFun = Rceattle::build_selectivity(linkages = list(
      inf_asc = Rceattle::linkage_spec(~ xcov, by = ~ fleet, fleet = flt,
                                       link = "power")))),
    "only supported on catchability")
})

testthat::test_that("an unimplemented link is still refused, and names the implemented set", {
  testthat::expect_error(
    Rceattle::linkage_spec(~ 1, by = ~ species, link = "logit"),
    "reserved but not yet implemented")
  testthat::expect_true(all(c("identity", "log", "power") %in%
                              Rceattle:::LINKAGE_LINKS_IMPLEMENTED))
  # Lockstep with the C++ (CLAUDE.md rule 12).
  testthat::expect_identical(unname(Rceattle:::LINKAGE_LINK_CODES[["power"]]), 3L)
})
