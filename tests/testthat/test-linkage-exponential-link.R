# `link = "exponential"` MULTIPLIES log q by exp(beta * x) rather than shifting it:
#
#   index_q_yr = exp((log_q + q_offset(yr)) * exp(sum beta * x(yr)) + q_dev(yr))
#
# Stock Synthesis's environmental link type 1, verified against the pinned
# v3.30.22.1 source and cited so nobody has to re-derive it:
#   SS_timevaryparm.tpl:50       parm_timevary = baseparm
#   SS_timevaryparm.tpl:206-211  `case 1: // exponential env link`, *= mfexp(beta * env)
#   SS_timevaryparm.tpl:215-220  type 2 adds it (this is Rceattle's `log`)
#   SS_expval.tpl:408            Svy_log_q = parm_timevary
#   SS_expval.tpl:413-419        q is its exponential ONLY for errtype >= 0
#                                (lognormal/t); errtype -1 (normal) reads it
#                                arithmetically, where type 1 is our `log`
#   SS_readcontrol_330.tpl:3289  env-var 1xx decodes to link type 1, for Q_parm
#
# Catchability only, and only where the base is estimable: SS3 holds M and growth
# naturally (SS_biofxn.tpl:1063, :265-275), where its type 1 is already our
# `log`. It does store the recruitment level as a log (`SR_LN(R0)`), so the form
# applies there in principle, but no accumulator consumes it yet.
# vignette("environmental-linkages-and-priors") has the rest.

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

testthat::test_that("an exponential-link q reproduces SS3's exponential env link exactly", {
  testthat::skip_if_not_installed("TMB")
  d <- make_test_data(nyrs = 20, nages = 5, seed = 42)
  flt <- d$fleet_control$Fleet_code[d$fleet_control$Fleet_type == "Survey"][1]
  testthat::skip_if(is.na(flt))
  d$fleet_control$Catchability[d$fleet_control$Fleet_code == flt] <- "Estimated"

  x <- as.numeric(scale(seq_len(length(d$styr:d$projyr))))
  d$env_data <- env_frame(d, x)

  fit <- fit3(d, qFun = Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ xcov, by = ~ fleet, fleet = flt, link = "exponential"))))

  # beta starts at 0, where the exponential and the additive forms coincide -- so set
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
  # and the log-multiplier tensor is what carries it
  testthat::expect_equal(as.numeric(rep$q_linkage_log_mult[flt, ]),
                         BETA * x[seq_len(nyr)], tolerance = 1e-10)
})

testthat::test_that("a zero multiplier tensor leaves the model arithmetically unchanged", {
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
    q = Rceattle::linkage_spec(~ xcov, by = ~ fleet, fleet = flt, link = "exponential"))))
  p <- fit$obj$env$last.par.best
  p[which(names(p) == "beta_linkage")[1]] <- 0.9   # a large effect, deliberately
  q <- as.numeric(fit$obj$report(p)$index_q[flt, ])
  testthat::expect_equal(q, rep(1, length(q)), tolerance = 1e-12)
})

testthat::test_that("an unimplemented link is still refused, and names the implemented set", {
  testthat::expect_error(
    Rceattle::linkage_spec(~ 1, by = ~ species, link = "logit"),
    "reserved but not yet implemented")
  testthat::expect_true(all(c("identity", "log", "exponential") %in%
                              Rceattle:::LINKAGE_LINKS_IMPLEMENTED))
  # Lockstep with the C++ (CLAUDE.md rule 12).
  testthat::expect_identical(unname(Rceattle:::LINKAGE_LINK_CODES[["exponential"]]), 3L)
})


testthat::test_that("a shared q linkage row is checked on every fleet", {
  # NA fleet is the shared sentinel and the cpp expands it to all fleets, so a
  # shared row must be checked against fleets that have no catchability at all.
  fc <- Rceattle::switch_check(Rceattle::clean_data(Rceattle::BS2017SS))$fleet_control
  tbl <- Rceattle:::bind_linkage(
    Rceattle:::linkage_row(process = "q", param = "q", X_col = 1L, fleet = 7L,
                           design_col = "xcov", link = "log"),
    Rceattle:::linkage_row(process = "q", param = "q", X_col = 2L,
                           design_col = "PDO", link = "log"))
  testthat::expect_error(Rceattle:::.check_q_linkage_support(tbl, fc),
                         "does not estimate q")
})


testthat::test_that("exponential is refused on every process but catchability", {
  # The correction an adversarial review forced. An earlier form of this link ran
  # on recruitment, M and growth on the reasoning that Rceattle stores all four
  # as logs; what matters is the scale SS3 stores, and on M -- where log M < 0 --
  # the form inverted the sign of the covariate effect.
  env <- data.frame(Year = 1980:1999,
                    xcov = as.numeric(scale(seq_len(20))))
  for (proc in c("recruitment", "M", "growth", "sel")) {
    prm <- switch(proc, recruitment = "R0", M = "M1", growth = "Linf",
                  sel = "inf_asc")
    testthat::expect_error(
      Rceattle:::materialize_linkage(
        Rceattle::linkage_spec(~ xcov, param = prm, link = "exponential"),
        proc, env),
      "only consumed by catchability",
      info = proc)
  }
  # and the message does not hand a recruitment user `log`, which is SS3's type 2
  err <- tryCatch(
    Rceattle:::materialize_linkage(
      Rceattle::linkage_spec(~ xcov, param = "R0", link = "exponential"),
      "recruitment", env),
    error = conditionMessage)
  testthat::expect_match(err, "SR_LN\\(R0\\)")
  testthat::expect_match(err, "NOT a substitute")
})


testthat::test_that("exponential is refused on a random-effect row", {
  # A deviation multiplied by exp(beta * x) has effective sd sigma * exp(beta * x)
  # while its density still scores it at a constant sigma -- the same defect the
  # consume site avoids for index_q_dev by keeping it outside the multiply.
  testthat::expect_error(
    Rceattle:::linkage_row(process = "q", param = "q", X_col = 1L, fleet = 1L,
                           link = "exponential", re_struct = "ar1",
                           re_group = "Year"),
    "cannot carry a random effect")
})


testthat::test_that("exponential and a random-effect q on one fleet are refused", {
  # Across specs: every log-link row accumulates into q_linkage_offset, which
  # sits inside the multiply, so the deviations would be scaled too.
  ex <- Rceattle:::linkage_row(process = "q", param = "q", X_col = 1L,
                               fleet = 3L, link = "exponential")
  re <- Rceattle:::linkage_row(process = "q", param = "q", X_col = 2L,
                               fleet = 3L, link = "log",
                               re_struct = "ar1", re_group = "Year")
  testthat::expect_error(Rceattle:::bind_linkage(ex, re),
                         "random-effect catchability linkage")
  # a shared (NA fleet) row reaches every fleet, so it clashes too
  re_all <- Rceattle:::linkage_row(process = "q", param = "q", X_col = 2L,
                                   link = "log", re_struct = "ar1",
                                   re_group = "Year")
  testthat::expect_error(Rceattle:::bind_linkage(ex, re_all),
                         "random-effect catchability linkage")
})


testthat::test_that("exponential warns on a natural-scale index family", {
  # SS3 stores q as a log only for a lognormal/t survey (SS_expval.tpl:413-419);
  # under MVN / Normal / TruncatedNormal it reads the same slot arithmetically,
  # where its type 1 is q * exp(beta * x) -- our `log` link. A warning rather
  # than a refusal: Rceattle holds q on the log scale whatever the index family,
  # so the model is well defined, it is just the wrong bridge.
  fc <- Rceattle::switch_check(Rceattle::clean_data(Rceattle::BS2017SS))$fleet_control
  # Intercept-bearing, so only the index-family check can fire here.
  tbl <- Rceattle:::bind_linkage(
    Rceattle:::linkage_row(process = "q", param = "q", X_col = 1L, fleet = 7L,
                           design_col = "(Intercept)", link = "exponential"),
    Rceattle:::linkage_row(process = "q", param = "q", X_col = 2L, fleet = 7L,
                           design_col = "xcov", link = "exponential"))
  fc$Index_distribution[fc$Fleet_code == 7L] <- "Normal"
  testthat::expect_warning(
    Rceattle:::.check_q_linkage_support(tbl, fc),
    "NATURAL scale")
  fc$Index_distribution[fc$Fleet_code == 7L] <- "Lognormal"
  testthat::expect_silent(Rceattle:::.check_q_linkage_support(tbl, fc))
})


testthat::test_that("exponential is refused where index_log_q is mapped out", {
  # It multiplies log q, so it needs a free base. A slope-only formula, or an
  # intercept fixed at 0, makes map_linkage_adjuster() mask index_log_q: log q is
  # then frozen at log(Catchability_init), and at 1 that is exactly 0, where beta
  # has an identically zero gradient for every value it could take.
  fc <- Rceattle::switch_check(Rceattle::clean_data(Rceattle::BS2017SS))$fleet_control
  slope_only <- Rceattle:::linkage_row(process = "q", param = "q", X_col = 1L,
                                       fleet = 7L, design_col = "xcov",
                                       link = "exponential")
  testthat::expect_error(Rceattle:::.check_q_linkage_support(slope_only, fc),
                         "needs an estimated base catchability")

  fixed_icept <- Rceattle:::bind_linkage(
    Rceattle:::linkage_row(process = "q", param = "q", X_col = 1L, fleet = 7L,
                           design_col = "(Intercept)", link = "exponential",
                           est_phase = 0L),
    Rceattle:::linkage_row(process = "q", param = "q", X_col = 2L, fleet = 7L,
                           design_col = "xcov", link = "exponential"))
  testthat::expect_error(Rceattle:::.check_q_linkage_support(fixed_icept, fc),
                         "needs an estimated base catchability")

  ok <- Rceattle:::bind_linkage(
    Rceattle:::linkage_row(process = "q", param = "q", X_col = 1L, fleet = 7L,
                           design_col = "(Intercept)", link = "exponential"),
    Rceattle:::linkage_row(process = "q", param = "q", X_col = 2L, fleet = 7L,
                           design_col = "xcov", link = "exponential"))
  testthat::expect_error(Rceattle:::.check_q_linkage_support(ok, fc), NA)
})


testthat::test_that("a q starting at exactly 1 warns that beta has no gradient", {
  fc <- Rceattle::switch_check(Rceattle::clean_data(Rceattle::BS2017SS))$fleet_control
  fc$Catchability_init[fc$Fleet_code == 7L] <- 1
  tbl <- Rceattle:::bind_linkage(
    Rceattle:::linkage_row(process = "q", param = "q", X_col = 1L, fleet = 7L,
                           design_col = "(Intercept)", link = "exponential"),
    Rceattle:::linkage_row(process = "q", param = "q", X_col = 2L, fleet = 7L,
                           design_col = "xcov", link = "exponential"))
  testthat::expect_warning(Rceattle:::.check_q_linkage_support(tbl, fc),
                           "Catchability_init is 1")
})
