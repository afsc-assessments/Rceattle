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
#   SS_readcontrol_330.tpl:3328  env-var 1xx decodes to link type 1, for Q_parm
#
# Catchability only, and only where the base is estimable: SS3 holds M and growth
# naturally (SS_biofxn.tpl:1074, :265-275), where its type 1 is already our
# `log`. It does store the recruitment level as a log (`SR_LN(R0)`), so the form
# applies there in principle, but no accumulator consumes it yet.
# vignette("environmental-linkages-and-priors") has the rest.

testthat::skip_on_cran()

.exp_env <- function(d) {
  yrs <- d$styr:d$projyr
  data.frame(Year = yrs, xcov = as.numeric(scale(seq_along(yrs))))
}

.exp_fit3 <- function(d, qFun = NULL) {
  suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = d, file = NULL, inits = NULL, estimateMode = 3,
    random_rec = FALSE, msmMode = 0, qFun = qFun,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0))))
}

# A survey q of 0.3: away from 1, so `log q` is a real number the covariate can
# scale, and negative, which is where the sign of the effect inverts.
.exp_data <- function(q_init = 0.3, nyrs = 20) {
  d <- make_test_data(nyrs = nyrs, nages = 5, seed = 42)
  flt <- d$fleet_control$Fleet_code[d$fleet_control$Fleet_type == "Survey"][1]
  d$fleet_control$Catchability[d$fleet_control$Fleet_code == flt] <- "Estimated"
  d$fleet_control$Catchability_init[d$fleet_control$Fleet_code == flt] <- q_init
  d$env_data <- .exp_env(d)
  list(d = d, flt = flt)
}

# Pick a coefficient by its design column. `beta_linkage` is aligned row-for-row
# with the linkage table, but TMB drops map = NA entries from last.par.best, so
# the parameter vector holds only the ESTIMATED rows -- an intercept-bearing q
# formula contributes one entry, not two. Subset the table by the map first.
.exp_beta_index <- function(fit, col) {
  est <- !is.na(fit$map$mapList$beta_linkage)
  which(names(fit$obj$env$last.par.best) == "beta_linkage")[
    which(fit$data_list$linkage_table$design_col[est] == col)]
}


testthat::test_that("an exponential q reproduces SS3's env link type 1 exactly", {
  testthat::skip_if_not_installed("TMB")
  s <- .exp_data()
  fit <- .exp_fit3(s$d, qFun = Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ xcov, by = ~ fleet, fleet = s$flt,
                               link = "exponential"))))

  # beta starts at 0, where every link coincides, so drive it and re-report.
  slope <- .exp_beta_index(fit, "xcov")
  testthat::expect_length(slope, 1L)
  BETA <- 0.35
  p <- fit$obj$env$last.par.best
  p[slope] <- BETA
  rep <- fit$obj$report(p)

  q    <- rep$index_q
  base <- as.numeric(fit$estimated_params$index_log_q)[s$flt]
  x    <- s$d$env_data$xcov[seq_len(ncol(q))]
  # The base has to be a real number, or the comparison below holds trivially.
  testthat::expect_lt(base, -1)

  # SS3: log q_y = log q_base * exp(beta * env_y); Rceattle exponentiates that.
  testthat::expect_equal(as.numeric(q[s$flt, ]), exp(base * exp(BETA * x)),
                         tolerance = 1e-10)
  # It is NOT the additive form, which is what `log` gives.
  testthat::expect_false(isTRUE(all.equal(as.numeric(q[s$flt, ]),
                                          exp(base + BETA * x),
                                          tolerance = 1e-6)))
  # and q genuinely moves, so the assertion above is not comparing constants
  testthat::expect_gt(stats::sd(as.numeric(q[s$flt, ])), 1e-6)
  # the log-multiplier tensor is what carries it
  testthat::expect_equal(as.numeric(rep$q_linkage_log_mult[s$flt, ]), BETA * x,
                         tolerance = 1e-10)

  # and the identifiability record runs on a real fit: q here is 0.3, well away
  # from 1, so it must stay silent rather than merely not erroring.
  testthat::expect_false("exponential_q_near_one" %in%
                           names(fit$convergence$checks))
})


testthat::test_that("the covariate effect inverts in sign below q = 1", {
  # SS3's property, not ours: beta multiplies a LOG, so for q < 1 (log q < 0) a
  # positive beta DECREASES q. Pinned so nobody reads it as a broken linkage.
  testthat::skip_if_not_installed("TMB")
  s <- .exp_data(q_init = 0.3)
  fit <- .exp_fit3(s$d, qFun = Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ xcov, by = ~ fleet, fleet = s$flt,
                               link = "exponential"))))
  p <- fit$obj$env$last.par.best
  p[.exp_beta_index(fit, "xcov")] <- 0.5
  q <- as.numeric(fit$obj$report(p)$index_q[s$flt, ])
  x <- s$d$env_data$xcov[seq_len(length(q))]
  testthat::expect_lt(stats::cor(q, x), -0.9)   # q falls as the covariate rises
  testthat::expect_true(all(q < 1))             # and stays below 1 throughout
})


testthat::test_that("a zero multiplier tensor leaves the model unchanged", {
  testthat::skip_if_not_installed("TMB")
  d <- make_test_data(nyrs = 20, nages = 5, seed = 42)
  # No linkage at all: exp(0) = 1, so the new multiply is the identity.
  testthat::expect_true(is.finite(.exp_fit3(d)$quantities$jnll))
})


testthat::test_that("the identified quantity is beta * log q, not beta", {
  # Simulate an index from a known beta under q^exp(beta * x) and refit. What
  # comes back is NOT beta: with a free base the fit is happy to cross q = 1,
  # where log q changes sign, and re-express the same curve with the opposite
  # sign of beta. Measured here: a truth of (q = 0.3, beta = +0.4) refits to
  # (q = 467.0, beta = -0.0775), preserving log(q) * beta to 1.1% -- to first
  # order log q_y = log q * (1 + beta * x), so log q * beta is the slope on log q
  # and is what the index informs. The assertions below pin that invariant, not
  # those values, which are one platform's optimizer path. Report beta with its
  # base, never alone.
  testthat::skip_if_not_installed("TMB")
  s <- .exp_data(q_init = 0.3, nyrs = 60)
  s$d$fleet_control$Index_sd[s$d$fleet_control$Fleet_code == s$flt] <- 0.1
  qfun <- Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ xcov, by = ~ fleet, fleet = s$flt,
                               link = "exponential")))
  fit <- .exp_fit3(s$d, qFun = qfun)

  BETA <- 0.4
  p <- fit$obj$env$last.par.best
  p[.exp_beta_index(fit, "xcov")] <- BETA
  q_true <- as.numeric(fit$obj$report(p)$index_q[s$flt, ])
  L_true <- as.numeric(fit$estimated_params$index_log_q)[s$flt]
  set.seed(11)
  sim <- fit$obj$simulate(p)
  i <- which(s$d$index_data$Fleet_code == s$flt)
  dsim <- s$d
  dsim$index_data$Observation[i] <- as.numeric(sim$index_hat)[i]

  ref <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = dsim, file = NULL, inits = NULL, estimateMode = 1,
    random_rec = FALSE, msmMode = 0, qFun = qfun,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0))))
  tbl  <- ref$data_list$linkage_table
  bhat <- as.numeric(ref$estimated_params$beta_linkage)[
    tbl$design_col == "xcov"]
  Lhat <- as.numeric(ref$estimated_params$index_log_q)[s$flt]
  testthat::expect_length(bhat, 1L)

  # The q series is recovered in shape (its level trades off against biomass, as
  # any catchability does).
  q_hat <- as.numeric(ref$quantities$index_q[s$flt, seq_along(q_true)])
  testthat::expect_gt(stats::cor(log(q_hat), log(q_true)), 0.95)

  # And the first-order effect on log q is recovered, while beta alone is not.
  testthat::expect_equal(Lhat * bhat, L_true * BETA, tolerance = 0.15)
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
        proc, env, strata = list(species = 1L, sex = 1L, fleet = 1L)),
      "only consumed by catchability",
      info = proc)
  }
  # and the message does not hand a recruitment user `log`, which is SS3's type 2
  err <- tryCatch(
    Rceattle:::materialize_linkage(
      Rceattle::linkage_spec(~ xcov, param = "R0", link = "exponential"),
      "recruitment", env, strata = list(species = 1L)),
    error = conditionMessage)
  testthat::expect_match(err, "SR_LN\\(R0\\)")
  testthat::expect_match(err, "rather than a substitute")
})


testthat::test_that("exponential is refused on a random-effect row", {
  # A deviation multiplied by exp(beta * x) has effective sd sigma * exp(beta * x)
  # while its density still scores it at a constant sigma -- the same defect the
  # consume site avoids for index_q_dev by keeping it outside the multiply. One
  # row that is both is caught by the same fleet-overlap check as two rows.
  testthat::expect_error(
    Rceattle:::linkage_row(process = "q", param = "q", X_col = 1L, fleet = 1L,
                           link = "exponential", re_struct = "ar1",
                           re_group = "Year"),
    "random-effect catchability linkage")
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


testthat::test_that("a q fitted at 1 is flagged on the fit, not just documented", {
  # The identifiability hazard is a silently-wrong-number shape, so it is a
  # convergence record. It reads the FITTED log q: Catchability_init is not where
  # index_log_q started once an intercept `init`, fit_mod(inits = ) or a shared
  # Catchability_index group is involved, so a start-based test false-positives.
  fc <- Rceattle::switch_check(Rceattle::clean_data(Rceattle::BS2017SS))$fleet_control
  tbl <- Rceattle:::bind_linkage(
    Rceattle:::linkage_row(process = "q", param = "q", X_col = 1L, fleet = 7L,
                           design_col = "(Intercept)", link = "exponential"),
    Rceattle:::linkage_row(process = "q", param = "q", X_col = 2L, fleet = 7L,
                           design_col = "xcov", link = "exponential"))
  mk <- function(q_mle) {
    lq <- rep(0, nrow(fc)); lq[7] <- log(q_mle)
    list(data_list = list(linkage_table = tbl, fleet_control = fc),
         estimated_params = list(index_log_q = lq))
  }
  r <- Rceattle:::.check_exponential_q_sign(mk(1 + 1e-6))
  testthat::expect_true("exponential_q_near_one" %in% names(r))
  testthat::expect_equal(r$exponential_q_near_one$severity, "WARN")
  testthat::expect_match(r$exponential_q_near_one$message, "within 0.1% of 1",
                         fixed = TRUE)
  # a q well away from 1 has nothing to say, whatever it started at
  testthat::expect_length(Rceattle:::.check_exponential_q_sign(mk(0.25)), 0L)
  # and a model with no exponential linkage is untouched
  testthat::expect_length(
    Rceattle:::.check_exponential_q_sign(list(data_list = list())), 0L)
})


testthat::test_that("an unimplemented link is still refused, and the codes are in lockstep", {
  testthat::expect_error(
    Rceattle::linkage_spec(~ 1, by = ~ species, link = "logit"),
    "reserved but not yet implemented")
  testthat::expect_true(all(c("identity", "log", "exponential") %in%
                              Rceattle:::LINKAGE_LINKS_IMPLEMENTED))
  # Lockstep with the C++ (CLAUDE.md rule 12).
  testthat::expect_identical(
    unname(Rceattle:::LINKAGE_LINK_CODES[["exponential"]]), 3L)
})
