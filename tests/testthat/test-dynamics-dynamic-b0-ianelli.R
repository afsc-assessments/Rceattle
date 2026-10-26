# Dynamic B0 is the fitted history with no fishing: the stock-recruit curve at
# the dynamic SSB times each year's realized deviation from the curve. Under the
# Ianelli form the hindcast is on mean recruitment until srr_mse_switchyr, so
# rec_dev there is measured from R0, not from the curve; the dynamic runs take
# log R - log(curve at the hindcast SSB) instead. With fishing near zero the
# dynamic SSB equals the hindcast's, so dynamic B0 equals hindcast biomass in
# every year, while numbers-at-age stay above the hindcast's 0.001 floor (the
# dynamic runs have none). Taking rec_dev directly (before 5.33.0) put the curve over a
# mean-based deviation, which on the Pacific hake operating model ran dynamic
# recruitment about 1.5 times too high in every year before the switch.

.ianelli_no_f_pars <- function(m, seed = 1) {
  set.seed(seed)
  pars <- m$obj$env$last.par
  pars[names(pars) == "log_F"] <- -30
  rid <- which(names(pars) == "rec_dev")
  testthat::expect_gt(length(rid), 0)
  pars[rid] <- stats::rnorm(length(rid), 0, 0.3)
  pars
}

test_that("under the Ianelli form dynamic B0 follows the realized deviation from the curve", {
  testthat::skip_on_cran()
  d  <- make_test_data()
  sw <- d$endyr - 5
  for (curve in c("BevertonHolt", "Ricker")) {
    m <- suppressMessages(suppressWarnings(fit_mod(
      data_list = d, inits = NULL, estimateMode = 3, random_rec = FALSE,
      recFun = build_srr(srr_fun = "mean", srr_pred_fun = curve,
                         srr_mse_switchyr = sw),
      fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE))))
    pars <- .ianelli_no_f_pars(m)
    r    <- m$obj$report(pars)
    nh   <- length(d$styr:d$endyr)
    expect_lt(max(abs(r$DynamicB0[1, 1:nh] / r$biomass[1, 1:nh] - 1)), 1e-8, label = curve)

    # Not vacuous: before the switch the curve is not R0, so rec_dev is not the
    # deviation from the curve there.
    pre <- which((d$styr:d$endyr) <= sw)[-1]
    expect_gt(max(abs(r$R_hat[1, pre] / r$R0[1] - 1)), 1e-3, label = curve)
  }
})

test_that("the same holds under predation, with the fitted suitability", {
  testthat::skip_on_cran()
  set.seed(1)
  d  <- make_msm_test_data()$data_list
  sw <- d$endyr - 5
  m  <- suppressMessages(suppressWarnings(fit_mod(
    data_list = d, inits = NULL, estimateMode = 3, msmMode = 1, suitMode = 0,
    initMode = "NonEquilibrium", random_rec = FALSE,
    recFun = build_srr(srr_fun = "mean", srr_pred_fun = "BevertonHolt",
                       srr_est_mode = "Estimated", srr_mse_switchyr = sw),
    fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE))))
  pars <- .ianelli_no_f_pars(m)
  r    <- m$obj$report(pars)
  nh   <- length(d$styr:d$endyr)

  # The hindcast floors numbers-at-age at 0.001 and the dynamic runs do not, so
  # compare only species-years whose every age has stayed well above it. At these
  # start values the prey reaches the floor in year 2, and the predator only after
  # the switch year, so the predator's pre-switch years -- where the change acts --
  # are all compared.
  pre   <- which((d$styr:d$endyr) <= sw)[-1]
  min_n <- apply(r$N_at_age[, , , 1:nh, drop = FALSE], c(1, 4), min)
  clear <- t(apply(min_n > 1, 1, cumprod)) == 1
  expect_true(all(clear[1, pre]))
  rel <- r$DynamicB0[, 1:nh] / r$biomass[, 1:nh] - 1
  expect_lt(max(abs(rel[clear])), 1e-8)

  expect_gt(max(abs(r$R_hat[1, pre] / exp(m$obj$env$parList(pars)$rec_pars[1, 1]) - 1)), 1e-3)
})

test_that("under predation and fishing, dynamic B0 is the biomass of a run with no fishing", {
  testthat::skip_on_cran()
  # Only a fished run shows dynamic B0 using unfished predators, the term pinned
  # here. The identity holds under LognormalWeight suitability (suitMode 4),
  # where parametric suitability does not read abundance; under Empirical
  # (suitMode 0, the model_config() default, so NOT the common configuration) it
  # does not, and that is by design -- dynamic B0 holds suitability at its
  # FITTED value, since suitability is a property of the predator-prey
  # relationship and must not change because of fishing. See "Deliberately not
  # changed" in ../Rceattle-dev/CLEANUP_BACKLOG.md before changing predation.hpp.
  set.seed(1)
  d  <- make_msm_test_data()$data_list
  nh <- d$endyr - d$styr + 1
  bld <- function(inits = NULL) suppressMessages(suppressWarnings(fit_mod(
    data_list = d, inits = inits, estimateMode = 3, msmMode = 1, suitMode = 4,
    initMode = "NonEquilibrium", random_rec = FALSE,
    fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE))))
  # The simulator's mean recruitment (100 and 1,000) keeps both runs off the
  # hindcast's 0.001 numbers-at-age floor, which dynamic B0 does not have.
  p <- bld()$estimated_params
  p$rec_pars[, 1] <- log(c(100, 1000))
  fishery <- d$fleet_control[["Fleet_type"]] == 1
  p$log_F[fishery, ] <- log(0.3)
  fished <- bld(p)$quantities
  p$log_F[fishery, ] <- -30
  unfished <- bld(p)$quantities

  expect_gt(min(fished$N_at_age[, 1, , 1:nh], unfished$N_at_age[, 1, , 1:nh]), 0.01)
  expect_equal(unname(fished$DynamicB0[, 1:nh]), unname(unfished$biomass[, 1:nh]),
               tolerance = 1e-8)
  # Not vacuous: fishing lowers biomass, and changes M2 by more than 0.1 somewhere.
  expect_lt(max(fished$biomass[, nh] / fished$DynamicB0[, nh]), 0.5)
  expect_gt(max(unfished$M2_at_age[, 1, , 1:nh] - fished$M2_at_age[, 1, , 1:nh]), 0.1)
})

test_that("with F near zero, dynamic SB0 equals hindcast SSB for every recruitment form", {
  # Fast, and run on CRAN: the bundled Bering Sea data, built but not optimized, with
  # mean recruitment and a curve fitted in the hindcast as controls.
  d   <- Rceattle::BS2017SS
  nh  <- d$endyr - d$styr + 1
  bld <- function(rf, inits = NULL) suppressMessages(suppressWarnings(fit_mod(
    data_list = d, inits = inits, estimateMode = 3, msmMode = 0, recFun = rf,
    fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE))))
  forms <- list(
    ianelli = build_srr(srr_fun = "mean", srr_pred_fun = "BevertonHolt"),
    curve   = build_srr(srr_fun = "BevertonHolt"),
    mean    = build_srr(srr_fun = "mean"))
  for (nm in names(forms)) {
    p <- bld(forms[[nm]])$estimated_params
    p$log_F[] <- -30
    set.seed(1)
    p$rec_dev[, 1:nh] <- stats::rnorm(length(p$rec_dev[, 1:nh]), 0, 0.3)
    q <- bld(forms[[nm]], p)$quantities
    expect_true(all(is.finite(q$ssb[, 1:nh])) && all(q$ssb[, 1:nh] > 0), info = nm)
    # Not vacuous: at these starting values the Ianelli curve is far from mean recruitment.
    if (nm == "ianelli") expect_gt(max(abs(q$R_hat[, 2:nh] / q$R0[, 1] - 1)), 1e-3)
    # Compare only species-years whose every real age (N_at_age is padded to the oldest
    # species) has stayed at 10x the hindcast's 0.001 floor, which the dynamic runs lack.
    min_n <- t(vapply(seq_len(d$nspp), function(sp) apply(
      q$N_at_age[sp, seq_len(d$nsex[sp]), seq_len(d$nages[sp]), 1:nh, drop = FALSE], 4, min),
      numeric(nh)))
    clear <- t(apply(min_n > 0.01, 1, cumprod)) == 1
    expect_true(any(rowSums(clear) == nh), label = nm)   # at least one whole species compared
    rel <- q$DynamicSB0[, 1:nh] / q$ssb[, 1:nh] - 1
    expect_lt(max(abs(rel[clear])), 1e-8, label = nm)
  }
})
