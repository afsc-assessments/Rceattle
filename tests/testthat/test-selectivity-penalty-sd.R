# =============================================================================
# Non-parametric selectivity penalties can be specified as standard deviations
# instead of the cryptic penalty WEIGHTS. Every such penalty is a Gaussian SSQ
# (weight * x^2 = x^2 / (2*sd^2)), so switch_check() converts the SD columns
# (Sel_shape_sd [+ Sel_shape_dir], Sel_curvature_sd, Sel_devmag_sd) into
# Sel_curve_pen1/2/3 via weight = 1/(2*sd^2). Legacy Sel_curve_pen columns are
# left untouched (so existing models are bit-identical); a fleet supplying the SD
# columns fits equivalently to the same model expressed as weights.
# =============================================================================

testthat::test_that("switch_check converts penalty SDs to Sel_curve_pen weights", {
  testthat::skip_if_not_installed("Rceattle")
  d <- Rceattle::BS2017SS
  np <- which(d$fleet_control$Selectivity == 2)
  testthat::skip_if(length(np) == 0)
  d$fleet_control$Sel_curve_pen1[np] <- NA_real_   # drop legacy so conversion fires
  d$fleet_control$Sel_curve_pen2[np] <- NA_real_
  d$fleet_control$Sel_shape_sd      <- NA_real_; d$fleet_control$Sel_shape_sd[np]      <- 1 / sqrt(2 * 20)
  d$fleet_control$Sel_shape_dir     <- NA;       d$fleet_control$Sel_shape_dir[np]     <- "Decreasing"
  d$fleet_control$Sel_curvature_sd  <- NA_real_; d$fleet_control$Sel_curvature_sd[np]  <- 1 / sqrt(2 * 12.5)

  out <- suppressMessages(Rceattle::switch_check(d))
  testthat::expect_equal(out$fleet_control$Sel_curve_pen1[np], rep(20, length(np)),  tolerance = 1e-8)
  testthat::expect_equal(out$fleet_control$Sel_curve_pen2[np], rep(12.5, length(np)), tolerance = 1e-8)

  # "Increasing" is refused on NonParametric: its shape term does not branch on
  # the sign, so a negative weight would reward a decreasing curve without bound.
  d2 <- d; d2$fleet_control$Sel_shape_dir[np] <- "Increasing"
  testthat::expect_error(suppressMessages(Rceattle::switch_check(d2)),
                         "does not read the sign")

  # NonParametricPM is the one form that does branch on it, so it still flips.
  d3 <- d; d3$fleet_control$Selectivity[np] <- "NonParametricPM"
  d3$fleet_control$Sel_shape_dir[np] <- "Increasing"
  out3 <- suppressMessages(Rceattle::switch_check(d3))
  testthat::expect_equal(out3$fleet_control$Sel_curve_pen1[np], rep(-20, length(np)), tolerance = 1e-8)
})

testthat::test_that("penalty SD columns are rejected on non-non-parametric forms", {
  testthat::skip_if_not_installed("Rceattle")
  d <- Rceattle::BS2017SS
  # fleet 4 is logistic (Selectivity == 1): a penalty SD there is meaningless and
  # (for AR1 forms) would corrupt the logit-rho reuse of Sel_curve_pen.
  logi <- which(d$fleet_control$Selectivity == 1)
  testthat::skip_if(length(logi) == 0)
  d$fleet_control$Sel_shape_sd <- NA_real_; d$fleet_control$Sel_shape_sd[logi[1]] <- 0.2
  testthat::expect_error(suppressMessages(Rceattle::switch_check(d)),
                         "NonParametric")
})

testthat::test_that("LogisticPM accepts Sel_shape_sd / Sel_devmag_sd (two-sided, positive)", {
  testthat::skip_if_not_installed("Rceattle")
  d <- Rceattle::BS2017SS
  f <- which(d$fleet_control$Selectivity == 1)[1]   # repurpose a logistic fleet
  testthat::skip_if(length(f) == 0)
  d$fleet_control$Selectivity[f] <- 11              # LogisticPM
  d$fleet_control$Sel_shape_sd  <- NA_real_; d$fleet_control$Sel_shape_sd[f]  <- 1 / sqrt(2 * 20)
  d$fleet_control$Sel_devmag_sd <- NA_real_; d$fleet_control$Sel_devmag_sd[f] <- 1 / sqrt(2 * 8)
  out <- suppressMessages(Rceattle::switch_check(d))
  testthat::expect_equal(out$fleet_control$Sel_curve_pen1[f], 20, tolerance = 1e-8)  # pen1 (shape)
  testthat::expect_equal(out$fleet_control$Sel_curve_pen3[f], 8,  tolerance = 1e-8)  # pen3 (dev-mag)
})

testthat::test_that("Sel_curvature_sd is rejected on LogisticPM (pen2 unused there)", {
  testthat::skip_if_not_installed("Rceattle")
  d <- Rceattle::BS2017SS
  f <- which(d$fleet_control$Selectivity == 1)[1]
  testthat::skip_if(length(f) == 0)
  d$fleet_control$Selectivity[f] <- 11
  d$fleet_control$Sel_curvature_sd <- NA_real_; d$fleet_control$Sel_curvature_sd[f] <- 0.2
  testthat::expect_error(suppressMessages(Rceattle::switch_check(d)), "NonParametric")
})

testthat::test_that("Sel_shape_dir = 'Increasing' is rejected on LogisticPM", {
  testthat::skip_if_not_installed("Rceattle")
  d <- Rceattle::BS2017SS
  f <- which(d$fleet_control$Selectivity == 1)[1]
  testthat::skip_if(length(f) == 0)
  d$fleet_control$Selectivity[f] <- 11
  d$fleet_control$Sel_shape_sd  <- NA_real_; d$fleet_control$Sel_shape_sd[f]  <- 0.2
  d$fleet_control$Sel_shape_dir <- NA;       d$fleet_control$Sel_shape_dir[f] <- "Increasing"
  testthat::expect_error(suppressMessages(Rceattle::switch_check(d)),
                         "does not read the sign")
})

testthat::test_that("a non-positive penalty SD is rejected", {
  testthat::skip_if_not_installed("Rceattle")
  d <- Rceattle::BS2017SS
  np <- which(d$fleet_control$Selectivity == 2)
  testthat::skip_if(length(np) == 0)
  d$fleet_control$Sel_curve_pen1[np] <- NA_real_
  d$fleet_control$Sel_shape_sd <- NA_real_; d$fleet_control$Sel_shape_sd[np[1]] <- 0
  testthat::expect_error(suppressMessages(Rceattle::switch_check(d)),
                         "positive standard deviation")
})

testthat::test_that("penalty SD columns fit equivalently to the legacy weights", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")
  ctl <- Rceattle::fit_control(phase = TRUE, verbose = 0)

  base <- suppressMessages(Rceattle::fit_mod(
    data_list = Rceattle::BS2017SS, estimateMode = 0, random_rec = FALSE,
    msmMode = 0, fit_control = ctl))

  d <- Rceattle::BS2017SS
  np <- which(d$fleet_control$Selectivity == 2)
  d$fleet_control$Sel_shape_sd     <- NA_real_; d$fleet_control$Sel_shape_sd[np]     <- 1 / sqrt(2 * d$fleet_control$Sel_curve_pen1[np])
  d$fleet_control$Sel_curvature_sd <- NA_real_; d$fleet_control$Sel_curvature_sd[np] <- 1 / sqrt(2 * d$fleet_control$Sel_curve_pen2[np])
  d$fleet_control$Sel_curve_pen1[np] <- NA_real_
  d$fleet_control$Sel_curve_pen2[np] <- NA_real_
  sdfit <- suppressMessages(Rceattle::fit_mod(
    data_list = d, estimateMode = 0, random_rec = FALSE, msmMode = 0, fit_control = ctl))

  # Equivalent up to the floating-point of the sd<->weight round-trip amplified
  # by the optimisation (~1e-7 on an objective ~1e4, i.e. relative ~1e-11).
  testthat::expect_equal(sdfit$opt$objective, base$opt$objective, tolerance = 1e-5)
})


testthat::test_that("mode 5 does not feed sel_dev_log_sd from unestimated deviates", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # Time_varying_sel = "RandomWalkAscending" varies the ascending limb only, so
  # build_map() estimates no descending deviate. Those deviates sit at 0, and the
  # descending random-walk penalty on them used to be accumulated anyway. With
  # random_sel = FALSE that is a constant (sel_dev_log_sd is mapped out); with
  # random_sel = TRUE the SD IS estimated, so the spurious term is
  # 2 * nyrs * nsex * log(sigma) and it biases the SD downward. Pin the gradient
  # so that term cannot come back unnoticed for either setting.
  testthat::skip_if_not(exists("GOA2018SS"))
  d <- Rceattle::GOA2018SS
  flt <- which(d$fleet_control$Time_varying_sel == 5)
  testthat::skip_if(length(flt) == 0)

  fit <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = d, file = NULL, inits = NULL, estimateMode = 3,
    random_rec = FALSE, random_sel = TRUE, msmMode = 0,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0))))

  # The SD is estimable under random_sel = TRUE -- otherwise this test is vacuous.
  testthat::expect_true("sel_dev_log_sd" %in% names(fit$obj$par))

  # Its gradient must come only from deviates the model actually estimates. The
  # descending deviates are mapped out at 0, so a descending penalty would add a
  # term with no data behind it; the ascending pair is what legitimately informs
  # the SD.
  g <- fit$obj$gr(fit$obj$par)
  gsd <- g[names(fit$obj$par) == "sel_dev_log_sd"]
  testthat::expect_true(all(is.finite(gsd)))
  testthat::expect_equal(length(gsd), sum(names(fit$obj$par) == "sel_dev_log_sd"))
})

# A Sel_curve_pen slot multiplies a squared deviation, so a negative weight
# rewards that deviation rather than penalizing it. Measured on BS2017SS fleet 1
# (form 2, N_sel_bins = 8, Sel_curve_pen2 = 12.5) by ramping sel_coff down a
# constant step per bin, the JNLL_SEL_NONPARAM row is exactly
#   [Sel_curve_pen1 * (N_sel_bins - 1) + Sel_curve_pen2] * step^2 + 8.6
# -- -29, -122, -502, -2032 at steps 0.5/1/2/4 with Sel_curve_pen1 = -20, and
# the whole objective -6.6e6 at step 256. The coefficient 1 on Sel_curve_pen2 is
# the single second difference where the ramp meets the coefficient repeated
# past N_sel_bins, so the curvature penalty DOES bind and the objective diverges
# only past |Sel_curve_pen1| > Sel_curve_pen2 / (N_sel_bins - 1) = 1.79: at
# Sel_curve_pen1 = -1 the row rises (+1417, +5641, +22537 at steps 16/32/64).
# A smaller negative weight is anti-shrinking rather than divergent, and is
# refused too. Only NonParametricPM branches on the sign; the AR1 forms reuse
# the columns as logit-scale correlations, and 13 never reads slot 3.
testthat::test_that("a negative Sel_curve_pen is refused on the forms that read it as a weight", {
  d0 <- Rceattle::BS2017SS
  np <- which(d0$fleet_control$Selectivity == 2)
  d0$fleet_control$Sel_curve_pen3 <- NA_real_

  # Distinguish "refused for the right reason" from "refused for some other
  # reason": an expect_false on a tryCatch that swallows every error would pass
  # on a missing column just as happily as on the exemption being honoured.
  run <- function(d) tryCatch({
    suppressWarnings(suppressMessages(
      Rceattle:::data_check(suppressMessages(Rceattle::switch_check(d)))))
    ""
  }, error = function(e) conditionMessage(e))
  refused <- function(d) grepl("rewards the deviation", run(d))
  allowed <- function(d) identical(run(d), "")
  # tv drives the slots charged on the DEVIATES rather than on the base curve;
  # BS2017SS ships Time_varying_sel = "Off", where those terms are identically 0.
  set_form <- function(form, col, val, tv = NULL) {
    d <- d0
    d$fleet_control$Selectivity[np]    <- form
    d$fleet_control$Sel_curve_pen1[np] <- 20
    d$fleet_control$Sel_curve_pen2[np] <- 12.5
    if (!is.null(tv)) {
      d$fleet_control$Time_varying_sel[np] <- tv
      # Estimated deviations are penalized at this sd; without it data_check()
      # stops on that instead and the case below would never be reached.
      d$fleet_control$Time_varying_sel_sd[np] <- 0.2
    }
    d$fleet_control[[col]][np] <- val
    d
  }

  # Slot 1 is charged on the BASE coefficients on 2 / 13, so it is read whatever
  # Time_varying_sel says. Allowed on 9, whose shape term reads the sign.
  testthat::expect_true(refused(set_form("NonParametric",           "Sel_curve_pen1", -20)))
  testthat::expect_true(refused(set_form("NonParametricIntegrable", "Sel_curve_pen1", -20)))
  testthat::expect_true(allowed(set_form("NonParametricPM",         "Sel_curve_pen1", -20)))

  # Slot 2: two-sided on every non-parametric form, and also on the base curve.
  testthat::expect_true(refused(set_form("NonParametric",           "Sel_curve_pen2", -5)))
  testthat::expect_true(refused(set_form("NonParametricIntegrable", "Sel_curve_pen2", -5)))

  # Slot 3 is read only by 9 and 11. Forms 2 and 13 also estimate deviates, but
  # score them with a density on Time_varying_sel_sd, so a sign check on the
  # slot they never read would be noise.
  testthat::expect_true(allowed(set_form("NonParametric",           "Sel_curve_pen3", -5)))
  testthat::expect_true(allowed(set_form("NonParametricIntegrable", "Sel_curve_pen3", -5)))

  # The three deviate-charged slots -- LogisticPM 1 and 3, NonParametricPM 3 --
  # are inert under Time_varying_sel = "Off" (objective bit-identical at +w, 0
  # and -w), so a negative weight there is allowed, and refused once the
  # deviates are estimated.
  # "RandomWalk" on both: it is the only time-varying structure LogisticPM
  # accepts, and NonParametricPM takes it too.
  testthat::expect_true(allowed(set_form("LogisticPM",      "Sel_curve_pen1", -20)))
  testthat::expect_true(allowed(set_form("NonParametricPM", "Sel_curve_pen3",  -5)))
  testthat::expect_true(refused(set_form("LogisticPM",      "Sel_curve_pen1", -20, tv = "RandomWalk")))
  testthat::expect_true(refused(set_form("NonParametricPM", "Sel_curve_pen3",  -5, tv = "RandomWalk")))

  # An "Off" fleet is charged no selectivity penalty, so its weight is not read.
  off <- set_form("NonParametric", "Sel_curve_pen1", -20)
  off$fleet_control$Fleet_type[np] <- "Off"
  testthat::expect_true(allowed(off))

  # The other half of the template's gate: flt_sel_lead == 1. Fleets sharing a
  # Selectivity_index are charged the penalty once, on the group's lead, so a
  # FOLLOWER's weight is never read and a negative one there is inert. The lead's
  # is read, so it is still refused.
  grp <- set_form("NonParametric", "Sel_curve_pen1", 20)
  two <- np[1:2]
  grp$fleet_control$Selectivity_index[two] <- grp$fleet_control$Fleet_code[two[1]]
  follower <- grp; follower$fleet_control$Sel_curve_pen1[two[2]] <- -20
  testthat::expect_true(allowed(follower))
  leader <- grp; leader$fleet_control$Sel_curve_pen1[two[1]] <- -20
  testthat::expect_true(refused(leader))

  # Positive weights are untouched.
  testthat::expect_true(allowed(set_form("NonParametric",           "Sel_curve_pen1", 20)))
})

testthat::test_that("Sel_devmag_sd is refused where Sel_curve_pen3 is not read", {
  # Neither "NonParametric" (2) nor "NonParametricIntegrable" (13) reads the
  # slot: both estimate selectivity deviates, but score them with a Gaussian
  # density on Time_varying_sel_sd. Converting an SD into it is a no-op.
  for (form in c("NonParametric", "NonParametricIntegrable")) {
    d  <- Rceattle::BS2017SS
    np <- which(d$fleet_control$Selectivity == 2)
    d$fleet_control$Selectivity[np]    <- form
    d$fleet_control$Sel_devmag_sd      <- NA_real_
    d$fleet_control$Sel_devmag_sd[np]  <- 0.5
    testthat::expect_error(suppressMessages(Rceattle::switch_check(d)),
                           "does not use that penalty slot")
  }
})

# "NonParametricPM" (9) is the only form whose shape term branches on the sign of
# Sel_curve_pen1, and only in its "Directional" mode; until 5.42.0 nothing
# checked that the increasing direction does what it says. Spec: ADMB/AMAK
# Selectivity_Likelihood sel_like(1), a one-sided SSQ over adjacent bin pairs --
# w * sum(max(+-d, 0)^2) with d = base(bin) - base(bin+1). Computable from
# sel_coff alone, the per-year mean-centring cancelling in the difference; bins
# past N_sel_bins repeat the last coefficient and contribute d = 0.
testthat::test_that("the form-9 shape penalty charges the direction Sel_shape_dir names", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  d   <- Rceattle::BS2017SS
  np  <- which(d$fleet_control$Selectivity == 2)
  flt <- np[1]
  d$fleet_control$Selectivity[np]       <- "NonParametricPM"
  d$fleet_control$Time_varying_sel[np]  <- "Off"

  bld <- function(pen1, inits) {
    dd <- d; dd$fleet_control$Sel_curve_pen1[np] <- pen1
    suppressMessages(suppressWarnings(Rceattle::fit_mod(
      data_list = dd, msmMode = 0, estimateMode = 3, inits = inits,
      fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
  }
  m0  <- bld(20, NULL)
  nsb <- m0$obj$env$data$flt_n_sel_bins[flt]
  nb  <- m0$obj$env$data$nages[d$fleet_control$Species[flt]]
  # Address the jnll row through the registry, never a bare integer: .JNLL_ROW_AXIS
  # is keyed by display name in JnllRow order (R/6-rename_output.R).
  .np_row <- which(names(Rceattle:::.JNLL_ROW_AXIS) == "Non-parametric selectivity")
  testthat::expect_length(.np_row, 1L)
  # The oracle sums one sex and writes sex 1 only, so it is exact only on a
  # one-sex fixture; the C++ loops every sex.
  testthat::expect_true(all(d$nsex == 1))
  testthat::expect_equal(dim(m0$initial_params$sel_coff)[3], nsb)
  # It also differences every adjacent pair from bin 1 and ignores the cap, so it
  # matches the C++ only on a fleet with no excluded leading bins, no cap and no
  # narrowed penalty range -- and it reads JNLL_SEL_NONPARAM, which the avgsel
  # penalty would also enter. Asserted rather than assumed: lifting this oracle
  # to a fleet like EBS pollock (excluded age-1 bin, AMAK caps) gives a wrong
  # answer, and these are what would make it wrong.
  testthat::expect_equal(as.integer(d$fleet_control$Bin_first_selected[flt]), 1L)
  testthat::expect_true(is.null(d$fleet_control$Sel_cap_bin) ||
                          is.na(d$fleet_control$Sel_cap_bin[flt]))
  testthat::expect_true(is.null(d$fleet_control$Sel_pen_first_bin) ||
                          is.na(d$fleet_control$Sel_pen_first_bin[flt]))
  testthat::expect_true(is.null(d$fleet_control$Sel_pen_last_bin) ||
                          is.na(d$fleet_control$Sel_pen_last_bin[flt]))
  testthat::expect_true(is.null(d$fleet_control$Sel_avgsel_pen) ||
                          is.na(d$fleet_control$Sel_avgsel_pen[flt]) ||
                          d$fleet_control$Sel_avgsel_pen[flt] == 0)

  oracle <- function(coff, w, increasing) {
    base <- coff[pmin(seq_len(nb), nsb)]        # tail-repeat past N_sel_bins
    dif  <- base[1:(nb - 1)] - base[2:nb]       # > 0 where selectivity decreases
    w * sum(pmax(if (increasing) -dif else dif, 0)^2)
  }
  # sel_curve_pen is a PARAMETER, so `inits` carry the weight and would override
  # the fleet_control column; set both or the sign never reaches the template.
  charged <- function(ramp, pen1) {
    p <- m0$initial_params
    p$sel_coff[flt, 1, ]     <- ramp
    p$sel_curve_pen[flt, 1]  <- pen1
    r <- bld(pen1, p)$obj$report()$jnll_comp
    r[.np_row, flt]
  }

  up   <- seq(0, 1.4, length.out = nsb)
  down <- rev(up)
  # Asymmetric on purpose: a symmetric tent charges both directions equally and
  # would not show that the sign selects which segments are scored.
  .k   <- ceiling(nsb * 2 / 3)
  mix  <- c(seq(0, 1.4, length.out = .k),
            seq(1.4, 0.7, length.out = nsb - .k + 1L)[-1])

  # The limiting cases: each direction is blind to a curve that never moves the
  # way it penalizes, and charges the mirror case.
  testthat::expect_equal(charged(up,   20), 0)
  testthat::expect_equal(charged(down, -20), 0)
  testthat::expect_gt(charged(up,   -20), 0)
  testthat::expect_gt(charged(down,  20), 0)

  # And the amount charged is the penalty it claims to be, on all three shapes.
  for (r in list(up, down, mix)) {
    testthat::expect_equal(charged(r,  20), oracle(r, 20, FALSE), tolerance = 1e-10)
    testthat::expect_equal(charged(r, -20), oracle(r, 20, TRUE),  tolerance = 1e-10)
  }
  # A mixed curve charges both directions, so the sign genuinely selects which
  # segments are scored rather than rescaling one set.
  testthat::expect_gt(charged(mix, 20), 0)
  testthat::expect_gt(charged(mix, -20), 0)
  testthat::expect_false(isTRUE(all.equal(charged(mix, 20), charged(mix, -20))))

  # The same through the column a user actually writes: Sel_shape_dir is what
  # produces the negative weight, so drive the whole chain rather than the
  # weight alone.
  via_dir <- function(dir) {
    dd <- d
    dd$fleet_control$Sel_curve_pen1[np] <- NA_real_
    dd$fleet_control$Sel_shape_sd  <- NA_real_
    dd$fleet_control$Sel_shape_sd[np]  <- 1 / sqrt(2 * 20)
    dd$fleet_control$Sel_shape_dir <- "Decreasing"
    dd$fleet_control$Sel_shape_dir[np] <- dir
    suppressMessages(Rceattle::switch_check(dd))$fleet_control$Sel_curve_pen1[flt]
  }
  testthat::expect_equal(via_dir("Decreasing"),  20)
  testthat::expect_equal(via_dir("Increasing"), -20)
  p <- m0$initial_params
  p$sel_coff[flt, 1, ]    <- up
  p$sel_curve_pen[flt, 1] <- via_dir("Increasing")
  testthat::expect_equal(bld(via_dir("Increasing"), p)$obj$report()$jnll_comp[.np_row, flt],
                         oracle(up, 20, TRUE), tolerance = 1e-10)
})

testthat::test_that("Sel_shape_mode = 'Smooth' is not sign-aware, so it takes no direction", {
  # "Smooth" applies the weight two-sided (pen * d^2), so a negative weight
  # there rewards the shape exactly as on the forms that never read a sign.
  d  <- Rceattle::BS2017SS
  np <- which(d$fleet_control$Selectivity == 2)
  d$fleet_control$Selectivity[np]     <- "NonParametricPM"
  d$fleet_control$Sel_curve_pen1[np]  <- NA_real_
  d$fleet_control$Sel_shape_mode      <- "Smooth"

  dir_run <- d
  dir_run$fleet_control$Sel_shape_sd      <- NA_real_
  dir_run$fleet_control$Sel_shape_sd[np]  <- 1 / sqrt(2 * 20)
  dir_run$fleet_control$Sel_shape_dir     <- "Decreasing"
  dir_run$fleet_control$Sel_shape_dir[np] <- "Increasing"
  testthat::expect_error(suppressMessages(Rceattle::switch_check(dir_run)),
                         "does not read the sign")

  col_run <- d
  col_run$fleet_control$Sel_curve_pen1[np] <- -20
  col_run$fleet_control$Sel_curve_pen2[np] <- 12.5
  testthat::expect_error(
    suppressMessages(Rceattle:::data_check(suppressMessages(Rceattle::switch_check(col_run)))),
    "rewards the deviation")

  # "Directional" is unaffected: it is the branch that reads the sign.
  ok <- dir_run; ok$fleet_control$Sel_shape_mode <- "Directional"
  testthat::expect_equal(
    suppressMessages(Rceattle::switch_check(ok))$fleet_control$Sel_curve_pen1[np[1]], -20)
})

testthat::test_that("the direction exemption survives a numeric Selectivity code", {
  # switch_check() canonicalizes Selectivity late, so this rule sees whatever the
  # workbook holds -- and every bundled data set stores the integer code.
  d  <- Rceattle::BS2017SS
  np <- which(d$fleet_control$Selectivity == 2)
  testthat::expect_true(is.numeric(d$fleet_control$Selectivity))
  d$fleet_control$Selectivity[np]    <- 9          # "NonParametricPM"
  d$fleet_control$Sel_curve_pen1[np] <- NA_real_
  d$fleet_control$Sel_shape_sd       <- NA_real_
  d$fleet_control$Sel_shape_sd[np]   <- 1 / sqrt(2 * 20)
  d$fleet_control$Sel_shape_dir      <- "Decreasing"
  d$fleet_control$Sel_shape_dir[np]  <- "Increasing"
  out <- suppressMessages(Rceattle::switch_check(d))
  testthat::expect_equal(out$fleet_control$Sel_curve_pen1[np], rep(-20, length(np)),
                         tolerance = 1e-8)

  # And the column check agrees when handed the raw codes directly.
  raw <- d; raw$fleet_control$Sel_curve_pen1[np] <- -20
  testthat::expect_length(Rceattle:::.rce_sel_pen_sign_errors(raw$fleet_control), 0L)
  raw$fleet_control$Selectivity[np] <- 2           # "NonParametric": not sign-aware
  testthat::expect_true(length(Rceattle:::.rce_sel_pen_sign_errors(raw$fleet_control)) > 0L)
})

testthat::test_that("a negative Sel_curve_pen in `inits` is refused, not silently used", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  # sel_curve_pen is a PARAMETER, so `inits` from a fit saved before 5.42.0
  # override the corrected column. Every refit path (retrospective(), profile(),
  # run_mse(), self_test()) passes estimated_params as inits, so without this
  # the refusal would not reach any of them: measured -3.81 on a decreasing
  # ramp, -28.99 at three times the ramp.
  d   <- Rceattle::BS2017SS
  flt <- which(d$fleet_control$Selectivity == 2)[1]
  bld <- function(inits) suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = d, msmMode = 0, estimateMode = 3, inits = inits,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
  m0 <- bld(NULL)
  testthat::expect_no_error(bld(m0$initial_params))     # a clean refit is untouched
  bad <- m0$initial_params
  bad$sel_curve_pen[flt, 1] <- -20
  # Names inits as the source, so this pins the parameter check rather than the
  # column check, which still sees a valid +20.
  testthat::expect_error(bld(bad), "in the supplied `inits`", fixed = TRUE)
})

# The check skips a fleet that follows another's Selectivity_index, because the
# shared block's penalty is charged once, on its lead. Which fleet leads is the
# question: ceattle.cpp reads flt_sel_lead, which rearrange_data() groups by
# Selectivity_index AND selectivity form, while the parameter map shares on the
# index alone. Borrowing the map's rule let a negative weight through on a fleet
# the template charges. Found reviewing the 5.34.0-5.42.0 release PR (#158).
testthat::test_that("the penalty lead follows the template's grouping, not the map's", {
  d  <- Rceattle::BS2017SS
  fc <- d$fleet_control
  # Fleet 5 is made to follow fleet 4's Selectivity_index and carry the negative
  # weight; only the FORM differs between the two cases.
  mk <- function(same_form) {
    x <- fc
    x$Selectivity_index[5] <- x$Selectivity_index[4]
    if (!same_form) x$Selectivity[5] <- 2
    x$Sel_curve_pen1[5] <- -20
    x
  }

  # Same form: one group, charged once on fleet 4. Fleet 5's weight is never
  # read, so a negative one there is inert and stays allowed.
  same <- mk(TRUE)
  testthat::expect_false(Rceattle:::.rce_sel_pen_lead(same)[5])
  testthat::expect_length(Rceattle:::.rce_sel_pen_sign_errors(same), 0L)

  # Different form: two groups to the template, and fleet 5 leads its own. The
  # weight IS read, so it must be refused. This is the case that slipped through.
  mixed <- mk(FALSE)
  testthat::expect_true(Rceattle:::.rce_sel_pen_lead(mixed)[5])
  testthat::expect_match(Rceattle:::.rce_sel_pen_sign_errors(mixed),
                         "Sel_curve_pen1 is negative", all = FALSE)

  # The rule the check used to borrow calls fleet 5 a follower in BOTH cases --
  # this is the divergence itself, asserted so a revert is loud.
  for (x in list(same, mixed)) {
    testthat::expect_false(
      is.na(Rceattle:::.shared_block_lead(list(fleet_control = x), 5L, "sel")))
  }

  # BS2017SS gives every fleet its own Selectivity_index, so every fleet leads.
  testthat::expect_true(all(Rceattle:::.rce_sel_pen_lead(fc)))

  # The lead of a group is its first fleet that is NOT Off: an Off fleet
  # estimates nothing, and the template gates on flt_type as well, so leading
  # with one would drop the group's penalty entirely. Asserted here because
  # .rce_sel_pen_sign_errors() has its own Off skip, so without this the `off`
  # argument to .group_lead() could be dropped and this file would stay green.
  off_fc <- fc
  off_fc$Selectivity_index[5] <- off_fc$Selectivity_index[4]
  off_fc$Fleet_type[4] <- "Off"
  testthat::expect_false(Rceattle:::.rce_sel_pen_lead(off_fc)[4])  # Off, never leads
  testthat::expect_true(Rceattle:::.rce_sel_pen_lead(off_fc)[5])   # so fleet 5 does

  # And it matches what rearrange_data() hands the template for the same table.
  dl <- Rceattle::BS2017SS; dl$fleet_control <- off_fc
  lead_tmb <- suppressMessages(suppressWarnings(
    Rceattle:::rearrange_data(Rceattle:::switch_check(dl))))$flt_sel_lead
  testthat::expect_identical(as.integer(Rceattle:::.rce_sel_pen_lead(off_fc)), lead_tmb)
})

# .canon_switch() trims and rearrange_data()'s .pull_int() does not, so a
# Selectivity the R side resolves and the template does not would group the
# fleet here and leave it leading there -- the guard would skip a weight
# ceattle.cpp charges. Such a value gets a key of its own, so the fleet leads
# here too. switch_check() normalizes the spelling first, so this guards a
# hand-built fleet_control rather than a workbook path.
testthat::test_that("a Selectivity the template cannot read leads, so its weight is checked", {
  fc <- Rceattle::BS2017SS$fleet_control
  fc$Selectivity <- as.character(fc$Selectivity)
  fc$Selectivity[4] <- " NonParametric"          # resolves here, NA to the template
  fc$Selectivity[5] <- "NonParametric"
  fc$Selectivity_index[5] <- fc$Selectivity_index[4]
  fc$Sel_curve_pen1[5] <- -20

  testthat::expect_true(Rceattle:::.rce_sel_pen_lead(fc)[5])
  testthat::expect_match(Rceattle:::.rce_sel_pen_sign_errors(fc),
                         "Sel_curve_pen1 is negative", all = FALSE)

  # Trailing whitespace is the same case, and a clean spelling still groups.
  fc$Selectivity[4] <- "NonParametric "
  testthat::expect_true(Rceattle:::.rce_sel_pen_lead(fc)[5])
  fc$Selectivity[4] <- "NonParametric"
  testthat::expect_false(Rceattle:::.rce_sel_pen_lead(fc)[5])
  testthat::expect_length(Rceattle:::.rce_sel_pen_sign_errors(fc), 0L)
})
