# fit_control(selectivity_se = TRUE) adds a delta-method standard error for the
# selectivity curve, asked for in issue #107 so a fitted curve can be drawn with
# an interval.
#
# Reported on the log scale rather than the logit that issue proposes. The
# non-parametric forms renormalize to mean selectivity 1, so sel_at_age is not
# bounded by 1: on Atka2022, 58% of entries exceed it and the largest is 3.06.
# logit() is undefined for those. Atka2022 is the fixture precisely because it
# is one of those models.
#
# The recurring hazard in this file is a zero. log(0) = -Inf on the AD tape
# makes sdreport return NaN for EVERY quantity, biomass and ssb included, not
# just selectivity, so each test below pins one route to a zero cell shut:
# Fixed selectivity, a length-based fleet, and Bin_first_selected.

testthat::test_that("selectivity_se adds a log-scale SE without moving the fit", {
  testthat::skip_on_cran()

  run <- function(se) suppressMessages(suppressWarnings(fit_mod(
    data_list    = Atka2022,
    msmMode      = 0,
    estimateMode = "Hindcast",
    fit_control  = fit_control(selectivity_se = se, verbose = 0))))

  off <- run(FALSE)
  on  <- run(TRUE)

  # An ADREPORT is a report, not a likelihood term.
  testthat::expect_identical(on$opt$objective, off$opt$objective)

  testthat::expect_false("log_sel_at_age" %in% names(off$sdrep$value))
  testthat::expect_true("log_sel_at_age" %in% names(on$sdrep$value))

  # One row per reported fleet x its own species' sexes x selected ages x
  # hindcast year. Iterating the padded array instead would report structural
  # zeros, and log(0) = -Inf NaNs every other quantity in the same sdreport.
  # Built from the same rearrange_data() switches the template read -- they live
  # on obj$env$data, not on data_list -- so this restates the skip rule rather
  # than re-deriving it from the fleet_control columns.
  # bin_first_selected is the 0-based count of ages the curve zeroes below.
  dat <- on$obj$env$data
  dl  <- on$data_list
  spp <- dl$fleet_control$Species
  lead <- which(dat$flt_type > 0 & dat$flt_sel_lead == 1 &
                  dat$flt_sel_type != 0 & !(dat$flt_sel_dim %in% 1L))
  testthat::expect_gt(length(lead), 0)
  nyrs_hind <- dl$endyr - dl$styr + 1
  nbin <- dl$nages[spp[lead]] - dat$bin_first_selected[lead]
  expected <- sum(dl$nsex[spp[lead]] * nbin) * nyrs_hind

  k <- names(on$sdrep$value) == "log_sel_at_age"
  testthat::expect_equal(sum(k), expected)

  # Columns 3 and 4 are the absolute age and the calendar year, not 0-based or
  # 1-based positions. nages counts bins, so age and age index differ wherever
  # minage > 1 -- which no bundled dataset has, so only this assertion catches
  # a regression to raw indices.
  idx <- on$quantities$log_sel_at_age_index
  testthat::expect_equal(nrow(idx), expected)
  testthat::expect_equal(ncol(idx), 4L)
  testthat::expect_setequal(unique(idx[, 4]), on$data_list$styr:on$data_list$endyr)
  minage_sp <- on$data_list$minage[spp[lead]]
  testthat::expect_gte(min(idx[, 3]), min(minage_sp))
  testthat::expect_true(all(idx[, 1] %in% lead))
  testthat::expect_true(all(idx[, 2] >= 1))

  # The index has to address the right cell, or an interval gets drawn against
  # the wrong age. Convert back to array positions to look it up.
  sel  <- on$quantities$sel_at_age
  back <- exp(on$sdrep$value[k])
  age_i <- idx[, 3] - on$data_list$minage[spp[idx[, 1]]] + 1L
  yr_i  <- idx[, 4] - on$data_list$styr + 1L
  at   <- sel[cbind(idx[, 1], idx[, 2], age_i, yr_i)]
  testthat::expect_equal(unname(back), unname(at), tolerance = 1e-12)
  testthat::expect_true(all(at > 0))   # a zero would have put -Inf on the tape

  # Finite everywhere, and the rest of the sdreport is untouched by it.
  sd_sel <- on$sdrep$sd[k]
  testthat::expect_true(all(is.finite(sd_sel)))
  testthat::expect_false(any(is.nan(on$sdrep$sd)))

  # A bin the model pins carries no uncertainty. Atka2022's survey normalizes at
  # Sel_norm_bin = 4, so that age is 1 by construction with an SE of 0.
  testthat::expect_true(any(sd_sel == 0))
})


# The error belongs to whichever sdreport the fit ends on. estimateMode
# "Projection" estimates no selectivity at all, so it has none to report, and
# reading a 0 as "the curve is certain" is the failure this warning stops.
testthat::test_that("estimateMode 'Projection' warns that it reports no usable error", {
  testthat::skip_on_cran()

  # projection_uncertainty turns the hindcast parameters back on under
  # "Estimate"; here build_map(debug = TRUE) has already mapped them off, so it
  # cannot help and the warning must fire anyway.
  testthat::expect_warning(
    m <- suppressMessages(fit_mod(
      data_list = Atka2022, msmMode = 0, estimateMode = "Projection",
      fit_control = fit_control(selectivity_se = TRUE,
                                projection_uncertainty = TRUE, verbose = 0))),
    "no usable standard error")

  # The plotter reaches the same verdict from the fit alone, so add_ci declines
  # rather than drawing a hairline band.
  testthat::expect_null(.rce_sel_se(m))
})


testthat::test_that("selectivity_se defaults to off and survives a data list without it", {
  # The TMB template declares adreport_sel as DATA_INTEGER, so MakeADFun() fails
  # outright if the field is missing. A data list built before the option existed
  # has no such field, so build_osa_data() -- the funnel every list passes
  # through on its way to MakeADFun() -- has to supply it. Asserted on that
  # funnel directly: fit_mod() sets the field itself, so a fit would pass whether
  # the default existed or not.
  testthat::expect_false(fit_control()$selectivity_se)

  d <- suppressMessages(suppressWarnings(
    rearrange_data(switch_check(Atka2022))))
  d$adreport_sel <- NULL
  testthat::expect_identical(
    build_osa_data(d, build_osa = FALSE)$adreport_sel, 0L)

  # And an explicit TRUE survives the funnel rather than being reset to the
  # default, which is what carries the setting through a refit.
  d$adreport_sel <- 1L
  testthat::expect_identical(
    build_osa_data(d, build_osa = FALSE)$adreport_sel, 1L)

  # adreport_catch is the same kind of DATA_INTEGER, off unless a peel's
  # correction sets it.
  d$adreport_catch <- NULL
  testthat::expect_identical(
    build_osa_data(d, build_osa = FALSE)$adreport_catch, 0L)
})


# BS2017SS is the regression case for the Fixed-selectivity skip: its
# EIT_Pollock survey is Selectivity = "Fixed", which leaves 0 in every year with
# no emp_sel_obs row -- 228 cells.
testthat::test_that("a Fixed-selectivity fleet cannot NaN the sdreport", {
  testthat::skip_on_cran()

  m <- suppressMessages(suppressWarnings(fit_mod(
    data_list = BS2017SS, msmMode = 0, estimateMode = "Hindcast",
    fit_control = fit_control(selectivity_se = TRUE, verbose = 0))))

  testthat::expect_false(any(is.nan(m$sdrep$sd)))
  testthat::expect_false(any(is.nan(m$sdrep$value)))

  k <- names(m$sdrep$value) == "log_sel_at_age"
  testthat::expect_gt(sum(k), 0)
  testthat::expect_true(all(is.finite(m$sdrep$value[k])))

  # The Fixed fleet contributes no rows.
  fc <- m$data_list$fleet_control
  fixed <- fc$Fleet_code[fc$Selectivity == "Fixed"]
  idx <- m$quantities$log_sel_at_age_index
  testthat::expect_length(intersect(idx[, 1], fixed), 0)
})


# Bin_first_selected is the second route to a zero, and the one a real
# assessment reaches first: selectivity.hpp zeroes the curve outright below it,
# for every estimated fleet, before normalizing. No bundled dataset sets it above
# minage, but 117 workbooks across the sibling assessment repositories do --
# GOA pollock 2024 (Age_first_selected = 3 on two fleets) and the Pacific hake
# MSE model among them -- so the fixture sets it here rather than waiting for a
# real fit to find it. Measured on Atka2022 with Bin_first_selected = 3: 92 of
# fleet 1's 506 hindcast sel_at_age cells are exactly 0.
testthat::test_that("Bin_first_selected cannot NaN the sdreport", {
  testthat::skip_on_cran()

  d <- Atka2022
  d$fleet_control$Bin_first_selected <- c(3, 1)

  m <- suppressMessages(suppressWarnings(fit_mod(
    data_list = d, msmMode = 0, estimateMode = "Hindcast",
    fit_control = fit_control(selectivity_se = TRUE, verbose = 0))))

  # The premise: the curve really is zero there, so this is not a vacuous pass.
  sel <- m$quantities$sel_at_age
  nyrs_hind <- m$data_list$endyr - m$data_list$styr + 1
  testthat::expect_gt(sum(sel[1, 1, seq_len(m$data_list$nages[1]),
                              seq_len(nyrs_hind)] == 0), 0)

  testthat::expect_false(any(is.nan(m$sdrep$sd)))
  testthat::expect_false(any(is.nan(m$sdrep$value)))

  k <- names(m$sdrep$value) == "log_sel_at_age"
  testthat::expect_true(all(is.finite(m$sdrep$value[k])))

  # Ages below Bin_first_selected are absent, not reported as -Inf. Ages are
  # absolute, so the first one reported is the column value itself.
  idx <- m$quantities$log_sel_at_age_index
  testthat::expect_equal(min(idx[idx[, 1] == 1, 3]), 3)
  testthat::expect_equal(min(idx[idx[, 1] == 2, 3]), 1)
})


# The fourth route is not structural: a reported cell can underflow to 0. A
# double-normal peaking above the oldest age with a narrow ascending limb puts
# exp(-x^2/2) below 1e-308, and log(0) = -Inf NaNs every quantity in the
# sdreport. The log is floored at 1e-300, which is exact above the floor.
testthat::test_that("a selectivity that underflows to zero is floored, not -Inf", {
  testthat::skip_on_cran()

  d <- Atka2022
  d$fleet_control$Selectivity <- c("DoubleNormal", d$fleet_control$Selectivity[2])
  build <- function(inits) suppressMessages(suppressWarnings(fit_mod(
    data_list = d, msmMode = 0, estimateMode = "DebugBuild", inits = inits,
    fit_control = fit_control(selectivity_se = TRUE, getsd = FALSE,
                              verbose = 0))))

  p <- build(NULL)$obj$env$parList()
  p$sel_inf[1, 1, ]     <- 60          # peak far above the oldest age
  p$log_sel_slp[1, 1, ] <- log(0.15)   # narrow ascending limb
  m <- build(p)

  # The premise: the curve really is 0 at reported ages, so this is not vacuous.
  s <- m$quantities$sel_at_age[1, 1, seq_len(d$nages[1]), 1]
  testthat::expect_gt(sum(s == 0), 0)

  ls <- m$quantities$log_sel_at_age
  testthat::expect_true(all(is.finite(ls)))
  testthat::expect_equal(min(ls), log(1e-300))
})


# The third route: a length-based fleet's sel_at_age is
# sum(growth_matrix * sel_at_length), which is 0 for an age overlapping no
# selected length bin. Such a fleet is fitted and plotted on sel_at_length, which
# carries no error, so it is skipped rather than guarded.
testthat::test_that("a length-based fleet reports no age-based standard error", {
  testthat::skip_on_cran()

  # Whole columns, not one element: Atka2022's workbook carries neither, and
  # assigning into an absent column recycles the value over every fleet.
  d <- Atka2022
  d$fleet_control$Selectivity_dimension <- c("Length", "Age")
  d$fleet_control$Selectivity <- c("Logistic", d$fleet_control$Selectivity[2])

  m <- suppressMessages(suppressWarnings(fit_mod(
    data_list = d, msmMode = 0, estimateMode = "Hindcast",
    fit_control = fit_control(selectivity_se = TRUE, verbose = 0))))
  testthat::expect_equal(m$obj$env$data$flt_sel_dim, c(1L, 0L))

  testthat::expect_false(any(is.nan(m$sdrep$sd)))
  idx <- m$quantities$log_sel_at_age_index
  testthat::expect_length(intersect(idx[, 1], 1L), 0)
  testthat::expect_gt(sum(idx[, 1] == 2L), 0)
})
