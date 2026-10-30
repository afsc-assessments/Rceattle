# A Fixed fleet's selectivity skipped the step that carries estimated selectivity
# through the projection, so its sel_at_age was 0 in every projection year: a Fixed
# survey predicted an index and comps of 0 there, and a Fixed fishery had no
# exploitable biomass and, under an HCR that reads SBF, the wrong F. Found reviewing
# the age-selectivity comp fix, 2026-10-30.

# One row at Year 0 applies a curve to every hindcast year; `years` supplies only those.
.fixed_sel <- function(d, flt, sel, years = 0) {
  d$fleet_control$Selectivity[flt] <- "Fixed"
  d$emp_sel <- cbind(
    data.frame(Fleet_name = as.character(d$fleet_control$Fleet_name[flt]), Fleet_code = flt,
               Species = 1L, Sex = 0L, Year = years),
    as.data.frame(matrix(sel, length(years), length(sel), byrow = TRUE,
                         dimnames = list(NULL, paste0("Comp_", seq_along(sel))))))
  d
}

testthat::test_that("a Fixed survey predicts its index and comps in projection years", {
  testthat::skip_if_not_installed("TMB")
  d <- make_test_data(nyrs = 8, nages = 5, seed = 1)
  sel <- c(0.1, 0.4, 0.8, 1, 1)
  d <- .fixed_sel(d, 1, sel)
  py <- d$endyr + 2

  comp <- d$comp_data[d$comp_data$Fleet_code == 1, ][1, ]
  comp$Year <- -py
  d$comp_data <- rbind(d$comp_data, comp)
  ind <- d$index_data[d$index_data$Fleet_code == 1, ][1, ]
  ind$Year <- -py
  d$index_data <- rbind(d$index_data, ind)

  q <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    d, file = NULL, estimateMode = "DebugBuild", msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0))))$quantities

  yr <- py - d$styr + 1
  testthat::expect_equal(as.numeric(q$sel_at_age[1, 1, , yr]), sel)
  testthat::expect_gt(tail(as.numeric(q$index_hat), 1), 0)

  # Survey numbers at age in that year, by hand; age comps carry the identity ageing error.
  x <- q$N_at_age[1, 1, , yr] * sel *
    exp(-d$fleet_control$Month[1] / 12 * q$Z_at_age[1, 1, , yr])
  testthat::expect_equal(as.numeric(q$comp_hat[nrow(d$comp_data), 1:5]),
                         as.numeric(x / sum(x)), tolerance = 1e-10)
})

testthat::test_that("a Fixed fishery projects as an estimated fishery with the same curve", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  d <- make_test_data(nyrs = 12, nages = 5, seed = 1)
  hcr <- Rceattle::build_hcr(HCR = "NPFMC", DynamicHCR = FALSE, Ftarget = 0.4, Flimit = 0.35)
  ctl <- Rceattle::fit_control(getsd = FALSE, verbose = 0, phase = FALSE,
                               nlminb_control = list(iter.max = 10))
  fit <- function(d, inits, mode) suppressWarnings(suppressMessages(Rceattle::fit_mod(
    d, inits = inits, file = NULL, estimateMode = mode, msmMode = 0, random_rec = FALSE,
    HCR = hcr, fit_control = ctl)))

  est <- fit(d, NULL, "Estimate")
  ny  <- length(d$styr:d$endyr)
  # The estimated curve is time-invariant here, so one row is exact.
  curve <- est$quantities$sel_at_age[2, 1, , seq_len(ny)]
  testthat::expect_lt(max(apply(curve, 1, function(s) diff(range(s)))), 1e-12)
  fixed <- fit(.fixed_sel(d, 2, curve[, 1]), est$estimated_params, "Projection")

  py <- seq(ny + 1, dim(est$quantities$sel_at_age)[4])
  testthat::expect_gt(max(est$quantities$proj_F), 0)
  testthat::expect_equal(as.numeric(fixed$quantities$sel_at_age[2, 1, , py]),
                         as.numeric(est$quantities$sel_at_age[2, 1, , py]), tolerance = 1e-12)
  testthat::expect_equal(as.numeric(fixed$quantities$Ftarget), as.numeric(est$quantities$Ftarget),
                         tolerance = 1e-5)
  testthat::expect_equal(as.numeric(fixed$quantities$proj_F), as.numeric(est$quantities$proj_F),
                         tolerance = 1e-5)
  testthat::expect_equal(as.numeric(fixed$quantities$ssb), as.numeric(est$quantities$ssb),
                         tolerance = 1e-5)
})

testthat::test_that("a Fixed fleet holds its last supplied year, not its last hindcast year", {
  testthat::skip_if_not_installed("TMB")
  # BS2017SS's EIT survey has this shape: emp_sel rows stop a year before endyr.
  d <- make_test_data(nyrs = 8, nages = 5, seed = 1)
  sel <- c(0.1, 0.4, 0.8, 1, 1)
  d <- .fixed_sel(d, 1, sel, years = d$styr:(d$endyr - 1))

  q <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    d, file = NULL, estimateMode = "DebugBuild", msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0))))$quantities
  ny <- length(d$styr:d$endyr)
  testthat::expect_equal(as.numeric(q$sel_at_age[1, 1, , ny - 1]), sel)
  testthat::expect_equal(sum(q$sel_at_age[1, 1, , ny]), 0)
  testthat::expect_equal(as.numeric(q$sel_at_age[1, 1, , ny + 1]), sel)
})

testthat::test_that("data_check() warns when a Fixed fishery has no row for endyr", {
  testthat::skip_if_not_installed("TMB")
  warns <- function(d) {
    d$HCR <- 0
    x <- suppressMessages(suppressWarnings(Rceattle::switch_check(Rceattle::clean_data(d))))
    w <- character(0)
    withCallingHandlers(suppressMessages(data_check(x)), warning = function(e) {
      w <<- c(w, conditionMessage(e))
      invokeRestart("muffleWarning")
    })
    any(grepl("no emp_sel row for endyr", w))
  }
  d <- make_test_data(nyrs = 8, nages = 5, seed = 1)
  sel <- c(0.1, 0.4, 0.8, 1, 1)

  testthat::expect_true(warns(.fixed_sel(d, 2, sel, years = d$styr:(d$endyr - 1))))
  testthat::expect_false(warns(.fixed_sel(d, 2, sel, years = d$styr:d$endyr)))
  testthat::expect_false(warns(.fixed_sel(d, 2, sel)))
  # A survey with this shape predicts nothing that year; BS2017SS's EIT has no data there.
  testthat::expect_false(warns(.fixed_sel(d, 1, sel, years = d$styr:(d$endyr - 1))))
})

testthat::test_that("the endyr warning keys on sex, and a blank Selectivity still gets its own error", {
  testthat::skip_if_not_installed("TMB")
  d <- Rceattle::GOAatf
  flt <- which(as.character(d$fleet_control$Fleet_type) %in% c("1", "Fishery"))[1]
  d$fleet_control$Selectivity[flt] <- "Fixed"
  d$comp_data <- d$comp_data[d$comp_data$Fleet_code != flt, ]
  row <- function(sex, year) {
    cbind(data.frame(Fleet_name = as.character(d$fleet_control$Fleet_name[flt]), Fleet_code = flt,
                     Species = 1L, Sex = sex, Year = year),
          as.data.frame(matrix(1, 1, d$nages[1], dimnames = list(NULL, paste0("Comp_", seq_len(d$nages[1]))))))
  }
  w <- function(es) {
    x <- d
    x$emp_sel <- es
    got <- character(0)
    withCallingHandlers(.warn_fixed_fishery_endyr(x), warning = function(e) {
      got <<- c(got, conditionMessage(e))
      invokeRestart("muffleWarning")
    })
    length(got) > 0
  }
  testthat::expect_true(w(row(1L, d$endyr)))
  testthat::expect_false(w(rbind(row(1L, d$endyr), row(2L, d$endyr))))
  testthat::expect_false(w(row(0L, d$endyr)))

  # An unset Selectivity must reach switch_check()'s message, not an `if (NA)`.
  blank <- make_test_data(nyrs = 8, nages = 5, seed = 1, growth = "vonBertalanffy")
  blank$fleet_control$Selectivity[2] <- NA
  blank$fleet_control$Selectivity_dimension <- "Length"
  testthat::expect_error(
    suppressWarnings(suppressMessages(Rceattle::fit_mod(
      blank, file = NULL, estimateMode = "DebugBuild", msmMode = 0,
      growthFun = Rceattle::build_growth(fun = "vonBertalanffy"),
      fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0)))),
    "Invalid 'Selectivity'")
})
