# NonParametricIntegrable (13), 5.40.0: the Ianelli base curve whose deviations
# have a proper density, with Time_varying_sel picking the structure.
# NonParametric (2) charges its shape penalties on each year's realized curve,
# so integrating its deviates tilts the density and biases the deviation SD.
# This form charges the decreasing, curvature and average-selectivity penalties
# once on the base coefficients and scores the deviates with
# dnorm(0, sel_dev_sd) alone.
# Atka2022's fishery ships as NonParametric with IID deviates (sd 0.35), the
# configuration this exists for. estimateMode = 3 evaluates at the starts.

np_data <- function(form, tv) {
  d <- Rceattle::Atka2022
  d$fleet_control$Selectivity[2]      <- form
  d$fleet_control$Time_varying_sel[2] <- tv
  d
}

np_build <- function(d, inits = NULL, random_sel = FALSE, mode = 3) {
  suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = d, inits = inits, msmMode = 0, estimateMode = mode,
    random_sel = random_sel,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
}

np_rows <- function(m) m$quantities$jnll_comp[c("Non-parametric selectivity", "Selectivity deviates"), 2]

testthat::test_that("with no deviates the form gives the NonParametric objective", {
  testthat::skip_on_cran()
  # A curved base, so the decreasing and curvature penalties are not zero.
  bent <- function(m) { ini <- m$initial_params; ini$sel_coff[2, 1, 1:10] <- log(seq(0.2, 1, length.out = 10)); ini }
  m2  <- np_build(np_data("NonParametric",    "Off")); m2  <- np_build(np_data("NonParametric",    "Off"), inits = bent(m2))
  m13 <- np_build(np_data("NonParametricIntegrable", "Off")); m13 <- np_build(np_data("NonParametricIntegrable", "Off"), inits = bent(m13))
  testthat::expect_gt(abs(np_rows(m2)[["Non-parametric selectivity"]]), 0)
  testthat::expect_equal(m13$obj$fn(), m2$obj$fn(), tolerance = 1e-12)
  testthat::expect_equal(np_rows(m13), np_rows(m2), tolerance = 1e-12)
  testthat::expect_equal(m13$quantities$sel_at_age[2, 1, , ], m2$quantities$sel_at_age[2, 1, , ])
})

testthat::test_that("the Off equivalence holds with a first selected bin above 1", {
  testthat::skip_on_cran()
  bfs <- function(form) { d <- np_data(form, "Off"); d$fleet_control$Bin_first_selected[2] <- 3; d }
  bent <- function(m) { ini <- m$initial_params; ini$sel_coff[2, 1, 3:10] <- log(seq(0.2, 1, length.out = 8)); ini }
  m2  <- np_build(bfs("NonParametric"));    m2  <- np_build(bfs("NonParametric"),    inits = bent(m2))
  m13 <- np_build(bfs("NonParametricIntegrable")); m13 <- np_build(bfs("NonParametricIntegrable"), inits = bent(m13))
  testthat::expect_equal(m13$obj$fn(), m2$obj$fn(), tolerance = 1e-12)

  # The walk construction itself: built under RandomWalk so the walk branch of
  # the template runs, with the increments at their zero starting values. The
  # curve must then be flat across years. Checked on the curve rather than the
  # objective, which differs by the increments' own density.
  dw <- bfs("NonParametricIntegrable"); dw$fleet_control$Time_varying_sel[2] <- "RandomWalk"
  mw <- np_build(dw); mw <- np_build(dw, inits = bent(mw))
  sw <- mw$quantities$sel_at_age[2, 1, , ]
  testthat::expect_equal(sw[, ncol(sw)], sw[, 1])
  testthat::expect_equal(sw[, 2], sw[, 1])
})

testthat::test_that("NonParametricIntegrable reads no increment at or before the fleet's start year", {
  testthat::skip_on_cran()
  d <- np_data("NonParametricIntegrable", "RandomWalk")
  d$fleet_control$Sel_start_year[2] <- 1990
  m <- np_build(d)
  start_idx <- 1990 - d$styr + 1
  testthat::expect_true(all(is.na(m$map$mapList$sel_coff_dev[2, 1, 1:10, 1:start_idx])))
  # Values written into those mapped-off cells (a warm start from an older
  # walk fit) neither move the curve nor enter the density.
  ini <- m$initial_params
  ini$sel_coff_dev[2, 1, 1:10, c(1, start_idx)] <- 0.7
  m2 <- np_build(d, inits = ini)
  testthat::expect_equal(m2$quantities$sel_at_age[2, 1, , ], m$quantities$sel_at_age[2, 1, , ])
  testthat::expect_equal(m2$obj$fn(), m$obj$fn())
})

testthat::test_that("a zero deviation SD is refused for the new forms as for NonParametric", {
  testthat::skip_on_cran()
  # Both deviation structures read the sd, so both must refuse a zero.
  for (tv in c("IID", "RandomWalk")) {
    d <- np_data("NonParametricIntegrable", tv)
    d$fleet_control$Time_varying_sel_sd[2] <- 0
    testthat::expect_error(np_build(d), "Time_varying_sel_sd", info = tv)
  }
})

testthat::test_that("NonParametricIntegrable scores the deviates alone and leaves the shape penalty on the base", {
  testthat::skip_on_cran()
  d <- np_data("NonParametricIntegrable", "IID")
  m <- np_build(d)
  cells <- which(!is.na(m$map$mapList$sel_coff_dev[2, 1, , ]))
  testthat::expect_gt(length(cells), 0)
  ini <- m$initial_params
  devs <- seq(-0.5, 0.5, length.out = length(cells))
  ini$sel_coff_dev[2, 1, , ][cells] <- devs
  m2 <- np_build(d, inits = ini)
  r0 <- np_rows(m); r1 <- np_rows(m2)
  testthat::expect_equal(r1[["Non-parametric selectivity"]], r0[["Non-parametric selectivity"]])
  testthat::expect_equal(r1[["Selectivity deviates"]], -sum(stats::dnorm(devs, 0, 0.35, log = TRUE)))
  # The curve does move with the deviates.
  testthat::expect_false(isTRUE(all.equal(m2$quantities$sel_at_age[2, 1, , ], m$quantities$sel_at_age[2, 1, , ])))

  # NonParametric charges the shape penalty on the realized curve, so the same
  # deviates move that row there: the difference these forms exist for.
  d2 <- np_data("NonParametric", "IID")
  n0 <- np_build(d2); ini2 <- n0$initial_params
  ini2$sel_coff_dev[2, 1, , ][cells] <- devs
  n1 <- np_build(d2, inits = ini2)
  testthat::expect_false(isTRUE(all.equal(np_rows(n1)[["Non-parametric selectivity"]],
                                          np_rows(n0)[["Non-parametric selectivity"]])))
})

testthat::test_that("NonParametricIntegrable keeps the base, fixes the first increment and scores the rest", {
  testthat::skip_on_cran()
  d <- np_data("NonParametricIntegrable", "RandomWalk")
  m <- np_build(d)
  ml <- m$map$mapList
  testthat::expect_true(all(!is.na(ml$sel_coff[2, 1, 1:10])))
  testthat::expect_true(all(is.na(ml$sel_coff_dev[2, 1, 1:10, 1])))
  testthat::expect_true(all(!is.na(ml$sel_coff_dev[2, 1, 1:10, 2])))
  nyh <- d$endyr - d$styr + 1
  ini <- m$initial_params
  inc <- seq(-0.4, 0.4, length.out = 10)
  ini$sel_coff_dev[2, 1, 1:10, 3] <- inc
  m2 <- np_build(d, inits = ini)
  s0 <- m$quantities$sel_at_age[2, 1, , ]; s1 <- m2$quantities$sel_at_age[2, 1, , ]
  # A step at year 3 leaves years 1-2 alone and moves every later year alike.
  testthat::expect_equal(s1[, 1:2], s0[, 1:2])
  testthat::expect_false(isTRUE(all.equal(s1[, 3], s0[, 3])))
  testthat::expect_equal(s1[, 3], s1[, nyh])
  # Increments are scored from the year after the start year, on the estimated bins.
  scored <- ini$sel_coff_dev[2, 1, 1:10, 2:nyh]
  testthat::expect_equal(np_rows(m2)[["Selectivity deviates"]], -sum(stats::dnorm(scored, 0, 0.35, log = TRUE)))
  testthat::expect_equal(np_rows(m2)[["Non-parametric selectivity"]], np_rows(m)[["Non-parametric selectivity"]])
})

testthat::test_that("random_sel = TRUE is accepted for the new form and still refused for NonParametric", {
  testthat::skip_on_cran()
  m13 <- np_build(np_data("NonParametricIntegrable", "IID"), random_sel = TRUE)
  testthat::expect_identical(sum(names(m13$obj$par) == "sel_dev_log_sd"), 1L)
  testthat::expect_true("sel_coff_dev" %in% m13$obj$env$.random)
  m13_walk <- np_build(np_data("NonParametricIntegrable", "RandomWalk"), random_sel = TRUE)
  testthat::expect_identical(sum(names(m13_walk$obj$par) == "sel_dev_log_sd"), 1L)
  testthat::expect_error(np_build(np_data("NonParametric", "IID"), random_sel = TRUE),
                         "NonParametricIntegrable")
  testthat::expect_error(np_build(np_data("NonParametric", "RandomWalk"), random_sel = TRUE),
                         "NonParametricIntegrable")
})

# The form supplies the density; Time_varying_sel supplies the structure, so all
# three of its non-parametric modes are legal and a mode with no density is not.
testthat::test_that("NonParametricIntegrable takes Off, IID or RandomWalk", {
  testthat::skip_on_cran()
  # `message =` on expect_no_error() is a regexp that scopes which error is
  # tolerated, not a label, so the assertion is left unscoped and the failing
  # call identifies the mode.
  for (mode in c("Off", "IID", "RandomWalk")) {
    testthat::expect_no_error(np_build(np_data("NonParametricIntegrable", mode)))
  }
  testthat::expect_error(np_build(np_data("NonParametricIntegrable", "Block")),
                         "'Off', 'IID' or 'RandomWalk'")
})

testthat::test_that("coefficients and deviates below the first selected bin are held at 0", {
  testthat::skip_on_cran()
  # They are mapped off, but the curve centres each year by the log mean over
  # every bin, so a value there shifts the whole curve with no density scoring
  # it. `inits` from a fit with a lower Bin_first_selected carry such values:
  # measured before the guard, 0.9 in those cells moved Atka2022's fishery
  # objective by 704 nats (116258.42 -> 115554.32) and year-1 selectivity by 0.21.
  d <- np_data("NonParametricIntegrable", "IID")
  d$fleet_control$Bin_first_selected[2] <- 3
  m   <- np_build(d)
  ini <- m$initial_params
  ini$sel_coff[2, 1, 1:2]      <- 0.9
  ini$sel_coff_dev[2, 1, 1:2, ] <- 0.9
  m2 <- np_build(d, inits = ini)
  testthat::expect_equal(m2$obj$fn(), m$obj$fn())
  testthat::expect_equal(m2$quantities$sel_at_age[2, 1, , ], m$quantities$sel_at_age[2, 1, , ])
  testthat::expect_true(all(m2$initial_params$sel_coff[2, , 1:2] == 0))
  testthat::expect_true(all(m2$initial_params$sel_coff_dev[2, , 1:2, ] == 0))
})


testthat::test_that("above a first selected bin of 1, the scored cells are exactly the estimated ones", {
  testthat::skip_on_cran()
  # build_map() frees bins bin_first_selected:N_sel_bins in 1-based R indexing;
  # the density loops from the template's bin_first_selected, which
  # rearrange_data() made 0-based. The two conventions have to agree: an
  # estimated deviate with no density is a free parameter the integration never
  # sees, and under random_sel the reported SD would absorb it. The Off
  # equivalence tests above cannot catch it, having no deviates at all.
  for (tv in c("IID", "RandomWalk")) {
    iid <- identical(tv, "IID")
    d <- np_data("NonParametricIntegrable", tv)
    d$fleet_control$Bin_first_selected[2] <- 3
    m   <- np_build(d)
    nyh <- d$endyr - d$styr + 1
    yrs <- if (iid) seq_len(nyh) else 2:nyh          # the walk fixes year 1 at 0
    ml  <- m$map$mapList$sel_coff_dev[2, 1, , ]
    testthat::expect_true(all(is.na(ml[1:2, ])), info = tv)         # below the first selected bin
    testthat::expect_true(all(!is.na(ml[3:10, yrs])), info = tv)    # estimated
    # Every estimated cell carries its density: set them all and read the row.
    amp  <- if (iid) 0.4 else 0.05                  # a walk accumulates
    devs <- seq(-amp, amp, length.out = length(3:10) * length(yrs))
    ini  <- m$initial_params
    ini$sel_coff_dev[2, 1, 3:10, yrs] <- devs
    m2 <- np_build(d, inits = ini)
    testthat::expect_equal(np_rows(m2)[["Selectivity deviates"]],
                           -sum(stats::dnorm(devs, 0, 0.35, log = TRUE)), info = tv)
  }
})


testthat::test_that("the deviation SD is estimated, not collapsed, when the deviates are integrated", {
  testthat::skip_on_cran()
  fit <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = np_data("NonParametricIntegrable", "IID"), msmMode = 0, estimateMode = "Hindcast",
    random_sel = TRUE,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
  testthat::expect_true(is.finite(fit$opt$objective))
  testthat::expect_lt(max(abs(fit$obj$gr())), 0.1)
  sd_hat <- exp(fit$estimated_params$sel_dev_log_sd[2])
  testthat::expect_gt(sd_hat, 0.02)
  testthat::expect_lt(sd_hat, 3)
})
