# Accuracy tests for the ggplot plotters: each returns a ggplot whose `$data`
# holds exactly the model's source quantities. We cross-check against
# as.data.frame.Rceattle() -- the canonical tidy output of the same quantities
# -- so the plot and the public data extractor are guaranteed consistent.
# This is the regression guarantee for the base->ggplot rewrite ("accuracy is
# key"), complementing the render-only smoke tests.
testthat::skip_on_cran()

# --- time-series family (plot_timeseries engine) ------------------------------
testthat::test_that("time-series plotters match as.data.frame() quantities", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("BS2017SS")
  ss <- suppressMessages(suppressWarnings(
    Rceattle::fit_mod(data_list = BS2017SS, estimateMode = 3, msmMode = 0,
                      fit_control = Rceattle::fit_control(verbose = 0))
  ))
  sp1 <- ss$data_list$spnames[1]
  ad <- as.data.frame(ss, which = c("biomass", "ssb", "R",
                                    "exploitable_biomass", "ssb_depletion",
                                    "F_spp"))

  # Compare a plotter's species-1 trajectory to the same quantity from
  # as.data.frame(), aligned on the years the plot actually shows. `scale` is
  # the display rescaling the plotter applies: numbers-at-age are in thousands
  # and weight-at-age in kg, so every biomass series is in mt (1e6 -> million
  # mt) and recruitment is in thousands of fish (1e3 -> millions of recruits).
  # Depletion and F are ratios and are plotted unscaled.
  check <- function(p, quantity, scale = 1) {
    testthat::expect_s3_class(p, "ggplot")
    pd <- p$data[p$data$Species == sp1, c("Year", "value")]
    pd <- pd[order(pd$Year), ]
    # Guard against an empty selection: if the Species filter matches nothing,
    # both pd$value and src are numeric(0) and expect_equal() passes vacuously.
    testthat::expect_gt(nrow(pd), 0)
    sd <- ad[ad$quantity == quantity & ad$species == sp1, c("year", "value")]
    src <- sd$value[match(pd$Year, sd$year)] / scale
    testthat::expect_equal(pd$value, src)
  }

  check(Rceattle::plot_biomass(ss),             "biomass",             scale = 1e6)
  check(Rceattle::plot_ssb(ss),                 "ssb",                 scale = 1e6)
  check(Rceattle::plot_recruitment(ss),         "R",                   scale = 1e3)
  check(Rceattle::plot_exploitable_biomass(ss), "exploitable_biomass", scale = 1e6)
  check(Rceattle::plot_ssb_depletion(ss),       "ssb_depletion")
  check(Rceattle::plot_f(ss),                   "F_spp")  # wraps plot_timeseries
})

testthat::test_that("minyr and maxyr narrow what the timeseries plotters draw", {
  # `minyr` was the argument the standardization was meant to make work, and it
  # is the one with no coverage: a plotter that ignores it still passes every
  # accuracy check above, because those align on the years the plot shows.
  testthat::skip_if_not_installed("TMB")

  data("BS2017SS")
  ss <- suppressMessages(suppressWarnings(
    Rceattle::fit_mod(data_list = BS2017SS, estimateMode = 3, msmMode = 0,
                      fit_control = Rceattle::fit_control(verbose = 0))
  ))
  styr  <- ss$data_list$styr
  endyr <- ss$data_list$endyr
  cut   <- styr + 5L

  for (fn in c("plot_biomass", "plot_ssb", "plot_recruitment", "plot_f",
               "plot_ssb_depletion", "plot_exploitable_biomass")) {
    f <- getExportedValue("Rceattle", fn)
    testthat::expect_equal(min(f(ss, minyr = cut)$data$Year), cut, label = fn)
    testthat::expect_equal(max(f(ss, maxyr = cut)$data$Year), cut, label = fn)
    testthat::expect_equal(range(f(ss, minyr = cut, maxyr = cut + 3L)$data$Year),
                           c(cut, cut + 3L), label = fn)
    # Unset, the window is the model's own.
    testthat::expect_equal(range(f(ss)$data$Year), c(styr, endyr), label = fn)
  }

  # Narrowing the data, not just the axis: the panel's y scale must follow the
  # window. Clipping the axis alone leaves it trained on the hidden years.
  full <- ggplot2::ggplot_build(Rceattle::plot_biomass(ss))
  cutp <- ggplot2::ggplot_build(Rceattle::plot_biomass(ss, minyr = endyr - 2L))
  testthat::expect_false(isTRUE(all.equal(full$layout$panel_params[[1]]$y.range,
                                          cutp$layout$panel_params[[1]]$y.range)))

  testthat::expect_error(Rceattle::plot_biomass(ss, minyr = endyr + 50L),
                         "No years left")
})

# The per-species reference-point frame a plot carries, ordered by species.
# plot_f() and plot_ssb_depletion() build it through the same helper, so one
# accessor reads both.
ref_lines <- function(p) {
  d <- Filter(function(l) "target" %in% names(l$data), p$layers)[[1]]$data
  d[order(as.character(d$Species)), ]
}

testthat::test_that("plot_f keys its reference lines to the species it drew", {
  # plot_f() is hand-written rather than a .ts_wrapper() child, and indexed
  # Ftarget and the facet labels with the raw `species`. That works for indices
  # and gives an NA facet key for a name, which puts the reference lines in the
  # wrong panel -- on the same argument NEWS advertises as newly accepting names.
  testthat::skip_if_not_installed("TMB")

  data("BS2017SS")
  ss <- suppressMessages(suppressWarnings(
    Rceattle::fit_mod(data_list = BS2017SS, estimateMode = 3, msmMode = 0,
                      fit_control = Rceattle::fit_control(verbose = 0))
  ))
  nm <- ss$data_list$spnames

  by_index <- ref_lines(Rceattle::plot_f(ss, species = c(1, 3)))
  by_name  <- ref_lines(Rceattle::plot_f(ss, species = nm[c(1, 3)]))
  testthat::expect_false(anyNA(by_name$Species))
  testthat::expect_equal(by_name, by_index, ignore_attr = TRUE)
  testthat::expect_setequal(as.character(by_index$Species), nm[c(1, 3)])
  testthat::expect_equal(by_index$target,
                         as.numeric(ss$quantities$Ftarget[c(1, 3)]))

  # `spnames` relabels the panels; the reference lines must follow.
  lab <- c("A", "B", "C")
  testthat::expect_setequal(
    as.character(
      ref_lines(Rceattle::plot_f(ss, species = 2, spnames = lab))$Species),
    "B")

  # lty was missing from plot_f's formals, so it stopped with `unused argument`.
  testthat::expect_equal(
    unique(ggplot2::ggplot_build(Rceattle::plot_f(ss, lty = 2))$data[[1]]$linetype),
    2)
})

# --- fit diagnostics (index/catch) --------------------------------------------
# The fit plots draw the MEDIAN of the fitted observation distribution, so the
# expected value here is hat * exp(-b*sd^2/2) rather than hat; before 5.64.3
# this pinned the mean. Recomputed from the reported sd and flag, not from the
# plotter's own helper.
testthat::test_that("index/catch fit plots plot the model's predictions", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("BS2017SS")
  ss <- suppressMessages(suppressWarnings(
    Rceattle::fit_mod(data_list = BS2017SS, estimateMode = 3, msmMode = 0,
                      fit_control = Rceattle::fit_control(verbose = 0))
  ))
  fc <- ss$data_list$fleet_control

  # survey index: plotted Predicted == quantities$index_hat on the survey rows
  pidx <- Rceattle::plot_index(ss)
  testthat::expect_s3_class(pidx, "ggplot")
  idx <- ss$data_list$index_data
  idx$hat <- ss$quantities$index_hat
  idx$sd  <- as.numeric(ss$quantities$index_sd)
  surv <- fc$Fleet_code[fc$Fleet_type == "Survey"]
  idx <- idx[idx$Fleet_code %in% surv & idx$Year <= ss$data_list$endyr, ]
  pd  <- pidx$data[order(pidx$data$Fleet, pidx$data$Year), ]
  src <- idx[order(idx$Fleet_code, abs(idx$Year)), ]
  ba  <- 1   # bias_adjust_obs defaults to TRUE
  testthat::expect_equal(pd$Predicted, src$hat * exp(-ba * src$sd^2 / 2))

  # fishery catch: plotted Predicted == quantities$catch_hat on the fishery rows
  pcat <- Rceattle::plot_catch(ss)
  testthat::expect_s3_class(pcat, "ggplot")
  cat_d <- ss$data_list$catch_data
  cat_d$hat <- ss$quantities$catch_hat
  cat_d$sd  <- as.numeric(ss$quantities$catch_sd)
  fish <- fc$Fleet_code[fc$Fleet_type == "Fishery"]
  cat_d <- cat_d[cat_d$Fleet_code %in% fish & cat_d$Year <= ss$data_list$endyr, ]
  pdc <- pcat$data[order(pcat$data$Fleet, pcat$data$Year), ]
  srcc <- cat_d[order(cat_d$Fleet_code, abs(cat_d$Year)), ]
  testthat::expect_equal(pdc$Predicted, srcc$hat * exp(-ba * srcc$sd^2 / 2))
})

# --- maturity -----------------------------------------------------------------
testthat::test_that("plot_maturity plots data_list$maturity at age", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("BS2017SS")
  ss <- suppressMessages(suppressWarnings(
    Rceattle::fit_mod(data_list = BS2017SS, estimateMode = 3, msmMode = 0,
                      fit_control = Rceattle::fit_control(verbose = 0))
  ))
  p <- Rceattle::plot_maturity(ss)
  testthat::expect_s3_class(p, "ggplot")
  mat <- ss$data_list$maturity[, -1]
  sp1 <- ss$data_list$spnames[1]
  d <- p$data[p$data$Species == sp1, ]
  d <- d[order(d$Age), ]
  src <- as.numeric(mat[1, seq_len(ss$data_list$nages[1])])
  testthat::expect_equal(d$Maturity, src)
})

# --- stock-recruit ------------------------------------------------------------
testthat::test_that("plot_stock_recruit plots SSB vs recruitment", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("BS2017SS")
  ss <- suppressMessages(suppressWarnings(
    Rceattle::fit_mod(data_list = BS2017SS, estimateMode = 3, msmMode = 0,
                      fit_control = Rceattle::fit_control(verbose = 0))
  ))
  p <- Rceattle::plot_stock_recruit(ss)
  testthat::expect_s3_class(p, "ggplot")
  sp1 <- ss$data_list$spnames[1]
  nyr <- ss$data_list$endyr - ss$data_list$styr + 1L
  d <- p$data[p$data$Species == sp1, ]   # points layer data (SSB, R)
  testthat::expect_equal(sort(d$SSB),
                         sort(as.numeric(ss$quantities$ssb[1, seq_len(nyr)] / 1e6)))
  # R is carried in thousands of fish and plotted in millions of recruits, so
  # this divisor is 1e3 where SSB's (mt -> million mt) is 1e6.
  testthat::expect_equal(sort(d$R),
                         sort(as.numeric(ss$quantities$R[1, seq_len(nyr)] / 1e3)))
})

# --- selectivity --------------------------------------------------------------
testthat::test_that("plot_selectivity plots sel_at_age", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("BS2017SS")
  ss <- suppressMessages(suppressWarnings(
    Rceattle::fit_mod(data_list = BS2017SS, estimateMode = 3, msmMode = 0,
                      fit_control = Rceattle::fit_control(verbose = 0))
  ))
  p <- Rceattle::plot_selectivity(ss)
  testthat::expect_s3_class(p, "ggplot")

  fc <- ss$data_list$fleet_control
  nyrshind <- ss$data_list$endyr - ss$data_list$styr + 1L
  nages1 <- ss$data_list$nages[fc$Species[fc$Fleet_code == 1]]
  flt1 <- as.character(fc$Fleet_name[fc$Fleet_code == 1])
  # `Bin` is the fleet's own selectivity dimension -- an age here, a length-bin
  # ordinal for a length-based fleet -- so it is not named `Age`.
  d <- p$data[p$data$Fleet == flt1 & p$data$Year == ss$data_list$endyr, ]
  d <- d[order(d$Bin), ]
  src <- as.numeric(ss$quantities$sel_at_age[1, 1, seq_len(nages1), nyrshind])
  testthat::expect_equal(d$Selectivity, src)
  testthat::expect_setequal(unique(d$Dimension), "Age")
  testthat::expect_equal(d$Bin, seq_len(nages1) - 1L + ss$data_list$minage[1])
})

# --- mortality-at-age ---------------------------------------------------------
testthat::test_that("the mortality and predation plotters draw their source arrays", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("BS2017MS")
  ms <- suppressMessages(suppressWarnings(
    Rceattle::fit_mod(data_list = BS2017MS, estimateMode = 3, msmMode = 1,
                      fit_control = Rceattle::fit_control(verbose = 0))
  ))
  p <- Rceattle::plot_mortality(ms)  # default M2 = TRUE -> M2_at_age
  testthat::expect_s3_class(p, "ggplot")
  sp1 <- ms$data_list$spnames[1]
  nages1 <- ms$data_list$nages[1]
  d <- p$data[p$data$Species == sp1 & p$data$Year == ms$data_list$styr, ]
  d <- d[order(d$Age), ]
  src <- as.numeric(ms$quantities$M2_at_age[1, 1, seq_len(nages1), 1])
  testthat::expect_equal(d$M, src)
  testthat::expect_equal(p$labels$y, "M2")

  # M2 = FALSE draws M1, not the sum. The help said "total M (M1 + M2)" for
  # this branch, and nothing here contradicted it.
  p1 <- Rceattle::plot_mortality(ms, M2 = FALSE)
  d1 <- p1$data[p1$data$Species == sp1 & p1$data$Year == ms$data_list$styr, ]
  d1 <- d1[order(d1$Age), ]
  testthat::expect_equal(
    d1$M, as.numeric(ms$quantities$M1_at_age[1, 1, seq_len(nages1), 1]))
  testthat::expect_equal(p1$labels$y, "M1")
  # The two branches draw different arrays -- neither is the sum, which is what
  # plot_m_at_age() draws from M_at_age. (M_at_age is defined as M1 + M2 in the
  # template, but the REPORTed arrays only satisfy that once the multispecies
  # iteration converges: M2_at_age is recomputed after M_at_age each pass, so on
  # an estimateMode = 3 build it lags by one. Not an identity to assert here.)
  testthat::expect_false(isTRUE(all.equal(d$M, d1$M)))

  # type: "heatmap" tiles the same series; the default and anything else lines.
  geom1 <- function(p) class(p$layers[[1]]$geom)[1]
  testthat::expect_equal(geom1(Rceattle::plot_mortality(ms, type = "heatmap")),
                         "GeomTile")
  testthat::expect_equal(geom1(Rceattle::plot_mortality(ms, type = 0)),
                         "GeomTile")
  testthat::expect_equal(geom1(p), "GeomLine")

  # maxage names an age, not a count of bins. At minage = 1 the two agree; a
  # shifted copy is what tells them apart -- ages 2..5, not bins 1..5.
  testthat::expect_equal(max(Rceattle::plot_mortality(ms, maxage = 5)$data$Age), 5)
  shifted <- ms
  shifted$data_list$minage <- rep(2L, ms$data_list$nspp)
  ds <- Rceattle::plot_mortality(shifted, maxage = 5)$data
  testthat::expect_equal(max(ds$Age), 5)
  testthat::expect_equal(min(ds$Age), 2)
  testthat::expect_error(Rceattle::plot_mortality(shifted, maxage = 1),
                         "below the youngest age")

  # plot_m_at_age: total M at a fixed age over time == M_at_age[sp, sex, age, ]
  pa <- Rceattle::plot_m_at_age(ms, age = 1)
  testthat::expect_s3_class(pa, "ggplot")
  nyr <- ms$data_list$endyr - ms$data_list$styr + 1L
  d2 <- pa$data[pa$data$Species == sp1, ]
  d2 <- d2[order(d2$Year), ]
  src2 <- as.numeric(ms$quantities$M_at_age[1, 1, 1, seq_len(nyr)])
  testthat::expect_equal(d2$M[seq_len(nyr)], src2)

  # plot_b_eaten: total biomass eaten as prey == sum of B_eaten_as_prey, in
  # million mt (the model reports it in mt).
  pb <- Rceattle::plot_b_eaten(ms)
  testthat::expect_s3_class(pb, "ggplot")
  be <- ms$quantities$B_eaten_as_prey
  tot1 <- apply(be[, , , seq_len(nyr), drop = FALSE], c(1, 4), sum)[1, ]
  db <- pb$data[pb$data$Species == sp1, ]
  db <- db[order(db$Year), ]
  testthat::expect_equal(db$value[seq_len(nyr)], as.numeric(tot1) / 1e6)

  # ... on the same scale as plot_b_eaten_prop(), which splits the same
  # quantity by predator. Summing that over predators recovers this.
  pp <- Rceattle::plot_b_eaten_prop(ms)
  by_pred <- stats::aggregate(value ~ Year, data = pp$data[pp$data$Prey == sp1, ],
                              FUN = sum)
  by_pred <- by_pred[order(by_pred$Year), ]
  testthat::expect_equal(by_pred$value[seq_len(nyr)],
                         db$value[seq_len(nyr)], tolerance = 1e-8)
})

testthat::test_that("depletion reference lines come from one model, per species", {
  # plot_timeseries() collected Ptarget / Plimit into a models x species matrix
  # and subset it with the species indices, which flattens column-major. On a
  # two-model overlay whose models carry different Ptarget, species 2's line was
  # drawn from model 2's species 1. Models in one figure normally share their
  # reference points, and then every row is identical, which is why it survived.
  testthat::skip_if_not_installed("TMB")

  data("BS2017SS")
  a <- suppressMessages(suppressWarnings(
    Rceattle::fit_mod(data_list = BS2017SS, estimateMode = 3, msmMode = 0,
                      fit_control = Rceattle::fit_control(verbose = 0))
  ))
  b <- a
  b$data_list$Ptarget <- c(0.11, 0.22, 0.33)
  b$data_list$Plimit  <- c(0.011, 0.022, 0.033)

  # One horizontal line per panel can only show one model's value: the first.
  for (mods in list(list(a, b), list(b, a))) {
    d <- ref_lines(Rceattle::plot_ssb_depletion(mods, model_names = c("x", "y")))
    d <- d[match(a$data_list$spnames, as.character(d$Species)), ]
    testthat::expect_equal(d$target, as.numeric(mods[[1]]$data_list$Ptarget))
    testthat::expect_equal(d$limit, as.numeric(mods[[1]]$data_list$Plimit))
  }

  # A single model is unaffected, which is every bundled example.
  d1 <- ref_lines(Rceattle::plot_ssb_depletion(a))
  testthat::expect_equal(d1$target[order(as.character(d1$Species))],
                         as.numeric(a$data_list$Ptarget)[
                           order(a$data_list$spnames)])
})

# --- catchability -------------------------------------------------------------
# plot_catchability() reads `quantities$index_q`, the realized catchability the
# model scaled each survey by, whichever route set it. There is no
# as.data.frame() quantity to cross-check against, so the source is index_q
# itself; what this pins is the fleet-to-row and column-to-year alignment, which
# is where a fleet-indexed array goes wrong (see the Fleet_code invariant).
testthat::test_that("plot_catchability matches quantities$index_q", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("BS2017SS")
  ss <- suppressMessages(suppressWarnings(
    Rceattle::fit_mod(data_list = BS2017SS, estimateMode = 3, msmMode = 0,
                      fit_control = Rceattle::fit_control(verbose = 0))
  ))

  p  <- Rceattle::plot_catchability(ss)
  testthat::expect_s3_class(p, "ggplot")
  fc <- ss$data_list$fleet_control
  q  <- ss$quantities$index_q

  # Only fleets carrying index_data are drawn -- a fishery with a CPUE series as
  # much as a survey -- and every one of them is.
  drawn <- levels(p$data$Fleet)
  expected <- fc$Fleet_name[match(
    Rceattle:::.fleets_with_index(ss$data_list, fitted_only = FALSE),
    fc$Fleet_code)]
  testthat::expect_equal(drawn, expected)

  for (nm in drawn) {
    pd <- p$data[p$data$Fleet == nm, ]
    pd <- pd[order(pd$Year), ]
    testthat::expect_equal(pd$value, unname(q[match(nm, fc$Fleet_name), ]))
    # Hindcast years only: the model does not project catchability.
    testthat::expect_equal(range(pd$Year),
                           c(ss$data_list$styr, ss$data_list$endyr))
  }

  # log = TRUE plots the same series on the scale a lognormal survey
  # likelihood and the q linkages work on.
  pl <- Rceattle::plot_catchability(ss, log = TRUE)
  testthat::expect_equal(pl$data$value, log(p$data$value))
})
