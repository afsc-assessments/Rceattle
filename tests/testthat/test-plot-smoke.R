# Smoke tests for the plotting / output functions (R/7-plot_*.R), the single
# largest untested module. These don't check pixels -- they fit a small model
# once and assert every exported plotter runs without error and draws to a
# device. Output is sent to a throwaway PDF so nothing is written to the repo.
# (Requested in discussion #75 by Cole-Monnahan-NOAA.)
testthat::skip_on_cran()

# Open a throwaway graphics device for the duration of `code` so base-graphics
# plotters don't spew files or pop up windows; always closes it again.
with_null_device <- function(code) {
  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  force(code)
}

# --- single-species plotters --------------------------------------------------
testthat::test_that("single-species plotters run on a fitted model", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("BS2017SS")
  ss <- suppressMessages(suppressWarnings(
    Rceattle::fit_mod(data_list = BS2017SS,
                      estimateMode = 3,  # build only, no optimisation (fast)
                      msmMode = 0,
                      fit_control = Rceattle::fit_control(verbose = 0))
  ))

  with_null_device({
    # time-series family (accept a single fit or a list)
    testthat::expect_no_error(Rceattle::plot_biomass(ss))
    testthat::expect_no_error(Rceattle::plot_ssb(ss))
    testthat::expect_no_error(Rceattle::plot_recruitment(ss))
    testthat::expect_no_error(Rceattle::plot_depletion(ss))
    testthat::expect_no_error(Rceattle::plot_exploitable_biomass(ss))
    testthat::expect_no_error(Rceattle::plot_index(ss))
    testthat::expect_no_error(Rceattle::plot_index(ss, log = TRUE))
    testthat::expect_no_error(Rceattle::plot_catchability(ss))
    testthat::expect_no_error(Rceattle::plot_catchability(ss, log = TRUE))
    testthat::expect_no_error(Rceattle::plot_catch(ss))
    testthat::expect_no_error(Rceattle::plot_f(ss))

    # structural / data plotters
    testthat::expect_no_error(Rceattle::plot_selectivity(ss))
    testthat::expect_no_error(Rceattle::plot_selectivity_vs_maturity(ss))
    testthat::expect_no_error(Rceattle::plot_maturity(ss))
    testthat::expect_no_error(Rceattle::plot_mortality(ss))
    testthat::expect_no_error(Rceattle::plot_stock_recruit(ss))
    testthat::expect_no_error(Rceattle::plot_data(ss))
    testthat::expect_no_error(Rceattle::plot_comp(ss))

    # S3 plot() dispatch
    testthat::expect_no_error(plot(ss, what = "biomass"))
    testthat::expect_no_error(plot(ss, what = "ssb"))
    testthat::expect_no_error(plot(ss, what = "recruitment"))
    testthat::expect_no_error(plot(ss, what = "selectivity"))
    testthat::expect_no_error(plot(ss, what = "data"))
  })
})

# --- multispecies / predation plotters ---------------------------------------
testthat::test_that("multispecies predation plotters run on a fitted model", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("BS2017MS")
  ms <- suppressMessages(suppressWarnings(
    Rceattle::fit_mod(data_list = BS2017MS,
                      estimateMode = 3,
                      msmMode = 1,   # multispecies (predation on)
                      fit_control = Rceattle::fit_control(verbose = 0))
  ))

  with_null_device({
    testthat::expect_no_error(Rceattle::plot_b_eaten(ms))
    testthat::expect_no_error(Rceattle::plot_b_eaten_prop(ms))
    testthat::expect_no_error(Rceattle::plot_ration(ms))
    testthat::expect_no_error(Rceattle::plot_diet_comp(ms))
    testthat::expect_no_error(Rceattle::plot_m_at_age(ms))
    testthat::expect_no_error(Rceattle::plot_m2_at_age_prop(ms))
    testthat::expect_no_error(plot(ms, what = "mortality"))
  })
})

# --- regression: NA observations in projection rows ---------------------------
# `clean_data()` gives every projection year an NA catch, so `.obs` is NA on
# those rows and the lognormal interval mask must not pass an NA to a subscript.
# No special fixture is needed -- plain BS2017SS carries 99 such rows. The bug
# survived because the smoke test above only ever called plot_catch() with the
# default incl_proj = FALSE, which filters them out.
testthat::test_that("catch plotter tolerates NA observations in projection rows", {
  testthat::skip_if_not_installed("TMB")

  data("BS2017SS")
  ss <- suppressMessages(suppressWarnings(
    Rceattle::fit_mod(data_list = BS2017SS,
                      estimateMode = 3,  # build only, no optimisation (fast)
                      msmMode = 0,
                      fit_control = Rceattle::fit_control(verbose = 0))
  ))
  testthat::expect_true(
    any(is.na(ss$data_list$catch_data$Catch[ss$data_list$catch_data$Year > 0])))

  with_null_device({
    testthat::expect_no_error(Rceattle::plot_catch(ss, incl_proj = TRUE))
    testthat::expect_no_error(Rceattle::plot_index(ss, incl_proj = TRUE))
  })
})

testthat::test_that("plot_indexresidual honours incl_proj", {
  testthat::skip_if_not_installed("TMB")

  data("BS2017SS")
  dat <- BS2017SS
  # A projection-year index row, carrying the placeholder a roll-forward writes.
  proj <- dat$index_data[nrow(dat$index_data), ]
  proj$Year <- dat$endyr + 1L
  proj$Observation <- 99999
  dat$index_data <- rbind(dat$index_data, proj)

  ss <- suppressMessages(suppressWarnings(
    Rceattle::fit_mod(data_list = dat, estimateMode = 3, msmMode = 0,
                      fit_control = Rceattle::fit_control(verbose = 0))))

  yrs_of <- function(...) sort(unique(Rceattle::plot_indexresidual(...)$data$Year))
  with_null_device({
    testthat::expect_false(any(yrs_of(ss) > ss$data_list$endyr))
    testthat::expect_true(any(yrs_of(ss, incl_proj = TRUE) > ss$data_list$endyr))
  })
})


# --- the whole family at minage != 1 ------------------------------------------
# `nages` is a COUNT of age bins, not the oldest age: ages run
# minage .. minage + nages - 1, and the array's 3rd dimension is a 1-based bin
# index, so age `a` sits at index `a - minage + 1`. Every bundled dataset and
# all three live assessments have minage = 1, which makes the two coincide -- so
# a bin index read as an age passes every other test in this suite and every
# real model, and only surfaces on the next minage != 1 species.
#
# test-plot-predation-args.R pins the five predation plotters against that
# specifically. This sweeps the REST of the exported family, so the class is
# closed rather than the five instances: a plotter that indexes an age array
# directly will either error or silently mislabel here.
testthat::test_that("every exported plotter runs on a minage != 1 model", {
  testthat::skip_if_not_installed("TMB")

  data("BS2017SS")
  dat <- BS2017SS
  dat$minage <- rep(3L, dat$nspp)   # ages 3 .. 3 + nages - 1

  fit <- suppressMessages(suppressWarnings(
    Rceattle::fit_mod(data_list = dat, estimateMode = 3, msmMode = 0,
                      fit_control = Rceattle::fit_control(verbose = 0))))

  testthat::expect_equal(unique(fit$data_list$minage), 3L)

  plotters <- c("plot_biomass", "plot_ssb", "plot_recruitment", "plot_depletion",
                "plot_exploitable_biomass", "plot_index", "plot_catch", "plot_f",
                "plot_selectivity", "plot_maturity", "plot_mortality",
                "plot_stock_recruit", "plot_comp", "plot_data",
                "plot_indexresidual", "plot_selectivity_vs_maturity",
                "plot_logindex")

  # A name that is not an exported plotter FAILS. Skipping it instead -- which
  # this did -- let two misspellings sit here unnoticed: `plot_catchresidual`,
  # which no version of the package has ever had, and `plot_sel_vs_mat` for
  # `plot_selectivity_vs_maturity`. The loop then ran 15 of the 17 it listed
  # and still passed its `> 10` floor.
  missing <- plotters[!vapply(plotters, exists, logical(1),
                              envir = asNamespace("Rceattle"))]
  testthat::expect_equal(missing, character(0))

  failed <- character(0)
  ran <- 0L
  with_null_device({
    for (fn in plotters) {
      f <- get(fn, envir = asNamespace("Rceattle"))
      ran <- ran + 1L
      tryCatch({
        p <- suppressMessages(suppressWarnings(f(fit)))
        # A ggplot that assembles but cannot build is not a pass.
        if (inherits(p, "ggplot")) invisible(ggplot2::ggplot_build(p))
      }, error = function(e) {
        failed <<- c(failed, paste0(fn, "(): ", conditionMessage(e)))
      })
    }
  })

  # Exact, not a floor: a floor with slack is what hid the two missing names.
  testthat::expect_equal(ran, length(plotters))
  testthat::expect_equal(failed, character(0))
})
