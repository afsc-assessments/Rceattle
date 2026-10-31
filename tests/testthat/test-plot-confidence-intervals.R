# Confidence-interval plumbing for the time-series plotters.
#
# Every other plotting test builds its fixture with estimateMode = 3, which
# skips sdreport() and leaves `sdrep` NULL -- so the whole `add_ci = TRUE`
# branch never runs. This file is the one that fits with getsd = TRUE and
# exercises it.
#
# What it pins: that the interval comes from the model's own standard errors,
# that the log-scale series use the reported log sd while the depletions fall
# back to the delta-method identity sd(log x) = sd(x)/x, and that a series the
# template never ADREPORTs draws no interval and says so.
testthat::skip_on_cran()

fit_with_sd <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    testthat::skip_if_not_installed("TMB")
    data("BS2017SS", package = "Rceattle")
    cached <<- suppressMessages(suppressWarnings(
      Rceattle::fit_mod(data_list = BS2017SS,
                        estimateMode = 1,      # hindcast only
                        msmMode = 0,
                        phase = FALSE,
                        fit_control = Rceattle::fit_control(getsd = TRUE,
                                                            verbose = 0))))
    cached
  }
})


testthat::test_that("add_ci draws an interval built from the model's standard errors", {
  fit <- fit_with_sd()
  testthat::skip_if(is.null(fit$sdrep), "fit carries no sdreport")

  p <- Rceattle::plot_biomass(fit, add_ci = TRUE)
  testthat::expect_true(all(c("lower95", "upper95") %in% names(p$data)))
  testthat::expect_false(all(is.na(p$data$lower95)))
  # A 95% interval brackets the point estimate and is right-skewed on the log
  # scale, so it cannot cross zero for a positive series.
  ok <- !is.na(p$data$lower95)
  testthat::expect_true(all(p$data$lower95[ok] <= p$data$value[ok]))
  testthat::expect_true(all(p$data$upper95[ok] >= p$data$value[ok]))
  testthat::expect_true(all(p$data$lower95[ok] > 0))

  # No interval is drawn without the flag.
  testthat::expect_true(all(is.na(Rceattle::plot_biomass(fit)$data$lower95)))
})

testthat::test_that("the interval matches exp(log(x) +/- z * sd_log) from sdrep", {
  # The figure must agree with the model, not merely look plausible. biomass is
  # ADREPORTed on both scales, so the log sd is read directly.
  fit <- fit_with_sd()
  testthat::skip_if(is.null(fit$sdrep), "fit carries no sdreport")

  nspp  <- fit$data_list$nspp
  nyrs  <- fit$data_list$endyr - fit$data_list$styr + 1L
  n_tot <- ncol(fit$quantities$biomass) * nspp
  sd_log <- .rce_series_sd(fit, "log_biomass", nyrs * nspp, n_tot)
  testthat::skip_if(is.null(sd_log), "log_biomass not reported by this fit")

  p  <- Rceattle::plot_biomass(fit, add_ci = TRUE)
  sp <- fit$data_list$spnames[1]
  d  <- p$data[p$data$Species == sp, ]
  d  <- d[order(d$Year), ]
  # Column-major: species varies fastest, so species 1 is every nspp-th value.
  sd_sp1 <- sd_log[seq(1, by = nspp, length.out = nyrs)]

  testthat::expect_equal(d$upper95,
                         exp(log(d$value) + stats::qnorm(0.975) * sd_sp1),
                         tolerance = 1e-8)
  testthat::expect_equal(d$lower95,
                         exp(log(d$value) - stats::qnorm(0.975) * sd_sp1),
                         tolerance = 1e-8)
})

testthat::test_that("a series with no log_ counterpart recovers sd_log as sd/x", {
  # ssb_depletion cannot be ADREPORTed on the log scale (it divides by SB0), so
  # the plotter uses the delta-method identity instead. The two routes must give
  # the same interval, which is what makes the fallback safe.
  fit <- fit_with_sd()
  testthat::skip_if(is.null(fit$sdrep), "fit carries no sdreport")

  nspp  <- fit$data_list$nspp
  nyrs  <- fit$data_list$endyr - fit$data_list$styr + 1L
  n_tot <- ncol(fit$quantities$ssb_depletion) * nspp
  sd_nat <- .rce_series_sd(fit, "ssb_depletion", nyrs * nspp, n_tot)
  testthat::skip_if(is.null(sd_nat), "ssb_depletion not reported by this fit")
  testthat::expect_null(.rce_series_sd(fit, "log_ssb_depletion",
                                       nyrs * nspp, n_tot))

  p  <- Rceattle::plot_ssb_depletion(fit, add_ci = TRUE)
  sp <- fit$data_list$spnames[1]
  d  <- p$data[p$data$Species == sp, ]
  d  <- d[order(d$Year), ]
  sd_sp1 <- sd_nat[seq(1, by = nspp, length.out = nyrs)]

  ok <- !is.na(d$upper95) & d$value > 0
  testthat::expect_equal(
    d$upper95[ok],
    exp(log(d$value[ok]) + stats::qnorm(0.975) * (sd_sp1[ok] / d$value[ok])),
    tolerance = 1e-8)
})

testthat::test_that("a REPORT-only series draws no interval and stays quiet", {
  # F_spp is not ADREPORTed, so warning about it would fire on every fit.
  fit <- fit_with_sd()
  testthat::skip_if(is.null(fit$sdrep), "fit carries no sdreport")

  testthat::expect_false(.rce_quantity_adreport("F_spp"))
  p <- testthat::expect_silent(
    Rceattle::plot_timeseries(fit, output = "F_spp", add_ci = TRUE))
  testthat::expect_true(all(is.na(p$data$lower95)))
})

testthat::test_that("a fit without sdreport warns once and draws no interval", {
  testthat::skip_if_not_installed("TMB")
  data("BS2017SS", package = "Rceattle")
  no_sd <- suppressMessages(suppressWarnings(
    Rceattle::fit_mod(data_list = BS2017SS, estimateMode = 3, msmMode = 0,
                      fit_control = Rceattle::fit_control(verbose = 0))))

  testthat::expect_warning(p <- Rceattle::plot_ssb(no_sd, add_ci = TRUE),
                           "no standard errors")
  testthat::expect_true(all(is.na(p$data$lower95)))
})

testthat::test_that("a peel whose sd rebuild failed draws no band and says why", {
  fit <- fit_with_sd()
  testthat::skip_if(is.null(fit$sdrep), "fit carries no sdreport")
  attr(fit, "sd_report_error") <- "peel 1 (through 2015): indefinite Hessian"

  testthat::expect_warning(p <- Rceattle::plot_ssb(fit, add_ci = TRUE),
                           "peel 1 \\(through 2015\\): indefinite Hessian")
  testthat::expect_true(all(is.na(p$data$lower95)))
  testthat::expect_null(Rceattle:::.rce_sel_se(fit))
  tab <- as.data.frame(fit)
  testthat::expect_true(all(is.na(tab$se)))
  testthat::expect_null(stats::vcov(fit))
  testthat::expect_true(any(grepl("peel 1 \\(through 2015\\)",
                                  utils::capture.output(print(summary(fit))))))
})

testthat::test_that("the standard errors read are the ones the block actually holds", {
  # The series are flattened column-major over the whole modelled period, so a
  # hindcast plot takes a leading slice. Taking that slice without checking the
  # block length would return another cell's standard errors.
  fit <- fit_with_sd()
  testthat::skip_if(is.null(fit$sdrep), "fit carries no sdreport")

  nspp  <- fit$data_list$nspp
  n_tot <- ncol(fit$quantities$biomass) * nspp
  testthat::expect_equal(sum(names(fit$sdrep$value) == "biomass"), n_tot)
  # A wrong claim about the block length yields nothing rather than a
  # misaligned interval.
  testthat::expect_null(.rce_series_sd(fit, "biomass", 10L, n_tot + 1L))
  testthat::expect_length(.rce_series_sd(fit, "biomass", 10L, n_tot), 10L)
})
