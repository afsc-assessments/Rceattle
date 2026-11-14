# Integration test (fits a CEATTLE model): skipped on CRAN to keep R CMD check
# fast; the coverage workflow runs it in full (NOT_CRAN=true). See README.md.
testthat::skip_on_cran()

# A fit averaged with itself must reproduce its own standard errors, since
# model_average() draws each fit's parameters from its asymptotic covariance.
# Before 5.67.6 it drew from the precision matrix without random effects (log
# SSB SE 1,320-4,380 times the fit's here) and, with them, set the first
# parameter's variance to 100 (pollock's log SSB SE up to 383 times).
testthat::test_that("model_average() of a fit with itself reproduces its standard errors", {
  testthat::skip_if_not_installed("TMB")

  data(BS2017SS, envir = environment())
  d <- BS2017SS
  d$projyr <- 2020

  inits <- NULL
  for (re in c(FALSE, TRUE)) {
    fit <- suppressMessages(suppressWarnings(fit_mod(
      data_list = d, inits = inits, file = NULL, estimateMode = 0,
      random_rec = re, msmMode = 0,
      fit_control = fit_control(phase = FALSE, getsd = TRUE, verbose = 0))))
    inits <- fit$estimated_params

    set.seed(1)
    avg <- model_average(list(fit, fit), uncertainty = TRUE, nboot = 400)

    # log_ssb lists every species for year 1, then year 2, ..., so the
    # hindcast is its first nspp * nyrs_hind entries.
    nyrs_hind <- fit$data_list$endyr - fit$data_list$styr + 1
    rows <- which(names(fit$sdrep$value) == "log_ssb")
    rows <- rows[seq_len(fit$data_list$nspp * nyrs_hind)]
    ratio <- avg$sdrep$sd[rows] / fit$sdrep$sd[rows]

    # 400 draws, and SSB's curvature in the parameters, let single years stray
    # from the delta-method SE (up to about 1.5 times); the median should not.
    label <- paste0("random_rec = ", re)
    testthat::expect_true(all(is.finite(ratio)), label = label)
    testthat::expect_true(all(ratio > 0.5 & ratio < 2), label = label)
    testthat::expect_lt(abs(stats::median(ratio) - 1), 0.1, label = label)
  }
})
