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
  fits <- list()
  for (re in c(FALSE, TRUE)) {
    fit <- suppressMessages(suppressWarnings(fit_mod(
      data_list = d, inits = inits, file = NULL, estimateMode = 0,
      random_rec = re, msmMode = 0,
      fit_control = fit_control(phase = FALSE, getsd = TRUE, verbose = 0))))
    inits <- fit$estimated_params
    fits[[length(fits) + 1]] <- fit

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
    # The average starts as a copy of the fit, so the ratio alone would also
    # pass if the SEs were never rewritten; they must be the draws' own SD.
    testthat::expect_equal(
      avg$sdrep$sd[rows],
      as.vector(apply(log(avg$asymptotic_samples$ssb[, seq_len(nyrs_hind), ]),
                      c(1, 2), stats::sd)),
      label = label)
  }

  # A fixed-effects fit has no joint precision matrix, but mixing the two kinds
  # is the cause to name.
  testthat::expect_error(
    model_average(fits, uncertainty = TRUE, nboot = 10),
    "purely fixed and random effects")
  # A random-effects fit without an sdreport reads as fixed effects, so the
  # missing sdreport must be named before the mixed-list check runs.
  re_nosd <- fits[[2]]
  re_nosd$sdrep <- NULL
  testthat::expect_error(
    model_average(list(fits[[2]], re_nosd), uncertainty = TRUE, nboot = 10),
    "no sdreport")
})


# Under any HCR but NoFishing, without projection_uncertainty, fit$obj is the
# projection's and frees only the rule's F: on BS2017SS under NPFMC every
# averaged hindcast log-SSB SE was 0. Under ConstantF, whose sdrep is the
# hindcast's, model_average() stopped on "incompatible arguments"; under
# estimateMode = "Projection" with projection_uncertainty, on "infinite or
# missing values in 'x'"; a getsd = FALSE fit on "'data' must be of a vector
# type, was 'NULL'".
testthat::test_that("model_average() refuses a fit it cannot draw the hindcast from", {
  testthat::skip_if_not_installed("TMB")

  data(BS2017SS, envir = environment())
  d <- BS2017SS
  d$projyr <- 2020
  d$fleet_control$Proj_F_proportion <- 1  # an HCR needs projected fishing
  fit <- function(hcr = build_hcr(), inits = NULL, estimateMode = 0, ...) {
    suppressMessages(suppressWarnings(fit_mod(
      data_list = d, inits = inits, file = NULL, estimateMode = estimateMode,
      msmMode = 0, HCR = hcr,
      fit_control = fit_control(phase = FALSE, verbose = 0, ...))))
  }
  avg <- function(m1, m2) model_average(list(m1, m2), uncertainty = TRUE, nboot = 400)

  hind   <- fit()
  inits  <- hind$estimated_params
  npfmc  <- fit(build_hcr(HCR = "NPFMC"), inits)
  constf <- fit(build_hcr(HCR = "ConstantF", Ftarget = 0.1), inits)
  # Under "Projection", projection_uncertainty frees `dummy` beside F, not the hindcast.
  proj   <- fit(build_hcr(HCR = "NPFMC"), inits, estimateMode = 2,
                projection_uncertainty = TRUE)
  refused <- "free no hindcast parameter"
  testthat::expect_error(avg(npfmc, npfmc), refused)
  testthat::expect_error(avg(constf, constf), refused)
  testthat::expect_error(avg(proj, proj), refused)
  # The cause is `obj`, so that is what a fit lacking both is told.
  npfmc_nosd <- npfmc
  npfmc_nosd$sdrep <- NULL
  testthat::expect_error(avg(npfmc_nosd, hind), refused)
  # Point estimates draw nothing, so they are not refused.
  testthat::expect_no_error(model_average(list(npfmc, npfmc)))

  # Both remedies the message names reproduce the hindcast fit's SEs.
  rows <- which(names(hind$sdrep$value) == "log_ssb")
  rows <- rows[seq_len(d$nspp * (d$endyr - d$styr + 1))]
  remedies <- list(
    projection_uncertainty = fit(build_hcr(HCR = "NPFMC"), inits,
                                 projection_uncertainty = TRUE),
    Hindcast = fit(build_hcr(HCR = "NPFMC"), inits, estimateMode = 1))
  for (nm in names(remedies)) {
    set.seed(1)
    m <- remedies[[nm]]
    ratio <- avg(m, m)$sdrep$sd[rows] / hind$sdrep$sd[rows]
    testthat::expect_lt(abs(stats::median(ratio) - 1), 0.1, label = nm)
  }

  nosd <- fit(inits = inits, getsd = FALSE)
  testthat::expect_error(avg(hind, nosd), "no sdreport")
})
