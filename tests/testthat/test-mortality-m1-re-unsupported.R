# Every M1 random-effect family in build_map_m1() nests its deviation writes
# inside an M1_model test, and only M1_model 1 and 2 have arms. Under M1_model
# 3, 4 or 5 no arm fires, so log_M1_dev stays entirely NA -- while the
# standard deviation and the AR1 correlations were freed anyway, outside that
# guard.
#
# Measured on GOA2018SS (nspp 3, nsex c(1, 2, 1)) before the fix, free
# parameter counts from the built map, for every M1_re 1-6:
#
#   M1_model 1   log_M1_dev 4736 / 2688   M1_dev_log_sd 6   M1_rho 0 / 6 / 12
#   M1_model 2   log_M1_dev 6364 / 1764   M1_dev_log_sd 6   M1_rho 0 / 6 / 12
#   M1_model 3   log_M1_dev    0          M1_dev_log_sd 6   M1_rho 0 / 6 / 12
#   M1_model 4   log_M1_dev    0          M1_dev_log_sd 6   M1_rho 0 / 6 / 12
#   M1_model 5   log_M1_dev    0          M1_dev_log_sd 6   M1_rho 0 / 6 / 12
#
# So 18 of the 35 (M1_model, M1_re) pairs asked for time-varying M and got
# constant M, and M1_model 4 raised no condition at all. The free sd is worse
# than inert: it scores N(0, sigma) against a vector of zeros.
# `Rceattle-models/EBS pollock/2024/06-time-varying-M.R` measured that at
# 56.06 nats on a 61-year hindcast (61 * log(2*pi)/2), minimised by driving
# sigma to its bound, "which makes the objective incomparable with anything",
# and worked around it by telling the reader to avoid the combination.
#
# Now a requested random effect with nothing to attach it to is refused. The
# refusal reads whether an arm fired rather than restating the supported pairs,
# so adding an arm legalises its combination with no second registry to keep in
# step.
#
# Gating the sd write itself was considered and dropped: the refusal fires on
# exactly the complementary per-species condition, so the gate could never
# decide anything -- dead code by construction.

m1_map <- function(fit, d) {
  mp <- fit$obj$env$map$log_M1_dev
  list(dev = sum(!is.na(mp)),
       sd  = sum(!is.na(fit$obj$env$map$M1_dev_log_sd)),
       rho = sum(!is.na(fit$obj$env$map$M1_rho)))
}

fit_m1 <- function(mm, re) {
  data("GOA2018SS", package = "Rceattle", envir = environment())
  suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = GOA2018SS, inits = NULL, file = NULL,
    estimateMode = "DebugBuild", msmMode = 0, niter = 3, random_rec = FALSE,
    M1Fun = Rceattle::build_M1(M1_model = mm, M1_re = re),
    fit_control = Rceattle::fit_control(verbose = 0))))
}


testthat::test_that("a random effect with no deviations to scale is refused", {
  testthat::skip_on_cran()
  for (mm in 3:5) {
    for (re in 1:6) {
      testthat::expect_error(fit_m1(mm, re),
        "is not implemented for M1_model",
        info = paste("M1_model", mm, "M1_re", re))
    }
  }
})


testthat::test_that("the supported M1_re models still free their deviations", {
  testthat::skip_on_cran()
  data("GOA2018SS", package = "Rceattle", envir = environment())
  d <- suppressMessages(Rceattle::switch_check(GOA2018SS))

  # rho: 0 for an IID family, 6 for a single AR1, 12 for the separable 2D-AR1.
  want_rho <- c(`1` = 0L, `2` = 0L, `3` = 0L, `4` = 6L, `5` = 6L, `6` = 12L)
  for (mm in 1:2) {
    for (re in 1:6) {
      got <- m1_map(fit_m1(mm, re), d)
      lab <- paste("M1_model", mm, "M1_re", re)
      testthat::expect_gt(got$dev, 0L)
      # One sd per species, shared across sexes.
      testthat::expect_equal(got$sd, 6L, info = lab)
      testthat::expect_equal(got$rho, want_rho[[as.character(re)]], info = lab)
    }
  }
})


testthat::test_that("M1_re = 0 frees no deviation, sd or correlation", {
  # Passes either side of the fix -- with M1_re = 0 none of the family blocks
  # executes, so the sd was never freed there. Here to pin that the refusal
  # does not fire on the default, and that the no-random-effects case stays
  # exactly empty.
  testthat::skip_on_cran()
  data("GOA2018SS", package = "Rceattle", envir = environment())
  d <- suppressMessages(Rceattle::switch_check(GOA2018SS))
  for (mm in 1:5) {
    got <- m1_map(fit_m1(mm, 0), d)
    lab <- paste("M1_model", mm, "M1_re 0")
    testthat::expect_equal(got$dev, 0L, info = lab)
    # The sd is what used to be freed with nothing attached to it.
    testthat::expect_equal(got$sd, 0L, info = lab)
    testthat::expect_equal(got$rho, 0L, info = lab)
  }
})
