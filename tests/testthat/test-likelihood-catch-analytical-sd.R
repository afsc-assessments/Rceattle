# Estimate_catch_sd = "Analytical" (2) concentrates the observation sd out of
# the catch likelihood instead of fixing or estimating it.
#
# Until 5.12.0 the option was documented in the schema, accepted by
# validate_switches(), and treated as a real mode by build_map() (which maps
# catch_log_sd out) -- but the template's dispatch had only case 0 and case 1,
# so the fit died with "Invalid 'Estimate_sigma_catch'" after passing every R
# check. est_sigma_index had implemented case 2 all along; this is its mirror.
#
# Rceattle fits catch to a mean-unbiased prediction
# (log(obs) ~ N(log(pred) - b*sigma^2/2, sigma), b = bias_adjust_obs, default 1),
# so the sd that minimises that density is
#   sigma^2 = 2*S / (sqrt(1 + b^2*S) + 1),   S = mean squared log residual,
# which reduces to the Walters and Ludwig (1994) sqrt(S) only at b = 0. The
# second test below is the one that pins this: it profiles the real catch
# likelihood over sigma rather than re-deriving the same expression the
# template computes, which is what an earlier draft did and could not have
# caught the bias term being left out.

test_that("the analytical catch sd equals the closed form", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  data("BS2017SS")
  d <- BS2017SS
  # NOTE: Fleet_type is integer-coded on the raw bundled data (1 = Fishery).
  # Subsetting on the string "Fishery" here matches nothing and silently
  # produces a mode-0 model -- which is how a first draft of this test passed
  # while exercising nothing.
  d$fleet_control$Estimate_catch_sd[d$fleet_control$Fleet_type == 1] <- 2

  fit <- suppressMessages(suppressWarnings(
    Rceattle::fit_mod(data_list = d, estimateMode = 3, msmMode = 0,
                      fit_control = Rceattle::fit_control(verbose = 0))))

  # The switch really reached the template.
  testthat::expect_true(all(fit$obj$env$data$est_sigma_fsh[1:3] == 2))

  dl <- fit$obj$env$data
  cc <- dl$catch_ctl; co <- dl$catch_obs
  hat <- fit$quantities$catch_hat
  got <- fit$quantities$catch_analytical_sd
  b   <- dl$bias_adjust_obs

  checked <- 0L
  for (flt in sort(unique(cc[, 1]))) {
    if (dl$flt_type[flt] != 1) next
    rows <- which(cc[, 1] == flt & cc[, 3] > 0 & cc[, 3] <= dl$endyr & co[, 1] > 0)
    if (!length(rows)) next
    checked <- checked + 1L
    S <- sum((log(co[rows, 1]) - log(hat[rows]))^2) / length(rows)
    want <- sqrt(2 * S / (sqrt(1 + b^2 * S) + 1))
    testthat::expect_equal(as.numeric(got[flt]), want, tolerance = 1e-10,
                           info = paste("fleet", flt))
  }
  testthat::expect_gt(checked, 0L)
})

test_that("the analytical catch sd minimises the catch likelihood", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # The estimator is only worth the name if it is the argmin of the density the
  # model actually evaluates. Profiling that density directly is independent of
  # how the template computes it: the plain Walters-Ludwig sqrt(S) fails this
  # by 4-17% on these three fisheries, because it ignores the -b*sigma^2/2 the
  # mean of the density carries.
  data("BS2017SS")
  d <- BS2017SS
  d$fleet_control$Estimate_catch_sd[d$fleet_control$Fleet_type == 1] <- 2

  fit <- suppressMessages(suppressWarnings(
    Rceattle::fit_mod(data_list = d, estimateMode = 1, msmMode = 0,
                      fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                                          verbose = 0))))
  dl <- fit$obj$env$data
  cc <- dl$catch_ctl; co <- dl$catch_obs
  hat <- fit$quantities$catch_hat
  b   <- dl$bias_adjust_obs
  got <- fit$quantities$catch_analytical_sd

  checked <- 0L
  for (flt in sort(unique(cc[, 1]))) {
    if (dl$flt_type[flt] != 1) next
    rows <- which(cc[, 1] == flt & cc[, 3] > 0 & cc[, 3] <= dl$endyr & co[, 1] > 0)
    if (!length(rows)) next
    checked <- checked + 1L
    nll <- function(s) {
      -sum(stats::dnorm(log(co[rows, 1]), log(hat[rows]) - b * s^2 / 2, s, log = TRUE))
    }
    mle <- stats::optimize(nll, c(1e-4, 5), tol = 1e-12)$minimum
    testthat::expect_equal(as.numeric(got[flt]), mle, tolerance = 1e-6,
                           info = paste("fleet", flt))
  }
  testthat::expect_gt(checked, 0L)
})

test_that("the analytical sd is the one the likelihood actually applies", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  data("BS2017SS")
  d <- BS2017SS
  d$fleet_control$Estimate_catch_sd[d$fleet_control$Fleet_type == 1] <- 2

  fit <- suppressMessages(suppressWarnings(
    Rceattle::fit_mod(data_list = d, estimateMode = 1, msmMode = 0,
                      fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                                          verbose = 0))))
  # Computing it and then not using it would look identical from the outside.
  applied <- fit$quantities$catch_sd
  applied <- sort(unique(round(applied[applied > 0], 8)))
  analytical <- sort(round(fit$quantities$catch_analytical_sd[1:3], 8))
  testthat::expect_equal(applied, as.numeric(analytical))

  testthat::expect_true(is.finite(fit$opt$objective))
})

test_that("a fishery with no fitted catch observation is refused, not given sd 0", {
  testthat::skip_if_not_installed("TMB")

  # A fishery whose hindcast catch is all zero has no residuals to concentrate,
  # so the analytical sd is undefined. It used to fall through as 0, which the
  # likelihood never reads but catch_sd reports and the SIMULATE draw uses --
  # rnorm(mean, 0) is a deterministic catch that sim_mod() would write back as
  # data. A zero-catch year is legal input, which is what makes this reachable.
  data("BS2017SS")
  d <- BS2017SS
  d$fleet_control$Estimate_catch_sd[d$fleet_control$Fleet_type == 1] <- 2
  # Zero the hindcast catch for fleet 1 only, leaving its rows in place.
  zeroed <- d$catch_data$Fleet_code == 1 & d$catch_data$Year <= d$endyr
  d$catch_data$Catch[zeroed] <- 0

  testthat::expect_error(
    suppressMessages(suppressWarnings(
      Rceattle::fit_mod(data_list = d, estimateMode = 3, msmMode = 0,
                        fit_control = Rceattle::fit_control(verbose = 0)))),
    "analytical")
})

test_that("a model that does not ask for it is bit-identical", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # The whole point: a new dispatch arm must not perturb the modes that were
  # already there. BS2017SS ships at Estimate_catch_sd = 0.
  data("BS2017SS")
  fit <- suppressMessages(suppressWarnings(
    Rceattle::fit_mod(data_list = BS2017SS, estimateMode = 1, msmMode = 0,
                      fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                                          verbose = 0))))
  # The pinned value for this configuration; see test-golden-regression.R for
  # the phased/polished references.
  testthat::expect_equal(fit$opt$objective, 10241.03042750, tolerance = 1e-6)
})
