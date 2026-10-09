# Everything downstream of the index likelihood has to know which SCALE the
# fleet is fitted on.
#
# `Lognormal` (0) is a log-scale family whose `Log_sd` is a CV / log-sd.
# `MVN` (1), `MVNORM` (2), `Normal` (3) and `TruncatedNormal` (4) are
# natural-scale families whose sd is ABSOLUTE, in the units of the index.
# Applying a log-scale formula to the second group does not error, it silently
# returns nonsense: `sigma^2 / 2` becomes a number the size of the index squared.
#
# These pin the two places that got it wrong -- the Pearson residual and the
# observation interval -- with fixtures where the wrong answer is not subtle.
#
# Guarded per test rather than per file: the registry-drift check below runs no
# fit and is the one test here that should run under a plain R CMD check.

.ns_fixture <- function(dist, sd, nyrs = 12) {
  dat <- make_test_data(nyrs = nyrs, nages = 5, seed = 123)
  srv <- dat$fleet_control$Fleet_name == "Survey"
  dat$fleet_control$Index_distribution[srv] <- dist
  if (dist %in% c("MVN", "MVNORM")) {
    dat$fleet_control$Catchability[srv] <- "AnalyticalArith"
    dat$index_cov <- list(Survey = diag(rep(sd^2, nyrs)))
  }
  dat$index_data$Log_sd[dat$index_data$Fleet_name == "Survey"] <- sd
  suppressMessages(suppressWarnings(Rceattle::fit_mod(
    dat, file = NULL, estimateMode = 1, msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0, phase = FALSE))))
}

testthat::test_that("Pearson index residuals use the fleet's own scale", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  fit <- .ns_fixture("Normal", sd = 150)
  r   <- stats::residuals(fit, type = "pearson", source = "index")
  srv <- r$Fleet_name == "Survey" & is.finite(r$Residual)
  res <- r$Residual[srv]

  # Natural scale: (obs - hat) / sd. With sd = 150 these are O(1).
  idx <- fit$data_list$index_data
  sel <- idx$Fleet_name == "Survey"
  want <- (idx$Observation[sel] - as.numeric(fit$quantities$index_hat)[sel]) / 150
  testthat::expect_equal(res, want[is.finite(want)][seq_along(res)],
                         tolerance = 1e-6)

  # The lognormal formula would return the same large constant for every row --
  # about +75 here, since -sd^2/2 dominates. Assert we are nowhere near it.
  testthat::expect_lt(max(abs(res)), 10)
  testthat::expect_gt(stats::sd(res), 0)

  # A lognormal fleet is untouched by the fix.
  ln  <- .ns_fixture("Lognormal", sd = 0.2)
  rl  <- stats::residuals(ln, type = "pearson", source = "index")
  sl  <- rl$Fleet_name == "Survey" & is.finite(rl$Residual)
  idl <- ln$data_list$index_data
  sel <- idl$Fleet_name == "Survey"
  sig <- as.numeric(ln$quantities$index_sd)[sel]
  wl  <- (log(idl$Observation[sel]) -
            (log(as.numeric(ln$quantities$index_hat)[sel]) - sig^2 / 2)) / sig
  testthat::expect_equal(rl$Residual[sl], wl[is.finite(wl)][seq_len(sum(sl))],
                         tolerance = 1e-6)
})

testthat::test_that("index observation intervals use the fleet's own scale", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("ggplot2")

  for (dist in c("Normal", "TruncatedNormal", "MVN")) {
    fit <- .ns_fixture(dist, sd = 150)
    df  <- Rceattle:::.fleet_fit_df(list(fit), kind = "index")
    d   <- df[df$Fleet == "Survey" & is.finite(df$Upper95), ]
    testthat::expect_gt(nrow(d), 0)

    # qlnorm(0.975, log(obs), 150) overflows to ~1e130; the natural-scale band is
    # obs +/- 1.96*150, so a few hundred at most on this fixture.
    testthat::expect_lt(max(d$Upper95), 1e5, label = dist)
    # An index cannot be negative, so the lower bound is clamped at zero.
    testthat::expect_gte(min(d$Lower95), 0)
  }

  # Lognormal keeps the multiplicative band.
  ln <- .ns_fixture("Lognormal", sd = 0.2)
  dl <- Rceattle:::.fleet_fit_df(list(ln), kind = "index")
  dl <- dl[dl$Fleet == "Survey" & is.finite(dl$Upper95), ]
  testthat::expect_gt(nrow(dl), 0)
  testthat::expect_equal(dl$Upper95,
                         stats::qlnorm(0.975, log(dl$Observation), 0.2),
                         tolerance = 1e-6)
})


testthat::test_that("every natural-scale family is recognized as one", {
  testthat::skip_if_not_installed("TMB")

  # The scale question is answered once, in .index_rows_natural_scale(). A family
  # added to index_distribution_map but not to that vector reverts silently to
  # the log-scale treatment -- the exact failure the rest of this file pins. This
  # enumerates the map rather than a hand-written list, so a new family fails
  # here rather than in a user's residual plot.
  d <- make_test_data(nyrs = 12, nages = 5, seed = 123)
  srv <- d$fleet_control$Fleet_name == "Survey"

  for (nm in names(index_distribution_map)) {
    dat <- d
    dat$fleet_control$Index_distribution[srv] <- nm
    nat <- Rceattle:::.index_rows_natural_scale(dat)
    rows <- dat$index_data$Fleet_name == "Survey"
    want <- nm != "Lognormal"
    testthat::expect_equal(unique(nat[rows]), want, label = nm)
    # ...and the code path agrees with the name path, since fleet_control may
    # hold either spelling.
    dat$fleet_control$Index_distribution[srv] <- index_distribution_map[[nm]]
    testthat::expect_equal(unique(Rceattle:::.index_rows_natural_scale(dat)[rows]),
                           want, label = paste(nm, "as code"))
  }
})

testthat::test_that("TruncatedNormal gets natural-scale Pearson residuals", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  fit <- .ns_fixture("TruncatedNormal", sd = 150)
  r   <- stats::residuals(fit, type = "pearson", source = "index")
  srv <- r$Fleet_name == "Survey" & is.finite(r$Residual)
  res <- r$Residual[srv]

  idx <- fit$data_list$index_data
  sel <- idx$Fleet_name == "Survey"
  want <- (idx$Observation[sel] - as.numeric(fit$quantities$index_hat)[sel]) / 150
  testthat::expect_equal(res, want[is.finite(want)][seq_along(res)],
                         tolerance = 1e-6)
  # Nowhere near the ~+75 the lognormal formula would return.
  testthat::expect_lt(max(abs(res)), 10)
})

testthat::test_that("the analytical sd is refused where the likelihood reads it", {
  testthat::skip_if_not_installed("TMB")

  # index_analytical_sd accumulates squared LOG residuals (Walters and Ludwig
  # 1994), so it is a log-scale sd. What that costs splits by family, and so does
  # the response:
  #   Normal / TruncatedNormal read the sd, so the LIKELIHOOD is on the wrong
  #     scale -- an error.
  #   MVN / MVNORM score through index_cov and never read the scalar sd, so the
  #     fit is fine and only the reported index_sd (and the Pearson residual and
  #     plot interval built from it) is wrong -- a warning, because refusing
  #     would reject a model that fits correctly.
  mk <- function(fam, est) {
    dat <- make_test_data(nyrs = 12, nages = 5, seed = 123)
    srv <- dat$fleet_control$Fleet_name == "Survey"
    dat$fleet_control$Index_distribution[srv] <- fam
    dat$fleet_control$Estimate_index_sd[srv]  <- est
    if (fam %in% c("MVN", "MVNORM")) {
      dat$fleet_control$Catchability[srv] <- "AnalyticalArith"
      dat$index_cov <- list(Survey = diag(rep(20^2, 12)))
    }
    dat
  }
  fit <- function(dat) suppressMessages(Rceattle::fit_mod(
    dat, file = NULL, estimateMode = 3, msmMode = 0,
    fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0)))

  for (fam in c("Normal", "TruncatedNormal")) {
    testthat::expect_error(suppressWarnings(fit(mk(fam, "Analytical"))),
                           "wrong scale", label = fam)
  }

  # The covariance families warn and still fit.
  for (fam in c("MVN", "MVNORM")) {
    testthat::expect_warning(m <- fit(mk(fam, "Analytical")),
                             "never read the scalar sd", label = fam)
    testthat::expect_false(is.null(m$obj), label = fam)
  }

  # Lognormal is the family the analytical route was derived for, and an "Off"
  # fleet is not checked at all -- or this would just be banning working setups.
  testthat::expect_no_error(suppressWarnings(fit(mk("Lognormal", "Analytical"))))
  d_off <- mk("Normal", "Analytical")
  d_off$fleet_control$Fleet_type[d_off$fleet_control$Fleet_name == "Survey"] <- "Off"
  testthat::expect_no_error(suppressWarnings(fit(d_off)))
})


testthat::test_that("plot_indexresidual and residuals() agree on sign and scale", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # Two user-facing views of one quantity. Until 5.9.0 the plotter drew
  # predicted - observed while residuals() returned observed - predicted, so the
  # two were mirror images and a reader could take the sign of a survey residual
  # backwards. Both are observed - predicted now; this pins that, on both scales.
  for (fam in c("Lognormal", "Normal")) {
    fit <- .ns_fixture(fam, sd = if (fam == "Lognormal") 0.2 else 20)

    plotted <- Rceattle::plot_indexresidual(fit)$data
    direct  <- stats::residuals(fit, type = "response", source = "index",
                                scale = if (fam == "Lognormal") "log" else "natural")
    testthat::skip_if(!nrow(plotted), "no index residuals")

    # Same rows, same order: both drop Year <= 0 and keep the fitted window.
    testthat::expect_equal(nrow(plotted), nrow(direct), label = fam)
    testthat::expect_equal(as.numeric(plotted$Residual),
                           as.numeric(direct$Residual),
                           tolerance = 1e-10, label = fam)
    # ...and it is the observed-minus-predicted orientation, not its negative.
    testthat::expect_false(
      isTRUE(all.equal(as.numeric(plotted$Residual),
                       -as.numeric(direct$Residual))),
      label = fam)
  }
})
