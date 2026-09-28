# OSA (one-step-ahead) residuals across ALL survey-index likelihood families.
#
# The index likelihood family is set by fleet_control$Index_distribution ->
# index_ll_type: 0 lognormal IID, 1 MVN (bare quadratic form), 2 MVNORM (full
# density), 3 natural-scale Normal. osa_residuals() must produce statistically
# correct residuals for every family:
#   * lognormal / Normal -> independent (log- or natural-scale) normal residuals;
#   * MVN / MVNORM -> the correlated survey block is whitened by the lower Cholesky
#     of its covariance Sigma = L L', so the residuals are the multivariate-Gaussian
#     one-step-ahead innovations u = L^-1 (obs - q*pred). For a Gaussian block this
#     is exactly the closed form that TMB::oneStepPredict() (the SAM/TMB OSA engine)
#     reproduces -- Thygesen, Albertsen, Berg, Kristensen & Nielsen (2017),
#     "Validation of ecological state space models using the Laplace approximation".
# With random_rec = FALSE the estimated model has no random effects, so the
# oneStepPredict residual equals this closed form to numerical precision -- which is
# the check below (osa_residuals() runs oneStepPredict internally).
#
# External cross-check (when installed): the SAM-author / WHAM OSA packages
#   compResidual::resmvnorm   -- devtools::install_github("fishfollower/compResidual/compResidual")
#   OSA_multivariate_dists    -- TMB:::install.contrib("https://github.com/vtrijoulet/OSA_multivariate_dists/archive/main.zip")
# implement the same multivariate-normal one-step-ahead residual; the last test
# below verifies Rceattle's MVN OSA residuals against compResidual::resmvnorm().

.osa_index_fit <- function(dist, nyrs = 8, nages = 5, seed = 42, rho = 0.3,
                           sd = 20, absolute_sd = FALSE) {
  dat <- make_test_data(nyrs = nyrs, nages = nages, seed = seed)
  srv <- dat$fleet_control$Fleet_name == "Survey"
  dat$fleet_control$Index_distribution[srv] <- dist
  dat$fleet_control$Catchability[srv]       <- "AnalyticalArith"
  Sigma <- NULL
  if (dist %in% c("MVN", "MVNORM")) {
    Rho   <- matrix(rho, nyrs, nyrs); diag(Rho) <- 1
    Sigma <- diag(rep(sd, nyrs)) %*% Rho %*% diag(rep(sd, nyrs))
    dat$index_cov <- list(Survey = Sigma)
  }
  if (absolute_sd) dat$index_data$Log_sd[dat$index_data$Fleet_name == "Survey"] <- sd
  fit <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    dat, file = NULL, estimateMode = 1, msmMode = 0, random_rec = FALSE,
    fit_control = fit_control(getsd = FALSE, verbose = 0, phase = FALSE))))
  list(fit = fit, Sigma = Sigma)
}

# Survey (Fleet_code 1) fitted index rows, in index_obs row order (== Sigma order).
.fitted_index <- function(fit) {
  td <- fit$obj$env$data
  which(td$index_ctl[, 1] == 1 & td$index_ctl[, 3] > 0 &
          td$index_ctl[, 3] <= td$endyr & td$index_obs[, 1] > 0)
}

testthat::test_that("MVN / MVNORM OSA residuals equal the Cholesky innovation L^-1 (obs - pred)", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  for (dist in c("MVN", "MVNORM")) {
    o <- .osa_index_fit(dist); fit <- o$fit; Sigma <- o$Sigma
    osa <- suppressWarnings(osa_residuals(fit, source = "index"))
    testthat::expect_s3_class(osa, "rceattle_osa")
    testthat::expect_true(all(is.finite(osa$residual)), info = dist)

    # Independent oracle: the closed-form multivariate-Gaussian one-step-ahead
    # residual = L^-1 (obs - index_hat), with L the lower Cholesky of Sigma. This
    # is what oneStepPredict() reproduces for a Gaussian block; osa_residuals()
    # runs oneStepPredict internally, so an exact match validates the whole path.
    sel <- .fitted_index(fit)
    rp  <- fit$obj$report(fit$obj$env$last.par.best)
    innov <- as.numeric(forwardsolve(t(chol(Sigma)),
                                     fit$obj$env$data$index_obs[sel, 1] - rp$index_hat[sel]))
    # osa rows are ordered by year (ascending) == sel order.
    testthat::expect_equal(osa$residual[order(osa$year)], innov, tolerance = 1e-6,
                           info = dist)
  }
})

testthat::test_that("MVN whitening is consistent: OSA jnll = fit quadratic form up to 0.5*logdet(Sigma)", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # The whitened OSA density (independent standard normals on z_k = (L^-1 obs)_k,
  # keep == 1) equals the fitting MVNORM density minus 0.5*logdet(Sigma) -- the
  # normalizing constant the whitening drops. An exact offset proves the R-side L
  # (chol) and the C++-side L (Eigen LLT) are the identical factor, i.e. the
  # observations and the predicted mean are whitened consistently.
  o <- .osa_index_fit("MVNORM"); fit <- o$fit; Sigma <- o$Sigma
  jc0  <- fit$obj$report(fit$obj$env$last.par.best)$jnll_comp[1, 1]   # MVNORM fit, fleet 1
  obj1 <- Rceattle:::.osa_build_obj(fit)                              # osa_mode = 1
  jc1  <- obj1$report(obj1$par)$jnll_comp[1, 1]                       # whitened, keep == 1
  logdet <- as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)
  testthat::expect_equal(jc1, jc0 - 0.5 * logdet, tolerance = 1e-6)
})

testthat::test_that("natural-scale Normal OSA residuals equal (obs - pred) / sd", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  sd <- 20
  o  <- .osa_index_fit("Normal", sd = sd, absolute_sd = TRUE); fit <- o$fit
  osa <- suppressWarnings(osa_residuals(fit, source = "index"))
  testthat::expect_true(all(is.finite(osa$residual)))
  sel <- .fitted_index(fit)
  rp  <- fit$obj$report(fit$obj$env$last.par.best)
  ana <- (fit$obj$env$data$index_obs[sel, 1] - rp$index_hat[sel]) / sd
  testthat::expect_equal(osa$residual[order(osa$year)], ana, tolerance = 1e-6)
})

testthat::test_that("the survey covariance actually shapes the residuals (whitening != naive standardize)", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # With a diagonal Sigma (rho = 0) whitening reduces to element-wise
  # standardizing, so the OSA residuals equal (obs - pred)/sd. With an
  # off-diagonal Sigma (rho = 0.3) the whitening mixes years and the residuals
  # depart from the naive per-year standardization -- proving the covariance is
  # genuinely used, not ignored.
  sd <- 20
  naive_std <- function(fit) {
    sel <- .fitted_index(fit); rp <- fit$obj$report(fit$obj$env$last.par.best)
    (fit$obj$env$data$index_obs[sel, 1] - rp$index_hat[sel]) / sd
  }

  o0 <- .osa_index_fit("MVNORM", rho = 0, sd = sd); fit0 <- o0$fit
  osa0 <- suppressWarnings(osa_residuals(fit0, source = "index"))
  testthat::expect_equal(osa0$residual[order(osa0$year)], naive_std(fit0), tolerance = 1e-6)

  oR <- .osa_index_fit("MVNORM", rho = 0.3, sd = sd); fitR <- oR$fit
  osaR <- suppressWarnings(osa_residuals(fitR, source = "index"))
  testthat::expect_false(isTRUE(all.equal(osaR$residual[order(osaR$year)], naive_std(fitR),
                                          tolerance = 1e-6)))
})

testthat::test_that("MVN OSA residuals run with random recruitment (Laplace re-conditioning)", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # With random_rec = TRUE the survey mean depends on random recruitment
  # deviations, so oneStepPredict() re-conditions the random effects as each
  # observation is added -- the whitening (a constant-L linear map of the mean)
  # stays valid but the residuals are no longer the fixed-parameter closed form.
  # There is no simple analytic oracle for the marginalized case, so pin only that
  # the machinery runs and returns finite, reasonably-calibrated residuals.
  nyrs <- 12; nages <- 5
  dat  <- make_test_data(nyrs = nyrs, nages = nages, seed = 42)
  srv  <- dat$fleet_control$Fleet_name == "Survey"
  dat$fleet_control$Index_distribution[srv] <- "MVNORM"
  dat$fleet_control$Catchability[srv]       <- "AnalyticalArith"
  Rho <- matrix(0.3, nyrs, nyrs); diag(Rho) <- 1
  dat$index_cov <- list(Survey = diag(rep(20, nyrs)) %*% Rho %*% diag(rep(20, nyrs)))

  fit <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    dat, file = NULL, estimateMode = 1, msmMode = 0, random_rec = TRUE,
    fit_control = fit_control(getsd = FALSE, verbose = 0, phase = FALSE))))
  osa <- suppressWarnings(osa_residuals(fit, source = "index"))

  # No analytic oracle for the marginalized (random-effect) case, and the fixture
  # is tiny synthetic data, so assert only that the machinery runs and returns
  # finite, non-exploded residuals (not the fixed-parameter closed form).
  testthat::expect_true(all(is.finite(osa$residual)))
  testthat::expect_equal(nrow(osa), nyrs)
  testthat::expect_true(all(abs(osa$residual) < 10))     # not blown up
})

testthat::test_that("MVN OSA residuals match compResidual::resmvnorm (SAM-author oracle)", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("compResidual")   # GitHub-only; local cross-check

  # Cross-check Rceattle's whitened MVN OSA residuals against the independent
  # multivariate-normal one-step-ahead residual from compResidual (Nielsen, the SAM
  # author). On the same fitted state both must return L^-1 (obs - index_hat).
  o <- .osa_index_fit("MVNORM"); fit <- o$fit; Sigma <- o$Sigma
  osa <- suppressWarnings(osa_residuals(fit, source = "index"))
  sel <- .fitted_index(fit)
  rp  <- fit$obj$report(fit$obj$env$last.par.best)
  obs <- fit$obj$env$data$index_obs[sel, 1]
  mu  <- rp$index_hat[sel]

  r_sam <- as.numeric(compResidual::resmvnorm(matrix(obs, ncol = 1),
                                              matrix(mu, ncol = 1), Sigma))
  testthat::expect_equal(osa$residual[order(osa$year)], r_sam, tolerance = 1e-6)
})
testthat::test_that("method = \"cdf\" reproduces the MVN whitened innovation", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # The whitened block is independent standard normals, so its conditional CDF is
  # pnorm of the innovation and method = "cdf" must return the same closed form
  # the Gaussian methods do. This is the only cover the MVN CDF call site has.
  for (dist in c("MVN", "MVNORM")) {
    o <- .osa_index_fit(dist); fit <- o$fit; Sigma <- o$Sigma
    osa <- suppressWarnings(suppressMessages(
      osa_residuals(fit, source = "index", method = "cdf", parallel = FALSE)))
    sel <- .fitted_index(fit)
    rp  <- fit$obj$report(fit$obj$env$last.par.best)
    innov <- as.numeric(forwardsolve(t(chol(Sigma)),
                                     fit$obj$env$data$index_obs[sel, 1] - rp$index_hat[sel]))
    testthat::expect_equal(osa$residual[order(osa$year)], innov, tolerance = 1e-5,
                           info = dist)
    testthat::expect_true(all(is.na(osa$predicted)), info = dist)
  }
})


testthat::test_that("a survey whose rows are out of year order is dropped, with a warning", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # The lower-triangular whitening conditions row k on rows 1..k-1, so it is the
  # one-step-ahead conditioning order only while the rows are chronological and
  # aligned with Sigma. Out of order, the residuals would be a valid
  # decomposition of the wrong conditioning sequence -- so the fleet is dropped
  # instead, and the warning has to say that the drop is from the residual
  # MODEL, not just the output: the other fleets' residuals move with it under
  # random effects.
  o  <- .osa_index_fit("MVNORM")
  dl <- o$fit$obj$env$data
  rows <- .fitted_index(o$fit)
  testthat::expect_gt(length(rows), 1)

  # Same rows, same Sigma, years no longer ascending.
  dl$index_ctl[rows, 3] <- rev(dl$index_ctl[rows, 3])
  testthat::expect_warning(
    osa_dat <- build_osa_data(dl, build_osa = TRUE),
    "out of year order")

  # Dropped means excluded from obsvec, which is what makes the residual model
  # differ from the fitted one -- not merely filtered out of the returned frame.
  testthat::expect_true(all(osa_dat$index_obsvec_idx[rows] == -1L))
})
