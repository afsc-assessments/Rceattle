# The diet (stomach content) likelihood, reconstructed in R and checked against
# the TMB tape for both families it supports. It is the only test of that
# likelihood, and the C++ side is awkward enough to want one: the observed prey
# proportions get an "other prey" category appended, both observed and predicted
# are offset and renormalized, and only then is the density taken on counts.
#
# Integration test (fits a CEATTLE model): skipped on CRAN to keep R CMD check
# fast; the coverage workflow runs it in full (NOT_CRAN=true). See README.md.
testthat::skip_on_cran()

# Fit a two-species predation model at fixed parameters (estimateMode = 3, so
# nothing is estimated and the reconstruction has a known answer to hit).
.diet_fixture <- function(diet_distribution = NULL) {
  nyrs <- 30
  nspp <- 2
  Fmort  <- c(seq(0.02, 0.3, length.out = nyrs / 2),
              seq(0.3, 0.05, length.out = nyrs / 2))
  Fmort2 <- seq(0.02, 0.3, length.out = nyrs)
  gam_a   <- c(1, 0.1)
  gam_b   <- rep(0.15, nspp)
  log_phi <- matrix(c(-5, 0.5, -10, -2), nspp, nspp, byrow = TRUE)

  set.seed(123)
  sim <- make_msm_test_data(
    years   = 1:nyrs,
    Fmort   = matrix(c(Fmort, Fmort2), nspp, nyrs, byrow = TRUE),
    gam_a   = gam_a,
    gam_b   = gam_b,
    log_phi = log_phi)

  simData <- sim$data_list
  if (!is.null(diet_distribution)) {
    simData$Diet_distribution <- rep(diet_distribution, nspp)
  }

  ss_run <- Rceattle::fit_mod(
    data_list = simData, estimateMode = 3,
    fit_control = fit_control(phase = FALSE, verbose = 0))
  inits <- ss_run$estimated_params

  inits$log_gam_a           <- log(gam_a)
  inits$log_gam_b           <- log(gam_b)
  inits$log_phi             <- log_phi
  inits$sel_inf[1, , 1]     <- c(3, 6, 2.5, 4)
  inits$log_sel_slp[1, , 1] <- log(c(2, 2.5, 2, 2.5))
  inits$log_F[2, ]          <- log(Fmort)
  inits$log_F[4, ]          <- log(Fmort2)
  inits$rec_pars[, 1]       <- log(c(1e2, 1e3))
  inits$index_log_q[]       <- log(1)
  inits$R_log_sd[]          <- log(1)
  inits$rec_dev[, 1:30]     <- sim$model_quantities$rec_devs
  inits$init_dev[, 1:14]    <- sim$model_quantities$init_devs

  mod <- Rceattle::fit_mod(
    data_list = simData, inits = inits, file = NULL,
    estimateMode = 3, random_rec = FALSE, msmMode = 1, suitMode = 4,
    niter = 5, initMode = "NonEquilibrium",
    fit_control = fit_control(phase = FALSE, verbose = 0))
  list(mod = mod, nspp = nspp)
}


# Walk the stomachs as ceattle.cpp section 13.6 does, apply `dens` to each one's
# observed counts and predicted proportions, and accumulate per predator species.
.diet_nll_r <- function(mod, nspp, dens) {
  out <- rep(0, nspp)
  re_data <- Rceattle::rearrange_data(mod$data_list)
  stomach_ids <- re_data$stomach_id + 1   # stored 0-indexed for TMB

  for (i in unique(stomach_ids)) {
    idx <- which(stomach_ids == i)
    rsp <- re_data$diet_ctl[idx[1], 1]
    if (mod$data_list$suitMode[rsp] <= 0) next

    N_s       <- re_data$diet_obs[idx[1], 1]
    obs_prop  <- re_data$diet_obs[idx, 2]
    pred_prop <- mod$quantities$diet_hat[idx, 2]

    # "Other prey" is the balance of the stomach, appended as a final bin.
    obs_prop  <- c(obs_prop,  1.0 - sum(obs_prop))
    pred_prop <- c(pred_prop, 1.0 - sum(pred_prop))

    # The same offset and renormalization the C++ applies before the density.
    obs_prop  <- (obs_prop  + 1e-5) / sum(obs_prop  + 1e-5)
    pred_prop <- (pred_prop + 1e-5) / sum(pred_prop + 1e-5)

    out[rsp] <- out[rsp] + dens(obs_prop * N_s, pred_prop, N_s, rsp)
  }
  out
}


testthat::test_that("the multinomial diet likelihood matches the R math", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  f <- .diet_fixture()
  tmb <- f$mod$quantities$jnll_comp["Stomach content data", 1:f$nspp]

  # Multinomial: weighted by Diet_comp_weights for the predator.
  r <- .diet_nll_r(f$mod, f$nspp, function(counts, pred, N_s, rsp) {
    f$mod$data_list$Diet_comp_weights[rsp] * calc_multinom_nll(counts, pred)
  })

  testthat::expect_equal(as.numeric(tmb), as.numeric(r), tolerance = 1e-5)
  testthat::expect_true(all(is.finite(as.numeric(tmb))))
  testthat::expect_gt(sum(abs(as.numeric(tmb))), 0)   # the slot is actually used
})


testthat::test_that("the Dirichlet-multinomial diet likelihood matches the R math", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  f <- .diet_fixture(diet_distribution = "DirichletMultinomial")
  tmb <- f$mod$quantities$jnll_comp["Stomach content data", 1:f$nspp]

  # The weight enters through the concentration here
  # (alpha = pred * N * exp(diet_comp_weights)) rather than as a multiplier on
  # the density, so unlike the multinomial branch there is no outer weight.
  theta <- exp(as.numeric(f$mod$estimated_params$diet_comp_weights))
  ddirmult <- function(obs, alpha) {
    N <- sum(obs); phi <- sum(alpha)
    lgamma(N + 1) + lgamma(phi) - lgamma(N + phi) +
      sum(-lgamma(obs + 1) + lgamma(obs + alpha) - lgamma(alpha))
  }
  r <- .diet_nll_r(f$mod, f$nspp, function(counts, pred, N_s, rsp) {
    -ddirmult(counts, pred * N_s * theta[rsp])
  })

  testthat::expect_equal(as.numeric(tmb), as.numeric(r), tolerance = 1e-5)
  testthat::expect_true(all(is.finite(as.numeric(tmb))))
})
