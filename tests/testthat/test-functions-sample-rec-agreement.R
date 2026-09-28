# sample_rec() projects recruitment two ways: resampling hindcast deviations
# (sample_rec = TRUE) or one deviation for every projection year (FALSE). The
# two must agree in expectation, on the natural scale, at both draw sites: the
# curve fitted in the hindcast (pool = rec_dev) and the Ianelli penalty (pool =
# log R - log R_hat over the penalty years). The FALSE value is the log of the
# mean of exp() over the pool, so exp(FALSE) equals the mean of the resampled
# multipliers exactly, which is what these tests read off the pool rather than
# from a finite number of draws. Both sites are exercised at estimateMode = 3,
# where the deviations are the build's starting values.

srr_agree_data <- function() {
  set.seed(123)
  make_msm_test_data()$data_list
}

srr_agree_build <- function(recFun) {
  suppressMessages(suppressWarnings(fit_mod(
    data_list = srr_agree_data(), estimateMode = 3, msmMode = 1, suitMode = 0,
    initMode = "NonEquilibrium", random_rec = FALSE, recFun = recFun,
    fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE))))
}

proj_devs <- function(m, resample) {
  out <- suppressMessages(suppressWarnings(
    sample_rec(m, sample_rec = resample, update_model = FALSE)))
  hind <- m$data_list$endyr - m$data_list$styr + 1
  out$estimated_params$rec_dev[, -seq_len(hind), drop = FALSE]
}

testthat::test_that("a curve fitted in the hindcast projects the mean multiplier", {
  testthat::skip_on_cran()
  m <- srr_agree_build(build_srr(srr_fun = "BevertonHolt",
                                 srr_alpha_init = c(1, 10), srr_beta_init = c(0.01, 0.01)))
  # Give the pool spread so the mean and the median differ.
  hind <- m$data_list$endyr - m$data_list$styr + 1
  m$estimated_params$rec_dev[, seq_len(hind)] <- matrix(seq(-1, 1, length.out = hind),
                                                        nrow = 2, ncol = hind, byrow = TRUE)
  fixed <- proj_devs(m, FALSE)
  pool  <- m$estimated_params$rec_dev[, seq_len(hind), drop = FALSE]
  for (sp in 1:2) {
    testthat::expect_equal(unique(fixed[sp, ]), log(mean(exp(pool[sp, ]))))
    # Resampling draws from that pool and nothing else.
    set.seed(1)
    drawn <- proj_devs(m, TRUE)
    testthat::expect_true(all(drawn[sp, ] %in% pool[sp, ]))
  }
})

testthat::test_that("the Ianelli penalty projects the mean of R / R_hat over the penalty years", {
  testthat::skip_on_cran()
  m <- srr_agree_build(build_srr(srr_fun = "mean", srr_pred_fun = "BevertonHolt",
                                 srr_alpha_init = c(1, 10), srr_beta_init = c(0.01, 0.01)))
  hind  <- m$data_list$endyr - m$data_list$styr + 1
  fixed <- proj_devs(m, FALSE)
  hat   <- Rceattle:::.srr_hat_cols(m$data_list, hind)
  ratio <- (m$quantities$R / m$quantities$R_hat)[, hat, drop = FALSE]
  for (sp in 1:2) {
    testthat::expect_equal(unique(fixed[sp, ]), log(mean(ratio[sp, ])))
  }
  # Resampling draws log(R / R_hat) from every hindcast year, the documented
  # difference from the penalty-year mean (inst/dev/TODO-srr-multispecies.md).
  set.seed(1)
  drawn <- proj_devs(m, TRUE)
  pool  <- (log(m$quantities$R) - log(m$quantities$R_hat))[, seq_len(hind), drop = FALSE]
  for (sp in 1:2) testthat::expect_true(all(drawn[sp, ] %in% pool[sp, ]))
})
