# The descending limb of the double logistic and descending logistic forms at a
# steep slope: selectivity at the oldest ages must stay strictly positive and
# match 1 / (1 + exp(slope * (age - a50))) to 1e-10 relative error. Written
# as 1 - 1 / (1 + exp(-slope * (age - a50))), the limb is exactly 0 in double
# precision once slope * (age - a50) passes about 37 (ages 13 to 15 below).
testthat::skip_on_cran()

testthat::test_that("a steep descending limb stays positive and matches the closed form", {
  testthat::skip_if_not_installed("TMB")

  data_list <- make_test_data(nages = 15, seed = 1)
  data_list$fleet_control$Selectivity <- c("DescendingLogistic", "DoubleLogistic")

  slp1 <- 2; inf1 <- 3; slp2 <- 8; inf2 <- 8
  mod0 <- suppressMessages(Rceattle::fit_mod(
    data_list = data_list, inits = NULL, estimateMode = 3, random_rec = FALSE,
    msmMode = 0, fit_control = fit_control(verbose = 0)))
  inits <- mod0$estimated_params
  inits$log_sel_slp[1, , ] <- log(slp1)
  inits$sel_inf[1, , ] <- inf1
  inits$log_sel_slp[2, , ] <- log(slp2)
  inits$sel_inf[2, , ] <- inf2

  fit <- suppressMessages(Rceattle::fit_mod(
    data_list = data_list, inits = inits, estimateMode = 3, random_rec = FALSE,
    msmMode = 0, fit_control = fit_control(verbose = 0)))

  # The curve is evaluated at bin ordinals 1..nages, which are the ages here
  # because make_test_data() has minage = 1. Sel_norm_bin is NA: no normalization.
  ages <- 1:15
  descending <- 1 / (1 + exp(slp2 * (ages - inf2)))
  ascending <- 1 / (1 + exp(-slp1 * (ages - inf1)))
  expected <- list(descending, ascending * descending)

  # Fleet 1 is the descending logistic, fleet 2 the double logistic.
  for (flt in 1:2) {
    sel <- as.numeric(fit$quantities$sel_at_age[flt, 1, ages, 1])
    testthat::expect_true(all(sel > 0), info = paste("fleet", flt))
    testthat::expect_lt(max(abs(sel / expected[[flt]] - 1)), 1e-10,
                        label = paste("max relative error, fleet", flt))
  }
})
