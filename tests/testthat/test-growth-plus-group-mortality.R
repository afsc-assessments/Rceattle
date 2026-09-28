# The plus group's mean length pools the ages it holds, each weighted by
# survival and its length interpolated from the oldest age towards L-infinity
# (growth.hpp, "Plus-Group Correction"). Until 5.35.0 the survival weight was a
# hard-coded exp(-0.2 a) whatever the species' natural mortality; it is now
# exp(-M1 a) at the species' base M1 for the oldest age.

testthat::skip_on_cran()

plus_group_build <- function(d) {
  suppressWarnings(suppressMessages(fit_mod(
    data_list = d, inits = NULL, estimateMode = 3, msmMode = 0, random_rec = FALSE,
    growthFun = build_growth(fun = "vonBertalanffy"),
    fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE))))
}

testthat::test_that("the plus-group length is the survival-weighted mean at the species' own M", {
  set.seed(3)
  d <- make_msm_test_data()$data_list
  m <- plus_group_build(d)
  q <- m$quantities
  na <- d$nages[1]; a_L1 <- m$data_list$growth_age_L1[1]
  gp <- q$growth_parameters[1, 1, 1, ]                   # K, L1, Linf, m in year 1
  cs <- gp[3] + (gp[2] - gp[3]) * exp(-gp[1] * ((na - 1 + d$minage[1]) - a_L1))
  M  <- exp(m$estimated_params$log_M1[1, 1, na])
  a  <- 0:na
  w  <- exp(-M * a)
  expected <- sum(w * (cs + (a / na) * (gp[3] - cs))) / sum(w)
  testthat::expect_equal(unname(q$length_hat[1, 1, na, 1]), unname(expected), tolerance = 1e-8)

  # A different M at the oldest age moves the plus-group length; the age below it is untouched.
  d2 <- d
  d2$M1_base[d2$M1_base$Species == 1, ncol(d2$M1_base)] <-
    2 * d2$M1_base[d2$M1_base$Species == 1, ncol(d2$M1_base)]
  m2 <- plus_group_build(d2)
  testthat::expect_false(isTRUE(all.equal(q$length_hat[1, 1, na, 1], m2$quantities$length_hat[1, 1, na, 1])))
  testthat::expect_equal(q$length_hat[1, 1, na - 1, 1], m2$quantities$length_hat[1, 1, na - 1, 1])
})
