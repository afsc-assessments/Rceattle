# The q and selectivity deviations in years past a model's data (a retrospective
# peel's forecast, the years run_mse() adds) keep their last year's values, as
# the projection holds selectivity; a selectivity form whose deviations are
# year-to-year increments takes an increment of 0. M1 is left as it is.

testthat::test_that(".hold_devs() keeps the last year's q and selectivity deviations", {
  n_yrs <- 5
  last <- 3
  fc <- data.frame(
    Selectivity      = c("Logistic", "NonParametric", "NonParametricPM", "NonParametricIntegrable"),
    Time_varying_sel = c("RandomWalk", "IID", "RandomWalk", "RandomWalk"))
  ramp <- function(dims) array(seq_len(prod(dims)) / 10, dims)
  p <- list(index_q_dev     = ramp(c(4, n_yrs)),
            log_sel_slp_dev = ramp(c(2, 4, 1, n_yrs)),
            sel_inf_dev     = ramp(c(2, 4, 1, n_yrs)),
            sel_coff_dev    = ramp(c(4, 1, 3, n_yrs)),
            log_M1_dev      = ramp(c(1, 1, 2, n_yrs + 2)))   # runs on to projyr
  h <- Rceattle:::.hold_devs(p, list(fleet_control = fc), last)
  later <- (last + 1):n_yrs

  # Every deviation stored as a level keeps its last value, IID included; the
  # increment forms take an increment of 0.
  for (yr in later) {
    testthat::expect_identical(h$index_q_dev[, yr], p$index_q_dev[, last])
    testthat::expect_identical(h$log_sel_slp_dev[, , , yr], p$log_sel_slp_dev[, , , last])
    testthat::expect_identical(h$sel_inf_dev[, , , yr], p$sel_inf_dev[, , , last])
    testthat::expect_identical(h$sel_coff_dev[1:2, , , yr], p$sel_coff_dev[1:2, , , last])
    testthat::expect_true(all(h$sel_coff_dev[3:4, , , yr] == 0))
  }
  testthat::expect_identical(h$log_M1_dev, p$log_M1_dev)
  # The years with data are untouched.
  testthat::expect_identical(h$index_q_dev[, 1:last], p$index_q_dev[, 1:last])
  testthat::expect_identical(h$sel_coff_dev[, , , 1:last], p$sel_coff_dev[, , , 1:last])
})
