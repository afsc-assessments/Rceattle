# Catchability, selectivity and M1 in years no data inform (a retrospective
# peel's forecast, the years run_mse() adds) keep their last year's values: each
# deviation keeps its last value, except in the selectivity forms that store
# each year's step, whose step is 0. Projection years are left as they are.

testthat::test_that(".hold_curves() keeps the last year's q, selectivity and M1", {
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
  h <- Rceattle:::.hold_curves(p, list(fleet_control = fc), last)
  later <- (last + 1):n_yrs

  # Every deviation stored as a level keeps its last value, IID included; the
  # step forms take a step of 0.
  for (yr in later) {
    testthat::expect_identical(h$index_q_dev[, yr], p$index_q_dev[, last])
    testthat::expect_identical(h$sel_inf_dev[, , , yr], p$sel_inf_dev[, , , last])
    testthat::expect_identical(h$sel_coff_dev[1:2, , , yr], p$sel_coff_dev[1:2, , , last])
    testthat::expect_true(all(h$sel_coff_dev[3:4, , , yr] == 0))
    testthat::expect_identical(h$log_M1_dev[, , , yr], p$log_M1_dev[, , , last])
  }
  testthat::expect_identical(h$log_M1_dev[, , , n_yrs + 1:2], p$log_M1_dev[, , , n_yrs + 1:2])
  # The years with data are untouched.
  testthat::expect_identical(h$index_q_dev[, 1:last], p$index_q_dev[, 1:last])
  testthat::expect_identical(h$sel_coff_dev[, , , 1:last], p$sel_coff_dev[, , , 1:last])
})
