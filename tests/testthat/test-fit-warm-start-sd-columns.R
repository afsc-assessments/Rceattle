# Time_varying_sel_sd and Time_varying_q_sd are read from fleet_control only on
# a fresh build: build_params() stores their log in sel_dev_log_sd /
# index_q_dev_log_sd, and a warm start from `inits` keeps the old value. Until
# 5.35.0 editing the column and refitting from a previous fit was a silent
# no-op, the same shape as the Comp_weights trap fit_mod() already warns about.

testthat::skip_on_cran()

testthat::test_that("editing a deviation-sd column under a warm start is warned about", {
  d <- make_test_data()
  if (!is.null(d$data_list)) d <- d$data_list
  build <- function(d, inits = NULL) fit_mod(
    data_list = d, inits = inits, estimateMode = 3, msmMode = 0, random_rec = FALSE,
    fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE))
  m0 <- suppressWarnings(suppressMessages(build(d)))
  inits <- m0$estimated_params

  d2 <- d
  d2$fleet_control$Time_varying_sel_sd[1] <- d$fleet_control$Time_varying_sel_sd[1] * 2
  testthat::expect_warning(suppressMessages(build(d2, inits)),
                           "Time_varying_sel_sd.*differs from the supplied `inits\\$sel_dev_log_sd`")

  d3 <- d
  d3$fleet_control$Time_varying_q_sd[1] <- d$fleet_control$Time_varying_q_sd[1] * 2
  testthat::expect_warning(suppressMessages(build(d3, inits)),
                           "Time_varying_q_sd.*differs from the supplied `inits\\$index_q_dev_log_sd`")

  # An unchanged column is silent: the comparison is on the parameter's own log scale.
  testthat::expect_no_warning(suppressMessages(build(d, inits)),
                              message = "Time_varying")
})
