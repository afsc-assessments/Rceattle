# run_mse(regenerate_past = TRUE) refits the EM to past data regenerated from the
# OM. Under simulate_data = TRUE those data carry observation error, drawn once
# under regenerate_seed and shared by every simulation. Before 5.70.4 they were
# the OM's expected values whatever simulate_data was.

testthat::test_that("regenerate_past draws observation error under simulate_data", {
  testthat::skip_on_cran()
  data(BS2017SS)
  d <- BS2017SS
  d$projyr <- d$endyr + 2L
  d$fleet_control$Proj_F_proportion <- rep(1, nrow(d$fleet_control))
  om <- suppressWarnings(suppressMessages(fit_mod(
    d, estimateMode = 1, msmMode = 0,
    fit_control = fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))

  hindcast_obs <- function(index_data) {
    index_data$Observation[index_data$Year > 0 & index_data$Year <= d$endyr]
  }
  regenerated <- function(regenerate_seed) {
    mse <- suppressWarnings(suppressMessages(run_mse(
      om = om, em = om, nsim = 1, assessment_period = 2, sampling_period = 1,
      simulate_data = TRUE, regenerate_past = TRUE, sample_rec = FALSE,
      regenerate_seed = regenerate_seed, cores = 1, seed = 1)))
    hindcast_obs(mse$Sim_1$EM[[1]]$data_list$index_data)
  }

  expected <- hindcast_obs(sim_mod(om, simulate = FALSE)$index_data)
  draw_1 <- regenerated(1)
  testthat::expect_length(draw_1, length(expected))
  # Drawn with observation error, not the OM's expected values ...
  testthat::expect_false(isTRUE(all.equal(draw_1, expected)))
  # ... and regenerate_seed sets the draw.
  testthat::expect_identical(regenerated(1), draw_1)
  testthat::expect_false(isTRUE(all.equal(draw_1, regenerated(2))))
})
