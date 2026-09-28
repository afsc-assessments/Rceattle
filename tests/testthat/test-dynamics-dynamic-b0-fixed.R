# Dynamic B0 keeps a species' input numbers-at-age (estDynamics > 0). Before
# 5.30.0 it projected them on placeholder recruitment, so a fixed predator
# collapsed (hake 4-spp: arrowtooth dynamic SSB 38.6 mt vs input 43,490 mt).
# fixed_dyn_build() is in helpers-fixed-natage.R.

test_that("a fixed-dynamics species keeps its input numbers in dynamic B0", {
  m <- fixed_dyn_build(c(2, 0))
  q <- m$quantities
  nh <- m$data_list$endyr - m$data_list$styr + 1

  # Same numbers, same weights: dynamic B0 equals biomass exactly.
  expect_equal(unname(q$DynamicB0[1, seq_len(nh)]), unname(q$biomass[1, seq_len(nh)]),
               tolerance = 1e-10)

  expect_true(all(is.finite(q$DynamicSB0)))
  expect_true(all(q$DynamicSB0[1, ] > 0))
  expect_true(is.finite(m$obj$fn()))
  expect_true(all(is.finite(m$obj$gr())))
})

test_that("an estimated species' dynamic B0 is untouched by the fixed-species branch", {
  set.seed(123)
  d <- make_msm_test_data()$data_list
  m <- suppressMessages(suppressWarnings(fit_mod(
    data_list = d, inits = NULL, estimateMode = 3, msmMode = 1, suitMode = 0,
    initMode = "NonEquilibrium", random_rec = FALSE,
    fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE))))
  q <- m$quantities
  # Year 1 is the shared starting state.
  expect_equal(unname(q$DynamicB0[, 1]), unname(q$biomass[, 1]), tolerance = 1e-10)
  expect_true(all(is.finite(q$DynamicSB0)))
})
