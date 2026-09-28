# estDynamics = 3 ("FixedScaledByAge") was retired in 5.35.0. It never estimated
# an age-specific multiplier: build_map() freed only the first age's scalar
# under predation and none in single-species mode, so every such fit was the
# estDynamics = 2 model. Measured on the fixed-numbers fixture before the
# retirement (5.34.0): codes 2 and 3 both gave objective 222899009.1243875 with
# one free scalar; code 1 gave the same objective with none. (5.35.0's plus-group
# survival weighting moved this parametric-growth fixture's objective, whose
# second species has M = 0.3, so the test pins codes 1 and 2 to each other.)

testthat::test_that("estDynamics = 3 is refused and names the code it fitted as", {
  d <- Rceattle::BS2017SS
  for (v in list(3, "3", "FixedScaledByAge")) {
    d$estDynamics <- c(v, 0, 0)
    testthat::expect_error(suppressMessages(Rceattle::switch_check(d)),
                           "retired in 5.35.0.*Set estDynamics = 2", info = v)
  }
})

testthat::test_that("an out-of-range estDynamics code is refused in R, not in the template", {
  d <- Rceattle::BS2017SS
  d$estDynamics <- c(7, 0, 0)
  testthat::expect_error(suppressMessages(Rceattle::switch_check(d)),
                         "Invalid 'estDynamics' value\\(s\\): 7")
})

testthat::test_that("log_pop_scalar is one value per species and code 2 fits unchanged", {
  testthat::skip_on_cran()
  m2 <- fixed_dyn_build(c(2, 0))
  m1 <- fixed_dyn_build(c(1, 0))
  testthat::expect_equal(m2$obj$fn(), m1$obj$fn(), tolerance = 1e-12)
  testthat::expect_equal(sum(names(m2$obj$par) == "log_pop_scalar"), 1L)
  testthat::expect_equal(sum(names(m1$obj$par) == "log_pop_scalar"), 0L)
  testthat::expect_null(dim(m2$quantities$pop_scalar))
  testthat::expect_named(m2$quantities$pop_scalar, m2$data_list$spnames)
})

testthat::test_that("an older fit's [nspp, nages] log_pop_scalar warm-starts as its first column", {
  testthat::skip_on_cran()
  m2 <- fixed_dyn_build(c(2, 0))
  inits <- m2$estimated_params
  nspp  <- m2$data_list$nspp
  old   <- matrix(0, nspp, max(m2$data_list$nages))
  old[1, 1] <- 0.25
  inits$log_pop_scalar <- old
  m_old <- suppressMessages(suppressWarnings(fit_mod(
    data_list = m2$data_list, inits = inits, estimateMode = 3, msmMode = 1,
    suitMode = 0, initMode = "NonEquilibrium", random_rec = FALSE,
    fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE))))
  testthat::expect_equal(unname(m_old$estimated_params$log_pop_scalar), c(0.25, 0))
  testthat::expect_equal(unname(m_old$quantities$pop_scalar[1]), exp(0.25))
})
