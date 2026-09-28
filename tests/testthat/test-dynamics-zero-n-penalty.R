# The "Zero n-at-age penalty" row is the per-species sum of the posfun()
# excursions at the 0.001 floors (numbers-at-age, the Ricker intercept, and the
# recruitment floors under an identity-link linkage). Until 5.35.0 the
# accumulator behind the numbers-at-age floor was reset once per iteration and
# added as a running total, so every cell of every species added every earlier
# excursion. Measured on BS2017SS at estimateMode = 3 with species 1's log R0
# moved from +8 to -20 and species 2 and 3 untouched: species 2's row went from
# 0 to 1.42e-3.

testthat::skip_on_cran()

testthat::test_that("a floor excursion in one species leaves the other species' penalty at zero", {
  build <- function(inits = NULL) suppressWarnings(suppressMessages(fit_mod(
    data_list = Rceattle::BS2017SS, inits = inits, estimateMode = 3, msmMode = 0,
    random_rec = FALSE, fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE))))
  m0 <- build()
  inits <- m0$estimated_params
  inits$rec_pars[1, 1] <- -20
  m1 <- build(inits)
  pen0 <- m0$quantities$jnll_comp["Zero n-at-age penalty", ]
  pen1 <- m1$quantities$jnll_comp["Zero n-at-age penalty", ]
  testthat::expect_gt(pen1[1], pen0[1])
  # Species 2 and 3 share nothing with species 1's R0, so their rows do not move
  # (at the starting values species 3 already touches the floor, 1.2e-6).
  testthat::expect_equal(unname(pen1[2:3]), unname(pen0[2:3]))
  rec <- Rceattle:::.check_zero_n_penalty(m1)
  testthat::expect_true(rec$zero_n_penalty$severity %in% c("WARN", "FAIL"))
  testthat::expect_match(rec$zero_n_penalty$message, "Pollock")
  # An excursion of over a thousand fish is a FAIL; a few fish is a WARN.
  fake <- m1; fake$quantities$jnll_comp["Zero n-at-age penalty", 1] <- 0.01 * 2^2
  testthat::expect_equal(Rceattle:::.check_zero_n_penalty(fake)$zero_n_penalty$severity, "FAIL")
})
