# Length-based suitability (suitMode 1, 3, 5) takes the predator:prey size ratio
# from length_hat, whose rows run as weight_hat's do: population length at row
# 2 * sp, spawning-time length at 2 * sp + 1 (0-based), then one row per fleet.
# Before 5.66.15 the template read row `sp`, so the second species of two took
# the first species' spawning-time lengths. On this fixture at its true growth
# that put suitability off by up to 0.38 and moved the objective by 0.4% to 5%.
# data_check() refuses these modes, so it is mocked out to reach the template.

.length_suit_fit <- function(suitMode) {
  set.seed(20261105)
  sim <- make_msm_test_data()
  d <- sim$data_list
  ctl <- list(data_list = d, estimateMode = 3, msmMode = 1, niter = 3,
              growthFun = build_growth(fun = "vonBertalanffy"),
              fit_control = fit_control(getsd = FALSE, verbose = 0))
  start <- suppressMessages(suppressWarnings(
    do.call(fit_mod, c(ctl, suitMode = 2))))$estimated_params
  # The fixture's own growth (K, L1, Linf; cm), so the species differ in length.
  start$log_growth_pars[, 1, 1:3] <- log(c(0.3, 0.35, 4.5, 4.5, 90, 50))
  suppressMessages(suppressWarnings(
    do.call(fit_mod, c(ctl, list(suitMode = suitMode, inits = start)))))
}

# Suitability rebuilt from the population-length rows (R row 2 * sp - 1).
.length_suit_oracle <- function(q, suitMode) {
  L <- q$length_hat
  S <- array(0, dim(q$suitability))
  nspp <- nrow(q$vulnerability)
  for (r in seq_len(nspp)) for (k in seq_len(nspp)) {
    a <- q$gam_a[r]; b <- q$gam_b[r]
    for (ra in seq_len(dim(L)[3])) {
      x <- sweep(1 / L[2 * k - 1, 1, , ], 2, L[2 * r - 1, 1, ra, ], `*`)
      if (suitMode %in% c(1, 3)) x <- log(x)
      dens <- if (suitMode == 1) {
        dgamma(x, a, scale = b) / dgamma((a - 1) * b, a, scale = b)
      } else {
        dnorm(x, a, b) / dnorm(a, a, b)
      }
      S[r, k, ra, , ] <- ifelse(x > 0, q$vulnerability[r, k] * dens, 0)
    }
  }
  S
}

testthat::test_that("data_check() still refuses length-based suitability", {
  set.seed(20261105)
  d <- make_msm_test_data()$data_list
  d$suitMode <- rep(1, d$nspp)
  testthat::expect_error(suppressWarnings(data_check(d)),
                         "Length-based suitability")
})

testthat::test_that("length-based suitability reads each species' population length", {
  testthat::local_mocked_bindings(data_check = function(...) invisible(TRUE),
                                  .package = "Rceattle")
  for (m in c(1, 3, 5)) {
    q <- .length_suit_fit(m)$quantities
    # Non-vacuous: species 2's population row (3) must differ from row 2, the
    # first species' spawning-time row that a species-number index reads.
    testthat::expect_gt(max(abs(q$length_hat[2, 1, , ] - q$length_hat[3, 1, , ])), 10)
    testthat::expect_equal(unname(q$suitability),
                           .length_suit_oracle(q, m), tolerance = 1e-12)
  }
})
