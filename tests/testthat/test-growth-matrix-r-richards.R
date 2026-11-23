# get_growth_matrix_r()'s Richards cohort step after year 1 must match the Richards
# branch of estimate_growth() in src/TMB/growth.hpp. It read L1 where K belongs and
# divided by m instead of taking the 1/m root; found by the PR #346 attribution review.

# One sex, 10 ages, L1 = 20 cm at age 1, 5-cm bins from 0 cm; parameters vary by year.
growth_r <- function(model, m, K = c(0.30, 0.20, 0.25), linf = c(90, 80, 85)) {
  nyrs <- length(K)
  gp <- array(0, dim = c(1, nyrs, 4))
  gp[1, , 1] <- K
  gp[1, , 2] <- 20
  gp[1, , 3] <- linf
  gp[1, , 4] <- m
  Rceattle:::get_growth_matrix_r(
    fracyr = 0, nsex_sp = 1, nages_sp = 10, nlengths_sp = 21, nyrs = nyrs,
    lengths_sp = seq(0, 100, by = 5), minage_sp = 1, maxage_sp = 10,
    growth_params_sp = gp, growth_log_sd_sp = matrix(log(c(2, 6)), 1, 2),
    growth_model_sp = model, M_plus_sp = 0.3
  )
}

testthat::test_that("Richards with m = 1 reproduces von Bertalanffy after year 1", {
  vb <- growth_r(model = 1, m = 1)
  ri <- growth_r(model = 2, m = 1)
  testthat::expect_equal(ri$length_at_age, vb$length_at_age)
  testthat::expect_equal(ri$growth_matrix, vb$growth_matrix)
})

testthat::test_that("Richards cohort step uses last year's K, Linf and m", {
  K <- c(0.30, 0.20, 0.25)
  linf <- c(90, 80, 85)
  m <- c(1.3, 0.8, 1.1)
  out <- growth_r(model = 2, m = m, K = K, linf = linf)

  # growth.hpp: L[a, y] = (L[a-1, y-1]^m + (L[a-1, y-1]^m - Linf^m) * (exp(-K) - 1))^(1/m),
  # with year y-1's K, Linf and m; year 1 is the curve anchored at L1 at age 1.
  step <- function(l, y) (l^m[y] + (l^m[y] - linf[y]^m[y]) * (exp(-K[y]) - 1))^(1 / m[y])
  curve1 <- (linf[1]^m[1] + (20^m[1] - linf[1]^m[1]) * exp(-K[1] * (0:7)))^(1 / m[1])
  yr2 <- step(curve1, 1)          # ages 2-9 from ages 1-8
  yr3 <- step(c(20, yr2[1:7]), 2) # age 1 is L1 = 20 cm on the ramp

  # Ages 2-9: age 1 is on the ramp and age 10 is the plus group.
  testthat::expect_equal(unname(out$length_at_age[1, 2:9, 2]), yr2)
  testthat::expect_equal(unname(out$length_at_age[1, 2:9, 3]), yr3)
})
