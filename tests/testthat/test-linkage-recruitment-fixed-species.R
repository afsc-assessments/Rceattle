# A recruitment linkage on a species with input numbers-at-age (estDynamics > 0)
# fits nothing: rec_pars is mapped out and N-at-age is read from NByageFixed.
# Measured on the fixed-numbers fixture before the refusal (5.34.0): an
# intercept prior on species 1's mean recruitment moved the objective from
# 222899009.1243875 to 222899496.1049893 (486.98 nats, a constant) with no new
# free parameter, and a covariate slope on the same species was a free
# parameter with gradient exactly 0 whose value left the objective unchanged.

testthat::skip_on_cran()

testthat::test_that("a recruitment linkage on a fixed-numbers species is refused", {
  d <- fixed_natage_data(c(2, 0))
  d$env_data$EnvData <- round(stats::rnorm(nrow(d$env_data)), 3)

  rec_prior <- build_srr(srr_fun = "mean", linkages = list(
    R0 = linkage_spec(~ 1, species = 1,
                      priors = list("(Intercept)" = prior_normal(5000, 100)))))
  testthat::expect_error(fixed_natage_build(d, rec_prior),
                         "species Species1.*input numbers-at-age")

  rec_slope <- build_srr(srr_fun = "mean", linkages = list(
    R0 = linkage_spec(~ 1 + EnvData, species = 1, data = d$env_data)))
  testthat::expect_error(fixed_natage_build(d, rec_slope),
                         "species Species1.*input numbers-at-age")

  # A spec with no species expands to one row per species, the fixed one included.
  rec_all <- build_srr(srr_fun = "mean", linkages = list(
    R0 = linkage_spec(~ 1 + EnvData, data = d$env_data)))
  testthat::expect_error(fixed_natage_build(d, rec_all),
                         "species Species1.*species = c\\(2\\)")
})

testthat::test_that("the same linkage on the estimated species builds and scores its prior", {
  d <- fixed_natage_data(c(2, 0))
  rec <- build_srr(srr_fun = "mean", linkages = list(
    R0 = linkage_spec(~ 1, species = 2,
                      priors = list("(Intercept)" = prior_normal(5000, 100)))))
  m <- fixed_natage_build(d, rec)
  testthat::expect_true(is.finite(m$obj$fn()))
  testthat::expect_gt(m$quantities$jnll_comp["Linkage-table priors", 2], 0)
})
