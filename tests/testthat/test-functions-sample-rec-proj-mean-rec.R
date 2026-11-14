# Under build_srr(proj_mean_rec = TRUE), the default, projected recruitment is
# avg_R and never reads the projection deviations sample_rec() writes. Before
# 5.67.8 sample_rec(update_model = TRUE) returned that unchanged projection
# silently. run_mse() calls it with update_model = FALSE, so it stays silent
# there.
testthat::test_that("sample_rec() warns only when its refit ignores the deviations", {
  testthat::skip_if_not_installed("TMB")

  d <- make_test_data(seed = 1)
  build <- function(proj_mean_rec) {
    suppressMessages(suppressWarnings(Rceattle::fit_mod(
      d, file = NULL, estimateMode = 3, msmMode = 0,
      recFun = Rceattle::build_srr(proj_mean_rec = proj_mean_rec),
      fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0))))
  }
  nyrs_hind <- d$endyr - d$styr + 1
  proj <- (nyrs_hind + 1):(d$projyr - d$styr + 1)

  # Default: warns, and projected recruitment does not move.
  m_mean <- build(TRUE)
  testthat::expect_warning(
    s_mean <- Rceattle::sample_rec(m_mean, rec_trend = -0.5), "do not drive it")
  testthat::expect_equal(s_mean$quantities$R[, proj], m_mean$quantities$R[, proj])

  # run_mse()'s call writes the deviations without warning.
  testthat::expect_no_warning(
    Rceattle::sample_rec(m_mean, rec_trend = -0.5, update_model = FALSE))

  # proj_mean_rec = FALSE: silent, and the 50% decline reaches projected
  # recruitment. The hindcast deviations are 0 here, so resampling adds none.
  m_dev <- build(FALSE)
  testthat::expect_no_warning(
    s_dev <- Rceattle::sample_rec(m_dev, rec_trend = -0.5))
  trend <- 1 - 0.5 * seq_along(proj) / length(proj)
  testthat::expect_equal(
    as.vector(s_dev$quantities$R[, proj] / m_dev$quantities$R[, proj]),
    rep(trend, each = d$nspp))
})
