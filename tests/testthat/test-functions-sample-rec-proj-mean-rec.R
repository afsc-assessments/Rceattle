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
    s_mean <- Rceattle::sample_rec(m_mean, rec_trend = -0.5),
    "do not drive it")
  testthat::expect_equal(s_mean$quantities$R[, proj],
                         m_mean$quantities$R[, proj])

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

# Before 5.69.2 the refit was build-only (estimateMode = 3), which never runs the
# harvest control rule, so projected F was 0 under every rule while the fit still
# named it: BS2017SS under ConstantF (0.2, 0.3, 0.05) projected F_spp 0, with no
# warning. The refit also rebuilt the map, freeing a hand-fixed cell.
testthat::test_that("sample_rec() projects under the fit's harvest control rule", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  d <- make_test_data(seed = 1)
  hind <- seq_len(d$endyr - d$styr + 1)
  proj <- (d$endyr + 1):d$projyr - d$styr + 1
  fit <- function(hcr, mode = 0, ...) {
    suppressMessages(suppressWarnings(Rceattle::fit_mod(
      d, file = NULL, estimateMode = mode, msmMode = 0, HCR = hcr, ...,
      recFun = Rceattle::build_srr(proj_mean_rec = FALSE),
      fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
  }
  cf_hcr <- Rceattle::build_hcr(HCR = "ConstantF", Ftarget = 0.2)

  # ConstantF: the input F every projection year, on the fit's own hindcast and map.
  cf <- fit(cf_hcr)
  set.seed(1)
  s_cf <- Rceattle::sample_rec(cf, rec_trend = -0.5)
  testthat::expect_equal(unname(s_cf$quantities$F_spp[1, proj]), rep(0.2, length(proj)),
                         tolerance = 1e-12)
  testthat::expect_equal(s_cf$quantities$ssb[, hind], cf$quantities$ssb[, hind], tolerance = 1e-12)
  testthat::expect_identical(s_cf$map, cf$map)
  testthat::expect_lt(s_cf$quantities$R[1, max(proj)], 0.5 * cf$quantities$R[1, max(proj)])

  # A projection-only fit reruns its rule too.
  m2 <- fit(cf_hcr, mode = 2, inits = cf$estimated_params)
  set.seed(1)
  s_m2 <- Rceattle::sample_rec(m2, rec_trend = -0.5)
  testthat::expect_equal(unname(s_m2$quantities$F_spp[1, proj]), rep(0.2, length(proj)),
                         tolerance = 1e-12)

  # A hindcast-only fit ran no projection, so it stays unfished, and its map keeps
  # a hand-fixed selectivity cell fixed.
  map <- cf$map
  cell <- which(!is.na(map$mapList$sel_inf))[1]
  map$mapList$sel_inf[cell] <- NA
  map$mapFactor$sel_inf <- factor(map$mapList$sel_inf)
  m1 <- fit(cf_hcr, mode = 1, inits = cf$estimated_params, map = map)
  set.seed(1)
  s_m1 <- Rceattle::sample_rec(m1, rec_trend = -0.5)
  testthat::expect_true(all(s_m1$quantities$F_spp[1, proj] == 0))
  testthat::expect_equal(s_m1$quantities$ssb[, hind], m1$quantities$ssb[, hind], tolerance = 1e-12)
  testthat::expect_identical(s_m1$map, m1$map)

  # Tier 3: F40% is per-recruit, so it is unchanged; projected F ramps down as lower
  # resampled recruitment, and the 50% decline, take SSB below B40%. A fit without
  # standard errors gets none from the refit.
  np <- fit(Rceattle::build_hcr(HCR = "NPFMC", Ftarget = 0.4, Flimit = 0.35))
  set.seed(1)
  s_np <- Rceattle::sample_rec(np, rec_trend = -0.5)
  testthat::expect_equal(s_np$quantities$Ftarget, np$quantities$Ftarget, tolerance = 1e-8)
  testthat::expect_equal(s_np$quantities$F_spp[1, proj[1]], np$quantities$F_spp[1, proj[1]],
                         tolerance = 1e-8)
  testthat::expect_gt(min(s_np$quantities$F_spp[1, proj]), 0)
  testthat::expect_lt(s_np$quantities$F_spp[1, max(proj)], np$quantities$F_spp[1, max(proj)])
  testthat::expect_null(s_np$sdrep)
})
