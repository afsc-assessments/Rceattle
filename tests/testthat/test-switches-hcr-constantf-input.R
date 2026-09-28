# ConstantF projects every species at its input Ftarget, recycled from a single
# value. Up to 5.33.0 the multispecies HCR loop reset log_Ftarget to 0 (F = 1),
# and a single Ftarget with nspp > 1 stopped on a map-size error.

.constantf_fit <- function(d, Ftarget, msmMode) {
  base <- suppressMessages(suppressWarnings(fit_mod(
    data_list = d, inits = NULL, estimateMode = 3, msmMode = msmMode, random_rec = FALSE,
    fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE))))
  fit <- suppressMessages(suppressWarnings(fit_mod(
    data_list = d, inits = base$estimated_params, estimateMode = 2, msmMode = msmMode,
    random_rec = FALSE, HCR = build_hcr(HCR = "ConstantF", Ftarget = Ftarget),
    fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE))))
  fit
}
.proj_cols <- function(d) (d$endyr + 1):d$projyr - d$styr + 1

test_that("multispecies ConstantF projects at the input F", {
  testthat::skip_on_cran()
  set.seed(123)
  d    <- make_msm_test_data()$data_list
  proj <- .proj_cols(d)
  f02  <- .constantf_fit(d, Ftarget = 0.2, msmMode = 1)          # single value, recycled
  expect_equal(unname(f02$quantities$F_spp[, proj]), matrix(0.2, 2, length(proj)), tolerance = 1e-8)
  f13  <- .constantf_fit(d, Ftarget = c(0.1, 0.3), msmMode = 1)  # one per species
  expect_equal(unname(f13$quantities$F_spp[, proj[1]]), c(0.1, 0.3), tolerance = 1e-8)
  # F-at-age, which removes fish, scales with each species' input F (0.1/0.2, 0.3/0.2).
  faa <- function(fit, sp) sum(fit$quantities$F_at_age[sp, 1, , proj[1]])
  expect_equal(faa(f13, 1) / faa(f02, 1), 0.5, tolerance = 1e-8)
  expect_equal(faa(f13, 2) / faa(f02, 2), 1.5, tolerance = 1e-8)
  expect_error(.constantf_fit(d, Ftarget = -0.1, msmMode = 1), "non-negative")
  # An input F of 0 is the model's no-fishing value on the log scale, not -Inf.
  f00  <- .constantf_fit(d, Ftarget = c(0, 0.2), msmMode = 1)
  expect_true(all(f00$quantities$F_spp[1, proj] == 0))
  expect_equal(unname(f00$quantities$F_spp[2, proj[1]]), 0.2, tolerance = 1e-8)
  expect_true(all(is.finite(f00$estimated_params$log_Ftarget)))
})

test_that("single-species ConstantF projects at the input F", {
  testthat::skip_on_cran()
  d <- make_test_data()
  f <- .constantf_fit(d, Ftarget = 0.2, msmMode = 0)
  expect_equal(unname(f$quantities$F_spp[1, .proj_cols(d)]), rep(0.2, length(.proj_cols(d))),
               tolerance = 1e-8)
  # The projection is unfished either way; the finite check is what guards the fix.
  f0 <- .constantf_fit(d, Ftarget = 0, msmMode = 0)
  expect_true(all(f0$quantities$F_spp[1, .proj_cols(d)] == 0))
  expect_true(is.finite(f0$estimated_params$log_Ftarget))
  # Its inits must not pin a rule that estimates Ftarget at F = 0.
  f40 <- suppressMessages(suppressWarnings(fit_mod(
    data_list = d, inits = f0$estimated_params, estimateMode = 2, msmMode = 0,
    random_rec = FALSE, HCR = build_hcr(HCR = "NPFMC", Ftarget = 0.4, Flimit = 0.35),
    fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE))))
  expect_gt(unname(f40$quantities$Ftarget[1]), 0.01)
})
