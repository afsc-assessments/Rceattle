# Projected F_flt is each fishery's share of its species' projected F; surveys are 0.
# Up to 5.33.0 section 6.7 wrote it by species index, in every fit.

test_that("projected F_flt is indexed by fleet", {
  testthat::skip_on_cran()
  set.seed(123)
  d    <- make_msm_test_data()$data_list
  base <- suppressMessages(suppressWarnings(fit_mod(
    data_list = d, inits = NULL, estimateMode = 3, msmMode = 1, random_rec = FALSE,
    fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE))))
  fit  <- suppressMessages(suppressWarnings(fit_mod(
    data_list = d, inits = base$estimated_params, estimateMode = 2, msmMode = 1,
    random_rec = FALSE, HCR = build_hcr(HCR = "ConstantF", Ftarget = c(0.1, 0.3)),
    fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE))))
  q    <- fit$quantities
  fc   <- d$fleet_control
  proj <- (d$endyr + 1):d$projyr - d$styr + 1
  fishery <- fc$Fleet_type %in% c("Fishery", 1)
  for (flt in fc$Fleet_code[fishery]) {
    sp  <- fc$Species[flt]
    exp <- fc$Proj_F_proportion[flt] * q$F_spp[sp, proj]
    expect_equal(unname(q$F_flt[flt, proj]), unname(exp), tolerance = 1e-10,
                 info = paste("fleet", flt))
  }
  expect_true(all(q$F_flt[fc$Fleet_code[!fishery], proj] == 0))
})
