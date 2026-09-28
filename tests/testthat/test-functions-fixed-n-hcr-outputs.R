# A species with input numbers-at-age (estDynamics > 0) has no harvest control rule:
# it is projected at F = 0 and its F-based reference points are NA (5.34.0).

.fixed_n_fit <- function() {
  set.seed(123)
  d <- make_msm_test_data()$data_list
  build <- function(d) suppressMessages(suppressWarnings(fit_mod(
    data_list = d, inits = NULL, estimateMode = 3, msmMode = 1, suitMode = 0,
    initMode = "NonEquilibrium", random_rec = FALSE,
    fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE))))
  N   <- build(d)$quantities$N_at_age
  yrs <- d$styr:(d$styr + dim(N)[4] - 1)
  nb  <- do.call(rbind, lapply(seq_along(yrs), function(i)
    data.frame(Species_name = "Species1", Species = 1, Sex = 0, Year = yrs[i],
               t(N[1, 1, , i]))))
  colnames(nb) <- c("Species_name", "Species", "Sex", "Year",
                    paste("Age", seq_len(dim(N)[3])))
  d$NByageFixed <- nb
  d$estDynamics <- c(1, 0)
  build(d)
}

test_that("reference points set by the unestimated F are NA for that species", {
  testthat::skip_on_cran()
  fit <- .fixed_n_fit()
  q   <- fit$quantities
  for (nm in c("Ftarget", "Flimit", "SPRtarget", "SPRlimit")) {
    expect_true(is.na(q[[nm]][1]), info = nm)
    expect_false(is.na(q[[nm]][2]), info = nm)
  }
  for (nm in c("SBF", "DynamicSBF")) {
    expect_true(all(is.na(q[[nm]][1, ])), info = nm)
    expect_true(all(is.finite(q[[nm]][2, ])), info = nm)
  }
  # Depletion is still reported for the fixed species.
  expect_true(all(is.finite(q$ssb_depletion[1, ])))
})

test_that("a species with input numbers-at-age is projected at F = 0", {
  testthat::skip_on_cran()
  fit0 <- .fixed_n_fit()
  d    <- fit0$data_list
  proj <- (d$endyr + 1):d$projyr - d$styr + 1
  # Projection only. PFMC with Ptarget near 1 resets projected F by stock status
  # (section 6.7) in every year between Plimit and Ptarget.
  hcrs <- list(ConstantF = build_hcr(HCR = "ConstantF", Ftarget = rep(0.2, d$nspp)),
               PFMC = build_hcr(HCR = "PFMC", Flimit = 0.45, Ptarget = 0.999,
                                Plimit = 0.001, Pstar = 0.45, Sigma = 0.5))
  for (nm in names(hcrs)) {
    fit <- suppressMessages(suppressWarnings(fit_mod(
      data_list = d, inits = fit0$estimated_params, estimateMode = 2, msmMode = 1,
      suitMode = 0, initMode = "NonEquilibrium", random_rec = FALSE, HCR = hcrs[[nm]],
      fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE))))
    expect_true(all(fit$quantities$F_spp[1, proj] == 0), info = nm)
    if (nm == "ConstantF") expect_gt(min(fit$quantities$F_spp[2, proj]), 1e-3)
    if (nm == "PFMC") {
      # Precondition: the fixed species' depletion sits where the reset fires.
      dep <- fit$quantities$ssb_depletion[1, proj - 1]
      expect_true(any(dep > 0.001 & dep < 0.999))
    }
  }
})

test_that("mse_summary() gives NA, not NaN, overfishing probabilities for that species", {
  testthat::skip_on_cran()
  fit  <- .fixed_n_fit()
  mse  <- list(Sim_1 = list(EM = list(EM = fit), use_sim = TRUE, failure = NA,
                            OM = fit, OM_no_F = suppressMessages(suppressWarnings(remove_F(fit)))))
  summ <- suppressWarnings(mse_summary(mse, om_only = TRUE))
  p    <- summ$species$om_p_overfishing
  expect_true(is.na(p[1]) && !is.nan(p[1]))
  expect_false(is.na(p[2]))
  for (cc in c("catch_iav", "p_closed")) {
    expect_true(is.na(summ$species[[cc]][1]) && !is.nan(summ$species[[cc]][1]), info = cc)
  }
  expect_true(is.finite(summ$species$om_terminal_depletion_dynamic[1]))
})
