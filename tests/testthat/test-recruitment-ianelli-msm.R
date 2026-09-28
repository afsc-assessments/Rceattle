# Stock-recruit curves under predation (msmMode > 0), where SPR is undefined.
# The Ianelli penalty and a hindcast curve (with a free R_init) are allowed; a BH
# steepness prior is refused. Before 5.30.0 all were refused.

msm_srr_data <- function() {
  set.seed(123)
  make_msm_test_data()$data_list
}

msm_srr_build <- function(recFun, d = msm_srr_data(), map = NULL) {
  suppressMessages(suppressWarnings(fit_mod(
    data_list = d, inits = NULL, map = map, estimateMode = 3, msmMode = 1,
    suitMode = 0, initMode = "NonEquilibrium", random_rec = FALSE,
    recFun = recFun,
    fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE))))
}

srr_penalty_row <- function(m) {
  jc <- m$quantities$jnll_comp
  jc[grep("Stock-recruit penalty", rownames(jc)), ]
}

test_that("the Ianelli Beverton-Holt penalty builds under predation", {
  m <- msm_srr_build(build_srr(srr_fun = "mean", srr_pred_fun = "BevertonHolt",
                               srr_est_mode = "Estimated"))
  q <- m$quantities

  expect_true(is.finite(m$obj$fn()))
  expect_true(all(is.finite(m$obj$gr())))

  # Year 1 sits on R_init, not the undefined SPRFinit.
  expect_true(all(is.finite(q$R_hat)))
  expect_equal(as.numeric(q$R_hat[, 1]), as.numeric(q$R_init), tolerance = 1e-12)

  # The deviation sim_mod() and retrospective() resample.
  nyrs_hind <- m$data_list$endyr - m$data_list$styr + 1
  dev <- (log(q$R) - log(q$R_hat))[, seq_len(nyrs_hind)]
  expect_true(all(is.finite(dev)))

  expect_true(all(q$steepness == 0))
  expect_equal(as.numeric(q$SPR0), rep(0, m$data_list$nspp))
})

test_that("a Ricker penalty with a prior on alpha builds under predation", {
  m <- msm_srr_build(build_srr(srr_fun = "mean", srr_pred_fun = "Ricker",
                               srr_est_mode = "LognormalPrior",
                               srr_prior = 4, srr_prior_sd = 1))
  expect_true(is.finite(m$obj$fn()))
  expect_true(all(is.finite(m$obj$gr())))
  expect_true(all(is.finite(m$quantities$R_hat)))
})

test_that("a steepness prior is refused under predation", {
  for (f in c("mean", "BevertonHolt")) {
    expect_error(msm_srr_build(build_srr(srr_fun = f, srr_pred_fun = "BevertonHolt",
                                         srr_est_mode = "LognormalPrior",
                                         srr_prior = 0.8, srr_prior_sd = 0.2)),
                 "steepness", info = f)
    expect_error(msm_srr_build(build_srr(srr_fun = f, srr_pred_fun = "BevertonHolt",
                                         srr_est_mode = "BetaPrior",
                                         srr_prior = 0.8, srr_prior_sd = 0.2)),
                 "steepness", info = f)
  }
})

test_that("data_check() resolves string switches before comparing them", {
  # "mean" >= 2 is TRUE in R; "2" is a code read in as text.
  d <- msm_srr_data()
  d$msmMode <- 1
  srr_msg <- function(d) {
    msg <- tryCatch({ Rceattle:::data_check(d); "" },
                    error = function(e) conditionMessage(e))
    grepl("stock-recruit|steepness|srr_fun", msg)
  }

  d$srr_fun <- "mean"; d$srr_pred_fun <- "BevertonHolt"; d$srr_est_mode <- "Estimated"
  expect_false(srr_msg(d))

  d$srr_est_mode <- "LognormalPrior"
  expect_true(srr_msg(d))

  d$srr_fun <- "BevertonHolt"; d$srr_est_mode <- "Estimated"
  expect_false(srr_msg(d))

  d$srr_fun <- "2"
  expect_false(srr_msg(d))

  d$srr_fun <- "0"; d$srr_pred_fun <- "2"; d$srr_est_mode <- "2"
  expect_true(srr_msg(d))
})

test_that("convergence_diagnostics() reads the curve against the data, not steepness, under predation", {
  m <- msm_srr_build(build_srr(srr_fun = "mean", srr_pred_fun = "BevertonHolt",
                               srr_est_mode = "Estimated"))
  sr <- convergence_diagnostics(m)$checks$stock_recruit
  expect_true(!is.null(sr))
  # Steepness is never the test here; the curve is read over the SSB range
  # (test-recruitment-srr-degenerate-check.R). At the default starts (alpha
  # e^3, beta 3) both curves are flat over the fixture's SSB.
  expect_identical(sr$severity, "WARN")
  expect_match(sr$message, "Species1: flat over the observed SSB range")
  expect_match(sr$message, "Species2: flat over the observed SSB range")
  expect_false(grepl("steepness", sr$message))
})

test_that("the Ianelli penalty fits under predation from on-scale starts", {
  # The 5.31.0 tests above are build-only. One fitted case: the penalty form
  # converges on the fixture when alpha and beta start on the stock's scale
  # (build_srr(srr_alpha_init =, srr_beta_init =)), and the curve check reads
  # the fitted curve rather than skipping it.
  skip_on_cran()
  d  <- msm_srr_data()
  m0 <- msm_srr_build(build_srr(), d)
  hind <- d$endyr - d$styr + 1
  ssb  <- m0$quantities$ssb[, seq_len(hind)]
  R    <- m0$quantities$R[, seq_len(hind)]
  s_med <- apply(ssb, 1, stats::median)
  fit <- suppressMessages(suppressWarnings(fit_mod(
    data_list = d, estimateMode = "Hindcast", msmMode = 1, suitMode = 0,
    initMode = "NonEquilibrium", random_rec = FALSE,
    recFun = build_srr(srr_fun = "mean", srr_pred_fun = "BevertonHolt",
                       srr_alpha_init = 2 * rowMeans(R) / s_med, srr_beta_init = 1 / s_med),
    fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE))))
  expect_true(is.finite(fit$opt$objective))
  expect_lt(max(abs(fit$obj$gr())), 0.1)
  expect_true(all(is.finite(fit$estimated_params$rec_pars[, 2:3])))
  expect_true(all(abs(fit$estimated_params$rec_pars[, 2:3]) < 30))
  # What this fit produces: species 1 bends inside the data (predicted /
  # asymptote 0.73 at the lowest SSB, 0.90 at the highest), species 2 runs to
  # the flat ridge (alpha 2.4e6, beta 1251, ratio 1.000 everywhere).
  sr <- convergence_diagnostics(fit)$checks$stock_recruit
  expect_identical(sr$severity, "WARN")
  expect_match(sr$message, "Species2: flat over the observed SSB range")
  expect_false(grepl("Species1:", sr$message))
  expect_lt(sr$data$Species1$dd_at_smin, 0.9)
  expect_gt(sr$data$Species2$dd_at_smin, 0.9)
})

test_that("a fixed-dynamics species carries no curve penalty", {
  # Species 1 is input (estDynamics = 2); ungated, its penalty was ~1100 nats.
  d  <- msm_srr_data()
  m0 <- msm_srr_build(build_srr(), d)
  N  <- m0$quantities$N_at_age
  yrs <- d$styr:d$endyr
  nb <- do.call(rbind, lapply(seq_along(yrs), function(i)
    data.frame(Species_name = "Species1", Species = 1, Sex = 0, Year = yrs[i],
               t(N[1, 1, , i]))))
  colnames(nb) <- c("Species_name", "Species", "Sex", "Year",
                    paste("Age", seq_len(dim(N)[3])))
  d$NByageFixed <- nb
  d$estDynamics <- c(2, 0)

  m <- msm_srr_build(build_srr(srr_fun = "mean", srr_pred_fun = "BevertonHolt",
                               srr_est_mode = "Estimated"), d)
  pen <- srr_penalty_row(m)
  expect_equal(unname(pen[1]), 0)
  expect_gt(unname(pen[2]), 0)
  expect_true(is.finite(m$obj$fn()))
  expect_true(all(is.finite(m$obj$gr())))
})

test_that("a supplied map that fixes the curve warns", {
  d <- msm_srr_data()
  mean_map <- msm_srr_build(build_srr(), d)$map
  expect_warning(
    suppressMessages(fit_mod(
      data_list = d, inits = NULL, map = mean_map, estimateMode = 3, msmMode = 1,
      suitMode = 0, initMode = "NonEquilibrium", random_rec = FALSE,
      recFun = build_srr(srr_fun = "mean", srr_pred_fun = "BevertonHolt"),
      fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE))),
    "fixes both stock-recruit parameters")
})

test_that("build_srr() refuses a Ricker beta prior and flags an ignored Bmsy_lim", {
  expect_error(build_srr(srr_fun = "mean", srr_pred_fun = "Ricker",
                         srr_est_mode = "BetaPrior",
                         srr_prior = 0.8, srr_prior_sd = 0.2),
               "no Ricker form")
  expect_warning(build_srr(srr_fun = "mean", srr_pred_fun = "BevertonHolt",
                           Bmsy_lim = 1e6),
                 "Bmsy_lim")
  expect_no_warning(build_srr(srr_fun = "mean", srr_pred_fun = "Ricker",
                              Bmsy_lim = 1e6))
  # -999 is the stored "off" value a refit passes back in.
  expect_no_warning(build_srr(srr_fun = "mean", srr_pred_fun = "BevertonHolt",
                              Bmsy_lim = -999))
})

test_that("a hindcast curve under predation builds in every initMode", {
  modes <- c("FreeParams", "Equilibrium", "NonEquilibrium", "FishedNonEquilibrium",
             "FishedNonEquilibriumScaled", "OffsetEquilibrium")
  cases <- rbind(data.frame(f = "BevertonHolt", im = modes),
                 data.frame(f = "Ricker", im = c("FreeParams", "NonEquilibrium")))
  for (k in seq_len(nrow(cases))) {
    lab <- paste(cases$f[k], cases$im[k])
    m <- suppressMessages(suppressWarnings(fit_mod(
      data_list = msm_srr_data(), inits = NULL, estimateMode = 3, msmMode = 1,
      suitMode = 0, initMode = cases$im[k], random_rec = FALSE,
      recFun = build_srr(srr_fun = cases$f[k]),
      fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE))))
    q <- m$quantities
    expect_true(is.finite(m$obj$fn()), info = lab)
    expect_true(all(is.finite(m$obj$gr())), info = lab)
    expect_true(all(is.finite(q$R)) && all(is.finite(q$R_hat)), info = lab)

    # R_init is the free level exp(rec_pars[, "R0"]), not an SPR equilibrium.
    expect_false(anyNA(m$map$mapList$rec_pars[, 1]), info = lab)
    expect_equal(as.numeric(q$R_init), as.numeric(exp(m$estimated_params$rec_pars[, 1])),
                 tolerance = 1e-12, info = lab)
    expect_true(all(q$steepness == 0), info = lab)

    # sample_rec() projects the mean ratio to the curve, not log(mean R / R0).
    s  <- sample_rec(m, sample_rec = FALSE, update_model = FALSE)
    nh <- m$data_list$endyr - m$data_list$styr + 1
    expect_equal(unname(s$estimated_params$rec_dev[, nh + 1]),
                 unname(log(rowMeans((q$R / q$R_hat)[, 1:nh]))), tolerance = 1e-12,
                 info = lab)
  }
})

test_that("the Ricker positivity penalty is skipped under predation", {
  # posfun(alpha * SPR0 - 1) with SPR0 = 0 would add ~0.02 per species.
  m <- msm_srr_build(build_srr(srr_fun = "Ricker"))
  zn <- m$quantities$jnll_comp[grep("Zero n-at-age", rownames(m$quantities$jnll_comp)), 1:2]
  expect_lt(max(zn), 0.005)
})
