# jitter() and self_test() test the model that was fitted, so they refit under its map.
#
# Provenance: the 2026 GOA multispecies script hand-rolled its jitter and skipped
# self_test() because both released cod's fixed coefficients. On BS2017SS with one
# catchability fixed, jitter() returned the free model's objective and printed status: OK.

test_that(".refit_map() keeps the stored map except for a projection-only fit", {
  fit <- list(data_list = list(estimateMode = 0), map = list(tag = "STORED"))
  expect_identical(.refit_map(fit), fit$map)
  fit$data_list$estimateMode <- 3
  expect_identical(.refit_map(fit), fit$map)
  fit$data_list$estimateMode <- 2
  expect_null(.refit_map(fit))
})

test_that("the redraw warning reads linkage deviates by their own slots", {
  # M1 deviates are penalized (slots 1-8) and R0 deviates integrated (slots 9-16).
  M1 <- linkage_spec(~ 1 + (1 | Year), by = ~ species, init = list(sigma = 0.2), integrate = FALSE)
  R0 <- linkage_spec(~ 1 + (1 | Year), by = ~ species)
  fit <- suppressWarnings(suppressMessages(fit_mod(
    data_list = make_test_data(), estimateMode = 3, M1Fun = build_M1(M1_model = 1, linkages = list(M1 = M1)),
    recFun = build_srr(srr_fun = 0, linkages = list(R0 = R0)),
    fit_control = fit_control(verbose = 0, phase = FALSE, getsd = FALSE))))
  integrate <- .re_slot_routing(fit$data_list$linkage_table)$integrate
  expect_true(any(integrate) && any(!integrate))
  hold_first <- function(nm) {
    fit$map$mapFactor[[nm]][1] <- NA
    fit
  }
  # Only the penalized (M1) slots were drawn.
  m_only <- list(list(beta_linkage_re = numeric(length(integrate)), beta_linkage_re_drawn = !integrate))
  expect_warning(.warn_redrawn_holds(hold_first("beta_linkage_re_pen"), m_only), "beta_linkage_re_pen")
  expect_no_warning(.warn_redrawn_holds(hold_first("beta_linkage_re"), m_only))
  # Only the integrated (R0) slots were drawn.
  r_only <- list(list(beta_linkage_re = numeric(length(integrate)), beta_linkage_re_drawn = integrate))
  expect_warning(.warn_redrawn_holds(hold_first("beta_linkage_re"), r_only), "beta_linkage_re")
  expect_no_warning(.warn_redrawn_holds(hold_first("beta_linkage_re_pen"), r_only))
})

test_that("jitter() and self_test() keep a fixed parameter fixed", {
  testthat::skip_on_cran()
  fit_model <- function(...) suppressWarnings(suppressMessages(fit_mod(
    data_list = make_test_data(seed = 1), estimateMode = 1, ...,
    fit_control = fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
  free <- fit_model()
  cell <- which(!is.na(free$map$mapList$rec_pars))[1]

  # Fix mean recruitment (log scale) away from its estimate, so a released one would move.
  inits <- free$estimated_params
  inits$rec_pars[cell] <- inits$rec_pars[cell] + 0.05
  map <- free$map
  map$mapList$rec_pars[cell] <- NA
  map$mapFactor$rec_pars <- factor(map$mapList$rec_pars)
  held <- fit_model(inits = inits, map = map)
  expect_gt(held$opt$objective, free$opt$objective + 1e-3)

  jit <- suppressWarnings(suppressMessages(jitter(held, njitter = 1, sd = 0.01, cores = 1)))
  expect_length(jit$Rceattle_list, 1L)
  refit <- jit$Rceattle_list[[1]]
  expect_equal(refit$estimated_params$rec_pars[cell], held$estimated_params$rec_pars[cell])
  expect_equal(length(refit$obj$par), length(held$obj$par))
  expect_equal(refit$opt$objective, held$opt$objective, tolerance = 1e-4)

  sims <- suppressWarnings(suppressMessages(self_test(held, nsim = 1, cores = 1, start = "estimated")))
  expect_length(sims, 1L)
  expect_equal(sims[[1]]$estimated_params$rec_pars[cell], held$estimated_params$rec_pars[cell])
  expect_equal(length(sims[[1]]$obj$par), length(held$obj$par))

  # Redrawing a fixed deviation is a contradiction the refit cannot show on its own.
  dev_map <- free$map
  dev_map$mapList$rec_dev[1, 2] <- NA
  dev_map$mapFactor$rec_dev <- factor(dev_map$mapList$rec_dev)
  held_dev <- fit_model(map = dev_map)
  expect_warning(
    suppressMessages(self_test(held_dev, nsim = 1, cores = 1, start = "estimated", process = "recruitment")),
    "redrew deviations the fit's map fixes, in: rec_dev")

  # TMB reads mapFactor, so a cell fixed only there must not be perturbed either.
  factor_only <- free$map
  factor_only$mapFactor$rec_pars[cell] <- NA
  held_factor <- fit_model(inits = inits, map = factor_only)
  jit <- suppressWarnings(suppressMessages(jitter(held_factor, njitter = 1, sd = 0.2, cores = 1)))
  expect_equal(jit$Rceattle_list[[1]]$estimated_params$rec_pars[cell],
               held_factor$estimated_params$rec_pars[cell])
})
