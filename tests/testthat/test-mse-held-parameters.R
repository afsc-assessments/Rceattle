# run_mse() builds a new map at each assessment, because endyr and the array sizes
# move. Parameters an author held in the EM's map have to stay held in that new map.
#
# Provenance: the 2026 GOA Pacific cod EM holds 31 linkage coefficients and 3 sel_dn6
# cells its data cannot inform. Every refit released them; a 2-replicate pilot with
# random_rec on had 5 of 24 assessments at status FAIL and none at OK.

.free_map <- function(n = 6L) {
  v <- seq_len(n)
  list(mapList = list(b = v), mapFactor = list(b = factor(v)))
}

.cells <- function(positions, n_cells) structure(positions, n_cells = n_cells)
.hold <- function(map, holds, params = list(b = as.numeric(map$mapList$b))) {
  .apply_map_holds(map, params, holds)
}

test_that(".apply_map_holds() holds the cells and renumbers nothing else", {
  map <- .free_map()
  out <- .hold(map, list(b = .cells(c(2L, 5L), 6L)))
  expect_equal(out$map$mapList$b, c(1L, NA, 3L, 4L, NA, 6L))
  expect_equal(levels(out$map$mapFactor$b), c("1", "3", "4", "6"))
  expect_true(all(is.na(out$map$mapFactor$b[c(2L, 5L)])))
})

test_that(".apply_map_holds() leaves a map alone when nothing is held", {
  map <- .free_map()
  expect_identical(.hold(map, NULL)$map, map)
  expect_identical(.hold(map, list())$map, map)
})

test_that("a held cell keeps its position when the array has grown", {
  # A later assessment has more years, added at the end of the last dimension.
  earlier <- matrix(seq_len(6L), 2L, 3L)
  earlier[1, 2] <- NA
  grown <- matrix(seq_len(8L), 2L, 4L)
  map <- list(mapList = list(b = grown), mapFactor = list(b = factor(grown)))
  out <- .hold(map, list(b = .cells(which(is.na(earlier)), length(earlier))))
  expect_true(is.na(out$map$mapList$b[1, 2]))
  expect_equal(sum(is.na(out$map$mapList$b)), 1L)
  expect_equal(dim(out$map$mapList$b), c(2L, 4L))
})

test_that("a time block held for every earlier year stays held in its new years", {
  # Two blocks over six years; two new years join the second block.
  block <- c(1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L)
  map <- list(mapList = list(b = block), mapFactor = list(b = factor(block)))
  out <- .hold(map, list(b = .cells(4:6, 6L)))
  expect_equal(out$map$mapList$b, c(1L, 1L, 1L, NA, NA, NA, NA, NA))
  # run_mse() starts a new year from the last year, here 0; it takes the block's value.
  values <- list(b = c(0.1, 0.1, 0.1, 0.5, 0.5, 0.5, 0, 0))
  expect_equal(.hold(map, list(b = .cells(4:6, 6L)), values)$params$b,
               c(0.1, 0.1, 0.1, 0.5, 0.5, 0.5, 0.5, 0.5))

  # A block with an earlier year still free is still estimated, new years included.
  out <- .hold(map, list(b = .cells(4L, 6L)))
  expect_equal(out$map$mapList$b, c(1L, 1L, 1L, NA, 2L, 2L, 2L, 2L))

  # Yearly deviations share no level, so a new year stays free.
  yearly <- seq_len(8L)
  map <- list(mapList = list(b = yearly), mapFactor = list(b = factor(yearly)))
  out <- .hold(map, list(b = .cells(1:6, 6L)))
  expect_equal(out$map$mapList$b, c(rep(NA, 6L), 7L, 8L))
})

test_that(".map_holds() reads a Laplace hyperparameter and a cell fixed only in mapFactor", {
  d <- make_test_data()
  d$fleet_control$Time_varying_sel[1] <- "IID"
  d$fleet_control$Sel_dev_sd[1] <- 0.2
  fit <- suppressWarnings(suppressMessages(fit_mod(
    data_list = d, estimateMode = 3, random_rec = TRUE, random_sel = TRUE,
    fit_control = fit_control(verbose = 0, phase = FALSE, getsd = FALSE))))
  expect_length(.map_holds(fit), 0L)
  fit$map$mapList$R_log_sd[1] <- NA
  fit$map$mapFactor$R_log_sd <- factor(fit$map$mapList$R_log_sd)
  cell <- which(!is.na(fit$map$mapFactor$sel_dev_log_sd))[1]
  fit$map$mapFactor$sel_dev_log_sd[cell] <- NA
  holds <- .map_holds(fit)
  expect_setequal(names(holds), c("R_log_sd", "sel_dev_log_sd"))
  expect_equal(as.vector(holds$R_log_sd), 1L)
  expect_equal(as.vector(holds$sel_dev_log_sd), unname(cell))
})

test_that(".refit_like(map_holds =) reaches the map a build makes", {
  fit <- suppressWarnings(suppressMessages(fit_mod(
    data_list = make_test_data(), estimateMode = 3,
    fit_control = fit_control(verbose = 0, phase = FALSE, getsd = FALSE))))
  free <- unname(which(!is.na(fit$map$mapFactor$rec_dev))[1:2])
  out <- suppressWarnings(suppressMessages(.refit_like(
    data_list = fit$data_list, inits = fit$estimated_params, estimateMode = 3,
    map_holds = list(rec_dev = .cells(free, length(fit$map$mapFactor$rec_dev))))))
  expect_true(all(is.na(out$map$mapFactor$rec_dev[free])))
  expect_equal(length(out$obj$par), length(fit$obj$par) - 2L)
})

test_that("fit_mod() starts a held time block's new year at the block's value", {
  d <- make_test_data()
  flt <- d$index_data$Fleet_code[1]
  d$fleet_control$Time_varying_sel[flt] <- "Block"
  fit <- suppressWarnings(suppressMessages(fit_mod(
    data_list = d, estimateMode = 3,
    fit_control = fit_control(verbose = 0, phase = FALSE, getsd = FALSE))))
  cells <- which(!is.na(fit$map$mapList$sel_inf_dev))
  last <- max(cells)
  inits <- fit$estimated_params
  inits$sel_inf_dev[cells] <- 0.3
  inits$sel_inf_dev[last] <- 0
  out <- suppressWarnings(suppressMessages(.refit_like(
    data_list = fit$data_list, inits = inits, estimateMode = 3,
    map_holds = list(sel_inf_dev = .cells(unname(setdiff(cells, last)), last - 1L)))))
  expect_true(is.na(out$map$mapFactor$sel_inf_dev[last]))
  expect_equal(out$estimated_params$sel_inf_dev[last], 0.3)
})

test_that(".map_holds() finds only the cells a map fixes beyond build_map()", {
  fit <- suppressWarnings(suppressMessages(fit_mod(
    data_list = make_test_data(), estimateMode = 3,
    fit_control = fit_control(verbose = 0, phase = FALSE, getsd = FALSE))))
  expect_length(.map_holds(fit), 0L)

  free <- which(!is.na(fit$map$mapList$rec_dev))[1:2]
  fit$map$mapList$rec_dev[free] <- NA
  fit$map$mapFactor$rec_dev <- factor(fit$map$mapList$rec_dev)
  expect_equal(names(.map_holds(fit)), "rec_dev")
  expect_equal(as.vector(.map_holds(fit)$rec_dev), unname(free))

  # A projection-only or debug fit's stored map is not read as holds.
  fit$data_list$estimateMode <- 2
  expect_length(.map_holds(fit), 0L)
  fit$data_list$estimateMode <- 4
  expect_length(.map_holds(fit), 0L)
  fit$data_list$estimateMode <- 3

  # An EM that run_mse() stored has had its map removed.
  fit$map <- NULL
  expect_length(.map_holds(fit), 0L)
})

test_that("run_mse() keeps a held parameter at the EM's value at every assessment", {
  testthat::skip_on_cran()
  data(BS2017SS)
  d <- BS2017SS
  d$projyr <- d$endyr + 4L
  d$fleet_control$Proj_F_proportion <- rep(1, nrow(d$fleet_control))
  hcr <- build_hcr(HCR = 5, Ftarget = 0.4, Flimit = 0.35, Plimit = 0.2, Alpha = 0.05)
  ctl <- fit_control(phase = TRUE, getsd = FALSE, verbose = 0)
  om <- suppressWarnings(suppressMessages(fit_mod(d, estimateMode = 1, msmMode = 0, fit_control = ctl)))

  # Hold one catchability that the fit estimates, slightly away from its estimate, so
  # that a released parameter would move.
  flt <- which(!is.na(om$map$mapList$index_log_q))[1]
  inits <- om$estimated_params
  inits$index_log_q[flt] <- inits$index_log_q[flt] + 0.05
  mp <- om$map
  mp$mapList$index_log_q[flt] <- NA
  mp$mapFactor$index_log_q <- factor(mp$mapList$index_log_q)
  # And one fishing mortality, whose array grows as endyr advances.
  f_cell <- which(!is.na(mp$mapList$log_F))[1]
  mp$mapList$log_F[f_cell] <- NA
  mp$mapFactor$log_F <- factor(mp$mapList$log_F)
  em <- suppressWarnings(suppressMessages(fit_mod(
    d, inits = inits, map = mp, estimateMode = 0, msmMode = 0, HCR = hcr,
    fit_control = fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
  held <- em$estimated_params$index_log_q[flt]
  expect_equal(unname(held), unname(inits$index_log_q[flt]))

  expect_message(
    mse <- suppressWarnings(run_mse(om = om, em = em, nsim = 1, assessment_period = 2, sampling_period = 1,
                                    simulate_data = FALSE, sample_rec = FALSE, cores = 1, seed = 1)),
    "index_log_q \\(1\\)")
  fits <- mse$Sim_1$EM[grepl("EM_yr_", names(mse$Sim_1$EM))]
  expect_gt(length(fits), 0L)
  for (m in fits) {
    expect_equal(unname(m$estimated_params$index_log_q[flt]), unname(held))
    expect_equal(m$estimated_params$log_F[f_cell], em$estimated_params$log_F[f_cell])
  }
  # The refit did run: early-year SSB moved, which a new log_F column alone would not show.
  expect_gt(ncol(fits[[1]]$estimated_params$log_F), ncol(em$estimated_params$log_F))
  expect_false(isTRUE(all.equal(fits[[1]]$quantities$ssb[1, 1:10], em$quantities$ssb[1, 1:10], tolerance = 1e-10)))

  # The refit to regenerated past data holds it too.
  mse <- suppressWarnings(suppressMessages(run_mse(
    om = om, em = em, nsim = 1, assessment_period = 2, sampling_period = 1, simulate_data = FALSE,
    regenerate_past = TRUE, sample_rec = FALSE, cores = 1, seed = 1)))
  expect_equal(unname(mse$Sim_1$EM[[1]]$estimated_params$index_log_q[flt]), unname(held))
})
