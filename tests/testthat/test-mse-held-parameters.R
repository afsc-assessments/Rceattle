# run_mse() builds a new map at each assessment, because endyr and the array sizes
# move. Parameters an author held in the EM's map have to stay held in that new map.
#
# Provenance: the 2026 GOA Pacific cod EM holds 31 linkage coefficients and 3 sel_dn6
# cells its data cannot inform. Every refit released them; a 2-replicate pilot had 5 of
# 24 assessments at status FAIL and none at OK.

.held_map <- function(cells, length_out = 6L) {
  v <- seq_len(length_out)
  v[cells] <- NA
  list(mapList = list(b = v), mapFactor = list(b = factor(v)))
}

test_that(".apply_map_holds() holds the cells and renumbers nothing else", {
  map <- .held_map(integer(0))
  out <- .apply_map_holds(map, list(b = c(2L, 5L)))
  expect_equal(out$mapList$b, c(1L, NA, 3L, 4L, NA, 6L))
  expect_equal(levels(out$mapFactor$b), c("1", "3", "4", "6"))
  expect_true(all(is.na(out$mapFactor$b[c(2L, 5L)])))
})

test_that(".apply_map_holds() leaves a map alone when nothing is held", {
  map <- .held_map(integer(0))
  expect_identical(.apply_map_holds(map, NULL), map)
  expect_identical(.apply_map_holds(map, list()), map)
})

test_that("a held cell keeps its position when the array has grown", {
  # A later assessment has more years, added at the end of the last dimension.
  grown <- matrix(seq_len(8L), 2L, 4L)
  map <- list(mapList = list(b = grown), mapFactor = list(b = factor(grown)))
  out <- .apply_map_holds(map, list(b = 3L))
  expect_true(is.na(out$mapList$b[1, 2]))
  expect_equal(sum(is.na(out$mapList$b)), 1L)
  expect_equal(dim(out$mapList$b), c(2L, 4L))
})

test_that(".map_holds() finds only what the author held", {
  fit <- suppressWarnings(suppressMessages(fit_mod(
    data_list = make_test_data(), estimateMode = 3,
    fit_control = fit_control(verbose = 0, phase = FALSE, getsd = FALSE))))
  expect_length(.map_holds(fit), 0L)

  free <- which(!is.na(fit$map$mapList$rec_dev))[1:2]
  fit$map$mapList$rec_dev[free] <- NA
  fit$map$mapFactor$rec_dev <- factor(fit$map$mapList$rec_dev)
  expect_equal(.map_holds(fit), list(rec_dev = free))

  # A projection-only fit's stored map is all NA, which is not a hold.
  fit$data_list$estimateMode <- 2
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
  em <- suppressWarnings(suppressMessages(fit_mod(
    d, inits = inits, map = mp, estimateMode = 0, msmMode = 0, HCR = hcr,
    fit_control = fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
  held <- em$estimated_params$index_log_q[flt]
  expect_equal(unname(held), unname(inits$index_log_q[flt]))

  expect_message(
    mse <- suppressWarnings(run_mse(om = om, em = em, nsim = 1, assessment_period = 2, sampling_period = 1,
                                    simulate_data = FALSE, sample_rec = FALSE, cores = 1, seed = 1)),
    "holding 1 parameter")
  fits <- mse$Sim_1$EM[grepl("EM_yr_", names(mse$Sim_1$EM))]
  expect_gt(length(fits), 0L)
  for (m in fits) expect_equal(unname(m$estimated_params$index_log_q[flt]), unname(held))
  # The refit did run: it estimated the new years' fishing mortality.
  expect_gt(ncol(fits[[1]]$estimated_params$log_F), ncol(em$estimated_params$log_F))
  expect_false(isTRUE(all.equal(fits[[1]]$quantities$ssb[1, 1:10], em$quantities$ssb[1, 1:10], tolerance = 1e-10)))

  # The refit to regenerated past data holds it too.
  mse <- suppressWarnings(suppressMessages(run_mse(
    om = om, em = em, nsim = 1, assessment_period = 2, sampling_period = 1, simulate_data = FALSE,
    regenerate_past = TRUE, sample_rec = FALSE, cores = 1, seed = 1)))
  expect_equal(unname(mse$Sim_1$EM[[1]]$estimated_params$index_log_q[flt]), unname(held))
})
