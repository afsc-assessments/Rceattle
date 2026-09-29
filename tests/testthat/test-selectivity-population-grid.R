# A length-based selectivity curve is built on the population length bins, so
# everything that walks the curve has to walk all of them. The normalizer, the
# below-first-bin zeroing and the projection copy were bounded by the DATA bin
# count instead: on a finer population grid the projection kept only the lowest
# lengths, where an ascending curve is ~0, so projected selectivity-at-age came
# back identically zero and with it projected F, catch and the reference points.
# Measured on the fixture below at 20 data bins against 96 population bins.

testthat::skip_on_cran()

e <- new.env(parent = asNamespace("Rceattle"))
for (h in c("helpers.R", "helpers-make-msm-data.R")) {
  sys.source(testthat::test_path(h), envir = e)
}

set.seed(11)
pg_data <- local({
  d <- e$make_msm_test_data()$data_list
  flt <- which(d$fleet_control$Species == 1 & d$fleet_control$Fleet_type == 1)[1]
  d$fleet_control$Selectivity_dimension[flt] <- "Length"
  d$fleet_control$Selectivity[flt] <- "Logistic"
  d$fleet_control$Time_varying_sel[flt] <- "Off"
  list(d = d, flt = flt)
})
pg_edges <- sort(unique(pg_data$d$caal_data$Length[pg_data$d$caal_data$Species == 1]))
pg_fine  <- sort(unique(c(pg_edges, seq(min(pg_edges), max(pg_edges), by = 1))))

pg_fit <- function(d, ...) suppressWarnings(suppressMessages(Rceattle::fit_mod(
  data_list = d, inits = NULL, estimateMode = 3, msmMode = 0, random_rec = FALSE,
  growthFun = Rceattle::build_growth(fun = "vonBertalanffy", ...),
  fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0, getsd = FALSE))))

testthat::test_that("a finer population grid fills every projection bin", {
  testthat::skip_if_not_installed("TMB")
  d <- pg_data$d; flt <- pg_data$flt
  m <- pg_fit(d, pop_lengths = list(pg_fine, pg_fine))

  sl <- m$quantities$sel_at_length
  testthat::expect_equal(dim(sl)[3], length(pg_fine))   # the premise: the finer grid
  nh <- d$endyr - d$styr + 1

  # Projection is the terminal hindcast year, on every bin, not just the first
  # length(pg_edges) of them.
  for (yr in c(nh + 1, dim(sl)[4])) {
    testthat::expect_equal(unname(sl[flt, 1, , yr]), unname(sl[flt, 1, , nh]))
  }
  testthat::expect_gt(max(sl[flt, 1, , nh + 1]), 0.99)   # normalized, so it reaches 1
  testthat::expect_equal(unname(m$quantities$sel_at_age[flt, 1, , nh + 1]),
                         unname(m$quantities$sel_at_age[flt, 1, , nh]))
})

testthat::test_that("the bin columns name DATA bins on either grid", {
  testthat::skip_if_not_installed("TMB")
  d <- pg_data$d; flt <- pg_data$flt
  # Bin_first_selected is a 1-based data-bin ordinal, so bins below the 5th data
  # bin's lower edge are zeroed whether or not a population grid is supplied.
  d$fleet_control$Bin_first_selected[flt] <- 5L
  cut_at <- pg_edges[5]

  coarse <- pg_fit(d)
  fine   <- pg_fit(d, pop_lengths = list(pg_fine, pg_fine))
  nh <- d$endyr - d$styr + 1

  zero_below <- function(m, edges) {
    s <- m$quantities$sel_at_length[flt, 1, seq_along(edges), nh]
    c(below = max(s[edges < cut_at]), at_or_above = max(s[edges >= cut_at]))
  }
  testthat::expect_equal(unname(zero_below(coarse, pg_edges)["below"]), 0)
  testthat::expect_equal(unname(zero_below(fine, pg_fine)["below"]), 0)
  testthat::expect_gt(unname(zero_below(fine, pg_fine)["at_or_above"]), 0)
})
