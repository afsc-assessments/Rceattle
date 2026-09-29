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

# Every column that names a selectivity bin is a DATA bin ordinal, so on a finer
# population grid it cannot address the grid the curve is built on. Translating
# one column and not the rest is worse than refusing: Bin_first_selected on the
# population grid against N_sel_bins on the data grid gives the non-parametric
# base curve a negative length, which aborts R inside MakeADFun.
testthat::test_that("a bin-naming column on a finer grid is refused", {
  d <- pg_data$d; flt <- pg_data$flt
  fine <- list(pg_fine, pg_fine)

  refused <- function(dd) {
    dd$pop_lengths <- fine
    inherits(tryCatch(suppressMessages(suppressWarnings(
      Rceattle:::data_check(Rceattle::switch_check(dd)))),
      error = function(e) e), "error")
  }
  # A parametric curve naming no bin is a function of length: allowed, and it is
  # what the population grid exists for.
  testthat::expect_false(refused(d))

  for (cl in c("Bin_first_selected", "Sel_norm_bin", "Sel_norm_bin_upper")) {
    dd <- d; dd$fleet_control[[cl]][flt] <- 5L
    testthat::expect_true(refused(dd), info = cl)
  }
  # ... and a bin-indexed form, whatever its columns say.
  dd <- d
  dd$fleet_control$Selectivity[flt] <- "NonParametric"
  dd$fleet_control$Sel_curve_pen1[flt] <- 10
  dd$fleet_control$Sel_curve_pen2[flt] <- 10
  testthat::expect_true(refused(dd))

  # The same fleets are fine when the two grids coincide.
  dd <- d; dd$fleet_control$Sel_norm_bin[flt] <- 5L
  testthat::expect_false(inherits(tryCatch(suppressMessages(suppressWarnings(
    Rceattle:::data_check(Rceattle::switch_check(dd)))),
    error = function(e) e), "error"))
})

# Bin_first_selected carries no schema default, so a workbook can arrive without
# it. The refusal above read it first and every flag inherited that length, so
# one absent column silently excused every fleet -- including the bin-indexed
# forms, which the column has nothing to do with.
testthat::test_that("the refusal survives a fleet_control with no Bin_first_selected", {
  d <- pg_data$d; flt <- pg_data$flt
  d$fleet_control$Selectivity[flt] <- "NonParametric"
  d$fleet_control$Sel_curve_pen1[flt] <- 10
  d$fleet_control$Sel_curve_pen2[flt] <- 10
  d$fleet_control$Bin_first_selected <- NULL
  d$pop_lengths <- list(pg_fine, pg_fine)

  sc <- suppressMessages(suppressWarnings(Rceattle::switch_check(d)))
  # The premise: nothing fills the column back in.
  testthat::expect_false("Bin_first_selected" %in% colnames(sc$fleet_control))

  testthat::expect_error(
    suppressMessages(suppressWarnings(Rceattle:::data_check(sc))),
    "indexed by bin")
})
