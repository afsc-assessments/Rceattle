# Fleets sharing a Selectivity_index share one parameter block, but the columns
# that shape the curve -- Selectivity, Selectivity_dimension,
# Age_first_selected, N_sel_bins, Sel_norm_bin, Sel_norm_bin_upper -- are read
# per fleet, so a group disagreeing on any of them does not share one curve.
# Time_varying_sel is the exception build_map() reconciles to the lead fleet's
# setting, leaving the curves matched. data_check() reports both, differently.

testthat::skip_on_cran()

# A fishery plus a mirror fleet built the documented way -- copy the fishery's
# whole row, change only identity and catchability -- then break one column.
.mirror_dat <- function(mutate = identity) {
  d   <- make_test_data(nyrs = 12, nages = 6, seed = 42)
  fc  <- d$fleet_control
  fsh <- which(fc$Fleet_type == "Fishery")[1]
  srv <- which(fc$Fleet_type == "Survey")[1]
  testthat::skip_if(is.na(fsh) || is.na(srv), "fixture lacks both fleet types")

  m <- fc[fsh, ]
  m$Fleet_name         <- "Fishery_CPUE"
  m$Fleet_code         <- 3L
  m$Fleet_type         <- "Survey"
  m$Catchability_index <- 3L
  m$Catchability       <- "Estimated"
  m$Catchability_init  <- 1
  m$Estimate_index_sd  <- 0
  d$fleet_control <- rbind(fc, mutate(m))

  add <- d$index_data[d$index_data$Fleet_code == fc$Fleet_code[srv], , drop = FALSE]
  add$Fleet_code <- 3L
  add$Fleet_name <- "Fishery_CPUE"
  d$index_data   <- rbind(d$index_data, add)
  d
}

# data_check()'s shared-Selectivity_index warnings only, as text.
.mirror_warnings <- function(d) {
  w <- character(0)
  # data_check() is internal: call it bare (the test env's parent is the package
  # namespace). Rceattle::data_check() errors, and try() would swallow that into
  # an empty warning list -- every assertion below would then pass for nothing.
  withCallingHandlers(
    suppressMessages(try(
      data_check(suppressMessages(switch_check(d))),
      silent = TRUE)),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    })
  grep("sharing Selectivity_index", w, value = TRUE)
}


testthat::test_that("a correctly copied mirror row raises nothing", {
  # The baseline the other cases are measured against: if this warned, the check
  # would be noise and users would learn to ignore it.
  testthat::expect_length(.mirror_warnings(.mirror_dat()), 0)
})


testthat::test_that("a curve-shaping column that differs within the group is reported", {
  # Each of these is read per fleet by the cpp, so a difference means the two
  # fleets get different curves despite sharing the block.
  breaks <- list(
    Selectivity_dimension = function(m) { m$Selectivity_dimension <- "Length"; m },
    Age_first_selected    = function(m) { m$Age_first_selected <- 3L; m },
    Sel_norm_bin          = function(m) { m$Sel_norm_bin <- 2L; m },
    Sel_norm_bin_upper    = function(m) { m$Sel_norm_bin_upper <- 4L; m }
  )
  for (col in names(breaks)) {
    w <- .mirror_warnings(.mirror_dat(breaks[[col]]))
    testthat::expect_length(w, 1)
    testthat::expect_match(w[1], col, fixed = TRUE, info = col)
    testthat::expect_match(w[1], "will not share one selectivity",
                           fixed = TRUE, info = col)
  }
})


testthat::test_that("a blank counts as a value among the shaping columns", {
  # Sel_norm_bin is blank on the fixture's fishery, so the case above is
  # blank-against-2, and it has to warn: blank means "do not normalize", a
  # different curve from normalizing at bin 2, not "inherit the lead's value".
  d <- .mirror_dat(function(m) { m$Sel_norm_bin <- 2L; m })
  # Blank canonicalizes to the word "Off", which is the same instruction.
  testthat::expect_equal(d$fleet_control$Sel_norm_bin[
    d$fleet_control$Fleet_name == "Fishery"][1], "Off")
  testthat::expect_length(.mirror_warnings(d), 1)
})


testthat::test_that("a blank first selected age agrees with naming minage", {
  # Blank means selected from minage, so the two fleets get the same curve.
  w <- .mirror_warnings(.mirror_dat(function(m) { m$Age_first_selected <- 1L; m }))
  testthat::expect_length(w, 0)
})


testthat::test_that("Time_varying_sel is reported as overridden, not as divergent", {
  # build_map() copies the lead's deviation map over the group, so the curves
  # still match and the wording has to differ.
  w <- .mirror_warnings(.mirror_dat(function(m) { m$Time_varying_sel <- "IID"; m }))
  testthat::expect_length(w, 1)
  testthat::expect_match(w[1], "Time_varying_sel", fixed = TRUE)
  testthat::expect_match(w[1], "the others are ignored", fixed = TRUE)
  testthat::expect_false(grepl("will not share one selectivity", w[1], fixed = TRUE))
})


testthat::test_that("an Off fleet in the group is not reported", {
  # An Off fleet never leads the group and its parameters are all mapped out, so
  # its columns cannot change what the estimated fleets get.
  d <- .mirror_dat(function(m) {
    m$Fleet_type            <- "Off"
    m$Selectivity_dimension <- "Length"
    m
  })
  testthat::expect_length(.mirror_warnings(d), 0)
})


testthat::test_that("the bundled datasets are clean", {
  # GOA2018SS is the worked mirroring example (fleets 9 and 10 share selectivity
  # and q), so a warning here would mean the check is wrong, not the data.
  for (nm in c("BS2017SS", "BS2017MS", "GOA2018SS")) {
    d <- get(nm, envir = asNamespace("Rceattle"))
    testthat::expect_length(.mirror_warnings(d), 0)
  }
})
