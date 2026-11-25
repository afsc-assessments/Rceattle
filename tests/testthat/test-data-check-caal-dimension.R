# CAAL is the age composition WITHIN a length bin, so the model predicts it as
# selectivity-at-length convolved with the growth matrix (ceattle.cpp, section
# 10.2). selectivity.hpp writes sel_at_length only for a length-dimensioned
# fleet and leaves it zero otherwise, so pred_CAAL on an age-dimensioned fleet
# is zero however the surrounding gates are written -- there is no
# selectivity-at-length to convolve with.
#
# The likelihood still runs. With pred_CAAL zero, comp_offset makes the
# predicted composition uniform, so the observations are scored against a flat
# composition no parameter can move: they add a term to the objective and inform
# nothing. Selectivity_dimension defaults to "Age", so this is the default
# outcome for anyone adding CAAL data rather than an unusual mistake.
testthat::skip_on_cran()

.caal_dim_prep <- function(x) {
  x$HCR <- 0                       # build_hcr() would set this
  suppressMessages(suppressWarnings(
    Rceattle::switch_check(Rceattle::clean_data(x))))
}

# Warnings raised by data_check(), as a character vector.
.caal_dim_warnings <- function(x) {
  w <- character(0)
  withCallingHandlers(suppressMessages(data_check(x)),
                      warning = function(e) {
                        w <<- c(w, conditionMessage(e))
                        invokeRestart("muffleWarning")
                      })
  w
}


testthat::test_that("data_check() flags CAAL data on an age-dimensioned fleet", {
  testthat::skip_if_not_installed("TMB")

  # The fixture carries CAAL only under parametric growth, and puts it on a
  # fleet whose Selectivity_dimension is the default "Age".
  d <- make_test_data(nyrs = 20, nages = 6, seed = 123, growth = "vonBertalanffy")
  testthat::expect_gt(nrow(d$caal_data), 0)
  d$caal_data$Sample_size <- 100
  testthat::expect_identical(
    as.character(d$fleet_control$Selectivity_dimension[d$fleet_control$Fleet_code == 1]),
    "Age")

  testthat::expect_true(
    any(grepl("Selectivity_dimension", .caal_dim_warnings(.caal_dim_prep(d)))))
})


testthat::test_that("a length-dimensioned fleet with CAAL data does not flag", {
  testthat::skip_if_not_installed("TMB")

  d <- make_test_data(nyrs = 20, nages = 6, seed = 123, growth = "vonBertalanffy")
  d$caal_data$Sample_size <- 100
  d$fleet_control$Selectivity_dimension[d$fleet_control$Fleet_code == 1] <- "Length"

  testthat::expect_false(
    any(grepl("Selectivity_dimension", .caal_dim_warnings(.caal_dim_prep(d)))))
})


testthat::test_that("CAAL rows that are not fitted do not flag", {
  testthat::skip_if_not_installed("TMB")

  # A row with no sample size, or a prediction-only year, is not in the
  # likelihood, so its selectivity dimension is nobody's problem.
  d <- make_test_data(nyrs = 20, nages = 6, seed = 123, growth = "vonBertalanffy")

  off <- d; off$caal_data$Sample_size <- 0
  testthat::expect_false(
    any(grepl("Selectivity_dimension", .caal_dim_warnings(.caal_dim_prep(off)))))

  proj <- d
  proj$caal_data$Sample_size <- 100
  proj$caal_data$Year <- -abs(proj$caal_data$Year)
  testthat::expect_false(
    any(grepl("Selectivity_dimension", .caal_dim_warnings(.caal_dim_prep(proj)))))
})


testthat::test_that("CAAL on a Fixed fleet labelled Length warns", {
  testthat::skip_if_not_installed("TMB")

  # Fixed selectivity is always at age, whatever Selectivity_dimension says, so
  # CAAL on such a fleet is unpredicted even though its label reads "Length".
  d <- make_test_data(nyrs = 20, nages = 6, seed = 123, growth = "vonBertalanffy")
  d$caal_data$Sample_size <- 100
  d$fleet_control$Selectivity_dimension <- "Length"
  testthat::expect_false(any(grepl("Fixed", .caal_dim_warnings(.caal_dim_prep(d)))))

  d$fleet_control$Selectivity[d$fleet_control$Fleet_code == 1] <- "Fixed"
  d$emp_sel <- cbind(data.frame(Fleet_name = "Survey", Fleet_code = 1L, Species = 1L,
                                Sex = 0L, Year = d$styr:d$endyr),
                     as.data.frame(matrix(1, length(d$styr:d$endyr), 6,
                                          dimnames = list(NULL, paste0("Comp_", 1:6)))))
  w <- .caal_dim_warnings(.caal_dim_prep(d))
  testthat::expect_true(any(grepl("whose Selectivity is 'Fixed'", w)))
  testthat::expect_false(any(grepl("Selectivity_dimension is not 'Length'", w)))
})
