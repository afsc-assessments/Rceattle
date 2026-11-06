# data_check()'s warnings describe the data_list, not the fit, so a diagnostic
# that refits the same data repeats them once per peel, jitter, or MSE
# iteration. `fit_mod(quiet_data_check = TRUE)` drops them; `.refit_like()` sets
# it, which covers all eight refit entry points. Errors are unaffected -- a
# data_list that cannot be fitted still stops -- and so are convergence and TMB
# warnings.

testthat::skip_on_cran()

# A fleet_control that raises a data_check WARNING and nothing else: two fleets
# in one Selectivity_index group disagreeing on a curve-shaping column.
.warned_dat <- function() {
  d   <- make_test_data(nyrs = 12, nages = 6, seed = 42)
  fc  <- d$fleet_control
  fsh <- which(fc$Fleet_type == "Fishery")[1]
  srv <- which(fc$Fleet_type == "Survey")[1]
  testthat::skip_if(is.na(fsh) || is.na(srv), "fixture lacks both fleet types")

  m <- fc[fsh, ]
  m$Fleet_name          <- "Fishery_CPUE"
  m$Fleet_code          <- 3L
  m$Fleet_type          <- "Survey"
  m$Catchability_index  <- 3L
  m$Catchability        <- "Estimated"
  m$Catchability_init   <- 1
  m$Estimate_index_sd   <- 0
  m$Age_first_selected  <- 3L          # differs from the fishery it mirrors
  d$fleet_control <- rbind(fc, m)

  add <- d$index_data[d$index_data$Fleet_code == fc$Fleet_code[srv], , drop = FALSE]
  add$Fleet_code <- 3L
  add$Fleet_name <- "Fishery_CPUE"
  d$index_data   <- rbind(d$index_data, add)
  d
}

.warnings_of <- function(expr) {
  w <- character(0)
  withCallingHandlers(suppressMessages(force(expr)),
                      warning = function(cnd) {
                        w <<- c(w, conditionMessage(cnd))
                        invokeRestart("muffleWarning")
                      })
  w
}


testthat::test_that("an ordinary fit still reports data_check warnings", {
  testthat::skip_if_not_installed("TMB")
  # The baseline: silence on a refit is only useful if the user heard it once.
  w <- .warnings_of(Rceattle::fit_mod(
    .warned_dat(), file = NULL, estimateMode = 3, msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0, phase = FALSE)))
  testthat::expect_true(any(grepl("sharing Selectivity_index", w)))
})


testthat::test_that("quiet_data_check drops them", {
  testthat::skip_if_not_installed("TMB")
  w <- .warnings_of(Rceattle::fit_mod(
    .warned_dat(), file = NULL, estimateMode = 3, msmMode = 0, random_rec = FALSE,
    quiet_data_check = TRUE,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0, phase = FALSE)))
  testthat::expect_false(any(grepl("sharing Selectivity_index", w)))
})


testthat::test_that("a refit through .refit_like() is quiet", {
  testthat::skip_if_not_installed("TMB")
  # This is what makes it reach retrospective(), jitter(), run_mse() and the
  # other five entry points without each of them opting in.
  fit <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    .warned_dat(), file = NULL, estimateMode = 3, msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0, phase = FALSE))))

  w <- .warnings_of(.refit_like(
    data_list    = fit$data_list,
    inits        = fit$estimated_params,
    estimateMode = 3))
  testthat::expect_false(any(grepl("sharing Selectivity_index", w)))
})


testthat::test_that("quiet_data_check does not swallow errors", {
  testthat::skip_if_not_installed("TMB")
  # Warnings only. A data_list data_check() rejects must still stop the fit --
  # otherwise the refits would run on input the model cannot use.
  d <- make_test_data(nyrs = 12, nages = 6, seed = 42)
  d$fleet_control$Fleet_code[1] <- 99L        # must equal the row number
  testthat::expect_error(
    suppressMessages(suppressWarnings(Rceattle::fit_mod(
      d, file = NULL, estimateMode = 3, msmMode = 0, random_rec = FALSE,
      quiet_data_check = TRUE,
      fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0,
                                          phase = FALSE)))))
})
