# Conditional age-at-length. CAAL is the age composition WITHIN a length bin, so
# it is predicted only for a fleet whose selectivity is length-dimensioned
# (selectivity.hpp writes sel_at_length for those alone). An age-dimensioned
# fleet predicts nothing and leaves caal_hat identically zero, so both the draw
# and the data_check() warning turn on Selectivity_dimension.
#
# Split out of test-functions-sim-mod-comp.R, which was 27% of the coverage
# suite in one file. Shared fixtures: helpers-comp-sim.R.
testthat::skip_on_cran()

testthat::test_that("CAAL is drawn from its own predicted composition", {
  testthat::skip_if_not_installed("TMB")

  # CAAL is only PREDICTED for a fleet whose selectivity is length-dimensioned:
  # pred_CAAL is filled under `flt_sel_dim(flt) == 1` (ceattle.cpp, section
  # 10.2), so an age-dimensioned fleet leaves caal_hat identically zero and there
  # is nothing to draw. Every fixture in this repo puts its CAAL rows on an
  # age-based fleet, which is why this path went unexercised; set the dimension
  # and it comes alive.
  d <- make_test_data(nyrs = 20, nages = 6, seed = 123, growth = "vonBertalanffy")
  testthat::expect_gt(nrow(d$caal_data), 0)
  d$fleet_control$Selectivity_dimension[d$fleet_control$Fleet_code == 1] <- "Length"
  d$caal_data$Sample_size <- 100

  # At the starting values, not a fit: these flat data leave growth free to collapse
  # P(L|age) onto one length bin, and a draw needs a prediction in every row.
  fit <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    d, file = NULL, estimateMode = "DebugBuild", msmMode = 0, random_rec = FALSE,
    growthFun = Rceattle::build_growth(fun = "vonBertalanffy"),
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0,
                                        phase = FALSE))))
  # The fixture guard: without a prediction the assertions below are vacuous.
  testthat::expect_equal(sum(rowSums(fit$quantities$caal_hat) > 0),
                         nrow(fit$quantities$caal_hat))

  cols <- grep("^CAAL_", names(fit$data_list$caal_data))
  set.seed(6)
  sim <- suppressWarnings(Rceattle::sim_mod(fit, simulate = TRUE))
  got <- as.matrix(sim$caal_data[, cols])

  testthat::expect_false(anyNA(got))
  testthat::expect_true(all(got >= 0))
  testthat::expect_equal(got, round(got))
  testthat::expect_equal(unname(rowSums(got)),
                         unname(fit$data_list$caal_data$Sample_size))
  # It sampled rather than returning the expectation.
  exp_dat <- suppressWarnings(Rceattle::sim_mod(fit, simulate = FALSE))
  testthat::expect_false(isTRUE(all.equal(
    got, unname(as.matrix(exp_dat$caal_data[, cols])))))

  # And the draw is centred on the predicted composition.
  N <- fit$data_list$caal_data$Sample_size[1]
  set.seed(7)
  reps <- replicate(400, suppressWarnings(
    as.numeric(Rceattle::sim_mod(fit, simulate = TRUE)$caal_data[1, cols])))
  p <- fit$quantities$caal_hat[1, seq_along(cols)]
  p <- p / sum(p)
  keep <- which(N * p > 5)
  testthat::expect_gt(length(keep), 1)
  z <- (rowMeans(reps)[keep] - N * p[keep]) /
    (sqrt(N * p[keep] * (1 - p[keep])) / sqrt(ncol(reps)))
  testthat::expect_lt(max(abs(z)), 4)
})


testthat::test_that("data_check() warns about CAAL data on an age-dimensioned fleet", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # CAAL is the age composition WITHIN a length bin, so it is predicted from
  # selectivity-at-length convolved with the growth matrix. selectivity.hpp
  # writes sel_at_length only for a length-dimensioned fleet, so an age-based one
  # predicts nothing -- and the likelihood still runs, scoring the observations
  # against a flat composition it cannot move. Selectivity_dimension defaults to
  # "Age", so this was the default outcome.
  d <- make_test_data(nyrs = 20, nages = 6, seed = 123, growth = "vonBertalanffy")
  testthat::expect_gt(nrow(d$caal_data), 0)
  d$caal_data$Sample_size <- 100
  d$HCR <- 0
  prep <- function(x) suppressMessages(suppressWarnings(
    Rceattle::switch_check(Rceattle::clean_data(x))))
  quiet_check <- function(x) suppressMessages(suppressWarnings(data_check(x)))

  # The fixture's CAAL fleet is age-dimensioned, so it warns. A warning rather
  # than an error only because such models currently fit -- the CAAL data in
  # them inform nothing, so this should tighten once configurations are checked.
  warn_of <- function(x) {
    w <- character(0)
    withCallingHandlers(suppressMessages(data_check(x)),
      warning = function(e) { w <<- c(w, conditionMessage(e))
                              invokeRestart("muffleWarning") })
    w
  }
  testthat::expect_true(any(grepl("Selectivity_dimension", warn_of(prep(d)))))

  # Length-dimensioned: silent on this point.
  ok <- d
  ok$fleet_control$Selectivity_dimension[ok$fleet_control$Fleet_code == 1] <- "Length"
  testthat::expect_false(any(grepl("Selectivity_dimension", warn_of(prep(ok)))))

  # And inactive CAAL rows are not the model's problem.
  off <- d
  off$caal_data$Sample_size <- 0
  testthat::expect_false(any(grepl("Selectivity_dimension", warn_of(prep(off)))))
})
