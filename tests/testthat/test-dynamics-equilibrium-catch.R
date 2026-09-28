# The initial equilibrium catch is a catch_data row at Year == styr - 1: the catch
# the stock yielded under the initial fishing mortality before the hindcast began.
#
# That year is NOT a free marker. Real data carries ordinary catch history there --
# GOA2018SS has 23 catch rows before styr, two of them on styr - 1 (1976) -- and an
# earlier version of this feature read those two as equilibrium observations. Under
# the default initMode the predicted equilibrium catch is 0, so the likelihood took
# log(0) and every GOA fit returned a non-finite objective, which surfaced as
# `optimHess: non-finite value supplied by optim` in the golden references.
#
# A negative sentinel cannot be used instead: run_mse() reserves negative Year for
# rows it splices back in as the next assessment's data, and its window filters are
# on abs(Year), so a -999 row survives them as year 999.
#
# The rule is therefore gated, not marked: a styr - 1 row is read as an equilibrium
# catch only under an initMode that estimates Finit (3, 4, 6), it is named in a
# message whenever it IS read, and it is otherwise dropped as history exactly as it
# was before the feature existed.

testthat::skip_on_cran()

fit3 <- function(d, im) suppressMessages(suppressWarnings(Rceattle::fit_mod(
  data_list = d, file = NULL, inits = NULL, estimateMode = 3, random_rec = FALSE,
  msmMode = 0, initMode = im,
  fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))

testthat::test_that("catch history at styr - 1 is not read as an equilibrium catch", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not(exists("GOA2018SS"))
  d <- Rceattle::GOA2018SS
  # The dataset really does carry the collision this guards.
  testthat::expect_gt(sum(d$catch_data$Year == d$styr - 1L, na.rm = TRUE), 0)

  f <- fit3(d, 2)                       # NonEquilibrium: Finit is held at 0
  testthat::expect_true(is.finite(f$quantities$jnll))
  testthat::expect_length(f$quantities$equil_catch_hat, 0)
})

testthat::test_that("a Finit-estimating mode does read it, and says so", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not(exists("GOA2018SS"))
  d <- Rceattle::GOA2018SS
  n <- sum(d$catch_data$Year == d$styr - 1L, na.rm = TRUE)

  testthat::expect_message(
    suppressWarnings(Rceattle::fit_mod(
      data_list = d, file = NULL, inits = NULL, estimateMode = 3,
      random_rec = FALSE, msmMode = 0, initMode = 3,
      fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))),
    "Initial equilibrium catch read from catch_data")

  f <- fit3(d, 3)
  testthat::expect_true(is.finite(f$quantities$jnll))
  testthat::expect_length(f$quantities$equil_catch_hat, n)
})

testthat::test_that("two equilibrium catches for one fleet are refused", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not(exists("GOA2018SS"))
  d <- Rceattle::GOA2018SS
  row <- d$catch_data[d$catch_data$Year == d$styr, ][1, ]
  row$Year <- d$styr - 1L
  e <- d; e$catch_data <- rbind(row, d$catch_data)
  testthat::expect_error(fit3(e, 3), "More than one initial equilibrium catch")
})

testthat::test_that("a non-positive equilibrium catch is refused", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not(exists("GOA2018SS"))
  d <- Rceattle::GOA2018SS
  e <- d
  i <- which(e$catch_data$Year == e$styr - 1L)
  testthat::skip_if(length(i) == 0)
  e$catch_data$Catch[i] <- 0
  testthat::expect_error(fit3(e, 3), "must be positive")
})
