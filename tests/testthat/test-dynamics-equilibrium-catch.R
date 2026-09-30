# The initial equilibrium catch is a catch_data row at Year == styr - 1, read only
# under initMode 6 (FishedNonEquilibriumSelected).
#
# That year is not a free marker: GOA2018SS carries 23 catch rows before styr, two
# of them on 1976, and it is the only bundled dataset that does. Reading those as
# equilibrium observations under a mode that holds Finit at 0 predicts 0, takes
# log(0), and returns a non-finite objective on both GOA golden references, which
# is how this surfaced. A negative Year cannot be the marker either --
# run_mse() reserves those for rows it splices in as the next assessment's data,
# and its window filters are on abs(Year).
#
# The rows are held in data_list$equil_catch_data rather than left in catch_data,
# because catch_hat, catch_sd, plot_catch(), residuals(), sim_mod() and run_mse()
# all index catch_data by position.

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

# Modes 3 and 4 estimate Finit but reach an age with it differently -- 3
# accumulates a flat Finit over ages, 4 applies it once -- so the equilibrium
# catch, which is Baranov at Finit * selectivity, would be scored against a
# population that mortality never produced.
testthat::test_that("the other fished modes do not read it", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not(exists("GOA2018SS"))
  d <- Rceattle::GOA2018SS
  for (im in c(3, 4)) {
    f <- fit3(d, im)
    testthat::expect_length(f$quantities$equil_catch_hat, 0)
    testthat::expect_true(is.finite(f$quantities$jnll))
    testthat::expect_false(any(is.na(f$map$mapFactor$log_Finit)))  # still estimated
  }
})

# The template scores an equilibrium catch on flt_type == 1 alone, so a survey's
# row would leave the catch history and be fitted by nothing.
testthat::test_that("an equilibrium catch on a survey is refused", {
  testthat::skip_if_not(exists("GOA2018SS"))
  d <- Rceattle::GOA2018SS; d$initMode <- 6
  srv <- d$fleet_control$Fleet_code[d$fleet_control$Fleet_type == 2][1]
  testthat::skip_if(is.na(srv))
  row <- d$catch_data[d$catch_data$Year == d$styr, ][1, ]
  row$Fleet_code <- srv
  row$Species <- d$fleet_control$Species[d$fleet_control$Fleet_code == srv][1]
  row$Year <- d$styr - 1L
  e <- d; e$catch_data <- rbind(row, d$catch_data)
  testthat::expect_error(
    suppressMessages(suppressWarnings(
      Rceattle:::data_check(Rceattle::switch_check(e)))),
    "not fisheries")
})

# catch_hat has one entry per catch_data row, and plot_catch(), residuals(),
# sim_mod() and run_mse() all pair the two by position. GOA2018SS is the only
# bundled dataset with a styr - 1 catch row, so it is the only one that can
# catch the equilibrium rows being left in catch_data.
testthat::test_that("catch_data stays aligned with catch_hat", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not(exists("GOA2018SS"))
  d <- Rceattle::GOA2018SS
  testthat::expect_gt(sum(d$catch_data$Year == d$styr - 1L, na.rm = TRUE), 0)

  for (im in c(2, 6)) {
    f <- fit3(d, im)
    testthat::expect_length(f$quantities$catch_hat, nrow(f$data_list$catch_data))
    testthat::expect_length(f$quantities$catch_sd,  nrow(f$data_list$catch_data))
    testthat::expect_false(any(f$data_list$catch_data$Year == d$styr - 1L, na.rm = TRUE))
    # The three consumers that pair them, each of which errored on the mismatch.
    testthat::expect_s3_class(Rceattle::plot_catch(f), "ggplot")
    testthat::expect_s3_class(stats::residuals(f, source = "catch"), "data.frame")
    testthat::expect_type(Rceattle::sim_mod(f, simulate = FALSE), "list")
  }
})

# clean_data() runs on every way into fit_mod(), so retrospective(), jitter(),
# self_test(), model_average() and run_mse()'s refits all clean a data_list that
# has been cleaned once already. If the second pass rebuilt equil_catch_data
# from catch_data alone it would come back empty, and the refit would quietly
# drop the observation: 51,873 nats on GOA2018SS at initMode 6.
testthat::test_that("cleaning twice keeps the equilibrium catch", {
  testthat::skip_if_not(exists("GOA2018SS"))
  d <- Rceattle::GOA2018SS
  n <- sum(d$catch_data$Year == d$styr - 1L, na.rm = TRUE)
  cl <- suppressMessages(suppressWarnings(
    Rceattle::clean_data(Rceattle::switch_check(d))))
  testthat::expect_equal(nrow(cl$equil_catch_data), n)
  testthat::expect_equal(nrow(suppressMessages(
    Rceattle::clean_data(cl))$equil_catch_data), n)

  testthat::skip_if_not_installed("TMB")
  f1 <- fit3(d, 6)
  f2 <- fit3(f1$data_list, 6)                 # what a refit path hands back in
  testthat::expect_equal(f2$quantities$equil_catch_hat,
                         f1$quantities$equil_catch_hat)
  testthat::expect_equal(f2$quantities$jnll, f1$quantities$jnll)
})

# A row added to catch_data after the split must still be found, or it is
# neither fitted as catch nor as an equilibrium catch.
testthat::test_that("a row added after the split is not masked", {
  testthat::skip_if_not(exists("GOA2018SS"))
  d <- Rceattle::GOA2018SS
  cl <- suppressMessages(suppressWarnings(
    Rceattle::clean_data(Rceattle::switch_check(d))))
  extra <- cl$catch_data[1, ]
  extra$Fleet_code <- 14L; extra$Year <- d$styr - 1L; extra$Catch <- 5000
  cl$catch_data <- rbind(extra, cl$catch_data)
  testthat::expect_equal(
    nrow(Rceattle:::.rce_equil_catch_candidates(cl)),
    nrow(cl$equil_catch_data) + 1L)
})

# The rows are out of catch_data by the time data_check() runs on the fit path,
# so the loops over catch_data no longer see them. The template indexes flt_type
# and flt_units by Fleet_code with no range test.
testthat::test_that("an equilibrium catch gets the checks a catch row gets", {
  testthat::skip_if_not(exists("GOA2018SS"))
  d <- Rceattle::GOA2018SS; d$initMode <- 6
  cl <- suppressMessages(suppressWarnings(
    Rceattle::clean_data(Rceattle::switch_check(d))))
  bad <- function(f) {
    e <- cl; e$equil_catch_data <- f(e$equil_catch_data)
    suppressMessages(suppressWarnings(Rceattle:::data_check(e)))
  }
  testthat::expect_error(bad(function(x) { x$Fleet_code[1] <- 99L; x }),
                         "not in\\s+fleet_control")
  testthat::expect_error(bad(function(x) { x$Log_sd[1] <- 0; x }), "Log_sd")
  testthat::expect_error(bad(function(x) { x$Species[1] <- 7L; x }), "Species")
})

# The row is an ordinary catch row, so write_data() puts it back in the catch
# sheet; a data_list element with no workbook support round-trips to nothing.
testthat::test_that("the equilibrium catch survives a workbook round trip", {
  testthat::skip_if_not_installed("openxlsx")
  testthat::skip_if_not(exists("GOA2018SS"))
  d <- Rceattle::GOA2018SS
  cl <- suppressMessages(suppressWarnings(
    Rceattle::clean_data(Rceattle::switch_check(d))))
  testthat::expect_equal(nrow(cl$equil_catch_data),
                         sum(d$catch_data$Year == d$styr - 1L, na.rm = TRUE))

  f <- tempfile(fileext = ".xlsx"); on.exit(unlink(f))
  suppressMessages(suppressWarnings(Rceattle::write_data(cl, file = f)))
  back <- suppressMessages(suppressWarnings(Rceattle::read_data(f)))

  want <- cl$equil_catch_data[order(cl$equil_catch_data$Fleet_code),
                              c("Fleet_code", "Year", "Catch")]
  got <- back$catch_data[back$catch_data$Year == d$styr - 1L,
                         c("Fleet_code", "Year", "Catch")]
  got <- got[order(got$Fleet_code), ]
  testthat::expect_equal(unname(as.matrix(got)), unname(as.matrix(want)))
})

testthat::test_that("a Finit-estimating mode does read it, and says so", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not(exists("GOA2018SS"))
  d <- Rceattle::GOA2018SS
  n <- sum(d$catch_data$Year == d$styr - 1L, na.rm = TRUE)

  testthat::expect_message(
    suppressWarnings(Rceattle::fit_mod(
      data_list = d, file = NULL, inits = NULL, estimateMode = 3,
      random_rec = FALSE, msmMode = 0, initMode = 6,
      fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))),
    "Initial equilibrium catch read from catch_data")

  f <- fit3(d, 6)
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
  testthat::expect_error(fit3(e, 6), "More than one initial equilibrium catch")
})

testthat::test_that("a non-positive equilibrium catch is refused", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not(exists("GOA2018SS"))
  d <- Rceattle::GOA2018SS
  e <- d
  i <- which(e$catch_data$Year == e$styr - 1L)
  testthat::skip_if(length(i) == 0)
  e$catch_data$Catch[i] <- 0
  testthat::expect_error(fit3(e, 6), "must be positive")
})
