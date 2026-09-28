# rec_pars alpha and beta are bounded at +/-30 on the log scale (5.39.0). A fit
# saved on the old unbounded ridge (log alpha reached 702) starts outside them,
# so build_bounds() stops -- and with it retrospective(), profile() and
# run_mse(), which refit from the saved values. The error has to name the block
# and the bound; "Initial parameter values are not within bounds" alone sends
# the reader looking at the data.

testthat::test_that("out-of-bounds rec_pars names the block and says where the bound came from", {
  d <- suppressMessages(Rceattle::switch_check(Rceattle::BS2017SS))
  p <- suppressWarnings(Rceattle::build_params(d))
  p$rec_pars[1, 2] <- 702                      # the measured flat-ridge value
  err <- tryCatch(suppressMessages(
    Rceattle::build_bounds(param_list = p, data_list = d)),
    error = function(e) conditionMessage(e))
  testthat::expect_type(err, "character")
  testthat::expect_match(err, "rec_pars", fixed = TRUE)
  testthat::expect_match(err, "5.39.0", fixed = TRUE)
  testthat::expect_match(err, "srr_alpha_init", fixed = TRUE)

  # A violation in another block keeps the plain message, with no curve advice.
  p2 <- suppressWarnings(Rceattle::build_params(d))
  p2$rec_dev[1, 1] <- 99                       # bounded at +/-15
  err2 <- tryCatch(suppressMessages(
    Rceattle::build_bounds(param_list = p2, data_list = d)),
    error = function(e) conditionMessage(e))
  testthat::expect_match(err2, "rec_dev", fixed = TRUE)
  testthat::expect_false(grepl("5.39.0", err2, fixed = TRUE))
})
