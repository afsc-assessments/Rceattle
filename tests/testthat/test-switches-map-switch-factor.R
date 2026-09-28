# .map_switch() returned a factor unchanged, and every downstream comparison of
# a factor against an integer code is FALSE rather than an error, so a factor
# switch fitted as its level index: factor("LognormalPrior") is level 1, which is
# code 1 ("Estimated"), and the prior was silently dropped. read.csv() with
# stringsAsFactors = TRUE hands back exactly that.

testthat::test_that("a factor switch resolves to the code its label names", {
  m <- Rceattle:::.map_switch(factor(c("FixedScaled", "Estimated")),
                              Rceattle:::estDynamics_map, "estDynamics")
  testthat::expect_equal(m, c(2, 0))
})

testthat::test_that("a numeric-looking string or factor is the code it names", {
  m <- Rceattle:::.map_switch(c("2", "Fixed", NA), Rceattle:::estDynamics_map, "estDynamics")
  testthat::expect_equal(m, c(2, 1, NA))
  f <- Rceattle:::.map_switch(factor(c("1", "0")), Rceattle:::estimate_sd_map, "Estimate_index_sd")
  testthat::expect_equal(f, c(1, 0))
  testthat::expect_error(Rceattle:::.map_switch("99", Rceattle:::suitMode_map, "suitMode"),
                         "Invalid 'suitMode' value\\(s\\): 99")
})

testthat::test_that("a factor Estimate_index_sd column reaches switch_check() as its code", {
  d <- Rceattle::BS2017SS
  d$fleet_control$Estimate_index_sd <- factor(rep("Fixed", nrow(d$fleet_control)))
  out <- suppressMessages(Rceattle::switch_check(d))
  testthat::expect_false(is.factor(out$fleet_control$Estimate_index_sd))
  testthat::expect_true(all(out$fleet_control$Estimate_index_sd %in% c("Fixed", 0)))
})
