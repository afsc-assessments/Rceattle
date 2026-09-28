# Sel_norm_scope is read only by the shared normalizer in selectivity.hpp. Hake
# normalizes each sex to its own maximum in its own block, and LogisticPM reads
# Sel_norm_bin as a penalty range, so on a two-sex fleet of either form the
# column changes nothing: measured on GOAatf fleet 3 (Hake, two sexes), the two
# scopes gave identical selectivity and objective (46.4710) to every digit.

testthat::test_that("data_check() says Sel_norm_scope is not read on a two-sex Hake or LogisticPM fleet", {
  for (sel in c("Hake", "LogisticPM")) {
    d <- Rceattle::GOAatf2023
    d$fleet_control$Selectivity[3] <- sel
    d <- suppressMessages(Rceattle::switch_check(d))
    testthat::expect_message(suppressWarnings(Rceattle:::data_check(d)),
                             "'Sel_norm_scope' is not read", info = sel)
  }
})

# The two forms reach that conclusion for opposite reasons, and saying Hake's
# reason for LogisticPM would tell the user the level is pinned when it is free.
testthat::test_that("the notice gives each form its own reason", {
  d <- Rceattle::GOAatf2023
  d$fleet_control$Selectivity[3] <- "Hake"
  d <- suppressMessages(Rceattle::switch_check(d))
  testthat::expect_message(suppressWarnings(Rceattle:::data_check(d)),
                           "normalizes each sex to its own maximum")

  d <- Rceattle::GOAatf2023
  d$fleet_control$Selectivity[3] <- "LogisticPM"
  d <- suppressMessages(Rceattle::switch_check(d))
  testthat::expect_message(suppressWarnings(Rceattle:::data_check(d)),
                           "does not normalize at all")
  # LogisticPM must not inherit Hake's consequence: it does not constrain level.
  testthat::expect_no_message(suppressWarnings(Rceattle:::data_check(d)),
                              message = "cannot differ in selectivity level")
})

testthat::test_that("the notice does not fire on a one-sex species", {
  d <- Rceattle::BS2017SS
  d$fleet_control$Selectivity[1] <- "Hake"
  d <- suppressMessages(Rceattle::switch_check(d))
  testthat::expect_no_message(suppressWarnings(Rceattle:::data_check(d)),
                              message = "Sel_norm_scope")
})
