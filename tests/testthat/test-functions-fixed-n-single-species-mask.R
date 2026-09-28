# In single-species mode a species with input numbers-at-age (estDynamics > 0)
# has rec_pars fixed, so its equilibrium SB0 and B0 are built on placeholder
# recruitment (R0 = exp(9)), and with DynamicHCR = FALSE ssb_depletion and
# biomass_depletion divide by them. Those four are NA for such a species, and
# its R is the input recruits (first-age NByageFixed). Under DynamicHCR = TRUE the depletions are DynamicSB0-based, which for input
# numbers is the numbers relative to themselves, and stay reported. Under
# predation MSSB0 replaces SB0 (test-functions-fixed-n-hcr-outputs.R).
# A model whose only species is fixed has nothing to estimate, so the fixture is
# the two-species one run without predation.

testthat::skip_on_cran()

testthat::test_that("placeholder SB0, B0 and depletion are NA in single-species mode", {
  d <- fixed_natage_data(c(1, 0))
  m <- fixed_natage_build(d, msmMode = 0)
  q <- m$quantities
  for (nm in c("SB0", "B0", "ssb_depletion", "biomass_depletion", "R0")) {
    testthat::expect_true(all(is.na(q[[nm]][1, ])), info = nm)
    testthat::expect_true(all(is.finite(q[[nm]][2, ])), info = nm)
  }
  # R is the input recruits, read back from the workbook's first age (pop_scalar = 1).
  nb <- d$NByageFixed[order(d$NByageFixed$Year), ]
  testthat::expect_equal(unname(q$R[1, ]), nb[["Age 1"]])
  testthat::expect_true(all(is.finite(q$R[2, ])))
  # The fixture's input recruits equal the exp(9) placeholder they replace, so
  # scale them: R must follow the input, and a year with no row is NA, not 0.
  d3 <- d
  d3$NByageFixed[["Age 1"]] <- 3 * d3$NByageFixed[["Age 1"]]
  d3$NByageFixed <- d3$NByageFixed[d3$NByageFixed$Year <= d3$endyr, ]
  q3 <- fixed_natage_build(d3, msmMode = 0)$quantities
  nh <- d$endyr - d$styr + 1
  testthat::expect_equal(unname(q3$R[1, seq_len(nh)]), 3 * nb[["Age 1"]][seq_len(nh)])
  testthat::expect_true(all(is.na(q3$R[1, -seq_len(nh)])))

  # Its R carries no standard error in the table either: give the fit a fake sdreport.
  n <- length(q$R)
  m_sd <- m
  m_sd$sdrep <- list(value = stats::setNames(rep(1, n), rep("R", n)), sd = rep(0.5, n))
  r <- as.data.frame(m_sd, which = "R")
  testthat::expect_gt(nrow(r), 0)
  testthat::expect_true(all(is.na(r$se[r$species == m$data_list$spnames[1]])))
  testthat::expect_true(all(r$se[r$species == m$data_list$spnames[2]] == 0.5))
  for (nm in c("R_init", "steepness", "SPR0")) testthat::expect_true(is.na(q[[nm]][1]), info = nm)
  testthat::expect_true(all(is.finite(q$ssb[1, ])))
  testthat::expect_true(all(is.finite(q$DynamicSB0[1, ])))

  # sample_rec() reads R and R0 for every species it projects; the fixed one is skipped.
  s <- sample_rec(m, sample_rec = FALSE, update_model = FALSE)
  testthat::expect_true(all(is.finite(s$estimated_params$rec_dev)))
})

testthat::test_that("estDynamics = 2 in single-species mode is announced as fitting as 1", {
  d <- fixed_natage_data(c(2, 0))
  d$msmMode <- 0
  d <- suppressMessages(Rceattle::switch_check(d))
  testthat::expect_message(suppressWarnings(Rceattle:::data_check(d)),
                           "estimated only under predation")
})

testthat::test_that("the fixed species' R carries no standard error for the band", {
  # sdrep series are [nspp, nyrs] flattened column-major: species 1 is every odd cell.
  fake <- list(sdrep = list(value = stats::setNames(1:6, rep("R", 6)), sd = 1:6),
               data_list = list(nspp = 2, estDynamics = c(1, 0)))
  testthat::expect_equal(Rceattle:::.rce_series_sd(fake, "R", 6), c(NA, 2, NA, 4, NA, 6))
  fake$data_list$estDynamics <- c(0, 0)
  testthat::expect_equal(Rceattle:::.rce_series_sd(fake, "R", 6), 1:6)
})

testthat::test_that("under DynamicHCR = TRUE the depletions stay reported", {
  m <- fixed_natage_build(fixed_natage_data(c(1, 0)), msmMode = 0,
                          HCR = build_hcr(HCR = "ConstantF", DynamicHCR = TRUE, Ftarget = 0.1))
  q <- m$quantities
  testthat::expect_true(all(is.na(q$SB0[1, ])))
  testthat::expect_true(all(is.finite(q$ssb_depletion[1, ])))
})
