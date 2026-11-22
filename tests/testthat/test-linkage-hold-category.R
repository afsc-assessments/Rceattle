# In a retrospective peel, a categorical fixed-effect linkage term (a time block,
# a regime, a flag) holds its value from the peel's terminal year, so the peel's
# forecast never reads a block that starts after it. A fit's own projection reads
# env_data as given.

testthat::test_that("a category holds its value after hold_after; numbers do not", {
  set.seed(3)
  env <- data.frame(Year = 2000:2015, temp = stats::rnorm(16),
                    reg = rep(1:4, each = 4), ch = rep(c("a", "b"), each = 8))
  f <- ~ temp + factor(reg) + cut(Year, c(-Inf, 2010, 2013, Inf)) + ch +
    I(Year > 2011) + temp:cut(Year, c(-Inf, 2010, Inf)) + poly(temp, 2)
  design <- function(h) Rceattle:::.linkage_design(f, env, hold_after = h)$X
  X0 <- design(NULL)
  X1 <- design(2009)
  ref <- stats::model.matrix(Rceattle:::.parse_linkage_formula(f)$fixed, data = env)
  testthat::expect_identical(X0, ref[, , drop = FALSE])

  # The same columns, including factor(reg)4, which no year keeps after the hold.
  testthat::expect_identical(colnames(X1), colnames(X0))
  pre <- env$Year <= 2009
  testthat::expect_identical(X1[pre, ], X0[pre, ])
  num <- c("temp", "poly(temp, 2)1", "poly(temp, 2)2")
  testthat::expect_identical(X1[!pre, num], X0[!pre, num])
  cats <- setdiff(colnames(X1), c("(Intercept)", num, grep("^temp:", colnames(X1), value = TRUE)))
  testthat::expect_true(all(t(X1[!pre, cats]) == X1[env$Year == 2009, cats]))
  testthat::expect_true(all(X1[, "factor(reg)4"] == 0))
  # An interaction holds the block, not the covariate.
  ix <- grep("^temp:", colnames(X1), value = TRUE)
  testthat::expect_equal(unname(X1[!pre, ix]),
                         env$temp[!pre] * unname(X0[env$Year == 2009, ix]))

  # No row for the hold year: nothing is held.
  testthat::expect_identical(design(1990), X0)
})

testthat::test_that("a peel holds its terminal-year block; a projection does not", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  d <- make_test_data(nyrs = 20, seed = 1)
  yrs <- d$styr:d$projyr
  m1_block <- function(brk) {
    d$env_data <- data.frame(Year = yrs)
    f <- stats::as.formula(sprintf("~ cut(Year, c(-Inf, %d, Inf), dig.lab = 4)", brk))
    nm <- colnames(stats::model.matrix(f, d$env_data))[2]
    suppressWarnings(suppressMessages(Rceattle::fit_mod(
      data_list = d, inits = NULL, file = NULL, estimateMode = 0,
      random_rec = FALSE, msmMode = 0,
      M1Fun = Rceattle::build_M1(M1_model = 1, linkages = list(
        M1 = Rceattle::linkage_spec(f, est_phase = 0, init = stats::setNames(list(0.3), nm)))),
      fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
  }
  m1 <- function(fit) fit$quantities$M1_at_age[1, 1, 1, ]

  # A break after endyr: the projection follows it.
  fit <- m1_block(d$endyr + 2)
  testthat::expect_gt(m1(fit)[yrs == d$endyr + 3], m1(fit)[yrs == d$endyr])

  # A break two years before endyr: a peel ending before it holds its own block.
  fit <- m1_block(d$endyr - 2)
  testthat::expect_gt(m1(fit)[yrs == d$endyr], m1(fit)[yrs == d$endyr - 2])
  r <- suppressWarnings(suppressMessages(Rceattle::retrospective(
    fit, peels = 3, cores = 1, getsd = FALSE)))
  deep <- m1(r$Rceattle_list[[paste0("Year_", d$endyr - 3)]])
  testthat::expect_true(all(deep[yrs > d$endyr - 3] == deep[yrs == d$endyr - 3]))
})

testthat::test_that("endyr_peel sets the hold year only inside a peel", {
  # retrospective() leaves endyr_peel and endyr_full on the models it returns; a
  # refit at another endyr (run_mse() advances it) is not a peel and holds nothing.
  d <- make_test_data(nyrs = 20, seed = 1)
  yrs <- d$styr:d$projyr
  d$env_data <- data.frame(Year = yrs)
  f <- ~ cut(Year, c(-Inf, 18, Inf), dig.lab = 4)
  nm <- colnames(stats::model.matrix(f, d$env_data))[2]
  block <- function(endyr_peel, endyr_full) {
    d$endyr_peel <- endyr_peel
    d$endyr_full <- endyr_full
    fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
      data_list = d, inits = NULL, file = NULL, estimateMode = 3,
      random_rec = FALSE, msmMode = 0,
      M1Fun = Rceattle::build_M1(M1_model = 1, linkages = list(
        M1 = Rceattle::linkage_spec(f, est_phase = 0, init = stats::setNames(list(0.3), nm)))),
      fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
    fit$data_list$linkage_X[yrs > 18, nm]
  }
  testthat::expect_true(all(block(17, d$endyr) == 0))   # a peel ending in 17
  testthat::expect_true(all(block(17, 18) == 1))        # left over; endyr is 20
  testthat::expect_true(all(block(d$endyr, d$endyr) == 1))  # the unpeeled model
})
