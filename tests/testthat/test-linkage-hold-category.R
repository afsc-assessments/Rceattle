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
  block <- function(endyr_peel, endyr_full, brk = 18) {
    d$endyr_peel <- endyr_peel
    d$endyr_full <- endyr_full
    f <- stats::as.formula(sprintf("~ cut(Year, c(-Inf, %d, Inf), dig.lab = 4)", brk))
    nm <- colnames(stats::model.matrix(f, d$env_data))[2]
    fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
      data_list = d, inits = NULL, file = NULL, estimateMode = 3,
      random_rec = FALSE, msmMode = 0,
      M1Fun = Rceattle::build_M1(M1_model = 1, linkages = list(
        M1 = Rceattle::linkage_spec(f, est_phase = 0, init = stats::setNames(list(0.3), nm)))),
      fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
    fit$data_list$linkage_X[yrs > brk, nm]
  }
  testthat::expect_true(all(block(17, d$endyr) == 0))   # a peel ending in 17
  testthat::expect_true(all(block(17, 18) == 1))        # left over; endyr is 20
  # The unpeeled model (endyr_peel = endyr): a break after endyr is not held.
  testthat::expect_true(all(block(d$endyr, d$endyr, brk = d$endyr + 2) == 1))
})

# A q or selectivity linkage deviation holds its level after the hold year too, as
# a fleet_control deviation does; a fixed covariate does not.
testthat::test_that("a held random effect keeps its level after hold_after", {
  set.seed(4)
  env <- data.frame(Year = 2000:2015, temp = stats::rnorm(16))
  f <- ~ temp + ar1(1 | Year)
  X0 <- Rceattle:::.linkage_design(f, env, hold_after = 2009)$X
  X1 <- Rceattle:::.linkage_design(f, env, hold_after = 2009, hold_re = TRUE)$X
  later <- env$Year > 2009
  fixed <- 1:2   # (Intercept), temp
  re <- seq(3, ncol(X1))
  testthat::expect_identical(colnames(X1), colnames(X0))
  testthat::expect_identical(X1[!later, ], X0[!later, ])
  testthat::expect_identical(X1[later, fixed], X0[later, fixed])
  testthat::expect_true(all(t(X1[later, re]) == X1[env$Year == 2009, re]))
  testthat::expect_false(identical(X1[later, re], X0[later, re]))
  # A formula with no random effect holds nothing and keeps its columns.
  testthat::expect_identical(
    Rceattle:::.linkage_design(~ temp, env, hold_after = 2009, hold_re = TRUE)$X,
    Rceattle:::.linkage_design(~ temp, env, hold_after = 2009)$X)
})

# A survey 50% higher in years 15-18, so the peel ending in year 18 has a nonzero
# terminal q deviation.
.bumped_survey <- function() {
  d <- make_test_data(nyrs = 20, seed = 1)
  d$fleet_control$Catchability[1] <- "Estimated"
  d$env_data <- data.frame(Year = d$styr:d$projyr)
  late <- d$index_data$Year > d$endyr - 6 & d$index_data$Year <= d$endyr - 2
  d$index_data$Observation[late] <- 1.5 * d$index_data$Observation[late]
  d
}

testthat::test_that("a peel holds a q linkage deviation at its terminal year", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  d <- .bumped_survey()
  endyr_peel <- d$endyr - 2
  k <- endyr_peel - d$styr + 1
  for (f in c(~ ar1(1 | Year), ~ (1 | Year))) {
    fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
      data_list = d, inits = NULL, file = NULL, estimateMode = 0,
      random_rec = FALSE, msmMode = 0,
      qFun = Rceattle::build_catchability(linkages = list(q = Rceattle::linkage_spec(
        f, by = ~ fleet, fleet = 1L, link = "log", init = list(sigma = 0.2)))),
      fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
    peel <- suppressWarnings(suppressMessages(Rceattle::retrospective(
      fit, peels = 2, cores = 1, getsd = FALSE, phase = FALSE)))$Rceattle_list[[paste0("Year_", endyr_peel)]]
    q <- peel$quantities$index_q[1, ]
    testthat::expect_equal(q[(k + 1):(d$endyr - d$styr + 1)], rep(q[[k]], 2),
                           ignore_attr = TRUE, info = deparse(f))
  }
})

# A held q effect and an unheld R0 effect on the same grouping column pool apart,
# under their own names: pooled by name, the peel stopped ("design column
# 'Year_re::19' has inconsistent values"), and renamed, an `init` on a level did.
testthat::test_that("a peel holds a q effect that shares its grouping column with R0", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  d <- .bumped_survey()
  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 0,
    random_rec = FALSE, msmMode = 0,
    qFun = Rceattle::build_catchability(linkages = list(q = Rceattle::linkage_spec(
      ~ (1 | Year), by = ~ fleet, fleet = 1L, link = "log",
      init = list(sigma = 0.2, `Year_re::5` = 0.1)))),
    recFun = Rceattle::build_srr(linkages = list(R0 = Rceattle::linkage_spec(~ (1 | Year)))),
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
  r <- suppressWarnings(suppressMessages(Rceattle::retrospective(
    fit, peels = 2, cores = 1, getsd = FALSE, phase = FALSE)))
  testthat::expect_length(r$Rceattle_list, 3L)
  peel <- r$Rceattle_list[[paste0("Year_", d$endyr - 2)]]
  q <- peel$quantities$index_q[1, ]
  k <- d$endyr - 2 - d$styr + 1
  testthat::expect_equal(q[[k + 1]], q[[k]])
  # The held q columns sit beside R0's under the same names.
  testthat::expect_setequal(colnames(peel$data_list$linkage_X),
                            colnames(fit$data_list$linkage_X))
})

# A q linkage with no random effect (a covariate or a prior on its mean) has
# nothing to hold, and every peel still runs.
testthat::test_that("a peel runs with a q linkage that has no random effect", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  d <- make_test_data(nyrs = 20, seed = 1)
  d$fleet_control$Catchability[1] <- "Estimated"
  set.seed(2)
  d$env_data <- data.frame(Year = d$styr:d$projyr, temp = stats::rnorm(30))
  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 0,
    random_rec = FALSE, msmMode = 0,
    qFun = Rceattle::build_catchability(linkages = list(q = Rceattle::linkage_spec(
      ~ temp, by = ~ fleet, fleet = 1L))),
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
  r <- suppressWarnings(suppressMessages(Rceattle::retrospective(
    fit, peels = 1, cores = 1, getsd = FALSE, phase = FALSE)))
  testthat::expect_length(r$Rceattle_list, 2L)
})

# A selectivity linkage that is a category alone (5.70.1's block hold) has no
# random effect to hold, and every peel still runs.
testthat::test_that("a peel runs with a selectivity block and no random effect", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  d <- make_test_data(nyrs = 20, seed = 1)
  d$env_data <- data.frame(Year = d$styr:d$projyr)
  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 0,
    random_rec = FALSE, msmMode = 0,
    selFun = Rceattle::build_selectivity(linkages = list(inf_asc = Rceattle::linkage_spec(
      ~ cut(Year, c(0, 10, 100)), by = ~ fleet, fleet = 2L))),
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
  r <- suppressWarnings(suppressMessages(Rceattle::retrospective(
    fit, peels = 2, cores = 1, getsd = FALSE, phase = FALSE)))
  testthat::expect_length(r$Rceattle_list, 3L)
})
