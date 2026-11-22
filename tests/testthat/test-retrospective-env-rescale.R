# A peel used to drop env_data rows after its terminal year under rescale = TRUE.
# An ar1(1 | Year) linkage has one deviation per env_data year, so each peel's
# model was shorter than the parent's `inits` and retrospective() stopped:
# "beta_linkage (inits 46, model 45); beta_linkage_re (inits 40, model 39)" on a
# Pacific hake model with krill observed as a QAR1 latent on R0. rescale = TRUE
# also turned a bioenergetics temperature into z-scores, so rescale now names
# columns. A column with no retained value, blanked, made M1 NaN. A column
# listed in projected_env_data_known (a climate projection) is kept whole. A
# logical column, which the model reads as 0/1, was not withheld in 5.69.0.

testthat::test_that("a peel keeps an observed AR1 latent and rescales only named columns", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  d <- make_test_data(nyrs = 20, seed = 1)
  d$fleet_control$Catchability[1] <- "Estimated"
  yrs <- d$styr:d$endyr
  set.seed(2)
  d$env_data <- data.frame(Year = yrs, qcov = 3 + 2 * stats::rnorm(length(yrs)),
                           temp = 5 + stats::rnorm(length(yrs)),
                           warm = stats::rnorm(length(yrs)) > 0)
  qf <- Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ ar1(1 | Year), by = ~ fleet, fleet = 1L,
                               observe = "qcov", obs_sd = 0.3)))
  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 0,
    random_rec = FALSE, msmMode = 0, qFun = qf,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))

  testthat::expect_error(
    Rceattle::retrospective(fit, peels = 2, rescale = "krill", cores = 1),
    "not in `env_data`: krill", fixed = TRUE)
  testthat::expect_error(
    Rceattle::retrospective(fit, peels = 2, rescale = NA, cores = 1),
    "must be TRUE, FALSE", fixed = TRUE)
  testthat::expect_error(
    Rceattle::retrospective(fit, peels = 2, rescale = "warm", cores = 1),
    "cannot be standardized: warm", fixed = TRUE)
  testthat::expect_error(
    Rceattle::retrospective(fit, peels = 2, projected_env_data_known = "krill", cores = 1),
    "`projected_env_data_known` names column(s) that are not in `env_data`: krill",
    fixed = TRUE)
  testthat::expect_error(
    Rceattle::retrospective(fit, peels = 2, projected_env_data_known = NA, cores = 1),
    "`projected_env_data_known` must be TRUE, FALSE", fixed = TRUE)

  deep_env <- function(r) {
    testthat::expect_length(r$Rceattle_list, 3L)   # both peels plus the parent
    deep <- r$Rceattle_list[[paste0("Year_", d$endyr - 2)]]
    testthat::expect_length(deep$estimated_params$beta_linkage_re,
                            length(fit$estimated_params$beta_linkage_re))
    deep$data_list$env_data
  }
  z_kept <- function(x, post) {
    testthat::expect_equal(mean(x[!post]), 0, tolerance = 1e-10)
    testthat::expect_equal(stats::sd(x[!post]), 1, tolerance = 1e-10)
  }

  # Named: the observed covariate is withheld after the peel and standardized on
  # the retained years; the unnamed columns are withheld but not standardized.
  env <- deep_env(suppressWarnings(suppressMessages(Rceattle::retrospective(
    fit, peels = 2, rescale = "qcov", cores = 1, getsd = FALSE))))
  post <- env$Year > d$endyr - 2
  testthat::expect_true(all(is.na(env$qcov[post])))
  z_kept(env$qcov, post)
  testthat::expect_true(all(is.na(env$temp[post])))
  testthat::expect_identical(env$temp[!post], d$env_data$temp[!post])
  testthat::expect_identical(env$warm, replace(d$env_data$warm, post, NA))

  # TRUE: every numeric column is standardized. A column known ahead keeps its
  # post-peel values, on the retained years' centre and scale.
  env <- deep_env(suppressWarnings(suppressMessages(Rceattle::retrospective(
    fit, peels = 2, rescale = TRUE, cores = 1, getsd = FALSE,
    projected_env_data_known = "temp"))))
  testthat::expect_true(all(is.na(env$qcov[post])))
  z_kept(env$qcov, post)
  z_kept(env$temp, post)
  ctr <- mean(d$env_data$temp[!post]); scl <- stats::sd(d$env_data$temp[!post])
  testthat::expect_equal(env$temp, (d$env_data$temp - ctr) / scl, tolerance = 1e-12)
  # A logical flag is withheld but never standardized.
  testthat::expect_identical(env$warm, replace(d$env_data$warm, post, NA))
})

testthat::test_that("a peel keeps fixed-effect, grouping and no-retained-value columns whole", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  d <- make_test_data(nyrs = 20, seed = 1)
  d$fleet_control$Catchability[1] <- "Estimated"
  yrs <- d$styr:d$endyr
  set.seed(3)
  n <- length(yrs)
  d$env_data <- data.frame(Year = yrs, temp = 5 + stats::rnorm(n),
                           regime = ifelse(yrs < d$styr + 12, 1, 2),
                           late = c(rep(NA, n - 2), 1, 2), x = stats::rnorm(n))
  qf <- Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ temp + (1 | regime), by = ~ fleet, fleet = 1L)))
  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 0,
    random_rec = FALSE, msmMode = 0, qFun = qf,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))

  testthat::expect_error(
    Rceattle::retrospective(fit, peels = 3, rescale = "regime", cores = 1),
    "cannot be standardized: regime", fixed = TRUE)

  ws <- character(0)
  r <- withCallingHandlers(suppressMessages(Rceattle::retrospective(
    fit, peels = 3, rescale = TRUE, cores = 1, getsd = FALSE)),
    warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") })
  testthat::expect_length(r$Rceattle_list, 4L)
  testthat::expect_true(any(grepl("Covariate(s) temp enter", ws, fixed = TRUE)))
  testthat::expect_true(any(grepl("grouping column(s) regime", ws, fixed = TRUE)))
  testthat::expect_true(any(grepl("column(s) late have no value", ws, fixed = TRUE)))
  testthat::expect_true(any(grepl("leaves late in their own units", ws, fixed = TRUE)))

  env <- r$Rceattle_list[[paste0("Year_", d$endyr - 3)]]$data_list$env_data
  kept <- env$Year <= d$endyr - 3
  # The fixed-effect covariate keeps its post-peel values, standardized on the
  # retained years; the grouping label and the column with no retained value
  # are untouched.
  ctr <- mean(d$env_data$temp[kept]); scl <- stats::sd(d$env_data$temp[kept])
  testthat::expect_equal(env$temp, (d$env_data$temp - ctr) / scl, tolerance = 1e-12)
  testthat::expect_identical(env$regime, d$env_data$regime)
  testthat::expect_identical(env$late, d$env_data$late)

  # TRUE declares every column known ahead: no peel's env_data changes (`x`,
  # unlinked, would be blanked in every peel), and no kept-column warning fires.
  ws <- character(0)
  r <- withCallingHandlers(suppressMessages(Rceattle::retrospective(
    fit, peels = 3, cores = 1, getsd = FALSE, projected_env_data_known = TRUE)),
    warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") })
  testthat::expect_length(r$Rceattle_list, 4L)
  testthat::expect_false(any(grepl(
    "enter the model as fixed-effect|grouping column|have no value", ws)))
  for (peel in r$Rceattle_list) {
    testthat::expect_identical(peel$data_list$env_data, fit$data_list$env_data)
  }
})

testthat::test_that("rescale refuses Year and a column that is also a grouping variable", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # `t` is both a fixed-effect term and the AR1 grouping: standardized, its
  # levels would no longer match and observations would drop silently.
  d <- make_test_data(nyrs = 20, seed = 1)
  d$fleet_control$Catchability[1] <- "Estimated"
  yrs <- d$styr:d$endyr
  set.seed(4)
  d$env_data <- data.frame(Year = yrs, t = seq_along(yrs), qcov = stats::rnorm(length(yrs)))
  qf <- Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ t + ar1(1 | t), by = ~ fleet, fleet = 1L,
                               observe = "qcov", obs_sd = 0.3)))
  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 3,
    random_rec = FALSE, msmMode = 0, qFun = qf,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))

  testthat::expect_identical(Rceattle:::.rce_linkage_fixed_covariates(
    fit$data_list, re_groups = TRUE, fixed = FALSE), "t")
  for (col in c("t", "Year")) {
    testthat::expect_error(
      Rceattle::retrospective(fit, peels = 2, rescale = col, cores = 1),
      paste0("cannot be standardized: ", col), fixed = TRUE)
  }
})

testthat::test_that("env_data rows with no Year or before styr do not stop a peel", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # `y` has values only before styr and after the peel: the model years hold
  # none, so the peel keeps it whole rather than blank it into a NaN M1.
  d <- make_test_data(nyrs = 20, seed = 1)
  yrs <- d$styr:d$endyr
  n <- length(yrs)
  d$env_data <- data.frame(Year = c(d$styr - 1, yrs, NA), x = c(0, seq_along(yrs), 0),
                           y = c(1, rep(NA, n - 1), 2, NA))
  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 0,
    random_rec = FALSE, msmMode = 0,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
  r <- suppressWarnings(suppressMessages(Rceattle::retrospective(
    fit, peels = 1, cores = 1, getsd = FALSE)))
  testthat::expect_length(r$Rceattle_list, 2L)
})
