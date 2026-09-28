# fit_mod(config =) replaced data_list$model_config wholesale, so a config whose
# model_config was at the defaults -- run_config(model_config(), ...) or any
# load_config() output, which reconstructs defaults for every unwritten field --
# silently dropped every linkage on the data (57 random effects became 0 on a
# GOAatf survey rw(1 | Year) fit). The overlay is now field by field.

testthat::skip_on_cran()

survey_code <- function(d) d$fleet_control$Fleet_code[d$fleet_control$Fleet_name == "Survey"]

config_fixture <- function() {
  d <- make_test_data()
  if (!is.null(d$data_list)) d <- d$data_list
  srv <- survey_code(d)
  d$fleet_control$Catchability[d$fleet_control$Fleet_code == srv] <- "Estimated"
  q <- build_catchability(linkages = list(
    q = linkage_spec(~ (1 | Year), by = ~ fleet, fleet = srv)))
  build_data(base = d, model_config = model_config(qFun = q))
}

build <- function(d, ...) suppressWarnings(suppressMessages(fit_mod(
  data_list = d, inits = NULL, estimateMode = 3, msmMode = 0, random_rec = FALSE, ...,
  fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE))))

testthat::test_that("a default-model_config config keeps the data's linkages", {
  d <- config_fixture()
  n_re <- function(m) sum(names(m$obj$env$par)[m$obj$env$random] == "beta_linkage_re")
  m0 <- build(d)
  testthat::expect_gt(n_re(m0), 0)
  cfg <- run_config(model_config(), estimateMode = "DebugBuild")
  testthat::expect_equal(n_re(build(d, config = cfg)), n_re(m0))
  cfg2 <- run_config(d, estimateMode = "DebugBuild")
  testthat::expect_equal(n_re(build(d, config = cfg2)), n_re(m0))
})

testthat::test_that("a field the config set wins, defaults included; an unset one keeps the data's", {
  d <- make_test_data()
  if (!is.null(d$data_list)) d <- d$data_list
  d <- build_data(base = d, model_config = model_config(initMode = "OffsetEquilibrium"))
  # Unset in the config: the data's value stands.
  m <- build(d, config = run_config(model_config(), estimateMode = "DebugBuild"))
  testthat::expect_equal(m$run_config$model_config$initMode, "OffsetEquilibrium")
  # Set in the config to the default: it wins, with a warning that the two differ.
  cfg <- run_config(model_config(initMode = "NonEquilibrium"), estimateMode = "DebugBuild")
  testthat::expect_warning(
    m2 <- suppressMessages(fit_mod(data_list = d, inits = NULL, estimateMode = 3, msmMode = 0,
                                   random_rec = FALSE, config = cfg,
                                   fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE))),
    "`initMode` differs between the data's model_config and `config`")
  testthat::expect_equal(m2$run_config$model_config$initMode, "NonEquilibrium")
  # A config saved from a fit sets every field, so it reproduces that fit on
  # a data object carrying a different structure.
  d0 <- make_test_data(); if (!is.null(d0$data_list)) d0 <- d0$data_list
  f0 <- build(d0)
  tf <- tempfile(fileext = ".yaml")
  suppressMessages(save_config(f0, tf))
  m3 <- suppressWarnings(build(d, config = load_config(tf)))
  testthat::expect_equal(m3$run_config$model_config$initMode, "NonEquilibrium")
})

testthat::test_that("a config saved from a fit and reloaded onto the same data does not warn", {
  # The reloaded linkage is the same spec with a new formula environment.
  d <- config_fixture()
  tf <- tempfile(fileext = ".yaml")
  suppressMessages(save_config(build(d), tf))
  testthat::expect_no_warning(
    suppressMessages(fit_mod(data_list = d, inits = NULL, estimateMode = 3, msmMode = 0,
                             random_rec = FALSE, config = load_config(tf),
                             fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE))),
    message = "differs between")
})

testthat::test_that("a field both set differently is taken from the config with a warning", {
  d <- config_fixture()
  # The config's q linkage is a fixed intercept with a prior, not a random walk.
  cfg <- run_config(model_config(qFun = build_catchability(linkages = list(
    q = linkage_spec(~ 1, by = ~ fleet, fleet = survey_code(d),
                     priors = list("(Intercept)" = prior_normal(0, 1)))))),
    estimateMode = "DebugBuild")
  testthat::expect_warning(
    m <- suppressMessages(fit_mod(data_list = d, inits = NULL, estimateMode = 3, msmMode = 0,
                                  random_rec = FALSE, config = cfg,
                                  fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE))),
    "`qFun` differs between the data's model_config and `config`")
  testthat::expect_true(all(is.na(m$data_list$linkage_table$re_index)))
  testthat::expect_equal(sum(names(m$obj$env$par) == "beta_linkage_re"), 0L)
})
