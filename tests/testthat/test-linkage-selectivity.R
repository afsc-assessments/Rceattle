# Environmental linkages on selectivity parameters (logistic family).

.sel_env_data <- function(d) {
  yrs <- d$styr:d$projyr
  set.seed(1)
  data.frame(Year = yrs, temp = as.numeric(scale(sin(seq_along(yrs) / 4))))
}

.logistic_fleet <- function(d) {
  fc <- Rceattle::switch_check(Rceattle::clean_data(d))$fleet_control
  which(fc$Selectivity == "Logistic")[1]
}

testthat::test_that("sel is a wired linkage process", {
  testthat::expect_true("sel" %in% Rceattle:::LINKAGE_PROCESSES_IMPLEMENTED)
})


testthat::test_that("build_selectivity() validates its linkages", {
  testthat::expect_null(Rceattle::build_selectivity()$linkages)
  ok <- Rceattle::build_selectivity(linkages = list(
    inf_asc = Rceattle::linkage_spec(~ temp, by = ~ fleet)))
  testthat::expect_s3_class(ok$linkages$inf_asc, "Rceattle_linkage_spec")
  testthat::expect_error(
    Rceattle::build_selectivity(linkages = list(
      nope = Rceattle::linkage_spec(~ temp))),
    "unknown")
})


testthat::test_that("the logistic slots map to the right parameter arrays", {
  # slope slots on log_sel_slp, inflection slots on sel_inf.
  m <- Rceattle:::.SEL_PARAM_TO_SLOT
  testthat::expect_equal(m$slp_asc$arr,  "log_sel_slp")
  testthat::expect_equal(m$inf_desc$arr, "sel_inf")
  testthat::expect_equal(m$inf_desc$slot, 2L)
  # DoubleNormal aliases share the slots.
  testthat::expect_equal(m$sigma_asc$slot, m$slp_asc$slot)
  testthat::expect_equal(m$peak$slot,      m$inf_asc$slot)
})


testthat::test_that("a model without a sel linkage is numerically unchanged", {
  testthat::skip_on_cran()
  fit <- Rceattle::fit_mod(
    data_list = Rceattle::BS2017SS, estimateMode = 1, msmMode = 0,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0))
  testthat::expect_equal(fit$opt$objective, 10241.0304274978,
                         tolerance = 1e-8)
})


testthat::test_that("a sel linkage on a logistic fleet is estimated and moves the fit", {
  testthat::skip_on_cran()
  d <- Rceattle::BS2017SS
  d$env_data <- .sel_env_data(d)
  flt <- .logistic_fleet(d)
  testthat::skip_if(is.na(flt))

  fit <- Rceattle::fit_mod(
    data_list = d, estimateMode = 1, msmMode = 0,
    selFun = Rceattle::build_selectivity(linkages = list(
      inf_asc = Rceattle::linkage_spec(~ temp, by = ~ fleet, fleet = flt))),
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0))

  tbl <- fit$data_list$linkage_table
  testthat::expect_equal(sum(tbl$process == "sel"), 2L)  # intercept + slope
  testthat::expect_lt(fit$opt$objective, 10241.0304274978)
})


testthat::test_that("a right_floor (log-link) linkage on a double-normal fits and stays bounded", {
  testthat::skip_on_cran()
  # An unphased double-normal fit of the whole GOA fleet set, and by far the most
  # expensive test in the suite under coverage: covr rebuilds the TMB model at
  # -O0, so this one block is 27 s here and ~26 minutes there -- 79% of a file
  # that was measured at 33 minutes and set the floor for the slowest shard. The
  # R it covers is the same linkage machinery the two fits above already run at a
  # tenth of the cost; what is unique to it is the numerical result, and that is
  # R-CMD-check's job.
  testthat::skip_on_covr()
  # right_floor is logit-parameterised (sel_inf(1) = logit(right_floor)), so a log-link
  # linkage offset must act on the logit, not multiply the probability (which could push
  # right_floor > 1 and corrupt the descending limb). Here the fit must stay finite and
  # selectivity a valid probability for every bin/year.
  data("GOA2018SS")
  d <- GOA2018SS
  rows_use <- which(d$fleet_control$Selectivity != 0 & d$fleet_control$Fleet_type != 0)
  d$fleet_control$Fleet_type[-rows_use] <- 0
  d$fleet_control$Selectivity        <- "DoubleNormal"
  d$fleet_control$Time_varying_sel   <- 0
  d$fleet_control$Selectivity_index  <- seq_len(nrow(d$fleet_control))
  d$fleet_control$Age_first_selected <- 1
  flt <- rows_use[1]
  yrs <- d$styr:d$projyr
  d$env_data <- data.frame(Year = yrs, temp = as.numeric(scale(seq_along(yrs))))

  fit <- Rceattle::fit_mod(
    data_list = d, estimateMode = 1, msmMode = 0,
    selFun = Rceattle::build_selectivity(linkages = list(
      right_floor = Rceattle::linkage_spec(~ temp, by = ~ fleet, fleet = flt, link = "log"))),
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))

  testthat::expect_equal(sum(fit$data_list$linkage_table$process == "sel"), 2L)  # intercept + slope
  testthat::expect_true(is.finite(fit$opt$objective))
  sel <- fit$quantities$sel_at_age
  testthat::expect_true(all(is.finite(sel)))
  testthat::expect_true(all(sel >= -1e-8 & sel <= 1 + 1e-8))
})


testthat::test_that("coff (non-parametric) is rejected as a normalization no-op", {
  testthat::skip_on_cran()
  d <- Rceattle::BS2017SS
  d$env_data <- .sel_env_data(d)

  # A per-year offset on mean-centred coefficients cancels exactly, so coff
  # is refused with an explanation rather than silently doing nothing.
  err <- tryCatch(
    Rceattle::fit_mod(
      data_list = d, estimateMode = 3, msmMode = 0,
      selFun = Rceattle::build_selectivity(linkages = list(
        coff = Rceattle::linkage_spec(~ temp, by = ~ fleet, fleet = 1))),
      fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                          verbose = 0)),
    error = function(e) conditionMessage(e))
  testthat::expect_match(err, "not supported")
  testthat::expect_match(err, "cancels exactly")
})


testthat::test_that("a sel linkage on a non-parametric fleet is rejected", {
  testthat::skip_on_cran()
  d <- Rceattle::BS2017SS
  d$env_data <- .sel_env_data(d)

  # Any fleet whose form is not one of the parametric wired forms.
  fc <- Rceattle::switch_check(Rceattle::clean_data(d))$fleet_control
  wired <- c("Logistic", "DoubleLogistic", "DescendingLogistic",
             "DoubleNormal", "LogisticPM")
  unwired <- which(!fc$Selectivity %in% wired & fc$Fleet_type != "Off")
  testthat::skip_if(length(unwired) == 0)
  testthat::expect_error(
    Rceattle::fit_mod(
      data_list = d, estimateMode = 3, msmMode = 0,
      selFun = Rceattle::build_selectivity(linkages = list(
        inf_asc = Rceattle::linkage_spec(~ temp, by = ~ fleet,
                                         fleet = unwired[1]))),
      fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                          verbose = 0)),
    "not yet")
})
