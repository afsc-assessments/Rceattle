# A deviation standard deviation that has been estimated to zero.
#
# Provenance. `Atka2022` under `random_sel = TRUE` with a non-parametric random
# walk returns an objective of 185.82 and `sel_dev_sd = 2.7e-08`: a
# time-invariant selectivity fitted by a model configured as time-varying. The
# battery flags that fit through `max_gradient`, which reports that the optimizer
# stopped rather than what went wrong; collapse can also occur at a CLEAN
# gradient, for data reasons -- few composition years, small sample sizes, or a
# composition weight competing with the process variance -- and nothing else
# catches that.
#
# What these tests do and do not establish. `fit_mod()` now refuses the
# configuration above, so the collapse cannot be produced from a supported
# setting and the third test below injects the value instead. That exercises the
# threshold and the message, NOT the wiring. The wiring was checked separately by
# running that fit on the pre-change package: `sel_dev_log_sd` reaches
# `.conv_hindcast$par` carrying 2.69e-08, so a collapsing fit does put a small
# value under a name this check reads. Every other name in `.CONV_PROCESS_SD` is
# wired by the same mechanism but has not been observed collapsing in a real fit;
# `R_log_sd` and `log_sigma_linkage` are confirmed present and healthy below.
testthat::skip_on_cran()

# The healthy and collapsed-sd tests read the same BS2017SS random-recruitment
# hindcast, which takes ~150 s, so it is fitted once on first use.
.healthy_fit <- local({
  m <- NULL
  function() {
    if (is.null(m)) {
      data("BS2017SS")
      m <<- suppressMessages(suppressWarnings(Rceattle::fit_mod(
        data_list = BS2017SS, inits = NULL, msmMode = 0,
        estimateMode = "Hindcast", random_rec = TRUE,
        fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0))))
    }
    m
  }
})

testthat::test_that("a healthy fit does not fire the check", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  m <- .healthy_fit()

  # The recruitment sds are O(1) here, which is what an assessment intends.
  p   <- m$.conv_hindcast$par
  sds <- exp(p[names(p) %in% Rceattle:::.CONV_PROCESS_SD])
  testthat::expect_gt(min(sds), 1e-3)
  testthat::expect_length(Rceattle:::.check_variance_collapse(m), 0)
})

testthat::test_that("an integrated linkage sd reaches the snapshot", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # The second of the two names confirmed to reach `$par` from a real fit. A
  # healthy value here is what the check must stay silent on.
  data("BS2017SS")
  d <- BS2017SS
  d$fleet_control$Catchability[7] <- "Estimated"
  q <- Rceattle::build_catchability(linkages = list(
         q = Rceattle::linkage_spec(~ (1 | Year), by = ~ fleet, fleet = 7)))
  m <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = d, inits = NULL, msmMode = 0, estimateMode = "Hindcast", qFun = q,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0))))

  testthat::expect_true("log_sigma_linkage" %in% names(m$.conv_hindcast$par))
  testthat::expect_gt(exp(m$.conv_hindcast$par[["log_sigma_linkage"]]), 1e-3)
  testthat::expect_length(Rceattle:::.check_variance_collapse(m), 0)
})

testthat::test_that("a collapsed sd is reported, and named", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # A local copy: the injection below does not reach the shared fit.
  m <- .healthy_fit()

  # The value Atka2022 reaches, injected so the assertion does not depend on a
  # configuration fit_mod() now refuses.
  m$.conv_hindcast$par["R_log_sd"] <- log(2.69e-8)
  r <- Rceattle:::.check_variance_collapse(m)
  testthat::expect_length(r, 1)
  # WARN, never FAIL: a variance at the boundary is a well-posed optimum and a
  # verdict about the model, not the optimizer, and report_tables() prints this
  # status as a fit's `converged` field.
  testthat::expect_equal(r[[1]]$severity, "WARN")
  testthat::expect_match(r[[1]]$message, "R_log_sd")
  testthat::expect_true("R_log_sd" %in% r[[1]]$data$param)

  # BS2017SS has three species, so the block is named once per element; the
  # message has to say WHICH one.
  testthat::expect_match(r[[1]]$message, "R_log_sd[1]", fixed = TRUE)
  testthat::expect_equal(r[[1]]$data$index, 1L)

  # And it reaches the assembled battery, not just the helper.
  testthat::expect_true("variance_collapse" %in%
    vapply(Rceattle::convergence_diagnostics(m)$checks,
           function(x) x$id, character(1)))
})

testthat::test_that("a fixed variance is a choice, not a collapse", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Only ESTIMATED parameters are in `$par`, so a deviation sd held at
  # Time_varying_sel_sd -- the ordinary penalized fit -- cannot fire this.
  data("Atka2022")
  m <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = Atka2022, inits = NULL, msmMode = 0, estimateMode = "Hindcast",
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0))))
  testthat::expect_false("sel_dev_log_sd" %in% names(m$.conv_hindcast$par))
  testthat::expect_length(Rceattle:::.check_variance_collapse(m), 0)
})
