# Convergence diagnostics (R/0-convergence.R). Exercised on minimal synthetic
# fit-like objects so the suite stays fast (no TMB fits).

make_fake_fit <- function(max_gradient = 1e-5, worst = "log_F", pdHess = TRUE) {
  structure(list(
    sdrep = NULL,
    .conv_hindcast = list(
      max_gradient = max_gradient,
      worst = list(param = worst, gradient = max_gradient),
      pdHess = pdHess)
  ), class = "Rceattle")
}

test_that("high gradient and non-PD Hessian are flagged", {
  fit <- make_fake_fit(max_gradient = 4e12, worst = "sel_inf", pdHess = FALSE)
  cv <- convergence_diagnostics(fit)
  expect_s3_class(cv, "Rceattle_convergence")
  expect_equal(cv$status, "FAIL")
  expect_equal(cv$checks$max_gradient$severity, "FAIL")
  expect_equal(cv$checks$pdHess$severity, "FAIL")
  expect_match(cv$checks$max_gradient$message, "sel_inf")
})

test_that("parameters print with the quantity they estimate", {
  expect_equal(.rce_par_display(c("log_M1", "rec_dev", "not_a_block")),
               c("log_M1 (M1)", "rec_dev (recruitment deviations)",
                 "not_a_block"))
})

# A hindcast snapshot whose gradient and index agree, so the checks can say
# where a parameter sits. Three parameters: M1 at ages 1-2 and one F.
.fake_located_fit <- function(gradient, cov = NULL) {
  nm  <- c("log_M1", "log_M1", "log_F")
  idx <- data.frame(par_index = 1:3, block = nm,
                    label = c("age 1", "age 2", "1990"),
                    stringsAsFactors = FALSE)
  fit <- make_fake_fit(max_gradient = max(abs(gradient)))
  fit$.conv_hindcast$gradient <- stats::setNames(gradient, nm)
  fit$.conv_hindcast$index    <- idx
  if (!is.null(cov)) {
    dimnames(cov) <- list(nm, nm)
    fit$sdrep <- list(cov.fixed = cov, pdHess = TRUE)
  }
  fit
}

test_that("the largest gradient is named by quantity and coordinate", {
  fit <- .fake_located_fit(c(1e-4, 2.2e-3, -5e-4))
  mg  <- convergence_diagnostics(fit)$checks$max_gradient
  expect_equal(mg$severity, "WARN")
  expect_match(mg$message, "on log_M1 (M1): age 2.", fixed = TRUE)
})

test_that("the distance to the optimum is reported in standard errors", {
  # The Newton step is cov %*% gradient = (0.04, 0.002, 0); over SEs (2, 1, 1)
  # that is (0.02, 0.002, 0), largest on the first parameter.
  fit <- .fake_located_fit(c(0.01, 0.002, 0), cov = diag(c(4, 1, 1)))
  mg  <- convergence_diagnostics(fit)$checks$max_gradient
  expect_equal(mg$data$newton_step_se, 0.02)
  expect_match(mg$message,
               "at most 0.02 standard errors (log_M1 (M1): age 1)", fixed = TRUE)
})

test_that("no step is reported when the sdreport is for another parameter vector", {
  # Under an estimating HCR the sdreport holds the projection's parameters.
  fit <- .fake_located_fit(c(0.01, 0.002, 0))
  fit$sdrep <- list(cov.fixed = matrix(1, 1, 1,
                                       dimnames = list("log_Ftarget", "log_Ftarget")),
                    pdHess = TRUE)
  mg <- convergence_diagnostics(fit)$checks$max_gradient
  expect_null(mg$data$newton_step_se)
  expect_false(grepl("standard errors", mg$message))
})

test_that("a scattered year set is counted against its span", {
  yrs <- c(1980:1985, 1990, 1995:1998, 2021)            # 12 years, not a run
  idx <- data.frame(par_index = seq_along(yrs), block = "log_F",
                    species = NA, fleet = "Hake_fishery", sex = NA, age = NA,
                    bin = NA, year = as.character(yrs), slot = NA,
                    stringsAsFactors = FALSE)
  out <- .rce_par_summary(idx$par_index, idx)
  expect_match(out, "log_F (F)", fixed = TRUE)
  expect_match(out, "12 years in 1980-2021  (12)", fixed = TRUE)
})

test_that("a converged fit is OK", {
  fit <- make_fake_fit(max_gradient = 1e-5, pdHess = TRUE)
  cv <- convergence_diagnostics(fit)
  expect_equal(cv$status, "OK")
})

# The conditioning check reads severity on the CORRELATION matrix, so the three
# tests below separate the two things a covariance condition number confounds:
# near-linear dependence between estimates (flagged) and a difference in the
# units the parameters are held in (not flagged). Provenance: the fixture here
# was previously diag(c(1, 1e8)) with a 1e3 covariance -- correlation 0.1, a
# scale difference and nothing more -- which the covariance reading called 1e8.
.corr_cov <- function(nm, rho, se) {
  p <- length(nm)
  cr <- matrix(rho, p, p); diag(cr) <- 1
  cov <- outer(se, se) * cr
  dimnames(cov) <- list(nm, nm)
  cov
}

test_that("Hessian eigen check flags a near-collinear parameter pair", {
  # Correlation 1 - 1e-8: a ridge. The standard errors differ by 1e4 to show the
  # verdict comes from the confounding and not from the units.
  cov <- .corr_cov(c("a", "b"), 1 - 1e-8, c(1, 1e4))
  fit <- make_fake_fit()
  fit$sdrep <- list(cov.fixed = cov, pdHess = TRUE)
  hc <- convergence_diagnostics(fit)$checks$hessian_conditioning
  expect_true(hc$severity %in% c("WARN", "FAIL"))
  expect_gt(hc$data$condition_number, 1e6)
  # se_ratio is the readable form of the same number.
  expect_equal(hc$data$se_ratio, sqrt(hc$data$condition_number))
  expect_match(hc$message, "standard error")
})

test_that("a mis-scaled but uncorrelated covariance is not flagged", {
  # diag(c(1, 1e8)) with a 1e3 covariance is correlation 0.1. Its covariance
  # condition number is 1e8; the estimates are not confounded, so the verdict is
  # OK and the covariance number is still recorded.
  cov <- diag(c(1, 1e8))
  cov[1, 2] <- cov[2, 1] <- 1e3
  dimnames(cov) <- list(c("a", "b"), c("a", "b"))
  fit <- make_fake_fit()
  fit$sdrep <- list(cov.fixed = cov, pdHess = TRUE)
  hc <- convergence_diagnostics(fit)$checks$hessian_conditioning
  expect_equal(hc$severity, "OK")
  expect_lt(hc$data$condition_number, 10)
  expect_gt(hc$data$covariance_condition_number, 1e6)
})

test_that("the conditioning verdict does not move when a parameter is rescaled", {
  # The point of reading the correlation: log_F sits near -2 and sel_inf near 10,
  # so a change of units must not change what the diagnostic says.
  cov <- .corr_cov(c("a", "b", "c"), 0.99, c(1, 1, 1))
  fit <- make_fake_fit(); fit$sdrep <- list(cov.fixed = cov, pdHess = TRUE)
  base <- convergence_diagnostics(fit)$checks$hessian_conditioning

  s <- c(1, 1e6, 1e-3)                       # rescale each parameter
  scaled <- diag(s) %*% cov %*% diag(s)
  dimnames(scaled) <- dimnames(cov)
  fit$sdrep <- list(cov.fixed = scaled, pdHess = TRUE)
  resc <- convergence_diagnostics(fit)$checks$hessian_conditioning

  expect_equal(resc$data$condition_number, base$data$condition_number,
               tolerance = 1e-6)
  expect_equal(resc$severity, base$severity)
  # The covariance reading is the one that moves, which is why it is not the
  # quantity severity is read on.
  expect_gt(resc$data$covariance_condition_number,
            base$data$covariance_condition_number * 1e3)
})

test_that("Hessian eigen check names a diffusely-loaded flat direction", {
  # Flat direction spread evenly over 40 rec_dev coefficients (each loading
  # ~0.16, below any single-coefficient threshold) plus a little ln_srv_sel.
  # The pre-fix check printed "loads on: ." with nothing after it.
  nm  <- c(rep("rec_dev", 40), rep("ln_srv_sel", 4))
  p   <- length(nm)
  v1  <- c(rep(1, 40), rep(1.5, 4)); v1 <- v1 / sqrt(sum(v1^2))
  Q   <- qr.Q(qr(cbind(v1, diag(p)[, -1])))
  cov <- Q %*% diag(c(1e8, rep(1, p - 1))) %*% t(Q)
  dimnames(cov) <- list(nm, nm)
  fit <- make_fake_fit(); fit$sdrep <- list(cov.fixed = cov, pdHess = TRUE)
  hc  <- convergence_diagnostics(fit)$checks$hessian_conditioning
  expect_true(hc$severity %in% c("WARN", "FAIL"))
  expect_match(hc$message, "loads on: [A-Za-z]")   # never blank after "loads on: "
  expect_match(hc$message, "rec_dev")              # the dominant block is named
  expect_true("rec_dev" %in% hc$data$loadings$param)
})

test_that("Hessian eigen check falls back to par.fixed names without dimnames", {
  nm  <- c("a", "b", "c")
  # 'a' and 'b' near-collinear, 'c' free: the flat direction is a against b.
  cov <- diag(3); cov[1, 2] <- cov[2, 1] <- 1 - 1e-9
  fit <- make_fake_fit()
  fit$sdrep <- list(cov.fixed = cov, pdHess = TRUE,
                    par.fixed = stats::setNames(c(0, 0, 0), nm))
  hc  <- convergence_diagnostics(fit)$checks$hessian_conditioning
  expect_match(hc$message, "loads on: [ab]: ")     # named from par.fixed, ...
  expect_false(grepl("p1", hc$message))            # ... not the "p1" placeholder
})

test_that("Hessian eigen check is OK on a well-conditioned covariance", {
  cov <- diag(c(1, 2)); dimnames(cov) <- list(c("a", "b"), c("a", "b"))
  fit <- make_fake_fit(); fit$sdrep <- list(cov.fixed = cov, pdHess = TRUE)
  cv <- convergence_diagnostics(fit)
  expect_equal(cv$checks$hessian_conditioning$severity, "OK")
  expect_equal(cv$status, "OK")
})

test_that("sdreport failure is a FAIL", {
  fit <- structure(list(.conv_hindcast = list(
    sd_requested = TRUE, sd_present = FALSE)), class = "Rceattle")
  expect_equal(convergence_diagnostics(fit)$checks$sdreport_failed$severity,
               "FAIL")
  fit2 <- structure(list(.conv_hindcast = list(
    sd_requested = FALSE, sd_present = FALSE)), class = "Rceattle")
  expect_false("sdreport_failed" %in%
                 names(convergence_diagnostics(fit2)$checks))
})

test_that("a parameter at a configured bound is flagged (WARN)", {
  fit <- structure(list(.conv_hindcast = list(
    par   = c(a = 0.5, b = 2.0),
    lower = c(0, -1), upper = c(1, 2))), class = "Rceattle")
  cv <- convergence_diagnostics(fit)
  expect_equal(cv$checks$parameters_on_bounds$severity, "WARN")
  expect_match(cv$checks$parameters_on_bounds$message, "b")
  fit2 <- structure(list(.conv_hindcast = list(
    par = c(log_F = -999), lower = c(-999), upper = c(10))), class = "Rceattle")
  expect_false("parameters_on_bounds" %in%
                 names(convergence_diagnostics(fit2)$checks))
})

test_that(".capture_opt_convergence aligns each MLE with its own bounds", {
  # Regression for the bug where the (mle, lower, upper) triple was assembled
  # from two independently-ordered sources, pairing one parameter's MLE with
  # another's bounds: rec_pars (unbounded, MLE 14.9) was reported "at" log_F's
  # [-1000, 10] upper bound. The capture must read the MLE back through the SAME
  # parameter-list shape as the bounds so rows can't drift apart.
  fake_obj <- list(
    env = list(
      last.par.best = c(rec_pars = 14.9, log_F = -999, catch_log_sd = 3),
      random = integer(0),
      # parList(x = par[lfixed()], par = last.par) expands the reduced vector
      # into the full param-list shape that `bounds` and `mapFactor` share.
      parList = function(x = NULL, par = NULL) list(rec_pars = 14.9,
                                                    log_F = -999,
                                                    catch_log_sd = 3)
    ),
    gr = function(p) rep(0, length(p))
  )
  bounds <- list(
    lower = list(rec_pars = -Inf, log_F = -1000, catch_log_sd = -10),
    upper = list(rec_pars =  Inf, log_F =    10, catch_log_sd =   3)
  )
  mapFactor <- list(rec_pars = 1, log_F = 1, catch_log_sd = 1)

  snap <- .capture_opt_convergence(opt = list(), obj = fake_obj,
                                   bounds = bounds, mapFactor = mapFactor,
                                   random_vars = NULL, getsd = FALSE)
  # bounds line up with their own parameter
  expect_equal(unname(snap$par["rec_pars"]), 14.9)
  expect_equal(snap$lower[match("rec_pars", names(snap$par))], -Inf)
  expect_equal(snap$upper[match("catch_log_sd", names(snap$par))], 3)

  fit <- structure(list(.conv_hindcast = snap), class = "Rceattle")
  cv  <- convergence_diagnostics(fit)
  # rec_pars is unbounded -> must NOT be flagged; log_F is a -999 sentinel ->
  # skipped; only catch_log_sd genuinely sits at its upper bound.
  ob <- cv$checks$parameters_on_bounds
  expect_equal(ob$severity, "WARN")
  expect_match(ob$message, "catch_log_sd")
  expect_false(grepl("rec_pars", ob$message))
  expect_false(grepl("log_F", ob$message))
})

test_that("a phase that ends with a high gradient is flagged (WARN)", {
  fit <- structure(list(.conv_phase = list(
    list(phase = 1, max_grad = 1e-4),
    list(phase = 2, max_grad = 5e3))), class = "Rceattle")
  cv <- convergence_diagnostics(fit)
  expect_equal(cv$checks$phasing$severity, "WARN")
  expect_match(cv$checks$phasing$message, "phase 2")
})

test_that("print method runs and is non-erroring", {
  fit <- make_fake_fit(max_gradient = 4e12, worst = "sel_inf", pdHess = FALSE)
  cv <- convergence_diagnostics(fit)
  expect_output(print(cv), "status: FAIL")
  # Captured, not printed. expect_invisible() does not sink output, so this
  # line used to write the fixture's "[FAIL] max_gradient ..." to stdout -- and
  # when a Windows testthat worker dies, R CMD check dumps whatever that file
  # had printed, which has been read as a real convergence regression three
  # times. The fixture is deliberately non-converged; it should stay quiet.
  invisible(capture.output(expect_invisible(print(cv))))
})

# getsd = FALSE leaves sdrep NULL, so the Hessian eigenvalue, sdreport, pdHess
# and estimability checks all return nothing and the battery used to report
# "OK" -- which report_tables() prints into a SAFE table as `converged`. A NOTE
# separates "every check passed" from "the strongest checks never ran".
test_that("getsd = FALSE is reported rather than passing silently", {
  fit <- make_fake_fit(max_gradient = 1e-5, pdHess = TRUE)
  fit$.conv_hindcast$sd_requested <- FALSE
  cv <- convergence_diagnostics(fit)
  expect_equal(cv$checks$hessian_not_run$severity, "NOTE")
  expect_match(cv$checks$hessian_not_run$message, "getsd = FALSE", fixed = TRUE)
  expect_equal(cv$status, "NOTE")

  # With an sdreport requested the record is absent, so a real battery is unchanged.
  fit$.conv_hindcast$sd_requested <- TRUE
  expect_null(convergence_diagnostics(fit)$checks$hessian_not_run)

  # A fit object that never recorded the flag says nothing either way.
  fit$.conv_hindcast$sd_requested <- NULL
  expect_null(convergence_diagnostics(fit)$checks$hessian_not_run)
  expect_equal(convergence_diagnostics(fit)$status, "OK")
})

# The two claims the Newton-step report rests on, neither of which had a test:
# that it changes no severity, and that it is withheld when the covariance
# cannot support it. Added reviewing PR #161 for the 5.43.0 release.
test_that("the Newton step reports without changing any severity", {
  g <- c(0.01, 0.002, 0)
  with_cov    <- convergence_diagnostics(.fake_located_fit(g, cov = diag(c(4, 1, 1))))
  without_cov <- convergence_diagnostics(.fake_located_fit(g))

  # The step is reported in one and absent from the other ...
  expect_false(is.null(with_cov$checks$max_gradient$data$newton_step_se))
  expect_null(without_cov$checks$max_gradient$data$newton_step_se)
  # ... and the severity is the same either way: it is read on the gradient.
  expect_identical(with_cov$checks$max_gradient$severity,
                   without_cov$checks$max_gradient$severity)
})

test_that("the Newton step is withheld when the covariance cannot carry it", {
  cov <- diag(c(4, 1, 1))

  # On a positive-definite Hessian it is reported.
  pd <- convergence_diagnostics(.fake_located_fit(c(0.01, 0.002, 0), cov = cov))
  expect_equal(pd$checks$max_gradient$data$newton_step_se, 0.02)

  # Not on a saddle: cov %*% gradient is the distance to a minimum only if cov
  # inverts one, and the sentence would read as reassurance beside the pdHess
  # check that just failed.
  saddle <- .fake_located_fit(c(0.01, 0.002, 0), cov = cov)
  saddle$sdrep$pdHess <- FALSE
  mg <- convergence_diagnostics(saddle)$checks$max_gradient
  expect_null(mg$data$newton_step_se)
  expect_false(grepl("standard errors", mg$message, fixed = TRUE))

  # A negative variance is dropped rather than square-rooted: this battery
  # reports through message() and must not emit "NaNs produced" from sqrt().
  neg <- cov; neg[3, 3] <- -1
  expect_no_warning(convergence_diagnostics(.fake_located_fit(c(0.01, 0.002, 0), cov = neg)))
})
