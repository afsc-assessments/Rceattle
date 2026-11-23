# TMBhelper::fit_tmb() keeps every Newton step, even one that raises the
# objective, and reports the last, while fit_mod() reads its estimates at TMB's
# best point (last.par.best). On a NOAA model-comparison replicate (C0 sim 211,
# rel.tol = 1e-4) that reported max_gradient 9.2e-5 for estimates whose own
# gradient was 0.131.

# (x^2 - 1)^2 + y^2: minima at x = +/-1, a local maximum in x at x = 0. From
# (0.3, 0.5), where the curvature in x is negative, Newton step 1 lowers the
# objective (1.078 to 0.989) and steps 2-3 climb to the maximum (1.0).
# iter.max = 0 stands in for an nlminb stop there. fn() keeps the best point
# the way TMB::MakeADFun() does for a model without random effects.
.saddle_obj <- function(start = c(x = 0.3, y = 0.5)) {
  env <- new.env()
  env$parameters <- list(x = 0, y = 0)
  env$last.par.best <- start
  env$value.best <- Inf
  fn <- function(p) {
    v <- (p[[1]]^2 - 1)^2 + p[[2]]^2
    if (is.finite(v) && v < env$value.best) {
      env$last.par.best <- stats::setNames(as.numeric(p), c("x", "y"))
      env$value.best <- v
    }
    v
  }
  gr <- function(p) c(4 * p[[1]] * (p[[1]]^2 - 1), 2 * p[[2]])
  list(par = start, fn = fn, gr = gr, env = env)
}
.early_stop <- list(iter.max = 0, eval.max = 1e9, trace = 0)
.fit_saddle <- function(obj, ..., control = .early_stop) {
  suppressMessages(Rceattle:::.fit_tmb(obj, newtonsteps = 3,
                                       control = control, ...))
}
.tmbhelper_saddle <- function(obj, control = .early_stop) {
  suppressMessages(TMBhelper::fit_tmb(
    obj, startpar = obj$par, loopnum = 1, newtonsteps = 3, getsd = FALSE,
    control = control, quiet = TRUE))
}

testthat::test_that("after an uphill Newton step the fit reports at the estimates", {
  testthat::skip_if_not_installed("TMBhelper")
  # Positive control: TMBhelper alone reports a near-zero gradient at the local
  # maximum, a point that is not TMB's best.
  raw_obj <- .saddle_obj()
  raw <- .tmbhelper_saddle(raw_obj)
  testthat::expect_lt(raw$max_gradient, 1e-4)
  testthat::expect_gt(raw$objective, raw_obj$env$value.best)

  obj <- .saddle_obj()
  opt <- .fit_saddle(obj, getsd = FALSE)
  best <- raw_obj$env$last.par.best
  testthat::expect_lt(best[["x"]], -0.05)        # step 1, not the maximum at 0
  testthat::expect_identical(unname(opt$par), unname(best))
  testthat::expect_identical(opt$objective, obj$env$value.best)
  testthat::expect_equal(opt$max_gradient, max(abs(obj$gr(best))))
  testthat::expect_gt(opt$max_gradient, 0.1)
  testthat::expect_equal(opt$diagnostics$final_gradient, obj$gr(best))
  testthat::expect_identical(opt$diagnostics$starting_value, c(0.3, 0.5))
  testthat::expect_identical(opt$Convergence_check,
                             "The model is likely not converged")
})

testthat::test_that("the re-evaluation at the estimates does not optimize again", {
  testthat::skip_if_not_installed("TMBhelper")
  # One nlminb iteration first, so only .fit_tmb()'s own iter.max = 0 stops the
  # re-evaluation from moving off TMBhelper's best point.
  ctl <- list(iter.max = 1, eval.max = 1e9, trace = 0)
  start <- c(x = 0.05, y = 0.5)
  raw_obj <- .saddle_obj(start)
  raw <- .tmbhelper_saddle(raw_obj, control = ctl)
  opt <- .fit_saddle(.saddle_obj(start), getsd = FALSE, control = ctl)
  testthat::expect_identical(unname(opt$par), unname(raw_obj$env$last.par.best))
  testthat::expect_identical(opt$iterations, raw$iterations)
  testthat::expect_identical(opt$evaluations, raw$evaluations)
})

testthat::test_that("a best point past a bound is reported there, not clamped", {
  testthat::skip_if_not_installed("TMBhelper")
  # Newton steps ignore bounds; step 1 lands at x = -0.074, below the bound at 0.
  obj <- .saddle_obj()
  opt <- .fit_saddle(obj, lower = c(0, -Inf), upper = c(Inf, Inf),
                     getsd = FALSE)
  testthat::expect_lt(opt$par[["x"]], 0)
  testthat::expect_identical(unname(opt$par), unname(obj$env$last.par.best))
  testthat::expect_identical(opt$diagnostics$Lower, c(0, -Inf))
})

testthat::test_that("getsd = TRUE judges the Hessian at the estimates, once", {
  testthat::skip_if_not_installed("TMBhelper")
  # The curvature in x is negative at step 1, so TMBhelper stops before
  # sdreport() and warns.
  obj <- .saddle_obj()
  w <- character(0)
  opt <- withCallingHandlers(
    .fit_saddle(obj, getsd = TRUE),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    })
  testthat::expect_length(grep("not positive definite", w), 1L)
  x <- obj$env$last.par.best[["x"]]
  testthat::expect_identical(unname(opt$par), unname(obj$env$last.par.best))
  # The Hessian TMBhelper judged is d2/dx2 = 12 x^2 - 4 at step 1 (-3.93), not
  # at the maximum (-4).
  testthat::expect_lt(abs(opt$hessian[1, 1] - (12 * x^2 - 4)), 1e-4)
  testthat::expect_identical(opt$Convergence_check,
                             "The model is definitely not converged")
})

testthat::test_that("downhill Newton steps report what TMBhelper reports", {
  testthat::skip_if_not_installed("TMBhelper")
  start <- c(x = 0.9, y = 0.1)
  raw <- .tmbhelper_saddle(.saddle_obj(start))
  opt <- .fit_saddle(.saddle_obj(start), getsd = FALSE)
  testthat::expect_identical(opt$par, raw$par)
  testthat::expect_identical(opt$objective, raw$objective)
  testthat::expect_identical(opt$max_gradient, raw$max_gradient)
})

testthat::test_that("the fallback without TMBhelper reports at the estimates", {
  # Run .fit_tmb() as if TMBhelper were not installed.
  fallback <- Rceattle:::.fit_tmb
  environment(fallback) <- list2env(
    list(requireNamespace = function(...) FALSE),
    parent = asNamespace("Rceattle"))
  obj <- .saddle_obj()
  opt <- fallback(obj, newtonsteps = 3, getsd = FALSE, control = .early_stop)
  best <- obj$env$last.par.best
  testthat::expect_lt(opt$objective, 0.99)       # step 1, not the maximum (1.0)
  testthat::expect_identical(unname(opt$par), unname(best))
  testthat::expect_equal(opt$max_gradient, max(abs(obj$gr(best))))
  # The fallback's own verdict, so this cannot pass on the TMBhelper path.
  testthat::expect_identical(opt$Convergence_check,
    "The maximum gradient is high; the model may not be converged")
})
