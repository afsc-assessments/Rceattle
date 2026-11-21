# The C++ and RTMB objects for one model, built from the pieces fit_mod() hands the C++
# MakeADFun(); shared by the test-rtmb-oracle-*.R files.
rtmb_pair <- function(data = BS2017SS, inits = NULL, msmMode = 0, suitMode = 0, niter = 3) {
  fit <- suppressMessages(suppressWarnings(fit_mod(
    data, inits = inits, file = NULL, msmMode = msmMode, estimateMode = 3, niter = niter,
    suitMode = rep(suitMode, data$nspp), random_rec = FALSE,
    fit_control = fit_control(getsd = FALSE, verbose = 0))))
  d <- fit$obj$env$data
  params <- fit$estimated_params[names(fit$obj$env$parameters)]
  list(cpp = fit$obj, data = d, params = params,
       rtmb = .rtmb_make_adfun(d, params, map = fit$map$mapFactor))
}

# fn, gr and jnll_comp under the gate's rules, and the tape against the report pass
expect_rtmb_matches <- function(o) {
  x <- o$cpp$par
  expect_identical(names(o$rtmb$par), names(x))
  fc <- o$cpp$fn(x); fr <- o$rtmb$fn(x)
  expect_lte(abs(fc - fr), 1e-12 * abs(fc))
  gc <- o$cpp$gr(x)
  expect_lte(max(abs(gc - o$rtmb$gr(x))), 1e-12 * max(abs(fc), max(abs(gc))))
  jc <- o$cpp$report(x)$jnll_comp; jr <- o$rtmb$report(x)$jnll_comp
  expect_identical(dim(jr), dim(jc))
  expect_true(all(abs(jc - jr) <= 1e-12 * abs(jc) + 1e-15 * max(abs(jc))))
  # The tape (fn) and the plain-R report pass compute the same objective
  expect_lte(abs(fr - sum(jr)), 1e-12 * abs(fr))
}
