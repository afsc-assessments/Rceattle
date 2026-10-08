# retrospective() fits each peel's hindcast with the map stored on the fitted
# object (R/9-retro_and_jitter.R, `map <- object$map`), so the peels reproduce
# the parameterisation the model was ORIGINALLY fitted with. That is
# deliberate -- Mohn's rho is a comparison against that fit -- and it is not
# changed here.
#
# The cost is that a later fix to how a map is built never reaches a saved fit.
# Worse, within one peel the two passes disagree: the hindcast uses the stored
# map while the forecast-catch refit calls build_map() afresh, so a peel can be
# fitted and reported under different parameter counts.
#
# So the divergence is reported instead of silently tolerated. The rebuild for
# the comparison uses the ORIGINAL data and parameters, not the peel's -- a
# difference then means the map-building code changed since the fit was saved,
# not that the peel has fewer years. Measured on a fresh fit: 49 map blocks
# compared, 0 differ.
#
# The check is hoisted ABOVE run_one_peel's definition on purpose. A warning()
# raised inside a .parallel_lapply() worker is discarded, and the default
# `cores` is detectCores() - 6, so the parallel path is the one anyone
# actually runs -- a per-peel warning would be invisible exactly where it
# matters. inst/dev/TRAPS.md records the class.

retro_fit <- function() {
  suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = Rceattle::BS2017SS, inits = NULL, file = NULL,
    estimateMode = "Hindcast", msmMode = 0, niter = 3, random_rec = FALSE,
    fit_control = Rceattle::fit_control(verbose = 0, getsd = FALSE,
                                        newtonsteps = 0))))
}

drift_warnings <- function(expr) {
  w <- character()
  withCallingHandlers(
    try(expr, silent = TRUE),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    },
    message = function(m) invokeRestart("muffleMessage"))
  grep("stored map differs", w, value = TRUE)
}


testthat::test_that("a current fit's stored map matches a fresh build", {
  testthat::skip_on_cran()
  fit <- retro_fit()
  fresh <- suppressWarnings(suppressMessages(Rceattle::build_map(
    data_list = fit$data_list, params = fit$estimated_params,
    debug = FALSE, random_rec = fit$data_list$random_rec)))
  shared <- intersect(names(fit$map$mapList), names(fresh$mapList))
  differing <- shared[!vapply(shared, function(nm) identical(
    as.integer(fit$map$mapList[[nm]]),
    as.integer(fresh$mapList[[nm]])), logical(1))]

  testthat::expect_gt(length(shared), 40L)   # 49 at the time of writing
  testthat::expect_equal(differing, character(0))
  # So retrospective() says nothing on a fit built by the current code.
  testthat::expect_length(drift_warnings(
    Rceattle::retrospective(fit, peels = 1, cores = 1)), 0L)
})


testthat::test_that("a stored map that no longer matches is reported", {
  testthat::skip_on_cran()
  fit <- retro_fit()
  stale <- fit
  # Stand in for a fit saved by code that built this block differently.
  stale$map$mapList$rec_dev[1, 1] <- NA

  hit <- drift_warnings(Rceattle::retrospective(stale, peels = 1, cores = 1))
  testthat::expect_length(hit, 1L)
  testthat::expect_match(hit, "rec_dev", fixed = TRUE)
  # It must say which map the peels actually used, since that is what Mohn's
  # rho is computed on.
  testthat::expect_match(hit, "STORED map", fixed = TRUE)
})


testthat::test_that("the warning survives the parallel peel path", {
  # The point of hoisting it. With cores > 1 the peels run through
  # .parallel_lapply(), which discards a worker's warnings -- so a check
  # placed inside run_one_peel() would be silent on the path the default
  # `cores` takes.
  testthat::skip_on_cran()
  stale <- retro_fit()
  stale$map$mapList$rec_dev[1, 1] <- NA

  for (cr in c(1L, 2L)) {
    # expect_warning() rather than the helper above: testthat installs its own
    # calling handlers, so a nested withCallingHandlers() does not reliably see
    # every warning inside a test_that() block.
    testthat::expect_warning(
      suppressMessages(Rceattle::retrospective(stale, peels = 2, cores = cr)),
      "stored map differs")
  }
})
