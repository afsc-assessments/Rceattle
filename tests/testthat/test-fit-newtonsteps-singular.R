# TMBhelper::fit_tmb() takes Newton steps with solve(H, g), unguarded, so a Hessian
# singular in one direction ("system is exactly singular") used to lose the whole fit
# from fit_mod(). Seen in the NOAA model comparison on a state-space replicate.

# An objective informed in `a` only: its Hessian is exactly singular in `b`.
.flat_obj <- function() {
  env <- new.env()
  env$parameters <- list(a = 0, b = 0)
  list(par = c(a = 0, b = 0),
       fn  = function(p) (p[[1]] - 1)^2,
       gr  = function(p) c(2 * (p[[1]] - 1), 0),
       env = env)
}

testthat::test_that("a singular Hessian stops Newton steps instead of the fit", {
  testthat::skip_if_not_installed("TMBhelper")
  obj <- .flat_obj()
  # Positive control: TMBhelper alone errors here.
  testthat::expect_error(
    suppressWarnings(TMBhelper::fit_tmb(obj, startpar = obj$par, newtonsteps = 1,
                                        getsd = FALSE, quiet = TRUE,
                                        control = list(trace = 0))),
    "singular")
  testthat::expect_warning(
    opt <- Rceattle:::.fit_tmb(obj, newtonsteps = 3, getsd = FALSE),
    "singular Hessian")
  testthat::expect_equal(unname(opt$par[["a"]]), 1, tolerance = 1e-6)
  # The refit is the fit newtonsteps = 0 gives.
  ref <- Rceattle:::.fit_tmb(.flat_obj(), newtonsteps = 0, getsd = FALSE)
  testthat::expect_identical(opt$par, ref$par)
  testthat::expect_identical(opt$objective, ref$objective)
})

testthat::test_that("only an error from solve() is retried", {
  testthat::skip_if_not_installed("TMBhelper")
  obj <- .flat_obj()
  obj$fn <- function(p) stop("a singular value, not from solve()")
  testthat::expect_error(Rceattle:::.fit_tmb(obj, newtonsteps = 3, getsd = FALSE),
                         "not from solve")
})

testthat::test_that("without Newton steps the same object fits with no warning", {
  testthat::skip_if_not_installed("TMBhelper")
  obj <- .flat_obj()
  testthat::expect_no_warning(opt <- Rceattle:::.fit_tmb(obj, newtonsteps = 0, getsd = FALSE))
  testthat::expect_equal(unname(opt$par[["a"]]), 1, tolerance = 1e-6)
})
