#' Prior distributions for Rceattle linkage coefficients
#'
#' Each prior is captured as a small object of class
#' `"Rceattle_prior"` holding a `family` name and two positional
#' parameters (`p1`, `p2`). The two-parameter shape is enforced so the
#' linkage table can store priors uniformly as four columns
#' (`prior_family`, `prior_p1`, `prior_p2`, plus reserved future use).
#'
#' Two surfaces are provided:
#'
#'   1. **Programmatic constructors**: exported with the `prior_`
#'      prefix (`prior_normal()`, `prior_lognormal()`, ...), safe to
#'      call anywhere without masking [base::gamma()]/[base::beta()].
#'   2. **Inline form** inside `linkage_spec(priors = ...)`. There the
#'      argument is captured unevaluated and resolved with a private
#'      data mask that makes `normal()`, `lognormal()`, `gamma()`, and
#'      `beta()` shorthand for the respective `prior_*` constructors,
#'      *only* inside that argument. Base R remains untouched at the
#'      package namespace.
#'
#' @name Rceattle_priors
#' @keywords internal
NULL


#' Allowed prior families for linkage coefficients
#'
#' @keywords internal
PRIOR_FAMILIES <- c("none", "normal", "lognormal", "gamma", "beta")


#' Internal: build a prior object after parameter validation.
#' @keywords internal
#' @noRd
new_prior <- function(family, p1, p2) {
  structure(
    list(family = family, p1 = as.numeric(p1), p2 = as.numeric(p2)),
    class = "Rceattle_prior"
  )
}


#' @keywords internal
#' @noRd
.assert_finite_scalar <- function(x, what) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    stop(sprintf("`%s` must be a finite numeric scalar", what), call. = FALSE)
  }
}


#' Normal prior on a linkage coefficient
#'
#' @param mean prior mean.
#' @param sd prior standard deviation (must be positive).
#' @return An `Rceattle_prior` of family `"normal"`.
#' @export
prior_normal <- function(mean, sd) {
  .assert_finite_scalar(mean, "mean")
  .assert_finite_scalar(sd, "sd")
  if (sd <= 0) stop("normal prior `sd` must be > 0", call. = FALSE)
  new_prior("normal", mean, sd)
}


#' Lognormal prior on a linkage coefficient
#'
#' A normal density on the log of the coefficient. Under
#' `fit_control(bias_adjust_proc = TRUE)`, the default, it is centred at
#' `meanlog - sdlog^2/2`, so `exp(meanlog)` is the prior mean; with `FALSE` it
#' matches [stats::dlnorm()] and `exp(meanlog)` is the prior median.
#'
#' @param meanlog log of the prior mean (of the median when `bias_adjust_proc = FALSE`).
#' @param sdlog prior standard deviation of the log (must be positive).
#' @return An `Rceattle_prior` of family `"lognormal"`.
#' @export
prior_lognormal <- function(meanlog, sdlog) {
  .assert_finite_scalar(meanlog, "meanlog")
  .assert_finite_scalar(sdlog, "sdlog")
  if (sdlog <= 0) stop("lognormal prior `sdlog` must be > 0", call. = FALSE)
  new_prior("lognormal", meanlog, sdlog)
}


#' Gamma prior on a linkage coefficient
#'
#' Shape-rate parameterization, matching [stats::dgamma()].
#'
#' @param shape positive shape parameter.
#' @param rate positive rate parameter.
#' @return An `Rceattle_prior` of family `"gamma"`.
#' @export
prior_gamma <- function(shape, rate) {
  .assert_finite_scalar(shape, "shape")
  .assert_finite_scalar(rate, "rate")
  if (shape <= 0 || rate <= 0) {
    stop("gamma prior `shape` and `rate` must both be > 0", call. = FALSE)
  }
  new_prior("gamma", shape, rate)
}


#' Beta prior on a linkage coefficient
#'
#' Standard Beta(shape1, shape2) on (0, 1), matching [stats::dbeta()].
#' For priors on stock-recruit steepness on the standard (0.2, 1)
#' interval, transform inside the model and use [prior_beta()] on the
#' rescaled quantity.
#'
#' @param shape1,shape2 positive shape parameters.
#' @return An `Rceattle_prior` of family `"beta"`.
#' @export
prior_beta <- function(shape1, shape2) {
  .assert_finite_scalar(shape1, "shape1")
  .assert_finite_scalar(shape2, "shape2")
  if (shape1 <= 0 || shape2 <= 0) {
    stop("beta prior `shape1` and `shape2` must both be > 0", call. = FALSE)
  }
  new_prior("beta", shape1, shape2)
}


#' Test whether an object is an Rceattle prior
#' @param x object to test.
#' @keywords internal
is_rceattle_prior <- function(x) inherits(x, "Rceattle_prior")


#' @export
print.Rceattle_prior <- function(x, ...) {
  cat(sprintf("<Rceattle prior: %s(%g, %g)>\n", x$family, x$p1, x$p2))
  invisible(x)
}


#' Internal: data mask exposing short prior names to NSE-aware callers
#'
#' Returned as a named list rather than an environment so it can be
#' passed straight to [rlang::eval_tidy()] as a `data` mask. Names are
#' bound to the corresponding `prior_*` constructors so user code such
#' as `priors = list(temp = normal(0, 1))` resolves correctly without
#' masking [base::gamma()] or [base::beta()] at the package level.
#'
#' @keywords internal
.prior_dispatch_mask <- function() {
  list(
    normal    = prior_normal,
    lognormal = prior_lognormal,
    gamma     = prior_gamma,
    beta      = prior_beta
  )
}
