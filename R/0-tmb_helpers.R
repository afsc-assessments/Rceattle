#' Fit a TMB object with optional TMBhelper enhancements
#'
#' Internal wrapper used by [fit_mod()]. Delegates to
#' \code{TMBhelper::fit_tmb} when the suggested package is installed,
#' and otherwise performs an equivalent core fit using \code{stats::nlminb}
#' followed by \code{TMB::sdreport}. The fallback path covers `loopnum`
#' restarts, optional `newtonsteps`, and the `getsd` /
#' `getJointPrecision` / `getReportCovariance` reporting options.
#'
#' @keywords internal
#' @noRd
.fit_tmb <- function(obj, lower = -Inf, upper = Inf,
                     control = list(eval.max = 1e9, iter.max = 1e9, trace = 0),
                     loopnum = 1, newtonsteps = 0,
                     getsd = TRUE, bias.correct = FALSE,
                     getJointPrecision = FALSE, getReportCovariance = FALSE,
                     quiet = TRUE) {

  if (requireNamespace("TMBhelper", quietly = TRUE)) {
    return(.normalize_fit_tmb(
      TMBhelper::fit_tmb(
        obj                 = obj,
        fn                  = obj$fn,
        gr                  = obj$gr,
        startpar            = obj$par,
        lower               = lower,
        upper               = upper,
        loopnum             = loopnum,
        newtonsteps         = newtonsteps,
        getsd               = getsd,
        control             = control,
        bias.correct        = bias.correct,
        bias.correct.control = list(sd = getsd),
        getJointPrecision   = getJointPrecision,
        getReportCovariance = getReportCovariance,
        quiet               = quiet
      )
    ))
  }

  # Fallback: plain nlminb + sdreport
  start_par <- obj$par
  best <- NULL
  for (i in seq_len(max(1, loopnum))) {
    init <- if (is.null(best)) start_par else best$par
    fit_i <- tryCatch(
      stats::nlminb(start = init, objective = obj$fn, gradient = obj$gr,
                    lower = lower, upper = upper, control = control),
      error = function(e) NULL
    )
    if (!is.null(fit_i) && (is.null(best) || fit_i$objective < best$objective)) {
      best <- fit_i
    }
  }
  if (is.null(best)) {
    stop("Optimization failed in .fit_tmb fallback (install TMBhelper for richer diagnostics).")
  }

  # Optional Newton steps to refine the gradient. These are UNCONSTRAINED, as
  # TMBhelper's own are, so a parameter nlminb parked on a bound can be pushed
  # past it; convergence_diagnostics() reports that as parameters_outside_bounds.
  if (newtonsteps > 0) {
    par <- best$par
    for (i in seq_len(newtonsteps)) {
      g <- obj$gr(par)
      h <- tryCatch(stats::optimHess(par, fn = obj$fn, gr = obj$gr),
                    error = function(e) NULL)
      if (is.null(h)) break
      step <- tryCatch(solve(h, as.numeric(g)), error = function(e) NULL)
      if (is.null(step)) break
      par <- par - step
    }
    best$par <- par
    best$objective <- obj$fn(par)
  }

  best$diagnostics <- data.frame(
    Param = names(best$par),
    starting_value = start_par,
    Lower = if (length(lower) == 1) rep(lower, length(start_par)) else lower,
    MLE = best$par,
    Upper = if (length(upper) == 1) rep(upper, length(start_par)) else upper,
    final_gradient = as.numeric(obj$gr(best$par))
  )
  best$max_gradient <- max(abs(best$diagnostics$final_gradient))
  best$Convergence_check <- if (best$max_gradient < 1e-4) {
    "There is no evidence that the model is not converged"
  } else {
    "The maximum gradient is high; the model may not be converged"
  }
  best$number_of_coefficients <- c(
    Total = length(start_par),
    Fixed_effects = length(start_par) - length(obj$env$random),
    Random_effects = length(obj$env$random)
  )

  if (isTRUE(getsd)) {
    best$SD <- tryCatch(
      TMB::sdreport(obj,
                    bias.correct        = bias.correct,
                    bias.correct.control = list(sd = getsd),
                    getJointPrecision   = getJointPrecision,
                    getReportCovariance = getReportCovariance),
      error = function(e) NULL
    )
  }

  best
}


#' Normalize `TMBhelper::fit_tmb()`'s non-positive-definite early return
#'
#' @description
#' With `getsd = TRUE`, `TMBhelper::fit_tmb()` factorizes the fixed-effect
#' Hessian before calling `sdreport`. If that `chol()` fails it warns and
#' returns `list(opt = <estimates>, h = <Hessian>)`, a *different shape* from
#' its documented return, with no `$objective`, `$max_gradient` or
#' `$Convergence_check` at the top level.
#'
#' Everything downstream expects the documented shape, so a non-positive-definite
#' fit otherwise arrives as a malformed `fit$opt`: `fit$opt$objective` is `NULL`,
#' and the convergence gate the diagnostic re-fitters share
#' (`.refit_converged()`) sees no verdict at all rather than the failure that
#' actually occurred.
#'
#' Unwrap it back to the estimates and record the verdict `fit_tmb()` itself
#' uses when `sdreport` returns `pdHess = FALSE`, what it would have said had
#' it not returned early. `$SD` stays absent, so `fit_mod()` still records
#' `sdrep = NULL` and `.check_sdreport_failed()` still fires.
#'
#' The Hessian is kept as `$hessian`, the name
#' `fit_tmb(getHessian = TRUE)` gives it. Nothing in the package reads it; it is
#' kept only so the unwrap does not discard something the malformed shape used
#' to expose (as `$opt$h`), since it is the one artifact that says *how* the
#' factorization failed.
#'
#' A no-op for every converged fit and for the whole `getsd = FALSE` path.
#'
#' Indexes with `[[` throughout, so the guard reads exactly the names it means
#', `$` partially matches on lists, and this discriminates between two shapes
#' by which names are present.
#'
#' @param x The value returned by [TMBhelper::fit_tmb()].
#' @return `x` in the documented `fit_tmb()` shape.
#' @keywords internal
#' @noRd
.normalize_fit_tmb <- function(x) {
  if (!is.list(x) || !is.null(x[["objective"]])) return(x)
  inner <- x[["opt"]]
  if (!is.list(inner) || is.null(inner[["objective"]])) return(x)
  inner[["hessian"]] <- x[["h"]]
  inner[["Convergence_check"]] <- "The model is definitely not converged"
  inner
}


#' Check parameter estimability (TMBhelper::check_estimability fallback)
#'
#' @keywords internal
#' @noRd
.check_estimability <- function(obj) {
  if (requireNamespace("TMBhelper", quietly = TRUE)) {
    return(TMBhelper::check_estimability(obj))
  }
  # Minimal in-package fallback: test which fixed effects move the gradient
  par <- obj$env$last.par.best
  fixed <- setdiff(seq_along(par), obj$env$random)
  g <- as.numeric(obj$gr(par[fixed]))
  bad <- which(abs(g) > 1)
  list(
    BadParams = data.frame(
      Param = names(par[fixed]),
      MLE = unname(par[fixed]),
      Param_check = ifelse(seq_along(fixed) %in% bad, 2L, 1L),
      Final_gradient = g
    ),
    WhichBad = bad
  )
}


#' Gap behind convergence warning (8)
#'
#' @description
#' The objective `nlminb` reported must survive a fresh evaluation of the object
#' it came from. Both sides are the MARGINAL objective. `obj$report()$jnll` is
#' the JOINT negative log-likelihood at the random-effect mode, a Laplace
#' correction away, so it cannot check a model holding random effects.
#'
#' Evaluated at `last.par.best`, the best point TMB saw. That is not always the
#' iterate `nlminb` returned, and stopping on a worse one is itself the
#' discontinuity this looks for. A fresh call rather than `obj$env$value.best`
#' because a `retape()` resets that to `Inf`.
#'
#' Only meaningful while `opt` came from `obj`. `fit_mod()` rebuilds `obj` for
#' the projection and again for `projection_uncertainty` without re-optimizing,
#' and those maps change which parameters are random. Measure beside each
#' optimization, never at the end.
#'
#' @param obj A TMB object from [TMB::MakeADFun()].
#' @param opt The `.fit_tmb()` result for that same object.
#' @return The absolute difference in negative log-likelihood units, or
#'   `NA_real_` when it cannot be formed.
#' @noRd
.rce_marginal_gap <- function(obj, opt) {
  if (is.null(obj) || is.null(opt) || is.null(opt$objective)) return(NA_real_)
  par <- obj$env$last.par.best
  if (is.null(par)) return(NA_real_)

  # -integer(0) would drop every parameter rather than none.
  fixed <- if (length(obj$env$random) > 0) par[-obj$env$random] else par

  # obj$fn() writes to the TMB environment, and that object is returned to the
  # caller, so the evaluation state is restored after the call.
  lpb <- obj$env$last.par.best
  lp  <- obj$env$last.par
  vb  <- obj$env$value.best
  val <- tryCatch(as.numeric(obj$fn(fixed)), error = function(e) NA_real_)
  obj$env$last.par.best <- lpb
  obj$env$last.par      <- lp
  obj$env$value.best    <- vb

  # A re-evaluation that failed is not evidence either way.
  if (!is.finite(val)) return(NA_real_)
  abs(as.numeric(opt$objective) - val)
}
