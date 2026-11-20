# .fit_tmb() ----

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

  # fit_mod() reads its estimates at last.par.best, the lowest objective TMB has
  # evaluated. A Newton step that raised the objective is not that point.
  best_fixed <- function() {
    p <- obj$env$last.par.best
    if (length(obj$env$random)) p[-obj$env$random] else p
  }

  if (requireNamespace("TMBhelper", quietly = TRUE)) {
    fit <- function(start, lo, hi, steps, sd, ctl, loops,
                    elapsed = as.difftime(0, units = "secs")) .normalize_fit_tmb(
      TMBhelper::fit_tmb(
        obj                 = obj,
        fn                  = obj$fn,
        gr                  = obj$gr,
        startpar            = start,
        lower               = lo,
        upper               = hi,
        loopnum             = loops,
        newtonsteps         = steps,
        getsd               = sd,
        control             = ctl,
        bias.correct        = bias.correct,
        bias.correct.control = list(sd = getsd),
        getJointPrecision   = getJointPrecision,
        getReportCovariance = getReportCovariance,
        quiet               = quiet,
        start_time_elapsed  = elapsed
      )
    )
    run <- function(steps) {
      if (!isTRUE(steps > 0)) {
        return(fit(obj$par, lower, upper, steps = 0, sd = getsd, ctl = control,
                   loops = loopnum))
      }
      opt <- fit(obj$par, lower, upper, steps = steps, sd = FALSE, ctl = control,
                 loops = loopnum)
      # TMBhelper reports the last Newton step even when it raised the objective.
      # Re-evaluate at the estimates: iter.max = 0 does not move, and no bounds,
      # because nlminb clamps a start to them.
      ctl <- control
      ctl[["iter.max"]] <- 0
      at_best <- fit(best_fixed(), -Inf, Inf, steps = 0, sd = getsd, ctl = ctl,
                     loops = 1, elapsed = opt[["time_for_run"]])
      keep <- c("iterations", "evaluations", "time_for_MLE")
      at_best[keep] <- opt[keep]
      cols <- c("starting_value", "Lower", "Upper")
      at_best[["diagnostics"]][cols] <- opt[["diagnostics"]][cols]
      at_best
    }
    # TMBhelper's Newton steps solve(H, g) without a guard, so a singular Hessian
    # (a direction the data do not inform) would lose the fit: refit without them.
    return(tryCatch(run(newtonsteps), error = function(e) {
      cl <- conditionCall(e)
      from_solve <- is.call(cl) && identical(cl[[1]], quote(solve.default))
      if (!(isTRUE(newtonsteps > 0) && from_solve)) stop(e)
      warning(sprintf(paste0(
        "newtonsteps = %s met a singular Hessian, so the fit was redone without ",
        "Newton steps (the same fit as newtonsteps = 0). A parameter the data do ",
        "not inform is the usual cause."), newtonsteps), call. = FALSE)
      # Reset TMB's best point so the inner optimizer starts as a fresh fit would.
      obj$env$last.par.best <- obj$env$par
      obj$env$value.best <- Inf
      run(0)
    }))
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
      if (is.null(step)) {
        warning(sprintf(paste0(
          "newtonsteps = %s stopped at a singular Hessian after %d step(s); ",
          "see fit$convergence."), newtonsteps, i - 1L), call. = FALSE)
        break
      }
      par <- par - step
      obj$fn(par)  # as TMBhelper does, so TMB keeps the step if it is the lowest
    }
    best$par <- best_fixed()
    best$objective <- obj$fn(best$par)
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


# Convergence and estimability checks ----

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
#' Indexes with `[[` throughout, so the guard reads exactly the names it means:
#' `$` partially matches on lists, and this discriminates between two shapes
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
  # Minimal in-package fallback: test which fixed effects move the gradient.
  # Param_check uses TMBhelper's own "Bad"/"OK" so one shape leaves this
  # function; the consumer reads that column, not WhichBad, which means
  # eigenvalue positions in TMBhelper and parameter positions here.
  par <- obj$env$last.par.best
  fixed <- setdiff(seq_along(par), obj$env$random)
  g <- as.numeric(obj$gr(par[fixed]))
  bad <- which(abs(g) > 1)
  list(
    BadParams = data.frame(
      Param = names(par[fixed]),
      MLE = unname(par[fixed]),
      Param_check = ifelse(seq_along(fixed) %in% bad, "Bad", "OK"),
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


# Suitability floor diagnostic ----

#' Warn when the empirical-suitability floor bound a runaway
#'
#' `predation.hpp` holds "other food" at a fraction of the modelled-prey term, so
#' a floored M2 is bounded rather than estimated. Names the ages it bound.
#' @param data_list The model's `data_list`.
#' @param quantities The renamed report list.
#' @noRd
.warn_suitability_floor <- function(data_list, quantities) {
  msm <- suppressWarnings(as.integer(as.character(data_list[["msmMode"]])))
  # suit_other is setZero() and written only under predation, so on a
  # single-species fit every cell is 0 and would read as floored.
  if (!length(msm) || !any(msm > 0L, na.rm = TRUE)) return(invisible(NULL))
  so <- quantities[["suit_other"]]
  if (is.null(so)) return(invisible(NULL))

  nspp  <- as.integer(data_list[["nspp"]])
  nsex  <- as.integer(data_list[["nsex"]])
  nages <- as.integer(data_list[["nages"]])
  # Each species' own sexes and ages: the rest of the array is padding at 0.
  if (!length(nspp) || !length(nsex) || !length(nages) ||
      length(dim(so)) != 4L || dim(so)[4L] < 1L) {
    return(invisible(NULL))
  }

  floored <- 1e-6 * 1.001   # within 0.1% of 1e-6/(1 + 1e-6)
  # The floor is in calculate_msvpa_suitability(), so only suitMode 0 reaches
  # it; elsewhere a small suit_other is a fitted vulnerability. NA warns.
  suit <- suppressWarnings(as.integer(as.character(data_list[["suitMode"]])))
  if (length(suit) == 1L) suit <- rep(suit, nspp)
  if (length(suit) != nspp) suit <- rep(NA_integer_, nspp)
  # Ages run minage .. minage + nages - 1 (rule 11), so report the age, not the
  # array index.
  minage <- suppressWarnings(as.integer(as.character(data_list[["minage"]])))
  if (length(minage) == 1L) minage <- rep(minage, nspp)
  if (length(minage) != nspp) minage <- rep(1L, nspp)

  hits <- list()
  for (sp in seq_len(nspp)) {
    if (!is.na(suit[sp]) && suit[sp] > 0L) next
    for (sx in seq_len(nsex[sp])) {
      for (ag in seq_len(nages[sp])) {
        # One value per predator-sex-age, written to every year by
        # predation.hpp, so one year counts each bound cell once.
        v <- so[sp, sx, ag, 1L]
        # > 0 excludes a cell the template never wrote.
        if (!is.finite(v) || v <= 0 || v > floored) next
        hits[[length(hits) + 1L]] <- data.frame(
          sp = sp, sex = sx, age = ag - 1L + minage[sp], value = v)
      }
    }
  }
  if (!length(hits)) return(invisible(NULL))
  h <- do.call(rbind, hits)
  w <- h[which.min(h$value), ]
  spname <- data_list[["spnames"]][w$sp]
  if (is.null(spname) || is.na(spname)) spname <- paste("species", w$sp)

  warning(
    "The 'other food' suitability floor bound ", nrow(h),
    " predator-sex-age cell(s); the lowest is ", spname, " age ", w$age,
    # Named only on a two-sex predator: on one sex it is noise, and the count
    # above already promises the sex axis.
    if (length(nsex) >= w$sp && isTRUE(nsex[w$sp] > 1L))
      paste0(", sex ", w$sex) else "",
    ". A prey age with (near) zero biomass and a positive diet ",
    "proportion has taken almost the whole suitability share, so available ",
    "food, M2 and consumption for that predator age are bounded rather than ",
    "estimated. Check `diet_data` against that prey's numbers at age, or widen ",
    "`suit_styr`/`suit_endyr`.",
    call. = FALSE)
  invisible(h)
}
