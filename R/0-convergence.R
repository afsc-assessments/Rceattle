# -----------------------------------------------------------------------------
# Convergence diagnostics
# -----------------------------------------------------------------------------
# Each check emits one record -- list(id, tier, severity, message, data),
# severity in c("OK","NOTE","WARN","FAIL") -- so the battery is a single object
# (fit$convergence) rather than scattered messages. fit_mod() attaches the
# result and surfaces non-OK checks via message() (never errors).
# -----------------------------------------------------------------------------

# Severity ordering (low -> high). Used to roll individual records up into an
# overall status (the worst severity present).
.CONV_SEVERITY <- c("OK", "NOTE", "WARN", "FAIL")

# Build one check record.
.conv_record <- function(id, tier, severity, message, data = NULL) {
  severity <- match.arg(severity, .CONV_SEVERITY)
  list(id = id, tier = tier, severity = severity,
       message = message, data = data)
}

# Worst severity across a list of records ("OK" when empty).
.conv_overall <- function(checks) {
  if (length(checks) == 0) return("OK")
  sev <- vapply(checks, function(x) x$severity, character(1))
  .CONV_SEVERITY[max(match(sev, .CONV_SEVERITY))]
}

# Snapshot of the hindcast optimizer result, captured in fit_mod() before the
# projection re-optimization overwrites `opt`. `obj` lets us recompute the
# gradient and MLE vector directly -- essential in the hard-failure case, where
# opt$SD/$max_gradient/$diagnostics are often absent. `getsd` flags whether an
# sdreport was requested.
#
# `bounds` (the full build_bounds() list, $lower/$upper), `mapFactor`
# (map$mapFactor) and `random_vars` are used to assemble the per-parameter
# (mle, lower, upper) table.
.capture_opt_convergence <- function(opt, obj = NULL, bounds = NULL,
                                     mapFactor = NULL, random_vars = NULL,
                                     getsd = NA, data_list = NULL) {
  if (is.null(opt) && is.null(obj)) return(NULL)

  # Fixed-effect MLE vector (named), recomputed from the object.
  par_fixed <- NULL
  if (!is.null(obj)) {
    par_fixed <- tryCatch({
      par   <- obj$env$last.par.best
      ran   <- obj$env$random
      fixed <- if (length(ran)) setdiff(seq_along(par), ran) else seq_along(par)
      stats::setNames(as.numeric(par[fixed]), names(par)[fixed])
    }, error = function(e) NULL)
  }

  diag  <- opt$diagnostics
  worst <- NULL
  gg <- NULL
  if (!is.null(diag) && !is.null(diag$final_gradient)) {
    # Named by block so the index can say where the largest gradient sits.
    gg <- stats::setNames(as.numeric(diag$final_gradient),
                          as.character(diag$Param))
    i <- which.max(abs(diag$final_gradient))
    if (length(i) == 1L) {
      worst <- list(param = as.character(diag$Param[i]),
                    gradient = as.numeric(diag$final_gradient[i]))
    }
  }
  mg <- opt$max_gradient
  if (is.null(mg) || !is.finite(mg)) {
    mg <- if (!is.null(worst)) abs(worst$gradient) else NA_real_
  }

  # Fallback: recompute the marginal gradient (and the worst parameter) from the
  # hindcast object when opt did not provide it.
  if ((is.na(mg) || is.null(worst)) && !is.null(par_fixed)) {
    gg <- tryCatch(stats::setNames(as.numeric(obj$gr(par_fixed)),
                                   names(par_fixed)),
                   error = function(e) NULL)
    if (!is.null(gg) && length(gg) > 0 && any(is.finite(gg))) {
      i <- which.max(abs(gg))
      worst <- list(param = names(gg)[i], gradient = unname(gg[i]))
      if (is.na(mg)) mg <- max(abs(gg))
    }
  }

  # Per-parameter (mle, lower, upper) for the bounds check. Build all three from
  # a single shared shape so they cannot drift out of alignment.
  bnd_par <- bnd_lo <- bnd_hi <- NULL
  if (!is.null(bounds) && !is.null(bounds$lower) && !is.null(bounds$upper) &&
      !is.null(mapFactor) && !is.null(obj)) {
    # parList(x = par[lfixed()], par = last.par): pass last.par.best as `par`
    # (the full vector) so `x` defaults to its fixed-effects subset. Passing it
    # positionally as `x` triggers `par[lfixed()] <- x` with a length mismatch.
    full <- tryCatch(obj$env$parList(par = obj$env$last.par.best),
                     error = function(e) NULL)
    if (!is.null(full)) {
      pv <- lv <- uv <- nv <- list()
      for (nm in names(mapFactor)) {
        if (nm %in% random_vars) next                 # random effects are unbounded
        if (is.null(full[[nm]]) || is.null(bounds$lower[[nm]])) next
        mf <- as.numeric(unlist(mapFactor[[nm]]))
        pb <- as.numeric(unlist(full[[nm]]))
        lb <- as.numeric(unlist(bounds$lower[[nm]]))
        ub <- as.numeric(unlist(bounds$upper[[nm]]))
        # one row per estimated coefficient (drop NA-mapped/fixed and duplicates)
        keep <- which(!is.na(mf) & !duplicated(mf))
        if (length(keep) == 0L) next
        pv[[nm]] <- pb[keep]; lv[[nm]] <- lb[keep]; uv[[nm]] <- ub[keep]
        nv[[nm]] <- rep(nm, length(keep))
      }
      if (length(pv)) {
        bnd_par <- stats::setNames(unlist(pv, use.names = FALSE),
                                   unlist(nv, use.names = FALSE))
        bnd_lo  <- unlist(lv, use.names = FALSE)
        bnd_hi  <- unlist(uv, use.names = FALSE)
      }
    }
  }

  # Where each hindcast parameter sits in the model. Built here because an HCR
  # projection replaces obj with one holding log_Ftarget / log_Flimit alone.
  index <- if (!is.null(obj) && !is.null(data_list)) {
    tryCatch(.rce_par_index(obj, data_list), error = function(e) NULL)
  } else NULL

  # Standard errors of the hindcast fixed effects, in obj$par order. Kept here
  # because under any estimating HCR `fit$sdrep` is the projection's.
  se_fixed <- NULL
  if (!is.null(opt$SD) && !is.null(opt$SD$cov.fixed)) {
    se_fixed <- tryCatch(
      stats::setNames(.conv_se_from_cov(opt$SD$cov.fixed), names(opt$SD$par.fixed)),
      error = function(e) NULL)
  }

  list(
    par                = if (!is.null(bnd_par)) bnd_par else par_fixed,
    index              = index,
    gradient           = gg,
    lower              = bnd_lo,
    upper              = bnd_hi,
    se_fixed           = se_fixed,
    sd_requested       = isTRUE(getsd),
    sd_present         = !is.null(opt$SD),
    max_gradient       = as.numeric(mg),
    worst              = worst,
    nlminb_convergence = opt$convergence,
    message            = opt$message,
    convergence_check  = opt$Convergence_check,
    pdHess             = if (!is.null(opt$SD)) isTRUE(opt$SD$pdHess) else NA
  )
}


#' Did a diagnostic re-fit converge well enough to keep?
#'
#' @description
#' The shared keep/drop gate for the re-fitting diagnostics,
#' [retrospective()], [jitter()], [self_test()] and [profile.Rceattle()], each
#' of which silently drops the runs that did not converge.
#'
#' These call sites used to test `opt$Convergence_check` against the string
#' `TMBhelper::fit_tmb()` uses for a non-invertible Hessian. `fit_tmb()` assigns
#' that particular string in exactly one place, when `sdreport` returns
#' `pdHess = FALSE`, and the test could not work in either direction:
#'
#' * with `getsd = TRUE`, `fit_tmb()` returns early when the Hessian fails
#'   `chol()`, so it never reaches that assignment, and the shape it returns
#'   instead holds no `Convergence_check` at all, the run was dropped by the
#'   enclosing `is.null()` guard, by accident rather than by the test;
#' * with `getsd = FALSE` the assignment is unreachable, so *nothing* was ever
#'   dropped, a run that ended with a maximum gradient of 1e13 counted as
#'   converged. (`Convergence_check` is still set, but to one of the two gradient
#'   verdicts, neither of which the test matched.)
#'
#' So judge the gradient, which is available either way.
#' `.capture_opt_convergence()` recomputes it from the objective function when
#' the optimizer did not report one, and `fit_mod()` captures that snapshot from
#' the *hindcast* optimization, before any projection re-optimization overwrites
#' `opt`.
#'
#' This judges the model it is handed, so a caller that re-fits more than once
#' has to call it more than once: [retrospective()] runs a second, F-only refit
#' whose gradient says nothing about whether the peeled hindcast converged, and
#' gates on both.
#'
#' The threshold is [convergence_diagnostics()]'s own `max_gradient` FAIL tier
#' and compares the same way (`>`, not `>=`), so the two cannot disagree at the
#' boundary. Deliberately the FAIL tier and not the tighter WARN tier (1e-3):
#' dropping a run is destroying it, so the bar for that is "this fit is broken",
#' not "this fit is worth a look".
#'
#' Two consequences worth being explicit about, since dropping is silent to
#' anything that only reads `length()`:
#'
#' * this is an OPTIMIZER gate, not the whole battery. A kept run can still
#'   hold a WARN (gradient between 1e-3 and 1) or even a FAIL from one of the
#'   other checks, a non-positive-definite Hessian, a non-identifiable
#'   parameter, a stock-recruit curve under the replacement line. Read
#'   `$convergence` on what comes back; do not treat "returned" as "clean";
#' * the one case that drops without a matching battery record is a non-finite
#'   gradient. `.check_optimizer()` emits nothing at all then (it guards on
#'   `is.finite`), so the battery can say OK for a run this drops. A `NaN`
#'   gradient is a diverged fit, which is exactly what the gate is for.
#'
#' This is the gate for the diagnostics that DROP runs. [run_mse()] also re-fits
#' through `.refit_like()` and is deliberately not gated: an MSE cannot skip an
#' assessment year, so a bad estimation-model fit has to propagate rather than
#' vanish.
#'
#' @param newmod A fitted `Rceattle`, or `NULL`.
#' @param max_grad Drop the run above this maximum absolute marginal gradient.
#'   Defaults to 1, `.check_optimizer()`'s FAIL tier.
#' @return `TRUE` to keep the run, `FALSE` to drop it.
#' @noRd
.refit_converged <- function(newmod, max_grad = 1) {
  if (is.null(newmod)) return(FALSE)

  ch <- newmod[[".conv_hindcast"]]
  if (is.null(ch)) {
    # No hindcast snapshot, so no gradient to judge: fit_mod() captures one only
    # for estimateMode 0 and 1. Defer to the optimizer's own verdict, which is
    # what these call sites did before. All four rewrite estimateMode to 0/1/2,
    # so in practice this branch is only reached at >= 3, where `opt` is never
    # stored either and the run drops. It is NOT a safe general rule: at mode 2
    # `opt` is attached from the PROJECTION optimization, so a fifth caller
    # would be keeping a run on a verdict about its reference points.
    cc <- newmod[["opt"]][["Convergence_check"]]
    return(!is.null(cc) && cc != "The model is definitely not converged")
  }

  mg <- ch[["max_gradient"]]
  if (length(mg) != 1L || !is.finite(mg) || mg > max_grad) return(FALSE)

  # A Hessian that is not positive definite -- the condition the old string test
  # was reaching for. Keep dropping on it, but only when an sdreport was actually
  # asked for, since otherwise there is nothing to judge.
  #
  # Both ways it can present: TMBhelper::fit_tmb() bails at its own chol() and
  # never returns an SD at all, while the in-package .fit_tmb() fallback calls
  # TMB::sdreport() directly -- which does NOT error on a singular Hessian, it
  # returns an object with pdHess = FALSE. Testing only for a missing SD would
  # keep on the fallback path exactly what it drops on the TMBhelper path.
  if (isTRUE(ch[["sd_requested"]]) &&
      (!isTRUE(ch[["sd_present"]]) || isFALSE(ch[["pdHess"]]))) {
    return(FALSE)
  }

  TRUE
}


#' Tell the user that a diagnostic dropped runs
#'
#' @description
#' The re-fitting diagnostics drop non-converged runs and return a short list.
#' While the keep/drop gate could not actually drop anything (see
#' `.refit_converged()`) that silence cost nothing; now that it can, a caller
#' who does not think to compare `length()` against what they asked for would
#' read a thinned list as a complete one, and for `jitter()` and `self_test()`
#' a thinned list is a biased sample, since the runs that failed are exactly the
#' ones that would have shown the spread.
#'
#' @param n_dropped,n_total Counts for the message.
#' @param what Singular noun for the unit, e.g. `"peel"`.
#' @param warn Raise a warning rather than a message. `retrospective()` sets it
#'   because Mohn's rho is averaged over the peels that survive, so a drop
#'   changes a reported number; `suppressMessages()` is common enough in the
#'   assessment scripts that a message would not reach the reader.
#' @return `invisible(NULL)`, called for the message.
#' @noRd
.report_dropped <- function(n_dropped, n_total, what, warn = FALSE) {
  if (!isTRUE(n_dropped > 0)) return(invisible(NULL))
  txt <- sprintf(
    "%d of %d %s%s dropped as non-converged (max|gradient| > 1, or sdreport failed); %d returned.",
    n_dropped, n_total, what, if (n_dropped == 1L) "" else "s",
    n_total - n_dropped)
  if (isTRUE(warn)) warning(txt, call. = FALSE) else message(txt)
  invisible(NULL)
}


#' Run one diagnostic re-fit under an optional wall-clock limit
#'
#' @description
#' `.fit_tmb()` optimizes with `eval.max = iter.max = 1e9`, so a re-fit that
#' wanders somewhere pathological has no bound and one replicate can stall a
#' whole `jitter()` or `self_test()` run, the failure this is for is a hang,
#' which no convergence check can reach because the fit never returns.
#'
#' The limit is approximate by construction: [setTimeLimit()] is checked when
#' control returns to R, so it fires between the optimizer's function
#' evaluations rather than inside one. That is enough here (`nlminb` re-enters R
#' every evaluation) but a single very long evaluation can overrun it.
#'
#' Errors, including the timeout, are returned rather than thrown, so one bad
#' replicate cannot abort the run and, under a cluster, take every other
#' replicate with it.
#'
#' @param expr The re-fit, evaluated lazily.
#' @param timeout Elapsed-second limit, or `Inf` for none.
#' @return The fitted model, or the `condition` if it errored or timed out.
#' @noRd
.refit_with_timeout <- function(expr, timeout = Inf) {
  if (is.finite(timeout) && timeout > 0) {
    # transient = FALSE plus an explicit reset: `transient = TRUE` would restore
    # the limit only at the end of the whole top-level call (the lapply over all
    # replicates), not at the end of this one, so every replicate after the first
    # would inherit a clock that started before it did.
    setTimeLimit(elapsed = timeout, transient = FALSE)
    on.exit(setTimeLimit(cpu = Inf, elapsed = Inf), add = TRUE)
  }
  tryCatch(expr, error = function(e) e)
}


#' Is this condition a `.refit_with_timeout()` timeout rather than a model error?
#' @noRd
.is_timeout <- function(e) {
  inherits(e, "condition") &&
    grepl("reached elapsed time limit|reached CPU time limit", conditionMessage(e))
}


#' Report replicates that errored or timed out
#'
#' Separate from `.report_dropped()` on purpose: a hang or a hard error is a
#' different problem from a bad gradient, and the fix is different too (raise
#' `timeout`, versus look at the model).
#'
#' @param errs Character vector of messages, `NA` where the replicate ran.
#' @param what Singular noun for the unit.
#' @return `invisible(NULL)`, called for the message.
#' @noRd
.report_errors <- function(errs, what) {
  bad <- errs[!is.na(errs)]
  if (length(bad) == 0L) return(invisible(NULL))
  n_to <- sum(grepl("reached elapsed time limit|reached CPU time limit", bad))
  if (n_to > 0L) {
    message(sprintf("%d %s%s exceeded `timeout` and were stopped; raise it to let them finish.",
                    n_to, what, if (n_to == 1L) "" else "s"))
  }
  if (length(bad) > n_to) {
    message(sprintf("%d %s%s errored; first: %s",
                    length(bad) - n_to, what, if (length(bad) - n_to == 1L) "" else "s",
                    bad[!grepl("reached elapsed time limit|reached CPU time limit", bad)][1]))
  }
  invisible(NULL)
}


# --- checks ------------------------------------------------------------------

# Where one parameter of a named vector sits: "log_M1 (M1): age 1", or the
# display name alone when no index matches the vector.
.conv_par_where <- function(i, nms, index) {
  txt <- .rce_par_display(nms[i])
  idx <- .conv_index_for(index, nms)
  if (!is.null(idx)) {
    lab <- idx$label[match(i, idx$par_index)]
    if (length(lab) == 1L && !is.na(lab) && nzchar(lab)) {
      txt <- paste0(txt, ": ", lab)
    }
  }
  txt
}

# How far the estimates still are from the optimum, in standard errors. A
# quadratic approximation puts the optimum a Newton step cov %*% gradient away;
# dividing each element by its own SE makes the size comparable across log,
# logit and natural-scale parameters, which the raw gradient is not. NULL
# unless cov.fixed describes the same parameter vector as the hindcast gradient
# (under an estimating HCR the sdreport is the projection's).
# Standard errors from a covariance diagonal, without warning on a negative
# variance. An indefinite Hessian gives one, and this battery reports through
# message() and must not emit "NaNs produced" from inside a diagnostic; the
# callers already drop non-finite entries.
.conv_se_from_cov <- function(cov) {
  v <- diag(cov)
  v[!is.finite(v) | v < 0] <- NA_real_
  sqrt(v)
}

.conv_newton_step_se <- function(object, gg) {
  cov <- tryCatch(object$sdrep$cov.fixed, error = function(e) NULL)
  if (is.null(cov) || !is.matrix(cov) || is.null(gg) || is.null(names(gg)) ||
      nrow(cov) != length(gg)) return(NULL)
  # Only reported on a positive-definite Hessian. The step is the distance to a
  # minimum only if `cov` inverts one; on a saddle it points away from it, and
  # "the optimum is a fraction of a standard error away" would read as
  # reassurance beside the pdHess check that just failed.
  if (!isTRUE(object$sdrep$pdHess)) return(NULL)
  nm <- rownames(cov)
  if (is.null(nm)) nm <- names(object$sdrep$par.fixed)
  if (!identical(unname(as.character(nm)), names(gg))) return(NULL)
  # A negative variance is dropped rather than square-rooted: this battery
  # reports through message() and never raises, and sqrt() would emit "NaNs
  # produced" from inside a diagnostic.
  se   <- unname(.conv_se_from_cov(cov))
  step <- as.numeric(cov %*% gg) / se
  if (!any(is.finite(step))) return(NULL)
  i <- which.max(abs(replace(step, !is.finite(step), NA)))
  list(max = abs(step[i]), i = i, step_se = stats::setNames(step, names(gg)))
}

# Optimizer convergence: max |gradient| (+ the parameter carrying it) and
# Hessian positive-definiteness. Reads the hindcast snapshot so the result is
# not clobbered by the projection re-optimization.
.check_optimizer <- function(object, index = NULL) {
  ch <- object$.conv_hindcast
  out <- list()
  if (is.null(ch)) return(out)

  mg <- ch$max_gradient
  if (!is.null(mg) && is.finite(mg)) {
    sev <- if (mg > 1) "FAIL" else if (mg > 1e-3) "WARN" else "OK"
    gg  <- ch$gradient
    # Coordinates cost a parList() pass, so a clean gradient does not force them.
    if (sev == "OK") index <- NULL
    worst_txt <- if (!is.null(gg) && !is.null(names(gg)) && any(is.finite(gg))) {
      sprintf(" on %s", .conv_par_where(which.max(abs(gg)), names(gg), index))
    } else if (!is.null(ch$worst)) {
      sprintf(" on %s", .rce_par_display(ch$worst$param))
    } else ""
    step <- .conv_newton_step_se(object, gg)
    step_txt <- if (!is.null(step)) {
      sprintf(" Reaching the optimum would move the estimates by at most %.2g standard errors (%s).",
              step$max, .conv_par_where(step$i, names(gg), index))
    } else ""
    out$max_gradient <- .conv_record(
      "max_gradient", "fit", sev,
      sprintf("Maximum absolute marginal gradient = %.3g%s.%s",
              mg, worst_txt, step_txt),
      c(ch, list(newton_step_se = step$max, step_se = step$step_se)))
  }

  if (!is.null(ch$pdHess) && !is.na(ch$pdHess)) {
    out$pdHess <- .conv_record(
      "pdHess", "fit",
      if (isTRUE(ch$pdHess)) "OK" else "FAIL",
      if (isTRUE(ch$pdHess)) "Hessian is positive definite."
      else "Hessian is not positive definite; standard errors are unavailable.")
  }
  out
}

# How nearly linearly dependent the estimates are. Reads the fixed-effect
# covariance (sdrep$cov.fixed), standardized to a correlation; its eigenvalues
# are variances along directions in parameter space, so kappa = max/min is the
# severity and the eigenvector of the largest is the least-determined direction.
# Complements (does not duplicate) TMBhelper::check_estimability: that gives a
# per-parameter verdict, this a continuous severity and the linear *combination*.
.check_hessian_eigen <- function(object, index = .conv_par_index(object)) {
  out <- list()
  cov <- tryCatch(object$sdrep$cov.fixed, error = function(e) NULL)
  if (is.null(cov) || !is.matrix(cov) || nrow(cov) < 2L) return(out)
  # Severity is read on the CORRELATION matrix. The covariance condition number
  # is not scale-invariant -- log_F sits near -2 and sel_inf near 10, so
  # rescaling a parameter moves it without changing the model. Standardising
  # leaves a measure of how nearly linearly dependent the estimates are.
  # The covariance number is still reported, as the numerical cost of inverting
  # the Hessian; see inst/dev/TRAPS.md for the measured values.
  se <- .conv_se_from_cov(cov)
  ok <- is.finite(se) & se > 0
  if (sum(ok) < 2L) return(out)
  keep_i <- which(ok)
  corr <- cov[ok, ok, drop = FALSE] / outer(se[ok], se[ok])

  cov_ev <- tryCatch(eigen(cov[ok, ok, drop = FALSE], symmetric = TRUE,
                           only.values = TRUE), error = function(e) NULL)
  cov_kappa <- NA_real_
  if (!is.null(cov_ev)) {
    cp <- cov_ev$values[is.finite(cov_ev$values) & cov_ev$values > 0]
    if (length(cp) >= 2L) cov_kappa <- max(cp) / min(cp)
  }

  ev <- tryCatch(eigen(corr, symmetric = TRUE), error = function(e) NULL)
  if (is.null(ev)) return(out)
  vals <- ev$values
  pos  <- vals[is.finite(vals) & vals > 0]
  if (length(pos) < 2L) return(out)

  kappa <- max(pos) / min(pos)
  # Correlation eigenvalues are variances of the STANDARDIZED estimates, so the
  # square root is a ratio of standard errors each measured against its own
  # parameter's -- the readable form of the same number.
  se_ratio <- sqrt(kappa)
  v  <- ev$vectors[, which.max(vals)]         # least-determined direction (unit norm)
  nm <- rownames(cov)
  if (is.null(nm) || length(nm) != length(cov[, 1])) nm <- names(object$sdrep$par.fixed)
  if (is.null(nm) || length(nm) != length(cov[, 1])) nm <- paste0("p", seq_len(nrow(cov)))
  nm_full <- nm                 # one name per row of cov, before the finite-se subset
  nm <- nm[keep_i]

  # The eigenvector is unit-norm, so v^2 partitions the direction across
  # coefficients. Aggregate by parameter block (a vector parameter shares one
  # name across its coefficients): each block's share is the sum of its squared
  # loadings, and the shares sum to 1. This names the block(s) the flat direction
  # lives in even when it is spread diffusely over many coefficients -- the case
  # a single-coefficient threshold missed, leaving the message blank.
  share <- tapply(v^2, factor(nm, levels = unique(nm)), sum)
  share <- sort(share, decreasing = TRUE)
  cum   <- cumsum(as.numeric(share))
  ntop  <- which(cum >= 0.90)[1]                          # blocks explaining >=90%
  if (is.na(ntop)) ntop <- length(share)
  ntop  <- max(1L, min(ntop, 5L))                         # always name >=1, cap at 5
  combo <- paste(sprintf("%s: %.0f%%", .rce_par_display(names(share)[seq_len(ntop)]),
                         100 * as.numeric(share)[seq_len(ntop)]),
                 collapse = " + ")
  top   <- data.frame(param = names(share), share = round(as.numeric(share), 3))

  sev <- if (kappa > 1e10) "FAIL" else if (kappa > 1e6) "WARN" else "OK"
  msg <- sprintf(
    "Condition number = %.2g: standard errors span %.2gx across parameter combinations%s",
    kappa, se_ratio,
    if (sev != "OK") sprintf("; loads on: %s.", combo) else ".")

  # Where the loading sits INSIDE each named block. The direction is spread over
  # coefficients rather than confined to a set, so the parameters carrying the
  # block's top 90% are summarised by coordinate and the count is printed: a
  # ridge over ten years and one over a single year read very differently.
  if (sev != "OK") {
    # cov.fixed's rows are the sdreport's parameter vector, which is the
    # projection fit's under an estimating HCR and the hindcast's otherwise.
    idx <- .conv_index_for(index, nm_full)
    if (!is.null(idx)) {
      lines <- character(0)
      # One width for every block printed below, so the coordinates align.
      .w <- max(16L, nchar(.rce_par_display(names(share)[seq_len(ntop)])))
      for (b in names(share)[seq_len(ntop)]) {
        inb <- which(nm == b)
        o   <- inb[order(v[inb]^2, decreasing = TRUE)]
        keep <- o[seq_len(max(1L, which(cumsum(v[o]^2) / sum(v[o]^2) >= 0.90)[1]))]
        # `v` and `nm` index the parameters kept above, so map back before
        # asking parameter_index() where they are.
        s <- .rce_par_summary(keep_i[keep], idx, max_lines = 2L, width = .w)
        if (length(s) > 0) {
          lines <- c(lines, sub("\\((\\d+)\\)$",
                                sprintf("(\\1 of %d parameters; %.0f%% of the direction)",
                                        length(inb), 100 * share[[b]]), s))
        }
      }
      if (length(lines) > 0) msg <- paste0(msg, "\n", paste(lines, collapse = "\n"))
    }
  }
  out$hessian_conditioning <- .conv_record(
    "hessian_conditioning", "fit", sev, msg,
    list(condition_number = kappa, se_ratio = se_ratio,
         covariance_condition_number = cov_kappa, loadings = top))
  out
}

# The template floors numbers-at-age, the Ricker intercept and (under an
# identity-link recruitment linkage) recruitment at 0.001 with a penalty, so a
# non-zero row means the fit is of a floored model, not the one specified.
.check_zero_n_penalty <- function(object) {
  jc <- object$quantities$jnll_comp
  if (is.null(jc) || !"Zero n-at-age penalty" %in% rownames(jc)) return(list())
  pen <- jc["Zero n-at-age penalty", ]
  # posfun() charges 0.01 * (x - 0.001)^2 per floored value, so sqrt(pen / 0.01)
  # is the root-sum-square excursion (thousands of fish for N and R). Under 1e-3
  # is numerical.
  excursion <- sqrt(pmax(pen, 0) / 0.01)
  bad <- !is.finite(pen)
  hit <- which(bad | excursion > 1e-3)
  if (!length(hit)) return(list())
  sp <- (object$data_list$spnames %||% seq_along(pen))[hit]
  severity <- if (any(bad) || max(excursion[is.finite(excursion)], 0) > 1) "FAIL" else "WARN"
  list(zero_n_penalty = .conv_record(
    "zero_n_penalty", "fit", severity,
    sprintf("Numbers-at-age, the Ricker intercept or recruitment sat on the 0.001 floor in species %s (root-sum-square excursion %s); the fit is of a floored model, not the one specified.",
            paste(sp, collapse = ", "),
            paste(ifelse(bad[hit], "not finite", signif(excursion[hit], 3)), collapse = ", ")),
    list(penalty = pen, excursion = excursion)))
}

# getsd = FALSE leaves sdrep NULL, so the Hessian eigenvalue, sdreport, pdHess
# and estimability checks all return nothing. Without this record the battery
# reports "OK" and report_tables() prints it into a SAFE table, which reads as
# "every check passed" rather than "the strongest checks never ran".
.check_hessian_not_run <- function(object) {
  ch <- object$.conv_hindcast
  # fit_mod() always records sd_requested as a logical, so a positive FALSE is
  # what marks getsd = FALSE. An absent field is an older or synthetic fit
  # object and says nothing either way, so it earns no record.
  if (is.null(ch) || !identical(ch$sd_requested, FALSE)) return(list())
  list(hessian_not_run = .conv_record(
    "hessian_not_run", "fit", "NOTE",
    paste0("Hessian checks not run: the fit was made with getsd = FALSE, so ",
           "the positive-definite Hessian, condition-number, sdreport and ",
           "estimability checks were all skipped. Refit with getsd = TRUE ",
           "before reading this status as convergence.")))
}

# sdreport failed: requested but did not return (Hessian not invertible). A
# strong non-convergence signal even when no gradient is available.
.check_sdreport_failed <- function(object) {
  ch <- object$.conv_hindcast
  if (is.null(ch) || !isTRUE(ch$sd_requested) || isTRUE(ch$sd_present)) {
    return(list())
  }
  list(sdreport_failed = .conv_record(
    "sdreport_failed", "fit", "FAIL",
    "sdreport failed (Hessian not invertible); standard errors unavailable."))
}

# Parameters at a configured bound. Optimization is unbounded in fit_mod(), so a
# parameter at/beyond its build_bounds() range means the MLE hit the edge of the
# plausible range -- often unidentified or mis-scaled.
.check_bounds <- function(object, index = .conv_par_index(object)) {
  ch <- object$.conv_hindcast
  if (is.null(ch) || is.null(ch$par) || is.null(ch$lower) || is.null(ch$upper)) {
    return(list())
  }
  par <- ch$par; lo <- ch$lower; hi <- ch$upper
  rng <- hi - lo
  tol <- pmax(1e-6, 1e-3 * ifelse(is.finite(rng) & rng > 0, rng, 1))
  at_lo <- is.finite(lo) & par <= lo + tol
  at_hi <- is.finite(hi) & par >= hi - tol
  hit <- which((at_lo | at_hi) & par > -900)   # skip -999 sentinels (e.g. log_F)
  if (length(hit) == 0) return(list())
  tab <- data.frame(param = names(par)[hit], mle = signif(par[hit], 4),
                    lower = signif(lo[hit], 4), upper = signif(hi[hit], 4),
                    bound = ifelse(at_lo[hit], "lower", "upper"))
  # `hit` counts positions in the hindcast bounds vector, so the index has to be
  # the one built from it.
  idx <- .conv_index_for(index, names(par))
  tab <- .conv_attach_label(tab, hit, idx)
  list(parameters_on_bounds = .conv_record(
    "parameters_on_bounds", "fit", "WARN",
    .conv_with_coords(
      sprintf("%d parameter(s) at a configured bound: %s.",
              length(hit), paste(unique(names(par)[hit]), collapse = ", ")),
      hit, idx),
    tab))
}

#' A process variance estimated to zero
#'
#' A deviation standard deviation at the floor means the process it governs is
#' not varying: a model configured as time-varying has fitted something
#' time-invariant. `Atka2022` under `random_sel = TRUE` with a non-parametric
#' random walk reaches `sel_dev_sd = 2.7e-08`. The battery flags that particular
#' fit through `max_gradient`, which reports that the optimizer stopped, not what
#' went wrong; and a collapse at a CLEAN gradient, a well-posed maximum at the
#' boundary, has nothing else to catch it.
#'
#' Scope is deliberately narrow. Only the standard deviations of a modelled
#' DEVIATION are read, all of which are log-scale and O(0.1)-O(1) in any
#' assessment, so one absolute threshold means the same thing for all of them.
#' Observation sds (`catch_log_sd`, `index_log_sd`, `log_obs_sd_linkage`) are
#' excluded: a concentrated observation sd at zero is a different diagnosis, it
#' already trips `parameters_on_bounds` because `build_bounds()` bounds those,
#' and it governs no time-varying process. `growth_log_sd` is excluded because it
#' is a length-at-age sd in CENTIMETRES, where 1e-3 means nothing. Correlations
#' and penalty weights (`sel_curve_pen`, `M1_rho`,
#' `trans_rho_linkage`) are not sds at all, and the catchability PRIOR sd
#' (`index_q_log_sd`) is a deliberate choice.
#'
#' `WARN`, never `FAIL`: a variance at the boundary is a well-posed optimum and a
#' verdict about the MODEL, not about the optimizer, and `report_tables()` prints
#' this status as a fit's `converged` field.
#'
#' Only estimated parameters appear in the snapshot, so a deviation sd held at
#' `Time_varying_sel_sd` cannot fire this.
#' @noRd
.CONV_PROCESS_SD <- c("R_log_sd", "sel_dev_log_sd", "index_q_dev_log_sd",
                      "M1_dev_log_sd", "log_sigma_linkage")

#' @param object A fitted `Rceattle` object.
#' @return A list holding zero or one check record.
#' @keywords internal
#' @noRd
.check_variance_collapse <- function(object) {
  ch <- object$.conv_hindcast
  if (is.null(ch) || is.null(ch$par) || is.null(names(ch$par))) return(list())

  par   <- ch$par
  is_sd <- names(par) %in% .CONV_PROCESS_SD
  if (!any(is_sd)) return(list())

  sd_val <- exp(par[is_sd])
  hit    <- which(sd_val < 1e-3)
  if (length(hit) == 0) return(list())

  # A block is named once per element, so say WHICH element: "R_log_sd[2]" is
  # actionable where "R_log_sd" on a three-species model is not.
  nm  <- names(sd_val)
  idx <- stats::ave(seq_along(nm), nm, FUN = seq_along)
  n_in_block <- as.integer(table(nm)[nm])
  label <- ifelse(n_in_block > 1L, sprintf("%s[%d]", nm, idx), nm)

  tab <- data.frame(param = nm[hit], index = idx[hit],
                    sd = signif(unname(sd_val[hit]), 3), row.names = NULL)

  list(variance_collapse = .conv_record(
    "variance_collapse", "fit", "WARN",
    sprintf(paste0("%d estimated deviation standard deviation(s) at zero: %s. ",
                   "The deviations they govern are not varying, so this fit is ",
                   "time-invariant in that process despite being configured ",
                   "otherwise."),
            length(hit),
            paste(sprintf("%s = %.2g", label[hit], sd_val[hit]), collapse = ", ")),
    tab))
}


# Phasing log: phases that ended with a high gradient localize which parameter
# block breaks. Reads the per-phase log captured from TMBphase().
.check_phasing <- function(object) {
  pl <- object$.conv_phase
  if (is.null(pl) || length(pl) == 0) return(list())
  df <- tryCatch(do.call(rbind, lapply(pl, as.data.frame)),
                 error = function(e) NULL)
  if (is.null(df) || is.null(df$max_grad)) return(list())
  bad <- df[is.finite(df$max_grad) & df$max_grad > 1, , drop = FALSE]
  if (nrow(bad) == 0) return(list())
  list(phasing = .conv_record(
    "phasing", "fit", "WARN",
    sprintf("Phase(s) ended with max|grad| > 1: %s.",
            paste(sprintf("phase %s (%.2g)", bad$phase, bad$max_grad),
                  collapse = ", ")),
    df))
}

# Surface TMBhelper::check_estimability (run by fit_mod() and stored as
# object$identified when sdreport fails). Binary per-parameter verdict; pairs
# with .check_hessian_eigen()'s continuous view.
.check_estimability_record <- function(object, index = .conv_par_index(object)) {
  out <- list()
  id <- object$identified
  if (is.null(id) || is.character(id)) return(out)
  bad <- tryCatch(id$WhichBad, error = function(e) NULL)
  if (is.null(bad)) {
    pc  <- tryCatch(id$BadParams$Param_check, error = function(e) NULL)
    bad <- if (!is.null(pc)) which(pc == 2L) else NULL
  }
  if (!is.null(bad) && length(bad) > 0) {
    nms <- tryCatch(unique(as.character(id$BadParams$Param[bad])),
                    error = function(e) NULL)
    # BadParams has one row per hindcast fixed parameter, in obj$par order, so
    # its Param column is the vector `bad` indexes into.
    idx <- .conv_index_for(index, tryCatch(as.character(id$BadParams$Param),
                                           error = function(e) NULL))
    out$estimability <- .conv_record(
      "estimability", "fit", "FAIL",
      .conv_with_coords(
        sprintf("check_estimability flagged %d non-identifiable fixed parameter(s)%s.",
                length(bad),
                if (!is.null(nms)) paste0(" (", paste(nms, collapse = ", "), ")")
                else ""),
        bad, idx),
      id)
  }
  out
}

# Stock-recruit curve sanity. Beverton-Holt steepness is derived as
# h = alpha * SPR0 / (4 + alpha * SPR0), and alpha carries no lower bound, so the
# optimizer can reach alpha < 1/SPR0 -- below the replacement line. Steepness
# then falls under 0.2, the stock cannot replace itself, and the implied unfished
# recruitment (alpha - 1/SPR0)/beta is negative, which carries into the initial
# age structure.
.check_stock_recruit <- function(object) {
  dl <- object[["data_list"]]
  q  <- object[["quantities"]]
  if (is.null(dl) || is.null(q)) return(list())
  # Read the curve whose parameters are being estimated. Normally that is
  # srr_fun, the curve the hindcast is fit with. Under the Ianelli
  # configuration (srr_fun < 2 with srr_pred_fun > 1) the hindcast is annual
  # deviates about mean recruitment and the curve enters as a penalty, but its
  # alpha is still estimated and can still fall below the replacement line --
  # so check srr_pred_fun there rather than skipping the fit.
  srr_fun      <- suppressWarnings(as.integer(dl[["srr_fun"]])[1])
  srr_pred_fun <- suppressWarnings(as.integer(dl[["srr_pred_fun"]])[1])
  ianelli <- is.na(srr_fun) || srr_fun < 2
  curve   <- if (ianelli) srr_pred_fun else srr_fun
  if (is.na(curve) || curve < 2) return(list())

  # Steepness needs spawning biomass per recruit, which is undefined under
  # predation, so the curve cannot be tested against the replacement line.
  # What can be tested without SPR is whether the data informed the curve.
  if (isTRUE(as.integer(dl[["msmMode"]] %||% 0L)[1] > 0L)) {
    return(.check_stock_recruit_msm(object, curve))
  }

  # Both are [nspp, nyrs]; the stock-recruit curve is summarised by its first
  # year, so read column 1.
  col1 <- function(x) {
    if (is.null(x)) return(numeric(0))
    if (is.matrix(x)) as.numeric(x[, 1]) else as.numeric(x)
  }
  h  <- suppressWarnings(col1(q[["steepness"]]))
  # Under Ianelli the reported R0 is the mean-recruitment level, not the value
  # the penalty curve implies, so it says nothing about that curve. Steepness
  # is the complete test in any case: for Beverton-Holt,
  # h = alpha*SPR0/(4 + alpha*SPR0), so alpha < 1/SPR0 is exactly h < 0.2.
  r0 <- if (ianelli) numeric(0) else suppressWarnings(col1(q[["R0"]]))
  if (!length(h) && !length(r0)) return(list())

  bad_r0 <- which(is.finite(r0) & r0 <= 0)
  bad_h  <- if (curve %in% c(2L, 3L)) which(is.finite(h) & h < 0.2) else integer(0)
  if (!length(bad_r0) && !length(bad_h)) {
    return(list(stock_recruit = .conv_record(
      "stock_recruit", "fit", "OK",
      sprintf("Stock-recruit curve is well posed (steepness %s).",
              paste(sprintf("%.3f", h), collapse = ", ")),
      list(steepness = h, R0 = r0))))
  }

  spp <- dl[["spnames"]]
  nm  <- function(i) if (length(spp) >= max(i)) paste(spp[i], collapse = ", ") else paste(i, collapse = ", ")
  msg <- character(0)
  if (length(bad_h)) {
    msg <- c(msg, sprintf(
      "Beverton-Holt steepness below 0.2 for %s (h = %s): alpha is under the replacement line 1/SPR0, so the stock cannot replace itself.",
      nm(bad_h), paste(sprintf("%.3f", h[bad_h]), collapse = ", ")))
  }
  if (length(bad_r0)) {
    msg <- c(msg, sprintf(
      "Implied unfished recruitment R0 is non-positive for %s (R0 = %s).",
      nm(bad_r0), paste(signif(r0[bad_r0], 4), collapse = ", ")))
  }
  msg <- c(msg, "Refit from stock-recruit starting values on the right scale -- see build_srr(srr_alpha_init =, srr_beta_init =).")

  list(stock_recruit = .conv_record(
    "stock_recruit", "fit", "FAIL", paste(msg, collapse = " "),
    list(steepness = h, R0 = r0)))
}


# Under predation steepness does not exist and nothing anchors alpha and beta,
# so an uninformed curve runs to a flat ridge (recruitment independent of SSB
# over the data) or a linear one (beta at 0), each with a positive-definite
# Hessian. Read the curve over the hindcast SSB range, the +/-30 bound and the
# standard errors instead.
.check_stock_recruit_msm <- function(object, curve) {
  dl <- object[["data_list"]]
  q  <- object[["quantities"]]
  rp <- object[["estimated_params"]][["rec_pars"]]
  ssb <- q[["ssb"]]
  if (is.null(rp) || is.null(ssb)) return(list())
  nyrs_hind <- dl[["endyr"]] - dl[["styr"]] + 1
  spp  <- dl[["spnames"]] %||% paste0("Species", seq_len(nrow(rp)))
  estd <- dl[["estDynamics"]] %||% rep(0, nrow(rp))
  ch   <- object[[".conv_hindcast"]]

  pos_of <- function(sp, slot) {
    if (is.null(ch$index)) return(integer(0))
    ch$index[ch$index$block == "rec_pars" & ch$index$slot == slot &
               ch$index$species == spp[sp], "par_index"]
  }
  # NA without an sdreport (getsd = FALSE), so the standard-error clause is
  # silent then.
  se_of <- function(sp, slot) {
    i <- pos_of(sp, slot)
    if (is.null(ch$se_fixed) || length(i) != 1L || i > length(ch$se_fixed)) return(NA_real_)
    unname(ch$se_fixed[i])
  }
  # A curve held at its inputs (srr_est_mode = "Fixed", or a linkage fixing
  # alpha and beta) is what the user asked for, not something the data failed
  # to inform; it is reported without a warning.
  fixed_curve <- function(sp) {
    !is.null(ch$index) && !length(pos_of(sp, "alpha")) && !length(pos_of(sp, "beta"))
  }

  msg <- character(0); dat <- list(); held <- character(0)
  for (sp in seq_len(nrow(rp))) {
    if (isTRUE(estd[sp] > 0)) next                       # input numbers, no curve
    if (fixed_curve(sp)) { held <- c(held, spp[sp]); next }
    la <- rp[sp, 2]; lb <- rp[sp, 3]
    a  <- exp(la);   b  <- exp(lb)
    s  <- ssb[sp, seq_len(nyrs_hind)]
    s  <- s[is.finite(s) & s > 0]
    if (!length(s)) next
    smin <- min(s); smax <- max(s)
    se_a <- se_of(sp, "alpha"); se_b <- se_of(sp, "beta")
    why <- character(0)
    # Density dependence at the two ends of the observed SSB range.
    if (curve %in% c(2L, 3L)) {
      # Beverton-Holt: predicted / asymptote = beta S / (1 + beta S).
      dd_smin <- b * smin / (1 + b * smin); dd_smax <- b * smax / (1 + b * smax)
      if (dd_smin > 0.9) why <- c(why, sprintf(
        "flat over the observed SSB range (predicted/asymptote %.3f at the lowest SSB, %.4g): recruitment does not depend on SSB in these data",
        dd_smin, smin))
      if (dd_smax < 0.1) why <- c(why, sprintf(
        "linear over the observed SSB range (predicted/asymptote %.3f at the highest SSB, %.4g): no density dependence, beta is not informed",
        dd_smax, smax))
    } else {
      # Ricker: R = alpha S exp(-beta S / 1e6); the density-dependence factor is
      # exp(-beta S / 1e6) and the curve peaks at S* = 1e6 / beta.
      dd_smin <- exp(-b * smin / 1e6); dd_smax <- exp(-b * smax / 1e6)
      speak   <- 1e6 / b
      if (dd_smax > 0.9) why <- c(why, sprintf(
        "linear over the observed SSB range (density-dependence factor %.3f at the highest SSB, %.4g): no density dependence, beta is not informed",
        dd_smax, smax))
      if (speak < smin) why <- c(why, sprintf(
        "peak at SSB %.4g, below the lowest observed SSB (%.4g): every observation is on the descending limb",
        speak, smin))
    }
    at_bound <- abs(c(la, lb)) >= 30 - 1e-3
    if (any(at_bound)) why <- c(why, sprintf(
      "%s at the +/-30 log-scale overflow bound",
      paste(c("alpha", "beta")[at_bound], collapse = " and ")))
    big_se <- c(alpha = se_a, beta = se_b)
    big_se <- big_se[is.finite(big_se) & big_se > 10 | is.nan(big_se)]
    if (length(big_se)) why <- c(why, sprintf(
      "log-scale standard error %s on %s",
      paste(signif(big_se, 3), collapse = " and "), paste(names(big_se), collapse = " and ")))
    dat[[spp[sp]]] <- list(alpha = a, beta = b, se_log_alpha = se_a, se_log_beta = se_b,
                           dd_at_smin = dd_smin, dd_at_smax = dd_smax,
                           ssb_range = c(smin, smax))
    if (length(why)) msg <- c(msg, sprintf("%s: %s.", spp[sp], paste(why, collapse = "; ")))
  }
  if (!length(dat)) {
    if (!length(held)) return(list(stock_recruit = .conv_record(
      "stock_recruit", "fit", "NOTE",
      "Stock-recruit curve not checked: no species carries an estimated curve with usable spawning biomass.",
      list())))
    return(list(stock_recruit = .conv_record(
      "stock_recruit", "fit", "NOTE",
      sprintf("Stock-recruit curve held at its inputs for %s; nothing to check.",
              paste(held, collapse = ", ")),
      list(held = held))))
  }
  if (!length(msg)) {
    return(list(stock_recruit = .conv_record(
      "stock_recruit", "fit", "OK",
      "Stock-recruit curve bends within the observed SSB range (steepness is undefined under msmMode > 0, so the replacement line is not checked).",
      dat)))
  }
  list(stock_recruit = .conv_record(
    "stock_recruit", "fit", "WARN",
    paste(c("Stock-recruit curve not informed by the data under predation.", msg,
            "Anchor it with an alpha prior through build_srr(linkages = ), or report the fit as a mean-recruitment model."),
          collapse = " "),
    dat))
}


#' Convergence diagnostics for a fitted Rceattle model
#'
#' Runs the post-fit convergence battery and returns a single structured
#' object. Each check yields a record with a common schema
#' (\code{id}, \code{tier}, \code{severity}, \code{message}, \code{data});
#' \code{severity} is one of \code{"OK"}, \code{"NOTE"}, \code{"WARN"},
#' \code{"FAIL"}. The object's \code{status} is the worst severity present.
#'
#' \code{fit_mod()} runs this automatically and attaches the result as
#' \code{fit$convergence}; call \code{convergence_diagnostics()} directly to
#' re-run it on any fit. Checks cover the optimizer gradient, a requested
#' \code{sdreport} that did not return, an \code{sdreport} that was never
#' requested, Hessian positive-definiteness and conditioning, parameters on
#' bounds, a deviation variance estimated to zero, phasing, parameter
#' estimability, a numbers-at-age, Ricker-intercept or recruitment floor that
#' was reached, and the stock-recruit curve.
#'
#' @param object An object of class \code{"Rceattle"} returned by [fit_mod()].
#' @param ... Currently unused.
#'
#' @return An object of class \code{"Rceattle_convergence"}: a list with
#'   \code{status} (overall worst severity) and \code{checks} (named list of
#'   records).
#' @export
convergence_diagnostics <- function(object, ...) {
  # Four checks name their parameters by coordinate. The hindcast index was
  # stored by fit_mod(); the one for the fit's final parameter vector means
  # pushing a tagged vector through TMB's parList(), so a promise builds the pair
  # at most once, and not at all on a fit where every check passes.
  delayedAssign("index", .conv_par_index(object))
  checks <- c(
    .check_phasing(object),
    .check_optimizer(object, index),
    .check_sdreport_failed(object),
    .check_hessian_not_run(object),
    .check_hessian_eigen(object, index),
    .check_bounds(object, index),
    .check_variance_collapse(object),
    .check_estimability_record(object, index),
    .check_zero_n_penalty(object),
    .check_stock_recruit(object)
  )
  structure(
    list(status = .conv_overall(checks), checks = checks),
    class = "Rceattle_convergence"
  )
}


#' @export
print.Rceattle_convergence <- function(x, all = FALSE, ...) {
  cat("<Rceattle convergence>  status:", x$status, "\n")
  if (length(x$checks) == 0) {
    cat("  (no checks run)\n")
    return(invisible(x))
  }
  show <- if (isTRUE(all)) x$checks else
    Filter(function(c) c$severity != "OK", x$checks)
  if (length(show) == 0) {
    cat("  all checks OK\n")
    return(invisible(x))
  }
  for (c in show) {
    cat(sprintf("  [%-4s] %-16s %s\n", c$severity, c$id, c$message))
  }
  invisible(x)
}
