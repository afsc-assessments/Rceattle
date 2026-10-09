#' Bundle the optimizer / sdreport / phasing controls for `fit_mod()`
#'
#' Collects the optimizer and reporting settings into one object, so a call to
#' [fit_mod()] stays about the model rather than how it is fit. Pass the result
#' through `fit_mod()`'s `fit_control` argument, where it overrides the
#' equivalent individual arguments, which are kept for back-compatibility.
#'
#' @details
#' # Selectivity standard errors
#'
#' `selectivity_se` needs `getsd = TRUE`, and reports the error belonging to
#' whichever `sdreport` the fit ends on. Under `estimateMode = "Estimate"` with
#' an HCR that re-optimizes, that is the *projection* fit, where every
#' selectivity parameter is mapped off, so every error comes back exactly 0;
#' `"Projection"` estimates no selectivity at all. Use `"Hindcast"`, or
#' `"Estimate"` with `projection_uncertainty = TRUE`, to get an error from a fit
#' that estimated the curve. [fit_mod()] warns in each of these cases.
#'
#' It is off by default because the delta method forms a Jacobian of every
#' reported value against every parameter, so the cost is their product: on
#' `Atka2022`, 1,012 values against 584 parameters. The error is on the log
#' scale rather than the logit, because the non-parametric forms normalize to
#' mean selectivity 1 instead of a maximum of 1 -- 58% of `Atka2022`'s
#' `sel_at_age` exceeds 1, up to 3.06 -- so a logit is undefined over most of
#' the array.
#'
#' Rows cover estimated, age-based lead fleets only, starting at each fleet's
#' first selected bin. Four kinds of cell are a structural zero (a `Fixed`
#' fleet's empirical curve, a length-based fleet's growth-matrix projection, a
#' bin below `Bin_first_selected`, and array padding), and one `log(0) = -Inf`
#' on the tape turns *every* quantity in the `sdreport` to `NaN`, biomass and
#' SSB included. All four are identified from the data, so the reported set
#' never depends on a parameter value, and a value that underflows to zero is
#' floored so it cannot reintroduce the `-Inf`. See [plot_selectivity()], which
#' draws the interval.
#'
#' @param bias.correct logical. Apply bias correction via [TMB::sdreport()].
#'   Default `FALSE`.
#' @param getsd logical. Run [TMB::sdreport()] after optimization. Default
#'   `TRUE`.
#' @param getJointPrecision logical. Return the full Hessian of fixed and
#'   random effects. Default `TRUE`.
#' @param getReportCovariance logical. Return the variance-covariance of
#'   `ADREPORT` variables. Default `FALSE`.
#' @param projection_uncertainty logical. Carry hindcast parameter uncertainty
#'   into an HCR projection, by refitting with the hindcast and
#'   biological-reference-point parameters turned on. Default `FALSE` for speed.
#' @param selectivity_se logical. Also return a standard error for log
#'   selectivity-at-age. Default `FALSE`; see Details.
#' @param comp_offset Numeric or `NULL`. Small proportion added to the observed
#'   and predicted composition and CAAL bins before the multinomial /
#'   Dirichlet-multinomial likelihood, so an empty bin is not `log(0)`. Stored
#'   on `data_list`, so the fitted likelihood, the OSA observation vector and
#'   any internal refit share it. `NULL` (default) inherits
#'   `data_list$comp_offset` if set, else `1e-5`. Does not apply to diet
#'   compositions.
#' @param bias_adjust_obs logical, default `TRUE`. Apply a bias adjustment
#'   (mean - sd^2/2) to lognormal data likelihoods.
#' @param bias_adjust_proc logical, default `TRUE`. Shift lognormal process
#'   likelihoods, lognormal priors and the Ianelli stock-recruit penalty by
#'   `-sd^2/2` on the log scale, making each prior value or curve a mean rather
#'   than a median. A value between 0 and 1 scales the shift.
#' @param use_gradient logical. Use the analytic gradient during phasing.
#'   Default `TRUE`.
#' @param rel_tol Numeric tolerance for flagging a discontinuous likelihood,
#'   comparing `nlminb`'s objective with a fresh evaluation of the object it
#'   came from. Default `1`.
#' @param loopnum Integer. Times to restart optimization; `3` sometimes reaches
#'   a lower final gradient than `1`. Default `5`.
#' @param newtonsteps Integer. Extra unconstrained Newton steps after
#'   optimization, an alternative to `loopnum`. Default `0`.
#' @param phase `TRUE`/`FALSE`, or a list of parameter names with their phases.
#'   Default `FALSE`.
#' @param TMBfilename Optional path (without `.cpp`) to an alternate TMB
#'   template for development. Default `NULL` uses the bundled `ceattle`.
#' @param verbose `0` silent, `1` model-fit updates, `2` model-fit and TMB
#'   progress. Default `1`.
#' @param nlminb_control Control parameters passed to [stats::nlminb()]; see
#'   `?nlminb`. Default `list(eval.max = 1e9, iter.max = 1e9, trace = 0)`.
#'
#' @return A list of class `"Rceattle_fit_control"`.
#' @export
#'
#' @examples
#' # Quick-and-dirty fit: skip sdreport, single optimizer pass
#' ctl <- fit_control(getsd = FALSE, loopnum = 1)
#'
#' # Production fit with bias correction and joint precision
#' ctl <- fit_control(bias.correct = TRUE, getJointPrecision = TRUE)
fit_control <- function(
  bias.correct        = FALSE,
  getsd               = TRUE,
  getJointPrecision   = TRUE,
  getReportCovariance = FALSE,
  projection_uncertainty = FALSE,
  selectivity_se      = FALSE,
  comp_offset         = NULL,
  bias_adjust_obs     = TRUE,
  bias_adjust_proc    = TRUE,
  use_gradient        = TRUE,
  rel_tol             = 1,
  loopnum             = 5,
  newtonsteps         = 0,
  phase               = FALSE,
  TMBfilename         = NULL,
  verbose             = 1,
  nlminb_control      = list(eval.max = 1e+09, iter.max = 1e+09, trace = 0)
) {
  ans <- list(
    bias.correct        = bias.correct,
    getsd               = getsd,
    getJointPrecision   = getJointPrecision,
    getReportCovariance = getReportCovariance,
    projection_uncertainty = projection_uncertainty,
    selectivity_se      = selectivity_se,
    comp_offset         = comp_offset,
    bias_adjust_obs     = bias_adjust_obs,
    bias_adjust_proc    = bias_adjust_proc,
    use_gradient        = use_gradient,
    rel_tol             = rel_tol,
    loopnum             = loopnum,
    newtonsteps         = newtonsteps,
    phase               = phase,
    TMBfilename         = TMBfilename,
    verbose             = verbose,
    nlminb_control      = nlminb_control
  )
  class(ans) <- c("Rceattle_fit_control", "list")

  # Which fields the caller actually named. The returned list holds every field,
  # defaults included, so the value alone cannot say whether it was asked for --
  # and `fit_control(getsd = TRUE)` is a request even though TRUE is the default.
  # `match.call()` resolves positional and partially-matched arguments to their
  # full names first. The `$<-`/`[[<-`/`[<-` methods below keep it current when
  # the bundle is edited afterwards. Read it with `.rce_fit_control_supplied()`.
  attr(ans, "supplied") <- names(as.list(match.call())[-1])
  ans
}


# Editing a bundle after the fact -- `ctl$getsd <- TRUE`, or the `[[<-` that
# modifyList() uses -- is a request just as much as naming the field in the
# fit_control() call, and the value cannot say so when it equals the default.
# These record it. Nothing else about assignment changes.
.rce_record_supplied <- function(x, y, i) {
  nm <- if (is.character(i)) i else names(y)[i]
  nm <- nm[!is.na(nm) & nzchar(nm)]
  attr(y, "supplied") <- union(attr(x, "supplied"), nm)
  y
}

#' @export
#' @noRd
`$<-.Rceattle_fit_control` <- function(x, name, value) {
  .rce_record_supplied(x, NextMethod(), name)
}

#' @export
#' @noRd
`[[<-.Rceattle_fit_control` <- function(x, i, value) {
  .rce_record_supplied(x, NextMethod(), i)
}

#' @export
#' @noRd
`[<-.Rceattle_fit_control` <- function(x, i, value) {
  .rce_record_supplied(x, NextMethod(), i)
}


#' Which `fit_control()` fields the caller asked for
#'
#' @description
#' A field counts as asked for if the caller named it in the `fit_control()`
#' call or assigned to it afterwards. Restricted to fields the bundle still
#' holds, so deleting one (`ctl$getsd <- NULL`) withdraws the request.
#'
#' A value that no longer matches the default also counts, whatever the record
#' says. That is the fallback for a bundle whose record did not survive, one
#' rebuilt by `load_config()` from the non-default fields, or by `structure()`.
#'
#' @param fit_control an `Rceattle_fit_control` object.
#' @return Character vector of field names.
#' @keywords internal
#' @noRd
.rce_fit_control_supplied <- function(fit_control) {
  defaults <- fit_control()
  nms      <- intersect(names(fit_control), names(defaults))
  edited   <- nms[!vapply(nms, function(n) identical(fit_control[[n]], defaults[[n]]),
                          logical(1))]
  union(intersect(attr(fit_control, "supplied"), nms), edited)
}


#' Read the refit knobs a diagnostic can actually honour
#'
#' @description
#' `retrospective()`, `jitter()` and `self_test()` refit through `.refit_like()`,
#' which forwards `phase` and `getsd` from its caller. The rest of
#' `fit_control()` either comes off the source `data_list` (`loopnum`,
#' `bias_adjust_obs`, `bias_adjust_proc`, `comp_offset`) or falls back to
#' `.refit_like()`'s own defaults (`newtonsteps`, `nlminb_control`, `rel_tol`,
#' `use_gradient`, ...). Accepting a whole `fit_control()` and quietly dropping
#' those would leave a user believing every peel had used them.
#'
#' So a field this path cannot reach is an error, not a silent no-op, but only
#' when the caller changed it from `fit_control()`'s default. `.refit_like()`
#' builds a fresh `fit_control()` naming only `phase`, `loopnum`, `getsd`,
#' `verbose` and the bias-adjustment flags, so an unreachable field left at its
#' default is what the refit uses anyway, and refusing it would be pedantry.
#'
#' **`phase` and `getsd` are read from what the caller SET**, via
#' `.rce_fit_control_supplied()`, named in the `fit_control()` call or
#' assigned to afterwards, not from whether the value differs from a
#' default. The two disagree: `fit_control()` defaults `phase` to `FALSE`, but
#' `retrospective()` phases its peels because otherwise the parameters do not
#' move from the full-model fit and Mohn's rho is biased towards zero. Keying on
#' the value would mean `fit_control(getsd = FALSE)`, a request about standard
#' errors, silently un-phasing every peel, while `fit_control(getsd = TRUE)`
#' would do nothing at all. Keying on the name gives each field exactly what was
#' asked for and nothing else.
#'
#' @param fit_control a `fit_control()` object, or `NULL`.
#' @param fn calling function, for the message.
#'
#' @return A named list of the honoured knobs the caller asked for. Empty if
#'   they asked for none; `NULL` if `fit_control` was `NULL`.
#' @keywords internal
#' @noRd
.rce_refit_control <- function(fit_control, fn) {
  if (is.null(fit_control)) return(NULL)
  if (!inherits(fit_control, "Rceattle_fit_control")) {
    stop(sprintf("%s(): `fit_control` must come from fit_control().", fn),
         call. = FALSE)
  }

  honoured <- c("phase", "getsd")
  defaults <- fit_control()
  changed  <- function(nms) {
    nms[!vapply(nms, function(n) identical(fit_control[[n]], defaults[[n]]),
                logical(1))]
  }

  unreachable <- changed(setdiff(names(fit_control), honoured))
  if (length(unreachable)) {
    stop(sprintf(
      "%s() refits through .refit_like(), which takes only %s from fit_control(). It cannot apply %s. Set %s in the fit_mod() call that produced the model -- note that .refit_like() recovers only `loopnum`, `comp_offset` and the bias-adjustment flags from it, and refits under the defaults for the rest.",
      fn, paste(sprintf("`%s`", honoured), collapse = " and "),
      paste(sprintf("`%s`", unreachable), collapse = ", "),
      if (length(unreachable) == 1) "it" else "them"), call. = FALSE)
  }

  fit_control[intersect(honoured, .rce_fit_control_supplied(fit_control))]
}
