#' Iteratively reweight composition data (McAllister-Ianelli)
#'
#' @description
#' Tunes multinomial composition weights by refitting: each pass sets a fleet's
#' weight to the McAllister & Ianelli (1997) weight implied by the previous fit,
#' and stops once the weights settle. The weight multiplies the input sample
#' size, so changing it changes the fit and hence the next implied weight,
#' which is why the weights are tuned iteratively rather than in one step.
#'
#' @details
#' Only multinomial fleets are tuned. A Dirichlet-multinomial fleet estimates
#' its own weight within the likelihood, so tuning it externally would compete
#' with that estimate; those fleets are named and skipped, as are any requested
#' through `fleets` that have no fitted composition data.
#'
#' Only the age/length composition weight (`Comp_weights`) is tuned. The
#' conditional age-at-length and diet weights have McAllister-Ianelli analogues
#' of their own (`CAAL_weights`, `Diet_comp_weights`), but they are not tuned
#' here; set them by hand if needed.
#'
#' Tuning stops when the largest relative change in any weight falls below
#' `tol`. Hitting `n_iter` first gives a warning, not an error, the partly
#' tuned fit is returned either way, and `$reweight$history` shows whether the
#' weights were still settling.
#'
#' The weight is a parameter, not a data input, so each pass holds it to the
#' next through `inits`; `fleet_control$Comp_weights` supplies the value only
#' when a model is built from scratch, and editing that column then refitting
#' from an existing fit has no effect (`fit_mod()` warns when it detects this).
#' Each pass writes the new weight to both, so the returned model's data list
#' rebuilds the tuned model from scratch.
#'
#' @param object An `Rceattle` model fitted at `estimateMode` `"Estimate"` (0) or `"Hindcast"` (1), the modes that optimize the hindcast.
#' @param fit deprecated name for `object`, still accepted so existing
#'   scripts keep working. Supplying both is an error.
#' @param n_iter Maximum number of iterations (default 10).
#' @param tol Relative change in the weights below which to stop (default 0.01).
#' @param fleets Fleets to tune, as `Fleet_code` values or as names matching
#'   `fleet_control$Fleet_name` (as in [linkage_spec()]). Give ids or names, not
#'   a mix, R coerces `c(1, "Survey")` to `c("1", "Survey")`. Default `NULL`
#'   tunes every eligible fleet.
#' @param verbose Print the weights each iteration (default `TRUE`).
#'
#' @return The model refitted with the final weights, holding a `reweight`
#'   element: `history` (one row per fleet per iteration), `iterations`, and
#'   `converged`. The last rows of `history` are the weights the returned model
#'   was fitted with.
#'
#' @references
#' McAllister, M.K. and Ianelli, J.N. (1997) Bayesian stock assessment using
#' catch-age data and the sampling-importance resampling algorithm.
#' \emph{Canadian Journal of Fisheries and Aquatic Sciences} 54:284-300.
#'
#' Francis, R.I.C.C. (2011) Data weighting in statistical fisheries stock
#' assessment models. \emph{Canadian Journal of Fisheries and Aquatic Sciences}
#' 68:1124-1138. (The main alternative; applied by hand rather than by this
#' function.)
#'
#' @examples
#' \dontrun{
#' fit  <- fit_mod(BS2017SS, estimateMode = "Hindcast")
#' tuned <- reweight_comps(fit)
#' tuned$reweight$history      # weight per fleet, per iteration
#' }
#' @export
reweight_comps <- function(object = NULL, n_iter = 10, tol = 0.01, fleets = NULL,
                           verbose = TRUE, fit = NULL) {
  # `fit` was the old name for `object`; see R/0-deprecate.R.
  if (!missing(fit))
    object <- .rce_deprecated_arg(fit, !missing(object), "fit", "object", "reweight_comps")


  if (!inherits(object, "Rceattle")) {
    stop("`object` must be a fitted Rceattle model.", call. = FALSE)
  }
  fc <- object$data_list$fleet_control

  # Fleets the loop can tune: they carry composition data that is actually
  # fitted, and their likelihood takes the weight as a multiplier rather than
  # estimating it. A negative year marks a row carried but not fitted, and the
  # McAllister-Ianelli weight averages over fitted rows alone, so a fleet with
  # none of them has no weight to tune towards.
  cd <- object$data_list$comp_data
  has_comp <- unique(cd$Fleet_code[cd$Fleet_code > 0 & cd$Year > 0])
  dist <- as.character(fc$Comp_distribution)
  multinom <- fc$Fleet_code[dist %in% c("Multinomial", "MultinomialAFSC")]
  eligible <- intersect(has_comp, multinom)

  skipped <- setdiff(has_comp, eligible)
  if (length(skipped) && verbose) {
    message("Not tuning fleet(s) ", paste(skipped, collapse = ", "),
            ": their composition likelihood estimates its own weight.")
  }
  if (!is.null(fleets)) {
    # Accept names as well as codes, matching linkage_spec(fleet = ). The shared
    # resolver reports an unknown name against the model's own fleet list; left
    # unresolved, a name would fall through as "not eligible" and surface as the
    # misleading "no fitted composition data" warning below.
    requested <- fleets
    if (is.character(fleets)) {
      pos <- .resolve_stratum_names(fleets, fc$Fleet_name, "fleets",
                                    "fleet names (fleet_control$Fleet_name)")
      # The resolver returns row positions; this function works in Fleet_code.
      fleets <- fc$Fleet_code[pos]
    }
    # A fleet the caller named but cannot be tuned is reported: quietly tuning
    # the remainder would return a weighting that was never asked for. Report it
    # the way the caller wrote it.
    drop_idx <- !fleets %in% eligible
    if (any(drop_idx)) {
      warning("Not tuning requested fleet(s) ",
              paste(requested[drop_idx], collapse = ", "),
              ": they have no fitted composition data, or their likelihood ",
              "estimates its own weight.", call. = FALSE)
    }
    eligible <- intersect(eligible, fleets)
  }
  if (!length(eligible)) {
    stop("No fleets to reweight: none have composition data fitted with a ",
         "multinomial likelihood.", call. = FALSE)
  }

  estimate_mode <- object$data_list$estimateMode
  if (!is.null(estimate_mode) && estimate_mode > 1) {
    stop("`object` was run at estimateMode ", estimate_mode, ", which does not ",
         "optimize the hindcast, so the weights cannot be tuned: every pass ",
         "would return the same value and report convergence. Refit at ",
         "\"Estimate\" or \"Hindcast\" first.", call. = FALSE)
  }
  # Refit the way the input fit was produced. Left to .refit_like()'s defaults
  # (map = NULL, phase = FALSE, getsd = FALSE) the loop would rebuild a
  # user-supplied map, drop phasing the model may need to reach its optimum, and
  # return the tuned fit without standard errors. `phase` may be a named list of
  # per-parameter phases, so it travels verbatim rather than as a flag.
  fc_cfg     <- object$run_config$fit_control
  fit_phase  <- if (is.null(fc_cfg$phase)) FALSE else fc_cfg$phase
  fit_getsd  <- if (is.null(fc_cfg$getsd)) FALSE else fc_cfg$getsd
  fit_map    <- object$map

  history <- NULL
  converged <- FALSE
  iter_run <- 0L

  for (iter in seq_len(n_iter)) {
    iter_run <- iter

    # Weight implied by the current fit, per fleet.
    implied <- object$data_list$fleet_control$Comp_weights_mcallister[
      match(eligible, object$data_list$fleet_control$Fleet_code)]
    if (all(is.na(implied))) {
      stop("The fitted model carries no McAllister-Ianelli weights; it must be fitted ",
           "with composition data before it can be reweighted.", call. = FALSE)
    }
    # The weight has to be finite and positive to be usable: zero drops the
    # fleet's composition data from the likelihood, and a NaN carries into every
    # later warm start, leaving the loop unable to converge.
    bad <- !is.finite(implied) | implied <= 0
    if (any(bad)) {
      stop("McAllister-Ianelli weight is ",
           paste(unique(implied[bad]), collapse = ", "), " for fleet(s) ",
           paste(eligible[bad], collapse = ", "),
           "; it cannot be applied. Check that those fleets have fitted ",
           "composition data with non-zero sample sizes.", call. = FALSE)
    }

    current <- object$estimated_params$comp_weights[eligible]
    # Relative move, guarding a zero denominator.
    rel <- abs(implied - current) / pmax(abs(current), .Machine$double.eps)
    history <- rbind(history, data.frame(
      iteration = iter, Fleet_code = eligible, weight = implied,
      rel_change = rel, row.names = NULL))

    if (verbose) {
      message("Reweight ", iter, ": ",
              paste0("fleet ", eligible, " = ", signif(implied, 4),
                     collapse = "; "))
    }

    converged <- all(is.finite(rel)) && max(rel) < tol

    # Apply and refit even on the final pass, so the model returned is the one
    # fitted with the converged weights rather than the one that implied them.
    # The weight the fit reads is the parameter, so it travels in `inits`; the
    # column is moved with it so the two never disagree. Keeping them in step
    # every pass -- rather than only writing the column back at the end -- means
    # a data_list taken from any iteration rebuilds that iteration's model, and
    # it keeps fit_mod()'s "this column is being ignored" warning quiet, which
    # here would fire once per fleet per pass on a weighting the loop is setting
    # deliberately.
    inits <- object$estimated_params
    inits$comp_weights[eligible] <- implied
    data_list <- object$data_list
    data_list$fleet_control$Comp_weights[eligible] <- implied

    object <- .refit_like(data_list = data_list, inits = inits,
                       map = fit_map, estimateMode = estimate_mode,
                       phase = fit_phase, getsd = fit_getsd)

    if (converged) break
  }

  if (!converged) {
    warning("Composition weights had not converged to within tol = ", tol,
            " after ", iter_run, " iteration(s); returning the last fit. ",
            "Inspect `$reweight$history`.", call. = FALSE)
  }

  # The loop carries the column along with the parameter, so the two already
  # agree; re-assert it from the returned fit rather than relying on that, since
  # this is the guarantee the docs make -- saving this data_list and refitting
  # from scratch reproduces the tuned model.
  object$data_list$fleet_control$Comp_weights[eligible] <-
    object$estimated_params$comp_weights[eligible]

  object$reweight <- list(history = history, iterations = iter_run,
                       converged = converged, tol = tol, fleets = eligible)
  object
}
