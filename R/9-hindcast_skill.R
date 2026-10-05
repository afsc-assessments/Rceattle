#' Forecast skill across retrospective peels
#'
#' Scores how well each retrospective peel's projection recovers the full
#' time-series model's estimate over the years the peel did not see, given the
#' known catch. The full model is taken as the best available estimate of what
#' happened, so this measures forecast skill rather than the estimation
#' consistency [retrospective()] reports through Mohn's rho.
#'
#' @details
#' Each peel drops the last `i` years of data, refits, and then solves fishing
#' mortality against the catch that was actually taken in those years -- so the
#' projection differs from the full model only in what it could not know:
#' the recruitment, and anything downstream of it. That makes this the natural
#' comparison for recruitment-projection assumptions, e.g.
#' `proj_mean_rec = TRUE` against `FALSE` against a DSEM.
#'
#' Skill is reported as the mean absolute scaled error (MASE), following the
#' hindcast cross-validation of Kell et al. (2016) but scored against the full
#' model's estimate rather than against held-out observations:
#'
#' \deqn{MASE = \frac{\overline{|forecast - reference|}}{\overline{|naive - reference|}}}
#'
#' The naive forecast is the peel's own terminal estimate carried forward,
#' i.e. "assume nothing changes". So:
#' \itemize{
#'   \item `MASE < 1` -- the projection beats assuming no change.
#'   \item `MASE > 1` -- it is worse than assuming no change, which usually means
#'     the recruitment assumption is pulling the projection the wrong way.
#'   \item `MASE = 1` -- no better than persistence.
#' }
#'
#' A MASE is undefined when the naive forecast happens to be exactly right
#' (a zero denominator); those rows come back `NA` rather than `Inf`.
#'
#' One row per horizon per species per quantity: the mean runs ACROSS PEELS at a
#' fixed steps-ahead, which is Kell et al. (2021) eq. 5, whose sums run over
#' \eqn{t = T-n \ldots T} at fixed \eqn{h}. `N` is how many peels
#' contributed, so `peels` buys terms inside each MASE rather than more rows.
#'
#' **Know what the baseline is under `reference = "model"`.** It is NOT a
#' persistence error of one series. `naive` is the PEEL's terminal estimate and
#' `reference` is the FULL model's, so each denominator term is a difference
#' between two different models, carrying the peel's own retrospective bias at
#' its terminal year. A peel whose terminal estimate happens to sit near the
#' full model's later value therefore contributes a tiny term: measured on
#' BS2017SS pollock SSB with `peels = c(2,3,5)`, the three one-year denominator
#' terms spanned 18,254 to 1,965,189 -- 108-fold -- because one peel's 2015 SSB
#' landed within 0.3% of the full model's 2016 value.
#'
#' Averaging across peels does tame the ratio: those same peels give per-peel
#' ratios of 1.02, 2.82 and 34.22, and a reported MASE of 1.485. That is the
#' point of scoring per horizon rather than per peel. But a horizon resting on
#' few peels inherits whatever bias cancellation they happen to carry, which is
#' a reason to read `N`, not a reason to distrust short horizons as such.
#'
#' Read `N` before reading a MASE, per row -- it is not constant down the
#' column. A horizon `h` can only be scored by a peel at least `h` years deep, so
#' on an ANNUAL series with a contiguous `peels` run it falls as the horizon
#' grows: with `peels = 2:4`, three peels at one year ahead and one at four. A
#' MASE resting on a single peel says as much about that peel as about the
#' projection.
#'
#' That pattern holds for an ANNUAL series. Under `reference = "observed"` an
#' irregular survey breaks it, because a fleet is scored at a horizon only where
#' it has an observation both at the peel's terminal year and after it -- so
#' `N` can rise with the horizon, and some horizons carry no rows at all.
#' On BS2017SS's acoustic pollock survey with `peels = c(2,3,5)`, only
#' `Forecast year` 2 and 4 are populated, with `N` 2 and 1.
#'
#' `$mase` is therefore ragged rather than rectangular: a horizon with no rows is
#' MISSING, not `NA`. Join two models' tables on `Forecast year` rather than
#' comparing them positionally, or a two-year row lines up against a four-year
#' one.
#'
#' The naive forecast differs slightly from Kell's by construction. Theirs is
#' \eqn{y_{t-h}}, the last OBSERVED value, which `reference = "observed"` uses
#' directly. Under `reference = "model"` there is no observation to carry
#' forward, so the baseline is the peel's own terminal estimate held flat -- the
#' model-based analogue of the same idea. Kell also averages the denominator over
#' a wider window (\eqn{n+1+h} terms against the numerator's \eqn{n+1}), because
#' \eqn{y_{t-h}} exists for more years; here both use the peels actually fitted.
#'
#' @param object A fitted Rceattle model (the full time series).
#' @param peels Which peels to fit. Passed to [retrospective()], so a single
#'   number `n` means every peel from 1 to `n`, and a vector names the depths
#'   exactly. Horizons run 1 to the DEEPEST depth asked for, so a deeper `peels`
#'   adds rows while more peels at the same depth only add terms to the mean
#'   inside each row -- `peels = 3:10` and `1:10` produce the same ten horizons,
#'   the second with more peels behind the short ones and two more model fits.
#'   Ignored when `retro` is supplied.
#' @param quantity Quantities to score against the full model. Any of `"ssb"`,
#'   `"biomass"`, `"R"`. Ignored when `reference = "observed"`.
#' @param reference What to score the projection against.
#'   `"model"` (default) compares each peel's projection to the full time-series
#'   model's estimate, taking that as the best available account of what
#'   happened. `"observed"` is classic hindcast cross-validation: it compares the
#'   peel's PREDICTED SURVEY INDEX to the index values actually observed in the
#'   held-out years. Both can be asked for. `"model"` scores the quantity a
#'   projection is used for (SSB, recruitment); `"observed"` scores the only
#'   thing that was really measured.
#'
#'   **`"observed"` is not numerically comparable to published MASE values.**
#'   Kell et al. (2021) eq. 5 averages its denominator over a WIDER window than
#'   its numerator -- \eqn{n+1+h} terms against \eqn{n+1}, because
#'   \eqn{y_{t-h}} exists for years the numerator cannot score -- and this uses
#'   the peels actually fitted for both. Computed from BS2017SS's index data,
#'   the two denominators differ by up to 2.6x, which is enough to move a MASE
#'   across 1 and so to flip the "beats persistence" verdict. The ordering of
#'   two models scored the same way here is unaffected; the absolute level is
#'   not the published statistic. `"model"` does not have this issue, since its
#'   baseline is defined from the peels themselves.
#' @param forecast_rec How the peeled years get their recruitment.
#'   `"model"` (default here) uses the model's own projection rule, in
#'   precedence order: `proj_mean_rec = TRUE` projects at mean recruitment,
#'   whatever process the model carries, because that is what the setting means;
#'   otherwise the LATENT STATES supply it where the deviations are random
#'   effects (`random_rec = TRUE`, or a DSEM), so an AR1's autocorrelation or a
#'   DSEM's lagged and covariate paths propagate into the forecast; otherwise
#'   recruitment comes off the stock-recruit curve. `"mean"` forces the
#'   historical mean for all of them, which is [retrospective()]'s default and
#'   the convention Mohn's rho is computed under. Applies only when this
#'   function fits the peels itself: a `retro` handed in already carries its
#'   forecast years, so supplying both is an error rather than a setting that
#'   does nothing.
#'
#'   The default differs from [retrospective()]'s on purpose. A peel's forecast
#'   years are hindcast years, so `proj_mean_rec` cannot reach them by itself,
#'   and with `"mean"` every model forecasts identically -- which makes this
#'   function unable to tell projection methods apart, the one thing it exists
#'   to do. Measured on the GOA arrowtooth model, `proj_mean_rec = TRUE` and
#'   `FALSE` returned byte-identical MASE under `"mean"`.
#' @param retro Optionally an already-computed [retrospective()] result, to
#'   avoid refitting when both are wanted. Its peels carry the forecast rule
#'   they were fitted under -- `forecast_rec = "mean"` with [retrospective()]'s
#'   defaults -- and that rule is what gets scored, so passing a conflicting
#'   `forecast_rec` alongside it is an error. To compare the two rules, fit one
#'   retrospective under each.
#' @param ... Passed to [retrospective()] (`cores`, `getsd`, `rescale`).
#'
#' @return A list with
#'   \describe{
#'     \item{`mase`}{one row per horizon x species x quantity. Keyed by
#'       `Object`, `Forecast year` and `species`, the column vocabulary
#'       [retrospective()]'s `$mohns` uses, so the two merge on those three;
#'       then `mase`, `mae_forecast`, `mae_naive`, `N` -- how many peels were
#'       averaged at that horizon -- and `peels_used`, which ones. A peel that
#'       could not be scored at a horizon (an NA in its forecast years, or a
#'       species with `estDynamics > 0`, whose numbers-at-age are input and so
#'       was never forecast) is left out of both rather than nulling the row.}
#'     \item{`by_year`}{the underlying series, in the same vocabulary:
#'       `Object`, `Forecast year`, `species`, `peel`, `year`, `forecast`,
#'       `reference`, `naive`. `year` is the calendar year; `Forecast year` is
#'       the horizon.}
#'   }
#'
#' @references
#' Kell, L.T., Sharma, R., Kitakado, T., Winker, H., Mosqueira, I., Cardinale,
#' M., Fu, D. (2021) Validation of stock assessment methods: is it me or my
#' model talking? \emph{ICES Journal of Marine Science} 78(6), 2244-2255.
#' \doi{10.1093/icesjms/fsab104}. Equation 5 is the MASE this function reports;
#' equation 4 is the naive baseline.
#'
#' Kell, L.T., Kimoto, A., Kitakado, T. (2016) Evaluation of the prediction
#' skill of stock assessment using hindcasting. \emph{Fisheries Research} 183,
#' 119-127.
#'
#' @export
hindcast_skill <- function(object = NULL, peels = 5,
                           quantity = c("ssb", "biomass", "R"),
                           reference = c("model", "observed"),
                           forecast_rec = c("model", "mean"),
                           retro = NULL, ...) {
  if (!inherits(object, "Rceattle")) {
    stop("`object` must be a fitted Rceattle model (from fit_mod()).",
         call. = FALSE)
  }
  quantity     <- match.arg(quantity, several.ok = TRUE)
  # match.arg(several.ok = TRUE) returns EVERY choice when the argument is left
  # at its default, so `reference` defaulted to both -- not to "model" as
  # documented. That is not a harmless extra: "observed" rebuilds every peel at
  # estimateMode = 3, and it refuses outright on a model with analytical
  # catchability, so the documented default call failed on models the function
  # can serve. Take the first when the caller did not ask.
  reference    <- if (missing(reference)) "model" else
    match.arg(reference, several.ok = TRUE)
  # Before match.arg(), which ASSIGNS to the formal and so makes missing() FALSE
  # from here on, whether or not the caller passed anything.
  .fr_supplied <- !missing(forecast_rec)
  forecast_rec <- match.arg(forecast_rec)

  # Only where `quantity` is actually read. Under reference = "observed" the
  # score is the predicted survey index, so refusing a fit that does not report
  # the default `biomass` would reject a model this function can serve.
  if ("model" %in% reference) {
    for (q in quantity) {
      if (is.null(object$quantities[[q]])) {
        stop("`quantity` '", q, "' is not reported by this fit.", call. = FALSE)
      }
    }
  }

  if (is.null(retro)) {
    retro <- retrospective(object, peels = peels,
                           forecast_rec = forecast_rec, ...)
  } else if (.fr_supplied) {
    # A supplied `retro` already holds its forecast years, computed under the
    # rule its own call was given, and nothing here recomputes them. So
    # `forecast_rec` cannot apply to it, and scoring the same peels twice under
    # two values of this argument would return the SAME numbers both times --
    # which reads as "the projection rule makes no difference" rather than as an
    # ignored argument. Refuse instead: the comparison needs one retrospective
    # per rule, because the rule is fixed when the peels are fitted.
    # `$forecast_rec` on the object, not `$mase$forecast_rec`: $mase is NULL when
    # no peel converged, and absent on a retro saved before it existed, which is
    # where ignoring the argument would do the most damage.
    .retro_rec <- retro[["forecast_rec"]] %||% unique(retro[["mase"]][["forecast_rec"]])
    if (length(.retro_rec) == 1L && !identical(.retro_rec, forecast_rec)) {
      stop("`forecast_rec` cannot be applied to a `retro` that was already ",
           "fitted: these peels carry their forecast under forecast_rec = \"",
           .retro_rec, "\", and nothing here recomputes it. To compare the two ",
           "rules, fit one retrospective under each -- ",
           "retrospective(object, peels = , forecast_rec = \"mean\") and the ",
           "same call with \"model\" -- then score each. See ?hindcast_skill.",
           call. = FALSE)
    }
  } else {
    # Not an error: reusing a retro WITHOUT naming forecast_rec is the documented
    # way to score one. But this function's default is "model" and
    # retrospective()'s is "mean", so a default-built retro is scored under a rule
    # that is not this function's default -- the same mismatch the error above
    # refuses, reached by saying nothing instead. Say which rule was used.
    .retro_rec <- retro[["forecast_rec"]] %||% unique(retro[["mase"]][["forecast_rec"]])
    if (length(.retro_rec) == 1L && !identical(.retro_rec, forecast_rec)) {
      message("Scoring under forecast_rec = \"", .retro_rec, "\", the rule this ",
              "`retro`'s peels were fitted with -- not this function's default ",
              "of \"", forecast_rec, "\".")
    }
  }
  # retrospective() returns rev(c(list(object), peels)), so the LAST element
  # is the input model and the peels run deepest-first. Score only the peels;
  # comparing the parent to itself would give MASE 0 and mean nothing.
  ml <- retro$Rceattle_list
  if (length(ml) < 2L) {
    stop("No peel survived, so there is nothing to score. Inspect ",
         "retrospective(..., peels = c(1, 1)) and read its $convergence.",
         call. = FALSE)
  }
  peel_models <- ml[-length(ml)]

  styr   <- object$data_list$styr
  endyr  <- object$data_list$endyr
  nspp   <- object$data_list$nspp
  spnames <- if (!is.null(object$data_list$spnames)) {
    object$data_list$spnames
  } else as.character(seq_len(nspp))

  by_year <- list()
  if ("model" %in% reference) {
    by_year <- .rce_mase_by_year_model(object, peel_models, quantity)
  }
  # --- reference = "observed": classic hindcast cross-validation ----------
  # A peel's index_data is filtered to endyr_peel, so its index_hat has no rows
  # for the held-out years. Rebuild each peel at estimateMode = 3 -- build, do
  # not optimize -- against the FULL index_data, with the peel's own parameters
  # and map. That yields the peel's PREDICTION at observations it never fitted,
  # which is the whole point, without letting it re-estimate anything.
  if ("observed" %in% reference) {
    idx_full <- object$data_list$index_data
    # Fleet-peel pairs with no observation at the peel's terminal year, so no
    # y_{t-h} and no scale for the MASE. Reported once, after the loop.
    .no_baseline <- NULL

    # Analytical catchability solves q INSIDE the template from log(obs/pred)
    # over every fitted index row, so handing the rebuild the held-out rows lets
    # them set the q that the prediction is then scored against -- measured, a
    # 4.2% shift in q on one fleet, moving the scored index_hat toward the very
    # observations being scored, and moving the retained years too. Refuse
    # rather than report a leaked skill score.
    # Which codes those are comes from the schema, not from a literal here: the
    # allowed set is `q_map`, and a second copy of it would drift the first time
    # a catchability form is added or retired.
    .q_map <- .rce_allowed_map("Catchability")
    .analytic_names <- c("Analytical", "AnalyticalArith")
    # Assert rather than intersect quietly: if either form is ever renamed, an
    # empty code set would take the integer arm dead, and the name arm below
    # reads the same literals, so BOTH would stop matching and the leak guard
    # would vanish without erroring.
    if (!all(.analytic_names %in% names(.q_map))) {
      stop("Internal: the Catchability map no longer names ",
           paste(setdiff(.analytic_names, names(.q_map)), collapse = ", "),
           "; hindcast_skill() cannot tell which fleets solve q analytically.",
           call. = FALSE)
    }
    .analytic_codes <- unname(.q_map[.analytic_names])
    .q <- object$data_list$fleet_control$Catchability
    .q_int <- suppressWarnings(as.integer(.q))
    .analytic <- (!is.na(.q_int) & .q_int %in% .analytic_codes) |
      (as.character(.q) %in% .analytic_names)
    if (any(.analytic, na.rm = TRUE)) {
      stop("hindcast_skill(reference = \"observed\") cannot be used on a model ",
           "with analytical catchability (fleet(s) ",
           paste(unique(object$data_list$fleet_control$Fleet_name[which(.analytic)]),
                 collapse = ", "),
           "): q is solved inside the model from the index observations, so a ",
           "peel asked to predict held-out rows would use those rows to set its ",
           "own q and the score would be leaked. Use reference = \"model\".",
           call. = FALSE)
    }

    for (m in peel_models) {
      ep <- m$data_list$endyr_peel
      if (is.null(ep) || ep >= endyr) next
      # Take the parent's SHAPE -- the peel's own data_list implies a different
      # parameter geometry and the rebuild fails the inits check -- but put the
      # PEEL's fixed inputs back. The peel carried weight-at-age, empirical
      # selectivity and ration forward from its last fitted year, and handing
      # the prediction the parent's copies would give it the true held-out
      # values for all three, which is the leak this is avoiding.
      dl <- object$data_list
      for (.nm in c("weight", "emp_sel", "ration_data")) {
        if (!is.null(m$data_list[[.nm]])) dl[[.nm]] <- m$data_list[[.nm]]
      }
      dl$index_data <- idx_full
      dl$endyr_peel <- ep
      pred <- try(suppressWarnings(suppressMessages(.refit_like(
        data_list = dl, inits = m$estimated_params, map = m$map,
        estimateMode = 3, phase = FALSE, getsd = FALSE))), silent = TRUE)
      # Say WHY a peel could not be scored. Dropping it silently leaves the
      # caller with "No peel had forecast years to score", which points at the
      # peel horizon when the real cause is a failed rebuild.
      if (inherits(pred, "try-error")) {
        warning("Could not predict the held-out index for the peel at ",
                "endyr_peel = ", ep, ": ",
                conditionMessage(attr(pred, "condition")), call. = FALSE)
        next
      }
      ih <- pred$quantities$index_hat
      if (is.null(ih)) {
        warning("The peel at endyr_peel = ", ep, " reported no index_hat.",
                call. = FALSE)
        next
      }
      ih <- as.numeric(ih)
      if (length(ih) != nrow(idx_full)) {
        warning("The peel at endyr_peel = ", ep, " predicted ", length(ih),
                " index rows but index_data has ", nrow(idx_full),
                "; cannot align them.", call. = FALSE)
        next
      }

      held <- which(idx_full$Year > ep & idx_full$Year <= endyr)
      for (flt in unique(idx_full$Fleet_code[held])) {
        rows <- held[idx_full$Fleet_code[held] == flt]
        if (!length(rows)) next
        # Naive = y_{t-h}, Kell et al. (2021) eq. 4: the observation h steps
        # before the one being scored. Every held-out row of this fleet has
        # t - h = ep by construction (years_ahead = Year - ep), so that is the
        # observation AT ep -- not merely the last one at or before it.
        #
        # The difference bites on an irregular survey. Taking the last prior
        # observation scores a one-step forecast against a two- or three-step
        # persistence baseline whenever the survey skipped year ep, inflating
        # the denominator and flattering the fleet.
        #
        # Without an observation at ep there is no y_{t-h}, so the BASELINE is
        # NA -- not the row. Dropping the row would throw away the forecast
        # error too, and $by_year is where a caller goes to recompute; an NA
        # denominator already means "undefined" everywhere else here. The
        # skipped fleet-peel pairs are collected and reported once below.
        at_ep <- which(idx_full$Fleet_code == flt & idx_full$Year == ep)
        if (length(at_ep)) {
          naive <- mean(idx_full$Observation[at_ep])
        } else {
          naive <- NA_real_
          .no_baseline <- rbind(.no_baseline, data.frame(
            fleet = flt, peel = endyr - ep, ep = ep, stringsAsFactors = FALSE))
        }
        sp <- idx_full$Species[rows][1]
        by_year[[length(by_year) + 1L]] <- data.frame(
          peel        = endyr - ep,
          species     = if (is.numeric(sp) && sp >= 1 && sp <= nspp) spnames[sp] else as.character(sp),
          quantity    = paste0("index_fleet_", flt),
          year        = idx_full$Year[rows],
          years_ahead = idx_full$Year[rows] - ep,
          forecast    = ih[rows],
          reference   = idx_full$Observation[rows],
          naive       = as.numeric(naive),
          stringsAsFactors = FALSE)
      }
    }

    # Say which fleet-peel pairs have no scale. Their forecast error is still in
    # $by_year; only the MASE is undefined, and a horizon whose every peel landed
    # here comes back NA rather than quietly averaging fewer peels.
    if (!is.null(.no_baseline)) {
      .fn <- object$data_list$fleet_control$Fleet_name
      warning(nrow(.no_baseline), " fleet-peel combination(s) have no survey ",
              "observation in the peel's terminal year, so there is no ",
              "y_{t-h} to scale against and their MASE is NA: ",
              paste0(.fn[.no_baseline$fleet] %||% .no_baseline$fleet,
                     " (peel ", .no_baseline$peel, ", needs ", .no_baseline$ep,
                     ")", collapse = "; "),
              ". An irregular survey will do this at most peel depths; the ",
              "forecast errors are still in $by_year.", call. = FALSE)
    }
  }

  if (!length(by_year)) {
    stop("No peel had forecast years to score. With reference = \"observed\" ",
         "check the survey years as well as the peel depths: a fleet is scored ",
         "only where it has observations after the peel.", call. = FALSE)
  }
  by_year <- do.call(rbind, by_year)

  # Grouped by HORIZON, averaging across peels -- Kell et al. (2021) eq. 5, whose
  # sums run over t = T-n .. T at a FIXED h. Each peel contributes one
  # |observed - predicted| at a given steps-ahead, and the mean is over those.
  #
  # This was grouped by peel and averaged over horizons instead, which is the
  # transpose: it reported one MASE per peel rather than per horizon, and made
  # the h = 1 denominator a SINGLE |naive - reference| instead of a mean over
  # peels. That is where the "a one-year MASE is dominated by its own
  # denominator" caveat came from -- a property of the mis-grouping, not of the
  # statistic.
  mase <- .rce_mase_aggregate(by_year)

  # Both tables leave in $mohns's column vocabulary; see .rce_as_mohns_names().
  list(mase = .rce_as_mohns_names(mase),
       by_year = .rce_as_mohns_names(by_year))
}


#' Score a peeled forecast against the full model, one row per peel-year.
#'
#' The `reference = "model"` half of [hindcast_skill()], shared with
#' [retrospective()] so the statistic has ONE implementation. Each peel's
#' estimate over the years it did not see is compared to the unpeeled model's,
#' with the peel's own terminal estimate held flat as the baseline.
#'
#' @param object the unpeeled fit.
#' @param peel_models the peels, excluding `object` itself.
#' @param quantity which reported quantities to score.
#' @return a list of data frames, one per peel x quantity x species, with
#'   `peel`, `species`, `quantity`, `year`, `years_ahead`, `forecast`,
#'   `reference` and `naive`. These are the INTERNAL names the Kell eq. 5
#'   arithmetic is written against; `.rce_as_mohns_names()` renames them on the
#'   way out of [hindcast_skill()]. Empty where no peel has forecast years.
#' @noRd
.rce_mase_by_year_model <- function(object, peel_models, quantity) {
  styr  <- object$data_list$styr
  endyr <- object$data_list$endyr
  nspp  <- object$data_list$nspp
  spnames <- object$data_list$spnames %||% as.character(seq_len(nspp))

  out <- list()
  for (m in peel_models) {
    ep <- m$data_list$endyr_peel
    # A peel with no forecast years (endyr_peel == endyr) scores nothing.
    if (is.null(ep) || ep >= endyr) next
    fyrs <- (ep + 1):endyr
    cols <- fyrs - styr + 1L
    last <- ep - styr + 1L

    for (q in quantity) {
      ref <- object$quantities[[q]]
      fc  <- m$quantities[[q]]
      if (is.null(fc) || is.null(ref) || ncol(fc) < max(cols)) next
      for (sp in seq_len(nspp)) {
        # A species with input numbers-at-age was not forecast: estDynamics > 0
        # makes the peel reproduce the parent bit for bit, so the forecast error
        # is 0 and the MASE would read as perfect skill for a series the model
        # copied from its input. Scored NA instead -- the convention here for a
        # quantity that cannot be judged.
        .fixed <- isTRUE((object$data_list$estDynamics %||% 0)[sp] > 0)
        out[[length(out) + 1L]] <- data.frame(
          peel        = endyr - ep,
          species     = spnames[sp],
          quantity    = q,
          year        = fyrs,
          years_ahead = seq_along(fyrs),
          forecast    = if (.fixed) NA_real_ else as.numeric(fc[sp, cols]),
          reference   = as.numeric(ref[sp, cols]),
          # "Assume nothing changes": the peel's own terminal estimate held flat.
          naive       = if (.fixed) NA_real_ else as.numeric(fc[sp, last]),
          stringsAsFactors = FALSE)
      }
    }
  }
  out
}


#' Mean absolute scaled error, per horizon, averaged across peels.
#'
#' Kell et al. (2021) eq. 5. The sums run over \eqn{t = T-n \ldots T} at a
#' FIXED \eqn{h}, so each peel contributes one term at a given steps-ahead and
#' the mean is over peels -- not over horizons within a peel, which is the
#' transpose and was what this computed before 5.23.0.9002.
#'
#' Shared by [hindcast_skill()] and [retrospective()]; eq. 5 is implemented
#' once, here.
#'
#' @param by_year rows from [.rce_mase_by_year_model()] or the observed path.
#' @return one row per `years_ahead` x `species` x `quantity`, with `n_peels`,
#'   `peels_used` (which peels the row rests on), `mae_forecast`, `mae_naive`
#'   and `mase`. Internal names; renamed on the way out, as above.
#' @noRd
.rce_mase_aggregate <- function(by_year) {
  key <- interaction(by_year$years_ahead, by_year$species, by_year$quantity,
                     drop = TRUE, lex.order = TRUE)
  mase <- do.call(rbind, lapply(split(by_year, key), function(z) {
    # ONE term per peel, then the mean across peels. Eq. 5 sums over t, each a
    # different peel's terminal year, so a peel contributes once. On the model
    # path that is already true. On `reference = "observed"` it is not: a fleet
    # may legally carry two index rows in the same Year at different Months
    # (data_check() tests duplicates on Fleet_code/Year/Month), and both land in
    # the same years_ahead from the same peel. Averaging within the peel first
    # stops it being weighted twice, and makes n_peels a count of PEELS.
    af <- tapply(abs(z$forecast - z$reference), z$peel, mean)
    an <- tapply(abs(z$naive    - z$reference), z$peel, mean)
    # A peel that could not be scored -- an NA anywhere in its forecast years,
    # or a species whose dynamics are fixed -- drops OUT rather than nulling the
    # horizon for every other peel. n_peels and peels_used then count what
    # actually contributed, which is what those columns claim to mean. With no
    # peel left the row is NA, not a mean of nothing.
    .ok <- !is.na(af) & !is.na(an)
    af <- af[.ok]
    an <- an[.ok]
    mae_f <- if (length(af)) mean(af) else NA_real_
    mae_n <- if (length(an)) mean(an) else NA_real_
    data.frame(
      years_ahead  = z$years_ahead[1], species = z$species[1],
      quantity     = z$quantity[1],
      n_peels      = length(af),
      # WHICH peels contributed, not just how many. For a contiguous `peels`
      # with nothing dropped this is just h:max(peels) and adds nothing -- it
      # earns its place when `peels` is a subset like c(2,5,9), when a peel was
      # dropped as non-converged, or when one was skipped as unscoreable above,
      # since then n_peels alone cannot say WHICH fits the row rests on.
      #
      # Named peels_used, not peel_depths: retrospective()'s own $peel_depths is
      # the depths REQUESTED, including ones that were dropped, and two columns
      # a level apart with the same name and different meanings read as a bug.
      peels_used   = paste(sort(as.integer(names(af))), collapse = ","),
      mae_forecast = mae_f,
      mae_naive    = mae_n,
      # NA rather than Inf when persistence was exactly right: an undefined
      # ratio is not infinitely bad skill, and Inf would poison any mean. NA
      # when the baseline itself is undefined -- a fleet with no observation at
      # the peel's terminal year has no y_{t-h}, so mae_n is NA and `NA > 0`
      # would be an error rather than a verdict.
      mase       = if (isTRUE(mae_n > 0)) mae_f / mae_n else NA_real_,
      stringsAsFactors = FALSE)
  }))
  rownames(mase) <- NULL
  mase[order(mase$quantity, mase$species, mase$years_ahead), ]
}


#' Put a skill table into `$mohns`'s column vocabulary
#'
#' @description
#' The public tables share the leading `Object`, `Forecast year`, `N`,
#' `species` columns with [retrospective()]'s `$mohns`, so the bias and skill
#' tables merge on those three and read in the same terms. `Forecast year` is a
#' HORIZON in both -- steps ahead of the peel's terminal year, not a calendar
#' year -- and `$mohns` additionally carries horizon 0, the terminal retrospective
#' bias, which has no forecast to score.
#'
#' Applied at the exit only. The internals keep the names the Kell eq. 5
#' arithmetic is written against, so renaming cannot move a number.
#'
#' @param d a `$mase` or `$by_year` frame.
#' @return `d` with the shared columns renamed and moved to the front.
#' @noRd
.rce_as_mohns_names <- function(d) {
  if (is.null(d) || !nrow(d)) return(d)
  # `mase` is deliberately NOT renamed. The other three had to move to key the
  # join against $mohns; `mase` collides with nothing, and `$mase$mase` is the
  # most natural read there is -- a rename would make it return NULL silently.
  # It also pairs with $mohns's lowercase `rho`.
  ren <- c(quantity = "Object", years_ahead = "Forecast year", n_peels = "N")
  hit <- intersect(names(ren), names(d))
  if (length(hit)) names(d)[match(hit, names(d))] <- unname(ren[hit])
  # The shared keys first, then the value, then its parts -- so a reader meets
  # `mase` before the two means it is a ratio of.
  lead <- intersect(c("Object", "Forecast year", "N", "species", "mase"),
                    names(d))
  d <- d[, c(lead, setdiff(names(d), lead)), drop = FALSE]
  # split() keyed the rows, so they come back in key order with its rownames.
  # Both tables leave through here, so they come back in ONE order -- quantity,
  # then horizon, then species. Also drops the rownames, which are split keys
  # rather than anything a caller can use.
  # method = "radix" so the order does not depend on the collation locale: plain
  # order() on a character column sorts "R" before or after "biomass" depending
  # on LC_COLLATE, which would make row order differ between a dev machine and a
  # C-collate R CMD check.
  ord <- intersect(c("Object", "Forecast year", "species"), names(d))
  if (length(ord)) {
    d <- d[do.call(order, c(unname(as.list(d[, ord, drop = FALSE])),
                            list(method = "radix"))), , drop = FALSE]
  }
  rownames(d) <- NULL
  d
}
