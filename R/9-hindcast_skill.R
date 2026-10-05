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
#' \eqn{t = T-n \ldots T} at fixed \eqn{h}. `n_peels` is how many peels
#' contributed. So `peels` buys precision in each MASE rather than more rows, and
#' a short horizon is not inherently noisier than a long one -- it is averaged
#' over the same peels.
#'
#' Read `n_peels` before reading a MASE. A horizon `h` can only be scored by a
#' peel at least `h` years deep, so `n_peels` FALLS as the horizon grows: with
#' `peels = 2:4` the one-year-ahead MASE averages three peels and the
#' four-year-ahead one averages a single peel. The deepest horizon in any run is
#' therefore the least precise, and a MASE resting on one peel says as much about
#' that peel as about the projection. Widen `peels` to score a long horizon, and
#' treat the last row or two of `$mase` with the same caution you would give any
#' n = 1 statistic.
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
#'   exactly. More peels do not add rows to `$mase` -- they add terms to the mean
#'   inside each one, since the scoring averages across peels at a fixed horizon.
#'   A vector is still useful for cost: `peels = 3:10` pays for two fewer model
#'   fits than `1:10`. Ignored when `retro` is supplied.
#' @param quantity Quantities to score against the full model. Any of `"ssb"`,
#'   `"biomass"`, `"R"`. Ignored when `reference = "observed"`.
#' @param reference What to score the projection against.
#'   `"model"` (default) compares each peel's projection to the full time-series
#'   model's estimate, taking that as the best available account of what
#'   happened. `"observed"` is classic hindcast cross-validation: it compares the
#'   peel's PREDICTED SURVEY INDEX to the index values actually observed in the
#'   held-out years. Both can be asked for. `"model"` scores the quantity a
#'   projection is used for (SSB, recruitment); `"observed"` scores the only
#'   thing that was really measured, and is the version comparable to MASE
#'   values published for other assessment platforms.
#' @param forecast_rec How the peeled years get their recruitment.
#'   `"model"` (default here) uses the model's own projection rule, in
#'   precedence order: `proj_mean_rec = TRUE` projects at mean recruitment,
#'   whatever process the model carries, because that is what the setting means;
#'   otherwise the LATENT STATES supply it where the deviations are random
#'   effects (`random_rec = TRUE`, or a DSEM), so an AR1's autocorrelation or a
#'   DSEM's lagged and covariate paths propagate into the forecast; otherwise
#'   recruitment comes off the stock-recruit curve. `"mean"` forces the
#'   historical mean for all of them, which is [retrospective()]'s default and
#'   the convention Mohn's rho is computed under.
#'
#'   The default differs from [retrospective()]'s on purpose. A peel's forecast
#'   years are hindcast years, so `proj_mean_rec` cannot reach them by itself,
#'   and with `"mean"` every model forecasts identically -- which makes this
#'   function unable to tell projection methods apart, the one thing it exists
#'   to do. Measured on the GOA arrowtooth model, `proj_mean_rec = TRUE` and
#'   `FALSE` returned byte-identical MASE under `"mean"`.
#' @param retro Optionally an already-computed [retrospective()] result, to
#'   avoid refitting when both are wanted. Note a `retro` computed with
#'   [retrospective()]'s defaults carries `forecast_rec = "mean"`.
#' @param ... Passed to [retrospective()] (`cores`, `getsd`, `rescale`).
#'
#' @return A list with
#'   \describe{
#'     \item{`mase`}{one row per years-ahead x species x quantity, with `mase`,
#'       `mae_forecast`, `mae_naive` and `n_peels`, the number of peels averaged
#'       over at that horizon.}
#'     \item{`by_year`}{the underlying series: peel, species, quantity, year,
#'       years-ahead, `forecast`, `reference`, `naive`.}
#'   }
#'
#' @references
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
  if ("model" %in% reference) for (m in peel_models) {
    ep <- m$data_list$endyr_peel
    # A peel with no forecast years (endyr_peel == endyr) scores nothing.
    if (is.null(ep) || ep >= endyr) next
    fyrs <- (ep + 1):endyr
    cols <- fyrs - styr + 1L
    last <- ep - styr + 1L

    for (q in quantity) {
      ref <- object$quantities[[q]]
      fc  <- m$quantities[[q]]
      if (is.null(fc) || ncol(fc) < max(cols)) next
      for (sp in seq_len(nspp)) {
        by_year[[length(by_year) + 1L]] <- data.frame(
          peel        = endyr - ep,
          species     = spnames[sp],
          quantity    = q,
          year        = fyrs,
          years_ahead = seq_along(fyrs),
          forecast    = as.numeric(fc[sp, cols]),
          reference   = as.numeric(ref[sp, cols]),
          # "Assume nothing changes": the peel's own terminal estimate held flat.
          naive       = as.numeric(fc[sp, last]),
          stringsAsFactors = FALSE)
      }
    }
  }
  # --- reference = "observed": classic hindcast cross-validation ----------
  # A peel's index_data is filtered to endyr_peel, so its index_hat has no rows
  # for the held-out years. Rebuild each peel at estimateMode = 3 -- build, do
  # not optimize -- against the FULL index_data, with the peel's own parameters
  # and map. That yields the peel's PREDICTION at observations it never fitted,
  # which is the whole point, without letting it re-estimate anything.
  if ("observed" %in% reference) {
    idx_full <- object$data_list$index_data

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
        # Naive = the last index value this fleet was observed at before the
        # peel, carried forward. Undefined if the fleet has no pre-peel
        # observation, in which case there is no persistence baseline.
        prior <- which(idx_full$Fleet_code == flt & idx_full$Year <= ep)
        if (!length(prior)) next
        naive <- idx_full$Observation[prior[which.max(idx_full$Year[prior])]]
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
  }

  if (!length(by_year)) {
    stop("No peel had forecast years to score.", call. = FALSE)
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
  key <- interaction(by_year$years_ahead, by_year$species, by_year$quantity,
                     drop = TRUE, lex.order = TRUE)
  mase <- do.call(rbind, lapply(split(by_year, key), function(z) {
    mae_f <- mean(abs(z$forecast - z$reference))
    mae_n <- mean(abs(z$naive    - z$reference))
    data.frame(
      years_ahead  = z$years_ahead[1], species = z$species[1],
      quantity     = z$quantity[1],
      n_peels      = nrow(z),
      mae_forecast = mae_f,
      mae_naive    = mae_n,
      # NA rather than Inf when persistence was exactly right: an undefined
      # ratio is not infinitely bad skill, and Inf would poison any mean.
      mase       = if (mae_n > 0) mae_f / mae_n else NA_real_,
      stringsAsFactors = FALSE)
  }))
  rownames(mase) <- NULL
  mase <- mase[order(mase$quantity, mase$species, mase$years_ahead), ]

  list(mase = mase, by_year = by_year)
}
