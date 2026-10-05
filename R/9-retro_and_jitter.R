#' Retrospective peels
#'
#' @description Calculate Mohn's rho and run retrospective peels for an Rceattle model. The function also evaluates retrospective forecast skill. To evaluate both retrospective bias and forecast skill, the function uses the map functionality of TMB to peel the model:
#' 1. Filters data, filters fixed inputs, and maps out time-varying parameters for the peeled years. All time-varying parameters for the peeled years are set to the terminal year of the model for that peel.
#' 2. Fits the peeled model.
#' 3. Turns off all hindcast parameters, turns on F for the peeled years, and fits to the peeled catch series to update the "forecast" dynamics given projection assumptions and observed catch from the peeled years.
#'
#' @inheritParams rceattle-refit-args
#' @param peels which retrospective peels to fit. A single number `n` fits every
#'   peel from 1 to `n`, as it always has. A VECTOR instead fits exactly the peel
#'   depths named, so `peels = 2:10` skips the one-year peel and `peels = c(1, 5,
#'   10)` fits three. Each element is a number of years removed from `endyr`, not
#'   a year. Depths are sorted and de-duplicated, so the returned list runs
#'   deepest-first whatever order they are given in. Mohn's rho is averaged over
#'   whichever peels are fitted, so a subset changes the number it reports --
#'   a rho from `peels = 2:10` is not comparable to one from `peels = 10`.
#'
#'   `peels = 5` and `peels = c(5)` are the same object in R, so a length-1
#'   value always means 1:5 and one depth on its own cannot be written that way.
#'   Name it twice -- `peels = c(5, 5)` de-duplicates to the 5-year peel alone --
#'   when you want to inspect a single peel without paying for the shallower ones.
#' @param rescale TRUE/FALSE whether to subset and rescale environmental predictors for the range of peel years.
#' @param nyrs_forecast Number of forecast years to calculate Mohn's Rho in addition to terminal year
#' @param forecast_rec How the peeled years get their recruitment. `"mean"`
#'   (default) projects them at the bias-adjusted historical mean, which is the
#'   convention Mohn's rho has always been computed under. `"model"` uses the
#'   model's own projection rule instead, in precedence order: `proj_mean_rec =
#'   TRUE` projects at mean recruitment whatever process the model carries;
#'   otherwise the latent states supply it where the deviations are random
#'   effects (`random_rec`, or a DSEM), so an AR1's autocorrelation or a DSEM's
#'   lagged and covariate paths propagate into the forecast; otherwise
#'   recruitment comes off the stock-recruit curve. Use `"model"` to compare
#'   projection methods -- see [hindcast_skill()], which defaults to it.
#' @param getsd whether each peel runs \code{TMB::sdreport} (standard errors).
#'   Costs an extra model build per peel; see Details.
#'   Mohn's rho uses only point estimates, so \code{FALSE} is faster with no
#'   effect on rho. Default \code{NULL} inherits the input model's setting
#'   (\code{TRUE} if it was fit with \code{getsd = TRUE}, i.e. carries an
#'   \code{sdrep}); the returned peel models then carry standard errors only
#'   when \code{getsd} is \code{TRUE}.
#' @param phase whether each peel is refitted in phases (default \code{TRUE}).
#'   A peel restarts from the unpeeled fit's starting values with a year removed;
#'   without phasing the parameters barely move, the peels sit on top of the full
#'   model, and Mohn's rho is biased towards zero. Change it only deliberately.
#'
#' @details
#' Each peel is fitted twice: a peeled hindcast, then a forecast refit that
#' estimates only the peeled years' F. The second holds every hindcast parameter
#' fixed, so on its own it reports a standard error of zero for the whole
#' hindcast. Under `getsd = TRUE` the peel is therefore rebuilt at those same
#' parameters with the hindcast free in the map and reported from there. Nothing
#' is re-estimated, so no point estimate moves.
#'
#' @return a list of 1. list of Rceattle models, 2. retrospective bias
#'   (`$mohns`, Mohn's rho) and 3. forecast skill (`$mase`).
#'
#'   Both are LONG and share their key columns -- `Object`, `Forecast year`,
#'   `species` -- so they merge on those three:
#'   `merge(retro$mohns, retro$mase, by = c("Object", "Forecast year", "species"))`.
#'   `Object` is the quantity (`"biomass"`, `"ssb"`, `"R"`, `"F_spp"`),
#'   `Forecast year` a HORIZON in steps ahead of the peel's terminal year rather
#'   than a calendar year, and `species` a name from `spnames`.
#'
#'   `N` means the same thing in both -- the number of observations the statistic
#'   was computed from -- but **join on the three keys, not on `N` as well.** The
#'   two are computed from different peels (rho scores a horizon the peel
#'   reached; MASE additionally needs a usable naive baseline), so they agree for
#'   a well-behaved species and diverge exactly where a term was unusable.
#'   Joining on `N` would silently DISCARD those rows; merging on the keys brings
#'   both columns through, and a disagreement between them is a finding.
#'
#'   They answer different questions and cover different horizons. Rho (`rho`)
#'   is a SIGNED relative error, so errors cancel across peels and it measures
#'   BIAS; its `Forecast year` runs `0:nyrs_forecast`, with 0 the terminal-year
#'   bias the +/-0.2 rule is about. MASE (`mase`) is an absolute error scaled by
#'   "assume nothing changes", so it measures SKILL and `mase < 1` beats that
#'   baseline. Read [hindcast_skill()]'s `@details` before leaning on the word
#'   "persistence": under `reference = "model"`, which is what this returns, the
#'   baseline is the PEEL's terminal estimate against the FULL model's, so it
#'   carries the peel's own retrospective bias and is not the persistence error
#'   of one series; it has no horizon 0, because a zero-step forecast is the
#'   terminal estimate compared with itself.
#'
#'   So the merge is an INNER join on purpose, and drops more than the horizons:
#'   rho covers `0:nyrs_forecast` and MASE `1:max(peels)`, so they share
#'   `1:min(nyrs_forecast, max(peels))`; and rho scores `F_spp` while MASE is
#'   scored only on `ssb`, `biomass` and `R`, so every `F_spp` row is rho-only.
#'   Use `all.x = TRUE` to keep them.
#'
#'   `$mase` additionally carries `mae_forecast` and `mae_naive` (the ratio's two
#'   halves), `peels_used` (which depths the row rests on, since `N` alone cannot
#'   say which), and `forecast_rec`, because it scores whatever that argument
#'   asked for. Under the default `"mean"` every model projects the peeled years
#'   at its own historical mean whatever process it carries, so two models'
#'   tables differ by their hindcast FITS, not their projection rules. Use
#'   `forecast_rec = "model"` to score the model's own rule, and see
#'   [hindcast_skill()] to compare the two on one fit. `N` falls as the horizon
#'   grows, because a horizon `h` is scored only by peels at least `h` deep, so
#'   the deepest rows rest on one or two fits. A species with `estDynamics > 0`
#'   is `NA` -- its numbers-at-age are input, so it was never forecast.
#'
#'   A peel that did not converge is dropped, so \code{Rceattle_list} can be
#'   shorter than one plus the number of peels asked for (a message reports how
#'   many); \code{peels_requested} is that number and \code{peel_depths} the
#'   depths themselves. Each entry is named for its own terminal year
#'   (\code{Year_2017}, ...) rather than by position, so index it by name --
#'   \code{Rceattle_list[[3]]} is not necessarily the 3-year peel, and with a
#'   subset such as \code{peels = 2:10} it never is. With no peel left, Mohn's
#'   rho is \code{NaN} and the function warns.
#'
#'   Each peel reports its own terminal year as \code{data_list$endyr}, so plots
#'   draw it only as far as it was fit and the peels fan out.
#'
#'   A peel still estimates the years it dropped: they are its retrospective
#'   forecast, fit to the observed catch with recruitment held at the peel's mean
#'   and the survey and composition data withheld. Three years therefore matter,
#'   and each peel carries all three:
#'   \describe{
#'     \item{\code{endyr}, \code{endyr_peel}}{the peel's terminal year -- what it
#'       was fit through. Equal to each other.}
#'     \item{\code{endyr_full}}{the unpeeled model's terminal year, where the
#'       retrospective forecast ends.}
#'     \item{\code{projyr}}{the end of the harvest-control-rule projection.}
#'   }
#'   So the forecast years are those after \code{endyr_peel} through
#'   \code{endyr_full}, and the projection follows through \code{projyr};
#'   \code{incl_proj = TRUE} plots both. Take the forecast years as
#'   \code{endyr_peel + seq_len(endyr_full - endyr_peel)}, which is empty for the
#'   unpeeled model, rather than \code{(endyr_peel + 1):endyr_full}, which counts
#'   \emph{down} there.
#'
#'   Mohn's rho is computed from \code{endyr_peel} and is unaffected by any of
#'   this.
#'
#'   Catchability is estimated only for a fleet that carries fitted index rows
#'   (see \code{\link{build_map}}), and a peel moves \code{endyr}. A survey whose
#'   index observations all fall in the peeled-off years therefore has no q
#'   estimated in that peel -- the parameter count is not constant across peels.
#'   That is deliberate: a q with no index to inform it is a flat direction in
#'   the likelihood. It does not affect Mohn's rho, which is computed from SSB,
#'   but it does mean \code{npar} and the reported catchability differ between a
#'   shallow and a deep peel for such a fleet.
#'
#' @examples
#' \donttest{
#' data(BS2017SS)
#' ss_run <- fit_mod(data_list = BS2017SS,
#'     inits = NULL, file = NULL,
#'     estimateMode = 0, random_rec = FALSE,
#'     msmMode = 0, avgnMode = 0,
#'     phase = FALSE, verbose = 0)
#' retro <- retrospective(ss_run, peels = 10)
#' }
#' @export
retrospective <- function(object = NULL, peels = 5, rescale = FALSE, nyrs_forecast = 3, cores = NULL, getsd = NULL, phase = TRUE, fit_control = NULL, forecast_rec = c("mean", "model"), Rceattle = NULL) {
  # `phase`, `fit_control` and `forecast_rec` are appended after the arguments
  # that predate them, so a positional call keeps its meaning. The deprecated
  # `Rceattle` formal sits last; see R/0-deprecate.R.
  if (!missing(Rceattle))
    object <- .rce_deprecated_arg(Rceattle, !missing(object), "Rceattle", "object", "retrospective")
  # Cleared once consumed: .parallel_lapply() exports every binding in this frame
  # to the PSOCK workers, and two names for one fitted model send it twice.
  Rceattle <- NULL

  if (!inherits(object, "Rceattle")) {
    stop("Object is not of class 'Rceattle'")
  }
  forecast_rec <- match.arg(forecast_rec)

  # A peel does not shorten the model: it sets endyr_peel and turns off DATA
  # after it, so a DSEM still spans every year. Its latent states stay free and
  # in `random`, and the Laplace approximation integrates the peeled-year states
  # out against the GMRF prior with no data informing them -- which is the
  # peeled marginal likelihood. Do NOT mirror what rec_dev does below (zero the
  # tail, map it out): pinning is inert for an independent deviate but not for
  # states coupled through the RAM, where the pinned zeros stay in the quadratic
  # form and shrink the terminal retained state.
  #
  # rescale = TRUE works under a DSEM and is the more defensible setting for
  # one. It standardizes env_data on styr:endyr_peel only, so the peel does not
  # centre its covariates using years it is supposed to have not seen -- with
  # the full-series mean and sd, post-peel information leaks in and the
  # retrospective understates error. It does not change the latent dimensions:
  # build_dsem_objects() full-joins env_data back onto styr:dsem_endyr, so the
  # trimmed years return as NA and the DSEM predicts those covariate values from
  # its own process instead of reading data the peel should not have.

  # Peels inherit the input model's sdreport setting unless overridden. Mohn's
  # rho reads only point estimates, so getsd = FALSE is faster and rho-neutral.
  if (is.null(getsd)) getsd <- !is.null(object$sdrep)

  ctl <- .rce_refit_control(fit_control, "retrospective")
  if (!is.null(ctl$phase)) phase <- ctl$phase
  if (!is.null(ctl$getsd)) getsd <- ctl$getsd

  # Get objects
  object$data_list$endyr_peel <- object$data_list$endyr
  # Terminal year of the model being peeled. Each peel reports its OWN terminal
  # year as `endyr` (see run_one_peel), which makes the plots fan out but leaves
  # `endyr` and `endyr_peel` holding the same value -- so without this the
  # unpeeled terminal year is unrecoverable from a peel, and with it the boundary
  # between the retrospective FORECAST years, (endyr_peel + 1):endyr_full, and
  # the true projection, (endyr_full + 1):projyr. Set once here: run_one_peel
  # copies this data_list, and extra fields survive the refits the same way
  # `endyr_peel` already does.
  object$data_list$endyr_full <- object$data_list$endyr
  data_list <- object$data_list # used by Mohn's rho block below
  endyr <- object$data_list$endyr
  styr <- object$data_list$styr
  nyrs <- length(styr:endyr)
  projyr <- object$data_list$projyr
  nyrs_proj <- projyr - styr + 1

  # Which peel depths to fit. A scalar `n` keeps the historical meaning, 1:n;
  # a vector names the depths exactly, so a caller can skip the shallow peels,
  # which cost a model fit each and contribute only to the shortest horizons in
  # hindcast_skill(). Validate here rather than let a
  # bad value reach run_one_peel(): `endyr - i` would silently build a model
  # ending before styr, or `(endyr_peel + 1):endyr` would count DOWN and fit the
  # retained years as a forecast.
  # `is.finite`, not just `anyNA`: Inf passes the whole-number check below
  # (`Inf != round(Inf)` is FALSE), then as.integer() turns it into NA and the
  # `< 1` test becomes `if (NA)`. So `peels = 1e10` -- a plausible typo -- used
  # to surface as a coercion warning plus "missing value where TRUE/FALSE
  # needed" rather than any of the messages here. A matrix is refused for the
  # same reason: as.integer() drops its `dim` and it would run as a vector.
  if (!is.numeric(peels) || length(peels) == 0L || !all(is.finite(peels)) ||
      !is.null(dim(peels))) {
    stop("`peels` must be a positive whole number, or a plain vector of peel ",
         "depths; got ",
         if (length(peels) == 0L) "an empty value"
         else if (!is.null(dim(peels)))
           paste0("a ", paste(dim(peels), collapse = "x"), " ", class(peels)[1])
         else paste0("a ", class(peels)[1], " of length ", length(peels),
                     # Both clauses guarded on is.numeric(): `is.finite()` has
                     # no list method and errors outright, and on a character or
                     # raw it returns FALSE for every element, which read as
                     # "containing Inf" on a plain `peels = "3"`.
                     if (is.numeric(peels) && anyNA(peels)) " containing NA"
                     else if (is.numeric(peels) && !all(is.finite(peels)))
                       " containing Inf"
                     else ""), ".",
         call. = FALSE)
  }
  if (any(peels != round(peels))) {
    stop("`peels` must be whole numbers of years: got ",
         paste(peels[peels != round(peels)], collapse = ", "), ".", call. = FALSE)
  }
  peels <- round(peels)
  # Both range checks run on the DOUBLES, before as.integer(). A finite value
  # above .Machine$integer.max -- `peels = 1e10`, the same typo -- coerces to NA
  # with only a coercion warning, and every test after that reads `if (NA)`.
  # The depth bound rejects it for what it is instead, since `max(peels)` is the
  # same number whether `peels` is a scalar `n` (whose deepest peel is n) or a
  # vector of depths.
  if (any(peels < 1)) {
    stop("`peels` must be at least 1: a depth of 0 is the unpeeled model, ",
         "which is already returned. Got ",
         paste(sort(unique(peels[peels < 1])), collapse = ", "), ".",
         call. = FALSE)
  }
  # Leave at least two years fitted. One year cannot inform a selectivity or a
  # recruitment deviation, so the peel would not be a shallower version of the
  # same model -- the comparison Mohn's rho rests on.
  if (max(peels) > endyr - styr - 1L) {
    stop("`peels` asks for a ", format(max(peels), scientific = FALSE),
         "-year peel, but the model spans ", styr, ":", endyr, " (", nyrs,
         " years); ",
         # On a 2-year model the bound is 0, and "the deepest peel is 0" would
         # contradict the >= 1 rule two checks up. Such a model cannot be peeled.
         if (endyr - styr - 1L < 1L)
           "a model this short cannot be peeled at all and leave two years fitted."
         else paste0("the deepest peel that leaves two years fitted is ",
                     endyr - styr - 1L, "."),
         call. = FALSE)
  }
  peels <- as.integer(peels)
  if (length(peels) == 1L) peels <- seq_len(peels)
  peel_seq <- sort(unique(peels))

  # Cross-platform parallel via parallel::parLapply on a PSOCK cluster
  # (same approach as run_mse). Respect the CRAN core limit
  # ('_R_CHECK_LIMIT_CORES_' is set during R CMD check;
  # parallel::makeCluster errors if we exceed 2 cores then).
  chk <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  cran_cap <- nzchar(chk) && !chk %in% c("false", "0", "no")
  if (is.null(cores)) {
    cores <- if (cran_cap) 2L else max(1L, parallel::detectCores() - 6L)
  } else {
    cores <- max(1L, as.integer(cores))
    if (cran_cap) cores <- min(cores, 2L)
  }
  use_parallel <- length(peel_seq) > 1L && cores > 1L

  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Once-per-call warnings ----
  # Both conditions read only the input model, so they are settled here rather
  # than inside run_one_peel(). Gating them on the peel index got both wrong
  # twice over: `peels = 2:10` never reaches peel 1 and lost them silently, and
  # a PSOCK worker's warning is discarded, so under the default `cores` they
  # never surfaced at all. The second of them is the one that says a DSEM
  # comparison is not actually comparing anything, which is the worst of these
  # to lose.
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # A fixed-effect linkage covariate may not be NA -- materialize_linkage()
  # rejects one, because model.matrix() would drop the row and misalign the
  # design -- so a peel keeps its post-peel values and Mohn's rho is
  # conditional on that covariate being known. Checked against the DEEPEST peel,
  # which is the one with the MOST post-peel rows: a deeper peel has a smaller
  # endyr_peel, so if any peel sees a post-peel covariate value the deepest one
  # does. Gating on the shallowest instead goes silent whenever env_data stops
  # short of endyr - 1 -- a real configuration, since rearrange_data() full-joins
  # a short series onto styr:projyr and mean-fills it -- while every deeper peel
  # still carries post-peel values.
  .env_cols_all <- setdiff(names(object$data_list$env_data), "Year")
  .keep_all <- intersect(.env_cols_all,
                         .rce_linkage_fixed_covariates(object$data_list))
  if (length(.keep_all) > 0 &&
      any(object$data_list$env_data$Year > endyr - max(peel_seq))) {
    warning("Covariate(s) ", paste(.keep_all, collapse = ", "), " enter the model ",
            "as fixed-effect linkage terms, which may not be NA, so each peel ",
            "still sees their post-peel values. Mohn's rho for a model with a ",
            "fixed-effect environmental covariate is conditional on that ",
            "covariate being known.", call. = FALSE)
  }

  # An unrecorded `random_vars` means the peels cannot tell a random effect from
  # a fixed one, so .rce_peel_map() pins every deviation. A pinned deviation is
  # still scored by its density, and "the deviation was exactly zero" is the
  # strongest possible evidence for a small process SD -- measured at -6.6% on
  # sigma at 5 peels and monotone in peel depth, so it becomes a trend in the
  # very quantity Mohn's rho reports.
  #
  # `is.null()`, not `length() == 0`: since 5.23.0 a fit that declared no random
  # effects records character(0), and pinning every deviation is then exactly
  # what the fit itself did -- nothing to report. Only an absent record is a
  # problem, and only a fit made before that is absent.
  if (is.null(object$random_vars)) {
    warning("This fit does not record which blocks were random effects, so ",
            "every peel pins every deviation. If any of them was a random ",
            "effect, that shrinks the estimated process SDs with peel depth ",
            "(-6.6% on sigma at 5 peels) and biases Mohn's rho. Refit with ",
            "fit_mod() at estimateMode 0 or 1 for an unbiased retrospective.",
            call. = FALSE)
  }

  # Say so when `forecast_rec = "model"` resolves to the mean anyway. A peel's
  # forecast years are HINDCAST years, so proj_mean_rec -- a projection switch --
  # reaching them is a surprise, and build_srr() defaults it to TRUE: every model
  # fitted without naming it gets Mohn's rho identical under both settings, and
  # hindcast_skill(), which defaults to "model" precisely to tell projection
  # methods apart, cannot. Note fit_mod() sets proj_mean_rec = 0 when a DSEM
  # carries estimate_projection = TRUE, so such a fit reads FALSE here and is
  # correctly left alone.
  if (identical(forecast_rec, "model") &&
      isTRUE(as.logical(object$data_list$proj_mean_rec))) {
    warning("forecast_rec = \"model\" is inert on this fit: proj_mean_rec = ",
            "TRUE, so the peeled years take mean recruitment and Mohn's rho ",
            "will match forecast_rec = \"mean\". Refit with ",
            "build_srr(proj_mean_rec = FALSE) for the model's own process to ",
            "supply the forecast.", call. = FALSE)
  }

  # Was this block a random effect in the HINDCAST? One definition, shared with
  # .rce_peel_map()'s local `pin()`, and read off `object$random_vars` for the
  # reason that function's @param gives: `obj$env$random` is empty on any fit
  # carrying an HCR, so it would call every block a fixed effect on exactly
  # those fits. An absent record pins everything, which is what the map does
  # too; `retrospective()` has already warned about it above.
  #
  # Unprefixed because it is a local in this frame, like `run_one_peel`, not one
  # of the package's unexported `.rce_*` helpers. The peel closure reaches it
  # either way: FORK inherits the frame, and PSOCK serialises the closure
  # together with its enclosing environment.
  pin_block <- function(nm) !(nm %in% object$random_vars)

  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Per-peel closure ----
  # Each peel only reads the original model, so peels are independent.
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  run_one_peel <- function(i) {
    # Warnings this peel raises that the CALLER has to see. A warning raised in
    # here is discarded by the FORK/PSOCK worker that runs the peel, so it is
    # collected on the returned model and re-raised after dispatch. Only for
    # conditions that are genuinely per-peel -- one that reads the input model
    # belongs above, with the other once-per-call warnings.
    .peel_warnings <- character(0)

    # * Get end year of peel ----
    data_list <- object$data_list
    endyr_peel <- endyr - i
    data_list$endyr_peel <- endyr_peel
    nyrs_peel <- endyr_peel - styr + 1
    peel_prj_yrs <- (endyr_peel+1):endyr
    nyrs_proj_peel <- length(peel_prj_yrs)


    # * Turn off data after endyr_peel ----
    data_list$index_data <- data_list$index_data |>
      dplyr::filter(Year <= endyr_peel)

    data_list$comp_data <- data_list$comp_data |>
      dplyr::filter(Year <= endyr_peel)

    data_list$caal_data <- data_list$caal_data |>
      dplyr::filter(Year <= endyr_peel)

    data_list$diet_data <- data_list$diet_data |>
      dplyr::filter(Year <= endyr_peel)

    peeled_catch_data <- data_list$catch_data |>
      dplyr::filter(Year > endyr_peel)
    data_list$catch_data$Catch[which(data_list$catch_data$Year > endyr_peel)] <- 0 # Set catch data in peeled years to 0 to avoid fitting


    # * Turn off fixed inputs after endyr_peel ----
    #FIXME ignores forecasted growth
    data_list$weight <- data_list$weight |>
      dplyr::filter(Year <= endyr_peel)

    data_list$emp_sel <- data_list$emp_sel |>
      dplyr::filter(Year <= endyr_peel)

    data_list$ration_data <- data_list$ration_data |>
      dplyr::filter(Year <= endyr_peel)

    # * Extend fixed inputs for "projection years"
    # - Assume weight/ration/empirical sel is same as last year of peel

    # - Weight
    #FIXME ignores forecasted growth
    proj_wt <- data_list$weight |>
      dplyr::filter(Year != 0)

    if(nrow(proj_wt) > 0){
      proj_wt <- proj_wt |>
        dplyr::group_by(Wt_index , Sex) |>
        dplyr::slice(rep(dplyr::n(),  nyrs_proj_peel)) |>
        dplyr::mutate(Year = peel_prj_yrs)
    }
    data_list$weight  <- rbind(data_list$weight, proj_wt) |>
      dplyr::arrange(Wt_index, Year)

    # - Empirical selectivity
    proj_emp_sel <- data_list$emp_sel |>
      dplyr::filter(Year != 0)

    if(nrow(proj_emp_sel) > 0){
      proj_emp_sel <- proj_emp_sel |>
        dplyr::group_by(Fleet_code, Sex) |>
        dplyr::slice(rep(n(),  nyrs_proj_peel)) |>
        dplyr::mutate(Year = peel_prj_yrs)
      data_list$emp_sel  <- rbind(data_list$emp_sel, proj_emp_sel) |>
        dplyr::arrange(Fleet_code, Year)
    }

    # - Ration data
    proj_ration_data <- data_list$ration_data |>
      dplyr::filter(Year != 0)

    if(nrow(proj_ration_data) > 0){
      proj_ration_data <- proj_ration_data |>
        dplyr::group_by(Species, Sex) |>
        dplyr::slice(rep(dplyr::n(),  nyrs_proj_peel)) |>
        dplyr::mutate(Year = peel_prj_yrs)
      data_list$ration_data  <- rbind(data_list$ration_data, proj_ration_data) |>
        dplyr::arrange(Species, Year)
    }

    # * Environmental predictors ----
    # Withhold post-peel covariate VALUES, always, not only when rescaling.
    # env_data is data the peel is supposed not to have: a QAR1 catchability
    # (est_index_q = 6) fits index_q_dev to env_index over every hindcast year,
    # and now that the peeled index_q_dev are free rather than pinned, leaving
    # env_data full-length lets a peel estimate its post-peel catchability from
    # post-peel environmental observations -- look-ahead in a diagnostic whose
    # whole purpose is to withhold the future.
    #
    # Blank the values; do NOT drop the rows. Linkage design matrices are built
    # POSITIONALLY from env_data (materialize_linkage() -> model.matrix()), so
    # dropping rows drops design columns and random-effect levels, and `inits`
    # -- carried from the parent fit -- no longer matches the model the peel
    # builds. On GOA pollock 2025, which carries ar1(1|Year) and rw(1|Year)
    # linkages, that refused every peel outright: beta_linkage 241 vs 237,
    # beta_linkage_re 110 vs 108.
    #
    # A FIXED-EFFECT linkage covariate is exempt: materialize_linkage() rejects
    # an NA in one by design, because model.matrix() would silently drop the row
    # and misalign the whole design. Those columns keep their values and the
    # peel is told which they are, once. Everything else -- a state-space
    # `observe` covariate, a DSEM covariate, an env_index column read by
    # Time_varying_q / Cindex / M1_indices -- is blanked, and the mean-fill in
    # rearrange_data() replaces it with the retained years' mean.
    .post_peel <- data_list$env_data$Year > endyr_peel
    .env_cols  <- setdiff(names(data_list$env_data), "Year")
    .keep      <- intersect(.env_cols, .rce_linkage_fixed_covariates(data_list))
    .blank     <- setdiff(.env_cols, .keep)
    if (any(.post_peel) && length(.blank) > 0) {
      data_list$env_data[.post_peel, .blank] <- NA_real_
    }
    # The caller was told about `.keep` once, before dispatch; see
    # "Once-per-call warnings" above.
    if(rescale && length(.env_cols) > 0){
      # Standardize on the RETAINED years only, so the peel does not centre its
      # covariates using years it is supposed not to have seen. The centre and
      # scale come from those years and are then applied to the whole column,
      # which keeps a kept fixed-effect covariate on the same scale as the years
      # the peel was fit to. A column with no retained variation is left alone
      # rather than turned into NaN.
      .kept <- as.matrix(data_list$env_data[!.post_peel, .env_cols, drop = FALSE])
      .ctr  <- colMeans(.kept, na.rm = TRUE)
      .scl  <- apply(.kept, 2, stats::sd, na.rm = TRUE)
      for (j in seq_along(.env_cols)) {
        if (!is.finite(.ctr[j]) || !is.finite(.scl[j]) || .scl[j] == 0) next
        data_list$env_data[[.env_cols[j]]] <-
          (data_list$env_data[[.env_cols[j]]] - .ctr[j]) / .scl[j]
      }
    }

    # * Adjust parameters ----
    #FIXME: adjust for forecasting via MVN
    inits <- object$estimated_params

    # Under a DSEM with family = "fixed" the covariate columns of dsem_x_tj ARE
    # the environmental data: the map holds them fixed, and their VALUES come
    # from inits. merge_dsem_params() keeps a caller-supplied dsem_x_tj rather
    # than the one just built, so warm-starting the peel from the parent carries
    # the parent's UNRESCALED covariate through and rescale = TRUE is silently
    # ignored -- measured: every path coefficient and the whole of x_tj came back
    # bit-identical with and without it. Drop the block so the rebuild supplies
    # covariates standardized on styr:endyr_peel. The recdev columns lose their
    # warm start, which the phased refit recovers and which is if anything the
    # more appropriate starting point for a peel.
    if (rescale && !is.null(inits$dsem_x_tj)) inits$dsem_x_tj <- NULL

    inits$rec_dev[, (nyrs_peel + 1):nyrs_proj] <- 0
    inits$log_M1_dev[,,,(nyrs_peel+1):nyrs_proj] <- inits$log_M1_dev[,,,nyrs_peel]
    inits$index_q_dev[,(nyrs_peel+1):nyrs] <- inits$index_q_dev[,nyrs_peel]
    inits$log_sel_slp_dev[,,,(nyrs_peel+1):nyrs] <- inits$log_sel_slp_dev[,,,nyrs_peel]
    inits$sel_inf_dev[,,,(nyrs_peel+1):nyrs] <- inits$sel_inf_dev[,,,nyrs_peel]
    inits$sel_coff_dev[,,,(nyrs_peel+1):nyrs] <- inits$sel_coff_dev[,,,nyrs_peel]

    # * Adjust map size ----
    # Which deviations the peel holds fixed. Extracted so the decision can be
    # exercised directly: it depends on convergence of nothing, but reaching it
    # through a fitted retrospective needs a model that converges, and the
    # selectivity forms it distinguishes are exactly the ones hardest to fit.
    map <- .rce_peel_map(object$map,
                         random_vars    = object$random_vars,
                         fleet_control  = object$data_list$fleet_control,
                         nyrs_peel      = nyrs_peel,
                         nyrs           = nyrs,
                         nyrs_proj      = nyrs_proj)

    # -- Map out Fdev for years with 0 catch to very low number
    zero_catch <- data_list$catch_data |>
      dplyr::filter(Year <= endyr &
                      Catch == 0) |>
      dplyr::mutate(Year = Year - styr + 1) |>
      dplyr::select(Fleet_code, Year) |>
      as.matrix()
    inits$log_F[zero_catch] <- -999
    map$mapList$log_F[zero_catch] <- NA
    map$mapFactor$log_F <- factor(map$mapList$log_F)

    # * Refit ----
    newmod <- suppressWarnings(
      suppressMessages(
        # Refit the peeled data_list, reusing its HCR / SR / M / growth
        # configuration; clamp the SR-switch, SR-fit-end, and suitability-end
        # years back to this peel's terminal year.
        .refit_like(
          data_list        = data_list,
          inits            = inits,
          map              = map,
          estimateMode     = ifelse(data_list$estimateMode < 3, 0, data_list$estimateMode),
          phase            = phase,  # see `phase` above: peels need phasing to move
          getsd            = getsd,
          srr_mse_switchyr = min(data_list$srr_mse_switchyr, endyr_peel),
          srr_hat_endyr    = min(data_list$srr_hat_endyr, endyr_peel),
          suit_endyr       = pmin(data_list$suit_endyr, endyr_peel))
      )
    )

    # Judge the peeled hindcast NOW. `newmod` is overwritten below by the
    # forecast-catch refit, whose map (build_map(debug = TRUE), plus the peeled
    # years' log_F) leaves only a handful of F parameters free -- so its gradient
    # says nothing about whether the hindcast converged. Mohn's rho reads the
    # hindcast quantities, which come from these parameters, so a peel that
    # diverged here has to be dropped even when the F refit lands cleanly.
    hindcast_converged <- .refit_converged(newmod)

    # Forecast ----
    peeled_pars <- newmod$estimated_params

    # - Add in peeled catch to fit to
    data_list$catch_data <- data_list$catch_data |>
      dplyr::filter(Year <= endyr_peel) |>
      rbind(peeled_catch_data) |>
      dplyr::arrange(Fleet_code, Year)

    # Keep the peeled hindcast's map before it is replaced: it is what frees the
    # hindcast again for the report-only pass after the forecast refit.
    map_hindcast <- map

    # - Update map
    # Only new parameter we are estimating in is the log_F of the peeled years
    map <- build_map(
      data_list = data_list,
      params = newmod$estimated_params,
      debug = TRUE,
      random_rec = newmod$data_list$random_rec)
    map$mapFactor$dummy <- as.factor(NA); map$mapList$dummy <- NA

    # build_map() knows nothing about the DSEM, and fit_mod() only fills the
    # DSEM entries in when it builds the map itself -- so without this the
    # dsem_* blocks would be ABSENT from a supplied map, which means unmapped,
    # which means estimated. Worse, `random` is derived from the map, so the
    # latent states would come back as fixed effects rather than being
    # integrated out. This refit only solves log_F against the peeled catch, so
    # hold every DSEM parameter at what the peeled hindcast estimated: that is
    # what debug = TRUE does to the rest of the hindcast, and Mohn's rho reads
    # those hindcast quantities. (Pinning is safe HERE, unlike in the peel
    # itself, because no peeled likelihood is being computed.)
    for (.nm in grep("^dsem_", names(newmod$estimated_params), value = TRUE)) {
      map$mapList[[.nm]]   <- rep(NA, length(newmod$estimated_params[[.nm]]))
      map$mapFactor[[.nm]] <- factor(map$mapList[[.nm]])
    }

    # - Turn on F for peeled years to fit to catch (matches full model)
    peeled_pars$log_F[,(nyrs_peel+1):nyrs] <- object$estimated_params$log_F[,(nyrs_peel+1):nyrs]
    map$mapList$log_F[,(nyrs_peel+1):nyrs] <- object$map$mapList$log_F[,(nyrs_peel+1):nyrs]
    map$mapFactor$log_F <-  factor(map$mapList$log_F)

    # The peeled years' recruitment, for the forecast-catch refit. Under a DSEM
    # the value goes into the LATENT STATE, not into rec_dev -- the model
    # derives rec_dev from the states, so writing rec_dev alone is a silent
    # no-op. Pinning is correct here, unlike in the peel itself: debug = TRUE
    # already holds the whole hindcast fixed and only log_F is solved, so no
    # peeled likelihood is computed and nothing is shrunk by fixing a state.
    #
    # forecast_rec chooses the rule, in precedence order:
    #   proj_mean_rec = TRUE                    mean recruitment, whatever
    #                                           process the model carries
    #   FALSE, deviations are random effects    the latent states supply it, so
    #                                           write nothing -- an AR1's
    #                                           autocorrelation or a DSEM's
    #                                           covariate paths propagate
    #   FALSE, deviations are fixed effects     off the stock-recruit curve,
    #                                           i.e. a zero deviation
    # "mean" is the default so Mohn's rho keeps the convention it has always had.
    .use_model_rec <- identical(forecast_rec, "model")
    .mean_rec      <- isTRUE(as.logical(newmod$data_list$proj_mean_rec))
    # Same source as the peel map's `pin()`, or the two disagree: with
    # random_rec = TRUE under an HCR the map would leave rec_dev free while this
    # said nothing supplies the forecast, and the "forecast" would be a
    # deterministic zero deviation inherited from inits.
    .states_supply <- .use_model_rec && !.mean_rec &&
      (.has_dsem(newmod) || !pin_block("rec_dev"))
    # The caller was told once, before dispatch, if `.mean_rec` makes
    # `forecast_rec = "model"` inert; see "Once-per-call warnings" above.
    if (!.states_supply) for(sp in 1:newmod$data_list$nspp){

      # -- where SR curve is estimated directly
      if(newmod$data_list$srr_fun == newmod$data_list$srr_pred_fun){
        rec_dev <- log(mean(newmod$quantities$R[sp,1:nyrs_peel]))  - log(newmod$quantities$R0[sp])
      }

      # -- OMs where SR curve is estimated as penalty (sensu Ianelli)
      if(newmod$data_list$srr_fun != newmod$data_list$srr_pred_fun){
        # Already a log-scale deviation, so take the mean directly.
        rec_dev <- mean((log(newmod$quantities$R) - log(newmod$quantities$R_hat))[sp, 1:nyrs_peel])

      }

      # Reached only with proj_mean_rec = FALSE and no recruitment process (the
      # random-effect cases returned above): the projection comes off the
      # stock-recruit curve, which is a zero deviation.
      if (.use_model_rec && !.mean_rec) {
        rec_dev <- 0
      }

      # - Update OM with devs
      peeled_pars$rec_dev[sp, (peel_prj_yrs - styr + 1)] <- replace(
        peeled_pars$rec_dev[sp, (peel_prj_yrs - styr + 1)],
        values =  rec_dev)

      # ... and, under a DSEM, into the state rec_dev is derived FROM, or the
      # line above is overwritten before it is ever read. rec_dev_col is 0-based
      # (it indexes the C++ column), hence the + 1.
      # Condition on the DSEM, not on dsem_x_tj: a non-DSEM fit carries a 0x0
      # dsem_x_tj (not NULL), so keying off the matrix enters this block on
      # every ordinary retrospective and is saved only by the row count being
      # zero. Assert the rows fit rather than filtering them: they always do
      # today (nyrs_dsem spans styr:endyr and the forecast years are inside it),
      # and the one future state where the filter WOULD bite -- a DSEM rebuilt
      # per peel over styr:endyr_peel -- is exactly where silently writing
      # nothing would be the bug this block exists to fix.
      if (!is.null(newmod$dsem) && !is.null(peeled_pars$dsem_x_tj)) {
        .col <- newmod$dsem$tmb_inputs$data$rec_dev_col[sp] + 1L
        .rows <- peel_prj_yrs - styr + 1L
        stopifnot(length(.col) == 1L, !is.na(.col),
                  max(.rows) <= nrow(peeled_pars$dsem_x_tj))
        # ADD BACK the lognormal bias correction. The template derives
        #   rec_dev = x_tj - bias_adjust_proc * margvar / 2
        # so writing the intended deviation straight into x_tj lands margvar/2
        # low -- measured at -22.1% in realised recruitment on BS2017SS with a
        # naive sem at the default bias_adjust_proc, and worse for a lagged sem,
        # where margvar is sigma^2/(1-rho^2) rather than sigma^2.
        .mv <- newmod$quantities$dsem_margvar_tj
        .adj <- if (!is.null(.mv) && nrow(.mv) >= max(.rows) && ncol(.mv) >= .col) {
          as.numeric(.mv[.rows, .col])
        } else rep(0, length(.rows))
        .bias <- as.numeric(newmod$data_list$bias_adjust_proc %||% 1)
        peeled_pars$dsem_x_tj[.rows, .col] <- rec_dev + .bias * .adj / 2
      }
    }

    newmod <- suppressMessages(
      suppressWarnings(
        # Second refit: same peeled configuration, now started from the
        # bias-adjusted peeled parameters with peeled-year F turned back on.
        .refit_like(
          data_list        = data_list,
          inits            = peeled_pars,
          map              = map,
          estimateMode     = ifelse(data_list$estimateMode < 3, 0, data_list$estimateMode),
          phase            = phase,  # see `phase` above: peels need phasing to move
          getsd            = getsd,
          srr_mse_switchyr = min(data_list$srr_mse_switchyr, endyr_peel),
          srr_hat_endyr    = min(data_list$srr_hat_endyr, endyr_peel),
          suit_endyr       = pmin(data_list$suit_endyr, endyr_peel))
      )
    )

    # * Standard errors for the peeled hindcast ----
    # The refit above holds the hindcast fixed, so its sdreport gives every
    # hindcast quantity a standard error of zero. Report instead from the same
    # parameters with the hindcast free: estimateMode 3 builds without
    # optimizing, so no point estimate moves.
    if (isTRUE(getsd) && !is.null(newmod)) {
      sd_map <- map_hindcast
      # The peeled years' F was estimated too, so it belongs in the report.
      sd_map$mapList$log_F[, (nyrs_peel + 1):nyrs] <-
        object$map$mapList$log_F[, (nyrs_peel + 1):nyrs]
      sd_map$mapFactor$log_F <- factor(sd_map$mapList$log_F)

      newmod$sdrep <- tryCatch({
        report_mod <- suppressWarnings(suppressMessages(
          .refit_like(
            data_list        = data_list,
            inits            = newmod$estimated_params,
            map              = sd_map,
            estimateMode     = 3,   # build only; do not re-estimate
            getsd            = FALSE,
            srr_mse_switchyr = min(data_list$srr_mse_switchyr, endyr_peel),
            srr_hat_endyr    = min(data_list$srr_hat_endyr, endyr_peel),
            suit_endyr       = pmin(data_list$suit_endyr, endyr_peel))))
        TMB::sdreport(report_mod$obj)
      }, error = function(e) {
        # A peel that cannot be reported is not a peel that failed: keep the fit
        # and its point estimates, and drop the band.
        #
        # NULL, not `newmod$sdrep`. At this point that is the FORECAST refit's
        # sdreport, built with the whole hindcast pinned -- every hindcast
        # standard error in it is exactly zero, which is the thing this report
        # pass exists to replace. Returning it hands back a zero-width band as
        # if it were real. NULL is the state `getsd = FALSE` already produces,
        # so vcov() and the plotters give NA bands rather than a false one.
        #
        # Collected rather than warned here: this runs inside the per-peel
        # closure, and a warning raised in a FORK/PSOCK worker is discarded, so
        # at the default `cores` nobody was ever told. Re-raised by the caller
        # after dispatch, as self_test() does with sim_warns.
        .peel_warnings <<- c(.peel_warnings, paste0(
          "Peel ", i, ": could not report hindcast standard errors (",
          conditionMessage(e), "), so this peel has no uncertainty band. ",
          "Point estimates are unaffected."))
        NULL
      })
    }

    # gc()
    #
    # map$mapFactor <- map$mapFactor[names(newmod$map$mapFactor)]
    # check <- c()
    # check_na <- c()
    # for(j in 1:length(map$mapList)){
    #   check[j] <- sum(map$mapFactor[[j]] != newmod$map$mapFactor[[j]], na.rm = TRUE)
    #   check_na[j] <- sum(is.na(map$mapFactor[[j]]) != is.na(newmod$map$mapFactor[[j]]), na.rm = TRUE)
    # }

    # * Report the peel's own terminal year ----
    # Set here, AFTER both refits, and never before: `endyr` sizes the model,
    # and the forecast refit above turns F back on over `(nyrs_peel+1):nyrs`
    # against the FULL nyrs, so peeling it earlier would index off the end of
    # log_F. At this point it is output metadata only.
    #
    # Every plot builds its year axis per model as `styr:endyr`
    # (`R/7-plot_ceattle.R`), and nothing outside this file reads `endyr_peel`.
    # Without this each peel was drawn to the full model's terminal year, so
    # the peels were indistinguishable -- the opposite of what a retrospective
    # plot is for.
    #
    # Mohn's rho is unaffected: it reads `endyr_peel` off each peel and the
    # full model's `endyr` from this function's enclosing scope, never a peel's
    # `data_list$endyr`.
    #
    # Note this makes the returned peel deliberately inconsistent: its
    # parameters, quantities, and `catch_data` still span the full hindcast,
    # because the peeled years are its retrospective FORECAST. `endyr` marks
    # what was fit, not what was estimated. Plot with `incl_proj = TRUE` to see
    # the forecast years, and read `endyr_full` (carried through from the source
    # model) for where those forecast years end.
    newmod$data_list$endyr <- endyr_peel

    # Return model only if BOTH refits converged, else NULL (dropped
    # post-dispatch)
    if (hindcast_converged && .refit_converged(newmod)) {
      if (length(.peel_warnings))
        attr(newmod, "peel_warnings") <- .peel_warnings
      return(newmod)
    }
    return(NULL)
  } # End run_one_peel closure


  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Dispatch peels (parallel via PSOCK or sequential) ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  if (use_parallel) {
    peel_results <- .parallel_lapply(peel_seq, run_one_peel,
                                     min(cores, length(peel_seq)), environment())
  } else {
    peel_results <- lapply(peel_seq, run_one_peel)
  }

  # Drop non-converged peels and prepend the original model
  peel_results <- peel_results[!vapply(peel_results, is.null, logical(1))]
  # A warning, not a message: rho below is averaged over the peels that
  # survive, so a drop changes the number that gets reported.
  .report_dropped(length(peel_seq) - length(peel_results), length(peel_seq),
                  "peel", warn = TRUE)

  # Re-raise whatever the peels collected. They run in a FORK/PSOCK worker,
  # whose warnings are discarded, so this is the only place the caller can hear
  # them -- the pattern self_test() uses for sim_warns. Each message already
  # names its peel, and the attribute is stripped so a returned fit carries no
  # trace of the plumbing.
  for (w in unique(unlist(lapply(peel_results, attr, "peel_warnings")))) {
    warning(w, call. = FALSE)
  }
  peel_results <- lapply(peel_results, function(m) {
    attr(m, "peel_warnings") <- NULL
    m
  })

  mod_list <- c(list(object), peel_results)

  # Mohn's rho averages the peels against the full model, so with none left the
  # sums below stay at their initialized zeros and every species column comes
  # nothing is accumulated and $mohns comes back with no rows. Say so, rather
  # than leave the caller to infer it from an empty table.
  if (length(peel_results) == 0L) {
    warning("No peel converged, so Mohn's rho cannot be computed and $mohns has no rows. ",
            "Inspect one with retrospective(..., peels = c(", peel_seq[1],
            ", ", peel_seq[1], ")) -- a repeated depth asks for that peel ",
            "alone -- and read its $convergence.", call. = FALSE)
  }


  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Calculate Mohs rho ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # * Data frame to save ----
  # LONG, one row per Object x Forecast year x species, sharing the leading
  # `Object`, `Forecast year`, `N`, `species` columns with $mase so the bias and
  # skill tables merge on those keys and read in the same terms. A species is a
  # ROW here rather than a column: $mase's N varies by species, which a column
  # per species cannot carry, and a long table is the same shape whatever nspp
  # is.
  objects <- c("biomass", "ssb", "R", "F_spp")
  spnames <- data_list[["spnames"]] %||% as.character(seq_len(data_list$nspp))

  # Accumulate the signed relative error and its count per Object x forecast
  # year x species. Keyed rather than pre-allocated by position: a peel too
  # shallow to reach a horizon contributes nothing, and a pre-allocated row it
  # never reached used to survive as Object "0" with N = 0 and a NaN rho, which
  # reads as a computed result.
  acc <- list()

  # * Loop through peels ----
  # seq_len, not 1:(n-1): with every peel dropped `mod_list` holds only the
  # original model and `1:0` would count DOWN to mod_list[[2]].
  for (i in seq_len(length(mod_list) - 1)) {
    endyr_peel <- mod_list[[i + 1]]$data_list$endyr_peel
    nyrs_peel <- mod_list[[i + 1]]$data_list$endyr_peel - styr + 1

    # * Loop output ----
    for (j in 1:length(objects)) {
      term_quantities <- mod_list[[1]]$quantities[[objects[j]]]
      retro_quantities <- mod_list[[i + 1]]$quantities[[objects[j]]]

      # * Loop forecast years ----
      for(yr in 0:nyrs_forecast){

        # If data exist for forecast (save)
        if(endyr_peel + yr <= endyr){

          # * Get full and peeled models ----
          base <- term_quantities[, nyrs_peel + yr]
          peel <- retro_quantities[, nyrs_peel + yr]
          rel_error <- ((peel - base)/base)

          # * Save and sum relative error ----
          # `N` is the number of observations rho is computed from, so a species
          # whose relative error is not finite does not count toward it. A zero
          # base -- an unfished year for `F_spp`, a collapsed stock for `ssb` or
          # `R` -- makes rel_error NaN or Inf, and summing that would make rho
          # NaN for every peel at this horizon rather than only for the peel that
          # could not be scored. Counted per SPECIES, for the same reason.
          ok <- is.finite(rel_error)
          key <- paste(objects[j], yr, sep = "\r")
          if (is.null(acc[[key]])) {
            acc[[key]] <- list(Object = objects[j], fyr = yr,
                               N = integer(data_list$nspp),
                               total = numeric(data_list$nspp))
          }
          acc[[key]]$N[ok]     <- acc[[key]]$N[ok] + 1L
          acc[[key]]$total[ok] <- acc[[key]]$total[ok] + rel_error[ok]
        }
      }
    }
  }

  # * Divide N ----
  mohns <- if (length(acc)) {
    out <- do.call(rbind, lapply(acc, function(a) data.frame(
      Object          = a$Object,
      "Forecast year" = a$fyr,
      N               = a$N,
      species         = spnames,
      # No contributing peel is no answer, not 0/0.
      rho             = ifelse(a$N > 0L, a$total / a$N, NA_real_),
      check.names = FALSE, stringsAsFactors = FALSE)))
    # Through the same exit as $mase, so both tables come back in one order.
    # `acc` is in first-touch order, which depends on which peels reached which
    # horizon -- so two runs differing only in a dropped peel would otherwise
    # return differently ordered tables.
    .rce_as_mohns_names(out)
  } else {
    # No peel reached any horizon. An empty table of the right shape, so a
    # caller's subset and merge still work rather than erroring on a missing
    # column.
    data.frame(Object = character(0), "Forecast year" = integer(0),
               N = integer(0), species = character(0), rho = numeric(0),
               check.names = FALSE, stringsAsFactors = FALSE)
  }



  # # * Beta coefficients ----
  # objects <- colnames(object$estimated_params$beta_rec_pars)
  # beta_mohns <- data.frame(matrix(0, nrow = length(objects), ncol = 3 + data_list$nspp))
  # colnames(beta_mohns) <- c("Object", "Forecast year", "N", data_list$spnames)

  #
  #   # * Loop through peels ----
  #   for (i in 1:(length(mod_list) - 1)) {
  #     endyr_peel <- mod_list[[i + 1]]$data_list$endyr_peel
  #     nyrs_peel <- mod_list[[i + 1]]$data_list$endyr_peel - styr + 1
  #     ind <- 1
  #
  #     # * Loop output ----
  #     for (j in 1:length(objects)) {
  #       base <- mod_list[[1]]$estimated_params$beta_rec_pars[,j]
  #       peel <- mod_list[[i + 1]]$estimated_params$beta_rec_pars[,j]
  #       rel_error <- ((peel - base)/base)
  #
  #       # * Save and sum relative error ----
  #       beta_mohns[j, 1] <- objects[j]        # Object
  #       beta_mohns[j, 2] <- 0                 # Year
  #       beta_mohns[j, 3] <- beta_mohns[j, 3] + 1   # N
  #       beta_mohns[j, 4:(data_list$nspp + 3) ] <- beta_mohns[j, 4:(data_list$nspp + 3)] + rel_error # Relative error
  #     }
  #   }

  # * Divide N ----
  # beta_mohns[, 4:(data_list$nspp + 3) ] <- beta_mohns[, 4:(data_list$nspp + 3)]/beta_mohns[, 3]

  mod_list <- rev(mod_list)
  # Name each entry from its own terminal year rather than assuming all `peels`
  # survived: a dropped peel then leaves a gap in the labels instead of shifting
  # every later one onto the wrong model (and `Year_...` running one short of the
  # vector, which errored). The original model carries endyr_peel = endyr, set
  # at the top of the function.
  names(mod_list) <- paste0(
    "Year_",
    vapply(mod_list, function(x) as.numeric(x$data_list$endyr_peel), numeric(1)))

  # Still the same list -- $Rceattle_list and $mohns are unchanged. The class
  # adds a print method that carries Mohn's reference band, which the bare
  # number never did; see print.Rceattle_retro().
  #
  # `peels_requested` is carried so print() can say how many were asked for. The
  # warning above is gone by the time anyone reads the object back off disk, and
  # the list length alone cannot show a drop. It is a COUNT, not the depths --
  # print() subtracts it from the number kept -- so `peel_depths` carries the
  # depths that were asked for, which a subset call cannot otherwise recover
  # once a peel has been dropped.

  # Forecast SKILL, alongside Mohn's rho's forecast BIAS ----
  # Both are per horizon, and they answer different questions: rho is a signed
  # relative error, so +/- cancel across peels and it measures bias; MASE is an
  # absolute error scaled by "assume nothing changes", so it measures skill and
  # MASE < 1 is a claim. Computed here rather than left to hindcast_skill()
  # because it costs nothing -- every input is a quantity these peels already
  # report -- and because reading bias without skill is how a model that is
  # unbiased on average and wrong every year looks fine.
  #
  # Always, and for all three quantities. On a short horizon SSB barely responds
  # to the recruitment the peels disagree about -- measured at 0.03% against
  # 46.8% in recruitment on a BS2017SS DSEM at a 3-year peel -- so scoring one
  # without the others invites reading "no difference" off the least responsive.
  #
  # reference = "model" only. The "observed" variant rebuilds every peel at
  # estimateMode = 3 and refuses on analytical catchability, which does not
  # belong in the default diagnostic; hindcast_skill() still offers it.
  #
  # IT SCORES WHATEVER `forecast_rec` ASKED FOR, and the default is "mean" --
  # which projects every peel at its own historical mean recruitment whatever
  # process the model carries. Two models scored that way forecast by the SAME
  # rule, so a difference between their tables is a difference between their
  # hindcast fits, not between their projection rules. That is a real quantity
  # and not the one most people want. The setting is therefore a COLUMN: two
  # tables can be rbind()ed and stay interpretable, and a reader cannot mistake
  # one for the other. To compare projection rules, score the same fit twice --
  # forecast_rec = "mean" against "model" -- see ?hindcast_skill.
  mase <- NULL
  if (length(peel_results)) {
    # intersect, not stop: hindcast_skill() refuses an unreported quantity
    # because the caller named it, but here the three are this function's own
    # choice and a fit that cannot report one should still return a retro. Say
    # so rather than drop it silently.
    .want <- c("ssb", "biomass", "R")
    .q <- intersect(.want, names(object$quantities))
    if (length(.q) < length(.want)) {
      message("Forecast skill scored on ", paste(.q, collapse = ", "),
              "; this fit does not report ",
              paste(setdiff(.want, .q), collapse = ", "), ".")
    }
    .by <- .rce_mase_by_year_model(object, peel_results, .q)
    if (length(.by)) {
      mase <- .rce_mase_aggregate(do.call(rbind, .by))
      # Which projection rule these peels were forecast under. A column rather
      # than an attribute because the number is unreadable without it: the
      # default is "mean", so a DSEM scored here is being told how well MEAN
      # recruitment forecasts. Last, so the key columns match $mohns.
      mase$forecast_rec <- forecast_rec
      mase <- .rce_as_mohns_names(mase)
      mase <- mase[, c(setdiff(names(mase), "forecast_rec"), "forecast_rec")]
    }
  }

  structure(list(Rceattle_list = mod_list, mohns = rbind(mohns),
                 mase = mase,
                 # The rule these peels were forecast under, recorded whether or
                 # not $mase could be computed. hindcast_skill() reads it to
                 # refuse a conflicting forecast_rec, and $mase is NULL whenever
                 # no peel converged -- exactly when a silent mismatch is most
                 # likely, so the record cannot live only on that table.
                 forecast_rec = forecast_rec,
                 peels_requested = length(peel_seq),
                 peel_depths = peel_seq),
            class = "Rceattle_retro")
}


#' Print method for a retrospective analysis
#'
#' @description Reports Mohn's rho against a reference band rather than as a
#' bare number. The default band is +/- 0.2 on SSB, the rule of thumb
#' `vignette("model-diagnostics")` states; Hurtado-Ferro et al. (2015) give the
#' asymmetric, life-history-dependent alternatives (-0.15 to 0.20 for
#' long-lived, -0.22 to 0.30 for short-lived), which is why the band is an
#' argument rather than a constant.
#'
#' Only the terminal-year peel (`Forecast year` 0) is judged. The forecast-skill
#' rows are reported for information: a rho computed over a forecast horizon is
#' not the quantity the +/- 0.2 rule was calibrated on.
#'
#' A rho averaged over a SUBSET of peels is not the quantity the band was
#' calibrated on either. The band assumes the standard contiguous set, and a
#' deep-peel-only relative error is systematically larger in magnitude than the
#' average over `1:n` in any model with real retrospective bias, so the
#' false-pass rate moves with the depth set. Whenever
#' `retrospective(peels = )` was given anything but `1:n`, the depths are
#' printed and the verdict carries a `NOTE`.
#'
#' @param x A `"Rceattle_retro"` object from [retrospective()].
#' @param band Symmetric reference band for Mohn's rho. Default `0.2`.
#' @param ... Currently unused.
#' @return `x`, invisibly.
#' @references Hurtado-Ferro, F., et al. 2015. Looking in the rear-view mirror:
#'   bias and retrospective patterns in integrated, age-structured stock
#'   assessment models. ICES J. Mar. Sci. 72:99-110.
#' @export
print.Rceattle_retro <- function(x, band = 0.2, ...) {
  m <- .rce_mohns_as_long(as.data.frame(x$mohns))
  term <- m[!is.na(m[["Forecast year"]]) & m[["Forecast year"]] == 0, , drop = FALSE]

  rho <- suppressWarnings(as.numeric(term[["rho"]]))
  # A verdict computed over no rho, from rows that DO exist, is a shape we do not
  # recognize -- not a clean bill of health. Refuse rather than print "OK".
  if (nrow(term) > 0L && length(rho) == 0L) {
    stop("This object's $mohns is not a shape print() recognizes (no `rho`, ",
         "and no species columns to read one from). Rebuild it with ",
         "retrospective().", call. = FALSE)
  }
  sev <- ifelse(is.na(rho), "OK", ifelse(abs(rho) > band, "WARN", "OK"))
  n_bad <- sum(sev == "WARN")

  # Rho is averaged over the peels that survived, so a drop changes it. The
  # list length alone cannot show one; `peels_requested` is absent on an object
  # saved before it was carried.
  n_kept <- length(x$Rceattle_list) - 1L          # the unpeeled model
  n_asked <- x$peels_requested
  n_dropped <- if (is.null(n_asked)) 0L else max(0L, n_asked - n_kept)
  if (n_dropped > 0L) sev <- c(sev, "NOTE")

  # A subset of depths is a NOTE in its own right: the band is calibrated on the
  # contiguous set, so a rho averaged over 3:10 is not the number it judges.
  # `peel_depths` is absent on an object saved before it was carried, which reads
  # as "not a subset" rather than guessing.
  depths <- x$peel_depths
  is_subset <- !is.null(depths) && !identical(as.integer(depths),
                                              seq_len(max(depths)))
  if (is_subset) sev <- c(sev, "NOTE")

  .rce_diag_header(
    "retrospective", .rce_worst(sev),
    paste0(if (n_dropped > 0L) paste0(n_kept, " of ", n_asked, " peel(s); ")
           else paste0(n_kept, " peel(s); "),
           .rce_n_of(n_bad, length(rho)),
           " terminal Mohn's rho outside +/-", band))

  if (n_dropped > 0L) {
    cat("  ", n_dropped, " peel(s) dropped as non-converged, so rho is ",
        if (n_kept > 0L) paste("averaged over", n_kept) else "undefined",
        "\n", sep = "")
  }
  if (is_subset) {
    cat("  peel depth", if (length(depths) > 1L) "s" else "", " ",
        paste(depths, collapse = ", "),
        " -- rho is averaged over these only; the +/-", band,
        " band assumes the standard contiguous set\n", sep = "")
  }

  if (nrow(term)) {
    # One row per quantity x species now that $mohns is long, so the tag is per
    # row rather than per quantity across species columns.
    show <- term
    show$.tag <- vapply(seq_len(nrow(show)), function(i) {
      r <- suppressWarnings(as.numeric(show[["rho"]][i]))
      .rce_sev_tag(if (!is.na(r) && abs(r) > band) "WARN" else "OK")
    }, character(1))
    .rce_diag_table(show, stats::setNames(
      c(".tag", "Object", "species", "N", "rho"),
      c(" ", "quantity", "species", "N", "rho")))
  }
  if (any(!is.na(m[["Forecast year"]]) & m[["Forecast year"]] > 0)) {
    cat("  rows above forecast year 0 are BIAS at that horizon; the +/-", band,
        "rule is for the terminal peel only\n")
  }
  # Bias without skill is how a model that is unbiased on average and wrong
  # every year reads as fine, so say the skill table is there. Naming
  # forecast_rec matters: under the default "mean" every model forecasts by the
  # same rule, so these numbers compare FITS, not projection rules.
  if (!is.null(x$mase) && nrow(x$mase)) {
    .fr <- unique(x$mase[["forecast_rec"]])
    cat("  forecast SKILL (MASE) is in $mase, per forecast year, ",
        "forecast_rec = \"", .fr[1], "\"\n", sep = "")
    if (identical(.fr[1], "mean")) {
      cat("    that is the MEAN-recruitment projection, the same rule for ",
          "every model:\n    re-run with forecast_rec = \"model\" to score ",
          "the model's own\n", sep = "")
    }
    .thin <- x$mase[!is.na(x$mase[["mase"]]) & x$mase[["N"]] < 2L, ]
    if (nrow(.thin)) {
      cat("    ", nrow(.thin), " row(s) rest on a single peel (the deepest ",
          "horizons); read N\n", sep = "")
    }
  }
  invisible(x)
}


#' Read a `$mohns` table of either shape as the long one
#'
#' @description
#' `$mohns` was WIDE before 5.23.0.9007 -- `Object`, `Forecast year`, `N`, then
#' one column per species -- and retrospectives are saved to disk and reloaded
#' months later, so `print()` still meets that shape. Reshaped here rather than
#' refused, because the alternative is a diagnostic that cannot read its own
#' saved output.
#'
#' A table already carrying `rho` is returned untouched.
#'
#' @param m a `$mohns` data frame, either shape.
#' @return the long form: `Object`, `Forecast year`, `N`, `species`, `rho`.
#' @noRd
.rce_mohns_as_long <- function(m) {
  if (!is.data.frame(m) || "rho" %in% names(m)) return(m)
  keys <- c("Object", "Forecast year", "N")
  spp <- setdiff(names(m), keys)
  if (!all(keys %in% names(m)) || length(spp) == 0L) return(m)
  out <- do.call(rbind, lapply(spp, function(sp) data.frame(
    m[, keys, drop = FALSE], species = sp,
    rho = suppressWarnings(as.numeric(m[[sp]])),
    check.names = FALSE, stringsAsFactors = FALSE)))
  rownames(out) <- NULL
  out
}




#' Jitter analysis
#'
#' @description Refits the Rceattle model from starting values perturbed by N(0, sd) around the model's initial (pre-fit) parameters, to check convergence robustness.
#'
#' @note Attaching Rceattle masks \code{\link[base]{jitter}}, so \code{jitter(x)}
#'   on a numeric vector reaches this function and reports a missing model rather
#'   than adding noise to \code{x}. Call \code{base::jitter()} explicitly for the
#'   base-graphics behaviour.
#'
#' @inheritParams rceattle-refit-args
#' @param njitter the number of jitters to run
#' @param sd standard deviation for jitter (default = 0.2)
#' @param phase as in \code{\link{fit_mod}} default = FALSE. Jitters restart from
#'   perturbed \emph{starting} values, so a model that needed phasing to fit its
#'   real data needs it here too; leave this at \code{FALSE} for such a model and
#'   the jitters end far from any optimum and are dropped as non-converged, which
#'   reads as multimodality rather than as an unphased fit.
#' @param seed random number seed. Each jitter \code{i} uses \code{seed + i}
#'   so results are reproducible under both sequential and parallel execution.
#' @param getsd whether each jitter runs \code{TMB::sdreport}. Jitter compares
#'   objectives and point estimates across starts, so \code{FALSE} is faster
#'   with no effect on that comparison. Default \code{NULL} inherits the input
#'   model's setting (\code{TRUE} only if it carries an \code{sdrep}).
#' @param timeout elapsed-second limit per jitter, \code{Inf} (default) for none.
#'   A jitter is a deliberately perturbed start and the optimizer runs with no
#'   iteration cap, so this is the diagnostic most likely to send one somewhere
#'   pathological and stall the whole run -- a hang no convergence check can
#'   catch, because the fit never returns. One that exceeds the limit is stopped,
#'   counted as non-converged and reported separately. Approximate: the limit is
#'   checked when control returns to R, so it fires between the optimizer's
#'   function evaluations rather than inside one.
#'
#' @return a list of 1. \code{Rceattle_list}, the converged jitters, and
#'   2. \code{nll}, their objective values. Non-converged (or timed-out) starts
#'   are dropped and reported in a message, so both can be shorter than
#'   \code{njitter} -- and that count is itself the result, since the whole point
#'   is what fraction of random starts reach the same optimum.
#'
#' @examples
#' \donttest{
#' data(BS2017SS)
#' ss_run <- fit_mod(data_list = BS2017SS,
#'     inits = NULL, file = NULL,
#'     estimateMode = 0, random_rec = FALSE,
#'     msmMode = 0, avgnMode = 0,
#'     phase = FALSE, verbose = 0)
#' jitters <- jitter(ss_run, njitter = 10)
#' }
#' @export
jitter <- function(object = NULL, njitter = 50, sd = 0.2, phase = FALSE, seed = 123, cores = NULL, getsd = NULL, timeout = Inf, fit_control = NULL, Rceattle = NULL) {
  # `Rceattle` was the old name for `object`; see R/0-deprecate.R.
  if (!missing(Rceattle))
    object <- .rce_deprecated_arg(Rceattle, !missing(object), "Rceattle", "object", "jitter")
  # Cleared once consumed: .parallel_lapply() exports every binding in this frame
  # to the PSOCK workers, and two names for one fitted model send it twice.
  Rceattle <- NULL

  if (!inherits(object, "Rceattle")) {
    stop("Object is not of class 'Rceattle'")
  }

  # A DSEM needs nothing special here, unlike the peel in retrospective(). The
  # perturbation below only touches entries whose map is not NA, so under a DSEM
  # it jitters the path coefficients and the latent recruitment-deviation
  # columns while leaving the covariate columns of dsem_x_tj alone -- those are
  # mapped out because, with family = "fixed", they ARE the environmental data
  # and jittering them would perturb the data rather than the starting values.
  # And no map is passed to .refit_like() below, so fit_mod() rebuilds one and
  # fills in the DSEM blocks itself.

  # Jitters inherit the input model's sdreport setting unless overridden;
  # multimodality is judged from objectives and point estimates, not sdrep.
  if (is.null(getsd)) getsd <- !is.null(object$sdrep)

  # fit_control() overrides `phase` and `getsd`, but only where the caller
  # changed them from fit_control()'s own defaults; see .rce_refit_control().
  ctl <- .rce_refit_control(fit_control, "jitter")
  if (!is.null(ctl$phase)) phase <- ctl$phase
  if (!is.null(ctl$getsd)) getsd <- ctl$getsd

  # Cross-platform parallel via parallel::parLapply on a PSOCK cluster
  # (same approach as run_mse). Respect the CRAN core limit
  # ('_R_CHECK_LIMIT_CORES_' is set during R CMD check;
  # parallel::makeCluster errors if we exceed 2 cores then).
  chk <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  cran_cap <- nzchar(chk) && !chk %in% c("false", "0", "no")
  if (is.null(cores)) {
    cores <- if (cran_cap) 2L else max(1L, parallel::detectCores() - 6L)
  } else {
    cores <- max(1L, as.integer(cores))
    if (cran_cap) cores <- min(cores, 2L)
  }
  use_parallel <- njitter > 1L && cores > 1L

  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Per-jitter closure ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  run_one_jitter <- function(i) {

    set.seed(seed + i) # unique seed per jitter for reproducibility under parallel

    # * Adjust initial values ----
    inits <- object$initial_params
    mapList <- object$map$mapList
    data_list <- object$data_list

    for(j in 1:length(inits)){
      par <- names(inits)[j]
      inits[[j]] <- replace(inits[[j]],
                            values = ifelse(is.na(as.numeric(mapList[[par]])),
                                            as.numeric(inits[[j]]),
                                            as.numeric(inits[[j]]) + stats::rnorm(length(as.numeric(inits[[j]])), 0, sd))
      )
    }


    # * Refit ----
    # Bounded and error-trapped: a jitter is a deliberately perturbed start, so
    # it is the diagnostic most likely to send the optimizer somewhere
    # pathological -- and one such start must not abort the other njitter - 1.
    newmod <- .refit_with_timeout(
      suppressMessages(
        suppressWarnings(
          # Refit from the jittered starting values (map rebuilt from scratch).
          .refit_like(
            data_list        = data_list,
            inits            = inits,
            estimateMode     = ifelse(data_list$estimateMode < 3, 0, data_list$estimateMode),
            phase            = phase,
            getsd            = getsd,
            srr_mse_switchyr = min(data_list$srr_mse_switchyr, data_list$endyr),
            suit_endyr       = pmin(data_list$suit_endyr, data_list$endyr))
        )
      ),
      timeout = timeout)

    # Report the verdict beside the model; the dispatcher filters below.
    if (inherits(newmod, "condition")) {
      return(list(model = NULL, converged = FALSE, error = conditionMessage(newmod)))
    }
    list(model = newmod, converged = .refit_converged(newmod))
  } # End run_one_jitter closure


  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Dispatch jitters (parallel via PSOCK or sequential) ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  if (use_parallel) {
    mod_list <- .parallel_lapply(1:njitter, run_one_jitter, min(cores, njitter), environment())
  } else {
    mod_list <- lapply(1:njitter, run_one_jitter)
  }

  # Drop non-converged
  converged <- vapply(mod_list, function(x) isTRUE(x$converged), logical(1))
  errs      <- vapply(mod_list, function(x) x$error %||% NA_character_, character(1))
  mod_list  <- lapply(mod_list, function(x) x$model)[converged]
  .report_dropped(sum(!converged), length(converged), "jitter")
  .report_errors(errs, "jitter")

  # vapply, not sapply: over an empty list sapply returns list(), and the
  # documented `min(jitters$nll)` / `hist(log(nll - min(nll)))` then errors
  # rather than reporting that nothing converged. numeric(0) keeps `nll`
  # numeric at every length.
  jnll <- vapply(mod_list, function(x) as.numeric(x$opt$objective), numeric(1))
  if (length(mod_list) > 0) {
    names(mod_list) <- paste0("Jitter_", seq_along(mod_list))
  }


  # Return ----
  # Still the same list -- $Rceattle_list and $nll are unchanged. `njitter` is
  # carried so the print method can report the fraction of STARTS that reached
  # the optimum: the converged runs alone cannot say how many were attempted,
  # and that fraction is the result jitter() exists to produce.
  structure(list(Rceattle_list = mod_list, nll = unname(jnll),
                 njitter = njitter),
            class = "Rceattle_jitter")
}


#' Print method for a jitter analysis
#'
#' @description Reports what the run was for: how many random starts reached the
#' best optimum found. The objective values alone cannot say that -- non-converged
#' starts are dropped before the result is returned, so the count of returned
#' fits is not the count attempted.
#'
#' @param x A `"Rceattle_jitter"` object from [jitter()].
#' @param tol Objective units within which a start counts as reaching the same
#'   optimum. Default `0.01`, which is far below any difference that would change
#'   management advice and far above optimizer noise.
#' @param ... Currently unused.
#' @return `x`, invisibly.
#' @export
print.Rceattle_jitter <- function(x, tol = 0.01, ...) {
  nll <- x$nll[is.finite(x$nll)]
  tried <- x$njitter %||% NA_integer_
  best <- if (length(nll)) min(nll) else NA_real_
  at_best <- if (length(nll)) sum(nll - best <= tol) else 0L

  # A start that did not converge is as informative as one that landed high:
  # both say the optimum is hard to reach from an arbitrary point.
  sev <- if (!length(nll)) "FAIL"
         else if (!is.na(tried) && at_best < 0.9 * tried) "WARN"
         else if (!is.na(tried) && length(nll) < tried) "NOTE"
         else "OK"

  .rce_diag_header(
    "jitter", sev,
    paste0(at_best, " of ", if (is.na(tried)) length(nll) else tried,
           " start(s) reached the best optimum (within ", tol, ")",
           if (!is.na(tried) && length(nll) < tried)
             paste0("; ", tried - length(nll), " did not converge") else ""))

  if (length(nll)) {
    cat("  best -log L : ", formatC(best, format = "f", digits = 4), "\n", sep = "")
    if (at_best < length(nll)) {
      worse <- sort(nll[nll - best > tol])
      cat("  higher optima found at +",
          paste(formatC(utils::head(worse, 5) - best, format = "g", digits = 3),
                collapse = ", +"),
          if (length(worse) > 5) ", ..." else "", "\n", sep = "")
      cat("  a higher optimum means the reported fit is start-dependent\n")
    }
  }
  invisible(x)
}


#' Which fleets' selectivity deviations are scored by an UNNORMALIZED kernel.
#'
#' A random effect's peeled tail can be freed only where the density scoring it is
#' normalized. The Laplace integral over a data-free tail is then exactly 1, so
#' freeing it removes the fabricated "the deviation was exactly zero" evidence and
#' changes nothing else. Several selectivity forms carry the ADMB/AMAK lineage's
#' bare sums of squares instead, with no dnorm constant: marginalizing one of
#' those contributes -k * log(sel_dev_sd) to the objective, which drives the
#' estimated SD UP with peel depth -- the same bias as pinning, in the other
#' direction. Those fleets stay pinned.
#'
#' From ceattle.cpp slots 4-5, by Selectivity form:
#'
#'   sel_coff_dev   normalized  Hake -- dnorm(dev, 0, sel_dev_sd);
#'                              2DAR1 -- SCALE(SEPARABLE(AR1, AR1));
#'                              3DAR1 -- GMRF(Q)
#'                  bare SSQ    NonParametric and NonParametricPM (decreasing,
#'                              curvature, random-walk and dev-magnitude
#'                              penalties); LogisticPM (weighted first difference
#'                              of realized log selectivity)
#'
#'   log_sel_slp_dev, sel_inf_dev
#'                  normalized  every form that estimates them -- the IID and
#'                              random-walk branches are all dnorm(..., true)
#'                  bare SSQ    LogisticPM only, whose free age-1 deviate random
#'                              walk is sel_curve_pen * (dev_t - dev_{t-1})^2
#'
#' Decided per Selectivity_index GROUP, not per fleet. Fleets sharing an index
#' share ONE parameter block through shared map levels, so pinning one member
#' while freeing another would leave the pinned fleet's cells at their starting
#' value while the shared level moved -- two fleets that must have identical
#' selectivity would silently stop sharing it in the peeled years.
#'
#' @param fleet_control the model's `fleet_control` table.
#' @return A list of integer fleet-row vectors: `sel_coff_dev` and `limb`.
#' @noRd
.rce_sel_unnormalized_rows <- function(fleet_control) {
  sel_map <- .rce_allowed_map("Selectivity")

  # The forms whose sel_coff_dev density is normalized, and the one form whose
  # limb deviates are not.
  coff_normalized <- c("Hake", "2DAR1", "3DAR1")
  limb_unnormalized <- "LogisticPM"

  # Assert rather than intersect quietly. A renamed or retired form would leave
  # an empty code set, and an empty set here reads as "nothing is unnormalized"
  # -- it would free a tail that must stay pinned, silently, which is the defect
  # this function exists to prevent.
  need <- c(coff_normalized, limb_unnormalized)
  if (!all(need %in% names(sel_map))) {
    stop("Internal: the Selectivity map no longer names ",
         paste(setdiff(need, names(sel_map)), collapse = ", "),
         "; retrospective() cannot tell which selectivity deviations are ",
         "scored by a normalized density.", call. = FALSE)
  }

  # fleet_control$Selectivity holds names on a workbook that has been through
  # revert_switches() and integer codes on one that has not, so read both.
  .is_form <- function(nms) {
    codes <- unname(sel_map[nms])
    s  <- fleet_control$Selectivity
    si <- suppressWarnings(as.integer(as.character(s)))
    (!is.na(si) & si %in% codes) | (as.character(s) %in% nms)
  }

  n   <- nrow(fleet_control)
  grp <- fleet_control$Selectivity_index
  # A fleet with no Selectivity_index shares with nobody, so it is its own
  # group. Same idiom as .sel_start_year_by_group() in rearrange_data().
  grp <- if (is.null(grp)) paste0("_row", seq_len(n)) else
    ifelse(is.na(grp), paste0("_row", seq_len(n)), as.character(grp))

  widen <- function(bad) which(grp %in% unique(grp[bad]))

  list(
    # Anything that is not one of the three normalized forms. A fleet whose form
    # does not use sel_coff_dev at all is included and is a no-op: build_map()
    # has already mapped those cells out.
    sel_coff_dev = widen(!.is_form(coff_normalized)),
    limb         = widen(.is_form(limb_unnormalized))
  )
}


#' The map a retrospective peel runs under: which deviations it holds fixed over
#' the years the peel did not see.
#'
#' Turn off forecasted parameters -- but NOT the ones that are random effects. A
#' pinned deviation is still scored by its density, and "the deviation was exactly
#' zero" is the strongest possible evidence for a small process SD; leaving a
#' random effect free lets the Laplace approximation integrate it out instead,
#' which is what no data should mean. The NEWS entry for 5.15.0 carries the
#' measured bias (-6.6% on sigma at 5 peels, monotone in peel depth, so it becomes
#' a trend in the quantity Mohn's rho measures).
#'
#' A block is freed only if it is BOTH a random effect and scored by a NORMALIZED
#' density: the Laplace integral over a data-free tail is exactly 1 for a
#' normalized Gaussian, so freeing it removes the fabricated term and nothing
#' else. Several selectivity forms are bare sums of squares with no dnorm
#' constant, and those stay pinned -- but PER FLEET, not per block, because the
#' exempt forms are a property of the fleet; see .rce_sel_unnormalized_rows().
#'
#' @param map the parent fit's `map` (`mapList` + `mapFactor`).
#' @param random_vars `fit$random_vars`, recorded by fit_mod() at the HINDCAST
#'   build. This is the one source of truth -- `obj$env$random` is not: under
#'   estimateMode = 0 with a harvest control rule that object is the PROJECTION
#'   object, whose map turns every hindcast entry off, so its `random` declaration
#'   is empty and every block reads as a fixed effect.
#'
#'   `character(0)` means the fit declared no random effects, and every
#'   deviation is then pinned because that is what the fit itself did. `NULL`
#'   means the fit never recorded it (made before 5.10.0, or at an
#'   `estimateMode` that does not store it) and is treated the same way, but it
#'   is NOT the same situation: there, pinning may be wrong and shrinks the
#'   process SDs with peel depth. Reporting that is `retrospective()`'s job --
#'   it checks `is.null()` once before dispatching its peels, because a warning
#'   raised in here reaches no one when the peels run on a parallel cluster.
#' @param fleet_control the model's fleet table, for the selectivity forms.
#' @param nyrs_peel,nyrs,nyrs_proj retained hindcast years, total hindcast years,
#'   and total years including the projection.
#' @return `map`, with the peeled tail of each pinned block set to NA.
#' @noRd
.rce_peel_map <- function(map, random_vars, fleet_control,
                          nyrs_peel, nyrs, nyrs_proj) {
  # Silent: `retrospective()` has already reported an unrecorded `random_vars`
  # (see the @param note). Warning here instead reached no one, because the
  # peels run on a parallel cluster at the default `cores` and a worker's
  # warnings are discarded.
  if (is.null(random_vars)) random_vars <- character(0)
  pin <- function(nm) !(nm %in% random_vars)
  has <- function(nm) !is.null(map$mapList[[nm]])
  relevel <- function(nm) factor(map$mapList[[nm]])

  peeled_hind <- (nyrs_peel + 1):nyrs
  peeled_all  <- (nyrs_peel + 1):nyrs_proj

  # rec_dev is [species, year]; log_M1_dev is [species, sex, age, year]; both run
  # through the projection. index_q_dev is [fleet, year], hindcast only.
  if (pin("rec_dev") && has("rec_dev")) {
    map$mapList$rec_dev[, peeled_all] <- NA
    map$mapFactor$rec_dev <- relevel("rec_dev")
  }
  if (pin("log_M1_dev") && has("log_M1_dev")) {
    map$mapList$log_M1_dev[, , , peeled_all] <- NA
    map$mapFactor$log_M1_dev <- relevel("log_M1_dev")
  }
  if (pin("index_q_dev") && has("index_q_dev")) {
    map$mapList$index_q_dev[, peeled_hind] <- NA
    map$mapFactor$index_q_dev <- relevel("index_q_dev")
  }

  # Selectivity is decided PER FLEET. One model can carry a 3DAR1 fleet, whose
  # deviations are a proper GMRF and must be freed, beside a NonParametricPM
  # fleet, whose penalty is a bare SSQ and must stay pinned; pinning the whole
  # block for the second fleet's sake reintroduces the shrinkage on the first.
  #
  # The selectivity dimension is indexed by fleet ROW -- Fleet_code equals the
  # row number, and build_params() dimensions these arrays by
  # nrow(fleet_control).
  sel_un <- .rce_sel_unnormalized_rows(fleet_control)
  # dim() is NULL for a block this model does not carry, which is not an error --
  # there is simply nothing to pin. Only a PRESENT block whose fleet dimension
  # disagrees with fleet_control is one, because the row indices would then pin
  # the wrong fleets.
  n_sel <- dim(map$mapList$sel_coff_dev)[1]
  n_flt <- nrow(fleet_control)
  if (!is.null(n_sel) && !identical(as.integer(n_sel), as.integer(n_flt))) {
    stop("sel_coff_dev has ", n_sel, " selectivity rows but fleet_control has ",
         n_flt, "; the peel cannot tell which fleet's deviations are scored by ",
         "which density.", call. = FALSE)
  }
  all_sel <- if (is.null(n_sel)) integer(0) else seq_len(n_sel)

  # The limb deviates are [limb, selectivity, sex, year], so the fleet rows index
  # the SECOND dimension.
  slp_rows <- if (pin("log_sel_slp_dev")) all_sel else sel_un$limb
  if (length(slp_rows) && has("log_sel_slp_dev")) {
    map$mapList$log_sel_slp_dev[, slp_rows, , peeled_hind] <- NA
    map$mapFactor$log_sel_slp_dev <- relevel("log_sel_slp_dev")
  }
  inf_rows <- if (pin("sel_inf_dev")) all_sel else sel_un$limb
  if (length(inf_rows) && has("sel_inf_dev")) {
    map$mapList$sel_inf_dev[, inf_rows, , peeled_hind] <- NA
    map$mapFactor$sel_inf_dev <- relevel("sel_inf_dev")
  }
  # sel_coff_dev is [selectivity, sex, bin, year] -- rows index the FIRST.
  coff_rows <- if (pin("sel_coff_dev")) all_sel else sel_un$sel_coff_dev
  if (length(coff_rows) && has("sel_coff_dev")) {
    map$mapList$sel_coff_dev[coff_rows, , , peeled_hind] <- NA
    map$mapFactor$sel_coff_dev <- relevel("sel_coff_dev")
  }

  map
}
