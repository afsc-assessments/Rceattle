#' Retrospective peels
#'
#' @description Refit an Rceattle model with its last 1, 2, ..., `peels` years
#'   of data removed, and calculate Mohn's rho for the terminal year and the
#'   `nyrs_forecast` years after it. Each peel is fitted twice: its hindcast,
#'   then F in the removed years, fitted to their observed catch with the
#'   hindcast held fixed.
#'
#' @inheritParams rceattle-refit-args
#' @param peels number of retrospective peels (default 5).
#' @param rescale which `env_data` columns to standardize on each peel's
#'   retained years: `FALSE` (default) none, `TRUE` every numeric column except
#'   `Year`, random-effect grouping columns and columns a linkage reads as a
#'   category, or a character vector of column names.
#' @param nyrs_forecast number of years after each peel's terminal year for
#'   which Mohn's rho is also calculated (default 3).
#' @param getsd whether each peel reports standard errors; `NULL` (default) does
#'   so only if `object` holds an `sdrep`.
#' @param phase whether each peel is refitted in phases (default `TRUE`); a peel
#'   starts from the unpeeled model's estimates.
#' @param projected_env_data_known `env_data` columns whose values after each
#'   peel's terminal year were known then, such as a climate projection, so the
#'   peels keep them (a categorical linkage term is held regardless): `FALSE`
#'   (default) none, `TRUE` all, or a character vector of column names.
#'
#' @details
#' **Hindcast and forecast years.** A peel's hindcast is the model as it would
#' have been fitted at its terminal year: it ends there, so no later index,
#' composition or catch, and no parameter for a later year, enters its
#' estimates. It starts from the unpeeled estimates for its own years and uses
#' the unpeeled map cut to them. Its removed years are then its retrospective
#' forecast: F is fitted to the observed catch with the hindcast held, and q,
#' selectivity, M1, weight-at-age, empirical selectivity and ration stay at the
#' terminal year (a selectivity form that stores each year's step takes a step
#' of 0). Recruitment takes the deviation [sample_rec()] uses with
#' `sample_rec = FALSE`, computed from the peel's fit (for a penalty-form curve
#' with no penalty years in the peel, the mean over its own years after the
#' first). Each peel's
#' `data_list` holds `endyr` and `endyr_peel` (its terminal year), `endyr_full`
#' (the unpeeled terminal year, where the forecast ends) and `projyr` (the end
#' of the projection). The forecast years are
#' `endyr_peel + seq_len(endyr_full - endyr_peel)`, which is empty for the
#' unpeeled model, where `(endyr_peel + 1):endyr_full` would count down.
#'
#' **Dropped and flagged peels.** A peel whose hindcast did not converge is
#' dropped, with a warning; under `getsd = TRUE` that includes a hindcast
#' Hessian that is not positive definite. A peel whose forecast refit did not
#' fit F to the catch (gradient check, or a peeled F at a bound) is kept but
#' left out of the `$mohns` rows with `Forecast year` above 0, and
#' `attr(peel, "forecast_refit_error")` gives the reason. Check it before
#' reading that peel's forecast years: only `plot_timeseries()` and its
#' wrappers warn about them.
#'
#' **Standard errors.** Mohn's rho reads only point estimates, so
#' `getsd = FALSE` is faster. Under `getsd = TRUE` each peel is rebuilt at its
#' estimates with the hindcast parameters treated as estimated, so no estimate
#' moves. A fitted year's band is the peeled hindcast's own, except a depletion
#' relative to last-projection-year biomass (`HCR = "NoFishing"`,
#' `msmMode > 0`, `DynamicHCR = FALSE`), which reads the peeled F. A forecast
#' year's band treats the catch as known and carries the hindcast's uncertainty
#' through the F needed to take it, without recruitment variability, so it is
#' narrower than a prediction interval. If the rebuild fails (an error, a
#' Hessian that is not positive definite, or no finite standard error), the
#' peel is kept, but plots, `as.data.frame()` and `vcov()` ignore its `sdrep`,
#' which may not hold valid standard errors; `attr(peel, "sd_report_error")`
#' gives the reason. `print()` counts the peels carrying either attribute;
#' `retrospective()` warns about neither.
#'
#' **Environmental data.** A peel's forecast should read only what was known at
#' its terminal year, so after that year a peel treats each `env_data` column by
#' how the model reads it:
#'
#' * an observed covariate (`linkage_spec(observe = )`) loses its observations;
#' * a column the model reads directly (a `Cindex` temperature, an environmental
#'   catchability or M1) takes the retained years' mean, as in any missing year;
#' * a numeric fixed-effect linkage covariate, a 0/1 block column included, takes
#'   the retained years' mean, with a warning; a term computed over the whole
#'   column, such as `scale()`, then changes in the retained years too;
#' * a fixed-effect category, such as `cut(Year, ...)` or `factor(regime)`, holds
#'   its terminal-year value ([linkage_spec()]), though one cut from the whole
#'   column, such as `cut(temp, 3)`, still sets its breaks from the whole column,
#'   post-peel years included;
#' * a random-effect grouping column, a column with no retained value, and a
#'   fixed-effect covariate that is also observed, also read as a category (the
#'   category term is still held), or logical keep their post-peel values, with a
#'   warning; Mohn's rho is then conditional on them;
#' * a column listed in `projected_env_data_known` keeps its values with no
#'   warning, though a category is still held.
#'
#' Refitting a returned peel, with [fit_mod()] or a diagnostic on it, no longer
#' holds its categories, while its numeric covariates stay mean-filled.
#'
#' `rescale` standardizes a column on the retained years and applies that centre
#' and scale to the whole column (a column with fewer than two distinct retained
#' values stays in its units, with a warning). Name the columns when `env_data`
#' also holds an input that is not a covariate, such as a `Cindex` temperature in
#' degrees C, which `TRUE` would turn into z-scores. The unpeeled model is not
#' rescaled, so standardize those columns the same way before fitting it.
#'
#' **Catchability.** A peel's hindcast uses the unpeeled model's map, cut to its
#' years, so a survey whose index falls only in the removed years keeps a free q
#' with no index to inform it, which under `getsd = TRUE` can leave the Hessian
#' not positive definite and drop the peel.
#'
#' @return An `Rceattle_retro` list:
#'   \describe{
#'     \item{\code{Rceattle_list}}{each kept peel, deepest first, then the
#'       unpeeled model, named for its terminal year (e.g. \code{Year_2016});
#'       index by name, since a dropped peel leaves a gap.}
#'     \item{\code{mohns}}{Mohn's rho by \code{Object} and
#'       \code{Forecast year}, one column per species, with \code{N} the number
#'       of peels averaged (\code{NaN} where \code{N} is 0; the call warns only
#'       if no peel converged).}
#'     \item{\code{peels_requested}}{the \code{peels} asked for.}
#'   }
#'
#' @examples
#' \donttest{
#' data(BS2017SS)
#' ss_run <- fit_mod(data_list = BS2017SS, inits = NULL, file = NULL,
#'                   estimateMode = 0, random_rec = FALSE, msmMode = 0,
#'                   fit_control = fit_control(phase = FALSE, verbose = 0))
#' retro <- retrospective(ss_run, peels = 10)
#' }
#' @export
retrospective <- function(object = NULL, peels = 5, rescale = FALSE, nyrs_forecast = 3, cores = NULL, getsd = NULL, phase = TRUE, fit_control = NULL, projected_env_data_known = FALSE, Rceattle = NULL) {
  # Newer arguments follow the older ones so positional calls keep working; the
  # deprecated `Rceattle` sits last (R/0-deprecate.R).
  if (!missing(Rceattle))
    object <- .rce_deprecated_arg(Rceattle, !missing(object), "Rceattle", "object", "retrospective")
  # Cleared once consumed: .parallel_lapply() exports every binding in this frame
  # to the PSOCK workers, and two names for one fitted model send it twice.
  Rceattle <- NULL

  if (!inherits(object, "Rceattle")) {
    stop("Object is not of class 'Rceattle'")
  }

  # By default, peels report SEs only if `object` holds an sdrep. FALSE is faster
  # (rho reads point estimates) and keeps a peel TRUE would drop for a non-PD
  # Hessian.
  if (is.null(getsd)) getsd <- !is.null(object$sdrep)

  ctl <- .rce_refit_control(fit_control, "retrospective")
  if (!is.null(ctl$phase)) phase <- ctl$phase
  if (!is.null(ctl$getsd)) getsd <- ctl$getsd

  # Get objects
  object$data_list$endyr_peel <- object$data_list$endyr
  # Where each peel's forecast ends, since a peel reports its own terminal year
  # as `endyr`. Set once: run_one_peel copies this data_list into every peel.
  object$data_list$endyr_full <- object$data_list$endyr
  data_list <- object$data_list # used by Mohn's rho block below
  endyr <- object$data_list$endyr
  styr <- object$data_list$styr
  nyrs <- length(styr:endyr)
  projyr <- object$data_list$projyr
  nyrs_proj <- projyr - styr + 1

  # FORK cluster where available, PSOCK on Windows (.parallel_lapply()); R CMD
  # check sets _R_CHECK_LIMIT_CORES_ and makeCluster() errors above 2 cores.
  chk <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  cran_cap <- nzchar(chk) && !chk %in% c("false", "0", "no")
  if (is.null(cores)) {
    cores <- if (cran_cap) 2L else max(1L, parallel::detectCores() - 6L)
  } else {
    cores <- max(1L, as.integer(cores))
    if (cran_cap) cores <- min(cores, 2L)
  }
  use_parallel <- peels > 1L && cores > 1L

  # Stored-map divergence ----
  # A peel's hindcast uses `object$map`, cut to its years, so rho compares like
  # with like. Warn if build_map() on the ORIGINAL inputs now differs (an edited
  # map or newer map code, not a shorter peel); a worker's warning is discarded.
  .map_drift <- tryCatch({
    # Pass random_sel: build_map() defaults it to FALSE with no data_list
    # fallback, so a fresh map would map out a random_sel fit's sel_dev_log_sd.
    .fresh <- suppressWarnings(suppressMessages(build_map(
      data_list  = object$data_list,
      params     = object$estimated_params,
      debug      = FALSE,
      random_rec = object$data_list$random_rec,
      random_sel = isTRUE(as.logical(object$data_list$random_sel)))))
    .shared <- intersect(names(object$map$mapList), names(.fresh$mapList))
    .shared[!vapply(.shared, function(nm) identical(
      as.integer(object$map$mapList[[nm]]),
      as.integer(.fresh$mapList[[nm]])), logical(1))]
  }, error = function(e) NULL)
  if (length(.map_drift)) {
    warning("The stored map differs from `build_map()` in: ",
            paste0("`", .map_drift, "`", collapse = ", "),
            ". Peels use the stored map. If you did not edit or supply the map ",
            "yourself, this fit predates a change to how maps are built; refit it ",
            "to use the current map.", call. = FALSE)
  }

  # env_data columns ----
  # How a peel treats each column, warned about here because a warning in a
  # .parallel_lapply() worker is discarded; judged on the deepest peel, which
  # retains the fewest years.
  .envd <- object$data_list$env_data
  .env_all <- setdiff(names(.envd), "Year")
  .is_num <- vapply(.env_all, function(v) is.numeric(.envd[[v]]), logical(1))
  # The model reads a logical column as 0/1, so a peel withholds it too.
  .env_read <- .env_all[.is_num |
    vapply(.env_all, function(v) is.logical(.envd[[v]]), logical(1))]
  .keep_fixed <- intersect(.env_all, .rce_linkage_fixed_covariates(object$data_list))
  .groups <- intersect(.env_all, .rce_linkage_fixed_covariates(
    object$data_list, re_groups = TRUE, fixed = FALSE))
  .keep_group <- setdiff(.groups, .keep_fixed)
  # Which columns the linkages read as a category, as a number, or as observed.
  .use <- .rce_linkage_column_use(object$data_list)
  .post_any <- any(.envd[["Year"]] > endyr - peels, na.rm = TRUE)
  .ret_deep <- !is.na(.envd[["Year"]]) & .envd[["Year"]] >= styr &
    .envd[["Year"]] <= endyr - peels
  .n_kept <- vapply(.env_read, function(v) {
    x <- .envd[[v]][.ret_deep]
    length(unique(x[is.finite(x)]))
  }, integer(1))

  # A grouping column is a level label and a logical a flag; neither has units.
  # A category would change its levels, cut points or which years fall in each.
  .can_scale <- setdiff(.env_all[.is_num], c(.groups, .use$category))
  if (is.character(rescale)) {
    .bad <- setdiff(rescale, names(.envd))
    if (length(.bad)) {
      stop("`rescale` names column(s) that are not in `env_data`: ",
           paste(.bad, collapse = ", "), ".", call. = FALSE)
    }
    .bad <- setdiff(rescale, .can_scale)
    if (length(.bad)) {
      stop("`rescale` names column(s) that cannot be standardized: ",
           paste(.bad, collapse = ", "), ". `Year`, a random-effect grouping ",
           "column and a logical or non-numeric column have no units to ",
           "standardize, and standardizing a column a linkage reads as a category ",
           "would change its levels, its cut points or which years fall in each.",
           call. = FALSE)
    }
    .rescale_cols <- intersect(.env_all, rescale)
  } else if ((is.logical(rescale) || is.numeric(rescale)) && length(rescale) == 1L &&
             !is.na(rescale)) {
    # rescale also takes names, so a number reads like a column position.
    if (is.numeric(rescale)) {
      warning("`rescale = ", rescale, "` is read as ", as.logical(rescale),
              ", not as a column position. Pass TRUE, FALSE, or the names of ",
              "the env_data columns to standardize.", call. = FALSE)
    }
    .rescale_cols <- if (rescale) .can_scale else character(0)
  } else {
    stop("`rescale` must be TRUE, FALSE, or a character vector of env_data ",
         "column names.", call. = FALSE)
  }

  # Columns known ahead of each peel's terminal year: kept whole, and left out
  # of the kept-column and mean-fill warnings.
  .known <- projected_env_data_known
  if (is.character(.known)) {
    .bad <- setdiff(.known, names(.envd))
    if (length(.bad)) {
      stop("`projected_env_data_known` names column(s) that are not in `env_data`: ",
           paste(.bad, collapse = ", "), ".", call. = FALSE)
    }
  } else if (is.logical(.known) && length(.known) == 1L && !is.na(.known)) {
    .known <- if (.known) .env_all else character(0)
  } else {
    stop("`projected_env_data_known` must be TRUE, FALSE, or a character vector ",
         "of env_data column names.", call. = FALSE)
  }

  # fit_mod() holds a category (a block or regime) from the peel's terminal year.
  # A number takes the retained mean unless it is also a category, observed, a
  # grouping column or known ahead.
  .held <- setdiff(.use$category, c(.use$number, .use$observed, .groups))
  .fill_fixed <- setdiff(intersect(.use$number, .env_all[.is_num]),
                         c(.use$category, .use$observed, .groups, .known))
  .warn_fixed <- setdiff(.keep_fixed, c(.held, .fill_fixed, .known))
  .warn_group <- setdiff(.keep_group, .known)

  .no_kept <- setdiff(.env_read[.n_kept == 0L],
                      c(.keep_fixed, .keep_group, .known))
  if (.post_any && length(.no_kept)) {
    warning("env_data column(s) ", paste(.no_kept, collapse = ", "), " have no ",
            "value in years ", styr, "-", endyr - peels, ", so a peel that retains ",
            "none of them keeps their post-peel values: a blanked column would have ",
            "no retained mean to fill it.", call. = FALSE)
  }
  .raw <- .rescale_cols[.n_kept[.rescale_cols] < 2L]
  if (length(.raw)) {
    warning("`rescale` leaves ", paste(.raw, collapse = ", "), " in their own ",
            "units in any peel that keeps fewer than two distinct values of them.",
            call. = FALSE)
  }
  if (.post_any && length(.fill_fixed)) {
    warning("Fixed-effect linkage covariate(s) ", paste(.fill_fixed, collapse = ", "),
            " cannot be NA, so each peel gives their post-peel years the mean of ",
            "its retained years. List a series known ahead, such as a climate ",
            "projection, in `projected_env_data_known`. A 0/1 block column takes ",
            "the mean too: list it there, or refit it as `factor(<column>)`, which ",
            "holds its terminal-year value under new coefficient names.",
            call. = FALSE)
  }
  if (.post_any && length(.warn_fixed)) {
    warning("Fixed-effect linkage covariate(s) ", paste(.warn_fixed, collapse = ", "),
            " are also observed, used as a random-effect grouping column, read as a ",
            "category, or are logical, so a peel cannot give them the mean and keeps ",
            "their post-peel values. List any known ahead in ",
            "`projected_env_data_known`; otherwise Mohn's rho is conditional on ",
            "those values being known.", call. = FALSE)
  }
  if (.post_any && length(.warn_group)) {
    warning("Random-effect grouping column(s) ", paste(.warn_group, collapse = ", "),
            " label which deviation each year takes, so each peel keeps their ",
            "post-peel values. If a label is derived from data (a warm/cold ",
            "regime), each peel knows which level its forecast years are in.",
            call. = FALSE)
  }

  # Per-peel closure ----
  # Each peel only reads the original model, so peels are independent.
  run_one_peel <- function(i) {

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
    # Blank post-peel values rather than drop rows: a linkage design is built by
    # row position, so a dropped row would drop a random-effect level.
    .yr <- data_list$env_data[["Year"]]
    .post_peel <- !is.na(.yr) & .yr > endyr_peel
    # Retained: the model years the peel fits, which rearrange_data() fills from.
    .retained <- !is.na(.yr) & .yr >= styr & .yr <= endyr_peel
    .blank <- setdiff(.env_read, c(.keep_fixed, .keep_group, .known))
    # A column with no retained value stays whole: blanked, rearrange_data() would
    # fill it with NaN for whatever reads it.
    .blank <- .blank[vapply(.blank, function(v)
      any(is.finite(data_list$env_data[[v]][.retained])), logical(1))]
    if (any(.post_peel) && length(.blank)) {
      data_list$env_data[.post_peel, .blank] <- NA
    }
    # A numeric fixed-effect covariate cannot be NA: it takes the retained mean.
    for (col in .fill_fixed) {
      data_list$env_data[[col]][.post_peel] <- mean(data_list$env_data[[col]][.retained])
    }
    # Centre and scale from the retained years only, applied to the whole column.
    # A column with no retained variation is left alone rather than made NaN.
    for (col in .rescale_cols) {
      kept <- data_list$env_data[[col]][.retained]
      ctr <- mean(kept, na.rm = TRUE)
      scl <- stats::sd(kept, na.rm = TRUE)
      if (!is.finite(scl) || scl == 0) next
      data_list$env_data[[col]] <- (data_list$env_data[[col]] - ctr) / scl
    }

    # * Hindcast ----
    # The model as it would have been fitted at the peel's terminal year: it ends
    # there, so no later index, composition or catch, and no parameter for a later
    # year, enters its estimates. It starts from the unpeeled estimates for its own
    # years and uses the unpeeled map cut to them.
    hind_list <- data_list
    hind_list$endyr <- endyr_peel
    hind_list$catch_data <- dplyr::filter(data_list$catch_data, Year <= endyr_peel)
    hind_defaults <- suppressWarnings(suppressMessages(build_params(hind_list)))
    inits <- .peel_params(object$estimated_params, hind_defaults, nyrs_peel)
    map <- .peel_map(object$map, hind_defaults, nyrs_peel)

    newmod <- suppressWarnings(
      suppressMessages(
        # Clamp the stock-recruit switch and fit-end years and the suitability
        # end year to the peel's terminal year.
        .refit_like(
          data_list        = hind_list,
          inits            = inits,
          map              = map,
          estimateMode     = ifelse(data_list$estimateMode == 2, 0, data_list$estimateMode), # Projection-only (2) is refitted at mode 0
          phase            = phase,
          getsd            = getsd,
          srr_mse_switchyr = min(data_list$srr_mse_switchyr, endyr_peel),
          srr_hat_endyr    = min(data_list$srr_hat_endyr, endyr_peel),
          suit_endyr       = pmin(data_list$suit_endyr, endyr_peel))
      )
    )

    # Judge the hindcast here: the forecast refit below replaces `newmod` and frees only F,
    # so its gradient cannot show whether the hindcast that Mohn's rho reads converged.
    if (!.refit_converged(newmod)) return(NULL)
    hind_sdrep <- newmod$sdrep

    # Forecast ----
    # Back to the full model: the peel's estimates in its own years, and in the
    # removed years F started from the unpeeled estimates and fitted to the
    # observed catch, with q, selectivity and M1 held at the terminal year.
    peeled_pars <- .unpeel_params(newmod$estimated_params, object$estimated_params)
    peeled_pars <- .hold_curves(peeled_pars, data_list, nyrs_peel)

    # - Add in peeled catch to fit to
    data_list$catch_data <- data_list$catch_data |>
      dplyr::filter(Year <= endyr_peel) |>
      rbind(peeled_catch_data) |>
      dplyr::arrange(Fleet_code, Year)

    # The hindcast's map at full length, the removed years fixed: it frees the
    # hindcast again for the report-only pass after the forecast refit.
    map_hindcast <- .unpeel_map(newmod$map, peeled_pars)

    # - Update map
    # Only new parameter we are estimating in is the log_F of the peeled years
    map <- build_map(
      data_list = data_list,
      params = peeled_pars,
      debug = TRUE,
      random_rec = newmod$data_list$random_rec)
    map$mapFactor$dummy <- as.factor(NA); map$mapList$dummy <- NA

    # - Turn on F for peeled years to fit to catch (matches full model)
    map$mapList$log_F[,(nyrs_peel+1):nyrs] <- object$map$mapList$log_F[,(nyrs_peel+1):nyrs]
    map$mapFactor$log_F <-  factor(map$mapList$log_F)

    # Forecast-year rec_dev. Under a penalty-form curve, the years to average
    # R/R_hat over; a peel ending before them uses its own (warned before the
    # peels run, below).
    hat_yrs <- integer(0)
    if (newmod$data_list$srr_fun != newmod$data_list$srr_pred_fun) {
      hat_yrs <- .srr_hat_cols(newmod$data_list, nyrs_peel)
      if (!length(hat_yrs)) hat_yrs <- seq(min(2, nyrs_peel), nyrs_peel)
    }
    for(sp in 1:newmod$data_list$nspp){
      # A species with input numbers-at-age has no rec_dev to set (its R is the input recruits).
      if (isTRUE((newmod$data_list$estDynamics %||% 0)[sp] > 0)) next

      # -- where SR curve is estimated directly
      if(newmod$data_list$srr_fun == newmod$data_list$srr_pred_fun){
        rec_dev <- log(mean(newmod$quantities$R[sp,1:nyrs_peel]))  - log(newmod$quantities$R0[sp])
        # A curve fitted in the hindcast takes its mean deviation from the curve, as sample_rec() does.
        if (newmod$data_list$srr_fun > 1) {
          rec_dev <- log(mean(exp(newmod$estimated_params$rec_dev[sp, 1:nyrs_peel])))
        }
      }

      # -- OMs where SR curve is estimated as penalty (sensu Ianelli)
      if(newmod$data_list$srr_fun != newmod$data_list$srr_pred_fun){
        # Log of the mean of R / R_hat over the penalty years, so projected
        # recruitment is the mean around the curve, not the median.
        rec_dev <- log(mean((newmod$quantities$R / newmod$quantities$R_hat)[sp, hat_yrs]))

      }

      # - Update OM with devs
      peeled_pars$rec_dev[sp, (peel_prj_yrs - styr + 1)] <- replace(
        peeled_pars$rec_dev[sp, (peel_prj_yrs - styr + 1)],
        values =  rec_dev)
    }

    newmod <- suppressMessages(
      suppressWarnings(
        # Second refit: same peeled configuration, now started from the
        # bias-adjusted peeled parameters with peeled-year F turned back on.
        .refit_like(
          data_list        = data_list,
          inits            = peeled_pars,
          map              = map,
          estimateMode     = ifelse(data_list$estimateMode == 2, 0, data_list$estimateMode),
          phase            = phase,
          getsd            = getsd,
          srr_mse_switchyr = min(data_list$srr_mse_switchyr, endyr_peel),
          srr_hat_endyr    = min(data_list$srr_hat_endyr, endyr_peel),
          suit_endyr       = pmin(data_list$suit_endyr, endyr_peel))
      )
    )

    # Judge gradient and bounds: at a bound F can miss the catch with a small gradient.
    # Not the Hessian (F alone). The hindcast is fixed, so terminal-year rho stands.
    if (!.refit_converged(newmod, hessian = FALSE) || length(.check_bounds(newmod))) {
      attr(newmod, "forecast_refit_error") <-
        paste0("peel ", i, " (through ", endyr_peel, "): the forecast refit did ",
               "not fit F to the peeled catch (gradient check failed, or F at a bound)")
    }

    # * Standard errors for the peeled hindcast ----
    # The forecast refit holds the hindcast fixed (SE 0), so rebuild with it free and, to stay
    # at its optimum, the peeled catch dropped and F held; .rce_catch_conditioned_sd() then
    # moves F with the hindcast parameters.
    if (isTRUE(getsd)) {
      sd_try <- tryCatch({
        report_mod <- suppressWarnings(suppressMessages(
          .refit_like(
            data_list        = data_list,
            inits            = newmod$estimated_params,
            map              = map_hindcast,
            estimateMode     = 3,   # build only; do not re-estimate
            getsd            = FALSE,
            srr_mse_switchyr = min(data_list$srr_mse_switchyr, endyr_peel),
            srr_hat_endyr    = min(data_list$srr_hat_endyr, endyr_peel),
            suit_endyr       = pmin(data_list$suit_endyr, endyr_peel))))
        obj <- report_mod$obj
        dat <- obj$env$data
        catch_flt <- dat$catch_ctl[, 1]
        catch_yr <- dat$catch_ctl[, 3]
        peeled <- catch_yr > endyr_peel
        # The rows the forecast refit fitted: fishery catch through endyr whose
        # F that refit estimated (`map` is its map).
        yr_i <- catch_yr - styr + 1
        f_free <- !is.na(map$mapList$log_F[cbind(catch_flt, pmin(pmax(yr_i, 1), nyrs))])
        fitted_rows <- which(peeled & catch_yr <= endyr & dat$flt_type[catch_flt] == 1 &
                               dat$catch_obs[, 1] > 0 & f_free)
        # A catch of 0 is read as missing; F, and so the removals, are unchanged.
        obj$env$data$catch_obs[peeled, 1] <- 0
        obj$retape()   # rebuild the model with the changed catch
        # The hindcast's own Hessian, not this full-length model's, which also
        # carries the removed years' penalties. The correction needs the joint
        # precision under random effects. Quiet: an indefinite Hessian warns "NaNs
        # produced", and the reason is recorded below.
        sdr <- suppressWarnings(TMB::sdreport(
          obj, hessian.fixed = solve(hind_sdrep$cov.fixed),
          getJointPrecision = length(obj$env$random) > 0))
        tryCatch({
          if (isTRUE(sdr$pdHess) && length(fitted_rows)) {
            sdr <- .rce_catch_conditioned_sd(obj, sdr, fitted_rows, data_list$styr)
          }
          list(sdrep = sdr)
        }, error = function(e) list(
          sdrep = sdr,
          error = paste("the catch-conditioned correction failed:",
                        conditionMessage(e))))
      }, error = function(e) list(error = conditionMessage(e)))

      # A rebuild can fail without throwing: TMB::sdreport() returns pdHess = FALSE
      # on a Hessian it cannot factor. Name which condition fired for the reader.
      .sdr <- sd_try[["sdrep"]]
      if (is.null(sd_try[["error"]])) {
        if (isFALSE(.sdr[["pdHess"]])) {
          sd_try[["error"]] <- "the rebuilt hindcast Hessian is not positive definite"
        } else if (!any(is.finite(.sdr[["sd"]]))) {
          sd_try[["error"]] <- "the rebuilt hindcast reported no finite standard error"
        }
      }

      # Keep the point estimates and what TMB returned (none if the rebuild threw
      # first). The reason stops plots and tables reading it as standard errors.
      newmod$sdrep <- .sdr
      if (!is.null(sd_try[["error"]])) {
        attr(newmod, "sd_report_error") <-
          paste0("peel ", i, " (through ", endyr_peel, "): ", sd_try[["error"]])
      }
    }


    # * Report the peel's own terminal year ----
    # Only after the refits: `endyr` sizes the model, and the forecast refit
    # indexes log_F over the full hindcast. Here it is metadata (the arrays
    # still span the full hindcast): plots draw styr:endyr by default, so a peel
    # shows only as far as it was fit. Rho reads endyr_peel.
    newmod$data_list$endyr <- endyr_peel

    newmod
  } # End run_one_peel closure


  # A penalty-form peel with no penalty years has no curve fitted to recruitment.
  # Warned here, before the peels run: warnings from peels run in parallel are lost.
  penalty_form <- isTRUE(suppressWarnings(
    .srr_fun_structural(data_list$srr_fun) != .srr_fun_structural(data_list$srr_pred_fun)))
  if (penalty_form) {
    hat_s <- data_list$srr_hat_styr %||% (styr + 1)
    hat_e <- data_list$srr_hat_endyr %||% endyr
    no_curve <- vapply(seq_len(peels), function(i) {
      e <- endyr - i
      !length(.srr_hat_cols(list(styr = styr, srr_hat_styr = hat_s,
                                 srr_hat_endyr = min(hat_e, e)), e - styr + 1))
    }, logical(1))
    if (any(no_curve)) {
      warning("Peel(s) ending ", paste(endyr - which(no_curve), collapse = ", "),
              " contain no stock-recruit penalty years, so their curve stays at the ",
              "unpeeled fit's estimates (fitted to years the peel withholds) unless a ",
              "prior or a reference-point penalty changes it; their projected recruitment deviation averages over the ",
              "peel's own years after the first (or over its one year, for a peel that keeps a single year). ",
              "With getsd = TRUE such a peel is usually dropped for a ",
              "non-positive-definite Hessian.", call. = FALSE)
    }
  }

  # Dispatch peels (parallel via PSOCK or sequential) ----
  if (use_parallel) {
    peel_results <- .parallel_lapply(1:peels, run_one_peel, min(cores, peels), environment())
  } else {
    peel_results <- lapply(1:peels, run_one_peel)
  }

  # Drop peels whose hindcast did not converge and prepend the original model
  peel_results <- peel_results[!vapply(peel_results, is.null, logical(1))]
  # A warning, not a message: rho below is averaged over the peels that
  # survive, so a drop changes the number that gets reported.
  .report_dropped(peels - length(peel_results), peels, "peel", warn = TRUE)

  mod_list <- c(list(object), peel_results)

  # With no peel left every rho below is 0/0; say so rather than return a table
  # of NaN that looks computed.
  if (length(peel_results) == 0L) {
    warning("No peel converged, so Mohn's rho is undefined (every value NaN). ",
            "Inspect a peel with retrospective(..., peels = 1) and read its ",
            "$convergence.", call. = FALSE)
  }


  # Calculate Mohn's rho ----
  # * Data frame to save ----
  objects <- c("biomass", "ssb", "R", "F_spp")
  mohns <- data.frame(matrix(0, nrow = length(objects) * (nyrs_forecast+1), ncol = 3 + data_list$nspp))
  colnames(mohns) <- c("Object", "Forecast year", "N", data_list$spnames)

  # * Loop through peels ----
  # seq_len, not 1:(n-1): with every peel dropped `mod_list` holds only the
  # original model and `1:0` would count DOWN to mod_list[[2]].
  for (i in seq_len(length(mod_list) - 1)) {
    endyr_peel <- mod_list[[i + 1]]$data_list$endyr_peel
    nyrs_peel <- mod_list[[i + 1]]$data_list$endyr_peel - styr + 1
    ind <- 1
    # F not fitted to the peeled catch: its forecast years are not a forecast of
    # that catch, while its terminal year comes from the hindcast.
    forecast_fitted <- is.null(
      attr(mod_list[[i + 1]], "forecast_refit_error", exact = TRUE))

    # * Loop output ----
    for (j in 1:length(objects)) {
      term_quantities <- mod_list[[1]]$quantities[[objects[j]]]
      retro_quantities <- mod_list[[i + 1]]$quantities[[objects[j]]]

      # * Loop forecast years ----
      for(yr in 0:nyrs_forecast){

        # If data exist for forecast (save)
        if(endyr_peel + yr <= endyr){

          mohns[ind, 1] <- objects[j]         # Object
          mohns[ind, 2] <- yr                 # Year

          if (yr == 0 || forecast_fitted) {
            # * Get full and peeled models ----
            base <- term_quantities[, nyrs_peel + yr]
            peel <- retro_quantities[, nyrs_peel + yr]
            rel_error <- ((peel - base)/base)

            # * Save and sum relative error ----
            mohns[ind, 3] <- mohns[ind, 3] + 1  # N
            mohns[ind, 4:(data_list$nspp + 3) ] <- mohns[ind, 4:(data_list$nspp + 3)] + rel_error # Relative error
          }
        }
        ind = ind+1
      }
    }
  }

  # * Divide N ----
  mohns[, 4:(data_list$nspp + 3) ] <- mohns[, 4:(data_list$nspp + 3)]/mohns[, 3]

  mod_list <- rev(mod_list)
  # Named from each model's own terminal year, so a dropped peel leaves a gap
  # rather than shifting the labels; the unpeeled model's endyr_peel is endyr.
  names(mod_list) <- paste0(
    "Year_",
    vapply(mod_list, function(x) as.numeric(x$data_list$endyr_peel), numeric(1)))

  # `peels_requested` lets print.Rceattle_retro() show a dropped peel after the
  # object is saved and read back, when the drop warning is gone.
  structure(list(Rceattle_list = mod_list, mohns = rbind(mohns),
                 peels_requested = peels),
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
#' @param x A `"Rceattle_retro"` object from [retrospective()].
#' @param band Symmetric reference band for Mohn's rho. Default `0.2`.
#' @param ... Currently unused.
#' @return `x`, invisibly.
#' @references Hurtado-Ferro, F., et al. 2015. Looking in the rear-view mirror:
#'   bias and retrospective patterns in integrated, age-structured stock
#'   assessment models. ICES J. Mar. Sci. 72:99-110.
#' @export
print.Rceattle_retro <- function(x, band = 0.2, ...) {
  m <- as.data.frame(x$mohns)
  spp <- setdiff(names(m), c("Object", "Forecast year", "N"))
  term <- m[!is.na(m[["Forecast year"]]) & m[["Forecast year"]] == 0, , drop = FALSE]

  rho <- suppressWarnings(as.numeric(unlist(term[, spp, drop = FALSE])))
  sev <- ifelse(is.na(rho), "OK", ifelse(abs(rho) > band, "WARN", "OK"))
  n_bad <- sum(sev == "WARN")

  # Rho averages the surviving peels, so a drop changes it. `peels_requested` is
  # absent on an object saved before it was carried.
  n_kept <- length(x$Rceattle_list) - 1L          # the unpeeled model
  n_asked <- x$peels_requested
  n_dropped <- if (is.null(n_asked)) 0L else max(0L, n_asked - n_kept)
  # A peel out of the forecast rows also changes reported numbers, so it is as visible.
  n_noforecast <- sum(vapply(
    x$Rceattle_list,
    function(m) !is.null(attr(m, "forecast_refit_error", exact = TRUE)),
    logical(1)))
  if (n_dropped > 0L || n_noforecast > 0L) sev <- c(sev, "NOTE")

  .rce_diag_header(
    "retrospective", .rce_worst(sev),
    paste0(if (n_dropped > 0L) paste0(n_kept, " of ", n_asked, " peel(s); ")
           else paste0(length(x$Rceattle_list), " peel(s); "),
           .rce_n_of(n_bad, length(rho)),
           " terminal Mohn's rho outside +/-", band))

  if (n_dropped > 0L) {
    cat("  ", n_dropped, " peel(s) dropped as non-converged, so rho is ",
        if (n_kept > 0L) paste("averaged over", n_kept) else "undefined",
        "\n", sep = "")
  }

  # A kept peel can still have no band and the call does not warn, so say it
  # here -- rho is unaffected either way.
  n_noband <- sum(vapply(
    x$Rceattle_list,
    function(m) !is.null(attr(m, "sd_report_error", exact = TRUE)),
    logical(1)))
  if (n_noband > 0L) {
    cat("  ", n_noband, " peel(s) have no usable standard errors",
        " (see attr(peel, \"sd_report_error\"))\n", sep = "")
  }

  if (n_noforecast > 0L) {
    cat("  ", n_noforecast, " peel(s) left out of the forecast-year rows of $mohns",
        " (see attr(peel, \"forecast_refit_error\"))\n", sep = "")
  }

  if (nrow(term)) {
    show <- term
    show$.tag <- vapply(seq_len(nrow(show)), function(i) {
      r <- suppressWarnings(as.numeric(show[i, spp]))
      .rce_sev_tag(if (any(!is.na(r) & abs(r) > band)) "WARN" else "OK")
    }, character(1))
    .rce_diag_table(show, stats::setNames(
      c(".tag", "Object", "N", spp), c(" ", "quantity", "N", spp)))
  }
  if (any(!is.na(m[["Forecast year"]]) & m[["Forecast year"]] > 0)) {
    cat("  forecast-skill peels are in $mohns; the +/-", band,
        "rule is for the terminal peel only\n")
  }
  invisible(x)
}


# A peel's hindcast ends at its terminal year, so its hindcast-year parameters
# (log_F and the catchability and selectivity deviations) are shorter. Year is
# their last dimension, and R fills an array column by column, so a shorter
# model's years are the leading cells of the longer model's array. rec_dev and
# log_M1_dev run to projyr in every model, so only their later years change role.

#' The unpeeled estimates as starting values for a peel's hindcast: its own years
#' from the unpeeled fit, later years of rec_dev and log_M1_dev from `defaults`
#' (the peel's own starting values).
#' @noRd
.peel_params <- function(params, defaults, nyrs_peel) {
  yrs <- seq_len(nyrs_peel)
  for (nm in intersect(names(defaults), names(params))) {
    if (nm == "rec_dev") {
      defaults[[nm]][, yrs] <- params[[nm]][, yrs]
    } else if (nm == "log_M1_dev") {
      defaults[[nm]][, , , yrs] <- params[[nm]][, , , yrs]
    } else {
      defaults[[nm]][] <- params[[nm]][seq_along(defaults[[nm]])]
    }
  }
  defaults
}

#' The unpeeled map cut to a peel's years, so a parameter fixed, freed or shared
#' by hand stays so. A later year of rec_dev or log_M1_dev is fixed unless it
#' shares its level with one of the peel's years, as an M1 deviation by age does.
#' @noRd
.peel_map <- function(map, defaults, nyrs_peel) {
  for (nm in intersect(names(defaults), names(map$mapList))) {
    x <- map$mapList[[nm]]
    if (nm %in% c("rec_dev", "log_M1_dev")) {
      later <- slice.index(x, length(dim(x))) > nyrs_peel
      x[later & !x %in% x[!later]] <- NA
    } else if (length(x) > length(defaults[[nm]])) {
      cut <- defaults[[nm]]
      cut[] <- x[seq_along(cut)]
      x <- cut
    } else {
      next
    }
    map$mapList[[nm]] <- x
    map$mapFactor[[nm]] <- factor(x)
  }
  map
}

#' The peel's estimates back in arrays as long as the unpeeled model's: a shortened
#' array takes the unpeeled values in the years the peel removed; rec_dev and
#' log_M1_dev, already full length, are the peel's.
#' @noRd
.unpeel_params <- function(peel, full) {
  for (nm in intersect(names(full), names(peel))) {
    full[[nm]][seq_along(peel[[nm]])] <- peel[[nm]]
  }
  full
}

#' The peel's map at full length, with the years it removed fixed.
#' @noRd
.unpeel_map <- function(map, full) {
  for (nm in names(map$mapList)) {
    x <- map$mapList[[nm]]
    if (length(x) < length(full[[nm]])) {
      y <- full[[nm]]
      y[] <- NA            # the full array's shape, every cell fixed
      y[seq_along(x)] <- x
      map$mapList[[nm]] <- y
      map$mapFactor[[nm]] <- factor(y)
    }
  }
  map
}


#' Jitter analysis
#'
#' @description Refits the Rceattle model from starting values perturbed by N(0, sd) around the model's initial (pre-fit) parameters, to check convergence robustness.
#'
#' @details Each start refits under the fit's own map, so a parameter fixed there keeps
#' its value and is not estimated. The objectives do not compare with the fit's if it
#' used \code{fit_mod(bounds = )} (a refit rebuilds the bounds) or is projection-only
#' (\code{estimateMode = "Projection"}). A projection-only fit refits under a rebuilt map,
#' and its starts are not perturbed unless it was fitted under a supplied map.
#'
#' @note Attaching Rceattle masks \code{\link[base]{jitter}}, because that
#'   function is not generic and this name is exported, so the mask cannot be
#'   removed by making this an S3 method. It is harmless for the spellings base
#'   callers use: a numeric first argument is handed straight to
#'   \code{base::jitter()}, so \code{jitter(1:5)} and
#'   \code{jitter(1:5, factor = 2)} behave as base does. \code{jitter(x = 1:5)}
#'   and a positional second argument do not. \code{x} is not a parameter here,
#'   so it is refused as an unused argument; a positional second argument binds
#'   \code{njitter}, so \code{jitter(1:5, 2)} jitters with base's defaults and
#'   silently drops the \code{2}. Name base's arguments, or call
#'   \code{base::jitter()} directly.
#'
#' @param ... for a numeric \code{object}, passed to \code{base::jitter()};
#'   for a fitted model, refused, so that a misspelled argument is not silently
#'   ignored.
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
#'   model's setting (\code{TRUE} only if it holds an \code{sdrep}).
#' @param timeout elapsed-second limit per jitter, \code{Inf} (default) for none.
#'   A jitter is a deliberately perturbed start and the optimizer runs with no
#'   iteration cap, so this is the diagnostic most likely to send one somewhere
#'   pathological and stall the whole run, a hang no convergence check can
#'   catch, because the fit never returns. One that exceeds the limit is stopped,
#'   counted as non-converged and reported separately. Approximate: the limit is
#'   checked when control returns to R, so it fires between the optimizer's
#'   function evaluations rather than inside one.
#'
#' @return a list of 1. \code{Rceattle_list}, the converged jitters, and
#'   2. \code{nll}, their objective values. Non-converged (or timed-out) starts
#'   are dropped and reported in a message, so both can be shorter than
#'   \code{njitter}, and that count is itself the result, since the whole point
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
jitter <- function(object = NULL, njitter = 50, sd = 0.2, phase = FALSE, seed = 123, cores = NULL, getsd = NULL, timeout = Inf, fit_control = NULL, ..., Rceattle = NULL) {
  # A numeric first argument is base::jitter()'s, not a fitted model. This name
  # masks base::jitter(), which is not generic and so cannot be dispatched on,
  # and S3 is no escape: UseMethod() errors before any method runs when the
  # dispatch argument is absent, so a generic on `x` could not also accept the
  # model under `object` or `Rceattle`. Forwarding here is what makes the mask
  # harmless.
  if (is.numeric(object)) return(base::jitter(object, ...))

  # Past the numeric path nothing reads `...`, so an entry there is a typo:
  # jitter(mod, njiter = 10) would otherwise silently run the default 50.
  if (...length()) {
    # Nine formals precede `...`, so an unnamed extra takes ten positional
    # arguments to reach; a mixed call gives "" beside a name.
    nm <- ...names()
    if (is.null(nm)) nm <- rep("", ...length())
    nm[!nzchar(nm)] <- "<unnamed>"
    stop("jitter() does not take these arguments: ", paste(nm, collapse = ", "),
         ". A fitted model takes njitter, sd, phase, seed, cores, getsd, ",
         "timeout and fit_control; for base::jitter()'s own arguments give a ",
         "numeric first argument, or call base::jitter() directly.",
         call. = FALSE)
  }

  # `Rceattle` was the old name for `object`; see R/0-deprecate.R.
  if (!missing(Rceattle))
    object <- .rce_deprecated_arg(Rceattle, !missing(object), "Rceattle", "object", "jitter")
  # Cleared once consumed: .parallel_lapply() exports every binding in this frame
  # to the PSOCK workers, and two names for one fitted model send it twice.
  Rceattle <- NULL

  if (!inherits(object, "Rceattle")) {
    stop("Object is not of class 'Rceattle'")
  }

  # Jitters inherit the input model's sdreport setting unless overridden;
  # multimodality is judged from objectives and point estimates, not sdrep.
  if (is.null(getsd)) getsd <- !is.null(object$sdrep)

  # fit_control() overrides `phase` and `getsd`, but only where the caller
  # changed them from fit_control()'s own defaults; see .rce_refit_control().
  ctl <- .rce_refit_control(fit_control, "jitter")
  if (!is.null(ctl$phase)) phase <- ctl$phase
  if (!is.null(ctl$getsd)) getsd <- ctl$getsd

  # Cross-platform parallel: a FORK cluster where available, PSOCK on
  # Windows (see .parallel_lapply())
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

  # Per-jitter closure ----
  run_one_jitter <- function(i) {

    set.seed(seed + i) # unique seed per jitter for reproducibility under parallel

    # * Adjust initial values ----
    inits <- object$initial_params
    mapFactor <- object$map$mapFactor
    data_list <- object$data_list

    # Perturb only what the map estimates, read from mapFactor because that is what TMB
    # fixes; a fixed cell keeps its starting value.
    for(j in 1:length(inits)){
      par <- names(inits)[j]
      inits[[j]] <- replace(inits[[j]],
                            values = ifelse(is.na(as.numeric(mapFactor[[par]])),
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
          # Refit from the jittered starting values, under the fit's own map.
          .refit_like(
            data_list        = data_list,
            inits            = inits,
            map              = .refit_map(object),
            # As in retrospective(): only a projection-only fit (2) refits at 0.
            estimateMode     = ifelse(data_list$estimateMode == 2, 0, data_list$estimateMode),
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


  # Dispatch jitters (parallel via PSOCK or sequential) ----
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
  # `njitter` is carried so the print method can report the fraction of STARTS
  # that reached the optimum: the converged runs alone cannot say how many were
  # attempted, and that fraction is the result jitter() exists to produce.
  structure(list(Rceattle_list = mod_list, nll = unname(jnll),
                 njitter = njitter),
            class = "Rceattle_jitter")
}


#' Print method for a jitter analysis
#'
#' @description Reports what the run was for: how many random starts reached the
#' best optimum found. The objective values alone cannot say that, non-converged
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


#' Catch-conditioned standard errors for a retrospective peel
#'
#' The rebuilt peel holds the peeled (forecast) years' F at the values fitted
#' to their catch. A catch-conditioned band instead lets that F move with the
#' hindcast parameters so the predicted catch stays at the observed catch, and
#' the delta method carries that into each standard error. Catch is treated as
#' known.
#'
#' With theta the hindcast's parameters (fixed and random effects), `j_q` holds
#' each reported quantity's derivative with respect to theta and the peeled F,
#' and `j_c` each peeled log predicted catch's. Holding catch fixed gives
#' `j_cc = j_q[theta] - j_q[F] j_c[F]^-1 j_c[theta]`, each quantity's
#' sensitivity to theta once F adjusts to keep the catch, and the covariance is
#' `j_cc cov_theta t(j_cc)`, with cov_theta the covariance of theta. For example, a
#' hindcast parameter that implies a smaller stock needs a higher F to take the same
#' catch, which lowers forecast-year SSB further; the band carries that.
#'
#' Both Jacobians come from one reverse-mode pass per reported quantity, which
#' is why predicted catch joins the ADREPORT set here (`adreport_catch`) and not
#' in a fit, where it would add to `sdrep$value`.
#'
#' @param obj the rebuilt peel's TMB object, with the peeled F held fixed.
#' @param sdrep its `TMB::sdreport()`, with the joint precision (inverse
#'   covariance of fixed and random effects) when `obj` has random effects.
#' @param rows rows of the catch table (`catch_obs`) for the peeled years that
#'   had a catch to fit.
#' @param styr first model year.
#' @return `sdrep` with `sd` and `cov` replaced and `jointPrecision` dropped;
#'   `par.fixed`, `cov.fixed` and `gradient.fixed` still describe the model
#'   with the peeled F held fixed.
#' @noRd
.rce_catch_conditioned_sd <- function(obj, sdrep, rows, styr) {
  dat <- obj$env$data
  p <- obj$env$last.par.best   # the peel's estimates, fixed and random together
  # The log_F cell (fleet, model-year index) that each peeled catch row determines
  flt <- dat$catch_ctl[rows, 1]
  yr <- dat$catch_ctl[rows, 3] - styr + 1
  cells <- unique(cbind(flt, yr))
  n_f <- nrow(cells)

  # The same model with those F freed, numbered after the estimated F so they come
  # last among log_F. ADreport = TRUE: ad$fn() returns the quantities sdreport gives
  # errors for (SSB, R, ...) plus predicted catch, and ad$gr() their exact derivatives.
  map <- obj$env$map
  # TMB keeps the full log_F array, held values included, as "shape".
  lf_shape <- attr(obj$env$parameters$log_F, "shape")
  lf <- array(as.integer(map$log_F), dim = dim(lf_shape))
  lf[cells] <- max(c(0L, lf), na.rm = TRUE) + seq_len(n_f)
  map$log_F <- factor(lf)
  dat$adreport_catch <- 1   # double: MakeADFun() has already checked this data
  ad <- TMB::MakeADFun(dat, obj$env$parList(par = p), map = map, ADreport = TRUE,
                       DLL = obj$env$DLL, silent = TRUE)

  # ad starts at the peel's estimates, the freed F at their held values.
  pa <- ad$par
  f_col <- utils::tail(which(names(pa) == "log_F"), n_f)
  value <- ad$fn(pa)
  is_catch <- names(value) == "catch_hat"

  # Dividing by predicted catch turns its derivative into that of its log (j_c).
  jac <- ad$gr(pa)
  j_q <- jac[!is_catch, , drop = FALSE]
  catch_at <- which(is_catch)[rows]
  j_c <- jac[catch_at, , drop = FALSE] / value[catch_at]

  # qr.solve() applies j_c[F]^-1 from the formula above (least squares when a fleet
  # has several catch rows in a year).
  keep <- setdiff(seq_along(pa), f_col)
  j_cc <- j_q[, keep, drop = FALSE] - j_q[, f_col, drop = FALSE] %*%
    qr.solve(j_c[, f_col, drop = FALSE], j_c[, keep, drop = FALSE])
  # Covariance of theta: fixed and random effects jointly, else the fixed.
  cov_theta <- if (length(obj$env$random)) solve(sdrep$jointPrecision) else sdrep$cov.fixed
  cov_q <- j_cc %*% cov_theta %*% t(j_cc)   # the delta method

  sdrep$sd <- sqrt(diag(cov_q))
  sdrep$cov <- cov_q
  # The joint precision describes the model with F held fixed, not this band.
  sdrep$jointPrecision <- NULL
  sdrep
}
