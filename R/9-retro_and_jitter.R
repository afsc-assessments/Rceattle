#' Retrospective peels
#'
#' @description Calculate Mohn's rho and run retrospective peels for an Rceattle model. The function also evaluates retrospective forecast skill. To evaluate both retrospective bias and forecast skill, the function uses the map functionality of TMB to peel the model:
#' 1. Filters data, filters fixed inputs, and maps out time-varying parameters for the peeled years. All time-varying parameters for the peeled years are set to the terminal year of the model for that peel.
#' 2. Fits the peeled model.
#' 3. Turns off all hindcast parameters, turns on F for the peeled years, and fits to the peeled catch series to update the "forecast" dynamics given projection assumptions and observed catch from the peeled years.
#'
#' @inheritParams rceattle-refit-args
#' @param peels the number of retrospective peels to use in the calculation of rho and for model estimation
#' @param rescale TRUE/FALSE whether to subset and rescale environmental predictors for the range of peel years.
#' @param nyrs_forecast Number of forecast years to calculate Mohn's Rho in addition to terminal year
#' @param getsd whether each peel runs \code{TMB::sdreport} (standard errors).
#'   Costs an extra model build per peel; see Details.
#'   Mohn's rho uses only point estimates, so \code{FALSE} is faster with no
#'   effect on rho. Default \code{NULL} inherits the input model's setting
#'   (\code{TRUE} if it was fit with \code{getsd = TRUE}, i.e. holds an
#'   \code{sdrep}); the returned peel models then hold standard errors only
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
#' @return a list of 1. list of Rceattle models and 2. vector of Mohn's rho for
#'   each species.
#'
#'   A peel that did not converge is dropped, so \code{Rceattle_list} can be
#'   shorter than \code{peels + 1} (a message reports how many). Each entry is
#'   named for its own terminal year (\code{Year_2017}, ...) rather than by
#'   position, so index it by name, \code{Rceattle_list[[3]]} is not
#'   necessarily the 3-year peel. With no peel left, Mohn's rho is \code{NaN}
#'   and the function warns.
#'
#'   Each peel reports its own terminal year as \code{data_list$endyr}, so plots
#'   draw it only as far as it was fit and the peels fan out.
#'
#'   A peel still estimates the years it dropped. They are its retrospective
#'   forecast, fit to the observed catch with the survey and composition data
#'   withheld. Their recruitment deviation is the one [sample_rec()] sets with
#'   `sample_rec = FALSE`, computed from the peel's own fit; a penalty-form peel
#'   with no penalty years averages over its own years after the first (or over
#'   its one year, for a peel that keeps a single year), with a warning. Three
#'   years therefore matter, and each peel has all three:
#'   \describe{
#'     \item{\code{endyr}, \code{endyr_peel}}{the peel's terminal year, what it
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
#'   Catchability is estimated only for a fleet that holds fitted index rows
#'   (see \code{\link{build_map}}), and a peel moves \code{endyr}. A survey whose
#'   index observations all fall in the peeled-off years therefore has no q
#'   estimated in that peel, the parameter count is not constant across peels.
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
retrospective <- function(object = NULL, peels = 5, rescale = FALSE, nyrs_forecast = 3, cores = NULL, getsd = NULL, phase = TRUE, fit_control = NULL, Rceattle = NULL) {
  # `phase` and `fit_control` are appended after the arguments that predate them,
  # so a positional call keeps its meaning. The deprecated `Rceattle` formal sits
  # last; see R/0-deprecate.R.
  if (!missing(Rceattle))
    object <- .rce_deprecated_arg(Rceattle, !missing(object), "Rceattle", "object", "retrospective")
  # Cleared once consumed: .parallel_lapply() exports every binding in this frame
  # to the PSOCK workers, and two names for one fitted model send it twice.
  Rceattle <- NULL

  if (!inherits(object, "Rceattle")) {
    stop("Object is not of class 'Rceattle'")
  }

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
  use_parallel <- peels > 1L && cores > 1L

  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Per-peel closure ----
  # Each peel only reads the original model, so peels are independent.
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
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

    # * Rescale environmental predictors ----
    if(rescale){
      data_list$env_data <- data_list$env_data |>
        dplyr::filter(Year <= endyr_peel)
      data_list$env_data[,2:ncol(data_list$env_data)]<-scale(data_list$env_data[,2:ncol(data_list$env_data)])
    }

    # * Adjust parameters ----
    #FIXME: adjust for forecasting via MVN
    inits <- object$estimated_params
    inits$rec_dev[, (nyrs_peel + 1):nyrs_proj] <- 0
    inits$log_M1_dev[,,,(nyrs_peel+1):nyrs_proj] <- inits$log_M1_dev[,,,nyrs_peel]
    inits$index_q_dev[,(nyrs_peel+1):nyrs] <- inits$index_q_dev[,nyrs_peel]
    inits$log_sel_slp_dev[,,,(nyrs_peel+1):nyrs] <- inits$log_sel_slp_dev[,,,nyrs_peel]
    inits$sel_inf_dev[,,,(nyrs_peel+1):nyrs] <- inits$sel_inf_dev[,,,nyrs_peel]
    inits$sel_coff_dev[,,,(nyrs_peel+1):nyrs] <- inits$sel_coff_dev[,,,nyrs_peel]

    # * Adjust map size ----
    # Turn off forecasted parameters
    map <- object$map
    map$mapList$rec_dev[, (nyrs_peel + 1):nyrs_proj] <- NA
    map$mapFactor$rec_dev <- factor(map$mapList$rec_dev)

    map$mapList$log_M1_dev[,,,(nyrs_peel+1):nyrs_proj] <- NA
    map$mapFactor$log_M1_dev <- factor(map$mapList$log_M1_dev)

    map$mapList$index_q_dev[,(nyrs_peel+1):nyrs] <- NA
    map$mapFactor$index_q_dev <- factor(map$mapList$index_q_dev)

    map$mapList$log_sel_slp_dev[,,,(nyrs_peel+1):nyrs] <- NA
    map$mapFactor$log_sel_slp_dev <- factor(map$mapList$log_sel_slp_dev)

    map$mapList$sel_inf_dev[,,,(nyrs_peel+1):nyrs] <- NA
    map$mapFactor$sel_inf_dev <- factor(map$mapList$sel_inf_dev)

    map$mapList$sel_coff_dev[,,,(nyrs_peel+1):nyrs] <- NA
    map$mapFactor$sel_coff_dev <- factor(map$mapList$sel_coff_dev)

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

    # - Turn on F for peeled years to fit to catch (matches full model)
    peeled_pars$log_F[,(nyrs_peel+1):nyrs] <- object$estimated_params$log_F[,(nyrs_peel+1):nyrs]
    map$mapList$log_F[,(nyrs_peel+1):nyrs] <- object$map$mapList$log_F[,(nyrs_peel+1):nyrs]
    map$mapFactor$log_F <-  factor(map$mapList$log_F)

    # Adjust forecased rec_dev in new mod for bias and refit
    # Penalty years to average R/R_hat over. A peel that ends before them uses its
    # own years instead; retrospective() warns about this before the peels run.
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
        # A peel that cannot be reported is not a peel that failed: keep the
        # fit and its point estimates, and say why the band is missing.
        warning("Peel ", i, ": could not report hindcast standard errors (",
                conditionMessage(e), "). Point estimates are unaffected.",
                call. = FALSE)
        newmod$sdrep
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
      return(newmod)
    }
    return(NULL)
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

  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Dispatch peels (parallel via PSOCK or sequential) ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  if (use_parallel) {
    peel_results <- .parallel_lapply(1:peels, run_one_peel, min(cores, peels), environment())
  } else {
    peel_results <- lapply(1:peels, run_one_peel)
  }

  # Drop non-converged peels and prepend the original model
  peel_results <- peel_results[!vapply(peel_results, is.null, logical(1))]
  # A warning, not a message: rho below is averaged over the peels that
  # survive, so a drop changes the number that gets reported.
  .report_dropped(peels - length(peel_results), peels, "peel", warn = TRUE)
  mod_list <- c(list(object), peel_results)

  # Mohn's rho averages the peels against the full model, so with none left the
  # sums below stay at their initialized zeros and every species column comes
  # back 0/0. Say so rather than return a table of NaN that looks computed.
  if (length(peel_results) == 0L) {
    warning("No peel converged, so Mohn's rho is undefined (every value NaN). ",
            "Inspect a peel with retrospective(..., peels = 1) and read its ",
            "$convergence.", call. = FALSE)
  }


  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # Calculate Mohs rho ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
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
          mohns[ind, 1] <- objects[j]         # Object
          mohns[ind, 2] <- yr                 # Year
          mohns[ind, 3] <- mohns[ind, 3] + 1  # N
          mohns[ind, 4:(data_list$nspp + 3) ] <- mohns[ind, 4:(data_list$nspp + 3)] + rel_error # Relative error
        }
        ind = ind+1
      }
    }
  }

  # * Divide N ----
  mohns[, 4:(data_list$nspp + 3) ] <- mohns[, 4:(data_list$nspp + 3)]/mohns[, 3]



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
  # the list length alone cannot show a drop.
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

  # Rho is averaged over the peels that survived, so a drop changes it. The
  # list length alone cannot show one; `peels_requested` is absent on an object
  # saved before it was carried.
  n_kept <- length(x$Rceattle_list) - 1L          # the unpeeled model
  n_asked <- x$peels_requested
  n_dropped <- if (is.null(n_asked)) 0L else max(0L, n_asked - n_kept)
  if (n_dropped > 0L) sev <- c(sev, "NOTE")

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
