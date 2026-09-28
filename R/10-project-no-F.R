#' Rerun with F = 0.
#'
#' @description
#' Refits the model with fishing mortality set to 0 from `styr` on, keeping
#' every other parameter. The projection after `endyr` is always unfished.
#' `run_mse()` uses it for the no-fishing run (`OM_no_F`) behind the collapse
#' metrics.
#'
#' @details
#' `styr` may be any year from the model's own `styr` to `endyr + 1`; the
#' projection is unfished whatever harvest control rule the model was fit under.
#' Under predation, empirical suitability (`suitMode = 0`) is derived from the
#' fitted abundance over each predator's `suit_styr:suit_endyr`, so `styr` must
#' fall after that window for every predator with `suitMode = 0` and diet data
#' in it.
#'
#' @param object A fitted Rceattle model object
#' @param styr First year with F = 0; default `endyr + 1`, which leaves the hindcast unchanged.
#' @param Rceattle deprecated name for `object`, still accepted so existing
#'   scripts keep working. Supplying both is an error.
#' @export
#'
remove_F <- function(object = NULL, styr = NULL, Rceattle = NULL){
  # `Rceattle` was the old name for `object`; see R/0-deprecate.R.
  if (!missing(Rceattle))
    object <- .rce_deprecated_arg(Rceattle, !missing(object), "Rceattle", "object", "remove_F")

  if (!inherits(object, "Rceattle")) {
    stop("`object` must be a fitted Rceattle model (from fit_mod()).", call. = FALSE)
  }

  dl <- object$data_list
  if (is.null(styr)) styr <- dl$endyr + 1
  # The projection is always unfished, so the no-F period starts by endyr + 1.
  if (!is.numeric(styr) || length(styr) != 1 || is.na(styr) ||
      styr != round(styr) || styr < dl$styr || styr > dl$endyr + 1) {
    stop("`styr` must be a single year from the first model year (", dl$styr,
         ") to the year after endyr (", dl$endyr + 1, "); the projection is always unfished.",
         call. = FALSE)
  }

  # Empirical suitability (suitMode 0) is derived from the fitted abundance over each
  # predator's suit_styr:suit_endyr, so removing F inside that window re-derives it. A
  # predator with no diet data in its window has an all-zero suitability slice; a
  # parametric form does not read abundance.
  if (isTRUE(.map_switch(dl$msmMode, msmMode_map, "msmMode") > 0)) {
    suit_mode <- rep_len(.map_switch(dl$suitMode, suitMode_map, "suitMode"), dl$nspp)
    suit <- object$quantities$suitability
    if (!is.null(suit)) {
      # A predator's rows are pred + nspp * (sex - 1), the r_idx of predation.hpp.
      has_suit <- vapply(seq_len(dl$nspp), function(p) {
        !isTRUE(all(suit[seq(p, dim(suit)[1], by = dl$nspp), , , , ] == 0))
      }, logical(1))
    } else {
      # A model stored by run_mse() keeps only the MSE quantities. Without the
      # suitability array, every predator with a diet record counts.
      has_suit <- seq_len(dl$nspp) %in% unique(as.integer(dl$diet_data$Pred))
    }
    emp       <- suit_mode == 0 & has_suit
    suit_end  <- rep_len(pmin(dl$suit_endyr, dl$endyr), dl$nspp)[emp]
    if (length(suit_end) && styr <= max(suit_end)) {
      stop("`styr` (", styr, ") must be after the empirical suitability window ",
           "(suit_endyr ", max(suit_end), "): removing fishing inside it would re-derive ",
           "the predation suitability the model was fit with.", call. = FALSE)
    }
  }

  # * Years for F = 0 ----
  proj_years <- styr:dl$projyr - dl$styr + 1
  fdevs_cols <- 1:ncol(object$estimated_params$log_F)
  fdevs_change <- which(fdevs_cols %in% proj_years)

  # * Set F to 0 ----
  object$estimated_params$log_F[,fdevs_change] <- replace(object$estimated_params$log_F[,fdevs_change], values = -999)

  # * Update fit ----
  # Build-only refit (projection unfished); SR-switch and suitability years clamp to endyr.
  estMode <- object$data_list$estimateMode
  object <- .refit_like(
    data_list        = object$data_list,
    inits            = object$estimated_params,
    estimateMode     = 3,
    getsd            = TRUE,
    srr_mse_switchyr = min(object$data_list$srr_mse_switchyr, object$data_list$endyr),
    suit_endyr       = pmin(object$data_list$suit_endyr, object$data_list$endyr))

  object$data_list$estimateMode <- estMode
  return(object)

}

