

# Display name -> returned column name for the mse_summary() performance
# metrics. The left-hand strings are what the metrics are called while they are
# being assembled (and what earlier versions returned); the right-hand names are
# what callers index. Prefixes say whose view a metric is: `om_` the operating
# model's truth, `em_` the estimation model's perception, unprefixed a quantity
# that has only one reading.
#
# Kept as one table so the two never drift, and attached to each returned frame
# as a "labels" attribute so a plot or table can print the long form.
.RCE_MSE_METRIC_NAMES <- c(
  "Average Catch"                                   = "avg_catch",
  "Catch IAV"                                       = "catch_iav",
  "P(Closed)"                                       = "p_closed",
  "Avg SSB Relative MSE"                            = "ssb_rmse_avg",
  "Avg terminal SSB Relative MSE"                   = "ssb_rmse_terminal",
  "EM: P(Fy > Flimit)"                              = "em_p_overfishing",
  "EM: P(SSB < SSBlimit)"                           = "em_p_overfished",
  "OM: P(Fy > Flimit)"                              = "om_p_overfishing",
  "OM: P(SSB < SSBlimit)"                           = "om_p_overfished",
  # Misclassification: the estimation model's status call disagrees with the
  # operating model's truth. "false_pos" = perceived overfishing/overfished when
  # it is not; "false_neg" = missed when it is.
  "EM: P(Fy > Flimit) but OM: P(Fy < Flimit)"       = "p_overfishing_false_pos",
  "EM: P(Fy < Flimit) but OM: P(Fy > Flimit)"       = "p_overfishing_false_neg",
  "EM: P(SSB < SSBlimit) but OM: P(SSB > SSBlimit)" = "p_overfished_false_pos",
  "EM: P(SSB > SSBlimit) but OM: P(SSB < SSBlimit)" = "p_overfished_false_neg",
  "OM: Terminal B"                                  = "om_terminal_biomass",
  "OM: Terminal SSB"                                = "om_terminal_ssb",
  "OM: Terminal Dynamic SB0"                        = "om_terminal_dynamic_sb0",
  "OM: Terminal SSB Depletion"                      = "om_terminal_depletion",
  "OM: Terminal SSB Depletion (Dynamic)"            = "om_terminal_depletion_dynamic",
  "OM: Average SSB Depletion"                       = "om_avg_depletion",
  # Counts of simulations in which SSB fell below the collapse cutoff, not
  # probabilities -- the names say `sims` so they are not read as the `p_` rows.
  "OM: SSB Collapse"                                = "om_sims_collapsed",
  "OM no F: SSB Collapse"                           = "om_no_f_sims_collapsed",
  "OM: SSB Collapse from F"                         = "om_sims_collapsed_from_f"
)

#' Rename the metric columns of an mse_summary() frame and record the old names
#'
#' Every metric column is renamed through `.RCE_MSE_METRIC_NAMES`; key columns
#' (`Species`, `Fleet_code`, `Fleet_name`) are already syntactic and pass
#' through. Errors on an unmapped metric rather than returning it under its
#' display name, so adding a metric without adding its entry cannot ship a
#' frame that is half-renamed.
#'
#' @param df A per-entity frame from mse_summary().
#' @return `df` with renamed columns and a `"labels"` attribute mapping each new
#'   name to its display string.
#' @keywords internal
#' @noRd
.rce_rename_mse_metrics <- function(df) {
  keys    <- intersect(c("Species", "Fleet_code", "Fleet_name"), names(df))
  metrics <- setdiff(names(df), keys)
  unmapped <- setdiff(metrics, names(.RCE_MSE_METRIC_NAMES))
  if (length(unmapped)) {
    stop("mse_summary(): no returned name for metric(s) ",
         paste(sQuote(unmapped), collapse = ", "),
         ". Add them to .RCE_MSE_METRIC_NAMES.", call. = FALSE)
  }
  names(df)[match(metrics, names(df))] <- .RCE_MSE_METRIC_NAMES[metrics]
  attr(df, "labels") <- stats::setNames(metrics, .RCE_MSE_METRIC_NAMES[metrics])
  df
}

#' Management strategy evaluation performance metric summary
#'
#' @param mse MSE runs from \code{\link{run_mse}} or \code{\link{load_mse}}
#' @param om_only If TRUE, report only operating-model (true) status and skip the
#'   estimation-model (EM) perception metrics.
#'
#' @return A named list, one element per entity dimension so each metric lives
#'   only where it applies (no NA padding):
#'   * `species`, a data.frame with one row per species (keyed by `Species`)
#'     of the conservation / status metrics: per-species average catch, catch
#'     inter-annual variability (IAV), and P(Closed) = the probability the
#'     fishery is closed (catch ~ 0); the relative mean-squared error of
#'     terminal and average SSB; the estimation-model- and operating-model-
#'     perceived overfishing / overfished probabilities (via
#'     \code{\link{build_hcr}}) and the probability each status is misclassified
#'     (estimation model disagrees with the operating-model truth); and terminal
#'     biomass, SSB, dynamic SB0, SSB depletion (equilibrium and dynamic),
#'     average SSB depletion, and SSB-collapse counts.
#'   * `fleet`, a data.frame with one row per fishery fleet (keyed by
#'     `Fleet_code` / `Fleet_name`) of average catch, catch IAV, and P(Closed).
#'   * `total`, a named numeric of the across-fleet total average catch and
#'     catch IAV.
#'   * `meta`, run provenance (`nsim`, `nspp`, `nflts`, `HCR`, projection-year
#'     range).
#'
#'   All metrics are averaged across projection years and simulations.
#'
#'   Metric columns hold syntactic names, so they can be typed without
#'   backticks: `avg_catch`, `catch_iav`, `p_closed`, `ssb_rmse_avg`,
#'   `ssb_rmse_terminal`, `em_p_overfishing`, `em_p_overfished`,
#'   `om_p_overfishing`, `om_p_overfished`, `p_overfishing_false_pos`,
#'   `p_overfishing_false_neg`, `p_overfished_false_pos`,
#'   `p_overfished_false_neg`, `om_terminal_biomass`, `om_terminal_ssb`,
#'   `om_terminal_dynamic_sb0`, `om_terminal_depletion`,
#'   `om_terminal_depletion_dynamic`, `om_avg_depletion`, `om_sims_collapsed`,
#'   `om_no_f_sims_collapsed`, `om_sims_collapsed_from_f`. An `om_` prefix is
#'   the operating model's truth and `em_` the estimation model's perception;
#'   the four `*_false_pos` / `*_false_neg` metrics are the probability the two
#'   disagree. The three `*_sims_collapsed` metrics are **counts of
#'   simulations**, not probabilities.
#'
#'   The three depletions divide by different references:
#'   - `om_terminal_depletion_dynamic`: terminal SSB over the OM's own
#'     `DynamicSB0`, its history with no fishing on the stock-recruit curve with
#'     the realized recruitment deviations and, under predation, the suitability
#'     fitted in the hindcast. Defined for single- and multispecies OMs. A
#'     species with input numbers-at-age (`estDynamics > 0`) keeps those numbers
#'     in the dynamic run, so its dynamic depletion is near 1 and says nothing
#'     about fishing.
#'   - `om_terminal_depletion`: terminal SSB over `SB0`, the equilibrium value
#'     for a single-species OM. Under `msmMode > 0` it is `MSSB0`, SSB at the end
#'     of a no-fishing projection, which `fit_mod()` derives only when a harvest
#'     control rule is present; without one (`HCR = "NoFishing"`) this is `NA`.
#'   - `om_avg_depletion`: the OM's own `ssb_depletion` averaged over the
#'     projection. Under predation with no HCR that series divides by SSB in the
#'     last projection year, not by an unfished reference.
#'
#'   The collapse metrics compare the OM with `OM_no_F`, the OM refit with no
#'   fishing after the original operating model's `endyr` ([remove_F()]).
#'
#'   Each frame holds a `"labels"` attribute mapping those names to the long
#'   display strings (e.g. `om_terminal_depletion_dynamic` ->
#'   `"OM: Terminal SSB Depletion (Dynamic)"`) for plots and tables:
#'   `attr(summ$species, "labels")`.
#'
#' @examples
#' \dontrun{
#' mse <- load_mse(dir = "mse_runs")
#' s <- mse_summary(mse)
#' s$species   # one row per species
#' s$fleet     # one row per fleet
#' s$total     # across the whole system
#' s$meta      # nsim, nspp, nflts, HCR
#' }
#' @export
#'
mse_summary <- function(mse, om_only = FALSE){

  ############################################
  ## Set up
  ############################################
  # HCR Switches (make length of nspp if not)
  extend_length <- function(x, nspp){
    if(length(x) == nspp){ return(x)}
    else {return(rep(x, nspp))}
  }

  # NA, not NaN, when every value is NA (e.g. no Flimit).
  mean_or_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

  ## Drop simulations that did not run to completion ----
  # run_mse() returns only a failure marker for those (no OM, no EMs), so there
  # is nothing to summarise and averaging over them would mix a partial
  # trajectory into every performance metric.
  failed <- vapply(mse, function(x) isFALSE(x$use_sim), logical(1))
  if (any(failed)) {
    warning(sum(failed), " of ", length(mse), " simulations did not complete ",
            "and are excluded: ",
            paste(names(mse)[failed], collapse = ", "), call. = FALSE)
    mse <- mse[!failed]
  }
  if (!length(mse)) {
    stop("No simulations completed; there is nothing to summarise.", call. = FALSE)
  }

  ## Simulations whose unfished reference run failed ----
  # The unfished run is a separate fit and can fail on its own, in which case
  # run_mse() keeps the simulation but leaves OM_no_F NULL. Such a simulation is
  # valid for everything except the no-F comparisons, so it is excluded from
  # those specifically -- numerator and denominator both. Left in,
  # `sum(NULL < 1000) > 0` answers FALSE, counting a missing reference run as
  # "did not collapse" and biasing the collapse metrics low.
  has_no_f <- vapply(mse, function(x) !is.null(x$OM_no_F), logical(1))
  mse_no_f <- mse[has_no_f]
  if (any(!has_no_f)) {
    warning(sum(!has_no_f), " of ", length(mse), " simulations have no unfished ",
            "reference run; the 'OM no F' and 'SSB Collapse from F' metrics are ",
            "over the remaining ", length(mse_no_f), ".", call. = FALSE)
  }

  ## OM dimensions ----
  # - determined from OM sim 1
  # - should be the same as for the EM
  msmMode <- mse[[1]]$OM$data_list$msmMode
  nspp <- mse[[1]]$OM$data_list$nspp
  nsex <- mse[[1]]$OM$data_list$nsex
  flt_type <- mse[[1]]$OM$data_list$fleet_control$Fleet_type
  flts <- mse[[1]]$OM$data_list$fleet_control$Fleet_code[which(flt_type == "Fishery")]
  nflts = length(flts)
  flt_spp <- mse[[1]]$OM$data_list$fleet_control$Species
  styr <- mse[[1]]$OM$data_list$styr
  endyr <- mse[[1]]$EM[[1]]$data_list$endyr
  projyr <- min(c(mse[[1]]$EM[[1]]$data_list$projyr, mse[[1]]$OM$data_list$projyr))
  projyrs <- (endyr+1):projyr

  ## HCR ----
  # - Normalize HCR to integer code. build_hcr() accepts either integer or
  #   string alias (e.g. "NPFMC"), and fit_mod()/run_mse() may store either
  #   form depending on processing path. All downstream comparisons in this
  #   function use integer codes (HCR == 5, etc.), so coerce once here.
  .normalize_hcr <- function(x) {
    if (is.numeric(x) || (is.character(x) && !is.na(suppressWarnings(as.integer(x))))) {
      return(as.integer(x))
    }
    idx <- match(as.character(x), names(hcr_map))
    if (is.na(idx)) {
      stop("Unrecognized HCR value: '", x, "'. Expected one of ",
           paste(names(hcr_map), collapse = ", "),
           " or integer code 0:", max(hcr_map), ".")
    }
    as.integer(unname(hcr_map[idx]))
  }
  HCR <- .normalize_hcr(mse[[1]]$EM[[1]]$data_list$HCR)
  DynamicHCR <- as.logical(mse[[1]]$EM[[1]]$data_list$DynamicHCR)
  Ptarget <- extend_length(mse[[1]]$EM[[1]]$data_list$Ptarget, nspp)
  Plimit <- extend_length(mse[[1]]$EM[[1]]$data_list$Plimit, nspp)
  Alpha <- extend_length(mse[[1]]$EM[[1]]$data_list$Alpha, nspp)
  # -- HCR = 0: No catch - Params off
  # -- HCR = 1: Constant catch - Params off
  # -- HCR = 2: Constant input F - Params off
  # -- HCR = 3: F that achieves X% of SSB0 in the end of the projection - Ftarget on
  # -- HCR = 4: Constant target Fspr - Ftarget on
  # -- HCR = 5: NPFMC Tier 3 - Flimit and Ftarget on
  # -- HCR = 6: PFMC Cat 1 - Flimit on
  # -- HCR = 7: SESSF Tier 1 - Flimit and Ftarget on

  ## MSE specifications
  nsim <- length(mse)

  ## MSE Output ----
  # - Catch is by fleet
  mse_summary <- data.frame(matrix(NA, nrow = nflts+nspp+1, ncol = 25))
  colnames(mse_summary) <- c("Species",
                             "Fleet_name",
                             "Fleet_code",
                             "Average Catch",
                             "Catch IAV",
                             "P(Closed)",
                             "Avg SSB Relative MSE",
                             "Avg terminal SSB Relative MSE",
                             "EM: P(Fy > Flimit)",
                             "EM: P(SSB < SSBlimit)",
                             "OM: P(Fy > Flimit)",
                             "OM: P(SSB < SSBlimit)",
                             "EM: P(Fy > Flimit) but OM: P(Fy < Flimit)",
                             "EM: P(Fy < Flimit) but OM: P(Fy > Flimit)",
                             "EM: P(SSB < SSBlimit) but OM: P(SSB > SSBlimit)",
                             "EM: P(SSB > SSBlimit) but OM: P(SSB < SSBlimit)",
                             # "OM: Recovery Time",
                             "OM: Terminal B",
                             "OM: Terminal SSB",
                             "OM: Terminal Dynamic SB0",
                             "OM: Terminal SSB Depletion",
                             "OM: Terminal SSB Depletion (Dynamic)",
                             "OM: Average SSB Depletion",
                             "OM no F: SSB Collapse",
                             "OM: SSB Collapse",
                             "OM: SSB Collapse from F")
  mse_summary$Fleet_name <- c(rep("Conservation metric", nspp), mse[[1]]$OM$data_list$fleet_control$Fleet_name[flts], "All")
  mse_summary$Fleet_code <- c(rep(NA, nspp), mse[[1]]$OM$data_list$fleet_control$Fleet_code[flts], "All")
  mse_summary$Species <- c(mse[[1]]$OM$data_list$spnames,
                           mse[[1]]$OM$data_list$spnames[mse[[1]]$OM$data_list$fleet_control$Species[flts]],
                           "All")


  ## Catch performance metrics by fleet ----
  # - Average catch
  # - Catch IAV
  # - P(Closed)
  # Species with input numbers-at-age are projected at F = 0: no catch to vary or close.
  fixed_n <- (mse[[1]]$OM$data_list$estDynamics %||% rep(0, nspp)) > 0
  for(i in 1:nflts){
    flt = flts[i]

    # * Mean catch ----
    mse_summary$`Average Catch`[i+nspp] <- mean(
      sapply(mse, function(x)
        x$OM$data_list$catch_data |>
          filter(Fleet_code == flt & Year %in% projyrs) |>
          pull(Catch)
      ), na.rm = TRUE)

    # * Catch IAV ----
    catch_list_tmp <- lapply(mse, function(x)
      x$OM$data_list$catch_data |>
        dplyr::filter(Fleet_code == flt & Year %in% projyrs) |>
        pull(Catch)
    )

    # -- Average across simulations by fleet
    mse_summary$`Catch IAV`[i+nspp] = 0
    for(sim in 1:nsim){
      iav_tmp <- sum((lag(catch_list_tmp[[sim]], 1) - catch_list_tmp[[sim]])^2, na.rm = TRUE)/(length(projyrs) - 1) # Squared difference
      iav_tmp <- sqrt(iav_tmp) / (sum(catch_list_tmp[[sim]], na.rm = TRUE)/ length(projyrs))/nsim # Divide by mean
      mse_summary$`Catch IAV`[i+nspp] <- mse_summary$`Catch IAV`[i+nspp] + iav_tmp
    }

    # * P(Closed) ----
    mse_summary$`P(Closed)`[i+nspp] <- sum(
      sapply(catch_list_tmp, function(x)
        length(which(x < 1)) # Using less than 1 here just in case super small catches and fishery is effectively close
        /length(x)))/nsim

    if (isTRUE(fixed_n[flt_spp[flt]])) {
      mse_summary$`Catch IAV`[i+nspp] <- mse_summary$`P(Closed)`[i+nspp] <- NA_real_
    }
  }


  ## Catch performance metrics by species ----
  # - Average catch
  # - Catch IAV
  # - P(Closed)
  #
  # A multispecies model can carry a species with no fishery at all -- a
  # predator included for its consumption, as arrowtooth and sablefish are in
  # the hake model. Such a species has no `catch_data` rows, so each quantity
  # below reduces to an operation on a zero-length vector: `pull(Catch)` returns
  # numeric(0), sapply() over those returns a LIST rather than a numeric, and
  # mean() of a list warns "argument is not numeric or logical: returning NA".
  # The two ratio metrics were worse -- they divide by length(x) == 0 and
  # returned NaN with no warning at all.
  #
  # Handled explicitly instead. Catch on an unfished species is zero, which is a
  # real answer; its catch variability and probability-of-closure are not
  # defined, because there is no fishery to vary or to close.
  #
  # A species counts as fished if `fleet_control` declares a fishery for it OR
  # `catch_data` carries any row for it. The fleet table is the authoritative
  # statement -- it distinguishes "has a fishery that landed nothing" from "has
  # no fishery" -- but it is not used alone: `Fleet_type` reaches here as the
  # normalized "Fishery"/"Survey" string, and were it ever to arrive as the raw
  # integer code the comparison would match nothing and EVERY species would be
  # silently reported as unfished. Falling back to the catch rows makes that
  # failure impossible.
  fished_spp <- union(unique(flt_spp[flt_type == "Fishery"]),
                      unique(mse[[1]]$OM$data_list$catch_data$Species))

  for(sp in 1:nspp){

    # - Per-simulation projection-year catch. lapply()+unlist() rather than
    #   sapply() so an empty result stays numeric(0) instead of degenerating to
    #   a list; with rows present the pooled values are identical to before.
    catch_by_sim <- lapply(mse, function(x)
      x$OM$data_list$catch_data |>
        filter(Species == sp & Year %in% projyrs) |>
        pull(Catch))

    # - Catch IAV ----
    catch_list_tmp <- suppressMessages(lapply(mse, function(x)
      x$OM$data_list$catch_data |>
        filter(Species == sp & Year %in% projyrs) |>
        group_by(Year) |>
        summarise(Catch = sum(Catch)) |>
        pull(Catch)
    )) # Sum catch across fleets within a year

    if (!(sp %in% fished_spp) || fixed_n[sp] || sum(lengths(catch_by_sim)) == 0L) {
      if (sp %in% fished_spp && !fixed_n[sp] && sum(lengths(catch_by_sim)) == 0L) {
        # Has a fishery, but nothing landed in the projection years: a data
        # problem rather than a modelling choice, so do not report it as zero.
        warning("Species ", sp, " has a fishery but no catch in the projection ",
                "years (", min(projyrs), "-", max(projyrs), "); its catch ",
                "metrics are NA.", call. = FALSE)
        mse_summary$`Average Catch`[sp] <- NA_real_
      } else {
        mse_summary$`Average Catch`[sp] <- 0
      }
      mse_summary$`Catch IAV`[sp] <- NA_real_
      mse_summary$`P(Closed)`[sp] <- NA_real_
      next
    }

    # * Mean catch ----
    mse_summary$`Average Catch`[sp] <- mean(unlist(catch_by_sim), na.rm = TRUE)

    # - Average across simulations
    mse_summary$`Catch IAV`[sp] <- 0 # Initialize
    for(sim in 1:nsim){
      iav_tmp <- sum((lag(catch_list_tmp[[sim]], 1) - catch_list_tmp[[sim]])^2, na.rm = TRUE)/(length(projyrs) - 1) # Squared difference
      iav_tmp <- sqrt(iav_tmp) / (sum(catch_list_tmp[[sim]], na.rm = TRUE)/ length(projyrs))/nsim # Divide by mean
      mse_summary$`Catch IAV`[sp] <- mse_summary$`Catch IAV`[sp] + iav_tmp
    }

    # * P(Closed) ----
    mse_summary$`P(Closed)`[sp] <- sum(
      sapply(catch_list_tmp, function(x)
        length(which(x < 1)) # Using less than 1 here just in case super small catches and fishery is effectively close
        /length(x)))/nsim
  }


  ## Catch performance metrics across species ----
  # - Average catch
  # - Catch IAV
  # - P(Closed)
  catch_list_tmp <- suppressMessages(lapply(mse, function(x)
    x$OM$data_list$catch_data |>
      filter(Year %in% projyrs) |>
      group_by(Year) |>
      summarise(Catch = sum(Catch)) |>
      pull(Catch)
  )
  ) # Sum catch across species

  # * Mean catch ----
  mse_summary$`Average Catch`[nspp + nflts + 1] <- mean(unlist(catch_list_tmp), na.rm = TRUE)

  # * Catch IAV ----
  # -- Average across simulations
  mse_summary$`Catch IAV`[nspp + nflts + 1] <- 0 # Initialize
  for(sim in 1:nsim){
    iav_tmp <- sum((lag(catch_list_tmp[[sim]], 1) - catch_list_tmp[[sim]])^2, na.rm = TRUE)/(length(projyrs) - 1) # Squared difference
    iav_tmp <- sqrt(iav_tmp) / (sum(catch_list_tmp[[sim]], na.rm = TRUE)/ length(projyrs))/nsim # Divide by mean
    mse_summary$`Catch IAV`[nspp + nflts + 1] <- mse_summary$`Catch IAV`[nspp + nflts + 1] + iav_tmp
  }


  ## Conservation performance metrics ----
  # - Avg terminal SSB MSE
  # - EM: P(Fy > Flimit)
  # - EM: P(SSB < SSBlimit)
  # - OM: P(Fy > Flimit)
  # - OM: P(SSB < SSBlimit)
  # - EM: P(Fy > Flimit) but OM: P(Fy < Flimit)
  # - EM: P(Fy < Flimit) but OM: P(Fy > Flimit)
  # - EM: P(SSB < SSBlimit) but OM: P(SSB > SSBlimit)
  # - EM: P(SSB > SSBlimit) but OM: P(SSB < SSBlimit)
  # - OM: Terminal Depletion Relative to equilibrium SB0
  # - OM: Terminal Depletion Relative to dynamic SB0
  # - EM: Average age-1 M_at_age
  # - EM: Variance of age-1 M_at_age

  # -- Tier 3 for single-species models
  # - Produces vectors of Flimits given depletion and input Flimit (Fspr)
  # - Note, it doesnt have Plimit because thats for cod
  flimit_tier3_fun <- function(ssb_depletion, ssb, SBF, plimit, alpha, Flimit){
    # No Flimit, no Tier 3 limit.
    if (length(Flimit) != 1 || is.na(Flimit)) return(rep(NA_real_, length(ssb)))
    tier3_flimit <- c()
    for(i in 1:length(ssb)){

      # Tier-3 HCR
      if(ssb[i] >= SBF[i]){
        tier3_flimit[i] = Flimit
      }else if(ssb[i] < SBF[i] & ssb[i] > alpha * SBF[i]){
        tier3_flimit[i] = Flimit * (ssb[i]/SBF[i] - alpha)/(1-alpha)
      }else{
        tier3_flimit[i] = 0
      }

      # If below 20%
      if(ssb_depletion[i] < plimit){
        tier3_flimit[i] = 0
      }
    }
    return(tier3_flimit)
  }


  # Helper: SSB-limit threshold (depletion or absolute) given HCR
  # - returns list(value, is_depletion) so callers can compare against
  #   either ssb_depletion or absolute ssb depending on HCR
  ssb_limit_thresh <- function(HCR, DynamicHCR, Ptarget_sp, Plimit_sp, SBF_val, DynamicSBF_val){
    if(HCR == 2)                                  list(val = 0.5 * SBF_val,        is_dep = FALSE)
    else if(HCR == 4 & !DynamicHCR)               list(val = 0.5 * SBF_val,        is_dep = FALSE)
    else if(HCR == 4 &  DynamicHCR)               list(val = 0.5 * DynamicSBF_val, is_dep = FALSE)
    else if(HCR == 5 & !DynamicHCR)               list(val = 0.5 * SBF_val,        is_dep = FALSE)
    else if(HCR == 5 &  DynamicHCR)               list(val = 0.5 * DynamicSBF_val, is_dep = FALSE)
    else if(HCR == 6 & Ptarget_sp == 0.4)         list(val = 0.25,                 is_dep = TRUE)
    else if(HCR == 6)                             list(val = 0.5 * Plimit_sp,      is_dep = TRUE)
    else                                          list(val = Plimit_sp,            is_dep = TRUE)
  }

  for(sp in 1:nspp){

    ## Perceived status (EM) and OM at matching assess years ----
    if(!om_only){
      em_f_flimit          <- c()
      em_sb_sblimit        <- c()
      om_f_flimit_assess   <- c()
      om_sb_sblimit_assess <- c()
      em_ssb_assess        <- c()
      om_ssb_assess        <- c()

      for(sim in 1:length(mse)){
        em_list <- mse[[sim]]$EM
        om_q    <- mse[[sim]]$OM$quantities

        for(em in 2:length(em_list)){ # First EM is "conditioned" model

          # Terminal year of intermediate assessment
          em_endyr_sim <- em_list[[em]]$data_list$endyr
          end_yr_col   <- em_endyr_sim - styr + 1
          em_q         <- em_list[[em]]$quantities

          # * EM: P(F > Flimit) ----
          # - Tier-3 aware: replace base Flimit with depletion-adjusted Flimit when HCR == 5
          em_Flimit <- em_q$Flimit[sp]
          if(HCR == 5){
            em_SBF_yr <- if(DynamicHCR) em_q$DynamicSBF[sp, end_yr_col] else em_q$SBF[sp, end_yr_col]
            em_Flimit <- flimit_tier3_fun(
              ssb_depletion = em_q$ssb_depletion[sp, end_yr_col],
              ssb           = em_q$ssb[sp, end_yr_col],
              SBF           = em_SBF_yr,
              plimit = Plimit[sp], alpha = Alpha[sp],
              Flimit = em_q$Flimit[sp]
            )
          }
          em_f_flimit <- c(em_f_flimit, em_q$F_spp[sp, end_yr_col] > em_Flimit)

          # * EM: P(SSB < SSBlimit) ----
          # HCR 2 ("ConstantF") fishes at a fixed rate and carries no depletion
          # target, so it takes the absolute 0.5 * SBF that ssb_limit_thresh()
          # already gives the OM -- read from the same helper so the two cannot
          # drift, and gated on msmMode as the OM arm is. The cross-tab then
          # compares the EM and the OM on one criterion. The remaining arms are
          # depletions the run was configured with.
          # Under msmMode > 0 both sides fall through to Plimit, which defaults
          # to 0 and so reports the stock as never overfished. That is the OM's
          # own multispecies rule (see the msmMode branch below), so the two
          # still agree; set Plimit explicitly for a multispecies ConstantF run.
          if(HCR == 2 & msmMode == 0){
            em_thr <- ssb_limit_thresh(HCR, DynamicHCR, Ptarget[sp], Plimit[sp],
                                       SBF_val        = em_q$SBF[sp, end_yr_col],
                                       DynamicSBF_val = em_q$DynamicSBF[sp, end_yr_col])
            # Branch on the helper's own scale flag, as the OM arm does, so a
            # later change to that arm cannot leave this one comparing an
            # absolute biomass against a depletion.
            em_below <- if(em_thr$is_dep) em_q$ssb_depletion[sp, end_yr_col] < em_thr$val
                        else              em_q$ssb[sp,            end_yr_col] < em_thr$val
          } else {
            em_thresh <- if(HCR == 4)                            0.5 * Ptarget[sp]
                         else if(HCR == 5)                       0.5 * Plimit[sp]
                         else if(HCR == 6 & Ptarget[sp] == 0.40) 0.25
                         else if(HCR == 6)                       0.5 * Ptarget[sp]
                         else                                    Plimit[sp]
            em_below <- em_q$ssb_depletion[sp, end_yr_col] < em_thresh
          }
          em_sb_sblimit <- c(em_sb_sblimit, em_below)

          # * OM at the SAME assess year (for cross-tab and SSB MSE) ----
          # - OM F > Flimit (Tier-3 aware for single-species mode)
          om_Flimit <- om_q$Flimit[sp]
          if(msmMode == 0 & HCR == 5){
            om_SBF_yr <- if(DynamicHCR) om_q$DynamicSBF[sp, end_yr_col] else om_q$SBF[sp, end_yr_col]
            om_Flimit <- flimit_tier3_fun(
              ssb_depletion = om_q$ssb_depletion[sp, end_yr_col],
              ssb           = om_q$ssb[sp, end_yr_col],
              SBF           = om_SBF_yr,
              plimit = Plimit[sp], alpha = Alpha[sp],
              Flimit = om_q$Flimit[sp]
            )
          }
          om_f_flimit_assess <- c(om_f_flimit_assess, om_q$F_spp[sp, end_yr_col] > om_Flimit)

          # - OM SSB < SSBlimit (uses absolute SBF for SPR-based HCRs, depletion otherwise)
          if(msmMode == 0){
            om_thr <- ssb_limit_thresh(HCR, DynamicHCR, Ptarget[sp], mse[[sim]]$OM$data_list$Plimit[sp],
                                       SBF_val        = om_q$SBF[sp, end_yr_col],
                                       DynamicSBF_val = om_q$DynamicSBF[sp, end_yr_col])
            om_below <- if(om_thr$is_dep) om_q$ssb_depletion[sp, end_yr_col] < om_thr$val
                        else              om_q$ssb[sp,            end_yr_col] < om_thr$val
          } else {
            om_below <- om_q$ssb_depletion[sp, end_yr_col] < mse[[sim]]$OM$data_list$Plimit[sp]
          }
          om_sb_sblimit_assess <- c(om_sb_sblimit_assess, om_below)

          # - Paired SSB values for relative MSE across assessments
          em_ssb_assess <- c(em_ssb_assess, em_q$ssb[sp, end_yr_col])
          om_ssb_assess <- c(om_ssb_assess, om_q$ssb[sp, end_yr_col])
        }
      }

      ## Perceived status (averages over all assess-year observations)
      mse_summary$`EM: P(Fy > Flimit)`[sp]    <- mean_or_na(em_f_flimit)
      mse_summary$`EM: P(SSB < SSBlimit)`[sp] <- mean_or_na(em_sb_sblimit)
    }


    ## Actual status across all projection years ----
    # * OM: P(F > Flimit) ----
    om_f_flimit <- lapply(mse, function(x) x$OM$quantities$F_spp[sp, (projyrs - styr + 1)] > (x$OM$quantities$Flimit[sp]))

    # - Tier 3
    if(msmMode == 0 & HCR == 5 & DynamicHCR == FALSE){
      om_f_flimit <- lapply(mse, function(x) x$OM$quantities$F_spp[sp, (projyrs - styr + 1)] >
                              flimit_tier3_fun(
                                ssb_depletion = x$OM$quantities$ssb_depletion[sp,(projyrs - styr + 1)],
                                ssb = x$OM$quantities$ssb[sp,(projyrs - styr + 1)],
                                SBF = x$OM$quantities$SBF[sp,(projyrs - styr + 1)],
                                Plimit[sp], Alpha[sp], Flimit = x$OM$quantities$Flimit[sp]
                              )
      )
    }

    # - Dynamic Tier 3
    if(msmMode == 0 & HCR == 5 & DynamicHCR == TRUE){
      om_f_flimit <- lapply(mse, function(x) x$OM$quantities$F_spp[sp, (projyrs - styr + 1)] >
                              flimit_tier3_fun(
                                ssb_depletion = x$OM$quantities$ssb_depletion[sp,(projyrs - styr + 1)],
                                ssb = x$OM$quantities$ssb[sp,(projyrs - styr + 1)],
                                SBF = x$OM$quantities$DynamicSBF[sp,(projyrs - styr + 1)],
                                Plimit[sp], Alpha[sp], Flimit = x$OM$quantities$Flimit[sp]
                              )
      )
    }

    om_f_flimit <- unlist(om_f_flimit)
    mse_summary$`OM: P(Fy > Flimit)`[sp] <- mean_or_na(om_f_flimit)


    # * OM: P(SSB < SSBlimit) ----
    # - default (multi-species or HCR with no SPR ref): depletion-based
    om_sb_sblimit <- lapply(mse, function(x) x$OM$quantities$ssb_depletion[sp, (projyrs - styr + 1)] < x$OM$data_list$Plimit[sp])

    if(msmMode == 0){
      proj_idx <- projyrs - styr + 1

      # - Avg F SPR based (HCR 2)
      if(HCR == 2){
        om_sb_sblimit <- lapply(mse, function(x) x$OM$quantities$ssb[sp, proj_idx] < 0.5 * x$OM$quantities$SBF[sp, proj_idx])
      }

      # - New England SPR based (HCR 4)
      if(HCR == 4 & !DynamicHCR){
        om_sb_sblimit <- lapply(mse, function(x) x$OM$quantities$ssb[sp, proj_idx] < 0.5 * x$OM$quantities$SBF[sp, proj_idx])
      } else if(HCR == 4 & DynamicHCR){
        om_sb_sblimit <- lapply(mse, function(x) x$OM$quantities$ssb[sp, proj_idx] < 0.5 * x$OM$quantities$DynamicSBF[sp, proj_idx])
      }

      # - Tier 3 SPR based (HCR 5)
      if(HCR == 5 & !DynamicHCR){
        om_sb_sblimit <- lapply(mse, function(x) x$OM$quantities$ssb[sp, proj_idx] < 0.5 * x$OM$quantities$SBF[sp, proj_idx])
      } else if(HCR == 5 & DynamicHCR){
        om_sb_sblimit <- lapply(mse, function(x) x$OM$quantities$ssb[sp, proj_idx] < 0.5 * x$OM$quantities$DynamicSBF[sp, proj_idx])
      }

      # - Cat 1 depletion based (HCR 6)
      if(HCR == 6){
          om_sb_sblimit <- lapply(mse, function(x) x$OM$quantities$ssb_depletion[sp, proj_idx] < Ptarget[sp] * 0.5)
        if(Ptarget[sp] == 0.4){
          om_sb_sblimit <- lapply(mse, function(x) x$OM$quantities$ssb_depletion[sp, proj_idx] < 0.25)
        }
      }
    }

    om_sb_sblimit <- unlist(om_sb_sblimit)
    mse_summary$`OM: P(SSB < SSBlimit)`[sp] <- mean_or_na(om_sb_sblimit)


    ## * Perceived status relative to actual status (paired at assess years) ----
    if(!om_only){
      # NA when either side is NA.
      both_f <- !is.na(em_f_flimit) & !is.na(om_f_flimit_assess)
      mse_summary$`EM: P(Fy > Flimit) but OM: P(Fy < Flimit)`[sp]      <- mean_or_na(ifelse(both_f, em_f_flimit == 1 & om_f_flimit_assess == 0, NA))
      mse_summary$`EM: P(Fy < Flimit) but OM: P(Fy > Flimit)`[sp]      <- mean_or_na(ifelse(both_f, em_f_flimit == 0 & om_f_flimit_assess == 1, NA))
      both_sb <- !is.na(em_sb_sblimit) & !is.na(om_sb_sblimit_assess)
      mse_summary$`EM: P(SSB < SSBlimit) but OM: P(SSB > SSBlimit)`[sp] <- mean_or_na(ifelse(both_sb, em_sb_sblimit == 1 & om_sb_sblimit_assess == 0, NA))
      mse_summary$`EM: P(SSB > SSBlimit) but OM: P(SSB < SSBlimit)`[sp] <- mean_or_na(ifelse(both_sb, em_sb_sblimit == 0 & om_sb_sblimit_assess == 1, NA))
    }


    # * Bias in SSB ----
    if(!om_only){
      # - terminal projection year (last EM vs OM)
      terminal_ssb_om <- sapply(mse, function(x) x$OM$quantities$ssb[sp, (projyr - styr + 1)])
      terminal_ssb_em <- sapply(mse, function(x) x$EM[[length(x$EM)]]$quantities$ssb[sp, (projyr - styr + 1)])
      mse_summary$`Avg terminal SSB Relative MSE`[sp] <- mean((terminal_ssb_em - terminal_ssb_om)^2 / terminal_ssb_om^2, na.rm = TRUE)

      # - across all intermediate assessments (EM terminal SSB vs OM SSB at the same year)
      mse_summary$`Avg SSB Relative MSE`[sp] <- mean((em_ssb_assess - om_ssb_assess)^2 / om_ssb_assess^2, na.rm = TRUE)
    }

    # * OM: Terminal B, SSB, depletion ----
    terminal_b_om <- sapply(mse, function(x) x$OM$quantities$biomass[sp, (projyr - styr + 1)])
    terminal_ssb_om <- sapply(mse, function(x) x$OM$quantities$ssb[sp, (projyr - styr + 1)])

    # Dynamic SB0: the OM's history with no fishing (fitted curve, realized deviations, suitability).
    terminal_dynamic_sb0_om <- sapply(mse, function(x) x$OM$quantities$DynamicSB0[sp, (projyr - styr + 1)])

    if(mse[[1]]$OM$data_list$msmMode == 0){
      terminal_sb0_om <- sapply(mse, function(x) x$OM$quantities$SB0[sp, (projyr - styr + 1)])
    }

    if(mse[[1]]$OM$data_list$msmMode > 0){
      # Terminal year, as in the single-species arm above. The multispecies
      # SB0 is the same in every year -- the model overwrites its own
      # derivation with the `MSSB0` input -- so this reads the same number a
      # bare `SB0[sp]` did, and keeps reading the right one if that changes.
      terminal_sb0_om <- sapply(mse, function(x) x$OM$quantities$SB0[sp, (projyr - styr + 1)])

      # An unfished reference is only derived for a run carrying a harvest
      # control rule, so under `HCR = "NoFishing"` MSSB0 is still the
      # placeholder. Dividing by it reports SSB/999 as a depletion -- on the
      # Pacific hake three-species model that read 2.7e3. Not defined here.
      #
      # Read the per-species flag fit_mod() recorded, not the value. Comparing
      # the reported SB0 to the placeholder constant would also null a
      # legitimately-derived 999 mt, and could not see a workbook that supplied
      # its own MSSB0. A fit from before the flag existed has no field, so fall
      # back to the value test for those.
      derived <- vapply(mse, function(x) {
        d <- x$OM$data_list$MSSB0_derived
        if (is.null(d)) NA else isTRUE(d[sp])
      }, logical(1))
      undefined <- if (all(is.na(derived))) {
        terminal_sb0_om == .RCE_MSSB0_PLACEHOLDER
      } else {
        !derived | is.na(derived)
      }
      terminal_sb0_om[undefined] <- NA_real_
    }

    mse_summary$`OM: Terminal B`[sp] <- mean(terminal_b_om)
    mse_summary$`OM: Terminal SSB`[sp] <- mean(terminal_ssb_om)
    mse_summary$`OM: Terminal Dynamic SB0`[sp] <- mean(terminal_dynamic_sb0_om)
    mse_summary$`OM: Terminal SSB Depletion`[sp] <- mean(terminal_ssb_om/terminal_sb0_om)
    mse_summary$`OM: Terminal SSB Depletion (Dynamic)`[sp] <- mean(terminal_ssb_om/terminal_dynamic_sb0_om)


    # - OM: Average SSB depletion
    sb_depletion <- lapply(mse, function(x) x$OM$quantities$ssb_depletion[sp, (projyrs - styr + 1)])
    sb_depletion <- unlist(sb_depletion)
    mse_summary$`OM: Average SSB Depletion`[sp] <- mean(sb_depletion)

    # * OM: Collapse ----
    mse_summary$`OM no F: SSB Collapse`[sp] <- if (length(mse_no_f)) {
      sum(sapply(mse_no_f, function(x) sum(x$OM_no_F$quantities$ssb[sp, (projyrs - styr + 1)] < 1000) > 0))
    } else NA_integer_
    mse_summary$`OM: SSB Collapse`[sp] <- sum(sapply(mse, function(x) sum(x$OM$quantities$ssb[sp, (projyrs - styr + 1)] < 1000) > 0))

    # -- OM without F is above cutoff, but OM with F is below
    mse_summary$`OM: SSB Collapse from F`[sp] <- if (length(mse_no_f)) {
      sum(sapply(mse_no_f,
                 function(x) sum((x$OM$quantities$ssb[sp, (projyrs - styr + 1)] < 1000) *
                                   (x$OM_no_F$quantities$ssb[sp, (projyrs - styr + 1)] > 1000)) > 0))
    } else NA_integer_
  }

  ############################################
  ## Reshape into per-entity tidy frames
  ############################################
  # The metrics above were accumulated into a single stacked frame (species
  # rows, then fishery-fleet rows, then an "All" row), so every column was
  # NA-padded outside the entity it applies to and the same column name (e.g.
  # "Average Catch") meant a species value, a fleet value, or a total depending
  # on the row. Return them split by entity instead: one frame per dimension,
  # each with a labelled key column and no padding.
  metric_cols   <- setdiff(colnames(mse_summary),
                           c("Species", "Fleet_name", "Fleet_code"))
  fleet_metrics <- c("Average Catch", "Catch IAV", "P(Closed)")
  species_rows  <- seq_len(nspp)
  fleet_rows    <- nspp + seq_len(nflts)
  total_row     <- nspp + nflts + 1L

  species <- data.frame(
    Species = mse_summary$Species[species_rows],
    mse_summary[species_rows, metric_cols, drop = FALSE],
    check.names = FALSE, row.names = NULL, stringsAsFactors = FALSE)

  fleet <- data.frame(
    # Fleet_code was character in the stacked frame only because of the "All"
    # total row (now in `total`); the per-fleet frame carries the integer codes.
    Fleet_code = as.integer(mse_summary$Fleet_code[fleet_rows]),
    Fleet_name = mse_summary$Fleet_name[fleet_rows],
    mse_summary[fleet_rows, fleet_metrics, drop = FALSE],
    check.names = FALSE, row.names = NULL, stringsAsFactors = FALSE)

  total <- unlist(mse_summary[total_row, c("Average Catch", "Catch IAV")])

  ## Syntactic metric names ----
  # The working names above are display strings ("OM: Terminal SSB Depletion
  # (Dynamic)"), which a caller can only reach through backticks or [[ ]].
  # Rename on the way out so the returned frames are ordinary data frames whose
  # columns can be typed, and keep the display strings on a "labels" attribute
  # so a plot or a SAFE table can still render them.
  species <- .rce_rename_mse_metrics(species)
  fleet   <- .rce_rename_mse_metrics(fleet)
  # `total` is a named vector, not a frame, so it goes through the same helper
  # as a one-row frame rather than being renamed by a bare lookup -- an unmapped
  # name there would otherwise become NA silently instead of erroring.
  total <- unlist(.rce_rename_mse_metrics(
    as.data.frame(as.list(total), check.names = FALSE)))

  # Still the same ragged list -- $species, $fleet, $total and $meta are
  # unchanged. The class only adds a print method, so the object says what it
  # holds instead of dumping four blocks of differing shape; as.data.frame() on
  # it still errors, which is the caller's bug and is now said out loud.
  structure(
    list(
      species = species,
      fleet   = fleet,
      total   = total,
      meta    = list(nsim = nsim, nspp = nspp, nflts = nflts, HCR = HCR,
                     proj_years = c(first = min(projyrs), last = max(projyrs)))
    ),
    class = "Rceattle_mse_summary")
}


#' Print method for an MSE summary
#'
#' @description Says what the summary holds. The object is deliberately ragged,
#' per-species, per-fleet and whole-system metrics have different shapes and
#' cannot share one frame, so it reports the blocks and their dimensions rather
#' than printing them end to end.
#'
#' @param x An `"Rceattle_mse_summary"` object from [mse_summary()].
#' @param ... Currently unused.
#' @return `x`, invisibly.
#' @export
print.Rceattle_mse_summary <- function(x, ...) {
  mt <- x$meta

  # An undefined depletion is the one thing here worth a severity: it means no
  # unfished reference was derived (HCR = "NoFishing" leaves MSSB0 at its
  # placeholder), so the column is NA rather than a ratio against 999 mt.
  dep <- suppressWarnings(as.numeric(
    x$species[["om_terminal_depletion"]] %||% NA_real_))
  sev <- if (length(dep) && all(is.na(dep))) "NOTE" else "OK"

  .rce_diag_header(
    "MSE summary", sev,
    paste0(mt$nsim, " simulation(s), ", mt$nspp, " species, ", mt$nflts,
           " fleet(s); HCR ", mt$HCR, "; projection ",
           mt$proj_years[["first"]], "-", mt$proj_years[["last"]]))

  dims <- function(v) if (is.null(v)) "-" else
    if (is.null(dim(v))) paste0(length(v), " value(s)") else
      paste0(nrow(v), " x ", ncol(v))
  cat("  $species  ", dims(x$species), "\n", sep = "")
  cat("  $fleet    ", dims(x$fleet), "\n", sep = "")
  cat("  $total    ", dims(x$total), "\n", sep = "")
  cat("  $meta     nsim, nspp, nflts, HCR, proj_years\n")

  if (identical(sev, "NOTE")) {
    cat("  terminal SSB depletion is NA: no unfished reference was derived for\n",
        "  this run, so there is nothing to divide by. The dynamic column is\n",
        "  unaffected.\n", sep = "")
  }
  cat("  the blocks have different shapes on purpose; as.data.frame() on the\n",
      "  whole object will not work -- take one block.\n", sep = "")
  invisible(x)
}


#' Check which saved MSE simulation files can be loaded
#'
#' @param dir Directory used to save files from \code{\link{run_mse}}
#' @param file file name used to save files from \code{\link{run_mse}}
#'
#' @return list of MSE simulations/run
#' @export
#'
check_mse <- function(dir = NULL, file = NULL){

  # - Get file names
  mse_files <- list.files(path = dir, pattern = paste0(file, "EMs_from_OM_Sim_"))
  mse_order <- as.numeric(gsub(".rds", "", sapply(strsplit(mse_files, "EMs_from_OM_Sim_"), "[[", 2)))
  mse_files <- mse_files[order(mse_order)]

  ### Set up parallel processing (PSOCK cluster -- cross-platform,
  ### avoids the foreach::%dopar% rlang deparser interaction).
  chk <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  cran_cap <- nzchar(chk) && !chk %in% c("false", "0", "no")
  cores <- if (cran_cap) 2L else max(1L, parallel::detectCores() - 6L)
  use_parallel <- length(mse_files) > 1L && cores > 1L

  read_one <- function(i) {
    mse_tmp <- readRDS(file = paste0(dir, "/", mse_files[i]))
    check_tmp <- data.frame(Dir = dir, File = mse_files[i], Use = NA, Ind = i)
    if (!is.null(mse_tmp$use_sim)) {
      check_tmp$Use <- mse_tmp$use_sim
    }
    check_tmp
  }

  if (use_parallel) {
    cl <- parallel::makeCluster(min(cores, length(mse_files)))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl, c("dir", "mse_files"), envir = environment())
    rows <- parallel::parLapply(cl, seq_along(mse_files), read_one)
  } else {
    rows <- lapply(seq_along(mse_files), read_one)
  }
  check_df <- do.call(rbind, rows)

  return(check_df)
}




#' Load saved MSE simulation runs
#'
#' @param dir Directory used to save files from \code{\link{run_mse}}
#' @param file file name used to save files from \code{\link{run_mse}}
#' @param exclude index of MSE simulations not to load
#' @param include_em whether the EMs should be loaded or not (default = TRUE)
#'
#' @return list of MSE simulations/run
#' @export
#'
load_mse <- function(dir = NULL, file = NULL, exclude = NULL, include_em = TRUE){

  # - Get file names
  mse_files <- list.files(path = dir, pattern = paste0(file, "EMs_from_OM_Sim_"))
  mse_order <- as.numeric(gsub(".rds", "", sapply(strsplit(mse_files, "EMs_from_OM_Sim_"), "[[", 2)))
  mse_files <- mse_files[order(mse_order)]


  if(!is.null(exclude)){
    mse_files <- mse_files[-exclude]
  }

  ### Set up parallel processing
  # cores = (detectCores()/2)-1
  # doParallel::registerDoParallel(cores)
  mse_tmp <- list()
  for(i in 1:length(mse_files)){
    # mse <- foreach(i = 1:length(mse_files),
    #                .combine = "c") %dopar% {
    # mse_tmp <- list(readRDS(file = paste0(dir,"/", mse_files[i])))
    mse_tmp[[i]] <- (readRDS(file = paste0(dir,"/", mse_files[i])))

    # - Remove EM
    if(!include_em){
      mse_tmp[[i]]$EM <- NULL
    }

    # - Clean EM data
    if(length(mse_tmp[[i]]$EM) > 1){
      for(em in 2:length(mse_tmp[[i]]$EM)){
        mse_tmp[[i]]$EM[[em]]$data_list$weight <- NULL
        mse_tmp[[i]]$EM[[em]]$data_list$emp_sel <- NULL
        mse_tmp[[i]]$EM[[em]]$data_list$age_trans_matrix <- NULL
        mse_tmp[[i]]$EM[[em]]$data_list$age_error <- NULL
        mse_tmp[[i]]$EM[[em]]$data_list$NByageFixed <- NULL
        mse_tmp[[i]]$EM[[em]]$data_list$aLW <- NULL
        mse_tmp[[i]]$EM[[em]]$data_list$diet_data <- NULL
        mse_tmp[[i]]$EM[[em]]$data_list$ration_data <- NULL
        mse_tmp[[i]]$EM[[em]]$data_list$aLW <- NULL
        mse_tmp[[i]]$EM[[em]]$estimated_params <- NULL
      }
    }


    # Only use these bits for OM
    mse_tmp[[i]]$OM$initial_params <- NULL
    mse_tmp[[i]]$OM$bounds <- NULL
    mse_tmp[[i]]$OM$map <- NULL
    mse_tmp[[i]]$OM$obj <- NULL
    mse_tmp[[i]]$OM$opt <- NULL
    mse_tmp[[i]]$OM$sdrep <- NULL
    mse_tmp[[i]]$OM$quantities[!names(mse_tmp[[i]]$OM$quantities) %in% .mse_keep_quantities] <- NULL


    # Only use these bits for OM no F.
    # Guarded: run_mse() leaves OM_no_F present and NULL when the unfished run
    # failed, and `NULL$field <- NULL` DELETES the element -- which would undo
    # that, so a simulation read back off disk would lose the marker that says
    # its reference run was attempted and failed.
    if (!is.null(mse_tmp[[i]]$OM_no_F)) {
      mse_tmp[[i]]$OM_no_F$initial_params <- NULL
      mse_tmp[[i]]$OM_no_F$bounds <- NULL
      mse_tmp[[i]]$OM_no_F$map <- NULL
      mse_tmp[[i]]$OM_no_F$obj <- NULL
      mse_tmp[[i]]$OM_no_F$opt <- NULL
      mse_tmp[[i]]$OM_no_F$sdrep <- NULL
      mse_tmp[[i]]$OM_no_F$quantities[!names(mse_tmp[[i]]$OM_no_F$quantities) %in% .mse_keep_quantities] <- NULL
    }

    # - Return
    mse_tmp[[i]]$name <- mse_files[i]
    #  mse_tmp
  }
  # mse <- lapply(mse_files, function(x) readRDS(file = paste0(dir,"/", x)))
  closeAllConnections()

  names(mse_tmp) <- paste0("Sim_", 1:length(mse_tmp))
  return(mse_tmp)
}

