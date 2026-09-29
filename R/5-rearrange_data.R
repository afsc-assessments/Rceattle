#' Flag the fleet that holds each shared parameter block's penalty
#'
#' Fleets sharing an index estimate one block of parameters, so its prior /
#' penalty must be accumulated once. Returns 1 for the first estimated fleet in
#' each group and 0 for the rest. An "Off" fleet estimates nothing, so it is
#' never chosen while an estimated fleet is available, the same donor rule
#' `adjust_map_shared_params()` uses. Fleets with no index (`NA`) share with
#' nobody and always lead.
#'
#' @keywords internal
#' @noRd
.group_lead <- function(key, is_off) {
  lead <- rep(0L, length(key))
  for (k in unique(key[!is.na(key)])) {
    rows <- which(!is.na(key) & key == k)
    est  <- rows[!is_off[rows]]
    lead[if (length(est)) est[1] else rows[1]] <- 1L
  }
  lead[is.na(key)] <- 1L
  as.integer(lead)
}

# Pull a `fleet_control` column as an integer vector. `.pull_int0()` additionally
# shifts a 1-based R index to the 0-based index the C++ template expects -- the
# single place that conversion happens, so an off-by-one can only be introduced
# (or fixed) once. `.pull_int0()` keeps `- 1` (double), matching the original
# `pull() %>% as.integer() - 1` result type exactly.
.pull_int  <- function(fc, col) as.integer(fc[[col]])
.pull_int0 <- function(fc, col) as.integer(fc[[col]]) - 1

# Scale each observation to proportions, as the composition and CAAL
# likelihoods expect. An all-zero or all-NA row divides by zero and is set to 0,
# i.e. no information from that row. Divides by the row sums rather than going
# through `t(apply())`, which returns a transposed 1 x nrow matrix when the
# input has a single bin column. The zero-row guard is load-bearing: a matrix
# with no rows divides by a zero-length vector and loses its dimensions.
.normalise_rows <- function(m) {
  if (nrow(m) == 0) return(m)
  storage.mode(m) <- "numeric"
  m <- m / rowSums(m, na.rm = TRUE)
  m[is.infinite(m)] <- 0
  m[is.na(m)] <- 0
  colnames(m) <- NULL   # as before: the bin columns reach TMB unnamed
  m
}

#' Effective selectivity start year, resolved across mirrored fleets
#'
#' Fleets sharing a `Selectivity_index` share one set of selectivity deviations,
#' so the shared block must start at the earliest year any of them has data.
#' Returns the per-fleet group minimum of `Sel_start_year` (NA = `styr`), keeping
#' the mask and the penalty anchor independent of `fleet_control` row order.
#'
#' @keywords internal
#' @noRd
.sel_start_year_by_group <- function(fleet_control, styr) {
  n   <- nrow(fleet_control)
  styr <- as.integer(styr)
  ssy <- fleet_control$Sel_start_year
  if (is.null(ssy)) return(rep(styr, n))

  eff <- suppressWarnings(as.integer(ssy))
  eff[is.na(eff)] <- styr                       # NA -> styr (no pre-start masking)

  grp <- fleet_control$Selectivity_index
  if (is.null(grp)) return(eff)
  # A fleet with no Selectivity_index shares with nobody -> its own group.
  grp <- ifelse(is.na(grp), paste0("_row", seq_len(n)), as.character(grp))

  mins <- tapply(eff, grp, min)
  as.integer(mins[match(grp, names(mins))])
}

#' Rearrange a data_list for TMB
#'
#' @description Function to rearrange a \code{data_list} object to be read into TMB
#'
#' @param data_list an Rceattle data_list
#' @param build_osa Logical. Passed to [build_osa_data()]; when `TRUE` the full
#'   one-step-ahead (OSA) observation data is assembled so [osa_residuals()] can
#'   be computed. Default `FALSE` (the fast path used by simulation testing). The
#'   composition proportion offset is read from `data_list$comp_offset` (default
#'   `1e-5`).
#'
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom dplyr n
#' @importFrom tidyselect contains
rearrange_data <- function(data_list, build_osa = FALSE){

  # Upgrade any deprecated fleet_control column names to canonical, so the
  # catchability/selectivity reads below see the canonical names even when
  # rearrange_data() is called directly (build_params()/build_map() self-upgrade
  # via switch_check(); this closes the same gap here). No-op in the fit_mod()
  # pipeline, where switch_check() has already run.
  data_list$fleet_control <-
    .rce_upgrade_fleet_control_aliases(data_list$fleet_control)

  # Fail clearly on a malformed fleet_control rather than via a cryptic dplyr
  # error deep in the reshaping below. These identity / structural columns are
  # read by name (and by `$`) throughout; a missing one has no sensible default.
  if (is.null(data_list$fleet_control) || nrow(data_list$fleet_control) == 0)
    stop("rearrange_data(): 'fleet_control' is missing or empty.", call. = FALSE)
  .required_fc <- c("Species", "Fleet_code", "Fleet_type", "Fleet_name",
                    "Selectivity", "Selectivity_index", "Sel_start_year",
                    "Catchability_init")
  .missing_fc <- setdiff(.required_fc, names(data_list$fleet_control))
  if (length(.missing_fc) > 0)
    stop("rearrange_data(): 'fleet_control' is missing required column(s): ",
         paste(.missing_fc, collapse = ", "),
         ". Run switch_check() first to fill defaulted columns.", call. = FALSE)

  # Both read the column as supplied: convert_switches() turns a blank and a type
  # it cannot map into the same NA, and an out-of-range code into itself.
  .rce_stop_blank_fleet_type(data_list$fleet_control)
  # convert_switches() resolves factor columns; do it here too so the refusal
  # below reads the label the user wrote rather than a level index.
  data_list$fleet_control <-
    .rce_defactor_fleet_control(data_list$fleet_control)
  .rce_stop_unreadable_fleet_type(data_list$fleet_control)

  # Convert text to integer for switches used in TMB
  data_list <- convert_switches(data_list)

  # Data dimensions
  max_sex <- max(data_list$nsex, na.rm = TRUE)
  max_age <- max(data_list$nages, na.rm = TRUE)
  max_length <- max(data_list$nlengths, na.rm = TRUE)
  yrs_hind <- data_list$styr:data_list$endyr
  yrs_proj <- data_list$styr:data_list$projyr
  nyrs_hind <- length(yrs_hind)
  nyrs_proj <- length(yrs_proj)

  # 1 - Fleet control ----
  # - 0) Vector to save  species
  data_list$flt_spp <- .pull_int0(data_list$fleet_control, "Species")

  # - 1) Fleet type; 0 = don't fit, 1 = fishery, 2 = survey
  data_list$flt_type <- .pull_int(data_list$fleet_control, "Fleet_type")

  # - 3) Month of observation
  data_list$flt_month <- data_list$fleet_control %>%
    dplyr::pull(.data$Month)

  # - 4) Selectivity type
  data_list$flt_sel_type <- .pull_int(data_list$fleet_control, "Selectivity")

  # - 4b) Selectivity penalty "lead": 1 for the one fleet that carries each shared
  #       block's penalty, 0 for the fleets mirroring it, so a mirrored selectivity
  #       is penalized once (matching ADMB). Fleets are grouped by Selectivity_index
  #       AND selectivity type, so a fleet with empirical selectivity (type 0)
  #       sharing an index with a parametric fleet does not suppress the parametric
  #       fleet's penalty. The lead is the first ESTIMATED fleet in the group: an
  #       "Off" fleet estimates nothing and the cpp penalty is gated on flt_type,
  #       so leading with one would drop the group's penalty entirely. This matches
  #       the map donor chosen in adjust_map_shared_params().
  {
    .off <- data_list$flt_type == 0
    data_list$flt_sel_lead <- .group_lead(
      paste(data_list$fleet_control$Selectivity_index, data_list$flt_sel_type), .off)
  }

  data_list$flt_sel_dim <- data_list$fleet_control %>%
    dplyr::mutate(Selectivity_dimension = dplyr::case_when(
      .data$Selectivity_dimension == "Age" ~ 0,
      .data$Selectivity_dimension == "Length" ~ 1,
      .default = NA
    )) %>%
    dplyr::pull(.data$Selectivity_dimension) %>% as.integer()

  # - 4b) Per-fleet offset converting a selectivity bin column to the 0-based C++
  #       index. Age-based fleets give absolute ages -> offset is the species'
  #       minage; length-based fleets give 1-based length-bin ordinals -> offset 1.
  sel_bin_offset <- ifelse(!is.na(data_list$flt_sel_dim) & data_list$flt_sel_dim == 1L,
                           1L,
                           data_list$minage[data_list$fleet_control$Species])

  # - 5) Number of ages/lengths for non-parametric selectivity
  data_list$flt_n_sel_bins <- data_list$fleet_control %>%
    dplyr::mutate(N_sel_bins = ifelse(is.na(.data$N_sel_bins), -999, .data$N_sel_bins)) %>%
    dplyr::pull(.data$N_sel_bins) %>% as.integer()

  # - 6) Time-varying selectivity type.
  data_list$flt_varying_sel <- .pull_int(data_list$fleet_control, "Time_varying_sel")

  # - 7) First age selected
  data_list$bin_first_selected <- data_list$fleet_control %>%
    dplyr::mutate(Bin_first_selected = .data$Bin_first_selected - 1, # R to C++ indexing
                  Bin_first_selected = ifelse(is.na(.data$Bin_first_selected), 0, .data$Bin_first_selected)) %>%
    dplyr::pull(.data$Bin_first_selected) %>% as.integer()

  # - 8) Age of max selectivity (used for normalization). If NA, does not normalize
  # First bin each fleet may normalize at, on the column's own scale: below it
  # the curve is zeroed, so a value there means "normalize by the maximum".
  .norm_bounds <- lapply(seq_len(nrow(data_list$fleet_control)), function(i) {
    age <- isTRUE(data_list$fleet_control$Selectivity_dimension[i] == "Age")
    nb  <- if (age) data_list$nages[data_list$fleet_control$Species[i]]
           else data_list$nlengths[data_list$fleet_control$Species[i]]
    .rce_sel_norm_bounds(data_list, i, age, nb)
  })
  .norm_lo <- vapply(.norm_bounds, function(b) as.numeric(b$lo), numeric(1))
  .norm_hi <- vapply(.norm_bounds, function(b) as.numeric(b$hi), numeric(1))
  # A Fixed curve is read from emp_sel_obs and never normalized, so its value is
  # allowed to be stale -- the same fleets switch_check() leaves alone.
  .norm_hi[data_list$flt_sel_type == 0] <- NA_real_

  .is_pm <- data_list$fleet_control$Selectivity %in% c(11, "LogisticPM")

  # Resolve the word first, so only a real bin is shifted to the 0-based index:
  # -99 normalizes by the maximum and -999 does not normalize. An unreadable
  # value is refused here too -- switch_check() does not run on every path in.
  .norm1 <- .rce_sel_norm_code(data_list$fleet_control$Sel_norm_bin,
                               lo = .norm_lo, allow_all = .is_pm)
  .rce_stop_bad_sel_norm(.norm1, "Sel_norm_bin", data_list$fleet_control, .norm_hi)
  data_list$sel_norm_bin1 <- as.integer(dplyr::case_when(
    is.na(.norm1) ~ -999,                        # do not normalize
    .norm1 <  0   ~ -99,                         # normalize by the maximum
    TRUE          ~ .norm1 - sel_bin_offset))

  # - 8b) Whether normalization pools its reference across sexes (WithinSex = 0,
  #       AcrossSexes = 1); already an integer code by this point.
  data_list$sel_norm_scope <- data_list$fleet_control %>%
    dplyr::pull(.data$Sel_norm_scope) %>% as.integer()

  # - 9) upper age of max selectivity (used for normalization). If NA, does not normalize
  # As above, but the upper bin keeps its own convention: selectivity.hpp reads
  # anything below 0 as "no age range", so a negative is passed through.
  .norm2 <- .rce_sel_norm_code(data_list$fleet_control$Sel_norm_bin_upper,
                               lo = .norm_lo, allow_all = .is_pm)
  .rce_stop_bad_sel_norm(.norm2, "Sel_norm_bin_upper", data_list$fleet_control, .norm_hi)
  data_list$sel_norm_bin2 <- as.integer(dplyr::case_when(
    is.na(.norm2) ~ -999,
    .norm2 <  0   ~ .norm2,
    TRUE          ~ .norm2 - sel_bin_offset))

  # - 9b) Per-fleet selectivity start year (0-based from styr). Selectivity
  #       penalties begin the year after this (excludes pre-survey years + the
  #       start-year boundary). NA -> styr (index 0). Used by LogisticPM (type 11).
  #       Resolved to the minimum across each Selectivity_index group, since
  #       mirrored fleets share one deviation block (see .sel_start_year_by_group).
  data_list$flt_sel_start_yr <- as.integer(
    .sel_start_year_by_group(data_list$fleet_control, data_list$styr) - data_list$styr)

  # - 9c) First bin (0-based) for the non-parametric shape penalty. NA -> -999
  #       (cpp falls back to bin_first_selected). Lets the ascending/descending
  #       shape constraint span a narrower range than the first selected bin
  #       (e.g. ATS: age-1 is selected but the shape penalty starts at mina_ats).
  data_list$flt_sel_pen_first_bin <- data_list$fleet_control %>%
    dplyr::mutate(Sel_pen_first_bin = .data$Sel_pen_first_bin - sel_bin_offset,
                  Sel_pen_first_bin = ifelse(is.na(.data$Sel_pen_first_bin), -999, .data$Sel_pen_first_bin)) %>%
    dplyr::pull(.data$Sel_pen_first_bin) %>% as.integer()

  # - 9d) Last (left) bin (0-based) of the shape-penalty adjacent pairs. NA -> -999
  #       (cpp falls back to nbins-2). With 9c, bounds the shape penalty to a sub-
  #       range (e.g. RTMB "rpm" fishery smoothness over ages 6-11, ATS 5-7).
  data_list$flt_sel_pen_last_bin <- data_list$fleet_control %>%
    dplyr::mutate(Sel_pen_last_bin = .data$Sel_pen_last_bin - sel_bin_offset,
                  Sel_pen_last_bin = ifelse(is.na(.data$Sel_pen_last_bin), -999, .data$Sel_pen_last_bin)) %>%
    dplyr::pull(.data$Sel_pen_last_bin) %>% as.integer()

  # - 9e) Non-parametric shape-penalty mode: "Directional" (0; sign of Sel_curve_pen1
  #       -> one-sided decreasing/increasing, ADMB/AMAK) or "Smooth" (1; two-sided
  #       d^2 over adjacent ages, RTMB "rpm"). NA/other -> 0.
  data_list$flt_sel_shape_mode <- data_list$fleet_control %>%
    dplyr::mutate(Sel_shape_mode = dplyr::case_when(
      .data$Sel_shape_mode %in% c(1, "Smooth", "smooth") ~ 1L,
      .default = 0L)) %>%
    dplyr::pull(.data$Sel_shape_mode) %>% as.integer()

  # - 9g) AMAK "avgsel" base-level penalty weight (type 9). ADMB adds
  #       10*square(avgsel_*), avgsel = log(mean(exp(base coffs))); 0 = off (default).
  data_list$flt_sel_avgsel_pen <- suppressWarnings(
    as.numeric(data_list$fleet_control$Sel_avgsel_pen))
  data_list$flt_sel_avgsel_pen[is.na(data_list$flt_sel_avgsel_pen)] <- 0

  # - 9f) NonParametricPM (type 9) bin cap (0-based): the realized selectivity is
  #       held flat at/after this bin (RTMB cap_old_age). NA -> -999 (no cap).
  data_list$flt_sel_cap_bin <- data_list$fleet_control %>%
    dplyr::mutate(Sel_cap_bin = .data$Sel_cap_bin - sel_bin_offset,
                  Sel_cap_bin = ifelse(is.na(.data$Sel_cap_bin), -999, .data$Sel_cap_bin)) %>%
    dplyr::pull(.data$Sel_cap_bin) %>% as.integer()

  # - 10) Index indicating whether to do dirichlet multinomial or a multinomial
  data_list$comp_ll_type <- .pull_int(data_list$fleet_control, "Comp_distribution")
  data_list$caal_ll_type <- .pull_int(data_list$fleet_control, "CAAL_distribution")
  data_list$diet_ll_type <- data_list$Diet_distribution %>%
    as.integer()

  # - 10b) Composition young/old-tail accumulation bins (AFSC ac_yng / ac_old).
  # Absent columns or NAs default to NO accumulation (young = 1; old = 0, which
  # the cpp reads as "use the last bin"), so existing models are bit-identical.
  .accum_col <- function(nm, default) {
    v <- if (nm %in% names(data_list$fleet_control)) {
      as.integer(dplyr::pull(data_list$fleet_control, nm))
    } else {
      rep(NA_integer_, nrow(data_list$fleet_control))
    }
    v[is.na(v)] <- default
    v
  }
  data_list$comp_accum_young <- .accum_col("Comp_accum_young", 1L)
  data_list$comp_accum_old   <- .accum_col("Comp_accum_old", 0L)

  # - 11) Index units (1 = weight, 2 = numbers)
  data_list$flt_units <- .pull_int(data_list$fleet_control, "Observation_units")

  # - 12) Dim1 of weight (what weight-at-age data set)
  data_list$flt_wt_index <- .pull_int0(data_list$fleet_control, "Weight_index")

  # - 13) Dim3 of age transition matrix (what ALK to use)
  data_list$flt_age_transition_index <- .pull_int0(data_list$fleet_control, "Age_transition_index")

  # - 13b) Ageing error matrix per fleet. Absent or NA means the fleet's own
  #        species, which reproduces the one-matrix-per-species behaviour
  #        exactly, so an existing fleet_control needs no new column.
  {
    .ae <- data_list$fleet_control[["Ageing_error_index"]]
    .sp <- as.integer(data_list$fleet_control$Species)
    .ae <- if (is.null(.ae)) .sp else {
      v <- suppressWarnings(as.integer(.ae)); ifelse(is.na(v), .sp, v)
    }
    data_list$flt_ageing_error_index <- as.integer(.ae) - 1L   # 0-based for C++
  }

  # - 14) Parametric form of q
  data_list$est_index_q <- .pull_int(data_list$fleet_control, "Catchability")

  # - 15) Time varying q type
  #
  # `Time_varying_q` is overloaded: for most fleets it holds a tv_q_map mode
  # code (0 Off / 1 IID / 2 AR1 / 4 RandomWalk), but for "AR1" and
  # "Environmental" catchability it holds `env_data` column indices, which
  # `convert_switches()` leaves unconverted. "Environmental" (5) may carry a
  # comma-separated list; those indices reach the model via the `index_q_beta`
  # map, so emit 0 ("not applicable") here -- it collides with no mode code
  # and is never read for Catchability = 5.
  .tv_q <- data_list$fleet_control %>% dplyr::pull(.data$Time_varying_q)
  .tv_q[data_list$est_index_q %in% 5L] <- "0"
  data_list$index_varying_q <- as.integer(.tv_q)

  # - 15b) Catchability "lead", the q analogue of flt_sel_lead. Fleets sharing a
  #        Catchability_index estimate ONE catchability (and one deviate vector),
  #        so the q prior and the deviate penalties are accumulated on the lead
  #        fleet only. Without this they were applied once per sharing fleet to
  #        the same parameter, e.g. counting a q prior twice for a mirrored pair.
  data_list$flt_q_lead <- .group_lead(data_list$fleet_control$Catchability_index,
                                      data_list$flt_type == 0)

  # - 16) Whether to estimate standard deviation of index time series
  data_list$est_sigma_index <- .pull_int(data_list$fleet_control, "Estimate_index_sd")

  # - 17) Whether to estimate standard deviation of fishery time series
  data_list$est_sigma_fsh <- .pull_int(data_list$fleet_control, "Estimate_catch_sd")

  # - 18) Survey/index biomass likelihood family (0 = lognormal IID, 1 = MVN covariance)
  data_list$index_ll_type <- .pull_int(data_list$fleet_control, "Index_distribution")

  # NOTE: index timing carries no data element. The cpp reads it off `flt_type`
  # (section 8.1): a survey's index is a snapshot at its observation month, a
  # fishery's is the year-average N*(1-exp(-Z))/Z that CPUE is proportional to.

  data_list$index_log_q_prior <- log(data_list$fleet_control$Catchability_init)

  # Species names
  data_list$spnames <- NULL

  # * F reference points ----
  data_list$Ftarget_percent <- data_list$Ftarget
  data_list$Flimit_percent <- data_list$Flimit


  # 2 -  Index data ----
  # - Seperate index metadata from observation
  data_list$index_ctl <- data_list$index_data %>%
    dplyr::select(Fleet_code, Species, Year) %>%
    dplyr::mutate_all(as.integer) %>%
    as.matrix()

  data_list$index_n <- as.matrix(data_list$index_data[,c("Month")])

  data_list$index_obs <- data_list$index_data %>%
    dplyr::select(Observation, Log_sd) %>%
    dplyr::mutate_all(as.numeric) %>%
    as.matrix()

  # - Survey-index covariance matrices (Index_distribution == "MVN" or "MVNORM").
  #   TMB evaluates TMB's density::MVNORM(Sigma) on each covariance fleet's fitted
  #   residual vector r = obs - q*pred, i.e. 0.5*(r' Sigma^-1 r + logdet(Sigma) +
  #   n*log(2*pi)). "MVN" reports the bare quadratic form 0.5 r' Sigma^-1 r (the
  #   AMAK/ebswp DoCovBTS value) by subtracting the fixed normalizing constant
  #   index_cov_const = 0.5*(logdet(Sigma) + n*log(2*pi)); "MVNORM" keeps the full
  #   density. We pass Sigma directly (MVNORM factorizes it internally -- more
  #   robust than an explicit inverse) plus the precomputed constant. Non-covariance
  #   fleets get a 1x1 inert dummy (const 0) so the list is length n_flt and TMB can
  #   index it by fleet code. Sigma must be square/symmetric with dimension equal to
  #   the number of fitted survey obs for the fleet, in index_data row order.
  n_flt_cov <- nrow(data_list$fleet_control)
  data_list$index_cov_mat   <- vector("list", n_flt_cov)
  data_list$index_cov_const <- numeric(n_flt_cov)
  for(flt in seq_len(n_flt_cov)){
    if(isTRUE(data_list$index_ll_type[flt] %in% c(1L, 2L))){
      flt_name <- data_list$fleet_control$Fleet_name[flt]
      flt_code <- data_list$fleet_control$Fleet_code[flt]
      Sigma <- data_list$index_cov[[as.character(flt_name)]]
      if(is.null(Sigma)) Sigma <- data_list$index_cov[[as.character(flt_code)]]
      Sigma <- as.matrix(Sigma)
      n_fit <- sum(data_list$index_data$Fleet_code == flt_code &
                     data_list$index_data$Year > 0 &
                     data_list$index_data$Year <= data_list$endyr &
                     data_list$index_data$Observation > 0)
      if(nrow(Sigma) != n_fit || ncol(Sigma) != n_fit){
        stop(sprintf(paste0("Index covariance matrix for fleet '%s' is %d x %d but the fleet has %d fitted ",
                            "survey observations. Provide a square Sigma matching the fitted survey years ",
                            "(Year in [styr, endyr], Observation > 0), in index_data row order."),
                     flt_name, nrow(Sigma), ncol(Sigma), n_fit))
      }
      Sigma <- (Sigma + t(Sigma)) / 2  # symmetrize away tiny numerical asymmetry from the source file
      logdet <- tryCatch(as.numeric(determinant(Sigma, logarithm = TRUE)$modulus),
                         error = function(e) stop(sprintf(
                           "Index covariance matrix for fleet '%s' is not invertible: %s", flt_name, conditionMessage(e))))
      data_list$index_cov_mat[[flt]]  <- Sigma
      data_list$index_cov_const[flt]  <- 0.5 * (logdet + nrow(Sigma) * log(2 * pi))
    } else {
      data_list$index_cov_mat[[flt]]  <- matrix(0, 1, 1)  # inert dummy (present-but-off)
      data_list$index_cov_const[flt]  <- 0
    }
  }


  # 3 -  Catch data ----
  # - The initial equilibrium catch is predicted from the equilibrium age
  #   structure rather than from a hindcast year, so it gets its own inputs.
  #   Only a mode that estimates Finit reads it; under any other the rows are
  #   catch history and are dropped.
  .equil_rows <- .rce_equil_catch_candidates(data_list)
  if (!is.null(.equil_rows) && nrow(.equil_rows)) {
    .equil_rows <- .equil_rows[.rce_equil_catch_rows(
      .equil_rows, data_list$styr, data_list$initMode), , drop = FALSE]
  }
  if (is.null(.equil_rows) || !nrow(.equil_rows)) {
    # TMB needs a matrix with the right column count even when there are no rows.
    data_list$equil_catch_ctl <- matrix(0L, 0, 2)
    data_list$equil_catch_obs <- matrix(0,  0, 2)
  } else {
    data_list$equil_catch_ctl <- .equil_rows %>%
      dplyr::select(Fleet_code, Species) %>%
      dplyr::mutate_all(as.integer) %>%
      as.matrix()
    data_list$equil_catch_obs <- .equil_rows %>%
      dplyr::select(Catch, Log_sd) %>%
      dplyr::mutate_all(as.numeric) %>%
      as.matrix()
  }
  # catch_data holds the fitted catch only. catch_obs is built from this frame,
  # and row names here would travel into obsvec's OSA bookkeeping.
  data_list$catch_data <- data_list$catch_data[
    is.na(data_list$catch_data$Year) |
      data_list$catch_data$Year != (data_list$styr - 1L), , drop = FALSE]
  rownames(data_list$catch_data) <- NULL

  # - Seperate catch metadata from observation
  data_list$catch_ctl <- data_list$catch_data %>%
    dplyr::select(Fleet_code, Species, Year) %>%
    dplyr::mutate_all(as.integer) %>%
    as.matrix()

  data_list$catch_n <- as.matrix(data_list$catch_data[,c("Month")])

  data_list$catch_obs <- data_list$catch_data %>%
    dplyr::select(Catch, Log_sd) %>%
    dplyr::mutate_all(as.numeric) %>%
    as.matrix()


  # 4 -  Comp data ----
  # - Seperate comp metadata from observation
  data_list$comp_ctl <- data_list$comp_data %>%
    dplyr::select(Fleet_code, Species, Sex, Age0_Length1, Year) %>%
    dplyr::mutate_all(as.integer) %>%
    as.matrix()

  data_list$comp_n <- data_list$comp_data %>%
    dplyr::select(Month, Sample_size) %>%
    dplyr::mutate_all(as.numeric) %>%
    as.matrix()

  data_list$comp_obs <- data_list$comp_data %>%
    dplyr::select(contains("Comp_")) %>%
    dplyr::mutate_all(as.numeric) %>%
    as.matrix()

  data_list$comp_obs <- .normalise_rows(data_list$comp_obs)
  data_list <- check_composition_data(data_list)

  # 5 - CAAL data ----
  data_list$caal_ctl <- data_list$caal_data %>%
    dplyr::mutate(Length_bin = factor(.data$Length)) %>%
    dplyr::select(Fleet_code, Species, Sex, Year, Length_bin) %>%
    dplyr::mutate_all(as.integer) %>%
    as.matrix()

  data_list$caal_n <- data_list$caal_data %>%
    dplyr::select(Sample_size) %>%
    dplyr::mutate_all(as.numeric) %>%
    as.matrix()

  data_list$caal_obs <- data_list$caal_data %>%
    dplyr::select(contains("CAAL_")) %>%
    dplyr::mutate_all(as.numeric) %>%
    as.matrix()

  data_list$caal_obs <- .normalise_rows(data_list$caal_obs)
  data_list <- check_caal_data(data_list)

  # * Get lengths from CAAL ----
  data_list$lengths <- matrix(rep(1:max(data_list$nlengths), data_list$nspp), nrow = data_list$nspp, byrow = TRUE)

  if(nrow(data_list$caal_data) > 0){
    caal_lengths <- data_list$caal_data %>%
      dplyr::distinct(.data$Species, .data$Length) %>%
      dplyr::arrange(.data$Species, .data$Length) %>%
      dplyr::group_by(.data$Species) %>%
      dplyr::mutate(Bin = paste0("Bin", 1:dplyr::n())) %>%
      tidyr::pivot_wider(names_from = "Bin", values_from = "Length")

    data_list$lengths[caal_lengths$Species, 1:ncol(caal_lengths[,-1])] <- as.matrix(caal_lengths[,-1])

    # - Check to make sure n-bins in data match nlengths
    n_bins_data <- data_list$caal_data %>%
      dplyr::distinct(.data$Species, .data$Length) %>%
      dplyr::count(.data$Species)

    if(any(data_list$nlengths[n_bins_data$Species] != n_bins_data$n)){
      stop("Number of length bins in CAAL data does not match nlengths in control file.
           If some lengths are missing from CAAL data, please add them to the data
           as rows of all 1s with Sample_size = 0")
    }
  }

  # * DoubleNormalSS3 ends ----
  # Both ends scaled unless fit_mod() finds SS3's -999 in the starting values.
  data_list$sel_dn6_ends <- matrix(1L, nrow(data_list$fleet_control), 2)

  # * Population length bins ----
  # The age-length key, weight-at-length and maturity-at-length are computed on
  # these finer bins (SS3's population length bins) and summed into the data
  # bins above. A species without its own grid uses its data bins.
  data_list <- .rce_pop_length_bins(data_list)

  # A length-based curve is built on the population bins, but the columns that
  # name a bin are written and checked as DATA bin ordinals, which is what the
  # user has. Translate them onto the population grid: the first selected bin
  # and a normalization reference take that data bin's first population bin, an
  # upper range bound its last, so each still covers the lengths it names. A
  # species with no separate grid is unchanged, every run being one bin long.
  .len_flt <- which(as.character(data_list$fleet_control$Selectivity_dimension) == "Length")
  for (i in .len_flt) {
    sp <- as.integer(data_list$fleet_control$Species[i])
    if (is.na(sp) || data_list$nlengths_pop[sp] == data_list$nlengths[sp]) next
    map <- data_list$pop_to_data_bin[sp, seq_len(data_list$nlengths_pop[sp])]
    to_pop <- function(b, last = FALSE) {
      if (is.na(b) || b < 0) return(b)             # -99 / -999 are not bins
      run <- which(map == b)
      if (!length(run)) return(b)
      as.integer(if (last) run[length(run)] else run[1]) - 1L
    }
    data_list$bin_first_selected[i] <- to_pop(data_list$bin_first_selected[i])
    data_list$sel_norm_bin1[i]      <- to_pop(data_list$sel_norm_bin1[i])
    data_list$sel_norm_bin2[i]      <- to_pop(data_list$sel_norm_bin2[i], last = TRUE)
  }

  # * Maturity-at-length ----
  # Logistic in length (cm); a species with no L50 / slope keeps the age-based
  # maturity sheet.
  L50   <- as.numeric(data_list$L50_mat_len   %||% rep(NA_real_, data_list$nspp))
  slope <- as.numeric(data_list$slope_mat_len %||% rep(NA_real_, data_list$nspp))
  data_list$mat_len_use  <- as.integer(!is.na(L50) & !is.na(slope))
  data_list$mat_len_pars <- cbind(ifelse(is.na(L50), 0, L50), ifelse(is.na(slope), 0, slope))


  # 6 -  Diet data ----
  # - Seperate diet metadata from observation
  if(!is.null(data_list$diet_data)){ # Add a check in case there's no diet data

    # Create a temporary diet data frame to work with
    diet_dat <- data_list$diet_data

    # Add n_stomach_obs and the stomach_id vector to the main data_list
    if(length(diet_dat$stomach_id) == 0){
      data_list$n_stomach_obs = 0
    } else{
      data_list$n_stomach_obs <- max(diet_dat$stomach_id) + 1
    }
    data_list$stomach_id <- diet_dat$stomach_id

    # Create diet_ctl (without the stomach_id column)
    data_list$diet_ctl <- diet_dat %>%
      dplyr::select(Pred, Prey, Pred_sex, Prey_sex, Pred_age, Prey_age, Year) %>%
      dplyr::mutate_all(as.integer)

    # The observations themselves, in the same row order as diet_ctl
    data_list$diet_obs <- diet_dat %>%
      dplyr::select(Sample_size, Stomach_proportion_by_weight) %>%
      dplyr::mutate_all(as.numeric)
  }


  # 7 -  Seperate empirical selectivity info from observation ----
  yrs <- data_list$styr:data_list$endyr
  if(nrow(data_list$emp_sel) > 0 ){
    for(i in 1:nrow(data_list$emp_sel)){
      # Fill in years
      if(data_list$emp_sel$Year[i] == 0){

        # Change first year
        data_list$emp_sel$Year[i] <- yrs[1]

        # Change the rest
        emp_sel_tmp <- data_list$emp_sel[i,]
        emp_sel_tmp$Year <- yrs[2]
        emp_sel <- emp_sel_tmp

        for(yr in 3:length(yrs)){
          emp_sel_tmp$Year <- yrs[yr]
          emp_sel <- rbind(emp_sel, emp_sel_tmp)
        }

        data_list$emp_sel <- rbind(data_list$emp_sel, emp_sel)
      }
    }
  }

  data_list$emp_sel_ctl <- data_list$emp_sel %>%
    dplyr::select(Fleet_code, Species, Sex, Year) %>%
    dplyr::mutate_all(as.integer)

  data_list$emp_sel_obs <- data_list$emp_sel %>%
    dplyr::select(contains("Comp_")) %>%
    dplyr::mutate_all(as.numeric)


  # 8 - Rearrange age-transition matrix ----
  age_trans_matrix <- data_list$age_trans_matrix
  unique_age_transition <- unique(as.character(age_trans_matrix$Age_transition_index))
  if(sum(!data_list$pop_age_transition_index %in% unique_age_transition) > 0){
    stop("Check population age_transition index, not in age_transition file")
  }
  age_transition <- array(0, dim = c(length(unique_age_transition), max_sex, max_age, max_length))


  for (i in 1:nrow(age_trans_matrix)) {
    age_transition_ind <- as.numeric(as.character(age_trans_matrix$Age_transition_index[i]))
    sex <- as.numeric(as.character(age_trans_matrix$Sex[i]))
    sp <- as.numeric(as.character(age_trans_matrix$Species[i]))
    age <- as.numeric(as.character(age_trans_matrix$Age[i])) - data_list$minage[sp] + 1

    if (age > data_list$nages[sp]) {
      message(paste0("Error: number of ages in age_trans_matrix for species: ", sp, " and index: ", age_transition_ind))
      message(paste0("is greater than the number of ages specified in the control"))
      message(paste0("Please remove or change nages in control"))
      stop()
    }

    # Handle sex == 0 case for 2-sex species
    sex_values <- if (sex == 0) 1:data_list$nsex[sp] else sex

    for(j in 1:length(sex_values)){
      age_transition[age_transition_ind, sex_values[j], age, 1:data_list$nlengths[sp]] <- as.numeric(as.character(age_trans_matrix[i, (1:data_list$nlengths[sp]) + 5]))

      # Normalize
      age_transition[age_transition_ind, sex_values[j], age, 1:data_list$nlengths[sp]] <- age_transition[age_transition_ind, sex_values[j], age, 1:data_list$nlengths[sp]] / sum(age_transition[age_transition_ind, sex_values[j], age, 1:data_list$nlengths[sp]], na.rm = TRUE)
    }
  }
  data_list$age_trans_matrix <- age_transition


  # 9 - Rearrange age_error matrices ----
  # Indexed by Ageing_error_index, not by species, so one species can carry
  # several matrices -- an ageing method that changed part way through a series
  # reads the same otolith differently before and after, and a fleet picks its
  # matrix through fleet_control$Ageing_error_index. Absent the column the index
  # IS the species, which is what the model did when there was one matrix each.
  data_list$age_error <- as.data.frame(data_list$age_error)
  if (is.null(data_list$age_error$Ageing_error_index)) {
    data_list$age_error$Ageing_error_index <- data_list$age_error$Species
  }
  n_ae <- max(as.numeric(as.character(data_list$age_error$Ageing_error_index)),
              data_list$nspp, na.rm = TRUE)
  arm <- array(0, dim = c(n_ae, max_age, max_age))

  # The observed-age columns were read positionally as "everything after the
  # first two". Name the metadata instead, so adding a column cannot silently
  # shift which columns are read as probabilities.
  .ae_obs  <- setdiff(colnames(data_list$age_error), .RCE_AGE_ERROR_META)

  for (i in seq_len(nrow(data_list$age_error))) {
    sp <- as.numeric(as.character(data_list$age_error$Species[i]))
    idx <- as.numeric(as.character(data_list$age_error$Ageing_error_index[i]))
    true_age <- as.numeric(as.character(data_list$age_error$True_age[i])) - data_list$minage[sp] + 1

    if (true_age > data_list$nages[sp]) {
      message(paste0("Error: number of true ages specified in age_error for species: ", sp))
      message(paste0("is greater than the number of ages specified in the control"))
      message(paste0("Please remove or change nages in control"))
      stop()
    }

    arm[idx, true_age, 1:data_list$nages[sp]] <-
      as.numeric(as.character(data_list$age_error[i, .ae_obs[1:data_list$nages[sp]]]))

    # Normalize
    arm[idx, true_age, 1:data_list$nages[sp]] <- arm[idx, true_age, 1:data_list$nages[sp]] / sum(arm[idx, true_age, 1:data_list$nages[sp]], na.rm = TRUE)
  }
  data_list$age_error <- arm


  # 11 - Set up weight-at-age ----
  # - Convert to array
  data_list$weight <- data_list$weight %>%
    dplyr::mutate(
      Wt_index = as.numeric(as.character(.data$Wt_index)),
      Species = as.numeric(as.character(.data$Species)),
      Sex = as.numeric(as.character(.data$Sex)),
      Year = as.numeric(as.character(.data$Year)),
      Year = ifelse(.data$Year == 0,
                    0,
                    Year - data_list$styr + 1)
    )

  # Pre-allocate the array
  unique_wt <- unique(as.numeric(data_list$weight$Wt_index))
  weight <- array(0, dim = c(max(unique_wt), max_sex, max_age, nyrs_hind))

  # Convert weight data to numeric once, outside the loop
  weight_matrix <- data_list$weight[, (1:max_age) + 5]
  weight_matrix <- apply(weight_matrix, 2, function(x) as.numeric(as.character(x)))
  weight_matrix <- matrix(weight_matrix, ncol = max_age)

  # Check for spaces in the weight data
  if (any(grepl("[[:space:]]", weight_matrix))) {
    stop("Space found in weight data")
  }

  # Main processing
  for (i in seq_len(nrow(data_list$weight))) {
    wt_ind <- data_list$weight$Wt_index[i]
    sp <- data_list$weight$Species[i]
    sex <- data_list$weight$Sex[i]
    yr <- data_list$weight$Year[i]

    # Check if the year is within the valid range
    if (yr <= data_list$endyr - data_list$styr + 1) {
      # Handle year == 0 case
      if(yr == 0) {
        yr <- seq_len(nyrs_hind)
      }

      # Handle sex == 0 case for 2-sex species
      sex_values <- if (sex == 0) 1:data_list$nsex[sp] else sex

      # Assign weights to the array
      weight[wt_ind, sex_values, 1:data_list$nages[sp], yr] <- rep(weight_matrix[i, 1:data_list$nages[sp]], each = length(sex_values))
    }
  }

  # Update the data_list with the new weight array
  data_list$weight_obs <- weight


  # 12 - Set up NByageFixed ----
  # - Convert to array
  data_list$NByageFixed <- data_list$NByageFixed %>%
    dplyr::mutate(Species = as.numeric(as.character(.data$Species)),
           Sex = as.numeric(as.character(.data$Sex)),
           Year = as.numeric(as.character(.data$Year)) - data_list$styr + 1)

  NByageFixed <- array(0, dim = c(data_list$nspp, max_sex, max_age, nyrs_proj))

  if(nrow(data_list$NByageFixed) > 0){
    for (i in 1:nrow(data_list$NByageFixed)) {

      sp <- data_list$NByageFixed$Species[i]
      sex <- data_list$NByageFixed$Sex[i]
      yr <- data_list$NByageFixed$Year[i]

      # Handle sex == 0 case for 2-sex species
      sex_values <- if(sex == 0) 1:data_list$nsex[sp] else sex

      for(j in 1:length(sex_values)){
        if(yr > 0){
          NByageFixed[sp, sex_values[j], 1:max_age, yr] <- as.numeric(data_list$NByageFixed[i,c(1:max_age)+4])
        }
      }
    }
  }

  data_list$NByageFixed <- NByageFixed


  # 13 - Set up environmental indices ----
  # - Fill in missing years with column mean
  #
  # env_index is read BY ROW POSITION -- row r is model year styr + r - 1, and
  # the Year label is dropped below -- so a row outside styr:projyr would shift
  # every later row and feed the previous year's temperature to consumption.
  data_list$env_index <- merge(data_list$env_data, data.frame(Year = data_list$styr:data_list$projyr), all = TRUE)
  data_list$env_index <- data_list$env_index[
    data_list$env_index$Year %in% data_list$styr:data_list$projyr, , drop = FALSE]
  data_list$env_index <- data_list$env_index[order(data_list$env_index$Year), , drop = FALSE]
  # A repeated Year is inside the range, so the filter above does not catch it,
  # and one extra row shifts every year after it. data_check() reports it by
  # name; this keeps the row-per-model-year contract even when that is bypassed.
  data_list$env_index <-
    data_list$env_index[!duplicated(data_list$env_index$Year), , drop = FALSE]
  data_list$env_index <-  data_list$env_index %>%
    dplyr::select(-Year) %>%
    dplyr::mutate_all(~ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x)) %>%
    as.matrix()


  # 14 - Set up ration_data ----
  # - Convert to array
  data_list$ration_data <- data_list$ration_data %>%
    dplyr::mutate(Species = as.numeric(as.character(.data$Species)),
                  Sex = as.numeric(as.character(.data$Sex)),
                  Year = as.numeric(as.character(.data$Year)),
                  Year = ifelse(Year == 0,
                                0,
                                Year - data_list$styr + 1)
    )

  ration_data <- array(0, dim = c(data_list$nspp, max_sex, max_age, nyrs_hind)) # FIXME: Change for forecast

  if( nrow(data_list$ration_data) > 0){
    for (i in 1:nrow(data_list$ration_data)) {
      sp <- data_list$ration_data$Species[i]
      sex <- data_list$ration_data$Sex[i]
      yr <- data_list$ration_data$Year[i]

      # Handle sex == 0 case for 2-sex species
      sex_values <- if(sex == 0) 1:data_list$nsex[sp] else sex

      if(yr <= (data_list$endyr - data_list$styr + 1)){

        # Handle year == 0 case
        if(yr == 0) {
          yr <- seq_len(nyrs_hind)
        }

        # Fill in years
        for(j in 1:length(sex_values)){
          ration_data[sp, sex_values[j], 1:data_list$nages[sp], yr] <- as.numeric(
            as.character(
              as.matrix(data_list$ration_data[i, (1:data_list$nages[sp]) + 3])
            )
          )
        }
      }
    }
  }
  data_list$ration_data <- ration_data


  # 15 - Remove species column from maturity and sex_ratio ----
  data_list$sex_ratio <- data_list$sex_ratio %>%
    dplyr::select(contains("Age")) %>%
    dplyr::mutate_all(as.numeric) %>%
    as.matrix()

  data_list$maturity <- data_list$maturity %>%
    dplyr::select(contains("Age")) %>%
    dplyr::mutate_all(as.numeric) %>%
    as.matrix()


  # 16 - Make data.frames into matrices ----
  df_to_mat <- which(sapply(data_list, function(x) class(x)[1]) == "data.frame")
  data_list[df_to_mat] <- lapply(data_list[df_to_mat], as.matrix)

  # model_config is a code-side configuration slot (a nested list of build_*()
  # spec objects), not TMB data; strip it here alongside the other list-of-spec
  # objects so it does not ride into obj$env$data (inert to the objective, but it
  # bloats the fit and relies on TMB's sanitizer tolerating an arbitrary list).
  items_to_remove <- c("emp_sel",  "fsh_comp",    "srv_comp",    "catch_data", "equil_catch_data", "index_data", "comp_data", "caal_data", "env_data", "spnames",
                       "aLW", "diet_data", "index_cov", "model_config", # "NByageFixed", "estDynamics", "Ceq",
                       "avgnMode", "minNByage", "weight", "fleet_control")
  data_list[items_to_remove] <- NULL

  # 17 - One-step-ahead (OSA) residual vector ----
  # Assemble the flat observation vector (obsvec), the position-map metadata
  # (obs_ctl), and the per-type obsvec index vectors that the TMB template uses
  # to compute OSA residuals (see build_osa_data()). Built from the *_ctl/*_obs
  # matrices above. The aggregate index/catch entries are always built (the
  # template reads them); the costly composition/caal/diet metadata is built
  # only when build_osa = TRUE (off by default for fast simulation refits). This
  # runs AFTER section 16 on purpose: obs_ctl is a mixed-type data frame (R-side
  # metadata) and must not be swept into the data.frame->matrix coercion above.
  data_list <- build_osa_data(data_list, build_osa = build_osa)

  return(data_list)
}



#' Check and Clean Composition Data
#'
#' This function checks the composition data for zero sum rows, handles NA values,
#' and verifies that the composition data spans the required range of ages/lengths.
#'
#' @param data_list A list containing the following components:
#'   - comp_obs: A matrix or data frame of composition observations.
#'   - comp_data: A data frame with metadata including species, sex, and age/length information.
#'   - nsex: A vector of sex counts by species.
#'   - nages: A vector of age counts by species.
#'   - nlengths: A vector of length counts by species.
#' throws An error if any rows of `comp_obs` sum to 0 or if the composition data does not span the range of ages/lengths.
#'
#' @return The modified `data_list` with NA values in `comp_obs` converted to 0.
#' @examples
#' data_list <- list(
#'   comp_obs = matrix(c(1, 2, 3, 4, 5, 6), nrow = 2),
#'   comp_data = data.frame(Species = c(1, 2), Sex = c(1, 1), Age0_Length1 = c(0, 0)),
#'   nsex = c(1, 1),
#'   nages = c(3, 3),
#'   nlengths = c(3, 3)
#' )
#' cleaned_data_list <- check_composition_data(data_list)
#'
#' @keywords internal
#' @noRd
check_composition_data <- function(data_list) {

  # If no data, convert to empty matrix
  if(is.null(dim(data_list$comp_obs))){
    data_list$comp_obs <- matrix(NA, ncol = 10, nrow = 0)
  } else if(nrow(data_list$comp_obs) > 0){

    # # Check for zero sum rows in composition data
    # if (any(rowSums(data_list$comp_obs, na.rm = TRUE) == 0)) {
    #   stop("Some rows of composition data sum to 0: please remove or set all to 1 and sample size to 0")
    # }

    # Calculate adjustments for sex and age/length
    joint_adjust <- ifelse(data_list$nsex[data_list$comp_data$Species] == 2 &
                             data_list$comp_data$Sex == 3, 2, 1)

    col_adjust <- ifelse(data_list$comp_data$Age0_Length1 == 0,
                         data_list$nages[data_list$comp_data$Species],
                         data_list$nlengths[data_list$comp_data$Species])

    col_adjust <- col_adjust * joint_adjust

    # Check for NAs in composition data and convert them to 0
    na_check <- apply(cbind(col_adjust, data_list$comp_obs), 1, function(x) sum(x[2:(x[1] + 1)]))

    if (any(is.na(na_check))) {
      na_rows <- which(is.na(na_check))
      message(sprintf("Composition data have NAs in row(s): %s. Converting to 0s.",
                      paste(na_rows, collapse = ", ")))
      data_list$comp_obs[is.na(data_list$comp_obs)] <- 0
    }

    # Verify composition data size
    if (ncol(data_list$comp_obs) < max(col_adjust)) {
      stop("Comp data does not span range of ages/lengths")
    }
  }

  return(data_list)
}


#' Check and Clean CAAL Data
#'
#' This function checks the CAAL data for zero sum rows, handles NA values,
#' and verifies that the CAAL data spans the required range of ages
#'
#' @param data_list A list containing the following components:
#'   - caal_obs: A matrix or data frame of composition observations.
#'   - caal_data: A data frame with metadata including species, sex, and age/length information.
#'   - nsex: A vector of sex counts by species.
#'   - nages: A vector of age counts by species.
#'   - nlengths: A vector of length counts by species.
#' throws An error if any rows of `caal_obs` sum to 0 or if the composition data does not span the range of ages
#'
#' @return The modified `data_list` with NA values in `caal_obs` converted to 0.
#' cleaned_data_list <- check_caal_data(data_list)
#'
#' @keywords internal
#' @noRd
check_caal_data <- function(data_list) {

  # If no data, convert to empty matrix
  if(is.null(dim(data_list$caal_obs)) | nrow(data_list$caal_obs) == 0){
    data_list$caal_obs <- matrix(NA, ncol = 10, nrow = 0)
  } else{

    # # Check for zero sum rows in composition data
    # if (any(rowSums(data_list$caal_obs, na.rm = TRUE) == 0)) {
    #   stop("Some rows of composition data sum to 0: please remove or set all to 1 and sample size to 0")
    # }

    # Calculate adjustments for joint sex and CAAL data
    joint_adjust <- ifelse(data_list$nsex[data_list$caal_data$Species] == 2 &
                             data_list$caal_data$Sex == 3, 2, 1)

    col_adjust <- data_list$nages[data_list$caal_data$Species]
    col_adjust <- col_adjust * joint_adjust

    # Check for NAs in composition data and convert them to 0
    na_check <- apply(cbind(col_adjust, data_list$caal_obs), 1, function(x) sum(x[2:(x[1] + 1)]))

    if (any(is.na(na_check))) {
      na_rows <- which(is.na(na_check))
      message(sprintf("CAAL data have NAs in row(s): %s. Converting to 0s.",
                      paste(na_rows, collapse = ", ")))
      data_list$caal_obs[is.na(data_list$caal_obs)] <- 0
    }

    # Verify CAAL data size
    if (ncol(data_list$caal_obs) < max(col_adjust)) {
      stop("CAAL data does not span range of ages")
    }
  }

  return(data_list)
}


#' @rdname rearrange_data
#' @description `rearrange_dat()` is a deprecated alias for `rearrange_data()`,
#'   kept so older scripts keep running. It is scheduled for removal in 6.0.0;
#'   use `rearrange_data()`.
#' @export
rearrange_dat <- function(data_list){
  .Deprecated("rearrange_data",
              msg = "rearrange_dat() is deprecated and will be removed in Rceattle 6.0.0; use rearrange_data().")
  rearrange_data(data_list)
}


#' Population length bins and their map onto the data length bins
#'
#' Adds `lengths_pop` (lower edges, cm; species x bins), `nlengths_pop` and
#' `pop_to_data_bin` (0-based data bin that each population bin sums into) to
#' `data_list`. The grid comes from `data_list$pop_lengths` (one vector for every
#' species, or a list per species); a species without one, or with empirical
#' growth, uses its data bins, so the map is the identity. A population bin below
#' the first data edge sums into the first data bin and one above the last edge
#' into the last, the data bins' minus and plus groups. Every data edge must also
#' be a population edge, so no population bin straddles two data bins.
#' @keywords internal
#' @noRd
.rce_pop_length_bins <- function(data_list) {
  nspp <- data_list$nspp
  pop  <- data_list$pop_lengths
  grids <- lapply(seq_len(nspp), function(sp) {
    data_edges <- as.numeric(data_list$lengths[sp, seq_len(data_list$nlengths[sp])])
    g <- if (is.list(pop)) (if (sp <= length(pop)) pop[[sp]] else NULL) else pop
    if (is.null(g) || isTRUE(data_list$growth_model[sp] == 0)) g <- data_edges
    missing_edge <- data_edges[vapply(data_edges, function(e) !any(abs(g - e) < 1e-8), logical(1))]
    if (length(missing_edge)) {
      stop(sprintf("pop_lengths for species %d must include every data length-bin lower edge; missing: %s",
                   sp, paste(missing_edge, collapse = ", ")), call. = FALSE)
    }
    bin <- pmax(findInterval(g + 1e-8, data_edges), 1L) - 1L
    list(edges = g, bin = as.integer(bin))
  })
  n_pop <- vapply(grids, function(x) length(x$edges), integer(1))
  data_list$nlengths_pop    <- as.integer(n_pop)
  data_list$lengths_pop     <- matrix(0, nspp, max(n_pop))
  data_list$pop_to_data_bin <- matrix(0L, nspp, max(n_pop))
  for (sp in seq_len(nspp)) {
    data_list$lengths_pop[sp, seq_len(n_pop[sp])]     <- grids[[sp]]$edges
    data_list$pop_to_data_bin[sp, seq_len(n_pop[sp])] <- grids[[sp]]$bin
  }
  data_list
}
