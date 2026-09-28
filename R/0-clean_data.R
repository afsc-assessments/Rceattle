# Placeholder multispecies unfished SSB / biomass, metric tons. `MSSB0` and
# `MSB0` are required DATA_VECTORs that no workbook supplies, so clean_data()
# seeds them with this and fit_mod() section 10.2 replaces it with a real
# no-fishing projection -- but only for a run that has a harvest control rule.
# A value still equal to this one means no unfished reference was ever derived,
# so anything reading SSB against it reports a ratio, not a depletion.
.RCE_MSSB0_PLACEHOLDER <- 999


#' Align survey-index covariance matrices to the current fitted observations
#'
#' `index_cov` Sigma matrices are positional: rows/cols follow a fleet's fitted
#' survey observations (Year in (0, endyr], Observation > 0) in `index_data`
#' order. Sigma dimnames are tagged with those years the first time it matches
#' the data, then re-keyed to the current fitted set on each `clean_data()` pass:
#' retained years keep their covariance block, new years (future / simulated
#' observations) get an independent diagonal `(Observation * Log_sd)^2`, and
#' cross terms between the two are 0. Non-MVN fleets pass through unchanged.
#' Assumes one fitted observation per year per fleet.
#'
#' @keywords internal
#' @noRd
.align_index_cov <- function(data_list) {
  ic <- data_list$index_cov
  if (is.null(ic) || length(ic) == 0) return(data_list)
  fc  <- data_list$fleet_control
  idx <- data_list$index_data
  if (is.null(fc) || is.null(idx) || !"Index_distribution" %in% colnames(fc)) return(data_list)
  endyr <- data_list$endyr
  mvn_codes <- c("MVN", "MVNORM", 1, 2, "1", "2")

  for (nm in names(ic)) {
    Sigma <- as.matrix(ic[[nm]])
    row <- which(fc$Fleet_name == nm | as.character(fc$Fleet_code) == nm)
    if (length(row) != 1) next
    if (!(fc$Index_distribution[row] %in% mvn_codes)) next        # only MVN/MVNORM fleets
    code <- fc$Fleet_code[row]

    fit <- idx[idx$Fleet_code == code & idx$Year > 0 &
                 idx$Year <= endyr & idx$Observation > 0, , drop = FALSE]
    now_yrs <- fit$Year

    have_yrs <- suppressWarnings(as.integer(rownames(Sigma)))
    if (length(have_yrs) == 0 || any(is.na(have_yrs))) have_yrs <- NULL

    if (is.null(have_yrs)) {
      # First sight: tag with the fitted years only if the Sigma already matches
      # the current fitted set; otherwise leave it for rearrange_data() to report.
      if (nrow(Sigma) == length(now_yrs) && length(now_yrs) > 0) {
        dimnames(Sigma) <- list(as.character(now_yrs), as.character(now_yrs))
      }
      ic[[nm]] <- Sigma
      next
    }

    # Re-key the tagged Sigma to the current fitted years (subset / reorder / grow).
    n   <- length(now_yrs)
    S2  <- matrix(0, n, n)
    pos <- match(now_yrs, have_yrs)                          # NA => a new (untagged) year
    keep <- !is.na(pos)
    if (any(keep)) S2[keep, keep] <- Sigma[pos[keep], pos[keep], drop = FALSE]
    new_i <- which(!keep)
    if (length(new_i)) {
      v <- (fit$Observation[new_i] * fit$Log_sd[new_i])^2
      S2[cbind(new_i, new_i)] <- v
    }
    if (n > 0) dimnames(S2) <- list(as.character(now_yrs), as.character(now_yrs))
    ic[[nm]] <- S2
  }
  data_list$index_cov <- ic
  data_list
}

#' Default-fill and tidy a data list before fitting
#'
#' Fills optional blocks with correctly-shaped empties, filters observations to
#' the model's year range, extends catch to the projection years, and re-keys
#' any survey-index covariance matrices. Whether a missing block is actually a
#' problem is left to the fit-time validation, which knows the model configuration.
#'
#' @param data_list Rceattle data list
#'
#' @export
#'
clean_data <- function(data_list){

  # --- 0. Default-fill optional data.frame fields ----
  # These fields can be NULL when the user is not using the corresponding feature
  # (e.g., comp_data in a model without composition likelihoods, ration_data /
  # diet_data in single-species mode, NByageFixed when estDynamics == 0).
  # Filling with empty data.frames that carry the metadata columns the
  # downstream code expects lets rearrange_data / build_params use uniform
  # `nrow > 0` checks without separate NULL guards. Whether the missing data is
  # actually a problem is enforced by data_check() based on the model
  # configuration (msmMode, growth_model, estDynamics, Selectivity, etc.).
  default_dfs <- list(
    comp_data   = data.frame(Fleet_code = integer(), Species = integer(),
                             Sex = integer(), Age0_Length1 = integer(),
                             Year = integer(), Month = integer(),
                             Sample_size = numeric()),
    caal_data   = data.frame(Fleet_code = integer(), Species = integer(),
                             Sex = integer(), Year = integer(),
                             Length = numeric(), Sample_size = numeric()),
    emp_sel     = data.frame(Fleet_code = integer(), Species = integer(),
                             Sex = integer(), Year = integer()),
    NByageFixed = data.frame(Species = integer(), Sex = integer(),
                             Year = integer()),
    ration_data = data.frame(Species = integer(), Sex = integer(),
                             Year = integer()),
    diet_data   = data.frame(Pred = integer(), Prey = integer(),
                             Pred_sex = integer(), Prey_sex = integer(),
                             Pred_age = integer(), Prey_age = integer(),
                             Year = integer(),
                             Sample_size = numeric(),
                             Stomach_proportion_by_weight = numeric())
  )
  for(df_name in names(default_dfs)){
    if(is.null(data_list[[df_name]])){
      data_list[[df_name]] <- default_dfs[[df_name]]
    }
  }

  # Selectivity_block: optional per-observation time-block id in index_data /
  # catch_data. It is READ ONLY for Time_varying_{sel,q} == "Block" (build_map);
  # every other configuration (Off / IID / RW / AR1, and the random-effect
  # linkages) ignores it. Default a missing column to 1 (a single block = no
  # time-blocking), so a user need only supply it for Block-mode fleets. (The
  # `Q_block` column is vestigial -- never read; the q time-blocking reuses
  # Selectivity_block -- see data_check() for its soft-deprecation.)
  for (.blk in c("index_data", "catch_data")) {
    if (!is.null(data_list[[.blk]]) && nrow(data_list[[.blk]]) > 0L &&
        !"Selectivity_block" %in% names(data_list[[.blk]])) {
      data_list[[.blk]]$Selectivity_block <- 1L
    }
  }

  # env_data: default to a Year-only data.frame so downstream code that does
  # `ncol(env_data) - 1` (number of indices) gets 0 when no environmental data
  # are supplied, and `merge(env_data, ...)` in rearrange_data() works.
  # data_check() still errors when a model feature (env-q catchability,
  # temperature-dependent consumption, env linkages, srr/M1 indices) needs
  # actual indices.
  if(is.null(data_list$env_data)){
    yrs <- if(!is.null(data_list$styr) && !is.null(data_list$projyr))
             data_list$styr:data_list$projyr else integer(0)
    data_list$env_data <- data.frame(Year = yrs)
  }

  # index_cov: named list of survey-index variance-covariance matrices, keyed by
  # Fleet_name, used only by fleets with Index_distribution == "MVN" (the AMAK/ebswp
  # DoCovBTS covariance survey likelihood). Default to an empty list so
  # rearrange_data() / data_check() can treat "not supplied" uniformly; a fleet
  # that is not MVN never needs one, and data_check() enforces that an MVN fleet
  # has a matrix of the correct dimension.
  if(is.null(data_list$index_cov)){
    data_list$index_cov <- list()
  }

  # --- 1. Filter Data by Year ----
  # Data in likelihood (use absolute Year)
  abs_year_data <- c("index_data", "comp_data", "caal_data")
  for(df_name in abs_year_data) {
    if(!is.null(data_list[[df_name]])) {
      data_list[[df_name]] <- data_list[[df_name]] |>
        dplyr::filter(abs(Year) >= data_list$styr & abs(Year) <= data_list$projyr)
    }
  }

  # catch_data admits one year more than the others: styr - 1 may hold the
  # initial equilibrium catch. Whether it does depends on initMode, which
  # rearrange_data() applies -- under a mode that holds Finit at 0 the row is
  # ordinary history and is dropped there, as it always was.
  if(!is.null(data_list$catch_data)) {
    data_list$catch_data <- data_list$catch_data |>
      dplyr::filter(abs(Year) >= (data_list$styr - 1L) & abs(Year) <= data_list$projyr)
  }

  # Fixed data (allow Year == 0)
  fixed_year_data <- c("diet_data", "weight", "emp_sel", "NByageFixed", "ration_data")
  for(df_name in fixed_year_data) {
    if(!is.null(data_list[[df_name]])) {
      data_list[[df_name]] <- as.data.frame(data_list[[df_name]]) |>
        dplyr::filter((Year >= data_list$styr & Year <= data_list$projyr) | Year == 0)
    }
  }

  # --- 1b. Re-key survey-index covariance matrices to the current fitted
  #         survey observations (see .align_index_cov above).
  data_list <- .align_index_cov(data_list)

  # --- 2. Multi-species unfished SSB and biomass (metric tons) ----
  # Required DATA_VECTORs, and no workbook can supply them: neither has
  # read_data()/write_data() support, so this is the only thing that creates
  # them. The model derives SB0/B0 itself and reads these only to overwrite
  # it under msmMode > 0, where fit_mod() fills them from a no-fishing
  # projection; .RCE_MSSB0_PLACEHOLDER stands in until it does.
  if(is.null(data_list$MSSB0)){
    data_list$MSSB0 <- rep(.RCE_MSSB0_PLACEHOLDER, data_list$nspp)
    data_list$MSB0 <- rep(.RCE_MSSB0_PLACEHOLDER, data_list$nspp)
  }
  # Per species: has an unfished reference actually been derived? Seeded FALSE
  # alongside the placeholder and set by fit_mod() section 10.2 when it projects
  # under no fishing. Carried rather than inferred from the value, because
  # recognising the placeholder by comparing a double to 999 would also null a
  # genuinely-derived 999 mt, and is blind to a workbook that supplied its own
  # MSSB0. Not a TMB input -- the model never sees it.
  if(is.null(data_list$MSSB0_derived)){
    data_list$MSSB0_derived <- rep(FALSE, data_list$nspp)
  }


  # --- 3. Extend catch data to proj year for projections ----
  if(data_list$projyr > data_list$endyr){
    for(flt in unique(data_list$catch_data$Fleet_code)){
      catch_data_sub <- data_list$catch_data |> dplyr::filter(Fleet_code == flt)

      yrs_proj <- (data_list$endyr + 1):data_list$projyr
      yrs_proj <- yrs_proj[which(!yrs_proj %in% catch_data_sub$Year)]
      nyrs_proj <- length(yrs_proj)

      if(nyrs_proj > 0) {
        proj_catch_data <- data.frame(
          Fleet_name = rep(catch_data_sub$Fleet_name[1], nyrs_proj),
          Fleet_code = rep(flt, nyrs_proj),
          Species = rep(catch_data_sub$Species[1], nyrs_proj),
          Year = yrs_proj,
          Month = rep(dplyr::last(catch_data_sub$Month), nyrs_proj),
          Selectivity_block = rep(dplyr::last(catch_data_sub$Selectivity_block), nyrs_proj),
          Catch = rep(NA, nyrs_proj),
          Log_sd = rep(dplyr::last(catch_data_sub$Log_sd), nyrs_proj)
        )
        data_list$catch_data <- rbind(data_list$catch_data, proj_catch_data)
      }
    }
  }


  # --- 4. Column names ----
  expected_cols <- c("Species", paste0("Age", 1:max(data_list$nages)))

  if(any(!colnames(data_list$sex_ratio) %in% expected_cols)){
    colnames(data_list$sex_ratio) <- expected_cols
    message("Renaming column names in 'sex_ratio' data to 'Age1', 'Age2', ....")
  }

  if(any(!colnames(data_list$maturity) %in% expected_cols)){
    colnames(data_list$maturity) <- expected_cols
    message("Renaming column names in 'maturity' data to 'Age1', 'Age2', ....")
  }


  # --- 5. Arrange diet data ----
  if(!is.null(data_list$diet_data)){
    data_list$diet_data <- data_list$diet_data |>
      dplyr::arrange(Pred, Pred_sex, Pred_age, Prey, Prey_sex, Prey_age, Year) |>
      dplyr::mutate(stratum_id = paste(Pred, Pred_sex, Pred_age, Year, sep = "_"),
                    stomach_id = as.numeric(as.factor(stratum_id)) - 1) |>
      dplyr::arrange(stomach_id)
  }

  return(.rce_as_data(data_list))
}


#' Locate the observed-proportion columns of a composition table
#'
#' @description
#' `comp_data` and `caal_data` each begin with identifying columns (fleet,
#' species, sex, year, sample size, ...) followed by the observed proportions
#' at each age or length bin (`Comp_1`, `Comp_2`, ... / `CAAL_1`, `CAAL_2`,
#' ...). The number of identifying columns differs between the two tables and
#' is not fixed, so the proportion columns must be found by name rather than by
#' counting from a fixed position. `rearrange_data()` does this with
#' `dplyr::contains()`; this helper is the base-R equivalent for subsetting and
#' assignment.
#'
#' @param x a composition table (`comp_data` or `caal_data`).
#' @param prefix column-name prefix identifying the proportion columns, e.g.
#'   `"Comp_"` or `"CAAL_"`.
#'
#' @return Integer vector of column indices, in table order. Throws if no
#'   column matches, rather than returning a silent zero-length result.
#'
#' @keywords internal
#' @noRd
.composition_cols <- function(x, prefix) {
  idx <- grep(paste0("^", prefix), names(x))
  if (length(idx) == 0L) {
    stop(sprintf(
      "no '%s*' columns found in the supplied table; available columns: %s",
      prefix, paste(names(x), collapse = ", ")), call. = FALSE)
  }
  idx
}
