

#' Label the derived quantities reported by a CEATTLE fit
#'
#' Attaches dimension names (species, sex, age, length bin, year, fleet) to the
#' report objects returned by the TMB model, and names the rows/columns of the
#' likelihood-component matrices, so the fitted quantities read as labeled
#' arrays rather than bare numbers.
#'
#' @param data_list an Rceattle data_list
#' @param quantities list of "report" objects from Rceattle.
#'
#' @export
#'
rename_output <- function(data_list = NULL, quantities = NULL){

  # Dimension attributed
  max_age <- max(data_list$nages, na.rm = TRUE)
  max_sex <- max(data_list$nsex, na.rm = TRUE)
  max_length <- max(data_list$nlengths, na.rm = TRUE)
  sex_labels <- c("Sex combined or females", "males")
  if(max_sex == 1){
    sex_labels <- "Sex combined"
  }
  yrs_hind <- data_list$styr:data_list$endyr
  yrs_proj <- data_list$styr:data_list$projyr
  nyrs_hind <- length(yrs_hind)
  nyrs_proj <- length(yrs_proj)


  # Rename ----
  # * Vectors ----
  #  Biological quantities
  names(quantities$avg_R) <- data_list$spnames
  names(quantities$Flimit) <- data_list$spnames
  names(quantities$Ftarget) <- data_list$spnames
  names(quantities$gam_a) <- data_list$spnames
  names(quantities$gam_b) <- data_list$spnames
  names(quantities$R_init) <- data_list$spnames
  names(quantities$SPR0) <- data_list$spnames
  names(quantities$SPRFinit) <- data_list$spnames
  names(quantities$SPRlimit) <- data_list$spnames
  names(quantities$SPRtarget) <- data_list$spnames
  names(quantities$steepness) <- data_list$spnames

  # Input numbers-at-age (estDynamics > 0): no HCR and F = 0 in projection, so these are NA.
  fixed_n <- (data_list$estDynamics %||% rep(0, data_list$nspp)) > 0
  if (any(fixed_n)) {
    # Its rec_pars are fixed, so the stock-recruit quantities are the build_params()
    # placeholder (R0 = exp(9)), not something the model or the user set.
    mask <- c("Ftarget", "Flimit", "SPRtarget", "SPRlimit", "SBF", "DynamicSBF",
              "R0", "R_init", "avg_R", "steepness", "SPR0")
    # R is its input recruits (thousands of fish): first-age N-at-age, which is
    # NByageFixed times pop_scalar, summed over sexes.
    for (sp in which(fixed_n)) {
      quantities$R[sp, ] <- apply(
        quantities$N_at_age[sp, seq_len(data_list$nsex[sp]), 1, , drop = FALSE], 4, sum)
      # A year with no NByageFixed row holds zeros, not recruits.
      yrs_in <- data_list$NByageFixed$Year[data_list$NByageFixed$Species == sp]
      quantities$R[sp, !yrs_proj %in% yrs_in] <- NA
    }
    # In single-species mode the equilibrium SB0 and B0 are built on that placeholder,
    # and with DynamicHCR = FALSE the depletions divide by them. Under predation MSSB0
    # replaces SB0; with DynamicHCR = TRUE the depletions are the input numbers
    # relative to themselves.
    if (isTRUE(as.integer(data_list$msmMode %||% 0L)[1] == 0L)) {
      mask <- c(mask, "SB0", "B0")
      if (!isTRUE(as.logical(data_list$DynamicHCR %||% FALSE)))
        mask <- c(mask, "ssb_depletion", "biomass_depletion")
    }
    for (nm in mask) {
      x <- quantities[[nm]]
      if (is.null(x)) next
      if (is.matrix(x)) x[fixed_n, ] <- NA else x[fixed_n] <- NA
      quantities[[nm]] <- x
    }
  }

  # * Fleets ----
  names(quantities$catch_sd) <- data_list$catch_data$Fleet_name
  names(quantities$index_sd) <- data_list$index_data$Fleet_name

  # Deprecated spellings, kept one release so existing scripts and saved
  # workflows keep running. Neither was ever a log: both hold the observation sd
  # the likelihood actually used, and for a natural-scale Index_distribution
  # index_sd is an ABSOLUTE sd in index units, which is exactly what the old
  # name talked people out of believing. Removed in the next minor release.
  quantities$log_catch_sd <- quantities$catch_sd
  quantities$log_index_sd <- quantities$index_sd


  # * 2D array ----
  # - Population quantities
  dimnames(quantities$biomass) <- list(data_list$spnames, yrs_proj)
  dimnames(quantities$biomass_depletion) <- list(data_list$spnames, yrs_proj)
  dimnames(quantities$B0) <- list(data_list$spnames, yrs_proj)
  dimnames(quantities$DynamicB0) <- list(data_list$spnames, yrs_proj)
  dimnames(quantities$exploitable_biomass) <- list(data_list$spnames, yrs_proj)
  dimnames(quantities$R0) <- list(data_list$spnames, yrs_proj)

  dimnames(quantities$ssb_depletion) <- list(data_list$spnames, yrs_proj)
  dimnames(quantities$ssb) <- list(data_list$spnames, yrs_proj)
  dimnames(quantities$SB0) <- list(data_list$spnames, yrs_proj)
  dimnames(quantities$DynamicSB0) <- list(data_list$spnames, yrs_proj)
  dimnames(quantities$DynamicSBF) <- list(data_list$spnames, yrs_proj)

  dimnames(quantities$vulnerability) <- list(paste0("Pred: ", data_list$spnames), paste0("Prey: ", data_list$spnames)) # Pred/prey?

  dimnames(quantities$R) <- list(data_list$spnames, yrs_proj)

  dimnames(quantities$F_spp) <- list(data_list$spnames, yrs_proj)
  dimnames(quantities$proj_F) <- list(data_list$spnames, yrs_proj)

  dimnames(quantities$fT) <- list(data_list$spnames, yrs_proj) # Temperature function of consumption

  names(quantities$pop_scalar) <- data_list$spnames


  # - Fleet quantities
  dimnames(quantities$F_flt) <- list(data_list$fleet_control$Fleet_name, yrs_proj) # Sex specific?
  dimnames(quantities$index_q) <- list(data_list$fleet_control$Fleet_name, yrs_hind)

  # * 4D array ----
  # - Biological
  dimnames(quantities$biomass_at_age) <- list(data_list$spnames, sex_labels, paste0("Age", 1:max_age), yrs_proj)
  dimnames(quantities$consumption_at_age) <- list(data_list$spnames, sex_labels, paste0("Age", 1:max_age), yrs_proj)
  dimnames(quantities$B_eaten_as_prey) <-  list(data_list$spnames, sex_labels, paste0("Age", 1:max_age), yrs_proj)
  dimnames(quantities$F_at_age) <- list(data_list$spnames, sex_labels, paste0("Age", 1:max_age), yrs_proj)
  dimnames(quantities$M_at_age) <- list(data_list$spnames, sex_labels, paste0("Age", 1:max_age), yrs_proj)
  dimnames(quantities$M1_at_age) <- list(data_list$spnames, sex_labels, paste0("Age", 1:max_age), yrs_proj)
  dimnames(quantities$M2_at_age) <- list(data_list$spnames, sex_labels, paste0("Age", 1:max_age), yrs_proj)
  dimnames(quantities$Z_at_age) <- list(data_list$spnames, sex_labels, paste0("Age", 1:max_age), yrs_proj)
  dimnames(quantities$N_at_age) <- list(data_list$spnames, sex_labels, paste0("Age", 1:max_age), yrs_proj)
  dimnames(quantities$NByage0) <- list(data_list$spnames, sex_labels, paste0("Age", 1:max_age), yrs_proj)
  dimnames(quantities$NByageF) <- list(data_list$spnames, sex_labels, paste0("Age", 1:max_age), yrs_proj)
  # 3D, and the only quantity whose first dimension is a reference point rather
  # than a species: numbers at age per recruit, in the order section 6.2 writes
  # them. Unlabelled it reads as NbyageSPR[3, 1, ], which requires knowing that
  # 3 is Finit. Ages follow the convention above -- the label is the age bin
  # index, so it is the age itself only when minage is 1.
  dimnames(quantities$NbyageSPR) <- list(c("F0", "Flimit", "Ftarget", "Finit"),
                                         data_list$spnames, paste0("Age", 1:max_age))
  # `growth_parameters` and `growth_linkage_offset` carry only the
  # mean-growth params (log_K, log_L1, log_Linf, log_m); the SD
  # endpoints in `GROWTH_LINKAGE_PARAMS` live on `growth_log_sd` and
  # have no year dim, so they are not part of this last dimension.
  dimnames(quantities$growth_parameters) <- list(data_list$spnames, sex_labels, yrs_proj, .GROWTH_MEAN_PARAMS)
  dimnames(quantities$growth_linkage_offset) <- list(data_list$spnames, sex_labels, yrs_proj, .GROWTH_MEAN_PARAMS)

  dimnames(quantities$weight_hat) <- dimnames(quantities$length_hat) <- list(
    c(paste(rep(data_list$spnames, each = 2), rep(c("biomass length", "spawn length"), data_list$nspp)), data_list$fleet_control$Fleet_name),
    sex_labels, paste0("Age", 1:max_age), yrs_proj)
  if (!is.null(quantities$spawn_output)) {
    dimnames(quantities$spawn_output) <- list(data_list$spnames, paste0("Age", 1:max_age), yrs_proj)
  }
  # The age-length key and selectivity-at-length are on the POPULATION length
  # bins. The count comes off the arrays, since `data_list` here is the
  # pre-rearrange_data() list and carries no nlengths_pop.
  .dim_or <- function(x, i, alt) if (!is.null(x) && length(dim(x)) >= i) dim(x)[i] else alt
  max_length_pop <- .dim_or(quantities$growth_matrix, 4L, max_length)
  dimnames(quantities$growth_matrix) <- list(
    c(paste(rep(data_list$spnames, each = 2), rep(c("biomass weight", "spawn weight"), data_list$nspp)), data_list$fleet_control$Fleet_name),
    sex_labels, paste0("Age", 1:max_age), paste0("PopBin", 1:max_length_pop), yrs_proj)

  # - Fleet
  dimnames(quantities$F_flt_age) <- list(data_list$fleet_control$Fleet_name, sex_labels, paste0("Age", 1:max_age), yrs_proj)
  dimnames(quantities$sel_at_age) <- list(data_list$fleet_control$Fleet_name, sex_labels, paste0("Age", 1:max_age), yrs_proj)
  dimnames(quantities$sel_at_length) <- list(data_list$fleet_control$Fleet_name, sex_labels,
    paste0("PopBin", 1:.dim_or(quantities$sel_at_length, 3L, max_length_pop)), yrs_proj)

  # * 5D arrays ----
  dimnames(quantities$B_eaten) <- list(paste("Pred:", data_list$spnames, rep(sex_labels, each = data_list$nspp)),
                                       paste("Prey:", data_list$spnames, rep(sex_labels, each = data_list$nspp)),
                                       paste0("Pred age", 1:max_age),
                                       paste0("Prey age", 1:max_age),
                                       yrs_proj)
  dimnames(quantities$suitability) <- list(paste("Pred:", data_list$spnames, rep(sex_labels, each = data_list$nspp)),
                                           paste("Prey:", data_list$spnames, rep(sex_labels, each = data_list$nspp)),
                                           paste0("Pred age", 1:max_age),
                                           paste0("Prey age", 1:max_age),
                                           yrs_proj)


  # Likelihood components ----
  S = paste0(data_list$spnames, " or ")
  X = data_list$fleet_control$Fleet_name
  llcolnames <- paste0(append(S, rep("", max(c(0, length(X)-length(S))))), append(X, rep(NA, max(c(0, length(S)-length(X))))))

  colnames(quantities$jnll_comp) <- llcolnames
  colnames(quantities$unweighted_jnll_comp) <- llcolnames

  rownames(quantities$jnll_comp) <- rownames(quantities$unweighted_jnll_comp) <- c(
    "Index data",
    "Catch data",
    "Composition data",
    "CAAL data",
    "Non-parametric selectivity",
    "Selectivity deviates",
    "Catchability prior",
    "Catchability deviates",
    "Stock-recruit prior",
    "Initial abundance deviates",
    "Recruitment deviates",
    "Stock-recruit penalty",
    "Reference point penalties",
    "Zero n-at-age penalty",
    "M prior",
    "M random effects",
    "Ration",
    "Ration penalties",
    "Stomach content data",
    "Linkage-table priors",
    "Linkage random effects",
    "Initial equilibrium catch"
  )

  return(quantities)
}


#' Calculate McAllister-Ianelli composition weights
#'
#' @param data_list an Rceattle data_list
#' @param data_list_reorganized reorganized data_list
#' @param quantities list of "report" objects from Rceattle.
#'
#' @keywords internal
#' @noRd
calc_mcall_ianelli <- function(data_list = NULL, data_list_reorganized = NULL, quantities = NULL){


  # - Calculate Mcallister-Iannelli coefficients
  # Effective sample size for the length data for year y
  eff_n_mcallister <- rowSums(quantities$comp_hat * (1 - quantities$comp_hat), na.rm = TRUE)/rowSums((data_list_reorganized$comp_obs - quantities$comp_hat)^2, na.rm = TRUE) # sum_length (p_hat * (1 - p_hat))/ sum_length ((p - p_hat) ^ 2)


  # Loop fleets and take harmonic mean
  data_list$fleet_control$Comp_weights_mcallister <- NA
  for(flt in unique(data_list$comp_data$Fleet_code)){
    comp_sub <- which(data_list$comp_data$Fleet_code == flt & data_list$comp_data$Year > 0)
    data_list$fleet_control$Comp_weights_mcallister[which(data_list$fleet_control$Fleet_code == flt)] <- ((1/length(comp_sub))*sum((eff_n_mcallister[comp_sub]/data_list$comp_data$Sample_size[comp_sub])^-1))^-1
  }

  # Back-compat: this output column was named `Est_weights_mcallister` before
  # v4.12 (renamed to `Comp_weights_mcallister` for parity with the diet output).
  # Mirror the old name so downstream reweighting scripts that read it off a
  # fitted object keep working. Deprecated; prefer `Comp_weights_mcallister`.
  data_list$fleet_control$Est_weights_mcallister <-
    data_list$fleet_control$Comp_weights_mcallister

  return(data_list)
}





#' Calculate McAllister-Ianelli weights for diet composition data
#'
#' @param data_list an Rceattle data_list
#' @param quantities list of "report" objects from Rceattle, including diet_hat predictions
#'
#' @keywords internal
#' @noRd
calc_mcall_ianelli_diet <- function(data_list = NULL, quantities = NULL){

  # - Calculate Mcallister-Iannelli coefficients for diet data
  diet_multiplier <- data_list$Diet_comp_weights

  # Small constant to avoid division by zero
  epsilon <- 1e-10


  # Calculate effective sample size for diet data (predator specific)
  # Using the same formula as for length: sum(p_hat * (1 - p_hat)) / sum((p - p_hat)^2)
  eff_n_mcallister <- data_list$diet_data |>
    dplyr::mutate(Diet_hat = quantities$diet_hat) |>
    dplyr::group_by(Pred, Pred_age) |>
    dplyr::summarise(
      Sample_size = dplyr::first(Sample_size), # Sample size should be the same across predators of the same age
      eff_n_mcallister = sum(Diet_hat * (1 - Diet_hat), na.rm = TRUE) /
        (sum((Stomach_proportion_by_weight - Diet_hat)^2, na.rm = TRUE) + epsilon)
    )

  # Take harmonic mean across predator ages
  data_list$Diet_weights_mcallister <- eff_n_mcallister |>
    dplyr::group_by(Pred) |>
    filter(eff_n_mcallister != 0) |>
    dplyr::summarise(Diet_weights_mcallister = (1/n() * sum((eff_n_mcallister /Sample_size)^-1))^-1 ) |>
    dplyr::pull(Diet_weights_mcallister)

  return(data_list)
}


#' Read an observation-sd vector under either spelling
#'
#' `index_sd` / `catch_sd` were named `log_index_sd` / `log_catch_sd` until
#' 5.9.0. Neither was ever a log, and `index_sd` is an ABSOLUTE sd for a
#' natural-scale `Index_distribution`, so the old name actively misled. Fits
#' saved before the rename hold only the old name; a fresh report holds only
#' the new one. Everything downstream reads through this.
#'
#' @param quantities A fit's `$quantities` list, or a TMB report.
#' @param what `"index"` or `"catch"`.
#' @return The sd vector, or `NULL` when neither name is present.
#' @noRd
.observation_sd <- function(quantities, what = c("index", "catch")) {
  what <- match.arg(what)
  got <- quantities[[paste0(what, "_sd")]]
  if (is.null(got)) got <- quantities[[paste0("log_", what, "_sd")]]
  got
}

#' Replace a covariance fleet's reported sd with its marginal sd
#'
#' An MVN / MVNORM fleet scores `index_cov` and leaves the reported `index_sd`
#' at whatever its `index_data` `Log_sd` column carried, which that likelihood
#' never reads, so standardizing those rows by it can be wrong by orders of
#' magnitude. Their sd is `sqrt(diag(index_cov))`, the
#' MARGINAL sd: it ignores the off-diagonal correlation the likelihood does use,
#' so a residual standardized this way is marginal, not fully whitened.
#'
#' @param data_list The fit's `data_list` (pre-`rearrange_data()`, so the input
#'   `index_cov`, not the `index_cov_mat` that `rearrange_data()` builds).
#' @param sigma `index_sd`, one value per `index_data` row.
#' @return `sigma`, with covariance-fleet rows replaced, and `NA` on such a row
#'   wherever the matrix cannot be resolved against the fitted years.
#' @noRd
.index_marginal_sd <- function(data_list, sigma) {
  cov_rows <- .index_rows_covariance(data_list)
  if (!length(cov_rows) || !any(cov_rows) || length(sigma) != length(cov_rows)) {
    return(sigma)
  }
  idx <- data_list$index_data
  fc  <- data_list$fleet_control
  for (code in unique(idx[["Fleet_code"]][cov_rows])) {
    i  <- match(code, fc[["Fleet_code"]])
    nm <- as.character(fc[["Fleet_name"]][i])
    Sigma <- data_list$index_cov[[nm]]
    if (is.null(Sigma)) Sigma <- data_list$index_cov[[as.character(code)]]
    Sigma <- tryCatch(as.matrix(Sigma), error = function(e) NULL)
    # Sigma's rows follow this fleet's FITTED observations in index_data order
    # and are year-dimnamed by clean_data(); see .align_index_cov().
    this <- cov_rows & idx[["Fleet_code"]] == code
    # !is.na() first: an NA observation would otherwise put NA in the mask, and
    # a logical subscript carrying NA errors on assignment.
    fitted <- this & idx[["Year"]] > 0 & idx[["Year"]] <= data_list$endyr &
      !is.na(idx[["Observation"]]) & idx[["Observation"]] > 0
    sigma[this] <- NA_real_
    if (is.null(Sigma) || !nrow(Sigma)) next
    # diag() of a non-square matrix, or of a bare vector passed as one, is
    # shorter than the fitted rows and would recycle silently.
    sd_diag <- sqrt(diag(Sigma))
    if (anyNA(sd_diag) || any(sd_diag <= 0)) next
    yrs <- suppressWarnings(as.integer(rownames(Sigma)))
    if (sum(fitted) == nrow(Sigma) && length(sd_diag) == nrow(Sigma)) {
      # Positional first: Sigma's rows ARE this fleet's fitted rows in
      # index_data order, which rearrange_data() enforces with a stop(), so a
      # mismatch is a fit-time error rather than a plot-time wrong number.
      sigma[fitted] <- sd_diag
    } else if (length(yrs) == nrow(Sigma) && !anyNA(yrs) &&
               !anyDuplicated(idx[["Year"]][fitted])) {
      # Year key for a Sigma whose row count no longer matches -- and only when
      # years are unique, since match() would give two rows in one year the
      # same entry. A legacy Index_loglike workbook reaches fit_mod() untagged,
      # because switch_check() renames the column after clean_data() runs.
      sigma[fitted] <- sd_diag[match(idx[["Year"]][fitted], yrs)]
    }
  }
  sigma
}

#' Read a fit's observation bias-adjustment flag
#'
#' The index and catch lognormal likelihoods fit `log(obs)` to a mean of
#' `log(hat) - bias_adjust_obs * sigma^2 / 2`, so anything comparing an
#' observation against `hat` needs the same flag. `data_list` carries it after
#' `fit_mod()`; older fits fall back to the TMB data, then to the default 1.
#'
#' @param object A fitted `Rceattle` object.
#' @return Length-one numeric: 1, 0, or a fraction where a partial correction
#'   was configured.
#' @noRd
.bias_adjust_obs <- function(object) {
  got <- object$data_list$bias_adjust_obs
  if (is.null(got) && !is.null(object$obj)) {
    got <- object$obj$env$data$bias_adjust_obs
  }
  if (is.null(got)) got <- 1
  as.numeric(got)[1]
}

#' Convert a reported mean observation to the median its likelihood centres on
#'
#' Log-scale rows only: a natural-scale family has no offset, and its `sigma` is
#' an ABSOLUTE sd, so the lognormal form would collapse the prediction to zero
#' (see `.observation_sd()`).
#'
#' @param hat Reported `index_hat` / `catch_hat`.
#' @param sigma The sd the likelihood used, from `.observation_sd()`.
#' @param ba_obs The flag, from `.bias_adjust_obs()`.
#' @param log_scale Logical, one per row: is this row fitted on the log scale?
#' @return `hat`, with the log-scale rows converted to medians.
#' @noRd
.hat_as_median <- function(hat, sigma, ba_obs, log_scale) {
  out <- hat
  # No usable sd keeps the reported MEAN rather than wiping the row, since
  # exp(-b * NA / 2) is NA and the figure can still draw a prediction.
  use <- log_scale & is.finite(sigma) & is.finite(hat)
  out[use] <- hat[use] * exp(-ba_obs * sigma[use]^2 / 2)
  out
}
