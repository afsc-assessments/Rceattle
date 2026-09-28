

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
  dimnames(quantities$growth_matrix) <- list(
    c(paste(rep(data_list$spnames, each = 2), rep(c("biomass weight", "spawn weight"), data_list$nspp)), data_list$fleet_control$Fleet_name),
    sex_labels, paste0("Age", 1:max_age), paste0("Bin", 1:max_length), yrs_proj)

  # - Fleet
  dimnames(quantities$F_flt_age) <- list(data_list$fleet_control$Fleet_name, sex_labels, paste0("Age", 1:max_age), yrs_proj)
  dimnames(quantities$sel_at_age) <- list(data_list$fleet_control$Fleet_name, sex_labels, paste0("Age", 1:max_age), yrs_proj)
  dimnames(quantities$sel_at_length) <- list(data_list$fleet_control$Fleet_name, sex_labels, paste0("Bin", 1:max_length), yrs_proj)

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
    "Linkage random effects"
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
