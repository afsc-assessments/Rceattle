#' Build parameter list from cpp file
#'
#' @description Construct the TMB parameter list, with every parameter set to its
#'   starting value, for an Rceattle model.
#'
#' @param data_list an Rceattle data_list
#'
#' @return a named list of model parameters at their starting values
#' @export
build_params <- function(data_list) {

  # Fill out switches if missing
  data_list <- Rceattle::switch_check(data_list)

  # - Dimensions
  param_list <- list()

  max_age <- max(data_list$nages, na.rm = TRUE)
  max_sex <- max(data_list$nsex, na.rm = TRUE)
  sex_labels <- c("Sex combined or females", "males")
  if(max_sex == 1){
    sex_labels <- "Sex combined"
  }
  yrs_hind <- data_list$styr:data_list$endyr
  yrs_proj <- data_list$styr:data_list$projyr
  nyrs_hind <- length(yrs_hind)
  nyrs_proj <- length(yrs_proj)


  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # 1. Population dynamics parameters ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

  param_list$dummy = 0  # Variable to test derived quantities given input parameters; n = [1]

  # * 1.0. Population scalar ----
  # Log multiplier on input numbers-at-age (estDynamics = 2); 0 is a multiplier of 1.
  param_list$log_pop_scalar = stats::setNames(rep(0, data_list$nspp), data_list$spnames)

  # * 1.1. Recruitment parameters ----
  # - Stock recruit parameters
  param_list$rec_pars = matrix(9, nrow = data_list$nspp, ncol = 3,
                               dimnames = list(data_list$spnames, c("R0", "Alpha", "Beta")))  # col 1 = mean rec, col 2 = alpha from srr curve, col 3 = beta from srr curve
  param_list$rec_pars[,3] <- log(3) # Starting low here for beta
  param_list$rec_pars[,2] <- 3

  # srr_prior is a prior on alpha for Ricker but on steepness for Beverton-Holt,
  # so it is a valid alpha starting value only in the former case
  # (see .srr_prior_is_alpha).
  if(!is.null(data_list$srr_prior) && .srr_prior_is_alpha(data_list)){
    param_list$rec_pars[,2] <- log(data_list$srr_prior)
  }

  # Explicit starting values from build_srr(). The defaults above carry no
  # information about the stock's scale; see ?build_srr.
  if(!is.null(data_list$srr_alpha_init)){
    param_list$rec_pars[,2] <- log(data_list$srr_alpha_init)
  }
  if(!is.null(data_list$srr_beta_init)){
    param_list$rec_pars[,3] <- log(data_list$srr_beta_init)
  }

  # - Recruitment deviations
  param_list$rec_dev = matrix(0, nrow = data_list$nspp, ncol = nyrs_proj,
                              dimnames = list(data_list$spnames, yrs_proj))  # Annual recruitment deviation; n = [nspp, nyrs_hind]

  param_list$R_log_sd = log(as.numeric(data_list$sigma_rec))  # Standard deviation of recruitment deviations; n = [1, nspp]
  names(param_list$R_log_sd) <- data_list$spnames

  # * 1.3. Initial age-structure parameters ----
  param_list$init_dev = matrix(0, nrow = data_list$nspp, ncol = max_age,
                               dimnames = list(data_list$spnames, paste0("Age", 1:max_age)))

  for(sp in 1:data_list$nspp){

    # Pad every column the cpp never reads with the -999 sentinel. init_dev is
    # indexed as init_dev(sp, age - 1) over age = 1:(nages - 1) by EVERY
    # initMode -- the free-parameter branch, the equilibrium / non-equilibrium
    # branches, and the init_dev penalty alike (ceattle.cpp sections 6.4/6.5 and
    # slot 9) -- so columns nages:max_age are unused for all six modes.
    #
    # The padding is therefore unconditional. Do not gate it on
    # `data_list$initMode`: by this point switch_check() has resolved that to a
    # canonical STRING, so `initMode > 0` is a lexicographic comparison
    # ("FreeParams" > "0" is TRUE) rather than a numeric one, and its result
    # would turn on the spelling of an alias rather than on the mode.
    param_list$init_dev[sp, data_list$nages[sp]:max_age] = -999
  }



  # * 1.4. Natural mortality ----
  # M1 = residual M, M2 = predation M
  # ** Fixed effects ----
  m1 <- array(1, dim = c(data_list$nspp, max_sex, max_age),
              dimnames = list(data_list$spnames, sex_labels, paste0("Age", 1:max_age))) # Set up array
  # Which (species, sex) an M1_base row actually covers
  m1_written <- matrix(FALSE, data_list$nspp, max_sex)

  # Initialize from inputs. The array is dimensioned to the WIDEST species
  # (max_sex, max_age), so a species with fewer sexes or ages leaves padding
  # cells behind. Write only that species' own sexes and ages: a workbook row
  # is legitimately blank past that species' last age, so reading 1:max_age
  # from it would put log(NA) in the padding. The padding keeps the 1 above and
  # reaches the template as log(1) = 0.
  for (i in 1:nrow(data_list$M1_base)) {
    sp <- as.numeric(as.character(data_list$M1_base$Species[i]))
    sex <- as.numeric(as.character(data_list$M1_base$Sex[i]))

    # Handle sex == 0 case for 2-sex species
    sex_values <- if (sex == 0) seq_len(data_list$nsex[sp]) else sex
    ages <- seq_len(data_list$nages[sp])

    # Fill in M1 array from fixed values for each sex
    for (j in seq_along(sex_values)) {
      m1[sp, sex_values[j], ages] <- as.numeric(data_list$M1_base[i, ages + 2])
      m1_written[sp, sex_values[j]] <- TRUE
    }
  }
  # A (species, sex) with no M1_base row keeps the 1 above, which is a residual
  # M of 1.0 per year -- inside the parameter bounds, so nothing downstream
  # catches it. Refuse it here rather than fit it.
  .rce_stop_if_M1_row_missing(m1_written, data_list)
  param_list$log_M1 <- log(m1)


  # ** Age and annual random effects ----
  param_list$log_M1_dev <- array(0, dim = c(data_list$nspp, max_sex, max_age, nyrs_proj),
                                dimnames = list(data_list$spnames, sex_labels, paste0("Age", 1:max_age), yrs_proj)) # Set up array

  # ** M1 fixed parameters ----
  # - Regression coefficients for environment-M1 linkage
  param_list$M1_beta = array(0, dim = c(data_list$nspp, max_sex, ncol(data_list$env_data) - 1),
                             dimnames = list(data_list$spnames, sex_labels, colnames(data_list$env_data)[-1]))

  # - Rho for AR1
  param_list$M1_rho = array(0, dim = c(data_list$nspp, max_sex, 2),
                            dimnames = list(data_list$spnames, sex_labels, c("Age", "Year")))

  # - SD for random effects
  param_list$M1_dev_log_sd = array(0, dim = c(data_list$nspp, max_sex),
                                  dimnames = list(data_list$spnames, sex_labels))


  # * 1.5. fishing mortality ----

  # Future fishing mortality limit
  param_list$log_Flimit = rep(0, data_list$nspp)
  names(param_list$log_Flimit) <- data_list$spnames

  # - Future fishing mortality target
  param_list$log_Ftarget = rep(0, data_list$nspp)
  names(param_list$log_Ftarget) <- data_list$spnames

  # - Initial F when population is not at equilibrium
  param_list$log_Finit = rep(-10, data_list$nspp)
  names(param_list$log_Finit) <- data_list$spnames

  # - Proportion of future fishing mortality for projections for each fleet
  param_list$proj_F_prop = data_list$fleet_control$Proj_F_proportion
  names(param_list$proj_F_prop) <- data_list$fleet_control$Fleet_name

  # - Annual fishing mortality deviations
  param_list$log_F = matrix(0, nrow = nrow(data_list$fleet_control), ncol = nyrs_hind,
                           dimnames = list(data_list$fleet_control$Fleet_name, yrs_hind))

  # -- Make log_F very low if the fleet is turned off or not a fishery
  for (i in 1:nrow(data_list$fleet_control)) {
    if (data_list$fleet_control$Fleet_type[i] != "Fishery") {
      param_list$log_F[i,] <- -999
    }
  }

  # -- Set Fdev for years with 0 catch to very low number
  zero_catch <- data_list$catch_data |>
    dplyr::filter(Year <= data_list$endyr &
                    Catch == 0) |>
    dplyr::mutate(Year = Year - data_list$styr + 1) |>
    dplyr::select(Fleet_code, Year) |>
    as.matrix()
  param_list$log_F[zero_catch] <- -999


  # * 1.6. Growth ----
  # - Mean growth parameters
  param_list$log_growth_pars <- array(0, dim = c(data_list$nspp, max_sex, 4),
                                     dimnames = list(data_list$spnames, sex_labels, c("log_K", "log_L1", "log_Linf", "log_m")))
  # - Initialize
  param_list$log_growth_pars[, , 1] <- log(0.3)

  if(nrow(data_list$caal_data) > 0){
    caal_lengths <- data_list$caal_data |>
      dplyr::distinct(Species, Length) |>
      dplyr::arrange(Species, Length) |>
      dplyr::group_by(Species) |>
      dplyr::mutate(Bin = paste0("Bin", 1:n())) |>
      dplyr::slice(c(1, n())) |>
      dplyr::ungroup() |>
      tidyr::pivot_wider(names_from = Bin, values_from = Length)
    param_list$log_growth_pars[caal_lengths$Species, 1, 2:3] <- as.matrix(log(caal_lengths[,-1]))
    if(max_sex == 2){
      param_list$log_growth_pars[caal_lengths$Species, 2, 2:3] <- as.matrix(log(caal_lengths[,-1]))
    }
  }

  param_list$growth_log_sd <- array(0, dim = c(data_list$nspp, max_sex, 2),
                                   dimnames = list(data_list$spnames, sex_labels, c("log_sd_minage", "log_sd_maxage")))
  # Endpoints read as CVs (sd_form = "CV") start at a CV of 0.1; as SDs, at 1 cm.
  cv_sp <- which(rep_len(data_list$growth_sd_form %||% 1L, data_list$nspp) == 2L)
  param_list$growth_log_sd[cv_sp, , ] <- log(0.1)
  param_list$weight_length_pars <- matrix(0, nrow = data_list$nspp, ncol = 2,
                                          dimnames = list(data_list$spnames, c("a", "b")))  # Weight-length parameters
  param_list$weight_length_pars[,1] <- data_list$alpha_wt_len
  param_list$weight_length_pars[,2] <- data_list$beta_wt_len

  # * 1.3b. Linkage-table coefficients ----
  # Aligned row-for-row with the pooled `data_list$linkage_table` (set
  # by `pool_linkages()` inside `fit_mod`). Initial values come from
  # the `init` column of the table; (Intercept) rows are forced to 0
  # because they're mapped out of estimation -- the base parameter
  # carries the level instead (see 1.3c below). Absent table =>
  # length-0 vector.
  if (!is.null(data_list$linkage_table) &&
      nrow(data_list$linkage_table) > 0L) {
    lt <- data_list[["linkage_table"]]
    init_vals <- as.numeric(lt[["init"]])
    init_vals[.is_pinned_intercept(lt)] <- 0
    # An `R_init` intercept keeps its own starting value: it has no base
    # parameter to re-target, so nothing else holds the level. Given on the
    # natural scale -- a multiplier on R0, as every other intercept's `init` is
    # natural-scale -- and stored logged, because the level is added inside the
    # exp() that builds the initial numbers-at-age.
    for (r in which(.is_level_intercept(lt))) {
      init_vals[r] <- .r_init_log_start(init_vals[r], lt[["init_supplied"]][r])
    }
    param_list$beta_linkage <- init_vals
  } else {
    param_list$beta_linkage <- numeric(0)
  }

  # Random-effect linkage machinery, sized from the RE registry that
  # pool_linkages() wrote onto the table. beta_linkage_re holds the deviation
  # coefficients that enter the Laplace approximation (one per RE row, indexed
  # by `re_index`); log_sigma_linkage is one log-SD per RE group (`sigma_index`);
  # trans_rho_linkage is one transformed correlation per autocorrelated (ar1)
  # group. All length-0 until a random linkage spec is supplied, so a model
  # without one is numerically unchanged.
  if (!is.null(data_list$linkage_table) &&
      nrow(data_list$linkage_table) > 0L &&
      any(!is.na(data_list$linkage_table$re_index))) {
    lt      <- data_list$linkage_table
    gt      <- .re_group_table(lt)
    n_re    <- sum(!is.na(lt$re_index))
    # Split the slot space across the two destination vectors. Only
    # beta_linkage_re joins TMB's `random` set; beta_linkage_re_pen stays a fixed
    # effect whose density is a plain penalty. Both are sized from the same
    # routing the encoder and build_map use.
    rt      <- .re_slot_routing(lt)
    # ar1 groups get a rho; us/rw groups do not. (rho estimation lands with
    # ar1() -- until then no group is ar1 and this stays length 0.)
    n_ar1   <- sum(gt$re_struct == "ar1")
    param_list$beta_linkage_re     <- numeric(sum(rt$integrate))   # integrated, init 0
    param_list$beta_linkage_re_pen <- numeric(sum(!rt$integrate))  # penalized, init 0
    # One log-SD per group; start from linkage_spec(init = list(sigma = )) when
    # supplied (fixed there via build_map), else a default. gt is ordered by
    # sigma_index so element g is group g - 1.
    param_list$log_sigma_linkage <- log(gt$sigma_start)
    # One transformed correlation per ar1 group (ar1 groups in gt order, matching
    # linkage_re_rho). trans is the rho_trans pre-image atanh(rho) of the start
    # correlation (default 0). rho_trans(x) = 2/(1+exp(-2x)) - 1, so x=atanh(rho).
    if (n_ar1 > 0L) {
      rho_start <- gt$rho_start[gt$re_struct == "ar1"]
      param_list$trans_rho_linkage <- atanh(pmin(pmax(rho_start, -0.999), 0.999))
    } else {
      param_list$trans_rho_linkage <- numeric(0)
    }
    # Rogers QAR1 effect size: one estimated beta per observed group (the latent
    # ar1 deviate enters the target as beta * deviate). Init 0 (no effect).
    param_list$beta_linkage_obs <- numeric(sum(gt$observed))
    # Rogers QAR1 observation SD: one log-SD per observed group, estimated by
    # default (matching the reference Estimate_q = 6 / GOApollock, which estimates
    # the measurement SD), started from the `obs_sd` value supplied on the spec.
    param_list$log_obs_sd_linkage <- log(gt$obs_sd[gt$observed])
  } else {
    param_list$beta_linkage_re     <- numeric(0)
    param_list$beta_linkage_re_pen <- numeric(0)
    param_list$log_sigma_linkage <- numeric(0)
    param_list$trans_rho_linkage <- numeric(0)
    param_list$beta_linkage_obs  <- numeric(0)
    param_list$log_obs_sd_linkage <- numeric(0)
  }

  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # 2. Observation model parameters ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

  # * 2.1. Catchability parameters ----
  # - Catchability on log scale
  param_list$index_log_q = log(data_list$fleet_control$Catchability_init)
  names(param_list$index_log_q) <- data_list$fleet_control$Fleet_name

  # - Regression coefficients for environment-q linkage
  param_list$index_q_beta = matrix(0, nrow = nrow(data_list$fleet_control), ncol = ncol(data_list$env_data) - 1,
                                   dimnames = list(data_list$fleet_control$Fleet_name, colnames(data_list$env_data)[-1]))

  # param_list$index_q_pow = rep(0, nrow(data_list$fleet_control))

  # - Annual index catchability deviations
  param_list$index_q_dev = matrix(0, nrow = nrow(data_list$fleet_control), ncol = nyrs_hind,
                                  dimnames = list(data_list$fleet_control$Fleet_name, yrs_hind))

  # - Log standard deviation prior on Q (maybe should be data...)
  param_list$index_q_log_sd <- log(data_list$fleet_control$Catchability_prior_sd)
  names(param_list$index_q_log_sd) <- data_list$fleet_control$Fleet_name

  # - Log SD of the time-varying survey catchability deviations
  param_list$index_q_dev_log_sd <- log(data_list$fleet_control$Time_varying_q_sd)
  names(param_list$index_q_dev_log_sd) <- data_list$fleet_control$Fleet_name


  # * 2.2. Selectivity parameters ----
  n_selectivities <- nrow(data_list$fleet_control)
  max_sel_bins <- max(c(1, as.numeric(data_list$fleet_control$N_sel_bins)), na.rm = TRUE)

  # - Non-parametric selectivity coefficients
  param_list$sel_coff =  array(0, dim = c(n_selectivities, max_sex, max_sel_bins),
                               dimnames = list(data_list$fleet_control$Fleet_name, sex_labels, paste0("Bin", 1:max_sel_bins)))

  # - Non-parametric selectivity penalties (sensu Ianelli / ADMB AMAK):
  #   col 1 = decreasing penalty weight (ADMB ctrl_flag(13));
  #   col 2 = curvature (2nd-difference) weight (ADMB ctrl_flag(11)/nch);
  #   col 3 = dev-magnitude weight on coefficient increments (ADMB ctrl_flag(10)/group_num).
  param_list$sel_curve_pen = matrix( c(data_list$fleet_control$Sel_curve_pen1, data_list$fleet_control$Sel_curve_pen2, data_list$fleet_control$Sel_curve_pen3), nrow = n_selectivities, ncol = 3)
  param_list$sel_curve_pen[is.na(param_list$sel_curve_pen)] <- 0

  # - Non-parametric selectivity coef annual deviates
  param_list$sel_coff_dev = array(0, dim = c(n_selectivities, max_sex, max_sel_bins, nyrs_hind),
                                  dimnames = list(data_list$fleet_control$Fleet_name, sex_labels, paste0("Bin", 1:max_sel_bins), yrs_hind))

  # - Selectivity slope parameters for logistic
  param_list$log_sel_slp = array(0.5, dim = c(2, n_selectivities, max_sex),
                                dimnames = list(c("Ascending" , "Descending"), data_list$fleet_control$Fleet_name, sex_labels))

  # - Selectivity asymptotic parameters for logistic
  param_list$sel_inf = array(0, dim = c(2, n_selectivities, max_sex),
                             dimnames = list(c("Ascending" , "Descending"), data_list$fleet_control$Fleet_name, sex_labels))
  param_list$sel_inf[1,,] <- 0
  param_list$sel_inf[2,,] <- 10

  # For length-based selectivity the ascending parameter sel_inf[1] is an
  # inflection *length*, so the age-scale default of 0 starts below the
  # smallest length bin; from there the optimizer can wander to nonsensical
  # (even negative) lengths. Start it near the middle of the species' length
  # range instead. Age-based selectivity is left at 0 (unchanged). The length
  # scale used here mirrors data_list$lengths as built in rearrange_data():
  # physical bin centres from caal_data when present, else 1:nlengths indices.
  sel_dim <- data_list$fleet_control$Selectivity_dimension
  if (!is.null(sel_dim)) {
    length_midpoint <- function(sp) {
      cd <- data_list$caal_data
      if (!is.null(cd) && nrow(cd) > 0 && all(c("Species", "Length") %in% names(cd))) {
        len <- cd$Length[cd$Species == sp]
        len <- len[is.finite(len)]
        if (length(len) > 0) return(mean(range(len)))
      }
      nl <- data_list$nlengths[min(sp, length(data_list$nlengths))]
      (1 + nl) / 2
    }
    for (flt in which(!is.na(sel_dim) & tolower(sel_dim) == "length")) {
      param_list$sel_inf[1, flt, ] <- length_midpoint(data_list$fleet_control$Species[flt])
    }
  }

  # LogisticPM (type 11) reuses the unused descending-limb slot sel_inf[2] as the
  # free first-bin (age-1) LOG-selectivity (AMAK sel_age_one_bts), not an
  # inflection age. The shared default of 10 would start age-1 selectivity at
  # exp(10) = 22026, swamping the index prediction, so start at 0 (selectivity 1).
  sel_type <- data_list$fleet_control$Selectivity
  if (!is.null(sel_type)) {
    logisticpm_flts <- which(sel_type %in% c(11, "11", "LogisticPM"))
    if (length(logisticpm_flts) > 0) {
      param_list$sel_inf[2, logisticpm_flts, ] <- 0
    }
  }

  # - DoubleNormalSS3 (SS3 size pattern 24): six parameters on SS3's own
  #   scales -- peak (cm or age), logit top width, log ascending and descending
  #   widths, logit initial and final selectivity. The start is a dome with
  #   neither end scaled (SS3's -999 for P5 and P6), so only the first four are
  #   estimated unless the ends are given values.
  param_list$sel_dn6 <- array(0, dim = c(6, n_selectivities, max_sex),
                              dimnames = list(c("peak", "top_logit", "ascend_se", "descend_se", "start_logit", "end_logit"),
                                              data_list$fleet_control$Fleet_name, sex_labels))
  dn6_flts <- which(data_list$fleet_control$Selectivity %in% c(15, "15", "DoubleNormalSS3"))
  for (flt in dn6_flts) {
    is_len <- isTRUE(tolower(data_list$fleet_control$Selectivity_dimension[flt]) == "length")
    sp <- data_list$fleet_control$Species[flt]
    param_list$sel_dn6[, flt, ] <- c(
      if (is_len) length_midpoint(sp) else (data_list$nages[sp] + 1) / 2,  # peak
      -5,                                        # top width, logit: narrow
      if (is_len) 5 else log(4),                 # ascending width, log
      if (is_len) 5 else log(4),                 # descending width, log
      -999, -999)                                # ends unscaled
  }

  # - Annual selectivity slope deviation for logistic
  param_list$log_sel_slp_dev = array(0, dim = c(2, n_selectivities, max_sex, nyrs_hind),
                                    dimnames = list(c("Ascending" , "Descending"), data_list$fleet_control$Fleet_name, sex_labels, yrs_hind))

  # - Annual selectivity asymptotic deviations for logistic
  param_list$sel_inf_dev = array(0, dim = c(2, n_selectivities, max_sex, nyrs_hind),
                                 dimnames = list(c("Ascending" , "Descending"), data_list$fleet_control$Fleet_name, sex_labels, yrs_hind))

  # - Per-sex apical height, log scale: the whole curve of one sex times
  #   exp(log_sel_apical), applied after the form and before normalization. 0 is
  #   no offset; estimated only through a selectivity linkage on `apical`.
  param_list$log_sel_apical = array(0, dim = c(n_selectivities, max_sex),
                                    dimnames = list(data_list$fleet_control$Fleet_name, sex_labels))

  # - Log standard deviation for selectivity random walk - used for logistic
  param_list$sel_dev_log_sd <- log(data_list$fleet_control$Time_varying_sel_sd)
  names(param_list$sel_dev_log_sd) <- data_list$fleet_control$Fleet_name


  # * 2.3. Variance of survey and fishery time series ----
  # - Log standard deviation of survey index time-series
  param_list$index_log_sd = log(data_list$fleet_control$Index_sd)
  names(param_list$index_log_sd) <- data_list$fleet_control$Fleet_name

  # - Log standard deviation of fishery catch time-series
  param_list$catch_log_sd = log(data_list$fleet_control$Catch_sd)
  names(param_list$catch_log_sd) <- data_list$fleet_control$Fleet_name

  # * 2.4. Comp weighting ----
  param_list$comp_weights = data_list$fleet_control$Comp_weights
  names(param_list$comp_weights) <- data_list$fleet_control$Fleet_name

  # * 2.5. CAAL weighting ----
  param_list$caal_weights = data_list$fleet_control$CAAL_weights
  names(param_list$caal_weights) <- data_list$fleet_control$Fleet_name

  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  # 3. Predation model parameters ----
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#


  # * 3.1. Diet composition weighting ----
  param_list$diet_comp_weights = data_list$Diet_comp_weights
  names(param_list$diet_comp_weights) <- paste0("Pred: ", data_list$spnames)

  # * 3.2. Suitability parameters ----
  param_list$log_gam_a = rep(0.5, data_list$nspp)  # Log predator selectivity;
  names(param_list$log_gam_a) = paste0("Pred: ", data_list$spnames)

  param_list$log_gam_b = rep(-.5, data_list$nspp)  # Log predator selectivity
  names(param_list$log_gam_b) = paste0("Pred: ", data_list$spnames)


  # * 3.3. Preference parameters ----
  param_list$log_phi = matrix(0.5, data_list$nspp, data_list$nspp,
                              dimnames = list(paste0("Pred: ", data_list$spnames), paste0("Prey: ", data_list$spnames)))



  param_list <- .push_linkage_intercept_inits(param_list, data_list)

  # Checked HERE, not where log_M1 is filled: an M1 linkage carrying an
  # `init` for its intercept writes the level over every real age above, so a
  # workbook may legitimately leave M1_base blank and let the linkage supply
  # it. Checking at the fill would refuse that.
  .rce_stop_if_nonfinite_M1(param_list$log_M1, "build_params()")

  param_list
}


# Push linkage (Intercept) inits onto the base parameters. `fixed_only = TRUE` pushes only
# est_phase = 0 rows; fit_mod() re-applies those over supplied `inits`, so a fixed value wins.
.push_linkage_intercept_inits <- function(param_list, data_list, fixed_only = FALSE) {
  # * Push (Intercept) inits to the base parameter ----
  # An intercept-bearing linkage formula (`~ 1`, `~ temp`, ...) emits an
  # "(Intercept)" row whose coefficient stays fixed at 0 (mapped out by
  # `build_map_linkages()`); the base parameter carries the level and remains
  # estimable. An `init` supplied for the intercept column therefore sets the
  # BASE parameter's starting value, so phasing and priors operate from a
  # sensible point. Without `init`, the base keeps its `build_params()` default.
  #
  # This runs at the end of build_params() because the base parameters it writes
  # span the whole function -- log_M1 and rec_pars are built in section 1,
  # index_log_q / sel_* / *_weights in section 2 -- and every process must be
  # able to reach its own.
  #
  # Values are supplied on the parameter's NATURAL scale. Base parameters that
  # are stored logged (log_M1, rec_pars, log_growth_pars, index_log_q, the DM
  # weights, log_sel_slp) take log(init); ones stored as-is (sel_inf) take it
  # directly.
  #
  # Slope-only formulas (`~ 0 + temp`) emit no (Intercept) row and are handled
  # entirely by `map_linkage_adjuster()`, which maps the base parameter out so
  # it stays at its default.
  if (!is.null(data_list$linkage_table) &&
      nrow(data_list$linkage_table) > 0L) {
    lt <- data_list$linkage_table
    fixed_now <- as.integer(lt$est_phase) == 0L
    # An `R_init` level fixed at phase 0 is re-pushed whether or not an `init`
    # was given: unlike a base-parameter intercept, which falls back on its own
    # build_params() default, the level's fixed value IS this row's `init`, and
    # the table default of 0 is a well-defined multiplier of 1. Gating it on
    # `init_supplied` left "pinned at no shift" to be overwritten by a warm
    # start, which is the same silent hold this push exists to prevent.
    int_rows <- which(lt$design_col == "(Intercept)" &
                        (lt$init_supplied | (.is_level_intercept(lt) & fixed_now)) &
                        (!fixed_only | fixed_now))
    if (any(is.na(lt$init[int_rows]))) {
      stop("Initial value provided for '(Intercept)' is NA.", call. = FALSE)
    }
    for (ri in int_rows) {
      row <- lt[ri, , drop = FALSE]
      init_val <- as.numeric(row$init)
      # `R_init` is the one intercept whose level lives in beta_linkage rather
      # than in a base parameter, so this is where a fixed (est_phase = 0)
      # level is re-applied over supplied `inits`. Without it a warm start
      # would overwrite the level and the mapped-out row would hold the wrong
      # value for the whole fit. Taken before .linkage_row_indices(), which
      # resolves a base parameter this row does not have.
      if (.is_level_intercept(row)) {
        param_list$beta_linkage[ri] <-
          .r_init_log_start(init_val, row[["init_supplied"]])
        next
      }
      idx <- .linkage_row_indices(row, data_list)
      switch(row$process,
        growth = {
          .stop_unless_positive(init_val, row$param, "the growth parameter")
          # Mean-growth params live on log_growth_pars[sp, sex, k];
          # SD endpoints live on growth_log_sd[sp, sex, k']. Both expose
          # the same intercept-init contract -- dispatch on the param
          # name once and write into the right tensor.
          mean_idx <- .GROWTH_PARAM_TO_INDEX[row$param]
          sd_idx   <- .GROWTH_SD_PARAM_TO_INDEX[row$param]
          if (!is.na(mean_idx)) {
            for (s in idx$species) {
              param_list$log_growth_pars[s, idx$per_sp[[as.character(s)]]$sex, mean_idx] <- log(init_val)
            }
          } else if (!is.na(sd_idx)) {
            for (s in idx$species) {
              param_list$growth_log_sd[s, idx$per_sp[[as.character(s)]]$sex, sd_idx] <- log(init_val)
            }
          }
        },
        M = {
          .stop_unless_positive(init_val, row$param, "log_M1")
          for (s in idx$species) {
            sx <- idx$per_sp[[as.character(s)]]$sex
            ag <- idx$per_sp[[as.character(s)]]$age
            param_list$log_M1[s, sx, ag] <- log(init_val)
          }
        },
        recruitment = {
          # `R_init` has no rec_pars column -- its starting value was logged
          # onto beta_linkage above, so there is nothing to push to a base
          # parameter and nothing here to check against one.
          par_idx <- .REC_PARAM_TO_INDEX[row$param]
          if (is.na(par_idx)) next
          .stop_unless_positive(init_val, row$param, "rec_pars")
          param_list$rec_pars[idx$species, par_idx] <- log(init_val)
        },
        q = {
          # index_log_q is log catchability, one entry per fleet.
          .stop_if_shared_block(data_list, row$fleet, "q", row$param)
          .stop_unless_positive(init_val, row$param, "index_log_q")
          param_list$index_log_q[idx$fleet] <- log(init_val)
        },
        comp = {
          # The Dirichlet-multinomial weight is exp(<weight parameter>), so an
          # init given as the weight itself is stored logged. theta_comp and
          # theta_caal are fleet-indexed, theta_diet predator(species)-indexed.
          #
          # The log is correct only for the DM reading of this slot: a
          # multinomial likelihood uses comp_weights raw, as a constant
          # multiplier on the sample size (ceattle.cpp JNLL_COMP). That case
          # cannot arrive here because .check_comp_linkage_support() rejects a
          # comp linkage on a non-DM fleet -- if that ever relaxes, this needs
          # to branch on the fleet's Comp_distribution.
          .stop_unless_positive(init_val, row$param, "the DM weight")
          switch(row$param,
            theta_comp = { param_list$comp_weights[idx$fleet]      <- log(init_val) },
            theta_caal = { param_list$caal_weights[idx$fleet]      <- log(init_val) },
            theta_diet = { param_list$diet_comp_weights[idx$species] <- log(init_val) })
        },
        sel = {
          slot <- .SEL_PARAM_TO_SLOT[[row$param]]
          if (is.null(slot)) next
          .stop_if_shared_block(data_list, row$fleet, "sel", row$param)
          if (identical(slot$arr, "log_sel_slp")) {
            .stop_unless_positive(init_val, row$param, "log_sel_slp")
            for (s in idx$species) {
              param_list$log_sel_slp[slot$slot, idx$fleet,
                                     idx$per_sp[[as.character(s)]]$sex] <- log(init_val)
            }
          } else if (identical(slot$arr, "sel_inf")) {
            # Slot 2 is an inflection only for the logistic family; DoubleNormal
            # and LogisticPM store a logit and a log there instead.
            for (f in idx$fleet) {
              .stop_unless_natural_sel_inf(data_list, f, slot$slot, row$param)
            }
            for (s in idx$species) {
              param_list$sel_inf[slot$slot, idx$fleet,
                                 idx$per_sp[[as.character(s)]]$sex] <- init_val
            }
          } else if (identical(slot$arr, "sel_dn6")) {
            # DoubleNormalSS3 holds each parameter on SS3's scale, so the init
            # is the value as an SS3 control file gives it.
            for (s in idx$species) {
              param_list$sel_dn6[slot$slot, idx$fleet,
                                 idx$per_sp[[as.character(s)]]$sex] <- init_val
            }
          } else if (identical(slot$arr, "log_sel_apical")) {
            # The init is the multiplier itself (1 = no offset), stored logged.
            .stop_unless_positive(init_val, row$param, "log_sel_apical")
            for (s in idx$species) {
              param_list$log_sel_apical[idx$fleet,
                                        idx$per_sp[[as.character(s)]]$sex] <- log(init_val)
            }
          }
          # `coff` is a per-bin vector with no single level to set, so an
          # intercept init has no well-defined target; rejected up front by
          # .check_sel_linkage_support() rather than guessed at here.
        }
      )
    }
  }

  return(param_list)
}


#' Refuse a non-finite starting value for log_M1
#'
#' @description
#' `log_M1` is dimensioned to the widest species, so a species with fewer sexes
#' or ages leaves padding cells that no `nages(sp)`-bounded loop in the template
#' reads. `NA` or `-Inf` there is never an intended natural mortality: the two
#' fill paths disagreed on the value through 5.49.4, a loop bound that widened
#' to the padding would read `exp(-Inf)` as an M1 of 0, and under
#' `M1_model >= 1` the padding sex cells are estimated rather than constant.
#'
#' @param log_M1 The filled `log_M1` array.
#' @param where Name of the calling path, for the message.
#' @noRd
.rce_stop_if_nonfinite_M1 <- function(log_M1, where) {
  bad <- which(!is.finite(log_M1), arr.ind = TRUE)
  if (!nrow(bad)) return(invisible(TRUE))
  stop(where, ": log_M1 holds ", nrow(bad),
       " non-finite starting value(s), at (species, sex, age) ",
       paste(apply(utils::head(bad, 5), 1, paste, collapse = ","),
             collapse = "; "),
       if (nrow(bad) > 5) ", ..." else "",
       ". M1_base must give a finite value for every age of every sex the ",
       "species has.", call. = FALSE)
}


#' Refuse an `M1_base` that does not cover every sex of every species
#'
#' @description
#' A `(species, sex)` the model has but `M1_base` does not name keeps the
#' array's initial 1, i.e. a residual M of 1.0 per year. That sits inside
#' `build_bounds()`'s `[log(0.001), log(2)]`, so nothing downstream catches it
#' and the fit simply runs with it.
#'
#' @param written `[nspp, max_sex]` logical: which cells an `M1_base` row set.
#' @param data_list The `data_list`, for `nsex` and `spnames`.
#' @noRd
.rce_stop_if_M1_row_missing <- function(written, data_list) {
  gaps <- character()
  for (sp in seq_len(data_list$nspp)) {
    for (sx in seq_len(data_list$nsex[sp])) {
      if (!written[sp, sx]) {
        gaps <- c(gaps, paste0(data_list$spnames[sp], " (species ", sp,
                               ", sex ", sx, ")"))
      }
    }
  }
  if (!length(gaps)) return(invisible(TRUE))
  stop("M1_base has no row for ", paste(gaps, collapse = "; "),
       ". Every sex of every species needs one, or its residual M starts at ",
       "1.0 per year, which is inside the parameter bounds and so is fit ",
       "rather than refused.", call. = FALSE)
}
