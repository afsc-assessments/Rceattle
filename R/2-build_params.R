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

  # Initialize from inputs
  for (i in 1:nrow(data_list$M1_base)) {
    sp <- as.numeric(as.character(data_list$M1_base$Species[i]))
    sex <- as.numeric(as.character(data_list$M1_base$Sex[i]))


    # Handle sex == 0 case for 2-sex species
    sex_values <- if (sex == 0) 1:data_list$nsex[sp] else sex

    # Fill in M1 array from fixed values for each sex
    for(j in 1:length(sex_values)){
      m1[sp, sex_values[j], 1:max_age] <- as.numeric(data_list$M1_base[i,(1:max(data_list$nages, na.rm = TRUE)) + 2])
    }
  }
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
    init_vals <- as.numeric(data_list$linkage_table$init)
    init_vals[data_list$linkage_table$design_col == "(Intercept)"] <- 0
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

  # # * 3.2. Kinzery predation function parameters ----
  # param_list$logH_1 = matrix(-8.5, nrow = data_list$nspp, ncol = data_list$nspp + 1)  # Predation functional form
  # param_list$logH_1a = rep(-3, data_list$nspp)  # Age adjustment to H_1; n = [1, nspp]; # FIXME: make matrix
  # param_list$logH_1b = rep(0, data_list$nspp)  # Age adjustment to H_1; n = [1, nspp]; # FIXME: make matrix
  #
  # param_list$logH_2 = matrix(-9, nrow = data_list$nspp, ncol = data_list$nspp)  # Predation functional form; n = [nspp, nspp]
  # param_list$logH_3 = matrix(-9, nrow = data_list$nspp, ncol = data_list$nspp)  # Predation functional form; n = [nspp, nspp]; bounds = LowerBoundH3,UpperBoundH3;
  # param_list$H_4 = matrix(1, nrow = data_list$nspp, ncol = data_list$nspp)  # Predation functional form; n = [nspp, nspp]; bounds = LowerBoundH4,UpperBoundH4;
  #
  #
  # * 3.3. Suitability parameters ----
  param_list$log_gam_a = rep(0.5, data_list$nspp)  # Log predator selectivity;
  names(param_list$log_gam_a) = paste0("Pred: ", data_list$spnames)

  param_list$log_gam_b = rep(-.5, data_list$nspp)  # Log predator selectivity
  names(param_list$log_gam_b) = paste0("Pred: ", data_list$spnames)


  # * 3.4. Preference parameters ----
  param_list$log_phi = matrix(0.5, data_list$nspp, data_list$nspp,
                              dimnames = list(paste0("Pred: ", data_list$spnames), paste0("Prey: ", data_list$spnames)))



  .push_linkage_intercept_inits(param_list, data_list)
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
    intercepts <- lt[lt$design_col == "(Intercept)" & lt$init_supplied &
                       (!fixed_only | as.integer(lt$est_phase) == 0L), , drop = FALSE]
    if (any(is.na(intercepts$init))) {
      stop("Initial value provided for '(Intercept)' is NA.", call. = FALSE)
    }
    for (i in seq_len(nrow(intercepts))) {
      row <- intercepts[i, , drop = FALSE]
      idx <- .linkage_row_indices(row, data_list)
      init_val <- as.numeric(row$init)
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
          .stop_unless_positive(init_val, row$param, "rec_pars")
          par_idx <- .REC_PARAM_TO_INDEX[row$param]
          if (is.na(par_idx)) next
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
