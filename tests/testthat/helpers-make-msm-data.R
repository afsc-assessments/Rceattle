# Adapted from Matt Cheng
#
#' Make a multi-species operating-model data set
#'
#' Simulates an age-structured, multi-species CEATTLE population (MSVPA-style
#' predation) and returns both the underlying truth and an Rceattle-shaped data
#' list, so estimation tests can compare a fitted model against known values.
#'
#' The population is projected forward `niter` times so predation mortality
#' (`M2`) converges, then observations (survey index, catch, age/length comps,
#' diet proportions) are drawn with the given sampling error / input sample
#' sizes. Selectivity, growth, suitability, and bioenergetics are all
#' configurable via the formals.
#'
#' @param nspp,years,ages Dimensions of the simulated system.
#' @param mean_Rec,sigma_R,M,rhoM,sigmaM Recruitment level/variability and
#'   natural-mortality level plus its AR(1) time-variation.
#' @param fish_sel,srv_sel,fish_size_sel,srv_size_sel,use_size_sel Age- or
#'   length-based selectivity for the fishery and survey.
#' @param growth_params,growth_model,nlengths,lengths Growth (K, L1, Linf, m)
#'   and the length bins used to build the age-length transition.
#' @param gam_a,gam_b,log_phi,other_food,ration_mult,normalize_suitability
#'   Predation/diet controls (size-preference, vulnerability, ration).
#' @param sigma_catch,sigma_srv,fish_ISS,srv_ISS,diet_ISS,fish_CAAL_ISS,srv_CAAL_ISS
#'   Observation error and input sample sizes for each data type.
#'
#' @return A list with `model_quantities` (the simulated truth: NAA, FAA, SSB,
#'   suitability, diet_prop, ...), `data_list` (the Rceattle data list ready for
#'   `fit_mod()`), and `base_data` (the data list before ration is added).
make_msm_test_data <- function(
    nspp = 2,
    years = 1:30,
    ages = 1:15,
    MatAA = matrix(c(1 / (1 + exp(-1 * (ages - 5))),
                     1 / (1 + exp(-0.8 * (ages - 5)))),
                   nspp, length(ages), byrow = TRUE),
    mean_Rec = c(1e2, 1e3),
    sigma_R = 1,
    sigma_catch = 0.001,
    sigma_srv = 0.05,
    diet_ISS = 1e5,
    fish_ISS = 1e5,
    srv_ISS = 1e5,
    M = c(0.2, 0.3),
    rhoM = 0,
    sigmaM = 0,
    fish_sel = matrix(c(1 / (1 + exp(-2.5 * (ages - 6))),
                        1 / (1 + exp(-2.5 * (ages - 4)))), nspp, length(ages), byrow = TRUE),
    srv_sel = matrix(c(1 / (1 + exp(-2 * (ages - 3))),
                       1 / (1 + exp(-2 * (ages - 2.5)))), nspp, length(ages), byrow = TRUE),
    Fmort = matrix(c(seq(0.02, 0.3, length.out = length(years)/2), seq(0.3, 0.05, length.out = length(years)/2),
                     seq(0.02, 0.3, length.out = length(years))),
                   nrow = nspp, ncol = length(years), byrow = TRUE),
    srv_q = rep(1, nspp),

    # Growth bits
    growth_params = matrix(c(0.3, 4.5, 90, 1.0,
                             0.35, 4.5, 50, 1.0), # K, L1, Linf, M
                           nrow = nspp, ncol = 4, byrow = TRUE),
    growth_model = 1,
    nlengths = 20,
    lengths = 5 * 1:nlengths,
    fish_CAAL_ISS = 0,
    srv_CAAL_ISS = 0,
    fish_size_sel = matrix(c(1 / (1 + exp(-2.5 * (lengths + 2.5 - 35))), # Sp 1
                             1 / (1 + exp(-2.5 * (lengths + 2.5 - 30)))), # Sp 2
                           nrow = nspp, ncol = nlengths, byrow = TRUE),
    srv_size_sel = matrix(c(1 / (1 + exp(-2 * (lengths + 2.5 - 20))), # +2.5 for bin midpoint
                            1 / (1 + exp(-2 * (lengths + 2.5 - 15)))),
                          nrow = nspp, ncol = nlengths, byrow = TRUE),
    use_size_sel = FALSE,

    # Multispecies bits
    niter = 5,
    gam_a = c(1, 0.1),
    gam_b =  rep(0.15, nspp),
    log_phi =  matrix(c(-5,0.5,-10,-2), nspp, nspp, byrow = TRUE),
    other_food =  rep(1e5, nspp),
    ration_mult = 50, # Multiplied against WAA for ration
    normalize_suitability = FALSE # NEW: Flag to force suitability to sum to 1.0

) {

  nyrs <- length(years)
  nages <- length(ages)

  # Initialize arrays ----
  NAA <- array(0, dim = c(nspp, nages, nyrs))  # Numbers at age
  avgNAA <- array(0, dim = c(nspp, nages, nyrs))  # Numbers at age
  ZAA <- array(0, dim = c(nspp, nages, nyrs))      # Total mortality
  MAA <- array(0, dim = c(nspp, nages, nyrs))      # Total mortality
  M2_at_age <- array(0, dim = c(nspp, nages, nyrs))
  CAA <- array(0, dim = c(nspp, nages, nyrs))      # Catch at age
  FAA <- array(0, dim = c(nspp, nages, nyrs))      # Fishing mortality at age

  # Population metrics
  SSB <- matrix(0, nspp, nyrs)
  Total_Biom <- matrix(0, nspp, nyrs)
  Catch <- matrix(0, nspp, nyrs)

  # Generate M deviations
  M_vec <- suppressWarnings(
    matrix(arima.sim(n = nyrs * nspp, list(order=c(1,0,0), ar=rhoM)
                     , sd = sigmaM), nspp, nyrs))
  Mtv = M * exp(M_vec)

  # Growth ----
  # - Growth params: (sp, sex, yr, param) -> params: K, L1, Linf, m
  WAA <- matrix(NA, nrow = nspp, ncol = nages)
  growth_matrix <- array(0, dim = c(nspp, nages, nlengths))
  length_at_age <- array(0, dim = c(nspp, nages))
  for(sp in 1:nspp){
    gp <- array(growth_params[sp,], dim=c(1, 1, 4))
    gsd <- array(log(3), dim=c(nspp, 2))

    # - Run Growth Matrix
    gm <-  Rceattle:::get_growth_matrix_r(fracyr=0,
                              nsex_sp=1,
                              nages_sp=nages,
                              nlengths_sp = nlengths,
                              nyrs = 1,
                              lengths_sp= lengths,
                              minage_sp= 1,
                              maxage_sp = nages,
                              growth_params_sp = gp,
                              growth_log_sd_sp = gsd,
                              growth_model_sp = growth_model,
                              M_plus_sp = M[sp])
    growth_matrix[sp, ,] = gm$growth_matrix
    length_at_age[sp, ] = gm$length_at_age

    # 4. Calculate Weight-at-Age
    lw_p <- array(c(0.00001, 3.0), dim=c(1, 1, 2)) # a=0.00001, b=3.0
    WAA[sp,] <- Rceattle:::get_weight_at_age_r(nsex_sp = 1, nages_sp = nages, nlengths_sp = nlengths, nyrs = 1,
                                    lengths_sp = lengths, growth_matrix = gm$growth_matrix, lw_params = lw_p)

    # Do any growth transformations
    if(use_size_sel){
      fish_sel[sp,] <- growth_matrix[sp, ,] %*% fish_size_sel[sp,]
      srv_sel[sp,] <- growth_matrix[sp, ,] %*% srv_size_sel[sp,]
    }
  }

  # Ration
  ration = ration_mult * WAA

  # Generate recruitment deviations
  rec_devs <- matrix(stats::rnorm(nyrs * nspp, 0, sigma_R), nspp, nyrs)
  init_devs <- matrix(stats::rnorm((nages-1) * nspp, 0, sigma_R), nspp, nages-1)

  # Vulnerability ----
  suitability = array(0, dim = c(nspp, nspp, nages, nages))
  vulnerability <- matrix(0, nspp, nspp)
  suit_other <- c()

  for(sp in 1:nspp) {
    vulnerability[sp,] = exp(log_phi[sp,])/(1+sum(exp(log_phi[sp,]))) # multinomial logistic transformation
    suit_other[sp] = 1 - sum(vulnerability[sp,])
  }

  # Calculate suitability ----
  for(sp in 1:nspp) { # Predator
    for(r_age in 1:nages){  # Pred age
      for(ksp in 1:nspp) {   # Prey loop
        for(k_age in 1:nages){ # Prey age

          # Weight-based lognormal suitability
          log_size_ratio = log(WAA[sp,r_age] / WAA[ksp,k_age]) # Log ratio of weights

          if(log_size_ratio > 0){
            suitability[sp, ksp, r_age, k_age] = vulnerability[sp, ksp] * dnorm(log_size_ratio, gam_a[sp], gam_b[sp]) / dnorm(gam_a[sp], gam_a[sp], gam_b[sp]) # Divide by mode to scale to 1
          }
        }
      }
    }
  }

  # NEW: Normalize suitability to sum to 1.0 across all prey for MSVPA testing
  suit_other_mat <- matrix(0, nspp, nages)
  for(sp in 1:nspp) suit_other_mat[sp, ] <- suit_other[sp]

  if(normalize_suitability) {
    for(sp in 1:nspp) {
      for(r_age in 1:nages) {
        S_tot = sum(suitability[sp, , r_age, ]) + suit_other_mat[sp, r_age]
        if(S_tot > 0) {
          suitability[sp, , r_age, ] = suitability[sp, , r_age, ] / S_tot
          suit_other_mat[sp, r_age] = suit_other_mat[sp, r_age] / S_tot
        }
      }
    }
  }

  # Initialize population
  for(sp in 1:nspp){
    init_age_idx <- 1:(nages - 2)
    NAA[sp, init_age_idx + 1, 1] <- mean_Rec[sp] * exp(- (init_age_idx * Mtv[sp, 1]))
    NAA[sp, nages, 1] <- mean_Rec[sp] * exp(-(nages - 1) * Mtv[sp, 1]) / (1 - exp(-Mtv[sp, 1]))
    NAA[sp,2:nages, 1] <- NAA[sp,2:nages, 1] * exp(init_devs[sp,])
  }

  # Project population forward
  for(iter in 1:niter){

    for(y in 1:nyrs) {
      for(sp in 1:nspp){
        # New recruits
        NAA[sp, 1, y] <- mean_Rec[sp] * exp(rec_devs[sp,y])

        # Calculate mortality
        FAA[sp, ,y] <- Fmort[sp, y] * fish_sel[sp,]
        ZAA[sp, ,y] <- FAA[sp,,y] + Mtv[sp, y] + M2_at_age[sp, , y]
        MAA[sp, ,y] <- Mtv[sp, y] + M2_at_age[sp, ,y]

        # Calculate catch
        CAA[sp, ,y] <- FAA[sp, ,y] / ZAA[sp, ,y] * NAA[sp, ,y] * (1 - exp(-ZAA[sp, ,y]))

        # Project survivors
        if(y < nyrs) {
          for(a in 1:(nages-1)) {
            NAA[sp, a+1, y+1] <- NAA[sp, a, y] * exp(-ZAA[sp, a, y])
          }
          # Plus group
          NAA[sp, nages, y+1] <- NAA[sp, nages, y+1] + NAA[sp, nages, y] * exp(-ZAA[sp, nages, y])
        }
        avgNAA[sp,,y] = NAA[sp, ,y] * (1-exp(-ZAA[sp, ,y]))/ZAA[sp,,y]

        # Calculate annual metrics
        Total_Biom[sp, y] <- sum(NAA[sp, ,y] * WAA[sp, ])
        SSB[sp, y] <- sum(NAA[sp, ,y] * WAA[sp, ] * MatAA[sp, ]) * 0.5
        Catch[sp, y] <- sum(CAA[sp, ,y] * WAA[sp, ])
      }
    }

    # Available food ----
    avail_food <- array(0, dim = c(nspp, nages, nyrs))
    for(rsp in 1:nspp) {    # Predator species loop
      for(r_age in 1:nages) { # Predator age loop
        for(y in 1:nyrs) {
          for(ksp in 1:nspp) {
            for(k_age in 1:nages) { # Prey age loop
              avail_food[rsp, r_age, y] = avail_food[rsp, r_age, y] + suitability[rsp, ksp, r_age, k_age] * avgNAA[ksp, k_age, y] * WAA[ksp, k_age]
            }
          }
          # Other food (NEW logic handling normalization switch)
          if(normalize_suitability) {
            avail_food[rsp, r_age, y] = avail_food[rsp, r_age, y] + other_food[rsp] * suit_other_mat[rsp, r_age]
          } else {
            avail_food[rsp, r_age, y] = avail_food[rsp, r_age, y] + other_food[rsp] * suit_other[rsp]
          }
        }
      }
    }

    # Predation mortality ----
    M2_at_age[] <- 0
    B_eaten_as_prey <- array(0, dim = c(nspp, nages, nyrs))
    diet_prop <- array(0, dim = c(nspp, nspp, nages, nages, nyrs))
    for(ksp in 1:nspp) {
      for(k_age in 1:nages) { # Prey age loop
        for(rsp in 1:nspp) {    # Predator species loop
          for(r_age in 1:nages) { # Predator age loop
            for(y in 1:nyrs) {

              if(avail_food[rsp, r_age, y] > 0){
                # MSVPA
                # - M2
                M2_at_age[ksp, k_age, y] = M2_at_age[ksp, k_age, y] + avgNAA[rsp, r_age, y] * ration[rsp, r_age] * suitability[rsp, ksp, r_age, k_age] / avail_food[rsp, r_age, y]

                # Biomass consumed as prey
                B_eaten_as_prey[ksp, k_age, y] = B_eaten_as_prey[ksp, k_age, y] + avgNAA[ksp, k_age, y] * WAA[ksp, k_age] * avgNAA[rsp, r_age, y] * ration[rsp, r_age] * suitability[rsp, ksp, r_age, k_age] / avail_food[rsp, r_age, y]

                # Diet
                diet_prop[rsp, ksp, r_age, k_age, y] = avgNAA[ksp, k_age, y] * WAA[ksp, k_age] * suitability[rsp, ksp, r_age, k_age] / avail_food[rsp, r_age, y]
              }
            }
          }
        }
      }
    }

    # End iterations
  }

  # Observation model
  ObsCatch <- Catch * rlnorm(nyrs * nspp, 0, sigma_catch)
  SrvIdx <- srv_q * Total_Biom * rlnorm(nyrs * nspp, 0, sigma_srv)

  # Age composition data (simplified multinomial)
  ObsFishAges <- array(0, dim=c(nspp, nages, nyrs))
  ObsSrvAges <- array(0, dim=c(nspp, nages, nyrs))

  for(sp in 1:nspp){
    for(y in 1:nyrs) {
      ObsFishAges[sp, ,y] <- rmultinom(1, fish_ISS, CAA[sp, ,y])
      ObsSrvAges[sp, ,y] <- rmultinom(1, srv_ISS, NAA[sp, ,y] * srv_sel[sp,])
    }
  }

  # Age composition data (simplified multinomial)
  ObsFishCAAL <- array(0, dim=c(nspp, nages, nlengths, nyrs))
  ObsSrvCAAL <- array(0, dim=c(nspp, nages, nlengths, nyrs))


  for(sp in 1:nspp){
    for(y in 1:nyrs) {
      for(ln in 1:nlengths) {
        if(sum(growth_matrix[sp, , ln]) > 0){
          # Fishery CAAL
          ObsFishCAAL[sp,,ln, y] <- rmultinom(1, fish_CAAL_ISS, CAA[sp, ,y] * growth_matrix[sp, , ln])
          # Survey CAAL
          ObsSrvCAAL[sp,,ln, y] <- rmultinom(1, srv_CAAL_ISS, NAA[sp, ,y] * srv_sel[sp,] * growth_matrix[sp, , ln])
        }
      }
    }
  }

  # Export the array instead of the vector if normalized
  if(normalize_suitability) {
    suit_other <- suit_other_mat
  }

  # Set up Rceattle data -------------------------------------------------------------
  simData <- list()

  # * Data controls
  simData$nspp <- nspp
  simData$styr <- 1
  simData$endyr <- nyrs
  simData$projyr <- nyrs+10
  simData$spnames <- paste0("Species",1:nspp)
  simData$nsex <- rep(1, nspp)
  simData$spawn_month <- rep(0, nspp)
  simData$nages <- rep(nages, nspp)
  simData$minage <- rep(1, nspp)
  simData$nlengths <- rep(nlengths, nspp)
  simData$estDynamics <- rep(0, nspp)
  simData$pop_wt_index <- 1:nspp
  simData$ssb_wt_index <- 1:nspp
  simData$alpha_wt_len = rep(0.00001, nspp)
  simData$beta_wt_len = rep(3, nspp)
  simData$pop_age_transition_index <- rep(1, nspp)
  simData$sigma_rec = rep(1, nspp)
  simData$other_food <- rep(1e5, nspp)
  simData$estDynamics = rep(0, nspp)

  # * Fleet control
  simData$fleet_control <- data.frame(
    Fleet_name = paste0(c("Survey", "Fishery"), rep(paste(" Species", 1:nspp), each = 2)),
    Fleet_code = 1:4,
    Fleet_type = rep(2:1, nspp),
    Species = rep(1:nspp, each = 2),
    Month = 0,
    Selectivity_index = 1:(nspp * 2),
    Selectivity = "Logistic", # age or length based logistic
    Selectivity_dimension = ifelse(use_size_sel, "Length", "Age"), # age or length based logistic
    N_sel_bins = NA,
    Sel_curve_pen1 = NA,
    Sel_curve_pen2 = NA,
    Time_varying_sel = 0,
    Time_varying_sel_sd = 1,
    Bin_first_selected = 1,
    Sel_norm_bin = NA,
    Sel_norm_bin_upper = NA,
    Comp_distribution = "Multinomial",
    Comp_weights = 1,
    CAAL_distribution = 0,
    Observation_units = 1,
    Weight_index = rep(1:nspp, each = 2),
    Age_transition_index = 1,
    Catchability_index = 1:(nspp * 2),
    Catchability = rep(c(0, NA), nspp),
    Catchability_init = rep(c(1, NA), nspp),
    Catchability_prior_sd = rep(c(0.2, NA), nspp),
    Time_varying_q = rep(c(0, NA), nspp),
    Time_varying_q_sd = rep(c(1, NA), nspp),
    Estimate_index_sd = rep(c(0, NA), nspp),
    Index_sd = rep(c(1, NA), nspp),
    Estimate_catch_sd = rep(c(NA, 0), nspp),
    Catch_sd = rep(c(NA, 1), nspp),
    Proj_F_proportion = rep(c(NA, 1), nspp),
    CAAL_weights = 1,
    Comp_weights_mcallister = 1
  )

  # * Index data
  simData$index_data <- data.frame(
    Fleet_name = rep(paste("Survey Species", 1:nspp), each = nyrs),
    Fleet_code = rep(seq(1, nspp * 2, by = 2), each = nyrs),
    Species = rep(1:nspp, each = nyrs),
    Year = rep(years, nspp),
    Month = 0,
    Selectivity_block = 1,
    Observation = as.numeric(t(SrvIdx)),
    Log_sd = sigma_srv
  )

  # * Catch data
  simData$catch_data <- data.frame(
    Fleet_name = rep(paste("Fishery Species", 1:nspp), each = nyrs),
    Fleet_code = rep(seq(2, nspp * 2, by = 2), each = nyrs),
    Species = rep(1:nspp, each = nyrs),
    Year = rep(years, nspp),
    Month = 0,
    Selectivity_block = 1,
    Catch = as.numeric(t(ObsCatch)),
    Log_sd = sigma_catch
  )

  # * Comp data
  # - Index
  tmp <- do.call(rbind, lapply(seq_len(dim(ObsSrvAges)[1]), function(i) t(ObsSrvAges[i, , ])))
  colnames(tmp) <- paste0("Comp_",ages)
  index_comp <- cbind(
    data.frame(
      Fleet_name = rep(paste("Survey Species", 1:nspp), each = nyrs),
      Fleet_code = rep(seq(1, nspp * 2, by = 2), each = nyrs),
      Species = rep(1:nspp, each = nyrs),
      Sex = 0,
      Age0_Length1 = 0,
      Year = rep(years, nspp),
      Month = 0,
      Sample_size = rowSums(tmp)),
    tmp
  )

  # - Fishery
  tmp <- do.call(rbind, lapply(seq_len(dim(ObsFishAges)[1]), function(i) t(ObsFishAges[i, , ])))
  colnames(tmp) <- paste0("Comp_",ages)
  fishery_comp <- cbind(
    data.frame(
      Fleet_name = rep(paste("Fishery Species", 1:nspp), each = nyrs),
      Fleet_code = rep(seq(2, nspp * 2, by = 2), each = nyrs),
      Species = rep(1:nspp, each = nyrs),
      Sex = 0,
      Age0_Length1 = 0,
      Year = rep(years, nspp),
      Month = 0,
      Sample_size = rowSums(tmp)),
    tmp
  )

  simData$comp_data <- rbind(index_comp, fishery_comp)


  # CAAL
  # - Fishery
  fish_caal_list <- list()
  ind = 1
  for(sp in 1:nspp){
    for(y in 1:nyrs) {
      tmp <- data.frame(t(ObsFishCAAL[sp,,, y]))
      colnames(tmp) <- paste0("CAAL_", 1:nages)
      tmp <- tmp |>
        dplyr::mutate(
          Fleet_name = paste("Fishery Species", sp),
          Fleet_code = 2*(sp-1)+2,
          Species = sp,
          Sex = 0,
          Year = y,
          Length = lengths,
          Sample_size = fish_CAAL_ISS
        ) |>
        dplyr::select("Fleet_name", "Fleet_code", "Species", "Sex", "Year", "Length", "Sample_size", paste0("CAAL_", 1:nages))
      fish_caal_list[[ind]] <- tmp
      ind = ind + 1
    }
  }
  fish_caal <- do.call("rbind", fish_caal_list)

  # - Survey
  srv_caal_list <- list()
  ind = 1
  for(sp in 1:nspp){
    for(y in 1:nyrs) {
      tmp <- data.frame(t(ObsSrvCAAL[sp,,, y]))
      colnames(tmp) <- paste0("CAAL_", 1:nages)
      tmp <- tmp |>
        dplyr::mutate(
          Fleet_name = paste("Survey Species", sp),
          Fleet_code = 2*(sp-1)+1,
          Species = sp,
          Sex = 0,
          Year = y,
          Length = lengths,
          Sample_size = srv_CAAL_ISS
        ) |>
        dplyr::select("Fleet_name", "Fleet_code", "Species", "Sex", "Year", "Length", "Sample_size", paste0("CAAL_", 1:nages))
      srv_caal_list[[ind]] <- tmp
      ind = ind + 1
    }
  }
  srv_caal <- do.call("rbind", srv_caal_list)
  simData$caal_data <- rbind(fish_caal, srv_caal)

  #  Empirical selectivity
  simData$emp_sel <- data.frame(matrix(NA, nrow = 0, ncol = 5 + nages))
  colnames(simData$emp_sel ) = c("Fleet_name", "Fleet_code", "Species", "Sex", "Year", paste("Comp_", 1:nages))

  # Input N-at-age
  simData$NByageFixed <- data.frame(matrix(NA, nrow = 0, ncol = 4 + nages))
  colnames(simData$NByageFixed ) = c("Species_name ", "Species", "Sex", "Year", paste("Age", 1:nages))

  # * Age transition matrix
  age_trans_list <- lapply(1:nspp, function(sp) {
    tmp <- matrix(0, nrow = nages, ncol = nlengths)
    colnames(tmp) <- paste0("Length_",1:nlengths)
    cbind(
      data.frame(Age_transition_name = paste("Species", sp),
                 Age_transition_index = sp,
                 Species = sp,
                 Sex = 0,
                 Age = 1:nages),
      tmp
    )
  })
  simData$age_trans_matrix <- do.call(rbind, age_trans_list)

  # * Age error
  age_error_list <- lapply(1:nspp, function(sp) {
    tmp <- as.data.frame(diag(1,nages))
    colnames(tmp) <- paste0("Obs_age",1:nages)
    cbind(data.frame(Species = sp,
                     True_age = 1:nages),
          tmp
    )
  })
  simData$age_error <- do.call(rbind, age_error_list)

  # * Weight-at-age
  weight_list <- lapply(1:nspp, function(sp) {
    weight_tmp <- as.data.frame(matrix(WAA[sp,], ncol = nages))
    colnames(weight_tmp) <- paste0("Age",1:nages)
    cbind(data.frame(Wt_name = paste("Species", sp),
                     Wt_index = sp,
                     Species = sp,
                     Sex = 0,
                     Year = 0),
          weight_tmp
    )
  })
  simData$weight <- do.call(rbind, weight_list)

  # * Maturity
  maturity_list <- lapply(1:nspp, function(sp) {
    mat_tmp <- as.data.frame(matrix(MatAA[sp,], ncol = nages))
    colnames(mat_tmp) <- paste0("Age",1:nages)
    cbind(data.frame(Species = sp),
          mat_tmp
    )
  })
  simData$maturity <- do.call(rbind, maturity_list)

  # * Sex ratio
  sexratio <- as.data.frame(matrix(0.5, nrow = nspp, ncol = nages))
  colnames(sexratio) <- paste0("Age",1:nages)
  simData$sex_ratio <- cbind(data.frame(Species = 1:nspp),
                             sexratio
  )

  # * Mortality
  mort_list <- lapply(1:nspp, function(sp) {
    mort_tmp <- as.data.frame(matrix(M[sp], nrow = 1, ncol = nages))
    colnames(mort_tmp) <- paste0("Age",1:nages)
    cbind(data.frame(Species = sp,
                     Sex = 0),
          mort_tmp
    )
  })
  simData$M1_base <- do.call(rbind, mort_list)

  # * Environmental data ----
  simData$env_data <- data.frame(Year = 1:nyrs, EnvData = 1)

  # * Base data ----
  simDataBase <- simData

  # * Ration ----
  ration_list <- lapply(1:nspp, function(sp) {
    ration_tmp <- as.data.frame(matrix(WAA[sp,] * 50, ncol = nages))
    colnames(ration_tmp) <- paste0("Age",1:nages)
    cbind(data.frame(
      Species = sp,
      Sex = 0,
      Year = 0),
      ration_tmp
    )
  })
  simData$ration_data <- do.call(rbind, ration_list)

  # * Bioenergetics
  simData$Ceq <- rep(4,nspp)
  simData$Cindex <- rep(1, nspp)
  simData$Pvalue <- rep(1, nspp)
  simData$fday <- rep(1, nspp)
  simData$CA <- rep(1, nspp)
  simData$CB <- rep(-1, nspp)
  simData$Qc <- rep(1,nspp)
  simData$Tco <- rep(1, nspp)
  simData$Tcm <- rep(1, nspp)
  simData$Tcl <- rep(1, nspp)
  simData$CK1 <- rep(1, nspp)
  simData$CK4 <- rep(1, nspp)
  simData$Diet_distribution <- rep(0,nspp)
  simData$Diet_comp_weights <- rep(1,nspp)

  # * Diet proportion ----
  # Get all combinations of indices
  idx <- which(!is.na(diet_prop), arr.ind = TRUE)  # or diet_prop != 0 for nonzero only

  # Build the data frame directly
  simData$diet_data <- data.frame(
    Year = idx[, 5],
    Pred = idx[, 1],
    Prey = idx[, 2],
    Pred_sex = 0,
    Prey_sex = 0,
    Pred_age = idx[, 3],
    Prey_age = idx[, 4],
    Sample_size = 1000,
    Stomach_proportion_by_weight = diet_prop[idx]
  )
  simData$diet_data <- simData$diet_data |>
    dplyr::filter(!is.na(Pred))

  # Return list of true and observed values ----
  return(list(
    model_quantities = list(
      NAA = NAA,
      avgNAA = avgNAA,
      CAA = CAA,
      FAA = FAA,
      MAA = MAA,
      SSB = SSB,
      Total_Biom = Total_Biom,
      B_eaten_as_prey = B_eaten_as_prey,
      Catch = Catch,
      ObsCatch = ObsCatch,
      SrvIdx = SrvIdx,
      ObsFishAges = ObsFishAges,
      ObsSrvAges = ObsSrvAges,
      fish_sel = fish_sel,
      srv_sel = srv_sel,
      WAA = WAA,
      growth_matrix = growth_matrix,
      lengths = lengths,
      length_at_age = length_at_age,
      fish_size_sel = fish_size_sel,
      srv_size_sel = srv_size_sel,
      MatAA = MatAA,
      M = M,
      Mtv = Mtv,
      M_vec = M_vec,
      srv_q = srv_q,
      rec_devs = rec_devs,
      init_devs = init_devs,
      M2_at_age = M2_at_age,
      vulnerability = vulnerability,
      suitability = suitability,
      suit_other = suit_other, # Automatically outputs as matrix if normalized
      avail_food = avail_food,
      diet_prop = diet_prop,
      ration = ration
    ),
    data_list = simData,
    base_data = simDataBase
  ))
}
