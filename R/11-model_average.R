#' Model average of derived quantities
#'
#' @param object list of Rceattle model objects
#' @param Rceattle deprecated name for `object`, still accepted so existing
#'   scripts keep working. Supplying both is an error.
#' @param weights vector of weights to be used for weighting models
#' @param uncertainty TRUE/FALSE Sample uncertainty across derived quantities using weighted bootstrap from the asymptotic  distribution of MLEs
#' @param nboot Number of bootstraps taken from asymptotic distribution of MLEs. Default = 10000
#'
#' @return an Rceattle object with derived quantities weighted by the specified weights. The length of the derived quantities spans the years which overlap across all models.
#' @examples
#' \dontrun{
#' # Equal weights across two candidate models.
#' model_average(list(fit1, fit2), weights = c(1, 1))
#' # Weighted by AIC, with bootstrapped uncertainty.
#' model_average(list(fit1, fit2), weights = c(0.7, 0.3), uncertainty = TRUE)
#' }
#' @export
#'
model_average <- function(object = NULL, weights = NULL, uncertainty = FALSE, nboot = 10000, Rceattle = NULL){
  # `Rceattle` was the old name for `object`; see R/0-deprecate.R.
  if (!missing(Rceattle))
    object <- .rce_deprecated_arg(Rceattle, !missing(object), "Rceattle", "object", "model_average")

  # A list of fits is also valid here, so check for either shape.
  if (!inherits(object, "Rceattle") && !is.list(object)) {
    stop("`object` must be a fitted Rceattle model, or a list of them.", call. = FALSE)
  }

  # --------------------------------------------------------------------------------------------
  # Average derived quantities of models
  # --------------------------------------------------------------------------------------------
  if(inherits(object, "Rceattle")){
    stop("Only one model provided")
  }
  if(is.null(weights)){
    weights = rep(1/length(object), length(object))
  }
  weights <- weights/sum(weights) # Normalize

  # Extract number of species of each model
  nspp <- sapply(object, function(x) x$data_list$nspp)
  nspp <- unique(nspp)
  if(length(nspp) > 1){stop("Number of species does not match across models")}

  # Extract timespan of each model
  years <- sapply(object, function(x) x$data_list$styr:x$data_list$projyr)
  endyrs <- sapply(object, function(x) x$data_list$endyr)
  styrs <- sapply(object, function(x) x$data_list$styr)
  projyrs <- sapply(object, function(x) x$data_list$projyr)
  nyrs <- sapply(years, length)

  min_projyr <- min(projyrs, na.rm = TRUE) # Find coverage across all models
  min_hindyr <- min(endyrs, na.rm = TRUE) # Find coverage across all models
  max_styr <- max(styrs, na.rm = TRUE)

  # Initialize model average object
  mod_avg <- object[[which(nyrs == min(nyrs))[1]]] # Copy a model of the smallest length to fill in
  mod_avg_rel_proj_yrs <- (max_styr - mod_avg$data_list$styr + 1) : (min_projyr - mod_avg$data_list$styr + 1) # Find years of overlap
  mod_avg_rel_hind_yrs <- (max_styr - mod_avg$data_list$styr + 1) : (min_hindyr - mod_avg$data_list$styr + 1) # Find years of overlap
  mod_avg$estimated_params <- mod_avg$initial_params <- mod_avg$opt <- mod_avg$run_time <- mod_avg$obj <- mod_avg$map <- mod_avg$bounds <- NULL
  if(!uncertainty){mod_avg$sdrep <- NULL}

  # -- Set all quantities to zero
  mod_avg$quantities <- sapply(mod_avg$quantities, function(x) replace(x, values = rep(0, length(x))))



  # Average out derived quantities (Note: dimensions may be different)
  for(i in 1:length(object)){
    sub_rel_proj_yrs <- (max_styr - object[[i]]$data_list$styr + 1) : (min_projyr - object[[i]]$data_list$styr + 1) # Find years of overlap
    sub_rel_hind_yrs <- (max_styr - object[[i]]$data_list$styr + 1) : (min_hindyr - object[[i]]$data_list$styr + 1) # Find years of overlap

    # Average quantities of interest
    # Biomass components
    mod_avg$quantities$R[, mod_avg_rel_proj_yrs] <- mod_avg$quantities$R[, mod_avg_rel_proj_yrs] + object[[i]]$quantities$R[, sub_rel_proj_yrs] * weights[i]
    mod_avg$quantities$avg_R <- mod_avg$quantities$avg_R + object[[i]]$quantities$avg_R * weights[i]

    mod_avg$quantities$N_at_age[,,,mod_avg_rel_proj_yrs] <- mod_avg$quantities$N_at_age[,,,mod_avg_rel_proj_yrs] + object[[i]]$quantities$N_at_age[,,,sub_rel_proj_yrs] * weights[i]
    mod_avg$quantities$avgN_at_age[,,,mod_avg_rel_proj_yrs] <- mod_avg$quantities$avgN_at_age[,,,mod_avg_rel_proj_yrs] + object[[i]]$quantities$avgN_at_age[,,,sub_rel_proj_yrs] * weights[i]
    #mod_avg$quantities$biomass_at_age[,,mod_avg_rel_proj_yrs] <- mod_avg$quantities$biomass_at_age[,,mod_avg_rel_proj_yrs] + object[[i]]$quantities$biomass_at_age[,,sub_rel_proj_yrs] * weights[i]
    mod_avg$quantities$ssb_at_age[,,mod_avg_rel_proj_yrs] <- mod_avg$quantities$ssb_at_age[,,mod_avg_rel_proj_yrs] + object[[i]]$quantities$ssb_at_age[,,sub_rel_proj_yrs] * weights[i]
    mod_avg$quantities$biomass[, mod_avg_rel_proj_yrs] <- mod_avg$quantities$biomass[, mod_avg_rel_proj_yrs] + object[[i]]$quantities$biomass[, sub_rel_proj_yrs] * weights[i]
    mod_avg$quantities$ssb[, mod_avg_rel_proj_yrs] <- mod_avg$quantities$ssb[, mod_avg_rel_proj_yrs] + object[[i]]$quantities$ssb[, sub_rel_proj_yrs] * weights[i]
    mod_avg$quantities$B_eaten_as_prey[,,,mod_avg_rel_proj_yrs] <- mod_avg$quantities$B_eaten_as_prey[,,,mod_avg_rel_proj_yrs] + object[[i]]$quantities$B_eaten_as_prey[,,,sub_rel_proj_yrs] * weights[i]
    mod_avg$quantities$sel_at_age[,,,mod_avg_rel_hind_yrs] <- mod_avg$quantities$sel_at_age[,,,mod_avg_rel_hind_yrs] + object[[i]]$quantities$sel_at_age[,,,sub_rel_hind_yrs] * weights[i]

    # Reference points
    mod_avg$quantities$SB0 <- mod_avg$quantities$SB0 + object[[i]]$quantities$SB0 * weights[i]

    # Mortality components
    mod_avg$quantities$M1_at_age <- mod_avg$quantities$M1_at_age + object[[i]]$quantities$M1_at_age * weights[i]
    mod_avg$quantities$M2_at_age[,,,mod_avg_rel_proj_yrs] <- mod_avg$quantities$M2_at_age[,,,mod_avg_rel_proj_yrs] + object[[i]]$quantities$M2_at_age[,,,sub_rel_proj_yrs] * weights[i]
    mod_avg$quantities$M_at_age[,,,mod_avg_rel_proj_yrs] <- mod_avg$quantities$M_at_age[,,,mod_avg_rel_proj_yrs] + object[[i]]$quantities$M_at_age[,,,sub_rel_proj_yrs] * weights[i]
    mod_avg$quantities$Z_at_age[,,,mod_avg_rel_proj_yrs] <- mod_avg$quantities$Z_at_age[,,,mod_avg_rel_proj_yrs] + object[[i]]$quantities$Z_at_age[,,,sub_rel_proj_yrs] * weights[i]
    mod_avg$quantities$S_at_age[,,,mod_avg_rel_proj_yrs] <- mod_avg$quantities$S_at_age[,,,mod_avg_rel_proj_yrs] + object[[i]]$quantities$S_at_age[,,,sub_rel_proj_yrs] * weights[i]
    mod_avg$quantities$M2_prop[,,,,mod_avg_rel_proj_yrs] <- mod_avg$quantities$M2_prop[,,,,mod_avg_rel_proj_yrs] + object[[i]]$quantities$M2_prop[,,,,sub_rel_proj_yrs] * weights[i]
    mod_avg$quantities$B_eaten_prop[,,,,mod_avg_rel_proj_yrs] <- mod_avg$quantities$B_eaten_prop[,,,,mod_avg_rel_proj_yrs] + object[[i]]$quantities$B_eaten_prop[,,,,sub_rel_proj_yrs] * weights[i]
  }

  # --------------------------------------------------------------------------------------------
  # Sample uncertainty across derived quantities via a weighted bootstrap from the asymptotic intervals
  # --------------------------------------------------------------------------------------------

  if(uncertainty){

    # Assuming asymptotic multivariate normal
    # - Must have estimated joint precision matrix
    # TRUE where the model's sdrep carries a random block; a getsd = FALSE fit
    # reads FALSE.
    length_ran <- sapply(object, function(x) !is.null(x$sdrep$par.random))
    joint_est <- sapply(object, function(x) is.null(x$sdrep$jointPrecision))
    if(sum(length_ran) > 0){
      if(sum(joint_est) > 0){
        stop(paste("Model(s):", which(joint_est > 0)," do not have joint precision matrices: re-estimate with getJointPrecision = TRUE"))
      }

      if(sum(length_ran) != length(object)){
        stop("Does not currently support purely fixed and random effects models together")
      }
    }



    # - Pointers
    draws <- round(weights * nboot)
    nrows <- 1:sum(draws)
    rowid <- rep(1:length(draws), draws)

    # List to save samples of SSB, B, and R
    recruitment <-
      array(NA, dim = c(nspp, length(mod_avg_rel_proj_yrs),  sum(draws)))
    ssb <-
      array(NA, dim = c(nspp, length(mod_avg_rel_proj_yrs),  sum(draws)))
    biomass <-
      array(NA, dim = c(nspp, length(mod_avg_rel_proj_yrs),  sum(draws)))


    # Loop across models
    for(i in 1:length(object)){

      # - Sample parameters from asymptotic normal distribution, whose covariance
      #   is the inverse joint precision with random effects, cov.fixed without.
      # - Models with random effects
      if(sum(length_ran) > 0){
        mle <- object[[i]]$obj$env$last.par.best # Includes fixed and random
        vcov <- solve(object[[i]]$sdrep$jointPrecision) # names(mle) == rownames(vcov)
      }


      # - Fixed effects models
      if(sum(length_ran) == 0){
        mle <- object[[i]]$obj$env$last.par.best # Includes fixed
        vcov <- object[[i]]$sdrep$cov.fixed # names(mle) == rownames(vcov)
      }
      samples <- MASS::mvrnorm(draws[i], mu = mle, Sigma = vcov)

      # - Get quantities
      quantities <- apply(samples, 1, function(x) object[[i]]$obj$report(x)[c("biomass", "ssb", "R")]) # Only want uncertainty in SSB, B, and R

      # - Subset years of interest and assign to objects
      sub_rel_proj_yrs <- (max_styr - object[[i]]$data_list$styr + 1) : (min_projyr - object[[i]]$data_list$styr + 1) # Find years of overlap

      # Extract R, B, and SSB and assign to objects
      recruitment[1:nspp, 1:length(mod_avg_rel_proj_yrs), nrows[which(rowid == i)]] <- unlist(lapply(quantities, function(x) x$R[, sub_rel_proj_yrs]))
      ssb[1:nspp, 1:length(mod_avg_rel_proj_yrs), nrows[which(rowid == i)]] <- unlist(lapply(quantities, function(x) x$ssb[, sub_rel_proj_yrs]))
      biomass[1:nspp, 1:length(mod_avg_rel_proj_yrs), nrows[which(rowid == i)]] <- unlist(lapply(quantities, function(x) x$biomass[, sub_rel_proj_yrs]))
    }

    # - Calculate SD
    # Both scales get the across-model sampling SD. The plotters prefer the
    # log-scale rows, so leaving log_* at its single-model value would draw
    # model-averaged trajectories with single-model intervals. A non-positive
    # sample gives NaN, which the plotters treat as "no log SD".
    sd_samples <- list(R           = recruitment,
                       biomass      = biomass,
                       ssb          = ssb,
                       log_R        = log(recruitment),
                       log_biomass  = log(biomass),
                       log_ssb      = log(ssb))
    for (nm in names(sd_samples)) {
      rows <- which(names(mod_avg$sdrep$value) == nm)
      if (length(rows)) {
        mod_avg$sdrep$sd[rows] <- sqrt(apply(sd_samples[[nm]], c(1, 2), var))
      }
    }

    # - Save samples
    mod_avg$asymptotic_samples <- list(recruitment = recruitment, biomass = biomass, ssb = ssb)
  }

  return(mod_avg)
}
