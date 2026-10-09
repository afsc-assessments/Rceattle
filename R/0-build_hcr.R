#' @title Specify the harvest control rule (HCR) used for Rceattle
#'
#' @description Defines the harvest control rule and its associated reference
#'   points. The rules are tabulated, with which ones work under predation, in
#'   \code{vignette("projections-and-reference-points")}; the formulations of
#'   the three ramped rules are below.
#'
#' @param HCR Harvest control rule, as an integer or string alias:
#'   \code{0}/\code{"NoFishing"} (no catch, hindcast only),
#'   \code{1}/\code{"CMSY"} (maximize joint catch, optionally holding depletion
#'   at or above \code{Plimit}), \code{2}/\code{"ConstantF"} (constant input F
#'   at \code{Ftarget}), \code{3}/\code{"ConstantFSSB"} (the F reaching
#'   \code{Ftarget} percent of SSB0 by the end of the projection),
#'   \code{4}/\code{"ConstantFSPR"} (constant Fspr, scalable by \code{Fmult}
#'   following NEFSC), \code{5}/\code{"NPFMC"}, \code{6}/\code{"PFMC"} or
#'   \code{7}/\code{"SESSF"}. Default \code{0}. **Only 0, 1, 2, 3 and 6 work in
#'   multispecies mode.**
#' @param DynamicHCR TRUE/FALSE. Use dynamic rather than static reference points
#'   (default FALSE).
#' @param Ftarget Target fishing mortality rate (yr^-1), SPR- or
#'   depletion-based, or the input F for projections. For SPR F40%, enter 0.40.
#' @param Flimit Limit fishing mortality rate (yr^-1), SPR- or depletion-based.
#'   For SPR F35%, enter 0.35.
#' @param Ptarget Target spawning-stock biomass as a percentage of static or
#'   dynamic spawning-stock biomass at F = 0, which accounts for recruitment.
#' @param Plimit Limit spawning-stock biomass on the same basis.
#' @param Alpha Parameter used in the NPFMC Tier 3 rule.
#' @param Pstar Quantile for the uncertainty buffer, given \code{Flimit} and
#'   \code{Sigma}.
#' @param Sigma Standard deviation of the normal uncertainty buffer, given
#'   \code{Flimit} and \code{Pstar}.
#' @param Fmult Multiplier on the target F (default 1), scaling it following
#'   NEFSC convention.
#' @param HCRorder For multispecies models, the order in which to project
#'   fishing (predators first, then prey, for instance).
#'
#' @details
#' \code{hcr = 5} or \code{"NPFMC"}, the North Pacific Fishery Management
#' Council Tier 3 spawner-per-recruit rule:
#' 	Stock status: \eqn{SB > SB at Ftarget}
#' 	\deqn{Fofl = Flimit}
#' 	\deqn{Fuse = Ftarget}
#' 	Stock status: \eqn{Alpha < SB / SB at Ftarget \le 1}
#' 	\deqn{Fofl = Flimit * (SB / SB_{Ftarget} - Alpha) / (1 - Alpha)}
#' 	\deqn{Fuse = Ftarget * (SB / SB_{Ftarget} - Alpha) / (1 - Alpha)}
#' 	Stock status: \eqn{SB / SB at Ftarget \le Alpha} or \eqn{SB < Plimit * SB0}
#' 	\deqn{Fofl = 0}
#' 	\deqn{Fuse = 0}
#'
#' \code{hcr = 6} or \code{"PFMC"}, the Pacific Fishery Management Council
#' category 1 40-10 annual catch limit rule (PFMC 2020), using Fspr in single
#' species or the F reaching X% of SSB0 under predation. The uncertainty buffer
#' is \code{qnorm(Pstar, Flimit, Sigma)}, which equals
#' \code{Flimit + qnorm(Pstar, 0, Sigma)} and so sits *below* \code{Flimit}
#' when \eqn{Pstar < 0.5}. The taper runs between \eqn{SB0 * Plimit} and
#' \eqn{SB0 * Ptarget}, so the 40-10 shape needs \code{Ptarget = 0.40} and
#' \code{Plimit = 0.10}; the default \code{Plimit = 0} gives a 40-0 rule.
#' 	Stock status: \eqn{SB > SB0 * Ptarget}
#' 	\deqn{Fofl = Flimit}
#' 	\deqn{Fuse = qnorm(Pstar, Flimit, Sigma)}
#' 	Stock status: \eqn{SB0 * Plimit < SB \le SB0 * Ptarget}
#' 	\deqn{Fofl = Flimit}
#' 	\deqn{Fuse = qnorm(Pstar, Flimit, Sigma) * \frac{SB0 * Ptarget * (SB - SB0 * Plimit)}{SB * SB0 * (Ptarget - Plimit)}}
#' 	Stock status: \eqn{SB < SB0 * Plimit}
#' 	\deqn{Fofl = 0}
#' 	\deqn{Fuse = 0}
#'
#' \code{hcr = 7} or \code{"SESSF"}, the Southern and Eastern Scalefish and
#' Shark Fishery Tier 1 spawner-per-recruit rule, with F_Limit = F_(20%) and
#' B_Limit = SB_20 (AFMA 2017):
#' 	Stock status: \eqn{SB > SB0 * Ptarget}
#' 	\deqn{Fofl = Flimit}
#' 	\deqn{Fuse = Ftarget}
#' 	Stock status: \eqn{SB0 * Ptarget > SB > SB0 * Plimit}
#' 	\deqn{Fofl = Flimit * (SB / (SB0 * Plimit) - 1)}
#' 	\deqn{Fuse = Ftarget * (SB / (SB0 * Plimit) - 1)}
#' 	Stock status: \eqn{SB < SB0 * Plimit}
#' 	\deqn{Fofl = 0}
#' 	\deqn{Fuse = 0}
#'
#' @return A \code{list} containing the harvest control rule and associated
#'   biological reference points.
#' @export
#' @examples
#' # Tier 3 NPFMC control rule: F40% as FABC, F35% as FOFL.
#' build_hcr(HCR = "NPFMC", DynamicHCR = FALSE, Ftarget = 0.4, Flimit = 0.35)
#'
#' # No fishing: the projection a B0 reference point is read from.
#' build_hcr(HCR = "NoFishing")
##'
build_hcr <- function(HCR = 0, DynamicHCR = FALSE, Ftarget = 0.40, Flimit = 0.35, Ptarget = 0.4, Plimit = 0.0, Alpha = 0.05, Pstar = 0.45, Sigma = 0.5, Fmult = 1, HCRorder = 1) {
  if(0 %in% Alpha & HCR %in% c(5, "NPFMC")){stop(paste0("Alpha = 0 for NPFMC Tier 3 HCR, divide by zero error"))}
  list(HCR = HCR, DynamicHCR = DynamicHCR, Ftarget = Ftarget, Flimit = Flimit, Ptarget = Ptarget, Plimit = Plimit, Alpha = Alpha, Pstar = Pstar, Sigma = Sigma, Fmult = Fmult, HCRorder = HCRorder)
}


#' Function to construct the TMB map argument for CEATTLE for projecting under alternative harvest control rules
#'
#' @description Reads a data list and map to update the map argument based on the HCR specified in \code{\link{build_hcr}}
#'
#' @param data_list an Rceattle data_list
#' @param map a map object created from \code{\link{build_map}}.
#' @param debug logical. If TRUE, turns off all parameters for debugging (default = FALSE).
#' @param all_params_on logical. If TRUE, leaves all hindcast parameters turned on (default = FALSE).
#' @param HCRiter for multi-species models, the order in which to project fishing (e.g. predators first, then prey)
#'
#' @return a list of map arguments for each parameter
#' @export
build_hcr_map <- function(data_list, map, debug = FALSE, all_params_on = FALSE, HCRiter = 1){

  # Turn off all population/fleet parameters ---- and turn on Fspr parameters
  params_on <- 1:data_list$nspp
  if(!all_params_on){
    map$mapList = sapply(map$mapList, function(x) replace(x, values = rep(NA, length(x))))
    yrs_proj = data_list$endyr:data_list$projyr - data_list$styr
    params_on <- c(1:data_list$nspp)[which(data_list$HCRorder <= HCRiter)]
  }

  # Turn on Fspr parameters depending on HCR ----
  # -- HCR = 0: No catch - Params off
  # -- HCR = 1: Constant catch - Params off
  # -- HCR = 2: Constant input F - Params off
  # -- HCR = 3: F that achieves X% of SSB0 in the end of the projection - Ftarget on
  # -- HCR = 4: Constant target Fspr - Ftarget on
  # -- HCR = 5: NPFMC Tier 3 - Flimit and Ftarget on
  # -- HCR = 6: PFMC Cat 1 - Flimit on
  # -- HCR = 7: SESSF Tier 1 - Flimit and Ftarget on
  # --- Dynamic BRPS - 1 value per species and year
  if(!debug){
    if(data_list$HCR == "CMSY"){ # CMSY
      map$mapList$log_Ftarget[params_on] <- params_on
    }

    if(data_list$HCR == "ConstantF"){ # Fixed F - still have Flimit for single-species
      map$mapList$log_Flimit[params_on] <- params_on
    }
    if(data_list$HCR == "ConstantFSSB"){
      map$mapList$log_Ftarget[params_on] <- params_on
    }
    if(data_list$HCR == "ConstantFSPR"){
      map$mapList$log_Ftarget[params_on] <- params_on
      map$mapList$log_Flimit[params_on] <- params_on
    }
    if(data_list$HCR %in% c("NPFMC", "SESSF")){
      map$mapList$log_Ftarget[params_on] <- params_on
      map$mapList$log_Flimit[params_on] <- params_on
    }
    if(data_list$HCR == "PFMC"){
      map$mapList$log_Flimit[params_on] <- params_on
    }


    # Turn off SPR parameters for special cases ----
    # -- Turn off SPR parameters for species with no fishing (sum(Proj_F_proportion) == 0)
    # -- Turn off SPR parameters for species with fixed Nbyage
    for(sp in 1:data_list$nspp){

      # Check proj F if proj F prop is all 0
      prop_check <- data_list$fleet_control$Proj_F_proportion[which(data_list$fleet_control$Species == sp & data_list$fleet_control$Fleet_type == "Fishery")]
      # Turn off future F only when *every* fishery for this species is inactive.
      if(sum(prop_check, na.rm = TRUE) == 0){
        message(paste("F_prop for species",sp,"sums to 0"))
        map$mapList$log_Ftarget[sp] <- NA
        map$mapList$log_Flimit[sp] <- NA
      }

      # Fixed n-at-age: Turn off parameters
      if(data_list$estDynamics[sp] > 0){
        map$mapList$log_Ftarget[sp] <- NA
        map$mapList$log_Flimit[sp] <- NA
      }
    }
  }
  if(debug){
    map$mapList$dummy = 1
  }

  # Convert to factor ----
  map$mapFactor <- sapply(map$mapList, factor)
  return(map)
}
