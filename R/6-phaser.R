#' Run TMB with ADMB-style phased estimation
#'
#' Fits a TMB model in phases, like ADMB: parameters are switched on in stages
#' rather than all at once, which stabilizes a difficult optimization. The
#' optimizer (nlminb by default) runs once per phase; at each phase the map holds
#' fixed any parameter whose phase is greater than the current phase, so
#' parameters turn on progressively as the phase counter advances. `phases` is a
#' tagged list giving each named parameter its integer phase; parameters not
#' listed are estimated from the first phase onward.
#'
#' @param  data A list to be passed to TMB
#' @param  parameters A list of parameters of the model
#' @param  map a list of map object from the model
#' @param  random A character vector of names of parameters that are random effects
#' @param  phases Tagged list assigning each named parameter its integer phase (as returned by [set_phases()]).
#' @param silent logical. If TRUE, suppresses output from TMB (default = TRUE).
#' @param use_gradient logical. If TRUE, uses gradient in optimization (default = TRUE).
#' @param control A list of control parameters. For details see \code{?nlminb}
#' @param  model_name A string describing the model name. Must be the name of your .cpp file
#' @return The parameter list estimated in the final phase, with a per-phase convergence log attached as the `phase_log` attribute.
#' @author Gavin Fay https://github.com/kaskr/TMB_contrib_R/blob/master/TMBphase/R/TMBphase.R
#' @export
#'
TMBphase <- function(data, parameters, map, random, phases, model_name,
                     silent, use_gradient = TRUE,
                     control = list(eval.max = 1e+09, iter.max = 1e+09, trace = 0)) {

  # function to fill list component with a factor
  fill_vals <- function(x,vals){rep(as.factor(vals), length(x))}

  #loop over phases
  phase_log <- list()
  for (phase_cur in 1:max(unlist(phases))) {
    #phase_cur <- 1

    # work out the map for this phase
    # if phases for parameters is less than the current phase
    # then map will contain a factor filled with NAs
    map_use <- map
    # Iterate by name: `phases` covers only parameters with a declared phase,
    # and its order is the phase table's, not the parameter list's.
    for (nm in names(phases)) {
      if (phases[[nm]] > phase_cur) {
        map_val <- which(names(map_use) %in% nm)
        if (length(map_val) == 1L) {
          map_use[[map_val]] <- fill_vals(map[[map_val]], NA)
        }
      }
    }

    # Phase random effects as PENALISED fixed effects. During phasing, do NOT pass
    # `random` -- the deviates (rec_dev, sel_coff_dev, index_q_dev, ...) are estimated
    # as ordinary fixed effects, constrained by their density/penalty -- and hold the
    # RE variance / correlation hyperparameters at their initial values. This builds up
    # realistic non-zero deviates before the Laplace approximation is switched on in the
    # final hindcast fit (which uses the full map + `random`, seeded from here). A flat
    # random field integrated from zero gives a NaN marginal objective, so it must be
    # warmed up first. Fits with no random effects keep the original behaviour.
    random_use <- random
    if (length(random) > 0) {
      random_use <- NULL
      re_var_pars <- c("R_log_sd", "sel_dev_log_sd", "sel_curve_pen",
                       "index_q_dev_log_sd", "index_q_log_sd",
                       "M1_dev_log_sd", "M1_rho",
                       "log_sigma_linkage", "trans_rho_linkage", "log_obs_sd_linkage",
                       "growth_log_sd")
      for (vp in intersect(re_var_pars, names(map_use)))
        map_use[[vp]] <- fill_vals(map_use[[vp]], NA)
    }

    # initialize the parameters at values in previous phase
    params_use <- parameters
    if (phase_cur>1) params_use <- obj$env$parList(opt$par)

    # Fit the model
    obj <- TMB::MakeADFun(data,parameters =  params_use,random=random_use,DLL=model_name,map=map_use, silent = silent)

    if(use_gradient){
      opt <- nlminb(obj$par,obj$fn,obj$gr, control = control)
    }else{
      opt <- nlminb(obj$par,obj$fn)
    }
    last_par = obj$env$parList(par = obj$env$last.par.best)

    # Per-phase convergence record (for the phasing diagnostic).
    phase_log[[phase_cur]] <- list(
      phase       = phase_cur,
      n_par       = length(opt$par),
      objective   = opt$objective,
      max_grad    = tryCatch(max(abs(obj$gr(opt$par))),
                             error = function(e) NA_real_),
      convergence = opt$convergence
    )
    #close phase loop
  }

  attr(last_par, "phase_log") <- phase_log
  return(last_par)
}

#' Default phasing order for CEATTLE parameters
#'
#' @returns Tagged list mapping each parameter name to its estimation phase.
#' @export
#'
set_phases <- function(){
  phaseList = list(
    dummy = 1,
    log_pop_scalar = 4, # Scalar for input numbers-at-age
    rec_pars = 1,      # Stock-recruit parameters or log(mean rec) if no stock-recruit relationship
    R_log_sd = 2,       # Variance for annual recruitment deviats
    rec_dev = 2,       # Annual recruitment deviats
    init_dev = 2,      # Age specific initial age-structure deviates or parameters
    # sex_ratio_log_sd = 3, # Variance of sex ratio (usually fixed)
    log_M1 = 4,         #  Estimated natural or residual mortality
    log_M1_dev = 5,     #  Estimated natural or residual mortality
    M1_beta = 5,       #  Estimated natural or residual mortality
    M1_rho = 5,        #  Estimated natural or residual mortality
    M1_dev_log_sd = 5,  #  Estimated natural or residual mortality
    log_Flimit = 3,     # Estimated F limit
    log_Ftarget = 3,    # Estimated F target
    log_Finit = 3,      # Estimated fishing mortality for non-equilibrium initial age-structure
    proj_F_prop = 1,   # Fixed fleet-specific proportion of Flimit and Ftarget apportioned within each species
    log_F = 1,          # Annual fleet specific fishing mortality
    index_log_q = 3,    # Survey catchability
    index_q_dev = 5,   # Annual survey catchability deviates (if time-varying)
    index_q_log_sd = 4, # Prior SD for survey catchability deviates
    index_q_beta = 4,  # Regression coefficients for environmental linkage
    index_q_dev_log_sd = 4, # SD for annual survey catchability deviates (if time-varying)
    sel_coff = 3,      # Non-parametric selectivity coefficients
    sel_coff_dev = 4,  # Annual deviates for non-parametric selectivity coefficients
    log_sel_slp = 3,    # Slope parameters for logistic forms of selectivity
    sel_inf = 3,       # Asymptote parameters for logistic forms of selectivity
    log_sel_apical = 3, # Per-sex apical height (log multiplier), via a selectivity linkage
    log_sel_slp_dev = 5,# Annual deviates for slope parameters for logistic forms of selectivity (if time-varying)
    sel_inf_dev = 5,   # Annual deviates for asymptote parameters for logistic forms of selectivity (if time-varying)
    sel_dev_log_sd = 4, # SD for annual selectivity deviates (if time-varying)
    sel_curve_pen = 4, # Penalty for non-parametric selectivity
    index_log_sd = 2,   # Log SD for survey lognormal index likelihood (usually input)
    catch_log_sd = 2,   # Log SD for lognormal catch likelihood (usually input)
    comp_weights = 5,  # Weights for comp likelihood
    caal_weights = 5,  # Weights for CAAL likelihood
    # ,logH_1 = 6,     # Functional form parameter (not used in MSVPA functional form)
    # logH_1a = 6,     # Functional form parameter (not used in MSVPA functional form)
    # logH_1b = 6,     # Functional form parameter (not used in MSVPA functional form)
    # logH_2 = 6,      # Functional form parameter (not used in MSVPA functional form)
    # logH_3 = 6,      # Functional form parameter (not used in MSVPA functional form)
    # H_4 = 6,         # Functional form parameter (not used in MSVPA functional form)
    log_gam_a = 5,     # Suitability parameter (size-preference mean)
    log_gam_b = 5,     # Suitability parameter (size-preference sd)
    log_phi = 5,       # Suitability parameter (pred-prey vulnerability)
    log_growth_pars = 4,# Mean growth parameters
    growth_log_sd = 4,   # SD in weight-at-age at youngest and oldest ages
    weight_length_pars = 5, # Weight-length parameters
    beta_linkage = 4 # Long-format linkage-table coefficients (env covariates on process params)

  )

  return(phaseList)
}
