#' Build parameter bounds
#'
#' Function to build parameter bounds based on Holsman et al 2015 and Kinzey and Punt 2010
#'
#' @param param_list Parameter list object built from \code{\link{build_params}}
#' @param data_list a Rceattle data object
#'
#' @return List of upper and lower bounds
#' @export
#'
build_bounds <- function(param_list = NULL, data_list) {

  upper_bnd <- param_list
  lower_bnd <- param_list

  # General bounds
  for (i in 1:length(param_list)) {
    upper_bnd[[i]] <- replace(upper_bnd[[i]], values = rep(Inf, length(upper_bnd[[i]])))
    lower_bnd[[i]] <- replace(lower_bnd[[i]], values = rep(-Inf, length(lower_bnd[[i]])))
  }

  # Predator selectivity bounds for gamma suitability
  lower_bnd$log_gam_b <- replace(lower_bnd$log_gam_b, values = rep(-10, length(lower_bnd$log_gam_b)))
  upper_bnd$log_gam_b <- replace(upper_bnd$log_gam_b, values = rep(20, length(upper_bnd$log_gam_b)))

  for(sp in 1:data_list$nspp){
    if(data_list$suitMode[sp] %in% c(1:2)) {
      lower_bnd$log_gam_a[sp] <- 1e-10
      upper_bnd$log_gam_a[sp] <- 19.9
    }
  }
  #
  # # Kinzey functional form
  # lower_bnd$logH_3 <- replace(lower_bnd$logH_3, values = rep(-30, length(lower_bnd$logH_3)))
  # upper_bnd$logH_3 <- replace(upper_bnd$logH_3, values = rep(-1e-06, length(upper_bnd$logH_3)))
  #
  # lower_bnd$H_4 <- replace(lower_bnd$H_4, values = rep(-0.1, length(lower_bnd$H_4)))
  # upper_bnd$H_4 <- replace(upper_bnd$H_4, values = rep(20, length(upper_bnd$H_4)))


  # Recruitment ----
  lower_bnd$rec_dev <- replace(lower_bnd$rec_dev, values = rep(-15, length(lower_bnd$rec_dev)))
  upper_bnd$rec_dev <- replace(upper_bnd$rec_dev, values = rep(15, length(upper_bnd$rec_dev)))

  lower_bnd$init_dev <- replace(lower_bnd$init_dev, values = rep(-1000, length(lower_bnd$init_dev)))
  upper_bnd$init_dev <- replace(upper_bnd$init_dev, values = rep(23, length(upper_bnd$init_dev)))


  # # Selectivity ----
  # lower_bnd$log_sel_slp <- replace(lower_bnd$log_sel_slp, values = rep(log(0.01), length(lower_bnd$log_sel_slp)))
  # upper_bnd$log_sel_slp <- replace(upper_bnd$log_sel_slp, values = rep(log(100), length(upper_bnd$log_sel_slp)))
  #
  # lower_bnd$sel_inf <- replace(lower_bnd$sel_inf, values = rep(0, length(lower_bnd$sel_inf)))
  # for(flt in 1:nrow(data_list$fleet_control)){
  #   upper_bnd$sel_inf[,flt,] <- replace(upper_bnd$sel_inf[,flt,], values = rep(data_list$nages[data_list$fleet_control$Species[flt]]+0.5, length(upper_bnd$sel_inf[,flt,])))
  # }


  # # Selectivity deviates ----
  # # - Slope
  # for(i in 1:nrow(data_list$fleet_control)){
  #   # If using blocks don't put bounds on deviates, as these are estimated
  #   if(data_list$fleet_control$Time_varying_sel[i] != 3){
  #     lower_bnd$log_sel_slp_dev[,i,,] <- replace(lower_bnd$log_sel_slp_dev[,i,,], values = rep(-5, length(lower_bnd$log_sel_slp_dev[,i,,])))
  #     upper_bnd$log_sel_slp_dev[,i,,] <- replace(upper_bnd$log_sel_slp_dev[,i,,], values = rep(5, length(upper_bnd$log_sel_slp_dev[,i,,])))
  #   }
  # }
  #
  # # - Asymptotic
  # for(i in 1:nrow(data_list$fleet_control)){
  #   # If using blocks don't put bounds on deviates, as these are estimated
  #   if(data_list$fleet_control$Time_varying_sel[i] != 3){
  #     lower_bnd$sel_inf_dev[,i,,] <- replace(lower_bnd$sel_inf_dev[,i,,], values = rep(-5, length(lower_bnd$sel_inf_dev[,i,,])))
  #     upper_bnd$sel_inf_dev[,i,,] <- replace(upper_bnd$sel_inf_dev[,i,,], values = rep(5, length(upper_bnd$sel_inf_dev[,i,,])))
  #   }
  # }

  # Survey variance ----
  lower_bnd$index_log_sd <- replace(lower_bnd$index_log_sd, values = rep(-10, length(lower_bnd$index_log_sd)))
  upper_bnd$index_log_sd <- replace(upper_bnd$index_log_sd, values = rep(10, length(upper_bnd$index_log_sd)))

  # Fishery variance ----
  lower_bnd$catch_log_sd <- replace(lower_bnd$catch_log_sd, values = rep(-10, length(lower_bnd$catch_log_sd)))
  upper_bnd$catch_log_sd <- replace(upper_bnd$catch_log_sd, values = rep(3, length(upper_bnd$catch_log_sd)))

  # F ----
  lower_bnd$log_F <- replace(lower_bnd$log_F, values = rep(-1000, length(lower_bnd$log_F)))
  upper_bnd$log_F <- replace(upper_bnd$log_F, values = rep(10, length(upper_bnd$log_F)))

  lower_bnd$log_M1 <- replace(lower_bnd$log_M1, values = rep(log(0.001), length(lower_bnd$log_M1)))
  upper_bnd$log_M1 <- replace(upper_bnd$log_M1, values = rep(log(2), length(upper_bnd$log_M1)))

  # Stock-recruit alpha and beta (log scale) ----
  # Overflow guard only: nothing anchors the curve under predation, and on a
  # flat ridge log alpha has reached 702, next to the double-precision limit.
  # +/-30 is 13 orders of magnitude either side of any stock's scale, so it
  # never binds a determined estimate. A linkage bound overrides it below.
  lower_bnd$rec_pars[, 2:3] <- -30
  upper_bnd$rec_pars[, 2:3] <- 30


  # Linkage-table coefficients ----
  # Honor the per-row `lower`/`upper` from data_list$linkage_table.
  # The default Inf/-Inf was set by the generic loop above; here we
  # only override rows that supplied a finite bound. Length-0 vectors
  # (no linkages) are a no-op.
  if (!is.null(data_list$linkage_table) &&
      nrow(data_list$linkage_table) > 0L) {
    tbl <- data_list$linkage_table
    lower_bnd$beta_linkage <- as.numeric(tbl$lower)
    upper_bnd$beta_linkage <- as.numeric(tbl$upper)

    # Intercept beta_linkage entries are mapped out (fixed at 0) and the
    # bound for those rows is enforced on the BASE param below. Loosen the
    # beta_linkage bound for those rows to [-Inf, Inf] so the "inits within
    # bounds" check downstream doesn't trip on the 0 init vs the
    # user-supplied natural-scale bound (e.g. Linf [70, 130]).
    is_int <- tbl$design_col == "(Intercept)"
    lower_bnd$beta_linkage[is_int] <- -Inf
    upper_bnd$beta_linkage[is_int] <- Inf

    # Mirror the init-push at 2-build_params.R:218-256: for (Intercept)
    # rows the beta_linkage coefficient is mapped out and the base
    # parameter (log_growth_pars / growth_log_sd / log_M1 / rec_pars)
    # carries the level. Bounds on the (Intercept) MUST therefore be
    # propagated to the base parameter, otherwise the user's natural-
    # scale bounds (e.g. K in [0.1, 1.0]) have no effect and the
    # optimizer can wander into degenerate regions (K -> 0.07 etc).
    # Inputs are natural-scale; base params are log-scale -> apply log.
    int_rows <- which(tbl$design_col == "(Intercept)" &
                        is.finite(as.numeric(tbl$lower)) &
                        is.finite(as.numeric(tbl$upper)))
    if (length(int_rows) > 0) {
      for (ri in int_rows) {
        row <- tbl[ri, , drop = FALSE]
        idx <- .linkage_row_indices(row, data_list)
        lo  <- as.numeric(row$lower); hi <- as.numeric(row$upper)
        sel_slot <- if (identical(row$process, "sel")) .SEL_PARAM_TO_SLOT[[row$param]] else NULL
        if (!is.null(sel_slot)) {
          .stop_if_shared_block(data_list, row$fleet, "sel", row$param)
          # sel_inf holds a natural-scale inflection only for the logistic
          # family; DoubleNormal stores logit(right_floor) there and LogisticPM
          # a log, so those are refused rather than bounded on the wrong scale.
          if (identical(sel_slot$arr, "sel_inf")) {
            for (f in idx$fleet) {
              .stop_unless_natural_sel_inf(data_list, f, sel_slot$slot, row$param)
            }
          }
        }
        if (identical(row$process, "q")) {
          .stop_if_shared_block(data_list, row$fleet, "q", row$param)
        }
        # Bounds are given on the natural scale. Every base parameter here is
        # stored logged except sel_inf, whose inflections are ages or length
        # midpoints held unlogged -- a lower bound of 0 is ordinary there, and
        # impossible to carry across anywhere else.
        logged   <- !(identical(sel_slot$arr, "sel_inf"))
        if (logged && !(lo > 0)) {
          warning(sprintf("Linkage (Intercept) bound lower = %g <= 0 for ",
                          lo),
                  "process '", row$process, "' param '", row$param,
                  "'; skipping push to log-scale base parameter.",
                  call. = FALSE)
          next
        }
        switch(row$process,
               growth = {
                 mean_idx <- .GROWTH_PARAM_TO_INDEX[row$param]
                 sd_idx   <- .GROWTH_SD_PARAM_TO_INDEX[row$param]
                 if (!is.na(mean_idx)) {
                   for (s in idx$species) {
                     sx <- idx$per_sp[[as.character(s)]]$sex
                     lower_bnd$log_growth_pars[s, sx, mean_idx] <- log(lo)
                     upper_bnd$log_growth_pars[s, sx, mean_idx] <- log(hi)
                   }
                 } else if (!is.na(sd_idx)) {
                   for (s in idx$species) {
                     sx <- idx$per_sp[[as.character(s)]]$sex
                     lower_bnd$growth_log_sd[s, sx, sd_idx] <- log(lo)
                     upper_bnd$growth_log_sd[s, sx, sd_idx] <- log(hi)
                   }
                 }
               },
               M = {
                 for (s in idx$species) {
                   sx <- idx$per_sp[[as.character(s)]]$sex
                   ag <- idx$per_sp[[as.character(s)]]$age
                   lower_bnd$log_M1[s, sx, ag] <- log(lo)
                   upper_bnd$log_M1[s, sx, ag] <- log(hi)
                 }
               },
               recruitment = {
                 par_idx <- .REC_PARAM_TO_INDEX[row$param]
                 if (!is.na(par_idx)) {
                   lower_bnd$rec_pars[idx$species, par_idx] <- log(lo)
                   upper_bnd$rec_pars[idx$species, par_idx] <- log(hi)
                 }
               },
               q = {
                 lower_bnd$index_log_q[idx$fleet] <- log(lo)
                 upper_bnd$index_log_q[idx$fleet] <- log(hi)
               },
               comp = {
                 # DM weight is exp(parameter), so natural-scale bounds log.
                 switch(row$param,
                   theta_comp = {
                     lower_bnd$comp_weights[idx$fleet] <- log(lo)
                     upper_bnd$comp_weights[idx$fleet] <- log(hi)
                   },
                   theta_caal = {
                     lower_bnd$caal_weights[idx$fleet] <- log(lo)
                     upper_bnd$caal_weights[idx$fleet] <- log(hi)
                   },
                   theta_diet = {
                     lower_bnd$diet_comp_weights[idx$species] <- log(lo)
                     upper_bnd$diet_comp_weights[idx$species] <- log(hi)
                   })
               },
               sel = {
                 if (!is.null(sel_slot)) {
                   for (s in idx$species) {
                     sx <- idx$per_sp[[as.character(s)]]$sex
                     if (identical(sel_slot$arr, "log_sel_slp")) {
                       lower_bnd$log_sel_slp[sel_slot$slot, idx$fleet, sx] <- log(lo)
                       upper_bnd$log_sel_slp[sel_slot$slot, idx$fleet, sx] <- log(hi)
                     } else if (identical(sel_slot$arr, "sel_inf")) {
                       lower_bnd$sel_inf[sel_slot$slot, idx$fleet, sx] <- lo
                       upper_bnd$sel_inf[sel_slot$slot, idx$fleet, sx] <- hi
                     } else if (identical(sel_slot$arr, "log_sel_apical")) {
                       lower_bnd$log_sel_apical[idx$fleet, sx] <- log(lo)
                       upper_bnd$log_sel_apical[idx$fleet, sx] <- log(hi)
                     }
                   }
                 }
               }
        )
      }
    }
  }


  # Combine bounds in list ----
  bounds <- list(upper = upper_bnd, lower = lower_bnd)

  # Make sure inits are within bounds ----
  if (sum(unlist(bounds$upper) < as.numeric(unlist(param_list)), na.rm = TRUE) > 0 | sum(as.numeric(unlist(param_list)) <
                                                                                         unlist(bounds$lower), na.rm = TRUE) > 0) {
    lower_check <- param_list
    upper_check <- param_list
    param_check <- data.frame(matrix(NA, nrow = length(param_list), ncol = 3))
    colnames(param_check) <- c("Parameter", "Lower", "Upper")
    param_check$Parameter <- names(param_list)

    for (i in 1:length(param_list)) {
      lower_check[[i]] <- param_list[[i]] < lower_bnd[[i]]
      upper_check[[i]] <- param_list[[i]] > upper_bnd[[i]]
      param_check$Lower[i] <- sum(lower_check[[i]], na.rm = TRUE)
      param_check$Upper[i] <- sum(upper_check[[i]], na.rm = TRUE)
    }

    message("Non-zero value indicates error in initial value")
    print(param_check)
    off <- param_check$Parameter[param_check$Lower > 0 | param_check$Upper > 0]
    # rec_pars is the one block whose bounds arrived after fits were being
    # saved (5.39.0), so a warm start from the old unbounded ridge lands here.
    hint <- if ("rec_pars" %in% off) paste0(
      " Stock-recruit alpha and beta are bounded at +/-30 on the log scale from",
      " 5.39.0, so a fit saved on the old unbounded ridge cannot be refit as it",
      " stands; restart it from values on the stock's scale with",
      " build_srr(srr_alpha_init = , srr_beta_init = ).") else ""
    stop("Initial parameter values are not within bounds: ",
         paste(off, collapse = ", "), ".", hint, call. = FALSE)
  }

  return(bounds)
}
