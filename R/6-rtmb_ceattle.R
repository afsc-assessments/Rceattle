# The RTMB model of CEATTLE: its MakeADFun() wrapper and objective, which holds the
# population dynamics, observation models and likelihood. Population quantities are per
# species, [age, year], in thousands of fish, kg and mt. C++ pointers are to this tree's
# src/TMB (5.67.11): each block ends with the ceattle.cpp section it ports and its lines,
# unchanged since 5.67.3 (the numbering repeats and skips). Only configurations
# .rtmb_check_scope() accepts are ported; the AD mechanics are in 6-rtmb_helpers.R.


#' Build the RTMB objective for a model the C++ template would fit
#'
#' @param data The data list the C++ `MakeADFun()` receives (`fit$obj$env$data`).
#' @param parameters Parameter list. RTMB keeps the order given, where TMB orders
#'   `obj$par` by the template's PARAMETER declarations.
#' @param map,random,silent As for `TMB::MakeADFun()`.
#' @param ... Passed to `RTMB::MakeADFun()`, e.g. `ADreport = TRUE`.
#' @noRd
.rtmb_make_adfun <- function(data, parameters, map = list(), random = NULL,
                             silent = TRUE, ...) {
  if (length(random) > 0) stop("The RTMB model does not yet support random effects.", call. = FALSE)
  pd <- .rtmb_prepare_data(data)
  # Catchability is computed where its starting log q is not NA, as the C++ does
  pd$q_fleets <- which(!is.na(parameters[["index_log_q"]]))
  RTMB::MakeADFun(function(p) .rtmb_objective(p, pd), parameters, map = map,
                  random = random, silent = silent, ...)
}

#' REPORT each element of a list under its own name
#' @noRd
.rtmb_report_list <- function(lst) {
  for (nm in names(lst)) eval(as.call(list(quote(RTMB::REPORT), as.name(nm))), lst)
}

#' Joint negative log-likelihood of the CEATTLE model, with its REPORTs and ADREPORTs
#'
#' @param p Parameter list, as RTMB passes it.
#' @param pd Prepared data, from `.rtmb_prepare_data()`.
#' @noRd
.rtmb_objective <- function(p, pd) {
  "[<-" <- RTMB::ADoverload("[<-")
  c <- RTMB::ADoverload("c")
  nspp <- pd$nspp; n_flt <- pd$n_flt; nages <- pd$nages; max_age <- pd$max_age
  # H hindcast years and Y in all; `hind` and `proj` index the hindcast and projection years
  H <- pd$nyrs_hind; Y <- pd$nyrs; hind <- seq_len(H); proj <- seq_len(Y - H) + H
  spawn_month <- pd$spawn_month
  fishery <- which(pd$flt_type == 1)
  # [species or fleet, sex, age, year] array from per-species or per-fleet [age, year] matrices
  to_array <- function(x, rows, ages = nages[pd$flt_spp]) {
    a <- array(0, c(rows, 1, max_age, Y))
    for (i in seq_len(rows)) if (!is.null(x[[i]])) a[i, 1, seq_len(ages[i]), ] <- x[[i]]
    a
  }

  # -- Parameter transforms; mean recruitment R0 has no linkage here, so it is constant over
  # years (cpp 4, 5.1, 5.6-5.7: 573, 655-656, 745-746, 762-767, 997-1029)
  R_sd <- exp(p[["R_log_sd"]])
  pop_scalar <- exp(p[["log_pop_scalar"]])
  gam_a <- exp(p[["log_gam_a"]])
  gam_b <- exp(p[["log_gam_b"]])
  Flimit <- exp(p[["log_Flimit"]])
  Ftarget <- exp(p[["log_Ftarget"]])
  R0_sp <- exp(p[["rec_pars"]][, 1])
  growth_parameters <- array(0, c(nspp, 1, Y, 4))
  for (k in 1:4) growth_parameters[, 1, , k] <- RTMB::matrix(exp(p[["log_growth_pars"]][, 1, k]), nspp, Y)

  # -- Catchability [fleet, hindcast year]; NA where log q is NA, the NaN the C++ reports
  # there (cpp 5.3, 897-933)
  index_q <- RTMB::matrix(NA_real_, n_flt, H)
  for (f in pd$q_fleets) index_q[f, ] <- exp(p[["index_log_q"]][f] + p[["index_q_dev"]][f, hind])

  # -- Selectivity-at-age (cpp 5.8, 1089-1163)
  selectivity <- .rtmb_selectivity(p, pd)
  sel <- selectivity$sel

  # -- Fishing mortality. Under HCR 0 projected F is 0; Flimit_at_age and Ftarget_at_age
  # scale each fishery's selectivity by its share proj_F_prop (cpp 5.12, 1226-1319)
  proj_F_prop <- p[["proj_F_prop"]]
  F_flt_age <- lapply(seq_len(n_flt), function(f) RTMB::matrix(0, nages[pd$flt_spp[f]], Y))
  F_at_age <- Flimit_at_age <- Ftarget_at_age <- lapply(nages, function(na) RTMB::matrix(0, na, Y))
  F_flt <- RTMB::matrix(0, n_flt, Y)
  F_spp <- RTMB::matrix(0, nspp, Y)
  for (f in fishery) {
    sp <- pd$flt_spp[f]; na <- nages[sp]
    F_hind <- exp(p[["log_F"]][f, hind])
    F_flt_age[[f]] <- sel[[f]][, hind, drop = FALSE] * rep(F_hind, each = na)
    if (Y > H) F_flt_age[[f]] <- cbind(F_flt_age[[f]], RTMB::matrix(0, na, Y - H))
    F_at_age[[sp]] <- F_at_age[[sp]] + F_flt_age[[f]]
    Flimit_at_age[[sp]] <- Flimit_at_age[[sp]] + sel[[f]] * proj_F_prop[f] * Flimit[sp]
    Ftarget_at_age[[sp]] <- Ftarget_at_age[[sp]] + sel[[f]] * proj_F_prop[f] * Ftarget[sp]
    F_flt[f, hind] <- F_hind
    F_spp[sp, hind] <- F_spp[sp, hind] + F_hind
  }

  # -- Natural mortality M1 (log M1 + annual deviation + environmental effect), and the
  # first pass's Z = M1 + F; later passes add predation M2 (cpp 6.1, 1513-1543)
  M1 <- Z1 <- vector("list", nspp)
  for (sp in seq_len(nspp)) {
    na <- nages[sp]; a <- seq_len(na)
    env_effect <- 0
    for (i in seq_len(ncol(pd$env))) {
      env_effect <- env_effect + pd$env[, i] * p[["M1_beta"]][sp, 1, i]
    }
    M1[[sp]] <- exp(RTMB::matrix(p[["log_M1"]][sp, 1, a], na, Y) +
                      RTMB::matrix(p[["log_M1_dev"]][sp, 1, a, ], na, Y) +
                      rep(env_effect, each = na))
    Z1[[sp]] <- M1[[sp]] + F_at_age[[sp]]
  }

  # -- Spawning biomass per recruit (kg) at F = 0, Flimit and Ftarget from the terminal
  # hindcast year's rates, and at Finit (0 under initMode 2) from the first's. Zero under
  # predation: M carries M2, so output per recruit is not the stock's own (cpp 6.2, 1546-1602)
  SPR0 <- SPRlimit <- SPRtarget <- SPRFinit <- numeric(nspp)
  NbyageSPR <- array(0, c(4, nspp, max_age))
  if (pd$msmMode == 0) for (sp in seq_len(nspp)) {
    a <- seq_len(nages[sp])
    Z_F0 <- M1[[sp]][, H]
    Z_Flimit <- Z_F0 + Flimit_at_age[[sp]][, H]
    Z_Ftarget <- Z_F0 + Ftarget_at_age[[sp]][, H]
    Z_Finit <- M1[[sp]][, 1]
    NPR_F0 <- .rtmb_per_recruit_survivors(Z_F0)   # survivors per recruit, by age
    NPR_Flimit <- .rtmb_per_recruit_survivors(Z_Flimit)
    NPR_Ftarget <- .rtmb_per_recruit_survivors(Z_Ftarget)
    NPR_Finit <- .rtmb_per_recruit_survivors(Z_Finit)
    so_H <- pd$spawn_output[[sp]][, H]   # spawning output (kg), last and first hindcast years
    so_1 <- pd$spawn_output[[sp]][, 1]
    SPR0[sp] <- .rtmb_spawning_biomass_per_recruit(NPR_F0, Z_F0, so_H, 1, spawn_month[sp])
    SPRlimit[sp] <- .rtmb_spawning_biomass_per_recruit(NPR_Flimit, Z_Flimit, so_H, 1,
                                                       spawn_month[sp])
    SPRtarget[sp] <- .rtmb_spawning_biomass_per_recruit(NPR_Ftarget, Z_Ftarget, so_H, 1,
                                                        spawn_month[sp])
    SPRFinit[sp] <- .rtmb_spawning_biomass_per_recruit(NPR_Finit, Z_Finit, so_1, 1,
                                                       spawn_month[sp])
    NbyageSPR[1, sp, a] <- NPR_F0
    NbyageSPR[2, sp, a] <- NPR_Flimit
    NbyageSPR[3, sp, a] <- NPR_Ftarget
    NbyageSPR[4, sp, a] <- NPR_Finit
  }
  R_init <- R0_sp
  steepness <- RTMB::matrix(0.99, nspp, Y)

  # -- Population dynamics, run niter times under predation: each pass's M and Z carry the
  # PREVIOUS pass's M2, so the last pass's M2 is reported but never applied, and with
  # niter = 1 none is (cpp 6.1, 6.4-6.11, 7: 1511-2445, 2453-2472)
  msm <- pd$msm
  suit <- predation <- M2 <- M2_dB0 <- M2_dBF <- NULL
  if (!is.null(msm) && pd$suitMode[1] > 0) {
    # Parametric suitability depends on parameters and weight only (cpp 5.10, 1208-1214)
    suit <- .rtmb_suit_parametric(msm, p[["log_phi"]], gam_a, gam_b)
  }
  # Ages one year on: recruits enter the first age, and the plus group keeps its survivors
  age_on <- function(recruits, survivors) {
    na <- length(survivors)
    c(recruits, survivors[seq_len(na - 2)], survivors[na - 1] + survivors[na])
  }
  for (iter in seq_len(pd$niter)) {
    if (is.null(M2)) {
      M <- M_dB0 <- M_dBF <- M1
      Z <- Z1
    } else {
      M <- M_dB0 <- M_dBF <- Z <- vector("list", nspp)
      for (sp in seq_len(nspp)) {
        M[[sp]] <- M1[[sp]] + M2[[sp]]
        M_dB0[[sp]] <- M1[[sp]] + M2_dB0[[sp]]
        M_dBF[[sp]] <- M1[[sp]] + M2_dBF[[sp]]
        Z[[sp]] <- M1[[sp]] + F_at_age[[sp]] + M2[[sp]]
      }
    }

    # Numbers-at-age, recruitment, and the unfished and target-F reference age structures,
    # species by species (cpp 6.4-6.8, 1744-2351)
    N <- avgN <- avgN_dB0 <- avgN_dBF <- NByage0 <- NByageF <- vector("list", nspp)
    R <- RTMB::matrix(0, nspp, Y)
    avg_R <- numeric(nspp)
    mort_sum <- RTMB::matrix(0, nspp, max_age)
    zero_N_pen <- numeric(nspp)   # the likelihood uses the last pass's penalty (cpp:1608)
    biomass <- ssb <- B0 <- SB0 <- SBF <- DynamicB0 <- DynamicSB0 <- DynamicSBF <-
      RTMB::matrix(0, nspp, Y)
    for (sp in seq_len(nspp)) {
      na <- nages[sp]
      mid <- seq_len(na - 2) + 1   # ages between the recruits and the plus group
      rec_dev <- p[["rec_dev"]][sp, ]
      R_hind <- R0_sp[sp] * exp(rec_dev[hind])   # mean recruitment times exp(deviation)

      # Initial age structure (initMode 2, Finit = 0): mean recruitment decayed by cumulative
      # M1, a log deviation per age, geometric plus group; no M2, though SSB reads Z (cpp 6.4)
      m1 <- M1[[sp]][, 1]
      cum_M1 <- c(0, cumsum(m1[seq_len(na - 1)]))
      mort_sum[sp, seq_len(na)] <- cum_M1
      init_dev <- p[["init_dev"]][sp, seq_len(na - 1)]
      N_sp <- RTMB::matrix(0, na, Y)
      N_sp[, 1] <- c(R_hind[1], R0_sp[sp] * exp(-cum_M1[mid] + init_dev[mid - 1]),
                     R0_sp[sp] * exp(-cum_M1[na] + init_dev[na - 1]) / (1 - exp(-m1[na])))

      # Hindcast cells floored at 0.001 by posfun; penalty to "Zero n-at-age penalty" (cpp 6.5)
      for (y in hind[-1]) {
        survivors <- N_sp[, y - 1] * exp(-Z[[sp]][, y - 1])
        floored <- .rtmb_posfun(age_on(R_hind[y], survivors), 0.001)
        N_sp[, y] <- floored$value
        zero_N_pen[sp] <- zero_N_pen[sp] + sum(floored$penalty)
      }
      # Mean recruitment over the first nyrs_srrmean years (cpp 6.6); the projection
      # recruits at that mean (proj_mean_rec = 1), with no floor (cpp 6.8)
      avg_R[sp] <- sum(R_hind[seq_len(pd$nyrs_srrmean)] / pd$nyrs_srrmean)
      for (y in proj) {
        survivors <- N_sp[, y - 1] * exp(-Z[[sp]][, y - 1])
        N_sp[, y] <- age_on(avg_R[sp], survivors)
      }
      R[sp, ] <- c(R_hind, rep(avg_R[sp], length(proj)))
      N[[sp]] <- N_sp
      biomass[sp, ] <- RTMB::colSums(N_sp * pd$W_pop[[sp]])
      ssb[sp, ] <- RTMB::colSums(N_sp * exp(-Z[[sp]] * spawn_month[sp] / 12) *
                                   pd$spawn_output[[sp]])
      avgN[[sp]] <- N_sp * (1 - exp(-Z[[sp]])) / Z[[sp]]

      # NByage0 (F = 0), NByageF (Ftarget): mean recruitment (proj_mean_rec = 1), realized M.
      # Dynamic dB0, dBF: hindcast R, then exp(log mean + deviation); M1 + own M2 (cpp:2035)
      M_sp <- M[[sp]]
      M_dB0_sp <- M_dB0[[sp]]
      M_dBF_sp <- M_dBF[[sp]]
      Ftarget_sp <- Ftarget_at_age[[sp]]
      NByage0_sp <- NByageF_sp <- N_dB0_sp <- N_dBF_sp <- RTMB::matrix(0, na, Y)
      NByage0_sp[, 1] <- N_sp[, 1]; NByageF_sp[, 1] <- N_sp[, 1]
      N_dB0_sp[, 1] <- N_sp[, 1]; N_dBF_sp[, 1] <- N_sp[, 1]
      for (y in seq_len(Y)[-1]) {
        R_dyn <- if (y <= H) R_hind[y] else exp(log(avg_R[sp]) + rec_dev[y])
        surv_F0 <- exp(-M_sp[, y - 1])   # annual survival
        surv_Ftarget <- exp(-M_sp[, y - 1] - Ftarget_sp[, y - 1])
        surv_dB0 <- if (is.null(M2)) surv_F0 else exp(-M_dB0_sp[, y - 1])
        surv_dBF <- if (is.null(M2)) surv_Ftarget else
          exp(-M_dBF_sp[, y - 1] - Ftarget_sp[, y - 1])
        NByage0_sp[, y] <- age_on(avg_R[sp], NByage0_sp[, y - 1] * surv_F0)
        NByageF_sp[, y] <- age_on(avg_R[sp], NByageF_sp[, y - 1] * surv_Ftarget)
        N_dB0_sp[, y] <- age_on(R_dyn, N_dB0_sp[, y - 1] * surv_dB0)
        N_dBF_sp[, y] <- age_on(R_dyn, N_dBF_sp[, y - 1] * surv_dBF)
      }
      NByage0[[sp]] <- NByage0_sp; NByageF[[sp]] <- NByageF_sp
      so_H <- pd$spawn_output[[sp]][, H]   # spawning output (kg) of the last hindcast year
      SB0[sp, ] <- RTMB::colSums(NByage0_sp * so_H * exp(-M_sp * spawn_month[sp] / 12))
      SBF[sp, ] <- RTMB::colSums(NByageF_sp * so_H *
                                   exp(-(M_sp + Ftarget_sp) * spawn_month[sp] / 12))
      DynamicSB0[sp, ] <- RTMB::colSums(N_dB0_sp * pd$spawn_output[[sp]] *
                                          exp(-M_dB0_sp * spawn_month[sp] / 12))
      DynamicSBF[sp, ] <- RTMB::colSums(N_dBF_sp * pd$spawn_output[[sp]] *
                                          exp(-(M_dBF_sp + Ftarget_sp) * spawn_month[sp] / 12))
      B0[sp, ] <- RTMB::colSums(NByage0_sp * pd$W_pop[[sp]][, H])
      DynamicB0[sp, ] <- RTMB::colSums(N_dB0_sp * pd$W_pop[[sp]])
      if (!is.null(msm)) {
        # Average numbers of the dynamic reference structures (cpp 6.11, 2440-2441)
        avgN_dB0[[sp]] <- N_dB0_sp * (1 - exp(-M_dB0_sp)) / M_dB0_sp
        avgN_dBF[[sp]] <- N_dBF_sp * (1 - exp(-M_dBF_sp - Ftarget_sp)) /
          (M_dBF_sp + Ftarget_sp)
      }
    }

    # Predation: suitability (empirical: from this pass's numbers), then the M2 that the
    # next pass applies (cpp 7, 2453-2472)
    if (!is.null(msm)) {
      avgN_cells <- do.call(rbind, avgN)   # [cell, year], cells species then age
      if (pd$suitMode[1] == 0) suit <- .rtmb_suit_empirical(msm, avgN_cells)
      predation <- .rtmb_predation(msm, suit$suitability, suit$suit_other, avgN_cells,
                                   do.call(rbind, avgN_dB0), do.call(rbind, avgN_dBF))
      M2 <- lapply(msm$cell_rows, function(i) predation$M2[i, , drop = FALSE])
      M2_dB0 <- lapply(msm$cell_rows, function(i) predation$M2_dB0[i, , drop = FALSE])
      M2_dBF <- lapply(msm$cell_rows, function(i) predation$M2_dBF[i, , drop = FALSE])
    }
  }
  if (!is.null(msm)) {
    # Under predation the static references are the input MSSB0 and MSB0 (cpp:2179-2182)
    SB0 <- RTMB::matrix(msm$MSSB0, nspp, Y)
    B0 <- RTMB::matrix(msm$MSB0, nspp, Y)
  }
  R0 <- RTMB::matrix(R0_sp, nspp, Y)
  R_hat <- R0   # expected recruitment from the mean, with no deviation (cpp 6.9)

  # -- Survey index: numbers at the survey month times selectivity and weight (kg), then
  # catchability (cpp 8.1, 8.3: 2481-2534, 2654-2671)
  index_hat <- numeric(pd$n_index)
  for (f in seq_len(n_flt)) {
    ix <- pd$index[[f]]
    if (is.null(ix)) next
    sp <- pd$flt_spp[f]; na <- nages[sp]; yrs <- ix$yr
    N_survey <- N[[sp]][, yrs, drop = FALSE] *
      exp(Z[[sp]][, yrs, drop = FALSE] * rep(-(ix$month / 12), each = na))
    survey_biomass <- RTMB::colSums(N_survey * sel[[f]][, yrs, drop = FALSE] *
                                      pd$W_flt[[f]][, yrs, drop = FALSE])
    index_hat[ix$rows] <- index_q[f, ix$q_yr] * survey_biomass
  }

  # -- Baranov catch by weight (mt); its concentrated sd (Walters and Ludwig 1994) is
  # reported, not fitted, while the input sd is used (cpp 9.1, 9.1b: 2705-2798)
  catch_hat <- max_catch_hat <- numeric(pd$n_catch)
  catch_analytical_sd <- numeric(n_flt)
  for (f in seq_len(n_flt)) {
    cx <- pd$catch[[f]]
    if (is.null(cx)) next
    sp <- pd$flt_spp[f]; yrs <- cx$yr
    Z_yrs <- Z[[sp]][, yrs, drop = FALSE]
    N_yrs <- N[[sp]][, yrs, drop = FALSE]
    W_yrs <- pd$W_flt[[f]][, yrs, drop = FALSE]
    catch_wt <- RTMB::colSums(F_flt_age[[f]][, yrs, drop = FALSE] / Z_yrs * (1 - exp(-Z_yrs)) *
                                N_yrs * W_yrs)
    catch_hat[cx$rows] <- catch_wt
    max_catch_hat[cx$rows] <- RTMB::colSums(N_yrs * W_yrs * sel[[f]][, yrs, drop = FALSE] *
                                              proj_F_prop[f])
    if (any(cx$analytical)) {
      resid <- cx$log_obs_analytical - log(catch_wt[cx$analytical])
      catch_analytical_sd[f] <- .rtmb_concentrated_lognormal_sd(
        sum(resid * resid) / sum(cx$analytical), pd$bias_adjust_obs)
    }
  }

  # -- Exploitable biomass (mt), weighted by each fishery's share (cpp 9.2, 2801-2818)
  exploitable_biomass <- RTMB::matrix(0, nspp, Y)
  for (f in fishery) {
    sp <- pd$flt_spp[f]
    exploitable_biomass[sp, ] <- exploitable_biomass[sp, ] +
      RTMB::colSums(N[[sp]] * pd$W_flt[[f]] * sel[[f]] * proj_F_prop[f])
  }

  # -- Compositions: numbers at true age through ageing error (age) or the age-length key
  # (length), normalized; multinomial likelihood (cpp 10.1, 13.4: 2825-3135, 3613-3692)
  age_hat <- age_obs_hat <- RTMB::matrix(0, pd$n_comp, pd$max_age_cols)
  comp_hat <- RTMB::matrix(0, pd$n_comp, ncol(pd$comp_obs))
  comp_nll <- comp_nll_unweighted <- numeric(n_flt)
  for (g in pd$comp) {
    f <- g$flt; sp <- g$sp; na <- nages[sp]; yrs <- g$yr
    N_yrs <- N[[sp]][, yrs, drop = FALSE]
    Z_yrs <- Z[[sp]][, yrs, drop = FALSE]
    # Catch (fishery) or survey numbers at true age, [age, row]
    n_true_age <- if (pd$flt_type[f] == 1) {
      F_flt_age[[f]][, yrs, drop = FALSE] / Z_yrs * (1 - exp(-Z_yrs)) * N_yrs
    } else {
      N_yrs * sel[[f]][, yrs, drop = FALSE] * rep(index_q[f, g$q_yr], each = na) *
        exp(Z_yrs * rep(-(g$month / 12), each = na))
    }
    n_obs_age <- g$age_error_t %*% n_true_age
    comp_pred <- if (g$type == 0) n_obs_age else g$age_length_t %*% n_true_age
    comp_pred <- comp_pred /
      rep(.rtmb_posfun(RTMB::colSums(comp_pred), 0.001)$value, each = g$nb)
    age_hat[g$rows, seq_len(na)] <- t(n_true_age)
    age_obs_hat[g$rows, seq_len(na)] <- t(n_obs_age)
    comp_hat[g$rows, seq_len(g$nb)] <- t(comp_pred)
    if (any(g$fit)) {
      # Multinomial with the offset on both proportions; counts are fractional
      hp <- comp_pred[, g$fit, drop = FALSE] + pd$comp_offset
      prob <- hp / rep(RTMB::colSums(hp) + 1e-300, each = g$nb)
      ll <- g$ll_const + RTMB::colSums(g$obs_n * log(prob))
      comp_nll[f] <- comp_nll[f] - sum(p[["comp_weights"]][f] * ll)
      comp_nll_unweighted[f] <- comp_nll_unweighted[f] - sum(ll)
    }
  }

  # -- Depletion against the final-year static reference; under predation with no harvest
  # rule, against the final (unfished) projection year (cpp 12.1, 3245-3270)
  if (is.null(msm)) {
    biomass_depletion <- biomass / B0[, Y]
    ssb_depletion <- ssb / SB0[, Y]
  } else {
    biomass_depletion <- biomass / biomass[, Y]
    ssb_depletion <- ssb / ssb[, Y]
  }

  # -- Negative log-likelihood components: rows are addressed by the JnllRow labels in
  # .JNLL_ROW_AXIS, and columns count fleets or species by row (cpp 13, 3279-4658)
  jrow <- stats::setNames(seq_along(.JNLL_ROW_AXIS), names(.JNLL_ROW_AXIS))
  jnll_comp <- RTMB::matrix(0, length(jrow), pd$n_col)
  for (f in seq_len(n_flt)) {
    # Lognormal index, mean bias-corrected by -sd^2/2 (cpp 13.1)
    ix <- pd$index[[f]]
    if (!is.null(ix) && any(ix$fit)) {
      mu <- log(index_hat[ix$rows[ix$fit]]) - ix$bias
      jnll_comp[jrow[["Index data"]], f] <- -sum(RTMB::dnorm(ix$obs, mu, ix$sd[ix$fit], log = TRUE))
    }
    # Lognormal catch (cpp 13.2)
    cx <- pd$catch[[f]]
    if (!is.null(cx) && any(cx$fit)) {
      mu <- log(catch_hat[cx$rows[cx$fit]]) - cx$bias
      jnll_comp[jrow[["Catch data"]], f] <- -sum(RTMB::dnorm(cx$obs, mu, cx$sd[cx$fit], log = TRUE))
    }
  }
  jnll_comp[jrow[["Composition data"]], seq_len(n_flt)] <- comp_nll
  # NonParametric shape penalties on the first year's log curve (Sel_curve_pen1 on any
  # decline between adjacent bins, Sel_curve_pen2 on curvature) and on its level avg_sel,
  # counted once per shared selectivity block (cpp 13.7, 3892-3960)
  for (f in which(pd$flt_sel_type == 2 & pd$flt_type > 0 & pd$flt_sel_lead == 1)) {
    log_sel <- selectivity$log_sel_yr1[[f]]
    nb <- length(log_sel)
    fall <- log_sel[-nb] - log_sel[-1]        # from each bin to the next; negative for a rise
    decline <- (abs(fall) + fall) / 2.0   # the fall where positive, else 0
    curvature <- diff(diff(log_sel))
    jnll_comp[jrow[["Non-parametric selectivity"]], f] <-
      sum(p[["sel_curve_pen"]][f, 1] * (decline * decline)) +
      sum(p[["sel_curve_pen"]][f, 2] * curvature * curvature) +
      2.0 * (selectivity$avg_sel[[f]] * selectivity$avg_sel[[f]])
  }
  # Initial-age and recruitment deviations, N(-sigma_R^2/2, sigma_R) under bias_adjust_proc
  # (cpp 13.9, 4253-4297); and the zero-N floor penalty (cpp 13.11, 4384-4390)
  for (sp in seq_len(nspp)) {
    mu <- -pd$bias_adjust_proc * (R_sd[sp] * R_sd[sp]) / 2.0
    jnll_comp[jrow[["Initial abundance deviates"]], sp] <-
      -sum(RTMB::dnorm(p[["init_dev"]][sp, seq_len(nages[sp] - 1)], mu, R_sd[sp], log = TRUE))
    jnll_comp[jrow[["Recruitment deviates"]], sp] <-
      -sum(RTMB::dnorm(p[["rec_dev"]][sp, hind], mu, R_sd[sp], log = TRUE))
    jnll_comp[jrow[["Zero n-at-age penalty"]], sp] <- zero_N_pen[sp]
  }
  # Predicted stomach content from the last pass, and its likelihood, where suitability is
  # estimated (cpp 11, 13.6: 3233-3237, 3753-3817)
  diet_hat <- stomach <- NULL
  if (!is.null(msm) && pd$suitMode[1] > 0) {
    diet_hat <- .rtmb_diet_hat(msm, suit$suitability, avgN_cells, predation$inv_avail)
    stomach <- .rtmb_stomach_nll(msm, diet_hat, p[["diet_comp_weights"]])
    jnll_comp[jrow[["Stomach content data"]], seq_len(nspp)] <- stomach$weighted
  }
  # Unweighted copy; the composition and stomach rows differ (cpp 13.15, 4646-4658)
  unweighted_jnll_comp <- jnll_comp
  unweighted_jnll_comp[jrow[["Composition data"]], seq_len(n_flt)] <- comp_nll_unweighted
  if (!is.null(stomach)) {
    unweighted_jnll_comp[jrow[["Stomach content data"]], seq_len(nspp)] <- stomach$unweighted
  }
  jnll <- sum(jnll_comp)

  # -- Reports, at the C++ names and shapes (cpp 14, 4661-4887)
  N_at_age <- to_array(N, nspp, nages)
  avgN_at_age <- to_array(avgN, nspp, nages)
  M1_at_age <- to_array(M1, nspp, nages)
  M_at_age <- to_array(M, nspp, nages)
  Z_at_age <- to_array(Z, nspp, nages)
  F_at_age <- to_array(F_at_age, nspp, nages)
  NByage0 <- to_array(NByage0, nspp, nages)
  NByageF <- to_array(NByageF, nspp, nages)
  biomass_at_age <- to_array(lapply(seq_len(nspp), function(sp) N[[sp]] * pd$W_pop[[sp]]),
                             nspp, nages)
  F_flt_age <- to_array(F_flt_age, n_flt)
  sel_at_age <- to_array(sel, n_flt)
  log_biomass <- log(biomass)
  log_ssb <- log(ssb)
  log_R <- log(R)
  log_index_hat <- log(index_hat)
  beta_linkage <- p[["beta_linkage"]]
  beta_linkage_re <- p[["beta_linkage_re"]]
  beta_linkage_re_pen <- p[["beta_linkage_re_pen"]]
  beta_linkage_obs <- p[["beta_linkage_obs"]]
  beta_linkage_re_all <- numeric(0)
  .rtmb_report_list(pd$fixed_reports)
  if (!is.null(msm) && !inherits(jnll, "advector")) {
    # Built in the report() pass only, which keeps the per-pair products off the tape
    vul <- if (is.null(suit$vulnerability)) matrix(0, nspp, nspp) else suit$vulnerability
    vul_other <- if (is.null(suit$vulnerability_other)) numeric(nspp) else
      suit$vulnerability_other
    .rtmb_report_list(.rtmb_predation_reports(msm, suit$suitability, suit$suit_other, vul,
                                              vul_other, predation, avgN_cells, diet_hat))
  }
  .rtmb_report_list(list(
    pop_scalar = pop_scalar, avg_R = avg_R, Z_at_age = Z_at_age, N_at_age = N_at_age,
    avgN_at_age = avgN_at_age, biomass_at_age = biomass_at_age,
    biomass_depletion = biomass_depletion, ssb_depletion = ssb_depletion,
    biomass = biomass, ssb = ssb, exploitable_biomass = exploitable_biomass, R_sd = R_sd,
    R0 = R0, R_init = R_init, R = R, R_hat = R_hat, M_at_age = M_at_age,
    log_biomass = log_biomass, log_ssb = log_ssb, log_R = log_R,
    NByageF = NByageF, NByage0 = NByage0, NbyageSPR = NbyageSPR, B0 = B0, SB0 = SB0,
    SBF = SBF, DynamicB0 = DynamicB0, DynamicSB0 = DynamicSB0, DynamicSBF = DynamicSBF,
    SPR0 = SPR0, SPRFinit = SPRFinit, SPRlimit = SPRlimit, SPRtarget = SPRtarget,
    steepness = steepness, Ftarget = Ftarget, Flimit = Flimit, sel_at_age = sel_at_age,
    index_hat = index_hat, index_q = index_q, log_index_hat = log_index_hat,
    beta_linkage = beta_linkage, beta_linkage_re = beta_linkage_re,
    beta_linkage_re_pen = beta_linkage_re_pen, beta_linkage_re_all = beta_linkage_re_all,
    beta_linkage_obs = beta_linkage_obs, F_spp = F_spp, F_flt = F_flt,
    F_flt_age = F_flt_age, F_at_age = F_at_age, catch_hat = catch_hat,
    max_catch_hat = max_catch_hat, catch_analytical_sd = catch_analytical_sd,
    age_obs_hat = age_obs_hat, age_hat = age_hat, comp_hat = comp_hat,
    jnll_comp = jnll_comp, unweighted_jnll_comp = unweighted_jnll_comp,
    M1_at_age = M1_at_age, gam_a = gam_a, gam_b = gam_b, mort_sum = mort_sum,
    growth_parameters = growth_parameters, jnll = jnll))
  # ADREPORTs in the C++ order, so the flat vectors line up
  RTMB::ADREPORT(pop_scalar)
  RTMB::ADREPORT(biomass)
  RTMB::ADREPORT(ssb)
  RTMB::ADREPORT(exploitable_biomass)
  RTMB::ADREPORT(biomass_depletion)
  RTMB::ADREPORT(ssb_depletion)
  RTMB::ADREPORT(log_biomass)
  RTMB::ADREPORT(log_ssb)
  RTMB::ADREPORT(log_R)
  RTMB::ADREPORT(R_sd)
  RTMB::ADREPORT(R)
  RTMB::ADREPORT(log_index_hat)
  RTMB::ADREPORT(beta_linkage)
  RTMB::ADREPORT(beta_linkage_obs)
  jnll
}
