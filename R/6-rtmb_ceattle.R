# The RTMB model of CEATTLE: its MakeADFun() wrapper and objective. Each block names the
# src/TMB/ceattle.cpp section it ports, with that section's line range on 5.67.5: the
# file's own numbering repeats and skips in places, so the line range is the reliable
# pointer. Only configurations .rtmb_check_scope() accepts are ported (stage 1).


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

#' Negative log-likelihood of the CEATTLE model
#'
#' @param p Parameter list, as RTMB passes it.
#' @param pd Prepared data, from `.rtmb_prepare_data()`.
#' @noRd
.rtmb_objective <- function(p, pd) {
  "[<-" <- RTMB::ADoverload("[<-")
  c <- RTMB::ADoverload("c")
  # Between operands of equal length RTMB takes dim() from the first, so a matrix leads each
  # product; RTMB's matrix() and colSums() replace base's only when attached, so name them.
  nspp <- pd$nspp; n_flt <- pd$n_flt; nages <- pd$nages; max_age <- pd$max_age
  H <- pd$nyrs_hind; Y <- pd$nyrs; hind <- seq_len(H); proj <- seq_len(Y - H) + H
  sm <- pd$spawn_month
  fishery <- which(pd$flt_type == 1)
  # [species or fleet, sex, age, year] array from per-species or per-fleet [age, year] matrices
  to_array <- function(x, rows, ages = nages[pd$flt_spp]) {
    a <- array(0, c(rows, 1, max_age, Y))
    for (i in seq_len(rows)) if (!is.null(x[[i]])) a[i, 1, seq_len(ages[i]), ] <- x[[i]]
    a
  }

  # -- 4 (cpp:573, 655-656, 745-746), 5.1 (cpp:762-767), 5.6-5.7 (cpp:997-1029): parameter
  # transforms. Mean recruitment R0 has no linkage here, so it is constant over years.
  R_sd <- exp(p[["R_log_sd"]])
  pop_scalar <- exp(p[["log_pop_scalar"]])
  gam_a <- exp(p[["log_gam_a"]])
  gam_b <- exp(p[["log_gam_b"]])
  Flimit <- exp(p[["log_Flimit"]])
  Ftarget <- exp(p[["log_Ftarget"]])
  R0_sp <- exp(p[["rec_pars"]][, 1])
  growth_parameters <- array(0, c(nspp, 1, Y, 4))
  for (k in 1:4) growth_parameters[, 1, , k] <- RTMB::matrix(exp(p[["log_growth_pars"]][, 1, k]), nspp, Y)

  # -- 5.3 (cpp:897-933): catchability [fleet, hindcast year]; NA where log q is NA, the
  # NaN the C++ reports there
  index_q <- RTMB::matrix(NA_real_, n_flt, H)
  for (f in pd$q_fleets) index_q[f, ] <- exp(p[["index_log_q"]][f] + p[["index_q_dev"]][f, hind])

  # -- 5.8 (cpp:1089-1163): selectivity-at-age
  S <- .rtmb_selectivity(p, pd)
  sel <- S$sel

  # -- 5.12 (cpp:1226-1319): fishing mortality. Under HCR 0 projected F is 0; Flimit_at_age
  # and Ftarget_at_age scale each fishery's selectivity by its share proj_F_prop.
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

  # -- 6.1 (cpp:1513-1543): natural mortality M1 (log M1 + annual deviation + environmental
  # effect), and the first pass's Z = M1 + F; later passes add predation M2 (below)
  M1 <- Z1 <- vector("list", nspp)
  for (sp in seq_len(nspp)) {
    na <- nages[sp]; a <- seq_len(na)
    mult <- 0
    for (i in seq_len(ncol(pd$env))) mult <- mult + pd$env[, i] * p[["M1_beta"]][sp, 1, i]
    M1[[sp]] <- exp(RTMB::matrix(p[["log_M1"]][sp, 1, a], na, Y) +
                      RTMB::matrix(p[["log_M1_dev"]][sp, 1, a, ], na, Y) + rep(mult, each = na))
    Z1[[sp]] <- M1[[sp]] + F_at_age[[sp]]
  }

  # -- 6.2 (cpp:1546-1602): spawning biomass per recruit at F = 0, Flimit and Ftarget from
  # the terminal hindcast year's rates, and at Finit (0 under initMode 2) from the first's.
  # Zero under predation: M carries M2, so output per recruit is not the stock's own.
  SPR0 <- SPRlimit <- SPRtarget <- SPRFinit <- numeric(nspp)
  NbyageSPR <- array(0, c(4, nspp, max_age))
  if (pd$msmMode == 0) for (sp in seq_len(nspp)) {
    a <- seq_len(nages[sp])
    Zu <- M1[[sp]][, H]
    Zl <- Zu + Flimit_at_age[[sp]][, H]
    Zt <- Zu + Ftarget_at_age[[sp]][, H]
    Zi <- M1[[sp]][, 1]
    nu <- .rtmb_per_recruit_survivors(Zu)
    nl <- .rtmb_per_recruit_survivors(Zl)
    nt <- .rtmb_per_recruit_survivors(Zt)
    ni <- .rtmb_per_recruit_survivors(Zi)
    w_term <- pd$so[[sp]][, H]
    SPR0[sp] <- .rtmb_spawning_biomass_per_recruit(nu, Zu, w_term, 1, sm[sp])
    SPRlimit[sp] <- .rtmb_spawning_biomass_per_recruit(nl, Zl, w_term, 1, sm[sp])
    SPRtarget[sp] <- .rtmb_spawning_biomass_per_recruit(nt, Zt, w_term, 1, sm[sp])
    SPRFinit[sp] <- .rtmb_spawning_biomass_per_recruit(ni, Zi, pd$so[[sp]][, 1], 1, sm[sp])
    NbyageSPR[1, sp, a] <- nu
    NbyageSPR[2, sp, a] <- nl
    NbyageSPR[3, sp, a] <- nt
    NbyageSPR[4, sp, a] <- ni
  }
  R_init <- R0_sp
  steepness <- RTMB::matrix(0.99, nspp, Y)

  # -- 6.1 and 6.4-6.11 (cpp:1511-2445), then 7 (cpp:2453-2472) under predation: the
  # dynamics run niter times. Each pass's M and Z carry the PREVIOUS pass's M2, so the last
  # pass's M2 is reported but never applied (and with niter = 1, no M2 is applied).
  msm <- pd$msm
  su <- pred <- M2 <- M2_dB0 <- M2_dBF <- NULL
  if (!is.null(msm) && pd$suitMode[1] > 0) {
    # 5.10 (cpp:1208-1214): parametric suitability depends on parameters and weight only
    su <- .rtmb_suit_parametric(msm, p[["log_phi"]], gam_a, gam_b)
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

    # 6.4-6.8 (cpp:1744-2351): numbers-at-age (thousands of fish), recruitment, and the
    # unfished and target-F reference age structures, species by species
    N <- avgN <- avgN_dB0 <- avgN_dBF <- NB0 <- NBF <- dB0 <- dBF <- vector("list", nspp)
    R <- RTMB::matrix(0, nspp, Y)
    avg_R <- numeric(nspp)
    mort_sum <- RTMB::matrix(0, nspp, max_age)
    zero_N_pen <- numeric(nspp)   # each pass starts again (cpp:1608)
    biomass <- ssb <- B0 <- SB0 <- SBF <- DynamicB0 <- DynamicSB0 <- DynamicSBF <-
      RTMB::matrix(0, nspp, Y)
    for (sp in seq_len(nspp)) {
      na <- nages[sp]
      mid <- seq_len(na - 2) + 1   # ages between the recruits and the plus group
      rec_dev <- p[["rec_dev"]][sp, ]
      # Hindcast recruitment, mean recruitment times exp(deviation)
      Rh <- R0_sp[sp] * exp(rec_dev[hind])

      # 6.4: initial age structure under initMode 2: mean recruitment decayed by
      # cumulative M1 (Finit = 0), with a log deviation per age; the plus group sums the
      # geometric tail. It ignores M2, while the first year's SSB below reads Z.
      m1 <- M1[[sp]][, 1]
      ms <- c(0, cumsum(m1[seq_len(na - 1)]))
      mort_sum[sp, seq_len(na)] <- ms
      idev <- p[["init_dev"]][sp, seq_len(na - 1)]
      Nm <- RTMB::matrix(0, na, Y)
      Nm[, 1] <- c(Rh[1], R0_sp[sp] * exp(-ms[mid] + idev[mid - 1]),
                   R0_sp[sp] * exp(-ms[na] + idev[na - 1]) / (1 - exp(-m1[na])))

      # 6.5: hindcast; every cell is floored at 0.001 by posfun, its penalty to ZERO_N
      for (y in hind[-1]) {
        s <- Nm[, y - 1] * exp(-Z[[sp]][, y - 1])
        fl <- .rtmb_posfun(c(Rh[y], s[mid - 1], s[na - 1] + s[na]), 0.001)
        Nm[, y] <- fl$value
        zero_N_pen[sp] <- zero_N_pen[sp] + sum(fl$penalty)
      }
      # 6.6: mean recruitment over the first nyrs_srrmean years; 6.8: the projection
      # recruits at that mean (proj_mean_rec = 1), with no floor
      avg_R[sp] <- sum(Rh[seq_len(pd$nyrs_srrmean)] / pd$nyrs_srrmean)
      for (y in proj) {
        s <- Nm[, y - 1] * exp(-Z[[sp]][, y - 1])
        Nm[, y] <- c(avg_R[sp], s[mid - 1], s[na - 1] + s[na])
      }
      R[sp, ] <- c(Rh, rep(avg_R[sp], length(proj)))
      N[[sp]] <- Nm
      biomass[sp, ] <- RTMB::colSums(Nm * pd$W_pop[[sp]])
      ssb[sp, ] <- RTMB::colSums(Nm * exp(-Z[[sp]] * sm[sp] / 12) * pd$so[[sp]])
      avgN[[sp]] <- Nm * (1 - exp(-Z[[sp]])) / Z[[sp]]

      # 6.6 (cpp:1999-2184), Option 1a: unfished (NByage0) and target-F (NByageF) age
      # structures recruiting at the mean, on the realized M; and their dynamic versions
      # recruiting at the hindcast R, then at exp(log mean + deviation) in the projection,
      # on M1 plus their own M2 (M_dB0, M_dBF)
      Mm <- M[[sp]]
      M0 <- M_dB0[[sp]]
      MF <- M_dBF[[sp]]
      Ft <- Ftarget_at_age[[sp]]
      b0 <- bF <- d0 <- dF <- RTMB::matrix(0, na, Y)
      b0[, 1] <- Nm[, 1]; bF[, 1] <- Nm[, 1]; d0[, 1] <- Nm[, 1]; dF[, 1] <- Nm[, 1]
      for (y in seq_len(Y)[-1]) {
        r_dyn <- if (y <= H) Rh[y] else exp(log(avg_R[sp]) + rec_dev[y])
        e0 <- exp(-Mm[, y - 1])
        eF <- exp(-Mm[, y - 1] - Ft[, y - 1])
        e0d <- if (is.null(M2)) e0 else exp(-M0[, y - 1])
        eFd <- if (is.null(M2)) eF else exp(-MF[, y - 1] - Ft[, y - 1])
        s <- b0[, y - 1] * e0; b0[, y] <- c(avg_R[sp], s[mid - 1], s[na - 1] + s[na])
        s <- bF[, y - 1] * eF; bF[, y] <- c(avg_R[sp], s[mid - 1], s[na - 1] + s[na])
        s <- d0[, y - 1] * e0d; d0[, y] <- c(r_dyn, s[mid - 1], s[na - 1] + s[na])
        s <- dF[, y - 1] * eFd; dF[, y] <- c(r_dyn, s[mid - 1], s[na - 1] + s[na])
      }
      NB0[[sp]] <- b0; NBF[[sp]] <- bF; dB0[[sp]] <- d0; dBF[[sp]] <- dF
      so_term <- pd$so[[sp]][, H]
      SB0[sp, ] <- RTMB::colSums(b0 * so_term * exp(-Mm * sm[sp] / 12))
      SBF[sp, ] <- RTMB::colSums(bF * so_term * exp(-(Mm + Ft) * sm[sp] / 12))
      DynamicSB0[sp, ] <- RTMB::colSums(d0 * pd$so[[sp]] * exp(-M0 * sm[sp] / 12))
      DynamicSBF[sp, ] <- RTMB::colSums(dF * pd$so[[sp]] * exp(-(MF + Ft) * sm[sp] / 12))
      B0[sp, ] <- RTMB::colSums(b0 * pd$W_pop[[sp]][, H])
      DynamicB0[sp, ] <- RTMB::colSums(d0 * pd$W_pop[[sp]])
      if (!is.null(msm)) {
        # 6.11 (cpp:2440-2441): average numbers of the dynamic reference structures
        avgN_dB0[[sp]] <- d0 * (1 - exp(-M0)) / M0
        avgN_dBF[[sp]] <- dF * (1 - exp(-MF - Ft)) / (MF + Ft)
      }
    }

    # 7 (cpp:2453-2472): predation, with suitability (empirical: from this pass's numbers),
    # giving the M2 that the next pass applies
    if (!is.null(msm)) {
      avgN_c <- do.call(rbind, avgN)
      if (pd$suitMode[1] == 0) su <- .rtmb_suit_empirical(msm, avgN_c)
      pred <- .rtmb_predation(msm, su$S, su$suit_other, avgN_c, do.call(rbind, avgN_dB0),
                              do.call(rbind, avgN_dBF))
      M2 <- lapply(msm$cell_rows, function(i) pred$M2[i, , drop = FALSE])
      M2_dB0 <- lapply(msm$cell_rows, function(i) pred$M2_dB0[i, , drop = FALSE])
      M2_dBF <- lapply(msm$cell_rows, function(i) pred$M2_dBF[i, , drop = FALSE])
    }
  }
  if (!is.null(msm)) {
    # cpp:2179-2182: under predation the static references are the input MSSB0 and MSB0
    SB0 <- RTMB::matrix(msm$MSSB0, nspp, Y)
    B0 <- RTMB::matrix(msm$MSB0, nspp, Y)
  }
  R0 <- RTMB::matrix(R0_sp, nspp, Y)
  R_hat <- R0   # 6.9: expected recruitment from the mean, with no deviation

  # -- 8.1, 8.3 (cpp:2481-2534, 2654-2671): survey index: numbers at the survey month
  # times selectivity and weight (kg), then catchability
  index_hat <- numeric(pd$n_index)
  for (f in seq_len(n_flt)) {
    ix <- pd$index[[f]]
    if (is.null(ix)) next
    sp <- pd$flt_spp[f]; na <- nages[sp]; cl <- ix$col
    n_index <- N[[sp]][, cl, drop = FALSE] *
      exp(Z[[sp]][, cl, drop = FALSE] * rep(-(ix$mo / 12), each = na))
    raw <- RTMB::colSums(n_index * sel[[f]][, cl, drop = FALSE] * pd$W_flt[[f]][, cl, drop = FALSE])
    index_hat[ix$rows] <- index_q[f, ix$qcol] * raw
  }

  # -- 9.1, 9.1b (cpp:2705-2798): Baranov catch by weight (mt); its concentrated sd
  # (Walters and Ludwig 1994) is reported, not fitted, while the input sd is used
  catch_hat <- max_catch_hat <- numeric(pd$n_catch)
  catch_analytical_sd <- numeric(n_flt)
  for (f in seq_len(n_flt)) {
    cx <- pd$catch[[f]]
    if (is.null(cx)) next
    sp <- pd$flt_spp[f]; cl <- cx$col
    Zc <- Z[[sp]][, cl, drop = FALSE]
    Nc <- N[[sp]][, cl, drop = FALSE]
    Wc <- pd$W_flt[[f]][, cl, drop = FALSE]
    ch <- RTMB::colSums(F_flt_age[[f]][, cl, drop = FALSE] / Zc * (1 - exp(-Zc)) * Nc * Wc)
    catch_hat[cx$rows] <- ch
    max_catch_hat[cx$rows] <- RTMB::colSums(Nc * Wc * sel[[f]][, cl, drop = FALSE] * proj_F_prop[f])
    if (any(cx$analytical)) {
      r <- cx$log_obs_analytical - log(ch[cx$analytical])
      catch_analytical_sd[f] <- .rtmb_concentrated_lognormal_sd(
        sum(r * r) / sum(cx$analytical), pd$bias_adjust_obs)
    }
  }

  # -- 9.2 (cpp:2801-2818): exploitable biomass, weighted by each fishery's share
  exploitable_biomass <- RTMB::matrix(0, nspp, Y)
  for (f in fishery) {
    sp <- pd$flt_spp[f]
    exploitable_biomass[sp, ] <- exploitable_biomass[sp, ] +
      RTMB::colSums(N[[sp]] * pd$W_flt[[f]] * sel[[f]] * proj_F_prop[f])
  }

  # -- 10.1 (cpp:2825-3135), 13.4 (cpp:3613-3692): age and length compositions, and
  # their multinomial likelihood. Predicted catch-at-true-age goes through the ageing-
  # error matrix (age comps) or the age-length key (length comps), then is normalized.
  age_hat <- age_obs_hat <- RTMB::matrix(0, pd$n_comp, pd$max_age_cols)
  comp_hat <- RTMB::matrix(0, pd$n_comp, ncol(pd$comp_obs))
  comp_nll <- comp_nll_unweighted <- numeric(n_flt)
  for (g in pd$comp) {
    f <- g$flt; sp <- g$sp; na <- nages[sp]; cl <- g$col
    Nc <- N[[sp]][, cl, drop = FALSE]
    Zc <- Z[[sp]][, cl, drop = FALSE]
    ah <- if (pd$flt_type[f] == 1) {
      F_flt_age[[f]][, cl, drop = FALSE] / Zc * (1 - exp(-Zc)) * Nc
    } else {
      Nc * sel[[f]][, cl, drop = FALSE] * rep(index_q[f, g$qcol], each = na) *
        exp(Zc * rep(-(g$mo / 12), each = na))
    }
    aoh <- g$ae_t %*% ah
    ch <- if (g$type == 0) aoh else g$alk_t %*% ah
    ch <- ch / rep(.rtmb_posfun(RTMB::colSums(ch), 0.001)$value, each = g$nb)
    age_hat[g$rows, seq_len(na)] <- t(ah)
    age_obs_hat[g$rows, seq_len(na)] <- t(aoh)
    comp_hat[g$rows, seq_len(g$nb)] <- t(ch)
    if (any(g$fit)) {
      # Multinomial with the offset on both proportions; counts are fractional
      hp <- ch[, g$fit, drop = FALSE] + pd$comp_offset
      prob <- hp / rep(RTMB::colSums(hp) + 1e-300, each = g$nb)
      ll <- g$ll_const + RTMB::colSums(g$obs_n * log(prob))
      comp_nll[f] <- comp_nll[f] - sum(p[["comp_weights"]][f] * ll)
      comp_nll_unweighted[f] <- comp_nll_unweighted[f] - sum(ll)
    }
  }

  # -- 12.1 (cpp:3245-3270): depletion against the final-year static reference; under
  # predation with no harvest rule, against the final (unfished) projection year
  if (is.null(msm)) {
    biomass_depletion <- biomass / B0[, Y]
    ssb_depletion <- ssb / SB0[, Y]
  } else {
    biomass_depletion <- biomass / biomass[, Y]
    ssb_depletion <- ssb / ssb[, Y]
  }

  # -- 13 (cpp:3279-4658): negative log-likelihood components. Rows are addressed by the
  # JnllRow labels in .JNLL_ROW_AXIS; columns count fleets or species by row.
  jrow <- stats::setNames(seq_along(.JNLL_ROW_AXIS), names(.JNLL_ROW_AXIS))
  jnll_comp <- RTMB::matrix(0, length(jrow), pd$n_col)
  for (f in seq_len(n_flt)) {
    # 13.1: lognormal index, mean bias-corrected by -sd^2/2
    ix <- pd$index[[f]]
    if (!is.null(ix) && any(ix$fit)) {
      mu <- log(index_hat[ix$rows[ix$fit]]) - ix$bias
      jnll_comp[jrow[["Index data"]], f] <- -sum(RTMB::dnorm(ix$obs, mu, ix$sd[ix$fit], log = TRUE))
    }
    # 13.2: lognormal catch
    cx <- pd$catch[[f]]
    if (!is.null(cx) && any(cx$fit)) {
      mu <- log(catch_hat[cx$rows[cx$fit]]) - cx$bias
      jnll_comp[jrow[["Catch data"]], f] <- -sum(RTMB::dnorm(cx$obs, mu, cx$sd[cx$fit], log = TRUE))
    }
  }
  jnll_comp[jrow[["Composition data"]], seq_len(n_flt)] <- comp_nll
  # 13.7 (cpp:3892-3960): NonParametric shape penalties on the first year's log curve
  # (Sel_curve_pen1 on any decline between adjacent bins, Sel_curve_pen2 on curvature)
  # and on its level, counted once per shared selectivity block
  for (f in which(pd$flt_sel_type == 2 & pd$flt_type > 0 & pd$flt_sel_lead == 1)) {
    lnp <- S$log_np[[f]]
    nb <- length(lnp)
    r <- lnp[-nb] - lnp[-1]
    h <- (abs(r) + r) / 2.0
    d2 <- diff(diff(lnp))
    jnll_comp[jrow[["Non-parametric selectivity"]], f] <-
      sum(p[["sel_curve_pen"]][f, 1] * (h * h)) + sum(p[["sel_curve_pen"]][f, 2] * d2 * d2) +
      2.0 * (S$avg_sel[[f]] * S$avg_sel[[f]])
  }
  # 13.9 (cpp:4253-4297): initial-age and recruitment deviations, N(-sigma_R^2/2, sigma_R)
  # under bias_adjust_proc; 13.11 (cpp:4384-4390): the zero-N floor penalty
  for (sp in seq_len(nspp)) {
    mu <- -pd$bias_adjust_proc * (R_sd[sp] * R_sd[sp]) / 2.0
    jnll_comp[jrow[["Initial abundance deviates"]], sp] <-
      -sum(RTMB::dnorm(p[["init_dev"]][sp, seq_len(nages[sp] - 1)], mu, R_sd[sp], log = TRUE))
    jnll_comp[jrow[["Recruitment deviates"]], sp] <-
      -sum(RTMB::dnorm(p[["rec_dev"]][sp, hind], mu, R_sd[sp], log = TRUE))
    jnll_comp[jrow[["Zero n-at-age penalty"]], sp] <- zero_N_pen[sp]
  }
  # 11 (cpp:3233-3237) and 13.6 (cpp:3753-3817): predicted stomach content from the last
  # pass, and its likelihood, where suitability is estimated
  diet_hat <- stomach <- NULL
  if (!is.null(msm) && pd$suitMode[1] > 0) {
    diet_hat <- .rtmb_diet_hat(msm, su$S, avgN_c, pred$inv_avail)
    stomach <- .rtmb_stomach_nll(msm, diet_hat, p[["diet_comp_weights"]])
    jnll_comp[jrow[["Stomach content data"]], seq_len(nspp)] <- stomach$weighted
  }
  # 13.15 (cpp:4646-4658): unweighted copy; the composition and stomach rows differ
  unweighted_jnll_comp <- jnll_comp
  unweighted_jnll_comp[jrow[["Composition data"]], seq_len(n_flt)] <- comp_nll_unweighted
  if (!is.null(stomach)) {
    unweighted_jnll_comp[jrow[["Stomach content data"]], seq_len(nspp)] <- stomach$unweighted
  }
  jnll <- sum(jnll_comp)

  # -- 14 (cpp:4661-4887): reports, at the C++ names and shapes
  N_at_age <- to_array(N, nspp, nages)
  avgN_at_age <- to_array(avgN, nspp, nages)
  M1_at_age <- to_array(M1, nspp, nages)
  M_at_age <- to_array(M, nspp, nages)
  Z_at_age <- to_array(Z, nspp, nages)
  F_at_age <- to_array(F_at_age, nspp, nages)
  NByage0 <- to_array(NB0, nspp, nages)
  NByageF <- to_array(NBF, nspp, nages)
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
    vul <- if (is.null(su$vulnerability)) matrix(0, nspp, nspp) else su$vulnerability
    vul_other <- if (is.null(su$vulnerability_other)) numeric(nspp) else su$vulnerability_other
    .rtmb_report_list(.rtmb_predation_reports(msm, su$S, su$suit_other, vul, vul_other, pred,
                                              avgN_c, diet_hat))
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
