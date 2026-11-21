# One-time data preparation for the RTMB model, in plain R and off the AD tape. It reads
# the same `data` list the C++ template receives (TMB stores every element as a double,
# with 0-based fleet indices where rearrange_data() made them so) and precomputes every
# quantity that depends on data alone.


#' Refuse a configuration the RTMB model does not implement yet
#'
#' @details The RTMB port reproduces the C++ only on the paths ported so far: the
#'   single-species core BS2017SS exercises (stage 1), and the predation BS2017MS
#'   exercises (stage 2). Any other configuration would return a plausible but wrong
#'   number, so it stops here and names the setting.
#' @noRd
.rtmb_check_scope <- function(d) {
  need <- function(ok, what) {
    if (!isTRUE(all(ok))) {
      stop("The RTMB model does not yet support ", what, ".", call. = FALSE)
    }
  }
  # Every switch read below must exist: need() passes on an empty input
  switches <- c("estimateMode", "msmMode", "suitMode", "nsex", "nages", "estDynamics",
                "growth_model", "mat_len_use", "initMode", "srr_fun", "srr_pred_fun",
                "proj_mean_rec", "HCR", "DynamicHCR", "M1_re", "M1_use_prior", "M2_use_prior",
                "osa_mode", "adreport_sel", "adreport_catch", "flt_type", "flt_sel_type",
                "flt_sel_dim", "flt_varying_sel", "bin_first_selected", "sel_norm_bin1",
                "flt_units", "est_index_q", "index_varying_q", "comp_offset", "niter",
                "diet_ll_type", "suit_styr", "suit_endyr", "MSSB0", "MSB0", "other_food")
  missing <- switches[!vapply(switches, function(n) length(d[[n]]) > 0, TRUE)]
  if (length(missing)) stop("The RTMB model needs data element(s) ", toString(missing), ".", call. = FALSE)
  idx_flt <- unique(d[["index_ctl"]][, 1])
  cat_flt <- unique(d[["catch_ctl"]][, 1])
  comp_flt <- unique(d[["comp_ctl"]][, 1])
  need(d[["estimateMode"]] < 4, "estimateMode 4 (a placeholder objective)")
  need(d[["msmMode"]] %in% 0:2, "msmMode other than 0, 1 or 2")
  need(d[["suitMode"]] %in% c(0, 4), "suitMode other than 0 (empirical) or 4 (lognormal on weight)")
  need(length(unique(d[["suitMode"]])) == 1, "suitMode differing between predators")
  # The C++ still scores stomachs for an estimated suitability without predation
  need(d[["msmMode"]] > 0 | d[["suitMode"]] == 0, "estimated suitability at msmMode 0")
  if (isTRUE(d[["msmMode"]] > 0)) {
    dc <- d[["diet_ctl"]]
    need(d[["niter"]] >= 1 & d[["niter"]] == round(d[["niter"]]), "niter other than a whole number >= 1")
    need(d[["minage"]] == 1, "a minimum age other than 1 under predation")
    # data_check() refuses 0; available food is positive only for positive other food
    need(d[["other_food"]] > 0, "other food of 0 or less")
    need(dc[, 7] == 0, "diet rows for a single year (Year other than 0)")
    need(dc[, 5] >= d[["minage"]][dc[, 1]] & dc[, 6] >= d[["minage"]][dc[, 2]],
         "diet-data types 2-4 (rows summed over ages)")
    need(d[["diet_ll_type"]][d[["suitMode"]] > 0] == 0, "a Dirichlet-multinomial diet likelihood")
  }
  need(d[["nsex"]] == 1, "two-sex species")
  need(d[["nages"]] >= 2, "a species with one age")
  need(d[["estDynamics"]] == 0, "estDynamics > 0 (fixed numbers-at-age)")
  need(d[["growth_model"]] == 0, "estimated growth")
  need(d[["mat_len_use"]] == 0, "maturity-at-length")
  need(d[["initMode"]] == 2, "initMode other than 2 (NonEquilibrium)")
  need(c(d[["srr_fun"]], d[["srr_pred_fun"]]) == 0, "a stock-recruit curve")
  need(d[["proj_mean_rec"]] == 1, "proj_mean_rec = 0")
  need(c(d[["HCR"]], d[["DynamicHCR"]]) == 0, "a harvest control rule")
  need(c(d[["M1_re"]], d[["M1_use_prior"]], d[["M2_use_prior"]]) == 0,
       "M1 random effects or M priors")
  need(length(d[["linkage_process"]]) == 0, "linkages")
  need(nrow(d[["caal_ctl"]]) == 0, "CAAL data")
  need(nrow(d[["equil_catch_ctl"]]) == 0, "initial equilibrium catch")
  need(c(d[["osa_mode"]], d[["adreport_sel"]], d[["adreport_catch"]]) == 0,
       "OSA builds or selectivity/catch standard errors")
  need(d[["flt_type"]] %in% c(1, 2), "fleets switched off (Fleet_type 0)")
  need(d[["flt_sel_type"]] %in% c(0, 1, 2),
       "selectivity forms other than Fixed, Logistic and NonParametric")
  need(d[["flt_sel_dim"]] == 0, "length-based selectivity")
  need(d[["flt_varying_sel"]] == 0, "time-varying selectivity")
  need(d[["bin_first_selected"]] == 0, "Bin_first_selected above the first bin")
  need(d[["sel_norm_bin1"]][d[["flt_sel_type"]] != 0] == -999, "selectivity normalization")
  need(d[["flt_units"]] == 1, "observations in numbers (Observation_units 2)")
  need(d[["flt_type"]][idx_flt] == 2, "a fishery-type index")
  need(d[["index_ll_type"]][idx_flt] == 0, "index likelihoods other than lognormal")
  need(d[["est_sigma_index"]][idx_flt] == 0, "an estimated index sd")
  # Every fleet, not only those with an index: the q prior and deviations are scored on any
  # fleet leading its catchability block, and survey compositions read q
  need(d[["est_index_q"]] %in% c(NA, 0, 1), "catchability forms other than Fixed or Estimated")
  need(is.na(d[["index_varying_q"]]) | d[["index_varying_q"]] == 0, "time-varying catchability")
  need(d[["est_sigma_fsh"]][cat_flt] == 0, "an estimated catch sd")
  need(d[["comp_ll_type"]][comp_flt] == 0, "composition likelihoods other than multinomial")
  need(d[["comp_offset"]] > 0, "comp_offset = 0")
  need(d[["comp_ctl"]][, 3] == 0, "sex-specific or joint compositions")
  need(d[["comp_ctl"]][, 4] %in% c(0, 1), "composition types other than age or length")
  # The C++ takes a row's species from its Species column; the port takes the fleet's
  row_sp_ok <- function(ctl) ctl[, 2] == d[["flt_spp"]][ctl[, 1]] + 1
  need(row_sp_ok(d[["index_ctl"]]), "an index row whose Species is not its fleet's")
  need(row_sp_ok(d[["catch_ctl"]]), "a catch row whose Species is not its fleet's")
  need(row_sp_ok(d[["comp_ctl"]]), "a composition row whose Species is not its fleet's")
  es <- d[["emp_sel_ctl"]][d[["flt_sel_type"]][d[["emp_sel_ctl"]][, 1]] == 0, , drop = FALSE]
  need(row_sp_ok(es), "an emp_sel row whose Species is not its fleet's")
  need(es[, 3] %in% c(0, 1), "emp_sel Sex 2 on a one-sex species")
  need(es[, 4] >= d[["styr"]], "emp_sel rows before the start year")
  comp_sp <- d[["flt_spp"]][d[["comp_ctl"]][, 1]] + 1
  nb <- ifelse(d[["comp_ctl"]][, 4] == 0, d[["nages"]][comp_sp], d[["nlengths"]][comp_sp])
  yng <- d[["comp_accum_young"]][d[["comp_ctl"]][, 1]]
  old <- d[["comp_accum_old"]][d[["comp_ctl"]][, 1]]
  need(yng <= 1 & (old < 1 | old >= nb), "composition tail accumulation")
  invisible(TRUE)
}

#' Data-only quantities and row groupings for the RTMB objective
#'
#' @param data The data list the C++ `MakeADFun()` receives (`fit$obj$env$data`).
#' @return A list read by `.rtmb_objective()`. Population matrices are per species,
#'   [age, year]; observation rows are grouped by fleet.
#' @noRd
.rtmb_prepare_data <- function(data) {
  d <- data
  .rtmb_check_scope(d)
  int <- function(x) as.integer(round(x))
  nspp <- int(d[["nspp"]])
  nages <- int(d[["nages"]])
  max_age <- max(nages)
  styr <- int(d[["styr"]])
  endyr <- int(d[["endyr"]])
  nyrs <- int(d[["projyr"]]) - styr + 1L
  nyrs_hind <- endyr - styr + 1L
  flt_spp <- int(d[["flt_spp"]]) + 1L
  n_flt <- length(flt_spp)
  yr_ind <- pmin(seq_len(nyrs), nyrs_hind)   # projection years read the last hindcast year
  # A data year is calendar; a negative year (a held-out row) counts from its magnitude.
  col_of <- function(y) int(ifelse(y > 0, y - styr, ifelse(y < 0, -y - styr, 0))) + 1L
  # A comparison with an NA observation is false in the C++ (NaN > 0), not NA.
  true_only <- function(x) !is.na(x) & x

  # 5.7 (growth.hpp, empirical weight): weight-at-age by species role and by fleet (kg)
  wo <- d[["weight_obs"]]
  w_of <- function(row, sp) wo[row, 1, seq_len(nages[sp]), yr_ind, drop = TRUE]
  W_pop <- lapply(seq_len(nspp), function(sp) w_of(int(d[["pop_wt_index"]][sp]), sp))
  W_ssb <- lapply(seq_len(nspp), function(sp) w_of(int(d[["ssb_wt_index"]][sp]), sp))
  W_flt <- lapply(seq_len(n_flt), function(f) w_of(int(d[["flt_wt_index"]][f]) + 1L, flt_spp[f]))
  weight_hat <- array(0, c(nspp * 2 + n_flt, 1, max_age, nyrs))
  for (sp in seq_len(nspp)) {
    weight_hat[2 * sp - 1, 1, seq_len(nages[sp]), ] <- W_pop[[sp]]
    weight_hat[2 * sp, 1, seq_len(nages[sp]), ] <- W_ssb[[sp]]
  }
  for (f in seq_len(n_flt)) weight_hat[2 * nspp + f, 1, seq_len(nages[flt_spp[f]]), ] <- W_flt[[f]]

  # 5.4 and 5.7: spawning output, female SSB weight times mature females (one-sex
  # species: maturity times the data sex ratio, before 6.4 resets the first age's ratio to 1)
  spawn_output <- array(0, c(nspp, max_age, nyrs))
  so <- vector("list", nspp)
  for (sp in seq_len(nspp)) {
    a <- seq_len(nages[sp])
    mature_females <- d[["maturity"]][sp, a] * d[["sex_ratio"]][sp, a]
    so[[sp]] <- W_ssb[[sp]] * mature_females
    spawn_output[sp, a, ] <- so[[sp]]
  }

  # 5.9 (bioenergetics.hpp): temperature function and annual ration (kg/yr); data only
  env <- d[["env_index"]]
  fT <- matrix(0, nspp, nyrs)
  consumption_at_age <- array(0, c(nspp, 1, max_age, nyrs))
  for (sp in seq_len(nspp)) {
    ci <- int(d[["Cindex"]][sp])
    eq <- if (ci >= 1 && ci <= ncol(env)) int(d[["Ceq"]][sp]) else 4L
    temp <- if (ci >= 1 && ci <= ncol(env)) env[, ci] else rep(0, nyrs)
    Qc <- d[["Qc"]][sp]; Tcm <- d[["Tcm"]][sp]; Tco <- d[["Tco"]][sp]
    fT[sp, ] <- switch(as.character(eq),
      "1" = exp(Qc * temp),
      "2" = {
        Yc <- log(Qc) * (Tcm - Tco + 2.0)
        Zc <- log(Qc) * (Tcm - Tco)
        Vc <- (Tcm - temp) / (Tcm - Tco)
        Xc <- Zc^2 * (1.0 + (1.0 + 40.0 / Yc)^0.5)^2 / 400.0
        Vc^Xc * exp(Xc * (1.0 - Vc))
      },
      "3" = {
        Tcl <- d[["Tcl"]][sp]; CK1 <- d[["CK1"]][sp]; CK4 <- d[["CK4"]][sp]
        G2 <- (1.0 / (Tcl - Tcm)) * log((0.98 * (1.0 - CK4)) / (CK4 * 0.02))
        L2 <- exp(G2 * (Tcl - temp))
        Kb <- (CK4 * L2) / (1.0 + CK4 * (L2 - 1.0))
        G1 <- (1.0 / (Tco - Qc)) * log((0.98 * (1.0 - CK1)) / (CK1 * 0.02))
        L1 <- exp(G1 * (temp - Qc))
        Ka <- (CK1 * L1) / (1.0 + CK1 * (L1 - 1.0))
        Ka * Kb
      },
      rep(1.0, nyrs))
    a <- seq_len(nages[sp])
    ration <- d[["ration_data"]][sp, 1, a, yr_ind, drop = TRUE]
    cons <- 0   # an equation code above 4 leaves consumption at 0, as in the C++
    if (d[["Ceq"]][sp] < 4) {
      cons <- d[["CA"]][sp] * (W_pop[[sp]] * 1000.0)^(1.0 + d[["CB"]][sp]) *
        rep(fT[sp, ], each = nages[sp]) * d[["fday"]][sp]
      cons <- cons * d[["Pvalue"]][sp] * ration / 1000.0
    }
    if (d[["Ceq"]][sp] == 4) cons <- d[["Pvalue"]][sp] * ration
    consumption_at_age[sp, 1, a, ] <- cons
  }

  # 5.8 (selectivity.hpp), form 0: selectivity read from emp_sel; zero in hindcast
  # years without a row, and the last year with a row held through the projection
  sel_fixed <- vector("list", n_flt)
  for (f in which(d[["flt_sel_type"]] == 0)) {
    na <- nages[flt_spp[f]]
    s <- matrix(0, na, nyrs)
    rows <- which(d[["emp_sel_ctl"]][, 1] == f)
    yrs <- int(d[["emp_sel_ctl"]][rows, 4]) - styr + 1L
    rows <- rows[yrs <= nyrs_hind]; yrs <- yrs[yrs <= nyrs_hind]
    for (i in seq_along(rows)) {
      v <- d[["emp_sel_obs"]][rows[i], seq_len(na)]
      skip <- is.na(v) & !is.nan(v)   # the C++ skips NA only (R_IsNA), not NaN
      s[!skip, yrs[i]] <- v[!skip]
    }
    if (length(yrs) > 0) s[, seq_len(nyrs - nyrs_hind) + nyrs_hind] <- s[, max(yrs)]
    sel_fixed[[f]] <- s
  }

  # 8.1: index rows by fleet; the sd is the row's log-scale sd and the observation the
  # log index (obsvec); a row is fitted in 0 < year <= endyr with a positive index
  ic <- d[["index_ctl"]]
  index <- lapply(seq_len(n_flt), function(f) {
    r <- which(ic[, 1] == f)
    if (length(r) == 0) return(NULL)
    yr <- ic[r, 3]
    col <- col_of(yr)
    pos <- int(d[["index_obsvec_idx"]][r])
    fit <- true_only(yr > 0 & yr <= endyr & d[["index_obs"]][r, 1] > 0 & pos >= 0)
    sd <- d[["index_obs"]][r, 2]
    list(rows = r, col = col, qcol = pmin(col, nyrs_hind), mo = d[["index_n"]][r, 1],
         fit = fit, sd = sd, obs = d[["obsvec"]][pos[fit] + 1L],
         bias = d[["bias_adjust_obs"]] * (sd[fit] * sd[fit]) / 2.0)
  })

  # 9.1: catch rows by fleet; fitted in 0 < year <= endyr with a positive catch
  cc <- d[["catch_ctl"]]
  catch <- lapply(seq_len(n_flt), function(f) {
    r <- which(cc[, 1] == f)
    if (length(r) == 0) return(NULL)
    yr <- cc[r, 3]
    pos <- int(d[["catch_obsvec_idx"]][r])
    pos_obs <- d[["catch_obs"]][r, 1] > 0
    fit <- true_only(yr > 0 & yr <= endyr & d[["flt_type"]][f] == 1 & pos_obs & pos >= 0)
    analytical <- true_only(yr > 0 & yr <= endyr & d[["flt_type"]][f] == 1 & pos_obs)
    sd <- d[["catch_obs"]][r, 2]
    list(rows = r, col = col_of(yr), fit = fit, sd = sd,
         obs = d[["obsvec"]][pos[fit] + 1L],
         bias = d[["bias_adjust_obs"]] * (sd[fit] * sd[fit]) / 2.0,
         analytical = analytical, log_obs_analytical = log(d[["catch_obs"]][r[analytical], 1]))
  })

  # 10.1 and 13.4: composition rows grouped by fleet and type (0 age, 1 length); the
  # multinomial's data-only term lgamma(sum(x) + 1) - sum(lgamma(x + 1)) is precomputed
  cmp <- d[["comp_ctl"]]
  comp_offset <- d[["comp_offset"]]
  keys <- unique(cmp[, c(1, 4), drop = FALSE])
  comp <- lapply(seq_len(nrow(keys)), function(k) {
    f <- int(keys[k, 1]); type <- int(keys[k, 2]); sp <- flt_spp[f]
    r <- which(cmp[, 1] == f & cmp[, 4] == type)
    yr <- cmp[r, 5]
    na <- nages[sp]
    nb <- if (type == 0) na else int(d[["nlengths"]][sp])
    col <- col_of(yr)
    fit <- true_only(yr <= endyr & yr > 0 & d[["flt_type"]][f] > 0 & d[["comp_n"]][r, 2] > 0)
    obs_n <- t((d[["comp_obs"]][r[fit], seq_len(nb), drop = FALSE] + comp_offset) *
               d[["comp_n"]][r[fit], 2])
    ae <- d[["age_error"]][int(d[["flt_ageing_error_index"]][f]) + 1L, seq_len(na), seq_len(na)]
    alk <- if (type == 1) {
      d[["age_trans_matrix"]][int(d[["flt_age_transition_index"]][f]) + 1L, 1,
                              seq_len(na), seq_len(nb)]
    }
    list(flt = f, type = type, sp = sp, rows = r, col = col, qcol = pmin(col, nyrs_hind),
         mo = d[["comp_n"]][r, 1], nb = nb, fit = fit, obs_n = obs_n,
         ll_const = lgamma(colSums(obs_n) + 1) - colSums(lgamma(obs_n + 1)),
         ae_t = t(matrix(ae, na, na)), alk_t = if (type == 1) t(matrix(alk, na, nb)))
  })

  # 14: reports that are data, or zero on every path ported so far (predation, growth,
  # linkages, CAAL), at the C++ names and shapes so rename_output() reads them unchanged
  max_nlengths <- max(int(d[["nlengths"]]))
  max_nlengths_pop <- max(int(d[["nlengths_pop"]]))
  pop4 <- c(nspp, 1, max_age, nyrs)
  pair5 <- c(nspp, nspp, max_age, max_age, nyrs)
  fixed_reports <- list(
    weight_hat = weight_hat, spawn_output = spawn_output, fT = fT,
    consumption_at_age = consumption_at_age,
    comp_obs = d[["comp_obs"]], caal_obs = d[["caal_obs"]], diet_obs = d[["diet_obs"]],
    caal_hat = matrix(0, nrow(d[["caal_obs"]]), ncol(d[["caal_obs"]])),
    diet_hat = matrix(0, nrow(d[["diet_obs"]]), ncol(d[["diet_obs"]])), equil_catch_hat = numeric(0),
    catch_sd = d[["catch_obs"]][, 2], index_sd = d[["index_obs"]][, 2],
    proj_F = matrix(0, nspp, nyrs), rec_srr_single_density =
      as.numeric(!(d[["srr_fun"]] == 0 && d[["srr_pred_fun"]] > 0)),
    pred_CAAL = array(0, c(n_flt, 1, max_age, max_nlengths, nyrs)),
    growth_matrix = array(0, c(nspp * 2 + n_flt, 1, max_age, max_nlengths_pop, nyrs)),
    length_hat = array(0, c(nspp * 2 + n_flt, 1, max_age, nyrs)),
    sel_at_length = array(0, c(n_flt, 1, max_nlengths_pop, nyrs)),
    M2_at_age = array(0, pop4), B_eaten_as_prey = array(0, pop4),
    suit_other = array(0, pop4), avail_food = array(0, pop4),
    suitability = array(0, pair5), M2_prop = array(0, pair5), B_eaten = array(0, pair5),
    stom_div_bio = array(0, pair5),
    vulnerability = matrix(0, nspp, nspp), vulnerability_other = numeric(nspp),
    q_linkage_offset = matrix(0, n_flt, nyrs_hind),
    q_linkage_offset_nat = matrix(0, n_flt, nyrs_hind),
    q_linkage_log_mult = matrix(0, n_flt, nyrs_hind),
    recruitment_linkage_offset = array(0, c(nspp, 4, nyrs)),
    recruitment_linkage_offset_nat = array(0, c(nspp, 4, nyrs)),
    M_linkage_offset = array(0, pop4), M_linkage_offset_nat = array(0, pop4),
    growth_linkage_offset = array(0, c(nspp, 1, nyrs, 4)),
    growth_linkage_offset_nat = array(0, c(nspp, 1, nyrs, 4)))

  # Under predation (msmMode > 0) the predation reports are computed, and the C++ runs
  # niter passes of the dynamics; at msmMode 0 it runs one (cpp:189)
  msmMode <- int(d[["msmMode"]])
  msm <- NULL
  if (msmMode > 0) {
    msm <- .rtmb_prepare_predation(d, nages, nyrs, nyrs_hind, styr, W_pop, consumption_at_age)
    fixed_reports[c("M2_at_age", "B_eaten_as_prey", "suit_other", "avail_food", "suitability",
                    "M2_prop", "B_eaten", "stom_div_bio", "vulnerability",
                    "vulnerability_other", "diet_hat")] <- NULL
  }

  list(
    msmMode = msmMode, niter = if (msmMode > 0) int(d[["niter"]]) else 1L,
    suitMode = int(d[["suitMode"]]), msm = msm,
    nspp = nspp, nages = nages, max_age = max_age, n_flt = n_flt, n_col = max(n_flt, nspp),
    styr = styr, endyr = endyr, nyrs = nyrs, nyrs_hind = nyrs_hind,
    nyrs_srrmean = min(int(d[["srr_mse_switchyr"]]) - styr + 1L, nyrs_hind),
    flt_spp = flt_spp, flt_type = int(d[["flt_type"]]), flt_sel_type = int(d[["flt_sel_type"]]),
    flt_n_sel_bins = int(d[["flt_n_sel_bins"]]), flt_sel_lead = int(d[["flt_sel_lead"]]),
    spawn_month = d[["spawn_month"]], bias_adjust_proc = d[["bias_adjust_proc"]],
    bias_adjust_obs = d[["bias_adjust_obs"]], comp_offset = comp_offset, env = env,
    W_pop = W_pop, W_flt = W_flt, so = so, sel_fixed = sel_fixed,
    index = index, catch = catch, comp = comp,
    fixed_reports = fixed_reports, comp_obs = d[["comp_obs"]],
    n_catch = nrow(cc), n_index = nrow(ic), n_comp = nrow(cmp),
    max_age_cols = max(ncol(d[["comp_obs"]]), nages[int(cmp[, 2])])
  )
}

#' Data-only pieces of the predation model, on one-sex (species, age) cells
#'
#' @details Cells run species, then age: the C++ cell order (predation.hpp:387-401) and
#'   the rows of `do.call(rbind, <per-species [age, year] list>)`. Weights and ration are
#'   data, because estimated growth is refused. The scope guard allows only Year 0 diet
#'   rows for one predator and one prey age, so the diet proportions are year-invariant.
#' @param W_pop Population weight-at-age (kg), per species [age, year].
#' @param consumption_at_age Ration (kg/predator/yr), [nspp, 1, max_age, nyrs].
#' @noRd
.rtmb_prepare_predation <- function(d, nages, nyrs, nyrs_hind, styr, W_pop, consumption_at_age) {
  int <- function(x) as.integer(round(x))
  nspp <- length(nages)
  max_age <- max(nages)
  cell_sp <- rep(seq_len(nspp), nages)
  cell_age <- sequence(nages)
  nc <- length(cell_sp)
  first <- cumsum(c(0L, nages))[seq_len(nspp)]
  suitMode <- int(d[["suitMode"]])
  W <- do.call(rbind, W_pop)
  # Positive available food under empirical suitability rests on non-negative prey biomass
  if (any(W < 0)) {
    stop("The RTMB model does not yet support negative weight-at-age under predation.", call. = FALSE)
  }
  Rat <- do.call(rbind, lapply(seq_len(nspp), function(sp) {
    matrix(consumption_at_age[sp, 1, seq_len(nages[sp]), ], nages[sp], nyrs)
  }))
  # Suitability window per predator, as year columns (cpp:139-143: calendar years in data)
  suit_cols <- lapply(seq_len(nspp), function(sp) {
    seq(int(d[["suit_styr"]][sp]) - styr + 1L, int(d[["suit_endyr"]][sp]) - styr + 1L)
  })

  # Diet proportion by (predator cell, prey cell) (diet_data.hpp:24-87); a repeated row
  # keeps its last value, as the C++ assignment does
  dc <- d[["diet_ctl"]]
  dob <- d[["diet_obs"]]
  minage <- int(d[["minage"]])
  dr <- first[int(dc[, 1])] + int(dc[, 5]) - minage[int(dc[, 1])] + 1L
  dk <- first[int(dc[, 2])] + int(dc[, 6]) - minage[int(dc[, 2])] + 1L
  dp <- matrix(0, nc, nc)
  for (i in seq_len(nrow(dc))) dp[dr[i], dk[i]] <- dob[i, 2]
  # Unattributed diet share per mt of other food (diet_data.hpp:97-133), subtracted in
  # prey-cell order
  other_food <- d[["other_food"]]
  ofdp <- vapply(seq_len(nc), function(r) {
    of <- other_food[cell_sp[r]]
    if (of > 0) Reduce(`-`, dp[r, ], 1) / of else 0
  }, 0)

  # Empirical suitability: the (predator, prey) pairs with a positive diet proportion,
  # ordered predator then prey, the order the C++ sums them in
  emp <- NULL
  emp_sp <- which(suitMode == 0)
  if (length(emp_sp)) {
    tr <- which(dp > 0 & (cell_sp %in% emp_sp)[row(dp)], arr.ind = TRUE)
    tr <- tr[order(tr[, 1], tr[, 2]), , drop = FALSE]
    cols <- sort(unique(unlist(suit_cols[emp_sp])))
    by_pred <- lapply(emp_sp, function(sp) {
      list(t = which(cell_sp[tr[, 1]] == sp), cols = match(suit_cols[[sp]], cols),
           nyrs = length(suit_cols[[sp]]))
    })
    emp <- list(r = tr[, 1], k = tr[, 2], d = dp[tr], cols = cols,
                by_pred = by_pred[lengths(lapply(by_pred, `[[`, "t")) > 0],
                G = Matrix::sparseMatrix(i = tr[, 1], j = seq_len(nrow(tr)), x = 1,
                                         dims = c(nc, nrow(tr))))
  }

  # Parametric suitability: the (predator, prey, hindcast year) cells whose log weight
  # ratio is positive (predation.hpp:115-123); data, because weight is
  ps <- NULL
  param_sp <- which(suitMode > 0)
  if (length(param_sp)) {
    g <- expand.grid(r = which(cell_sp %in% param_sp), k = seq_len(nc), y = seq_len(nyrs_hind))
    sr <- log(W[cbind(g$r, g$y)] / W[cbind(g$k, g$y)])
    keep <- !is.na(sr) & sr > 0
    g <- g[keep, ]
    ps <- list(sp_r = cell_sp[g$r], sp_k = cell_sp[g$k], sr = sr[keep],
               lin = g$r + nc * (g$k - 1L) + nc * nc * (g$y - 1L))
  }

  # Predicted stomach content (diet_data.hpp:160-235): each row's (predator cell, year)
  # and (prey cell, year) positions over its predator's window, [row, year] column-major
  dh <- st <- NULL
  if (length(param_sp)) {
    dh <- lapply(param_sp, function(sp) {
      rows <- which(dc[, 1] == sp)
      if (!length(rows)) return(NULL)
      yoff <- rep(suit_cols[[sp]] - 1L, each = length(rows))
      list(rows = rows, nyrs = length(suit_cols[[sp]]),
           ry = rep(dr[rows], length(suit_cols[[sp]])) + nc * yoff,
           ky = rep(dk[rows], length(suit_cols[[sp]])) + nc * yoff,
           rky = rep(dr[rows] + nc * (dk[rows] - 1L), length(suit_cols[[sp]])) + nc * nc * yoff)
    })
    dh <- Filter(Negate(is.null), dh)
    # Stomachs (cpp:3753-3817): runs of stomach_id, scored when their predator's
    # suitability is estimated. The observed side is data: proportions plus the other-prey
    # balance, offset and renormalized, times the sample size; and dmultinom's constant.
    sid_all <- int(d[["stomach_id"]])
    groups <- split(seq_len(nrow(dc)), factor(sid_all, levels = seq_len(int(d[["n_stomach_obs"]])) - 1L))
    groups <- groups[lengths(groups) > 0]
    gsp <- vapply(groups, function(j) int(dc[j[1], 1]), 1L)
    groups <- groups[gsp %in% param_sp]
    gsp <- gsp[gsp %in% param_sp]
    off <- d[["comp_offset"]]
    obs <- lapply(groups, function(j) {
      s <- min(Reduce(`+`, dob[j, 2]), 1)
      o <- c(dob[j, 2], 1 - s) + off
      x <- o / Reduce(`+`, o) * dob[j[1], 1]
      list(x = x, const = lgamma(Reduce(`+`, x) + 1) - Reduce(`+`, lgamma(x + 1)))
    })
    rows <- unlist(groups, use.names = FALSE)
    sid <- rep(seq_along(groups), lengths(groups))
    st <- list(rows = rows, sid = sid, sp = gsp,
               x = unlist(lapply(obs, function(o) o$x[-length(o$x)]), use.names = FALSE),
               x_other = vapply(obs, function(o) o$x[length(o$x)], 0),
               const = vapply(obs, `[[`, 0, "const"),
               G = Matrix::sparseMatrix(i = sid, j = seq_along(rows), x = 1,
                                        dims = c(length(groups), length(rows))))
  }

  # Positions of each cell pair in the C++ [nspp, nspp, max_age, max_age, nyrs] reports,
  # and of each cell in [nspp, 1, max_age, nyrs]
  pair_lin <- outer(cell_sp + nspp * nspp * (cell_age - 1L),
                    nspp * (cell_sp - 1L) + nspp * nspp * max_age * (cell_age - 1L), "+")
  lin4 <- outer(cell_sp + nspp * (cell_age - 1L), nspp * max_age * (seq_len(nyrs) - 1L), "+")
  list(nspp = nspp, nc = nc, nyrs = nyrs, nyrs_hind = nyrs_hind, msmMode = int(d[["msmMode"]]),
       cell_sp = cell_sp, cell_rows = split(seq_len(nc), cell_sp), W = W, Rat = Rat,
       other_food_c = other_food[cell_sp], ofdp = ofdp, emp = emp, param_sp = param_sp,
       ps = ps, dh = dh, st = st, n_diet = nrow(dc), comp_offset = d[["comp_offset"]],
       pair_lin = pair_lin, n_pair4 = nspp * nspp * max_age * max_age,
       dim5 = c(nspp, nspp, max_age, max_age, nyrs), dim4 = c(nspp, 1, max_age, nyrs),
       lin4 = lin4, MSSB0 = d[["MSSB0"]], MSB0 = d[["MSB0"]])
}
