# Multispecies predation for the RTMB model: Type II and Type III MSVPA (Holsman et al.
# 2015). C++ pointers are to this tree's src/TMB (5.67.11): predation.hpp, diet_data.hpp.
#
# Predators and prey are the same one-sex (species, age) "cells", stacked species then age
# into rows, so a quantity by cell and year is a [cell, year] matrix, and suitability is
# [predator cell, prey cell], or [predator cell, prey cell, year] where it varies by year.
# Numbers are thousands of fish and weight kg, so numbers x weight is mt; ration is kg per
# predator per year. Pairs with no interaction hold 0. The AD mechanics (0/1 indicators,
# `[<-` overloads, the summation helpers) are in 6-rtmb_helpers.R.


#' Size-based suitability, lognormal on the predator-to-prey weight ratio (suitMode 4)
#'
#' @details `suitability = vulnerability[predator, prey] x` a normal curve on the log
#'   predator-to-prey weight ratio, mean `gam_a` and sd `gam_b`, scaled to peak at 1; 0 where
#'   the ratio is not above 1. Vulnerability is a multinomial logit of `log_phi` whose
#'   remainder is other food (Trijoulet et al. 2020). Weights vary by hindcast year, so
#'   suitability does too; projection years reuse the last (predation.hpp:20-45, 76-142).
#' @return List: `suitability` `[cell, cell, hindcast year]`, `suit_other` `[cell]`,
#'   `vulnerability` `[nspp, nspp]`, `vulnerability_other` `[nspp]`.
#' @noRd
.rtmb_suit_parametric <- function(msm, log_phi, gam_a, gam_b) {
  "[<-" <- RTMB::ADoverload("[<-")
  nspp <- msm$nspp
  vulnerability <- matrix(0, nspp, nspp)
  vulnerability_other <- numeric(nspp)
  for (rsp in msm$param_sp) {
    phi <- exp(log_phi[rsp, ])
    vulnerability[rsp, ] <- phi / (1 + sum(phi))
    vulnerability_other[rsp] <- 1 - sum(vulnerability[rsp, ])
  }
  # One element per (predator cell, prey cell, hindcast year) with a positive log ratio
  pairs <- msm$size_pairs
  suit <- vulnerability[cbind(pairs$pred_sp, pairs$prey_sp)] *
    RTMB::dnorm(pairs$log_wt_ratio, gam_a[pairs$pred_sp], gam_b[pairs$pred_sp]) /
    RTMB::dnorm(gam_a, gam_a, gam_b)[pairs$pred_sp]
  suitability <- array(0, c(msm$ncell, msm$ncell, msm$nyrs_hind))
  suitability[pairs$index] <- suit
  list(suitability = suitability, suit_other = vulnerability_other[msm$cell_sp],
       vulnerability = vulnerability, vulnerability_other = vulnerability_other)
}

#' Empirical (MSVPA) suitability from observed diet and prey biomass (suitMode 0)
#'
#' @details `suitability[pred, prey] = mean over the predator's window years of
#'   (diet_prop / B_prey) / (sum over its prey of diet_prop / B_prey + other)`, where B_prey
#'   is N x W, prey biomass in mt (N^2 x W under Type III, not a biomass), and `other` the
#'   other-food diet share per mt, floored at 1e-6 of that sum. Held for every year
#'   (predation.hpp:179-308).
#' @param avgN Average numbers-at-age (thousands of fish), `[cell, year]`.
#' @return List: `suitability` `[cell, cell]`, `suit_other` `[cell]`.
#' @noRd
.rtmb_suit_empirical <- function(msm, avgN) {
  "[<-" <- RTMB::ADoverload("[<-")
  emp <- msm$empirical   # the (predator cell, prey cell) pairs with a positive diet proportion
  yrs <- emp$years       # year columns inside any predator's window
  prey_N_response <- avgN[, yrs, drop = FALSE]
  if (msm$msmMode == 2) prey_N_response <- prey_N_response * prey_N_response
  prey_denom <- msm$W[, yrs, drop = FALSE] * prey_N_response   # mt under Type II only
  # 1 / prey_denom where it is positive, else 0; NaN passes and reaches suitability
  has_prey <- 1 - .rtmb_ind_le(prey_denom, 0)   # (predation.hpp:227-229)
  inv_prey_denom <- has_prey / .rtmb_safe_den(prey_denom, has_prey)
  # Diet proportion / prey_denom, [pair, year], and its sum over each predator's prey
  stom_div_bio <- inv_prey_denom[emp$prey, , drop = FALSE] * emp$diet_prop
  sum_stom_div_bio <- .rtmb_group_sum(emp$sum_by_pred, stom_div_bio)
  # A predator age whose diet prey all have zero biomass gets suitability 0 in that year
  diet_prop_sum <- .rtmb_group_sum(emp$sum_by_pred,
                                   has_prey[emp$prey, , drop = FALSE] * emp$diet_prop)
  other_term <- .rtmb_max2(matrix(msm$other_stom_div_bio, msm$ncell, length(yrs)),
                           1e-6 * sum_stom_div_bio)
  inv_pred_denom <- .rtmb_ind_gt(diet_prop_sum, 0) /
    .rtmb_max2(sum_stom_div_bio + other_term, 1e-10)
  suit_by_year <- stom_div_bio * inv_pred_denom[emp$pred, , drop = FALSE]
  suitability <- matrix(0, msm$ncell, msm$ncell)
  for (w in emp$windows) {
    suitability[cbind(emp$pred[w$pairs], emp$prey[w$pairs])] <-
      RTMB::rowSums(suit_by_year[w$pairs, w$cols, drop = FALSE]) / w$nyrs
  }
  list(suitability = suitability, suit_other = 1 - RTMB::rowSums(suitability))
}

#' Available food and predation mortality M2 (per year), MSVPA (Holsman et al. 2015)
#'
#' @details `avail_food[predator] = sum over prey of suitability x N x W + other_food x
#'   suit_other`, in mt (prey N^2 x W under Type III, so not a biomass); then
#'   `M2[prey] = sum over predators of suitability x N x ration / avail_food`, unsquared.
#'   The dynamic B0 (F = 0) and BF (F = Ftarget) populations get the same with the
#'   fitted suitability. A predator age whose available food is not positive adds no M2;
#'   NaN available food passes that test and gives NaN (predation.hpp:349-528).
#' @param suitability `[cell, cell]`, or `[cell, cell, hindcast year]`.
#' @param avgN,avgN_dB0,avgN_dBF Average numbers-at-age (thousands of fish), `[cell, year]`:
#'   fished, dynamic B0 and dynamic BF.
#' @return List of `[cell, year]` matrices: `avail_food`, `inv_avail` (its reciprocal),
#'   `pred_mort_scalar` (N x ration / available food), `M2`, `M2_dB0`, `M2_dBF`.
#' @noRd
.rtmb_predation <- function(msm, suitability, suit_other, avgN, avgN_dB0, avgN_dBF) {
  nyrs <- msm$nyrs; nyrs_hind <- msm$nyrs_hind
  # For population k (fished, dB0, dBF): suit_sum(x, "prey")[[k]][pred, y] = sum over prey of
  # suitability[pred, prey, y] * x[[k]][prey, y]; "predators" sums over predators instead.
  suit_sum <- function(x, over = c("prey", "predators")) {
    over <- match.arg(over)
    orient <- function(S) if (over == "prey") S else t(S)
    X <- do.call(cbind, x)   # side by side, so one product covers all three: taping is faster
    if (length(dim(suitability)) == 2) {
      out <- .rtmb_matmul(orient(suitability), X)
    } else {
      cols <- function(yrs) c(yrs, nyrs + yrs, 2 * nyrs + yrs)   # yrs in each population
      # One product per hindcast year, plus one for all projection years (last hindcast's)
      blocks <- lapply(seq_len(nyrs_hind), cols)
      if (nyrs > nyrs_hind) blocks[[nyrs_hind + 1]] <- cols(seq(nyrs_hind + 1, nyrs))
      out <- do.call(cbind, lapply(seq_along(blocks), function(i) {
        .rtmb_matmul(orient(suitability[, , min(i, nyrs_hind)]),
                     X[, blocks[[i]], drop = FALSE])
      }))
      out <- out[, order(unlist(blocks)), drop = FALSE]   # one population after another
    }
    yr_cols <- function(k) out[, k * nyrs + seq_len(nyrs), drop = FALSE]
    list(fished = yr_cols(0), dB0 = yr_cols(1), dBF = yr_cols(2))
  }
  prey_N_response <- function(N) if (msm$msmMode == 2) N * N else N   # N, or N^2 (Type III)
  prey_food <- suit_sum(list(prey_N_response(avgN) * msm$W, prey_N_response(avgN_dB0) * msm$W,
                             prey_N_response(avgN_dBF) * msm$W))
  suitable_other_food <- rep(msm[["other_food"]] * suit_other, nyrs)
  avail_food <- prey_food[["fished"]] + suitable_other_food
  avail_food_dB0 <- prey_food[["dB0"]] + suitable_other_food
  avail_food_dBF <- prey_food[["dBF"]] + suitable_other_food
  # N (thousands) x ration (kg/predator/yr) is mt consumed, divided by available food
  if (is.null(msm$size_pairs)) {
    # suit_other >= 1e-6/(1 + 1e-6) and other food > 0, so available food is never <= 0; the
    # dB0, dBF lines write N x ration / avail_food in the C++'s order, to round alike
    inv_avail <- 1 / avail_food
    pred_mort_scalar <- avgN * msm$ration * inv_avail
    pred_mort_scalar_dB0 <- (avgN_dB0 * msm$ration) / avail_food_dB0
    pred_mort_scalar_dBF <- (avgN_dBF * msm$ration) / avail_food_dBF
  } else {
    # suit_other = 1 - sum(vulnerability) rounds to <= 0 at large log_phi, so food can be 0:
    # inv_avail = 1 / food where > 0, else 0; NaN passes and reaches M2 (predation.hpp:475-477)
    has_food <- 1 - .rtmb_ind_le(avail_food, 0)
    has_food_dB0 <- (1 - .rtmb_ind_le(avail_food_dB0, 0)) * has_food   # and fished food > 0
    has_food_dBF <- (1 - .rtmb_ind_le(avail_food_dBF, 0)) * has_food
    inv_avail <- has_food / .rtmb_safe_den(avail_food, has_food)
    pred_mort_scalar <- avgN * msm$ration * inv_avail
    pred_mort_scalar_dB0 <- has_food_dB0 *
      ((avgN_dB0 * msm$ration) / .rtmb_safe_den(avail_food_dB0, has_food_dB0))
    pred_mort_scalar_dBF <- has_food_dBF *
      ((avgN_dBF * msm$ration) / .rtmb_safe_den(avail_food_dBF, has_food_dBF))
  }
  M2 <- suit_sum(list(pred_mort_scalar, pred_mort_scalar_dB0, pred_mort_scalar_dBF),
                 "predators")
  list(avail_food = avail_food, inv_avail = inv_avail, pred_mort_scalar = pred_mort_scalar,
       M2 = M2[["fished"]], M2_dB0 = M2[["dB0"]], M2_dBF = M2[["dBF"]])
}

#' Predicted diet proportion of each diet row's prey age in its predator age
#'
#' @details `diet_prop_hat = N_prey x W_prey x suitability / avail_food`, the predicted share
#'   of this prey age in the diet (not a proportion under Type III, N being unsquared:
#'   predation.hpp:512), averaged unweighted over the predator's window years, 0 in a year
#'   without predators (diet_data.hpp:160-235, 286; diet-data type 1).
#' @return Vector over the diet rows (column 2 of the C++ `diet_hat`).
#' @noRd
.rtmb_diet_hat <- function(msm, suitability, avgN, inv_avail) {
  "[<-" <- RTMB::ADoverload("[<-")
  prey_B <- avgN * msm$W
  diet_hat <- numeric(msm$n_diet)
  for (g in msm$diet_rows) {   # one predator species' rows, over its window years
    diet_prop_hat <- prey_B[g$prey_pos] * suitability[g$pair_pos] * inv_avail[g$pred_pos]
    pred_N <- avgN[g$pred_pos]
    has_pred <- .rtmb_ind_gt(pred_N, 0)
    # (diet_prop_hat x N) / N is the C++'s numbers-weighted mean over sexes, here one sex
    term <- has_pred * (diet_prop_hat * pred_N) / .rtmb_safe_den(pred_N, has_pred)
    dim(term) <- c(length(g$rows), g$nyrs)
    diet_hat[g$rows] <- RTMB::rowSums(term) / g$nyrs
  }
  diet_hat
}

#' Stomach-content negative log-likelihood: multinomial on prey proportions plus other prey
#'
#' @details Per stomach sample, the predicted proportions of its prey and of other prey (1
#'   minus their sum, floored at 1e-5 by posfun, whose penalty does not enter the objective)
#'   each get `comp_offset` and are renormalized; the multinomial log-likelihood is weighted
#'   by the raw `diet_comp_weights` (`Diet_distribution` 0; cpp 13.6, 3753-3817).
#' @return List: `weighted` and `unweighted` negative log-likelihood by predator species.
#' @noRd
.rtmb_stomach_nll <- function(msm, diet_hat, diet_comp_weights) {
  "[<-" <- RTMB::ADoverload("[<-")
  stom <- msm$stomach
  off <- msm$comp_offset
  prey_hat <- diet_hat[stom$rows]   # each sample's prey rows, sample by sample
  other_hat <- .rtmb_posfun(1 - .rtmb_group_sum(stom$by_sample, prey_hat), 1e-5)$value
  total <- .rtmb_group_sum(stom$by_sample, prey_hat + off) + (other_hat + off)
  ll <- stom$ll_const +
    (.rtmb_group_sum(stom$by_sample, stom$obs_n * log((prey_hat + off) / total[stom$sample])) +
       stom$obs_n_other * log((other_hat + off) / total))
  weighted <- unweighted <- numeric(msm$nspp)
  for (sp in unique(stom$pred_sp)) {
    j <- stom$pred_sp == sp
    weighted[sp] <- -sum(diet_comp_weights[sp] * ll[j])
    unweighted[sp] <- -sum(ll[j])
  }
  list(weighted = weighted, unweighted = unweighted)
}

#' The predation REPORTs at the C++ names and shapes, from the last pass
#'
#' @details Plain R, run only in the `report()` pass. Each predator-prey product is formed
#'   as the C++ forms it (predation.hpp:496-526). At a NaN objective the C++ report can
#'   hold NaN in more cells, from 0 x NaN products this port never forms.
#' @noRd
.rtmb_predation_reports <- function(msm, suitability, suit_other, vulnerability,
                                    vulnerability_other, predation, avgN, diet_hat) {
  ncell <- msm$ncell; nyrs <- msm$nyrs; nyrs_hind <- msm$nyrs_hind
  # Suitability in every year, [cell, cell, year]; projection years repeat the last hindcast
  suit_all <- if (length(dim(suitability)) == 2) {
    array(suitability, c(ncell, ncell, nyrs))
  } else {
    suitability[, , c(seq_len(nyrs_hind), rep(nyrs_hind, nyrs - nyrs_hind)), drop = FALSE]
  }
  # Predator cell, prey cell and year of each element of suit_all, and its position in the
  # C++ 5-d pair reports (as5); as4 places [cell, year] into [nspp, 1, max_age, nyrs]
  pred <- rep(seq_len(ncell), ncell * nyrs)
  prey <- rep(rep(seq_len(ncell), each = ncell), nyrs)
  yr <- rep(seq_len(nyrs), each = ncell * ncell)
  pair <- rep(msm$pair_lin, nyrs) + msm$n_pair4 * (yr - 1)
  as5 <- function(v) { a <- array(0, msm$dim5); a[pair] <- v; a }
  as4 <- function(m) { a <- array(0, msm$dim4); a[msm$lin4] <- m; a }
  prey_B <- avgN * msm$W
  # M2_prop: a predator age's part of its prey's M2; B_eaten: the prey biomass (mt) it eats
  M2_prop <- predation$pred_mort_scalar[cbind(pred, yr)] * suit_all
  B_eaten <- prey_B[cbind(prey, yr)] * M2_prop
  eaten <- B_eaten
  dim(eaten) <- c(ncell, ncell * nyrs)
  B_eaten_as_prey <- matrix(colSums(eaten), ncell, nyrs)   # summed over predators
  stom_div_bio <- array(0, msm$dim5)
  if (!is.null(msm$empirical)) {
    emp <- msm$empirical
    prey_N_response <- if (msm$msmMode == 2) avgN * avgN else avgN
    prey_denom <- msm$W * prey_N_response
    has_prey <- 1 - .rtmb_ind_le(prey_denom, 0)
    inv_prey_denom <- has_prey / .rtmb_safe_den(prey_denom, has_prey)
    sdb <- inv_prey_denom[emp$prey, , drop = FALSE] * emp$diet_prop
    stom_div_bio[msm$pair_lin[cbind(emp$pred, emp$prey)] + msm$n_pair4 * (col(sdb) - 1)] <- sdb
  }
  dh <- matrix(0, msm$n_diet, 2)
  if (!is.null(diet_hat)) dh[, 2] <- diet_hat
  list(suitability = as5(as.vector(suit_all)), M2_prop = as5(as.vector(M2_prop)),
       B_eaten = as5(as.vector(B_eaten)), stom_div_bio = stom_div_bio,
       B_eaten_as_prey = as4(B_eaten_as_prey), M2_at_age = as4(predation$M2),
       suit_other = as4(matrix(suit_other, ncell, nyrs)),
       avail_food = as4(predation$avail_food),
       vulnerability = vulnerability, vulnerability_other = vulnerability_other, diet_hat = dh)
}
