# Multispecies predation for the RTMB model: Type II and Type III MSVPA (Holsman et al.
# 2015), ported from src/TMB/predation.hpp and diet_data.hpp on 5.67.11 (line ranges cited
# per function). Quantities are on the one-sex (species, age) cells of
# .rtmb_prepare_predation(), [cell, year]; a cell the C++ never writes stays a constant 0.


#' Group sum `G %*% x` that keeps the AD class on the tape
#'
#' @param G Sparse 0/1 matrix, one row per group.
#' @param x Vector, or matrix with one row per summed term; terms add in row order.
#' @noRd
.rtmb_group_sum <- function(G, x) {
  out <- G %*% x
  if (!inherits(out, "advector")) out <- as.matrix(out)
  if (is.null(dim(x))) out <- out[, 1]
  out
}

#' Matrix product summed term by term in column order of `A`, as the C++ loops sum
#'
#' @details `%*%` sums in its own order: BLAS on plain numbers (the `report()` pass), and
#'   on the tape an order that differs from the C++ in the last bit. Either is amplified
#'   where it matters: near the 1e-6 other-food floor, or through exp(-sum(Z)) when M2 is
#'   large (4e-10 in a projected SSB of 4e-22 mt). A product with a constant zero folds
#'   on the tape, so this records no more operations than `%*%`.
#' @noRd
.rtmb_matmul <- function(A, B) {
  out <- 0
  for (k in seq_len(ncol(A))) {
    out <- out + A[, k, drop = FALSE] %*% B[k, , drop = FALSE]
  }
  out
}

#' Parametric suitability, lognormal on the predator-to-prey weight ratio (suitMode 4)
#'
#' @details predation.hpp:20-45 (vulnerability, a multinomial logit of `log_phi`) and
#'   :76-142 (suitability). A pair whose log weight ratio is not positive keeps
#'   suitability 0. Weights vary by hindcast year, so suitability is [cell, cell, year];
#'   projection years repeat the last hindcast year, as their weights do.
#' @return List: `S` [nc, nc, nyrs_hind], `suit_other` [nc], `vulnerability`,
#'   `vulnerability_other`.
#' @noRd
.rtmb_suit_parametric <- function(pm, log_phi, gam_a, gam_b) {
  "[<-" <- RTMB::ADoverload("[<-")
  nspp <- pm$nspp
  vul <- matrix(0, nspp, nspp)
  vul_other <- numeric(nspp)
  for (rsp in pm$param_sp) {
    e <- exp(log_phi[rsp, ])
    vul[rsp, ] <- e / (1 + sum(e))
    vul_other[rsp] <- 1 - sum(vul[rsp, ])
  }
  ps <- pm$ps
  dens <- vul[cbind(ps$sp_r, ps$sp_k)] * RTMB::dnorm(ps$sr, gam_a[ps$sp_r], gam_b[ps$sp_r]) /
    RTMB::dnorm(gam_a, gam_a, gam_b)[ps$sp_r]
  S <- array(0, c(pm$nc, pm$nc, pm$nyrs_hind))
  S[ps$lin] <- dens
  suit_other <- vul_other[pm$cell_sp]
  list(S = S, suit_other = suit_other, vulnerability = vul, vulnerability_other = vul_other)
}

#' Empirical (MSVPA) suitability from stomach proportions and prey biomass (suitMode 0)
#'
#' @details predation.hpp:179-308. Diet proportion over prey biomass (squared prey
#'   numbers under Type III), normalized over the modelled prey and "other food", the
#'   latter floored at 1e-6 of the prey term, then averaged over each predator's
#'   suitability window and held for every year.
#' @param Nc Average numbers-at-age (thousands of fish), [nc, nyrs].
#' @return List: `S` [nc, nc], `suit_other` [nc].
#' @noRd
.rtmb_suit_empirical <- function(pm, Nc) {
  "[<-" <- RTMB::ADoverload("[<-")
  e <- pm$emp
  u <- e$cols
  prey_N <- Nc[, u, drop = FALSE]
  if (pm$msmMode == 2) prey_N <- prey_N * prey_N
  prey_denom <- pm$W[, u, drop = FALSE] * prey_N
  ok <- 1 - .rtmb_ind_le(prey_denom, 0)    # NaN passes (predation.hpp:229)
  inv <- ok / .rtmb_safe_den(prey_denom, ok)
  sdb <- inv[e$k, , drop = FALSE] * e$d   # stom_div_bio on the diet triplets
  suma <- .rtmb_group_sum(e$G, sdb)
  # The diet proportion is summed only where the prey test passed, as in the C++
  dps <- .rtmb_group_sum(e$G, ok[e$k, , drop = FALSE] * e$d)
  other <- .rtmb_max2(matrix(pm$ofdp, pm$nc, length(u)), 1e-6 * suma)
  inv_pred <- .rtmb_ind_gt(dps, 0) / .rtmb_max2(suma + other, 1e-10)
  prod <- sdb * inv_pred[e$r, , drop = FALSE]
  S <- matrix(0, pm$nc, pm$nc)
  for (w in e$by_pred) {
    S[cbind(e$r[w$t], e$k[w$t])] <- RTMB::rowSums(prod[w$t, w$cols, drop = FALSE]) / w$nyrs
  }
  list(S = S, suit_other = 1 - RTMB::rowSums(S))
}

#' Available food and predation mortality M2 (MSVPA)
#'
#' @details predation.hpp:349-528. Available food (mt) is suitability-weighted prey
#'   biomass plus other food; under Type III the prey numbers are squared there but not
#'   in M2. The same is done for the dynamic B0 (F = 0) and BF (F = Ftarget) numbers, with
#'   the fit's suitability. A predator age whose available food is not positive adds no
#'   M2; NaN available food passes that test and gives NaN.
#' @param S Suitability, [nc, nc] or [nc, nc, nyrs_hind].
#' @param Nc,N0c,NFc Average numbers-at-age, [nc, nyrs]: fished, dynamic B0, dynamic BF.
#' @return List of [nc, nyrs] matrices: `avail`, `inv_avail`, `pms` (predator scalar),
#'   `M2`, `M2_dB0`, `M2_dBF`.
#' @noRd
.rtmb_predation <- function(pm, S, suit_other, Nc, N0c, NFc) {
  Y <- pm$nyrs; H <- pm$nyrs_hind
  adj <- function(N) if (pm$msmMode == 2) N * N else N
  B <- cbind(adj(Nc) * pm$W, adj(N0c) * pm$W, adj(NFc) * pm$W)
  # S %*% B over the three copies; a year-varying S multiplies year by year, and the
  # projection years share the last hindcast year's suitability
  by_year <- function(Sm, X, tr) {
    if (length(dim(Sm)) == 2) return(if (tr) .rtmb_matmul(t(Sm), X) else .rtmb_matmul(Sm, X))
    out <- vector("list", Y)
    for (y in seq_len(H)) {
      Sy <- Sm[, , y]
      out[[y]] <- if (tr) .rtmb_matmul(t(Sy), X[, c(y, Y + y, 2 * Y + y), drop = FALSE]) else
        .rtmb_matmul(Sy, X[, c(y, Y + y, 2 * Y + y), drop = FALSE])
    }
    res <- do.call(cbind, out[seq_len(H)])
    res <- res[, c(seq(1, 3 * H, 3), seq(2, 3 * H, 3), seq(3, 3 * H, 3)), drop = FALSE]
    if (Y > H) {
      pj <- seq_len(Y - H) + H
      Sy <- Sm[, , H]
      P <- if (tr) .rtmb_matmul(t(Sy), X[, c(pj, Y + pj, 2 * Y + pj), drop = FALSE]) else
        .rtmb_matmul(Sy, X[, c(pj, Y + pj, 2 * Y + pj), drop = FALSE])
      n <- Y - H
      res <- cbind(res[, seq_len(H), drop = FALSE], P[, seq_len(n), drop = FALSE],
                   res[, H + seq_len(H), drop = FALSE], P[, n + seq_len(n), drop = FALSE],
                   res[, 2 * H + seq_len(H), drop = FALSE], P[, 2 * n + seq_len(n), drop = FALSE])
    }
    res
  }
  A <- by_year(S, B, FALSE)
  other <- rep(pm$other_food_c * suit_other, Y)
  avail <- A[, seq_len(Y), drop = FALSE] + other
  avail0 <- A[, Y + seq_len(Y), drop = FALSE] + other
  availF <- A[, 2 * Y + seq_len(Y), drop = FALSE] + other
  # pms: avgN (thousands) x ration (kg/predator/yr) is mt consumed, over available food (mt)
  if (is.null(pm$ps)) {
    # Empirical suitability floors suit_other at 1e-6/(1 + 1e-6) and other food is
    # positive, so available food is never <= 0 and the C++ test always passes
    inv_avail <- 1 / avail
    pms <- Nc * pm$Rat * inv_avail
    pms0 <- (N0c * pm$Rat) / avail0
    pmsF <- (NFc * pm$Rat) / availF
  } else {
    # Parametric: suit_other is 1 - sum(vulnerability), which rounds to 0 or just below at
    # large log_phi, leaving a predator age with no smaller prey no available food
    m <- 1 - .rtmb_ind_le(avail, 0)
    m0 <- (1 - .rtmb_ind_le(avail0, 0)) * m
    mF <- (1 - .rtmb_ind_le(availF, 0)) * m
    inv_avail <- m / .rtmb_safe_den(avail, m)
    pms <- Nc * pm$Rat * inv_avail
    pms0 <- m0 * ((N0c * pm$Rat) / .rtmb_safe_den(avail0, m0))
    pmsF <- mF * ((NFc * pm$Rat) / .rtmb_safe_den(availF, mF))
  }
  M2 <- by_year(S, cbind(pms, pms0, pmsF), TRUE)
  list(avail = avail, inv_avail = inv_avail, pms = pms, M2 = M2[, seq_len(Y), drop = FALSE],
       M2_dB0 = M2[, Y + seq_len(Y), drop = FALSE], M2_dBF = M2[, 2 * Y + seq_len(Y), drop = FALSE])
}

#' Predicted stomach proportions on the diet rows (diet-data type 1)
#'
#' @details diet_data.hpp:160-235: the predicted diet proportion of each row's prey age
#'   in its predator age, averaged over the predator's suitability window.
#' @return Vector over the diet rows (column 2 of the C++ `diet_hat`).
#' @noRd
.rtmb_diet_hat <- function(pm, S, Nc, inv_avail) {
  "[<-" <- RTMB::ADoverload("[<-")
  NW <- Nc * pm$W
  out <- numeric(pm$n_diet)
  for (g in pm$dh) {
    dph <- NW[g$ky] * S[g$rky] * inv_avail[g$ry]
    Nr <- Nc[g$ry]
    gt <- .rtmb_ind_gt(Nr, 0)
    term <- gt * (dph * Nr) / .rtmb_safe_den(Nr, gt)
    dim(term) <- c(length(g$rows), g$nyrs)
    out[g$rows] <- RTMB::rowSums(term) / g$nyrs
  }
  out
}

#' Stomach-content likelihood: multinomial on prey proportions plus "other prey"
#'
#' @details ceattle.cpp section 13.6 (3753-3817), `Diet_distribution` 0. The predicted
#'   other-prey share is floored at 1e-5 by posfun, whose penalty the C++ discards; both
#'   sides get `comp_offset` and are renormalized, and the log-likelihood is TMB's
#'   `dmultinom()`, weighted by the raw `diet_comp_weights`.
#' @return List: `weighted` and `unweighted` negative log-likelihood by predator species.
#' @noRd
.rtmb_stomach_nll <- function(pm, diet_hat, diet_comp_weights) {
  "[<-" <- RTMB::ADoverload("[<-")
  st <- pm$st
  off <- pm$comp_offset
  q <- diet_hat[st$rows]
  q_other <- .rtmb_posfun(1 - .rtmb_group_sum(st$G, q), 1e-5)$value
  tot <- .rtmb_group_sum(st$G, q + off) + (q_other + off)
  ll <- st$const + (.rtmb_group_sum(st$G, st$x * log((q + off) / tot[st$sid])) +
                      st$x_other * log((q_other + off) / tot))
  weighted <- unweighted <- numeric(pm$nspp)
  for (sp in unique(st$sp)) {
    j <- st$sp == sp
    weighted[sp] <- -sum(diet_comp_weights[sp] * ll[j])
    unweighted[sp] <- -sum(ll[j])
  }
  list(weighted = weighted, unweighted = unweighted)
}

#' The predation REPORTs at the C++ names and shapes, from the last iteration
#'
#' @details Plain R, run only in the `report()` pass. Each predator-prey product is formed
#'   as the C++ forms it (predation.hpp:496-526). At a NaN objective the C++ report can
#'   hold NaN in more cells, from 0 x NaN products this port never forms.
#' @noRd
.rtmb_predation_reports <- function(pm, S, suit_other, vul, vul_other, pred, Nc, diet_hat) {
  nc <- pm$nc; Y <- pm$nyrs; H <- pm$nyrs_hind
  yrs <- if (length(dim(S)) == 3) c(seq_len(H), rep(H, Y - H)) else NULL
  S_full <- if (is.null(yrs)) array(S, c(nc, nc, Y)) else S[, , yrs, drop = FALSE]
  rr <- rep(seq_len(nc), nc * Y)
  kk <- rep(rep(seq_len(nc), each = nc), Y)
  yy <- rep(seq_len(Y), each = nc * nc)
  pair <- rep(pm$pair_lin, Y) + pm$n_pair4 * (yy - 1)
  as5 <- function(v) { a <- array(0, pm$dim5); a[pair] <- v; a }
  as4 <- function(m) { a <- array(0, pm$dim4); a[pm$lin4] <- m; a }
  NW <- Nc * pm$W
  M2_prop <- pred$pms[cbind(rr, yy)] * S_full
  B_eaten <- NW[cbind(kk, yy)] * M2_prop
  eaten <- B_eaten
  dim(eaten) <- c(nc, nc * Y)
  B_eaten_as_prey <- matrix(colSums(eaten), nc, Y)
  stom_div_bio <- array(0, pm$dim5)
  if (!is.null(pm$emp)) {
    e <- pm$emp
    prey_N <- if (pm$msmMode == 2) Nc * Nc else Nc
    prey_denom <- pm$W * prey_N
    ok <- 1 - .rtmb_ind_le(prey_denom, 0)
    inv <- ok / .rtmb_safe_den(prey_denom, ok)
    sdb <- inv[e$k, , drop = FALSE] * e$d
    stom_div_bio[pm$pair_lin[cbind(e$r, e$k)] + pm$n_pair4 * (col(sdb) - 1)] <- sdb
  }
  dh <- matrix(0, pm$n_diet, 2)
  if (!is.null(diet_hat)) dh[, 2] <- diet_hat
  list(suitability = as5(as.vector(S_full)), M2_prop = as5(as.vector(M2_prop)),
       B_eaten = as5(as.vector(B_eaten)), stom_div_bio = stom_div_bio,
       B_eaten_as_prey = as4(B_eaten_as_prey), M2_at_age = as4(pred$M2),
       suit_other = as4(matrix(suit_other, nc, Y)), avail_food = as4(pred$avail),
       vulnerability = vul, vulnerability_other = vul_other, diet_hat = dh)
}
