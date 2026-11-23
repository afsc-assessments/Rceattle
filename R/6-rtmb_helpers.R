# Helpers for the RTMB model, and the one place its automatic-differentiation (AD)
# mechanics are explained. Most mirror a function in src/TMB/helper_functions.hpp or
# spr.hpp in its evaluation order, so the objective matches the C++ to rounding.
#
# RTMB records the objective's arithmetic on estimated quantities once, as a "tape", and
# differentiates the tape. Three rules follow from that, and none of them is ecology:
# - A branch on an estimated value is a 0/1 indicator (.rtmb_ind_le/gt/ge) times each
#   branch, never `if` or `ifelse()`: an `if` errors at taping, or under the taped setting
#   silently takes one branch whatever the value. The tape re-decides the indicator at
#   every point, and .rtmb_safe_den() keeps the branch it switches off finite.
# - A function that writes into a container or calls c() on the tape declares `"[<-"` and
#   `c` from RTMB::ADoverload(): package code reaches base's versions before RTMB's, and
#   they drop the AD class, c(1, <advector>) included.
# - .rtmb_matmul() sums in the C++ order, which `%*%` does not; .rtmb_group_sum() sums
#   terms into groups, which an indexed `+=` cannot do in R; .rtmb_max2() and
#   .rtmb_posfun() are the C++'s branch-free max and floor.
# Also: on the tape RTMB takes dim() from the first of two equal-length operands, so a matrix
# leads each product; base's matrix(), colSums(), rowSums() come first, so RTMB:: is written.


#' 0/1 indicator of `x >= a`, recorded on the AD tape
#'
#' @details Turns on RTMB's taped comparisons for this one comparison and restores the
#'   previous setting, so a stray `if` anywhere else in the objective still errors.
#'   On plain numbers (the `report()` pass) it returns a logical.
#' @noRd
.rtmb_ind_ge <- function(x, a) {
  old <- RTMB::TapeConfig()
  on.exit(RTMB::TapeConfig(old))
  RTMB::TapeConfig(comparison = "tape")
  x >= a
}

#' 0/1 indicator of `x <= a`, recorded on the AD tape
#'
#' @details As `.rtmb_ind_ge()`. NaN compares as 0 on the tape, and as FALSE here in the
#'   `report()` pass, so `1 - .rtmb_ind_le(x, 0)` is the C++ `!(x <= 0)`, which lets NaN
#'   through.
#' @noRd
.rtmb_ind_le <- function(x, a) {
  if (!inherits(x, "advector")) return(!is.na(x) & x <= a)
  old <- RTMB::TapeConfig()
  on.exit(RTMB::TapeConfig(old))
  RTMB::TapeConfig(comparison = "tape")
  x <= a
}

#' 0/1 indicator of `x > a`, recorded on the AD tape; FALSE for NaN, as `.rtmb_ind_le()`
#' @noRd
.rtmb_ind_gt <- function(x, a) {
  if (!inherits(x, "advector")) return(!is.na(x) & x > a)
  old <- RTMB::TapeConfig()
  on.exit(RTMB::TapeConfig(old))
  RTMB::TapeConfig(comparison = "tape")
  x > a
}

#' Denominator that is `den` where `m` is 1 and 1 where `m` is 0
#'
#' @details Keeps the branch an indicator switches off finite, so 0 times that branch
#'   is 0 rather than NaN.
#' @noRd
.rtmb_safe_den <- function(den, m) den * m + (1 - m)

#' Smooth floor at `eps` (helper_functions.hpp `posfun`)
#'
#' @details Returns `x` at or above `eps`, otherwise `eps / (2 - x / eps)`, plus the
#'   penalty `0.01 (x - eps)^2` below `eps`. Vectorized; both elements match the C++
#'   CondExpGe form, evaluation order included.
#' @return List: `value`, and `penalty` per element.
#' @noRd
.rtmb_posfun <- function(x, eps) {
  m <- .rtmb_ind_ge(x, eps)
  b <- 1 - m
  # The indicator multiplies (x - eps) before it is squared, so a huge x gives a penalty of
  # 0, not 0 * Inf = NaN; below eps it is the C++'s 0.01 (x - eps) (x - eps) exactly
  d <- (x - eps) * b
  list(value = m * x + b * (eps / .rtmb_safe_den(2 - x / eps, b)),
       penalty = 0.01 * d * d)
}

#' Branch-free maximum, exact with a kink at x = y (helper_functions.hpp `max2`)
#' @noRd
.rtmb_max2 <- function(x, y) 0.5 * (abs(x - y) + x + y)

#' Group sum `G %*% x` that keeps the AD class on the tape
#'
#' @details Sums terms into their groups, which an indexed `+=` cannot do in R: a
#'   repeated index keeps only its last write.
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

#' Log of the mean of `exp(x)` over the rows of each column (helper_functions.hpp
#' `log_mean_exp`)
#'
#' @details The shift is the running `max2` down each column, as in the C++.
#' @param x Matrix; one log-mean-exp per column.
#' @noRd
.rtmb_log_mean_exp <- function(x) {
  n <- nrow(x)
  mx <- x[1, ]
  for (i in seq_len(n)[-1]) mx <- .rtmb_max2(mx, x[i, ])
  total <- RTMB::colSums(exp(x - rep(mx, each = n)))
  mx + log(total / n)
}

#' Survivors per recruit at age (spr.hpp `per_recruit_survivors`)
#'
#' @details One recruit at the first age, surviving at total mortality `Z`; the last
#'   age is a plus group.
#' @noRd
.rtmb_per_recruit_survivors <- function(Z) {
  c <- RTMB::ADoverload("c")
  na <- length(Z)
  n <- cumprod(c(1, exp(-Z[seq_len(na - 2)])))
  c(n, n[na - 1] * exp(-Z[na - 1]) / (1 - exp(-Z[na])))
}

#' Spawning biomass per recruit (spr.hpp `spawning_biomass_per_recruit`)
#'
#' @param spawn_month Month of spawning; mortality before it discounts the spawners.
#' @noRd
.rtmb_spawning_biomass_per_recruit <- function(n, Z, weight, mature_female, spawn_month) {
  sum(n * weight * mature_female * exp(-Z * spawn_month / 12))
}

#' Lognormal-sd estimate concentrated out of the likelihood (helper_functions.hpp
#' `concentrated_lognormal_sd`; Walters and Ludwig 1994)
#' @noRd
.rtmb_concentrated_lognormal_sd <- function(mean_sq_resid, bias_adjust) {
  sqrt(2 * mean_sq_resid / (sqrt(1 + (bias_adjust * bias_adjust) * mean_sq_resid) + 1))
}
