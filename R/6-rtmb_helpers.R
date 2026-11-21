# Helpers for the RTMB model. Most mirror a function in src/TMB/helper_functions.hpp or
# spr.hpp in its evaluation order, so the objective matches the C++ to rounding; the
# .rtmb_ind_*() indicators and .rtmb_safe_den() are the RTMB form of the C++ branches.
#
# A branch on an estimated quantity is a 0/1 indicator times each branch, never `if` or
# `ifelse()`: an `if` on an AD value errors under RTMB's default comparison setting, and
# under the taped setting silently takes one branch whatever the value.
#
# A function that writes into a container or calls c() on the tape declares `"[<-"` and
# `c` from RTMB::ADoverload(): package code reaches base's versions before RTMB's, and
# they drop the AD class, c(1, <advector>) included.


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
