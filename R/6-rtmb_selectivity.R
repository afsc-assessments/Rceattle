# Selectivity-at-age for the RTMB model; mirrors src/TMB/selectivity.hpp
# (calculate_selectivity and normalize_and_project_selectivity) for the forms ported so
# far: 0 Fixed, 1 Logistic and 2 NonParametric (Ianelli et al. 2018, AMAK), age-based,
# one sex, without normalization (.rtmb_check_scope() refuses the rest). The AD mechanics
# are in 6-rtmb_helpers.R.


#' Selectivity-at-age by fleet (unitless, 0 to 1 for Logistic)
#'
#' @details An age-based curve is evaluated at the bin number (bin + 1), not the age, so
#'   `sel_inf` is on the bin scale; it equals age only when minage = 1. Hindcast years are
#'   computed; projection years hold the last hindcast year (for Fixed, the last year with
#'   an emp_sel row).
#' @return List: `sel`, per fleet a matrix `[age, year]`; and per NonParametric fleet, for
#'   the penalties, the first year's log curve `log_sel_yr1` and `avg_sel`, AMAK's avgsel:
#'   the log mean of exp(coefficients) over the estimated bins, penalized as 2 avg_sel^2.
#' @noRd
.rtmb_selectivity <- function(p, pd) {
  "[<-" <- RTMB::ADoverload("[<-")
  H <- pd$nyrs_hind
  proj <- seq_len(pd$nyrs - H) + H
  sel <- log_sel_yr1 <- avg_sel <- vector("list", pd$n_flt)
  for (f in seq_len(pd$n_flt)) {
    na <- pd$nages[pd$flt_spp[f]]
    type <- pd$flt_sel_type[f]
    if (type == 0) {
      sel[[f]] <- pd$sel_fixed[[f]]
      next
    }
    if (type == 1) {
      # Logistic: 1 / (1 + exp(-slope (bin - inflection))), per hindcast year
      slope <- exp(p[["log_sel_slp"]][1, f, 1] + p[["log_sel_slp_dev"]][1, f, 1, seq_len(H)])
      inflection <- p[["sel_inf"]][1, f, 1] + p[["sel_inf_dev"]][1, f, 1, seq_len(H)]
      z <- (RTMB::matrix(seq_len(na), na, H) - rep(inflection, each = na)) *
        rep(slope, each = na)
      sel_hind <- RTMB::matrix(1 / (1 + exp(-z)), na, H)
    }
    if (type == 2) {
      # NonParametric: one log coefficient per estimated bin, the last carried to the
      # older bins, centred so the curve's mean on the natural scale is 1
      n_bins <- pd$flt_n_sel_bins[f]
      log_coef <- RTMB::matrix(p[["sel_coff"]][f, 1, seq_len(n_bins)], n_bins, H) +
        RTMB::matrix(p[["sel_coff_dev"]][f, 1, seq_len(n_bins), seq_len(H)], n_bins, H)
      log_curve <- log_coef[pmin(seq_len(na), n_bins), , drop = FALSE]
      log_sel <- log_curve - rep(.rtmb_log_mean_exp(log_curve), each = na)
      sel_hind <- exp(log_sel)
      log_sel_yr1[[f]] <- log_sel[, 1]
      avg_sel[[f]] <- .rtmb_log_mean_exp(log_coef[, 1, drop = FALSE])
    }
    sel[[f]] <- if (length(proj)) {
      cbind(sel_hind, RTMB::matrix(sel_hind[, H], na, length(proj)))
    } else {
      sel_hind
    }
  }
  list(sel = sel, log_sel_yr1 = log_sel_yr1, avg_sel = avg_sel)
}
