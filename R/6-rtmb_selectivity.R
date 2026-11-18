# Selectivity-at-age for the RTMB model; mirrors src/TMB/selectivity.hpp
# (calculate_selectivity and normalize_and_project_selectivity) for the forms ported so
# far: 0 Fixed, 1 Logistic and 2 NonParametric (Ianelli et al. 2018, AMAK), age-based,
# one sex, without normalization (.rtmb_check_scope() refuses the rest).


#' Selectivity-at-age by fleet
#'
#' @details An age-based curve is evaluated at the bin number (bin + 1), not the age.
#'   Hindcast years are computed; projection years hold the last hindcast year (for
#'   Fixed, the last year with an emp_sel row).
#' @return List: `sel`, per fleet a matrix [age, year]; and per NonParametric fleet the
#'   first year's log curve `log_np` and level `avg_sel`, which the penalties read.
#' @noRd
.rtmb_selectivity <- function(p, pd) {
  "[<-" <- RTMB::ADoverload("[<-")
  H <- pd$nyrs_hind
  proj <- seq_len(pd$nyrs - H) + H
  sel <- log_np <- avg_sel <- vector("list", pd$n_flt)
  for (f in seq_len(pd$n_flt)) {
    na <- pd$nages[pd$flt_spp[f]]
    type <- pd$flt_sel_type[f]
    if (type == 0) {
      sel[[f]] <- pd$sel_fixed[[f]]
      next
    }
    if (type == 1) {
      # Logistic: 1 / (1 + exp(-slope (bin - inflection))), per hindcast year
      slp <- exp(p[["log_sel_slp"]][1, f, 1] + p[["log_sel_slp_dev"]][1, f, 1, seq_len(H)])
      inf <- p[["sel_inf"]][1, f, 1] + p[["sel_inf_dev"]][1, f, 1, seq_len(H)]
      z <- (RTMB::matrix(seq_len(na), na, H) - rep(inf, each = na)) * rep(slp, each = na)
      s <- RTMB::matrix(1 / (1 + exp(-z)), na, H)
    }
    if (type == 2) {
      # NonParametric: one log coefficient per estimated bin, the last carried to the
      # older bins, centred so the curve's mean on the natural scale is 1
      nsb <- pd$flt_n_sel_bins[f]
      coef <- RTMB::matrix(p[["sel_coff"]][f, 1, seq_len(nsb)], nsb, H) +
        RTMB::matrix(p[["sel_coff_dev"]][f, 1, seq_len(nsb), seq_len(H)], nsb, H)
      log_curve <- coef[pmin(seq_len(na), nsb), , drop = FALSE]
      lnp <- log_curve - rep(.rtmb_log_mean_exp(log_curve), each = na)
      s <- exp(lnp)
      log_np[[f]] <- lnp[, 1]
      avg_sel[[f]] <- .rtmb_log_mean_exp(coef[, 1, drop = FALSE])
    }
    sel[[f]] <- if (length(proj)) cbind(s, RTMB::matrix(s[, H], na, length(proj))) else s
  }
  list(sel = sel, log_np = log_np, avg_sel = avg_sel)
}
