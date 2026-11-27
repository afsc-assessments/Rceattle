#' @references
#' \itemize{
#'   \item Adams, G.D., Holsman, K.K., Barbeaux, S.J., Dorn, M.W., Ianelli, J.N.,
#'         Spies, I., Stewart, I.J., Punt, A.E. (2022). An ensemble approach to
#'         understand predation mortality for groundfish in the Gulf of Alaska.
#'         Fisheries Research, 251, 106303. doi:10.1016/j.fishres.2022.106303
#'   \item Holsman, K.K., Ianelli, J., Aydin, K., Punt, A.E., Moffitt, E.A. (2016).
#'         A comparison of fisheries biological reference points estimated from
#'         temperature-specific multi-species and single-species climate-enhanced
#'         stock assessment models. Deep Sea Research Part II: Topical Studies in
#'         Oceanography, 134, 360-378.
#'   \item Wassermann, S.N., Adams, G.D., Haltuch, M.A., Kaplan, I.C., Marshall,
#'         K.N., Punt, A.E. (2025). Even low levels of cannibalism can bias
#'         population estimates for Pacific hake. ICES Journal of Marine Science,
#'         82(1), fsae064.
#'   \item Kristensen, K., Nielsen, A., Berg, C.W., Skaug, H., Bell, B.M. (2016).
#'         TMB: automatic differentiation and Laplace approximation. Journal of
#'         Statistical Software, 70(5). doi:10.18637/jss.v070.i05.
#' }
#' @rawNamespace useDynLib(Rceattle, .registration=TRUE); useDynLib(ceattle)
#' @keywords internal
#' @importFrom dplyr filter mutate select arrange group_by summarise summarize rename left_join right_join full_join bind_rows bind_cols distinct pull n lag vars matches
#' @importFrom plyr rbind.fill
#' @importFrom stats sd median nlminb dnorm qlnorm aggregate rgamma rmultinom coef logLik residuals vcov simulate
#' @importFrom graphics mtext axis text symbols filled.contour persp curve points image box
#' @importFrom grDevices col2rgb rgb
#' @importFrom utils write.csv type.convert packageVersion
#' @importFrom ggplot2 geom_tile scale_y_continuous coord_equal scale_x_continuous element_rect element_text scale_fill_viridis_c geom_contour geom_contour_filled scale_fill_viridis_d geom_line scale_color_viridis_c
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL

utils::globalVariables(c(
  # Data column names used in dplyr/ggplot2 non-standard evaluation
  "Year", "Fleet_code", "Time_varying_q", "Catch", "Species", "Length", "Bin", "Sex",
  "Pred", "Pred_age", "Pred_sex", "Sample_size", "Diet_hat",
  "Stomach_proportion_by_weight", "Diet_weights_mcallister",
  "Prey", "Prey_sex", "Prey_age", "stratum_id", "stomach_id",
  "Selectivity_index", "Catchability", "Q_index", "Catchability_index", "Fleet_type",
  "Selectivity", "Comp_loglike", "CAAL_loglike", "Index_loglike",
  "Comp_distribution", "CAAL_distribution", "Index_distribution", "Wt_index", "Tmp_ind",
  "Wt_name", "Age_transition_name", "Age_transition_index",
  "Age", "True_age", "M",
  "Observation", "Log_sd", "Age0_Length1", "Month", "Length_bin",
  "Estimate_survey_sd", "Survey_sd_prior", "Nselages", "Sel_sd_prior",
  "Estimate_q", "Age_first_selected", "Age_max_selected",
  "Age_max_selected_upper",
  "Time_varying_sel", "Time_varying_sel_sd_prior", "Time_varying_sel_sd",
  "Time_varying_q_sd_prior", "Time_varying_q_sd",
  "Q_prior", "Q_init", "Catchability_init", "Q_sd_prior", "Catchability_prior_sd",
  "Index_sd_prior", "Index_sd", "Catch_sd_prior", "Catch_sd",
  "Sel_curve_pen1", "Sel_curve_pen2",
  # Loop iterators and anonymous function variables
  "i", "x", "spp", "sim", "var"
))
