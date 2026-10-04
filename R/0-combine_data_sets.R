#' Combine two Rceattle data lists
#'
#' @description
#' Concatenate two data lists into one. `env_data` is taken from `data_list1`,
#' and diet data must be updated afterward.
#'
#' @param data_list1 Rceattle data_list 1
#' @param data_list2 Rceattle data_list 2
#'
#' @export
#'
combine_data <- function(data_list1 = NULL, data_list2 = NULL){

  # Fold the deprecated `est_M1` into `M1_model` on both inputs first, so a
  # legacy (est_M1) list combines correctly with a modern (M1_model) one -- both
  # then concatenate under the single `M1_model` name below.
  data_list1 <- .alias_est_M1(data_list1)
  data_list2 <- .alias_est_M1(data_list2)

  # Upgrade deprecated fleet_control column names on both inputs, so a legacy
  # (e.g. Q_index) list combines correctly with a modern (Catchability_index)
  # one under the single canonical names read below.
  if (!is.null(data_list1))
    data_list1$fleet_control <- .rce_upgrade_fleet_control_aliases(data_list1$fleet_control)
  if (!is.null(data_list2))
    data_list2$fleet_control <- .rce_upgrade_fleet_control_aliases(data_list2$fleet_control)
  # ...and the control / bioenergetics element names (sigma_rec, Diet_distribution)
  # so both inputs concatenate under the single canonical vec_names below.
  data_list1 <- .rce_upgrade_data_list_aliases(data_list1)
  data_list2 <- .rce_upgrade_data_list_aliases(data_list2)

  data_list_new <- data_list1

  dat_names <- names(data_list_new)

  names_not_used <- c("nspp", "styr", "endyr", "projyr")

  vec_names <- c("spnames", "nsex", "spawn_month", "nages", "minage", "nlengths","pop_wt_index", "ssb_wt_index", "M1_model", "M1_re", "pop_alk_index", "sigma_rec", "other_food", "estDynamics", "Ceq", "Cindex","Pvalue", "fday", "CA","CB", "Qc", "Tco",  "Tcm",  "Tcl",  "CK1", "CK4", "beta_wt_len", "alpha_wt_len", "Diet_distribution", "Diet_comp_weights") # Object names of per-species vectors

  # equil_catch_data is here so a cleaned data set keeps its initial equilibrium
  # catch, and gets the same Fleet_code and Species reindexing as catch_data.
  mat_names <- c("fleet_control", "index_data", "catch_data", "equil_catch_data", "comp_data", "caal_data", "emp_sel", "NByageFixed", "age_trans_matrix", "age_error", "weight",   "maturity", "sex_ratio", "M1_base", "ration_data", "diet_data") # Object names of matrices

  # Get index from data_set1 of the 4 indices
  fleet_index1 <- max(data_list1$fleet_control$Fleet_code, na.rm = TRUE)
  weight_index1 <- max(data_list1$weight$Wt_index, na.rm = TRUE)
  alk_index1 <- max(data_list1$age_trans_matrix$Age_transition_index, na.rm = TRUE)
  nspp1 <- data_list1$nspp
  q_index1 <- max(data_list1$fleet_control$Catchability_index, na.rm = TRUE)
  sel_index1 <- max(data_list1$fleet_control$Selectivity_index, na.rm = TRUE)

  # Update vector indices
  data_list2[["pop_wt_index"]] <- data_list2[["pop_wt_index"]] + weight_index1
  data_list2[["ssb_wt_index"]] <- data_list2[["ssb_wt_index"]] + weight_index1
  data_list2[["pop_age_transition_index"]] <- data_list2[["pop_age_transition_index"]] + alk_index1

  # Update fleet control
  data_list2$fleet_control$Age_transition_index <- data_list2$fleet_control$Age_transition_index + alk_index1
  data_list2$fleet_control$Selectivity_index <- data_list2$fleet_control$Selectivity_index + sel_index1
  data_list2$fleet_control$Catchability_index <- data_list2$fleet_control$Catchability_index + q_index1
  data_list2$fleet_control$Species <- data_list2$fleet_control$Species + nspp1
  data_list2$fleet_control$Weight_index <- data_list2$fleet_control$Weight_index + weight_index1
  data_list2$fleet_control$Fleet_code <- data_list2$fleet_control$Fleet_code + fleet_index1

  # Update fleet and spp indices of matrices
  # Skip fleet_control (index 1); re-indexed explicitly above.
  for(i in seq_along(mat_names)[-1]){
    if(!is.null(data_list2[[mat_names[i]]]$Species)){
      data_list2[[mat_names[i]]]$Species <- data_list2[[mat_names[i]]]$Species + nspp1
    }
    if(!is.null(data_list2[[mat_names[i]]]$Fleet_code)){
      data_list2[[mat_names[i]]]$Fleet_code <- data_list2[[mat_names[i]]]$Fleet_code + fleet_index1
    }
  }

  # Update stomach pred sp and prey sp
  data_list2$diet_data$Pred <- data_list2$diet_data$Pred + nspp1
  data_list2$diet_data$Prey <- data_list2$diet_data$Prey + nspp1

  # Update weight and alk indices
  data_list2$weight$Wt_index <- data_list2$weight$Wt_index + weight_index1
  data_list2$age_trans_matrix$Age_transition_index <- data_list2$age_trans_matrix$Age_transition_index + alk_index1



  # Combine vectors
  for(i in vec_names){
    data_list_new[[i]] <- c(data_list1[[i]], data_list2[[i]])
  }

  # Optional per-species vectors: a data set without one uses NA for each of its
  # species (for maturity-at-length, "use the age-based maturity sheet"), so the
  # combined vector stays aligned with the species.
  for(i in c("L50_mat_len", "slope_mat_len")){
    if(!is.null(data_list1[[i]]) || !is.null(data_list2[[i]])){
      data_list_new[[i]] <- c(data_list1[[i]] %||% rep(NA_real_, data_list1$nspp),
                              data_list2[[i]] %||% rep(NA_real_, data_list2$nspp))
    }
  }

  # Combine matrices. An optional one -- equil_catch_data, which only a cleaned
  # data set carries -- may be absent from both, and stays absent.
  for(i in mat_names){
    if(is.null(data_list1[[i]]) && is.null(data_list2[[i]])) next
    data_list_new[[i]] <- plyr::rbind.fill(data_list1[[i]], data_list2[[i]])
  }

  data_list_new$env_data <- merge(data_list1$env_data, data_list2$env_data,
                                  by = "Year", all = TRUE)


  # Add new species
  data_list_new$nspp <- data_list1$nspp + data_list2$nspp

  return(.rce_as_data(data_list_new))
}
