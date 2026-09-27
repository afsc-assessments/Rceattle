
# Helpers for tests
# Minimal test-data factories and small utilities used across the suite:
#   - expect_all_true(): custom expectation (below)
#   - make_test_data():  single-species data fixture (further below)
#   - calc_*_nll():      independent R reference likelihoods for cross-checking
#                        the TMB jnll_comp components

# Custom expectation: every element of `object` must be TRUE; NA counts as a
# failure. Used by the TMB map / quantity checks in the selectivity, growth,
# and recruitment integration tests, where a whole sub-array of the parameter
# map must be entirely (un)mapped. Reports how many elements failed so a
# mismatch points straight at the offending sub-array. testthat auto-loads
# this file, so test files call it unqualified (like the other helpers here).
expect_all_true <- function(object) {
  lab <- paste(deparse(substitute(object)), collapse = " ")
  # Guard against an empty selection: `all(logical(0))` is TRUE, so without this
  # an index expression that selects nothing would pass vacuously instead of
  # flagging that the sub-array it was meant to check does not exist.
  testthat::expect(
    length(object) > 0L,
    sprintf("%s is empty: expected a non-empty logical vector.", lab)
  )
  ok  <- !is.na(object) & object
  testthat::expect(
    all(ok),
    sprintf("%s is not all TRUE: %d of %d element(s) are FALSE or NA.",
            lab, sum(!ok), length(ok))
  )
  invisible(object)
}

# Reference NLL of a 1-D AR(1) deviation vector `x` with marginal SD `sd` and
# lag-1 correlation `rho`, evaluated as a full multivariate normal. Used to
# cross-check TMB's AR1() density (e.g. selectivity / recruitment deviations).
calc_nll_ar1_1d <- function(x, sd, rho) {
  n <- length(x)
  Sigma_M <- sqrt(sd^2 / (1 - rho^2))

  # 2. Build the Correlation Matrix for an AR1 process
  times <- 1:n
  distance_matrix <- abs(outer(times, times, "-"))
  Correlation_Matrix <- rho^distance_matrix

  # 3. Scale Correlation to Covariance using the Marginal Variance
  Covariance_Matrix <- (Sigma_M^2) * Correlation_Matrix

  # 4. Calculate Negative Log-Likelihood
  # (mean is a vector of 0s because these are mean-zero deviations)
  nll <- -mvtnorm::dmvnorm(x, mean = rep(0, n), sigma = Covariance_Matrix, log = TRUE)

  return(nll)
}

#' Make single-species test data set
#'
#' @param nyrs Number of hindcast years
#' @param nprojyrs Number of projection years
#' @param nages Number of age bins
#' @param seed Optional RNG seed
#' @param minage Minimum age of the population. Default 1 to match the
#'   historical convention; pass 0 to exercise the age = 0 code paths in
#'   the C++ growth and recruitment routines.
#' @param growth One of `"empirical"` (default), `"vonBertalanffy"`, or
#'   `"Richards"`. Controls whether `caal_data` is populated: empirical
#'   leaves it empty (required by `data_check` which rejects empirical
#'   growth + CAAL); parametric variants populate one CAAL row per
#'   length bin so `data_check` is satisfied for tests that later pass
#'   `growthFun = build_growth(fun = ...)` to `fit_mod()`. The helper
#'   does NOT set the growth model itself -- callers control that via
#'   `growthFun` -- it only aligns the data fixture with the choice.
make_test_data <- function(nyrs = 8, nprojyrs = 10, nages = 5, seed = NULL,
                           minage = 1,
                           growth = c("empirical", "vonBertalanffy",
                                      "Richards")) {
  growth <- match.arg(growth)
  if (!requireNamespace("Rceattle", quietly = TRUE)) {
    stop("Rceattle package required for test helpers")
  }

  simData <- list()

  if (!is.null(seed)) set.seed(seed)
  nspp = 1
  simData$nspp <- nspp
  simData$styr <- 1
  simData$endyr <- nyrs
  simData$projyr <- nyrs + nprojyrs
  simData$spnames <- "Sim"
  simData$nsex = 1
  simData$spawn_month = 0
  simData$nages = nages
  simData$minage = minage
  simData$nlengths = nages
  simData$pop_wt_index = 1
  simData$ssb_wt_index = 1
  simData$alpha_wt_len = 0.0001
  simData$beta_wt_len = 3
  simData$pop_age_transition_index = 1
  simData$sigma_rec = 1
  simData$other_food = 1e6
  simData$estDynamics = 0


  years <- seq(simData$styr, simData$endyr)

  # Fleet control
  # - Two simple fleets (survey + fishery)
  simData$fleet_control <- data.frame(
    Fleet_name = c("Survey", "Fishery"),
    Fleet_code = 1:2,
    Fleet_type = c("Survey", "Fishery"),
    Species = 1,
    Month = 0,
    Selectivity_index = 1:2,
    Selectivity = "Logistic",
    Selectivity_dimension = "Age",
    N_sel_bins = NA,
    Sel_curve_pen1 = NA,
    Sel_curve_pen2 = NA,
    Time_varying_sel = 0,
    Time_varying_sel_sd = 1,
    Bin_first_selected = 1,
    Sel_norm_bin = NA,
    Sel_norm_bin_upper = NA,
    Comp_distribution = "Multinomial",
    Comp_weights = 1,
    CAAL_distribution = 0,
    Observation_units = 1,
    Weight_index = 1,
    Age_transition_index = 1,
    Catchability_index = c(1, NA),
    Catchability = c("Fixed", NA),
    Catchability_init = c(1, NA),
    Catchability_prior_sd = c(0.2, NA),
    Time_varying_q = c(0, NA),
    Time_varying_q_sd = c(1, NA),
    Estimate_index_sd = c(0, NA),
    Index_sd = c(1, NA),
    Estimate_catch_sd = c(NA, 0),
    Catch_sd = c(NA, 1),
    Proj_F_proportion = c(NA, 1),
    CAAL_weights = 1,
    Comp_weights_mcallister = 1
  )

  # Deterministic simple observations for fast tests
  total_biom <- rep(100.0, length(years))
  simData$index_data <- data.frame(
    Fleet_name = "Survey",
    Fleet_code = 1,
    Species = 1,
    Year = years,
    Month = 0,
    Selectivity_block = 1,
    Observation = total_biom,
    Log_sd = 0.1
  )

  simData$catch_data <- data.frame(
    Fleet_name = "Fishery",
    Fleet_code = 2,
    Species = 1,
    Year = years,
    Month = 0,
    Selectivity_block = 1,
    Catch = total_biom * 0.1,
    Log_sd = 0.05
  )

  # Minimal age-comp: one row per fleet so the data_check guard
  # ("Fleet(s) with estimated Selectivity but no comp_data or caal_data")
  # passes regardless of `growth`. Proportions are flat 1/nages; tests
  # that don't fit comps just ignore them, while tests that do see a
  # uniform prior.
  comp_cols <- c("Fleet_name", "Fleet_code", "Species", "Sex",
                 "Age0_Length1", "Year", "Month", "Sample_size",
                 paste0("Comp_", 1:nages))
  simData$comp_data <- data.frame(matrix(NA, nrow = 0, ncol = length(comp_cols)))
  colnames(simData$comp_data) <- comp_cols
  flat_comp <- matrix(1 / nages, nrow = 1, ncol = nages)
  simData$comp_data <- rbind(
    simData$comp_data,
    setNames(cbind(data.frame(Fleet_name = "Survey",  Fleet_code = 1L,
                              Species = 1L, Sex = 0L, Age0_Length1 = 0L,
                              Year = years[1], Month = 0L, Sample_size = 1L,
                              stringsAsFactors = FALSE),
                   as.data.frame(flat_comp)),
             comp_cols),
    setNames(cbind(data.frame(Fleet_name = "Fishery", Fleet_code = 2L,
                              Species = 1L, Sex = 0L, Age0_Length1 = 0L,
                              Year = years[1], Month = 0L, Sample_size = 1L,
                              stringsAsFactors = FALSE),
                   as.data.frame(flat_comp)),
             comp_cols)
  )

  # CAAL: empty when growth == "empirical" (data_check rejects empirical
  # growth + CAAL); populated with one row per length bin when growth
  # is parametric so tests passing `growthFun = build_growth(fun = "...")`
  # to fit_mod() have a CAAL likelihood to evaluate.
  caal_cols <- c("Fleet_name", "Fleet_code", "Species", "Sex", "Year",
                 "Length", "Sample_size", paste0("CAAL_", 1:nages))
  if (growth == "empirical") {
    simData$caal_data <- data.frame(matrix(NA, nrow = 0, ncol = length(caal_cols)))
    colnames(simData$caal_data) <- caal_cols
  } else {
    caal_obs_init <- matrix(1 / nages, nrow = nages, ncol = nages)
    colnames(caal_obs_init) <- paste0("CAAL_", 1:nages)
    simData$caal_data <- cbind(
      data.frame(Fleet_name  = "Survey",
                 Fleet_code  = 1L,
                 Species     = 1L,
                 Sex         = 0L,
                 Year        = years[1],
                 Length      = seq_len(nages),
                 Sample_size = 1L),
      caal_obs_init
    )
  }

  # Empirical selectivity stays empty -- both fleets use "Logistic"
  # selectivity above (parametric, estimated). Tests that want flat /
  # fixed selectivity override `Selectivity` and supply emp_sel rows.
  simData$emp_sel <- data.frame(matrix(NA, nrow = 0, ncol = 5 + nages))
  colnames(simData$emp_sel ) = c("Fleet_name", "Fleet_code", "Species", "Sex", "Year", paste("Comp_", 1:nages))

  # Input N-at-age
  simData$NByageFixed <- data.frame(matrix(NA, nrow = 0, ncol = 4 + nages))
  colnames(simData$NByageFixed ) = c("Species_name ", "Species", "Sex", "Year", paste("Age", 1:nages))

  # Age-transition (Age column spans minage..minage+nages-1)
  age_transition <- as.data.frame(diag(nages))
  colnames(age_transition) <- paste0("Length_", seq_len(nages))
  simData$age_trans_matrix <- cbind(data.frame(Age_transition_name = "Base",
                                               Age_transition_index = 1,
                                               Species = 1,
                                               Sex = 0,
                                               Age = minage:(minage + nages - 1L)),
                                    age_transition)

  # Age-error (True_age spans minage..minage+nages-1)
  age_error <- as.data.frame(diag(nages))
  colnames(age_error) <- paste0("Obs_age", seq_len(nages))
  simData$age_error <- cbind(data.frame(
    Species = 1,
    True_age = minage:(minage + nages - 1L)),
    age_error)

  # Weight-at-age
  WAA <- rep(1, nages)
  WAA_df <- as.data.frame(matrix(WAA, nrow = 1))
  colnames(WAA_df) <- paste0("Age", seq_len(nages))
  simData$weight <- cbind(data.frame(Wt_name = "Base", Wt_index = 1, Species = 1, Sex = 0, Year = 0), WAA_df)

  # Maturity
  MatAA_df <- as.data.frame(matrix(1, nrow = 1, ncol = nages))
  colnames(MatAA_df) <- paste0("Age", seq_len(nages))
  simData$maturity <- cbind(data.frame(Species = 1), MatAA_df)

  # Sex ratio
  sexratio <- as.data.frame(matrix(0.5, nrow = 1, ncol = nages))
  colnames(sexratio) <- paste0("Age", seq_len(nages))
  simData$sex_ratio <- cbind(data.frame(Species = 1), sexratio)

  # Mortality
  mort <- as.data.frame(matrix(0.2, nrow = 1, ncol = nages))
  colnames(mort) <- paste0("Age", seq_len(nages))
  simData$M1_base <- cbind(data.frame(Species = 1, Sex = 0), mort)

  # Bioenergetics
  simData$Ceq = 1
  simData$Cindex = 1
  simData$Pvalue = 1
  simData$fday = 1
  simData$CA = 1
  simData$CB = 1
  simData$Qc = 1
  simData$Tco = 1
  simData$Tcm = 1
  simData$Tcl = 1
  simData$CK1 = 1
  simData$CK4 = 1
  simData$Diet_distribution = 1
  simData$Diet_comp_weights = 1

  # Environmental data
  simData$env_data <- data.frame(
    Year = years,
    Index1 = stats::rnorm(nyrs)
  )

  # Diet information Pyrs (relative foraging rate) ----
  simData$ration_data <- simData$weight |>
    dplyr::select(Species, Sex, Year, contains("Age"))


  # Diet proportion ----
  simData$diet_data <- as.data.frame(matrix(NA, nrow = 0, ncol = 9))
  colnames(simData$diet_data) = c("Pred", "Prey", "Pred_sex", "Prey_sex", "Pred_age", "Prey_age",
                                  "Year", "Sample_size", "Stomach_proportion_by_weight")


  # Clean and normalize switches/defaults, then return. switch_check() fills the
  # fleet_control columns that fit_mod()/rearrange_data() rely on (Sel_start_year,
  # Sel_pen_*, Index_distribution, ...) so the fixture can be passed straight to
  # rearrange_data() in tests, not only through fit_mod(). switch_check() is
  # idempotent, so fit_mod() re-running it on this data is a no-op.
  simData <- Rceattle::clean_data(simData)
  simData <- suppressMessages(Rceattle::switch_check(simData))
  return(simData)
}

compile_tmb_if_needed <- function() {
  compile_script <- file.path("src", "TMB", "compile.R")
  if (file.exists(compile_script)) {
    tryCatch(source(compile_script), error = function(e) stop("TMB compile failed: ", e$message))
  }
}

with_loaded_dll <- function(lib, code) {
  if (!file.exists(lib)) stop("DLL not found: ", lib)
  dll <- dyn.load(lib)
  on.exit(if (!is.null(dll)) dyn.unload(dll[["path"]]), add = TRUE)
  force(code)
}





# Reference NLL of a 2-D separable AR(1) x AR(1) deviation field (age x year),
# the R equivalent of TMB's SEPARABLE(AR1(age), AR1(yr)). Cross-checks the
# 2-D random-effects densities (e.g. time- and age-varying selectivity).
calc_nll_ar1_2d <- function(x_matrix, sigma_innov, rho_a, rho_y) {
  n_age <- nrow(x_matrix)
  n_yr  <- ncol(x_matrix)

  # 1. Calculate the Marginal Standard Deviation
  Sigma_M <- sqrt(sigma_innov^2 / ((1 - rho_a^2) * (1 - rho_y^2)))

  # 2. Build the Age Correlation Matrix (Rows)
  dist_age <- abs(outer(1:n_age, 1:n_age, "-"))
  R_age <- rho_a^dist_age

  # 3. Build the Year Correlation Matrix (Columns)
  dist_yr <- abs(outer(1:n_yr, 1:n_yr, "-"))
  R_yr <- rho_y^dist_yr

  # 4. Combine them using the Kronecker Product
  # TMB's SEPARABLE(AR1(age), AR1(yr)) assumes rows=age, cols=yr.
  Correlation_2D <- kronecker(R_yr, R_age)

  # 5. Scale to Covariance
  Covariance_2D <- (Sigma_M^2) * Correlation_2D

  # 6. Flatten the 2D matrix of deviations into a 1D vector
  x_vec <- as.vector(x_matrix)

  # 7. Calculate Negative Log-Likelihood
  nll <- -mvtnorm::dmvnorm(x_vec, mean = rep(0, length(x_vec)), sigma = Covariance_2D, log = TRUE)

  return(nll)
}


# Reference NLL of multinomial composition data: observed counts `obs_num`
# against expected proportions `hat_prop`. Cross-checks the "Multinomial"
# Comp_distribution branch of jnll_comp.
calc_multinom_nll <- function(obs_num, hat_prop) {
  p <- hat_prop / sum(hat_prop)
  # TMB uses the continuous lgamma instead of factorial: x! = gamma(x+1)
  ll <- lgamma(sum(obs_num) + 1) - sum(lgamma(obs_num + 1)) + sum(obs_num * log(p))
  return(-ll)
}

# Reference NLL of Dirichlet-multinomial composition data: observed counts
# `obs_num` with concentration parameters `alpha`. Cross-checks the
# "Dirichlet-multinomial" Comp_distribution branch of jnll_comp.
calc_dirmultinom_nll <- function(obs_num, alpha) {
  N <- sum(obs_num)
  sum_alpha <- sum(alpha)
  ll <- lgamma(N + 1) - sum(lgamma(obs_num + 1)) +
    lgamma(sum_alpha) - lgamma(N + sum_alpha) +
    sum(lgamma(obs_num + alpha) - lgamma(alpha))
  return(-ll)
}

# Every warning an expression raises, not just the first: build_map() raises
# several and the one under test is rarely the earliest.
collect_warnings <- function(expr) {
  w <- character(0)
  withCallingHandlers(expr, warning = function(cond) {
    w <<- c(w, conditionMessage(cond))
    invokeRestart("muffleWarning")
  })
  w
}
