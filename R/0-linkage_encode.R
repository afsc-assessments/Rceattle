#' Encode an Rceattle linkage table into TMB-friendly inputs
#'
#' The template reads the linkage table as parallel integer and numeric
#' vectors plus a design matrix. Text columns (`process`, `param`, `link`,
#' `prior_family`) become 0-based integer codes, an `NA` stratum id
#' (`species`, `sex`, `age_bin`, `fleet`) becomes `0` ("every level"), and
#' `X_col` becomes 0-based.
#'
#' @keywords internal
#' @name linkage_encode
NULL


# Process codes ----

#' Integer codes for the `process` column
#'
#' Stable across versions: a new process is appended, so no code shifts.
#'
#' @keywords internal
LINKAGE_PROCESS_CODES <- c(
  recruitment = 0L,
  M           = 1L,
  growth      = 2L,
  q           = 3L,
  sel         = 4L,
  comp        = 5L
)

# Row predicates ----

#' Rows of a linkage table that are the recruitment `R_init` parameter
#'
#' `R_init` scales the initial age structure and has no base parameter in
#' `rec_pars`. `.check_r_init_linkage()` and, through `.is_pinned_intercept()`
#' and `.is_level_intercept()`, the builders all read it, so they agree.
#'
#' @param tbl A pooled `linkage_table`.
#' @keywords internal
#' @noRd
.is_r_init_linkage_row <- function(tbl) {
  tbl[["process"]] == "recruitment" & tbl[["param"]] == "R_init"
}


#' Intercept rows whose coefficient is pinned at `NA`
#'
#' The base parameter (`rec_pars`, `log_M1`, `index_log_q`, ...) already sets
#' the level, so estimating the intercept too would be a flat ridge. `R_init`
#' has no base parameter, so its intercept IS the level and stays estimable.
#'
#' @param tbl A pooled `linkage_table`.
#' @keywords internal
#' @noRd
.is_pinned_intercept <- function(tbl) {
  tbl[["design_col"]] == "(Intercept)" & !.is_r_init_linkage_row(tbl)
}


#' The one intercept row that carries a level instead of pinning one
#'
#' The `R_init` intercept, whose `init` and `bounds` are natural-scale
#' multipliers on R0 that the builders log. A slope row is on the link scale
#' and is not included.
#'
#' @param tbl A pooled `linkage_table`.
#' @keywords internal
#' @noRd
.is_level_intercept <- function(tbl) {
  tbl[["design_col"]] == "(Intercept)" & .is_r_init_linkage_row(tbl)
}


#' The log-scale starting coefficient for a natural-scale `R_init` level
#'
#' Used by `build_params()` and by `fit_mod()` (which rewrites a fixed,
#' `est_phase = 0` level over a warm start), so both use the same scale.
#'
#' @param init_val The `init` column value, a natural-scale multiplier on R0.
#' @param supplied The row's `init_supplied` flag.
#' @keywords internal
#' @noRd
.r_init_log_start <- function(init_val, supplied) {
  if (init_val < 0 || (isTRUE(supplied) && init_val == 0)) {
    stop("a recruitment `R_init` linkage takes a natural-scale `init` -- a ",
         "multiplier on R0, so it must be greater than 0. Got ", init_val, ".",
         call. = FALSE)
  }
  # An unset `init` arrives as the table default 0 and means a multiplier of 1.
  if (init_val == 0) return(0)
  log(init_val)
}


# Link, structure, prior and param codes ----

# Process names as a message reads them.
.LINKAGE_PROCESS_LABELS <- c(
  recruitment = "recruitment", M = "natural mortality", growth = "growth",
  q = "catchability", sel = "selectivity", comp = "composition-weight"
)


#' Integer codes for the `link` column
#' @keywords internal
LINKAGE_LINK_CODES <- c(
  identity    = 0L,
  log         = 1L,
  logit       = 2L,
  exponential = 3L   # q ^ exp(beta * x): SS3's environmental link type 1
)


#' Integer codes for the random-effect covariance `re_struct`
#'
#' Stable across versions; the C++ density dispatches on these. Covariance
#' structures: `us` (IID), `rw` (random walk), `ar1` (first-order autoregressive).
#'
#' @keywords internal
LINKAGE_STRUCT_CODES <- c(
  us  = 0L,
  rw  = 1L,
  ar1 = 2L
)


#' Integer codes for the `prior_family` column
#' @keywords internal
LINKAGE_PRIOR_CODES <- c(
  none      = 0L,
  normal    = 1L,
  lognormal = 2L,
  gamma     = 3L,
  beta      = 4L
)


#' Integer codes for the `param` column, per process
#'
#' Each process numbers its own parameters, so a name is looked up within its
#' process. An unknown name is an error.
#'
#' @keywords internal
LINKAGE_PARAM_CODES <- list(
  # 0-3 are year-varying mean-growth offsets. The SD endpoints 4-5 take
  # intercept rows only; keep them >= RCEATTLE_N_GROWTH_PARAMS (linkage.hpp).
  growth      = c(K = 0L, L1 = 1L, Linf = 2L, m = 3L,
                  sd_L1 = 4L, sd_Linf = 5L),
  M           = c(M1 = 0L),
  # R_init: the initial recruitment level, read at year 0 only, no deviate penalty.
  recruitment = c(R0 = 0L, alpha = 1L, beta = 2L, R_init = 3L),
  q           = c(q = 0L),
  # Selectivity codes are parameter slots shared across the forms:
  #   0/1  log_sel_slp asc/desc (log slope, or log sigma for DoubleNormal)
  #   2/3  sel_inf asc/desc     (inflection/peak; 3 is DoubleNormal's logit floor
  #                             and LogisticPM's log first-bin selectivity)
  #   4    sel_coff             (non-parametric, every bin)
  #   5    log_sel_apical       (log-scale per-sex multiplier on the curve)
  #   6-11 DoubleNormalSS3 P1-P6 (SS3 pattern 24), each on SS3's own scale
  # Aliases (sigma_*, peak, right_floor) name the same slot in a form's terms.
  # linkage.hpp defines each primary name as RCEATTLE_SEL_<NAME> and comp as
  # RCEATTLE_COMP_<NAME>; LINKAGE_LINK_CODES as RCEATTLE_LINK_<NAME>, and
  # LINKAGE_STRUCT_CODES as RCEATTLE_RE_IID / _RW / _AR1. test-linkage-encode.R
  # checks the values agree and that each sel / comp code has a consume site.
  sel         = c(slp_asc     = 0L, slp_desc = 1L,
                  inf_asc     = 2L, inf_desc = 3L,
                  coff        = 4L,
                  apical      = 5L,
                  # DoubleNormal aliases
                  sigma_asc   = 0L, sigma_desc = 1L,
                  peak        = 2L, right_floor = 3L,
                  # DoubleNormalSS3 P1-P6 under the SS3 manual's names
                  # (P1 is dn_peak, as peak is taken above).
                  dn_peak     = 6L,
                  top_logit   = 7L,
                  ascend_se   = 8L,
                  descend_se  = 9L,
                  start_logit = 10L,
                  end_logit   = 11L),
  # Dirichlet-multinomial theta = exp(weight), prior only, no year-varying
  # offset: the intercept re-targets comp_weights, caal_weights (per fleet)
  # or diet_comp_weights (per predator).
  comp        = c(theta_comp = 0L, theta_caal = 1L, theta_diet = 2L)
)


#' Sentinel meaning "applies to every level of this stratum"
#'
#' Stored where the user passed `NA` for `species`, `sex`, `age_bin` or
#' `fleet`; the offsets loop over every level, but the intercept-prior loop
#' reads 0 as level 1 (which `.drop_off_fleet_sel_priors()` relies on).
#'
#' @keywords internal
LINKAGE_STRATUM_ALL <- 0L


# encode_linkage_for_tmb() ----

#' Encode a linkage table into TMB-friendly parallel vectors
#'
#' @param table an `Rceattle_linkage_table` (typically the output of
#'   `pool_linkages()`).
#' @param X the global design matrix from `pool_linkages()` (passed
#'   through unchanged so callers can stash it alongside).
#'
#' @return A named list with components:
#' \describe{
#'   \item{n_linkage}{`integer(1)`. Number of coefficients
#'     (`nrow(table)`).}
#'   \item{linkage_process, linkage_param, linkage_species,
#'     linkage_sex, linkage_age_bin, linkage_fleet, linkage_X_col, linkage_link,
#'     linkage_re_index, linkage_prior_family}{Parallel `integer(n_linkage)`
#'     vectors. `linkage_re_index` is `-1` on fixed rows, else the 0-based
#'     deviation across `beta_linkage_re` and `beta_linkage_re_pen`.}
#'   \item{linkage_prior_p1, linkage_prior_p2}{Parallel
#'     `numeric(n_linkage)` vectors.}
#'   \item{linkage_X}{The design matrix as passed in.}
#'   \item{linkage_design_names}{`character` colnames of `linkage_X`.}
#' }
#' @keywords internal
encode_linkage_for_tmb <- function(table, X) {
  if (!is_linkage_table(table)) {
    stop("`table` must be an Rceattle_linkage_table", call. = FALSE)
  }
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("`X` must be a numeric matrix", call. = FALSE)
  }
  n <- nrow(table)
  if (n == 0L) {
    return(.empty_tmb_encoding(X))
  }

  # process: known-good (validated upstream); fail loudly on surprises
  proc_int <- unname(LINKAGE_PROCESS_CODES[table$process])
  if (anyNA(proc_int)) {
    stop("encode_linkage_for_tmb: unknown process(es) in table",
         call. = FALSE)
  }

  # param: encoded against per-process namespace
  param_int <- vapply(seq_len(n), function(i) {
    codes <- LINKAGE_PARAM_CODES[[table$process[i]]]
    code  <- codes[table$param[i]]
    if (is.na(code)) {
      stop(sprintf(
        "encode_linkage_for_tmb: unknown param '%s' for process '%s'",
        table$param[i], table$process[i]), call. = FALSE)
    }
    unname(code)
  }, integer(1))

  link_int  <- unname(LINKAGE_LINK_CODES[table$link])
  prior_int <- unname(LINKAGE_PRIOR_CODES[table$prior_family])

  # Random-effect slot per row: -1 on a fixed row, else the 0-based deviation
  # the template reads in place of beta_linkage(i).
  re_index_int <- ifelse(is.na(table$re_index), -1L, as.integer(table$re_index))
  is_re <- re_index_int >= 0L
  # Per deviation: its sigma group, sigma = exp(log_sigma_linkage(group)).
  re_sigma_int <- integer(0)
  if (any(is_re)) {
    # Each deviation belongs to exactly one row, so it enters the model once
    # and is scored once.
    n_re <- sum(is_re)
    if (!setequal(re_index_int[is_re], 0:(n_re - 1L)) ||
        anyDuplicated(re_index_int[is_re])) {
      stop("encode_linkage_for_tmb: re_index is not a bijection onto ",
           "0:(n_re-1); the RE registry is inconsistent.", call. = FALSE)
    }
    # The fixed-effect prior loop would read such a prior on beta_linkage(i),
    # which is held at 0 on a random-effect row.
    if (any(prior_int[is_re] != 0L)) {
      stop("encode_linkage_for_tmb: random-effect rows must not carry a ",
           "fixed-effect prior (prior on the deviation SD goes through ",
           "`priors = list(sigma = ...)`).", call. = FALSE)
    }
    re_rows <- which(is_re)
    re_sigma_int <- integer(n_re)
    re_sigma_int[re_index_int[re_rows] + 1L] <- as.integer(table$sigma_index[re_rows])
    # Observed covariate per deviation (Rogers QAR1), with a mask: a year missing
    # from env_data keeps its latent deviation but adds no observation term.
    re_obs_value <- numeric(n_re)
    re_obs_mask  <- integer(n_re)
    ov <- table$re_obs_value[re_rows]
    obs_finite <- is.finite(ov)
    ov[!obs_finite] <- 0
    slots <- re_index_int[re_rows] + 1L
    re_obs_value[slots] <- as.numeric(ov)
    re_obs_mask[slots]  <- as.integer(obs_finite)
  } else {
    re_obs_value <- numeric(0)
    re_obs_mask  <- integer(0)
  }

  # Per random-effect group, ordered by sigma_index: structure, priors on its
  # SD and correlation, and its QAR1 observation.
  gt <- .re_group_table(table)
  n_re_group <- if (is.null(gt)) 0L else nrow(gt)
  if (is.null(gt)) {
    re_struct_codes    <- integer(0)
    re_rho_idx         <- integer(0)
    re_sigma_prior_fam <- integer(0)
    re_sigma_prior_p1  <- numeric(0)
    re_sigma_prior_p2  <- numeric(0)
    re_rho_prior_fam   <- integer(0)
    re_rho_prior_p1    <- numeric(0)
    re_rho_prior_p2    <- numeric(0)
    re_obs_idx         <- integer(0)
    re_obs_sd          <- numeric(0)
  } else {
    re_struct_codes    <- unname(LINKAGE_STRUCT_CODES[gt$re_struct])
    # 0-based slot in trans_rho_linkage for an ar1 group, -1 otherwise.
    re_rho_idx <- rep(-1L, nrow(gt))
    is_ar1 <- gt$re_struct == "ar1"
    re_rho_idx[is_ar1] <- seq_len(sum(is_ar1)) - 1L
    re_sigma_prior_fam <- unname(LINKAGE_PRIOR_CODES[gt$prior_family])
    re_sigma_prior_p1  <- as.numeric(gt$prior_p1)
    re_sigma_prior_p2  <- as.numeric(gt$prior_p2)
    re_rho_prior_fam   <- unname(LINKAGE_PRIOR_CODES[gt$rho_prior_family])
    re_rho_prior_p1    <- as.numeric(gt$rho_prior_p1)
    re_rho_prior_p2    <- as.numeric(gt$rho_prior_p2)
    # Observed (QAR1) groups: 0-based slot in beta_linkage_obs and the
    # observation SD; -1 and 0 for an unobserved group.
    re_obs_idx <- rep(-1L, nrow(gt))
    re_obs_idx[gt$observed] <- seq_len(sum(gt$observed)) - 1L
    re_obs_sd <- ifelse(gt$observed, as.numeric(gt$obs_sd), 0)
  }

  # Where each deviation is stored: beta_linkage_re if Laplace-integrated,
  # beta_linkage_re_pen if penalized, and its 0-based position there.
  rt <- .re_slot_routing(table)
  if (is.null(rt)) {
    re_integrate_int <- integer(0)
    re_slot_int      <- integer(0)
  } else {
    re_integrate_int <- as.integer(rt$integrate)
    re_slot_int      <- as.integer(rt$pos)
    # Positions must cover 0:(k-1) once each, or two deviations share a
    # parameter or one is never used.
    for (want in c(TRUE, FALSE)) {
      p <- rt$pos[rt$integrate == want]
      if (length(p) && !setequal(p, seq_along(p) - 1L)) {
        stop("internal error: a random-effect deviation is stored twice, or not ",
             "at all. Please report this with the linkage specification that ",
             "produced it.", call. = FALSE)
      }
    }
  }

  # NA stratum ids => sentinel 0 ("applies to all"); else 1-based
  to_stratum <- function(v) {
    out <- as.integer(v)
    out[is.na(out)] <- LINKAGE_STRATUM_ALL
    out
  }

  list(
    n_linkage             = n,
    linkage_process       = proc_int,
    linkage_param         = param_int,
    linkage_species       = to_stratum(table$species),
    linkage_sex           = to_stratum(table$sex),
    linkage_age_bin       = to_stratum(table$age_bin),
    linkage_fleet         = to_stratum(table$fleet),
    linkage_X_col         = as.integer(table$X_col) - 1L,   # 0-based for TMB
    linkage_link          = link_int,
    linkage_re_index      = re_index_int,
    linkage_re_sigma      = re_sigma_int,
    linkage_re_integrate  = re_integrate_int,
    linkage_re_slot       = re_slot_int,
    n_re_group            = n_re_group,
    linkage_re_struct     = re_struct_codes,
    linkage_re_rho        = re_rho_idx,
    linkage_re_sigma_prior_family = re_sigma_prior_fam,
    linkage_re_sigma_prior_p1     = re_sigma_prior_p1,
    linkage_re_sigma_prior_p2     = re_sigma_prior_p2,
    linkage_re_rho_prior_family   = re_rho_prior_fam,
    linkage_re_rho_prior_p1       = re_rho_prior_p1,
    linkage_re_rho_prior_p2       = re_rho_prior_p2,
    linkage_re_obs        = re_obs_idx,
    linkage_re_obs_sd     = re_obs_sd,
    linkage_re_obs_value  = re_obs_value,
    linkage_re_obs_mask   = re_obs_mask,
    linkage_is_intercept  = as.integer(table$design_col == "(Intercept)"),
    linkage_prior_family  = prior_int,
    linkage_prior_p1      = as.numeric(table$prior_p1),
    linkage_prior_p2      = as.numeric(table$prior_p2),
    linkage_X             = X,
    linkage_design_names  = colnames(X) %||% character(ncol(X))
  )
}


#' @keywords internal
#' @noRd
.empty_tmb_encoding <- function(X) {
  list(
    n_linkage             = 0L,
    linkage_process       = integer(0),
    linkage_param         = integer(0),
    linkage_species       = integer(0),
    linkage_sex           = integer(0),
    linkage_age_bin       = integer(0),
    linkage_fleet         = integer(0),
    linkage_X_col         = integer(0),
    linkage_link          = integer(0),
    linkage_re_index      = integer(0),
    linkage_re_sigma      = integer(0),
    linkage_re_integrate  = integer(0),
    linkage_re_slot       = integer(0),
    n_re_group            = 0L,
    linkage_re_struct     = integer(0),
    linkage_re_rho        = integer(0),
    linkage_re_sigma_prior_family = integer(0),
    linkage_re_sigma_prior_p1     = numeric(0),
    linkage_re_sigma_prior_p2     = numeric(0),
    linkage_re_rho_prior_family   = integer(0),
    linkage_re_rho_prior_p1       = numeric(0),
    linkage_re_rho_prior_p2       = numeric(0),
    linkage_re_obs        = integer(0),
    linkage_re_obs_sd     = numeric(0),
    linkage_re_obs_value  = numeric(0),
    linkage_re_obs_mask   = integer(0),
    linkage_is_intercept  = integer(0),
    linkage_prior_family  = integer(0),
    linkage_prior_p1      = numeric(0),
    linkage_prior_p2      = numeric(0),
    linkage_X             = X,
    linkage_design_names  = colnames(X) %||% character(ncol(X))
  )
}
