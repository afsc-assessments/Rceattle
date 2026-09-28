#' Encode an Rceattle linkage table into TMB-friendly inputs
#'
#' The TMB template consumes the pooled linkage table as a set of
#' parallel `IVECTOR` / `VECTOR` inputs plus a dense design matrix.
#' R encodes the string-valued columns (`process`, `param`, `link`,
#' `prior_family`) as 0-based integer codes; converts the
#' `NA` stratum ids (`species`, `sex`, `age_bin`, `fleet`) to a sentinel `0`
#' meaning "applies to all levels" (TMB-side dispatch expands the
#' shared rows over the relevant 1-based levels); and converts the 1-based
#' `X_col` to 0-based.
#'
#' @keywords internal
#' @name linkage_encode
NULL


#' Integer codes for the `process` column
#'
#' Stable across versions. New processes append at the end so existing
#' codes never shift.
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

# Process names as a message reads them.
.LINKAGE_PROCESS_LABELS <- c(
  recruitment = "recruitment", M = "natural mortality", growth = "growth",
  q = "catchability", sel = "selectivity", comp = "composition-weight"
)


#' Integer codes for the `link` column
#' @keywords internal
LINKAGE_LINK_CODES <- c(
  identity = 0L,
  log      = 1L,
  logit    = 2L
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
#' Each process has its own independent parameter namespace (`K`
#' means something different to growth than `alpha` does to
#' recruitment), so the param column is encoded against a per-process
#' table. Unknown parameter names error.
#'
#' @keywords internal
LINKAGE_PARAM_CODES <- list(
  # Codes 0..3 index `growth_linkage_offset` (year-varying mean-growth
  # offsets); codes 4..5 are the SD endpoints which live on
  # `growth_log_sd` and only consume intercept rows. Keep the SD codes
  # >= RCEATTLE_N_GROWTH_PARAMS so the env-offset loop in linkage.hpp
  # short-circuits on them (see the param-range guard there).
  growth      = c(K = 0L, L1 = 1L, Linf = 2L, m = 3L,
                  sd_L1 = 4L, sd_Linf = 5L),
  M           = c(M1 = 0L),
  recruitment = c(R0 = 0L, alpha = 1L, beta = 2L),
  q           = c(q = 0L),
  # Selectivity codes index the underlying parameter slots, which are shared
  # across the parametric forms:
  #   0/1 = log_sel_slp[asc/desc]  (log-scale: slope, or log-sigma in DoubleNormal)
  #   2/3 = sel_inf[asc/desc]      (natural: inflection, or peak/logit-floor)
  #   4   = sel_coff               (non-parametric per-bin coefficients)
  #   5   = log_sel_apical         (log-scale per-sex multiplier on the curve)
  # The form-specific aliases (sigma_*, peak, right_floor) resolve to the same
  # slot, so a user names the quantity their form actually has.
  sel         = c(slp_asc     = 0L, slp_desc = 1L,
                  inf_asc     = 2L, inf_desc = 3L,
                  coff        = 4L,
                  apical      = 5L,
                  # DoubleNormal aliases
                  sigma_asc   = 0L, sigma_desc = 1L,
                  peak        = 2L, right_floor = 3L),
  # Dirichlet-multinomial composition-weighting overdispersion. Prior-only
  # (no year-varying accumulator): the intercept re-targets the log DM scalar
  # (comp_weights / caal_weights per fleet, diet_comp_weights per predator;
  # theta = exp()), so a prior on the natural theta falls out of the shared
  # intercept prior loop. theta_comp = age/length comps, theta_caal = CAAL,
  # theta_diet = diet composition.
  comp        = c(theta_comp = 0L, theta_caal = 1L, theta_diet = 2L)
)


#' Sentinel meaning "applies to every level of this stratum"
#'
#' Stored where the user passed `NA` for `species`, `sex`, or
#' `age_bin`. The TMB template expands these by iterating over
#' `1:nspp` / `1:nsex` / `1:nages` rather than indexing a single cell.
#'
#' @keywords internal
LINKAGE_STRATUM_ALL <- 0L


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
#'     `beta_linkage_re` slot.}
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

  # Random-effect coefficient index (per row): -1 on fixed rows, else the
  # 0-based slot in `beta_linkage_re`. The C++ accumulators read
  # `beta_linkage_re(re_index)` for these rows instead of `beta_linkage(i)`.
  re_index_int <- ifelse(is.na(table$re_index), -1L, as.integer(table$re_index))
  is_re <- re_index_int >= 0L
  # Per-slot sigma-group index (length n_re): slot g's log_sigma_linkage group.
  # The C++ density reads sigma = exp(log_sigma_linkage(linkage_re_sigma(g))).
  re_sigma_int <- integer(0)
  if (any(is_re)) {
    # Bijection: every beta_linkage_re slot is referenced by exactly one row,
    # so each has one data path (via the accumulators) and one density term.
    n_re <- sum(is_re)
    if (!setequal(re_index_int[is_re], 0:(n_re - 1L)) ||
        anyDuplicated(re_index_int[is_re])) {
      stop("encode_linkage_for_tmb: re_index is not a bijection onto ",
           "0:(n_re-1); the RE registry is inconsistent.", call. = FALSE)
    }
    # A stray prior on an RE indicator column would leak into the fixed-beta
    # prior loop (which reads beta_linkage(i), frozen at 0 for RE rows).
    if (any(prior_int[is_re] != 0L)) {
      stop("encode_linkage_for_tmb: random-effect rows must not carry a ",
           "fixed-effect prior (prior on the deviation SD goes through ",
           "`priors = list(sigma = ...)`).", call. = FALSE)
    }
    re_rows <- which(is_re)
    re_sigma_int <- integer(n_re)
    re_sigma_int[re_index_int[re_rows] + 1L] <- as.integer(table$sigma_index[re_rows])
    # Per-slot observed covariate value (Rogers QAR1) + a mask flagging which
    # slots actually carry an observation. A slot is un-observed when its year is
    # absent from env_data (`re_obs_value` NA there): the latent still exists (the
    # AR1 runs over all years) but no observation term is added. Store 0 for the
    # value (never read when masked) and 0 in the mask.
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

  # Per-group sigma-prior triple (length n_re_group), from the reserved `sigma`
  # key in linkage_spec(priors = ). Groups are ordered by sigma_index.
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
    # Per group: its 0-based slot in trans_rho_linkage (ar1 groups only, in
    # group order), -1 for non-ar1 groups.
    re_rho_idx <- rep(-1L, nrow(gt))
    is_ar1 <- gt$re_struct == "ar1"
    re_rho_idx[is_ar1] <- seq_len(sum(is_ar1)) - 1L
    re_sigma_prior_fam <- unname(LINKAGE_PRIOR_CODES[gt$prior_family])
    re_sigma_prior_p1  <- as.numeric(gt$prior_p1)
    re_sigma_prior_p2  <- as.numeric(gt$prior_p2)
    re_rho_prior_fam   <- unname(LINKAGE_PRIOR_CODES[gt$rho_prior_family])
    re_rho_prior_p1    <- as.numeric(gt$rho_prior_p1)
    re_rho_prior_p2    <- as.numeric(gt$rho_prior_p2)
    # Per group (Rogers QAR1): 0-based slot in beta_linkage_obs (observed groups,
    # in group order), -1 otherwise; and the fixed observation SD (0 otherwise).
    re_obs_idx <- rep(-1L, nrow(gt))
    re_obs_idx[gt$observed] <- seq_len(sum(gt$observed)) - 1L
    re_obs_sd <- ifelse(gt$observed, as.numeric(gt$obs_sd), 0)
  }

  # Per-slot destination routing: whether a slot is Laplace-integrated (and so
  # stored in beta_linkage_re) or penalized (beta_linkage_re_pen), plus its dense
  # 0-based position within that vector. Same shape as linkage_re_rho above -- a
  # dense sub-index over a subset. For every model without a penalized spec these
  # are all 1 and 0:(n_re - 1), i.e. the identity, so the C++ read is unchanged.
  rt <- .re_slot_routing(table)
  if (is.null(rt)) {
    re_integrate_int <- integer(0)
    re_slot_int      <- integer(0)
  } else {
    re_integrate_int <- as.integer(rt$integrate)
    re_slot_int      <- as.integer(rt$pos)
    # Each destination's positions must be a bijection onto 0:(k-1); otherwise a
    # deviation would share a parameter with another, or leave one unreferenced.
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
