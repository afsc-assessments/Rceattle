#' Linkage table: long-format coefficient registry for Rceattle processes
#'
#' A *linkage table* is the unified data structure used to describe every
#' estimated coefficient that connects a process (recruitment, natural
#' mortality, growth, ...) to a covariate or stratum. Each row of the
#' table corresponds to exactly one scalar coefficient (a `beta`) and
#' identifies (a) which process and parameter it modifies, (b) which
#' species/sex/age subset it applies to, (c) which column of the shared
#' design matrix `X` it multiplies, and (d) its prior, bounds, link, and
#' estimation phase.
#'
#' The motivation is that the historical Rceattle pattern of one wide
#' parameter array per process (e.g. `beta_rec_pars[nspp, n_env]`,
#' `M1_beta[nspp, nsex, n_env]`) does not generalize cleanly across
#' processes that need different stratifications. A long-format table
#' grows only with the linkages the user actually requests, regardless
#' of the underlying species/sex/age dimensionality.
#'
#' This file defines the schema and basic helpers; user-facing
#' formula-driven constructors live in `R/0-build_linkage.R`.
#'
#' @keywords internal
#' @name linkage_table
NULL


#' Column schema for a linkage table
#'
#' Named character vector mapping column name -> R type. Used by
#' [new_linkage_table()] and [validate_linkage_table()].
#'
#' @keywords internal
LINKAGE_COLS <- c(
  process       = "character",  # "recruitment", "M", "growth", "q", "sel"
  param         = "character",  # e.g. "R0", "alpha", "K"
  species       = "integer",    # 1-based species id; NA = shared
  sex           = "integer",    # 1 or 2; NA = shared
  age_bin       = "integer",    # 1-based age index; NA = shared
  fleet         = "integer",    # 1-based Fleet_code; NA = shared
  X_col         = "integer",    # column of the global design matrix
  design_col    = "character",  # name of the design matrix column
  link          = "character",  # "identity", "log", "logit"
  init          = "numeric",    # initial value on the linear predictor scale
  init_supplied = "logical",    # TRUE iff user explicitly supplied init
  lower         = "numeric",    # lower bound (-Inf for unbounded)
  upper         = "numeric",    # upper bound ( Inf for unbounded)
  prior_family  = "character",  # "none" or one of PRIOR_FAMILIES
  prior_p1      = "numeric",    # family param 1 (mean/meanlog/shape/shape1)
  prior_p2      = "numeric",    # family param 2 (sd/sdlog/rate/shape2)
  re_group      = "character",  # random-effect group name (NA = fixed)
  re_struct     = "character",  # covariance structure (us/rw/ar1/...; NA = fixed)
  est_phase     = "integer",    # phase ordinal; 0 = fixed
  # Random-effect registry (filled by pool_linkages(); NA on fixed rows).
  re_index      = "integer",    # 0-based slot in beta_linkage_re (NA = fixed row)
  sigma_index   = "integer",    # 0-based slot in log_sigma_linkage (NA = fixed row)
  re_time       = "numeric",    # numeric grouping value, for rw/ar1 time order (NA = fixed)
  # Per-group RE-SD (sigma) routing from linkage_spec(); identical across a
  # group's rows, deduped per sigma_index at encode time. NA on fixed rows.
  re_sigma_init = "numeric",       # start / fixed SD on natural scale (NA = default)
  re_sigma_prior_family = "character", # prior on the SD; "none"/NA = no prior
  re_sigma_prior_p1     = "numeric",
  re_sigma_prior_p2     = "numeric",
  # Per-group ar1 correlation (rho) routing from linkage_spec(); natural (-1,1)
  # scale. NA on non-ar1 rows.
  re_rho_init   = "numeric",       # start / fixed correlation (NA = default 0)
  re_rho_prior_family = "character",   # prior on rho; "none"/NA = no prior
  re_rho_prior_p1     = "numeric",
  re_rho_prior_p2     = "numeric",
  # State-space (Rogers QAR1) observation: this row's latent ar1 deviate is also
  # observed as `re_obs_value` (an env_data column) with fixed SD `re_obs_sd`,
  # and enters the target through an estimated effect size. NA when not observed.
  re_obs_value  = "numeric",       # observed covariate value at this row's time
  re_obs_sd     = "numeric",       # measurement SD start/fixed value (NA = not observed)
  re_obs_est    = "logical",       # TRUE = estimate the obs SD; FALSE/NA = hold it fixed
  # Laplace-integrated vs penalized fixed effect, from linkage_spec(integrate =).
  # Constant within a sigma group (the registry key includes it, so a group cannot
  # straddle both). NA on fixed rows.
  re_integrate  = "logical",       # TRUE/NA = integrated; FALSE = penalized fixed effect
  # Where this row's deviation actually sits. `re_index` is the GLOBAL slot; the
  # deviations are stored in two vectors, so the position within the one that
  # holds it is a different number. Supplied so a caller setting `inits` by hand
  # can write beta_linkage_re[re_pos + 1] / beta_linkage_re_pen[re_pos + 1]
  # without having to reconstruct the split. NA on fixed rows.
  re_pos        = "integer"        # 0-based position within its own parameter vector
)


#' Allowed values for fixed enum-like columns
#' @keywords internal
LINKAGE_PROCESSES <- c("recruitment", "M", "growth", "q", "sel", "comp")

#' Processes with a C++ accumulator behind them
#'
#' `LINKAGE_PROCESSES` is the *reserved* set, shared with
#' `src/TMB/linkage.hpp` and `R/0-linkage_encode.R::LINKAGE_PROCESS_CODES`.
#' Only this subset is wired end-to-end; the rest are rejected up front so a
#' reserved process cannot be estimated without affecting the model.
#'
#' @keywords internal
#' @noRd
LINKAGE_PROCESSES_IMPLEMENTED <- c("recruitment", "M", "growth", "q", "sel", "comp")


#' Error on a reserved-but-unwired process
#'
#' @param process character scalar, already matched against
#'   [LINKAGE_PROCESSES].
#' @return `process` invisibly, on success. Throws otherwise.
#' @keywords internal
#' @noRd
.check_process_implemented <- function(process) {
  if (!process %in% LINKAGE_PROCESSES_IMPLEMENTED) {
    stop(sprintf(
      paste0("linkages on process \"%s\" are reserved but not yet ",
             "implemented: no C++ accumulator consumes them, so the ",
             "coefficients would be estimated without affecting the ",
             "model.\n  Implemented processes: %s"),
      process, paste(LINKAGE_PROCESSES_IMPLEMENTED, collapse = ", ")),
      call. = FALSE)
  }
  invisible(process)
}

#' @rdname LINKAGE_PROCESSES
#' @keywords internal
LINKAGE_LINKS <- c("identity", "log", "logit")

#' Link functions with a C++ accumulator behind them
#'
#' Every accumulator in `src/TMB/linkage.hpp` gates on `linkfn == 1` (log) or
#' `linkfn == 0` (identity). `"logit"` stays reserved, the code is referenced
#' by the C++ header, but is rejected until an accumulator implements it;
#' the processes wired today expose only log-scale parameters.
#'
#' @keywords internal
#' @noRd
LINKAGE_LINKS_IMPLEMENTED <- c("identity", "log")


#' Error on a reserved-but-unimplemented link function
#'
#' @param link character scalar, already matched against [LINKAGE_LINKS].
#' @return `link` invisibly, on success. Throws otherwise.
#' @keywords internal
#' @noRd
.check_link_implemented <- function(link) {
  if (!link %in% LINKAGE_LINKS_IMPLEMENTED) {
    stop(sprintf(
      paste0("link = \"%s\" is reserved but not yet implemented: no C++ ",
             "accumulator consumes it, so the linkage would be estimated ",
             "and contribute nothing to the model.\n  Implemented links: %s"),
      link, paste(LINKAGE_LINKS_IMPLEMENTED, collapse = ", ")),
      call. = FALSE)
  }
  invisible(link)
}


#' Construct an empty linkage table with the canonical schema
#'
#' @return An empty `data.frame` with the columns and types defined in
#'   [LINKAGE_COLS], holding class `c("Rceattle_linkage_table", "data.frame")`.
#' @keywords internal
new_linkage_table <- function() {
  cols <- lapply(LINKAGE_COLS, function(type) {
    switch(type,
      character = character(0),
      integer   = integer(0),
      numeric   = numeric(0),
      logical   = logical(0),
      stop("Unknown column type: ", type)
    )
  })
  df <- as.data.frame(cols, stringsAsFactors = FALSE)
  class(df) <- c("Rceattle_linkage_table", "data.frame")
  df
}


#' Test whether an object is a linkage table
#'
#' @param x object to test.
#' @return logical scalar.
#' @keywords internal
is_linkage_table <- function(x) {
  inherits(x, "Rceattle_linkage_table")
}


#' Validate a linkage table against the schema
#'
#' Checks that all required columns are present with the correct R type,
#' that no required values are `NA`, and that enum-like columns
#' (`process`, `link`) only contain allowed values.
#'
#' @param x object to validate.
#' @return `x` invisibly, on success. Throws an error otherwise.
#' @keywords internal
validate_linkage_table <- function(x) {
  if (!is.data.frame(x)) {
    stop("linkage table must be a data.frame; got ", class(x)[1])
  }
  missing_cols <- setdiff(names(LINKAGE_COLS), names(x))
  if (length(missing_cols) > 0) {
    stop("linkage table missing columns: ",
         paste(missing_cols, collapse = ", "))
  }
  for (col in names(LINKAGE_COLS)) {
    expected <- LINKAGE_COLS[[col]]
    actual <- typeof(x[[col]])
    ok <- switch(expected,
      character = is.character(x[[col]]),
      integer   = is.integer(x[[col]]),
      numeric   = is.numeric(x[[col]]),
      logical   = is.logical(x[[col]])
    )
    if (!ok) {
      stop(sprintf("linkage column '%s' must be %s; got %s",
                   col, expected, actual))
    }
  }
  if (nrow(x) == 0) return(invisible(x))

  required_non_na <- c("process", "param", "X_col", "link", "init", "est_phase")
  for (col in required_non_na) {
    if (anyNA(x[[col]])) {
      stop(sprintf("linkage column '%s' contains NA", col))
    }
  }
  bad_proc <- setdiff(unique(x$process), LINKAGE_PROCESSES)
  if (length(bad_proc) > 0) {
    stop("unknown process(es) in linkage table: ",
         paste(bad_proc, collapse = ", "),
         "; allowed: ", paste(LINKAGE_PROCESSES, collapse = ", "))
  }
  # Backstop for tables assembled without going through
  # materialize_linkage(): a reserved-but-unwired process must not
  # reach TMB as a silently inert (q) or cryptically failing (sel) row.
  for (pr in setdiff(unique(x$process), LINKAGE_PROCESSES_IMPLEMENTED)) {
    .check_process_implemented(pr)
  }
  bad_link <- setdiff(unique(x$link), LINKAGE_LINKS)
  if (length(bad_link) > 0) {
    stop("unknown link(s) in linkage table: ",
         paste(bad_link, collapse = ", "),
         "; allowed: ", paste(LINKAGE_LINKS, collapse = ", "))
  }
  # Backstop for tables assembled without going through linkage_spec():
  # a reserved-but-unimplemented link must not reach TMB silently.
  for (lk in setdiff(unique(x$link), LINKAGE_LINKS_IMPLEMENTED)) {
    .check_link_implemented(lk)
  }
  bad_fam <- setdiff(unique(x$prior_family), PRIOR_FAMILIES)
  if (length(bad_fam) > 0) {
    stop("unknown prior family in linkage table: ",
         paste(bad_fam, collapse = ", "),
         "; allowed: ", paste(PRIOR_FAMILIES, collapse = ", "))
  }
  has_prior <- x$prior_family != "none"
  if (any(has_prior & (is.na(x$prior_p1) | is.na(x$prior_p2)))) {
    stop("rows with prior_family != 'none' must have non-NA ",
         "prior_p1 and prior_p2")
  }
  if (any(x$lower > x$upper)) {
    stop("linkage table has rows where lower > upper")
  }
  invisible(x)
}


#' Build a single linkage row
#'
#' Convenience constructor that returns a one-row linkage table with
#' default values for optional columns. Useful in tests and for
#' incremental table assembly.
#'
#' @param process,param,X_col required identifying fields.
#' @param species,sex,age_bin,fleet stratum ids; `NA` = shared across the
#'   dimension. `fleet` is a 1-based `Fleet_code`, used by catchability and
#'   selectivity linkages; the process-level linkages leave it `NA`.
#' @param design_col name of the design matrix column.
#' @param link link function; one of [LINKAGE_LINKS].
#' @param init initial value (default `0`).
#' @param lower,upper bounds (default `-Inf`, `Inf`).
#' @param prior_family one of [PRIOR_FAMILIES]. `"none"` = no prior.
#' @param prior_p1,prior_p2 family-specific prior parameters; ignored when
#'   `prior_family == "none"`.
#' @param re_group random-effect grouping label; `NA` = fixed.
#' @param re_struct random-effect covariance structure (`"us"`/`"rw"`/`"ar1"`);
#'   `NA` = fixed.
#' @param est_phase estimation phase ordinal; `0` = fix at `init`.
#' @param re_index,sigma_index,re_time random-effect registry fields filled by
#'   `pool_linkages()`; `NA` on fixed rows. `re_index` is the 0-based slot in
#'   `beta_linkage_re`, `sigma_index` the 0-based slot in `log_sigma_linkage`,
#'   and `re_time` the numeric grouping value used to order `rw()`/`ar1()`
#'   deviations in real elapsed time.
#' @param re_sigma_init,re_sigma_prior_family,re_sigma_prior_p1,re_sigma_prior_p2
#'   per-group RE-SD routing from `linkage_spec(init = list(sigma = ), priors =
#'   list(sigma = ))`; identical across a group's rows, `NA` on fixed rows.
#'   `re_sigma_init` is the start (or, when supplied without a prior, fixed) SD
#'   on the natural scale; the prior triple places a prior on that SD.
#' @param re_rho_init,re_rho_prior_family,re_rho_prior_p1,re_rho_prior_p2
#'   per-group `ar1` correlation routing from `linkage_spec(init = list(rho = ),
#'   priors = list(rho = ))`; natural `(-1, 1)` scale, `NA` on non-`ar1` rows.
#' @param re_obs_value,re_obs_sd state-space (Rogers QAR1) observation from
#'   `linkage_spec(observe = , obs_sd = )`: the observed covariate value at this
#'   row's time and the fixed measurement SD. `NA` when the group is unobserved.
#' @param re_obs_est `TRUE` to estimate the QAR1 measurement SD; `FALSE`/`NA`
#'   holds it fixed.
#' @param re_integrate `FALSE` when the deviations are estimated as a penalized
#'   fixed effect rather than integrated out by the Laplace approximation, from
#'   `linkage_spec(integrate = FALSE)`. `TRUE`/`NA` = integrated.
#' @param re_pos 0-based position of this row's deviation within the parameter
#'   vector that holds it (`beta_linkage_re` when `re_integrate`, else
#'   `beta_linkage_re_pen`). Distinct from `re_index`, which is the global slot.
#' @return A one-row `Rceattle_linkage_table`.
#' @keywords internal
linkage_row <- function(process, param, X_col,
                        species       = NA_integer_,
                        sex           = NA_integer_,
                        age_bin       = NA_integer_,
                        fleet         = NA_integer_,
                        design_col    = NA_character_,
                        link          = "identity",
                        init          = 0,
                        init_supplied = FALSE,
                        lower         = -Inf,
                        upper         =  Inf,
                        prior_family  = "none",
                        prior_p1      = NA_real_,
                        prior_p2      = NA_real_,
                        re_group      = NA_character_,
                        re_struct     = NA_character_,
                        est_phase     = 1L,
                        re_index      = NA_integer_,
                        sigma_index   = NA_integer_,
                        re_time       = NA_real_,
                        re_sigma_init = NA_real_,
                        re_sigma_prior_family = NA_character_,
                        re_sigma_prior_p1     = NA_real_,
                        re_sigma_prior_p2     = NA_real_,
                        re_rho_init   = NA_real_,
                        re_rho_prior_family = NA_character_,
                        re_rho_prior_p1     = NA_real_,
                        re_rho_prior_p2     = NA_real_,
                        re_obs_value  = NA_real_,
                        re_obs_sd     = NA_real_,
                        re_obs_est    = NA,
                        re_integrate  = NA,
                        re_pos        = NA_integer_) {
  out <- new_linkage_table()
  out[1L, ] <- list(
    process       = as.character(process),
    param         = as.character(param),
    species       = as.integer(species),
    sex           = as.integer(sex),
    age_bin       = as.integer(age_bin),
    fleet         = as.integer(fleet),
    X_col         = as.integer(X_col),
    design_col    = as.character(design_col),
    link          = as.character(link),
    init          = as.numeric(init),
    init_supplied = as.logical(init_supplied),
    lower         = as.numeric(lower),
    upper         = as.numeric(upper),
    prior_family  = as.character(prior_family),
    prior_p1      = as.numeric(prior_p1),
    prior_p2      = as.numeric(prior_p2),
    re_group      = as.character(re_group),
    re_struct     = as.character(re_struct),
    est_phase     = as.integer(est_phase),
    re_index      = as.integer(re_index),
    sigma_index   = as.integer(sigma_index),
    re_time       = as.numeric(re_time),
    re_sigma_init = as.numeric(re_sigma_init),
    re_sigma_prior_family = as.character(re_sigma_prior_family),
    re_sigma_prior_p1     = as.numeric(re_sigma_prior_p1),
    re_sigma_prior_p2     = as.numeric(re_sigma_prior_p2),
    re_rho_init   = as.numeric(re_rho_init),
    re_rho_prior_family = as.character(re_rho_prior_family),
    re_rho_prior_p1     = as.numeric(re_rho_prior_p1),
    re_rho_prior_p2     = as.numeric(re_rho_prior_p2),
    re_obs_value  = as.numeric(re_obs_value),
    re_obs_sd     = as.numeric(re_obs_sd),
    re_obs_est    = as.logical(re_obs_est),
    re_integrate  = as.logical(re_integrate),
    re_pos        = as.integer(re_pos)
  )
  validate_linkage_table(out)
  out
}


#' Resolve a linkage row's stratum ids against per-species dimensions.
#'
#' A row's `species`, `sex`, and `age_bin` are integer ids or `NA` (which
#' means "shared / broadcast across this dimension"). Both
#' `map_linkage_adjuster()` and the intercept-zeroing pass in
#' `build_params()` need to translate that into concrete index vectors.
#'
#' @param row a one-row slice of an `Rceattle_linkage_table`.
#' @param data_list the data list (used for `nspp`, `nsex`, `nages`).
#' @return a list with components `species`, `sex`, `age`, each a list
#'   keyed by species id, giving the sex/age index vectors to apply for
#'   that species.
#' @keywords internal
#' @noRd
.linkage_row_indices <- function(row, data_list) {
  spp <- if (is.na(row$species)) seq_len(data_list$nspp) else as.integer(row$species)
  per_sp <- vector("list", length(spp))
  names(per_sp) <- as.character(spp)
  for (s in spp) {
    sx <- if (is.na(row$sex))     seq_len(data_list$nsex[s])  else as.integer(row$sex)
    ag <- if (is.na(row$age_bin)) seq_len(data_list$nages[s]) else as.integer(row$age_bin)
    per_sp[[as.character(s)]] <- list(sex = sx, age = ag)
  }
  # Fleet is not nested inside species the way sex and age are: a fleet
  # already implies its species via fleet_control$Species, so it resolves
  # once against the full fleet set rather than per species.
  flt <- if (is.na(row$fleet)) {
    seq_len(nrow(data_list$fleet_control))
  } else {
    as.integer(row$fleet)
  }

  list(species = spp, per_sp = per_sp, fleet = flt)
}


#' Is a `sel_inf` slot held on the natural scale for this selectivity form?
#'
#' `sel_inf` is dual-purpose. Slot 1 (`inf_asc` / `peak`) is always an age or
#' length midpoint. Slot 2 is an inflection for the logistic family, but
#' DoubleNormal reuses it as `logit(right_floor)` and LogisticPM as a log
#' age-1 selectivity override. A value written on the wrong one of those is
#' silently wrong, `right_floor = 0.2` would become `plogis(0.2) = 0.55`,
#' so the transformed slots are refused rather than guessed at.
#'
#' @param param linkage parameter name; `slot` its `.SEL_PARAM_TO_SLOT` index.
#' @param form the fleet's `Selectivity` value.
#' @keywords internal
#' @noRd
.sel_inf_is_natural <- function(slot, form) {
  if (identical(as.integer(slot), 1L)) return(TRUE)
  as.character(form) %in% c("Logistic", "DoubleLogistic", "DescendingLogistic")
}


#' Is this fleet a follower in a shared selectivity / catchability block?
#'
#' `Selectivity_index` / `Catchability_index` are group keys, not fleet codes:
#' fleets holding the same value estimate ONE parameter block, and
#' `adjust_map_shared_params()` copies the group's donor slice over the rest.
#' The donor is the first estimated fleet in the group, the same rule
#' `.group_lead()` applies for `flt_sel_lead` / `flt_q_lead`, and the reason an
#' `Off` fleet (whose slice is all NA) never leads. A value set on the donor is
#' what the whole group uses; one set on any other member is overwritten.
#'
#' A group of one is not shared, whatever its key happens to be, a survey
#' catchability counter runs 1..n_survey and rarely matches the fleet code.
#'
#' @return `NA_integer_` if `flt` is not a follower, otherwise the fleet_control
#'   ROW of the donor whose value would win. `data_check()` requires `Fleet_code`
#'   to equal the row number, so the two coincide, but this is a row index.
#' @keywords internal
#' @noRd
.shared_block_lead <- function(data_list, flt, process) {
  col <- switch(process, q = "Catchability_index", sel = "Selectivity_index", NULL)
  if (is.null(col)) return(NA_integer_)
  fc  <- data_list$fleet_control
  idx <- fc[[col]]
  if (is.null(idx) || is.na(idx[flt])) return(NA_integer_)

  rows <- which(!is.na(idx) & idx == idx[flt])
  if (length(rows) < 2L) return(NA_integer_)

  # Read through the map: a workbook that has not been through switch_check()
  # still holds the integer code, and 0 == "Off" is FALSE.
  off  <- if (is.null(fc$Fleet_type)) rep(FALSE, nrow(fc)) else
    vapply(fc$Fleet_type, function(x) identical(.canon_switch(x, fleet_map), "Off"), logical(1))
  est  <- rows[!off[rows]]
  lead <- if (length(est)) est[1] else rows[1]
  if (identical(as.integer(lead), as.integer(flt))) NA_integer_ else as.integer(lead)
}


#' Guards shared by the intercept init and bounds pushes
#'
#' @keywords internal
#' @noRd
.stop_if_shared_block <- function(data_list, flt, process, param) {
  # Only a fleet the user NAMED is checked. An unstratified linkage expands to
  # every fleet, where setting the whole group is unambiguous and the followers
  # take the donor's value anyway.
  if (length(flt) != 1L || is.na(flt)) return(invisible())
  lead <- .shared_block_lead(data_list, flt, process)
  if (is.na(lead)) return(invisible())
  stop(sprintf(
    paste0("linkage `%s` names fleet %d, which mirrors fleet %d's %s. One ",
           "block is estimated for the group and fleet %d's value is the one ",
           "used, so a value set here would be overwritten. Put it on fleet %d."),
    param, as.integer(flt), lead,
    if (identical(process, "q")) "catchability" else "selectivity", lead, lead),
    call. = FALSE)
  invisible()
}

.stop_unless_positive <- function(value, param, target) {
  if (!(value > 0)) {
    stop(sprintf(
      paste0("linkage `%s` intercept value %g is not > 0, and %s is stored on ",
             "the log scale, so it has no logarithm."),
      param, value, target), call. = FALSE)
  }
  invisible()
}

.stop_unless_natural_sel_inf <- function(data_list, flt, slot, param) {
  form <- data_list$fleet_control$Selectivity[flt]
  if (!.sel_inf_is_natural(slot, form)) {
    stop(sprintf(
      paste0("linkage `%s` on fleet %d (%s selectivity) targets a sel_inf slot ",
             "that is not on the natural scale -- DoubleNormal stores ",
             "logit(right_floor) there and LogisticPM a log age-1 selectivity. ",
             "Setting it from a natural-scale value would silently mis-scale ",
             "it; set that parameter through `inits` instead."),
      param, flt, as.character(form)), call. = FALSE)
  }
  invisible()
}


#' Row-bind one or more linkage tables, preserving the schema
#'
#' Wraps `rbind` and re-applies validation and class. Accepts either a
#' list of tables or several tables passed as separate arguments. Empty
#' inputs are skipped.
#'
#' @param ... linkage tables, or a single list of them.
#' @return An `Rceattle_linkage_table`.
#' @keywords internal
bind_linkage <- function(...) {
  args <- list(...)
  if (length(args) == 1L && is.list(args[[1L]]) && !is.data.frame(args[[1L]])) {
    args <- args[[1L]]
  }
  args <- args[vapply(args, NROW, integer(1)) > 0]
  if (length(args) == 0L) return(new_linkage_table())
  for (i in seq_along(args)) validate_linkage_table(args[[i]])
  out <- do.call(rbind, c(args, list(make.row.names = FALSE)))
  class(out) <- c("Rceattle_linkage_table", "data.frame")
  validate_linkage_table(out)
  out
}


#' @export
print.Rceattle_linkage_table <- function(x, ...) {
  cat(sprintf("<Rceattle linkage table: %d coefficient(s)>\n", nrow(x)))
  if (nrow(x) == 0L) return(invisible(x))
  show <- c("process", "param", "species", "sex", "age_bin", "fleet",
            "design_col", "link", "init", "prior_family", "est_phase")
  print(format(x[, show, drop = FALSE]), row.names = FALSE)
  invisible(x)
}
