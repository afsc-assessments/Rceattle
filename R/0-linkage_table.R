#' Linkage table: long-format coefficient registry for Rceattle processes
#'
#' A *linkage table* lists every coefficient (`beta`) that ties a process
#' (recruitment, natural mortality, growth, ...) to a covariate, time block or
#' stratum. Each row is one coefficient and records the process and parameter
#' it changes, the species/sex/age/fleet it applies to, the column of the
#' design matrix `X` it multiplies, and its link, prior, bounds and phase.
#'
#' One row per coefficient, rather than one array per process, means the table
#' grows only with the linkages requested. The user-facing constructors are in
#' `R/0-build_linkage.R`.
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
  process       = "character",  # "recruitment", "M", "growth", "q", "sel", "comp"
  param         = "character",  # e.g. "R0", "alpha", "K"
  species       = "integer",    # 1-based species id; NA = every level
  sex           = "integer",    # 1 or 2; NA = every level
  age_bin       = "integer",    # 1-based age index; NA = every level
  fleet         = "integer",    # 1-based Fleet_code; NA = every level
  X_col         = "integer",    # column of the global design matrix
  design_col    = "character",  # name of the design matrix column
  link          = "character",  # "identity", "log", "logit", "exponential"
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
  # Random effects, filled by pool_linkages(); NA on fixed rows.
  re_index      = "integer",    # 0-based deviation slot, across both RE vectors
  sigma_index   = "integer",    # 0-based slot in log_sigma_linkage
  re_time       = "numeric",    # grouping value that orders rw/ar1 deviations in time
  # Per-group SD (sigma), from linkage_spec(init = , priors = list(sigma = ));
  # the same on every row of a group.
  re_sigma_init = "numeric",       # start / fixed SD on natural scale (NA = default)
  re_sigma_prior_family = "character", # prior on the SD; "none"/NA = no prior
  re_sigma_prior_p1     = "numeric",
  re_sigma_prior_p2     = "numeric",
  # Per-group ar1 correlation rho on (-1, 1), from linkage_spec(init = ,
  # priors = list(rho = )); NA on non-ar1 rows.
  re_rho_init   = "numeric",       # start / fixed correlation (NA = default 0)
  re_rho_prior_family = "character",   # prior on rho; "none"/NA = no prior
  re_rho_prior_p1     = "numeric",
  re_rho_prior_p2     = "numeric",
  # Rogers QAR1, from linkage_spec(observe = , obs_sd = , obs_sd_est = ): the
  # latent ar1 deviation is also observed as an env_data column, and acts
  # through an estimated effect size. NA when not observed.
  re_obs_value  = "numeric",       # observed covariate value at this row's time
  re_obs_sd     = "numeric",       # measurement SD start/fixed value (NA = not observed)
  re_obs_est    = "logical",       # TRUE = estimate the obs SD; FALSE/NA = hold it fixed
  # From linkage_spec(integrate = ); the same on every row of a group.
  re_integrate  = "logical",       # TRUE/NA = integrated; FALSE = penalized fixed effect
  # Position within beta_linkage_re or beta_linkage_re_pen (re_index is the
  # slot across both), for setting `inits` by hand.
  re_pos        = "integer"        # 0-based position within its own parameter vector
)


# Reserved sets and their checks ----

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
LINKAGE_LINKS <- c("identity", "log", "logit", "exponential")

#' Link functions with a C++ accumulator behind them
#'
#' The accumulators read link codes 1 (log), 0 (identity) and 3 (exponential).
#' `"logit"` is reserved and refused until one implements it.
#'
#' `"exponential"` multiplies the parameter on its log scale,
#' `exp((log_base + log_offset) * exp(beta * x))`: Stock Synthesis's
#' environmental link type 1 (`SS_timevaryparm.tpl:206`). It applies only where
#' SS3 stores the parameter as a log, so it is wired for catchability under a
#' lognormal/t index error; under a normal index SS3 keeps `Q_parm` arithmetic
#' (`SS_expval.tpl:413-419`) and its type 1 is Rceattle's `"log"`. SS3's
#' `SR_LN(R0)` is a log too, but no recruitment accumulator reads this link yet.
#'
#' @keywords internal
#' @noRd
LINKAGE_LINKS_IMPLEMENTED <- c("identity", "log", "exponential")


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


# Constructor and validator ----

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


# `exponential` is SS3's environmental link type 1, which multiplies a parameter
# on whatever scale SS3 stores it on (SS_timevaryparm.tpl:206-211).
.check_exponential_link <- function(tbl) {
  if (is.null(tbl) || !is.data.frame(tbl) || nrow(tbl) == 0L) {
    return(invisible(tbl))
  }
  if (is.null(tbl[["link"]])) return(invisible(tbl))
  is_exp <- tbl[["link"]] == "exponential"
  if (!any(is_exp)) return(invisible(tbl))

  bad <- which(is_exp & tbl[["process"]] != "q")
  if (length(bad)) {
    stop("link = \"exponential\" is only consumed by catchability (process ",
         "\"q\"), and was given on: ",
         paste(unique(tbl[["process"]][bad]), collapse = ", "), ".\n",
         "  It reproduces Stock Synthesis's environmental link type 1, which ",
         "multiplies the parameter on\n  whatever scale SS3 stores it on, so ",
         "the right substitute depends on the process:\n",
         "  - M, growth: SS3 stores these on the natural scale ",
         "(SS_biofxn.tpl:1074, :265-275), where its\n    type 1 is ",
         "parm * exp(beta * x) -- exactly what link = \"log\" gives. Use ",
         "\"log\".\n",
         "  - recruitment R0 and R_init: SS3 stores both levels as logs ",
         "(SR_LN(R0), and the initial\n    level as an SR_regime block), so ",
         "the form is right there, but no accumulator\n    consumes it yet, ",
         "and \"log\" is SS3's type 2 rather than a substitute.\n",
         "  See vignette(\"environmental-linkages-and-priors\") for the other ",
         "processes.",
         call. = FALSE)
  }

  # Every RE row lands in q_linkage_offset, which this link multiplies, so its
  # deviations would be scored at a constant sigma while varying by covariate.
  q_re <- which(!is.na(tbl[["re_struct"]]) & tbl[["process"]] == "q")
  if (length(q_re)) {
    exp_flt <- tbl[["fleet"]][is_exp]
    re_flt  <- tbl[["fleet"]][q_re]
    any_shared <- anyNA(exp_flt) || anyNA(re_flt)
    both <- intersect(exp_flt, re_flt)
    if (any_shared || length(both)) {
      shared <- if (any_shared) {
        "all fleets (a spec with no `fleet =` is shared)"
      } else {
        paste(sort(unique(both)), collapse = ", ")
      }
      stop("an \"exponential\" catchability linkage and a random-effect ",
           "catchability linkage both target\n  fleet(s) ", shared,
           ". Every log-link row accumulates into q_linkage_offset, which the\n",
           "  exponential link multiplies, so the deviations would be scaled ",
           "by exp(beta * x) while\n  their density still scores them at a ",
           "constant sigma.\n",
           "  Put the environmental effect and the deviations on different ",
           "fleets, or drop one. For both on\n  one fleet, express the ",
           "deviations as Time_varying_q, which the template keeps outside\n",
           "  the multiply, as SS3 does its dev blocks.",
           call. = FALSE)
    }
  }
  invisible(tbl)
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
  # Backstop for a table not built by materialize_linkage(): a reserved process
  # must not reach the template, where nothing would read it.
  for (pr in setdiff(unique(x$process), LINKAGE_PROCESSES_IMPLEMENTED)) {
    .check_process_implemented(pr)
  }
  bad_link <- setdiff(unique(x$link), LINKAGE_LINKS)
  if (length(bad_link) > 0) {
    stop("unknown link(s) in linkage table: ",
         paste(bad_link, collapse = ", "),
         "; allowed: ", paste(LINKAGE_LINKS, collapse = ", "))
  }
  # The same backstop for a link not set through linkage_spec().
  for (lk in setdiff(unique(x$link), LINKAGE_LINKS_IMPLEMENTED)) {
    .check_link_implemented(lk)
  }
  # After the type checks, as it reads `link`, `process` and `re_struct`.
  .check_exponential_link(x)
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


# linkage_row() ----

#' Build a single linkage row
#'
#' Returns a one-row linkage table with defaults for the optional columns, for
#' tests and for assembling a table row by row.
#'
#' @param process,param,X_col required identifying fields.
#' @param species,sex,age_bin,fleet stratum ids, `NA` = every level; `fleet`
#'   is a 1-based `Fleet_code`, used by catchability and selectivity linkages.
#' @param design_col name of the design matrix column.
#' @param link link function; one of [LINKAGE_LINKS].
#' @param init initial value (default `0`).
#' @param init_supplied `TRUE` if the user supplied `init` (default `FALSE`).
#' @param lower,upper bounds (default `-Inf`, `Inf`).
#' @param prior_family one of [PRIOR_FAMILIES]. `"none"` = no prior.
#' @param prior_p1,prior_p2 family-specific prior parameters; ignored when
#'   `prior_family == "none"`.
#' @param re_group random-effect grouping label; `NA` = fixed.
#' @param re_struct random-effect covariance structure (`"us"`/`"rw"`/`"ar1"`);
#'   `NA` = fixed.
#' @param est_phase estimation phase ordinal; `0` = fix at `init`.
#' @param re_index,sigma_index,re_time random-effect fields filled by
#'   `pool_linkages()` (deviation slot, SD slot, time order); `NA` on fixed rows.
#' @param re_sigma_init,re_sigma_prior_family,re_sigma_prior_p1,re_sigma_prior_p2
#'   the group's natural-scale SD start value (held fixed if given without a
#'   prior) and prior, from `linkage_spec(init = , priors = )`; `NA` on fixed rows.
#' @param re_rho_init,re_rho_prior_family,re_rho_prior_p1,re_rho_prior_p2
#'   the `ar1` group's correlation on `(-1, 1)`, start value (held fixed if
#'   given without a prior) and prior, from `linkage_spec()`; `NA` otherwise.
#' @param re_obs_value,re_obs_sd the Rogers QAR1 observed covariate value and
#'   measurement SD, from `linkage_spec(observe = , obs_sd = )`; `NA` when unobserved.
#' @param re_obs_est `TRUE` to estimate the QAR1 measurement SD
#'   (`linkage_spec(obs_sd_est = )`); `FALSE`/`NA` holds it fixed.
#' @param re_integrate `FALSE` for a penalized fixed effect rather than a
#'   Laplace-integrated random effect, from `linkage_spec(integrate = )`;
#'   `TRUE`/`NA` = integrated.
#' @param re_pos 0-based position of the deviation within `beta_linkage_re` or
#'   `beta_linkage_re_pen`, whichever holds it.
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


# Row indices and slot scales ----

#' Resolve a linkage row's stratum ids against per-species dimensions.
#'
#' Turns the row's `species`, `sex`, `age_bin` and `fleet` ids (`NA` = every
#' level) into index vectors, for the parameter, map and bounds builders.
#'
#' @param row a one-row slice of an `Rceattle_linkage_table`.
#' @param data_list the data list (used for `nspp`, `nsex`, `nages`, `fleet_control`).
#' @return a list: `species` (ids), `per_sp` (the `sex` and `age` vectors for
#'   each species, keyed by its id as a string) and `fleet` (fleet rows).
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
  # A fleet already implies its species, so fleets resolve once, not per species.
  flt <- if (is.na(row$fleet)) {
    seq_len(nrow(data_list$fleet_control))
  } else {
    as.integer(row$fleet)
  }

  list(species = spp, per_sp = per_sp, fleet = flt)
}


#' Is a `sel_inf` slot held on the natural scale for this selectivity form?
#'
#' Slot 1 (`inf_asc` / `peak`) is always an age or length. Slot 2 is an
#' inflection for Logistic, DoubleLogistic and DescendingLogistic, but
#' DoubleNormal stores `logit(right_floor)` there and LogisticPM a log age-1
#' selectivity, where a natural-scale value would be wrong (`right_floor = 0.2`
#' would read as 0.55), so the caller refuses those slots.
#'
#' @param slot the parameter's `.SEL_PARAM_TO_SLOT` index (1 or 2).
#' @param form the fleet's `Selectivity` value.
#' @keywords internal
#' @noRd
.sel_inf_is_natural <- function(slot, form) {
  if (identical(as.integer(slot), 1L)) return(TRUE)
  as.character(form) %in% c("Logistic", "DoubleLogistic", "DescendingLogistic")
}


# Shared and mirrored blocks ----

#' Is this fleet a follower in a shared selectivity / catchability block?
#'
#' Fleets with the same `Selectivity_index` / `Catchability_index` estimate ONE
#' parameter block, copied from the donor: the group's first fleet that is not
#' `Off` (the `.group_lead()` rule). A value set on any other member is
#' overwritten. A group of one is not shared.
#'
#' @return `NA_integer_` if `flt` is not a follower, else the donor's
#'   fleet_control row (equal to its `Fleet_code`, which `data_check()` requires).
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

  # Read through the map: before switch_check() Fleet_type may be the code 0.
  off  <- if (is.null(fc$Fleet_type)) rep(FALSE, nrow(fc)) else
    vapply(fc$Fleet_type, function(x) identical(.canon_switch(x, fleet_map), "Off"), logical(1))
  est  <- rows[!off[rows]]
  lead <- if (length(est)) est[1] else rows[1]
  if (identical(as.integer(lead), as.integer(flt))) NA_integer_ else as.integer(lead)
}


#' Every member of the blocks the named fleets belong to
#'
#' A shared block starts at `TMB:::updateMap()`'s mean over all its members,
#' and `fit_mod()` takes one bound per map level by first occurrence, so a start
#' value or bound written to the donor alone would not survive. Groups may span
#' species, so members are not filtered by species; a second-sex write to a
#' one-sex member lands in a mapped-off padding cell (harmless: it is never read).
#'
#' @keywords internal
#' @noRd
.shared_block_members <- function(data_list, flts, process) {
  col <- switch(process, q = "Catchability_index", sel = "Selectivity_index", NULL)
  idx <- if (is.null(col)) NULL else data_list$fleet_control[[col]]
  if (is.null(idx) || !length(flts)) return(flts)
  out <- as.integer(flts)
  for (f in out[!is.na(out)]) {
    if (!is.na(idx[f])) out <- c(out, which(!is.na(idx) & idx == idx[f]))
  }
  sort(unique(as.integer(out)))
}


#' Guards shared by the intercept init and bounds pushes
#'
#' @keywords internal
#' @noRd
.stop_if_shared_block <- function(data_list, flt, process, param) {
  # Only a fleet the user named is checked; a linkage on every fleet sets the
  # whole group, which is unambiguous.
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

#' Refuse a selectivity / catchability linkage a mirrored block cannot honour
#'
#' A mirrored group shares one coefficient but each member gets its own offset,
#' so the members stay identical only if each has exactly one row and the rows
#' agree. Refuses partial coverage, a random effect, a duplicate row, and specs
#' disagreeing on `link`/`bounds`/`init`; pinned intercepts are exempt.
#' Measured numbers in `inst/dev/TRAPS.md`, "Shared parameter blocks".
#'
#' @param process `"sel"` or `"q"`.
#' @keywords internal
#' @noRd
.stop_if_mirrored_block_linkage <- function(linkage_table, fleet_control, process) {
  col <- switch(process, q = "Catchability_index", sel = "Selectivity_index", NULL)
  if (is.null(col) || is.null(linkage_table) || nrow(linkage_table) == 0L) {
    return(invisible())
  }
  idx <- fleet_control[[col]]
  if (is.null(idx)) return(invisible())
  tbl <- linkage_table[linkage_table[["process"]] == process, , drop = FALSE]
  if (nrow(tbl) == 0L) return(invisible())

  what <- if (identical(process, "q")) "catchability" else "selectivity"
  nm   <- fleet_control[["Fleet_name"]]
  # Read through the map, as `.shared_block_lead()` does, so the two agree.
  off  <- vapply(fleet_control[["Fleet_type"]],
                 function(x) identical(.canon_switch(x, fleet_map), "Off"),
                 logical(1))
  # Coverage is owed only by members whose parameter is estimated, as only they
  # can carry a row. q is read raw, as `.canon_switch()` maps an absent switch
  # to "<blank>", not NA.
  linkable <- if (identical(process, "q")) {
    qf <- as.character(fleet_control[["Catchability"]])
    !is.na(qf) & !qf %in% c(.Q_LINKAGE_SELFBUILT_FORMS,
                            .Q_LINKAGE_UNESTIMATED_FORMS)
  } else {
    .canon_switch(fleet_control[["Selectivity"]], sel_map) %in%
      .SEL_LINKAGE_WIRED_FORMS
  }

  for (g in unique(idx[!is.na(idx)])) {
    members <- which(!is.na(idx) & idx == g)
    # The map ties any group of two or more, `Off` members included; only
    # `covers` owes a row.
    if (length(members) < 2L) next
    covers <- members[!off[members] & linkable[members]]

    grp <- tbl[!is.na(tbl[["fleet"]]) & tbl[["fleet"]] %in% members, , drop = FALSE]
    if (nrow(grp) == 0L) next

    # (a) Each named fleet would get its own deviations and SD: separate sigma
    # groups would score the same deviations twice and fit two SDs to them.
    re <- grp[!is.na(grp[["re_struct"]]), , drop = FALSE]
    if (nrow(re) > 0L) {
      stop(sprintf(paste0(
        "random-effect %s linkage `%s` on fleet(s) %s, which share a %s and so ",
        "estimate ONE parameter block: each named fleet would get its own ",
        "deviation series and its own deviation SD, and the fleets would not ",
        "mirror each other. Place the random effect on a parameter the group ",
        "does not share, or give the fleets separate %s values."),
        what, paste(unique(re[["param"]]), collapse = ", "),
        paste(sprintf("'%s'", nm[unique(re[["fleet"]])]), collapse = ", "),
        col, col), call. = FALSE)
    }

    # One offset per design column, so every rule below is per column.
    key <- paste(grp[["param"]], grp[["X_col"]], grp[["design_col"]],
                 grp[["species"]], grp[["sex"]], grp[["age_bin"]], sep = "\r")
    for (k in unique(key)) {
      rows <- grp[key == k, , drop = FALSE]
      # A pinned intercept has no per-fleet offset; `.stop_if_shared_block()`
      # governs where its value may be set.
      if (all(.is_pinned_intercept(rows))) next
      covered <- unique(rows[["fleet"]])

      # (b) The offset is added once per ROW, so a member with two rows moves
      # twice as far as the rest on the same coefficient.
      dup <- covered[vapply(covered, function(f)
        sum(rows[["fleet"]] == f) > 1L, logical(1))]
      if (length(dup) > 0L) {
        stop(sprintf(paste0(
          "%s linkage `%s` has %d rows for design column `%s` on fleet(s) %s, ",
          "which share %s %d and so estimate ONE coefficient. The offset is ",
          "added once per row, so that fleet would move twice as far as the ",
          "rest of the group on the same coefficient. Give the column one row ",
          "per fleet."),
          what, rows[["param"]][1], max(table(rows[["fleet"]])),
          rows[["design_col"]][1],
          paste(sprintf("'%s'", nm[dup]), collapse = ", "), col, g),
          call. = FALSE)
      }

      # (c) The named members get an offset the rest never get.
      missed <- setdiff(covers, covered)
      if (length(missed) > 0L && length(intersect(covered, covers)) > 0L) {
        stop(sprintf(paste0(
          "%s linkage `%s` names fleet(s) %s but not %s, which share %s %d and ",
          "so estimate ONE parameter block. The linkage offset is applied per ",
          "fleet, so the named fleet(s) would move and the rest would not. ",
          "Name every fleet in the group, or pass `by = NULL` so the linkage ",
          "is one row that reaches all of them."),
          what, rows[["param"]][1],
          paste(sprintf("'%s'", nm[intersect(covered, covers)]), collapse = ", "),
          paste(sprintf("'%s'", nm[missed]), collapse = ", "),
          col, g), call. = FALSE)
      }

      # (d) One coefficient, so link, bounds and init must agree: the surviving
      # bound would be the first row's, and two inits would start at their mean.
      # est_phase may differ (a fixed donor holds the group). Defaults are labelled.
      dflt <- function(v, is_default) ifelse(is_default, paste(v, "(default)"), v)
      shown <- list(
        link   = as.character(rows[["link"]]),
        bounds = dflt(sprintf("[%s, %s]", format(rows[["lower"]], trim = TRUE),
                              format(rows[["upper"]], trim = TRUE)),
                      !is.finite(rows[["lower"]]) & !is.finite(rows[["upper"]])),
        init   = dflt(format(rows[["init"]], trim = TRUE),
                      !rows[["init_supplied"]]))
      for (f in names(shown)) {
        v <- unique(shown[[f]])
        if (length(v) > 1L) {
          stop(sprintf(paste0(
            "two %s linkage specs for `%s` on fleet(s) %s, which share %s %d ",
            "and so estimate ONE coefficient, disagree on `%s`: %s. The group ",
            "would take one of them by table position, so give both specs the ",
            "same value."),
            what, rows[["param"]][1],
            paste(sprintf("'%s'", nm[covered]), collapse = ", "), col, g, f,
            paste(v, collapse = " vs ")), call. = FALSE)
        }
      }
    }
  }
  invisible()
}

# Argument validators ----

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


# bind_linkage() and print ----

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
