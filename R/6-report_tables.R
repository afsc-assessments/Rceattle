# Assessment reporting tables.
#
# The sections follow the AFSC Alaska Groundfish Stock Assessment Guidelines
# (2024) for what a SAFE chapter reports: the executive-summary quantities
# (section 4.2.2), the fit and likelihood decomposition, and the required
# diagnostics -- retrospective (4.9.7), likelihood profiles (4.9.6), jitter, and
# residuals (one-step-ahead and/or Pearson). The standard harvest scenarios of
# section 4.11.3 need a projection module Rceattle does not have yet, so there
# is no `scenarios` section; it is the one required element this cannot fill.


# Year -> era, in the vocabulary the NOAA standardized output uses: "time" is
# the hindcast a model was fit to, "fore" the projection. Kept as one helper so
# every section splits the two the same way.
.rce_era <- function(year, endyr) {
  ifelse(year > endyr, "fore", "time")
}


# Accept a single fit or a list of them, and give every model a name. Unnamed
# models are numbered rather than left blank, so the `model` column is always a
# usable key for a join or a facet.
.rce_named_models <- function(object, model_names = NULL) {
  models <- if (inherits(object, "Rceattle")) list(object) else object
  if (!is.list(models) || !length(models)) {
    stop("`object` must be an Rceattle fit or a non-empty list of them.",
         call. = FALSE)
  }
  bad <- !vapply(models, inherits, logical(1), "Rceattle")
  if (any(bad)) {
    stop("Every element of `object` must be an Rceattle fit; element(s) ",
         paste(which(bad), collapse = ", "), " are not.", call. = FALSE)
  }
  nm <- model_names %||% names(models)
  if (is.null(nm) || !all(nzchar(nm))) nm <- paste0("Model_", seq_along(models))
  if (length(nm) != length(models)) {
    stop("`model_names` must have one name per model (", length(models), ").",
         call. = FALSE)
  }
  names(models) <- nm
  models
}


# Line up a diagnostic argument with the models. A single object is used for a
# single model; a list must be one entry per model. NULL means the section is
# simply absent -- report_tables() never runs a diagnostic itself, because a
# retrospective or a jitter is tens to hundreds of refits and that is not
# something a table-building call should spend on the user's behalf.
.rce_align_diag <- function(x, models, what, class) {
  if (is.null(x)) return(vector("list", length(models)))
  if (inherits(x, class)) x <- list(x)
  if (!is.list(x)) {
    stop("`", what, "` must be a ", class, " object, a list of them, or NULL.",
         call. = FALSE)
  }
  if (length(x) != length(models)) {
    stop("`", what, "` must be one ", class, " object per model (",
         length(models), "), or NULL.", call. = FALSE)
  }

  nm <- names(models)
  if (!is.null(names(x)) && all(nzchar(names(x)))) {
    # Named: match on the name, so the pairing cannot depend on order.
    unknown <- setdiff(names(x), nm)
    if (length(unknown)) {
      stop("`", what, "` is named for model(s) not in `object`: ",
           paste(unknown, collapse = ", "), ". Models are: ",
           paste(nm, collapse = ", "), call. = FALSE)
    }
    missing <- setdiff(nm, names(x))
    if (length(missing)) {
      stop("`", what, "` has no entry for model(s): ",
           paste(missing, collapse = ", "), call. = FALSE)
    }
    return(x[nm])
  }

  if (length(x) > 1) {
    # Unnamed and ambiguous: one model's diagnostic stored as a list of parts
    # (an osa_residuals() result split by source) would otherwise be attributed
    # one part per model. State the pairing rather than assume it was meant.
    message("`", what, "` is unnamed, so it is paired with the models in order: ",
            paste(sprintf("[%d] -> %s", seq_along(nm), nm), collapse = ", "),
            ". Name the list to be explicit.")
  }
  x
}


#' Assessment reporting tables from one or more Rceattle fits
#'
#' @description Collects the quantities a stock assessment reports into one set
#' of tidy tables, so a SAFE chapter or a model comparison is built from a
#' single call rather than from a dozen ad-hoc extractions. Every table holds
#' a `model` column, so passing several fits gives a like-for-like comparison.
#'
#' @details
#' The sections follow the AFSC Alaska Groundfish Stock Assessment Guidelines
#' for what a chapter reports:
#' \describe{
#'   \item{`model`}{One row per fit: dimensions, switches, the marginal and
#'     joint negative log-likelihoods, the number of estimated parameters, AIC,
#'     the maximum gradient, whether the Hessian was positive definite, and the
#'     run time.}
#'   \item{`parameters`}{Every estimated parameter with its standard error, and
#'     the natural-scale name and process from [parameter_dictionary()]. Where
#'     `sigma_R` and an estimated M are found. Estimates are on the parameter's
#'     own scale, so a `log_` name needs `exp()`; a **fixed** M is not here at
#'     all, because it was never estimated, read it off `M_at_age`.}
#'   \item{`likelihood`}{The negative log-likelihood by component and fleet or
#'     species, weighted and unweighted.}
#'   \item{`timeseries`}{Biomass, female spawning-stock biomass, recruitment,
#'     depletion and fishing mortality by species and year, with standard errors
#'     and confidence intervals, split into hindcast (`era = "time"`) and
#'     projection (`era = "fore"`).}
#'   \item{`reference_points`}{The executive-summary quantities: the SPR-based
#'     F proxies, unfished and target female spawning-stock biomass, the biomass
#'     proxies implied by `Ptarget` / `Plimit`, and terminal status. A `basis`
#'     column says whether each was estimated, and if not, why, see below.}
#'   \item{`fits`}{Observed against predicted index and catch, with the standard
#'     deviation of normalized residuals (SDNR) per fleet.}
#'   \item{`retrospective`, `jitter`, `osa`}{Present only when the corresponding
#'     object is supplied.}
#' }
#'
#' Nothing here refits. The three diagnostics are tens to hundreds of
#' optimizations each, so they are supplied as already-computed objects rather
#' than run inside a table-building call; a section whose object is `NULL` is
#' simply absent from the result.
#'
#' The standard harvest scenarios of guideline section 4.11.3 are **not**
#' produced, they need a standard projection module, which Rceattle does not
#' have. Projected biomass under the model's own harvest control rule is in
#' `timeseries` with `era = "fore"`.
#'
#' # Reference points a fit does not define
#'
#' A CEATTLE fit leaves a *number* in the reported array for several reference
#' points it never estimated, so reading the array directly puts a plausible
#' wrong figure into the one table that becomes the executive summary. Each is
#' returned as `NA` with the reason in `basis`:
#' \itemize{
#'   \item `Ftarget` / `Flimit` are estimated only under a harvest control rule
#'     that defines them, and are switched off for a species with no projected
#'     fishery or with fixed numbers-at-age. Unestimated, they sit at their
#'     initial value of `exp(0) = 1`, which reads as an F of 1.0/yr.
#'   \item `SB0` / `B0` under `msmMode > 0` come from the `MSSB0` / `MSB0`
#'     inputs, which stand at a placeholder until [fit_mod()] derives them from
#'     a no-fishing projection. `B_target` and `B_limit` are fractions of `SB0`
#'     and go with it.
#'   \item The per-recruit quantities (`SPR0`, `SPRtarget`, `SPRlimit`) are
#'     computed only under `msmMode = 0`.
#' }
#' The depletions are **not** blanked alongside `SB0`: under a no-fishing
#' harvest control rule in multispecies mode the model divides by biomass in the
#' last projection year, which is the equilibrated unfished reference, so the
#' series is meaningful there.
#'
#' @section Two negative log-likelihoods:
#' `model` reports both, and **both are minimized**: a smaller value is the
#' better fit. `marginal_nll` is the negative log marginal likelihood the
#' optimizer minimized, random effects integrated out by the Laplace
#' approximation, and is what `AIC` is built from. `joint_nll` is what the
#' template evaluated at the conditional modes, so it is what `likelihood` sums
#' to, on the same scale as `jnll_comp`. They are equal when `n_random` is 0 and
#' differ by the Laplace correction otherwise.
#'
#' @section Supplying diagnostics for several models:
#' A diagnostics list is matched to models **by name**, so `list(alt = ..., base
#' = ...)` pairs correctly whatever the order. An unnamed list is paired
#' positionally and says so in a message. Names that are not model names are an
#' error, which catches the realistic mistake of passing one model's
#' [osa_residuals()] result stored as a list of parts.
#'
#' @param object An Rceattle fit from [fit_mod()], or a list of them.
#' @param model_names Names for the models, defaulting to the names of `object`
#'   or to `Model_1`, `Model_2`, ...
#' @param retro A [retrospective()] result per model, or `NULL` (default) to
#'   omit that section.
#' @param jitter A [jitter()] result per model, or `NULL` (default) to omit that
#'   section.
#' @param osa An [osa_residuals()] result per model, summarized here with
#'   [osa_diagnostics()], or `NULL` (default) to omit that section.
#' @param ci_level Confidence level for the `timeseries` interval, default
#'   `0.95`.
#' @param quantities Time-series quantities to include, from the set
#'   [as.data.frame.Rceattle()] accepts; [quantity_dictionary()] says what each
#'   one means.
#'
#' @return A list of data frames with class `"rceattle_report"`, one element per
#'   section described above. Each holds a `model` column.
#'
#' @seealso [quantity_dictionary()] for what each quantity means and its units,
#'   [as.data.frame.Rceattle()] for the time series alone, and
#'   [standard_output()] to relabel the result into the NOAA standardized
#'   assessment output that `stockplotr` and `asar` consume.
#'
#' @examples
#' \dontrun{
#' fit <- fit_mod(data_list = BS2017SS, estimateMode = 1)
#' tabs <- report_tables(fit)
#' tabs$reference_points
#'
#' # Compare two models, with diagnostics attached
#' tabs <- report_tables(list(base = fit0, alt = fit1),
#'                       retro = list(retro0, retro1))
#' }
#' @export
report_tables <- function(object,
                          model_names = NULL,
                          retro = NULL,
                          jitter = NULL,
                          osa = NULL,
                          ci_level = 0.95,
                          quantities = c("biomass", "ssb", "R",
                                         "biomass_depletion", "ssb_depletion",
                                         "F_spp")) {
  models <- .rce_named_models(object, model_names)
  retro  <- .rce_align_diag(retro,  models, "retro",  "Rceattle_retro")
  jitter <- .rce_align_diag(jitter, models, "jitter", "Rceattle_jitter")
  osa    <- .rce_align_diag(osa,    models, "osa",    "rceattle_osa")

  bind <- function(lst) {
    lst <- lst[!vapply(lst, is.null, logical(1))]
    lst <- lst[vapply(lst, nrow, 1L) > 0]
    if (!length(lst)) return(NULL)
    do.call(rbind, c(unname(lst), list(make.row.names = FALSE)))
  }

  nm <- names(models)
  out <- list(
    model            = bind(lapply(nm, function(m) .rce_tab_model(models[[m]], m))),
    parameters       = bind(lapply(nm, function(m) .rce_tab_parameters(models[[m]], m))),
    likelihood       = bind(lapply(nm, function(m) .rce_tab_likelihood(models[[m]], m))),
    timeseries       = bind(lapply(nm, function(m)
                         .rce_tab_timeseries(models[[m]], m, quantities, ci_level))),
    reference_points = bind(lapply(nm, function(m) .rce_tab_refpoints(models[[m]], m))),
    fits             = bind(lapply(nm, function(m) .rce_tab_fits(models[[m]], m))),
    retrospective    = bind(lapply(seq_along(nm), function(i)
                         .rce_tab_retro(retro[[i]], nm[i]))),
    jitter           = bind(lapply(seq_along(nm), function(i)
                         .rce_tab_jitter(jitter[[i]], nm[i]))),
    osa              = bind(lapply(seq_along(nm), function(i)
                         .rce_tab_osa(osa[[i]], nm[i])))
  )

  out <- out[!vapply(out, is.null, logical(1))]
  structure(out, class = c("rceattle_report", "list"))
}


# -- section builders --------------------------------------------------------

# One row per fit. Everything is guarded: a model built with estimateMode = 3
# carries no $opt, $sdrep, $convergence or $.conv_hindcast at all, and a fit run
# with fit_control(getsd = FALSE) has no sdreport, so each of these is NA rather
# than an error.
.rce_tab_model <- function(fit, model) {
  d <- fit$data_list
  npar <- length(fit$opt$par %||% numeric(0))
  obj  <- fit$opt$objective %||% NA_real_
  aic  <- if (!is.null(fit$opt$objective) && npar > 0)
    2 * npar + 2 * as.numeric(fit$opt$objective) else NA_real_

  # Convergence is judged on the hindcast fit, so read the `.conv_hindcast`
  # snapshot first: the projection re-optimizes and overwrites `opt` and `sdrep`.
  ch <- fit$.conv_hindcast
  mg <- ch$max_gradient %||% fit$opt$max_gradient %||% NA_real_
  pd <- if (!is.null(ch$pdHess)) ch$pdHess
        else if (!is.null(fit$sdrep)) isTRUE(fit$sdrep$pdHess) else NA

  secs <- tryCatch(as.numeric(fit$run_time, units = "secs"),
                   error = function(e) NA_real_)

  # Marginal (what nlminb minimized, and what AIC uses) and joint (what
  # jnll_comp decomposes). They differ by the Laplace correction whenever there
  # are random effects; see ?report_tables.
  n_random <- length(fit$obj$env$random %||% integer(0))
  joint <- fit$quantities$jnll %||% NA_real_

  data.frame(
    model         = model,
    rceattle      = tryCatch(as.character(utils::packageVersion("Rceattle")),
                             error = function(e) NA_character_),
    nspp          = d$nspp %||% NA_integer_,
    species       = paste(d$spnames, collapse = ", "),
    styr          = d$styr %||% NA_integer_,
    endyr         = d$endyr %||% NA_integer_,
    projyr        = d$projyr %||% NA_integer_,
    msmMode       = unname(.rce_alias_show(d$msmMode, msmMode_map)),
    estimateMode  = unname(tryCatch(.rce_alias_show(d$estimateMode, estimateMode_map),
                                    error = function(e) NA_character_)),
    HCR           = unname(tryCatch(.rce_alias_show(d$HCR, hcr_map),
                                    error = function(e) NA_character_)),
    n_parameters  = npar,
    n_random      = as.integer(n_random),
    marginal_nll  = as.numeric(obj),
    joint_nll     = as.numeric(joint),
    AIC           = as.numeric(aic),
    max_gradient  = as.numeric(mg),
    pdHess        = as.logical(pd),
    converged     = fit$convergence$status %||% NA_character_,
    run_time_secs = as.numeric(secs),
    stringsAsFactors = FALSE
  )
}


# The estimated parameters and their standard errors, named through the
# parameter dictionary. Where sigma_R and an estimated M are found; a FIXED M
# was never estimated, so read that off `M_at_age`.
.rce_tab_parameters <- function(fit, model) {
  est <- tryCatch(stats::coef(fit), error = function(e) NULL)
  if (is.null(est) || !length(est)) return(NULL)

  # Fixed effects only. NA when getsd = FALSE, and NA for a negative variance
  # too: a SAFE table must not carry sqrt(|variance|) as an uncertainty.
  se <- tryCatch({
    v <- stats::vcov(fit)
    if (is.null(v)) NULL else .conv_se_from_cov(v)
  }, error = function(e) NULL)
  se <- if (!is.null(se) && length(se) == length(est)) unname(se)
        else rep(NA_real_, length(est))

  # coef() repeats a name per element of a block, so number within the block.
  nm <- names(est)
  idx <- stats::ave(seq_along(nm), nm, FUN = seq_along)
  n_in_block <- as.integer(table(nm)[nm])

  info <- .par_info(unique(nm))
  m <- match(nm, info$internal)

  data.frame(
    model      = model,
    parameter  = nm,
    index      = ifelse(n_in_block > 1L, idx, NA_integer_),
    natural    = info$natural[m],
    process    = info$process[m],
    estimate   = unname(est),
    std_error  = se,
    stringsAsFactors = FALSE
  )
}


# The likelihood decomposition. jnll_comp's columns count FLEETS on rows 1-8,
# SPECIES on rows 9-20 and neither on row 21, so the column key comes from
# .JNLL_ROW_AXIS -- reading every column as a fleet would report a species'
# recruitment penalty against a survey.
.rce_tab_likelihood <- function(fit, model) {
  j <- fit$quantities$jnll_comp
  if (is.null(j)) return(NULL)
  u <- fit$quantities$unweighted_jnll_comp
  d <- fit$data_list

  rn <- rownames(j) %||% paste("Row", seq_len(nrow(j)))
  axis <- .JNLL_ROW_AXIS[rn]
  fleet_names <- d$fleet_control$Fleet_name
  sp_names    <- d$spnames

  grid <- expand.grid(row = seq_len(nrow(j)), col = seq_len(ncol(j)),
                      KEEP.OUT.ATTRS = FALSE)
  # `.JNLL_ROW_AXIS` has three values: a row's columns count fleets, or species,
  # or -- for the linkage random effects -- neither, because that term belongs to
  # the model as a whole. Only the first two name a column; "model" and an
  # unrecognised row both get NA rather than a guessed key.
  ax <- unname(axis[grid$row])
  name_for <- function(a, nms) {
    ok <- !is.na(ax) & ax == a & grid$col <= length(nms)
    out <- rep(NA_character_, length(ax))
    out[ok] <- nms[grid$col[ok]]
    out
  }
  key <- name_for("fleet", fleet_names)
  sp_key <- name_for("species", sp_names)
  key[is.na(key)] <- sp_key[is.na(key)]

  res <- data.frame(
    model      = model,
    component  = rn[grid$row],
    axis       = ax,
    unit       = key,
    weighted   = as.numeric(j[cbind(grid$row, grid$col)]),
    unweighted = if (is.null(u)) NA_real_ else
      as.numeric(u[cbind(grid$row, grid$col)]),
    stringsAsFactors = FALSE
  )
  # A cell that scored nothing is structurally zero, not a small number, so
  # carrying every empty cell would triple the table for no information.
  res <- res[is.finite(res$weighted) & res$weighted != 0, , drop = FALSE]

  total <- data.frame(model = model, component = "Total", axis = NA_character_,
                      unit = NA_character_,
                      weighted = sum(j, na.rm = TRUE),
                      unweighted = if (is.null(u)) NA_real_ else sum(u, na.rm = TRUE),
                      stringsAsFactors = FALSE)
  rbind(res, total, make.row.names = FALSE)
}


# The derived time series. Built by as.data.frame.Rceattle() rather than by a
# second reading of the same arrays, so this table and the figures cannot
# disagree about a confidence interval.
.rce_tab_timeseries <- function(fit, model, quantities, ci_level) {
  df <- as.data.frame(fit, which = quantities, ci_level = ci_level)
  if (!nrow(df)) return(NULL)
  endyr <- fit$data_list$endyr

  # Both depletions divide by an unfished reference, and which one depends on
  # the HCR -- see .rce_refpoint_availability(). Where that reference is the
  # underived MSSB0 / MSB0 placeholder the series is a ratio to a made-up
  # denominator, so it is blanked rather than shipped as a depletion.
  avail <- .rce_refpoint_availability(fit)
  sp <- fit$data_list$spnames %||% as.character(seq_along(avail$depletion))
  undefined <- sp[!avail$depletion]
  if (length(undefined)) {
    hit <- df$quantity %in% c("ssb_depletion", "biomass_depletion") &
      df$species %in% undefined
    df$value[hit] <- NA_real_
    df$se[hit] <- NA_real_
    df$lwr[hit] <- NA_real_
    df$upr[hit] <- NA_real_
  }

  data.frame(
    model             = model,
    species           = df$species,
    year              = df$year,
    era               = .rce_era(df$year, endyr),
    age               = df$age,
    sex               = df$sex,
    quantity          = df$quantity,
    estimate          = df$value,
    uncertainty       = df$se,
    uncertainty_label = "sd",
    lwr               = df$lwr,
    upr               = df$upr,
    stringsAsFactors  = FALSE
  )
}


# Which reference points a fit actually defines, one logical per species.
# Where it defines none the model still leaves a number behind, so see
# ?report_tables ("Reference points a fit does not define") for what each
# condition means and why it matters.
.rce_refpoint_availability <- function(fit) {
  d <- fit$data_list
  nspp <- d$nspp
  ss <- identical(as.integer(d$msmMode %||% 0L), 0L)

  # build_hcr_map() holds the gating, so it is called rather than re-read here;
  # the fit's own `map` cannot stand in, because build_map() sets both to NA in
  # the hindcast map whatever the harvest control rule. An unresolvable fit
  # counts as NOT estimated: a missing reference point is recoverable, a
  # fabricated one is not.
  f_est <- tryCatch({
    # build_hcr_map() matches HCR by name, so an integer code would match no
    # rule and blank an F40% that was estimated.
    d$HCR <- unname(.rce_alias(d$HCR, hcr_map))
    hm <- suppressMessages(build_hcr_map(d, fit$map))
    list(target = !is.na(hm$mapList$log_Ftarget),
         limit  = !is.na(hm$mapList$log_Flimit))
  }, error = function(e) list(target = rep(FALSE, nspp), limit = rep(FALSE, nspp)))

  # Single-species derives SB0/B0 itself, so the placeholder never applies.
  sb0_ok <- if (ss) rep(TRUE, nspp) else {
    dv <- d$MSSB0_derived
    if (is.null(dv)) rep(FALSE, nspp) else as.logical(dv)
  }

  # Section 12.1 of the template divides the depletions by DynamicSB0 under a
  # dynamic rule, and by biomass in the LAST projection year -- the equilibrated
  # unfished reference -- under no fishing in multispecies mode. Only the
  # remaining case reads the SB0 array, so only it inherits the placeholder.
  dyn <- isTRUE(as.logical(d$DynamicHCR %||% FALSE))
  hcr_none <- identical(unname(tryCatch(.rce_alias_show(d$HCR, hcr_map),
                                        error = function(e) NA_character_)),
                        "NoFishing")
  depl_ok <- if (dyn || (hcr_none && !ss)) rep(TRUE, nspp) else sb0_ok

  list(
    ftarget   = rep_len(f_est$target, nspp),
    flimit    = rep_len(f_est$limit, nspp),
    sb0       = rep_len(sb0_ok, nspp),
    depletion = rep_len(depl_ok, nspp),
    per_rec   = rep(ss, nspp)
  )
}


# The executive-summary reference points, per species. Biomass proxies are built
# from the model's own Ptarget / Plimit rather than hardcoded 0.40 / 0.35, so a
# model configured to another fraction reports the fraction it actually used.
#
# A quantity the fit does not define is returned as NA with a `basis` saying
# why, rather than as the number the array happens to hold.
.rce_tab_refpoints <- function(fit, model) {
  q <- fit$quantities
  d <- fit$data_list
  nspp <- d$nspp
  if (is.null(nspp) || !nspp) return(NULL)
  sp <- d$spnames %||% as.character(seq_len(nspp))

  avail <- .rce_refpoint_availability(fit)
  hcr_txt <- unname(tryCatch(.rce_alias_show(d$HCR, hcr_map),
                             error = function(e) NA_character_))

  # Reference points are taken at the terminal hindcast year, which is the year
  # the SPR schedules themselves are built from.
  yrs <- d$styr:d$projyr
  term <- match(d$endyr, yrs)

  col <- function(x, i) if (is.null(x)) NA_real_ else {
    if (is.matrix(x)) as.numeric(x[i, term]) else as.numeric(x[i])
  }

  rows <- lapply(seq_len(nspp), function(i) {
    sb0  <- col(q$SB0, i)
    ssb  <- col(q$ssb, i)
    ptar <- if (is.null(d$Ptarget)) NA_real_ else as.numeric(d$Ptarget)[i]
    plim <- if (is.null(d$Plimit))  NA_real_ else as.numeric(d$Plimit)[i]

    vals <- c(
      Ftarget            = col(q$Ftarget, i),
      Flimit             = col(q$Flimit, i),
      SPRtarget          = col(q$SPRtarget, i),
      SPRlimit           = col(q$SPRlimit, i),
      SPR0               = col(q$SPR0, i),
      B0                 = col(q$B0, i),
      SB0                = sb0,
      SBF                = col(q$SBF, i),
      B_target           = ptar * sb0,
      B_limit            = plim * sb0,
      terminal_ssb       = ssb,
      terminal_depletion = col(q$ssb_depletion, i),
      terminal_F         = col(q$F_spp, i)
    )

    # Blank anything this fit does not define, and say why. The basis text is
    # what stops a reader treating an absent reference point as a missing
    # number they could look up elsewhere.
    basis <- rep("estimated", length(vals))
    blank <- function(which_q, why) {
      hit <- names(vals) %in% which_q
      vals[hit] <<- NA_real_
      basis[hit] <<- why
    }
    if (!avail$per_rec[i]) {
      blank(c("SPR0", "SPRtarget", "SPRlimit"),
            "not defined under multispecies (M carries predation)")
    }
    if (!avail$ftarget[i]) {
      blank(c("Ftarget", "SPRtarget", "SBF"),
            paste0("not estimated under HCR = ", hcr_txt))
    }
    if (!avail$flimit[i]) {
      blank(c("Flimit", "SPRlimit"),
            paste0("not estimated under HCR = ", hcr_txt))
    }
    if (!avail$sb0[i]) {
      blank(c("SB0", "B0", "SBF", "B_target", "B_limit"),
            "unfished reference not derived (MSSB0 placeholder)")
    }
    if (!avail$depletion[i]) {
      blank("terminal_depletion",
            "unfished reference not derived (MSSB0 placeholder)")
    }
    if (isTRUE((d$estDynamics[i] %||% 0) > 0)) {
      blank(c("Ftarget", "Flimit", "SPRtarget", "SPRlimit", "SBF"),
            "input numbers-at-age: no harvest control rule, projected at F = 0")
    }

    # A quantity this fit does define can still come back non-finite from the
    # model. Saying "estimated" beside an empty cell would send the reader
    # looking for a number that was never produced.
    gone <- !is.finite(vals) & basis == "estimated"
    basis[gone] <- "not available (model returned NA)"

    data.frame(
      model    = model,
      species  = sp[i],
      year     = d$endyr,
      quantity = names(vals),
      estimate = as.numeric(vals),
      basis    = basis,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, c(rows, list(make.row.names = FALSE)))
}


# Observed against predicted for the index and catch series, with the standard
# deviation of normalized residuals per fleet -- the statistic the assessment
# guidelines ask for on fits to indices. SDNR is computed from the Pearson
# residual, which standardizes on the scale each fleet's own likelihood uses.
.rce_tab_fits <- function(fit, model) {
  res <- tryCatch(
    stats::residuals(fit, type = "pearson", source = c("index", "catch")),
    error = function(e) NULL)
  if (is.null(res) || !nrow(res)) return(NULL)

  key <- paste(res$Source, res$Fleet_code)
  sdnr <- stats::ave(res$Residual, key,
                     FUN = function(z) stats::sd(z[is.finite(z)]))

  # residuals() keys species by its integer code; every other section keys it by
  # name. Resolve to the name here so the sections join, and so a species filter
  # written against one table does not silently empty the other.
  sp <- fit$data_list$spnames
  species <- if (is.null(sp)) as.character(res$Species) else
    sp[pmin(pmax(as.integer(res$Species), 1L), length(sp))]

  data.frame(
    model      = model,
    source     = res$Source,
    fleet      = res$Fleet_code,
    fleet_name = res$Fleet_name,
    species    = species,
    year       = res$Year,
    era        = .rce_era(res$Year, fit$data_list$endyr),
    observed   = res$Observed,
    predicted  = res$Fitted,
    residual   = res$Residual,
    sdnr       = as.numeric(sdnr),
    stringsAsFactors = FALSE
  )
}


# Mohn's rho, reshaped from the wide species-per-column frame retrospective()
# returns into one row per species and forecast year.
.rce_tab_retro <- function(retro, model) {
  if (is.null(retro)) return(NULL)
  m <- retro$mohns
  if (is.null(m) || !nrow(m)) return(NULL)
  keys <- c("Object", "Forecast year", "N")
  sp_cols <- setdiff(names(m), keys)
  if (!length(sp_cols)) return(NULL)

  do.call(rbind, c(lapply(sp_cols, function(s) data.frame(
    model         = model,
    species       = s,
    quantity      = m[["Object"]],
    forecast_year = m[["Forecast year"]],
    n_peels       = m[["N"]],
    mohns_rho     = as.numeric(m[[s]]),
    stringsAsFactors = FALSE)), list(make.row.names = FALSE)))
}


# What the jitter run is for: the fraction of random starts that reached the
# best optimum found. The returned fits alone cannot say that -- non-converged
# starts are dropped before jitter() returns -- so njitter carries the
# denominator.
.rce_tab_jitter <- function(jit, model, tol = 0.01) {
  if (is.null(jit)) return(NULL)
  nll <- jit$nll[is.finite(jit$nll)]
  tried <- jit$njitter %||% NA_integer_
  best <- if (length(nll)) min(nll) else NA_real_
  data.frame(
    model          = model,
    n_started      = as.integer(tried),
    n_converged    = length(nll),
    n_at_best      = if (length(nll)) sum(nll - best <= tol) else 0L,
    best_nll       = best,
    worst_nll      = if (length(nll)) max(nll) else NA_real_,
    nll_range      = if (length(nll)) max(nll) - best else NA_real_,
    tolerance      = tol,
    stringsAsFactors = FALSE
  )
}


# One-step-ahead residual diagnostics, summarized per source and fleet.
.rce_tab_osa <- function(osa, model) {
  if (is.null(osa)) return(NULL)
  d <- tryCatch(osa_diagnostics(osa), error = function(e) NULL)
  if (is.null(d) || !nrow(d)) return(NULL)
  cbind(model = model, as.data.frame(d), stringsAsFactors = FALSE)
}


#' Print method for a set of assessment reporting tables
#'
#' @param x An `"rceattle_report"` from [report_tables()].
#' @param ... Currently unused.
#' @return `x`, invisibly.
#' @export
print.rceattle_report <- function(x, ...) {
  mods <- unique(x$model$model)
  cat("<Rceattle report tables>\n")
  cat("  models : ", paste(mods, collapse = ", "), "\n", sep = "")
  cat("  species: ", paste(unique(x$model$species), collapse = "; "), "\n",
      sep = "")
  cat("  sections\n")
  nms <- names(x)
  kw <- max(nchar(nms))
  for (i in seq_along(nms)) {
    # Box-drawing glyphs as \u escapes: a portable package keeps its R code
    # ASCII, and the printed characters are identical either way.
    glyph <- if (i == length(nms)) "\u2514\u2500" else "\u251c\u2500"
    cat(sprintf("  %s %s : %d rows\n", glyph,
                formatC(nms[i], width = kw, flag = "-"), nrow(x[[nms[i]]])))
  }
  missing <- setdiff(c("retrospective", "jitter", "osa"), nms)
  if (length(missing)) {
    cat("  not supplied: ", paste(missing, collapse = ", "),
        " (pass retro = / jitter = / osa =)\n", sep = "")
  }
  invisible(x)
}


# The NOAA standardized assessment output, in column order. `stockplotr` and
# `asar` read this schema; every column is present even when empty, because
# their filters index by name and a missing column is an error there rather than
# an absence. These are the 34 names `stockplotr::example_data` carries;
# `test-report-tables.R` checks the emitted names against that package.
.RCE_STANDARD_COLS <- c(
  "module_name", "label", "time", "era", "year", "month", "season",
  "subseason", "birthseas", "initial", "estimate", "uncertainty",
  "uncertainty_label", "likelihood", "fleet", "platoon", "area", "age", "sex",
  "growth_pattern", "bio_pattern", "settlement", "morph", "beg_mid", "type",
  "factor", "part", "kind", "nsim", "bin", "age_a", "length_bins", "count",
  "block"
)


# Pad a frame out to the standard schema and put the columns in its order. Extra
# columns (`species`, `model`) are kept on the end rather than dropped -- see
# standard_output() on why species cannot be folded into the standard itself.
.rce_standard_frame <- function(df) {
  extra <- setdiff(names(df), .RCE_STANDARD_COLS)
  # rep() rather than a bare NA, so a zero-row frame pads to an empty table of
  # the right shape instead of erroring on the recycling.
  for (nm in setdiff(.RCE_STANDARD_COLS, names(df)))
    df[[nm]] <- rep(NA, nrow(df))
  df[, c(.RCE_STANDARD_COLS, extra), drop = FALSE]
}


#' Rceattle output in the NOAA standardized assessment format
#'
#' @description Relabels [report_tables()] output into the standardized
#' assessment output that the `stockplotr` and `asar` packages consume, so
#' Rceattle results can be plotted and written into a report by the same tools
#' used for SS3, BAM, WHAM and FIMS.
#'
#' @details
#' Quantity names are translated through the `standard_label` column of
#' [quantity_dictionary()], so `ssb` becomes `spawning_biomass`, `R` becomes
#' `recruitment`, `F_spp` becomes `fishing_mortality`, and so on. A quantity the
#' standard has no name for keeps its Rceattle name, so nothing is silently
#' dropped.
#'
#' **The standard has no species dimension.** It describes one stock, so a
#' multispecies CEATTLE fit cannot be represented in it as a whole. A `species`
#' column is held alongside the standard columns and `species` selects one
#' stock; with several species in the fit and no selection, this errors rather
#' than returning a frame in which two stocks' biomass share a year.
#'
#' `era` follows the standard's own vocabulary: `"time"` for the hindcast the
#' model was fit to, `"fore"` for the projection.
#'
#' @param x An `"rceattle_report"` from [report_tables()], or an Rceattle fit,
#'   which is passed through [report_tables()] first.
#' @param species Which stock to emit: a species name, or an index into the
#'   species this report holds, required when it holds more than one.
#' @param model Which model to emit when `x` holds several, defaulting to the
#'   first.
#'
#' @return A data frame with the standardized columns, in the standard's own
#'   order, plus `species` and `model`.
#'
#' @seealso [report_tables()] for the native tables and
#'   [quantity_dictionary()] for the name crosswalk.
#'
#' @examples
#' \dontrun{
#' fit <- fit_mod(data_list = BS2017SS, estimateMode = 1)
#' std <- standard_output(fit, species = 1)
#' # ...then, with the stockplotr package installed:
#' # stockplotr::plot_spawning_biomass(std)
#' }
#' @export
standard_output <- function(x, species = NULL, model = NULL) {
  if (inherits(x, "Rceattle")) x <- report_tables(x)
  if (!inherits(x, "rceattle_report")) {
    stop("`x` must be an Rceattle fit or an rceattle_report from report_tables().",
         call. = FALSE)
  }

  mods <- unique(x$model$model)
  if (is.null(model)) {
    model <- mods[1]
    if (length(mods) > 1) {
      message("Emitting model '", model, "'; pass `model` to choose another.")
    }
  }
  if (!model %in% mods) {
    stop("`model` must be one of: ", paste(mods, collapse = ", "), call. = FALSE)
  }

  # Species are read back off the tables, not by splitting the `model` row's
  # comma-joined summary string: a stock whose own name contains ", " ("Pollock,
  # GOA") would otherwise split into pieces and become unselectable.
  all_species <- unique(unlist(lapply(
    x[c("timeseries", "reference_points", "fits")], function(tb) {
      if (is.null(tb) || !"species" %in% names(tb)) return(NULL)
      tb$species[tb$model == model]
    }), use.names = FALSE))
  all_species <- all_species[!is.na(all_species)]
  if (!length(all_species)) {
    stop("No species found for model '", model, "' in this report.", call. = FALSE)
  }
  if (is.null(species)) {
    if (length(all_species) > 1) {
      stop("The standardized output describes ONE stock and has no species ",
           "dimension, but this fit has ", length(all_species), " (",
           paste(all_species, collapse = ", "), "). Pass `species` to choose one.",
           call. = FALSE)
    }
    species <- all_species[1]
  }
  if (is.numeric(species)) species <- all_species[species]
  if (!species %in% all_species) {
    stop("`species` must be one of: ", paste(all_species, collapse = ", "),
         call. = FALSE)
  }

  # Rceattle name -> standard label, falling back to the Rceattle name so a
  # quantity the standard has no word for is still emitted under its own.
  dict <- quantity_dictionary()
  relabel <- function(q) {
    i <- match(q, dict$quantity)
    lab <- dict$standard_label[i]
    ifelse(is.na(lab), q, lab)
  }

  keep <- function(df, sp_col = "species") {
    if (is.null(df) || !nrow(df)) return(NULL)
    df <- df[df$model == model, , drop = FALSE]
    if (sp_col %in% names(df)) {
      df <- df[is.na(df[[sp_col]]) | df[[sp_col]] == species, , drop = FALSE]
    }
    if (!nrow(df)) NULL else df
  }

  parts <- list()

  ts <- keep(x$timeseries)
  if (!is.null(ts)) {
    parts$timeseries <- .rce_standard_frame(data.frame(
      module_name = "timeseries", label = relabel(ts$quantity),
      year = ts$year, era = ts$era, estimate = ts$estimate,
      uncertainty = ts$uncertainty, uncertainty_label = ts$uncertainty_label,
      age = ts$age, sex = ts$sex, species = ts$species, model = model,
      stringsAsFactors = FALSE))
  }

  rp <- keep(x$reference_points)
  if (!is.null(rp)) {
    # `basis` is not part of the standard, so it rides along as an extra column:
    # without it an NA estimate here looks like a number that went missing
    # rather than one this fit never defined.
    parts$reference_points <- .rce_standard_frame(data.frame(
      module_name = "reference_points", label = relabel(rp$quantity),
      year = rp$year, estimate = rp$estimate,
      species = rp$species, model = model, basis = rp$basis,
      stringsAsFactors = FALSE))
  }

  # `unit` names a fleet on rows 1-8 and a species on rows 9-20, so only species
  # rows are filtered here. Fleet rows are kept whatever the species and carried
  # in the standard's `fleet` column -- they are most of the decomposition a
  # model comparison is built from (31 of 38 rows on a 16-fleet assessment).
  lk <- x$likelihood
  if (!is.null(lk) && nrow(lk)) {
    lk <- lk[lk$model == model, , drop = FALSE]
    is_sp <- !is.na(lk$axis) & lk$axis == "species"
    lk <- lk[!is_sp | lk$unit == species, , drop = FALSE]
  }
  if (!is.null(lk) && nrow(lk)) {
    is_flt <- !is.na(lk$axis) & lk$axis == "fleet"
    parts$likelihood <- .rce_standard_frame(data.frame(
      module_name = "likelihood", label = lk$component,
      likelihood = lk$weighted, estimate = lk$weighted,
      fleet = ifelse(is_flt, lk$unit, NA_character_),
      species = species, model = model,
      stringsAsFactors = FALSE))
  }

  ft <- keep(x$fits)
  if (!is.null(ft)) {
    obs <- .rce_standard_frame(data.frame(
      module_name = "fits", label = paste0(ft$source, "_observed"),
      year = ft$year, era = ft$era, estimate = ft$observed,
      fleet = ft$fleet_name, species = ft$species, model = model,
      stringsAsFactors = FALSE))
    pred <- .rce_standard_frame(data.frame(
      module_name = "fits", label = paste0(ft$source, "_predicted"),
      year = ft$year, era = ft$era, estimate = ft$predicted,
      fleet = ft$fleet_name, species = ft$species, model = model,
      stringsAsFactors = FALSE))
    parts$fits <- rbind(obs, pred, make.row.names = FALSE)
  }

  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (!length(parts)) {
    return(.rce_standard_frame(data.frame(species = character(),
                                          model = character(),
                                          stringsAsFactors = FALSE)))
  }
  # Sections carry different extra columns beyond the standard (only the
  # reference points have a `basis`), so pad every part to their union before
  # binding rather than requiring them to agree.
  all_cols <- unique(unlist(lapply(parts, names), use.names = FALSE))
  parts <- lapply(parts, function(p) {
    for (nm in setdiff(all_cols, names(p))) p[[nm]] <- NA
    p[, all_cols, drop = FALSE]
  })
  do.call(rbind, c(unname(parts), list(make.row.names = FALSE)))
}
