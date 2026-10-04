#' Fit the CEATTLE assessment model
#' @description Estimate CEATTLE population parameters by maximum likelihood, and
#'   optionally project the stock and apply a harvest control rule.
#'
#' @param data_list A data list read in via \code{\link{read_data}} or built
#'   directly in R; see \code{vignette("data-without-excel", package = "Rceattle")}.
#' @param inits (Optional) A named list of initial parameter values, as returned by
#'   \code{\link{build_params}} or extracted from a previous fit
#'   (\code{model$estimated_params}). If \code{NULL}, parameters are initialized
#'   from scratch via \code{\link{build_params}}.
#' @param map (Optional) A map object from \code{\link{build_map}}.
#' @param bounds (Optional) A bounds object from \code{\link{build_bounds}}.
#' @param file (Optional) Filename where files will be saved. If NULL, no file is saved.
#' @param estimateMode What to fit, given as a string alias or the integer code:
#'   \code{"Estimate"} (0) = fit the hindcast model and the HCR projection
#'   (\code{HCR}); \code{"Hindcast"} (1) = fit the hindcast only (no fitting BRPs/HCR/projection);
#'   \code{"Projection"} (2) = fit the BRPs/HCR/projection only, from the initial
#'   parameters in \code{inits}; \code{"DebugBuild"} (3) = build through
#'   \code{MakeADFun} but not \code{nlminb} -- the returned \code{obj} carries the
#'   real objective and gradient, so \code{obj$fn()} / \code{obj$gr()} are usable
#'   for diagnosing a model before committing to a fit; \code{"DebugOptimize"}
#'   (4) = optimize with all parameters mapped out, so the objective is a
#'   placeholder (\code{dummy^2}), not a likelihood. Defaults to \code{"Estimate"}.
#' @param random_rec logical. If TRUE, treats recruitment deviations as random effects using the Laplace approximation. The default is FALSE.
#' @param random_q logical. If TRUE, treats annual catchability deviations as random effects using the Laplace approximation, and estimates their standard deviation rather than fixing it at `Time_varying_q_sd`. The default is FALSE.
#' @param random_sel logical. If TRUE, treats annual selectivity deviations as random effects using the Laplace approximation, and estimates their standard deviation rather than fixing it at `Time_varying_sel_sd`. The default is FALSE.
#' @param HCR HCR list object from \code{\link{build_hcr}}
#' @param niter Number of iterations for multispecies model
#' @param recFun The stock recruit-relationship parameterization from \code{\link{build_srr}}.
#' @param M1Fun M1 parameterizations and priors. Use \code{build_M1}.
#' @param growthFun The weight-at-age parameterization from \code{\link{build_growth}}.
#' @param qFun Catchability specification from \code{\link{build_catchability}},
#'   carrying any environmental linkages on q.
#' @param selFun Selectivity specification from \code{\link{build_selectivity}}, carrying any environmental linkages on selectivity parameters.
#' @param compFun Composition-weighting specification from \code{\link{build_composition}}, carrying any priors on the Dirichlet-multinomial weights.
#' @param dsem Optional dynamic structural equation model on the recruitment
#'   deviations, from \code{\link{build_DSEM}}. The recruitment deviations
#'   become latent states of a GMRF, so environmental covariates and
#'   lagged/simultaneous paths drive recruitment. \code{NULL} (the default)
#'   fits the standard recruitment-deviation parameterization and is a complete
#'   no-op. A DSEM and a recruitment linkage (\code{recFun = build_srr(linkages
#'   = ...)}) both structure the recruitment deviations and cannot be combined.
#' @param msmMode The predation-mortality mode, as a string alias or integer code:
#'   \code{"SingleSpecies"} (0, the default, no predation), \code{"MSVPA"}
#'   (1, the Type-II MSVPA predation of Holsman et al. 2015) or
#'   \code{"TypeIIIMSVPA"} (2). Higher integer
#'   codes (Kinzey-Punt, Holling forms) are declared but not implemented, and are
#'   rejected by the data check.
#' @param avgnMode the average abundance-at-age approximation used in the predation-mortality equations. Only mode 0, \eqn{N/Z(1 - exp(-Z))} (the MSVPA form), is implemented; the alternatives are declared but have no effect.
#' @param initMode how the population is initialized, as a string alias or integer
#'   code: \code{"FreeParams"} (0), \code{"Equilibrium"} (1),
#'   \code{"NonEquilibrium"} (2, the default), \code{"FishedNonEquilibrium"} (3),
#'   \code{"FishedNonEquilibriumScaled"} (4), \code{"OffsetEquilibrium"} (5). See
#'   the \strong{Initial age structure} section below for what each one estimates.
#' @param suitMode how predator-prey suitability is derived, per predator (a single value or a vector of length \code{nspp}): 0 = empirical from diet data (Holsman et al. 2015), 2 = weight-based gamma, 4 = weight-based lognormal, 6 = weight-based normal. The length-based forms (1, 3, 5) are declared but not implemented and are rejected by the data check.
#' @param suit_styr The first year used to calculate mean suitability. A single integer is applied to every predator, or a vector of length `nspp` sets a distinct start year per predator. Defaults to `styr` in `data_list`. Used when diet data were sampled from a subset of years.
#' @param suit_endyr The last year used to calculate mean suitability. A single integer is applied to every predator, or a vector of length `nspp` sets a distinct end year per predator. Defaults to `endyr` in `data_list`. Used when diet data were sampled from a subset of years.
#' @param fit_control A list returned by [fit_control()] that bundles the
#'   optimizer / sdreport / phasing knobs (`phase`, `bias.correct`, `getsd`,
#'   `getJointPrecision`, `getReportCovariance`, `use_gradient`, `rel_tol`,
#'   `loopnum`, `newtonsteps`, `TMBfilename`, `verbose`, `nlminb_control`).
#'   Defaults to `fit_control()`. See [fit_control()] for the meaning and
#'   defaults of each field.
#' @param config (Optional) An `Rceattle_run_config` from [load_config()] (or
#'   [run_config()]). Its stored `model_config` structure and estimation controls
#'   (`estimateMode`, `random_rec`/`random_q`/`random_sel`, `suit_styr`/
#'   `suit_endyr`, `fit_control`) overlay only the arguments the caller did *not*
#'   pass -- an explicit argument always wins. `NULL` (default) applies no
#'   configuration. Example: `fit_mod(data_list, config = load_config("run.yaml"))`.
#' @param quiet_data_check Drop the warnings the fit-time validation raises (errors still
#'   stop the fit). `FALSE` (default) for an ordinary fit. The diagnostic refits
#'   -- [retrospective()], [jitter()], [self_test()], [profile()], [run_mse()],
#'   [remove_F()], [sample_rec()], [reweight_comps()] -- set it, since they
#'   re-validate a `data_list` the caller has already fitted once and would
#'   otherwise repeat the same warnings per peel, jitter, or MSE iteration.
#'   Convergence and TMB warnings are unaffected.
#' @param ... Deprecated optimizer / sdreport / phasing arguments
#'   (e.g. `phase`, `getsd`, `bias.correct`, `use_gradient`, `rel_tol`,
#'   `control`, `getJointPrecision`, `getReportCovariance`, `loopnum`,
#'   `newtonsteps`, `verbose`, `TMBfilename`). These are forwarded into
#'   `fit_control` with a deprecation warning; pass them via
#'   [fit_control()] instead.
#'
#' @details
#' CEATTLE is an age-structured population dynamics model that can be fit with or without predation mortality. The default is to exclude predation mortality by setting \code{msmMode} to 0. Predation mortality can be included by setting \code{msmMode} with the following options:
#' \itemize{
#' \item{0. Single species mode}
#' \item{1. Holsman et al. 2015 predation based on multi-species virtual population analysis (MSVPA) based predation formation.}
#' \item{2. MSVPA Holling Type III}
#' }
#'
#' Values 3 through 9 (Kinzey & Punt 2009 functional responses --
#' Holling Type I/II/III, predator interference, predator preemption,
#' Hassell-Varley, Ecosim) are blocked at runtime by \code{data_check()}
#' because the implementations have not been validated against the
#' current parameter set. See \code{src/TMB/predation.hpp}.
#'
#'
#' @section Initial age structure:
#' What \code{initMode} estimates, and from what:
#' \describe{
#'   \item{\code{"FreeParams"} (0)}{The initial age-structure is estimated directly, one free
#'     parameter per age.}
#'   \item{\code{"Equilibrium"} (1)}{Unfished (\eqn{F_{init} = 0}) equilibrium age-structure,
#'     carried out from \eqn{R_0} and residual natural mortality \eqn{M1}.}
#'   \item{\code{"NonEquilibrium"} (2)}{As (1), plus estimated initial population deviates, so
#'     the first year need not sit at equilibrium. The default.}
#'   \item{\code{"FishedNonEquilibrium"} (3)}{As (2), with an estimated initial fishing
#'     mortality \eqn{F_{init}} added to the mortality that shapes the first year.}
#'   \item{\code{"FishedNonEquilibriumScaled"} (4)}{As (3), but \eqn{F_{init}} scales
#'     \eqn{R_0} rather than entering the mortality directly.}
#'   \item{\code{"OffsetEquilibrium"} (5)}{Unfished equilibrium seeded by the first year's
#'     recruitment, \code{R_init * exp(rec_dev[1])}, decayed by \eqn{M1} and closed with the
#'     usual geometric plus group. Initial deviates are turned off and no init-dev penalty is
#'     applied. This is the Cole Monnahan / AFSC GOA pollock convention.}
#' }
#'
#' Modes 1 and 5 differ by exactly one term: both start from the initial equilibrium
#' recruitment \eqn{R_{init}}, but (1) carries it forward unchanged while (5) seeds the first
#' year with the realized recruitment \code{R_init * exp(rec_dev[1])}. On a stock whose first
#' year was not average, that is not a small difference.
#'
#' \eqn{R_{init}} is not always \eqn{R_0}: under a random-about-mean stock-recruit relationship
#' the two are equal, but under Beverton-Holt or Ricker \eqn{R_{init}} is the equilibrium
#' recruitment implied by the curve.
#'
#' **In multispecies mode the decay uses \eqn{M1} only**, so predation mortality (\eqn{M2}) does
#' not enter the initial age structure under any of these modes.
#'
#' @return A list of class "Rceattle" including:
#'
#' \itemize{
#'  \item{data_list: List of data inputs}
#'  \item{initial_params: List of starting parameters}
#'  \item{bounds: Parameter bounds used for estimation}
#'  \item{map: List of map used in TMB}
#'  \item{obj: TMB model object}
#'  \item{opt: Optimized model object from `nlminb`}
#'  \item{sdrep: Object of class `sdreport` exported by TMB including the standard errors of estimated parameters}
#'  \item{estimated_params: List of estimated parameters}
#'  \item{quantities: Derived quantities from CEATTLE}
#'  \item{run_time: Model run time}
#'  }
#'
#'
#'
#'
#' @examples
#' \donttest{
#' data(BS2017SS)
#' ss_run <- fit_mod(
#'   data_list    = BS2017SS,
#'   estimateMode = 0,
#'   msmMode      = 0,
#'   fit_control  = fit_control(phase = FALSE, verbose = 0)
#' )
#' }
#'
#' @export
fit_mod <-
  function(
    data_list = NULL,
    inits = NULL,
    map = NULL,
    bounds = NULL,
    file = NULL,
    estimateMode = 0,
    random_rec = FALSE,
    random_q = FALSE,
    random_sel = FALSE,
    HCR = build_hcr(),
    niter = 3,
    recFun = build_srr(),
    M1Fun = build_M1(),
    growthFun = build_growth(),
    qFun = build_catchability(),
    selFun = build_selectivity(),
    compFun = build_composition(),
    dsem = NULL,
    msmMode = 0,
    avgnMode = 0,
    initMode = "NonEquilibrium",
    suitMode = 0,
    suit_styr = NULL,
    suit_endyr = NULL,
    fit_control = NULL,
    config = NULL,
    quiet_data_check = FALSE,
    ...){

    # Whether the caller supplied fit_control -- captured before the default is
    # filled below, since assigning to the argument makes missing() unreliable.
    # (The config= overlay needs this to know if it may fill fit_control.)
    fc_supplied <- !missing(fit_control)

    # Default control bundle. Built inside the body to avoid the
    # "recursive default argument reference" error that occurs when a
    # parameter and the function used in its default share a name.
    if (is.null(fit_control)) fit_control <- Rceattle::fit_control()

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # Deprecation catch ----
    # Older versions of fit_mod() accepted the optimizer / sdreport /
    # phasing knobs as individual arguments. They now live exclusively
    # on `fit_control`. Catch the old names from `...`, warn, and
    # forward into the supplied fit_control so existing scripts keep
    # working.
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    .deprecated_ctl_args <- c(
      "phase", "getsd", "bias.correct", "use_gradient", "rel_tol",
      "control", "getJointPrecision", "getReportCovariance",
      "loopnum", "newtonsteps", "verbose", "TMBfilename",
      "projection_uncertainty"
    )
    .extra  <- list(...)
    .legacy <- intersect(names(.extra), .deprecated_ctl_args)
    if (length(.legacy)) {
      warning(
        sprintf(
          paste0(
            "Passing %s directly to fit_mod() is deprecated and will be ",
            "removed in a future release. Bundle these into fit_control() ",
            "instead, e.g. fit_control(%s). Forwarding for now."
          ),
          paste(sQuote(.legacy), collapse = ", "),
          paste(.legacy, "= ...", collapse = ", ")
        ),
        call. = FALSE
      )
      for (.nm in .legacy) {
        .field <- if (.nm == "control") "nlminb_control" else .nm
        fit_control[[.field]] <- .extra[[.nm]]
      }
      .extra <- .extra[setdiff(names(.extra), .deprecated_ctl_args)]
    }
    if (length(.extra)) {
      stop(
        "Unused arguments to fit_mod(): ",
        paste(sQuote(names(.extra)), collapse = ", "),
        call. = FALSE
      )
    }

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # Unpack fit_control ----
    # The optimizer / sdreport / phasing knobs all live in fit_control;
    # expand them to local variables so the rest of the function body
    # can refer to them by name without rewriting every reference.
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    if (!inherits(fit_control, "Rceattle_fit_control")) {
      stop("`fit_control` must be created with fit_control().")
    }
    bias.correct        <- fit_control$bias.correct
    getsd               <- fit_control$getsd
    getJointPrecision   <- fit_control$getJointPrecision
    getReportCovariance <- fit_control$getReportCovariance
    projection_uncertainty <- fit_control$projection_uncertainty
    use_gradient        <- fit_control$use_gradient
    rel_tol             <- fit_control$rel_tol
    loopnum             <- fit_control$loopnum
    newtonsteps         <- fit_control$newtonsteps
    phase               <- fit_control$phase
    TMBfilename         <- fit_control$TMBfilename
    verbose             <- fit_control$verbose
    control             <- fit_control$nlminb_control

    # ---------------------------------------------------------------------
    # Pipeline overview (file prefixes match this execution order):
    #   0-clean_data.R            clean_data() / 0-switches.R switch_check()
    #   1-data_check.R            data_check()      validate inputs
    #   2-build_params.R          build_params()    starting parameter list
    #   3-build_map.R             build_map()       TMB map (fixed vs estimated)
    #   4-build_parameter_bounds.R build_bounds()   lower/upper bounds
    #   5-rearrange_data.R        rearrange_data()  reshape data for TMB
    #   6-fit_mod.R               (this file)       MakeADFun + nlminb + sdreport
    #   6-rename_output.R         rename_output()   label derived quantities
    # HCR map (0-build_hcr.R build_hcr_map()) is applied during projection below.
    # ---------------------------------------------------------------------

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # 0 - Start ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    mod_objects <- list()
    start_time  <- Sys.time()

    # Not every estimateMode / HCR combination runs an optimization, and the
    # reporting block below reads `opt$SD`.
    opt <- NULL

    extend_length <- function(x) {
      # A scalar recycles to every species; a length-nspp vector passes through.
      # A multi-element vector of any other length is a misconfiguration (it
      # would otherwise be silently re-recycled, e.g. length 3 with nspp 2 -> 6).
      if (length(x) > 1L && length(x) != data_list$nspp) {
        stop("Expected a scalar or a length-", data_list$nspp,
             " (nspp) vector; got length ", length(x), call. = FALSE)
      }
      if (length(x) == data_list$nspp) x else rep(x, data_list$nspp)
    }

    # Prefer function arg > existing data_list field > fallback year
    resolve_yr <- function(arg, data_val, fallback) {
      if (!is.null(arg)) arg else if (!is.null(data_val)) data_val else fallback
    }


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # 1 - Load data and switches ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    if (is.null(data_list)) {
      stop("Missing data_list object")
    }

    data_list <- Rceattle::clean_data(data_list)

    # Overlay a run configuration (from load_config()) onto arguments the caller
    # did not supply. Its model_config structure is attached to the data list so
    # the model_config overlay below applies it to the omitted structure args;
    # its estimation controls fill the omitted scalar args. An explicitly-passed
    # argument always wins (missing() / fc_supplied is FALSE). config = NULL
    # (default) is a complete no-op, so the fit path is unchanged.
    if (!is.null(config)) {
      if (!inherits(config, "Rceattle_run_config"))
        stop("`config` must be an Rceattle_run_config (from load_config() / run_config()).",
             call. = FALSE)
      if (!is.null(config$model_config)) data_list$model_config <- config$model_config
      if (missing(estimateMode) && !is.null(config$estimateMode)) estimateMode <- config$estimateMode
      if (missing(random_rec)   && !is.null(config$random_rec))   random_rec   <- config$random_rec
      if (missing(random_q)     && !is.null(config$random_q))     random_q     <- config$random_q
      if (missing(random_sel)   && !is.null(config$random_sel))   random_sel   <- config$random_sel
      if (missing(suit_styr)    && !is.null(config$suit_styr))    suit_styr    <- config$suit_styr
      if (missing(suit_endyr)   && !is.null(config$suit_endyr))   suit_endyr   <- config$suit_endyr
      if (!fc_supplied          && !is.null(config$fit_control))  fit_control  <- config$fit_control
    }

    # Overlay a stored model_config (if any) onto arguments the caller did not
    # supply. An explicitly-passed argument always wins (missing() is FALSE); an
    # argument the caller omitted falls back to the stored field. With no
    # model_config slot this is a no-op.
    # Detected with missing() rather than a sentinel so 0 / FALSE / "" stored
    # values are honoured. NOTE: a call that passes one of these arguments --
    # even at its default -- overrides the slot; omit it to let the config win.
    cfg <- data_list$model_config
    if (!is.null(cfg)) {
      if (missing(msmMode)   && !is.null(cfg$msmMode))   msmMode   <- cfg$msmMode
      if (missing(initMode)  && !is.null(cfg$initMode))  initMode  <- cfg$initMode
      if (missing(avgnMode)  && !is.null(cfg$avgnMode))  avgnMode  <- cfg$avgnMode
      if (missing(suitMode)  && !is.null(cfg$suitMode))  suitMode  <- cfg$suitMode
      if (missing(niter)     && !is.null(cfg$niter))     niter     <- cfg$niter
      if (missing(HCR)       && !is.null(cfg$HCR))       HCR       <- cfg$HCR
      if (missing(recFun)    && !is.null(cfg$recFun))    recFun    <- cfg$recFun
      if (missing(M1Fun)     && !is.null(cfg$M1Fun))     M1Fun     <- cfg$M1Fun
      if (missing(growthFun) && !is.null(cfg$growthFun)) growthFun <- cfg$growthFun
      if (missing(qFun)      && !is.null(cfg$qFun))      qFun      <- cfg$qFun
      if (missing(selFun)    && !is.null(cfg$selFun))    selFun    <- cfg$selFun
      if (missing(compFun)   && !is.null(cfg$compFun))   compFun   <- cfg$compFun
      if (missing(dsem)      && !is.null(cfg$dsem))      dsem      <- cfg$dsem
    }

    # Validate the mode switches HERE, before any work happens. model_config()
    # re-checks them when the run configuration is recorded, but that call runs
    # after the optimization and sdreport -- a throw there would discard a
    # converged fit over a typo in an argument the model never reads.
    .check_switch(initMode, initMode_map, "initMode")
    .check_switch(suitMode, suitMode_map, "suitMode")
    if (!is.null(avgnMode) && length(avgnMode) &&
        !all(is.na(avgnMode) | avgnMode %in% 0:2)) {
      stop("`avgnMode` must be 0, 1, or 2 (only 0 is implemented).", call. = FALSE)
    }

    # Resolve string aliases for the model-mode switches (e.g. "Hindcast" -> 1,
    # "SingleSpecies" -> 0); integers pass through unchanged and an unknown string
    # errors clearly. Done here, after the config overlay, so every downstream
    # integer test on estimateMode / msmMode still applies.
    estimateMode <- .map_switch(estimateMode, estimateMode_map, "estimateMode")
    msmMode      <- .map_switch(msmMode,      msmMode_map,      "msmMode")
    .check_switch(msmMode, msmMode_map, "msmMode")

    # Add switches from function call
    data_list$random_rec  <- as.numeric(random_rec)
    data_list$random_q    <- as.numeric(random_q)
    data_list$random_sel  <- as.numeric(random_sel)
    data_list$estimateMode <- estimateMode
    data_list$niter       <- niter
    data_list$avgnMode    <- avgnMode
    data_list$initMode    <- initMode
    data_list$loopnum     <- loopnum
    data_list$msmMode     <- msmMode
    data_list$suitMode    <- extend_length(suitMode)

    # * Suitability switches ----
    # suit_styr/suit_endyr are per-predator: a scalar is recycled to all species,
    # a length-nspp vector sets a distinct suitability-averaging window per predator.
    data_list$suit_styr  <- extend_length(resolve_yr(suit_styr,  data_list$suit_styr,  data_list$styr))
    data_list$suit_endyr <- extend_length(resolve_yr(suit_endyr, data_list$suit_endyr, data_list$endyr))


    # * Recruitment switches ----
    data_list$srr_fun       <- recFun$srr_fun
    data_list$srr_pred_fun  <- recFun$srr_pred_fun
    data_list$proj_mean_rec <- recFun$proj_mean_rec

    data_list$srr_mse_switchyr <- resolve_yr(recFun$srr_mse_switchyr, data_list$srr_mse_switchyr, data_list$endyr)
    data_list$srr_hat_styr     <- resolve_yr(recFun$srr_hat_styr,     data_list$srr_hat_styr,     data_list$styr + 1)
    data_list$srr_hat_endyr    <- resolve_yr(recFun$srr_hat_endyr,    data_list$srr_hat_endyr,    data_list$endyr)

    data_list$srr_est_mode <- recFun$srr_est_mode
    data_list$srr_prior    <- extend_length(recFun$srr_prior)
    data_list$srr_prior_sd <- extend_length(recFun$srr_prior_sd)
    data_list$srr_alpha_init <- if (is.null(recFun$srr_alpha_init)) NULL else extend_length(recFun$srr_alpha_init)
    data_list$srr_beta_init  <- if (is.null(recFun$srr_beta_init))  NULL else extend_length(recFun$srr_beta_init)
    data_list$srr_linkages <- recFun$linkages
    data_list$Bmsy_lim     <- extend_length(recFun$Bmsy_lim)

    # * M switches ----
    # Fold the deprecated `est_M1` into `M1_model` here -- BEFORE the M1Fun
    # reconciliation below -- so it gets the same treatment as `M1_model`
    # (the switch_check() alias runs later, after M1_model is already resolved
    # from M1Fun, so it cannot help the fit path).
    data_list <- .alias_est_M1(data_list)
    if (!is.null(data_list$M1_model)) {
      if (sum(data_list$M1_model != extend_length(M1Fun$M1_model))) {
        warning("M1_model in data is different than in call `fit_mod`, using switch from 'fit_mod'")
      }
    }

    # FIXME: may want to pull from data here too??
    data_list$M1_model     <- extend_length(M1Fun$M1_model)
    data_list$M1_model     <- ifelse(data_list$nsex == 1 & data_list$M1_model == 2, 1, data_list$M1_model) # sex-specific -> sex-invariant for 1-sex model
    data_list$M1_re        <- extend_length(M1Fun$M1_re)
    updateM1               <- M1Fun$updateM1
    data_list$M1_use_prior <- extend_length(M1Fun$M1_use_prior) * (data_list$M1_model > 0) # 0 when M1 is fixed
    data_list$M2_use_prior <- extend_length(M1Fun$M2_use_prior) * (msmMode > 0)            # 0 in single-species mode
    data_list$M_prior      <- extend_length(M1Fun$M_prior)
    data_list$M_prior_sd   <- extend_length(M1Fun$M_prior_sd)
    data_list$M1_indices   <- M1Fun$M1_indices
    data_list$M1_linkages  <- M1Fun$linkages


    # * Growth switches ----
    data_list$growth_fun      <- growthFun$fun
    data_list$growth_linkages <- growthFun$linkages
    data_list$q_linkages      <- qFun$linkages
    data_list$sel_linkages    <- selFun$linkages
    data_list$comp_linkages   <- compFun$linkages
    data_list$growth_model    <- extend_length(growthFun$growth_model)
    # Plus-group SD-at-age style, resolved like growth_age_L1 below:
    # build_growth() value (if given) > existing data_list$growth_sd_style
    # (so a refit that rebuilds growth via build_growth(fun=) keeps the
    # original SS3/WHAM choice) > WHAM (1) fallback.
    gsd <- extend_length(growthFun$growth_sd_style)
    if (!is.null(data_list$growth_sd_style)) {
      gsd[is.na(gsd)] <- extend_length(data_list$growth_sd_style)[is.na(gsd)]
    }
    gsd[is.na(gsd)] <- 1L   # WHAM
    data_list$growth_sd_style <- gsd
    # VB anchor age per species (= SS3 Growth_Age_for_L1). Resolution
    # order: build_growth() user arg > data_list$growth_age_L1 (e.g. from
    # ss3_to_rceattle converter) > max(0.5, minage[sp]) fallback. The
    # fallback keeps minage >= 1 models backwards-compatible and gives
    # minage = 0 models an SS3-style half-year anchor.
    gal1 <- extend_length(growthFun$growth_age_L1)
    if (!is.null(data_list$growth_age_L1)) {
      gal1[is.na(gal1)] <- extend_length(data_list$growth_age_L1)[is.na(gal1)]
    }
    gal1[is.na(gal1)] <- pmax(0.5, as.numeric(extend_length(data_list$minage)))[is.na(gal1)]
    data_list$growth_age_L1 <- gal1


    # * HCR Switches ----
    data_list$HCR        <- HCR$HCR
    data_list$DynamicHCR <- HCR$DynamicHCR
    if (!HCR$HCR %in% c(2, "ConstantF")) { # Ftarget is also used for fixed F (so may be of length nflts)
      data_list$Ftarget <- extend_length(HCR$Ftarget)
    } else {
      data_list$Ftarget <- HCR$Ftarget
    }
    data_list$Flimit   <- extend_length(HCR$Flimit)
    data_list$Ptarget  <- extend_length(HCR$Ptarget)
    data_list$Plimit   <- extend_length(HCR$Plimit)
    data_list$Alpha    <- extend_length(HCR$Alpha)
    data_list$Pstar    <- extend_length(HCR$Pstar)
    data_list$Sigma    <- extend_length(HCR$Sigma)
    data_list$Fmult    <- extend_length(HCR$Fmult)
    data_list$HCRorder <- extend_length(HCR$HCRorder)
    data_list$QnormHCR <- stats::qnorm(data_list$Pstar, 0, data_list$Sigma)

    # Fill out switches if missing. `quiet_data_check` covers switch_check()'s
    # deprecation messages for the same reason it covers data_check()'s
    # warnings: the diagnostic refits re-validate the caller's own data_list, so
    # a retrospective would otherwise repeat them once per peel.
    if (isTRUE(quiet_data_check)) {
      data_list <- suppressMessages(Rceattle::switch_check(data_list))
    } else {
      data_list <- Rceattle::switch_check(data_list)
    }

    # A Dirichlet-multinomial reads its composition weight as a log, so the
    # workbook default of 1 is a starting weight of e. Reported here rather
    # than inside switch_check(), which build_params() and build_map() also call
    # -- that would print it three times per fit. Reports only; changes nothing.
    if (isTRUE(quiet_data_check)) suppressMessages(.rce_flag_dm_weight_scale(data_list))
    else .rce_flag_dm_weight_scale(data_list)

    # Check for data error. `quiet_data_check` drops the WARNINGS only -- errors
    # still stop the fit. The diagnostic refits set it: they re-validate a
    # data_list the caller already fitted once, so every warning it raises has
    # been seen, and a retrospective or an MSE would otherwise repeat it once per
    # peel or per iteration.
    if (isTRUE(quiet_data_check)) suppressWarnings(data_check(data_list))
    else data_check(data_list)

    # * Pool process linkages into a global table + design matrix ----
    # No-op when no build_*() supplied a `linkages` list. When
    # linkages are present, this materializes each spec against
    # data_list$env_data and unions design columns by name.
    #
    # Linkages consume env_data POSITIONALLY: row r is applied to model year
    # styr + r - 1 (years beyond env_data get a zero offset). So a `Year`
    # column must be sorted, start at styr, and be contiguous, or a covariate /
    # deviate lands on the wrong year. Validate up front rather than misalign.
    .has_linkage <- any(vapply(
      list(data_list$growth_linkages, data_list$M1_linkages,
           data_list$srr_linkages, data_list$q_linkages, data_list$sel_linkages,
           data_list$comp_linkages),
      function(x) !is.null(x) && length(x) > 0L, logical(1)))
    if (.has_linkage) {
      # Prepend/gap-fill env_data to start at styr (NA for missing years) so a
      # later-starting `observe` covariate aligns; the check then validates the
      # extended table. NA years are skipped in the QAR1 observation (masked).
      data_list$env_data <- .extend_env_data(data_list$env_data, data_list$styr)
      .check_env_data_years(data_list$env_data, data_list$styr)
    }
    .message_auto_fleet_linkages(list(q   = data_list$q_linkages,
                                      sel = data_list$sel_linkages))
    .linkage_pool <- pool_linkages(
      spec_groups = list(growth      = data_list$growth_linkages,
                         M           = data_list$M1_linkages,
                         recruitment = data_list$srr_linkages,
                         q           = data_list$q_linkages,
                         sel         = data_list$sel_linkages,
                         comp        = data_list$comp_linkages),
      env_data    = data_list$env_data,
      # The fleet / species levels carry the model's own names, so a spec built
      # with `linkage_spec(fleet = "Shelikof")` can be resolved to an id.
      strata      = list(
        fleet   = .label_strata(seq_len(nrow(data_list$fleet_control)),
                                data_list$fleet_control$Fleet_name),
        species = .label_strata(seq_len(data_list$nspp), data_list$spnames),
        sex     = if (length(data_list$nsex) > 1L &&
                      length(unique(data_list$nsex)) > 1L) {
          stats::setNames(lapply(seq_len(data_list$nspp),
                                 function(sp) seq_len(data_list$nsex[sp])),
                          as.character(seq_len(data_list$nspp)))
        } else {
          seq_len(max(data_list$nsex))
        },
        age_bin = if (length(data_list$nages) > 1L &&
                      length(unique(data_list$nages)) > 1L) {
          stats::setNames(lapply(seq_len(data_list$nspp),
                                function(sp) seq_len(data_list$nages[sp])),
                         as.character(seq_len(data_list$nspp)))
        } else {
          seq_len(max(data_list$nages))
        }
      )
    )
    data_list$linkage_table <- .linkage_pool$table
    data_list$linkage_X     <- .linkage_pool$X

    # A q linkage is applied per fleet and is not reconciled across a shared
    # Catchability_index, so it can break a group data_check() has already
    # passed. Checked here because the table above is what it reads.
    if (!isTRUE(quiet_data_check)) .warn_q_linkage_shared_group(data_list)
    # (Fixed-effect covariates with missing years are rejected earlier, in
    # materialize_linkage(), before model.matrix() can silently drop NA rows.)

    # Selectivity linkages are only consumed by the parametric forms wired in
    # the TMB template. Reject a sel linkage on a fleet whose selectivity form
    # is not yet wired, and the non-parametric `coff` param, so the effect is
    # never silently dropped. (Empirical and the RPM random walk cannot carry
    # a covariate offset at all.)
    .check_sel_linkage_support(data_list$linkage_table, data_list$fleet_control)
    .check_q_linkage_support(data_list$linkage_table, data_list$fleet_control)
    .check_comp_linkage_support(data_list$linkage_table, data_list)

    # Random-effect linkage rows (IID `~ (1 | group)`) are now consumed by the
    # TMB template: beta_linkage_re carries the deviations and jnll_comp row 20
    # the N(0, sigma) density. Correlated structures (rw()/ar1()) are still
    # rejected upstream in .materialize_re_design() until their densities land.

    # * DSEM ----
    # A dynamic structural equation model on the recruitment deviations: the
    # recdevs become latent states of a GMRF, so environmental covariates and
    # lagged/simultaneous paths drive recruitment. Built after the recruitment
    # switches (build_dsem_objects() reads random_rec and proj_mean_rec) and
    # after the linkage pool, keeping the structure specs contiguous.
    #
    # dsem = NULL is the default and must stay completely inert: with no DSEM
    # attached, every path below is skipped and the model is textually the
    # non-DSEM model. That is what lets /golden-check gate this change.
    if (is.null(dsem)) {
      # A DSEM fit's own data_list carries its specification, and re-fitting
      # from that data_list with dsem = NULL is the ordinary way to compare a
      # model with and without one. Clear the specification, or the fit reports
      # a DSEM it does not have: every diagnostic guard reads dsem_settings, and
      # save_config() would write out a `dsem:` block that load_config() then
      # replays as a different model (measured at 86.2 nll on BS2017SS). This is
      # the data_list counterpart of the `inits` scrub further down.
      data_list$dsem_settings <- NULL
      data_list$model_config$dsem <- NULL
    } else {
      # The build_*() specs are plain unclassed lists in this package, so tell a
      # specification from already-built objects structurally: built objects
      # carry $tmb_inputs, a spec carries the build_DSEM() fields. Anything else
      # is a mistake worth naming -- without this, a bare list with a $sem
      # element would be silently accepted and its family / estimate_projection
      # quietly defaulted. Checked first so a bad `dsem` is named here rather
      # than by a base R error a few lines further down.
      .dsem_built <- is.list(dsem) && !is.null(dsem$tmb_inputs)
      .dsem_spec  <- is.list(dsem) && "sem" %in% names(dsem)
      if (!.dsem_built && !.dsem_spec) {
        stop("`dsem` must be a specification from build_DSEM() (or built ",
             "objects from build_dsem_objects()); got ",
             paste(class(dsem), collapse = "/"), ".", call. = FALSE)
      }

      # estimate_projection = TRUE says the SEM supplies projected recruitment,
      # so it overrides proj_mean_rec. The two are otherwise contradictory:
      # under mean-recruitment projection the model sets R = avg_R past
      # endyr and the estimated projection states reach only the dynamic B0/BF
      # reference series, so the SEM would be fitted over the projection and
      # then ignored there. Overriding before the DSEM is built keeps one
      # answer -- build_dsem_objects() reads proj_mean_rec too.
      .ep <- if (.dsem_built) {
        isTRUE(dsem$dsem_settings$estimate_projection)
      } else isTRUE(dsem$estimate_projection)
      if (.ep && isTRUE(as.logical(data_list$proj_mean_rec))) {
        if (verbose > 0) {
          message("DSEM `estimate_projection = TRUE`: projecting recruitment ",
                  "through the SEM rather than at mean recruitment ",
                  "(`proj_mean_rec` set to FALSE). Use ",
                  "`estimate_projection = FALSE` to keep mean-recruitment ",
                  "projection.")
        }
        data_list$proj_mean_rec <- 0L
      }

      if (.dsem_built) {
        mod_objects$dsem <- dsem
        # Record the specification alongside the built objects. .refit_like()
        # and the diagnostic guards read data_list$dsem_settings, so a fit built
        # from pre-built objects would otherwise look DSEM-free to them.
        data_list$dsem_settings <- dsem$dsem_settings
      } else {
        data_list$dsem_settings <- dsem
        mod_objects$dsem <- build_dsem_objects(
          dsem_settings = dsem,
          debug         = estimateMode %in% c(2, 4),   # match build_map()
          data_list     = data_list)
        data_list$dsem_settings$sem <- mod_objects$dsem$sem
      }

      # DSEM and a recruitment linkage both claim the recruitment deviations.
      # Reject the combination up front rather than silently applying both.
      # linkage_table$process holds the process NAME (see LINKAGE_COLS), not the
      # integer code -- comparing against LINKAGE_PROCESS_CODES coerces to "0"
      # and never matches, which would let the two structures stack silently
      # once the C++ lands.
      if (!is.null(data_list$linkage_table) &&
          any(data_list$linkage_table$process == "recruitment")) {
        stop("A DSEM and a recruitment linkage both structure the recruitment ",
             "deviations; they cannot be combined. Supply either `dsem` or ",
             "`recFun = build_srr(linkages = ...)`, not both.", call. = FALSE)
      }

    }

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # 2: Load/build parameters ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    if (is.character(inits) | is.null(inits)) {
      start_par <- suppressWarnings(Rceattle::build_params(data_list = data_list))
    } else {
      start_par <- inits

      # Guard: catch `inits` that would make TMB::MakeADFun() segfault in
      # getParameterOrder() instead of raising an R error, by comparing the
      # supplied `inits` to the parameter template `build_params(data_list)`
      # implies. Flagged when a declared parameter is absent, or when any
      # parameter's length differs. Both are fatal: build_map() keeps the map
      # consistent with `inits`, but the C++ still reads several parameters at
      # data-driven bounds (linkage-table sizes for `*_linkage`; year bounds
      # nyrs_hind / nyrs_proj for e.g. `index_q_dev`, `log_M1_dev`, `rec_dev`),
      # so a too-short array reads out of bounds. The usual causes are `inits`
      # from a fit that used build_catchability / build_selectivity /
      # build_composition without the matching *Fun re-supplied (e.g. an older
      # .refit_like()), or a warm start not extended to a later `endyr`.
      .skel    <- suppressWarnings(Rceattle::build_params(data_list = data_list))
      .missing <- setdiff(names(.skel), names(start_par))
      .shared  <- intersect(names(.skel), names(start_par))
      .badlen  <- .shared[vapply(.shared, function(nm)
        length(unlist(start_par[[nm]])) != length(unlist(.skel[[nm]])), logical(1))]
      if (length(.missing) || length(.badlen)) {
        stop("`inits` do not match the parameters implied by `data_list`; ",
             "TMB would crash on this. ",
             if (length(.missing))
               paste0("Missing: ", paste(.missing, collapse = ", "), ". ") else "",
             if (length(.badlen))
               paste0("Length mismatch: ",
                      paste(sprintf("%s (inits %d, model %d)", .badlen,
                        vapply(.badlen, function(nm) length(unlist(start_par[[nm]])), 1L),
                        vapply(.badlen, function(nm) length(unlist(.skel[[nm]])), 1L)),
                      collapse = "; "), ". ") else "",
             "Rebuild `inits` for this `data_list` (inits = NULL builds them ",
             "fresh); if these came from a fit at an earlier `endyr`, re-fit to ",
             "extend the year-varying parameters, and if that fit used ",
             "build_catchability() / build_selectivity() / build_composition(), ",
             "pass the matching qFun / selFun / compFun.", call. = FALSE)
      }
      # Drop anything the model no longer has a parameter for, in skeleton
      # order. `inits` from an older fit can carry a retired block (e.g.
      # log_growth_par_devs, removed in 5.9.0). MakeADFun drops names the
      # template does not declare, but build_map() runs on start_par first and
      # produces a map entry per element, so the stale name reached MakeADFun
      # through the map and stopped with "Names in map must correspond to
      # parameter names" -- a refit that failed on a parameter the model no
      # longer has. Extra names are dropped silently because they are inert by
      # definition; a MISSING one is the error above.
      # ...EXCEPT the DSEM blocks, which build_params() does not declare (they
      # come from build_dsem_objects() and are merged in below), so a bare
      # start_par[names(.skel)] would drop them. They are not inert: dropped,
      # merge_dsem_params() refills them from the freshly built template, i.e.
      # the START values, discarding the warm start -- which matters wherever
      # the DSEM is PINNED rather than estimated (retrospective()'s forecast
      # refit, estimateMode = 2).
      # .DSEM_PARAM_NAMES, not grep("^dsem_"), so this stays the exact
      # complement of the scrub below that removes these blocks when the fit has
      # no DSEM.
      start_par <- start_par[c(names(.skel),
                               intersect(names(start_par), .DSEM_PARAM_NAMES))]
      rm(.skel, .missing, .shared, .badlen)

      # Set F for years with 0 catch to very low number
      zero_catch <- data_list$catch_data |>
        dplyr::filter(.data$Year <= data_list$endyr &
                        .data$Catch == 0) |>
        dplyr::mutate(Year = Year - data_list$styr + 1) |>
        dplyr::select("Fleet_code", "Year") |>
        as.matrix()
      start_par$log_F[zero_catch] <- -999
      rm(zero_catch)

      # Update proj F prop
      start_par$proj_F_prop <- data_list$fleet_control$Proj_F_proportion
    }

    # Drop DSEM parameters from a warm start when this fit has no DSEM. They are
    # not in build_params(), so the shape guard above passes them through as
    # unflagged extras; build_map() then propagates the extra names into
    # mapFactor and TMB rejects the map ("Names in map must correspond to
    # parameter names"). Reached whenever `inits` come from a DSEM fit
    # (including one produced on the dev-DSEM branch) and `dsem` is left NULL.
    if (is.null(mod_objects$dsem) && !is.null(start_par)) {
      .stale_dsem <- intersect(names(start_par), .DSEM_PARAM_NAMES)
      if (length(.stale_dsem) > 0) {
        if (verbose > 0) {
          message("Dropping DSEM parameters from `inits` (no `dsem` supplied): ",
                  paste(.stale_dsem, collapse = ", "))
        }
        start_par[.stale_dsem] <- NULL
      }
      rm(.stale_dsem)
    }

    # DSEM parameters. Added after build_params() rather than inside it, so they
    # cannot collide with the "(Intercept)" init push that runs at the end of
    # that function. Zero-length when there is no DSEM, which is what makes the
    # non-DSEM objective identical.
    # Add the DSEM blocks that `inits` does not already carry, and leave the ones
    # it does alone. Appending with c() duplicated them for a warm start from a
    # DSEM fit ("Map and parameter objects are not the same size"); overwriting
    # them unconditionally was worse -- it silently reset x_tj to 0 and beta_z to
    # its start value, so estimateMode = 2 reprojected from a different model
    # than the one it was handed (measured: 47% SSB error, no warning).
    # merge_dsem_params() has exactly the right semantics.
    start_par <- merge_dsem_params(start_par, .dsem_par_shape(mod_objects$dsem))

    # A warm start whose DSEM blocks do not match this model's shape is a
    # different model, not a starting value -- e.g. a retrospective peel that
    # shortens nyrs_dsem while inits carries the full-length x_tj. Fail rather
    # than let TMB reinterpret it.
    start_par <- .dsem_check_shape(start_par, mod_objects$dsem)

    if (verbose > 0) { message("Step 1: Parameter build complete") }


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # 3: Load/build map ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    if (is.null(map)) {
      # build_map() warns about configurations the caller asked for and will not
      # get: an M1 model wanting sex-specific mortality on a single-sex species,
      # a Time_varying_sel the selectivity form ignores, a Laplace request the
      # map cannot honour. Each one changes what is estimated, so it has to be
      # seen. They are de-duplicated rather than passed straight through because
      # build_map() raises the selectivity ones once per fleet and per sex, and
      # because .refit_like() re-enters fit_mod() once per retrospective peel,
      # jitter and MSE simulation -- unfiltered, one real warning becomes
      # hundreds of identical lines and stops being read.
      seen <- character(0)
      map <- withCallingHandlers(
        build_map(data_list, start_par,
                  debug = estimateMode %in% c(2, 4), # turn off hindcast parameters in projection / debug mode
                  random_rec = random_rec, random_sel = random_sel,
                  random_q = random_q),
        warning = function(w) {
          msg <- conditionMessage(w)
          if (msg %in% seen) invokeRestart("muffleWarning")
          seen <<- c(seen, msg)
        }
      )

      # DSEM map entries, if any. Assigned by name rather than concatenated:
      # this block can run more than once in a fit, and c() would append a second
      # copy of every DSEM block. The bounds loop iterates mapFactor, so a
      # duplicate silently produces more bounds than obj$par has parameters.
      .dsem_map <- .dsem_map_entries(mod_objects$dsem)
      for (nm in names(.dsem_map$mapList))   map$mapList[[nm]]   <- .dsem_map$mapList[[nm]]
      for (nm in names(.dsem_map$mapFactor)) map$mapFactor[[nm]] <- .dsem_map$mapFactor[[nm]]

      # Under a DSEM the model OVERWRITES rec_dev from the latent states and
      # R_sd from beta_z, so both blocks have exactly zero gradient. Left
      # estimated they would be free parameters that cannot move the objective:
      # flat directions, a singular Hessian, and no standard errors. Map them
      # out. (They stay declared as parameters either way -- that is what keeps
      # the non-DSEM path textually unchanged.)
      if (!is.null(mod_objects$dsem)) {
        for (nm in c("rec_dev", "R_log_sd")) {
          map$mapList[[nm]][]   <- NA
          map$mapFactor[[nm]]   <- factor(rep(NA, length(map$mapList[[nm]])))
        }
      }
    }
    if (verbose > 0) { message("Step 2: Map build complete") }

    # Now that the map is known, drop any composition-weighting prior whose DM
    # weight is fixed in this configuration -- it would only add a constant.
    data_list$linkage_table <- .neutralize_inert_comp_priors(
      data_list$linkage_table, map, verbose = verbose)


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # 4: Get bounds ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # A DSEM standard deviation is a Cholesky diagonal, so its sign is not
    # identified and build_bounds() bounds it at 0 unless the DSEM was built
    # with bound_sd = FALSE. A model carried forward from before that bound can
    # sit at the negative root, which would put a warm start outside its own
    # bound -- and build_bounds() rejects that outright, so the flip has to
    # happen BEFORE it rather than against the bounds it returns. Flipping is
    # exact where the bound applies: the variable has no other two-headed path,
    # so the likelihood sees only the square.
    #
    # Gated on bound_sd, so `FALSE` leaves a start the caller chose where it is.
    # Rewriting it would defeat the one job that argument has -- reproducing an
    # old fit parameter for parameter, or showing the two optima are mirrored.
    if (!is.null(mod_objects$dsem) && !is.null(start_par$dsem_beta_z) &&
        .dsem_bound_sd(mod_objects$dsem)) {
      .sdz <- .dsem_sd_indices(mod_objects$dsem)
      .k <- .sdz$sd[.sdz$sd <= length(start_par$dsem_beta_z)]
      # which(), not a bare logical: an NA start value would index with NA and
      # write an NA back into the parameter list.
      .k <- .k[which(start_par$dsem_beta_z[.k] < 0)]
      if (length(.k)) {
        start_par$dsem_beta_z[.k] <- abs(start_par$dsem_beta_z[.k])
        if (verbose > 0) {
          message("Flipped ", length(.k), " DSEM standard deviation(s) to the ",
                  "positive root; the sign is not identified and is bounded at 0.")
        }
      }
      if (length(.sdz$entangled) && verbose > 0) {
        message("DSEM variable(s) ", paste(.sdz$entangled, collapse = ", "),
                " carry a cross-covariance or lagged variance path, so their ",
                "standard deviations keep an unidentified sign and no bound. ",
                "See ?build_bounds.")
      }
    }
    if (is.null(bounds)) {
      bounds <- Rceattle::build_bounds(param_list = start_par, data_list,
                                       dsem = mod_objects$dsem)
    }
    if (verbose > 0) { message("Step 3: Parameter bounds complete") }


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # 5: Setup random effects ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # Turns on laplace approximation.
    #
    # What the two "is it free?" guards below are and are NOT for. They are NOT
    # protecting MakeADFun: TMB resolves `random` AFTER applying `map`, so a
    # block whose map is entirely NA contributes no indices and TMB drops it
    # (`random <- NULL`). Naming a fully-mapped parameter is a no-op there --
    # verified: `random = "init_dev"` on an initMode = "Equilibrium" fit gives
    # the same objective and gradient, and the same zero-length random vector,
    # as `random = NULL`. It does not yield NaN.
    #
    # What DOES read this vector's emptiness is TMBphase() (R/6-phaser.R): when
    # `length(random) > 0` it pins the twelve RE variance / correlation
    # hyperparameters to NA in every phase. So whether a fully-mapped block is
    # listed here decides how a phased fit treats those hyperparameters, and
    # that -- not a NaN -- is why the guards are kept.
    #
    # build_map() returns a list with $mapList (the raw per-parameter maps,
    # NA = fixed); that is what the guards read.
    random_vars <- c()

    # Under a DSEM the recruitment deviations are DERIVED from the latent states,
    # so it is dsem_x_tj that gets integrated out, not rec_dev -- rec_dev is
    # overwritten in the model and carries no density of its own. Appended,
    # never assigned over the top: doing the latter would silently drop
    # init_dev, index_q_dev, log_M1_dev and beta_linkage_re from the Laplace
    # approximation.
    if (!is.null(mod_objects$dsem)) {
      # Take the random declaration from the DSEM itself rather than naming
      # x_tj here: dsem::dsem() decides which blocks are integrated out (it
      # declares mu_j random alongside x_tj), and hardcoding one of them
      # estimates the other as a fixed effect -- a different model, which shows
      # up as a non-positive-definite Hessian rather than an error.
      .dsem_re <- paste0("dsem_", mod_objects$dsem$tmb_inputs$random)
      for (nm in .dsem_re) {
        if (any(!is.na(map$mapList[[nm]]))) random_vars <- c(random_vars, nm)
      }
    }

    if (random_rec) {
      if (is.null(mod_objects$dsem)) random_vars <- c(random_vars, "rec_dev")
      # init_dev is integrated out only where it carries a density. The initial
      # deviate penalty dnorm(init_dev, -sigma^2/2, R_sd) applies when
      # initMode > 1 and != 5 (ceattle.cpp, JNLL_INIT_DEV); keep this in lockstep
      # with that gate. Under "FreeParams" the initial age structure is estimated
      # as fixed effects, and the equilibrium modes hold init_dev at 0.
      #
      # initMode is a canonical string here, so resolve it to its integer code
      # before comparing -- a character/numeric comparison sorts lexicographically.
      init_mode_int <- if (is.character(data_list$initMode)) {
        unname(initMode_map[data_list$initMode])
      } else {
        as.integer(data_list$initMode)
      }
      init_dev_has_density <- !is.na(init_mode_int) &&
        init_mode_int > 1 && init_mode_int != 5
      if (init_dev_has_density && any(!is.na(map$mapList$init_dev))) {
        random_vars <- c(random_vars, "init_dev")
      }
    }
    if (random_q) {
      random_vars <- c(random_vars, "index_q_dev")
    }
    if (random_sel) {
      # Selectivity deviates are integrated out only where the model scores
      # them. The densities on log_sel_slp_dev / sel_inf_dev are gated on
      # Time_varying_sel being IID, AR1, RandomWalk or RandomWalkAscending
      # (ceattle.cpp, JNLL_SEL_DEV). "Block" is deliberately unscored -- one
      # fixed parameter per block, no penalty -- and sel_dev_log_sd is mapped
      # out for it too, so there is not even a variance to estimate. Integrating
      # a flat block against nothing is an improper random effect: the marginal
      # objective is NaN, and TMB reports it as "NA/NaN gradient evaluation",
      # which names neither selectivity nor random_sel.
      blocked <- .rce_unscored_sel_dev_fleets(data_list$fleet_control, map$mapList)
      if (length(blocked)) {
        stop("Cannot use `random_sel = TRUE` with `Time_varying_sel = \"Block\"` ",
             "on fleet(s) ", paste(blocked, collapse = ", "), ". Selectivity ",
             "blocks are fixed effects carrying no penalty, so there is no ",
             "density to integrate them against. Either set ",
             "`Time_varying_sel` to \"IID\" or \"RandomWalk\" on those fleets, ",
             "or fit with `random_sel = FALSE` (blocks are then estimated as ",
             "the fixed effects they are).")
      }

      # Listed whether or not the maps are empty. TMB drops a fully-mapped name
      # itself, but TMBphase() reads `length(random) > 0` to decide whether to
      # pin the RE variance hyperparameters in every phase (see the note above),
      # so dropping an empty array here would change a phased fit rather than
      # tidy it.
      random_vars <- c(random_vars, "log_sel_slp_dev", "sel_inf_dev", "sel_coff_dev")
    }
    if (sum(data_list$M1_re) > 0) {
      random_vars <- c(random_vars, "log_M1_dev")
    }
    # Random-effect linkage deviations enter the Laplace approximation. Unlike
    # the init_dev guard above this one is load-bearing: a linkage random effect
    # can be the ONLY random effect, so listing a fully-mapped `beta_linkage_re`
    # would make `random_vars` non-empty and send a phased fit down TMBphase()'s
    # hyperparameter-pinning branch for a model that has no random effects at
    # all. Changing it therefore changes results -- it is not cosmetic.
    # beta_linkage_re_pen is deliberately absent: it holds the deviations a spec
    # asked to keep OUT of the Laplace approximation (integrate = FALSE), so it
    # is an ordinary fixed effect carrying a penalty. Adding it here would
    # integrate exactly what the user asked not to integrate. A purely penalized
    # model leaves beta_linkage_re length 0, so the guard below is FALSE and
    # `random` is correctly empty.
    if (any(!is.na(map$mapList$beta_linkage_re))) {
      random_vars <- c(random_vars, "beta_linkage_re")
    }

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # 6: Reorganize data ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    if (!is.null(TMBfilename)) {
      TMB::compile(paste0(TMBfilename, ".cpp"),
                   PKG_CXXFLAGS = "",
                   framework = "TMBad",
                   safebounds = FALSE, safeunload = FALSE)
      dyn.load(TMB::dynlib(TMBfilename))
      TMBfilename <- basename(TMBfilename)
    }
    if (is.null(TMBfilename)) {
      TMBfilename <- "ceattle"
    }

    # Composition proportion offset. It lives on `data_list` so that every
    # internal re-fit (projections, retrospective, jitter, run_mse, ...) inherits
    # the same value without threading it through each fit_mod() call. An explicit
    # fit_control(comp_offset=) overrides. Defaults to 1e-5.
    if (!is.null(fit_control$comp_offset)) data_list$comp_offset <- fit_control$comp_offset
    if (is.null(data_list$comp_offset))    data_list$comp_offset <- 1e-5

    # Bias-adjustment flags for lognormal likelihoods (default TRUE/1). Stored on
    # data_list as numeric (0/1) so the DATA_SCALAR is a clean double, so the TMB
    # template can toggle the -sigma^2/2 correction.
    data_list$bias_adjust_obs <- data_list$bias_adjust_proc <- 1
    # TODO(review): a user-supplied NA (e.g. fit_control(bias_adjust_obs = NA))
    # becomes NaN via as.numeric() and propagates a NaN objective. Add a finite
    # check to reject it with a clear error.
    if(!is.null(fit_control$bias_adjust_obs)) data_list$bias_adjust_obs <- as.numeric(fit_control$bias_adjust_obs)
    if(!is.null(fit_control$bias_adjust_proc)) data_list$bias_adjust_proc <- as.numeric(fit_control$bias_adjust_proc)

    # Carried like the bias-adjustment flags above, so a refit can recover how
    # the model was fitted rather than take fit_control()'s default.
    data_list$projection_uncertainty <- isTRUE(projection_uncertainty)

    # Reorganize data for .cpp file. The OSA observation vector is built on
    # demand by osa_residuals(), so only the fast aggregate metadata is assembled
    # here (build_osa = FALSE, the rearrange_data() default).
    data_list_reorganized <- Rceattle::rearrange_data(data_list)
    data_list_reorganized <- c(list(model = TMBfilename), data_list_reorganized)

    # DSEM inputs. The template declares these unconditionally, so every fit has
    # to supply them; with no DSEM they are zero-length and dsem_on = 0.
    .dsem_in <- if (is.null(mod_objects$dsem)) .dsem_null_inputs() else
      .dsem_tmb_inputs(mod_objects$dsem)
    data_list_reorganized <- c(data_list_reorganized, .dsem_in$data)
    data_list_reorganized$forecast <- rep(0, data_list_reorganized$nspp) # hindcast switch

    # Scrub fields that TMB's dataSanitize cannot recurse into (data
    # frames with character columns, list-of-spec objects, character
    # vectors). They live on data_list for downstream R code; they
    # are re-injected below in TMB-friendly form.
    data_list_reorganized$linkage_table   <- NULL
    data_list_reorganized$growth_linkages <- NULL
    data_list_reorganized$growth_fun      <- NULL
    data_list_reorganized$M1_linkages     <- NULL
    data_list_reorganized$srr_linkages    <- NULL
    data_list_reorganized$q_linkages      <- NULL
    data_list_reorganized$sel_linkages    <- NULL
    data_list_reorganized$comp_linkages   <- NULL

    # OSA residual metadata: obs_ctl maps each obsvec element back to its
    # fleet/species/year/age. It is an R-side data frame (TMB's dataSanitize
    # cannot recurse into it), so stash it for the returned object and scrub it
    # from the TMB data list. The numeric obsvec / *_obsvec_idx / osa_mode
    # fields are TMB-friendly and stay.
    .obs_ctl <- data_list_reorganized$obs_ctl
    data_list_reorganized$obs_ctl <- NULL

    # * Inject linkage-table encoding into the TMB DATA ----
    # Empty when no build_*() supplied a `linkages` list. TMB's
    # DATA_MATRIX can be touchy about 0-dim shapes during ADFun
    # construction, so we substitute a single-cell sentinel when
    # there are no linkages -- the cpp loop short-circuits on
    # n_linkage == 0 without ever indexing into `linkage_X`.
    .linkage_table_now <- data_list$linkage_table %||% new_linkage_table()
    .linkage_X_now     <- if (nrow(.linkage_table_now) > 0L) {
      data_list$linkage_X
    } else {
      matrix(0, nrow = 1L, ncol = 1L)
    }
    .linkage_enc <- encode_linkage_for_tmb(
      table = .linkage_table_now,
      X     = .linkage_X_now
    )
    # The TMB cpp reads these by exact name; do not rename.
    data_list_reorganized$linkage_process      <- .linkage_enc$linkage_process
    data_list_reorganized$linkage_param        <- .linkage_enc$linkage_param
    data_list_reorganized$linkage_species      <- .linkage_enc$linkage_species
    data_list_reorganized$linkage_sex          <- .linkage_enc$linkage_sex
    data_list_reorganized$linkage_age_bin      <- .linkage_enc$linkage_age_bin
    data_list_reorganized$linkage_fleet        <- .linkage_enc$linkage_fleet
    data_list_reorganized$linkage_X_col        <- .linkage_enc$linkage_X_col
    data_list_reorganized$linkage_link         <- .linkage_enc$linkage_link
    data_list_reorganized$linkage_re_index     <- .linkage_enc$linkage_re_index
    data_list_reorganized$linkage_re_sigma      <- .linkage_enc$linkage_re_sigma
    data_list_reorganized$linkage_re_integrate  <- .linkage_enc$linkage_re_integrate
    data_list_reorganized$linkage_re_slot       <- .linkage_enc$linkage_re_slot
    data_list_reorganized$linkage_re_struct     <- .linkage_enc$linkage_re_struct
    data_list_reorganized$linkage_re_rho        <- .linkage_enc$linkage_re_rho
    data_list_reorganized$linkage_re_sigma_prior_family <- .linkage_enc$linkage_re_sigma_prior_family
    data_list_reorganized$linkage_re_sigma_prior_p1     <- .linkage_enc$linkage_re_sigma_prior_p1
    data_list_reorganized$linkage_re_sigma_prior_p2     <- .linkage_enc$linkage_re_sigma_prior_p2
    data_list_reorganized$linkage_re_rho_prior_family   <- .linkage_enc$linkage_re_rho_prior_family
    data_list_reorganized$linkage_re_rho_prior_p1       <- .linkage_enc$linkage_re_rho_prior_p1
    data_list_reorganized$linkage_re_rho_prior_p2       <- .linkage_enc$linkage_re_rho_prior_p2
    data_list_reorganized$linkage_re_obs        <- .linkage_enc$linkage_re_obs
    data_list_reorganized$linkage_re_obs_sd     <- .linkage_enc$linkage_re_obs_sd
    data_list_reorganized$linkage_re_obs_value  <- .linkage_enc$linkage_re_obs_value
    data_list_reorganized$linkage_re_obs_mask   <- .linkage_enc$linkage_re_obs_mask
    data_list_reorganized$linkage_is_intercept <- .linkage_enc$linkage_is_intercept
    data_list_reorganized$linkage_prior_family <- .linkage_enc$linkage_prior_family
    data_list_reorganized$linkage_prior_p1     <- .linkage_enc$linkage_prior_p1
    data_list_reorganized$linkage_prior_p2     <- .linkage_enc$linkage_prior_p2
    data_list_reorganized$linkage_X            <- .linkage_enc$linkage_X

    # The composition weights are estimated parameters -- the Dirichlet-
    # multinomial likelihood fits them, which is why they can carry a prior --
    # so they warm-start from `inits` like every other parameter. Their
    # fleet_control columns are their STARTING values and are read by
    # build_params() on a fresh build; re-reading them here would discard the
    # estimate a refit was handed, which is not what a starting value means.
    # A reweighting workflow therefore updates the parameter, not the column.
    #
    # That makes editing a column and re-fitting from an existing fit a silent
    # no-op, which is how composition weights were tuned by hand for years. Warn
    # where the edit cannot possibly take effect: a weight the map holds fixed
    # under a multinomial likelihood is a pure data-weighting input, so `inits`
    # and the column disagreeing means the column is being ignored and nothing
    # in the fit will show it. A Dirichlet-multinomial weight is deliberately
    # excluded -- the likelihood estimates it, so a refit's `inits` differing
    # from the column is the normal, correct state and warning on it would fire
    # on every diagnostic refit.
    if (!is.null(inits)) {
      .weight_blocks <- list(
        comp_weights      = list(col = data_list$fleet_control$Comp_weights,
                                 dist = data_list$fleet_control$Comp_distribution,
                                 label = "fleet_control$Comp_weights"),
        caal_weights      = list(col = data_list$fleet_control$CAAL_weights,
                                 dist = data_list$fleet_control$CAAL_distribution,
                                 label = "fleet_control$CAAL_weights"),
        diet_comp_weights = list(col = data_list$Diet_comp_weights,
                                 dist = data_list$Diet_distribution,
                                 label = "Diet_comp_weights"))
      for (.nm in names(.weight_blocks)) {
        .b   <- .weight_blocks[[.nm]]
        .par <- start_par[[.nm]]
        .fix <- map$mapList[[.nm]]
        if (is.null(.b$col) || is.null(.par) || is.null(.fix)) next
        if (length(.b$col) != length(.par) || length(.fix) != length(.par)) next
        # NA == NA counts as agreement; only a real numeric difference is a
        # discarded edit.
        # Coerce only a non-numeric column (read.csv with stringsAsFactors can
        # hand back a factor, which would compare as NA and poison `if()`
        # below). A numeric column is compared as-is: routing it through
        # as.character() would round it to 15 significant digits on one side of
        # the comparison only, making a value that IS equal read as different.
        .col <- .b$col
        if (is.factor(.col) || is.character(.col)) {
          .col <- suppressWarnings(as.numeric(as.character(.col)))
        }
        .differs <- !is.na(.col) & !is.na(.par) & .col != as.numeric(.par)
        # A Dirichlet-multinomial weight is exempt. The likelihood estimates it,
        # so a fit's parameter having moved away from its column is the normal
        # state -- and under the debug map run_mse() / retrospective() supply,
        # every weight reads as fixed, so without this the warning would fire
        # once per fleet per assessment year of every simulation.
        .stuck <- .differs & is.na(.fix) &
          !as.character(.b$dist) %in% c("1", "DirichletMultinomial")
        if (isTRUE(any(.stuck, na.rm = TRUE))) {
          warning("`", .b$label, "` differs from the supplied `inits$", .nm,
                  "` at ", if (identical(.nm, "diet_comp_weights")) "species " else "fleet ",
                  paste(which(.stuck), collapse = ", "),
                  ", and the fit reads `inits`. The column is only read when a ",
                  "model is built from scratch (`inits = NULL`), so this edit ",
                  "has no effect. Set `inits$", .nm, "` instead, or see ",
                  "?reweight_comps.", call. = FALSE)
        }
      }
      rm(.weight_blocks)
    }
    # Proportion of projected F to each fleet
    start_par$proj_F_prop <- data_list$fleet_control$Proj_F_proportion
    # Fixed fishing mortality for projections for each species
    if (!is.null(HCR$Ftarget) & HCR$HCR %in% c(2, "ConstantF")) {
      start_par$log_Ftarget <- log(HCR$Ftarget)
    }

    # Update M1 parameter object from data if initial parameter values input
    if (updateM1) {
      m1 <- array(0, dim = c(data_list$nspp,
                             max(data_list$nsex, na.rm = TRUE),
                             max(data_list$nages, na.rm = TRUE)))

      for (i in 1:nrow(data_list$M1_base)) {
        sp  <- as.numeric(as.character(data_list$M1_base$Species[i]))
        sex <- as.numeric(as.character(data_list$M1_base$Sex[i]))

        # Handle sex == 0 case for 2-sex species
        sex_values <- if (sex == 0) 1:data_list$nsex[sp] else sex

        for (j in 1:length(sex_values)) {
          m1[sp, sex_values[j], 1:max(data_list$nages, na.rm = TRUE)] <- as.numeric(data_list$M1_base[i, (1:max(data_list$nages, na.rm = TRUE)) + 2])
        }
      }
      start_par$log_M1 <- log(m1)
    }

    # Fix alpha at the prior mean, where srr_prior is an alpha rather than a
    # steepness (see .srr_prior_is_alpha). build_params() applies the same rule;
    # this covers the case where the caller supplied `inits` instead.
    if (data_list$srr_est_mode %in% c(0, 2) & data_list$srr_pred_fun > 1 &&
        .srr_prior_is_alpha(data_list) && !is.null(data_list$srr_prior)) {
      start_par$rec_pars[, 2] <- log(data_list$srr_prior)
    }

    # Explicit stock-recruit starting values from build_srr(), which override
    # both the defaults and the prior mean.
    if (!is.null(data_list$srr_alpha_init)) {
      start_par$rec_pars[, 2] <- log(data_list$srr_alpha_init)
    }
    if (!is.null(data_list$srr_beta_init)) {
      start_par$rec_pars[, 3] <- log(data_list$srr_beta_init)
    }

    # Starting parameters as the model uses them: the blocks above set
    # proj_F_prop, log_Ftarget, log_M1 and the stock-recruit alpha / beta.
    # retrospective() and jitter() reuse this as their refit starting values.
    # Taken before TMBphase() replaces start_par with a fitted state.
    mod_objects$initial_params <- start_par

    if (verbose > 0) { message("Step 4: Data rearrange complete") }


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # 7: Set up parameter bounds ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    L <- c()
    U <- c()
    L_block <- c()   # which parameter each bound belongs to, for the reorder below
    for (i in 1:length(map$mapFactor)) {
      nm <- names(map$mapFactor)[i]
      if (!nm %in% random_vars) { # no bounds for random effects
        mf   <- unlist(map$mapFactor[[i]])
        if (!is.factor(mf)) mf <- factor(mf)
        # One bound per distinct map index, in FACTOR-LEVEL order -- the order
        # TMB collapses a mapped block in (TMB:::updateMap uses tapply over the
        # map factor). First-occurrence order diverges wherever mirrored fleets
        # share indices out of order (adjust_map_shared_params), swapping their
        # bounds.
        keep <- match(levels(droplevels(mf)), as.character(mf))
        L <- c(L, unlist(bounds$lower[[nm]])[keep])
        U <- c(U, unlist(bounds$upper[[nm]])[keep])
        L_block <- c(L_block, rep(nm, length(keep)))
      }
    }

    # Dimension check
    start_par <- start_par[names(map$mapFactor), drop = FALSE]
    dim_check <- sapply(start_par, function(x) length(unlist(x))) == sapply(map$mapFactor, function(x) length(unlist(x)))
    if (sum(dim_check) != length(dim_check)) {
      stop(paste0("Map and parameter objects are not the same size for: ", names(dim_check)[which(dim_check == FALSE)]))
    }


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # 8: Phase hindcast ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    if (is.logical(phase)) {
      if (phase) {
        phaseList <- set_phases()
      }
    }

    if (!is.logical(phase)) {
      warning("Using input phase. Please set phase = TRUE if using defaults.")
      phaseList <- phase
      phase     <- TRUE
    }


    step <- 5
    if (phase & estimateMode %in% c(0, 1)) {
      if (verbose > 0) { message(paste0("Step ", step, ": Phasing begin")) }; step <- step + 1
      phase_pars <- Rceattle::TMBphase(
        data         = data_list_reorganized,
        parameters   = start_par,
        map          = map$mapFactor,
        random       = random_vars,
        phases       = phaseList,
        model_name   = TMBfilename,
        silent       = verbose != 2,
        use_gradient = use_gradient,
        control      = control
      )

      # Pull the per-phase convergence log (attached by TMBphase) for the
      # phasing diagnostic, then strip it so it doesn't ride along in start_par.
      mod_objects$.conv_phase <- attr(phase_pars, "phase_log")
      attr(phase_pars, "phase_log") <- NULL

      mod_objects$phase_params <- phase_pars
      start_par <- phase_pars

      if (verbose > 0) { message(paste0("Step ", step, ": Phasing complete")) }
      step <- step + 1
    }


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # 9: Fit hindcast ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # * Build ----
    if (sum(as.numeric(unlist(map$mapFactor)), na.rm = TRUE) == 0) { stop("Map of length 0: all NAs") }
    obj <- .rce_make_adfun(
      data_list_reorganized,
      parameters = start_par,
      DLL        = TMBfilename,
      map        = map$mapFactor,
      random     = random_vars,
      silent     = verbose != 2
    )

    # Align the bounds with obj$par. They were assembled in build_params() order,
    # but TMB orders obj$par by the sequence the PARAMETER_* macros appear in the
    # template -- and the two disagree: build_params() lists the linkage
    # coefficients after log_F, ceattle.cpp declares them before it. The lengths
    # match either way, so a mismatch is silent: box constraints simply land on
    # the wrong parameters (log_F losing its upper rail to a linkage coefficient,
    # say). Reorder blockwise, then assert, because getting this wrong is
    # undetectable downstream.
    if (length(L)) {
      blocks <- rle(names(obj$par))$values
      stopifnot(setequal(blocks, unique(L_block)))
      L <- unlist(split(L, factor(L_block, levels = blocks))[blocks], use.names = FALSE)
      U <- unlist(split(U, factor(L_block, levels = blocks))[blocks], use.names = FALSE)
      if (length(L) != length(obj$par)) {
        stop("Parameter bounds could not be aligned with obj$par (",
             length(L), " vs ", length(obj$par), ").")
      }
      # A short bounds entry would leave NAs here, and nlminb accepts an NA bound
      # silently -- returning convergence = 0 with objective = Inf and par = NA
      # rather than erroring. Also catches split() having recycled a short vector
      # up to the right length, which the check above cannot see.
      if (anyNA(L) || anyNA(U)) {
        stop("Parameter bounds contain NA after alignment; a bounds entry is ",
             "shorter than its parameter.")
      }
      # Keep the aligned vectors, not just the per-block lists. This is the only
      # place the alignment is known to be right -- it happens after obj exists,
      # against obj$par's own block order, with the assertions above -- and it is
      # what a sampler needs. tmbstan::tmbstan(lower =, upper =) takes a vector
      # of exactly this length and expands it across the random effects itself,
      # so `bounds$lower` (a list keyed by parameter block) is the wrong object
      # to hand it: indexing that list by names(obj$par) repeats whole blocks and
      # silently returns a vector orders of magnitude too long.
      bounds$par_lower <- L
      bounds$par_upper <- U
    }

    mod_objects <- c(
      list(
        TMBfilename = TMBfilename,
        bounds      = bounds,
        map         = map
      ),
      mod_objects)

    if (verbose > 0) { message(paste0("Step ", step, ": Hindcast build complete")) }
    step <- step + 1


    # * Optimize hindcast ----
    if (estimateMode %in% c(0, 1, 4)) {
      opt <- suppressMessages(
        .fit_tmb(obj                 = obj,
                 lower               = L,
                 upper               = U,
                 loopnum             = loopnum,
                 newtonsteps         = newtonsteps,
                 getsd               = getsd,
                 control             = control,
                 bias.correct        = bias.correct,
                 getJointPrecision   = getJointPrecision,
                 getReportCovariance = getReportCovariance,
                 quiet               = verbose < 2)
      )

      if (verbose > 0 & estimateMode != 4) {
        message("Step ", step, ": Hindcast optimization complete")
        step <- step + 1
      }
      if (verbose > 0 & estimateMode == 4) {
        message("Step ", step, ": 'dummy' optimization complete")
        step <- step + 1
      }

      # Convergence warnings
      if (estimateMode %in% c(0, 1)) {
        if (is.null(opt$SD) & getsd) {

          if (verbose > 0) {
            message("#################################################")
            message("Model did not converge, check 'identified'")
            message("#################################################")
          }

          identified <- tryCatch(
            {
              # check_estimability() print()s its table to stdout, which
              # suppressMessages() cannot catch. Capture it, re-emit under
              # verbose only; the verdict is kept on `identified` regardless.
              res <- NULL
              tbl <- utils::capture.output(
                res <- suppressMessages(.check_estimability(obj)))
              if (verbose > 0 && length(tbl)) {
                message(paste(tbl, collapse = "\n"))
              }
              res
            },
            error = function(e) {
              "Some gradients are high, please improve optimization and only then use `Check_Identifiable`"
            }
          )

          # Make into list if gradients were low for diagnostics
          if (!is.character(identified)) {
            identified_param_list <- obj$env$parList(identified$BadParams$Param_check)
            identified_param_list <- rapply(identified_param_list, function(x) ifelse(x == 0, "Not estimated", x), how = "replace")
            identified_param_list <- rapply(identified_param_list, function(x) ifelse(x == 1, "OK",            x), how = "replace")
            identified_param_list <- rapply(identified_param_list, function(x) ifelse(x == 2, "BAD",           x), how = "replace")
            identified$param_list <- identified_param_list
          }
          mod_objects$identified <- identified
        }
      }

      # Record which parameter blocks the HINDCAST declared random. Reading this
      # back off `$obj` later does not work: under estimateMode = 0 with a
      # harvest control rule, `$obj` is the PROJECTION object, whose map comes
      # from build_hcr_map(all_params_on = FALSE) and sets every hindcast entry
      # to NA -- so TMB drops the whole `random` declaration and
      # obj$env$random is empty. Anything asking "was this block a random
      # effect?" would then get FALSE for every block, on exactly the fits that
      # carry an HCR. retrospective() asks precisely that.
      #
      # Stored as `character(0)`, never NULL, when the fit declared no random
      # effects at all. `random_vars` starts life as `c()`, and assigning NULL
      # into a list DELETES the element -- so a perfectly ordinary
      # `random_rec = FALSE` fit came back with no `random_vars` at all, which is
      # indistinguishable from a fit made before 5.10.0 recorded it. The two
      # need opposite handling: with no random effects, pinning every deviation
      # in a peel is correct, while an unrecorded one has to be reported. An
      # empty character vector says "recorded, and there were none".
      if (estimateMode %in% c(0, 1)) {
        mod_objects$random_vars <-
          if (is.null(random_vars)) character(0) else random_vars
      }

      # Capture the hindcast optimizer convergence snapshot now, before any
      # projection re-optimization overwrites `opt` (see R/0-convergence.R).
      if (estimateMode %in% c(0, 1)) {
        mod_objects$.conv_hindcast <- .capture_opt_convergence(
          opt, obj, bounds = bounds, mapFactor = map$mapFactor,
          random_vars = random_vars, getsd = getsd)
      }
    }


    # * Get MLEs ----
    if (estimateMode > 1) { # debugging / projection-only: use initial parameters
      last_par <- start_par
    } else {
      last_par <- obj$env$parList(par=obj$env$last.par.best)
    }

    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # 10: Run projection ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    if (estimateMode %in% c(0, 2, 4)) {
      if (data_list$HCR != "NoFishing") { # all HCRs except no F and fixed F

        # * 10.1: Single species mode ----
        if (msmMode == 0) {

          # Turn BRP estimation on within likelihood
          data_list_reorganized$forecast <- rep(1, data_list_reorganized$nspp)

          hcr_map <- Rceattle::build_hcr_map(data_list, map, debug = estimateMode > 3)
          if (sum(!is.na(unlist(hcr_map$mapFactor))) == 0) { stop("HCR map of length 0: all NAs") }

          obj <- TMB::MakeADFun(
            data_list_reorganized,
            parameters = last_par,
            DLL        = TMBfilename,
            map        = hcr_map$mapFactor,
            random     = random_vars,
            silent     = verbose != 2
          )

          if (verbose > 0) { message(paste0("Step ", step, ": Projection build complete")) }
          step <- step + 1

          if (data_list$HCR != "ConstantF") { # fixed F does not need estimation
            opt <- suppressMessages(
              .fit_tmb(obj               = obj,
                       loopnum           = loopnum,
                       getsd             = getsd,
                       control           = control,
                       bias.correct      = bias.correct,
                       getJointPrecision = FALSE,
                       quiet             = verbose < 2)
            )
          }
        }


        # * 10.2: Multi-species mode ----
        if (msmMode > 0) {

          for (HCRiter in 1:max(data_list$HCRorder)) {

            hcr_map <- Rceattle::build_hcr_map(data_list, map,
                                               debug         = estimateMode > 3,
                                               all_params_on = FALSE,
                                               HCRiter       = HCRiter)
            if (sum(as.numeric(unlist(hcr_map$mapFactor)), na.rm = TRUE) == 0) { stop("HCR map of length 0: all NAs") }

            # Get SB0: SSB when model is projected forward under no fishing
            data_list_reorganized$forecast <- data_list$HCRorder <= HCRiter # species to include in likelihood
            params_on  <- c(1:data_list$nspp)[which(data_list$HCRorder == HCRiter)] # only update the species being estimated
            quantities <- obj$report(obj$env$last.par.best)
            SB0        <- quantities$ssb[, ncol(quantities$ssb)]
            B0         <- quantities$biomass[, ncol(quantities$biomass)]
            # Mark the species whose unfished reference is now a real
            # no-fishing projection rather than the placeholder. Recorded
            # explicitly because the alternative -- recognising the placeholder
            # by its value -- is a float comparison against 999 that would also
            # null a legitimately-derived 999 mt, and cannot see a workbook that
            # supplied its own MSSB0.
            data_list_reorganized$MSSB0_derived[params_on] <- TRUE
            data_list_reorganized$MSSB0[params_on] <- SB0[params_on]
            data_list_reorganized$MSB0[params_on]  <- B0[params_on]

            # Adjust Ftarget inits
            params_off <- c(1:data_list$nspp)[which(data_list$HCRorder > HCRiter)]
            last_par$log_Ftarget[params_on]  <- 0
            last_par$log_Ftarget[params_off] <- -999

            obj <- TMB::MakeADFun(
              data_list_reorganized,
              parameters = last_par,
              DLL        = TMBfilename,
              map        = hcr_map$mapFactor,
              random     = random_vars,
              silent     = verbose != 2
            )

            if (verbose > 0) { message(paste0("Step ", step, " - HCRiter ", HCRiter, ": Projection build complete. Optimizing.")) }
            step <- step + 1

            if (data_list$HCR != "ConstantF") { # fixed F does not need estimation
              opt <- suppressMessages(
                .fit_tmb(obj               = obj,
                         loopnum           = loopnum,
                         getsd             = getsd,
                         bias.correct      = bias.correct,
                         control           = control,
                         getJointPrecision = FALSE,
                         quiet             = verbose < 2)
              )

              # Update F from opt
              last_par$log_Ftarget[params_on] <- opt$par[1:length(params_on)]
            }
          }
        }


        if (verbose > 0) { message(paste0("Step ", step, ": Projection optimization complete")) }
        step <- step + 1


        # * Update MLEs ----
        if (estimateMode > 2) { # debugging: use initial parameters
          last_par <- start_par
        } else {
          last_par <- obj$env$parList(par=obj$env$last.par.best)
        }


        # * Projection uncertainty ----
        # Updates the model with all hindcast and BRP parameters "turned on"
        # to get uncertainty estimates in the projection.
        if (projection_uncertainty) {

          hcr_map_proj <- Rceattle::build_hcr_map(data_list, map, debug = estimateMode > 3, all_params_on = TRUE)
          if (sum(as.numeric(unlist(hcr_map_proj$mapFactor)), na.rm = TRUE) == 0) { stop("HCR projection map of length 0: all NAs") }

          obj <- TMB::MakeADFun(
            data_list_reorganized,
            parameters = last_par,
            DLL        = TMBfilename,
            map        = hcr_map_proj$mapFactor,
            random     = random_vars,
            silent     = verbose != 2
          )

          opt$SD <- TMB::sdreport(obj)
        }

      } # End estimable BRP/HCR projections
    } # End projection


    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    # 11: Save output ----
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    mod_objects$estimated_params <- last_par
    mod_objects$obj              <- obj

    quantities <- obj$report(obj$env$last.par.best)

    # Warning for discontinuous likelihood
    if (estimateMode %in% c(0:2)) {
      if (!(estimateMode == 2 & data_list$HCR == "ConstantF")) { # no optimization of projections with fixed F
        if (!is.null(opt) && !is.null(opt$SD) && random_rec == FALSE) {
          if (abs(opt$objective - quantities$jnll) > rel_tol) {
            message("#################################################")
            message("Convergence warning (8): discontinuous likelihood")
            message("#################################################")
          }
        }
      }
    }

    mod_objects$quantities <- Rceattle::rename_output(data_list = data_list, quantities = quantities)

    mod_objects$data_list <- calc_mcall_ianelli(data_list = data_list, data_list_reorganized = data_list_reorganized, quantities = quantities)
    mod_objects$data_list <- calc_mcall_ianelli_diet(data_list = mod_objects$data_list, quantities = quantities)

    # Unfished SSB and biomass under predation (metric tons). Section 10.2
    # derives these by projecting under no fishing and writes them into the
    # reorganized copy it refits from, so the fitted object has to carry them
    # forward: anything that refits from `data_list` alone -- `.refit_like()`,
    # `remove_F()`, every `run_mse()` projection -- would otherwise re-enter the
    # template at the 999 mt placeholder and read SSB/999 as depletion.
    mod_objects$data_list$MSSB0_derived <- data_list_reorganized$MSSB0_derived
    mod_objects$data_list$MSSB0 <- data_list_reorganized$MSSB0
    mod_objects$data_list$MSB0  <- data_list_reorganized$MSB0

    # OSA residual metadata for the aggregate (index/catch) series, mapping
    # obsvec positions to fleet/species/year/age. osa_residuals() rebuilds the
    # full composition / caal / diet metadata on demand, so only the aggregate
    # map is stored here.
    mod_objects$obs_ctl <- .obs_ctl

    mod_objects$run_time <- (Sys.time() - start_time)

    # Record the resolved run configuration so save_config(fit) / run_config(fit)
    # reproduce this fit without re-deriving it. Built from the arguments as
    # actually used (after any config= / model_config overlay).
    # Recording the configuration must never cost the caller a converged fit:
    # this runs after the optimization and sdreport, so a failure here is
    # reported and the run_config dropped, not raised.
    mod_objects$run_config <- tryCatch(
      .rce_run_config(
        mc = model_config(msmMode = msmMode, initMode = initMode, avgnMode = avgnMode,
                          suitMode = suitMode, niter = niter, HCR = HCR, recFun = recFun,
                          M1Fun = M1Fun, growthFun = growthFun, qFun = qFun,
                          # The SPECIFICATION, never the built objects. A run
                          # config is written back out as a build_DSEM() call
                          # (.RCE_CONFIG_BUILDERS), and built objects share only
                          # `sem` with that function's arguments -- so storing
                          # them would drop `family`, `estimate_projection` and
                          # `sigmaR_prior_sd`, and load_config() would rebuild a
                          # different model without saying so.
                          selFun = selFun, compFun = compFun,
                          dsem = data_list$dsem_settings),
        estimateMode = estimateMode, random_rec = random_rec, random_q = random_q,
        random_sel = random_sel, suit_styr = suit_styr, suit_endyr = suit_endyr,
        fc = fit_control),
      error = function(e) {
        warning("Could not record the run configuration on this fit: ",
                conditionMessage(e), call. = FALSE)
        NULL
      })

    if (estimateMode < 3) {
      if (!(estimateMode == 2 & data_list$HCR == "ConstantF")) { # no optimization of projections with fixed F
        mod_objects$opt   <- opt
        mod_objects$sdrep <- opt$SD
      }
    }

    class(mod_objects) <- "Rceattle"

    # Convergence diagnostics (R/0-convergence.R): attach the structured object
    # and surface non-OK checks via message(). message() + tryCatch so a
    # non-converged fit is never turned into an error and is always returned.
    if (estimateMode %in% c(0, 1)) {
      mod_objects$convergence <- tryCatch(
        convergence_diagnostics(mod_objects),
        error = function(e) NULL
      )
      tryCatch({
        if (!is.null(mod_objects$convergence) &&
            mod_objects$convergence$status %in% c("WARN", "FAIL")) {
          message("Convergence diagnostics flagged (status: ",
                  mod_objects$convergence$status,
                  "); inspect fit$convergence.")
          for (ch in mod_objects$convergence$checks) {
            if (ch$severity %in% c("WARN", "FAIL")) {
              message(sprintf("  [%s] %s: %s", ch$severity, ch$id, ch$message))
            }
          }
        }
      }, error = function(e) NULL)
    }

    if (!is.null(file)) {
      save(mod_objects, file = paste0(file, ".Rdata"))
    }

    # Free up memory
    if (estimateMode %in% 0:1) {
      TMB::FreeADFun(obj)
    }

    return(mod_objects)
  }
