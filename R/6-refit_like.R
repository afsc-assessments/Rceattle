#' Re-fit a model reusing a source `data_list`'s structural configuration
#'
#' @description
#' Internal engine behind the repeated [fit_mod()] re-invocations in
#' [retrospective()], [jitter()], [self_test()], [profile.Rceattle()],
#' [run_mse()], [remove_F()], [sample_rec()], and [reweight_comps()]. Each of
#' those re-fits a model that shares a source model's structure but changes a few
#' things (a data peel, a jittered start, a projection year, a set of
#' composition weights).
#'
#' This helper reconstructs the harvest-control-rule, stock-recruit,
#' natural-mortality and growth specifications and the shared
#' process switches (`msmMode`, `suitMode`, `random_rec`, ...) from `data_list`
#' (the "...like" part), then calls [fit_mod()]. Every field that genuinely
#' varies across the callers is a named argument that defaults to the
#' `data_list` value, so each call site shows only the handful of things it
#' actually changes.
#'
#' Console suppression is otherwise left to the caller: the call sites wrap this
#' differently (`suppressMessages()`, `suppressWarnings()`, both, or neither),
#' and suppression is console-only, it does not affect the fit. The one
#' exception is `data_check()`'s warnings, dropped here via
#' `fit_mod(quiet_data_check = TRUE)`, because every caller re-validates a
#' `data_list` the user has already fitted. Convergence and TMB warnings still
#' reach whatever the call site allows through.
#'
#' @param data_list Source `data_list`. Drives both the data to fit and the
#'   reconstructed HCR / recFun / M1Fun / growthFun (unless a field is overridden
#'   below).
#' @param inits Starting parameters (an `estimated_params`-shaped list).
#' @param estimateMode [fit_mod()] estimate mode for this refit. Required and
#'   passed explicitly by every caller (retro/jitter refit at 0, profile at 1,
#'   the MSE OM rebuild at 3, the input-F EM at 2, ...), so it has no default.
#' @param map Parameter map, or `NULL` to let [fit_mod()] build it.
#'
#'   There is deliberately no `bounds` counterpart. Bounds are always rebuilt
#'   from the source `data_list`, which reproduces them exactly for every bound
#'   the schema or the linkage table drives (`linkage_spec(bounds = )`), the
#'   documented way to set them. A raw `fit_mod(bounds = )` hand-override is
#'   therefore the one thing a refit does not hold, and holding it would not be
#'   safe: [run_mse()] grows `log_F` and re-dimensions the selectivity deviation
#'   blocks at every assessment, so the source fit's bounds no longer line up
#'   with the parameters they would be indexed against.
#' @param HCR A pre-built [build_hcr()] object to use instead of reconstructing
#'   one from `data_list` (used for [run_mse()]'s input-F EM refit, whose HCR is
#'   a freshly computed average F rather than the stored rule).
#' @param phase,getsd,loopnum [fit_control()] knobs. `getsd` defaults to `FALSE`
#'   (diagnostics rarely need `sdreport`); `remove_F()` overrides it to `TRUE`.
#'   The bias-adjustment settings are recovered from `data_list` rather than
#'   taken from a caller (see below); every other [fit_control()] field falls
#'   back to its default, so a source model fitted with non-default optimizer
#'   settings (`newtonsteps`, `nlminb_control`, ...) refits under the defaults.
#' @param proj_mean_rec Recruitment projection mode; defaults to the source's.
#'   [run_mse()]'s OM refit pins it to `TRUE` (project on the mean across sims).
#' @param srr_mse_switchyr,srr_hat_styr,srr_hat_endyr Stock-recruit reference
#'   fields; default to the source's values. Callers clamp `srr_mse_switchyr` /
#'   `srr_hat_endyr` to a peel or assessment year, or pin them to the pristine OM.
#' @param suit_styr,suit_endyr Suitability window; default to the source's,
#'   clamped to a peel year or pinned to the pristine OM by some callers.
#' @param projection_uncertainty Whether the refit reports standard errors with
#'   the hindcast parameters on. Defaults to the source's, held on
#'   `data_list`; under an HCR a refit without it reports zero for the whole
#'   hindcast.
#'
#' @return The fitted `Rceattle` object returned by [fit_mod()].
#' @noRd
.refit_like <- function(data_list,
                        inits,
                        estimateMode,
                        map              = NULL,
                        HCR              = NULL,
                        phase            = FALSE,
                        getsd            = FALSE,
                        loopnum          = data_list$loopnum,
                        proj_mean_rec    = data_list$proj_mean_rec,
                        srr_mse_switchyr = data_list$srr_mse_switchyr,
                        srr_hat_styr     = data_list$srr_hat_styr,
                        srr_hat_endyr    = data_list$srr_hat_endyr,
                        suit_styr        = data_list$suit_styr,
                        suit_endyr       = data_list$suit_endyr,
                        projection_uncertainty =
                          isTRUE(data_list$projection_uncertainty)) {
  dl <- data_list

  # Reconstruct the harvest control rule from the source unless the caller
  # supplied one (the input-F EM refit builds its own from a computed average F).
  if (is.null(HCR)) {
    HCR <- build_hcr(
      HCR        = dl$HCR,
      DynamicHCR = dl$DynamicHCR,
      Ftarget    = dl$Ftarget,
      Flimit     = dl$Flimit,
      Ptarget    = dl$Ptarget,
      Plimit     = dl$Plimit,
      Alpha      = dl$Alpha,
      Pstar      = dl$Pstar,
      Sigma      = dl$Sigma,
      Fmult      = dl$Fmult,
      HCRorder   = dl$HCRorder)
  }

  # A fit stored with srr_fun 1|3|5 refits as the structural form it fitted. The warning
  # escapes the suppressWarnings() below, but retrospective(), jitter(), profile() and self_test() hide it.
  srr_fun_refit      <- .srr_fun_structural(dl$srr_fun)
  srr_pred_fun_refit <- .srr_fun_structural(
    dl$srr_pred_fun, penalty = isTRUE(as.integer(dl$srr_fun)[1] == 0L))

  # These inits hold the fitted (or, in profile(), the grid) value, so fit_mod() must not
  # put a linkage-fixed intercept back to its init.
  if (is.list(inits)) attr(inits, "rceattle_refit") <- TRUE

  fit_mod(
    data_list    = dl,
    inits        = inits,
    map          = map,
    bounds       = NULL,
    file         = NULL,
    estimateMode = estimateMode,
    HCR          = HCR,
    # suppressWarnings: build_srr()'s advisory warnings only.
    recFun = suppressWarnings(build_srr(
      srr_fun          = srr_fun_refit,
      srr_pred_fun     = srr_pred_fun_refit,
      proj_mean_rec    = proj_mean_rec,
      srr_mse_switchyr = srr_mse_switchyr,
      srr_hat_styr     = srr_hat_styr,
      srr_hat_endyr    = srr_hat_endyr,
      srr_est_mode     = dl$srr_est_mode,
      srr_prior        = dl$srr_prior,
      srr_prior_sd     = dl$srr_prior_sd,
      srr_alpha_init   = dl$srr_alpha_init,
      srr_beta_init    = dl$srr_beta_init,
      Bmsy_lim         = dl$Bmsy_lim,
      linkages         = dl$srr_linkages)),
    # suppressWarnings: legacy M1_indices may travel via data_list.
    M1Fun = suppressWarnings(build_M1(
      M1_model     = dl$M1_model,
      M1_re        = dl$M1_re,
      updateM1     = FALSE,
      M1_use_prior = dl$M1_use_prior,
      M2_use_prior = dl$M2_use_prior,
      M_prior      = dl$M_prior,
      M_prior_sd   = dl$M_prior_sd,
      M1_indices   = dl$M1_indices,
      linkages     = dl$M1_linkages)),
    growthFun = build_growth(fun      = dl$growth_fun,
                             linkages = dl$growth_linkages),
    # fit_mod() takes these three specifications as the source of truth for
    # data_list$q_linkages / sel_linkages / comp_linkages, so they are rebuilt
    # from the source model like the HCR / recruitment / M / growth specs above.
    # That holds the linkage table -- and with it the size of the beta_linkage
    # block the warm-start `inits` carry -- fixed across the refit. Each
    # build_*(NULL) is a no-op for a model that uses no such linkage.
    qFun       = build_catchability(linkages = dl$q_linkages),
    selFun     = build_selectivity(linkages = dl$sel_linkages),
    compFun    = build_composition(linkages = dl$comp_linkages),
    random_rec = dl$random_rec,
    random_q   = isTRUE(as.logical(dl$random_q)),
    random_sel = isTRUE(as.logical(dl$random_sel)),
    niter      = dl$niter,
    msmMode    = dl$msmMode,
    avgnMode   = dl$avgnMode,
    suitMode   = dl$suitMode,
    suit_styr  = suit_styr,
    suit_endyr = suit_endyr,
    initMode   = dl$initMode,
    fit_control = fit_control(
      phase   = phase,
      loopnum = loopnum,
      getsd   = getsd,
      verbose = 0,
      # Bias adjustment defines the LIKELIHOOD, not the optimizer, so a model
      # fitted with it off has to refit with it off or the diagnostic compares
      # two different objectives. It cannot be left to the default: fit_mod()
      # resets data_list's copies to 1 and then re-applies fit_control's, whose
      # defaults are TRUE, so a freshly built control silently switches it back
      # on (worth ~880 jnll units on BS2017SS). The resolved values ride on the
      # data_list, so read them back from there -- the same way `loopnum` above
      # recovers its value -- which covers every call site without one of them
      # having to remember to pass it.
      bias_adjust_obs  = .refit_bias_adjust(dl$bias_adjust_obs),
      bias_adjust_proc = .refit_bias_adjust(dl$bias_adjust_proc),
      # Same reason as the bias-adjustment flags: a freshly built control would
      # reset this to FALSE, and under an HCR that costs the hindcast its
      # standard errors. Gated on `getsd` because fit_mod() runs the extra
      # sdreport() either way, and run_mse() refits here with getsd = FALSE.
      projection_uncertainty = isTRUE(projection_uncertainty) && isTRUE(getsd)),
    # Every caller here re-validates a data_list the user has already fitted, so
    # data_check()'s warnings have been seen; repeating them once per peel,
    # jitter, or MSE iteration is noise. Errors still stop the refit, and
    # convergence / TMB warnings are untouched.
    quiet_data_check = TRUE)
}


#' Recover a bias-adjustment multiplier from a fitted `data_list`.
#'
#' Returned unchanged rather than coerced to a logical. The cpp reads these as
#' `DATA_SCALAR` and uses them as a plain multiplier on the correction
#' (`bias_adjust_obs * sigma^2 / 2`), so a fractional value is a partial
#' bias-adjustment ramp, not a malformed flag, `as.logical()` would quantize
#' 0.5 up to 1 and apply a full bias correction where a half one was configured.
#' `DATA_SCALAR` also means a vector would fail in `MakeADFun`, so taking the
#' first element cannot silently pick the wrong one.
#'
#' `NULL` means the `data_list` holds no resolved value, it never went
#' through [fit_mod()], or it came from a fit predating the field being recorded.
#' In the second case a model fitted with the correction off refits with it on,
#' which is not recoverable from the object; [fit_control()]'s own default is the
#' only available answer.
#'
#' @keywords internal
#' @noRd
.refit_bias_adjust <- function(x) {
  if (is.null(x) || !length(x) || anyNA(x)) return(TRUE)
  x[[1]]
}
