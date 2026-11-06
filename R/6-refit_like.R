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
#' @param map Parameter map, or `NULL` to let [fit_mod()] build it, which releases anything
#'   fixed by hand; a caller that must keep the fitted model passes the stored map.
#'
#'   There is deliberately no `bounds` counterpart. Bounds are always rebuilt
#'   from the source `data_list`, which reproduces them exactly for every bound
#'   the schema or the linkage table drives (`linkage_spec(bounds = )`), the
#'   documented way to set them. A raw `fit_mod(bounds = )` hand-override is
#'   therefore the one thing a refit does not hold, and holding it would not be
#'   safe: [run_mse()] grows `log_F` and re-dimensions the selectivity deviation
#'   blocks at every assessment, so the source fit's bounds no longer line up
#'   with the parameters they would be indexed against.
#' @param map_holds Cells from `.map_holds()` to keep fixed on the map `fit_mod()` builds,
#'   for a caller whose `endyr` moves.
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
#' @param dsem DSEM specification to carry into the refit. Defaults to the
#'   source `data_list`'s stored `dsem_settings`, so a refit inherits it.
#' @return The fitted `Rceattle` object returned by [fit_mod()].
#' @noRd
.refit_like <- function(data_list,
                        inits,
                        estimateMode,
                        map              = NULL,
                        map_holds        = NULL,
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
                        dsem             = data_list$dsem_settings,
                        projection_uncertainty =
                          isTRUE(data_list$projection_uncertainty)) {
  dl <- data_list

  # The DSEM travels as the SPECIFICATION, never as built objects. Its latent
  # states x_tj (and eps_tj / y_tj / the map) are dimensioned over the model's
  # year span, and every caller here changes that span: retrospective() peels
  # endyr, run_mse() shortens projyr. Rebuilding from the spec against each
  # refit's own data_list re-derives the right dimensions; forwarding a built
  # object would hand MakeADFun a latent-state matrix of the wrong length.
  #
  # This is also why the DSEM is not in the `if (is.null(HCR))` style block
  # below: build_dsem_objects() is called by fit_mod(), after it has resolved
  # the refit's own styr/endyr/projyr.


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
  if (is.list(inits)) {
    attr(inits, "rceattle_refit") <- TRUE
    # fit_mod() has no argument for held parameters, so they travel on inits like the flag
    # above; it applies them only to a map it builds itself.
    attr(inits, "rceattle_map_holds") <- map_holds
  }

  fit_mod(
    data_list    = dl,
    inits        = inits,
    map          = map,
    bounds       = NULL,
    file         = NULL,
    estimateMode = estimateMode,
    HCR          = HCR,
    dsem = dsem,
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


# True when a model carries a DSEM. A DSEM can live in three places, and a check
# that looks in only one is silent on exactly the objects it exists for:
#   fit$data_list$dsem_settings  -- the specification, set by fit_mod()
#   fit$dsem                     -- the built objects, the canonical slot on a
#                                   fitted object (what summary.Rceattle reads),
#                                   and all that is set when fit_mod() is handed
#                                   pre-built objects
#   fit$data_list$model_config$dsem -- carried by a stored run configuration
.has_dsem <- function(x) {
  # `[[`, not `$`: `x$dsem` PARTIAL-MATCHES `dsem_settings`, so an object
  # carrying only a specification reported built objects and every caller below
  # took the wrong branch (measured). The schema's partial-match trap, live.
  !is.null(x[["dsem"]]) ||
    !is.null(x$data_list$dsem_settings) ||
    !is.null(x$data_list$model_config$dsem)
}

# True when a DSEM is attached AND it drives `process`. `.has_dsem()` answers
# "is there a SEM"; this answers "does that SEM own this process's deviations",
# and the two stop being the same question as soon as a SEM can drive something
# other than recruitment. Read it off the built objects, which is where the
# template's own flag comes from, so R and C++ cannot disagree about which
# processes a fit's SEM drives.
#
# A fit built before the flag existed could only have been a recruitment DSEM,
# so its absence reads as recruitment -- the same default .dsem_tmb_inputs()
# applies, and the two have to match or a saved fit would refit as a different
# model.
#
# Accepts a fit, anything carrying a `data_list`, or the built dsem object
# itself: fit_mod() holds `mod_objects$dsem` while the diagnostics hold a
# fitted model, and both ask this question.
#' @noRd
.dsem_drives <- function(x, process = "recruitment") {
  known <- names(.DSEM_LATENT_PREFIXES)
  if (length(process) != 1L || !process %in% known) {
    stop("`process` must be one of: ", paste(known, collapse = ", "), ".",
         call. = FALSE)
  }
  # Built objects carry $tmb_inputs; a fit carries them under $dsem. `[[` for
  # the same partial-match reason as .has_dsem().
  obj <- if (!is.null(x[["tmb_inputs"]])) x else x[["dsem"]]
  if (is.null(obj)) {
    # A specification with no built objects beside it. That is not only a
    # pre-fit state -- profile() handles a fitted model in exactly this shape
    # and says so. A growth-only sem is legal, so a spec alone no longer
    # implies recruitment is driven; the honest answer is read off the sem text,
    # which is the only thing available here.
    if (!.has_dsem(x)) return(FALSE)
    # The two spec locations .has_dsem() accepts that are not the built object.
    # There is no top-level `dsem_settings` on a model: fit_mod() stores it at
    # data_list$dsem_settings, and the built object carries $tmb_inputs and so
    # takes the branch above.
    spec <- x$data_list$dsem_settings %||% x$data_list$model_config$dsem
    sem_txt <- paste(as.character(spec$sem %||% ""), collapse = "\n")
    if (!nzchar(sem_txt)) {
      # build_DSEM(sem = NULL) synthesises recruitment self-loops only.
      return(identical(process, "recruitment"))
    }
    stem <- if (identical(process, "recruitment")) {
      .DSEM_LATENT_PREFIXES[["recruitment"]]
    } else {
      .DSEM_LATENT_PREFIXES[[process]]
    }
    return(grepl(paste0("\\b", stem, "[0-9]+[fm]?\\b"), sem_txt))
  }
  # ONE NAMED BRANCH PER PROCESS, and a loud default. The accepted set comes from
  # .DSEM_LATENT_PREFIXES, so deriving the non-growth case as "not a growth
  # target" hands back the RECRUITMENT answer for any process added to that
  # table later -- silently wrong. The stop() below is what makes growing the
  # table fail loudly instead of answering about the wrong process.
  if (identical(process, "recruitment")) {
    return(isTRUE(as.integer(obj$tmb_inputs$data$drives_rec %||% 1L) == 1L))
  }
  if (process %in% names(.DSEM_GROWTH_TARGETS)) {
    # The index map carries one entry per named growth latent, so the question is
    # whether any entry claims this parameter's slot. An object built before
    # growth was drivable has no map and answers FALSE.
    par <- unname(LINKAGE_PARAM_CODES$growth[[sub("^growth_", "", process)]])
    return(isTRUE(par %in% (obj$tmb_inputs$data$growth_dev_par %||% integer(0))))
  }
  if (identical(process, "natural_mortality")) {
    # M's map has no parameter slot to test -- the age dimension takes it -- so
    # a non-empty column map is the whole answer. An object built before M was
    # drivable has no map and answers FALSE.
    return(length(obj[["tmb_inputs"]][["data"]][["mdev_col"]] %||% integer(0)) > 0L)
  }
  stop("`.dsem_drives()` has no branch for process '", process,
       "'. It is registered in `.DSEM_LATENT_PREFIXES` but nothing here answers ",
       "for it; add a branch rather than letting it fall through.", call. = FALSE)
}

# TRUE when the SEM drives ANY growth parameter. The gate in sim_mod() asks the
# question that way because the GMRF draw is joint: it does not matter which
# growth parameter is driven, only whether the field has growth states in it.
#' @noRd
.dsem_drives_growth <- function(x) {
  any(vapply(names(.DSEM_GROWTH_TARGETS),
             function(p) isTRUE(.dsem_drives(x, p)), logical(1)))
}


# There is deliberately no shared "does not support a DSEM" refusal. Only three
# places refuse at all -- run_mse(), sample_rec() on a hindcast-only field, and
# profile() on `rec_dev` -- and each says something specific about why, which one
# generic message could not. Everything else runs on a DSEM; see
# vignette("model-diagnostics") for what differs.

