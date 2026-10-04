#' Print method for fitted Rceattle models
#'
#' Provides a compact summary so that auto-printing inside R Markdown /
#' knitr / RStudio does not recurse into the (very deep) data and TMB
#' objects stored on the fit. Only structural metadata, convergence
#' status, headline derived quantities, and the package / TMB-DLL
#' versions used to produce the fit are printed.
#'
#' For operational use, the package version line is meant to make it
#' obvious which version of `Rceattle` produced an archived fit so that
#' results can be reproduced even if `main` has moved on. Tag a
#' release (`devtools::install_github("afsc-assessments/Rceattle@X.Y.Z")`,
#' bare version, no leading `v`)
#' and the same version string will reappear here on a fresh run.
#'
#' @param x An object of class \code{"Rceattle"} returned by [fit_mod()].
#' @param ... Currently unused.
#'
#' @return The input `x`, invisibly.
#' @export
print.Rceattle <- function(x, ...) {
  dat <- x$data_list
  pkg_ver <- tryCatch(
    as.character(utils::packageVersion("Rceattle")),
    error = function(e) NA_character_
  )

  cat("<Rceattle model>\n")
  cat("  Rceattle  :", pkg_ver, "\n")

  # Model specification as an indented spec tree (dimensions -> fleets ->
  # processes -> linkages -> config), shared with print.Rceattle_data(). Guarded
  # so a rendering edge case can never break auto-printing a fitted model.
  tree <- tryCatch(.rce_spec_tree(dat),
                   error = function(e) paste0("  (spec tree unavailable: ",
                                              conditionMessage(e), ")"))
  cat(paste(tree, collapse = "\n"), "\n", sep = "")

  # Assembled as rows and emitted together, so the tree closes with the right
  # glyph and every value sits in one column. Written line-by-line it drifted:
  # initMode carried a mid-tree "|-" with nothing after it, the fit statistics
  # then left the tree entirely with a different indent, and cat()'s separator
  # left a trailing space on every line.
  #
  # Switch codes are aliased to their string forms so the block reads the same
  # whatever the underlying code, matching the spec tree above.
  rows <- list()
  add_row <- function(k, v) rows[[length(rows) + 1L]] <<- c(k, v)

  add_row("msmMode", if (is.null(dat$msmMode)) NA else .rce_alias_show(dat$msmMode, msmMode_map))
  # HCR = 0 is the code for "NoFishing", which as a bare label reads as a
  # management choice. Under an estimateMode that runs no projection there is no
  # harvest control rule in play at all, so say that instead.
  # unname() throughout: .rce_alias_show() carries the map's name through, so a
  # bare identical() against the string is FALSE and the guard below silently
  # never fires. cat() hides that, which is how it went unnoticed.
  hcr_txt <- if (is.null(dat$HCR)) "(none)" else
    unname(.rce_alias_show(dat$HCR, hcr_map))
  est_mode <- unname(tryCatch(.rce_alias_show(dat$estimateMode, estimateMode_map),
                              error = function(e) NA_character_))
  if (identical(est_mode, "Hindcast") && identical(hcr_txt, "NoFishing")) {
    hcr_txt <- "(not applicable -- hindcast only)"
  }
  add_row("HCR", hcr_txt)
  add_row("initMode", if (is.null(dat$initMode)) NA else .rce_alias_show(dat$initMode, initMode_map))

  if (!is.null(x$opt) && !is.null(x$opt$objective)) {
    add_row("-log L", signif(x$opt$objective, 6))
  }
  if (!is.null(x$opt$max_gradient)) {
    add_row("max |grad|", signif(x$opt$max_gradient, 4))
  } else if (!is.null(x$obj)) {
    g <- tryCatch(max(abs(as.numeric(x$obj$gr()))), error = function(e) NA_real_)
    if (is.finite(g)) add_row("max |grad|", signif(g, 4))
  }
  if (!is.null(x$run_time)) {
    # signif to match the objective above, rather than format()'s seven digits
    # on a wall-clock time nobody reads to the microsecond.
    secs <- tryCatch(as.numeric(x$run_time, units = "secs"), error = function(e) NA_real_)
    add_row("run time", if (is.finite(secs)) paste(signif(secs, 3), "secs")
                        else format(x$run_time))
  }
  if (!is.null(x$convergence)) {
    add_row("converged", paste0(x$convergence$status,
                                if (x$convergence$status != "OK")
                                  "  (see fit$convergence)" else ""))
  }

  if (length(rows)) {
    cat("  fit\n")
    kw <- max(vapply(rows, function(r) nchar(r[1]), 1L))
    for (i in seq_along(rows)) {
      glyph <- if (i == length(rows)) "\u2514\u2500" else "\u251c\u2500"
      cat(paste0("  ", glyph, " ", formatC(rows[[i]][1], width = kw, flag = "-"),
                 " : ", rows[[i]][2], "\n"))
    }
  }
  invisible(x)
}


#' Compact summary method for Rceattle fits
#'
#' @param object An object of class \code{"Rceattle"} returned by [fit_mod()].
#' @param ... Currently unused.
#'
#' @export
summary.Rceattle <- function(object, ...) {
  # For a fitted assessment the conventional answer to summary() is the
  # estimates and their uncertainty plus the likelihood decomposition, joined
  # from coef(), vcov() and quantities$jnll_comp. print() gives the spec tree;
  # this gives the numbers.
  est <- tryCatch(stats::coef(object), error = function(e) NULL)
  se  <- tryCatch({
    v <- stats::vcov(object)
    # An indefinite Hessian inverts to a negative variance, which has no
    # standard error; NA reports it as unavailable, not as sqrt(|variance|).
    if (is.null(v)) NULL else .conv_se_from_cov(v)
  }, error = function(e) NULL)

  coefs <- NULL
  if (!is.null(est) && length(est)) {
    # vcov() covers the fixed effects only, and getsd = FALSE leaves sdrep NULL,
    # so the standard error is NA rather than absent -- the column has to line up
    # with the estimates either way.
    s <- rep(NA_real_, length(est))
    if (!is.null(se) && length(se) == length(est)) s <- unname(se)
    coefs <- data.frame(parameter = names(est), estimate = unname(est),
                        std_error = s, stringsAsFactors = FALSE)
  }

  jnll <- tryCatch(object$quantities$jnll_comp, error = function(e) NULL)

  structure(
    list(spec = object, coefficients = coefs, jnll_comp = jnll,
         objective = tryCatch(object$opt$objective, error = function(e) NULL),
         convergence = object$convergence),
    class = "summary.Rceattle")
}


#' Print method for an Rceattle model summary
#'
#' @param x A `"summary.Rceattle"` object from [summary.Rceattle()].
#' @param n Number of parameters to show, in `coef()` order. Default 10; use
#'   `Inf` for all, or take `x$coefficients`.
#' @param ... Currently unused.
#' @return `x`, invisibly.
#' @export
print.summary.Rceattle <- function(x, n = 10, ...) {
  print(x$spec)

  if (!is.null(x$coefficients) && nrow(x$coefficients)) {
    cf <- x$coefficients
    cat("\n  estimated parameters (", nrow(cf), ")\n", sep = "")
    show <- utils::head(cf, n)
    show$estimate  <- formatC(show$estimate, format = "g", digits = 4)
    show$std_error <- ifelse(is.na(show$std_error), "-",
                             formatC(show$std_error, format = "g", digits = 3))
    .rce_diag_table(show, c("parameter" = "parameter", "estimate" = "estimate",
                            "std error" = "std_error"))
    if (nrow(cf) > n) {
      cat("  ... ", nrow(cf) - n, " more; all of them are in $coefficients\n",
          sep = "")
    }
    if (all(is.na(cf$std_error))) {
      # Two reasons, and the user needs to know which: no sdreport at all, or one
      # whose variances are unusable because the Hessian is not positive definite.
      # A print method must not stop, so an unreadable spec falls back to the
      # commoner reason.
      v <- tryCatch(stats::vcov(x$spec), error = function(e) NULL)
      if (is.null(v)) {
        cat("  standard errors are unavailable -- this fit has no sdreport",
            "(fit_control(getsd = FALSE))\n")
      } else {
        cat("  standard errors are unavailable -- no variance in the sdreport is",
            "usable; see convergence_diagnostics() 'pdHess'\n")
      }
    }
  }

  if (!is.null(x$jnll_comp)) {
    tot <- sum(x$jnll_comp, na.rm = TRUE)
    cat("\n  likelihood components summing to ",
        formatC(tot, format = "f", digits = 4), "\n", sep = "")
    rn <- rownames(x$jnll_comp)
    by_row <- rowSums(x$jnll_comp, na.rm = TRUE)
    keep <- which(abs(by_row) > 0)
    if (length(keep)) {
      d <- data.frame(component = if (is.null(rn)) keep else rn[keep],
                      value = formatC(by_row[keep], format = "f", digits = 3),
                      stringsAsFactors = FALSE)
      .rce_diag_table(d, c("component" = "component", "-log L" = "value"))
    }
    cat("  the full fleet-by-component matrix is in $jnll_comp\n")
  }
  invisible(x)
}


#' Plot method for fitted Rceattle models
#'
#' Thin S3 dispatcher around the package's existing `plot_*()` functions
#' so that `plot(fit)` works the way users expect. Pick the panel with
#' `what`; everything in `...` is forwarded to the underlying function.
#'
#' @param x An object of class \code{"Rceattle"} returned by [fit_mod()].
#' @param what Character. One of `"biomass"` (default), `"ssb"`,
#'   `"recruitment"`, `"depletion"` (total biomass / B0), `"ssb_depletion"`
#'   (female spawning biomass / SB0, the quantity a Tier 3 HCR compares
#'   against B40%), `"index"`, `"catch"`, `"selectivity"`, `"mortality"`, or
#'   `"data"`.
#' @param ... Passed to the underlying plotting function.
#'
#' @return Invisibly returns `NULL`. Called for the side effect of
#'   producing a plot.
#' @export
plot.Rceattle <- function(x, what = "biomass", ...) {
  what <- match.arg(
    what,
    c("biomass", "ssb", "recruitment", "depletion", "ssb_depletion",
      "index", "catch", "selectivity", "mortality", "data")
  )
  switch(
    what,
    biomass      = plot_biomass(x, ...),
    ssb          = plot_ssb(x, ...),
    recruitment  = plot_recruitment(x, ...),
    depletion    = plot_depletion(x, ...),
    # SSB depletion against SB0 is the management quantity under the
    # Amendment-56 proxies, and this dispatcher could not reach it: "depletion"
    # is plot_depletion(), which draws TOTAL biomass depletion. Note the current
    # plotter is plot_depletionSSB(); plot_ssb_depletion() is its deprecated
    # alias.
    ssb_depletion = plot_depletionSSB(x, ...),
    index        = plot_index(x, ...),
    catch        = plot_catch(x, ...),
    selectivity  = plot_selectivity(x, ...),
    mortality    = plot_mortality(x, ...),
    data         = plot_data(x, ...)
  )
  invisible(NULL)
}


#' Extract estimated parameters from an Rceattle fit
#'
#' @param object An object of class \code{"Rceattle"} returned by [fit_mod()].
#' @param ... Currently unused.
#'
#' @return A named numeric vector of estimated fixed-effect parameters
#'   (the optimizer's `par`). Returns `NULL` if the model was not
#'   estimated.
#' @export
coef.Rceattle <- function(object, ...) {
  if (is.null(object$opt) || is.null(object$opt$par)) return(NULL)
  object$opt$par
}


#' Variance-covariance matrix for an Rceattle fit
#'
#' Returns the fixed-effect covariance matrix produced by
#' [TMB::sdreport()]. Random-effect covariance is not returned here,
#' use `object$sdrep` for the full report.
#'
#' @param object An object of class \code{"Rceattle"} returned by [fit_mod()].
#' @param ... Currently unused.
#'
#' @return A numeric matrix, or `NULL` if `sdreport` was not run
#'   (i.e. the fit was produced with `getsd = FALSE`).
#' @export
vcov.Rceattle <- function(object, ...) {
  sdr <- object$sdrep
  if (is.null(sdr) || is.null(sdr$cov.fixed)) return(NULL)
  sdr$cov.fixed
}


#' Log-likelihood of an Rceattle fit
#'
#' Returns the joint negative-objective from `nlminb`, sign-flipped, as
#' a `"logLik"` object. `df` is the number of estimated fixed-effect
#' parameters. The `nobs` attribute is intentionally omitted: counting
#' "observations" in a stock assessment likelihood (with composition
#' cells, indices, catches, and priors) is not well-defined, so
#' [stats::AIC()] works (uses `df`) while [stats::BIC()] does not.
#'
#' @param object An object of class \code{"Rceattle"} returned by [fit_mod()].
#' @param ... Currently unused.
#'
#' @return An object of class `"logLik"`, or `NULL` if the model was
#'   not estimated.
#' @export
logLik.Rceattle <- function(object, ...) {
  if (is.null(object$opt) || is.null(object$opt$objective)) return(NULL)
  ll <- -as.numeric(object$opt$objective)
  attr(ll, "df") <- length(object$opt$par)
  class(ll) <- "logLik"
  ll
}


#' Residuals from an Rceattle fit
#'
#' Returns a long-format data frame of residuals, following the convention of
#' [stats::residuals.glm()] where `type` selects the *kind* of residual. By
#' default residuals are returned for every applicable data source, survey
#' indices, fishery catches, age/length composition (`comp`), and conditional
#' age-at-length (`caal`); use `source` to restrict to particular ones.
#'
#' Residual kinds (`type`):
#' \describe{
#'   \item{`"response"`}{Observed minus fitted. For `index` / `catch` this is on
#'     the log scale by default (matching the lognormal likelihood; set
#'     `scale = "natural"` for the arithmetic difference); for `comp` / `caal`
#'     it is the difference in proportions, observed minus fitted.}
#'   \item{`"pearson"`}{Standardized residuals. For `catch` and a log-scale
#'     survey index,
#'     \eqn{(\log o - (\log\hat{o} - b\,\sigma^2/2))/\sigma} using the model's
#'     realized observation log-SD \eqn{\sigma} and the observation
#'     bias-adjustment flag \eqn{b} (`bias_adjust_obs`, default 1). A
#'     natural-scale index fleet (`Index_distribution` `"MVN"`, `"MVNORM"`,
#'     `"Normal"` or `"TruncatedNormal"`) holds an ABSOLUTE \eqn{\sigma} and is
#'     standardized as \eqn{(o - \hat{o})/\sigma} instead. Two caveats there: a
#'     covariance fleet gets its marginal residual, not the whitened one (use
#'     `type = "osa"` for that), and `"TruncatedNormal"` is standardized on the
#'     untruncated moments, so it is approximate where \eqn{\hat{o}/\sigma} is
#'     small enough for truncation to shift the mean. For `comp` / `caal`,
#'     \eqn{(p - \hat{p})/\sqrt{\hat{p}(1 - \hat{p})/N}} with input sample size
#'     N.}
#'   \item{`"osa"`}{One-step-ahead residuals via [osa_residuals()], which builds
#'     the composition observation data on demand from any fit. Pass its
#'     arguments through `...`; note that `method = "cdf"` forms no conditional
#'     mean, so `Fitted` is `NA` for every row under it.}
#'   \item{`"process"`}{Process residuals via [process_residuals()] for the
#'     model's random-effect deviations; `source` does not apply.}
#' }
#'
#' Composition rows are returned in long form (one row per observation x
#' age/length bin) and hold the `Age0_Length1` flag from `comp_data` (`0` age,
#' `1` length); CAAL rows hold both the conditioning `Length` and the age `Bin`.
#'
#' Where a fleet uses tail accumulation (`Comp_accum_young` /
#' `Comp_accum_old`), composition residuals describe the bins the likelihood
#' actually fit: the tails are folded into their boundary bin and the bins
#' outside the window are not reported, because the model never fit them
#' separately. `Bin` names the age or length each residual belongs to and
#' `Accumulated` marks the boundary bins, which stand for a range rather than
#' the single bin they are named for. Fleets without accumulation are unchanged.
#'
#' @param object An object of class \code{"Rceattle"} returned by [fit_mod()].
#' @param type Residual kind: one of `"response"` (default), `"pearson"`,
#'   `"osa"`, or `"process"`.
#' @param source Data source(s) to include: any of `"index"`, `"catch"`,
#'   `"comp"`, `"caal"`, or `"all"` (default). `"diet"` (predator
#'   stomach-content composition) is also accepted but, because it uses a
#'   predator/prey schema rather than the fleet/bin layout, must be requested on
#'   its own. Ignored when `type = "process"`.
#' @param scale `"log"` (default) or `"natural"`. Only affects `"response"`
#'   residuals for `index` / `catch`.
#' @param species Optional species code(s) to include (matched against the
#'   `Species` column). Default `NULL` keeps all species. Mirrors the `species`
#'   argument of [plot.rceattle_osa()].
#' @param ... Passed to [osa_residuals()] (e.g. `method`, `seed`) when
#'   `type = "osa"`, or to [process_residuals()] when `type = "process"`.
#'
#' @return A `data.frame` with columns `Source`, `Fleet_code`, `Fleet_name`,
#'   `Species`, `Sex`, `Year`, `Length`, `Bin`, `Age0_Length1`, `Sample_size`,
#'   `Accumulated`, `Observed`, `Fitted`, `Residual`. Columns are `NA` where they
#'   do not apply; `Accumulated` is `FALSE` except on a composition bin that
#'   absorbed a folded tail.
#' @export
residuals.Rceattle <- function(object, type = "response", source = "all",
                               scale = "log", species = NULL, ...) {
  # Following stats::residuals.glm(), `type` selects the residual *kind* and
  # `source` the data source(s).
  type   <- match.arg(type, c("response", "pearson", "osa", "process"))
  source <- match.arg(source, c("all", "index", "catch", "comp", "caal", "diet"),
                      several.ok = TRUE)
  # TODO(review): osa path -- osa_residuals("all") includes diet, but this "all"
  # expansion (forwarded when type = "osa") omits it, and source = c("all",
  # "diet") silently drops diet too. Add "diet" here or document the asymmetry.
  if ("all" %in% source) source <- c("index", "catch", "comp", "caal")
  scale  <- match.arg(scale, c("log", "natural"))

  # Optional species filter, applied to whichever data frame is returned. Keep
  # rows with NA species (e.g. catchability process residuals, which belong to a
  # fleet, not a species) so a species filter does not silently drop them.
  .sp_filter <- function(df) {
    if (is.null(species)) df
    else df[df$Species %in% species | is.na(df$Species), , drop = FALSE]
  }

  # Diet (predator stomach-content) composition uses a predator/prey schema that
  # does not fit the fleet/bin layout of the other sources, so it is returned on
  # its own. (For type = "osa" it flows through to osa_residuals() below.)
  if ("diet" %in% source && type %in% c("response", "pearson")) {
    if (!all(source == "diet")) {
      stop("source = 'diet' uses a predator/prey schema and must be requested ",
           "on its own (not combined with index/catch/comp/caal).")
    }
    dd <- object$data_list$diet_data
    if (is.null(dd) || nrow(dd) == 0 || is.null(object$quantities$diet_hat)) {
      stop("This model has no fitted diet composition (need msmMode > 0 with ",
           "diet_data).")
    }
    obs <- dd$Stomach_proportion_by_weight
    hat <- object$quantities$diet_hat[, 2]
    res <- if (type == "pearson") .rce_diet_pearson(object, dd, obs, hat) else
      obs - hat
    return(.sp_filter(data.frame(
      Source = "diet", Species = dd$Pred, Pred_sex = dd$Pred_sex,
      Prey = dd$Prey, Prey_sex = dd$Prey_sex, Pred_age = dd$Pred_age,
      Prey_age = dd$Prey_age, Year = dd$Year, Sample_size = dd$Sample_size,
      Observed = obs, Fitted = hat, Residual = res, stringsAsFactors = FALSE)))
  }

  # Process residuals (recruitment / random-effect deviations) are computed
  # separately and reshaped into the common schema; `source` does not apply.
  if (type == "process") {
    pr <- process_residuals(object, ...)
    return(.sp_filter(data.frame(
      Source       = pr$source,
      Fleet_code   = NA_integer_,
      Fleet_name   = NA_character_,
      Species      = pr$species,
      Sex          = pr$sex,
      Year         = pr$year,
      Length       = rep(NA_real_, nrow(pr)),
      Bin          = pr$age_length_bin,
      Age0_Length1 = rep(NA_integer_, nrow(pr)),
      Sample_size  = rep(NA_real_, nrow(pr)),
      Observed     = pr$observed,
      Fitted       = pr$predicted,
      Residual     = pr$residual,
      stringsAsFactors = FALSE)))
  }

  # One-step-ahead residuals are computed separately (oneStepPredict) and
  # reshaped into the common schema. They are standard-normal, so `scale` does
  # not apply. `source` selects the observation sources; extra args (method,
  # seed) flow through `...` to osa_residuals().
  if (type == "osa") {
    osa <- osa_residuals(object, source = source, ...)
    fc  <- object$data_list$fleet_control
    return(.sp_filter(data.frame(
      Source       = osa$source,
      Fleet_code   = osa$fleet,
      Fleet_name   = fc$Fleet_name[match(osa$fleet, fc$Fleet_code)],
      Species      = osa$species,
      Sex          = osa$sex,
      Year         = osa$year,
      Length       = osa$length,
      Bin          = osa$age_length_bin,
      Age0_Length1 = ifelse(osa$index_label == "length", 1L,
                            ifelse(osa$index_label == "age", 0L, NA_integer_)),
      Sample_size  = rep(NA_real_, nrow(osa)),
      Observed     = osa$observed,
      Fitted       = osa$predicted,
      Residual     = osa$residual,
      stringsAsFactors = FALSE)))
  }

  q <- object$quantities
  d <- object$data_list

  empty_row <- function(n) {
    data.frame(
      Source       = character(n),
      Fleet_code   = integer(n),
      Fleet_name   = character(n),
      Species      = integer(n),
      Sex          = rep(NA_integer_, n),
      Year         = integer(n),
      Length       = rep(NA_real_, n),
      Bin          = rep(NA_real_, n),
      Age0_Length1 = rep(NA_integer_, n),
      Sample_size  = rep(NA_real_, n),
      # TRUE where tail accumulation folded neighbouring bins into this one, so
      # `Bin` names a range rather than the single age or length it appears to.
      Accumulated  = rep(FALSE, n),
      Observed     = numeric(n),
      Fitted       = numeric(n),
      Residual     = numeric(n),
      # For type = "pearson" on a composition source, the sd the fleet's own
      # likelihood assumes for the observed proportion -- the denominator the
      # residual was divided by. NA elsewhere.
      Sd           = rep(NA_real_, n),
      stringsAsFactors = FALSE
    )
  }

  # Observation bias-adjustment flag (default 1): the index/catch lognormal
  # likelihood fits to mean log(hat) - bias_adjust_obs * sigma^2/2, so the Pearson
  # residual must subtract the same offset.
  ba_obs <- d$bias_adjust_obs
  if (is.null(ba_obs) && !is.null(object$obj)) ba_obs <- object$obj$env$data$bias_adjust_obs
  if (is.null(ba_obs)) ba_obs <- 1
  ba_obs <- as.numeric(ba_obs)[1]

  out <- list()

  if ("index" %in% source &&
      !is.null(d$index_data) && !is.null(q$index_hat)) {
    idx <- d$index_data
    obs <- idx$Observation
    hat <- as.numeric(q$index_hat)
    if (type == "pearson") {
      sigma <- .observation_sd(q, "index")
      sigma <- if (!is.null(sigma)) as.numeric(sigma) else idx$Log_sd
      # Standardize on the scale the fleet's own likelihood uses. The lognormal
      # form below is wrong for a natural-scale family by orders of magnitude,
      # not by a little: sigma there is ABSOLUTE, so -sigma^2/2 is a number the
      # size of the index squared and every residual comes back at the same
      # large constant (an absolute sd of 150 gives ~+75 for every row,
      # regardless of fit).
      nat <- .index_rows_natural_scale(d)
      if (!length(nat)) nat <- rep(FALSE, length(obs))
      res <- rep(NA_real_, length(obs))
      res[!nat] <- (log(obs[!nat]) -
                      (log(hat[!nat]) - ba_obs * sigma[!nat]^2 / 2)) / sigma[!nat]
      res[nat]  <- (obs[nat] - hat[nat]) / sigma[nat]
    } else {
      res <- if (scale == "log") log(obs) - log(hat) else obs - hat
    }
    df <- empty_row(nrow(idx))
    df$Source     <- "index"
    df$Fleet_code <- idx$Fleet_code
    df$Fleet_name <- idx$Fleet_name
    df$Species    <- idx$Species
    df$Year       <- idx$Year
    df$Observed   <- obs
    df$Fitted     <- hat
    df$Residual   <- res
    # Keep only genuine observations: drop projection / NA rows and non-positive
    # values, which carry no (lognormal) residual.
    # TODO(review): held-out rows (Year <= 0) with a positive observation still
    # pass this filter and appear with a negative Year, though the C++ likelihood
    # excludes them (Year>0 & Year<=endyr & flt_type>0). Gate on that same rule;
    # the "drop projection rows" comment above is not fully accurate for them.
    df <- df[is.finite(df$Observed) & is.finite(df$Fitted) &
             df$Observed > 0 & df$Fitted > 0, , drop = FALSE]
    out$index <- df
  }

  if ("catch" %in% source &&
      !is.null(d$catch_data) && !is.null(q$catch_hat)) {
    ctc <- d$catch_data
    obs <- ctc$Catch
    hat <- as.numeric(q$catch_hat)
    if (type == "pearson") {
      sigma <- .observation_sd(q, "catch")
      sigma <- if (!is.null(sigma)) as.numeric(sigma) else ctc$Log_sd
      res   <- (log(obs) - (log(hat) - ba_obs * sigma^2 / 2)) / sigma
    } else {
      res <- if (scale == "log") log(obs) - log(hat) else obs - hat
    }
    df <- empty_row(nrow(ctc))
    df$Source     <- "catch"
    df$Fleet_code <- ctc$Fleet_code
    df$Fleet_name <- ctc$Fleet_name
    df$Species    <- ctc$Species
    df$Year       <- ctc$Year
    df$Observed   <- obs
    df$Fitted     <- hat
    df$Residual   <- res
    # Keep only genuine observations: drop projection / NA rows and non-positive
    # values, which carry no (lognormal) residual.
    # TODO(review): held-out rows (Year <= 0) with a positive observation still
    # pass this filter and appear with a negative Year, though the C++ likelihood
    # excludes them (Year>0 & Year<=endyr & flt_type>0). Gate on that same rule;
    # the "drop projection rows" comment above is not fully accurate for them.
    df <- df[is.finite(df$Observed) & is.finite(df$Fitted) &
             df$Observed > 0 & df$Fitted > 0, , drop = FALSE]
    out$catch <- df
  }

  if ("comp" %in% source &&
      !is.null(d$comp_data) && !is.null(q$comp_hat) &&
      nrow(d$comp_data) > 0) {
    cd <- d$comp_data
    bin_cols <- grep("^Comp_", colnames(cd), value = TRUE)
    obs_mat <- as.matrix(cd[, bin_cols, drop = FALSE])
    row_tot <- rowSums(obs_mat, na.rm = TRUE)
    row_tot[row_tot == 0] <- NA_real_
    obs_prop <- obs_mat / row_tot
    hat_prop <- as.matrix(q$comp_hat)
    bin_vals <- suppressWarnings(as.numeric(sub("^Comp_", "", bin_cols)))
    n_obs <- nrow(cd); n_bin <- length(bin_cols)
    a0l1 <- if (!is.null(cd$Age0_Length1)) cd$Age0_Length1 else NA_integer_

    # Fold the tails onto the bins the likelihood fit (see .fold_comp_props()).
    # NULL means no fleet accumulates, so the rectangular path below is used and
    # every existing model is unchanged.
    folded <- .fold_comp_props(d, cd, obs_prop, hat_prop, a0l1, n_bin)

    if (is.null(folded)) {
      row_id <- rep(seq_len(n_obs), times = n_bin)
      df <- empty_row(n_obs * n_bin)
      df$Source       <- "comp"
      df$Fleet_code   <- rep(cd$Fleet_code, times = n_bin)
      df$Fleet_name   <- rep(cd$Fleet_name, times = n_bin)
      df$Species      <- rep(cd$Species,    times = n_bin)
      df$Sex          <- rep(if (!is.null(cd$Sex)) cd$Sex else NA_integer_,
                             times = n_bin)
      df$Year         <- rep(cd$Year,       times = n_bin)
      df$Bin          <- rep(bin_vals, each = n_obs)
      df$Age0_Length1 <- rep(a0l1, times = n_bin)
      df$Sample_size  <- rep(cd$Sample_size, times = n_bin)
      df$Observed     <- as.numeric(obs_prop)
      df$Fitted       <- as.numeric(hat_prop)
    } else {
      row_id <- folded$row
      df <- empty_row(length(folded$obs))
      df$Source       <- "comp"
      df$Fleet_code   <- cd$Fleet_code[folded$row]
      df$Fleet_name   <- cd$Fleet_name[folded$row]
      df$Species      <- cd$Species[folded$row]
      df$Sex          <- if (!is.null(cd$Sex)) cd$Sex[folded$row] else NA_integer_
      df$Year         <- cd$Year[folded$row]
      df$Bin          <- folded$bin
      df$Age0_Length1 <- a0l1[folded$row]
      df$Sample_size  <- cd$Sample_size[folded$row]
      df$Accumulated  <- folded$acc
      df$Observed     <- folded$obs
      df$Fitted       <- folded$hat
    }
    if (type == "pearson") {
      fw <- .rce_comp_family(object, df$Fleet_code, "Comp_distribution",
                             "Comp_weights", "comp_weights")
      pr <- .rce_comp_pearson(df$Observed, df$Fitted, df$Sample_size,
                              row_id, fw$family, fw$weight,
                              .rce_comp_offset(d))
      df$Residual <- as.numeric(pr)
      df$Sd       <- attr(pr, "sd")
    } else {
      df$Residual <- df$Observed - df$Fitted
    }
    # Drop zero-padded phantom bins (ragged multispecies comps have
    # Observed == Fitted == 0), which otherwise give 0/0 = NaN Pearson residuals.
    df <- df[!is.na(df$Observed) & !is.na(df$Fitted) &
             !(df$Observed == 0 & df$Fitted == 0), , drop = FALSE]
    # The TMB composition guard is Neff > 0, so a row with Sample_size 0 enters
    # no likelihood and must not get a residual. The OSA path already excludes
    # them; this is the same rule on the Pearson side.
    df <- df[is.finite(df$Sample_size) & df$Sample_size > 0, , drop = FALSE]
    out$comp <- df
  }

  if ("caal" %in% source &&
      !is.null(d$caal_data) && !is.null(q$caal_hat) &&
      nrow(d$caal_data) > 0) {
    cd <- d$caal_data
    bin_cols <- grep("^CAAL_", colnames(cd), value = TRUE)
    obs_mat <- as.matrix(cd[, bin_cols, drop = FALSE])
    row_tot <- rowSums(obs_mat, na.rm = TRUE)
    row_tot[row_tot == 0] <- NA_real_
    obs_prop <- obs_mat / row_tot
    hat_prop <- as.matrix(q$caal_hat)
    if (ncol(hat_prop) >= length(bin_cols)) {
      hat_prop <- hat_prop[, seq_along(bin_cols), drop = FALSE]
    }
    bin_vals <- suppressWarnings(as.numeric(sub("^CAAL_", "", bin_cols)))
    n_obs <- nrow(cd); n_bin <- length(bin_cols)
    df <- empty_row(n_obs * n_bin)
    df$Source      <- "caal"
    df$Fleet_code  <- rep(cd$Fleet_code, times = n_bin)
    df$Fleet_name  <- rep(cd$Fleet_name, times = n_bin)
    df$Species     <- rep(cd$Species,    times = n_bin)
    df$Sex         <- rep(if (!is.null(cd$Sex)) cd$Sex else NA_integer_,
                          times = n_bin)
    df$Year        <- rep(cd$Year,       times = n_bin)
    df$Length      <- rep(if (!is.null(cd$Length)) cd$Length else NA_real_,
                          times = n_bin)
    df$Bin         <- rep(bin_vals, each = n_obs)
    df$Sample_size <- rep(cd$Sample_size, times = n_bin)
    df$Observed    <- as.numeric(obs_prop)
    df$Fitted      <- as.numeric(hat_prop)
    if (type == "pearson") {
      fw <- .rce_comp_family(object, df$Fleet_code, "CAAL_distribution",
                             "CAAL_weights", "caal_weights")
      pr <- .rce_comp_pearson(df$Observed, df$Fitted, df$Sample_size,
                              rep(seq_len(n_obs), times = n_bin),
                              fw$family, fw$weight, .rce_comp_offset(d))
      df$Residual <- as.numeric(pr)
      df$Sd       <- attr(pr, "sd")
    } else {
      df$Residual <- df$Observed - df$Fitted
    }
    # Drop zero-padded phantom bins (ragged multispecies comps have
    # Observed == Fitted == 0), which otherwise give 0/0 = NaN Pearson residuals.
    df <- df[!is.na(df$Observed) & !is.na(df$Fitted) &
             !(df$Observed == 0 & df$Fitted == 0), , drop = FALSE]
    df <- df[is.finite(df$Sample_size) & df$Sample_size > 0, , drop = FALSE]
    out$caal <- df
  }

  if (length(out) == 0) return(empty_row(0))
  .sp_filter(do.call(rbind, out))
}


#' The composition proportion offset the fit used
#'
#' Added to both observed and fitted proportions before every composition
#' density. `switch_check()` fills it at 1e-5 when unset, so a saved `data_list`
#' from before it existed still resolves.
#' @param d The model's `data_list`.
#' @noRd
.rce_comp_offset <- function(d) {
  off <- d$comp_offset
  if (is.null(off) || !is.finite(off[1])) 1e-5 else as.numeric(off)[1]
}


#' Pearson residual for an observed proportion
#'
#' Divides by the variance the fleet's own likelihood assumed. `n` is the
#' EFFECTIVE count total, not the input sample size: under a multinomial the
#' weight column multiplies the log-likelihood, so it acts as an effective
#' sample size (the model's own simulator draws at `N * Comp_weights`). `conc`
#' is the Dirichlet-multinomial concentration total, i.e. the sum of the alphas.
#' @param observed,fitted Observed and fitted proportions, on the offset scale
#'   the likelihood works in (see `.rce_comp_pearson()`).
#' @param n Effective count total for the observation.
#' @param conc Dirichlet-multinomial concentration total, or `NULL` for the
#'   multinomial variance.
#' @noRd
.pearson_proportion <- function(observed, fitted, n, conc = NULL) {
  v <- if (is.null(conc)) {
    fitted * (1 - fitted) / n
  } else {
    # A Dirichlet-multinomial proportion is overdispersed relative to the
    # multinomial by (n + conc)/(1 + conc). conc -> Inf recovers the multinomial;
    # conc -> 0 leaves a single effective sample.
    fitted * (1 - fitted) * (n + conc) / (n * (1 + conc))
  }
  (observed - fitted) / sqrt(v)
}


#' Composition Pearson residuals on the likelihood's own scale
#'
#' Rebuilds the proportions, count total and concentration the likelihood used.
#' `comp_offset` is added to both proportions before the density; the weight is
#' an effective sample size under a multinomial and a log concentration under a
#' Dirichlet-multinomial. Totals are summed over the bins the fleet actually
#' fit, so tail-folded rows are handled with no special case.
#'
#' @param obs,hat Observed and fitted proportions, long, one element per bin.
#' @param n_input Input sample size, recycled to the length of `obs`.
#' @param row Row grouping, so a per-observation total sums its own bins.
#' @param family Integer likelihood code per element: 1 = Dirichlet-multinomial,
#'   0 / -1 = multinomial.
#' @param weight Per-element weight parameter, on the scale the family reads it
#'   (natural for a multinomial, log for a Dirichlet-multinomial).
#' @param offset The model's `comp_offset`.
#' @param S Optional per-element normalizing total. `NULL` (the default) sums
#'   each row's own bins, which is right for comp and CAAL. Diet must pass its
#'   own, because it renormalizes over the prey bins PLUS an "other prey" bin
#'   that never reaches the residual frame.
#' @param n_total Optional count total the density sees, per element. `NULL`
#'   uses `n_input * S`, which is comp and CAAL's construction (they scale the
#'   offset proportions by N, so the counts sum to `N * S`). Diet normalizes
#'   BEFORE scaling by N, so its total is `N` and it passes that.
#' @param alpha_scale Optional factor holding the second appearance of the
#'   total in the Dirichlet-multinomial concentration, per element. `NULL` uses
#'   `S`: comp and CAAL build `alpha = (N*S) * hat_offset * theta` where
#'   `hat_offset` sums to `S`. Diet's alphas are built on a renormalized vector
#'   summing to one, so it passes 1.
#' @noRd
.rce_comp_pearson <- function(obs, hat, n_input, row, family, weight, offset,
                              S = NULL, n_total = NULL, alpha_scale = NULL) {
  # Total over the row's own bins, on the offset scale the density sees.
  if (is.null(S)) S <- stats::ave(hat + offset, row, FUN = sum)
  S     <- as.numeric(S)
  p_obs <- (obs + offset) / S
  p_hat <- (hat + offset) / S
  n_row <- if (is.null(n_total)) n_input * S else as.numeric(n_total)
  a_scl <- if (is.null(alpha_scale)) S else as.numeric(alpha_scale)

  res <- sdv <- rep(NA_real_, length(obs))
  dm  <- !is.na(family) & family == 1L

  # Comp_weights = 0 multiplies a multinomial fleet's log-likelihood by zero, so
  # it is not fit and has no effective sample size: NA, not a residual of 0.
  i <- which(!dm & is.finite(weight) & weight > 0 & n_row > 0)
  if (length(i)) {
    sdv[i] <- sqrt(p_hat[i] * (1 - p_hat[i]) / (weight[i] * n_row[i]))
  }
  if (any(dm)) {
    i <- dm
    cc <- n_row[i] * a_scl[i] * exp(weight[i])
    sdv[i] <- sqrt(p_hat[i] * (1 - p_hat[i]) * (n_row[i] + cc) /
                     (n_row[i] * (1 + cc)))
  }
  res <- (p_obs - p_hat) / sdv

  # Travels with the residual so the aggregated-composition band sums the same
  # variances the residuals were standardized by.
  attr(res, "sd") <- sdv
  res
}


#' Diet (stomach-content) Pearson residuals on the likelihood's own scale
#'
#' Each stomach holds an "other prey" balance that enters the density but not
#' the residual frame, so the normalizing total is rebuilt from it. Proportions
#' are renormalized before scaling by the stomach sample size
#' (`ceattle.cpp:4837`), so the count total is `N_s`. The predicted other-prey
#' bin is floored at 1e-5, which is within 1e-5 of the C++ `posfun()` it stands
#' in for, `posfun()` decays towards zero once the prey proportions sum past
#' one, where a flat floor holds. A difference at the scale of the offset itself,
#' on a bin that is a residual balance rather than an observation.
#'
#' @param object The fitted model.
#' @param dd Its `diet_data`.
#' @param obs,hat Observed and fitted prey proportions, one element per row.
#' @noRd
.rce_diet_pearson <- function(object, dd, obs, hat) {
  d   <- object$data_list
  off <- .rce_comp_offset(d)

  # Each stomach is scored separately, on a contiguous run of stomach_id
  # (ceattle.cpp, section 13.2). Without that column the predator/sex/age/year
  # key is one stomach per group, which is every bundled dataset.
  sid <- dd$stomach_id
  if (is.null(sid) || length(sid) != nrow(dd)) {
    sid <- paste(dd$Pred, dd$Pred_sex, dd$Pred_age, dd$Year)
  }

  n_prey <- as.numeric(stats::ave(rep(1, length(obs)), sid, FUN = sum))
  sum_o  <- as.numeric(stats::ave(obs, sid, FUN = sum))
  sum_h  <- as.numeric(stats::ave(hat, sid, FUN = sum))

  # The "other prey" balance, exactly as the C++ forms it.
  other_o <- 1 - pmin(sum_o, 1)
  other_h <- pmax(1 - sum_h, 1e-5)

  # Each vector is normalized by its own offset-inflated total over n_prey + 1
  # bins, so observed and fitted carry different totals here.
  n_bin <- n_prey + 1
  S_o   <- sum_o + other_o + off * n_bin
  S_h   <- sum_h + other_h + off * n_bin

  # Both switches are per-predator vectors on the data_list itself, indexed by
  # species code. `bioenergetics_control` is a workbook sheet, not a data_list
  # element -- reading one off it resolves to NULL and silently scores every
  # predator as a multinomial.
  fw  <- .rce_resolve_family_weight(
    .rce_family_code(.rce_switch_column(d, "Diet_distribution"),
                     "Diet_distribution")[dd$Pred],
    .rce_weight_value(object$estimated_params$diet_comp_weights,
                      d$Diet_comp_weights)[dd$Pred],
    dd$Pred, "predator species")
  fam <- fw$family
  w   <- fw$weight

  # Observed and fitted each divide by their own total, so this does not go
  # through .rce_comp_pearson() (which applies one total to both).
  p_obs <- (obs + off) / S_o
  p_hat <- (hat + off) / S_h
  n     <- as.numeric(dd$Sample_size)

  res <- rep(NA_real_, length(obs))
  dm  <- !is.na(fam) & fam == 1L
  if (any(!dm)) {
    res[!dm] <- .pearson_proportion(p_obs[!dm], p_hat[!dm], w[!dm] * n[!dm])
  }
  if (any(dm)) {
    # alpha = normalized fitted proportions * N_s * theta, so the concentration
    # total is N_s * theta -- no second factor, unlike comp and CAAL.
    res[dm] <- .pearson_proportion(p_obs[dm], p_hat[dm], n[dm],
                                   conc = n[dm] * exp(w[dm]))
  }
  res
}


#' Per-fleet composition family and fitted weight
#'
#' The family column may hold either the switch name or its integer code; the
#' weight is the FITTED parameter, which a refit moves away from the
#' `fleet_control` starting value.
#' @param object The fitted model.
#' @param fleet_code Fleet codes to look up, one per residual row.
#' @param dist_col,weight_col Names of the `fleet_control` switch and starting
#'   weight columns.
#' @param par Name of the `estimated_params` weight vector.
#' @return A list of `family` and `weight`, both aligned to `fleet_code`.
#' @noRd
.rce_comp_family <- function(object, fleet_code, dist_col, weight_col, par) {
  fc  <- object$data_list$fleet_control
  idx <- match(fleet_code, fc$Fleet_code)

  fam <- .rce_family_code(.rce_switch_column(fc, dist_col), dist_col)[idx]
  w   <- .rce_weight_value(object$estimated_params[[par]],
                           .rce_switch_column(fc, weight_col))[idx]

  .rce_resolve_family_weight(fam, w, fleet_code, "fleet")
}


#' Read a switch column by its canonical name or any deprecated spelling
#'
#' `switch_check()` upgrades an alias in place at build time, but a fit SAVED
#' before a rename still holds the old spelling, and a fresh `residuals()`
#' call on one would otherwise find nothing and fall back to the schema default
#', the same silent wrong-family failure the canonical lookup exists to avoid.
#' Works on a `fleet_control` data frame and on the `data_list` itself, both
#' being `[[`-indexable.
#'
#' @param x The `fleet_control` or `data_list` to read from.
#' @param col Canonical schema column name.
#' @return The column, or `NULL` if neither it nor an alias is present.
#' @noRd
.rce_switch_column <- function(x, col) {
  if (!is.null(x[[col]])) return(x[[col]])
  row <- .rce_column_schema()[[col]]
  for (old in row$aliases) if (!is.null(x[[old]])) return(x[[old]])
  NULL
}


#' Normalize a likelihood-family column to its integer switch code
#'
#' The column may hold the switch name or the code itself. `NULL` (an older or
#' hand-built `data_list`) returns `NA`, which the callers resolve to the
#' multinomial.
#' @param x The column, or `NULL`.
#' @param col The schema column name, for its allowed-value map.
#' @noRd
.rce_family_code <- function(x, col) {
  if (is.null(x)) return(NA_integer_)
  if (is.numeric(x)) return(as.integer(x))
  as.integer(.rce_allowed_map(col)[as.character(x)])
}


#' The fitted composition weight, falling back to its starting value
#'
#' A refit moves the weight away from the `fleet_control` column, so the fitted
#' parameter is the one the likelihood used. A model built but never optimized
#' has no `estimated_params`.
#' @param fitted,start The `estimated_params` vector and the control column.
#' @noRd
.rce_weight_value <- function(fitted, start) {
  if (!is.null(fitted)) return(as.numeric(fitted))
  if (!is.null(start))  return(as.numeric(start))
  NA_real_
}


#' Resolve a family/weight pair, demoting a Dirichlet-multinomial with no theta
#'
#' A concentration that cannot be recovered has no safe substitute, so the row
#' falls back to the multinomial variance, too small for an overdispersed
#' composition, and warned about rather than left silent. Leaving the row on the
#' Dirichlet-multinomial with a substituted weight would be worse: `exp(1)` is a
#' concentration of e, which inflates the variance by a factor that looks like a
#' considered choice and is not one.
#'
#' An unresolved family (a hand-built `data_list` holding no switch column) is
#' the schema default, the multinomial, and needs no warning.
#'
#' @param family,weight Resolved family codes and weights, aligned to `ids`.
#' @param ids,label The units to name in the message.
#' @return A list of `family` and `weight` with no `NA`s.
#' @noRd
.rce_resolve_family_weight <- function(family, weight, ids, label) {
  bad <- !is.na(family) & family == 1L & is.na(weight)
  if (any(bad)) {
    warning("No fitted Dirichlet-multinomial weight found for ", label, "(s) ",
            paste(utils::head(unique(ids[bad]), 10L), collapse = ", "),
            "; Pearson residuals there use the multinomial variance, which is ",
            "too small for an overdispersed composition.", call. = FALSE)
  }
  family[bad]            <- 0L
  family[is.na(family)]  <- 0L
  weight[is.na(weight)]  <- 1
  list(family = family, weight = weight)
}


#' Fold composition proportions onto the bins the likelihood actually fit
#'
#' Tail accumulation (`Comp_accum_young` / `Comp_accum_old`) folds a
#' composition's tails into a boundary bin and fits only `[yng, old]`, per fleet
#' and per sex block. A residual taken on the unfolded row therefore describes
#' bins the model never fit separately.
#'
#' Reuses [build_osa_data()]'s fold, so the Pearson and OSA residuals of a fleet
#' cover the same bins with the same labels.
#'
#' @param d The model's `data_list`.
#' @param cd Its `comp_data`.
#' @param obs_prop,hat_prop Observed / predicted proportions, `[row, bin]`.
#' @param a0l1 Per-row `Age0_Length1` (0 = age bins, 1 = length bins).
#' @param n_bin Number of bin columns, including any joint-sex doubling.
#' @return `NULL` when no fleet accumulates (the caller keeps its vectorized
#'   path); otherwise a list of equal-length vectors `row` (index into `cd`),
#'   `bin` (the ordinal each value stands for), `obs`, `hat`, and `acc` (whether
#'   the bin absorbed a folded tail).
#' @keywords internal
#' @noRd
.fold_comp_props <- function(d, cd, obs_prop, hat_prop, a0l1, n_bin) {
  acc_y <- d$comp_accum_young
  acc_o <- d$comp_accum_old
  if (is.null(acc_y)) acc_y <- d$fleet_control$Comp_accum_young
  if (is.null(acc_o)) acc_o <- d$fleet_control$Comp_accum_old
  if (is.null(acc_y) && is.null(acc_o)) return(NULL)

  nages    <- d$nages
  nlengths <- d$nlengths
  fleets   <- cd$Fleet_code
  sexes    <- if (!is.null(cd$Sex)) cd$Sex else rep(NA_integer_, nrow(cd))

  # Per-row window, clamped exactly as build_osa_data() and the cpp do.
  win <- lapply(seq_len(nrow(cd)), function(r) {
    flt <- fleets[r]
    sp  <- cd$Species[r]
    nb  <- if (isTRUE(a0l1[r] == 1)) nlengths[sp] else nages[sp]
    if (is.na(nb) || nb < 1L) nb <- n_bin
    yng <- if (!is.null(acc_y)) acc_y[flt] else 1L
    old <- if (!is.null(acc_o)) acc_o[flt] else 0L
    if (is.na(yng) || yng < 1L) yng <- 1L
    if (is.na(old) || old < 1L || old > nb) old <- nb
    if (yng > old) yng <- old
    list(nb = nb, yng = yng, old = old,
         blk = if (isTRUE(sexes[r] == 3)) 2L else 1L)
  })
  if (!any(vapply(win, function(w) w$yng > 1L || w$old < w$nb, logical(1)))) {
    return(NULL)   # nothing folds: leave the caller's fast path alone
  }

  parts <- lapply(seq_len(nrow(cd)), function(r) {
    w <- win[[r]]
    n_use <- w$blk * w$nb
    # A row declaring more bins than it has columns is inconsistent, and
    # data_check() is where that is judged. Report its raw bins rather than
    # index past the end.
    if (n_use > n_bin) {
      o <- as.numeric(obs_prop[r, seq_len(n_bin)])
      h <- as.numeric(hat_prop[r, seq_len(n_bin)])
      return(list(row = rep(r, n_bin), bin = seq_len(n_bin),
                  obs = o, hat = h, acc = rep(FALSE, n_bin)))
    }
    o <- as.numeric(obs_prop[r, seq_len(n_use)])
    h <- as.numeric(hat_prop[r, seq_len(n_use)])
    does_fold <- w$yng > 1L || w$old < w$nb
    if (does_fold) {
      o <- .fold_comp_bins(o, w$nb, w$blk, w$yng, w$old)
      h <- .fold_comp_bins(h, w$nb, w$blk, w$yng, w$old)
    }
    keep <- seq.int(w$yng, w$old)
    bin  <- rep((seq_len(w$blk) - 1L) * w$nb, each = length(keep)) +
      rep(keep, times = w$blk)
    # Only the boundary bins absorbed a tail, so only they cover more than the
    # age they are named for.
    acc1 <- rep(FALSE, length(keep))
    if (does_fold) {
      if (w$yng > 1L)     acc1[1] <- TRUE
      if (w$old < w$nb)   acc1[length(acc1)] <- TRUE
    }
    list(row = rep(r, length(o)), bin = bin, obs = o, hat = h,
         acc = rep(acc1, times = w$blk))
  })

  list(row = unlist(lapply(parts, `[[`, "row"), use.names = FALSE),
       bin = unlist(lapply(parts, `[[`, "bin"), use.names = FALSE),
       obs = unlist(lapply(parts, `[[`, "obs"), use.names = FALSE),
       hat = unlist(lapply(parts, `[[`, "hat"), use.names = FALSE),
       acc = unlist(lapply(parts, `[[`, "acc"), use.names = FALSE))
}


# Confidence bounds for a derived series, shared by as.data.frame.Rceattle()
# and plot_timeseries() so the table and the figure cannot disagree.
#
# A strictly positive series takes its interval on the log scale,
# exp(log(x) +/- z * sd_log). These series are built multiplicatively
# (R = R0 * exp(rec_dev); n-at-age is a product of survivals), so log(x) is
# close to linear in the estimated parameters where x itself is exponential --
# which is the approximation sdreport() actually makes. The interval is
# right-skewed and cannot cross zero the way the symmetric one does for a weak
# year class or a depleted stock.
#
# `log_sd` is the model's own sd(log x) where it ADREPORTed one. Where it did
# not -- exploitable biomass and the depletions cannot be logged on the tape,
# and fits predating v5.8.0 have no log_* rows -- it is recovered from the
# delta-method identity sd(log x) = sd(x) / x. Where x is not positive the
# quotient is undefined and the symmetric interval stands; the series is
# degenerate there anyway.
.rce_ci_bounds <- function(value, sd, log_sd = NA_real_, z = stats::qnorm(0.975)) {
  n <- length(value)
  sd     <- rep_len(sd, n)
  log_sd <- rep_len(log_sd, n)

  positive <- !is.na(value) & value > 0
  derive   <- is.na(log_sd) & !is.na(sd) & positive
  log_sd[derive] <- sd[derive] / value[derive]

  lwr <- value - z * sd
  upr <- value + z * sd

  use_log <- !is.na(log_sd) & positive
  if (any(use_log)) {
    mu <- log(value[use_log])
    lwr[use_log] <- exp(mu - z * log_sd[use_log])
    upr[use_log] <- exp(mu + z * log_sd[use_log])
  }
  list(lwr = as.numeric(lwr), upr = as.numeric(upr))
}


# Catalogue of derived quantities that as.data.frame.Rceattle() can extract,
# keyed by REPORT name, with their shape: "sy" = matrix(nspp, nyrs);
# "ssay" = array(nspp, max_sex, max_age, nyrs). Quantities that are also
# ADREPORT'd in the TMB template are flagged so we can pull standard errors
# out of sdrep$value / sdrep$sd by name.
.RCEATTLE_QUANTITIES <- list(
  biomass             = list(shape = "sy",   adreport = TRUE),
  ssb                 = list(shape = "sy",   adreport = TRUE),
  R                   = list(shape = "sy",   adreport = TRUE),
  biomass_depletion   = list(shape = "sy",   adreport = TRUE),
  ssb_depletion       = list(shape = "sy",   adreport = TRUE),
  B0                  = list(shape = "sy",   adreport = FALSE),
  SB0                 = list(shape = "sy",   adreport = FALSE),
  DynamicB0           = list(shape = "sy",   adreport = FALSE),
  DynamicSB0          = list(shape = "sy",   adreport = FALSE),
  DynamicSBF          = list(shape = "sy",   adreport = FALSE),
  exploitable_biomass = list(shape = "sy",   adreport = TRUE),
  F_spp               = list(shape = "sy",   adreport = FALSE),
  N_at_age            = list(shape = "ssay", adreport = FALSE),
  biomass_at_age      = list(shape = "ssay", adreport = FALSE),
  Z_at_age            = list(shape = "ssay", adreport = FALSE),
  M_at_age            = list(shape = "ssay", adreport = FALSE),
  M1_at_age           = list(shape = "ssay", adreport = FALSE),
  M2_at_age           = list(shape = "ssay", adreport = FALSE),
  F_at_age            = list(shape = "ssay", adreport = FALSE),
  consumption_at_age  = list(shape = "ssay", adreport = FALSE),
  B_eaten_as_prey     = list(shape = "ssay", adreport = FALSE),
  NByage0             = list(shape = "ssay", adreport = FALSE),
  NByageF             = list(shape = "ssay", adreport = FALSE)
)


#' Tidy long-format derived quantities from an Rceattle fit
#'
#' Returns the model's derived population quantities in long form so that
#' custom plots and post-processing don't have to walk the deeply nested
#' `quantities` list or inherit the dimnames decisions in
#' [rename_output()]. Two shapes are flattened into one tidy frame:
#' species-by-year quantities (e.g. `biomass`, `ssb`, `R`, `F_spp`) and
#' species-by-sex-by-age-by-year quantities (e.g. `N_at_age`,
#' `biomass_at_age`, `M_at_age`). For the species-level shapes, `sex`
#' and `age` are returned as `NA`. Cells of the 4D arrays that are
#' padded out to `max(nsex)` / `max(nages)` for species with fewer
#' sexes or ages are dropped.
#'
#' Standard errors (`se`) and confidence intervals (`lwr`, `upr`) are
#' populated from the TMB `sdreport` for any quantity that was
#' `ADREPORT`'d (`biomass`, `ssb`, `R`, `exploitable_biomass`,
#' `biomass_depletion`, `ssb_depletion`); other quantities and fits produced
#' with `getsd = FALSE` get `NA` for `se` / `lwr` / `upr`. Set `ci_level` to
#' widen or narrow the band.
#'
#' A strictly positive series takes its interval on the log scale, so it is
#' right-skewed and cannot reach zero. This is the same construction
#' [plot_timeseries()] draws, so the table and the figure agree.
#'
#' @param x An object of class \code{"Rceattle"} returned by [fit_mod()].
#' @param row.names,optional Ignored; present for the [as.data.frame()] generic.
#' @param which Character vector of quantity names to extract, or
#'   `"all"` for every quantity with a known shape. Defaults to a common
#'   population-level summary. See `names(Rceattle:::.RCEATTLE_QUANTITIES)`
#'   for the full list.
#' @param ci_level Confidence level for `lwr` and `upr`. Default `0.95`.
#' @param ... Currently unused.
#'
#' @return A `data.frame` with columns `year`, `species`, `sex`, `age`,
#'   `quantity`, `value`, `se`, `lwr`, `upr`. `species` is the character
#'   species name from `data_list$spnames`. Rows are sorted in the order
#'   `which` was given.
#' @export
as.data.frame.Rceattle <- function(x,
                                   row.names = NULL,
                                   optional = FALSE,
                                   which = c("biomass", "ssb", "R",
                                             "biomass_depletion",
                                             "ssb_depletion", "F_spp"),
                                   ci_level = 0.95,
                                   ...) {
  known <- names(.RCEATTLE_QUANTITIES)
  if (length(which) == 1 && identical(which, "all")) which <- known
  bad <- setdiff(which, known)
  if (length(bad)) {
    stop("Unknown quantity name(s): ", paste(bad, collapse = ", "),
         ".\nKnown: ", paste(known, collapse = ", "))
  }
  if (!is.numeric(ci_level) || length(ci_level) != 1 ||
      ci_level <= 0 || ci_level >= 1) {
    stop("`ci_level` must be a single number in (0, 1).")
  }

  d <- x$data_list
  q <- x$quantities
  yrs   <- d$styr:d$projyr
  nyrs  <- length(yrs)
  nspp  <- d$nspp
  nages <- d$nages
  nsex  <- d$nsex
  minage  <- if (!is.null(d$minage)) d$minage else rep(1L, nspp)
  spnames <- if (!is.null(d$spnames)) d$spnames else as.character(seq_len(nspp))

  z <- stats::qnorm(0.5 + ci_level / 2)

  # sdreport vectors are stored in the same column-major order TMB used to
  # stamp the matrix/array, so reshaping `sd` against the report's dim
  # lines up cell-for-cell. A fit run with getsd=FALSE has no $sdrep.
  sdrep <- x$sdrep
  sd_lookup <- function(name, dims) {
    if (is.null(sdrep) || is.null(sdrep$value)) return(NULL)
    idx <- which(names(sdrep$value) == name)
    if (!length(idx) || length(idx) != prod(dims)) return(NULL)
    array(sdrep$sd[idx], dim = dims)
  }

  empty <- data.frame(
    year     = integer(),
    species  = character(),
    sex      = integer(),
    age      = integer(),
    quantity = character(),
    value    = numeric(),
    se       = numeric(),
    lwr      = numeric(),
    upr      = numeric(),
    stringsAsFactors = FALSE
  )

  out <- vector("list", length(which))
  names(out) <- which

  for (qn in which) {
    spec <- .RCEATTLE_QUANTITIES[[qn]]
    arr  <- q[[qn]]
    if (is.null(arr)) next

    if (spec$shape == "sy") {
      mat <- as.matrix(arr)
      if (length(dim(mat)) != 2 ||
          dim(mat)[1] != nspp || dim(mat)[2] != nyrs) next
      sd_mat  <- if (spec$adreport) sd_lookup(qn, dim(mat)) else NULL
      lsd_mat <- if (spec$adreport) sd_lookup(paste0("log_", qn), dim(mat)) else NULL
      # A species with input numbers-at-age reports its input recruits as R
      # (rename_output()); the sdreport's R is the placeholder curve's, so no SE.
      fixed_sp <- which((x$data_list$estDynamics %||% 0) > 0)
      if (qn == "R" && length(fixed_sp)) {
        if (!is.null(sd_mat))  sd_mat[fixed_sp, ]  <- NA
        if (!is.null(lsd_mat)) lsd_mat[fixed_sp, ] <- NA
      }
      grid <- expand.grid(species_idx = seq_len(nspp),
                          year_idx    = seq_len(nyrs),
                          KEEP.OUT.ATTRS = FALSE,
                          stringsAsFactors = FALSE)
      cell <- cbind(grid$species_idx, grid$year_idx)
      val  <- mat[cell]
      sdv  <- if (!is.null(sd_mat)) sd_mat[cell] else NA_real_
      lsd  <- if (!is.null(lsd_mat)) lsd_mat[cell] else rep(NA_real_, length(val))
      bnds <- .rce_ci_bounds(val, sdv, lsd, z)
      out[[qn]] <- data.frame(
        year     = yrs[grid$year_idx],
        species  = spnames[grid$species_idx],
        sex      = NA_integer_,
        age      = NA_integer_,
        quantity = qn,
        value    = as.numeric(val),
        se       = as.numeric(sdv),
        lwr      = bnds$lwr,
        upr      = bnds$upr,
        stringsAsFactors = FALSE
      )
    } else if (spec$shape == "ssay") {
      a4 <- arr
      if (length(dim(a4)) != 4 || dim(a4)[1] != nspp || dim(a4)[4] != nyrs) next
      max_sex <- dim(a4)[2]
      max_age <- dim(a4)[3]
      sd4 <- if (spec$adreport) sd_lookup(qn, dim(a4)) else NULL
      grid <- expand.grid(species_idx = seq_len(nspp),
                          sex_idx     = seq_len(max_sex),
                          age_idx     = seq_len(max_age),
                          year_idx    = seq_len(nyrs),
                          KEEP.OUT.ATTRS = FALSE,
                          stringsAsFactors = FALSE)
      keep <- grid$sex_idx <= nsex[grid$species_idx] &
              grid$age_idx <= nages[grid$species_idx]
      grid <- grid[keep, , drop = FALSE]
      cell <- cbind(grid$species_idx, grid$sex_idx,
                    grid$age_idx, grid$year_idx)
      val  <- a4[cell]
      sdv  <- if (!is.null(sd4)) sd4[cell] else NA_real_
      out[[qn]] <- data.frame(
        year     = yrs[grid$year_idx],
        species  = spnames[grid$species_idx],
        sex      = as.integer(grid$sex_idx),
        age      = as.integer(grid$age_idx + minage[grid$species_idx] - 1L),
        quantity = qn,
        value    = as.numeric(val),
        se       = as.numeric(sdv),
        lwr      = as.numeric(val - z * sdv),
        upr      = as.numeric(val + z * sdv),
        stringsAsFactors = FALSE
      )
    }
  }

  out <- out[!vapply(out, is.null, logical(1))]
  if (!length(out)) return(empty)
  do.call(rbind, c(unname(out), list(make.row.names = FALSE)))
}


#' Simulate data from a fitted Rceattle model
#'
#' @description The [stats::simulate()] method for CEATTLE fits: draws `nsim`
#' replicate data sets from the fitted observation model, and optionally from
#' the process model as well.
#'
#' @details
#' A wrapper on [sim_mod()], which documents the observation model in full. For
#' expected values rather than draws, call `sim_mod(simulate = FALSE)`.
#'
#' Draws are taken by the TMB model, so `simulate()` needs a live `$obj`: a
#' model loaded from disk has one, a [model_average()] result does not, see
#' [sim_mod()].
#'
#' @param object An object of class \code{"Rceattle"} returned by [fit_mod()].
#' @param nsim Number of replicate data sets to draw. Default 1.
#' @param seed Optional seed. The caller's random state is restored afterwards,
#'   per [stats::simulate()]. Recorded as `attr(, "seed")`.
#' @param process Which process error to redraw alongside the observations.
#'   `FALSE` (default) keeps the fitted deviations, so replicates differ only in
#'   observation error. See [sim_mod()] for the alternatives.
#' @param ... Currently unused.
#'
#' @return A list of `nsim` `data_list` objects, always a list, including at
#'   `nsim = 1`, so callers do not have to special-case the length. When
#'   `process` redrew something, each element holds the deviations that
#'   generated it as `attr(, "process_sim")`: a named list of whichever of
#'   `rec_dev`, `init_dev`, `log_M1_dev` and `beta_linkage_re` were drawn, each
#'   with a `_drawn` logical of the same shape marking the cells the draw
#'   touched. Compare estimates against those rather than against `object`,
#'   see [sim_mod()].
#'
#' @seealso [sim_mod()] for the observation model and the `process` options,
#'   [self_test()] for simulating and refitting in one step.
#'
#' @examples
#' \dontrun{
#' fit  <- fit_mod(data_list = BS2017SS, estimateMode = 1)
#' reps <- simulate(fit, nsim = 10, seed = 1)
#' # ...and with recruitment redrawn as well:
#' reps <- simulate(fit, nsim = 10, seed = 1, process = "recruitment")
#' truth <- attr(reps[[1]], "process_sim")$rec_dev
#' }
#' @export
simulate.Rceattle <- function(object, nsim = 1, seed = NULL,
                              process = FALSE, ...) {
  if (!inherits(object, "Rceattle")) {
    stop("`object` must be an Rceattle model fit.", call. = FALSE)
  }
  nsim <- as.integer(nsim)
  if (length(nsim) != 1L || is.na(nsim) || nsim < 1L) {
    stop("`nsim` must be a single positive integer.", call. = FALSE)
  }

  # The stats::simulate() convention: seeding is local, so a caller mid-stream
  # is not silently displaced by having simulated.
  if (!is.null(seed)) {
    if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      stats::runif(1)
    }
    old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
    set.seed(seed)
  }

  out <- lapply(seq_len(nsim), function(i)
    sim_mod(object, simulate = TRUE, process = process))
  names(out) <- paste0("Sim_", seq_len(nsim))
  if (!is.null(seed)) attr(out, "seed") <- seed
  out
}
