# =============================================================================
# Fit diagnostics: survey-index and fishery-catch fits, and index residuals.
# All plotters build a tidy data frame, render it with the shared ggplot theme
# (theme_classic) + colourblind palette, and return the ggplot object.
# =============================================================================

# Assemble a tidy observed-vs-predicted data frame for a fleet "fit" plot.
# `kind` selects the index (survey) or catch (fishery) data/quantities. Returns
# one row per (model, fleet, year) with the observation, its lognormal 95%
# band, and the model prediction.
.fleet_fit_df <- function(Rceattle, kind = c("index", "catch"),
                          species = NULL, incl_proj = FALSE, maxyr = NULL,
                          model_names = NULL) {
  kind <- match.arg(kind)
  cfg <- switch(kind,
    index = list(data = "index_data", obs = "Observation",
                 hat = "index_hat", sd = "index", type = "Survey"),
    catch = list(data = "catch_data", obs = "Catch",
                 hat = "catch_hat", sd = "catch", type = "Fishery"))

  model_names_use <- .model_labels(Rceattle, model_names)
  fc <- Rceattle[[1]]$data_list$fleet_control
  # Catch belongs to fisheries by definition; an index belongs to whichever fleet
  # carries index_data. Selecting surveys here dropped a fishery's CPUE rows from
  # the figure without comment.
  codes <- if (kind == "index") {
    .fleets_with_index(Rceattle[[1]]$data_list, fitted_only = FALSE)
  } else {
    fc$Fleet_code[fc$Fleet_type == cfg$type]
  }
  if (is.null(species)) species <- seq_len(Rceattle[[1]]$data_list$nspp)

  out <- list()
  for (i in seq_along(Rceattle)) {
    dat <- Rceattle[[i]]$data_list[[cfg$data]]
    dat$.obs <- dat[[cfg$obs]]
    dat$.sd  <- .observation_sd(Rceattle[[i]]$quantities, cfg$sd)
    dat$.hat <- Rceattle[[i]]$quantities[[cfg$hat]]
    # Year < 0 flags a prediction-only row (excluded from the likelihood; see
    # ceattle.cpp "Likelihood (yr > 0) vs prediction (yr < 0)"). Drop
    # them so they are not drawn as fitted observations.
    keep <- dat$Species %in% species & dat$Fleet_code %in% codes & dat$Year > 0
    if (!incl_proj) keep <- keep & dat$Year <= Rceattle[[i]]$data_list$endyr
    if (!is.null(maxyr)) keep <- keep & dat$Year <= maxyr
    dat <- dat[keep, , drop = FALSE]
    if (nrow(dat) == 0) next
    # Projection rows carry no observation, so exclude NAs as well as
    # non-positives; those rows keep the NA_real_ interval and draw no error bar.
    pos <- !is.na(dat$.obs) & dat$.obs > 0
    upr <- lwr <- rep(NA_real_, nrow(dat))
    # Interval on the scale the fleet is actually fitted on. A natural-scale
    # survey family (MVN / MVNORM / Normal) carries an ABSOLUTE sd, so the
    # lognormal form treats it as a log-sd: an sd of 150 on an index of ~100
    # draws a band of [0, 1e130] and the panel is unreadable. Catch is lognormal
    # throughout, so it always takes the log branch.
    #
    # The natural-scale interval is the untruncated normal's, clamped at zero.
    # For a "TruncatedNormal" fleet that is not the truncated quantile, so the
    # band is approximate wherever truncation carries real mass -- which is also
    # where sim_mod() warns.
    nat <- rep(FALSE, nrow(dat))
    if (kind == "index") {
      nat_all <- .index_rows_natural_scale(Rceattle[[i]]$data_list)
      if (length(nat_all)) nat <- nat_all[keep]
    }
    lg <- pos & !nat
    nt <- pos & nat
    upr[lg] <- stats::qlnorm(0.975, log(dat$.obs[lg]), dat$.sd[lg])
    lwr[lg] <- stats::qlnorm(0.025, log(dat$.obs[lg]), dat$.sd[lg])
    upr[nt] <- dat$.obs[nt] + stats::qnorm(0.975) * dat$.sd[nt]
    # An index cannot be negative; clamp rather than draw a bar through zero.
    lwr[nt] <- pmax(0, dat$.obs[nt] + stats::qnorm(0.025) * dat$.sd[nt])
    out[[length(out) + 1L]] <- data.frame(
      Model       = model_names_use[i],
      Fleet       = as.character(fc$Fleet_name[match(dat$Fleet_code,
                                                     fc$Fleet_code)]),
      Year        = abs(dat$Year),
      Observation = dat$.obs,
      Lower95     = lwr,
      Upper95     = upr,
      Predicted   = dat$.hat,
      stringsAsFactors = FALSE)
  }
  if (length(out) == 0L) {
    stop("No data to plot: the selected species/fleet(s) have no ",
         "likelihood observations (all rows are projections or excluded).",
         call. = FALSE)
  }
  df <- do.call(rbind, out)
  df$Model <- factor(df$Model, levels = unique(model_names_use))
  df$Fleet <- factor(df$Fleet,
                     levels = unique(fc$Fleet_name[match(codes, fc$Fleet_code)]))
  df
}

# Shared renderer for index/catch fit plots.
.render_fleet_fit <- function(df, ylab, error = TRUE, single_model = FALSE) {
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Year, colour = .data$Model,
                                        fill = .data$Model))
  if (error) {
    p <- p +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$Lower95,
                                          ymax = .data$Upper95),
                             width = 0, colour = "grey60") +
      ggplot2::geom_point(ggplot2::aes(y = .data$Observation),
                          shape = 21, fill = "white", colour = "grey30")
  }
  p <- p + ggplot2::geom_line(ggplot2::aes(y = .data$Predicted), linewidth = 1)
  p <- .rceattle_scale(
    p +
      ggplot2::facet_wrap(~ Fleet, scales = "free_y") +
      ggplot2::labs(x = "Year", y = ylab) +
      .rceattle_theme())
  if (single_model) p <- p + ggplot2::guides(colour = "none", fill = "none")
  p
}


#' Survey index fits
#'
#' Plots fitted index series: observed points with 95% observation intervals and
#' the model-predicted index, faceted by fleet. Every fleet holding
#' `index_data` is drawn, a fishery with a CPUE series as much as a survey.
#'
#' The interval follows the fleet's `Index_distribution`: lognormal quantiles on
#' a `"Lognormal"` fleet, and `observed +/- 1.96 sd` truncated at zero on a
#' natural-scale one (`"MVN"`, `"MVNORM"`, `"Normal"`, `"TruncatedNormal"`),
#' whose sd is absolute and in the units of the index.
#'
#' @param Rceattle A single [fit_mod()] object or a list of them (overlaid).
#' @param file Optional file stem; if given the figure is written to
#'   `<file>_survey_indices.png`.
#' @param model_names Legend labels for the models.
#' @param species Species (indices) to include. Default all.
#' @param incl_proj Include projection years.
#' @param width,height Saved figure size (inches).
#' @param error Draw observed points and error bars.
#' @param log Plot the index on the log scale.
#' @param line_col,right_adj,top_adj,single.plots Deprecated base-graphics
#'   arguments, retained for back-compatibility and ignored.
#' @return A `ggplot` object.
#' @export
plot_index <- function(Rceattle,
                       file = NULL,
                       model_names = NULL,
                       species = NULL,
                       incl_proj = FALSE,
                       width = 7,
                       height = 6.5,
                       error = TRUE,
                       log = FALSE,
                       line_col = NULL,
                       right_adj = 0,
                       top_adj = 0.05,
                       single.plots = FALSE) {
  Rceattle <- .as_plot_models(Rceattle)
  df <- .fleet_fit_df(Rceattle, kind = "index", species = species,
                      incl_proj = incl_proj, model_names = model_names)
  ylab <- "Index"
  if (log) {
    for (col in c("Observation", "Lower95", "Upper95", "Predicted")) {
      df[[col]] <- log(df[[col]])
    }
    ylab <- "log(Index)"
  }
  p <- .render_fleet_fit(df, ylab = ylab, error = error,
                         single_model = nlevels(df$Model) < 2L)
  .save_ggplot(p, file = file, suffix = "survey_indices",
               width = width, height = height)
}


#' Survey index fits on the log scale (deprecated)
#'
#' @description Deprecated. Use `plot_index(..., log = TRUE)`, which draws the
#'   survey index on the log scale.
#' @param ... Passed through to [plot_index()].
#' @return A `ggplot` object.
#' @keywords internal
#' @export
plot_logindex <- function(...) {
  .Deprecated("plot_index")
  plot_index(..., log = TRUE)
}


# Tidy the fleet x year catchability matrix into one row per (model, fleet, year).
# `index_q` is the realized q, whichever route set it (see ?plot_catchability).
.catchability_df <- function(Rceattle, species = NULL, spnames = NULL,
                             model_names = NULL) {
  labels <- .model_labels(Rceattle, model_names)
  fc     <- Rceattle[[1]]$data_list$fleet_control
  codes  <- .fleets_with_index(Rceattle[[1]]$data_list, fitted_only = FALSE)
  sp     <- .resolve_species(Rceattle, species, spnames)$index
  codes  <- codes[fc$Species[match(codes, fc$Fleet_code)] %in% sp]
  if (length(codes) == 0L) {
    stop("No fleet to plot: the selected species have no fleet carrying ",
         "`index_data`, so no catchability is estimated for them.",
         call. = FALSE)
  }

  out <- list()
  for (i in seq_along(Rceattle)) {
    q <- Rceattle[[i]]$quantities$index_q
    if (is.null(q)) next
    dl   <- Rceattle[[i]]$data_list
    # Fleets are identified by model 1's fleet_control. Models whose fleets
    # differ would be drawn under the wrong names, silently, so say so instead.
    if (!identical(as.character(dl$fleet_control$Fleet_name),
                   as.character(fc$Fleet_name))) {
      stop("Model ", i, " has different fleets from model 1, so their ",
           "catchabilities cannot be drawn on shared panels.", call. = FALSE)
    }
    # index_q is hindcast-only, so its columns are styr:endyr in order.
    yrs  <- dl$styr:(dl$styr + ncol(q) - 1L)
    rows <- match(codes, fc$Fleet_code)
    out[[i]] <- data.frame(
      Model = labels[i],
      Fleet = rep(fc$Fleet_name[rows], each = length(yrs)),
      Year  = rep(yrs, times = length(rows)),
      value = as.vector(t(q[rows, , drop = FALSE])),
      stringsAsFactors = FALSE)
  }
  if (length(out) == 0L) {
    stop("No catchability to plot: `quantities$index_q` is absent from every ",
         "model. It is written by a fitted model; a `data_list` alone has none.",
         call. = FALSE)
  }
  df <- do.call(rbind, out)
  df$Model <- factor(df$Model, levels = unique(labels))
  df$Fleet <- factor(df$Fleet, levels = fc$Fleet_name[match(codes, fc$Fleet_code)])
  df
}


#' Survey catchability over time
#'
#' Plots the fitted survey catchability `q` by year, faceted by fleet. Every
#' fleet holding `index_data` is drawn, a fishery with a CPUE series as much as
#' a survey.
#'
#' The value plotted is the realized catchability the model scaled the survey by,
#' whichever route produced it: an estimated or fixed mean, annual deviations
#' under `Time_varying_q`, an environmental regression (`Catchability =
#' "Environmental"`), a `q` linkage such as `rw(1 | Year)` or `ar1(1 | Year)`, or
#' the closed-form value under `"Analytical"` / `"AnalyticalArith"`. A fleet with
#' a time-invariant catchability draws a flat line, which is the honest picture
#' of a constant `q` rather than an empty panel.
#'
#' **Hindcast years only.** The model applies catchability over the hindcast and
#' does not project it, so the series stops at `endyr` regardless of the
#' projection horizon. There is no `incl_proj`.
#'
#' `line_col` and `lty` separate the **models**, so overlaying several fits
#' compares their catchability series panel by panel.
#'
#' @inheritParams rceattle-plot-args
#' @param log Plot catchability on the log scale.
#' @param mse Is `Rceattle` an MSE object? Its operating models are drawn.
#' @return A `ggplot` object.
#' @examples
#' \donttest{
#' data(BS2017SS)
#' fit <- Rceattle::fit_mod(data_list = BS2017SS, msmMode = 0)
#' plot_catchability(fit)
#' plot_catchability(fit, log = TRUE)
#' }
#' @export
plot_catchability <- function(Rceattle,
                              file = NULL,
                              model_names = NULL,
                              species = NULL,
                              spnames = NULL,
                              mse = FALSE,
                              minyr = NULL,
                              maxyr = NULL,
                              incl_mean = FALSE,
                              log = FALSE,
                              line_col = NULL,
                              lwd = 3,
                              lty = 1,
                              width = 7,
                              height = 6.5,
                              right_adj = 0,
                              top_adj = 0.05,
                              single.plots = FALSE) {
  Rceattle <- .as_plot_models(Rceattle, mse = mse, OM = TRUE)
  df <- .catchability_df(Rceattle, species = species, spnames = spnames,
                         model_names = model_names)
  df <- .rce_year_filter(df, minyr = minyr, maxyr = maxyr)

  ylab <- "Catchability (q)"
  if (log) {
    df$value <- log(df$value)
    ylab <- "log(q)"
  }

  n_mod <- nlevels(df$Model)
  lines <- .rce_line_params(lwd = lwd, lty = lty,
                            lwd_by = "model", lty_by = "model",
                            lwd_n = n_mod, lty_n = n_mod)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Year, y = .data$value,
                                        colour = .data$Model)) +
    do.call(ggplot2::geom_line, lines$args) +
    .rce_mean_line(df, incl_mean = incl_mean, by = c("Fleet", "Model"),
                   value = "value",
                   hind_endyr = .rce_model_endyr(Rceattle, levels(df$Model)),
                   colour_by = "Model") +
    ggplot2::facet_wrap(~ Fleet, scales = "free_y") +
    ggplot2::labs(x = "Year", y = ylab) +
    .rceattle_theme()
  for (s in lines$scales) p <- p + s
  p <- .rceattle_scale(p, line_col = line_col)
  if (n_mod < 2L) p <- p + ggplot2::guides(colour = "none", fill = "none")

  .save_ggplot(p, file = file, suffix = "catchability",
               width = width, height = height)
}


#' Fishery catch fits
#'
#' Plots fitted fishery catch: observed points with lognormal 95% intervals and
#' the model-predicted catch, faceted by fishery fleet. For MSE objects the
#' projection period is summarized with 50% / 95% ribbons across simulations.
#'
#' @inheritParams plot_index
#' @param maxyr Last year to plot.
#' @param mse Is `Rceattle` an MSE object (operating models)?
#' @param fleets Fishery fleets to include (indices into the sorted fleet
#'   codes). Default all.
#' @param alpha,lwd,ymax Deprecated base-graphics arguments, ignored.
#' @return A `ggplot` object.
#' @export
plot_catch <- function(Rceattle,
                       file = NULL,
                       model_names = NULL,
                       species = NULL,
                       incl_proj = FALSE,
                       width = 7,
                       height = 6.5,
                       error = TRUE,
                       maxyr = NULL,
                       mse = FALSE,
                       fleets = NULL,
                       line_col = NULL,
                       right_adj = 0,
                       top_adj = 1.2,
                       single.plots = FALSE,
                       alpha = 0.4,
                       lwd = 2,
                       ymax = NULL) {
  if (mse) {
    Rceattle <- .as_plot_models(Rceattle, mse = TRUE, OM = TRUE)
    incl_proj <- TRUE
  } else {
    Rceattle <- .as_plot_models(Rceattle)
  }

  df <- .fleet_fit_df(Rceattle, kind = "catch", incl_proj = incl_proj,
                      maxyr = maxyr, model_names = model_names)

  if (!is.null(fleets)) {
    lev <- levels(df$Fleet)
    df <- df[df$Fleet %in% lev[fleets], , drop = FALSE]
    df$Fleet <- droplevels(df$Fleet)
  }

  if (mse) {
    # Summarise predicted catch across operating models per fleet/year
    agg <- stats::aggregate(Predicted ~ Fleet + Year, data = df,
      FUN = function(v) c(mean = mean(v),
                          l95 = stats::quantile(v, 0.025, names = FALSE),
                          u95 = stats::quantile(v, 0.975, names = FALSE),
                          l50 = stats::quantile(v, 0.25, names = FALSE),
                          u50 = stats::quantile(v, 0.75, names = FALSE)))
    mse_df <- data.frame(Fleet = agg$Fleet, Year = agg$Year,
                         agg$Predicted)
    p <- ggplot2::ggplot(mse_df, ggplot2::aes(x = .data$Year)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$l95, ymax = .data$u95),
                           alpha = 0.25, fill = "grey40") +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$l50, ymax = .data$u50),
                           alpha = 0.4, fill = "grey40") +
      ggplot2::geom_line(ggplot2::aes(y = .data$mean), linewidth = 1) +
      ggplot2::facet_wrap(~ Fleet, scales = "free_y") +
      ggplot2::labs(x = "Year", y = "Catch") +
      .rceattle_theme()
    return(.save_ggplot(p, file = file, suffix = "fishery_catch",
                        width = width, height = height))
  }

  p <- .render_fleet_fit(df, ylab = "Catch", error = error,
                         single_model = nlevels(df$Model) < 2L)
  .save_ggplot(p, file = file, suffix = "fishery_catch",
               width = width, height = height)
}


#' Survey index residuals
#'
#' Plots survey index residuals by year, faceted by survey fleet
#' (`residual_type = "pearson"`, the default), or one-step-ahead (OSA) residual
#' diagnostics for a single fit (`residual_type = "osa"`, via [osa_residuals()] /
#' [plot.rceattle_osa()]).
#'
#' The residual is **observed minus predicted**, taken on the scale the fleet is
#' fitted on and read from its `Index_distribution`:
#' `log(observed) - log(predicted)` for a log-scale (`"Lognormal"`) fleet, and
#' `observed - predicted` for a natural-scale one (`"MVN"`, `"MVNORM"`,
#' `"Normal"`, `"TruncatedNormal"`), whose sd is absolute. A positive residual
#' means the survey saw more than the model predicted. Panels hold different
#' units where a model mixes the two families, which is why the y scale is free
#' per fleet.
#'
#' Before 5.9.0 this plotted `predicted - observed`, the negative of what
#' [residuals.Rceattle()] returns for the same fleet. Plots made with an earlier
#' version are mirrored about zero relative to these.
#'
#' @inheritParams plot_index
#' @param residual_type `"pearson"` (index residuals on the fleet's own scale)
#'   or `"osa"`.
#' @return A `ggplot` object.
#' @export
plot_indexresidual <- function(Rceattle,
                               file = NULL,
                               model_names = NULL,
                               species = NULL,
                               width = 7,
                               height = 6.5,
                               residual_type = c("pearson", "osa"),
                               line_col = NULL,
                               right_adj = 0,
                               top_adj = 0.05,
                               incl_proj = FALSE,
                               single.plots = FALSE) {
  residual_type <- match.arg(residual_type)
  if (residual_type == "osa") {
    fit <- if (inherits(Rceattle, "Rceattle")) Rceattle else Rceattle[[1]]
    osa <- osa_residuals(fit, source = c("index", "catch"))
    return(plot(osa))
  }

  Rceattle <- .as_plot_models(Rceattle)
  model_names_use <- .model_labels(Rceattle, model_names)
  fc <- Rceattle[[1]]$data_list$fleet_control
  if (is.null(species)) species <- seq_len(Rceattle[[1]]$data_list$nspp)

  out <- list()
  for (i in seq_along(Rceattle)) {
    dat <- Rceattle[[i]]$data_list$index_data
    dat$.hat <- Rceattle[[i]]$quantities$index_hat
    # Drop prediction-only rows (Year < 0); they are not in the likelihood and
    # have no residual to plot. Projection years are excluded unless asked for,
    # matching `.fleet_fit_df()` -- a projection row's "observation" is whatever
    # placeholder its workbook carries, which plots as a spurious residual.
    # Same fleet set as plot_index(). Without this an Off fleet's index rows were
    # drawn as residuals indistinguishably from fitted ones -- GOA2018SS fleet 7
    # is the case: 9 index rows on a fleet the likelihood never reads.
    keep <- dat$Species %in% species & dat$Year > 0 &
      dat$Fleet_code %in% .fleets_with_index(Rceattle[[i]]$data_list,
                                             fitted_only = FALSE)
    if (!incl_proj) keep <- keep & dat$Year <= Rceattle[[i]]$data_list$endyr
    dat <- dat[keep, , drop = FALSE]
    if (nrow(dat) == 0) next
    out[[length(out) + 1L]] <- data.frame(
      Model    = model_names_use[i],
      Fleet    = as.character(fc$Fleet_name[match(dat$Fleet_code,
                                                  fc$Fleet_code)]),
      Year     = abs(dat$Year),
      # OBSERVED minus PREDICTED, matching residuals.Rceattle() and the usual
      # assessment convention (WHAM, SS3): a positive residual means the survey
      # saw more than the model predicted.
      #
      # Log ratio for a log-scale fleet, plain difference for a natural-scale one
      # -- log() of a natural-scale residual is not a residual. Assigned by
      # subset rather than with ifelse(), which evaluates both arms and so would
      # take log() of a natural-scale observation it then discards.
      Residual = {
        nat <- .index_rows_natural_scale(Rceattle[[i]]$data_list)
        nat <- if (length(nat)) nat[keep] else rep(FALSE, nrow(dat))
        res <- rep(NA_real_, nrow(dat))
        res[nat]  <- dat$Observation[nat] - dat$.hat[nat]
        res[!nat] <- log(dat$Observation[!nat]) - log(dat$.hat[!nat])
        res
      },
      stringsAsFactors = FALSE)
  }
  if (length(out) == 0L) {
    stop("No index residuals to plot: the selected species have no ",
         "likelihood observations (all rows are projections or excluded).",
         call. = FALSE)
  }
  df <- do.call(rbind, out)
  df$Model <- factor(df$Model, levels = unique(model_names_use))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Year, y = .data$Residual,
                                        colour = .data$Model)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
    ggplot2::geom_segment(ggplot2::aes(xend = .data$Year, yend = 0),
                          position = ggplot2::position_dodge(width = 0.6)) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.6))
  p <- .rceattle_scale(
    p +
      # Panels are per fleet with a free y scale, which is what lets the two
      # scales share a figure: a log-scale fleet's panel is dimensionless, a
      # natural-scale fleet's is in the units of the index.
      ggplot2::facet_wrap(~ Fleet, scales = "free_y") +
      ggplot2::labs(x = "Year", y = "Index residual") +
      .rceattle_theme(),
    aesthetics = "colour")
  if (nlevels(df$Model) < 2L) p <- p + ggplot2::guides(colour = "none")
  .save_ggplot(p, file = file, suffix = "index_residuals",
               width = width, height = height)
}
