#' Make a vector of colors.
#'
#' A subset of rich.colors by Arni Magnusson from the gplots package, with the
#' addition of alpha transparency (which is now available in the gplots version
#' as well)
#'
#'
#' @param n Number of colors to generate.
#' @param alpha Alpha transparency value for all colors in vector. Value is
#' passed to rgb function.
#' @author Arni Magnusson, Ian Taylor
rich.colors.short <- function(n,alpha=1){
  x <- seq(0, 1, length = n)
  r <- 1/(1 + exp(20 - 35 * x))
  g <- pmin(pmax(0, -0.8 + 6 * x - 5 * x^2), 1)
  b <- dnorm(x, 0.25, 0.15)/max(dnorm(x, 0.25, 0.15))
  rgb.m <- matrix(c(r, g, b), ncol = 3)
  apply(rgb.m, 1, function(v) rgb(v[1], v[2], v[3], alpha = alpha))
}

#' Plot time-series
#'
#' @description Function that plots the time-series (SSB/B/R/Depletion) 95% CI trends as estimated from Rceattle
#'
#' @details
#' # Units
#'
#' The model holds numbers-at-age in **thousands** and weight-at-age in
#' **kg**, so every biomass series (`biomass`, `ssb`, `exploitable_biomass`)
#' comes out of the model in **mt** and recruitment comes out in **thousands of
#' fish**. For display these are divided by 1e6 (million mt) and 1e3 (millions
#' of recruits) respectively; depletion is a ratio and is not rescaled. Supply
#' the model's inputs on that convention, catch and index in mt, weight-at-age
#' in kg, or the axis labels will not describe what is plotted.
#'
#' # Confidence intervals
#'
#' `add_ci = TRUE` draws a 95% interval. Every strictly positive series takes it
#' on the log scale as `exp(log(x) +/- qnorm(0.975) * sd_log)`. These quantities are
#' built multiplicatively, so `log(x)` is close to linear in the estimated
#' parameters where `x` is exponential in them: the interval is both a better
#' linearization and right-skewed, and it cannot cross zero the way a symmetric
#' natural-scale interval does for weak year classes and depleted stocks.
#'
#' `sd_log` comes from the model's own `log_biomass` / `log_ssb` / `log_R` where
#' those are reported, and is otherwise recovered as `sd(x) / x`, the delta
#' method's own identity, which matches the reported values to machine
#' precision. That covers `exploitable_biomass` and the two depletions, which
#' cannot be reported on the log scale (`exploitable_biomass` is identically 0
#' without projection F), and models fit before the `log_*` series existed. A
#' non-positive or unreported value keeps the symmetric interval.
#'
#' @inheritParams rceattle-plot-args
#' @param output derived quantity of interest: recruitment, biomass, ssb, depletion, or ssb_depletion. Uses same name as ".cpp" file.
#' @param ylab Y-axis label. `NULL` (default) derives one from `output` and the
#'   model's `minage`.
#' @param save Write the plotted series to CSV alongside the figure.
#' @param mse Is an MSE object from \code{\link{load_mse}} or \code{\link{run_mse}}
#' @param OM if mse == TRUE, use the OM (TRUE) or EM (FALSE) for plotting?
#' @param mod_avg TRUE/FALSE
#' @param zero_y Anchor the y-axis at zero. TRUE for the absolute series
#'   (biomass, SSB, recruitment), where a truncated axis exaggerates change;
#'   FALSE for the depletions, which are already on a relative scale.
#' @param reference A model to overlay for comparison, drawn at 1.5x `lwd` and
#'   labelled "Reference". It takes the next colour from the palette, or black
#'   if `line_col` is supplied.
#' @param ref_lines Internal. `function(models, sp_sel)` returning per-species
#'   reference-point layers, added before the figure is saved. Supplied by
#'   [plot_f()] and [plot_depletionSSB()] through `.ts_wrapper()`; `sp_sel` is
#'   the species resolution the panels were built from.
#' @param suffix Internal. Overrides the `<file>_<suffix>.png` stem, which
#'   otherwise names the series.
#'
#' @importFrom stats quantile
#' @export
#'
#' @return Returns and saves a figure with the population trajectory.
plot_timeseries <- function(Rceattle,
                            output = "biomass",
                            ylab = NULL,
                            file = NULL,
                            model_names = NULL,
                            line_col = NULL,
                            species = NULL,
                            spnames = NULL,
                            add_ci = FALSE,
                            lwd = 3,
                            save = FALSE,
                            legend.pos = "topright",
                            right_adj = 0,
                            width = 7,
                            height = 6.5,
                            minyr = NULL,
                            maxyr = NULL,
                            incl_proj = FALSE,
                            mod_cex = 1,
                            lty = rep(1, length(Rceattle)),
                            alpha = 0.4,
                            mse = FALSE,
                            OM = TRUE,
                            reference = NULL,
                            zero_y = FALSE,
                            ref_lines = NULL,
                            suffix = NULL,
                            mod_avg = rep(FALSE, length(Rceattle))) {

  .rce_check_alpha(alpha)

  ## Model object manipulation ----
  # Convert mse object to Rceattle list
  if(mse){
    if(OM){
      Rceattle <- lapply(Rceattle, function(x) x$OM)
    } else {
      Rceattle <- lapply(Rceattle, function(x) x$EM[[length(x$EM)]])
    }
    nmse = length(Rceattle)
    add_ci = TRUE
    incl_proj = TRUE
  }


  # Convert single one into a list
  if(inherits(Rceattle, "Rceattle")){
    Rceattle <- list(Rceattle)
  }


  ## Line width / colour ----
  lwd <- rep(lwd, length(Rceattle))
  # line_col stays NULL unless the user supplies colours, in which case it falls
  # through to the default colourblind-safe scale (.rceattle_scale) at render
  # time. When supplied, recycle to one colour per model.
  if (!is.null(line_col)) {
    line_col <- rep(line_col, length.out = length(Rceattle))
  }

  ## Add reference model
  if(!is.null(reference)){
    Rceattle <- c(Rceattle, list(reference))
    mod_avg = c(mod_avg, FALSE)
    lty = c(lty, 1)
    if (!is.null(line_col)) line_col <- c(line_col, "black")
    lwd <- c(lwd, lwd[1] * 1.5)
    # Give the reference its own label so .model_labels() does not recycle the
    # user's model_names (which would collapse the reference into the first
    # model's factor level and drop its distinct colour / width / type).
    if (!is.null(model_names)) model_names <- c(model_names, "Reference")
  }

  # Extract years
  Endyrs <-  sapply(Rceattle, function(x) x$data_list$endyr)
  years <- lapply(Rceattle, function(x) x$data_list$styr:x$data_list$endyr)
  if(incl_proj){
    years <- lapply(Rceattle, function(x) x$data_list$styr:x$data_list$projyr)
  }
  if(!is.null(maxyr)){
    years <- lapply(Rceattle, function(x) x$data_list$styr:min(c(maxyr, x$data_list$projyr)))
  }

  max_endyr <- max(unlist(Endyrs), na.rm = TRUE)
  nyrs_vec <- sapply(years, length)
  nyrs <- max(nyrs_vec)
  if(is.null(minyr)){minyr <- min((sapply(years, min)))}
  if(is.null(maxyr)){maxyr <- max((sapply(years, max)))}

  # Extract species information
  nspp <- Rceattle[[1]]$data_list$nspp
  minage <- Rceattle[[1]]$data_list$minage
  estDynamics <- Rceattle[[1]]$data_list$estDynamics

  # `species` selects (indices, names, a logical mask, or "all"); `spnames`
  # labels. Shared with every other plotter so the pair means one thing.
  sp_sel  <- .resolve_species(Rceattle, species, spnames)
  spp     <- sp_sel$index
  spnames <- sp_sel$spnames

  # Built here rather than in .ts_wrapper() because the age it names is
  # per-species data. An explicit `ylab` always wins.
  if (is.null(ylab)) ylab <- .rce_ts_ylab(output, minage)


  ## Get quantity of interest ----
  # Save objects
  quantity <-
    array(NA, dim = c(nspp, nyrs,  length(Rceattle)))
  quantity_sd <-
    array(NA, dim = c(nspp, nyrs,  length(Rceattle)))
  # Delta-method SD of log(quantity), when the model ADREPORTed a matching
  # `log_<output>` series. Distinct from log_quantity_sd below, which is the
  # across-sample SD used by the model-averaging path.
  quantity_log_sd <-
    array(NA, dim = c(nspp, nyrs,  length(Rceattle)))
  log_quantity_sd <-
    array(NA, dim = c(nspp, nyrs,  length(Rceattle)))
  log_quantity_mu <-
    array(NA, dim = c(nspp, nyrs,  length(Rceattle)))


  for (i in 1:length(Rceattle)) {
    # - Get quantities
    quantity[, 1:nyrs_vec[i] , i] <- Rceattle[[i]]$quantities[[output]][,1:nyrs_vec[i]]

    # Get SD of quantity
    if(add_ci & !mse){
      # The series is ADREPORTed over the whole modelled period, flattened
      # column-major with species varying fastest, so the first
      # nspp * nyrs_hindcast values are exactly the years being plotted.
      # n_total states that full shape, so the slice is checked rather than
      # assumed -- an unexpected block length draws no interval instead of one
      # belonging to different cells.
      n_need  <- nyrs_vec[i] * nspp
      n_total <- ncol(Rceattle[[i]]$quantities[[output]]) * nspp
      sd_vals <- .rce_series_sd(Rceattle[[i]], output, n_need, n_total)
      if (is.null(sd_vals)) {
        # Leave the SD array as NA and draw no ribbon. Only warn for a series
        # the model could have reported: F_spp and the other REPORT()-only
        # quantities never carry standard errors, so warning on them would fire
        # on every fit.
        .rce_no_ci(TRUE, paste0("`", output, "` in model ", i),
                   paste("this fit carries no standard errors for it",
                         "(not reported by the model, or getsd = FALSE)"),
                   warn = .rce_quantity_adreport(output))
      } else {
        quantity_sd[, 1:nyrs_vec[i], i] <-
          replace(quantity_sd[, 1:nyrs_vec[i], i], values = sd_vals)

        # Log-scale SD where the model reported one; the rest are recovered
        # from sd/x below.
        log_vals <- .rce_series_sd(Rceattle[[i]], paste0("log_", output),
                                   n_need, n_total)
        if (!is.null(log_vals)) {
          quantity_log_sd[, 1:nyrs_vec[i], i] <-
            replace(quantity_log_sd[, 1:nyrs_vec[i], i], values = log_vals)
        }
      }
    }

    # - Model average
    if(mod_avg[i]){
      log_quantity_sd[,  1:nyrs_vec[i], i] <- apply(Rceattle[[i]]$asymptotic_samples[quantity][,1:nyrs_vec[i],], c(1,2), function(x) sd(as.vector(log(x))))
      log_quantity_mu[,  1:nyrs_vec[i], i] <- apply(Rceattle[[i]]$asymptotic_samples[quantity][,1:nyrs_vec[i],], c(1,2), function(x) mean(as.vector(log(x))))
    }
  }

  ## Get confidence intervals ----
  if(!mse){
    # - Single model. .rce_ci_bounds() is shared with as.data.frame.Rceattle(),
    #   so the figure and the table cannot disagree: a strictly positive series
    #   takes exp(log(x) +/- z * sd_log), recovering sd(log x) = sd(x)/x where
    #   the model reported no log_* row, and everything else stays symmetric.
    ci95 <- .rce_ci_bounds(quantity, quantity_sd, quantity_log_sd,
                           z = stats::qnorm(0.975))
    ci50 <- .rce_ci_bounds(quantity, quantity_sd, quantity_log_sd,
                           z = stats::qnorm(0.75))

    quantity_upper95 <- array(ci95$upr, dim = dim(quantity))
    quantity_lower95 <- array(ci95$lwr, dim = dim(quantity))
    quantity_upper50 <- array(ci50$upr, dim = dim(quantity))
    quantity_lower50 <- array(ci50$lwr, dim = dim(quantity))
  } else {
    # - MSE objects
    # -- Get quantiles and mean across simulations
    quantity_upper95 <- apply( quantity[,,1:nmse], c(1,2), function(x) quantile(x, probs = 0.975) )
    quantity_lower95 <- apply( quantity[,,1:nmse], c(1,2), function(x) quantile(x, probs = 0.025) )
    quantity_upper50 <- apply( quantity[,,1:nmse], c(1,2), function(x) quantile(x, probs = 0.75) )
    quantity_lower50 <- apply( quantity[,,1:nmse], c(1,2), function(x) quantile(x, probs = 0.25) )

    # -- Put back in array for indexing below
    if(is.null(reference)){
      quantity <- array(apply( quantity[,,1:nmse], c(1,2), mean ), dim = c(nspp, nyrs,  1))
      quantity_upper95 <- array(quantity_upper95, dim = c(nspp, nyrs,  1))
      quantity_lower95 <- array(quantity_lower95, dim = c(nspp, nyrs,  1))
      quantity_upper50 <- array(quantity_upper50, dim = c(nspp, nyrs,  1))
      quantity_lower50<- array(quantity_lower50, dim = c(nspp, nyrs,  1))
    } else {
      quantity_upper95 <- array(c(quantity_upper95, quantity[,,nmse+1]), dim = c(nspp, nyrs,  2))
      quantity_lower95 <- array(c(quantity_lower95, quantity[,,nmse+1]), dim = c(nspp, nyrs,  2))
      quantity_upper50 <- array(c(quantity_upper50, quantity[,,nmse+1]), dim = c(nspp, nyrs,  2))
      quantity_lower50<- array(c(quantity_lower50, quantity[,,nmse+1]), dim = c(nspp, nyrs,  2))
      quantity <- array(c(apply( quantity[,,1:nmse], c(1,2), mean ), quantity[,,nmse+1]), dim = c(nspp, nyrs,  2))
    }
  }

  # - Model Average
  for (i in 1:length(Rceattle)) {
    if(mod_avg[i]){
      quantity[,,i] <- qlnorm(0.5, meanlog = log_quantity_mu[,,i], sdlog = log_quantity_sd[,,i])
      quantity_upper95[,,i] <- qlnorm(0.975, meanlog = log_quantity_mu[,,i], sdlog = log_quantity_sd[,,i])
      quantity_lower95[,,i] <- qlnorm(0.025, meanlog = log_quantity_mu[,,i], sdlog = log_quantity_sd[,,i])
    }
  }

  ## Rescale ----
  # Numbers-at-age are in thousands and weight-at-age in kg, so biomass is in mt
  # and recruitment in thousands of fish. Depletion is a ratio and is not
  # rescaled.
  rescale <- .RCE_TS_RESCALE[output]
  if (!is.na(rescale)) {
    quantity <- quantity / rescale
    quantity_upper95 <- quantity_upper95 / rescale
    quantity_lower95 <- quantity_lower95 / rescale
    quantity_upper50 <- quantity_upper50 / rescale
    quantity_lower50 <- quantity_lower50 / rescale
  }


  # Legend / CSV column labels for the models, defaulting to "Model 1", ...
  model_names_use <- .model_labels(Rceattle, model_names)

  ## Save ----
  if (save) {
    if (is.null(file)) {
      stop("`save = TRUE` needs a `file` stem to write to.", call. = FALSE)
    }
    # The CSV holds what the figure holds, so it is cut to the same window.
    save_idx <- which(years[[1]] >= minyr & years[[1]] <= maxyr)
    for (i in spp) {
      dat_new <- data.frame(Year = years[[1]][save_idx])
      for (j in seq_along(Rceattle)) {
        block <- data.frame(quantity[i, save_idx, j],
                            quantity_lower95[i, save_idx, j],
                            quantity_upper95[i, save_idx, j])
        colnames(block) <- paste0(model_names_use[j],
                                  c("", "_lower95", "_upper95"))
        dat_new <- cbind(dat_new, block)
      }
      utils::write.csv(dat_new, row.names = FALSE,
                       file = paste0(file, "_", output, "_",
                                     gsub("[^A-Za-z0-9]+", "_", spnames[i]),
                                     ".csv"))
    }
  }


  ## Assemble tidy data ----
  # Build one row per (species, year, model) holding exactly the plotted
  # quantities (the same arrays, reshaped) so the returned ggplot's `$data` can
  # be tested against the model's source quantities.
  df_list <- list()
  for (k in seq_len(dim(quantity)[3])) {
    yk  <- years[[k]]
    nyk <- length(yk)
    for (sp in spp) {
      df_list[[length(df_list) + 1L]] <- data.frame(
        Species = spnames[sp],
        Year    = yk,
        Model   = model_names_use[k],
        value   = quantity[sp, seq_len(nyk), k],
        lower95 = quantity_lower95[sp, seq_len(nyk), k],
        upper95 = quantity_upper95[sp, seq_len(nyk), k],
        lower50 = quantity_lower50[sp, seq_len(nyk), k],
        upper50 = quantity_upper50[sp, seq_len(nyk), k],
        show_ci = (estDynamics[sp] == 0),
        stringsAsFactors = FALSE
      )
    }
  }
  plot_df <- do.call(rbind, df_list)
  # Narrow to the requested window, rather than clipping the axis over it. The
  # species panels use scales = "free_y", so a scale trained on the hidden years
  # squeezes the window the caller asked for into the bottom of the panel.
  # `maxyr` has already truncated `years` above; `minyr` is applied here.
  plot_df <- .rce_year_filter(plot_df, minyr, maxyr)
  plot_df$Species <- factor(plot_df$Species, levels = spnames[spp])
  plot_df$Model   <- factor(plot_df$Model, levels = unique(model_names_use))

  ## Render ----
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data$Year, y = .data$value,
                 colour = .data$Model, fill = .data$Model))

  # Credible-interval ribbons (only for estimated species)
  if (add_ci) {
    ci_df <- plot_df[plot_df$show_ci, , drop = FALSE]
    if (nrow(ci_df) > 0) {
      p <- p + ggplot2::geom_ribbon(
        data = ci_df,
        ggplot2::aes(ymin = .data$lower95, ymax = .data$upper95),
        alpha = alpha / 2, colour = NA)
      if (mse) {
        p <- p + ggplot2::geom_ribbon(
          data = ci_df,
          ggplot2::aes(ymin = .data$lower50, ymax = .data$upper50),
          alpha = alpha, colour = NA)
      }
    }
  }

  # Line width / type: honour per-model lwd / lty, keyed to Model when they
  # vary so the uniform default adds no extra legend. The lwd / 3 bridge lives
  # in the helper, which every other plotter now shares.
  .rce_check_line_col(line_col, nlevels(plot_df$Model), "models")
  line_p <- .rce_line_params(lwd = lwd, lty = lty,
                             lwd_by = "Model", lty_by = "Model")
  p <- .rce_add_line(p, line_p)

  # Projection divider at the terminal hindcast year
  p <- p + .rce_proj_divider(Rceattle, incl_proj, minyr, maxyr)

  # A single-species model needs no species strip: the facet label would name
  # the only thing on the plot. Faceting is what carries that label, so drop it
  # rather than blanking the strip and leaving its gap.
  if (length(spp) > 1L) {
    p <- p + ggplot2::facet_wrap(~ Species, scales = "free_y", ncol = 1,
                                 strip.position = "top")
  }
  p <- p +
    ggplot2::labs(x = "Year", y = ylab) +
    .rceattle_theme()
  # A quantity on an absolute scale is read against zero, so an axis starting at
  # the series minimum overstates how much it moves. Applied per panel, so
  # `scales = "free_y"` still fits each species.
  if (isTRUE(zero_y)) p <- p + ggplot2::expand_limits(y = 0)
  # `line_col = NULL` keeps the default Okabe-Ito palette; supplied colours are
  # applied to the model levels in plotting order.
  p <- .rceattle_scale(p, line_col = line_col)

  # A single model needs no colour legend
  if (nlevels(plot_df$Model) < 2L) {
    p <- p + ggplot2::guides(colour = "none", fill = "none",
                             linewidth = "none", linetype = "none")
  }

  # Per-series reference points (F, depletion), supplied by the wrapper rather
  # than switched on `output` here. They are handed the same `sp_sel` the panels
  # were built from: resolving `species` a second time is what put plot_f()'s
  # Ftarget lines in the wrong panel.
  if (!is.null(ref_lines)) p <- p + ref_lines(Rceattle, sp_sel)

  .save_ggplot(p, file = file,
               suffix = if (is.null(suffix)) paste0(output, "_trajectory")
                        else suffix,
               width = width, height = height)
}


# Display divisor per derived series; a series absent here is plotted unscaled.
# Any entry added needs its matching unit in .rce_ts_ylab(): the two drifting
# apart is how recruitment came to be plotted in thousands under a "million"
# label.
.RCE_TS_RESCALE <- c(
  biomass             = 1e6,   # mt              -> million mt
  ssb                 = 1e6,   # mt              -> million mt
  exploitable_biomass = 1e6,   # mt              -> million mt
  R                   = 1e3    # thousands of fish -> millions of recruits
)


# Default y-axis label for a derived series. The minimum age is per-species
# data, so it is interpolated at plot time; where species disagree on minage the
# age prefix is dropped, since one shared y axis cannot name two ages.
.rce_ts_ylab <- function(output, minage) {
  # Only recruitment names an age -- it is an age class. SSB is female spawning
  # biomass, not the minage+ stock, so it carries no age prefix.
  ages <- unique(as.integer(minage[!is.na(minage)]))
  age  <- if (length(ages) == 1L) paste0("Age-", ages, " ") else ""
  switch(output,
         biomass             = "Biomass (million mt)",
         ssb                 = "SSB (million mt)",
         exploitable_biomass = "Max exploitable biomass (million mt)",
         R                   = paste0(age, "recruits (millions)"),
         ssb_depletion       = "SSB depletion",
         biomass_depletion   = "Biomass depletion",
         F_spp               = "Fishing mortality (F)",
         output)
}


# Per-species reference lines: target in blue, limit in red, dashed. Shared by
# the F and depletion hooks below, which differ only in where the two values
# come from.
#
# `sp_sel` is the resolution plot_timeseries() built its panels from, so the
# lines are keyed to the same species and labelled the same way. Both vectors
# are indexed by the species index, never by the raw `species` argument.
.rce_ref_hlines <- function(sp_sel, target, limit) {
  if (is.null(target) || is.null(limit)) return(NULL)
  target <- as.numeric(target)
  limit  <- as.numeric(limit)
  if (max(sp_sel$index) > min(length(target), length(limit))) return(NULL)
  ref_df <- data.frame(
    Species = factor(sp_sel$labels, levels = sp_sel$labels),
    target  = target[sp_sel$index],
    limit   = limit[sp_sel$index])
  # No line where Ftarget or Flimit is NA.
  list(
    ggplot2::geom_hline(data = ref_df[!is.na(ref_df$target), , drop = FALSE], inherit.aes = FALSE,
                        ggplot2::aes(yintercept = .data$target),
                        colour = "blue", linetype = 2),
    ggplot2::geom_hline(data = ref_df[!is.na(ref_df$limit), , drop = FALSE], inherit.aes = FALSE,
                        ggplot2::aes(yintercept = .data$limit),
                        colour = "red", linetype = 2))
}


# Ptarget / Plimit for the depletion plots.
#
# One horizontal line per panel can only show one model's value, so it reads the
# first model, and the axis label says so.
#
# Take the value per model and index it directly. Collecting the models into a
# models x species matrix and subsetting that with species indices flattens
# column-major, so on an overlay whose models carry different Ptarget, species 2
# is drawn from model 2's species 1. It is invisible whenever the models share
# their reference points, which is the normal case.
.depletion_reference_lines <- function(models, sp_sel) {
  .rce_ref_hlines(sp_sel, models[[1]]$data_list$Ptarget,
                  models[[1]]$data_list$Plimit)
}


# Ftarget / Flimit for plot_f(). These are estimated, so they come from
# `quantities` rather than `data_list`.
.f_reference_lines <- function(models, sp_sel) {
  .rce_ref_hlines(sp_sel, models[[1]]$quantities$Ftarget,
                  models[[1]]$quantities$Flimit)
}


# Build a plot_timeseries() wrapper that pins the derived series (`output`) for
# one quantity while exposing plot_timeseries()'s full argument list unchanged.
# The timeseries plotters below differ only in that string and, for two of them,
# a set of reference lines, so they share one body through this factory. The
# y-axis label is resolved inside plot_timeseries() from .rce_ts_ylab(), which
# needs the model's minage.
#
# `ref_lines` is a function(models, sp_sel) returning per-species reference
# layers; `suffix` names the saved file when it differs from the series. Both
# are handed to plot_timeseries(), which adds the layers before it saves, so the
# reference lines reach the written PNG and the CSV still has the `file` it
# needs. They exist so every timeseries plotter goes through this factory: F and
# SSB depletion both draw reference points, and writing either by hand outside
# it is how plot_f() came to take a different argument list from its siblings
# and to resolve `species` its own way.
.ts_wrapper <- function(output, zero_y = FALSE, ref_lines = NULL,
                        suffix = NULL) {
  force(output); force(zero_y); force(ref_lines); force(suffix)
  function(Rceattle,
           file = NULL,
           model_names = NULL,
           line_col = NULL,
           species = NULL,
           spnames = NULL,
           add_ci = FALSE,
           lwd = 3,
           save = FALSE,
           right_adj = 0,
           legend.pos = "topright",
           width = 7,
           height = 6.5,
           minyr = NULL,
           maxyr = NULL,
           incl_proj = FALSE,
           mod_cex = 1,
           lty = rep(1, length(Rceattle)),
           alpha = 0.4,
           mod_avg = rep(FALSE, length(Rceattle)),
           mse = FALSE,
           OM = TRUE,
           reference = NULL,
           ylab = NULL) {

    plot_timeseries(Rceattle,
                    output = output,
                    ylab = ylab,
                    zero_y = zero_y,
                    ref_lines = ref_lines,
                    suffix = suffix,
                    file = file,
                    model_names = model_names,
                    line_col = line_col,
                    species = species,
                    spnames = spnames,
                    add_ci = add_ci,
                    lwd = lwd,
                    save = save,
                    right_adj = right_adj,
                    legend.pos = legend.pos,
                    width = width,
                    height = height,
                    minyr = minyr,
                    maxyr = maxyr,
                    incl_proj = incl_proj,
                    mod_cex = mod_cex,
                    lty = lty,
                    alpha = alpha,
                    mod_avg = mod_avg,
                    mse = mse,
                    OM = OM,
                    reference = reference)
  }
}

#' Plot biomass
#'
#' @description Plots the mean minage+ biomass (million mt) and 95\% CI trends as estimated from Rceattle.
#'
#' @inheritParams plot_timeseries
#'
#' @export
#'
#' @return Returns and saves a figure with the biomass trajectory.
plot_biomass <- .ts_wrapper("biomass", zero_y = TRUE)


#' Plot recruitment
#'
#' @description Plots the mean minage recruitment (millions) and 95\% CI trends as estimated from Rceattle.
#'
#' @inheritParams plot_timeseries
#'
#' @export
#'
#' @return Returns and saves a figure with the recruitment trajectory.
plot_recruitment <- .ts_wrapper("R", zero_y = TRUE)


#' Plot spawning stock biomass (SSB)
#'
#' @description Plots the mean SSB (million mt) and 95\% CI trends as estimated from Rceattle.
#'
#' @inheritParams plot_timeseries
#'
#' @export
#'
#' @return Returns and saves a figure with the SSB trajectory.
plot_ssb <- .ts_wrapper("ssb", zero_y = TRUE)


#' Plot exploitable biomass
#'
#' @description Plots the mean exploitable biomass (million mt) and 95\% CI trends as estimated from Rceattle.
#'
#' @inheritParams plot_timeseries
#'
#' @export
#'
#' @return Returns and saves a figure with the exploitable biomass trajectory.
plot_exploitable_biomass <- .ts_wrapper("exploitable_biomass", zero_y = TRUE)

#' Plot SSB depletion
#'
#' @description Plots the mean SSB depletion and 95\% CI trends as estimated from Rceattle.
#'   Depletion reference lines for Ptarget and Plimit are drawn in blue and red respectively.
#'
#' @inheritParams plot_timeseries
#'
#' @export
#'
#' @return Returns and saves a figure with the SSB depletion trajectory.
plot_depletionSSB <- .ts_wrapper("ssb_depletion", zero_y = TRUE,
                                 ref_lines = .depletion_reference_lines)

#' Plot SSB depletion (deprecated name)
#'
#' @description Deprecated alias for \code{\link{plot_depletionSSB}}. Please use
#'   \code{plot_depletionSSB()} instead.
#'
#' @inheritParams plot_timeseries
#' @export
plot_ssb_depletion <- plot_depletionSSB

#' Plot biomass depletion
#'
#' @description Plots the mean biomass depletion and 95\% CI trends as estimated from Rceattle.
#'
#' @inheritParams plot_timeseries
#'
#' @export
#'
#' @return Returns and saves a figure with the biomass depletion trajectory.
plot_depletion <- .ts_wrapper("biomass_depletion", zero_y = TRUE)


#' Plot selectivity
#'
#' @description Selectivity over the hindcast, one panel per fleet. Each fleet
#'   is drawn on the dimension it was estimated on: age for an age-based fleet,
#'   length bin for a length-based one (`Selectivity_dimension`). A model mixing
#'   the two gets one figure per dimension.
#'
#' @details
#' # What the colours mean
#'
#' With one model, colour is the year, so time-varying selectivity reads as a
#' fan and a time-invariant fleet collapses to a single line. With several
#' models, colour separates the models and the year fan moves to transparency,
#' faintest at the earliest year drawn. The transparency scale spans the years
#' shown across all models, so the same year is the same shade everywhere and a
#' short retrospective peel stops short of solid. `colour_by` overrides the
#' choice either way, and `alpha` sets the faintest end.
#'
#' Line type separates the sexes, and `lty` supplies its values. Panels are
#' fleets, so `spnames` does not label anything here, it only lets `species`
#' select by name.
#'
#' # Confidence intervals
#'
#' `add_ci = TRUE` draws `exp(log(sel) +/- 1.96 * sd)`, so the band is positive
#' and right-skewed. It needs a fit run with `fit_control(selectivity_se = TRUE)`.
#'
#' Only estimated, age-based fleets hold one. A length-based fleet is drawn on
#' `sel_at_length`, which is what was fitted and holds no error; a
#' `Selectivity = "Fixed"` fleet estimates nothing; and no fleet gets a band
#' below its first selected bin (`Bin_first_selected`, a 1-based bin ordinal,
#' not an absolute age), where selectivity is 0 by construction. A
#' fleet mirroring another's `Selectivity_index` has no errors of its own and
#' borrows its lead's, but only where the two curves agree, `data_check()` only
#' warns when a shared group differs in a shaping column such as `Sel_norm_bin`.
#'
#' Every year drawn gets its own band, so a time-varying fleet is dense. Pair it
#' with `minyr` / `maxyr` on a single year, or use it on a time-invariant fleet.
#'
#' The interval is **relative to the normalization anchor**, so quote the anchor
#' with it. Selectivity is identified only up to a scalar, jointly with F for a
#' fishery and q for a survey, so `Sel_norm_bin` fixes the level rather than
#' measuring it: the bin it pins has a standard error of exactly 0 and the band
#' widens away from it. On a penalized fit (`random_sel = FALSE`) the band is also
#' conditional on the smoothing, and it is marginal rather than joint. See
#' `vignette("model-options-and-functionality")` for both, and for the
#' `estimateMode` that determines whether the fit reports an error at all.
#'
#' @inheritParams rceattle-plot-args
#' @param colour_by What colour separates: `"year"` (a fan), `"model"`, or
#'   `"auto"` (the default) for year with a single model and model with several.
#' @param alpha Transparency of the faintest year when the fan is drawn by
#'   transparency, i.e. under `colour_by = "model"`.
#'
#' @return A `ggplot`, or for a model mixing age- and length-based fleets a
#'   named list of them, one per dimension.
#' @export
plot_selectivity <-
  function(Rceattle,
           file = NULL,
           model_names = NULL,
           line_col = NULL,
           width = 7,
           height = 6.5,
           species = NULL,
           lwd = 3,
           spnames = NULL,
           lty = 1,
           minyr = NULL,
           maxyr = NULL,
           alpha = 0.25,
           colour_by = c("auto", "year", "model"),
           add_ci = FALSE) {

    colour_by <- match.arg(colour_by)
    .rce_check_alpha(alpha)
    models <- .as_plot_models(Rceattle)

    # One model without the standard errors would draw a band on some panels and
    # not others, which reads as "this fleet is certain" rather than "this fit
    # was not asked for them", so the request is declined for the whole figure.
    if (isTRUE(add_ci)) {
      se_all  <- lapply(models, .rce_sel_se)
      lacking <- vapply(se_all, is.null, logical(1))
      if (any(lacking)) {
        add_ci <- .rce_no_ci(
          add_ci, "selectivity",
          paste0("model(s) ", paste(which(lacking), collapse = ", "),
                 " report no selectivity standard error. Refit with ",
                 "fit_control(selectivity_se = TRUE, getsd = TRUE)"))
      }
    }
    # No usable error, for one of two reasons that need different advice. All
    # zero: the sdreport held no selectivity free, which is what an HCR
    # projection returns. None finite: the Hessian was not invertible.
    if (isTRUE(add_ci)) {
      zero <- vapply(se_all, function(s)
        any(is.finite(s$sd)) && all(s$sd == 0, na.rm = TRUE), logical(1))
      nonf <- vapply(se_all, function(s) !any(is.finite(s$sd)), logical(1))
      if (any(zero)) {
        add_ci <- .rce_no_ci(
          add_ci, "selectivity",
          paste0("model(s) ", paste(which(zero), collapse = ", "),
                 " report a standard error of 0 for every age, so the sdreport ",
                 "came from the HCR projection, where the selectivity ",
                 "parameters are fixed. Refit with estimateMode = \"Hindcast\", ",
                 "or fit_control(projection_uncertainty = TRUE)"))
      }
      if (any(nonf)) {
        add_ci <- .rce_no_ci(
          add_ci, "selectivity",
          paste0("model(s) ", paste(which(nonf), collapse = ", "),
                 " report no finite selectivity standard error, so the Hessian ",
                 "was not invertible. See fit$convergence -- the ",
                 "`sdreport_failed` and `hessian_conditioning` checks"))
      }
    }
    model_names_use <- .model_labels(models, model_names)
    sp_sel  <- .resolve_species(models, species, spnames)
    keep_sp <- sp_sel$index

    # Every model contributes: reading only the first silently dropped the rest
    # of an overlay. Each is read with its own fleet_control and dimensions,
    # since a comparison run may differ in both.
    se_list <- if (isTRUE(add_ci)) se_all else vector("list", length(models))

    df_list <- list()
    for (k in seq_along(models)) {
      dl     <- models[[k]]$data_list
      fc     <- dl$fleet_control
      nages  <- dl$nages
      minage <- dl$minage
      nsex   <- dl$nsex
      nlen   <- dl$nlengths
      hindyears <- dl$styr:dl$endyr
      # "Age" / "Length" per fleet; the column defaults to Age and may be blank
      # on an older workbook.
      dimn <- as.character(fc$Selectivity_dimension)
      dimn[is.na(dimn) | !nzchar(dimn)] <- "Age"
      # Empirical selectivity is read straight into sel_at_age -- the model
      # skips sel_type 0 when it fills sel_at_length (selectivity.hpp, "ESTIMATED
      # SELECTIVITY") -- so a Fixed fleet is age-based whatever the column says,
      # and reading the length array would draw it as identically zero. Scripts
      # set Selectivity_dimension across every fleet at once, so this happens.
      is_fixed <- as.character(fc$Selectivity) %in% c("Fixed", "0")
      dimn[is_fixed] <- "Age"

      for (i in seq_len(nrow(fc))) {
        flt <- fc$Fleet_code[i]
        sp  <- fc$Species[i]
        if (is.na(sp) || !(sp %in% keep_sp)) next
        is_len <- identical(dimn[i], "Length")
        # A length-based fleet estimates selectivity on length bins; sel_at_age
        # is that curve pushed through the growth matrix, not what was fitted.
        sel  <- if (is_len) models[[k]]$quantities$sel_at_length
                else        models[[k]]$quantities$sel_at_age
        if (is.null(sel)) {
          warning("Model ", model_names_use[k], " reports no ",
                  if (is_len) "sel_at_length" else "sel_at_age",
                  "; skipping fleet ", fc$Fleet_name[i], ".", call. = FALSE)
          next
        }
        nbin <- if (is_len) nlen[sp] else nages[sp]
        if (is.na(nbin) || nbin < 1L) next
        # Length bins are 1-based ordinals (see the column schema); ages are
        # offset by the species' minage.
        bins <- if (is_len) seq_len(nbin) else seq_len(nbin) - 1L + minage[sp]

        # A mirrored fleet shares the lead's parameter block and reports no rows
        # of its own, so it borrows the lead's band.
        sel_lead <- .rce_sel_lead(dl, i)

        for (sex in seq_len(nsex[sp])) {
          sex_lab <- if (nsex[sp] == 1) "Combined" else c("Female", "Male")[sex]
          for (yr in seq_along(hindyears)) {
            est <- as.numeric(sel[flt, sex, seq_len(nbin), yr])
            # The standard error is on sel_at_age, so a length-based fleet gets
            # NA and geom_ribbon drops it. NA rather than the estimate, so an
            # absent band never reads as a zero-width one.
            lo <- up <- rep(NA_real_, length(est))
            if (isTRUE(add_ci) && !is_len) {
              ci <- .rce_sel_ci(se_list[[k]], flt, sex, bins, hindyears[yr],
                                lead = sel_lead, curve = est)
              if (!is.null(ci)) { lo <- ci$lower; up <- ci$upper }
            }
            df_list[[length(df_list) + 1L]] <- data.frame(
              Model = model_names_use[k],
              Fleet = as.character(fc$Fleet_name[i]),
              Dimension = if (is_len) "Length" else "Age",
              Sex   = sex_lab,
              Bin   = bins,
              Year  = hindyears[yr],
              Selectivity = est,
              Lower = lo,
              Upper = up,
              stringsAsFactors = FALSE)
          }
        }
      }
    }
    if (length(df_list) == 0L) {
      stop("No fleets to plot: `species` selected none of the fleets in this ",
           "model.", call. = FALSE)
    }

    plot_df <- do.call(rbind, df_list)
    plot_df <- .rce_year_filter(plot_df, minyr, maxyr)
    # droplevels: a species filter can leave a model with no rows, and an
    # unused level would still make the figure look like an overlay.
    plot_df$Model <- droplevels(
      factor(plot_df$Model, levels = unique(model_names_use)))
    plot_df$Fleet <- factor(plot_df$Fleet, levels = unique(plot_df$Fleet))

    if (identical(colour_by, "auto")) {
      colour_by <- if (nlevels(plot_df$Model) > 1L) "model" else "year"
    }
    if (identical(colour_by, "model")) {
      .rce_check_line_col(line_col, nlevels(plot_df$Model), "models")
    } else if (nlevels(plot_df$Model) > 1L) {
      warning("`colour_by = \"year\"` with ", nlevels(plot_df$Model),
              " models draws them superimposed with nothing to tell them ",
              "apart. Use colour_by = \"model\", or vary `lwd`.",
              call. = FALSE)
    }

    # One figure per dimension: ages and length bins do not share an axis.
    dims <- unique(plot_df$Dimension)
    out <- lapply(dims, function(d) {
      dd <- droplevels(plot_df[plot_df$Dimension == d, , drop = FALSE])
      .plot_selectivity_one(dd, dimension = d, colour_by = colour_by,
                            line_col = line_col, lwd = lwd, lty = lty,
                            alpha = alpha, add_ci = add_ci,
                            file = if (is.null(file)) NULL else
                              paste0(file, if (length(dims) > 1L)
                                paste0("_", tolower(d)) else ""),
                            width = width, height = height)
    })
    names(out) <- tolower(dims)
    if (length(out) == 1L) out[[1]] else out
  }


# Render one selectivity panel grid, for a single Selectivity_dimension.
# Split out because a model mixing age- and length-based fleets needs one of
# these per dimension and they differ only in the x axis label.
.plot_selectivity_one <- function(dd, dimension, colour_by, line_col, lwd, lty,
                                  alpha, add_ci, file, width, height) {
  xlab <- if (identical(dimension, "Length")) "Length bin" else "Age"
  # A length figure, or a dimension whose fleets are all Fixed, has no band to
  # draw even when the request stands for the age figure beside it.
  draw_ci <- isTRUE(add_ci) && any(is.finite(dd$Lower))

  if (identical(colour_by, "model")) {
    # Colour is the model, so the year fan moves to transparency. The scale is
    # trained on the years shown across all models, so a short retrospective
    # peel stops short of solid rather than being rescaled to its own end --
    # the same year reads as the same shade in every model.
    p <- ggplot2::ggplot(
      dd,
      ggplot2::aes(x = .data$Bin, y = .data$Selectivity,
                   group = interaction(.data$Year, .data$Sex, .data$Model,
                                       .data$Fleet),
                   colour = .data$Model, alpha = .data$Year,
                   linetype = .data$Sex))
    lp <- .rce_line_params(lwd = lwd, lty = lty,
                           lwd_by = "Model",
                           lwd_n = nlevels(dd$Model),
                           lty_by = "Sex",
                           lty_n = length(unique(dd$Sex)),
                           lty_in_aes = TRUE)
    if (draw_ci) {
      # Under the ribbon, so a band never hides the curve it belongs to. The
      # fill follows the model and inherits the year alpha, so a band is shaded
      # like the line above it; show.legend keeps a second key off the figure.
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$Lower, ymax = .data$Upper, fill = .data$Model),
        colour = NA, alpha = 0.15, show.legend = FALSE, na.rm = TRUE)
      p <- .rceattle_scale(p, aesthetics = "fill", line_col = line_col)
    }
    p <- .rce_add_line(p, lp)
    p <- .rceattle_scale(p, aesthetics = "colour", line_col = line_col) +
      ggplot2::scale_alpha_continuous(range = c(alpha, 1), guide = "none")
  } else {
    p <- ggplot2::ggplot(
      dd,
      ggplot2::aes(x = .data$Bin, y = .data$Selectivity,
                   group = interaction(.data$Year, .data$Sex, .data$Model,
                                       .data$Fleet),
                   colour = .data$Year, linetype = .data$Sex))
    lp <- .rce_line_params(lwd = lwd, lty = lty,
                           lwd_by = "Model",
                           lwd_n = nlevels(dd$Model),
                           lty_by = "Sex",
                           lty_n = length(unique(dd$Sex)),
                           lty_in_aes = TRUE)
    if (draw_ci) {
      # Colour is a continuous year here, so the band takes a single neutral
      # grey rather than a second year ramp -- two ramps on one panel read as
      # two variables.
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$Lower, ymax = .data$Upper),
        fill = "grey40", colour = NA, alpha = 0.15, na.rm = TRUE)
    }
    p <- .rce_add_line(p, lp)
    # Colour is the year, a magnitude, so `line_col` gives the ramp anchors.
    p <- .rceattle_scale(p, discrete = FALSE, aesthetics = "colour",
                         line_col = line_col)
  }

  p <- p +
    ggplot2::facet_wrap(~ Fleet) +
    ggplot2::labs(x = xlab, y = "Selectivity") +
    .rceattle_theme()
  if (length(unique(dd$Sex)) < 2L) p <- p + ggplot2::guides(linetype = "none")
  if (nlevels(dd$Model) < 2L && identical(colour_by, "model")) {
    p <- p + ggplot2::guides(colour = "none")
  }

  .save_ggplot(p, file = file, suffix = "selectivity",
               width = width, height = height)
}




#' Plot functional form
#'
#' @description Function to plot the functional form estimated or specified by \code{Rceattle}
#'
#' @param params Parameter list object from \code{\link{build_params}} or \code{Rceattle}
#' @param pred Predator index
#' @param pred_age Predator age
#' @param prey Prey index
#' @param msmMode Multispecies mode integer specifying functional form
#' @export
plot_form <- function( params = NULL, pred = 1, pred_age = 1, prey = 1, msmMode = 3){
  # The Kinzey & Punt (2009) functional-response parameters are commented out
  # in the TMB template, build_params(), build_map() and build_bounds(), and
  # data_check() blocks msmMode 3:9. Body retained below so this can be revived
  # alongside the C++.
  stop(paste0(
    "plot_form() plots the Kinzey & Punt (2009) functional responses ",
    "(msmMode 3-9), which are not available in this version of Rceattle: ",
    "the functional-response parameters are not estimated, and data_check() ",
    "blocks these modes. Use msmMode = 1 (Holsman et al. 2015 MSVPA) or ",
    "2 (Holling Type III MSVPA)."), call. = FALSE)

  # ---- retained for future development ----
  #
  #   # Get indices
  #   rsp = pred
  #   ksp = prey
  #
  #   # Get parameter values
  #   H_1 <- exp(params$logH_1)
  #   H_1a <- exp(params$logH_1a)
  #   H_1b <- exp(params$logH_1b)
  #
  #   H_2 <- exp(params$logH_2)
  #   H_3 <- exp(params$logH_3)
  #
  #   H_4 <- params$H_4
  #
  #   # Set up ratios
  #   Pred_r <- seq(from = 0.001, to = 5, length.out = 100) # Pred biomass relative to equilibrium
  #   Prey_r <- seq(from = 0.001, to = 5, length.out = 100) # Prey biomass relative to equilibrium
  #
  #   # Calculate functional form
  #   Term = H_1[rsp, ksp] * (1 + H_1a[rsp] * H_1b[rsp] / (pred_age + H_1b[rsp]))
  #   response <- matrix(NA, ncol = length(Prey_r), nrow = length(Pred_r))
  #   rownames(response) <- Pred_r
  #   colnames(response) <- Prey_r
  #
  #   for(i in 1:length(Pred_r)){
  #     for(j in 1:length(Prey_r)){
  #       response[i, j] <- Prey_r[j] * switch (
  #         as.character(msmMode),
  #         "2" = { # Holling Type I (linear)
  #           Term},
  #         "3" = { #  Holling Type II
  #           Term * (1 + H_2[rsp, ksp] ) / ( 1 + H_2[rsp, ksp] * Prey_r[j] )
  #         },
  #         "4" = { #  Holling Type III
  #           Term * (1 + H_2[rsp, ksp]) * ((Prey_r[j] ) ^ H_4[rsp, ksp]) / (1 + H_2[rsp, ksp] * ((Prey_r[j] ) ^ H_4[rsp, ksp])  )},
  #         "5" = { #  predator interference
  #           Term * (1 + H_2[rsp, ksp] ) / ( 1 + H_2[rsp, ksp] * Prey_r[j] + H_3[rsp, ksp] * (Pred_r[i] - 1) )
  #         },
  #         "6" = { # predator preemption
  #           Term * (1 + H_2[rsp, ksp] ) / ( (1 + H_2[rsp, ksp] * Prey_r[j]) * (1 + H_3[rsp, ksp] * (Pred_r[i] - 1)) )
  #         },
  #         "7" = { # Hassell-Varley
  #           Term * (2 + H_2[rsp, ksp] ) / (1 + H_2[rsp, ksp] * Prey_r[j] + ((Prey_r[j] ) ^ H_4[rsp, ksp]))
  #         },
  #         "8" = { #  Ecosim
  #           Term / (1 + H_3[rsp, ksp] * (Pred_r[i] - 1 ))},
  #         {
  #           stop("msmMode not implemented: ", msmMode)
  #         }
  #       )
  #     }
  #   }
  #
  #   df <- reshape2::melt(response)
  #   colnames(df) <- c("Pred_r", "Prey_r", "Response")
  #   df$Pred_r <- as.numeric(as.character(df$Pred_r))
  #   df$Prey_r <- as.numeric(as.character(df$Prey_r))
  #
  #   if (msmMode %in% c(2, 3, 4, 7)) {
  #     # functional response depends on the prey ratio only
  #     d <- stats::aggregate(Response ~ Prey_r, data = df, FUN = mean)
  #     p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$Prey_r, y = .data$Response)) +
  #       ggplot2::geom_line(linewidth = 1) +
  #       ggplot2::labs(x = "Prey ratio", y = "Functional response")
  #   } else if (msmMode == 8) {
  #     d <- stats::aggregate(Response ~ Pred_r, data = df, FUN = mean)
  #     p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$Pred_r, y = .data$Response)) +
  #       ggplot2::geom_line(linewidth = 1) +
  #       ggplot2::labs(x = "Predator ratio", y = "Functional response")
  #   } else if (msmMode %in% c(5, 6)) {
  #     p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Prey_r, y = .data$Pred_r,
  #                                           fill = .data$Response)) +
  #       ggplot2::geom_raster() +
  #       ggplot2::scale_fill_viridis_c("Functional\nresponse") +
  #       ggplot2::labs(x = "Prey ratio", y = "Predator ratio")
  #   } else {
  #     stop("msmMode not implemented: ", msmMode)
  #   }
  #   p + .rceattle_theme()
}





#' Plot M1 or M2 at age
#'
#' @description Mortality-at-age over the hindcast for one model: predation
#'   mortality (`M2`, the default) or residual natural mortality (`M1`). One
#'   component at a time, never their sum, [plot_m_at_age()] draws total M
#'   (M1 + M2) as a time series.
#'
#' @param file name of a file to identified the files exported by the
#'   function.
#' @param Rceattle A single [fit_mod()] object. A list of more than one is an
#'   error; call it per model.
#' @param incl_proj Include the projection years (TRUE/FALSE)
#' @param type `"lines"` (default) draws M-at-age with one line per year;
#'   `"heatmap"` (or `0`) draws the same series as an age-by-year tile plot.
#'   Any other value gives the lines.
#' @param width Plot width when saved "inches"
#' @param height Plot height when saved "inches"
#' @param title Additional title to add. Will also add species names if not NULL
#' @param species Species to plot. Plots all if null.
#' @param log TRUE/FALSE plot the series on a log scale
#' @param minyr First year to plot
#' @param maxage Oldest age to draw, on the species' own age scale rather than
#'   as a count of age bins. `NULL` draws every age.
#' @param M2 Draw predation mortality M2 (TRUE, the default) or residual
#'   natural mortality M1 (FALSE). Neither is the sum of the two.
#' @param zlim Ignored. Base-graphics leftover: the fill range of the tile
#'   plot. Use `p + ggplot2::scale_fill_viridis_c(limits = ...)`.
#' @param theta Ignored. Base-graphics leftover: viewing angle of a `persp`
#'   surface, which is no longer drawn.
#' @param title_cex Ignored. Base-graphics leftover: title font size. Use
#'   `p + ggplot2::theme(plot.title = ggplot2::element_text(size = ...))`.
#'
#' @return A `ggplot`.
#' @export
plot_mortality <-
  function(Rceattle,
           file = NULL,
           incl_proj = FALSE,
           zlim = NULL,
           type = "lines",
           width = 8,
           height = 5.5,
           title = NULL,
           log = FALSE,
           minyr = NULL,
           theta = 155,
           species = NULL,
           maxage = NULL,
           title_cex = 10,
           M2 = TRUE) {

    Rceattle <- .as_plot_models(Rceattle)
    if (length(Rceattle) > 1) stop("Can only plot one model")
    dl <- Rceattle[[1]]$data_list

    if (is.null(minyr)) minyr <- dl$styr
    endyr <- if (incl_proj) dl$projyr else dl$endyr
    years <- minyr:endyr
    nyrs  <- length(years)

    nspp        <- dl$nspp
    spnames     <- dl$spnames
    estdynamics <- dl$estDynamics
    nages       <- dl$nages
    minage      <- dl$minage
    # `maxage` is an age, so the number of bins it leaves is how many of the
    # species' ages are at or below it -- pmin(nages, maxage) treats it as a bin
    # count, which truncates in the wrong place for any minage but 1.
    if (!is.null(maxage)) {
      nages <- pmin(nages, pmax(0L, as.integer(maxage) - as.integer(minage) + 1L))
      if (all(nages < 1L)) {
        stop("`maxage` = ", maxage, " is below the youngest age of every ",
             "species.", call. = FALSE)
      }
    }
    nsex   <- dl$nsex
    if (is.null(species)) species <- seq_len(nspp)
    spp <- intersect(species, seq_len(nspp))

    # M-at-age over the selected years: M2 (M2 = TRUE) or M1 (M2 = FALSE).
    # One component or the other -- plot_m_at_age() is what sums them.
    yr_idx <- seq_len(nyrs) + (minyr - dl$styr)
    Msrc <- if (M2) Rceattle[[1]]$quantities$M2_at_age
            else    Rceattle[[1]]$quantities$M1_at_age

    df_list <- list()
    for (sp in spp) {
      if (estdynamics[sp] != 0) next
      # A species whose youngest age is already above `maxage` has no bins left.
      if (nages[sp] < 1L) next
      ages <- seq_len(nages[sp]) - 1 + minage[sp]
      for (sex in seq_len(nsex[sp])) {
        sex_lab <- if (nsex[sp] == 1) "Combined" else c("Female", "Male")[sex]
        for (yi in seq_along(years)) {
          val <- as.numeric(Msrc[sp, sex, seq_len(nages[sp]), yr_idx[yi]])
          df_list[[length(df_list) + 1L]] <- data.frame(
            Species = spnames[sp], Sex = sex_lab,
            Age = ages, Year = years[yi],
            M = if (log) log(val) else val,
            stringsAsFactors = FALSE)
        }
      }
    }
    plot_df <- do.call(rbind, df_list)
    plot_df$Species <- factor(plot_df$Species, levels = spnames[spp])
    # Msrc is M2_at_age when M2 = TRUE, else M1_at_age; label matches the source.
    ylab <- if (M2) "M2" else "M1"
    if (log) ylab <- paste0("log(", ylab, ")")

    if (identical(type, "heatmap") || identical(type, 0)) {
      # Age x year heatmap, faceted by species
      p <- ggplot2::ggplot(plot_df,
                           ggplot2::aes(x = .data$Year, y = .data$Age,
                                        fill = .data$M)) +
        ggplot2::geom_tile() +
        .facet_species(plot_df) +
        ggplot2::scale_fill_viridis_c(ylab) +
        ggplot2::labs(x = "Year", y = "Age") +
        .rceattle_theme()
    } else {
      # Year-coloured M-at-age curves (standard), faceted by species
      p <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(x = .data$Age, y = .data$M,
                     group = interaction(.data$Year, .data$Sex),
                     colour = .data$Year, linetype = .data$Sex)) +
        ggplot2::geom_line() +
        .facet_species(plot_df, scales = "free_y") +
        ggplot2::scale_colour_viridis_c("Year") +
        ggplot2::labs(x = "Age", y = ylab) +
        .rceattle_theme()
      if (length(unique(plot_df$Sex)) < 2L) {
        p <- p + ggplot2::guides(linetype = "none")
      }
    }
    if (!is.null(title)) p <- p + ggplot2::ggtitle(title)

    .save_ggplot(p, file = file, suffix = "mortality_at_age",
                 width = width, height = height)
  }







#' Plot maturity
#'
#' @description Function that plots the maturity of each species
#'
#' @param file name of a file to identified the files exported by the
#'   function.
#' @param Rceattle Single or list of Rceattle model objects exported from \code{Rceattle}
#' @param model_names Names of models to be used in legend
#' @param line_col Colors of models to be used for line color
#' @param species Species names for legend
#' @param lwd Line width as specified by user
#' @param width Figure width in inches
#' @param height Figure height in inches
#'
#' @export
plot_maturity <-
  function(Rceattle,
           file = NULL,
           model_names = NULL,
           line_col = NULL,
           species = NULL,
           width = 4,
           height = 5.5,
           lwd = 3) {

    Rceattle <- .as_plot_models(Rceattle)
    model_names_use <- .model_labels(Rceattle, model_names)
    nspp  <- Rceattle[[1]]$data_list$nspp
    nages <- Rceattle[[1]]$data_list$nages
    if (is.null(species)) species <- Rceattle[[1]]$data_list$spnames

    # Tidy maturity-at-age (data_list$maturity, species column dropped)
    df_list <- list()
    for (k in seq_along(Rceattle)) {
      mat <- Rceattle[[k]]$data_list$maturity[, -1, drop = FALSE]
      for (sp in seq_len(nspp)) {
        df_list[[length(df_list) + 1L]] <- data.frame(
          Model    = model_names_use[k],
          Species  = species[sp],
          Age      = seq_len(nages[sp]),
          Maturity = as.numeric(mat[sp, seq_len(nages[sp])]),
          stringsAsFactors = FALSE)
      }
    }
    plot_df <- do.call(rbind, df_list)
    plot_df$Model   <- factor(plot_df$Model, levels = unique(model_names_use))
    plot_df$Species <- factor(plot_df$Species, levels = species)

    p <- ggplot2::ggplot(plot_df,
                         ggplot2::aes(x = .data$Age, y = .data$Maturity,
                                      colour = .data$Model)) +
      ggplot2::geom_line(linewidth = 1) +
      .facet_species(plot_df, ncol = 1) +
      ggplot2::coord_cartesian(ylim = c(0, 1.1)) +
      ggplot2::labs(x = "Age", y = "Maturity")
    p <- .rceattle_scale(p + .rceattle_theme(), aesthetics = "colour")
    if (nlevels(plot_df$Model) < 2L) p <- p + ggplot2::guides(colour = "none")

    .save_ggplot(p, file = file, suffix = "maturity",
                 width = width, height = height)
  }


#' Plot biomass eaten by predation
#'
#' @description Total biomass of each species consumed by its predators, summed
#'   over sex and age. For an MSE object the projection is summarized across
#'   simulations with a mean and a 95% interval.
#'
#' @inheritParams rceattle-plot-args
#' @param save Ignored. Only the time-series plotters write their series to CSV.
#' @param mse Is an MSE object from \code{\link{load_mse}} or \code{\link{run_mse}}
#' @param OM if mse == TRUE, use the OM (TRUE) or EM (FALSE) for plotting?
#' @param mod_avg TRUE/FALSE
#'
#' @export
#'
plot_b_eaten <- function(Rceattle,
                         file = NULL,
                         model_names = NULL,
                         line_col = NULL,
                         species = NULL,
                         spnames = NULL,
                         add_ci = FALSE,
                         lwd = 3,
                         save = FALSE,
                         right_adj = 0,
                         width = 7,
                         height = 6.5,
                         minyr = NULL,
                         incl_proj = FALSE,
                         mod_cex = 1,
                         alpha = 0.4,
                         mod_avg = rep(FALSE, length(Rceattle)),
                         mse = FALSE,
                         OM = TRUE,
                         maxyr = NULL,
                         lty = 1,
                         incl_mean = FALSE,
                         top_adj = 0.15) {

  .rce_check_alpha(alpha)
  models <- .as_plot_models(Rceattle, mse = mse, OM = OM)
  if (mse) incl_proj <- TRUE
  model_names_use <- .model_labels(models, model_names)
  sp_sel  <- .resolve_species(models, species, spnames)
  species <- sp_sel$index
  spnames <- sp_sel$spnames
  # B_eaten_as_prey is REPORTed but not ADREPORTed, so there is no standard
  # error to draw from.
  add_ci <- .rce_no_ci(add_ci, "biomass eaten",
                       "the model reports it without standard errors",
                       warn = !mse)

  # Total biomass eaten as prey per species/year: sum of B_eaten_as_prey over
  # sex and age. The model reports it in mt; the /1e6 puts it in million mt,
  # the display unit every other biomass axis in the package uses -- including
  # plot_b_eaten_prop(), which breaks this same quantity down by predator.
  df_list <- list()
  for (k in seq_along(models)) {
    dl  <- models[[k]]$data_list
    yrs <- dl$styr:(if (incl_proj) dl$projyr else dl$endyr)
    be  <- models[[k]]$quantities$B_eaten_as_prey
    tot <- apply(be[, , , seq_along(yrs), drop = FALSE], c(1, 4), sum)
    for (sp in species) {
      df_list[[length(df_list) + 1L]] <- data.frame(
        Model = model_names_use[k], Species = spnames[sp],
        Year = yrs, value = as.numeric(tot[sp, ]) / 1e6,
        stringsAsFactors = FALSE)
    }
  }
  plot_df <- do.call(rbind, df_list)
  plot_df <- .rce_year_filter(plot_df, minyr, maxyr)
  plot_df$Model   <- factor(plot_df$Model, levels = unique(model_names_use))
  plot_df$Species <- factor(plot_df$Species, levels = spnames[species])

  if (mse) {
    # Simulations are summarized into one band, so there are no model series to
    # colour or to tell apart by line type.
    if (!is.null(line_col)) {
      warning("`line_col` is not used with `mse = TRUE`: the simulations are ",
              "summarized into a single band, not drawn as separate series.",
              call. = FALSE)
    }
    agg <- stats::aggregate(value ~ Species + Year, plot_df,
      FUN = function(v) c(m = mean(v),
                          l95 = stats::quantile(v, 0.025, names = FALSE),
                          u95 = stats::quantile(v, 0.975, names = FALSE)))
    mdf <- data.frame(Species = agg$Species, Year = agg$Year, agg$value)
    # One summarized series, so no keying variable -- but lwd / lty are still
    # validated and applied through the same helper as everywhere else.
    lp <- .rce_line_params(lwd = lwd, lty = lty)
    p <- ggplot2::ggplot(mdf, ggplot2::aes(x = .data$Year)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$l95, ymax = .data$u95),
                           alpha = alpha, fill = "grey40") +
      do.call(ggplot2::geom_line,
              c(list(mapping = ggplot2::aes(y = .data$m)), lp$args))
    p <- p +
      .facet_species(mdf, scales = "free_y") +
      .rce_mean_line(mdf, incl_mean, by = "Species", value = "m",
                     hind_endyr = max(.rce_model_endyr(models, model_names_use),
                                      na.rm = TRUE)) +
      ggplot2::labs(x = "Year", y = "Biomass eaten as prey (million mt)") +
      .rceattle_theme()
    return(.save_ggplot(p, file = file, suffix = "biomass_eaten",
                        width = width, height = height))
  }

  .rce_check_line_col(line_col, length(unique(plot_df$Model)), "models")
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data$Year, y = .data$value, colour = .data$Model))
  p <- .rce_add_line(p, .rce_line_params(
    lwd = lwd, lty = lty,
    lwd_by = "Model", lwd_n = length(unique(plot_df$Model)),
    lty_by = "Model", lty_n = length(unique(plot_df$Model))))
  p <- p +
    .facet_species(plot_df, scales = "free_y") +
    .rce_proj_divider(models, incl_proj, minyr, maxyr) +
    .rce_mean_line(plot_df, incl_mean, by = c("Species", "Model"),
                   hind_endyr = .rce_model_endyr(models, model_names_use),
                   colour_by = "Model") +
    ggplot2::labs(x = "Year", y = "Biomass eaten as prey (million mt)")
  p <- .rceattle_scale(p + .rceattle_theme(), aesthetics = "colour",
                       line_col = line_col)
  if (nlevels(plot_df$Model) < 2L) p <- p + ggplot2::guides(colour = "none")
  .save_ggplot(p, file = file, suffix = "biomass_eaten",
               width = width, height = height)
}



#' Plot biomass consumed of each prey species by predator
#'
#' @description Function that plots the biomass consumed trends as estimated from Rceattle. Returns and saves a figure with the biomass eaten trajectory.
#'
#' @details Colour separates predators; line type separates models. Panels are
#'   prey species.
#'
#' @inheritParams rceattle-plot-args
#' @param mohns Ignored. Formerly annotated the figure with Mohn's rho from
#'   [retrospective()]; add it with `ggplot2::labs(subtitle = ...)` instead.
#'
#' @export
#'
plot_b_eaten_prop <-
  function(Rceattle,
           file = NULL,
           model_names = NULL,
           line_col = NULL,
           spnames = NULL,
           species = NULL,
           lwd = 3,
           right_adj = 0,
           top_adj = 0.15,
           minyr = NULL,
           mohns = NULL,
           width = 7,
           height = 6.5,
           incl_proj = FALSE,
           incl_mean = FALSE,
           add_ci = FALSE,
           mod_cex = 1,
           maxyr = NULL,
           lty = 1) {

    models <- .as_plot_models(Rceattle)
    model_names_use <- .model_labels(models, model_names)
    sp_sel  <- .resolve_species(models, species, spnames)
    species <- sp_sel$index
    spnames <- sp_sel$spnames
    nspp    <- models[[1]]$data_list$nspp
    max_sex <- max(models[[1]]$data_list$nsex)
    add_ci <- .rce_no_ci(add_ci, "biomass eaten by predator",
                         "the model reports it without standard errors")

    # Biomass of prey ksp eaten by predator rsp (summed over sex/age), in
    # million mt. B_eaten indexing kept verbatim from the original.
    df_list <- list()
    for (k in seq_along(models)) {
      dl  <- models[[k]]$data_list
      yrs <- dl$styr:(if (incl_proj) dl$projyr else dl$endyr)
      Be  <- models[[k]]$quantities$B_eaten
      for (rsp in 1:nspp) {
        for (ksp in species) {
          val <- vapply(seq_along(yrs), function(yr)
            sum(Be[c(rsp, (rsp + nspp) * (max_sex - 1)),
                   c(ksp, (ksp + nspp) * (max_sex - 1)), , , yr, drop = FALSE]),
            numeric(1))
          df_list[[length(df_list) + 1L]] <- data.frame(
            Model = model_names_use[k], Prey = spnames[ksp],
            Predator = spnames[rsp], Year = yrs, value = val / 1e6,
            stringsAsFactors = FALSE)
        }
      }
    }
    plot_df <- do.call(rbind, df_list)
    plot_df <- .rce_year_filter(plot_df, minyr, maxyr)
    plot_df$Predator <- factor(plot_df$Predator, levels = spnames)
    plot_df$Prey     <- factor(plot_df$Prey, levels = spnames[species])
    .rce_check_line_col(line_col, length(unique(plot_df$Predator)), "predators")

    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = .data$Year, y = .data$value,
                   colour = .data$Predator, linetype = .data$Model))
    p <- .rce_add_line(p, .rce_line_params(
      lwd = lwd, lty = lty,
      lwd_by = "Model", lwd_n = length(unique(plot_df$Model)),
      lty_by = "Model", lty_n = length(unique(plot_df$Model)),
      lty_in_aes = TRUE))
    p <- p +
      ggplot2::facet_wrap(~ Prey, scales = "free_y") +
      .rce_proj_divider(models, incl_proj, minyr, maxyr) +
      .rce_mean_line(plot_df, incl_mean, by = c("Prey", "Predator", "Model"),
                     hind_endyr = .rce_model_endyr(models, model_names_use),
                     colour_by = "Predator") +
        ggplot2::labs(x = "Year", y = "Biomass eaten (million mt)")
    p <- .rceattle_scale(p + .rceattle_theme(), aesthetics = "colour",
                         line_col = line_col)
    if (length(unique(plot_df$Model)) < 2L) p <- p + ggplot2::guides(linetype = "none")
    .save_ggplot(p, file = file, suffix = "biomass_eaten_by_predator",
                 width = width, height = height)
  }




#' Plot natural mortality by age
#'
#' @description Function that plots the natural mortality at age (M1 + M2) as estimated from Rceattle. Returns and saves a figure with the M-at-age trajectory.
#'
#' @details Colour separates the models; line type separates the sexes. A
#'   sex-combined model has one sex, so a varying `lty` has nothing to key on
#'   and warns.
#'
#' @inheritParams rceattle-plot-args
#' @param age Age to plot M at, on the species' own age scale (so `age = 1` is
#'   age 1, not the first age bin, and means nothing to a species whose `minage`
#'   is 2). A species that has no such age is dropped, with a warning.
#'
#' @export
#'
plot_m_at_age <-
  function(Rceattle,
           file = NULL,
           age = 1,
           model_names = NULL,
           line_col = NULL,
           spnames = NULL,
           species = NULL,
           lwd = 3,
           lty = 1,
           right_adj = 0,
           minyr = NULL,
           width = 7,
           height = 6.5,
           incl_proj = FALSE,
           incl_mean = FALSE,
           add_ci = FALSE,
           maxyr = NULL,
           top_adj = 0.15) {

    Rceattle <- .as_plot_models(Rceattle)
    model_names_use <- .model_labels(Rceattle, model_names)
    sp_sel  <- .resolve_species(Rceattle, species, spnames)
    species <- sp_sel$index
    spnames <- sp_sel$spnames
    nsex <- Rceattle[[1]]$data_list$nsex
    # `age` is an age, not a bin index: the arrays run 1 .. nages while the ages
    # run minage .. minage + nages - 1, so it is resolved per species.
    age_sel <- .rce_age_index(age, species, Rceattle[[1]]$data_list$minage,
                              Rceattle[[1]]$data_list$nages, spnames)
    species <- age_sel$species
    # M_at_age is REPORTed but its ADREPORT is commented out in the model,
    # so the fit carries no standard errors for it.
    add_ci <- .rce_no_ci(add_ci, "M at age",
                         "the model reports it without standard errors")

    # Total natural mortality (M1 + M2) at a single age, as a time series.
    df_list <- list()
    for (k in seq_along(Rceattle)) {
      dl  <- Rceattle[[k]]$data_list
      yrs <- dl$styr:(if (incl_proj) dl$projyr else dl$endyr)
      Ma  <- Rceattle[[k]]$quantities$M_at_age
      for (j in seq_along(species)) {
        sp  <- species[j]
        bin <- age_sel$index[j]
        for (sex in seq_len(nsex[sp])) {
          sex_lab <- if (nsex[sp] == 1) "Combined" else c("Female", "Male")[sex]
          df_list[[length(df_list) + 1L]] <- data.frame(
            Model   = model_names_use[k],
            Species = spnames[sp],
            Sex     = sex_lab,
            Year    = yrs,
            M       = as.numeric(Ma[sp, sex, bin, seq_along(yrs)]),
            stringsAsFactors = FALSE)
        }
      }
    }
    plot_df <- do.call(rbind, df_list)
    plot_df <- .rce_year_filter(plot_df, minyr, maxyr)
    plot_df$Model   <- factor(plot_df$Model, levels = unique(model_names_use))
    plot_df$Species <- factor(plot_df$Species, levels = spnames[species])
    .rce_check_line_col(line_col, length(unique(plot_df$Model)), "models")

    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = .data$Year, y = .data$M,
                   colour = .data$Model, linetype = .data$Sex))
    # Line type separates the sexes here, so `lty` supplies the sex values.
    p <- .rce_add_line(p, .rce_line_params(
      lwd = lwd, lty = lty,
      lwd_by = "Model", lwd_n = length(unique(plot_df$Model)),
      lty_by = "Sex",   lty_n = length(unique(plot_df$Sex)),
      lty_in_aes = TRUE))
    p <- p +
      .facet_species(plot_df, scales = "free_y") +
      .rce_proj_divider(Rceattle, incl_proj, minyr, maxyr) +
      .rce_mean_line(plot_df, incl_mean, by = c("Species", "Model", "Sex"),
                     value = "M",
                     hind_endyr = .rce_model_endyr(Rceattle, model_names_use),
                     colour_by = "Model") +
        ggplot2::labs(x = "Year", y = paste0("M at age ", age))
    p <- .rceattle_scale(p + .rceattle_theme(), aesthetics = "colour",
                         line_col = line_col)
    if (length(unique(plot_df$Sex)) < 2L) p <- p + ggplot2::guides(linetype = "none")
    if (nlevels(plot_df$Model) < 2L)      p <- p + ggplot2::guides(colour = "none")

    .save_ggplot(p, file = file, suffix = paste0("m_at_age", age),
                 width = width, height = height)
  }


#' Plot predation mortality by age and predator
#'
#' @description Share of the predation mortality (M2) on each prey age that is
#'   attributable to each predator species. The shares sum to 1 across
#'   predators for every prey age and year; a year with no predation on that
#'   age leaves them undefined and draws nothing.
#'
#' @details Colour separates predators; line type separates models. Panels are
#'   prey species (and sex, where the prey is sexed).
#'
#' @inheritParams rceattle-plot-args
#' @param age Prey age to plot the M2 proportions at, on the prey species' own
#'   age scale rather than as an age-bin index. A prey species that has no such
#'   age is dropped, with a warning.
#'
#' @export
#'
plot_m2_at_age_prop <-
  function(Rceattle,
           file = NULL,
           age = 1,
           model_names = NULL,
           line_col = NULL,
           spnames = NULL,
           species = NULL,
           lwd = 3,
           right_adj = 0,
           top_adj = 0.15,
           minyr = NULL,
           width = 7,
           height = 6.5,
           incl_proj = FALSE,
           incl_mean = FALSE,
           add_ci = FALSE,
           maxyr = NULL,
           lty = 1) {

    Rceattle <- .as_plot_models(Rceattle)
    model_names_use <- .model_labels(Rceattle, model_names)
    sp_sel  <- .resolve_species(Rceattle, species, spnames)
    species <- sp_sel$index
    spnames <- sp_sel$spnames
    years <- lapply(Rceattle, function(x) x$data_list$styr:x$data_list$endyr)
    if (incl_proj) years <- lapply(Rceattle, function(x) x$data_list$styr:x$data_list$projyr)
    nyrs_vec <- sapply(years, length)
    nyrs <- max(nyrs_vec)
    nsex <- Rceattle[[1]]$data_list$nsex
    nspp <- Rceattle[[1]]$data_list$nspp
    # `age` is a prey age, not a bin index; resolve it per prey species.
    age_sel <- .rce_age_index(age, species, Rceattle[[1]]$data_list$minage,
                              Rceattle[[1]]$data_list$nages, spnames)
    species <- age_sel$species
    # A ratio of REPORT-only arrays: no standard errors to draw from.
    add_ci <- .rce_no_ci(add_ci, "the M2 proportions",
                         "they are a ratio of series the model reports without standard errors")

    # Share of the predation mortality on each prey age attributable to each
    # predator species.
    #
    # M2_prop holds each predator's CONTRIBUTION to M2, not a share of it:
    # summing it over predators reproduces M2_at_age exactly. Plotting that sum
    # gave a "proportion" that reached 1500 on BS2017MS. Dividing by the total
    # over predators gives the share the axis has always claimed, which sums to
    # 1 across predators for every prey age and year.
    m2_at_age_prop <- array(NA, dim = c(nspp, nspp, 2, nyrs, length(Rceattle)))
    for (i in 1:length(Rceattle)) {
      for (j in seq_along(species)) {
        ksp <- species[j]
        # `age` is a PREY age here, so it resolves against the prey species'
        # own minage -- the 4th dimension of M2_prop is the prey age bin.
        k_bin <- age_sel$index[j]
        for (k_sex in 1:nsex[ksp]) {
          for (yr in 1:nyrs_vec[i]) {
            by_pred <- vapply(1:nspp, function(rsp)
              sum(Rceattle[[i]]$quantities$M2_prop[
                c(rsp, (rsp + nspp) * (max(nsex) - 1)),
                ksp + (nspp * (k_sex - 1)), , k_bin, yr]),
              numeric(1))
            total <- sum(by_pred)
            # No predation on this prey age in this year leaves the shares
            # undefined rather than 0/0.
            m2_at_age_prop[, ksp, k_sex, yr, i] <-
              if (total > 0) by_pred / total else NA_real_
          }
        }
      }
    }

    # Tidy: one row per (model, prey species/sex panel, predator, year)
    df_list <- list()
    for (k in seq_along(Rceattle)) {
      yrs <- years[[k]]
      for (ksp in species) {            # prey
        for (sex in seq_len(nsex[ksp])) {
          panel <- if (nsex[ksp] == 1) spnames[ksp]
                   else paste(spnames[ksp], c("Female", "Male")[sex])
          for (rsp in seq_len(nspp)) {  # predator
            df_list[[length(df_list) + 1L]] <- data.frame(
              Model     = model_names_use[k],
              Prey      = panel,
              Predator  = spnames[rsp],
              Year      = yrs,
              Proportion = as.numeric(m2_at_age_prop[rsp, ksp, sex, seq_along(yrs), k]),
              stringsAsFactors = FALSE)
          }
        }
      }
    }
    plot_df <- do.call(rbind, df_list)
    plot_df <- .rce_year_filter(plot_df, minyr, maxyr)
    plot_df$Predator <- factor(plot_df$Predator, levels = spnames)
    # Panels follow the order `species` asked for, as in the other plotters.
    plot_df$Prey     <- factor(plot_df$Prey, levels = unique(plot_df$Prey))
    .rce_check_line_col(line_col, length(unique(plot_df$Predator)), "predators")

    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = .data$Year, y = .data$Proportion,
                   colour = .data$Predator, linetype = .data$Model))
    p <- .rce_add_line(p, .rce_line_params(
      lwd = lwd, lty = lty,
      lwd_by = "Model", lwd_n = length(unique(plot_df$Model)),
      lty_by = "Model", lty_n = length(unique(plot_df$Model)),
      lty_in_aes = TRUE))
    p <- p +
      ggplot2::facet_wrap(~ Prey, scales = "free_y") +
      .rce_proj_divider(Rceattle, incl_proj, minyr, maxyr) +
      .rce_mean_line(plot_df, incl_mean, by = c("Prey", "Predator", "Model"),
                     value = "Proportion",
                     hind_endyr = .rce_model_endyr(Rceattle, model_names_use),
                     colour_by = "Predator") +
        ggplot2::labs(x = "Year",
                      y = paste0("Share of M2 at age ", age, " by predator"))
    p <- .rceattle_scale(p + .rceattle_theme(), aesthetics = "colour",
                         line_col = line_col)
    if (length(unique(plot_df$Model)) < 2L) p <- p + ggplot2::guides(linetype = "none")

    .save_ggplot(p, file = file, suffix = paste0("m2_at_age_prop", age),
                 width = width, height = height)
  }


#' plot F
#'
#' @description Fishing mortality over time, one panel per species, with the
#'   Ftarget (blue) and Flimit (red) reference points drawn per panel.
#'
#' @inheritParams plot_timeseries
#' @export
#'
#' @return Returns and saves a figure with the fishing mortality trajectory.
plot_f <- .ts_wrapper("F_spp", ref_lines = .f_reference_lines,
                      suffix = "f_trajectory")


#' Plot ration
#'
#' @description Population-level consumption for ages `minage`+: the individual
#'   annual ration (kg/yr) multiplied by average numbers-at-age and summed over
#'   age, in million mt. This is how the model forms total consumption
#'   (`avgN_at_age * ration`, `predation.hpp`), so it is the consumption that
#'   generates the predation mortality in [plot_m2_at_age_prop()] and the
#'   biomass in [plot_b_eaten()], plus the other-food term.
#'
#' @details Colour separates the models; line type separates the sexes.
#'
#' @inheritParams rceattle-plot-args
#' @param minage Youngest age to sum consumption over, so the figure is
#'   "age `minage`+". An age on the species' own scale, not an age-bin index; a
#'   species with no age that old is dropped, with a warning.
#'
#' @export
#'
plot_ration <-
  function(Rceattle,
           file = NULL,
           minage = 1,
           model_names = NULL,
           line_col = NULL,
           spnames = NULL,
           species = NULL,
           lwd = 3,
           lty = 1,
           right_adj = 0,
           minyr = NULL,
           width = 7,
           height = 6.5,
           incl_proj = FALSE,
           incl_mean = FALSE,
           add_ci = FALSE,
           maxyr = NULL,
           top_adj = 0.15) {

    models <- .as_plot_models(Rceattle)
    model_names_use <- .model_labels(models, model_names)
    sp_sel  <- .resolve_species(models, species, spnames)
    species <- sp_sel$index
    spnames <- sp_sel$spnames
    nsex  <- models[[1]]$data_list$nsex
    nages <- models[[1]]$data_list$nages
    # `minage` names an age ("age 3+"), not a bin, so the bins summed over are
    # resolved against each species' own age vector rather than assumed to start
    # at the same index for all of them.
    age_sel <- .rce_age_plus_index(minage, species,
                                   models[[1]]$data_list$minage, nages, spnames)
    species <- age_sel$species
    # Consumption is the product of two REPORT-only arrays, so its standard
    # error would need a covariance the fit does not carry.
    add_ci <- .rce_no_ci(add_ci, "consumption", paste(
      "it is a product of two series the model reports without standard",
      "errors"))

    # Population-level consumption: individual ration x average numbers-at-age,
    # summed over age minage+. consumption_at_age is the annual ration of ONE
    # fish (kg/yr) and numbers-at-age are in thousands, so the product is mt.
    #
    # avgN_at_age, not N_at_age: the model multiplies the ration by the
    # year's average numbers (predation.hpp, `avgN_at_age(rsp, ...) * pred_rat`),
    # and ration_hat is B_eaten_as_pred / avgN_at_age. Start-of-year numbers
    # would overstate consumption by 1 / ((1 - exp(-Z)) / Z) and would not
    # reconcile with plot_b_eaten().
    #
    # Multiplying the ration by biomass instead weights the age-sum by
    # weight-at-age and is not a quantity in any unit.
    df_list <- list()
    for (k in seq_along(models)) {
      dl  <- models[[k]]$data_list
      yrs <- dl$styr:(if (incl_proj) dl$projyr else dl$endyr)
      Caa <- models[[k]]$quantities$consumption_at_age
      Naa <- models[[k]]$quantities$avgN_at_age
      for (j in seq_along(species)) {
        sp   <- species[j]
        bins <- age_sel$index[[j]]
        for (sex in seq_len(nsex[sp])) {
          sex_lab <- if (nsex[sp] == 1) "Combined" else c("Female", "Male")[sex]
          val <- vapply(seq_along(yrs), function(yr)
            sum(Caa[sp, sex, bins, yr] * Naa[sp, sex, bins, yr]),
            numeric(1)) / 1e6
          df_list[[length(df_list) + 1L]] <- data.frame(
            Model = model_names_use[k], Species = spnames[sp], Sex = sex_lab,
            Year = yrs, value = val, stringsAsFactors = FALSE)
        }
      }
    }
    plot_df <- do.call(rbind, df_list)
    plot_df <- .rce_year_filter(plot_df, minyr, maxyr)
    plot_df$Model   <- factor(plot_df$Model, levels = unique(model_names_use))
    plot_df$Species <- factor(plot_df$Species, levels = spnames[species])
    .rce_check_line_col(line_col, length(unique(plot_df$Model)), "models")

    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = .data$Year, y = .data$value,
                   colour = .data$Model, linetype = .data$Sex))
    # Line type separates the sexes here, so `lty` supplies the sex values.
    p <- .rce_add_line(p, .rce_line_params(
      lwd = lwd, lty = lty,
      lwd_by = "Model", lwd_n = length(unique(plot_df$Model)),
      lty_by = "Sex",   lty_n = length(unique(plot_df$Sex)),
      lty_in_aes = TRUE))
    p <- p +
      .facet_species(plot_df, scales = "free_y") +
      .rce_proj_divider(models, incl_proj, minyr, maxyr) +
      .rce_mean_line(plot_df, incl_mean, by = c("Species", "Model", "Sex"),
                     hind_endyr = .rce_model_endyr(models, model_names_use),
                     colour_by = "Model") +
        ggplot2::labs(x = "Year",
                    y = paste0("Consumption (million mt), age ", minage, "+"))
    p <- .rceattle_scale(p + .rceattle_theme(), aesthetics = "colour",
                         line_col = line_col)
    if (length(unique(plot_df$Sex)) < 2L) p <- p + ggplot2::guides(linetype = "none")
    if (nlevels(plot_df$Model) < 2L)      p <- p + ggplot2::guides(colour = "none")
    .save_ggplot(p, file = file, suffix = paste0("ration", minage, "plus"),
                 width = width, height = height)
  }

