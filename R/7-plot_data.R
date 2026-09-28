#' Timeline of data used in the model likelihoods
#'
#' Plots a timeline of presence/absence (and optionally relative quantity) of
#' the data sources contributing to an Rceattle model's likelihood, by year and
#' fleet. Modeled after Stock Synthesis's `r4ss::SSplotData()`.
#'
#' Years with `Year > 0` contribute to the likelihood; rows with `Year < 0`
#' (kept for fitting comparison only) are shown when `ghost = TRUE`. Diet data
#' uses `Year == 0` for the fixed reference stomach composition; those entries
#' are plotted at the model `endyr` and labeled per predator species.
#'
#' @param Rceattle Either a single Rceattle model object exported from
#'   \code{Rceattle::fit_mod()}, or a raw `data_list` (e.g. one of the
#'   bundled datasets such as `BS2017MS`).
#' @param file Path/prefix used to save the figure. If `NULL` (default) the
#'   plot is drawn to the active device only. When supplied, the figure is
#'   written to `paste0(file, "_data_plot.png")` (and `_data_plot2.png` for
#'   the bubble subplot).
#' @param subplots Integer vector controlling which subplots are produced:
#'   \itemize{
#'     \item 1, equal-size points showing presence/absence by year/fleet
#'     \item 2, points scaled to relative quantity / precision within each
#'       data type (catch tonnage, 1/SE for indices, sample size for comps)
#'   }
#' @param datatypes Either `"all"` or a subset of
#'   `c("catch", "index", "agecomp", "lencomp", "caal", "diet")`.
#' @param fleets Either `"all"` or an integer vector of `Fleet_code`s to
#'   include. Diet "fleets" are encoded after the regular fleets (one per
#'   predator species).
#' @param species Either `"all"` or an integer vector of species indices to
#'   include.
#' @param ghost Logical. If `TRUE`, also include rows excluded from the
#'   likelihood (`Year < 0`). Defaults to `FALSE`.
#' @param fleetcol Either `"default"` or a vector of colors (one per fleet
#'   shown after subsetting). Diet entries get an additional color appended.
#' @param width,height Figure dimensions in inches.
#' @param res Resolution (dpi) for the saved PNG.
#' @param ptsize Pointsize passed to `png()`.
#' @param margins `par("mar")` for the plot. Increase the right margin if
#'   long fleet names are clipped.
#' @param cex Character expansion for points showing isolated years.
#' @param lwd Line width for runs of consecutive years (subplot 1).
#' @param maxsize Max bubble radius (in plot units) for subplot 2.
#' @param alphasize Bubble fill transparency (0-1).
#' @param mainTitle Logical; if `TRUE` add a default title.
#' @param cex.main Title character expansion.
#'
#' @return Invisibly, a list with `typetable`, the long data frame underlying
#'   the plot (year, fleet, data type, relative size).
#'
#' @export
plot_data <- function(Rceattle,
                      file = NULL,
                      subplots = 1:2,
                      datatypes = "all",
                      fleets = "all",
                      species = "all",
                      ghost = FALSE,
                      fleetcol = "default",
                      width = 7,
                      height = 5,
                      res = 200,
                      ptsize = 10,
                      margins = c(5.1, 2.1, 2.1, 10.1),
                      cex = 2,
                      lwd = 12,
                      maxsize = 1,
                      alphasize = 1,
                      mainTitle = FALSE,
                      cex.main = 1) {

  # Accept either a fitted Rceattle object or a raw data_list
  if (inherits(Rceattle, "Rceattle")) {
    data_list <- Rceattle$data_list
  } else if (is.list(Rceattle) &&
             all(c("fleet_control", "nspp", "endyr") %in% names(Rceattle))) {
    data_list <- Rceattle
  } else if (is.list(Rceattle) && is.list(Rceattle$data_list)) {
    data_list <- Rceattle$data_list
  } else {
    stop("`Rceattle` must be either an Rceattle model object or a data_list ",
         "(a list containing `fleet_control`, `nspp`, `endyr`, ...).")
  }
  fleet_control <- data_list$fleet_control
  nspp <- data_list$nspp
  endyr <- data_list$endyr

  if (length(species) == 1 && species == "all") {
    species <- seq_len(nspp)
  }

  # Species names for labelling diet entries
  spnames <- data_list$spnames
  if (is.null(spnames) || length(spnames) < nspp) {
    spnames <- paste0("Sp", seq_len(nspp))
  }

  # Catalogue of data types: name, label, source df
  typetable_def <- data.frame(
    type     = c("catch", "index", "agecomp", "lencomp", "caal", "diet"),
    label    = c("Catch", "Abundance index", "Age composition",
                 "Length composition", "Conditional age-at-length",
                 "Diet (stomach contents)"),
    stringsAsFactors = FALSE
  )

  if (length(datatypes) == 1 && datatypes == "all") {
    datatypes <- typetable_def$type
  }
  typetable_def <- typetable_def[typetable_def$type %in% datatypes, , drop = FALSE]

  # Helper: aggregate a fleet/year data frame, return year + size vectors
  agg_year <- function(df, value_col, fun = sum) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    out <- stats::aggregate(df[[value_col]],
                            by = list(Yr = abs(df$Year)),
                            FUN = fun, na.rm = TRUE)
    names(out) <- c("yr", "size")
    out[is.finite(out$size) & out$size > 0, , drop = FALSE]
  }

  # Build the long table of (year, fleet, type, size) entries
  rows <- list()
  push <- function(yr, size, fleet, type) {
    if (length(yr) == 0) return(invisible())
    rows[[length(rows) + 1L]] <<- data.frame(
      yr = yr, size = size, fleet = fleet, type = type,
      stringsAsFactors = FALSE
    )
  }

  # Optional Year-sign filter (drop excluded rows unless ghost = TRUE)
  filt_used <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(df)
    if (ghost) df else df[df$Year > 0, , drop = FALSE]
  }

  for (this_type in typetable_def$type) {

    if (this_type == "catch") {
      df <- filt_used(data_list$catch_data)
      if (!is.null(df) && nrow(df) > 0) {
        df <- df[df$Species %in% species, , drop = FALSE]
        for (flt in unique(df$Fleet_code)) {
          agg <- agg_year(df[df$Fleet_code == flt, ], "Catch", sum)
          if (!is.null(agg)) push(agg$yr, agg$size, flt, this_type)
        }
      }
    }

    if (this_type == "index") {
      df <- filt_used(data_list$index_data)
      if (!is.null(df) && nrow(df) > 0) {
        df <- df[df$Species %in% species, , drop = FALSE]
        for (flt in unique(df$Fleet_code)) {
          # precision = 1 / mean Log_sd within year
          agg <- agg_year(df[df$Fleet_code == flt, ], "Log_sd", mean)
          if (!is.null(agg)) push(agg$yr, 1 / agg$size, flt, this_type)
        }
      }
    }

    if (this_type %in% c("agecomp", "lencomp")) {
      df <- filt_used(data_list$comp_data)
      if (!is.null(df) && nrow(df) > 0) {
        flag <- if (this_type == "agecomp") 0 else 1
        df <- df[df$Age0_Length1 == flag & df$Species %in% species, , drop = FALSE]
        for (flt in unique(df$Fleet_code)) {
          agg <- agg_year(df[df$Fleet_code == flt, ], "Sample_size", sum)
          if (!is.null(agg)) push(agg$yr, agg$size, flt, this_type)
        }
      }
    }

    if (this_type == "caal") {
      df <- filt_used(data_list$caal_data)
      if (!is.null(df) && nrow(df) > 0 && "Sample_size" %in% names(df)) {
        df <- df[df$Species %in% species, , drop = FALSE]
        for (flt in unique(df$Fleet_code)) {
          agg <- agg_year(df[df$Fleet_code == flt, ], "Sample_size", sum)
          if (!is.null(agg)) push(agg$yr, agg$size, flt, this_type)
        }
      }
    }

    if (this_type == "diet") {
      df <- data_list$diet_data
      # Diet enters the likelihood only in multispecies mode (msmMode > 0) and
      # only for predators with non-zero composition weight. Fall back to
      # "show if present" when those switches are absent.
      if (!is.null(df) && nrow(df) > 0) {
        msm <- if (is.null(data_list$msmMode)) 1L else data_list$msmMode
        wts <- data_list$Diet_comp_weights
        if (is.null(wts)) wts <- rep(1, nspp)
        active_pred <- if (msm > 0) which(wts > 0) else integer(0)
        df <- df[df$Pred %in% intersect(active_pred, species), , drop = FALSE]
        if (nrow(df) > 0) {
          for (pred in unique(df$Pred)) {
            sub <- df[df$Pred == pred, , drop = FALSE]
            # Year == 0 means fixed reference: pin to endyr for display
            disp_yrs <- ifelse(sub$Year == 0, endyr, abs(sub$Year))
            agg <- stats::aggregate(sub$Sample_size,
                                    by = list(yr = disp_yrs),
                                    FUN = sum, na.rm = TRUE)
            names(agg) <- c("yr", "size")
            # Encode predator as a pseudo-fleet after the regular fleets
            push(agg$yr, agg$size, max(fleet_control$Fleet_code) + pred, this_type)
          }
        }
      }
    }
  }

  if (length(rows) == 0) {
    message("No data of the requested types/fleets/species found in this model.")
    return(invisible(NULL))
  }
  typetable <- do.call(rbind, rows)

  # Subset by user-requested fleets
  if (!(length(fleets) == 1 && fleets == "all")) {
    typetable <- typetable[typetable$fleet %in% fleets, , drop = FALSE]
    if (nrow(typetable) == 0) {
      message("No data remaining after fleet subset.")
      return(invisible(NULL))
    }
  }

  # Order data types by the catalogue (so layout is consistent)
  typetable$itype <- match(typetable$type, typetable_def$type)
  typetable <- typetable[order(typetable$itype, typetable$fleet, typetable$yr), ]

  # Build display fleet list (regular fleets + diet pseudo-fleets per species)
  fleets_used <- sort(unique(typetable$fleet))
  reg_max <- max(fleet_control$Fleet_code)
  fleet_label <- function(f) {
    if (f <= reg_max) {
      lab <- fleet_control$Fleet_name[match(f, fleet_control$Fleet_code)]
      if (is.na(lab)) lab <- paste0("Fleet ", f)
      lab
    } else {
      paste0("Diet: ", spnames[f - reg_max])
    }
  }

  # Colors: one per displayed fleet
  nflt <- length(fleets_used)
  if (length(fleetcol) == 1 && fleetcol == "default") {
    fleetcol <- if (nflt == 1) "grey40"
                else if (nflt == 2) c("blue", "red")
                else if (nflt == 3) c("blue", "red", "green3")
                else rich.colors.short(nflt + 1)[-1]
  } else if (length(fleetcol) < nflt) {
    fleetcol <- rep(fleetcol, length.out = nflt)
  }

  # Build display labels and relative point size (within each data type, max
  # -> 1), then render a faceted ggplot: year on x, fleet on y, one facet row
  # per data type, point area scaled to the relative quantity/precision.
  typetable$Fleet_label <- vapply(typetable$fleet, fleet_label, character(1))
  typetable$Type_label  <- factor(
    typetable_def$label[match(typetable$type, typetable_def$type)],
    levels = typetable_def$label)

  typetable <- do.call(rbind, lapply(split(typetable, typetable$type),
    function(d) {
      mx <- max(d$size, na.rm = TRUE)
      d$rel <- if (is.finite(mx) && mx > 0) d$size / mx else 0
      d
    }))
  fl_order <- unique(typetable$Fleet_label[order(typetable$itype,
                                                 typetable$fleet)])
  typetable$Fleet_label <- factor(typetable$Fleet_label, levels = rev(fl_order))

  p <- ggplot2::ggplot(
    typetable,
    ggplot2::aes(x = .data$yr, y = .data$Fleet_label,
                 colour = .data$Type_label, size = .data$rel)) +
    ggplot2::geom_vline(xintercept = endyr + 0.5, linetype = 2,
                        colour = "grey50") +
    ggplot2::geom_point(alpha = alphasize) +
    ggplot2::facet_grid(rows = ggplot2::vars(.data$Type_label),
                        scales = "free_y", space = "free_y") +
    ggplot2::scale_size_area(max_size = 6, guide = "none") +
    ggplot2::labs(x = "Year", y = NULL)
  p <- .rceattle_scale(p + .rceattle_theme() +
                         ggplot2::theme(legend.position = "none"),
                       aesthetics = "colour")

  .save_ggplot(p, file = file, suffix = "data_plot",
               width = width, height = height)
}
