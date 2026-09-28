# =============================================================================
# Shared helpers for the Rceattle ggplot2 plotting functions
# -----------------------------------------------------------------------------
# The exported plot_*() functions follow one standardized contract:
#   1. coerce the `Rceattle` argument to a list of fits      (.as_plot_models)
#   2. assemble a tidy data frame of exactly what is plotted (per function)
#   3. render it with the shared theme + colourblind palette  (.rceattle_theme,
#      .rceattle_scale)
#   4. optionally write the figure to `file`                  (.save_ggplot)
#   5. RETURN the ggplot object (so plots are composable and the plotted data
#      is testable as `p$data`).
#
# Keeping the data assembly separate from rendering is deliberate: tests assert
# the assembled data equals the model's source quantities, so a ggplot plotter
# is verified to draw the *correct numbers*, not merely to render.
# =============================================================================

# Largest |residual| the bubble figures draw at full size. Every residual-bubble
# panel is pinned to [0, this] so a well-fitting fleet and a badly-fitting one
# cannot be made to look alike by a free scale, and so the OSA and Pearson
# panels drawn from the same fit stay comparable. Beyond it the bubble stops
# growing (.rce_truncate_resid) -- a standard normal residual past 6 is already
# far outside anything a size scale can usefully distinguish.
.RCE_BUBBLE_MAX <- 6

#' Coerce the `Rceattle` plotting argument to a list of fits
#'
#' Accepts a single `Rceattle` fit, a list of fits (multi-model overlay), or an
#' MSE object, and returns a plain list of `Rceattle` fits.
#'
#' @param Rceattle A fit, list of fits, or (when `mse = TRUE`) MSE object.
#' @param mse,OM When `mse = TRUE`, pull the operating model (`OM = TRUE`) or
#'   the terminal estimation model from each MSE element.
#' @return A list of `Rceattle` fits.
#' @noRd
.as_plot_models <- function(Rceattle, mse = FALSE, OM = TRUE) {
  if (mse) {
    if (OM) {
      return(lapply(Rceattle, function(x) x$OM))
    }
    return(lapply(Rceattle, function(x) x$EM[[length(x$EM)]]))
  }
  if (inherits(Rceattle, "Rceattle")) {
    return(list(Rceattle))
  }
  Rceattle
}

#' Model display names for a list of fits
#'
#' Returns `model_names` if supplied, otherwise `"Model 1"`, `"Model 2"`, ... so
#' the colour/legend mapping always has labels.
#'
#' `model_names` is often built as a `list()`, the package's own vignettes do
#', so it is flattened to character here. Left as a list it becomes a
#' one-element list per model and the plot frame fails to bind.
#'
#' Too few names would be recycled, drawing two models as one series under one
#' legend key, so that is a warning rather than a silent merge.
#'
#' @param model_list List of fits.
#' @param model_names Labels, or `NULL` for positional ones.
#' @keywords internal
.model_labels <- function(model_list, model_names = NULL) {
  n <- length(model_list)
  if (!is.null(model_names)) {
    model_names <- as.character(unlist(model_names))
    if (length(model_names) != n) {
      warning("`model_names` has ", length(model_names), " name",
              if (length(model_names) == 1L) "" else "s", " for ", n,
              " models; recycling. Models sharing a name are drawn as one ",
              "series.", call. = FALSE)
    }
    return(rep(model_names, length.out = n))
  }
  if (n == 1L) {
    return("Model")
  }
  paste("Model", seq_len(n))
}

#' Standard Rceattle ggplot theme
#'
#' Clean, legible, and colorblind-friendly: `theme_classic` with bold facet
#' strips, a light panel border, and the legend along the bottom.
#' @param base_size Base font size.
#' @keywords internal
.rceattle_theme <- function(base_size = 12) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text       = ggplot2::element_text(face = "bold", size = base_size),
      panel.border     = ggplot2::element_rect(colour = "black", fill = NA,
                                               linewidth = 0.5),
      legend.position  = "bottom",
      legend.title     = ggplot2::element_blank()
    )
}

#' Okabe-Ito colorblind-safe qualitative palette
#'
#' Eight fixed hues in a fixed order, reordered so the strongest-contrast,
#' best-separated hues come first (2-4 model overlays are the common case) and
#' the low-contrast yellow comes last. Unlike a sequential ramp (viridis), this
#' encodes *identity* (which model), not magnitude: the colours hold no implied
#' ordering and every hue holds a legible contrast against a white panel.
#'
#' Reference: Okabe & Ito (2008), "Color Universal Design".
#' @keywords internal
.okabe_ito <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7",
                "#56B4E9", "#E69F00", "#000000", "#F0E442")

#' Discrete palette function: Okabe-Ito for up to 8 series, else interpolate
#'
#' Mirrors the arity contract of `scale_*_viridis_d` (accepts any `n`) so the
#' scale never errors on a many-model overlay; beyond 8 series it interpolates
#' the Okabe-Ito anchors (rare enough that CVD-optimality is not guaranteed).
#' @keywords internal
.okabe_ito_pal <- function(n) {
  if (n <= length(.okabe_ito)) return(.okabe_ito[seq_len(n)])
  grDevices::colorRampPalette(.okabe_ito)(n)
}

#' Normalize a user colour specification to something ggplot2 accepts
#'
#' `line_col` is a base-graphics argument and users pass base-graphics colours:
#' `line_col = 1` and `line_col = c(1, 2, 4)` are as common as `"black"` or
#' `"#0072B2"`. Base wraps indices past the end of the palette; ggplot2 does not,
#' so they are resolved here.
#'
#' Invalid colours fail here, where the argument still has a name, rather than
#' silently, an `NA` colour draws nothing, so a bad `line_col` would otherwise
#' give a blank panel and no message.
#'
#' @param x Colours as a character vector, or integers indexing
#'   [grDevices::palette()], or `NULL`.
#' @return A character vector of colours, or `NULL` for `NULL`/empty input.
#' @keywords internal
#' @noRd
.as_colour <- function(x) {
  if (is.null(x) || length(x) == 0L) return(NULL)
  if (is.logical(x) && all(is.na(x))) {
    stop("`line_col` is NA. Pass colours, palette indices, or NULL for the ",
         "default palette.", call. = FALSE)
  }
  if (is.numeric(x)) {
    if (anyNA(x) || any(x < 0)) {
      stop("`line_col` must be non-negative palette indices; got ",
           paste(utils::head(x, 5), collapse = ", "), ".", call. = FALSE)
    }
    pal <- grDevices::palette()
    # Base wraps the palette by (i - 1) %% n and reads 0 as the background,
    # i.e. transparent.
    return(ifelse(x == 0, NA_character_,
                  pal[(as.integer(x) - 1L) %% length(pal) + 1L]))
  }
  if (!is.character(x) && !is.factor(x)) {
    stop("`line_col` must be colour names, hex codes, or palette indices; got ",
         class(x)[1], ".", call. = FALSE)
  }
  x <- as.character(x)
  # Reject unrenderable colours now. grDevices::col2rgb() is the same parser
  # ggplot2 eventually uses, so anything it accepts will draw.
  bad <- vapply(x, function(one) {
    if (is.na(one)) return(FALSE)   # NA is a legitimate "draw nothing"
    inherits(try(grDevices::col2rgb(one), silent = TRUE), "try-error")
  }, logical(1))
  if (any(bad)) {
    stop("`line_col` is not a valid colour: ",
         paste(unique(x[bad]), collapse = ", "), ".", call. = FALSE)
  }
  x
}


#' Warn when a colour vector is too short for what it is about to colour
#'
#' `line_col` recycles, so too few colours merge series: `line_col = 1` on
#' `plot_b_eaten_prop()`, whose colour separates predators, draws every predator
#' in black. The message names the variable the colours landed on, since a
#' caller passing one colour usually meant one model.
#'
#' @param line_col The user's colours, or `NULL`.
#' @param n Number of levels of the variable mapped to colour.
#' @param what Name of that variable, for the message.
#' @keywords internal
#' @noRd
.rce_check_line_col <- function(line_col, n, what) {
  if (is.null(line_col) || n <= 1L || length(line_col) >= n) return(invisible())
  warning("`line_col` has ", length(line_col), " colour",
          if (length(line_col) > 1L) "s" else "", " but colour separates ", n,
          " ", what, " on this figure, so they are recycled. Supply ", n,
          " colours to tell them apart.", call. = FALSE)
  invisible()
}


#' Validate a transparency value
#'
#' `alpha` reaches the device unchecked otherwise: a value above 1 saturates the
#' ribbon, and an `NA` or a vector draws nothing, in both cases without a
#' message. Checked here, where the argument still has a name, like `line_col`
#' and `lwd`.
#'
#' @param alpha The user's transparency, or `NULL` where the figure has no
#'   ribbon or fan to apply it to.
#' @return `alpha`, invisibly.
#' @keywords internal
#' @noRd
.rce_check_alpha <- function(alpha) {
  if (is.null(alpha)) return(invisible(alpha))
  if (length(alpha) != 1L || !is.numeric(alpha) || is.na(alpha) ||
      alpha < 0 || alpha > 1) {
    stop("`alpha` must be a single number between 0 and 1.", call. = FALSE)
  }
  invisible(alpha)
}


#' Add the standard colorblind-safe scales to a plot
#'
#' Discrete aesthetics (series identity, e.g. model) use the Okabe-Ito
#' qualitative palette; continuous aesthetics (ordered magnitude, e.g. year)
#' use viridis. Both are colorblind-safe.
#'
#' `line_col` overrides the default, supplying the palette for whichever
#' variable the plot maps to colour (see `?rceattle-plot-args`). On a continuous
#' colour aesthetic the colours are ramp anchors instead.
#'
#' @param p A ggplot.
#' @param discrete Use the discrete (Okabe-Ito) or continuous (viridis) scale.
#' @param aesthetics Which of `"colour"`/`"fill"` to add.
#' @param line_col User-supplied colours, or `NULL` for the package default.
#'   Recycled over the mapped variable's levels in plotting order.
#' @keywords internal
.rceattle_scale <- function(p, discrete = TRUE,
                            aesthetics = c("colour", "fill"),
                            line_col = NULL) {
  line_col <- .as_colour(line_col)

  # One closure for both aesthetics, so a user's colours cannot reach the line
  # but not its ribbon.
  scale_for <- function(aes_name) {
    if (is.null(line_col)) {
      return(if (discrete)
               ggplot2::discrete_scale(aes_name, palette = .okabe_ito_pal)
             else if (aes_name == "colour") ggplot2::scale_colour_viridis_c()
             else ggplot2::scale_fill_viridis_c())
    }
    if (discrete) {
      # Applied in level order by a palette function. Naming the values instead
      # renders every series grey50 whenever the names miss the data's levels,
      # and level order is the documented contract anyway.
      return(ggplot2::discrete_scale(
        aes_name,
        palette = function(n) rep(line_col, length.out = n),
        na.value = "grey50"))
    }
    # Continuous: ramp anchors. A single colour gives a flat ramp, which is what
    # `line_col = "black"` is asking for.
    if (length(line_col) == 1L) {
      return(if (aes_name == "colour")
               ggplot2::scale_colour_gradient(low = line_col, high = line_col)
             else ggplot2::scale_fill_gradient(low = line_col, high = line_col))
    }
    if (aes_name == "colour") ggplot2::scale_colour_gradientn(colours = line_col)
    else ggplot2::scale_fill_gradientn(colours = line_col)
  }

  for (a in c("colour", "fill")) {
    if (a %in% aesthetics) p <- p + scale_for(a)
  }
  p
}

#' Save a ggplot to `file` if a file stem is supplied; always return the plot
#'
#' @param p A ggplot.
#' @param file File stem (no extension), or `NULL` to skip saving.
#' @param suffix Appended to `file` to form the PNG name.
#' @param width,height Figure size in inches.
#' @keywords internal
.save_ggplot <- function(p, file = NULL, suffix = "plot",
                         width = 7, height = 6.5) {
  if (!is.null(file)) {
    ggplot2::ggsave(paste0(file, "_", suffix, ".png"), plot = p,
                    width = width, height = height, dpi = 300)
  }
  p
}


#' Facet by species, but only when there is more than one
#'
#' With a single species the strip label just names the only thing on the plot.
#' Faceting is what draws that strip, so the single-species case omits the facet
#' entirely rather than blanking the strip and leaving its gap.
#'
#' Returns `NULL` when there is nothing to facet, which ggplot2 accepts as a
#' no-op inside a `+` chain, so call sites read the same either way.
#'
#' @param df The plot data frame; its `Species` column supplies the count.
#' @param ... Passed to [ggplot2::facet_wrap()] (e.g. `scales`, `ncol`).
#' @return A `facet_wrap` layer, or `NULL` for a single-species model.
#' @keywords internal
#' @noRd
.facet_species <- function(df, ...) {
  n <- if (is.null(df$Species)) 0L else length(unique(stats::na.omit(df$Species)))
  if (n > 1L) ggplot2::facet_wrap(~ Species, ...) else NULL
}


# =============================================================================
# Standardized argument handling
# -----------------------------------------------------------------------------
# One resolver per shared argument, so every plot_*() reads `species`, `lwd`,
# `lty`, `line_col`, `minyr`/`maxyr`, `incl_proj`, `incl_mean` and `add_ci` the
# same way. Before these existed, plot_timeseries() was the only plotter that
# honoured any of them and the rest declared them as formals and ignored them.
# =============================================================================

#' Resolve the `species` / `spnames` pair to indices and labels
#'
#' `species` selects species by **index**; `spnames` supplies their **labels**.
#' Every plotter resolves the pair here, so `species` accepts each of the
#' spellings the assessment scripts use and normalizes it:
#'
#' * `NULL` or `"all"`, every species, in model order.
#' * numeric, indices, validated against `1:nspp`.
#' * logical, a mask over `1:nspp`.
#' * character matching `spnames`, selection by name.
#' * character matching *nothing*, back-compatibility: these are labels, not a
#'   selection, which is how `plot_selectivity()` and `plot_maturity()` used the
#'   argument. Emits a message pointing at `spnames` and keeps every species.
#' * character matching *some* names, a selection with typos. Warns, naming
#'   the unmatched entries, and keeps the matches.
#'
#' The returned `index` preserves the order the caller asked for, so
#' `species = c(3, 1)` yields facets in that order rather than model order.
#'
#' @param models A list of `Rceattle` fits (from `.as_plot_models()`).
#' @param species Species selection; see above.
#' @param spnames Species labels, length `nspp`. `NULL` takes the model's own.
#' @return `list(index = <integer>, spnames = <character, nspp>,
#'   labels = <character, length(index)>)`.
#' @keywords internal
#' @noRd
.resolve_species <- function(models, species = NULL, spnames = NULL) {
  dl   <- models[[1]]$data_list
  nspp <- dl$nspp
  if (is.null(nspp) || is.na(nspp)) nspp <- length(dl$spnames)
  if (!length(nspp) || nspp < 1L) {
    # Report the model, not the argument: a malformed fit would otherwise give
    # "`species` selected no species" for an argument the caller never passed.
    stop("This model reports no species (`data_list$nspp` is ",
         format(dl$nspp), ").", call. = FALSE)
  }

  # Default labels come from the model; a model without them gets positional
  # ones so the legend is never blank.
  default_names <- dl$spnames
  if (is.null(default_names) || length(default_names) != nspp) {
    default_names <- paste("Species", seq_len(nspp))
  }
  user_named <- !is.null(spnames)
  if (!user_named) spnames <- default_names
  # Recycling a short `spnames` would label species 3 with species 1's name.
  # Refuse instead.
  spnames <- as.character(spnames)
  if (length(spnames) != nspp) {
    stop("`spnames` has ", length(spnames), " name(s) but the model has ",
         nspp, " species. Supply one name per species.", call. = FALSE)
  }

  idx <- NULL
  if (is.null(species)) {
    idx <- seq_len(nspp)
  } else if (is.factor(species)) {
    # A factor of species names is a plausible input (e.g. carried out of a
    # data frame); treat it as the character vector it prints as.
    return(.resolve_species(models, as.character(species), spnames))
  } else if (is.character(species)) {
    if (length(species) == 1L && identical(tolower(species), "all")) {
      idx <- seq_len(nspp)
    } else {
      # Match against the model's own names as well as any relabelling, so
      # renaming species for display does not make them unselectable.
      hit <- match(species, spnames)
      hit[is.na(hit)] <- match(species[is.na(hit)], default_names)
      legacy_shape <- length(species) == nspp && !user_named
      if (all(!is.na(hit))) {
        idx <- hit
      } else if (all(is.na(hit)) && legacy_shape) {
        # The legacy spelling: `species` held one display label per species,
        # which is how plot_selectivity() and plot_maturity() used it. Only
        # this shape -- one string per species, no separate `spnames` -- is
        # read that way.
        spnames <- as.character(species)
        idx <- seq_len(nspp)
        message("`species` matched no species name and gives one string per ",
                "species, so it is being read as display labels. Pass labels ",
                "as `spnames` and a selection as `species`.")
      } else if (legacy_shape) {
        # One string per species, some matching and some not: it could be a
        # relabelling with a collision or a selection with a typo, and the two
        # give different figures. Say so rather than pick one.
        stop("`species` gives one string per species but only ",
             paste(species[!is.na(hit)], collapse = ", "), " match a species ",
             "name. Read as labels it would rename all ", nspp,
             "; read as a selection it would plot ", sum(!is.na(hit)),
             ". Pass labels as `spnames`, or fix ",
             paste(species[is.na(hit)], collapse = ", "), ".", call. = FALSE)
      } else if (all(is.na(hit))) {
        stop("`species` matched no species name: ",
             paste(species, collapse = ", "), ". This model has: ",
             paste(default_names, collapse = ", "), ".", call. = FALSE)
      } else {
        warning("`species` did not match: ",
                paste(species[is.na(hit)], collapse = ", "),
                ". Plotting the ", sum(!is.na(hit)), " that did.", call. = FALSE)
        idx <- hit[!is.na(hit)]
      }
    }
  } else if (is.logical(species)) {
    idx <- which(rep(species, length.out = nspp))
  } else if (is.numeric(species)) {
    idx <- as.integer(species)
  } else {
    stop("`species` must be species indices, names, a logical mask, or \"all\".",
         call. = FALSE)
  }

  if (length(idx) == 0L) {
    stop("`species` selected no species.", call. = FALSE)
  }
  bad <- idx[is.na(idx) | idx < 1L | idx > nspp]
  if (length(bad)) {
    stop("`species` out of range: ", paste(bad, collapse = ", "),
         ". This model has ", nspp, " species (valid indices 1:", nspp, ").",
         call. = FALSE)
  }

  list(index = idx, spnames = spnames, labels = spnames[idx])
}


#' Resolve an age to its bin index in a species' at-age array
#'
#' The at-age arrays are indexed by BIN, 1 .. `nages[sp]`, while a species' ages
#' run `minage[sp] .. minage[sp] + nages[sp] - 1`. The two coincide only at
#' `minage = 1`, so an `age` argument passed straight through as a subscript
#' plots a different age than the axis label claims for any other species, and
#' runs off the end of the array for an age the species does not hold.
#'
#' A species that has no such age is dropped rather than plotted wrong, since
#' the ages differ between species in one figure. That matches how
#' `.resolve_species()` handles a partly-matching selection.
#'
#' @param age The age to plot, on the species' own age scale.
#' @param species Species indices being drawn.
#' @param minage,nages The model's `minage` and `nages` vectors.
#' @param spnames Species labels, for the messages.
#' @param arg Name of the calling argument, for the messages.
#' @return `list(species = <the species that hold this age>,
#'   index = <their bin indices, same order>)`.
#' @keywords internal
#' @noRd
.rce_age_index <- function(age, species, minage, nages, spnames, arg = "age") {
  if (length(age) != 1L || !is.numeric(age) || !is.finite(age)) {
    stop("`", arg, "` must be a single age.", call. = FALSE)
  }
  idx  <- as.integer(age) - as.integer(minage[species]) + 1L
  keep <- idx >= 1L & idx <= nages[species]
  if (!any(keep)) {
    stop("No species has age ", age, ". Age ranges: ",
         paste(sprintf("%s %d-%d", spnames[species], minage[species],
                       minage[species] + nages[species] - 1L),
               collapse = "; "), ".", call. = FALSE)
  }
  if (any(!keep)) {
    warning("Age ", age, " is outside the age range of ",
            paste(sprintf("%s (%d-%d)", spnames[species[!keep]],
                          minage[species[!keep]],
                          minage[species[!keep]] + nages[species[!keep]] - 1L),
                  collapse = ", "), "; ",
            if (sum(!keep) == 1L) "it is" else "they are",
            " not drawn.", call. = FALSE)
  }
  list(species = species[keep], index = idx[keep])
}


#' Bin indices of the ages at or above `minage`, per species
#'
#' The `minage` argument of [plot_ration()] names an age ("age 3+"), not a bin,
#' so it is resolved against each species' own age vector. A species with no age
#' that old is dropped rather than summed over a shifted range.
#'
#' @param minage The lower age, on the species' own age scale.
#' @param species Species indices being drawn.
#' @param model_minage,nages The model's `minage` and `nages` vectors.
#' @param spnames Species labels, for the messages.
#' @return `list(species = <kept>, index = <list of bin-index vectors>)`.
#' @keywords internal
#' @noRd
.rce_age_plus_index <- function(minage, species, model_minage, nages, spnames) {
  if (length(minage) != 1L || !is.numeric(minage) || !is.finite(minage)) {
    stop("`minage` must be a single age.", call. = FALSE)
  }
  idx <- lapply(species, function(sp) {
    ages <- seq_len(nages[sp]) - 1L + as.integer(model_minage[sp])
    which(ages >= minage)
  })
  keep <- lengths(idx) > 0L
  if (!any(keep)) {
    stop("No species has an age at or above ", minage, ". Age ranges: ",
         paste(sprintf("%s %d-%d", spnames[species], model_minage[species],
                       model_minage[species] + nages[species] - 1L),
               collapse = "; "), ".", call. = FALSE)
  }
  if (any(!keep)) {
    warning("No age at or above ", minage, " in ",
            paste(spnames[species[!keep]], collapse = ", "),
            "; not drawn.", call. = FALSE)
  }
  list(species = species[keep], index = idx[keep])
}


#' Line width / type / transparency parameters for a `geom_line()`
#'
#' A value that does not vary is passed as a fixed geom parameter, so the
#' uniform default adds no spurious legend; one that varies is mapped to the
#' keying variable and given a scale.
#'
#' `lwd` keeps the base-graphics convention, where the default of 3 renders as a
#' standard-weight ggplot line, the `lwd / 3` below. Do not change that ratio:
#' every default figure in the package and in the assessment scripts is drawn at
#' `linewidth = 1`, so `lwd = 3` must keep producing exactly that.
#'
#' @param lwd Line width(s), base-graphics scale (default 3 == linewidth 1).
#' @param lty Line type(s).
#' @param alpha Fixed transparency, or `NULL` to leave it unset.
#' @param lwd_by,lty_by Name of the data column keying each aesthetic when it
#'   varies, or `NULL` if the plot has nothing to key it to. Values are applied
#'   to that column's levels in plotting order.
#' @param lwd_n,lty_n Number of levels the keying column actually has. A column
#'   with one level cannot hold a varying value, so saying so lets the caller
#'   be told rather than having the extra values dropped in silence, line type
#'   keys on sex, and most models here are sex-combined.
#' @param lty_in_aes The plot's own `aes()` already maps line type (as
#'   `plot_ration()` maps it to sex). A vector then supplies one value per
#'   level; a single non-default value is applied to every level through a
#'   one-value palette, since passing it as a fixed geom parameter would
#'   override the mapping and drop its legend. The default `lty = 1` is what the
#'   first level already draws, so it is left alone and no default figure moves.
#' @return `list(mapping = <aes or NULL>, args = <fixed geom params>,
#'   scales = <list of scale objects>)`.
#' @keywords internal
#' @noRd
.rce_line_params <- function(lwd = 3, lty = 1, alpha = NULL,
                             lwd_by = NULL, lty_by = NULL,
                             lwd_n = NULL, lty_n = NULL,
                             lty_in_aes = FALSE) {
  args   <- list()
  scales <- list()

  # An NA or negative width does not draw: the line silently disappears, or the
  # device errors deep inside the render. Catch it where the argument has a name.
  if (is.null(lwd) || !length(lwd)) lwd <- 3
  if (is.null(lty) || !length(lty)) lty <- 1
  if (!is.numeric(lwd) || anyNA(lwd) || any(lwd < 0)) {
    stop("`lwd` must be non-negative numbers; got ",
         paste(utils::head(format(lwd), 5), collapse = ", "), ".", call. = FALSE)
  }
  if (anyNA(lty)) stop("`lty` must not be NA.", call. = FALSE)

  # An aesthetic can only be mapped if the plot has a variable to key it to AND
  # that variable separates more than one series. Otherwise the first value is
  # used throughout, and the caller is told which of the two it hit.
  has_lwd_key <- !is.null(lwd_by) && (is.null(lwd_n) || lwd_n > 1L)
  has_lty_key <- !is.null(lty_by) && (is.null(lty_n) || lty_n > 1L)
  map_lwd <- length(unique(lwd)) > 1L && has_lwd_key
  map_lty <- length(unique(lty)) > 1L && has_lty_key

  no_key_msg <- function(arg, by, n, value) {
    if (is.null(by)) {
      paste0("`", arg, "` varies but this plot draws one line ",
             if (arg == "lwd") "width" else "type", "; using ", arg, " = ",
             value, " throughout.")
    } else {
      paste0("`", arg, "` varies but this plot keys line ",
             if (arg == "lwd") "width" else "type", " on ", by,
             ", which has ", n, " level; using ", arg, " = ", value,
             " throughout.")
    }
  }
  if (length(unique(lwd)) > 1L && !map_lwd) {
    warning(no_key_msg("lwd", lwd_by, lwd_n, lwd[1]), call. = FALSE)
  }
  if (length(unique(lty)) > 1L && !map_lty) {
    warning(no_key_msg("lty", lty_by, lty_n, lty[1]), call. = FALSE)
  }

  # Values are applied in level order by a palette function rather than by
  # naming them: a named scale_*_manual() whose names miss the data's levels
  # renders every line blank, and the caller cannot be expected to know the
  # factor's levels.
  #
  # A fixed value would override an aesthetic the plot's own aes() maps -- on
  # plot_ration() that flattens female and male onto one line -- so where the
  # mapping already exists, a single value leaves it alone.
  if (map_lwd) {
    scales <- c(scales, list(ggplot2::discrete_scale(
      "linewidth", palette = function(n) rep(lwd, length.out = n) / 3)))
  } else {
    args$linewidth <- lwd[1] / 3
  }
  if (map_lty) {
    scales <- c(scales, list(ggplot2::discrete_scale(
      "linetype", palette = function(n) rep(lty, length.out = n))))
  } else if (!lty_in_aes) {
    args$linetype <- lty[1]
  } else if (!isTRUE(all.equal(as.vector(lty[1]), 1))) {
    # The plot maps line type itself and `lty` is a single non-default value, so
    # it applies to every level. Leaving it unset here is what silently dropped
    # `plot_ration(fit, lty = 2)`: the argument had no effect and said nothing.
    # Collapsing a key that does separate something is worth a word, since the
    # figure then draws those levels alike.
    if (!is.null(lty_n) && lty_n > 1L) {
      warning("`lty` is one value but this plot keys line type on ", lty_by,
              ", which has ", lty_n, " levels; they are drawn alike. Supply ",
              lty_n, " line types to tell them apart.", call. = FALSE)
    }
    scales <- c(scales, list(ggplot2::discrete_scale(
      "linetype", palette = function(n) rep(lty[1], length.out = n))))
  }
  if (!is.null(alpha)) args$alpha <- alpha

  # Built explicitly rather than by deleting elements, so the quosures capture
  # the column names from this frame.
  mapping <- if (map_lwd && map_lty) {
    ggplot2::aes(linewidth = .data[[lwd_by]], linetype = .data[[lty_by]])
  } else if (map_lwd) {
    ggplot2::aes(linewidth = .data[[lwd_by]])
  } else if (map_lty) {
    ggplot2::aes(linetype = .data[[lty_by]])
  } else {
    NULL
  }

  list(mapping = mapping, args = args, scales = scales)
}


#' Add a `geom_line()` built from `.rce_line_params()`, plus its scales
#'
#' Adds the geom and its scales, which every converted plotter needs together.
#' @param p A ggplot.
#' @param lp The list returned by `.rce_line_params()`.
#' @param ... Further arguments to [ggplot2::geom_line()].
#' @keywords internal
#' @noRd
.rce_add_line <- function(p, lp, ...) {
  p <- p + do.call(ggplot2::geom_line,
                   c(list(mapping = lp$mapping), lp$args, list(...)))
  for (s in lp$scales) p <- p + s
  p
}


#' Restrict a plot frame to a year range
#'
#' Narrows a plot to `minyr:maxyr` by dropping rows, so that a `free_y` panel
#' rescales to what is actually shown. Clipping the axis with
#' `coord_cartesian()` leaves the y scale trained on the hidden years, which on
#' a series that has grown by two orders of magnitude squeezes the requested
#' window into a few pixels at the bottom of the panel.
#'
#' Also used to thin a year fan (selectivity, mortality-at-age), where Year is
#' the colour or grouping variable rather than the x axis.
#'
#' An empty frame is an error, not a blank panel; a window that only empties
#' *some* panels just loses those rows.
#'
#' @param df A plot data frame.
#' @param minyr,maxyr First / last year to keep, or `NULL` for no limit.
#' @param year Name of the year column.
#' @keywords internal
#' @noRd
.rce_year_filter <- function(df, minyr = NULL, maxyr = NULL, year = "Year") {
  if (is.null(df[[year]])) return(df)
  if (!is.null(minyr) && (length(minyr) != 1L || is.na(minyr))) {
    stop("`minyr` must be a single year.", call. = FALSE)
  }
  if (!is.null(maxyr) && (length(maxyr) != 1L || is.na(maxyr))) {
    stop("`maxyr` must be a single year.", call. = FALSE)
  }
  keep <- rep(TRUE, nrow(df))
  if (!is.null(minyr)) keep <- keep & df[[year]] >= minyr
  if (!is.null(maxyr)) keep <- keep & df[[year]] <= maxyr
  # A comparison against an NA year gives NA, which subsets to a fabricated
  # all-NA row rather than dropping it.
  keep[is.na(keep)] <- FALSE
  if (!any(keep)) {
    stop("No years left to plot between minyr = ",
         if (is.null(minyr)) "NULL" else minyr, " and maxyr = ",
         if (is.null(maxyr)) "NULL" else maxyr, ".", call. = FALSE)
  }
  df[keep, , drop = FALSE]
}


#' Dashed divider at the last hindcast year
#'
#' Marks where the projection starts. Consolidates the identical `geom_vline()`
#' that six plotters held inline.
#'
#' Takes the **latest** hindcast year across the models rather than the last
#' model's. On a retrospective peel list the models end in different years, and
#' keying on whichever happens to be last puts the divider mid-hindcast and
#' labels real data as projection.
#'
#' Drawn only when the boundary is inside the window being plotted. A
#' `geom_vline()` trains the x scale like any other layer, so a divider outside
#' a `minyr`/`maxyr` window would stretch the axis back to it and reintroduce
#' the empty span the window was asked for to remove.
#'
#' @param models A list of `Rceattle` fits.
#' @param incl_proj Whether projection years are being drawn.
#' @param minyr,maxyr The plotted year window, or `NULL` for no limit.
#' @return A `geom_vline` layer, or `NULL`.
#' @keywords internal
#' @noRd
.rce_proj_divider <- function(models, incl_proj = FALSE,
                              minyr = NULL, maxyr = NULL) {
  if (!isTRUE(incl_proj)) return(NULL)
  endyr <- vapply(models, function(m) {
    e <- m$data_list$endyr
    if (is.null(e) || !length(e)) NA_real_ else as.numeric(e[1])
  }, numeric(1))
  if (all(is.na(endyr))) return(NULL)
  at <- max(endyr, na.rm = TRUE)
  if (!is.null(minyr) && !is.na(minyr) && at < minyr) return(NULL)
  if (!is.null(maxyr) && !is.na(maxyr) && at > maxyr) return(NULL)
  ggplot2::geom_vline(xintercept = at, linetype = 2, colour = "grey50")
}


#' Last hindcast year of each model, named by its plot label
#'
#' Feeds `.rce_mean_line()`, which has to cut each model at its own `endyr`.
#'
#' @param models A list of `Rceattle` fits.
#' @param labels Plot labels for those models, from [.model_labels()].
#' @keywords internal
#' @noRd
.rce_model_endyr <- function(models, labels) {
  stats::setNames(vapply(models, function(m) {
    e <- m$data_list$endyr
    if (is.null(e) || !length(e)) NA_real_ else as.numeric(e[1])
  }, numeric(1)), labels)
}


#' Horizontal long-term mean line(s)
#'
#' `incl_mean` was documented as "include time series mean as horizontal line"
#' on four plotters and implemented on none. The mean is taken over **hindcast
#' years only**: a mean that folded in the projection would not be a historical
#' reference, and would move when the projection horizon changed. `hind_endyr`
#' is therefore required, not optional, a caller that forgot it would silently
#' get the projection-contaminated mean this exists to avoid.
#'
#' The mean line takes a colour only when the caller says colour encodes the
#' model. Inferring that from `by` is wrong: this helper's own callers map
#' colour to predator or to species, and adding a `colour = Model` layer to such
#' a plot trains model names into that scale's legend.
#'
#' Models in one figure can end in different years, a retrospective peel is
#' the usual case, so `hind_endyr` may be a vector named by model. Cutting
#' every model at the first one's `endyr` averages a peel over years it never
#' fitted, and the answer then depends on the order of the list.
#'
#' @param df The plot data frame.
#' @param incl_mean Draw it?
#' @param by Columns defining a group (a facet and/or a series), e.g.
#'   `c("Species", "Model")`.
#' @param value Name of the value column.
#' @param hind_endyr Last hindcast year. A single year, or one per model named
#'   by the levels of `model`. Pass `NA` to average every row deliberately.
#' @param year Name of the year column.
#' @param model Name of the model column, used to look up a per-model
#'   `hind_endyr`.
#' @param colour_by Column the plot maps to colour, or `NULL` if the mean line
#'   should be drawn in a neutral grey.
#' @return A `geom_hline` layer, or `NULL`.
#' @keywords internal
#' @noRd
.rce_mean_line <- function(df, incl_mean = FALSE, by = NULL, value = "value",
                           hind_endyr, year = "Year", model = "Model",
                           colour_by = NULL) {
  if (!isTRUE(incl_mean)) return(NULL)
  if (missing(hind_endyr)) {
    stop("`hind_endyr` is required: the mean is over hindcast years only.",
         call. = FALSE)
  }
  hind <- df
  if (!is.null(hind_endyr) && !all(is.na(hind_endyr)) && !is.null(hind[[year]])) {
    cut_at <- if (length(hind_endyr) > 1L && !is.null(hind[[model]])) {
      # One endyr per model, matched by name.
      unname(hind_endyr[as.character(hind[[model]])])
    } else {
      hind_endyr[1]
    }
    keep <- !is.na(hind[[year]]) & !is.na(cut_at) & hind[[year]] <= cut_at
    hind <- hind[keep, , drop = FALSE]
  }
  hind <- hind[!is.na(hind[[value]]), , drop = FALSE]
  if (nrow(hind) == 0L) return(NULL)

  # A grouping column must survive into the layer's data for faceting to place
  # the line, so drop any the frame does not carry.
  by <- by[by %in% names(hind)]
  by <- setdiff(by, ".mean")
  agg <- if (length(by)) {
    stats::aggregate(hind[[value]], by = hind[by], FUN = mean)
  } else {
    data.frame(x = mean(hind[[value]]))
  }
  names(agg)[ncol(agg)] <- ".mean"

  colour_by <- if (!is.null(colour_by) && colour_by %in% names(agg)) colour_by
  args <- list(data = agg, inherit.aes = FALSE, linetype = 3)
  if (is.null(colour_by)) {
    args$mapping <- ggplot2::aes(yintercept = .data$.mean)
    args$colour  <- "grey40"
  } else {
    args$mapping <- ggplot2::aes(yintercept = .data$.mean,
                                 colour = .data[[colour_by]])
  }
  do.call(ggplot2::geom_hline, args)
}


#' Standard errors for an ADREPORTed series, by name
#'
#' Factored out of `plot_timeseries()` so any plotter drawing a confidence
#' interval reads `sdrep` the same way.
#'
#' The block length must match `n_total` **exactly**. `sdrep$value` holds the
#' whole flattened series, so taking the first `n_need` of a block of unverified
#' length silently returns the standard errors of different cells, an interval
#' that is wrong rather than absent.
#'
#' Taking a leading slice is legitimate only when the caller knows the block's
#' full shape, which is why `n_total` has to be stated rather than inferred. The
#' species-by-year series are flattened column-major with the hindcast years
#' first, so the first `nspp * nyrs_hindcast` values are exactly the hindcast,
#' that is the one prefix any caller here needs.
#'
#' @param model An `Rceattle` fit.
#' @param name Name of the ADREPORTed series.
#' @param n_need How many leading values the caller wants.
#' @param n_total Length the whole block must have. Defaults to `n_need`, i.e.
#'   no slicing.
#' @return `n_need` standard errors, or `NULL` when the fit does not hold a
#'   block of exactly `n_total` (no `sdreport`, or a `REPORT()`-only series).
#' @keywords internal
#' @noRd
.rce_series_sd <- function(model, name, n_need, n_total = n_need) {
  sdrep <- model$sdrep
  if (is.null(sdrep) || is.null(sdrep$value)) return(NULL)
  rows <- which(names(sdrep$value) == name)
  if (length(rows) != n_total || n_need > n_total) return(NULL)
  out <- sdrep$sd[rows]
  # A species with input numbers-at-age reports its input recruits as R
  # (rename_output()), but the sdreport's R is the placeholder curve's: no SE to show.
  nspp  <- model$data_list$nspp
  fixed <- which((model$data_list$estDynamics %||% 0) > 0)
  if (name %in% c("R", "log_R") && length(fixed) && length(out) %% nspp == 0) {
    out[((seq_along(out) - 1) %% nspp + 1) %in% fixed] <- NA
  }
  out[seq_len(n_need)]
}


#' Can this quantity hold a confidence interval at all?
#'
#' Reads the same `.RCEATTLE_QUANTITIES` registry that `as.data.frame.Rceattle()`
#' uses, so the figure and the table agree on which series have standard errors.
#' A quantity absent from the registry is a derived combination (ration is
#' consumption x biomass; the M2 proportions are ratios) and has none.
#'
#' @param quantity REPORT name of the plotted quantity.
#' @keywords internal
#' @noRd
.rce_quantity_adreport <- function(quantity) {
  isTRUE(.RCEATTLE_QUANTITIES[[quantity]]$adreport)
}


#' Decline an `add_ci` request that the model cannot support
#'
#' Says so once rather than ignoring the request. Kept quiet for a quantity the
#' model was never going to report: `F_spp` and the other `REPORT()`-only series
#' would otherwise warn on every fit.
#'
#' @param add_ci The user's request.
#' @param quantity Name of the plotted quantity, for the message.
#' @param reason Why it cannot be drawn.
#' @param warn Emit the warning? Pass `FALSE` where the quantity is
#'   `REPORT()`-only by design.
#' @return `FALSE`, always, so the caller can write `add_ci <- .rce_no_ci(...)`.
#' @keywords internal
#' @noRd
.rce_no_ci <- function(add_ci, quantity, reason, warn = TRUE) {
  if (isTRUE(add_ci) && isTRUE(warn)) {
    warning("`add_ci = TRUE` is not available for ", quantity, ": ", reason,
            ". Drawing no confidence interval.", call. = FALSE)
  }
  FALSE
}


#' Shared plotting arguments
#'
#' A common vocabulary, documented here once and inherited by the plotters that
#' take it: [plot_timeseries()] and its wrappers ([plot_biomass()], [plot_ssb()],
#' [plot_recruitment()], the depletions, [plot_exploitable_biomass()],
#' [plot_f()]), the predation plotters ([plot_b_eaten()], [plot_b_eaten_prop()],
#' [plot_m_at_age()], [plot_m2_at_age_prop()], [plot_ration()]), and
#' [plot_selectivity()]. Each argument means the same thing wherever it appears,
#' but not every plotter takes every one, `incl_mean` is on the predation
#' plotters, `add_ci` only where the quantity holds standard errors, and
#' `alpha` only where the figure has a ribbon or a fan. The remaining
#' `plot_*()` functions still take their own arguments; see each one's help.
#'
#' @section How `line_col` and `lty` are applied:
#'
#' They supply values for whichever **discrete variable the plot already
#' encodes with that aesthetic**, matched in level order, which is not always
#' the model. Each function's help says what its own figure separates.
#'
#' Where colour encodes a continuous variable, the year fan in
#' [plot_selectivity()], `line_col` supplies the ramp anchors instead: one
#' colour draws the fan in that colour, several interpolate between them.
#'
#' @param Rceattle A single [fit_mod()] object or a list of them (overlaid).
#' @param file Optional file stem; the figure is written to
#'   `<file>_<suffix>.png` if given.
#' @param model_names Legend labels for the models.
#' @param species Species to include, as indices (`c(1, 3)`), names, a logical
#'   mask, or `"all"`. Default `NULL` plots every species. Species **labels**
#'   belong in `spnames`; a character `species` that matches no species name is
#'   read as labels for back-compatibility.
#' @param spnames Species labels, length `nspp`. Default: the model's own.
#' @param line_col Line colours; names, hex, or base-graphics integers. `NULL`
#'   uses the colorblind-safe Okabe-Ito palette. Applied to whichever variable
#'   the figure separates by colour, in legend order. Too few colours are
#'   recycled, with a warning.
#' @param lwd Line width on the base-graphics scale: the default `3` renders as
#'   a standard-weight ggplot line. A vector varies it across series.
#' @param lty Line type. A vector varies it across the levels of whatever the
#'   figure separates by line type.
#' @param alpha Transparency of confidence ribbons and shaded areas, between 0
#'   and 1.
#' @param add_ci Add a 95% confidence interval. Only available where the
#'   plotted quantity holds standard errors; warns and draws none otherwise.
#' @param minyr,maxyr First / last year to plot.
#' @param incl_proj Include the projection years, with a dashed divider at the
#'   last hindcast year.
#' @param incl_mean Add a horizontal line at the hindcast mean of each series.
#' @param width,height Saved figure size in inches.
#' @param right_adj Ignored. Base-graphics leftover: the figure widened its
#'   right margin to fit the legend. Set margins on the returned ggplot instead.
#' @param top_adj Ignored. Base-graphics leftover; see `right_adj`.
#' @param mod_cex Ignored. Base-graphics leftover: legend text size. Use
#'   `p + ggplot2::theme(legend.text = ggplot2::element_text(size = ...))`.
#' @param legend.pos Ignored. Base-graphics leftover: legend placement. Use
#'   `p + ggplot2::theme(legend.position = "right")`.
#' @param single.plots Ignored. Base-graphics leftover: one device per panel.
#'   The ggplot facets instead.
#' @param theta Ignored. Base-graphics leftover: viewing angle of a `persp`
#'   surface, which is no longer drawn.
#' @param ymax Ignored. Base-graphics leftover: y-axis maximum. Use
#'   `p + ggplot2::coord_cartesian(ylim = c(0, ymax))`.
#' @param cex Ignored. Base-graphics leftover: point expansion.
#' @return A `ggplot` object.
#' @name rceattle-plot-args
NULL


#' The selectivity standard errors a fit reports, or NULL
#'
#' `fit_control(selectivity_se = TRUE)` adds `log_sel_at_age` to the sdreport
#' and `log_sel_at_age_index` to the reported quantities. Both are needed: the
#' index says which (fleet, sex, age, year) each sdreport row belongs to, and
#' its length depends on the model, since Fixed fleets and mirrored fleets are
#' left out.
#'
#' @param model One fitted `Rceattle` object.
#' @return A data frame of `Fleet`, `Sex`, `Age`, `Year`, `log_sel`, `sd`, or
#'   `NULL` where the fit has no sdreport or was not asked for the errors.
#' @keywords internal
#' @noRd
.rce_sel_se <- function(model) {
  sdrep <- model$sdrep
  idx   <- model$quantities$log_sel_at_age_index
  if (is.null(sdrep) || is.null(sdrep$value) || is.null(idx)) return(NULL)
  k <- names(sdrep$value) == "log_sel_at_age"
  if (!any(k) || sum(k) != nrow(idx)) return(NULL)
  out <- data.frame(Fleet = as.integer(idx[, 1]), Sex = as.integer(idx[, 2]),
                    Age   = as.integer(idx[, 3]), Year = as.integer(idx[, 4]),
                    log_sel = unname(sdrep$value[k]), sd = unname(sdrep$sd[k]),
                    stringsAsFactors = FALSE)
  # The plotter asks for one fleet, sex and year at a time, so index once here.
  attr(out, "by") <- split(seq_len(nrow(out)),
                           paste(out$Fleet, out$Sex, out$Year))
  out
}


#' The fleets whose reported selectivity a given fleet may borrow
#'
#' The template reports one row set per selectivity group, on the group's lead
#' fleet, so a mirrored fleet has no rows of its own. The lead is resolved with
#' `rearrange_data()`'s own `.group_lead()`, on the same key it uses,
#' `Selectivity_index` and `Selectivity` together, with an `Off` fleet never
#' leading, so the plotter and the template cannot disagree about who leads.
#' Recomputed rather than read off `data_list`, which is the pre-`rearrange_data()`
#' list and holds no `flt_sel_lead`. Reading `Selectivity_index` as a fleet code
#' instead would miss a group keyed on a value that is nobody's `Fleet_code`, and
#' silently drop the band.
#'
#' @param dl One model's `data_list`.
#' @param i Row index of the fleet, which is also its `Fleet_code`.
#' @return Fleet codes in the same selectivity group, the lead first.
#' @keywords internal
#' @noRd
.rce_sel_lead <- function(dl, i) {
  fc  <- dl$fleet_control
  flt <- fc$Fleet_code[i]
  if (is.null(fc$Selectivity_index) || is.null(fc$Selectivity)) return(flt)
  off  <- as.character(fc$Fleet_type) %in% c("Off", "0")
  lead <- .group_lead(paste(fc$Selectivity_index, fc$Selectivity), off)
  key  <- paste(fc$Selectivity_index, fc$Selectivity)
  unique(c(fc$Fleet_code[key == key[i] & lead == 1L], flt))
}


#' Rows of the `.rce_sel_se()` frame for one fleet, sex and year
#'
#' Uses the lookup `.rce_sel_se()` built, so a 16-fleet model does not rescan a
#' 30,000-row frame once per panel and year.
#'
#' @param se The `.rce_sel_se()` frame.
#' @param flt,sex,year Fleet code, sex index and calendar year.
#' @return Integer row positions, possibly none.
#' @keywords internal
#' @noRd
.rce_sel_rows <- function(se, flt, sex, year) {
  by <- attr(se, "by")
  if (!is.null(by)) return(by[[paste(flt, sex, year)]])
  which(se$Fleet == flt & se$Sex == sex & se$Year == year)
}


#' One fleet-sex-year's selectivity interval, on the natural scale
#'
#' Returned in the order of `bins`, with `NA` for any age the fit did not
#' report, so the caller can bind it to the curve without checking lengths.
#' `exp()` of a log-scale interval, so it is positive and right-skewed.
#'
#' A mirrored fleet holds no rows of its own, so it borrows its lead's, but
#' only where the two curves actually agree. `data_check()` merely warns when
#' fleets sharing a `Selectivity_index` differ in a shaping column, and a
#' differing `Sel_norm_bin` alone rescales the curve, so sharing the parameter
#' block is not on its own enough to share the band.
#'
#' @param se The `.rce_sel_se()` frame, or NULL.
#' @param flt,sex Fleet code and sex index.
#' @param bins Absolute ages, in plotting order.
#' @param year Calendar year.
#' @param lead Fleet codes whose rows this fleet may borrow (`.rce_sel_lead()`).
#' @param curve This fleet's own selectivity over `bins`, to check a borrowed
#'   band against. `NULL` borrows unchecked.
#' @return A list of `lower` and `upper`, or `NULL` if no band applies.
#' @keywords internal
#' @noRd
.rce_sel_ci <- function(se, flt, sex, bins, year, lead = flt, curve = NULL) {
  if (is.null(se)) return(NULL)
  band <- function(hit) {
    m <- match(bins, se$Age[hit])
    list(lower = exp(se$log_sel[hit][m] - 1.96 * se$sd[hit][m]),
         upper = exp(se$log_sel[hit][m] + 1.96 * se$sd[hit][m]),
         fit   = exp(se$log_sel[hit][m]))
  }
  own <- .rce_sel_rows(se, flt, sex, year)
  if (length(own)) return(band(own)[c("lower", "upper")])

  for (ld in setdiff(lead[!is.na(lead)], flt)) {
    hit <- .rce_sel_rows(se, ld, sex, year)
    if (!length(hit)) next
    b <- band(hit)
    # Compared where both are present; an age the lead did not report stays NA
    # in the band and is not evidence the curves differ.
    if (!is.null(curve)) {
      both <- is.finite(b$fit) & is.finite(curve)
      if (!any(both)) next
      if (!isTRUE(all.equal(b$fit[both], curve[both], tolerance = 1e-8))) next
    }
    return(b[c("lower", "upper")])
  }
  NULL
}


#' Truncate residuals onto the fixed bubble-size scale
#'
#' `plot_comp()` and `.osa_bubble_plot()` both draw residual bubbles on a fixed
#' `[0, .RCE_BUBBLE_MAX]` size scale so two figures compare by eye. Truncation
#' has to happen here rather than being left to the scale's `limits`, because
#' `scale_size_continuous()` sets an out-of-bounds value to `NA` and drops the
#' point silently, the largest residuals on the panel would be the ones that
#' vanished.
#'
#' The warning names the count and the largest magnitude rather than every
#' offending value: an upweighted fleet's residuals scale by `sqrt(w)`, so this
#' can fire on hundreds of points at once.
#'
#' @param x Residuals, signed.
#' @param what Figure name, to open the warning.
#' @return `x` with magnitudes above the cap pulled back to it, sign kept.
#' @keywords internal
#' @noRd
.rce_truncate_resid <- function(x, what) {
  big <- which(abs(x) > .RCE_BUBBLE_MAX)
  if (length(big)) {
    warning(what, ": ", length(big), " residual(s) beyond +/-", .RCE_BUBBLE_MAX,
            " truncated for plotting (largest |residual| ",
            sprintf("%.2f", max(abs(x[big]))), "). The bubble area understates ",
            "the misfit at those points.", call. = FALSE)
    x[big] <- .RCE_BUBBLE_MAX * sign(x[big])
  }
  x
}
