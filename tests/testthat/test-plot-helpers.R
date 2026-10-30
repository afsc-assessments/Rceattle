# Unit tests for the shared plotting-argument helpers (R/7-plot_helpers.R).
#
# These are what make `line_col`, `lwd`, `lty`, `alpha` and `species` behave the
# same way in every plot_*() function. They need no fitted model, so they stay
# unguarded (fast) while the plotters that use them are skip_on_cran.
#
# The load-bearing assertion is the `lwd / 3` bridge: every plotter converted to
# .rce_line_params() previously hardcoded `linewidth = 1`, so the default
# `lwd = 3` must keep producing exactly that or every existing figure changes.

# A minimal stand-in for a fitted model: .resolve_species() only reads
# `data_list$nspp` and `data_list$spnames`.
fake_models <- function(nspp = 3L, spnames = c("Pollock", "Cod", "ATF")) {
  list(list(data_list = list(nspp = nspp, spnames = spnames)))
}

# Colours of the first layer of a built plot.
built_colours <- function(p) {
  unique(ggplot2::ggplot_build(p)$data[[1]]$colour)
}

demo_df <- data.frame(
  x     = 1:4,
  y     = 1:4,
  Model = factor(rep(c("a", "b"), each = 2)),
  Year  = 1:4
)


# --- .as_colour ---------------------------------------------------------------
testthat::test_that(".as_colour resolves base-graphics colour specifications", {
  # `line_col = 1` and `line_col = c(1, 2, 4)` are how the real assessment
  # scripts pass colours.
  testthat::expect_equal(.as_colour(1), "black")
  testthat::expect_equal(.as_colour(c(1, 2, 4)),
                         grDevices::palette()[c(1, 2, 4)])
  testthat::expect_equal(.as_colour("black"), "black")
  testthat::expect_equal(.as_colour("#0072B2"), "#0072B2")
  testthat::expect_null(.as_colour(NULL))
  # Base treats 0 as transparent rather than as palette entry 0.
  testthat::expect_true(is.na(.as_colour(0)))
  # Recycling matches base: index 9 wraps around an 8-colour palette.
  testthat::expect_equal(.as_colour(length(grDevices::palette()) + 1L),
                         .as_colour(1))
  testthat::expect_equal(.as_colour(factor("red")), "red")
})

testthat::test_that(".as_colour refuses colours that would silently not draw", {
  # An NA colour renders as a missing line, so without these guards a bad
  # `line_col` produces a blank panel and no message.
  testthat::expect_null(.as_colour(numeric(0)))
  testthat::expect_null(.as_colour(character(0)))
  testthat::expect_error(.as_colour(-1), "non-negative")
  testthat::expect_error(.as_colour(NA), "NA")
  testthat::expect_error(.as_colour("notacolour"), "not a valid colour")
  testthat::expect_error(.as_colour(TRUE), "colour names, hex codes")
})


# --- .resolve_species ---------------------------------------------------------
testthat::test_that(".resolve_species accepts every spelling of a selection", {
  m <- fake_models()

  testthat::expect_equal(.resolve_species(m)$index, 1:3)
  testthat::expect_equal(.resolve_species(m, "all")$index, 1:3)
  testthat::expect_equal(.resolve_species(m, c(TRUE, FALSE, TRUE))$index, c(1L, 3L))
  testthat::expect_equal(.resolve_species(m, c("Cod", "ATF"))$index, c(2L, 3L))

  # The caller's order is preserved, so `species = c(3, 1)` faceted in that
  # order rather than silently sorting back to model order.
  res <- .resolve_species(m, c(3, 1))
  testthat::expect_equal(res$index, c(3L, 1L))
  testthat::expect_equal(res$labels, c("ATF", "Pollock"))
})

testthat::test_that(".resolve_species reads a full label vector as labels", {
  # Back-compatibility: plot_selectivity() and plot_maturity() used `species`
  # for display labels, and the assessment scripts still call them that way.
  # Only this exact shape -- one string per species, no separate `spnames` --
  # is read that way.
  m <- fake_models()
  ebs <- c("Walleye pollock", "Pacific cod", "Arrowtooth flounder")

  testthat::expect_message(res <- .resolve_species(m, ebs), "display labels")
  testthat::expect_equal(res$index, 1:3)      # nothing is dropped
  testthat::expect_equal(res$labels, ebs)     # the strings become the labels
})

testthat::test_that(".resolve_species refuses to relabel every species from a typo", {
  # The dangerous case: a mistyped selection must not silently plot every
  # species under the user's misspelling. A figure captioned with the wrong
  # species name is worse than no figure.
  m <- fake_models()

  testthat::expect_error(.resolve_species(m, "Codd"), "matched no species name")
  testthat::expect_error(.resolve_species(m, c("Pollok", "Codd")),
                         "matched no species name")
  # Partial match is a selection with a typo: keep what matched, say what did not.
  testthat::expect_warning(res <- .resolve_species(m, c("Cod", "Halibut")),
                           "did not match")
  testthat::expect_equal(res$index, 2L)
})

testthat::test_that(".resolve_species will not guess between relabelling and selecting", {
  # One string per species with SOME matching is ambiguous: read as labels it
  # renames all three, read as a selection it plots one. Both are plausible and
  # they give different figures, so neither is chosen silently.
  m <- fake_models()
  testthat::expect_error(.resolve_species(m, c("Pollock", "Zzz", "Yyy")),
                         "one string per species")
})

testthat::test_that(".resolve_species can still select by the model's own names after relabelling", {
  # Renaming species for display must not make them unselectable, and must not
  # fail by silently plotting everything.
  m <- fake_models()
  res <- .resolve_species(m, species = "Cod", spnames = c("A", "B", "C"))
  testthat::expect_equal(res$index, 2L)
  testthat::expect_equal(res$labels, "B")
})

testthat::test_that(".resolve_species rejects a wrong-length spnames", {
  # Recycling would label species 3 with species 1's name -- a plausible-looking
  # wrong answer on an assessment figure.
  m <- fake_models()
  testthat::expect_error(.resolve_species(m, spnames = c("A", "B")),
                         "one name per species")
  testthat::expect_error(.resolve_species(m, spnames = c("A", "B", "C", "D")),
                         "one name per species")
})

testthat::test_that(".resolve_species errors on bad input", {
  m <- fake_models()
  testthat::expect_error(.resolve_species(m, 5), "out of range")
  testthat::expect_error(.resolve_species(m, 0), "out of range")
  testthat::expect_error(.resolve_species(m, integer(0)), "selected no species")
  testthat::expect_error(.resolve_species(m, list(1)), "must be species indices")

  # A malformed model must blame itself, not an argument the caller never passed.
  testthat::expect_error(
    .resolve_species(list(list(data_list = list(nspp = 0L, spnames = NULL)))),
    "reports no species")
})

testthat::test_that(".resolve_species accepts a factor of species names", {
  m <- fake_models()
  testthat::expect_equal(.resolve_species(m, factor(c("Cod", "ATF")))$index,
                         c(2L, 3L))
})

testthat::test_that(".resolve_species falls back to positional labels", {
  m <- list(list(data_list = list(nspp = 2L, spnames = NULL)))
  testthat::expect_equal(.resolve_species(m)$spnames,
                         c("Species 1", "Species 2"))
  # A model whose own spnames are the wrong length gets positional ones rather
  # than recycled wrong ones.
  m2 <- list(list(data_list = list(nspp = 3L, spnames = c("a", "b"))))
  testthat::expect_equal(.resolve_species(m2)$spnames,
                         c("Species 1", "Species 2", "Species 3"))
})


# --- .rce_line_params ---------------------------------------------------------
testthat::test_that("the default lwd renders as linewidth 1", {
  # The invariant that keeps every converted plotter's default figure identical.
  lp <- .rce_line_params(lwd = 3)
  testthat::expect_equal(lp$args$linewidth, 1)
  testthat::expect_null(lp$mapping)
  testthat::expect_length(lp$scales, 0L)

  testthat::expect_equal(.rce_line_params(lwd = 6)$args$linewidth, 2)
  testthat::expect_equal(.rce_line_params(lwd = 2)$args$linewidth, 2 / 3)
})

testthat::test_that("a varying lwd/lty is mapped when there is a key for it", {
  lp <- .rce_line_params(lwd = c(3, 6), lwd_by = "Model")
  testthat::expect_false(is.null(lp$mapping))
  testthat::expect_length(lp$scales, 1L)
  testthat::expect_null(lp$args$linewidth)   # mapped, not fixed

  lp2 <- .rce_line_params(lty = c(1, 2), lty_by = "Sex")
  testthat::expect_false(is.null(lp2$mapping))
  testthat::expect_null(lp2$args$linetype)
})

testthat::test_that("a mapped lwd applies lwd/3 in level order, like the fixed one", {
  # The mapped branch must use the same base-graphics bridge as the fixed
  # branch, and must not depend on the caller knowing the factor's levels: a
  # named manual scale whose names miss the data renders every line blank.
  p <- ggplot2::ggplot(demo_df, ggplot2::aes(x = .data$x, y = .data$y,
                                             colour = .data$Model))
  b <- ggplot2::ggplot_build(
    .rce_add_line(p, .rce_line_params(lwd = c(3, 6), lwd_by = "Model")))
  testthat::expect_equal(sort(unique(b$data[[1]]$linewidth)), c(1, 2))
  testthat::expect_false(anyNA(b$data[[1]]$linewidth))
})

testthat::test_that("lty reaches a plot that maps line type in its own aes()", {
  # plot_ration() and plot_m_at_age() map line type to sex in their own aes().
  # A fixed geom parameter would override that mapping and drop its legend, so
  # the value goes through a scale instead.

  # The default is what the first level already draws: leave the mapping alone,
  # so no default figure moves.
  lp <- .rce_line_params(lty = 1, lty_by = "Sex", lty_in_aes = TRUE)
  testthat::expect_null(lp$args$linetype)
  testthat::expect_length(lp$scales, 0L)

  # A vector overrides, one value per level.
  lp2 <- .rce_line_params(lty = c(1, 3), lty_by = "Sex", lty_in_aes = TRUE)
  testthat::expect_length(lp2$scales, 1L)
  testthat::expect_null(lp2$args$linetype)

  # A single non-default value applies to EVERY level. Dropping it here is what
  # made `plot_ration(fit, lty = 2)` silently do nothing on a sex-combined fit.
  lp3 <- .rce_line_params(lty = 2, lty_by = "Sex", lty_n = 1L, lty_in_aes = TRUE)
  testthat::expect_null(lp3$args$linetype)
  testthat::expect_length(lp3$scales, 1L)
  testthat::expect_equal(lp3$scales[[1]]$palette(3), rep(2, 3))

  # ... and says so when the key it collapses does separate something.
  testthat::expect_warning(
    .rce_line_params(lty = 2, lty_by = "Sex", lty_n = 2L, lty_in_aes = TRUE),
    "drawn alike")

  # Where the plot does NOT map line type, a scalar is still applied.
  testthat::expect_equal(.rce_line_params(lty = 2, lty_by = "Model")$args$linetype, 2)
})

testthat::test_that(".rce_line_params rejects widths that would not draw", {
  # NA silently removes the line; a negative width errors deep inside the device.
  testthat::expect_error(.rce_line_params(lwd = NA), "non-negative")
  testthat::expect_error(.rce_line_params(lwd = -3), "non-negative")
  testthat::expect_error(.rce_line_params(lwd = "3"), "non-negative")
  testthat::expect_error(.rce_line_params(lty = NA), "must not be NA")
  # A character line type is legitimate.
  testthat::expect_equal(.rce_line_params(lty = "dashed")$args$linetype, "dashed")
})

testthat::test_that("a varying lwd/lty without a key warns and uses the first", {
  testthat::expect_warning(lp <- .rce_line_params(lwd = c(3, 6)),
                           "draws one line width")
  testthat::expect_equal(lp$args$linewidth, 1)

  testthat::expect_warning(lp2 <- .rce_line_params(lty = c(1, 2)),
                           "draws one line type")
  testthat::expect_equal(lp2$args$linetype, 1)
})

testthat::test_that("alpha is passed through only when supplied", {
  testthat::expect_null(.rce_line_params()$args$alpha)
  testthat::expect_equal(.rce_line_params(alpha = 0.25)$args$alpha, 0.25)
})


# --- .rce_check_alpha ---------------------------------------------------------
testthat::test_that("a transparency outside [0, 1] stops, naming the argument", {
  testthat::expect_error(.rce_check_alpha(1.5), "between 0 and 1")
  testthat::expect_error(.rce_check_alpha(-0.5), "between 0 and 1")
  testthat::expect_error(.rce_check_alpha(NA), "between 0 and 1")
  # A vector would silently apply only its first element.
  testthat::expect_error(.rce_check_alpha(c(0.2, 0.5)), "between 0 and 1")
  testthat::expect_error(.rce_check_alpha("0.4"), "between 0 and 1")
})

testthat::test_that("a valid alpha, and NULL, pass through unchanged", {
  testthat::expect_silent(.rce_check_alpha(0.4))
  testthat::expect_silent(.rce_check_alpha(0))
  testthat::expect_silent(.rce_check_alpha(1))
  # NULL is the "this figure has no ribbon to shade" case, not an error.
  testthat::expect_silent(.rce_check_alpha(NULL))
})


# --- .rce_age_index / .rce_age_plus_index -------------------------------------
# Ages run minage .. minage + nages - 1; the arrays are indexed 1 .. nages. The
# two coincide only at minage = 1, which is every bundled dataset, so these are
# the tests that hold the distinction.
age_model <- list(minage = c(1L, 2L, 0L), nages = c(12L, 10L, 5L),
                  spnames = c("Pollock", "Cod", "ATF"))

testthat::test_that("an age resolves to its bin on each species' own scale", {
  r <- with(age_model, .rce_age_index(3, 1:3, minage, nages, spnames))
  testthat::expect_equal(r$species, 1:3)
  # age 3: bin 3 at minage 1, bin 2 at minage 2, bin 4 at minage 0.
  testthat::expect_equal(r$index, c(3L, 2L, 4L))
})

testthat::test_that("a species without the age is dropped, with a warning", {
  testthat::expect_warning(
    r <- with(age_model, .rce_age_index(1, 1:3, minage, nages, spnames)),
    "Cod \\(2-11\\)")
  testthat::expect_equal(r$species, c(1L, 3L))
  testthat::expect_equal(r$index, c(1L, 2L))
})

testthat::test_that("an age no species carries is an error, not an empty plot", {
  testthat::expect_error(
    with(age_model, .rce_age_index(40, 1:3, minage, nages, spnames)),
    "No species has age 40")
  testthat::expect_error(
    with(age_model, .rce_age_index(c(1, 2), 1:3, minage, nages, spnames)),
    "single age")
})

testthat::test_that("`minage`+ keeps the bins at or above that age", {
  r <- with(age_model, .rce_age_plus_index(3, 1:3, minage, nages, spnames))
  # Pollock ages 1-12 -> bins 3:12; Cod 2-11 -> bins 2:10; ATF 0-4 -> bins 4:5.
  testthat::expect_equal(r$index[[1]], 3:12)
  testthat::expect_equal(r$index[[2]], 2:10)
  testthat::expect_equal(r$index[[3]], 4:5)

  # A species with nothing that old drops out; all of them is an error.
  testthat::expect_warning(
    r2 <- with(age_model, .rce_age_plus_index(11, 1:3, minage, nages, spnames)),
    "ATF")
  testthat::expect_equal(r2$species, c(1L, 2L))
  testthat::expect_error(
    with(age_model, .rce_age_plus_index(99, 1:3, minage, nages, spnames)),
    "No species has an age at or above")
})

testthat::test_that(".rce_add_line reaches the built plot", {
  p <- ggplot2::ggplot(demo_df, ggplot2::aes(x = .data$x, y = .data$y,
                                             colour = .data$Model))
  b <- ggplot2::ggplot_build(.rce_add_line(p, .rce_line_params(lwd = 6, lty = 2)))
  testthat::expect_equal(unique(b$data[[1]]$linewidth), 2)
  testthat::expect_equal(unique(b$data[[1]]$linetype), 2)
})


# --- .rceattle_scale ----------------------------------------------------------
testthat::test_that("line_col = NULL keeps the package default palettes", {
  # Guards against the extension changing any existing figure.
  disc <- ggplot2::ggplot(demo_df, ggplot2::aes(.data$x, .data$y,
                                                colour = .data$Model)) +
    ggplot2::geom_line()
  testthat::expect_equal(built_colours(.rceattle_scale(disc, aesthetics = "colour")),
                         .okabe_ito[1:2])

  # Pin viridis by value, not merely "some hex codes" -- the continuous default
  # is exactly what the refactor could have broken.
  cont <- ggplot2::ggplot(demo_df, ggplot2::aes(.data$x, .data$y,
                                                colour = .data$Year,
                                                group = .data$Model)) +
    ggplot2::geom_line()
  testthat::expect_equal(
    built_colours(.rceattle_scale(cont, discrete = FALSE, aesthetics = "colour")),
    built_colours(cont + ggplot2::scale_colour_viridis_c()))
})

testthat::test_that("line_col overrides a discrete colour scale", {
  p <- ggplot2::ggplot(demo_df, ggplot2::aes(.data$x, .data$y,
                                             colour = .data$Model)) +
    ggplot2::geom_line()
  testthat::expect_equal(
    built_colours(.rceattle_scale(p, aesthetics = "colour",
                                  line_col = c("red", "blue"))),
    c("red", "blue"))

  # Integer colours must not error.
  testthat::expect_equal(
    built_colours(.rceattle_scale(p, aesthetics = "colour", line_col = 1)),
    "black")

  # Recycled when fewer colours than series, never grey50.
  testthat::expect_equal(
    built_colours(.rceattle_scale(p, aesthetics = "colour", line_col = "red")),
    "red")

  # An NA level keeps the grey the named scales used to give it; without an
  # explicit na.value it would draw nothing at all.
  na_df <- data.frame(x = 1:4, y = 1:4,
                      Model = factor(c("A", "A", NA, NA)))
  na_p <- ggplot2::ggplot(na_df, ggplot2::aes(.data$x, .data$y,
                                              colour = .data$Model)) +
    ggplot2::geom_line()
  testthat::expect_true("grey50" %in%
    built_colours(.rceattle_scale(na_p, aesthetics = "colour",
                                  line_col = c("red", "blue"))))
})

testthat::test_that("line_col reaches the fill aesthetic too", {
  # A ribbon and its line must not disagree, which is the reason the two
  # aesthetics share one code path.
  df <- data.frame(x = 1:4, y = 1:4, lo = 0:3, hi = 2:5,
                   Model = factor(rep(c("a", "b"), each = 2)))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y,
                                        colour = .data$Model,
                                        fill = .data$Model)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lo, ymax = .data$hi))
  b <- ggplot2::ggplot_build(
    .rceattle_scale(p, line_col = c("red", "blue")))
  testthat::expect_equal(unique(b$data[[1]]$fill), c("red", "blue"))

  # Continuous fill must dispatch to a fill scale, not a colour one.
  pc <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y,
                                         fill = .data$x)) +
    ggplot2::geom_raster()
  testthat::expect_equal(
    unique(ggplot2::ggplot_build(
      .rceattle_scale(pc, discrete = FALSE, aesthetics = "fill",
                      line_col = "black"))$data[[1]]$fill),
    "#000000")
})

testthat::test_that("line_col supplies ramp anchors on a continuous colour scale", {
  p <- ggplot2::ggplot(demo_df, ggplot2::aes(.data$x, .data$y,
                                             colour = .data$Year,
                                             group = .data$Model)) +
    ggplot2::geom_line()
  # One colour: the whole year fan drawn in it -- "line_col = 'black'" meaning
  # what a base-graphics user expects it to mean.
  testthat::expect_equal(
    built_colours(.rceattle_scale(p, discrete = FALSE, aesthetics = "colour",
                                  line_col = "black")),
    "#000000")
  # Two or more: interpolated between them, so the fan still reads as a ramp.
  many <- built_colours(.rceattle_scale(p, discrete = FALSE,
                                        aesthetics = "colour",
                                        line_col = c("red", "blue")))
  testthat::expect_gt(length(many), 1L)
  testthat::expect_true("#FF0000" %in% many)
})


# --- year window --------------------------------------------------------------
testthat::test_that(".rce_year_filter drops rows outside the window", {
  df <- data.frame(Year = 1990:1999, v = 1:10)
  testthat::expect_equal(nrow(.rce_year_filter(df, 1995, 1997)), 3L)
  testthat::expect_equal(nrow(.rce_year_filter(df, minyr = 1998)), 2L)
  testthat::expect_equal(nrow(.rce_year_filter(df)), 10L)
  testthat::expect_error(.rce_year_filter(df, 2050), "No years left")
})

testthat::test_that(".rce_year_filter drops NA years instead of fabricating a row", {
  # A comparison against an NA year gives NA, which subsets to an all-NA row.
  df <- data.frame(Year = c(1990, NA, 1992, 1995), v = 1:4)
  out <- .rce_year_filter(df, 1990, 1993)
  testthat::expect_equal(nrow(out), 2L)
  testthat::expect_false(anyNA(out$Year))

  testthat::expect_error(.rce_year_filter(df, NA), "single year")
  testthat::expect_error(.rce_year_filter(df, maxyr = c(1, 2)), "single year")
})


# --- projection divider / mean line -------------------------------------------
testthat::test_that(".rce_proj_divider draws at the latest hindcast year", {
  m <- list(list(data_list = list(endyr = 2017)))
  testthat::expect_null(.rce_proj_divider(m, incl_proj = FALSE))

  lyr <- .rce_proj_divider(m, incl_proj = TRUE)
  # Check what it draws, not merely that it is a layer -- geom_point() would
  # satisfy an expect_s3_class(., "ggproto").
  testthat::expect_equal(lyr$data$xintercept, 2017)
  testthat::expect_equal(lyr$aes_params$linetype, 2)
  testthat::expect_equal(lyr$aes_params$colour, "grey50")

  # On a retrospective peel the models end in different years; keying on
  # whichever is last would put the divider mid-hindcast.
  peels <- list(list(data_list = list(endyr = 2017)),
                list(data_list = list(endyr = 2015)))
  testthat::expect_equal(.rce_proj_divider(peels, TRUE)$data$xintercept, 2017)

  # A model without an endyr must not produce a vline at NA.
  testthat::expect_null(
    .rce_proj_divider(list(list(data_list = list(endyr = NULL))), TRUE))
})

testthat::test_that(".rce_proj_divider stays out of a window it falls outside", {
  # A vline trains the x scale, so a divider at 2017 on a minyr = 2020 window
  # would stretch the axis back to 2017 -- the empty span the window removed.
  m <- list(list(data_list = list(endyr = 2017)))
  testthat::expect_null(.rce_proj_divider(m, TRUE, minyr = 2020))
  testthat::expect_null(.rce_proj_divider(m, TRUE, maxyr = 2010))
  testthat::expect_equal(
    .rce_proj_divider(m, TRUE, minyr = 2000, maxyr = 2050)$data$xintercept, 2017)
})

testthat::test_that(".rce_mean_line averages hindcast years only", {
  df <- data.frame(Species = "a", Model = "m", Year = 1:10,
                   value = c(rep(1, 5), rep(100, 5)))
  testthat::expect_null(.rce_mean_line(df, incl_mean = FALSE, hind_endyr = 5))

  lyr <- .rce_mean_line(df, incl_mean = TRUE, by = c("Species", "Model"),
                        hind_endyr = 5)
  # Mean over years 1:5 only; the projection's 100s must not drag it up.
  testthat::expect_equal(lyr$data$.mean, 1)

  # hind_endyr is required, so a caller cannot silently get the
  # projection-contaminated mean the helper exists to avoid.
  testthat::expect_error(.rce_mean_line(df, incl_mean = TRUE, by = "Species"),
                         "hind_endyr` is required")
  # NA opts into averaging everything, deliberately.
  testthat::expect_equal(
    .rce_mean_line(df, TRUE, by = c("Species", "Model"),
                   hind_endyr = NA)$data$.mean,
    mean(df$value))
})

testthat::test_that(".rce_mean_line takes a colour only when colour encodes it", {
  # Mapping colour = Model on a plot whose colour encodes species trains model
  # names into the species legend.
  df <- data.frame(Species = rep(c("s1", "s2"), each = 4),
                   Model = rep(c("m1", "m2"), 4),
                   Year = rep(1:4, 2), value = 1:8)

  neutral <- .rce_mean_line(df, TRUE, by = c("Species", "Model"), hind_endyr = 4)
  testthat::expect_false("colour" %in% names(neutral$mapping))
  testthat::expect_equal(neutral$aes_params$colour, "grey40")

  keyed <- .rce_mean_line(df, TRUE, by = c("Species", "Model"), hind_endyr = 4,
                          colour_by = "Model")
  testthat::expect_true("colour" %in% names(keyed$mapping))

  # A colour_by column the aggregate does not carry falls back to neutral
  # rather than erroring at render time.
  testthat::expect_false(
    "colour" %in% names(.rce_mean_line(df, TRUE, by = "Species", hind_endyr = 4,
                                       colour_by = "Model")$mapping))
})


# --- confidence-interval plumbing ---------------------------------------------
testthat::test_that(".rce_series_sd requires an exact length match", {
  testthat::expect_null(.rce_series_sd(list(sdrep = NULL), "ssb", 5))
  fit <- list(sdrep = list(value = stats::setNames(1:3, rep("ssb", 3)),
                           sd = c(0.1, 0.2, 0.3)))
  testthat::expect_null(.rce_series_sd(fit, "ssb", 5))          # too few
  testthat::expect_equal(.rce_series_sd(fit, "ssb", 3), c(0.1, 0.2, 0.3))
  testthat::expect_null(.rce_series_sd(fit, "biomass", 1))      # wrong name

  # A peel keeps an sdreport whose rebuild failed, but it is not a band.
  attr(fit, "sd_report_error") <- "peel 1 (through 2015): indefinite"
  testthat::expect_null(.rce_series_sd(fit, "ssb", 3))

  # Taking the first n of a longer block would return another cell's standard
  # errors -- a confidence interval that is wrong rather than absent.
  two_spp <- list(sdrep = list(value = stats::setNames(1:10, rep("ssb", 10)),
                               sd = c(101:105, 201:205)))
  testthat::expect_null(.rce_series_sd(two_spp, "ssb", 5))
  testthat::expect_equal(.rce_series_sd(two_spp, "ssb", 10), c(101:105, 201:205))
})

testthat::test_that(".rce_no_ci warns only when the interval was asked for", {
  testthat::expect_silent(res <- .rce_no_ci(FALSE, "ration", "no standard errors"))
  testthat::expect_false(res)
  testthat::expect_warning(res2 <- .rce_no_ci(TRUE, "ration", "no standard errors"),
                           "not available for ration")
  testthat::expect_false(res2)

  # A REPORT()-only quantity was never going to carry standard errors, so
  # warning about it would fire on every fit.
  testthat::expect_silent(.rce_no_ci(TRUE, "F_spp", "REPORT only", warn = FALSE))
})

testthat::test_that(".rce_quantity_adreport matches the shared registry", {
  # The figure and as.data.frame.Rceattle() must agree on which series have
  # standard errors, so both read one registry.
  testthat::expect_true(.rce_quantity_adreport("ssb"))
  testthat::expect_true(.rce_quantity_adreport("biomass"))
  # M_at_age and B_eaten_as_prey have their ADREPORT commented out in
  # ceattle.cpp; they are REPORT-only and carry no standard errors.
  testthat::expect_false(.rce_quantity_adreport("M_at_age"))
  testthat::expect_false(.rce_quantity_adreport("B_eaten_as_prey"))
  testthat::expect_false(.rce_quantity_adreport("F_spp"))
  # A derived combination is absent from the registry entirely.
  testthat::expect_false(.rce_quantity_adreport("ration"))
})
