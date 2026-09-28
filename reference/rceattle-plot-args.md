# Shared plotting arguments

A common vocabulary, documented here once and inherited by the plotters
that take it:
[`plot_timeseries()`](https://afsc-assessments.github.io/Rceattle/reference/plot_timeseries.md)
and its wrappers
([`plot_biomass()`](https://afsc-assessments.github.io/Rceattle/reference/plot_biomass.md),
[`plot_ssb()`](https://afsc-assessments.github.io/Rceattle/reference/plot_ssb.md),
[`plot_recruitment()`](https://afsc-assessments.github.io/Rceattle/reference/plot_recruitment.md),
the depletions,
[`plot_exploitable_biomass()`](https://afsc-assessments.github.io/Rceattle/reference/plot_exploitable_biomass.md),
[`plot_f()`](https://afsc-assessments.github.io/Rceattle/reference/plot_f.md)),
the predation plotters
([`plot_b_eaten()`](https://afsc-assessments.github.io/Rceattle/reference/plot_b_eaten.md),
[`plot_b_eaten_prop()`](https://afsc-assessments.github.io/Rceattle/reference/plot_b_eaten_prop.md),
[`plot_m_at_age()`](https://afsc-assessments.github.io/Rceattle/reference/plot_m_at_age.md),
[`plot_m2_at_age_prop()`](https://afsc-assessments.github.io/Rceattle/reference/plot_m2_at_age_prop.md),
[`plot_ration()`](https://afsc-assessments.github.io/Rceattle/reference/plot_ration.md)),
and
[`plot_selectivity()`](https://afsc-assessments.github.io/Rceattle/reference/plot_selectivity.md).
Each argument means the same thing wherever it appears, but not every
plotter takes every one, `incl_mean` is on the predation plotters,
`add_ci` only where the quantity holds standard errors, and `alpha` only
where the figure has a ribbon or a fan. The remaining `plot_*()`
functions still take their own arguments; see each one's help.

## Arguments

- Rceattle:

  A single
  [`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)
  object or a list of them (overlaid).

- file:

  Optional file stem; the figure is written to `<file>_<suffix>.png` if
  given.

- model_names:

  Legend labels for the models.

- species:

  Species to include, as indices (`c(1, 3)`), names, a logical mask, or
  `"all"`. Default `NULL` plots every species. Species **labels** belong
  in `spnames`; a character `species` that matches no species name is
  read as labels for back-compatibility.

- spnames:

  Species labels, length `nspp`. Default: the model's own.

- line_col:

  Line colours; names, hex, or base-graphics integers. `NULL` uses the
  colorblind-safe Okabe-Ito palette. Applied to whichever variable the
  figure separates by colour, in legend order. Too few colours are
  recycled, with a warning.

- lwd:

  Line width on the base-graphics scale: the default `3` renders as a
  standard-weight ggplot line. A vector varies it across series.

- lty:

  Line type. A vector varies it across the levels of whatever the figure
  separates by line type.

- alpha:

  Transparency of confidence ribbons and shaded areas, between 0 and 1.

- add_ci:

  Add a 95% confidence interval. Only available where the plotted
  quantity holds standard errors; warns and draws none otherwise.

- minyr, maxyr:

  First / last year to plot.

- incl_proj:

  Include the projection years, with a dashed divider at the last
  hindcast year.

- incl_mean:

  Add a horizontal line at the hindcast mean of each series.

- width, height:

  Saved figure size in inches.

- right_adj:

  Ignored. Base-graphics leftover: the figure widened its right margin
  to fit the legend. Set margins on the returned ggplot instead.

- top_adj:

  Ignored. Base-graphics leftover; see `right_adj`.

- mod_cex:

  Ignored. Base-graphics leftover: legend text size. Use
  `p + ggplot2::theme(legend.text = ggplot2::element_text(size = ...))`.

- legend.pos:

  Ignored. Base-graphics leftover: legend placement. Use
  `p + ggplot2::theme(legend.position = "right")`.

- single.plots:

  Ignored. Base-graphics leftover: one device per panel. The ggplot
  facets instead.

- theta:

  Ignored. Base-graphics leftover: viewing angle of a `persp` surface,
  which is no longer drawn.

- ymax:

  Ignored. Base-graphics leftover: y-axis maximum. Use
  `p + ggplot2::coord_cartesian(ylim = c(0, ymax))`.

- cex:

  Ignored. Base-graphics leftover: point expansion.

## Value

A `ggplot` object.

## How `line_col` and `lty` are applied

They supply values for whichever **discrete variable the plot already
encodes with that aesthetic**, matched in level order, which is not
always the model. Each function's help says what its own figure
separates.

Where colour encodes a continuous variable, the year fan in
[`plot_selectivity()`](https://afsc-assessments.github.io/Rceattle/reference/plot_selectivity.md),
`line_col` supplies the ramp anchors instead: one colour draws the fan
in that colour, several interpolate between them.
