# plot F

Fishing mortality over time, one panel per species, with the Ftarget
(blue) and Flimit (red) reference points drawn per panel.

## Usage

``` r
plot_f(
  Rceattle,
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
  ylab = NULL
)
```

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

- line_col:

  Line colours; names, hex, or base-graphics integers. `NULL` uses the
  colorblind-safe Okabe-Ito palette. Applied to whichever variable the
  figure separates by colour, in legend order. Too few colours are
  recycled, with a warning.

- species:

  Species to include, as indices (`c(1, 3)`), names, a logical mask, or
  `"all"`. Default `NULL` plots every species. Species **labels** belong
  in `spnames`; a character `species` that matches no species name is
  read as labels for back-compatibility.

- spnames:

  Species labels, length `nspp`. Default: the model's own.

- add_ci:

  Add a 95% confidence interval. Only available where the plotted
  quantity holds standard errors; warns and draws none otherwise.

- lwd:

  Line width on the base-graphics scale: the default `3` renders as a
  standard-weight ggplot line. A vector varies it across series.

- save:

  Write the plotted series to CSV alongside the figure.

- right_adj:

  Ignored. Base-graphics leftover: the figure widened its right margin
  to fit the legend. Set margins on the returned ggplot instead.

- legend.pos:

  Ignored. Base-graphics leftover: legend placement. Use
  `p + ggplot2::theme(legend.position = "right")`.

- width, height:

  Saved figure size in inches.

- minyr, maxyr:

  First / last year to plot.

- incl_proj:

  Include the projection years, with a dashed divider at the last
  hindcast year.

- mod_cex:

  Ignored. Base-graphics leftover: legend text size. Use
  `p + ggplot2::theme(legend.text = ggplot2::element_text(size = ...))`.

- lty:

  Line type. A vector varies it across the levels of whatever the figure
  separates by line type.

- alpha:

  Transparency of confidence ribbons and shaded areas, between 0 and 1.

- mod_avg:

  TRUE/FALSE

- mse:

  Is an MSE object from
  [`load_mse`](https://afsc-assessments.github.io/Rceattle/reference/load_mse.md)
  or
  [`run_mse`](https://afsc-assessments.github.io/Rceattle/reference/run_mse.md)

- OM:

  if mse == TRUE, use the OM (TRUE) or EM (FALSE) for plotting?

- reference:

  A model to overlay for comparison, drawn at 1.5x `lwd` and labelled
  "Reference". It takes the next colour from the palette, or black if
  `line_col` is supplied.

- ylab:

  Y-axis label. `NULL` (default) derives one from `output` and the
  model's `minage`.

## Value

Returns and saves a figure with the fishing mortality trajectory.
