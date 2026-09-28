# Plot time-series

Function that plots the time-series (SSB/B/R/Depletion) 95% CI trends as
estimated from Rceattle

## Usage

``` r
plot_timeseries(
  Rceattle,
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
  mod_avg = rep(FALSE, length(Rceattle))
)
```

## Arguments

- Rceattle:

  A single
  [`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)
  object or a list of them (overlaid).

- output:

  derived quantity of interest: recruitment, biomass, ssb, depletion, or
  ssb_depletion. Uses same name as ".cpp" file.

- ylab:

  Y-axis label. `NULL` (default) derives one from `output` and the
  model's `minage`.

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

- legend.pos:

  Ignored. Base-graphics leftover: legend placement. Use
  `p + ggplot2::theme(legend.position = "right")`.

- right_adj:

  Ignored. Base-graphics leftover: the figure widened its right margin
  to fit the legend. Set margins on the returned ggplot instead.

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

- zero_y:

  Anchor the y-axis at zero. TRUE for the absolute series (biomass, SSB,
  recruitment), where a truncated axis exaggerates change; FALSE for the
  depletions, which are already on a relative scale.

- ref_lines:

  Internal. `function(models, sp_sel)` returning per-species
  reference-point layers, added before the figure is saved. Supplied by
  [`plot_f()`](https://afsc-assessments.github.io/Rceattle/reference/plot_f.md)
  and
  [`plot_depletionSSB()`](https://afsc-assessments.github.io/Rceattle/reference/plot_depletionSSB.md)
  through `.ts_wrapper()`; `sp_sel` is the species resolution the panels
  were built from.

- suffix:

  Internal. Overrides the `<file>_<suffix>.png` stem, which otherwise
  names the series.

- mod_avg:

  TRUE/FALSE

## Value

Returns and saves a figure with the population trajectory.

## Units

The model holds numbers-at-age in **thousands** and weight-at-age in
**kg**, so every biomass series (`biomass`, `ssb`,
`exploitable_biomass`) comes out of the model in **mt** and recruitment
comes out in **thousands of fish**. For display these are divided by 1e6
(million mt) and 1e3 (millions of recruits) respectively; depletion is a
ratio and is not rescaled. Supply the model's inputs on that convention,
catch and index in mt, weight-at-age in kg, or the axis labels will not
describe what is plotted.

## Confidence intervals

`add_ci = TRUE` draws a 95% interval. Every strictly positive series
takes it on the log scale as `exp(log(x) +/- qnorm(0.975) * sd_log)`.
These quantities are built multiplicatively, so `log(x)` is close to
linear in the estimated parameters where `x` is exponential in them: the
interval is both a better linearization and right-skewed, and it cannot
cross zero the way a symmetric natural-scale interval does for weak year
classes and depleted stocks.

`sd_log` comes from the model's own `log_biomass` / `log_ssb` / `log_R`
where those are reported, and is otherwise recovered as `sd(x) / x`, the
delta method's own identity, which matches the reported values to
machine precision. That covers `exploitable_biomass` and the two
depletions, which cannot be reported on the log scale
(`exploitable_biomass` is identically 0 without projection F), and
models fit before the `log_*` series existed. A non-positive or
unreported value keeps the symmetric interval.
