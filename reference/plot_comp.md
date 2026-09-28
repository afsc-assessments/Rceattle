# Plot composition fits and residuals

Diagnostic plots for age / length composition data from a fitted
[Rceattle](https://afsc-assessments.github.io/Rceattle/reference/Rceattle-package.md)
model, drawn with ggplot2 for a consistent look with
[`plot.rceattle_osa()`](https://afsc-assessments.github.io/Rceattle/reference/plot.rceattle_osa.md).
With the default `residual_type = "pearson"` three figures are produced
per fleet x composition type:

- **Pearson residual bubbles** by year and bin, faceted by fleet (and,
  for joint-sex data, by sex); red = positive, blue = negative, sized by
  magnitude. The Pearson residual is the same form used by
  [`residuals.Rceattle()`](https://afsc-assessments.github.io/Rceattle/reference/residuals.Rceattle.md),
  standardized by the variance the fleet's own likelihood assumes, so
  the weight column and any Dirichlet-multinomial overdispersion are
  already in the denominator.

- **Annual composition**: observed (shaded area) vs fitted (line)
  proportion at age / length, one panel per year. Joint-sex data are
  mirrored (females up, males down).

- **Aggregated composition**: observed vs fitted, counts pooled over the
  fitted hindcast years. The interval holds 95% of the data the model
  predicts (not a confidence interval on the mean), so observations
  outside it indicate misfit. Poor where the expected count is below
  about 10, and slightly narrow throughout, because it treats the fitted
  proportions as known rather than estimated from these same data.

The shaded area and fitted line span only the observed bins (they do not
extend past the first/last bin), and bins with zero observed proportion
are retained (only `NA` bins are dropped), so the curves are not
interpolated across gaps.

## Usage

``` r
plot_comp(
  Rceattle,
  file = NULL,
  model_names = NULL,
  species = NULL,
  cex = 3,
  lwd = 3,
  right_adj = 0,
  residual_type = c("pearson", "osa"),
  add_agg_ci = TRUE,
  add_agg_n = TRUE
)
```

## Arguments

- Rceattle:

  A single fitted `Rceattle` model.

- file:

  Optional filename stem; when supplied, each figure is also saved as a
  PNG.

- model_names:

  Unused (kept for back-compatibility).

- species:

  Optional species code(s) to plot (matched against the `Species`
  column); `NULL` (default) plots all. Mirrors the `species` argument of
  [`residuals.Rceattle()`](https://afsc-assessments.github.io/Rceattle/reference/residuals.Rceattle.md)
  and
  [`plot.rceattle_osa()`](https://afsc-assessments.github.io/Rceattle/reference/plot.rceattle_osa.md).

- cex:

  Unused (kept for back-compatibility).

- lwd:

  Width of the fitted-composition line. Default `3`.

- right_adj:

  Unused (kept for back-compatibility).

- residual_type:

  `"pearson"` (default) for the ggplot2 Pearson-residual and
  composition-fit figures drawn here, or `"osa"` to instead draw the
  one-step-ahead residual diagnostics via
  [`osa_residuals()`](https://afsc-assessments.github.io/Rceattle/reference/osa_residuals.md)
  and
  [`plot.rceattle_osa()`](https://afsc-assessments.github.io/Rceattle/reference/plot.rceattle_osa.md),
  a Q-Q plot (with SDNR / tail annotation) alongside signed OSA- and
  Pearson-residual bubbles. The `"osa"` path builds its observation data
  on demand, so it works with any fit.

- add_agg_ci:

  Logical. Draw the 95% prediction interval on the aggregated
  composition figure. Default `TRUE`.

- add_agg_n:

  Logical. Annotate the aggregated composition figure with the input
  sample size, the effective sample size the likelihood assumed, and the
  McAllister-Ianelli effective sample size this fit implies. Default
  `TRUE`.

## Value

Invisibly, a named list of the `ggplot` objects. Called for its side
effect of drawing (and optionally saving) the figures.
