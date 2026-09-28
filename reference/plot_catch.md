# Fishery catch fits

Plots fitted fishery catch: observed points with lognormal 95% intervals
and the model-predicted catch, faceted by fishery fleet. For MSE objects
the projection period is summarized with 50% / 95% ribbons across
simulations.

## Usage

``` r
plot_catch(
  Rceattle,
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
  ymax = NULL
)
```

## Arguments

- Rceattle:

  A single
  [`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)
  object or a list of them (overlaid).

- file:

  Optional file stem; if given the figure is written to
  `<file>_survey_indices.png`.

- model_names:

  Legend labels for the models.

- species:

  Species (indices) to include. Default all.

- incl_proj:

  Include projection years.

- width, height:

  Saved figure size (inches).

- error:

  Draw observed points and error bars.

- maxyr:

  Last year to plot.

- mse:

  Is `Rceattle` an MSE object (operating models)?

- fleets:

  Fishery fleets to include (indices into the sorted fleet codes).
  Default all.

- line_col, right_adj, top_adj, single.plots:

  Deprecated base-graphics arguments, retained for back-compatibility
  and ignored.

- alpha, lwd, ymax:

  Deprecated base-graphics arguments, ignored.

## Value

A `ggplot` object.
