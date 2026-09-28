# Plot stock recruit function

Plots the stock-recruit relationship estimated by Rceattle: spawning
stock biomass (x) against recruitment (y) as points, with the fitted
stock-recruit curve (mean recruitment, Beverton-Holt, or Ricker)
overlaid, faceted by species. A 95% normal data ellipse of the SSB-R
cloud is added when `add_ci = TRUE`.

## Usage

``` r
plot_stock_recruit(
  Rceattle,
  file = NULL,
  model_names = NULL,
  line_col = NULL,
  width = 7,
  height = 6.5,
  species = NULL,
  spnames = NULL,
  lwd = 3,
  lty = 1,
  incl_proj = FALSE,
  plot_env = FALSE,
  mod_cex = 1,
  add_ci = TRUE
)
```

## Arguments

- Rceattle:

  A single
  [`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)
  object or a list of them (overlaid).

- file:

  Optional file stem; the figure is written to
  `<file>_stock_recruit.png` if given.

- model_names:

  Legend labels for the models.

- line_col, lwd, lty, plot_env, mod_cex:

  Deprecated base-graphics arguments, retained for back-compatibility
  and ignored.

- width, height:

  Saved figure size (inches).

- species:

  Species (indices) to include. Default all.

- spnames:

  Species names for the facet labels.

- incl_proj:

  Currently unused (kept for back-compatibility).

- add_ci:

  Add a 95% normal data ellipse of the SSB-R points.

## Value

A `ggplot` object.
