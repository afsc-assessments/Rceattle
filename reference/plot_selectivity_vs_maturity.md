# Plot fishery selectivity and maturity

Overlays terminal-year fishery selectivity-at-age against input
maturity-at-age, faceted by species. Useful for debugging SPR-based
reference points (the two curves drive spawning-potential ratio).

## Usage

``` r
plot_selectivity_vs_maturity(
  Rceattle,
  file = NULL,
  model_names = NULL,
  line_col = NULL,
  width = 7,
  height = 6.5,
  species = NULL,
  lwd = 3
)
```

## Arguments

- Rceattle:

  A single
  [`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)
  object or a list of them (the first is used).

- file:

  Optional file stem; the figure is written to
  `<file>_selectivity_vs_maturity.png` if given.

- model_names, line_col, lwd:

  Deprecated base-graphics arguments, retained for back-compatibility
  and ignored.

- width, height:

  Saved figure size (inches).

- species:

  Species names for the facet labels. Default: model species names.

## Value

A `ggplot` object.
