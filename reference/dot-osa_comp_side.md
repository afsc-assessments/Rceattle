# One bin-side (age or length) composition column: Q-Q + OSA + Pearson bubbles

One bin-side (age or length) composition column: Q-Q + OSA + Pearson
bubbles

## Usage

``` r
.osa_comp_side(comp, pear, side, add_sdnr_ci = TRUE, add_qq_quantiles = TRUE)
```

## Arguments

- comp:

  Composition rows with `source` and `.side` columns.

- pear:

  Reshaped Pearson rows with `source` and `.side` columns, or `NULL`.

- side:

  `"age"` or `"length"`.

- add_sdnr_ci, add_qq_quantiles:

  Passed to
  [`.osa_qqplot()`](https://afsc-assessments.github.io/Rceattle/reference/dot-osa_qqplot.md).

## Value

A stacked `cowplot`/`ggplot` object, or `NULL` if no rows on that side.
