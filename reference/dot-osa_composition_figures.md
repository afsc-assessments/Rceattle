# Composition OSA figure(s): Q-Q + OSA bubbles + Pearson bubbles (age / length)

Composition OSA figure(s): Q-Q + OSA bubbles + Pearson bubbles (age /
length)

## Usage

``` r
.osa_composition_figures(
  comp,
  pearson = NULL,
  combine = TRUE,
  nages = NULL,
  nlengths = NULL,
  add_sdnr_ci = TRUE,
  add_qq_quantiles = TRUE
)
```

## Arguments

- comp:

  The composition rows of an `rceattle_osa` object.

- pearson:

  Matching Pearson residuals (the `"pearson"` attribute of the
  `rceattle_osa` object), or `NULL`.

- combine:

  When `TRUE`, return a single `composition` figure with age bins in the
  left column and length bins in the right; when `FALSE`, return
  separate `composition_age` and `composition_length` figures.

- nages, nlengths:

  Per-species bin counts (the `rceattle_osa` `"nages"` / `"nlengths"`
  attributes), used to split joint-sex (Sex == 3) bins onto a single
  age/length axis, matching
  [`plot_comp()`](https://afsc-assessments.github.io/Rceattle/reference/plot_comp.md).

- add_sdnr_ci, add_qq_quantiles:

  Passed to
  [`.osa_qqplot()`](https://afsc-assessments.github.io/Rceattle/reference/dot-osa_qqplot.md).

## Value

A named list of `cowplot`/`ggplot` objects.
