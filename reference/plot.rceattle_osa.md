# Plot one-step-ahead (OSA) residual diagnostics

Diagnostic plots for an `rceattle_osa` object (from
[`osa_residuals()`](https://afsc-assessments.github.io/Rceattle/reference/osa_residuals.md)),
following the recommendations of Stewart and Monnahan (2025) and the
styling of the NOAA-AFSC `afscOSA` package. Under a correctly specified
model OSA residuals are iid standard normal, so the headline diagnostic
is statistical (the Q-Q plot with its annotated SDNR / tail statistics).

Up to two separate figures are drawn, depending on which data are
present:

1.  **Aggregate** (`index` / `catch`): a Q-Q panel faceted by data
    source. These series have no age/length bin, so no bubble plots are
    drawn.

2.  **Composition** (`comp` / `caal`): a Q-Q panel, a signed
    OSA-residual bubble panel, and a signed Pearson-residual bubble
    panel (the Pearson residuals held on the `rceattle_osa` object). By
    default age-based bins (age composition and conditional
    age-at-length) are shown in the left column and length-based bins in
    the right column, each with its own bin axis; set `combine = FALSE`
    to draw the age and length composition as two separate figures
    instead (useful when a model has many fleets/species).

Panel headers use the fleet name from `fleet_control`. Process residuals
(from
[`process_residuals()`](https://afsc-assessments.github.io/Rceattle/reference/process_residuals.md))
are drawn as a Q-Q panel plus a residual-by-year panel.

## Usage

``` r
# S3 method for class 'rceattle_osa'
plot(
  x,
  source = "all",
  species = NULL,
  combine = TRUE,
  add_sdnr_ci = TRUE,
  add_qq_quantiles = TRUE,
  ...
)
```

## Arguments

- x:

  An `rceattle_osa` object from
  [`osa_residuals()`](https://afsc-assessments.github.io/Rceattle/reference/osa_residuals.md)
  or
  [`process_residuals()`](https://afsc-assessments.github.io/Rceattle/reference/process_residuals.md).

- source:

  Data source(s) to plot: any of `"index"`, `"catch"`, `"comp"`,
  `"caal"`, `"diet"`, or `"all"` (default). Mirrors the `source`
  argument of
  [`residuals.Rceattle()`](https://afsc-assessments.github.io/Rceattle/reference/residuals.Rceattle.md);
  filters which figures are produced.

- species:

  Optional species code(s) to include (matched against the `species`
  column). Default `NULL` keeps all species.

- combine:

  Logical. When `TRUE` (default), age and length composition share one
  figure (age in the left column, length in the right). When `FALSE`,
  they are drawn as separate `composition_age` / `composition_length`
  figures.

- add_sdnr_ci:

  Logical. Show the chi-square null interval beside the SDNR annotation
  on each Q-Q panel. Default `TRUE`.

- add_qq_quantiles:

  Logical. Annotate each Q-Q panel with the tail order statistics and
  their exact null intervals. Default `TRUE`.

- ...:

  Unused.

## Value

Invisibly, a named list of the assembled `ggplot` / `cowplot` objects
(some of `aggregate`, `composition` / `composition_age` /
`composition_length`, and `process`). Called for its side effect of
drawing the plot(s).

## References

Stewart, I.J., and Monnahan, C.C. 2025. Can. J. Fish. Aquat. Sci.
82:1-13.

## See also

[`osa_residuals()`](https://afsc-assessments.github.io/Rceattle/reference/osa_residuals.md),
[`osa_diagnostics()`](https://afsc-assessments.github.io/Rceattle/reference/osa_diagnostics.md)
