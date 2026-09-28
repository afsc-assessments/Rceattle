# Plot method for fitted Rceattle models

Thin S3 dispatcher around the package's existing `plot_*()` functions so
that `plot(fit)` works the way users expect. Pick the panel with `what`;
everything in `...` is forwarded to the underlying function.

## Usage

``` r
# S3 method for class 'Rceattle'
plot(x, what = "biomass", ...)
```

## Arguments

- x:

  An object of class `"Rceattle"` returned by
  [`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md).

- what:

  Character. One of `"biomass"` (default), `"ssb"`, `"recruitment"`,
  `"depletion"` (total biomass / B0), `"ssb_depletion"` (female spawning
  biomass / SB0, the quantity a Tier 3 HCR compares against B40%),
  `"index"`, `"catch"`, `"selectivity"`, `"mortality"`, or `"data"`.

- ...:

  Passed to the underlying plotting function.

## Value

Invisibly returns `NULL`. Called for the side effect of producing a
plot.
