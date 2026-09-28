# Print method for an MSE summary

Says what the summary holds. The object is deliberately ragged,
per-species, per-fleet and whole-system metrics have different shapes
and cannot share one frame, so it reports the blocks and their
dimensions rather than printing them end to end.

## Usage

``` r
# S3 method for class 'Rceattle_mse_summary'
print(x, ...)
```

## Arguments

- x:

  An `"Rceattle_mse_summary"` object from
  [`mse_summary()`](https://afsc-assessments.github.io/Rceattle/reference/mse_summary.md).

- ...:

  Currently unused.

## Value

`x`, invisibly.
