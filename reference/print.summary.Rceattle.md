# Print method for an Rceattle model summary

Print method for an Rceattle model summary

## Usage

``` r
# S3 method for class 'summary.Rceattle'
print(x, n = 10, ...)
```

## Arguments

- x:

  A `"summary.Rceattle"` object from
  [`summary.Rceattle()`](https://afsc-assessments.github.io/Rceattle/reference/summary.Rceattle.md).

- n:

  Number of parameters to show, in
  [`coef()`](https://rdrr.io/r/stats/coef.html) order. Default 10; use
  `Inf` for all, or take `x$coefficients`.

- ...:

  Currently unused.

## Value

`x`, invisibly.
