# Variance-covariance matrix for an Rceattle fit

Returns the fixed-effect covariance matrix produced by
[`TMB::sdreport()`](https://rdrr.io/pkg/TMB/man/sdreport.html).
Random-effect covariance is not returned here, use `object$sdrep` for
the full report.

## Usage

``` r
# S3 method for class 'Rceattle'
vcov(object, ...)
```

## Arguments

- object:

  An object of class `"Rceattle"` returned by
  [`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md).

- ...:

  Currently unused.

## Value

A numeric matrix, or `NULL` if `sdreport` was not run (i.e. the fit was
produced with `getsd = FALSE`).
