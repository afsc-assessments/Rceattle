# Extract estimated parameters from an Rceattle fit

Extract estimated parameters from an Rceattle fit

## Usage

``` r
# S3 method for class 'Rceattle'
coef(object, ...)
```

## Arguments

- object:

  An object of class `"Rceattle"` returned by
  [`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md).

- ...:

  Currently unused.

## Value

A named numeric vector of estimated fixed-effect parameters (the
optimizer's `par`). Returns `NULL` if the model was not estimated.
