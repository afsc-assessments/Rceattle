# Log-likelihood of an Rceattle fit

Returns the joint negative-objective from `nlminb`, sign-flipped, as a
`"logLik"` object. `df` is the number of estimated fixed-effect
parameters. The `nobs` attribute is intentionally omitted: counting
"observations" in a stock assessment likelihood (with composition cells,
indices, catches, and priors) is not well-defined, so
[`stats::AIC()`](https://rdrr.io/r/stats/AIC.html) works (uses `df`)
while [`stats::BIC()`](https://rdrr.io/r/stats/AIC.html) does not.

## Usage

``` r
# S3 method for class 'Rceattle'
logLik(object, ...)
```

## Arguments

- object:

  An object of class `"Rceattle"` returned by
  [`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md).

- ...:

  Currently unused.

## Value

An object of class `"logLik"`, or `NULL` if the model was not estimated.
