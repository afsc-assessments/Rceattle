# Convergence diagnostics for a fitted Rceattle model

Runs the post-fit convergence battery and returns a single structured
object. Each check yields a record with a common schema (`id`, `tier`,
`severity`, `message`, `data`); `severity` is one of `"OK"`, `"NOTE"`,
`"WARN"`, `"FAIL"`. The object's `status` is the worst severity present.

## Usage

``` r
convergence_diagnostics(object, ...)
```

## Arguments

- object:

  An object of class `"Rceattle"` returned by
  [`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md).

- ...:

  Currently unused.

## Value

An object of class `"Rceattle_convergence"`: a list with `status`
(overall worst severity) and `checks` (named list of records).

## Details

[`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)
runs this automatically and attaches the result as `fit$convergence`;
call `convergence_diagnostics()` directly to re-run it on any fit.
Checks cover the optimizer gradient, a requested `sdreport` that did not
return, an `sdreport` that was never requested, Hessian
positive-definiteness and conditioning, parameters on bounds, parameters
past them, a deviation variance estimated to zero, phasing, parameter
estimability, a numbers-at-age, Ricker-intercept or recruitment floor
that was reached, and the stock-recruit curve.
