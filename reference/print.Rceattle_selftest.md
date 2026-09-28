# Print method for a self-test

Reports what a self-test is run to find out: how many simulations the
estimator brought back to an optimum, and under what the replicates were
generated. The returned list cannot say the first on its own ,
non-converged runs are dropped before it is returned, so the number of
fits returned is not the number attempted.

## Usage

``` r
# S3 method for class 'Rceattle_selftest'
print(x, rate = 0.9, ...)
```

## Arguments

- x:

  A `"Rceattle_selftest"` object from
  [`self_test()`](https://afsc-assessments.github.io/Rceattle/reference/self_test.md).

- rate:

  Fraction of simulations that must converge before the run counts as
  clean. Default `0.9`, matching
  [`jitter()`](https://afsc-assessments.github.io/Rceattle/reference/jitter.md)'s.

- ...:

  Currently unused.

## Value

`x`, invisibly.

## Details

The status is the convergence RATE, since that is the quantity the run
produces: `FAIL` if nothing converged, `WARN` below `rate`, `NOTE` if
any replicate was dropped, `OK` otherwise. The table beneath tallies
each returned fit's own `$convergence$status`, which is a separate
question, a replicate can reach a zero gradient and still hold a `NOTE`
or `WARN`, and those are not folded into the header for that reason.

The last line says what generated the replicates, because it decides
what the spread means. With processes held fixed (the default)
[`sim_mod()`](https://afsc-assessments.github.io/Rceattle/reference/sim_mod.md)
redraws the observations alone, so the spread reflects observation error
only and is a lower bound on estimation uncertainty, do not read it
against the model's own uncertainty bands, which include process error.
See **Interpreting the spread** in
[`self_test()`](https://afsc-assessments.github.io/Rceattle/reference/self_test.md).
