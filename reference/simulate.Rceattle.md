# Simulate data from a fitted Rceattle model

The [`stats::simulate()`](https://rdrr.io/r/stats/simulate.html) method
for CEATTLE fits: draws `nsim` replicate data sets from the fitted
observation model, and optionally from the process model as well.

## Usage

``` r
# S3 method for class 'Rceattle'
simulate(object, nsim = 1, seed = NULL, process = FALSE, ...)
```

## Arguments

- object:

  An object of class `"Rceattle"` returned by
  [`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md).

- nsim:

  Number of replicate data sets to draw. Default 1.

- seed:

  Optional seed. The caller's random state is restored afterwards, per
  [`stats::simulate()`](https://rdrr.io/r/stats/simulate.html). Recorded
  as `attr(, "seed")`.

- process:

  Which process error to redraw alongside the observations. `FALSE`
  (default) keeps the fitted deviations, so replicates differ only in
  observation error. See
  [`sim_mod()`](https://afsc-assessments.github.io/Rceattle/reference/sim_mod.md)
  for the alternatives.

- ...:

  Currently unused.

## Value

A list of `nsim` `data_list` objects, always a list, including at
`nsim = 1`, so callers do not have to special-case the length. When
`process` redrew something, each element holds the deviations that
generated it as `attr(, "process_sim")`: a named list of whichever of
`rec_dev`, `init_dev`, `log_M1_dev` and `beta_linkage_re` were drawn,
each with a `_drawn` logical of the same shape marking the cells the
draw touched. Compare estimates against those rather than against
`object`, see
[`sim_mod()`](https://afsc-assessments.github.io/Rceattle/reference/sim_mod.md).

## Details

A wrapper on
[`sim_mod()`](https://afsc-assessments.github.io/Rceattle/reference/sim_mod.md),
which documents the observation model in full. For expected values
rather than draws, call `sim_mod(simulate = FALSE)`.

Draws are taken by the TMB model, so
[`simulate()`](https://rdrr.io/r/stats/simulate.html) needs a live
`$obj`: a model loaded from disk has one, a
[`model_average()`](https://afsc-assessments.github.io/Rceattle/reference/model_average.md)
result does not, see
[`sim_mod()`](https://afsc-assessments.github.io/Rceattle/reference/sim_mod.md).

## See also

[`sim_mod()`](https://afsc-assessments.github.io/Rceattle/reference/sim_mod.md)
for the observation model and the `process` options,
[`self_test()`](https://afsc-assessments.github.io/Rceattle/reference/self_test.md)
for simulating and refitting in one step.

## Examples

``` r
if (FALSE) { # \dontrun{
fit  <- fit_mod(data_list = BS2017SS, estimateMode = 1)
reps <- simulate(fit, nsim = 10, seed = 1)
# ...and with recruitment redrawn as well:
reps <- simulate(fit, nsim = 10, seed = 1, process = "recruitment")
truth <- attr(reps[[1]], "process_sim")$rec_dev
} # }
```
