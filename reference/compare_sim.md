# Evaluate simulation performance

Function to evaluate the simulation performance with regard to bias
using the median relative error (MRE) and precision using the
coefficient of variation.

## Usage

``` r
compare_sim(operating_mod, simulation_mods, object = "quantities")
```

## Arguments

- operating_mod:

  CEATTLE model object exported from `Rceattle` to be used as the
  operating model

- simulation_mods:

  List of CEATTLE model objects exported from `Rceattle` fit to
  simulated data

- object:

  character string specifying which part of the model to compare
  (default = "quantities")

## Value

A data frame summarizing simulation performance metrics

## Note

Compares against `operating_mod` only, so it is valid when the
replicates redrew the observations alone. If they were produced with
`sim_mod(process = )` / `self_test(process = )`, the operating model's
deviations are no longer what generated the data and the bias this
reports is an artefact. Compare against `attr(sims, "process_sim")` in
that case. Passing a
[`self_test`](https://afsc-assessments.github.io/Rceattle/reference/self_test.md)
result that holds process draws warns, since the attribute makes that
case detectable rather than only documented.
