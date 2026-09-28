# Save a model run configuration to a documented YAML file

Round-trips a full run configuration, the
[`model_config()`](https://afsc-assessments.github.io/Rceattle/reference/model_config.md)
structure plus the estimation controls and
[`fit_control()`](https://afsc-assessments.github.io/Rceattle/reference/fit_control.md)
bundle, to a git-diffable YAML file, with each field's documentation
emitted as a comment and a spec-tree + provenance header. Only fields
that differ from their defaults are written, so two configurations diff
to just their real differences. The parameter values
(`inits`/`map`/`bounds`) are NOT stored; pair the config with a saved
fit for those.

## Usage

``` r
save_config(x, file = "Rceattle_config.yaml", ...)
```

## Arguments

- x:

  A fitted Rceattle object, a data list holding `$model_config`, an
  `Rceattle_run_config`, or an `Rceattle_model_config`.

- file:

  Output path for the `.yaml` file.

- ...:

  Estimation controls / `fit_control` to record (passed to
  [`run_config()`](https://afsc-assessments.github.io/Rceattle/reference/run_config.md)).

## Value

Invisibly, the `Rceattle_run_config` that was written.

## See also

[`load_config()`](https://afsc-assessments.github.io/Rceattle/reference/load_config.md),
[`run_config()`](https://afsc-assessments.github.io/Rceattle/reference/run_config.md),
[`model_config()`](https://afsc-assessments.github.io/Rceattle/reference/model_config.md),
[`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md).

## Examples

``` r
cfg <- model_config(msmMode = 1, initMode = "OffsetEquilibrium")
f <- file.path(tempdir(), "run.yaml")
save_config(cfg, f)
identical(load_config(f)$model_config$msmMode, 1)
#> [1] TRUE
```
