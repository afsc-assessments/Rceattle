# Extract the run configuration from a fit, data list, or config object

Returns the
[`model_config()`](https://afsc-assessments.github.io/Rceattle/reference/model_config.md)
structure plus the estimation controls and
[`fit_control()`](https://afsc-assessments.github.io/Rceattle/reference/fit_control.md)
bundle as a single `Rceattle_run_config`. Accepts a fitted Rceattle
object, a data list holding `$model_config`, an `Rceattle_run_config`,
or an `Rceattle_model_config`. Estimation controls and `fit_control`
supplied via `...` override any found on the object.

## Usage

``` r
run_config(x, ...)
```

## Arguments

- x:

  A fitted Rceattle object, a data list, an `Rceattle_run_config`, or an
  `Rceattle_model_config`.

- ...:

  Estimation controls / `fit_control` to override (`estimateMode`,
  `random_rec`, `random_q`, `random_sel`, `suit_styr`, `suit_endyr`,
  `fit_control`).

## Value

An `Rceattle_run_config`.

## See also

[`save_config()`](https://afsc-assessments.github.io/Rceattle/reference/save_config.md),
[`load_config()`](https://afsc-assessments.github.io/Rceattle/reference/load_config.md),
[`model_config()`](https://afsc-assessments.github.io/Rceattle/reference/model_config.md).

## Examples

``` r
# Coerce a model_config() into a full run configuration, overriding the
# estimation controls. A fitted model or a data list works the same way.
rc <- run_config(model_config(msmMode = 1), estimateMode = "Hindcast")
rc$estimateMode
#> [1] "Hindcast"

# Optimizer and uncertainty settings belong inside fit_control(), not here;
# an unrecognized field is an error rather than a silent no-op.
rc <- run_config(rc, fit_control = fit_control(getsd = FALSE))
```
