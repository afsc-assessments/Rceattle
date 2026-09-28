# Load a model run configuration from a YAML file

The inverse of
[`save_config()`](https://afsc-assessments.github.io/Rceattle/reference/save_config.md):
reads a run-config YAML file and rebuilds the
[`model_config()`](https://afsc-assessments.github.io/Rceattle/reference/model_config.md)
structure, estimation controls, and
[`fit_control()`](https://afsc-assessments.github.io/Rceattle/reference/fit_control.md)
bundle (reconstructing linkage formulas and priors). Attach the result
to a fit with `fit_mod(data_list, config = load_config("run.yaml"))`, or
read `$model_config` off it to attach to a data list.

## Usage

``` r
load_config(file)
```

## Arguments

- file:

  Path to a YAML file written by
  [`save_config()`](https://afsc-assessments.github.io/Rceattle/reference/save_config.md).

## Value

An `Rceattle_run_config`.

## See also

[`save_config()`](https://afsc-assessments.github.io/Rceattle/reference/save_config.md),
[`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md).

## Examples

``` r
f <- file.path(tempdir(), "run.yaml")
save_config(model_config(msmMode = 1, initMode = "Equilibrium"), f)

cfg <- load_config(f)
cfg$model_config$msmMode
#> [1] 1
# Apply with fit_mod(data_list, config = cfg).
unlink(f)
```
