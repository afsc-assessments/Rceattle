# Build a model-configuration slot for a data list

Bundles the model-structure arguments of
[`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md),
predation mode, initialization, harvest control rule, and the
`build_*()` process specifications, into a single validated object that
can be stored on a data list (`data_list$model_config`, e.g. via
`build_data(..., model_config = model_config(...))`). A configuration
then travels with the data instead of living only in the
[`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)
call.

## Usage

``` r
model_config(
  msmMode = 0,
  initMode = "NonEquilibrium",
  avgnMode = 0,
  suitMode = 0,
  niter = 3,
  HCR = build_hcr(),
  recFun = build_srr(),
  M1Fun = build_M1(),
  growthFun = build_growth(),
  qFun = build_catchability(),
  selFun = build_selectivity(),
  compFun = build_composition()
)
```

## Arguments

- msmMode:

  Predation-mortality mode (see
  [`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)).
  Default `0`.

- initMode:

  Population initialization (see
  [`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)).
  Default `"NonEquilibrium"`.

- avgnMode:

  Average-N mode. Default `0`.

- suitMode:

  Suitability mode. Default `0`.

- niter:

  Number of predation iterations. Default `3`.

- HCR:

  A harvest-control-rule specification from
  [`build_hcr()`](https://afsc-assessments.github.io/Rceattle/reference/build_hcr.md).

- recFun:

  A stock-recruit specification from
  [`build_srr()`](https://afsc-assessments.github.io/Rceattle/reference/build_srr.md).

- M1Fun:

  A natural-mortality specification from
  [`build_M1()`](https://afsc-assessments.github.io/Rceattle/reference/build_M1.md).

- growthFun:

  A growth specification from
  [`build_growth()`](https://afsc-assessments.github.io/Rceattle/reference/build_growth.md).

- qFun:

  A catchability specification from
  [`build_catchability()`](https://afsc-assessments.github.io/Rceattle/reference/build_catchability.md).

- selFun:

  A selectivity specification from
  [`build_selectivity()`](https://afsc-assessments.github.io/Rceattle/reference/build_selectivity.md).

- compFun:

  A composition specification from
  [`build_composition()`](https://afsc-assessments.github.io/Rceattle/reference/build_composition.md).

## Value

An object of class `"Rceattle_model_config"` (a named list) to be stored
as `data_list$model_config`.

## Details

The defaults are exactly
[`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)'s
own argument defaults, so a data list holding `model_config()` fits
identically to one with no slot at all. When a data list has a
`model_config`,
[`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)
reads each field only for arguments the caller did **not** pass
explicitly, an argument passed to
[`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)
always wins, even when passed at its default. Omit an argument to let
the stored configuration take effect for that field.

## Persistence

The configuration is code-side model structure, not one of the workbook
data sheets, so it is **not** written by
[`write_data()`](https://afsc-assessments.github.io/Rceattle/reference/write_data.md)
and does not survive an xlsx round-trip,
`build_data(base = x, model_config = cfg)` piped through
[`write_data()`](https://afsc-assessments.github.io/Rceattle/reference/write_data.md)
then
[`read_data()`](https://afsc-assessments.github.io/Rceattle/reference/read_data.md)
returns without the slot. Re-attach it in code, store it alongside the
data, or persist it as a documented, git-diffable YAML with
[`save_config()`](https://afsc-assessments.github.io/Rceattle/reference/save_config.md)
/
[`load_config()`](https://afsc-assessments.github.io/Rceattle/reference/load_config.md)
and apply it with
`fit_mod(data_list, config = load_config("run.yaml"))`.

## See also

[`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md),
[`build_data()`](https://afsc-assessments.github.io/Rceattle/reference/build_data.md),
[`save_config()`](https://afsc-assessments.github.io/Rceattle/reference/save_config.md),
[`load_config()`](https://afsc-assessments.github.io/Rceattle/reference/load_config.md).

## Examples

``` r
cfg <- model_config(msmMode = 1, initMode = "NonEquilibrium")
dat <- build_data(base = BS2017MS, model_config = cfg)
# fit_mod(dat) would then fit as multispecies without passing msmMode.
```
