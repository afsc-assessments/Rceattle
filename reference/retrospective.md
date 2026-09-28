# Retrospective peels

Calculate Mohn's rho and run retrospective peels for an Rceattle model.
The function also evaluates retrospective forecast skill. To evaluate
both retrospective bias and forecast skill, the function uses the map
functionality of TMB to peel the model:

1.  Filters data, filters fixed inputs, and maps out time-varying
    parameters for the peeled years. All time-varying parameters for the
    peeled years are set to the terminal year of the model for that
    peel.

2.  Fits the peeled model.

3.  Turns off all hindcast parameters, turns on F for the peeled years,
    and fits to the peeled catch series to update the "forecast"
    dynamics given projection assumptions and observed catch from the
    peeled years.

## Usage

``` r
retrospective(
  object = NULL,
  peels = 5,
  rescale = FALSE,
  nyrs_forecast = 3,
  cores = NULL,
  getsd = NULL,
  phase = TRUE,
  fit_control = NULL,
  Rceattle = NULL
)
```

## Arguments

- object:

  an Rceattle model fit using
  [`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md).

- peels:

  the number of retrospective peels to use in the calculation of rho and
  for model estimation

- rescale:

  TRUE/FALSE whether to subset and rescale environmental predictors for
  the range of peel years.

- nyrs_forecast:

  Number of forecast years to calculate Mohn's Rho in addition to
  terminal year

- cores:

  Number of cores for the parallel refits. Default `NULL` picks
  `parallel::detectCores() - 6`, capped at 2 under `R CMD check` (which
  sets `_R_CHECK_LIMIT_CORES_`). Set to 1 to force sequential execution.

- getsd:

  whether each peel runs
  [`TMB::sdreport`](https://rdrr.io/pkg/TMB/man/sdreport.html) (standard
  errors). Costs an extra model build per peel; see Details. Mohn's rho
  uses only point estimates, so `FALSE` is faster with no effect on rho.
  Default `NULL` inherits the input model's setting (`TRUE` if it was
  fit with `getsd = TRUE`, i.e. holds an `sdrep`); the returned peel
  models then hold standard errors only when `getsd` is `TRUE`.

- phase:

  whether each peel is refitted in phases (default `TRUE`). A peel
  restarts from the unpeeled fit's starting values with a year removed;
  without phasing the parameters barely move, the peels sit on top of
  the full model, and Mohn's rho is biased towards zero. Change it only
  deliberately.

- fit_control:

  optional
  [`fit_control()`](https://afsc-assessments.github.io/Rceattle/reference/fit_control.md)
  bundle for the refits. Only `phase` and `getsd` are read; see **What
  [`fit_control()`](https://afsc-assessments.github.io/Rceattle/reference/fit_control.md)
  reaches**.

- Rceattle:

  deprecated name for `object`, still accepted so existing scripts keep
  working. Supplying both is an error.

## Value

a list of 1. list of Rceattle models and 2. vector of Mohn's rho for
each species.

A peel that did not converge is dropped, so `Rceattle_list` can be
shorter than `peels + 1` (a message reports how many). Each entry is
named for its own terminal year (`Year_2017`, ...) rather than by
position, so index it by name, `Rceattle_list[[3]]` is not necessarily
the 3-year peel. With no peel left, Mohn's rho is `NaN` and the function
warns.

Each peel reports its own terminal year as `data_list$endyr`, so plots
draw it only as far as it was fit and the peels fan out.

A peel still estimates the years it dropped. They are its retrospective
forecast, fit to the observed catch with the survey and composition data
withheld. Their recruitment deviation is the one
[`sample_rec()`](https://afsc-assessments.github.io/Rceattle/reference/sample_rec.md)
sets with `sample_rec = FALSE`, computed from the peel's own fit; a
penalty-form peel with no penalty years averages over its own years
after the first (or over its one year, for a peel that keeps a single
year), with a warning. Three years therefore matter, and each peel has
all three:

- `endyr`, `endyr_peel`:

  the peel's terminal year, what it was fit through. Equal to each
  other.

- `endyr_full`:

  the unpeeled model's terminal year, where the retrospective forecast
  ends.

- `projyr`:

  the end of the harvest-control-rule projection.

So the forecast years are those after `endyr_peel` through `endyr_full`,
and the projection follows through `projyr`; `incl_proj = TRUE` plots
both. Take the forecast years as
`endyr_peel + seq_len(endyr_full - endyr_peel)`, which is empty for the
unpeeled model, rather than `(endyr_peel + 1):endyr_full`, which counts
*down* there.

Mohn's rho is computed from `endyr_peel` and is unaffected by any of
this.

Catchability is estimated only for a fleet that holds fitted index rows
(see
[`build_map`](https://afsc-assessments.github.io/Rceattle/reference/build_map.md)),
and a peel moves `endyr`. A survey whose index observations all fall in
the peeled-off years therefore has no q estimated in that peel, the
parameter count is not constant across peels. That is deliberate: a q
with no index to inform it is a flat direction in the likelihood. It
does not affect Mohn's rho, which is computed from SSB, but it does mean
`npar` and the reported catchability differ between a shallow and a deep
peel for such a fleet.

## Details

Each peel is fitted twice: a peeled hindcast, then a forecast refit that
estimates only the peeled years' F. The second holds every hindcast
parameter fixed, so on its own it reports a standard error of zero for
the whole hindcast. Under `getsd = TRUE` the peel is therefore rebuilt
at those same parameters with the hindcast free in the map and reported
from there. Nothing is re-estimated, so no point estimate moves.

## Examples

``` r
# \donttest{
data(BS2017SS)
ss_run <- fit_mod(data_list = BS2017SS,
    inits = NULL, file = NULL,
    estimateMode = 0, random_rec = FALSE,
    msmMode = 0, avgnMode = 0,
    phase = FALSE, verbose = 0)
#> Warning: Passing ‘phase’, ‘verbose’ directly to fit_mod() is deprecated and will be removed in a future release. Bundle these into fit_control() instead, e.g. fit_control(phase = ..., verbose = ...). Forwarding for now.
#> `age_trans_matrix` data does not span range of age for species 1 will fill with 0s
retro <- retrospective(ss_run, peels = 10)
# }
```
