# Retrospective peels

Calculate Mohn's rho and run retrospective peels for an Rceattle model,
and evaluate retrospective forecast skill. Each peel uses TMB's map to:

1.  filter the data and fixed inputs and map out time-varying parameters
    for the peeled years, holding them at that peel's terminal year;

2.  fit the peeled model;

3.  turn the hindcast parameters off and F for the peeled years on, and
    fit to the peeled catch series, updating the "forecast" dynamics
    given the projection assumptions and the observed catch.

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

  number of retrospective peels used for rho and estimation

- rescale:

  TRUE/FALSE whether to subset and rescale environmental predictors over
  the peel years.

- nyrs_forecast:

  forecast years over which to calculate Mohn's rho, in addition to the
  terminal year

- cores:

  Number of cores for the parallel refits. Default `NULL` picks
  `parallel::detectCores() - 6`, capped at 2 under `R CMD check` (which
  sets `_R_CHECK_LIMIT_CORES_`). Set to 1 to force sequential execution.

- getsd:

  whether each peel reports standard errors through
  [`TMB::sdreport`](https://rdrr.io/pkg/TMB/man/sdreport.html), at the
  cost of extra model builds (see Details); `NULL` (default) inherits
  the input model's setting (`TRUE` only if it holds an `sdrep`).

- phase:

  whether each peel is refitted in phases, default `TRUE`. A peel
  restarts from the unpeeled fit's starting values with a year removed,
  so without phasing the parameters barely move, the peels sit on top of
  the full model, and Mohn's rho is biased toward zero. Change it only
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

An `Rceattle_retro` list:

- `Rceattle_list`:

  the unpeeled model and each kept peel, named for its own terminal year
  (`Year_2017`, ...). Index by name: the list runs from the deepest peel
  to the unpeeled model, and a dropped peel leaves a gap, so it can be
  shorter than `peels + 1` and `Rceattle_list[[3]]` is not in general
  the 3-year peel. Each peel reports its terminal year as
  `data_list$endyr`, so plots draw it only as far as it was fit.

- `mohns`:

  Mohn's rho by `Object` and `Forecast year`, one column per species,
  with `N` the number of peels averaged. With no peel left it is `NaN`,
  and the function warns.

- `peels_requested`:

  the `peels` asked for.

A peel can carry `attr(peel, "forecast_refit_error")` or
`attr(peel, "sd_report_error")`, each a one-line reason naming the peel
and its terminal year (see Details); each is absent on every other peel.

## Details

Each peel is fitted twice: a peeled hindcast, then a forecast refit
estimating only the peeled years' F. The second holds every hindcast
parameter fixed, so on its own it reports a standard error of zero
across the hindcast. Under `getsd = TRUE` the peel is therefore rebuilt
at those same parameters with the hindcast parameters treated as
estimated and reported from there; nothing is re-estimated, so no point
estimate moves.

A band in a fitted year is the peeled hindcast's own: the uncertainty of
its parameters given the data through the peel's terminal year. The
exception is a depletion that divides by last-projection-year biomass
(`HCR = "NoFishing"` under `msmMode > 0` and `DynamicHCR = FALSE`),
which reads the peeled F, so its fitted-year bands move too. A band in a
peeled (forecast) year treats that year's catch as known and its F as
inferred from it: as the hindcast parameters vary, F moves, to first
order, to keep the predicted catch at the observed catch, and the delta
method carries that into the band. The peeled years' recruitment
deviations are fixed at the projected deviation, so the band leaves out
recruitment variability there and is narrower than a forecast's
prediction interval.

Where that rebuild fails – it errors, its Hessian is not positive
definite, it reports no finite standard error, or moving F with the
hindcast parameters fails (its forecast years would then leave out F's
response) – the peel keeps its point estimates, `peel$sdrep` keeps what
[`TMB::sdreport()`](https://rdrr.io/pkg/TMB/man/sdreport.html) returned
(`NULL` if the rebuild errored first; under an indefinite Hessian its
numbers are not valid standard errors), and
`attr(peel, "sd_report_error")` gives the reason, so plots,
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) and
[`vcov()`](https://rdrr.io/r/stats/vcov.html) read none of it.
[`print()`](https://rdrr.io/r/base/print.html) counts such peels; no
warning is raised, since the point estimates and Mohn's rho are
unaffected.

Mohn's rho reads point estimates, so `getsd = FALSE` is faster, but
`TRUE` also drops a peel whose hindcast Hessian is not positive
definite, so the two can average over different peels.

A peel whose hindcast did not converge is dropped, with a warning. A
peel whose forecast refit did not fit F to the peeled catch (gradient
check, or a peeled F at a bound) is kept without a warning: its
terminal-year Mohn's rho stands, but it is left out of the `$mohns` rows
with `Forecast year` above 0 (`N` counts the peels used). Check
`attr(peel, "forecast_refit_error")` before reading its forecast years
or their standard errors; [`print()`](https://rdrr.io/r/base/print.html)
counts such peels, and
[`plot_timeseries()`](https://afsc-assessments.github.io/Rceattle/reference/plot_timeseries.md)
and its wrappers warn when they draw those years.

A peel still estimates the years it dropped: that is its retrospective
forecast, fit to observed catch with the survey and composition data
withheld. Its recruitment deviation is the one
[`sample_rec()`](https://afsc-assessments.github.io/Rceattle/reference/sample_rec.md)
sets with `sample_rec = FALSE`, computed from the peel's own fit. Three
years matter, and every peel's `data_list` carries all three:

- `endyr`, `endyr_peel`:

  the peel's terminal year, what it was fit through. Equal to each
  other.

- `endyr_full`:

  the unpeeled model's terminal year, where the retrospective forecast
  ends.

- `projyr`:

  the end of the harvest-control-rule projection.

Take the forecast years as
`endyr_peel + seq_len(endyr_full - endyr_peel)`, which is empty for the
unpeeled model, rather than `(endyr_peel + 1):endyr_full`, which counts
*down* there. Mohn's rho is computed from `endyr_peel` and is unaffected
by any of this.

Catchability is estimated only for a fleet holding fitted index rows
(see
[`build_map`](https://afsc-assessments.github.io/Rceattle/reference/build_map.md)),
and a peel moves `endyr`, so a survey whose index observations all fall
in the peeled-off years has no q estimated in that peel: the parameter
count is not constant across peels. That is deliberate, since a q with
no index to inform it is a flat direction in the likelihood. Rho is
computed from SSB and is unaffected, but `npar` and the reported
catchability do differ between a shallow and a deep peel for such a
fleet.

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
