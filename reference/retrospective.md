# Retrospective peels

Refit an Rceattle model with its last 1, 2, ..., `peels` years of data
removed, and calculate Mohn's rho (Mohn 1999), averaged over peels as in
Hurtado-Ferro et al. (2015), for the terminal year and the
`nyrs_forecast` years after it. Each peel is fitted twice: its hindcast,
then F in the removed years, fitted to their observed catch with the
hindcast held fixed.

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
  projected_env_data_known = FALSE,
  Rceattle = NULL
)
```

## Arguments

- object:

  an Rceattle model fit using
  [`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md).

- peels:

  number of retrospective peels (default 5).

- rescale:

  which `env_data` columns to standardize on each peel's retained years:
  `FALSE` (default) none, `TRUE` every numeric column except `Year`,
  random-effect grouping columns and columns a linkage reads as a
  category, or a character vector of column names.

- nyrs_forecast:

  number of years after each peel's terminal year for which Mohn's rho
  is also calculated (default 3).

- cores:

  Number of cores for the parallel refits. Default `NULL` picks
  `parallel::detectCores() - 6`, capped at 2 under `R CMD check` (which
  sets `_R_CHECK_LIMIT_CORES_`). Set to 1 to force sequential execution.

- getsd:

  whether each peel reports standard errors; `NULL` (default) does so
  only if `object` holds an `sdrep`.

- phase:

  whether each peel is refitted in phases (default `TRUE`); a peel
  starts from the unpeeled model's estimates.

- fit_control:

  optional
  [`fit_control()`](https://afsc-assessments.github.io/Rceattle/reference/fit_control.md)
  bundle for the refits. Only `phase` and `getsd` are read; see **What
  [`fit_control()`](https://afsc-assessments.github.io/Rceattle/reference/fit_control.md)
  reaches**.

- projected_env_data_known:

  `env_data` columns whose values after each peel's terminal year were
  known then, such as a climate projection, so the peels keep them (a
  categorical linkage term or a q or selectivity random effect keeps its
  terminal-year value regardless): `FALSE` (default) none, `TRUE` all,
  or a character vector of column names.

- Rceattle:

  deprecated name for `object`, still accepted so existing scripts keep
  working. Supplying both is an error.

## Value

An `Rceattle_retro` list:

- `Rceattle_list`:

  each kept peel, deepest first, then the unpeeled model, named for its
  terminal year (e.g. `Year_2016`); index by name, since a dropped peel
  leaves a gap.

- `mohns`:

  Mohn's rho, the mean over peels of (peel - unpeeled) / unpeeled, by
  `Object` (`biomass`, `ssb`, `R`, `F_spp`) and `Forecast year` (0 = the
  peel's terminal year), one column per species, with `N` the number of
  peels averaged (`NaN` where `N` is 0, with no warning unless no peel
  converged).

- `peels_requested`:

  the `peels` asked for.

## Details

**Hindcast and forecast years.** A peel's hindcast is the model as it
would have been fitted at its terminal year: it ends there, so no later
index, composition or catch, and no parameter for a later year, enters
its estimates. It starts from the unpeeled estimates for its own years
and uses the unpeeled map cut to them. Its removed years are then its
retrospective forecast: F is fitted to the observed catch with the
hindcast fixed at its estimates. Weight-at-age, empirical selectivity
and ration stay at the terminal year, and the q and selectivity
deviations keep their terminal-year values (a form whose deviations are
year-to-year increments takes an increment of 0). M1 deviations take the
peel's projection values: 0 for one estimated by year, an input series
as given. A q or selectivity linkage random effect keeps its
terminal-year level too, as the projection holds selectivity. Any other
linkage effect takes the peel's own estimate for each removed year, as
the projection does: a random walk its terminal value, an AR1 its decay
toward 0, an IID effect 0. Recruitment takes the deviation
[`sample_rec()`](https://afsc-assessments.github.io/Rceattle/reference/sample_rec.md)
uses with `sample_rec = FALSE`, computed from the peel's fit (for a
penalty-form curve with no penalty years in the peel, the mean over its
own years after the first). Each peel's `data_list` holds `endyr` and
`endyr_peel` (its terminal year), `endyr_full` (the unpeeled terminal
year, where the forecast ends) and `projyr` (the end of the projection).
The forecast years are `endyr_peel + seq_len(endyr_full - endyr_peel)`,
which is empty for the unpeeled model, where
`(endyr_peel + 1):endyr_full` would count down.

**Dropped and flagged peels.** A peel whose hindcast did not converge is
dropped, with a warning; under `getsd = TRUE` that includes a hindcast
Hessian that is not positive definite, as when a category level begins
after the peel's terminal year and so leaves its estimated coefficient
with no data (`getsd = FALSE` keeps that peel). A peel whose forecast
refit did not fit F to the catch (gradient check, or a peeled F at a
bound) is kept but left out of the `$mohns` rows with `Forecast year`
above 0, and `attr(peel, "forecast_refit_error")` gives the reason.
Check it before reading that peel's forecast years: only
[`plot_timeseries()`](https://afsc-assessments.github.io/Rceattle/reference/plot_timeseries.md)
and its wrappers warn about them.

**Standard errors.** Mohn's rho reads only point estimates, so
`getsd = FALSE` is faster. Under `getsd = TRUE` each peel is rebuilt at
its estimates with the hindcast parameters treated as estimated, so no
estimate moves. A fitted year's band is the peel's hindcast's own,
except a depletion relative to last-projection-year biomass
(`HCR = "NoFishing"`, `msmMode > 0`, `DynamicHCR = FALSE`), which reads
the peeled F. A forecast year's band treats the catch as known and
carries the hindcast's uncertainty through the F needed to take it,
without recruitment variability, so it is narrower than a prediction
interval. If the rebuild fails (an error, or no finite standard error),
the peel is kept, but plots,
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) and
[`vcov()`](https://rdrr.io/r/stats/vcov.html) ignore its `sdrep`, which
may not hold valid standard errors; `attr(peel, "sd_report_error")`
gives the reason. [`print()`](https://rdrr.io/r/base/print.html) counts
the peels carrying either attribute; `retrospective()` warns about
neither.

**Environmental data.** A peel's forecast should read only what was
known at its terminal year, so after that year a peel treats each
`env_data` column by how the model reads it:

- an observed covariate (`linkage_spec(observe = )`) loses its
  observations;

- a column the model reads directly (a `Cindex` temperature, an
  environmental catchability or M1) takes the retained years' mean, as
  in any missing year;

- a numeric fixed-effect linkage covariate, a 0/1 block column included,
  takes the retained years' mean, with a warning; a term computed over
  the whole column, such as
  [`scale()`](https://rdrr.io/r/base/scale.html), then changes in the
  retained years too;

- a fixed-effect category, such as `cut(Year, ...)`, `factor(regime)` or
  a logical column, holds its terminal-year value
  ([`linkage_spec()`](https://afsc-assessments.github.io/Rceattle/reference/linkage_spec.md)),
  though one cut from the whole column, such as `cut(temp, 3)`, still
  sets its breaks from the whole column, post-peel years included;

- a random-effect grouping column, a column with no retained value, and
  a fixed-effect covariate that is also observed, also read as a
  category (the category term is still held), or logical but read as a
  number (`as.numeric(flag)`) keep their post-peel values, with a
  warning; Mohn's rho is then conditional on them (a q or selectivity
  random effect takes its terminal-year level whatever its grouping
  column says);

- a column listed in `projected_env_data_known` keeps its values with no
  warning, though a category is still held.

Refitting a returned peel's `data_list` with
[`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)
does not hold its categories, while its numeric covariates stay
mean-filled; a diagnostic such as
[`jitter()`](https://afsc-assessments.github.io/Rceattle/reference/jitter.md)
or [`profile()`](https://rdrr.io/r/stats/profile.html) does not work on
a returned peel.

`rescale` standardizes a column on the retained years and applies that
centre and scale to the whole column (a column with fewer than two
distinct retained values stays in its units, with a warning). Name the
columns when `env_data` also holds an input that is not a covariate,
such as a `Cindex` temperature in degrees C, which `TRUE` would turn
into z-scores. The unpeeled model is not rescaled, so standardize those
columns the same way before fitting it.

**Catchability.** A peel's hindcast uses the unpeeled model's map, cut
to its years, so a survey whose index falls only in the removed years
keeps a free q with no index to inform it, which under `getsd = TRUE`
can leave the Hessian not positive definite and drop the peel.

## References

Mohn, R. 1999. The retrospective problem in sequential population
analysis: an investigation using cod fishery and simulated data. ICES J.
Mar. Sci. 56:473-488.

Hurtado-Ferro, F., et al. 2015. Looking in the rear-view mirror: bias
and retrospective patterns in integrated, age-structured stock
assessment models. ICES J. Mar. Sci. 72:99-110.

## Examples

``` r
# \donttest{
data(BS2017SS)
ss_run <- fit_mod(data_list = BS2017SS, inits = NULL, file = NULL,
                  estimateMode = 0, random_rec = FALSE, msmMode = 0,
                  fit_control = fit_control(phase = FALSE, verbose = 0))
#> `age_trans_matrix` data does not span range of age for species 1 will fill with 0s
retro <- retrospective(ss_run, peels = 10)
# }
```
