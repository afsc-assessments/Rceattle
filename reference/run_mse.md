# Run a management strategy evaluation

Runs a forward-projecting management strategy evaluation (MSE).
Projected selectivity, catchability, foraging days and weight-at-age are
held at the operating model's terminal hindcast year, survey SD is the
average over the historical series, and composition sample size is held
at the last year. There is no implementation error and no observation
error on catch. Fixed numbers-at-age (`estDynamics > 0`) are an input
rather than a forecast, so they are NOT held forward: supply
`NByageFixed` rows for every projection year, or the uncovered years
read 0 and that species contributes nothing, with a warning naming it.

## Usage

``` r
run_mse(
  om,
  em,
  nsim = 10,
  start_sim = 1,
  assessment_period = 1,
  sampling_period = 1,
  simulate_data = TRUE,
  regenerate_past = FALSE,
  sample_rec = TRUE,
  rec_trend = 0,
  fut_sample = 1,
  cap = NULL,
  catch_mult = NULL,
  seed = 666,
  regenerate_seed = seed,
  loopnum = 1,
  file = NULL,
  dir = NULL,
  timeout = 999,
  endyr = NA,
  cores = NULL
)
```

## Arguments

- om:

  CEATTLE model object exported from `Rceattle`

- em:

  CEATTLE model object exported from `Rceattle`

- nsim:

  Number of simulations to run (default 10)

- start_sim:

  First simulation number to start at, useful when a run stops at a
  particular seed (default 1).

- assessment_period:

  Assessment schedule. A single number is the period in years between
  assessments, counted from the operating model's terminal year; a
  vector is the explicit set of assessment years, which must name two or
  more years inside the projection horizon. Default 1.

- sampling_period:

  Period of years data sampling is conducted. Single value or vector the
  same length as the number of fleets.

- simulate_data:

  Include simulated random error proportional to that estimated or
  provided for the data from the OM.

- regenerate_past:

  Refit the EM to historical conditioning data before the MSE, generated
  from the OM with or without sampling error according to
  `simulate_data`.

- sample_rec:

  Include resampled hindcast recruitment deviations in the OM
  projection. Resampled rather than drawn from N(0, sigmaR) because the
  initial deviations bias R0 low. `FALSE` uses the single deviation
  described in
  [`sample_rec()`](https://afsc-assessments.github.io/Rceattle/reference/sample_rec.md).

- rec_trend:

  Linear change in mean recruitment from `endyr` to `projyr`, as the
  terminal multiplier
  `mean rec * (1 + (rec_trend/projection years) * 1:projection years)`.
  Length 1 or `nspp`.

- fut_sample:

  Future sampling effort relative to the last year:
  `Log_sd * 1 / fut_sample` for the index and `Sample_size * fut_sample`
  for comps.

- cap:

  A ceiling on projected catch. A single number is an ANNUAL ceiling on
  total removals summed across species, shared between them in
  proportion to the catch the control rule recommended; a vector of
  length `nspp` is a separate ceiling per species. Default NULL.

- catch_mult:

  A multiplier on the catch the control rule recommends, applied before
  `cap` and before the exploitable-biomass limit. A single number or a
  vector of length `nspp` applies in every projection year; a
  `data.frame` with columns `Year`, `Species` and `mult` applies only to
  the pairs it lists, and any pair it omits is multiplied by 1. `Year`
  must be a year the schedule covers. Default NULL.

- seed:

  seed for the simulation

- regenerate_seed:

  seed for regenerating data

- loopnum:

  Number of times to restart optimization; `loopnum = 3` sometimes
  reaches a lower final gradient than `1`.

- file:

  (Optional) Filename where each OM simulation with EMs is saved. If
  NULL, no files are saved.

- dir:

  (Optional) Directory where each OM simulation is saved

- timeout:

  Minutes estimation runs before a sim is stopped (default 999)

- endyr:

  Terminal year of the MSE projection. Default NA uses `projyr` from the
  operating model.

- cores:

  Cores for the parallel simulations. `NULL` (default) picks
  `parallel::detectCores() - 6`, capped at 2 under `R CMD check`; 1
  forces sequential.

## Value

A list of operating models (differing by simulated recruitment, per
`nsim`) and estimation models fit to each (differing by terminal year).

## Details

Designing a schedule, comparing monitoring scenarios, and what the
common random numbers do and do not share across schedules are in
[`vignette("hcrs-and-mses")`](https://afsc-assessments.github.io/Rceattle/articles/hcrs-and-mses.md).

One piece of arithmetic to get right when buffering a missed assessment:
the assessment in year `Y` sets catch for `Y + 1` onward. So dropping
2031 from a biennial cycle leaves **2032 and 2033** on 2029's advice,
and 2031 itself is set by the 2029 assessment either way – buffering
2031 changes a year the missed assessment never touched.

Catch recorded past a model's terminal year is blanked to `NA` at setup,
with a warning naming the years: the MSE sets those years from the
control rule, and the likelihood never scored them, since it fits only
`Year <= endyr`. That is what a workbook looks like when `endyr` has
fallen behind the catch series, and it is worth resolving first, because
conditioning on those years is a different question from projecting over
them.

Each assessment rebuilds the map as `endyr` advances. A parameter that
`em`'s map fixes and
[`build_map()`](https://afsc-assessments.github.io/Rceattle/reference/build_map.md)
would estimate stays fixed at `em`'s value, and `run_mse()` reports how
many. In years after `em`'s `endyr`, a yearly deviation fixed for every
earlier year can still be estimated, under its usual penalty; a year
joining a time block fixed for every earlier year stays fixed at that
block's value. A new year's block comes from the OM's index rows on a
survey and the EM's catch rows on a fishery. A block fixed on only one
of the fleets sharing it stays estimated on the others and in its new
years. Nothing else in a hand-edited map carries forward: shared values,
and parameters freed where
[`build_map()`](https://afsc-assessments.github.io/Rceattle/reference/build_map.md)
would fix them, return to the
[`build_map()`](https://afsc-assessments.github.io/Rceattle/reference/build_map.md)
default. Nothing carries over from an `em` fitted with
`estimateMode = "Projection"`, even under a supplied map, or from a
re-assessment stored by `run_mse()`, which has no map.

## Examples

``` r
if (FALSE) { # \dontrun{
data(BS2017SS)
om <- fit_mod(BS2017SS, estimateMode = "Estimate", HCR = build_hcr("NPFMC"))
# Closed loop: assess every 2 years, survey every year, 10 simulations.
mse <- run_mse(om = om, em = om, nsim = 10,
               assessment_period = 2, sampling_period = 1)
mse_summary(mse)$species
} # }
```
