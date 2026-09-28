# Run a management strategy evaluation

Runs a forward-projecting management strategy evaluation (MSE).
Projected selectivity, catchability, foraging days, and weight-at-age
are held at the operating model's terminal hindcast year. Survey SD is
set to the average over the historical time series, and composition
sample size is held at the last year. There is no implementation error
and no observation error on catch.

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

  First simulation number to start at. Useful if the code stops at
  specific seed/sim (default = 1).

- assessment_period:

  Assessment schedule. A single number is the period in years between
  assessments, counted from the operating model's terminal year. A
  vector is the explicit set of years an assessment is completed.
  Default = 1.

- sampling_period:

  Period of years data sampling is conducted. Single value or vector the
  same length as the number of fleets.

- simulate_data:

  Include simulated random error proportional to that estimated/provided
  for the data from the OM.

- regenerate_past:

  Refits the EM to historical/conditioning data prior to the MSE, where
  the data are generated from the OM with `simulate_data = TRUE` or
  without `simulate_data = FALSE` sampling error.

- sample_rec:

  Include resampled recruitment deviations from the hindcast in the OM
  projection. Resampled deviations are used rather than drawing from
  N(0, sigmaR) because the initial deviations bias R0 low. If FALSE,
  uses the single deviation described in
  [`sample_rec()`](https://afsc-assessments.github.io/Rceattle/reference/sample_rec.md).

- rec_trend:

  Linear increase or decrease in mean recruitment from `endyr` to
  `projyr`. This is the terminal multiplier
  `mean rec * (1 + (rec_trend/projection years) * 1:projection years)`.
  Can be of length 1 or of length nspp. If length 1, all species get the
  same trend.

- fut_sample:

  future sampling effort relative to last year.
  ` Log_sd * 1 / fut_sample` for index and ` Sample_size * fut_sample`
  for comps

- cap:

  A cap on the catch in the projection. Can be a single number applied
  to all species (proportional to recommended catch) or vector of length
  `nspp` applied to each species. Default = NULL

- catch_mult:

  A multiplier applied to the catch the control rule recommends. A
  single number or a vector of length `nspp` applies in every projection
  year; a `data.frame` with columns `Year`, `Species` and `mult` applies
  only in the year and species pairs it lists, where `Year` is a year
  the assessment schedule covers. Default = NULL

- seed:

  seed for the simulation

- regenerate_seed:

  seed for regenerating data

- loopnum:

  number of times to re-start optimization (where `loopnum=3` sometimes
  achieves a lower final gradient than `loopnum=1`)

- file:

  (Optional) Filename where each OM simulation with EMs will be saved.
  If NULL, no files are saved.

- dir:

  (Optional) Directory where each OM simulation is saved

- timeout:

  length of time (minutes) estimation will run before stopping a sim
  (default 999 minutes)

- endyr:

  Terminal year of the MSE projection. Default = NA uses `projyr` from
  the operating model.

- cores:

  Number of cores to use for parallel simulations. Default `NULL` picks
  `parallel::detectCores() - 6`, capped at 2 when running under
  `R CMD check` (which sets `_R_CHECK_LIMIT_CORES_`). Set to 1 to force
  sequential execution.

## Value

A list of operating models (differ by simulated recruitment determined
by `nsim`) and estimation models fit to each operating model (differ by
terminal year).

## Assessment schedule

A single `assessment_period` is a fixed cycle. A vector is the schedule
itself, for a design whose years are not evenly spaced, one assessment
missed inside an otherwise biennial cycle, for instance:

    biennial <- seq(om$data_list$endyr + 2, om$data_list$projyr, by = 2)
    run_mse(om, em, assessment_period = setdiff(biennial, 2031))

Every year must be after the operating model's terminal year and within
the projection horizon; a year outside that window is an error rather
than being dropped. One missed assessment and a permanently longer cycle
are different questions, and the second does not answer the first.

A schedule must name two or more years. One year on its own cannot be
told from a period, and the two readings are a world apart, so it is
refused rather than guessed at.

Two schedules run on the same `seed` share their recruitment deviations,
which are drawn once per replicate before the first assessment, and each
assessment's observation draw is seeded on that assessment's own year
rather than on wherever earlier assessments left the stream. So an
assessment in year `Y` starts from the same place in every schedule that
assesses `Y`, and a divergence at one assessment no longer displaces
every assessment after it.

Common random numbers are not complete, and the gap is worth knowing
before designing a comparison. One
[`sim_mod()`](https://afsc-assessments.github.io/Rceattle/reference/sim_mod.md)
call draws every year in the assessment interval, under that
assessment's seed — so a year sitting *inside* a longer interval is
drawn under a different seed than the same year in a schedule that
assessed it directly. Two schedules that differ in which years they
assess therefore realize different observation error in the years
between, even where the stock is in the same state. Per-observation-year
streams would close it; `inst/dev/TODO-mse-horizon.md` records what that
would take.

Past the point where two schedules genuinely diverge — different catch
taken, different years surveyed — the draws necessarily diverge too.
That is a real difference in the runs, not an artefact.

Make the last assessment year reach the projection horizon. Catch is
only ever filled up to the last assessment, so a schedule that stops
short leaves the trailing years at `NA`, and
[`mse_summary()`](https://afsc-assessments.github.io/Rceattle/reference/mse_summary.md)
summarises over the whole projection — understating Average Catch, Catch
IAV and P(Closed) with nothing in the table to say which years it
covered. A biennial cycle over an odd number of projection years lands
here every time, so this is warned about, for a period and a schedule
alike.

## Catch multiplier

`catch_mult` multiplies the catch the estimation model's control rule
recommends, before `cap` and before the exploitable-biomass limit. Given
as a `data.frame` it applies only in the years and species it lists,
which is how a buffer held for the years advice is stale is expressed.

Mind which years those are. The assessment in year `Y` sets catch for
`Y+1` onward, so a missed assessment in 2031 leaves **2032 and 2033** on
2029's advice — 2031 itself is set by the 2029 assessment either way,
and cutting it changes a year the missed assessment never touched:

    buffer <- expand.grid(Year = 2032:2033, Species = seq_len(om$data_list$nspp))
    buffer$mult <- 0.90
    run_mse(om, em, assessment_period = setdiff(biennial, 2031),
            catch_mult = buffer)

`Species` is the species number, matching the catch data's own column,
and any (year, species) pair the table omits is multiplied by 1. `Year`
must be a year the assessment schedule actually covers — from the
operating model's terminal year to the last assessment, which is where
catch is filled, not the whole projection.

Note that this reduces catch, not ABC. Where realized catch sits well
below ABC, GOA arrowtooth flounder, for one, reducing ABC changes
removals only to the extent the fishery attains it, while reducing catch
changes them in full. Either scale the multiplier by recent attainment,
`1 - (1 - mult) * attainment`, or report the unscaled result as an upper
bound on the effect of the reduction.

## Projection catch

Catch recorded past a model's terminal year is blanked to `NA` at setup,
with a warning naming the years. Those years are the projection, and the
MSE sets their catch from the control rule; the likelihood never scored
them either, since it fits only `Year <= endyr`. This is what a workbook
looks like when `endyr` has fallen behind the catch series, catch
through 2023 with `endyr` still 2019, and it is worth resolving before
running the MSE, because conditioning the assessment on those years is a
different question from projecting over them.

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
