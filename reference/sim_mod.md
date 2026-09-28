# Simulate Rceattle data

Simulates the data an Rceattle model would have produced, either as
expected values or as a random draw. Every observation type is covered:
survey biomass (under the fleet's own `Index_distribution`, lognormal,
natural-scale normal, or the correlated MVN/MVNORM draw from its
covariance), total catch (lognormal), age/length composition and
conditional age-at-length (multinomial or Dirichlet-multinomial), and
stomach contents (multinomial or Dirichlet-multinomial).

## Usage

``` r
sim_mod(object = NULL, simulate = FALSE, process = FALSE, Rceattle = NULL)
```

## Arguments

- object:

  A CEATTLE model object exported from `Rceattle`.

- simulate:

  Logical. If `TRUE`, simulates data from distributions. If `FALSE`,
  returns the expected values (hats).

- process:

  Which process error to redraw alongside the observations. `FALSE`
  (default) or `"none"` keeps the fitted deviations; `TRUE` or `"all"`
  redraws every process; `"dynamics"` covers recruitment, natural
  mortality and growth, `"observation"` covers catchability and
  selectivity. Any subset of `"recruitment"`, `"M"`, `"growth"`,
  `"catchability"` and `"selectivity"` may be given instead. Ignored
  when `simulate = FALSE`.

- Rceattle:

  deprecated name for `object`, still accepted so existing scripts keep
  working. Supplying both is an error.

## Value

A `data_list` object containing the simulated or expected data values,
formatted for use in `Rceattle`. When `process` redrew something, the
deviations that generated the data are attached as
`attr(x, "process_sim")`, a named list holding whichever of `rec_dev`,
`init_dev`, `log_M1_dev` and `beta_linkage_re` were drawn. Those are the
truth a refit has to recover; without them the only comparison available
is against the original fitted deviations, which are no longer the
values that generated the data.

Each is accompanied by a logical of the same shape named with a `_drawn`
suffix (`rec_dev_drawn`, ...), `TRUE` where the draw touched that cell.
The draws cover the hindcast only, and `beta_linkage_re` is one vector
over every random-linkage slot whether or not its process was asked for,
so the arrays hold fitted values alongside simulated ones. Restrict any
recovery statistic to the `_drawn` cells, over the full array it reports
perfect recovery on the cells that were never redrawn.

## Details

Every draw is taken by the TMB model itself, in a `SIMULATE` block
beside the likelihood that defines it, so the two are edited together. A
simulator that has drifted from its likelihood does not error, it makes
[`self_test`](https://afsc-assessments.github.io/Rceattle/reference/self_test.md)
report recovery against a process the likelihood never assumed.

Consequently `simulate = TRUE` needs a model to evaluate. A model loaded
from an `.Rdata`/`.rds` file has one, and a fit whose `$obj` was dropped
to save space is rebuilt from its `data_list` and estimates, provided
the rebuild reproduces the fit's own expected values.
[`model_average`](https://afsc-assessments.github.io/Rceattle/reference/model_average.md)
output cannot be simulated from at all: its quantities are an average
over models rather than any one model's fit, so no parameters produced
them. `simulate = FALSE` reads only `$quantities` and draws no random
numbers, so it works on any model.

Rows the model predicts nothing for are left as they are, and each is
reported by a warning: a composition for a fleet with no catch that year
comes back empty; a stomach whose predator has an empirical suitability,
and a covariance survey fleet's observations outside its fitted window,
keep their observed values. A composition, CAAL or diet row whose sample
size times its weight rounds below one observation also comes back
empty, and its `Sample_size` is dropped to zero so the refit does not
score a row with no data behind it.

One deliberate exception to draw-equals-density: compositions are drawn
from the predicted proportions before `comp_offset` (the constant added
to observed and predicted alike to keep `log(0)` out of the density,
default 1e-5). A bin the model predicts at zero is therefore drawn at
zero while the density scores it at `comp_offset`. The difference is
about 4e-4 relative on a 20-bin composition and does not bias the round
trip, since the offset is applied to both sides on the next fit.

Simulating leaves the drawn values in the object's report environment,
under names ending `_sim`. The estimates, the data and the objective
function are untouched, so a later
[`osa_residuals()`](https://afsc-assessments.github.io/Rceattle/reference/osa_residuals.md)
or [`vcov()`](https://rdrr.io/r/stats/vcov.html) on the same model is
unaffected.

## Examples

``` r
if (FALSE) { # \dontrun{
data(BS2017SS)
fit <- fit_mod(BS2017SS, estimateMode = "Hindcast")
# Expected values only, no observation error drawn.
sim_mod(fit)
# One stochastic replicate, drawn from each fleet's own likelihood.
sim_mod(fit, simulate = TRUE)
} # }
```
