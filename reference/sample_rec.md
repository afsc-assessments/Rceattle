# Sample historical recruitment deviations into the projection

Sample historical recruitment deviations into the projection

## Usage

``` r
sample_rec(
  object = NULL,
  sample_rec = TRUE,
  update_model = TRUE,
  rec_trend = 0,
  Rceattle = NULL
)
```

## Arguments

- object:

  CEATTLE model object exported from `Rceattle`

- sample_rec:

  Include resampled recruitment deviations from the hindcast in the OM
  projection. Resampled deviations are used rather than drawing from
  N(0, sigmaR) because the initial deviations bias R0 low. If FALSE,
  uses one deviation (plus the log of the `rec_trend` multiplier) for
  every projection year: `log(mean(R)) - log(R0)` under mean
  recruitment, and `log(mean(R / R_hat))` under a curve, over every
  hindcast year when the curve is fitted in the hindcast and over the
  penalty years for the penalty form.

- update_model:

  Update model dynamics. Default = TRUE

- rec_trend:

  Linear increase or decrease in mean recruitment from `endyr` to
  `projyr`. This is the terminal multiplier
  `mean rec * (1 + (rec_trend/projection years) * 1:projection years)`.
  Can be of length 1 or of length nspp. If length 1, all species get the
  same trend.

- Rceattle:

  deprecated name for `object`, still accepted so existing scripts keep
  working. Supplying both is an error.

## Value

Rceattle model

## Details

Projected recruitment uses the sampled deviations only if the model was
fitted with `build_srr(proj_mean_rec = FALSE)`; under the default it is
mean recruitment over the hindcast, and `update_model = TRUE` warns. In
[`run_mse()`](https://afsc-assessments.github.io/Rceattle/reference/run_mse.md)
each projection year becomes an operating-model hindcast year, where the
deviations apply either way.

`update_model = TRUE` holds the hindcast estimates fixed. A fit run with
its harvest control rule (`estimateMode` `"Estimate"` or `"Projection"`)
reruns that rule, so an F the rule estimates is re-solved for the new
recruitment. Under predation (`msmMode > 0`) a fishing rule's unfished
SSB and biomass (`MSSB0`, `MSB0`) are re-derived too, so SB0, B0 and
depletion in every year follow the new recruitment, as in a fresh
[`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)
from those deviations.
