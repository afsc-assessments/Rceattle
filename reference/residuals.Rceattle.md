# Residuals from an Rceattle fit

Returns a long-format data frame of residuals, following the convention
of
[`stats::residuals.glm()`](https://rdrr.io/r/stats/glm.summaries.html)
where `type` selects the *kind* of residual. By default residuals are
returned for every applicable data source, survey indices, fishery
catches, age/length composition (`comp`), and conditional age-at-length
(`caal`); use `source` to restrict to particular ones.

## Usage

``` r
# S3 method for class 'Rceattle'
residuals(
  object,
  type = "response",
  source = "all",
  scale = "log",
  species = NULL,
  ...
)
```

## Arguments

- object:

  An object of class `"Rceattle"` returned by
  [`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md).

- type:

  Residual kind: one of `"response"` (default), `"pearson"`, `"osa"`, or
  `"process"`.

- source:

  Data source(s) to include: any of `"index"`, `"catch"`, `"comp"`,
  `"caal"`, or `"all"` (default). `"diet"` (predator stomach-content
  composition) is also accepted but, because it uses a predator/prey
  schema rather than the fleet/bin layout, must be requested on its own.
  Ignored when `type = "process"`.

- scale:

  `"log"` (default) or `"natural"`. Only affects `"response"` residuals
  for `index` / `catch`.

- species:

  Optional species code(s) to include (matched against the `Species`
  column). Default `NULL` keeps all species. Mirrors the `species`
  argument of
  [`plot.rceattle_osa()`](https://afsc-assessments.github.io/Rceattle/reference/plot.rceattle_osa.md).

- ...:

  Passed to
  [`osa_residuals()`](https://afsc-assessments.github.io/Rceattle/reference/osa_residuals.md)
  (e.g. `method`, `seed`) when `type = "osa"`, or to
  [`process_residuals()`](https://afsc-assessments.github.io/Rceattle/reference/process_residuals.md)
  when `type = "process"`.

## Value

A `data.frame` with columns `Source`, `Fleet_code`, `Fleet_name`,
`Species`, `Sex`, `Year`, `Length`, `Bin`, `Age0_Length1`,
`Sample_size`, `Accumulated`, `Observed`, `Fitted`, `Residual`. Columns
are `NA` where they do not apply; `Accumulated` is `FALSE` except on a
composition bin that absorbed a folded tail.

## Details

Residual kinds (`type`):

- `"response"`:

  Observed minus fitted. For `index` / `catch` this is on the log scale
  by default (matching the lognormal likelihood; set `scale = "natural"`
  for the arithmetic difference); for `comp` / `caal` it is the
  difference in proportions, observed minus fitted.

- `"pearson"`:

  Standardized residuals. For `catch` and a log-scale survey index,
  \\(\log o - (\log\hat{o} - b\\\sigma^2/2))/\sigma\\ using the model's
  realized observation log-SD \\\sigma\\ and the observation
  bias-adjustment flag \\b\\ (`bias_adjust_obs`, default 1). A
  natural-scale index fleet (`Index_distribution` `"MVN"`, `"MVNORM"`,
  `"Normal"` or `"TruncatedNormal"`) holds an ABSOLUTE \\\sigma\\ and is
  standardized as \\(o - \hat{o})/\sigma\\ instead. Two caveats there: a
  covariance fleet gets its marginal residual, not the whitened one (use
  `type = "osa"` for that), and `"TruncatedNormal"` is standardized on
  the untruncated moments, so it is approximate where \\\hat{o}/\sigma\\
  is small enough for truncation to shift the mean. For `comp` / `caal`,
  \\(p - \hat{p})/\sqrt{\hat{p}(1 - \hat{p})/N}\\ with input sample size
  N.

- `"osa"`:

  One-step-ahead residuals via
  [`osa_residuals()`](https://afsc-assessments.github.io/Rceattle/reference/osa_residuals.md),
  which builds the composition observation data on demand from any fit.
  Pass its arguments through `...`; note that `method = "cdf"` forms no
  conditional mean, so `Fitted` is `NA` for every row under it.

- `"process"`:

  Process residuals via
  [`process_residuals()`](https://afsc-assessments.github.io/Rceattle/reference/process_residuals.md)
  for the model's random-effect deviations; `source` does not apply.

Composition rows are returned in long form (one row per observation x
age/length bin) and hold the `Age0_Length1` flag from `comp_data` (`0`
age, `1` length); CAAL rows hold both the conditioning `Length` and the
age `Bin`.

Where a fleet uses tail accumulation (`Comp_accum_young` /
`Comp_accum_old`), composition residuals describe the bins the
likelihood actually fit: the tails are folded into their boundary bin
and the bins outside the window are not reported, because the model
never fit them separately. `Bin` names the age or length each residual
belongs to and `Accumulated` marks the boundary bins, which stand for a
range rather than the single bin they are named for. Fleets without
accumulation are unchanged.
