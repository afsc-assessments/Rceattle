# Build the flat observation vector and metadata for OSA residuals

Internal helper used by
[`rearrange_data()`](https://afsc-assessments.github.io/Rceattle/reference/rearrange_data.md)
to assemble the inputs the TMB template needs to compute one-step-ahead
(OSA) residuals via
[`TMB::oneStepPredict()`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html).

OSA residuals require every observation that enters the likelihood to
live in a single flat vector (`obsvec`) with a companion `keep`
indicator. This function walks the already-built control/observation
matrices and produces:

- `obsvec`, the flat vector of observations. Aggregate (catch and index)
  observations are stored as log(observation) because their likelihood
  is lognormal; composition (comp) and conditional-age-at-length (caal)
  observations are stored as bin counts, `(proportion + 1e-5) * N`,
  matching exactly what the TMB likelihood forms during fitting.

- `obs_ctl`, a data frame with one row per element of `obsvec`, mapping
  each position back to its data type, fleet, species, year, age/length
  bin, etc., so residuals stay interpretable. R-side metadata only;
  removed before the data list is passed to TMB.

- `*_obsvec_idx`, per observation-row index vectors. For aggregate
  series this is the 0-based `obsvec` position of each row; for
  composition and caal it is the 0-based `obsvec` position of the row's
  FIRST bin (the template reads the row's bins as
  `obsvec.segment(start, n_bins)`). `-1` marks rows excluded from the
  likelihood.

The inclusion rules and per-row bin counts below must match the guards
in the TMB model exactly, so that every observation the model evaluates
has a valid `obsvec` position. All fitted observation types are handled:
aggregate catch and index, comp and caal composition, and predator diet
composition.

## Usage

``` r
build_osa_data(data_list, build_osa = FALSE)
```

## Arguments

- data_list:

  A `data_list` whose `*_ctl`/`*_obs`/`*_n` matrices have already been
  built (and coerced to matrices) by
  [`rearrange_data()`](https://afsc-assessments.github.io/Rceattle/reference/rearrange_data.md).

- build_osa:

  Logical. When `TRUE`, build the full OSA observation data for every
  type (aggregate, composition, caal, and diet) so
  [`osa_residuals()`](https://afsc-assessments.github.io/Rceattle/reference/osa_residuals.md)
  can be computed. When `FALSE` (the default), only the aggregate
  index/catch entries the TMB template always reads are built and the
  (much larger) composition/caal/diet metadata is skipped, this is the
  fast path for simulation testing (e.g.
  [`run_mse()`](https://afsc-assessments.github.io/Rceattle/reference/run_mse.md)),
  where the fitted objective is identical but OSA composition residuals
  are not produced.

## Value

The input `data_list` with `obsvec`, `obs_ctl`, `osa_mode`,
`comp_offset`, and the per-type `*_obsvec_idx` vectors added.

## Details

The composition proportion offset is read from `data_list$comp_offset`
(defaulting to `1e-5`), so the comp/caal bin counts are
`(proportion + comp_offset) * Neff`, the same offset the TMB likelihood
applies when fitting. Set it via `fit_control(comp_offset = )` or on
`data_list` directly.
