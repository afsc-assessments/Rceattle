# Tidy long-format derived quantities from an Rceattle fit

Returns the model's derived population quantities in long form so that
custom plots and post-processing don't have to walk the deeply nested
`quantities` list or inherit the dimnames decisions in
[`rename_output()`](https://afsc-assessments.github.io/Rceattle/reference/rename_output.md).
Two shapes are flattened into one tidy frame: species-by-year quantities
(e.g. `biomass`, `ssb`, `R`, `F_spp`) and species-by-sex-by-age-by-year
quantities (e.g. `N_at_age`, `biomass_at_age`, `M_at_age`). For the
species-level shapes, `sex` and `age` are returned as `NA`. Cells of the
4D arrays that are padded out to `max(nsex)` / `max(nages)` for species
with fewer sexes or ages are dropped.

## Usage

``` r
# S3 method for class 'Rceattle'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  which = c("biomass", "ssb", "R", "biomass_depletion", "ssb_depletion", "F_spp"),
  ci_level = 0.95,
  ...
)
```

## Arguments

- x:

  An object of class `"Rceattle"` returned by
  [`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md).

- row.names, optional:

  Ignored; present for the
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
  generic.

- which:

  Character vector of quantity names to extract, or `"all"` for every
  quantity with a known shape. Defaults to a common population-level
  summary. See `names(Rceattle:::.RCEATTLE_QUANTITIES)` for the full
  list.

- ci_level:

  Confidence level for `lwr` and `upr`. Default `0.95`.

- ...:

  Currently unused.

## Value

A `data.frame` with columns `year`, `species`, `sex`, `age`, `quantity`,
`value`, `se`, `lwr`, `upr`. `species` is the character species name
from `data_list$spnames`. Rows are sorted in the order `which` was
given.

## Details

Standard errors (`se`) and confidence intervals (`lwr`, `upr`) are
populated from the TMB `sdreport` for any quantity that was `ADREPORT`'d
(`biomass`, `ssb`, `R`, `exploitable_biomass`, `biomass_depletion`,
`ssb_depletion`); other quantities and fits produced with
`getsd = FALSE` get `NA` for `se` / `lwr` / `upr`. Set `ci_level` to
widen or narrow the band.

A strictly positive series takes its interval on the log scale, so it is
right-skewed and cannot reach zero. This is the same construction
[`plot_timeseries()`](https://afsc-assessments.github.io/Rceattle/reference/plot_timeseries.md)
draws, so the table and the figure agree.
