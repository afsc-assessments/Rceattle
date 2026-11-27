# Generate Length-at-Age Transition Matrix

Calculates the age-to-length transition matrix under von Bertalanffy or
Richards growth with a Stock Synthesis plus-group correction; like
`growth.hpp`, adapted from Giancarlo Correa's growth code in WHAM's
growth branch.

## Usage

``` r
get_growth_matrix_r(
  fracyr,
  nsex_sp,
  nages_sp,
  nlengths_sp,
  nyrs,
  lengths_sp,
  minage_sp,
  maxage_sp,
  growth_params_sp,
  growth_log_sd_sp,
  growth_model_sp,
  M_plus_sp
)
```

## Arguments

- fracyr:

  Numeric fraction of the year (`0` = Jan 1).

- nsex_sp:

  Integer number of sexes for the species.

- nages_sp:

  Integer number of age classes.

- nlengths_sp:

  Integer number of length bins.

- nyrs:

  Integer number of years in the simulation.

- lengths_sp:

  Numeric vector of the `nlengths_sp` length-bin lower edges (cm); the
  first bin is open below and the last open above.

- minage_sp:

  Numeric reference age (`L1`) for growth estimation.

- maxage_sp:

  Unused by the body; `helpers-make-msm-data.R` passes it.

- growth_params_sp:

  Numeric `[sex, yr, 4]` array of `K`, `L1`, `Linf`, and Richards `m`.

- growth_log_sd_sp:

  Numeric `[sex, 2]` array of log length SD at `minage` and `maxage`.

- growth_model_sp:

  Integer growth code (`1` = von Bertalanffy, `2` = Richards).

- M_plus_sp:

  Numeric natural mortality at the oldest age (scalar or per sex),
  weighting ages pooled in the plus group by survival to match the TMB
  growth model.

## Value

A list: `length_at_age`, mean length (cm) by (sex, age, year), and
`growth_matrix`, the probability of each length bin by (sex, age,
length, year).
