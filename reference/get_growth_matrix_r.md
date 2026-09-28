# Generate Length-at-Age Transition Matrix

This function calculates a probability transition matrix that defines
the probability of a fish of a given age belonging to specific length
bins. It supports Von Bertalanffy and Richards growth models and
includes a Stock Synthesis (SS) style plus-group correction.

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

  Numeric. Fraction of the year (0 = Jan 1st).

- nsex_sp:

  Integer. Number of sexes for the species.

- nages_sp:

  Integer. Number of age classes.

- nlengths_sp:

  Integer. Number of length bins.

- nyrs:

  Integer. Number of years in the simulation.

- lengths_sp:

  Vector. Boundaries of the length bins.

- minage_sp:

  Numeric. The reference age (L1) for growth estimation.

- maxage_sp:

  Numeric. The age at which growth enters the asymptotic phase.

- growth_params_sp:

  Array. Dimensions (sex, yr, 4). Params: K, L1, Linf, Richards m.

- growth_log_sd_sp:

  Array. Dimensions (sex, 2). Log-SD of length: 1st param is SD at
  minage, 2nd param is SD at maxage.

- growth_model_sp:

  Integer. 1 = Von Bertalanffy, 2 = Richards.

- M_plus_sp:

  Numeric. Natural mortality at the oldest age (one value, or one per
  sex), which weights the ages pooled in the plus group by survival, as
  the template does.

## Value

A 4D array of probabilities with dimensions (sex, age, length, year).
