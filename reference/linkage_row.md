# Build a single linkage row

Returns a one-row linkage table with defaults for the optional columns,
for tests and for assembling a table row by row.

## Usage

``` r
linkage_row(
  process,
  param,
  X_col,
  species = NA_integer_,
  sex = NA_integer_,
  age_bin = NA_integer_,
  fleet = NA_integer_,
  design_col = NA_character_,
  link = "identity",
  init = 0,
  init_supplied = FALSE,
  lower = -Inf,
  upper = Inf,
  prior_family = "none",
  prior_p1 = NA_real_,
  prior_p2 = NA_real_,
  re_group = NA_character_,
  re_struct = NA_character_,
  est_phase = 1L,
  re_index = NA_integer_,
  sigma_index = NA_integer_,
  re_time = NA_real_,
  re_sigma_init = NA_real_,
  re_sigma_prior_family = NA_character_,
  re_sigma_prior_p1 = NA_real_,
  re_sigma_prior_p2 = NA_real_,
  re_rho_init = NA_real_,
  re_rho_prior_family = NA_character_,
  re_rho_prior_p1 = NA_real_,
  re_rho_prior_p2 = NA_real_,
  re_obs_value = NA_real_,
  re_obs_sd = NA_real_,
  re_obs_est = NA,
  re_integrate = NA,
  re_pos = NA_integer_
)
```

## Arguments

- process, param, X_col:

  required identifying fields.

- species, sex, age_bin, fleet:

  stratum ids, `NA` = every level; `fleet` is a 1-based `Fleet_code`,
  used by catchability, selectivity and composition-weight linkages.

- design_col:

  name of the design matrix column.

- link:

  link function; one of
  [LINKAGE_LINKS](https://afsc-assessments.github.io/Rceattle/reference/LINKAGE_PROCESSES.md).

- init:

  initial value (default `0`).

- init_supplied:

  `TRUE` if the user supplied `init` (default `FALSE`).

- lower, upper:

  bounds (default `-Inf`, `Inf`).

- prior_family:

  one of
  [PRIOR_FAMILIES](https://afsc-assessments.github.io/Rceattle/reference/PRIOR_FAMILIES.md).
  `"none"` = no prior.

- prior_p1, prior_p2:

  family-specific prior parameters; ignored when
  `prior_family == "none"`.

- re_group:

  random-effect grouping label; `NA` = fixed.

- re_struct:

  random-effect covariance structure (`"us"`/`"rw"`/`"ar1"`); `NA` =
  fixed.

- est_phase:

  estimation phase ordinal; `0` = fix at `init`.

- re_index, sigma_index, re_time:

  random-effect fields filled by `pool_linkages()` (deviation slot, SD
  slot, time order); `NA` on fixed rows.

- re_sigma_init, re_sigma_prior_family, re_sigma_prior_p1,
  re_sigma_prior_p2:

  the group's natural-scale SD start value (held fixed if given without
  a prior) and prior, from `linkage_spec(init = , priors = )`; `NA` on
  fixed rows.

- re_rho_init, re_rho_prior_family, re_rho_prior_p1, re_rho_prior_p2:

  the `ar1` group's correlation on `(-1, 1)`, start value (held fixed if
  given without a prior) and prior, from
  [`linkage_spec()`](https://afsc-assessments.github.io/Rceattle/reference/linkage_spec.md);
  `NA` otherwise.

- re_obs_value, re_obs_sd:

  the Rogers QAR1 observed covariate value and measurement SD, from
  `linkage_spec(observe = , obs_sd = )`; `NA` when unobserved.

- re_obs_est:

  `TRUE` to estimate the QAR1 measurement SD
  (`linkage_spec(obs_sd_est = )`); `FALSE`/`NA` holds it fixed.

- re_integrate:

  `FALSE` for a penalized fixed effect rather than a Laplace-integrated
  random effect, from `linkage_spec(integrate = )`; `TRUE`/`NA` =
  integrated.

- re_pos:

  0-based position of the deviation within `beta_linkage_re` or
  `beta_linkage_re_pen`, whichever holds it.

## Value

A one-row `Rceattle_linkage_table`.
