# Build a single linkage row

Convenience constructor that returns a one-row linkage table with
default values for optional columns. Useful in tests and for incremental
table assembly.

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

  stratum ids; `NA` = shared across the dimension. `fleet` is a 1-based
  `Fleet_code`, used by catchability and selectivity linkages; the
  process-level linkages leave it `NA`.

- design_col:

  name of the design matrix column.

- link:

  link function; one of
  [LINKAGE_LINKS](https://afsc-assessments.github.io/Rceattle/reference/LINKAGE_PROCESSES.md).

- init:

  initial value (default `0`).

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

  random-effect registry fields filled by `pool_linkages()`; `NA` on
  fixed rows. `re_index` is the 0-based slot in `beta_linkage_re`,
  `sigma_index` the 0-based slot in `log_sigma_linkage`, and `re_time`
  the numeric grouping value used to order `rw()`/`ar1()` deviations in
  real elapsed time.

- re_sigma_init, re_sigma_prior_family, re_sigma_prior_p1,
  re_sigma_prior_p2:

  per-group RE-SD routing from
  `linkage_spec(init = list(sigma = ), priors = list(sigma = ))`;
  identical across a group's rows, `NA` on fixed rows. `re_sigma_init`
  is the start (or, when supplied without a prior, fixed) SD on the
  natural scale; the prior triple places a prior on that SD.

- re_rho_init, re_rho_prior_family, re_rho_prior_p1, re_rho_prior_p2:

  per-group `ar1` correlation routing from
  `linkage_spec(init = list(rho = ), priors = list(rho = ))`; natural
  `(-1, 1)` scale, `NA` on non-`ar1` rows.

- re_obs_value, re_obs_sd:

  state-space (Rogers QAR1) observation from
  `linkage_spec(observe = , obs_sd = )`: the observed covariate value at
  this row's time and the fixed measurement SD. `NA` when the group is
  unobserved.

- re_obs_est:

  `TRUE` to estimate the QAR1 measurement SD; `FALSE`/`NA` holds it
  fixed.

- re_integrate:

  `FALSE` when the deviations are estimated as a penalized fixed effect
  rather than integrated out by the Laplace approximation, from
  `linkage_spec(integrate = FALSE)`. `TRUE`/`NA` = integrated.

- re_pos:

  0-based position of this row's deviation within the parameter vector
  that holds it (`beta_linkage_re` when `re_integrate`, else
  `beta_linkage_re_pen`). Distinct from `re_index`, which is the global
  slot.

## Value

A one-row `Rceattle_linkage_table`.
