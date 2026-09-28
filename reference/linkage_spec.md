# Capture a linkage specification

Capture a linkage specification

## Usage

``` r
linkage_spec(
  formula,
  param = NULL,
  data = NULL,
  by = ~species,
  species = NULL,
  sex = NULL,
  fleet = NULL,
  link = "log",
  init = NULL,
  bounds = NULL,
  priors = NULL,
  re_group = NA_character_,
  est_phase = 1L,
  observe = NULL,
  obs_sd = NULL,
  obs_sd_est = FALSE,
  integrate = TRUE
)
```

## Arguments

- formula:

  one-sided R formula whose RHS describes the linear predictor for
  `param` (e.g. `~ 1`, `~ temp`, `~ temp + PDO`).

- param:

  target parameter name on the natural scale (e.g. `"alpha"`, `"M1"`,
  `"K"`). May be `NULL` when the spec is built inside a `build_*()` call
  that infers the parameter name from the enclosing list key (see
  [`build_growth()`](https://afsc-assessments.github.io/Rceattle/reference/build_growth.md)).

- data:

  Optional data frame for formula validation; validation is performed at
  fit time inside
  [`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md).

- by:

  one-sided formula naming stratifying factors that should each get
  their own coefficients. Allowed names are `species`, `sex`, `age_bin`,
  and `fleet` (`fleet` for catchability and selectivity linkages).
  **When omitted, `by` defaults to the base stratum of whichever process
  the spec is attached to**, `~fleet` for catchability, selectivity, and
  the fleet composition weights (`theta_comp` / `theta_caal`), and
  `~species` for recruitment, M, growth, and the diet weight
  (`theta_diet`), so you rarely need to spell it out for the base case.
  Pass it explicitly to override: e.g. `~species + sex` for
  per-(species, sex) coefficients, or `NULL` to share a single
  coefficient across every stratum. An explicit `by` (including `NULL`)
  is always kept as given.

- species:

  optional vector of species that this spec applies to, given either as
  1-based species ids (`c(1L, 2L)`) or as species **names** matching
  `data_list$spnames` (`c("Pollock", "Cod")`). Names are matched
  exactly, after trimming whitespace, when the model is assembled in
  [`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md);
  an unrecognized name is an error that lists the model's species. Give
  ids or names, not a mix, R coerces `c(1, "Cod")` to `c("1", "Cod")`.
  `NULL` (default) means every species in `strata$species` at
  materialization time. Use this to give different species different
  formulas, e.g. by registering multiple specs against the same
  parameter, see
  [`build_growth()`](https://afsc-assessments.github.io/Rceattle/reference/build_growth.md)
  for the multi-spec syntax.

- sex:

  optional vector of sex ids that this spec applies to. May be supplied
  as integers (`1L` = female, `2L` = male) or as character strings
  (`"Females"`/`"Males"`, case-insensitive; `"female"`, `"male"`, `"f"`,
  `"m"` are also accepted). `NULL` (default) means every sex in
  `strata$sex` at materialization time. `by` must include `sex` for it
  to apply; otherwise it warns and has no effect. Use this to register
  separate specs per sex (e.g. one prior on females, another on males)
  against the same parameter.

- fleet:

  optional vector of fleets this spec applies to, given either as
  1-based `Fleet_code`s (`c(1L, 3L)`) or as fleet **names** matching
  `fleet_control$Fleet_name` (`c("Shelikof", "Summer BT")`). Names are
  matched exactly, after trimming whitespace, when the model is
  assembled in
  [`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md);
  an unrecognized name, or one that is not unique in `fleet_control`, is
  an error that lists the model's fleets. Prefer names: a `Fleet_code`
  that is wrong but in range attaches the linkage to a different fleet
  and still fits, whereas a misspelled name cannot. Give ids or names,
  not a mix, R coerces `c(7, "Pollock")` to `c("7", "Pollock")`. `NULL`
  (default) means every fleet in `strata$fleet` at materialization time.
  `by` must include `fleet` for it to apply; otherwise it warns and has
  no effect. Used by catchability and selectivity linkages to give
  different fleets different formulas.

- link:

  link function relating the linear predictor to the natural-scale
  target parameter. One of `"log"` (default) or `"identity"`. With
  `link = "log"`, `log(param) = X * beta`, slope contributions are
  multiplicative on the natural-scale parameter. With
  `link = "identity"`, `param = X * beta`, slope contributions are
  additive on the natural scale. The linkage targets are estimated on
  the log scale, so `"log"` is the default.

- init:

  optional named list of initial values keyed by the design-matrix
  column name (e.g. `` list(`(Intercept)` = 4, temp = 0) ``), an
  intercept on the parameter's natural scale and a slope on the link
  scale; missing entries default to `0`.

- bounds:

  optional named list of `c(lower, upper)` keyed the same way as `init`.

- priors:

  optional named list of
  [Rceattle_priors](https://afsc-assessments.github.io/Rceattle/reference/Rceattle_priors.md)
  objects, keyed by design-matrix column name. Inside this argument you
  may write `normal()`, `lognormal()`,
  [`gamma()`](https://rdrr.io/r/base/Special.html), or
  [`beta()`](https://rdrr.io/r/base/Special.html) directly, e.g.
  `priors = list(temp = normal(0, 1))`, equivalent to
  `priors = list(temp = prior_normal(0, 1))`.

- re_group:

  optional character: name of a random-effect grouping for these
  coefficients. `NA` (default) means fixed.

- est_phase:

  optional integer estimation phase. Default `1L`; `0` fixes the
  coefficient at its `init`, which then holds over any `inits` given to
  [`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)
  (an estimated intercept's `init` is only a starting value, and `inits`
  win). Applies to **fixed-effect** rows only, the coefficients in
  `beta_linkage`. A random-effect term's deviations are held in a
  separate vector that `est_phase` does not reach, so `est_phase < 1` on
  a formula containing one is an error rather than a silent no-op; drop
  the term, or fix a small SD via `init = list(sigma = )`, to remove the
  time variation. Values above `1` are currently inert for every linkage
  row.

- observe:

  optional character: for an `ar1(1 | group)` term, the name of an
  `env_data` column that measures the AR1 latent (a state-space
  covariate, sensu Rogers et al. 2024). The latent enters the linked
  parameter through an estimated effect size and is observed against
  this column. `NULL` (default) leaves the AR1 as a plain random effect.
  The latent is zero-mean (no estimated level), so standardize the
  observed covariate to mean 0; a non-zero-mean column confounds its
  level with the intercept and warns.

- obs_sd:

  optional positive numeric: the measurement SD for the `observe`
  covariate (one per observed group). Required with `observe`, unused
  otherwise. Held **fixed** at this value by default
  (`obs_sd_est = FALSE`); it is the *starting* value when
  `obs_sd_est = TRUE`.

- obs_sd_est:

  optional single `TRUE`/`FALSE` (default `FALSE`): estimate the
  `observe` measurement SD instead of holding it fixed, as the
  state-space survey-catchability (GOA pollock) model does. **Caveat:**
  the effect size and `obs_sd` are only jointly identified when the
  observed covariate is informative; on a smooth series the AR1 latent
  can track it exactly and the freely-estimated `obs_sd` collapses
  toward 0. Keep it fixed unless the covariate is informative. Only used
  with `observe`.

- integrate:

  single `TRUE`/`FALSE` (default `TRUE`): whether the random effect's
  deviations are integrated out by the Laplace approximation.
  `integrate = FALSE` instead estimates them as a **penalized fixed
  effect**, the deviations stay in the objective as a plain penalty and
  are reported with standard errors like any other fixed effect. This
  reproduces the ADMB/AMAK convention behind the legacy
  `Time_varying_sel` / `Time_varying_q` switches, which a
  Laplace-integrated `rw()` cannot match (the marginal likelihood holds
  a log-determinant term the penalized form has no counterpart for).
  Permitted **only with a fixed SD**, `init = list(sigma = )` and no
  `sigma` prior, plus a fixed `rho` for `ar1`, because estimating
  deviations and their SD jointly as fixed effects is degenerate. Cannot
  be combined with `observe`: an observed latent state must stay
  integrated.

## Value

An `Rceattle_linkage_spec` object.

## Details

The reserved keys `sigma` and `rho` in `init` / `priors` route the
random-effect deviation SD and (for `ar1`) the correlation: e.g.
`init = list(sigma = 0.1)` fixes the SD,
`priors = list(rho = normal(0, 0.3))` places a prior on the correlation.
`sigma` means different things by structure: for `rw()` it is the
innovation (per-step) SD; for `ar1()` it is the marginal (stationary)
SD. The two are not directly comparable across structures, see
[`vignette("environmental-linkages-and-priors")`](https://afsc-assessments.github.io/Rceattle/articles/environmental-linkages-and-priors.md).

## Examples

``` r
# A covariate on a growth parameter, one coefficient per species. Attach to a
# parametric build_growth(); the default "empirical" form consumes no linkages.
linkage_spec(~ temp, param = "K")
#> <Rceattle linkage spec>
#>   param:   K
#>   formula: ~temp
#>   link:    log

# A covariate on catchability for one fleet, named rather than numbered.
# `by` defaults to ~ fleet for catchability, so it need not be given.
linkage_spec(~ temp, fleet = "Pollock_survey_1_shelikof_acoustic")
#> <Rceattle linkage spec>
#>   param:   NA
#>   formula: ~temp
#>   fleet:   Pollock_survey_1_shelikof_acoustic
#>   link:    log

# A separate coefficient per year, the fixed-effect alternative to rw()/ar1().
linkage_spec(~ factor(Year), fleet = "Pollock_survey_1_shelikof_acoustic")
#> <Rceattle linkage spec>
#>   param:   NA
#>   formula: ~factor(Year)
#>   fleet:   Pollock_survey_1_shelikof_acoustic
#>   link:    log

# IID annual deviations; rw() for a random walk, ar1() for AR1.
linkage_spec(~ (1 | Year))
#> <Rceattle linkage spec>
#>   param:   NA
#>   formula: ~(1 | Year)
#>   link:    log

# A random walk estimated as a penalized fixed effect, which requires a
# fixed SD, the ADMB/AMAK convention behind the legacy Time_varying_*
# switches.
linkage_spec(~ rw(1 | Year), init = list(sigma = 0.05), integrate = FALSE)
#> <Rceattle linkage spec>
#>   param:   NA
#>   formula: ~rw(1 | Year)
#>   link:    log
#>   estimate: penalized fixed effect (integrate = FALSE)

# Intercept-only: a prior on the base parameter itself. The bare normal()
# resolves to prior_normal() inside `priors`.
linkage_spec(~ 1, priors = list(`(Intercept)` = normal(0, 1)))
#> <Rceattle linkage spec>
#>   param:   NA
#>   formula: ~1
#>   prior:    (Intercept) ~ normal(0, 1)
#>   link:    log
```
