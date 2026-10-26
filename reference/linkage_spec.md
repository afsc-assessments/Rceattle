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
  `param` (e.g. `~ 1`, `~ temp`, `~ temp + (1 | Year)`), where a
  `(1 | g)`, `rw(1 | g)` or `ar1(1 | g)` term adds a random effect with
  one deviation per level of `g`.

- param:

  target parameter name on the natural scale (e.g. `"alpha"`, `"M1"`,
  `"K"`). May be `NULL` when the spec is built inside a `build_*()` call
  that infers it from the enclosing list key; see
  [`build_growth()`](https://afsc-assessments.github.io/Rceattle/reference/build_growth.md).

- data:

  Optional data frame for formula validation; validation happens at fit
  time inside
  [`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md).

- by:

  one-sided formula naming stratifying factors that each get their own
  coefficients: `species`, `sex`, `age_bin`, or `fleet`. **Omitted, it
  defaults to the base stratum of the process the spec is attached to**
  – `~fleet` for catchability, selectivity and the fleet composition
  weights, `~species` for recruitment, M, growth and the diet weight –
  so the base case needs no `by`. Pass `~species + sex` for
  per-(species, sex) coefficients, or `NULL` to share one coefficient
  across every stratum; an explicit `by`, `NULL` included, is kept as
  given. In a **selectivity or catchability** linkage under `~fleet`,
  fleets sharing a `Selectivity_index` or `Catchability_index` still
  take one coefficient between them, because they share one parameter
  block.

- species:

  optional species this spec applies to, as 1-based ids (`c(1L, 2L)`) or
  as names matching `data_list$spnames` (`c("Pollock", "Cod")`). Names
  are matched exactly after trimming whitespace, and an unrecognized one
  is an error listing the model's species. Give ids or names, not a mix,
  since R coerces `c(1, "Cod")` to character. `NULL` (default) means
  every species. Register several specs against one parameter to give
  species different formulas.

- sex:

  optional sex ids, as integers (`1L` female, `2L` male) or strings
  (`"Females"`/`"Males"`, case-insensitive; `"f"`/`"m"` accepted).
  `NULL` (default) means every sex. `by` must include `sex` or it warns
  and has no effect.

- fleet:

  optional fleets, as 1-based `Fleet_code`s or as names matching
  `fleet_control$Fleet_name`. **Prefer names:** a `Fleet_code` that is
  wrong but in range attaches the linkage to a different fleet and still
  fits, where a misspelled name cannot. Names are matched exactly after
  trimming; one that is unrecognized or not unique is an error listing
  the fleets. Give ids or names, not a mix. `NULL` (default) means every
  fleet, and `by` must include `fleet` or it warns and has no effect.

- link:

  how the linear predictor reaches the natural-scale parameter: `"log"`
  (default), so `log(param) = X * beta` and slopes are multiplicative;
  `"identity"`, so `param = X * beta` and slopes are additive; or
  `"exponential"`, where the predictor multiplies the log rather than
  shifting it, `param^exp(X * beta)`. The exponential form reproduces
  Stock Synthesis's environmental link type 1 and is accepted on
  catchability only, for a lognormal index, with an estimated intercept.

- init:

  optional named list of initial values keyed by design-matrix column
  (e.g. `` list(`(Intercept)` = 4, temp = 0) ``), an intercept on the
  parameter's natural scale and a slope on the link scale; missing
  entries default to `0`. An intercept's value is logged onto the
  parameter's base (`rec_pars`, `log_M1`, ...), or onto the coefficient
  itself for recruitment `R_init`, which has no base.

- bounds:

  optional named list of `c(lower, upper)`, keyed and scaled as `init`.

- priors:

  optional named list of
  [Rceattle_priors](https://afsc-assessments.github.io/Rceattle/reference/Rceattle_priors.md)
  objects keyed by design-matrix column. Inside this argument
  `normal()`, `lognormal()`,
  [`gamma()`](https://rdrr.io/r/base/Special.html) and
  [`beta()`](https://rdrr.io/r/base/Special.html) may be written
  directly, so `priors = list(temp = normal(0, 1))` is the same as
  `prior_normal(0, 1)`.

- est_phase:

  optional integer estimation phase, default `1L`. `0` fixes the
  coefficient at its `init`, which then holds over any `inits` given to
  [`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)
  (an estimated intercept's `init` is only a starting value, and `inits`
  win). It applies to the fixed-effect rows in `beta_linkage` only: a
  random-effect term's deviations live in a separate vector it cannot
  reach, so `est_phase < 1` on a formula containing one is an error
  rather than a silent no-op. Drop the term, or fix a small SD with
  `init = list(sigma = )`, to remove the time variation. Values above
  `1` are inert for every linkage row.

- observe:

  optional `env_data` column measuring the AR1 latent of an
  `ar1(1 | group)` term, making it a state-space covariate (sensu Rogers
  et al. 2024): the latent reaches the linked parameter through an
  estimated effect size and is observed against this column. `NULL`
  (default) leaves the AR1 a plain random effect. The latent is
  zero-mean, so standardize the column to mean 0; a non-zero-mean column
  confounds its level with the intercept and warns.

- obs_sd:

  optional positive measurement SD for the `observe` covariate, one per
  observed group. Required with `observe` and unused otherwise. Held
  fixed by default, or the starting value when `obs_sd_est = TRUE`.

- obs_sd_est:

  estimate the `observe` measurement SD instead of fixing it, as the
  state-space survey-catchability (GOA pollock) model does. Default
  `FALSE`. **Caveat:** the effect size and `obs_sd` are jointly
  identified only when the observed covariate is informative; on a
  smooth series the AR1 latent tracks it exactly and a free `obs_sd`
  collapses toward 0.

- integrate:

  whether the random effect's deviations are integrated out by the
  Laplace approximation, default `TRUE`. `FALSE` estimates them as a
  penalized fixed effect instead; see Details for when that is
  permitted.

## Value

An `Rceattle_linkage_spec` object.

## Details

The reserved keys `sigma` and `rho` in `init` / `priors` route the
random-effect deviation SD and, for `ar1`, the correlation. `sigma`
means different things by structure – the innovation (per-step) SD for
`rw()`, the marginal (stationary) SD for `ar1()` – so the two are not
comparable across structures; see
[`vignette("environmental-linkages-and-priors")`](https://afsc-assessments.github.io/Rceattle/articles/environmental-linkages-and-priors.md).

`integrate = FALSE` leaves the deviations in the objective as a plain
penalty and reports them with standard errors. That reproduces the
ADMB/AMAK convention behind the legacy `Time_varying_sel` /
`Time_varying_q` switches, which a Laplace-integrated `rw()` cannot
match, because the marginal likelihood carries a log-determinant term
the penalized form has no counterpart for. Permitted only with a fixed
SD (`init = list(sigma = )`, no `sigma` prior, plus a fixed `rho` for
`ar1`), since estimating deviations and their SD jointly as fixed
effects is degenerate, and never with `observe`, which must stay
integrated.

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

# IID annual deviations; rw() for a random walk, ar1() for AR1.
linkage_spec(~ (1 | Year))
#> <Rceattle linkage spec>
#>   param:   NA
#>   formula: ~(1 | Year)
#>   link:    log

# A random walk estimated as a penalized fixed effect, which requires a fixed
# SD: the ADMB/AMAK convention behind the legacy Time_varying_* switches.
linkage_spec(~ rw(1 | Year), init = list(sigma = 0.05), integrate = FALSE)
#> <Rceattle linkage spec>
#>   param:   NA
#>   formula: ~rw(1 | Year)
#>   link:    log
#>   estimate: penalized fixed effect (integrate = FALSE)

# Intercept-only: a prior on the base parameter itself.
linkage_spec(~ 1, priors = list(`(Intercept)` = normal(0, 1)))
#> <Rceattle linkage spec>
#>   param:   NA
#>   formula: ~1
#>   prior:    (Intercept) ~ normal(0, 1)
#>   link:    log
```
