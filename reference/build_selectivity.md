# Selectivity specification

Holds environmental linkages on selectivity parameters. An effect is
written as a formula and composes additively with any `Time_varying_sel`
process error on the same fleet, the two being separate mechanisms: a
covariate effect against a deviation.

The `linkages` keys are the shape parameters of the parametric forms, so
which apply depends on the fleet's `Selectivity`:

- `slp_asc`, `slp_desc`:

  ascending / descending logistic slope (log scale); the double-normal
  widths, aliased `sigma_asc` / `sigma_desc`.

- `inf_asc`, `inf_desc`:

  ascending / descending inflection age or length (natural scale); the
  double-normal peak and logit right-floor, aliased `peak` /
  `right_floor`.

- `coff`:

  non-parametric selectivity-at-bin coefficients.

- `apical`:

  a log-scale multiplier on one sex's whole curve, applied after the
  form and before normalization, so every estimated form takes it. See
  Details.

`DoubleNormalSS3` (code 15) carries its own six under Stock Synthesis's
names – `dn_peak`, `top_logit`, `ascend_se`, `descend_se`,
`start_logit`, `end_logit` – and a linkage is the **only** way to vary
that form over time, since `Time_varying_sel` is refused on it. A time
block is `~ cut(Year, ...)`; annual deviations are a random-effect term.

Every parameter but `apical` takes `link = "log"` (multiplicative on the
natural parameter) or `link = "identity"` (additive); `apical` is
already a multiplier and takes `"log"` only.

## Usage

``` r
build_selectivity(linkages = NULL)
```

## Arguments

- linkages:

  Optional named list of
  [`linkage_spec()`](https://afsc-assessments.github.io/Rceattle/reference/linkage_spec.md)
  objects keyed by selectivity parameter. Coefficients are per fleet by
  default (`by = ~ fleet`), except that fleets sharing a
  `Selectivity_index` estimate one parameter block and so take ONE
  coefficient for the whole group, which is what makes them mirror each
  other. Use
  [`linkage_spec()`](https://afsc-assessments.github.io/Rceattle/reference/linkage_spec.md)'s
  `fleet` argument to restrict a spec to particular fleets.

## Value

A list of selectivity settings for
[`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md).

## Details

**The `apical` offset.** Fishing mortality is one `log_F` per fleet and
year shared by the sexes, and no form has a height parameter, so a sex
difference in F can only come from selectivity. `apical` multiplies one
sex's curve by `exp(log_sel_apical)` bin by bin, as Stock Synthesis's
male-offset option does. Name the fleet and the sex that holds it
(`fleet = 3`, `sex = "male"`; `by` defaults to `~ fleet + sex` here) and
the other sex is the reference. It equals the ratio of the sexes' peak
heights only where their shapes peak equally, so for a sex-specific dome
read it as the multiplier and take the peak ratio from
`fit$quantities$sel_at_age`.

Only the contrast between sexes is identified, the common level being
`log_F`, so one sex holds it and the fit is refused if both do, if no
fleet or sex is named, if the species has one sex, on a `Fixed`, AR1 or
`Fleet_type = "Off"` fleet, on a fleet sharing another's
`Selectivity_index`, under `link = "identity"` (which could drive the
multiplier negative), or where `Sel_norm_scope = "WithinSex"` would
divide it straight back out. The contrast is informed only by joint
composition (`comp_data$Sex = 3`); with single-sex compositions it rests
on its prior, and
[`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)
warns when a fleet has neither. An intercept prior is on the
multiplier's natural scale, so a `lognormal()` centred on 1 means no
offset. See
[`vignette("environmental-linkages-and-priors")`](https://afsc-assessments.github.io/Rceattle/articles/environmental-linkages-and-priors.md).

**Priors on a selectivity parameter.** An intercept-only formula (`~ 1`)
with a `priors` entry puts a prior on the parameter itself, adding no
year-to-year offset. Read it on the parameter's own scale: the slopes
are log scale (`lognormal()`), the inflections natural scale
(`normal()`). A prior targets one parameter, so in a two-sex model an
unstratified `~ 1` prior constrains sex 1 only – use `by = ~ sex`. An
`init` on a selectivity intercept has no effect, since the starting
value comes from the data, and a prior on the double-normal
`right_floor` is unsupported. Fleets sharing a `Selectivity_index`
estimate one block, so place the prior on the group's lead fleet (its
first that is not `Off`): a prior on a follower would penalize the
shared block once per sharing fleet. On an `Off` fleet a prior is
dropped with a warning and the fleet's linkage coefficients are held at
their starting values, since the fleet has no data; on an `Off` member
of a shared block it moves to the block's lead instead.

Like every selectivity linkage, a covariate acts in the hindcast years;
projection years hold the last hindcast year's curve.

## Examples

``` r
# \donttest{
# A cold-pool effect on the ascending inflection of a logistic fleet
build_selectivity(linkages = list(
  inf_asc = linkage_spec(~ cold_pool, by = ~ fleet)))
#> $linkages
#> $linkages$inf_asc
#> <Rceattle linkage spec>
#>   param:   inf_asc
#>   formula: ~cold_pool
#>   link:    log
#> 
#> 

# A normal prior on the ascending inflection (intercept-only formula)
build_selectivity(linkages = list(
  inf_asc = linkage_spec(~ 1, priors = list(`(Intercept)` = normal(0, 3)))))
#> $linkages
#> $linkages$inf_asc
#> <Rceattle linkage spec>
#>   param:   inf_asc
#>   formula: ~1
#>   prior:    (Intercept) ~ normal(0, 3)
#>   link:    log
#> 
#> 
# }
```
