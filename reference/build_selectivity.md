# Selectivity specification

Holds environmental linkages on selectivity parameters. The effect on a
parameter is written as a formula and composes additively with any
`Time_varying_sel` process error on the same fleet (the two are separate
mechanisms: a covariate effect versus a deviation).

The parameter names are the shape parameters of the parametric
selectivity forms:

- `slp_asc`, `slp_desc`:

  ascending / descending logistic slope (log scale); for a double-normal
  the ascending / descending width, aliased `sigma_asc` / `sigma_desc`.

- `inf_asc`, `inf_desc`:

  ascending / descending inflection age/length (natural scale); for a
  double-normal the peak and the logit right-floor, aliased `peak` /
  `right_floor`.

- `coff`:

  non-parametric selectivity-at-bin coefficients.

- `apical`:

  a multiplier on one sex's whole curve (log scale), applied after the
  form and before normalization, so every estimated form takes it. Name
  the fleet and the sex that holds it (`by = ~ fleet + sex`,
  `fleet = 3`, `sex = "male"`), as Stock Synthesis's male-offset option
  does; the other sex is the reference. See Details.

Every parameter but `apical` accepts `link = "log"` (multiplicative on
the natural parameter) or `link = "identity"` (additive), like the other
processes; `apical` is a multiplier already and takes `"log"` only.

## Usage

``` r
build_selectivity(linkages = NULL)
```

## Arguments

- linkages:

  Optional named list of
  [`linkage_spec()`](https://afsc-assessments.github.io/Rceattle/reference/linkage_spec.md)
  objects keyed by selectivity parameter. Coefficients are per fleet by
  default (`by = ~ fleet`); use the `fleet` argument of
  [`linkage_spec()`](https://afsc-assessments.github.io/Rceattle/reference/linkage_spec.md)
  to restrict a spec to particular fleets.

## Value

A list of selectivity settings for
[`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md).

## Details

**The `apical` offset.** Fishing mortality is one `log_F` per fleet and
year shared by the sexes, so a sex difference in F can only come from
selectivity, and no form has a height parameter: the logistic family and
DoubleNormal peak at 1 for every sex, the non-parametric forms re-centre
each sex, Hake normalizes each sex by its own maximum. `apical`
multiplies one sex's curve by `exp(log_sel_apical)`, bin by bin. That
equals the ratio of the sexes' peak heights only where their shapes peak
equally (the logistic family on an age axis); for a dome with
sex-specific shape, read it as the multiplier on that sex's curve and
take the peak ratio from `fit$quantities$sel_at_age`. Only the contrast
between the sexes is identified (the common level is `log_F`), so one
sex holds it and the fit is refused if both do, if no fleet or no sex is
named, if the species has one sex, on a `Fixed`, AR1 or
`Fleet_type = "Off"` fleet, on a fleet that shares another's
`Selectivity_index` block, or under `link = "identity"`, which could
drive the multiplier negative; use the default `link = "log"`. It is
also refused where `Sel_norm_scope = "WithinSex"` normalization would
divide it straight back out; use `"AcrossSexes"`, under which the
more-selected sex peaks at 1, or turn `Sel_norm_bin` off. The contrast
is informed only by joint composition (`comp_data$Sex = 3`); with
single-sex compositions it rests on its prior, and
[`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)
warns when a fleet has neither. Read the fitted multiplier with
`exp(fit$estimated_params$log_sel_apical[fleet, sex])`, and the realized
ratio of the sexes' maxima from `fit$quantities$sel_at_age`. Naming
`fleet` and `sex` is enough: `by` defaults to `~ fleet + sex` for this
parameter. An intercept prior is on the multiplier's natural scale
(`lognormal()` centred on 1 means no offset). Like every selectivity
linkage, a covariate on it acts in the hindcast years; projection years
hold the last hindcast year's curve.

**Priors on a selectivity parameter.** An intercept-only formula (`~ 1`)
with a `priors` entry places a prior on the selectivity parameter itself
(no year-to-year offset is added). Read the prior on the parameter's own
scale: the slopes (`slp_asc` / `slp_desc`) are on the log scale (use
`lognormal()`), the inflections (`inf_asc` / `inf_desc`) on the natural
scale (use `normal()`). See Examples for a normal prior on the ascending
inflection. This mirrors the prior-only
[`build_composition()`](https://afsc-assessments.github.io/Rceattle/reference/build_composition.md)
path.

A selectivity prior targets one parameter, so in a two-sex model an
unstratified `~ 1` prior constrains sex 1 only, use `by = ~ sex` for a
per-sex prior. An `init` on a selectivity intercept has no effect (the
starting value comes from the data), and a prior on the double-normal
`right_floor` is not supported. Fleets sharing a `Selectivity_index`
estimate one parameter block, so place the prior on the group's lead
fleet (its first fleet that is not `Off`); a prior on a follower would
penalize the shared block once per sharing fleet. A prior on an `Off`
fleet is refused: its selectivity is not estimated.

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
