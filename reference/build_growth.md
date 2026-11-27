# Specify the growth model for Rceattle

Specify the growth model for Rceattle

## Usage

``` r
build_growth(
  fun = "empirical",
  growth_age_L1 = NA,
  sd_plus_group = NA,
  linkages = NULL,
  sd_form = NA,
  plus_group_length = NA,
  plus_group_decay = NA,
  pop_lengths = NULL
)
```

## Arguments

- fun:

  Growth function: a string
  ([GROWTH_FUNS](https://afsc-assessments.github.io/Rceattle/reference/GROWTH_FUNS.md):
  `"empirical"` (default), `"vonBertalanffy"`, `"Richards"`) or the
  equivalent integer code (`0`, `1`, `2`). The canonical string is
  stored on the returned object.

- growth_age_L1:

  Von Bertalanffy / Richards anchor age, where mean length equals `L1`;
  SS3's `Growth_Age_for_L1`. Falls back to `max(0.5, minage[sp])`, which
  keeps `minage >= 1` models backwards-compatible and gives `minage = 0`
  models an SS3-consistent half-year anchor.

- sd_plus_group:

  How the oldest age class's SD-at-age is treated, and only for
  estimated growth. `"WHAM"` (the fallback) pins the plus-group SD to
  the upper anchor `exp(sd_Linf)`, the WHAM SDAA convention; `"SS3"`
  interpolates it by length like any interior age. Integer codes `1`/`2`
  are accepted. Where SS3's `Growth_Age_for_L2` is 999 it pins the plus
  group to the upper anchor, which is `"WHAM"` here.

- linkages:

  Optional named list of
  [`linkage_spec()`](https://afsc-assessments.github.io/Rceattle/reference/linkage_spec.md)
  objects keyed by parameter name (one of
  [GROWTH_LINKAGE_PARAMS](https://afsc-assessments.github.io/Rceattle/reference/GROWTH_LINKAGE_PARAMS.md)).
  The mean-growth keys (`K`, `L1`, `Linf`, `m`) take arbitrary one-sided
  formulas and make that parameter year-varying, as a per-year offset
  around its mean. The SD-endpoint keys (`sd_L1`, `sd_Linf`) honour only
  intercept-bearing formulas, typically `~ 1`, threading `init`,
  `bounds` and `priors` onto the growth SD-at-age so the SDs get the
  same prior, fix and starting-value contract as the means. A slope row
  on an SD spec warns and has no effect; a slope-only formula
  (`~ 0 + temp`) is an error.

- sd_form:

  What the two variability endpoints (`sd_L1`, `sd_Linf`) are: `"SD"`
  (the fallback), standard deviations of length-at-age in cm, SS3's
  `CV_Growth_Pattern` 2; or `"CV"`, coefficients of variation, so the SD
  is CV x mean length, SS3 pattern 0.

- plus_group_length:

  How the plus group's mean length is set: `"M1"` (the fallback),
  `"none"`, `"SS3.24"` or `"decay"`. See Details.

- plus_group_decay:

  Annual decay rate for `plus_group_length = "decay"`, SS3's positive
  `Linf_decay`, roughly the plus group's total mortality. Required for,
  and used only by, `"decay"`.

- pop_lengths:

  Population length bins (lower edges, cm) on which the age-length key,
  weight-at-length and maturity-at-length are computed before being
  summed into the data length bins: SS3's population length bins. Every
  data-bin lower edge must also be a population-bin edge. A vector
  applies to every species, a list gives one per species, and the
  fallback is the data bins.

## Value

A list of switches defining the growth model.

## Details

Every switch above takes a scalar, recycled across species, or a
length-`nspp` vector for per-species values. Each defaults to `NA`,
which inherits the matching field from `data_list` when it is set – so a
refit, or a model built by the SS3 converter, keeps the original choice
– and otherwise uses the fallback named with it.

**Plus-group mean length.** Fish older than the oldest age are pooled,
so the plus group's mean length lies between the growth curve at the
oldest age, L_A, and L-infinity:

- `"M1"`: mean of L_A + (a/n)(Linf - L_A) over a = 0..n (n = `nages`),
  weighted by survival at the oldest age's base M1.

- `"none"`: L_A (SS3 `Linf_decay = -998`).

- `"SS3.24"`: as `"M1"` with weights exp(-0.2 a) and a = 0..A, A the
  oldest age (SS3 `Linf_decay = -999`).

- `"decay"`: L_A and 2A further ages, each one more year along the von
  Bertalanffy curve, weighted by exp(-d a) with d = `plus_group_decay`
  (SS3 `Linf_decay = d`).

Under `"none"`, `"SS3.24"` and `"decay"` the plus group also grows
within the year like every other age, as SS3 does; under `"M1"` it keeps
its Jan-1 length through the year.

**Maturity-at-length** is set in the data, not here: the per-species
control columns `L50_mat_len` and `slope_mat_len`.

## References

Richards, F.J. 1959. A flexible growth function for empirical use. J.
Exp. Bot. 10:290-301.

## Examples

``` r
# \donttest{
# Sex-specific von Bertalanffy with temperature on K, by species + sex
build_growth(
  fun = "vonBertalanffy",   # or fun = 1
  linkages = list(
    K = linkage_spec(
      formula = ~ temp,
      by      = ~ species + sex,
      priors  = list(temp = normal(0, 1))
    )
  )
)
#> $fun
#> [1] "vonBertalanffy"
#> 
#> $linkages
#> $linkages$K
#> <Rceattle linkage spec>
#>   param:   K
#>   formula: ~temp
#>   prior:    temp ~ normal(0, 1)
#>   link:    log
#> 
#> 
#> $growth_model
#> [1] 1
#> 
#> $sd_plus_group
#> [1] NA
#> 
#> $growth_sd_style
#> [1] NA
#> 
#> $growth_age_L1
#> [1] NA
#> 
#> $sd_form
#> [1] NA
#> 
#> $growth_sd_form
#> [1] NA
#> 
#> $plus_group_length
#> [1] NA
#> 
#> $growth_plus_length
#> [1] NA
#> 
#> $plus_group_decay
#> [1] NA
#> 
#> $pop_lengths
#> NULL
#> 
# }
```
