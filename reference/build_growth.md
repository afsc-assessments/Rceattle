# Specify the growth model for Rceattle

Specify the growth model for Rceattle

## Usage

``` r
build_growth(
  fun = "empirical",
  growth_age_L1 = NA,
  sd_plus_group = NA,
  linkages = NULL
)
```

## Arguments

- fun:

  Growth function. Either a string
  ([GROWTH_FUNS](https://afsc-assessments.github.io/Rceattle/reference/GROWTH_FUNS.md):
  `"empirical"` (default), `"vonBertalanffy"`, `"Richards"`) or the
  equivalent integer code (`0`, `1`, `2`). The canonical string form is
  stored on the returned object.

- growth_age_L1:

  Von Bertalanffy / Richards anchor age (the age at which mean length
  equals `L1`). Matches SS3's `Growth_Age_for_L1` control input. Scalar
  (recycled across species) or a length-`nspp` vector for per-species
  values. Default `NA` inherits `data_list$growth_age_L1` if supplied
  (e.g. from the SS3 converter), otherwise falls back to
  `max(0.5, minage[sp])` so `minage >= 1` models stay
  backwards-compatible and `minage = 0` models pick up an SS3-consistent
  half-year anchor.

- sd_plus_group:

  How the oldest age class's SD-at-age is treated (only affects
  estimated growth, `fun != "empirical"`). `"WHAM"` pins the plus-group
  SD to the upper anchor `exp(sd_Linf)` (the WHAM SDAA convention);
  `"SS3"` instead interpolates it by length like any interior age.
  Accepts a string or the integer code (`1`/`2`), scalar or a
  length-`nspp` vector. Default `NA` inherits
  `data_list$growth_sd_style` if present (so a refit keeps the original
  choice), otherwise `"WHAM"`.

- linkages:

  Optional named list of
  [`linkage_spec()`](https://afsc-assessments.github.io/Rceattle/reference/linkage_spec.md)
  objects keyed by parameter name (must be one of
  [GROWTH_LINKAGE_PARAMS](https://afsc-assessments.github.io/Rceattle/reference/GROWTH_LINKAGE_PARAMS.md)).
  The mean-growth keys (`K`, `L1`, `Linf`, `m`) accept arbitrary
  one-sided formulas and make that growth parameter year-varying (a
  per-year offset around its mean). The SD-endpoint keys (`sd_L1`,
  `sd_Linf`) only honor intercept-bearing formulas (typically `~ 1`),
  they thread `init`, `bounds`, and `priors` onto the growth SD-at-age,
  giving the SDs the same prior/fix/initial-value contract as the mean
  parameters. Slope rows on SD specs raise a warning and have no effect;
  slope-only formulas (`~ 0 + temp`) error.

## Value

A list of switches defining the growth model.

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
# }
```
