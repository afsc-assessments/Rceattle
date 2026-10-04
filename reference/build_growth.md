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
  choice), otherwise `"WHAM"`. When SS3's `Growth_Age_for_L2` is 999,
  SS3 pins the plus group to the upper anchor, which is `"WHAM"` here.

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

- sd_form:

  What the two growth-variability endpoints (`sd_L1`, `sd_Linf`) are:
  `"SD"`, standard deviations of length-at-age in cm (SS3
  `CV_Growth_Pattern` 2), or `"CV"`, coefficients of variation so the SD
  is CV x mean length (SS3 pattern 0). Scalar or length-`nspp`; default
  `NA` inherits `data_list$growth_sd_form`, otherwise `"SD"`.

- plus_group_length:

  How the plus group's mean length is set: `"M1"`, `"none"`, `"SS3.24"`
  or `"decay"` (see Details). Scalar or length-`nspp`; default `NA`
  inherits `data_list$growth_plus_length`, otherwise `"M1"`.

- plus_group_decay:

  Annual decay rate (per year) for `plus_group_length = "decay"`, SS3's
  positive `Linf_decay`, roughly the plus group's total mortality.
  Required for, and only used by, `"decay"`.

- pop_lengths:

  Population length bins (lower edges, cm) on which the age-length key,
  weight-at-length and maturity-at-length are computed before being
  summed into the data length bins: SS3's population length bins. A
  vector applies to every species, a list gives one per species. Every
  data-bin lower edge must also be a population-bin edge. Default `NULL`
  inherits `data_list$pop_lengths`, otherwise uses the data bins.

## Value

A list of switches defining the growth model.

## Details

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
