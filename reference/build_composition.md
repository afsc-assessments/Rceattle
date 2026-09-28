# Composition-weighting specification

Holds **priors** on the Dirichlet-multinomial composition-weighting
overdispersion. The DM weight (the "theta" that scales the effective
sample size) is otherwise an unpenalized free parameter; a linkage lets
you put a prior on it through the same grammar as every other parameter.
The three parameters target the three DM likelihoods:

- `theta_comp`:

  age / length composition weight (`comp_weights`), per fleet.

- `theta_caal`:

  conditional-age-at-length weight (`caal_weights`), per fleet.

- `theta_diet`:

  diet composition weight (`diet_comp_weights`), per predator.

This is a **prior-only** process: the DM weight is a scalar, not a
year-varying quantity, so only intercept formulas (`~ 1`) with a
`priors` entry are meaningful. The prior is placed on the natural-scale
DM weight `theta = exp(weight)`, so a
[`gamma()`](https://rdrr.io/r/base/Special.html) prior reads naturally.
A linkage on a fleet whose `Comp_distribution` (or `CAAL_distribution`)
is not `"DirichletMultinomial"` errors, since the weight has no effect
there.

## Usage

``` r
build_composition(linkages = NULL)
```

## Arguments

- linkages:

  Optional named list of
  [`linkage_spec()`](https://afsc-assessments.github.io/Rceattle/reference/linkage_spec.md)
  objects keyed by `theta_comp` / `theta_caal` (per fleet by default,
  `by = ~ fleet`) or `theta_diet` (per predator by default,
  `by = ~ species`).

## Value

A list of composition-weighting settings for
[`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md).

## Examples

``` r
# \donttest{
# A weak gamma prior on the DM overdispersion of every fleet's age comps
build_composition(linkages = list(
  theta_comp = linkage_spec(~ 1, by = ~ fleet,
                            priors = list(`(Intercept)` = gamma(2, 0.5)))))
#> $linkages
#> $linkages$theta_comp
#> <Rceattle linkage spec>
#>   param:   theta_comp
#>   formula: ~1
#>   prior:    (Intercept) ~ gamma(2, 0.5)
#>   link:    log
#> 
#> 
# }
```
