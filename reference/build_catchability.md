# Catchability specification

Holds environmental linkages on survey/index catchability `q`. The
effect of an `env_data` covariate is written as a formula and can hold
priors, bounds, and an estimation phase like any other linkage.

## Usage

``` r
build_catchability(linkages = NULL)
```

## Arguments

- linkages:

  Optional named list of
  [`linkage_spec()`](https://afsc-assessments.github.io/Rceattle/reference/linkage_spec.md)
  objects keyed by catchability parameter. The only parameter is `q`.
  Coefficients are per fleet by default (`by = ~ fleet`); use the
  `fleet` argument of
  [`linkage_spec()`](https://afsc-assessments.github.io/Rceattle/reference/linkage_spec.md)
  to restrict a spec to particular fleets.

## Value

A list of catchability settings for
[`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md).

## Details

Catchability is the one process that accepts `link = "exponential"`,
Stock Synthesis's environmental link type 1: `q^exp(beta * x)`, which
MULTIPLIES `log q` where `"log"` shifts it. It needs an estimated base
`q` and a lognormal or t index likelihood, and cannot share a fleet with
a random-effect q linkage; each of those is refused with the reason. See
[`vignette("environmental-linkages-and-priors")`](https://afsc-assessments.github.io/Rceattle/articles/environmental-linkages-and-priors.md).

## Examples

``` r
# \donttest{
# One temperature effect on q per fleet. `by = ~ fleet` is the default, and a
# fleet without an estimated survey catchability is not linkable, so at fit time
# restrict the linkage to the estimated-q fleets (here fleet 7 of BS2017SS).
build_catchability(linkages = list(
  q = linkage_spec(~ temp, fleet = 7)))
#> $linkages
#> $linkages$q
#> <Rceattle linkage spec>
#>   param:   q
#>   formula: ~temp
#>   fleet:   7
#>   link:    log
#> 
#> 

# Restrict it to fleets 1 and 3, with a prior on the slope
build_catchability(linkages = list(
  q = linkage_spec(~ temp, fleet = c(1, 3),
                   priors = list(temp = prior_normal(0, 1)))))
#> $linkages
#> $linkages$q
#> <Rceattle linkage spec>
#>   param:   q
#>   formula: ~temp
#>   prior:    temp ~ normal(0, 1)
#>   fleet:   1, 3
#>   link:    log
#> 
#> 
# }
```
