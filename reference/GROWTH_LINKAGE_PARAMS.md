# Allowed growth-parameter names for `linkages` in [`build_growth()`](https://afsc-assessments.github.io/Rceattle/reference/build_growth.md)

Natural-scale names of the underlying growth-function parameters. Von
Bertalanffy uses `K`, `L1`, `Linf`; Richards adds `m`. `sd_L1` /
`sd_Linf` are the standard deviations of length-at-age anchored at `L1`
and `Linf` (the SD-at-age interpolation endpoints). Only intercept-only
specs (`~ 1`) are honored on the SD endpoints, they thread through
`init` / `bounds` / `priors` onto the growth SD-at-age but do not vary
by year. The empirical weight-at-age model admits no linkages.

## Usage

``` r
GROWTH_LINKAGE_PARAMS
```
