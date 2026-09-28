# Allowed recruitment-parameter names for `linkages` in [`build_srr()`](https://afsc-assessments.github.io/Rceattle/reference/build_srr.md)

Natural-scale names of the underlying recruitment parameters that the
linkage system can address. Linkages on `R0` act under mean recruitment
(the offset is added to log mean recruitment with the default log link);
under a hindcast curve a single-species `R0` is derived from alpha and
beta, so an `R0` linkage is refused there, and a multispecies one takes
an intercept only. Linkages on `alpha` and `beta` only do work when the
model has a curve (Beverton-Holt, Ricker).

## Usage

``` r
RECRUITMENT_LINKAGE_PARAMS
```
