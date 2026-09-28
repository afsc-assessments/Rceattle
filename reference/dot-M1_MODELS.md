# String\<-\>integer mapping for `M1_model` in [`build_M1()`](https://afsc-assessments.github.io/Rceattle/reference/build_M1.md)

Either form is accepted by
[`build_M1()`](https://afsc-assessments.github.io/Rceattle/reference/build_M1.md);
the canonical integer code is what the TMB template ultimately consumes.
The env-driven integer codes 4 and 5 still work with a deprecation
warning, their structural part is identical to 1 and 2 respectively, and
the env effect is expressed via the `linkages` argument to
[`build_M1()`](https://afsc-assessments.github.io/Rceattle/reference/build_M1.md)
(see
[`vignette("environmental-linkages-and-priors")`](https://afsc-assessments.github.io/Rceattle/articles/environmental-linkages-and-priors.md)).
No string alias is offered for 4 or 5 to discourage their use in new
code.

## Usage

``` r
.M1_MODELS
```
