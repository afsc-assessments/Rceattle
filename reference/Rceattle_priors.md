# Prior distributions for Rceattle linkage coefficients

Each prior is captured as a small object of class `"Rceattle_prior"`
holding a `family` name and two positional parameters (`p1`, `p2`). The
two-parameter shape is enforced so the linkage table can store priors
uniformly as four columns (`prior_family`, `prior_p1`, `prior_p2`, plus
reserved future use).

## Details

Two surfaces are provided:

1.  **Programmatic constructors**: exported with the `prior_` prefix
    ([`prior_normal()`](https://afsc-assessments.github.io/Rceattle/reference/prior_normal.md),
    [`prior_lognormal()`](https://afsc-assessments.github.io/Rceattle/reference/prior_lognormal.md),
    ...), safe to call anywhere without masking
    [`base::gamma()`](https://rdrr.io/r/base/Special.html)/[`base::beta()`](https://rdrr.io/r/base/Special.html).

2.  **Inline form** inside `linkage_spec(priors = ...)`. There the
    argument is captured unevaluated and resolved with a private data
    mask that makes `normal()`, `lognormal()`,
    [`gamma()`](https://rdrr.io/r/base/Special.html), and
    [`beta()`](https://rdrr.io/r/base/Special.html) shorthand for the
    respective `prior_*` constructors, *only* inside that argument. Base
    R remains untouched at the package namespace.
