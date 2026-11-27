# Allowed growth functions for `fun` in [`build_growth()`](https://afsc-assessments.github.io/Rceattle/reference/build_growth.md)

Currently:

- `"empirical"`: use the empirical weight-at-age input directly, no
  parameters estimated.

- `"vonBertalanffy"`: sex-specific von Bertalanffy (parameters `K`,
  `L1`, `Linf`).

- `"Richards"`: sex-specific Richards (1959) curve (adds `m`).

## Usage

``` r
GROWTH_FUNS
```

## References

Richards, F.J. 1959. A flexible growth function for empirical use. J.
Exp. Bot. 10:290-301.
