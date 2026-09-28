# Allowed growth functions for `fun` in [`build_growth()`](https://afsc-assessments.github.io/Rceattle/reference/build_growth.md)

Currently:

- `"empirical"`: use the empirical weight-at-age input directly, no
  parameters estimated.

- `"vonBertalanffy"`: sex-specific von Bertalanffy (parameters `K`,
  `L1`, `Linf`).

- `"Richards"`: sex-specific Richards (adds `m`).

## Usage

``` r
GROWTH_FUNS
```
