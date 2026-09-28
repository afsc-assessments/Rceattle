# Helper to set map for linkage-table parameters

Maps `beta_linkage` (one entry per row of `data_list$linkage_table`).
Rows whose `est_phase == 0` are fixed at their initial values via `NA`;
`(Intercept)` rows are also fixed (their value stays at 0, the base
parameter holds the level). Everything else is estimated. Phased
estimation honoring nonzero phase ordinals can layer on later via the
`phase` argument to
[`fit_control()`](https://afsc-assessments.github.io/Rceattle/reference/fit_control.md).

## Usage

``` r
build_map_linkages(map_list, data_list)
```

## Arguments

- map_list:

  The current TMB map list.

- data_list:

  an Rceattle data_list (with the pooled `linkage_table` from
  `pool_linkages()`).

## Value

Updated `map_list`.
