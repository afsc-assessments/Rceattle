# Helper to turn off base linked parameters

Maps the base parameter (`rec_pars`, `log_M1`, `log_growth_pars`) out of
estimation only for stratum groups whose linkage formula holds *no*
intercept. With an intercept in the formula (`~ 1`, `~ temp`, ...) and a
nonzero `est_phase` the base parameter holds the level and stays
estimable; the linkage `(Intercept)` row is fixed at 0 instead. When an
intercept row has `est_phase == 0` (the documented "fix at init"
contract), the base parameter is mapped off as well, so a fixed
intercept truly fixes the parameter at its
[`build_params()`](https://afsc-assessments.github.io/Rceattle/reference/build_params.md)
initial value. Slope-only formulas (`~ 0 + temp`) emit no intercept row,
so we mask the base parameter to keep it at its
[`build_params()`](https://afsc-assessments.github.io/Rceattle/reference/build_params.md)
default and let the slope rows define the year-by-year offset.

## Usage

``` r
map_linkage_adjuster(map_list, data_list)
```

## Arguments

- map_list:

  The current TMB map list.

- data_list:

  an Rceattle data_list (with the pooled `linkage_table` from
  `pool_linkages()`).

## Value

Updated `map_list`.
