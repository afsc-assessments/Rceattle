# Data inputs for CEATTLE model with WHAM-estimated growth

A data list containing inputs for a CEATTLE model that uses growth
estimated from the Woods Hole Assessment Model (WHAM). See
[`BS2017SS`](https://afsc-assessments.github.io/Rceattle/reference/BS2017SS.md)
for format details. Its two 2DAR1 fleets ship with
`Time_varying_sel_sd = 0` and no `N_sel_bins`, so
[`data_check()`](https://afsc-assessments.github.io/Rceattle/reference/data_check.md)
refuses it as shipped; the vignettes set a different selectivity before
fitting.

## Usage

``` r
whamGrowthData
```

## Format

A list with the same structure as
[`BS2017SS`](https://afsc-assessments.github.io/Rceattle/reference/BS2017SS.md).
