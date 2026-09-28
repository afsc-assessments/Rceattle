# Rearrange a data_list for TMB

Function to rearrange a `data_list` object to be read into TMB

`rearrange_dat()` is a deprecated alias for `rearrange_data()`, kept so
older scripts keep running. It is scheduled for removal in 6.0.0; use
`rearrange_data()`.

## Usage

``` r
rearrange_data(data_list, build_osa = FALSE)

rearrange_dat(data_list)
```

## Arguments

- data_list:

  an Rceattle data_list

- build_osa:

  Logical. Passed to
  [`build_osa_data()`](https://afsc-assessments.github.io/Rceattle/reference/build_osa_data.md);
  when `TRUE` the full one-step-ahead (OSA) observation data is
  assembled so
  [`osa_residuals()`](https://afsc-assessments.github.io/Rceattle/reference/osa_residuals.md)
  can be computed. Default `FALSE` (the fast path used by simulation
  testing). The composition proportion offset is read from
  `data_list$comp_offset` (default `1e-5`).
