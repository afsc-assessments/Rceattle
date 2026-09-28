# Main function to construct the TMB map argument for CEATTLE

Builds the TMB map, which tells TMB which parameters to estimate and
which to hold fixed (a fixed parameter is mapped to `NA`), and which
parameters share a single estimated value. One helper handles each
process block (recruitment, M1, predation, selectivity, catchability,
...).

## Usage

``` r
build_map(
  data_list,
  params,
  debug = FALSE,
  random_rec = FALSE,
  random_sel = FALSE,
  random_q = NULL
)
```

## Arguments

- data_list:

  an Rceattle data_list

- params:

  A parameter list created from
  [`build_params`](https://afsc-assessments.github.io/Rceattle/reference/build_params.md).

- debug:

  Logical. If TRUE, fixes every parameter except the dummy (maps all to
  `NA`), so the model runs with no parameters estimated.

- random_rec:

  Logical. If TRUE, treats recruitment deviations as random effects,
  meaning the variance parameter (`R_log_sd`) is estimated.

- random_sel:

  Logical. If TRUE, treats selectivity deviations as random effects,
  meaning the variance parameter (`sel_dev_log_sd`) is estimated.

- random_q:

  Logical. If TRUE, treats catchability deviations as random effects,
  meaning the variance parameter (`index_q_dev_log_sd`) is estimated.
  Defaults to `data_list$random_q`, which is what
  [`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)
  stores there.

## Value

A list containing the factorized TMB map (`mapFactor`) and the original
map matrix/array list (`mapList`).

## Details

Fleets sharing a `Selectivity_index` or a `Catchability_index` have the
lead fleet's map slice copied over the rest of the group, so a per-fleet
setting that differs within a group is resolved here rather than
reported. Those disagreements are checked in
[`data_check()`](https://afsc-assessments.github.io/Rceattle/reference/data_check.md),
which
[`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)
runs first; `build_map()` does not call it, so a caller invoking
`build_map()` directly should run
[`data_check()`](https://afsc-assessments.github.io/Rceattle/reference/data_check.md)
too.
