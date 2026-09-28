# Function to construct the TMB map argument for CEATTLE for projecting under alternative harvest control rules

Reads a data list and map to update the map argument based on the HCR
specified in
[`build_hcr`](https://afsc-assessments.github.io/Rceattle/reference/build_hcr.md)

## Usage

``` r
build_hcr_map(
  data_list,
  map,
  debug = FALSE,
  all_params_on = FALSE,
  HCRiter = 1
)
```

## Arguments

- data_list:

  an Rceattle data_list

- map:

  a map object created from
  [`build_map`](https://afsc-assessments.github.io/Rceattle/reference/build_map.md).

- debug:

  logical. If TRUE, turns off all parameters for debugging (default =
  FALSE).

- all_params_on:

  logical. If TRUE, leaves all hindcast parameters turned on (default =
  FALSE).

- HCRiter:

  for multi-species models, the order in which to project fishing (e.g.
  predators first, then prey)

## Value

a list of map arguments for each parameter
