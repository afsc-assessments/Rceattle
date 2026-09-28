# Helper to set map for Catchability parameters

Maps catchability base parameters (`index_log_q`), time-varying
deviations (`index_q_dev`), and environmental linkages (`index_q_beta`)
for every fleet that holds fitted `index_data`, a fishery with a CPUE
series as much as a survey. A fleet with no index rows gets none of
them, whatever its `Catchability` says, since a q with no index to
inform it is a flat direction. Sharing overrides this:
[`adjust_map_shared_params()`](https://afsc-assessments.github.io/Rceattle/reference/adjust_map_shared_params.md)
then copies the lead fleet's slice across each `Catchability_index`
group.

## Usage

``` r
build_map_catchability(map_list, data_list, nyrs_hind, random_q = FALSE)
```

## Arguments

- map_list:

  The current TMB map list.

- data_list:

  The data list containing model settings.

- nyrs_hind:

  Number of historical years.

- random_q:

  Logical. If TRUE, the deviation sd (`index_q_dev_log_sd`) is estimated
  rather than held at `Time_varying_q_sd`.

## Value

Updated `map_list`.
