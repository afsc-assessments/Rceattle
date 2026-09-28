# Helper to set map for growth parameters

Maps the growth parameters (`log_growth_pars`, `growth_log_sd`,
`weight_length_pars`) according to `growth_model`. Everything starts
mapped off and each growth function turns on the parameters it uses.

Time-varying growth comes from the linkage grammar
(`build_growth(linkages = )`), whose random effects hold their own
density and map.

## Usage

``` r
build_map_growth(map_list, data_list, nyrs_hind)
```

## Arguments

- map_list:

  The current TMB map list.

- data_list:

  The data list containing model settings.

- nyrs_hind:

  Number of historical years.

## Value

Updated `map_list`.
