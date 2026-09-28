# Helper to set map for Selectivity parameters

Maps base selectivity parameters (`log_sel_slp`, `sel_inf`, `sel_coff`)
and time-varying deviations, based on `Selectivity` and
`Time_varying_sel` settings in `fleet_control`.

`Selectivity` in `fleet_control` of the data determines shape of
selectivity curve: `0` = "Fixed" empirical selectivity provided in
`emp_sel` in the data `1` = "Logistic" `2` = "NonParametric" selectivity
sensu Ianelli et al 2018 `3` = "DoubleLogistic" `4` =
"DescendingLogistic" `5` = "Hake" non-parametric selectivity sensu
Taylor et al 2014 (Hake) `6` = "2DAR1" across age x year `7` = "3DAR1"
across age x cohort x year (Cheng et al 2024) `8` = "DoubleNormal"
Gaussian ascending/descending limbs blended at peak (analogous to SS3
pattern 24). Parameters: `sel_inf[1]` = peak; `sel_inf[2]` =
logit(right_floor) (right-tail floor, SS3 P6/end_logit);
`log_sel_slp[1]` = log(sigma_asc); `log_sel_slp[2]` = log(sigma_desc).
right_floor-\>0: dome-shaped; right_floor-\>1: logistic ascending only.

`N_sel_bins` Number of age/length bins to estimate for non-parametric
and AR1 selectivity (Selectivity = 2, 5, 6, 7, 9, or 13). Not used
otherwise

`Time_varying_sel` determines if time-varying selectivity should be
estimated for logistic, double logistic selectivity, descending logistic
, non-parametric, or hake (`Selectivity = 1, 2, 3, 4, or 5`). `0` =
'Off' `1` = 'IID' penalized deviates given `sel_sd_prior` `3` = 'Block'
time blocks with no penalty `4` = 'RandomWalk' random walk `5` =
'RandomWalkAscending' random walk on ascending portion of double
logistic only. `random_sel` in `fit_mod` treats random deviates and
random walk parameters as random effects, estimating the variance.

## Usage

``` r
build_map_selectivity(map_list, data_list, nyrs_hind, random_sel)
```

## Arguments

- map_list:

  The current TMB map list.

- data_list:

  The data list containing model settings.

- nyrs_hind:

  Number of historical years.

- random_sel:

  Logical indicating if selectivity deviations are random effects.

## Value

Updated `map_list`.
