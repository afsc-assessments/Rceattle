# Look up what a CEATTLE parameter is

The TMB parameter vector uses transformed, abbreviated names, `log_M1`,
`R_log_sd`, `index_log_q`. This returns the table mapping each one to
the quantity it represents on its natural scale, which process it
belongs to, what it means, and its dimensions.

## Usage

``` r
parameter_dictionary(internal = NULL, process = NULL)
```

## Arguments

- internal:

  Character vector of internal parameter names, as they appear in
  `names(fitted$estimated_params)` or in an `sdreport()`. Default `NULL`
  returns every parameter.

- process:

  Character vector restricting the table to one or more parts of the
  model: `"recruitment"`, `"mortality"`, `"growth"`, `"fishing"`,
  `"catchability"`, `"selectivity"`, `"observation"`, `"predation"`,
  `"linkage"` or `"internal"`.

## Value

A data frame with one row per parameter and columns `internal`,
`natural`, `process`, `meaning` and `dims`.

## Details

A name beginning `log_`, or ending `_log_sd`, is estimated on the log
scale, so [`exp()`](https://rdrr.io/r/base/Log.html) turns it into the
quantity named in the `natural` column. Most fitted values are already
back-transformed in `fitted$quantities`, so this is mainly needed when
reading `estimated_params` or an `sdreport()` directly.

Dimensions use the model's own notation: `nspp` species, `nsex` sexes,
`nages` age bins, `nyrs` years including the projection, `nyrs_hind`
hindcast years only, `n_flt` fleets and `n_sel` selectivity blocks.
`n_sel` equals `n_flt`: fleets sharing a `Selectivity_index` share one
block through the map, not through a narrower array.

## See also

[`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)
for the fitted object these names come from.

## Examples

``` r
# The whole table
head(parameter_dictionary())
#>         internal                natural     process
#> 1          dummy                  dummy    internal
#> 2 log_pop_scalar             pop_scalar    internal
#> 3       rec_pars      R0 / alpha / beta recruitment
#> 4        rec_dev recruitment deviations recruitment
#> 5       R_log_sd                sigma_R recruitment
#> 6       init_dev initial-age deviations recruitment
#>                                                                                                                                              meaning
#> 1                                                                             Placeholder parameter; the only free parameter under estimateMode = 4.
#> 2                                            Multiplier on user-supplied numbers-at-age; estimated for estDynamics = 2 under predation, 1 otherwise.
#> 3                                                                    Stock-recruit parameters: mean recruitment (or R0), and the SRR alpha and beta.
#> 4                                                                                                           Annual log-scale recruitment deviations.
#> 5                                                           Standard deviation of the recruitment deviations; estimated only when random_rec = TRUE.
#> 6 Deviations defining the initial (first-year) age structure; one per age from minage + 1 upward, since age minage in the first year is recruitment.
#>              dims
#> 1             [1]
#> 2          [nspp]
#> 3       [nspp, 3]
#> 4    [nspp, nyrs]
#> 5          [nspp]
#> 6 [nspp, nages-1]

# What is R_log_sd?
parameter_dictionary("R_log_sd")
#>   internal natural     process
#> 1 R_log_sd sigma_R recruitment
#>                                                                                    meaning
#> 1 Standard deviation of the recruitment deviations; estimated only when random_rec = TRUE.
#>     dims
#> 1 [nspp]

# Everything governing selectivity
parameter_dictionary(process = "selectivity")
#>          internal                                        natural     process
#> 1     log_sel_slp                              selectivity slope selectivity
#> 2         sel_inf                         selectivity inflection selectivity
#> 3  log_sel_apical                      selectivity apical height selectivity
#> 4 log_sel_slp_dev                               slope deviations selectivity
#> 5     sel_inf_dev                          inflection deviations selectivity
#> 6        sel_coff                       selectivity coefficients selectivity
#> 7    sel_coff_dev                         coefficient deviations selectivity
#> 8  sel_dev_log_sd                                      sigma_sel selectivity
#> 9   sel_curve_pen selectivity penalty weights / AR1 correlations selectivity
#>                                                                                                                                                                                                                                                                                                                                        meaning
#> 1                                                                                                                                                                                                                                                                        Logistic-family selectivity slope; row 1 ascending, row 2 descending.
#> 2                                                                                                                                                                                                                                                              Logistic-family age/length at 50% selection; row 1 ascending, row 2 descending.
#> 3                                                                                                                                                                           Log multiplier on one sex's whole curve, applied after the form and before normalization; 0 = no offset. Estimated only through a selectivity linkage on `apical`.
#> 4                                                                                                                                                                                                                                                                                                  Annual deviations on the selectivity slope.
#> 5                                                                                                                                                                                                                                                                                       Annual deviations on the selectivity inflection point.
#> 6                                                                                                                                                                                                                                                                                              Non-parametric selectivity-at-bin coefficients.
#> 7                                                                                                                                                                                                                                                                            Annual deviations on the non-parametric selectivity coefficients.
#> 8                                                                                                                                                                                                                                                                               Standard deviation of the time-varying selectivity deviations.
#> 9 Meaning depends on the fleet's Selectivity. For the non-parametric forms these are weights on the shape and curvature penalties, supplied via fleet_control and not estimated. For '2DAR1' and '3DAR1' the same slots hold estimated correlations: slot 1 across selectivity BINS, slot 2 across YEARS, and for 3DAR1 slot 3 across COHORTS.
#>                                   dims
#> 1                     [2, n_sel, nsex]
#> 2                     [2, n_sel, nsex]
#> 3                        [n_sel, nsex]
#> 4          [2, n_sel, nsex, nyrs_hind]
#> 5          [2, n_sel, nsex, nyrs_hind]
#> 6            [n_sel, nsex, n_sel_bins]
#> 7 [n_sel, nsex, n_sel_bins, nyrs_hind]
#> 8                              [n_sel]
#> 9                           [n_sel, 3]

# Search the meanings by keyword
dict <- parameter_dictionary()
dict[grep("catchability", dict$meaning, ignore.case = TRUE), ]
#>              internal       natural      process
#> 27        index_log_q             q catchability
#> 29        index_q_dev  q deviations catchability
#> 30     index_q_log_sd sigma_q_prior catchability
#> 31 index_q_dev_log_sd   sigma_q_dev catchability
#>                                                                        meaning
#> 27                                                  Survey/index catchability.
#> 29                   Annual deviations on catchability when q is time-varying.
#> 30 Standard deviation of the prior on catchability (used when Estimate_q = 2).
#> 31             Standard deviation of the time-varying catchability deviations.
#>                  dims
#> 27            [n_flt]
#> 29 [n_flt, nyrs_hind]
#> 30            [n_flt]
#> 31            [n_flt]
```
