# Build an Rceattle data list in R

Assemble (or edit) an Rceattle data list in R code rather than from an
xlsx workbook. Supply only the data blocks a model uses, dimensions,
biology (`weight`, `maturity`, `sex_ratio`, `M1_base`), a
`fleet_control`, and the observation tables (`catch_data`, `index_data`,
`comp_data`, ...), and the optional blocks a single-species model does
not need (`caal_data`, `emp_sel`, `diet_data`, ...) are default-filled
by
[`clean_data()`](https://afsc-assessments.github.io/Rceattle/reference/clean_data.md).
The result is the same bare list
[`read_data()`](https://afsc-assessments.github.io/Rceattle/reference/read_data.md)
returns and round-trips through
[`write_data()`](https://afsc-assessments.github.io/Rceattle/reference/write_data.md)
unchanged.

## Usage

``` r
build_data(base = NULL, file = NULL, ..., .check = TRUE)
```

## Arguments

- base:

  Optional existing Rceattle data list to start from and override.

- file:

  Optional path to an Rceattle xlsx workbook; read via
  [`read_data()`](https://afsc-assessments.github.io/Rceattle/reference/read_data.md)
  to form the starting object. Supply at most one of `base` / `file`.

- ...:

  Named top-level data-list elements to set or override (e.g. `nspp`,
  `styr`, `endyr`, `fleet_control`, `weight`, `maturity`, `catch_data`,
  `index_data`, `comp_data`). Every argument must be named.

- .check:

  Logical; run the presence pre-check (default `TRUE`). Set `FALSE` to
  assemble a deliberately partial object.

## Value

An Rceattle data list (a bare `list`, as from
[`read_data()`](https://afsc-assessments.github.io/Rceattle/reference/read_data.md)).

## Details

Three entry points, which may be combined:

- **from blocks**: pass the elements as named arguments:
  `build_data(nspp = 1, styr = 1977, ..., fleet_control = fc, catch_data = catch)`.

- **from a file**: `build_data(file = "model.xlsx", projyr = 2060)`
  reads the workbook via
  [`read_data()`](https://afsc-assessments.github.io/Rceattle/reference/read_data.md),
  then applies the overrides.

- **from an existing object**:
  `build_data(base = BS2017SS, projyr = 2060)` starts from a data list
  (e.g. a bundled dataset or a
  [`combine_data()`](https://afsc-assessments.github.io/Rceattle/reference/combine_data.md)
  result) and overrides the named blocks. This is the common
  copy-and-edit / combine-and-restamp workflow.

Element names are checked against the recognized schema: an override
whose name is a near-miss of a known element (e.g. `maturty`) errors
with a suggestion, so a typo is caught here rather than surfacing much
later in a fit. Legacy top-level names (`fsh_biom`, `srv_biom`, `wt`,
`pmature`, `Pyrs`) are mapped to their canonical equivalents.

Full validation runs at fit time, inside
[`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md).
`build_data()` runs only a light presence pre-check (`.check = TRUE`) so
a missing *required* block is reported at construction with a clear
message. The pre-check reads an attached
[`model_config()`](https://afsc-assessments.github.io/Rceattle/reference/model_config.md),
so a configuration held on the object is accounted for. Requirements
that depend on fit-time settings passed directly to
[`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)
and stored nowhere on the data list are not knowable here and are left
to that fit-time validation; see
[`data_requirements()`](https://afsc-assessments.github.io/Rceattle/reference/data_requirements.md)
to preview them.

## See also

[`read_data()`](https://afsc-assessments.github.io/Rceattle/reference/read_data.md),
[`clean_data()`](https://afsc-assessments.github.io/Rceattle/reference/clean_data.md),
[`data_requirements()`](https://afsc-assessments.github.io/Rceattle/reference/data_requirements.md),
[`combine_data()`](https://afsc-assessments.github.io/Rceattle/reference/combine_data.md),
[`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md).

## Examples

``` r
# Copy-and-edit a bundled dataset.
dat <- build_data(base = BS2017SS, projyr = 2060)

# Preview what a configuration needs before building it.
data_requirements(msmMode = 0)
#>                     element    category   status
#> 1                   M1_base     biology Required
#> 2                 age_error     biology Required
#> 3          age_trans_matrix     biology Required
#> 4                  maturity     biology Required
#> 5  pop_age_transition_index     biology Required
#> 6              pop_wt_index     biology Required
#> 7                 sex_ratio     biology Required
#> 8              ssb_wt_index     biology Required
#> 9                    weight     biology Required
#> 10                    endyr  dimensions Required
#> 11                   minage  dimensions Required
#> 12                    nages  dimensions Required
#> 13                 nlengths  dimensions Required
#> 14                     nsex  dimensions Required
#> 15                     nspp  dimensions Required
#> 16               other_food  dimensions Required
#> 17                   projyr  dimensions Required
#> 18              spawn_month  dimensions Required
#> 19                  spnames  dimensions Required
#> 20                     styr  dimensions Required
#> 21               catch_data     fishery Required
#> 22            fleet_control     fishery Required
#> 23                caal_data composition Optional
#> 24                comp_data composition Optional
#> 25                 env_data environment Optional
#> 26               index_data     fishery Optional
#> 27              NByageFixed  dimensions  Ignored
#> 28                  emp_sel     fishery  Ignored
#> 29                index_cov     fishery  Ignored
#> 30            bioenergetics   predation  Ignored
#> 31                diet_data   predation  Ignored
#> 32              ration_data   predation  Ignored
#>                                          condition
#> 1                                           always
#> 2                                           always
#> 3                                           always
#> 4                                           always
#> 5                                           always
#> 6                                           always
#> 7                                           always
#> 8                                           always
#> 9                                           always
#> 10                                          always
#> 11                                          always
#> 12                                          always
#> 13                                          always
#> 14                                          always
#> 15                                          always
#> 16                                          always
#> 17                                          always
#> 18                                          always
#> 19                                          always
#> 20                                          always
#> 21                                          always
#> 22                                          always
#> 23         any growth_model > 0 (estimated growth)
#> 24                     optional (used if supplied)
#> 25 any Ceq > 1 (temperature-dependent consumption)
#> 26                     optional (used if supplied)
#> 27      any estDynamics > 0 (fixed numbers-at-age)
#> 28                any fleet Selectivity == 'Fixed'
#> 29           any fleet Index_distribution == 'MVN'
#> 30                      msmMode > 0 (multispecies)
#> 31                      msmMode > 0 (multispecies)
#> 32                      msmMode > 0 (multispecies)
#>                            default
#> 1                                 
#> 2                                 
#> 3                                 
#> 4                                 
#> 5                                 
#> 6                                 
#> 7                                 
#> 8                                 
#> 9                                 
#> 10                                
#> 11                                
#> 12                                
#> 13                                
#> 14                                
#> 15                                
#> 16                                
#> 17                                
#> 18                                
#> 19                                
#> 20                                
#> 21                                
#> 22                                
#> 23    empty caal_data (clean_data)
#> 24    empty comp_data (clean_data)
#> 25 Year-only env_data (clean_data)
#> 26          no survey index fitted
#> 27                                
#> 28                                
#> 29                                
#> 30                                
#> 31                                
#> 32                                
```
