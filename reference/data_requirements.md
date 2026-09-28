# Report which data inputs a model configuration requires

Introspects the Rceattle data-element catalogue and reports, for a given
model configuration, which top-level `data_list` elements are
**Required**, **Optional** (used if supplied, otherwise default-filled
by
[`clean_data()`](https://afsc-assessments.github.io/Rceattle/reference/clean_data.md)),
or **Ignored** (not consulted because the feature that would use them is
switched off). It answers "what do I actually need to supply for *this*
model?" without having to read the validation code or the switch tables,
the same conditions enforced at fit time (they share one declarative
table).

## Usage

``` r
data_requirements(
  data_list = NULL,
  msmMode = 0,
  growth_model = 0,
  estDynamics = 0,
  Selectivity = NULL,
  Index_distribution = NULL,
  Ceq = NULL
)
```

## Arguments

- data_list:

  Optional. An existing Rceattle data list (e.g. from
  [`read_data()`](https://afsc-assessments.github.io/Rceattle/reference/read_data.md)
  /
  [`build_data()`](https://afsc-assessments.github.io/Rceattle/reference/build_data.md),
  or a bundled dataset). When supplied, the configuration is read from
  it and the convenience arguments are ignored.

- msmMode:

  Predation mode used to build the configuration when `data_list` is
  `NULL`. `0` single-species (default), `> 0` multispecies.

- growth_model:

  Growth mode (scalar or per-species). `> 0` estimates growth and
  requires `caal_data`.

- estDynamics:

  Numbers-at-age mode (scalar or per-species). `> 0` fixes
  numbers-at-age and requires `NByageFixed`.

- Selectivity:

  Optional character vector of per-fleet selectivity forms (used only
  when `data_list` is `NULL`) so the `emp_sel` requirement
  (`Selectivity = "Fixed"`) can be evaluated. Named for the
  `fleet_control` column it stands in for.

- Index_distribution:

  Optional character vector of per-fleet index likelihood families (used
  only when `data_list` is `NULL`) so the `index_cov` requirement
  (`Index_distribution = "MVN"`) can be evaluated. Named for the
  `fleet_control` column it stands in for.

- Ceq:

  Optional per-species consumption-equation codes (used only when
  `data_list` is `NULL`); `> 1` requires an environmental temperature
  index.

## Value

A `data.frame` with one row per data element and columns:

- `element`:

  the `data_list` element name.

- `category`:

  grouping (dimensions / biology / fishery / composition / predation /
  environment).

- `status`:

  `"Required"`, `"Optional"`, or `"Ignored"`.

- `condition`:

  the condition under which the element is required (`"always"` for the
  core backbone).

- `default`:

  for Optional elements, the
  [`clean_data()`](https://afsc-assessments.github.io/Rceattle/reference/clean_data.md)
  default used when the element is absent.

Rows are ordered Required, then Optional, then Ignored.

## Details

The configuration can be given either as an existing (possibly partial)
`data_list`, its switches are normalized through
[`clean_data()`](https://afsc-assessments.github.io/Rceattle/reference/clean_data.md)
/
[`switch_check()`](https://afsc-assessments.github.io/Rceattle/reference/switch_check.md)
so the conditions evaluate against filled defaults, or, when no
`data_list` is supplied, built from the convenience arguments.

Requirements are *conditional*: e.g. `diet_data`, `ration_data` and the
bioenergetics scalars are Ignored under single-species (`msmMode = 0`)
but Required under multispecies (`msmMode > 0`); `caal_data` is Optional
unless `growth_model > 0`; `NByageFixed` is Ignored unless
`estDynamics > 0`; `emp_sel` is Required only for fleets with
`Selectivity = "Fixed"`; `index_cov` only for
`Index_distribution = "MVN"`.

## See also

[`clean_data()`](https://afsc-assessments.github.io/Rceattle/reference/clean_data.md),
[`build_data()`](https://afsc-assessments.github.io/Rceattle/reference/build_data.md).

## Examples

``` r
# Single-species: predation/diet inputs are Ignored, comp_data Optional.
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

# Multispecies: diet_data and bioenergetics become Required.
data_requirements(msmMode = 1)
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
#> 23            bioenergetics   predation Required
#> 24                diet_data   predation Required
#> 25              ration_data   predation Required
#> 26                caal_data composition Optional
#> 27                comp_data composition Optional
#> 28                 env_data environment Optional
#> 29               index_data     fishery Optional
#> 30              NByageFixed  dimensions  Ignored
#> 31                  emp_sel     fishery  Ignored
#> 32                index_cov     fishery  Ignored
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
#> 23                      msmMode > 0 (multispecies)
#> 24                      msmMode > 0 (multispecies)
#> 25                      msmMode > 0 (multispecies)
#> 26         any growth_model > 0 (estimated growth)
#> 27                     optional (used if supplied)
#> 28 any Ceq > 1 (temperature-dependent consumption)
#> 29                     optional (used if supplied)
#> 30      any estDynamics > 0 (fixed numbers-at-age)
#> 31                any fleet Selectivity == 'Fixed'
#> 32           any fleet Index_distribution == 'MVN'
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
#> 23                                
#> 24                                
#> 25                                
#> 26    empty caal_data (clean_data)
#> 27    empty comp_data (clean_data)
#> 28 Year-only env_data (clean_data)
#> 29          no survey index fitted
#> 30                                
#> 31                                
#> 32                                

# From a bundled dataset.
data_requirements(BS2017SS)
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
#> 21                 env_data environment Required
#> 22               catch_data     fishery Required
#> 23                  emp_sel     fishery Required
#> 24            fleet_control     fishery Required
#> 25            bioenergetics   predation Required
#> 26                diet_data   predation Required
#> 27              ration_data   predation Required
#> 28                caal_data composition Optional
#> 29                comp_data composition Optional
#> 30               index_data     fishery Optional
#> 31              NByageFixed  dimensions  Ignored
#> 32                index_cov     fishery  Ignored
#>                                          condition                      default
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
#> 21 any Ceq > 1 (temperature-dependent consumption)                             
#> 22                                          always                             
#> 23                any fleet Selectivity == 'Fixed'                             
#> 24                                          always                             
#> 25                      msmMode > 0 (multispecies)                             
#> 26                      msmMode > 0 (multispecies)                             
#> 27                      msmMode > 0 (multispecies)                             
#> 28         any growth_model > 0 (estimated growth) empty caal_data (clean_data)
#> 29                     optional (used if supplied) empty comp_data (clean_data)
#> 30                     optional (used if supplied)       no survey index fitted
#> 31      any estDynamics > 0 (fixed numbers-at-age)                             
#> 32           any fleet Index_distribution == 'MVN'                             
```
