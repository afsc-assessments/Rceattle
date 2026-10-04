# Look up what a CEATTLE reported quantity is

`fit$quantities` uses the model's own abbreviated names. This returns
the table mapping each one to what it means, the units it is in, how it
is shaped, whether it holds a standard error, and what the same quantity
is called in the NOAA standardized assessment output.

## Usage

``` r
quantity_dictionary(quantity = NULL, process = NULL)
```

## Arguments

- quantity:

  Report names as they appear in `names(fit$quantities)`, default `NULL`
  for every documented quantity.

- process:

  Character vector restricting the table to one or more parts of the
  model: `"population"`, `"recruitment"`, `"mortality"`, `"fishing"`,
  `"reference_points"`, `"catchability"`, `"selectivity"`, `"growth"`,
  `"composition"`, `"predation"`, `"likelihood"`, `"linkage"` or
  `"internal"`.

## Value

A data frame with one row per quantity and columns `quantity`,
`process`, `meaning`, `units`, `dims`, `se` and `standard_label`.

## Details

Units follow the model throughout: numbers-at-age are **thousands of
fish** and weight-at-age is **kg**, so their product is **mt**. A
quantity whose units read `"mt or thousands of fish"` follows its
fleet's `Observation_units` column.

`se = TRUE` means the TMB template `ADREPORT`s the quantity, so
`fit$sdrep` holds a standard error for it and
[`as.data.frame.Rceattle()`](https://afsc-assessments.github.io/Rceattle/reference/as.data.frame.Rceattle.md)
can fill `se`, `lwr` and `upr`. `se = FALSE` means no standard error
exists anywhere on the fit for that quantity. Nothing has a standard
error when the fit was produced with `fit_control(getsd = FALSE)`, which
leaves `sdrep` NULL.

`standard_label` gives the name the quantity takes in the NOAA
standardized assessment output that `stockplotr` and `asar` consume, so
a table can be relabelled for those tools; `NA` marks a quantity that
standard has no equivalent for. Rceattle's own names are kept as the
primary vocabulary because several quantities (predation, multispecies)
have no standard name.

Every per-recruit reference point (`SPR0`, `SPRlimit`, `SPRtarget`,
`SPRFinit`, `NbyageSPR`) is computed only under `msmMode = 0` and is
exactly **zero on a multispecies fit**: M there holds predation
mortality, which scales with predator abundance, so spawning output per
recruit is not a property of the prey stock alone.

## See also

[`parameter_dictionary()`](https://afsc-assessments.github.io/Rceattle/reference/parameter_dictionary.md)
for the estimated parameters,
[`as.data.frame.Rceattle()`](https://afsc-assessments.github.io/Rceattle/reference/as.data.frame.Rceattle.md)
for these quantities in tidy form, and
[`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)
for the fitted object they come from.

## Examples

``` r
# The whole table
head(quantity_dictionary())
#>              quantity    process
#> 1             biomass population
#> 2         log_biomass population
#> 3                 ssb population
#> 4             log_ssb population
#> 5 exploitable_biomass population
#> 6   biomass_depletion population
#>                                                                                                                              meaning
#> 1                                                                                   Total stock biomass, summed over sexes and ages.
#> 2              Total stock biomass on the log scale; its standard error is the CV of biomass, which is the form an ABC buffer wants.
#> 3                                                                               Female spawning-stock biomass at the spawning month.
#> 4                                               Female spawning-stock biomass on the log scale; its standard error is the CV of SSB.
#> 5 Biomass selected by the fisheries, summed over fishery fleets; exactly zero for a survey-only species or when proj_F_prop is zero.
#> 6                                                                          Total biomass relative to unfished biomass, biomass / B0.
#>        units         dims   se   standard_label
#> 1         mt [nspp, nyrs] TRUE          biomass
#> 2     log mt [nspp, nyrs] TRUE             <NA>
#> 3         mt [nspp, nyrs] TRUE spawning_biomass
#> 4     log mt [nspp, nyrs] TRUE             <NA>
#> 5         mt [nspp, nyrs] TRUE             <NA>
#> 6 proportion [nspp, nyrs] TRUE             <NA>

# What is ssb_depletion, and what units is it in?
quantity_dictionary("ssb_depletion")
#>        quantity    process
#> 1 ssb_depletion population
#>                                                                                                                                                                                                                                                                                           meaning
#> 1 Female spawning biomass relative to unfished, ssb / SB0; the quantity a Tier 3 harvest control rule compares against B40%. NA for a species with input numbers-at-age in single-species mode unless DynamicHCR = TRUE, where it is those numbers relative to themselves (1 at spawn_month = 0).
#>        units         dims   se            standard_label
#> 1 proportion [nspp, nyrs] TRUE relative_spawning_biomass

# Everything that holds a standard error
dict <- quantity_dictionary()
dict[dict$se, c("quantity", "meaning")]
#>                quantity
#> 1               biomass
#> 2           log_biomass
#> 3                   ssb
#> 4               log_ssb
#> 5   exploitable_biomass
#> 6     biomass_depletion
#> 7         ssb_depletion
#> 11                    R
#> 12                log_R
#> 17                 R_sd
#> 53        log_index_hat
#> 58       log_sel_at_age
#> 89         beta_linkage
#> 93     beta_linkage_obs
#> 103          pop_scalar
#>                                                                                                                                                                                                                                                                                                            meaning
#> 1                                                                                                                                                                                                                                                                 Total stock biomass, summed over sexes and ages.
#> 2                                                                                                                                                                                            Total stock biomass on the log scale; its standard error is the CV of biomass, which is the form an ABC buffer wants.
#> 3                                                                                                                                                                                                                                                             Female spawning-stock biomass at the spawning month.
#> 4                                                                                                                                                                                                                             Female spawning-stock biomass on the log scale; its standard error is the CV of SSB.
#> 5                                                                                                                                                                               Biomass selected by the fisheries, summed over fishery fleets; exactly zero for a survey-only species or when proj_F_prop is zero.
#> 6                                                                                                                                                                                                                                                        Total biomass relative to unfished biomass, biomass / B0.
#> 7                  Female spawning biomass relative to unfished, ssb / SB0; the quantity a Tier 3 harvest control rule compares against B40%. NA for a species with input numbers-at-age in single-species mode unless DynamicHCR = TRUE, where it is those numbers relative to themselves (1 at spawn_month = 0).
#> 11                                                                                                                                                   Recruitment: numbers entering at the youngest age bin. For a species with input numbers-at-age (estDynamics > 0), the input recruits, with no standard error.
#> 12                                                                                                                                                                                                                                      Recruitment on the log scale; its standard error is the CV of recruitment.
#> 17                                                                                                                                                                                                                                Standard deviation of the recruitment deviations, sigma_R, on the natural scale.
#> 53                                                                                                                                                                                                                                                                        Predicted survey index on the log scale.
#> 58  Log selectivity at age, reported only under fit_control(selectivity_se = TRUE), for a delta-method interval on the curve. One entry per estimated age-based lead fleet, sex, age at or above the fleet's first selected bin, and hindcast year; log_sel_at_age_index says which. Use exp(value +/- 1.96 * sd).
#> 89                                                                                                                                                                                                                    Fitted coefficients of the environmental linkage formulas, one per row of the linkage table.
#> 93                                                                                                                                                                                                QAR1 effect size scaling the latent deviate into the linked parameter; length 0 without a state-space covariate.
#> 103                                                                                                                                                                                                        Multiplier on user-supplied numbers-at-age; estimated for estDynamics = 2 under predation, 1 otherwise.

# Every reference point
quantity_dictionary(process = "reference_points")
#>      quantity          process
#> 1      Flimit reference_points
#> 2     Ftarget reference_points
#> 3          B0 reference_points
#> 4         SB0 reference_points
#> 5         SBF reference_points
#> 6   DynamicB0 reference_points
#> 7  DynamicSB0 reference_points
#> 8  DynamicSBF reference_points
#> 9        SPR0 reference_points
#> 10   SPRlimit reference_points
#> 11  SPRtarget reference_points
#> 12   SPRFinit reference_points
#> 13  NbyageSPR reference_points
#> 14    NByage0 reference_points
#> 15    NByageF reference_points
#>                                                                                                                                  meaning
#> 1           Limit fishing mortality, the FOFL proxy (F35% under Tier 3) used in projections; NA for a species with input numbers-at-age.
#> 2  Target fishing mortality, the maximum FABC proxy (F40% under Tier 3) used in projections; NA for a species with input numbers-at-age.
#> 3                                                             Total biomass at F = 0, carrying the estimated stock-recruit relationship.
#> 4                                                         Female spawning biomass at F = 0; the B100% the Tier 3 proxies are taken from.
#> 5                                                    Female spawning biomass at F = Ftarget; NA for a species with input numbers-at-age.
#> 6               Total biomass under the realized recruitment history with F set to zero. A species with input numbers-at-age keeps them.
#> 7     Female spawning biomass under the realized recruitment history with F set to zero. A species with input numbers-at-age keeps them.
#> 8      Female spawning biomass under the realized recruitment history with F set to Ftarget; NA for a species with input numbers-at-age.
#> 9                                                                         Spawning biomass per recruit at F = 0. Zero under msmMode > 0.
#> 10                       Spawning biomass per recruit at F = Flimit. Zero under msmMode > 0; NA for a species with input numbers-at-age.
#> 11                      Spawning biomass per recruit at F = Ftarget. Zero under msmMode > 0; NA for a species with input numbers-at-age.
#> 12                                                                    Spawning biomass per recruit at F = Finit. Zero under msmMode > 0.
#> 13    Survivors per recruit behind the four SPR schedules; the first dimension is F = 0, Flimit, Ftarget, Finit. Zero under msmMode > 0.
#> 14                                                    Numbers at age at mean recruitment and F = 0, the age structure behind B0 and SB0.
#> 15                                                     Numbers at age at mean recruitment and F = Ftarget, the age structure behind SBF.
#>                units                      dims    se        standard_label
#> 1              yr^-1                    [nspp] FALSE                  <NA>
#> 2              yr^-1                    [nspp] FALSE                  <NA>
#> 3                 mt              [nspp, nyrs] FALSE                  <NA>
#> 4                 mt              [nspp, nyrs] FALSE spawning_biomass_zero
#> 5                 mt              [nspp, nyrs] FALSE                  <NA>
#> 6                 mt              [nspp, nyrs] FALSE                  <NA>
#> 7                 mt              [nspp, nyrs] FALSE                  <NA>
#> 8                 mt              [nspp, nyrs] FALSE                  <NA>
#> 9     kg per recruit                    [nspp] FALSE                  <NA>
#> 10    kg per recruit                    [nspp] FALSE                  <NA>
#> 11    kg per recruit                    [nspp] FALSE                  <NA>
#> 12    kg per recruit                    [nspp] FALSE                  <NA>
#> 13  fish per recruit          [4, nspp, nages] FALSE                  <NA>
#> 14 thousands of fish [nspp, nsex, nages, nyrs] FALSE                  <NA>
#> 15 thousands of fish [nspp, nsex, nages, nyrs] FALSE                  <NA>

# Search the meanings by keyword
dict[grep("depletion", dict$meaning, ignore.case = TRUE), ]
#> [1] quantity       process        meaning        units          dims          
#> [6] se             standard_label
#> <0 rows> (or 0-length row.names)
```
