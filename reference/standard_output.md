# Rceattle output in the NOAA standardized assessment format

Relabels
[`report_tables()`](https://afsc-assessments.github.io/Rceattle/reference/report_tables.md)
output into the standardized assessment output that the `stockplotr` and
`asar` packages consume, so Rceattle results can be plotted and written
into a report by the same tools used for SS3, BAM, WHAM and FIMS.

## Usage

``` r
standard_output(x, species = NULL, model = NULL)
```

## Arguments

- x:

  An `"rceattle_report"` from
  [`report_tables()`](https://afsc-assessments.github.io/Rceattle/reference/report_tables.md),
  or an Rceattle fit, which is passed through
  [`report_tables()`](https://afsc-assessments.github.io/Rceattle/reference/report_tables.md)
  first.

- species:

  Which stock to emit: a species name, or an index into the species this
  report holds, required when it holds more than one.

- model:

  Which model to emit when `x` holds several, defaulting to the first.

## Value

A data frame with the standardized columns, in the standard's own order,
plus `species` and `model`.

## Details

Quantity names are translated through the `standard_label` column of
[`quantity_dictionary()`](https://afsc-assessments.github.io/Rceattle/reference/quantity_dictionary.md),
so `ssb` becomes `spawning_biomass`, `R` becomes `recruitment`, `F_spp`
becomes `fishing_mortality`, and so on. A quantity the standard has no
name for keeps its Rceattle name, so nothing is silently dropped.

**The standard has no species dimension.** It describes one stock, so a
multispecies CEATTLE fit cannot be represented in it as a whole. A
`species` column is held alongside the standard columns and `species`
selects one stock; with several species in the fit and no selection,
this errors rather than returning a frame in which two stocks' biomass
share a year.

`era` follows the standard's own vocabulary: `"time"` for the hindcast
the model was fit to, `"fore"` for the projection.

## See also

[`report_tables()`](https://afsc-assessments.github.io/Rceattle/reference/report_tables.md)
for the native tables and
[`quantity_dictionary()`](https://afsc-assessments.github.io/Rceattle/reference/quantity_dictionary.md)
for the name crosswalk.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- fit_mod(data_list = BS2017SS, estimateMode = 1)
std <- standard_output(fit, species = 1)
# ...then, with the stockplotr package installed:
# stockplotr::plot_spawning_biomass(std)
} # }
```
