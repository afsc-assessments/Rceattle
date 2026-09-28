# Read a CEATTLE excel data file

Read a CEATTLE excel data file

## Usage

``` r
read_data(file = "Rceattle_data.xlsx")
```

## Arguments

- file:

  Filename to read. Must end in '.xlsx'

## See also

[`build_data()`](https://afsc-assessments.github.io/Rceattle/reference/build_data.md)
to assemble or edit a data list in R (and to read a workbook then
override blocks in one call, `build_data(file = ...)`),
[`data_requirements()`](https://afsc-assessments.github.io/Rceattle/reference/data_requirements.md)
to see which inputs a configuration needs.

## Examples

``` r

library(Rceattle)
data(BS2017SS)
out_file <- file.path(tempdir(), "BS2017SS.xlsx")
write_data(data_list = BS2017SS, file = out_file)
data_list <- read_data(file = out_file)
# Or read + override blocks in one call:
data_list <- build_data(file = out_file, projyr = 2060)
file.remove(out_file)
#> [1] TRUE
```
