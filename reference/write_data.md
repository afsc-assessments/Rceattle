# Write a data list to a CEATTLE xlsx workbook

Write a data list to a CEATTLE xlsx workbook

## Usage

``` r
write_data(data_list, file = "Rceattle_data.xlsx")
```

## Arguments

- data_list:

  Rceattle data_list object

- file:

  Filename to write. Must end in '.xlsx'

## See also

[`build_data()`](https://afsc-assessments.github.io/Rceattle/reference/build_data.md)
to assemble or edit a data list in R,
[`read_data()`](https://afsc-assessments.github.io/Rceattle/reference/read_data.md)
to read one back. A
[`model_config()`](https://afsc-assessments.github.io/Rceattle/reference/model_config.md)
slot is not written to the workbook and does not survive the round-trip
(a warning is issued); persist it separately with
[`save_config()`](https://afsc-assessments.github.io/Rceattle/reference/save_config.md)
/
[`load_config()`](https://afsc-assessments.github.io/Rceattle/reference/load_config.md).

## Examples

``` r

library(Rceattle)
data(BS2017SS)
out_file <- file.path(tempdir(), "BS2017SS.xlsx")
write_data(data_list = BS2017SS, file = out_file)
file.remove(out_file)
#> [1] TRUE
```
