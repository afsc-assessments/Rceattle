# Print method for an Rceattle data list

Shows the model specification as an indented tree, dimensions, fleets
and their selectivity / catchability forms, configured processes, active
linkages, and any attached
[`model_config()`](https://afsc-assessments.github.io/Rceattle/reference/model_config.md),
rather than dumping the full data list.

## Usage

``` r
# S3 method for class 'Rceattle_data'
print(x, config = TRUE, ...)

# S3 method for class 'Rceattle_data'
summary(object, config = FALSE, ...)
```

## Arguments

- x:

  An `"Rceattle_data"` object from
  [`build_data()`](https://afsc-assessments.github.io/Rceattle/reference/build_data.md).

- config:

  Show the attached
  [`model_config()`](https://afsc-assessments.github.io/Rceattle/reference/model_config.md)
  block. Defaults to `TRUE` for
  [`print()`](https://rdrr.io/r/base/print.html) and `FALSE` for
  [`summary()`](https://rdrr.io/r/base/summary.html). Has no effect on
  an object that holds no configuration.

- ...:

  Currently unused.

- object:

  An `"Rceattle_data"` object.

## Value

`x`, invisibly.

## Details

[`print()`](https://rdrr.io/r/base/print.html) includes the
configuration block; [`summary()`](https://rdrr.io/r/base/summary.html)
omits it and reports the data structure alone. Either can be overridden
with `config`.
