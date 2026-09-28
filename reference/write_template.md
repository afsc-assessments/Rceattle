# Write a minimal starter CEATTLE data workbook

Emits a small, structurally complete single-species workbook, one survey
and one fishery, flat placeholder data, that a user can open and edit as
a starting point. fleet_control is built on the canonical column names
with schema defaults filled by
[`switch_check`](https://afsc-assessments.github.io/Rceattle/reference/switch_check.md),
so the template is always in sync with the current schema. The template
round-trips through
[`read_data`](https://afsc-assessments.github.io/Rceattle/reference/read_data.md)
and
[`data_check`](https://afsc-assessments.github.io/Rceattle/reference/data_check.md)
and builds under `fit_mod(estimateMode = 3)`; replace the placeholder
observations with real data before fitting.

## Usage

``` r
write_template(
  file = "Rceattle_data_template.xlsx",
  nages = 10,
  nyrs = 30,
  nprojyrs = 12,
  minage = 1
)
```

## Arguments

- file:

  Output path for the `.xlsx` workbook.

- nages:

  Number of ages to template (default 10).

- nyrs:

  Number of hindcast years (default 30).

- nprojyrs:

  Number of projection years beyond the hindcast (default 12).

- minage:

  Minimum age (recruitment age; default 1).

## Value

Invisibly, the minimal `data_list` that was written.

## Examples

``` r
f <- file.path(tempdir(), "template.xlsx")
dat <- write_template(f, nages = 8, nyrs = 20)

# The workbook is correctly shaped but holds placeholder observations;
# data_requirements() reports what a given configuration still needs.
head(data_requirements(dat), 4)
#>            element category   status condition default
#> 1          M1_base  biology Required    always        
#> 2        age_error  biology Required    always        
#> 3 age_trans_matrix  biology Required    always        
#> 4         maturity  biology Required    always        
unlink(f)
```
