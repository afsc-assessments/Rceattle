# Calculate Predicted Weight-at-Age

Converts a growth matrix (length-at-age probabilities) into mean
weight-at-age using a length-weight relationship (W = a \* L^b).

## Usage

``` r
get_weight_at_age_r(
  nsex_sp,
  nages_sp,
  nlengths_sp,
  nyrs,
  lengths_sp,
  length_at_age,
  growth_matrix,
  lw_params
)
```

## Arguments

- nsex_sp:

  Integer. Number of sexes.

- nages_sp:

  Integer. Number of age classes.

- nlengths_sp:

  Integer. Number of length bins.

- nyrs:

  Integer. Number of years.

- lengths_sp:

  Vector. Lower edges of the length bins (cm), as in
  [`get_growth_matrix_r()`](https://afsc-assessments.github.io/Rceattle/reference/get_growth_matrix_r.md).

- length_at_age:

  Array. Mean length at age from get_growth_matrix_r.

- growth_matrix:

  Array. 4D array (sex, age, length, year) from get_growth_matrix_r.

- lw_params:

  Array. Dimensions (sex, yr, 2). Params: 1st is alpha (a), 2nd is beta
  (b).

## Value

A 3D array of mean weights with dimensions (sex, age, year).

## Details

The function calculates midpoints for length bins to avoid bias. For the
first bin, it assumes the width is equal to the second bin's width. The
final weight-at-age is the expected value across all length bins for
that age.
