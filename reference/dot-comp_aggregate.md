# Pool composition counts across the fitted hindcast years

Multiplies each year's proportions by its input sample size, sums the
counts over years (for one fleet x type panel), and puts the total back
on the proportion scale; joint-sex groups keep their shared
normalization (females + males sum to 1). Pooling counts rather than
averaging proportions is what stops a year with 20 otoliths counting as
much as one with 2000.

## Usage

``` r
.comp_aggregate(d, endyr = NULL)
```

## Arguments

- d:

  One panel's rows from
  [`.comp_resid_long()`](https://afsc-assessments.github.io/Rceattle/reference/dot-comp_resid_long.md).

- endyr:

  Last year the composition likelihood fits, or `NULL` to keep every
  hindcast row. Default `NULL`.

## Value

A data frame with `bin`, `sex_grp`, the pooled `obs` / `hat` proportions
and their counts `o_n` / `e_n`, the pooled count variance `v_n`, the 95%
prediction limits `lwr` / `upr`, the summed input and assumed effective
sample sizes `ISS` / `ESS`, and the sex-mirrored drawing columns
`y_obs`, `y_hat`, `y_lwr`, `y_upr`.

## Details

The pooled count is a sum of independent draws, so its variance is the
exact sum of the per-year variances the fleet's own likelihood assumes
(the `Sd` column). The resulting interval treats the fitted proportions
as known, so it is slightly narrower than one that held the estimation
uncertainty in `p_hat` as well.
