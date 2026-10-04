# Split joint-sex (Sex == 3) composition bins onto a single age/length axis

Joint-sex compositions stack females in bins `1..nbin` and males in bins
`nbin+1..2*nbin` (where `nbin` is `nages` or `nlengths` for the
species). This re-bases the male bins to `1..nbin` and tags the source
label by sex so males and females face the same bin axis, matching
[`plot_comp()`](https://afsc-assessments.github.io/Rceattle/reference/plot_comp.md).
Rows with Sex != 3 (single-sex or combined) are returned unchanged.

## Usage

``` r
.osa_jointsex(df, nages, nlengths)
```

## Arguments

- df:

  A data frame with `species`, `sex`, `index_label`, `age_length_bin`,
  and `source` columns.

- nages, nlengths:

  Per-species bin counts (or `NULL` to skip the split).

## Details

The label says "joint" because the two series are halves of ONE
decomposition unit, not two: a joint row is appended as a single density
over the stacked vector with one cell dropped for the sum-to-one
constraint, not one per sex. An assessment author read the previous " -
male" / " - female" labels as two separate distributions and asked why
the residuals did not match the likelihood; they do, and the label now
says so.
