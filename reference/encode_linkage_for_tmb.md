# Encode a linkage table into TMB-friendly parallel vectors

Encode a linkage table into TMB-friendly parallel vectors

## Usage

``` r
encode_linkage_for_tmb(table, X)
```

## Arguments

- table:

  an `Rceattle_linkage_table` (typically the output of
  `pool_linkages()`).

- X:

  the global design matrix from `pool_linkages()` (passed through
  unchanged so callers can stash it alongside).

## Value

A named list with components:

- n_linkage:

  `integer(1)`. Number of coefficients (`nrow(table)`).

- linkage_process, linkage_param, linkage_species, linkage_sex,
  linkage_age_bin, linkage_fleet, linkage_X_col, linkage_link,
  linkage_re_index, linkage_prior_family:

  Parallel `integer(n_linkage)` vectors. `linkage_re_index` is `-1` on
  fixed rows, else the 0-based `beta_linkage_re` slot.

- linkage_prior_p1, linkage_prior_p2:

  Parallel `numeric(n_linkage)` vectors.

- linkage_X:

  The design matrix as passed in.

- linkage_design_names:

  `character` colnames of `linkage_X`.
