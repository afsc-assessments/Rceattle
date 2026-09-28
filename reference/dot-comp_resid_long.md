# Tidy long-format composition observed / fitted / Pearson table

Reuses
[`residuals.Rceattle()`](https://afsc-assessments.github.io/Rceattle/reference/residuals.Rceattle.md)
(`type = "pearson"`, `source = "comp"`) for the observed and fitted
proportions and the Pearson residual, the single source of truth, then
adds the plotting-only columns: joint-sex bins re-based onto a single
age/length axis (males stored as bins `nbin+1 .. 2*nbin` are mapped to
`1 .. nbin` and tagged `sex_grp = "male"`) and the facet labels. Zero
observed proportions are kept; only `NA` observed/fitted are dropped (by
[`residuals()`](https://rdrr.io/r/stats/residuals.html)).

## Usage

``` r
.comp_resid_long(Rceattle, species = NULL)
```

## Arguments

- Rceattle:

  A fitted `Rceattle` model.

- species:

  Optional species code(s) to keep.

## Value

A data frame, or `NULL` when there are no composition residuals.
