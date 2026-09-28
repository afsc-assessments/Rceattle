# Locate every estimated parameter in the model's own coordinates

One row per estimated (fixed-effect) parameter, giving the fleet,
species, sex, age, bin, year or slot it occupies.

## Usage

``` r
parameter_index(object)
```

## Arguments

- object:

  A fitted Rceattle model from
  [`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md).

## Value

A data frame with one row per estimated parameter and columns
`par_index` (position in `obj$par`), `block`, `species`, `fleet`, `sex`,
`age`, `bin`, `year`, `slot`, `n_cells` and `label`.

## Details

TMB names each element of `obj$par` after its parameter block, so a
convergence diagnostic can report that `sel_coff_dev` is
non-identifiable but not which of its elements.

The mapping is recovered by pushing a tagged parameter vector back
through TMB's own `parList()`, so it reflects
[`build_map()`](https://afsc-assessments.github.io/Rceattle/reference/build_map.md)
exactly: parameters mapped off are absent, and fleets sharing a
`Selectivity_index` or `Catchability_index` appear once, as the single
parameter they are, with `n_cells` recording how many array cells they
drive.

Ages are absolute. `init_dev` holds one deviation per age from
`minage + 1` upward, since age `minage` in the first year is
recruitment.

Coordinate columns are `NA` where the axis does not apply to a block.
The linkage and environmental-covariate blocks (`beta_linkage`,
`M1_beta`, ...) are indexed by a linkage-table row rather than a model
coordinate, so they hold no coordinates and are reported by element
number.

A selectivity `slot` is named from the fleet's `Selectivity`: slot 2 of
`sel_inf` is a descending inflection for the double-logistic family, the
right-tail floor under `DoubleNormal`, and the age-1 selectivity under
`LogisticPM`.

Under `estimateMode = "Estimate"` with any HCR but `"NoFishing"`,
`object$obj` is the projection object, whose only free parameters are
`log_Ftarget` / `log_Flimit`, so the index describes those, not the
hindcast. `fit$convergence` reads the hindcast index that
[`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)
captured before the projection remapped the model, and checks it against
the vector it is labelling either way.

## See also

[`parameter_dictionary()`](https://afsc-assessments.github.io/Rceattle/reference/parameter_dictionary.md)
for what each block means.
