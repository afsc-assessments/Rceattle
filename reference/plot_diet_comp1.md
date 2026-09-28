# Plot diet composition fits (bubble/grid diagnostic)

`plot_diet_comp1()` is an alias of
[`plot_diet_comp()`](https://afsc-assessments.github.io/Rceattle/reference/plot_diet_comp.md):
it produces the predator-age x prey-age bubble grids (observed /
estimated / Pearson residual) for each predator-prey-year interaction.
It is provided as a named entry point alongside
[`plot_diet_comp2()`](https://afsc-assessments.github.io/Rceattle/reference/plot_diet_comp2.md)
(the aggregation-aware line / bar / bubble plots) for scripts that call
both by name.

## Usage

``` r
plot_diet_comp1(Rceattle, file = NULL, species = NULL)
```

## Arguments

- Rceattle:

  Single Rceattle model object exported from `Rceattle`.

- file:

  Optional file-name prefix for saved figures.

- species:

  Optional species names for the legend.

## Value

Invisibly returns a list of the printed plot grids (see
[`plot_diet_comp()`](https://afsc-assessments.github.io/Rceattle/reference/plot_diet_comp.md)).

## See also

[`plot_diet_comp()`](https://afsc-assessments.github.io/Rceattle/reference/plot_diet_comp.md),
[`plot_diet_comp2()`](https://afsc-assessments.github.io/Rceattle/reference/plot_diet_comp2.md)

## Examples

``` r
if (FALSE) { # \dontrun{
plot_diet_comp1(my_msm_fit)
} # }
```
