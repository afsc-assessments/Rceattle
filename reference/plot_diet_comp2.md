# Plot diet composition fits (aggregation-aware)

Diagnostic plots for diet-composition fits that adapt to how each
predator-prey interaction is aggregated (see
[`plot_diet_comp()`](https://afsc-assessments.github.io/Rceattle/reference/plot_diet_comp.md)
for the aggregation conventions):

- prey-age aggregated (predator age resolved): line plot of observed vs
  estimated diet proportion against predator age, with 95% CI ribbons;

- predator aggregated (prey age resolved): line plot against prey age;

- both aggregated: dodged bar plot of observed vs estimated proportion;

- fully disaggregated: predator-age x prey-age bubble grids (observed /
  estimated / Pearson residual).

The observed proportions, fitted proportions and Pearson residuals come
from the diet data path (`quantities$diet_hat` /
`residuals(source = "diet")` / `data_list$diet_data`) via the internal
`.diet_plot_data()` helper. Observed 95% CI ribbons are a normal
approximation to the binomial proportion from `Sample_size`. Estimated
95% CI ribbons are drawn only when the `sdreport` exposes a `diet_hat`
standard error; the C++ template `REPORT()`s but does not `ADREPORT()`
`diet_hat`, so the estimated ribbon is unavailable unless the template
is changed to `ADREPORT(diet_hat)` (the code path is retained so it
activates automatically if that ever happens, or when `sdrep` is
`NULL`).

## Usage

``` r
plot_diet_comp2(Rceattle, file = NULL, species = NULL)
```

## Arguments

- Rceattle:

  A single Rceattle model object.

- file:

  Optional file-name prefix for saved figures.

- species:

  Optional species names for the legend.

## Value

Invisibly returns a list of the printed plot objects.

## See also

[`plot_diet_comp()`](https://afsc-assessments.github.io/Rceattle/reference/plot_diet_comp.md),
[`plot_diet_comp1()`](https://afsc-assessments.github.io/Rceattle/reference/plot_diet_comp1.md)

## Examples

``` r
if (FALSE) { # \dontrun{
plot_diet_comp2(my_msm_fit)
} # }
```
