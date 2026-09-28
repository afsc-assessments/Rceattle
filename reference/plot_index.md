# Survey index fits

Plots fitted index series: observed points with 95% observation
intervals and the model-predicted index, faceted by fleet. Every fleet
holding `index_data` is drawn, a fishery with a CPUE series as much as a
survey.

## Usage

``` r
plot_index(
  Rceattle,
  file = NULL,
  model_names = NULL,
  species = NULL,
  incl_proj = FALSE,
  width = 7,
  height = 6.5,
  error = TRUE,
  log = FALSE,
  line_col = NULL,
  right_adj = 0,
  top_adj = 0.05,
  single.plots = FALSE
)
```

## Arguments

- Rceattle:

  A single
  [`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)
  object or a list of them (overlaid).

- file:

  Optional file stem; if given the figure is written to
  `<file>_survey_indices.png`.

- model_names:

  Legend labels for the models.

- species:

  Species (indices) to include. Default all.

- incl_proj:

  Include projection years.

- width, height:

  Saved figure size (inches).

- error:

  Draw observed points and error bars.

- log:

  Plot the index on the log scale.

- line_col, right_adj, top_adj, single.plots:

  Deprecated base-graphics arguments, retained for back-compatibility
  and ignored.

## Value

A `ggplot` object.

## Details

The interval follows the fleet's `Index_distribution`: lognormal
quantiles on a `"Lognormal"` fleet, and `observed +/- 1.96 sd` truncated
at zero on a natural-scale one (`"MVN"`, `"MVNORM"`, `"Normal"`,
`"TruncatedNormal"`), whose sd is absolute and in the units of the
index.
