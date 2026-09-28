# Survey index residuals

Plots survey index residuals by year, faceted by survey fleet
(`residual_type = "pearson"`, the default), or one-step-ahead (OSA)
residual diagnostics for a single fit (`residual_type = "osa"`, via
[`osa_residuals()`](https://afsc-assessments.github.io/Rceattle/reference/osa_residuals.md)
/
[`plot.rceattle_osa()`](https://afsc-assessments.github.io/Rceattle/reference/plot.rceattle_osa.md)).

## Usage

``` r
plot_indexresidual(
  Rceattle,
  file = NULL,
  model_names = NULL,
  species = NULL,
  width = 7,
  height = 6.5,
  residual_type = c("pearson", "osa"),
  line_col = NULL,
  right_adj = 0,
  top_adj = 0.05,
  incl_proj = FALSE,
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

- width, height:

  Saved figure size (inches).

- residual_type:

  `"pearson"` (index residuals on the fleet's own scale) or `"osa"`.

- line_col, right_adj, top_adj, single.plots:

  Deprecated base-graphics arguments, retained for back-compatibility
  and ignored.

- incl_proj:

  Include projection years.

## Value

A `ggplot` object.

## Details

The residual is **observed minus predicted**, taken on the scale the
fleet is fitted on and read from its `Index_distribution`:
`log(observed) - log(predicted)` for a log-scale (`"Lognormal"`) fleet,
and `observed - predicted` for a natural-scale one (`"MVN"`, `"MVNORM"`,
`"Normal"`, `"TruncatedNormal"`), whose sd is absolute. A positive
residual means the survey saw more than the model predicted. Panels hold
different units where a model mixes the two families, which is why the y
scale is free per fleet.

Before 5.9.0 this plotted `predicted - observed`, the negative of what
[`residuals.Rceattle()`](https://afsc-assessments.github.io/Rceattle/reference/residuals.Rceattle.md)
returns for the same fleet. Plots made with an earlier version are
mirrored about zero relative to these.
