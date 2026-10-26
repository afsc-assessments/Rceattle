# Survey index residuals

Plots survey index residuals by year, faceted by survey fleet
(`residual_type = "pearson"`, the default), or one-step-ahead (OSA)
residual diagnostics for a single fit (`residual_type = "osa"`, via
[`osa_residuals()`](https://afsc-assessments.github.io/Rceattle/reference/osa_residuals.md)
/
[`plot.rceattle_osa()`](https://afsc-assessments.github.io/Rceattle/reference/plot.rceattle_osa.md)).

## Usage

``` r
plot_index_residual(
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

  `"pearson"` (standardized index residuals) or `"osa"`.

- line_col, right_adj, top_adj, single.plots:

  Deprecated base-graphics arguments, retained for back-compatibility
  and ignored.

- incl_proj:

  Include projection years.

## Value

A `ggplot` object.

## Details

The residual is observed minus predicted, standardized by the sd the row
was fitted with – the same quantity
[`residuals.Rceattle()`](https://afsc-assessments.github.io/Rceattle/reference/residuals.Rceattle.md)
returns for `type = "pearson"`. On a log-scale fleet it is measured
against the likelihood's centre,
`log(hat) - bias_adjust_obs * sd^2 / 2`; on a natural-scale fleet
(`"MVN"`, `"MVNORM"`, `"Normal"`, `"TruncatedNormal"` in
`Index_distribution`), whose sd is absolute, it is `(obs - hat) / sd`. A
covariance fleet (`"MVN"`, `"MVNORM"`) takes that sd from
`sqrt(diag(index_cov))`, which its likelihood actually scores, so the
residual is marginal rather than whitened. A positive residual means the
survey saw more than the model predicted. The y scale stays free per
fleet, since residual spread differs by fleet.
