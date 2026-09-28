# Plot M1 or M2 at age

Mortality-at-age over the hindcast for one model: predation mortality
(`M2`, the default) or residual natural mortality (`M1`). One component
at a time, never their sum,
[`plot_m_at_age()`](https://afsc-assessments.github.io/Rceattle/reference/plot_m_at_age.md)
draws total M (M1 + M2) as a time series.

## Usage

``` r
plot_mortality(
  Rceattle,
  file = NULL,
  incl_proj = FALSE,
  zlim = NULL,
  type = "lines",
  width = 8,
  height = 5.5,
  title = NULL,
  log = FALSE,
  minyr = NULL,
  theta = 155,
  species = NULL,
  maxage = NULL,
  title_cex = 10,
  M2 = TRUE
)
```

## Arguments

- Rceattle:

  A single
  [`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)
  object. A list of more than one is an error; call it per model.

- file:

  name of a file to identified the files exported by the function.

- incl_proj:

  Include the projection years (TRUE/FALSE)

- zlim:

  Ignored. Base-graphics leftover: the fill range of the tile plot. Use
  `p + ggplot2::scale_fill_viridis_c(limits = ...)`.

- type:

  `"lines"` (default) draws M-at-age with one line per year; `"heatmap"`
  (or `0`) draws the same series as an age-by-year tile plot. Any other
  value gives the lines.

- width:

  Plot width when saved "inches"

- height:

  Plot height when saved "inches"

- title:

  Additional title to add. Will also add species names if not NULL

- log:

  TRUE/FALSE plot the series on a log scale

- minyr:

  First year to plot

- theta:

  Ignored. Base-graphics leftover: viewing angle of a `persp` surface,
  which is no longer drawn.

- species:

  Species to plot. Plots all if null.

- maxage:

  Oldest age to draw, on the species' own age scale rather than as a
  count of age bins. `NULL` draws every age.

- title_cex:

  Ignored. Base-graphics leftover: title font size. Use
  `p + ggplot2::theme(plot.title = ggplot2::element_text(size = ...))`.

- M2:

  Draw predation mortality M2 (TRUE, the default) or residual natural
  mortality M1 (FALSE). Neither is the sum of the two.

## Value

A `ggplot`.
