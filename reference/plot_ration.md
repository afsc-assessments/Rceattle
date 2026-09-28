# Plot ration

Population-level consumption for ages `minage`+: the individual annual
ration (kg/yr) multiplied by average numbers-at-age and summed over age,
in million mt. This is how the model forms total consumption
(`avgN_at_age * ration`, `predation.hpp`), so it is the consumption that
generates the predation mortality in
[`plot_m2_at_age_prop()`](https://afsc-assessments.github.io/Rceattle/reference/plot_m2_at_age_prop.md)
and the biomass in
[`plot_b_eaten()`](https://afsc-assessments.github.io/Rceattle/reference/plot_b_eaten.md),
plus the other-food term.

## Usage

``` r
plot_ration(
  Rceattle,
  file = NULL,
  minage = 1,
  model_names = NULL,
  line_col = NULL,
  spnames = NULL,
  species = NULL,
  lwd = 3,
  lty = 1,
  right_adj = 0,
  minyr = NULL,
  width = 7,
  height = 6.5,
  incl_proj = FALSE,
  incl_mean = FALSE,
  add_ci = FALSE,
  maxyr = NULL,
  top_adj = 0.15
)
```

## Arguments

- Rceattle:

  A single
  [`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)
  object or a list of them (overlaid).

- file:

  Optional file stem; the figure is written to `<file>_<suffix>.png` if
  given.

- minage:

  Youngest age to sum consumption over, so the figure is "age
  `minage`+". An age on the species' own scale, not an age-bin index; a
  species with no age that old is dropped, with a warning.

- model_names:

  Legend labels for the models.

- line_col:

  Line colours; names, hex, or base-graphics integers. `NULL` uses the
  colorblind-safe Okabe-Ito palette. Applied to whichever variable the
  figure separates by colour, in legend order. Too few colours are
  recycled, with a warning.

- spnames:

  Species labels, length `nspp`. Default: the model's own.

- species:

  Species to include, as indices (`c(1, 3)`), names, a logical mask, or
  `"all"`. Default `NULL` plots every species. Species **labels** belong
  in `spnames`; a character `species` that matches no species name is
  read as labels for back-compatibility.

- lwd:

  Line width on the base-graphics scale: the default `3` renders as a
  standard-weight ggplot line. A vector varies it across series.

- lty:

  Line type. A vector varies it across the levels of whatever the figure
  separates by line type.

- right_adj:

  Ignored. Base-graphics leftover: the figure widened its right margin
  to fit the legend. Set margins on the returned ggplot instead.

- minyr, maxyr:

  First / last year to plot.

- width, height:

  Saved figure size in inches.

- incl_proj:

  Include the projection years, with a dashed divider at the last
  hindcast year.

- incl_mean:

  Add a horizontal line at the hindcast mean of each series.

- add_ci:

  Add a 95% confidence interval. Only available where the plotted
  quantity holds standard errors; warns and draws none otherwise.

- top_adj:

  Ignored. Base-graphics leftover; see `right_adj`.

## Details

Colour separates the models; line type separates the sexes.
