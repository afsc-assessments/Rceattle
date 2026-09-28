# Plot predation mortality by age and predator

Share of the predation mortality (M2) on each prey age that is
attributable to each predator species. The shares sum to 1 across
predators for every prey age and year; a year with no predation on that
age leaves them undefined and draws nothing.

## Usage

``` r
plot_m2_at_age_prop(
  Rceattle,
  file = NULL,
  age = 1,
  model_names = NULL,
  line_col = NULL,
  spnames = NULL,
  species = NULL,
  lwd = 3,
  right_adj = 0,
  top_adj = 0.15,
  minyr = NULL,
  width = 7,
  height = 6.5,
  incl_proj = FALSE,
  incl_mean = FALSE,
  add_ci = FALSE,
  maxyr = NULL,
  lty = 1
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

- age:

  Prey age to plot the M2 proportions at, on the prey species' own age
  scale rather than as an age-bin index. A prey species that has no such
  age is dropped, with a warning.

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

- right_adj:

  Ignored. Base-graphics leftover: the figure widened its right margin
  to fit the legend. Set margins on the returned ggplot instead.

- top_adj:

  Ignored. Base-graphics leftover; see `right_adj`.

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

- lty:

  Line type. A vector varies it across the levels of whatever the figure
  separates by line type.

## Details

Colour separates predators; line type separates models. Panels are prey
species (and sex, where the prey is sexed).
