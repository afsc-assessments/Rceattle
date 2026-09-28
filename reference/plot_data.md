# Timeline of data used in the model likelihoods

Plots a timeline of presence/absence (and optionally relative quantity)
of the data sources contributing to an Rceattle model's likelihood, by
year and fleet. Modeled after Stock Synthesis's `r4ss::SSplotData()`.

## Usage

``` r
plot_data(
  Rceattle,
  file = NULL,
  subplots = 1:2,
  datatypes = "all",
  fleets = "all",
  species = "all",
  ghost = FALSE,
  fleetcol = "default",
  width = 7,
  height = 5,
  res = 200,
  ptsize = 10,
  margins = c(5.1, 2.1, 2.1, 10.1),
  cex = 2,
  lwd = 12,
  maxsize = 1,
  alphasize = 1,
  mainTitle = FALSE,
  cex.main = 1
)
```

## Arguments

- Rceattle:

  Either a single Rceattle model object exported from
  [`Rceattle::fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md),
  or a raw `data_list` (e.g. one of the bundled datasets such as
  `BS2017MS`).

- file:

  Path/prefix used to save the figure. If `NULL` (default) the plot is
  drawn to the active device only. When supplied, the figure is written
  to `paste0(file, "_data_plot.png")` (and `_data_plot2.png` for the
  bubble subplot).

- subplots:

  Integer vector controlling which subplots are produced:

  - 1, equal-size points showing presence/absence by year/fleet

  - 2, points scaled to relative quantity / precision within each data
    type (catch tonnage, 1/SE for indices, sample size for comps)

- datatypes:

  Either `"all"` or a subset of
  `c("catch", "index", "agecomp", "lencomp", "caal", "diet")`.

- fleets:

  Either `"all"` or an integer vector of `Fleet_code`s to include. Diet
  "fleets" are encoded after the regular fleets (one per predator
  species).

- species:

  Either `"all"` or an integer vector of species indices to include.

- ghost:

  Logical. If `TRUE`, also include rows excluded from the likelihood
  (`Year < 0`). Defaults to `FALSE`.

- fleetcol:

  Either `"default"` or a vector of colors (one per fleet shown after
  subsetting). Diet entries get an additional color appended.

- width, height:

  Figure dimensions in inches.

- res:

  Resolution (dpi) for the saved PNG.

- ptsize:

  Pointsize passed to [`png()`](https://rdrr.io/r/grDevices/png.html).

- margins:

  `par("mar")` for the plot. Increase the right margin if long fleet
  names are clipped.

- cex:

  Character expansion for points showing isolated years.

- lwd:

  Line width for runs of consecutive years (subplot 1).

- maxsize:

  Max bubble radius (in plot units) for subplot 2.

- alphasize:

  Bubble fill transparency (0-1).

- mainTitle:

  Logical; if `TRUE` add a default title.

- cex.main:

  Title character expansion.

## Value

Invisibly, a list with `typetable`, the long data frame underlying the
plot (year, fleet, data type, relative size).

## Details

Years with `Year > 0` contribute to the likelihood; rows with `Year < 0`
(kept for fitting comparison only) are shown when `ghost = TRUE`. Diet
data uses `Year == 0` for the fixed reference stomach composition; those
entries are plotted at the model `endyr` and labeled per predator
species.
