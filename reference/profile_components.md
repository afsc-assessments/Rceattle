# Likelihood components along a profile

Pulls the per-fleet and per-species negative log-likelihood components
out of every fit in a
[`profile()`](https://rdrr.io/r/stats/profile.html) and returns them as
one long data frame, one row per grid point per component.

The total is the least informative curve on a profile: a well-behaved
quadratic total can hide two data sources pulling the parameter in
opposite directions, and their conflict is visible only once the
components are drawn separately. This is the extractor behind
[`plot_profile()`](https://afsc-assessments.github.io/Rceattle/reference/plot_profile.md);
use it directly to tabulate, or to draw the figure yourself.

## Usage

``` r
profile_components(
  object,
  weighted = TRUE,
  relative = c("own", "scaled", "minimum", "none"),
  minfraction = 0,
  include_total = TRUE
)
```

## Arguments

- object:

  An `"Rceattle_profile"` object from
  [`profile.Rceattle()`](https://afsc-assessments.github.io/Rceattle/reference/profile.Rceattle.md).

- weighted:

  Report the weighted components the optimizer minimized (`TRUE`, the
  default, `quantities$jnll_comp`) or the unweighted ones (`FALSE`,
  `quantities$unweighted_jnll_comp`). Conflict is normally read off the
  weighted components, since those are what moved the fit.
  `unweighted_jnll_comp` exists so Francis and McAllister-Ianelli can
  read a composition likelihood without its `Comp_weights` multiplier,
  so only the rows that hold such a multiplier are filled: composition,
  CAAL, stomach content and the two linkage rows. Every other row is
  zero there and is dropped as unfitted, so `weighted = FALSE` returns a
  much smaller set of series, the index, catch, selectivity,
  catchability and penalty components are absent, not flat.

- relative:

  How to place each series on the y axis.

  `"own"`

  :   (default) subtract each series' own minimum, so every curve starts
      at zero and its minimum marks the value that component prefers.
      Objective units, so depth still says how strongly.

  `"scaled"`

  :   as `"own"`, then divide each series by its own change over the
      grid, so every curve runs 0 to 1. Use it when one component dwarfs
      the rest and flattens the others onto the axis. It **discards
      magnitude**: a component moving 0.02 draws like one moving 40, so
      raise `minfraction` with it. A series that does not move stays at
      zero.

  `"minimum"`

  :   subtract each series' value at the grid point where the total is
      lowest, so the curves show each component's change away from the
      fitted optimum.

  `"none"`

  :   the raw negative log-likelihoods.

  `minfraction` and the series ordering are always computed on the raw
  change over the grid, before any of this is applied, so `"scaled"`
  does not quietly disable the filter by making every span 1.

- minfraction:

  Drop components whose change over the grid is less than this fraction
  of the total's change, as in `r4ss::SSplotProfile()`. Default `0`
  keeps everything;
  [`plot_profile()`](https://afsc-assessments.github.io/Rceattle/reference/plot_profile.md)
  uses `0.01`. `"Total"` is never dropped.

- include_total:

  Include the `"Total"` series. Default `TRUE`; `FALSE` also disables
  `minfraction`, which is defined against the total.

## Value

A data frame with one row per grid point per retained component: the
profile's `grid` columns (`slot_1`, ...) holding the value profiled
over, then `fit` (grid row index), `component` (the `jnll_comp` row),
`unit` (fleet or species name, `NA` for model-wide rows), `axis`
(`"fleet"`, `"species"` or `"model"`), `series` (the plotting label),
and `value` (the re-zeroed negative log-likelihood). Series are ordered
by decreasing change over the grid, with `"Total"` first. The profile's
`param`, `alias` and the `relative` used are held as attributes.

## Details

**Which cells are reported.** `jnll_comp` is a component-by-column
matrix whose columns mean different things on different rows, fleets on
the data, selectivity and catchability rows, species on the priors,
penalties and predation rows. Each cell is labelled from the axis its
row uses, so a cell becomes e.g. `"Shelikof acoustic: Index data"`. The
unit is dropped from the label when the model has only one fleet (or one
species) to distinguish. Cells that are zero at every grid point are
dropped: they are components the model does not fit.

**What `Total` is.** The `"Total"` series is `object$nll`, the objective
each grid fit actually minimized. Under `random_rec = TRUE` that is the
Laplace-approximated marginal likelihood while the components are the
inner joint negative log-likelihood, so the components will not sum to
it; the function says so when they differ. Compare the shapes, not the
sums.

**Non-converged grid points** keep their row with a value of `NA`, so a
failed fit leaves a gap in the curve rather than a straight segment
drawn across it.

## See also

[`plot_profile()`](https://afsc-assessments.github.io/Rceattle/reference/plot_profile.md)
to draw it,
[`profile.Rceattle()`](https://afsc-assessments.github.io/Rceattle/reference/profile.Rceattle.md)
to produce the profile.

## Examples

``` r
# \donttest{
data(BS2017SS)
ss_run <- fit_mod(data_list = BS2017SS, inits = NULL, file = NULL,
    estimateMode = 0, random_rec = FALSE, msmMode = 0, avgnMode = 0,
    phase = FALSE, verbose = 0)
#> Warning: Passing ‘phase’, ‘verbose’ directly to fit_mod() is deprecated and will be removed in a future release. Bundle these into fit_control() instead, e.g. fit_control(phase = ..., verbose = ...). Forwarding for now.
#> `age_trans_matrix` data does not span range of age for species 1 will fill with 0s

prof <- profile(ss_run, param = "M1", slots = list(c(1, 1, 1)),
                values = list(seq(0.2, 0.5, by = 0.05)))

comps <- profile_components(prof)
head(comps)
#>   slot_1 fit component unit  axis series     value
#> 1   0.20   1     Total <NA> model  Total 28.479708
#> 2   0.25   2     Total <NA> model  Total 23.399627
#> 3   0.30   3     Total <NA> model  Total 18.451669
#> 4   0.35   4     Total <NA> model  Total 13.636791
#> 5   0.40   5     Total <NA> model  Total  8.955937
#> 6   0.45   6     Total <NA> model  Total  4.410035
# }
```
