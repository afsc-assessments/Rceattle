# Likelihood components along a profile

Pulls the per-fleet and per-species negative log-likelihood components
out of every fit in a
[`profile()`](https://rdrr.io/r/stats/profile.html) and returns them as
one long data frame, one row per grid point per component.

The total is the least informative curve on a profile: a well-behaved
quadratic total can hide two data sources pulling the parameter in
opposite directions, and their conflict shows only once the components
are drawn separately. This is the extractor behind
[`plot_profile()`](https://afsc-assessments.github.io/Rceattle/reference/plot_profile.md);
use it to tabulate, or to draw the figure yourself.

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
  default) or the unweighted ones (`FALSE`). See Details.

- relative:

  How to place each series on the y axis: `"own"` (default), `"scaled"`,
  `"minimum"` or `"none"`. See Details.

- minfraction:

  Drop components whose change over the grid is less than this fraction
  of the total's change, as in `r4ss::SSplotProfile()`. Default `0`
  keeps everything;
  [`plot_profile()`](https://afsc-assessments.github.io/Rceattle/reference/plot_profile.md)
  uses `0.01`, so a total moving 60 units cuts at 0.6. `"Total"` is
  never dropped.

- include_total:

  Include the `"Total"` series. Default `TRUE`; `FALSE` also disables
  `minfraction`, which is defined against the total.

## Value

A data frame with one row per grid point per retained component: the
profile's `grid` columns (`slot_1`, ...), then `fit` (grid row index),
`component` (the `jnll_comp` row), `unit` (fleet or species name, `NA`
for model-wide rows), `axis` (`"fleet"`, `"species"` or `"model"`),
`series` (the plotting label) and `value` (the re-zeroed negative
log-likelihood). Series are ordered by decreasing change over the grid,
`"Total"` first. The profile's `param`, `alias` and the `relative` used
are attributes.

## Details

**Which cells are reported.** `jnll_comp`'s columns mean different
things on different rows – fleets on the data, selectivity and
catchability rows, species on the priors, penalties and predation rows –
so each cell is labelled from the axis its row uses, giving e.g.
`"Shelikof acoustic: Index data"`. The unit is dropped from the label
when the model has only one fleet (or species) to distinguish. Cells
that are zero at every grid point are dropped as components the model
does not fit.

**What `Total` is.** The `"Total"` series is `object$nll`, the objective
each grid fit minimized. Under `random_rec = TRUE` that is the
Laplace-approximated marginal likelihood while the components are the
inner joint negative log-likelihood, so they will not sum to it, and the
function says so when they differ. Compare the shapes, not the sums.

**Non-converged grid points** keep their row with a value of `NA`, so a
failed fit leaves a gap in the curve rather than a segment drawn across
it.

**Weighted or unweighted components** (`weighted`). Conflict is normally
read off the weighted components, since those are what moved the fit.
`unweighted_jnll_comp` scores the composition, CAAL and stomach-content
rows without their data weights (`Comp_weights`, `CAAL_weights`,
`Diet_comp_weights`); every other row equals `jnll_comp`, so
`weighted = FALSE` changes only those three components. A
Dirichlet-multinomial row is scored at theta = 1, an effective sample
size of about (1 + N)/2 rather than N, so it is not comparable with a
multinomial row.

**Placing each series on the y axis** (`relative`).

- `"own"`:

  (default) subtract each series' own minimum, so every curve starts at
  zero and its minimum marks the value that component prefers. Objective
  units, so depth still says how strongly.

- `"scaled"`:

  as `"own"`, then divide by each series' own change over the grid, so
  every curve runs 0 to 1. Use it where one component dwarfs the rest
  and flattens them onto the axis. It **discards magnitude** – a
  component moving 0.02 draws like one moving 40 – so raise
  `minfraction` with it. A series that does not move stays at zero.

- `"minimum"`:

  subtract each series' value at the grid point where the total is
  lowest, showing each component's change away from the fitted optimum.

- `"none"`:

  the raw negative log-likelihoods.

`minfraction` and the series ordering are always computed on the raw
change over the grid, before any of this, so `"scaled"` cannot quietly
disable the filter by making every span 1.

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
#> 1   0.20   1     Total <NA> model  Total 28.479709
#> 2   0.25   2     Total <NA> model  Total 23.399627
#> 3   0.30   3     Total <NA> model  Total 18.451669
#> 4   0.35   4     Total <NA> model  Total 13.636791
#> 5   0.40   5     Total <NA> model  Total  8.955937
#> 6   0.45   6     Total <NA> model  Total  4.410035
# }
```
