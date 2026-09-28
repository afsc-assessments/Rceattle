# Likelihood profile across one or more parameter cells

Re-fits an Rceattle model while holding selected cells of a parameter
fixed at user-specified values. Supports profiling a single cell (e.g.
`R_log_sd[species = 1]`) and arbitrary N-dimensional cross-profiles over
multiple cells, e.g. `log_M1[1, 1, 1]` and `log_M1[1, 2, 1]` jointly, to
profile residual M for males against females. For each grid point the
targeted cells are fixed in the TMB map and the remaining parameters are
re-estimated; the result is a grid of Rceattle models for downstream NLL
surfaces.

## Usage

``` r
# S3 method for class 'Rceattle'
profile(
  fitted = NULL,
  param = NULL,
  slots = NULL,
  values = NULL,
  transform = "log",
  cores = NULL,
  getsd = NULL,
  joint = c("none", "multiply", "add", "value"),
  ...
)
```

## Arguments

- fitted:

  an Rceattle model fit using
  [`fit_mod`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)

- param:

  Name of the parameter to profile. Two ways to specify it:

  Raw parameter slot

  :   any name in `Rceattle$estimated_params`; tested for `"R_log_sd"`,
      `"rec_pars"`, and `"log_M1"`. `slots` must index into the full
      array and `transform` controls the scale.

  Natural-scale alias

  :   convenience shortcut for the three documented parameters. Aliases
      imply `transform = "log"` (values are taken in natural units and
      log'd before being substituted) and, for `rec_pars`, fill in the
      column from the alias name so `slots` only needs the species
      index:

      - `"sigmaR"`, `"R_sd"` -\> `R_log_sd`

      - `"M1"` -\> `log_M1`

      - `"R0"` -\> `rec_pars[, 1]`

      - `"alpha"` -\> `rec_pars[, 2]`

      - `"beta"` -\> `rec_pars[, 3]`

      - `"q"` -\> `index_log_q`, base catchability, one entry per fleet.
        `slots` takes a fleet index or a **fleet name**. Fleets sharing
        a `Catchability_index` share one q, so naming any of them
        profiles the whole group: the rest are added with a message and
        `joint` becomes `"value"`.

      If `transform` is supplied with an alias it is ignored (with a
      warning).

- slots:

  A list whose entries are integer index vectors, one entry per cell to
  fix. Each entry's length must equal the number of dimensions of the
  resolved parameter, 1 for vectors (`R_log_sd`), 2 for matrices
  (`rec_pars`), 3 for 3-D arrays (`log_M1`). When using the
  `"R0"`/`"alpha"`/`"beta"` aliases, supply only the species index
  (length 1); the column is filled in from the alias. E.g.
  `list(c(1, 2, 1))` fixes `log_M1[1, 2, 1]`;
  `list(c(1, 1, 1), c(1, 2, 1))` fixes both sex cells for a
  males-vs-females cross-profile of species 1; `list(1, 2)` with
  `param = "sigmaR"` cross-profiles species 1 and 2. If omitted,
  defaults to a single species-1 slot shaped to match the resolved
  parameter (e.g. `list(1)` for `R_log_sd`, `list(c(1, 1, 1))` for
  `log_M1`, `list(1)` for the `rec_pars` aliases) and emits a warning;
  pass `slots` explicitly to silence the warning. Defaulting requires
  `length(values) == 1L` (otherwise the user must explicitly say which
  cell each grid targets).

- values:

  A list of numeric vectors, one per entry of `slots`. The full grid of
  fits is `expand.grid(values)`, so a single slot gives a 1-D profile
  and *k* slots give a *k*-D cross-profile.

- transform:

  How to map user values onto the internal parameter scale before
  substituting them into `inits`. Either `"log"` (default),
  `"identity"`, or a unary function (e.g. `qlogis`). Applied
  element-wise to every grid value. Aliases override this with `"log"`.

- cores:

  Number of cores to use for parallel fits. Default `NULL` picks
  `parallel::detectCores() - 6`, capped at 2 when running under
  `R CMD check` (which sets `_R_CHECK_LIMIT_CORES_`). Set to 1 to force
  sequential execution.

- getsd:

  whether each grid fit runs
  [`TMB::sdreport`](https://rdrr.io/pkg/TMB/man/sdreport.html). The
  profile reads only the objective (`nll`), so `FALSE` is faster with no
  effect on the profile. Default `NULL` inherits the input model's
  setting (`TRUE` only if it holds an `sdrep`).

- joint:

  How the grid moves the cells in `slots`. Use a joint mode to profile a
  parameter supplied as a vector, such as an empirical age-based M:
  under `"none"` ten ages over 13 values is `13^10` fits, not 13.

  `"none"`

  :   (default) each slot gets its own grid and they cross, giving a
      *k*-dimensional profile for *k* slots.

  `"multiply"`

  :   one grid, multiplying every cell's current value. The schedule
      keeps its shape and moves up or down as a whole. `1` is the fitted
      model.

  `"add"`

  :   one grid, added to every cell. The schedule keeps the differences
      between cells instead. `0` is the fitted model.

  `"value"`

  :   one grid; every cell takes that value.

  A joint mode needs `values` to be one vector. `"multiply"` and `"add"`
  work in natural units, so `transform` must be `"log"` or `"identity"`.
  A grid value that would take a log-scale parameter to zero or below is
  an error naming the cell, not a failed fit.

- ...:

  Unused; present for consistency with the
  [`stats::profile`](https://rdrr.io/r/stats/profile.html) generic.

## Value

A list with elements:

- Rceattle_list:

  list of fitted Rceattle models, one per grid row; entries for
  non-converged fits are `NULL` so positions stay aligned with `grid`.

- grid:

  data frame of grid values on the user scale (before `transform`); one
  column per profiled cell, named `slot_1`, `slot_2`, ...

- nll:

  numeric vector of joint negative log-likelihoods (`opt$objective`);
  `NA` where the fit did not converge.

- param:

  the profiled parameter name (echoed).

- slots:

  the slots list (echoed for downstream plotting).

- alias:

  the name you asked for, when it was one of the natural-scale aliases.
  Profiling `"M1"` returns `param = "log_M1"`, because the model
  estimates M on the log scale, while `grid` holds the M values you
  supplied. `alias` keeps `"M1"`, so a figure can label the axis in the
  units profiled rather than in log units. `NA` if you named the
  parameter slot directly.

- joint:

  the `joint` mode used, so a figure and
  [`print()`](https://rdrr.io/r/base/print.html) can say that `grid`
  holds a multiplier or an offset rather than a parameter value.

Has class `"Rceattle_profile"`, so printing it reports whether the grid
brackets the minimum; see
[`print.Rceattle_profile`](https://afsc-assessments.github.io/Rceattle/reference/print.Rceattle_profile.md).
Every element indexes exactly as before.

## Examples

``` r
# \donttest{
data(BS2017SS)
ss_run <- fit_mod(data_list = BS2017SS,
    inits = NULL, file = NULL,
    estimateMode = 0, random_rec = FALSE,
    msmMode = 0, avgnMode = 0,
    phase = FALSE, verbose = 0)
#> Warning: Passing ‘phase’, ‘verbose’ directly to fit_mod() is deprecated and will be removed in a future release. Bundle these into fit_control() instead, e.g. fit_control(phase = ..., verbose = ...). Forwarding for now.
#> `age_trans_matrix` data does not span range of age for species 1 will fill with 0s

# 1-D profile of sigmaR for species 1 (alias form, natural scale)
p1 <- profile(ss_run,
    param  = "sigmaR",
    slots  = list(1),
    values = list(seq(0.1, 1.5, by = 0.1)))

# Equivalent raw form (log scale, user does the transform)
p1_raw <- profile(ss_run,
    param     = "R_log_sd",
    slots     = list(1),
    values    = list(log(seq(0.1, 1.5, by = 0.1))),
    transform = "identity")

# 2-D cross-profile of M1 across species 1 and 2 (sex 1, age 1).
# BS2017SS is single-sex; with a multi-sex model the same form
# (e.g. c(1, 1, 1), c(1, 2, 1)) would cross-profile males vs females.
p2 <- profile(ss_run,
    param  = "M1",
    slots  = list(c(1, 1, 1), c(2, 1, 1)),
    values = list(seq(0.1, 0.4, length.out = 3),
                  seq(0.1, 0.4, length.out = 3)))

# 1-D profile of SRR alpha for species 1 (alias drops the rec_pars column)
p3 <- profile(ss_run,
    param  = "alpha",
    slots  = list(1),
    values = list(seq(2, 80, length.out = 20)))

# An empirical age-based M, moved as a whole. One slot per age bin of
# species 1, all on ONE grid of multipliers, so the shape of the schedule is
# preserved and 1 is the fitted model. Ages run minage .. minage + nages - 1.
sp <- 1
p4 <- profile(ss_run,
    param  = "M1",
    slots  = lapply(seq_len(ss_run$data_list$nages[sp]),
                    function(a) c(sp, 1, a)),
    values = list(seq(0.6, 1.4, by = 0.05)),
    joint  = "multiply")

# The same schedule shifted by a constant instead of scaled; 0 is the fit.
p5 <- profile(ss_run,
    param  = "M1",
    slots  = lapply(seq_len(ss_run$data_list$nages[sp]),
                    function(a) c(sp, 1, a)),
    values = list(seq(-0.1, 0.1, by = 0.02)),
    joint  = "add")

# Base catchability for a fleet, named rather than counted. If it shares a
# Catchability_index the whole group moves with it.
p6 <- profile(ss_run,
    param  = "q",
    slots  = list("BT_Pollock"),
    values = list(seq(0.5, 2.0, by = 0.1)))
# }
```
