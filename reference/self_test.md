# Self-test simulation analysis

Simulates data from an Rceattle model and refits the model to the
simulated data, to check that the fitting procedure recovers the
operating-model parameters. Pass `process` to redraw a process as well
as the observations.

## Usage

``` r
self_test(
  object = NULL,
  nsim = 50,
  simulate = TRUE,
  seed = 123,
  cores = NULL,
  getsd = NULL,
  phase = NULL,
  start = c("initial", "estimated"),
  debug = FALSE,
  timeout = Inf,
  process = FALSE,
  fit_control = NULL,
  Rceattle = NULL
)
```

## Arguments

- object:

  an Rceattle model fit using
  [`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md).

- nsim:

  number of simulations

- simulate:

  passed to
  [`sim_mod`](https://afsc-assessments.github.io/Rceattle/reference/sim_mod.md).
  `TRUE` (default) simulates data with observation error; `FALSE` uses
  the model's expected values.

- seed:

  random number seed. Simulation `i` uses `seed + i`, so results are
  reproducible under sequential and parallel execution alike.

- cores:

  Number of cores for the parallel refits. Default `NULL` picks
  `parallel::detectCores() - 6`, capped at 2 under `R CMD check` (which
  sets `_R_CHECK_LIMIT_CORES_`). Set to 1 to force sequential execution.

- getsd:

  whether each refit runs
  [`TMB::sdreport`](https://rdrr.io/pkg/TMB/man/sdreport.html).
  Self-test compares refit point estimates to the operating model, so
  `FALSE` is faster with no effect on that comparison. `NULL` (default)
  inherits the input model's setting (`TRUE` only if it holds an
  `sdrep`).

- phase:

  as in
  [`fit_mod`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md).
  Under the default `start = "initial"` each refit covers the ground the
  original fit did, so a model that needed phasing for its real data
  needs it for every simulated one; without it such a model's refits can
  end many orders of magnitude from a zero gradient and be dropped as
  non-converged. `NULL` (default) reads what
  [`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)
  recorded on the source fit, so a model fitted under the package
  default `phase = FALSE` is refitted unphased.

- start:

  Which parameter set each refit starts from: `"initial"` (default) or
  `"estimated"`. See Details.

- debug:

  return every simulation rather than only the converged ones. The
  dropped runs are the interesting ones when a self-test comes back
  short, and each holds its own `$convergence` diagnostics. See
  **Value**.

- timeout:

  elapsed-second limit per simulation, `Inf` (default) for none. The
  optimizer runs with no iteration cap, so a replicate that wanders
  somewhere pathological can stall the run – a hang no convergence check
  can catch, because the fit never returns. One that exceeds the limit
  is stopped, counted as non-converged and reported separately. It is
  approximate: the limit is checked when control returns to R, so it
  fires between the optimizer's function evaluations rather than inside
  one.

- process:

  passed to
  [`sim_mod`](https://afsc-assessments.github.io/Rceattle/reference/sim_mod.md).
  `FALSE` (default) keeps the fitted process deviations, so the test
  measures whether the estimator recovers its own parameters from new
  observations. Naming a process (`"recruitment"`, `"M"`, `"growth"`,
  `"dynamics"`, `TRUE`) redraws it, so the test instead measures
  recovery of a process the estimator has not been shown.

- fit_control:

  optional
  [`fit_control()`](https://afsc-assessments.github.io/Rceattle/reference/fit_control.md)
  bundle for the refits. Only `phase` and `getsd` are read; see **What
  [`fit_control()`](https://afsc-assessments.github.io/Rceattle/reference/fit_control.md)
  reaches**.

- Rceattle:

  deprecated name for `object`, still accepted so existing scripts keep
  working. Supplying both is an error.

## Value

A list of Rceattle models named `Sim_1`, `Sim_2`, ... By default only
the converged simulations, renumbered contiguously, with a message
reporting how many were dropped. The list has class
`"Rceattle_selftest"` and carries the number attempted in
`attr(, "nsim")`, so printing it reports the convergence rate; it is
otherwise the plain list it always was.

With `debug = TRUE`, every simulation, `Sim_i` being simulation `i` so
it pairs with seed `seed + i`, plus the convergence verdicts in
`attr(, "converged")`. A refit that errored outright comes back as the
condition object rather than a model, so it cannot abort the run.

Where `process` redrew something, `attr(, "process_sim")` holds the
deviations that generated each replicate, keyed by the same `Sim_i`
names and subset alongside the models. Each entry names whichever of
`rec_dev`, `init_dev`, `log_M1_dev` and `beta_linkage_re` were drawn,
each with a same-shaped `_drawn` logical marking the cells the draw
touched. Restrict any recovery statistic to those cells, and compare
estimates against these rather than against the operating model, whose
fitted deviations are no longer what generated the data.

## Details

Each simulation refits under the fit's own map, so a parameter fixed
there keeps its value and is not estimated. A projection-only fit
(`estimateMode = "Projection"`) refits under a rebuilt map.
`self_test()` warns if `process` redraws a deviation that the fit's map
fixes and
[`build_map()`](https://afsc-assessments.github.io/Rceattle/reference/build_map.md)
would estimate: the refit holds it at the fit's value, so cannot recover
the draw.

**What `start` costs.** `"initial"` uses `initial_params`, so the
estimator travels the same distance to the optimum on simulated data as
it did on the real data. `"estimated"` starts from the generating values
instead – the fixed effects, and under `random_rec = TRUE` the inner
Laplace problem too – so it is much faster and far likelier to converge,
but on a multimodal or weakly identified surface the optimizer never
leaves the basin containing them, so recovery is close to guaranteed by
construction. Read it as optimistic about recovery, not merely less
powerful. It is not a complete warm start either –
[`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)
resets `log_Ftarget`, `proj_F_prop` and the stock-recruit
\\\alpha\\/\\\beta\\ under both settings – and non-identifiability stays
visible either way through `$convergence`'s Hessian conditioning and
estimability checks. What a warm start stops testing is *reachability*.

## Interpreting the spread

[`sim_mod`](https://afsc-assessments.github.io/Rceattle/reference/sim_mod.md)
redraws the observations only, and warns about each row it leaves alone:
a predator with empirical rather than estimated suitability has no
predicted diet to draw from, and a covariance (MVN/MVNORM) survey fleet
has no covariance outside its fitted years. Those data are identical in
every replicate, so recovery of whatever they inform is optimistic.

By default recruitment is not redrawn, so under `random_rec = TRUE`
every replicate shares the operating model's single realization – its
shrunk empirical-Bayes modes, not a draw from N(0, sigmaR). The spread
across replicates is then observation error alone and a lower bound on
estimation uncertainty in SSB and recruitment, so do not read it against
the model's own uncertainty bands, which include process error; and
sigmaR is re-estimated from deviations shrunk toward zero the same way
every time, a downward bias averaging does not remove.
`process = "recruitment"` removes both, at the cost of asking a
different question. See
[`vignette("model-diagnostics")`](https://afsc-assessments.github.io/Rceattle/articles/model-diagnostics.md).

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
sims <- self_test(ss_run, nsim = 10)
# }
```
