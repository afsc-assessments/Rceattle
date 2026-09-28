# Jitter analysis

Refits the Rceattle model from starting values perturbed by N(0, sd)
around the model's initial (pre-fit) parameters, to check convergence
robustness.

## Usage

``` r
jitter(
  object = NULL,
  njitter = 50,
  sd = 0.2,
  phase = FALSE,
  seed = 123,
  cores = NULL,
  getsd = NULL,
  timeout = Inf,
  fit_control = NULL,
  Rceattle = NULL
)
```

## Arguments

- object:

  an Rceattle model fit using
  [`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md).

- njitter:

  the number of jitters to run

- sd:

  standard deviation for jitter (default = 0.2)

- phase:

  as in
  [`fit_mod`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)
  default = FALSE. Jitters restart from perturbed *starting* values, so
  a model that needed phasing to fit its real data needs it here too;
  leave this at `FALSE` for such a model and the jitters end far from
  any optimum and are dropped as non-converged, which reads as
  multimodality rather than as an unphased fit.

- seed:

  random number seed. Each jitter `i` uses `seed + i` so results are
  reproducible under both sequential and parallel execution.

- cores:

  Number of cores for the parallel refits. Default `NULL` picks
  `parallel::detectCores() - 6`, capped at 2 under `R CMD check` (which
  sets `_R_CHECK_LIMIT_CORES_`). Set to 1 to force sequential execution.

- getsd:

  whether each jitter runs
  [`TMB::sdreport`](https://rdrr.io/pkg/TMB/man/sdreport.html). Jitter
  compares objectives and point estimates across starts, so `FALSE` is
  faster with no effect on that comparison. Default `NULL` inherits the
  input model's setting (`TRUE` only if it holds an `sdrep`).

- timeout:

  elapsed-second limit per jitter, `Inf` (default) for none. A jitter is
  a deliberately perturbed start and the optimizer runs with no
  iteration cap, so this is the diagnostic most likely to send one
  somewhere pathological and stall the whole run, a hang no convergence
  check can catch, because the fit never returns. One that exceeds the
  limit is stopped, counted as non-converged and reported separately.
  Approximate: the limit is checked when control returns to R, so it
  fires between the optimizer's function evaluations rather than inside
  one.

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

a list of 1. `Rceattle_list`, the converged jitters, and 2. `nll`, their
objective values. Non-converged (or timed-out) starts are dropped and
reported in a message, so both can be shorter than `njitter`, and that
count is itself the result, since the whole point is what fraction of
random starts reach the same optimum.

## Note

Attaching Rceattle masks [`jitter`](https://rdrr.io/r/base/jitter.html),
so `jitter(x)` on a numeric vector reaches this function and reports a
missing model rather than adding noise to `x`. Call
[`base::jitter()`](https://rdrr.io/r/base/jitter.html) explicitly for
the base-graphics behaviour.

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
jitters <- jitter(ss_run, njitter = 10)
# }
```
