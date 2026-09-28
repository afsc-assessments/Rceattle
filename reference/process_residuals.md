# Process residuals for an Rceattle model's random-effect processes

One-sample process residuals for the model's random-effect (or
penalized) process deviations, in the style of SAM's `procres()`
(Nielsen and Berg 2014). They validate the *process* model, whether the
deviations behave like their assumed iid normal process, a complement to
the observation-based
[`osa_residuals()`](https://afsc-assessments.github.io/Rceattle/reference/osa_residuals.md).

Each set of deviations holds a Gaussian process prior in the model. The
posterior *mode* of the deviations is shrunk toward that prior, so,
following SAM, a single draw is taken from the joint posterior of the
deviations (from the joint precision when they are random effects, or
the fixed-effect covariance when they are penalized fixed effects) and
standardized by the process standard deviation. Under a correctly
specified process these are approximately iid N(0, 1). Because a random
draw is used the residuals are stochastic; set `seed` for
reproducibility.

Supported processes and their deviations:

- `"recruitment"`:

  `rec_dev`, prior `N(-bias_adjust_proc * R_sd^2/2, R_sd^2)` per
  species.

- `"initial"`:

  `init_dev`, prior `N(-bias_adjust_proc * R_sd^2/2, R_sd^2)` per
  species.

- `"catchability"`:

  `index_q_dev`, prior `N(0, q_dev_sd^2)` per index.

`"all"` returns every supported process present in the fit. Selectivity
and natural-mortality deviations (which use random-walk / 2D-AR1 priors)
are not yet supported. Recruitment and initial-abundance residuals are
exact (their priors are iid). Catchability residuals are exact only for
the iid deviate prior (`Time_varying_q = 1` or `2`); for a random-walk
or AR1 catchability prior the marginal-SD standardization ignores the
prior correlation, so those residuals are approximate and a warning is
emitted.

## Usage

``` r
process_residuals(
  object = NULL,
  process = c("recruitment", "initial", "catchability", "all"),
  seed = 123,
  fit = NULL
)
```

## Arguments

- object:

  A fitted `Rceattle` model. The targeted deviations must be estimated,
  as random effects (e.g. `random_rec = TRUE`) or as penalized fixed
  effects, with a usable covariance.

- process:

  One of `"recruitment"`, `"initial"`, `"catchability"`, or `"all"`.

- seed:

  Seed for the posterior draw. Default 123.

- fit:

  deprecated name for `object`, still accepted so existing scripts keep
  working. Supplying both is an error.

## Value

A data frame of class `rceattle_osa` (so it can be passed to
[`osa_diagnostics()`](https://afsc-assessments.github.io/Rceattle/reference/osa_diagnostics.md)
and
[`plot.rceattle_osa()`](https://afsc-assessments.github.io/Rceattle/reference/plot.rceattle_osa.md))
with one row per deviation: columns `source` (the process), `fleet`,
`species`, `year`, `age_length_bin`, and `residual`.

## References

Nielsen, A., and Berg, C.W. 2014. Estimation of time-varying selectivity
in stock assessments using state-space models. Fish. Res. 158:96-101.

## See also

[`osa_residuals()`](https://afsc-assessments.github.io/Rceattle/reference/osa_residuals.md),
[`osa_diagnostics()`](https://afsc-assessments.github.io/Rceattle/reference/osa_diagnostics.md)
