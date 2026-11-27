# One-step-ahead (OSA) residuals for an Rceattle model

One-step-ahead residuals, also called forecast or quantile residuals
(Thygesen et al. 2017), computed via
[`TMB::oneStepPredict()`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html).
Unlike Pearson residuals they are iid standard normal under a correctly
specified model even when observations are correlated through
composition bins or the model holds random effects, so they support
objective goodness-of-fit testing (Trijoulet et al. 2023; Stewart and
Monnahan 2025). The residualization runs inside the assessment, so it
also accounts for correlation induced by the model's own random effects.

## Usage

``` r
osa_residuals(
  object = NULL,
  source = c("ecov", "index", "catch", "comp", "caal"),
  method = "oneStepGaussianOffMode",
  discrete = NULL,
  parallel = TRUE,
  seed = 123,
  trace = FALSE,
  ...,
  fit = NULL
)
```

## Arguments

- object:

  A fitted `Rceattle` object, from
  [`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md).

- source:

  Observation sources to residualize: any of `"ecov"`, `"index"`,
  `"catch"`, `"comp"`, `"caal"`, `"diet"` or `"all"`, defaulting to the
  five non-diet sources. Mirrors the `source` argument of
  [`residuals.Rceattle()`](https://afsc-assessments.github.io/Rceattle/reference/residuals.Rceattle.md);
  a source with no observations is skipped silently, `"diet"` is opt-in
  because it needs a multispecies model and is expensive, and `"ecov"`
  is the state-space covariate, residualized first against its own
  series as in WHAM's `make_osa_residuals()`.

- method:

  How the conditional distribution is read: one of
  `"oneStepGaussianOffMode"` (default, the WHAM/SAM choice),
  `"oneStepGaussian"`, `"fullGaussian"`, `"oneStepGeneric"` or `"cdf"`,
  passed to
  [`TMB::oneStepPredict()`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html).

- discrete:

  Whether to treat composition observations as the discrete counts they
  are. `NULL` (default) picks `TRUE` under `method = "cdf"`, where the
  residual is not standard normal without it, and `FALSE` otherwise.
  `TRUE` makes composition residuals randomized quantile residuals (Dunn
  and Smyth 1996) and so stochastic – set `seed`. A Gaussian `method`
  cannot score a discrete observation, so that group falls back to
  `"oneStepGeneric"`; `"cdf"` with `discrete = FALSE` is allowed and
  says in a message that those residuals are biased up.

- parallel:

  Compute the per-observation loop with
  [`mclapply`](https://rdrr.io/r/parallel/mclapply.html), default
  `TRUE`, which is the main speedup for a model with random effects; set
  `options(mc.cores = )` to choose how many. Where a forked worker
  aborts instead of returning, the loop recomputes serially and prints
  the worker's own "irrecoverable exception" message, which comes from C
  and cannot be suppressed – it does not mean the call failed.

- seed:

  Seed passed to
  [`TMB::oneStepPredict()`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html)
  for reproducible randomized-quantile residuals. Default `123`.

- trace:

  Print
  [`TMB::oneStepPredict()`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html)
  progress. Default `FALSE`.

- ...:

  Further arguments passed to
  [`TMB::oneStepPredict()`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html).

- fit:

  Deprecated name for `object`, still accepted; supplying both is an
  error.

## Value

A data frame of class `rceattle_osa`, one row per residualized
observation, with columns `source`, `fleet`, `fleet_name`, `species`,
`sex`, `year`, `age_length_bin`, `accumulated` (the bin folds
neighbouring ages, so it stands for a range rather than the age named),
`length` (the conditioning length bin for caal, `NA` otherwise),
`index_label`, `observed`, `predicted`, `sd` and `residual`.

`observed` and `predicted` are on the residualization scale: log for
lognormal catch and index, natural scale for `"Normal"` and
`"TruncatedNormal"`, whitened by the lower Cholesky of the survey
covariance for `"MVN"`/`"MVNORM"`, and bin counts for compositions.
Three cases read differently, and none is an error:

- Under `method = "cdf"`, `predicted` is `NA` on every row because the
  method forms no conditional mode. `sd` is `NA` under the default
  method too; only `"oneStepGaussian"` returns one.

- On a `"TruncatedNormal"` fleet under a Gaussian method, `predicted` is
  the truncated conditional mean `E[x | x > 0]`, which sits above the
  fitted index (163.8 against a fitted 100.0 on the package fixture),
  and `sd` is `NA`. The residual is unaffected.

- A composition `predicted` can go slightly negative where a bin holds
  almost no fish, and the function warns, naming the count and the
  years. The bin's own count drives it, not the composition's sample
  size – median observed count 0.05 on those rows against 4.9 elsewhere
  – so the warning is not by itself evidence of thin data. Those
  residuals are biased positive.

Adding a `"TruncatedNormal"` fleet can move the *other* fleets'
residuals on a random-effects model, because its rows interleave with
the other index fleets by year instead of forming a contiguous block: a
contiguous group reproduces a single call to 1.5e-14, while an
interleaved split moved a residual by 5.8e-2. Both orderings are valid
probability-integral-transform sequences, so neither is wrong.
Fixed-effect models are unaffected.

Attributes carry `seed` and `method` – the string that was passed, or a
named vector where a group of rows ran under its own, such as
`TruncatedNormal = "oneStepGeneric"`. Where composition types are
present a `"pearson"` attribute holds the matching Pearson residuals, on
proportions with the sample size in `sample_size` rather than the bin
counts above, so the two are not comparable column by column.

## Details

Computed *post hoc* and expensive – TMB re-optimizes the random effects
as each observation is added – so
[`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)
does not produce them, and the fit must have been optimized at
`estimateMode < 3`. Composition data are decomposed into univariate
conditional residuals (binomial / beta-binomial; Trijoulet et al. 2023),
using code adapted from WHAM (Stock and Miller 2021). The last bin of
each composition is fixed by the sum-to-N constraint, so it has no
residual (`NA`).

Choosing between the methods, what each costs, how they behave under
random effects, the ceiling on `|residual|` under `"cdf"`, and where
compositions at scale defeat it are in
[`vignette("model-diagnostics")`](https://afsc-assessments.github.io/Rceattle/articles/model-diagnostics.md),
with the measurements behind each.

## References

Thygesen, U.H., et al. 2017. Validation of ecological state space models
using the Laplace approximation. Environ. Ecol. Stat. 24:317-339.

Trijoulet, V., et al. 2023. Model validation for compositional data in
stock assessment models. Fish. Res. 257:106487.

Stock, B.C., and Miller, T.J. 2021. The Woods Hole Assessment Model
(WHAM): a general state-space assessment framework that incorporates
time- and age-varying processes via random effects and links to
environmental covariates. Fish. Res. 240:105967.

Stewart, I.J., and Monnahan, C.C. 2025. Diagnosing common sources of
lack of fit to composition data in fisheries stock assessment models
using one-step-ahead (OSA) residuals. Can. J. Fish. Aquat. Sci. 82:1-13.

Dunn, P.K., and Smyth, G.K. 1996. Randomized quantile residuals. J.
Comput. Graph. Stat. 5:236-244.

## See also

[`osa_diagnostics()`](https://afsc-assessments.github.io/Rceattle/reference/osa_diagnostics.md),
[`plot.rceattle_osa()`](https://afsc-assessments.github.io/Rceattle/reference/plot.rceattle_osa.md),
[`process_residuals()`](https://afsc-assessments.github.io/Rceattle/reference/process_residuals.md),
[`vignette("model-diagnostics")`](https://afsc-assessments.github.io/Rceattle/articles/model-diagnostics.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(BS2017SS)
fit <- fit_mod(BS2017SS, estimateMode = "Hindcast")
osa <- osa_residuals(fit, source = c("index", "comp"))
plot(osa)
} # }
```
