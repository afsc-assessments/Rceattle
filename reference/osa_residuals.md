# One-step-ahead (OSA) residuals for an Rceattle model

Computes one-step-ahead (OSA) residuals, also called forecast or
quantile residuals (Thygesen et al. 2017), for a fitted
[Rceattle](https://afsc-assessments.github.io/Rceattle/reference/Rceattle-package.md)
model via
[`TMB::oneStepPredict()`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html).
Unlike Pearson residuals, OSA residuals are distributed iid standard
normal under a correctly specified model even when observations are
correlated (through composition bins) or when the model contains random
effects, so they support objective goodness-of-fit testing (Trijoulet et
al. 2023; Stewart and Monnahan 2025).

These are *internal* OSA residuals: the residualization is integrated
into the assessment via TMB, so it also accounts for correlation induced
by the model's random effects across years, the gold standard relative
to the *external* `compResidual` approach (Stewart and Monnahan 2025).

OSA residuals are computed *post hoc* and are expensive (TMB
re-optimizes the random effects as each observation is added), so they
are not produced during
[`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md).
The model must have been optimized with `estimateMode < 3`.

Supported observation types are the aggregate `"catch"` and `"index"`
series, the `"comp"` (age/length composition) and `"caal"` (conditional
age-at-length) compositions, and `"diet"` (predator stomach-content
composition, for multispecies models with estimated suitability). Diet
is opt-in (not in the default `types`) because it applies only to
multispecies models and can be expensive.

For composition data the multivariate multinomial /
Dirichlet-multinomial is decomposed into a sequence of univariate
conditional residuals (binomial / beta-binomial; Trijoulet et al. 2023).
The final bin of each composition is fixed by the sum-to-N constraint
and so has no residual (returned as `NA`). Composition OSA uses an
internal model rebuild with unweighted, proper densities (the `osa_mode`
switch); fleets fit with the AFSC `MultinomialAFSC` pseudo-likelihood
are residualized under the full multinomial.

Survey-index OSA residuals are supported for every index likelihood
family (`Index_distribution`). Lognormal IID (`"Lognormal"`)
residualizes on the log scale, and the natural-scale `"Normal"` and
`"TruncatedNormal"` on the natural scale. The correlated covariance
families (`"MVN"` / `"MVNORM"`) are whitened by the lower Cholesky of
the fleet's survey covariance Sigma = L L', so the residuals are the
multivariate-Gaussian one-step-ahead innovations L^-1 (obs - q\*pred),
the closed form
[`TMB::oneStepPredict()`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html)
reproduces for a Gaussian block.

Under a Gaussian `method`, `"TruncatedNormal"` rows are residualized in
their own
[`TMB::oneStepPredict()`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html)
call, with `method = "oneStepGeneric"` and a range starting at zero.
(Under `method = "cdf"` none of the rest of this applies: the model
supplies the truncated CDF in closed form, so the family is residualized
in the main call, exactly.) Its density differs from `"Normal"` only by
`log Phi(mu/sd)`, which is a function of the prediction and not of the
observation, so a Gaussian method, which reads the curvature of the
density in the observation, cannot see the truncation at all and returns
the untruncated `(obs - mu)/sd`. Integrating over the family's own
support instead gives the truncated CDF
`F(x) = [Phi((x - mu)/sd) - Phi(-mu/sd)] / Phi(mu/sd)`, so `qnorm(F(x))`
is standard normal by the probability integral transform however hard
the truncation bites. The upper limit is finite rather than `Inf`, ten
standard deviations past the largest fitted index in the group, which
leaves under 1e-23 of the mass outside while keeping the Laplace inner
problem in a region it can solve. That group also runs with
`splineApprox = FALSE`, because the spline shortcut integrates over
whatever range its profile slice covered. The range is a property of the
family, so it cannot be shared with the other fleets: `"Normal"` is
genuinely untruncated, and a lognormal fleet's stored observation is
`log(obs)`, which is negative for a small index.

The size of the correction is the truncated mass: on a fleet predicting
100 with an absolute sd of 150 (a quarter of the density below zero), an
observation exactly at the prediction has an untruncated residual of 0
and a truncated one of -0.44.

Three consequences of that group being residualized separately, worth
knowing before reading the output:

- **On a model with random effects the exact integration can fail.** It
  evaluates the Laplace marginal at arbitrary values of the observation,
  and the inner problem does not always converge across the whole
  support. Rather than return `NA` for those rows, which would shrink
  the sample
  [`osa_diagnostics()`](https://afsc-assessments.github.io/Rceattle/reference/osa_diagnostics.md)
  passes verdict on, without saying so, the fleet is recomputed under
  TMB's spline approximation and a warning says the residuals for it are
  approximate. Fixed-effect models are unaffected.

- **`sd` is `NA` and `predicted` means something different for this
  group.** `oneStepGeneric` returns neither a standard deviation nor the
  fitted value; `predicted` is the truncated conditional mean
  `E[x | x > 0]`, which sits above the fitted index (163.8 against a
  fitted 100.0 on the fixture above), not the fitted index itself. The
  residual is unaffected.

- **Adding a `"TruncatedNormal"` fleet can still move the other fleets'
  residuals on a random-effects model.** Each
  [`TMB::oneStepPredict()`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html)
  call is given everything earlier in the sequence as `conditional`, so
  a group that is a contiguous block reproduces a single call exactly
  (measured at 1.5e-14 on a 21-random-effect fixture; without it the
  same split moved a residual by 5.8e-2). A truncated fleet's rows are
  interleaved with the other index fleets by year rather than
  contiguous, so the rows falling inside its span are still marked
  unconditional, which zeroes their data terms and changes the
  conditional distribution of the latent states. Both orderings are
  valid probability-integral-transform sequences, so no value is wrong,
  but catch and lognormal-index residuals need not match those from an
  otherwise identical model with no truncated fleet. Fixed-effect models
  are unaffected, since their observations are independent.

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

  A fitted object of class `Rceattle` (from
  [`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)).

- source:

  Character vector of observation sources to residualize: any of
  `"ecov"`, `"index"`, `"catch"`, `"comp"`, `"caal"`, `"diet"`, or
  `"all"`. Defaults to the five non-diet sources (`diet` is opt-in
  because it applies only to multispecies models and can be expensive);
  pass `"all"` to include `diet`. `"ecov"` is the state-space covariate
  (QAR1 `observe=` term), residualized first against its own series, as
  in WHAM's `make_osa_residuals()`. Sources with no observations in the
  model are silently skipped. Mirrors the `source` argument of
  [`residuals.Rceattle()`](https://afsc-assessments.github.io/Rceattle/reference/residuals.Rceattle.md)
  and
  [`plot.rceattle_osa()`](https://afsc-assessments.github.io/Rceattle/reference/plot.rceattle_osa.md).

- method:

  One of `"oneStepGaussianOffMode"` (default; the WHAM/SAM choice),
  `"oneStepGaussian"`, `"fullGaussian"`, `"oneStepGeneric"` or `"cdf"`,
  passed to
  [`TMB::oneStepPredict()`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html).
  See the section below on choosing between the Gaussian methods and
  `"cdf"`.

- discrete:

  Logical; whether to treat *composition* (comp / caal / diet)
  observations as the discrete counts they are. `NULL` (the default)
  chooses per method: `TRUE` under `method = "cdf"`, where it is needed
  for the residual to be standard normal at all, and `FALSE` otherwise,
  matching how CEATTLE fits the composition likelihood with
  effective-sample-size-scaled counts. Passing `FALSE` under `"cdf"` is
  allowed, and says in a message that those composition residuals are
  biased up. When `TRUE`, composition residuals are randomized quantile
  residuals (Dunn and Smyth 1996) and so are stochastic; set `seed` for
  reproducibility. The aggregate index/catch series are always
  continuous (lognormal); the
  [`TMB::oneStepPredict()`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html)
  call is split by observation type so `discrete` is applied correctly
  per type. A Gaussian `method` cannot score a discrete observation, so
  that group falls back to `"oneStepGeneric"`; `method = "cdf"` already
  reads a CDF and keeps it.

- parallel:

  Logical; compute the per-observation OSA loop in parallel via
  [`mclapply`](https://rdrr.io/r/parallel/mclapply.html). Default
  `TRUE`. This is the main speedup for models with random effects, where
  each observation triggers a Laplace re-evaluation, it gives a
  near-linear speedup across cores (set `options(mc.cores = )` to choose
  how many; forking falls back to serial on Windows). Some models, heavy
  random-effect structures such as a time-varying catchability, abort
  the forked worker instead of returning; the loop then recomputes
  serially, after rebuilding, and prints the worker's own "irrecoverable
  exception" message, which comes from C and cannot be suppressed. That
  message does not mean the call failed. Pass `FALSE` to skip the
  attempt. `method = "cdf"` is parallelized whether or not `discrete` is
  `TRUE`: the workers only evaluate the objective, and
  [`TMB::oneStepPredict()`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html)
  draws the randomizing uniforms under `seed` once they return, so the
  residuals are the serial ones. A discrete group under
  `"oneStepGeneric"` runs serially, where that has not been measured.

- seed:

  Random seed passed to
  [`TMB::oneStepPredict()`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html)
  for reproducibility of randomized-quantile residuals. Default `123`.

- trace:

  Logical; print
  [`TMB::oneStepPredict()`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html)
  progress. Default `FALSE`.

- ...:

  Further arguments passed to
  [`TMB::oneStepPredict()`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html).

- fit:

  deprecated name for `object`, still accepted so existing scripts keep
  working. Supplying both is an error.

## Value

A data frame of class `rceattle_osa` with one row per residualized
observation and columns `source` (the data source:
index/catch/comp/caal/ diet), `fleet`, `fleet_name`, `species`, `sex`,
`year`, `age_length_bin` (the age or length bin the value stands for),
`accumulated` (`TRUE` where tail accumulation folded neighbouring bins
into this one, so it covers a range of ages rather than the one named),
`length` (the conditioning length bin for caal; `NA` otherwise),
`index_label` (`"age"`/`"length"`/`NA`), `observed`, `predicted`, `sd`,
and `residual`. For aggregate series `observed` and `predicted` are on
the residualization scale, log for lognormal catch/index, natural scale
for a `"Normal"` or `"TruncatedNormal"` index, and the whitened (`L^-1`)
scale for an `"MVN"`/`"MVNORM"` index; for compositions they are bin
counts. `predicted` is `NA` for every row under `method = "cdf"`, which
forms no conditional mode (see above). `sd` is `NA` under the default
method too, only `method = "oneStepGaussian"` returns one. It holds
`method` and `seed` attributes, `method` is the string that was passed,
or a named vector `c(default = <method>, ...)` when a group of rows was
residualized with its own. Two names are likelihood families:
`TruncatedNormal = "oneStepGeneric"` for a truncated index fleet under a
Gaussian method, and `DirichletMultinomial = "oneStepGaussianOffMode"`
for a D-M composition under `"cdf"`. The third is not a family but the
rows `discrete = TRUE` selects:
`DiscreteComposition = "oneStepGeneric"`, written only under a Gaussian
`method`, which cannot score a discrete observation because it is
continuous-only (`"cdf"` already is a CDF method, so it needs no
override and writes no entry). It also holds (when composition types are
present) a `"pearson"` attribute holding the matching Pearson residuals
so
[`plot.rceattle_osa()`](https://afsc-assessments.github.io/Rceattle/reference/plot.rceattle_osa.md)
can show both. The attribute uses this data frame's column names rather
than the data-sheet names
[`residuals.Rceattle()`](https://afsc-assessments.github.io/Rceattle/reference/residuals.Rceattle.md)
returns, so the two halves of one object read alike. Note the shared
names do not mean a shared scale: in the attribute `observed` and
`predicted` are proportions summing to one within a fleet-year, with the
sample size in `sample_size`, because composition Pearson residuals are
defined on proportions, not the bin counts the columns above hold. Do
not compare the two directly. Both describe the bins the likelihood fit,
so a fleet with tail accumulation reports the folded window in each,
with one asymmetry: the one-step-ahead decomposition drops each group's
last bin (it is fixed by sum-to-N), and under an *old*-tail accumulation
that dropped bin is the upper accumulated one. Such a fleet therefore
shows its upper boundary bin in the Pearson residuals and not in the OSA
residuals. Summarize it with
[`osa_diagnostics()`](https://afsc-assessments.github.io/Rceattle/reference/osa_diagnostics.md)
and plot it with
[`plot.rceattle_osa()`](https://afsc-assessments.github.io/Rceattle/reference/plot.rceattle_osa.md).

## Choosing a method

The default `"oneStepGaussianOffMode"` approximates each conditional
distribution as Gaussian: it treats the observation as a free variable,
finds the mode of the conditional density in that coordinate, and
standardises the observation against it. That is fast, and it is what
WHAM and SAM use.

`"cdf"` instead asks the model for the conditional CDF and returns
`qnorm(F(x))`, the probability integral transform, which is standard
normal whatever shape the conditional has. No conditional mean is
formed, so the method has three properties the Gaussian ones do not:

- **`predicted` is `NA` on every row** (and so is `sd`, which
  `"oneStepGaussianOffMode"` does not return either).
  [`TMB::oneStepPredict()`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html)
  gives `Fx`, `px` and `nll` for this method and no fitted value, and
  there is no conditional mode to put in the column; a marginal fitted
  value in its place would mean the conditional mode under one method
  and the marginal fit under another. Fitted values are in
  `fit$quantities` and in `residuals(fit, type = "pearson")`. The column
  is blanked on the Dirichlet-multinomial rows too, which run under a
  Gaussian method, so it means one thing across the object.

- **The conditional mean cannot leave the support**, because it is never
  computed. That removes the negative composition `predicted` values
  described below, and the positive bias they pass into the residual on
  those rows.

- **`Index_distribution = "TruncatedNormal"` is exact and needs no
  separate call**, so none of the three consequences listed above
  applies to that family under `"cdf"`.

What it costs, in three places.

- **Not available for a Dirichlet-multinomial composition**: the
  conditional is a beta-binomial, which has no closed-form CDF and
  cannot be summed at a fractional count. Those fleets are residualized
  with `"oneStepGaussianOffMode"`, announced in a message and recorded
  in the `method` attribute.

- **`|residual|` is censored at 8.04, in both directions.** The upper
  end is forced:
  [`TMB::oneStepPredict()`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html)
  recovers `F` as `1 / (1 + exp(.))`, which saturates at the last double
  below one. The lower end is not, that same expression takes a small
  `F` down to a residual of -37, and is censored to match anyway,
  because an asymmetric ceiling would show as a long left tail against a
  wall on the right, which is what skewness in the residuals looks like.
  This is a ceiling, not a large number standing in for a larger one:
  [`osa_diagnostics()`](https://afsc-assessments.github.io/Rceattle/reference/osa_diagnostics.md)
  computes SDNR and the tail statistics on the censored values, so it
  bites hardest on a short series where one observation drives the
  statistic, and the function warns when any residual sits there.

  For an observation past the ceiling, **reach for `"oneStepGaussian"`
  specifically** and on the fleet in question, on a 12-year survey with
  one observation multiplied by 200 it reports 38.98 uncensored, where
  the package default returns `NaN` and `"oneStepGeneric"` compresses
  the same row to 3.33. It costs an `nlminb` and an `optimHess` per
  observation. The comparison is in
  [`vignette("model-diagnostics")`](https://afsc-assessments.github.io/Rceattle/articles/model-diagnostics.md)
  and reproduced by `tools/verify/verify-osa-cdf-accuracy.R`.

- **Compositions are residualized in their own
  [`TMB::oneStepPredict()`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html)
  call**, because `discrete` differs between them and the aggregate
  series and TMB takes one setting per call. Everything earlier in the
  sequence is passed as `conditional`, so the CONDITIONING is what a
  single call would have used. The residual values are not identical to
  a single call under `discrete = TRUE`: oneStepPredict re-seeds and
  draws `nrow(subset)` randomizing uniforms per call, so a composition
  residual depends on how many rows share its call. Set `seed` for
  reproducibility of a given call, and do not compare individual
  randomized residuals across different `source` selections.
  Fixed-effect models are unaffected by the conditioning.

**Which is right for compositions.** Residualizing simulated data at the
parameters that generated it makes the answer exactly standard normal,
so the methods can be scored rather than argued about
(`tools/verify/verify-osa-cdf.R`, BS2017SS, 20 replicates, 4538
composition bins each):

|  |  |  |  |  |
|----|----|----|----|----|
| method | mean | sd | lag-1 acf within a composition | KS rejects |
| `oneStepGaussianOffMode` | +0.103 | 0.918 | +0.060 | every replicate |
| `cdf`, `discrete = FALSE` | +0.610 | 1.262 | +0.434 | every replicate |
| `cdf`, `discrete = TRUE` | -0.007 | 0.995 | +0.001 | none |

The null standard error on the autocorrelation is 0.015. A composition
bin holds a count, so its conditional CDF is a step function and
`qnorm(F(x))` inherits the step; only the randomized quantile residual
`qnorm(F(x) - U f(x))` (Dunn and Smyth 1996, the construction Trijoulet
et al. 2023 prescribe) removes it. That is why `discrete` defaults to
`TRUE` under this method, and why `cdf` with `discrete = FALSE` is the
worst of the three. On the aggregate index and catch series, which are
genuinely Gaussian, all three agree (to 5e-5 wherever `|residual| < 8`,
i.e. everywhere the `"cdf"` ceiling above does not bind) and all three
pass.

## Random effects, where the choice reverses for the aggregate series

`"cdf"` is exact only on a fixed-effect model. With random effects
[`TMB::oneStepPredict()`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html)
integrates the CDF over the latent states by Laplace, and the integrand
there is a Gaussian times a sigmoid rather than a density, so the
approximation is not exact, where for a Gaussian observation the
Gaussian methods integrate a density and are. Against the exact Kalman
innovations of a linear-Gaussian state space model the Gaussian methods
are exact to 1e-14, while `"cdf"` errs by 7e-4 to 4e-2 as the latent
state becomes more informative relative to the observation
(`tools/verify/verify-osa-cdf-accuracy.R`, which compiles that model).

**That result is about a LINEAR-Gaussian model, and does not transfer
wholesale.** Those methods are exact when the one-step-ahead
*predictive* is Gaussian, which needs the model to be linear in the
random effects. Rceattle's index and catch are
[`exp()`](https://rdrr.io/r/base/Log.html) of cumulated log recruitment
deviations pushed through the population dynamics, and are not.
`fullGaussian` and `oneStepGaussian` cannot differ for a Gaussian
conditional, so their disagreement measures the departure, on a
17-deviation fixture they differ by 0.091 on both index and catch, while
`"cdf"` differs from `oneStepGaussian` by 0.017 on index. **No method is
exact for index or catch under random effects**, and the choice there is
not settled by this package.

**`"ecov"` is the exception, and there a Gaussian method is right.** Its
conditional, a Gaussian measurement of an AR1 latent, every other data
term unconditional, genuinely is linear-Gaussian: `fullGaussian` and
`oneStepGaussian` agree to 4e-14 on the QAR1 fixture in
`test-likelihood-osa-cdf.R`, while `"cdf"` sits 0.139 away, about a
quarter of the residual standard deviation.

**It does not reverse for compositions.** Their conditional is a
discrete, skewed binomial, which is what the Gaussian methods get wrong,
and by far more than the Laplace error above. Simulating from a
22-random-effect model with the recruitment deviations redrawn and
residualizing at the generating parameters (1680 residuals, 120
replicates; null standard errors 0.024 and 0.017),
`"oneStepGaussianOffMode"` gives mean +0.513 and sd 0.404 with
Kolmogorov-Smirnov rejecting all 120, against mean +0.006, sd 1.002 and
6 of 120, the nominal 5%, for `"cdf"` with `discrete = TRUE`. The
Gaussian default is not merely biased there; it is under-dispersed by a
factor of two and a half. What limits `"cdf"` on compositions is scale,
not random effects, see the section below.

## Known limitation, compositions at scale under random effects

On a random-effects model with a large composition data set,
`method = "cdf"` returns non-finite residuals in bulk and is very slow.
Measured on `BS2017SS` with `random_rec = TRUE` (159 random effects,
4538 composition bins): **1879 of 4538 residuals non-finite**, against 0
for `"oneStepGaussianOffMode"` on the same fit, and hours rather than
minutes.

The failures are a contiguous tail, and the same 1880 rows residualized
on their own return 1 failure, so the observations are not the problem.
What separates the two runs is the depth of the conditioning (2658 prior
observations against none): the Laplace inner problem fails on the
conditioning itself, and redoing the tail on a fresh call does not
recover it (1879 before, 1879 after).

What binds is the depth, not the presence of random effects: the same
method residualizes 1680 composition bins on a 22-random-effect model
correctly (the section above). So try `"cdf"` and read the warning it
issues, when it returns non-finite residuals in bulk, fall back to a
Gaussian `method` for that source, remembering that its composition
residuals are under-dispersed by about a factor of two and a half.
`"cdf"` is sound on fixed-effect models, and on random-effects models
for the aggregate and covariate series, which are few enough not to
reach the depth where this bites.

One caveat that applies to every method, not just this one:
`Estimate_index_sd = "Analytical"` / `Estimate_catch_sd = "Analytical"`
and the analytical `Catchability` forms concentrate their parameter out
of the likelihood using **all** the observations, including the one
being residualized. The conditioning is then not strictly one-step-ahead
and the residuals are approximate, by an amount that shrinks as the
series lengthens.

## Negative composition `predicted` values

This section describes the Gaussian methods. Under `method = "cdf"` no
conditional mean is formed and `predicted` is `NA`, so none of it
applies.

A composition `predicted` is an expected bin count and cannot truly be
negative, but it goes slightly negative where a bin holds almost no
fish. Composition observations enter as counts,
`(proportion + comp_offset) * N`, and
[`TMB::oneStepPredict()`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html)'s
conditional mean is a numerical step away from the observation, which
overshoots below zero when the count is near it. The function warns,
naming the count and the years.

It is the *bin's* count that drives this, not the composition's sample
size: on EBS pollock the negative rows have a median observed count of
0.05 against 4.9 for the rest, while their sample sizes span the same 1
to 821 as everything else (69 of 353 occur above a sample size of 100).
A rare age in a well-sampled year does it as readily as a poorly sampled
year, so the warning is not by itself evidence of thin data.

The values are reported rather than clamped, because a negative expected
count is the signal that the bin is too sparse for the decomposition to
describe and clamping would hide it. Treat `predicted` on those rows as
uninformative; their `residual` standardises the observation against
that same mean and is therefore biased positive.

## References

Thygesen, U.H., et al. 2017. Validation of ecological state space models
using the Laplace approximation. Environ. Ecol. Stat. 24:317-339.

Trijoulet, V., et al. 2023. Model validation for compositional data in
stock assessment models. Fish. Res. 257:106487.

Stewart, I.J., and Monnahan, C.C. 2025. Diagnosing common sources of
lack of fit to composition data using one-step-ahead residuals. Can. J.
Fish. Aquat. Sci. 82:1-13.

## See also

[`osa_diagnostics()`](https://afsc-assessments.github.io/Rceattle/reference/osa_diagnostics.md),
[`plot.rceattle_osa()`](https://afsc-assessments.github.io/Rceattle/reference/plot.rceattle_osa.md),
[`process_residuals()`](https://afsc-assessments.github.io/Rceattle/reference/process_residuals.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(BS2017SS)
fit <- fit_mod(BS2017SS, estimateMode = "Hindcast")
osa <- osa_residuals(fit, source = c("index", "comp"))
plot(osa)
} # }
```
