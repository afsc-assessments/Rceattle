# Fit the CEATTLE assessment model

Estimate CEATTLE population parameters by maximum likelihood, and
optionally project the stock and apply a harvest control rule.

## Usage

``` r
fit_mod(
  data_list = NULL,
  inits = NULL,
  map = NULL,
  bounds = NULL,
  file = NULL,
  estimateMode = 0,
  random_rec = FALSE,
  random_q = FALSE,
  random_sel = FALSE,
  HCR = build_hcr(),
  niter = 3,
  recFun = build_srr(),
  M1Fun = build_M1(),
  growthFun = build_growth(),
  qFun = build_catchability(),
  selFun = build_selectivity(),
  compFun = build_composition(),
  msmMode = 0,
  avgnMode = 0,
  initMode = "NonEquilibrium",
  suitMode = 0,
  suit_styr = NULL,
  suit_endyr = NULL,
  fit_control = NULL,
  config = NULL,
  quiet_data_check = FALSE,
  ...
)
```

## Arguments

- data_list:

  A data list read in via
  [`read_data`](https://afsc-assessments.github.io/Rceattle/reference/read_data.md)
  or built directly in R; see
  [`vignette("data-without-excel", package = "Rceattle")`](https://afsc-assessments.github.io/Rceattle/articles/data-without-excel.md).

- inits:

  (Optional) A named list of initial parameter values, as returned by
  [`build_params`](https://afsc-assessments.github.io/Rceattle/reference/build_params.md)
  or extracted from a previous fit (`model$estimated_params`). If
  `NULL`, parameters are initialized from scratch via
  [`build_params`](https://afsc-assessments.github.io/Rceattle/reference/build_params.md).

- map:

  (Optional) A map object from
  [`build_map`](https://afsc-assessments.github.io/Rceattle/reference/build_map.md).

- bounds:

  (Optional) A bounds object from
  [`build_bounds`](https://afsc-assessments.github.io/Rceattle/reference/build_bounds.md).

- file:

  (Optional) Filename where files will be saved. If NULL, no file is
  saved.

- estimateMode:

  What to fit, given as a string alias or the integer code:
  `"Estimate"` (0) = fit the hindcast model and the HCR projection
  (`HCR`); `"Hindcast"` (1) = fit the hindcast only (no fitting
  BRPs/HCR/projection); `"Projection"` (2) = fit the BRPs/HCR/projection
  only, from the initial parameters in `inits`; `"DebugBuild"` (3) =
  build through `MakeADFun` but not `nlminb`, the returned `obj` holds
  the real objective and gradient, so `obj$fn()` / `obj$gr()` are usable
  for diagnosing a model before committing to a fit;
  `"DebugOptimize"` (4) = optimize with all parameters mapped out, so
  the objective is a placeholder (`dummy^2`), not a likelihood. Defaults
  to `"Estimate"`.

- random_rec:

  logical. If TRUE, treats recruitment deviations as random effects
  using the Laplace approximation. The default is FALSE.

- random_q:

  logical (default FALSE); if TRUE the `Time_varying_q` deviations are
  integrated as random effects with one estimated sd per
  `Catchability_index` group, not fixed at `Time_varying_q_sd` (linkage
  random effects are integrated either way).

- random_sel:

  logical (default FALSE); if TRUE the `Time_varying_sel` deviations are
  integrated as random effects with one estimated sd per
  `Selectivity_index` group, not fixed at `Time_varying_sel_sd` (linkage
  random effects are integrated either way).

- HCR:

  HCR list object from
  [`build_hcr`](https://afsc-assessments.github.io/Rceattle/reference/build_hcr.md)

- niter:

  Number of iterations for multispecies model

- recFun:

  The stock recruit-relationship parameterization from
  [`build_srr`](https://afsc-assessments.github.io/Rceattle/reference/build_srr.md).

- M1Fun:

  M1 parameterizations and priors. Use `build_M1`.

- growthFun:

  The weight-at-age parameterization from
  [`build_growth`](https://afsc-assessments.github.io/Rceattle/reference/build_growth.md).

- qFun:

  Catchability specification from
  [`build_catchability`](https://afsc-assessments.github.io/Rceattle/reference/build_catchability.md),
  holding any environmental linkages on q.

- selFun:

  Selectivity specification from
  [`build_selectivity`](https://afsc-assessments.github.io/Rceattle/reference/build_selectivity.md),
  holding any environmental linkages on selectivity parameters.

- compFun:

  Composition-weighting specification from
  [`build_composition`](https://afsc-assessments.github.io/Rceattle/reference/build_composition.md),
  holding any priors on the Dirichlet-multinomial weights.

- msmMode:

  The predation-mortality mode, as a string alias or integer code:
  `"SingleSpecies"` (0, the default, no predation), `"MSVPA"` (1, the
  Type-II MSVPA predation of Holsman et al. 2015) or `"TypeIIIMSVPA"`
  (2). Higher integer codes (Kinzey-Punt, Holling forms) are declared
  but not implemented, and are rejected by the data check.

- avgnMode:

  the average abundance-at-age approximation used in the
  predation-mortality equations. Only mode 0, \\N/Z(1 - exp(-Z))\\ (the
  MSVPA form), is implemented; the alternatives are declared but have no
  effect.

- initMode:

  how the population is initialized, as a string alias or integer code:
  `"FreeParams"` (0), `"Equilibrium"` (1), `"NonEquilibrium"` (2, the
  default), `"FishedNonEquilibrium"` (3), `"FishedNonEquilibriumScaled"`
  (4), `"OffsetEquilibrium"` (5). See the **Initial age structure**
  section below for what each one estimates.

- suitMode:

  how predator-prey suitability is derived, per predator (a single value
  or a vector of length `nspp`): 0 = empirical from diet data (Holsman
  et al. 2015), 2 = weight-based gamma, 4 = weight-based lognormal, 6 =
  weight-based normal. The length-based forms (1, 3, 5) are declared but
  not implemented and are rejected by the data check.

- suit_styr:

  The first year used to calculate mean suitability. A single integer is
  applied to every predator, or a vector of length `nspp` sets a
  distinct start year per predator. Defaults to `styr` in `data_list`.
  Used when diet data were sampled from a subset of years.

- suit_endyr:

  The last year used to calculate mean suitability. A single integer is
  applied to every predator, or a vector of length `nspp` sets a
  distinct end year per predator. Defaults to `endyr` in `data_list`.
  Used when diet data were sampled from a subset of years.

- fit_control:

  A list returned by
  [`fit_control()`](https://afsc-assessments.github.io/Rceattle/reference/fit_control.md)
  that bundles the optimizer / sdreport / phasing knobs (`phase`,
  `bias.correct`, `getsd`, `getJointPrecision`, `getReportCovariance`,
  `use_gradient`, `rel_tol`, `loopnum`, `newtonsteps`, `TMBfilename`,
  `verbose`, `nlminb_control`). Defaults to
  [`fit_control()`](https://afsc-assessments.github.io/Rceattle/reference/fit_control.md).
  See
  [`fit_control()`](https://afsc-assessments.github.io/Rceattle/reference/fit_control.md)
  for the meaning and defaults of each field.

- config:

  (Optional) An `Rceattle_run_config` from
  [`load_config()`](https://afsc-assessments.github.io/Rceattle/reference/load_config.md)
  or
  [`run_config()`](https://afsc-assessments.github.io/Rceattle/reference/run_config.md)
  whose stored settings overlay the ones you did not pass; `NULL`
  (default) applies no configuration. See Details for what it overlays.

- quiet_data_check:

  Drop the warnings the fit-time validation raises (errors still stop
  the fit). `FALSE` (default) for an ordinary fit. The diagnostic refits
  ,
  [`retrospective()`](https://afsc-assessments.github.io/Rceattle/reference/retrospective.md),
  [`jitter()`](https://afsc-assessments.github.io/Rceattle/reference/jitter.md),
  [`self_test()`](https://afsc-assessments.github.io/Rceattle/reference/self_test.md),
  [`profile()`](https://rdrr.io/r/stats/profile.html),
  [`run_mse()`](https://afsc-assessments.github.io/Rceattle/reference/run_mse.md),
  [`remove_F()`](https://afsc-assessments.github.io/Rceattle/reference/remove_F.md),
  [`sample_rec()`](https://afsc-assessments.github.io/Rceattle/reference/sample_rec.md),
  [`reweight_comps()`](https://afsc-assessments.github.io/Rceattle/reference/reweight_comps.md),
  set it, since they re-validate a `data_list` the caller has already
  fitted once and would otherwise repeat the same warnings per peel,
  jitter, or MSE iteration. Also drops a linkage filter's "has no
  effect" warning; a filter that drops its whole spec still warns.
  Convergence and TMB warnings are unaffected.

- ...:

  Deprecated optimizer / sdreport / phasing arguments (e.g. `phase`,
  `getsd`, `bias.correct`, `use_gradient`, `rel_tol`, `control`,
  `getJointPrecision`, `getReportCovariance`, `loopnum`, `newtonsteps`,
  `verbose`, `TMBfilename`). These are forwarded into `fit_control` with
  a deprecation warning; pass them via
  [`fit_control()`](https://afsc-assessments.github.io/Rceattle/reference/fit_control.md)
  instead.

## Value

A list of class "Rceattle" including:

- data_list: List of data inputs

- initial_params: List of starting parameters

- bounds: Parameter bounds used for estimation

- map: List of map used in TMB

- obj: TMB model object

- opt: Optimized model object from `nlminb`

- sdrep: Object of class `sdreport` exported by TMB including the
  standard errors of estimated parameters

- estimated_params: List of estimated parameters

- quantities: Derived quantities from CEATTLE

- run_time: Model run time

## Details

CEATTLE is an age-structured population dynamics model that can be fit
with or without predation mortality. The default is to exclude predation
mortality by setting `msmMode` to 0. Predation mortality can be included
by setting `msmMode` with the following options:

- 0\. Single species mode

- 1\. Holsman et al. 2015 predation based on multi-species virtual
  population analysis (MSVPA) based predation formation.

- 2\. MSVPA Holling Type III

Values 3 through 9 (Kinzey & Punt 2009 functional responses, Holling
Type I/II/III, predator interference, predator preemption,
Hassell-Varley, Ecosim) are blocked at runtime by
[`data_check()`](https://afsc-assessments.github.io/Rceattle/reference/data_check.md)
because the implementations have not been validated against the current
parameter set. See `src/TMB/predation.hpp`.

**What `config` overlays.** Two overlays happen, at different levels.
The estimation controls (`estimateMode`, `random_rec` / `random_q` /
`random_sel`, `suit_styr` / `suit_endyr`, `fit_control`) overlay only
the arguments you did not pass, so an explicit argument always wins. The
stored `model_config` is merged into the data object's **field by
field**, not wholesale: only the fields the config actually set are
imposed, and the data object keeps the rest. Since 5.36.0 a config built
with
[`model_config()`](https://afsc-assessments.github.io/Rceattle/reference/model_config.md)
therefore no longer drops the linkages held on the data object, which it
did when the whole structure was replaced. Where a field is set on both
and the two disagree, the config's value is used and the difference is
reported as a warning naming the field. A config written before that
field record existed is treated as having set its non-default fields.

## Initial age structure

What `initMode` estimates, and from what:

- `"FreeParams"` (0):

  The initial age-structure is estimated directly, one free parameter
  per age.

- `"Equilibrium"` (1):

  Unfished (\\F\_{init} = 0\\) equilibrium age-structure, projected from
  \\R_0\\ and residual natural mortality \\M1\\.

- `"NonEquilibrium"` (2):

  As (1), plus estimated initial population deviates, so the first year
  need not sit at equilibrium. The default.

- `"FishedNonEquilibrium"` (3):

  As (2), with an estimated initial fishing mortality \\F\_{init}\\
  added to the mortality that shapes the first year.

- `"FishedNonEquilibriumScaled"` (4):

  As (3), but \\F\_{init}\\ scales \\R_0\\ rather than entering the
  mortality directly.

- `"OffsetEquilibrium"` (5):

  Unfished equilibrium seeded by the first year's recruitment,
  `R_init * exp(rec_dev[1])`, decayed by \\M1\\ and closed with the
  usual geometric plus group. Initial deviates are turned off and no
  init-dev penalty is applied. This is the Cole Monnahan / AFSC GOA
  pollock convention.

Modes 1 and 5 differ by exactly one term: both start from the initial
equilibrium recruitment \\R\_{init}\\, but (1) projects it forward
unchanged while (5) seeds the first year with the realized recruitment
`R_init * exp(rec_dev[1])`. On a stock whose first year was not average,
that is not a small difference.

\\R\_{init}\\ is not always \\R_0\\: under a random-about-mean
stock-recruit relationship the two are equal, but under Beverton-Holt or
Ricker \\R\_{init}\\ is the equilibrium recruitment implied by the
curve.

**In multispecies mode the decay uses \\M1\\ only**, so predation
mortality (\\M2\\) does not enter the initial age structure under any of
these modes.

## Examples

``` r
# \donttest{
data(BS2017SS)
ss_run <- fit_mod(
  data_list    = BS2017SS,
  estimateMode = 0,
  msmMode      = 0,
  fit_control  = fit_control(phase = FALSE, verbose = 0)
)
#> `age_trans_matrix` data does not span range of age for species 1 will fill with 0s
# }
```
