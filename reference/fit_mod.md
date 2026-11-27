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
  or taken from a previous fit (`model$estimated_params`). `NULL`
  initializes from scratch.

- map:

  (Optional) A map object from
  [`build_map`](https://afsc-assessments.github.io/Rceattle/reference/build_map.md).

- bounds:

  (Optional) A bounds object from
  [`build_bounds`](https://afsc-assessments.github.io/Rceattle/reference/build_bounds.md).

- file:

  (Optional) Filename to save to. If NULL, no file is saved.

- estimateMode:

  What to fit, as a string alias or integer code: `"Estimate"` (0) fits
  the hindcast and the `HCR` projection, `"Hindcast"` (1) the hindcast
  alone, `"Projection"` (2) the reference points and projection alone
  from `inits`, `"DebugBuild"` (3) builds through `MakeADFun` without
  optimizing, leaving a usable `obj$fn()` / `obj$gr()`, and
  `"DebugOptimize"` (4) maps every parameter out, so its objective is a
  placeholder rather than a likelihood. Default `"Estimate"`.

- random_rec:

  logical; treat recruitment deviations as random effects via the
  Laplace approximation. Default FALSE.

- random_q:

  logical; integrate the `Time_varying_q` deviations as random effects
  with one estimated sd per `Catchability_index` group instead of fixing
  them at `Time_varying_q_sd`. Default FALSE. Linkage random effects are
  integrated either way.

- random_sel:

  logical; the same for `Time_varying_sel` deviations, one estimated sd
  per `Selectivity_index` group. Default FALSE.

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

  Predation-mortality mode (string or integer): `"SingleSpecies"` (0,
  default; no predation), `"MSVPA"` (1; Type-II MSVPA, Holsman et al.
  2016), or `"TypeIIIMSVPA"` (2); codes 3–9 (Kinzey & Punt 2009
  functional responses: Holling forms through Ecosim) are declared but
  refused by
  [`data_check()`](https://afsc-assessments.github.io/Rceattle/reference/data_check.md)
  until validated against the current parameter set.

- avgnMode:

  the average abundance-at-age approximation in the predation equations.
  Only mode 0, \\N/Z(1 - exp(-Z))\\ (the MSVPA form), is implemented;
  the alternatives are declared but have no effect.

- initMode:

  how the population is initialized, as a string alias or integer code:
  `"FreeParams"` (0), `"Equilibrium"` (1), `"NonEquilibrium"` (2, the
  default), `"FishedNonEquilibrium"` (3), `"FishedNonEquilibriumScaled"`
  (4), `"OffsetEquilibrium"` (5), `"FishedNonEquilibriumSelected"` (6).
  What each one estimates, and how \\F\_{init}\\ reaches an age under 3,
  4 and 6, is tabulated in
  [`vignette("model-options-and-functionality")`](https://afsc-assessments.github.io/Rceattle/articles/model-options-and-functionality.md).
  In multispecies mode the initial decay uses \\M1\\ only, so predation
  does not shape the initial age structure under any mode.

- suitMode:

  Predator-prey suitability mode per predator (scalar or length-`nspp`
  vector): 0 = empirical from diet data (Holsman et al. 2016), 2 =
  weight-based gamma, 4 = weight-based lognormal, 6 = weight-based
  normal (length-based codes 1, 3, 5 are declared but refused by
  [`data_check()`](https://afsc-assessments.github.io/Rceattle/reference/data_check.md)).

- suit_styr:

  The first year used to calculate mean suitability, as one integer for
  every predator or a vector of length `nspp`. Defaults to `styr`. Set
  it where diet data were sampled from a subset of years.

- suit_endyr:

  The last year used to calculate mean suitability, in the same form.
  Defaults to `endyr`.

- fit_control:

  A list from
  [`fit_control()`](https://afsc-assessments.github.io/Rceattle/reference/fit_control.md)
  bundling the optimizer, sdreport and phasing settings. Defaults to
  [`fit_control()`](https://afsc-assessments.github.io/Rceattle/reference/fit_control.md),
  which documents each field.

- config:

  (Optional) An `Rceattle_run_config` from
  [`load_config()`](https://afsc-assessments.github.io/Rceattle/reference/load_config.md)
  or
  [`run_config()`](https://afsc-assessments.github.io/Rceattle/reference/run_config.md)
  whose stored settings overlay the ones you did not pass; `NULL`
  (default) applies none. See Details.

- quiet_data_check:

  Drop the warnings fit-time validation raises; errors still stop the
  fit. `FALSE` (default) for an ordinary fit. The diagnostic refits
  ([`retrospective()`](https://afsc-assessments.github.io/Rceattle/reference/retrospective.md),
  [`jitter()`](https://afsc-assessments.github.io/Rceattle/reference/jitter.md),
  [`self_test()`](https://afsc-assessments.github.io/Rceattle/reference/self_test.md),
  [`profile()`](https://rdrr.io/r/stats/profile.html),
  [`run_mse()`](https://afsc-assessments.github.io/Rceattle/reference/run_mse.md),
  [`remove_F()`](https://afsc-assessments.github.io/Rceattle/reference/remove_F.md),
  [`sample_rec()`](https://afsc-assessments.github.io/Rceattle/reference/sample_rec.md),
  [`reweight_comps()`](https://afsc-assessments.github.io/Rceattle/reference/reweight_comps.md))
  set it, since they re-validate a `data_list` the caller has already
  fitted and would otherwise repeat each warning per peel or iteration.
  It also drops a linkage filter's "has no effect" warning, though a
  filter that drops its whole spec still warns. Convergence and TMB
  warnings are unaffected.

- ...:

  Deprecated optimizer / sdreport / phasing arguments, forwarded into
  `fit_control` with a deprecation warning. Pass them via
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

- sdrep: Object of class `sdreport` from TMB, with standard errors

- estimated_params: List of estimated parameters

- quantities: Derived quantities from CEATTLE

- run_time: Model run time

## Details

CEATTLE is an age-structured population model, fit with or without
predation mortality; `msmMode` selects which.

**What `config` overlays**, at two levels. The estimation controls
(`estimateMode`, `random_rec` / `random_q` / `random_sel`, `suit_styr` /
`suit_endyr`, `fit_control`) overlay only arguments you did not pass, so
an explicit argument always wins. The stored `model_config` merges into
the data object's **field by field**: only fields the config actually
set are imposed, so a config built with
[`model_config()`](https://afsc-assessments.github.io/Rceattle/reference/model_config.md)
keeps the linkages held on the data object. Where both set a field and
disagree, the config wins and the difference is warned about by name. A
config written before that field record existed counts as having set its
non-default fields.

## References

Holsman, K.K., Ianelli, J., Aydin, K., Punt, A.E., and Moffitt, E.A.
2016. A comparison of fisheries biological reference points estimated
from temperature-specific multi-species and single-species
climate-enhanced stock assessment models. Deep Sea Res. II 134:360-378.

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
