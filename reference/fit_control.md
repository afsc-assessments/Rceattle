# Bundle the optimizer / sdreport / phasing controls for `fit_mod()`

Collects the optimizer and reporting settings into one object, so a call
to
[`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)
stays about the model rather than how it is fit. Pass the result through
[`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)'s
`fit_control` argument, where it overrides the equivalent individual
arguments, which are kept for back-compatibility.

## Usage

``` r
fit_control(
  bias.correct = FALSE,
  getsd = TRUE,
  getJointPrecision = TRUE,
  getReportCovariance = FALSE,
  projection_uncertainty = FALSE,
  selectivity_se = FALSE,
  comp_offset = NULL,
  bias_adjust_obs = TRUE,
  bias_adjust_proc = TRUE,
  use_gradient = TRUE,
  rel_tol = 1,
  loopnum = 5,
  newtonsteps = 0,
  phase = FALSE,
  TMBfilename = NULL,
  verbose = 1,
  nlminb_control = list(eval.max = 1e+09, iter.max = 1e+09, trace = 0)
)
```

## Arguments

- bias.correct:

  Logical (`FALSE` by default): apply TMB's epsilon bias correction
  (Thorson and Kristensen 2016) in
  [`TMB::sdreport()`](https://rdrr.io/pkg/TMB/man/sdreport.html), as
  distinct from the lognormal bias adjustment of `bias_adjust_obs` /
  `bias_adjust_proc`.

- getsd:

  logical. Run
  [`TMB::sdreport()`](https://rdrr.io/pkg/TMB/man/sdreport.html) after
  optimization. Default `TRUE`.

- getJointPrecision:

  logical. Return the full Hessian of fixed and random effects. Default
  `TRUE`.

- getReportCovariance:

  logical. Return the variance-covariance of `ADREPORT` variables.
  Default `FALSE`.

- projection_uncertainty:

  logical. Carry hindcast parameter uncertainty into an HCR projection,
  by refitting with the hindcast and biological-reference-point
  parameters turned on. Default `FALSE` for speed.

- selectivity_se:

  logical. Also return a standard error for log selectivity-at-age.
  Default `FALSE`; see Details.

- comp_offset:

  Numeric or `NULL`. Small proportion added to the observed and
  predicted composition and CAAL bins before the multinomial /
  Dirichlet-multinomial likelihood, so an empty bin is not `log(0)`.
  Stored on `data_list`, so the fitted likelihood, the OSA observation
  vector and any internal refit share it. `NULL` (default) inherits
  `data_list$comp_offset` if set, else `1e-5`. Diet compositions take
  the same offset.

- bias_adjust_obs:

  logical, default `TRUE`. Apply a bias adjustment (mean - sd^2/2) to
  lognormal data likelihoods.

- bias_adjust_proc:

  logical, default `TRUE`. Shift lognormal process likelihoods,
  lognormal priors and the Ianelli stock-recruit penalty by `-sd^2/2` on
  the log scale, making each prior value or curve a mean rather than a
  median. A value between 0 and 1 scales the shift.

- use_gradient:

  logical. Use the analytic gradient during phasing. Default `TRUE`.

- rel_tol:

  Numeric tolerance for flagging a discontinuous likelihood, comparing
  `nlminb`'s objective with a fresh evaluation at the estimates (under
  `newtonsteps > 0` the reported objective is that evaluation, so
  nothing is flagged). Default `1`.

- loopnum:

  Integer. Times to restart optimization; `3` sometimes reaches a lower
  final gradient than `1`. Default `5`.

- newtonsteps:

  Integer. Extra unconstrained Newton steps after optimization, an
  alternative to `loopnum`; a singular Hessian (usually a parameter the
  data do not inform) ends them with a warning, and with TMBhelper
  installed the fit is redone without Newton steps. Default `0`.

- phase:

  `TRUE`/`FALSE`, or a list of parameter names with their phases.
  Default `FALSE`.

- TMBfilename:

  Optional path (without `.cpp`) to an alternate TMB template for
  development. Default `NULL` uses the bundled `ceattle`.

- verbose:

  `0` silent, `1` model-fit updates, `2` model-fit and TMB progress.
  Default `1`.

- nlminb_control:

  Control parameters passed to
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html); see
  [`?nlminb`](https://rdrr.io/r/stats/nlminb.html). Default
  `list(eval.max = 1e9, iter.max = 1e9, trace = 0)`.

## Value

A list of class `"Rceattle_fit_control"`.

## Selectivity standard errors

`selectivity_se` needs `getsd = TRUE`, and reports the error belonging
to whichever `sdreport` the fit ends on. Under
`estimateMode = "Estimate"` with an HCR that re-optimizes, that is the
*projection* fit, where every selectivity parameter is mapped off, so
every error comes back exactly 0; `"Projection"` estimates no
selectivity at all. Use `"Hindcast"`, or `"Estimate"` with
`projection_uncertainty = TRUE`, to get an error from a fit that
estimated the curve.
[`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)
warns in each of these cases.

It is off by default because the delta method forms a Jacobian of every
reported value against every parameter, so the cost is their product: on
`Atka2022`, 1,012 values against 584 parameters. The error is on the log
scale rather than the logit, because the non-parametric forms normalize
to mean selectivity 1 instead of a maximum of 1 – 58% of `Atka2022`'s
`sel_at_age` exceeds 1, up to 3.06 – so a logit is undefined over most
of the array.

Rows cover estimated, age-based lead fleets only, starting at each
fleet's first selected age (`Age_first_selected`, or a deprecated
`Bin_first_selected`). Four kinds of cell are a structural zero (a
`Fixed` fleet's empirical curve, a length-based fleet's growth-matrix
projection, an age below the first selected age, and array padding), and
one `log(0) = -Inf` on the tape turns *every* quantity in the `sdreport`
to `NaN`, biomass and SSB included. All four are identified from the
data, so the reported set never depends on a parameter value, and a
value that underflows to zero is floored so it cannot reintroduce the
`-Inf`. See
[`plot_selectivity()`](https://afsc-assessments.github.io/Rceattle/reference/plot_selectivity.md),
which draws the interval.

## References

Thorson, J.T., and Kristensen, K. 2016. Implementing a generic method
for bias correction in statistical models using random effects, with
spatial and population dynamics examples. Fish. Res. 175:66-74.

## Examples

``` r
# Quick-and-dirty fit: skip sdreport, single optimizer pass
ctl <- fit_control(getsd = FALSE, loopnum = 1)

# Production fit with bias correction and joint precision
ctl <- fit_control(bias.correct = TRUE, getJointPrecision = TRUE)
```
