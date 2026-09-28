# Bundle the optimizer / sdreport / phasing controls for `fit_mod()`

[`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)
holds roughly a dozen optimizer- and reporting-related arguments
(`bias.correct`, `getsd`, `loopnum`, `newtonsteps`, ...). That is a lot
of surface area when the user mostly cares about "what model am I
fitting" rather than "how is it being fit." `fit_control()` collects
those knobs into a single object so calls to
[`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)
can stay focused on the model spec:

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

  logical. If `TRUE`, applies bias correction via
  [`TMB::sdreport()`](https://rdrr.io/pkg/TMB/man/sdreport.html).
  Default `FALSE`.

- getsd:

  logical. If `TRUE`, run
  [`TMB::sdreport()`](https://rdrr.io/pkg/TMB/man/sdreport.html) after
  optimization. Default `TRUE`.

- getJointPrecision:

  logical. Return the full Hessian of fixed and random effects. Default
  `TRUE` (matches
  [`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)
  default).

- getReportCovariance:

  logical. Return the variance-covariance of `ADREPORT` variables.
  Default `FALSE`.

- projection_uncertainty:

  logical. If `TRUE`, accounts for hindcast parameter uncertainty in
  projections when using an HCR (refits with all hindcast and
  biological-reference-point parameters turned on). Default `FALSE` for
  speed.

- selectivity_se:

  logical. If `TRUE`,
  [`TMB::sdreport()`](https://rdrr.io/pkg/TMB/man/sdreport.html) also
  returns a standard error for log selectivity-at-age. Default `FALSE`.

- comp_offset:

  Numeric or `NULL`. Small proportion offset added to the observed and
  predicted age/length composition and conditional-age-at-length (caal)
  bins before the multinomial / Dirichlet-multinomial likelihood (to
  avoid `log(0)` for empty bins). Stored on `data_list` so the fitted
  likelihood and the OSA observation vector use the same offset, and so
  internal re-fits inherit it. `NULL` (the default) inherits
  `data_list$comp_offset` if set, otherwise `1e-5`; a number overrides
  it. Does not apply to diet (stomach-content) compositions.

- bias_adjust_obs:

  logical with default TRUE. Whether to apply a bias adjustment
  (mean-sd^2/2) to lognormal data likelihoods

- bias_adjust_proc:

  logical with default TRUE. Whether lognormal process likelihoods,
  lognormal priors and the Ianelli stock-recruit penalty are shifted by
  `-sd^2/2` on the log scale, making each prior value or curve a mean
  rather than a median; a value between 0 and 1 scales the shift.

- use_gradient:

  logical. Use the analytic gradient during phasing. Default `TRUE`.

- rel_tol:

  Numeric tolerance used to flag discontinuous likelihood warnings,
  comparing `nlminb`'s objective with a fresh evaluation of the object
  it came from. Default `1`.

- loopnum:

  Integer. Number of times to re-start optimization (`loopnum = 3`
  sometimes achieves a lower final gradient than `loopnum = 1`). Default
  `5`.

- newtonsteps:

  Integer. Number of extra unconstrained Newton steps to take after
  optimization (alternative to `loopnum`). Default `0`.

- phase:

  `TRUE`/`FALSE` or a list. If `FALSE`, the model is not phased. If
  `TRUE`, default phasing is used. Can also accept a list of parameter
  object names with corresponding phase. Default `FALSE`.

- TMBfilename:

  Optional character. Path (without `.cpp`) to an alternate TMB template
  for development. Default `NULL` (use the bundled `ceattle`).

- verbose:

  `0` = silent, `1` = print updates of model fit, `2` = print updates of
  model fit and TMB estimation progress. Default `1`.

- nlminb_control:

  A list of control parameters passed to
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html). See
  [`?nlminb`](https://rdrr.io/r/stats/nlminb.html). Default
  `list(eval.max = 1e9, iter.max = 1e9, trace = 0)`.

## Value

A list of class `"Rceattle_fit_control"`.

## Details

    fit <- fit_mod(
      data_list   = BS2017SS,
      msmMode     = 0,
      fit_control = fit_control(loopnum = 1, getsd = FALSE)
    )

Pass the result via the `fit_control` argument to
[`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md).
When supplied, the values in the `fit_control` object override the
corresponding individual arguments to
[`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md).
Individual arguments are kept for backward compatibility.

## Selectivity standard errors

`selectivity_se` needs `getsd = TRUE`, and the error it reports is the
one belonging to whichever `sdreport` the fit ends on. Under
`estimateMode = "Estimate"` with any HCR that re-optimizes, that is the
*projection* fit, in which every selectivity parameter is mapped off, so
every error comes back exactly 0. `estimateMode = "Projection"`
estimates no selectivity at all and reports none where it runs no
`sdreport`. Use `estimateMode = "Hindcast"`, or `"Estimate"` with
`projection_uncertainty = TRUE`, to get an error from a fit that
estimated the curve.
[`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)
warns before fitting in each of these cases.

Off by default because the delta method forms a Jacobian of every
reported value against every parameter, so its cost is the product of
the two: on `Atka2022` it adds 1,012 values against 584 parameters.

The error is on the log scale, not the logit, because the non-parametric
forms normalize to mean selectivity 1 rather than a maximum of 1, 58% of
`Atka2022`'s `sel_at_age` exceeds 1, to 3.06, so a logit is undefined
over most of the array.

Rows cover estimated, age-based lead fleets only, and start at each
fleet's first selected bin. Four kinds of cell hold a structural zero, a
`Fixed` fleet's empirical curve, a length-based fleet's growth-matrix
projection, a bin below `Bin_first_selected`, and array padding, and one
`log(0) = -Inf` on the tape turns *every* quantity in the `sdreport` to
`NaN`, biomass and SSB included. All four are identified from the data,
so the reported set never depends on a parameter value; a value that
underflows to zero is floored, so it cannot reintroduce the `-Inf`. See
[`plot_selectivity()`](https://afsc-assessments.github.io/Rceattle/reference/plot_selectivity.md),
which draws the interval.

## Examples

``` r
# Quick-and-dirty fit: skip sdreport, single optimizer pass
ctl <- fit_control(getsd = FALSE, loopnum = 1)

# Production fit with bias correction and joint precision
ctl <- fit_control(bias.correct = TRUE, getJointPrecision = TRUE)
```
