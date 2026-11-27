# Specify the stock-recruit relationship (SRR) for Rceattle

Sets the stock-recruit curve and how recruitment is estimated.

**Relationships currently implemented:**

- `srr_fun = 0` or `"mean"`: no stock-recruit relationship. Recruitment
  is \\R0\\ (log scale) plus annual deviates, and the reported
  `steepness` holds its placeholder 0.99 rather than being estimated.
  \$\$R_y = exp(R0 + R\_{dev,y})\$\$

- `srr_fun = 2` or `"BevertonHolt"` (Beverton and Holt 1957): \$\$R_y =
  \frac{\alpha\_{srr} \* SB\_{y-minage}}{1+\beta\_{srr} \*
  SB\_{y-minage}}\$\$

- `srr_fun = 4` or `"Ricker"` (Ricker 1954): \$\$R_y = \alpha\_{srr} \*
  SB\_{y-minage} \* exp(-\beta\_{srr} \* SB\_{y-minage})\$\$

Both curves are the deterministic mean; realized recruitment applies the
annual log deviation, \\R_y \cdot exp(R\_{dev,y})\\. For numerical
stability the Ricker \\\beta\_{srr}\\ is estimated divided by 1,000,000,
so the fitted `beta` is 1e6 times the coefficient above. `Bmsy_lim` is
on the natural scale: the template divides it by 1e6 before comparing,
so supply it in metric tons.

With `srr_pred_fun > 0` and `srr_fun = 0`, the hindcast estimates
recruitment as in the mean form while the curve is fitted between
`srr_hat_styr` and `srr_hat_endyr` as a penalty, then used in the
projection.

**Multispecies models.** Spawning biomass per recruit is undefined once
mortality includes predation, so under `msmMode > 0` a curve fitted in
the hindcast estimates its initial recruitment level (`R0`) rather than
deriving it, and under the fished initModes (3, 4) that level trades off
against the initial F. A steepness prior is refused and `steepness` is
reported as 0; priors on alpha or beta go through `linkages`.

Two refusals worth knowing before you write a linkage: under a curve
fitted in the hindcast a single-species `R0` is derived from alpha and
beta, so an `R0` linkage is refused (a multispecies model accepts an
intercept-only one), and a Ricker lognormal linkage prior on alpha
duplicates the deprecated `srr_est_mode = "LognormalPrior"`, so using
both is refused.

## Usage

``` r
build_srr(
  srr_fun = 0,
  srr_pred_fun = srr_fun,
  proj_mean_rec = TRUE,
  srr_mse_switchyr = NULL,
  srr_hat_styr = NULL,
  srr_hat_endyr = NULL,
  srr_est_mode = 1,
  srr_prior = 4,
  srr_prior_sd = 1,
  srr_alpha_init = NULL,
  srr_beta_init = NULL,
  srr_indices = NA,
  Bmsy_lim = NA,
  linkages = NULL
)
```

## Arguments

- srr_fun:

  Stock-recruit function used in the hindcast estimation (see below).
  Default = 0

- srr_pred_fun:

  Stock-recruit function used for projection, reference points and
  penalties. When `srr_fun == 0` the curve is added as a penalty on the
  annually estimated hindcast recruitment (the AMAK / Jim Ianelli
  pollock form); when `srr_fun > 0` it equals `srr_fun` and no penalty
  is added.

- proj_mean_rec:

  Recruitment used in the projection: `TRUE`/1 (default) is mean
  recruitment over the hindcast, `FALSE`/0 follows `srr_pred_fun`.
  Equilibrium and dynamic reference points follow the curve whenever
  `srr_pred_fun` is a stock-recruit form, whatever this is set to.

- srr_mse_switchyr:

  Year at which an MSE switches from the annual recruitment-penalty
  estimate to the stock-recruit function.

- srr_hat_styr, srr_hat_endyr:

  First and last year used to estimate the recruitment-penalty function
  (active when `srr_pred_fun > 0` and `srr_fun = 0`), defaulting to
  `styr + 1` and `endyr`. Set them where the environmental data
  conditioning the curve do not span the full series but projections are
  still wanted.

- srr_est_mode:

  The curve's built-in prior, as a code or string: 1 / `"Estimated"`
  (default, no prior), or 2 / `"LognormalPrior"` or 3 / `"BetaPrior"` on
  Beverton-Holt steepness, single-species only. 0 / `"Fixed"` and a
  Ricker `"LognormalPrior"` are deprecated in favour of an `alpha`
  linkage.

- srr_prior:

  Natural-scale prior mean of Beverton-Holt steepness under modes 2 and
  3 (the median of a lognormal prior when `bias_adjust_proc = FALSE`).
  Its other uses, as a Ricker alpha prior, a fixed alpha, or alpha's
  starting value, are deprecated.

- srr_prior_sd:

  Prior standard deviation: log scale for the lognormal prior (mode 2),
  natural scale for the beta prior (mode 3).

- srr_alpha_init, srr_beta_init:

  Optional natural-scale starting values for alpha and beta, one per
  species, used only where the curve estimates them.

- srr_indices:

  Defunct: supplying it is an error; express an environmental effect
  through `linkages`.

- Bmsy_lim:

  Upper limit on Ricker SSB-MSY (1/beta) in metric tons, adding a
  likelihood penalty once estimated SSB-MSY rises above it. Default `NA`
  is unused.

- linkages:

  Named list of
  [`linkage_spec()`](https://afsc-assessments.github.io/Rceattle/reference/linkage_spec.md)
  objects keyed by `"R0"`, `"alpha"`, `"beta"` or `"R_init"` – note the
  key is a unitless log-scale MULTIPLIER on the initial age structure,
  not the reported quantity `R_init`, which is equilibrium recruitment
  at F = F_init in thousands of fish. It is the recommended way to put a
  prior on, fix, or add an environmental effect to those parameters. See
  [`vignette("environmental-linkages-and-priors")`](https://afsc-assessments.github.io/Rceattle/articles/environmental-linkages-and-priors.md).

## Value

A `list` containing the stock recruitment relationship settings

## Starting values

None of the starting values knows the stock's scale, and a bad one fails
the fit rather than fitting badly. Mean recruitment (`R0`) starts at
\\e^9 = 8103\\ thousand fish; alpha starts at `srr_prior` (default 4)
wherever that is an alpha and at \\e^3\\ otherwise; beta starts at 3.
Set them with `srr_alpha_init` / `srr_beta_init` or a linkage `init`.

\\\beta\\ sets the density dependence in \\R = \alpha S / (1 + \beta
S)\\, so it must be on the order of \\(\alpha - 1/\phi_0) / R_0\\,
typically \\10^{-3}\\ or smaller for a stock measured in tonnes.
Starting three orders of magnitude away drives predicted recruitment to
near zero and the optimizer returns `NA/NaN gradient evaluation`. From a
steepness \\h\\ and unfished spawning biomass per recruit \\\phi_0\\:
\$\$\alpha = \frac{4h}{\phi_0 (1 - h)}, \qquad \beta = \frac{\alpha -
1/\phi_0}{R_0}.\$\$ Under predation, where \\\phi_0\\ is undefined, seed
them from a curve fitted to an earlier fit's SSB-recruitment pairs.

## References

Beverton, R.J.H., and Holt, S.J. 1957. On the dynamics of exploited fish
populations. Fishery Investigations Series II, Vol. 19. UK Ministry of
Agriculture, Fisheries and Food, London.

Ricker, W.E. 1954. Stock and recruitment. J. Fish. Res. Board Can.
11:559-623.

## Examples

``` r
# Mean recruitment: no stock-recruit relationship fitted.
build_srr(srr_fun = "mean")
#> $srr_fun
#> [1] 0
#> 
#> $srr_pred_fun
#> [1] 0
#> 
#> $proj_mean_rec
#> [1] TRUE
#> 
#> $srr_mse_switchyr
#> NULL
#> 
#> $srr_hat_styr
#> NULL
#> 
#> $srr_hat_endyr
#> NULL
#> 
#> $srr_est_mode
#> [1] 1
#> 
#> $srr_prior
#> [1] 4
#> 
#> $srr_prior_sd
#> [1] 1
#> 
#> $srr_alpha_init
#> NULL
#> 
#> $srr_beta_init
#> NULL
#> 
#> $srr_indices
#> [1] NA
#> 
#> $Bmsy_lim
#> [1] -999
#> 
#> $linkages
#> NULL
#> 

# Beverton-Holt with a lognormal prior on alpha (mean 5, log-scale SD 0.5),
# species 1 only.
build_srr(srr_fun = "BevertonHolt",
          linkages = list(alpha = linkage_spec(~ 1, species = 1,
            priors = list(`(Intercept)` = prior_lognormal(log(5), 0.5)))))
#> $srr_fun
#> [1] 2
#> 
#> $srr_pred_fun
#> [1] 2
#> 
#> $proj_mean_rec
#> [1] TRUE
#> 
#> $srr_mse_switchyr
#> NULL
#> 
#> $srr_hat_styr
#> NULL
#> 
#> $srr_hat_endyr
#> NULL
#> 
#> $srr_est_mode
#> [1] 1
#> 
#> $srr_prior
#> [1] 4
#> 
#> $srr_prior_sd
#> [1] 1
#> 
#> $srr_alpha_init
#> NULL
#> 
#> $srr_beta_init
#> NULL
#> 
#> $srr_indices
#> [1] NA
#> 
#> $Bmsy_lim
#> [1] -999
#> 
#> $linkages
#> $linkages$alpha
#> <Rceattle linkage spec>
#>   param:   alpha
#>   formula: ~1
#>   prior:    (Intercept) ~ lognormal(1.60944, 0.5)
#>   species: 1
#>   link:    log
#> 
#> 

# Mean recruitment with the curve as a penalty (Ianelli form).
build_srr(srr_fun = "mean", srr_pred_fun = "BevertonHolt")
#> $srr_fun
#> [1] 0
#> 
#> $srr_pred_fun
#> [1] 2
#> 
#> $proj_mean_rec
#> [1] TRUE
#> 
#> $srr_mse_switchyr
#> NULL
#> 
#> $srr_hat_styr
#> NULL
#> 
#> $srr_hat_endyr
#> NULL
#> 
#> $srr_est_mode
#> [1] 1
#> 
#> $srr_prior
#> [1] 4
#> 
#> $srr_prior_sd
#> [1] 1
#> 
#> $srr_alpha_init
#> NULL
#> 
#> $srr_beta_init
#> NULL
#> 
#> $srr_indices
#> [1] NA
#> 
#> $Bmsy_lim
#> [1] -999
#> 
#> $linkages
#> NULL
#> 

# A prior on Beverton-Holt steepness, the one prior linkages cannot express
# (single-species only).
build_srr(srr_fun = "BevertonHolt", srr_est_mode = "LognormalPrior",
          srr_prior = 0.8, srr_prior_sd = 0.2)
#> $srr_fun
#> [1] 2
#> 
#> $srr_pred_fun
#> [1] 2
#> 
#> $proj_mean_rec
#> [1] TRUE
#> 
#> $srr_mse_switchyr
#> NULL
#> 
#> $srr_hat_styr
#> NULL
#> 
#> $srr_hat_endyr
#> NULL
#> 
#> $srr_est_mode
#> [1] 2
#> 
#> $srr_prior
#> [1] 0.8
#> 
#> $srr_prior_sd
#> [1] 0.2
#> 
#> $srr_alpha_init
#> NULL
#> 
#> $srr_beta_init
#> NULL
#> 
#> $srr_indices
#> [1] NA
#> 
#> $Bmsy_lim
#> [1] -999
#> 
#> $linkages
#> NULL
#> 
```
