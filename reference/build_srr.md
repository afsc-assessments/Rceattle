# Specify the stock-recruit relationship (SRR) for Rceattle

Sets the stock-recruit curve and how recruitment is estimated. Priors,
fixed values and environmental effects on `R0`, alpha, beta and `R_init`
go through `linkages`; see **Priors, fixed values and covariates**
below.

**Stock recruitment relationships currently implemented in Rceattle:**

- `srr_fun = 0` or `"mean"`: No stock recruit relationship. Recruitment
  is a function of \\R0\\ (on the log scale) and annual deviates (i.e.
  steepness = 0.99). \$\$R_y = exp(R0 + R\_{dev,y})\$\$

- `srr_fun = 2` or `"BevertonHolt"`: Beverton-Holt stock-recruitment
  relationship \$\$R_y = \frac{\alpha\_{srr} \*
  SB\_{y-minage}}{1+\beta\_{srr} \* SB\_{y-minage}}\$\$

- `srr_fun = 4` or `"Ricker"`: Ricker stock-recruitment relationship
  \$\$R_y = \alpha\_{srr} \* SB\_{y-minage} \* exp(-\beta\_{srr} \*
  SB\_{y-minage})\$\$

The Beverton-Holt and Ricker curves above are the deterministic mean;
realized recruitment applies the annual log deviation, \\R_y \cdot
exp(R\_{dev,y})\\, as in the mean form. For numerical stability the
Ricker \\\beta\_{srr}\\ is estimated on a scale divided by 1,000,000, so
the fitted `beta` is 1e6 times the density-dependence coefficient in the
equation above; `Bmsy_lim` (\\\approx 1/\beta\_{srr}\\) holds the same
scaling.

When `srr_pred_fun > 0` and `srr_fun = 0` recruitment in the hindcast is
estimated as in `srr_fun = 0` \$\$R_y = exp(R0 + R\_{dev,y})\$\$, but an
additional stock recruitment relationship defined by `srr_pred_fun` is
estimated between `srr_hat_styr` and `srr_hat_endyr` and treated as an
additional penalty. The stock recruitment relationship defined by
`srr_pred_fun` is then used in the projection.

**Multispecies models.** Spawning biomass per recruit is undefined when
mortality includes predation, so under `msmMode > 0` a curve fitted in
the hindcast estimates its initial recruitment level (`R0`) rather than
deriving it; under the fished initModes (3, 4) that level trades off
against the initial F. The curve can also enter as the penalty above. A
steepness prior is refused and `steepness` is reported as 0; priors on
alpha or beta go through `linkages`.

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

  Stock-recruit function used in the hindcast estimation (see the list
  below). Default = 0

- srr_pred_fun:

  Stock-recruit function used for projection, reference points, and
  penalties (see below). When `srr_fun == 0`, the stock-recruit curve is
  added as a penalty on the annually estimated hindcast recruitment
  (following AMAK and Jim Ianelli's pollock model). If `srr_fun > 0`,
  then `srr_pred_fun = srr_fun` and no extra penalty is added.

- proj_mean_rec:

  Recruitment used in the projection: `TRUE`/1 (default) = mean
  recruitment, the average R over the hindcast; `FALSE`/0 = the
  stock-recruit relationship given by `srr_pred_fun`. Equilibrium and
  dynamic reference points follow the curve whenever `srr_pred_fun` is a
  stock-recruit form, regardless of this switch.

- srr_mse_switchyr:

  Year at which an MSE switches from the annual recruitment-penalty
  estimate to the stock-recruit function (the `srr_fun = 0`,
  `srr_pred_fun > 0` case).

- srr_hat_styr:

  Integer. First year used to estimate the recruitment-penalty function
  (the AMAK/Ianelli penalty, active when `srr_pred_fun > 0` and
  `srr_fun = 0`), starting at `styr + 1`. Defaults to `styr + 1` in
  `data_list`. Useful when the environmental data conditioning the
  stock-recruit relationship is not available until the terminal year
  but projections are still wanted.

- srr_hat_endyr:

  Integer. Last year used to estimate the recruitment-penalty function
  (the AMAK/Ianelli penalty, active when `srr_pred_fun > 0` and
  `srr_fun = 0`). Defaults to `endyr` in `data_list`. Useful when the
  environmental data conditioning the stock-recruit relationship does
  not span the full time series but projections are still wanted.

- srr_est_mode:

  The curve's built-in prior, as a code or string: 1 / `"Estimated"`
  (default, no prior), or 2 / `"LognormalPrior"` or 3 / `"BetaPrior"` on
  Beverton-Holt steepness (single-species only); 0 / `"Fixed"` (alpha
  held at `srr_prior`) and `"LognormalPrior"` on a Ricker curve (a prior
  on alpha) are deprecated in favour of an `alpha` linkage.

- srr_prior:

  Natural-scale prior mean (the median of a lognormal prior when
  `bias_adjust_proc = FALSE`) of Beverton-Holt steepness under modes 2
  and 3; its other uses, as a Ricker alpha prior, a fixed alpha or
  alpha's starting value, are deprecated (see **Starting values**).

- srr_prior_sd:

  Prior standard deviation: log scale for the lognormal prior (mode 2),
  natural scale for the beta prior (mode 3).

- srr_alpha_init, srr_beta_init:

  Optional starting values for alpha and beta, natural scale, one per
  species, used only when the curve estimates them (see **Starting
  values**).

- srr_indices:

  Defunct: supplying it is an error; express an environmental effect
  through `linkages`.

- Bmsy_lim:

  Upper limit for Ricker based SSB-MSY (e.g 1/Beta). Will add a
  likelihood penalty if beta is estimated above this limit. Default `NA`
  is not used.

- linkages:

  Named list of
  [`linkage_spec()`](https://afsc-assessments.github.io/Rceattle/reference/linkage_spec.md)
  objects keyed by `"R0"`, `"alpha"`, `"beta"` or `"R_init"`: the
  recommended way to put a prior on, fix, or add an environmental effect
  to those parameters (see **Priors, fixed values and covariates**).

## Value

A `list` containing the stock recruitment relationship settings

## Priors, fixed values and covariates

Use `linkages` for `R0`, alpha, beta and `R_init`. Each entry is a
[`linkage_spec()`](https://afsc-assessments.github.io/Rceattle/reference/linkage_spec.md),
and an intercept-only formula (`~ 1`) acts on the parameter itself:

- **Prior:**
  `` priors = list(`(Intercept)` = prior_lognormal(log(m), s)) `` is
  lognormal with mean `m` (median `m` when `bias_adjust_proc = FALSE`)
  and log-scale SD `s`;
  [`prior_normal()`](https://afsc-assessments.github.io/Rceattle/reference/prior_normal.md)
  is normal on the natural scale.

- **Fixed value:** `` init = list(`(Intercept)` = v), est_phase = 0 ``
  holds the parameter at `v`, over supplied `inits` too.

- **One species:** add `species = 1`; the default applies to every
  species.

- **Environmental effect:** a covariate formula such as `~ temp` adds a
  log-scale effect by year; see
  [`vignette("environmental-linkages-and-priors")`](https://afsc-assessments.github.io/Rceattle/articles/environmental-linkages-and-priors.md).

A linkage on `R0` acts under mean recruitment, penalty form included.
Under a curve fitted in the hindcast a single-species `R0` is derived
from alpha and beta, so an `R0` linkage is refused; in a multispecies
model `R0` is the initial recruitment level, and only an intercept-only
linkage is accepted.

For a Ricker curve the lognormal linkage prior on alpha gives the same
objective as `srr_est_mode = "LognormalPrior"`, which is deprecated;
using both is refused. `srr_est_mode` and `srr_prior` remain for the one
prior a linkage cannot express, on Beverton-Holt steepness, which needs
spawning biomass per recruit and so exists only in single-species
models.

## Starting values

Mean recruitment (`R0`) starts at \\e^9 = 8103\\ thousand fish; under a
curve fitted in the hindcast of a multispecies model the same slot is
the free initial recruitment level `R_init`, with the same start. Alpha
starts at `srr_prior` (default 4) wherever that is an alpha, and at
\\e^3\\ otherwise; beta starts at 3. None of them knows the stock's
scale. Set them with `srr_alpha_init` / `srr_beta_init` or a linkage
`init`; supplying `srr_prior` as alpha's starting value is deprecated.
\\\beta\\ sets the density dependence in \\R = \alpha S / (1 + \beta
S)\\, so it must be on the order of \\(\alpha - 1/\phi_0) / R_0\\,
typically \\10^{-3}\\ or smaller for a stock measured in tonnes;
starting three orders of magnitude away drives predicted recruitment to
near zero and the optimizer returns `NA/NaN gradient evaluation`. From a
steepness \\h\\ and unfished spawning biomass per recruit \\\phi_0\\:
\$\$\alpha = \frac{4h}{\phi_0 (1 - h)}, \qquad \beta = \frac{\alpha -
1/\phi_0}{R_0}.\$\$ Under predation, where \\\phi_0\\ is undefined, seed
them from a curve fitted to an earlier fit's SSB-recruitment pairs.

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

# Beverton-Holt with alpha fixed at 20.
build_srr(srr_fun = "BevertonHolt",
          linkages = list(alpha = linkage_spec(~ 1, est_phase = 0,
            init = list(`(Intercept)` = 20))))
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
#>   link:    log
#> 
#> 

# A temperature effect on alpha (log scale; BTempC is a column of env_data).
build_srr(srr_fun = "BevertonHolt",
          linkages = list(alpha = linkage_spec(~ BTempC)))
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
#>   formula: ~BTempC
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
