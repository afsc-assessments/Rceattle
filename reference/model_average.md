# Model average of derived quantities

Model average of derived quantities

## Usage

``` r
model_average(
  object = NULL,
  weights = NULL,
  uncertainty = FALSE,
  nboot = 10000,
  Rceattle = NULL
)
```

## Arguments

- object:

  list of Rceattle model objects

- weights:

  vector of weights to be used for weighting models

- uncertainty:

  TRUE/FALSE Sample uncertainty across derived quantities using weighted
  bootstrap from the asymptotic distribution of MLEs

- nboot:

  Number of bootstraps taken from asymptotic distribution of MLEs.
  Default = 10000

- Rceattle:

  deprecated name for `object`, still accepted so existing scripts keep
  working. Supplying both is an error.

## Value

an Rceattle object with derived quantities weighted by the specified
weights. The length of the derived quantities spans the years which
overlap across all models.

## Details

With `uncertainty = TRUE` each fit's parameters are drawn from its
asymptotic covariance, so every fit needs an sdreport
(`fit_control(getsd = TRUE)`) and a `fit$obj` in which the hindcast
parameters are free. An `estimateMode = "Estimate"` fit under any
harvest control rule but `"NoFishing"` leaves only the rule's F free in
`fit$obj` unless it used `projection_uncertainty`, and a `"Projection"`
fit with the default map holds every hindcast parameter fixed, so such
fits are refused. Refit with `estimateMode = "Hindcast"`, which projects
with no fishing, or, for a fixed-effects fit under any rule but
`"ConstantF"`, with `estimateMode = "Estimate"` and
`fit_control(projection_uncertainty = TRUE)`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Equal weights across two candidate models.
model_average(list(fit1, fit2), weights = c(1, 1))
# Weighted by AIC, with bootstrapped uncertainty.
model_average(list(fit1, fit2), weights = c(0.7, 0.3), uncertainty = TRUE)
} # }
```
