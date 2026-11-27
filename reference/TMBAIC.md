# Calculate marginal AIC for a fitted model

`TMBAIC` calculates AIC for a given model fit

## Usage

``` r
TMBAIC(opt, p = 2, n = Inf)
```

## Arguments

- opt:

  the output from `nlminb` or `optim`

- p:

  the penalty on additional fixed effects (default=2, for AIC)

- n:

  the sample size, for use in AICc calculation (default=Inf, for which
  AICc=AIC)

## Value

AIC, where a parsimonious model has a AIC relative to other candidate
models

## Author

James Thorson; copied from
[`TMBhelper::TMBAIC()`](https://rdrr.io/pkg/TMBhelper/man/TMBAIC.html)
(<https://github.com/kaskr/TMB_contrib_R>, GPL-3).
