# Print method for a retrospective analysis

Reports Mohn's rho against a reference band rather than as a bare
number. The default band is +/- 0.2 on SSB, the rule of thumb
[`vignette("model-diagnostics")`](https://afsc-assessments.github.io/Rceattle/articles/model-diagnostics.md)
states; Hurtado-Ferro et al. (2015) give the asymmetric,
life-history-dependent alternatives (-0.15 to 0.20 for long-lived, -0.22
to 0.30 for short-lived), which is why the band is an argument rather
than a constant.

Only the terminal-year peel (`Forecast year` 0) is judged. The
forecast-skill rows are reported for information: a rho computed over a
forecast horizon is not the quantity the +/- 0.2 rule was calibrated on.

## Usage

``` r
# S3 method for class 'Rceattle_retro'
print(x, band = 0.2, ...)
```

## Arguments

- x:

  A `"Rceattle_retro"` object from
  [`retrospective()`](https://afsc-assessments.github.io/Rceattle/reference/retrospective.md).

- band:

  Symmetric reference band for Mohn's rho. Default `0.2`.

- ...:

  Currently unused.

## Value

`x`, invisibly.

## References

Hurtado-Ferro, F., et al. 2015. Looking in the rear-view mirror: bias
and retrospective patterns in integrated, age-structured stock
assessment models. ICES J. Mar. Sci. 72:99-110.
