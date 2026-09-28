# Print method for fitted Rceattle models

Provides a compact summary so that auto-printing inside R Markdown /
knitr / RStudio does not recurse into the (very deep) data and TMB
objects stored on the fit. Only structural metadata, convergence status,
headline derived quantities, and the package / TMB-DLL versions used to
produce the fit are printed.

## Usage

``` r
# S3 method for class 'Rceattle'
print(x, ...)
```

## Arguments

- x:

  An object of class `"Rceattle"` returned by
  [`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md).

- ...:

  Currently unused.

## Value

The input `x`, invisibly.

## Details

For operational use, the package version line is meant to make it
obvious which version of `Rceattle` produced an archived fit so that
results can be reproduced even if `main` has moved on. Tag a release
(`devtools::install_github("afsc-assessments/Rceattle@X.Y.Z")`, bare
version, no leading `v`) and the same version string will reappear here
on a fresh run.
