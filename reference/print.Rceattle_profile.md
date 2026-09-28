# Print method for a likelihood profile

Reports whether the grid actually brackets the minimum. A profile whose
lowest point is its first or last grid value has not found the optimum,
it has run out of grid, and the curve drawn from it understates how far
the parameter can move. That is the failure the numbers alone hide,
since a partial profile plots as a perfectly ordinary line.

## Usage

``` r
# S3 method for class 'Rceattle_profile'
print(x, cutoff = 1.92, ...)
```

## Arguments

- x:

  A `"Rceattle_profile"` object from
  [`profile.Rceattle()`](https://afsc-assessments.github.io/Rceattle/reference/profile.Rceattle.md).

- cutoff:

  Objective units above the minimum that bound the reported interval.
  Default `1.92`, the 95% profile-likelihood cutoff for one parameter.

- ...:

  Currently unused.

## Value

`x`, invisibly.

## Details

For a one-dimensional profile the interval reported is the usual
profile-likelihood one: the grid values within `cutoff` of the minimum,
where the default 1.92 is \\\chi^2_1(0.95)/2\\. It is read off the grid,
so it is no finer than the spacing of `values`, and it is reported as
open on either side the grid does not close. It is also referenced to
the best GRID point rather than to the unconstrained MLE, which the
object does not hold: the grid minimum sits at or above the MLE, so the
interval errs wide. And it is reported as a range, so a profile with a
second basin, or a failed point inside the range, is called out as not
contiguous rather than left to read as one interval. No interval is
given for a cross-profile over two or more cells, the cutoff would be
\\\chi^2_k(0.95)/2\\ and the region is not an interval.

Under `random_rec = TRUE` the objective is the Laplace-approximated
marginal likelihood, so this is a profile of that, with the usual caveat
that the approximation is what is being profiled.
