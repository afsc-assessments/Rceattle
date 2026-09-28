# Print method for a jitter analysis

Reports what the run was for: how many random starts reached the best
optimum found. The objective values alone cannot say that, non-converged
starts are dropped before the result is returned, so the count of
returned fits is not the count attempted.

## Usage

``` r
# S3 method for class 'Rceattle_jitter'
print(x, tol = 0.01, ...)
```

## Arguments

- x:

  A `"Rceattle_jitter"` object from
  [`jitter()`](https://afsc-assessments.github.io/Rceattle/reference/jitter.md).

- tol:

  Objective units within which a start counts as reaching the same
  optimum. Default `0.01`, which is far below any difference that would
  change management advice and far above optimizer noise.

- ...:

  Currently unused.

## Value

`x`, invisibly.
