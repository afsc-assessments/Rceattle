# Prior mean/SD and labels for the estimated elements of a deviation parameter

Maps the estimated (non-fixed) elements of `par_name`, in the
column-major order used by TMB and the covariance matrices, to their
species/fleet, year/age labels, and the process prior mean and SD.

## Usage

``` r
.process_prior_spec(fit, process, par_name, n)
```
