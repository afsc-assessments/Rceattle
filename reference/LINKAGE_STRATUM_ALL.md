# Sentinel meaning "applies to every level of this stratum"

Stored where the user passed `NA` for `species`, `sex`, `age_bin` or
`fleet`; the offsets loop over every level, but the intercept-prior loop
reads 0 as level 1 (which `.drop_off_fleet_sel_priors()` relies on).

## Usage

``` r
LINKAGE_STRATUM_ALL
```
