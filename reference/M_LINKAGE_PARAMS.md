# Allowed M-parameter names for `linkages` in [`build_M1()`](https://afsc-assessments.github.io/Rceattle/reference/build_M1.md)

Natural-scale names of the underlying natural-mortality parameters that
the linkage system can address. Currently just `M1`, with the default
log link the offset is added to M1 on the log scale (applied across all
ages unless the linkage row pins a specific `age_bin`).

## Usage

``` r
M_LINKAGE_PARAMS
```

## Details

Note: predation mortality `M2` is a derived quantity in CEATTLE (a
function of predator abundance, suitability, and ration), not a
parameter. There is no `M2` linkage target; environmental effects on
predation are mediated upstream via recruitment, growth, suitability, or
ration inputs.
