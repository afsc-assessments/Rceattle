# Model display names for a list of fits

Returns `model_names` if supplied, otherwise `"Model 1"`, `"Model 2"`,
... so the colour/legend mapping always has labels.

## Usage

``` r
.model_labels(model_list, model_names = NULL)
```

## Arguments

- model_list:

  List of fits.

- model_names:

  Labels, or `NULL` for positional ones.

## Details

`model_names` is often built as a
[`list()`](https://rdrr.io/r/base/list.html), the package's own
vignettes do , so it is flattened to character here. Left as a list it
becomes a one-element list per model and the plot frame fails to bind.

Too few names would be recycled, drawing two models as one series under
one legend key, so that is a warning rather than a silent merge.
