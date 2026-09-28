# Rerun with F = 0.

Refits the model with fishing mortality set to 0 from `styr` on, keeping
every other parameter. The projection after `endyr` is always unfished.
[`run_mse()`](https://afsc-assessments.github.io/Rceattle/reference/run_mse.md)
uses it for the no-fishing run (`OM_no_F`) behind the collapse metrics.

## Usage

``` r
remove_F(object = NULL, styr = NULL, Rceattle = NULL)
```

## Arguments

- object:

  A fitted Rceattle model object

- styr:

  First year with F = 0; default `endyr + 1`, which leaves the hindcast
  unchanged.

- Rceattle:

  deprecated name for `object`, still accepted so existing scripts keep
  working. Supplying both is an error.

## Details

`styr` may be any year from the model's own `styr` to `endyr + 1`; the
projection is unfished whatever harvest control rule the model was fit
under. Under predation, empirical suitability (`suitMode = 0`) is
derived from the fitted abundance over each predator's
`suit_styr:suit_endyr`, so `styr` must fall after that window for every
predator with `suitMode = 0` and diet data in it.
