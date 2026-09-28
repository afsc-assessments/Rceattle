# Okabe-Ito colorblind-safe qualitative palette

Eight fixed hues in a fixed order, reordered so the strongest-contrast,
best-separated hues come first (2-4 model overlays are the common case)
and the low-contrast yellow comes last. Unlike a sequential ramp
(viridis), this encodes *identity* (which model), not magnitude: the
colours hold no implied ordering and every hue holds a legible contrast
against a white panel.

## Usage

``` r
.okabe_ito
```

## Details

Reference: Okabe & Ito (2008), "Color Universal Design".
