# Gulf of Alaska 2018 SAFE report reference values

A list containing biomass, spawning stock biomass, and recruitment
reference values from the 2018 Gulf of Alaska Stock Assessment and
Fishery Evaluation (SAFE) report, plus the 2020 Pacific cod assessment.

## Usage

``` r
GOAsafe2018
```

## Format

A list with components, each a data frame of `Year` and one column per
assessment (`2018_Pollock`, `2018_Cod`, `2018_ATF`, `2020_Cod`):

- biomass:

  Total biomass time series

- ssb:

  Spawning stock biomass time series

- recruitment:

  Recruitment time series

## Source

The 2018 columns are from the chapters cited in
[`GOApollock`](https://afsc-assessments.github.io/Rceattle/reference/GOApollock.md),
[`GOAcod`](https://afsc-assessments.github.io/Rceattle/reference/GOAcod.md)
and
[`GOAatf`](https://afsc-assessments.github.io/Rceattle/reference/GOAatf.md);
`2020_Cod` is from Barbeaux, S., Ferriss, B., Palsson, W., Shotwell, K.,
Spies, I., Wang, M., and Zador, S. 2020. Assessment of the Pacific cod
stock in the Gulf of Alaska. In: Stock assessment and fishery evaluation
report for the groundfish resources of the Gulf of Alaska. North Pacific
Fishery Management Council, Anchorage, AK.
