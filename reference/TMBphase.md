# Run TMB with ADMB-style phased estimation

Fits a TMB model in phases, like ADMB: parameters are switched on in
stages rather than all at once, which stabilizes a difficult
optimization. The optimizer (nlminb by default) runs once per phase; at
each phase the map holds fixed any parameter whose phase is greater than
the current phase, so parameters turn on progressively as the phase
counter advances. `phases` is a tagged list giving each named parameter
its integer phase; parameters not listed are estimated from the first
phase onward.

## Usage

``` r
TMBphase(
  data,
  parameters,
  map,
  random,
  phases,
  model_name,
  silent,
  use_gradient = TRUE,
  control = list(eval.max = 1e+09, iter.max = 1e+09, trace = 0)
)
```

## Arguments

- data:

  A list to be passed to TMB

- parameters:

  A list of parameters of the model

- map:

  a list of map object from the model

- random:

  A character vector of names of parameters that are random effects

- phases:

  Tagged list assigning each named parameter its integer phase (as
  returned by
  [`set_phases()`](https://afsc-assessments.github.io/Rceattle/reference/set_phases.md)).

- model_name:

  A string describing the model name. Must be the name of your .cpp file

- silent:

  logical. If TRUE, suppresses output from TMB (default = TRUE).

- use_gradient:

  logical. If TRUE, uses gradient in optimization (default = TRUE).

- control:

  A list of control parameters. For details see
  [`?nlminb`](https://rdrr.io/r/stats/nlminb.html)

## Value

The parameter list estimated in the final phase, with a per-phase
convergence log attached as the `phase_log` attribute.

## Author

Gavin Fay
https://github.com/kaskr/TMB_contrib_R/blob/master/TMBphase/R/TMBphase.R
