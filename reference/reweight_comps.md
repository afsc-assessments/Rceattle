# Iteratively reweight composition data (McAllister-Ianelli)

Tunes multinomial composition weights by refitting: each pass sets a
fleet's weight to the McAllister & Ianelli (1997) weight implied by the
previous fit, and stops once the weights settle. The weight multiplies
the input sample size, so changing it changes the fit and hence the next
implied weight, which is why the weights are tuned iteratively rather
than in one step.

## Usage

``` r
reweight_comps(
  object = NULL,
  n_iter = 10,
  tol = 0.01,
  fleets = NULL,
  verbose = TRUE,
  fit = NULL
)
```

## Arguments

- object:

  An `Rceattle` model fitted at `estimateMode` `"Estimate"` (0) or
  `"Hindcast"` (1), the modes that optimize the hindcast.

- n_iter:

  Maximum number of iterations (default 10).

- tol:

  Relative change in the weights below which to stop (default 0.01).

- fleets:

  Fleets to tune, as `Fleet_code` values or as names matching
  `fleet_control$Fleet_name` (as in
  [`linkage_spec()`](https://afsc-assessments.github.io/Rceattle/reference/linkage_spec.md)).
  Give ids or names, not a mix, R coerces `c(1, "Survey")` to
  `c("1", "Survey")`. Default `NULL` tunes every eligible fleet.

- verbose:

  Print the weights each iteration (default `TRUE`).

- fit:

  deprecated name for `object`, still accepted so existing scripts keep
  working. Supplying both is an error.

## Value

The model refitted with the final weights, holding a `reweight` element:
`history` (one row per fleet per iteration), `iterations`, and
`converged`. The last rows of `history` are the weights the returned
model was fitted with.

## Details

Only multinomial fleets are tuned. A Dirichlet-multinomial fleet
estimates its own weight within the likelihood, so tuning it externally
would compete with that estimate; those fleets are named and skipped, as
are any requested through `fleets` that have no fitted composition data.

Only the age/length composition weight (`Comp_weights`) is tuned. The
conditional age-at-length and diet weights have McAllister-Ianelli
analogues of their own (`CAAL_weights`, `Diet_comp_weights`), but they
are not tuned here; set them by hand if needed.

Tuning stops when the largest relative change in any weight falls below
`tol`. Hitting `n_iter` first gives a warning, not an error, the partly
tuned fit is returned either way, and `$reweight$history` shows whether
the weights were still settling.

The weight is a parameter, not a data input, so each pass holds it to
the next through `inits`; `fleet_control$Comp_weights` supplies the
value only when a model is built from scratch, and editing that column
then refitting from an existing fit has no effect
([`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)
warns when it detects this). Each pass writes the new weight to both, so
the returned model's data list rebuilds the tuned model from scratch.

## References

McAllister, M.K. and Ianelli, J.N. (1997) Bayesian stock assessment
using catch-age data and the sampling-importance resampling algorithm.
*Canadian Journal of Fisheries and Aquatic Sciences* 54:284-300.

Francis, R.I.C.C. (2011) Data weighting in statistical fisheries stock
assessment models. *Canadian Journal of Fisheries and Aquatic Sciences*
68:1124-1138. (The main alternative; applied by hand rather than by this
function.)

## Examples

``` r
if (FALSE) { # \dontrun{
fit  <- fit_mod(BS2017SS, estimateMode = "Hindcast")
tuned <- reweight_comps(fit)
tuned$reweight$history      # weight per fleet, per iteration
} # }
```
