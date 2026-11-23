# 5. Single- vs. multi-species models

When to fit a single-species or a multispecies (predation) model, and
how to switch between them, with Bering Sea and Gulf of Alaska examples.
The switch is `msmMode` in
[`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md):

| `msmMode` | Predation mortality (M2)            |
|-----------|-------------------------------------|
| `0`       | Off, single-species; total M = M1   |
| `1`       | Type II MSVPA (Holsman et al. 2015) |
| `2`       | Type III MSVPA                      |

In single-species mode the model is a standard statistical catch-at-age
model. In multi-species mode, a bioenergetics-based predation component
(M2) is added to natural mortality so that predator abundance directly
affects prey survival. The two modes share the same data format and most
[`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)
arguments: switch `msmMode`, and lower a fixed M1 where predation now
supplies part of it (see `BS2017MS`).

## When to use each mode

**Single-species** is appropriate when:

- Predation interactions are unknown or negligible
- The species of interest has no major predators or prey in the system
- You are building an initial model or debugging data inputs

**Multi-species** is appropriate when:

- Predation is an important source of mortality (e.g., groundfish in the
  Bering Sea or Gulf of Alaska)
- You want to evaluate how climate-driven changes in predator or prey
  abundance propagate through the food web
- Reference points should account for ecosystem effects

A common workflow is to fit single-species models first and use those
estimates as starting values (`inits`) for the multi-species run.

## Bering Sea example

``` r

library(Rceattle)

# BS2017SS and BS2017MS use the same data structure but differ in M1_base:
# the multi-species dataset has lower residual mortality for pollock ages 1-2,
# where predation (M2) explains part of it in multi-species mode.
data(BS2017SS)
data(BS2017MS)
BS2017SS$projyr <- 2060
BS2017MS$projyr <- 2060

plot_data(BS2017SS)
```

``` r

# Single-species
ss_run <- fit_mod(
  data_list    = BS2017SS,
  inits        = NULL,
  file         = NULL,
  estimateMode = 0,
  random_rec   = FALSE,
  msmMode      = 0, # <-- single-species
  fit_control = fit_control(
    phase        = TRUE,
    verbose      = 1))

# Multi-species, start from single-species MLEs for stability
ms_run <- fit_mod(
  data_list    = BS2017MS,
  inits        = ss_run$estimated_params,  # warm start
  file         = NULL,
  estimateMode = 0,
  niter        = 5,        # iterations between population and predation dynamics
  random_rec   = FALSE,
  msmMode      = 1,        # <-- Type II MSVPA
  suitMode     = 0, # empirical suitability from diet data
  fit_control = fit_control(
    verbose      = 1))
```

The `niter` argument sets how many times the population dynamics and
predation mortality equations are iterated in each evaluation of the
objective; the count is fixed, with no convergence test. Five iterations
is typically sufficient.

## Gulf of Alaska example

`GOA2018SS` is a three-species data list (pollock, arrowtooth flounder,
Pacific cod) for the 2018 GOA assessment, used for both the single- and
multi-species runs.

``` r

data("GOA2018SS")
plot_data(GOA2018SS)

# Single-species
goa_ss <- fit_mod(
  data_list    = GOA2018SS,
  inits        = NULL,
  file         = NULL,
  estimateMode = 0,
  random_rec   = FALSE,
  msmMode      = 0,
  fit_control = fit_control(
    phase        = TRUE,
    verbose      = 1))

# Single-species with estimated M
goa_ss_M <- fit_mod(
  data_list    = GOA2018SS,
  inits        = goa_ss$estimated_params,
  file         = NULL,
  estimateMode = 0,
  M1Fun        = build_M1(M1_model = c(1, 2, 1)),  # Estimate M: single M for spp 1 & 3, sex-specific for spp 2
  random_rec   = FALSE,
  msmMode      = 0,
  fit_control = fit_control(
    phase        = TRUE,
    verbose      = 1))

# Multi-species
goa_ms <- fit_mod(
  data_list    = GOA2018SS,
  inits        = goa_ss$estimated_params,
  file         = NULL,
  estimateMode = 0,
  M1Fun        = build_M1(M1_model = c(1, 2, 1)),
  niter        = 3,
  random_rec   = FALSE,
  msmMode      = 1,
  suitMode     = 0,
  fit_control = fit_control(
    phase        = TRUE,
    verbose      = 1))
```

## Predator-prey suitability (`suitMode`)

`suitMode` sets how prey preference is estimated for each predator
species: one value for all predators, or a vector of length `nspp`.

| `suitMode` | String | Description | Available? |
|----|----|----|:--:|
| `0` | `"Empirical"` | Empirical proportions from observed stomach content data (MSVPA-style; Holsman et al. 2015) | Yes |
| `1` | `"GammaLength"` | Gamma suitability on prey length-at-age | Not yet |
| `2` | `"GammaWeight"` | Gamma suitability on prey weight-at-age | Yes |
| `3` | `"LognormalLength"` | Log-normal suitability on prey length-at-age | Not yet |
| `4` | `"LognormalWeight"` | Log-normal suitability on prey weight-at-age | Yes |
| `5` | `"NormalLength"` | Normal suitability on prey length-at-age | Not yet |
| `6` | `"NormalWeight"` | Normal suitability on prey weight-at-age | Yes |

The length-based modes are blocked pending validation of the
growth-model integration:
[`data_check()`](https://afsc-assessments.github.io/Rceattle/reference/data_check.md)
errors with “Length-based suitability (suitMode 1, 3, or 5) is not yet
implemented”.

## Comparing outputs

``` r

mod_list  <- list(ss_run, ms_run)
mod_names <- c("Single-species", "Multi-species")

# Population trajectories
plot_biomass(Rceattle = mod_list, model_names = mod_names)
plot_ssb(Rceattle = mod_list, model_names = mod_names)
plot_recruitment(Rceattle = mod_list, model_names = mod_names, add_ci = TRUE)
plot_ssb_depletion(Rceattle = mod_list, model_names = mod_names)

# Mortality breakdown (M1 vs M2), these accept model lists
plot_m_at_age(Rceattle = mod_list, model_names = mod_names, age = 2)
plot_m2_at_age_prop(Rceattle = mod_list, model_names = mod_names)

# plot_mortality() takes one model and draws one component at a time, never the
# sum; plot_m_at_age() above gives total M (M1 + M2).
plot_mortality(ms_run)                    # M2-at-age, one line per year
plot_mortality(ms_run, M2 = FALSE)        # M1-at-age instead
plot_mortality(ms_run, type = "heatmap")  # the same series as an age x year tile

# Predation (multi-species only), these accept model lists
plot_b_eaten(Rceattle = mod_list, model_names = mod_names)
plot_b_eaten_prop(Rceattle = mod_list, model_names = mod_names)
```

For model averaging across single- and multi-species variants, see
[`?model_average`](https://afsc-assessments.github.io/Rceattle/reference/model_average.md).
