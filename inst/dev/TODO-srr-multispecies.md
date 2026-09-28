# TODO: stock-recruit curves in multispecies models

Status: **open**, pruned 2026-09-20 after 5.39.0 closed items 3 and 4.
Opened 2026-09-11 with 5.30.0 (`35fbc19b`) and 5.31.0, after
two adversarial reviews of each. **5.39.0 closed items 3 (bounds), 4, 9 (the
harness) and the tests of 12**: `rec_pars[, 2:3]` are bounded at +/-30 (log);
`.check_stock_recruit_msm()` reads the curve over the observed SSB range (flat /
linear ridge, Ricker peak below the data, bound, standard error);
`tools/verify/verify-sim-recovery-srr-msm.R` is the recovery harness. The Ricker
criteria (density-dependence factor above 0.9 at the highest SSB = linear; peak
`1e6 / beta` below the lowest SSB = every observation on the descending limb)
were chosen without a reviewer; revisit if a real Ricker fit trips them wrongly.

**Item 2 was attempted in 5.39.0 and reverted.** Decaying the initial ages with
`M_at_age(..., 0)` (M1 + lagged M2) inside the predation iterations makes year-1
N and year-1 M2 mutually dependent: M2 is consumption over prey numbers, so a
smaller initial N raises M2, which shrinks the next iteration's initial N. On
`BS2017MS` from default starts the year-1 numbers of pollock ages 2-3 fell from
8020 / 5940 to 3e-25 / 3e-28 by `niter = 2` and the objective was `Inf` at
`niter = 10`; the golden `ms` fit (warm-started, `niter = 5`) happened to
converge, to 10231.72 against 10267.25. It also changes the initial structure of
`Equilibrium` and `OffsetEquilibrium`, whose `init_dev` is mapped off, and breaks
the equivalence tests against the R simulator and CEATTLE classic, both of which
decay with M1 only. A stable form needs M2 for the pre-`styr` years that does
not feed back through the initial numbers: an M2 taken from the M1-only initial
structure and held (a one-shot correction), or a switch that reads a
user-supplied pre-`styr` M2. Either is a design decision, not a patch.

Spawning biomass per recruit (SPR) is undefined under predation: total mortality
carries M2, which scales with predator abundance, so the template leaves `SPR0`
and `SPRFinit` at 0 when `msmMode > 0`. Everything below follows from that.

## Where it stands

- **5.30.0**
  - The Ianelli penalty is allowed under predation (`srr_fun = "mean"`, curve as
    `srr_pred_fun`).
  - A Beverton-Holt steepness prior is refused. An alpha prior goes through a
    linkage, which acts on `rec_pars` itself: `normal()` is on natural-scale
    alpha; `prior_lognormal()` is `dnorm(log alpha)`, centred so `exp(meanlog)`
    is alpha's mean under `bias_adjust_proc` (5.33.0) and its median without.
  - Species with input numbers-at-age (`estDynamics > 0`) carry no curve terms.
  - Dynamic B0 keeps those species at their input numbers.
  - `fit_mod()` warns when a supplied map fixes the curve.
- **5.31.0**
  - A curve fitted in the hindcast is allowed under predation, in every
    initMode (see item 1).
  - `R_init` is the free level `exp(rec_pars[, "R0"])`; steepness is 0; the
    Ricker `posfun(alpha * SPR0 - 1)` penalty is skipped.
  - `sample_rec(sample_rec = FALSE)` and `retrospective()` project
    `log(mean(R / R_hat))` = `log(mean(exp(rec_dev)))`, which is
    mean-unbiased.

## Measured on hake (4 species, `MSE_yr2024.R` setup, MSVPA, 2026-09-11)

| Configuration | Objective | Hake alpha, beta | Terminal SSB vs baseline |
|---|---|---|---|
| NonEquilibrium, mean recruitment (baseline) | 2663.81 | — | — |
| Ianelli BH penalty, alpha/beta free | 2732.20 | 18.8, 2.7e-6 | +0.04% |
| FreeParams, mean recruitment (control) | 2520.14 | — | −34.7% |
| FreeParams, BH in the hindcast | 2520.08 | 4.5e8, 57 (flat) | −34.7% |
| NonEquilibrium, BH in the hindcast, free `R_init` | 2654.13 | 5.2e8, 58.5 (flat) | −24.9% |

Measured on 5.31.0. 5.33.0 mean-centred the lognormal priors and the Ianelli
penalty, so these objectives cannot be compared with 5.33.0 fits.

- **The penalty row's objective** excludes the fixed predators' curve penalty
  (5,137.9 nats), which 5.30.0 removed. Objectives are not comparable across
  initModes: FreeParams drops the initial-abundance density (−28.8 nats on
  hake).
- **The hake data do not inform a curve fitted in the hindcast.** Both hindcast
  fits run to the flat ridge, recruitment ≈ alpha/beta ≈ 8e6 at any SSB (the
  predicted/asymptote ratio is 1.000 at median SSB). Their Hessians are
  positive definite but near-singular (condition 1.5e14 and 6.5e13).
- **FreeParams drives the shift, not the curve.** It drops M1 from 0.265 to
  0.190 and mean recruitment by 43%, and moves SSB −56% in 1980 and −35% at the
  end. The 1980 ages 16–20 fall to 259–450 fish, with log-scale standard errors
  of 2.7–4.0. The composition fit improves only 6.7 nats.

## Open items, most important first

1. **Done in 5.31.0: every initMode is allowed.** Under a hindcast curve under
   predation, `R_init = R0` for every initMode (`ceattle.cpp` 6.3), and none
   reads SPR. Still open:
   - Under 3/4, `R_init` is confounded with Finit.
   - Under FreeParams or OffsetEquilibrium with `minage = 1`, `R_init` enters only through
     `R(0) = R_init * exp(dev0)`, so one deviation's density pins it. No fitted, Hessian or
     simulation-recovery test covers the free level.
   - NonEquilibrium converged on the fixture when started from the
     mean-recruitment fit (gradient 0.05), but collapsed from default starts.
     Add a fitted test and document warm-starting.
2. **NonEquilibrium decays the initial ages with M1 only**, so under predation
   the initial deviations absorb the missing M2. Hake's mean `init_dev` is −1.6,
   against age-1 M2 of 0.81 plus σ²/2 of 0.77.
   - **Built and measured during 5.39.0, then REVERTED.** Decaying with
     `M_at_age(..., 0)` moved the golden multispecies fits as expected (`ms` to
     10231.72, `goa_ms` to 12915.74), but year-1 numbers and M2 feed back across
     the predation iterations: `BS2017MS` from its default starts reached an
     infinite objective at `niter = 10`. It also breaks the CEATTLE-classic and
     R-simulator equivalence tests, which decay with M1.
   - So the fix needs a stable formulation, not just the substitution: either
     hold year-1 M2 at the previous iteration's value, or converge the initial
     age structure inside the iteration. Whoever picks this up should start from
     the divergence, not from the substitution.
3. **Anchor the curve — bounds done in 5.39.0.** Under predation nothing ties
   alpha and beta to an equilibrium, so they run to the flat ridge.
   `rec_pars[, 2:3]` now carries a default `[-30, 30]` on the log scale, applied
   before the linkage override so a user bound still wins. It is an overflow
   guard, not an anchor: Ricker's beta is divided by 1e6 in the template, so the
   two forms sit 13.8 log units apart and no tighter default is defensible. A
   prior through a linkage remains the documented way to anchor a curve.
   - On the 2-species fixture log alpha reached 702, next to double-precision
     overflow at 709.8, and a 1-peel retrospective failed with NaN gradients.
     Holding alpha at its single-species value cost 2.5 nats.
   - `rec_pars` has no bounds.
   - Options:
     - bounds on `rec_pars[, 2:3]`;
     - an alpha prior through a linkage as the documented configuration;
     - item 4's diagnostic.
4. **Done in 5.39.0: the degenerate-curve check.** `.check_stock_recruit_msm()`
   reports a flat or linear Beverton-Holt, a linear or descending-limb Ricker, a
   parameter at the 5.39.0 bound, and a standard error above 10, reading the
   hindcast `sdreport` through `fit$.conv_hindcast` — under an estimating HCR
   `fit$sdrep` is the projection's and reports 0 rather than NaN. It ships with
   the bounds deliberately: a bound without it converts a loud ridge into a fit
   that looks converged.
5. **Projected-deviation conventions: done in 5.33.0.** The Ianelli branch now
   projects `log(mean(R / R_hat))` over the penalty years, and its
   penalty is mean-centred under `bias_adjust_proc`. The median form ran hake
   projections at about 45% of the curve's mean. A single-species curve fitted
   in the hindcast now takes `log(mean(exp(rec_dev)))`, as the multispecies one
   does (`test-functions-sample-rec-curve.R`). Still open: a test that expected
   projected R agrees with and without resampling.
6. **Done in 5.33.0: the Ianelli dynamic B0 uses the realized deviation from
   the curve.** Before `srr_mse_switchyr` it takes `log R - log R_hat` rather
   than `rec_dev`, which is measured from R0 there. On the hake operating model
   `R_hat / R0` was about 1.5 in those years. The dynamic runs keep the
   suitability fitted in the hindcast; a full no-fishing refit would re-derive
   empirical suitability from its own abundances, which is why `remove_F()`
   removes fishing only after `suit_endyr`.
7. **Multispecies SB0 (`MSSB0`).**
   - It is derived only when `HCR != "NoFishing"` (`fit_mod.R`); otherwise it is
     the 999 mt placeholder.
   - It comes from a mean-recruitment projection to `projyr`, which is not an
     equilibrium: at M ≈ 0.27, cohorts alive at `endyr` still carry about 33%
     of SSB seven years on.
   - Verify whether HCRs 3/6 can reach the objective with the placeholder
     (reviewer: `ceattle.cpp` ~2136–2183, 4543–4548). The other session's
     `test-switches-hcr-multispecies.R` covers the same area.
   - Options:
     - a converged no-fishing projection;
     - reporting depletion against dynamic B0;
     - an in-model multispecies equilibrium (fixed-point M2), which item 8
       needs.
8. **Steepness under predation.**
   - `h = 0.2(1 + beta S0)/(1 + 0.2 beta S0)` is steepness only when `S0` is
     the curve's own equilibrium.
   - If it is needed, define it from an in-model equilibrium S0 and report it
     as conditional steepness.
   - Alpha (recruits per unit female SSB, the slope at the origin) is what
     carries over from a single-species fit, not h. It must be restated in the
     multispecies model's SSB units and definition.
9. **Simulation recovery.**
   - No test shows the hindcast curve or the Ianelli penalty recovering a known
     alpha/beta. The fixture's SSB spans only 1.7×, too little contrast; build
     one with real contrast.
   - Also measure the σR inflation the Ianelli double density implies under
     `random_rec` (up to √2 by the algebra; `review4/re.R` did not finish).
10. **Per-species `srr_fun` / `srr_pred_fun`**, so hake can carry a curve while
    the predators stay on mean recruitment. Both are global today.
11. **MSE resampling.** `sample_rec()` resamples independent hindcast
    deviations, including year 1 and years outside the `srr_hat` window. A
    block bootstrap within the window would keep autocorrelated recruitment.
12. **Smaller.**
    - `retrospective()`'s predation guard is untested, and the peel refit
      overwrites the column it writes.
    - A year-varying R0 linkage is inert on the hindcast-curve path.
    - Under mean recruitment with `minage > 1`, the hindcast takes `R_init` for
      the first `minage - 1` years (`ceattle.cpp` section 6.5), so a year-varying R0
      linkage has no effect there. Dynamic B0 copies the hindcast's R.
    - `R_init` starts at exp(9) = 8103, and the alpha/beta defaults are
      off-scale for tonnes; document seeding with `srr_alpha_init` /
      `srr_beta_init`.
    - The 5.31.0 tests are build-only (`estimateMode = 3`); add one fitted
      case.
13. **What to report with a multispecies stock-recruit fit:**
    - the SSB-recruitment pairs with the curve;
    - the predicted/asymptote ratio;
    - alpha and beta with their standard errors;
    - sensitivity to initMode;
    - how σR was configured, and the projection convention;
    - the dynamic B0 series;
    - a note that objectives are not comparable across these configurations.

    Document what the Ianelli penalty fits: recruitment centred halfway between
    the mean and the curve, at half the variance.
14. **Recruitment floor under an identity-link offset (5.35.0) is incomplete.**
    An identity-link offset adds to R0, alpha and Beta on the natural scale
    (`ceattle.cpp` 5.6), so it can make them non-positive. 5.35.0 floors hindcast
    R, R_hat, the penalty curve and the derived year-0 R0/R_init (under
    `rec_floor_on`), but not:
    - the equilibrium recursion `NByage0`/`NByageF` (6.6), so SB0, SBF and B0,
      and the depletion reference HCRs 5 and 6 read, can go negative;
    - dynamic-B0 recruitment `N_at_age_dB0`/`dBF` (6.6);
    - projected recruitment (6.8).

    Each spot carries a `TODO` in the template. Flooring dynamic-B0 recruitment
    has moved fits before, so run `/golden-check` and the hake MSE after the change.
    `test-recruitment-curve-floor.R` asserts only `is.finite(DynamicSB0)`, which
    a negative value passes; assert `> 0` on R0, SB0 and DynamicSB0 with the fix.
