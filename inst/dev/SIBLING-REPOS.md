# The sibling assessment repos

The repos below are live consumers of this package's API, and this is the one list of them. A
breaking change here breaks scripts that produce federal catch advice.

- **`../Rceattle-models`** — EBS/GOA pollock, sablefish, arrowtooth, plaice, POP, hake.
- **`../GOA-ATF-ESP`** — GOA arrowtooth and its multispecies (cannibalism) run: **the only live
  two-sex, `suitMode = 0` model**, so it is what exercises the sexed and predation paths.
- **`../Climate_MSE`** — GOA climate-linked multispecies MSE: pollock, arrowtooth and cod, with
  SSP126/245/585 operating models. Brought to the current API on 2026-09-11 but **not yet refit**;
  see its section below.
- **`../GOA-multispecies-assessment`** — the GOA multispecies assessment; `run_all.R` and
  `R/02_fit_models.R` call `Rceattle::` directly.
- **Ignore `EBS_CEATTLE_TMB`** — a vendored fork, not a consumer.

Fitted `*.rds` are ~50 MB each. Keep them out of git.

## Sweeping

```
grep -rn "<symbol>" --include=*.R "../Rceattle-models" "../GOA-ATF-ESP" "../Climate_MSE" "../GOA-multispecies-assessment"
```

### `Climate_MSE`

Entry point `R/Climate_MSE_GOA_runs.R`: it sources the OM and EM conditioning scripts, then
`run_climate_mse()`. Both conditioning scripts start their fits from
`Models/GOA_20_1_1_mod_list.RData`, which `Models/GOA_23.1.1. fit models.R` writes.

It will not run until the saved 2024 fits are regenerated. They predate the current parameter
set, so `fit_mod()` stops on them as `inits`; rerun the fit script first.

Two data checks that 2024 Rceattle did not have also refused its workbook, and the scripts now
handle both after `read_data()`. `Pcod_spawn_srv` and `Pcod_seine_srv` estimated selectivity with
no composition data, so they are turned off. 252 of 4,096 `diet_data` rows carried cod at ages
11–12 against a cod model of ages 1–10; Climate_MSE's own helper `fold_diet_plus_group()`
(`R/Functions/` there, sourced per script -- it is not an Rceattle export) folds them into age 10.

The port had to catch three silent changes. Any 2024-era script carries the same risk:

- **`initMode = 1` meant unfished equilibrium *with* initial deviates in 2024. That is now `2`.**
  Leaving it at `1` drops the deviates without an error.
- **`srr_fun = 1|3|5` with `srr_indices` was inert from 4.4.0 through 5.31.0, and is an error from
  5.32.0.** Commit `862ad197` removed the term, although NEWS 4.4.0 says both still work. The objective and parameter count are
  identical to the non-environmental model. Climate-driven recruitment is now a linkage on `R0`
  (mean recruitment) or `alpha` (Ricker). The old `srr_env_indices` counted `env_data` columns
  *after* `Year`, so Climate_MSE's `c(2,3,4)` meant winter SST, SST squared and zooplankton. It
  did not include bottom temperature.
- **Projected recruitment ignores an `R0` linkage under the default `proj_mean_rec = TRUE`.** It is
  the hindcast mean, where 2024 multiplied that mean by the environmental term. The ported
  mean-recruitment climate OMs therefore lose the effect in every projected year. `run_mse()`
  moves assessed years into the OM's hindcast, where the effect applies. The Ricker OMs set
  `proj_mean_rec = FALSE` and keep it. Restoring the 2024 behaviour was declined (2026-09-11), so
  treat those OMs' projected years as climate-naive recruitment.

Two limits, and both need saying out loud when you report a sweep:

- **It catches removed or renamed API. It does not catch behavioural drift.** A script that
  still parses can still produce different numbers.
- **Some models there are partially implemented and do not run regardless**, so not every hit
  needs chasing.

The v4.7.0 ggplot migration is the cautionary case: `plot_*()` began returning ggplot objects,
so every `plot_x(...); mtext(...)` chain in those scripts failed with "plot.new has not been
called yet". A grep for the function names would have found them; nobody ran one.

**Sweep the workbooks, not only the scripts.** Three of the four behaviour changes in 5.25.0
are settings a workbook holds and a script never names, so a grep over `*.R` sees none of
them. `readxl::read_excel()` over every `.xlsx` in the four repos is the net: 374 workbooks,
183 with a `fleet_control` sheet, about two minutes. It is what turned "183 sheets have
`Accumatation_age_*`" into the true 147, and what showed that the pre-`styr` `env_data` row
shift misses the hake operating model. Note the `control` sheet is TRANSPOSED -- the first
column holds the element names and each species is a column -- so `"nsex" %in% names(control)`
is always FALSE and reads as a missing `nsex` on every workbook in the ecosystem.

## Running them, when a sweep is not enough

Neither repo caches a fitted object, so verifying against a real assessment means refitting.
That is cheaper than it sounds: the terminal fit is under a minute for either pollock model,
~13 min for the ATF chain.

| Model | Entry point |
|---|---|
| GOA pollock 2025 | `../Rceattle-models/GOA pollock/2025/04-fit-and-diagnostics.R` |
| EBS pollock 2024 | `../Rceattle-models/EBS pollock/2024/04-fit-and-diagnostics.R` |
| GOA arrowtooth | `../GOA-ATF-ESP/R/2026 assessment w HCR projection.R` (2025: `R/Run_2025_ceattle.R`) |
| Pacific hake MSE | `../Rceattle-models/Pacific hake/04-mse.R` |
| Pacific hake MSE, 2024 | `../Rceattle-models/Pacific hake/MSE_yr2024.R` |

**The hake MSE is the one script that runs `run_mse()` end to end**, and the only routine
exercise of three-species predation with estimated suitability, of `suitMode` differing per
predator, and of Dirichlet-multinomial comps with a prior on their own weight. Golden
covers none of that: its four models are single- and multi-species Bering Sea and Gulf of
Alaska hindcasts, with no MSE and no estimated suitability. Run it after touching predation,
suitability, the DM likelihood, `sim_mod()`, or `run_mse()`.

Its four fits take ~3.5 min together on an M-series Mac, plus ~2 min for `nsim = 2, cores = 2`.
Reference objectives on 5.33.0 (2026-09-11), after its lognormal priors became mean-centred
under `bias_adjust_proc`. The 5.32.1 values are in the right-hand column; the change comes from the DM
`prior_lognormal(0, 2)` weights and the M prior. Survey DM theta fell from 35 to 25
(single-species) and 32 to 23 (MSVPA); hake terminal SSB changed by at most 0.71%. Every fit kept a
positive-definite Hessian. The 5.32.1 column is measured on dev `cff500c7`.

| Stage | -log L (5.33.0) | 5.32.1 |
|---|---|---|
| single-species | 2136.8588547522 | 2133.8207228717 |
| single-species + category-1 HCR | 2137.5094597505 | 2134.4713944220 |
| MSVPA, estimated M | 2140.4295989555 | 2137.4433306648 |
| estimated suitability | 2267.4725502601 | 2260.7063099168 |

**Those four are STALE against the script as it stands.** `MSE_yr2024.R` changed on 2026-09-15
(`bc9596e`, Rceattle-models): it added lognormal intercept priors and an M prior to the existing
fits and appended an Ianelli Beverton-Holt section, so the script now fits different models than
it did when the table above was recorded on 2026-09-11. No input data changed. Re-measured
2026-09-21 on the current script, running it twice against different package versions:

| Fit (in script order) | 5.41.0 | 5.33.0 |
|---|---|---|
| single-species | 2440.0942 | 2440.0942 |
| single-species + category-1 HCR | 2440.6633 | 2440.6633 |
| MSVPA, estimated M | 2447.0049 | 2447.0049 |
| estimated suitability | 2669.3776 | 2669.3776 |

All six fits the script reports, not just these four, agree to every digit printed, and neither
run errored. **So 5.34.0-5.41.0 moved nothing here**; the gap against the older table is the
script, and a reader comparing to it would see a 300-nat regression that does not exist.
Re-record against whichever script revision you ran, and say which.

Re-run on 5.25.0 (2026-09-01), against that day's references (stage 2 2134.4713926593, stage 4
2260.7063099135): stages 1, 3 and 4 bit-identical, and stage 2 higher
by 1.8e-06 (8.3e-10 relative, below the optimizer's own tolerance). Stage 2 is the only one that
runs the reference-point penalty, so it is the only one that touches the SPR sum, whose factors
5.24.1 reordered -- floating-point addition is not associative, so the last bits move and the F
solve propagates it into the objective. **A delta of that size on stage 2 alone is the expected
result of an SPR change, not a numeric regression.** Any of the other three moving is.

All four were bit-identical across the two branches (delta 0.000e+00), as was the vulnerability
matrix: 0.8172 (arrowtooth to hake) and 0.7686 (sablefish to hake). **Record which quantity that
is next time it is measured** -- `exp(log_phi)` on the stage-4 fit prints 4.4709 and 3.3213, so
0.8172/0.7686 is some other transform and the note does not say which. The objective above is
bit-identical, so the parameters are unchanged either way. The script's own
inline comments give the first three ~5 higher and the fourth as 2262.318: those are Rceattle
5.6.1 numbers that still included the `theta_diet` prior constants. `README.md` in that folder
records the "clean" values, which are what the current package reproduces.

### `MSE_yr2024.R` — the four-species MSE

`MSE_yr2024.R` is the newer run and the one to check first: California sea lion, sablefish,
arrowtooth and hake from `MSE_hake_yr24_final.xlsx`, `endyr` 2023, projected to 2030, with
Dirichlet-multinomial age and diet composition and a lognormal prior on every DM weight. Seven
fits plus `run_mse(nsim = 2, cores = 2)`, about 13 min.

Baseline objectives on 5.33.0 (2026-09-11), after the lognormal DM and M priors and the
Ianelli penalty became mean-centred; 5.32.1 on the right. Hake terminal SSB changed by at most
0.86%, survey DM theta fell from 43-49 to 31-35 and diet theta from 94.2 to 56.5 and 18.4 to 11.8 (the M-estimated fits), and every fit kept a
positive-definite Hessian. The script ran end to end, `run_mse()` included.

| Fit | -log L (5.33.0) | 5.32.1 |
|---|---|---|
| `ss_run_DM_CSL` | 2440.0941615088 | 2436.8886423469 |
| `ss_run_DM_hcr_CSL` | 2440.6633379056 | 2437.4573180484 |
| `ms_run_DM_CSL` | 2447.0048917469 | 2443.8538668697 |
| `run_ms_CSL_Mest_prior_DM_CSL` | 2669.3775502006 | 2663.8053181169 |
| `run_ms_CSL_Mest_prior_DM_CSL_BH` | 2735.6792676724 | n/a (see below) |
| `run_ms_CSL_Mest_prior_DM_CSL_stable` | 2669.3775502006 | 2663.8053181169 |
| `ss_run_DM_hcr_B0` | 2440.0941615088 | 2436.8886423469 |

The `_BH` row is re-recorded on 5.33.0 (2026-09-12) for section 5 of the script as it now stands:
a `prior_lognormal(log(100), 0.05)` on hake alpha, chosen for strong density dependence. The
earlier 2737.74 and the 5.32.1 value came from drafts of that section without this prior. The fit
has a positive-definite Hessian and max gradient 5.9e-4; hake alpha is 99.87 and beta 7.63e-6, so
1/beta (1.3e5 t) lies well below observed SSB (0.80-3.29 million t) and R/R_max is 0.86 / 0.93 /
0.96 at the minimum / median / maximum. With an SD of 0.05 the prior, not the data, sets alpha.
`run_mse()` takes the `_stable` fit as its OM, not this one, so the prior leaves the MSE unchanged.

Two of those equalities are structural, not coincidences: the `_stable` refit starts from its
parent's `data_list` and returns to the same optimum, and `ss_run_DM_hcr_B0` matches
`ss_run_DM_CSL` because a `ConstantF` HCR never re-optimizes the projection, so `fit$opt` stays
the hindcast's. `ss_run_DM_hcr_CSL` differs because its HCR does estimate.

Its linkage table is **composition-only** — 6 rows, all `process = "comp"`, none carrying a sex
stratum — so a selectivity or per-sex linkage change cannot reach it.

Traps:

- **`run_mse(cores > 1)` works under `pkgload::load_all()`** because `.parallel_lapply()` forks;
  a PSOCK fallback would not see the loaded tree. Do not assume a parallel MSE failure is the
  model.
- `mse_summary()` returns a ragged list (`species`, `fleet`, `total`, `meta`), not a data frame.
  `as.data.frame()` on it errors -- that is the caller's bug, not a broken MSE.

- The pollock scripts' `Data/` paths are relative to the **project** root, not the year folder.
- **`../GOA-multispecies-assessment` runs entirely off saved fits.** `R/03`–`R/07` load
  `models/GOA_26_mod_list.RData` rather than refitting, so a breaking change surfaces there as a
  script error, not as a moved number, and the saved objects can lag the package. Fits are
  ~1 min each; `run_all.R` rebuilds them from `R/02`. Its `R/07_figures_tables.R` renders every
  figure the chapter uses, so it is a cheap end-to-end check of the `plot_*()` surface: source it
  and confirm the manifest it prints reports no new skips.
- **The ATF entry point is small and the sourcing caveat below is stale.**
  `2026 assessment w HCR projection.R` is 68 lines, so the three unassigned-object line numbers
  this note used to carry (`:364`, `:480`, `:570`) cannot refer to it; in the 2025 fallback
  `Run_2025_ceattle.R` the objects at those lines are all assigned earlier in the same file.
  Whatever revision that described is gone. Re-derive against the script you actually run before
  repeating the claim.
- Its `file =` arguments write **into the assessment repo**. Run it from a sandbox that
  symlinks `Data/`.
- Force plots through `ggplot2::ggplot_build()`. A figure that assembles but cannot render is
  not a pass.
