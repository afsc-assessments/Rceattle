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

**A workbook sweep is only a result next to a BASELINE.** Reading every consumer workbook through
`read_data()` + `data_check()` on the release tree reports **258 of 380 failing** -- old data sets,
dropped columns, `Fleet_code` mismatches -- and essentially all of it predates the release. Run the
same sweep on `main` and diff, or the number says nothing. At 5.48.0 (2026-10-02): main 122 ok /
258 error, dev 122 ok / 258 error, **0 newly broken, 0 newly passing**. One workbook
(`BSAI pop/Data/bsai_pop_single_species_2024.xlsx`) changed its message without changing its
verdict -- it already failed on `main` for an unrelated reason and dev additionally reports the new
ageing-error validation. The sweep needs no compiled DLL, since both functions are pure R, so a
`main` worktree with any `.so` dropped in will do.

**Re-run for the 5.49.0 review** (2026-10-03), `dev` `527a4919` against the review tree, over
every `.xlsx` under the four consumer repos: **0 newly broken, 0 newly passing, and not one
workbook's verdict TEXT changed** on 375 files. That is the gate that matters for a `data_check()`
change, and 5.49.0 adds two refusals and loosens two. Note the absolute counts from that run (17
ok / 358 error) are **not** comparable to the 122 / 258 above: this sweep globbed every `.xlsx` in
the trees, including report and output files with no `control` sheet, whereas the earlier figure
came from a script that is not committed. Compare a sweep only with another sweep run the same
way, which is the whole point of the paragraph above.

**The cheapest real API check is the NAMESPACE and the formals**, not a grep: diff `export(...)`
between the two trees, then parse both trees' `R/` and compare each exported function's formal
argument names. At 5.48.0 that was 90 exports either side with none removed or renamed, and the
only signature change was `build_growth` GAINING `sd_form`, `plus_group_length`,
`plus_group_decay`, `pop_lengths`. Nothing a consumer calls could break.

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

**Two things about running these.** Each script opens with `library(Rceattle)`, which attaches the
INSTALLED package -- on this machine 15 feature versions stale -- so load the release tree first
(`pkgload::load_all(<tree>)`) or the run silently verifies the wrong code. And the pollock scripts
resolve their relative paths from the MODEL root, where the `.Rproj` sits, not from the year
subdirectory the table names: run them with the working directory at `GOA pollock/`, not
`GOA pollock/2025/`.

**Which live assessments the 5.48.0 length-comp ageing-error fix reaches** (a non-identity
`age_error` AND length comps on the same species), measured 2026-10-02 by fitting each under
`main` 5.45.3 and `dev` 5.48.0:

| Model | Reached? | Objective | Terminal SSB |
|---|---|---|---|
| GOA arrowtooth 2026, `mod_26_0` (M fixed) | yes, all 31 length rows joint-sex | 38.6514 -> 38.0811 | 429260 -> 429746 mt (**+0.113%**) |
| GOA arrowtooth 2026, `mod_26_1` (M estimated) | same | 0.1475 -> -0.3021 | 331426 -> 330639 mt (**-0.238%**) |
| GOA pollock 2025 | yes -- `max|A-I|` 0.3550, 74 length rows, all COMBINED-sex | unmeasured, the script cannot run (below) | |
| EBS pollock 2024 | no -- identity matrix, 0 length rows | 713.6765 unchanged in kind | |
| Pacific hake (MSE_yr2024) | no | all four references reproduce exactly | |

Every objective that moves, falls. **Note the SSB sign is not universal** -- the same assessment
moves +0.11% with M fixed and -0.24% with M estimated, so do not generalise a direction from one
model. The bundled-data figures in `NEWS.md` (`GOAatf` +0.223%) are a different, older workbook
than the live 2026 one and are not a substitute for it. A reminder that bundled `GOApollock` has
no length comps while the live 2025 pollock workbook has 74: the bundled blast radius is not the
live one.

**`GOA pollock 2025` cannot run as committed** (checked 2026-10-02). Its line 9 loads
`Data/2024pollock_mfix_estSigR.Rdata` for the parameter skeleton, and that file was DELETED from
`Rceattle-models` in commit `9d0d5c1` ("pollock profile M") -- it is not in HEAD, although
`GOA pollock/.gitignore:15` still un-ignores it by name. So this entry point verifies nothing for
anyone until the file is restored there or the script is pointed elsewhere.

**The hake MSE is the one script that runs `run_mse()` end to end**, and the only routine
exercise of three-species predation with estimated suitability, of `suitMode` differing per
predator, and of Dirichlet-multinomial comps with a prior on their own weight. Golden
covers none of that: its four models are single- and multi-species Bering Sea and Gulf of
Alaska hindcasts, with no MSE and no estimated suitability. Run it after touching predation,
suitability, the DM likelihood, `sim_mod()`, or `run_mse()`.

**The two scripts have separate reference tables, and they are 300 nats apart.** `04-mse.R` is
three species and lands near 2136-2267; `MSE_yr2024.R` is four and lands near 2440-2669. Read the
table under the heading for the script you ran. Mixing them shows a 300-nat regression that does
not exist.

### `04-mse.R` — the three-species MSE

Hake, arrowtooth and sablefish. Unchanged since 2026-08-17 (`72bd887`, Rceattle-models), so the
table below is live rather than historical. Four fits plus `run_mse(nsim = 2, cores = 2)` take
~3.75 min together on an M-series Mac.

Reference objectives on 5.33.0 (2026-09-11), after its lognormal priors became mean-centred
under `bias_adjust_proc`. The 5.32.1 values are in the right-hand column; the change comes from the DM
`prior_lognormal(0, 2)` weights and the M prior. Survey DM theta fell from 35 to 25
(single-species) and 32 to 23 (MSVPA); hake terminal SSB changed by at most 0.71%. Every fit kept a
positive-definite Hessian. The 5.32.1 column is measured on dev `cff500c7`.

| Stage | -log L (5.33.0) | 5.48.0 | 5.32.1 |
|---|---|---|---|
| single-species | 2136.8588547522 | 2136.8588547554 | 2133.8207228717 |
| single-species + category-1 HCR | 2137.5094597505 | 2137.5094602853 | 2134.4713944220 |
| MSVPA, estimated M | 2140.4295989555 | 2140.4295989804 | 2137.4433306648 |
| estimated suitability | 2267.4725502601 | 2267.4725502653 | 2260.7063099168 |

**Confirmed at 5.48.0** (2026-10-02): all four reproduce the 5.33.0 references to between 3.2e-09
and 5.3e-07 absolute, 2.5e-10 relative at worst, and `run_mse()` ran end to end for both sims.
**Re-run at 5.49.0** after the release review: every one of the four reproduces the 5.48.0 column
to all 16 digits, as do both vulnerabilities, so none of that pass's fixes reaches this model --
which is the point of running it, since hake is the only routine exercise of three-species
predation with estimated suitability and the only model here with `estDynamics > 0`, the
configuration 5.49.0's new equilibrium-catch refusal names.
5.48.0 is the release that owed this check, because it added a `SIMULATE` draw for the initial
equilibrium catch and a `sim_mod()` write-back for it. Hake also confirms the length-comp
ageing-error fix does not reach it, and that the new equilibrium-catch validation only reports
here: `catch_data` carries rows at 1979 against a 1980 `styr`, and under
`initMode = "NonEquilibrium"` they are dropped with a message, exactly as they were before.

Re-run on 5.25.0 (2026-09-01), against that day's references (stage 2 2134.4713926593, stage 4
2260.7063099135): stages 1, 3 and 4 bit-identical, and stage 2 higher
by 1.8e-06 (8.3e-10 relative, below the optimizer's own tolerance). Stage 2 is the only one that
runs the reference-point penalty, so it is the only one that touches the SPR sum, whose factors
5.24.1 reordered -- floating-point addition is not associative, so the last bits move and the F
solve propagates it into the objective. **A delta of that size on stage 2 alone is the expected
result of an SPR change, not a numeric regression.** Any of the other three moving is. Stage 2
still carries the largest delta at 5.48.0, for the same reason.

**The vulnerability quantity is `fit$quantities$vulnerability[predator, prey]`**, a 3x3 matrix,
which answers the question this note used to leave open. On the stage-4 fit at 5.48.0 it is 0.8189
(arrowtooth to hake, `[2, 1]`) and 0.7717 (sablefish to hake, `[3, 1]`); `exp(log_phi)` on the same
fit prints 4.5215 and 3.3794, so it is not that transform. The 0.8172 / 0.7686 recorded on the
5.25.0 re-run is **not** a comparison for these: it was measured against the 5.32.1-era objective
(stage 4 2260.7063), and the mean-centring in 5.33.0 moved that fit by 6.8 nats, so a different
optimum is the expected result rather than drift. The script's own inline comments give the first
three ~5 higher and the fourth as 2262.318: those are Rceattle 5.6.1 numbers that still included
the `theta_diet` prior constants. `README.md` in that folder records the "clean" values, which are
what the current package reproduces.

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
