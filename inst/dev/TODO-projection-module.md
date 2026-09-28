# TODO: a native projection module (retire the SPM ADMB code)

Status: **open, scoped, not started.** Scoped 2026-09-09 to 2026-09-11.

The goal is an integrated forecast in the style of Stock Synthesis's forecast
file and WHAM's `project_wham()`: the standard AFSC harvest projections run
inside the Rceattle model, not by exporting to the ADMB Standard Projection
Model (SPM, `afsc-assessments/spmR`). `report_tables()` names the gap in its own
header (`R/6-report_tables.R:1-9`): Guidelines §4.11.3, the seven Tier 1–3
harvest scenarios, is the one required SAFE element Rceattle cannot produce.

## What exists today

### spmR is an ADMB wrapper; its R engine is a stub

spmR 0.3.0 (GitHub only, not on CRAN) exports `runSPM()`, `dat2list()`,
`list2dat()`, `as_spm_result()`, `tier3_scenario_table()`, `plotSPM()`,
`plotSPMx()`. The only real engine is `inst/admb/spm.tpl` (2,560 lines). The
package ships no binary; it needs ADMB 13+ to compile, or a prebuilt `spm.exe`
on Windows.

`runSPM(engine = "rtmb")` is **not** an implementation. On both the installed
build and `main` (`R/spmR.R`), `runSPM_rtmb()` writes `SSB = mean(SSB)`,
`Tot_biom = 2 * mean(SSB)`, `F = 0` and `Catch = ABC = OFL = fixed catch` —
the right columns filled with placeholders. Do not use it as a reference.

`spm.tpl` is the operational definition of the seven scenarios. Write the
native version against it, not against the Guidelines' prose.

### A working, validated bridge lives in GOA-ATF-ESP

`../GOA-ATF-ESP/R/Functions/spm_bridge.R` (589 lines) writes SPM's three input
files from an Rceattle fit, compiles and runs SPM, and builds the §4.2.2
executive table and the §4.11.3 scenario table. It also builds the executive
table a second way from Rceattle's own projection, and the two share no code.
Driver: `../GOA-ATF-ESP/R/2026 assessment projections.R`.

Measured on the 2026 GOA arrowtooth SAFE models
(`../GOA-ATF-ESP/2026/projections/mod_26_*/comparison.csv`):

| 2027 | 26.0 Rceattle | 26.0 SPM | 26.1 Rceattle | 26.1 SPM |
|---|---|---|---|---|
| Female SSB (t) | 874,501 | 874,501 | 709,459 | 709,459 |
| FOFL | 0.1179 | 0.1179 | 0.2891 | 0.2891 |
| max FABC | 0.0994 | 0.0994 | 0.2354 | 0.2354 |
| OFL (t) | 157,261 | 157,261 | 284,071 | 284,071 |
| ABC (t) | 134,083 | 134,083 | 237,731 | 237,731 |
| B100% (t) | 1,071,948 | 1,073,870 | 828,565 | 829,454 |
| Total biomass (t) | 1,437,469 | 1,437,490 | 1,980,285 | 1,980,340 |

SSB, FOFL, FABC, OFL and ABC agree exactly. B100% differs by 0.18% / 0.11%
because the two average recruitment over slightly different year sets. **These
are the acceptance numbers for the native module.**

The bridge's header records four file-format traps that fail silently: SPM
wants *female* recruits; SPM's spawn month is Rceattle's `spawn_month + 1`;
spawning weight is one row but M, maturity, fishery weight, selectivity and N
are two rows each (females then males); units go across in Rceattle's own
(thousands of fish, kg, t) with SPM's `scalars = 1`.

It hard-stops outside what it was written for: `nsex == 2`, `minage == 1`,
exactly one fishery, species index 1, and `spawn_month == 0` on the Rceattle-side
table. It writes `spr_abc`/`spr_msy` as literal `0.4`/`0.35`
(`spm_bridge.R:164-165`) instead of reading `data_list$Ftarget`/`Flimit`, so a
fit built with `build_hcr(Ftarget = 0.45)` would still project at F40%. It has
no `msmMode` guard: under `msmMode > 0` it would average M1 + M2 into a static M
and silently turn a multispecies fit into a single-species projection.

### CEATTLE already projects inside TMB

Section 6.7–6.9 of `src/TMB/ceattle.cpp` projects numbers-at-age inside the AD
tape, with the HCR setting F. Structurally it already does what SPM does:

| SPM | Rceattle | where |
|---|---|---|
| Tier 3 ramp on B/B40% | `HCR = "NPFMC"`, static or dynamic reference points | `ceattle.cpp`, the projection's HCR block |
| F apportioned across gears (`Fratio`) | `proj_F_prop` | `F_flt_age(...) = sel_at_age(...) * proj_F_prop(flt) * proj_F(sp, yr)` |
| Baranov catch, plus group, sex ratio, spawn month | same | the projection's catch and numbers-at-age loop |
| B100%, B40%, B35% | `SB0` | `Calculate Dynamic SB0 and SB at F target` |
| Recruitment from the historical mean or an SRR | `proj_mean_rec` | `if(proj_mean_rec == 1)` in the projection |

Two things it does that SPM cannot: propagate parameter uncertainty into the
projection (`fit_control(projection_uncertainty = TRUE)`, the WHAM approach),
and carry predation mortality forward under `msmMode > 0`.

## The gaps

**Gap 0 — averaging window.** The projection freezes biology at the terminal
year: `sel_at_age(..., nyrs_hind - 1)` in the projection's F (marked
`// FIXME using last year of selectivity`), and `weight_hat(..., nyrs_hind - 1)`
in the SB0/SBF sums. SPM, SS3 and WHAM average a window (WHAM `avg.yrs`,
default the last 5); the bridge averages 5 (`spm_bridge.R:104-110`). Identical
for arrowtooth, whose biology is time-invariant; different on any stock with
time-varying selectivity or weight-at-age. Moves numbers, so `/golden-check`.

**Two sites hardcode the terminal year; three more inherit it.** The projection's F
(`ceattle.cpp:2275`, `FIXME using last year of selectivity`) and the reference-point
schedules (`:1595`, marked at `:1602`, `FIXME: time-vary sel in the forecast`) index
`nyrs_hind - 1` themselves. The two `max_catch_hat` sums and `exploitable_biomass` index the
loop year and only see the terminal curve because `normalize_and_project_selectivity()` copies
it forward into every projection year (`selectivity.hpp:120`). **That fill is the single
lever**: an averaging window placed there covers those three automatically, and the two
hardcoded sites still have to be changed by hand -- a window that reaches only the projection's
F leaves the reference points on a different basis from the projection they are compared
against.

**Gap 1 — catch-conditioned projection years (the keystone).** There is no way
to specify catch in a projection year; `proj_F` comes from the HCR only. SS3
and WHAM solve F from a specified catch inside the tape (WHAM `proj.catch`, and
`proj_F_opt` to set each year independently). The bridge works around it with
`uniroot()` outside the model (`spm_bridge.R:525`), which propagates no
uncertainty and duplicates the Baranov equation.

Implement as a fixed-iteration Newton solve on Baranov per species-year (fixed
count, so the tape stays static), apportioned across fleets by `proj_F_prop`,
capped at attainable catch (`max_catch_hat` already models the ceiling). Needed
in practice: GOA arrowtooth projects on 12,408.7 t a year — the ABC in force
times the five-year yield ratio — because removing the full ABC understates the
second-year population.

**Gap 2 — the seven scenarios.** Six are orchestration over existing HCRs:
max permissible ABC (`NPFMC`), author's F (`ConstantF`, or `Fmult`), average
recent F (`ConstantF`; window `endyr-5 .. endyr-1` per §4.11.3, see
`spm_bridge.R:112-117`), alternative SPR (`ConstantFSPR`), no fishing
(`NoFishing`), OFL determination (`ConstantF` at `Flimit`). The
status-determination scenario changes F partway through the projection, which
is the same per-year F primitive as Gap 1. Take the exact definitions of 6 and 7
from `spm.tpl` and §4.11.3.

**Gap 3 — the distribution over future recruitment.** SPM reports means and
percentiles over 1,000 recruitment trajectories at the MLE; Rceattle carries
one. This is a different quantity from `projection_uncertainty` (parameter
uncertainty at fixed recruitment), and the AFSC executive table uses SPM's.
SPM draws recruitment from an **inverse Gaussian** fitted to the historical
arithmetic and harmonic means, not a lognormal; Rceattle's `SIMULATE{}` blocks
draw lognormal `rec_dev`. SAFE parity needs the inverse Gaussian as an explicit
option. Mind the RNG-stream traps in `TRAPS.md` if this goes through `sim_mod()`.

## Proposed API

A forecast specification orthogonal to `build_hcr()`. `build_hcr()` keeps
saying how F responds to the stock; the new object says what the forecast
assumes. SS3's forecast file and WHAM's `proj.opts` are the models.

```r
build_forecast(
  n_years   = 14,
  catch     = c("2027" = 12408.7, "2028" = 12408.7),  # or F =, or NULL for the HCR
  avg_years = 5,            # window for selectivity, weight-at-age, M
  rec_years = 1978:2025,    # recruitment the projection draws from; also sets B100%
  rec_dist  = "invgauss",   # or "lognormal"
  n_sims    = 1000
)

fit  <- fit_mod(..., hcr = build_hcr("NPFMC"), forecast = build_forecast(...))
scen <- harvest_scenarios(fit)            # seven scenarios, spm_result-shaped
report_tables(fit, scenarios = scen)      # fills §4.11.3 and §4.2.2
```

Keep SPM reachable behind the same call, `harvest_scenarios(fit, engine =
"spm")` — spmR's own adapter pattern — so the ADMB path can be deleted the day
the native one matches. Names are provisional; `build_*()` in this repo means an
input list for `fit_mod()`, so the bridge's `build_spm()` (which compiles ADMB)
would become `compile_spm()`. Every default above is a placeholder until the
decisions below are made (CLAUDE.md rule 9).

## Sequencing

| step | work | effort |
|---|---|---|
| 1 | Move the bridge into Rceattle (`R/12-projection_spm.R`), generalized: `nsex` 1 and 2, `minage != 1`, several fisheries, species argument, SPR targets read from `data_list`, `msmMode > 0` guard | ~2 wk |
| 2 | Gap 0, `avg_years`; `/golden-check`, NEWS | ~3 d |
| 3 | Gap 1, F from catch and per-year F in TMB; `build_forecast()` | 1–2 wk |
| 4 | Gap 2, scenario runner, executive table, `report_tables(scenarios =)` | ~1 wk |
| 5 | Gap 3, recruitment draws and percentiles | 1–1.5 wk |
| 6 | Bridge becomes the acceptance harness; retire ADMB; vignette, NEWS, `_pkgdown.yml` | ~1 wk |

About 6–8 weeks to a native module a Tier 3 SAFE could run on. Step 1 is worth
doing first even though the bridge is eventually retired: it is the reference
the native module is tested against.

## Verification

- `tools/verify/verify-spm-export.R`: write inputs, read back with
  `spmR::dat2list()`, assert every field against the fit. No ADMB needed. SPM
  reads its files positionally, so a dropped row shifts every field after it
  with no error; this is the only net for that.
- Unit tests on `make_test_data()` per shape: `nsex` 1 and 2, `minage != 1`,
  two fisheries, `msmMode > 0` refused.
- Bridge test: the table above, native vs SPM on GOA arrowtooth 26.0 and 26.1.
  SSB, FOFL, FABC, OFL, ABC identical; B100%/B40%/B35% within 0.2%.
  `skip_on_cran()` and skip when no `spm` binary is present.
- Gaps 0 and 1 change the TMB model, so `/golden-check`, then `/verify` for the
  projection and simulation paths golden does not cover.

## Open decisions

1. **`avg_years` default** — terminal year (today) or a 5-year mean (SPM, SS3,
   WHAM)? Changing it moves every existing projection.
2. **Recruitment draw** — inverse Gaussian for SAFE parity, or the model's own
   lognormal `rec_dev`? Changes B100%, B40%, B35% and every percentile.
3. **Bridge first**, or straight to native, validated offline against the
   arrowtooth numbers?
4. **Multispecies** — in scope, or deferred behind an error? Under
   `msmMode > 0` the SPR block does not run, so `SPR0`, `SPRlimit` and
   `SPRtarget` are exactly zero and F40%/F35% do not exist
   (`R/0-quantity_dictionary.R:151-155`). A multispecies Tier 3 executive table
   has no defined reference points today. That is a research question about what
   a multispecies SPR proxy means, not a port.
5. **Home and dependency** — spmR is not on CRAN. The bridge needs it only for
   `system.file("admb", "spm.tpl")`, so Rceattle can take it as `Suggests` with
   `Additional_repositories`, or no dependency at all with `compile_spm(tpl =)`.

## Where to resume

Read `spm_bridge.R` and the arrowtooth driver first; they are the specification
of what works. Then settle decisions 1–3, and start at step 1.
