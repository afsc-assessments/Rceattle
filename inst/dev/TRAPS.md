# Verified traps

Each entry is something that has actually bitten, with the evidence. `CLAUDE.md` carries a
one-line version of the ones you need before touching the code; this file holds the detail and
the measured numbers so a claim can be checked rather than believed.

Verified against source, 2026-08.

## Build and test

**Dev builds are `-O2`, not pkgbuild's `-O0`.** The repo `.Rprofile` sets
`options(pkg.build_extra_flags = FALSE)`, so `load_all()` compiles the TMB model with the same
optimization as a production `R CMD INSTALL` — `fit_mod()` runs ~10x faster than an
unoptimized debug build, bit-identically (measured). A normal install was always `-O2`; only
`load_all` was `-O0`. To debug the C++ line-by-line (gdb/lldb), start R with
`RCEATTLE_DEBUG_CPP=1`. **Any absolute fit timing must state its build** — an `-O0` number
overstates real cost ~10x.

**A slow fit is the model, not a regression.** A full `BS2017SS` `estimateMode = 0,
phase = TRUE` fit is ~55 s: ~14 s for one optimization, ~3x for phasing, +~14 s for
`sdreport`. That single fit has needed ~500–700 `nlminb` iterations since at least 2023
(bisected back through `main`), so don't hunt a recent cause. For dev loops
`phase = FALSE, getsd = FALSE` gives ~14 s.

**`make_test_data()` cannot produce a positive-definite Hessian**, so it cannot exercise
anything behind `getsd = TRUE`: no `sdrep`, no confidence intervals, no
`.refit_converged()` PD check. `test-functions-retrospective.R:597` says so in passing ("the
small fixture loses every peel to it"), and `test-plot-confidence-intervals.R` reaches for
`BS2017SS` instead.

The bundled datasets that DO converge with a PD Hessian, measured 2026-08-30 with
`fit_mod(estimateMode = 1, msmMode = 0, phase = FALSE, fit_control(getsd = TRUE))`:

| dataset | nspp | nsex | sdrep | pdHess | convergence |
|---|---|---|---|---|---|
| `GOAcod` | 1 | 1 | yes | yes | OK |
| `GOAatf` | 1 | 2 | yes | yes | OK |
| `NorthernRockfish2022` | 1 | 1 | yes | yes | WARN |
| `GOAatf2023` | 1 | 2 | yes | yes | FAIL |
| `GOApollock` | 1 | 1 | no | -- | FAIL |
| `Atka2022` | 1 | 1 | -- | -- | fits from 5.25.0; `data_check` refused its IID non-parametric selectivity before |

`GOAcod` and `GOAatf` are the two clean ones, and `GOAatf` is the only bundled **two-sex**
fit that converges -- which makes it the fixture for anything sex-structured, including the
SPR-per-recruit paths.

**The parallel-test DLL failure is a toolchain failure, not a test failure.**
`DESCRIPTION` sets `Config/testthat/parallel: true`, and the workers cannot load a freshly
rebuilt DLL, so `devtools::test()` aborts before running anything with
`testthat subprocess failed to start … getDLLRegisteredRoutines.DLLInfo`. After a `.cpp`/`.hpp`
edit run `pkgload::load_all(".")` first, then test with `TESTTHAT_PARALLEL=false`.
(`test-coverage.yaml` forces serial for a related reason.)

**`load_all()` does not notice a header-only edit.** TMB's `compile()` compares `ceattle.cpp`'s
mtime against `ceattle.so`'s and nothing else, so editing only a `.hpp` — where most of the
process logic lives — recompiles nothing and every measurement afterwards is against the old
model. It is silent: `load_all()` reports success. Cost the OSA CDF work half an hour of
debugging a fix that was already correct and simply had never been built. **After a `.hpp`-only
edit, `touch src/TMB/ceattle.cpp` before `load_all()`**, and confirm the `.so` mtime moved. Same
failure mode as the `safebounds` trap below, from the same cause.

**Roxygen is pinned at 8.1.0** (`DESCRIPTION: Config/roxygen2/version`). A *different* local
roxygen2 rewrites *every* `man/*.Rd`; an older one also swaps the key back to the legacy
`RoxygenNote:`, leaving both present and contradicting each other. After `document()`, check
`git diff DESCRIPTION` **first**: if the version key moved, the `man/` churn is the version, not
your change. Either regenerate everything under one version deliberately, or `git checkout` the
unrelated churn and keep only the `.Rd` you meant to change. A matching version is necessary but
not sufficient — the build environment can rewrite `man/` and `NAMESPACE` too.

## Shared parameter blocks

**Fleets sharing a `Selectivity_index` / `Catchability_index` share ONE parameter block.**
(`Q_index` is a deprecated alias of `Catchability_index`; the schema upgrades it on read.)
`adjust_map_shared_params()` copies the donor fleet's map slice over the rest of the group, so
any per-fleet setting that differs within a group is silently taken from the donor
(`Sel_start_year`, `Bin_first_selected`, `N_sel_bins`, …). The donor is the first *estimated*
fleet — an `Off` fleet's slice is all `NA` and must never lead. Penalties and priors on a shared
block are accumulated once, on the lead fleet (`flt_sel_lead` / `flt_q_lead`); without that gate
they are counted once per sharing fleet.

**Worked example: GOA2018SS.** Fleets 1 and 7 share selectivity; fleets 9 and 10 share
selectivity *and* q.

**There are TWO lead rules and they do not agree.** Which one is right depends on whether you
are asking about the parameter block or about the penalty:

| rule | grouping key | who uses it |
|---|---|---|
| `.shared_block_lead()` (`R/0-linkage_table.R`) | `Selectivity_index` **alone** | the map — `adjust_map_shared_params()` shares on this, so it is the rule for "does a linkage or prior here free anything?" |
| `flt_sel_lead` (`R/5-rearrange_data.R`, via `.group_lead()`) | `Selectivity_index` **and** `Selectivity` | the template — `ceattle.cpp` gates the penalty block on it, so it is the rule for "does this fleet's penalty weight get read?" |

They coincide on every group that shares one form, which is every group in every bundled data
set and in all 375 consumer-repo workbooks, so the divergence hides. It opens when two live
fleets share an index with **different** forms: the template makes them two groups and charges
both, while the map still copies one block over both. `data_check()` only *warns* about that
(`R/1-data_check.R`, the `.sel_shaping_cols` check), so the configuration reaches a fit.

5.42.0's negative-`Sel_curve_pen` guard borrowed the map's rule to decide which fleets to
check, and so skipped a fleet whose weight `ceattle.cpp` does read — a negative weight there
rewarded the deviation without bound while `data_check()` reported clean. Fixed at 5.42.1 with
`.rce_sel_pen_lead()`. **If you write anything that predicts the template's behaviour, group
by index AND form; if you write anything about the map, group by index alone.**

Still open, and pre-existing: in that same mixed-form state the one shared `sel_coff` block is
penalized **twice**, once per group, because the map shares more widely than the penalty gate
groups. `CLEANUP_BACKLOG.md` carries it.

**A fixed-width parameter slot does not have a fixed meaning.** Several blocks are declared
`[…, 2]` or `[…, 3]`, and what the slot holds depends on a switch, so a static label names the
wrong quantity:

- `sel_curve_pen` is mapped out entirely by `build_map_selectivity()` and re-enabled only under
  `Selectivity` 6 (2DAR1, slots 1–2) and 7 (3DAR1, slots 1–3), where the slots are logit-scale
  AR1 correlations across **bins, years and cohorts**. The `fleet_control` reading of the same
  column — the shape/curvature penalty weights for the non-parametric and LogisticPM forms — is
  never estimated, so an element that reaches `obj$par` is always a correlation.
- `sel_inf[2]` is a descending inflection only for the double-logistic family. Under
  `DoubleNormal` (8) it is `logit(right_floor)`; under `LogisticPM` (11) it is the free age-1
  log-selectivity. `log_sel_slp` under 8 is `log(sigma_asc)` / `log(sigma_desc)`.
- `rec_pars(sp, 1)` is the SRR **alpha** on the log scale (`alpha = exp(rec_pars(sp,1) + …)`,
  `ceattle.cpp`), not steepness. Steepness is derived from alpha and `SPR0`.

`.PAR_SEL_SLOTS` (`R/0-parameter_index.R`) is where the switch-dependent ones are resolved.

**`parameter_dictionary()$dims` is prose, and it has been wrong.** A rank-only check passes a
token that names the wrong axis of the right arity. Three were wrong until 5.26.0: `log_F`
declared `[n_fsh, nyrs]` against an `[n_flt, nyrs_hind]` array, `log_pop_scalar` `[nspp, nyrs]`
against an **age** axis (which would have printed ages as calendar years), and `M1_dev_log_sd`
`[nspp, nsex, 2]` against `[nspp, nsex]`. `n_sel` and `n_fsh` both equal `nrow(fleet_control)`;
sharing happens through the map, not through a narrower array. `init_dev` is the one block whose
declared extent is deliberately the estimated portion (`nages-1`) of a wider (`nages`) array.
`test-schema-parameter-index.R` now checks rank *and* extent.

**A shortened axis also starts somewhere.** `nages-1` does not say which end it drops, and
`init_dev` drops the youngest: `ceattle.cpp` reads `init_dev(sp, age - 1)` over
`age = 1 .. nages-1`, so its ages are `minage + 1 .. minage + nages - 1` — age `minage` in the
first year is recruitment. `minage = 1` everywhere bundled hides the shift, so reading the axis
from `minage` named every `init_dev` one age too young (verified on `BS2017SS`: labelled 1–11
and 1–20, actually 2–12 and 2–21) and never named the oldest. `.PAR_AXIS_OFFSET`
(`R/0-parameter_index.R`) is the registry, and `test-schema-parameter-index.R` requires an entry
for every shortened token the dictionary uses.

## `fit$data_list` is the pre-`rearrange_data()` list

**Nothing `rearrange_data()` derives survives onto a fitted object.** Its output goes to
`data_list_reorganized`, and `mod_objects$data_list` (`R/6-fit_mod.R:1664`) is the *input* list
plus `calc_mcall_ianelli()`'s two numeric `fleet_control` columns. So `fit$data_list` has no
`flt_sel_lead`, no `flt_sel_type`, no `bin_first_selected` — those exist only on
`fit$obj$env$data`.

Reading one off `fit$data_list` gets `NULL`, and a helper that early-returns on `NULL` then
degrades **silently**. `.rce_sel_lead()` was briefly rewritten to resolve a mirrored fleet's
selectivity lead from `flt_sel_lead`; it returned the fleet itself on every real fit, so
mirrored fleets drew no confidence band beside a banded lead — the exact "some panels certain,
some not" reading the whole-figure decline exists to prevent. Measured on `Atka2022` with fleet
2 mirroring fleet 1: band at ages 1-3 went `0.0015, 0.0383, 0.4467` -> `NULL`, and a real
`GOApollock` hindcast resolved `Pollock_survey_1_shelikof_acoustic_BS` to lead `7` instead of
`1, 7`. Bundled datasets with shared selectivity groups: `GOApollock`, `GOA2018SS`,
`GOAatf2023`, `GeorgesBank3spp`.

**The unit test passed throughout**, because it was written against a hand-built list carrying
`flt_sel_lead`. A fixture that supplies the field under test cannot detect its absence — assert
on a fitted object. Either recompute from `fleet_control` (with `.group_lead()`, so the key
stays the template's) or read `fit$obj$env$data`.

## `fit$obj` is the projection's under any HCR but `NoFishing`

**A fit run at the default `estimateMode = "Estimate"` does not end on the object that
estimated the hindcast.** `build_hcr_map()` (`R/0-build_hcr.R`) replaces every entry of
`map$mapList` with `NA` and turns on `log_Ftarget` / `log_Flimit` alone, `fit_mod()` rebuilds
`obj` against that map, and for every HCR but `NoFishing` and `ConstantF` it re-optimizes — so
`mod_objects$obj` and `mod_objects$sdrep <- opt$SD` both belong to the projection. Measured on
`Atka2022` with `HCR = 5` (2026-09-02):

| HCR | `fit$obj$par` | `nrow(sdrep$cov.fixed)` |
|---|---|---|
| `NPFMC` (Tier 3) | **2** (`log_Flimit`, `log_Ftarget`) | 2 |
| `ConstantF` | **1** (`log_Flimit`) | **584** — built but never re-optimized |
| `NoFishing` | 584 | 584 |
| `NPFMC`, `projection_uncertainty = TRUE` | 586 | 586 |

`fit$.conv_hindcast$par` is 584 in every row. Two things fall out of it, and both are silent:

- **Anything indexing `obj$par` by position describes the wrong parameters.**
  `fit$identified` and `fit$.conv_hindcast` are captured *before* the projection, so labelling
  one from an index built on the other named `log_Flimit` where `rec_pars` was flagged.
  `.conv_index_ok()` / `.conv_index_for()` (`R/0-parameter_index.R`) now check the index against
  the block names of the vector being labelled, and `fit_mod()` stores the hindcast index on
  `.conv_hindcast$index` because that is the last point holding the hindcast `obj`.
- **Every delta-method standard error of a hindcast quantity is exactly 0**, since the
  projection `sdreport` holds those parameters fixed. On the same fit all 1,012
  `log_sel_at_age` errors are 0. `fit_mod()` warns, and `plot_selectivity(add_ci = TRUE)`
  declines rather than drawing a zero-width band that reads as certainty.
  `projection_uncertainty = TRUE` re-enables the hindcast parameters and is the fix that keeps
  the HCR — but **not** under `estimateMode = "Projection"`, where `build_map(debug = TRUE)`
  has already mapped the hindcast off and nothing turns it back on. `fit_mod()` warns there too.

`ConstantF` is the awkward middle: the projection `obj` is built but never re-optimized, so
`sdrep` stays the hindcast's beside a projection `obj`. That is why the checks verify rather
than switch on `estimateMode`.

## Silent-wrong-number traps

**Which fleet leads a `Selectivity_index` group is row-order dependent, so a group's penalty
weights can change meaning when rows move.** `.group_lead()` picks the group's first fleet that
is not `Off`, and `Fleet_code` must equal the row number, so inserting or reordering a fleet —
or switching the lead `Off` — promotes a different row. Only the lead's `Sel_curve_pen1/2/3` are
read (`ceattle.cpp:4051` gates on `flt_sel_lead(flt) == 1`), and a follower's are neither read
nor checked against the lead's: `.sel_shaping_cols` deliberately excludes them, and the negative-
weight refusal skips followers for the same reason. So a stale or wrong weight sitting on a
follower is inert until a reordering makes that fleet the lead, at which point it is charged
silently. Keep the whole group's penalty columns in agreement even though nothing enforces it.

**A `NonParametricPM` (9) `RandomWalk` fleet at the default `Sel_curve_pen3 = 0` has an exactly
flat direction.** Slot 3 is the only term charged on the RAW `sel_coff_dev`
(`ceattle.cpp:4233`); every other term — shape, curvature, the random walk, and the data —
reads `log_non_par_sel`, which `selectivity.hpp:449-459` mean-centres per year. A shift common
to all bins in a year is therefore invisible everywhere except slot 3, verified by objectives
bit-identical across `Time_varying_sel_sd` of 0.1, 0.2, 1 and 3 while the deviates moved. With
the schema default of 0 nothing scores that direction at all, so under `random_sel = TRUE` the
Laplace inner Hessian is singular. Set `Sel_curve_pen3` (or `Sel_devmag_sd`) on a form-9
random-walk fleet. Pre-existing; found reviewing 5.42.0.

**A single-species PFMC fit's `Flimit` depends on the `log_Ftarget` start value, a parameter
PFMC never estimates.** Measured 2026-09-14 on `make_test_data()`, `HCR = "PFMC"`,
`estimateMode = 2`: fitted `Flimit` is 0.4812 from the `build_params()` start of 0, 0.5853 from
3, and 0.7379 from -999 (the value a `ConstantF, Ftarget = 0` fit stores). The cause is in
`ceattle.cpp` section 5.12: `Ftarget(sp) = Flimit(sp) + QnormHCR(sp)` is assigned inside the
fleet/sex/age/year loop at the first projection year, so for the first fishery, first sex and
age `minage` the hindcast years of `Ftarget_at_age` still carry `exp(log_Ftarget)`. `SPRtarget`
reads a hindcast year, and the single-species reference-point penalty pulls `Flimit` toward
whatever that start implies. Every single-species PFMC fit, the hake EM included, carries the
start-0 version. Fix: assign `Ftarget` from `Flimit` before the F loop; see
`CLEANUP_BACKLOG.md`.

**A Pearson residual must divide by the effective sample size the likelihood used.** Fixed in
5.29.0; the structure that caused it is still there, so a fourth composition-like source would
repeat it.

- `Comp_weights` multiplies the multinomial log-likelihood (`ceattle.cpp:3790`), so it is an
  effective sample size — the model's own simulator says so, drawing at
  `n_nom * comp_weights(flt)` (`:3870`). Dividing by the raw `Sample_size` gives a residual
  `1/sqrt(w)` times the right one: too small on an upweighted fleet, too large on a downweighted
  one. `MultinomialAFSC` is the schema default, so this was the *common* case, and
  `reweight_comps()` tunes either side of 1.
- A Dirichlet-multinomial is additionally overdispersed by `(n + conc)/(1 + conc)`, so those
  residuals were inflated into apparent misfit.
- **Three sources, three different alpha constructions.** comp (`ceattle.cpp:3776`) and CAAL
  (`:3931`) build the concentration on the offset-inflated total, so `sum(alpha) = N·S²·theta`
  with `S = 1 + comp_offset·nbins`; diet (`:5009`) renormalizes first, so it is the clean
  `N_s·theta`. Sharing one formula across all three is wrong for two of them — measured up to
  1% in the sd on the smallest bins from the offset alone, at the default `comp_offset = 1e-5`.

`.rce_comp_pearson()` and `.rce_diet_pearson()` are the only places this is resolved.
`test-likelihood-pearson-effective-n.R` checks it by simulation — data drawn under the family
the likelihood assumes must return residuals with `sd = 1` — because a test that restates the
formula passes on a wrong variance.

**The OSA tail statistic and its null must be the same estimator.** Also 5.29.0. The observed
value was `quantile()`'s type-7 interpolation; the null was simulated from type-7 quantiles, so
the pair was consistent but seed-dependent. Replacing only the null with the closed-form order
statistic looks like a strict improvement and is not: measured coverage of the nominal 95%
interval falls to **0.855 at n ≈ 50**, and oscillates non-monotonically in n because
`round(q(n+1))` steps while type-7's `1 + (n-1)q` moves continuously. Both sides are now
`sort(resid)[r]` against `Beta(r, n-r+1)`; measured coverage 0.949–0.951 for n = 5..800.
`r` must be clamped to `[1, n]`: unclamped, `qbeta(p, n+1, 0)` returns 1 and `qnorm(1)` is `Inf`,
so `upper_ok` is TRUE for every series with **n ≤ 19** — a silent pass, not a visible gap.

**A reference point CEATTLE never estimated is a NUMBER, not a gap.** Three cases, all
verified against the 2026 GOA three-species assessment (2026-08-29):

- `Ftarget` / `Flimit` are estimated only under an HCR that defines them (`build_hcr_map()`
  turns them on for CMSY / ConstantF / ConstantFSSB / ConstantFSPR / NPFMC / SESSF / PFMC only),
  and are switched off per species when `Proj_F_proportion` sums to 0 or `estDynamics > 0`.
  Unestimated they sit at `exp(0) = 1`, so that assessment reports **Ftarget = Flimit = 1.0/yr**
  for all three species — an enormous F that reads as an estimate. Gate on
  `build_hcr_map(data_list, fit$map)`. **Do NOT read `fit$map$mapList$log_Ftarget` instead**:
  `R/3-build_map.R:1536` sets both to `NA` in the hindcast map whatever the HCR, so that test says
  "never estimated" for every model.
- Under `msmMode > 0` the template overwrites `SB0` / `B0` with the `MSSB0` / `MSB0` inputs
  (`ceattle.cpp:2189`), which stand at `.RCE_MSSB0_PLACEHOLDER = 999` mt until `fit_mod()`
  derives them from a no-fishing projection. Measured on that fit: `SB0 = 999` for all three
  species against true terminal SSB of 1,136,070 / 564,108 / 96,504 mt, so a
  `B_target = Ptarget * SB0` proxy comes out at **399.6 mt**. `data_list$MSSB0_derived` is the
  per-species flag — carried precisely because recognising the placeholder by comparing a double
  to 999 would also null a genuine 999 mt.
- The per-recruit quantities (`SPR0`, `SPRlimit`, `SPRtarget`, `SPRFinit`, `NbyageSPR`) are
  inside `if(msmMode == 0)` and are **exactly zero** on a multispecies fit.

`report_tables()` returns all three as `NA` with the reason in `basis`.

**The depletions do NOT simply divide by `SB0`, so do not blank them alongside it.** Section
12.1 divides by `DynamicSB0` under a dynamic rule, and by biomass in the LAST projection year —
the equilibrated unfished reference — when `HCR == 0 & msmMode > 0`. Only the remaining case
reads the `SB0` array. Blanking depletion whenever `SB0` is a placeholder destroys a valid
series on exactly the multispecies fits the package exists for; measured range on that fit is
0.169–1.724, which is visibly not `ssb / 999`.

**`sex_ratio` enters SPR differently for a one-sex and a two-sex species.** SPR is spawning
output per TOTAL recruit, so the female fraction belongs in it for every species -- once. A
one-sex species' schedule is sex-combined, so the age-varying ratio applies at every age (that
is what `mature_females` folds in, section 5.4). A two-sex species' sex-0 schedule is already
female, so only the recruitment split `sex_ratio(sp, 0)` applies, the same split 6.6 uses for
`NByage0` / `NByageF`.

Until 5.24.1 section 6.2 applied the age-varying ratio to both. On a two-sex species that
re-applied a split already in the schedule, and returned NaN wherever the table stopped short of
the species' own `nages` -- measured: GOA arrowtooth (`nsex = 2`, `nages = 21`, table filled to
10) returned NA for `SPR0`, `SPRtarget`, `SPRlimit` and `SPRFinit`, so no Tier 3 proxy could be
formed. Fixed value 1.2300187 kg/recruit.

**Both wrong answers are easy to reach.** Applying the age-varying ratio double-counts;
dropping it entirely reports per FEMALE recruit, which is not comparable with the one-sex
species in the same model. `test-dynamics-brps.R:229` and `test-dynamics-spr.R:168` pin the
per-total-recruit convention on fixtures whose sex_ratio is 0.5 at every age, so the two
branches must agree there -- they are what catches either mistake.

After the fix a two-sex species reads `sex_ratio` only at age 1, so `data_check()` requires the
full schedule for a one-sex species only. **SPR reaches the objective only under
`HCR = "ConstantF"` or an SPR-based rule (`HCR > 3`)**, so most fits cannot move at all.

**A grep for `REPORT(` over `ceattle.cpp` over-counts.** `mature_females`, `sex_ratio_hat`,
`N_at_age_dB0` / `N_at_age_dBF` and the whole Kinzey functional-response block sit behind
comments. A fit reports **99** quantities, identical single- and multi-species; enumerate from
`names(fit$quantities)`, not from the source.

**A missing OSA CDF term returns 0, not an error.** `oneStepPredict(method = "cdf")` reads
`Fx = 1 / (1 + exp(nlcdf.lower - nlcdf.upper))`. A likelihood that supplies no
`keep.cdf_lower` / `keep.cdf_upper` term makes those two evaluations identical, so `Fx = 0.5`
and every residual on that data source is exactly 0 — a clean-looking Q-Q plot from a
likelihood that was never asked anything. The Dirichlet-multinomial composition families are in
exactly that position (no elementary beta-binomial CDF), so `osa_residuals()` routes those
fleets to a Gaussian method and says so; `test-likelihood-osa-cdf.R` asserts the routing and
that the residuals are not degenerate. **Any new observation likelihood owes a CDF term the way
it owes a `SIMULATE` draw.**

**Squeezing an OSA CDF at exactly one machine epsilon lands on a rounding tie.** `Fx` above is
recovered in double precision, and at a tail of `DBL_EPSILON/2 = 1.11e-16` the denominator is
`1 + 1.11e-16` — the round-half-to-even tie above 1 — so whether it rounds to 1 depends on the
last bits of an objective summed over thousands of observations. When it does, `Fx = 1` and the
residual is `+Inf`. Measured on BS2017SS: 6 of 4538 composition residuals infinite at one
epsilon, 0 at four (`osa_squeeze_cdf()` in `comp_osa.hpp`). The reference implementation TMB
ships uses `squeeze()`, i.e. one epsilon, and has the same latent failure. The cost of four is
where `|residual|` saturates: 8.04 rather than 8.13.

**`Index_distribution` has a second hand-synced registry.** A family added to
`index_distribution_map` (`R/0-switches.R`) must also be classified in
`.index_rows_natural_scale()` in the same file, which is what `residuals(type = "pearson")`,
`plot_index()`'s observation interval and `plot_indexresidual()` read to pick a scale. A
natural-scale family that misses it **does not error** — it silently gets the log-scale formula,
whose `sigma^2/2` is then a number the size of the index squared.
`test-index-natural-scale-downstream.R` enumerates the map and is the net.

**`nages` is a COUNT of age bins, not the oldest age.** Ages run
`minage .. minage + nages - 1`, and the array's 3rd dimension is a 1-based bin index, so age `a`
is at index `a - minage + 1`. `minage = 1` makes the two coincide — which is every bundled
dataset and all three live assessments — so a bin/age confusion passes every test and every real
model, and only shows up on the next `minage != 1` species. Write
`seq_len(nages[sp]) - 1 + minage[sp]` and mean it.

**`jnll_comp` columns mean different things on different rows.** The matrix is dimensioned
`max(n_flt, nspp)`, and the column index is whatever loop variable wrote the cell: rows 1-8
(index, catch, composition, CAAL, non-parametric selectivity, selectivity deviates, catchability
prior, catchability deviates) are written inside `for(flt = 0; flt < n_flt; flt++)`, rows 9-20 by
`sp`/`rsp`/`slot_col`, and "Linkage random effects" always into column 1. So `jnll_comp[3, 2]` is
fleet 2's composition likelihood while `jnll_comp[11, 2]` is species 2's recruitment deviates,
and `rowSums()` pools across two different axes. `.JNLL_ROW_AXIS` (`R/9-profile.R`) is the
registry; `test-schema-jnll-rows.R` parses every `jnll_comp(JNLL_*, col)` write in the template
and asserts each row's declared axis matches the column it is actually indexed by. Verified
2026-08-26: 124 writes, all 21 rows covered; 133 writes at 5.41.0, same 21 rows.

**`unweighted_jnll_comp` is populated for 5 of its 21 rows.** It exists so Francis and
McAllister-Ianelli can read a composition likelihood without its `Comp_weights` multiplier, so
only the rows carrying such a multiplier are written: composition, CAAL, stomach content, and the
two linkage rows. Index, catch, selectivity, catchability and every penalty are **structurally
zero** there, not small. Anything that treats the two matrices as the same quantity at different
weights — a reweighting diagnostic, a component profile, a "which term dominates" table — reports
those components as absent rather than unweighted. `profile_components(weighted = FALSE)` says so
in its `@param`.

**`retrospective(getsd = TRUE)` can return fewer peels than `getsd = FALSE`, so Mohn's rho
differs.** `.refit_converged()` drops a peel on a non-positive-definite Hessian, and that check
can only run when an `sdreport` was asked for (`R/0-convergence.R:231`, deliberate). A peel that
converges on gradient but has a singular Hessian is therefore kept without standard errors and
dropped with them. On `make_test_data()` every peel is lost this way. Rho is computed over the
peels that survive, so **two runs of the same model can report different rho** depending only on
`getsd`. `test-functions-retrospective.R` compares the two runs peel by peel over the shared
names for exactly this reason.

**~~`fit_mod(d, config = cfg)` replaces `d$model_config`, linkages and all.~~ Fixed in
5.36.0** (`65259384`): `R/6-fit_mod.R:333-363` overlays only the fields the config's `"set"`
attribute names, keeps the data's for the rest, and warns where the two differ. Kept because the
measurement is what the fix was written against, and because a config saved before that release
still carries no `"set"` attribute -- those set their non-default fields, which is the old
behaviour for any field the author did change. Measured on 5.29.0, `GOAatf` with a survey
`slp_asc ~ rw(1 | Year)` linkage:

| call | random effects |
|---|---|
| `fit_mod(d)` | 57 (`beta_linkage_re`) |
| `fit_mod(d, config = run_config(d, random_sel = FALSE))` | 57 |
| `fit_mod(d, config = run_config(model_config(), random_sel = FALSE))` | **0** |

Before 5.36.0 there was no warning: the fit ran and returned a model without the time-varying
survey selectivity. Building the config from the data (`run_config(d, ...)`) or from a fit is
still the habit worth keeping, and it is the only thing that works for a config written before
the attribute existed.

`random_sel` does not reach the linkage REs either way — it gates only the `fleet_control`
`Time_varying_sel` deviations (`R/6-fit_mod.R:852`); linkage REs are integrated whenever
present, unless the spec sets `integrate = FALSE` (`:923`). `random_q` is the same: it gates
only `index_q_dev` (`:849`). The config description that misled an external user into reading
`random_sel` as a switch over linkage REs was corrected in 5.36.0 and now says
"Integrate the Time_varying_sel deviations and estimate their sd (linkages integrate either
way)" (`R/0-save_config.R:309`), so the trap is the behaviour, not the wording.

**`fit_mod(initMode =)` never reads the data object's own `initMode`.** The argument defaults to
`"NonEquilibrium"` (`R/6-fit_mod.R:185`) and is written to the data list unconditionally
(`:424`), so a value set on the data is overwritten before anything reads it. `BS2017MS$initMode`
is 1 (`Equilibrium`), and the golden `ms` fit -- which passes no `initMode` -- initializes as
`NonEquilibrium`. A
`model_config` slot on the data is honoured (`:385`, since 5.36.0); the bare field is not.

**A `data_list` element with no `write_data()`/`read_data()` support round-trips to nothing.**
The feature is then silently lossy through the standard xlsx format. This is how `index_cov` was
lost.

**A `Comp_weights` of 1 under a Dirichlet-multinomial is a starting weight of e.** That
likelihood reads the column as a log (`DM_pars_comp = exp(comp_weights)`); a multinomial reads it
as a natural-scale multiplier. `write_template()` seeds 1. `fit_mod()` reports this since 5.11.0.

**Golden reference: `BS2017SS = 10241.0304272585`,** measured on `dev` 2026-08-21 with the
`/golden-check` recipe (every fit at `newtonsteps = 3`). Other values have been quoted in this
repo's tooling and were wrong. The constants are platform-independent (verified across Windows
and macOS/arm64 to 16 digits), but on a new machine capture the baseline **before** the change
you want to test — a failure you cannot attribute is worse than no check.

## The refit path, and the eight callers

**Eight entry points re-invoke `fit_mod()` through `.refit_like()`** (`R/6-refit_like.R`), which
rebuilds the HCR / SR / M1 / growth specs from a source `data_list` and exposes each per-caller
divergence as a named override:

`retrospective()` · `jitter()` · `self_test()` · `profile()` · `run_mse()` · `remove_F()` ·
`sample_rec()` · `reweight_comps()`

`tools/verify/verify-refit-like.R` covers **six**. `sample_rec(update_model = TRUE)` and
`reweight_comps()` are **not in it** and must be checked by hand.

**`run_mse()` pins the OM's stock-recruit and suitability reference windows to the pristine
`om$`**, not the advancing `om_use$`, so the hindcast does not drift through the projection —
essential for multispecies, whose predation suitability must stay fixed. In `.refit_like()` these
are the `srr_mse_switchyr` / `srr_hat_styr` / `srr_hat_endyr` / `suit_styr` / `suit_endyr`
overrides. The EM instead advances `srr_mse_switchyr` to its current assessment `endyr` each
iteration. `tools/verify/verify-mse-hindcast-invariant.R` checks the invariant.

## The SIMULATE contract

Every observation and process error is drawn in a `SIMULATE{}` block **beside the density that
scores it** (`ceattle.cpp` sections 5.12b and 5.13, plus one per likelihood slot; the
multivariate draws live in `comp_sim.hpp`). `sim_mod()` implements no observation model in R — it
calls `obj$simulate()` once and writes the result back. **A new likelihood family therefore owes
a draw.** Three rules, with their reasons:

1. **Draw what the density assumes** — bias-correction convention and scale included. A draw that
   does not match its own density makes a self-test look biased when the model is fine.
2. **REPORT under a `*_sim` name.** TMB never clears the report environment, so a draw REPORTed
   under the observed object's own name reads back as the data. A process draw also reports a
   `*_drawn_sim` mask, since the deviation arrays are REPORTed whole.
3. **Don't draw what the model does not define.** Two densities on one latent — the AMAK/Ianelli
   stock-recruit penalty — have no distribution to draw from. Leave it, and warn.

## MSE draws and the assessment schedule

`sim_mod()` draws **once per observation row**, so anything that changes how many rows the
operating model carries changes how far the random stream advances — and therefore every draw
after it. That makes the random stream sensitive to things that look like pure bookkeeping.

**The trap it already caused.** 5.14.0 refit the operating model only as far as `assess_yrs[k+1]`,
the *next* assessment year, because that is all the loop reads (`max_catch_hat`, for the
exploitable-biomass cap). `clean_data()` filters every data frame to `styr:projyr`, so the
shortened model carried fewer projection rows and consumed fewer draws — a count set by **when the
next assessment fell**. Two runs differing only in their assessment schedule then drew different
observation error from the *first* assessment onward, before diverging in any other way.

Measured on BS2017SS (endyr 2017, projyr 2020, `nsim = 1`, `seed = 666`), dropping the 2019
assessment moved 2019 catch by **2.1%** — a year whose advice is identical by construction, set by
the same 2018 assessment under both schedules. Per-year seeding took it to 0.53%; refitting over
the whole projection took it to 0.003%, which is optimizer noise rather than the draws. Both
landed in 5.18.0, and the 9% the shortening was worth on a Bering Sea multispecies MSE went with
it — see `inst/dev/TODO-mse-horizon.md` for how to get it back.

**What is still open.** Seeding is keyed to the **assessment** year, but one `sim_mod()` call
draws every year in the interval. A year sitting inside a longer interval is drawn under that
interval's seed, while the same year in a schedule that assessed it directly gets its own. So two
schedules realize different observation error in the years *between* assessments even where the
stock is in the same state, and that enters the next assessment as noise rather than as the
schedule. Closing it needs the draw keyed to the observation's year, not the assessment's.

**The rule.** Before changing anything that alters the operating model's row count, its horizon,
or the order of draws, ask what it does to a comparison of two schedules — not just to a single
run's reproducibility. `tools/verify/verify-mse-schedule.R`'s `gap_crn_clean` check is the gate:
it asserts a year both schedules take from the same assessment comes out the same. Note it
measures **catch**, so it does not see the between-assessment observation gap above.

## The `nages` age-vs-index class, and what closes it

`plot_ration(minage=)`, `plot_m_at_age(age=)` and `plot_m2_at_age_prop(age=)` once indexed the
age array directly while labelling the axis with an age. They now resolve the argument per
species through `.rce_age_index()` / `.rce_age_plus_index()`.

The class is closed by fixtures, not by those three fixes. `test-plot-predation-args.R` shifts a
fitted model to `minage = 2` and requires `age = 2` to draw what `age = 1` drew before;
`test-plot-smoke.R` runs the rest of the exported plotters against a `minage = 3` model. Without
a `minage != 1` fixture the whole suite passes on a bin index read as an age, because every
bundled dataset and all three live assessments have `minage = 1`.

**The same class in prose: `Bin_first_selected` is a 1-based BIN ORDINAL, not an age.**
`rearrange_data()` converts it with a bare `- 1`, while `Sel_norm_bin` beside it gets
`sel_bin_offset = minage[species]` — the two columns use opposite conventions, and the alias
`Age_first_selected` argues for the wrong one. On a `minage = 3` stock `Bin_first_selected = 3`
is age 5. `minage = 1` hides it everywhere it is currently written, so the column description in
`R/0-column_schema.R` now says which convention it is, and the selectivity documentation says
"first selected bin" rather than naming an age.

## The guards are not themselves guarded

Three of this repo's most valuable checks have been guarded off in a way nothing verified. All
the same shape: a speed optimization, or a green-looking job, that silently removes coverage.

**The bounds check measured nothing from before 5.20.0 until 2026-08-27, and reported success.**
`deep-checks`' `safebounds` job is `continue-on-error: true`, so the run's conclusion was
`success` while the job failed. It failed at the DLL, not on a violation: `NAMESPACE` loads
`Rceattle` (the package SHLIB) as well as `ceattle` (the TMB model), `verify-safebounds.R` built
only the second and then called `load_all(compile = FALSE)`, and on a fresh checkout
`src/Rceattle.so` had never existed. All five cases errored with "not found in the list of loaded
DLLs" and were printed as `BOUNDS VIOLATIONS or errors in 5 case(s)` -- reading as five
violations when the truth was zero measurements. **It passed locally only because a developer's
tree already had a `src/Rceattle.so` lying around**, which is why "five configurations clean on
macOS" in the 5.20.0 notes was not the clearance it looked like. Fixed by building through make
(`tmblib` is a `.PHONY` prerequisite of `$(SHLIB)`, so `compile.R` re-runs with the flag still
set), asserting both DLLs loaded before any case runs, and reporting "could not run" separately
from "violation". The restore path had the same hole and left the tree unloadable.

**And the step named for the crash was not bounds-checked either.** `deep-checks` ran
`test-selectivity-catchability.R` -- the file the Windows worker died in -- as a *second* step,
after `verify-safebounds.R` had restored the unchecked build on exit. `TMB::compile()` is
incremental, so that step's `load_all(compile = TRUE)` no-op'd: measured 2026-08-27,
`src/ceattle.so`'s mtime did not move, so it ran against the ordinary model. **A build flag is
not in force just because the job env sets it** -- check that the artifact was actually rebuilt.
It is a case inside the harness now, sharing the one checked build, which also made it cheaper
than the second compile a separate step needs (7.8s against ~90s).

**The golden regression ran in NO automated job.** `test-golden-regression.R` carries both
`skip_on_cran()` and `skip_on_covr()`. `R-CMD-check` sets `NOT_CRAN=false` deliberately (to keep
the multi-OS matrix fast), so the first fires; `test-coverage` runs under covr, so the second
does. The committed backstop for "a numeric change needs golden equivalence" therefore only ever
ran when someone typed `/golden-check` locally. Its measured cost of `0` in
`tools/ci/coverage-costs.tsv` is the tell — four phased fits do not cost nothing. The
`deep-checks` workflow now runs it nightly, and **asserts the run produced assertions**, because
a skipped file otherwise reports as a pass.

**And golden was only the worst case of it: no `skip_on_cran()` test could fail any job.** That
is 173 of the 240 test files — every test that runs a real `fit_mod()` optimization. `R-CMD-check`
sets `NOT_CRAN=false`, so it skips all of them. `test-coverage` does run them, but
`tools/ci/coverage-shard.R` passes `stop_on_failure = FALSE` (deliberately: a failure must not
throw away the shard's coverage before it is written) and no later step reads the result — the
"Upload test results" step is `if: failure()`, which never fires because nothing failed. So a
broken regression test printed into a shard log and the run went green. Found while reviewing
5.25.1, whose own PR body claimed a new test "fails if anyone restores the old comparison"; it
would not have. `deep-checks` now runs the whole suite nightly with `NOT_CRAN=true`, asserting
the mode **before** the fits rather than after, since `NOT_CRAN` going false presents as a fast
green run of nothing but the unguarded unit tests. **Writing a regression test is not the same as
having one run** — check which job executes it, and whether that job can go red.

**`NOT_CRAN=false` was set through `$GITHUB_ENV` and did not hold reliably.** On `7e16fa73` the
skip fired (testthat recorded `test-selectivity-catchability.R:3:1`); on `7611bd1c`, differing
only in a Markdown file, the same file executed and the worker died with exit code
`-1073741819` = `0xC0000005`, a Windows **access violation**. Set it as step-level `env:` on the
`check-r-package` step instead — step env beats both the job env and anything a prior step
exported, and it is printed in the step's own env block, so the log says which mode ran. 1 of 7
recent Windows runs slipped; when it does not slip you learn nothing.

**What the crash DUMPS is not a finding either — fixed at the source in 5.43.0.** When the
worker dies, `R CMD check` prints whatever that test file had written to stdout, under a
`── Test failures ──` heading. `test-convergence.R` used `expect_invisible(print(cv))`, which
unlike `expect_output()` does not sink output, so the file wrote three lines of a **deliberately
non-converged synthetic fixture** on every run:

```
<Rceattle convergence>  status: FAIL
  [FAIL] max_gradient     Maximum absolute marginal gradient = 4e+12 on sel_inf (...)
  [FAIL] pdHess           Hessian is not positive definite; standard errors are unavailable.
```

Those lines were read as a real convergence regression three separate times, most recently on
#161 -- and 5.43.0 naming the quantity (`sel_inf (selectivity inflection)`) made them read
*more* like a genuine finding, not less. The call is captured now and the file prints nothing,
so a dead-worker dump is quieter. **Diagnose from the traceback**
(`parallel_event_loop_chunky` -> `handle_error` -> `cli_abort` means a worker died), never from
the text it dumped.

**The file testthat names for this crash carries NO information — measured 2026-09-21.** Tests
run in parallel (`Config/testthat/parallel: true`, two workers on the runner), so when a worker
dies testthat reports whichever file that worker was holding. Three occurrences have named three
unrelated files: `test-selectivity-catchability.R`, `test-switch-string-aliases.R` (on `main`)
and `test-data-input-validation.R`. **The last of those executes nothing in CI** -- it opens with
`skip_on_cran()` and an unconditional `skip()`, and CI runs `NOT_CRAN=false` -- so a file that
ran zero lines was blamed for the fault. Do not investigate the named file, and do not read
`verify-safebounds.R`'s "the CI crash config" case as targeting anything established; that
config was chosen from one such attribution.

**It is not specific to any release line.** `main` at 5.33.0, released and unchanged, crashed
with the same exit code on 2026-09-21. Rate over the 30 most recent `R-CMD-check` runs: 2
failures, both that day, while Windows also PASSED on that day on another branch -- so
intermittent, not an image change. Everything from 2026-09-16 to 2026-09-20 passed.

**An access violation is memory corruption, not a bad optimum.** A fit in a different basin
gives a huge gradient and a failed `sdreport` — not a fault. The model is built
`safebounds = FALSE`, so an out-of-range access is a silent write into adjacent memory; that is
precisely what the 5.15.0 `age_hat` overflow was. Build with `RCEATTLE_SAFEBOUNDS=true` to turn
that into an R error naming the object, and run `tools/verify/verify-safebounds.R`. That harness
carries a **positive control** — it deletes the artifacts, rebuilds, and asserts
`-DTMB_SAFEBOUNDS` is in the compile line — because `pkgload` only recompiles when sources
change, so setting the flag alone guarantees nothing and a clean result against a stale
unchecked `.so` is meaningless. As of 5.16.0 it reports no violation over five configurations
(ragged comps both directions, joint sex, predation arrays, the CI crash config, `sim_mod()`
draws) on macOS. **That is not a clearance** — the fault is intermittent and was seen on Windows.

**Editing only a `.hpp` does not rebuild the model.** `src/TMB/compile.R` calls `TMB::compile()`,
which builds through `R CMD SHLIB`, and that compares `ceattle.cpp` against `ceattle.o` alone —
the headers it `#include`s are not in the dependency check. So a change confined to
`selectivity.hpp`, `bioenergetics.hpp`, `predation.hpp` or any other header is *not* picked up:
`pkgload::load_all(".")` prints nothing unusual and loads the previous `.so`, so the run that
follows tests the code you did not change. It fails the way a real defect does — a header guard
added to `bioenergetics.hpp` in 5.25.0 appeared not to work, and the segfault it was meant to
close was still there, because the binary predated it. Touch `src/TMB/ceattle.cpp`, or delete
`src/TMB/ceattle.o` and `src/ceattle.so`, before trusting a header-only result.

**A workflow that lives only on a feature branch does not run.** `deep-checks.yaml` triggers on
`schedule` and `workflow_dispatch`, and GitHub reads both from the **default branch** only:
`gh run list --workflow=deep-checks.yaml` answers `not found on the default branch`, the cron
never fires, and a manual dispatch has nothing to dispatch. It was added on `dev` as the guard
for the Windows access violation, so between adding it and merging to `main` the guard measured
nothing at all. Do not read a quiet nightly as a clean nightly until the workflow is on `main`;
check `gh run list --workflow=<file>` says something other than "not found".

Two further sharp edges found while building that harness, both worth knowing:
`src/TMB/compile.R` resolves `ceattle.cpp` against the **working directory**, so running it from
the repo root skips the compile, prints its messages, and exits 0; and `on.exit()` registered at
a script's top level has no function frame and never fires, which is how the harness first left
a slow bounds-checked `.so` in the tree for every later fit to use.

## The conditioning check reads correlation, not covariance

**Severity is read on the correlation of the estimates, and the thresholds were not
recalibrated.** The covariance condition number is not scale-invariant — `log_F` sits near -2,
`sel_inf` near 10 — so it mixes confounding with the units a parameter is held in. Measured
(2026-09-02):

| model | covariance | correlation | se ratio |
|---|---|---|---|
| BS2017SS (golden) | 4.6e4 | 2.0e4 | 140 |
| GOA pollock NonParametric | 1.9e5 | 1.6e4 | 126 |
| GOA pollock 2DAR1 | 2.5e9 | 1.9e8 | 13634 |

Standardising lowers every model by 0.4–1.1 orders, so at unchanged `WARN` 1e6 / `FAIL` 1e10 the
check is **less** likely to fire than before 5.26.0. Three models is not a calibration; if a
model that should have been flagged comes back `OK`, the thresholds are the first thing to
question. `data$covariance_condition_number` keeps the old number, and `pdHess` /
`sdreport_failed` are what catch an inversion that actually fails.

`fit$convergence$checks$hessian_conditioning$data$condition_number` therefore means something
different from 5.25.1 and earlier. `../Rceattle-models/Model comparison/run_rceattle.R` reads it
(`get_condition_number()`), and saved `condition_number_<scenario>.RDS` files are not comparable
across the boundary.

## Coverage gaps in the golden check

`/golden-check` fits four models and diffs the result. It does not reach:

- **The refit paths.** Eight entry points re-invoke `fit_mod()` through `.refit_like()`. Use
  `tools/verify/verify-refit-like.R`, which covers six of them —
  `sample_rec(update_model = TRUE)` and `reweight_comps()` are **not** in it.
- **The simulation draws.** `tools/verify/verify-sim-*.R`.
- **Any figure.** A plotter can change every number it draws and stay green. The net is
  `tests/testthat/test-plot-*.R` plus a before/after `ggplot_build()` diff at the merge base.
- **Linkage rows.** The four golden models carry none. The constructive linkage tests
  (`test-dynamics-recruitment-linkage.R`, `test-linkage-random-effects.R`) are the net.
- **Selectivity normalization.** No bundled dataset enters the normalization block in
  `selectivity.hpp` — every one has `Sel_norm_bin` all-`NA` except one-sex `Atka2022` — so the
  four reference fits pass trivially through any change to it. `GOA2018SS` *is* reachable
  (**species 2**, arrowtooth, is the two-sex one -- `nsex = c(1, 2, 1)`; 9 is the `Fleet_code` of
  `ATF_bottom_trawl`, not a species index) once you set a named `Sel_norm_bin`; that is exactly how
  `681f3ba0` shipped with `test-selectivity-normalization.R` failing.
- **Growth.** `growth_fun` is unset on BS2017SS and GOA2018SS, so `growth.hpp`'s SDAA path —
  `sd_plus_group` included — is never exercised. Covered constructively by the
  `test-growth-*.R` files.

**`goa_ms` (fixed-M GOA multispecies) sits on a flat likelihood ridge:** the same objective at
different `par`/`ssb` across *different* code, though deterministic on same-code re-runs. Judge
it on `obj`/`jnll`, not `par`/`ssb`.

**`goa_ss` has a second local minimum of the negative log-likelihood 52.9 units higher, and a
one-ULP change in one gradient element is enough to send `nlminb` there.** Recorded on the
5.34.0 branch (PR #144, 2026-09-14): adding code the objective never evaluates left every
objective bit-identical and bounds-checked builds clean, but changed one `log_F` gradient
element by 3e-16 (presumably summation order), and from there `nlminb` reached the higher
minimum with `newtonsteps = 3` in place. HEAD reproduces the reference (12867.9902664788). The
52.9 that `golden-check.md` attributes to tolerance-stopping (commit `1a172677`) is the same
gap; polishing did not remove it. A `goa_ss` delta of 52.9 with the other three models
bit-identical is this, not a numeric regression. Diagnose it from the gradient at the reference
`par`, not from the objective. A robustness fix (a warm start from the reference `par`, or a
second start keeping the lower minimum) is open; see `CLEANUP_BACKLOG.md`.

## Prior centring shares `bias_adjust_proc` with the recruitment deviations

From 5.33.0 every lognormal prior and the Ianelli penalty are centred at `-sd^2/2` when
`bias_adjust_proc = TRUE`, so a prior value is a mean; with `FALSE` it is a median, but the
recruitment deviations lose their centring too. To reproduce a pre-5.33 fit's priors with the
flag on, shift the inputs instead:

| Prior | Pre-5.33 input `v` becomes |
|---|---|
| `prior_lognormal(p1, s)` | `prior_lognormal(p1 + s^2/2, s)` |
| q prior | `Catchability_init = v * exp(s^2/2)` (also moves the starting value) |
| Ricker alpha prior | `srr_prior = v * exp(s^2/2)` |
| M prior (old centre `log M + s^2/2`) | `M_prior = v * exp(s^2)` |

The one case no input reproduces is the Ianelli penalty with centred deviations: the penalty has
no input of its own and shares `R_sd` with them. The GOA northern rockfish bridges are not this
case: they fit no stock-recruit curve, and urm's median-centred q and M priors return with
`Catchability_init = 1.107` (1.0 * exp(0.45^2/2)) and `M_prior = 0.06008` (0.06 * exp(0.05^2/2)). A `bias_adjust_proc` between 0 and 1 gives priors (DM weights included) a centre that is
neither mean nor median.

**In the penalty years the Ianelli form penalizes recruitment twice, and from 5.33.0
`bias_adjust_proc` centres both.** `rec_dev` is penalized around `R0` (`JNLL_REC_DEV`) and log R
around `R_hat` (`JNLL_SRR_PENALTY`). With flag value b, their product in log R is
N((log R0 + log R_hat)/2 - b*sigma_R^2/2, sigma_R/sqrt(2)). Before 5.33.0 the penalty lacked its
`-b*sigma_R^2/2`, so the recruitment the pair favours now falls by `exp(-b*sigma_R^2/4)`. At b = 1
that is 8.6% at sigma_R = 0.6 and 22% at 1.0, for both the prior mean (from `sqrt(R0 * R_hat)`)
and the median a penalized fit shrinks towards (from `sqrt(R0 * R_hat) * exp(-sigma_R^2/4)`).

Only centring exactly one term gives a mean of `sqrt(R0 * R_hat)`; with the flag off the mean is
`exp(+sigma_R^2/4)` high. The curve term alone still treats `R_hat` as the mean of R, which the
reference points and a `proj_mean_rec = FALSE` projection read. The data dampen the shift on
hindcast R. Accepted for 5.33.0 (2026-09-12); removing it would move the hake baselines again.
`test-likelihood-prior-bias-adjust.R` pins the formula, not this property.

## Porting an ADMB model

Absorbed from `ADMB-CONVERSION.md`, which this file replaces; nothing referenced it.

**An ADMB `dev_vector`'s sum-to-zero constraint cannot be replicated in TMB, and a straight
port is unidentified.** The deviations and the parameter they deviate from trade off freely and
the optimizer wanders the ridge. Fix: map the first element to `NA` in `R/3-build_map.R`, which
pins the level. **Exception:** if the vector already carries its own penalty -- a normal
penalty or a random-effects density -- that pins it, and turning off the first element as well
over-constrains it.

**AMAK conventions are reproduced deliberately, not corrected.** `src/TMB/selectivity.hpp:609`
evaluates the logistic at mid-age (`age_vector(j) = j + 0.5`, so `bin + 1.5`) for `LogisticPM`
(11); the standard `Logistic` (1) uses `bin + 1`. `NonParametricPM` (9) and `LogisticPM` exist
to match ADMB AMAK's "pm" parameterizations, penalties included, and
`NonParametricIntegrable` (13) came from the same lineage in 5.40.0.

The literature citations through `src/TMB/` are the specification for those blocks, not
historical notes; `CLAUDE.md` says so under Comments.

## CLAUDE.md's traps in full

`CLAUDE.md` carries each of these as one line. This is the fuller text it held until 2026-09-14,
kept verbatim (only a stale line number replaced by the call it pointed at). Some repeat a
section above.

- **`Index_distribution` has a second hand-synced registry** — a family added to
  `index_distribution_map` must also be classified in `.index_rows_natural_scale()`
  (`R/0-switches.R`), or it silently gets the log-scale residual formula.
- **`jnll_comp` columns count fleets on rows 1–8 and species on rows 9–20**, so `rowSums()`
  pools across two different axes. Row 21 (linkage random effects) is model-wide.
  `.JNLL_ROW_AXIS` (`R/9-profile.R`) is the registry.
- **A reference point CEATTLE never estimated is a number, not a gap** — `Ftarget`/`Flimit` sit
  at `exp(0) = 1` unless the HCR estimates them (gate on `build_hcr_map()`, never on `fit$map`),
  `SB0` under `msmMode > 0` is the 999 mt `MSSB0` placeholder until `MSSB0_derived` is TRUE, and
  the per-recruit quantities are zero outside `msmMode = 0`.
- **The depletions do not simply divide by `SB0`** — under `HCR = 0 & msmMode > 0` they divide
  by biomass in the last projection year, so blanking them alongside a placeholder `SB0`
  discards a valid series.
- **A grep for `REPORT(` over `ceattle.cpp` over-counts** — several sit behind comments. A fit
  reports 101 quantities at 5.41.0 (99 at 5.24.0); enumerate from `names(fit$quantities)`. `quantity_dictionary()` is the
  registry, and `test-schema-quantity-dictionary.R` holds the two together.
- **`retrospective(getsd = TRUE)` can drop peels `getsd = FALSE` keeps** — the non-PD Hessian
  check only runs when an `sdreport` exists — so Mohn's rho can differ between the two.
- **`unweighted_jnll_comp` is written for 5 of its 21 rows** — composition, CAAL, stomach and the
  two linkage rows. Everything else is structurally zero there, not small.
- **`fit_mod(d, config = cfg)` replaces `d$model_config` with the config's** — a config from
  `run_config(model_config(), ...)` silently drops every linkage on `d` (57 REs → 0). Build it
  with `run_config(d, ...)`. `random_sel` never gates linkage REs.
- **`bias_adjust_proc` centres the lognormal priors and the recruitment deviations together**
  (5.33.0) — `FALSE` gives median priors *and* uncentred deviations. Reproduce an old prior by
  shifting its input; the Ianelli penalty with centred deviations cannot be reproduced. With the
  flag on, in the penalty years the pair favours `exp(-sigma_R^2/4)` times the pre-5.33.0
  recruitment.
- **A `data_list` element with no `write_data()`/`read_data()` support round-trips to nothing** —
  this is how `index_cov` was lost.
- **A `Comp_weights` of 1 under a Dirichlet-multinomial is a starting weight of e** — that
  likelihood reads the column as a log.
- **A Pearson residual divides by the effective sample size the likelihood used, not the input
  N** — `Comp_weights` multiplies the multinomial log-likelihood, so it *is* an effective sample
  size (the comp `SIMULATE` draw is `rmultinom_rce(n_nom * comp_weights(flt), p)`), and a DM is
  overdispersed by `(n + conc)/(1 + conc)`. Comp, CAAL and diet each have their own switch,
  weight and DM parameter block, and three *different* alpha constructions — only diet's is the
  clean `p·N·theta`. `.rce_comp_pearson()` is the one place that resolves this.
- **`fit_mod(estimateMode=)`** takes a string or the integer behind it: `"Estimate"` (0) =
  hindcast + HCR projection, `"Hindcast"` (1) = hindcast only, `"Projection"` (2) =
  projection-only from `inits`, `"DebugBuild"` (3) = build without optimizing,
  `"DebugOptimize"` (4). Prefer the strings.
- **`fit_mod(estimateMode = 4)` returns a placeholder objective** (`dummy*dummy`), because
  `build_map()` maps out every hindcast parameter. Don't read anything into a mode-4 objective,
  gradient, or Hessian. **Mode 3 returns the real objective**, so `obj$fn()` / `obj$gr()` are
  usable for diagnosing a model before fitting it — the analogue of WHAM's
  `fit_wham(do.fit = FALSE)` and SAM's `sam.fit(run = FALSE)`.
- **`fit$obj` is the PROJECTION's under any HCR but `NoFishing`; `fit$sdrep` too unless the HCR
  is also `ConstantF`** — at the default `estimateMode = "Estimate"`, `build_hcr_map()` maps
  every hindcast parameter off, so `obj$par` is `log_Ftarget`/`log_Flimit` alone (2 against the
  hindcast's 584 on `Atka2022`) and every delta-method SE of a hindcast quantity is exactly 0.
  Anything indexing `obj$par` by position must verify it against the vector it is labelling;
  `fit$identified` and `fit$.conv_hindcast` are the hindcast's.
- **`fit$data_list` is the PRE-`rearrange_data()` list** — it carries no `flt_sel_lead`,
  `flt_sel_type` or any other `rearrange_data()` output. Those live on
  `data_list_reorganized`, which `fit_mod()` does not keep. Recompute from `fleet_control`.
- **`Bin_first_selected` is a 1-based bin ordinal; `Sel_norm_bin` is an absolute age** — opposite
  conventions in adjacent columns, and `minage = 1` hides it.
- **`init_dev`'s ages start at `minage + 1`** — age `minage` in the first year is recruitment, so
  a `nages-1` axis is shifted, not just short. `.PAR_AXIS_OFFSET` (`R/0-parameter_index.R`) is the
  registry; `minage = 1` hides the shift everywhere bundled.
- **The conditioning check reads the correlation matrix, not the covariance** — `condition_number`
  changed meaning at 5.26.0; `covariance_condition_number` carries the old value, and the
  1e6/1e10 thresholds now fire less readily.
- **`fit_control()` bundles the optimizer and uncertainty knobs** — `getsd`, `bias.correct`,
  `loopnum`, `newtonsteps`, `getJointPrecision`, `nlminb_control`, and the bias-adjustment
  flags. `getsd = FALSE` leaves `sdrep` NULL, so `vcov()` returns NULL and uncertainty bands
  are NA. The refit diagnostics forward `phase`, `getsd`, the bias-adjust flags and
  `projection_uncertainty` — the last two are read back off `data_list`, because a freshly
  built `fit_control()` would silently reset them to its own defaults.
- **`run_mse()` pins the OM's stock-recruit and suitability windows to the pristine `om$`**, not
  the advancing `om_use$`, so the hindcast does not drift through the projection — essential for
  multispecies, whose predation suitability must stay fixed.
  `tools/verify/verify-mse-hindcast-invariant.R` checks it.
- **Every observation and process error is drawn in a `SIMULATE{}` block beside the density that
  scores it.** `sim_mod()` implements no observation model in R — it calls `obj$simulate()` once
  and writes the result back, so a new likelihood family owes a draw. Draw what the density
  assumes (bias-correction convention and scale included), REPORT under a `*_sim` name, and
  don't draw what the model does not define. `tools/verify/verify-sim-*.R` is the net.
- **An MSE draw is per observation row, so the row count is part of the RNG stream.** Anything
  changing the operating model's horizon or row count changes every draw after it — that is how a
  refit horizon set by the *next* assessment year made observation error depend on the assessment
  schedule (2.1% on a year whose advice was identical by construction). Before touching the
  horizon, the row count, or the draw order, ask what it does to a comparison of two schedules,
  not just to one run's reproducibility. Common random numbers are still incomplete between
  assessments; `inst/dev/TODO-mse-horizon.md` has the design and the measured numbers.
- **The guards are not themselves guarded.** `test-golden-regression.R` is `skip_on_cran()` AND
  `skip_on_covr()`, so until 5.16.0 it ran in no CI job at all; the `deep-checks` workflow now
  runs it nightly and asserts it produced assertions. `NOT_CRAN=false` must be step-level `env:`
  on `check-r-package`, never `$GITHUB_ENV` — through `$GITHUB_ENV` it did not hold reliably, and
  when it slipped Windows died with `0xC0000005`.
- **An access violation is memory corruption, not a bad optimum.** The model builds
  `safebounds = FALSE`, so an out-of-range access writes silently into adjacent memory. Build
  `RCEATTLE_SAFEBOUNDS=true` and run `tools/verify/verify-safebounds.R`, which asserts
  `-DTMB_SAFEBOUNDS` actually reached the compile line — `pkgload` only recompiles when sources
  change, so a clean result against a stale `.so` means nothing.
- **A slow fit is the model, not a regression** — `BS2017SS` has needed ~500–700 `nlminb`
  iterations since at least 2023.
