# Verified traps

Each entry is something that has actually bitten, with the evidence. `CLAUDE.md` carries a
one-line version of the ones you need before touching the code; this file holds the detail and
the measured numbers so a claim can be checked rather than believed.

Verified against source, 2026-08.

## Build and test

**A header-only edit under `src/TMB/` can leave `pkgload::load_all(".")` loading the OLD binary,
with nothing printed amiss.** Observed three times: a one-line `selectivity.hpp` fix left the
defect it fixed still reproducing until `ceattle.cpp` was touched; `growth.hpp` edited at
21:18:18 ran against a `ceattle.o` of 21:02:45; and reviewing PR #263 (2026-10-21) the stable
descending-limb fix read as FAILING its own new test, 4 of 4, across two reload attempts — the
second after deleting `src/TMB/ceattle.so` and `src/TMB/ceattle.o`, which is not enough because
the loaded library is `src/ceattle.so`. It passed 4 of 4, matching the closed form to 4.8e-25 at
age 15, once a reload rebuilt the `.o`. **The failure direction is not fixed**: that reading
accused a correct fix, where the first observation above excused a broken one. **It silently voids a before/after comparison** — a
revert-and-retest ran the same binary twice and showed "no difference", which reads as a test that
does not discriminate when the test was never given the reverted code. Anything measured after a
header-only edit is suspect until the rebuild is confirmed.

**The recipe: delete `src/TMB/ceattle.o`, `src/TMB/ceattle.so`, `src/ceattle.so` and
`src/Rceattle.so`, then check the `ceattle.o` mtime is newer than the header.** That is what
`src/Makevars`' own `clean:` removes, and a missing `.o` is always remade. Verified on the growth
anchor gate (5.57.3): one test file went 10 failures to 0 across two forced rebuilds, each `.o`
about a minute newer than the header, with no `touch` of `ceattle.cpp`. Treat a rebuild that
finishes in seconds as proof nothing compiled — a real one is ~60–140 s.

**Check the `.o`, never an `.so`.** `src/TMB/compile.R` copies `ceattle.so` up to `../`
unconditionally, so `src/ceattle.so` gets a fresh mtime whether or not anything was compiled.

**Why it happens at all is not established**, so do not reason from a mechanism. `TMB::compile()`
holds no mtime logic and `src/TMB/compile.R` declares no header dependencies, which would explain
it — but on pkgbuild 1.4.8 `pkgbuild:::headers()` globs `\\.h.*$` recursively and so does see
`src/TMB/*.hpp`, which should make `needs_clean()` true and force `--preclean`. Both observations
above are real; the chain between them is not understood. Run the recipe regardless.

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

**Compiling a header standalone without `-DTMBAD_FRAMEWORK` tests a translation unit the package
never builds.** `src/TMB/compile.R` passes `framework = "TMBad"`, so the shipped model is built
with TMBad, which undefines `CPPAD_FRAMEWORK`, `KEEP_COL`, `KEEP_ROW` and `eval` relative to
CppAD. A quick syntax or instantiation check of one `.hpp` that omits the flag therefore exercises
a different AD scalar type from the one the model uses, and whatever it reports — including a
"bug" — is about code that is never compiled. Pass the flag, or check through a real
`pkgload::load_all(".")`.

**Reading a randomized cross-check's bare correlation is not reading agreement.** `compResidual`
(Trijoulet et al. 2023) draws its composition residuals randomly, so the best score any method
could get against it is its own seed-to-seed correlation. Report the RATIO to that ceiling.
Measured on `BS2017SS` 2026-10-12 by `tools/verify/verify-osa-compresidual.R`: fleet 2 scores
0.648, which looks poor until the ceiling there is 0.661 — a ratio of 0.98 — and the ceiling is
low only because 41% of that fleet's cells hold less than one fish. Fleets 1 and 4 sit at 1.00.
Fleet 3 returns NaN from `resMulti()`, an instability on their side rather than a disagreement.
The bare number produced a false "fleet 2 diverges" on the first run.

## Shared parameter blocks

**Fleets sharing a `Selectivity_index` / `Catchability_index` share ONE parameter block.**
(`Q_index` is a deprecated alias of `Catchability_index`; the schema upgrades it on read.)
`adjust_map_shared_params()` copies the donor fleet's map slice over the rest of the group, so
any per-fleet setting that differs within a group is silently taken from the donor
(`Sel_start_year`, `Bin_first_selected`, `N_sel_bins`, …). The donor is the first *estimated*
fleet — an `Off` fleet's slice is all `NA` and must never lead. Penalties and priors on a shared
block are accumulated once, on the lead fleet (`flt_sel_lead` / `flt_q_lead`); without that gate
they are counted once per sharing fleet.

**A shared block STARTS at the geometric mean of its members' starting values, not at the
lead's.** TMB collapses a shared parameter to `tapply(par, map, mean)`, and these are all held
on the log scale, so no fleet keeps the value in its own row. It bites hardest on catchability,
where `index_log_q` is seeded from `log(Catchability_init)`: a shared q at the mean scales a
survey's whole predicted index by a constant factor, and no residual pattern distinguishes that
from a real change in abundance. Measured on an SS3 bridge for GOA Pacific cod, where a
converter-created fleet kept a default init and pulled the survey's q to
`sqrt(1.496398) = 1.223270` — 18% low across all 16 index observations, 6.23 nats, with the
standard deviations right to 5e-07 and both composition components agreeing to under 0.001.

**A selectivity or catchability LINKAGE on a shared block used to estimate one coefficient per
fleet, not one per group** (fixed 5.57.0). The base block is reconciled by
`adjust_map_shared_params()` across 15 by-fleet map slices; `beta_linkage` is a flat
`PARAMETER_VECTOR` with no fleet dimension and was never among them, so one design column on a
mirrored group became one FREE coefficient per member and the copies diverged. This was the
DEFAULT path: `by` defaults to `~ fleet` for both `sel` and `q`, and omitting `fleet =` does not
avoid it, because that argument is a filter whose `NULL` default means *every* fleet — the
shared-block spelling is an explicit `by = NULL`. Measured on `GOA2018SS` fleets 9 and 10, two
`Logistic` surveys sharing both indices: a `~ cut(Year, 2)` column gave coefficients -0.232615
and -0.352109 and selectivity differing by 0.123 at its worst age. On the SS3-bridged GOA Pacific
cod model, 18 paired `Srv`/`Srv_ae1` coefficients differed by up to 20.89 on the log scale, with
`index_log_q` non-identifiable and a non-invertible Hessian; the fix took free parameters from 344
to 328 against SS3's 330 and cleared both. **A held donor holds the whole group** — the first fix
tied a follower only where the donor's map entry was non-NA, which left the group untied wherever
SS3 fixes a block replacement (`est_phase = 0`), so 2 of those 18 rows stayed divergent. The
collapse also means a refit of a saved fit starts the shared coefficient at the MEAN of the old
per-fleet values, as the entry above describes.

**NOTHING fitted by that configuration is pinnable -- not the objective, not the coefficient --
and the reason is NOT "a flat surface".** The mirrored-pair regression fit (`GOA2018SS` fleets 9
and 10, `phase = FALSE`) measures, same source and data throughout:

| where | objective | `beta_linkage` | iters | max&#124;grad&#124; |
|---|---|---|---|---|
| macOS arm64, cold | 12879.997040 | -0.240043 | 1,759 | 8.56e-4 |
| macOS arm64, 3 jittered starts | 12879.997040 | -0.240042/4 | | ~1e-3 |
| ubuntu CI, **2 runs** (39345727752, 39449526638) | 13807.6 | -0.240043 | | |
| ubuntu CI, **2 runs** (39624390458, 39710679542) | **13818.2** | **-0.216** | | |
| macOS, `phase = TRUE` | **12862.289542** | -0.240042 | **31** | 9.39e-4 |

**Four ubuntu runs gave exactly TWO answers, each reproduced twice**, two of them on *different*
runner images, on byte-identical model code (`git diff` over the intervening `dev` commits touches
only `.github/`, `tools/ci/`, `inst/dev/`, `.claude/` and test files), the same R 4.6.1 / TMB
1.9.25, and with **no OpenMP** (no `_OPENMP` or `parallel_accumulator` in `src/TMB/`), **no Laplace
block** (`random_rec = FALSE`), **no RNG** (`inits = NULL`) and no clock in the path. Deterministic
code cannot return two answers, so the earlier explanation in this entry -- that the optimizer's
arithmetic path decides where it lands -- **does not hold**. Linux has two reproducible attractors
~930 nats worse than macOS, and that is the anomaly to explain.

**The cold fit is also simply a bad fit**, which is a separate point worth keeping: 1,759
iterations to land 17.7 nats WORSE than a phased fit reaches in 31. Both stop near 1e-3 in max
gradient, so phasing wins on objective, not on convergence -- do not describe it as
"better converged".

**What the platforms actually differ in, measured 2026-10-22.** A probe on `ubuntu-latest`
(`tools/ci/machine-fingerprint.sh` now records this on every golden and coverage run):

| | macOS arm64 (local) | ubuntu-latest |
|---|---|---|
| CPU | Apple M3 Pro | AMD EPYC 9V45, AVX-512 |
| compiler | Apple clang 21.0.0 | gcc 13.3.0 |
| R / TMB | 4.5.1 / 1.9.21 | 4.6.1 / 1.9.25 |
| BLAS | `libRblas` (reference, **single-threaded**) | `libopenblasp-r0.3.26.so` (**pthread**) |

Three things that probe established. **The model is deterministic per machine**: the same fit twice
in ONE process is bitwise identical in objective, coefficient and iteration count -- so this is
arithmetic between machines, NOT memory corruption, and not something non-deterministic inside
TMB. **There are at least THREE Linux answers, not two**: that run gave objective 13817.0067341985
with `beta_linkage` **-0.030399468** after **13,581** iterations, against 13807.6/-0.240043 and
13818.2/-0.216 from the coverage shards. The coefficient therefore ranges 8-fold. **And the
landing point is locally stable on macOS**: perturbing the converged parameters by a relative
1e-14, 1e-12 and 1e-10 returns to 12879.997040 / -0.240043 in 5 iterations every time, so the
macOS point is a real attractor rather than wherever the optimizer stopped.

**BLAS is the LEAST likely culprit, despite being the loudest difference in that table.** This
template does not route Eigen through an external BLAS -- no `EIGEN_USE_BLAS` anywhere in `src/`,
and one BLAS-level operation in `ceattle.cpp`. Its arithmetic is transcendental: **246 `exp()`,
105 `log()`, 55 `pow()`, 11 `pnorm()`** across `src/TMB/`. The better candidate is therefore
**glibc's libm, which dispatches `exp`/`log`/`pow` by CPU feature through ifunc resolvers** -- a
different kernel on a different runner is a different last bit, at every one of those call sites,
many of which sit inside the species-age-year loops. That is a candidate, not a conclusion; an earlier draft of this entry asserted BLAS
and was wrong, which is the second mechanism this entry has had to withdraw.

**None of it is a defect in TMB, glibc or OpenBLAS.** Floating-point addition is not associative
and libm implementations are not required to agree bit-for-bit across platforms; both are working
as specified. **What makes it matter is local: the fit does not converge**, so there is no optimum
to absorb the difference. 1,759 iterations on macOS and 13,581 on Linux, both stopping near 1e-3
in max gradient, while a phased fit reaches a BETTER objective in **31**. Fix the convergence, not
the arithmetic.

**What has been ruled out, and what has not.** A bounds-checked build
(`RCEATTLE_SAFEBOUNDS=true`) driven over this exact configuration reports **no violation** on
macOS and reproduces the ordinary build's numbers to all digits, so an out-of-range Eigen access
is not the explanation there -- but the divergence is on Linux, and a bounds-checked run on
`ubuntu` is the outstanding test (note safebounds catches out-of-RANGE accesses, not uninitialized
reads). `TMBhelper` is not the explanation: the **coverage** shard logs -- which is the workflow
that produced these numbers, not deep-checks -- show `TMBhelper 1.4.0` in sessioninfo, and local
is 1.4.0 too, so `.fit_tmb()` (`R/0-tmb_helpers.R:21`) dispatches the same way on both. That does
NOT clear the macOS-vs-Linux half, where local TMB 1.9.21 / R 4.5.1 / arm64 faces CI's 1.9.25 /
4.6.1 / x86-64; it clears only the Linux-vs-Linux anomaly, which is the one that cannot be
arithmetic. **Do not reason from a mechanism.**

**Pin history, corrected.** The objective literal and the coefficient pin shipped **together** in
`54c3ff32` (5.57.0, 2026-10-09), before any Linux evidence existed -- an earlier draft of this
entry said they were added in sequence, each in response to the last failing, and only the
ceiling was. All three were then falsified: the literal and the ceiling by the objective, the
coefficient by run 3. The file now asserts only the structural facts -- one estimated level, two
identical curves -- which held on Linux in **all four** runs.

`../Rceattle-dev/CLEANUP_BACKLOG.md` carries this as an open reproducibility defect rather than a
test nuisance: a user can run this configuration and get a 10% different coefficient by machine,
and `convergence` does not flag a max gradient of 1e-3. **Do not pin any fitted number from an
unphased fit that stops this flat** -- and note it failed invisibly on Linux from 5.57.0 (on `dev`
2026-10-10 UTC) until 2026-10-20, because `tools/ci/coverage-shard.R` discarded the shard result
(see "The guards are not themselves guarded"). It surfaced on the first run that read one.

**One coefficient is not one curve: the linkage OFFSET is per fleet.** The base parameter is
shared, but the offset the coefficient scales accumulates into `inf_offset(param, flt, sex, yr)`
(`linkage.hpp`), indexed by fleet — so a row naming one member of a mirrored block moves that
member alone. On the same pair, fitted: naming both fleets gives one coefficient and selectivity
identical to 0; naming the lead alone gives one coefficient and 0.265 of divergence; naming the
follower alone, 0.248. **So "put it on the lead fleet" — right only for a PRIOR, which is
re-targeted onto the shared base parameter — is WRONG for a design column, and the pre-5.57.0
refusal messages said it anyway.** That includes the apical height: every param code 0–11 and the
q offset write a per-fleet slot, `apical_offset(flt, sex, yr)` among them, so the apical refusal's
"the fleets sharing the index inherit it" was false and is corrected at 5.57.0.
Partial coverage is refused since 5.57.0, along with three others
`.stop_if_mirrored_block_linkage()` carries: a random-effect linkage on a mirrored block (each
named fleet gets its own deviation series and its own SD — 84 slots in two sigma groups of 42 with
2 free `log_sigma_linkage` levels, where the lead alone gives 42 in one group with 1; they cannot
be map-tied, because separate sigma groups would each score the same 42 deviations and fit two SDs
to identical data), and two specs for one parameter that disagree on `link`, `bounds` or `init`
(merged, the surviving bound was whichever row sat FIRST in the table — `match()`'s
first-occurrence rule in `fit_mod()`'s bounds reduction, so writing the follower's spec first
silently won — and two `init` values started at their mean, 2.0 and -2.0 beginning at 0).
`est_phase` is deliberately exempt: a held donor holding the group is the rule, and the SS3 bridge
needs it.
`build_map()` warns since 5.44.0 (`.warn_shared_block_start()`), for both deviation sds and the
catchability. Two edges: a member whose `Catchability_init` is blank, zero or negative seeds the
WHOLE group at `NA`/`-Inf` and it cannot fit — `data_check()` requires that column positive only
on fleets carrying index rows, and exempts `Analytical`/`AnalyticalArith` — and under
`Estimated-with-prior` the prior centre, and the prior SD, are the **lead's** alone, so a
non-lead fleet's values there are simply never read. **Inject or set per-fleet starting values by
BLOCK, never by fleet name.**

**A linkage intercept `init` / `bounds` is injected per block since 5.58.0, and was per FLEET
before.** `.stop_if_shared_block()` refuses a value set on a follower and names the donor, so the
spelling it recommends was the diluted one: the push wrote the donor's cell and the rest of the
group kept its build defaults, so the level started at their mean. Measured on `BS2017SS` with
fleets 4 and 5 sharing a `Selectivity_index` and an intercept `init` of 0.5 on `slp_asc`, the
level started at **-0.096574** where `log(0.5)` is -0.693147 — a slope of 0.908 for a requested
0.5. The bounds half was worse, because `fit_mod()` takes one bound per level by first occurrence,
so a donor that was not first in array order left an unwritten default as the surviving bound.
Both pushes now write every member, `.shared_block_members()` being the expansion.

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

**A NON-INTEGER `growth_age_L1` made Jan-1 length-at-age year-varying.** Fixed in 5.57.3; the
reason it hid is the trap. The Jan-1 cohort recursion advances a slot from the same cohort one year
younger, so it is valid only once that predecessor is on the growth curve. The gate testing for
"not yet" compared the slot's age to `max(minage, 1)`, where the first slot above the anchor is at
`max(minage, floor(age_L1) + 1)`.

Those agree for every anchor anyone uses, by two different mechanisms — which is why no fit was
wrong and no test fired. At an INTEGER anchor the linear ramp ends exactly at `l1`, so the
recursion inherits a curve value and is right regardless of the gate. At a sub-year anchor such as
0.5 the first slot above it is Jan-1 age 1, which `max(minage, 1)` names exactly when
`minage <= 1`; above that, slot 0 is already past the anchor and the old gate caught it because
`minage == max(minage, 1)`. Measured over all 188 SS3 control files reachable from this machine
(values 1 x113, 0 x44, 0.5 x30, 3 x1), the 5.57.3 gate reroutes **nothing** at any `minage` 0-4 or
`nages` 2-20.

What it cost where it did fire: on a `make_test_data` fit at `minage = 1` with anchor 2.5, the slot
above the anchor came back **2.2959 cm against a correct 1.6965 (+35.3%)**, with a 0.5994 cm spread
across years that constant growth parameters forbid. Transcribed at realistic parameters
(`Lmin = 2, l1 = 10, Linf = 60, K = 0.25`, `minage` 0-4 x anchor 0.25-5.00, `nages` 8, `nyrs` 6)
the worst error is **+59.8% / 8.19 cm** on Richards and +54.1% / 7.05 cm on von Bertalanffy; both
grow with the series length, so quote the grid with any such figure. It propagates into every
fleet's weight-at-age through `id_pop`, and it put the two growth routines in direct contradiction:
the same cohort at the same instant read 15.8752 from the within-year advance and 19.8139 from
Jan-1, a 3.939 cm gap that is now exactly 0.

**The lesson is about the test, not the gate.** `test-growth-vonbert-year-invariance.R` pinned
exactly this invariant and could never fire, because it used the default anchor and every default
is an integer or 0.5. **A year-invariance test is only as good as its parameter grid**: pin a
convention at a value where the two readings differ, not at the default.

**A joint-sex composition is predicted in TWO loops, and the row is normalized across both.**
`ceattle.cpp` builds a `flt_sex == 3` composition as females over the first `nlengths` (or
`nages`) bins and males over the next block, in a second loop under `if(flt_sex == 3)`, then
divides the whole row by one shared sum. So changing how a composition is predicted in one loop
and not the other is **worse than changing neither**: the normalization carries the
inconsistency into the half that was corrected, and the row still sums to 1, so nothing
complains. That is how the 5.48.0 length-comp ageing-error fix shipped half-applied through a
review, a green suite and a golden re-pin — and every length-comp row in `GOAatf`,
`GOAatf2023` and `GOA2018SS` species 2 is `Sex == 3`, so the joint loop was the ONLY one those
fits exercise. Grep both loops, and pin the pair with an invariance test rather than a value:
`test-likelihood-length-comp-ageing-error.R` asserts a length comp is unchanged by scrambling
`age_error`, with a finiteness guard, because an alternative model that comes back `NaN`
differs from the base everywhere and passes a difference check for free.

**A map built from a SKELETON model freezes both double-normal ends before the real values land.**
`build_params()` defaults `sel_dn6[5:6]` -- `start_logit` and `end_logit` -- to the sentinel
**-999**, and `build_map()` masks those two slots on that VALUE (`R/3-build_map.R`, the
DoubleNormalSS3 branch). Both are correct in isolation: an unscaled end has no meaningful default,
and masking a sentinel is right. The hazard is the build ORDER when a caller stages it. A common
bridging pattern builds one model with `inits = NULL` purely to get the parameter skeleton, derives
a map from it, then passes that map to the real fit alongside real starting values — and the map
still carries the mask, so the ends are held at values the caller did supply. Nothing warns: the
parameters hold the right numbers and simply never move.

Measured on `make_test_data()` with both fleets on `Selectivity = "DoubleNormalSS3"` and
`Selectivity_dimension = "Age"`, `estimateMode = 3` throughout. The skeleton map frees **8 of 12**
`sel_dn6` entries and masks **both** `end_logit` cells. Injecting `inits$sel_dn6[6, , 1] <- -0.5`
and building twice:

| map | `end_logit` value | map entry | free `sel_dn6` |
|---|---|---|---|
| reused from the skeleton | -0.5 | `NA` | **8** |
| rebuilt from the same `inits` | -0.5 | 5 | **10** |

Same model, same starting values; the only difference is where the map came from. **Rebuild the
map from the inits you intend to fit, or inject `sel_dn6` before building the skeleton.** The
fixture above is the durable instance: it needs no consumer, no bridge and no SS3 report, and it
holds for as long as `build_params()` defaults the ends to the sentinel.

**What it looks like from the outside, which is why it costs a day to find.** The parameters hold
exactly the values the caller asked for, the forward pass matches, the fit converges, and the only
symptom is a free-parameter count two short of the reference model. Every mechanism a reader
suspects first is individually correct — the per-cell mask on the sentinel value,
`.rce_dn6_ends()`, `adjust_map_shared_params()`, `map_linkage_adjuster()` — so the search goes
through all of them before reaching the build order. Found in the SS3-bridged GOA Pacific cod
work, on two `end_logit`s that SS3 estimated and Rceattle held; stated here as the behaviour of a
staged caller rather than as a property of that model, because once such a caller rebuilds its map
the instance goes away and the trap does not. Note this is a different code path from the
group-level sentinel patch, which keys on `Selectivity_index` rather than on the cell's value.

**A linkage parameter has THREE registries, and the rule-12 pair is only two of them.**
`LINKAGE_PARAM_CODES` (`R/0-linkage_encode.R`) must match `linkage.hpp`, which is rule 12 — but
`build_srr()` keeps its own whitelist in `RECRUITMENT_LINKAGE_PARAMS` (`R/0-build_srr.R`), and
`.REC_PARAM_TO_INDEX` beside it maps a linkage param to its `rec_pars` column. Add a parameter to
the rule-12 pair alone and `build_srr()` rejects it with "unknown recruitment linkage
parameter(s)"; the other processes have the equivalent (`.SEL_PARAM_TO_SLOT`). This is the same
shape as `Index_distribution`'s second registry. Grep every registry before adding one.

**Every `(Intercept)` linkage coefficient but recruitment's `R_init` is pinned at `NA`, on the premise that the process's
own base parameter carries the level.** `R/3-build_map.R` masks `beta_linkage` for intercept rows,
`build_params` forces their starting value to 0, and `build_parameter_bounds` loosens their bound
to ±Inf and propagates the caller's bound to the base parameter instead. That premise held for
all six processes until recruitment `R_init` (5.48.0), which multiplies the initial age-structure
and has NO base parameter — so `~ 1` produced zero estimable `beta_linkage` entries and moved
nothing, while every builder reported success. The R sites now share `.is_pinned_intercept()` and
its complement `.is_level_intercept()` (`R/0-linkage_encode.R`), both keyed off
`.is_r_init_linkage_row()`. **Don't trust a count — this one was published as "four", then
"five", then "seven", then "eleven", and every published figure has been wrong within a release.
So the grep is the authority, not any number here, and not a line list: the line numbers in the
earlier version of this entry were all stale within two commits.** Grep
`.is_pinned_intercept|.is_level_intercept|.is_r_init_linkage_row|.REC_PARAM_TO_INDEX\[` and read
every hit. At 5.49.0 that is 13 conditions in R — `R/2-build_params.R` (5),
`R/4-build_parameter_bounds.R` (4), `R/3-build_map.R` (2), `R/1-data_check.R` (1),
`R/6-fit_mod.R` (1) — plus the definitions in `R/0-linkage_encode.R` (6) and a comment in
`R/0-build_srr.R`, which are hits but not conditions, and one guard in `src/TMB/ceattle.cpp`
(`param != RCEATTLE_REC_R_INIT`, next entry). Three of them read the `.REC_PARAM_TO_INDEX` `NA`
instead of the predicate: that is correct in `map_linkage_adjuster()` and in `build_bounds()`'s
base-parameter push (there is no base parameter to mask or to bound), and was a defect
in `.push_linkage_intercept_inits()`, which is where `fit_mod()` re-applies an `est_phase = 0`
intercept over supplied `inits` so a FIXED value beats a warm start. `R_init` is the one
parameter whose fixed value lives in `beta_linkage`, so it was the one fixed level a warm start
silently overwrote — and the row is mapped `NA`, so the wrong value was then held for the whole
fit. The golden recipe and the GOA cod bridge both warm-start, so that is the normal path, not an
exotic one. Note too that a linkage TABLE row is not an estimated parameter: `print()` on the
table said "1 coefficient(s)" for the inert case. Count `beta_linkage` in
`obj$env$last.par.best` to tell the difference.

**The C++ linkage-prior block re-targets an `(Intercept)` prior onto the base parameter, and will
read off the end of the matrix for a parameter that has none.** `ceattle.cpp` slot 19 sets
`b = rec_pars(sp_idx, param)` for every recruitment intercept; `rec_pars` is `nspp x 3` and
`R_init` is code 3, so a prior on an `R_init` intercept read one element past a
`PARAMETER_MATRIX` and promoted adjacent heap memory to an AD variable — no crash guaranteed, no
R-side refusal. Guarded at 5.48.0 so the prior stays on `beta_linkage(i)`, which IS the level, so
the density reads on the multiplier (`lognormal(0, sd)` centres on no shift off R0 — the same
contract as `log_sel_apical`). The growth branch beside it already had the pattern
(`if (param < RCEATTLE_N_GROWTH_PARAMS)`). **A natural-scale `init`/`bounds` on such a parameter
has to be logged onto the coefficient too**, since there is no base parameter to log it onto;
read raw, `init = 0.5` ("half of R0") starts the level at `exp(0.5)` = 1.65x R0.

**A free initial recruitment level has a FLAT TAIL at the low end, and a few percent of fits
converge into it with a clean gradient.** `R_init` is informed only by the first year's
observations, so once the level is low enough that the initial cohorts are effectively
annihilated the data cannot distinguish one tiny level from another: at `beta = -30` the
objective is flat to `1.7e-08` while still finite. `tools/verify/verify-sim-recovery-r-init.R`
measures, at a true multiplier of 0.5 over 60 replicates on `make_test_data(nyrs = 40)`: a few
replicates run away to 5e-05 or below with `max|grad|` around 1e-04, so a gradient filter does
NOT find them and `convergence` reports success. Over the rest the level recovers with a low bias
of about 0.12 log units (z around -3.5, empirical sd 0.27, ~two thirds of replicates below
truth). **Don't quote those figures as constants** — the harness optimizes its own base fit, so
`ini`, the operating model and every draw inherit a machine-dependent starting vector, and the
runaway count has come out at both 3 and 6 of 60 on identical arguments. Re-measure.

The runaways are **not** a local optimum a multi-start escapes: four starting levels (truth, 0,
+1, -3) all return the same estimate to four decimals at the same objective, so the data really
do prefer the annihilated initial state and `jitter()` will not find it. The bias is a genuine
finite-sample property, not an artifact of the fixture's unconverged base — the score on the
level at the truth is centred (z = 0.5 over 60 replicates, implied first-order bias 0.003) —
but its MAGNITUDE belongs to this fixture, where the level is identified against four penalised
`init_dev` on one year's data. Don't read 0.12 as a universal correction.

**Bound the level on a real assessment**: `bounds` on this intercept are natural-scale, so
`bounds = list(intercept = c(0.05, 20))` keeps a fit out of the tail (`intercept` is the alias
that saves nesting backticks round the `(Intercept)` key). A lower bound of exactly 0 does NOT:
it logs to `-Inf`, which is the same unbounded floor, so it reads as protection and gives none.
Read the estimated multiplier, never just the convergence flag.

**`N_eq` is an equilibrium only if every age of it carries the same recruitment level.** Its one
consumer is `equil_catch_hat`, which integrates Baranov over all ages including age 0, so an
age-0 term left at `R_init` while ages 1+ carry a level makes the "equilibrium" sit at two levels
at once. Age 0 of `N_at_age` is deliberately NOT scaled — year-1 recruitment is a recruitment
year with its own deviate, not part of the pre-hindcast level — so the two arrays diverge here on
purpose. Live only under `initMode = 6`, the one mode that reads an equilibrium catch.

**A process linkage assigned straight onto a `data_list` is silently discarded.**
`fit_mod()` overwrites `data_list$srr_linkages` from `recFun` (`R/6-fit_mod.R:549`), and does the
same for the other processes from their own `build_*()` objects. So
`d$srr_linkages <- list(...)` followed by `fit_mod(d)` fits a model with no linkage at all and
says nothing. Build the spec through `build_srr(linkages = ...)`. Same shape as
`fit_mod(d, config = cfg)` dropping every linkage unless `cfg` came from `run_config()`.

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
reads `log_non_par_sel`, which `selectivity.hpp` (case 9) mean-centres per year. A shift common
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
`N_at_age_dB0` / `N_at_age_dBF` sit behind
comments. Enumerate from `names(fit$quantities)`, not from the source.

**And the count depends on `selectivity_se`, which is why three totals are on record
(99, 101, 104).** Measured on `BS2017SS` at `estimateMode = 3`: a default fit reports **102**
quantities, `fit_control(selectivity_se = TRUE)` reports **104**, and
`nrow(quantity_dictionary())` is **104** — the two it adds are `log_sel_at_age` and
`log_sel_at_age_index`, which exist only when that flag is on. Multispecies is also 102, so the
count is identical single- and multi-species; it is the flag, not the mode, that moves it. Quote
the figure with the `selectivity_se` it was counted under, or it will disagree with the next
person's.

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

**`Index_distribution` has TWO more hand-synced registries.** A family added to
`index_distribution_map` (`R/0-switches.R`) must also be classified in
`.index_rows_natural_scale()` in the same file, which is what `residuals(type = "pearson")`,
`plot_index()`'s observation interval and `plot_index_residual()` read to pick a scale. A
natural-scale family that misses it **does not error** — it silently gets the log-scale formula,
whose `sigma^2/2` is then a number the size of the index squared. Since 5.64.3 a family that
carries a COVARIANCE rather than a scalar sd goes in `.index_fleets_covariance()` too, which
those same three consumers read to take their sd from `sqrt(diag(index_cov))`; miss it and the
fleet is standardized by a `Log_sd` its likelihood never read (measured 10,180-18,439x wrong on
one EBS pollock BTS configuration). The MVN/MVNORM code set is additionally spelled out by hand
in `0-clean_data.R`, `1-data_check.R` and `5-rearrange_data.R`.
`test-index-natural-scale-downstream.R` enumerates the map and is the net — but note its
fixture passed `cov_sd` equal to `sd` until 5.64.3, so the right and wrong denominators were
equal by construction and nothing there could see the covariance defect.

**`nages` is a COUNT of age bins, not the oldest age.** Ages run
`minage .. minage + nages - 1`, and the array's 3rd dimension is a 1-based bin index, so age `a`
is at index `a - minage + 1`. `minage = 1` makes the two coincide — which is every bundled
dataset and the live assessments — so a bin/age confusion passes every test and every real
model, and only shows up on the next `minage != 1` species. Write
`seq_len(nages[sp]) - 1 + minage[sp]` and mean it.

**`jnll_comp` columns mean different things on different rows.** The matrix is dimensioned
`max(n_flt, nspp)`, and the column index is whatever loop variable wrote the cell: rows 1-8
(index, catch, composition, CAAL, non-parametric selectivity, selectivity deviates, catchability
prior, catchability deviates) are written inside `for(flt = 0; flt < n_flt; flt++)`, rows 9-20 by
`sp`/`rsp`/`slot_col`, "Linkage random effects" always into column 1, and row 22 ("Initial
equilibrium catch", added 5.46.0) back into a FLEET column -- so the fleet axis is not
contiguous and stopping an axis guess at row 20 is wrong. So `jnll_comp[3, 2]` is
fleet 2's composition likelihood while `jnll_comp[11, 2]` is species 2's recruitment deviates,
and `rowSums()` pools across two different axes. `.JNLL_ROW_AXIS` (`R/9-profile.R`) is the
registry; `test-schema-jnll-rows.R` parses every `jnll_comp(JNLL_*, col)` write in the template
and asserts each row's declared axis matches the column it is actually indexed by. Verified
2026-08-26: 124 writes, all 21 rows covered; 133 writes at 5.41.0, same 21 rows; 134 writes at
5.48.0 over 22 rows, the new one being the initial equilibrium catch.

**`unweighted_jnll_comp` is populated for 5 of its 22 rows.** It exists so Francis and
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

**And since 5.58.2 a KEPT peel can also have no `sdrep`, so `is.null(peel$sdrep)` no longer
implies `getsd = FALSE`.** Under `getsd = TRUE` each peel's standard errors are rebuilt from a
free-hindcast `estimateMode = 3` build, because the forecast refit's own `sdreport` reports
exactly 0 across the span the peel is drawn over. Where that rebuild throws, cannot factor its
Hessian, or reports nothing finite, the peel keeps its point estimates and carries no `sdrep` at
all — previously it carried the refit's zero-error object and drew a zero-width band as a measured
one. Read `attr(peel, "sd_report_error")`, or `print()`'s count, to tell this apart from
`getsd = FALSE`.

**The failure is deliberately NOT recorded on `.conv_hindcast$sd_present`, and must not be.**
`.refit_converged()` drops a peel on `sd_requested && !sd_present` (`R/0-convergence.R`), and the
forecast refit records both as TRUE — so writing `sd_present <- FALSE` to make the convergence
battery notice would convert a missing band into a **dropped peel** and change the denominator
Mohn's rho averages over. A reviewer proposed exactly that; it was measured as wrong before it
shipped.

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

**A bundled maturity above 1 is a data defect, not a sex-ratio convention, and on a TWO-SEX
species it doubles SSB.** `GOA2018SS` and `GOAcod` carry Pacific cod maturity = **2.0 at
Age1-Age12** — every age those models use (`nages[cod] = 12`). It looks harmless only because cod
is `nsex = 1` with `sex_ratio = 0.5` and `ceattle.cpp:1025` multiplies by the sex ratio **for
one-sex species only**; `mature_females` is otherwise a plain copy of `maturity`
(`ceattle.cpp:1015`), so a two-sex species takes the value at face value. Measured by doubling
two-sex Arrowtooth's row in `GOA2018SS` at the golden reference parameters: `ssb[endyr]`
1,044,496.91 -> 2,088,993.82, `SB0` 1,047,684.48 -> 2,095,368.97, `SPR0` 1.230019 -> 2.460037 —
all exactly x2 — with the objective **unchanged**.

Why it is a defect rather than a convention, measured across all 11 bundled datasets: 14 one-sex
species rows carry maturity <= 1 and exactly 2 carry 2.0, both the same stock. The control is
decisive — **BS2017SS/MS Cod** is the same species, also `nsex = 1`, also `sex_ratio = 0.5`, and
carries an ordinary ogive topping at **0.99728**. The row also mixes two scales: Age1-12 read a
flat 2, while the *unused* padding at Age13-20 reads 0.99828 -> 0.99965, an ogive tail on the
undoubled scale. GOA Pcod's real maturity is a **length** ogive (the SS3 model the bridge targets
sets `maturity_option: 1 = length logistic`), so the age sheet was never a proportion-at-age, and
`Rceattle-models/GOA cod/2024_pcod.R` already overrides it with `maturity[,-1] <- 1`.

Three reasons nothing caught it. `data_check()` has a lower bound only —
`R/1-data_check.R` reads `if(any(mat_vals < 0, na.rm = TRUE)){   # | mat_vals > 1`, the upper half
**written and commented out in the commit that introduced it** (`ab9872bf`) and never live.
`plot_maturity()` sets `coord_cartesian(ylim = c(0, 1.1))`, so the one figure that would show a
2.0 clips it out of the panel. And no test would move if the data were corrected: the five
`GOA2018SS` tests touching SSB or SPR either overwrite maturity first or recompute the expectation
from `d$maturity` itself. Across every workbook in the ecosystem, 70 of 528 maturity rows exceed
1; **all 70 are exactly 2.0**, every one a cod row in the GOA/AI lineage, including the live
`GOA-multispecies-assessment/data/raw/GOA_25_data_1977_2025.xlsx`.

**Documented, not fixed (Grant's call, 2026-10-12.)** Enabling the upper bound would refuse two
bundled datasets and 70 rows of live assessment workbooks, and the correct cod ogive has to be
recovered from the length-based model rather than guessed. Do not "fix" the data before that is
settled; the maturity note in the maintainer repository carries the full diagnosis.

**Mohn's rho has two live conventions differing by the PEEL COUNT, and r4ss returns both from one
function.** `r4ss::SSmohnsrho()` gives `$SSB` as `sum(rel_error)` over peels and
`$AFSC_Hurtado_SSB` as `sum(rel_error)/n`; SAM's `mohn.samset()` uses `colMeans`. **Rceattle is
the MEAN form** — it accumulates `rel_error` and divides by `N` (`R/9-retro_and_jitter.R`) — so an
Rceattle rho compares against r4ss's `AFSC_Hurtado_*` and SAM's output, **never against r4ss's
plain `$SSB`**, which is `n_peels` times larger. At the default `peels = 5` that turns a clean
0.04 into 0.20, exactly the band edge.

Verified 2026-10-12 on three peels of `BS2017SS` (terminal years 2014-2016, `estimateMode = 0`,
`phase = FALSE`), recomputing both forms from `$Rceattle_list` and comparing to `$mohns` at
`Forecast year == 0`:

| quantity | Rceattle (Pollock) | mean form | max abs diff | sum form / Rceattle |
|---|---|---|---|---|
| `ssb` | 0.08874500767 | 0.08874500767 | 3.5e-18 | exactly 3 |
| `biomass` | 0.058518029015 | 0.058518029015 | 1.4e-17 | exactly 3 |
| `R` | 0.5109766466 | 0.5109766466 | 5.6e-17 | exactly 3 |

Cod and Arrowtooth agree to the same precision, and the ratio is `n_peels` on every cell. The
absolute values are not an assessment result — `phase = FALSE` biases rho toward zero — only the
agreement between forms is being measured. Pinned by
`test-functions-retrospective.R` ("Mohn's rho is the mean over peels, not the sum"), which also
asserts the sum form is `n` times larger so a convention switch fails loudly.

**This is what makes the +/-0.2 band legitimate.** r4ss labels the mean form `AFSC_Hurtado`, and
the Hurtado-Ferro et al. (2015) ranges `print(retro, band = )` cites (-0.15 to 0.20 long-lived,
-0.22 to 0.30 short-lived) are on the mean — so Rceattle applies the band to the quantity those
ranges were derived for. Two further differences, neither a defect: r4ss also offers a
"Wood's Hole" all-years rho (averaging the relative error over `styr+1 .. peel endyr`, then summing
over peels) which **Rceattle does not implement** — measured 0.05167 against the terminal 0.08875
on the fit above; and SAM's `lagR` option changes only a column label, not a value.

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
scores it** (`ceattle.cpp` sections 5.2 and 5.13, plus one per likelihood slot; the
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
it — see issue #237 for how to get it back.

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
bundled dataset and the live assessments have `minage = 1`.

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

**A registry guard that builds its check-set WITH a regex and then floors it at
`expect_gt(length(.), 0)` cannot fail for the drift it exists to catch.** The floor is satisfied
by whatever survives, so a rename, a respelling or a reformat shrinks the set and the test still
passes — reporting green while measuring less. Measured 2026-10-12 by mutating a scratch copy of
`ceattle.cpp` and running the three real guard files verbatim against it: **9 of 15 mutations were
invisible**, every one of them a change the guard was written to catch, and in every case the pass
count was identical to pristine, so `tools/ci/source-guards.R`'s own per-block "asserted nothing"
check was satisfied too.

The sharpest three, all since closed:

* `test-dynamics-sex-index-bounds.R` anchored its pattern at `^\s*`, so
  `if(true) N_at_age(sp, 1, 0, 0) = Type(0.0);` and a respelled `N_at_age(sp,1,0,yr)` both left the
  set at 10 hits. Its array list was also hand-written and **missed `N_eq` entirely** — 4 of the
  real 15 sex-index-1 writes were outside the census permanently.
* `test-composition-age-hat-width.R`: **deleting one line** (the per-row widening
  `if(comp_row_cols > max_age_cols){ ... }`) reinstated the out-of-bounds write and the set merely
  fell 4 to 3.
* `test-dynamics-recruitment-minage.R`: rebinding a lag through a differently-named local took the
  set 4 to 3, and the site that left it stopped being checked at all.

**The fix is to pin the set by name and count, never to floor it.** `tools/ci/source-guards.R`'s
own `EXPECTED` vector is the in-repo pattern: a `setequal()` plus a message naming what appeared
and what disappeared. After pinning, 14 of the 15 mutations are caught; the one that is not
(removing a lag local's clamp, which reinstates the original negative-column read) is outside a
line-scoped grep's reach and the test now says so rather than implying it is covered.

**Two implementation traps, both measured.** Count per array with a SEPARATE `grepl()` each: R's
default engine is leftmost-longest, so one alternation capture attributes a line holding two of
these arrays to only the last — `N_eq(sp, 1, age) = N_at_age(sp, 1, age, 0);` (`ceattle.cpp:2019`)
is such a line, which is why **15 writes sit on 14 lines** and both numbers need pinning. And do
NOT strip comment lines before a positional lookback: the window is counted in lines, so removing
some lets a neighbouring block's guard drift into range and satisfy it.

**And on a pull request these three have no second line of defence.**
`.github/workflows/R-CMD-check.yaml` sets step-level `NOT_CRAN: "false"`, so every `skip_on_cran()`
fit block in the same files skips; `test-coverage` runs at `NOT_CRAN=true` and since PR #261
`tools/ci/coverage-shard.R` fails the shard on a failing test, but its shard contexts are not
among the six the ruleset requires, so that red advises rather than blocks;
`deep-checks` runs golden per-PR since 5.64.0 and the whole suite nightly, gating no merge.
`source-guards.yaml` is the cover, and since the ruleset below it is a required one. For the
sex-index and age-width guards the consequence of a gap is an unchecked write, which is not
observable from R at all — so the static pin is the only cover there is.

**`dev` and `main` are governed by a RULESET, which the branch-protection endpoint cannot see.**
`gh api repos/afsc-assessments/Rceattle/branches/dev/protection` answers 404 "Branch not
protected", which reads as ungated and is wrong; `gh api repos/.../rulesets` shows "Require the
guards before merging", active, requiring **six** checks by name: `registry and documentation
guards` plus all five R-CMD-check legs (`ubuntu-latest` release/devel/oldrel-1, `macos-latest
(release)`, `windows-latest (release)`). Two consequences. It sets `strict: true`, so a branch
**`BEHIND` its base cannot merge** and every PR pays another full cycle after anything lands under
it — `mergeable = MERGEABLE` reports conflict-freedom only, not ruleset satisfaction. And because
the Windows and macOS legs are required, the intermittent `0xC0000005` and a cancelled macOS runner
each block the merge until re-run. Renaming a required job un-gates it silently.

**Anything short of an ERROR passes R-CMD-check here, so a green check can be naming your
mistake.** `.github/workflows/R-CMD-check.yaml:127` sets `error-on: '"error"'` — one step looser
than the action's own `"warning"` default — so NOTEs *and* WARNINGs conclude SUCCESS. It is
deliberate (the comment above it cites a spurious rtools45 + GCC 14.3 `-Warray-bounds` inside SSE2
intrinsics), so don't tighten it without checking what that would reintroduce.

**The baseline to compare against, measured on a clean PR (#229): four runners `Status: OK`, and
`macos-latest (release)` at `Status: 1 WARNING`.** So any NOTE at all is new, and a second warning
would be new. Note the macOS warning is NOT the Windows one the setting's comment describes: it is
`checking whether package 'Rceattle' can be installed ... WARNING` over three
`-Wbitwise-instead-of-logical` diagnostics in `ceattle.cpp` (see `CLEANUP_BACKLOG.md`). `--as-cran`
is also expected to flag an installed-size NOTE once `data/` grows, per the `args:` comment, so
re-measure the baseline rather than assuming this one.

How it bit: PR #226 shipped a scratch test-runner at the repository root, and all five runners
printed `checking top-level files ... NOTE / Non-standard file/directory found at top level:
'suite226.R' / Status: 1 NOTE` while reporting pass. It was read as 12 of 12 green and merged, and
the file reached `dev`, where it shipped — in neither `.gitignore` nor `.Rbuildignore`. **Read the
`Status:` line against that baseline, not the check's conclusion**, and keep scratch runners out of
the worktree: write them to a temp directory, since `git add -A` cannot tell them from the change.

**A job that checks "did this guard measure anything" must check per BLOCK, not per file.**
Building `source-guards.yaml` (5.49.5) produced the defect twice in a row. First a heuristic
exempted any file whose `test_that` blocks were all followed within two lines by
`skip_on_cran()` — one line in the suite's dominant idiom would have exempted eight of the
sixteen guards. Replacing it with a per-file `passed > 0` floor had the identical hole: most of
these files hold several blocks, so gating one block kept the file's total high.
`test-schema-cpp-dispatch.R` dropped from 114 passing assertions to 81 with the job still green.
Both were found by adversarial review, not by the job. The per-block check then immediately found
a **third** instance already in the tree: that file's not-yet-implemented inventory — the pins
that non-parametric growth is stubbed and that Kinzey & Punt predation has no live `msmMode`
branch — sat behind a `skip_on_cran()` that existed only for the fixture half of the block.

**Pin a guard set by NAME, not by count.** `length(targets) == 16L` is satisfied by deleting a
real guard and adding any test that happens to match the discovery pattern. Only 5 of the 16 are
named anywhere else in the repo, so 11 could have been deleted silently.

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

**The same cause, a different trigger: a GIT BRANCH SWITCH does not invalidate the incremental
build either.** A cross-branch numeric comparison must delete `src/TMB/ceattle.o`,
`src/ceattle.so` and `src/Rceattle.so` and then **assert** the rebuild happened, or it runs one
branch's SOURCE against the other's COMPILED template. Measured 2026-10-14 on the rule-15 hake
MSE check: the `dev` leg reported objectives "identical to the last bit" against a predation
branch because it used that branch's DLL. The tell was a number that was too good -- the `dev`
leg reported `min suit_other` = 9.99999e-07, exactly the branch's floor value
`1e-6 / (1 + 1e-6)`, on a tree whose source contained no floor. After a forced clean rebuild the
same leg gave 5.16e-09 and `max M2` = 502,536 /yr, and the real difference between the branches
turned out to be 3.4e-09 to 4.3e-08 nats on the four multispecies fits.
`tools/verify/verify-safebounds.R` already does delete-then-assert; copy it rather than trusting
`pkgbuild::compile_dll()` to notice.

**And on Windows the harness handed the flag to its rebuild by a route that does not exist
there.** `system2(..., env = "RCEATTLE_SAFEBOUNDS=true")` serves only a command that takes
assignments on its command line -- `R` or `make`, not `Rscript` (`?system2`, Details) -- so the
`windows-latest` arm's first nightly run, 2026-10-16, compiled nothing in 0.1 s, captured an empty
build log, and the positive control refused the result. The guard working as intended: it reported
no measurement rather than a green run of none. The flag reaches the child through the
`Sys.setenv()` above the call instead, which both platforms inherit. The two routes are
indistinguishable on a Unix-alike, which is why review on Linux could not see it -- measured on
macOS with `env =` dropped, the flag off gives 81 lines of compile output and no
`-DTMB_SAFEBOUNDS`, the flag on gives 82 and the define.

**The first real `windows-latest` result, 2026-10-16: all 8 configurations clean.** The rebuild
took 2m39s where the broken call took 0.1s, and `test-selectivity-catchability.R` -- the file the
0xC0000005 died in -- ran 14.9s with no violation. Evidence, not clearance: `safebounds` range-
tests Eigen array access, so a clean run does not rule out the other ways that fault arises, and
the crash was intermittent to begin with.

**The golden regression ran in NO automated job.** `test-golden-regression.R` carries both
`skip_on_cran()` and `skip_on_covr()`. `R-CMD-check` sets `NOT_CRAN=false` deliberately (to keep
the multi-OS matrix fast), so the first fires; `test-coverage` runs under covr, so the second
does. The committed backstop for "a numeric change needs golden equivalence" therefore only ever
ran when someone typed `/golden-check` locally. Its measured cost of `0` in
`tools/ci/coverage-costs.tsv` is the tell — four phased fits do not cost nothing. The
`deep-checks` workflow runs it nightly and, since 5.64.0, on every pull request that touches
`src/**`, `R/**`, `data/**`, the golden test or its fixture — and **asserts the run produced
assertions**, by name, because a skipped file otherwise reports as a pass. That per-PR run is
still not one of the six contexts the ruleset requires, so a red golden advises rather than
blocks.

**And golden was only the worst case of it: no `skip_on_cran()` test could fail any job.** That
is 195 of the 271 test files, 95 of them carrying the skip at file level — every test that runs a
real `fit_mod()` optimization. `R-CMD-check`
sets `NOT_CRAN=false`, so it skips all of them. `test-coverage` does run them, and
`tools/ci/coverage-shard.R` still passes `stop_on_failure = FALSE` (deliberately: a failure must
not throw away the shard's coverage before it is written) — but since PR #261 it reads the result
afterwards and stops, so the "Upload test results" step now fires and the shard goes red. Before
that, a broken regression test printed into a shard log and the run went green; the first run
with the check on found `test-linkage-shared-block.R` had been failing on Linux for 11 days.
Found while reviewing 5.25.1, whose own PR body claimed a new test "fails if anyone restores the
old comparison"; it would not have. `deep-checks` runs the whole suite nightly with
`NOT_CRAN=true`, and `coverage-shard.R` asserts the mode in its own subprocess since PR #261 —
both **before** the fits rather than after, since `NOT_CRAN` going false presents as a fast green
run of nothing but the unguarded unit tests. **The nightly suite has two gaps of its own.**
First, the `schedule` fires on the DEFAULT branch, which is `main`, so `dev` never gets a
nightly suite at all — on `dev` the per-PR coverage shards are the only job that runs the
`skip_on_cran()` tests. Second, it has mostly stopped completing: measured over the last eight
scheduled runs to 2026-10-21, 10-13 cancelled, **10-14 success**, 10-15 cancelled, 10-16
**failure**, and 10-17 through 10-21 cancelled — five consecutive, and six of the last eight not
completing. A `cancelled` conclusion here is the 300-minute `timeout-minutes` killing the job,
not someone cancelling it: on run 39403748904 the `suite` job ran 09:51:40 to 14:52:40, exactly
301 minutes, while `golden` and both `safebounds` legs succeeded. **So its last completed run
was 2026-10-14, and the right reading of a green deep-checks badge is "golden and safebounds
passed", not "the suite passed".**
**Writing a regression test is not the same as having one run** — check which job executes it,
and whether that job can go red.

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

**The same fault also shows up as a broken worker message, not a crash code -- measured
2026-09-27 on #158.** The parent died in
`test_check() -> rs__read_message() -> rs__parse_header()` with `Internal callr error, invalid
message header`, and the Windows log carried **no `Failure (` line and no `FAIL n` count**: no
test failed, the IPC stream to a parallel worker broke. Treat it as the same dead-worker class as
the access violation, and check the same two things -- whether any assertion actually failed, and
whether the same commit passed elsewhere. Here **commit `44864a79` passed Windows in the
`pull_request` run and failed in the `push` run**, one commit with both outcomes on one platform,
which is the strongest evidence available that a Windows red is infrastructure. A red with a
`Failure (` line in it is a different thing and is worth investigating.

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

The committed gate (`test-golden-regression.R`) EVALUATES four reference models at stored
parameters rather than refitting them — it has not refitted since 5.45.1 — and pins the
objective, the gradient and the derived quantities below. The `/golden-check` recipe is a looser
local read at a 1e-6 absolute tolerance, so the two are not the same gate. Neither reaches:

- **The refit paths.** Eight entry points re-invoke `fit_mod()` through `.refit_like()`. Use
  `tools/verify/verify-refit-like.R`, which covers six of them —
  `sample_rec(update_model = TRUE)` and `reweight_comps()` are **not** in it.
- **The simulation draws.** `tools/verify/verify-sim-*.R`.
- **Any figure.** A plotter can change every number it draws and stay green. The net is
  `tests/testthat/test-plot-*.R` plus a before/after `ggplot_build()` diff at the merge base.
- **Anything the likelihood does not read.** `GOA2018SS` sets no `srr_fun` and no `srr_pred_fun`,
  the corner where SSB is pure derived output, so `maturity` never enters the objective there:
  every maturity perturbation gives `dObj = +0.0000e+00` on both `goa_ss` and `goa_ms` (positive
  controls on the same probe: an index observation x1.5 moves it +3.7664, a fitted catch row x1.5
  moves it +32.674). Correcting the cod maturity defect above would leave both objectives
  bit-unchanged. Maturity reaches the objective only where `srr_pred_fun > 0` or `srr_fun > 1`,
  through `R_hat` into `jnll_comp(JNLL_SRR_PENALTY)`. The gate therefore also pins `ssb`,
  `biomass`, `R`, single-species `SB0` and `B0`, `DynamicSB0`, `DynamicB0` and `M2_at_age` at the reference
  parameters (relative 1e-8), so a change to maturity, spawn timing or dynamic-B0 predation fails
  through those while the objective stays put.
  **Read the element count honestly**: 24,280 of the 32,976 asserted elements (73.6%) are exact
  zeros, leaving 8,696 informative numbers. `M2_at_age` is identically zero on BOTH
  single-species references, where it asserts "predation has not leaked into `msmMode = 0`"
  rather than pinning a predation mortality — the real M2 pin is the multispecies pair. A zero
  element is not a free pass, though: its relative band is exactly 0, so 0 to 1e-300 fails the
  element check (`expect_equal` alone would pass it).
  `SB0`/`B0` ARE stored for the multispecies pair but are not asserted — they are the 999
  MSSB0/MSB0 inputs there — so 876 stored elements go unchecked by design.
- **Harvest control rules.** All four references run `HCR = NoFishing`, so `proj_F` is 0,
  `Ftarget`/`Flimit` are the placeholder 1, and every fleet's `Proj_F_proportion` is 0, so `SBF`
  equals `SB0`. None of them is pinned. `test-dynamics-brps.R`, `test-switches-hcr-*.R`,
  `test-functions-fixed-n-hcr-outputs.R`, `test-functions-mse.R` and
  `test-mse-cap-and-hcr2-threshold.R` exercise the rules, against no golden value.
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

**The `golden` job's red on `main` at the 5.45.0 release was NOT the 52.9 basin — measured
2026-09-28, run 36433121293.** It failed on the MVN/Normal index-family configuration pin, at
the GRADIENT: `max|gradient|` 3.5e-04 on `ubuntu-latest` against a 1e-4 gate, where local
macOS/arm64 gives 2.5e-05. The objective it pins reproduced to its 1e-6 tolerance on the same
run, and the other 19 assertions passed, so nothing about the model moved. That configuration
sits on the q-versus-population-scale ridge, softened 11.5x by compound symmetry in `Sigma`, so
its gradient records where the optimizer stopped rather than what the model is; the gate is 1e-3
from 5.45.1. Two lessons: the four integration references were not the only platform-dependent
thing in that file, and `TMBhelper` **is** installed on the runner (1.4.0 from
`kaskr/TMB_contrib_R@d275e52`), so a different optimizer is not the explanation for any of it.

**`goa_ms` (fixed-M GOA multispecies) sits on a flat likelihood ridge:** the same objective at
different `par`/`ssb` across *different* code, though deterministic on same-code re-runs. Judge
it on `obj`/`jnll`, not `par`/`ssb`.

**`goa_ss` has a second local minimum of the negative log-likelihood 52.9 units higher, and a
one-ULP change in one gradient element is enough to send `nlminb` there.** Recorded on the
5.34.0 branch (PR #144, 2026-09-14): adding code the objective never evaluates left every
objective bit-identical and bounds-checked builds clean, but changed one `log_F` gradient
element by 3e-16 (presumably summation order), and from there `nlminb` reached the higher
minimum with `newtonsteps = 3` in place. A cold phased fit reached it again at 5.48.0
(12920.1030998153, `max|gradient|` 5.7e-11 — it polishes as well as the lower one), but **on a
half-applied length-comp fix, not on coherent code**: under the preceding release's behaviour a
cold fit reproduces the then-current reference exactly (12867.990267), and with both halves of that fix in
place it lands at 12866.845728, the same basin. An incoherent likelihood is itself a way into
the upper minimum. `goa_ms` warm-started from those upper-basin MLEs still reached its own
reference objective to ten decimals, so it does not always inherit `goa_ss`'s basin.
**`regenerate-golden-reference.R` therefore fits each reference from two starts — the recipe
that created it and the committed reference — and pins the lower**: the recipe alone can pin a
52.9-up minimum as a model change, and the reference alone would hide a change that moves the
optimizer on the cold path assessment scripts take. The
52.9 that `golden-check.md` attributes to tolerance-stopping (commit `1a172677`) is the same
gap; polishing did not remove it. A `goa_ss` delta of 52.9 with the other three models
bit-identical is this, not a numeric regression. Diagnose it from the gradient at the reference
`par`, not from the objective. **Fixed at 5.45.1**, not by making the optimizer land in one
basin but by taking it out of the gate: `test-golden-regression.R` now evaluates the likelihood
AT the committed reference parameters (`tests/testthat/fixtures/golden-reference.rds`,
`estimateMode = 3`) and asserts they are still a stationary point, with the cold-start fits kept
as `tools/verify/verify-golden-cold-start.R`. The second minimum is still there, and the gate no
longer sees it at all: a change that moves the OPTIMIZER into the other basin without touching
the likelihood leaves both the objective and the gradient at the reference exactly where they
were. That class -- phases, starting values, bounds inactive at the reference -- is what the
cold-start harness is for, and it is the price of a deterministic gate.

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

**AMAK conventions are reproduced deliberately, not corrected.** `src/TMB/selectivity.hpp` (case 11)
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
- **`jnll_comp` columns count fleets on rows 1–8 and row 22, species on rows 9–20**, so
  `rowSums()` pools across two different axes, and the fleet axis is NOT contiguous. Row 21
  (linkage random effects) is model-wide; row 22 is the initial equilibrium catch, by fleet.
  `.JNLL_ROW_AXIS` (`R/9-profile.R`) is the registry, so an axis guess that stops at row 20 is
  wrong.
- **A reference point CEATTLE never estimated is a number, not a gap** — `Ftarget`/`Flimit` sit
  at `exp(0) = 1` unless the HCR estimates them (gate on `build_hcr_map()`, never on `fit$map`),
  `SB0` under `msmMode > 0` is the 999 mt `MSSB0` placeholder until `MSSB0_derived` is TRUE, and
  the per-recruit quantities are zero outside `msmMode = 0`.
- **The depletions do not simply divide by `SB0`** — under `HCR = 0 & msmMode > 0` they divide
  by biomass in the last projection year, so blanking them alongside a placeholder `SB0`
  discards a valid series.
- **A grep for `REPORT(` over `ceattle.cpp` over-counts** — several sit behind comments. A default
  fit reports **102** quantities at 5.58.0 and **104** under `selectivity_se = TRUE`, which is the
  whole of `quantity_dictionary()`; enumerate from `names(fit$quantities)` and say which flag you
  counted under. `quantity_dictionary()` is the registry, and
  `test-schema-quantity-dictionary.R` holds the two together.
- **`retrospective(getsd = TRUE)` can drop peels `getsd = FALSE` keeps** — the non-PD Hessian
  check only runs when an `sdreport` exists — so Mohn's rho can differ between the two.
- **`unweighted_jnll_comp` is written for 5 of its 22 rows** — composition, CAAL, stomach and the
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
  are NA — but a NULL `sdrep` does not imply `getsd = FALSE` on a retrospective peel (above).
  The refit diagnostics forward `phase`, `getsd`, the bias-adjust flags and
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
  assessments; issue #237 has the design and the measured numbers.
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
- **`newtonsteps > 0` can return a parameter OUTSIDE its bounds.** `fit_mod()` passes
  `build_bounds()`'s range to `nlminb`, which respects it, but the Newton refinement afterwards
  is a plain unconstrained step in both `.fit_tmb()` paths — TMBhelper's own included, which the
  package cannot clamp. Measured on a GOA Pacific cod bridge at `newtonsteps = 3`: a growth
  parameter came back at 5.08e-05 against a lower bound of 1e-3, gradient 0.263, Hessian not
  invertible; at `newtonsteps = 0`, same code, it sat on the bound at 1.00e-03, gradient
  0.00246, `sdreport` clean. The three steps bought 0.0017 nats. `convergence` reports it as
  `parameters_outside_bounds` (FAIL) since 5.45.0; before that it read as `parameters_on_bounds`
  (WARN, "at a configured bound") because `par <= lo + tol` is also true below `lo`. The default
  is 0, but the golden references all run 3.
- **A slow fit is the model, not a regression** — `BS2017SS` has needed ~500–700 `nlminb`
  iterations since at least 2023.
- **The reproducible-install check passes on a stale install.** `inst/RELEASE-CHECKLIST.md`
  section 4 used `withr::with_temp_libpaths()`, which *prepends* a temporary library rather
  than isolating: the working library stays on `.libPaths()`. If the install fails — and
  `quiet = TRUE` hides that it did — `library(Rceattle)` loads whatever is already installed
  and `packageVersion()` reports that version. Caught 2026-09-25 verifying the 5.42.1 branch:
  the check reported **`packageVersion: 5.33.0`**, the copy sitting in
  `/Library/Frameworks/R.framework/.../library/Rceattle`, and printed `INSTALL CHECK OK`. A
  release could be signed off against a version it never installed. The recipe now installs
  into an explicit `lib`, reads the version out of that directory's `DESCRIPTION` rather than
  from a loaded namespace, asserts `getNamespaceInfo()$path` is under it, and drops `quiet`.
  Use `ref = "X.Y.Z"` rather than `"repo@X.Y.Z"` so a branch name containing a `/` survives.
