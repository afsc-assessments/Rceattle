---
description: Pick and run the tools/verify harnesses that actually cover the current diff
argument-hint: "[optional: a harness name, to run just that one]"
---

`/golden-check` fits four models and diffs the result. It does **not** reach the diagnostic refit
paths, the simulation draws, or any figure. `tools/verify/` holds the harnesses that do.

**They are not all the same shape.** Three are before/after digests you run twice; eleven are
self-contained checks you run once. Running a self-contained check on two commits and comparing
the two "PASS" lines verifies nothing — that is the failure mode this section exists to prevent.

## Before running

Work out which harnesses the diff needs, and say so before running anything.

### Before/after digests — run twice, compare

These take an output path and an optional `compare` argument.

| Harness | Covers |
|---|---|
| `verify-refit-like.R` | the `.refit_like()` refit path, via 6 of its callers |
| `verify-mse-repro.R` | fixed-seed multi-year MSE trajectory fingerprint |
| `verify-data-check-conditions.R` | every `data_check()` condition over the bundled datasets |

```
git worktree add /tmp/rce-base <merge-base>
export PATH=/usr/bin:$PATH
cd /tmp/rce-base && NOT_CRAN=true Rscript tools/verify/<harness>.R /tmp/before.rds
cd -            && NOT_CRAN=true Rscript tools/verify/<harness>.R /tmp/after.rds compare /tmp/before.rds
git worktree remove /tmp/rce-base --force
```

### Self-contained checks — run once, on your branch

These take **no path argument** and assert an invariant or a recovery property directly. Do not
pass them a `.rds` path.

`verify-mse-hindcast-invariant.R` (the MSE must not perturb the hindcast) ·
`verify-mse-schedule.R` (assessment schedules, year-indexed `catch_mult`, common random
numbers -- it takes NO path, so running it twice and diffing the PASS lines compares nothing) ·
`verify-sim-centering.R` · `verify-sim-diet.R` · `verify-sim-index-families.R` ·
`verify-sim-linkage-re.R` · `verify-sim-local-optimum.R` · `verify-sim-process-error.R` ·
`verify-sim-recovery.R` · `verify-sim-recovery-process.R`

`verify-rtmb-oracle.R` (the RTMB port against the C++) is self-contained too, but takes two
optional arguments that are not an `.rds`: a reference C++ build compiled without fused
multiply-adds, and `installed=<lib>` to load a byte-compiled install. Build the reference by
copying `src/TMB/*.hpp` and `ceattle.cpp` (renamed `ceattle_nofma.cpp`) into a scratch folder,
then `Rscript -e 'TMB::compile("ceattle_nofma.cpp", PKG_CXXFLAGS = "-ffp-contract=off", framework = "TMBad", safebounds = FALSE)'`.

```
export PATH=/usr/bin:$PATH && NOT_CRAN=true Rscript tools/verify/<harness>.R
```

Six harnesses take one argument that is a **count**, not a path. Passing a filename
yields `NA_integer_` in five of them, with a non-fatal "NAs introduced by coercion"
warning; `verify-osa-cdf.R` falls back to its default instead, which is quieter still: `verify-sim-recovery-M.R` (`nrep`, default 8),
`verify-sim-recovery-apical.R` (40), `verify-sim-recovery-np-integrable.R` (12),
`verify-sim-recovery-r-init.R` (20), `verify-sim-recovery-srr-msm.R` (30) and
`verify-osa-cdf.R` (20).

### Which to run

| If the diff touches | Run |
|---|---|
| `R/6-refit_like.R` or a caller | `verify-refit-like.R` (digest) |
| `run_mse()`, the OM/EM loop, projection horizons | `verify-mse-hindcast-invariant.R` + `verify-mse-schedule.R` (once), `verify-mse-repro.R` (digest) |
| a `SIMULATE{}` block, `sim_mod()`, `simulate.Rceattle()` | the `verify-sim-*.R` set (once each) |
| a `plot_*()` function | **none of these** — see below |
| `R/6-rtmb_*.R`, or C++ on a path the RTMB port covers | `verify-rtmb-oracle.R` (once per fixture: `fixture=bs2017ss`, and `fixture=bs2017ms` for predation), with a no-FMA C++ build |

## What nothing here covers

- **`sample_rec(update_model = TRUE)` and `reweight_comps()`.** Both refit through
  `.refit_like()`; neither is in any harness. A change to that helper needs them checked by hand.
  Say so rather than implying the harnesses were exhaustive.
- **Figures.** `/golden-check` and `test-golden-regression.R` diff the fit, so a plotter can
  change every number it draws and stay green. The net is `tests/testthat/test-plot-*.R` plus a
  before/after `ggplot_build()` diff across the affected calls at the merge base.

## After running

Report per harness: which shape it was, PASS/FAIL, and for a FAIL the section and the largest
deviation. State which harnesses you judged not applicable and why, and name the gaps above.
**Do not report a change as verified on the strength of harnesses that do not reach it.**
