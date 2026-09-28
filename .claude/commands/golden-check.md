---
description: Golden-reference equivalence — fit BS2017SS/MS and diff vs a saved snapshot
---

Verify that a numeric-touching change is behavior-preserving: fit the four reference models
(Bering Sea + Gulf of Alaska, each single- and multi-species) and compare against a stored
snapshot.

**Prefer the committed `tests/testthat/test-golden-regression.R`.** It pins the same four
objectives, additionally asserts each is a *converged* optimum rather than merely a
reproducible number, and carries `skip_on_cran()` / `skip_on_covr()`. Note it also runs a
second, larger "golden reference by configuration" block, so it is heavier than these four fits. Run it with
`NOT_CRAN=true`. Use the recipe below when that test fails and you need to see *which*
quantity moved -- it diffs full fit objects (`par`, `jnll_comp`, SSB, R), not just the
objective.

Scratch snapshot path for that diff: `dev/golden-ref.rds` (gitignored; create `dev/` if
missing).

**On a new machine, capture the baseline BEFORE the change you want to test.** The constants
are platform-independent (verified across Windows and macOS/arm64 to all 16 digits), but a
failure you cannot attribute is worse than no check. A `git worktree` on the merge base is
the reliable way to do this.

**Recipe note (important):** every multispecies fit is warm-started from its single-species
MLEs (`inits = <ss>$estimated_params`). The predation likelihood is non-convex, so the start
point selects the local optimum — a fresh-init MS fit lands ~37 units away and is NOT
comparable. Keep these exact recipes; the pinned reference was generated from them.

Wrap everything in one `export PATH=/usr/bin:$PATH && NOT_CRAN=true Rscript -e '...'`:

1. Recompile + load: `suppressMessages(pkgload::load_all(".", quiet = TRUE))`.
   **Every fit uses `newtonsteps = 3`.** Without it `nlminb` stops on an
   objective-*relative* tolerance, i.e. wherever that tolerance happens to bite
   rather than at a stationary point, leaving the GOA fits ~1e-3 in gradient. The
   reference then moves under changes that cannot alter the model -- adding a
   constant to the objective once shifted `goa_ss` by 52.9 units. Polishing pins a
   true optimum (gradients ~1e-11), so the comparison means something. It does not
   remove that 52.9: `goa_ss` has a second local minimum that far above the reference,
   and a one-ULP gradient change can select it (`inst/dev/TRAPS.md`, "Coverage gaps in
   the golden check").
2. Fit the Bering Sea single-species reference (mirrors `examples/Simulation_testing.R`):
   ```r
   ss <- Rceattle::fit_mod(data_list = BS2017SS, file = NULL, inits = NULL,
       estimateMode = 0, random_rec = FALSE, msmMode = 0,
       fit_control = fit_control(phase = TRUE, verbose = 0, newtonsteps = 3))
   ```
3. Fit the Bering Sea multispecies reference from the SS MLEs:
   ```r
   ms <- Rceattle::fit_mod(data_list = BS2017MS, inits = ss$estimated_params,
       file = NULL, estimateMode = 0, niter = 5, random_rec = FALSE,
       msmMode = 1, suitMode = 0, fit_control = fit_control(verbose = 0, newtonsteps = 3))
   ```
4. Fit the Gulf of Alaska single- and multi-species references. The GOA multispecies
   fit uses **fixed M** (default `M1`), warm-started from the single-species MLEs,
   `niter = 3`, `phase = TRUE`. (The published Adams et al. 2022 recipe estimates M via
   `M1_model = c(1, 2, 1)`. That fit **does converge** on the bundled data — finite
   objective ~13440.28, max gradient ~7e-4 — but its Hessian is non-invertible because M1
   is confounded with fleet-8 (GOA_pollock_fishery) selectivity (an inherent M-vs-selectivity
   ridge, *not* predation-M2 ill-conditioning and *not* a code bug). With `getsd = TRUE`
   the failed `sdreport` makes `TMBhelper::fit_tmb` return an `opt` **without** `$objective`,
   which reads back as the historical "`NA` objective"; `getsd = FALSE` returns the finite
   13440.28 cleanly. A truly PD-Hessian/identified M-estimating fit is unreachable without
   pinning M (e.g. an M prior). The reference therefore keeps **fixed M** — stable, PD-Hessian,
   reproducible. Resolved 2026-07-24; kept fixed-M deliberately, not a bug.):
   ```r
   goa_ss <- Rceattle::fit_mod(data_list = GOA2018SS, file = NULL, inits = NULL,
       estimateMode = 0, random_rec = FALSE, msmMode = 0,
       fit_control = fit_control(phase = TRUE, verbose = 0, newtonsteps = 3))
   goa_ms <- Rceattle::fit_mod(data_list = GOA2018SS, inits = goa_ss$estimated_params,
       file = NULL, estimateMode = 0, niter = 3, random_rec = FALSE,
       msmMode = 1, suitMode = 0, fit_control = fit_control(phase = TRUE, verbose = 0, newtonsteps = 3))
   ```
5. Capture per fit (`cap`): `opt$par`, `opt$objective`, `quantities$jnll_comp`, and the
   SSB / R series. Assemble `list(ss = cap(ss), ms = cap(ms), goa_ss = cap(goa_ss),
   goa_ms = cap(goa_ms))`.
6. If `dev/golden-ref.rds` does **not** exist: save these as the reference and report
   "pinned new golden reference" — only do this on a known-good commit.
   If it **does** exist: load it and report the MAX absolute deviation per field per model.
   **PASS** if every field is within `1e-6` -- the same tolerance the committed
   `test-golden-regression.R` uses, so the two gates cannot disagree; otherwise **FAIL** and show the largest-deviating
   fields so I can investigate. NOTE: `goa_ms` (fixed-M GOA multispecies) has a flat
   likelihood ridge — its `par`/`ssb` are ridge-sensitive across *different* code, though
   deterministic on same-code re-runs; judge it primarily on `obj` (and `jnll`). The
   reference values (this branch): ss=10241.0304272585, ms=10267.2478324443,
   goa_ss=12867.9902664788, goa_ms=12932.7902167145 (5.33.0; the GOA values changed when the
   q prior on GOA2018SS fleet 2 became mean-centred under `bias_adjust_proc`).

Report a clear PASS/FAIL summary, per model. Don't change source to force a pass.
