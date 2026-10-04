# Session handoff

State, not policy. Policy lives in `CLAUDE.md` and changes rarely; this file changes every
session. Maintained by `/handoff`.

## Now

**`cdca670a` is pushed to `dsem-v5-integration` (2026-10-03), at `5.23.0.9001`.** It makes
`retrospective(peels = )` take a vector of peel depths and fixes three reporting defects in the
same function. Written up in full in the commit message; the summary is under Done & verified.

Everything below about PR #111 predates it — **#111 was CLOSED unmerged on 2026-09-09** and this
branch is parked, now 316 commits behind `dev` and 123 ahead. Treat the #111 sections as the
record of what the branch contains, not as an open review.

**The version keeps colliding, and this session proved the guidance below is load-bearing.**
Both lines have independently used 5.19.0, 5.21.0 and 5.22.0, and each time the branch had to
move up after the fact and renumber a NEWS section already written. `dev`'s number always wins,
because `dev` is what releases. Before bumping this branch again, read `dev`'s DESCRIPTION — or
better, leave the bump until the PR is about to merge, since anything chosen earlier is a guess
about what `dev` will not use.

This session bumped to `5.24.0` anyway, against exactly that advice, and an adversarial review
caught it: `dev` has published **5.24.0, 5.24.1 and 5.24.2** (`git show dev:NEWS.md`), and
`dev`'s DESCRIPTION is at **5.48.0**. A bare `5.24.0` here would have made an install from this
branch report a real release holding entirely different code — reversing `95153bbc`, the commit
whose whole purpose was to stop that, and the GOA arrowtooth DSEM script installs from here.
**The NEWS file's collide-and-annotate convention does NOT extend to `DESCRIPTION`**; `95153bbc`
deliberately separates them. The section is now headed `# Rceattle (development version)` and
`DESCRIPTION` carries a `.9001` development suffix. Do not take a bare number on this branch.

The last block of work answered a report that **`self_test(process = TRUE)` did not simulate
process error** — on a DSEM or without one. It did not, for five separate reasons; all five are
fixed, plus two features and two crash/consistency fixes that came out of the same conversation.

## Done & verified in `cdca670a` (2026-10-03)

**Golden regression FAIL 0 | WARN 4 | SKIP 0 | PASS 20.** `SKIP 0` is the point — the test
carries both `skip_on_cran()` and `skip_on_covr()`, so a green line means nothing without it.
The four warnings are the GOA fixture's own `Time_varying_q` / `Time_varying_sel`
soft-deprecation notices. Reference objectives unmoved: `ss = 10241.0304272585`,
`ms = 10267.2478324443`, `goa_ss = 12868.0052289274`, `goa_ms = 12932.7931701136`.

**Affected suite FAIL 0 | WARN 0 | SKIP 0 | PASS 368** across `test-functions-retrospective`,
`test-retrospective-selectivity-pinning`, `test-diagnostics-print-contract`,
`test-hindcast-skill`, `test-dsem-retrospective`, `test-dsem-naive-equivalence`,
`test-vignette-api`.

**`verify-refit-like.R`: all nine shared sections bit-identical** — the refit-path check golden
does not cover (rule 3). The "before" digest came from `95153bbc` in a `git worktree`, so this
validates the commit rather than only seeding a baseline: `base`, `retro`, `jitter`,
`self_test`, `profile`, `remove_F`, `mse_normal`, `mse_regen` and `mse_ms` (the multispecies one,
which is the only section exercising the pinned predation-suitability window) all `identical()`.
The sole reported difference is `DIFF: retro_sub`, the section `cdca670a` adds.

What shipped:

1. **`retrospective(peels = )` takes a vector of depths.** Scalar `n` still means `1:n`;
   `peels = 2:10` fits those nine; `peels = c(5, 5)` de-duplicates to the 5-year peel alone,
   which is the only way to name one depth because `peels = 5` and `peels = c(5)` are the same
   object in R. `$peel_depths` records the request and `print()` attaches a `NOTE` when the set
   is not `1:n`.
2. **`forecast_rec = "model"` errored on every non-DSEM fit** — it called a `.pin()` that was
   never defined anywhere. `hindcast_skill()` *defaults* to `"model"`, so it was unusable on
   exactly the models built for the comparison it exists to make. Now `pin_block()`, defined once
   in `retrospective()`'s dispatching frame.
3. **Three once-per-call warnings never reached a caller** — see Known flags.
4. **A no-random-effect fit no longer claims a biased retrospective.** `fit_mod()` stores
   `random_vars` as `character(0)`; previously `$<-` with a NULL deleted the element.

**Measured, and it changes how the comparison should be run.** On a BS2017SS DSEM (lag-1 on the
recruitment deviations plus a temperature covariate), 3-year peel: `forecast_rec = "mean"` vs
`"model"` leave the peeled hindcast **bit-identical** (max difference in hindcast recruitment
exactly 0) and move forecast recruitment **46.8%** — `"mean"` flat at 35.2 million, the SEM
decaying 35.2 -> 25.8 -> 18.7 million — but move forecast **SSB only 0.03%**, because the
recruits the two rules disagree about are not mature inside three years. The SSB MASE agreed to
four significant figures for two of three species. **Compare projection rules on
`quantity = "R"`, not the default SSB.** Peel cost is ~1 min each (they warm-start from the
parent fit); the parent DSEM fit itself took 6.1 min.

Both `proj_mean_rec` settings returned the **same objective, 10452.7**, under `estimateMode = 1`
— the projection is not in the likelihood there — so the two-fit comparison buys nothing over
scoring one fit both ways. `model-diagnostics.Rmd` now says so.

## Done & verified on #111

Full suite **198 files, 0 failures, 4 skips** (all pre-existing), run after the `dev` merge that
brought component profiles, on the pushed tree.

## Done & verified on #111

Full suite **196 files, 0 failures, 4 skips** (all pre-existing: three `Skipping` guards and
`test-dsem-recruitment.R:23`, whose reference objectives predate v5.0). Verified after the `dev`
merge, on the pushed tree.

- **Golden reference green** — the four fits reproduce their pinned objectives. Every new draw is
  gated on `dsem_on == 1`, so non-DSEM numerics are untouched by design and by measurement.
- **`tools/verify/verify-dsem-equivalence.R` PASSES — for the first time ever on this branch.**
  See "Known flags"; it was dead.
- **End to end on a `family = "normal"` DSEM**: fit, redraw, write back, rebuild, refit, both
  replicates converged.
- Golden, full suite and the equivalence harness were all re-run **after** the merge, not only
  before it.

## What the five defects were

Worth reading before touching `sim_mod()` or `dsem.hpp`; the detail and the measured numbers are
in `inst/dev/TRAPS.md` under "The SIMULATE contract".

1. **`calculate_dsem()` had its own `do_simulate` branch that discarded the R draw.** It assigns
   the whole latent field (`x_tj = draw_tj`), and it ran *after* `sim_mod()` substituted its
   conditional draw — adding 100 to every latent state changed nothing. It consulted neither the
   map nor `dsem_cond_k`, so it also redrew the covariate columns; under `family = "fixed"` those
   columns **are** the `env_data` series, so each replicate was generated under an environment the
   refit never saw (a scaled covariate with sd 1 moved by up to 2.96). The draw is now taken only
   in R and `dsem_do_sim` is pinned to 0.
2. **The mask stayed zero, so no truth came back.** `rec_dev_drawn_sim` is written inside the
   template's IID block, which is gated off under a DSEM, so `attr(x, "process_sim")` was `NULL`
   on exactly the models whose process *had* been redrawn — indistinguishable from nothing
   happening, and enough to stop `compare_sim()` warning it was measuring bias against the wrong
   baseline. The R draw now carries its own mask and hands back `dsem_x_tj` as well as `rec_dev`.
3. **`init_dev` was gated on `dsem_on == 0`**, so a replicate's initial age structure kept its
   fitted values while the hindcast was re-realized — two histories in one data set.
4. **The AMAK/Ianelli gate did not reach the SEM draw.** `srr_fun = 0` with `srr_pred_fun > 0`
   scores the deviations twice; the template has always refused the standard draw for that, the
   SEM draw did not.
5. **`self_test()` swallowed every `sim_mod()` warning.** Replicates run through `parLapply` and a
   worker's warnings never reach the caller, so at the default `cores` a `process` request that
   redrew nothing was completely silent. The unique set is now re-emitted once, and asking for
   process error and getting none back warns in its own right.

## The two features from the same conversation

- **A DSEM covariate with a measurement family is simulated, state and observation.**
  `family = "fixed"` is data and is untouched; any other family makes the covariate a latent state
  observed with error, so the state is drawn with the field and the observation is drawn from its
  own density (all seven families) and written back into `env_data`. Checked as a moment: over 200
  draws, `sd(y_sim - mu)` was 0.9995 against a fitted measurement sd of 1.
- **`build_bounds(dsem = )` bounds the SEM's standard deviations at 0.** A lag-0 two-headed self
  path is a Cholesky diagonal whose sign is unidentified, so the surface is exactly symmetric
  about 0 — harmless at the MLE, fatal for MCMC. Only diagonals that stand alone in their row
  qualify; a cross-covariance or a lagged `<->` keeps unbounded support and `fit_mod()` names those
  variables. A fit at the negative root is flipped rather than rejected (10401.62 either way).
  `vignette("dsem")` shows the `tmbstan` call.

## For Grant's review

- **`NEWS.md` 5.22.0 folds two release blocks into one.** The branch had filed its DSEM work as
  5.19.0 while separate from `dev`, and `dev` released a different 5.19.0. Both existed after the
  merge. `dev`'s keeps its number; the branch's is merged into 5.21.0 under the four normal
  headings (39 bullets in, 39 out). Read it once as a single release — the wording of the two
  halves has not been harmonized.
- **The PR body is out of date in the worst direction.** It still says `sim_mod(process = TRUE)`
  refuses on a DSEM and that a `self_test()` does not re-realize the process. Both are now the
  opposite, and that is the headline of the change.

## Known flags

- **A `warning()` raised inside a `.parallel_lapply()` worker is DISCARDED**, and the default
  `cores` is `detectCores() - 6`, so the parallel path is what anyone actually runs — the warning
  is visible only at `cores = 1`, which is what the suite passes. This hid three of
  `retrospective()`'s warnings completely; all three are now settled before dispatch
  (`R/9-retro_and_jitter.R:245-300`). Two had also been gated on `i == 1L`, which `peels = 2:10`
  never reaches; the third sat in `.rce_peel_map()`, which has no peel index, so it fired once
  per peel instead of once per call — **so an old `cores = 1` log shows that one N times, not
  once.** Hoist when the condition reads only the input model; collect and re-raise as
  `self_test()` does with `sim_warns` when it is genuinely per-item. Now in `TRAPS.md`.
- **`retrospective()`'s per-peel "could not report hindcast standard errors" is STILL lost**, and
  it is the one that cannot be hoisted (it depends on whether *that peel's* `sdreport()` threw).
  Worse than a missing band: the handler falls back to `newmod$sdrep`, which at that point is the
  **forecast refit's** sdreport — built with the whole hindcast pinned, so every hindcast standard
  error in it is **zero**. A peel whose `sdreport()` failed presents a zero band as real, with
  nothing saying so. `R/9-retro_and_jitter.R:700-703`. Point estimates unaffected. Needs
  collect-and-re-raise; not done.
- **A saved pre-`cdca670a` fit with no random effects still draws the bias warning.**
  `random_vars` is absent on it, which is indistinguishable from a pre-5.10.0 fit that never
  recorded it. Correct for the latter, a false alarm for the former; refitting resolves either.
- **`tools/verify/verify-refit-like.R` gained a `retro_sub` section in `cdca670a`**, so its
  "bit-identical across all sections" verdict is unreachable against any baseline captured before
  that commit — the compare reports `DIFF: retro_sub` and nothing else. Its verdict is a
  whole-object `identical()`, then a name-based detail loop, so nothing mis-pairs. **Re-baseline
  on `cdca670a` or later**; `dev/refit-after.rds` from the run below is one, but `dev/` is
  untracked scratch and does not survive a clone.
- **A `tools/verify/` harness runs in NO CI job, so a broken one is invisible.**
  `verify-dsem-equivalence.R` — the only numerical check on the vendored `dsem.hpp` — could not
  run for the whole time it sat on this branch: `cond_k` was added to `calculate_dsem()` and to
  `dsem_standalone.cpp`, but the harness's `dsem_data()` was never taught to supply it, so every
  invocation died in `getParameterOrder()` before the first comparison. **Adding an argument to
  `calculate_dsem()` means updating three places, not two** — the header, `ceattle.cpp`, and
  `tools/verify/dsem_standalone.cpp` *plus* the R harness feeding it. Now in `TRAPS.md`.
- **A file split defeats git's rename detection and merges cleanly wrong.** `dev` moved
  `self_test()` and `profile.Rceattle()` out of `R/9-retro_and_jitter.R` into new files. Those
  arrived *already staged as adds*, never appeared as conflicts, and a `git commit` of the merge
  took `dev`'s copies — dropping every change this branch had made to both functions, silently,
  leaving a package that loads and functions that exist. Fixed in `6ff40df5` by taking the
  branch's version of each and re-applying `dev`'s two small deltas on top. **After any merge that
  moves a function between files, assert each one is defined exactly once and grep for your own
  change inside the new file.**
- `sim_mod(simulate = FALSE)` returns `env_data` as supplied, not as fitted values — the one data
  type that does not return its hat. Deliberate and documented in `?sim_mod`.
- An `env_data` column that is also a `linkage_spec` covariate is an errors-in-variables case
  under a measurement family: the replicate is generated from the series as supplied while the
  refit sees the redrawn one. Faithful to the real assessment, documented, not a defect.
- **`jnll_comp` columns count fleets on rows 1–8 and species on rows 9–20**, so `rowSums()` pools
  two different axes. `.JNLL_ROW_AXIS` (`R/9-profile.R`) is now a third hand-synced partner to the
  `JnllRow` enum and `R/6-rename_output.R`; all three must move together.
- **`unweighted_jnll_comp` is written for 5 of its 21 rows** — composition, CAAL, stomach and the
  two linkage rows, the ones carrying a `Comp_weights` multiplier. Everything else is
  structurally zero there, not small, so `profile_components(weighted = FALSE)` returns a much
  smaller set of series.
- **The QAR1 block in `ceattle.cpp:4269` is dead code.** `data_check()` refuses
  `Catchability = 6` (`R/1-data_check.R:92`) and `est_index_q` is only ever set from that column;
  the live QAR1 form is a q linkage. Its AR1 density was scoring into the "Catchability prior"
  row and now scores into "Catchability deviates" — a reporting fix that no fit can reach.
  `R/6-process_residuals.R:204` still branches on `est_index_q == 6`, which is why the block was
  fixed rather than removed. **Removing the dead QAR1 path — the C++ block and that R branch
  together — is Grant's call and wants its own change.**
- **A model carrying GOA numbers forward needs a refit.** Result-changing changes on this line not
  labelled breaking: the mode-5 selectivity penalty fix (GOA Pacific cod SSB 2050 −14.1%),
  parameter bounds previously applied to the wrong parameters, `remove_F()` zeroing fitted hindcast
  F when `suit_endyr < endyr`, composition weights warm-starting from `inits`, failed `run_mse()`
  simulations returning only a marker, the `mse_summary()` reshape, the recruitment fixes, and
  `sim_mod()` drawing the index under the fleet's own `Index_distribution`.
- **5.10.0 moved three predation figures' numbers** — `plot_m2_at_age_prop()` (a share, not a
  contribution), `plot_ration()` (× average numbers-at-age) and `plot_b_eaten()` (million mt).
  `plot_selectivity()` also renames `p$data$Age` to `Bin`.
- **MSE projection statistics are not comparable across 5.13.0.** `sim_mod()` draws the index
  under the fleet's own `Index_distribution` now, which shifts the RNG stream. At `nsim = 2` the
  hake summary swung `catch_iav` 0.25 vs 0.74 between branches on identical OM and EM fits.
- The golden reference on `dev` is `ss = 10241.0304272585` (`ms = 10267.2478324443`,
  `goa_ss = 12868.0052289274`, `goa_ms = 12932.7931701136`), pinned in
  `test-golden-regression.R`.

## Blocked

Nothing.

**Cleared 2026-08-27: the bounds check runs.** `verify-safebounds.R` now builds the package
SHLIB as well as the TMB model, asserts both DLLs loaded before any case runs, and reports
"could not run" apart from "violation". All **six** configurations pass locally from a clean tree
(all three `.so` deleted first), which is the CI condition — the sixth is
`test-selectivity-catchability.R`, the file that crashed on Windows, which used to be a second CI
step and was never bounds-checked at all. 206 s end to end.

**`continue-on-error` is off now**, so a bounds failure fails the run. It gates nothing —
`deep-checks` is nightly and dispatch-only, never a PR check, and `main` is not protected — but
the failure is visible, which it was not before.

The history below is kept because it is the reason not to trust a green `deep-checks` run
without reading the job.

<details><summary>What was wrong</summary>

**The bounds check had measured nothing since at least 2026-08-25, and reported green.**
`deep-checks`' `safebounds` job is `continue-on-error: true`, so the workflow's overall
conclusion is `success` while the job fails. It fails at the DLL, not on a violation:

```
confirmed: -DTMB_SAFEBOUNDS in the compile line
Failed to load at least one DLL.
  ERROR: 'ceattle' was not found in the list of loaded DLLs.   x 5 cases
BOUNDS VIOLATIONS or errors in 5 case(s):
```

All five configurations error before running, so **zero cases are bounds-checked**. The runs on
2026-08-25 (5.20.0 dispatch) and the 2026-08-26 nightly failed identically, so it is not new.

Cause: `NAMESPACE` loads two DLLs (`useDynLib(Rceattle, .registration=TRUE); useDynLib(ceattle)`).
`verify-safebounds.R` rebuilds only the TMB lib, then calls `pkgload::load_all(compile = FALSE)`
— deliberately, since a recompile would overwrite its bounds-checked model — so nothing ever
builds `src/Rceattle.so`, and `src/Makevars` has `$(SHLIB): tmblib`. On a fresh CI checkout that
file has never existed. **It passes locally only because a developer's tree already has one**,
which is why "five configurations clean on macOS" in the 5.20.0 notes was not the clearance it
looked like.

Both are fixed. It is the only net for the Windows `0xC0000005`, which is still uncured — see the
5.21.0 release notes' known limitation.

</details>

## Resume here

**Three concrete follow-ons from `cdca670a`, in order of value:**

1. **Run the projection comparison on a real assessment.** Everything in `cdca670a` was verified
   on `BS2017SS` and synthetic fixtures; the question it was built to answer — does a DSEM's
   correlation structure improve the projection — has not been asked of a GOA model. The recipe
   is in `model-diagnostics.Rmd` under "Does a DSEM's correlation structure improve the
   projection?": one fit at `build_srr(proj_mean_rec = FALSE)`, then `retrospective()` twice at
   `forecast_rec = "mean"` and `"model"`, scored with `hindcast_skill(quantity = "R")`.
   `peels = 3:10` is about 30 min at ~1 min per peel.
2. **Fix the per-peel `sdreport` warning** (Known flags). It is the last of the four
   worker-swallowed warnings and the only one that needs collect-and-re-raise rather than
   hoisting. It currently presents **zero standard errors as real**, which is the worst failure
   shape of the four.
3. **Re-baseline `verify-refit-like.R`** on `cdca670a` or later (Known flags).

**Then decide what this branch is for.** It is parked, #111 is closed, and it sits 316 behind
`dev` / 123 ahead. `cdca670a` is self-contained — the `peels` vector, the `.pin()` crash fix and
the `random_vars` fix are all independent of DSEM — so **the `.pin()` fix in particular is worth
checking against any other branch that carries that call**, since it makes `hindcast_skill()`
unusable on a non-DSEM fit wherever it appears.

Earlier #111 notes, now historical: the PR body described only the original DSEM work, not
`bound_sd`, the zero-variance guard, or `$bounds$par_lower`. `run_mse()` and two
`process_residuals()` processes still refuse on a DSEM and are the remaining gap in that suite.

Loose ends inherited from `dev`, none blocking:

1. **`inst/dev/SIBLING-REPOS.md` has an uncommitted edit that predates these sessions** — the
   `../GOA-multispecies-assessment` entry. Coherent and complete; left out of the 5.21.0 and
   5.22.0 commits as unrelated, not as wrong.
2. **`../Rceattle-models/GOA pollock/2025/04-fit-and-diagnostics.R` gained an M-at-age-6
   component profile** (uncommitted, separate repo, not yet run). `minage = 1`, `nages = 10`, and
   `M1_base` runs 1.39 at age 1 to 0.30 from age 6 on, so age 6 is bin 6 and 0.30 is the value a
   0.18-0.42 grid brackets. It fixes the age-6 cell **alone**, leaving ages 7-10 at base --
   narrower than goa_pk's own M profile, which scales the whole vector. **5.22.0's
   `joint = "multiply"` now expresses the goa_pk version**; switching it is a two-line change.
3. **The dead QAR1 path.** `ceattle.cpp:4269` was fixed in place rather than removed because
   `R/6-process_residuals.R:204` still branches on `est_index_q == 6`. Removing both together is
   Grant's call.

After #111: `inst/dev/CLEANUP_BACKLOG.md` has 64 items left (Tier 0 and Tier 2 cleared);
`inst/dev/BACKLOG-PLAN.md` sequences them by who is exposed. Queued and not started:
`inst/dev/CONTRIBUTOR-EXPERIENCE.md` — **item 0 comes first and is a conversation, not code**: ask
the sibling-repo authors where they actually stopped.

## Older paused work

**The `accessibility-and-code-review` refactor is superseded, not outstanding.** The branch is
gone from local and remote, its plan no longer exists, and no commit in any branch mentions it.
Three of its four locked "chosen extras" landed by other routes — the `JnllRow` enum, the repaired
cpp section index, and Doxygen `@file`/`@brief` blocks on the previously bare headers. Its one
concrete leftover, splitting `R/0-build_srr_and_M.R`, is done in 5.14.0. Nothing of that plan is
outstanding. Do not resume from the branch name or the handoff path; both are dead references.
