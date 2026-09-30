# Session handoff

State, not policy. Policy lives in `CLAUDE.md` and changes rarely; this file changes every
session. Maintained by `/handoff`.

## Now

**5.45.0 is released.** Tagged `5.45.0` on `main`'s merge commit `b4506079` and published
2026-09-28; the `release: published` event fired pkgdown (the 5.21.0 silent failure did not
recur), `deep-checks` was dispatched, and **checklist section 4 passed**: installed from the tag
into a temporary library in a clean session, `packageVersion()` reads 5.45.0, `citation()` and
the `fit_mod` example resolve. `dev` carries 5.45.1, the `golden` robustness fix.

**Read this before trusting `deep-checks` `golden` again.** It is reworked at 5.45.1, so the
expected signature has changed: a `goa_ss` delta of 52.9 is no longer the thing to look for,
because the gate no longer re-optimizes. And the red seen at the 5.45.0 release (run
36433121293) was **not** the 52.9 at all -- it was the MVN/Normal configuration pin failing its
GRADIENT assertion, 3.5e-04 on ubuntu against a 1e-4 gate where local macOS gives 2.5e-05, with
its objective reproducing and the other 19 assertions passing. That gate is 1e-3 from 5.45.1,
with the reason recorded in the test. `TMBhelper` is installed on the runner, so a different
optimizer explains neither. See the 5.45.1 NEWS entry and `TRAPS.md`.

**`cod-bridge` is at 5.46.0 and open as PR #178 into `dev`.** It carries seven features from the
SS3 cod bridge: `initMode 6`, the SS3 growth / maturity / length-bin options,
`Selectivity = "DoubleNormalSS3"` (code 15), length-based selectivity on the population bins,
the initial equilibrium catch, and a per-fleet ageing error matrix.

**Its review fixes are on `fix/cod-bridge-blockers`, not yet merged into it.** That branch merges
`dev` (5.45.2) and fixes three blockers, each reproduced before it was fixed:

- **`catch_data` lost its alignment with `catch_hat`.** `clean_data()` kept the `styr - 1` row and
  `fit_mod()` stores the pre-`rearrange_data()` list, so on `GOA2018SS` `catch_data` was 372 rows
  against a `catch_hat` of 370. `plot_catch()`, `residuals(source = "catch")` and
  `sim_mod(simulate = FALSE)` threw; `run_mse()` indexes `catch_hat` by `catch_data` row position,
  so it shifted all 160 projection rows by two years and wrote `NA` into the last two. The rows now
  live in `data_list$equil_catch_data`.
- **Length selectivity was normalized and projected on the DATA bins** while the curve had moved to
  the population bins. On a 20-data / 96-population fixture, projected `sel_at_length` came back
  zero on every bin, and with it projected F, catch and the reference points.
- **The equilibrium catch was read under `initMode` 3 and 4**, whose initial age structure is not
  built at the `Finit * selectivity` the prediction assumes. Restricted to mode 6.

Measured at `14040f5d`, serial, R 4.5.1: **suite 9,976 / 0 / 0** (3 skips), **golden 21 / 0**,
**Pacific hake MSE reproduces all four reference objectives exactly** (2440.0942, 2440.6633,
2447.0049, 2669.3776), `verify-mse-hindcast-invariant` 0.000e+00, `verify-sim-recovery` and
`verify-sim-centering` both clean.

**Two things found on the way that outlive this PR.** `TMB::compile()` tracks no header
dependency, so a `.hpp`-only edit leaves the old object and `load_all()` reports success while
running the previous model -- see `TRAPS.md`, and treat any number measured after a header-only
edit as suspect. And `combine_data()` errors on the bundled datasets (`plyr::rbind.fill` on
`GOA2018SS`'s `maturity` / `sex_ratio`, a different error on `BS2017SS`) before reaching any of
this code; pre-existing, worth its own issue.

**The thing to know before touching the equilibrium catch.** It is a `catch_data` row at
`styr - 1`, and that year is NOT a free marker: `GOA2018SS` carries 23 catch rows before `styr`,
two of them on 1976, and it is the only bundled dataset that does. An earlier version read those
as equilibrium observations, predicted 0 under an `initMode` that holds `Finit` at 0, took
`log(0)`, and returned a non-finite objective on both GOA golden references. **That failure was
recorded for most of a session as a pre-existing `goa_ss` problem. It was not; it was this.** A
negative sentinel cannot be used instead -- `run_mse()` reserves negative `Year` for rows it
splices in as the next assessment's data, and its window filters are on `abs(Year)`, so `-999`
survives as year 999. What ships is the year, read only under `initMode 6`, with `data_check()`
naming the fleets whose rows it reads and saying when a mode will not read them.
`.rce_equil_catch_rows()` and `.rce_equil_catch_candidates()` hold the rule so `clean_data()`,
`rearrange_data()` and `data_check()` cannot drift.

**Two lessons from the merge worth carrying.** The branch's selectivity codes were stale --
it had 13/14 as two integrable forms where `dev` now has a single 13 -- and the merged C++ had
already resolved to dev's design, so the R maps were made to match the template rather than the
other way round. And five registries owed entries for this branch's features
(`not_a_column`, `.QUANT_INFO`, the jnll axis scanner, `R/data.R`'s `@format`, the pinned
fleet_control defaults); every one was caught by a schema test rather than by reading the diff.

**Read the suite from testthat, not from the log.** `grep` for failures missed golden's error
three times in one session, because testthat writes `── 1. Error (...)` and the pattern looked
for `^ERROR`. `as.data.frame(testthat::test_local(reporter = "silent"))` gives the counts
directly. A capped run also aborts on max-failures, so "N failures" from a capped log is a floor.


**`dev` is at 5.45.2 and `main` at 5.45.0.** The next step is one `dev` -> `main` release
covering 5.45.1 and 5.45.2, per `inst/RELEASE-CHECKLIST.md`. Read that file's pkgdown note
before tagging: the `release: published` event has silently failed to fire once already.
**The tag is the DESCRIPTION version, so read it off `DESCRIPTION` at the moment you tag; it
has moved four times during this release (5.41.0 -> 5.42.1 -> 5.43.0 -> 5.45.0) as review and
follow-up work landed on `dev`, and `README.md:43` has had to follow it each time.** #160 (5.42.0),
#169 (5.44.0, 5.45.0) and the
review of #158 (5.42.1) both landed after the release PR was written.

**Checklist state at 5.45.0.** Two of the four measurements have been re-taken at this head,
and the other two are argued rather than re-run:

- **Full suite, re-run 2026-09-28 at `be207905`, `dev` at 5.45.0**
  (`NOT_CRAN=true TESTTHAT_PARALLEL=false`, serial, R 4.5.1 on macOS with every Suggests
  installed): **9,765 assertions / 0 failures / 0 errors**, 3 skips, 244 files. Supersedes
  9,702 (round 3), 9,613 (5.42.1) and 9,506 (2026-09-21). **Re-take it after anything lands on
  `dev`**: this figure has gone stale four times in this release, most recently when #169 merged
  minutes after the 9,702 run. The failure count is
  the load-bearing number; **the skip count is environment-specific** -- 459 `skip_on_cran()`,
  763 `skip_if_not_installed()` and 93 `skip_if()` guards mean a clean machine will skip far
  more, so do not treat 3 as a target.
- **Ecosystem sweep, re-run and widened to the 5.42.0 refusals**: 375 workbooks across the
  four consumer repos, 183 with a `fleet_control` sheet (the count excludes `~$` Excel lock
  files; including them gives 423). Exactly one workbook, one column, two fleets carries a
  negative weight: EBS pollock 2024's `Sel_curve_pen1` on AVO and ATS, `NonParametricPM` with
  no `Sel_shape_mode` column at all, which the directional exemption keeps legal -- and ATS
  additionally follows AVO's `Selectivity_index`, so the template never reads its weight.
  `Sel_curve_pen1`/`2` exist in only 20 of the 183 workbooks and `pen3` in 1. No workbook
  sets `Sel_shape_dir` or `Sel_devmag_sd` -- those columns are absent everywhere, so the
  present-but-NA case never arises. No selectivity prior or apical linkage sits on an `Off`
  fleet or a shared-block follower -- the GOA pollock 2025 prior fleets are each their
  group's lead, and GOA cod bridging builds double-normal linkages with no priors. No
  group anywhere mixes selectivity forms, read raw or canonicalized through `sel_map`.
- **Hake `MSE_yr2024.R`: the 2026-09-23 run still stands**, and re-running it would prove
  nothing new. `MSE_hake_yr24_final.xlsx` is **`nspp = 4`** (Hake, ATF, Sablefish, CSL) with
  only two fleets, both `Selectivity = 5` (`Hake`), on separate `Selectivity_index` values,
  with no negative penalty weight and no `Sel_shape_dir` or `Sel_devmag_sd` column.
  Form 5 appears in **no** slot of `.RCE_SEL_PEN_POSITIVE`, so no
  5.42.0 refusal can fire on it, and with no shared group the 5.42.1 lead rule is a no-op
  there. Nothing in 5.42.0 or 5.42.1 touches predation, suitability, the DM likelihood,
  `sim_mod()` or `run_mse()`'s numerics (rule 15); the only MSE-visible change is that the
  estimation fits now report `NOTE` instead of `OK` under `getsd = FALSE`, a status.
  Do not read the two-row `fleet_control` as a single-species model: only hake has fishery and
  survey data, the other three are diet-only predators, so the script's four-element
  `suitMode` / `suit_styr` / `suit_endyr` vectors and `msmMode = "MSVPA"` are correct. A review
  pass misread this as a stale script; it is not.
- **Reproducible install: still owed at this head.** It was driven at 5.41.0. Run it as part
  of checklist section 4 once the tag is pushed.
- **`urlchecker::url_check()` and `devtools::spell_check()`: run 2026-09-24, both clear.**
  Neither is run by CI, so they are only ever done by hand. `url_check()` reports exactly one
  404, the `adding-a-selectivity-form.html` canary in step 4 below, which is the stale
  published site rather than a bad link and resolves when pkgdown rebuilds. `spell_check()`
  returns several hundred domain terms (`acf`, `ADMB`, `ADREPORT`, `al`, `Ageing`); there is
  no `inst/WORDLIST`, so it is advisory noise, not a gate. Adding a WORDLIST so this becomes
  a real check is a `CLEANUP_BACKLOG.md`-sized job, not a release one.

**PR #159 merged into `dev` on 2026-09-24**, #160 on 2026-09-25 (`dev` head `c01ea717`), and
the review of #158 after it. #159 asked four questions of #158: does the language read as
AI-written, is the API frictionless, are the docs concise, can a developer find and change the
model. Full suite green (235 files, 0 failures), golden unchanged to ~1e-11.

What a reviewer should still go at hardest:

- **The selectivity form collapse.** `NonParametricIID` (13) and `NonParametricRW` (14) became
  one `NonParametricIntegrable` (13), with `Time_varying_sel` picking the structure. Code 14 is
  free. **Golden does not cover it**: the four reference models use forms 0-4 (verified: the
  union of their `Selectivity` columns is {0,1,2,3,4}), so golden passing only shows the new
  `sel_case` dispatch is inert for the other forms. What covers the merge is
  `test-selectivity-nonparametric-integrable.R` and its independent `dnorm` oracle.
- **The collapse left stale text in six shipped schema descriptions and in NEWS 5.40.0.** Fixed
  on `fix/release-doc-corrections`; `test-docs-anchors.R` now fails if any schema description
  names a selectivity code `sel_map` does not accept. A schema `doc` string is written verbatim
  into `meta_data_names.xlsx`, so it is user-facing, not a comment.
- **The language sweep touches 69 files** and is isolated in one commit. It recasts clauses
  rather than transliterating dashes; the risk to look for is a `carry` that meant *propagate*
  being flattened to `hold`. Six such were caught in roxygen; assume more exist.

## Round 3 of the #158 review (branch `fix/pr158-review-round3`, 4 commits)

Folded into 5.43.0, no bump: `dev` was already there and 5.45.1 is reserved. Two adversarial
passes ran over it; between them they found nine and eleven items, of which these mattered.

- **A negative variance printed as a standard error** in `summary()` and `report_tables()`'s
  parameter table, which is the table a SAFE executive summary is built from. `vcov()` is
  `sdreport()`'s covariance with no `pdHess` gate and an indefinite Hessian inverts without
  being positive definite, so `sqrt(abs(variance))` reported a meaningless number unflagged.
  Measured on `diag(2, -3)`: `chol()` fails, the diagonal is `0.5, -0.333`, the old code
  reported `0.707` and `0.577`. Both now return `NA`, with a rounding tolerance of `1e-10` so a
  converged fit on a flat ridge still reports `0` rather than `NA`.
- **A factor switch column was read by its level index.** `revert_switches()` resolves eight of
  the nine columns `convert_switches()` handles; **`Time_varying_q` is the ninth**, so a factor
  there reached the template by level index through `fit_mod()` -- `factor("Off")` became `IID`,
  estimating catchability deviations nobody asked for, with `data_check()` clean. One shared
  rule (`.rce_defactor_fleet_control()`) now runs in `switch_check()` and `convert_switches()`.
  Nothing in the ecosystem supplies a factor (375 workbooks, 15 bundled `.rda`, every
  `data.frame()`/`read.csv()` in the four consumer repos), and R >= 4.1 defaults it off.
- **The blank-`Fleet_type` guard ran after `convert_switches()` had destroyed the evidence**, so
  a mistyped type was refused as a blank cell. It now reads the column as supplied, and a type
  outside the allowed set is refused by fleet and value -- closing `3`, `-1` (which was eligible
  to lead a `Selectivity_index` group whose penalty is gated on `flt_type > 0`, so the group's
  penalty went uncharged) and `2.7`.
- **Two claims the second reviewer made did not survive checking**, and one had already been
  written into NEWS on its word: `revert_switches()` de-factors `Fleet_type` four lines before
  the `"Off"` assignment that was said to corrupt it, so neither the `switch_check()` NA
  corruption nor a fit-moving factor `Fleet_type` through `fit_mod()` is real. Verify a
  reviewer's mechanism before it reaches a release note.
- **Deliberately left open**: `R/1-data_check.R:1255` still reads `Selectivity != "Fixed"` raw,
  so raw `GOA2018SS` fleets 4 and 5 are still false positives of the class 5.43.0 fixed -- and
  they are what satisfies the positive control in `test-switches-fleet-type-integer-off.R`, so
  fixing the line turns that assertion red. The class fix (canonicalise every schema `switch`
  column at `data_check()`'s entry) wants its own PR and a golden run; it is the
  `CLEANUP_BACKLOG.md` row "schema `switch` not enforced at the boundary".
- The two transfer notes no longer ship in the tarball. `TRAPS.md`, `CLEANUP_BACKLOG.md` and the
  other `TODO-*` notes still do, because `man/run_mse.Rd` cites two of them.
- The repository `homepage` setting pointed at `grantdadams.github.io/Rceattle/index.html` (404);
  it is now the live site. That was a GitHub setting, not a file, so no sweep would have caught it.
- **A green local suite did not establish the tests were sound.** The width-alignment assertion
  located its column by searching for `"ages "` in the printed line, and which elements the 90%
  loading cut keeps moves with rounding in `eigen()`: a consecutive set prints `ages 1-14`, a
  scattered one `13 ages in 1-14`, so the token sits three characters further along. It passed
  locally on both lines and failed the `R-CMD-check` `oldrel-1` leg (`FAIL 1 | PASS 3050`). It now
  measures the first non-space past the padded block name. Twice in this review the tests were
  the weak link rather than the code -- the 5.42.1 round found three that passed for the wrong
  reason -- so when a new assertion reads printed output, check what in that output is
  platform-dependent before trusting one machine's green.

**The release sequence, from here:**

1. **Done.** Every review branch is in: `fix/release-doc-corrections`, #162, #164, #165, #166,
   and #161 (merged at `eafcece3` with no further bump, since `dev` already read 5.43.0; its
   NEWS entries were filed under 5.42.0 and have been moved to 5.43.0).
   `fix/pr158-review-round3` carries the third review of #158 -- see "Round 3" below. **5.45.1
   is spoken for by the `golden` robustness fix**, so a further review round folds into 5.43.0.
2. Merge the `dev` -> `main` release PR #158. Its body must say what forces a refit, what
   breaks and what is new, and must cover 5.42.0, 5.42.1 and 5.43.0; do not paste `NEWS.md`. Suite and
   sweep are already re-taken at this head (above); the install is checklist section 4.
3. Tag the MERGE COMMIT on `main` with the DESCRIPTION version, then publish a GitHub Release
   from the tag.
4. **Confirm pkgdown rebuilt with the canary, not by eye**, then
   `gh workflow run deep-checks.yaml --ref main`. The canary is
   `https://afsc-assessments.github.io/Rceattle/articles/adding-a-selectivity-form.html`, which
   **404s today** and must return 200 after the Release is published. It 404s because that
   article landed at 5.37.0 on `dev` and does not exist on `main`, which is what pkgdown
   builds from -- so it is a live test of the exact silent failure the checklist warns about,
   with a known-bad starting state. Still 404 after publishing means the `release: published`
   event did not fire: `gh workflow run pkgdown.yaml --ref main`.
5. Tell the consumer repos to pin the tag rather than track `main`.

**The last installable tag is `5.28.0`, not 5.33.0.** `main` carried 5.29.0, 5.30.0, 5.31.0,
5.32.0, 5.32.1 and 5.33.0 without a tag being pushed for any of them. So a consumer who pins
tags, which is what step 5 asks for, moves **5.28.0 -> 5.45.0**, seventeen minor versions (29 through 45
inclusive), not eight. Say that in the release body, and treat step 3 as the fragile step it has proven to be:
the pkgdown `release: published` miss at 5.21.0 is the same step failing in a different way.

**Three things are red before the release starts, and none is from this batch.** None is
evidence against the tag. But "known" is not "ignore": read each post-merge run against the
signature below, because that is the only thing separating a known red from a new one:

- **`deep-checks` `golden` fails on `main` at 5.33.0.** `goa_ss` lands in the second local
  minimum and `goa_ms` inherits it through its warm start. Last fully green 2026-09-13.
  **Before dispatching `deep-checks` on `main`, expect exactly this signature: a `goa_ss` delta
  of 52.9 with the other three models bit-identical.** Any other pattern is a real regression
  and stops the release. Diagnose from the gradient at the reference `par`, not the objective
  (`TRAPS.md`). The robustness fix is the first job after the release and ships as 5.45.1; until
  it lands this guard cannot gate anything.
- **`deep-checks` `suite` never finishes.** It is `cancelled` in every recent run at 5h00-5h01
  wall clock, a timeout rather than a pass. So the one job that runs the 140 `skip_on_cran()`
  test files proves nothing about `main` right now. Either raise the ceiling or shard it; until
  then do not cite a `deep-checks` run as suite coverage.
- **`R-CMD-check` fails intermittently, and it is mostly macOS, not Windows.** Over the last 100
  runs: 79 completed, 8 failures, of which 7 were macOS-only and 2 touched Windows. The Windows
  access violation is real and reproduces on `main` (the file the framework names carries no
  information, see `TRAPS.md`), but it is ~2 in 79, not 2 in 30. **It shows up as a dead
  testthat worker** (`parallel_event_loop_chunky` -> `handle_error` -> `cli_abort`), and the
  log then dumps whatever that worker had printed. In the 2026-09-25 run on #158 that was
  `test-convergence.R`'s `[FAIL] max_gradient = 4e+12 (largest on 'sel_inf')` and
  `[FAIL] pdHess` -- which are `print()` output from a deliberately non-converged synthetic
  fixture (`make_fake_fit()`; the file runs no TMB fit and passes 46/0/0). **Those lines are
  not a convergence regression.** Check for the dead worker before reading a Windows red as
  one.

macOS was red 2026-09-20 to 09-22 for two unrelated upstream reasons and recovered on its own;
branch `ci/macos-libomp` holds an unmerged remedy if the OpenMP one recurs.

The 2026-09-14 backlog plan is finished. Nine branches across eight versions, listed below in
version order (#150 merged before #149; 5.37.0 took two branches), each reviewed
adversarially before commit and again by a second session before merge:

| Version | PR | What landed |
|---|---|---|
| 5.34.0 | #144 | MSE dynamic-SB0 and fixed-numbers reference-point masking |
| 5.35.0 | #145 | Silent-wrong-number fixes; `estDynamics = 3` retired |
| 5.36.0 | #146 | Config overlay by field; OSA outliers flagged per panel |
| 5.37.0 | #147, #148 | QAR1 path removed; stored-map guard; `CONTRIBUTING.md`, the Doxygen build and `adding-a-selectivity-form.Rmd` |
| 5.38.0 | #149 | Per-sex apical selectivity offset (`log_sel_apical`) |
| 5.39.0 | #150 | Multispecies stock-recruit bounds and a degenerate-curve check |
| 5.40.0 | #151 | Two integrable non-parametric forms (13, 14), collapsed to `NonParametricIntegrable` (13) in #159 |
| 5.41.0 | #152 | `osa_residuals(method = "cdf")` |
| 5.42.0 | #160 | Selectivity-penalty sign refusals; apical/prior `Off` and shared-block gates |
| 5.42.1 | review of #158 | Penalty lead keyed as the template keys it; the code-14 guard widened |
| 5.43.0 | review of #158 | A blank `Fleet_type` is refused, in `switch_check()` and `rearrange_data()` |

## After the release, in order

1. **Make `golden` robust**, so `deep-checks` can gate a release. Warm-start the reference fits
   from the pinned parameters, or take the lower of two starts. It is a harness change and
   cannot move a fitted number. **This is not just the next cleanup: it gates the NOAA
   transfer** (`PLAN-adoption-and-NOAA-transfer.md` section 0, item 5) and it has a release
   vehicle already chosen, 5.45.1 (`TODO-pre-transfer.md` B3). Do it before anything below.
   While doing it, fix the `deep-checks` `suite` timeout too; a guard that cannot finish is
   the same problem in a different job.
2. **Decide on the three inert test guards** (`CLEANUP_BACKLOG.md`): restore or delete. The
   multispecies one is hiding an unexplained disagreement with the old EBS CEATTLE.
3. **Work `SIMPLIFY-LOG.md`.** Seventeen rows, three struck through, so **fourteen open**: six
   change behaviour (two of them needing a deprecation path or a shim), one moves a golden
   reference, four are internal, two additive, one doc. Every row is logged rather than done,
   by standing rule; Grant picks which become PRs.

Done for this cycle, so do not repeat them: the ecosystem sweep of the four consumer repos, and
the hake `MSE_yr2024.R` run. Both are recorded above with their results.

## Open work, by where it is recorded

- `TODO-selectivity.md` — `Hake` ignores `Sel_norm_scope` (and `Sel_norm_bin_upper`), with the
  fix sketch and the double-normalization trap in it; the parked `sel-penalty-form` branch; and
  whether a bias-corrected Laplace or `tmbstan` should be the recommended way to report the
  non-parametric deviation SD.
- `CONTRIBUTOR-EXPERIENCE.md` — items A (two of three recipes), B, E, G and H are open, and so
  is item 0, which is still the one that should reorder the rest.
- `TODO-srr-multispecies.md` — initial ages do not decay with M1 + M2 under predation. Built
  and measured during #150, then reverted: year-1 N and M2 feed back across predation
  iterations, and `BS2017MS`'s default starts reached an infinite objective at `niter = 10`.
- `CLEANUP_BACKLOG.md` — everything found and deliberately not fixed, in tiers. Absorbed
  `TODO-5.34-followups.md`: the PFMC `Ftarget` assignment and the `goa_ss` second minimum
  live in `TRAPS.md`, the unbounded `log_Ftarget` and the dead average-F branch here.
- `TODO-maturity.md` — new. Whether `maturity` should carry a bin column and an age/length flag
  as `comp_data` does, rather than an age sheet with length scalars beside it (recommendation:
  yes, own PR, reuse `Age0_Length1` rather than inventing a fourth spelling). And `GOA2018SS`
  Cod maturity reading 2.0 at ages 1-12 with no range check in `data_check()` — found, not
  diagnosed, and it feeds SSB directly.
- `TODO-projection-module.md`, `TODO-mse-horizon.md` — unchanged by this batch.
- `TRAPS.md` — verified traps with the measured numbers behind them.
- `PLAN-adoption-and-NOAA-transfer.md` — moving the package off a personal account to a NOAA
  org, and the two adoption barriers behind it. **Section 0 holds five decisions reserved for
  Grant** (destination org, license, co-maintainer, whether `Rceattle-models` moves, timing);
  agents do not pick these.
- `TODO-pre-transfer.md` — the execution checklist for that plan, stages A-F with owner tags.
  Stage B is this release. **B3 is the `golden` robustness fix**, to ship as 5.45.1 if it lands
  after the tag.

## Tagged snapshots of the DSEM lines (2026-09-24)

Neither line is a release, and both report a version that must not be mistaken for one. Pin
the tag, never the branch: the branches move, and an assessment refit from a moving branch
does not reproduce.

- **`dsem-v5-2026-08-27`** -> `95153bbc`, annotated. The `dsem-v5-integration` snapshot the
  **GOA arrowtooth 2026 assessment** runs against. 121 commits divergent from `dev` and
  missing 5.24.0 through 5.42.1, so it does **not** carry that range's silent-wrong-number
  fixes. Pinning makes that run reproducible, not current; whether it should run on this line
  at all is a scientific call, not a tooling one.
  The branch previously declared `Version: 5.23.0`, which is a **published release tag**, so
  an install from it reported `packageVersion("Rceattle") == "5.23.0"` and any provenance
  record built from it named the wrong Rceattle. Bumped to `5.23.0.9000` in `95153bbc`.
- **`dev-DSEM-v4.5-archive`** -> `a82c99f5`, lightweight, and **already equal to
  `dev-DSEM`'s head**, so it needs no new tag. Declares `Version: 4.5.0`, and no `4.5.0`
  release tag exists (the line goes 4.4.1 -> 4.6.0), so there is no collision here. 796
  commits behind `dev`; superseded by `dsem-v5-integration`.

**Outside this repo, and still owed:** `GOA-ATF-ESP/R/2026 assessment model-DSEM.R` was
repointed from `@dsem-v5-integration` to `@dsem-v5-2026-08-27` but **the edit is uncommitted
in that repo**. Until it is committed, that assessment still installs from the moving branch.
The `@dev-DSEM` pins in `Rceattle-models` (GOA pollock 2025, EBS pollock 2024 and its
README), the 2025 ATF script and `GOA_circlulation_study` are all commented out, so they bite
only whoever uncomments one; they would each need a `dev-DSEM` tag, which is a different and
older line.

## Parked branches

- `sel-penalty-form` (`Sel_penalty_form`, 5 commits) — parked by decision, not by defect.
- `dsem-v5-integration` — PR #111 closed unmerged 2026-09-09. Head `95153bbc`, tagged
  `dsem-v5-2026-08-27`; see the snapshots section above. Grant's plan is to bring `dev`'s
  updates onto it later.
- `reporting-tables` — **local only, never pushed**, so it is not on the remote to triage. Its
  one stray doc commit, `4716968c`, reached `dev` as `3255fb49` via PR #132 (which merged from
  `docs/minfraction`, so the content was re-applied rather than merged from this branch).

## Resume here

**Finish the release.** Grant is doing it in a session after this one, so this is where to
start rather than `SIMPLIFY-LOG.md`.

1. **Decide #161** -- renumber and include, or hold. It is code, not docs, and its bump is
   stale (see the release sequence, step 1). Nothing else blocks the merge.
2. **Merge PR #158.** As of 2026-09-25 it is `CLEAN` / `MERGEABLE`, and CI is green on all
   five platforms -- Windows passed on both runs, which is worth noting given the intermittent
   access violation. Re-check before merging; the branch has moved since.
3. **Tag the merge commit with whatever `DESCRIPTION` reads then (5.45.0 today)**, bare, no `v`
   prefix, then publish a GitHub Release from
   it. **This is the step that has silently not happened five times** (5.29.0 through 5.33.0
   are all untagged), so do not defer it or hand it on.
4. **Run the canary** (release sequence, step 4). It 404s now and must return 200 after.
5. **Dispatch `deep-checks` on `main`** and read `golden` against the 52.9 signature above
   before concluding anything from it.
6. **Then 5.45.1: make `golden` robust**, which gates the NOAA transfer.

Two loose ends that are not release-blocking: commit the `GOA-ATF-ESP` pin change in that
repo, and the `deep-checks` `suite` 5h timeout, which belongs with the 5.45.1 work.
