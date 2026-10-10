# Session handoff

State, not policy. Policy lives in `CLAUDE.md` and changes rarely; this file changes every
session. Maintained by `/handoff`.

## Now

**`dev` is at 5.57.0** (`07b9f0a3`). Six PRs landed after 5.54.3; **two** carried a version.

| PR | version | what it did |
|---|---|---|
| #217 | 5.56.0 | DoubleNormal (type 8) started at peak 0 and a logit floor of 10, giving a curve flat at ~1 and no gradient. Derived starts: peak at the mid-range BIN ORDINAL `(nages + 1) / 2`, floor at logit 0. Objective 3085.98 -> **914.10** on the GOApollock fishery |
| #218 | -- | three doctrines at the top of `CLAUDE.md`; `AGENTS.md` and a new `.github/copilot-instructions.md`; `CONTRIBUTING.md` split |
| #219 | 5.57.0 | a sel/q linkage on a shared block estimated one FREE coefficient per member, so mirrored fleets fitted different curves. Now one coefficient, plus four refusals |
| #221 | -- | `source-guards.R` said 17 where `EXPECTED` held 18 |
| #223 | -- | `TRAPS.md`: a map built from a skeleton model freezes both DoubleNormalSS3 ends |
| #224 | -- | the three-file doctrine contract is now enforced by a test, not declared |

**#222 is open and reviewed, not merged.** Titled "condense the documentation" at +3962/-5579, but
it touches the schema and the template, so the numeric surface was measured rather than read:
**all six `src/TMB/` files are comment-only** (strip `//` and `/* */`, diff the remainder), **20 of
23 `R/` files are comment/roxygen only** (compare deparsed `parse()`), and the schema's **86
columns, their order and every non-`doc` argument are identical**. The whole code delta is three
string literals, one of them a message now pointing at `vignette("model-diagnostics")`. It overlaps
#219 on four files, so it wants a real merge, not a fast-forward.

**The numbers a next session should not re-derive:**

  * **Open Tier 0 rows: 9 under the strict rule, 19 counting every `**Open` variant.** Say which
    you mean: `**Open**` matches 9 in-table rows; `**Open, low**`, `**Open.**`, `**Open, latent**`
    and the rest bring it to 19. Count IN-TABLE rows only -- the section's own prose line is a hit
    and is not a row. (Was 10/20; #219 closed the shared-block row.)
  * **An intercept `init` on the DONOR of a shared block does not reach the optimizer's starting
    value** -- new Tier 0 row, found writing #219 and deliberately not fixed there. The push writes
    the donor's cell, the others keep their build default, and `TMB:::updateMap()` starts the level
    at their mean: `log(0.5) = -0.693147` asked for, **-0.096574** actual. Pinned as current
    behaviour in `test-linkage-intercept-base-param.R`. Fixing it moves any fit that sets one.
  * **A linkage offset is applied PER FLEET** (`inf_offset(param, flt, sex, yr)`), the apical height
    included. Only a PRIOR is re-targeted onto the shared base. So "put it on the lead fleet" is
    wrong for a design column, and the pre-5.57.0 apical refusal said it anyway.
  * **The suite is 266 files / 1651 blocks / 10,567 assertions** at `NOT_CRAN=true` serial, ~2h50m.
    Golden is inside it; it is 6 blocks / 21 assertions and must show **skip=0**.

**The method that found everything today, and the only one worth carrying forward:** every probe
needs a POSITIVE CONTROL -- a case that must come back non-empty -- because every failure this
session was an under-measurement reported green. Flatten prose with `tr '\n' ' '` before counting
(a line-based grep of a wrapping sentence returns a confident 0: this bit four times, once inside
a control experiment that therefore proved nothing), parse rather than grep for code, and scope the
probe to the CLAIM rather than to the file you edited. Two real defects in #219's own fix were found
this way after a green suite; neither would have surfaced otherwise.

## Next, in the order I would take them

The refactor plan's steps 1-6 are done: dead code deleted, the file->role table (67 of 67 `R/`
files named in the developer guide, was 21), the partial-match net, the guard job, the
`CLAUDE.md`/`CONTRIBUTING.md` split plus `AGENTS.md`, and one look per figure family. What remains
in that plan is documentation, cosmetics, or a two-week decomposition -- and `CLEANUP_BACKLOG.md`
says catch advice is wrong. So the order below is the backlog, not the plan.

1. **`// Input SB0 (if running in multi-species mode)` -- `msmMode = 0` with `estDynamics > 0`,
   NARROWED rather than what the row first claimed.** Measured 2026-10-08: the F and
   catch-advice half is **wrong** -- `proj_F = 0.0` is applied AFTER the HCR switch at both
   template sites, and that gate is HCR- and msmMode-independent, so `F_spp` is 0 under every HCR.
   The reported quantities are masked to NA by `rename_output()` (the 11 always-NA ones in 5.34.0,
   `SB0`/`B0` and the depletions in 5.35.0), so no user reading a fit the normal way sees a
   placeholder number. **What is left is real but smaller:** `ceattle.cpp` still accumulates
   `SB0 += NByage0(...)` and `SBF += NByageF(...)` from placeholder recruitment, so a consumer
   reading `fit$obj$report()` or `quantities` before `rename_output()` still gets those numbers.
   Closing it means either overriding `NByage0`/`NByageF` in the template for a fixed-numbers
   species, or deciding the R-side mask is the contract and documenting it where a direct-report
   consumer reads. `/golden-check` is green either way -- no reference model sets
   `estDynamics > 0`.
2. **The two DoubleNormal rows, as one PR.** From the default starts the curve fits flat at 1.000
   at every age (`GOApollock` fishery, static selectivity, phased: objective 3085.98 against
   914.10 from `inits`), and DoubleNormal + `RandomWalkAscending` silently fits 220 parameters
   instead of 316 with no message. A new user following `adding-a-selectivity-form.Rmd` hits
   both, and `test-selectivity-double-normal.R` sets its own starts, which is why the suite does
   not see the first one.
3. **Whether the GOA multispecies assessment's own fit moves is still unmeasured**, and the
   decision not to refit was taken while this file said phasing made it moot. It does not:
   `R/02_fit_models.R:76` sets `phase = FALSE` on the headline multispecies fit. One fit, and the
   highest information per hour in the repo. Maintainer call.
4. **Branch protection for the guards job** (above). Maintainer call.
5. **`CONTRIBUTOR-EXPERIENCE.md` item 0** -- one hour with the `Rceattle-models` and
   `GOA-ATF-ESP` authors. That file still says item 0 "is still the one that should reorder the
   rest", and steps 1-6 of the refactor plan were done without it.

## 5.49.0 review pass on PR #184 (2026-10-02/03)

Four reviewers over the whole `main...dev` delta -- C++/likelihood accuracy, R-pipeline accuracy,
language/documentation, and concise-clean-code -- plus the release-checklist items the PR body
listed as outstanding. What was fixed is in `NEWS.md` under 5.49.0. What follows is the state of
the gates and what was deliberately NOT done.

**Two reviewer findings were checked and REJECTED. Do not re-open them.**

* `inst/RELEASE-CHECKLIST.md`'s `.libPaths(c(lib, .Library))` was reported as unable to run,
  on the reasoning that `remotes` lives in a user library the line drops. **On this machine the
  user's packages are IN `.Library`** (`/Library/Frameworks/R.framework/.../library`), so
  `remotes`, `TMB` and `dplyr` all resolve after that line and the recipe runs. It would break on
  a machine with a personal library, which is worth knowing, but the recipe is correct as
  written here.
* `selectivity.hpp`'s scalar `binwidth` was reported as an edge difference where SS3 uses a
  midpoint difference. **SS3 uses an edge difference too**: `binwidth(z) = len_bins(z + 1) -
  len_bins(z)` and `binwidth2 = binwidth(nlength / 2)`, at `SS_readdata_330.tpl:1612` and `:1640`
  on the pinned v3.30.22.1. Rceattle matches exactly, non-uniform grids included.

**Gates run for this release, with results.**

| Gate | Result |
|---|---|
| Full suite, serial, `NOT_CRAN=true` | see the figure in the release note on #184 |
| `golden-regression` | inside the suite; `skip_on_cran()` + `skip_on_covr()` means a run that skips it measures nothing |
| `devtools::document()` | zero drift in `man/` and `NAMESPACE`; roxygen2 8.1.0 matches `Config/roxygen2/version` |
| `urlchecker::url_check()` | all 22 URLs OK |
| `devtools::spell_check()` | 414 words flagged, **no genuine typos** -- domain terms, British spellings, index notation, and three NEWS entries QUOTING historical typos (`selecitivty`, `specificed`, `dont`) |
| `pkgdown::build_reference_index()` | exit 0. No `man/*.Rd` added or removed and `NAMESPACE` unchanged, so `_pkgdown.yml` correctly needed nothing |
| `devtools::check(--as-cran)` | **1 WARNING, 1 NOTE** -- see below |
| Vignettes, `RCEATTLE_EVAL_VIGNETTES=true` | 13 of 13 execute, 37.5 min. `hcrs-and-mses` is 24.6 of those minutes |
| `Pacific hake/04-mse.R` (rule 15) | all four stages reproduce the 5.33.0 references; `run_mse()` end to end; bit-identical across the review fixes |
| `tools/verify/verify-safebounds.R` | **no bounds violation in 8 configurations**, `-DTMB_SAFEBOUNDS` confirmed on the compile line, including the file that crashed Windows CI and both of 5.46.0's new array extents. This is the right net for the one-bin `lengths_pop(sp, -1)` fix |
| Consumer-workbook sweep vs `dev` | 0 newly broken, 0 newly passing, 0 changed verdicts on 375 files |

**The as-cran WARNING is pre-existing and was deliberately left.** Three or more
`-Wbitwise-instead-of-logical` warnings: `if((forecast(sp) == 0) | (estDynamics(sp) > 0))` and
friends use bitwise `|` / `&` on boolean operands. **The same lines are on `main`** (e.g.
`2b0306f3:src/TMB/ceattle.cpp:1345` is character-for-character the same), and `main` carries 23
instances of the pattern, so this release does not regress it and the checklist's "0 warnings"
bar has been failing quietly for some time. Fixing it is a ~23-site sweep of the model source;
that was judged the wrong thing to do hours before a tag. **It is the first follow-up.** The NOTE
is `.git` in the package directory, an artefact of checking a worktree in place.

**Deferred findings, ranked. All are reviewed and confirmed; none blocks the release.**

1. **The bitwise-boolean as-cran WARNING**, above. `R CMD check` names only three lines
   (`ceattle.cpp` `if((forecast(sp) == 0) | (estDynamics(sp) > 0))` twice, and
   `if((sp == flt_spp(flt)) & (flt_type(flt) == 1))`), but it truncates its "significant
   warnings" list, so that is not the whole set: grepping `src/TMB` for a single `|` or `&`
   between two parenthesised comparisons finds **20 candidates** across `ceattle.cpp` and
   `diet_data.hpp`. **Enumerate them with a direct compile rather than another `R CMD check`** --
   capture the compiler's own stderr, since the check only shows a few. Every operand is a pure
   scalar comparison with no side effect, so `|` -> `||` and `&` -> `&&` is provably
   semantically identical (short-circuiting changes nothing), but it is the model source, so
   `/golden-check` after. Left out of 5.49.0 deliberately: the same lines are on `main`, so it
   is not a regression, and a 20-site rewrite of `ceattle.cpp` hours before a tag is the wrong
   trade.
2. **`Ageing_error_index`'s refusal is a false positive on a partial `age_error`.** The check
   applies to EVERY fleet, defaulting a missing index to the fleet's species, so a model whose
   species 2 has no `age_error` rows at all -- a predator carried for diet and an index -- is
   hard-errored although `rearrange_data()` would size the array and nothing would read the
   empty slice. All 11 bundled datasets have full coverage, which is why the suite cannot see
   it. Fix: restrict the refusal to fleets that actually read a matrix (age comps or CAAL) and
   leave the rest at message level.
3. **The three SS3 growth switches are unvalidated on the inherit path.** `growth_sd_form`,
   `growth_plus_length` and `plus_group_decay` have no schema row, no `switch_check()` allowed
   set and no `data_check()` rule, so an out-of-range code falls through to a different branch
   in silence: `growth.hpp` branches on 1 / 3 / 4 and anything else reads as "no plus-group
   correction", and `length_sd_at_age` reads anything but 2 as SD-in-cm. `plus_group_length = 4`
   inherited without a decay rate gives `exp(0) = 1`, i.e. an UNWEIGHTED mean length over `2A`
   further ages -- `build_growth()` refuses that pairing, the inherit path does not. 5.49.0
   warns that these do not round-trip; it does not validate them.
4. **The equilibrium-catch prediction has no independent numeric test.**
   `test-dynamics-equilibrium-catch.R` asserts finiteness and length, never
   `equil_catch_hat` against a hand-computed `sum_a Finit*s*w*N_eq*(1-exp(-Z))/Z`. The
   `initMode = 6` plus-group is untested too: `test-initmode-fished-selected.R` loops
   `2:(nage - 1)`, stopping one age short of the divisor this release changed. And there is no
   `Finit` recovery harness, so the claim that the observation makes `Finit` estimable is
   unmeasured.
5. **The equilibrium catch always uses the row's `Log_sd`, ignoring `Estimate_catch_sd`.** This
   matches SS3, which reads `catch_se` per row, but on a fleet with `Estimate_catch_sd` 1 or 2
   the same fleet's two catch observations are fitted at two silently different variances, and
   no `equil_catch_sd` is REPORTed so a user cannot see which was used. Needs a decision, then
   either a REPORT plus a vignette sentence, or `est_sigma_fsh` honoured.
6. **A warning raised in `sim_mod()` is lost in a `.parallel_lapply()` worker.** The new
   `.sim_warn_unusable(..., "initial equilibrium catch")` is on `self_test()`'s dispatched path,
   so it is visible only at `cores = 1`, which is what the suite passes. Note `CLAUDE.md`'s
   mitigation sentence cites `sim_warns` "as `self_test()` does" -- **`sim_warns` exists nowhere
   in the tree**, and `R/9-self_test.R` has no `withCallingHandlers`, so either build the
   collection or correct that sentence.
7. **C++ duplication the release paid for.** Four verbatim copies of the population-bin CAAL
   loop (`ceattle.cpp`), three formulations of the population-bin midpoint (`growth.hpp`,
   `selectivity.hpp`, `ceattle.cpp` -- now all three guarded), and `nlengths` / `lengths` left
   dead in both growth functions beside their live replacements OF THE SAME TYPE, where a future
   argument reorder would swap them silently. A `GrowthSpec<Type>` struct would cut
   `estimate_growth()` from 27 parameters to about 15. Bit-identical refactors, so
   `/golden-check` covers them completely -- but not before a release.
8. **The TMB comment review's remaining items.** A late reviewer went through every new C++
   comment. Four statements that were wrong about the code are fixed in 5.49.1; these are not:
   * **Three new helpers have no `@param` block at all** -- `length_sd_at_age()` (10 arguments)
     and `fill_age_length_key()` (11) in `growth.hpp`, against `recruitment.hpp`'s standard of a
     documented unit per argument. `calculate_weight()` is missing `@param` for five new ones,
     `plus_group_decay` (per year) most importantly, and `calculate_selectivity()` for six --
     where `nlengths` / `lengths` are now ambiguous against `nlengths_pop` / `lengths_pop` and
     nothing says which is the DATA grid. These are the arguments a wrong unit would silently
     corrupt.
   * `ceattle.cpp:225` says the length compositions sum each data bin's run of population bins.
     They do not: `pop_bin_lo` / `pop_bin_hi` are read only inside `pred_CAAL`, and a length comp
     takes selectivity through `sel_at_age` and lengths through the transition matrix, already on
     the data bins.
   * Three comments narrate the diff rather than the model ("rather than on a bin average", "so
     such a model is numerically unchanged"), one states bug history in the past tense
     (`ceattle.cpp`'s equilibrium-catch `SIMULATE` block, 9 lines where rule 8 allows 2), and one
     carries an unsourced SS3 claim ("SS3 shrinks at sigmaR/ave_age") where every neighbour cites
     a file and line -- under hard rule 9 that is exactly the shape not to guess.
   * `spawn_output` has no section number of its own but is cited as "section 5.7", which is
     `GROWTH`; and a block numbered 2.3c sits between 2.1 and 2.2.

   **One of that reviewer's findings was checked and REJECTED**: it read `spr.hpp`'s header as
   inverting the female-fraction convention. It does not. `spawn_output` carries `sex_ratio` for a
   one-sex species (`ceattle.cpp` `(nsex == 1) ? sex_ratio : 1.0`, and `mature_females = maturity *
   sex_ratio`) and not for a two-sex one, where `female_split = sex_ratio(sp, 0)` supplies it. The
   header's "the female fraction not already in it" describes the second argument and is accurate,
   if terse, and the `@param` block below it is explicit.

9. **Smaller ones.** A non-numeric `Ageing_error_index` (`"1a"`, or a name, which
   `Ageing_error_name` invites) reads as `NA` and silently becomes the species; `bin_indexed_forms`
   in `R/1-data_check.R` is a fourth hard-coded copy of the same form set; `.rce_equil_catch_candidates()`
   and `.rce_equil_catch_rows()` are behavioural functions parked at the top of the schema registry
   and belong in `R/0-clean_data.R`.

## In flight above `dev` (5.47.0)

One branch. PR #181 is merged, and the `R_init` work is rebased onto that result rather than
sitting beside it. The conflicts it actually raised were `DESCRIPTION`, `NEWS.md`, `README.md`
and this file -- all version numbers or section ordering, no code -- so the prediction that
`NEWS.md` and `DESCRIPTION` were "the real conflict" held for the code and undercounted the
prose. Every code hunk did auto-merge.

**`feat/srr-init-level`, 5.48.0** -- recruitment `R_init`, a fourth recruitment linkage parameter
(code 3). A log-scale multiplier on the initial age structure, read at year 0 only, carrying no
deviate penalty. It exists because `init_dev` is penalised per age, so a stock starting away from
`R0` paid a recruitment-deviate penalty to sit where its data say it sits -- and the optimiser
retires that penalty by LOWERING `R0`. On GOA cod that is 0.37 log units, `R0` 31% below what the
same data support, while terminal SSB moves only -1.35%; `SB0`, `B40%` and depletion all scale
with `R0`, so a penalty on the initial state was displacing the quantities that set the catch
limit. Measured: moving the level out of `init_dev` drops the objective 49.9133 nats, all of it
in `Initial abundance deviates` (54.97530 -> 5.06184), every fitted component unchanged to
1.9e-04.

Two things for anyone merging across it:

- `RCEATTLE_N_REC_PARAMS` goes 3 -> 4 in `linkage.hpp`; it dimensions
  `recruitment_linkage_offset`, so a textual merge taking one side mis-sizes the tensor.
- **The shared-intercept machinery changed.** `build_map` pinned EVERY `(Intercept)`
  `beta_linkage` at `NA`, `build_params` zeroed its starting value, and
  `build_parameter_bounds` loosened its bound and propagated it to the base parameter -- all on
  the premise that a base parameter exists to carry the level. True for all six processes until
  `R_init`, which has none, so `~ 1` estimated nothing and moved nothing while every builder
  reported success. The R sites now share `.is_pinned_intercept()` and its complement
  `.is_level_intercept()` (`R/0-linkage_encode.R`); those predicates must stay in each condition.
  **Don't trust a site count** -- it was published as four, five and seven before anyone counted,
  and it is ELEVEN in R plus one in the C++. TRAPS.md carries the grep and the file:line list.
  The C++ one is slot 19's prior block, which re-targeted an intercept prior onto
  `rec_pars(sp, param)`; for code 3 that read past a `nspp x 3` `PARAMETER_MATRIX` and
  segfaulted. Guarded, so the prior stays on `beta_linkage(i)`.
- Three review rounds found, beyond that: `init`/`bounds` on the level were read raw on a log
  coefficient though documented natural-scale; `N_eq` age 0 was unscaled, making the
  `initMode = 6` equilibrium catch 4.5% wrong; the one-coefficient rule keyed on a stratum the
  C++ discards, so `by = ~ species + age_bin` built five aliased coefficients; a fixed level was
  overwritten by a warm start, with and then without an `init`; and, pre-existing and affecting
  every process, a pinned intercept was never re-zeroed over `inits` (M1 held at 6.69x).

**The combination is now verified, not predicted.** PR #181 changed the q computation and this
work changes `build_params` / `build_map` / `build_parameter_bounds`, which every model goes
through, so neither branch's own suite proved it. Golden and the full suite run on the rebased
tree; see the "Now" section.

**The SS3 parity numbers `NEWS.md` points here for.** Measured on the Aleutian Islands
Pacific cod bridge (SS3 3.30.22.1, model M24_1) with SS3's MLE injected: the largest
difference from SS3 in length-at-age fell from 1.5% to 4.5e-6, in the Jan-1 age-length key
to 4.3e-7, and in fecundity-at-age (mature ages) from a 5-6% Jensen gap to 2.8e-6 -- all
within `Report.sso`'s printed precision. SSB is within 1.5%, the rest being selectivity.
`initMode 6` was found the same way: injecting SS3's MLE and inverting mode 4 left the
initial deviates differing from SS3's `Early_InitAge` by exactly
`const - Finit * cumsum(sel)`, residual 0.00000 at all 13 ages.

**Its first round of review fixes is IN, at `88e7233f`** (what was `fix/cod-bridge-blockers`,
which merged `dev` 5.45.2). Three blockers, each reproduced before it was fixed:

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

**A second review round then added five commits, `23a4d265`..`80e825d3`, pushed 2026-09-29.**
`cod-bridge` on `afsc` was the PR's head branch; `origin` is `grantdadams`, so a push there did
NOT update PR #178 -- which is still how every PR on this repo behaves. Seven blockers, all silent-wrong-number or memory-safety:

- **The population-grid refusal turned itself off.** It read `Bin_first_selected` first and that
  column has no schema default, so `as.integer(NULL) > 1L` is `logical(0)`, `logical(0) | <n>`
  stays length 0, and `which()` excused **every** fleet -- including the non-parametric forms and
  time-varying penalties, which that column has nothing to do with. `Selectivity_dimension`,
  `Selectivity` and `Time_varying_sel` collapsed the same way; the last is now read with `[[`
  because `$` returns `Time_varying_sel_sd` by partial match.
- **`max_bin` was still sized on the data bins** (`ceattle.cpp:239`) while selectivity cases 2, 9
  and 13 build their curve over `nlengths_pop`. On the 20-data / 96-population fixture that is 76
  out-of-bounds writes per fleet, sex and year into `non_par_sel` / `sel_coff_off`, reachable only
  because the refusal above had disabled itself. Inert with one grid, so golden is unchanged.
- **`age_error` coverage was checked per SPECIES but the array is filled per
  `Ageing_error_index`**, so a complete matrix vouched for an incomplete sibling whose missing
  true ages stayed 0 and were renormalized away -- a silently wrong age composition. Now per
  matrix, and a matrix spanning two species is refused.
- **A fleet could borrow another species' ageing error.** The index was checked for existence
  only, and with no index column on `age_error` the index IS the species, so `2` means "species
  2's matrix" to the model and "the second matrix" to the user.
- **`sel_dn6_ends` was read from sex 1 alone** while `build_map()` fixes the end per sex. A
  sex-scoped linkage (`linkage_spec()` takes `sex`) made them disagree: one sex fitted with the
  other's curve shape, its own end parameter free at zero gradient. Refused.
- **A shared `Selectivity_index` was un-shared after the fact.** `adjust_map_shared_params()`
  shares the map at `R/3-build_map.R:76`; the `-999` pass at `:84-85` then re-fixes a follower left
  at the default. Reordering would NOT fix it -- `sel_dn6_ends` reads the starting VALUES, not the
  map -- so the block must agree, and a disagreement is refused.
- **`.SEL_DN6_PARAMS` existed to state the form/parameter rule and was referenced only by a
  test.** `peak` resolves to `sel_inf`, `dn_peak` to `sel_dn6`, and the header and
  `parameter_dictionary()` both call form 15's P1 "peak", so the wrong name was the natural one
  and wrote a slot the curve ignores. Now wired, both directions, **scoped to `sel_dn6`**: the same
  gap on the older forms is real (cases 2/9/13 contain zero references to `log_sel_slp` or
  `sel_inf`) but those configurations are in released scripts, so it needs a sweep -- Open 3 in
  `TODO-selectivity.md`.

Measured at `c8ecff79`, serial, R 4.5.1, after `touch src/TMB/ceattle.cpp`: **suite 10,237 / 0
failures / 0 errors** (4 skips, 251 files), **golden 25 / 0 / 0**. Golden is gated only on
`skip_on_cran()`, so `NOT_CRAN=true` runs it in-suite. The hake MSE and the `verify-*` harnesses
were NOT re-run this round -- nothing touched predation, suitability, the DM likelihood,
`sim_mod()` or `run_mse()`.

**Three flags on that round.** (1) `document()` regenerates `man/dot-comp_aggregate.Rd` and
`man/dot-osa_jointsex.Rd`, and `git diff DESCRIPTION` is empty, so this is NOT the roxygen-version
trap: PR #179 updated the roxygen in `R/7-plot_comp.R` / `R/7-plot_osa.R` without regenerating
`man/`. Left uncommitted deliberately, so the next `document()` will surface it again. (2) Golden
does not run in the PR workflows (`deep-checks` only, `NOT_CRAN=false` at step level), so a green
PR run is not golden clearance. (3) `GOA2018SS` Cod maturity reads 2.0 at ages 1-12 -- see
`TODO-maturity.md`, Open 2; pre-existing, feeds SSB, and golden pins it rather than catching it.

**What the behaviour change costs a live assessment, measured 2026-09-30.** Nobody had put a
number on it. Refitting the 2024 GOA Pacific cod approximation (`Rceattle-models/GOA cod/
2024_pcod.R`, 224 params, length selectivity + estimated vB growth, `maturity` forced to 1, no
`pop_lengths`) on `dev` 5.45.3 against `cod-bridge` 5.46.0: objective 6415.8172 -> 6416.2875
(+0.470 nats), terminal SSB 1,004,925 -> 935,106 mt (**-6.9%**), SSB series mean **-7.1%**,
largest -12.9%, year-1 SSB -9.3%, convergence OK -> WARN (max gradient 0.0013 on R0, <= 0.00085
SE -- benign). With maturity at 1 and no population grid, the maturity-at-length and pop-grid
paths are the identity, so the movement is the selected-body-weight change.

**Do NOT read the SAFE comparison as validation either way.** That script is an *approximation*:
its SSB sits 63% below the accepted 2024 SAFE spawning biomass on 5.45.3 and 66% below here,
correlation 0.68, and it RISES across a series over which SAFE declines. A 3-point move against
a 63% gap is not evidence. The validation reference is the SS3 bridge in `GOA cod/Bridging/` and
`AI cod - Dev/Bridging/`, not this script. Recorded because a first pass at this comparison
printed "moves AWAY FROM SAFE", which is a meaningless verdict on a baseline that far off.

**A third review round, this one about legibility rather than numbers** (uncommitted at the
time of writing). Nothing in it can move a fit; the C++ change removes a parameter no body
read.

- **`pop_to_data_bin` was threaded through four C++ signatures and dereferenced in none.**
  `estimate_growth()`, `estimate_growth_within_yr()`, `calculate_weight()` and
  `calculate_selectivity()` all took it; `ceattle.cpp` §2.3c builds `pop_bin_lo`/`pop_bin_hi`
  from the data array itself. Removed, with the two call-site arguments.
- **Eleven linkage names for six `DoubleNormalSS3` parameters.** `top_logit`/`dn_top` and
  four more pairs, so the refusal message offered eleven options for six slots and the
  vignette and `NEWS.md` documented different sets. Cut to the SS3 manual's six (`dn_peak`
  keeps its prefix because `peak` is `DoubleNormal`'s). Unreleased, so nothing is owed a
  deprecation.
- **Two blocks lifted out of `data_check()`**, which was one 2,330-line function: the
  117-line equilibrium-catch check is now `.check_equil_catch()` and the population-grid
  refusal `.check_pop_grid_bins()`. Both use `.rce_has_data()`, since `has_data()` is a
  closure inside `data_check()` and does not reach a helper.
- **`initMode 6` on a species with several fisheries now warns.** Only the *catch row* was
  refused; the mode itself was silent, and there `Finit` is applied at the mean fishery
  selectivity and is not apical. The initial state sets the SSB scale.
- The population-grid refusal read `v >= 0` on the bin columns, counting a literal `0` as a
  bin. `> 0` now, with a test. No bundled dataset carries one (checked all 11), and
  `N_sel_bins` cannot take it -- a pre-existing check refuses anything outside `1:nbins`.
- `.rce_pop_length_bins()` read `growth_model[sp]` where `data_check()` uses
  `rep_len(..., nspp)`; a scalar from `build_growth()` made species 2 read `NA`, and
  `isTRUE(NA == 0)` is `FALSE`, so it took the population grid. Reachable only on a direct
  `rearrange_data()` call, since `fit_mod()` extends the vector first.
- Docs: `quantity_dictionary()` gave `sel_at_length` and `growth_matrix` on `nlengths` while
  `rename_output()` labels them `PopBin`; `parameter_dictionary()` did not mention that
  `-999` on a `sel_dn6` end is the switch; `plot_selectivity()` plotted population-bin
  ordinals under an axis reading "Length bin" (now "Population length bin" where the grids
  differ, unchanged otherwise, which is what the two existing label tests assert).
- `SPEC-equilibrium-catch.md` said "Status: proposed, not implemented" in the PR that
  implements it, and its §3.1 still described reading the row under every `Finit` mode.
  Both corrected; the user-facing half lives in the vignette, so the note could still go.

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


**At the time of that release `dev` was at 5.45.2 and `main` at 5.45.0.** The next step was one `dev` -> `main` release
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

**`dev` is at 5.54.3 (`433f02eb`).** The repo holds **264** `test-*.R` files. Check `gh pr list`
before assuming nothing is in flight -- several PRs were open when this was written.

**Tier 0 item 1 is narrowed, not open as first written** -- read its row before picking it up.
The working order is `CLEANUP_BACKLOG.md`'s own, and its first rule matters most here:

1. **Reproduce first, in a test that fails -- and NOT at `estimateMode = 3`.**
   `tests/testthat/helpers-fixed-natage.R` builds a fixed-numbers fixture from a fit's `N_at_age`,
   but `fixed_natage_build()` hardcodes `estimateMode = 3`, which leaves the reorganized
   `forecast` at `0,0` and zeroes `proj_F` for EVERY species through the first clause of the same
   gate. A measurement taken there shows F = 0 for the estimated species too and proves nothing
   about the HCR. Use `estimateMode = 0` so `forecast` is `1,1`, and read `forecast` off
   `obj$env$data`, never `fit$data_list`, which never carries it.
2. **Expect `/golden-check` to be green either way** -- none of the four references sets
   `estDynamics > 0`, so it measures nothing here. `tools/verify/verify-golden-cold-start.R` is
   the harness to run for a starting-value or bounds change; neither covers this path.
3. **Decide what the right number IS before changing one.** The row says equilibrium `SB0`/`SBF`
   are built on placeholder recruitment while the projection correctly uses `NByageFixed`. For a
   species whose numbers-at-age are input, what `SB0` *means* is a modelling question, not a
   code question -- `TRAPS.md` already records that a reference point CEATTLE never estimated is
   a number rather than a gap (`Ftarget`/`Flimit` = 1, `MSSB0` = 999 mt), and that under
   `HCR = 0 & msmMode > 0` the depletions deliberately divide by last-projection-year biomass.
   So ask rather than infer: refusing the configuration, reporting `NA`, and computing a
   fixed-numbers equilibrium are three different answers with different consequences for advice.

**Two items are maintainer calls and have been carried for several sessions** (items 3 and 4
above): the GOA multispecies refit, whose "phasing makes it moot" justification has been
retracted, and branch protection to make the guards job required.

**What NOT to pick up.** The refactor plan's remaining steps are deliberately deprioritised: the
save-path unification (invisible to a user, no numeric exposure -- the net in
`test-plot-save-paths.R` was worth building because it found two defects, but the refactor is
not worth a cycle), and the `data_check()` decomposition, which needs
`tools/verify/verify-data-check-conditions.R` run twice and diffed to prove its condition capture
is deterministic before anyone commits two weeks to it.
