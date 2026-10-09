# CLAUDE.md — Rceattle

Rceattle fits **CEATTLE**, a single- and multi-species, climate-linked, age-structured stock
assessment model. Its output sets **US federal catch limits under Magnuson-Stevens**. A
silently-wrong number here becomes a wrong quota, and it will not announce itself as a crash.

The likelihood is a **TMB / C++** model (`src/TMB/`); everything around it — data prep, fitting,
projection, MSE, diagnostics, plotting — is R.

## The three that govern every change

These are what goes wrong most often here, so they come before any other instinct about style or
structure, including a correct one. They do not displace the hard rules below, and never rules 2
and 9. Accuracy beats speed, every time. The rules and sections below cite them by number, as
doctrine 1, 2 and 3.

**1. Don't change an API or a workflow without asking.** Rule 1 forbids deleting one; this
forbids changing one. An argument name or default, a switch code or its meaning, a schema column,
the names or shape of anything in `fit$`, a message an assessment script reads, how a figure
looks: if a user's existing script would behave differently, stop and ask — even when the change
is plainly an improvement, and even when you are already in the file for another reason. File the
idea as one row in `inst/dev/SIMPLIFY-LOG.md` and finish the task you were given. Rule 1 makes a
*new* name permanent too, so agree the name of a new exported function, argument or column before
it ships; a column then goes in with `/new-column` (rule 4).

What doctrine 1 is **not**. A refactor behind an unchanged API is free (rule 1), and so is a
defect fix: refusing a configuration the model never actually fitted (rule 5), correcting a
number the model computed wrongly, or freeing a parameter it fused, mapped out or started badly.
Make those — measure them under rule 2, and name in `NEWS.md` whose fits move, with
`/ecosystem-sweep` for which sibling scripts stop. Adding a message is additive; changing or
removing one is not. And a silent wrong number is never a `SIMPLIFY-LOG.md` row: that file holds
design wishes, while a known defect belongs in `inst/dev/CLEANUP_BACKLOG.md` or an issue.

Where a trigger and a carve-out both fit — a bad starting value is both a default and a badly
started parameter; a message at the wrong severity is both a message change and a defect — fix
the behaviour and leave the interface alone: warn rather than move the default, as PR #168 did
for the shared-block start, and raise a severity without rewriting the text.

**2. Write the simplest thing that works, and nothing for later.** No helper, wrapper, class,
config layer, option flag or *new* registry this change does not need now, and no second grammar
for what the linkage formula already expresses (rule 12) — but registering in the registries that
already exist (the `JnllRow` partners, `.JNLL_ROW_AXIS`, `.index_rows_natural_scale()`) is not
optional. A helper earns its place at two callers, or when it names a concept the reader needs —
`.is_timeout()` and `.par_description()` have no caller in `R/` at all, and only their own tests
keep them alive. **A guard that cannot fire is not safety:** name the input that reaches it, or
leave it out — PR #181 deleted five such guards, four refusals and a defensive early return,
because the check above had already thrown for every input that could arrive. Name that input and
the guard has earned its place: the failure mode here is a plausible wrong number, so a guard
that *can* fire is cheap. Adding to `data_check()` (2,254 lines) or `fit_mod()` (1,668) means a
named helper beside the existing ones, not another inline block. If the structure of a change
needs a paragraph to justify it, propose it before writing it.

**3. Write for a fisheries scientist, briefly.** A comment is one or two lines and a `@param` is
one sentence; anything longer belongs in `@details`, a vignette, or `inst/dev/`. Comment and
roxygen lines are already 39% of `R/` (16,141 of 41,059), so match a file's *idiom*, never its
longest block — 16 roxygen blocks run past 80 lines and one reaches 395. Give the assessment
reason, the units and the convention, then stop; what stays however long is under "Comments", and
a literature citation is the specification — never strip one. **You may shorten what is already
there:** tightening a comment or roxygen block in a file your change is already editing belongs
in that change, which is what `SIMPLIFY-LOG.md` asks for — unless that change can move a fit,
where rule 2's reviewer must see only the lines that move it. But shortening can make a comment
**false** — check that every function and argument it names still exists (`fit_control()` has no
`osa` argument and no `...`, so a comment naming one describes a call that errors), and keep any
clause the code below still depends on.

Doctrine 1 is not licence to ask about everything: settle ordinary implementation choices
yourself — a local name, which of two equivalent forms, where a helper goes. Wide freedom over
how this package is written; very little over what it means.

## Read first

1. `inst/dev/SESSION_HANDOFF.md` — what is in flight, what is verified, where to resume.
2. `R/0-column_schema.R` — the source of truth for every workbook column and switch.
3. `vignettes/articles/developer-guide.Rmd` — the pipeline, the switch system, the schema, the
   linkage grammar.
4. `inst/dev/TRAPS.md` — verified traps with the measured numbers behind them.
5. `inst/RELEASE-CHECKLIST.md` — the release and tag process.
6. `CONTRIBUTING.md` — the same rules written for a human contributor, plus
   `vignettes/articles/adding-a-selectivity-form.Rmd`, one extension end to end.

---

## Hard rules

1. **Preserve the public API. Deprecate, never delete** — and don't change one without asking
   (doctrine 1). Rceattle ships and has users (`cran-comments.md`). Exported `build_*()`
   arguments carry deprecation paths. Internal refactors are free as long as golden-reference
   equivalence holds.
2. **A numeric change needs `/golden-check`.** Any edit that can move a fit must keep the four
   reference models within tolerance and the suite green.
3. **`/golden-check` does not cover the refit paths, the simulation draws, or any figure.**
   Use `/verify` to pick the right `tools/verify/*.R` harness. For a plotting change the net is
   `test-plot-*.R` plus a before/after `ggplot_build()` diff. See `inst/dev/TRAPS.md`.
4. **The column schema is the source of truth** for switch values, defaults, and column order.
   Add a column with `/new-column`; consume it by its canonical name. Don't hardcode a default,
   an allowed-value set, or a column order anywhere else.
5. **Behaviour, API, or doc change ⇒ `NEWS.md` + `DESCRIPTION` `Version:` + the affected
   vignette, in the same commit** — plus `_pkgdown.yml` when a documented topic appears or
   disappears. `/doc-sync` checks this. Exempt from it: the agent instructions (`CLAUDE.md`,
   `AGENTS.md`, `CONTRIBUTING.md`, `.claude/`, `.github/`), repo tooling (`tools/`), and
   developer notes (`inst/dev/`). The first five are in `.Rbuildignore`; `tools/ci/` and all but
   three files of `inst/dev/` *do* ship, so don't put anything there a user must not read.
   Correcting a figure in an already-released `NEWS.md` entry is not a new change and owes no
   bump (precedent: 389cd4a6). **"Breaking" means no back-compat path**: a
   deprecation that keeps old fits working, or refusing a configuration that never fitted the
   model it described, is a minor bump — see `inst/RELEASE-CHECKLIST.md` for the examples.
6. **Never hand-edit `man/*.Rd` or `NAMESPACE`.** Run `/document`. Check `git diff DESCRIPTION`
   **first** — if the roxygen version key moved, the `man/` churn is the version, not your change.
7. **TMB source is inert until `pkgload::load_all(".")`.** Then test with
   `TESTTHAT_PARALLEL=false`; the parallel workers cannot load a freshly rebuilt DLL.
8. **Comments state current behaviour, units, and the assessment reason — never bug history.**
   See "Comments" below.
9. **Don't invent a switch code, a default, or a unit. Ask.** If it is not stated in
   `R/0-column_schema.R` or a vignette it is not established, and a plausible placeholder that
   survives into a fit is the worst failure mode this repo has.
10. **Fleet invariants.** `fleet_control$Fleet_code` must equal the row number — arrays are
    dimensioned by `nrow()` but indexed by `Fleet_code`, so read columns by row index `i`, never
    by `flt`. Fleets sharing a `Selectivity_index` / `Catchability_index` share ONE parameter block.
    Selectivity bin columns are read on the fleet's own `Selectivity_dimension`, but not on one
    convention: `Bin_first_selected` is a 1-based bin ordinal, while `Sel_norm_bin`,
    `Sel_norm_bin_upper`, `Sel_pen_first_bin`, `Sel_pen_last_bin` and `Sel_cap_bin` are absolute
    AGES on an age-based fleet. Check `rearrange_data()`'s offset before reading one.
11. **`nages` is a count of age bins, not the oldest age.** Ages run
    `minage .. minage + nages - 1`; age `a` sits at index `a - minage + 1`. `minage = 1` hides
    every confusion, and that is all 11 bundled datasets and the live assessments — so write
    `seq_len(nages[sp]) - 1 + minage[sp]` and mean it. A plotter taking an `age`/`minage`
    argument resolves it with `.rce_age_index()` / `.rce_age_plus_index()`, never by indexing
    the array directly.
12. **`linkage.hpp` and `R/0-linkage_encode.R` are in lockstep.** Their process and param codes
    must match — change one, change both. This is the seam DSEM extends. In a linkage formula the
    fixed part goes straight to `model.matrix()`, so a time block is `~ cut(Year, ...)` — there
    is **no bespoke `block()` helper**, and adding one would be a second grammar.
13. **Commits: plain messages, no `Co-Authored-By` trailer.** Imperative subject, ≤72 chars; the
    body says *why*, and gives the numbers that changed. Keep it to that — a few short
    paragraphs. The full account of what was tried and ruled out goes in the PR body or
    `inst/dev/`, where someone looking for it will find it.
14. **The repos listed in `inst/dev/SIBLING-REPOS.md` consume this API.** Sweep them after a
    breaking change (`/ecosystem-sweep`), and refit a real assessment when a sweep is not enough.
15. **The Pacific hake MSEs are the MSE and predation check** — the only end-to-end `run_mse()`,
    and the only routine exercise of estimated suitability and of DM comps with a prior on their
    own weight, none of which `/golden-check` touches. Run
    `../Rceattle-models/Pacific hake/MSE_yr2024.R` after changing predation, suitability, the DM
    likelihood, `sim_mod()` or `run_mse()`; reference objectives are in `SIBLING-REPOS.md`.
16. **Review each stage adversarially, before you commit it** — don't wait to be asked. After
    any self-contained unit (a merge resolution, a refactor, a wired-up feature, a C++ change)
    send a
    SET of read-only reviewers with different attack angles — a claim-by-claim fact check, the
    mechanics, numeric safety, a coherence audit of what you wrote — not one generalist. Each
    angle finds what the others miss, and reviewing before the commit folds the fixes into it.
    `.claude/agents/density-reviewer.md` is the one for a likelihood, AD-taping or `jnll_comp`
    change. Ask each reviewer to say what it verified as *correct*, not only what it found: the
    failure mode here is plausible output, which a green suite does not catch.

---

## Dev workflow

This is a TMB package, so the C++ must be compiled before the R can run.

```r
pkgload::load_all(".", quiet = TRUE)   # recompiles the TMB DLL + loads R/ (after any .cpp/.hpp edit)
devtools::document(quiet = TRUE)       # regenerate man/*.Rd + NAMESPACE after roxygen changes
NOT_CRAN=true TESTTHAT_PARALLEL=false Rscript -e 'devtools::test()'   # full suite, serial (rule 7)
rcmdcheck::rcmdcheck()                 # what CI runs (slow; usually backgrounded)
```

- **Toolchain:** prefix R compile/check commands with `export PATH=/usr/bin:$PATH` — system
  toolchain first, so a Homebrew clang/gfortran does not shadow the TMB build.
- `load_all()` recompiles via `src/TMB/compile.R`; add `compile = FALSE` for R-only changes.
  Compiled artifacts are gitignored — never commit them. `ceattle.o` is the big one and it
  grows with the template (89 MB at 5.45.3, 94 MB at 5.46.0), so treat any figure here as a
  snapshot; `*.so` is a few MB.
- **To run one test file**, make the env's parent the package namespace so internal helpers
  resolve: `e <- new.env(parent = asNamespace("Rceattle"))`, then source the shared helpers into
  it. A plain `new.env()` fails with `could not find function "data_check"`.
- **A test that runs a real `fit_mod()` optimization needs `testthat::skip_on_cran()`** so plain
  `R CMD check` stays fast. Leave fast unit tests unguarded.
- **CI:** `.github/workflows/R-CMD-check.yaml` (multi-OS) + `pkgdown.yaml` + `test-coverage.yaml`
  + `source-guards.yaml` (the registry/doc guards, per PR) + `vignettes.yaml` (weekly,
  non-blocking) + `deep-checks.yaml` (nightly: golden, safebounds, the full `NOT_CRAN` suite).
  No lint config, no coverage gate. **`pkgdown.yaml` triggers on `main` only**, so a PR to `dev`
  gets no pkgdown CI — run `/pkgdown-check` yourself. **Adding or removing a test that reads
  `R/*.R` or `src/TMB/*` means updating `EXPECTED` in `tools/ci/source-guards.R`**, which pins
  the guard set by name so the diff says which one moved.
- **Slash commands:** `/recompile`, `/test [file]`, `/document`, `/check`, `/golden-check`,
  `/verify`, `/new-column`, `/doc-sync`, `/pkgdown-check`, `/ecosystem-sweep`, `/handoff`.

## Layout

- **`R/`** — numbered by pipeline order: `0-*` build/prep helpers, `1-*` data checks,
  `2..5-*` params/map/bounds/rearrange, `6-*` fit + rename output, `7-*` plotting, `8-*` sim,
  `9-*` retro/jitter, `10-*` MSE, `11-*` model averaging. **The numeric prefixes are meaningful
  — don't renumber or rename wholesale.**
- **`src/TMB/`** — `ceattle.cpp` is the main model (numbered section index); process logic lives
  in headers (`recruitment.hpp`, `selectivity.hpp`, `predation.hpp`, `growth.hpp`, `linkage.hpp`,
  `spr.hpp`, `comp_osa.hpp`, `comp_sim.hpp`, `helper_functions.hpp`, `bioenergetics.hpp`,
  `diet_data.hpp`).
  `jnll_comp` rows are addressed by the **`JnllRow` enum** — refer to a row by its constant,
  never a bare integer. The enum has **two hand-synced partners**: display names in
  `R/6-rename_output.R`, and `.JNLL_ROW_AXIS` in `R/9-profile.R`, which records whether a row's
  columns count fleets or species. Adding or reordering a component means updating all three.
  `test-schema-jnll-rows.R` reads the template and asserts they agree.
- **`tests/testthat/`** — **flat**: every test is a top-level `test-<area>-<topic>.R`. Shared
  `helpers-*.R` / `fixtures/` sit alongside. Fast fixtures: `make_test_data()` (single-species)
  or `make_msm_test_data()` (multispecies, incl. diet) with `estimateMode = 3` build a
  non-optimized object. `tests/comparison/` holds WHAM cross-checks (not part of `test_check`).
- **`vignettes/`** do not execute their code by default — several chunks fit real models, so
  running them is far too slow for `R CMD check`. Set `RCEATTLE_EVAL_VIGNETTES=true` to execute
  them, which is what the weekly `vignettes.yaml` job does. On a PR the guard is
  `test-vignette-api.R`, which parses every chunk and checks each Rceattle call names an
  exported function with arguments it has; that catches renames, not return-shape drift.
  `data/` has the bundled example datasets.
- **`inst/dev/`** — committed developer notes (handoff, traps, sibling repos, backlog).
  The ADMB porting notes are a section of `TRAPS.md`.
  The untracked `dev/` is scratch and does not survive a clone.

## Plotting

- **The shared argument vocabulary lives in `R/7-plot_helpers.R`** and is documented once, in
  `?"rceattle-plot-args"` (`@inheritParams rceattle-plot-args`). `line_col`, `lwd`, `lty`,
  `alpha`, `species`/`spnames`, `minyr`/`maxyr`, `incl_proj`, `incl_mean`, `add_ci`,
  `model_names` each go through one resolver — `.as_colour()`, `.rce_line_params()`,
  `.rce_check_alpha()`, `.resolve_species()`, `.rce_year_filter()`, `.rce_proj_divider()`,
  `.rce_mean_line()`. Adding or converting a plotter means calling those, not writing a second
  reading of the same argument; that divergence is what `plot_f()` and `plot_selectivity()` were.
- **`line_col` and `lty` supply values for whatever the figure separates** — predators in
  `plot_b_eaten_prop()`, sex in `plot_ration()`, the year fan in `plot_selectivity()` — not
  always the model. Say which in the function's own `@details`.
- **Base-graphics arguments from before the ggplot migration are accepted and ignored**
  (`right_adj`, `top_adj`, `mod_cex`, `legend.pos`, `single.plots`, `theta`, `ymax`, `cex`), so
  the assessment scripts keep running. Keep them; document them as ignored, with the ggplot
  equivalent.

## Code style

Two external guides settle what this file doesn't — the
[tidyverse style guide](https://style.tidyverse.org) for R, the
[Google C++ style guide](https://google.github.io/styleguide/cppguide.html) for C++ — as
tie-breakers only: the surrounding file wins, and so does TMB/Eigen idiom. Both already describe
most of what is here (snake_case, two-space R indents, ~80 columns: 37,911 of 41,059 R lines), so
what matters is where they don't:

- **C++ free functions stay snake_case** — `calculate_ration()`, `first_difference()` — not
  Google's `CamelCase`. The headers are one template-on-`Type` function per process; Google
  governs spacing, braces, `const`-correctness and where a comment sits, never a rename.
- **R internal helpers keep the `.` prefix** (`.rce_year_filter()`, `.as_colour()`). That prefix
  is how a reader knows it is unexported, and `@noRd` goes with it.
- **`[[` over `$` for a new read.** `$` partial-matches on a list *and* on a data.frame, silently
  and with no warning (checked, R 4.5.1): where `Time_varying_sel` is absent from a
  `fleet_control`, `fleet_control$Time_varying_sel` hands back `Time_varying_sel_sd`. Ten pairs
  among the schema's columns have that shape, and three of the shorter names
  (`Time_varying_sel`, `Time_varying_q`, `Sel_norm_bin`) have exactly one longer sibling — the
  case that resolves silently, since an ambiguous prefix returns NULL. So `$` is safe only while
  its column is guaranteed present; `.pull_int()` in `5-rearrange_data.R` has it right,
  `fc[[col]]`. Write `[[` in new code, and always when the name is computed; don't sweep the
  ~8,500 existing `$`.
  **What keeps the ~50 existing `$` reads of those three names safe is which accessor comes
  FIRST**: `switch_check()` and `data_check()` reach the column through rlang's `.data$` pronoun,
  which *errors* on a missing column instead of partial-matching, so the pipeline dies loudly
  before any bare `$` is reached — and where it does not, `switch_check()` has already supplied
  the column. Rewriting a `.data$` read as a bare `$` trades that loud failure for the silent one,
  which is why the idiom rule below is not cosmetic. `test-schema-partial-match.R` pins it.
- **Match the file's idiom; never translate between them.** dplyr and the pipe are load-bearing
  in the `7-*` plotters, the `0-*` data prep and `5-rearrange_data.R` (105 lines of it); most of
  `1-*` to `6-*` is base R. Neither is more correct here, and rewriting one as the other inside
  the fit pipeline is how a fit moves unannounced: `.pull_int0()` is pinned to the exact result
  type of the `pull() %>% as.integer() - 1` it replaced, because the C++ template reads that
  type. A verb that reorders rows, drops a dimname or returns a tibble does not error.

**Scope discipline.** One concern per commit, and never a formatting-only hunk inside a numeric
change — whoever reviews a fit change should see only the lines that can move it. Don't refactor
next to the fix. Tests and comments are part of the change, not a follow-up (rules 2, 5, 8). No
new dependency without asking; `DESCRIPTION` already imports 20.

## What review keeps asking for

Mined from the pull requests that did this work, so the next change arrives with it already done.
Doctrine 2 covers the unreachable guard and rule 8 the bug history; these are the rest.

- **Ask the predicate's actual question.** In PR #181 a `log()` and a `pmax()` testing whether a
  number is 1 became `abs(q - 1) < 1e-8`, and a group key of `param|fleet`, where `q` has one
  param, became the fleet. (Both were later deleted with the warning they served.)
- **Split a long function by moving a block out, and say in the body that it is a move**
  (8843c506: 117 lines of `data_check()` into `.check_equil_catch()` beside the other trailing
  helpers — same checks, same order, same messages). Watch for a closure the helper cannot reach.
- **Delete dead code rather than guarding it** (PR #192: `src/TMB/Dev`, the Kinzey & Punt
  scaffolding).
- **Derive from the fields that exist before adding a column** — `Fleet_type` already answers
  most per-fleet questions, and a new `fleet_control` column is a schema change (doctrine 1).
- **Prove a defect by reverting the hunk and quoting the number it produces.** Reasoning about
  what the old code "would have done" falsified 3 of 8 claims on one sibling-repo review.
- **Count what executed, not that it passed.** Golden, the `NOT_CRAN` suite, the vignette chunks
  and pkgdown have each reported green while measuring less than they claimed.

## Comments: write for a fisheries scientist, not a programmer

The reader knows fish and knows management. They may not know R or C++ idiom, and they were not
in the room when the decision was made. So:

- **Explain the assessment reason, not the code.** Why this bin edge, this year window, this
  constant.
- **Always give units and the convention.** "metric tons", "log scale", "female SSB".
- **State current behaviour, never bug history.** A comment is read by someone deciding how the
  code behaves *now*. In PR #196's sweep, two of the comments that narrated a past defect turned
  out to be false about current behaviour, not merely dated — and three of the rewrites came out
  false in the other direction, which is why doctrine 3 asks you to check the claim first.
- **One or two lines** — doctrine 3.
- **Exceptions that stay:** a comment explaining why an old input path still executes is
  *behaviour*, not history — phrase it as behaviour. Literature citations (AMAK, Ianelli, ADMB,
  Punt, Holsman, Francis, Methot, Kinzey & Punt) are the **specification**; never strip them.
  In a regression test, provenance is what stops the next person deleting the test — put it in
  one header block above `test_that()`, never repeated inline.

```r
# Bad -- narrates the code
comp_weights[flt] <- 0        # set weight to 0

# Good -- states the assessment reason
# A Dirichlet-multinomial estimates its own weight inside the likelihood, so
# tuning it externally would compete with that estimate. 0 on the log scale
# is a weight of 1, i.e. no external multiplier.
comp_weights[flt] <- 0
```

**Match the surrounding form, not the surrounding volume** (doctrine 3). The codebase favours
explanatory section headers in R and Doxygen on the C++, so write those. Canonical
references: `src/TMB/spr.hpp` for a fully Doxygen-documented header — one of four with an
`@file` block (`comp_osa.hpp`, `comp_sim.hpp`
and `helper_functions.hpp` are the others), where `recruitment.hpp` has `@brief`/`@param`/
`@return` per function but no file block. And any `R/*.R` + its `tests/testthat/*` pair.

**Roxygen:** markdown, regenerated with `/document`. A `@param` is one sentence: what the
argument means, its allowed values, its default. Anything longer belongs in `@details` or a
vignette. Give internal helpers `@noRd` (not just `@keywords internal`) so they generate no
`.Rd`. **Never insert a helper between a function's roxygen block and its definition** —
contiguous `#'` lines are ONE block and bind to whichever object follows, so the helper silently
steals the `@export` and the `@importFrom` tags, and the original loses them. Tests won't catch
it (NAMESPACE isn't regenerated); only the next `document()` will. **Put helpers above the block
or after the function.**

## Domain vocabulary (use these exact terms in plots, docs and messages)

Match this vocabulary in axis labels, documentation, and console messages; don't substitute lay
phrasing.

- **Reference points:** Amendment-56 SPR proxies — F40% = max FABC, F35% = FOFL, B40% = BMSY
  proxy (Tier 3); Tier 1 uses estimated FMSY/BMSY. Don't write "MSY" where an SPR proxy is meant.
- **SSB** = female spawning-stock biomass. **"Recruitment deviations"** (log-scale), not
  "recruitment error".
- **Selectivity:** name the form — logistic / double-normal / gamma / nonparametric /
  semi-parametric. Don't call every dome shape "double-normal".
- **Composition:** age comps, length comps, conditional age-at-length (CAAL). An ageing-error
  matrix applies only where age/CAAL data are fit; length-only stocks have none.
- **Data weighting:** Francis (2011), McAllister–Ianelli, or Dirichlet-multinomial.
- **Diagnostics:** Mohn's rho (retrospective), OSA residuals, likelihood profiles.

## Reference implementations

What to consult when documenting a switch, shaping a workflow, or naming a process argument:

- [**WHAM**](https://github.com/timjmiller/wham) — `basic_info`/`input` argument documentation;
  how a large option surface is presented to an assessment author.
- [**SAM**](https://github.com/fishfollower/SAM) — the `conf` table: a compact, complete
  configuration reference.
- [**Stock Synthesis**](https://github.com/nmfs-ost/ss3-source-code) — control-file reference
  style; exhaustive per-switch documentation with allowed values.
- [**dsem**](https://github.com/James-Thorson-NOAA/dsem) — the DSEM-linked models: formula/path
  grammar, and how linkage structure is specified and reported. DSEM lives on the
  `dsem-v5-integration` branch, not here.
- [**FIMS**](https://github.com/NOAA-FIMS/FIMS) — *its* `.github/copilot-instructions.md`, where
  the two style guides and the scope discipline above come from; `tests/` as the pattern an agent
  copies rather than reinventing.

## Known traps

One line each; the fuller text and the measured numbers are in `inst/dev/TRAPS.md` (the last
section holds every entry below in full).

- **`newtonsteps > 0` can return a parameter OUTSIDE its bounds**: nlminb respects them, the Newton refinement after it does not, in either `.fit_tmb()` path. `convergence` reports `parameters_outside_bounds` (FAIL); default is 0 but golden runs 3.
- **A shared parameter block starts at the GEOMETRIC MEAN of its members' starting values**, not the lead's; inject per BLOCK, never per fleet. Cost a GOA cod survey 18% of its index.
- **`Index_distribution` has a second registry**: a new family must also be classified in `.index_rows_natural_scale()`, or it gets the log-scale residual.
- **`jnll_comp` columns count fleets on rows 1–8 AND row 22, species on 9–20, and neither on row 21** — the fleet axis is NOT contiguous (row 21 is the model-wide linkage REs, row 22 the initial equilibrium catch, by fleet); `.JNLL_ROW_AXIS` is the registry, so `rowSums()` mixes axes and an axis guess that stops at row 20 is wrong.
- **A reference point CEATTLE never estimated is a number, not a gap**: `Ftarget`/`Flimit` = 1, `MSSB0` = 999 mt, per-recruit quantities 0 under `msmMode > 0`.
- **Under `HCR = 0 & msmMode > 0` the depletions divide by last-projection-year biomass**, not `SB0`; don't blank them with a placeholder `SB0`.
- **A fit reports 99 quantities**: enumerate `names(fit$quantities)`, not a `REPORT(` grep; `quantity_dictionary()` is the registry.
- **`retrospective(getsd = TRUE)` can drop peels `getsd = FALSE` keeps**, so Mohn's rho can differ.
- **`unweighted_jnll_comp` is written for 5 of its 22 rows**; the rest are structurally zero.
- **`fit_mod(d, config = cfg)` replaces `d$model_config`**: build `cfg` with `run_config(d, ...)` or every linkage is dropped.
- **`bias_adjust_proc` centres the lognormal priors and the recruitment deviations together** (5.33.0).
- **A `data_list` element without `write_data()`/`read_data()` support round-trips to nothing.**
- **Under a Dirichlet-multinomial `Comp_weights` is a log**: 1 is a starting weight of e.
- **A Pearson residual divides by the effective sample size the likelihood used**; `.rce_comp_pearson()` resolves it.
- **`estimateMode`: prefer the strings.** Mode 4's objective is a placeholder; mode 3's is real and usable before fitting.
- **`fit$obj` (and `fit$sdrep` unless `ConstantF`) is the projection's under any HCR but `NoFishing`**; `fit$identified` and `fit$.conv_hindcast` are the hindcast's.
- **`fit$data_list` is the pre-`rearrange_data()` list**; recompute rearranged fields from `fleet_control`.
- **`Bin_first_selected` is a 1-based bin ordinal; `Sel_norm_bin` is an absolute age** (rules 10, 11).
- **`init_dev`'s ages start at `minage + 1`**; `.PAR_AXIS_OFFSET` is the registry.
- **`condition_number` reads the correlation matrix since 5.26.0**; `covariance_condition_number` is the old value.
- **`getsd = FALSE` leaves `sdrep` NULL** (no `vcov()`, NA bands); the refit diagnostics read the bias-adjust flags and `projection_uncertainty` off `data_list`.
- **`run_mse()` pins the OM's stock-recruit and suitability windows to the pristine `om$`** (`verify-mse-hindcast-invariant.R`).
- **Every error is drawn in a `SIMULATE{}` block beside its density**; a new likelihood family owes a draw (`verify-sim-*.R`).
- **An MSE draw is per observation row**, so changing the OM horizon or row count changes every later draw (`TODO-mse-horizon.md`).
- **The guards are not themselves guarded**: golden runs only in `deep-checks`; keep `NOT_CRAN=false` a step-level `env:`.
- **An access violation is memory corruption**: build `RCEATTLE_SAFEBOUNDS=true` and run `verify-safebounds.R`.
- **A slow fit is the model**: `BS2017SS` takes ~500–700 `nlminb` iterations.
- **A fixed-numbers species (`estDynamics > 0`) reports its input recruits as `R`** but `NA` R0/steepness/SPR0, and in single-species mode `NA` SB0/B0/depletion; `estDynamics = 2` fits as 1 under `msmMode = 0`.
- **An identity-link recruitment linkage turns on `rec_floor_on`**, changing the AD tape; its floors miss projection, SB0 and dynamic-B0 recruitment (`TODO-srr-multispecies.md` item 14).
- Scratch outputs (`Rplots.pdf`, `*_osa.png`, `*.RDS` under `tests/comparison/`) are gitignored.
