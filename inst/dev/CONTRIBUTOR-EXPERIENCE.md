# Contributor experience — plan

Make the codebase navigable, editable and extensible by a **fisheries scientist or ecologist**,
not a software engineer. State and plan, not policy; policy lives in `CLAUDE.md`.

**Shipped in 5.37.0 (PRs #147, #148):** C, D, the A-guard, and the first of A's three recipes.
Open: A's other two recipes, B, E, G, H. Item F was deleted — its premise was wrong (the column
it proposed to fix already holds a canonical string; see "Where item F went" below). Item 0 has
still not been done, and it is still the one that should reorder the rest.

**Success measure for the whole effort, so it can fail:** someone outside the maintainer group
lands a working selectivity form, or a new likelihood family, without a maintainer walking them
through it. One observable event. If that has not happened a year after these items land, the
items were the wrong items.

**Every item is documentation or tooling except G, which changes the API.** Rule 5's
doc-sync obligations apply as written in `CLAUDE.md`; they are not restated per item, only where
an item is *exempt* or needs a version bump.

---

## 0. First: find out where people actually get stuck

**Do this before any of A–H, and let the answer reorder them.**

Everything below rests on the claim that a fisheries scientist cannot navigate this codebase.
The evidence for that claim is a comparison against FIMS's documentation — which is evidence
about FIMS, not about where Rceattle's users stop. Nobody has been asked.

Ask the people who have actually tried to extend this: the `../Rceattle-models` and
`../GOA-ATF-ESP` authors (`inst/dev/SIBLING-REPOS.md`), and anyone who has asked how to add a
form. One question: *where did you stop, and what did you have to ask a person for?*

If the answer is "the workbook, not the C++", then B and the schema docs matter and A does not.
If it is "I could not get it to compile", H jumps to the front. Three conversations invalidate
or confirm eight items.

**Size:** an hour of conversation. **Acceptance:** A–H reordered, in writing, with the reason.

---

## Where the rest came from

A review of [NOAA-FIMS/FIMS](https://github.com/NOAA-FIMS/FIMS) on 2026-08-22 (v0.10.0.9000,
HEAD `f0c5937`). On model content Rceattle is well ahead — FIMS has no reference points, no
ageing error, no composition weighting, no working OSA residuals, one selectivity pair, one
growth form, single sex. **None of that is worth copying.**

What FIMS does better is the contributor path:

| FIMS has | Rceattle has | Verified by |
|---|---|---|
| `vignettes/adding-new-module.Rmd`, 739 lines: ordered file list, validation checklist, troubleshooting | a 465-line `developer-guide.Rmd` plus one per-task recipe, `adding-a-selectivity-form.Rmd` (5.37.0); two more recipes open | `wc -l` both |
| `inst/include/common/glossary.md` — every symbol, meaning, units | nothing | `find` |
| Doxygen with rendered LaTeX, published on the pkgdown site | Doxygen comments in all 11 headers, built by `Doxyfile` and published from the pkgdown workflow (5.37.0) | `ls src/TMB/*.hpp` = 11, all with `@brief`; `.github/workflows/pkgdown.yaml` |
| `CONTRIBUTING.md`, 101 lines | `CONTRIBUTING.md`, 124 lines (5.37.0); `CLAUDE.md` keeps the agent-facing policy | `ls` |
| `use_gtest_template()` / `use_testthat_template()` scaffolding | none | `grep -rn "usethis::" R/` = none |
| `.devcontainer/` + a setup vignette + `setup_fims.sh` | the `export PATH=/usr/bin:$PATH` toolchain trap, documented only in `CLAUDE.md` | `ls` |
| `tidy()` / `glance()` / `augment()` via {generics} | none | `grep` R/ + NAMESPACE |

The concrete gap, as a task: **adding a selectivity form.** In FIMS it is one new file under
`functors/` plus registration, every step written down. In Rceattle it is a new `case` in the
`switch (sel_type)` at `src/TMB/selectivity.hpp:361` plus edits across the R pipeline — and that
list is written down nowhere a human will find it.

### The file list is a hypothesis, not a finding

Assembled by reading the dispatch and the pipeline order, **not** by tracing one form end to
end. Marked always / conditional on that basis. **The first job of item A is to confirm or
correct this table**, and the prompt says so.

| File | When | Why |
|---|---|---|
| `R/0-switches.R` | always | `sel_map` (`:73`–`:87`) — the name the user types. `"DoubleNormal" = 8` is at `:82`; the 5.40.0 forms close the list at `:85`–`:86` |
| `R/0-column_schema.R` | always | the `Selectivity` description (`:118`) is the user-facing switch documentation |
| `src/TMB/selectivity.hpp` | always | the `case` in `switch (sel_type)` at `:361` |
| `R/1-data_check.R` | usually | per-form checks: `N_sel_bins` at `:867`–`:889`; form x `Time_varying_sel` for the non-parametric forms at `:960`–`:984`, Hake at `:986`–`:989`, LogisticPM at `:1002`–`:1005` |
| `R/2-build_params.R` | if the form has new parameters | |
| `R/3-build_map.R` | if the form has new parameters | the densest file in the list: `grep -ci sel` returns 302 lines |
| `R/4-build_parameter_bounds.R` | if the form has new parameters | |
| `R/5-rearrange_data.R` | if the form needs a column passed to TMB | |
| `src/TMB/ceattle.cpp` + the `JnllRow` enum | if the form carries its own penalty | |
| `R/6-rename_output.R` | with the above | display names are hand-synced to `JnllRow`; the selectivity rows are at `:188`–`:189` |
| `tests/testthat/` | always | `test-schema-cpp-dispatch.R` binds map values to `case` labels and will go red until both exist |

**A trap the recipe must state:** `Selectivity` is compared against **canonical strings**
everywhere in the R pipeline. Integers or strings on input become strings at `switch_check()`,
which calls `revert_switches()` before it returns (`R/0-switches.R:1003`, converting at `:1152`;
nothing after that call re-converts the column). They become integers only inside
`rearrange_data()`, on the copy handed to TMB
(`convert_switches()`, `R/0-switches.R:1471`). So a form added to `sel_map` but not to the
per-form string checks passes validation and fails later. Some sites defend both spellings, for
example `R/0-switches.R:944`.

---

## The work, A–H

### Recommended order, and the case against it

**0 → B → C → A → D → H → E → G.**

A has the highest ceiling, which is why I first put it first — but B is cheaper, C is nearly
free, and **A will spend its whole length referring to names that B defines**, so writing A
first means writing it twice. D and H are small and independent. E only pays off once A exists
to be enforced.

The case for A first: if item 0 says people stop at "how do I add a thing", B and C are
navigation aids for a wall they never reach. **Overturn this with item 0's answer, not with my
ordering.**

Sizes are rough and assume familiarity with the pipeline.

---

### A. Task recipes — highest ceiling · 3–5 days for the first, 1–2 each after

Three articles under `vignettes/articles/`:

- `adding-a-selectivity-form.Rmd`
- `adding-a-likelihood-family.Rmd`
- `adding-a-linkage-process.Rmd`

Each carries, in order: the confirmed file table, marked always/conditional; a worked example
that compiles and fits, with real pasted numbers; which drift guards catch a half-finished job
and what their failure looks like; and troubleshooting — "it compiled and nothing changed", "the
parameter is not estimated", "the value does not come back in the output object", each with the
cause traced in the code.

The knowledge already exists in `.claude/commands/new-column.md`, written for an agent. These
are the same content for a person. The linkage article must state rule 12 (`linkage.hpp` and
`R/0-linkage_encode.R` in lockstep).

Fitting chunks are gated the way `vignettes/*.Rmd` do it — nine of them open with
`eval = identical(Sys.getenv("RCEATTLE_EVAL_VIGNETTES"), "true")`, as do both articles. Do not copy
`vignettes/articles/developer-guide.Rmd` for this; it has one chunk.

**Every recipe ships with its drift guard (below). A recipe without one is not done.**

**Acceptance:** a reader who has never opened `src/TMB/` adds a working form using only the
article. Test it on a named person, not on yourself — pick them at item 0.
**Doc-sync:** `_pkgdown.yml`. No version bump.

### A-guard. Make the recipes self-guarding — 1 day, and the item that decides whether A survives

A nine-file recipe rots on the first pipeline change, and a stale recipe is **worse than none**:
it sends a newcomer to the wrong file with confidence. FIMS's 739-line guide is a standing
maintenance liability and the plan must not import it.

Add `tests/testthat/test-docs-anchors.R`, modelled on `test-schema-cpp-dispatch.R` — which
already does exactly this shape for schema-versus-C++, with pinned exemptions in both
directions. It must assert that every file path, function name and switch code quoted in the
recipes still exists and still resolves in the relevant map.

**Acceptance:** renaming a file the recipe names turns the suite red. Mutation-test it, the way
the existing drift guards were.
**Doc-sync:** tests, exempt.

### B. Symbol glossary — 2–3 days · cheapest thing with broad payoff

`vignettes/articles/symbol-glossary.Rmd` — one table, five columns:
**symbol · code name · units · where computed · where reported.**

FIMS's version maps symbol → meaning. Rceattle's must map **symbol → variable name**, because
that is the actual barrier: someone opening the 5,671-line `ceattle.cpp` cannot tell what
`biomassSSB` or `NByage` hold or in what units. Cover at minimum the quantities named in
`R/6-rename_output.R` and every `JnllRow` row.

Use the `CLAUDE.md` domain vocabulary exactly — SSB is female spawning-stock biomass, F40% is a
max-FABC SPR proxy, not "MSY".

**Acceptance:** every `JnllRow` constant and every reported derived quantity appears; covered by
the A-guard test so it cannot drift.
**Doc-sync:** `_pkgdown.yml`. No version bump.

### C. Publish the Doxygen — half a day · nearly free, immediately visible

A `Doxyfile` plus a step in `.github/workflows/pkgdown.yaml`. The `@brief` blocks already exist
in all eleven headers; this is a build and a link, and it makes the C++ browsable without an editor.

`pkgdown.yaml` triggers on `main` only, so a PR to `dev` gets no CI for it — run
`/pkgdown-check` locally.

**Acceptance:** the site renders `recruitment.hpp` (the canonical Doxygen header per
`CLAUDE.md`) with its equations, and B is reachable from it.
**Doc-sync:** `_pkgdown.yml`. No version bump.

### D. `CONTRIBUTING.md`, split out of `CLAUDE.md` — half a day

`CLAUDE.md` is already the contributor guide, written to an agent. Split it: a human-facing
`CONTRIBUTING.md` (setup, running the suite, branch and PR conventions, the hard rules that
constrain a contribution, links to A), and keep `CLAUDE.md` as the agent's operating manual.
Cross-reference; do not duplicate.

**Acceptance:** no rule stated in full in both files.
**Doc-sync:** repo tooling, exempt.

### E. Scaffolding — 2–3 days · only worth it after A exists

An internal `use_switch_code()` stamping, in one call: the schema row in `R/0-column_schema.R`,
the map entry in `R/0-switches.R`, the `case` stub in the right header, the map/bounds entries,
and a test file. Plus `use_test_template()` writing a `test-<area>-<topic>.R` skeleton with the
`new.env(parent = asNamespace("Rceattle"))` boilerplate already correct — that trap is in
`CLAUDE.md` and costs everyone their first hour.

Internal, `@noRd`.

**Acceptance:** the stamped stub compiles, `test-schema-cpp-dispatch.R` **passes** (the stub
supplies both the map value and the `case`, which is what that guard checks), and the *generated
behaviour test for the new form* fails until the case body is written. A scaffold that leaves
the whole suite green has hidden the remaining work.
**Doc-sync:** internal helpers, no user-visible change.

### G. broom generics — 2–3 days · additive · ⚠ API

`tidy.Rceattle()` — one row per parameter with estimate, SE, bounds. `glance.Rceattle()` — one
row per fit: objective, max gradient, convergence, npar, run time. `augment.Rceattle()` —
observations with fits and residuals.

Ecologists already know {broom}; this lets them use dplyr and ggplot on a fit without learning
the object's structure. `R/0-convergence.R` and `R/6-osa_residuals.R` are the likely sources for
`glance()` and `augment()` — **unverified; confirm before designing against them.**

**Acceptance:** `tidy()` row count equals the estimated-parameter count the existing convergence
code reports, on all four golden models.
**Doc-sync:** NEWS + `DESCRIPTION` `Version:` + `_pkgdown.yml`. Minor bump.

### H. `.devcontainer/` + setup article — 1–2 days

A devcontainer putting the system toolchain first, so the Homebrew clang/gfortran shadowing trap
cannot happen, plus `vignettes/articles/developer-setup.Rmd` for people not using it. This
removes a barrier rather than documenting one — for a collaborator who is not a developer it is
often the difference between contributing and not.

**Acceptance:** a clean container runs `pkgload::load_all(".")` and
`NOT_CRAN=true Rscript -e 'devtools::test()'` green with no PATH intervention.
**Doc-sync:** repo tooling, exempt.

---

## Where item F went

Item F proposed reporting switch names rather than codes, on the premise that
`fit$data_list$fleet_control$Selectivity` holds an integer after `switch_check()` and that
messages therefore quote `8` rather than `"DoubleNormal"`. That premise was wrong, and had been
since before the item was written: `switch_check()` calls `revert_switches()` before returning,
so the stored column holds the canonical string and the messages already quote it. The integer
exists only on the copy `rearrange_data()` hands to TMB.

The item is deleted rather than rewritten, because what it asked for is what the code does.
What survives is the trap above, stated the right way round.

## Explicitly not doing

- **A functor/module rewrite of the C++.** It is the FIMS feature people point at, and it is not
  what makes code legible to an ecologist — naming and documentation are, and those are
  independent of it. Four golden references and three live assessments ride on the current
  template, and FIMS is the cautionary tale: 4.5 years in, "refactor model interface" and
  "refactor model families to dynamic model composition framework" are still open issues against
  a 1.0 milestoned for December 2026.
- **An Rcpp/S4 module interface.** FIMS carries an open issue titled "High Memory Overhead in
  Rcpp Layer" and is planning an XPtr rewrite to fix it. Plain TMB plus the workbook-and-schema
  configuration is simpler and already the better story for this audience.
- **FIMS's 21 CI workflows.** Take the PR checklist; leave the rest.

---

## Prompt to start the work

Item 0 is a conversation, not a session. This opens the first build item; retarget the first
paragraph if item 0 reorders things.

```
Read inst/dev/CONTRIBUTOR-EXPERIENCE.md, then CLAUDE.md and
vignettes/articles/developer-guide.Rmd.

Goal: a fisheries scientist or ecologist with no software background can navigate, edit and
add to this codebase. Build item A -- vignettes/articles/adding-a-selectivity-form.Rmd -- and
its drift guard. Do that one item properly rather than starting several.

FIRST, confirm or correct the file table in the plan's "The file list is a hypothesis"
section. It was assembled by reading, not by tracing, and it is probably wrong somewhere.
Trace one existing form end to end -- Selectivity = 8, "DoubleNormal" (R/0-switches.R:82) --
and report every file I named that is not actually involved and every file I missed. Correct
the always/conditional marking from what you find. Do not take any switch code, name or line
number in that plan on trust; rule 9 applies to the plan as much as to the code.

The article must contain, in order:
  1. the corrected file table, each row with what changes and why, marked always or
     conditional;
  2. the string-vs-integer trap: the column is compared against strings before switch_check()
     (R/1-data_check.R:633) and integers after it (R/0-switches.R:947). Verify this is stated
     correctly and show where a half-finished addition falls through;
  3. a worked example that compiles and fits on a bundled dataset -- write it, run it, paste
     the real numbers. Do not invent output;
  4. which drift guards catch a half-finished job (test-schema-cpp-dispatch.R,
     test-schema-registries.R) and what their failure actually looks like -- run one red;
  5. troubleshooting: "it compiled and nothing changed", "the parameter is not estimated",
     "the value does not come back in the output object". Trace the real cause of each in the
     code; do not guess.

THEN build the guard: tests/testthat/test-docs-anchors.R, modelled on
test-schema-cpp-dispatch.R, asserting every file path, function name and switch code the
article quotes still exists and still resolves in its map. Mutation-test it the way the
existing drift guards were -- rename a file the article names and confirm the suite goes red.
The article is not done without this.

Write for a fisheries scientist, per the Comments section of CLAUDE.md: assessment reasons,
units and conventions, no code narration. Use the domain vocabulary exactly.

Constraints: gate fitting chunks the way vignettes/*.Rmd do
(eval = identical(Sys.getenv("RCEATTLE_EVAL_VIGNETTES"), "true")) and verify with that set to
true -- not developer-guide.Rmd, which has one chunk. Add the article to _pkgdown.yml and run
/pkgdown-check. No version bump; this changes no behaviour. Do not touch src/ or R/ for this
item -- if the trace turns up a real defect, add it to inst/dev/CLEANUP_BACKLOG.md rather than
fixing it unasked.

When it is done, run an adversarial review: hand the article to a fresh agent with no context,
have it try to add a selectivity form using only the article, and report every point where it
had to guess or read source the article did not name. Fix those, then stop and report before
moving on.
```

Adapt the same shape for B–H: name the item, require the facts be derived from the code rather
than from this file, state the acceptance criterion, and end with the adversarial check.
