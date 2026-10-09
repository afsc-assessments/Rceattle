# Copilot instructions for Rceattle

Rceattle fits CEATTLE, a single- and multi-species, climate-linked, age-structured stock
assessment model. Its output sets US federal catch limits under Magnuson-Stevens. A
silently-wrong number here becomes a wrong quota, and it will not announce itself as a crash.
Accuracy beats speed, every time.

This file is a summary. The guidance in full is `AGENTS.md` at the repository root (tool-neutral,
read it first), `CONTRIBUTING.md` (setup, tests, what constrains a change) and `CLAUDE.md` (the
hard rules, the verified traps with their measured numbers, the domain vocabulary to use in
plots, documentation and messages).

## The three that govern every change

1. **Do not change an API or a workflow without asking.** Deprecate an argument; never delete
   one; and do not rename, re-default or re-mean one either. An argument name or default, a
   switch code or its meaning, a column in `R/0-column_schema.R`, the names or shape of anything
   in `fit$`, a message an assessment script reads, how a figure looks: if a user's existing
   script would behave differently, say so and ask rather than doing it — even when the
   change is
   plainly an improvement, and even when you are already in the file for another reason. File the
   idea as one row in `inst/dev/SIMPLIFY-LOG.md` and finish the task you were given — though a
   silent wrong number is a defect, not a wish, and belongs in
   `inst/dev/CLEANUP_BACKLOG.md` or an issue. The package ships, and the repositories in
   `inst/dev/SIBLING-REPOS.md` consume it. What this rule is
   *not*: a refactor behind an unchanged API is free, and so is a defect fix — refusing a
   configuration the model never actually fitted, correcting a number the model computed
   wrongly, or freeing a parameter it fused, mapped out or started badly. Make those, and say
   whose fits move. Deprecating never deletes, so a *new* name is permanent too: agree the name
   of a new exported function, argument or column before it ships.
2. **Write the simplest thing that works, and nothing for later.** No helper, wrapper, class,
   config layer, option flag or *new* registry the change does not need now — but registering in
   the registries that already exist (the `JnllRow` partners, `.JNLL_ROW_AXIS`,
   `.index_rows_natural_scale()`) is not optional, or the new family gets the wrong
   residual. A helper earns its place at two callers, or when it names a concept the
   reader needs. A guard that cannot fire is
   not safety: name the input that reaches it, or leave it out. Keep a change to one concern, and
   never put a formatting-only hunk in a change that can move a fit.
3. **Write for a fisheries scientist, briefly.** A comment is one or two lines — the
   assessment reason, the units, the convention — and a `@param` is one sentence. Anything
   longer belongs in
   `@details`, a vignette, or `inst/dev/`. A literature citation is the specification, so never
   strip one. State current behaviour, never bug history, and check that every function and
   argument a comment names still exists before you shorten it.

The first rule is not licence to ask about everything. Settle ordinary implementation choices
yourself, and stop at anything that changes what the code *means* or what a user's script does.

## Do not guess

If a switch code, a default or a unit is not stated in `R/0-column_schema.R` or a vignette, it is
not established: ask. A plausible placeholder that survives into a fit is the worst failure mode
this repository has. Two facts cause wrong answers on their own — `nages` counts age bins, not
the oldest age (ages run `minage .. minage + nages - 1`), and `$` partial-matches silently on a
list *and* a data frame, so use `[[ ]]` for a new read.

## Before suggesting a change is finished

R follows the [tidyverse style guide](https://style.tidyverse.org), C++ the
[Google C++ style guide](https://google.github.io/styleguide/cppguide.html), as tie-breakers only
— the surrounding file and TMB/Eigen idiom win. Never hand-edit `man/*.Rd` or `NAMESPACE`; run
`devtools::document()`. A change that can move a fit needs the golden regression
(`tests/testthat/test-golden-regression.R`) plus the `tools/verify/` harness that covers what
golden does not. A behaviour, API or documentation change updates `NEWS.md`, the `DESCRIPTION`
version and the affected vignette in the same commit.
