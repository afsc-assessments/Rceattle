# AGENTS.md

Rceattle fits CEATTLE, a single- and multi-species, climate-linked, age-structured stock
assessment model. Its output sets **US federal catch limits under Magnuson-Stevens**. A
silently-wrong number here becomes a wrong quota, and it will not announce itself as a crash.

This file is the tool-neutral entry point. The full guidance lives in two places and is not
repeated here:

- **`CONTRIBUTING.md`** — how to set up, run the tests, and what constrains a change.
- **`CLAUDE.md`** — the hard rules in full, the known traps with the measured numbers behind
  them, and the domain vocabulary to use in plots, docs and messages.

Read `CONTRIBUTING.md` first. If you are changing anything numeric, read `CLAUDE.md`'s "Hard
rules" and "Known traps" before you start.

`.github/copilot-instructions.md` is a summary of this file at the path GitHub Copilot reads.
This file, that one and `CLAUDE.md` each state the three doctrines in full rather than pointing
at one another, because each is injected by its own tool. They are therefore hand-synced, and
each carries its own examples and emphasis — what they must share is the set of operative clauses
pinned in `tests/testthat/test-docs-doctrine-sync.R`, which fails if one of them loses a clause
the others keep. Change a doctrine and that test tells you which files are now out of step.
`CONTRIBUTING.md` is different: a human can follow a link, so it points here instead.

## The five that matter most

1. **Do not invent a switch code, a default, or a unit. Ask.** If it is not stated in
   `R/0-column_schema.R` or a vignette, it is not established. A plausible placeholder that
   survives into a fit is the worst failure mode this repository has.
2. **Do not change an API or a workflow without asking.** Deprecate an argument; never delete
   one; and do not rename, re-default or re-mean one either. An argument name or default, a
   switch code or its meaning, a schema column, the names or shape of anything in `fit$`, a
   message an assessment script reads, how a figure looks: if a user's existing script would
   behave differently, stop and ask, even when the change is plainly an improvement. Adding a
   message is additive; changing or removing one is not. File the idea as one row in
   `inst/dev/SIMPLIFY-LOG.md` and finish the task you were given. This package ships and has
   users, and the repositories in `inst/dev/SIBLING-REPOS.md` consume it. Deprecating never
   deletes, so a *new* name is permanent too: agree the name of a new exported function,
   argument or column before it ships.
   What this rule is **not**: a refactor behind an unchanged API is free, and so is a defect fix
   — refusing a configuration the model never actually fitted, correcting a number the model
   computed wrongly, or freeing a parameter it fused, mapped out or started badly. Make those,
   and say which fits move and which scripts stop. A silent wrong number is never a
   `SIMPLIFY-LOG.md` row: that file holds design wishes, while a known defect belongs in
   `inst/dev/CLEANUP_BACKLOG.md` or an issue.
3. **A change that can move a fit needs the golden regression**, plus the `tools/verify/`
   harness covering what golden does not — and golden does not cover the refit paths, the
   simulation draws, or any figure.
4. **The column schema is the source of truth** for every switch value, default and column
   order. Read them from it; never hardcode one elsewhere.
5. **Verify what ran, not just that it passed.** Guards here have repeatedly reported green
   while measuring less than they claimed — a skipped file, a floor with slack, a regex whose
   match set silently shrank. Count what executed.

## Two facts that cause wrong answers

- **`nages` is a count of age bins, not the oldest age.** Ages run
  `minage .. minage + nages - 1`, and age `a` sits at index `a - minage + 1`. Every bundled
  dataset uses `minage = 1`, and so do the live assessments, which hides the confusion.
- **`$` partial-matches silently, on lists and on data frames.** Where `Time_varying_sel` is
  absent from a `fleet_control`, `fleet_control$Time_varying_sel` returns `Time_varying_sel_sd`.
  Use `[[ ]]` for a new read.

## Scope and style

Match the file you are editing. Keep a change to one concern, and keep reformatting out of any
change that can move a number — a reviewer checking a fit should not have to read past
whitespace to find it. Ask before adding a dependency.

**Write the simplest thing that works, and nothing for later.** No helper, wrapper, class,
config layer, option flag or *new* registry this change does not need now — but registering in
the registries that already exist (the `JnllRow` partners, `.JNLL_ROW_AXIS`,
`.index_rows_natural_scale()`) is not optional, or the new family gets the wrong
residual. A helper earns its place at two callers, or
when it names a concept the reader needs. A guard that cannot fire is not
safety: name the input that reaches it, or leave it out. If the structure of a change needs a
paragraph to justify it, propose it before writing it.

**Write for a fisheries scientist, briefly.** A comment is one or two lines — the assessment
reason, the units, the convention — and a `@param` is one sentence. Anything longer belongs in
`@details`, a vignette, or `inst/dev/`. A literature citation is the specification, so never
strip one, and a regression test keeps its provenance. State current behaviour, never bug
history; but check that every function and argument a comment names still exists before you
shorten it, because dropping the past tense can assert an API that no longer works.

Settle ordinary implementation choices yourself and stop at anything that changes what the code
*means* or what a user's script does. Commit messages are plain and imperative, subject ≤72
characters, with the body saying *why* and giving the numbers that changed, and no
`Co-Authored-By` trailer.
