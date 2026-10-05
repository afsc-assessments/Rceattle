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

## The five that matter most

1. **Do not invent a switch code, a default, or a unit. Ask.** If it is not stated in
   `R/0-column_schema.R` or a vignette, it is not established. A plausible placeholder that
   survives into a fit is the worst failure mode this repository has.
2. **A change that can move a fit needs the golden regression**, plus the `tools/verify/`
   harness covering what golden does not — and golden does not cover the refit paths, the
   simulation draws, or any figure.
3. **The column schema is the source of truth** for every switch value, default and column
   order. Read them from it; never hardcode one elsewhere.
4. **Preserve the public API.** Deprecate an argument; do not delete it. This package ships and
   has users, and three live assessments consume it.
5. **Verify what ran, not just that it passed.** Guards here have repeatedly reported green
   while measuring less than they claimed — a skipped file, a floor with slack, a regex whose
   match set silently shrank. Count what executed.

## Two facts that cause wrong answers

- **`nages` is a count of age bins, not the oldest age.** Ages run
  `minage .. minage + nages - 1`, and age `a` sits at index `a - minage + 1`. Every bundled
  dataset and all three live assessments use `minage = 1`, which hides the confusion.
- **`$` partial-matches silently, on lists and on data frames.** Where `Time_varying_sel` is
  absent from a `fleet_control`, `fleet_control$Time_varying_sel` returns `Time_varying_sel_sd`.
  Use `[[ ]]` for a new read.

## Scope and style

Match the file you are editing. Keep a change to one concern, and keep reformatting out of any
change that can move a number — a reviewer checking a fit should not have to read past
whitespace to find it. Ask before adding a dependency.

Commit messages are plain and imperative, subject ≤72 characters, with the body saying *why* and
giving the numbers that changed.
