---
description: Sweep the sibling assessment repos for a symbol, argument, or fleet_control column
argument-hint: "<symbol, argument name, or column>, e.g. rearrange_dat or Sel_norm_scope"
---

Rceattle's API has live consumers. A breaking change here breaks scripts that produce federal
catch advice. Sweep before merging one.

`$ARGUMENTS` is what to look for. If empty, ask.

## Before running

Scope, in order of importance:

- `../Rceattle-models` — EBS/GOA pollock, sablefish, arrowtooth, plaice, POP.
- `../GOA-ATF-ESP` — GOA arrowtooth and its cannibalism run: the only live two-sex,
  `suitMode = 0` model, so it is what exercises the sexed and predation paths.
- `../GOA-multispecies-assessment` — the GOA multispecies assessment (`run_all.R`,
  `R/02_fit_models.R`). `inst/dev/SIBLING-REPOS.md` is the full list.
- The other sibling directories under `../` (`CEATTLE`, `Climate_MSE*`, `hake-CEATTLE`,
  `GOAceattle`, `Rceattle_MSE`, …) — real runs, lower priority. `../Rceattle-models` also
  holds Pacific hake.
- **Ignore `EBS_CEATTLE_TMB`** — a vendored fork, not a consumer.

## Invocation

Grep the three consumers explicitly. Do not `cd` to a parent and filter -- the repo sits one
level below the ecosystem root, and a wrong level sweeps a dozen unrelated repositories while
the `Rceattle/` self-exclusion silently stops matching, so every internal call site is reported
as a consumer hit.

```
grep -rn "$ARGUMENTS" --include="*.R" "../Rceattle-models" "../GOA-ATF-ESP" "../GOA-multispecies-assessment"
```

To widen to the other sibling run directories, add them by name (`../CEATTLE`,
`../Climate_MSE`, `../hake-CEATTLE`, `../GOAceattle`, `../Rceattle_MSE`, ...).

Also count `.xlsx` workbooks if the change touches a schema column or its validation — they are
inputs the change is about to be applied to, and a grep over `.R` will not see them.

## After running

Two things a static sweep cannot do, and you must say so:

- **It catches removed or renamed API. It does not catch behavioural drift.** A script that
  still parses can still produce different numbers.
- **Some models there are partially implemented and do not run regardless**, so not every hit
  needs chasing.

When the sweep is not enough, refit one real assessment and diff. Neither repo caches a fitted
object, so this means a real fit: the terminal fit is under a minute for either pollock model
(`../Rceattle-models/{GOA pollock/2025, EBS pollock/2024}/04-fit-and-diagnostics.R` — their
`Data/` paths are relative to the *project* root, not the year folder), and ~13 min for the ATF
chain (`../GOA-ATF-ESP/R/Run_2025_ceattle.R`, which cannot be sourced straight through on any
version: it references three objects it never assigns, and its `file =` arguments write into
the assessment repo, so run it from a sandbox that symlinks `Data/`).

Report hit counts per repo with file:line, and say explicitly whether each hit still works,
breaks, or needs a human decision.
