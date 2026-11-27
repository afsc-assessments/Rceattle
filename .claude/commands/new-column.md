---
description: Add a fleet_control / control column end to end, so it cannot round-trip to nothing
argument-hint: "[column name and what it controls, e.g. Sel_cap_bin: cap on NonParametricPM bins]"
---

Add a new column to the model's input schema. This is the operation most often left half-done:
a column that works in your session but is silently lost through `write_data()`/`read_data()`,
or absent from `write_template()`, or documented in a vignette that disagrees with the schema.
This is how `index_cov` was lost (CLAUDE.md, hard rule on the schema).

`$ARGUMENTS` names the column and what it controls. If it does not, ask before writing anything.

## Before running

Stop and ask rather than guess if any of these is unclear:

- **Can the behaviour be derived from an existing field instead?** `Fleet_type` already
  distinguishes fisheries, surveys and Off fleets. Grant pushes back on column proliferation;
  a derived value needs no schema row, no round-trip, and no default. Say what you considered.
- **What is the unit and the scale?** State it. A weight that is a log under one distribution
  and natural-scale under another (`Comp_weights`) needs that in the `doc` string, not just in
  a vignette.
- **Is it a bin index? The two conventions are opposite, so say which (rule 10).**
  `Age_first_selected` is an absolute **age** on any fleet, and `Sel_norm_bin`,
  `Sel_norm_bin_upper`, `Sel_pen_first_bin`, `Sel_pen_last_bin` and `Sel_cap_bin` are an
  absolute **age** on an age-based fleet, so `rearrange_data()` subtracts the `minage` offset
  (`sel_bin_offset`); on a length-based fleet they are 1-based length-bin ordinals and the
  offset is 1. `Bin_first_selected` is a **1-based length-bin ordinal** on a length-based fleet
  and the deprecated spelling of `Age_first_selected` on an age-based fleet, resolved via
  `.rce_first_selected()`. `N_sel_bins` is neither; it is
  a count. The cpp penalties loop over `nbins`, not `nages`.
- **Is it read per fleet inside a shared parameter block?** Fleets sharing a
  `Selectivity_index` / `Catchability_index` share ONE block, and `adjust_map_shared_params()` copies the
  donor fleet's slice over the rest — so a per-fleet setting that differs within a group is
  silently taken from the donor. If the column can differ within a group, say what happens.

## The order to do it in

1. **`R/0-column_schema.R`** — one `.rce_col()` row. This is the single source of truth: it
   drives `switch_check()` defaults, `write_data()`/`read_data()` ordering, the meta sheet, the
   field dictionary, and the alias upgrade of older spellings. Set `doc` (a real sentence, with
   the unit), `type`, `has_default`/`default`/`default_msg`, `aliases` for any older spelling,
   and `tmb_target` if it reaches the template. **Do not hardcode the default anywhere else.**
2. **`R/0-switches.R`** — if it is a `type = "switch"`, add its map, set `allowed` on the schema
   row, and wire it through `validate_switches()`. Apply the default with
   `.rce_apply_default(..., .sch)`, reading the schema — not a literal.
3. **`R/5-rearrange_data.R`** — the R→TMB rename, if it reaches the template. Then the
   `DATA_*` declaration in `src/TMB/ceattle.cpp` and the code that reads it. Remember TMB source
   is inert until `pkgload::load_all(".")`.
4. **`R/0-read_write_excel_data.R`** — confirm `write_data()` and `read_data()` round-trip it,
   and that `write_template()` emits it. A new element without these is silently lossy.
5. **`R/data.R`** — the `@format` `\item{}` entry. `test-schema-canonical.R` compares these
   against the schema `doc`.
6. **`data-raw/regenerate_meta_xlsx.R`** — re-run if the column is `meta = TRUE`.
7. **The vignette switch table** in `vignettes/model-options-and-functionality.Rmd`, and
   `model-parameterizations.Rmd` if it changes an equation.
8. **`NEWS.md` + `DESCRIPTION` `Version:`** — minor for a new column.

## After running

- `/test` — `test-schema-canonical.R` must pass; it is the drift guard for steps 1 and 5.
- Round-trip it for real: `write_data()` → `read_data()` → `expect_identical()`. Do not take
  "it appeared in the workbook" as proof.
- `/document` — check `git diff DESCRIPTION` FIRST; if the roxygen version key moved, the
  `man/` churn is the version, not your change.
- If the column can change a fit: `/golden-check`. If it can reach a refit path,
  `tools/verify/verify-refit-like.R` as well — golden does not cover `.refit_like()`.
- `/ecosystem-sweep` for the column name, in case an assessment already uses that spelling for
  something else.

Report which of the eight steps applied and which did not, and why. A step skipped silently is
the failure mode this command exists to prevent.
