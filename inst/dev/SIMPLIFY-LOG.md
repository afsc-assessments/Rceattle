# Simplification log

Opportunities to clarify or simplify a workflow, an API or a document, found
while working in the code. **Add here; do not act on an entry unasked** -- a struck row records
that something else closed it, not that this file authorised the work. Comment
and roxygen tightening in a file a change already edits goes into that change;
anything that alters an API, a switch, a default or a workflow waits for Grant
to pick it from this list. One line each: where, what, the simpler form, and
whether behaviour changes.

| Where | What | Simpler form | Behaviour change |
|---|---|---|---|
| `fit_mod()` / `run_config()` | `random_sel` and `random_q` are global flags whose only job is a per-fleet decision `fleet_control` already determines (`.rce_np_unintegrable_fleets()`) | derive per fleet; deprecate the flags | yes, needs a deprecation path |
| `save_config()` | writes no `fleet_control`, so the YAML is not the whole model | serialise it | additive |
| `Catchability` / `Time_varying_q` | `Time_varying_q` is a mode, or comma-separated `env_data` column indices under `Catchability = "Environmental"` (`R/3-build_map.R`, `# FIXME: use formula`) | the q linkage, one grammar | yes, needs a shim |
| `fit_mod(initMode =)` | overwrites `data_list$initMode` unconditionally, so a switch stored on the data is never read (`BS2017MS$initMode = 1` never reached the golden fit) | read the data's value when the argument is missing, as `model_config` does | yes; moves `ms` golden |
| `estDynamics` | under `msmMode = 0` codes 1 and 2 estimate nothing, so two codes describe one model | one code in single-species mode, or estimate the scalar there | yes |
| `.map_switch()` vs `.conv()` | two switch normalisers; PR 1 gave `.map_switch()` the same factor and numeric-string handling | one should call the other | internal |
| `switch_check()` `.bad_ed` | the `estDynamics` range check repeats `.map_switch()`'s error, because `.map_switch()` passes numeric codes through for every switch (`msmMode`, `suitMode`, `Catchability`, the sd switches) | validate numerics in `.map_switch()` once each map's legacy numeric codes are known | yes, may reject old numeric codes |
| `build_map_fixed_natage()` | three blocks fix everything a fixed-numbers species would estimate; the linkage rows are refused in `fit_mod()` instead | one header comment saying it is the single place for "what a fixed species estimates" | internal |
| `.osa_bubble_plot()` | one function serves OSA and Pearson residuals, whose null distributions differ | a `type` argument carrying the outlier rule (PR 2) | internal |
| `linkage_spec(species =)` | a spec with no `species =` expands to one row per species, including species whose recruitment is input (`estDynamics > 0`), and is then refused; the estimated species have to be listed by hand | drop the fixed species from the expansion with a message | yes |
| ~~`Sel_norm_scope`, `Sel_norm_bin`~~ | inert unless the fleet is two-sex, normalized at a bin or the max, and not Hake/LogisticPM; three conditions the reader has to assemble | **done in 5.35.0**: one sentence in the schema (PR 1) and a `data_check()` note (PR 1) | doc |
| ~~`R/0-osa_data.R:79-85`~~ | five-line comment plus a FIXME for a one-line fact | **done in 5.36.0**: four lines, no FIXME, all three fill sites named | doc |
| ~~`CONTRIBUTOR-EXPERIENCE.md` item F~~ | proposed work the code already does (`revert_switches()`) | **done**: deleted, with a section recording why the premise was wrong | doc |
| `R/10-run_mse.R:901`, `:1126` | every deviation array but `log_M1_dev` is carried into the projection, in BOTH the operating model and the estimation model's refit | carry the terminal year, as `index_q_dev` is; see `CLEANUP_BACKLOG.md` for both sites | yes, `M1_re` models only |
| `switch_check()`, `data_check()` | the per-form lists (which forms read `Sel_curve_pen*`, which take `Sel_norm_bin`, which allow which `Time_varying_sel`) are inline literals repeated at each site; no predicate function says "this form is parametric / non-parametric / takes deviates" (found tracing DoubleNormal for `adding-a-selectivity-form.Rmd`) | one registry per form (parameters read, deviate modes, penalty columns) that `build_map_selectivity()`, `data_check()` and `.PAR_SEL_SLOTS` all read | internal |
| `R/6-osa_residuals.R` | `method` and `discrete` are resolved per observation inside one call, and the resolution is reported only in `attr(osa, "method")` and `attr(osa, "discrete")` | a `method` column on the returned rows, so each residual says how it was computed | additive |
| `R/` internal helpers | 27 `.helper()` roxygen blocks use `@keywords internal` where CLAUDE.md asks for `@noRd`, so each generates a `man/dot-*.Rd` nobody links to (pkgdown drops them by keyword, so nothing breaks; the churn is the cost) | one sweep to `@noRd`, deleting the 27 pages | doc |
