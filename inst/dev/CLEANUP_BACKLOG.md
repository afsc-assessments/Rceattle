# Cleanup backlog

The `TODO` / `FIXME` markers in the source, triaged. **Add to this file; don't fix these
unasked.** Fix the one in the file you were already asked to touch, in the same commit.

**Cite the marker text, not a line number.** These references have gone stale three times --
twice when roxygen was added above them, once from an unrelated `ceattle.cpp` edit. Grep for the
quoted FIXME text instead; it moves with the code.

Three tiers: a **known defect** is a wrong answer waiting for the right input and should become a
GitHub issue; a **design note** is a wish, not a bug; `TODO(review)` is a deliberate convention
marking a judgement call for Grant, and is never resolved by an agent.

**47 remain**, re-derived with the command below rather than carried forward — the old running
total had drifted. The chain: 57 at 5.41.0, **56** on `dev` at 5.49.1 (`b066015f`), then **−7** for
`src/TMB/Dev/` and **−2** for the Kinzey & Punt blocks, both deleted as dead code
(`REMOVED-kinzey-and-dev-cpp.md`). Note `src/TMB/Dev/` held 7 markers, not the 5 this file
previously attributed to `caal.hpp` alone.
Counts by area: `src/TMB/ceattle.cpp` 23 ·
`src/TMB/predation.hpp` 5 · `R/10-run_mse.R` 4 ·
`R/3-build_map.R` 3 · `R/9-retro_and_jitter.R` 3 · `R/0-rceattle_class.R` 3 ·
`R/6-fit_mod.R` 2 · `src/TMB/growth.hpp` 1 · rest 1. Re-derive with
`grep -rnE 'TODO|FIXME' R/ src/TMB/ | grep -v 'todo <-' | grep -v 'TODO-'` -- the `-E` is
needed for the alternation, the first filter drops a variable in `R/6-process_residuals.R`, and
the second drops pointers to `inst/dev/TODO-*.md` notes, which are not markers.


## How to work one

Absorbed from `BACKLOG-PLAN.md`, which this file replaces.

1. **Reproduce first, in a test that fails.** The marker names the triggering input; build the
   fixture that reaches it. A fix whose test passes before the fix is not a test.
2. **Check what actually covers it.** `/golden-check` will be green either way for almost every
   item here -- none of the four reference models reaches these inputs. Use `/verify` to pick the
   right `tools/verify/*.R` harness, and remember no harness reaches
   `sample_rec(update_model = TRUE)`, `reweight_comps()`, or any figure.
3. **For a C++ change**, recompile before testing (`pkgload::load_all(".")`) and run the suite
   serially (`TESTTHAT_PARALLEL=false`). Golden is required even when you expect no movement.
4. **Adversarially review the diff before committing.** Across this work every review found
   something real, including two changes that moved fitted numbers and would otherwise have
   shipped.
5. **NEWS + `DESCRIPTION` + the affected vignette, same commit.** `/doc-sync` checks it.
6. **Move the item to "Deliberately not changed" or delete it** -- a backlog that only grows
   stops being read.

**The tiers are not a priority order.** Sequence by who is exposed: a defect a live assessment
can reach outranks a tidier one a fixture cannot. Two of these produce a wrong number rather
than an error, which is the failure mode this package cannot afford, and they come first
whatever tier they sit in.

**A Tier 0 claim is a claim about behaviour, so check it against the code, not just the
comment.** One entry named the wrong switch (`Time_varying_q` instead of `Catchability`) for a
whole session, and the source comment goes out of its way to warn about that confusion.

---

## Tier 0 — known defects, with the input that triggers them

These say, in the source, that the code is wrong under a stated condition.

Three further defects of the same class were found reviewing the fixes below, and are resolved
in 5.13.0 alongside them. None carried a marker, which is why none appeared in this file: they
are what the markers pointed *near*, not what they said.

Rows marked **Open** were found across several reviews, from PR #143 (2026-09-12)
onward. Most carry a measurement from a real fit; the few that do not say so in the row.

| Where | Condition | Consequence |
|---|---|---|
| ~~`R/3-build_map.R` (`# Age-independent scalar`)~~ | ~~`estDynamics = 3` in a multispecies model (`msmMode > 0`)~~ | **Resolved in 5.35.0**: code 3 is retired (it always fitted as 2) and `log_pop_scalar` is one value per species. `test-switches-estdynamics-retired.R`. Original note: Gated on `estDynamics[sp] == 2 \| msmMode != 0`, so under `msmMode > 0` columns 2..`nages` of `log_pop_scalar` are mapped out for every species, and the age-specific scalar of `estDynamics = 3` is never estimated (it stays `exp(init)`, 1 without `inits`). Gating on `estDynamics[sp] == 2` alone is safe for codes 0-2. Check the Hessian before freeing them: they are informed only through predation and that species' index. The `# Age-dependent scalar` branch has the same `\|` but maps only padding beyond `nages`, which is harmless. |
| ~~`R/3-build_map.R` (`# Don't estimate the scalar`)~~ | ~~`msmMode = 0`~~ | **Resolved in 5.35.0**: the schema, `?BS2017SS` and `data_check()` now say `estDynamics = 2` fits as 1 in single-species mode. Whether to merge the two codes there is in `SIMPLIFY-LOG.md`. |
| ~~`JNLL_Q_PRIOR` against a q linkage intercept prior~~ | ~~`Catchability = "Estimated-with-prior"` plus a q linkage `(Intercept)` prior on the same fleet~~ | **Resolved in 5.33.0**: both penalized that fleet's log q, so the prior counted twice. `.check_q_linkage_support()` now refuses the pair, taking fleet 1 for a row with no fleet, as the template does. `test-linkage-double-prior-guards.R`. |
| ~~`JNLL_M_PRIOR` against an M1 linkage intercept prior~~ | ~~`M1_use_prior = TRUE` with `M2_use_prior = FALSE`, plus an M1 linkage `(Intercept)` prior~~ | **Resolved in 5.33.0**: both penalized that species' log M1. `.check_M_linkage_prior()`, run from `fit_mod()`, now refuses the pair, taking species 1 for a row with no species. `test-linkage-double-prior-guards.R`. |
| `R/3-build_map.R` (`if (sel_type == "DoubleNormal")`) | **Open** (found tracing the form for `adding-a-selectivity-form.Rmd`, 2026-09-16). `Selectivity = "DoubleNormal"` with `Time_varying_sel = "RandomWalkAscending"` (5), or any `Time_varying_sel` the branch does not name | The branch handles `IID`, `AR1`, `RandomWalk` and `Block` only, and `data_check()` restricts `Time_varying_sel` per form for NonParametric, NonParametricPM, Hake and LogisticPM but not DoubleNormal, so the deviates are silently mapped out and a static curve fits where the workbook asked for a time-varying one. Measured on `GOApollock` with the fishery (shipped `DoubleLogistic` + `RandomWalkAscending`, 316 parameters) switched to `DoubleNormal`: 220 parameters, no message. Refuse the combination in `data_check()` next to the per-form checks, or implement the ascending-only walk for the peak and ascending width. |
| `R/2-build_params.R` (`sel_inf` starting values) | **Open** (found 2026-09-16, same trace). `Selectivity = "DoubleNormal"` fitted from the default starting values | DoubleNormal reuses the logistic slots, so its peak starts at `sel_inf[1] = 0` (below the first age) and its right-tail floor at `sel_inf[2] = 10` on the logit scale (a floor of 1): the starting curve is flat at 1 for every age, the ascending width has no gradient, and the optimizer stays on that ridge. `GOApollock` fishery, static selectivity, phased fit: objective 3085.98 with selectivity 1.000 at every age from the defaults, against 914.10 (AIC 2268 vs the double-logistic's 2276) from `inits` with the peak at age 4, a logit floor of 0 and widths of 2 ages. `test-selectivity-double-normal.R` sets its own starts, which is why the suite does not see this. The only form-specific start in `build_params()` is LogisticPM's; add DoubleNormal's (peak mid-range, floor near 0). |
| `ceattle.cpp` 5.13 (`SIMULATE PROCESS ERROR`) | **Open** (found 2026-09-17 while adding `NonParametricIntegrable`). `sim_mod(process = "selectivity")` on any fleet whose `Time_varying_sel` deviates are scored in `JNLL_SEL_DEV` (`sel_coff_dev`, `log_sel_slp_dev`, `sel_inf_dev`; forms 1, 2, 3, 5, 8, 13) | Slot 4 of `simulate_state` is only consumed by the linkage random effects (5.12b): the `Time_varying_sel` deviates have no `SIMULATE` draw beside their density, so a "redraw selectivity" request keeps the fitted deviates and a self-test measures recovery of those deviates, not of the process. `tools/verify/verify-sim-recovery-np-integrable.R` draws them in R instead. Add the draws in 5.13 gated on `simulate_state(4)`, per form (iid about 0 for IID; increments for the walks), and report them as `*_sim` for `attr(x, "process_sim")`. |
| `R/0-column_schema.R` (`type = "switch"`) not enforced at the boundary | **Open** (found reviewing the 5.34.0-5.43.0 release, #158). Any `fleet_control` that has not been through `switch_check()` -- `data_check()` is callable on one, and `rearrange_data()` is exported | The schema types **thirteen** columns as `switch` with an `allowed` map (`Fleet_type`, `Selectivity`, `Time_varying_sel`, `Catchability`, `Comp_distribution`, `estDynamics`, ...), but nothing applies those maps on entry, so every comparison written against the canonical spelling is wrong on the integer form the workbook stores -- and every bundled data set stores the integer form. `0 != "Off"` is `TRUE`; `0 == "Off"` is `FALSE`; `0 != "Fishery"` is `TRUE`. The `Fleet_type` half of one line was fixed at 5.43.0 (`est_sel_flts`). **Its `Selectivity` half is still raw, on the same line, and still has a demonstrated effect** (measured reviewing the 5.43.0 delta, #158): `R/1-data_check.R:1255` reads `fc$Selectivity != "Fixed"`, and `0 != "Fixed"` is `TRUE`, so on raw `GOA2018SS` fleets 4 and 5 -- `Selectivity = 0` (Fixed), each on its own `Selectivity_index`, no comp or CAAL rows at `Year > 0` -- are still named by "estimated Selectivity but no comp_data". Canonicalizing both columns names nobody. Note before fixing it: those two fleets are what satisfies the positive control in `test-switches-fleet-type-integer-off.R`, so that assertion needs a fleet with a genuinely estimated form and no comps. The other sites are unreachable **only because their callers canonicalize first**, which is a property of the call graph, not of the code: `R/3-build_map.R:1425` (`flt_off <- Fleet_type == "Off"`) picks the DONOR ROW for a shared selectivity block and its own comment says getting it wrong "would silently stop estimating their selectivity/catchability"; `R/3-build_map.R:1570` maps out comp/CAAL weights the same way; `R/2-build_params.R:157` and `R/3-build_map.R:1554` use `!= "Fishery"`, which is TRUE for every fleet on an integer column. **Fix the class, not the instances:** canonicalize every schema `switch` column once at `data_check()`'s entry through its declared `allowed` map. That changes which errors fire for raw workbooks across thirteen columns, so it wants its own PR and a golden run. Converting sites one at a time was tried and rejected in #158: it left one file with two conventions and installed a third resolver disagreeing with `.canon_switch()` on `" 0 "`, `"0.0"` and `"00"`. |
| ~~`R/1-data_check.R` (`est_sel_flts <- ...`)~~ / `Fleet_type` read raw | ~~An `NA` `Fleet_type`~~ / **Open:** an integer-coded `Fleet_type` (`0` for Off) read before `switch_check()` canonicalizes | **The NA half is resolved in 5.43.0**: `switch_check()` now refuses a blank `Fleet_type`, naming the fleet, before anything reads the column, so the all-`NA` row that killed `data_check()`'s `vapply` with `missing value where TRUE/FALSE needed` cannot form. `test-switches-fleet-type-blank.R`. **The integer half is open**: `0 != "Off"` coerces to `"0" != "Off"`, which is `TRUE`, so a fleet the workbook marks Off reads as LIVE at every bare `!= "Off"` comparison on a `fleet_control` that has not been through `switch_check()` -- `data_check()` is callable on one, and `validate_switches()` says so in as many words and canonicalizes via `.canon_switch()` first. The remaining raw comparisons do not. Audit them with `grep -n '!= "Off"' R/` and route each through `.canon_switch()` or `%in% c(0, "0", "Off")`. |
| `.group_lead()` (`R/5-rearrange_data.R`) | **Open** (found reviewing the 5.34.0-5.42.1 release, #158). A fleet with an NA `Selectivity_index` | `data_check()` accepts NA there and the column has no schema default, but `.group_lead()` is fed `paste(Selectivity_index, form)`, which turns NA into the string `"NA"` -- so its own `lead[is.na(key)] <- 1L` guard is dead code, and two NA-index fleets sharing a form are grouped together with only the first leading. Measured on `BS2017SS` with `Selectivity_index[1:2] <- NA`: `data_check()` passes, `flt_sel_lead` is `1,0,1,...`, yet `.shared_block_lead()` returns NA for fleet 2, so fleet 2 estimates its OWN `sel_coff` block whose shape and curvature penalty is never charged -- the mirror image of the double-penalty row above. Either give `Selectivity_index` a schema default of `Fleet_code`, or refuse NA in `data_check()`. |
| `rearrange_data()` `flt_type` (`.off`) | **Open, low** (same review). `Fleet_type` held as a factor rather than a canonical string or integer | `rearrange_data()` computes `.off <- flt_type == 0` on the factor's LEVEL INDICES, so an `"Off"` fleet reads as non-zero and can lead its group, while every R-side check reads it correctly through `.canon_switch()`. The two then disagree about which fleet is the lead. `switch_check()` converts the column before either is reached, so this needs a hand-built `fleet_control`; the fix is to canonicalize in `rearrange_data()` as well. |
| `adjust_map_shared_params()` vs `flt_sel_lead` | **Open** (found reviewing the 5.34.0-5.42.1 release, #158). Two live fleets sharing a `Selectivity_index` with **different** `Selectivity` forms | The map shares one `sel_coff` block across the group (`adjust_map_shared_params()` keys on the index alone), while the template makes them two penalty groups (`flt_sel_lead` keys on index AND form) and sets the lead flag on both -- so the single shared block is charged its shape and curvature penalty **twice**, once per fleet. `data_check()` only warns about a mixed-form group (the `.sel_shaping_cols` check), so the configuration fits. No bundled data set and no consumer-repo workbook has one (11 data sets, 375 workbooks, 183 with a `fleet_control` sheet), so nothing measures it today. The 5.42.1 fix to the negative-weight guard (`.rce_sel_pen_lead()`) made the guard agree with the template but did not reconcile the two grouping rules; see `TRAPS.md`, "Shared parameter blocks". Decide which rule is right and make the other follow it, or refuse a mixed-form group outright. |
| `srr_terms_on` | **Open, low.** A recruitment linkage intercept prior on an `estDynamics > 0` species | The gate covers the stock-recruit prior and the curve penalty, not the linkage-prior loop. `rec_pars` is mapped out for such a species, so the prior only adds a constant: it moves the objective and `JNLL_LINKAGE_PRIOR` (likelihood tables, AIC), but no estimate. Not yet checked: slope rows on such a species, which `build_map_fixed_natage()` does not map out. |
| `// Input SB0 (if running in multi-species mode)` | **Open.** `msmMode = 0` with `estDynamics > 0` | Equilibrium `SB0`/`SBF` are still built on placeholder recruitment; from 5.30.0 only the dynamic runs keep the input numbers. HCRs 5, 6 and 7 read `SB0` when `DynamicHCR = FALSE` (5 also reads `SBF`), and with `DynamicHCR = FALSE` `ssb_depletion` and `biomass_depletion` divide by `SB0` and `B0` under every HCR. Projected numbers come from `NByageFixed` and are right; the reported depletion, and F and catch advice under those HCRs, are not. |
| `ceattle.cpp` (`Type R_curve = calculate_recruitment(`) | **Open, low** (found in a later review, 2026-09-13). Ianelli form with an identity-link linkage offset on alpha or beta that drives the curve to zero or below in a year before `srr_mse_switchyr` | The dynamic-B0 deviation `log R - log R_curve` is then NaN, and the recursion carries it into `DynamicSB0` for every later year; under `DynamicHCR = TRUE` it reaches depletion, the HCR and the objective. In penalty years the stock-recruit penalty is already NaN, so this is new only outside `srr_hat_styr`..`srr_hat_endyr`. Guard the curve with `posfun()`, or refuse identity-link offsets that can make it non-positive. |
| ~~`int spawn_yr = yr - minage(sp);`, `int rp_yr = yr - minage(sp);`~~ | ~~`minage = 0` with a stock-recruit curve~~ | **Resolved in 5.33.0** by refusal: each lag read the same year's `ssb`/`SB0`/`DynamicSB0` before it was accumulated, so the curve gave R = 0 in the hindcast (`srr_fun >= 2`), the reference points (`srr_pred_fun >= 2`) and a curve-based projection. `data_check()` now refuses `minage = 0` with any curve; computing SSB before recruitment would allow it, but the GOA and AI cod `minage = 0` models fit no curve. `test-data-check-srr-guards.R`. |
| ~~`src/TMB/ceattle.cpp` (male slot writes)~~ | ~~every species one-sex (`max_sex == 1`)~~ | **Resolved in 5.13.0**: ten lines wrote sex index 1 unconditionally, but arrays are dimensioned `max_sex`, so that index does not exist when no species has two sexes. Value written is 0 (`sex_ratio` is set to 1 for a one-sex species first), but the write is out of range and lands on `(sp, 0, age + 1, yr)` — the next age — surviving only because the age loop overwrites it. Fires on BS2017SS and BS2017MS every evaluation. Reproduced with `TMB::compile(safebounds = TRUE)`, which raises Eigen's range assertion; guarded, the fit is clean at an unchanged 1537036.287629372. `test-dynamics-sex-index-bounds.R`. |
| ~~`R/1-data_check.R` (no `comp_data$Sex` check)~~ | ~~`Sex` 2 or 3 on a one-sex species~~ | **Resolved in 5.13.0**: `M1_base`, `weight` and `ration_data` are all checked against `nsex`; composition was not. Two registries disagree on what "joint" means — `check_composition_data()` uses `nsex == 2 & Sex == 3`, the template uses `flt_sex == 3` alone — so a joint row on a one-sex species was sized at `nages` and written to `nages * 2`, corrupting the NEXT observation's predicted composition and its likelihood. Refused at the boundary rather than reconciled in the template. `test-data-check-comp-sex.R`. |
| ~~`src/TMB/ceattle.cpp` (reference-point recruitment arms)~~ | ~~a stock-recruit curve with `proj_mean_rec = TRUE` (the default)~~ | **Resolved in 5.13.0**: the mean-rec arm required `proj_mean_rec == 1 & srr_pred_fun < 2` and the curve arm required `proj_mean_rec == 0`, so that combination matched neither and reference-point recruitment stayed 0 after year 1. `SB0` became the initial cohort decaying (3.344 → 1.230 over six years), and `SB0` in the terminal year is what HCR 5 and 6 read as the depletion reference — so perceived depletion and the resulting catch advice were both wrong. The curve arm now fires whenever a curve exists; the projection switch is read separately and is unchanged. `build_srr(proj_mean_rec =)` was also documented backwards. `test-dynamics-refpoint-mean-rec.R`. |
| ~~`R/3-build_map.R` (`build_map_m1`, the `map_list$log_M1[sp, , 1:nages_sp]` writes)~~ | **Resolved** (found reviewing PR #197, 2026-10-04). `M1_model >= 1` on a species set with ragged `nsex`, i.e. any one-sex species in a model whose `max(nsex)` is 2. | **The starting M1 is diluted toward 1.0 per year by the padding sex cells.** `build_map_m1()` writes `[sp, , 1:nages_sp]` -- a bare `,`, so ALL sexes -- putting a one-sex species' padding sex-2 cells in the same map level as its real cells. `TMB:::updateMap()` is `tapply(parameter.entry, map.entry, mean)`, so the shared start is the MEAN over the level, and the padding sits at `log(1) = 0`. Measured on `GOA2018SS` at `M1_model = 1`: pollock starts at 0.637454 against an intended 0.406348 (**+57%**, and exactly `sqrt(0.406348)`), cod at 0.710571 against 0.504911 (**+41%**). `M1_model = 3` is the same, giving `1.178983, 0.830662, 0.692820` instead of `1.39, 0.69, 0.48`. Species 2 is correct: its two REAL sexes legitimately share one level under a sex-invariant model. **Pre-existing and unchanged by #197** -- `build_params()` has always left `1` in those cells, so the dilution is identical on `dev`. It is a START value, not the likelihood (the template reads `sex < nsex(sp)`), so the MLE stays identified; but `test-golden-regression.R`'s own header records `goa_ss` reaching a second minimum 52.9 nats up from a one-ULP perturbation, so a 40-57% wrong start is not cosmetic. **Fix:** write `[sp, seq_len(nsex_sp), 1:nages_sp]` and leave the single-sex fallbacks at `NA`; the reviewer verified that gives exact starts with unchanged free-parameter counts (3 / 4 / 64), and all four golden references run `M1_model = 0`, so it is free. **RESOLVED on `fix/m1-map-padding-dilution` (5.51.0).** The effect on an UNPHASED fit was then measured: fitting `GOA2018SS` at `M1_model = 1` cold with `phase = FALSE, newtonsteps = 0`, the two starts reach different points -- objective 12849.030601 vs 12838.508102 (**10.5 nats**), M1 species 1 0.38983395 vs 0.28725155 (**-26%**), SSB at `endyr` 2018 362,066 vs 330,402 mt (**-8.75%**). **Two numbers in an earlier version of this row were wrong and were quoted onward:** it said terminal SSB **+16.0%**, which read `ssb[, ncol(ssb)]` -- the year 2050 projection column, never optimised under `estimateMode = 1`, and the opposite sign from the `endyr` value; and it said both runs passed the max-gradient criterion, when `max\|gradient\|` is 2.2e-03 and 1.7e-03 against an `OK` tier of 1e-03, so both are `WARN`. The real consequence is that which optimum an unphased fit reaches is arbitrary -- at least four local minima over this range -- not any single pair of numbers. Under `phase = TRUE` the two land bit-identically. `M1_model = 2` moved by 5.7e-06. The live GOA multispecies assessment is in the affected configuration (`nsex` `c(1, 2, 1)`, `M1_model = c(1, 2, 1)`) for species 1 and 3, and **its headline multispecies fit is `phase = FALSE`** (`R/02_fit_models.R:76`), so phasing does not settle it there. Whether ITS fit moves is unmeasured; Grant declined a refit while this row still said phasing made the point moot, so revisiting it is his call. |
| `R/6-fit_mod.R` (no `log_M1` finiteness check on `start_par`) | **Open** (found reviewing PR #197, 2026-10-05). Any warm start -- `inits = <a saved fit>`, `retrospective()`, `self_test()`, `refit_like()`, `run_mse()` -- whose `log_M1` holds a non-finite REAL cell. | **The refusal 5.50.0 adds does not cover the path a warm start takes.** `.rce_stop_if_nonfinite_M1()` is called in `build_params()` and inside the `if (updateM1)` block, but **not** on `start_par$log_M1` immediately before `MakeADFun()` -- and on the `inits` branch `start_par` comes from `inits`, so `build_params()`'s array is used only as a name/order skeleton and thrown away. Measured: `tests/testthat/fixtures/golden-reference.rds` holds non-finite `log_M1` cells -- 18 (`ss`), 18 (`ms`), 11 (`goa_ss`), 11 (`goa_ms`) -- and `test-golden-regression.R` block 1 passes them straight into the fit. Harmless there, and the fixture proves why: a non-finite constant in a **mapped-off** `PARAMETER_ARRAY` is inert (block 1 asserts the objective to 1e-10 relative and `max\|gradient\| < 1e-8` at exactly those parameters, which a `NaN` on the tape could not pass). **The fix is cheap and was measured**: all 58 of those fixture cells are PADDING, zero are real, so a check restricted to real `(species, sex, age)` cells passes golden unchanged and closes the path. Left out of 5.50.0 deliberately -- a whole-array check breaks golden, and a real-cells check on its own would refuse a pre-5.50.0 saved fit of a short-`M1_base` workbook whose mapped-off NaN is provably inert. **Decide which**: real-and-estimated cells only (needs the map, so it must sit after section 7), or real cells with a warning rather than a refusal. |
| `R/3-build_map.R` (`M1_rho[sp, , k]`, `M1_dev_log_sd[sp, ]`, `M1_beta[sp, 2, ]`) | **Open, latent** (found reviewing PR #199, 2026-10-05, by two reviewers independently). A sex-specific starting value for any random-effect M1 hyperparameter. | **The same bare-sex-comma pattern 5.51.0 removed from `log_M1` is still in the random-effect half of `build_map_m1()`**, at five writes. It is inert *today* and only by coincidence: `build_params()` starts `log_M1_dev`, `M1_rho` and `M1_dev_log_sd` at a uniform 0 (`R/2-build_params.R`), so `TMB:::updateMap()`'s `tapply(..., mean)` over a level holding a padding cell returns exactly the real value. The moment any of those starts becomes sex-specific -- a build default, a hand-edited `inits`, a linkage supplying the level -- the dilution 5.51.0 fixed comes back silently, in the random-effect block instead. **Fix:** write `[sp, seq_len(nsex_sp), ...]` there too, which is free while the starts are uniform, OR leave it and state in a comment that the bare comma is safe only because the start is uniform. Do not leave it unstated. |
| `R/3-build_map.R` (`seq_len(nsex_sp)`) | **Open, low, unreached** (found reviewing PR #199, 2026-10-05). A `data_list$nsex` that is a **factor**. | 5.51.0 makes the M1 map depend on `nsex`'s TYPE, where `dev`'s bare `[sp, , ]` did not. Measured: `seq_len(2L)` and `seq_len("2")` both give `1 2`, but `seq_len(factor("2", levels = "2"))` gives **`1`** -- so a factor `nsex` of 2 would silently leave a two-sex species' MALE M1 fixed at its start. Rated unreached rather than fixed: `nsex` is `type = "integer"` in the schema (`R/0-column_schema.R:162`), it is `numeric` on all five bundled datasets checked both raw and after `switch_check()`, and `readxl` does not produce factors. **A local coercion would be false comfort** -- `as.integer(factor("2"))` is also 1, and if `nsex` were ever a factor the `nsex_sp == 2` comparisons throughout this file would misbehave too. So the fix, if one is wanted, is to assert the type once at the boundary (`switch_check()` or `data_check()`), not to coerce at each use. |
| ~~`R/3-build_map.R` (`build_map_growth`, `(sp - 1) * 4 + 5:7`)~~ | **Resolved** (found sweeping for the M1 padding-dilution class, 2026-10-04). A multispecies model where species *sp* has `nsex == 2` and `growth_model > 0`, AND species *sp+1* also has `growth_model > 0`. Not reachable by any current configuration -- see below. | **SPECIES *sp*'s MALE GROWTH PARAMETERS ARE SPECIES *sp+1*'s FEMALE GROWTH PARAMETERS.** The map stride is 4 per species but a two-sex species consumes 8 slots: females get `(sp-1)*4 + 1:3`, males `(sp-1)*4 + 5:7`, so sp1 males are `5,6,7` and sp2 females are also `5,6,7`. Verified identical by arithmetic and by inspecting the built map. Richards (`1:4` / `5:8`) has it too. This is **structural, not a start error**: K, L1 and L-infinity are fused into one TMB parameter across a species AND sex boundary, so the model cannot estimate them separately and 3 (vonB) or 4 (Richards) free parameters vanish with no message. The diluted start is only a side effect -- on an nlengths-equalised GOA2018SS, arrowtooth MALE L1 starts at `sqrt(10*5) = 7.0711` against an intended 10 cm (-29%) and cod L1 at the same value against 5 cm (+41%); L-infinity at `sqrt(80*120) = 97.98` against 80 and 120. **Not reached today:** no bundled dataset sets `growth_model` (schema default 0), all four golden models run 0, every multispecies sibling uses `build_growth(fun = "empirical")`, and every parametric-growth consumer (GOA cod, the SS3 bridge, deprecated hake) is single-species, which is immune since there is no species *sp+1*. `make_msm_test_data()` forces `nsex = 1`, so the suite cannot hit it either. **It arms the moment parametric growth is enabled on the GOA multispecies model**, where arrowtooth is two-sex at species 2 of 3. That assessment runs `build_growth(fun = "empirical")` today and its "Weight at age" section says the temperature-dependent von Bertalanffy formulation is **not** used there -- it is named as a candidate for the next cycle and as a data need, so the defect is latent rather than pending. (An earlier draft of this row called it "documented as the intent", which overstates it, and cited a numbered section the document does not have -- the Rmd sets `number_sections: false`.) A trailing two-sex species is harmless (`c(1,1,0)`, `c(1,0,1)` verified clean). **Fix:** a running counter, as `build_map_m1()` uses, rather than a fixed stride. **RESOLVED on `fix/growth-map-stride-and-caal` (5.52.0)**, with a running counter. Pinned at two levels: `test-growth-map-stride.R` drives `build_map()` over `growth_model` 1-2 x `nsex` `c(1,2,1)` / `c(2,2,2)` / `c(2,1,1)`, and one block fits `make_msm_test_data(nspp = 2)` with `nsex = c(2, 1)` through `fit_mod(estimateMode = "DebugBuild")` and pins `length(obj$par)` at **135** -- the old stride gave **132**, with 6 distinct growth levels instead of 9. The objective at the STARTING values is identical either way on that fixture, because every species starts at the same K/L1/Linf and the mean-collapse is then a no-op; that is the mechanism of the silence, and why a start-value comparison is not a test here. |
| ~~`R/2-build_params.R` (`pivot_wider(names_from = Bin)` in the growth block)~~ | **Resolved** (same sweep). Parametric growth on a multispecies model whose species have DIFFERENT CAAL bin counts. | `build_params()` hard-errors with "number of items to replace is not a multiple of replacement length": the pivot keys columns on the per-species bin ordinal, so species with different `nlengths` produce a wider matrix than the `log_growth_pars[, 1, 2:3]` target. GOA2018SS ships `nlengths = c(7, 26, 117)`, so this error currently **blocks** the growth stride collision above -- the reviewer had to equalise `nlengths` to reach it. Loud rather than silent, but fix the two together: curing this one unblocks the other. **RESOLVED on the same branch (5.52.0)** with a per-species `min`/`max` over `caal_data`, which is what the arrange-and-slice was approximating. Note the old code did NOT merely mis-size a one-bin species: `slice(c(1, n()))` duplicated the single row, `pivot_wider` warned and returned a list-column, and `log()` then errored with "non-numeric-alike variable(s) in data frame". So both the ragged and the one-bin case were hard errors. The new code warns instead, naming the species, because L1 == L-infinity leaves the growth K with exactly zero gradient at the starting values -- see the one-bin row below. |
| `R/3-build_map.R` (`build_map_m1`, the `* 2. Random Effects` block) | **Open** (found reviewing PR #203, 2026-10-05). `M1_re > 0` together with `M1_model` 3, 4 or 5 -- and `M1_re` 2, 3, 5 or 6 together with `M1_model` 2. | **The requested time-varying M is not estimated, and you get free standard-deviation and correlation parameters with nothing attached to them.** Every random-effect arm is nested inside an `M1_model` test: `M1_re %in% c(1, 4)` is handled for `M1_model` 1 and 2 (`:390`, `:402`), and `M1_re %in% c(2, 5)` and `%in% c(3, 6)` for `M1_model` 1 only (`:425`, `:464`). Nothing else matches, so `log_M1_dev` stays entirely `NA`. Measured on `GOA2018SS` (`nspp = 3`), free-parameter counts from the built map: at `M1_model` 1 `log_M1_dev` frees 4736 (`M1_re` 1), 2688 (2), 2688 (3); at `M1_model` **3, 4 and 5** it frees **0 for every `M1_re` 1-6** -- while `M1_dev_log_sd` frees **6** in all of them, and `M1_rho` frees 6 under AR1 (`M1_re` 4/5) and 12 under 2D-AR1 (`M1_re` 6). At `M1_re = 0` all three are correctly 0, so the stray sd is driven by `M1_re` alone and does not depend on the deviations existing. **`M1_model = 4` raises no condition at all**; 3 and 5 warn, but about being sex-specific on a single-sex species, which is a different thing. So a user asking for time-varying M under an environmentally-driven or age-specific M gets constant M plus up to 18 unidentified parameters in the Laplace approximation. **Fix:** refuse the combination in `switch_check()` or `data_check()` with a message naming both switches, OR implement the missing arms -- decide which, because the deviations are meaningful under `M1_model` 3 (age-specific) and the refusal would be a user-visible break. Either way map the sd and rho off when the deviations are mapped off. |
| `R/9-retro_and_jitter.R:249` (`map <- object$map`) | **Open** (found reviewing PR #203, 2026-10-05). `retrospective()` on a fit saved before any map-construction fix. | **A map fix does not reach a warm start, and the two passes of one peel can disagree.** `retrospective()` takes the peel's map from the stored object (`:249`) and hands it to `fit_mod(map = map)` (`:288`); `fit_mod()` builds a map only `if (is.null(map))` (`:850`). So a peel refit of a model fitted on older code inherits that code's map -- including the fused growth block this release fixed, and the diluted `log_M1` levels 5.51.0 fixed. The later report-only pass in the same function **rebuilds** it with `build_map()` (`:321`), so within one `retrospective()` call the peel is fitted under the old map and reported under the new one, with different free-parameter counts. Not specific to growth: it is the general shape of every map fix. **Decide the contract**: either rebuild the map for a peel (and accept that a peel of an old fit stops reproducing), or keep inheriting it and say so in `?retrospective`, and make the second pass inherit too. `/golden-check` does not cover any refit path. |
| ~~`R/2-build_params.R` (one-bin CAAL)~~ | ~~Parametric growth on a species with exactly one distinct CAAL length~~ | **RESOLVED on `fix/growth-map-stride-and-caal` (5.52.0)** by a refusal. Found reviewing PR #203, 2026-10-05, and the first diagnosis was wrong: I filed it as the growth `K` having zero gradient, which a warning would be adequate for. It is worse. `length_sd_at_age()` (`src/TMB/growth.hpp:26`) interpolates the length-at-age SD as `sd0 + (sd1 - sd0) / (linf - l1) * (len - l1)`, and `build_params()` starts both `growth_log_sd` at 0, so with `l1 == linf` that is **`0/0` = NaN** for every age above `age_L1` -- measured directly -- and the NaN reaches the age-length key and the likelihood. `data_check()` (`R/1-data_check.R`, the `length(sp_lengths) != nlengths[sp]` rule) already refuses one CAAL length unless `nlengths[sp]` is also 1, so a species with a single length bin is the only shape that got here; measured, that is exactly what it refuses. Now refused in `build_params()` with a message naming the species and pointing at `build_growth(fun = "empirical")`. |
| `build_map()` is exported but not callable on a `switch_check()`'d data_list | **Open, low** (found 2026-09, re-confirmed reviewing PR #203, 2026-10-05). Any direct `build_map(data_list, params)` call outside `fit_mod()`. | `fit_mod()` populates several switches before it calls `build_map()`, so a caller who does not repeat that work gets a zero-length switch and an error whose text names neither the switch nor the fix. Two confirmed: `suitMode` is a scalar on a bundled data_list and `fit_mod()` extends it to `nspp`; `growth_model` is **NULL** after `switch_check()` and `fit_mod()` sets it from `growthFun`, so `build_map_growth()` dies on `if (growth_model %in% c(1, "vonBertalanffy"))` with "argument is of length zero" (identical before and after the 5.52.0 stride fix -- not a regression, the same `if` on the same zero-length value). `test-growth-map-stride.R` works around both by setting them by hand, which is a note that the API is wrong, not that the test is. **Fix:** either have `build_map()` apply the same defaults `fit_mod()` does, or give it a one-line refusal naming the missing switch. Do not quietly default them -- a wrong `suitMode` is a wrong predation model. |
| `src/TMB/ceattle.cpp` (`Fixme: denominator is zero somewhere`) | **Open, and unreproduced** (found sweeping bug-history comments, 2026-10-04). Unknown -- the note states no condition. | The marker sits inside the section 16 CHANGE LOG, stranded between items 19 and 20, and reads in full: "denominator is zero somewhere. Log of negative number. Check suitability. Make other prey a very large number. Look at M2_at_age: suitability: and consumption. Make sure positive." It dates from the 2017->2018 ADMB conversion and change-log item 23 ("Fixed suitability estimation (use hindcast only)") may well be its resolution -- but nothing says so, and a log of a negative number in suitability would be a silently-wrong M2 rather than a crash. **Triage before deleting it:** decide whether it is live, and if it is, give it a reproducing input. Not removed with the other bug-history comments precisely because it may not be history. **Note the marker is lowercase `Fixme:`, so the re-derive command at the top of this file -- which is case-sensitive -- does not count it; `grep -i` gives 49 rather than 47.** |
| `plot_index()` `suffix = "survey_indices"` (`R/7-plot_diagnostics.R`) | **Open** (found covering the save path, 2026-10-08). `plot_index(file = <stem>)` and `plot_index(file = <stem>, log = TRUE)`, or `plot_logindex(file = <stem>)`, in the same script. | **Both scales write the same PNG, so only the last one survives.** The suffix does not depend on `log`, so the natural-scale and log-scale survey-index figures are both `<stem>_survey_indices.png`. `plot_logindex()` forwards with `log = TRUE`, and `GOA-ATF-ESP/R/Run_2025_ceattle.R:274` calls it with `file =`. Not a crash and not a wrong number -- a missing figure, or the wrong one embedded in a SAFE chapter. Pinned as current behaviour in `test-plot-save-paths.R`; the fix needs names chosen deliberately, since changing a written filename is user-visible. |
| `plot_diet_comp()` / `plot_diet_comp1()` filename construction (`R/7-plot_comp.R`, `paste("Pred-", species[pred])`) | **Open, low** (same sweep). Any `plot_diet_comp(file = <stem>)` where a species name contains a space. | **The written filenames contain spaces and a stray `"Pred- "`.** Measured on `BS2017SS`: `<stem>_aggregated_diet_comps_year1_Pred- Arrowtooth flounder_prey_Cod.png`. Spaces break a shell glob and a LaTeX `\includegraphics` path. `plot_timeseries()` already sanitises a species name for its CSV with `gsub("[^A-Za-z0-9]+", "_", spnames[i])` (`R/7-plot_ceattle.R`); the diet plotters do not use it. Pinned as current behaviour; the fix renames a written file, so it is user-visible. |
| `.stop_if_shared_block()` (`R/0-linkage_table.R`) / `adjust_map_shared_params()` (`R/3-build_map.R`) | **Open** (reported by the GOA cod bridge session 2026-10-08, mechanism verified here the same day). A selectivity or catchability linkage stratified by fleet -- which is the DEFAULT for `sel` and `q`, so any spec without an explicit `by = NULL` on two or more fleets that share a `Selectivity_index` or `Catchability_index`. | **The base block is shared and the linkage offsets are not, so fleets declared to mirror each other get different realised selectivity with no warning.** `adjust_map_shared_params()` reconciles 15 parameter slices across the group, every one indexed BY FLEET -- including `sel_dn6`, `sel_inf`, `log_sel_slp`, `sel_coff`, `log_sel_apical` and `index_q_beta`. The linkage coefficients are not among them: they live in `beta_linkage`, a flat `PARAMETER_VECTOR` (`ceattle.cpp`) with no fleet dimension, and the linkage table keys its rows on `fleet` ("1-based Fleet_code; NA = shared"), so one design column on a shared group becomes one FREE coefficient per member fleet. Nothing in `R/` reconciles `beta_linkage` across a group. `.stop_if_shared_block()` refuses a linkage that NAMES a single follower, but returns early on `length(flt) != 1L`, with the comment "an unstratified linkage expands to every fleet, where setting the whole group is unambiguous and the followers take the donor's value anyway" -- that last clause is the false premise. **Measured by the reporter on GOA cod (5.52.0.9001 and identically on afsc/dev 5.53.0), Srv (fleet 4) and Srv_ae1 (10) sharing `Selectivity_index` = 4:** 47 distinct selectivity block columns become 59 linkage rows, the 12 `s4*` columns appearing on both fleets; estimated freely the copies diverge on the log scale, `s4p5_blk1996` -3.82671 against 17.05918 (+20.89), `s4p1_blk2006` 9.13890 against -3.11320, largest absolute difference **20.89**. In that model it also explains 12 surplus free parameters, `check_estimability` flagging `index_log_q` non-identifiable on Srv_ae1, `sdreport` failing on a non-invertible Hessian, and the fit reaching ~33 nats below SS3's optimum with FEWER declared parameters. **Applies to q linkages too**, by the same table and the same flat vector -- the reporter did not test q; the legacy per-fleet `index_q_beta` route is safe because it IS copied. **Interim guidance -- and note this is the DEFAULT path, not a rare spelling:** `by` defaults to `~ fleet` for both `sel` and `q` (`.default_stratum()`), so ANY selectivity or catchability linkage written without an explicit `by` is stratified by fleet and emits one free coefficient per fleet. `fleet =` is a FILTER, not the stratification -- its roxygen says `NULL` (the default) means "every fleet in `strata$fleet`" and that "`by` must include `fleet` for it to apply" -- so OMITTING `fleet =` makes it worse, not better: the linkage then applies to every fleet in the model. The shared-block spelling is an explicit **`by = NULL`**, which `expand_linkage_strata()` turns into a single row with `fleet = NA` and one coefficient, and which `linkage_spec()` keeps as given. So: pass `by = NULL` for a parameter the group shares, or do not share the index. **The fix is a design choice and not an agent's:** reconcile `beta_linkage` on the sharing index the way the 15 slices are reconciled, or extend `.stop_if_shared_block()`'s refusal to the stratified case. Full write-up with the table in Rceattle-models `GOA cod/Bridging/Estimation_Differences.md` entry #21 (`39808bd`). |

---

**All nine are resolved as of 5.13.0.** The table is kept struck through rather than deleted:
each row records what the defect actually was once reproduced, which in six of the nine differed
from what its FIXME claimed. Add new rows above it.

| Where | Condition | Consequence |
|---|---|---|
| ~~`src/TMB/ceattle.cpp` (`will bomb if minage > 1`)~~ | ~~`minage > 1`~~ | **Resolved in 5.13.0**: it did not bomb, it read adjacent memory. `ssb(sp, yr - minage(sp))` went negative for the first `minage - 1` years and Eigen does not bounds-check in a release build. Measured (BevertonHolt, nages 5): R came back `8.6e-314` with `is.finite()` TRUE, objective 14407.38 → 15162.60 → 17532.15 across minage 1/2/3. Those years now take `R_init * exp(rec_dev)` -- equilibrium recruitment at `F = Finit`, already what year 0 uses -- following Stock Synthesis's equilibrium-plus-early-devs treatment of the pre-start period (WHAM fixes the lag at one year instead). NOT `R0`: under a stock-recruit hindcast `build_map()` maps the mean-recruit parameter out and only `R0[, 1]` is overwritten with the derived `(alpha - 1/SPR0)/Beta`, so a guard reading `R0[, yr]` gets `exp(9) = 8103.08`. The first fix did exactly that, giving 14407.38 &rarr; 20421.68 &rarr; 30920.21; corrected it is 14407.38 &rarr; 13959.97 &rarr; 13587.68. Expected recruitment (`R_hat`) shares the anchor, or the stock-recruit penalty reads the mean-versus-curve level gap as signal (0.83 nats/year under the Ianelli configuration). All four sites guarded; `minage = 1` cannot fire, golden unmoved. `test-dynamics-recruitment-minage.R`. |
| ~~`src/TMB/ceattle.cpp` (`will blow up if nlengths is less than nages`)~~ | ~~`nlengths < nages`~~ | **Resolved in 5.13.0**: `age_hat`/`age_obs_hat` are written at AGE indices (to `nages*2` joint-sex) but were sized `= comp_obs`, whose width is the workbook's `Comp_` columns — `nlengths` for a length-only model. Silent out-of-bounds WRITE, not a crash. Both now sized from the widest age index the model's own `comp_ctl` says will be written -- `nages`, or `nages*2` for a joint-sex row -- and never narrower than the observations. Reserving `max_age*2` unconditionally (the first fix) also closes the overrun but widens every model's REPORTed `age_hat`: BS2017SS went 25 &rarr; 42 columns, 17 permanently zero, and assessment scripts `cbind` that array. No end-to-end reproduction is possible from R (an unchecked write is not observable), so `test-composition-age-hat-width.R` pins the invariant and the absence of `= comp_obs`. |
| ~~`R/10-run_mse.R` (`does not work for assessments that don't occur annually`)~~ | ~~assessment interval ≠ 1 year~~ | **Resolved in 5.13.0**: only the scalar-`cap` branch was affected. `dat_fill_ind` spans the whole interval, so `sum(Catch[dat_fill_ind]) > cap` held a multi-year total to a one-year ceiling — at `assessment_period = 2`, roughly halving projected catch. Now applied per projection year. The same line carried a second and larger defect: `ifelse()` returns the shape of its length-1 test, so whichever branch was taken was truncated to its first element and recycled across every row. Two species at 80/20 t came back 40/40 against a 50 t cap (80 t total, over the ceiling) and 80/80 against a non-binding 500 t cap (160 t total, catch invented). **Any stored MSE using a scalar `cap` moves, at every `assessment_period` including 1, binding or not** -- only an exactly equal split was safe. The species-specific vector branch was always per row. `test-mse-cap-and-hcr2-threshold.R`. |
| ~~`R/10-mse_summary.R` (`will bug if not survey`)~~ | ~~a species whose fleets are not all surveys~~ | **Resolved in 5.13.0**: not a defect. `spp_rows <- which(flt_spp == sp)` was assigned once and read by nothing, in this file or anywhere else in the package. The FIXME speculated about code that did nothing; both lines are gone. No behaviour change. |
| ~~`R/3-build_map.R` (`QAR1 is inert`)~~ | ~~**`Catchability = "AR1"`** (the QAR1 form, Rogers et al. 2024)~~ | **Resolved in 5.12.0**: `data_check()` now errors on it, so the branch is unreachable. It was inert — the deviate map is gated on `Time_varying_q %in% c("IID","AR1","RandomWalk")`, but under `Catchability = "AR1"` that column holds an `env_data` **column index**, not a mode, so `index_q_dev` stayed mapped out and q was constant. Not repaired: the Rogers form is implemented correctly by a q linkage (`ar1(1 \| Year)` with `observe`), which GOA pollock 2025 runs. The dead `build_map()` branch is deleted in 5.13.0. Code 6 **stays in `q_map`**: `validate_switches()` runs before `data_check()`, so dropping it would replace that migration message with a generic "invalid value" for exactly the workbooks that need the recipe (GOA pollock 2024/2025 still carry a 6). These are two different switches sharing a string; an earlier draft of this file named the wrong one. |
| ~~`src/TMB/ceattle.cpp` (`caal_ll_type`)~~ | ~~`CAAL_distribution = "MultinomialAFSC"`~~ | **Resolved in 5.13.0**: implemented, as the catch-sigma case of this shape was in 5.12.0. The AFSC multinomial pseudo-likelihood is a published AMAK form already implemented for age comps, so extending it to CAAL was mechanical rather than inventing a likelihood. Verified against the form computed by hand from the reported CAAL proportions (2.8e-14). The `test-schema-cpp-dispatch.R` exemption is removed. `test-likelihood-caal-afsc.R`. |
| ~~`R/3-build_map.R` (`will fail if random_sel = TRUE?`)~~ | ~~`random_sel = TRUE` + `Time_varying_sel = "Block"`~~ | **Resolved in 5.13.0**: confirmed real, and worse than "will fail". The block parameters live in `log_sel_slp_dev`/`sel_inf_dev`, and `fit_mod()` declared those arrays random unconditionally — but the template scores selectivity deviates only for `IID`/`AR1`/`RandomWalk`/`RandomWalkAscending`, so blocks were Laplace-integrated against **no density**, `sel_dev_log_sd` mapped out so there was no variance either. Measured: 8 parameters random, `JNLL_SEL_DEV` identically 0, objective `NaN`, real fit dead with TMB's `NA/NaN gradient evaluation`. Now refused with a message naming the fleets and the way out. `test-selectivity-random-sel-block.R` carries the reproduction and a drift guard pinning `Block` as the only mode the template leaves unscored. |
| ~~`R/6-fit_mod.R` (`swallows EVERY warning build_map() raises`)~~ | ~~any~~ | **Resolved in 5.13.0**: the comment named the wrong warnings — the shared-block ones it cited are raised by `data_check()`, not `build_map()`. What was actually swallowed was `build_map()`'s own set (M1 sex mismatch, selectivity-form incompatibilities), each of which changes what is estimated. Now de-duplicated via `withCallingHandlers()` and passed through, so `.refit_like()`'s per-peel re-entry prints each distinct warning once instead of hundreds of times. |
| ~~`R/10-mse_summary.R` (`EM uses fixed-depletion proxy for HCR 2`)~~ | ~~HCR 2~~ | **Resolved in 5.13.0**: in single-species mode the HCR 2 arm now reads `ssb_limit_thresh()`, the helper the operating model already uses, so both sides of the cross-tab score one criterion (absolute `0.5 * SBF`) and branch on the same scale flag. `Plimit` is NOT the answer -- `build_hcr()` defaults it to 0, so reading it reports a default ConstantF run as never overfished; the first fix did that. Under `msmMode > 0` both sides still fall through to `Plimit`, which is the operating model's own multispecies rule, so they agree and it is left alone. `test-mse-cap-and-hcr2-threshold.R`. |


## Tier 1 — stated limitations, currently by design

Not bugs, but they bound what the model can be asked. Worth documenting in a vignette rather
than fixing.

- **An observation the R inclusion set keeps but the template skips returns a residual of
  exactly 0 under `method = "cdf"`, where a Gaussian method returns `NaN`.** A row scored by
  neither density nor CDF leaves `nlcdf.lower == nlcdf.upper`, so `oneStepPredict()` recovers
  `F = 0.5`. The two inclusion sets are built independently — `build_osa_data()` in R, the
  `pos >= 0` and year/type conditions in `ceattle.cpp` — and agree for every fitted row today.
  A future divergence would show as a clean-looking zero rather than a visible gap. A REPORTed
  count of CDF-scored positions, asserted in R against the residual count, would make it loud.
- **Forecast growth is ignored** by the retrospective and MSE projection paths
  (`ignores forecasted growth`, twice in `R/9-retro_and_jitter.R`, once in `R/10-run_mse.R`) —
  the terminal-year growth is carried forward.
- **Projection quantities are held at the terminal hindcast year** (`R/10-run_mse.R`,
  `assuming same as terminal year of hindcast`).
- **`ration_data` is sized for the hindcast only** (`R/5-rearrange_data.R`,
  `Change for forecast`).
- **SPR reference points**: `sex_ratio` is an input rather than estimated for two-sex models,
  and the M used is the terminal-year value. **Neither has a marker** -- `rates for a reference
  point are the terminal hindcast year's` (`src/TMB/ceattle.cpp:1595`) is an ordinary comment, so
  the grep recipe above will not find it, and no `TODO`/`FIXME` mentions `sex_ratio` at all.
  5.24.1 corrected SPR to apply only the recruitment split `sex_ratio(sp, 0)` to a two-sex
  species (`female_split`, `ceattle.cpp:1605`); the ratio itself is still read from data.
- **Linkage random-effect priors are penalties, not proper densities** (`src/TMB/ceattle.cpp`,
  `FIXME(jacobian)`, twice). The sigma and rho priors sit on the natural scale without the
  Jacobian of the `log` / `rho_trans` transform. Fine as a penalty under maximum likelihood;
  under a Bayesian (`tmbstan`) run the stated density is not the prior actually applied. A
  `lognormal` sigma prior is exempt; rho has only normal and beta families.

## Tier 2 — design notes and refactor wishes

Cleared in 5.14.0 except where noted. As in Tier 0, three of these were not what their marker
said, so each struck row records what it actually turned out to be.

Found reviewing the 5.34.0-5.41.0 release (PR #158) and recorded rather than fixed. The
three defects that review found are fixed in 5.42.0; these are what it left:

- **`log_sel_apical` is unbounded.** `build_parameter_bounds()` gives it the default
  `+-Inf`, while `.check_sel_apical_rows()` only *warns* when a fleet has neither
  joint-sex composition nor a prior -- i.e. when nothing informs the sexes' ratio and the
  estimate is "whatever the optimizer leaves". That is the flat ridge `rec_pars[, 2:3]`
  got `+-30` for in 5.39.0. A `+-10` bound costs nothing: `exp(10)` is already absurd for
  a selectivity multiplier. So the package now bounds some blocks and not this one --
  the same inconsistency the `log_Ftarget` note above records.

- **Four shape-penalty columns are form-9-only and silently inert on forms 2 and 13.**
  `Sel_shape_mode`, `Sel_pen_first_bin`, `Sel_pen_last_bin` and `Sel_avgsel_pen` are read
  only inside `if(flt_sel_type(flt) == 9)` in `ceattle.cpp`, but `data_check()`
  range-validates all four on every fleet and the schema `doc` strings -- which ship
  verbatim into `meta_data_names.xlsx` -- say "non-parametric" without qualification. A
  user narrowing the shape penalty on a `NonParametricIntegrable` fleet gets the full
  range with no message. 5.40.0 added `test-selectivity-norm-scope-inert-forms.R` for
  exactly this class on `Sel_norm_scope`; these four are owed the same treatment.
  (5.42.0 fixed the fifth member of the set, `Sel_devmag_sd`, because that one was being
  actively converted rather than merely accepted.)

- **`.check_stock_recruit_msm()` reports a false clean when `spnames` is `NULL`.** It
  falls back to `paste0("Species", seq_len(...))` while `parameter_index()` falls back to
  `as.character(idx)` -- `"1"`, `"2"`. The two disagree, so `pos_of()` matches nothing for
  every species, `fixed_curve()` returns `TRUE`, and the check reports "Stock-recruit
  curve held at its inputs for Species1; nothing to check" on a curve that is in fact
  being estimated. One fallback should call the other.

- **The Dirichlet-multinomial OSA fallback is announced with `message()`.** Substituting
  a method documented as failing the package's own KS self-test on composition data is
  the most consequential of the four announcements in `osa_residuals()`, and it is the
  only one that is not a `warning()`. A `message()` is erased by `suppressMessages()`, by
  a knitr chunk with `message = FALSE` (how the vignettes and most assessment scripts
  run), and by any log that keeps only warnings.

- **`discrete = TRUE` feeds fractional composition counts into TMB's integer lattice.**
  Counts enter as `(proportion + comp_offset) * N` with `comp_offset = 1e-5`, so none is
  an integer, while TMB's `discrete = TRUE` path replaces `integrate()` with a sum over
  `ceiling(lower):floor(upper)`. For `obs = 10.001` the observation's own mass is in
  neither tail sum and the mass at 11 is dropped, yet `px` is still the density at the
  fractional `obs`. It returns a finite, plausible residual rather than erroring.
  Pre-existing; 5.41.0 re-documented and re-defaulted around it.

- **`.rce_sel_norm_code(allow_all = )` is called two ways.** `R/0-build_selectivity.R`
  passes `allow_all = TRUE` unconditionally; `R/1-data_check.R` passes
  `fc$Selectivity[rows] %in% c(11, "LogisticPM")`. The two therefore disagree on what
  counts as "normalization is on" for a non-`LogisticPM` fleet, which can flip the
  `WithinSex` apical refusal the wrong way. They should share one rule.

- **The non-parametric form set is a literal in five places** (`R/6-fit_mod.R`,
  `R/1-data_check.R` twice, `R/3-build_map.R`, `R/0-switches.R`) with no predicate on
  `sel_map`. All five are correct today, and a sixth site omits
  `NonParametricIntegrable` deliberately. But 5.40.0 -> 5.41.0 churned this set twice
  inside one release, and the `fit_mod()` copy failing silently would mean coefficients
  below `Bin_first_selected` quietly stop being held at 0 -- a selectivity change, so an
  SSB change. 5.42.0 added `.RCE_SEL_PEN_POSITIVE` beside `sel_map` as a start; the form
  sets themselves are still literals.

- **A warm start silently overrides a changed `Sel_curve_pen` column.** `sel_curve_pen` is a
  `PARAMETER_MATRIX` mapped off, not a `DATA_` object, so its value comes from `inits` when
  `inits` are supplied and from the `Sel_curve_pen1/2/3` columns only otherwise. Editing the
  column and refitting from a stored fit therefore keeps the OLD weight, with no message.
  Found building the form-9 directional test in 5.42.0, where setting the column to -20 while
  passing `inits` produced the +20 penalty. 5.42.0's `fit_mod()` guard closes the SIGN only;
  the magnitude is still silently overridden, measured on `BS2017SS` fleet 1 with the column
  reading +20 throughout: `JNLL_SEL_NONPARAM` is 52.96 from the column and 464.39 when
  `inits$sel_curve_pen[flt, 1]` is 200, with no message. This is the same class as the known
  `Time_varying_sel_sd`-inert-on-a-warm-start trap in `TRAPS.md`, and the fix is the same
  shape: either reseed the parameter from the column in `fit_mod()`, or warn when they
  disagree. Note this also bounds the blast radius of the sign defect 5.42.0 fixed -- a
  refit from a stored fit kept whatever weight was taped.

- **`Sel_shape_dir = "Increasing"` has no FITTED recovery check.** 5.42.0 closed the
  specification half: `test-selectivity-penalty-sd.R` pins the limiting cases (a strictly
  increasing curve charges 0 under `"Decreasing"` and the mirror under `"Increasing"`)
  and matches three shapes against the ADMB/AMAK `sel_like(1)` SSQ recomputed from
  `sel_coff`, driven through the `Sel_shape_dir` column itself. What is still open is
  simulation self-consistency: nobody has simulated from an increasing-selectivity stock
  and checked the penalized fit returns that shape. That gap is why 5.42.0 refused the
  direction on the forms that do not read the sign rather than teaching them the branch
  `"NonParametricPM"` (9) has.

Found during the 5.34.0-5.41.0 batch and recorded rather than fixed:

- **An estimated `log_Ftarget` is a single unbounded `nlminb` start.** The single-species
  projection starts it at its `inits` value with no bounds (`R/4-build_parameter_bounds.R`
  mentions it nowhere). 5.34.0 resets a no-fishing start (-999 or non-finite) to 0, but a large
  start also sticks: NPFMC on `make_test_data()` from `log_Ftarget = 3` (F = 20) stays at 20.08
  with objective 83886.99, against 83877.42 from 0. The realistic route in is a stored `CMSY` or
  `ConstantFSSB` fit with a large estimated `Ftarget` handed on as `inits`. Bound it, or start
  every estimated `log_Ftarget` at 0 as the multispecies loop already does. Note 5.39.0 bounded
  `rec_pars[, 2:3]`, so this file now bounds some blocks and not this one.
- **`run_mse(regenerate_past = TRUE)`'s average-F refit never runs on the string path.**
  `R/10-run_mse.R:647` tests `em$data_list$HCR == 2`, but the rule reaches it under its name
  (`"ConstantF"`) there, so the branch is dead. `.normalize_hcr()`'s own comment says either
  form may be stored depending on the processing path, so it is dead on that path rather than
  universally; `.normalize_hcr()` (`R/10-mse_summary.R`) is the existing way to compare either
  spelling. If it is revived, `Ftarget` needs a full per-species vector: `avg_F$avg_F` covers
  only species with a fleet in `fleet_control`, and `extend_length()` stops on any other length.
- **Three tests read as guards and never run.** Each calls `testthat::skip()` unconditionally,
  so a full `NOT_CRAN=true` suite reports them as skips among 9,506 passing assertions and
  nobody notices. Found running the release suite 2026-09-21.
  - `test-dynamics-multi-spp-model.R:119`, "Equilibrium MSVPA suitability dynamics match" --
    the only one with a stated reason, inline: a minor unexplained difference, possibly bias in
    diet weighting, where the old EBS CEATTLE still matches. **That is an open numerical
    discrepancy against the reference implementation and it is recorded nowhere else.**
  - `test-dynamics-fit-sanity-model.R:5`, "key quantities match baseline" -- no reason given.
  - `test-data-input-validation.R:4` -- no reason given, and it also carries `skip_on_cran()`,
    so it is doubly inert.
  Either restore them or delete them; a skipped test that names a baseline is worse than no
  test, because the suite reports a guard that is not guarding.

- **`run_mse()` carries every deviation array into the projection except `log_M1_dev`.** The
  carry is commented out at `R/10-run_mse.R:901` for the operating model, under the
  `#FIXME - simulate` marker, and again at `:1126` for the estimation model's refit, so a model
  with `M1_re` projects at zero M deviation. Carry the terminal year, as `index_q_dev` is.
- **`fit_mod(initMode =)` overwrites `data_list$initMode` unconditionally** (the argument
  defaults to `"NonEquilibrium"` at `R/6-fit_mod.R:185` and is written to the data list at
  `:424`), so the data object's own `initMode` field is never read. `BS2017MS$initMode = 1` has
  therefore never reached the golden `ms` fit. A `model_config` slot on the data IS read, since
  5.36.0 (`:385`) -- it is the bare field that is not. Also in `TRAPS.md`.
  **`random_sel` has the same shape** and is worth knowing before anyone attempts the "one
  object carries every switch" simplification: `fit_mod()` writes `data_list$random_sel` from
  its own argument at `R/6-fit_mod.R:420`, so setting the field on the data object is silently
  ignored. A `config` slot is read first (`:369`), as for `initMode`.
- **The AMAK selectivity start is not the package's** (`src/TMB/ceattle.cpp:4073`,
  `FIXME: AMAK starts at nbins/2`). A formulation divergence, not a defect; record it where a
  bridging exercise will find it.
- **Two places name selectivity codes that no form uses.** The shared normalizer's gate is
  `sel_type != 5 && sel_type != 12 && sel_type != 11` (`src/TMB/selectivity.hpp:74`) though
  `sel_map` has no 12, and the invalid-`Selectivity` error offers `range(sel_map)` as the allowed
  set (`R/0-switches.R:1326`), i.e. "0:14", which includes the unused 10 and 12. Neither is a
  defect -- `switch_check()` refuses a fleet carrying 10 or 12 -- but the gate reads as though
  three forms normalize per sex when two do, and the error names codes it would reject.

- ~~**Split `R/0-build_srr_and_M.R`**~~ — **Done in 5.14.0**, but not as described. The file was
  1,497 lines and **52** top-level objects, not 29, and the three-way srr/M1/growth split named
  here would have stranded 612 lines (41%): catchability, selectivity and composition linkage
  machinery that touches none of those three. Split six ways instead, one file per process, all
  keeping the `0-` prefix. `.coerce_switch_arg` went to `0-switches.R` beside `.canon_switch()`
  (it is generic over `srr_fun`/`M1_model`/`sd_plus_group`); `.default_stratum`,
  `.resolve_auto_by` and `.stamp_param` went to `0-build_linkage.R`, which already called the
  last of them. The roxygen block for `.check_q_linkage_support` was bound to
  `.message_auto_fleet_linkages` and is reattached. Pure relocation: 441 object bodies unchanged,
  and the multiset of non-blank lines across `R/` is identical.
- ~~`R/5-rearrange_data.R` (empty comp/CAAL frames, `age_error`)~~ — **Done in 5.14.0.** The comp
  and CAAL blocks were the same five lines twice; now one `.normalise_rows()`. The zero-row guard
  stays — `t(apply())` does not preserve the shape of a matrix with no rows. `age_error` keeps its
  `as.data.frame()` coercion, which is load-bearing because the loop mixes `$` and positional
  `[i, ]` access, and its `1:nrow()` becomes `seq_len()`.
- ~~`R/2-build_params.R` (`variance and AR1 parameters`)~~ — **Not work.** Both already exist
  directly above the marker (`log_sigma_linkage`, `trans_rho_linkage`), and the `else` branch
  zeroes every name the `if` branch sets. A leftover placeholder at a section boundary; deleted.
- ~~`src/TMB/ceattle.cpp` (`can probably outside iter loop`)~~ — **Not work.** Section 5.12 is
  already ~260 lines above the only `iter` loop. The comment now states what the placement means,
  including that the forecast arm is provisional and recomputed in section 6.7.
- ~~`R/0-clean_data.R` (`may be redundant now?`)~~ — **Half right, and not the half that matters.**
  The template does derive `SB0`/`B0` itself, and reads `MSSB0`/`MSB0` only to overwrite that under
  `msmMode > 0`. But neither has `read_data()`/`write_data()` support, so no workbook can supply
  them and this default is the only thing that creates the required `DATA_VECTOR`s. Kept, with the
  reason and with `999` named as the placeholder `fit_mod()` fills in. **Correction, 5.15.0:**
  section 10.2 filled it into `data_list_reorganized` only, so the *returned* `data_list` kept the
  999 and every refit off a fitted object re-entered the template with it. Fixed by carrying
  `MSSB0`/`MSB0` onto `mod_objects$data_list`.
- ~~`R/3-build_map.R` (`add checks for surveys sel sigma`)~~ — **Done in 5.14.0.** Real: fleets
  sharing a `Selectivity_index` estimate one `sel_dev_log_sd` between them, and a differing
  `Time_varying_sel_sd` was reconciled silently. Not to the first member's value, which is the
  intuition to unlearn: TMB's `updateMap()` collapses a shared parameter with
  `tapply(par, map, mean)`, and this one is held on the log scale, so the group starts at the
  **geometric mean** of its estimated members' values — `sqrt(0.3 * 0.7) = 0.4583` for a
  two-fleet group at 0.3 and 0.7, a value neither row asks for. `.warn_shared_block_start()` reports
  it, once per group, and runs at the END of `build_map()`: `build_map_f_and_data_weights()`
  maps the parameter out for `Off` fleets and `build_map_fixed_natage()` for a fixed-dynamics
  species, both after the sharing pass, so a check placed inside
  `adjust_map_shared_params()` counts fleets that end up estimating nothing.
- ~~`R/3-build_map.R` (`add checks for surveys q sigma`)~~ — **Done in 5.14.0, after fixing the
  reason there was nothing to check.** `index_q_dev_log_sd` was mapped out for every fleet and
  never turned back on, so no shared group could discard it. That was itself the defect:
  `random_sel` frees `sel_dev_log_sd` alongside the selectivity deviates it integrates out, but
  `random_q` integrated the catchability deviates out and left their sd fixed at
  `Time_varying_q_sd`. Now symmetric, so `random_q = TRUE` estimates it — **and any fit using that
  flag moves.** With the sd estimable the shared-group copy is meaningful, so a shared
  `Catchability_index` goes through the same `.warn_shared_block_start()` check, on the same
  geometric-mean footing described in the row above.

  Caveat, measured rather than assumed: on a 40-year index with the observation sd FIXED and q
  deviations injected at a true sd of 0.4, the marginal MLE pins to its lower bound at observation
  sd 0.1, 0.3 and 1.0, and reaches only 0.06 at 0.05. A short or noisy index can therefore return
  an sd that reads as a constant q. That is a diagnostic to check, not a reason to withhold the
  parameter -- the estimate and its gradient show it, and whether the series informs it is the
  assessor's call. `index_q_log_sd`, the prior sd on q itself, stays fixed: estimating the width of
  one's own prior is not meaningful.

Still open. No user-visible consequence; do them opportunistically.

- `R/10-run_mse.R` (`extract run_one_sim as a top-level internal helper`) — 455 lines, 18 free
  variables. Three hazards, all verified: `%!in%` is defined **inside `run_mse()`** and is not a
  package-level object; the `<<-` in the OM-no-F handler currently resolves to the closure's own
  `sim_list` and would walk to the namespace if the handler moves out; and `estimate_mode_base`
  and `sim_dat` exist as both `run_mse()` locals and closure-local re-assignments, so the closure
  does not depend on the outer copies. Note also that the closure's environment is passed to
  `.parallel_lapply()` for the PSOCK export, so the free-variable set is load-bearing for
  parallelism, and that `test-mse-cap-and-hcr2-threshold.R` greps this file for literal source
  strings and must move in lockstep. Needs `verify-mse-repro.R` and `verify-mse-om-horizon.R` as
  before/after digests.
- `R/3-build_map.R` (`use formula`) — retiring the `Time_varying_q` overload that holds
  comma-separated `env_data` column indices in favour of the q linkage. Deprecates a public switch
  path, so it needs a shim, not an opportunistic edit.
- `R/3-build_map.R` (`may want sex-varying?? Hard to estimate`) — `M1_rho` already has a sex
  dimension, but the mode-6 branch offsets by `nspp` to avoid colliding with the mode-4/5
  `sp`-valued rho, and going sex-varying needs that counter scheme reworked. The markers say the
  obstacle is estimability; decide that first.
- `src/TMB/ceattle.cpp` (`fit the window form directly when a fleet needs it`) — needs `t1` and
  `D` per fleet plumbed through to the template, **plus seasonal dynamics that do not exist**: the
  annual recursion spreads F evenly across the year, so the window predictor would be inconsistent
  with the dynamics fitted against it. The derivation above the marker already quantifies the
  current approximation at ~1.5% trend error over F 0.05–0.8, against −29% to +33% for the
  snapshot it replaced.
- The `logH_*` and `log_gam_*` markers belong to the stubbed Kinzey-Punt predation forms
  (`H_4` is NOT one: it is a plain declaration comment inside the commented-out block)
  (`msmMode` 3–9) and the gamma predator selectivity. They are pinned as stubbed in
  `tests/testthat/test-schema-registries.R`; leave them until that work is picked up.
- `src/TMB/ceattle.cpp` (`penalize every selectivity deviation rather than a sub-range`) —
  would pin the unidentified directions and drop the year/bin indexing of the deviation penalty.
  It would **not** retire the four columns the marker names: `Sel_cap_bin` holds the
  NonParametricRPM curve flat past a bin, `Sel_start_year` builds the curve from the base
  coefficients through that year and pins the random walk's level in `build_map()`, and
  `Sel_pen_first_bin` / `Sel_pen_last_bin` bound the shape penalty, not the deviation penalty.
  Moves every fit with penalized deviations, so it needs `/golden-check`.
- ~~`R/0-osa_data.R:80` — the comment names only `switch_check()` as what fills
  `comp_offset`.~~ **Resolved in 5.36.0**: it names all three fill sites (`:79`-`:82`). Never a
  marker, so it never appeared in the counts above; the similar-sounding
  `switch_check() does not run` sits at `R/5-rearrange_data.R:202` and is about `Sel_norm_bin`.

- ~~**Single-species hindcast-curve projection double-counts the SSB drop.**~~ **Resolved in
  5.33.0**: `sample_rec(sample_rec = FALSE)` and `retrospective()` take the hindcast's mean
  deviation from the curve, `log(mean(exp(rec_dev)))`, for every curve fitted in the hindcast.
  The old `log(mean R) - log(R0)` was off by the hindcast-average R/R0 the curve implies: 8.6%
  low for Beverton-Holt at h = 0.8 and 40% of SB0, often high for Ricker.
  `test-functions-sample-rec-curve.R`; `test-functions-retrospective.R` covers `retrospective()`
  only for a multispecies curve, and there is no single-species peel test.
- ~~**Dynamic B0 under the penalty form**~~ **Resolved in 5.33.0** (fed8cff8, from a18331b1):
  before `srr_mse_switchyr` the unfished run takes `log R - log R_hat`, not `rec_dev` (R/R0). On
  the hake operating model the curve sat at about 1.5 × R0 there, so no-fishing recruitment ran
  about 1.5 times too high; hake dynamic SB0 falls to 0.60-0.83 of its old value.
  `test-dynamics-dynamic-b0-ianelli.R` checks dynamic B0 equals the hindcast with F near zero.
- ~~**`.map_switch()` passes a factor through**, so a factor `srr_est_mode` skips
  `build_srr()`'s checks and fits as its level code.~~ **Resolved in 5.35.0**: it coerces a
  factor to character first (`R/0-switches.R:378`), with a comment naming this failure.
- **The refit warning for retired srr codes is hidden** by the `suppressWarnings()` wrapped
  around `.refit_like()` in `retrospective()`, `jitter()`, `profile()` and `self_test()`.
- **A one-year retrospective peel** averages over that year, though its warning says "after the
  first".

- **`M1_mult.sum()` is marked "LEGACY (scheduled removal: v4.5.0)" and is still live**
  (`src/TMB/ceattle.cpp`, `M1_mult.sum() is LEGACY`). One path, not two -- the marker text itself
  is the second grep hit. The sum is added unconditionally inside the species/sex/age/year loop
  just below the marker. The NEWS target (`## Scheduled removal (v4.5.0)`, in the 4.1.0 section)
  still exists, so the pointer resolves -- but a reader at 5.49.5 will reasonably conclude this
  path is gone. Either retire it behind a deprecation message or restate the comment as current
  behaviour with the condition under which it executes. Found sweeping bug-history comments,
  2026-10-04.
- **A second, unrelated `LEGACY (scheduled removal: v4.5.0)` marks the dynamic reference-point
  recruitment split** (`src/TMB/ceattle.cpp`, `LEGACY (scheduled removal: v4.5.0). Use
  relationship below`). Inside the Option 1a arm (`proj_mean_rec == 1 & srr_pred_fun < 2`), the
  dynamic RPs take observed `R(sp, yr)` for hindcast years and `exp(log(avg_R) + rec_dev)` for
  projection years; the marker says to use the Option 2 relationship below instead. This is NOT
  the defect resolved in 5.13.0 -- that was Option 1a and the curve arm matching **neither**
  condition. It feeds `N_at_age_dB0`, hence `SB0(sp, nyrs-1)`, which HCR 5 and 6 read as the
  depletion reference, so decide it deliberately rather than by sweep. Found 2026-10-04.

- **Three registry guards build their check-set *with* a regex and then floor it at
  `expect_gt(length(.), 0)`, so a reformat shrinks the set silently instead of going red.**
  `source-guards.yaml` made these a merge gate in 5.49.5, so the weakness now costs more than it
  did. In each case the fix is to pin the hit *count* rather than a floor, with a one-line comment
  saying what the number is, so a new or reformatted site goes red and names itself. Found by the
  adversarial review of PR #195, 2026-10-04, which demonstrated the first one.
  - `test-dynamics-sex-index-bounds.R` (`every write to sex index 1 is guarded`) builds
    `pat <- paste0("^\\s*(", paste(arrays, collapse = "|"), ")\\(sp, 1,")` -- line-anchored, with
    the exact spacing `(sp, 1,`. **Demonstrated escape:** an unguarded
    `if(true) N_at_age(sp, 1, 0, 0) = Type(0.0);` inserted in `ceattle.cpp` was invisible because
    the array name is not at line start, and the guard still passed on the remaining hits. Any
    write placed after a condition on the same line, or spelled `(sp,1,` or `(sp, 1 ,`, escapes
    it. It also looks back only 6 lines for the guarding condition.
  - `test-dynamics-recruitment-minage.R` (`no recruitment site indexes the spawning year without
    guarding it`) does `uses <- grep("yr\\s*-\\s*minage\\(sp\\)", src, value = TRUE)` then
    `expect_true(all(grepl(named, uses)))` -- `all()` over a set the regex defined. Two assertions
    total, and the test's own comment concedes "a new site still needs reading".
  - `test-composition-age-hat-width.R` does
    `sizing <- grep("max_age_cols", src, value = TRUE); expect_gt(length(sizing), 0)` then
    `expect_false(any(grepl("max_age", gsub("max_age_cols", "", sizing))))` -- renaming the
    variable empties the set and the assertion passes on nothing. Its third assertion, a
    `fixed = TRUE` match on a whole C++ line, is the strong one and fails closed.
- **`source-guards.yaml` makes the registry check APPEAR on a pull request but cannot make it
  REQUIRED.** That is branch protection on `afsc-assessments/Rceattle` for `main` and `dev`, and
  this repo holds no branch-protection-as-code (`find .github -type f` is the workflows and
  nothing else). Until the setting is changed, a red result there does not stop a merge. Grant's
  call; noted here so it is not lost.

- **`fit_mod(updateM1 = TRUE)` bypasses the bounds check entirely.** `build_bounds()` runs at
  `R/6-fit_mod.R:915`; the `updateM1` refill of `log_M1` is ~380 lines later, so nothing it
  writes is ever bounds-checked. Measured on `BS2017SS` with `M1_base` set to 5.0/yr against an
  upper bound of `log(2)`: `updateM1 = FALSE` correctly errors "Initial parameter values are not
  within bounds: log_M1", while `updateM1 = TRUE` builds happily with `log_M1 = 1.609438`, handing
  `nlminb` a start outside its own `upper`. 5.50.0 added a non-finiteness refusal on that path but
  not a bounds check. Found reviewing PR #197, 2026-10-04.
- **`ceattle.cpp`'s linkage intercept-prior re-target reads `log_M1` with no clamp.**
  `b = log_M1(sp_idx, sx_idx, ab_idx)` under `RCEATTLE_PROC_M` takes the linkage row's stratum ids
  directly, unlike `rceattle_stratum_range()` (`linkage.hpp`), which clamps `hi` to `n_levels`. It
  stays in range only because `fit_mod()` builds per-species `sex` / `age_bin` strata and
  `linkage_spec()` exposes no `age_bin` argument, so it is not reachable through the public API
  today. Worth a one-line clamp anyway: before 5.50.0 an out-of-range read returned `NA` and gave
  a NaN objective (loud); now it returns `log(1) = 0`, so a prior would be evaluated silently at
  M = 1.0 per year. Found reviewing PR #197, 2026-10-04.

Added 2026-10-08, from three occurrences in one session:

- **A markdown table in this file silently truncates on an unescaped pipe, and nothing checks
  it.** In GitHub-flavoured markdown a cell separator is a `|` NOT preceded by a backslash, so
  `` `max|gradient|` `` inside a cell splits the row into five cells against a three-column header
  and **GitHub discards everything after the third** -- on one Tier 0 row that hid ~700 characters,
  including the gradient half of a retraction and the paragraph naming the live GOA multispecies
  assessment. The raw text an agent `Read`s is complete, so this is invisible to every reader who
  does not open the rendered file, which is the one a human uses. It has now happened three times:
  two pre-existing rows (fixed in 5.54.3), and once more while filing the shared-block linkage row
  above. The inverse also bites -- counting RAW pipes rather than unescaped ones reports a
  correctly-escaped row as broken: under the unescaped rule the pre-fix file has **2** malformed
  rows, under a raw-pipe count **4**. Neither is the "3" that a review and my own check both
  reported, so that figure was wrong twice over and is not reproducible from either rule.
  **The check is about fifteen lines in `tools/ci/source-guards.R`**, which already runs per PR
  against a real checkout: for every contiguous block of `|`-leading lines in `inst/dev/*.md`,
  assert each row has the same count of `(?<!\\)\|` as the block's first row. Key it on the
  block, NOT on a `| Where |` header -- that header appears in only 2 of the ~31 markdown tables
  under `inst/dev/`, and would miss `SESSION_HANDOFF.md`'s own PR table, `SIBLING-REPOS.md`'s five
  and `TRAPS.md`'s six. It belongs there rather than in an agent's script,
  because the two things that caught it both times were assertions a human would not have written
  by hand. Cheap, and it ends a footgun that has already cost one retraction's visibility.

## `TODO(review)` — Grant's calls, not an agent's

Six, each a judgement about what the right behaviour *is*:

- `R/0-rceattle_class.R` (`osa_residuals("all") includes diet`) — whether
  `residuals(source = "all")` should include diet too.
- `R/0-rceattle_class.R` (`held-out rows (Year <= 0) with a positive observation`, twice) — how
  those rows should be treated.
- `R/6-fit_mod.R` (`a user-supplied NA`) — what an `NA` bias-adjustment should mean.
- `R/7-plot_osa.R` (`process-residual objects`) — how `process_residuals()` output should be
  plotted.
- `src/TMB/growth.hpp` (`this branch (and its Richards mirror below) tests`) — see the file.

A seventh, on multispecies SBF, was restated as a known limitation in `5d423172`; it is settled
under "Deliberately not changed".

## Deliberately not changed

- **`data_check()` decomposition: the harness is built, the coverage is not there yet**
  (2026-10-05). The refactor plan's Step 8 says build
  `tools/verify/verify-data-check-conditions.R` FIRST and prove it deterministic before touching
  the 2,212-line function. Built, and determinism proven the way the plan asks: two separate R
  processes give byte-identical digests, and a selftest runs the capture twice in one process and
  refuses to pass if they differ. **But the coverage is 12 distinct conditions against 36
  stop/warning/message sites in `R/1-data_check.R`** -- roughly a third. Demonstrated
  insufficient: silencing the first `warning()` in the file leaves the digest IDENTICAL, because
  no case reaches it. So the net cannot gate the decomposition yet, and starting one behind it
  would be the worst outcome the plan names -- half the nets, none of the content. To close the
  gap: cases for the unreached sites, or a capture over the 375-workbook release corpus (183 with
  a `fleet_control` sheet). The harness reports its own coverage on every run so the number
  cannot quietly rot.
- **The three mechanical loop/apply conversions, measured and declined** (2026-10-05, counts
  corrected 2026-10-05 after the PR #203 review). The refactor plan listed `1:n` -> `seq_len()`,
  `sapply` -> `vapply` and `=` -> `<-` as conversions that "remove live bug classes". Measured
  against `dev` at 5.49.7 **by parsing the `for` call's sequence argument**, which is what makes
  the count trustworthy: an earlier grep said 103 and was wrong.
  - Bare `T`/`F`: **0 sites**. Nothing to do.
  - `=` as assignment: **111**, against `<-` everywhere else. (I have quoted a precise
    `<-` count before; two parse-based methods give 7,131 and 7,679 depending on whether
    `<<-` and function-argument defaults are counted, so the denominator is not worth a
    figure. The ratio is what matters and it is overwhelming.) Cosmetics.
  - `for (x in 1:expr)`: **98**. **55** bounds are dimension counts the schema or `data_check()`
    guarantees positive (`nspp`, `nsex`, `nages`, `nsim`, `nlengths`, `nyrs`); **43** could in
    principle be 0. Of those 43, the three that read an optional sheet are each already inside an
    explicit `if (nrow(...) > 0)` -- `emp_sel` (`R/5-rearrange_data.R:593`), `NByageFixed`
    (`:765`), `ration_data` (`:820`) -- so the `1:` is unreachable rather than safe, and the
    author clearly knew the hazard. Two are unguarded, both `nrow(data_list$M1_base)`
    (`R/2-build_params.R:99`, `R/6-fit_mod.R:1291`); a zero-row `M1_base` there fails **loudly**
    with "missing value where TRUE/FALSE needed" after writing nothing, and `data_check()`
    refuses a zero-row `M1_base` upstream. `nrow(age_trans_matrix)`
    (`R/5-rearrange_data.R:633`) is covered by the `stop()` four lines above it. The `yr_ind`
    case is harmless: `1:0` gives `c(1, 0)` and `m[NA, NA] <- NA` leaves the matrix
    completely unchanged (checked -- no cell written, no NA introduced). It is in
    `build_map_f_and_data_weights()` (`R/3-build_map.R:1676`), not
    `build_map_selectivity()` as first filed.
  - `sapply`: **34**. The one in the riskiest place, the bounds dimension check at
    `R/6-fit_mod.R`, behaves identically under `vapply` including on an empty list.

  So the conclusion stands but the reason is different from the first write-up: the zero-capable
  bounds are guarded or fail loudly, not structurally positive. The 243 edits across the fitting
  pipeline, where `CLAUDE.md` warns a rewrite is how a fit moves silently, buys nothing
  measurable. Convert opportunistically in a file being edited for another reason; do not sweep.


- **Multispecies SBF sits on the projection's realized M2** (`src/TMB/ceattle.cpp`, `Multispecies:
  M_at_age carries the projection's realized M2`). Under `msmMode > 0` it is reported but nothing
  live reads it. The one rule that does, NPFMC (5), is refused there along with 4 and 7;
  ConstantFSSB tunes realized SSB against SB0, CMSY reads depletion, and `mse_summary()` reads SBF
  only when `msmMode == 0`. Allowing HCR 5 in multispecies mode would make this a defect;
  `test-switches-hcr-multispecies.R` pins the refusal.
  **Measured, so the size is known before anyone reopens it** (`BS2017MS`, `estimateMode = 4`):
  pointing `NByageF` at `M_at_age_dBF` moves `SBF` by 5.27, and pointing `NByage0` at
  `M_at_age_dB0` moves summed `NByage0` by 16808.7. The scope is six new arrays and six new
  solver parameters inside the `iter` loop. Settle what the reference point should MEAN under
  predation before writing any of it -- an unfished equilibrium whose M2 comes from a fished
  projection is not one definition or the other. **The proposed answer**, for whoever picks it
  up: `NByage0` is not a strict equilibrium anyway (mean recruitment and terminal-year weights,
  but year-specific M1 and R0), so "equilibrium M2" here means dynamic-B0 without the
  recruitment deviations.
- **Non-parametric growth** is declared and calls `error("not yet implemented")`.
- **The `msmMode` 3–9 Kinzey-Punt branches are not declared at all** -- the whole block in
  `predation.hpp` is inside a `/* ... */`, so there is no dispatch, live or erroring. The live
  modes are handled by `if (msmMode == ...)` in `ceattle.cpp`.
  `test-schema-cpp-dispatch.R` pins both, and pins the absence of the switch.
- ~~`flt_sel_ind`~~ — removed in 5.12.0. It was computed from `Fleet_code` on every fit and read
  by nothing.
- **The `dmultinom_osa()` renormalization under `Comp_distribution` case 0.** Fitting routes
  through `dmultinom_osa()`, which renormalizes `p`, so the *reported* multinomial NLL carries a
  per-row constant the old `dmultinom()` did not. The gradient and the MLE are unchanged, so this
  is a reporting discrepancy, not a wrong fit. Correcting it would move the golden reference
  numbers for a cosmetic gain; reviewed and left in place 2026-08-23.
