# SS3's initial equilibrium catch, and how Rceattle fits it

**Status: SHIPPED in 5.46.0.** Written as a proposal while bridging the 2024 AI
Pacific cod SS3 assessment, and kept for its SS3 reading (section 2) and the
`styr - 1` marker trap. What a *user* needs is in
`vignette("model-options-and-functionality")`, "The initial equilibrium catch";
sections 3-5 below record why the shipped design is what it is, and where they
disagree with the code, the code is right.

Before 5.46.0 Rceattle had no equilibrium-catch concept anywhere in `R/` or
`src/TMB/`. SS3 has one, most AFSC SS3 assessments supply it, and it was the
single largest term in the AI cod gradient mismatch.

All SS3 references are pinned to tag **v3.30.22.1**, the version the AI cod model
was run with. The unversioned `ss3-source-code-main` checkout under `GOA cod/`
reports `#V3.30.xx.yy` and must not be cited.


## 1. Why it is worth doing

An initial equilibrium catch is the catch the stock would have yielded in the year
before the hindcast, under the initial fishing mortality that produced the starting
age structure. It is what makes `Finit` an estimable quantity rather than a shape
parameter of the initial numbers.

Measured on AI cod at SS3's MLE (`SS3-bridge/ss3_finite_difference.R`, which runs
SS3 as a pure forward pass and central-differences every reported component):

```
  d(component)/d(log ...)      K       L1     Linf   Richards   log InitF
  Equil_catch               3.861    0.836   6.617    0.903       1.270
```

Its **value** is 0.0035 nats, which is why dropping it looked free. Its **gradient**
is up to 6.6. Restoring it would move Rceattle's worst growth gradient at SS3's MLE
from 5.22 to about 1.40.

The `log InitF` column is the sharpest evidence. SS3's `Equil_catch` contributes
+1.270 to `d(objective)/d(log InitF)`; Rceattle's leftover `log_Finit` gradient at
that same point is -1.2700. Dropping the observation leaves exactly its own gradient
behind. That is the flat `Finit`/`init_dev` ridge G3 walks along: `log_Finit` moves
7.5x for 0.17 nats because nothing else in the likelihood is holding it.

**A component can be negligible in value and decisive in gradient.** That is the
general lesson and it is worth keeping.


## 2. What SS3 does

### 2.1 The observation

A catch row at year -999 in the data file. AI cod `data.ss` line 44:

```
#_Catch data: yr, seas, fleet, catch, catch_se
-999 1 1 7654 0.05
```

7654 metric tons with a log SD of 0.05, on fishery fleet 1. One value per fleet per
season; `nseas = 1` for both cod stocks.

### 2.2 The prediction

`SS_popdyn.tpl:2066`, inside `Do_Equil_Calc()`:

```cpp
equ_catch_fleet(3, s, f) += Hrate(f, t)
  * elem_prod(equ_numbers(s, p, g)(0, nages), sel_ret_bio(s, f, g)) * Zrate2(p, g);
```

with `Zrate2 = (1 - exp(-seasdur * equ_Z)) / equ_Z` (`SS_popdyn.tpl:2053`). That is
the Baranov catch equation on the equilibrium age structure:

> C_eq = sum_a Finit * s_a * w_a * N_eq,a * (1 - exp(-Z_a)) / Z_a,
> with Z_a = M_a + Finit * s_a

The third index is chosen by `i = 3 * catchunits(f)` in the likelihood, and the
inline comments give the meaning: **3 = retained catch in biomass**, 6 = retained
catch in numbers. AI cod fleet 1 has `catch_units = 1`, so biomass.

Two properties of `equ_numbers` matter and are easy to get wrong:

- **It carries no recruitment deviations.** `Do_Equil_Calc(equ_Recr)` seeds it from
  a scalar `equ_Recr` (`SS_popdyn.tpl:1842`). SS3's `Early_InitAge` deviates are
  applied afterwards to produce year-1 `natage`. So the equilibrium catch is *not*
  the catch implied by the year-1 numbers Rceattle already has; it is the catch
  implied by the deviation-free equilibrium behind them.
- **The plus group uses the selected F.** `equ_numbers(1,p,g,a+1) = Survivors /
  (1 - exp(-equ_Z(nseas,p,g,nages)))` (`SS_popdyn.tpl:2023`), and `equ_Z` at the
  oldest age is `M + Finit * s_nages`.

`equ_Recr` for the initial equilibrium is set by the control-file switch
`init_equ_steepness` (`SS_popdyn.tpl:486-521`):

- `0` -- `R1 = R0` (times any regime shift). No stock-recruit feedback.
- `1` -- `R1` is solved from the spawner-recruit curve at the SPR implied by
  `init_F`, via `Equil_Spawn_Recr_Fxn`.

AI cod `control.ss:121` sets it to 0.

**Rceattle already mirrors this switch.** `ceattle.cpp:1770/1776/1783` set
`R_init = R0` for `srr_fun` 0 and 1, and `ceattle.cpp:1801` sets
`R_init = (alpha - 1/SPRFinit)/Beta` for Beverton-Holt. Nothing new is needed here.

### 2.3 The likelihood

There is **no separate equilibrium-catch likelihood**. The block at
`SS_objfunc.tpl:696-708` that looks like one is commented out in v3.30.22.1, as is
the `est_equ_catch` assignment at `SS_popdyn.tpl:539-554`. The live path is the
ordinary catch loop, `SS_objfunc.tpl:711-734`:

```cpp
for (y = styr-1; y <= endyr; y++)
{
  temp = 0.5 * square((log(1.1 * catch_ret_obs(f, t))
        - log(catch_fleet(t, f, i) * catch_mult(y, f) + 0.1 * catch_ret_obs(f, t)))
        / catch_se(t, f));
  if (y == styr - 1) {equ_catch_like(f) += temp;} else {catch_like(f) += temp;}
}
```

The loop simply starts one year early and books that iteration to a different
accumulator. So the equilibrium catch uses **the same density as the hindcast
catch**, including SS3's 10% offset (`log(1.1*obs) - log(pred + 0.1*obs)`, which
weights catch 1.21x lighter than a plain lognormal -- see `NEWS.md` and
`SS3-bridge/ss3_to_rceattle.R` for how the converter carries that into the SDs).

Its lambda is `init_equ_lambda` (`SS_objfunc.tpl:1073`), like-component code 9. AI
cod leaves it unset, so 1.


## 3. The design, as shipped in 5.46.0

One thing below was overtaken by the implementation: the row is read **only** under
`initMode = "FishedNonEquilibriumSelected"` (6), not under every mode that
estimates `Finit`. Reading it under 3 or 4 fits `Finit` to a catch the population
was never subject to, and under a mode holding `Finit` at 0 the predicted catch is
0 and the likelihood takes `log(0)`.

### 3.1 Where the observation lives

**Recommended: a row in `catch_data` at `Year = styr - 1`.** It adds no columns.
`catch_data` already has exactly the right fields -- `Fleet_code`, `Year`, `Catch`,
`Log_sd` -- an equilibrium catch is a catch observation with a year, and
`write_data()`/`read_data()` round-trip it for free (a `data_list` element without
that support round-trips to nothing).

The hazard is that it is silently droppable: anything filtering the hindcast window
loses it without a word, which is exactly the failure that hid this for so long. The
converter's `build_catch_data()` drops it today at
`cat_raw$year >= datlist$styr`. Mitigations, both required:

- `data_check()` reports the equilibrium catch rows it found, per fleet, by name.
- A row at `Year < styr - 1`, or more than one per fleet, is refused rather than
  ignored.

Note that negative `Year` values are already reserved: `run_mse()` splices
negative-Year catch rows back in as the next assessment's data. `styr - 1` is a real
year and does not collide with that.

> **RESOLVED, 5.46.0, and this paragraph was half right.** `styr - 1` does not collide
> with the negative-Year reservation, but it does collide with ordinary catch history:
> `GOA2018SS` carries 23 catch rows before `styr`, two of them on 1976. Reading those
> as equilibrium observations under the default `initMode` predicted 0 and took
> `log(0)`, so every GOA fit returned a non-finite objective and the golden references
> failed with `optimHess: non-finite value supplied by optim`.
>
> A `-999` sentinel was tried and reverted: `abs(-999) = 999` passes `run_mse()`'s
> `abs(Year) <= endyr` window filters, so the row survives as year 999 — exactly the
> reservation this paragraph warns about.
>
> What shipped is `styr - 1` **gated on `initMode`**: the row is read only under a
> mode that estimates `Finit`, `data_check()` names the fleets whose rows it read
> whenever it reads any, and two rows for one fleet or a non-positive value are
> refused. Under any other mode the row is history and is dropped as before, which is
> what keeps `GOA2018SS` and the golden references unchanged.

**Alternative: `fleet_control$Equilibrium_catch` + `Equilibrium_catch_sd`.** Two new
columns, but immune to the silent-drop hazard, and defensible because the prediction
path is separate from the hindcast catch loop anyway. Rejected as the default
because it files an observation among settings and because it proliferates columns.
This is a schema decision and is Grant's call, not mine.

### 3.2 The prediction in TMB

The equilibrium numbers are already computed, minus the deviations. In `ceattle.cpp`
section 6.5 the initial numbers are

```cpp
N_at_age(sp, sex, age, 0) = R_init(sp) * exp(-mort_sum(sp, age) + init_dev(sp, age-1)) * sex_ratio
```

so the deviation-free equilibrium is the same expression with `init_dev` dropped.
It should be accumulated in the same loop rather than recomputed, so the two cannot
drift apart:

```
N_eq(sp, sex, age) = R_init(sp) * exp(-mort_sum(sp, age)) * sex_ratio(sp, sex)
```

with the existing plus-group geometric series, then

```
equil_catch_hat(flt) = sum_{sex, age} Finit(sp) * sel_init_flt(flt, sex, age)
                       * weight_hat(sp, sex, age, 0) * N_eq(sp, sex, age)
                       * (1 - exp(-Z_eq)) / Z_eq
Z_eq(sex, age) = M1_at_age(sp, sex, age, 0) + Finit(sp) * sel_init_flt(flt, sex, age)
```

Units: metric tons, matching `catch_data$Catch`. The weight-at-age is the fleet's
own `Wt_index`, as the hindcast catch uses.

Two things to be careful of:

- **Selectivity must be the fleet's own, not the `sel_init` average.** `sel_init`
  (added for `initMode = "FishedNonEquilibriumSelected"`) is the mean over the
  species' fishery fleets, because the initial age structure sees one pooled initial
  F. The *catch* is attributed per fleet, so each fleet's own `sel_at_age(flt,...)`
  in year 0 belongs in its own prediction. Using `sel_init` for both would be wrong
  wherever a species has more than one fishery fleet.
- **`Finit` is apical.** Selectivity is normalised to a maximum of 1, so `Finit` is
  the initial F at full selection, and the per-fleet split of a pooled `Finit` is an
  assumption, not a derivation. With one fishery fleet (both cod stocks) it does not
  arise. With several, refuse rather than guess -- see 3.5.

### 3.3 The likelihood

Same density as the hindcast catch, on the fleet's `Catch_distribution`, reading
`Log_sd` off the row. No new distribution, no new weight column.

A new `jnll_comp` row, **appended** as

```cpp
JNLL_EQUIL_CATCH = 21,   // Initial equilibrium catch
JNLL_N_ROWS      = 22
```

Appending rather than inserting after `JNLL_CATCH` matters: inserting renumbers rows
2-20 and silently invalidates every saved fit and every integer index. The enum has
two hand-synced partners that must be updated in the same commit -- the display name
in `R/6-rename_output.R` and `.JNLL_ROW_AXIS` in `R/9-profile.R`, where the new row
is registered as counting **fleets**. `test-schema-jnll-rows.R` reads the template
and asserts all three agree.

Folding it into `JNLL_CATCH` instead would avoid that, but it would also make the
bridge unverifiable against SS3, which reports the two separately, and would hide a
term whose whole lesson is that it is invisible in the total.

### 3.4 Interaction with `initMode`

**Read only under `initMode = 6` (`FishedNonEquilibriumSelected`).** The
prediction is Baranov at `Finit * s_a`, and mode 6 is the only mode that builds
the initial age structure with that mortality: mode 3 accumulates a flat
`M1 + Finit` over ages and mode 4 adds `Finit` once outside the cumulative loop.
Scoring the observation under 3 or 4 fits `Finit` to a catch the population was
never subject to, and it converges, so nothing says so.

That was the shipped behaviour for one review round and is recorded here because
`SS3-bridge/ss3_to_rceattle.R` and `run_g3.R` still set modes 3 and 4: they now
read no equilibrium catch, and `data_check()` says so. The AI-cod parity run in
section 4.3 needs `initMode = 6` to reproduce.

### 3.5 Refusals

Per the house rule that a plausible placeholder surviving into a fit is the worst
failure mode here, refuse rather than default:

- an equilibrium catch on a species whose `initMode` starts unfished (1, 2, 5) --
  `Finit` is 0, so the prediction is identically 0 and the row cannot be fitted;
- an equilibrium catch on a non-fishery fleet;
- more than one equilibrium catch row per fleet;
- an equilibrium catch on a species with more than one fishery fleet, until the
  per-fleet split of a pooled `Finit` is specified. Neither cod stock hits this,
  and inventing a split is exactly what rule 9 forbids.

### 3.6 Paperwork

Behaviour change, so `NEWS.md` + `DESCRIPTION` `Version:` + the affected vignette in
the same commit, and `_pkgdown.yml` if a documented topic appears. Minor bump, not
breaking: a `data_list` without an equilibrium catch row fits exactly as before.


## 4. Verification

1. **Unit.** A one-fleet, flat-selectivity, known-`Finit` fixture where the Baranov
   sum is computable by hand.
2. **Zero case.** No equilibrium catch row leaves `jnll` bit-identical and
   `jnll_comp[JNLL_EQUIL_CATCH, ] == 0`. This is the golden-reference guard; all
   four reference models must stay within tolerance.
3. **The real gate -- SS3.** Re-run `SS3-bridge/run_parity.R "AI cod - Dev"` after
   restoring the year -999 row in `build_catch_data()`. Expected, from the finite
   differences already measured:
   - `Equil_catch` residual after the density constant near 0 against SS3's 0.0035;
   - `log_Finit` gradient from -1.2700 to near 0;
   - worst growth gradient from 5.22 to about 1.40.

   The third is the honest prediction and the one most likely to be wrong; the first
   two are close to arithmetic.
4. **Identification.** With the row in, re-run G3. `log_Finit` should stop wandering
   7.5x for 0.17 nats. `verify-` harness, since `/golden-check` covers none of the
   refit paths.
5. **Round-trip.** `write_data()` then `read_data()` must preserve the row.


## 5. A defect found alongside this, and FIXED with it

`initMode = "FishedNonEquilibriumSelected"` (6) decays the initial numbers with
`M1 + Finit * sel`, and its plus-group geometric series now divides by the same
selected F. `ceattle.cpp:2097-2099` builds
`Z_plus = M1 + Finit * sel_init(sp, sex, nages-1)` under mode 6 and charges the
full `Finit` under every other mode, which is SS3's convention
(`SS_popdyn.tpl:2023`) and leaves modes 1-5 bit-identical. Before the fix it
divided by `1 - exp(-M1 - Finit)` -- the **unselected** `Finit`.

On AI cod the fishery's selectivity at the plus group is 0.9999, so the plus-group
factor moved from 2.621575 to 2.621589, **+0.001%** -- real but not what was moving
anything here. It matters for a dome-shaped fishery whose descending limb leaves
the oldest ages lightly selected. Fixed alongside the equilibrium catch because
both touch the same block; it was not the cause of anything measured above.


## 6. Future work: per-fleet initial F, to match SS3 on a multi-fishery species

**Status: AGREED, NOT STARTED** (Grant, 2026-10-02). Deliberately deferred — it
changes the length of a parameter vector, so it wants its own PR and a
deprecation path rather than a slot in a release that already moves SSB.

SS3 gives each fleet its own `init_F` and **sums** them over the species'
fisheries:

```
for (int ff = 1; ff <= N_catchfleets(p); ff++) {
  f = fish_fleet_area(p, ff);
  equ_Z(s, p, g, a1) += sel_dead_num(s, f, g, a1) * Hrate(f, t);
}
```

(`SS_popdyn.tpl:1990-1994` in the local v3.30.25.1 checkout.) One parameter per
fleet x season, created only where that fleet's equilibrium catch is positive
(`SS_readcontrol_330.tpl:2837-2847`, FATAL if the catch is positive and no
`init_F` exists).

Rceattle has `PARAMETER_VECTOR(log_Finit)` of length `nspp` (`ceattle.cpp:514`)
and applies it to the unweighted MEAN selectivity across the species' fisheries
(`ceattle.cpp:1536`, `sel_init = acc / n_fsh`). **The two cannot disagree
today**, because `.check_equil_catch()` refuses a species with more than one
fishery — so this is about lifting that refusal, not about correcting a number
in the current bridge. The AI cod bridge already reproduces SS3 at a year-1
`SSB_rce/SSB_ss3` of 1.0000 with 89 free parameters against SS3's 89 active,
using a combined `FshComb` fishery.

What it costs, and why it is not folded in here:

- `log_Finit` becomes per-FLEET, so the parameter vector changes length for
  every `initMode` 3/4/6 model. Stored fits stop refitting and an
  `inits$log_Finit` of length `nspp` no longer maps — rule 1 wants a
  deprecation path (broadcast the old per-species value to the species'
  fisheries, or refuse it with a message naming the refit).
- `build_map()`, `build_parameter_bounds()`, `rename_output()` and the
  parameter dictionary all carry the shape.
- It also closes a latent inconsistency worth recording: `equil_catch_hat`
  already uses the FLEET's own `sel_at_age` in the Baranov numerator
  (`ceattle.cpp:2962`) while `N_eq` was depleted at the MEAN selectivity
  (`:2088`). Harmless while >1 fishery is refused; wrong the moment it is not.
- Rceattle has no retention concept, so `sel_at_age` stands in for SS3's
  `sel_dead_num = sel * (retain + (1-retain)*discmort)` (`SS_param.tpl:532`).
  Identical with no discards; worth a line in the vignette either way.

Accept per-fleet `init_F` only where that fleet has a positive equilibrium
catch, as SS3 does, or the parameter is unidentified.
