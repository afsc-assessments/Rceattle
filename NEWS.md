<!--
Version-numbering note. Three gaps in this file are deliberate, not lost entries:

  * 4.14.0 was a real DESCRIPTION version whose entries were folded into 5.0.0.
  * 5.2.0-5.2.4 were likewise folded into 5.3.0.
  * main's 4.9.0 / 4.9.1 are the same recruitment changes this line has as
    5.5.0 / 5.5.1, applied to the two lines separately.

No tag existed above 4.8.0 while these were in flight, so nobody could have installed an
intermediate. They were folded rather than renumbered because renumbering a section breaks
every (x.y.z) cross-reference pointing at it, and the entries below cite each other by
version throughout.
-->

# Rceattle 5.43.0

## Breaking changes
* **The package moved to `afsc-assessments/Rceattle`, and every link follows
  it.** `DESCRIPTION`'s `URL` and `BugReports`, `_pkgdown.yml`'s `url` and
  navbar, the README's install commands and badges, `?print.Rceattle`'s install
  hint, `CONTRIBUTING.md`, `examples/Install_Rceattle.R` and two vignettes all
  named `grantdadams`. GitHub redirects the repository path, so the install
  commands kept working -- but **the documented website did not**:
  `https://grantdadams.github.io/Rceattle/` returns 404 and
  `https://afsc-assessments.github.io/Rceattle/` returns 200, so `URL`,
  `_pkgdown.yml` and the README all pointed users at a dead page. Pin
  `afsc-assessments/Rceattle@5.43.0`: a redirect is not a permanent address, and
  it breaks the moment anything is created at the old path.


* **An integer-coded `Fleet_type` no longer reads as estimated in
  `data_check()`.** `data_check()` is callable on a list straight from
  `read_data()`, where the switch columns are still the integer codes the
  workbook stores -- and every bundled data set that carries a `fleet_control`
  stores them that way (`GOA2018SS` is `2, 2, 2, 2, 2, 2, 0, 1, ...`), so the
  integer is the shipped representation. `0 != "Off"` coerces to
  `"0" != "Off"`, which is `TRUE`, so an `Off` fleet read as estimated: on
  `GOA2018SS` with fleet 7 given its own `Selectivity_index`, it was named in
  the "estimated Selectivity but no `comp_data`" error, and is not now. The
  subset reads the column through `.canon_switch()`, as the same file already
  did eleven lines further down, rather than through a fourth spelling of the
  test -- `.canon_switch()` also resolves `" 0 "`, `"0.0"` and `"00"`, which a
  `%in%` list does not.

  **No fitted number moves.** `fit_mod()`, `build_map()` and `build_params()`
  each call `switch_check()` first, and on the canonical strings the old and
  new readings are `identical()` -- asserted for five bundled data sets. Only
  the un-canonicalized path changes.

  The rest of the class is open, in `inst/dev/CLEANUP_BACKLOG.md`: the schema
  types thirteen columns as `switch` with an `allowed` map, but nothing
  enforces that at the boundary, so any comparison written against the
  canonical spelling is wrong for a `fleet_control` that has not been through
  `switch_check()`. Those sites are unreachable today because their callers
  canonicalize first; the fix is to enforce the schema once on entry rather
  than to convert them one at a time.

* **The distance-to-optimum report is withheld when the covariance cannot
  carry it.** `max_gradient` reports the Newton step `cov.fixed %*% gradient`
  in standard errors, which measures the distance to a minimum only if
  `cov.fixed` inverts a positive-definite Hessian. On a saddle the step points
  away from one, so "reaching the optimum would move the estimates by at most
  0.0009 standard errors" read as reassurance printed directly beside a failing
  `pdHess`. It is reported only when `pdHess` is `TRUE`. Severity is still read
  on the gradient, so no fit changes status -- now asserted, where neither that
  claim nor the laziness of the coordinate lookup had a test.

  A negative variance no longer escapes as `NaNs produced`. Three sites took
  `sqrt(diag(cov))` unguarded, in a battery documented never to raise; they now
  share `.conv_se_from_cov()`. The callers already dropped non-finite entries --
  it was the `sqrt()` itself that warned.

  `hessian_conditioning`'s coordinates line up again. The column width was
  computed inside `.rce_par_summary()`, which `.check_hessian_eigen()` calls
  once per block, so each line sized to its own block name and the aligned
  output `vignette("model-diagnostics")` illustrates was unreachable. The caller
  now passes one width, and the vignette's example is regenerated from the
  format `print()` emits rather than hand-written.

* **A blank `Fleet_type` is refused.** The column has no schema default, so
  nothing filled it, and a blank one is not "unset, take the default" -- it is
  a fleet whose role in the likelihood nobody stated. `switch_check()` now
  stops, naming the fleet and the values to choose from, and
  `rearrange_data()` does the same: it is exported and reachable without
  `switch_check()`, and a blank reaching it was handed to the template as
  `flt_type = NA`.

  `NA` was the damaging case, and the two halves of the package read it
  differently. `Fleet_type != "Off"` is `NA` rather than `FALSE`, so
  `data_check()`'s estimated-selectivity subset kept an all-`NA` row and died
  inside its `vapply` on `missing value where TRUE/FALSE needed`, naming no
  fleet and no column -- measured on `Atka2022` with `Fleet_type[1] <- NA` --
  while `build_map_selectivity()` treated `NA` as estimated
  (`.on[is.na(.on)] <- TRUE`). An empty string or whitespace already produced
  a correctly named error from `validate_switches()`; those are refused here
  too, earlier and from one rule, rather than because they were broken in the
  same way.

  Every spelling of a stated type -- canonical name, integer code, character
  code, and `Off` -- is accepted exactly as before. No bundled data set and
  none of the 183 consumer-repository workbooks carrying a `Fleet_type` column
  has a blank one, so nothing that fits today stops fitting.

  Still open, and tracked in `inst/dev/CLEANUP_BACKLOG.md`: an
  **integer-coded** `Fleet_type` read before `switch_check()` canonicalizes it
  reads as live at every bare `!= "Off"` comparison, because `0 != "Off"` is
  `TRUE`. `validate_switches()` documents the trap and canonicalizes first;
  fourteen other raw comparisons in `R/` do not.

# Rceattle 5.42.1

Found reviewing the release PR that carries 5.34.0 through 5.42.0 (#158), over
several adversarial passes. Corrections to 5.42.0's own guard, to the tests
meant to hold it, and to documentation it left behind.

## Bug fixes

* **The negative-`Sel_curve_pen` refusal no longer skips a fleet whose penalty
  the template does charge.** The check skipped any fleet that followed another
  fleet's `Selectivity_index`, on the grounds that a shared block is penalized
  once, on its lead. But the template's lead (`flt_sel_lead`, built in
  `rearrange_data()`) groups fleets by `Selectivity_index` **and** selectivity
  form, while the rule the check borrowed -- the one the parameter map shares
  on -- groups by the index alone. Two live fleets sharing an index with
  different forms are therefore two groups to the template, each charged, and
  the check called the second one a follower and passed over it: a negative
  weight there reached `ceattle.cpp` and rewarded the deviation it names,
  without bound, with `data_check()` reporting clean. The check now reads the
  lead the same way the template does. No bundled data set and no workbook in
  the consumer repositories has a group of mixed form (every bundled
  `fleet_control`, including the three carried inside fitted example objects,
  and 183 workbooks with a `fleet_control` sheet; the answer is the same read
  raw or canonicalized through `sel_map`), so no existing model changes;
  this closes the guard rather than moving a number.

  The lead it builds also **fails closed on a `Selectivity` the template cannot
  read.** `.canon_switch()` trims and `rearrange_data()`'s `.pull_int()` does
  not, so `" NonParametric"` resolves here and reaches the template as `NA` --
  a group of its own there, with the fleet leading and its weight charged.
  Such a value now gets a key of its own here too, so the fleet leads and is
  checked. `switch_check()` normalizes the spelling before either is reached,
  so this is a guard against a hand-built `fleet_control`, not a path a
  workbook takes.

* **`test-docs-anchors.R` now checks the schema column that shipped the stale
  code.** The guard added at 5.42.0, so that no schema description names a
  selectivity code `sel_map` does not accept, matched only a slash-separated
  run (`2/9/13`). Of the eight columns it names, four matched nothing at all --
  including `Selectivity` itself, which enumerates every code one per line and
  is the column form 14 was advertised in. It now also reads a comma-or-`or`
  list, a per-line `13 = ...` enumeration, and a code parenthesized after a
  form name, quoted or not -- `Sel_curve_pen1` writes `"type 2/9/13"` and
  `LogisticPM (11)` in one sentence, so requiring the quotes missed 11 and
  would have missed a stale code in that position. All eight columns now yield
  codes. It is a shape-matcher, not a parser: a code written in a shape none of
  the four patterns covers (`14 - ...`, `forms 2 and 14`) still slips, so it
  narrows the gap rather than closing it.

* **A refusal test that passed for the wrong reason.** The check that an
  `apical` linkage is refused on an AR1 selectivity form matched the string
  `"AR1"`, and `GOAatf` ships `Sel_curve_pen2 = 200` on the fleet the fixture
  switches: on the AR1 forms that slot is a correlation, so `data_check()`
  stopped first with its own out-of-range message, which also contains `"AR1"`.
  Deleting the refusal left the suite green. The fixture now zeroes that column
  and the assertion matches the apical clause, so removing the refusal fails the
  test.

* **`.rce_sel_pen_lead()`'s `Off` handling was untested.** Replacing its `off`
  argument with `rep(FALSE, n)` left `test-selectivity-penalty-sd.R` green,
  because the sign check has its own independent `Off` skip. The file now pins
  that an `Off` fleet never leads its group, that the next fleet leads instead,
  and that the result is `identical()` to the `flt_sel_lead` `rearrange_data()
  hands the template for the same table.

## Documentation

* `README.md`'s operational pinning example named 5.41.0, a version this line
  never releases, so an assessor pinning a version for management advice was
  sent to a reference `install_github()` can never resolve. It now names the
  version this release tags. Like every version named there it resolves only
  once that tag is pushed -- checklist section 3, and `TODO-pre-transfer.md`
  B5 says to re-check the line whenever the version moves.

* `vignette("model-parameterizations")` listed three cases the negative-weight
  refusal does not catch; there are four. The fleet following another's
  `Selectivity_index` was added to the check and to `NEWS.md` at 5.42.0 but not
  to the reference table.

* `vignette("adding-a-selectivity-form")` said code 14 was free and the next
  form should take 15, in a paragraph whose own rule is that a retired code is
  not free. 14 was advertised as a second integrable form in
  `meta_data_names.xlsx` -- the template every workbook is built from -- before
  being collapsed into 13, so a `Selectivity` column written in that window
  holds it. It is now listed with 10 and 12 as taken. The same article named
  three R sites that key on the selectivity form when there are at least eight,
  attributed two of them to the wrong check, and described
  `.RCE_SEL_PEN_POSITIVE` as driving the penalty-SD conversion; that conversion
  is a separate registry, `.sd_specs`, which fails closed where
  `.RCE_SEL_PEN_POSITIVE` fails open.

* `vignette("developer-guide")` gave three switch-code facts the code does not
  support: the `selectivity.hpp` dispatch list omitted the Ianelli
  non-parametric form (2), which every bundled reference model uses, and
  descending logistic (4); `sel_map` skips 10, 12 and 14, not 10 alone; and the
  `Sel_curve_pen1` / `Sel_curve_pen2` slot map omitted form 13.

* The 5.39.0 entry now says the stock-recruit bound is unconditional. Predation
  is why it was added, and it sits under a predation heading, but
  `build_bounds()` applies +/-30 to every model.

# Rceattle 5.42.0

Found reviewing the 5.34.0-5.41.0 release (`# Rceattle 5.41.0` and the versions
below it). The negative-weight defect is much older than that review: the
one-sided shape penalty has had no sign branch since it landed in December 2024,
and 4.10.0 added a second way in through `Sel_shape_dir`. 5.40.0 only extended it
to the integrable forms.

## Breaking changes

* **A negative `Sel_curve_pen1`/`2`/`3` is refused wherever the template reads
  it as a weight.** A penalty slot multiplies a squared deviation, so a negative
  weight rewards the deviation it names instead of penalizing it. Use a positive
  weight, or the matching `Sel_shape_sd` / `Sel_curvature_sd` / `Sel_devmag_sd`
  column, which is a standard deviation and cannot go negative. The one
  exemption is `"NonParametricPM"` (9) under `Sel_shape_mode = "Directional"`,
  the only branch in the template that reads a sign -- there it switches the
  penalty from the decreasing side to the increasing one rather than negating
  it. That form's `"Smooth"` mode applies the weight two-sided and is refused
  like the rest. `"2DAR1"` (6) and `"3DAR1"` (7) reuse these columns as AR1
  correlations and are untouched.

  Slot 3 is checked only on `"NonParametricPM"` (9) and `"LogisticPM"` (11):
  `"NonParametric"` (2) and `"NonParametricIntegrable"` (13) also estimate
  selectivity deviates, but score them with a Gaussian density on
  `Time_varying_sel_sd`, so neither reads the slot. The rule throughout is to
  refuse a weight only where the fleet's own wiring reaches it: the template
  gates the penalty block on `flt_type(flt) > 0 && flt_sel_lead(flt) == 1`, so a
  `Fleet_type = "Off"` fleet is skipped, and so is one that follows another
  fleet's `Selectivity_index` -- the group is charged once, on its lead.
  (**Narrowed in 5.42.1**: the template groups by `Selectivity_index` *and*
  selectivity form, so only a follower sharing the lead's form is skipped. A
  follower with a different form leads its own group, is charged, and is
  refused. Read that entry, not this sentence, for the rule in force.)

  Measured on `BS2017SS` fleet 1 (`NonParametric`, `N_sel_bins = 8`, with the
  shipped `Sel_curve_pen2 = 12.5` active): ramping the fleet's `sel_coff`
  downward by a constant step per bin takes the `JNLL_SEL_NONPARAM` row to -29,
  -122, -502 and -2032 at steps of 0.5, 1, 2 and 4, and the whole objective to
  -6.6e6 at a step of 256 -- quadratic, and falling without bound, while the
  composition likelihood flattens and cannot counteract it. The row is exactly
  `[Sel_curve_pen1 * (N_sel_bins - 1) + Sel_curve_pen2] * step^2` up to a
  constant, so the **curvature penalty does bind** -- on the single second
  difference where the ramp meets the coefficient repeated past `N_sel_bins` --
  and the objective diverges only once `|Sel_curve_pen1|` exceeds
  `Sel_curve_pen2 / (N_sel_bins - 1)`, 1.79 here. A smaller negative weight is
  merely anti-shrinking rather than divergent; it is refused too, because
  rewarding a deviation is wrong at any magnitude. `data_check()` now stops,
  naming the fleet, the slot and the standard-deviation column that sets the
  same weight safely.

* **`Sel_shape_dir = "Increasing"` is refused on every form but
  `"NonParametricPM"` (9) under `Sel_shape_mode = "Directional"`.** It used to
  negate `Sel_curve_pen1` on every non-parametric form. `"NonParametric"` (2)
  and `"NonParametricIntegrable"` (13) have no increasing-direction penalty at
  all -- the template hard-codes `max(d, 0)^2`, `d` being the log-selectivity
  drop from one bin to the next -- so `"Increasing"` never did what it said
  there; it fitted the reward above. There is no like-for-like migration: only 9
  in `"Directional"` mode implements the direction, and it is not a drop-in (it
  charges its average-selectivity term only when `Sel_avgsel_pen > 0`, defaults
  its penalty range to `Bin_first_selected` rather than the first bin, reads
  `Sel_curve_pen3`, and refuses `Time_varying_sel = "IID"`). If you did not mean
  an increasing penalty, drop the column -- the default is `"Decreasing"`.
  Unlike the sign check above, this one is not waived on a `Fleet_type = "Off"`
  fleet: the direction states what the form can express, and `write_data()`
  persists the weight it sets.

* **`Sel_devmag_sd` is refused on `"NonParametric"` (2) and
  `"NonParametricIntegrable"` (13).** It writes `Sel_curve_pen3`, which neither
  form reads: both estimate selectivity deviates, but score them with a Gaussian
  density on `Time_varying_sel_sd`. The column was a silent no-op on both.

* **A negative `sel_curve_pen` carried in `inits` is refused too.**
  `sel_curve_pen` is a parameter, not data, so `inits` from a stored fit
  override the `Sel_curve_pen` columns. Without this a fit saved before 5.42.0
  would keep its negative weight through `retrospective()`, `profile()`,
  `run_mse()` and `self_test()` on a workbook the user had already corrected --
  measured at -3.81 on a decreasing ramp and -28.99 at three times that ramp,
  with the column reading a valid +20. `fit_mod()` now applies the same rule to
  the parameter actually in use and names `inits` as the source. This closes the
  **sign** only: a stored `inits$sel_curve_pen` still supersedes an edited
  `Sel_curve_pen` column at any magnitude, silently (a column retuned from 200
  to 20 keeps 200 on a refit). That is tracked in `inst/dev/CLEANUP_BACKLOG.md`.

* **An `apical` selectivity linkage on a fleet with `Fleet_type = "Off"` is
  refused.** No data are fit to such a fleet, so `log_sel_apical` was a free
  parameter in a flat direction: a singular Hessian and a failed `getsd` with
  nothing naming the cause. It was the only selectivity parameter `build_map()`
  freed without a `Fleet_type` gate.

* **A selectivity *prior* on a fleet with `Fleet_type = "Off"` is refused**, for
  the same reason: `build_map_selectivity()` maps that fleet's `log_sel_slp` and
  `sel_inf` off, so the prior was evaluated against a fixed value and added a
  constant to the objective while constraining nothing. The apical linkage was
  given this gate above; the prior path had been left without one.

## Bug fixes

* **The negative-weight refusal no longer fires on a slot the fleet's
  configuration makes inert.** `LogisticPM`'s slots 1 and 3 and
  `"NonParametricPM"`'s slot 3 are charged on the time-varying deviates, so
  under `Time_varying_sel = "Off"` the term is identically zero and the weight is
  never read -- verified by objectives bit-identical at `+w`, `0` and `-w` on
  `BS2017SS`. A negative value there is inert rather than wrong, and is now
  allowed; it is still refused once the deviates are estimated. The slots charged
  on the base curve (1 and 2 on the non-parametric forms) are checked whatever
  `Time_varying_sel` says.

* **`data_check()` no longer requires `Sel_curve_pen2` on a time-varying
  `LogisticPM` fleet.** The template says in as many words that
  `sel_curve_pen(flt,1)` -- the 0-based C++ spelling of the same column -- is
  unused in that branch, so the column was required and
  then ignored -- and the schema, `?BS2017SS` and
  `vignette("model-parameterizations")` all say `LogisticPM` does not use it. The
  two weights it does read, `Sel_curve_pen1` (the random walk on realized
  log-selectivity) and `Sel_curve_pen3` (the walk on the age-1 deviates), are
  still required, and the message now names what each one weights.

* **An `apical` selectivity linkage is now refused on the fleets that actually
  share a block, and allowed on the ones that do not.** `Selectivity_index` is a
  group key, not a fleet code: the lead is the group's first fleet that is not
  `"Off"`, and a group of one shares nothing. The check compared the key to
  `Fleet_code`, which both refused offsets that were perfectly identifiable
  (`GOAatf` fleet 3 and `GOA2018SS` fleet 11, each the sole member of its group,
  and `GOA2018SS` fleet 9, the lead of the group `{9, 10}`) and, where a
  group's key equalled its *second* member's code, let the follower
  through: `build_map_selectivity()` freed the cell and
  `adjust_map_shared_params()` then mapped it off, pinning the offset at
  `exp(0) = 1` with the fit converging and nothing reported. Both now follow
  `.shared_block_lead()`, the rule the map itself applies, and the error names
  the lead fleet. The same comparison governed selectivity **priors** and is
  fixed with it. This tightens as well as loosens: a fleet whose `Fleet_code`
  happened to equal its group's key was always allowed before, so an `apical`
  linkage or selectivity prior on one that a lower-numbered fleet leads (group
  `{2, 5}` keyed 5, say) now stops. That block was the lead's, so the linkage
  was freeing nothing.

* **`osa_residuals()` records the `discrete = TRUE` method override.** A
  Gaussian method cannot score a discrete observation, so those composition rows
  fall back to `oneStepGeneric`. The fallback happened but was neither announced
  nor written to the returned object's `method` attribute, so the attribute --
  and `print()` -- named a method no composition row had used.

* **The non-finite-residual messages now name the `"cdf"` limitation where it
  applies.** Both the headline warning and `.osa_retry_tail()`'s message sent
  the analyst to re-check convergence and sample sizes without mentioning the
  measurement recorded in the same file: on a random-effects fit with a large
  composition data set the Laplace inner problem fails on the depth of
  conditioning, and the retry recovers nothing (1879 before, 1879 after). Both
  now add that reading, and only where it can hold -- the headline warning shows
  it under `method = "cdf"` on a fit with random effects, since a Gaussian
  method returns every residual finite on the same fit. Convergence and the
  sparsest compositions remain the first thing to check otherwise.

* **`plot()` on an `rceattle_osa` object warns when it drops residuals.** It
  filtered to the finite ones and warned only when *every* residual was
  non-finite. Under `method = "cdf"` the failures are a contiguous tail, so the
  survivors are a time-biased subset -- at the documented 1879-of-4538 loss the
  figure drew a clean Q-Q panel, with an SDNR annotation, and said nothing.

* **`convergence_diagnostics()` reports a `NOTE` under `getsd = FALSE`.** With
  no `sdreport` the Hessian eigenvalue, `pdHess`, sdreport and estimability
  checks all return nothing, and the battery reported `"OK"` --
  `report_tables()$model$converged` then printed `OK` into a SAFE table.
  Nothing distinguished "every check passed" from "the strongest checks never
  ran". The fit is unchanged. Only `estimateMode` `"Estimate"` and `"Hindcast"`
  run the battery, and each refitting diagnostic resolves its own `getsd`:
  `retrospective()`, `jitter()`, `self_test()` and `profile()` follow whether
  the source fit kept an `sdreport`; `reweight()` follows that fit's
  `fit_control$getsd`, so a default fit still reads `OK`; `run_mse()` refits at
  `getsd = FALSE`, so an MSE's estimation fits now read `NOTE` -- a status only,
  changing no result.

* **The `"NonParametricPM"` (9) directional shape penalty now has a
  limiting-case and specification check.** Nothing verified that
  `Sel_shape_dir = "Increasing"` penalizes an increasing curve on the one form
  that implements it; the existing tests
  covered the `weight = 1/(2*sd^2)` arithmetic and one `"Decreasing"` fit. On a
  strictly increasing curve the decreasing direction now charges 0 and the
  increasing direction 5.6, mirrored on a strictly decreasing curve, and all
  three shapes tested match the ADMB/AMAK `sel_like(1)` one-sided SSQ recomputed
  from `sel_coff` to 1e-10, driven both through `Sel_curve_pen1` and through the
  `Sel_shape_dir` column a user actually writes. This is why the direction is
  refused elsewhere rather than implemented there: on 9 it is now checked.

## Documentation

* The `estDynamics` list in `vignette("data-without-excel")` still offered code
  `3`, retired in 5.35.0 and refused by `switch_check()` since. The same file
  listed the codes correctly 344 lines earlier.

* `Sel_shape_dir` is documented as `NonParametricPM`-only in the column schema
  (which ships verbatim into `meta_data_names.xlsx`), in `?BS2017SS` and in
  `vignette("model-parameterizations")`. All three described the directional
  sign as working on every non-parametric form.

## Convergence messages name the quantity and the coordinate

* **`fit$convergence` names what a flagged parameter estimates.** The
  `max_gradient` message, the `hessian_conditioning` loadings, and the
  per-coordinate lines under `hessian_conditioning`, `parameters_on_bounds` and
  `estimability` now follow each block name with the natural-scale quantity
  from `parameter_dictionary()`: `log_M1 (M1)`, `rec_dev (recruitment
  deviations)`, `log_F (F)`. The block name stays first, since it is what
  `map` and the parameter list are keyed on.
* **`max_gradient` says where the largest gradient sits.** It gave the block
  alone (`'log_M1'`); it now gives the species, fleet, sex, age or year from
  `parameter_index()`, the way `estimability`, `parameters_on_bounds` and
  `hessian_conditioning` already did. The coordinate is resolved only when the
  check is not `OK`, so a clean fit still skips building the index.
* **`max_gradient` reports the distance to the optimum in standard errors**,
  as `$data$newton_step_se` and in the message. A quadratic approximation puts
  the optimum a Newton step `cov.fixed %*% gradient` away; dividing each
  element by its standard error makes the size comparable across log, logit
  and natural-scale parameters. It is reported only when the `sdreport`
  describes the hindcast parameters (not under an estimating HCR), and the
  severity is still read on the gradient, so no fit changes status.
* **Scattered years and ages read as a count against their span**: "38 years in
  1980-2021" rather than "38 of 1980-2021", and the `hessian_conditioning`
  count reads "(38 of 44 parameters; 67% of the direction)".

# Rceattle 5.41.0

## One-step-ahead residuals from the conditional CDF

* **`osa_residuals(method = "cdf")`.** The one-step-ahead residual is defined
  through the conditional CDF, `qnorm(F(x))`, and is standard normal by the
  probability integral transform whatever shape the conditional has. Until now
  the model could not supply that CDF, so only the Gaussian methods were
  available; those approximate the conditional as normal and standardize the
  observation against its conditional mode. The template now supplies it: the
  continuous binomial `1 - I_p(x + 1, n - x)` for a composition bin, which is
  defined at the fractional counts composition data hold, and `pnorm` for the
  aggregate index, catch and covariate series (`src/TMB/comp_osa.hpp`, gated by
  `keep.cdf_lower` / `keep.cdf_upper` under a new `osa_mode = 2`). The fitted
  objective is untouched: the gates are zero except inside a
  [TMB::oneStepPredict()] call, and the terms are not even evaluated outside it.

  **The composition residuals it returns are the only ones that pass a
  self-test.** Residualizing simulated observations at the parameters that
  generated them makes the answer exactly iid N(0, 1), so the methods can be
  scored rather than argued about. On BS2017SS over 20 replicates of 4538
  composition bins (`tools/verify/verify-osa-cdf.R`):

  | method | mean | sd | lag-1 acf within a composition | KS rejects |
  |---|---|---|---|---|
  | `oneStepGaussianOffMode` | +0.103 | 0.918 | +0.060 | every replicate |
  | `cdf`, `discrete = FALSE` | +0.610 | 1.262 | +0.434 | every replicate |
  | `cdf`, `discrete = TRUE` | -0.007 | 0.995 | +0.001 | none |

  The null standard error on the autocorrelation is 0.015. On the aggregate
  series, which are genuinely Gaussian, all three agree to 5e-5 and all three
  pass.

* **`discrete` now defaults per method:** `TRUE` under `"cdf"`, `FALSE`
  otherwise, which is what every method that existed before this one already
  did. The default changed from `FALSE` to `NULL` to express that; passing
  `TRUE` or `FALSE` explicitly still does exactly what it did, and `FALSE` under
  `"cdf"` now says in a message that those composition residuals are biased up.
  It has to be
  `TRUE` under `"cdf"`: a composition bin holds a count, so its conditional CDF
  is a step function and `qnorm(F(x))` inherits the step, which is why the
  middle row of the table above is the worst of the three rather than the best.
  Randomizing over the step removes it: `qnorm(F(x) - U f(x))`, Dunn and Smyth
  (1996), the construction Trijoulet et al. (2023) prescribe. The
  randomization is drawn serially under `seed` after the per-observation loop
  returns, so this method stays bit-reproducible with `parallel = TRUE`, and
  `attr(osa, "discrete")` records what was used.

* **`predicted` is `NA` on every row under `method = "cdf"`.**
  [TMB::oneStepPredict()] returns `Fx`, `px` and `nll` for this method and no
  fitted value, and nothing is substituted: the Gaussian methods' `predicted` is
  a conditional mode, and a marginal fitted value in its place would give the
  column two meanings depending on `method`. The Dirichlet-multinomial rows
  below run under a Gaussian method and are blanked with the rest, so the column
  means one thing across the object. Fitted values remain in `fit$quantities`
  and in `residuals(fit, type = "pearson")`. (`sd` is `NA` under the package
  default too: only `method = "oneStepGaussian"` ever returns one.)

  This is also what resolves the negative composition `predicted` values
  reported as issue #108 point 1: 404 of 4538 rows on BS2017SS, minimum
  -10.86, each with a positive-biased residual. Under `"cdf"` no expected
  count is formed at all. The Gaussian methods still report and warn about them,
  unchanged; WHAM does the same thing on its own example, so this is not an
  Rceattle defect.

* **`Index_distribution = "TruncatedNormal"` is exact under `"cdf"`, in the
  main call.** The template supplies the truncated CDF
  `[Phi(z) - Phi(-m)] / Phi(m)` in closed form, so the family no longer needs
  its own `"oneStepGeneric"` call over `(0, Inf)`. That removes the numerical
  integration, its cost, its failure to converge on some random-effects models
  (which fell back to an approximate spline), and its side effect of moving the
  other fleets' residuals.

* **Known limitation: compositions at scale under random effects.** On a
  random-effects model with a large composition data set, `method = "cdf"`
  returns non-finite residuals in bulk and is very slow. On `BS2017SS` with
  `random_rec = TRUE` (159 random effects, 4538 composition bins), **1879 of
  4538 residuals are non-finite** against 0 for `"oneStepGaussianOffMode"` on
  the same fit. The failures are a contiguous tail and the same rows come back
  clean when residualized on their own, so it is the depth of the conditioning,
  not the observations; redoing the tail on a fresh call recovers nothing.
  What binds is that depth, not the presence of random effects: the same
  method residualizes 1680 composition bins on a 22-random-effect model
  correctly. **So try `"cdf"`, and when the warning reports non-finite residuals
  in bulk, fall back to a Gaussian `method` for that source**, taking its
  composition residuals as under-dispersed by about a factor of two and a half
  (the numbers are in the bullet below). `"cdf"` is sound on fixed-effect
  models, and on random-effects models for the aggregate and covariate series.

* **A Dirichlet-multinomial composition cannot use `"cdf"`** and is residualized
  with `"oneStepGaussianOffMode"` instead, announced in a message and recorded
  as `attr(osa, "method")["DirichletMultinomial"]`. The conditional is a
  beta-binomial, which has no elementary CDF; the reference implementation TMB
  ships (`contrib/OSA_multivariate_dists-main/distr.hpp`) sums the pmf over
  `0..floor(x)`, which is a step function of a fractional count and costs `O(x)`
  beta functions per bin. This is said rather than approximated quietly,
  because a missing CDF term does not fail loudly. It makes both tails equal,
  giving `Fx = 0.5` and a residual of exactly 0 for every bin.

* **`|residual|` is censored at 8.04 under `"cdf"`, in both directions, and
  `osa_residuals()` warns when any residual sits there.** The upper end is
  forced: `Fx` is recovered from `1 / (1 + exp(nlcdf.lower - nlcdf.upper))` in
  double precision, which saturates at the last double below one, so nothing
  reading a CDF can report past 8.21 on that side. The lower end is *not*
  forced: that expression takes a small `F` down to about 1e-308, a residual of
  -37. It is censored to match anyway, because an asymmetric ceiling would show as a
  long left tail against a wall on the right, which is what skewness in the
  residuals looks like. The cost is real: `osa_diagnostics()` computes SDNR and
  the tail statistics on the censored values, and it bites hardest on a short
  series where one observation drives the statistic. Which method to reach for
  then is measured rather than asserted
  (`tools/verify/verify-osa-cdf-accuracy.R`, a 12-year survey with one
  observation multiplied by 200): `"oneStepGaussian"` reports it uncensored at
  38.98 (SDNR 12.89), `"cdf"` censors to 8.04 (SDNR 4.59), `"oneStepGeneric"`
  compresses it to 3.33 (SDNR 2.21), and the package default
  `"oneStepGaussianOffMode"` returns `NaN` there, making its SDNR unusable
  rather than merely large. That last failure is magnitude-dependent: at a x20
  outlier the default is finite and matches `"oneStepGaussian"`, so **reach
  for `"oneStepGaussian"` specifically**, not for "a Gaussian method", and on
  the fleet in question rather than a whole composition source, since it costs
  an `nlminb` per observation. The template shrinks the CDF away from 0 and 1 by
  four machine epsilons rather than one, so the ceiling does not land on that
  recovery's rounding tie, which produced 6 infinite residuals on BS2017SS.

* **`"cdf"` is exact only without random effects, and the recommendation splits
  by data type there.** With random effects [TMB::oneStepPredict()] integrates
  the CDF over the latent states by Laplace, and that integrand is a Gaussian
  times a sigmoid rather than a density. Against the exact Kalman innovations of
  a linear-Gaussian state space model, `fullGaussian` and `oneStepGaussian` are
  exact to machine precision while `"cdf"` errs by 7e-4 to 4e-2 as the latent
  state becomes more informative. **That is a result about a LINEAR-Gaussian
  model and does not transfer wholesale**: the Gaussian methods are exact when
  the one-step-ahead predictive is Gaussian, which needs the model to be linear
  in the random effects, and Rceattle's index and catch are `exp()` of cumulated
  log recruitment deviations through the population dynamics. `fullGaussian` and
  `oneStepGaussian` cannot differ for a Gaussian conditional, and on a
  17-deviation fixture they differ by 0.091 on both index and catch, where
  `"cdf"` differs from `oneStepGaussian` by 0.017 on index, so **no method is
  exact for index or catch under random effects** and this release does not
  claim one. **`"ecov"` is the exception**: its conditional genuinely is
  linear-Gaussian, the two Gaussian methods agree there to 4e-14, and `"cdf"`
  sits 0.139 away, about a quarter of the residual sd, so prefer a Gaussian
  method for that source. It does not reverse for
  compositions, whose conditional is discrete and skewed, exactly what the
  Gaussian methods get wrong, and by much more. Simulating from a
  22-random-effect model with the recruitment deviations redrawn (1680
  residuals): `oneStepGaussianOffMode` gives mean +0.513, sd 0.404 and KS
  rejection in 120 of 120 replicates; `"cdf"` with `discrete = TRUE` gives mean
  +0.006, sd 1.002 and 6 of 120, the nominal 5%.

  Both halves of that split are what the SAM authors do across two packages.
  `stockassessment::residuals.sam()` takes the Gaussian default with
  `discrete = FALSE`, their composition package `compResidual::resMulti()`
  hardcodes `method = "cdf", discrete = TRUE`. Rceattle fits both kinds of data,
  so it chooses per observation type inside one call
  (`vignette("model-diagnostics")`).

* **Every [TMB::oneStepPredict()] call now conditions on the observations that
  precede its group** (`conditional =`), instead of discarding them. The
  composition and aggregate rows must be split into separate calls under
  `"cdf"` because they need different `discrete` settings, and without this the
  split zeroed the aggregate data terms while the compositions were
  residualized, which on a random-effects model moved the composition
  residuals by up to 0.99 on a 21-random-effect fixture. It also improves the
  pre-existing `"TruncatedNormal"` split. Fixed-effect models are unaffected,
  and `verify-refit-like.R` is bit-identical.

## Ease of use

* **`osa_residuals()` says when composition residuals are taken at the default
  method.** The default is biased on composition data by the scoring table in
  `?osa_residuals`, and stays the default because `"cdf"` returns non-finite
  residuals in bulk on a deeply nested random-effects model. Naming any method,
  the default included, is taken as a choice and stays silent.
* **The install commands in the README and `?print.Rceattle` work as written.**
  The version pin pointed at 4.3.0 and the tag convention was documented as
  `@vX.Y.Z`; releases since 5.0.0 are tagged bare. The example links pointed at
  `blob/master`, and there is no `master` branch.
* **`print()` on a `model_config` lists the fields it imposes.** `fit_mod(config = )`
  overlays only the fields the config set, so this is what a config will change
  on a data object, readable before fitting rather than from the warnings the
  fit raises.

# Rceattle 5.40.0

## New features

* **A non-parametric selectivity form whose deviations integrate.**
  `NonParametric` charges its shape penalties on each year's realized curve, so
  under `random_sel = TRUE` the density the Laplace approximation integrates is
  tilted and the reported deviation SD is not the SD of the deviations;
  `fit_mod()` refuses that combination, and the random-walk mode for its own
  reason (`inst/dev/TODO-selectivity.md`). `Selectivity =
  "NonParametricIntegrable"` (code 13) keeps `NonParametric` and
  `NonParametricPM` exactly as they are and adds the missing density: it is the
  Ianelli base curve with the decreasing, curvature and average-selectivity
  penalties charged once on the base coefficients, and the deviations scored by
  `dnorm(0, sel_dev_sd)` on the estimated bins, so `random_sel = TRUE` estimates
  the SD from a complete density.

  `Time_varying_sel` picks the structure, as it does for `NonParametric`:
  `"Off"` estimates no deviations, `"IID"` gives independent annual deviations,
  and `"RandomWalk"` gives increments from the base curve with the start-year
  increment fixed at 0. Under `"Off"` the form reproduces the `NonParametric`
  objective to the last digit, at any `Bin_first_selected`. It reads no
  `Sel_cap_bin`, and `fit_mod()`'s refusals name it as the alternative. The
  per-year mean of a year's deviations is removed by the curve's centring, so
  the data never see it: under `random_sel = TRUE` it integrates out exactly,
  under `random_sel = FALSE` those directions are pure prior.
* **A non-parametric fleet's coefficients below `Bin_first_selected` are held at
  0.** They are mapped off, but the curve reads them: each year is centred by the
  log mean over every bin, so a value there shifted the whole curve while no
  density scored it. `inits` from a fit with a lower `Bin_first_selected` hold
  such values: 0.9 in those cells moved `Atka2022`'s fishery objective by 704
  nats and year-1 selectivity by 0.21. This affected `NonParametric` (2) and the
  new form 13. `NonParametricPM` (9) was never affected, because its
  branch already zeroed those cells before centring, so no `NonParametricPM` fit
  needs revisiting. A fit started from the build defaults, the golden fits
  included, is unchanged.
* **What the estimated SD is worth.** `tools/verify/verify-sim-recovery-np-integrable.R`
  draws deviates at a known SD on `Atka2022`'s fishery (multinomial age
  compositions, input sample sizes 2 to 236), simulates the observations and
  refits with `random_sel = TRUE`. The estimate is biased low: 0.24 to 0.26
  against a true 0.35 across two runs of eight and ten replicates, every
  replicate below the truth (reported SE on the log scale 0.12), 13% low at
  0.70, and 9% low at 0.35 with the sample sizes multiplied by ten. Every scored cell has its density and no estimated cell
  is unscored (checked cell by cell), so this is the Laplace marginal
  likelihood's known downward bias for variance components on small
  multinomial samples (Breslow and Lin 1995), not a scoring defect. Read the
  reported SD as a lower bound on the process SD unless the compositions are
  well sampled; MCMC through `tmbstan` on the fitted object removes the
  approximation.

# Rceattle 5.39.0

## Stock-recruit curves under predation

* **Alpha and beta are bounded at +/-30 on the log scale.** Under predation
  nothing ties a stock-recruit curve to an equilibrium, and on a flat ridge
  log alpha has reached 702, next to the double-precision limit; +/-30 is 13
  orders of magnitude either side of any stock's scale, so it never binds a
  determined estimate. A linkage bound on alpha or beta overrides it.
  Predation is the reason the bound was added, but the bound itself is
  unconditional: `build_bounds()` applies it to every model, single-species
  included, and the refit stop below reaches any fit whose `rec_pars` sit
  outside it, estimated or fixed.
  **A fit saved on the old unbounded ridge, with log alpha or log beta beyond
  +/-30, no longer refits**: `build_bounds()` stops because its starting values
  are outside the bounds, which takes `retrospective()`, `profile()` and
  `run_mse()` on that fit with it. The error names the block and says so.
  Restart such a fit from values on the stock's scale
  (`build_srr(srr_alpha_init =, srr_beta_init =)`); the curve it came from was
  the degenerate one this release exists to report.
* **The convergence battery now reads the curve under predation instead of
  skipping it.** Steepness needs spawning biomass per recruit, which is
  undefined under `msmMode > 0`, so `check_convergence()` returned a NOTE
  without looking at the curve, and the flat Pacific hake curve passed. It now
  reports a WARN when the curve is flat over the observed SSB range
  (Beverton-Holt predicted/asymptote above 0.9 at the lowest SSB), linear (below
  0.1 at the highest; for Ricker, a density-dependence factor above 0.9), when a
  Ricker peaks below the lowest observed SSB, when log alpha or log beta sits at
  the overflow bound, or (with `getsd = TRUE`) when a log-scale standard error
  exceeds 10; the record holds alpha, beta, their standard errors and the
  density dependence at both ends of the SSB range per species. A curve held
  at its inputs is a NOTE, not a warning. The hindcast standard errors it reads
  are now kept in the fit's convergence snapshot, since under an estimating HCR
  `fit$sdrep` is the projection's.
* **Simulation recovery.** `tools/verify/verify-sim-recovery-srr-msm.R`
  imposes a Beverton-Holt curve with its bend inside the observed SSB range on
  the two-species fixture (fished with a strong pulse, SSB spanning about 7x)
  and refits it from the true values after redrawing the observations and the
  recruitment process; a `dispersed` option starts a log unit away. Over 30
  replicates from the truth, species 2 recovers cleanly (log alpha 2.51
  against a true 2.55, empirical SD 0.19, mean reported SE 0.20; log beta
  -5.64 against -5.66); species 1, whose recruitment is noisier, recovers in
  the replicates that converge and runs to the ridge in those that do not.
  From dispersed starts species 2 still recovers (2.58 against 2.54) and
  species 1 reaches the ridge in half the replicates, which is the behaviour
  the new check reports.
* Tests: `sample_rec()` projects the mean multiplier at both draw sites
  (`test-functions-sample-rec-agreement.R`); one fitted Ianelli case under
  predation; the degenerate-curve check on both ridges and both forms.
  `?build_srr` says what the `R0` slot starts at under a multispecies curve.

# Rceattle 5.38.0

## New features

* **A per-sex apical selectivity offset.** Fishing mortality is one `log_F`
  per fleet and year shared by the sexes, so a sex difference in F can only
  come from selectivity, and no form had a height parameter: the logistic
  family and `DoubleNormal` peak at 1 for every sex, the non-parametric forms
  re-centre each sex, `Hake` normalizes each sex by its own maximum. The new
  selectivity linkage parameter `apical` (base block `log_sel_apical`,
  `[n_fleets, nsex]`, log scale) multiplies one sex's whole curve by
  `exp(log_sel_apical)`, applied after the form and before the shared
  normalizer, so it works for every estimated form. Name the fleet and the sex
  that holds it, as Stock Synthesis's male-offset option does:
  `build_selectivity(linkages = list(apical = linkage_spec(~ 1, by = ~ fleet + sex, fleet = 3, sex = "male", priors = list(intercept = lognormal(0, 0.5)))))`.
  The multiplier equals the ratio of the sexes' peak heights only where their
  shapes peak equally (the logistic family on an age axis); for a dome with
  sex-specific shape it is the pointwise multiplier on that sex's curve. Only
  the male:female contrast is identified (the common level is `log_F`), so
  `fit_mod()` refuses an `apical` linkage that names no fleet, no sex or both
  sexes, one on a one-sex species, a `Fixed` or AR1 fleet, a mirror fleet, one
  with `link = "identity"`, and one on a fleet whose
  `Sel_norm_scope = "WithinSex"` normalization would divide the offset out.
  The contrast is informed only by joint compositions (`comp_data$Sex = 3`), and
  `fit_mod()` warns when the named fleet has none and the offset has no
  prior, since nothing then informs the ratio. Naming `fleet` and `sex` is
  enough: `by` defaults to `~ fleet + sex` for this parameter.
  The multiplier is compiled into the template only when a linkage names it,
  so every existing model keeps its AD tape and the golden fits are unchanged;
  `inits` and a stored `map` from an older fit are filled with the block fixed
  at 0, so saved fits still refit. `tools/verify/verify-sim-recovery-apical.R`
  is the simulation-recovery harness: on `GOAatf`'s fishery with a true log
  offset of 0.693, 120 replicates give a mean of 0.685 (empirical SD 0.28,
  mean reported SE 0.35), so the estimate is unbiased and its standard error
  conservative.

# Rceattle 5.37.0

## Breaking changes

* **The dead QAR1 catchability path is removed.** `Catchability = "AR1"` has
  been refused since 5.12.0 and the live form is a q linkage
  (`ar1(1 | Year)` with `observe`), but two template blocks keyed on the
  retired code stayed: one that overwrote `index_q` (discarding the linkage
  offsets) and its AR1 density, plus the `index_q_rho` parameter only they
  read. Both blocks and the parameter are gone; `index_q_rho` drops out of
  `parameter_dictionary()`, `set_phases()` and the map. An older fit's `inits`
  and stored `map` naming it are accepted (the block is dropped as retired),
  and a stored `map` sizing `log_pop_scalar` by age (before 5.35.0) is
  collapsed as `inits` already were, levels included, so a map that estimated
  an age-specific scalar collapses to one per species rather than stopping in
  TMB, so `retrospective()`, `profile()` and `run_mse()` on a saved fit still
  run. A `map` name the model has no parameter for is dropped with a warning
  (a retired block is dropped silently), so a misspelling no longer fixes
  nothing; a `map` missing a parameter the model has stops with a message
  naming it, where TMB used to fail on the template read; a map level no cell
  holds is dropped rather than becoming a parameter without a start value;
  and a `map` that estimated an age-specific `log_pop_scalar` warns that only
  the first age is kept. A script that sets `inits$index_q_rho` keeps running,
  with that assignment now inert. No reachable fit changes; the golden fits
  are unchanged.

## Documentation

* **A contributor path.** `CONTRIBUTING.md` (setup, tests, branches, what a
  pull request owes) replaces the branch table and commit convention in the
  developer guide. A new article, *Adding a selectivity form*, traces
  `Selectivity = "DoubleNormal"` through every file, fits it on `GOApollock`,
  and shows what the drift guards report on a half-finished form;
  `test-docs-anchors.R` checks every path, function and code it quotes. The
  C++ template is now published as a Doxygen reference from the site's
  Contributing menu (`Doxyfile`, built by `pkgdown.yaml`); the growth header's
  equations render, and stale `@param` names in the growth, selectivity and
  predation headers are corrected. The forms table no longer calls the double
  normal six-parameter: it estimates four (peak, two widths, right-tail floor).
  Two defects the trace found are recorded in `inst/dev/CLEANUP_BACKLOG.md`,
  not fixed here: a `DoubleNormal` fleet with `Time_varying_sel =
  "RandomWalkAscending"` silently drops its deviates, and the form's default
  starting values describe a flat curve the optimizer does not leave.

# Rceattle 5.36.0

## Results change

* **`fit_mod(config =)` overlays only the fields the config set.** It replaced
  `data_list$model_config` wholesale, so a config from `model_config()` or
  `load_config()` silently dropped every linkage attached with
  `build_data(model_config = )` (57 random effects to 0 on a GOAatf survey
  `rw(1 | Year)` fit). Now only the fields the config set, defaults included,
  replace the data's, with a warning where they differ (a `build_*()` field
  compared as `save_config()` writes it). A config saved from a
  fit sets every field, so it still reproduces that fit. A script that relied
  on the wholesale reset must name every field it means to set.

## Bug fixes

* **A `species =`, `sex =` or `fleet =` filter on a linkage spec now warns
  when it does nothing**: when `by` does not include the term, and when it
  matches none of the model's levels (the whole spec was dropped, so a male M
  prior on a one-sex species vanished silently). Refits stay quiet.
* **OSA outlier symbols no longer scale with the panel's size.** A fixed
  `|resid| > 3` expected 13.5 flags in a 5,000-residual panel and 0.3 in a
  100-residual one. The OSA panel now flags above `qnorm(1 - 0.05 / (2 n))`
  (Bonferroni per panel), so under the model a panel expects 0.05 flags
  whatever its size. The Pearson panel keeps 3: those residuals are not N(0, 1).

## Documentation

* `fit_mod()` and `run_config()` say what `random_sel` and `random_q` gate: the
  `Time_varying_*` deviations only, with one estimated sd per
  `Selectivity_index` group; linkage random effects integrate regardless.
* `build_srr()`'s starting values name the `exp(9)` mean recruitment start.
* `build_osa_data()`'s note on `comp_offset` names its three fill sites.

# Rceattle 5.35.0

## Results change

* **The plus group's mean length is weighted by the species' own natural
  mortality.** Under parametric growth (`growth_model > 0`) the ages pooled in
  the plus group were weighted by survival at a hard-coded M = 0.2 whatever the
  species' M; they now use the base M1 at the oldest age (fixed or estimated,
  without year-varying offsets). This is still survival at M, not at total
  mortality: fishing and predation are excluded, so a heavily fished or preyed
  plus group is still somewhat too long. With `M1_model > 0` the plus group's
  length and weight now depend on the estimated M1. Length, weight and spawning biomass of the
  plus group move for every parametric-growth fit whose M is not 0.2: on the
  2024 GOA Pacific cod model at its starting values (M = 0.49) the plus-group
  mean length goes from 99.89 to 98.80 cm, its weight from 12.28 to 11.48 kg,
  and first-year female spawning biomass down 0.66%; refitted, its objective
  goes from 6416.12 to 6415.82 and terminal female spawning biomass from
  48,360 to 48,219 t (-0.29%). Empirical weight-at-age fits, and so the golden
  references, are unchanged.

## Breaking changes

None of the refused configurations below fitted the model it described, so this
release stays a minor version.

* **`estDynamics = 3` (`"FixedScaledByAge"`) is retired.** It never estimated
  an age-specific multiplier on the input numbers-at-age: `build_map()` freed
  only the first age's scalar under predation and none in single-species mode,
  so every such fit was the `estDynamics = 2` model. Measured on the
  fixed-numbers fixture, codes 2 and 3 both gave objective 222899009.1243875
  with one free scalar. `switch_check()` now refuses the code and names `2`;
  rebuild a stored fit with `estDynamics = 2`, unchanged. `log_pop_scalar` and
  the reported `pop_scalar` are one value per species: a script reading
  `quantities$pop_scalar[sp, 1]` now reads `pop_scalar[sp]` (the GOA CEATTLE
  figure scripts are updated), and `fit_mod()` collapses an older fit's matrix
  in `inits` to its first column. The multiplier is estimated only under
  predation (`msmMode > 0`); in single-species mode nothing informs it, so `2`
  fits as `1`, which the schema, `?BS2017SS` and `data_check()` now say.
* **A recruitment linkage on a species with input numbers-at-age
  (`estDynamics > 0`) is refused.** Its recruitment is read from `NByageFixed`
  and `rec_pars` is fixed, so the linkage fitted nothing: on the fixed-numbers
  fixture an intercept prior added 486.98 nats to the objective (and so to
  AIC) with no free parameter, and a covariate slope was a free parameter with
  a zero gradient, which leaves the Hessian singular. A spec with no
  `species =` expands to one row per species, so it is refused too; the message
  names the estimated species to put in `species = c(...)`.
* **`random_sel = TRUE` is refused for `Selectivity = "NonParametric"` with
  `Time_varying_sel = "IID"` at every `Sel_curve_pen` setting.** Until now
  `Sel_curve_pen1 = 0` lifted the refusal, but the average-selectivity penalty
  is always charged on each year's realized curve and does not scale with the
  deviation sd, so the sd that fit reported was the sd of a tilted density
  (about 5% of the precision low at sd 0.35, more as the sd grows). Fit with
  `random_sel = FALSE`, the penalized AMAK formulation. (From 5.40.0 the
  integrable form `NonParametricIntegrable` takes
  `random_sel = TRUE` instead; this refusal does not apply to it.)

## Bug fixes

* **A non-positive stock-recruit curve no longer gives a NaN objective.** An
  identity-link offset on alpha, beta or R0 can drive the curve to or below
  zero; on the single-species fixture an alpha offset of -100 per unit covariate
  gave `R_hat` of -89 and a NaN objective, stock-recruit penalty and dynamic B0
  under Beverton-Holt, and a NaN `log(alpha * SPR0)` in `R_hat` and steepness
  under Ricker. When the model has an identity-link recruitment linkage,
  hindcast recruitment, R0, R_init, the penalty curve, `R_hat` (first year
  included) and the Ricker log arguments are kept positive by `posfun()` (a
  0.001 barrier: a badly negative curve returns a value well below 0.001, not
  one fish) with the excursion charged to the "Zero n-at-age penalty" row, and
  `fit_mod()` warns on the offset. Projected recruitment and the recruitment
  behind SB0 and dynamic B0 are not yet floored and can still be negative under
  such a linkage (`inst/dev/TODO-srr-multispecies.md` item 14); check those
  series. Without such a linkage the template is exactly as before, so no
  existing fit moves.
* **The "Zero n-at-age penalty" row was a running total across cells and
  species.** The accumulator behind the numbers-at-age floor (and the Ricker
  intercept floor) was reset once per iteration, so each cell added every
  earlier excursion in every species: with species 1's log R0 at -20 on
  BS2017SS, species 2's row read 1.4e-3 from a parameter it does not share, and
  species 1's own value was inflated by its cell count. Each floor now adds its
  own excursion only. Zero for every fit that never touches a floor (the golden
  fits are unchanged); a fit that does gets a smaller, per-species penalty.
* **`check_convergence()` reports a non-zero or non-finite "Zero n-at-age
  penalty" row**, naming the species and the root-sum-square excursion below
  the floor (WARN under 1, FAIL above or when the row is not finite; thousands
  of fish for numbers and recruitment, unitless for the Ricker intercept): a fit
  whose numbers-at-age, Ricker intercept or recruitment sat on the 0.001 floor
  is not the model as specified, and nothing else showed it.
* **A factor switch fitted as its level index.** `.map_switch()` passed a
  factor through (`read.csv(stringsAsFactors = TRUE)`), and every downstream
  comparison against an integer code was FALSE rather than an error, so
  `srr_est_mode = factor("LognormalPrior")` fitted as code 1 with the prior
  dropped. Factors and numeric-looking strings now resolve to the code they
  name, and a numeric-looking string outside the map (`"99"`) is an error, as
  an unknown name is.
* **Editing `Time_varying_sel_sd` or `Time_varying_q_sd` and refitting from
  `inits` was a silent no-op.** The columns are read on a fresh build only;
  `fit_mod()` now warns, as it does for `Comp_weights`, unless `random_sel` /
  `random_q` estimates the sd.
* **A species with input numbers-at-age reports its input recruits as `R`,
  and `NA` for its stock-recruit quantities.** Its `rec_pars` are fixed, so
  `R`, `R0`, `R_init`, `avg_R`, `steepness` and `SPR0` held the
  `build_params()` placeholder (`R0 = exp(9)`, a flat 8103 that neither the
  model nor the user set); in single-species mode the equilibrium `SB0` and
  `B0` were built on it, and with `DynamicHCR = FALSE` the depletions divided
  by it. `R` is now the first-age input numbers (`NByageFixed` times
  `pop_scalar`, sexes summed) in the years `NByageFixed` covers and `NA` in
  any other (a year with no row holds zeros), with no confidence band; the rest
  are `NA`, and `plot_stock_recruit()` draws its points but no curve. `sample_rec()` and
  `retrospective()` skip such a species (its `rec_dev` is mapped out). Under `DynamicHCR = TRUE` the depletions are the input numbers
  relative to themselves and stay reported; under predation `MSSB0` replaces
  `SB0`.
* **`data_check()` says when `Sel_norm_scope` is not read.** On a two-sex
  `Hake` or `LogisticPM` fleet the column changed nothing, without saying so
  (measured identical to every digit on GOAatf fleet 3). Neither form goes
  through the shared normalizer, for different reasons: `Hake` normalizes each
  sex to its own maximum in its own block, while `LogisticPM` does not normalize
  at all and reuses `Sel_norm_bin` as a penalty bin range. The notice names the
  reason that applies, and the schema and the vignette say it too.
* `retrospective()`'s note on a penalty-form peel with no penalty years now
  covers a peel that keeps a single year.

# Rceattle 5.34.0

## Results change

* **Multispecies `ConstantF` projects at the input F.** The multispecies
  projection loop in `fit_mod()` reset every species' `log_Ftarget` to 0, so a
  `ConstantF` run under predation projected each species at F = 1 whatever
  `Ftarget` was, including a `ConstantF` estimation model refit inside
  `run_mse()`. Projected catch, SSB and depletion change for those runs; with
  `HCRorder > 1`, so do later species' multispecies SB0. A single `Ftarget` now
  recycles to every species, as `?build_hcr` documents (with more than one
  species it stopped with a map-size error), and a missing or negative
  `Ftarget` under `ConstantF` is an error. An `Ftarget` of 0 is stored as
  log F = -999, the value `build_params()` gives a fleet with no catch, not as
  -Inf; the projection was already unfished, only the stored parameter changes.
  A single-species projection under a rule that estimates `Ftarget` starts it
  at log F = 0 when `inits` hold that no-fishing value, as the multispecies
  loop already did; before, inits from a `ConstantF` fit with `Ftarget = 0`
  left the rule's `Ftarget` at 0, where its gradient is exactly 0.

* **Projected `F_flt` and `F_flt_age` are indexed by fleet, fisheries only.**
  The projection wrote `F_flt` by species index, so a fishery's row reported
  another species' F, and survey rows took `NA` times F or a full share of F.
  Both are reported quantities only, so no fit changes. Hindcast rows were
  already correct.

* **A species with input numbers-at-age (`estDynamics > 0`) is projected at
  F = 0 and takes no harvest control rule.** `build_hcr_map()` already left
  its `log_Ftarget` / `log_Flimit` unestimated and the reference-point
  penalties already skipped it, but the projection still fished it at those
  start values (F = 1 under most rules). Its projected catch is now 0. Its
  numbers are input, but a fixed species with a fishery now has a higher
  within-year mean abundance in the projection, which changes predation
  mortality on and by it, and its own SSB changes when `spawn_month > 0`.

  `Ftarget`, `Flimit`, `SPRtarget`, `SPRlimit`, `SBF` and `DynamicSBF`, all set
  by that unestimated F, are `NA` in `fit$quantities`; `report_tables()` blanks
  the five it reports with that reason, and `plot_f()` draws no `Ftarget` or
  `Flimit` line for the species. `mse_summary()`
  reports the species' `P(Fy > Flimit)` metrics, and the `P(SSB < SSBlimit)`
  metrics that read `SBF`, as `NA`; a fixed species with a fishery reports
  catch IAV and P(Closed) as `NA`, like an unfished one. Its depletion is still
  reported.

  In the four-species hake model this is arrowtooth, sablefish and California
  sea lions; none has a fishery and `spawn_month` is 0, so its dynamics and
  catches do not change, but its summary reports `NA` where it reported 0.

* **`mse_summary()` takes multispecies dynamic depletion from the operating
  model's `DynamicSB0`**, as it already did for single-species models.
  `om_terminal_dynamic_sb0` and `om_terminal_depletion_dynamic` were read from
  the no-fishing refit (`OM_no_F`). Dynamic SB0 is the OM's own history with no
  fishing: the stock-recruit curve, the realized recruitment deviations and the
  predation suitability fitted in the hindcast. On the Pacific hake MSE
  (`MSE_yr2024.R`, Beverton-Holt operating model, two simulations) hake's
  terminal dynamic depletion goes from 0.99 to 0.78.

* **`remove_F()` sets F to 0 from the year after `endyr` by default**, not from
  the year after the latest `suit_endyr`. When the suitability window ended
  before `endyr` it removed fishing inside the hindcast: on the hake MSE,
  2020–2023 of a 2023 hindcast. A new `styr` argument gives the first year
  fished at F = 0; under predation it must fall after the empirical-suitability
  window (the `suit_endyr` of every predator with `suitMode = 0` and non-zero
  fitted suitability), since removing fishing inside it would change the
  suitability the model was fit with. As before, the
  projection is unfished whatever harvest control rule the model was fit under,
  so `styr` can be no later than the year after `endyr`.
  `run_mse()` now builds `OM_no_F` with no fishing
  after the original operating model's terminal year, so it matches the OM
  through that year. The `OM no F: SSB Collapse` and `OM: SSB Collapse from F`
  metrics change for runs whose suitability window ended before `endyr`; an MSE
  saved under 5.33.0 keeps its old `OM_no_F`, so rerun `run_mse()` to update
  them.

# Rceattle 5.33.0

## Breaking changes

None of the refused configurations below fitted the model it described, so this
release stays a minor version. Stored fits with them no longer refit.

* **`data_check()` refuses an `R0` linkage under a stock-recruit curve fitted
  in the hindcast.** Beverton-Holt and Ricker recruitment read only alpha, beta
  and SSB, so in a single-species model `R0` is derived from them and the
  linkage did not change recruitment. In a multispecies model `R0` sets only
  the initial recruitment level (`R_init`), so an intercept-only linkage
  (`~ 1`, a prior or fixed value) is accepted and a covariate or random effect
  is refused. Put an environmental effect on productivity through an `alpha`
  linkage, or use `srr_fun = "mean"` with the `R0` linkage.
* **A stored fit with `srr_est_mode = "BetaPrior"` on a Ricker curve does not
  refit**; `build_srr()` has refused the pair since 5.30.0. The prior is on
  Beverton-Holt steepness and never applied to a Ricker curve. Rebuild with
  `srr_est_mode = "Estimated"`, the model it fitted.
* **`srr_est_mode = "Fixed"` together with an alpha linkage fixed at its init
  (`est_phase = 0`) is refused.** Each fixed alpha at a different value. Keep
  the linkage alone.
* **A stored fit with stock-recruit penalty years outside the hindcast no
  longer refits** (see Bug fixes). Rebuild with the model's full `build_srr()`
  call, with the penalty years inside `styr` to `endyr`:
  `build_srr(srr_hat_endyr = )` on its own resets the curve and its priors to
  the defaults (mean recruitment, `srr_pred_fun = 0`), which drops the penalty.
  `retrospective()` still runs when only `srr_hat_endyr` is past `endyr`,
  truncating it to each peel.
* **A stored fit with a Ricker `srr_est_mode = "LognormalPrior"` and an alpha
  linkage intercept prior no longer refits** (see Bug fixes): one density
  counted twice. Keep one.
* **`Catchability = "Estimated-with-prior"` with a prior on that fleet's q
  linkage intercept is refused**: both are priors on the same q. Set
  `Catchability = "Estimated"`, or drop the intercept prior.
* **`M1_use_prior = TRUE` with a prior on that species' M1 linkage intercept is
  refused** unless `M2_use_prior = TRUE`, which puts the M prior on total M.
  Both are priors on the same M1. Set `M1_use_prior = FALSE`, or drop the
  intercept prior.
* **`minage = 0` with a Beverton-Holt or Ricker curve is refused.** The curve
  read that year's spawning biomass before the model computed it, so it gave
  zero recruitment in the hindcast (`srr_fun`) or in `SB0`, the reference points
  and a `proj_mean_rec = FALSE` projection (`srr_pred_fun`, the Ianelli form
  included). Use mean recruitment.

## Deprecated

* **Alpha priors, fixed values and starting values outside a linkage.** Each
  still works and warns once from `build_srr()`; refits are silent.
  - `srr_est_mode = "Fixed"`: use ``linkages = list(alpha = linkage_spec(~ 1,
    est_phase = 0, init = list(`(Intercept)` = <alpha>)))``.
  - `srr_est_mode = "LognormalPrior"` on a Ricker curve: use an alpha linkage
    prior, which gives the same objective, with `init` at the prior mean for the
    same start.
  - `srr_prior` as alpha's starting value (supplied, not the default 4): use
    `srr_alpha_init` or `init` on an alpha linkage.

  `srr_prior` stays for the Beverton-Holt steepness prior. `Bmsy_lim` is
  unchanged.

## Results change

* **Lognormal priors and the Ianelli stock-recruit penalty are mean-centred
  under `bias_adjust_proc`.** With `fit_control(bias_adjust_proc = TRUE)`, the
  default, each is centred at `-sd^2/2` on the log scale, so the value given is
  the prior mean of the natural-scale quantity. The Beverton-Holt steepness
  prior and the recruitment deviations already worked this way. With `FALSE`
  the value is the median, as before, but the recruitment deviations lose their
  centring too, so a fit with both a lognormal prior and recruitment deviations
  cannot be reproduced as it was by this flag. Reproduce an old prior with the
  flag on by shifting its input: `prior_lognormal(p1 + s^2/2, s)`,
  `Catchability_init * exp(s^2/2)`, a Ricker `srr_prior * exp(s^2/2)`, or
  `M_prior * exp(s^2)` for the old M centre (`inst/dev/TRAPS.md`). The
  Ianelli penalty with centred deviations has no input to shift. A value
  between 0 and 1 gives priors a centre that is neither the mean nor the
  median. This applies to
  `prior_lognormal()` in every linkage (intercepts, slopes and random-effect
  SDs), the Ricker alpha prior (`srr_est_mode = "LognormalPrior"`) and the
  catchability prior (`Catchability = "Estimated-with-prior"`). The penalty's
  curve term now treats the curve as the mean of recruitment, as in Ianelli's
  EBS pollock model (`pm.tpl`, Dorn 2002); `pm.tpl` penalizes its recruitment
  deviations without a sigma, which no setting reproduces.
  `bias_adjust_proc = FALSE` gives the AMAK form, where the curve is the
  median. Every fit that uses these features with the flag on changes. With
  the flag off, only fits with an M prior (below) and penalty-form projections
  (next item) change. Wide priors change most: for `prior_lognormal(0, 2)` the
  mean falls from 7.39 to 1 and the median from 1 to 0.14.

  In the penalty years the Ianelli form penalizes recruitment twice, as a
  deviation from `R0` and as a departure from the curve. With the flag on,
  centring the penalty lowers the recruitment the two favour by
  `exp(-sigma_R^2/4)`: 8.6% at `sigma_R = 0.6`. Their mean was the geometric
  mean of `R0` and the curve; it and the median fall alike. Recruitment the
  data inform moves less.

  Measured changes:
  - **Pacific hake MSEs.** Their Dirichlet-multinomial weights have
    `prior_lognormal(0, 2)` priors, and they use an M prior. Objectives rise
    by 3.0 to 6.8, and hake terminal SSB changes by at most 0.86%. The survey
    DM weight falls from 43-49 to 31-35 (`MSE_yr2024.R`) and from 32-35 to
    23-25 (`04-mse.R`), and the Ianelli Beverton-Holt alpha, fitted then with
    no alpha prior, rises from 18.8 to 41.0. Every fit still has a
    positive-definite Hessian.
  - **GOA2018SS**, which has a q prior (SD 0.1). The golden objectives change
    by -0.015 (single-species) and -0.003 (multispecies).
  - **BSAI northern rock sole bridging** (a q prior; M priors and a Ricker
    penalty in the later models). At the 5.32.1 estimates, the first model's
    objective changes by 5.7, all from the q prior. None of the refits
    converges (maximum gradient 1.6 to 11.5), and they end at different
    parameter values, with objectives from 95 lower to 1,779 higher. Compare
    them only once they converge.
  - **Unchanged:** the Bering Sea golden models (to 3e-11), and any
    mean-recruitment fit run with `bias_adjust_proc = FALSE` and no M prior,
    such as the GOA pollock 2025 assessment (not refit).

* For a single-species curve fitted in the hindcast, `sample_rec(sample_rec =
  FALSE)` and `retrospective()` now set the projected recruitment deviation to
  the hindcast's mean deviation from the curve, `log(mean(exp(rec_dev)))`,
  which equals the multispecies form's `log(mean(R / R_hat))`. They used
  `log(mean R) - log(R0)`, which applied the curve's SSB effect twice: projected
  recruitment was off by the hindcast-average R/R0 the curve implies, 8.6% low
  for a Beverton-Holt curve at steepness 0.8 and a hindcast averaging 40% of
  SB0 (27% at steepness 0.5), and often high for a Ricker curve. This changes
  an operating model under `run_mse(sample_rec = FALSE)`, retrospective
  forecasts, and dynamic reference points in projection years (and projected F
  under `DynamicHCR = TRUE`). In the model `sample_rec()` returns, projected
  recruitment changes only with `proj_mean_rec = FALSE`. `run_mse()`'s
  default, `sample_rec = TRUE`, is unchanged.

* `sample_rec(sample_rec = FALSE)` and `retrospective()` now set the projected
  recruitment deviation for the Ianelli penalty form (`srr_fun = 0` with a
  curve in `srr_pred_fun`) to `log(mean(R / R_hat))` over the penalty years
  (`srr_hat_styr` to `srr_hat_endyr`), the mean ratio of recruitment to the
  curve. They used the mean log ratio over every hindcast year, which projects
  the median recruitment around the curve. This applies whatever
  `bias_adjust_proc` is. The deviation sets recruitment in a retrospective's
  peeled years and, under `run_mse(sample_rec = FALSE)`, in the operating
  model's added years. In the model `sample_rec()` returns, the dynamic
  reference points (`DynamicB0`, `DynamicSB0`, `DynamicSBF`) always change;
  projected recruitment changes only with `proj_mean_rec = FALSE`; and
  projected F changes only under HCRs 5, 6 and 7, through projected SSB with
  `proj_mean_rec = FALSE` and through the dynamic reference points with
  `DynamicHCR = TRUE`. `sample_rec()` stops if no penalty years fall in the
  hindcast; a retrospective peel that ends before the penalty years uses its
  own years instead, with a warning.

## Bug fixes

* **The M prior is centred on `M_prior`.** It was centred at
  `log(M_prior) + sd^2/2` for either flag setting, so `M_prior` was neither its
  mean nor its median. It is now the mean under `bias_adjust_proc = TRUE` and
  the median under `FALSE`. Every fit with `M1_use_prior` or `M2_use_prior`
  changes, including fits with the flag off.

* **Dynamic B0 uses the hindcast's realized recruitment for cohorts spawned
  before the first model year** (`minage > 1`). Under a stock-recruit curve,
  including the Ianelli form, those years used the curve at the first year's
  spawning biomass. Under mean recruitment with `proj_mean_rec = FALSE` they
  used that year's R0, which differs only with a year-varying R0 linkage.
  Neither is the recruitment the hindcast used. `DynamicB0`, `DynamicSB0` and
  `DynamicSBF` change. With `DynamicHCR = TRUE` so do the projections, the
  HCR 3 target, and `biomass_depletion` and `ssb_depletion`. Equilibrium `SB0`
  and `SBF` do not change. Models with `minage = 1`, including every bundled
  dataset, are unchanged.

* **Under the Ianelli penalty form, dynamic B0 uses each year's realized
  deviation from the curve.** Before `srr_mse_switchyr` the hindcast is on mean
  recruitment, so `rec_dev` there is measured from R0, and dynamic B0 multiplied
  the curve by it. On the Pacific hake operating model the curve sat at about
  1.5 times R0 in those years, so no-fishing recruitment ran about 1.5 times too
  high in every year before the switch. It now takes `log R - log R_hat`, and
  with fishing near zero dynamic B0 equals hindcast biomass, while
  numbers-at-age stay above the model's 0.001 floor. Suitability is
  unchanged: the dynamic runs use the suitability fitted in the hindcast.
  `DynamicB0`, `DynamicSB0` and `DynamicSBF` change for every penalty-form model
  (`srr_fun = 0` with a curve in `srr_pred_fun`), and so does `mse_summary()`'s
  dynamic depletion for a single-species penalty-form operating model. With
  `DynamicHCR = TRUE` so do the projections, the HCR 3 target, and
  `biomass_depletion` and `ssb_depletion`. Models with a curve fitted in the
  hindcast, or on mean recruitment, keep their hindcast dynamic B0.

* **A linkage intercept fixed at its `init` (`est_phase = 0`) holds that value
  when `fit_mod()` is given `inits`.** The init was applied only to freshly
  built parameters, so a warm start kept the `inits` value and mapped it off:
  a refit at alpha = 500 warm-started from a fit at 1170 stayed at 1170. It
  also wins over `srr_alpha_init`. This applies to every process with an
  intercept linkage (recruitment, M, q, selectivity, growth). Refits through
  `retrospective()`, `profile()` and the other diagnostics keep their fitted,
  or profiled, value, and skip the `srr_prior`, `srr_alpha_init` and
  `srr_beta_init` starting values too; `profile()` over alpha or beta on a fit
  with those was flat. A hand-written profile loop through
  `fit_mod(inits = , map = )` over a linkage-fixed parameter is now held at the
  init: use `profile()`, or give the loop's linkage `est_phase = 1`. On a fit
  with `srr_alpha_init` or `srr_beta_init`, `jitter()` and `self_test()` now
  perturb alpha and beta (they were reset after jittering), and retrospective
  peels start from the fitted alpha.

* **Stock-recruit penalty years must lie in the hindcast.** `data_check()`
  refuses `srr_hat_styr < styr` and `srr_hat_endyr > endyr` under the Ianelli
  penalty. Past `endyr` the penalty applied to projected recruitment (BS2017SS,
  `endyr + 5`: objective +1,215), and `sample_rec()` averaged over different
  years than the template. A stored `srr_hat_endyr` is kept when `endyr` is
  lowered; Breaking changes gives the rebuild.

* **A Ricker `srr_est_mode = "LognormalPrior"` together with an alpha linkage
  intercept prior is refused.** They are the same density, which the objective
  counted twice.

* **The warning that a supplied `map` fixes both stock-recruit parameters no
  longer fires on a refit's own map**: not for a debug or projection map, nor
  for a species whose alpha or beta linkage estimates a covariate or random
  effect. An intercept-only linkage still counts as fixed.

* `data_check()` resolves a string `msmMode` (e.g. `"SingleSpecies"`) before
  its multispecies stock-recruit checks; `"SingleSpecies" > 0` is TRUE in R.

## Documentation

* `?build_srr`: a linkage on `R0` acts under mean recruitment only. Under a
  curve fitted in the hindcast a single-species `R0` linkage is refused; in a
  multispecies model only an intercept-only one is accepted, as the initial
  recruitment level. The "Starting values" section gave alpha's default
  start as `e^3`; it is `srr_prior` (default 4) wherever that is an alpha.
* `?linkage_spec`: an intercept's `init` is on the parameter's natural scale,
  a slope's on the link scale.

# Rceattle 5.32.1

## Documentation

* **`?build_srr` now leads with `linkages` for priors, fixed values and
  covariates on `R0`, alpha and beta.** An intercept-only `linkage_spec()` acts
  on the parameter itself: a prior (`prior_lognormal()`, `prior_normal()`), a
  fixed value (`init` with `est_phase = 0`), or one species (`species =`). For a
  Ricker curve the lognormal linkage prior on alpha gives the same objective as
  `srr_est_mode = "LognormalPrior"`. `srr_est_mode` and `srr_prior` remain for a
  Beverton-Holt steepness prior, the one a linkage cannot express. The old
  example labelled "with a prior on steepness" set `srr_est_mode =
  "Estimated"`, which applies no prior and turns `srr_prior = 0.8` into alpha's
  starting value; it is replaced. The linkages vignette's Recruitment section
  gains the same recipes.

* `?linkage_spec`: `init` is a named list, not a numeric vector. The old
  example, `c(...)`, is refused by the function.

# Rceattle 5.32.0

## Breaking changes

* **`srr_fun` / `srr_pred_fun` codes 1, 3 and 5, and `srr_indices`, are now
  errors.** These were the environment-driven stock-recruit forms. 4.4.0
  removed their environmental term from the template, but `build_srr()` kept
  accepting them with only a soft-deprecation warning. A model built with them
  was fitted without its covariate and reported nothing wrong: its objective and
  parameter count were identical to the non-environmental model's. The 4.4.0
  entry below, which says they "continue to work", was wrong.

  `build_srr()` now stops and names the linkage that replaces them:
  `srr_fun = 0` with `linkages = list(R0 = linkage_spec(~ <covariate>))` for
  code 1, and `srr_fun` 2 or 4 with an `alpha` linkage for 3 or 5.
  `srr_indices = k` meant the k-th `env_data` column after `Year`. A fit made
  with code 1, 3 or 5 still refits, through `retrospective()`, `run_mse()` and
  the other refitting diagnostics, as code 0, 2 or 4; that is the model it
  actually fitted. `run_mse()` shows the warning, but `retrospective()`,
  `jitter()`, `profile()` and `self_test()` suppress it. The one exception is
  `srr_fun = 0` with `srr_pred_fun = 1`, which scored the penalty around mean
  recruitment: it refits without that penalty, so its objective changes. A code
  stored in a data object's `model_config` is not checked and fits as 0, 2 or 4
  without a warning, and a stored `srr_indices` is ignored. Results produced with
  these codes from 4.4.0 through 5.31.0 have no environmental effect and should
  be refit with a linkage. No fit without them changes.

# Rceattle 5.31.0

## New features

* **A multispecies model can fit a stock-recruit curve in the hindcast.**
  Spawning biomass per recruit is undefined under predation, so the curve's
  initial recruitment level (`R0`), which seeds the first year and the initial
  age structure, is estimated rather than derived as `(alpha - 1/SPR)/beta`;
  every year with modelled spawners follows the curve. Any `initMode` is
  accepted; under the fished modes (3, 4) `R0` trades off against the initial
  F. The initial state matters: on the Pacific hake four-species model,
  `"FreeParams"` changed terminal SSB by -35% and `"NonEquilibrium"` by -25%
  against the mean-recruitment fit, while the hindcast curve itself came out
  flat (see `inst/dev/TODO-srr-multispecies.md`). Steepness is reported as 0,
  the Ricker positivity penalty on `alpha * SPR0 - 1` (which added to the
  objective at every evaluation) is skipped, and `sample_rec(sample_rec =
  FALSE)` and `retrospective()` scale the projection by the mean ratio of
  recruitment to the curve, since there is no unfished recruitment to scale
  against. Single-species fits
  are unchanged.

# Rceattle 5.30.0

## New features

* **A multispecies model can have a stock-recruit curve as a recruitment
  penalty.** `data_check()` refused every Beverton-Holt or Ricker curve under
  `msmMode > 0` (5.12.0), including the Ianelli configuration
  (`build_srr(srr_fun = "mean", srr_pred_fun = "BevertonHolt")`). That
  configuration never needed spawning biomass per recruit: the penalty, the
  dynamic reference points and -- with `proj_mean_rec = FALSE` -- the
  projection read alpha, beta and SSB only. It is now allowed, and so are
  estimating alpha and beta, fixing alpha (`srr_est_mode = "Fixed"`) and a
  Ricker prior on alpha. Under the default `proj_mean_rec = TRUE` the
  projection, and so the multispecies `MSSB0` derived from it, runs on mean
  recruitment while dynamic B0 follows the curve.

  Two uses still read spawning biomass per recruit, which is undefined under
  predation, and are still refused: a curve fitted in the hindcast
  (`srr_fun` Beverton-Holt or Ricker), whose initial recruitment derives from
  `SPRFinit`, and a prior on Beverton-Holt steepness (`srr_est_mode =
  "LognormalPrior"` or `"BetaPrior"`).

  Under predation the penalty curve's first-year `R_hat` is `R_init`, where it
  was `-Inf`; `sim_mod()` and `retrospective()` take `log(R_hat)` over that
  year. `steepness` is reported as 0, and `convergence_diagnostics()` records a
  NOTE rather than testing the curve against the replacement line. No earlier
  fit changes.

* **A species with input numbers-at-age (`estDynamics > 0`) has no
  stock-recruit prior, Bmsy penalty or curve penalty.** Its recruitment is a
  placeholder, so those terms scored nothing real, and with an estimated
  population scalar (`estDynamics` 2 or 3) the curve penalty pulled on that
  scalar. Only a model combining a stock-recruit curve with such a species
  changes. In a single-species model the population scalar is never estimated,
  so no estimate changes; only the objective changes, by a constant.

* `fit_mod()` warns when a supplied `map` fixes both stock-recruit parameters
  while the curve still shapes recruitment, as a map reused from a
  mean-recruitment fit does. `build_srr()` stops on `srr_est_mode =
  "BetaPrior"` with a Ricker curve, which has no Ricker form and was estimated
  with no prior, and warns when `Bmsy_lim` is given for a non-Ricker curve,
  where it is ignored.

## Bug fixes

* **Dynamic B0 holds a species with input numbers-at-age (`estDynamics > 0`)
  at those numbers.** The dynamic reference runs (`DynamicB0`, `DynamicSB0`,
  `DynamicSBF`) projected such a species from its first-year numbers on its
  placeholder recruitment, so a fixed predator collapsed in the unfished run
  and took its prey's predation mortality there with it. In the Pacific hake
  four-species model, with arrowtooth, sablefish and California sea lion fixed,
  arrowtooth dynamic SSB was 38.6 mt in 2023 against an input 43,490 mt. Hake
  dynamic SB0 in 2023 falls from 2.73e6 to 1.95e6 mt, and its depletion
  against it rises from 0.516 to 0.723. Only multispecies models with a
  fixed-dynamics species change: their dynamic reference quantities, and the
  objective when `DynamicHCR = TRUE` fits the projection against them.

# Rceattle 5.29.0

## Bug fixes

* **Composition Pearson residuals divide by the effective sample size the
  likelihood used, not the raw input N.** `residuals(type = "pearson")` always
  used the multinomial variance at `Sample_size`, for every composition source.
  The likelihood does not. Under a multinomial (and the `MultinomialAFSC`
  default) the weight column multiplies the log-likelihood, so it *is* an
  effective sample size -- the model's own simulator draws at
  `n_nom * comp_weights(flt)` (`ceattle.cpp:3758`). The reported residual was
  therefore `1/sqrt(w)` times the right one: too small on an upweighted fleet,
  too large on a downweighted one, and `reweight_comps()` routinely tunes either
  way. Under a Dirichlet-multinomial the proportions are additionally
  overdispersed by `(n + conc)/(1 + conc)`, so a DM fleet's residuals were
  inflated and read as systematic misfit that was not there.

  Age/length composition, conditional age-at-length and diet are each rebuilt on
  their own likelihood's scale, because all three differ: comp and CAAL form
  their Dirichlet-multinomial concentration on the offset-inflated total, while
  diet renormalizes over the prey bins **plus** its "other prey" balance before
  scaling by the stomach sample size. The composition proportion offset
  (`comp_offset`, default 1e-5) is carried onto that scale too, which moves the
  smallest bins by ~1e-4 relative even at weight 1.

  Each of the three reads its family and weight from the `data_list` element
  that actually holds them: `fleet_control$Comp_distribution` /
  `CAAL_distribution` for the two composition sources, and the per-species
  `Diet_distribution` for diet. `bioenergetics_control` is a workbook sheet
  name rather than a `data_list` element, and a diet lookup routed through one
  resolves to `NULL`, silently scoring every predator as a multinomial and then
  reading `Diet_comp_weights` -- a **log** under a Dirichlet-multinomial -- as a
  natural-scale multiplier on the stomach sample size. Where a
  Dirichlet-multinomial weight genuinely cannot be recovered, the row falls back
  to the multinomial variance and says so; it is not left on the
  Dirichlet-multinomial with a substituted concentration, which would inflate
  the variance by a factor that looks deliberate and is not.

  Each switch is also read through its deprecated spellings (`Comp_loglike`,
  `CAAL_loglike`, `Diet_loglike`). `switch_check()` upgrades those in place when
  a model is built, but a fit **saved** before a rename still carries the old
  name, and `residuals()` runs on the saved object -- so without the fallback an
  old fit would quietly resolve to the schema default rather than to the family
  it was actually fitted under.

  This moves Pearson residuals on essentially every fit. It does not move any
  likelihood, parameter estimate or reference point -- `residuals()`,
  `plot_comp()` and the Pearson panel of `plot.rceattle_osa()` are the only
  outputs affected. Index and catch residuals are untouched (they have no
  composition family), and `report_tables()` reads only those.
  `test-likelihood-pearson-effective-n.R` checks the variance by simulation --
  data drawn under the family the likelihood assumes must return residuals with
  `sd = 1` -- rather than by restating the formula.

* **The OSA tail statistics and their null intervals now come from one
  estimator.** `osa_diagnostics()` reported `quantile()`'s type-7 interpolated
  tail against a *simulated* null. Both sides are now the `r`-th order
  statistic, whose exact distribution is `Beta(r, n - r + 1)` on the probability
  scale, so the interval is closed-form and covers at its nominal rate (measured
  0.94-0.96 for n = 10 to 200). Mixing the two -- an interpolated statistic
  against a closed-form order-statistic null -- drops coverage to about 0.86
  near n = 50, non-monotonically in n; `test-diagnostics-osa-tail-nulls.R` pins
  the coverage so that pairing cannot come back.

  The reported tail is now genuinely the `r`-th order statistic rather than "the
  2.5% quantile", so `osa_diagnostics()` gains `lower_r` / `upper_r` and
  `lower_p` / `upper_p` naming which order statistic was used and the nominal
  probability it sits at. Short series are handled: `r` is clamped to `[1, n]`,
  without which `qbeta(p, n + 1, 0)` returns 1 and `qnorm()` returns `Inf`,
  making the upper tail check pass for every series with `n <= 19`.

* **A composition fleet switched off returns no Pearson residual, rather than
  zero.** `Comp_weights = 0` multiplies a multinomial fleet's log-likelihood by
  zero, so the fleet is not fit and has no effective sample size. The residual
  divided by it anyway, giving `sd = Inf` and a residual of exactly `0` in every
  bin -- `plot_comp()` drew a fleet the model never saw as a perfect fit. Those
  rows are now `NA`. A Dirichlet-multinomial reads the same column as a log, so
  `0` there is a weight of 1 and is unaffected.

* **Diet Pearson residuals group by `stomach_id` where the data carry it.** The
  grouping was read off `fit$data_list`, which is the pre-`rearrange_data()`
  list and never holds that element, so the predator/sex/age/year fallback ran
  even when the table had an id. Where a dataset holds more than one stomach per
  predator-age-year, that pooled them into one normalizing total; the C++ scores
  each stomach separately (section 13.2). No bundled dataset carries the column,
  so nothing in the package changes.

## Deprecated

* **`osa_diagnostics(nsim=, seed=)` are ignored** and warn when supplied. The
  tail null intervals are exact, so nothing in that function simulates. Note
  the OSA residuals themselves remain randomized-quantile residuals: pass a seed
  to `osa_residuals()` to control those.

## Plotting

* **OSA and Pearson residual bubbles are drawn on a fixed `[0, 6]` size scale**
  in both `plot.rceattle_osa()` and `plot_comp()`, so two figures can be
  compared by eye -- a free scale made a well-fitting fleet and a badly-fitting
  one look alike, and the two figures scaled the same residuals differently.
  Residuals beyond 6 are truncated with a warning naming their count and largest
  magnitude; the truncation precedes the limit because `scale_size_continuous()`
  sets out-of-bounds values to `NA` and drops them silently, which would hide
  the very points worth looking at. The warning does not list every offending
  value: dividing by the effective sample size scales residuals by `sqrt(w)`, so
  an upweighted fleet can put thousands of them past the cap at once. Bin axes
  take whole-number breaks while there are fewer than 20 bins.

* **The OSA Q-Q panels annotate the tail statistics** and their exact null
  intervals in the lower right, alongside SDNR in the upper left, following
  `afscOSA`. `plot.rceattle_osa()` gains `add_sdnr_ci` and `add_qq_quantiles`
  (both `TRUE`) to suppress either annotation.

* **`plot_comp()`'s aggregated composition pools counts across years** rather
  than averaging proportions, so a year with 20 otoliths no longer carries the
  same weight as one with 2000, and gains a 95% interval and the input /
  effective sample sizes (`add_agg_ci`, `add_agg_n`, both `TRUE`). The interval
  is the range holding 95% of the data the fitted model predicts, not a
  confidence interval on the mean, so observations outside it indicate misfit.
  Its variance is the exact sum of the per-year variances under each fleet's own
  likelihood -- the same variances the Pearson residuals are divided by -- with
  a normal approximation for the interval itself, which is poor below about 10
  expected counts. It also treats the fitted proportions as known, so it is
  slightly narrower than one carrying the estimation uncertainty in `p_hat`.
  Rows past `endyr` are excluded, matching the C++ likelihood's own gate, so a
  projection row no longer inflates the input sample size.

  The annotation names **two** effective sample sizes rather than one, because
  they answer different questions and the gap between them is the reweighting
  decision. `ESS (likelihood)` is what the model assumed, recovered from the
  same per-year variances the band is built on -- `sd^2 = p(1-p)/N_eff` holds
  for every family, so it needs no family branching and covers a
  Dirichlet-multinomial's overdispersion as readily as a multinomial weight.
  `ESS (McAllister-Ianelli)` is the tuning target this fit's own residuals
  imply, the harmonic mean across years `fit_mod()` already computes, which is
  the unbiased scale to average a ratio estimator on (averaging `Neff` directly
  runs about 50% high). Labelling either one "ESS" alone invited it to be read
  as the other. The McAllister-Ianelli line is drawn on multinomial fleets only:
  `reweight_comps()` names and skips a Dirichlet-multinomial fleet, because it
  estimates its own weight inside the likelihood, so printing an external tuning
  target beside one would invite the adjustment the package refuses to make.

  The effective sample size belongs to the observation, so it is summed over
  `(fleet, species, sex, year)` -- `data_check()`'s own uniqueness key -- not
  over years. A fleet may record a female-only and a male-only row in the same
  year, which are two observations; pooling them would report half the effective
  sample size beside an `ISS` that correctly counted both, reading as a fleet
  downweighted by half. Joint-sex rows (`Sex = 3`) remain one observation,
  because one multinomial spans both sexes.

  `afscOSA`'s constructions were deliberately **not** ported here. Its band is a
  binomial at the pooled sample size, and its aggregate effective sample size a
  single ratio taken from the pooled composition; measured against the installed
  package on data simulated with no misspecification, that ESS returns a median
  of 1.16 x ISS with a range of 0.34-9.16, because pooling first discards the
  between-year replication.

## Output format

* **`residuals(type = "pearson")` gains an `Sd` column** on composition sources:
  the standard deviation the fleet's own likelihood assumes for the observed
  proportion, i.e. the denominator the residual was divided by. `NA` on index
  and catch. The aggregated composition band reads it, so the band and the
  residuals cannot disagree about the assumed variance.

  On the copy carried as `attr(osa_residuals(), "pearson")`, which renames every
  column into that object's style, it is `assumed_sd` -- not `sd`, which on the
  OSA frame beside it is the conditional standard deviation of the OSA residual,
  a different quantity.

# Rceattle 5.28.1

## Documentation

* **Which selectivity forms can give the two sexes different levels.**
  `vignette("model-options-and-functionality")` described `Sel_norm_bin` and
  `Sel_norm_scope` as though they decided whether males and females could be
  selected at different levels, and its `"Max"` / `"AcrossSexes"` row promised
  that "one sex peaks at 1, the other keeps its relative level". Normalization
  only *preserves* a ratio the form already has -- `"AcrossSexes"` divides both
  sexes by one pooled scalar -- so on a logistic both sexes measure 1.00 and the
  promise is empty.

  The section now says which forms can carry a level difference and which
  cannot, with the ratios measured on `GOAatf`: `Logistic` 1.00,
  `DoubleLogistic` 1.97, `2DAR1` 4.53. Going non-parametric does not help
  either, which is the counter-intuitive part: those forms re-centre each sex to
  a mean of one every year, so the ratio is a by-product of shape rather than
  something the data inform. `inst/dev/two-sex-selectivity-example.R` fits the
  comparison and is now tracked.

## Bug fixes

* **A per-sex selectivity linkage no longer fixes the other sex too.** When a
  linkage supplies the base level -- a fixed intercept (`est_phase = 0`), or a
  slope-only formula -- `build_map()` masks the base parameter it replaces. The
  selectivity branch masked it across **every** sex of the fleet, so a linkage
  stratified on one sex (`by = ~ fleet + sex, sex = 2`) also fixed the other
  sex's inflection or slope at its starting value. That is exactly the sex a
  Stock Synthesis-style offset parameterization needs left estimated: the
  reference the offset is measured from. The `M` and growth branches beside it
  always honoured the row's own sex; selectivity was the outlier and now
  matches them.

  A linkage with no sex stratum -- `by = ~ fleet`, the default -- still masks
  both sexes, because the offset then applies to the whole fleet. Only a
  *selectivity* linkage naming a single sex changes, and nothing in the package
  or the reference fits does that: the golden models' only linkages are on
  composition, which this branch never touches.
  `test-linkage-selectivity-per-sex.R` covers both cases, which nothing did
  before.

# Rceattle 5.28.0

## Input format

* **`Sel_norm_bin` takes a word.** It says where selectivity is normalized to 1,
  and encoded three intents in one number with nothing in the value to say
  which: blank meant "do not normalize", any negative meant "by the maximum",
  and a positive number was an absolute age. It now also accepts `Max`, `Off`
  (or `None`), and `All` on a `LogisticPM` fleet, where the column is a penalty
  age-range rather than a normalization reference. Matching is case-insensitive,
  and `switch_check()` writes the word back, so a saved workbook says what it
  does. `Sel_norm_bin_upper` takes `Off` the same way.

  **Every existing value keeps its meaning.** Blank, `0`, `-1`, `-999` and a
  positive bin reach the model as exactly the code they always did; a test pins
  all of them.

* **A normalization bin the fleet is not selected over is now an error.**
  Previously `switch_check()` silently *overwrote* it: anything above `nages`
  was clamped to `nages`, which is a bin count used as a ceiling for an absolute
  age, so on a `minage = 3` stock ages 11 and 12 were valid and became 10. The
  clamp also ran before `data_check()`, so nothing downstream could see the
  original value. The bounds are now the fleet's own selected range -- below
  `Bin_first_selected` the curve is zeroed, so a reference taken there divides
  by nothing -- and a value below it is read as `Max`, which is what a negative
  has always meant. Note `Bin_first_selected` is a 1-based bin ordinal while
  `Sel_norm_bin` is an absolute age; the conversion is handled for you.

  On a stock recruiting at age 0, `Sel_norm_bin = 0` is the first age and is
  kept as a reference bin rather than read as a flag.

* **An unreadable value is refused rather than read as "do not normalize".**
  `rearrange_data()` is exported and does not always run behind
  `switch_check()`, so it checks the column itself; a typo such as `Maxx` used
  to resolve to blank and silently turn normalization off.

## Bug fixes

* **An `env_data` year outside the model years is dropped on the linkage path,
  rather than refusing the fit.** Rows are matched to model years by position
  (row `r` is model year `styr + r - 1`), so a row before `styr` shifts every
  later row and feeds the wrong covariate to consumption and, under
  multispecies, to predation mortality. `rearrange_data()` has dropped such rows
  from `env_index` since 5.25.0, so a fit without a linkage already tolerated
  them; a fit with one stopped in `.check_env_data_years()` with "env_data$Year
  must start at the model start year". The same workbook read two ways depending
  on whether it used a linkage.

  Dropping a row is a **warning**: covariate rows are assessment input, and 15
  workbooks in the sibling repositories carry pre-`styr` rows. A row whose `Year`
  is `NA` is kept for the year check to reject by name -- it is unlabelled, not
  out of range -- and a table lying entirely outside the model years is an error.

* **Rows after `projyr` are dropped too, which can move a linkage fit.** They
  shift nothing, since alignment runs from the front. But the fixed part of a
  linkage formula goes to `model.matrix()`, so `cut()`, `poly()` and `scale()`
  were computed partly over years the model never fits, and the basis they build
  depends on the rows supplied. Dropping them keeps the design matrix on the
  model years, and moves such a fit: measured at +434.9 on `~ cut(Year, 3)` and
  +823.8 on `~ scale(temp)` at identical parameters, with `~ (1 | Year)`
  unchanged because a random effect on an unused year integrates out.

  **No assessment is affected.** Of the 171 workbooks carrying `env_data` across
  `../Rceattle-models`, `../GOA-ATF-ESP`, `../GOA-multispecies-assessment` and
  `../GOA_circlulation_study`, 15 have rows before `styr` and **none** has a row
  after `projyr`. No bundled dataset has either. Refit only a model whose
  `env_data` runs past its `projyr` *and* carries a linkage whose fixed formula
  transforms a covariate.

# Rceattle 5.27.0

## New features

* **A standard error for the selectivity curve.**
  `fit_control(selectivity_se = TRUE)` makes [TMB::sdreport()] return
  `log_sel_at_age`, so a fitted selectivity can be drawn with a confidence
  interval instead of as a bare line. `exp(value +/- 1.96 * sd)` stays positive
  and is right-skewed, the same treatment `log_biomass`, `log_ssb` and `log_R`
  have had since 5.20.0. `log_sel_at_age_index` is reported alongside it and
  gives each row's fleet, sex, **absolute age and calendar year** -- `minage +
  index` and `styr + index`, not array positions.

  **On the log scale, not the logit.** A logit needs selectivity in (0, 1), and
  the non-parametric forms renormalize to a *mean* of one rather than a maximum
  of one: on `Atka2022`, 58% of `sel_at_age` exceeds 1 and the largest entry is
  3.06, so a logit is undefined across most of the array.

  Off by default. The delta method forms a Jacobian of every reported value
  against every parameter, so the cost is the product of the two: `Atka2022`
  adds 1,012 values against 584 parameters, and a 16-fleet three-species
  assessment around 31,000. The objective is untouched either way -- an
  `ADREPORT` is a report, not a likelihood term.

  **Only cells that are estimated and strictly positive are reported**, because
  one `log(0) = -Inf` on the tape turns *every* quantity in that `sdreport` into
  `NaN`, biomass and ssb included, not just selectivity. So the report skips
  `Selectivity = "Fixed"` fleets (read from `emp_sel_obs`, which leaves 0 wherever
  a year has no empirical row -- 228 cells on `BS2017SS`'s `EIT_Pollock`),
  length-based fleets (`sel_at_age` there is `sum(growth_matrix * sel_at_length)`,
  0 for an age overlapping no selected length bin), mirrored fleets and array
  padding (neither adds a curve), and every bin below `Bin_first_selected`, where
  `selectivity.hpp` zeroes the curve outright. The skipped set is fixed by the
  data, never by a parameter value. Nothing is lost: a `Fixed` curve and a mirror
  both have their lead's standard error or none by construction.

  A reported *value* can still underflow, which the skip rule cannot anticipate:
  a double-normal peaking above the oldest age with a narrow ascending limb puts
  `exp(-x^2/2)` below 1e-308, and 11 cells of `Atka2022`'s first year go to
  exactly 0 that way. The log is therefore floored at 1e-300, using a `max2()`
  that is exact above the floor -- on a fitted `Atka2022` all 1,012 reported
  values are bit-identical to `log(sel_at_age)` taken without it.

  A bin the normalization pins carries an SE of exactly 0, the honest answer for
  a value that is 1 by construction.

  **The error belongs to whichever `sdreport` the fit ends on**, which is not
  always the fit that estimated the curve. Under `estimateMode = "Estimate"` with
  an estimating HCR, `fit_mod()` re-optimizes the projection with every hindcast
  parameter mapped off, so the standard error of every selectivity comes back
  exactly 0. Use `estimateMode = "Hindcast"`, or
  `fit_control(projection_uncertainty = TRUE)`. `estimateMode = "Projection"`
  estimates no selectivity at all, and `projection_uncertainty` cannot change
  that, because `build_map(debug = TRUE)` has already mapped the hindcast off.
  `fit_mod()` warns on each of these before fitting, and warns again if
  `getsd = FALSE` leaves nothing to report.

* **`plot_selectivity(add_ci = TRUE)`** draws that interval under the curve it
  belongs to. It needs a fit run with `fit_control(selectivity_se = TRUE)` and
  declines with a warning naming the option otherwise, rather than drawing a bare
  line that would read as certainty -- and declines for the whole figure if any
  model in an overlay lacks the errors, since a band on some panels and not
  others says "this fleet is certain" rather than "this fit was not asked".

  A fleet mirroring another's `Selectivity_index` reports no rows of its own, so
  it borrows its lead's band, but only where the two curves agree: `data_check()`
  merely warns when a shared group differs in a shaping column, and a differing
  `Sel_norm_bin` alone rescales the curve. Every year drawn gets its own band, so
  pair `add_ci` with `minyr` / `maxyr` on a time-varying fleet.

  A fit reporting an error of 0 at *every* age is declined the same way, since
  the band would be a hairline reading as certainty. That is what an HCR
  projection returns; the warning names the `estimateMode` to refit with.

## Bug fixes

* **`build_map()` names the column and fleet when a switch reaches it as `NA`.**
  `Selectivity` or `Time_varying_sel` arriving unset on an active fleet failed
  at the first `==` against it, with "missing value where TRUE/FALSE needed" and
  no indication of which column or fleet. The message names `data_check()`, which
  is the validator -- `build_map()` already calls `switch_check()`, and that
  canonicalizes rather than validates. It keeps the validator's exemption: a
  fleet that is `Off` may leave both unset, which now also survives
  `data_check()`'s own selectivity-form checks and the random-effects sigma
  loop, neither of which was gated on `Fleet_type`. `Catchability`, where `NA` is
  legal on any fleet, is no longer compared as though it were a string.

* **A 2DAR1 or 3DAR1 fleet no longer warns about `Time_varying_sel = "Off"`.**
  The AR1 forms supply their own annual deviations, so any other
  `Time_varying_sel` mode is overridden and the builder says so -- but `"Off"`
  asks for no extra deviation layer, which is what an AR1 fleet already has, so
  it is not a conflict. Same class of guard as the `"Hake"` fix in 5.3.0.

* **`report_tables()` names its objectives for what they are.** `model` gains
  `marginal_nll` / `joint_nll` in place of `marginal_objective` /
  `joint_objective`, and `jitter` gains `best_nll` / `worst_nll` / `nll_range`
  in place of `best_objective` / `worst_objective` / `objective_range`. All six
  are negative log-likelihoods on the `jnll_comp` scale, so a **smaller** value
  is the better fit; `_objective` said nothing about which, and `_likelihood`
  would have argued for the wrong direction. `fit$opt$objective` is `nlminb()`'s
  own name and is unchanged.

* **`fit_mod()` no longer prints the "Model did not converge" banner.** It fired
  when an `sdreport` was requested and did not return, which is a Hessian that
  would not invert -- not an optimizer that failed to converge, and a fit can
  reach a good optimum with a non-positive-definite Hessian. The convergence
  battery reports the same condition as `sdreport_failed`, says what it means,
  and prints at `verbose = 0` as well. `fit$identified` is unchanged.

* **Two long-standing `R CMD check --as-cran` warnings are cleared.**
  `print.rceattle_report()`'s box-drawing glyphs are written as `\u2514\u2500`
  / `\u251c\u2500` escapes, so the R code is ASCII and the package is
  portable -- the printed tree is byte-identical. `tibble` is declared in
  `Suggests`, which `test-data-check-age-coverage.R` has used (guarded) for
  several versions. CI runs `error-on: "error"`, so neither had ever failed a
  build.

* **`quantity_dictionary()` described `sel_at_age` as "normalized to a maximum
  of one".** It is normalized per `Sel_norm_bin`, which for the non-parametric
  forms -- the default -- is a mean of one, so values above one are ordinary. It
  is also 0 below the fleet's first selected bin. The units column said "proportion" for
  the same reason and now says "unitless". `log_sel_at_age` and
  `log_sel_at_age_index` are documented alongside them.

# Rceattle 5.26.0

## New features

* **`parameter_index()` locates every estimated parameter in the model's own
  coordinates.** TMB names each element of `obj$par` after its parameter block,
  so a diagnostic could report that `sel_coff_dev` was non-identifiable but not
  which of its 392 elements. Returns one row per estimated parameter with
  `species`, `fleet`, `sex`, `age`, `bin`, `year` and `slot` columns plus a
  rendered `label`. Recovered through TMB's own `parList()`, so it tracks
  `build_map()` exactly: mapped-off parameters are absent, and fleets sharing a
  `Selectivity_index` or `Catchability_index` appear once, with `n_cells`
  counting the cells they drive. An axis constant across the model is left out
  of the label; the structured columns carry it either way.

  Ages are absolute, and an axis that omits the youngest ages says so:
  `init_dev` reports ages `minage + 1` upward, because `ceattle.cpp` reads
  `init_dev(sp, age - 1)` over `age = 1 .. nages-1` and age `minage` in the
  first year is recruitment. Reading that axis from `minage` would have named
  every initial-age deviation one age too young and never named the oldest.

  It describes `fit$obj`, which under `estimateMode = "Estimate"` with any HCR
  but `"NoFishing"` is the *projection* object -- `log_Ftarget` and `log_Flimit`
  alone, since `build_hcr_map()` maps every hindcast parameter off. Pass a fit
  run with `estimateMode = "Hindcast"` to browse the hindcast parameters.

## Bug fixes

* **The convergence battery names the parameters it flags.** `estimability`,
  `parameters_on_bounds` and `hessian_conditioning` reported the parameter block
  and left the reader to find the element. They now summarise by coordinate,
  ordinal axes collapsed to ranges:

  ```
  [FAIL] estimability  49 non-identifiable fixed parameter(s)
    sel_coff_dev  GOA_pollock_fishery, bins 1-8, 2013-2018  (41)
  ```

  `parameters_on_bounds` gains a `where` column.

  The three do not all count positions in the same vector: `estimability` and
  `parameters_on_bounds` are captured at the hindcast optimization, while
  `hessian_conditioning` reads the `sdreport`, which under an estimating HCR
  belongs to the projection fit. `fit_mod()` now stores the hindcast index
  alongside the hindcast parameter vector, and each check verifies that the index
  it uses names that vector before labelling anything -- a coordinate naming the
  wrong parameter is worse than no coordinate, so a mismatch prints block names
  alone.

* **`hessian_conditioning` is read on the correlation matrix.** The covariance
  condition number is not scale-invariant (`log_F` near -2, `sel_inf` near 10),
  so rescaling a parameter moved it without changing the model. The message now
  leads with the square root, a ratio of standard errors -- and, because the
  matrix is scaled by those errors, one measured against each parameter's own
  rather than in the parameters' own units.

  **`data$condition_number` therefore means something different from 5.25.1.**
  Standardising lowers it by 0.4-1.1 orders on the models measured, so at
  unchanged `WARN` 1e6 / `FAIL` 1e10 the check fires less readily; a badly
  scaled but unconfounded fit is no longer flagged. `data$se_ratio` and
  `data$covariance_condition_number` are new, the latter carrying the old value.
  Code comparing condition numbers across this boundary needs the new name.

* **Selectivity slot labels follow the fleet's `Selectivity`.** Slot 2 of
  `sel_inf` is a descending inflection only for the double-logistic family: it
  is `logit(right_floor)` under `DoubleNormal` and the free age-1
  log-selectivity under `LogisticPM`.

* **`sel_curve_pen` slots are AR1 correlations.** The block is estimated only
  under 2DAR1 / 3DAR1, where its slots hold logit-scale correlations across
  bins, years and cohorts -- not the `fleet_control` shape and curvature penalty
  weights, which are never estimated.

* **`rec_pars` slot 2 is the SRR alpha**, not steepness; steepness is derived
  from alpha and `SPR0`. `growth_log_sd` is the standard deviation of
  length-at-age at the minimum and maximum age (`L1` / `Linf`), and
  `log_growth_pars` now names its four slots rather than numbering them.

* **Three `parameter_dictionary()` dimensions were wrong.** `log_F` declared
  `[n_fsh, nyrs]` against an `[n_flt, nyrs_hind]` array; `log_pop_scalar`
  `[nspp, nyrs]` against an age axis; `M1_dev_log_sd` `[nspp, nsex, 2]` against
  `[nspp, nsex]` (also corrected in `ceattle.cpp`). Documentation only, no fit
  changes. `test-schema-parameter-index.R` now checks every block's declared
  rank *and* extent against the built array.

* **`fit_mod()` no longer prints the `check_estimability` table at
  `verbose = 1`.** It is one row per parameter, each named after its block, and
  the convergence record now carries the same verdict by coordinate. Still at
  `verbose > 1`; `fit$identified` is unchanged.

# Rceattle 5.25.1

## Bug fixes

* **The discontinuous-likelihood warning compares two marginal objectives.**
  Convergence warning (8) tested `nlminb`'s objective against
  `obj$report()$jnll`. The first is the Laplace **marginal** negative
  log-likelihood whenever the model carries a random effect; the second is the
  **joint** negative log-likelihood at the random-effect mode. They differ by
  the Laplace correction, which is far larger than the default `rel_tol` of 1 --
  150.96 on a GOA pollock fit carrying two catchability linkages -- so the
  warning fired on every such model regardless of whether its likelihood was
  smooth. It is now tested against a fresh `obj$fn()` evaluation, which agrees
  to machine precision on all four configurations checked (an `ar1`/`rw` q
  linkage, no random effect at all, `random_rec = TRUE`, and `2DAR1`
  selectivity with `random_sel = FALSE`).

  The comparison is made beside each optimization rather than at the end.
  `fit_mod()` rebuilds the TMB object for the projection, and again under
  `projection_uncertainty = TRUE`, without re-optimizing, and `build_hcr_map()`
  maps every hindcast parameter off, so those rebuilds disagree with `opt` on
  which parameters are random. Comparing across a rebuild is the same units
  mismatch: on a `BS2017SS` fit with one `~ (1 | Year)` catchability linkage,
  `opt$objective` of 10207.488519 against 10228.940923 from the
  all-parameters-on rebuild, a gap of 21.45 with nothing wrong with the fit.
  The same mismatch reached `HCR = "ConstantF"` under `estimateMode = 0`,
  which rebuilds for the projection and never re-optimizes.

  The `random_rec == FALSE` guard is removed with it: it exempted recruitment
  random effects from a check that could not read them, and the comparison now
  holds for any random effect, so `random_rec = TRUE` models get the check
  rather than skipping it. No fit changes -- the objective, the joint negative
  log-likelihood and the parameter vector are identical either side of this.

# Rceattle 5.25.0

## Breaking changes

No exported argument is removed or renamed. Two configurations that previously
ran are now refused, and **neither produced a number anyone should use**. A
workbook with the same column under two spellings is also refused where
it previously read, silently keeping one of them -- one file in the sibling
repositories does (`GOA cod/Data/2024_cod_pcod_test.xlsx`, whose
`Bin_first_selected` and `Age_first_selected` disagree on every row):

* **`random_sel = TRUE` with non-parametric selectivity**, in two cases, each
  for its own reason. Fit these with `random_sel = FALSE`, which is the
  penalized formulation the AMAK model intends and is unaffected.

  Under `Time_varying_sel = "RandomWalk"` the walk is scored on the per-year
  renormalized selectivity, which never sees the level of each year's
  coefficients. Those directions are an **improper** random effect, and the
  deviation sd collapses: 2.7e-8 on `Atka2022`, a time-invariant selectivity
  reported as a time-varying one, reached through 17 `NA/NaN function
  evaluation` warnings.

  Under `"IID"` the density itself is proper, but the AMAK shape penalty beside
  it is one-sided (`max(d, 0)^2`), whose second derivative is a step. The
  Laplace correction is a log-determinant of that second derivative, so the
  marginal objective is only piecewise smooth and the optimizer halts at a kink:
  `Atka2022` stops with a maximum gradient of **6.8** and reports a deviation sd
  **27%** away from the value the same model reaches with `Sel_curve_pen1 = 0`
  (maximum gradient 4e-4). Setting that penalty to 0 integrates the deviates
  against the curvature penalty alone and is allowed.

* **`Ceq` outside 1-4.** The template switches on it with no other branch, so an
  unrecognized value left the temperature function at zero, and with it
  consumption and -- under multispecies -- predation mortality.

## Bug fixes

* **Non-parametric selectivity no longer returns `NaN` when a coefficient grows
  large.** The forms renormalize their log-scale coefficients against
  `log(mean(exp(coefficients)))` once per year, and the realized curve was
  stored on the natural scale and re-logged by every penalty. Both halves
  overflowed: the sum reached infinity once a coefficient passed about 709, and
  the round trip underflowed to `log(0) = -Inf` from about -745, either way
  giving `Inf - Inf = NaN` with no gradient to recover from. A new
  `log_mean_exp()` computes the mean with a max shift -- algebraically identical,
  and the shift cancels out of the derivative too -- and the centred curve is now
  held on the log scale rather than recovered from `exp()`. Measured on
  `BS2017SS`, the objective and gradient are finite to a coefficient of 1e8 in
  both the level and the spread direction, where they previously failed at 710
  and 730. The four reference models reproduce their optima to 1e-11 relative.

* **The non-parametric shape penalty cannot read past the last bin.** Its
  pairs are `(bin, bin + 1)`, so the left bin cannot be the last one, but
  `data_check()` bounds `Sel_pen_last_bin` by the number of bins. Set to the
  final bin, the penalty read one bin beyond the curve -- an entry that exists
  because the selectivity arrays are sized by the widest species, but that is
  never written -- and charged a shape penalty against a phantom
  fully-selected bin. The model clamps the range.

* **A model with no environmental covariate runs.** `Ceq` 1 (Stewart et al.
  1983), 2 (Kitchell et al. 1977) and 3 (Thornton and Lessem 1979) each read
  `env_index`, which `ceattle.cpp` evaluates for every species on every fit, but
  `data_check()` only guarded `Ceq > 1` and described 1 as
  temperature-independent. With no covariate column that read went outside the
  matrix -- observed both segfaulting and returning plausible numbers on
  different runs, the model being built without bounds checking. Consumption now
  falls back to the constant `fT = 1` of `Ceq = 4` where there is no column to
  read. `data_check()` errors on it under multispecies, where consumption drives
  predation mortality and so the fit, and warns under single-species, where the
  predation section does not run and `fT` reaches only the reported `ration` and
  `consumption_at_age`. The missing-value defaults become `Ceq = 4` (the only
  equation that reads no index) and `Cindex = 1` (it is 1-based; the documented
  default of 0 was outside the matrix).

* **An `env_data` year outside `styr:projyr` no longer shifts the series.**
  `env_index` is read by row position -- row `r` is model year `styr + r - 1` --
  so an observation from before `styr` displaced every later row, feeding an
  earlier year's covariate to consumption and, under multispecies, to predation
  mortality. **This moves real assessments.** No bundled dataset is affected, but
  15 workbooks in the sibling assessment repositories have pre-`styr`
  environmental rows and need refitting: the GOA CEATTLE 2018-2023
  configurations (13-14 years), AI cod and a single-species Pacific hake input
  (12 years each), and three GOA circulation-study assessments (2-9 years). The
  Pacific hake operating model that `04-mse.R` runs is **not** among them --
  `env_data` starts at its `styr` of 1980 -- so the MSE and predation reference
  objectives are unchanged. A repeated `Year`, which sits inside the range and so
  survived the filter, is now refused by `data_check()` and de-duplicated before
  the model sees it.

* **An integer `Time_varying_sel` on a non-parametric fleet is a mode again.**
  A pre-4.4 workbook stored the two AMAK penalty weights in the time-varying
  columns, and the upgrade that moves them across was triggered by the *value*
  in `Time_varying_sel` -- which cannot tell a shape weight of 4 from the
  `RandomWalk` code 4. `Time_varying_sel = 4` on a modern workbook therefore
  turned the walk off and overwrote the penalty weights (on `GOAcod`, 20 and
  12.5 became 4 and the deviation sd), while the string `"RandomWalk"` worked.
  The upgrade is now triggered by an empty `Sel_curve_pen1`, which is what
  actually distinguishes the old format.

* **`data_check()` refuses a non-parametric fleet whose `Bin_first_selected` is
  above its `N_sel_bins`**, which left it with no estimated coefficients at all.

* **`convergence_diagnostics()` reports a variance parameter estimated to
  zero.** A deviation standard deviation at the floor means the process it
  governs is not varying: a model configured as time-varying has fitted
  something time-invariant. `Atka2022` under `random_sel = TRUE` with a
  non-parametric random walk reaches `sel_dev_sd = 2.7e-08`; the battery flagged
  that fit through `max_gradient`, which says the optimizer did not converge
  rather than naming what went wrong. The new `variance_collapse` check names
  the parameter and the element within it, and fires whether or not the gradient
  is clean -- a collapse at a converged optimum is the case with nothing else to
  catch it. It reads the recruitment, selectivity, catchability, M1 and linkage
  random-effect standard deviations, which are log-scale and O(0.1)-O(1) so one
  threshold of 1e-3 means the same thing for each; observation sds,
  `growth_log_sd` and the catchability prior sd are excluded. `WARN` rather than
  `FAIL`, because a variance at the boundary is a well-posed optimum and a
  verdict about the model rather than the optimizer.

* **Two composition-accumulation columns were never being read.**
  `Accumatation_age_lower` and `Accumatation_age_upper` are the pre-4.4 names of
  `Comp_accum_young` and `Comp_accum_old`, and were missing from the schema's
  alias list -- so a workbook with them got the default of no accumulation
  and the setting was ignored. **147 of the 183 fleet_control sheets in the
  sibling assessment repositories have them, and 79 set a young-tail value
  that folds bins**; the GOA pollock configuration is the AFSC one, folding ages 1-2 into
  age 3 for the Shelikof survey and age 1 into age 2 for the fishery. On the
  2023 GOA pollock workbook, applying it moves the composition negative
  log-likelihood from 16333.05 to 16039.59 at identical parameters. Both are
  1-based ordinals on the fleet's composition dimension, which equals the age
  only when `minage` is 1 -- it is 1 on every affected workbook. No bundled
  dataset has either column, and none of the 147 also has the
  canonical spelling, so nothing had to be reconciled.

  `Est_weights_mcallister` is deliberately **not** treated the same way. Despite
  looking like a third case, it is the pre-4.12 name of an OUTPUT --
  [fit_mod()] computes `Comp_weights_mcallister` from the fitted object and
  already mirrors the old name for back-compatibility -- so a copy of it in a
  workbook is a saved result, not a setting, and reading it as an input would
  imply an effect on the fit that it does not have.

* **A deprecated column name is no longer silently discarded when its canonical
  name is also present.** `read_data()` canonicalises `fleet_control`, so a
  script that then assigns a deprecated name creates the old column beside the
  new one -- and the upgrade deleted it and kept the canonical value, behind a
  routine deprecation message. The two must now hold the same setting, with an
  integer switch code and its string counting as the same; anything else is an
  error naming both values. Nothing is merged and no blank is filled in from the
  other column: `NA` is a real setting in several of these (`Sel_norm_bin` and
  `Sel_cap_bin` mean "do not normalize" and "no cap", `Proj_F_proportion` means
  "no F apportioned"), so filling could not express clearing a value and would
  move a number without saying so.

  **Check your own scripts before upgrading**, because this has been hiding real
  settings and the error will now find them. One script across the sibling
  assessment repositories stops on it, at exactly such a setting: the Jack
  mackerel 2024 bridging script assigns `Sel_sd_prior = 0.35` on four fisheries
  after `read_data()` has canonicalised the workbook's own value of 40 to
  `Time_varying_sel_sd`, so the model was fitting effectively free deviations
  where the script asked for penalized ones. The fix is to assign the canonical
  name instead of the deprecated one; the error names both.

## New features

* **`plot_catchability()`** draws fitted survey catchability by year, faceted by
  fleet, for every fleet with `index_data` -- a fishery with a CPUE series as
  much as a survey. It reads `quantities$index_q`, the realized `q` the model
  scaled the survey by, so one figure covers an estimated or fixed mean, annual
  deviations under `Time_varying_q`, an environmental regression, a `q` linkage
  such as `rw(1 | Year)`, and the closed-form `"Analytical"` value. Hindcast
  years only: the model does not project catchability, so there is no
  `incl_proj`. A time-invariant `q` draws a flat line, which is the honest
  picture of a constant `q` rather than an empty panel.

* **Non-parametric selectivity accepts `Time_varying_sel = "IID"`.** The model
  scores `dnorm(sel_coff_dev, 0, Time_varying_sel_sd)` on each estimated
  coefficient, so the annual curves are independent deviations about an
  estimated base curve rather than a walk. `NonParametricPM` still takes
  `"RandomWalk"` only -- its deviates *are* the walk increments, so an
  independent-deviate reading of them would not describe the curve the model
  builds, and its error message says so and points at `NonParametric`.

  This is what bundled **`Atka2022`** is configured as, so that dataset -- and
  the ADMB bridge it comes from -- could not be fitted as shipped. It now fits.

# Rceattle 5.24.2

## Bug fixes

* **`standard_output()` spells the length-bin column `length_bins`, and emits
  `beg_mid`.** The schema `stockplotr` and `asar` read is the 34 columns of
  `stockplotr::example_data`; `convert_output()` spells it `len_bins` internally
  and renames on its last step, so 5.24.0 emitted the internal name. Those
  packages index by column name, so the wrong spelling is an error in the
  consumer rather than an empty column. `test-report-tables.R` now checks the
  emitted names against that package rather than against the constant itself.

* **`report_tables()` reports convergence from the hindcast fit.** `max_gradient`
  and `pdHess` came from `$opt` and `$sdrep`, which the harvest-control-rule
  projection re-optimizes and overwrites, with a fallback to a bare `obj$gr()`
  evaluated at whatever parameter vector was left in place. Both now read the
  `.conv_hindcast` snapshot [fit_mod()] takes before the projection -- the same
  source `convergence_diagnostics()` uses -- so a maximum gradient quoted in an
  executive summary describes the assessment.

* **The reference-point gating reads a harvest control rule held as an integer.**
  `build_hcr_map()` matches `HCR` by name, so a `data_list` holding the integer
  code matched no rule and `report_tables()` reported `NA` for an F40% that was
  estimated, with a `basis` saying it was not. `HCR` is normalized before the
  lookup.

# Rceattle 5.24.1

## Bug fixes

* **Spawning biomass per recruit is no longer wrong for a two-sex species.** SPR
  is spawning output per **total** recruit, so the female fraction belongs in it
  for every species -- but it must enter once. A one-sex species' schedule is
  sex-combined, so the age-varying `sex_ratio` applies at every age, which is
  what `mature_females` folds in (section 5.4). A two-sex species' sex-0
  schedule is already female, so only the recruitment split `sex_ratio(sp, 0)`
  applies -- the same split section 6.6 uses to divide recruitment between the
  sexes. Section 6.2 applied the age-varying ratio to both, which re-applied a
  split already in the schedule and returned `NaN` wherever that table stopped
  short of the species' own `nages`.

  Found on the 2026 GOA three-species assessment, where arrowtooth flounder
  (`nsex = 2`, `nages = 21`) returned `NA` for every per-recruit quantity and so
  could not be given a Tier 3 SPR reference point. **One-sex species are
  unchanged in every configuration**: `mature_females` already applies the
  age-varying ratio for them, so the new formula reduces to the old one.

  **Which two-sex fits move.** The change rescales all four SPR schedules
  wherever `sex_ratio` varies with age, and it reaches the objective two ways:

  1. The reference-point penalty (section 13) scores `SPRlimit / SPR0` and
     `SPRtarget / SPR0`, so it is live only under `HCR = "ConstantF"` or an
     SPR-based rule (`HCR > 3`).
  2. The stock-recruit derivation (section 6.3) is live under **any** harvest
     control rule: with `srr_fun = "BevertonHolt"` or `"Ricker"`, `SPR0` sets
     `steepness` and `R0` and `SPRFinit` sets `R_init`, which scales the initial
     age structure (section 6.5). `srr_fun = "mean"` ignores `SPR0`, which is why
     the four golden references reproduce their pinned objectives.

  Route 1 also needs the rescale to be **non-uniform**. It scores a ratio of two
  SPR schedules, so where `sex_ratio` is the same at every age with
  `maturity > 0` the old and new formulas give the same ratio and `Ftarget` /
  `Flimit` do not move -- only the reported per-recruit quantities do. Route 1
  moves a fit only where `sex_ratio` varies across the *mature* ages, or where
  the schedule stops short of `nages` and the old formula returned `NaN`.

  That is the split across the live assessments, all two-sex arrowtooth under
  `HCR = "NPFMC"` with `srr_fun` at its `"mean"` default:

  - The **2026 GOA three-species assessment** had `sex_ratio` filled only to age
    10 of 21, so every per-recruit quantity was `NaN`. It gains reference points
    it could not previously form. **Refit.**
  - The **GOA arrowtooth ESP** (`../GOA-ATF-ESP`) has the schedule fully
    populated at 0.5 for age 1 and 0.353955 for ages 2-21, and `maturity` of 0
    at age 1 -- a uniform 0.708. Its `Ftarget` and `Flimit` are unchanged; the
    `SPR0`, `SPRlimit`, `SPRtarget` and `SPRFinit` it reports were 0.708 of the
    correct value and now rise by 1.412. **No refit; recheck any quoted SPR.**
  - The **Pacific hake operating model** is two-sex with a flat 0.5, so the MSE
    and predation reference is unchanged in both respects.

  Bundled `GOAatf` meets the condition -- `sex_ratio` is 0.5 at age 1 and 0.354
  over the mature ages -- so the old formula rescaled every SPR by a clean factor
  of 0.708 with nothing in the output to give it away. Built with
  `srr_fun = "BevertonHolt"` and no harvest control rule, the old schedules give
  steepness 0.4655 against 0.5516 and `R_init` 0.9505 against 1.0623 at the
  starting values: -15.6% and -10.5% on a fit route 1 never touches.
  `test-dynamics-spr-two-sex.R` checks the result against the per-total-recruit
  rule.

  `ssb`, `biomass`, `R`, `F_spp`, `M_at_age`, `N_at_age`, `SB0` and the objective
  are unchanged for a one-sex species, and for a two-sex species on neither
  route.

* **`data_check()` validates `nsex`.** The model has exactly two sex schedules --
  index 0 (females, or the single combined sex) and index 1 (males) -- and reads
  `nsex` straight into loop bounds and array dimensions, so a value outside
  `{1, 2}` was a silent out-of-range read, and an `NA` aborted `data_check()`
  with R's "missing value where TRUE/FALSE needed", naming no table. `nsex` is
  now checked for length and value, and the `M1_base` / `weight` / `ration_data`
  sex-consistency checks look it up by the table's own `Species` code rather
  than by row position, so an unusable or short `nsex` is reported alongside
  every other error instead of thrown or silently recycled.

* **`data_check()` no longer demands a full `sex_ratio` schedule for a two-sex
  species.** Every remaining use of `sex_ratio` on a two-sex species is
  `sex_ratio(sp, 0)` — the age-1 recruitment split — so requiring values at every
  age rejected a workbook the model never reads past the first column. A one-sex
  species still needs the whole schedule, and age 1 is still required for both.

# Rceattle 5.24.0

## Reading a fitted model

* **`quantity_dictionary()` says what every reported quantity is.** `fit$quantities`
  holds 99 derived quantities under the model's own abbreviated names, and nothing
  said what they meant, what units they were in, or which had a standard error.
  The new table gives each one a plain-language meaning, its units, its dimensions,
  whether `sdreport()` provides an SE for it, and what the same quantity is called in
  the NOAA standardized assessment output:

  ```r
  quantity_dictionary("ssb_depletion")
  quantity_dictionary(process = "reference_points")
  quantity_dictionary()[quantity_dictionary()$se, ]   # everything with an SE
  ```

  It is checked against the template and against a real fit by
  `test-schema-quantity-dictionary.R`, so a quantity added, renamed or dropped in
  `ceattle.cpp` fails a test rather than leaving the table quietly wrong.

  Two things it records that were easy to get wrong. Units are **mt** and
  **thousands of fish** throughout (numbers-at-age are thousands, weight-at-age is
  kg, so their product is mt) — several C++ comments still said `kg`. And every
  per-recruit reference point (`SPR0`, `SPRlimit`, `SPRtarget`, `SPRFinit`,
  `NbyageSPR`) is computed only under `msmMode = 0` and is **exactly zero on a
  multispecies fit**, which reads as an estimate of zero unless you know.

* **`report_tables()` collects what an assessment reports into one set of tables.**
  Previously the likelihood decomposition, the time series with uncertainty, the
  reference points, the fits, Mohn's rho, the jitter and the OSA diagnostics each
  came out of a different call in a different shape, so comparing two models meant
  assembling a dozen ad-hoc extractions by hand. Sections follow the AFSC Alaska
  Groundfish Stock Assessment Guidelines for what a SAFE chapter reports, and every
  table has a `model` column so several fits stack for comparison:

  ```r
  tabs <- report_tables(list(base = fit0, alt = fit1),
                        retro = list(retro0, retro1), osa = osa0)
  tabs$reference_points
  ```

  It **never refits**. A retrospective or a jitter is tens to hundreds of
  optimizations, so they are passed in already computed; a section whose object is
  `NULL` is simply absent. The standard harvest scenarios of guideline section
  4.11.3 are *not* produced — they need a standard projection module Rceattle does
  not have. Projected biomass under the model's own harvest control rule is in
  `timeseries` with `era = "fore"`.

* **`report_tables()` gains a `parameters` section, and reports both
  objectives.** Asked for directly: where are sigma_R, the estimates of M, and
  the likelihood values? `parameters` joins `coef()` and `vcov()` to
  [parameter_dictionary()], so every estimated parameter gets its
  natural-scale name, its process and its standard error. Estimates are on the
  parameter's own scale, so a `log_` name needs `exp()`; sigma_R appears only
  when `random_rec = TRUE`, and a **fixed** M is not there at all because it was
  never estimated — read that off `M_at_age`.

  `model` now reports `marginal_objective` and `joint_objective` rather than one
  ambiguous `objective`, plus `n_random`. The marginal is what the optimizer
  minimized and what `AIC` uses; the joint is what `jnll_comp` decomposes, so it
  is what the `likelihood` table sums to. They are equal without random effects
  and differ by the Laplace correction with them — on `GOAatf` with random
  recruitment, 744.33 against 816.99. Reporting only one made the two tables
  look inconsistent.

* **`standard_output()` emits the NOAA standardized assessment format.** Relabels
  `report_tables()` output into the schema that the `stockplotr` and `asar` packages
  consume, so Rceattle results can be plotted and written into a report by the same
  tooling used for SS3, BAM, WHAM and FIMS. Names are translated through the
  dictionary's `standard_label`, so `ssb` becomes `spawning_biomass` and `F_spp`
  becomes `fishing_mortality`; a quantity the standard has no word for keeps its own
  name rather than being dropped.

  That standard describes **one stock and has no species dimension**, so a
  multispecies fit errors unless `species` selects one, rather than returning a frame
  in which two stocks' biomass share a year. The likelihood decomposition is the
  exception: `jnll_comp` counts fleets on rows 1-8 and species on rows 9-20, so
  fleet rows go in the standard's `fleet` column whatever the species
  rather than filtered away with it — on a 16-fleet assessment they are 31 of
  the 38 rows.

  Verified against the 2026 GOA three-species assessment (3 species, 16 fleets,
  one two-sex stock, a live `sdreport`, and its real retrospective, jitter and
  OSA objects), not only against the package fixtures. `quantity_dictionary()`
  covers exactly the 99 quantities both its single- and multi-species fits
  report.

* **A reference point a fit never estimated is `NA` with a stated `basis`, not
  the number the array happens to hold.** CEATTLE leaves a value behind in three
  cases where there is no reference point, and each one reads as an estimate:

  - `Ftarget` / `Flimit` are estimated only under a harvest control rule that
    defines them, and are switched off for a species with no projected fishery.
    Unestimated they sit at `exp(0) = 1`, so the GOA assessment reported a
    target F of **1.0/yr**. The gating is taken from `build_hcr_map()` rather
    than by reading the HCR switch a second time; the fit's own `map` cannot be
    used, because `build_map()` sets both to `NA` in the hindcast map whatever
    the HCR.
  - Under `msmMode > 0`, `SB0` / `B0` are overwritten by the `MSSB0` / `MSB0`
    inputs, which stand at a 999 mt placeholder until `fit_mod()` derives them.
    `B_target = Ptarget * SB0` was therefore reported as **399.6 mt** against a
    true scale of 1e5-1e6 mt. `MSSB0_derived` is the flag that distinguishes a
    placeholder from a genuinely derived value.
  - The per-recruit quantities are computed only under `msmMode = 0`.

  The depletions are deliberately **not** blanked alongside `SB0`: under a
  no-fishing rule in multispecies mode the model divides by biomass in the last
  projection year, the equilibrated unfished reference, so that series is
  meaningful.

* **A diagnostics list is matched to models by name.** `retro`, `jitter` and
  `osa` pair with `object` by name whatever the order; an unnamed list is paired
  positionally and says so; a name that is not a model name is an error. That
  last one catches passing a single model's `osa_residuals()` result stored as a
  list of parts (`index` / `catch` / `comp`), which would otherwise be
  attributed one part per model.

## Corrections

* **Five stale or swapped comments in `ceattle.cpp`.** No behaviour change; they were
  found while documenting the quantities, and each one described the code as doing
  the opposite of what it does. `Flimit` / `Ftarget` and `Flimit_at_age` /
  `Ftarget_at_age` had their descriptions the wrong way round (`Flimit` is the FOFL
  proxy, `Ftarget` the maximum FABC proxy); `NByageF` said `F = Flimit` where the
  recursion uses `Ftarget_at_age`; and `index_hat` and `catch_obs` were annotated
  `(kg)` where they are mt or thousands of fish, per the fleet's
  `Observation_units`.

# Rceattle 5.23.0

## Reading a fitted model

* **`parameter_dictionary()` is now exported.** CEATTLE's parameter vector uses
  transformed, abbreviated names (`log_M1`, `R_log_sd`, `index_log_q`), and the
  table explaining them existed only as an internal object with `@noRd`, so it
  generated no help page and nobody outside the package could find it. It now
  has a user-facing accessor that returns the internal name, the quantity on its
  natural scale, the process it belongs to, a one-sentence meaning and the
  dimensions, filterable by name or by process:

  ```r
  parameter_dictionary("R_log_sd")            # -> sigma_R
  parameter_dictionary(process = "selectivity")
  ```

  An unknown `process` is an error; an unknown parameter name only warns, so
  passing every name from a fit still works when that fit carries parameters
  from a branch the dictionary predates. `test-schema-parameter-dictionary.R`
  asserts that every parameter `build_params()` creates is documented, so the
  table cannot silently fall behind the model.

* **`AIC()` replaces `fit$opt$AIC` in the introduction vignette.** That slot is
  not always there. `.fit_tmb()` normalizes what `TMBhelper::fit_tmb()` returns,
  and when TMBhelper hands back a nested list the normalizer keeps the inner
  `opt` — where `AIC` never lived — so the slot silently disappears. Without
  TMBhelper the plain `nlminb` fallback has no `AIC` either. `AIC(fit)` works in
  every case, because `logLik.Rceattle()` builds it from `opt$objective` and
  `length(opt$par)`, both of which survive all three paths. The vignette does not
  execute by default and `test-vignette-api.R` checks call signatures rather than
  return shapes, so the stale line was invisible to CI.

* **The introduction vignette now covers the fitted object itself** — the S3
  methods an `Rceattle` fit responds to (`summary()`, `plot()`, `coef()`,
  `vcov()`, `logLik()`, `residuals()`, `simulate()`, `profile()`,
  `as.data.frame()`), the classes the diagnostics return, and where the numbers
  live in `$quantities` versus `$estimated_params`.

# Rceattle 5.22.0

## Diagnostics

* **`retrospective()` warns when a peel is dropped, and the object remembers how
  many were asked for.** Mohn's rho is averaged over the peels that survive, so a
  drop changes a number that gets reported — but it was announced by a
  `message()`, which `suppressMessages()` hides, and nothing on the returned
  object showed it afterwards. It now warns, carries `$peels_requested`, and
  `print()` reads "3 of 5 peel(s)" with a `NOTE` and says what rho was averaged
  over. `jitter()` and `self_test()` still message.

  Worth knowing: `getsd = TRUE` also turns on the non-positive-definite Hessian
  check in `.refit_converged()`, so it can drop peels that `getsd = FALSE`
  keeps, and rho can differ between the two runs of the same model.

* **A retrospective peel now carries uncertainty on its hindcast.** Each peel is
  fitted twice — a peeled hindcast, then a forecast refit that estimates only the
  peeled years' F — and the second holds the hindcast fixed, so it reported a
  standard error of **zero** for every hindcast year. Peels plotted with no
  confidence band while the unpeeled model had one.

  Under `getsd = TRUE` the peel is now reported from the same parameters with the
  hindcast free. Nothing is re-estimated: on BS2017SS under an NPFMC HCR, Mohn's
  rho, every peel's SSB and every peel's objective are unchanged. SSB standard
  errors in a common year now fall as data are retained — 248,693 with two years
  peeled, 220,476 with one, 179,416 with none — instead of reading zero. Costs
  one extra model build per peel.

* **`fit_control(projection_uncertainty)` survives a refit.** It is carried on
  `data_list` and read back by `.refit_like()`, so a peel, jitter or MSE refit
  inherits how the model was fitted rather than taking `fit_control()`'s `FALSE`
  default — the failure the bias-adjustment flags already guard against. Under an
  HCR it decides whether a fit reports 6 free parameters and no hindcast
  uncertainty, or 316 and real errors. Forwarded only when `getsd = TRUE`.

* **`profile(joint =)` moves a whole set of cells on one grid.** `slots` cross
  through `expand.grid()`, so moving a ten-age M schedule together over 13
  values was a ten-dimensional factorial — `13^10` fits, not 13. That made an
  empirical age-based M, the usual case when M is supplied as a vector rather
  than estimated, impossible to profile the way it is normally reported.

  `joint = "multiply"` applies one grid to every named cell as a multiplier on
  its current natural-scale value, preserving the **shape** of the schedule
  (ratios between ages) and moving it up or down as a whole; `1` is the fitted
  model. `joint = "add"` adds the grid value instead, preserving the
  **increments** between cells; `0` is the fitted model. `joint = "value"` sets
  every cell to the same value. `"none"` is the default and unchanged.

  Both act on the natural scale, so they need `transform = "log"` or
  `"identity"` — applying a multiplier to a stored log is a different model,
  not a different parameterisation, and an arbitrary `transform` cannot be
  inverted. On a log-scale parameter a grid value that would drive a cell to
  zero or below is an error naming the cell and the value, not an `NaN` fit.

* **`profile(param = "q")` profiles a fleet's base catchability.**
  `index_log_q` is one entry per fleet, named by `Fleet_name`, so `slots` takes
  a **fleet name** as well as an index.

  Fleets sharing a `Catchability_index` share one q parameter — the map copies
  the lead fleet's slice across the group — so fixing one member's cell alone
  would leave the rest estimating a common q, which is not that fleet's q
  profiled. Naming any member now profiles the whole group: the others are
  added to `slots` with a message, and `joint` becomes `"value"` so the group
  takes one q rather than a factorial over its members.

* `profile()` returns `$joint` alongside `$alias`, so `print()` and the figure
  axis can say the grid holds a multiplier or an offset rather than a parameter
  value. A grid running 0.7 to 1.3 under an axis reading `M1[1, 1, 1]` would be
  read as an M of 0.7. A catchability axis names the fleet, not its number.

* **`plot_profile(relative = "scaled")` puts every component on its own scale.**
  On a real profile one component routinely dwarfs the rest — on a Bering Sea
  pollock M profile the bottom-trawl composition moves 60.5 objective units
  while the bottom-trawl index moves 0.697, a factor of 87 — so the large one
  sets the y axis and every other curve is drawn as a flat line along the
  bottom. Where those components prefer the parameter cannot be read at all,
  which is the one thing the figure exists to show.

  `"scaled"` re-zeroes each series at its own minimum and then divides by its
  own change over the grid, so every curve runs 0 to 1 and the minima can be
  compared directly.

  **It discards magnitude**: a component moving 0.02 and one moving 40 draw
  identically, so the figure stops saying whether a component matters and says
  only where it points. Raise `minfraction` alongside it. `minfraction` and the
  legend ordering are computed on the **raw** change over the grid, before any
  re-zeroing or rescaling — computed afterwards every span would be 1, the
  filter would drop nothing, and scaling would have quietly removed the guard
  at exactly the moment it is needed.

  `relative` keeps `"own"` (the default, unchanged), `"minimum"` and `"none"`.

## Bug fixes

* **`run_mse()` no longer expands a time-invariant `weight` or `ration_data`
  series into the projection.** A series supplied at `Year 0` is time-invariant
  — `rearrange_data()` fills every hindcast year from the single row — so it
  needs no projection rows. Expanding it replaced one legal row with a series
  naming the projection years and no year before them, and `data_check()` reads
  any index carrying more than one year as time-varying and requires it to span
  `styr..endyr`. The first assessment advances `endyr`, the check runs, and the
  operating-model refit dies with

  ```
  Weight data for index = 4 & sex = 1 does not span all hindcast years
  ```

  recorded as a bare `"OM"` failure with the message discarded, so every
  simulation returned `use_sim = FALSE` and `mse_summary()` had nothing to
  summarise. Hit on a Pacific hake MSE whose workbook supplies weight indices
  4–6 at `Year 0`; a workbook that dates every weight row never reached it. Both
  expansions now skip a group whose rows are all `Year 0` and are shared by the
  operating and estimation models through one helper.

  Two further defects in the same expansion, found reviewing the fix. It read
  the group's **last table row** rather than its last row at or before `endyr`,
  so a series running past its own terminal year projected a value the hindcast
  never fitted — not the terminal hindcast year `?run_mse` documents holding.
  And it added a row for every projection year without checking whether one was
  already there: `rearrange_data()` assigns row by row into the weight array, so
  the appended duplicate sorted last and its value silently replaced the
  observation. A workbook with catch and weight through 2023 against `endyr`
  2019 hit both, losing four years of weight-at-age. `clean_data()` already
  skips an already-present projection year for catch, and the `NByageFixed`
  expansion beside this one does the same.

* **`run_mse()` refuses an operating and estimation model whose terminal or
  projection years disagree.** The assessment loop takes the row positions of
  the catch rows to fill from the estimation model's table and indexes the
  operating model's `max_catch_hat` with them, so the two tables have to
  describe the same rows in the same order. Mismatched years read the
  exploitable-biomass limit off the wrong years, or off `NA`, and the run
  completed either way — catch advice that is wrong without looking wrong. The
  requirement was always there; now it is checked.

* **`run_mse()` no longer skips every assessment whose year already carries
  catch.** The fill that hands the control rule's advice to the operating model
  selects projection rows on `is.na(Catch)`. A projection year arriving with a
  number in it was passed over: no catch was written, the interval's total came
  out `0`, and the operating model took the no-catch path — rebuilt at
  `estimateMode = 3` with its map dropped, rather than advanced. The next
  section then overwrote that year with the rebuilt model's realized catch, so
  the recorded value did not survive either; the assessment simply did not
  happen.

  This is what a workbook looks like when `endyr` has fallen behind the catch
  series. On a Pacific hake MSE workbook carrying catch through 2023 with
  `endyr` still 2019, four of six assessments (2020–2023) ran that way, because
  `clean_data()` creates an `NA` projection row only for a year with no row
  already — leaving 2024 and 2025 as the only fillable years. Nothing in the run
  reported it.

  Catch past either model's terminal year is now blanked to `NA` at setup, with
  a warning naming the model, the years, and the terminal year it read against
  — the three things needed to tell a stale `endyr` from data that belongs in
  the projection. Nothing the model uses is
  lost: the likelihood scores only `Year <= endyr`, and no control rule reads
  observed catch (HCR 1, constant catch, sums `catch_hat`). The blanking runs
  ahead of the `Proj_F_proportion` check, which compares recorded catch against
  projected exploitable biomass and was tripped by the same rows into forcing a
  rebuild the model did not need.

## Documentation

* The diet age-coverage warning is about half its former length. It led with
  two clauses saying the same thing and closed with three sentences of
  qualification, which buried the list of species and ages — the only part that
  differs between one fit and the next. The gap lines keep their consequence in
  a parenthesis (`never eaten`, `exerts no predation`) rather than a trailing
  clause. Same rule, same species reported.

# Rceattle 5.21.0

## Diagnostics

* **`profile_components()` and `plot_profile()`: read a profile for conflict,
  not just for precision.** `profile()` returned the total negative
  log-likelihood, which answers how well the data determine a parameter but
  not whether the data sources agree about it. A smooth, well-behaved total
  can be two surveys pulling in opposite directions, and nothing in the total
  says so.

  `profile_components()` pulls `quantities$jnll_comp` out of every grid fit and
  returns it as one long data frame — one row per grid point per component,
  each cell labelled with the fleet or species it belongs to. That labelling is
  the work: the columns of `jnll_comp` count fleets on the data, selectivity
  and catchability rows and species on the priors, penalties and predation
  rows, so a cell read off the wrong axis would attribute a survey's likelihood
  to a species.

  `plot_profile()` draws it, following `r4ss::SSplotProfile()`: change in
  negative log-likelihood on the y axis, every series re-zeroed at its own
  minimum with a point marking it, the total overlaid in black, and components
  moving less than `minfraction` of the total dropped. The spread of the points
  along the x axis is the disagreement. `plot()` on a profile is the same
  thing. `relative = "minimum"` re-zeroes at the fitted optimum instead, and
  passing a list of profiles facets them side by side.

  Under `random_rec = TRUE` the total is the Laplace-approximated marginal
  likelihood while the components are the inner joint negative
  log-likelihood; they will not sum, and `profile_components()` says so rather
  than reconciling them.

* **The QAR1 process error is reported as a deviate, not as a prior.** The AR1
  density on `index_q_dev` (Rogers et al. 2024) was scored into the
  "Catchability prior" row of `jnll_comp`, where it read as a prior on log q.
  It now scores into "Catchability deviates", and accumulates rather than
  assigns, so it cannot erase anything already in the cell. **No fit moves**:
  the objective is `jnll_comp.sum()`, so relocating a term between rows leaves
  it bit-identical, and the branch is in any case unreachable —
  `data_check()` refuses `Catchability = 6` and the live QAR1 form is a q
  linkage, `ar1(1 | Year)` with `observe`, which scores under "Linkage random
  effects". Fixed rather than left in place because a component profile reads
  these row labels to say which term is in conflict.

  `tests/testthat/test-schema-jnll-rows.R` is new and pins the three-way
  registry the `jnll_comp` rows depend on — the `JnllRow` enum, the row labels
  in `R/6-rename_output.R`, and the column axis each row is indexed by. Nothing
  linked them at build time before, so a row added to one and not the others
  was silent.

* **`profile()` reports the parameter in the units you profiled it in.**
  `profile(fit, param = "M1", ...)` resolves `param` to `"log_M1"`, so
  `print()` reported a minimum of, say, 0.28 next to the name `log_M1` when
  0.28 was an M and not a log M. The alias is now kept in the returned object
  as `$alias`, and both `print()` and the figure's axis use it. `$param` is
  unchanged — it is still the internal slot the fit mapped off.

# Rceattle 5.20.0

## Diagnostics

* **`self_test()` and `profile()` join the display contract.** 5.16.0 gave
  `convergence_diagnostics()`, `osa_diagnostics()`, `retrospective()` and
  `jitter()` one way of reporting — an overall status in the same four words,
  a one-line verdict, then a compact table. These two were left out, so
  `vignette("model-diagnostics")` claimed a contract the family did not all
  keep.

  `self_test()` now reports the fraction of simulations that reached an
  optimum, which the returned list could not say on its own: non-converged
  runs are dropped before it is returned, so the number of fits returned is
  not the number attempted. The number attempted is carried in
  `attr(, "nsim")`. Beneath it, a tally of each returned fit's own
  `$convergence$status` — a separate question, since a replicate can reach a
  zero gradient and still carry a `NOTE` — and a line naming what generated
  the replicates, because with processes held fixed the spread across them
  carries observation error only and is a lower bound on estimation
  uncertainty.

  `profile()` now reports whether the grid **brackets** the minimum. A profile
  whose lowest point is its first or last grid value has run out of grid rather
  than found the optimum, and it plots as a perfectly ordinary curve, so
  nothing in the numbers says so. For a one-dimensional profile it also gives
  the values within 1.92 of the minimum — the 95% profile-likelihood cutoff,
  \eqn{\chi^2_1(0.95)/2} — marked open on either side the grid does not close.
  A cross-profile gets the edge test but no interval.

  **Return values are unchanged.** Both objects are the list they always were:
  `sims[["Sim_1"]]`, `length(sims)`, `prof$grid$slot_1` and
  `prof$nll - min(prof$nll)` all index as before, `c(sims, list(fit))` and
  `sims[i]` return a plain list of fits, and `plot_biomass()` and
  `compare_sim()` are unaffected. Only `class()` gains an entry and `print()`
  gets an opinion.

## Other changes

* `R/9-retro_and_jitter.R` was split into `9-diagnostic_helpers.R`
  (`rceattle-refit-args`, `.parallel_lapply()`), `9-retro_and_jitter.R`,
  `9-self_test.R` and `9-profile.R`. No code moved between functions and no
  behaviour changed.

# Rceattle 5.19.1

## Bug fixes

* **A `maturity` or `sex_ratio` held as a matrix no longer aborts
  `data_check()`.** The per-species age-coverage check added in 5.19.0 read the
  `Species` column with `$`, which on a matrix is an error rather than `NULL`:
  the check died with `$ operator is invalid for atomic vectors` and took every
  error accumulated before it with it. `read_data()` returns both tables as
  data frames, but a hand-built `data_list` may hold either as a matrix, and
  every other check tolerated one. Read by column name now.

* **A `Species` column that is not a species number is reported.** It was
  compared with `as.integer()`, so a column of species *names* gave `NA`, `NA
  != i` gave `NA`, and `which()` dropped it — leaving the row order unchecked,
  which is the one thing that block exists to check. A factor column is now read
  by its label rather than its level code: the two agree only while the numbers
  run `1..nspp`, so a table carrying species 1, 3, 4 read as level codes 1, 2, 3
  and passed.

* **Every `nages` `NA` reports rather than aborting.** 5.19.0 added `na.rm` to
  the age-coverage comparisons, which covers one `NA` among several; with
  *every* value `NA` the maximum is `-Inf`, and two checks downstream still
  consumed it — the CAAL column list as `1:-Inf`, which errors with "result
  would be too long a vector", and `NByageFixed` as a demand for `-Inf`
  columns. Both are gated on a finite maximum now, so `nages` is reported by
  the check that owns it.

## Documentation

* `vignette("hcrs-and-mses")` attributed the per-fleet `Index_distribution`
  index draw to 5.13.0; it landed in **5.9.0**. It also had 5.13.0 breaking
  reproducibility for every seeded `run_mse()`, where it breaks it for two
  kinds of run — one using a scalar `cap`, and one whose control rule reads
  reference points from an estimated stock-recruit curve. Both the vignette and
  the 5.13.0 `NEWS.md` entry now say which.

* The 5.18.0 entry gives the runtime the common-random-numbers fix cost:
  refitting the operating model over its whole projection was worth **9% on a
  Bering Sea multispecies MSE projecting to 2040**, growing with the length of
  the projection. `inst/dev/TODO-mse-horizon.md` carries the design for getting
  it back.

* The version-numbering note at the top of this file had a sentence with words
  missing.

# Rceattle 5.19.0

## Bug fixes

* **`data_check()` checks `maturity` and `sex_ratio` age coverage per species,
  not only against the widest one.** The existing checks compare the table's
  column count with `max(nages)`, so a table wide enough for the longest-lived
  species passes even when a shorter-column species is left short of its own
  bins — and the `[0, 1]` range check then reads over the `NA`s with `na.rm`.

  Both tables are summed over age. `mature_females` is `maturity`, times
  `sex_ratio` where a species is modelled one-sex (`ceattle.cpp` 5.4), and feeds
  hindcast `ssb`; `spawning_biomass_per_recruit()` (`src/TMB/spr.hpp`) sums the
  same schedule for SPR. So a `maturity` gap makes SSB `NA` for any species, and
  a `sex_ratio` gap does the same for a **one-sex** species.

  On a **two-sex** species `sex_ratio` reaches only SPR — which is how this one
  hid: the model fits cleanly, and the failure waits until a harvest control
  rule asks for reference points, when `nlminb` reports `NA/NaN gradient
  evaluation`, naming neither the table nor the species.

  Found on the GOA multispecies workbook, where arrowtooth `sex_ratio` was
  filled to Age10 against 21 age bins. That same data now reports:

  ```
  sex_ratio is missing values for species 2, ages 11-21 of 21. Spawning
  biomass and spawner-per-recruit both sum over every age, so a gap leaves
  SSB or SPR0 NA.
  ```

  Ages past a species' own `nages` are padding and are not checked, so a ragged
  multispecies workbook still passes — asserted against the bundled datasets in
  `test-data-check-age-coverage.R`.

  Rows are read by **position**, since `rearrange_data()` drops the `Species`
  column and hands the model a matrix whose row *i* is species *i*. A
  `Species` column that disagrees with the row order, and a table with fewer
  rows than species, are both reported.

* **An `NA` in `nages` no longer aborts `data_check()`.** Six age-coverage
  comparisons (`weight`, `ration_data`, `age_error`, `maturity`/`sex_ratio`
  widths, `NByageFixed`, CAAL columns) evaluated `any(... < nages)` or
  `max(nages)` without `na.rm`, so a single `NA` made the `if()` condition `NA`
  and raised "missing value where TRUE/FALSE needed" — discarding every error
  accumulated up to that point, including `nages`' own "not included for all
  species".


# Rceattle 5.18.1

## Documentation

* **The workbook column descriptions are reference entries, not essays.** The
  `meta_data` sheet's `Description` column is what an assessment author reads
  while filling in a workbook, and several cells had grown past 1,500 characters
  — `Catchability` and `Catchability_index` the longest. Each now gives what the
  column means, its allowed values and its default, then points at
  `vignette("model-options-and-functionality")`, which already carried the
  sharing rules and the lead-fleet resolution in full. No allowed value, default
  or switch code was dropped, and `inst/extdata/meta_data_names.xlsx` is
  regenerated from the schema. The same trim is applied to the matching
  `?BS2017SS` field dictionary in `R/data.R`.

* **"The template" is called "the model" in user-facing text.** It is TMB build
  jargon for `ceattle.cpp`, and it collided with the *Excel* template
  `write_template()` writes — two different objects under one word, in help
  pages and error messages. `data_check()` and `sim_mod()`'s messages, the
  vignettes and the roxygen now say "the model"; `write_template()`'s own
  documentation still says "template", meaning the workbook.

* Code comments state current behaviour rather than what an earlier release did,
  per the repository's comment standard. No behaviour, return value or number
  changes anywhere in this release.

* The `NEWS.md` entries for the unreleased 5.9.0–5.18.0 line are condensed:
  every breaking change, deprecation, moved number and reproducibility boundary
  is kept, while derivations and design rationale defer to the vignettes that
  carry them.

* The `SB0`/`SBF`-under-predation limitation in `ceattle.cpp` is marked
  `KNOWN LIMITATION` and names `DynamicSB0` / `DynamicSBF` as the consistent
  alternative, rather than carrying a review-process `TODO` tag.

# Rceattle 5.18.0

## MSE

* **Two `run_mse()` schedules run on the same `seed` are now on common random
  numbers**, so a paired comparison of two assessment schedules measures the
  schedules. Each assessment's observation draws are seeded on its own **year**,
  from a table built once per replicate, and the operating model is refit over
  its whole projection rather than only as far as the next assessment. On
  BS2017SS, contamination in a year whose advice is identical by construction
  falls from **2.1% to 0.003%**; the remainder is optimizer noise.

  **Common random numbers are still incomplete.** One `sim_mod()` call draws
  every year in an assessment interval under that assessment's seed, so a year
  sitting inside a longer interval is drawn under a different seed than the same
  year in a schedule that assessed it directly. `?run_mse`,
  `vignette("hcrs-and-mses")` and `inst/dev/TODO-mse-horizon.md` say so.

  **Every seeded `run_mse(simulate_data = TRUE)` result changes.** Runs stay
  reproducible from a given `seed` within a version, and `simulate_data = FALSE`
  is unchanged to the bit.

  **It costs runtime.** Refitting the operating model over its whole projection
  rather than only as far as the next assessment was worth **9% on a Bering Sea
  multispecies MSE projecting to 2040**, growing with the length of the
  projection; single-species models saw little change. `projyr` sizes the AD
  tape, so the saving was real -- it just made the draws depend on when the next
  assessment fell, which is the one thing a comparison of two schedules has to
  hold fixed. `inst/dev/TODO-mse-horizon.md` carries the design for recovering
  it without that dependence.

## Performance

* **`run_mse()` no longer refits the operating model on a shortened horizon.**
  That saving was 9% on a Bering Sea multispecies MSE projecting to 2040; it
  goes because the shortened horizon was what made the draw count depend on the
  schedule. `inst/dev/TODO-mse-horizon.md` records how to recover it.
  `tools/verify/verify-mse-om-horizon.R` is replaced by the `gap_crn_clean`
  check in `tools/verify/verify-mse-schedule.R`.

# Rceattle 5.17.0

## MSE

* **`run_mse(assessment_period =)` takes an explicit vector of assessment
  years.** A single number is the period it always was; a vector is the exact
  years an assessment is completed, which is how a *skipped* assessment is
  expressed — a period cannot represent one gap in a regular cycle.

  ```r
  biennial <- seq(om$data_list$endyr + 2, om$data_list$projyr, by = 2)
  run_mse(om, em, assessment_period = setdiff(biennial, 2031))
  ```

  Every year must fall after the operating model's terminal year and within the
  projection horizon; a year outside it is an error rather than being dropped.
  A schedule must name two or more years, since one year cannot be told from a
  period. A scalar must be a whole year of at least 1 that leaves room for one
  assessment inside the projection. Nothing that produced a usable number is
  refused.

* **A schedule that stops short of the projection horizon now warns.** Catch is
  filled only up to the last assessment, so trailing years keep the `NA` their
  projection rows were created with, while `mse_summary()` summarises over the
  whole projection — understating `Average Catch`, `Catch IAV` and `P(Closed)`
  with nothing in the table to say so. A biennial cycle over an odd number of
  projection years lands here every time. The warning names the years and both
  fixes: extend the schedule, or set `projyr` to the last assessment year.

* **`run_mse(catch_mult =)` takes a `data.frame` of `Year`, `Species` and
  `mult`.** A single number or a vector of length `nspp` is unchanged and
  applies in every projection year; a `data.frame` applies only to the year and
  species pairs it lists, and every omitted pair is multiplied by 1. This
  expresses a buffer held only while advice is stale, which the vector form
  cannot.

  ```r
  buffer <- expand.grid(Year = 2032:2033, Species = seq_len(om$data_list$nspp))
  buffer$mult <- 0.90
  run_mse(om, em, assessment_period = setdiff(biennial, 2031),
          catch_mult = buffer)
  ```

  Mind which years are stale: the assessment in year `Y` sets catch for `Y+1`
  onward, so a missed 2031 assessment leaves **2032 and 2033** on 2029's advice.
  `Species` is the species number. An uncovered year, an out-of-range species, a
  non-finite or negative multiplier, and two rows for the same year and species
  are all errors — each would otherwise read as *no reduction* in a run that
  finishes and looks ordinary. The year bound is the last assessment year, not
  `projyr`. A vector `catch_mult` is now checked for finiteness and sign.

  Note this multiplies **catch**, not the ABC the control rule produces. Where
  realized catch sits well below ABC — GOA arrowtooth flounder, for one —
  reducing ABC changes removals only to the extent the fishery attains it.
  Scale by recent attainment, `1 - (1 - mult) * attainment`, or report the
  result as an upper bound. `?run_mse` and `vignette("hcrs-and-mses")` carry the
  caveat.

* No stored MSE result changes. `assessment_period = 1` and a scalar
  `catch_mult` take the same paths they did, and
  `tools/verify/verify-mse-schedule.R` checks that a schedule given as years
  reproduces the equivalent period exactly.

* **Documented: two schedules on the same `seed` are not on common random
  numbers.** Measured at **2.1%** on a year whose advice was identical by
  construction. Pre-existing, and **fixed in 5.18.0**.

## Documentation

* `vignette("hcrs-and-mses")` gains a worked missed-assessment scenario, and its
  reproducibility section now names **5.13.0** alongside 5.9.0. Both moved the
  random-number stream, so a seeded `run_mse()` reproduces across neither.

# Rceattle 5.16.0

## Breaking changes

* **`summary()` on a fit returns a summary object, not the fit.** It was
  `print(object)`. It now prints the spec tree, a parameter table with standard
  errors, and the likelihood components and their total, returning a
  `"summary.Rceattle"` with `$coefficients`, `$jnll_comp`, `$objective` and
  `$convergence`.

  `x <- summary(fit)` followed by `x$quantities` no longer works; use `fit`
  directly. Swept across `../Rceattle-models` and `../GOA-ATF-ESP`: every call
  is a bare `summary(mod)` at top level, so nothing downstream is affected.

* **`Time_varying_sel = "AR1"` and `Time_varying_q = "AR1"` (value `2`) are
  removed.** `data_check()` errors on them, naming the fleet and the
  replacement.

  Neither was ever an AR1: the model scored value `2` with the same independent
  normal penalty as value `1`, and neither deviation block has a correlation
  parameter. **Nothing that produced a usable number is refused** — every fit
  under value `2` was an IID fit, and `"IID"` reproduces it exactly. Of 173
  `fleet_control` sheets across the ecosystem, **zero** set `"AR1"` on either
  switch.

  `Time_varying_q` is exempt where it is not a mode: under
  `Catchability = "Environmental"` the column holds an `env_data` column index,
  so a literal `2` there is a series number, not a switch value.

  The replacement is the linkage grammar, which fits the stationary AR1 these
  names promised — `sigma` is the marginal sd, reducing to IID at `rho = 0`:

  ```r
  build_catchability(linkages = list(
    q = linkage_spec(~ ar1(1 | Year), by = ~ fleet)))

  build_selectivity(linkages = list(
    inf_asc = linkage_spec(~ ar1(1 | Year), by = ~ fleet)))
  ```

  It is also more expressive: per selectivity *parameter* rather than per fleet,
  with the deviation sd estimated or fixed, priors, bounds and a phase.

## Bug fixes

* **The unfished-reference placeholder is recognised by a flag, not by its
  value.** `mse_summary()` compared the reported `SB0` to the 999 mt
  placeholder — a float equality test that would also null a legitimately
  derived 999 mt, and that was blind to a workbook supplying its own `MSSB0`.
  `clean_data()` now seeds a per-species `MSSB0_derived`, which `fit_mod()` sets
  where it projects under no fishing. A fit saved before the flag existed falls
  back to the value test.

* **`print()` on a fit closes its tree and lines up its values.** `Run time` is
  `signif()` to match the objective rather than seven digits of wall clock, and
  an `estimateMode` that runs no projection reports
  `HCR : (not applicable -- hindcast only)` rather than `NoFishing`, which read
  as a management choice nobody had made.

* **`Time_varying_sel_sd` is validated, as `Time_varying_q_sd` already was.**
  `build_params()` takes its log, so a blank or non-positive value gave a
  non-finite objective from inside `MakeADFun`, naming neither the fleet nor the
  column.

  Required only where the model reads the sd, which is a property of the
  (`Selectivity`, `Time_varying_sel`) pair: the logistic family under `IID` /
  `RandomWalk` / `RandomWalkAscending`, `Hake` under `IID`, and
  `NonParametric`/`NonParametricPM` under `RandomWalk`. `LogisticPM`, `Block`,
  `Fixed`, `2DAR1` and `3DAR1` are excluded.

  Of 173 `fleet_control` sheets across the ecosystem, **one is refused** —
  `Pacific hake/Data/Old files/hake_from_ss.xlsx`, whose `Hake_fishery` pairs
  `Selectivity = "Hake"` and `Time_varying_sel = "IID"` with a sd of 0. That is
  `dnorm(dev, 0, 0)` and cannot produce a finite objective. Nothing that fitted
  is refused.

* **The model spec tree shows a parameter carrying more than one linkage.** The
  grammar lets one parameter hold a list of specs — a shared prior plus a
  fleet-specific random walk, which is how GOA pollock 2025 is configured — and
  `print()` rendered those rows as `NULL`. Each spec now gets its own line,
  tagged `[i/n]` when stacked.

## New features

* **`plot(fit, what = "ssb_depletion")`.** The S3 dispatcher offered only
  `"depletion"`, which is *biomass* depletion. In a package whose vocabulary is
  the Amendment-56 proxies, spawning-biomass depletion against B40% is the
  management quantity, and the generic could not draw it.

* **The diagnostics share one display contract.** `osa_diagnostics()`,
  `retrospective()`, `jitter()` and `mse_summary()` now open with the header
  `convergence_diagnostics()` uses — the object, a status, and a one-line
  verdict — before the detail.

  **No return value changed.** `retro$mohns`, `jit$nll`, `diag$sdnr`,
  `summ$species` and every column index as before; only `class()` gains an entry
  and `print()` gets an opinion.

  - `osa_diagnostics()` leads with how many sources failed and the overall SDNR
    against its null interval, then a severity-tagged line per source, worst
    first.
  - `retrospective()` judges the terminal peel against `print(retro, band = )`,
    default `+/- 0.2`. Forecast-skill peels are reported but not judged: a rho
    over a forecast horizon is not the quantity that rule was calibrated on.
  - `jitter()` gains `njitter`, so the fraction of starts reaching the best
    optimum has a knowable denominator — non-converged starts are dropped, so
    the count of returned fits is not the count attempted.
  - `mse_summary()` reports its four blocks and their dimensions, and says that
    `as.data.frame()` on the whole object will not work.

## Internals

* **The golden regression now runs in CI, and the multi-OS `NOT_CRAN` guard is
  deterministic.** `test-golden-regression.R` carries both `skip_on_cran()` and
  `skip_on_covr()`, so it fired in neither existing job and the committed
  backstop for golden equivalence was never automated. A new `deep-checks`
  workflow runs it nightly and asserts the file produced assertions, since a
  skipped file otherwise reports as a pass.

  The `NOT_CRAN=false` override moved to step-level `env:` on
  `check-r-package`, where GitHub's precedence makes it unshadowable. Through
  `$GITHUB_ENV` it did not hold reliably, and when it slipped the Windows worker
  died with `0xC0000005`.

* **The model can be built with array bounds checking.**
  `RCEATTLE_SAFEBOUNDS=true` builds with TMB's `safebounds`, turning an
  out-of-range access from a silent write into adjacent memory into an R error
  naming the object. `tools/verify/verify-safebounds.R` drives it, with a
  positive control asserting `-DTMB_SAFEBOUNDS` reached the compile line. The
  default is unchanged and no fit gets slower.

## Documentation

* **The ten vignette chunks the weekly job could not see now run.**
  `RCEATTLE_EVAL_VIGNETTES` flips each vignette's chunk *default*, but a chunk
  carrying its own `eval = FALSE` overrides it — so ten chunks across four
  vignettes stayed unexecuted in the one job built to execute them, including
  the OSA and process-residual workflows and three `run_mse()` variants. Exactly
  one chunk keeps an explicit `eval = FALSE`, the install block, and says so.

  Cost on an M-series Mac: `hcrs-and-mses` 2.7 -> 19.2 minutes,
  `introduction` 2.7 -> 3.8, `model-diagnostics` 1.9 -> 2.3. The weekly job's
  timeout went from 120 to 180 minutes, since a timeout loses every vignette's
  result rather than one.

* **`?mse_summary` says when `om_terminal_depletion` is `NA`** — any
  multispecies run without a harvest control rule, which derives no unfished
  reference. Points at `om_terminal_depletion_dynamic`, which is unaffected.

* **`vignette("model-parameterizations")` gives the fishery index predictor.**
  It still showed the survey snapshot for CPUE, which 5.9.0 replaced for a
  fishery with the year-average `N̄ = N (1 − e^{−Z}) / Z`. The predictor is now
  split by `Fleet_type`, with the Baranov derivation and the SS3
  correspondence. The survey half was wrong independently: it wrote
  `e^{-Month·Z}` where the model computes `exp(-(Month/12)·Z)`, an exponent out
  by a factor of twelve.

* Same vignette: the `Selectivity` list gained `DoubleNormal` (8),
  `NonParametricPM` (9) and `LogisticPM` (11) and the note that there is no 10;
  `N_sel_bins` is described in bins rather than ages; and the four
  "`Time_varying_sel = 2`, `sigma` is estimated" equation blocks now name what
  actually does that, `fit_mod(random_sel = TRUE)`.

* **Three switch-table rows in `vignette("model-options-and-functionality")` no
  longer break their table.** A `|` inside a table cell is a column separator
  even inside a code span, so `linkage_spec(~ ar1(1 | Year))` rendered a
  spurious fourth column. The replacement code moved out of the cells into
  fenced blocks below each table. Nothing catches this today: `R CMD check` does
  not execute the vignettes, `test-vignette-api.R` parses only the R chunks, and
  pandoc renders a ragged row without complaint.

* **`vignette("introduction")` plots SSB depletion where it says it does.** The
  "SSB depletion" example called `plot_depletion()` — biomass depletion, the
  same call as the line above it — so the two examples drew identical figures
  under different labels. It calls `plot_depletionSSB()` now.

* The `Time_varying_sel` soft-deprecation message now names `ar1(1 | Year)`
  alongside `rw(1 | Year)`, as the `Time_varying_q` message already did.

# Rceattle 5.15.0

## Behavior changes

* **`DynamicSB0` and `DynamicSBF` decay to spawning on their own mortality.**
  Both are carried forward on predation mortality solved at their own fishing
  rate (`M_at_age_dB0`, `M_at_age_dBF`) but then applied `M_at_age` — the
  *fished* mortality — over the fraction of the year to spawning. They now use
  the same M they were propagated on.

  Only a model with a non-zero `spawn_month` is affected, and only under
  predation: with `spawn_month = 0` the term is `exp(0)`, and with `msmMode = 0`
  the two mortalities are identical. No bundled dataset carries both, so
  `/golden-check` is bit-identical.

  The equilibrium pair, `NByage0`/`SB0` and `NByageF`/`SBF`, is **not** fixed
  here. `SB0` escapes it under `msmMode > 0`, being overwritten with the
  `MSSB0` input, but `SBF` does not: it applies the realized fished-predator
  `M2` to an F-target *equilibrium* numbers-at-age array. HCR 5 (NPFMC Tier 3)
  reads `SBF` and `SB0` together, so under predation the two legs of that rule
  sit on different mortality bases. `msmMode = 0` is unaffected — there is no
  `M2`. `DynamicHCR = TRUE` uses `DynamicSB0` / `DynamicSBF`, which are
  internally consistent. Marked `KNOWN LIMITATION` in `ceattle.cpp`.

* **A multispecies model can report depletion against `DynamicSB0`.** With
  `HCR = 0` the projection runs unfished, so SSB in the last projection year is
  multispecies SB0 once `projyr` is far enough out to equilibrate, and that is
  the reference `ssb_depletion` uses. That arm ran unconditionally, so it also
  overwrote the `DynamicHCR = TRUE` result, leaving no way to get a dynamic-B0
  depletion out of a multispecies model. It is now skipped under `DynamicHCR`.
  Runs with `DynamicHCR = FALSE` are unchanged.

## Bug fixes

* **A multispecies fit now carries its unfished SSB (`MSSB0`) on the returned
  object.** Under `msmMode > 0` the model discards its own equilibrium SB0
  and reads the `MSSB0` `DATA_VECTOR` instead, so `ssb_depletion` is
  `ssb / MSSB0`. No workbook can supply it — `clean_data()` seeds it at a 999 mt
  placeholder and `fit_mod()` derives the real value by projecting under no
  fishing. That derived value was written only into the reorganized copy the
  projection refits from, so the **returned** `data_list` kept the 999, and
  everything that refits from a fitted object — `.refit_like()`, `remove_F()`,
  and every `run_mse()` projection — re-entered the model with the
  placeholder.

  **Any multispecies refit carrying a harvest control rule moves.** `SB0` is not
  only the depletion denominator: it is also the HCR threshold and the
  `posfun` floor on `ssb_depletion - Plimit`, so a refit was comparing spawning
  biomass against 999 mt. A single `fit_mod()` call was already correct — the
  fit itself used the derived value — so `/golden-check` is bit-identical and
  no hindcast changes.

* **`mse_summary()` no longer reports SSB/999 as a terminal depletion.** An
  unfished reference is only derived for a run carrying a harvest control rule,
  so under `HCR = "NoFishing"` `MSSB0` is still the placeholder. The
  multispecies arm divided by it anyway: on the Pacific hake three-species
  model `OM: Terminal SSB Depletion` read **2.68e3** while
  `OM: Terminal SSB Depletion (Dynamic)`, which uses the model's own
  `DynamicSB0`, read a sane 0.96. It is `NA` now, which is what an undefined
  reference point should report. The dynamic column is unchanged.

  That arm also indexed `SB0[sp]` on an `nspp × nyrs` matrix, reading year 1
  rather than the terminal year. Masked today because the multispecies
  overwrite makes every year identical; corrected so it stays right if that
  changes.

# Rceattle 5.14.1

## Documentation

* **Two misspellings corrected where they described a switch value.**
  `vignette("model-parameterizations")` wrote "specificed" for the
  `Time_varying_sel = 3` block form in three places, and both that vignette and
  `build_map_selectivity()` wrote "selecitivty" for the
  `Selectivity = "NonParametric"` form of Ianelli et al. (2018). Text only; no
  switch, default or fit changes.

# Rceattle 5.14.0

## Behavior changes

* **`random_q = TRUE` now estimates the catchability deviation standard
  deviation.** `random_sel` has always freed `sel_dev_log_sd` alongside the
  selectivity deviates it integrates out; `random_q` integrated the catchability
  deviates out but left `index_q_dev_log_sd` mapped, so the deviations were
  smoothed at whatever `Time_varying_q_sd` the workbook happened to carry. The
  two are symmetric now: integrating the deviates out is what turns their sd into
  a marginal variance rather than one half of a degenerate joint mode.
  **Any fit using `random_q = TRUE` moves**; every other fit is unchanged, since
  the parameter stays mapped without the flag. `index_q_log_sd`, the prior sd the
  assessor sets on q itself, is still never estimated.

  How well the sd is informed depends on the series. On a 40-year index with the
  observation sd held fixed and q deviations injected at a true sd of 0.4, the
  estimate reaches its lower bound when the observation sd is 0.1 or larger, and
  recovers 0.06 at 0.05 -- i.e. a short or noisy index can return a value that
  reads as a constant q. That is visible in the estimate and its gradient, so
  check both rather than assuming the sd is informed.

## Minor improvements

* **A shared selectivity or catchability group now says that no fleet keeps its
  own deviation sd.** Fleets sharing a `Selectivity_index` estimate one
  `sel_dev_log_sd` between them. TMB collapses a shared parameter to the mean of
  its members' starting values, and that sd is held on the log scale, so the
  group starts at the **geometric mean** of their `Time_varying_sel_sd` — a
  value none of the rows asks for, reached in silence. `build_map()` now warns
  once per group, naming its fleets, their differing values and the geometric
  mean the group starts from, and does the same for `Time_varying_q_sd` across a
  shared `Catchability_index`. Both read the finished map, so they cover only
  fleets whose sd is actually estimated: with it fixed, or on an `Off` fleet,
  each fleet keeps its own value and there is nothing to report. No fit changes.

* **`rearrange_data()` no longer iterates once over an empty `age_error`.** The
  loop used `1:nrow()`, which runs for `i = 1` and then `i = 0` on a frame with
  no rows. It is `seq_len()` now. No bundled dataset ships an empty `age_error`,
  so no fit that works today changes.

## Internal

* **`R/0-build_srr_and_M.R` is split into one file per process.** It had reached
  1,497 lines and 52 top-level objects carrying six unrelated constructors, only
  three of which its name described: 612 lines were catchability, selectivity and
  composition linkage machinery. It is now `0-build_srr.R`, `0-build_m1.R`,
  `0-build_growth.R`, `0-build_catchability.R`, `0-build_selectivity.R` and
  `0-build_composition.R`. `.coerce_switch_arg` moves to `0-switches.R` beside
  the other switch coercion helpers, and the three linkage stratum helpers to
  `0-build_linkage.R`. A pure relocation: every one of the 441 top-level object
  bodies is unchanged, and the multiset of non-blank lines across `R/` is
  identical before and after. Exports, `NAMESPACE` and the pkgdown reference
  index are untouched, since all three index topics rather than files.

* **The composition and CAAL row normalization is one helper.** The two blocks
  were the same five lines twice; `.normalise_rows()` replaces both. It divides
  by the row sums instead of going through `t(apply())`, which returned a
  transposed matrix for a single-bin composition; every multi-bin result is
  unchanged. The zero-row guard stays.

* **Four stale or incorrect source comments corrected**, in `build_params()`,
  `clean_data()`, the model's fishing-mortality section, and the two
  shared-parameter blocks in `build_map()` that still claimed `fit_mod()`
  suppresses their warnings. Documented in `inst/dev/CLEANUP_BACKLOG.md`.

# Rceattle 5.13.0

## Bug fixes

* **Reference points ignored an estimated stock-recruit curve unless the
  projection also used it.** Reference-point recruitment is set by two arms:
  one for mean recruitment, requiring `proj_mean_rec = TRUE` *and* no curve, and
  one for the curve, requiring `proj_mean_rec = FALSE`. A model with an
  estimated Beverton-Holt or Ricker and `proj_mean_rec` left at its default of
  `TRUE` matched neither, so `NByage0` and `NByageF` stayed 0 for every year
  after the first and `SB0` was the initial cohort decaying away — on a five-age
  fixture, 3.344 falling to 1.230 over six years instead of holding at 3.344.
  `SB0` in the terminal year is the depletion reference HCR 5 and HCR 6 read, so
  this inflated perceived depletion and the catch advice that follows from it.
  The curve arm now applies whenever `srr_pred_fun` is a stock-recruit form.
  Projection recruitment is unaffected — that switch is read separately, so a
  mean-recruitment projection stays mean. No bundled dataset or golden reference
  reaches the gap, because the default `srr_pred_fun` is mean recruitment.

* **Ten writes to the male slot ran past the end of the sex dimension.** Arrays
  carrying a sex are dimensioned by `max_sex`, so where *every* species is
  one-sex that dimension has length 1 and sex index 1 does not exist. The "males
  take `1 - sex_ratio`" half of each recruitment and sex-ratio assignment wrote
  it unconditionally. The value is 0, but the index is out of range and the
  model is compiled without bounds checking, so it was a silent write into the
  next age of the same species and year. Reproduced on BS2017SS with
  `safebounds = TRUE`; guarded, the same fit is clean at an unchanged objective.
  Affects BS2017SS, BS2017MS, and any model whose species are all one-sex.

* **`comp_data$Sex` was never checked against `nsex`.** `M1_base`, `weight` and
  `ration_data` are all validated; composition was not. A joint row (`Sex = 3`)
  on a one-sex species is the damaging case, because the two registries
  disagree: `check_composition_data()` sizes it at `nages` while the model
  writes it out to `nages * 2`, so the surplus landed in the next observation's
  predicted composition and silently changed that observation's likelihood.
  `Sex = 2` on a one-sex species is the same class — the model reads a sex
  index that does not exist. Both are now refused with a message naming the
  species. `Sex = 0` and `Sex = 1` are unaffected.

* **`build_srr(proj_mean_rec =)` was documented backwards.** The roxygen said
  0 = mean recruitment and 1 = the stock-recruit relationship; the model
  reads 1 as mean recruitment and 0 as the relationship, and the argument
  defaults to `TRUE`. Anyone who set it from the documentation got the opposite
  of what they asked for. Corrected, with the reference-point behaviour stated.

* **Recruitment at `minage > 1` read memory instead of spawning biomass.**
  Recruits arriving at age `minage` in year `yr` were spawned in `yr - minage`,
  so for the first `minage - 1` hindcast years the stock-recruit relationship
  indexed `ssb()` at a negative column. Eigen does not bounds-check in a release
  build, so this did not error -- it returned whatever was in adjacent memory.
  On a `nages = 5` Beverton-Holt fixture, recruitment in those years came back
  as `8.6e-314`, `is.finite()` reported TRUE, and the objective moved from
  14407.38 (`minage = 1`) to 15162.60 (`2`) and 17532.15 (`3`) purely from the
  values read. It was only reachable with an active stock-recruit relationship:
  `srr_fun = 0` never calls the relationship with a spawning biomass.

  Those years now take `R_init`, equilibrium recruitment at `F = Finit`, with
  the recruitment deviation already estimated for them — the same
  `R_init * exp(rec_dev)` the first year is built from, so the whole
  pre-start-spawned block rests on one equilibrium. This follows Stock
  Synthesis, which covers the pre-start period with an equilibrium age
  composition plus early recruitment deviations rather than extending the
  relationship backwards; WHAM avoids the boundary by fixing the recruitment lag
  at one year regardless of recruit age. Not `R0`, which in a later year is
  still the starting value. On the same fixture the objective now moves
  14407.38 / 13959.97 / 13587.67 across `minage` 1/2/3.

  All four sites are guarded. Expected recruitment takes the same `R_init`
  anchor, so the stock-recruit penalty scores the deviation alone — there is no
  spawning biomass to carry information about the curve's level. The reference
  points keep the relationship and only clamp the lag into range.

  **No fit anyone runs today changes**: every bundled dataset and all three live
  assessments are `minage = 1`, where the guard cannot fire. The four
  golden-reference models are unmoved.

* **`CAAL_distribution = "MultinomialAFSC"` now works.** `CAAL_distribution`
  shares `comp_loglike_map` with `Comp_distribution`, so the value passed
  `validate_switches()` and then died inside the model with
  `Invalid 'caal_ll_type'` -- the CAAL dispatch had only cases 0 and 1. The AFSC
  multinomial pseudo-likelihood is now implemented for conditional age-at-length
  in the same form already used for age composition, reading the CAAL row's
  proportions with its own sample size and weight. As for age comps it is a
  pseudo-likelihood rather than a density, so OSA residuals are taken under the
  full multinomial and the simulation draw is multinomial.

* **A model with `nlengths < nages` wrote past the end of two matrices.**
  `age_hat` and `age_obs_hat` are indexed by age -- up to `nages * 2` for
  joint-sex composition -- but were sized from `comp_obs`, whose width is the
  number of `Comp_` columns in the workbook. For a model whose composition data
  are all length-based with fewer length bins than ages, that width is
  `nlengths` and the writes ran off the end. Unchecked in a release build, so
  the effect was a silent write into adjacent memory rather than a crash. Both
  are now sized from the widest age index the model's own composition rows will
  be written at -- `nages`, or `nages * 2` for a joint-sex row -- and never
  narrower than the observations. Sizing them at `nages * 2` unconditionally
  would also close the overrun, but `age_hat` is REPORTed and assessment scripts
  read it, so a single-sex model must not gain trailing all-zero columns it
  never had.

## MSE

* **A scalar `cap` is now an annual ceiling, not an assessment-interval one.**
  `run_mse()` filled catch for every year of the assessment interval and then
  tested the cap against the sum over all of them, so at
  `assessment_period = 2` a one-year cap was applied to a two-year total,
  roughly halving projected catch. The cap is now applied within each projection
  year.

  The same line carried a second and larger defect, so **any stored MSE result
  using a scalar `cap` changes, at every `assessment_period` including 1, and
  whether or not the cap ever bound**. It was
  `ifelse(sum(x) > cap, cap * x / sum(x), x)` with a length-1 test, and
  `ifelse()` returns the shape of its test -- so whichever branch was taken was
  truncated to its first element and then recycled across every row. Two species
  landing 80 t and 20 t came out as:

  | scalar `cap` | old | new |
  |---|---|---|
  | 50 (binds)         | 40 t, 40 t -- 80 t total, over a 50 t ceiling | 40 t, 10 t |
  | 500 (never binds)  | 80 t, 80 t -- 160 t total, catch invented     | 80 t, 20 t |

  The non-binding case is the more serious one: a cap set generously enough to
  be inert still rewrote every species' catch to the first species' value. Only
  an exactly equal split across species was unaffected. The species-specific
  vector form of `cap` was always per row and is untouched.

* **`mse_summary()` scored HCR 2 against a hardcoded depletion.** The estimation
  model's `P(SSB < SSBlimit)` used a literal `0.5 * 0.35` under HCR 2, on the
  depletion scale, while the operating-model arm scored the same rule as an
  absolute `0.5 * SBF` -- so the EM/OM cross-tab compared two different
  criteria. In single-species mode the EM now reads `ssb_limit_thresh()`, the
  helper the OM already uses, so both sides apply one rule and it cannot drift.
  `Plimit` is not that rule: it defaults to 0, which would report a default
  single-species ConstantF run as never overfished.

  Under `msmMode > 0` both sides still fall through to `Plimit`, so a
  multispecies ConstantF run left at the default `Plimit = 0` reports
  `P(SSB < SSBlimit)` as 0. That is the operating model's own long-standing
  multispecies rule and the two sides agree, so it is left alone here -- but set
  `Plimit` explicitly on a multispecies ConstantF run.

* **A reproducibility boundary, for two kinds of run.** A seeded `run_mse()`
  from 5.12.0 does not reproduce here if it used a scalar `cap` — the rewrite
  above changes the catch taken, cap binding or not — or if its control rule
  read reference points from an estimated stock-recruit curve, which the
  `NByage0` / `SB0` fix above moves. A run with neither is unaffected by this
  release. That is unlike 5.9.0 and 5.18.0, which move the random stream for
  every seeded run; `vignette("hcrs-and-mses")` lists all three.

## Other changes

* **`fit_mod()` no longer suppresses every warning `build_map()` raises.** The
  suppression hid messages that exist to be seen -- an M1 model asking for
  sex-specific mortality on a single-sex species, a `Time_varying_sel` the
  selectivity form ignores, a Laplace request the map cannot honour. Each
  changes what is estimated. They are now de-duplicated and passed through, so
  a retrospective or MSE that re-enters `fit_mod()` once per peel or simulation
  prints each distinct warning once rather than hundreds of times.

* **`random_sel = TRUE` is now refused on a fleet with `Time_varying_sel =
  "Block"`,** instead of producing a NaN objective. Selectivity blocks are fixed
  effects — one slope and inflection per block, with no penalty — but
  `fit_mod()` added the block parameters to TMB's `random` regardless, and the
  model scores selectivity deviates only under `"IID"`, `"RandomWalk"` and
  `"RandomWalkAscending"`. They were integrated against no density, with
  `sel_dev_log_sd` mapped out so there was no variance to estimate either, and
  the fit died with TMB's `NA/NaN gradient evaluation` — which names neither
  selectivity nor `random_sel`. The error now says which fleets and what to
  change.

  No bundled dataset or reference model reaches this, and no fit that works
  today changes. `"Block"` remains fully supported with `random_sel = FALSE`.

# Rceattle 5.12.0

## Breaking changes

* **A Beverton-Holt or Ricker stock-recruit curve can no longer be combined
  with `msmMode > 0`.** `data_check()` now refuses it. Those curves anchor
  steepness, `R0` and `R_init` on spawning biomass per recruit, and SPR is not
  defined under predation: total mortality carries `M2`, which scales with
  predator abundance, so per-recruit spawning output is not a property of the
  prey stock alone. The model has always declined to compute it under
  multispecies mode -- but the code that consumes it was never gated to match,
  so `SPR0 = 0` reached it. With `srr_fun >= 2` that is a division by zero:
  `R0` and `R_init` came back `-Inf` and the objective `NaN`. With the Ianelli
  configuration (`srr_fun` mean, `srr_pred_fun` Beverton-Holt or Ricker) it was
  worse, because **the fit ran to completion** and reported a finite objective
  with `steepness` silently 0 and `R_hat` `-Inf`.

  Nothing that produced a usable fit is refused. Use
  `build_srr(srr_fun = "mean", srr_pred_fun = "mean")` under predation, or fit
  the stock-recruit curve in single-species mode.

* **Three `fleet_control` columns are now validated: `Selectivity_dimension`,
  `Sel_shape_dir` and `Sel_shape_mode`.** A typo in any of them previously
  resolved to `NA` rather than erroring -- `Selectivity_dimension` became a
  missing selectivity dimension, the two `Sel_shape_*` columns a missing penalty
  mode -- and nothing downstream re-checked it, so the model fitted and reported
  a number on an input nobody had accepted. A workbook carrying an invalid value
  will now refuse to load.

  Each column is checked against exactly what its consumer implements, so no
  working spelling is refused: `Sel_shape_mode` accepts `"Smooth"`, `"smooth"`
  and `1`; `Sel_shape_dir` accepts `"Increasing"`, `"increasing"` and `"-1"`
  (the ADMB sign convention). `Selectivity_dimension` accepts only `"Age"` and
  `"Length"`, the only values `rearrange_data()` matches.

  A blank `Selectivity_dimension` cell takes the schema default (`"Age"`) rather
  than erroring, so the partial-assignment idiom keeps working, and on a model
  that estimates growth the fill names the fleets it applied to. **This can move
  a fit** on a model that carried blanks: those fleets previously reached the
  model as a missing dimension, which sizes selectivity by `nlengths` rather
  than `nages` and so changes a max-normalized curve. No bundled dataset or
  golden model carries the column at all.

  Pre-flighted over every workbook in the ecosystem: **196** carry a
  `fleet_control` sheet and not one value the new check rejects.

* **`write_data()` writes `fleet_control` in schema column order**, as the
  control and bioenergetics sheets already did. Values are unchanged and
  round-tripping is identical; only the column order in the workbook moves.
  The schema order was also tidied so the workbook reads sensibly: `Month` now
  sits with the fleet identity columns rather than behind 26 selectivity
  columns, and `Index_distribution` heads the index-observation block with the
  standard-error settings it governs.

* **`Catchability = "AR1"` (the QAR1 form of Rogers et al. 2024) is now an
  error.** It never worked. `build_map()` gates the log-q deviates on
  `Time_varying_q %in% c("IID", "AR1", "RandomWalk")`, but under this form
  `Time_varying_q` holds an `env_data` column index rather than a mode -- so
  the deviates were never estimated and q came back constant. The fit ran and
  reported a time-invariant catchability where a time-varying one was asked
  for, silently.

  It errors rather than warns because a warned fit still returns a `summary()`
  that looks ordinary, and nothing downstream can tell that its q is constant or
  its objective inflated. On BS2017SS fleet 7 the `Catchability deviates`
  likelihood row accumulates **54.8** from deviates that are identically zero,
  so the reported objective is not comparable with any other model's, and
  `index_q_dev_log_sd` is left free with nothing opposing it, driving its sigma
  to zero.

  The Rogers form is available as a linkage. `observe` is what makes it QAR1
  rather than a free AR1 on q — it names the environmental series the deviates
  are observed against — and `obs_sd` is that series' measurement SD, which the
  legacy switch never asked for and must now be supplied:

  ```r
  build_catchability(linkages = list(q = linkage_spec(
    ~ ar1(1 | Year), by = ~ fleet, fleet = <Fleet_code>,
    observe = "<that fleet's env_data column>", obs_sd = <its measurement SD>)))
  ```

  and pass it to `fit_mod(qFun = )`. `GOA pollock/2025/04-fit-and-diagnostics.R`
  in `../Rceattle-models` is a worked example: it runs exactly this form.

  The schema marks the code removed too, so `?BS2017SS` and the workbook meta
  sheet say so before a model is built rather than only at fit time.

  Note this is **not** `Time_varying_q = "AR1"`, a different switch sharing the
  same string, which puts an AR1 structure on an ordinary `"Estimated"` q and
  works correctly.

## Bug fixes

* **`Estimate_catch_sd = "Analytical"` now works.** The option was documented in
  the schema, accepted by `validate_switches()`, and treated as a real mode by
  `build_map()` (which maps `catch_log_sd` out) -- but the model's dispatch
  had only cases 0 and 1, so a model using it passed every R-side check and then
  died inside the fit with `Invalid 'Estimate_sigma_catch'`. It is now
  implemented, mirroring `Estimate_index_sd = "Analytical"`, and the fitted
  value is reported per fleet as `catch_analytical_sd`.

  Note this concentrates sigma out of the catch likelihood, so the catch data
  no longer pins F as tightly as a fixed sd does, and the objective is not
  comparable to a fixed-sd fit. A model that does not ask for it is unaffected.

* **Both analytical observation sds now account for the lognormal bias
  adjustment**, so `Estimate_index_sd = "Analytical"` changes fits that use it.
  Rceattle fits an index or catch series to a mean-unbiased prediction,
  `log(obs) ~ N(log(pred) - b*sigma^2/2, sigma)` with `b = bias_adjust_obs`
  (default 1). The Ludwig and Walters (1994) form
  `sigma = sqrt(mean((log(obs) - log(pred))^2))` is the sd that minimises that
  density only when `b = 0`, so on a default fit the option did not return the
  concentrated estimate it advertised. The sd used is now

  ```
  sigma^2 = 2*S / (sqrt(1 + b^2*S) + 1),   S = mean((log(obs) - log(pred))^2)
  ```

  which reduces to the Ludwig-Walters form at `b = 0`. On the three `BS2017SS`
  fisheries the old expression sat 4-17% high and left up to 1.6 units of
  negative log-likelihood on the table. No bundled dataset, reference model or
  live assessment sets either switch to `"Analytical"`, so no fit anyone runs
  today moves; a model that does use it will. The one workbook in the ecosystem
  that sets it (`Pacific hake/Data/Bridging/hake_bridging6_2022.xlsx`, read only
  by a deprecated script) returned a `NaN` objective before this change too --
  its second index fleet has a zero predicted index in every hindcast year.

* **A fishery asking for an analytical catch sd with nothing to estimate it
  from is now refused by `data_check()`**, as is the index equivalent. With no
  fitted positive observation the sd is undefined, and it used to fall through
  as 0 -- which the likelihood never reads, but the reported `index_sd` /
  `catch_sd` and the `sim_mod()` draw both do, and `rnorm(mean, 0)` is a
  deterministic observation. A zero-catch year is legal input, so this is
  reachable on the catch side; on the index side `data_check()` already refused
  a non-positive observation, and the matching guard added to the index
  estimator is parity rather than a live fix.

* **`write_data()` no longer fails on a workbook that predates a control or
  bioenergetics switch.** Both sheets were assembled from the full schema
  object list -- the control sheet by `rbind()`, which drops a `NULL`, and the
  bioenergetics sheet by hard-coded row index into a fixed-height matrix. An
  absent object therefore aborted the whole write with `arguments imply
  differing number of rows` or `number of items to replace is not a multiple of
  replacement length`. The live Pacific hake workbook hits this: it predates
  `alpha_wt_len` / `beta_wt_len`, so `read_data()` then `write_data()` could not
  round-trip it. Both sheets now write the objects the `data_list` actually
  carries, in schema order, as `fleet_control` already did.

  Absent objects are dropped rather than written at their schema default: the
  default belongs to the model, and baking one into a workbook would turn a
  value `switch_check()` announces at fit time into one the file asserts.
  Keying the bioenergetics rows by name also retires the duplicate row-order
  list that had to be kept in sync with the schema by hand.

* **`switch_check()` accepted the selectivity spelling `"Non-parametric"` while
  `validate_switches()` rejected it**, so a model written that way loaded, was
  normalised to nothing, and then failed its own data check. It is now upgraded
  to `"NonParametric"` on the way in.

## New features

* **`write_template()` writes every column the schema defines.** The
  hand-written list had fallen 14 columns behind, and a column missing from the
  template is how a user never learns an option exists. The penalty weights
  `Sel_curve_pen1/2/3` are deliberately left blank rather than seeded with their
  schema default of 0: `switch_check()` converts `Sel_shape_sd` /
  `Sel_curvature_sd` / `Sel_devmag_sd` into a weight only where the raw weight
  is still blank, so a seeded 0 would silently disable that interface and leave
  a non-parametric fit with no shape penalty at all.

* **The model-level switches have one table**, `.rce_model_switch_schema()`, and
  the comments `save_config()` writes are projected from it. `estimateMode`'s
  `DebugOptimize` was a real mode the comments never mentioned, and `msmMode`
  was described in prose naming none of its canonical values. Each comment now
  gives the integer code beside the readable name
  (`SingleSpecies (0) / MSVPA (1) / TypeIIIMSVPA (2)`), since a saved config
  writes the switch as the number.

## Minor improvements

* **`quantities$NbyageSPR` now carries dimension names.** Its first dimension is
  a reference point rather than a species -- `c("F0", "Flimit", "Ftarget",
  "Finit")` -- so reading it no longer means knowing that slot 3 is `Finit`.
  Species and age bins are labelled to match the other numbers-at-age arrays.
  Positional indexing is unchanged.

# Rceattle 5.11.1

## Documentation

* **The vignettes can now be executed.** They still render without running their
  code by default -- several chunks fit real models -- but
  `RCEATTLE_EVAL_VIGNETTES=true` turns evaluation on, and a new weekly
  `vignettes.yaml` CI job runs them for real. Nothing executed them before, so an
  API break in a vignette was invisible until a user copied the code. That is how
  `hcrs-and-mses.Rmd` kept calling `knitr::kable()` on `mse_summary()` for several
  releases after it started returning a list; that code is fixed, and the
  per-entity structure (`$species`, `$fleet`, `$total`, `$meta`) is now explained.

* **The switch tables document every value the model accepts.** `Selectivity`
  stopped at 7, omitting `DoubleNormal` (8), `NonParametricPM` (9) and
  `LogisticPM` (11); `Catchability` stopped at 6, omitting `AnalyticalArith` (7).
  `PowerEquation` is now marked as not implemented, and the 5.8.1 deprecation of
  `"Environmental"` (deprecated in 4.9.0) is recorded with the linkage that replaces it.

* **`Index_distribution` is documented.** It had no vignette coverage at all,
  despite deciding whether an index is scored on the log or the natural scale --
  and despite a second, hand-synced registry (`.index_rows_natural_scale()`) that
  a new natural-scale family must also be added to, or it silently gets the
  log-scale residual formula.

* `reweight_comps(fleets =)` is documented, as is the log scale `Comp_weights`
  takes under a Dirichlet-multinomial. `fit_mod()`'s `initMode` moves from an
  849-character `@param` into an **Initial age structure** section, and
  `suitMode` / `avgnMode` now say "declared but not implemented" in one voice
  rather than three. Examples added to `build_srr()`, `build_hcr()`,
  `mse_summary()`, `sim_mod()`, `model_average()`, `osa_residuals()` and
  `run_mse()`.

* The vignettes and `build_data()`'s help no longer tell users to call
  `data_check()`, which is internal; they point at `build_data(.check = TRUE)`.

## Internal

* Removed `flt_sel_ind`. `rearrange_data()` computed it from `Fleet_code` on
  every fit and nothing read it -- it was declared in no `DATA_` object and
  referenced nowhere in the package or the assessment repos.

* **The C++ dispatch branches are pinned to the R maps that select them**
  (`test-schema-cpp-dispatch.R`). The pinned exemptions are a machine-checked
  inventory of what is implemented in the model but unreachable from R
  (non-parametric growth; `msmMode` 3-9) and what R can encode but the dispatch
  never sees. The schema's `tmb_target` field, previously set on no rows, now
  carries the rename for 27 columns and is checked against both ends of it.

* `R/data.R`'s field dictionary is pinned against the switch maps: a documented
  code the map does not define, or a map value the dictionary omits, now fails.
  It also may no longer present a deprecated alias as a canonical column.

* `.rce_config_schema()` reads `estimateMode_map` and `msmMode_map` instead of
  hardcoding them, so the comments `save_config()` writes list every valid value.
  `estimateMode` was missing `DebugOptimize`.

# Rceattle 5.11.0

## New features

* **The diagnostics all take the fitted model as `object`.** The same argument
  was previously spelled `Rceattle` in `retrospective()`, `jitter()`,
  `self_test()`, `sim_mod()`, `sample_rec()`, `remove_F()` and
  `model_average()`, and `fit` in `reweight_comps()`, `osa_residuals()` and
  `process_residuals()` -- so which name to use depended on which diagnostic you
  reached for, and `Rceattle` collided with the package name
  (`Rceattle::retrospective(Rceattle = mod)`). All ten now take `object`, matching
  `convergence_diagnostics()` and R's own `summary()`/`coef()`/`vcov()` methods.

  **Existing scripts keep working unchanged.** The old names are still accepted
  and are silent in this release; they begin warning in a future version. To find
  your own call sites early, set
  `options(Rceattle.warn_deprecated_args = TRUE)`. Supplying both names for one
  model is an error rather than a guess.

  Positional calls are unaffected: the old name is the last formal on every one
  of the ten, and the new `phase` and `fit_control` arguments below are appended
  after the arguments that predate them rather than inserted among them.
  `profile()` keeps `fitted` because it is an S3 method for `stats::profile()`,
  and `mse_summary(mse =)` and `compare_sim(operating_mod =)` keep their names
  because they take an MSE result and a simulation set, not a fitted model.

* **`sim_mod()`, `sample_rec()`, `remove_F()` and `model_average()` now say what
  is wrong when they are given something that is not a fitted model.** They
  previously failed further in with `argument is of length zero` or
  `invalid 'type' (list) of argument`.

* **`fit_mod()` reports a composition weight of 1 on a Dirichlet-multinomial
  fleet.** That likelihood reads `Comp_weights`,
  `CAAL_weights` and `Diet_comp_weights` as the LOG of the starting weight, so
  the value `write_template()` seeds (1) is a starting weight of e, and a weight
  of 1 is written as 0. A model built from the model and switched to a
  Dirichlet-multinomial previously started at e with nothing saying so. Only an
  exact 1 is reported, so a deliberate value is left alone, and Off fleets and
  fleets carrying no composition data are skipped. It fires once per fit, and
  the diagnostic refits suppress it the way they already suppress
  `data_check()`. No value and no fit changes -- this is a message.

* **`retrospective()`, `jitter()` and `self_test()` accept
  `fit_control = fit_control(...)`,** matching `fit_mod()`. Only `phase` and
  `getsd` are read, because those are what these diagnostics forward to each
  refit; setting any other field is an error rather than a silent no-op.
  Supplying it replaces the defaults they otherwise infer from the fitted model.
  Omitting it -- the default -- leaves every existing call behaving exactly as
  before.

  Which fields you asked for is read from what you **set** -- named in the
  `fit_control()` call, or assigned to afterwards (`ctl$getsd <- TRUE`) -- not
  from whether the value differs from a default. So `fit_control(getsd = TRUE)`
  and `fit_control(phase = FALSE)` do what they say even though `TRUE` and
  `FALSE` are those fields' defaults. A field you never touch keeps the
  diagnostic's own default, which is not always `fit_control()`'s:
  `retrospective()` phases its peels where `fit_control()` does not, so
  `fit_control(getsd = FALSE)` asks about standard errors and cannot also
  flatten Mohn's rho.

* **`?rceattle-refit-args` documents the vocabulary `retrospective()`,
  `jitter()` and `self_test()` share** -- `object`, `cores`, `fit_control`, and
  what `fit_control()` reaches through a refit -- the way
  `?rceattle-plot-args` already does for the plotters. `phase`, `getsd` and
  `timeout` stay on each function, because their defaults and what they mean for
  that diagnostic differ.

* **`retrospective()` takes `phase`,** the value it previously fixed internally.
  The default is `TRUE`, which is what it always used: a peel restarts from the
  unpeeled fit's starting values with a year removed, so without phasing the
  parameters barely move, the peels sit on top of the full model, and Mohn's rho
  is biased towards zero. No peel changes unless you set it.

* **`?jitter` now records that attaching Rceattle masks `base::jitter()`.** The
  function keeps its name; call `base::jitter()` for the base-graphics
  behaviour.

## Bug fixes

* **A parallel `retrospective()`, `jitter()` or `self_test()` called under the
  deprecated `Rceattle =` name no longer sends the fitted model to each Windows
  worker twice.** The PSOCK path exports the caller's whole frame, and both
  argument names were bound to the same model. macOS and Linux use FORK and were
  never affected. No result changes -- this is transfer time and worker memory.

## Deprecations

* `Rceattle` and `fit` as the fitted-model argument of the ten diagnostics
  above. Accepted silently now, warning from 5.13.0, removed in 6.0.0. (5.11.0 and
  5.12.0 are both unreleased on this line, so the silent grace period has not
  yet reached a user; the warning moves with it.)
* `rearrange_dat()` now names its removal version (6.0.0) rather than
  deprecating open-endedly. Use `rearrange_data()`.

# Rceattle 5.10.0

## New features

* **The timeseries, predation and selectivity plotters now share one set of
  arguments, and use them.** Only `plot_timeseries()` ever honoured `line_col`,
  `lwd`, `lty` and `alpha`; the others declared them and ignored them, so
  colours and line widths silently did nothing. They are now resolved in one
  place, documented once in `?"rceattle-plot-args"`, and applied consistently.
  `line_col` accepts colour names, hex codes, or base-graphics palette indices
  (`line_col = 1`), and supplies the palette for whichever variable the figure
  maps to colour; on `plot_selectivity()`'s year fan it gives the ramp anchors
  instead. `lwd` keeps the base-graphics scale, where the default `3` is a
  standard-weight line. The remaining plotters -- `plot_mortality()`,
  `plot_maturity()`, `plot_comp()`, `plot_data()`, `plot_stock_recruit()`, the
  index, catch and diet families -- still take their own arguments.

* **The predation plotters honour the shared arguments.** `plot_b_eaten()`,
  `plot_b_eaten_prop()`, `plot_m_at_age()`, `plot_m2_at_age_prop()` and
  `plot_ration()` now use `line_col`, `lwd`, `lty`, `minyr`, `maxyr` and
  `incl_mean` (and `alpha`, where the figure has a ribbon), all of which they
  previously declared and ignored. They also accept `maxyr`, `lty`, `incl_mean`
  and `top_adj` where those were missing entirely, so scripts passing them no
  longer stop with `unused argument`. `species` and `spnames` worked before and
  now additionally take names, a logical mask and `"all"`, and validate.

  `line_col` and `lty` follow the figure, not the model: colour separates
  predators in `plot_b_eaten_prop()` and `plot_m2_at_age_prop()`, and line type
  separates the sexes in `plot_ration()` and `plot_m_at_age()`. Each function's
  help says which. Too few colours are recycled, now with a warning naming what
  they coloured, and a varying `lty` whose key has one level warns rather than
  being dropped in silence.

* **`add_ci = TRUE` says when it cannot draw an interval.** None of the
  predation quantities carry standard errors -- `M_at_age` and
  `B_eaten_as_prey` are `REPORT`ed but not `ADREPORT`ed, and consumption and the
  M2 proportions are products and ratios of such series -- so the argument was
  silently doing nothing. It now warns once and draws no ribbon.

* **`plot_selectivity()` draws every model, on the right dimension, and uses its
  arguments.** It previously read only the first fit, so a list of models
  silently lost all but one; `model_names`, `line_col`, `lwd` and `species` were
  declared and ignored; and `species` defaulted to three hard-coded Bering Sea
  species names.

  It now takes `line_col`, `lwd`, `lty` (which separates the sexes), `species`,
  `spnames`, `minyr`, `maxyr` and `alpha`. With one model colour is still the
  year, so the figure is unchanged apart from line width; with several, colour
  separates the models and the year fan moves to transparency, keeping the
  curves superimposed for comparison. `colour_by` forces either.

* **Length-based fleets are drawn on length bins.** `plot_selectivity()` read
  `sel_at_age` for every fleet and labelled the axis "Age", so a fleet whose
  `Selectivity_dimension` is `"Length"` showed the growth-matrix conversion of
  its curve rather than the curve that was fitted. Each fleet is now drawn on its
  own dimension, and a model mixing the two returns one figure per dimension (a
  named list) rather than putting ages and length bins on one axis.

* **`species` and `spnames` mean the same thing across these plotters.**
  `species` selects -- by index, name, logical mask, or `"all"`, in the order
  given -- and `spnames` labels. Several plotters previously read `species` as
  display labels; a character vector giving one label per species is still read
  that way, with a message. `plot_timeseries()` and its wrappers
  (`plot_biomass()`, `plot_ssb()`, `plot_recruitment()`, the depletions,
  `plot_exploitable_biomass()`, `plot_f()`) gained selection by name.

* **`data_check()` warns when `diet_data` does not cover every age of an
  empirical-suitability predator.** Under `suitMode = 0` suitability is read
  straight out of the diet data, so an age with no diet row is switched off
  rather than estimated: a predator age with no rows gets `suit_other = 1` and
  exerts no predation mortality, and a prey age with no rows is never eaten.
  Neither raises an error or moves the likelihood, so a diet table truncated at
  the wrong age silently drops part of the predation. The warning names the
  species, sex and age range for each role.

  Only prey-at-age-in-predator-at-age rows count toward coverage: the
  aggregated diet formats feed the diet likelihood but are skipped when the
  suitability array is filled, so they cannot close a gap. Parametric
  suitability (`suitMode > 0`) and single-species models are unaffected. A
  species with no prey rows at any age is not reported — nothing eating it is a
  modelling choice, not a truncated table.

## Bug fixes -- figures whose numbers change

Three predation plotters drew quantities that did not match their axis labels.
All are corrected, so figures regenerated from them will differ from earlier
runs of the same model.

* **`plot_m2_at_age_prop()` draws a share, not a contribution.** `M2_prop` holds
  each predator's contribution to M2, which sums over predators to `M2_at_age`,
  so the plotted "proportion" reached 1564 on `BS2017MS`. The contributions are
  now divided by their total, giving shares in [0, 1] that sum to 1 across
  predators for each prey age and year; a prey age with no predation in a year
  leaves them undefined and draws nothing. The y axis reads "Share of M2 at age
  `<age>` by predator".

* **`plot_ration()` multiplies the ration by average numbers-at-age, not
  biomass-at-age.** `consumption_at_age` is one fish's annual ration in kg and
  numbers-at-age are in thousands, so the product is mt -- the way the model
  forms total consumption (`avgN_at_age * ration`, `predation.hpp`). Multiplying
  by biomass instead weighted the age-sum by weight-at-age, so "million mt"
  described nothing it computed. Average numbers rather than start-of-year
  numbers, so the series reconciles with `plot_b_eaten()`; under the default
  `avgnMode = 0`, `N_at_age` would overstate it by `1 / ((1 - exp(-Z)) / Z)`. On
  a fitted `BS2017MS` the first year drops 38.3% for pollock, 20.1% for cod and
  13.5% for arrowtooth flounder.

* **`plot_b_eaten()` is in million mt.** It plotted `B_eaten_as_prey` in the mt
  the model reports it in, while `plot_b_eaten_prop()` -- the same quantity
  broken down by predator -- was in million mt, so the two could not be read
  side by side. Both are now in million mt, the display unit the timeseries
  plotters use. `p$data` moves by the same factor of 1e6.

## Bug fixes

* **The `diet_data` age-bound check compares each species against its own
  `nages`.** It matched them by position, so a species missing from the `Pred`
  or `Prey` column shifted every later species onto the wrong age limit and
  could reject a valid table. Out-of-range ages are still caught.

* `model_names` given as a `list()` works again, in every plotter that labels
  models. The package's own vignettes build it that way, and a list produced a
  one-element list per model that the plot frame could not bind. Supplying fewer
  names than models now warns instead of silently drawing two models as one
  series.

* A fleet with `Selectivity = "Fixed"` is drawn on ages whatever
  `Selectivity_dimension` says. Empirical selectivity is read into `sel_at_age`
  only, so such a fleet on a `"Length"` dimension would have been drawn as an
  identically-zero curve -- and scripts commonly set `Selectivity_dimension`
  across every fleet at once.

* The projection divider is drawn at the latest hindcast year across the models
  plotted, not at whichever model came last in the list. On a retrospective peel
  the peels end in different years, so the divider could land mid-hindcast and
  label real data as projection. Figures overlaying models with different `endyr`
  -- including any using `reference =` -- will show the divider in a new place.
  It is also omitted when the last hindcast year falls outside a `minyr`/`maxyr`
  window, rather than stretching the axis back to it.

* An invalid `line_col`, `lwd` or `alpha` now stops with a message naming the
  argument. Previously an `NA` colour or width drew nothing, and a transparency
  outside `[0, 1]` saturated the ribbon or errored inside the device, in both
  cases with no explanation.

* A `spnames` of the wrong length now stops rather than recycling, which had
  labelled one species with another's name. A `species` string matching no
  species now stops rather than silently plotting everything.

* `minyr` and `maxyr` narrow the data rather than only the axis, so a panel with
  a free y scale rescales to the window. Clipping the axis alone left the scale
  trained on the hidden years, which on a series spanning orders of magnitude
  squeezed the requested window into the bottom of the panel. They also now work
  on the predation plotters, which declared them and ignored them, and
  `plot_timeseries(save = TRUE)` writes the same window it plots.

* `plot_f()` keys its Ftarget and Flimit reference lines to the species it
  drew. It indexed `Ftarget` and the facet labels with the raw `species`
  argument, which works for indices but gives an `NA` facet key for a name --
  on the same argument that newly accepts names.

* `plot_depletionSSB()` draws the Ptarget and Plimit lines of **one** model, per
  species. The two were collected into a models-by-species matrix and then
  subset with the species indices, which flattens column-major: on a two-model
  overlay whose models carry different `Ptarget`, species 2's line came from
  model 2's species 1. Models in one figure normally share their reference
  points, and then every row of that matrix is identical, which is why it went
  unseen. The values now come from the first model.

* A single `lty` reaches the figures that map line type themselves --
  `plot_ration()`, `plot_m_at_age()`, `plot_b_eaten_prop()`,
  `plot_m2_at_age_prop()` and `plot_selectivity()`. It is applied to every level
  of whatever the figure keys line type on, and warns when that key has more
  than one level, since they are then drawn alike. The default `lty = 1` leaves
  the figure's own line types alone.

* `incl_mean` averages each model over **its own** hindcast, not the first
  model's. On a retrospective peel the answer otherwise depended on the order of
  the list.

* `plot_m2_at_age_prop()` orders its prey panels the way `species` asked for,
  like the other plotters, instead of alphabetically.

* `plot_b_eaten(mse = TRUE)` draws its projection ribbon at `alpha` (default
  0.4) rather than a fixed 0.3, honours `incl_mean`, validates `lwd` and `lty`
  as the other paths do, and says that `line_col` does not apply when the
  simulations are summarized into one band.

* `plot_timeseries(save = TRUE)` writes the years alongside the values, names
  the columns after the models, names the file after the species, and writes
  only the species plotted. It previously wrote unlabelled columns for every
  species, with no year column.

* The pkgdown site builds again. `simulate.Rceattle` (added in 5.9.0) was not
  in `_pkgdown.yml`, and pkgdown stops on a documented topic missing from its
  reference index. `?"rceattle-plot-args"` is listed there too, so the topic the
  plotter help and the vignettes point at has a page on the site.

## Behavior changes

* **`age`, `minage` and `maxage` are ages, not age-bin indices.**
  `plot_m_at_age(age =)`, `plot_m2_at_age_prop(age =)`, `plot_ration(minage =)`
  and `plot_mortality(maxage =)` were passed to the at-age arrays as subscripts.
  The arrays are indexed 1 to `nages`, while a species' ages run `minage` to
  `minage + nages - 1`, so on any species with `minage != 1` the figure drew a
  different age from the one its axis named -- and an age past the plus group
  ran off the end of the array rather than being reported. Each is now resolved
  against each species' own age vector, which is also what lets one figure hold
  species with different age ranges. Only a model with `minage != 1` changes:
  every bundled dataset is `minage = 1`, where the two readings coincide.

  A species that has no such age is dropped with a warning naming it, rather
  than drawn at a shifted age; an age no species carries is an error. Requesting
  an age past the plus group of some but not all species therefore now draws the
  species that have it instead of failing.

* `plot_selectivity()`'s `species` used to be an ignored label argument whose
  default was `c("Walleye pollock", "Pacific cod", "Arrowtooth flounder")`. It
  now selects. Passing species names that only partly match the model's own
  stops with a message naming both readings, rather than guessing -- which is
  what a call copied from the old default does on a model whose species are
  named differently. Pass labels as `spnames`, or the model's own names as
  `species`. A call passing one label per species, none of which match, is still
  read as labels.

* `plot_timeseries(save = TRUE)` needs a `file` stem, and stops without one. It
  previously wrote to a file called `NULL_...csv`.

* `plot_f()` is now built by the same factory as the other timeseries plotters,
  so it takes their full argument list -- it gains `lty`, `save`, `reference`,
  `legend.pos` and `ylab`, none of which it accepted before. `plot_timeseries()`
  gains two internal arguments, `ref_lines` and `suffix`, which the factory uses
  to attach the F and depletion reference points; the reference-point layers are
  consequently later in `p$layers` on the depletion plots than they were. The
  rendered figures are unchanged.

* `plot_selectivity()`'s `p$data` names the x variable `Bin`, not `Age`, and
  carries a `Dimension` column (`"Age"` or `"Length"`); the column holds an age
  or a length-bin ordinal depending on the fleet, so it is no longer named for
  one of them. It also gains a `Model` column, now that every model is drawn.
  `plot_b_eaten()`'s `value` column is in million mt (see above).

## Dependencies

* `ggplot2` now requires >= 3.5.0.

# Rceattle 5.9.0

## New features

* **Every observation type is now simulated by the TMB model, in a `SIMULATE`
  block beside the density that scores it** — survey index (under each fleet's
  own `Index_distribution`, including the correlated `MVN`/`MVNORM` draw), total
  catch, age/length compositions, conditional age-at-length, and, for the first
  time, stomach contents. `sim_mod()` no longer re-implements any observation
  model in R, so the draw and the density cannot be changed on one side only.
  The two copies had diverged: every survey was drawn as independent lognormal
  whatever its `Index_distribution`, and diet was not drawn at all, so a
  multispecies `self_test()` recovered suitability from stomachs that never
  varied.

  Compositions and CAAL are drawn in raw bin space, before tail accumulation and
  before `comp_offset`. Both families are closed under merging categories, so
  drawing raw and letting the refit fold again is exact.

  `Comp_weights` / `CAAL_weights` / `Diet_comp_weights` now enter the draw as
  they enter the density: as an effective sample size for the multinomial
  families (`N * weight`), and through the concentration for the
  Dirichlet-multinomial. **The old R draw used the nominal `Sample_size`
  regardless of the weight**, so a `self_test()` on a re-weighted model was
  handed data more informative than the estimator treats them as being. **Any
  model with a composition weight other than 1 gives a different self-test.**

  A row whose sample size times its weight rounds below one observation comes
  back empty, with `Sample_size` dropped to zero. Rows the model cannot draw —
  a predator under empirical suitability, a covariance fleet outside its fitted
  window — keep their observed values. `sim_mod()` warns for each case.

  Simulated quantities are reported under names ending `_sim` (`catch_obs_sim`,
  `index_obs_sim`, `comp_obs_sim`, `caal_obs_sim`, `diet_obs_sim`), since TMB
  never clears its report environment.

  **This changes results.** A seeded `self_test()` or
  `run_mse(simulate_data = TRUE)` will differ from earlier releases: the draws
  moved into `obj$simulate()` and the random-number stream moved with them.
  Simulating is also slower, since it evaluates the compiled model rather than
  reading `$quantities`.

* **`sim_mod()` can redraw process error as well as observations, via the new
  `process` argument.** Observations are always drawn; process error is a
  choice, because redrawing it changes what `self_test()` measures — from
  recovering parameters to recovering a process. Pass `"recruitment"`, `"M"`,
  `"growth"`, `"catchability"` or `"selectivity"`, the groupings `"dynamics"` or
  `"observation"`, or `TRUE` for all. The default `FALSE` keeps the previous
  behaviour.

  Recruitment covers the initial age structure, since the initial
  numbers-at-age are recruitment from before `styr`. Natural mortality covers
  all six `M1_re` modes. Growth, catchability and selectivity are drawn wherever
  they are written as a random linkage.

  Three cases are not drawn, each with a warning: the equilibrium
  initial-condition modes, which fix `init_dev`; observed AR1 (Rogers QAR1)
  linkages, whose latent state is measured by a covariate series; and the legacy
  `Time_varying_q` / `Time_varying_sel` deviations and the AR1 selectivity
  forms, which carry densities but no draw yet.

  When a process is redrawn, the deviations that generated the data come back as
  `attr(x, "process_sim")`. **Compare estimates against those**, not against the
  source model's fitted deviations, which would report bias by construction.
  Each deviation arrives with a same-shaped `_drawn` logical marking the cells
  the draw touched. `compare_sim()` compares against the operating model and is
  therefore not valid for `process`-drawn replicates.

* **`Index_distribution = "TruncatedNormal"`** fits and simulates a
  natural-scale normal left-truncated at zero, so the density is renormalized as
  `log f(x) = log phi(x; mu, sd) - log Phi(mu/sd)` and the draw is taken by
  inverse CDF on `(0, Inf)`. It is the only natural-scale family whose simulator
  and likelihood are the same distribution: `"Normal"` and the covariance
  families have to redraw the non-positive draws `data_check()` would refuse,
  which samples from a truncated normal while the likelihood scores an
  untruncated one. Prefer it unless an exact ADMB comparison is needed —
  `"Normal"` is unchanged and still reproduces AMAK `avo_like` / `cpue_like`
  term for term.

* **`TruncatedNormal` is registered as a natural-scale family in the
  diagnostics**, so `residuals(type = "pearson")`, `plot_index()`'s observation
  interval and `plot_indexresidual()` give it the natural-scale treatment. A
  test enumerates `index_distribution_map`, so a future family that misses the
  predicate fails in the suite rather than in a residual plot.

* **`simulate()` works on a fitted model**, as the `stats` generic.
  `simulate(fit, nsim = 10)` returns a list of `nsim` data sets — a list at
  `nsim = 1` too, so callers never special-case the length. `seed` follows the
  `stats::simulate()` convention. Use `sim_mod(simulate = FALSE)` for expected
  values rather than draws.

* **`self_test(process = )`** passes through to `sim_mod()`, so a self-test can
  ask whether the estimator recovers a process it was not shown. Each
  replicate's true deviations are returned as `attr(result, "process_sim")[[k]]`.

* **`sim_mod()` warns when a simulated observation is one the model cannot be
  refit on.** A non-finite or negative draw would otherwise be rejected by
  `data_check()` and counted by `self_test()` as a convergence problem rather
  than a data one. The usual cause is a fleet whose observation standard
  deviation never got a value.

## Deprecations

* **`quantities$log_index_sd` and `log_catch_sd` are renamed `index_sd` and
  `catch_sd`.** Neither was ever a log: both hold the observation standard
  deviation the likelihood used for that row. Its scale depends on the family —
  log-scale for `Lognormal`, an absolute sd in the units of the index for a
  natural-scale `Index_distribution`.

  The old spellings are still returned this release and will be removed in the
  next. Fits saved before the rename keep working. The TMB parameter
  `index_log_sd` is a different object and is unchanged.

## Breaking changes

* **A fishery carrying `index_data` now gets an estimable catchability.** The
  model has always fitted such a row, but `build_map_catchability()` entered its
  block only for `Fleet_type == "Survey"`, so a fishery's `index_log_q`,
  `index_q_dev` and `index_log_sd` stayed mapped out and
  `Catchability = "Estimated"` did nothing.

  In reverse, a **survey** with no index rows no longer gets an estimated
  catchability — a q with no index to inform it is a flat direction in the
  likelihood. Sharing takes precedence and is unchanged.

  No bundled dataset is affected, so `/golden-check` is silent. To keep the old
  behaviour set `Catchability = "Fixed"`, or give the fleet its own
  `Catchability_index` to estimate it separately. `print()` marks the case:
  `q: Estimated (fixed: no index data)`.

* **A fishery's index is predicted from the year-averaged numbers.** A survey
  index remains a snapshot at its observation month, `N exp(-(Month/12) Z)`. A
  fishery index is CPUE, which integrates over the year alongside the catch:

  ```
  C/E = q * sum_a sel_a * Nbar_a * w_a,   Nbar_a = N_a (1 - exp(-Z_a)) / Z_a
  ```

  the same mean-numbers term the catch equation and the fishery's age
  composition already use, so a fleet's index, catch and comps are now
  consistent. `Month` is not read for a fishery's index rows.

  This affects trend as well as scale, so catchability cannot absorb it; the
  size over a range of Z is tabulated in
  `vignette("model-options-and-functionality")`. Stock Synthesis splits the same
  way on its per-fleet survey timing.

  Only a fishery carrying `index_data` is affected, so `/golden-check` is
  bit-identical. For a snapshot-timed index on a fishery's selectivity — the
  AMAK/ebswp form the EBS pollock ADMB bridge needs — put the index on its own
  `Survey` fleet and point that fleet's `Selectivity_index` at the fishery.

* **`sim_mod(simulate = TRUE)` no longer works on an averaged model.**
  Simulating evaluates the compiled model, so it needs one. A fit whose `$obj`
  was dropped is rebuilt from its `data_list` and estimates. `model_average()`
  output cannot be rebuilt — its `quantities` are an average over models, so
  there is no parameter vector to draw around. Simulate from one of the
  underlying fits. `sim_mod(simulate = FALSE)` is unaffected.

* **`sim_mod()` errors instead of recycling when `catch_data` no longer lines up
  with the fitted model.** The write-back is a row-position copy, and the old
  draw passed mismatched vectors to `rnorm()`, which recycled the shorter one
  and returned a full-length, wrong answer. The error names both row counts.

* **The `growth_re` switch, `growth_indices`, and the `log_growth_par_devs`
  parameter array are removed.** Nothing consumed `growth_re`: the deviation
  array was mapped off in every configuration and given no density, and the
  parameters its documentation claimed to map were never declared. Setting it
  changed no fit, and removing it leaves every fit bit-identical. A `data_list`
  still carrying either name gets a deprecation message. Time-varying growth
  goes through `build_growth(linkages = )`. `fit_mod()` drops retired parameter
  blocks from `inits`, so warm starts from older fits keep working.

## Bug fixes

* **`data_check()` requires the catchability columns each switch actually
  reads** — `Index_sd` under `Estimate_index_sd = "Estimated"`,
  `Catchability_prior_sd` under `Estimated-with-prior` and `AR1`,
  `Time_varying_q_sd` under a penalized `Time_varying_q`, and
  `Catchability_init` for every form that reads it. All are taken as `log()`
  when the parameter list is built, so a blank or non-positive entry gave a
  non-finite objective reported by TMB, naming neither the fleet nor the column.
  Each condition matches `build_map()`'s own gate, so settings that read no
  starting value are not asked for one.

* **`data_check()` requires the catchability settings wherever an index is
  fitted.** `Catchability`, `Catchability_init` and `Estimate_index_sd` have no
  schema default, so on a fishery they arrived `NA` and the fit died inside
  `build_map()`. The error now names the fleet and the column. Those
  configurations did not run before either.

  Two limits worth knowing. An index needing a selectivity **different** from
  the fishery's still has to be its own fleet: `sel_at_age` is indexed by fleet.
  The same is true of `flt_units`, so a CPUE series in numbers cannot sit on a
  fleet whose catch is in weight.

* **A missing `Ceq` is reported instead of crashing `data_check()`.**

* **`data_check()` reports fleets that share a `Selectivity_index` but disagree
  on the columns that shape the curve**, which are read per fleet and so do not
  give one shared curve despite the shared parameter block. A blank counts as a
  value: an empty `Sel_norm_bin` means "do not normalize".

* **`data_check()` reports a shared `Catchability_index` that does not share a
  catchability.** `Analytical` and `AnalyticalArith` solve q from each fleet's
  own index observations, bypassing the group parameter, so it is reported on
  the form rather than only on a disagreement.

* **A catchability linkage that names only part of a shared `Catchability_index`
  group is reported.** The offset is added per fleet and not reconciled across
  the group, so the fleets end up with different catchabilities while
  `fleet_control` still says they share one.

* **A catchability linkage on a fleet whose `Catchability` is `"Environmental"`
  or `"AR1"` is now refused.** Those forms overwrite `index_q`, so the linkage
  offset never reached the likelihood and its parameters were free with an
  identically-zero gradient. Express the whole relationship as the linkage and
  set `Catchability = "Estimated"`.

* **The diagnostic refits no longer repeat `data_check()`'s warnings**, which
  describe the `data_list` and would otherwise appear per peel, jitter or MSE
  iteration. `fit_mod()` gains `quiet_data_check` (default `FALSE`), set by
  every refit. Errors still stop the fit.

* **A fishery's index now appears in the index diagnostics.** `plot_index()`
  selected `Fleet_type == "Survey"`; both it and `plot_indexresidual()` now
  select the fleets carrying `index_data`. `plot_indexresidual()` previously
  applied no fleet filter at all, so it also drew residuals for `Off` fleets —
  `GOA2018SS` fleet 7 has nine such rows.

* **Diagnostics for a natural-scale survey index no longer use log-scale
  formulas.** `Lognormal` carries a log-sd; `MVN`, `MVNORM`, `Normal` and
  `TruncatedNormal` carry an ABSOLUTE sd. Applying the lognormal form to those
  does not error, it returns nonsense, because `sigma^2 / 2` is a number the
  size of the index squared.

  - `residuals(type = "pearson", source = "index")` returned the same large
    constant for every row of a natural-scale fleet (about `+75` for an absolute
    sd of 150, whatever the fit). It now standardizes as `(obs - hat) / sigma`.
  - `plot_index()` drew its observation interval with `qlnorm()`, giving bands
    like `[0, 1e130]`. Natural-scale fleets now get `obs +/- 1.96 * sigma`,
    clamped at zero.
  - `plot_indexresidual()` now uses the plain difference where the fleet is
    fitted on the natural scale, labelled `"Index residual"`.

* **`plot_indexresidual()` plotted the negative of a residual.** It drew
  `predicted - observed` while `residuals()` returns `observed - predicted`.
  Both are now `observed - predicted`, the usual assessment convention (WHAM,
  SS3). **Index residual plots made before 5.9.0 are mirrored about zero**; a
  figure carried forward needs redrawing and any written interpretation of its
  sign needs rereading. `Lognormal` fleets are unchanged.

* **`Estimate_index_sd = "Analytical"` no longer passes silently on a
  natural-scale index family.** The analytical sd (Ludwig and Walters 1994) is
  accumulated from squared *log* residuals. `Normal` and `TruncatedNormal` read
  it as an absolute value, so the likelihood was evaluated on the wrong scale —
  `data_check()` now refuses that combination. `MVN` and `MVNORM` score through
  `index_cov` and never read the scalar sd, so their fits are unaffected, but
  the diagnostics divide by it, so those warn rather than refuse. `Lognormal` is
  untouched.

* **`TruncatedNormal` OSA residuals were the untruncated family's.** The
  truncation constant is a function of the prediction, not the observation, so
  it drops out of any method reading the curvature of the density in the
  observation — including `oneStepPredict()`'s default. `osa_residuals()` now
  residualizes this family in its own call with `method = "oneStepGeneric"` and
  a range starting at zero. The correction is the size of the truncated mass: on
  a fleet predicting 100 with an absolute sd of 150, an observation exactly at
  the prediction moves from 0 to **-0.44**. `"Normal"` is untouched; both are
  pinned against their closed forms in
  `tests/testthat/test-likelihood-index-truncated-normal.R`.

  Three consequences, all in `?osa_residuals`: exact integration does not always
  converge on a random-effects model, in which case the fleet falls back to
  TMB's spline approximation with a warning; `sd` is `NA` for the group and
  `predicted` is the truncated conditional mean; and adding a truncated fleet
  shifts the other fleets' residuals on a random-effects model — a different
  conditioning sequence, not a wrong value. `attr(x, "method")` records what ran.

* **`compare_sim()` warns when it is handed `process`-drawn replicates**, since
  every statistic it reports is a deviation from `operating_mod`, which is the
  truth only when the replicates redrew the observations alone.

* **`osa_residuals()` reported `observed = NA` for any group residualized with
  `oneStepGeneric`.** Affected the discrete composition path as well as the new
  truncated one.

* **`M1_re = 3` and `6` estimated far fewer deviations than intended.** The
  age-by-year map index was built from a vector of length `nyrs_hind` and
  assigned into an `nages x nyrs_hind` block, which R recycled silently — on
  `GOA2018SS`, **42 free deviations where 882 were meant**. Fits in these two
  modes change; the other `M1_re` modes recycle deliberately and are unaffected.

  These deviations are integrated out rather than estimated as fixed effects, so
  the model class is sound, but a free age-by-year mortality field with a single
  standard deviation and no prior on M is strongly confounded with selectivity
  and recruitment in a single-species assessment, and the latent dimension is now
  over twenty times larger. Read `fit$convergence` and the estimability table
  before trusting such a fit, and consider a prior on M.

* **The 2D-AR1 correlations acted on the wrong dimensions.** `SEPARABLE(f, g)`
  applies `f` to the outermost array dimension and `g` to the fastest-running
  one, and both 2D-AR1 densities pass their fields as `(bin, year)` or
  `(age, year)`, so the two correlations were exchanged. Affected `M1_re = 6`
  and `Selectivity = "2DAR1"`. The 3D-AR1 form was already bin/year/cohort.

  **Most fits do not move.** Both correlations start at 0 and `TMBphase()` holds
  them there, so SSB, F and reference points are unchanged and only the two
  reported estimates exchange names. A fit moves only where the two differ as
  *starting* values — for selectivity, unequal `Sel_curve_pen1` /
  `Sel_curve_pen2`.

  **Warm starts across this release are the case to watch.** `inits` from an
  older fit carry `M1_rho` and `sel_curve_pen` under the previous convention, so
  a refit starts from a mirrored point. That covers `retrospective()`,
  `jitter()`, `run_mse()` and any two-stage bootstrap. Refit from
  `inits = NULL`, or transpose the two slots, if the starting values differ.

* **A composition or CAAL row with no sample size came back holding the
  predicted proportions.** `comp_hat` rows are normalized to sum to one, so a
  row with `Sample_size = 0` was returned as noise-free proportions —
  indistinguishable from a real full-weight composition. `run_mse()` reaches
  this state directly, so the next assessment would have been handed a
  perfectly-observed composition for a year that was never sampled. Such a row
  now comes back empty.

* **`obj$simulate()` returned the drawn compositions under the observed data's
  own name.** `comp_obs` and `caal_obs` are now drawn into copies, as the catch,
  index and diet draws already were. `fit$quantities` was never affected.

* **The natural-scale survey draw did not update `obsvec`**, leaving the OSA
  path scoring different data from the draw.

* **Recruitment is no longer redrawn when two densities score it.** Under
  `srr_fun = 0` with `srr_pred_fun > 0` — the AMAK/Ianelli configuration —
  `rec_dev` is scored by both `JNLL_REC_DEV` and `JNLL_SRR_PENALTY`. Two terms
  on one latent do not compose into a distribution to draw from, so
  `self_test(process = "recruitment")` would have measured recovery against a
  process the model does not assume. The deviations are left at their fitted
  values and `sim_mod()` says why; a random linkage on a recruitment parameter
  is a separate latent and is still drawn.

* **The correlated survey draw could write an unchecked non-positive index.**
  The rejection loop tested its budget before testing the draw, so exhausting
  the budget wrote out whatever had last been drawn.

* **`data_check()` rejects a `diet_data` that is not sorted by `stomach_id`.**
  The diet likelihood walks `diet_ctl` with a single forward cursor, so the ids
  have to run 0, 1, 2, ... in order and with no gaps. Out of that order whole
  stomachs drop out of the likelihood with no warning and a lower objective:
  re-sorting a cleaned `BS2017MS` diet table by predator age leaves **3 of its
  45 stomachs** in the fit, and reversing the rows leaves 1. `clean_data()`
  sorts by `stomach_id`; this catches a hand-built or re-sorted table.

* **`data_check()` allows a stomach summing to 1 within floating-point
  tolerance** (`> 1 + 1e-12` rather than `> 1`). Real excesses are still
  rejected.

* **`sim_mod()` sizes the truncation warning correctly, per family.** For the
  independent `"Normal"` family the warning reports the gap as
  `P(draw <= 0) = Phi(-mu/sd)` on the worst row. The previous test counted
  rejections, but each row is drawn once per call, so the ratio could only be
  0, 1/2, 2/3, ... — an indicator rather than a rate, against an unreachable 2%
  threshold. A covariance fleet is judged on the JOINT rate, since it rejects
  the whole vector whenever any row is non-positive: an 8-row fleet whose worst
  margin is 31% is rejected about **92%** of the time.

## Documentation

* The developer guide covers the `SIMULATE` layer: what a new likelihood term
  owes its simulator, the `*_sim` reporting convention, and the two registries a
  new or retired parameter block has to be added to.

* `vignette("model-diagnostics")` documents `process = `, the natural-scale
  index families, and what a redrawn M can and cannot be read as. On `BS2017SS`
  with a year-varying M random effect, refits recover the simulated deviations
  (correlation around 0.5, positive in every replicate) but recover their
  standard deviation poorly — roughly 70% low, collapsing to zero in a third of
  replicates. That is an ordinary variance component estimated from little
  information, not a fault in the simulation;
  `tools/verify/verify-sim-recovery-M.R` reproduces the numbers.

* `vignette("model-options-and-functionality")` records which `Sel_curve_pen`
  slot is which correlation under the AR1 selectivity forms, and the
  estimability caveat on `M1_re = 3` / `6`. `vignette("growth-estimation")`
  gains a time-varying growth section. `vignette("hcrs-and-mses")` notes that
  seeded MSE results are not comparable across 5.8.x to 5.9.0.

# Rceattle 5.8.1

Documentation-only release. Four changes that landed earlier in the 5.x line
reached `main` without a NEWS entry; they are recorded here so the release
notes describe what actually shipped. No code changes.

## Changes that move results, recorded late

* **The Dirichlet-multinomial composition likelihood uses a different effective
  sample size.** The concentration is now built from the observed-count total,
  `sum(comp_obs)` = `N * (1 + offset * nbins)`, rather than the raw input `N`.
  That is the same `N` `ddirmultinom()` reads back from `obs.sum()`, so the
  alpha and the density are now consistent, and it matches the CAAL
  Dirichlet-multinomial and the AFSC linear-DM parameterization. **Any fleet
  with `Comp_distribution = "DirichletMultinomial"` will move**; multinomial
  fleets are unaffected. Refit before carrying a DM-weighted model forward.

* **`build_growth(sd_plus_group = )` defaults to `NA` ("inherit") rather than
  `"WHAM"`.** `NA` takes `data_list$growth_sd_style` when the data list carries
  one, falling back to `"WHAM"` otherwise, so a fresh fit is unchanged. The
  fix is for refits: a diagnostic that rebuilt growth through
  `build_growth(fun = )` previously reset an SS3-style plus-group SD back to
  WHAM without saying so. Retrospectives, jitters and self-tests of an SS3
  growth model change accordingly.

## Deprecations

* **`Catchability = "Environmental"` is deprecated.** It still fits with its
  existing numerics and results are unchanged, but it names covariates by
  position in `Time_varying_q` rather than by formula, and cannot carry priors,
  bounds or an estimation phase. Use a catchability linkage instead:
  `build_catchability(linkages = list(q = linkage_spec(~ temp, by = ~ fleet)))`.
  `switch_check()` now names the affected fleets.

## Known issues

* **The reported case-0 (full multinomial) composition NLL is not comparable
  with 4.9.x.** Composition fitting routes through `dmultinom_osa()`, which
  renormalizes the predicted proportions; the previous `dmultinom()` did not.
  The reported per-row NLL therefore shifts by an additive constant,
  `Neff * log(1 + n_comp * offset)`. The gradient, the MLE and every derived
  quantity are unchanged -- only the printed likelihood value moves. Compare
  likelihoods within a version, not across this boundary.

# Rceattle 5.8.0

## New features

* **`fleet_control$Sel_norm_scope` controls whether selectivity normalization
  pools its reference across sexes.** Normalization makes two independent
  decisions -- *where* the reference is taken and *whose* scale it sets -- and
  until now the second was a silent passenger on the first: a named
  `Sel_norm_bin` always used a per-sex reference, and max-normalization always
  pooled across sexes, so two of the four combinations were unreachable. The new
  column is orthogonal to `Sel_norm_bin`, so the two combine rather than
  multiplying into more columns:

  * `"AcrossSexes"` (default) -- one reference pooled over the sexes, so the
    less-selected sex stays below 1 and **relative sex-specific selectivity is
    retained**. This is what a dimorphic stock usually wants.
  * `"WithinSex"` -- each sex divided by its own reference, so both reach 1 and
    only the *shape* differs by sex.

  The newly reachable combination worth knowing about is max-normalization with
  `"WithinSex"`: each sex is scaled at its own plateau, wherever that falls. For
  a stock whose sexes plateau at different ages that is more robust than naming
  a bin, because a named bin stops being the plateau once selectivity is
  time-varying while the maximum tracks it.

  The scope has no effect on a one-sex species, or where `Sel_norm_bin` is `NA`
  and nothing is normalized. **The one behaviour change** is a *two-sex* fleet
  normalizing at a *named bin*: it previously used a per-sex reference and now
  pools, and `switch_check()` emits a message naming the fix (set
  `Sel_norm_scope = "WithinSex"`) when it sees that configuration. Max-normalized
  fleets -- including the operational arrowtooth configuration -- already pooled
  and are unaffected. No bundled dataset is in the affected case.

* **Biomass, SSB and recruitment confidence intervals are now taken on the log
  scale.** The model ADREPORTs `log_biomass` and `log_R` alongside the existing
  `log_ssb`, and `plot_timeseries()` builds `exp(log(x) +/- 1.92 * sd_log)` from
  the delta-method SD of `log(x)` whenever a model reports one. `sdreport()`
  linearizes once about the MLE, so its SD is exact only for a linear function
  of the parameters; these three series are built multiplicatively
  (`R = R0 * exp(rec_dev)`, and n-at-age is a product of survivals), so `log(x)`
  is close to linear in the estimated parameters where `x` itself is
  exponential. The interval is therefore both better approximated and
  right-skewed, and cannot cross zero the way the symmetric natural-scale
  interval did for weak year classes and depleted stocks. The log-scale SD is
  also the CV, the form the ABC / OFL buffer calculations want. Natural-scale
  `biomass`, `ssb` and `R` are still ADREPORTed, so existing callers of
  `sdrep$value` are unaffected. Models fit before this release have no `log_*`
  rows, but still get the log-scale interval -- see below.

* **The `plot_*()` timeseries wrappers take a `ylab` argument.** It is appended
  to the argument list, not inserted, so positional calls are unaffected. Left
  `NULL` (the default) the axis label is derived from the series and the model's
  `minage`, so a model whose recruitment is at age 3 no longer gets an axis
  labelled "Age-1 recruits". Only recruitment names an age now -- it is an age
  class, so the age is information. The biomass and SSB axes drop their `Age-1+`
  prefix and read `Biomass (million mt)` and `SSB (million mt)`; the prefix was
  noise on an aggregate, and wrong on SSB, which is mature females rather than
  the minage+ stock.

* **`plot_exploitable_biomass()` and the two depletion plotters accept
  `add_ci = TRUE`, on the log scale.** `exploitable_biomass`, `ssb_depletion`
  and `biomass_depletion` are now ADREPORTed, so an interval can be computed at
  all. They are ADREPORTed on the **natural scale only**, deliberately:
  `exploitable_biomass` sums over fisheries alone, so it is exactly 0 for any
  model without projection F (`Proj_F_proportion = 0`, which is every bundled
  dataset -- all 216 species-years of `BS2017SS`), and the depletions divide by
  `B0` / `SB0`; `log()` of either would put `-Inf` on the AD tape and turn the
  entire `sdreport`, not merely that row, into `NaN`. They still get a
  log-scale interval, because the delta method defines `sd(log x) = sd(x) / x`
  exactly -- the plotters recover the log-scale SD from the natural-scale one
  wherever the series is positive. Where it is 0 the quotient is undefined and
  the (degenerate) symmetric interval stands. sdreport grows about 20% on
  `BS2017SS` for the three extra series.

* **Models fit before this release get log-scale intervals too.** The same
  `sd(log x) = sd(x) / x` recovery applies to any strictly positive series with
  a natural-scale SD, so an existing fit with no `log_*` rows plots a
  right-skewed interval without being refit. Where the model does report
  `log_biomass` / `log_ssb` / `log_R`, those are used directly; the two agree to
  machine precision (2e-16 on `BS2017SS`).

## Bug fixes

* **`as.data.frame(fit)` and the plotters now report the same interval.** The
  extractor still built a symmetric natural-scale interval for `biomass`, `ssb`
  and `R` while `plot_biomass()` drew the log-scale one, and it returned `NA`
  standard errors for `exploitable_biomass` and the two depletions -- the three
  series this release just ADREPORTed. Both paths now go through one helper.
  The plotters also switch from the rounded `1.92` to `qnorm(0.975)`, so the
  band labelled 95% is a 95% band: **ribbons are about 2% wider than in 5.7.0**,
  and any interval read off a figure or table should be regenerated. The point
  estimates are unchanged.

* **`add_ci = TRUE` on a fit made with `getsd = FALSE` no longer errors.** It
  warned that it would draw no interval and then failed on the next line
  indexing an absent `sdreport`. It now warns and draws the series without a
  ribbon, as the message says. The warning is also gated on whether the series
  is ADREPORTed at all, so `plot_f(add_ci = TRUE)` no longer warns on every fit
  about standard errors `F_spp` can never have.

* **`rearrange_data()` and `data_check()` work again on a `fleet_control` that
  has not been through `switch_check()`.** `Sel_norm_scope` was read with no
  missing-column guard, so a data list read straight from a workbook failed with
  a cryptic "column not found". It and `Index_distribution` are now filled from
  the column schema. A blank `Sel_norm_scope` cell is filled too -- an `NA`
  reached TMB as `NA_integer_`, which the C++ reads as `WithinSex`, the opposite
  of the default.

* **The `Sel_norm_scope` behaviour-change notice fires only where it applies.**
  It was skippable (leaving one cell blank flipped that fleet silently) and it
  fired on Hake and LogisticPM fleets, which use `Sel_norm_bin` as a penalty
  age-range rather than a normalization reference, and on `Off` fleets. The gate
  now matches the one in `selectivity.hpp`.

* **A one-species plot drops the redundant species strip, and absolute y-scales
  start at zero** so panel heights are comparable. Depletion and other relative
  series are unaffected.

* **Diet defaults are announced only when the model reads them**, so a
  single-species fit no longer reports settings it never uses.

* **`model_config()` names the mistake when handed a `data_list`.** It takes
  model settings, not data, but almost every other entry point takes a
  `data_list` first -- so `model_config(my_data)` bound the list to `msmMode` and
  reported only "`msmMode` must be a single value". It now points at
  `build_data(base = , model_config = )`. The check itself was correct:
  `msmMode` is a model-wide scalar (unlike the per-species `suitMode`), and a
  single-species `msmMode > 0` cannibalism configuration was never blocked.

* **A confidence interval that cannot be drawn now warns instead of vanishing.**
  Requesting `add_ci = TRUE` for a series with no standard errors indexed out of
  range, filled the interval with `NA` and rendered an invisible ribbon -- no
  error, no warning, just a missing interval. That is how the three plotters
  above hid their missing `ADREPORT`. `plot_timeseries()` now says which series
  and which model lacks standard errors.

* **Recruitment was plotted 1000x too high, and the stock-recruit panel 1000x
  too low.** The model carries numbers-at-age in thousands and weight-at-age in
  kg, so biomass comes out in mt and recruitment in thousands of fish.
  `plot_recruitment()` never applied the matching `/1e3`, plotting thousands of
  recruits under an axis reading "Age-1 recruits (million)", while
  `plot_stock_recruit()` divided by `1e6` under an axis reading "Recruitment
  (millions)" -- so the two panels disagreed with each other by a factor of a
  million. Both now plot millions of recruits. `plot_exploitable_biomass()` was
  labelled "million mt" with no rescaling applied at all and is now divided by
  `1e6` like the other biomass series. Divisors and axis units are now held in
  one table (`.RCE_TS_RESCALE` / `.rce_ts_ylab()`) that `plot_stock_recruit()`
  reads too, and are asserted against each other in the test suite. **Any figure
  or number read off these three plotters needs regenerating.** Model results
  are unchanged -- this is display-only.

## Documentation

* **How sex-specific selectivity works is now documented.** Every selectivity
  form is sex-specific by default in a two-sex model. `Sel_norm_bin` says only
  *where* the normalization reference is taken; whether the *relative* scale
  between sexes survives is `Sel_norm_scope` (above). Relative sex selectivity
  is only informed where the composition is joint (`comp_data$Sex = 3`);
  sex-specific rows each sum to 1 and carry no sex-ratio information, and
  `Selectivity = "Hake"` always normalizes within sex regardless of either
  column. New "Sex structure and relative selectivity" section in
  `vignette("model-options-and-functionality")`, with the details on the
  `Sel_norm_bin` and `Sel_norm_scope` field dictionary entries.

* **The mt / thousands-of-fish convention is now stated consistently.** The
  `biomass` / `ssb` / `catch_hat` declarations in `ceattle.cpp`, the
  `index_data` / `catch_data` / `Observation_units` field dictionary entries and
  the Stock Synthesis catch-conversion table all said kg. Section 12 of
  `vignette("model-options-and-functionality")` now covers the display units,
  `ylab` and the log-scale intervals, and names the discrete palette correctly
  (Okabe-Ito, not viridis).

* `Sel_norm_bin` and `Sel_norm_bin_upper` now say that they are bin indices on
  the fleet's own `Selectivity_dimension` -- an absolute age for an age-based
  fleet (`6` means age 6, not the sixth bin), a 1-based length-bin ordinal for a
  length-based one.

* `fleet_control$Sex` is documented as inert. It is read nowhere in the model;
  sex is set per observation on `comp_data$Sex`. Retained only so older
  workbooks still read.

* `vignette("model-parameterizations")` no longer claims the model is fit to
  sex-ratio data. There is no sex-ratio likelihood component -- the text
  described the ADMB implementation. Sex-ratio information enters through joint
  composition data.

* `plot_timeseries()` documents the unit convention it assumes: numbers-at-age
  in thousands and weight-at-age in kg, hence biomass in mt and recruitment in
  thousands of fish.


# Rceattle 5.7.0

## Breaking changes

* **The remaining bundled datasets use the canonical `fleet_control` column
  names.** The 4.10.0 renames (`Q_prior` -> `Catchability_init`,
  `Comp_loglike` -> `Comp_distribution`, ...) had reached `BS2017SS`,
  `BS2017MS` and `GOA2018SS` only; `Atka2022`, `GOA2018SS`, `GOAatf`,
  `GOAatf2023`, `GOAcod`, `GOApollock`, `GeorgesBank3spp`,
  `NorthernRockfish2022` and `whamGrowthData` still shipped the legacy
  spellings, so loading one printed a dozen upgrade messages. A script reading a
  legacy name off one of these objects directly -- e.g.
  `GOApollock$fleet_control$Comp_loglike` -- now gets `NULL`. Use the canonical
  name; a workbook on disk still reads fine, since the aliases are upgraded on
  read.

* **`mse_summary()`'s metric columns now have syntactic names.** They were
  display strings held together by `check.names = FALSE`, so reading one meant
  `summ$species[["OM: Terminal SSB Depletion (Dynamic)"]]`. They are now
  ordinary names -- `avg_catch`, `catch_iav`, `p_closed`, `om_terminal_ssb`,
  `om_avg_depletion` and so on -- with an `om_` prefix for the operating model's
  truth and `em_` for the estimation model's perception. The four
  misclassification metrics become `p_overfishing_false_pos` / `_false_neg` and
  `p_overfished_false_pos` / `_false_neg`, and the three collapse metrics become
  `om_sims_collapsed`, `om_no_f_sims_collapsed` and `om_sims_collapsed_from_f`
  -- named `sims` because they are **counts of simulations**, not probabilities,
  which the old names did not say. The long display strings are kept on a
  `"labels"` attribute of each frame, so plots and SAFE tables can still print
  them: `attr(summ$species, "labels")`. `Species`, `Fleet_code` and `Fleet_name`
  are unchanged. Values are unchanged. Done in the same release as the
  per-entity reshape so `mse_summary()` breaks once rather than twice; no script
  in `Rceattle-models` or `GOA-ATF-ESP` reads a metric by name.

* **`data_requirements()` argument `selectivity` is now `Selectivity`, and
  `index_distribution` is now `Index_distribution`.** They stand in for the
  `fleet_control` columns of those names, and the rest of the package's
  arguments already mirror their source exactly (`estDynamics`, `Ceq`,
  `growth_model`, and `model_config()`'s mirror of `fit_mod()`). No deprecation
  path: `data_requirements()` is new in this release line and has never been on
  a released version.

## Bug fixes

* **`self_test()` no longer fails on every replicate for a model with a
  natural-scale survey index.** `sim_mod()` drew `Normal` and `MVN` index
  observations from the untruncated natural-scale normal, so a fleet with an
  absolute sd near its index -- an AVO-type index with sd 0.25-0.80 against
  observations from 0.70 -- drew non-positive most replicates. `data_check()`
  rejects those, so the refit failed and the replicate was counted *not
  converged*, reading as a convergence problem rather than a simulation one.
  Non-positive draws are now redrawn, correlated fleets jointly. On EBS pollock
  `self_test()` goes from 0 of 50 replicates to all of them, recovering SSB with
  a median bias of -0.4%.

  Redrawing samples from a normal **truncated at zero**, not the normal the
  likelihood maximises, so a fleet that needs it is self-testing against a
  different data-generating process. At index 0.70 with sd 0.80 a fifth of draws
  are rejected and the truncated mean is 39% high. Two warnings now fire: the
  existing one when a fleet still cannot draw positive, and a new one when any
  row is redrawn more than 2% of the time. Either means the absolute sd is too
  large relative to the index for a natural-scale normal.

  Lognormal fleets are unaffected, including the seeded all-lognormal fast path,
  and seeded results are bit-identical for any `Normal`/`MVN` fleet that never
  rejects. A fleet that does reject consumes extra random numbers, shifting that
  replicate's comps and catch and every later replicate in a seeded loop.

* **`data_check()` no longer errors on a `data_list` that carries no `HCR`.**
  `HCR` is a `fit_mod()` argument, not a data field, so a list read straight
  from a workbook has none. `NULL %in% ...` is `logical(0)`; `TRUE && logical(0)`
  is `NA` and `x & logical(0)` is `logical(0)`, so `data_check()` died with
  `missing value where TRUE/FALSE needed` and `validate_switches()` with
  `argument is of length zero`. `fit_mod()` sets `HCR` before checking and was
  unaffected.

* **`run_mse()` no longer errors when a survey fleet has `NA`
  `Proj_F_proportion`.** The check that some fleet takes projected F summed the
  column without `na.rm`, so the `NA` that fleets taking no catch legitimately
  carry made the sum `NA` and the `if` failed with `missing value where
  TRUE/FALSE needed` before the MSE started.

* **`plot_indexresidual()` now honours `incl_proj`.** The argument was accepted
  and then ignored, so projection years were always plotted -- the opposite of
  its documented default, and of `plot_index()` / `plot_catch()`. A projection
  row's "observation" is whatever placeholder its workbook carries (the
  roll-forward scripts in `Rceattle-models` write 99999), which drew as an
  enormous spurious residual. **This removes rows from the default output** for
  any model whose `index_data` runs past `endyr`; pass `incl_proj = TRUE` to
  keep them. No bundled dataset has such rows -- `clean_data()` extends
  `catch_data` to the projection horizon but not `index_data`, and `run_mse()`
  marks its projection index rows with a negative `Year`, which was already
  filtered -- so in practice this bites only workbooks that carry their own.
  Applies to the default `residual_type = "pearson"`; the `"osa"` path returns
  `osa_residuals()` directly and takes none of the plot arguments.

* **`plot_catch(incl_proj = TRUE)` no longer errors on projection rows.**
  `clean_data()` gives projection years an `NA` catch, and the lognormal
  error-bar mask in `.fleet_fit_df()` passed that `NA` to a subscript. The mask
  now excludes `NA` observations, and those rows draw the fitted line with no
  error bar, as non-positive observations already did. `plot_index()` shares the
  code path and is fixed with it, though only a workbook that supplies `NA`
  index observations can reach it.

* **`data_check()` now warns when CAAL data sit on a fleet whose
  `Selectivity_dimension` is not `"Length"`.** Conditional age-at-length is the
  age composition within a length bin, so the model predicts it from
  selectivity-at-length convolved with the growth matrix. An age-dimensioned
  fleet has no selectivity-at-length, so the prediction is zero -- and the
  likelihood still runs, scoring the observations against the flat composition
  the offset leaves behind. Those data then add a term to the objective that no
  parameter can move: on the test fixture the CAAL component costs 64.6 with
  uniform observations and 325.8 with realistic ones, and the whole objective
  shifts by exactly that amount. `Selectivity_dimension` defaults to `"Age"`, so
  this was the default outcome for anyone adding CAAL data. It warns rather than
  errors because such models currently fit; the sibling case (empirical growth,
  which zeroes the prediction for the same reason) has always errored, and this
  should tighten to match once configurations have been checked.

* **`data_check()` now rejects a `diet_data` that is not sorted by
  `stomach_id`.** The TMB diet likelihood walks `diet_ctl` with a single forward
  cursor, taking stomach *i*'s prey as the rows where `stomach_id == i`, so the
  ids have to run 0, 1, 2, ... in order and with no gaps. Out of that order the
  cursor runs past whole stomachs, which then drop out of the likelihood with no
  warning and a lower objective: re-sorting a cleaned `BS2017MS` diet table by
  predator age leaves 3 of its 45 stomachs in the fit, and reversing the rows
  leaves 1. `clean_data()` sorts by `stomach_id`, so anything that came through
  it already satisfies this; the check catches a hand-built or re-sorted table.

* **`mse_summary()` now reports the catch metrics for a species that has no
  fishery.** A multispecies model can carry a predator purely for its
  consumption -- arrowtooth and sablefish in the hake model do exactly this --
  and such a species has no `catch_data` rows at all. Every per-species catch
  metric then operated on a zero-length vector: `Average Catch` warned
  "argument is not numeric or logical: returning NA" and returned `NA`, while
  `Catch IAV` and `P(Closed)` divided by a length of zero and returned `NaN`
  silently. They now report `Average Catch = 0` -- an unfished species really
  does have zero catch -- and `NA` for the two ratios, which are not defined
  when there is no fishery to vary or to close. A species that *does* have a
  fishery but no landings in the projection years is reported as `NA` with a
  warning naming the species and the years, because that is a data problem
  rather than a modelling choice. Metrics for fished species are unchanged
  (verified bit-identical on the 3-species hake MSE), as are `$fleet` and
  `$total`.

## Documentation

* **`?osa_residuals` now explains the negative composition `predicted` values it
  warns about**, in a new section: why a numerical conditional mean overshoots
  below zero on a near-empty bin, that it follows the *bin's* count rather than
  the composition's sample size (so a rare age in a well-sampled year does it),
  and that those residuals are biased positive. The warning itself is now two
  sentences pointing at that section rather than carrying the whole explanation,
  because it fires once per call on real assessments -- 269 rows on GOA-ATF-ESP
  2025, 98 on GOA pollock 2025.

* Corrected the `comp_osa.hpp` note claiming the `pUsed` squeeze in
  `dmultinom_osa()` matches WHAM. WHAM does not squeeze `pUsed`; it squeezes
  `1 - pUsed` at point of use, which Rceattle also does. The extra squeeze is
  redundant and leaves an O(A * eps) difference from WHAM on an A-bin
  composition. Comment only -- no value changes.

* Fixed two doubled words in user-facing text: "the The Pacific Fishery
  Management Council" in `?build_hcr`, and "the the diet" in the `diet_data`
  schema description, which is written verbatim into the bundled
  `meta_data_names.xlsx` (regenerated to match).

# Rceattle 5.6.1

## Bug fixes

* **`fit_control(verbose = 0)` now silences the estimability table printed when
  a fit does not converge.** `TMBhelper` emits it with `print()`, which
  `suppressMessages()` cannot catch, so it appeared whatever `verbose` said --
  one row per fixed parameter, repeated for every refit a `self_test()`,
  `jitter()` or `retrospective()` run makes. It and the "Model did not converge"
  banner now appear only under `verbose > 0`. The verdict is unchanged on
  `fit$identified` and in `fit$convergence`.

* **A covariance (`MVN`/`MVNORM`) survey `Sigma` that is symmetric but not
  positive definite is now caught by `data_check()`.** The covariance index
  likelihood factorizes `Sigma`, so symmetry is not enough. An indefinite or
  singular matrix used to clear every check -- presence, squareness, dimension,
  symmetry -- and then fail inside the TMB objective, in a message that named
  neither the fleet nor `Sigma`. The error now names the fleet and reports the
  smallest eigenvalue, so it is clear how far off the matrix is.

* **An `index_cov` entry for a fleet that is not using the covariance likelihood
  now warns.** It was ignored in both `data_check()` and the internal aligner, so
  supplying a covariance but leaving `Index_distribution` unset gave a silent
  lognormal fit with no indication the matrix had been dropped. Entries keyed by
  either `Fleet_name` or `Fleet_code` are recognised, so a correctly configured
  fleet does not warn.

# Rceattle 5.6.0

## New features

* **`self_test(debug = TRUE)` returns every simulation, not just the converged
  ones.** When a self-test comes back short, the dropped runs are the ones worth
  looking at, and each already carries its own `$convergence` diagnostics --
  they were simply discarded. Under `debug`, `Sim_i` is simulation `i` (so it
  pairs with the seed `seed + i` that produced it) and the verdicts are in
  `attr(sims, "converged")`. The default return is unchanged. A refit that
  errors outright is now caught and returned as its condition rather than
  aborting the whole call (and, under a cluster, every other replicate with
  it) -- previously the hardest failures were the ones you could not inspect.

* **`jitter(timeout = )` and `self_test(timeout = )` bound a single run.**
  `.fit_tmb()` optimizes with `eval.max = iter.max = 1e9`, so a replicate that
  wanders somewhere pathological has no bound and one stall blocks the whole
  run -- a failure no convergence check can reach, because the fit never
  returns. A run past the limit is stopped, counted as non-converged, and
  reported separately from the gradient drops, since the fix is different
  (raise the limit, versus look at the model). The limit is approximate:
  `setTimeLimit()` is checked when control returns to R, so it fires between the
  optimizer's function evaluations rather than inside one. Both functions now
  also trap errors per run, so one bad replicate cannot abort the rest.

* **`self_test(start = )` chooses the starting values.** `"initial"` (the
  default, and the previous behavior) starts each refit from
  `initial_params` -- what the original fit started from -- so the estimator
  covers the same ground on simulated data that it did on the real data.
  `"estimated"` warm starts from `estimated_params`, which is faster and far
  more likely to converge but asks the estimator to travel no distance: it tests
  that the optimum is stable under resampling, not that it is reachable.

* **`reweight_comps(fleets = )` accepts fleet names as well as codes,** matching
  `linkage_spec(fleet = )`. Passing a name previously did not error: it fell
  through as ineligible and surfaced as a warning about having no fitted
  composition data, followed by "No fleets to reweight". An unknown name is now
  reported against the model's own fleet list, and mixing ids and names in one
  vector is caught with the coercion hint.

## Bug fixes

* **`self_test()` can now phase its refits, and inherits the setting by
  default.**
  Every refit starts from the source model's *starting* values, so it covers the
  same ground the original fit did -- but `self_test()` had no way to phase it.
  `retrospective()` fixes `phase = TRUE` ("phasing, or the parameters dont wanna
  move") and `jitter()` exposes the argument; `self_test()` took
  `.refit_like()`'s `phase = FALSE` default with no argument to change it. A
  model that needed phasing to fit its real data therefore had to reach the same
  optimum in one unphased pass for every simulated one. On a GOA pollock
  configuration that fits at a maximum gradient of 2e-4, 5 of 6 simulations
  ended between 4e9 and 2e13; with phasing all 6 converged to ~1e-4. (That model
  carries an `sdrep`, so the self-test ran at `getsd = TRUE` and those five were
  dropped through the non-positive-definite path below -- under `getsd = FALSE`
  the old gate would have returned all six and counted them as converged.)
  `phase` reads the setting `fit_mod()` recorded on the source fit, so a model
  fitted under the package default is still refitted unphased; pass
  `phase = TRUE` for a model that needs phasing but was not fitted with it.

* **The re-fitting diagnostics' convergence gate now tests the gradient.**
  `retrospective()`, `jitter()`, `self_test()` and `profile()` each drop the runs
  that did not converge, by comparing `opt$Convergence_check` against the string
  `TMBhelper::fit_tmb()` uses for a non-invertible Hessian. That comparison could
  not work in either direction. With `getsd = TRUE`, `fit_tmb()` returns *before*
  assigning that string when the Hessian fails `chol()` -- so runs were dropped
  by the enclosing `is.null()` guard, incidentally rather than by the test, and
  by a criterion (positive-definiteness) that is not the one intended. With
  `getsd = FALSE` that assignment is unreachable, so nothing was ever dropped
  and a run ending at a maximum gradient of 1e13 counted as converged.
  (`Convergence_check` itself is always set, but to one of the two gradient
  verdicts, neither of which the comparison matched.)
  Both paths now go through one gate that reads the hindcast maximum gradient,
  using `convergence_diagnostics()`'s own FAIL threshold. Because the gate can
  now actually drop, all three report how many runs they dropped rather than
  returning a thinned list in silence -- for `jitter()` and `self_test()` a
  thinned list is a biased sample, since the failures are exactly the runs that
  would have shown the spread.

* **`retrospective()` survives a dropped peel.** It named its output
  `Year_(endyr - peels):endyr` -- always `peels + 1` labels -- and looped
  `1:(length(mod_list) - 1)`, both of which assume every peel converged. One
  dropped peel errored with `'names' attribute [3] must be the same length as
  the vector [2]`; all of them dropped made `1:0` count down and index past the
  end. The assumption held only because the old gate could not drop anything
  (above). Entries are now labelled from each model's own `endyr_peel`, so a
  dropped peel leaves a gap rather than shifting every later label onto the
  wrong model.

* **`retrospective()` now judges the peeled hindcast, not just the F refit.**
  Each peel fits twice: the peeled hindcast, then a refit with only the peeled
  years' `log_F` free. Only the second was gated, and its gradient says nothing
  about whether the hindcast converged -- so a peel whose hindcast diverged was
  kept as long as the handful of F parameters settled, and its parameters went
  into Mohn's rho. Both refits are now gated.

* **A non-positive-definite fit no longer returns a malformed `fit$opt`.**
  On that early return `TMBhelper::fit_tmb()` hands back
  `list(opt = , h = )` rather than the estimates, so `fit$opt$objective`,
  `$max_gradient` and `$Convergence_check` were all `NULL`. `fit_mod()` now
  unwraps it, keeps the Hessian as `$hessian`, and records the verdict
  `fit_tmb()` gives when `sdreport` reports `pdHess = FALSE`. `fit$sdrep` stays
  `NULL`, so the `sdreport_failed` diagnostic is unaffected. A no-op for every
  converged fit.

* **`data_requirements()` and `build_data()`'s pre-check read an attached
  `model_config()`.** They classified from the top-level switches only, so an
  object carrying `model_config(msmMode = 1)` reported the predation inputs as
  Ignored while `print()` on the same object showed it as multispecies.
  `fit_mod()` resolves the configuration onto the data list, and the requirement
  layer now sees the same thing. The overlay covers what `model_config()` carries:
  `estDynamics` and `Ceq` have no field there, so `NByageFixed` and the
  bioenergetics inputs still classify from the top-level slots.

* **The pre-check classifies an integer-coded `fleet_control`.** It runs before
  `switch_check()`, so a `Selectivity` still holding the code `0` never matched
  the string `"Fixed"` and `emp_sel` was not required. The predicates now accept
  the code, the canonical string, and the numeric-looking string a workbook read
  can leave behind.

* **A mistyped mode switch is caught before the fit, not after it.**
  `model_config()` validates `initMode`, `avgnMode` and `suitMode` against the
  switch maps, and `fit_mod()` checks them before doing any work. `fit_mod()`
  rebuilds a `model_config()` after the optimization to record the run
  configuration, so validating only there meant `avgnMode = 3` -- a switch the
  model never reads -- threw away a converged fit and its `sdreport`. That call
  is now guarded as well, so recording the configuration can never cost a caller
  a fit. The checks accept a per-species vector and an integer code, since
  `.refit_like()` feeds resolved values back.

* **`run_config()` errors on an unrecognized field** instead of silently dropping
  it, naming the valid set and pointing at `fit_control()`, where the mistake
  usually belongs -- `run_config(fit, getsd = FALSE)` was a no-op. Names are
  validated, not values.

* **The linkage grammar documented a term that does not exist.** A time block was
  described as `~ block(Year)` in the roxygen, the developer guide and the
  vignettes, but there is no `block()`: the fixed part of a linkage formula is
  passed to `model.matrix()`, so a model carrying that term failed at fit time
  with `could not find function "block"`. Documented as `~ cut(Year, ...)`, with
  a note that anything `model.matrix()` accepts will work.

## Behavior changes

* **`summary()` on a data object no longer prints the `model_config` block.**
  `print()` still does, and both take `config =` to override, so
  `summary(dat, config = TRUE)` and `print(dat, config = FALSE)` give the other
  behavior. An object with no configuration prints identically either way.

* **`read_data()`, `clean_data()` and `combine_data()` return an
  `"Rceattle_data"` object.** A data list keeps its spec-tree printer through a
  `write_data()` / `read_data()` round trip instead of only when it came from
  `build_data()`. The visible consequence is that these now print an indented
  specification tree at the prompt rather than dumping the ~40-element list. The
  tag is inert: every consumer treats the object as a plain list, and it is
  stripped before anything reaches `MakeADFun()`, so no fit changes.

* **`data_requirements(index_loglike = )` is now `index_distribution = `,**
  matching the canonical `fleet_control` column. The function has not been
  released, so there is no deprecation path.

## Documentation

* Runnable examples for `linkage_spec()`, `write_template()`, `load_config()`
  and `run_config()`, none of which had one. `linkage_spec()` is the entry point
  to the whole linkage grammar; its examples show the covariate, per-year,
  penalized random walk and intercept-prior forms.

* `reweight_comps()` is in the reference index (its absence aborted the pkgdown
  build), is referenced from the composition-diagnostics workflow in
  `vignette("model-diagnostics")`, and cites Francis (2011) as the alternative.
  `examples/McAllister-Ianelli-reweighting.R` now uses it rather than teaching a
  hand-rolled tuning loop.

* Reference documentation describes current behavior rather than narrating past
  refactors, and the release notes for 5.3.0 and 5.4.0 are condensed to what
  changed, who is affected, and what to do about it.

* **A failed parallel OSA loop no longer crashes the R session.**
  `osa_residuals()` runs the one-step-ahead loop in parallel by default, via
  forking, and retried serially if that failed -- but reusing the object the
  failed attempt had used, which ended the session with "An irrecoverable
  exception occurred" rather than recovering. The retry now builds a fresh
  object first. Found on GOA pollock 2025, where `osa_residuals(mod)` was
  unrunnable; it now returns its 1100 residuals, the same as
  `parallel = FALSE` gave.

  Forked workers abort on that model at any width above one core, and this is
  usually absorbed: a direct one-step-ahead call over the same observations
  returned all of them at 1, 2, 4, 8 and 11 cores. Occasionally the absorption
  fails, and that case did not reproduce on demand -- it appears to depend on
  machine load rather than on the model -- so the fix makes the retry safe
  rather than preventing the underlying failure. The aborting worker still
  prints its own message, which comes from C and cannot be suppressed; it does
  not mean the call has failed. `parallel = FALSE` avoids the attempt entirely.

* **A negative expected composition count now warns** (issue #108). `predicted`
  is an expected bin count and cannot be negative, but it goes slightly negative
  where a bin holds almost no fish: compositions enter as counts
  (`(proportion + comp_offset) * N`), and the one-step-ahead conditional mean is
  a numerical step from the observation that overshoots below zero when the
  count is near it. This follows the count in the *bin*, not the sample size of
  the composition -- on EBS pollock the affected rows have a median observed
  count of 0.05 against 4.9 elsewhere, while their sample sizes span the same
  range as everything else, so a rare age in a well-sampled year does it too.
  The affected years are named, and the warning notes that the residual on those
  rows standardises the observation against the same mean, so it is biased
  positive -- every affected row was positive in both Rceattle and WHAM. Nothing
  is dropped or clamped: the negative value is the signal that the bin is too
  sparse, and which rows to set aside is the analyst's call. This is not specific
  to Rceattle; WHAM produces negative expected counts on its own example, with
  the same default one-step-ahead method and composition counts on the same
  scale.

* **OSA residuals on an accumulated composition are labelled by the age they
  represent** (reported as issue #108). Tail accumulation folds the tails into a
  boundary bin, so a fleet with `Comp_accum_young = 3` on a 10-age model is fit
  on 8 bins -- but those were numbered 1 to 8, which labels the boundary residual
  "age 1" when it stands for ages 1-3 combined and puts every other bin two ages
  out. The residuals themselves were right; only the labels were positional.
  `age_length_bin` now gives the age or length each residual belongs to, and a
  new `accumulated` column marks the bins that carry a folded tail. A fleet
  without accumulation is unchanged.

* **`residuals()` reports composition residuals on the bins the likelihood
  fit.** Tail accumulation folds the tails into a boundary bin and restricts the
  likelihood to that window, but `residuals()` knew nothing about it -- so a
  fleet folding ages 1-2 into age 3 still returned residuals for ages 1 and 2,
  against a model that never fit them separately, and the boundary bin's
  observed proportion excluded the tail it had absorbed. Both the composition
  and the `"pearson"` attribute of `osa_residuals()` now fold first, so the
  Pearson and OSA views of a fleet describe the same bins. The residual formula
  is unchanged, and a fleet without accumulation is bit-identical. A new
  `Accumulated` column marks the bins carrying a folded tail.

* **`residuals()` no longer returns composition residuals for observations the
  model did not fit** (issue #108). A row with `Sample_size` 0 enters no
  likelihood -- the TMB guard is `Neff > 0` -- but the Pearson path reported one
  anyway, so a diagnostic plot showed residuals for data the model never used.
  On the 2025 GOA pollock assessment that was 531 of 1568 rows: every length
  composition (all carried at `Sample_size` 0) plus 8 unfitted age rows. The OSA
  path already applied the guard, which is why the two disagreed. Same rule now
  on both, for `comp` and `caal`.

* **The `"pearson"` attribute of an `osa_residuals()` result uses the same
  column names as the result itself.** It came from `residuals()`, which names
  columns in the data-sheet style (`Fleet_code`, `Year`, `Observed`), so one
  object carried two conventions and a reader had to know which half they were
  holding. Shared names do not mean a shared scale, and `?osa_residuals` now
  says so: the attribute's `observed`/`predicted` are proportions with the
  sample size alongside, since composition Pearson residuals are defined on
  proportions, where the OSA columns carry bin counts.

* **`sim_mod()` draws the survey index under the fleet's own
  `Index_distribution`.** Every fleet was drawn as an independent lognormal
  whatever it was set to, so an `MVN`/`MVNORM` survey was simulated on the wrong
  scale and without the correlation its covariance carries: a survey with an
  absolute sd of 20 and correlation 0.6 came back with a log-scale sd of 0.1 and
  no correlation. `Normal` was likewise drawn on the log scale. Nothing errored,
  so `self_test()` reported recovery against a data-generating process the
  likelihood never assumed. **Affects `self_test()`, `jitter()` and
  `run_mse(simulate_data = TRUE)` for MVN/MVNORM/Normal surveys only** --
  all-lognormal models draw exactly as before. A natural-scale draw can come out
  non-positive, which `data_check()` rejects; that now warns, since otherwise
  the refit simply fails and the run is reported as not converged.

* **A parallel worker that dies no longer aborts `retrospective()`, `jitter()`
  or `run_mse()`.** Each worker fits whole models, so running several at once can
  exhaust memory and the machine kills one. The call then stopped with
  `Error in unserialize(node$con) : error reading from connection`, discarding
  every peel or simulation finished up to that point. The tasks are now
  recomputed sequentially instead -- slower, but the call returns -- and the
  message suggests a lower `cores`. An error raised inside a worker still
  surfaces as itself. These functions parallelize differently on Windows than on
  macOS and Linux, so the same model can hit this on one and not the other.

* **Each `retrospective()` peel reports its own terminal year, so the peels show
  up in a plot.** Every peel carried the terminal year of the unpeeled model, and
  plots take their year axis from it, so each peel was drawn across the full
  series and they all lay on top of one another.
  `plot_biomass(retro$Rceattle_list)` now gives the usual fan with no extra
  arguments. Mohn's rho is unchanged.

  A peel still estimates the years it dropped -- they are its retrospective
  forecast, fit to the observed catch -- so `endyr` marks what the peel was fit
  through, not how far it was estimated. `incl_proj = TRUE` plots those years,
  and `data_list$endyr_full` gives the unpeeled terminal year where they end.

* **`?osa_residuals` no longer carries a broken link to `parallel::mclapply()`.**
  The `parallel` argument's cross-reference had been generated with a
  Windows-specific file anchor (`parallel:mcdummies`, where Windows documents the
  forking dummies); `mclapply` lives in `mclapply.Rd` on Unix, so the link
  resolved nowhere on macOS or Linux and `R CMD check` reported it as a missing
  link. It now points at the topic rather than a platform's file, so it resolves
  everywhere.

# Rceattle 5.5.1

*Released from `main` as 4.9.1, before this development line was merged back.
The two sections carry the same changes under different version numbers.*

## Bug fixes

* **The stock-recruit convergence check now runs under the Ianelli
  configuration.** `.check_stock_recruit()` was keyed on `srr_fun`, so with
  `srr_fun = 0` and `srr_pred_fun > 1` it returned no result at all -- the case
  where 5.5.0 made \eqn{\alpha} estimable against a steepness prior, and so
  exactly the case most in need of checking. It now reads `srr_pred_fun` there.
  Steepness is the complete test: for Beverton-Holt
  \eqn{h = \alpha \phi_0 / (4 + \alpha \phi_0)}, so \eqn{\alpha < 1/\phi_0} is
  precisely \eqn{h < 0.2}. The reported \eqn{R_0} is not checked in that
  configuration, because there it is the mean-recruitment level rather than a
  value the penalty curve implies.

* **`build_srr()` warns when a steepness is passed as `srr_est_mode = 0`'s
  \eqn{\alpha}.** `srr_prior` means steepness for `srr_est_mode` 2 and 3 but
  \eqn{\alpha} for mode 0. Since 5.5.0 fixes \eqn{\alpha} at that value under
  the Ianelli configuration, a steepness supplied by mistake is now pinned
  rather than merely used as a starting value, putting the curve under the
  replacement line with a negative implied \eqn{R_0}. A Beverton-Holt
  `srr_prior` in (0, 1) under mode 0 now warns and gives the conversion
  \eqn{\alpha = 4h / (\phi_0 (1 - h))}. It warns rather than errors because
  \eqn{1/\phi_0} is not known until the model is built, and a small
  \eqn{\alpha} is legitimate for a stock with a large \eqn{\phi_0}.

# Rceattle 5.5.0

*Released from `main` as 4.9.0, before this development line was merged back.
The two sections carry the same changes under different version numbers.*

## New features

* **`build_srr(srr_alpha_init = , srr_beta_init = )` sets the stock-recruit
  starting values.** The package defaults (\eqn{\alpha = e^3}, \eqn{\beta = 3})
  are placeholders that carry no knowledge of the stock's scale. \eqn{\beta} sets
  the density dependence in \eqn{R = \alpha S / (1 + \beta S)} and needs to be of
  order \eqn{(\alpha - 1/\phi_0)/R_0} -- typically \eqn{10^{-3}} or smaller for a
  stock in tonnes -- so the default sits three orders of magnitude away and
  Beverton-Holt fits returned `NA/NaN gradient evaluation` before reaching a
  first useful step. Both arguments take natural-scale values, one per species.
  For a Beverton-Holt seeded from steepness \eqn{h} and unfished spawning biomass
  per recruit \eqn{\phi_0}:

  ```r
  alpha <- 4 * h / (phi0 * (1 - h))
  build_srr(srr_fun = "BevertonHolt",
            srr_alpha_init = alpha,
            srr_beta_init  = (alpha - 1 / phi0) / R0)
  ```

## Bug fixes

* **`init_dev` is no longer declared a random effect under
  `initMode = "FreeParams"` (0).** `ceattle.cpp` applies the initial-deviate
  density `dnorm(init_dev, -sigma^2/2, R_sd)` only when
  `initMode > 1 && initMode != 5`, but `fit_mod()` added `init_dev` to the
  Laplace random block whenever any element was free. Under mode 0 that made the
  approximation integrate over an improper (flat) prior rather than estimate the
  initial age structure as fixed effects, which is what mode 0 means and what
  comparable platforms (WHAM's `N1_model = 0`, FIMS's `log_init_naa`) do.

  Only `initMode = 0` combined with `random_rec = TRUE` is affected, but **for
  that combination results move substantially** -- on a small test model the
  objective went from -91.6 to -83.3 and SSB by ~59%, because a flat prior and a
  fixed effect are different models. No assessment in the ecosystem uses that
  combination (every `initMode = 0` model sets `random_rec = FALSE`), so nothing
  in practice changes; refit if you do use it.

* **`build_params()` no longer seeds the stock-recruit \eqn{\alpha} from
  `srr_prior` when `srr_prior` is a steepness.** `srr_prior` is not the same
  quantity for every curve: for Ricker it is a prior on \eqn{\alpha}, but for
  Beverton-Holt it is a prior on **steepness** and must lie in (0, 1). Seeding
  \eqn{\alpha} with `log(steepness)` started the optimizer at a near-zero
  recruits-per-spawner. The seeding is now applied only where `srr_prior` is
  genuinely an \eqn{\alpha}, matching the prior gates in `ceattle.cpp`.

  **This changes results for Beverton-Holt models with `srr_est_mode` 2 or 3**,
  which are the only configurations where `srr_prior` is a steepness. Six of the
  24 `srr_est_mode` x `srr_pred_fun` combinations move; the other 18 are
  byte-identical. In the ecosystem that is BSAI Atka mackerel
  (`srr_est_mode = 2`, `srr_prior = 0.8`), where the starting \eqn{\alpha} goes
  from \eqn{\log(0.8) = -0.22} to the package default and the fit moves by
  ~6e-4 in relative SSB. Refit affected models rather than assuming the old
  numbers carry over.

* **Removed the spurious "alpha was not initialized to `srr_prior`" message.**
  `fit_mod()` copies `recFun$srr_prior` into the `data_list` only after the
  caller has already built parameters, so the message fired on the documented
  workflow and was not actionable.

* **A prior on steepness is now applied under the Ianelli configuration.**
  Steepness belongs to the stock-recruit curve, but it was only ever derived
  inside the `switch(srr_fun)` block, so with `srr_fun = 0` and
  `srr_pred_fun > 1` -- the AMAK/Ianelli setup, where the curve is estimated as a
  recruitment penalty -- it stayed at the mean-recruitment constant 0.99. The
  stock-recruit prior (`srr_est_mode` 2 or 3) is evaluated on steepness, so it
  was being applied to a constant: no gradient, no effect on \eqn{\alpha}, and a
  large fixed offset added to the objective. Steepness is now derived from the
  penalty curve's \eqn{\alpha} in that configuration. \eqn{R_0} is deliberately
  left alone, since recruitment there is \eqn{R_0 \exp(\mathrm{rec\_dev})}.

  **This changes results for any model combining `srr_fun = 0`,
  `srr_pred_fun` 2 / 3, and `srr_est_mode` 2 / 3** -- in the ecosystem that is
  BSAI Atka mackerel (`srr_prior = 0.8`, `srr_prior_sd = 0.0001`), whose
  steepness prior previously contributed a constant ~2.27e6 to the objective and
  did nothing else. It now constrains \eqn{\alpha} as intended. Refit and expect
  the reported objective to change substantially, mostly through the removal of
  that offset. Models with `srr_est_mode = 1` (no prior) are unaffected, as are
  models whose hindcast already uses the curve (`srr_fun > 1`), where steepness
  was always derived.

* **`srr_est_mode = 0` ("fix alpha to prior mean") now works for Beverton-Holt.**
  The gate in `fit_mod()` was Ricker-only, so a Beverton-Holt fit that supplied
  `inits` -- the normal warm-start workflow -- had \eqn{\alpha} mapped out by
  `build_map()` but never set to the prior mean, leaving it pinned at
  `build_params()`' placeholder \eqn{e^3}. `build_params()` and `fit_mod()` now
  share one rule (`.srr_prior_is_alpha()`) so the two paths cannot disagree.
  `build_map()` now keys the \eqn{\alpha} mapping on `srr_pred_fun`, the curve
  \eqn{\alpha} belongs to, so `srr_est_mode = 0` also fixes it under the Ianelli
  configuration. \eqn{R_0} stays keyed on `srr_fun` and remains estimated there.
  No model in the ecosystem currently uses `srr_est_mode = 0`.

* **`srr_alpha_init` / `srr_beta_init` survive a refit.** `.refit_like()`, which
  backs `retrospective()`, `jitter()`, `self_test()`, `profile()`, `run_mse()`,
  `remove_F()`, `sample_rec()` and `reweight_comps()`, rebuilds the
  stock-recruit specification from the stored `data_list` and now carries both
  starting values through. Note they override a supplied `inits`, so
  warm-starting from a previous fit's `estimated_params` restarts
  \eqn{\alpha}/\eqn{\beta} at the specified values rather than the fitted ones;
  leave them `NULL` to warm-start from the fit.

* **`build_srr()` validates the steepness prior instead of failing silently.**
  For a Beverton-Holt curve, `srr_est_mode` 2 and 3 put the prior on steepness,
  so `srr_prior` must lie in (0, 1) -- this is now checked. The beta prior
  (`srr_est_mode = 3`) additionally converts (mean, sd) to shape parameters by
  moments, which are positive only when \eqn{sd^2 < \mu(1-\mu)}; outside that
  range the shapes go negative and the prior is meaningless. It failed silently
  rather than loudly, because TMB's `lgamma`-based `dbeta` returns a finite
  value for negative shapes where R's `dbeta` returns `NaN`. Note the package
  default `srr_prior_sd = 1` is never valid for a steepness and now errors with
  the largest permissible value.

* **`fit$initial_params` now records the parameters the model actually started
  from.** It was captured before the blocks that overwrite `proj_F_prop`,
  `log_Ftarget`, `log_M1` and the stock-recruit alpha / beta, so it reported
  values no fit ever used. This was not only cosmetic: `retrospective()` and
  `jitter()` reuse `initial_params` as their refit starting values
  (`R/9-retro_and_jitter.R`). Those overrides are deterministic functions of the
  `data_list` / `HCR` and are re-applied on every refit, so no fit changes --
  confirmed bit-identical across `retrospective`, `jitter`, `self_test`,
  `profile`, `remove_F` and three `run_mse` variants via
  `tools/verify/verify-refit-like.R`.

* **`steepness` is reported for every year, not just the first.** It is declared
  `[nspp, nyrs]` but only column 1 was ever assigned, so
  `fit$quantities$steepness` was a value followed by structural zeros. It is now
  filled from the per-year alpha, which also generalizes to a time-varying
  recruitment linkage.

  `R0` is deliberately left as it was -- column 1 carries the curve-derived value
  under Beverton-Holt and the later columns the `exp(rec_pars)` level. Filling it
  per-year makes `R0` an AD function of alpha and beta across the whole series,
  which feeds the reference-point and projection recruitment calls and costs
  roughly 6x in fit time (a random-effects fit on a 30-year, 12-age model went
  from 14 s to 91 s, with an identical objective) for a quantity that
  Beverton-Holt recruitment does not use.

* **New `stock_recruit` convergence check.** Rceattle parameterizes
  Beverton-Holt by \eqn{\alpha} and \eqn{\beta} and derives steepness as
  \eqn{h = \alpha \phi_0/(4 + \alpha \phi_0)}, with no lower bound -- unlike
  WHAM, which builds \eqn{h} on (0.2, 1) by construction. Since `rec_pars`
  carries no default bounds, the optimizer can reach \eqn{\alpha < 1/\phi_0},
  where the stock cannot replace itself and the implied unfished recruitment
  \eqn{R_0 = (\alpha - 1/\phi_0)/\beta} turns **negative**, propagating into the
  initial age structure as `NaN`. `fit$convergence` now reports this as a `FAIL`
  naming the species and the offending steepness / \eqn{R_0}.

# Rceattle 5.4.0

## New features

* **`linkage_spec(integrate = FALSE)` estimates a random-effect linkage term as a
  penalized fixed effect.** A `~ (1|Year)` / `rw(1|Year)` / `ar1(1|Year)` term is
  normally integrated out by the Laplace approximation. Many reference
  assessments -- ADMB/AMAK, `goa_pk`, and Rceattle's own legacy
  `Time_varying_sel` / `Time_varying_q` switches -- instead treat such deviations
  as *penalized fixed effects*: the deviations sit in the objective as a plain
  penalty with a fixed SD and are not integrated. These are different models, not
  a reparametrization, so a `rw()` could not reproduce them. On the GOA pollock
  fishery selectivity walk, integrating moved the fit by ~119 jnll units and
  shifted SSB ~14%.

  ```r
  build_selectivity(linkages = list(
    slp_asc = linkage_spec(~ rw(1 | Year), fleet = "GOA_pollock_fishery",
                           init = list(sigma = 0.05), integrate = FALSE)))
  ```

  Permitted only with a fixed SD, and rejected with `observe = `. See
  `?linkage_spec` for the restrictions and
  `vignette("environmental-linkages-and-priors")` for when to prefer it.

* **The linkage table records where each deviation is stored (`re_pos`).**
  `re_index` is the deviation's slot in the *global* random-effect numbering, but
  the deviations are split across `beta_linkage_re` and `beta_linkage_re_pen`, so
  the position within the vector that actually holds one is a different number as
  soon as a model uses both treatments. `re_pos` gives that position, alongside
  `re_integrate` for which vector, so setting `inits` by hand is
  `beta_linkage_re_pen[re_pos + 1]` rather than a reconstruction of the split.
  Indexing a mixed model by `re_index` writes to the wrong deviation without
  erroring.

## Bug fixes

* **`osa_residuals()` falls back to serial when the parallel loop fails.**
  `TMB::oneStepPredict(parallel = TRUE)` forks via `mclapply`, and on some
  model/observation combinations a worker aborts rather than returning. The dead
  worker left an error object where a gradient belonged, which surfaced from deep
  inside TMB as `non-numeric argument to mathematical function` -- a message with
  no visible connection to the cause. The parallel loop is now retried serially,
  reporting that it did so, so the call returns residuals instead of failing.
  Seen on the GOA pollock model's index residuals; the other sources, and every
  model tested with a linkage on catchability, run in parallel unaffected.

* **A selectivity prior on a limb the fleet's curve does not have is now
  rejected.** A Logistic fleet has no descending limb and a DescendingLogistic
  fleet no ascending one, so those slots never reach selectivity-at-age -- they
  sit at their `build_params()` defaults. A prior on them was still added to the
  objective, shifting the reported likelihood by a constant that moved with an
  unrelated default while doing nothing to the fit. `fit_mod()` now names the
  fleet and parameter and says which limb to prior instead. This is how a
  reconciliation against another model picks up an unexplained offset: on the GOA
  pollock bridge it was worth 13.19 likelihood units.

* **Parameter bounds are applied to the parameter they were written for.**
  `fit_mod()` assembled `lower`/`upper` in `build_params()`'s parameter order,
  but TMB orders `obj$par` by the sequence the `PARAMETER_*` macros appear in the
  template -- and `build_params()` lists the linkage coefficients after `log_F`
  while `ceattle.cpp` declares them before it. Both vectors have the same length
  either way, so nothing downstream could notice: the box constraints landed on
  the wrong parameters. The orders disagree in four places, so this reaches well
  beyond linkages: `rec_dev`/`R_log_sd`, the growth block against
  `log_Flimit`...`log_F`, `index_q_beta`/`index_q_rho`, and
  `sel_curve_pen`/`sel_coff_dev`. **Fits of growth-estimating and
  linkage-carrying models can therefore change**, since bounds that previously
  bound the wrong parameter now bind the right one; refit them. The reference
  fits are unaffected (their misplaced bounds never bound). Bounds are now
  aligned against `names(obj$par)` after `MakeADFun()` and asserted non-`NA`,
  since `nlminb` accepts an `NA` bound silently and returns `objective = Inf`
  with `par = NA` rather than erroring.

* **`est_phase` no longer silently fails to fix a random-effect term.**
  `est_phase` reaches only `beta_linkage`, the fixed-effect coefficient vector; a
  random effect's deviations live in a separate vector it never touches, so
  `est_phase = 0` on a formula containing `(1|g)` / `rw()` / `ar1()` left every
  deviation estimated -- the opposite of the documented "fix at `init`" contract.
  It is now an error naming the supported alternatives (drop the term, or fix a
  small SD). Values above `1` remain inert for every linkage row.

# Rceattle 5.3.0

*Also covers the unreleased 5.2.0-5.2.4 development increments, whose entries
were folded into this section.*

## New features

* **`reweight_comps()` runs the iterative McAllister-Ianelli tuning loop.**
  Every fit reports the implied weights in
  `fleet_control$Comp_weights_mcallister`, but tuning on them has meant copying
  them across and refitting by hand, over and over. `reweight_comps(fit)` does
  that loop: it refits until the largest relative change in any weight falls
  below `tol`, returns the model fitted with the settled weights, and attaches
  the per-iteration history as `fit$reweight`. Dirichlet-multinomial fleets
  estimate their own weight, so they are named and skipped.

* **`intercept` is accepted wherever `` `(Intercept)` `` is.** The key naming a
  linkage's intercept came straight from `model.matrix()`, so setting an init,
  bounds or a prior on it meant writing `` list(`(Intercept)` = ...) `` --
  backticks, a capital I, and parentheses, with no error if any of it was wrong.
  `intercept` now means the same thing in `init`, `bounds` and `priors` on every
  process. The original spelling still works, and a covariate actually named
  `intercept` still refers to itself.

* **An `init`, `bounds` or `priors` key that names nothing is now an error.** A
  misspelled key was ignored without complaint, and the model went on to fit from
  a starting value, or under a prior, that was not the one asked for. Keys are
  now checked against the linkage's own design columns (plus `sigma` / `rho`
  where the structure has them); an unrecognized key stops the fit and lists the
  keys that were available.

* **An intercept `init` or `bounds` reaches the base parameter for all six
  processes.** Setting the intercept of a recruitment, natural-mortality or
  growth linkage set that parameter's starting value and bounds, but the same key
  on a catchability, selectivity or composition linkage was accepted and then
  dropped. All six now behave alike. Values are given on the parameter's natural
  scale and stored on whichever scale the parameter uses, so a catchability of
  0.8 or a von Bertalanffy K of 0.25 is written as given.

  Three cases cannot be honoured, and each now stops the fit with a message
  naming the linkage: a value of zero or less on a parameter held as a log; a
  natural-scale value for a `sel_inf` slot that holds a logit or a log instead of
  an inflection (DoubleNormal, LogisticPM); and a fleet that mirrors another's
  selectivity or catchability, where the value would be overwritten by the fleet
  leading the shared block.

* **`osa_residuals()` gains an `"ecov"` source** for the state-space covariate of
  a QAR1 catchability (`build_catchability(..., observe=)`, Rogers et al. 2024).
  It one-step-ahead residualizes the covariate observation with
  `TMB::oneStepPredict()`, first and against its own series, as in WHAM's
  `make_osa_residuals()`; the fit itself is unchanged. The AR1 latent is
  zero-mean, so standardize the covariate -- `build_catchability()` warns
  otherwise.


* **`linkage_spec()` selects species and fleets by name.** `species =` now accepts
  species names matching `data_list$spnames`, and `fleet =` accepts fleet names
  matching `fleet_control$Fleet_name`, in place of 1-based ids:

  ```r
  q = linkage_spec(~ ar1(1 | Year), by = ~ fleet, fleet = "Shelikof")
  M1 = linkage_spec(~ temp, by = ~ species, species = c("Pollock", "Cod"))
  ```

  Ids continue to work unchanged. Names are resolved when the model is assembled
  in `fit_mod()`; an unrecognized name -- or one that is not unique in
  `fleet_control` -- errors and lists the model's own names, whereas a
  `Fleet_code` that is wrong but in range silently attaches the linkage to a
  different fleet and still fits. Matching is exact after trimming whitespace, and
  a misspelled name is rejected even when the corresponding filter would be a
  no-op. Named specs round-trip through `save_config()` / `load_config()`.

## Performance

* **`run_mse()` refits the operating model only as far ahead as the next
  assessment needs.** Between assessments the operating model has to reach one
  assessment step past its terminal year -- far enough for the exploitable-biomass
  cap on the next TAC -- but it was rebuilt over the whole projection every time,
  and `projyr` sizes the AD tape, which costs roughly a fixed amount per model
  year. The refit now runs on the shortened horizon and the operating model is
  restored to the full projection afterwards; the final assessment keeps the full
  horizon, so the returned operating model is unchanged. Measured at a 9%
  saving on a Bering Sea multispecies MSE projecting to 2040, growing with the
  length of the projection. Single-species models see little change.

  **A run with `simulate_data = TRUE` will not reproduce results generated before
  this version**, because the operating model carries fewer projection rows while
  the refit runs and `sim_mod()` therefore draws a different number of random
  values. Runs remain fully reproducible from a given `seed` within a version.
  With `simulate_data = FALSE` results are unchanged to the bit.

## Bug fixes

* **`Time_varying_sel = "RandomWalkAscending"` no longer penalizes the descending
  limb it never estimates.** The random-walk selectivity penalty was gated on the
  selectivity *type* alone, so on a double-logistic fleet the descending-limb term
  was accumulated under mode 5 as well as mode 4. Mode 5 varies the ascending limb
  only, so with `random_sel = FALSE` those deviates sat at their init of 0 and the
  term contributed a pure constant -- on `GOA2018SS` (fleet 8, 42 years, one sex)
  exactly 113.459 objective units. Models with no mode-5 fleet -- including
  `BS2017SS` / `BS2017MS` -- are bit-identical.

  **A constant shift can still move the fit, and on `GOA2018SS` it does.**
  `nlminb` stops on an objective-*relative* tolerance, so an additive constant
  changes where it halts and which local optimum it reaches. The fit now lands on
  a converged optimum 52.9 units better (12868.005 rather than 12920.897 on the
  same surface), and the derived quantities move with it:
  **GOA Pacific cod SSB in the final projection year (2050) falls 14.1%**
  (395,627 to 339,897 t), while SSB, recruitment and F differ by up to 40%, 39%
  and 230% across the hindcast. **Anyone carrying GOA numbers forward should
  refit.**

  With `random_sel = TRUE` the removed term is **not** constant: `sel_dev_log_sd`
  is estimated for a mode-5 fleet, so removing the term moves the estimated
  selectivity-deviation SD upward and rescales shrinkage on the ascending walk.
  The direction is right -- penalizing a deviate pinned at 0 is a spurious
  `-log sigma` -- but it is a behavior change for such models.

  All four pinned reference objectives are restated: the reference recipe now
  Newton-polishes each fit to a stationary point (gradients ~1e-11) instead of
  stopping on `nlminb`'s relative tolerance, and `test-golden-regression.R`
  asserts convergence alongside each pinned value. Note this certifies
  convergence, not uniqueness -- use a jitter or multi-start check for that.
  `fit_control()`'s `newtonsteps = 0` default is unchanged; this affects the
  reference recipe only, not any user fit.

* **`retrospective()` and `run_mse()` run again on a model with an observed
  (state-space covariate) linkage.** The map for the QAR1 observation SD
  (`log_obs_sd_linkage`) was stored as a factor, where every other entry of
  `mapList` is a raw integer. Both diagnostics build a debug map for each peel /
  assessment, which maps the whole list out at once, and on a factor that step
  aborted with `argument "list" is missing, with no default` -- so every peel of
  such a model failed, reported under the parallel dispatch only as `5 nodes
  produced errors`. Models with no `linkage_spec(observe = )` group were never
  affected, and no fit changes: the map TMB receives is identical.

* **A diagnostic refit keeps the source model's bias adjustment.** `fit_mod()`
  resets `data_list`'s `bias_adjust_obs` / `bias_adjust_proc` and then re-applies
  the values from `fit_control()`, whose defaults are `TRUE`. `.refit_like()`
  built a fresh `fit_control()`, so a model fitted with bias adjustment off
  refitted with it back on -- worth ~880 jnll units on `BS2017SS`. Bias
  adjustment defines the likelihood rather than the optimizer, so this made
  `retrospective()`, `jitter()`, `self_test()`, `profile()`, `run_mse()`,
  `remove_F()` and `sample_rec()` compare two different objectives: a Mohn's rho
  or a jitter spread computed this way was not measuring what it reported. The
  resolved settings are now recovered from the `data_list`, which covers every
  caller. Only models fitted with a non-default `bias_adjust_*` are affected;
  `comp_offset` was already carried correctly and is now pinned by a test.

  The optimizer and `sdreport` fields -- `newtonsteps`, `use_gradient`,
  `rel_tol`, `nlminb_control`, `getJointPrecision` and the rest -- still fall
  back to their defaults in a refit. Those choose how the optimizer runs rather
  than what it optimizes, so they do not change the model, but a source fit that
  needed non-default optimizer settings will refit without them. `TMBfilename` is
  the exception worth naming: it selects the compiled template, so a model fitted
  against an alternate template refits against the bundled `ceattle` and is not
  the same likelihood. Refitting such a model through any of these diagnostics is
  not supported.

* **`sim_mod()` simulates with the model's own observation bias adjustment.** It
  hardcoded the `-sigma^2/2` offset on the simulated index and catch, which
  agreed with the estimator only because a refit used to force bias adjustment
  back on. With the refit now keeping the source setting, a model fitted at
  `bias_adjust_obs = FALSE` would be simulated with the offset and then fitted by
  a likelihood expecting none -- a systematic bias in scale, and so in
  catchability, that no number of simulations averages away. `self_test()` and
  `run_mse(simulate_data = TRUE)` on such a model would report that artifact as
  estimation error. The simulator now reads the flag the same way
  `residuals.Rceattle()` already did. Models on the default bias adjustment are
  unaffected.

* **A simulation whose unfished reference run fails keeps an `OM_no_F` entry.**
  `run_mse()` assigned the failure result with `sim_list$OM_no_F <- NULL`, which
  *removes* a list element rather than setting it to `NULL`, so the simulation
  came back without the name at all. `is.null()` tests were unaffected, but code
  keying on `names()` saw something different from what the documentation
  described. The element is now present and `NULL`, as stated.

* **A failed re-assessment now stops its simulation instead of being absorbed.**
  When an estimation-model refit failed, `run_mse()` could not assign to
  `em_use`, so it still held the *previous* year's assessment; the guard that
  followed tested `is.null(em_use)`, which such a failure never produces. The
  simulation carried on, stored that stale assessment under the failed year's
  name, handed its catch advice to the next iteration, and was reported as
  successful once a later iteration reset the failure flag. A simulation whose
  assessment fails is now stopped and marked like one whose operating-model
  refit fails. Any earlier run is worth re-checking: a failed assessment left no
  trace in `use_sim`.

* **A simulation whose unfished reference run fails is no longer counted as
  "did not collapse".** `remove_F()` is a separate fit and can fail on its own.
  The simulation is kept, since its assessments are the catch-advice record, but
  `OM_no_F` is then `NULL` -- and `sum(NULL < 1000) > 0` is `FALSE` in R, so
  `mse_summary()` scored it as a non-collapse while keeping it in the
  denominator. Because `remove_F()` fails preferentially at low stock sizes, that
  biased the collapse metrics low, in the direction that flatters a harvest rule.
  Such simulations are now excluded from the `OM no F` and `SSB Collapse from F`
  metrics -- numerator and denominator both -- which report `NA` when none has an
  unfished run; every other metric still uses them. Under `msmMode > 0` a single
  such simulation previously threw and took the whole summary with it.

* **`run_mse(dir = )` accepts an absolute path.** The save built its target as
  `file.path(getwd(), dir)`, which turns an absolute `dir` into
  `<cwd>/<abs path>`, so the directory was never created and the run died with
  "cannot open the connection". Relative paths were unaffected.

* **`sampling_period` is now honoured per fleet when more than one year advances
  per assessment.** `run_mse()` selected the newly sampled rows by testing the
  year set and the fleet set separately, which is the same as matching
  fleet-year pairs only while `years_include` spans a single year. With
  `assessment_period > 1` and fleets on different sampling periods, the two
  tests together admitted every fleet in every year of the step: a fleet
  surveyed every second year was given observations in its off years as well.
  The estimation model was therefore fit to more survey and composition data
  than the sampling design specifies, understating terminal-year uncertainty and
  making a data-poor design perform like a data-rich one. Rows are now matched on
  the fleet and year together. Runs with `assessment_period = 1` -- the default
  -- are unaffected; results move for any run combining a longer assessment
  period with per-fleet sampling periods.

* **`sample_rec(update_model = TRUE)` works on models carrying catchability,
  selectivity, or composition linkages.** `.refit_like()` reconstructs the
  `build_catchability()` / `build_selectivity()` / `build_composition()`
  specifications from the source model, but `sample_rec()` still re-invoked
  `fit_mod()` through its own hand-written copy of that block, and so omitted all
  three. Because `fit_mod()` treats those arguments as the source of truth for
  `q_linkages` / `sel_linkages` / `comp_linkages`, the rebuilt model had no
  linkage table and therefore no `beta_linkage` parameters, which no longer
  matched the parameters carried over from the source model. The `inits` check in
  `fit_mod()` caught that and stopped with a length-mismatch error, so
  `sample_rec(update_model = TRUE)` failed outright on any model using an
  environmental linkage or a prior on catchability, on a selectivity parameter,
  or on a Dirichlet-multinomial composition weight — it did not return quietly
  wrong dynamics. `sample_rec()` now routes through `.refit_like()` like every
  other refit.

  Routing through `.refit_like()` also means `sample_rec(update_model = TRUE)`
  now carries the source model's `random_q` and `random_sel` settings into the
  rebuild, where the hand-written block left them at `FALSE`. For a model that
  estimates time-varying catchability or selectivity as random effects, the
  rebuilt object now keeps those random effects instead of dropping them.


* **The "invalid switch" errors for `initMode` and `HCR` now name the offending
  value.** They previously read `Invalid 'initMode' specified: .` with the value
  dropped; they now echo what was passed (e.g. a since-renamed alias), alongside
  the list of accepted names and integer codes.


* **The diet Dirichlet-multinomial weight is no longer estimated where there is
  no diet likelihood.** `build_map()` freed `diet_comp_weights` for any predator
  with `Diet_distribution = "DirichletMultinomial"` whenever `msmMode > 0`, but
  the model only fits the diet composition for predators with
  `suitMode > 0` -- under empirical suitability (`suitMode = 0`) the diet
  proportions are taken as given. A multispecies run with empirical suitability
  therefore carried one free, wholly uninformed overdispersion parameter per
  predator, inflating the parameter count and risking a rank-deficient Hessian
  under `getsd = TRUE`. The weight is now estimated only where the diet
  composition is fit. Fits that did not use a DM diet are unchanged.

* **A composition-weighting prior on a fixed DM weight is now ignored instead of
  adding a constant.** A `theta_comp` / `theta_caal` / `theta_diet` prior from
  `build_composition(linkages = )` is prior-only, so when the map fixes the
  weight it targets -- a single-species fit, empirical suitability, an `"Off"`
  fleet, a fleet with no composition data, or a user-supplied `map` -- the prior
  contributed a constant that shifted the reported `jnll` without moving any
  estimate, making likelihoods non-comparable across configurations. Such priors
  are now dropped with a message. The linkage rows themselves are kept, so
  `beta_linkage` keeps its length and `inits` remain portable between fits that
  share one `compFun`. `linkage_spec()` still errors up front on a prior that
  can never apply to the data at hand (a non-DM `Comp_distribution` /
  `Diet_distribution`).

* **`Hake` selectivity no longer warns about a valid `Time_varying_sel = 0`.**
  `build_map_selectivity()` warned "Time_varying_sel for fleet N is not
  compatible ... Current value: Off" on every fleet that was not `"IID"` --
  including `"Off"` (0), which its own message and `data_check()` both list as
  valid for the `"Hake"` (Taylor et al. 2014) form. Because `fit_mod()` wraps
  `build_map()` in `suppressWarnings()`, the spurious warning surfaced only in
  the callers that build the map directly, notably `run_mse()` and
  `retrospective()`, where it appeared once per iteration. The map itself was
  always correct (no coefficient deviates are estimated for `"Off"`), so fits
  are unchanged. The warning now fires only for a mode the `"Hake"` form cannot
  represent, and names `'Off'` rather than the non-existent `'None'`.


* **The diagnostic refit paths preserve catchability / selectivity / composition
  linkages.** `.refit_like()` -- the engine behind `run_mse()`,
  `retrospective()`, `jitter()`, `self_test()`, `profile()`, and `remove_F()` --
  rebuilt the recruitment, M1, and growth specifications from the source model
  but silently dropped the `build_catchability()` / `build_selectivity()` /
  `build_composition()` linkages. On a model using any of them (e.g. a
  Dirichlet-multinomial composition-weight prior), the reconstructed model
  lacked the linkages' `beta_linkage` parameters while the warm-start values
  still carried them, and `TMB::MakeADFun()` segfaulted on the length mismatch.
  The three linkage constructors are now rebuilt from the stored `data_list`
  fields (`q_linkages` / `sel_linkages` / `comp_linkages`), so these models run
  through the MSE / retrospective / jitter loops.

* **`run_mse()` extends the selectivity-deviation parameters at the model's own
  sex dimension.** The projection-year extension of `log_sel_slp_dev` /
  `sel_inf_dev` / `sel_coff_dev` hardcoded the sex dimension to 2, which for a
  single-sex model both recycled the fitted values into a phantom second sex and
  produced a warm-start whose length disagreed with the parameter template. Now
  taken from the fitted arrays (`max_sex`).

* **`fit_mod()` raises a clear error instead of crashing on mismatched
  `inits`.** When the supplied `inits` omit a parameter the model declares, or
  any parameter's length disagrees with the model implied by `data_list`
  (a dropped linkage, or a warm start not extended to a later `endyr`),
  `fit_mod()` now stops with a message naming the parameter and its lengths,
  rather than letting `TMB::MakeADFun()` segfault in `getParameterOrder()`.

* **The diagnostic refit paths preserve the `random_q` / `random_sel` switches.**
  These are now stored on `data_list` (like `random_rec`) and forwarded by
  `.refit_like()`, so a random-catchability or random-selectivity model keeps its
  random-effect structure across MSE / retrospective / jitter refits instead of
  silently reverting to fixed effects.

* **`build_params()` no longer pads `init_dev` behind a string comparison.** The
  `-999` padding of the unused `init_dev` columns was partly gated on
  `initMode > 0`, evaluated against the canonical `initMode` *string* that
  `switch_check()` has already resolved — `"FreeParams" > "0"` is `TRUE`, so the
  gate was always open. The padding is unconditional now, which is what the C++
  requires (every mode reads only columns `1:(nages - 1)`). No fit changes:
  `build_params()` is bit-identical for all six modes on `BS2017SS` and
  `BS2017MS`, from both string and integer input. This was the only place whose
  behavior depended on the lexicographic value of a switch alias, so an alias
  beginning with a digit can no longer flip it.

## Behavior changes

* **Composition weights now warm-start from `inits` like every other
  parameter.** `fit_mod()` re-read `comp_weights`, `caal_weights` and
  `diet_comp_weights` from their `fleet_control` columns on every fit, even when
  `inits` was supplied, so a weight handed to a refit was discarded. These are
  estimated parameters -- the Dirichlet-multinomial likelihood fits them, which
  is why they can carry a prior -- and their columns are starting values, read
  when a model is built from scratch. A refit now keeps the estimate it was
  given, which is what makes iterative reweighting possible and what every other
  parameter already did.

  Which models this moves depends on whether the weight is estimated. Under a
  multinomial likelihood it is a fixed multiplier that never leaves its column
  value, so those models -- including the bundled examples -- are unaffected, and
  an MSE over them is bit-identical. Under a Dirichlet-multinomial it is
  estimated and can move a long way: fitting `BS2017SS` with DM composition puts
  the weights between -0.6 and 12.7 (log scale), starting from the default column
  value of 1. **Results therefore move for any DM model** under every
  `run_mse()`, `retrospective()`, `jitter()`, `self_test()` and `profile()`
  refit, each of which previously discarded the fitted estimate and restarted
  from the column.

  Editing a `Comp_weights` column and re-fitting from an existing fit no longer
  has an effect. Either build the model afresh (`inits = NULL`), or set the
  parameter the fit actually reads -- `inits$comp_weights[flt] <- w` (likewise
  `caal_weights` / `diet_comp_weights`) -- which is what `reweight_comps()` does.
  `fit_mod()` now warns when supplied `inits` disagree with the columns, so an
  edit that would once have taken effect silently is no longer silent.


* **A simulation that does not run to completion now returns only a marker.**
  `run_mse()` previously returned the partially advanced operating and estimation
  models alongside `use_sim = FALSE`. Those models describe a state the
  simulation never reached, and nothing downstream filtered on `use_sim`, so they
  could be averaged into performance metrics as though the simulation had
  finished. A failed simulation is now `list(use_sim = FALSE, failure = ...)`
  with no models attached, and `mse_summary()` drops such simulations up front
  with a warning naming them, or errors if none completed. Code that reached into
  `mse$Sim_n$OM`, `$EM` or `$OM_no_F` without checking `use_sim` should now check
  it -- filter first with
  `mse <- mse[vapply(mse, function(x) isTRUE(x$use_sim), logical(1))]`.

## Documentation

* **Help pages, code comments, and vignettes were rewritten for clarity.** No
  code or model behavior changed.

* **Four help pages that disagreed with the code were corrected.** `jitter()`
  now states that starting values are perturbed around the model's initial
  (pre-fit) parameters, not the fitted values; the `build_params()` and
  `TMBphase()` `@return` descriptions were fixed (they had described a map object
  and standard errors respectively); and a mislabeled comment now correctly
  attributes the time-varying survey catchability SD.

# Rceattle 5.1.1

## New features

* **`Diet_distribution` accepts string aliases.** Like `Comp_distribution` and
  `CAAL_distribution`, the per-predator diet composition likelihood switch now takes
  `"Multinomial"` (0) or `"DirichletMultinomial"` (1) in addition to the integer code
  -- previously it accepted only the integer. The diet likelihood has no AFSC variant,
  so `"MultinomialAFSC"` / `-1` is rejected with a clear message.

* **Clearer spec-tree printing for data objects and fits.** `print()`/`summary()` of an
  `Rceattle_data` object or a fitted model now shows each switch by its string alias
  (e.g. `srr_fun = mean`, `M1_model = fixed`, `msmMode = MSVPA`) rather than its raw
  integer code, so the outline reads the same even if a code is renumbered. The fleet
  section is laid out as aligned, `│`-separated columns with `sel:` / `q:` labels and a
  `[shared: ...]` tag on mirrored selectivity/catchability blocks.

## Bug fixes

* **Optional-input default messages fire only when the model uses the input.** The
  "not specified, assuming ..." messages for `alpha_wt_len`, `beta_wt_len`, and
  `Selectivity_dimension` now print only when growth is estimated (`growth_model > 0`);
  `CAAL_distribution` / `CAAL_weights` only when the model has CAAL data; and
  `Sel_norm_bin_upper` only when a fleet normalizes selectivity at a specific bin. The
  default value is still applied silently in every case -- this only stops the messages
  nagging about inputs the configuration never consumes (e.g. weight-length parameters
  in a single-species age-based model).

# Rceattle 5.1.0

## New features

* **`linkage_spec()` fills in `by` from the process it is attached to.** Omitting
  `by` now defaults to the base stratum of the process -- `~ fleet` for
  catchability, selectivity, and the fleet composition weights (`theta_comp` /
  `theta_caal`), and `~ species` for recruitment, M, growth, and the diet weight
  (`theta_diet`) -- so a base-case linkage no longer needs an explicit
  `by = ~ fleet` / `by = ~ species`. An explicit `by` (a formula, or `NULL` for a
  single shared coefficient) is always kept as given, so code that passes `by`
  explicitly -- as the bundled examples and real assessments do -- is unchanged.
  A linkage that *omitted* `by` on a fleet-keyed process does change: it now spans
  every fleet (with a separate coefficient each) rather than collapsing to fleet 1.
  In particular a `~ 1` selectivity/catchability **prior** now applies per fleet, so
  name the target fleet(s) via `linkage_spec(fleet = ...)` -- and note it will error
  on a model with mirrored selectivity fleets unless a fleet is named. Restrict a
  covariate or prior to the fleets you mean with the `fleet` argument. When `by` is
  omitted on a catchability/selectivity linkage, `fit_mod()` now prints a one-time
  message noting that it spans every eligible fleet, so the per-fleet expansion is
  not a surprise.

## Bug fixes

* **A catchability linkage on a fleet with no survey index now errors.** A `q`
  linkage only makes sense on a fleet whose catchability is estimated. Fleets that
  hold q fixed or solve it analytically were already rejected; a fleet with **no
  index data** (so `Catchability` is `NA` and there is no q to estimate) was not, so
  an omitted-`by` linkage silently attached meaningless q coefficients to, e.g.,
  fishery fleets. Such fleets are now rejected with the same clear error, naming them
  and pointing to `linkage_spec(fleet = ...)` to restrict the linkage.

# Rceattle 5.0.6

## Bug fixes

* **A `right_floor` selectivity linkage now applies its offset on the logit scale.**
  For a double-normal fleet, `sel_inf(1)` parameterises the right-tail floor as
  `logit(right_floor)`, but a log-link linkage offset was applied by multiplying the
  already-transformed probability, which could push `right_floor` above 1 and corrupt the
  descending limb. The offset now acts on the logit (additive and multiplicative), like the
  peak and the other inflection slots, so `right_floor` stays in `[0, 1]` for any offset.
  A fleet with no `right_floor` linkage is bit-identical (the four reference models are
  unchanged).

# Rceattle 5.0.5

## Bug fixes

* **Random-effect deviate SDs are now held during phasing.** The 5.0.4 phased
  random-effect warm-up (deviates estimated as penalized fixed effects before the Laplace
  step) held every process's variance hyperparameter fixed except three observation/deviate
  SDs: `log_sigma_linkage` and `log_obs_sd_linkage` (the `~ (1 | Year)` / `rw()` / `ar1()`
  and Rogers-2024 QAR1 SDs) and `index_q_log_sd` (the legacy `Catchability = "AR1"` q
  observation SD). Left free, each was estimated jointly with its own penalized deviates and
  could collapse toward 0, so a *phased* fit silently converged to the random-effect-*off*
  solution -- the estimated time-variation vanished -- instead of the correct optimum. All
  three are now held during phasing like every other process SD. Fits without these random
  effects are unaffected (the four reference models stay bit-identical).

## Internal

* **The TMB model is now `src/TMB/ceattle.cpp`** (was `ceattle_v01_11.cpp`), and its
  compiled DLL is `ceattle`. Internal rename only -- the model and results are
  unchanged (the four reference models stay bit-identical), and `fit_mod()` loads the
  renamed DLL automatically. A fit built with a custom `fit_control(TMBfilename = ...)`
  is unaffected. Caveat: an object fitted by an earlier version stores
  `TMBfilename = "ceattle_v01_11"`, so re-running `osa_residuals()` on such a saved fit
  needs it re-fit under this version first (the old DLL name is no longer built).

## Documentation

* Moved the `?build_selectivity` prior example into `@examples` (it previously rendered
  a broken `\verb{}`); corrected the `?build_hcr` PFMC notation to use the normal quantile
  `qnorm()` rather than `Phi` (which denotes the CDF); and regrouped the 4.10.0 linkage
  additions under New features. Minor vignette / README fixes.

# Rceattle 5.0.4

## Bug fixes

* **Random-effect models phase without a flat-field NaN.** During phasing, the
  deviates of a random effect (recruitment, selectivity, catchability, ...) are now
  estimated as *penalized fixed effects* with their variance / correlation
  hyperparameters held, and the Laplace approximation is switched on only in the
  final hindcast fit. Previously the deviates were integrated out from the first
  phase, so a random field started flat (all-zero); for a 2D AR1 selectivity
  (`random_sel = TRUE`), integrating a flat field gives a `NaN` marginal objective
  and the fit could not start without an external warm start. Phasing the deviates
  as penalized effects first builds up realistic non-zero values, so the final
  Laplace begins near the mode. Gated on the presence of random effects, so fits
  without any are unchanged (the four reference models remain bit-identical).

# Rceattle 5.0.3

## New features

* **OSA residuals now support every survey-index likelihood family.**
  `osa_residuals()` previously produced one-step-ahead residuals only for the
  lognormal IID index (`Index_distribution` `"Lognormal"`). It now covers all
  families: natural-scale `"Normal"` residualizes as an independent normal on the
  natural scale, and the correlated covariance families `"MVN"` / `"MVNORM"` are
  whitened by the lower Cholesky of the fleet's survey covariance
  `Sigma = L L'`, so the residuals are the multivariate-Gaussian one-step-ahead
  innovations `L^-1 (obs - q*pred)` — the closed form `TMB::oneStepPredict()`
  reproduces for a Gaussian block (Thygesen et al. 2017). This supersedes the
  5.0.2 exclusion of these fleets. The ordinary model fit is unchanged (the
  whitening applies only to the post-hoc OSA observation vector).

* **String aliases for the `fit_mod()` mode switches.** `estimateMode` and
  `msmMode` now accept readable names alongside their integer codes:
  `estimateMode = "Estimate"` (0), `"Hindcast"` (1), `"Projection"` (2),
  `"DebugBuild"` (3), `"DebugOptimize"` (4); `msmMode = "SingleSpecies"` (0),
  `"MSVPA"` (1), `"TypeIIIMSVPA"` (2). The integer codes still work, and an
  unrecognized string errors with the valid options listed.

# Rceattle 5.0.2

## Bug fixes

* **OSA residuals no longer silently mis-residualize non-lognormal index fleets.**
  `osa_residuals()` / `build_osa_data()` previously laid MVN-covariance
  (`Index_distribution` "MVN"/"MVNORM") and natural-scale "Normal" survey
  observations into the one-step-ahead observation vector as if they were
  independent lognormals, even though the C++ likelihood for those families does
  not read the OSA observation vector — producing invalid residuals with no
  warning. Those fleets are now excluded from the OSA residuals with a warning
  naming them; lognormal (IID) index fleets in the same model are unaffected.
  (Covariance-aware OSA residuals for the MVN families are a separate follow-up.)

# Rceattle 5.0.1

## New features

* **OSA residuals now support composition tail accumulation.** When a fleet folds
  its age/length composition tails via `Comp_accum_young` / `Comp_accum_old` (AFSC
  `ac_yng`/`ac_old`), `osa_residuals()` now residualizes the folded composition
  that was actually fit -- one residual per fitted (folded) bin per sex block,
  rather than refusing the combination. `build_osa_data()` applies the same
  per-sex-block young/old fold to the OSA observation vector, so the one-step-ahead
  decomposition matches the fitted composition likelihood bin-for-bin. Fleets
  without accumulation are unchanged.

# Rceattle 5.0.0

*Also covers the unreleased 4.14.0 development increment, whose entries were
folded into this section.*

## Breaking changes

* **`mse_summary()` now returns a per-entity list instead of one stacked
  data.frame.** The result is `list(species, fleet, total, meta)`: `species`
  (one row per species, keyed by `Species`) holds the conservation/status
  metrics; `fleet` (one row per fishery fleet) holds average catch / catch IAV
  / P(Closed); `total` is the across-fleet totals; `meta` carries run
  provenance. This removes the NA padding and the ambiguity of the old shape,
  where species, fleet, and "All" rows were stacked in one frame and a column
  like `Average Catch` meant different things per row. The metric values are
  unchanged. Update `summ$"<metric>"` call sites to `summ$species$"<metric>"`
  or `summ$fleet$"<metric>"`.

## Bug fixes

* Fixed the `ConstantFSPR` harvest control rule (`HCR = 4`) applying `Fmult`
  twice, which set the projected F to `Ftarget * Fmult^2`. It is now
  `Ftarget * Fmult`. Only affects projections with `Fmult != 1` (the default
  `Fmult = 1` is unchanged).
* Corrected the documented reference-point formulas in `?build_hcr`: the SESSF
  (Tier 1) and NPFMC (Tier 3) fishing-mortality ramps and the PFMC 40-10
  buffer (the normal quantile function, not the CDF).
* Restored `plot_logindex()` as a deprecated shim for `plot_index(log = TRUE)`.
  It forwards to `plot_index()` with a `.Deprecated()` notice.

## New features

* **Plus-group SD-at-age convention (`build_growth(sd_plus_group=)`).** For
  estimated growth, choose how the oldest age class's standard deviation of
  length-at-age is set: `"WHAM"` (default) pins it to the upper anchor
  `exp(sd_Linf)` (the WHAM SDAA convention, unchanged from prior releases);
  `"SS3"` instead interpolates it by length like any interior age. Per-species
  (scalar or length-`nspp`), round-trips through `save_config()`/`load_config()`,
  and affects only `growth_model > 0` fits. Existing models are bit-identical
  under the default.

# Rceattle 4.13.1

## Bug fixes

* The NA-handling diagnostic in `check_caal_data()` now reports "CAAL data
  have NAs ..." instead of mislabeling them as "Composition data" (the
  message was copy-pasted from `check_composition_data()`). Message text only;
  no change to the data or fit.

* `model_average(uncertainty = TRUE)` no longer opens a plot device partway
  through the computation. A stray `plot_ssb()` call fired before the averaged
  `sdrep$sd` was populated, so it drew stale, pre-averaging confidence
  intervals as an unrequested side effect; the averaged object it returns is
  unchanged.

# Rceattle 4.13.0

## New features

* **Save / load a run configuration (`save_config()` / `load_config()`).** A full
  run configuration -- the [model_config()] structure plus the estimation controls
  and `fit_control()` bundle -- round-trips to a documented, git-diffable YAML file
  (each field carries its doc string as a comment; fields at their default are
  omitted so two runs diff to only their real differences). Apply a saved config
  with `fit_mod(data_list, config = load_config("run.yaml"))`; it fills only the
  arguments the caller did not pass (an explicit argument always wins), and every
  fit now records the run configuration it used (`fit$run_config`, also reachable
  via `run_config(fit)`). Follows SAM's `saveConf`/`loadConf` and the
  config-separate-from-weights idiom; adds a `yaml` dependency.

* **Composition tail accumulation (AFSC `ac_yng` / `ac_old`).** Two new
  `fleet_control` columns, `Comp_accum_young` and `Comp_accum_old`, fold the
  young and old tails of a fleet's age/length composition into a boundary bin
  before the composition likelihood (per sex block for joint-sex comps), for
  every composition family including the default `MultinomialAFSC`. They are
  1-based bin ordinals on the fleet's composition dimension; leaving them unset
  (or `young = 1` / `old` at the last bin) applies no accumulation, so every
  existing model is bit-identical. `data_check()` rejects out-of-range or
  inverted (`young > old`) bins, and OSA residuals are not available for a fleet
  with active accumulation (the residuals are built on the un-accumulated bins).

## Bug fixes

* **`self_test()` no longer errors with `object 'getsd' not found`.** Its
  per-simulation refit closure referenced a `getsd` value that the function --
  unlike `retrospective()`, `jitter()`, and `profile()` -- never defined, so
  every simulation died on both the sequential and the parallel-cluster path.
  `self_test()` now takes a `getsd` argument (default `NULL`, inheriting the
  input model's setting) like its sibling functions.

# Rceattle 4.12.0

## New features

* **Intuitive `fleet_control` / data column names.** The workbook and data-list
  columns were renamed to say what they hold, each keeping a back-compat alias
  that is upgraded on read, so existing data lists, `.rda`, workbooks, and
  scripts keep working (with a one-time deprecation message) and fit
  identically: `Q_index`/`Q_init`/`Q_sd_prior` →
  `Catchability_index`/`Catchability_init`/`Catchability_prior_sd`;
  `Comp_loglike`/`CAAL_loglike`/`Index_loglike` →
  `Comp_distribution`/`CAAL_distribution`/`Index_distribution`;
  `Sel_norm_bin1`/`Sel_norm_bin2` → `Sel_norm_bin`/`Sel_norm_bin_upper`;
  `weight1_Numbers2` → `Observation_units`; `proj_F_prop` → `Proj_F_proportion`;
  the control scalar `sigma_rec_prior` → `sigma_rec`; the bioenergetics
  `Diet_loglike` → `Diet_distribution`; and the post-fit output column
  `Est_weights_mcallister` → `Comp_weights_mcallister` (the old name is still
  written for downstream readers).

* **`write_template()`.** A new export that writes a minimal, structurally
  complete single-species starter workbook on the canonical column names; it
  round-trips through `read_data()` and builds under `fit_mod(estimateMode = 3)`.

* **Single-source workbook schema.** The column dictionary now lives once in the
  package (`R/0-column_schema.R`) and drives `switch_check()` defaults,
  `write_data()` object order, the embedded `meta_data` documentation sheet, and
  the roxygen field dictionary, which are kept in sync by guard tests. ~16
  previously used-but-undocumented columns are now documented.

## Bug fixes

* **`read_data()` names the offending sheet or cell instead of failing
  cryptically.** A required sheet (`control`, `fleet_control`) is named when
  absent, and optional sheets are skipped, so a minimal single-species workbook
  reads cleanly. A non-numeric cell in the control or bioenergetics rows is named
  rather than silently becoming `NA`. `rearrange_data()` likewise fails clearly
  on a malformed `fleet_control` rather than via a cryptic `dplyr` error.

* **Removed the dead accumulation-age feature.** `Accumulation_age_lower/upper`
  were only range-validated and never applied to composition data; they are
  dropped (with a soft-deprecation message if an old workbook carries them).

# Rceattle 4.11.0

## New features

* **`initMode = "OffsetEquilibrium"` (5).** A new initial-age-structure mode: an
  F = 0 equilibrium seeded by the *first-year* recruitment
  (`R_init * exp(rec_dev[year 1])`) decayed by residual natural mortality `M1`,
  with initial deviates turned off and no init-dev penalty. Unlike
  `initMode = "Equilibrium"`, which sits at the initial equilibrium recruitment
  `R_init`, the initial numbers track the year-1 recruitment deviation, matching
  the AFSC GOA pollock (Cole Monnahan) convention — the first-year cohort and
  the initial numbers-at-age share one deviation. The name refers to that
  offset; the equilibrium itself is unfished. Note the deviation enters as a
  single scalar on every age, so it moves the level of the initial numbers-at-age
  without changing their proportions-at-age. Every other mode is unchanged
  (golden references bit-identical).

* **Priors on base selectivity parameters through the linkage grammar.** A
  selectivity linkage with an intercept-only formula and a `priors` entry now
  places a prior on the base logistic parameter that carries the level — the
  ascending/descending slope (`slp_asc` / `slp_desc`, on `log_sel_slp`, log
  scale) or inflection (`inf_asc` / `inf_desc`, on `sel_inf`, natural scale),
  e.g. `build_selectivity(linkages = list(inf_asc = linkage_spec(~ 1, priors =
  list(\`(Intercept)\` = normal(0, 3)))))`. Use `lognormal()` for the log-scale
  slopes and `normal()` for the natural-scale inflections. This mirrors the
  prior-only `build_composition()` path and reproduces the AFSC GOA pollock
  selectivity priors. (Previously such an intercept prior was silently inert for
  selectivity — the `sel` process was missing from the prior re-targeting.)

* **Per-predator suitability reference years** (`suit_styr` / `suit_endyr`).
  These `fit_mod()` arguments (and the underlying `data_list` fields) now accept
  a vector of length `nspp` so each predator can average its suitability over a
  different set of years — e.g. a California Current model with hake 1980–2019,
  arrowtooth 2013–2018, and sablefish 2005–2008. A single scalar is recycled to
  every predator, exactly reproducing the previous global-window behavior
  (`BS2017MS` and the golden references are unchanged to numerical tolerance).
  Internally `suit_styr` / `suit_endyr` became per-predator `DATA_IVECTOR`s and
  the suitability-averaging and stomach-content prediction loops index them by
  predator.

* **`plot_diet_comp2()` and `plot_diet_comp1()`.** `plot_diet_comp2()` adds
  aggregation-aware diet-composition diagnostics (line plots when prey- or
  predator-age is aggregated, dodged bars when both are, bubble grids when fully
  disaggregated), built on `residuals(source = "diet")`. `plot_diet_comp1()` is
  an alias of `plot_diet_comp()` (the bubble/grid diagnostic).

## Bug fixes

* **`read_data()` tolerates trailing empty age columns in `NByageFixed`.** Older
  writers pad the fixed numbers-at-age sheet to a wider age range than
  `max(nages)`; the all-`NA` trailing columns are now dropped on read instead of
  tripping the `data_check()` column-count validation.

* **Parallel workers now run the in-session package.** `run_mse()`,
  `retrospective()`, and `jitter()` use a FORK cluster on non-Windows platforms,
  so workers run the parent's loaded namespace via copy-on-write. This fixes
  silently running a stale *installed* package on the workers during
  `pkgload::load_all()` development sessions, and removes the per-worker startup
  cost. Parallel and serial runs stay bit-identical for a given seed; Windows
  keeps the previous PSOCK path as a fallback.

# Rceattle 4.10.0

## New features

* **Random-effect linkages (`~ (1 | group)`).** A linkage formula may now carry
  an IID random effect, so a parameter can vary year to year (or over any
  grouping) as a set of deviations damped by an estimated
  `N(0, sigma)` density, rather than fixed covariate slopes:

    ```r
    fit_mod(...,
      qFun = build_catchability(linkages = list(
               q = linkage_spec(~ (1 | Year), by = ~ fleet))))
    ```

  Each level of the grouping variable gets one deviation (`beta_linkage_re`);
  each distinct group estimates its own log-SD (`log_sigma_linkage`). The
  density is reported in the new `jnll_comp` row *"Linkage random effects"*.
  The deviation SD is routed through the same `linkage_spec()` arguments as
  every other parameter: `init = list(sigma = v)` **fixes** it at an input
  value (reproducing the legacy `Time_varying_*_sd_prior` fixed input), and
  `priors = list(sigma = lognormal(...))` places a **prior** on it and
  estimates it — the first prior on a deviation SD anywhere in the model.

* **`data_requirements()` — see which inputs a model configuration needs.** A new
  exported reader reports, for a given model spec, which top-level data inputs are
  **Required**, **Optional** (used if supplied, otherwise default-filled by
  `clean_data()`), or **Ignored** (not consulted because the feature is switched
  off) — the same conditions `data_check()` enforces at fit time, surfaced up
  front instead of buried in the validator:

    ```r
    data_requirements(msmMode = 1)          # multispecies: diet/ration/bioenergetics Required
    data_requirements(BS2017SS, msmMode = 0) # preview a data object as single-species
    ```

  It accepts either an existing (possibly partial) data list or the convenience
  switch arguments; an explicit switch argument overrides the data list's stored
  value (matching `fit_mod()` precedence). Internally, `data_check()`'s conditional
  presence-requirement gates were refactored to consume one declarative
  requirement table, so the reader and the validator can never drift apart.

* **`build_data()` — assemble a data list in R.** A code-first constructor
  complementing `read_data()`: supply only the blocks a model uses and the
  optional blocks a single-species model does not need are default-filled by
  `clean_data()`. Three combinable entry points cover the workflows real
  assessments use:

    ```r
    build_data(base = BS2017SS, projyr = 2060)          # copy-and-edit a dataset
    build_data(file = "model.xlsx", fleet_control = fc) # read a workbook, override
    build_data(nspp = 1, styr = 1977, endyr = 2023,     # assemble from blocks
               fleet_control = fc, catch_data = catch, index_data = survey)
    ```

  Overrides are checked against the recognized schema, so a typo (`maturty`) is
  caught at construction with a suggestion rather than surfacing later in a fit;
  legacy names (`fsh_biom`, `srv_biom`, `wt`, `pmature`, `Pyrs`) are mapped to
  their canonical equivalents. Validation is deferred to `data_check()` at fit
  time (one source of truth); `build_data()` runs only a light presence
  pre-check so a missing required block is reported early. The result is the
  same bare list `read_data()` returns and round-trips through `write_data()`
  unchanged — a `build_data(base = X)` object fits bit-identically to `X`.

* **`model_config()` — a model configuration that travels with the data.** The
  model-structure arguments of `fit_mod()` (`msmMode`, `initMode`, the HCR and
  the `build_*()` process specifications) can now be bundled into a slot on the
  data list, so a data object records how it is meant to be fit:

    ```r
    dat <- build_data(base = BS2017MS, model_config = model_config(msmMode = 1))
    fit_mod(dat)                     # fits as multispecies without passing msmMode
    ```

  `fit_mod()`'s signature and defaults are unchanged; when a data list carries a
  `model_config`, `fit_mod()` reads each field only for arguments the caller did
  **not** pass (detected with `missing()`), and an explicitly-passed argument
  always overrides the slot. With no slot present the behavior is byte-identical
  (a `BS2017SS` fit is bit-identical). A call that passes an argument — even at
  its default — overrides the slot, so omit the argument to let the configuration
  take effect. The slot is code-side structure, not a workbook sheet, so it does
  not persist through a `write_data()`/`read_data()` round-trip (a warning fires).

* **Spec-tree `print()` / `summary()` for data objects and fits.** A
  `build_data()` object now carries the class `"Rceattle_data"` and prints as an
  indented specification tree — dimensions → fleets (with their selectivity /
  catchability forms and mirroring) → configured processes → active linkages →
  any attached `model_config()` — instead of dumping the ~40-element list. The
  same tree is shown by `print()` on a fitted model above its fit statistics, so
  "read 600 lines of switch tables" becomes "read the printout". The class is a
  thin tag: every consumer still treats the object as a plain list, so it does
  not change a fit.

* **Random-walk linkages (`rw(1 | group)`).** A deviation that follows a random
  walk — `N(0, sigma)` on successive first differences, the first deviate
  pinned so the walk's level stays with the base parameter — reproducing the
  Dorn-style `Time_varying_* = "RandomWalk"` process through the grammar. The
  grouping variable must be numeric (a real elapsed-time lag).

* **AR1 linkages (`ar1(1 | group)`).** A stationary first-order autoregressive
  deviation, `SCALE(AR1(rho), sigma)` with `sigma` the marginal SD and `rho`
  the correlation (the glmmTMB convention, and the same form as the Rogers et
  al. (2024) QAR1 catchability process). The correlation is routed through the
  same grammar as `sigma`: `init = list(rho = 0.7)` fixes it,
  `priors = list(rho = normal(0, 0.3))` places a prior on it, else it is
  estimated free. Reduces to the IID density at `rho = 0`.

* **State-space environmental covariate (Rogers et al. 2024 QAR1).** An
  `ar1(1 | Year)` term can be turned into a measured latent covariate via
  `observe` / `obs_sd`: the AR1 latent is observed as an `env_data` column with
  a fixed measurement SD, and enters the linked parameter through an estimated
  effect size. This reproduces the Rogers-2024 QAR1 catchability model
  (`Estimate_q = 6`) through the grammar:

    ```r
    build_catchability(linkages = list(
      q = linkage_spec(~ ar1(1 | Year), by = ~ fleet, fleet = 1,
                       observe = "QcovPol", obs_sd = 0.1)))
    ```

  The AR1 process, the covariate observation, and the effect size are all
  carried by the one term; `observe` requires an `ar1` term and a positive
  fixed `obs_sd`.

* **Priors on Dirichlet-multinomial data weights (`build_composition()`).** The
  DM likelihood self-tunes each composition dataset's effective sample size, but
  that weight can be poorly identified when a fleet has few comp years. A new
  `build_composition()` linkage attaches a prior to the DM weight, keeping the
  implied effective sample size in a believable range without reverting to
  Francis / McAllister–Ianelli hand-tuning:

    ```r
    fit_mod(...,
      compFun = build_composition(linkages = list(
                  theta_comp = linkage_spec(~ 1, by = ~ fleet, fleet = c(1, 2),
                                 priors = list(`(Intercept)` = gamma(2, 0.5))))))
    ```

  Keys `theta_comp` / `theta_caal` / `theta_diet` cover age/length composition,
  conditional age-at-length, and predator diet. The process is prior-only (the
  weight is a scalar, not year-varying): a covariate, an `init`, or an
  `est_phase` on a composition spec, or a linkage on a fleet/predator not fit
  with the `DirichletMultinomial` likelihood, errors up front.

* **Formula-linkage effect sizes are now reported.** The estimated linkage
  coefficients (`beta_linkage`), the random-effect deviations (`beta_linkage_re`),
  and the Rogers-2024 QAR1 effect size (`beta_linkage_obs`) appear in
  `fit$quantities`, and the coefficients / effect size are `ADREPORT`'d so they
  carry a standard error in the `sdreport`. Previously the QAR1 effect size — the
  quantity the model exists to estimate — was buried in the raw parameter vector
  with no readable exposure or uncertainty.

* **QAR1 covariates need not span the whole series.** A state-space `observe`
  covariate (e.g. a survey index that began mid-series) is now handled without
  hand-padding: `env_data` is auto-extended to start at `styr` with `NA` for
  missing years, the latent AR1 spans all hindcast years, and the observation is
  applied **only** where the covariate is present (un-observed years are masked
  out — no fabricated observations, unlike the legacy mean-fill). Fixed-effect
  covariates must still be finite over the model range: a missing year in a
  fixed-effect (non-`observe`) covariate now errors clearly rather than silently
  producing an NaN.

* **QAR1 observation SD can be fixed or estimated.** The Rogers-2024 `observe` /
  `obs_sd` state-space covariate holds the measurement SD **fixed** by default;
  pass `linkage_spec(..., obs_sd_est = TRUE)` to **estimate** it (one per observed
  group, started from `obs_sd`), as the reference `Estimate_q = 6` / GOApollock
  model does. Estimation is opt-in because a freely-estimated `obs_sd` collapses
  toward 0 on a smooth covariate (the effect-size / `obs_sd` identifiability
  degeneracy) — keep it fixed unless the covariate is informative (a prior on
  `obs_sd` to regularise this is future work).

* **Selectivity penalties can be given as standard deviations.** The cryptic
  penalty *weights* in `Sel_curve_pen1/2/3` can instead be supplied as SDs via
  `Sel_shape_sd` (+ `Sel_shape_dir` `"Decreasing"`/`"Increasing"` for the
  directional sign), `Sel_curvature_sd`, and `Sel_devmag_sd`. Each penalty is a
  Gaussian SSQ, so `switch_check()` converts `weight = 1/(2*sd^2)`. Legacy
  `Sel_curve_pen` values are never overwritten (existing models are bit-identical);
  a fleet supplying the SD columns fits equivalently. `Sel_shape_sd` /
  `Sel_devmag_sd` apply to `NonParametric` (2/9) and `LogisticPM` (11);
  `Sel_curvature_sd` is `NonParametric`-only (LogisticPM does not use
  `Sel_curve_pen2`). Setting an SD column on a form that doesn't use that slot as a
  weight (e.g. 2D/3D-AR1, which reuse `Sel_curve_pen` for logit-scale correlations),
  or a non-positive SD, errors clearly.

* **Readable string aliases for integer-coded switches.** Following the CIE
  review's "rename options to improve interpretability" (e.g. *constant → not
  estimated*), integer-only switches now also accept self-explanatory strings,
  resolved by `switch_check()` (or `build_srr()`) to the same integer codes (so
  fits are identical). New aliases follow a consistent convention — the
  not-estimated value is `"Fixed"`, estimated is `"Estimated"`:
    - `estDynamics`: `"Estimated"` (0) / `"Fixed"` (1) / `"FixedScaled"` (2) /
      `"FixedScaledByAge"` (3);
    - `Estimate_index_sd` / `Estimate_catch_sd`: `"Fixed"` (0) / `"Estimated"` (1) /
      `"Analytical"` (2);
    - `srr_est_mode` (via `build_srr()`): `"Fixed"` (0) / `"Estimated"` (1) /
      `"LognormalPrior"` (2) / `"BetaPrior"` (3);
    - `suitMode`: the documented string map (`"Empirical"`, `"GammaWeight"`, …) is
      now actually applied (previously defined but never wired).

## Bug fixes

* **`M1_re = 6` (separable age × year 2D-AR1 on M) now estimates its correlations.**
  A gate bug in `build_map()` (`if(M1_re_model == 5)` inside a block reachable only
  when the mode is 3 or 6) plus a reference to an undefined index left mode 6's age
  and year AR1 correlations unmapped, so the separable AR1 silently collapsed to IID
  (identical to `M1_re = 3`). Mode 6 now frees both `rho` hyperparameters, so it is a
  genuine separable 2D-AR1. This changes results only for models using `M1_re = 6`
  (none of the bundled examples do); golden references are unaffected.

* **The deprecated `est_M1` name is now recognized everywhere `M1_model` is.**
  The natural-mortality estimation switch was renamed `est_M1` → `M1_model`, but
  the data dictionary still documented `est_M1` and no alias existed, so a data
  list carrying `est_M1` was silently ignored — it fell through to the `M1_model`
  default with no warning. `est_M1` is now folded into `M1_model` in `fit_mod()`
  (before the `build_M1()` reconciliation), `switch_check()`, and `combine_data()`,
  with a deprecation message. As with `M1_model`, the value set on the data list
  is a default that an explicit `build_M1(M1_model = ...)` argument overrides — but
  a data-list value that differs from the `build_M1()` setting now *warns* rather
  than being dropped silently. The recommended way to request M1 estimation remains
  `fit_mod(..., M1Fun = build_M1(M1_model = ...))`. The dictionary now documents
  `M1_model`.

* **A catchability linkage on a non-estimated q now errors.** A `q` linkage
  (environmental or random-effect) on a fleet whose `Catchability` is `"Fixed"`
  or an analytical form (`"Analytical"` / `"AnalyticalArith"`) is rejected up
  front, naming the fleet. Previously only the analytical forms were caught, so
  a linkage on a `"Fixed"` fleet silently turned a fixed catchability
  time-varying — contrary to the assessor's setting. Set `Catchability` to
  `"Estimated"` / `"Estimated-with-prior"` to link q.
  This is the same statistical model as the existing `Time_varying_q`/`_sel`
  deviate processes, expressed through the linkage grammar.

## Deprecations

* **`Selectivity_block` is now optional; `Q_block` is deprecated and ignored.** The
  per-observation time-block columns in `index_data`/`catch_data` are read only for
  `Time_varying_{sel,q} = "Block"`; every other configuration (including the
  random-effect linkages) ignores them. `clean_data()` now default-fills a missing
  `Selectivity_block` with `1` (a single block), so you need only supply it for
  Block-mode fleets. `Q_block` was never read (q time-blocking reuses
  `Selectivity_block`) — supplying it now warns, and it has been dropped from all
  bundled example datasets (as has the unused `Selectivity_block`). Existing models
  are bit-identical.

* **`Time_varying_q`, `Time_varying_sel`, and `M1_re` are soft-deprecated in favor
  of random-effect linkages.** `fit_mod()` now warns (naming the fleets/species and
  the grammar equivalent) when a model uses these legacy time-variation switches,
  pointing at `build_catchability()` / `build_selectivity()` / `build_M1()` with
  `(1 | Year)` / `rw(1 | Year)` / `ar1(1 | Year)` — which additionally allow a prior
  on, or free estimation of, the deviation SD. **The legacy switches still fit with
  their exact numerics** (they keep their own C++ path); this is only a nudge. The
  warning is not raised where no grammar equivalent exists yet: the environmental /
  Rogers-AR1 catchability modes (which overload `Time_varying_q` to name env
  columns), the non-parametric selectivity forms, and the separable `M1_re = 6`
  (age × year).

* **`fleet_control` columns `Q_prior` → `Q_init`, `Index_sd_prior` → `Index_sd`,
  `Catch_sd_prior` → `Catch_sd`.** These are start/input values, not priors (the
  prior on q lives in `Q_sd_prior`), so the names misled. As with the other
  renames, the old names are accepted and upgraded in place by `switch_check()` /
  `read_data()`, and the bundled datasets were regenerated — existing scripts keep
  fitting identically.

* **`fleet_control` columns `Time_varying_q_sd_prior` / `Time_varying_sel_sd_prior`
  renamed to `Time_varying_q_sd` / `Time_varying_sel_sd`.** They hold the input
  *value* of the time-varying deviate SD, not a prior on it (no density is placed
  on the SD), so the `_prior` suffix was misleading. The old names are still
  accepted — `switch_check()` and `read_data()` upgrade them in place with a
  one-time message — so existing data lists, saved xlsx files, and the bundled
  example datasets keep fitting identically. Update scripts to the new names at
  your convenience.

# Rceattle 4.9.0

## New features

* **Environmental linkages on catchability and selectivity.** The
  formula-driven linkage grammar (`linkage_spec()`) now extends to survey
  catchability and to the parametric selectivity forms, both indexed by a new
  `fleet` stratum (`by = ~ fleet`):

    ```r
    fit_mod(...,
      qFun   = build_catchability(linkages = list(
                 q = linkage_spec(~ temp, by = ~ fleet))),
      selFun = build_selectivity(linkages = list(
                 inf_asc = linkage_spec(~ cold_pool, by = ~ fleet))))
    ```

  Selectivity keys are the shared parameter slots `slp_asc`/`slp_desc`/
  `inf_asc`/`inf_desc`, with DoubleNormal aliases `sigma_asc`/`sigma_desc`/
  `peak`/`right_floor`. Supported forms: `Logistic`, `DoubleLogistic`,
  `DescendingLogistic`, `DoubleNormal`, `LogisticPM`. Every parameter accepts
  `link = "log"` (multiplicative) or `link = "identity"` (additive). A linkage
  on an unsupported form or the non-parametric `coff` errors at fit time
  naming the fleet, rather than being estimated to no effect.

* **`build_catchability()` and `build_selectivity()`** exported, mirroring
  `build_growth()` / `build_M1()` / `build_srr()`, each carrying a `linkages`
  argument for its process.

* Linkage formulas are now parsed with the **`reformulas`** package (shared by
  lme4 and glmmTMB), so `(1 | Year)` / `ar1(Year + 0 | fleet)` syntax is
  recognized and an unknown covariance-structure wrapper errors instead of
  silently degrading to unstructured.

## Behavior changes

* **`Catchability = "Environmental"` (`Estimate_q = 5`) is deprecated.** It
  still fits with its existing numerics but emits a warning pointing at
  `build_catchability(linkages = ...)`, which names covariates by formula and
  carries priors, bounds, and an estimation phase.

# Rceattle 4.8.0

## New features

* **AMAK "avgsel" base-level selectivity penalty** (`fleet_control$Sel_avgsel_pen`).
  Non-parametric (type 9 / `NonParametricPM`) fleets can now carry the AMAK/ebswp
  base-level regulariser `weight * (log(mean(exp(base coefficients))))^2` — ADMB's
  `10 * square(avgsel_*)` — which mildly pins the overall level of the base
  selectivity coefficients that the per-year mean-centering otherwise leaves free.
  The per-fleet weight defaults to `0` (off), so existing models and the `BS2017SS`
  golden reference are unchanged; set `Sel_avgsel_pen = 10` to match AMAK. The
  penalty is accumulated once per shared-selectivity block (on the lead fleet).

## Bug fixes

* **Non-parametric (AMAK "pm", type 9 / `NonParametricPM`) selectivity deviate
  penalties.** Two corrections so the type-9 selectivity and its deviate penalty
  reproduce ADMB/AMAK exactly when the fleet excludes a first bin (e.g. the
  acoustic-survey age-1) and/or starts after the model start year:
    - Excluded bins (below `Bin_first_selected`) are now held at 0 before each
      year's mean-centering instead of being carried through the random walk.
      Previously their log-selectivity accumulated the per-year centering offset
      and drifted, inflating the curvature / random-walk penalty on a bin that is
      zeroed out of the fit anyway.
    - For years up to a fleet's `Sel_start_year`, the curve is now rebuilt from the
      base coefficients each year rather than iterating the running random walk
      over the (data-free) pre-survey years. Iterating instead converged the
      excluded first bin to a different fixed point and perturbed the start-year
      base selectivity, inflating the deviate penalty on the first change year.
  A fleet that starts at `styr` with no excluded bins is unaffected (the reset
  reduces to the original start-year behaviour); the `BS2017SS` golden reference is
  unchanged.
* **`hessian_conditioning` diagnostic now always names the flat direction.** When
  the Hessian's least-identified direction was spread diffusely over many
  coefficients (no single coefficient above the reporting threshold), the message
  read `loads on: .` with nothing after it. The check now aggregates the
  eigenvector's squared loadings by parameter block and reports the block(s) making
  up the direction with their percentage share (e.g. `loads on: rec_dev (69%) +
  ln_srv_sel (31%)`), falling back to `par.fixed` names when `cov.fixed` carries no
  dimnames.

* **Mohn's rho was computed from the wrong row for forecast years beyond the
  terminal year.** `retrospective()` accumulated relative error into
  `mohns[ind, ]` while reading from `mohns[j, ]`; the two indices coincide only
  at forecast year 0, so retrospective bias was wrong for every additional
  forecast year. The observation counter on the adjacent line was correctly
  indexed, which masked the problem.

* **Projection-year conditional age-at-length was missing from every MSE.**
  `run_mse()` built `proj_caal` and then appended a different object
  (`proj_comp`) to a non-existent `proj_caal` field instead of appending the
  CAAL rows to `caal_data`. When the composition branch had not run,
  `proj_comp` did not exist at all.

* **`sim_mod()` reproduced the previous observation's draw instead of
  simulating.** The CAAL branch tested `Comp_loglike` (rather than
  `CAAL_loglike`) against integer codes, and neither the composition nor the
  CAAL branch handled `"MultinomialAFSC"` — the default `Comp_loglike`. In both
  cases no branch fired, leaving the simulation variable holding its value from
  the previous row. Unrecognised likelihood codes now raise an error rather
  than falling through.

* **`sample_rec()` and retrospective forecasts produced `NaN` recruitment
  deviations.** Both took `log()` of the mean of `log(R) - log(R_hat)`, a
  log-scale quantity centred near zero whose mean is routinely negative. The
  trend adjustment is now applied additively on the log scale, matching the
  sibling sampling branch.

* **`combine_data()` discarded merged environmental data** (it assigned the
  merge to its input rather than to the returned object) and appended two
  `NA`-named elements by iterating past the end of its column-name vector.

* **`read_data()` overwrote a user-supplied `fleet_control$Month`** with zero.
  The guard was inverted relative to its own message, so per-fleet survey
  months set in the workbook were silently reset.

* **`build_hcr_map()` disabled SPR reference points too eagerly.** It turned off
  `Ftarget`/`Flimit` when *any* fishery for a species had zero `proj_F_prop`,
  though its comment and message both describe the all-zero case.

* **`set_phases()` declared `log_M1` twice**, misaligning the positional
  phase-to-parameter pairing in `TMBphase()` for every subsequent parameter.

* **`fit_mod()` could fail with `object 'opt' not found`** for
  `estimateMode = 2` combined with `HCR = "NoFishing"`.

* **`Time_varying_q` environmental-index lists were destroyed.** For fleets with
  `Catchability = "Environmental"`, a comma-separated list of `env_data` column
  indices (e.g. `"1,3"`) was coerced with `as.integer()` and became `NA`, both in
  `fleet_control` itself and in the `index_varying_q` vector passed to TMB.
  `process_residuals()` tests that vector against a mode code, so the `NA`
  silently poisoned a branch rather than failing.

* **Natural-mortality environmental linkages used the wrong coefficients.**
  Assigning a scalar to the vector `beta_M1_tmp` broadcast it, so the dot
  product collapsed to `beta_last * sum(env)` rather than `sum_j env_j * beta_j`.
  Dormant in practice because `M1_beta` is mapped out by default.

* **Observation years earlier than `styr` produced a large negative array
  index.** Two consecutive `if` statements (rather than `if`/`else if`) applied
  both year transformations in sequence. Five sites.

* **Selectivity max-normalisation branched on an AD type**, freezing the
  location of the maximum at its initial-parameter position for the whole
  optimisation. If the peak moved during fitting, the curve was normalised by
  the wrong bin and the gradient with respect to the true peak was absent. Now
  folded with `max2()`, matching the neighbouring single-bin normalisation.

* **`plot_form()` errored on every call.** It plots the Kinzey & Punt (2009)
  functional responses, whose parameters (`logH_1`, `logH_1a`, `logH_1b`,
  `logH_2`, `logH_3`, `H_4`) are commented out in the TMB template,
  `build_params()`, `build_map()` and `build_bounds()`, and whose modes
  (`msmMode` 3-9) `data_check()` blocks. `params$logH_1` was therefore always
  `NULL` and the first line failed with "non-numeric argument to mathematical
  function". The function remains exported but now raises a clear message
  naming the supported alternatives; its body is retained, commented, so it can
  be revived with the C++ if these formulations are validated.

* **Predation suitability could become infinite.** The zero-guard tested
  `diet_prop_sum`, but the division is by `suma_suit + other_food_diet_prop`,
  which can reach zero or go negative independently.

## Behavior changes

* **`estimateMode = 3` now returns the real objective and gradient.** It
  previously shared mode 4's placeholder objective (`dummy^2`), so `obj$fn()`
  returned zero and every gradient element was identically zero — a build-only
  object appeared healthy while carrying no information. Mode 3 is the
  build-without-optimizing entry point (the analogue of WHAM's
  `fit_wham(do.fit = FALSE)` and SAM's `sam.fit(run = FALSE)`) and can now be
  used to inspect a model before committing to a fit. **`estimateMode = 4` is
  unchanged** and still returns the placeholder: it maps out every hindcast
  parameter, so it is a plumbing smoke test rather than a likelihood.

## Internal

* **Observed-proportion columns are now found by name rather than by counting
  from a fixed position.** `comp_data` and `caal_data` begin with identifying
  columns (fleet, species, sex, year, sample size) followed by the observed
  proportion at each age or length bin. `sim_mod()` located those proportions
  by a fixed offset (`[, 9:ncol(x)]`, `[, 7:ncol(x)]`), so adding or reordering
  an identifying column would have written simulated values into the wrong
  columns without any error. Both tables now resolve the proportion columns by
  name, as `rearrange_data()` already did.

* **Added an internal parameter dictionary** (`R/0-parameter_dictionary.R`)
  mapping each of the 43 TMB parameters to its natural-scale name, its process
  group, a one-sentence meaning, and its dimensions — so diagnostics and error
  messages can describe `log_sel_slp_dev` or `index_q_rho` in terms a reader
  can act on. A test enforces that the dictionary and the parameter list stay
  in exact correspondence.

* The non-parametric selectivity curvature penalty is computed in one pass
  rather than recomputing the full second difference for each age (O(n) instead
  of O(n^2)). Numerically identical.

# Rceattle 4.7.0

## New features

* **Natural-scale normal survey-index likelihood.** `fleet_control$Index_loglike`
  gains `"Normal"`: the index residual `obs - q*pred` is normal with an *absolute*
  observation SD (the `index_data$Log_sd` column is read as a natural-scale SD, not
  a log-scale CV), i.e. `0.5*(obs - q*pred)^2 / sd^2`. This matches the AMAK/ebswp
  `avo_like` / `cpue_like` survey likelihoods (which use an absolute sigma) so those
  indices can be reproduced exactly, rather than approximated by the default
  lognormal. Fixed alongside a latent crash: the covariance (MVN/MVNORM) block was
  gated on `Index_loglike >= 1`, which now also matched `"Normal"` and applied
  `MVNORM()` to a 1x1 dummy Sigma; it is now restricted to `"MVN"` / `"MVNORM"`.
* **Multivariate-normal (covariance) survey-index likelihood.** A survey/index
  fleet can now use a full variance-covariance matrix for its biomass index
  instead of independent lognormal errors, via the new `fleet_control` column
  `Index_loglike` (`"Lognormal"`, the default, or `"MVN"` / `"MVNORM"`). Supply
  the covariance matrix (e.g. a VAST-derived Sigma) as a named element of the new
  `index_cov` data list, keyed by `Fleet_name`. The likelihood uses TMB's native
  `density::MVNORM(Sigma)` on the natural-scale residual `r = obs - q*pred`:
  `"MVN"` reports the bare quadratic form `0.5 * r' Sigma^-1 r` (the AMAK/ebswp
  `DoCovBTS` bottom-trawl survey value, matching ADMB's reported likelihood),
  while `"MVNORM"` reports the full normalized density
  `0.5 * (r' Sigma^-1 r + logdet(Sigma) + n*log(2*pi))` — the two give an
  identical fit and differ only by a fixed constant. A companion catchability
  option `Catchability = "AnalyticalArith"` gives the arithmetic-mean analytical q
  (`mean(obs)/mean(pred)`) that the AMAK covariance survey uses (the existing
  `"Analytical"` q remains the geometric mean). Defaults are fully back-compatible:
  existing models get `Index_loglike = "Lognormal"` and are numerically unchanged.
* **AMAK-style non-parametric and logistic selectivity forms.** Two new
  `fleet_control$Selectivity` options reproduce the ADMB AMAK ("pm") selectivity:
  `NonParametricPM` (Ianelli coefficient selectivity with the decreasing,
  curvature, and deviation-magnitude penalties, whose weights are set by
  `Sel_curve_pen1` / `Sel_curve_pen2` / `Sel_curve_pen3`) and `LogisticPM`
  (logistic + a free age-1 selectivity). Both take `Time_varying_sel =
  "RandomWalk"`. Four new per-fleet `fleet_control` columns tune the shape
  penalty: `Sel_pen_first_bin` and `Sel_pen_last_bin` (the bin range of the
  adjacent-bin shape penalty, letting it span a narrower range than the first
  selected bin), `Sel_shape_mode` (`"Directional"` one-sided decreasing/increasing,
  the AMAK default, or `"Smooth"` two-sided curvature), and `Sel_cap_bin`
  (hold the realized non-parametric curve flat at/after a bin). Each is given on
  the fleet's own selectivity dimension — an age for age-based fleets, a 1-based
  length-bin ordinal for length-based fleets — and is range-checked accordingly.
  All default to the previous behavior when unset.
* The plotting functions have been overhauled to **ggplot2**. Every exported
  `plot_*()` function now builds its figure with ggplot2 (colourblind-safe
  palettes — the Okabe-Ito qualitative palette for series identity and viridis
  for ordered magnitude such as year — and `theme_classic`) and **returns the
  `ggplot` object**, so
  plots print when called interactively and can be further customised
  (`plot_biomass(fit) + ggplot2::ggtitle(...)`) or saved with `ggplot2::ggsave()`.
  Pass `file = "stem"` to also write the figure to `stem_<plot>.png`. Plot
  conventions follow r4ss / WHAM / SAM (line + 95% ribbon time series; observed
  points + error bars + predicted line for index/catch fits; year-coloured
  at-age curves for selectivity and mortality surfaces).
* `plot_index()` gains a `log` argument for the log-scale survey-index fit.
* The time-series plotters (`plot_biomass()`, `plot_ssb()`, `plot_recruitment()`,
  and the other `plot_timeseries()` wrappers) again honour user-supplied
  `line_col`, `lwd`, and `lty`: pass a colour, line width, and/or line type per
  model to override the defaults
  (e.g. `plot_biomass(list(m1, m2), line_col = c("black", "red"), lty = c(1, 2))`).
  `lwd` keeps the base-graphics convention where the default (3) renders as a
  standard-weight line.
* `plot_stock_recruit()` adds a 95% data ellipse of the SSB–recruitment cloud
  (`add_ci`, default `TRUE`).
* The test suite is reorganised into a flat, navigable `tests/testthat/` (see
  the new `tests/testthat/README.md`), runs fully under continuous integration,
  and has automated code-coverage reporting.

## Breaking changes

* `plot_*()` functions now **return a `ggplot` object** instead of drawing as a
  base-graphics side effect (returning `NULL`). Scripts that only called them
  for the side effect still work (the object prints); scripts that depended on
  the base-graphics device state or on a `NULL` return may need updating.
* `plot_logindex()` has been **removed**; use `plot_index(..., log = TRUE)`.
* The `gplots` and `oce` dependencies have been dropped (no longer used).
* **Time-varying non-parametric selectivity now uses `"RandomWalk"`, not
  `"IID"`.** The Ianelli non-parametric form (`Selectivity = "NonParametric"`)
  previously fired its time-varying coefficient deviations on `Time_varying_sel =
  "IID"`; it now requires `"RandomWalk"` (matching the random-walk structure the
  penalty implements) and rejects `"IID"` with an error at `build_map()`. A model
  using `NonParametric` + `IID` must switch to `RandomWalk`.
* **A selectivity shared across fleets is now penalized once.** When two or more
  fleets share a `Selectivity_index` *and* selectivity type (a mirrored curve),
  the selectivity shape/deviation penalty is accumulated only on the lead fleet
  rather than once per fleet, matching ADMB. This changes the objective only for
  models that both mirror a selectivity **and** put a penalty on it (non-parametric
  or time-varying); models with a unique selectivity index per fleet, and all
  bundled examples, are numerically unaffected.

## Bug fixes

* **`Catchability = "AnalyticalArith"` left an unused free catchability parameter,
  making the Hessian singular.** The arithmetic-mean analytical q solves q from the
  data (like the geometric `"Analytical"`), so its `index_log_q` is never used —
  but `build_map()` excluded only `"Analytical"` from estimation, so the
  `AnalyticalArith` fleet's `index_log_q` was still freed. That parameter never
  entered the objective, leaving a zero-gradient flat direction that prevented
  `sdreport()` from inverting the Hessian (`pdHess = FALSE`). It is now mapped
  out, and such models converge with an invertible Hessian.
* **The catchability prior and deviate penalties were counted once per fleet
  sharing a `Q_index`, not once per estimated parameter.** Fleets sharing a
  `Q_index` estimate one catchability and one deviate vector, but the model
  looped over every fleet, so a mirrored pair applied the `Q_prior` twice to the
  same parameter — tightening an intended prior SD of 0.2 to 0.2/sqrt(2) — and
  penalized the shared `index_q_dev` vector once per fleet for the IID, random
  walk and AR1 forms. A new `flt_q_lead` (the catchability analogue of
  `flt_sel_lead`) accumulates them on one fleet per group. Models whose fleets
  all have distinct `Q_index` are unchanged.
* **`Catchability = "PowerEquation"` is now rejected as not yet implemented.** It
  was accepted as a valid switch, but the power coefficient (`index_q_pow`) is
  not built as a parameter and the model does not apply it, so the fleet
  silently got a plain estimated q. `data_check()` now errors, matching how
  length-based `suitMode` values are handled.
* **`flt_sel_lead` could put the selectivity penalty on an `Off` fleet.** The
  lead was the first fleet in a `Selectivity_index` group by row order. When that
  fleet was `Fleet_type = "Off"` the model's `flt_type > 0` gate then skipped
  the penalty for the whole group, leaving the shared selectivity unpenalized.
  The lead is now the first *estimated* fleet in the group, matching the map
  donor.
* **A mirrored group led by an `Off` fleet stopped estimating selectivity and
  catchability.** `adjust_map_shared_params()` copied the first sharing fleet's
  map slice onto the rest. When that fleet was `Fleet_type = "Off"` its slice is
  all `NA`, so every fleet sharing the index silently had its selectivity /
  catchability parameters fixed at their starting values. The donor is now the
  first *estimated* fleet in the group, and the copy is skipped when the group
  has none. Groups led by an estimated fleet are unchanged.
* **An explicitly set `Sel_start_year` was not shared across mirrored fleets,
  making the fit depend on `fleet_control` row order.** The default derived from
  the data is already the earliest first-observation year across each
  `Selectivity_index` group, but a value the user sets directly was used
  per fleet. Since fleets sharing an index share one deviation block,
  `adjust_map_shared_params()` then overwrote the mirrored fleet's mask with the
  lead fleet's, so whichever fleet appeared first governed the group: when that
  fleet started later, a sharing fleet with earlier data silently lost those
  deviations (12 years in a 1982/1994 pair). `Sel_start_year` now resolves to the
  group minimum for both the map mask and the model's penalty anchor, however
  it was set. `data_check()` warns when a mirrored group has differing
  `Sel_start_year`, and `build_map()` warns when `Bin_first_selected` or
  `N_sel_bins` differ within a group (those are likewise taken from the lead
  fleet). Unmirrored fleets and derived defaults are unchanged.
* **`fleet_control$Fleet_code` is now required to equal the row number.** It is
  used directly as the fleet slot of the per-fleet parameter and map arrays,
  which are built in `fleet_control` row order; a mismatch silently attached
  parameters to the wrong fleet. `data_check()` now rejects it, and the
  remaining places that read a `fleet_control` column by `Fleet_code` instead of
  row index were corrected.
* **Selectivity bin columns were converted to model indices using the species'
  `minage`, which is wrong for length-based fleets.** `Sel_norm_bin1`,
  `Sel_norm_bin2`, and the shape-penalty range/cap columns subtracted
  `minage` to get the 0-based template index. That is correct for an age-based
  fleet (the value is an age), but a length-based fleet's value is a 1-based
  length-bin ordinal and must be offset by 1 — so those columns silently pointed
  at the wrong length bin whenever `minage != 1`. The offset is now chosen per
  fleet from `Selectivity_dimension`. Age-based fits, and length-based fits with
  `minage == 1` (where the two offsets coincide), are unchanged.
* **Non-parametric selectivity penalties now span length bins for length-based
  fleets.** The shape, curvature, and random-walk penalties for the
  `NonParametric` (type 2), `NonParametricPM` (type 9), and `LogisticPM`
  (type 11) forms were hard-coded to the number of ages, so a length-based fleet
  (`Selectivity_dimension = "Length"`) with more length bins than ages left the
  upper length bins unpenalized (and read the wrong array). The penalties now run
  over the fleet's own selectivity dimension (`nlengths` for length-based,
  `nages` for age-based), so every parametric and non-parametric selectivity form
  works on both age and length. Age-based fits are numerically unchanged.
* **Survey-index covariance matrices were not re-aligned when the fitted year
  range changed.** An `index_cov` (MVN/MVNORM) Sigma is positionally keyed to a
  fleet's fitted survey observations, so any workflow that changes that set —
  a `retrospective()` peel, an `endyr` / `styr` subset, or a `run_mse()`
  assessment step that appends survey observations — left the Sigma at its
  original dimension and tripped `rearrange_data()`'s dimension check
  (`"N x N but the fleet has M fitted survey observations"`). `clean_data()` now
  tags each Sigma with its fitted years the first time it is seen and, on every
  subsequent pass, re-keys it to the current fitted set: retained years keep
  their full covariance block, and new (future/simulated) years are added as an
  independent diagonal block with variance `(Observation * Log_sd)^2`. Because
  every re-fit routes through `clean_data()`, retrospective, MSE, and jitter now
  all work with covariance-survey models; fresh fits and non-MVN fleets are
  numerically unchanged.
* **Time-varying selectivity deviations were estimated before a fleet had any
  data.** `build_map()` never consulted `fleet_control$Sel_start_year`, so a fleet
  with time-varying (`"RandomWalk"`) selectivity had deviations estimated across
  *every* hindcast year — including years before its first observation. Those
  deviations are informed by no data and constrained by no penalty (every
  selectivity penalty in the objective is anchored at `Sel_start_year`), leaving
  unidentified flat directions that stall the optimizer. In the EBS pollock model
  this left ~54 free pre-survey deviations on a bottom-trawl survey starting in
  1982 and ~240 on an acoustic survey starting in 1994 — a total parameter count
  of 1483 against 1225 for the equivalent ADMB (AMAK) model, which never declares
  them. Deviations before `Sel_start_year` are now fixed at 0, giving 1249. The
  two selectivity parameterizations differ in where the base curve lives and are
  handled accordingly: `LogisticPM` (and other curve-based forms) estimate a
  separate base, so deviations are fixed *through* the start year; the
  non-parametric random walk maps its mean (`sel_coff`) off and lets the start-year
  deviation carry the base, so only the deviations strictly *before* it are fixed.
* **`Sel_start_year` now defaults to the fleet's first year of data** rather than
  `styr`. It is an optional `fleet_control` column, so the fix above only took
  effect for users who knew to set it — a model with a late-starting survey would
  silently carry unidentified deviations. The default is derived from
  `catch_data` / `index_data` / `comp_data` / `caal_data`, consistent with how
  `switch_check()` already auto-`"Off"`s fleets with no observations. Fleets
  sharing a `Selectivity_index` share one selectivity curve, so the start year is
  the *earliest* first-observation year across the whole group: a fleet whose own
  data starts late but which mirrors an earlier fleet's curve (e.g. an AVO index
  starting in 2006 mirroring an acoustic survey starting in 1994) must not drop
  the deviations the mirrored fleet's data informs. Set the column explicitly to
  override. Only models with time-varying selectivity on a fleet whose data starts
  after `styr` are affected; for those, the previous behaviour is recovered by
  setting `Sel_start_year = styr`.
* **`LogisticPM` selectivity started with an unusable age-1 selectivity.**
  `build_params()` initializes `sel_inf[2]` to `10`, which is correct for its
  usual meaning (the descending-limb inflection *age*), but `LogisticPM`
  (type 11) repurposes that slot as the free first-bin (age-1) **log**-selectivity.
  The shared default therefore started age-1 selectivity at `exp(10)` = 22026
  (for reference, the ADMB AMAK "pm" estimate is `sel_age_one_bts` = -3.19, i.e.
  0.04). When age-1 was selected (`Bin_first_selected = 1`) this swamped the
  survey-index prediction and the optimizer diverged with a gradient blow-up on
  `sel_inf`; when it was not selected the bad value was silently masked by the
  zeroed first bin. `sel_inf[2]` now defaults to `0` (age-1 selectivity = 1) for
  `LogisticPM` fleets.
* `plot_maturity()` read a non-existent `pmature` field and errored on real
  fits; it now reads `data_list$maturity`.
* `plot_ration()` failed for single-sex models (a dropped array dimension); the
  sex dimension is now indexed explicitly.
* `plot_stock_recruit()` drew the mean-recruitment reference line a factor of
  1e6 too high (recruitment points are in millions but the line was not scaled),
  which under ggplot's free y-scale collapsed the SSB–recruitment cloud to the
  axis; the line is now scaled to match the points.
* `plot_index()`, `plot_catch()`, and `plot_indexresidual()` drew
  prediction-only rows (`Year < 0`, excluded from the likelihood) as if they
  were fitted observations; these rows are now omitted.
* `plot_mortality(M2 = TRUE)` labelled the y-axis "M1 + M2" while plotting M2;
  the label now reads "M2".
* `osa_residuals()` now includes fleet in the one-step-ahead conditioning order
  (source, year, **fleet**, then bin), making the sequence fully deterministic
  when several fleets report in the same year. Previously the within-year,
  multi-fleet order fell to the incidental row order of the observation table.
  This does not change results for fixed-effects fits (`random_rec = FALSE`),
  where OSA residuals are order-invariant, and ascending fleet order matches
  WHAM's within-stage sequencing so the WHAM cross-check is unaffected; for
  random-effects fits it makes individual residuals reproducible (the overall
  N(0,1) properties were already order-invariant).

# Rceattle 4.6.0

## New features

* `osa_residuals()` now builds the composition / conditional-age-at-length /
  diet one-step-ahead observation data on demand, so OSA residuals can be
  computed from any fit. Previously the composition OSA data had to be assembled
  during the fit via `fit_control(osa = TRUE)`; that pre-build is no longer
  required (the reconstruction reads only the `*_ctl` / `*_obs` arrays already
  carried by the fitted object, and is identical to an `osa = TRUE` fit).
  Aggregate index/catch OSA residuals are unaffected.
* **Breaking:** the `osa` argument to `fit_control()` has been removed; it is no
  longer needed (see above), and passing it now raises an "unused argument"
  error. The `$osa` element is no longer stored on fitted `Rceattle` objects.
* `fit_control()` gains `bias_adjust_obs` and `bias_adjust_proc` (both default
  `TRUE`) to toggle the lognormal bias correction (`-sigma^2/2`) applied to the
  observation (index / catch) and process (recruitment, initial abundance, and
  the Beverton-Holt steepness prior) likelihoods, respectively. The defaults give
  the standard mean-unbiased behaviour (`E[R] = R0`); set a flag to `FALSE` to
  drop the correction for that likelihood group. The values are passed to the TMB
  template as `DATA_SCALAR`s.

* Added one-step-ahead (OSA) residuals for model validation via the new
  exported `osa_residuals()` (Thygesen et al. 2017; Trijoulet et al. 2023).
  Unlike Pearson residuals, OSA residuals are iid standard normal under a
  correctly specified model even for correlated composition data. All fitted
  observation types are supported: survey indices, fishery catch, age/length
  composition, conditional age-at-length, and predator diet (stomach-content)
  composition. Composition data are residualized with the conditional binomial
  / beta-binomial decomposition (multinomial / Dirichlet-multinomial), and the
  `MultinomialAFSC` likelihood is residualized under the full multinomial. The
  fitted objective is unchanged: a new `obsvec`/`keep`/`osa_mode` machinery
  feeds `TMB::oneStepPredict()` while leaving normal fitting bit-for-bit
  identical. The `oneStepPredict()` call is split by observation type so the
  continuous (lognormal) index/catch series and the composition series can be
  residualized with the correct settings; `discrete = TRUE` treats composition
  as discrete (randomized quantile residuals; Dunn and Smyth 1996) while the
  aggregate series stay continuous.
* `fit_control()` gains `comp_offset`, the small proportion offset added to the
  observed and predicted age/length composition, conditional-age-at-length, and
  predator diet (stomach-content) bins before the multinomial /
  Dirichlet-multinomial likelihood (and to the matching OSA observation vector,
  so fitting and OSA residuals use the same offset). It applies to every
  composition likelihood, including the "Martin's" (`comp_ll_type = -1`) form. It
  defaults to `1e-5` (the historical CEATTLE value, which avoids `log(0)` for
  empty bins); set `comp_offset = 0` for a standard WHAM-style multinomial. The
  value is stored on `data_list$comp_offset` (filled by `switch_check()`), so
  internal re-fits (projections, `retrospective()`, `jitter()`, `run_mse()`)
  inherit it; `fit_control(comp_offset = ...)` overrides the stored value.
* `osa_residuals()` gains a `parallel` argument (default `TRUE`) that computes
  the per-observation one-step-ahead loop with `parallel::mclapply()` -- a
  near-linear speedup for models with random effects, where each observation
  triggers a Laplace re-evaluation (set `options(mc.cores = )` to choose cores;
  serial fallback on Windows). The discrete randomized-quantile path stays serial
  so it remains reproducible given `seed`.
* Added `osa_diagnostics()` for the Stewart and Monnahan (2025) statistical
  diagnostics -- the standard deviation of the normalized residuals (SDNR) and
  the lower/upper tail statistics, each with its standard-normal null interval
  -- and a `plot()` method for `rceattle_osa` objects, styled after the
  NOAA-AFSC `afscOSA` package. The plot draws up to two figures: an *aggregate*
  Q-Q figure for the index/catch series (which have no age/length bin, so no
  bubbles), and a *composition* figure pairing the Q-Q plot with signed
  OSA-residual and Pearson-residual bubble plots, with age-based bins (age
  composition and conditional age-at-length) in the left column and length-based
  bins in the right column, each on its own bin axis. Panel headers use the
  fleet name from `fleet_control`. The `plot()` method takes `source` and
  `species` arguments to subset the data shown (mirroring `residuals()`), and
  `combine = FALSE` to draw the age and length composition as separate figures.
  `osa_residuals()` now also carries a `fleet_name` column and (for composition)
  the matching Pearson residuals.
* Added `process_residuals()` for SAM-style process residuals on the model's
  random-effect deviations (recruitment, initial abundance, and catchability),
  validating the process model as a complement to the observation residuals.
* `residuals()` gains `type = "osa"` and `type = "process"`; `plot_comp()` and
  `plot_indexresidual()` gain a `residual_type = "osa"` option that draws the
  OSA diagnostics through the familiar plotting functions.
* `plot_comp()` was re-implemented in ggplot2 for a consistent look with the OSA
  plots: composition Pearson-residual bubbles plus observed-vs-fitted annual and
  aggregated composition figures. The observed area and fitted line now span only
  the observed bins (they no longer extend past the first / last bin), zero
  observed proportions are retained (only `NA` is dropped), joint-sex (Sex == 3)
  data are drawn on a single age/length axis with females above and males below
  zero, and a `species` argument subsets the species shown.

## Bug fixes

* Fixed the unweighted conditional-age-at-length (CAAL) log-likelihood being
  recorded in the composition slot of `unweighted_jnll_comp` instead of the
  CAAL slot. This affects the reported diagnostic matrix only; the fitted
  objective was unaffected.

## API

* `projection_uncertainty` moved from `fit_mod()` into `fit_control()`,
  consolidating it with the other optimizer / reporting controls. Passing it
  directly to `fit_mod()` still works but emits a deprecation warning and
  forwards the value into `fit_control()` -- the same backward-compatible path
  used by the other former `fit_mod()` control arguments (`phase`, `getsd`,
  ...).
* `residuals()` now follows the `stats::residuals.glm()` convention: `type`
  selects the residual *kind* -- `"response"` (default), `"pearson"`, `"osa"`,
  or `"process"` -- and a `source` argument selects the data source(s)
  (`"index"`, `"catch"`, `"comp"`, `"caal"`, or `"all"`, the default), so by
  default residuals for every applicable source are returned. A `species`
  argument subsets to particular species. Pearson residuals are now available
  for the aggregate index/catch series (standardized by the realized observation
  log-SD) and for predator diet via `source = "diet"` (returned on its own in a
  predator/prey schema); `plot_diet_comp()` now draws its diet residuals from
  this single `residuals()` path. `type` selects the residual kind only; data
  sources are selected with `source`.
* The `source` argument is shared across `residuals()`, `osa_residuals()`, and
  `plot()` (it replaces the earlier `types` argument of `osa_residuals()`),
  accepting `"index"`, `"catch"`, `"comp"`, `"caal"`, `"diet"`, and `"all"`, so
  the three entry points select data sources with one consistent vocabulary.

## Documentation

* Reviewed and revised the user vignettes for accuracy, clarity, and
  consistency, targeting assessment scientists rather than developers. Trimmed
  developer-oriented internals and repetition, and corrected the option tables
  so they match the code: `estimateMode`, `suitMode`, selectivity, catchability,
  harvest-control-rule, `M1_model`, and bioenergetics (`Ceq`) codes now agree
  across `model-options-and-functionality`, `single-vs-multispecies`,
  `projections-and-reference-points`, `data-without-excel`, and
  `model-parameterizations`.
* Clarified that the length-based suitability modes (`suitMode = 1/3/5`) are not
  yet enabled; use a weight-based mode (`2/4/6`) for parametric suitability.
* Standardized data-structure terminology across vignettes: weight-at-age is in
  kilograms, and the diet / stomach-content input is `diet_data` (removing the
  legacy `UobsWtAge` / `UobsAgeWt` / `stom_prop_data` names).
* Corrected the `residuals()` examples in the README to the `type` / `source`
  convention (e.g. `residuals(fit, type = "pearson", source = "comp")`).
* Added a website-only "Developer guide" article
  (`vignettes/articles/developer-guide.Rmd`) documenting the fit pipeline, the
  switch system, the TMB module layout, and recipes for adding parameters, data
  inputs, switch options, and likelihood components. Consolidates and updates
  the GitHub wiki developer notes.

# Rceattle 4.5.0

## New features

* Added a general post-fit convergence diagnostics framework via the new
  exported `convergence_diagnostics()` function. It runs a battery of checks
  covering the optimizer gradient, Hessian positive-definiteness and
  conditioning, parameters on bounds, phasing, and parameter estimability, and
  returns a structured `"Rceattle_convergence"` object whose `status` reflects
  the worst severity found (`"OK"`, `"NOTE"`, `"WARN"`, or `"FAIL"`).
* `fit_mod()` now runs the convergence battery automatically and attaches the
  result as `fit$convergence`; `convergence_diagnostics()` can also be called
  directly to re-run it on any fit.
* Added a model-diagnostics vignette section and accompanying unit tests for
  the new framework.

## Bug fixes

* Fixed a parameters-on-bounds misalignment in the convergence diagnostics.
* The joint negative log-likelihood (jnll) is now taken from the optimizer
  objective for consistency.

## Internal / R CMD check

* Moved the suppression of the spurious GCC/Eigen `-Warray-bounds` false
  positive from a `-Wno-array-bounds` compiler flag to a source-level
  `#pragma GCC diagnostic ignored "-Warray-bounds"`, so
  `R CMD check --as-cran` no longer warns about non-portable,
  warning-suppressing compilation flags. Real diagnostics still surface.
* Additional C++ cleanup to remove compiler warnings.

# Rceattle 4.4.2

## Code organisation (no change to fitted results)

The pre-fit pipeline files in `R/` were reorganised so they are easier to
navigate. None of these changes alter model output.

* **File prefixes now follow execution order.** `fit_mod()` runs its stages
  as `clean_data()` -> `data_check()` -> `build_params()` -> `build_map()` ->
  `build_bounds()` -> `rearrange_data()` -> fit -> `rename_output()`, so the
  files were renumbered to match (`data_check` is now `1-`, `build_params`
  `2-`, `build_map` `3-`, `build_parameter_bounds` `4-`). A pipeline map was
  added to the top of `fit_mod()`.
* **Switch lifecycle consolidated** into a single `R/0-switches.R`: the
  string<->integer maps (formerly `0-constants.R`) plus `switch_check()`,
  `revert_switches()`, `validate_switches()`, and `convert_switches()`, with a
  header documenting the order in which they run.
* **HCR helpers co-located.** `build_hcr_map()` moved into `R/0-build_hcr.R`
  alongside `build_hcr()` (the separate `2-build_hcr_map.R` was removed).

## Rename / deprecation

* `rearrange_dat()` is renamed **`rearrange_data()`**. The old name still works
  as a deprecated alias (emits a one-time `.Deprecated()` warning) and will be
  removed in a future release.

## Export hygiene

* `check_composition_data()`, `check_caal_data()`, `calc_mcall_ianelli()`, and
  `calc_mcall_ianelli_diet()` are no longer exported; they are internal
  helpers called only from within the package.

## Internal / R CMD check

* Removed `Rceattle:::` self-references in `build_bounds()` (a package should
  not use `:::` for its own objects).
* `profile.Rceattle()` gained `...` for S3 consistency with the
  `stats::profile` generic.

# Rceattle 4.4.1

## Rename

`profile_param()` was renamed `profile()` and turned into an S3 class.

# Rceattle 4.4.0

## New features

* **Double Normal Selectivity (Type 8)**: Added support for the four-parameter Double Normal selectivity curve (Peak, Ascending SD, Descending SD, and Floor). This includes full support for annual deviates on all four parameters and is compatible with the age- and length-based selectivity engines.
* **Growth SD Control**: The linkage system now supports `sd_L1` and `sd_Linf` parameters. This allows users to specify priors, initial values, and bounds for growth SD endpoints using the same formula-driven interface used for mean growth parameters.
* **Natural-scale linkage API**: Renamed the linkage parameter keys from `log_*` to their natural-scale counterparts (`K` / `L1` / `Linf` / `m`; `M1`; `R0` / `alpha` / `beta`). Internal parameters remain on the log scale.
* **Natural-scale priors**: Standardized prior evaluation so that probability densities are applied to parameters on their natural scale by default, unless a lognormal family is explicitly requested.
* **Natural-scale inits**: Standardized init evaluation so initial values for `"(Intercept)"` parameters are applied on their natural scale.
* **Linkage link functions**: Fully implemented dual-path linkage offsets. `link = "log"` (the new default) applies the offset multiplicatively on the natural scale (additive on the log scale); `link = "identity"` applies it additively on the natural scale.
* **Per-species VB anchor age**: `build_growth()` gains a `growth_age_L1` argument (scalar or length-`nspp` vector) for the age at which mean length equals `L1`. Matches SS3's `Growth_Age_for_L1`. Default `NA` inherits `data_list$growth_age_L1` if set, else falls back to `max(0.5, minage[sp])` so `minage = 0` models pick up an SS3-consistent half-year anchor while `minage >= 1` models stay backwards-compatible.
* **Self-test simulation**: New `self_test()` simulates `nsim` datasets from a fitted model and re-fits the model to each simulated dataset, returning the list of refits. Runs in parallel by default (PSOCK cluster, capped at 2 cores under `R CMD check`) with a per-simulation `seed + i` so results are reproducible under both sequential and parallel execution.
* **Likelihood profile**: New `profile_param()` generalises the legacy `profile_rsigma()` example helper to any parameter slot in `Rceattle$estimated_params`. Supports arbitrary N-D cross-profiles via `expand.grid()` over a list of per-cell value grids (e.g. cross-profiling `log_M1` across sex, or `R_log_sd` across multiple species). Natural-scale aliases `"sigmaR"` / `"R_sd"`, `"M1"`, and `"R0"` / `"alpha"` / `"beta"` apply the implicit `log()` transform and (for the `rec_pars` aliases) auto-fill the column, so `slots` only needs the species index. `slots` defaults to species 1 with a warning. Fits run in parallel on the same PSOCK harness as `jitter()` / `retrospective()` / `self_test()`.
* **Parallel `retrospective()` and `jitter()`**: Both diagnostics now run their independent peels / starts on a PSOCK cluster (same approach as `run_mse()`). New `cores` argument on each (default `parallel::detectCores() - 6`, capped at 2 when `_R_CHECK_LIMIT_CORES_` is set); pass `cores = 1` to force sequential execution.
* **Standard errors in `as.data.frame.Rceattle()`**: The tidy long-format frame now carries a `se` column alongside `value` / `lwr` / `upr`, populated from the TMB `sdreport` for any `ADREPORT`'d quantity. Set to `NA` for non-ADREPORT'd quantities and for fits produced with `getsd = FALSE`.

## Bug fixes

* **SRR logic**: Fixed a bug in `build_srr()` where the `Bmsy_lim` penalty was incorrectly disabled for current Ricker implementations due to an index mismatch.
* **Selectivity RW prior scaling**: Corrected the random walk prior scaling in the TMB template to ensure consistent $4 \times$ SD multipliers for both ascending and descending limb slope/SD parameters.
* **`last_par` returned wrong vector**: Fixed `fit_mod()` so the value stored on the returned fit is the optimizer's last parameter vector rather than a stale prior reference, removing the need for the surrounding `try()` guards in downstream callers.
* **Growth at `minage = 0`**: Fixed a segfault in `growth.hpp` when `minage = 0` by guarding `current_age`, `age_L1`, and the cohort-boundary `age_L1_ceil` against the zero-age anchor. Also corrected length-at-age at `minage = 0` so the closed-form anchor at `L1` is honored on both the within-year and cohort recursion paths.
* **Length-bin midpoint for weight-at-age**: `calculate_weight()` now computes each bin's midpoint as `(lengths[ln] + lengths[ln+1]) / 2` (with the plus-group extended by half the final interior width) rather than `lengths[ln] + (lengths[1] - lengths[0]) / 2`. The previous formula assumed uniform bin widths and silently mis-located the midpoint for non-uniform length bins; for uniform bins (the common case) the two are algebraically identical, so existing fits are unchanged.
* **Plus-group SD-at-age pinned to the upper anchor**: both growth builders in `growth.hpp` now pin the oldest age class's `length_sd` to `exp(growth_log_sd(sp, sex, 1))` (SD at `Linf`), matching WHAM (`SDAA` plus group `= SD_len(1)`), instead of the length-based interpolation between SD(`L1`) and SD(`Linf`) used previously. This makes the plus-group convention consistent across von Bertalanffy and Richards growth and across both builders; expect a small numerical change in `length_sd` at the oldest age relative to prior versions.
* **`fit_mod()` bounds ordering**: `fit_mod()` now indexes parameter bounds by name rather than positional order when assembling `L` / `U` for `nlminb`. Previously, when `map$mapFactor` and `bounds$lower` were not in identical order, parameters could be paired with another parameter's bounds, producing silently wrong constraints. `start_par` is now also subset by name with `drop = FALSE`.
* **`mse_summary()` Tier-3 Flimit**: The internal `flimit_tier3_fun()` returned `Flimit` (its argument) instead of the depletion-adjusted `tier3_flimit` it had just computed, so the Tier-3 (HCR = 5) branch of P(F > Flimit) reduced to the base-Flimit check. Now returns the adjusted vector.
* **`mse_summary()` HCR coercion**: `HCR` is now normalized to its integer code before downstream comparisons (`HCR == 5`, etc.). `build_hcr()` accepts either an integer or a string alias (e.g. `"NPFMC"`); `mse_summary()` previously assumed integer form and silently produced wrong status flags when fits carried the string form.
* **`mse_summary()` OM status at assessment years**: P(F > Flimit) and P(SSB < SSBlimit) are now reported for the OM evaluated at the same assessment years as the EM (previously only the EM's perceived status was returned), and the SSB-limit threshold dispatch is consolidated in one helper so the Tier-3 / Category-1 / dynamic-vs-static cases stay aligned across the EM and OM paths.
* **`clean_data()` inactive-fleet handling**: The auto-Off branch no longer nulls out `proj_F_prop` and `Catchability` on fleets it flips to `"Off"`. The downstream TMB code already ignores those columns for off fleets, and clearing them lost information users needed when re-enabling a fleet.
* **Selectivity block indexing**: Renamed the within-loop `biom_yrs` index vector to `block_yrs` in the selectivity and catchability block branches of `build_map_*()`. The previous name was a leftover from the index-data path and shadowed nothing, but read as if it referred to biomass-observation years.
* **`plot_index()` / `plot_catch()` / `plot_logindex()` warnings**: Wrapped the internal `gplots::plotCI()` calls in `suppressWarnings()` so plotting a fit no longer prints the recurring `"In arrows(...): zero-length arrow is of indeterminate angle and so skipped"` noise when CI half-widths are zero.

## Data checks

* **Empirical growth + CAAL**: `data_check()` now errors when `growth_model == 0` (empirical weight-at-age) is combined with non-empty `caal_data` for a given species. The C++ growth matrix is not populated from the age-transition matrix in the empirical branch, so `pred_CAAL` collapses to ~0 and the multinomial NLL becomes uninformative.
* **Selectivity identifiability**: Fleets with estimated `Selectivity` and `Fleet_type != "Off"` now require at least one positive-sample `comp_data` or `caal_data` row. Either provide composition / CAAL data, mark the fleet as `Selectivity = "Fixed"` with `emp_sel`, or set `Fleet_type = "Off"`.
* **Auto-Off inactive fleets**: `clean_data()` automatically flips `Fleet_type` to `"Off"` for fleets that carry no catch or index observations, preventing the optimizer from drifting on unconstrained selectivity / catchability blocks.
* **`minage` guard**: `data_check()` errors when any species has `minage < 0`.

## Documentation
* Added `vignette("environmental-linkages-and-priors")` (and updated `_pkgdown.yml`) to cover the new linkage intercept behavior, link-function semantics, growth SD endpoints, and Double Normal selectivity.
* Updated all cross-references in `build_srr()` / `build_M1()` / `build_growth()` (deprecation warnings, soft-deprecated arg docs, and the `vignette()` pointers in the model-options vignette) from the old `environmental-linkages` slug to the renamed `environmental-linkages-and-priors` vignette so the soft-deprecation warnings now resolve.
* Expanded the Double Normal (selectivity type 8) doxygen / roxygen so the four estimable parameters (peak, ascending SD, descending SD, and logit right-floor) and their TV deviates are documented in one place; `sel_inf(1)` is the right-tail floor (analogous to SS3 P6 / `end_logit`), not a fixed-by-map placeholder.

## Deprecations

* The soft-deprecated `srr_indices` / `M1_indices` arguments and the
  legacy env-driven integer codes (`srr_fun %in% c(1, 3, 5)`,
  `M1_model %in% c(4, 5)`) continue to work in 4.4.0 with a one-time
  warning that points users at the linkage table. **Removal has been
  rescheduled from v4.2.0 to v4.5.0** to extend the migration window;
  see the "Scheduled removal" section under 4.1.0 below for the
  unchanged cleanup checklist.

# Rceattle 4.3.1

## Bug fixes

* Fixed a segfault in `MakeADFun` triggered by Beverton-Holt / Ricker
  fits with a recruitment linkage. Section 6.9 (expected recruitment)
  was indexing `R0`, `alpha`, and `Beta` with the function-scope `yr`
  variable, which the preceding forecast loop left equal to `nyrs` --
  one column past the end of the `(nspp, nyrs)` matrices. The year-0
  block now uses an explicit `first_yr = 0` constant. Mean-recruitment
  fits (`srr_fun = 0`) happened to dodge the segfault because
  `calculate_recruitment()` doesn't dereference `alpha`/`Beta` for
  that case.
* Fixed the recruitment-parameter offset formula at the top of section
  5.6: `alpha` and `Beta` now apply the linkage offset on the log
  scale (`exp(rec_pars + offset)`) to match the documented contract
  and the `log_R0` formula on the same line. Previously the offset
  was added on the linear scale, which silently corrupted BH/Ricker
  recruitment whenever a non-zero `log_alpha` or `log_beta` linkage
  offset was active.
* Added a `R0/alpha/Beta` column-count assertion at the entry of
  section 6.9 so future stale-loop / sizing regressions surface as a
  clean R-level error from TMB rather than an opaque segfault.

## Reparameterised intercept handling for the linkage system

Linkage rows whose `design_col == "(Intercept)"` no longer carry the
parameter level themselves. Instead:

* The base parameter (`rec_pars`, `log_M1`, `log_growth_pars`)
  remains estimable and holds the level. Phasing and the per-process
  map machinery operate on the base parameter as they would without
  any linkages.
* The `(Intercept)` row's `beta_linkage` slot is fixed at `0` and
  mapped NA. It exists in the table for bookkeeping plus as a hook
  for `init` and `priors`.
* `init = list(`(Intercept)` = X)` on the spec **flows to the base
  parameter** -- it sets the base parameter's starting value rather
  than the linkage row.
* Priors attached to an `(Intercept)` row are **re-targeted to the
  base parameter** in the slot 19 contribution; the prior density
  evaluates against `rec_pars(sp, 0)` / `log_M1(sp, sx, ag)` /
  `log_growth_pars(sp, sx, par)` rather than the (zero) linkage
  value.

For slope-only formulas (`~ 0 + temp`) the behaviour is unchanged:
the base parameter is still mapped NA at its `build_params()`
default, and the linkage row carries the year-by-year offset.

### Recruitment offset semantics

Year 0 of the recruitment block no longer bakes the year-0 covariate
contribution into `R0`. `R0` is computed from `rec_pars(sp, 0)`
alone, and the linkage offset multiplies against `R0` each year
(including year 0):

```
R(yr) = R0 * exp(rec_dev(yr) + linkage_offset(yr))
```

This makes the legacy `srr_fun = 1 / 3 / 5` quirk (which double-counted
`Temp[0]`) obsolete. Users migrating from the legacy paths should now
get clean log-linear behaviour without surprise offsets.

### Schema additions

* New `init_supplied` (logical) column on `Rceattle_linkage_table`
  tracks whether the user explicitly supplied an `init` for that row.
  Used by `build_params()` to decide whether to push a base-parameter
  init.
* New `linkage_is_intercept` IVECTOR in the TMB encoding (set from
  `design_col == "(Intercept)"`) used by the slot-19 prior dispatch
  to evaluate intercept priors against the base parameter.

### Tests

`tests-Dynamics/test-linkage-auto-map.R`,
`tests-Dynamics/test-recruitment-linkage.R`, and
`tests-Dynamics/test-growth-linkage-species.R` were updated to assert
the new contract (base parameters estimable, intercept rows fixed at
0, slope-only offsets in the year-by-year tensor).

# Rceattle 4.3.0

## Tidy long-format extraction: `as.data.frame.Rceattle()`

A new S3 method on `as.data.frame()` flattens derived population
quantities into a long data.frame with columns
`year, species, sex, age, quantity, value, lwr, upr` so that custom
plotting and post-processing don't have to walk the nested
`quantities` list or rely on the dimnames decisions in
`rename_output()`. Two shapes are supported and combined into one
frame:

* **Species-by-year** (default `which`): `biomass`, `ssb`, `R`,
  `biomass_depletion`, `ssb_depletion`, `F_spp`. Other species/year
  series (`B0`, `SB0`, `DynamicB0`, `DynamicSB0`, `DynamicSBF`,
  `exploitable_biomass`, `proj_F`, `fT`) are available by name; pass
  `which = "all"` to get every known quantity present on the fit.

* **Species-by-sex-by-age-by-year**: `N_at_age`, `biomass_at_age`,
  `Z_at_age`, `M_at_age`, `M1_at_age`, `M2_at_age`, `F_at_age`,
  `consumption_at_age`, `B_eaten_as_prey`, `NByage0`, `NByageF`. The
  `age` column is biological age (offset by `data_list$minage`), and
  cells padded out to `max(nsex)` / `max(nages)` for species with
  fewer sexes or ages are dropped rather than returned as `NA`.

`lwr` / `upr` are populated from the TMB `sdreport` for any quantity
that was `ADREPORT`'d (currently `biomass`, `ssb`, `R`); other
quantities and fits produced with `getsd = FALSE` get `NA` for the
band. The `ci_level` argument (default `0.95`) controls width.

## Optional data fields, continued (Phases B, C, D)

Continuing the Phase A work from 4.2.0, three more classes of inputs
that were previously required as non-NULL can now be omitted, with
`data_check()` enforcing them only when the model actually needs them:

* **Phase B: bioenergetics scalars.** `Ceq`, `Cindex`, `Pvalue`,
  `fday`, `CA`, `CB`, `Qc`, `Tco`, `Tcm`, `Tcl`, `CK1`, `CK4` may be
  `NULL` in single-species mode. `switch_check()` fills them with
  safe sentinels so TMB's length-`nspp` `DATA_VECTOR` requirements
  are satisfied. When `msmMode > 0` the scalars are required;
  `data_check()` reports which ones are missing or wrong-length in
  a single grouped error.

* **Phase C: `env_data`.** May be `NULL`. `clean_data()` defaults it
  to a Year-only `data.frame(Year = styr:projyr)` with zero indices.
  Existing checks still error when a feature actually needs an
  index (env-dependent catchability, temperature-dependent
  consumption, env linkages, `srr_indices`, `M1_indices`).

* **Phase D: `emp_sel`.** New requirement check: when any fleet has
  `Selectivity = "Fixed"`, `emp_sel` must be supplied. Other fleets
  do not need it.

* **Tests.** A new `tests-Data-processing/test-optional-fields.R`
  file exercises 25 NULL / requirement scenarios across the four
  phases.

# Rceattle 4.2.0

## Optional data fields & data_check cleanup

Several fields in `data_list` that were previously required as non-NULL
data.frames are now truly optional. Users who do not need composition
data, conditional age-at-length, empirical selectivity, fixed
numbers-at-age, ration data, or diet data can omit them entirely;
`clean_data()` default-fills the missing fields with empty data.frames
that carry the metadata columns the downstream code expects, and
`data_check()` enforces the field only under the conditions where the
model actually needs it.

* **Phase A (this release): `comp_data`, `caal_data`, `emp_sel`,
  `NByageFixed`, `ration_data`, `diet_data` may be `NULL`.** Conditional
  requirements are still enforced (`caal_data` when
  `any(growth_model > 0)`; `NByageFixed` when `any(estDynamics > 0)`;
  `diet_data` when `msmMode > 0`).

* **`data_check()` reorganisation.** The validation function has been
  reorganised into eight topical sections (top-level scalars; per-species
  dimensions; biology; fleet control; observation tables; diet &
  predation; environmental data; switches), with shared `has_data()` /
  `fc_num()` helpers and consolidated duplicate guards. New checks were
  added for year-scalar ordering, lognormal SDs, sample sizes,
  probability ranges, observation values, fleet referential integrity,
  selectivity bin bounds, predation cross-checks, duplicate observations,
  and probability-matrix row sums. Several pre-existing dead branches and
  matrix-`$` access bugs were fixed at the same time.

* **`transpose_fleet_control()` removed.** The deprecated long-format
  fleet_control transposer has been removed from `clean_data()`,
  `read_data()`, and the package namespace.

# Rceattle 4.1.0

## Environmental linkages: a unified, formula-driven API

A new long-format **linkage table** lets users express how process
parameters depend on environmental covariates and on stratifying
factors (species, sex, age) through a single formula-driven helper,
`linkage_spec()`. Each row of the table corresponds to exactly one
estimated coefficient. `fit_mod()` pools every spec into a shared
design matrix `X` and a per-row parameter vector
`beta_linkage`; the TMB template iterates the table once and
accumulates per-process offsets on the linear predictor of the
underlying parameter.

* **New constructor: `linkage_spec()`.** Captures
  `(formula, by, species, link, init, bounds, priors, est_phase)`
  for one process parameter. Anything `model.matrix()` understands
  works: `~ 1`, `~ temp + PDO`, `~ poly(temp, 4)`, `~ I(temp^2)`,
  `~ splines::ns(temp, df = 4)`, `~ temp * PDO`, etc.

* **Per-species formulas.** Register multiple specs against the
  same parameter via `linkages = list(log_K = list(spec_a, spec_b))`
  with each spec's optional `species = ...` argument scoping it
  to a subset of stocks. The pooler unions the design columns
  across specs so there's no duplication when species share
  covariates.

* **Priors.** First-class via `prior_normal()`,
  `prior_lognormal()`, `prior_gamma()`, `prior_beta()`. The same
  constructors are available unprefixed (`normal()` /
  `lognormal()` / ...) **only inside** the `priors = ...` argument
  via a private NSE data mask, so user code stays close to
  mathematical notation without masking
  `base::gamma()` / `base::beta()` at the package level. Priors
  can be a single value applied to every species, or a named list
  keyed by species id (and shortly, by `(species, sex)`).

* **Bounds.** Per-row `lower` / `upper` flow into
  `build_bounds()$lower$beta_linkage` /
  `build_bounds()$upper$beta_linkage`.

* **Growth** (von Bertalanffy / Richards) is the first process
  fully wired to the new pipeline. `build_growth()` gains a
  `linkages` argument and a string-named `fun`
  (`"empirical"` / `"vonBertalanffy"` / `"Richards"`); integer
  codes still work (`fun = 1` is shorthand for
  `fun = "vonBertalanffy"`) so existing scripts don't need to be
  rewritten apart from substituting `fun =` for `growth_model =`.

  ```r
  build_growth(
    fun = "vonBertalanffy",
    linkages = list(
      log_K = linkage_spec(
        formula = ~ temp,
        by      = ~ species + sex,
        priors  = list(temp = normal(0, 1))
      )
    )
  )
  ```

* **TMB plumbing.** New `src/TMB/linkage.hpp` accumulator;
  `ceattle_v01_11.cpp` reads parallel `DATA_IVECTOR(linkage_*)`
  inputs plus a `DATA_MATRIX(linkage_X)` and writes a
  `growth_linkage_offset` tensor that is added (additively, on the
  log scale) to `growth_parameters`. Per-row prior densities
  contribute to slot 19 of the joint NLL ("Linkage-table priors"
  in `fit$quantities$jnll_comp`).

* **Documentation.** New vignette
  `vignette("environmental-linkages", package = "Rceattle")` walks
  through the API, prior families, species-keyed priors,
  per-species formulas, basis-expansion formulas, and the
  underlying pipeline.

* **Natural mortality** is the second process wired to the
  pipeline. `build_M1()` gains a `linkages` argument keyed by
  `log_M1`; the offset is added on the log scale to `log_M1` inside
  the `M1_at_age` compute. A row's `age_bin == NA` broadcasts the
  offset across ages; specific values pin it to that age slice.
  `build_M1()` also gains string-form acceptance for `M1_model`
  and `M1_re` (parity with `build_growth(fun)`):

  ```r
  build_M1(M1_model = "sex_age_invariant",   # or 1
           M1_re    = "ar1_age",             # or 4
           linkages = list(
             log_M1 = linkage_spec(formula = ~ temp, by = ~ species)
           ))
  ```

  Growth and M can be linked in the same fit; their rows share the
  same global linkage table and the same `beta_linkage`
  parameter vector.

* **Per-(species, sex) priors.** In addition to scalar and
  species-keyed priors, each `priors[[col]]` value may be a
  two-level nested list keyed first by species id then by sex id:

  ```r
  priors = list(
    temp = list(
      `1` = list(`1` = normal(0, 0.1),
                 `2` = normal(0, 0.2)),    # sp 1 by sex
      `2` = normal(0, 0.5)                  # sp 2, both sexes
    )
  )
  ```

  Missing keys at either level resolve to "no prior" for that
  stratum. The validator checks the structure recursively and
  emits actionable error messages keyed by
  `priors$<col>$<species>[$<sex>]` paths.

* **Default `by = ~ species`.** `linkage_spec()` now defaults the
  `by` argument to `~ species`, so each linkage produces one
  coefficient per species without the user having to spell it out.
  Pass `~ species + sex` for per-(species, sex) coefficients, or
  `by = NULL` to share a single coefficient across every
  species/sex (the prior default). This matches the typical
  multispecies assessment use case where each stock has its own
  environmental sensitivity.

* **Recruitment** is the third process wired to the pipeline.
  `build_srr()` gains a `linkages` argument keyed by `log_R0`,
  `log_alpha`, or `log_beta`; the offset is added on the log scale
  to the corresponding parameter at every recruitment compute
  call site (hindcast, BRPs, projections, expected R). `log_R0`
  is meaningful for any `srr_fun`; `log_alpha` and `log_beta`
  only do work for SRRs that consume them
  (`srr_fun in c(2, 3, 4, 5)` -- Beverton-Holt and Ricker).

  ```r
  build_srr(
    srr_fun  = 2,
    linkages = list(
      log_alpha = linkage_spec(formula = ~ temp,
                               priors  = list(temp = normal(0, 0.5)))
    )
  )
  ```

  Growth, M, and recruitment can be linked in the same fit; their
  rows share the same global linkage table and the same
  `beta_linkage` parameter vector. End-to-end tests in
  `tests/testthat/tests-Dynamics/test-linkage-auto-map.R` verify that
  the base parameter (e.g., `log_growth_pars`, `log_M1`, `rec_pars`)
  is automatically mapped out (set to `NA`) when a linkage is active
  for that parameter, allowing the linkage intercept to define the base value.
  
  `tests/testthat/tests-Dynamics/test-recruitment-linkage.R`
  cover the analytical relations
  `R = R0 * exp(beta * temp[yr])` (mean R) and the
  `growth + M + recruitment` composition.

* **Soft deprecation in `build_M1()`.** The legacy column-index
  argument `M1_indices` and the env-driven structural integer
  codes `M1_model %in% c(4, 5)` are subsumed by the new
  `linkages = list(log_M1 = ...)` argument. Both still work for one
  release cycle, but emit a one-time warning that points users at
  the equivalent linkage-table call. The string aliases
  `"env_sex_invariant"` / `"env_sex_specific"` (added briefly on
  the dev branch) are removed; they were never released. The cpp's
  `M1_beta` / `M1_mult` code path stays in place for now -- both
  paths add additively to `log_M1` on the log scale -- but do not
  use both for the same coefficient or you will double-count.

* **Roadmap.** Recruitment is next on the same pipeline; then
  random-effects pooling on `re_group` for hierarchical
  shrinkage. The legacy `M1_indices` / `M1_model = 4|5` paths
  retire when recruitment migrates.

## Scheduled removal (v4.5.0)

> **Schedule update (v4.4.0):** The removal originally targeted for
> v4.2.0 has been pushed to **v4.5.0** to give downstream users a
> longer migration window for the natural-scale linkage API rolled
> out in v4.4.0. The soft-deprecation warnings continue to point
> users at the equivalent linkage-table call. The cleanup checklist
> below is unchanged.

The soft-deprecated API surfaces below remain functional and emit
one-time warnings pointing users at the linkage table. They will be
**removed entirely in 4.5.0**. To migrate, replace:

| Legacy                                | New                                                      |
|---------------------------------------|----------------------------------------------------------|
| `build_srr(srr_indices = ...)`        | `build_srr(linkages = list(log_R0 = linkage_spec(...)))` |
| `build_srr(srr_fun = 1)`              | `build_srr(srr_fun = 0)` + linkage on `log_R0`           |
| `build_srr(srr_fun = 3)`              | `build_srr(srr_fun = 2)` + linkage on `log_alpha`        |
| `build_srr(srr_fun = 5)`              | `build_srr(srr_fun = 4)` + linkage on `log_alpha`        |
| `build_M1(M1_indices = ...)`          | `build_M1(linkages = list(log_M1 = linkage_spec(...)))`  |
| `build_M1(M1_model = 4)`              | `build_M1(M1_model = 1)` + linkage on `log_M1`           |
| `build_M1(M1_model = 5)`              | `build_M1(M1_model = 2)` + linkage on `log_M1`           |

**Cpp cleanup checklist** (search for `LEGACY` in
`src/TMB/ceattle_v01_11.cpp`):

* Drop `PARAMETER_MATRIX(beta_rec_pars)`, `PARAMETER_ARRAY(M1_beta)`,
  and the scratch vectors `srr_mult`, `beta_rec_tmp`, `env_rec_tmp`,
  `M1_mult`, `beta_M1_tmp`, `env_M1_tmp`.
* Delete the five `srr_env_mult` blocks (hindcast, BRPs, dynamic
  BRPs, projection, R_hat) and the `M1_mult.sum()` term inside the
  M1_at_age compute.
* Pass `Type(0.0)` for `srr_env_mult` at each
  `calculate_recruitment()` call site (or drop the parameter from
  the function signature in `recruitment.hpp` if no caller still
  needs it).

**R-side cleanup**:

* Remove `srr_indices` and `M1_indices` arguments from
  `build_srr()` and `build_M1()`.
* Reject `srr_fun %in% c(1, 3, 5)` and
  `M1_model %in% c(4, 5)` as unknown integer values in
  `.coerce_srr_fun()` and `.coerce_M1_arg()` (drop the
  `.SRR_DEPRECATED_FUNS` / `.M1_DEPRECATED_MODELS` constants).
* Remove the `suppressWarnings()` wrappers around internal
  `build_srr()` / `build_M1()` re-callers in `sim_mod()`,
  `retrospective()`, `jitter()`, `run_mse()`, `project_no_F()`.

# Rceattle 4.0.3

## API

* New `fit_control()` constructor bundles the optimizer / sdreport /
  phasing knobs that previously cluttered `fit_mod()`'s signature
  (`phase`, `getsd`, `bias.correct`, `use_gradient`, `rel_tol`,
  `loopnum`, `newtonsteps`, `getJointPrecision`, `getReportCovariance`,
  `verbose`, `TMBfilename`, `nlminb_control`). Pass the result via the
  new `fit_control` argument:

  ```r
  fit <- fit_mod(
    data_list   = BS2017SS,
    msmMode     = 0,
    fit_control = fit_control(phase = TRUE, getsd = FALSE, loopnum = 1)
  )
  ```

  `fit_mod()`'s visible argument list shrinks from ~33 to ~22 args, so
  calls now read as model spec rather than a pile of optimizer flags.

* `fit_mod()` emits a deprecation warning if any of the legacy
  control args are passed directly and forwards them into
  `fit_control` for the duration of the deprecation window. Truly
  unknown arguments still error with `Unused arguments to fit_mod():
  ...` (no silent drops).

* Internal callers (`run_mse()`, `retrospective()`, `jitter()`,
  `sim_mod()`, `project_no_F()`) now wrap their control args in
  `fit_control(...)` rather than passing them positionally.

## New methods

* S3 methods on the `"Rceattle"` class so a fit behaves like an R
  model object: `plot()`, `coef()`, `vcov()`, `logLik()`,
  `residuals()`. With `df` set on `logLik`, [stats::AIC()] also
  works without a dedicated method. `nobs()` is intentionally
  *not* defined: counting "observations" in a stock-assessment
  likelihood (composition cells, indices, catches, priors) is not
  well-defined, so [stats::BIC()] does not work — use AIC or
  domain-specific information criteria.

* `plot.Rceattle()` is a thin dispatcher: `plot(fit, what = "biomass")`
  / `"ssb"` / `"recruitment"` / `"depletion"` / `"index"` / `"catch"` /
  `"selectivity"` / `"mortality"` / `"data"`. `...` is forwarded to
  the underlying `plot_*()` function.

* `residuals.Rceattle(type = ...)` returns a long-format data frame
  with rows from one or more of the four fitted data sources:
  `"index"` and `"catch"` (log-scale by default; switch with
  `scale = "natural"`), `"comp"` (Pearson on fitted proportions, with
  the `Age0_Length1` flag preserved), and `"caal"` (Pearson on
  fitted proportions, with both the conditioning `Length` and the
  age `Bin`). `type = "all"` returns all four stacked.

## Documentation

* README now has a self-contained *Getting started* block that fits a
  bundled model and exercises every new S3 method, so first-time
  users on the pkgdown site / CRAN no longer have to bounce to the
  GitHub wiki to see a working example.

* Vignette 8 ("Model parameterizations") is being expanded to fill in
  coverage gaps surfaced during the 4.0.2 release audit:
  * M1 random-effects modes (`M1_re = 0..6`).
  * Full stock-recruit section (Mean / BevertonHolt / Ricker, env-linked
    forms, Beta prior on steepness via `srr_est_mode = 3`).
  * Composition likelihoods (Multinomial, Dirichlet-multinomial, CAAL).
  * Catchability equations for `Catchability = 4` (Power),
    `5` (Environmental), and `6` (AR1).
  * Selectivity equations for `Selectivity = 6` (2DAR1) and
    `7` (3DAR1).
  * Internal growth model (`growth_model = 1` von Bertalanffy /
    `2` Richards).
  * Initial-age-structure modes (`initMode = 0..4`).
  * Harvest control rules (HCR = 0..7) — possibly cross-linking
    vignette 0 §9 rather than duplicating.

# Rceattle 4.0.2

## Bug fixes

* `data_check()` now stops with a clear message if a user requests
  `msmMode = 3:9` (Kinzey & Punt 2009 functional responses --
  Holling I/II/III, predator interference, predator preemption,
  Hassell-Varley, Ecosim). Those code paths are unvalidated against
  the current parameter set and should not be used for advice.
* Plotting functions (`plot_timeseries`, `plot_selectivity`,
  `plot_mortality`, `plot_maturity`, `plot_b_eaten`, `plot_b_eaten_prop`,
  `plot_m_at_age`, `plot_m2_at_age_prop`, `plot_f`, `plot_ration`,
  `plot_index`, `plot_catch`, `plot_logindex`, `plot_indexresidual`,
  `plot_comp`, `plot_selectivity_vs_maturity`, `plot_stock_recruit`)
  now restore graphics `par()` on exit instead of mutating the
  user's device permanently.
* Replaced `T` / `F` shortcuts with `TRUE` / `FALSE` throughout the
  package source (~60 occurrences).
* Replaced `.data$Bin` / `.data$Length` inside `tidyr::pivot_wider`
  arguments with quoted strings (tidyselect 1.2.0 deprecation).
* `examples/Georges_bank_example.R` now calls `plot_mortality()`
  instead of the long-removed `plot_mort()`.
* `_pkgdown.yml` "Get started" / overview navbar links now point at
  the actual generated `articles/Rceattle-overview.html` (was
  `Rceattle-overview-4_17_2025.html`, which 404'd).
* `README.md` example links now reference the correct
  underscore-separated filenames (`Fit_2018_GOA_multi-species_model.R`
  etc.) -- the previous space-encoded URLs 404'd on GitHub.
* Removed duplicate `\seealso{}` block in `?Rceattle-package` by
  consolidating `R/Rceattle.R` into `R/Rceattle-package.R`. Both
  files declared `_PACKAGE`, so roxygen2 emitted the auto-generated
  links twice.
* Added `graphics::box` to package imports (cleared the lone
  `R CMD check` NOTE for `plot_data`).
* TMB: terminal-age geometric series now includes `Finit` in the
  denominator for fished initial-equilibrium modes
  (`initMode = 3` or `4`), correcting a bias in the plus-group
  N-at-age when `Finit > 0`.

## Documentation

* Added Wassermann et al. (2025) cannibalism / Pacific hake
  reference to `inst/CITATION` and `?Rceattle-package`.
* `initMode` accepts integer codes or string aliases.
* `HCR` accepts integer codes or string aliases.

## Tests

* `tests/testthat/test-subdir-folders.R` no longer calls
  `testthat::test_dir()` for each subdirectory. Nested `test_dir()`
  inside `test_check()` triggered an "evaluation nested too deeply:
  infinite recursion" abort inside `rlang`'s trace deparser whenever
  any test failed, masking the real failure. Subdirectory test files
  are now discovered with `list.files()` and pulled in via
  `source()` so they register against the outer reporter directly.
* `tests/testthat.R` now wraps `library()` calls in
  `suppressPackageStartupMessages()` so transient build-version
  notices (e.g. "package 'dplyr' was built under R version 4.5.2")
  do not get captured as test warnings whose backtraces then crash
  rlang's expr_deparse at end-of-run.
* `tests/testthat/test-Dynamics/test-initial-dynamics` evaluates the different starting conditions.

## Parallelism

* `run_mse()` and `check_mse()` now use `parallel::parLapply` on a
  PSOCK cluster instead of `foreach::foreach(...) %dopar%`.
  - The `%dopar%` path triggered `rlang::expr_deparse` infinite
    recursion under nested `test_that` backtraces because each
    `foreach` invocation captures call frames that recurse during
    error formatting. PSOCK workers are clean R processes with no
    captured promise chains, so the issue does not occur.
  - PSOCK clusters work identically on Windows and macOS/Linux.
  - `run_mse()` gains a `cores` argument (default `NULL` picks
    `parallel::detectCores() - 6`); both functions cap at 2 cores
    when `_R_CHECK_LIMIT_CORES_` is set so they comply with CRAN's
    R CMD check limit.
  - `foreach` and `doParallel` removed from `Imports:`.

## Installation / dependencies

* `TMBhelper` moved from `Imports:` to `Suggests:`. Rceattle now uses an
  internal `.fit_tmb()` wrapper that delegates to `TMBhelper::fit_tmb()`
  when the (optional, GitHub-only) package is installed and otherwise
  falls back to a `stats::nlminb()` + `TMB::sdreport()` path. This
  removes the largest install-friction barrier for new users.
* `dplyr` moved from `Depends:` to `Imports:`. The package no longer
  attaches `dplyr` to the user's search path on load (so it no longer
  masks `stats::filter` / `stats::lag`).
* `kableExtra` dropped from `Suggests:`; vignettes now use
  `knitr::kable()` for table rendering.
* `quarto` dropped as a vignette engine; all vignettes are now `.Rmd`.

## API

* `run_mse()` no longer has `om = ms_run`, `em = ss_run` defaults.
  Both arguments are now required and validated as objects of class
  `"Rceattle"` before the MSE loop runs. Calling `run_mse()` with no
  arguments previously produced a confusing "object 'ms_run' not found"
  error; it now stops with a clear message.

## New methods

* `print.Rceattle()` and `summary.Rceattle()`. Auto-printing a fit
  inside knitr / RStudio / R Markdown previously dumped tens of MB of
  nested data and could trigger deep recursion errors during vignette
  rendering.

## Build / packaging

* Tightened `.Rbuildignore` (excludes `examples/`, `R/dev/`,
  `src/TMB/Dev/`, `.Rhistory`, `.claude/`, `.DS_Store`, build tarballs,
  `.Rcheck` directories).
* Tightened `.gitignore` to catch all `*.o` / `*.so` / `*.dll` and
  `*.Rcheck/` directories.
* Suppressed a benign clang `-Wfixed-enum-extension` warning by
  scoping the diagnostic pragma around `#include <TMB.hpp>` rather
  than via global Makevars flags.

# Rceattle 4.0.1

The 4.0.1 development cycle reorganized several `data_list` columns
and `fit_mod` / `build_*` arguments. Models or data files saved
against earlier 4.x revisions may need updating; see the renames
below. Compiled from
`inst/Running_list_of_updates.qmd` plus the `dev` branch commit log.

## Data renames

* `Pyrs` -> `ration_data` (the old name is still accepted on read,
  but is silently renamed).
* `UobsWtAge` -> `stom_prop_data`.
* `fsh_biom` -> `catch_data`.
* `srv_biom` -> `index_data`.
* `Nselages` -> `N_sel_bins` (in `fleet_control`).
* `Sel_norm_bin1` / `Sel_norm_bin2` <- `Age_max_selected` /
  `Age_max_selected_upper` (selectivity normalization bins).
* `Age_first_selected` -> `Bin_first_selected` (in `fleet_control`).
* `sel` -> `sel_at_age` (model report).
* `fleet_control` now carries a `Month` column (month of observation
  for indices / fisheries).

## API renames

* `build_M1`: `M1_prior_mean` -> `M_prior`,
  `M1_prior_sd` -> `M_prior_sd`.
* `build_srr`: `srr_prior_mean` -> `srr_prior`;
  `R_hat_endyr` replaced by `srr_hat_styr` / `srr_hat_endyr`.
* `fit_mod`: `suit_meanyr` replaced by `suit_styr` / `suit_endyr`.
* `initMode` semantics revised: 0 = free-parameter N-at-age,
  1 = unfished equilibrium with no devs, 2 (default) = unfished
  equilibrium with initial devs, 3 = fished equilibrium with initial
  devs. Type 4 ("non-equilibrium scaled") added later.

## New features -- composition and diet likelihoods

* Dirichlet-multinomial composition likelihood. Selected per fleet via
  `fleet_control$Comp_loglike = 1` (or `"DirichletMultinomial"`).
* Conditional age-at-length (CAAL) data path, with `CAAL_loglike` /
  `CAAL_weights` controls in `fleet_control`. CAAL data also flow
  through `sim_mod()` for simulation testing.
* `Diet_loglike` switch on the bioenergetics control sheet selects
  between multinomial (0) and Dirichlet-multinomial (1) for diet
  composition.
* Other-food diet proportion estimates added to the model report.
* Weighted-mean diet data path (annual proportion of prey-at-age in
  predator-at-age averaged across years).

## New features -- selectivity, catchability, growth

* Hake non-parametric selectivity (`Selectivity = "Hake"` or `5`),
  after Taylor et al.
* `2DAR1` (`= 6`) and `3DAR1` (`= 7`) selectivity
  parameterizations, after Cheng et al. (2024).
* `Catchability = 6` ("AR1"): annual AR1 catchability deviates fit
  to an environmental index, after Rogers et al. (2024) for the GOA
  pollock model. Environmental q-link (`Catchability = 5`) also
  exposed.
* Internal growth model. See `build_growth()` and the `growthFun`
  argument to `fit_mod()`. `alpha_wt_len` / `beta_wt_len` added to
  the data control sheet. Length-based suitability (`suitMode = 1` /
  `2` / `3` / `4` / `5` / `6`) wired through to use the estimated
  growth model. Comparison with WHAM growth implemented under
  `tests/comparison/`.
* Predator-specific suitability mode (different `suitMode` per
  predator).
* Suitability calculation now uses configurable year ranges
  (`suit_styr` / `suit_endyr`) instead of "mean year".

## New features -- recruitment and reference points

* Beta-distributed prior on Beverton-Holt steepness, available via
  `srr_est_mode = 3`.
* M1 random effects with optional environmental linkage; `M_prior` /
  `M_prior_sd` priors carried through `build_M1()`.
* `remove_F()` function returns a fitted model with F set to 0 --
  used internally for dynamic reference point calculation.
* `DynamicHCR = TRUE` in `build_hcr()` to switch from static to
  dynamic SB0 reference points.
* CMSY harvest control rule (`HCR = 1`): maximize joint catch across
  species, optionally constrained to keep depletion above `Plimit`.
* PFMC Category 1 40-10 ABC HCR (`HCR = 6`) using `Pstar` /
  `Sigma` uncertainty buffer.
* SESSF Tier 1 HCR (`HCR = 7`).
* Iterative multi-species HCRs: `HCRorder` controls the order in
  which species F is solved (e.g. predators before prey) inside
  `build_hcr()`.

## New features -- MSE and projection

* `run_mse()` now writes per-simulation `.rds` files when `dir` is
  specified, for streaming-friendly long runs. `load_mse()` reads
  those back.
* `check_mse()` validates which OM/EM simulations converged.
* `mse_summary()` produces a per-fleet performance-metric table
  (mean catch, IAV, P(closed), MSE on SSB, P(F > Flimit),
  P(SSB < SSBlimit), terminal depletion, ...).
* MSE function now supports `cap` (catch cap), `catch_mult` (catch
  multiplier), `rec_trend` (linear projected recruitment trend),
  `fut_sample` (future sampling effort), per-fleet
  `assessment_period` / `sampling_period`, `regenerate_past` (refit
  EM to OM-simulated past data), and `timeout`/`try`-error handling
  per simulation.
* `Recruitment_and_fixed_F_projections.R` and `Simulation_testing.R`
  examples added.

## New features -- diagnostics and tooling

* `jitter()` function to perturb starting values and re-fit, for
  global-vs-local-minimum diagnostics.
* `retrospective()` peels with optional `nyrs_forecast`.
* `model_average()` for averaging derived quantities across multiple
  fitted models, with optional bootstrap uncertainty.
* `compare_sim()` and `sim_mod()` for parametric simulation testing.
* `McAllister-Ianelli-reweighting.R` example for composition
  reweighting.
* TMB log-likelihood pieces (unweighted) added to the report for
  composition diagnostics.
* `Selectivity = "Fixed"` (`= 0`) for empirically supplied selectivity
  blocks via the `emp_sel` data sheet.
* `TMBfilename` argument to `fit_mod()` to point at an alternate
  `.cpp` during development.

## Behavior changes

* Removed accumulation-age switches in `fleet_control`. Selectivity
  normalization is now controlled via `Age_max_selected` (i.e.
  `Sel_norm_bin1`) on a per-fleet basis instead of always
  normalizing by the maximum-selectivity age.
* `NA` values inside the valid age/length range of composition data
  are now coerced to 0 with a warning (previously silently dropped
  or errored).
* Selectivity dimensioning switched from age- to bin-indexed for the
  non-parametric and 2D/3D AR1 forms (driven by
  `Selectivity_dimension` and `N_sel_bins`).
* Age-error and age-transition matrices are now dimension-checked
  against `nages` at `data_check()` time.

# Rceattle 4.0.0

* CEATTLE TMB version 4.0.0. See Adams et al. (2022),
  *Fisheries Research*, 251, 106303.
