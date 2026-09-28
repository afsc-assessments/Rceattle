# Selectivity: what is open, and what shipped

Replaces the four selectivity notes that stood before 5.41.0, whose designs are now either
shipped or recorded here. Measured numbers are kept; the design prose that led to a shipped
feature is not.

---

## Open 1 — `Hake` selectivity ignores `Sel_norm_scope`

**Status: proposed, not implemented.** `selectivity.hpp` carries no `PROPOSED` block for this
(the five in `ceattle.cpp` belong to the superseded non-parametric penalty design), but it does
carry a comment naming this file at `src/TMB/selectivity.hpp:531`, sitting where the fix goes.
That comment is the anchor; grep for the filename, not for a marker.

`Hake` (type 5) takes its normalization reference inside the sex loop, so the reference is
always that sex's own maximum, and it is excluded from the shared normalizer
-- `normalize_and_project_selectivity()` is called for EVERY fleet (`selectivity.hpp:660`), but
its normalization block is gated on `sel_type` not in 5, 11, 12 (`:74`; no form maps to 12
today). Its other two blocks do run for Hake: zeroing bins below `Bin_first_selected` (`:55`)
and, load-bearing for the fix below, copying the terminal hindcast curve into every projection
year (`:120`). `Sel_norm_scope` is therefore inert on a Hake fleet. Measured on `GOAatf` fleet 3 (`nsex = 2`,
`Time_varying_sel = "IID"`, `Sel_norm_bin = 0`), the two settings agree to every digit reported:

| `Sel_norm_scope` | female max | male max | ratio | jnll |
|---|---|---|---|---|
| `"AcrossSexes"` (default) | 1.0000 | 1.0000 | 1.0000 | 46.4710 |
| `"WithinSex"` | 1.0000 | 1.0000 | 1.0000 | 46.4710 |

Two consequences. Setting the column changes nothing. And because F is shared across sexes
(`F_flt_age = sel_at_age * exp(log_F)`), a Hake fleet cannot express sex-specific fishing
mortality at all, whatever the data say.

**Pooling the reference is not sufficient on its own, and the reason is specific to this form.**
Hake's curve is a cumulative sum of `sel_coff` on the log scale, and the first coefficient is
mapped off at 0 for BOTH sexes (`R/3-build_map.R:1141`, "first parameter is not-identifiable and
is not estimated"; `sel_coff` inits to 0). So each sex's raw log-curve starts at exactly 0, and
the normalization decides **where the two sexes are welded together**:

| | reference subtracted | sexes forced equal at | free |
|---|---|---|---|
| `"WithinSex"` (today) | each sex's own max | the **peak** -- both are 1 | the ratio at young ages |
| `"AcrossSexes"` | one pooled max `M` | the **first selected bin** -- both `exp(-M)` | the ratio at the peak |

Neither frees the male:female level outright; pooling only relocates the pin. **Across-sex
relocates it to the better place.** For a dimorphic stock the sexes are similar in size where
they first recruit to the gear and diverge as females grow, so equal selectivity at the first
selected bin is defensible, where equal selectivity at the peak asserts they are equally
available exactly where they differ most. Measured on `GOAatf` under today's within-sex rule,
age-1 selectivity is 1.83e-4 for females and **2.25e-47** for males: with the peak welded to 1,
a near-step curve is the only way left for the model to say males are less available overall.
Suggestive rather than proof, but it is the shape the constraint predicts. This is also why
freeing `sel_coff`'s first bin for the second sex is tempting and still wrong -- see the fix
sketch.

Since 5.35.0 the silence is partly gone: `data_check()` messages a two-sex Hake or `LogisticPM`
fleet whose `Sel_norm_scope` is `AcrossSexes`, set or defaulted, that the column is not read and the sexes cannot differ
in level (`R/1-data_check.R:993`). A fleet explicitly set to `WithinSex` gets nothing, because
the behaviour happens to match what was asked for. The capability gap is untouched either way.

**A second, smaller gap sits in the same block.** Hake honours a single named bin
(`sel_norm_bin1 >= 0 && sel_norm_bin2 < 0`) but falls through to the maximum when a bin RANGE is
given, so `Sel_norm_bin_upper` is ignored too. The shared normalizer averages over the range.

**Fix sketch, kept because this one is not shipped.** Hake normalizes on the LOG scale
(`exp(val - max_sel)`) because its pre-normalization curve is a cumulative sum of `sel_coff`,
while the shared normalizer divides on the natural scale. Since `exp(a - m) == exp(a) / exp(m)`,
the pooled reference is `max2()` over sexes of the per-sex LOG reference, not a natural-scale
divisor. **So do not simply drop `sel_type != 5` from the normalizer's gate** -- that runs the
natural-scale division on top of Hake's own log-scale normalization and normalizes twice.
Instead move Hake's normalization out of the `switch` into a block after the sex loop closes but
still inside the year loop, where both sexes' un-normalized curves are in `sel_at_age` /
`sel_at_length`: compute `max_sel(sex)` per sex as now; if `Sel_norm_scope` is within-sex,
replace both with `max2()` of the two; then apply `exp(val - max_sel(sex))` per sex. Fold rather
than branch on the AD type, as the surrounding code does. While there, decide whether
`sel_norm_bin2` means a range mean -- on the log scale that is the geometric mean of the curve,
which is NOT what the shared normalizer computes -- and state the choice in
`R/0-column_schema.R`.

**Do not also free `sel_coff`'s first bin for the second sex to carry the level.** That is the
apical offset's job, and doing both gives two parameters for one quantity. Pooling here only
relocates where the sexes are welded together; it is not a substitute for the offset.

**What done looks like.** Without these the fix can ship leaving a vignette that contradicts the
code and no test that would have been red beforehand:

- No bundled golden model uses Hake (`BS2017SS` / `BS2017MS` are `Selectivity` 2/2/2/1 and
  `GOA2018SS` has no 5), so `/golden-check` should be bit-identical. Confirm that rather than
  assume it.
- The fit DOES move for any two-sex Hake fleet left on the default `AcrossSexes`, so this is a
  behaviour change: NEWS, `DESCRIPTION` `Version:`, and the `Sel_norm_scope` row in
  `vignettes/model-options-and-functionality.Rmd`, which currently states the limitation as
  permanent and would become wrong.
- Add the two-scope pair to `tests/testthat/` as a regression. Today they are identical; after
  the fix `AcrossSexes` must differ from `WithinSex` on `GOAatf` fleet 3. **The table above is
  the test.**

## Open 2 — `sel-penalty-form`, parked

Branch `sel-penalty-form` (`Sel_penalty_form`, 5 commits, 638 insertions) is parked by
decision, not by defect. It generalizes which penalty a non-parametric fleet is charged.

---

## Shipped — kept for the numbers and the caveats

### Per-sex apical offset (5.38.0, PR #149)

Selectivity is the only route to a sex difference in fishing mortality, and before this no form
had an apical height parameter: logistic forms asymptote to 1 for every sex, the non-parametric
forms re-centre each sex to mean 1 every year, Hake normalizes within sex, and `DoubleNormal`
peaks at exactly 1 for both.

Shipped as the selectivity linkage parameter `apical` (`log_sel_apical`, one row per
selectivity BLOCK not per fleet: `[n_selectivities, max_sex]`, which is why mirroring fleets
share it and a linkage on one is refused).
`.check_sel_apical_rows()` refuses no fleet, no sex, both sexes across rows, a one-sex species,
Fixed / AR1 / mirror fleets, an identity link, and within-sex normalization.

Two things to carry forward:

- **The multiplier is bin-wise.** It equals the ratio of the sexes' peak heights only where
  their shapes peak equally. Logistic on an age axis measured exactly 0.5; a double logistic
  with sex-specific descending inflections gave a peak ratio of 0.33 at a multiplier of 0.5.
- **The ratio is informed only by joint (`Sex = 3`) composition.** On a fishery with partial
  sex-structured sampling it is informed only over the years that have it. That is the case the
  parameterization exists for, and the reason it beats two free curves, but say so.

On `GOAatf`'s fishery the male multiplier fits at 2.15 (log 0.767, SE 0.213, positive-definite
Hessian, objective 356.56 with the `lognormal(0, 0.5)` prior against 364.01 without the offset).
Recovery over 120 replicates is consistent with unbiased (mean 0.685 against a true
0.693, empirical SD 0.28), which is not the same as demonstrating it.
Harness: `tools/verify/verify-sim-recovery-apical.R`.

**A prior on `inf_desc` is not a substitute**, though it was suggested as a stopgap. On
`GOAatf`, free per-sex parameters gave a male:female max ratio of 1.9687 and a tight
`normal(8, 0.1)` prior gave 1.9294: the prior constrained the shape, it barely moved the ratio.
A loose `normal(8, 1)` reached 2.80 with no sdreport.

**And read that 1.9687 against its own standard errors.** The same `DoubleLogistic` +
`AcrossSexes` fit on `GOAatf` fleet 3 (female max 0.5080, male 1.0000, jnll 348.979) returned
`log_sel_slp` SE 5391.6 and `sel_inf` SE 8878.7, with `hessian_conditioning` putting 96% of the
least-determined direction on `log_F`. **Quote a ratio from the sex whose parameters are
determined, and treat a standard error in the thousands as a boundary rather than an estimate.**
That confounding is what the apical parameterization sits inside, and it does not go away
because one parameter now carries the level.

**A naming trap worth stating before anyone wires this up.** Rceattle's `peak` is an alias of
`inf_asc`, the ascending inflection AGE; Stock Synthesis's "peak selectivity" is how HIGH the
curve is. Describing the existing `peak` linkage as the thing the apical offset does gives a
user a shifted curve, not a lowered one.

### Integrable non-parametric deviations (5.40.0, PR #151)

`NonParametric` charges its shape penalties on each year's realized curve, so under
`random_sel = TRUE` the Laplace approximation integrates a tilted density and the reported
deviation SD is not the SD of the deviations. Rather than move the penalty and change every
penalized AMAK fit, one form was added that carries a proper density:
`NonParametricIntegrable` (13), with `Time_varying_sel` picking the structure it covers.
`NonParametric` and `NonParametricPM` are bit-identical and still refuse `random_sel = TRUE`;
the refusal names the new form. **The two refusals have different reasons**, and the pointers
in `NEWS.md` and in `build_map()`'s roxygen land here for both: under `IID` the shape penalties are
charged on the realized curve, so the SD absorbs a term belonging to the normalizer; under
`RandomWalk` the per-year renormalization additionally leaves the LEVEL of each year's
coefficients improper, which no penalty change reaches. Form 13 under
`Time_varying_sel = "RandomWalk"` fixes the second by constructing the walk as a base plus a
running sum of increments and centring only for output.

**The superseded design is still in the template as comments.** Five `PROPOSED` blocks in
`ceattle.cpp` sketch moving the penalty onto `sel_coff` in place. That approach was rejected in
favour of new forms because it would move every penalized AMAK fit. Read them as history -- and
with the two approaches that were tried against the same problem and rejected, since the blocks
invite both:

- **Smoothing the hinge is not the fix.** Replacing `CppAD::abs(d)` with `sqrt(d*d + eps*eps)`
  gives a C-infinity penalty and cures the kink, but leaves the unnormalized prior untouched. It
  buys a converged fit with a wrong variance, which is worse than a refusal because nothing
  announces it.
- **Adding `log Z(sd, w)` to the objective is closed form for only half of it.** For the
  curvature term alone, a Gaussian times `exp(-w u'Du)` gives
  `log Z = -0.5 * log det(I + 2 w sd^2 D)`. Term 1's hinge has no closed-form normalizer, so
  that route needs the smoothing above AND a numerical constant.

**The guard was necessary but not sufficient, which is why new forms were the right answer.**
`.rce_np_unintegrable_fleets()` tested `Sel_curve_pen1` only. Zeroing it removes the kink and
the fit converges, but the average-selectivity term is hardcoded at weight 2.0 with no column
to disable it, charged every year on the realized curve. It is smooth, so it adds no kink and
all bias. `Sel_curve_pen2` behaves the same way. A fit that passed the guard still reported a
deviation SD maximizing a marginal missing its `log Z(sd, w)` term, with nothing announcing it.

**What the estimated SD is worth.** Measured with
`tools/verify/verify-sim-recovery-np-integrable.R` on `Atka2022`'s fishery (multinomial
compositions, input sample sizes 2 to 236): the estimate is biased **low**, 0.24 to 0.26
against a true 0.35 over two runs with every replicate below the truth, 13% low at 0.70, and 9%
low with the sample sizes multiplied by ten. It does not depend on the start. The scored set
equals the estimated set cell for cell and the normalizing constant is complete, so this is the
Laplace marginal likelihood's known downward bias for a variance component on small multinomial
samples (Breslow and Lin 1995), not a defect. Read the reported SD as a lower bound unless the
compositions are well sampled. Whether a bias-corrected Laplace or `tmbstan` should be the
recommended route for reporting it is **open**.

A related trap, fixed in the same PR: non-parametric coefficients below `Bin_first_selected` are
mapped off, but the curve centres each year by the log mean over every bin, so a value there
shifted the whole curve unscored. `inits` from a fit with a lower `Bin_first_selected` carry
exactly such values: 0.9 in those cells moved the Atka fishery objective 704 nats.
`fit_mod()` now holds them at 0.

### OSA outlier threshold (5.36.0)

Flagging `|resid| > 3` made a panel's expected outlier count grow linearly with its size, so a
length-composition panel always looked worse than an age panel on the same model. The OSA panel
now flags on `qnorm(1 - 0.05 / (2n))` for the panel's own `n`, which expects about 0.05 marked
points whatever the size. The Pearson panel keeps 3, documented as a heuristic, because those
residuals are sum-constrained and overdispersed rather than standard normal.

---

## From an external two-sex trial, settled

A user building a two-sex, random-effects assessment entirely in R raised four items. Items 1
and 2 became the apical offset and the integrable forms above; item 4 became the OSA threshold.
Item 3 was "combine the fleet control and `fit_mod` sections", and most of it already existed:
`build_data()` takes `fleet_control` and `model_config` as named blocks, and `run_config()` /
`save_config()` put the estimation controls in one YAML. The `config =` trap that silently
dropped a data object's linkages was fixed in 5.36.0, and the `random_sel` / `random_q`
documentation was corrected with it. What remains of item 3 is in `SIMPLIFY-LOG.md`: deriving
`random_sel` per fleet, and serializing `fleet_control` into the YAML.

Two things settled earlier that should not be redone: per-sex selectivity linkages fixing both
sexes was fixed in 5.28.1 (`564475b8`), and `random_sel` does not gate linkage random effects,
which enter through `beta_linkage_re` and integrate independently.

The correspondence itself is not reproduced here. The note that quoted it was build-ignored for
that reason, and this file is not.
