<!--
Version-numbering note. Three gaps in this file are deliberate, not lost entries:

  * 4.14.0 was a real DESCRIPTION version whose entries were folded into 5.0.0.
  * 5.2.0-5.2.4 were likewise folded into 5.3.0.
  * main's 4.9.0 / 4.9.1 are the same recruitment changes this line carries as
    5.5.0 / 5.5.1, applied to the two lines separately.

No tag existed above 4.8.0 while these were in flight, so nobody could have installed an
intermediate. They were folded rather than renumbered because renumbering a section breaks
every (x.y.z) cross-reference pointing at it, and the entries below cite each other by
version throughout.
-->

# Rceattle (development version)

Unreleased, on `dsem-v5-integration`. Deliberately not numbered: `dev` has
already published 5.24.0, 5.24.1 and 5.24.2, so a bare `5.24.0` here would make
an install from this branch report a version that exists and holds entirely
different code -- the mislabelling `95153bbc` fixed when it moved this branch to
a `.9000` suffix, and the GOA arrowtooth DSEM script installs from here. The
5.23.0 section below explains the same divergence for 5.19.0 and 5.22.0;
`DESCRIPTION` carries the development suffix until this work is given a release
number of its own.

## Diagnostics

* **A retrospective peel whose `sdreport()` failed presented a band of zeros as
  if it were real.** The handler returned `newmod$sdrep`, which at that point is
  the FORECAST refit's sdreport -- built with the whole hindcast pinned, so every
  hindcast standard error in it is exactly zero. That is the thing the report
  pass exists to replace, so the fallback handed back a zero-width band rather
  than no band. It now returns `NULL`, the state `getsd = FALSE` already
  produces, so `vcov()` and the plotters give NA bands.

  **A peel that previously drew a zero-width ribbon now draws none.** That looks
  like a loss in a figure and is not: the zero was never an estimate.

  The warning that should have said so was also unreachable. It was raised
  inside the per-peel closure, and a warning raised in a FORK/PSOCK worker is
  discarded, so at the default `cores` nobody was ever told. It is now collected
  on the peel and re-raised by the caller after dispatch -- the pattern
  `self_test()` uses for `sim_warns` -- and names the peel it came from. This was
  the last of the four worker-swallowed warnings in `retrospective()`; the other
  three were hoisted because their conditions read only the input model, while
  this one is genuinely per-peel and needed collecting instead.

  Point estimates are unaffected, and `tools/verify/verify-refit-like.R` is
  bit-identical across all ten sections before and after.

* **`hindcast_skill()` computed MASE over the wrong axis.** It grouped by PEEL
  and averaged over horizons; Kell et al. (2021) eq. 5 groups by HORIZON and
  averages over peels -- its sums run over \eqn{t = T-n \ldots T} at fixed
  \eqn{h}, which is why their Table 2 reports one MASE per steps-ahead rather
  than one per peel. The two are transposes, so every number the function
  returned was a mean over the wrong set.

  `$mase` is now keyed by `years_ahead` rather than `peel`, and carries
  `n_peels` -- how many peels were averaged at that horizon -- in place of
  `n_years`. **This changes the shape of an exported function's return value**;
  anything indexing `$mase$peel` or `$mase$n_years` needs updating. `$by_year`
  is unchanged. Species and quantity still split the table as before.

  Two pieces of documented advice went with the defect, because they were
  compensating for it rather than describing the statistic. The warning that a
  one-year MASE "is dominated by its own denominator" was true only because the
  mis-grouping made that denominator a SINGLE absolute error instead of a mean
  over peels; so was the advice to prefer deeper peels or a median across them.
  The real caution runs the other way: a horizon `h` can only be scored by a
  peel at least `h` deep, so `n_peels` FALLS as the horizon grows and the
  DEEPEST horizon in a run is the least precise. Read `n_peels` beside each
  MASE. A claim that the definition matched `ss3diags` was removed rather than
  restated, having never been checked against their source.

  More peels now buy precision inside each MASE rather than more rows, so
  `peels = 3:10` is a statement about cost -- two fewer model fits -- not about
  noise.

* **`retrospective(peels = )` takes a vector of peel depths, so a comparison can
  skip the shallow peels.** A single number still means every peel from 1 to `n`,
  as it always has; `peels = 2:10` now fits exactly those nine depths, and
  `peels = c(1, 5, 10)` three. Each element is a number of years removed from
  `endyr`, not a year, and depths are sorted and de-duplicated. This is worth
  having because the shallow peels are the ones that cost the most to interpret:
  a one-year peel's MASE in `hindcast_skill()` has a single
  \eqn{|naive - reference|} in its denominator, so it says more about the
  denominator than about the projection. `peels = 3:10` drops those and pays for
  two fewer model fits. `$peel_depths` records the depths asked for, alongside
  the existing `$peels_requested` count, and `print()` names them and attaches a
  `NOTE` whenever the set is not `1:n` — a rho averaged over a subset is not the
  quantity the +/-0.2 band was calibrated on.

  `peels = 5` and `peels = c(5)` are the same object in R, so a length-1 value
  always means `1:5`. Write `peels = c(5, 5)`, which de-duplicates, to fit one
  depth on its own.

  A bad `peels` is now refused before any peel is fitted rather than after. A
  depth of 0 made `(endyr_peel + 1):endyr` count *down*, so the retained years
  were fitted as a forecast; an over-deep one built a model ending before
  `styr`; and a non-finite or out-of-integer-range value (`peels = 1e10`)
  coerced to `NA` and surfaced as "missing value where TRUE/FALSE needed". All
  were silent, and all were paid for first.

* **Three once-per-call retrospective warnings never actually reached the
  caller.** Each was raised inside the per-peel dispatch, where a warning raised
  in a parallel worker is discarded — so at the default `cores` none of them
  surfaced at all. Two were additionally gated on `i == 1L`, which a subset such
  as `peels = 2:10` never reaches; the third sat in `.rce_peel_map()`, which has
  no peel index in scope, and so fired once per peel instead of once per call. They are now settled once
  before dispatch, where they belong — every one of the three conditions reads
  only the input model. `inst/dev/TRAPS.md` records the general trap.

  The three: **`forecast_rec = "model"` is inert when `proj_mean_rec = TRUE`**
  (`proj_mean_rec` takes precedence and `build_srr()` defaults it to `TRUE`, so
  a model fitted without naming it projects the peeled years at mean recruitment
  whatever process it carries — and `hindcast_skill()`, which defaults to
  `"model"` precisely to tell projection methods apart, silently cannot, which
  means a DSEM comparison could report "the correlation structure makes no
  difference" when the setting never engaged); **an unrecorded `random_vars`**,
  which announces a biased Mohn's rho and is the subject of the next entry; and
  **a fixed-effect linkage covariate** keeping its post-peel values, so rho is
  conditional on that covariate being known. That last one is gated on the
  DEEPEST peel, which is the one with the most post-peel rows; gating on the
  shallowest would go silent whenever `env_data` stops short of `endyr - 1` even
  though the deeper peels still see post-peel values.

* **A fit with no random effects no longer claims its retrospective is biased.**
  `fit_mod()` built `random_vars` as `c()` and stored it with
  `mod_objects$random_vars <- random_vars`, and assigning `NULL` into a list
  *deletes the element* — so an ordinary `random_rec = FALSE` fit came back with
  no `random_vars` at all, indistinguishable from a fit made before 5.10.0
  recorded it. `retrospective()` then told its owner that every peel pins every
  deviation and shrinks the process SDs with peel depth, and to refit. With no
  random effects, pinning every deviation is exactly what the fit itself did and
  there is no bias, so the warning was false on the most common single-species
  configuration.

  `fit_mod()` now stores `character(0)` for "recorded, and there were none", and
  only an absent record warns. **A fit returned by `estimateMode` 0 or 1 that
  declares no random effects therefore carries `random_vars` as `character(0)`
  where it previously had no such element** — observable on the returned object,
  though nothing in the package or the sibling repos reads it other than
  `retrospective()`. Nothing numeric moves: the `random=` declaration passed to
  `MakeADFun()`, the bounds and the convergence snapshot all read a local
  variable, not the stored one.

  A fit **saved** by an earlier version still has no `random_vars`, so it still
  draws the warning. That is correct for a pre-5.10.0 fit and a false alarm for
  an older no-random-effect fit; refitting resolves it either way.

* **`retrospective(forecast_rec = "model")` errored on any fit without a DSEM.**
  It called a helper `.pin()` that was never defined, so the branch reached
  whenever `forecast_rec = "model"` met `proj_mean_rec = FALSE` on a non-DSEM
  fit died with "could not find function `.pin`". Since `hindcast_skill()`
  *defaults* to `forecast_rec = "model"`, that function was unusable on exactly
  the non-DSEM models built for the comparison it exists to make — those with
  `proj_mean_rec = FALSE`. With the default `TRUE` it took the mean-recruitment
  arm and worked. The predicate is
  now a `pin_block()` defined once in the dispatching frame and reading
  `object$random_vars`, the same source `.rce_peel_map()` uses.

# Rceattle 5.23.0

This release carries the DSEM integration that was filed as 5.19.0 while this
branch was separate from `dev`, together with the simulation and sampling
fixes below. `dev` released its own, different 5.19.0, 5.21.0 and
5.22.0 in the meantime; none of them is this work.

## Bug fixes

* **An estimated SEM variance starting at 0 crashed R instead of erroring.**
  `x <-> x, 0, sigma, 0` asks dsem to estimate a variance from a start of zero.
  The exogenous covariance is singular there, and `MakeADFun()` inside
  `dsem::dsem()` does not error -- it segfaults, and R aborts with no condition
  for anything to catch. `check_dsem_spec()` could not report it either: it
  screens `build_dsem_objects()`'s output, so the crash happened before it was
  reachable. `build_dsem_objects()` now refuses that shape up front and says
  what to write instead.

  The parameter NAME is the discriminator, not the start value. `x <-> x, 0,
  NA, 0` pins the variance at zero deliberately -- the variable becomes
  deterministic, dsem handles it, and `check_dsem_spec()`'s
  `covariate_variance` record is what reports it. Measured against dsem 3.0.0:
  a named parameter starting at 0 aborts R; `NA` at 0 returns a usable model on
  a covariate column and on a recruitment-deviation column alike. A missing
  start, a zero cross-covariance and a lag-1 two-headed path are all untouched.

* **A DSEM's process draw was discarded, and took the covariate series with
  it.** `calculate_dsem()` carried its own `do_simulate` branch, and it assigned
  the whole latent field (`x_tj = draw_tj` in `dsem.hpp`), so it overwrote the
  conditional draw `sim_mod()` had just substituted -- adding 100 to every
  latent state changed nothing. It consulted neither the map nor `dsem_cond_k`,
  so it also redrew the covariate columns. Under `family = "fixed"` those
  columns **are** the `env_data` series, so every replicate was generated under
  an environment the refit was never shown: on a scaled covariate with sd 1 the
  drawn column moved by up to 2.96. The draw is now taken only in R, where it is
  conditional on the cells the map pins, and `family = "fixed"` covariates come
  back untouched.

* **`sim_mod(process = )` and `self_test(process = )` returned no truth under a
  DSEM.** The `rec_dev_drawn_sim` mask is written inside the template's own IID
  draw, which is gated off under a DSEM, so it stayed zero however much of the
  field was redrawn; `attr(x, "process_sim")` therefore came back `NULL` on
  exactly the models whose process **had** been redrawn. That is
  indistinguishable from nothing happening -- it is how the behaviour was
  reported -- and it stopped `compare_sim()` warning that it was measuring bias
  against an operating model that is no longer the truth. The draw now carries
  its own mask, and `process_sim` gains `dsem_x_tj` / `dsem_x_tj_drawn`: the
  latent states, not the deviations, are what says whether the SEM structure is
  identified.

* **`init_dev` is redrawn under a DSEM.** It was gated on `dsem_on == 0`, so a
  replicate's initial age structure kept its fitted values while the hindcast
  deviations were re-realized -- two different histories in one data set. The
  latent field starts at `styr` and section 13 scores `init_dev` with its own
  `N(-bias*R_sd^2/2, R_sd)` whatever `dsem_on` is, so this is the same rule the
  rest of the draw follows: draw what the density scores.

* **A DSEM under the AMAK/Ianelli recruitment penalty now refuses the draw, as
  the standard path already did.** `srr_fun = 0` with `srr_pred_fun > 0` scores
  the deviations a second time through the stock-recruit penalty, so there is no
  single distribution to draw from. The template's own draw has always been
  gated on that; the SEM draw was not, so it redrew a process the model scores
  twice. `sim_mod()` says why, as before.

* **`self_test()` no longer swallows the draw's warnings.** Replicates run
  through `parLapply`, and a worker's warnings never reach the caller, so at the
  default `cores` a `process` request that redrew nothing was completely silent
  -- no warning and no `process_sim`, with nothing to say why. The unique set is
  now collected from the replicates and re-emitted once (not `nsim` times), and
  asking for process error and getting none back warns in its own right.

* **A retrospective peel no longer shrinks a random-effect SELECTIVITY SD.**
  5.15.0 stopped pinning the peeled tail of a random-effect block -- a pinned
  deviation is still scored by its density, and "the deviation was exactly zero"
  is the strongest possible evidence for a small process SD -- but exempted
  selectivity wholesale, so `random_sel = TRUE` kept the full bias (-6.6% on
  sigma at 5 peels, monotone in peel depth, and therefore a trend in what Mohn's
  rho measures).

  The exemption had a real reason and the wrong scope. Several selectivity forms
  carry the ADMB/AMAK lineage's bare sums of squares rather than a normalized
  density, and marginalizing one of those contributes `-k * log(sigma)` to the
  objective, driving the SD UP instead -- the same bias in the other direction.
  But which kernel applies is a property of the FLEET, not of the parameter
  block, so one fleet's penalty was pinning every other fleet's deviations.

  The peel now decides per fleet, from the fleet's own `Selectivity`:
  `Hake` (`dnorm`), `2DAR1` (`SCALE(SEPARABLE(AR1, AR1))`) and `3DAR1` (`GMRF`)
  are freed; `NonParametric`, `NonParametricPM` and `LogisticPM` stay pinned.
  The limb deviates (`log_sel_slp_dev`, `sel_inf_dev`) are `dnorm` in every
  branch that estimates them except `LogisticPM`, which is the only form that
  pins them. A model can now carry a 3DAR1 fleet beside a NonParametricPM fleet
  and get the right treatment for each.

  Decided per `Selectivity_index` GROUP, because fleets sharing an index share
  one parameter block through shared map levels: pinning one member while
  freeing another would leave the pinned fleet's cells at their starting value
  while the shared level moved, so two fleets that must have identical
  selectivity would silently stop sharing it in the peeled years.

  `random_sel = FALSE` is unaffected -- nothing is a random effect, so every tail
  pins exactly as before. The map decision is now `.rce_peel_map()`, which
  `test-retrospective-selectivity-pinning.R` drives directly.

* **`process_residuals(process = "recruitment")` scored projection years on a
  DSEM.** With `build_DSEM(estimate_projection = TRUE)` -- which `sample_rec()`
  requires, so it is not an exotic setting -- the latent field runs to `projyr`,
  and every estimated recruitment state was whitened and returned. Those states
  are scored by no data, so their residuals are N(0, 1) *by construction*, and
  `osa_diagnostics()` pooled them into the SDNR: on BS2017SS with a five-year
  projection, 44 rows per species instead of 39, pulling the SDNR toward 1 for a
  reason that has nothing to do with fit. Hindcast years only now; the
  projection states stay in the conditioning set, which is the treatment an
  unobserved covariate cell already got. The function's own documentation always
  said one residual per hindcast year.

* **`hindcast_skill(reference = )` defaulted to both settings, not to `"model"`
  as documented.** `match.arg(several.ok = TRUE)` returns every choice when the
  argument is left alone. `"observed"` rebuilds each peel at
  `estimateMode = 3`, and it refuses outright on a model with analytical
  catchability -- so the documented default call did roughly twice the work it
  needed to, and failed outright on models the function can serve.

* **`retrospective(forecast_rec = "model")` says when it is inert.**
  `proj_mean_rec` takes precedence and `build_srr()` defaults it to `TRUE`, so a
  model fitted without naming it takes mean recruitment in the peeled years
  whatever process it carries -- measured on the GOA arrowtooth DSEM, Mohn's rho
  was identical to every digit under both settings. That is the documented
  precedence, but it silently defeats `hindcast_skill()`, whose whole purpose is
  telling projection methods apart. It now emits a message naming the switch and
  the way out.

* **A retrospective peel no longer refuses a model with a linkage on
  `env_data`.** Withholding post-peel covariate values was done by dropping the
  rows, and linkage design matrices are built POSITIONALLY from `env_data`
  (`materialize_linkage()` -> `model.matrix()`), so dropping rows dropped design
  columns and random-effect levels -- and `inits`, carried from the parent fit,
  no longer matched the model the peel built. On the GOA pollock 2025
  configuration, which carries `ar1(1 | Year)` and `rw(1 | Year)` linkages, every
  peel stopped outright: `beta_linkage` 241 against 237, `beta_linkage_re` 110
  against 108. The values are now blanked in place and the row count is
  preserved. A fixed-effect linkage covariate is exempt -- `materialize_linkage()`
  rejects an NA in one by design, because `model.matrix()` would silently drop
  the row and misalign the whole design -- and the peel says once which
  covariates it therefore still sees.

* **`check_dsem_spec()` no longer reports OK on a specification nothing has
  looked at.** It reads `sem_full`, which only a built object carries, so a
  `build_DSEM()` specification -- the natural thing to pass, given the argument
  is named `dsem` and that is what `fit_mod(dsem = )` takes -- returned an empty
  pass. It now names the call to make instead.


* **A retrospective peel no longer shrinks the estimated process SDs.** A peel
  does not shorten the model -- it turns data off after `endyr_peel` -- and it
  pinned the random-effect deviations in the peeled years: `rec_dev` at zero,
  and `log_M1_dev`, `index_q_dev` and the three selectivity blocks carried
  forward from the last fitted year. Those pinned values were still scored by
  their densities, which is fabricated rather than missing data, and maximally
  informative fabricated data: "the deviation was exactly zero", or for the
  carried-forward random-walk blocks "the increment was exactly zero", is the
  strongest evidence available for a small process SD. The estimate shrank as
  `sqrt((N - k)/N)` in the number of peeled years `k` -- on a 39-year model,
  -1.3% at one peel, **-6.6% at the default `peels = 5` and -13.8% at 10** --
  and because it grew with peel depth it was a trend across peels, in the same
  quantity Mohn's rho is computed from, leaving every peeled forecast
  overconfident.

  Deviations that are random effects are now left free in the peeled years, so
  the Laplace approximation integrates them out, which is what no data should
  mean. Deviations that are not random effects in a given model are still
  pinned; free would leave them unidentified. This was already how
  `beta_linkage_re` behaved -- it was never pinned -- and how a DSEM's latent
  states behave, and neither carried the bias, which is what identified it.

  Measured on BS2017SS: the deepest peel's recruitment SD now matches a DSEM
  peel of the same model to 1e-6 at both one and five peels, against ratios of
  0.985 and a predicted 0.934 before. Retrospectives and Mohn's rho on any model
  with estimated random effects will change; that difference is the bias being
  removed. Models with no random effects are unaffected.

* **A refit silently discarded a DSEM's fitted parameters.** `build_params()`
  does not declare the `dsem_*` blocks -- they come from the DSEM builder -- so
  `fit_mod()` dropped them out of any supplied `inits` and then refilled them
  from a freshly built template, i.e. their START values. This was invisible
  wherever the DSEM was re-estimated (any start reaches the same optimum) and
  wrong wherever it was held fixed. `retrospective()` holds the whole DSEM fixed
  in its forecast refit, so every peel ran with the recruitment SD at its start
  value: 0.707107 for every species regardless of the fit. On BS2017SS with a
  naive sem that put an 80% error in the peeled hindcast and flipped the sign of
  Mohn's rho (-0.49 against +0.13). Only models ESTIMATING the recruitment SD
  were affected; a sem that fixes it reads the fixed value and never saw this.
  `profile()` on a DSEM fit was affected the same way -- every grid point reset
  the model -- as was `estimateMode = 2`, which reprojected from the reset
  values.

* **`retrospective()`'s forecast bias adjustment did nothing on a DSEM.** It
  wrote `rec_dev`, which the template overwrites from the latent states on the
  next evaluation, so the retrospective forecast silently kept the GMRF's own
  projection instead of the bias-adjusted mean -- a 44% difference in forecast
  recruitment against the non-DSEM path. It now writes the latent state that
  `rec_dev` is derived from.

  Together these two restore the contract a naive (IID) DSEM is judged by: it
  now reproduces a non-DSEM retrospective (hindcast and forecast recruitment
  within 1e-3, Mohn's rho within 0.004, against 0.80 and 3.40 before).

* **`sim_mod(process = )` ignored its DSEM guard for every spelling but `TRUE`.**
  The guard tested `isTRUE(process)`, which is `FALSE` for `"all"`,
  `"dynamics"` and `"recruitment"`, so those returned a data set whose
  recruitment process had not been redrawn, with no error -- and under `"all"`
  the accompanying warning listed the other un-drawn processes, which read as
  confirmation that recruitment had been drawn. It now tests the resolved
  process vector. `self_test(process = )` had the same hole.

* **`retrospective(rescale = TRUE)` had no effect on a DSEM.** With
  `family = "fixed"` the covariate columns of the DSEM's latent state matrix ARE
  the environmental data, held fixed by the map, and their values reach the
  model as parameter starting values -- which a refit inherits from the parent
  fit rather than rebuilding. So the peel carried the parent's unrescaled
  covariate and the rescaling was silently discarded: every path coefficient and
  the whole latent state matrix came back bit-identical with and without it. The
  stale block is now dropped so the rebuild supplies covariates standardized on
  the retained years.

* **`build_DSEM(estimate_projection = FALSE)` now leaves the projection years out
  of the DSEM likelihood.** The latent states are built over `styr:endyr` only.
  They were previously built to `projyr` and held fixed, which still leaves them
  in the GMRF, where lagged paths tie them to the last hindcast years and pull on
  the terminal recruitment deviations -- and so on terminal SSB and the ABC. Fits
  with `projyr > endyr` will change; that difference is the bias being removed.

* **An `env_data` column with no values at all now warns.** There is no mean to
  fill the missing years with, so the column reaches TMB as `NaN` for every
  year, and any index pointed at it (`Cindex`, `M1_indices`, an
  environmental `Time_varying_q`, a linkage covariate) makes the objective
  `NaN` with nothing to say why.

* **A compiled model that does not match the R code now says so.** `MakeADFun()`
  reports a template variable the data and parameter lists do not supply as
  "Error when reading the variable: '<name>'. Please check data and
  parameters.", which reads like a problem with the data. The usual cause is a
  build mismatch -- the C++ was not rebuilt after an update, or an older DLL is
  still loaded in the session from before a reinstall -- and `fit_mod()` now
  says that alongside the original message.

## New features

* **A fit carries parameter bounds aligned to `obj$par`, in
  `$bounds$par_lower` / `$bounds$par_upper`.** This is what a sampler needs:
  `tmbstan::tmbstan(lower =, upper =)` takes a vector of exactly that length and
  expands it across the random effects itself, so the DSEM standard-deviation
  bound reaches an MCMC run without any change to the template.

  `$bounds$lower` / `$upper` are lists keyed by parameter BLOCK, one entry per
  block rather than one per estimated parameter, so the obvious construction --
  `unlist(bounds$lower[names(obj$par)])` -- repeats whole blocks and returns a
  vector orders of magnitude too long (33412 against an `obj$par` of 154, on
  `BS2017SS`). `fit_mod()` already had to solve this for `nlminb`, because
  `build_params()` order and the template's `PARAMETER_*` declaration order
  disagree; the aligned vectors are now kept rather than discarded.

* **A DSEM covariate observed with error is now simulated, state and
  observation.** `family = "fixed"` makes a covariate data: no measurement
  density, pinned in the map, held fixed -- unchanged. Any other family makes it
  a latent state observed with error, and neither half was drawn. Now the state
  is drawn with the rest of the field (it is process), and the observation is
  drawn from the family that scores it -- normal, Bernoulli, Poisson, Gamma,
  normal with known SD, lognormal or Tweedie -- and written back into
  `env_data`, so the refit is fitted to it. Previously a replicate carried the
  original environmental series beside freshly drawn recruitment, which
  overstates how much the SEM's paths can be recovered. Verified as a moment:
  over 200 draws the residual `y_sim - mu` had sd 0.9995 against a fitted
  measurement sd of 1.

* **`build_bounds(dsem = )` bounds a DSEM's standard deviations at 0, so the
  model can be sampled.** A lag-0 two-headed self path (`x <-> x`) is a diagonal
  of the Cholesky factor of the exogenous covariance. The likelihood sees only
  \eqn{\Gamma'\Gamma} and the template reads it as `sqrt(square(beta_z))`, so
  its sign is not identified and the surface is exactly symmetric about 0.
  Harmless for an MLE -- a mirrored pair of optima -- and fatal for MCMC, where
  the posterior is bimodal by construction, chains do not mix, and R-hat and any
  interval on sigma are meaningless. Bounding the diagonal below at 0 is the
  standard identifying restriction for a Cholesky factor and rules out nothing
  the likelihood could distinguish; pass `$lower` / `$upper` on to
  `tmbstan::tmbstan()` for the same constraint there. A variable that also
  carries a cross-covariance (`A <-> B`) or a lagged two-headed path keeps
  unbounded support, because identifying its sign means flipping a whole row of
  \eqn{\Gamma} rather than one element; `fit_mod()` names those variables. A fit
  carried forward from the negative root is flipped to the positive one rather
  than rejected as out of bounds -- exact where the bound applies, and it
  reproduces the same optimum. `dsem` is optional, so an existing
  `build_bounds()` call is unchanged.

  **`build_DSEM(bound_sd = )` turns it off.** `TRUE` by default, and a DSEM
  built before the argument existed reads as `TRUE`, so no fit changes
  parameterization silently. `FALSE` is for the two cases that want the
  unbounded surface: reproducing a fit made before the bound existed, parameter
  for parameter, and demonstrating that the two optima really are mirrored.
  It also switches off the sign flip, since rewriting a starting value the
  caller chose would defeat exactly those cases. The setting travels on
  `data_list$dsem_settings`, so `.refit_like()` and every diagnostic that
  refits inherit it.

* **`dsem_z_tj` is REPORTed**, so the latent field the model actually used is
  readable from R. `rec_dev` is derived from it, so it is what says whether a
  substituted draw reached the dynamics.

* **`summary()` reports a DSEM's path coefficients alongside the fixed effects.**
  5.17.0 made `summary()` return a `"summary.Rceattle"` object whose
  `$coefficients` is the fixed-effects table; a DSEM fit needs the SEM path table
  as well, and the two had the same name.

  On a DSEM fit `$coefficients` is the **path table** — the meaning the DSEM
  analysis scripts already read — and the fixed effects are `$fixed_coefficients`.
  On every other fit `$coefficients` is the fixed effects, unchanged from 5.17.0.
  `$dsem_coefficients`, `$fixed_coefficients` and `$recruitment_sd` each mean one
  thing on every fit and are what new code should read. `print()` shows the fixed
  effects, then the path table and the per-species recruitment SD.

* **`build_DSEM(constant_variance = )` chooses what a `<->` term means** once the
  SEM carries lagged paths, passed through to `dsem::dsem_control()`. `dsem`
  inherits RAM notation from the `sem` package, where a `<->` self-loop is the
  DISTURBANCE variance -- conditional on the incoming paths -- so the marginal
  variance grows along the causal graph, reaching `sigma^2 / (1 - rho^2)` for a
  first-order self-path. `"marginal"` and `"diagonal"` rescale `Gamma` and
  `(I - Rho)` to hold the marginal variance constant instead, which is the
  time-series convention `ar1(1 | Year)` already uses here, and WHAM's `Ecov`
  AR(1) and `glmmTMB`'s `ar1()` use elsewhere.

  **The default is unchanged** -- `"conditional"`, which is `dsem`'s -- so no
  existing fit moves. The argument only makes the choice reachable, and the
  three settings are identical for a SEM with no lags: measured on `GOApollock`,
  all three return an objective of 24246.37659. On a lagged SEM with rho = 0.6,
  the ratio of marginal variance to `R_sd^2` is 1.5625 under `"conditional"` --
  exactly `1/(1 - 0.6^2)` -- and 1.0000 under `"marginal"`.

  Two consequences of the default are now documented in `vignette("dsem")`
  rather than left implicit. The lognormal bias correction uses the MARGINAL
  variance (reported as `quantities$dsem_margvar_tj`), which is what makes
  `E[exp(rec_dev)] = 1` and is not `R_sd^2` once the SEM is lagged. And
  `sigmaR_prior_sd` centres on `data_list$sigma_rec`, an assessment's marginal
  recruitment SD, so under `"conditional"` that prior lands on the innovation
  SD -- the same quantity for an IID SEM, and different by `1/(1 - rho^2)` once
  a self-path is added.


* **A SEM written with its columns aligned now means what it reads as.**
  `dsem::make_dsem_ram()` decides whether a variable already has a variance of
  its own by comparing the path string it would add -- `paste(v, "<->", v)`,
  single-spaced -- against the paths already in the model, with no whitespace
  normalization. So `recdevs1  <->  recdevs1` failed to match, dsem appended a
  second variance row for the same node, and `MakeADFun()` then died inside the
  dsem DLL: a **segfault** on a multi-variable SEM, an error on a
  single-variable one, either way with a traceback that said nothing about
  spacing. `build_dsem_objects()` now puts exactly one space either side of
  every arrow before handing the SEM to dsem, and rejects a genuine duplicate
  variance term by name. Aligning a SEM's columns is the natural way to write
  one, and every SEM in the package's own tests, examples and vignette happened
  to be single-spaced, so nothing here could catch it.

* **`process_residuals(process = "recruitment")` works on a DSEM**, by whitening
  the latent states with the process precision instead of dividing by a scalar
  SD. Conditioning on the states the model was given -- the observed covariates,
  which are fixed in the map -- the estimated states have precision `Q[E, E]`,
  and `r = chol(Q[E, E]) %*% (draw - conditional mean)` is approximately iid
  N(0, 1) under a correctly specified SEM. That is the same claim the scalar
  path makes and exact in the same sense; it is the generalization of SAM's
  `procres()` to a correlated process. `"initial"` now reads the recruitment SD
  the density actually used (`quantities$R_sd`, which the template sets from the
  SEM) rather than `exp(R_log_sd)`, a mapped-out placeholder under a DSEM --
  measured on BS2017SS at 0.7071 against the 0.9888 / 1.3321 / 0.7976 in force.
  `"all"` therefore returns all three processes.

* **`sim_mod(process = TRUE)` and `self_test(process = TRUE)` redraw a DSEM's
  recruitment process**, by drawing the whole latent field from the same
  precision the density uses -- as `dsem::simulate()` does -- and substituting it
  into the parameter vector, so the template derives `rec_dev` from the drawn
  states and applies the lognormal bias correction once, where it always is. The
  states the model was given stay put: under `family = "fixed"` the covariate
  columns are the environmental data, and redrawing them would simulate
  different data rather than a different process. Drawing IID instead -- what
  the template's own process draw would have done -- would have discarded the
  SEM's lagged paths, covariate effects and cross-species structure, and used
  the GMRF's conditional SD, understating the marginal variance by
  `1 / (1 - rho^2)`.

* **`profile(param = "sigmaR")` profiles a DSEM's recruitment SD** instead of
  refusing. The alias now resolves to `dsem_beta_z[rec_sd_idx[sp]]` -- the SEM's
  variance path, on the natural scale, since `R_sd` is its absolute value --
  so the alias means the same thing on both parameterizations and `slots` still
  indexes by species. Profiling `rec_dev` is still refused: under a DSEM it is
  derived from the latent states on every evaluation, carries no gradient, and
  has no single parameter standing for it.


* **`retrospective(forecast_rec = )` chooses how the peeled years get their
  recruitment.** `"mean"` (the default, and what a retrospective has always
  done) projects them at the bias-adjusted historical mean. `"model"` uses the
  model's own rule instead: `proj_mean_rec = TRUE` projects at mean recruitment
  whatever process the model carries; otherwise the latent states supply it
  wherever the deviations are random effects (`random_rec`, or a DSEM), so an
  AR1's autocorrelation or a DSEM's lagged and covariate paths propagate into
  the forecast; otherwise recruitment comes off the stock-recruit curve.

  This exists because a peel's forecast years are HINDCAST years, so
  `proj_mean_rec` -- which the model only reads past `endyr` -- could not reach
  them, and every peel forecast at the historical mean whatever the model said.
  Measured on the GOA arrowtooth model, `proj_mean_rec = TRUE` and `FALSE`
  returned byte-identical skill scores before this. `hindcast_skill()` defaults
  to `"model"` for that reason; `retrospective()` keeps `"mean"` so Mohn's rho
  stays comparable to published values.

* **`hindcast_skill()` scores forecast skill across retrospective peels.** Each
  peel projects the years it did not see, given the catch actually taken, and is
  scored by mean absolute scaled error against either the full time-series
  model's estimate (`reference = "model"`, the default) or the survey index
  actually observed (`reference = "observed"`, classic hindcast
  cross-validation). The naive baseline is the peel's own terminal estimate held
  flat, so `MASE < 1` means the projection beats assuming nothing changes. This
  answers a different question from Mohn's rho, which measures how an estimate
  of a year moves as data accumulate rather than how well a year was predicted
  -- so it is the diagnostic for comparing recruitment projection assumptions,
  e.g. `proj_mean_rec = TRUE` against `FALSE` against a DSEM. One row per
  horizon per species, averaging across peels, so `n_peels` says how many peels
  a given steps-ahead rests on.

* **`fit_mod(dsem = build_DSEM(...))` fits a dynamic structural equation model
  on the recruitment deviations.** The deviations become the latent states of a
  GMRF, so environmental covariates and lagged or simultaneous paths drive
  recruitment. `summary()` reports the SEM path coefficients and the recruitment
  SD; the specification round-trips through `save_config()` / `load_config()`.
  `dsem = NULL` (the default) is unchanged and bit-identical to previous
  versions.

  **The whole diagnostic suite runs on a DSEM.** `retrospective()`, `jitter()`,
  `self_test()`, `profile()`, `process_residuals()`, `sim_mod()`,
  `osa_residuals()`, `remove_F()` and `reweight_comps()` all work, as does
  `sample_rec()` where the latent states span the projection. Each is described
  below where it needed something specific. `run_mse()` is the one exception: it
  stops with a message rather than project a model without the recruitment
  structure being tested.

* **`check_dsem_spec()` screens a specification before it is fitted**,
  reporting sparse predictors, a collinear lag-aligned design, covariates
  spanning orders of magnitude in SD, and any variable with no exogenous
  variance of its own. That last one matters more than it reads: a variable
  whose lag-0 two-headed path is pinned at zero (`BT <-> BT, 0, NA, 0`) is
  deterministic, so the model computes it from its incoming paths and the
  `env_data` series supplied for it is never read -- scaling that column leaves
  the objective bit-identical. It also drops out of the reported precision, so
  `sample_rec()` cannot draw a projection. Omitting the two-headed line entirely
  is a different thing and is fine: `dsem` then adds a freely estimated
  variance. `fit_mod()` does not run this check for you.

* **`sample_rec()` draws a DSEM's projected recruitment from the fitted
  process.** Recruitment deviations under a DSEM are not independent -- a
  self-path makes them autocorrelated, covariate paths make them respond to the
  environment, and a cross-species path couples the stocks -- so resampling the
  hindcast projects a different process from the one estimated. The projection
  is now drawn from the GMRF that scored the fit, conditional on everything the
  fit already knows, as `dsem::simulate()` does. All columns are drawn jointly,
  so cross-species and shared-covariate paths carry through.

  A covariate supplied over the projection is a known future value, not
  something to draw. Under `family = "fixed"` those states are pinned by the
  map, and the draw conditions on them, so projected recruitment responds to the
  climate scenario as specified.

  `sample_rec = FALSE` gives the conditional mean, matching what the argument
  means without a DSEM: projected recruitment at `R0 * exp(mean deviation)`,
  not at the lognormal median.

  This needs `build_DSEM(estimate_projection = TRUE)`, so the latent states
  span the years being drawn; with the default `FALSE` they stop at `endyr`,
  and a projection draw would extrapolate past the likelihood's support.
  `sample_rec()` says so rather than extrapolating.

* **`build_DSEM(estimate_projection = TRUE)` now overrides
  `build_srr(proj_mean_rec = )`.** The two are contradictory -- under
  mean-recruitment projection the model projects `avg_R` past `endyr`, so the
  SEM's projection states are fitted and then reach only the dynamic B0/BF
  reference series -- and `proj_mean_rec` defaults to `TRUE`, so asking for a
  SEM-projected recruitment was inert unless you also changed the stock-recruit
  setup. Estimating the projection now wins, and `fit_mod(verbose > 0)` reports
  it. Set `estimate_projection = FALSE` to keep mean-recruitment projection.

* **Lognormal bias correction is now applied to a DSEM's recruitment
  deviations, so `R0` is comparable to a non-DSEM fit again.** It was not
  applied before, and while SSB absorbed the offset `R0` did not: measured on
  BS2017SS with a naive sem, `R0` came out 20-51% below the non-DSEM fit while
  terminal SSB agreed to 0.4%, and dynamic B0 and the Tier-3 B40% proxy are
  keyed to `R0`.

  The correction is `-Var/2` using the variance of each recruitment state,
  taken from the GMRF and reported as `dsem_margvar_tj`. It is NOT `-R_sd^2/2`:
  `R_sd` is the GMRF's conditional (innovation) SD, and the two agree only for a
  state with no incoming lagged path -- for a first-order self path they differ
  by `1/(1-rho^2)`, so the naive form would under-correct exactly the lagged
  models a DSEM exists to fit.

  The variance is conditional on the environmental values the model was given,
  which are data rather than random; counting their prior variance instead
  inflated it by `beta^2 Var(covariate)/(1-rho^2)`, measured at +67%. **What
  counts as given is decided per year, not per series**, because `env_data` is
  padded to the model's full span: a covariate that starts after `styr`, has a
  gap, or carries no projection scenario is a free latent state in exactly those
  years. Treating the whole series as known under-corrects there -- measured on
  BS2017SS with a hindcast-only covariate, projected recruitment 17.1% high.
  Verified in `tools/verify/verify-dsem-recovery.R` (matches `sigma^2/(1-rho^2)`
  to +0.00%, collapses to `sigma^2` at `rho = 0`) and in
  `test-dsem-bias-correction.R`, which checks it against
  `diag(solve(Q[unknown, unknown]))` from the model's own reported precision for
  covariates observed everywhere, observed partly, and autoregressive. A naive
  DSEM reproduces a non-DSEM fit's objective, `R_sd` and `R0` at the default
  `bias_adjust_proc`.

* **`retrospective()`, `jitter()` and `self_test()` run on a DSEM.** A peel does
  not shorten the model -- it turns off data after `endyr_peel` -- so the latent
  states stay free and in `random`, and the Laplace approximation integrates the
  peeled-year states out against the GMRF prior. That is the peeled likelihood.
  Pinning them the way `rec_dev` peels would be wrong: pinning is inert for an
  independent deviation, but DSEM states are coupled, so pinned zeros stay in
  the quadratic form, shrink the terminal retained state, and bias Mohn's rho.
  Checked against a model that genuinely ends at `endyr_peel`: terminal
  recruitment agrees to 1.2e-4, while the peel differs from the unpeeled fit by
  19%.

  `retrospective(rescale = TRUE)` is the better setting for a DSEM. It
  standardizes `env_data` on `styr:endyr_peel` and withholds the covariate
  values after it, so the peel does not centre or condition on years it is
  supposed to have not seen.

  `jitter()` perturbs the SEM path coefficients and the latent recruitment
  deviations, but not the covariate columns of the latent state matrix: with
  `family = "fixed"` those are the environmental data, so jittering them would
  perturb the data rather than the starting values.

  A `self_test()` of a DSEM re-realizes the recruitment process when asked to:
  `process = TRUE` draws the latent field from the precision that scored the fit
  (see the `sim_mod()` entry above), while the default `process = FALSE` leaves
  the states at their estimates, as it does without a DSEM, so a clean self test
  says the model recovers itself under fresh observation error -- not that the
  SEM structure is identified across recruitment realizations.

* **`summary()` reports the SEM path coefficients and the recruitment SD** for a
  DSEM fit: `coefficients` (one row per path, with `Estimate`, `Std_Error`,
  `z_value`, `p_value`) and `recruitment_sd` (per species, with whether the SD
  was estimated or fixed). The three uncertainty columns are always present and
  `NA` without an `sdreport`, so tables from several fits can be row-bound. Note
  that a two-headed (`<->`) variance row's sign is not identified and its Wald
  p-value is a boundary test, so neither is interpretable in the usual way. Fits
  without a DSEM are unaffected.

## Breaking changes

* **`env_data` is now read as numbers.** It is coerced to a data.frame of
  numeric columns with an integer `Year` on the way in -- in `read_data()`,
  `clean_data()` and `rearrange_data()` -- so a table supplied as a matrix, or
  carrying a column that arrived as text (a workbook column formatted as text,
  a factor, a character matrix), is usable rather than fatal. A column of text
  previously survived to `MakeADFun()`, which stopped with "Only numeric
  matrices, vectors, arrays, factors, lists or length-1-characters can be
  interfaced" after the fill of missing years had already given up on it
  ("argument is not numeric or logical: returning NA") and left it `NA` for
  every year. A cell that is genuinely not a number now names its column and
  its value instead of becoming an `NA` that the fill would replace with a
  column mean -- a fabricated observation. Blank cells and the literal `NA` /
  `NaN` still mean "no value" and are still filled.

## Dependencies

* **`dsem` is a new Suggests, pinned to `== 3.0.0`.** `src/TMB/dsem.hpp` is
  vendored from that version and `build_dsem_objects()` harvests
  `dsem::dsem(run_model = FALSE)`'s TMB inputs directly, so a different `dsem`
  can emit a data layout the compiled template cannot read. The pin is exact
  rather than a minimum for that reason, and `build_dsem_objects()` checks it
  and says so. Only DSEM models need it; every other fit is unaffected.

* `utils` is now imported (`utils::head`, `utils::packageVersion`).

**This changes results** for a DSEM simulation: `sim_mod()`, `self_test()` and
anything built on them draw differently from 5.19.0, and a covariate with a
measurement density now varies across replicates. Non-DSEM models are
unaffected -- every new draw is gated on `dsem_on == 1`, and the four golden
reference fits are unchanged.

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
