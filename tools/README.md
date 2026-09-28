# tools/

`verify/` holds ad-hoc verification harnesses, run by hand. `ci/` holds scripts a
GitHub Actions workflow calls.

The harnesses are not part of `test_check()`: they run real fits, take minutes to
hours, and answer questions a unit test cannot. Run them by hand
(`Rscript tools/verify/<script>.R`) after a change to the simulation, the
likelihood, or the MSE loop.

Each script prints a pass/fail digest and leaves its numbers on stdout.

## verify/

**Equivalence — must stay bit-identical**

| Script | Proves |
|---|---|
| `verify-refit-like.R` | The `.refit_like()` refit path is unchanged across all nine entry points. Golden covers none of them. |
| `verify-mse-hindcast-invariant.R` | An MSE does not perturb the hindcast, independent of the refit refactor. |
| `verify-mse-schedule.R` | An explicit vector of assessment years reproduces the equivalent period, a year-indexed `catch_mult` moves only the years it names, and two schedules stay on common random numbers. |
| `verify-mse-repro.R` | Full-MSE reproduction digest, single- and multispecies. |

**Simulation — is the draw the distribution the likelihood assumes?**

| Script | Proves |
|---|---|
| `verify-sim-centering.R` | Simulated observations are centred where the likelihood expects. An sd check alone passes a mis-centred draw. |
| `verify-sim-index-families.R` | Each survey is drawn from its own `Index_distribution`, not from lognormal for everything. |
| `verify-sim-diet.R` | Stomach contents are drawn from the distribution the diet likelihood assumes. |
| `verify-sim-linkage-re.R` | Recruitment, M and growth written as a random linkage are drawn from the covariance they were given. |
| `verify-sim-process-error.R` | The process draws, and the two `M1_re` defects fixed alongside them. |

**Recovery — does a refit get back what was simulated?**

| Script | Proves |
|---|---|
| `verify-sim-recovery.R` | Which observation type costs a `self_test()` its recovery, starting from deterministic recovery. |
| `verify-sim-recovery-process.R` | Simulate a process, refit, and see what comes back. |
| `verify-sim-recovery-M.R` | Whether the data actually inform M. Reproduces the documented limit: deviations recover, their SD does not. |
| `verify-sim-local-optimum.R` | When recovery is poor, whether that basin is the MLE for the simulated data or a local optimum. |
| `verify-sim-recovery-apical.R` | The per-sex apical selectivity offset recovers from joint compositions (`GOAatf` fishery, known male offset), with the mean, empirical SD and reported SE over `n` replicates. |
| `verify-sim-recovery-np-integrable.R` | Under `NonParametricIntegrable` with `random_sel = TRUE`, the reported selectivity-deviation SD recovers a known one (deviates drawn in R on `Atka2022`'s fishery, observations from `sim_mod()`). |
| `verify-osa-cdf.R` | `osa_residuals(method = "cdf")` returns iid N(0, 1) residuals when observations simulated from `BS2017SS` are residualized at the parameters that generated them; the Gaussian default is scored alongside. |
| `verify-osa-cdf-accuracy.R` | The two `"cdf"` trade-offs `?osa_residuals` states: its error under random effects against exact Kalman innovations, and what the 8.04 censoring costs a reported statistic. |
| `verify-sim-recovery-srr-msm.R` | Whether the data can inform a Beverton-Holt curve under predation, where no equilibrium anchors it: a curve bending inside a 7x SSB range is imposed on the two-species fixture and refitted after redrawing observations and recruitment. |

## ci/

| Script | Does |
|---|---|
| `coverage-shard.R` | Runs one weight-balanced slice of `tests/testthat` under `covr` and writes `cobertura.xml`. Called once per shard by `.github/workflows/test-coverage.yaml`. |

Coverage has to run the full `NOT_CRAN` suite serially (covr's line counters are
per-process, so anything run in a testthat parallel worker goes uncounted), and
the suite outgrew a single runner's 120-minute budget. `coverage-shard.R` splits
it across four runners and Codecov merges the uploads. Preview a split without
running any tests (the only mode that does not need `covr`):

```sh
COV_PLAN_ONLY=1 COV_SHARD=1 COV_SHARDS=4 Rscript tools/ci/coverage-shard.R
```

Balancing uses a static cost proxy read off each file (an optimization outweighs
a `MakeADFun` build; a `.refit_like()` caller outweighs an optimization), not
measured seconds. If the shards drift out of balance, the workflow's per-file
elapsed times are the input for retuning the weights.
