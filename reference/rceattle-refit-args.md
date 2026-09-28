# Shared refit-diagnostic arguments

A common vocabulary, documented here once and inherited by the
diagnostics that refit the model:
[`retrospective()`](https://afsc-assessments.github.io/Rceattle/reference/retrospective.md),
[`jitter()`](https://afsc-assessments.github.io/Rceattle/reference/jitter.md)
and
[`self_test()`](https://afsc-assessments.github.io/Rceattle/reference/self_test.md).
Each argument means the same thing wherever it appears. `phase`, `getsd`
and `timeout` are deliberately **not** here, their defaults and their
consequences for the diagnostic differ, so each function documents its
own. [`profile()`](https://rdrr.io/r/stats/profile.html) refits too but
takes only `cores` and a bare `getsd`.

## Arguments

- object:

  an Rceattle model fit using
  [`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md).

- Rceattle:

  deprecated name for `object`, still accepted so existing scripts keep
  working. Supplying both is an error.

- cores:

  Number of cores for the parallel refits. Default `NULL` picks
  `parallel::detectCores() - 6`, capped at 2 under `R CMD check` (which
  sets `_R_CHECK_LIMIT_CORES_`). Set to 1 to force sequential execution.

- fit_control:

  optional
  [`fit_control()`](https://afsc-assessments.github.io/Rceattle/reference/fit_control.md)
  bundle for the refits. Only `phase` and `getsd` are read; see **What
  [`fit_control()`](https://afsc-assessments.github.io/Rceattle/reference/fit_control.md)
  reaches**.

## What [`fit_control()`](https://afsc-assessments.github.io/Rceattle/reference/fit_control.md) reaches

These diagnostics refit through an internal helper that forwards `phase`
and `getsd` and nothing else. Every other field either comes off the
fitted model's own `data_list` (`loopnum`, `comp_offset`, the
bias-adjustment flags) or takes
[`fit_control()`](https://afsc-assessments.github.io/Rceattle/reference/fit_control.md)'s
default. Setting one of those is an error rather than a silent no-op, so
a field you set is a field the refit used.

Only the fields you set are applied, named in the
[`fit_control()`](https://afsc-assessments.github.io/Rceattle/reference/fit_control.md)
call, or assigned to afterwards (`ctl$getsd <- TRUE`). One you never
touch keeps the diagnostic's own default, which is not always
[`fit_control()`](https://afsc-assessments.github.io/Rceattle/reference/fit_control.md)'s,
[`retrospective()`](https://afsc-assessments.github.io/Rceattle/reference/retrospective.md)
phases its peels where
[`fit_control()`](https://afsc-assessments.github.io/Rceattle/reference/fit_control.md)
does not, so `fit_control(getsd = FALSE)` asks about standard errors and
changes nothing else. Setting a field to the value that is already
[`fit_control()`](https://afsc-assessments.github.io/Rceattle/reference/fit_control.md)'s
default still counts as setting it.
