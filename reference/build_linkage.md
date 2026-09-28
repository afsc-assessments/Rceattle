# Formula-driven linkage specifications for Rceattle processes

These helpers let users describe how a process parameter (e.g. `alpha`
for the Beverton-Holt SRR, `M1` for natural mortality, `K` for von
Bertalanffy growth) depends on environmental covariates and on
stratifying factors (species, sex). They produce:

## Details

1.  an `Rceattle_linkage_spec` object that captures the user's intent
    (formula + grouping) without committing to a global column index,
    and

2.  a `materialize_linkage()` step that, given the env data and stratum
    levels, expands the spec into the canonical long-format
    linkage-table rows consumed by TMB.

Splitting capture from materialization lets
[`fit_mod()`](https://afsc-assessments.github.io/Rceattle/reference/fit_mod.md)
pool specs from every process, build a single shared design matrix, and
assign globally consistent `X_col` indices in one place.
