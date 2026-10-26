# Linkage table: long-format coefficient registry for Rceattle processes

A *linkage table* lists every coefficient (`beta`) that ties a process
(recruitment, natural mortality, growth, ...) to a covariate, time block
or stratum. Each row is one coefficient and records the process and
parameter it changes, the species/sex/age/fleet it applies to, the
column of the design matrix `X` it multiplies, and its link, prior,
bounds and phase.

## Details

One row per coefficient, rather than one array per process, means the
table grows only with the linkages requested. The user-facing
constructors are in `R/0-build_linkage.R`.
