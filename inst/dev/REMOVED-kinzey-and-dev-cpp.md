# Removed: the Kinzey & Punt scaffolding and `src/TMB/Dev/`

State, not policy. Written when the commented-out code was deleted, so the specification and the
reason survive the deletion. Nothing here is live; this file is the record of what was taken out and
where to find it if it is ever revived.

## 1. The Kinzey & Punt (2009) functional-response parameters

**What they were.** Six parameter blocks for the Kinzey & Punt predation functional responses
(`msmMode` 3–9): `logH_1`, `logH_1a`, `logH_1b`, `logH_2`, `logH_3`, `H_4`.

**Specification.** Kinzey, D. and Punt, A.E. (2009). *Multispecies virtual population analysis and
statistical catch-at-age models: a comparison.* The MSVPA forms that Rceattle does implement are
Holsman et al. (2015) (`msmMode = 1`) and Holling Type III (`msmMode = 2`).

**Why they were commented out, not deleted, before now.** They were commented out in the TMB
template *and* on the R side, consistently, and `data_check()` refuses `msmMode` 3–9. `NEWS.md`
(see the 4.x entry on `plot_form()`) records that `plot_form()` errored because of it. So the blocks
were dead in every direction — the R side could not supply a parameter the template does not
declare, and the template does not declare it.

**What was deleted, and the starting values / bounds they carried** (the only content worth
keeping, since the structure is recoverable from the citation):

| Where | Content |
|---|---|
| `R/2-build_params.R` (was section 3.2) | `logH_1 = matrix(-8.5, nspp, nspp + 1)`; `logH_1a = rep(-3, nspp)`; `logH_1b = rep(0, nspp)`; `logH_2 = matrix(-9, nspp, nspp)`; `logH_3 = matrix(-9, nspp, nspp)`; `H_4 = matrix(1, nspp, nspp)` |
| `R/4-build_parameter_bounds.R` | `logH_3` bounded `[-30, -1e-06]`; `H_4` bounded `[-0.1, 20]` |
| `R/6-phaser.R` | all six at phase 6, annotated "not used in MSVPA functional form" |
| `R/3-build_map.R` (`build_map_predation`, was sections 2 and 3) | 56 lines of per-`msmMode` map logic — see "The map logic, stated carefully" below |
| `R/7-plot_ceattle.R` (`plot_form()`) | 82 lines computing and plotting the response surfaces from `H_1`…`H_4`, on a mode encoding that **disagreed** with the rest of the package — see below |
| `src/TMB/ceattle.cpp` | the `PARAMETER_*` declarations and the `H_1 = exp(logH_1.array())` derivations, inside `/* */` blocks. **Still present in the template** — see "Not done here" below |

### The map logic, stated carefully

This is the only part of the deleted code that is not recoverable from the citation, so read the
convention before the content: **in a TMB map, `NA` means NOT estimated.** The deleted lines were of
the form `map_list$logH_2 <- map_list$logH_2 * NA`, which **fixes** that block. An earlier draft of
this file said Ecosim "freed `logH_2` and `H_4`" — that is exactly backwards, and it was the one
non-mechanical fact the file existed to preserve.

Stated correctly, for `msmMode = 9` (Ecosim): `logH_2` and `H_4` were **fixed**; `logH_1`,
`logH_1a`, `logH_1b` and `logH_3` were left **free**. MSVPA (`msmMode` 1 and 2) fixed all six. The
intermediate Holling and interference forms each fixed a different subset; recover those from
`git log -p` on `R/3-build_map.R` rather than from prose, since the pattern is per-mode and there is
no rule to summarise.

### The two deleted mode tables disagreed with each other

Worth recording because it is a defect, not a design. `build_map_predation()` encoded:

> 3 = Holling Type 1 · 4 = Holling Type 2 · 5 = Holling Type 3 · 6 = predator interference ·
> 7 = predator preemption · 8 = Hassell-Varley · 9 = Ecosim

`plot_form()`'s own `switch(as.character(msmMode))` encoded the same seven forms **shifted down by
one** (2 = Holling Type I … 8 = Ecosim), and had no case for 9 at all, so `msmMode = 9` fell through
to its own `stop("msmMode not implemented")`. `R/1-data_check.R` and
`vignettes/model-options-and-functionality.Rmd` both agree with `build_map_predation()`, so
`plot_form()` was the wrong one — and since its default argument was `msmMode = 3`, calling it bare
would have plotted Holling Type II for the model's Holling Type I. **Any revival must pick one
encoding and make `plot_form()` follow `data_check()`, not the reverse.**

**`plot_form()` still exists and still errors**, with a message naming the available modes. That
error is behaviour, not history: it is what a user who asks for `msmMode` 3–9 sees. It was kept.

**To revive this:** uncomment the **three** `/* */` blocks listed under "Not done here" — two in
`ceattle.cpp` and the 216-line implementation in `predation.hpp` — then restore the five R-side
pieces above from this file and from git history, then lift the `data_check()` refusal of
`msmMode` 3–9. Two parts are not mechanical and both are written out above: the per-mode map logic,
and reconciling `plot_form()`'s mode encoding with `data_check()`'s.

## 2. `src/TMB/Dev/`

**Deleted:** `wham_v0.cpp` (2,715 lines), `caal.hpp` (468), `osa.hpp` (71),
`CAAL-simulation.R` (237), `CAAL simulation.R` (349).

**Proof it was dead, four independent ways** (re-verified at `b066015f`, 5.49.1):

1. `src/TMB/compile.R` compiles exactly one translation unit — `tmb_name <- "ceattle"`, then
   `TMB::compile(file = paste0(tmb_name, ".cpp"))`. No glob, no second `compile()` call.
2. `src/Makevars` and `src/Makevars.win` each have a single `tmblib` recipe,
   `(cd TMB; Rscript --no-save --no-restore compile.R)`. Nothing else descends into `TMB/`, and
   `clean` removes only build artefacts (`*.so *.o TMB/*.so TMB/*.o`).
3. `grep -rn 'Dev/' src/` returns nothing, and none of `ceattle.cpp`'s `#include` lines names it.
   (A substring search for `osa.hpp` appears to hit many files; every one of those is
   `comp_osa.hpp`, which is live and unrelated.)
4. `.Rbuildignore` carried `^src/TMB/Dev$`, so the directory never shipped in the tarball. (That
   line is deleted by the same commit, so do not look for it by number -- grep the pattern.)

**What was in it.** `wham_v0.cpp` was a vendored copy of WHAM's template, kept for the OSA-residual
cross-check. That comparison now lives in `tests/comparison/WHAM-OSA-comparison.R`, which does not
reference it. `caal.hpp` and `osa.hpp` were prototypes superseded by the shipped `comp_osa.hpp` and
the CAAL support in `ceattle.cpp`.

**Backlog effect.** `src/TMB/Dev/` carried **7** TODO/FIXME markers -- `caal.hpp` 5 and
`wham_v0.cpp` 2 -- not the 5 the backlog had attributed to `caal.hpp` alone. The backlog's total
and per-area breakdown were re-derived in the same commit.

## Not done here — three C++ remnants remain

All still commented, all still in the shipped template:

| File | Block | Size |
|---|---|---|
| `src/TMB/ceattle.cpp` | `/*` 551 – `*/` 558 | the six `PARAMETER_*` declarations |
| `src/TMB/ceattle.cpp` | `/*` 790 – `*/` 808 | the `H_1 = exp(logH_1.array())` derivations and the Kinzey arrays |
| **`src/TMB/predation.hpp`** | **`/*` 521 – `*/` 736** | **the whole functional-response implementation — `// 8.2. KINZEY PREDATION EQUATIONS`, `if (msmMode > 2)`, the `switch` over the forms, `H_2`/`H_3`/`H_4`. 216 lines, and it holds the 5 `predation.hpp` TODO/FIXME markers `CLEANUP_BACKLOG.md` still counts.** |

**Why they were left.** Deleting them is a C++ edit: it needs `/recompile` and a golden run, which
is a different kind of change from an R-side comment sweep. **It is not because of the
source-parsing guards** — an earlier draft of this file said so and was wrong.
`test-schema-jnll-rows.R:41-47` and `test-schema-cpp-dispatch.R` both strip `/* */` blocks before
scanning, so deleting commented code is invisible to them. (The line-anchored concern is real but
different: it applies to *reformatting* live code, e.g. clang-format putting two enumerators of the
`JnllRow` enum on one line, which would break the `^\s*JNLL_` scan.)

Doing all three in one C++ commit would take `predation.hpp` from 738 lines to ~522 and drop the
backlog's marker count by a further 5.
