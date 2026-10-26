---
description: Recompile the TMB C++ model and reload Rceattle
---

Recompile the TMB model and reload the package, then report the result cleanly.

Run from the repo root:

```
export PATH=/usr/bin:$PATH && Rscript -e 'pkgload::load_all(".", quiet = TRUE)'
```

Notes:
- The `export PATH=/usr/bin:$PATH` prefix puts the system toolchain first — this repo's
  sessions do it consistently so the TMB DLL compiles cleanly (avoids a Homebrew
  clang/gfortran shadowing the build). Keep it on every compile/check command.
- `load_all()` recompiles `src/TMB/*.cpp` (slow, ~60-140 s).
- **A `.hpp`-ONLY edit can leave it loading the OLD binary, silently** — and every
  process lives in a header (`selectivity.hpp`, `predation.hpp`,
  `recruitment.hpp`, `growth.hpp`, `linkage.hpp`, `spr.hpp`), so this is the
  normal case, not an edge. Seen again on 2026-10-21: a `selectivity.hpp` fix
  read as failing its own new test 4 times over two runs, and passed 4/4 once
  the rebuild was forced. **Do not reason about why** — `inst/dev/TRAPS.md`,
  "Build and test", records that the mechanism is not established.
  **Use its recipe:** delete `src/TMB/ceattle.o`, `src/TMB/ceattle.so`,
  `src/ceattle.so` and `src/Rceattle.so`, reload, then confirm
  `file.info("src/TMB/ceattle.o")$mtime` is newer than the header you edited.
  **Check the `.o`, never an `.so`** — `compile.R` copies the `.so` up
  unconditionally, so its mtime moves whether or not anything compiled. A reload
  that finishes in seconds is proof nothing built.
- If `$ARGUMENTS` contains `r-only` or `no-compile`, skip the rebuild by passing
  `compile = FALSE` to `load_all()` (use this when you only changed `R/`).
- If compilation fails, surface the compiler error with file:line and STOP — do not
  continue as if it built.
