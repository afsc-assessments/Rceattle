#!/usr/bin/env bash
# What machine produced this number?
#
# Two jobs here pin numeric results and can go red for reasons that have nothing
# to do with the change under test: `golden` compares four objectives against
# literals, and `test-coverage` runs regression tests that pin fitted values.
# Both have done so -- `NEWS.md` (5.45.1) records `deep-checks` red on `main`
# for weeks over a summation order, and `inst/dev/TRAPS.md` the 5.45.0 red, a
# gradient of 3.5e-4 on ubuntu against 2.5e-5 on macOS.
#
# The prompting case (one linkage fit, three different ubuntu objectives on
# byte-identical code), the measured platform table and why libm rather than
# BLAS is the candidate are all in `inst/dev/TRAPS.md`, under "BLAS is the LEAST
# likely culprit". None of the data needed to attribute such a disagreement was
# being recorded: the runner's own logs carry only "Platform: x86_64-pc-linux-gnu".
#
# So this prints the things that plausibly move a last bit, cheaply (~1 s), into
# every log that pins a number:
#
#   * CPU model and the FMA/AVX flags. glibc's libm dispatches exp/log/pow by
#     CPU feature through ifunc resolvers, and this model has 246 `exp()` and
#     105 `log()` call sites across `src/TMB/`, many of them inside the
#     species-age-year loops. The runner fleet is heterogeneous.
#   * glibc, which is that libm.
#   * the compiler, which chooses whether to contract a*b+c into an FMA.
#   * R's BLAS/LAPACK, for completeness. Note this model does NOT route Eigen
#     through them (no EIGEN_USE_BLAS, and one BLAS-level operation in
#     `ceattle.cpp`), so they are the least likely culprit here -- recorded so
#     the next person does not have to re-establish that, as happened once.
#   * any threading override, since a multithreaded BLAS reduces in a
#     nondeterministic order.
#
# It is a REPORT, not a check: nothing here fails a job. The point is that a
# cross-platform disagreement becomes attributable instead of a mystery.
set -u

echo "=== machine fingerprint ==="
echo "  runner image : ${ImageOS:-unknown} ${ImageVersion:-unknown}"
echo "  arch         : $(uname -m)  $(uname -s)"

if command -v lscpu >/dev/null 2>&1; then
  echo "  cpu          : $(lscpu | sed -n 's/^Model name: *//p' | head -1)"
  echo "  cpu id       : family $(lscpu | sed -n 's/^CPU family: *//p' | head -1)"\
" model $(lscpu | sed -n 's/^Model: *//p' | head -1)"\
" stepping $(lscpu | sed -n 's/^Stepping: *//p' | head -1)"
elif [ "$(uname -s)" = "Darwin" ]; then
  echo "  cpu          : $(sysctl -n machdep.cpu.brand_string 2>/dev/null || echo unknown)"
fi

# The flags that change which libm kernel an ifunc resolver picks.
if [ -r /proc/cpuinfo ]; then
  flags=""
  for f in fma avx avx2 avx512f avx512dq; do
    if grep -qw "$f" /proc/cpuinfo 2>/dev/null; then flags="$flags $f"; fi
  done
  echo "  cpu features :${flags:- none detected}"
fi

if command -v ldd >/dev/null 2>&1; then
  echo "  glibc (libm) : $(ldd --version 2>/dev/null | head -1)"
fi
echo "  compiler     : $(${CC:-cc} --version 2>/dev/null | head -1 || echo unknown)"

# Threading overrides, if a caller has pinned them.
thr=""
for v in OPENBLAS_NUM_THREADS OMP_NUM_THREADS MKL_NUM_THREADS; do
  eval "val=\${$v:-}"
  if [ -n "${val}" ]; then thr="$thr $v=$val"; fi
done
echo "  thread limits:${thr:- unset (BLAS picks its own)}"

# R-side: version, and which BLAS/LAPACK it actually linked.
if command -v Rscript >/dev/null 2>&1; then
  Rscript -e 'si <- sessionInfo()
    cat("  R            :", R.version.string, "\n")
    cat("  R platform   :", R.version$platform, "\n")
    cat("  BLAS         :", if (is.null(si$BLAS)) "unknown" else si$BLAS, "\n")
    cat("  LAPACK       :", if (is.null(si$LAPACK)) "unknown" else si$LAPACK, "\n")
    for (p in c("TMB", "Matrix", "TMBhelper")) {
      v <- tryCatch(as.character(utils::packageVersion(p)), error = function(e) "not installed")
      cat(sprintf("  %-13s %s\n", p, v))
    }' 2>/dev/null || echo "  R            : present but sessionInfo failed"
fi
echo "=== end fingerprint ==="
