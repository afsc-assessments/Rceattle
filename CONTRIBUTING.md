# Contributing to Rceattle

Rceattle fits CEATTLE, a single- and multispecies age-structured stock assessment. Its output
sets US federal catch limits. A wrong number here does not crash; it becomes a quota. Every
convention below exists so that a change cannot move a fit without someone noticing.

This guide is written for a fisheries scientist or ecologist who wants to fix, extend or
understand the package. `CLAUDE.md` holds the same rules as an operating manual for the coding
agent; where a rule is stated in full there, this guide points to it rather than repeating it.

## Setting up

You need R 4.1 or later, a C++ compiler and the packages in `DESCRIPTION`. The likelihood is a
TMB template under `src/TMB/`; it compiles the first time you load the package.

```r
install.packages(c("devtools", "pkgload", "TMB"))
pkgload::load_all(".")     # compiles src/TMB/ceattle.cpp (a few minutes), then loads R/
```

On macOS put the system toolchain first, or a Homebrew clang or gfortran can shadow the one
TMB expects:

```sh
export PATH=/usr/bin:$PATH
```

If clang refuses with "You have not agreed to the Xcode license", set
`DEVELOPER_DIR=/Library/Developer/CommandLineTools` in the same shell. Compiled artifacts
(`*.o`, `*.so`, `*.dll`) are gitignored. Never commit them.

Edits to R code are picked up by `pkgload::load_all(".", compile = FALSE)`. Edits to
`src/TMB/*.cpp` or `*.hpp` are inert until `pkgload::load_all(".")` recompiles.

## Running the tests

```sh
NOT_CRAN=true TESTTHAT_PARALLEL=false Rscript -e 'devtools::test()'
```

`NOT_CRAN=true` turns on the tests that fit real models; without it most of the suite is
skipped. `TESTTHAT_PARALLEL=false` is required after a C++ change, because the parallel workers
cannot load a freshly rebuilt DLL.

To run one file with the shared helpers loaded:

```r
pkgload::load_all(".", quiet = TRUE)
e <- new.env(parent = asNamespace("Rceattle"))   # plain new.env() cannot see internal helpers
for (f in list.files("tests/testthat", "^helper", full.names = TRUE)) sys.source(f, envir = e)
testthat::test_file("tests/testthat/test-selectivity-logisticpm.R", env = e)
```

Three checks sit outside the ordinary suite:

- **The golden regression** (`tests/testthat/test-golden-regression.R`) fits four reference
  models and pins their objectives. Any change that can move a fit must leave them unchanged,
  or re-pin them with the numbers in `NEWS.md`. It runs only under `NOT_CRAN=true`.
- **`tools/verify/`** holds harnesses for what the golden fits cannot see: refit paths,
  simulation draws, MSE reproducibility. `tools/README.md` says which to run for which change.
- **`inst/dev/TRAPS.md`** lists the ways a change has silently gone wrong before, with the
  measured numbers. Read it before touching anything it names.

## Where to read first

- `vignettes/articles/developer-guide.Rmd`: the fit pipeline, the switch system, the column
  schema and the linkage grammar.
- `vignettes/articles/adding-a-selectivity-form.Rmd`: one extension traced end to end, with
  the files it touches and the tests that catch a half-finished job.
- `R/0-column_schema.R`: the source of truth for every workbook column and switch value.
- The C++ reference on the package site (built from the Doxygen comments in `src/TMB/`):
  the model equations, function by function. `src/TMB/recruitment.hpp` is the header to
  emulate when you document C++.

## Making a change

**Branches.** Cut a branch from `dev` and open a pull request into `dev`. The maintainer
merges `dev` into `main` as a release, following `inst/RELEASE-CHECKLIST.md`. Operational
assessments pin a tagged release, so `main` moves only at a release.

**What a pull request owes.** A behaviour, API or documentation change updates
`NEWS.md`, the `DESCRIPTION` version and the affected vignette together (what counts as
breaking, and the `_pkgdown.yml` and `man/` obligations, are rules 5 and 6 in `CLAUDE.md`).
Add a test named `tests/testthat/test-<area>-<topic>.R`; a test that runs a real
optimization starts with `testthat::skip_on_cran()`.

**Commit messages.** Plain text, imperative subject, a body that says why; rule 13 in
`CLAUDE.md` is the full statement.

**Continuous integration.** A pull request into `dev` runs `R CMD check` on macOS, Windows and
three Ubuntu R versions, plus six test-coverage shards. The package site (`pkgdown.yaml`)
builds only for `main`, so check a documentation change locally with
`pkgdown::build_reference_index(pkgdown::as_pkgdown("."))` and, for an article,
`pkgdown::build_article("<name>")`.

## The rules that constrain a change

Stated in full under "Hard rules" in `CLAUDE.md`. In brief:

1. Preserve the public API. Deprecate an argument; do not delete it.
2. A change that can move a fit needs the golden regression, and the `tools/verify/` harness
   that covers what golden cannot.
3. The column schema defines every switch value, default and column order. Read them from
   it; never hardcode one elsewhere.
4. Do not invent a switch code, a default or a unit. If it is not in the schema or a
   vignette, ask.
5. `fleet_control$Fleet_code` equals the row number, and selectivity bin columns follow two
   conventions (ordinal versus absolute age). `nages` is a count of bins, not the oldest age.
6. `src/TMB/linkage.hpp` and `R/0-linkage_encode.R` hold the same process and parameter
   codes; change one, change both.
7. The repositories in `inst/dev/SIBLING-REPOS.md` consume this API. Sweep them after a
   breaking change.

## Writing comments and documentation

Write for a fisheries scientist who was not in the room: the assessment reason, the units
and the convention, in one or two lines. The "Comments" and "Domain vocabulary" sections of
`CLAUDE.md` give the rules and a before-and-after example; `src/TMB/recruitment.hpp` is the
C++ header to emulate.

## Getting help

Open an issue at <https://github.com/afsc-assessments/Rceattle/issues> so the answer is
searchable. The maintainer address in `DESCRIPTION` is the fallback.
