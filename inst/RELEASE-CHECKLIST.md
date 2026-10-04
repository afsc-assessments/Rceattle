# Pre-release checklist

Steps to run *before* tagging a new version of Rceattle. This is a
research / management-facing package, so the priority is reproducibility
and clear release notes — not speed.

## 1. Source state

- [ ] All in-progress changes committed; working tree clean
      (`git status`).
- [ ] `DESCRIPTION` `Version:` field bumped (semver: bump patch for
      bug fixes / docs, minor for new features, major for breaking
      API changes). "Breaking" means no back-compat path. A removal with a
      deprecation message that keeps old fits working is a minor bump:
      `growth_re` was removed with a `switch_check()` message and a
      `fit_mod()` guard dropping retired blocks from `inits`, and shipped as
      a minor. Refusing a configuration that never fitted the model it
      described (inert, self-contradictory, double-counted, reading values
      not yet computed, penalizing years outside the hindcast) is also minor,
      though stored fits with it stop refitting; list each under
      `## Breaking changes` with the rebuild (5.33.0, 5.35.0).
- [ ] `NEWS.md` top section heading matches the new version. Convert any
      "Unreleased" placeholder to the version number. Headings carry the
      version alone, with no date, so that one entry can cite another as
      `# Rceattle X.Y.Z` without going stale.

## 2. Verify

```r
# Regenerate man/ and NAMESPACE from roxygen comments
roxygen2::roxygenise()

# Locally
devtools::document()
devtools::test()

# Vignettes (rendering, not just chunk parsing)
devtools::build_vignettes()

# CRAN-equivalent check
devtools::check(args = c("--no-manual", "--as-cran"))

# URL rot
urlchecker::url_check()

# Spell check (DESCRIPTION + man/)
devtools::spell_check()
```

All four should be 0 errors / 0 warnings. The `installed size`
NOTE is acceptable until we move large `data/` objects to
`inst/extdata/`.

## 3. Tag

Tags are the bare version, **no `v` prefix**: `4.4.0`, `4.8.0`, `5.8.1`,
`5.20.0`. Only the very first tag, `v4.3.0`, carried one. Consumers pin
against these, so the shape has to stay put.

```bash
# Replace X.Y.Z with the DESCRIPTION version. Tag the MERGE COMMIT on
# main -- the commit CI actually validated -- not a local branch head.
git fetch origin main
git tag -a X.Y.Z origin/main -m "Rceattle X.Y.Z"
git push origin X.Y.Z
```

**Read `DESCRIPTION` on the merge commit and tag that version, not the one you
remember.** While a release PR is open, anything merged into `dev` lands in the release
and moves the target: #158 went 5.41.0 → 5.42.0 → 5.42.1 while under review. A tag
naming the wrong version cannot be quietly fixed, because consumers pin against it.
`git show origin/main:DESCRIPTION | grep '^Version:'` is the check.

The `pkgdown` GitHub Actions workflow rebuilds the website on the
`release` event, so a GitHub Release must be published from the tag —
drafting one is not enough, the event is `release: published`.

**Check that a release actually rebuilt the site rather than assuming the
event fired.** It has silently not fired: 5.20.0 has a `release`-triggered
pkgdown run and 5.21.0, published the same way, got none. The recovery is
a manual dispatch, `gh workflow run pkgdown.yaml --ref main`.

Write the body; do not paste `NEWS.md`. A release that folds a dozen
versions spans well over a thousand lines there, and the reader needs the
short answer: what forces a refit, what breaks, what is new. Follow the
5.8.1 and 5.20.0 bodies — a "results change" table of change against
effect, then breaking changes, then new features — and link `NEWS.md` for
the detail.

```bash
gh release create X.Y.Z --title "Rceattle X.Y.Z" --notes-file notes.md --latest
```

## 3b. Dispatch the deep checks

`deep-checks` runs the golden regression and the bounds-checked build, and
**nothing else does** — it is `skip_on_cran()` and `skip_on_covr()`, so it
fires in neither `R-CMD-check` nor `test-coverage`. GitHub reads `schedule`
and `workflow_dispatch` from the default branch only, so a release is the
first moment it can run against the new code:

```bash
gh workflow run deep-checks.yaml --ref main
```

## 4. Confirm reproducible install

From a clean R session on a different machine (or in a `renv` sandbox):

**This step passes on a stale install if you let it.** `withr::with_temp_libpaths()`
*prepends* the temporary library, it does not isolate — the working library stays on
`.libPaths()`. So if the install fails, `library(Rceattle)` loads whatever version is
already installed and `packageVersion()` reports that, which is how this check once
returned `5.33.0` while claiming to verify a 5.42.1 branch. Install into an explicit
`lib`, read the version back **out of that directory**, and do not pass `quiet = TRUE`,
which hides the failure this is meant to catch.

```r
# An explicit temporary library, so verifying a release neither overwrites the
# working install mid-assessment nor silently reads it instead.
lib <- file.path(tempdir(), "relcheck"); dir.create(lib, showWarnings = FALSE)
.libPaths(c(lib, .Library))          # temp library + base/recommended ONLY
remotes::install_github("afsc-assessments/Rceattle", ref = "X.Y.Z",
                        lib = lib, upgrade = "never", force = TRUE)
# Read the version from the installed DESCRIPTION, not from a loaded namespace:
# that is the assertion, and it cannot be satisfied by another copy.
read.dcf(file.path(lib, "Rceattle", "DESCRIPTION"))[1, "Version"]   # must be X.Y.Z
library(Rceattle, lib.loc = lib)
dirname(getNamespaceInfo("Rceattle", "path"))   # must be `lib`, not the user library
citation("Rceattle")                            # should list all references
example("fit_mod", package = "Rceattle", give.lines = TRUE)
```

`ref = "X.Y.Z"` rather than `"...@X.Y.Z"` so a branch name containing a `/` is passed
through untouched.

## 5. Operational consumers

For users running Rceattle inside an assessment / MSE pipeline, communicate
that they should pin to the new tag in their lockfile (`renv.lock`) or
`DESCRIPTION` `Remotes:` field rather than tracking `main`. The
maintainer email in `DESCRIPTION` is the formal contact.

`../Rceattle-models`, `../GOA-ATF-ESP` and `../GOA-multispecies-assessment`
consume this API directly; see `inst/dev/SIBLING-REPOS.md`.

## 6. (Optional) CRAN submission

When ready:

```r
devtools::release()
```

This reads `cran-comments.md` and uploads to CRAN via the official
submission form.
