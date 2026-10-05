# Run the tests that read the package's own SOURCE, on a pull request.
#
# Sixteen test files assert that two hand-synced copies of something agree: a
# switch map against the `case` labels in the template, the column schema
# against `R/data.R`'s field list, the `JnllRow` enum against its two R-side
# partners, the HCR-2 threshold against the literal source of `run_mse()`, the
# developer guide against the files and C++ switches it names. They do that by
# reading `R/*.R` and `src/TMB/*.cpp` off disk.
#
# WHERE THEY RUN TODAY, precisely -- two earlier versions of this comment got
# the mechanism wrong, so it is spelled out.
#
#   * `R CMD check` copies the tests to `<pkg>.Rcheck/tests/` and runs them with
#     the working directory at `<pkg>.Rcheck/tests/testthat/`. Every guard here
#     builds its path as `test_path("..", "..", "R", ...)` or
#     `test_path("..", "..", "src", "TMB")`, which resolves to
#     `<pkg>.Rcheck/R` and `<pkg>.Rcheck/src/TMB`. NEITHER DIRECTORY EXISTS.
#     (The lazy-load database is a level deeper, at
#     `<pkg>.Rcheck/Rceattle/R/Rceattle.rdb`, and is never on the path these
#     tests construct -- so "R/ holds Rceattle.rdb" is not the reason.) Each
#     guard `skip_if()`s, and a skip is not a failure, so the job is green.
#
#   * `test-coverage` DOES run them against the real source. `coverage-shard.R`
#     resolves `tests/testthat` to an absolute path in the source checkout, so
#     `../../R/data.R` exists. But it passes `stop_on_failure = FALSE` and no
#     later step reads the result, so a failure prints into a shard log and the
#     job stays green.
#
#   * `deep-checks.yaml` runs them properly -- `load_all()` on a source
#     checkout, NOT_CRAN=true, failures fatal -- but nightly, gating no merge.
#
# So the gap was never "these never run on a pull request". It is that where
# they run, their result is discarded or arrives a day late. This job runs them
# on a source checkout and makes the result fatal.
#
# IT FAILS ON THREE CONDITIONS, not one, because two of the three are the ways
# a guard dies quietly:
#   * a guard fails;
#   * a guard returns NO RESULT ROWS -- a file-level `skip_on_cran()`, at column
#     0 outside any `test_that()`, aborts the file and yields zero rows, so it
#     reports neither a pass nor a skip and vanishes from the counts. 94 of the
#     255 test files are written that way;
#   * a BLOCK that reads package source skipped, or asserted nothing.
#
# That last condition is per BLOCK, deliberately. Two earlier versions checked
# per file -- first by a `skip_on_cran()` heuristic, then by a per-file
# `passed > 0` floor -- and an adversarial review demonstrated the same hole in
# both: 10 of the 16 files hold several `test_that()` blocks, so one block can
# be switched off while its siblings keep the file's total well above zero.
# Adding a `skip_on_cran()` is the single most likely accidental edit in this
# repo (it is the idiom in 181 of 255 files), and doing it to the block that
# guards the C++ dispatch map took 33 assertions dark with the job still green.
#
# WHAT THE PER-BLOCK CHECK DOES NOT COVER, because the first version of it
# implied more than it delivered. A block is checked only when its own body
# matches a discovery pattern. Three files read source through a TOP-LEVEL
# helper instead -- `cpp_source()` in test-schema-jnll-rows.R is the pattern --
# so no block body matches and only the zero-row check applies to them. The
# per-file `blocks=` count is printed for exactly this reason: a 0 there means
# the file has no per-block cover, not that it has nothing to cover. Requiring
# that no block in such a file skips is not viable -- most of these files
# legitimately skip their fit blocks when NOT_CRAN=false.
#
# Discovery is per LINE. Applied to a whole file as one string, `.` matches a
# newline in R's default regex engine, which matched two files that read no
# package source at all.

root <- normalizePath(".")
test_dir <- file.path(root, "tests", "testthat")
stopifnot(dir.exists(test_dir))

# The guards need no fit. With NOT_CRAN unset or "true" the fit tests in these
# files attempt real optimizations against a DLL this job does not build, and
# the job goes red for the wrong reason -- so assert the mode rather than trust
# it. `inst/dev/TRAPS.md` records that this value "did NOT hold reliably" when
# it was set through $GITHUB_ENV, and `R-CMD-check.yaml` records that
# `setup-r-dependencies` writes NOT_CRAN=true there, which is why the workflow
# sets it as STEP-level env (GitHub resolves step env last) and why this check
# exists at all.
if (!identical(Sys.getenv("NOT_CRAN"), "false")) {
  stop("source-guards: NOT_CRAN is \"", Sys.getenv("NOT_CRAN"),
       "\", expected \"false\". The fit tests in these files would attempt ",
       "real fits against a DLL this job does not build. Set it as ",
       "step-level env.")
}

# How a test reaches package source. Patterns 3 and 4 match nothing today and
# are kept as cheap insurance for a guard written with a direct `readLines()`.
patterns <- c(
  "test_path[(].*[.][.]",   # test_path("..", "..", "R", ...)
  "[.]docs_root[(][)]",
  "readLines[(].*R/[0-9]",
  "readLines[(].*src/TMB"
)

# The expected SET, by name. A count is not a set: an earlier version asserted
# `length(targets) == 16L`, and a review demonstrated that deleting a real guard
# while adding one decorative test that happens to match a pattern keeps the
# count at 16 and the job green. Only 5 of these 16 are named anywhere else
# (the contributor recipe), so 11 could be deleted silently. Pinning the names
# makes the diff say which guard left.
EXPECTED <- c(
  "test-composition-age-hat-width.R",
  "test-docs-anchors.R",
  "test-dynamics-recruitment-minage.R",
  "test-dynamics-sex-index-bounds.R",
  "test-likelihood-caal-afsc.R",
  "test-linkage-encode.R",
  "test-linkage-selectivity-apical.R",
  "test-linkage-srr-r-init-level.R",
  "test-mse-cap-and-hcr2-threshold.R",
  "test-schema-canonical.R",
  "test-schema-cpp-dispatch.R",
  "test-schema-jnll-rows.R",
  "test-schema-quantity-dictionary.R",
  "test-schema-registries.R",
  "test-selectivity-random-sel-block.R",
  "test-vignette-api.R"
)

# Blocks whose source read is INCIDENTAL to a fit test, so they legitimately
# skip without a compiled DLL. Named `file::test_that description` so the
# exemption is per BLOCK, not per file -- a file-wide exemption is the hole
# this job exists to close. Each needs a reason.
# (Two entries run past 80 characters: they are exact block labels, and
# trimming one would stop it matching.)
INCIDENTAL_BLOCKS <- c(
  # Reads src/TMB only to find the JNLL_CAAL row index, for assertions about a
  # fitted model. All three blocks need a fit, and this job loads the package
  # with compile = FALSE, so they skip here.
  "test-likelihood-caal-afsc.R::CAAL_distribution = 'MultinomialAFSC' fits instead of erroring",
  "test-likelihood-caal-afsc.R::the CAAL AFSC value is the AFSC form, computed from its own inputs",
  "test-likelihood-caal-afsc.R::the AFSC CAAL family is a different likelihood from the multinomial"
)

all_tests <- list.files(test_dir, pattern = "^test-.*[.][rR]$",
                        full.names = TRUE)
reads_source <- vapply(all_tests, function(f) {
  ln <- readLines(f, warn = FALSE)
  any(vapply(patterns, function(p) any(grepl(p, ln)), logical(1)))
}, logical(1))
targets <- all_tests[reads_source]

if (!setequal(basename(targets), EXPECTED)) {
  stop(sprintf(
    paste0("source-guards: the set of source-reading tests changed.\n",
           "  no longer discovered: %s\n  newly discovered:     %s\n",
           "Update EXPECTED, and say in the commit which guard was added or ",
           "removed and why."),
    paste(setdiff(EXPECTED, basename(targets)), collapse = ", "),
    paste(setdiff(basename(targets), EXPECTED), collapse = ", ")))
}

# Which `test_that()` blocks read package source. Attributing lines to the
# block they sit in lets the measured-nothing check run per block. A read in a
# top-level helper, outside any block, is attributed to the file as a whole.
source_reading_blocks <- function(f) {
  ln <- readLines(f, warn = FALSE)
  # Both spellings. Matching only the bare form missed every block in the four
  # guard files written entirely as `testthat::test_that(` -- 29 of them -- so
  # the measured-nothing check below covered none of those files. A skip added
  # to one of them left this job reporting "all passed".
  starts <- grep("^\\s*(testthat::)?test_that\\(", ln)
  if (!length(starts)) return(character())
  # The block label, as testthat reports it in the result's `test` column.
  labels <- sub("^\\s*(testthat::)?test_that\\(\\s*[\"'](.*?)[\"']\\s*,.*$",
                "\\2", ln[starts])
  ends <- c(starts[-1] - 1L, length(ln))
  hits <- character()
  for (i in seq_along(starts)) {
    body <- ln[starts[i]:ends[i]]
    if (any(vapply(patterns, function(p) any(grepl(p, body)), logical(1)))) {
      hits <- c(hits, labels[i])
    }
  }
  hits
}


# Blocks in a file that reads source only in a TOP-LEVEL helper. The helper's
# output is what they all assert on, but this job runs at NOT_CRAN=false
# against no DLL, so a `skip_on_cran()` fit block in one of them skips
# legitimately. They are therefore held to a weaker rule than a body-matching
# block: no result row, or ran-but-asserted-nothing, is a failure; a skip is
# not. The whole-file rule below is what catches a blackout.
helper_read_blocks <- function(f) {
  if (length(source_reading_blocks(f))) return(character())
  block_labels(f)
}


# Every label in the allow-list must name a real block in a real target, or
# the exemption is silently dead: a renamed block drops out of the per-block
# check and the job still reports it as exempt.
block_labels <- function(f) {
  ln <- readLines(f, warn = FALSE)
  starts <- grep("^\\s*(testthat::)?test_that\\(", ln)
  sub("^\\s*(testthat::)?test_that\\(\\s*[\"'](.*?)[\"']\\s*,.*$", "\\2",
      ln[starts])
}

# A target with no parsed `test_that` block means the block regex has gone
# stale, which silently disables the per-block check for that whole file.
# That is how the namespaced spelling hid 29 blocks, so assert it rather than
# trust it.
unparsed <- basename(targets)[vapply(targets, function(f) {
  !length(grep("^\\s*(testthat::)?test_that\\(", readLines(f, warn = FALSE)))
}, logical(1))]
if (length(unparsed)) {
  stop("source-guards: no test_that block parsed in ",
       paste(unparsed, collapse = ", "),
       ". The block regex is stale, so the measured-nothing check would skip ",
       "these files entirely.")
}

known <- unlist(lapply(targets, function(f) {
  paste0(basename(f), "::", block_labels(f))
}), use.names = FALSE)
if (length(setdiff(INCIDENTAL_BLOCKS, known))) {
  stop("source-guards: these INCIDENTAL_BLOCKS entries name no block that ",
       "exists:\n  ", paste(setdiff(INCIDENTAL_BLOCKS, known),
                            collapse = "\n  "),
       "\nA dead exemption drops its block out of the per-block check while ",
       "the job still counts it as exempt. Fix the label or drop the entry.")
}

cat(sprintf("source-guards: %d source-reading test files, %d exempt blocks\n",
            length(targets), length(INCIDENTAL_BLOCKS)))

pkgload::load_all(root, quiet = TRUE, compile = FALSE)
# A missing .so makes library.dynam2() return quietly and the registration
# failure surfaces as a warning, not an error -- so the namespace can be
# half-built with load_all() reporting success. Assert an R-level registry is
# actually reachable, or every guard below would error for the wrong reason.
if (!is.function(Rceattle:::.rce_column_schema)) {
  stop("source-guards: the package namespace loaded but its registries are ",
       "not reachable, so the guards would fail for the wrong reason.")
}

failed <- character()
no_rows <- character()
dark_blocks <- character()
no_block_cover <- character()

for (f in targets) {
  nm <- basename(f)
  res <- as.data.frame(testthat::test_file(f, package = "Rceattle",
                                           reporter = "silent"))

  if (nrow(res) == 0L) {
    cat(sprintf("  %-46s NO RESULT ROWS\n", nm))
    no_rows <- c(no_rows, nm)
    next
  }

  n_fail <- sum(res$failed) + sum(res$error)
  blocks <- source_reading_blocks(f)
  helper <- helper_read_blocks(f)
  cat(sprintf("  %-46s fail=%-3d pass=%-5d skip=%-3d blocks=%d%s\n",
              nm, n_fail, sum(res$passed), sum(res$skipped), length(blocks),
              if (length(helper)) sprintf(" (+%d via helper)", length(helper))
              else ""))
  if (n_fail > 0) failed <- c(failed, nm)

  # Per-block: a block whose own body reads source must have asserted
  # something and must not have skipped.
  for (lab in blocks) {
    key <- paste0(nm, "::", lab)
    if (key %in% INCIDENTAL_BLOCKS) next
    row <- res[res$test == lab, , drop = FALSE]
    if (nrow(row) == 0L) {
      dark_blocks <- c(dark_blocks, paste0(key, "  (no result row)"))
    } else if (sum(row$skipped) > 0) {
      dark_blocks <- c(dark_blocks, paste0(key, "  (skipped)"))
    } else if (sum(row$passed) == 0L &&
               sum(row$failed) + sum(row$error) == 0L) {
      dark_blocks <- c(dark_blocks, paste0(key, "  (asserted nothing)"))
    }
  }
  # Weaker rule for a helper-read file: a skip is allowed, asserting nothing
  # while not skipping is not.
  for (lab in helper) {
    key <- paste0(nm, "::", lab)
    if (key %in% INCIDENTAL_BLOCKS) next
    row <- res[res$test == lab, , drop = FALSE]
    if (nrow(row) && sum(row$skipped) == 0L && sum(row$passed) == 0L &&
        sum(row$failed) + sum(row$error) == 0L) {
      dark_blocks <- c(dark_blocks, paste0(key, "  (ran, asserted nothing)"))
    }
  }

  # Whole-file rule, and the one that closes the demonstrated hole: a target
  # whose EVERY block skipped has exercised its source read not at all, and
  # the per-block loop above cannot say so for a helper-read file. Adding a
  # `skip()` to the one block that still ran used to leave the job green.
  # A file whose every block is allow-listed is exempt by declaration -- that
  # is caal-afsc, whose three blocks all need a fit this job does not build.
  all_exempt <- all(paste0(nm, "::", block_labels(f)) %in% INCIDENTAL_BLOCKS)
  if (sum(res$passed) == 0L && n_fail == 0L && !all_exempt) {
    no_block_cover <- c(no_block_cover,
                        paste0(nm, "  (every block skipped)"))
  }
}

report <- function(what, items) {
  if (!length(items)) return(invisible())
  cat("\nsource-guards: ", what, ":\n  ", paste(items, collapse = "\n  "),
      "\n", sep = "")
}
report("guards that returned no result rows (file-level skip?)", no_rows)
report("source-reading BLOCKS that measured nothing", dark_blocks)
report("targets whose every block skipped", no_block_cover)

if (length(failed)) {
  stop("source-guards: failures in ", paste(failed, collapse = ", "))
}
# A target that asserted nothing at all is a failure, not a note: its source
# read was never exercised, whichever block the read sits in. That is the case
# a `skip()` added to the last running block of a helper-read file produced,
# and the job used to report it as "all passed".
if (length(no_rows) || length(dark_blocks) || length(no_block_cover)) {
  stop("source-guards: ",
       length(no_rows) + length(dark_blocks) + length(no_block_cover),
       " guard(s) measured nothing, which is the failure this job exists to ",
       "catch")
}
cat("source-guards: all passed\n")
