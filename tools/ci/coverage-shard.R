# -----------------------------------------------------------------------------
# One shard of the coverage run, for .github/workflows/test-coverage.yaml.
#
# Coverage has to run the full NOT_CRAN suite serially: covr's line counters are
# per-process, so anything executed in a testthat parallel worker is never merged
# back (see the workflow header). Serial + covr-instrumented, the suite outgrew
# the job's 120-minute budget on PR #109 -- it was cancelled mid-run at ~116
# minutes of `test_check()` with nothing having failed.
#
# So the files are split across COV_SHARDS runners instead. Each shard installs
# and instruments the package itself, runs its own slice, and uploads its own
# cobertura report; Codecov merges uploads that carry the same commit SHA, so the
# per-PR diff and the badge are unchanged. Wall-clock becomes the slowest shard
# rather than the sum of all of them.
#
# Env:
#   COV_SHARD     1-based index of this shard           (required)
#   COV_SHARDS    how many shards the run is split into (required)
#   COV_PLAN_ONLY print the split and exit, running no tests -- for retuning the
#                 weights below, and the only mode that needs no covr
#   RUNNER_TEMP   where the instrumented package is installed (GitHub sets it)
#
# Run one shard locally with:
#   COV_SHARD=1 COV_SHARDS=4 Rscript tools/ci/coverage-shard.R
# -----------------------------------------------------------------------------

shard  <- as.integer(Sys.getenv("COV_SHARD",  "1"))
nshard <- as.integer(Sys.getenv("COV_SHARDS", "1"))
stopifnot(!is.na(shard), !is.na(nshard), nshard >= 1, shard >= 1, shard <= nshard)

test_dir <- normalizePath(file.path("tests", "testthat"), winslash = "/", mustWork = TRUE)
files    <- sort(list.files(test_dir, pattern = "^test-.*\\.[rR]$"))
stopifnot(length(files) > 0)


# --- Cost: measured seconds --------------------------------------------------
# coverage-costs.tsv holds measured seconds per test file, regenerated from the
# coverage-timing-shard-* artifacts the workflow uploads.
#
# It replaced a static proxy that counted fit_mod() / refit calls in the source.
# That proxy was not merely imprecise, it was mis-ranked: it scored
# test-functions-retrospective.R at 1042 for a file that takes 105 s and
# test-functions-sim-mod-comp.R at 118 for one that takes 412 s. Every shard
# carried identical proxy weight (1653) while their real runtimes were 18, 26,
# 94 and 12 minutes, and the 94-minute one eventually overran the job timeout.
# It could not have been corrected from the job logs either, because covr runs
# the suite in a subprocess whose output appears only when it returns -- the
# shard that overran was the one that reported nothing.
#
# A file the table does not know gets COST_UNMEASURED, deliberately larger than
# the median file so a new test is not assumed free and packed alongside a heavy
# one. Counting calls in the source is NOT used as the fallback: that is the
# estimator this replaced, and a wrong number in the right units is harder to
# spot than an obvious placeholder. The next timing run measures it for real.
#
# A recorded 0 is trusted, because it is real: the 36 files carrying one measure
# 0.003 s each, 0.1 minutes for all of them together. They are the trivial
# validation tests. What made those zeros dangerous was never their value but
# the tie-break below.
COST_UNMEASURED <- 60

COST_TABLE <- local({
  p <- file.path("tools", "ci", "coverage-costs.tsv")
  if (!file.exists(p)) return(NULL)
  tab <- utils::read.delim(p, stringsAsFactors = FALSE)
  stats::setNames(as.numeric(tab$seconds), tab$file)
})

file_cost <- function(path) {
  # Single-bracket, NOT [[: COST_TABLE is a named ATOMIC vector, and `[[` on a
  # name it does not carry throws "subscript out of bounds" rather than
  # returning NULL. That made the COST_UNMEASURED fallback below unreachable and
  # turned any new test file into a hard CI failure across every shard. `[`
  # yields NA for an absent name, and NULL[nm] is length 0 when the tsv is
  # missing entirely, so both fall through to the placeholder as intended.
  measured <- unname(COST_TABLE[basename(path)])
  if (length(measured) == 1L && !is.na(measured)) return(measured)
  COST_UNMEASURED
}

costs <- vapply(file.path(test_dir, files), file_cost, numeric(1))
names(costs) <- files

# Tracked separately rather than by testing `costs == COST_UNMEASURED`, which
# would also count a file that genuinely measured 60 s.
is_unmeasured <- vapply(files, function(f) {
  m <- unname(COST_TABLE[f])
  !(length(m) == 1L && !is.na(m))
}, logical(1))


# --- Longest-processing-time bin packing -------------------------------------
# Heaviest file first, each one onto whichever shard is currently lightest. For a
# spread like this suite's it lands within a few percent of a perfect split, and
# unlike round-robin it is stable: adding a light test file does not reshuffle
# the heavy ones onto different shards.
#
# Ties in load go to the shard holding the fewest files. `which.min()` alone
# always answers with the lowest index, and a 0-cost file leaves the load it was
# placed on unchanged, so the same bin stays lightest and takes every one of
# them: 36 zero-cost files put 59 of 171 on a single runner, which was then
# cancelled at the 120-minute timeout. Any run of equal-cost files does this;
# equal-and-zero is only the case that occurs here. File count is the tie-break
# because it is what actually differs then.
assign_shards <- function(costs, nshard) {
  ord    <- order(costs, decreasing = TRUE)
  load   <- numeric(nshard)
  count  <- integer(nshard)
  bin    <- integer(length(costs))
  for (i in ord) {
    light    <- which(load == min(load))
    b        <- light[which.min(count[light])]
    bin[i]   <- b
    load[b]  <- load[b] + costs[i]
    count[b] <- count[b] + 1L
  }
  list(bin = bin, load = load)
}

plan <- assign_shards(costs, nshard)
mine <- files[plan$bin == shard]
stopifnot(length(mine) > 0)


# --- Report the whole plan, so the log explains itself ------------------------
message(sprintf("Coverage shard %d of %d -- %d of %d test files",
                shard, nshard, length(mine), length(files)))
for (b in seq_len(nshard)) {
  message(sprintf("  shard %d: %3d files, %5.1f min%s",
                  b, sum(plan$bin == b), plan$load[b] / 60,
                  if (b == shard) "   <- this runner" else ""))
}
n_unmeasured <- sum(is_unmeasured)
if (n_unmeasured > 0) {
  message(sprintf(
    "\n%d of %d files have no measured cost and are counted at %d s. Regenerate\ntools/ci/coverage-costs.tsv from the coverage-timing-shard-* artifacts to\nreplace the placeholders with real times.",
    n_unmeasured, length(files), COST_UNMEASURED))
}

message("\nFiles in this shard (estimated seconds, from coverage-costs.tsv):")
for (f in mine) message(sprintf("  %-52s %5.0f", f, costs[[f]]))

if (nzchar(Sys.getenv("COV_PLAN_ONLY"))) {
  message("\nCOV_PLAN_ONLY set -- not running the tests.")
  quit(save = "no", status = 0)
}


# --- Run it under covr -------------------------------------------------------
# testthat's `filter` is matched against the file name with the `test-` prefix
# and the extension stripped, so anchor an alternation of exactly this shard's
# names. Every test file in the package is [a-z0-9_-] only, so none of them carry
# a regex metacharacter (test-schema-canonical.R would be the place to pin that
# if it ever changes).
context_names <- sub("\\.[rR]$", "", sub("^test-", "", mine))
filter <- paste0("^(", paste(context_names, collapse = "|"), ")$")

# Per-file timings, written as the run goes.
#
# covr::package_coverage(code = ) evaluates `code` in a subprocess and surfaces
# its output only when that returns, so on a shard killed by `timeout-minutes`
# the progress reporter's per-file times are never flushed -- the job log jumps
# from "* DONE (Rceattle)" straight to "The operation was canceled". The shard
# that overruns is the one that reports nothing.
#
# CoverageTimingReporter appends each file's time to this path as that file
# ends, so the table is on disk even when the job is cancelled; the workflow
# uploads it with `if: always()`. Regenerate coverage-costs.tsv from those
# artifacts whenever the suite's shape changes.
timing_log <- Sys.getenv("COV_TIMING_LOG", unset = "")
if (!nzchar(timing_log)) {
  timing_log <- file.path(getwd(), sprintf("coverage-timing-shard-%d.tsv", shard))
}
timing_src <- normalizePath(file.path("tools", "ci", "coverage-timing.R"),
                            winslash = "/", mustWork = TRUE)
message(sprintf("\nPer-file timings -> %s", timing_log))

# `load_package = "installed"` picks up covr's instrumented build rather than the
# source tree, and `package =` makes the test environment's parent the Rceattle
# namespace so the internal helpers the suite calls resolve (see CLAUDE.md).
#
# stop_on_failure = FALSE so a failing test does not throw away the shard's
# coverage. Instead the subprocess writes the passing-expectation count, the
# failing and erroring tests, and the files that produced no result at all to
# `failures_log`, and this script stops after the coverage report is written.
# This is the only per-PR job that runs the skip_on_cran() suite, since
# R-CMD-check sets NOT_CRAN=false.
#
# The subprocess asserts NOT_CRAN itself, BEFORE the fits: this whole job's
# value is that the skip_on_cran() suite runs, and NOT_CRAN going false presents
# as a fast green run of nothing but the unguarded unit tests. 195 of 271 files
# carry the skip, 95 of them at file level, where a skipped file reports no row
# at all rather than a skip -- so every shard would still pass expectations and
# report green. Asserting the mode is the only way to see it.
#
# `[[` not `$` on the result columns: if testthat ever stops returning `error`,
# `results$error` is NULL, `results$failed > 0 | NULL` is logical(0), and a
# logical(0) row index selects NO rows whatever `failed` holds -- a shard with
# real failures would write an empty list and pass. The stopifnot() pins the
# column contract that the indexing depends on (DESCRIPTION floors testthat at
# 3.2.0, which does not promise it).
#
# MultiReporter leaves the progress reporter's console output exactly as it was
# and writes the timing rows alongside it.
failures_log <- tempfile(sprintf("coverage-failures-shard-%d-", shard), fileext = ".txt")
code <- sprintf(
  'source(%s)
   if (!identical(Sys.getenv("NOT_CRAN"), "true"))
     stop("NOT_CRAN must be true for this shard: every skip_on_cran() test ",
          "would otherwise leave the run, 195 of 271 files, and the shard ",
          "would still report green")
   results <- as.data.frame(testthat::test_dir(
     %s, package = "Rceattle", load_package = "installed",
     filter = %s, stop_on_failure = FALSE,
     reporter = testthat::MultiReporter$new(reporters = list(
       testthat::ProgressReporter$new(),
       CoverageTimingReporter$new(%s)))))
   stopifnot(all(c("file", "test", "passed", "failed", "error") %%in%%
                 names(results)))
   failed <- results[results[["failed"]] > 0 | results[["error"]],
                     c("file", "test")]
   writeLines(c(as.character(sum(results[["passed"]])),
                sprintf("FAIL %%s: %%s", failed[["file"]],
                        ifelse(is.na(failed[["test"]]), "(outside test_that)",
                               failed[["test"]])),
                sprintf("NORESULT %%s",
                        setdiff(%s, unique(results[["file"]]))),
                "--shard-complete--"), %s)',
  deparse1(timing_src), deparse1(test_dir), deparse1(filter),
  deparse1(timing_log), deparse1(mine), deparse1(failures_log)
)

install_path <- Sys.getenv("RUNNER_TEMP", unset = "")
install_path <- if (nzchar(install_path)) {
  file.path(normalizePath(install_path, winslash = "/"), "package")
} else {
  tempfile("Rceattle-covr-")
}
dir.create(install_path, recursive = TRUE, showWarnings = FALSE)

cov <- covr::package_coverage(
  quiet        = FALSE,
  clean        = FALSE,
  type         = "none",
  code         = code,
  install_path = install_path
)

covr::to_cobertura(cov)
message(sprintf("\nShard %d coverage: %s", shard, format(covr::percent_coverage(cov))))

# covr keeps the test output in Rceattle-commands.Rout, which
# test-coverage.yaml uploads as the coverage-test-failures-shard-N artifact, so
# name the failing tests here.
#
# The sentinel is written last, so a half-written file cannot read as a clean
# shard. Be honest about its reach: a subprocess that DIES is already caught
# upstream, because covr's run_commands() stops on a non-zero exit before this
# line is reached. What is left is a child that exits 0 having written only
# part of the file, which has not been observed. The sentinel is kept because
# it turns that case, and a missing or empty file, into a named error instead
# of a bare readLines() failure.
shard_results <- readLines(failures_log)
if (length(shard_results) == 0L ||
    shard_results[length(shard_results)] != "--shard-complete--") {
  stop("Shard ", shard, " did not finish writing its result: the test ",
       "subprocess died before its last line. Read the ",
       "coverage-test-failures-shard-", shard, " artifact.")
}
shard_passed <- as.integer(shard_results[1])
failed_tests <- grep("^FAIL ", shard_results, value = TRUE)
no_result    <- sub("^NORESULT ", "",
                    grep("^NORESULT ", shard_results, value = TRUE))

# Count what executed, not that it passed: a green shard should say how much it
# measured. The files with no result are the ones this job cannot see at all --
# a file-level skip_on_covr() reports no row rather than a skip.
message(sprintf("\nShard %d ran %d passing expectations in %d of %d files",
                shard, shard_passed, length(mine) - length(no_result),
                length(mine)))
if (length(no_result) > 0) {
  message(sprintf("  no result from %d file%s (skipped at file level): %s",
                  length(no_result), if (length(no_result) == 1) "" else "s",
                  paste(no_result, collapse = ", ")))
}

if (length(failed_tests) > 0) {
  stop("Shard ", shard, " has ", length(failed_tests), " failing or erroring tests ",
       "(full output in the coverage-test-failures-shard-", shard, " artifact):\n",
       paste0("  ", sub("^FAIL ", "", failed_tests), collapse = "\n"))
}
if (shard_passed == 0) {
  stop("Shard ", shard, " passed no expectations, so its tests were skipped ",
       "(NOT_CRAN false, or every file skipped at file level under covr)")
}
