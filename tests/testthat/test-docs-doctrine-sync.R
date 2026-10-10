# CLAUDE.md, AGENTS.md and .github/copilot-instructions.md each state the three
# doctrines in full, because each is injected by its own tool rather than read
# through a pointer. That makes them a hand-synced registry, the same shape as
# the JnllRow partners (test-schema-jnll-rows.R): nothing linked them, and
# AGENTS.md's own contract -- "a change to one of the three is a change in all
# three files" -- was already false in eight clauses when it was written. These
# tests make the contract enforced rather than declared.
#
# Pinned by CLAUSE, not by count: a count is satisfied by deleting one clause
# and adding another. And matched on FLATTENED text, because every clause below
# wraps in at least one file, and a line-based grep of a wrapping sentence finds
# nothing and reports a confident 0.
#
# None of the three ships (all are Rbuildignored), so this skips on a tarball
# and runs from the source tree.

.DOCTRINE_FILES <- c("CLAUDE.md", "AGENTS.md", ".github/copilot-instructions.md")

# The operative clauses, not the prose around them. Each is the sentence a
# reader has to act on, and each is currently verbatim in all three files.
# Doctrine 1's headline is deliberately absent: CLAUDE.md writes "Don't change
# an API" where the other two write "Do not", so its operative test is pinned
# instead of its title.
.DOCTRINE_CORE <- c(
  "if a user's existing script would behave differently",   # doctrine 1's test
  "SIMPLIFY-LOG.md",                                        # where an idea goes
  "Write the simplest thing that works, and nothing for later",
  "A helper earns its place at two callers",
  "A guard that cannot fire is not safety",
  "is not optional",                     # the registries that already exist
  ".index_rows_natural_scale",           # the one with a wrong-residual path
  "Write for a fisheries scientist",
  "never strip one"                      # a citation is the specification
)

.doctrine_root <- function() {
  cands <- c(".", testthat::test_path("..", ".."))
  ok <- file.exists(file.path(cands, "DESCRIPTION")) &
    file.exists(file.path(cands, "CLAUDE.md"))
  cands <- cands[ok]
  if (!length(cands)) testthat::skip("source tree not available")
  normalizePath(cands[1])
}

# Flattened to one line, runs of whitespace collapsed, so a clause that wraps
# still matches.
.doctrine_flat <- function(path) {
  gsub("[[:space:]]+", " ", paste(readLines(path, warn = FALSE), collapse = " "))
}


testthat::test_that("the three agent-instruction files state one doctrine core", {
  root <- .doctrine_root()
  paths <- file.path(root, .DOCTRINE_FILES)
  testthat::expect_true(all(file.exists(paths)))
  flat <- vapply(paths, .doctrine_flat, character(1))

  # CLAUDE.md is the canonical source, and checking it first is the positive
  # control: if a clause were renamed in all three, or mistyped here, these
  # assertions fail rather than the file-to-file comparison passing on an
  # absence found everywhere.
  for (cl in .DOCTRINE_CORE) {
    testthat::expect_true(grepl(cl, flat[[1]], fixed = TRUE), info = cl)
  }
  for (i in seq_along(.DOCTRINE_FILES)[-1]) {
    for (cl in .DOCTRINE_CORE) {
      testthat::expect_true(grepl(cl, flat[[i]], fixed = TRUE),
                            info = paste(.DOCTRINE_FILES[i], "|", cl))
    }
  }
})


testthat::test_that("CONTRIBUTING.md points at the doctrines rather than restating them", {
  # Item D's acceptance criterion, which governs the CLAUDE.md / CONTRIBUTING.md
  # pair: "no rule stated in full in both files". A human reads CONTRIBUTING.md
  # on GitHub and can follow a link, so it carries the minimum it must and
  # points for the rest -- unlike the three machine files above, where a pointer
  # is a hope. The registry carve-out is exempt because missing it silently
  # yields a log-scale residual, so a human must meet it where they are.
  root <- .doctrine_root()
  flat <- .doctrine_flat(file.path(root, "CONTRIBUTING.md"))
  for (cl in c("Write the simplest thing that works, and nothing for later",
               "A helper earns its place at two callers",
               "A guard that cannot fire is not safety")) {
    testthat::expect_false(grepl(cl, flat, fixed = TRUE), info = cl)
  }
  # What it must still reach a human with, and the pointer that gets them there.
  testthat::expect_true(grepl(".index_rows_natural_scale", flat, fixed = TRUE))
  testthat::expect_true(grepl("CLAUDE.md", flat, fixed = TRUE))
})
