# The jnll_comp rows are a three-way hand-synced registry:
#
#   src/TMB/ceattle.cpp   the JnllRow enum -- which row each likelihood scores into
#   R/6-rename_output.R   the row LABELS a user reads
#   R/9-profile.R         .JNLL_ROW_AXIS -- what each row's COLUMNS count
#
# Nothing links them at build time, so a row added or reordered in one and not
# the others is silent. The cost is not cosmetic: profile_components() names each
# cell from .JNLL_ROW_AXIS, so a row classified on the wrong axis reports a
# survey's likelihood against a species. These tests read the template and assert
# the three agree.

cpp_source <- function(strip_comments = FALSE) {
  cpp <- c("src/TMB/ceattle.cpp",
           testthat::test_path("..", "..", "src", "TMB", "ceattle.cpp"))
  cpp <- cpp[file.exists(cpp)]
  testthat::skip_if(length(cpp) == 0, "src/TMB/ceattle.cpp not available")
  src <- readLines(cpp[1], warn = FALSE)
  if (strip_comments) src <- .strip_cpp_comments(src)
  src
}

# Blank out `//` tails and `/* ... */` blocks, keeping the line count so any
# reported position still matches the file. Without this a scan for
# `jnll_comp(JNLL_*, col)` matches the Kinzey & Punt ration likelihood, which
# lives inside a block comment -- so the guard would report rows 16 and 17 as
# written and never take its own "a row nothing writes yet" branch. The guard
# then measures less than it claims, which is the failure mode it exists to
# prevent elsewhere.
.strip_cpp_comments <- function(src) {
  src <- sub("//.*$", "", src)
  out <- src
  in_block <- FALSE
  for (i in seq_along(src)) {
    line <- src[i]
    res <- ""
    pos <- 1L
    n <- nchar(line)
    while (pos <= n) {
      if (!in_block) {
        j <- regexpr("/\\*", substring(line, pos))
        if (j < 0) { res <- paste0(res, substring(line, pos)); break }
        res <- paste0(res, substring(line, pos, pos + j - 2))
        pos <- pos + j + 1L
        in_block <- TRUE
      } else {
        j <- regexpr("\\*/", substring(line, pos))
        if (j < 0) break
        pos <- pos + j + 1L
        in_block <- FALSE
      }
    }
    out[i] <- res
  }
  out
}


testthat::test_that("the row registry and the JnllRow enum are the same length", {
  src <- cpp_source()
  n <- regmatches(src, regexpr("JNLL_N_ROWS\\s*=\\s*[0-9]+", src))
  n <- as.integer(sub("\\D+", "", n[nzchar(n)][1]))
  testthat::expect_equal(length(.JNLL_ROW_AXIS), n)
})


testthat::test_that("the row registry uses the labels a user actually sees", {
  # .JNLL_ROW_AXIS is keyed by the display name, so a label reworded in
  # rename_output() and not here would silently stop matching and every cell
  # would fall through to the "not in the row registry" warning.
  f <- c("R/6-rename_output.R",
         testthat::test_path("..", "..", "R", "6-rename_output.R"))
  f <- f[file.exists(f)]
  testthat::skip_if(length(f) == 0, "R/6-rename_output.R not available")
  src <- paste(readLines(f[1], warn = FALSE), collapse = "\n")

  block <- regmatches(src, regexpr(
    "rownames\\(quantities\\$jnll_comp\\)[^\\n]*<-\\s*c\\((?s).*?\\)",
    src, perl = TRUE))
  testthat::expect_length(block, 1)
  labels <- regmatches(block, gregexpr('"[^"]+"', block))[[1]]
  labels <- gsub('"', "", labels, fixed = TRUE)

  testthat::expect_equal(labels, names(.JNLL_ROW_AXIS))
})


testthat::test_that("each row's declared axis matches the column the template writes", {
  src <- cpp_source(strip_comments = TRUE)

  # position -> enum constant, in declaration order, which is the row order of
  # both the matrix and .JNLL_ROW_AXIS.
  decl <- grep("^\\s*JNLL_[A-Z_]+\\s*=\\s*[0-9]+", src, value = TRUE)
  const <- sub("^\\s*(JNLL_[A-Z_]+)\\s*=.*$", "\\1", decl)
  const <- setdiff(const, "JNLL_N_ROWS")
  testthat::expect_equal(length(const), length(.JNLL_ROW_AXIS))

  # Every write, as constant -> column expression.
  writes <- regmatches(src, gregexpr(
    "jnll_comp\\(\\s*(JNLL_[A-Z_]+)\\s*,\\s*([A-Za-z_0-9]+)\\s*\\)", src))
  writes <- unlist(writes)
  testthat::expect_gt(length(writes), 0)
  who <- sub(".*\\((JNLL_[A-Z_]+)\\s*,.*", "\\1", writes)
  col <- sub(".*,\\s*([A-Za-z_0-9]+)\\s*\\)$", "\\1", writes)

  # The loop variables each axis is indexed by. `0` is the column-1 fallback a
  # model-wide term uses, and is legal on any row -- reference-point penalties
  # score both per species and once for the model.
  legal <- list(
    fleet   = c("flt", "index", "0"),
    species = c("sp", "rsp", "slot_col", "0"),
    model   = "0")

  for (i in seq_along(const)) {
    seen <- unique(col[who == const[i]])
    if (!length(seen)) next          # a row nothing writes yet
    axis <- unname(.JNLL_ROW_AXIS[i])
    testthat::expect_true(
      all(seen %in% legal[[axis]]),
      info = paste0(const[i], " (", names(.JNLL_ROW_AXIS)[i], ") is registered as '",
                    axis, "' but the template indexes it by ",
                    paste(setdiff(seen, legal[[axis]]), collapse = ", "),
                    ". Either the registry in R/9-profile.R or the template is wrong."))
  }
})


testthat::test_that("the QAR1 catchability block is gone from the template", {
  # Catchability = "AR1" (QAR1) was refused in 5.12.0 and its template block,
  # unreachable since then, deleted in 5.37.0 along with index_q_rho. The live
  # form is a q linkage, ar1(1 | Year) with `observe`.
  src <- cpp_source()
  testthat::expect_length(grep("est_index_q(flt) == 6", src, fixed = TRUE), 0)
  testthat::expect_length(grep("index_q_rho", src, fixed = TRUE), 0)
})
