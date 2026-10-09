# The contributor recipe (vignettes/articles/adding-a-selectivity-form.Rmd)
# names files, functions, R objects, C++ symbols, switch codes and two lines of
# other tests. A recipe that names a file that has moved sends a newcomer to the
# wrong place with confidence, which is worse than no recipe. So every anchor
# the article quotes is checked here, the way test-schema-cpp-dispatch.R checks
# the maps against the template: read the text, resolve each anchor, fail on the
# first that does not resolve.
#
# The article is not in the built package (vignettes/articles is Rbuildignored),
# so this file skips under R CMD check on a tarball and runs from the source tree.

.docs_root <- function() {
  cands <- c(".", testthat::test_path("..", ".."))
  cands <- cands[file.exists(file.path(cands, "DESCRIPTION")) &
                   file.exists(file.path(cands, "vignettes", "articles",
                                         "adding-a-selectivity-form.Rmd"))]
  if (!length(cands)) testthat::skip("source tree not available")
  normalizePath(cands[1])
}

# Inline-code tokens: `...` on one line. Fenced blocks are excluded by the
# no-newline rule, so pasted test output is not scanned for anchors.
.docs_inline_tokens <- function(txt) {
  m <- regmatches(txt, gregexpr("`[^`\n]+`", txt))[[1]]
  unique(gsub("^`|`$", "", m))
}

.docs_cpp_source <- function(root) {
  files <- list.files(file.path(root, "src", "TMB"), pattern = "\\.(cpp|hpp)$",
                      full.names = TRUE)
  paste(vapply(files, function(f) paste(readLines(f, warn = FALSE), collapse = "\n"),
               character(1)), collapse = "\n")
}

test_that("every file the recipe names exists", {
  root <- .docs_root()
  txt  <- paste(readLines(file.path(root, "vignettes", "articles",
                                    "adding-a-selectivity-form.Rmd"), warn = FALSE),
                collapse = "\n")
  tok  <- .docs_inline_tokens(txt)
  paths <- grep("^(R|src/TMB|tests/testthat|vignettes|inst)/[A-Za-z0-9_./-]+\\.(R|hpp|cpp|Rmd|md)$",
                tok, value = TRUE)
  # A bare test file name resolves under tests/testthat/; a bare root file
  # (CONTRIBUTING.md) at the root.
  bare  <- grep("^(test-[A-Za-z0-9_-]+\\.R|[A-Z]+\\.md)$", tok, value = TRUE)
  paths <- c(paths, ifelse(grepl("^test-", bare), file.path("tests", "testthat", bare), bare))
  testthat::expect_gt(length(paths), 20)
  missing <- paths[!file.exists(file.path(root, paths))]
  testthat::expect_equal(missing, character(0),
                         info = "the article names a file that no longer exists")
})

test_that("every function the recipe names is defined, in R or in the template", {
  root <- .docs_root()
  txt  <- paste(readLines(file.path(root, "vignettes", "articles",
                                    "adding-a-selectivity-form.Rmd"), warn = FALSE),
                collapse = "\n")
  tok  <- .docs_inline_tokens(txt)
  fns  <- sub("\\(\\)$", "", grep("^\\.?[A-Za-z][A-Za-z0-9_.]*\\(\\)$", tok, value = TRUE))
  testthat::expect_gt(length(fns), 12)
  cpp  <- .docs_cpp_source(root)
  ns   <- asNamespace("Rceattle")
  resolves <- vapply(fns, function(f) {
    (exists(f, envir = ns, inherits = TRUE) && is.function(get(f, envir = ns))) ||
      grepl(paste0("(^|[^A-Za-z0-9_])", f, "\\s*\\("), cpp, perl = TRUE)
  }, logical(1))
  testthat::expect_equal(fns[!resolves], character(0),
                         info = "the article names a function that no longer exists")
})

test_that("the R objects and C++ symbols the recipe relies on exist, and are named", {
  root <- .docs_root()
  txt  <- paste(readLines(file.path(root, "vignettes", "articles",
                                    "adding-a-selectivity-form.Rmd"), warn = FALSE),
                collapse = "\n")
  tok  <- .docs_inline_tokens(txt)
  ns   <- asNamespace("Rceattle")
  cpp  <- .docs_cpp_source(root)

  # Pinned rather than scraped: a bare identifier in backticks can be a column
  # name or a string as easily as an object. Each must be in the text too, so a
  # pin cannot outlive the sentence that needed it.
  r_objects   <- c("sel_map", ".PAR_SEL_SLOTS", ".SEL_LINKAGE_WIRED_FORMS",
                   ".SEL_PARAM_TO_SLOT", ".JNLL_ROW_AXIS")
  cpp_symbols <- c("flt_sel_type", "sel_at_age", "sel_at_length", "is_length_based",
                   "JnllRow", "JNLL_SEL_DEV", "REPORT",
                   "Logistic selectivity penalties")
  tok_bare <- sub("\\(\\)$", "", tok)
  testthat::expect_equal(setdiff(c(r_objects, cpp_symbols), tok_bare), character(0),
                         info = "a pinned anchor is no longer named in the article")
  for (o in r_objects) testthat::expect_true(exists(o, envir = ns), info = o)
  for (s in cpp_symbols) testthat::expect_true(grepl(s, cpp, fixed = TRUE), info = s)

  # What the article says those objects contain.
  testthat::expect_true("DoubleNormal" %in% names(get(".PAR_SEL_SLOTS", ns)$sel_inf))
  testthat::expect_true("DoubleNormal" %in% get(".SEL_LINKAGE_WIRED_FORMS", ns))
  testthat::expect_true(all(c("peak", "right_floor", "sigma_asc", "sigma_desc") %in%
                              names(get(".SEL_PARAM_TO_SLOT", ns))))

  # Workbook columns the article names resolve in the schema.
  cols <- c("Selectivity", "Selectivity_index", "Time_varying_sel",
            "Bin_first_selected", "Sel_norm_bin")
  testthat::expect_equal(setdiff(cols, tok), character(0))
  testthat::expect_true(all(cols %in% names(get(".rce_column_schema", ns)())))

  # Source text the article sends the reader to, by file.
  src_of <- function(...) paste(readLines(file.path(root, ...), warn = FALSE), collapse = "\n")
  testthat::expect_match(src_of("R", "1-data_check.R"), "sel_para <-", fixed = TRUE)
  testthat::expect_match(src_of("R", "0-parameter_index.R"), '"8" = "DoubleNormal"', fixed = TRUE)
  # Four hard-coded type lists carry code 8 in the deviate densities.
  main <- src_of("src", "TMB", "ceattle.cpp")
  testthat::expect_match(main, "Logistic selectivity penalties", fixed = TRUE)
  testthat::expect_length(gregexpr("flt_sel_type(flt) == 8", main, fixed = TRUE)[[1]], 4L)
})

test_that("the switch codes and modes the recipe quotes still hold", {
  root <- .docs_root()
  sel  <- paste(readLines(file.path(root, "src", "TMB", "selectivity.hpp"), warn = FALSE),
                collapse = "\n")
  testthat::expect_identical(unname(sel_map[["DoubleNormal"]]), 8)
  testthat::expect_identical(unname(sel_map[["Fixed"]]), 0)
  testthat::expect_false("Fake" %in% names(sel_map))
  # "The next form takes 16": 10 retired, 12 still named by the normalizer, 14
  # freed when the two integrable forms collapsed into 13, and 15 taken in
  # 5.46.0 by DoubleNormalSS3.
  testthat::expect_false(any(c(10, 12, 14, 16) %in% sel_map))
  testthat::expect_identical(unname(sel_map[["DoubleNormalSS3"]]), 15)
  testthat::expect_identical(unname(sel_map[["NonParametricIntegrable"]]), 13)
  testthat::expect_match(sel, "sel_type != 12", fixed = TRUE)
  testthat::expect_match(sel, "case 8:")
  # The switch dispatches on sel_case, not sel_type, so NonParametricIntegrable
  # can pick its construction from Time_varying_sel. Keep the no-default check
  # pointed at the name actually switched on, or it passes vacuously.
  testthat::expect_match(sel, "switch (sel_case)", fixed = TRUE)
  testthat::expect_false(grepl("switch \\(sel_case\\)[^}]*default:", sel, perl = TRUE))
  testthat::expect_true(all(c("NonParametric", "NonParametricPM", "Hake", "LogisticPM",
                              "DoubleLogistic") %in% names(sel_map)))
  testthat::expect_true(all(c("Off", "IID", "AR1", "RandomWalk", "Block",
                              "RandomWalkAscending") %in% names(tv_sel_map)))
  # The DoubleNormal map block frees deviates for exactly these modes; the
  # article's "deviates vanished" entry depends on RandomWalkAscending not being one.
  bm <- paste(readLines(file.path(root, "R", "3-build_map.R"), warn = FALSE), collapse = "\n")
  # From the DoubleNormal branch to the next form's branch. The block still
  # omits RandomWalkAscending deliberately: 5.55.0 refused that combination in
  # data_check() rather than implementing a second mode-5 form, so the
  # article's "deviates vanished" entry still describes this block.
  from  <- regexpr('if \\(sel_type == "DoubleNormal"\\)', bm)
  testthat::expect_gt(from, 0)
  rest  <- substr(bm, from + 30, nchar(bm))
  block <- substr(rest, 1, regexpr('if \\(sel_type == "', rest) - 1)
  testthat::expect_match(block, 'c\\("IID", "AR1", "RandomWalk"\\)')
  testthat::expect_match(block, '"Block"')
  testthat::expect_false(grepl("RandomWalkAscending", block))
})

# A schema `doc` string is written verbatim into meta_data_names.xlsx, so a
# selectivity code named there is a code an assessment author will try to use.
# Retiring form 14 left six of them advertising it, caught only by hand while
# reviewing the 5.34.0-5.41.0 release PR.
test_that("no schema description names a selectivity code sel_map does not accept", {
  schema <- .rce_column_schema()
  sel_cols <- c("N_sel_bins", "Sel_curve_pen1", "Sel_curve_pen2", "Sel_curve_pen3",
                "Sel_shape_sd", "Sel_curvature_sd", "Sel_devmag_sd", "Selectivity")
  docs <- vapply(schema[names(schema) %in% sel_cols], function(r) r$doc, character(1))
  testthat::expect_gt(length(docs), 0)
  # The three shapes a schema description writes a code in. A slash run alone
  # misses the column that matters most: `Selectivity` enumerates every code one
  # per line, so it was the one column shipping code 14 and the only one this
  # check never looked at. Years and sd values are in none of these shapes.
  .named_codes <- function(txt) {
    grab <- function(rx) {
      runs <- unlist(regmatches(txt, gregexpr(rx, txt, perl = TRUE)))
      if (!length(runs)) return(numeric(0))
      as.numeric(unlist(regmatches(runs, gregexpr("\\d{1,2}", runs))))
    }
    # A code parenthesized after a form name, quoted or not -- the schema writes
    # both ("NonParametricPM" (9), and LogisticPM (11) in the same sentence as
    # "type 2/9/13"). Taken in two steps because "2DAR1" and "3DAR1" carry a
    # digit in the name itself, and only the parenthesized number is the code.
    named <- unlist(regmatches(txt, gregexpr(
      "(\"[^\"]+\"|[A-Za-z][A-Za-z0-9_.-]*)\\s*\\(\\d{1,2}\\)", txt, perl = TRUE)))
    paren <- if (length(named)) as.numeric(sub(".*\\((\\d{1,2})\\)$", "\\1", named)) else numeric(0)
    unique(c(grab("\\b\\d{1,2}(\\s*/\\s*\\d{1,2})+"),            # "2/9/13"
             grab("\\b\\d{1,2}(\\s*,\\s*(?:or\\s+)?\\d{1,2})+"), # "2, 5, 6, 7, 9, or 13"
             grab("(?m)^\\s*\\d{1,2}(?=\\s*=)"),                 # "13 = non-parametric ..."
             paren))                                             # "\"NonParametricPM\" (9)"
  }
  for (nm in names(docs)) {
    codes <- .named_codes(docs[[nm]])
    testthat::expect_true(all(codes %in% sel_map),
                          info = paste0(nm, " names selectivity code(s) ",
                                        paste(setdiff(codes, sel_map), collapse = ", "),
                                        ", which sel_map does not accept"))
  }
})

test_that("the two test lines whose failures the recipe pastes still carry those assertions", {
  root <- .docs_root()
  dispatch  <- readLines(file.path(root, "tests", "testthat", "test-schema-cpp-dispatch.R"),
                         warn = FALSE)
  canonical <- readLines(file.path(root, "tests", "testthat", "test-schema-canonical.R"),
                         warn = FALSE)
  # The pasted output names these lines, so the check reads the number OUT of
  # the article rather than hard-coding it: an assertion that moves then fails
  # here until the pasted output is refreshed, and renumbering the article
  # alone cannot satisfy it. Re-run the mutation and paste the new output.
  art <- paste(readLines(file.path(root, "vignettes", "articles",
                                   "adding-a-selectivity-form.Rmd"), warn = FALSE),
               collapse = "\n")
  cited <- function(file) {
    m <- regmatches(art, regexpr(paste0(file, ":\\d+"), art))
    testthat::expect_length(m, 1L)
    as.integer(sub(".*:", "", m))
  }
  testthat::expect_match(dispatch[cited("test-schema-cpp-dispatch.R")],
                         "expect_setequal(setdiff(r, cpp)", fixed = TRUE)
  testthat::expect_match(canonical[cited("test-schema-canonical.R")],
                         "missing_from_docs", fixed = TRUE)
})


# The README's pinning example is the one install command an assessor is told to
# use for management advice, and it names a version. It has gone stale four times
# in the 5.34.0-5.45.0 release alone -- 5.41.0, 5.42.1, 5.43.0 and 5.45.0 -- each
# time because a version landed on dev after the paperwork was written, and each
# time it shipped pointing at a tag that does not exist. A reviewer caught three
# of those; this catches the rest.
test_that("the README pins the version in DESCRIPTION", {
  root <- .docs_root()
  readme <- paste(readLines(file.path(root, "README.md"), warn = FALSE),
                  collapse = "\n")
  ver <- as.character(read.dcf(file.path(root, "DESCRIPTION"), "Version")[1, 1])

  pins <- regmatches(readme,
                     gregexpr("Rceattle@[0-9]+\\.[0-9]+\\.[0-9]+", readme))[[1]]
  # If the example stops naming a version there is nothing to keep in step, but
  # silently passing would retire the check without anyone deciding to.
  expect_gt(length(pins), 0)
  expect_equal(unique(sub("^Rceattle@", "", pins)), ver)
})


# The developer guide's file map claims to cover every file in R/, one row
# each. A map that silently stops covering a new file is worse than no map: a
# reader trusts it and concludes the file does not exist. And a map whose rows
# are merely PRESENT is not enough -- an earlier version of this guard matched
# filenames anywhere in the document, which let three mutations through:
# deleting the `data.R` row (its name is a substring of six others), swapping
# two descriptions, and emptying every description cell. So these tests anchor
# on the table row itself, and check the description says something about the
# file it names.
.dg_file_rows <- function(root) {
  g <- readLines(file.path(root, "vignettes", "articles",
                           "developer-guide.Rmd"), warn = FALSE)
  h <- grep("^## Every file in", g)
  if (!length(h)) return(NULL)
  nxt <- grep("^## ", g)
  nxt <- nxt[nxt > h[1]][1]
  sec <- g[h[1]:(nxt - 1L)]
  m <- regmatches(sec, regexec("^\\| `([^`]+\\.R)` \\| (.+?) \\|\\s*$", sec))
  m <- Filter(function(x) length(x) == 3L, m)
  stats::setNames(vapply(m, `[`, character(1), 3L),
                  vapply(m, `[`, character(1), 2L))
}

test_that("the developer guide has one file-map row per file in R/", {
  root <- .docs_root()
  rows <- .dg_file_rows(root)
  testthat::expect_false(is.null(rows))

  files <- basename(list.files(file.path(root, "R"), pattern = "[.]R$"))
  expect_gt(length(files), 50)   # not vacuous on an empty R/

  expect_setequal(names(rows), files)
  # One row each: a duplicated row would make setequal pass but the map
  # ambiguous.
  expect_equal(length(rows), length(files))

  # A row present but empty, or a placeholder, is a map that lies by omission.
  short <- names(rows)[nchar(trimws(rows)) < 20L]
  expect_true(
    length(short) == 0L,
    info = paste("file-map rows with no real description:",
                 paste(short, collapse = ", ")))
})

# NOT checked here, deliberately: whether a description is CORRECT. Word
# overlap between a row and the file it names was tried and rejected -- the
# large files (6-fit_mod.R, 0-column_schema.R) mention nearly every concept, so
# 17 of 67 accurate descriptions scored higher against some other file than
# against their own, and a guard that fires on correct rows trains people to
# ignore it. Swapping two descriptions between adjacent topics therefore passes
# these tests. That one is on review, and it is why both reviewers of the
# original map read every row against the code.

# The converse: a filename the guide names that no longer exists sends a reader
# to a file that is not there. test-plot-smoke.R carried two function names for
# functions that had never existed, so this class of staleness is real here.
test_that("every R/ filename the developer guide names exists", {
  root <- .docs_root()
  guide <- readLines(file.path(root, "vignettes", "articles",
                               "developer-guide.Rmd"), warn = FALSE)
  files <- basename(list.files(file.path(root, "R"), pattern = "[.]R$"))
  tests <- basename(list.files(file.path(root, "tests", "testthat"),
                               pattern = "[.]R$"))
  raw <- basename(list.files(file.path(root, "data-raw"), pattern = "[.]R$"))
  known <- c(files, tests, raw, "compile.R")

  mentioned <- unique(unlist(regmatches(
    guide, gregexpr("[0-9A-Za-z._-]+[.]R\\b", guide))))
  expect_gt(length(mentioned), 20)
  expect_true(all(mentioned %in% known),
              info = paste("named in the guide but absent:",
                           paste(setdiff(mentioned, known), collapse = ", ")))
})

# Symbol anchors, not just filenames. The guide told readers to find
# `switch (sel_type)` in selectivity.hpp for three releases; the switch is on
# `sel_case`, and the string the guide gave appears nowhere in the file. The
# filename guards above could not see that, because the filename was right.
test_that("every C++ switch the developer guide names exists in that file", {
  root <- .docs_root()
  guide <- paste(readLines(file.path(root, "vignettes", "articles",
                                     "developer-guide.Rmd"), warn = FALSE),
                 collapse = "\n")

  # `switch (var)` quoted in the guide, and the .hpp/.cpp named near it.
  switches <- unique(unlist(regmatches(
    guide, gregexpr("switch \\(([A-Za-z_][A-Za-z_0-9]*)\\)", guide))))
  vars <- unique(sub("^switch \\((.*)\\)$", "\\1", switches))
  # `sel_type` is quoted once as the name that is NOT the switch; allow it.
  vars <- setdiff(vars, "sel_type")
  expect_gt(length(vars), 0)

  src <- paste(unlist(lapply(
    list.files(file.path(root, "src", "TMB"), pattern = "[.](hpp|cpp)$",
               full.names = TRUE),
    readLines, warn = FALSE)), collapse = "\n")

  missing <- vars[!vapply(vars, function(v)
    grepl(paste0("switch\\s*\\(\\s*", v), src), logical(1))]
  expect_true(
    length(missing) == 0L,
    info = paste0("the guide names a C++ switch that does not exist in ",
                  "src/TMB/: ", paste(missing, collapse = ", ")))
})
