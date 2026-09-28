# Every switch code exists three times: as a name in an R map (R/0-switches.R or
# a R/0-build_*.R constructor), as a string literal in build_map() or
# build_params(), and as a `case` label in the C++. Nothing binds them. A code
# added to a map with no branch to run, or a branch with no name to select it,
# is invisible until a user picks the option and gets silence or the wrong arm.
#
# This test binds the first and third: it reads the `case` labels out of the
# template and requires them to be exactly the map's values, minus a PINNED
# exemption in each direction. The exemptions are the point -- they are the
# machine-checked inventory of "implemented in C++ but unreachable from R" and
# "reachable from R but handled before dispatch", which otherwise lives as prose
# in three different phrasings.

# Strip comments before scanning, for two reasons. A retired feature is commented
# out rather than deleted here -- `// DATA_INTEGER(est_diet);` and
# `// DATA_INTEGER(avgnMode);` both sit at the top of ceattle.cpp -- so a scanner
# reading the raw text calls a dead declaration live. And the whole Kinzey
# predation block in predation.hpp is inside a /* ... */, so its `switch
# (msmMode)` would otherwise be reported as live dispatch; it is also what
# unbalances the brace count for every file concatenated after it.
strip_cpp_comments <- function(src) {
  src <- gsub("(?s)/\\*.*?\\*/", "", src, perl = TRUE)
  gsub("//[^\n]*", "", src, perl = TRUE)
}

# Read `case N:` labels belonging to `switch (var)`, at brace depth 1 only. A
# nested switch -- switch(flt_sex) inside switch(flt_type) -- must not leak its
# labels into the outer one.

cpp_switch_cases <- function(src, var) {
  starts <- gregexpr(paste0("switch\\s*\\(\\s*", var, "\\s*(\\([^)]*\\))?\\s*\\)"),
                     src, perl = TRUE)[[1]]
  if (starts[1] == -1L) return(integer(0))

  # Brace depth over the whole file, computed once and vectorised. A per-
  # character R loop over the concatenated template costs ~30 s; this is
  # milliseconds, and the drift check has to be cheap enough to always run.
  ob <- gregexpr("{", src, fixed = TRUE)[[1]]
  cb <- gregexpr("}", src, fixed = TRUE)[[1]]
  ob <- ob[ob > 0L]; cb <- cb[cb > 0L]
  br <- c(ob, cb)
  step <- c(rep(1L, length(ob)), rep(-1L, length(cb)))
  o <- order(br); br <- br[o]; step <- step[o]
  depth <- cumsum(step)                     # depth immediately AFTER each brace

  out <- integer(0)
  for (st in starts) {
    open_i <- which(br > st & step == 1L)
    if (!length(open_i)) next
    k <- open_i[1]; body_depth <- depth[k]   # depth just inside the switch body
    shut <- which(seq_along(br) > k & depth == body_depth - 1L)
    if (!length(shut)) next
    from <- br[k]; to <- br[shut[1]]

    block <- substr(src, from, to)
    m <- gregexpr("case\\s+-?\\d+\\s*:", block, perl = TRUE)[[1]]
    if (m[1] == -1L) next
    abs_pos <- from + m - 1L
    # A case belongs to THIS switch only at depth 1 of its body; a nested
    # switch -- switch(flt_sex) inside switch(flt_type) -- sits deeper.
    at <- findInterval(abs_pos, br)
    keep <- depth[pmax(at, 1L)] == body_depth
    txt <- regmatches(block, gregexpr("case\\s+-?\\d+\\s*:", block, perl = TRUE))[[1]]
    # Keep the sign: comp_ll_type dispatches `case -1:` (MultinomialAFSC).
    out <- c(out, as.integer(sub("case\\s+(-?\\d+)\\s*:", "\\1", txt[keep])))
  }
  sort(unique(out))
}

test_that("every C++ dispatch branch matches the R map that selects it", {
  dir <- c("src/TMB", testthat::test_path("..", "..", "src", "TMB"))
  dir <- dir[dir.exists(dir)]
  testthat::skip_if(length(dir) == 0, "src/TMB not available")
  files <- list.files(dir[1], pattern = "\\.(cpp|hpp)$", full.names = TRUE)
  src <- paste(vapply(files, function(f) paste(readLines(f, warn = FALSE), collapse = "\n"),
                      character(1)), collapse = "\n")
  src <- strip_cpp_comments(src)

  srr_all <- c(.SRR_FUNS, .SRR_DEPRECATED_FUNS)

  # var = the C++ switch; values = what R can encode; exemptions carry a reason,
  # and the reason is why the exemption is allowed to exist at all.
  spec <- list(
    # The selectivity switch dispatches on `sel_case`, not `sel_type`:
    # NonParametricIntegrable (13) picks its construction from Time_varying_sel,
    # so the case is computed just above the switch. Every user-selectable code
    # still reaches it unchanged; the one internal label is pinned below.
    list(var = "sel_case", values = sel_map,
         map_only = c(Fixed = 0),
         map_only_why = "Fixed selectivity is applied before the dispatch, not by it",
         cpp_only = integer(0), cpp_only_why = character(0)),

    list(var = "flt_type", values = fleet_map,
         map_only = c(Off = 0),
         map_only_why = "an Off fleet is dropped from the likelihood before dispatch",
         cpp_only = integer(0), cpp_only_why = character(0)),

    list(var = "comp_ll_type", values = comp_loglike_map,
         map_only = integer(0), map_only_why = character(0),
         cpp_only = integer(0), cpp_only_why = character(0)),

    list(var = "caal_ll_type", values = comp_loglike_map,
         map_only = integer(0), map_only_why = character(0),
         cpp_only = integer(0), cpp_only_why = character(0)),

    list(var = "est_sigma_index", values = estimate_sd_map,
         map_only = integer(0), map_only_why = character(0),
         cpp_only = integer(0), cpp_only_why = character(0)),

    # No exemption: `case 2` (the Ludwig and Walters analytical sigma) was added
    # in 5.12.0, so this now matches est_sigma_index exactly. The exemption that
    # stood here recorded the gap; this test failing when it was fixed is the
    # guard working, not a regression.
    list(var = "est_sigma_fsh", values = estimate_sd_map,
         map_only = integer(0), map_only_why = character(0),
         cpp_only = integer(0), cpp_only_why = character(0)),

    list(var = "diet_ll_type", values = diet_loglike_map,
         map_only = integer(0), map_only_why = character(0),
         cpp_only = integer(0), cpp_only_why = character(0)),

    list(var = "HCR", values = hcr_map,
         map_only = integer(0), map_only_why = character(0),
         cpp_only = integer(0), cpp_only_why = character(0)),

    list(var = "srr_fun", values = srr_all,
         map_only = integer(0), map_only_why = character(0),
         cpp_only = integer(0), cpp_only_why = character(0)),

    list(var = "srr_pred_fun", values = srr_all,
         map_only = integer(0), map_only_why = character(0),
         cpp_only = integer(0), cpp_only_why = character(0)),

    list(var = "growth_model", values = .GROWTH_FUN_TO_INT,
         map_only = c(empirical = 0),
         map_only_why = "empirical growth reads the weight-at-age input; it has no growth curve to dispatch",
         cpp_only = 3L,
         cpp_only_why = "non-parametric growth is stubbed and calls error('not yet implemented')"),

    # msmMode is deliberately absent from this list. The only `switch (msmMode)`
    # in the template sits inside the commented-out Kinzey predation block
    # (predation.hpp), so there is NO live dispatch to compare against -- the
    # implemented modes are handled by `if (msmMode == ...)` in ceattle.cpp.
    # A row here would have to be justified by dead code; see the separate test
    # below, which pins the absence instead.
    NULL
  )

  spec <- Filter(Negate(is.null), spec)
  for (s in spec) {
    cpp <- cpp_switch_cases(src, s$var)
    testthat::expect_gt(length(cpp), 0)
    r <- sort(unique(unname(s$values)))

    # A C++ branch R cannot select. Anything beyond the pinned set is a code
    # someone implemented and no user can reach.
    testthat::expect_setequal(setdiff(cpp, r), sort(unique(as.integer(s$cpp_only))))

    # A value R can produce with no branch to run it -- either handled before
    # dispatch, or a silent fall-through to `default`.
    testthat::expect_setequal(setdiff(r, cpp), sort(unique(as.integer(s$map_only))))
  }
})

test_that("the not-yet-implemented inventory is stated, not just tolerated", {
  # Each stubbed branch above must actually be stubbed, so the exemption cannot
  # quietly become a real-but-unreachable implementation.
  #
  # skip_on_cran: the msmMode half builds a fixture, which is the whole runtime
  # of this file. The dispatch comparison above stays static and sub-second.
  testthat::skip_on_cran()
  dir <- c("src/TMB", testthat::test_path("..", "..", "src", "TMB"))
  dir <- dir[dir.exists(dir)]
  testthat::skip_if(length(dir) == 0, "src/TMB not available")

  growth <- paste(readLines(file.path(dir[1], "growth.hpp"), warn = FALSE), collapse = "\n")
  testthat::expect_match(growth, "Non-parametric growth not yet implemented", fixed = TRUE)

  # The Kinzey-Punt predation forms are not "stubbed" in the sense growth is --
  # the whole block is commented out, so msmMode 3-9 have no branch at all, live
  # or erroring. Pin that, so the day someone uncomments it this test says so.
  pred <- paste(readLines(file.path(dir[1], "predation.hpp"), warn = FALSE), collapse = "\n")
  testthat::expect_match(pred, "switch (msmMode)", fixed = TRUE)
  testthat::expect_length(cpp_switch_cases(strip_cpp_comments(pred), "msmMode"), 0)

  # msmMode 3-9 must be refused in R, not silently accepted and dispatched.
  d <- make_test_data(nyrs = 4, nages = 3)
  d$msmMode <- 5
  testthat::expect_error(suppressWarnings(suppressMessages(data_check(d))),
                         "Kinzey")
})


test_that("every schema tmb_target is really produced and really consumed", {
  # The third leg of the triangle. test-schema-canonical.R pins the schema
  # against the docs, the test above pins the C++ dispatch against the maps, and
  # this pins the rename in between: a column whose `tmb_target` names an object
  # rearrange_data() never assigns, or that the template never declares, is a
  # column that silently reaches the model as nothing.
  rd  <- c("R/5-rearrange_data.R", testthat::test_path("..", "..", "R", "5-rearrange_data.R"))
  cpp <- c("src/TMB/ceattle.cpp",  testthat::test_path("..", "..", "src", "TMB", "ceattle.cpp"))
  rd  <- rd[file.exists(rd)]; cpp <- cpp[file.exists(cpp)]
  testthat::skip_if(length(rd) == 0 || length(cpp) == 0, "source not available")

  rearrange <- paste(readLines(rd[1],  warn = FALSE), collapse = "\n")
  # Commented out is retired, not declared; see strip_cpp_comments() above.
  template  <- strip_cpp_comments(paste(readLines(cpp[1], warn = FALSE), collapse = "\n"))

  schema  <- .rce_column_schema()
  targets <- vapply(schema,
                    function(r) if (is.null(r$tmb_target)) NA_character_ else r$tmb_target,
                    character(1))
  targets <- targets[!is.na(targets)]
  testthat::expect_gt(length(targets), 20)

  for (col in names(targets)) {
    tgt <- targets[[col]]

    # rearrange_data() must assign it -- unless the column already carries the
    # template's name, in which case it passes through untouched and there is
    # no rename to find.
    if (!identical(tgt, col)) {
      testthat::expect_match(
        rearrange,
        paste0("data_list(\\$", tgt, "|\\[\\[\"", tgt, "\"\\]\\])\\s*<-"), perl = TRUE,
        info = paste0(col, " -> ", tgt, ": rearrange_data() never assigns it"))
    }

    # ...and the template must declare it, or nothing reads it.
    testthat::expect_match(
      template, paste0("DATA_[A-Z]+\\(\\s*", tgt, "\\s*\\)"), perl = TRUE,
      info = paste0(col, " -> ", tgt, ": ceattle.cpp declares no such DATA_ object"))
  }
})


test_that("every renamed template input is either a schema tmb_target or pinned", {
  # The loop above walks only the rows that already carry a `tmb_target`, so on
  # its own it cannot see a rename nobody recorded. This closes that end: every
  # object rearrange_data() assigns AND the template declares must be accounted
  # for -- as some column's `tmb_target`, or by name in the exemption below.
  #
  # A column reaching the model under a name the schema never mentions is the
  # same failure `tmb_target` exists to catch, one level up.
  rd  <- c("R/5-rearrange_data.R", testthat::test_path("..", "..", "R", "5-rearrange_data.R"))
  cpp <- c("src/TMB/ceattle.cpp",  testthat::test_path("..", "..", "src", "TMB", "ceattle.cpp"))
  rd  <- rd[file.exists(rd)]; cpp <- cpp[file.exists(cpp)]
  testthat::skip_if(length(rd) == 0 || length(cpp) == 0, "source not available")

  rearrange <- paste(readLines(rd[1],  warn = FALSE), collapse = "\n")
  # Commented out is retired, not declared; see strip_cpp_comments() above.
  template  <- strip_cpp_comments(paste(readLines(cpp[1], warn = FALSE), collapse = "\n"))

  declared <- regmatches(template,
                         gregexpr("DATA_[A-Z_]+\\(\\s*[A-Za-z0-9_.]+\\s*\\)", template, perl = TRUE))[[1]]
  declared <- unique(trimws(sub("^DATA_[A-Z_]+\\(", "", sub("\\)$", "", declared))))

  # Both spellings of the assignment, so writing one as data_list[["x"]] <- does
  # not quietly drop it out of the guard's view.
  assigned <- regmatches(rearrange, gregexpr(
    "data_list(\\$[A-Za-z0-9_.]+|\\[\\[\"[A-Za-z0-9_.]+\"\\]\\])\\s*<-",
    rearrange, perl = TRUE))[[1]]
  assigned <- unique(gsub("^data_list(\\$|\\[\\[\")|(\"\\]\\])?\\s*<-$", "", assigned,
                          perl = TRUE))

  built <- intersect(assigned, declared)

  # Template inputs that are NOT a single workbook column, so no schema row can
  # own them. Each is either a whole data sheet assembled by rearrange_data()
  # (the *_obs / *_ctl / *_n triples, the biological sheets) or a value that
  # arrives from somewhere other than the workbook -- Ftarget_percent and
  # Flimit_percent come from build_hcr(), env_index from the environmental
  # covariates, index_cov_const from index_cov.
  not_a_column <- c(
    "age_error", "age_trans_matrix", "caal_ctl", "caal_n", "caal_obs",
    "catch_ctl", "catch_n", "catch_obs", "comp_ctl", "comp_n", "comp_obs",
    "diet_ctl", "diet_obs", "emp_sel_ctl", "emp_sel_obs", "env_index",
    "Flimit_percent", "Ftarget_percent", "index_cov_const", "index_ctl",
    "index_n", "index_obs", "lengths", "maturity", "n_stomach_obs",
    "NByageFixed", "ration_data", "sex_ratio", "stomach_id", "weight_obs")

  schema  <- .rce_column_schema()
  targets <- vapply(schema,
                    function(r) if (is.null(r$tmb_target)) NA_character_ else r$tmb_target,
                    character(1))
  targets <- unname(targets[!is.na(targets)])

  testthat::expect_setequal(setdiff(built, targets), not_a_column)

  # And the exemption may not quietly outlive the thing it exempts.
  testthat::expect_setequal(intersect(not_a_column, built), not_a_column)
})


# The selectivity switch carries one label the dispatch scan above cannot see,
# because that scan reads integer `case` literals and this one is named. It is
# the internal construction NonParametricIntegrable (13) uses when its
# deviations are a random walk. Pin it here: it must exist, it must be negative
# so it can never collide with a Selectivity code, and no user-facing map may
# offer it.
testthat::test_that("the internal selectivity dispatch label is out of the user's range", {
  dir <- c("src/TMB", testthat::test_path("..", "..", "src", "TMB"))
  dir <- dir[dir.exists(dir)]
  testthat::skip_if(length(dir) == 0, "src/TMB not available")
  sel <- paste(readLines(file.path(dir[1], "selectivity.hpp"), warn = FALSE),
               collapse = "\n")

  testthat::expect_match(sel, "switch (sel_case)", fixed = TRUE)
  testthat::expect_match(sel, "case SEL_CASE_NP_INTEGRABLE_WALK:", fixed = TRUE)

  val <- sub(".*SEL_CASE_NP_INTEGRABLE_WALK\\s*=\\s*(-?\\d+).*", "\\1",
             regmatches(sel, regexpr("SEL_CASE_NP_INTEGRABLE_WALK\\s*=\\s*-?\\d+", sel)))
  testthat::expect_lt(as.integer(val), 0L)
  testthat::expect_false(as.integer(val) %in% as.integer(sel_map))

  # It is reached only from Selectivity 13 under the random-walk mode.
  testthat::expect_match(
    sel, "sel_type == 13 && flt_varying_sel(flt) == 4", fixed = TRUE)
})
