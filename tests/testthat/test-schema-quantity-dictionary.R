# The quantity dictionary is a hand-maintained registry describing what the TMB
# template REPORTs. Nothing links the two at build time, so a quantity added,
# renamed or dropped in src/TMB/ceattle.cpp is silent -- and a cheat sheet that
# quietly stops matching the model is worse than none, because it is trusted.
#
# These tests read the template and a real fit and assert the two agree: every
# documented quantity is actually reported, every reported quantity is
# documented, and the `se` column matches the template's ADREPORT calls.

cpp_source <- function() {
  cpp <- c("src/TMB/ceattle.cpp",
           testthat::test_path("..", "..", "src", "TMB", "ceattle.cpp"))
  cpp <- cpp[file.exists(cpp)]
  testthat::skip_if(length(cpp) == 0, "src/TMB/ceattle.cpp not available")
  readLines(cpp[1], warn = FALSE)
}

# Names in an ACTIVE (non-commented) REPORT()/ADREPORT() call. A commented-out
# REPORT is the reason a bare grep over the template over-counts:
# mature_females and sex_ratio_hat sit behind comments and are not on any fit.
reported_names <- function(src, macro = "REPORT") {
  src <- trimws(src)
  src <- src[!grepl("^(//|/\\*|\\*)", src)]
  pat <- paste0("(^|[^A-Za-z_])", macro, "\\(\\s*([A-Za-z0-9_]+)\\s*\\)")
  m <- regmatches(src, regexpr(pat, src))
  m <- m[nzchar(m)]
  unique(sub(paste0(".*", macro, "\\(\\s*([A-Za-z0-9_]+)\\s*\\)"), "\\1", m))
}


testthat::test_that("the dictionary is structurally well formed", {
  d <- quantity_dictionary()
  testthat::expect_s3_class(d, "data.frame")
  testthat::expect_named(d, c("quantity", "process", "meaning", "units",
                              "dims", "se", "standard_label"))
  testthat::expect_false(any(duplicated(d$quantity)))
  testthat::expect_true(all(nzchar(d$meaning)))
  testthat::expect_true(all(nzchar(d$units)))
  testthat::expect_true(all(nzchar(d$dims)))
  testthat::expect_type(d$se, "logical")
  testthat::expect_false(any(is.na(d$se)))
  # Every meaning is one sentence a scientist can act on, not a paragraph.
  testthat::expect_true(all(nchar(d$meaning) < 400))
})


testthat::test_that("every documented process is one of the declared set", {
  d <- quantity_dictionary()
  known <- c("population", "recruitment", "mortality", "fishing",
             "reference_points", "catchability", "selectivity", "growth",
             "composition", "predation", "likelihood", "linkage", "internal")
  testthat::expect_true(all(d$process %in% known))
})


testthat::test_that("filtering and lookup behave", {
  testthat::expect_equal(nrow(quantity_dictionary("ssb")), 1L)
  testthat::expect_equal(quantity_dictionary("ssb")$standard_label,
                         "spawning_biomass")
  rp <- quantity_dictionary(process = "reference_points")
  testthat::expect_true(all(rp$process == "reference_points"))
  testthat::expect_true("Ftarget" %in% rp$quantity)
  testthat::expect_error(quantity_dictionary(process = "not_a_process"),
                         "Unknown process")
  testthat::expect_warning(quantity_dictionary("not_a_quantity"),
                           "Not in the quantity dictionary")
})


testthat::test_that("`se` matches the template's ADREPORT calls", {
  src <- cpp_source()
  adrep <- reported_names(src, "ADREPORT")
  d <- quantity_dictionary()

  # Everything the dictionary claims has a standard error must be ADREPORT'd.
  claimed <- d$quantity[d$se]
  testthat::expect_true(all(claimed %in% adrep),
                        info = paste("claimed but not ADREPORT'd:",
                                     paste(setdiff(claimed, adrep),
                                           collapse = ", ")))

  # ...and every ADREPORT'd quantity that the dictionary covers must say so.
  # catch_hat is ADREPORTed only inside a retrospective peel's band calculation
  # (adreport_catch), so no fit's sdrep ever holds it.
  covered <- setdiff(intersect(adrep, d$quantity), "catch_hat")
  testthat::expect_true(all(d$se[match(covered, d$quantity)]),
                        info = paste("ADREPORT'd but se = FALSE:",
                                     paste(covered[!d$se[match(covered, d$quantity)]],
                                           collapse = ", ")))
})


testthat::test_that("`se` agrees with the as.data.frame() quantity registry", {
  # .RCEATTLE_QUANTITIES records the same fact for the 23 quantities
  # as.data.frame.Rceattle() tidies, and fills `se`/`lwr`/`upr` from it. Two
  # hand-maintained registries of one fact drift.
  d <- quantity_dictionary()
  reg <- Rceattle:::.RCEATTLE_QUANTITIES
  shared <- intersect(names(reg), d$quantity)
  testthat::expect_gt(length(shared), 0)
  testthat::expect_equal(
    d$se[match(shared, d$quantity)],
    vapply(reg[shared], function(x) isTRUE(x$adreport), logical(1)),
    ignore_attr = TRUE)
})


testthat::test_that("the dictionary covers exactly what a fit reports", {
  # Needs a built model, so it is the one slow test here.
  #
  # Fitted with selectivity_se on so log_sel_at_age and log_sel_at_age_index are
  # reported: they are the only conditional quantities, and with the option off
  # this check could not see them in either direction. It costs nothing here,
  # since getsd = FALSE means no sdreport runs.
  testthat::skip_on_cran()
  set.seed(123)
  sim <- make_msm_test_data(years = 1:30)
  fit <- suppressMessages(fit_mod(
    data_list = sim$data_list, estimateMode = 3, msmMode = 1,
    growthFun = build_growth(fun = "vonBertalanffy"),
    fit_control = fit_control(getsd = FALSE, selectivity_se = TRUE, verbose = 0)))

  reported <- names(fit$quantities)
  documented <- quantity_dictionary()$quantity

  testthat::expect_equal(sort(setdiff(reported, documented)), character(0),
                         info = "reported but undocumented; add them to .QUANT_INFO")
  testthat::expect_equal(sort(setdiff(documented, reported)), character(0),
                         info = "documented but not reported; remove them from .QUANT_INFO")
})


testthat::test_that("documented dimensions match the shapes a fit actually has", {
  # The dims column is the part a user acts on when indexing an array, so a
  # wrong rank there sends them to the wrong margin.
  testthat::skip_on_cran()
  set.seed(123)
  sim <- make_msm_test_data(years = 1:30)
  fit <- suppressMessages(fit_mod(
    data_list = sim$data_list, estimateMode = 3, msmMode = 1,
    growthFun = build_growth(fun = "vonBertalanffy"),
    fit_control = fit_control(getsd = FALSE, selectivity_se = TRUE, verbose = 0)))

  d <- quantity_dictionary()
  # Rank implied by the dims string: count the comma-separated entries. Commas
  # inside a nested call are not dimension separators -- "[21, max(n_flt, nspp)]"
  # is rank 2, not 3 -- so parenthesised groups are removed first.
  implied <- vapply(d$dims, function(s) {
    s <- gsub("^\\[|\\]$", "", s)
    s <- gsub("\\([^()]*\\)", "", s)
    length(strsplit(s, ",")[[1]])
  }, 1L)

  actual <- vapply(d$quantity, function(nm) {
    z <- fit$quantities[[nm]]
    if (is.null(z)) return(NA_integer_)
    dm <- dim(z)
    if (is.null(dm)) 1L else length(dm)
  }, 1L)

  bad <- which(!is.na(actual) & actual != implied)
  testthat::expect_equal(
    length(bad), 0L,
    info = paste("dims rank disagrees for:",
                 paste(sprintf("%s (documented %d, actual %d)",
                               d$quantity[bad], implied[bad], actual[bad]),
                       collapse = "; ")))
})
