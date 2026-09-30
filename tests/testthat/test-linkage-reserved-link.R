# Regression tests for the reserved-but-unimplemented "logit" link.
#
# `LINKAGE_LINKS` shares its codes with src/TMB/linkage.hpp, but every
# accumulator there gates on linkfn == 1 (log) or linkfn == 0 (identity).
# A logit row therefore used to encode cleanly, be estimated, and
# contribute exactly zero to every offset tensor -- while
# map_linkage_adjuster() still masked its base parameter, fixing the
# process parameter at its init with nothing replacing it. It must be
# rejected up front until an accumulator implements it.

testthat::test_that("logit is reserved in the schema but not implemented", {
  testthat::expect_true("logit" %in% Rceattle:::LINKAGE_LINKS)
  testthat::expect_false("logit" %in% Rceattle:::LINKAGE_LINKS_IMPLEMENTED)
  testthat::expect_setequal(Rceattle:::LINKAGE_LINKS_IMPLEMENTED,
                            c("identity", "log", "exponential"))
})


testthat::test_that("linkage_spec() rejects link = 'logit'", {
  testthat::expect_error(
    Rceattle::linkage_spec(~temp, param = "M1", link = "logit"),
    "reserved but not yet implemented"
  )
})


testthat::test_that("linkage_spec() still accepts the implemented links", {
  for (lk in Rceattle:::LINKAGE_LINKS_IMPLEMENTED) {
    spec <- Rceattle::linkage_spec(~temp, param = "M1", link = lk)
    testthat::expect_s3_class(spec, "Rceattle_linkage_spec")
    testthat::expect_equal(spec$link, lk)
  }
})


testthat::test_that("validate_linkage_table() backstops hand-built logit rows", {
  # linkage_row() bypasses linkage_spec(), so the schema validator has to
  # catch it too -- otherwise a hand-assembled table reaches TMB silently.
  testthat::expect_error(
    Rceattle:::linkage_row(process = "M", param = "M1", X_col = 1L,
                           link = "logit"),
    "reserved but not yet implemented"
  )
})
