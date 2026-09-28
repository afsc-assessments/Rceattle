# A species =, sex = or fleet = filter on a linkage spec keeps only those levels
# of a term `by` stratifies on. On a spec whose `by` lacks the term there is
# nothing to filter, and until 5.36.0 that was a silent no-op.

testthat::test_that("a filter on a term absent from `by` warns instead of doing nothing", {
  d <- make_test_data()
  if (!is.null(d$data_list)) d <- d$data_list
  spec <- linkage_spec(~ 1, by = ~ fleet, species = 1L, param = "q")
  testthat::expect_warning(
    Rceattle:::materialize_linkage(spec, "q", d$env_data, strata = list(fleet = 1:2)),
    "`species =` on the catchability linkage for `q` has no effect: `by` does not include `species`")
  spec2 <- linkage_spec(~ 1, by = ~ fleet, fleet = 1L, param = "q")
  testthat::expect_no_warning(
    Rceattle:::materialize_linkage(spec2, "q", d$env_data, strata = list(fleet = 1:2)))
})
