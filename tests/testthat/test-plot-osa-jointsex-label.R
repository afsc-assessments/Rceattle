# A joint-sex (Sex == 3) composition is ONE decomposition unit: the female and
# male bins are stacked into a single vector, scored under one density, with one
# cell dropped for the sum-to-one constraint rather than one per sex. The OSA
# figure re-bases the male bins onto the age axis so both can be read, which
# makes one unit look like two series -- an assessment author read the previous
# " - male" / " - female" labels as two separate distributions and asked why the
# residuals did not match the likelihood. They do; the label now says so.
# Reported 2026-09-28.

test_that("joint-sex OSA series are labelled as halves of one unit", {
  df <- data.frame(
    species = 1L, sex = c(3L, 3L, 1L, 0L),
    index_label = "age",
    age_length_bin = c(2L, 12L, 2L, 2L),   # bin 12 is male under nages = 10
    source = "comp",
    stringsAsFactors = FALSE)

  out <- Rceattle:::.osa_jointsex(df, nages = c(10L), nlengths = c(10L))

  # The joint pair says so, and says which half.
  expect_equal(out$source[1], "comp - joint, female")
  expect_equal(out$source[2], "comp - joint, male")
  # The male half is re-based onto the same axis as the female half.
  expect_equal(out$age_length_bin[2], 2L)

  # A single-sex row is its own density, so it is NOT labelled joint and is
  # returned untouched -- that distinction is the whole point of the label.
  expect_equal(out$source[3], "comp")
  expect_equal(out$source[4], "comp")
  expect_equal(out$age_length_bin[3:4], c(2L, 2L))
})
