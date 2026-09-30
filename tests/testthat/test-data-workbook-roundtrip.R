# write_data() sources the control-sheet and bioenergetics-sheet object names
# and order from the canonical schema (R/0-column_schema.R). These guards assert
# the write -> read round-trip preserves every control and bioenergetics value on
# the bundled datasets, so the schema-driven assembly is faithful.

# Control scalars and bioenergetics scalars are single values recycled across
# species (or per-species vectors); read_data returns them at nspp width, so
# compare the original recycled to the read-back length.
expect_roundtrip_field <- function(orig, back, key) {
  a <- orig[[key]]
  if (is.null(a)) return(invisible())
  b <- back[[key]]
  expect_equal(as.numeric(b), rep(as.numeric(a), length.out = length(b)),
               info = key)
}

roundtrip <- function(d) {
  f <- tempfile(fileext = ".xlsx")
  on.exit(unlink(f))
  suppressWarnings(suppressMessages(write_data(d, f)))
  suppressWarnings(suppressMessages(read_data(f)))
}

CONTROL_KEYS <- c("nspp", "styr", "endyr", "projyr", "nsex", "spawn_month",
                  "nages", "minage", "nlengths", "pop_wt_index", "ssb_wt_index",
                  "alpha_wt_len", "beta_wt_len", "pop_age_transition_index",
                  "sigma_rec", "other_food", "estDynamics")
BIO_KEYS <- c("Ceq", "Cindex", "Pvalue", "fday", "CA", "CB", "Qc", "Tco", "Tcm",
              "Tcl", "CK1", "CK4", "Diet_distribution", "Diet_comp_weights")

# GOA2018MS is not bundled (the multispecies GOA reference reuses GOA2018SS), so
# only the three shipped data lists are exercised.
for (nm in Filter(exists, c("BS2017SS", "BS2017MS", "GOA2018SS"))) {
  test_that(paste("control + bioenergetics round-trip on", nm), {
    d <- get(nm)
    back <- roundtrip(d)
    for (k in c(CONTROL_KEYS, BIO_KEYS)) expect_roundtrip_field(d, back, k)
    # The fleet_control sheet is unchanged by this step but must still survive.
    expect_equal(nrow(back$fleet_control), nrow(d$fleet_control))
    expect_true(all(d$fleet_control$Fleet_code %in% back$fleet_control$Fleet_code))
  })
}

test_that("Diet_distribution / Diet_comp_weights round-trip (bundled data leave them empty)", {
  # Both are NULL in the shipped datasets, so set them explicitly to exercise the
  # two bioenergetics rows added to the schema for the write_data label vector.
  d <- BS2017MS
  d$Diet_distribution <- c(0, 1, 0)
  d$Diet_comp_weights <- c(2, 3, 4)
  back <- roundtrip(d)
  expect_equal(as.numeric(back$Diet_distribution), c(0, 1, 0))
  expect_equal(as.numeric(back$Diet_comp_weights), c(2, 3, 4))
})

test_that("read_data parses the control sheet by name, not position", {
  # The schema (hence write_data) orders alpha/beta before the age-transition
  # index; read_data must recover the same values regardless of row order. Write,
  # reverse the non-dimension control rows, and confirm the read-back is
  # unchanged -- the property that makes the schema-driven row order safe.
  d <- BS2017MS
  f1 <- tempfile(fileext = ".xlsx"); on.exit(unlink(f1), add = TRUE)
  suppressWarnings(suppressMessages(write_data(d, f1)))
  sheets <- readxl::excel_sheets(f1)
  book <- lapply(sheets, function(s) as.data.frame(readxl::read_xlsx(f1, sheet = s)))
  names(book) <- sheets
  ctrl <- book$control
  scal <- ctrl[1:4, , drop = FALSE]                    # nspp/styr/endyr/projyr
  rest <- ctrl[nrow(ctrl):5, , drop = FALSE]           # reverse rows 5+
  book$control <- rbind(scal, rest)
  f2 <- tempfile(fileext = ".xlsx"); on.exit(unlink(f2), add = TRUE)
  writexl::write_xlsx(book, f2)
  d_ord <- suppressWarnings(suppressMessages(read_data(f1)))
  d_rev <- suppressWarnings(suppressMessages(read_data(f2)))
  for (k in c(CONTROL_KEYS, BIO_KEYS))
    if (!is.null(d_ord[[k]]))
      expect_equal(as.numeric(d_rev[[k]]), as.numeric(d_ord[[k]]), info = k)
})


# A data_list need not carry every control object the schema knows: read_data()
# returns only the rows the workbook had, and switch_check() is what fills the
# schema defaults. write_data() assembled the control sheet from the FULL schema
# name list and rbind()-ed the values, and rbind() drops a NULL -- so one absent
# object left more labels than rows and aborted the write with "arguments imply
# differing number of rows". The live Pacific hake workbook hits this: it
# predates alpha_wt_len / beta_wt_len.
test_that("write_data() handles a control block missing a schema object", {
  d <- BS2017SS
  # Control objects BS2017SS does not carry. The weight-length pair has never
  # been on it; L50_mat_len / slope_mat_len are the 5.46.0 maturity-at-length
  # options, which an age-based model has no use for. All four must be dropped
  # rather than written at a schema default -- see the comment below.
  missing_two <- c("alpha_wt_len", "beta_wt_len", "L50_mat_len", "slope_mat_len")
  for (k in missing_two) d[[k]] <- NULL

  f <- tempfile(fileext = ".xlsx"); on.exit(unlink(f))
  expect_no_error(suppressWarnings(suppressMessages(write_data(d, f))))

  ctrl <- as.data.frame(readxl::read_xlsx(f, sheet = "control"))
  # The absent objects are dropped, not written at their schema default: the
  # default belongs to the model, and baking one into the workbook would turn a
  # value switch_check() announces at fit time into one the file asserts.
  expect_false(any(missing_two %in% ctrl$Object))
  expect_equal(nrow(ctrl), length(.rce_schema_names("control")) - length(missing_two))

  # Everything the data_list did carry still round-trips, and nothing is lost.
  back <- suppressWarnings(suppressMessages(read_data(f)))
  for (k in setdiff(CONTROL_KEYS, missing_two))
    expect_roundtrip_field(d, back, k)
  for (k in BIO_KEYS) expect_roundtrip_field(d, back, k)
  expect_false(any(missing_two %in% names(back)))
})


# Same defect, other sheet: the bioenergetics rows were assigned into a
# fixed-height matrix by hard-coded index, so an absent object aborted the write
# with "number of items to replace is not a multiple of replacement length".
# Rows are now keyed by schema name, which also removes the second copy of the
# row order that had to be kept in sync by hand.
test_that("write_data() handles a bioenergetics block missing a schema object", {
  d <- BS2017SS
  d$Tcm <- NULL

  f <- tempfile(fileext = ".xlsx"); on.exit(unlink(f))
  expect_no_error(suppressWarnings(suppressMessages(write_data(d, f))))

  bio <- as.data.frame(readxl::read_xlsx(f, sheet = "bioenergetics_control"))
  expect_false("Tcm" %in% bio$Object)
  # Order still follows the schema, and every other object survives.
  expect_equal(bio$Object, setdiff(.rce_schema_names("bioenergetics_control"), "Tcm"))

  back <- suppressWarnings(suppressMessages(read_data(f)))
  for (k in setdiff(BIO_KEYS, "Tcm")) expect_roundtrip_field(d, back, k)
})
