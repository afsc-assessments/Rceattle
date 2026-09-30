# switch_check() sources its per-column fleet_control defaults from the canonical
# schema (R/0-column_schema.R) via .rce_apply_default(). These guards assert the
# schema-driven path reproduces the historical defaulting behavior: the same
# values, and the same messages (including the non-parametric/Hake-gated ones).

# Collect the messages emitted while evaluating `expr`.
capture_messages <- function(expr) {
  msgs <- character(0)
  withCallingHandlers(
    force(expr),
    message = function(m) {
      msgs[[length(msgs) + 1]] <<- conditionMessage(m)
      invokeRestart("muffleMessage")
    })
  msgs
}

test_that("schema pins the exact fleet_control defaults switch_check applies", {
  # Non-self-referential guard: the expected values are the historical literals
  # that lived in switch_check() before Phase 2.1 (git HEAD~N), NOT re-derived
  # from the schema. This pins both the SET of defaulted columns (so a column
  # added to one side but not the other is caught) and each VALUE (so a wrong
  # schema default cannot slip through a same-source comparison).
  expected <- list(
    Sel_norm_bin = "Off", Sel_norm_bin_upper = "Off",
    Sel_curve_pen1 = 0, Sel_curve_pen2 = 0, Sel_curve_pen3 = 0,
    Sel_start_year = NA, Sel_pen_first_bin = NA, Sel_pen_last_bin = NA,
    Sel_shape_mode = NA, Sel_avgsel_pen = 0, Sel_cap_bin = NA,
    Sel_norm_scope = "AcrossSexes",
    Selectivity_dimension = "Age", Comp_distribution = "MultinomialAFSC",
    CAAL_distribution = "Multinomial", Index_distribution = "Lognormal",
    CAAL_weights = 1, Comp_accum_young = NA, Comp_accum_old = NA, Month = 0,
    Ageing_error_index = NA)

  schema <- .rce_column_schema()
  fc_defaulted <- vapply(
    Filter(function(r) r$sheet == "fleet_control" && isTRUE(r$has_default), schema),
    function(r) r$name, character(1))

  # Exact set of wired fleet_control defaults.
  expect_setequal(unname(fc_defaulted), names(expected))
  # Each default value matches the historical literal.
  for (nm in names(expected))
    expect_identical(schema[[nm]]$default, expected[[nm]], info = nm)
})

test_that("schema-driven defaults fill the documented values", {
  d <- BS2017SS
  for (col in c("Sel_norm_bin", "Sel_norm_bin_upper", "Sel_start_year",
                "Sel_pen_first_bin", "Sel_pen_last_bin", "Sel_shape_mode",
                "Sel_avgsel_pen", "Sel_cap_bin", "Selectivity_dimension",
                "Comp_distribution", "CAAL_distribution", "Index_distribution", "CAAL_weights",
                "Comp_accum_young", "Comp_accum_old", "Month"))
    d$fleet_control[[col]] <- NULL

  out <- suppressMessages(switch_check(d))
  schema <- .rce_column_schema()
  nf <- nrow(out$fleet_control)

  for (col in c("Sel_avgsel_pen", "Selectivity_dimension", "Comp_distribution",
                "CAAL_distribution", "CAAL_weights", "Month",
                "Sel_norm_bin", "Sel_norm_bin_upper")) {
    expect_equal(unname(out$fleet_control[[col]]),
                 rep(schema[[col]]$default, nf), info = col)
  }
  # NA defaults come back all-NA (including the composition-accumulation columns,
  # which switch_check() must materialize per the schema's has_default contract).
  for (col in c("Sel_shape_mode", "Sel_cap_bin",
                "Comp_accum_young", "Comp_accum_old"))
    expect_true(all(is.na(out$fleet_control[[col]])), info = col)

  # Index_distribution is defaulted then per-NA-filled to Lognormal.
  expect_true(all(out$fleet_control$Index_distribution == "Lognormal"))
})


# A BLANK Selectivity_dimension cell is a different event from an absent column:
# .rce_apply_default() fills only the latter. The blank fill is what keeps the
# partial-assignment idiom working, but it also decides the fleet's selectivity
# dimension, so on a model that estimates growth -- the case where "Length" was
# plausibly meant -- it has to say which fleets it filled.
test_that("a blank Selectivity_dimension takes the default, announced when growth is estimated", {
  d <- BS2017SS
  d$fleet_control$Selectivity_dimension <- NA_character_
  d$fleet_control$Selectivity_dimension[1] <- "Length"

  # growth_model == 0: filled, silently.
  d0 <- d; d0$growth_model <- rep(0L, d0$nspp)
  msgs0 <- capture_messages(out0 <- switch_check(d0))
  expect_false(any(grepl("Selectivity_dimension", msgs0, fixed = TRUE)))
  expect_equal(out0$fleet_control$Selectivity_dimension[1], "Length")
  expect_true(all(out0$fleet_control$Selectivity_dimension[-1] == "Age"))

  # growth_model > 0: same fill, but named.
  d1 <- d; d1$growth_model <- rep(1L, d1$nspp)
  msgs1 <- capture_messages(out1 <- switch_check(d1))
  hit <- grep("'Selectivity_dimension' is blank for fleet\\(s\\)", msgs1, value = TRUE)
  expect_length(hit, 1)
  expect_true(grepl("assuming 'Age'", hit, fixed = TRUE))
  # Exactly the blank fleets are named -- split the list rather than
  # substring-match, or "Pollock" is satisfied by "BT_Pollock".
  named <- trimws(strsplit(
    sub("^.*fleet\\(s\\) (.*) in 'fleet_control'.*$", "\\1", hit), ",", fixed = TRUE)[[1]])
  expect_setequal(named, d1$fleet_control$Fleet_name[-1])
  expect_equal(out1$fleet_control$Selectivity_dimension[1], "Length")

  # No blanks -> no message, whatever growth is doing.
  d2 <- d1; d2$fleet_control$Selectivity_dimension <- "Age"
  msgs2 <- capture_messages(switch_check(d2))
  expect_false(any(grepl("is blank for fleet", msgs2, fixed = TRUE)))
})

test_that("default messages fire for ungated messaged columns and stay silent otherwise", {
  d <- BS2017SS   # single-species, no estimated growth, no CAAL data
  for (col in c("Selectivity_dimension", "Comp_distribution", "Month",
                "Sel_avgsel_pen", "Sel_shape_mode"))
    d$fleet_control[[col]] <- NULL

  msgs <- capture_messages(switch_check(d))

  # Ungated messaged columns always announce their default.
  expect_true(any(grepl("'Comp_distribution' not specified", msgs, fixed = TRUE)))
  expect_true(any(grepl("'Month' not specified", msgs, fixed = TRUE)))

  # Selectivity_dimension's message is gated on estimated growth (growth_model > 0);
  # BS2017SS estimates none, so it stays silent (the default is still applied).
  expect_false(any(grepl("'Selectivity_dimension' not specified", msgs, fixed = TRUE)))

  # Silent-default columns (default_msg = NULL) never announce.
  expect_false(any(grepl("Sel_avgsel_pen", msgs)))
  expect_false(any(grepl("Sel_shape_mode", msgs)))
})

test_that(".rce_apply_default gates optional-input messages on named conditions", {
  schema <- .rce_column_schema()
  gated <- list(
    Selectivity_dimension = "growth_estimated",
    CAAL_distribution     = "has_caal",
    CAAL_weights          = "has_caal",
    Sel_norm_bin_upper    = "sel_norm_upper"
  )
  for (col in names(gated)) {
    cond <- gated[[col]]
    # condition FALSE -> value still filled, but no message
    expect_silent(v_off <- .rce_apply_default(NULL, col, schema,
                                              conditions = stats::setNames(list(FALSE), cond)))
    expect_equal(v_off, schema[[col]]$default)
    # condition TRUE -> the default message fires
    msg <- capture_messages(.rce_apply_default(NULL, col, schema,
                                               conditions = stats::setNames(list(TRUE), cond)))
    expect_true(any(grepl(col, msg, fixed = TRUE)),
                info = paste("expected default message for", col, "when", cond, "is TRUE"))
  }
})

test_that(".rce_apply_default gates the Sel_curve_pen message on np_hake", {
  # The Sel_curve_pen columns carry default_msg_when = "np_hake": the default
  # message only fires when a non-parametric/Hake/LogisticPM fleet is present,
  # matching switch_check()'s historical `if(.np_hake)` guard. NOTE: in a real
  # fit this default path is unreachable -- the penalty-SD conversion block in
  # switch_check() materializes Sel_curve_pen1/2/3 (as all-NA vectors when
  # absent) before the default runs, so the value is never NULL there. This
  # test covers the helper's contract directly, not a live fit path.
  schema <- .rce_column_schema()

  # np_hake = FALSE -> fill 0 silently.
  expect_silent(
    v0 <- .rce_apply_default(NULL, "Sel_curve_pen1", schema, np_hake = FALSE))
  expect_equal(v0, 0)

  # np_hake = TRUE -> fill 0 and announce.
  msg <- capture_messages(
    v1 <- .rce_apply_default(NULL, "Sel_curve_pen1", schema, np_hake = TRUE))
  expect_equal(v1, 0)
  expect_true(any(grepl("'Sel_curve_pen1' not specified", msg, fixed = TRUE)))

  # A present value is returned unchanged, no message, either way.
  expect_silent(
    expect_equal(.rce_apply_default(5, "Sel_curve_pen1", schema, np_hake = TRUE), 5))
})

test_that("Diet_distribution accepts string and integer aliases (like Comp/CAAL)", {
  skip_if_not_installed("Rceattle")
  d <- Rceattle::BS2017MS
  res <- function(v) {
    d$Diet_distribution <- v
    suppressMessages(Rceattle:::switch_check(Rceattle:::clean_data(d)))$Diet_distribution
  }
  # string aliases resolve to the integer codes the C++ diet likelihood reads
  expect_equal(as.numeric(res(rep("DirichletMultinomial", d$nspp))), rep(1, d$nspp))
  expect_equal(as.numeric(res(rep("Multinomial", d$nspp))),          rep(0, d$nspp))
  # integer codes pass through unchanged
  expect_equal(as.numeric(res(rep(1L, d$nspp))), rep(1, d$nspp))
  # the diet likelihood has no AFSC variant, so -1 / "MultinomialAFSC" is rejected
  d$Diet_distribution <- rep("MultinomialAFSC", d$nspp)
  expect_error(suppressMessages(Rceattle:::switch_check(Rceattle:::clean_data(d))),
               "Diet_distribution")
})

test_that("the Type-II MSVPA msmMode alias is 'MSVPA'", {
  expect_equal(unname(.map_switch("MSVPA", msmMode_map, "msmMode")), 1)
  expect_equal(unname(.map_switch("SingleSpecies", msmMode_map, "msmMode")), 0)
  # renamed from the pre-release "TypeIIMSVPA", which is no longer accepted
  expect_false("TypeIIMSVPA" %in% names(msmMode_map))
  expect_error(.map_switch("TypeIIMSVPA", msmMode_map, "msmMode"), "msmMode")
})

# -----------------------------------------------------------------------------
# rearrange_data() is exported and validate_switches() runs from data_check(),
# so both can see a fleet_control that never went through switch_check() -- a
# hand-built fixture, or a workbook predating a newly added switch column. Those
# columns must be schema-defaulted, not read blind: the .data pronoun fails with
# a cryptic "column not found", and an NA reaching TMB as sel_norm_scope is read
# by the C++ as WithinSex, the opposite of the documented default.
#
# (These callers still require the .required_fc identity columns, Sel_start_year
# among them, which the bundled .rda do not carry -- hence the switch_check()
# base below rather than a raw BS2017SS.)
# -----------------------------------------------------------------------------

# A fleet_control that is complete except for the newer switch columns.
fc_without_new_switches <- function() {
  d <- suppressMessages(Rceattle::switch_check(Rceattle::BS2017SS))
  d$fleet_control$Sel_norm_scope <- NULL
  d$fleet_control$Index_distribution <- NULL
  d
}

test_that("rearrange_data() defaults the newer switch columns when absent", {
  skip_if_not_installed("Rceattle")
  d <- fc_without_new_switches()
  out <- suppressMessages(Rceattle::rearrange_data(d))
  # "AcrossSexes" is code 1.
  expect_identical(out$sel_norm_scope, rep(1L, nrow(d$fleet_control)))
  expect_false(anyNA(out$sel_norm_scope))
})

test_that("a blank Sel_norm_scope cell is filled, not passed through as NA", {
  skip_if_not_installed("Rceattle")
  d <- fc_without_new_switches()
  d$fleet_control$Sel_norm_scope <- "WithinSex"
  d$fleet_control$Sel_norm_scope[2] <- NA

  out <- suppressMessages(Rceattle::rearrange_data(d))
  expect_identical(out$sel_norm_scope[2], 1L)   # filled with the default
  expect_identical(out$sel_norm_scope[1], 0L)   # supplied value preserved
})

test_that("validate_switches() skips the newer switches rather than failing", {
  skip_if_not_installed("Rceattle")
  d <- fc_without_new_switches()
  # No error, and no complaint about the columns that are simply not there.
  errs <- Rceattle:::validate_switches(d)
  expect_false(any(grepl("Sel_norm_scope|Index_distribution", errs)))
})
