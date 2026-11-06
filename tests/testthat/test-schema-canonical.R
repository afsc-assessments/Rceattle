# Canonical workbook-column schema (R/0-column_schema.R): faithful-reproduction
# and integrity guards. These prove the single-source schema reproduces the
# hand-maintained copies BEFORE any consumer is switched over to it.

# Compare descriptions on their meaningful text: collapse CR/LF and runs of
# spaces so a drift guard flags wording changes, not cosmetic whitespace.
norm_desc <- function(x) {
  x <- gsub("[\r\n]+", "\n", x)
  x <- gsub("[ \t]+", " ", x)
  x <- gsub(" ?\n ?", "\n", x)
  trimws(x)
}

test_that("schema is structurally sound", {
  schema <- .rce_column_schema()
  nms <- vapply(schema, function(r) r$name, character(1))

  # Unique canonical names; names() keyed correctly.
  expect_equal(anyDuplicated(nms), 0L)
  expect_equal(names(schema), unname(nms))

  # Valid sheets.
  sheets <- vapply(schema, function(r) r$sheet, character(1))
  expect_true(all(sheets %in% c("control", "fleet_control", "bioenergetics_control")))

  # No alias collides with a canonical name (would make a rename ambiguous).
  aliases <- unlist(lapply(schema, function(r) r$aliases))
  expect_length(intersect(aliases, nms), 0)

  # No alias is shared by two columns (would make an upgrade ambiguous).
  expect_equal(anyDuplicated(aliases), 0L)

  # The orphan set is pinned: these are documented historically but read
  # nowhere and slated for removal. Promoting one to "live" (or adding a new
  # orphan) must be a deliberate, visible change -- pinning catches it, so the
  # file's "explicit, testable record of the drop" claim actually holds.
  orphan <- vapply(Filter(function(r) r$status == "orphan", schema),
                   function(r) r$name, character(1))
  expect_setequal(unname(orphan), c("R_sexr", "Catch_units", "Sex"))

  # Orphans never leak into the live accessors used by the generators.
  expect_length(intersect(.rce_schema_names("fleet_control"), orphan), 0)
  expect_length(intersect(.rce_schema_names("control"), orphan), 0)
})

test_that("committed meta_data_names.xlsx matches the schema generator exactly", {
  # The bundled xlsx is generated from the schema (data-raw/regenerate_meta_xlsx.R);
  # this drift guard fails if the two diverge (i.e. the schema was edited without
  # re-running the generator, or the asset was hand-edited).
  gen <- .rce_build_meta_data_df()
  xlsx <- as.data.frame(readxl::read_xlsx(
    system.file("extdata", "meta_data_names.xlsx", package = "Rceattle")))

  expect_equal(dim(xlsx), dim(gen))
  expect_equal(colnames(xlsx), colnames(gen))

  blank <- function(x) { x <- as.character(x); x[is.na(x)] <- ""; x }
  expect_equal(blank(xlsx[["Sheet name"]]), blank(gen[["Sheet name"]]))
  expect_equal(blank(xlsx[["Column/row name"]]), blank(gen[["Column/row name"]]))
  # Descriptions compared on meaningful text (whitespace-normalized), row-aligned.
  # USE.NAMES = FALSE so the comparison is on the normalized text, not names
  # carrying the raw strings (which would defeat the whitespace tolerance).
  expect_equal(vapply(blank(xlsx$Description), norm_desc, character(1), USE.NAMES = FALSE),
               vapply(blank(gen$Description), norm_desc, character(1), USE.NAMES = FALSE))
})

test_that("schema pins the known legacy -> canonical column aliases", {
  # These are the deprecated names accepted on read and upgraded in place
  # (read_data / switch_check). Pinning them here guards the alias table the
  # generator and the rename cascades rely on.
  expected <- list(
    Estimate_index_sd   = "Estimate_survey_sd",
    Index_sd            = c("Survey_sd_prior", "Index_sd_prior"),
    Catch_sd            = "Catch_sd_prior",
    Catchability_index    = "Q_index",
    Catchability_init     = c("Q_prior", "Q_init"),
    Catchability_prior_sd = "Q_sd_prior",
    Comp_distribution     = "Comp_loglike",
    CAAL_distribution     = "CAAL_loglike",
    Index_distribution    = "Index_loglike",
    N_sel_bins          = "Nselages",
    Time_varying_sel_sd = c("Sel_sd_prior", "Time_varying_sel_sd_prior"),
    Time_varying_q_sd   = "Time_varying_q_sd_prior",
    Catchability        = "Estimate_q",
    Bin_first_selected  = character(0),
    Sel_norm_bin        = c("Age_max_selected", "Sel_norm_bin1"),
    Sel_norm_bin_upper  = c("Age_max_selected_upper", "Sel_norm_bin2"),
    Observation_units   = c("Weight1_Numbers2", "weight1_Numbers2"),
    Proj_F_proportion   = "proj_F_prop",
    Sel_pen_first_bin   = "Sel_pen_first_age",
    Sel_pen_last_bin    = "Sel_pen_last_age",
    Sel_cap_bin         = "Sel_cap_age"
  )
  schema <- .rce_column_schema()
  for (canon in names(expected))
    expect_setequal(schema[[canon]]$aliases, expected[[canon]])
})

test_that("R/data.R @format documents exactly the fleet_control schema columns", {
  # The roxygen field dictionary in R/data.R must stay in lock-step with the
  # schema: every live fleet_control column that appears in the meta sheet is
  # documented, and no stale/orphan item lingers. Guards against a rename or a
  # newly surfaced column that forgets its roxygen entry.
  candidates <- c("R/data.R",
                  testthat::test_path("..", "..", "R", "data.R"))
  data_r <- candidates[file.exists(candidates)]
  testthat::skip_if(length(data_r) == 0, "R/data.R source not available")
  lines <- readLines(data_r[1])

  # Bound the fleet_control @format block: from its `\describe` header to the
  # `"BS2017SS"` object line that closes the roxygen block.
  start <- grep("fleet_control: controls", lines)
  end   <- grep('^"BS2017SS"', lines)
  end   <- min(end[end > start])
  block <- lines[start:end]

  item_lines <- grep("\\\\item\\{", block, value = TRUE)
  labs <- sub(".*\\\\item\\{([^}]*)\\}.*", "\\1", item_lines)
  labs <- trimws(unlist(strsplit(labs, ",\\s*")))   # split combined {a, b, c}

  sc <- .rce_schema_rows("fleet_control")
  sc_meta <- vapply(Filter(function(r) isTRUE(r$meta), sc),
                    function(r) r$name, character(1))

  testthat::expect_setequal(labs, sc_meta)
})


test_that("R/data.R documents the switch codes the maps actually define", {
  # The @format block is the field dictionary a user reads before filling in a
  # workbook, so a code documented there that the map does not define is an
  # invented switch -- and a map value the block omits is an option nobody can
  # find. The prose itself is deliberately NOT compared against the schema
  # `doc`: the two are written in different registers (Rd markup here, plain
  # workbook text there), and forcing them equal would strip the \code{}
  # markup that makes the help page readable. The codes are the part that must
  # agree, and they are machine-checkable.
  candidates <- c("R/data.R", testthat::test_path("..", "..", "R", "data.R"))
  data_r <- candidates[file.exists(candidates)]
  testthat::skip_if(length(data_r) == 0, "R/data.R source not available")
  lines <- readLines(data_r[1])

  start <- grep("fleet_control: controls", lines)
  end   <- min(grep('^"BS2017SS"', lines))
  block <- lines[start:end]
  item  <- grep("\\\\item\\{", block, value = TRUE)
  labs  <- sub(".*\\\\item\\{([^}]*)\\}.*", "\\1", item)
  text  <- sub("^#' \\\\item\\{[^}]*\\}\\{(.*)\\}\\s*$", "\\1", item)

  # column -> the map that defines its valid values
  documented_maps <- list(
    Fleet_type       = fleet_map,
    Selectivity      = sel_map,
    Time_varying_sel = tv_sel_map,
    Time_varying_q   = tv_q_map
  )

  for (col in names(documented_maps)) {
    i <- match(col, labs)
    testthat::expect_false(is.na(i), info = paste(col, "has no @format item"))
    map <- documented_maps[[col]]

    # Pairs are written both ways: `1 = "Logistic"` and `0 or 'Off' = ...`.
    pairs <- regmatches(text[i],
      gregexpr('[0-9]+ *(?:=|or) *["\']([A-Za-z0-9_-]+)["\']', text[i], perl = TRUE))[[1]]
    testthat::expect_gt(length(pairs), 0)

    code <- as.integer(sub('^([0-9]+).*', '\\1', pairs))
    name <- sub('.*["\']([A-Za-z0-9_-]+)["\'].*', '\\1', pairs)

    # Every documented pair must be exactly what the map says.
    for (k in seq_along(code)) {
      testthat::expect_true(name[k] %in% names(map),
        info = paste0(col, ": '", name[k], "' is documented but not in the map"))
      if (name[k] %in% names(map)) {
        testthat::expect_equal(unname(map[[name[k]]]), code[k],
          info = paste0(col, ": '", name[k], "' documented as ", code[k],
                        " but the map says ", unname(map[[name[k]]])))
      }
    }

    # And every live map value must be findable in the docs. `Time_varying_q`
    # is exempt: it is overloaded, holding an env_data column index rather than
    # a mode when Catchability is Environmental or AR1.
    if (col != "Time_varying_q") {
      missing_from_docs <- setdiff(names(map), name)
      testthat::expect_length(missing_from_docs, 0)
    }
  }
})


test_that("R/data.R does not present a deprecated alias as a canonical column", {
  # Hard rule: consume a column by its canonical name. The field dictionary is
  # the one place a user learns which name to type, so a legacy spelling here
  # teaches the wrong one -- and read_data() would silently upgrade it anyway.
  candidates <- c("R/data.R", testthat::test_path("..", "..", "R", "data.R"))
  data_r <- candidates[file.exists(candidates)]
  testthat::skip_if(length(data_r) == 0, "R/data.R source not available")
  lines <- readLines(data_r[1])
  block <- lines[grep("fleet_control: controls", lines):min(grep('^"BS2017SS"', lines))]

  labs <- sub(".*\\\\item\\{([^}]*)\\}.*", "\\1",
              grep("\\\\item\\{", block, value = TRUE))
  labs <- trimws(unlist(strsplit(labs, ",\\s*")))

  schema  <- .rce_column_schema()
  aliases <- unlist(lapply(schema, function(r) r$aliases), use.names = FALSE)
  aliases <- aliases[!is.na(aliases)]

  testthat::expect_length(intersect(labs, aliases), 0)
})
