testthat::test_that("encode_linkage_for_tmb empty table -> empty encoding", {
  empty <- Rceattle:::new_linkage_table()
  X <- matrix(numeric(0), nrow = 5, ncol = 0)
  out <- Rceattle:::encode_linkage_for_tmb(empty, X)
  testthat::expect_equal(out$n_linkage, 0L)
  testthat::expect_length(out$linkage_process, 0L)
  testthat::expect_length(out$linkage_X_col, 0L)
  testthat::expect_length(out$linkage_prior_p1, 0L)
  testthat::expect_equal(out$linkage_X, X)
})


testthat::test_that("encode_linkage_for_tmb maps strings to expected codes", {
  env <- data.frame(Year = 2000:2004, temp = stats::rnorm(5))
  s_K <- Rceattle::linkage_spec(
    formula = ~ temp,
    param   = "K",
    by      = ~ species,
    link    = "log",
    priors  = list(temp = normal(0, 0.5))
  )
  pool <- Rceattle:::pool_linkages(
    spec_groups = list(growth = list(K = s_K)),
    env_data    = env,
    strata      = list(species = 1:2)
  )
  out <- Rceattle:::encode_linkage_for_tmb(pool$table, pool$X)

  testthat::expect_equal(out$n_linkage, 4L)   # 2 cols x 2 species
  # All rows are growth -> code 2
  testthat::expect_true(all(out$linkage_process ==
                            Rceattle:::LINKAGE_PROCESS_CODES[["growth"]]))
  # All rows are log_K -> code 0 (per growth namespace)
  testthat::expect_true(all(out$linkage_param ==
                            Rceattle:::LINKAGE_PARAM_CODES$growth[["K"]]))
  # link = "log" -> code 1
  testthat::expect_true(all(out$linkage_link ==
                            Rceattle:::LINKAGE_LINK_CODES[["log"]]))
  # X_col is 0-based on output
  testthat::expect_setequal(out$linkage_X_col, c(0L, 1L))
  # species = 1:2; sex / age_bin are NA -> sentinel 0
  testthat::expect_setequal(out$linkage_species, 1:2)
  testthat::expect_true(all(out$linkage_sex == 0L))
  testthat::expect_true(all(out$linkage_age_bin == 0L))

  # Prior on `temp` propagates to its rows; intercept rows = none/NA
  temp_col <- match("temp", out$linkage_design_names) - 1L  # 0-based
  is_temp <- out$linkage_X_col == temp_col
  testthat::expect_true(all(out$linkage_prior_family[is_temp] ==
                            Rceattle:::LINKAGE_PRIOR_CODES[["normal"]]))
  testthat::expect_true(all(out$linkage_prior_p1[is_temp] == 0))
  testthat::expect_true(all(out$linkage_prior_p2[is_temp] == 0.5))
  testthat::expect_true(all(out$linkage_prior_family[!is_temp] ==
                            Rceattle:::LINKAGE_PRIOR_CODES[["none"]]))
  testthat::expect_true(all(is.na(out$linkage_prior_p1[!is_temp])))
})


testthat::test_that("encode_linkage_for_tmb is deterministic on row order", {
  env <- data.frame(Year = 2000:2004, temp = 1:5)
  s_K <- Rceattle::linkage_spec(
    formula = ~ temp,
    param   = "K",
    by      = ~ species
  )
  pool <- Rceattle:::pool_linkages(
    spec_groups = list(growth = list(K = s_K)),
    env_data    = env,
    strata      = list(species = 1:3)
  )
  out1 <- Rceattle:::encode_linkage_for_tmb(pool$table, pool$X)
  out2 <- Rceattle:::encode_linkage_for_tmb(pool$table, pool$X)
  testthat::expect_equal(out1, out2)
})


testthat::test_that("encode_linkage_for_tmb rejects bad inputs", {
  testthat::expect_error(
    Rceattle:::encode_linkage_for_tmb(list(a = 1), matrix(0)),
    "Rceattle_linkage_table"
  )
  empty <- Rceattle:::new_linkage_table()
  testthat::expect_error(
    Rceattle:::encode_linkage_for_tmb(empty, X = "not a matrix"),
    "numeric matrix"
  )
})


testthat::test_that("LINKAGE_PROCESS_CODES has stable assignments", {
  testthat::expect_equal(
    Rceattle:::LINKAGE_PROCESS_CODES,
    c(recruitment = 0L, M = 1L, growth = 2L, q = 3L, sel = 4L, comp = 5L)
  )
})


testthat::test_that("*_LINKAGE_PARAMS stays in sync with LINKAGE_PARAM_CODES", {
  # The user-facing per-process allowed-key vectors live beside their own
  # constructor, one R/0-build_*.R per process, while the cpp-side code map
  # lives in R/0-linkage_encode.R. They must agree on which params
  # exist for each process or the encoder errors at fit time. This
  # test catches drift in either direction.
  testthat::expect_setequal(
    Rceattle:::GROWTH_LINKAGE_PARAMS,
    names(Rceattle:::LINKAGE_PARAM_CODES$growth)
  )
  testthat::expect_setequal(
    Rceattle:::M_LINKAGE_PARAMS,
    names(Rceattle:::LINKAGE_PARAM_CODES$M)
  )
  testthat::expect_setequal(
    Rceattle:::RECRUITMENT_LINKAGE_PARAMS,
    names(Rceattle:::LINKAGE_PARAM_CODES$recruitment)
  )
  testthat::expect_setequal(
    Rceattle:::SEL_LINKAGE_PARAMS,
    names(Rceattle:::LINKAGE_PARAM_CODES$sel)
  )
})


testthat::test_that("the selectivity param codes match the template's consume sites", {
  # linkage.hpp routes each code to a tensor; ceattle.cpp's prior re-target
  # reads the same codes. A code added on one side only is silent at fit time.
  src <- testthat::test_path("..", "..", "src", "TMB")
  testthat::skip_if(!dir.exists(src), "src/TMB not available")
  lk <- paste(readLines(file.path(src, "linkage.hpp"), warn = FALSE), collapse = "\n")
  codes <- sort(unique(unname(Rceattle:::LINKAGE_PARAM_CODES$sel)))
  testthat::expect_equal(codes, 0:5)
  for (k in codes) {
    testthat::expect_match(lk, paste0("param == ", k, "\\b"), perl = TRUE,
                           info = paste("sel code", k, "has no consume site in linkage.hpp"))
  }
})
