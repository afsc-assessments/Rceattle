testthat::test_that("linkage_spec(sex = ...) filters the level grid", {
  env <- data.frame(
    Year = 2000:2004,
    temp = stats::rnorm(5)
  )
  spec <- Rceattle::linkage_spec(
    formula = ~ temp,
    param   = "M1",
    by      = ~ species + sex,
    sex     = 1L
  )
  rows <- Rceattle:::materialize_linkage(
    spec, process = "M",
    env_data = env,
    strata   = list(species = 1:2, sex = 1:2)
  )
  testthat::expect_setequal(rows$sex, 1L)
  testthat::expect_false(2L %in% rows$sex)
  # 2 design cols x 2 species x 1 sex -> 4 rows
  testthat::expect_equal(nrow(rows), 4L)
  testthat::expect_setequal(rows$species, 1:2)
})


testthat::test_that("sex filter that excludes everything yields zero rows", {
  env <- data.frame(Year = 2000:2002, temp = 1:3)
  spec <- Rceattle::linkage_spec(
    formula = ~ temp,
    param   = "M1",
    by      = ~ species + sex,
    sex     = 99L
  )
  # Since 5.36.0 a filter that matches nothing says so, since the whole spec drops.
  testthat::expect_warning(
    rows <- Rceattle:::materialize_linkage(
      spec, process = "M",
      env_data = env,
      strata   = list(species = 1:2, sex = 1:2)),
    "`sex = 99` on the natural mortality linkage for `M1` matches none")
  testthat::expect_equal(nrow(rows), 0L)
  testthat::expect_setequal(attr(rows, "design_colnames"),
                            c("(Intercept)", "temp"))
})


testthat::test_that("sex filter warns, and has no effect, when `by` does not include sex", {
  env <- data.frame(Year = 2000:2002, temp = 1:3)
  spec <- Rceattle::linkage_spec(
    formula = ~ temp,
    param   = "K",
    by      = ~ species,
    sex     = 1L
  )
  testthat::expect_warning(
    rows <- Rceattle:::materialize_linkage(
      spec, process = "growth",
      env_data = env,
      strata   = list(species = 1:2)),
    "`sex =` on the growth linkage for `K` has no effect")
  # 2 design cols x 2 species -> 4 rows; sex column is NA throughout.
  testthat::expect_equal(nrow(rows), 4L)
  testthat::expect_true(all(is.na(rows$sex)))
})


testthat::test_that("sex accepts integer and character forms equivalently", {
  s_int   <- Rceattle::linkage_spec(formula = ~ 1, sex = 1L)
  s_chr_F <- Rceattle::linkage_spec(formula = ~ 1, sex = "Females")
  s_chr_f <- Rceattle::linkage_spec(formula = ~ 1, sex = "female")
  s_chr_m <- Rceattle::linkage_spec(formula = ~ 1, sex = "M")

  testthat::expect_identical(s_int$sex,   1L)
  testthat::expect_identical(s_chr_F$sex, 1L)
  testthat::expect_identical(s_chr_f$sex, 1L)
  testthat::expect_identical(s_chr_m$sex, 2L)

  s_both <- Rceattle::linkage_spec(formula = ~ 1, sex = c("Females", "Males"))
  testthat::expect_identical(s_both$sex, c(1L, 2L))
})


testthat::test_that("character sex form filters the level grid like the integer form", {
  env <- data.frame(Year = 2000:2002, temp = 1:3)
  spec_int <- Rceattle::linkage_spec(
    formula = ~ temp, param = "M1",
    by = ~ species + sex, sex = 2L
  )
  spec_chr <- Rceattle::linkage_spec(
    formula = ~ temp, param = "M1",
    by = ~ species + sex, sex = "Males"
  )
  rows_int <- Rceattle:::materialize_linkage(
    spec_int, process = "M", env_data = env,
    strata = list(species = 1:2, sex = 1:2)
  )
  rows_chr <- Rceattle:::materialize_linkage(
    spec_chr, process = "M", env_data = env,
    strata = list(species = 1:2, sex = 1:2)
  )
  testthat::expect_setequal(rows_int$sex, 2L)
  testthat::expect_setequal(rows_chr$sex, 2L)
  testthat::expect_equal(rows_int$species, rows_chr$species)
})


testthat::test_that("linkage_spec validates sex argument", {
  testthat::expect_error(
    Rceattle::linkage_spec(formula = ~ temp, sex = c(0L, 1L)),
    "positive 1-based sex ids"
  )
  testthat::expect_error(
    Rceattle::linkage_spec(formula = ~ temp, sex = NA_integer_),
    "positive 1-based sex ids"
  )
  testthat::expect_error(
    Rceattle::linkage_spec(formula = ~ temp, sex = "neither"),
    "Females.*Males|1/2"
  )
  testthat::expect_error(
    Rceattle::linkage_spec(formula = ~ temp, sex = c("Females", "bogus")),
    "bogus"
  )
})


testthat::test_that("per-sex multi-spec form pools into one table", {
  env <- data.frame(Year = 2000:2002, temp = 1:3)
  specs <- list(
    Rceattle::linkage_spec(
      formula = ~ 1, by = ~ species + sex, sex = "Females",
      init    = list(`(Intercept)` = log(0.12)),
      priors  = list(`(Intercept)` = Rceattle::prior_normal(log(0.12), 0.1))
    ),
    Rceattle::linkage_spec(
      formula = ~ 1, by = ~ species + sex, sex = "Males",
      init    = list(`(Intercept)` = log(0.20)),
      priors  = list(`(Intercept)` = Rceattle::prior_normal(log(0.20), 0.1))
    )
  )
  pooled <- Rceattle:::pool_linkages(
    spec_groups = list(M = list(M1 = specs)),
    env_data    = env,
    strata      = list(species = 1L, sex = 1:2)
  )
  tbl <- pooled$table
  # 1 design col x 1 species x 1 sex per spec x 2 specs -> 2 rows
  testthat::expect_equal(nrow(tbl), 2L)
  testthat::expect_setequal(tbl$sex, c(1L, 2L))
  fem <- tbl[tbl$sex == 1L, ]
  mal <- tbl[tbl$sex == 2L, ]
  testthat::expect_equal(fem$init, log(0.12))
  testthat::expect_equal(mal$init, log(0.20))
  testthat::expect_setequal(pooled$design_names, "(Intercept)")
})
