testthat::test_that("linkage_spec(species = ...) filters the level grid", {
  env <- data.frame(
    Year = 2000:2004,
    temp = stats::rnorm(5),
    PDO  = stats::rnorm(5)
  )
  spec <- Rceattle::linkage_spec(
    formula = ~ temp,
    param   = "K",
    by      = ~ species,
    species = c(1L, 3L)
  )
  rows <- Rceattle:::materialize_linkage(
    spec, process = "growth",
    env_data = env,
    strata   = list(species = 1:3)
  )
  testthat::expect_setequal(rows$species, c(1L, 3L))
  testthat::expect_false(2L %in% rows$species)
  # 2 design cols x 2 species -> 4 rows
  testthat::expect_equal(nrow(rows), 4L)
})


testthat::test_that("species filter that excludes everything yields zero rows", {
  env <- data.frame(Year = 2000:2002, temp = 1:3)
  spec <- Rceattle::linkage_spec(
    formula = ~ temp,
    param   = "K",
    by      = ~ species,
    species = 99L
  )
  # Since 5.36.0 a filter that matches nothing says so, since the whole spec drops.
  testthat::expect_warning(
    rows <- Rceattle:::materialize_linkage(
      spec, process = "growth",
      env_data = env,
      strata   = list(species = 1:3)),
    "`species = 99` on the growth linkage for `K` matches none")
  testthat::expect_equal(nrow(rows), 0L)
  testthat::expect_setequal(attr(rows, "design_colnames"),
                            c("(Intercept)", "temp"))
})


testthat::test_that("linkage_spec validates species argument", {
  testthat::expect_error(
    Rceattle::linkage_spec(formula = ~ temp, species = c(0L, 1L)),
    "positive 1-based species ids"
  )
  testthat::expect_error(
    Rceattle::linkage_spec(formula = ~ temp, species = NA_integer_),
    "positive 1-based species ids"
  )
})


testthat::test_that("build_growth accepts list of specs for per-species formulas", {
  g <- Rceattle::build_growth(
    fun = "vonBertalanffy",
    linkages = list(
      K = list(
        Rceattle::linkage_spec(formula = ~ temp,
                               by      = ~ species,
                               species = 1L),
        Rceattle::linkage_spec(formula = ~ temp + PDO,
                               by      = ~ species,
                               species = 2L)
      )
    )
  )
  testthat::expect_length(g$linkages$K, 2L)
  testthat::expect_equal(g$linkages$K[[1]]$param, "K")
  testthat::expect_equal(g$linkages$K[[2]]$param, "K")
  testthat::expect_equal(g$linkages$K[[1]]$species, 1L)
  testthat::expect_equal(g$linkages$K[[2]]$species, 2L)
})


testthat::test_that("pool_linkages stitches per-species formulas into one table", {
  env <- data.frame(
    Year = 2000:2004,
    temp = stats::rnorm(5),
    PDO  = stats::rnorm(5)
  )
  specs <- list(
    K = list(
      Rceattle::linkage_spec(formula = ~ temp,
                             by      = ~ species,
                             species = 1L),
      Rceattle::linkage_spec(formula = ~ temp + PDO,
                             by      = ~ species,
                             species = 2L)
    )
  )
  out <- Rceattle:::pool_linkages(
    spec_groups = list(growth = specs),
    env_data    = env,
    strata      = list(species = 1:2)
  )
  testthat::expect_setequal(out$design_names,
                            c("(Intercept)", "temp", "PDO"))

  # Species 1 contributes 2 rows (Intercept, temp).
  # Species 2 contributes 3 rows (Intercept, temp, PDO).
  testthat::expect_equal(nrow(out$table), 5L)
  testthat::expect_equal(sum(out$table$species == 1L), 2L)
  testthat::expect_equal(sum(out$table$species == 2L), 3L)

  # Verify the global X_col mapping: species 1 should never reference PDO.
  pdo_global <- match("PDO", out$design_names)
  testthat::expect_false(any(out$table$species == 1L &
                             out$table$X_col == pdo_global))
})


testthat::test_that("validate_growth_linkages rejects malformed list entries", {
  testthat::expect_error(
    Rceattle::build_growth(
      fun = "vonBertalanffy",
      linkages = list(K = list(Rceattle::linkage_spec(formula = ~ 1),
                                   "not a spec"))
    ),
    "must be a linkage_spec"
  )
})
