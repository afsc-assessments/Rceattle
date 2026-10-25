# unweighted_jnll_comp is jnll_comp without data weights: in ordinary fitting the
# composition, CAAL and stomach rows are scored unweighted, and every other row
# equals jnll_comp.
#
# Provenance: through 5.64.7 the copy in ceattle.cpp section 13.15 overwrote the
# unweighted CAAL row with the weighted one (CAAL_weights = 0.5 reported 2.957
# where 5.913 was due) and stopped short of the equilibrium-catch row, leaving it 0.

caal_fit <- function(family, w) {
  d <- make_test_data(nyrs = 8, nages = 5, seed = 1, growth = "vonBertalanffy")
  d$fleet_control$CAAL_distribution <- family
  d$fleet_control$CAAL_weights <- w
  suppressMessages(suppressWarnings(fit_mod(
    data_list = d, inits = NULL, estimateMode = 3, msmMode = 0,
    random_rec = FALSE, growthFun = build_growth(fun = "vonBertalanffy"),
    fit_control = fit_control(verbose = 0))))
}

row_of <- function(m, label) {
  i <- match(label, rownames(m))
  testthat::expect_false(is.na(i), info = label)
  m[i, ]
}

test_that("the unweighted CAAL row drops the CAAL weight", {
  testthat::skip_on_cran()
  for (family in c("Multinomial", "MultinomialAFSC")) {
    m <- caal_fit(family, 0.5)
    w <- row_of(m$quantities$jnll_comp, "CAAL data")
    u <- row_of(m$quantities$unweighted_jnll_comp, "CAAL data")
    expect_gt(sum(abs(w)), 0)          # the CAAL term is reached
    expect_equal(u, w / 0.5, info = family)
  }
})

test_that("the unweighted Dirichlet-multinomial CAAL row does not move with the weight", {
  testthat::skip_on_cran()
  # Under a DM the weight is log theta, and the unweighted row is scored at theta = 1.
  a <- caal_fit("DirichletMultinomial", 0)
  b <- caal_fit("DirichletMultinomial", 1)
  expect_false(isTRUE(all.equal(row_of(a$quantities$jnll_comp, "CAAL data"),
                                row_of(b$quantities$jnll_comp, "CAAL data"))))
  expect_equal(row_of(a$quantities$unweighted_jnll_comp, "CAAL data"),
               row_of(b$quantities$unweighted_jnll_comp, "CAAL data"))
})

test_that("rows without a data weight are the same in both matrices", {
  testthat::skip_on_cran()
  m <- caal_fit("Multinomial", 0.5)
  j <- m$quantities$jnll_comp
  u <- m$quantities$unweighted_jnll_comp
  own <- c("Composition data", "CAAL data", "Stomach content data")
  keep <- setdiff(rownames(j), own)
  expect_true(all(c("Index data", "Catch data") %in% keep))
  expect_gt(sum(abs(j[c("Index data", "Catch data"), ])), 0)
  expect_identical(u[keep, ], j[keep, ])
})

test_that("the initial equilibrium catch row is copied", {
  testthat::skip_on_cran()
  testthat::skip_if_not(exists("GOA2018SS"))
  f <- suppressMessages(suppressWarnings(fit_mod(
    data_list = Rceattle::GOA2018SS, file = NULL, inits = NULL, estimateMode = 3,
    random_rec = FALSE, msmMode = 0, initMode = 6,
    fit_control = fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
  label <- grep("equilibrium", rownames(f$quantities$jnll_comp), ignore.case = TRUE, value = TRUE)
  expect_length(label, 1)
  w <- row_of(f$quantities$jnll_comp, label)
  expect_gt(sum(abs(w)), 0)
  expect_identical(row_of(f$quantities$unweighted_jnll_comp, label), w)
})
