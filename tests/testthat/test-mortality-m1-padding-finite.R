# log_M1 is dimensioned to the WIDEST species -- [nspp, max(nsex), max(nages)]
# -- so a species with fewer sexes or fewer ages leaves padding cells. Through
# 5.49.4 the two fill paths disagreed on what went in them:
#
#   build_params()              array(1, ...)  -> padding = log(1) = 0
#   fit_mod(updateM1 = TRUE)    array(0, ...)  -> padding = log(0) = -Inf
#
# Both now write only the real cells and then set each padding cell from its
# own species, so the array is finite and uniform per species. The uniformity
# matters because a parameter sharing a map level with a padding cell starts at
# the mean over that level: a 1 there pulled an estimated M1 toward 1.0/yr.
#
# Measured on GOA2018SS (nspp 3, nsex c(1, 2, 1), nages c(10, 21, 12), so
# max_sex 2 and max_age 21): the updateM1 path produced 53 non-finite cells of
# 126 -- 42 -Inf in the padding sex-2 cells of the two one-sex species, plus 11
# NA. The 11 NA appeared on BOTH paths: species 1's M1_base row is blank past
# its 10th age, which is correct for a 10-age species, and both loops read
# 1:max_age from it.
#
# Three of the template's four log_M1 readers are bounded by the species' own
# dimensions -- the M1_at_age block by sex < nsex(sp) and age < nages(sp), the
# M prior by its nsex_tmp / nage_tmp, and growth.hpp by
# log_M1(sp, sex, nages(sp) - 1). The fourth, the linkage-intercept prior at
# ceattle.cpp:5188, indexes by the linkage row's stratum and stays in range
# because fit_mod() builds per-species strata -- an R-side bound, not a
# template one. So no fit moved, which is why this needs a structural test:
# an objective cannot show it. Verified inert across msmMode 0/1 cold fits,
# M1_model 1/2/3, and all six golden blocks.

testthat::test_that("build_params() log_M1 has no non-finite start value", {
  data("GOA2018SS", package = "Rceattle", envir = environment())

  # The dataset that exposes it: ragged in BOTH dimensions.
  testthat::expect_true(max(GOA2018SS$nsex) > min(GOA2018SS$nsex))
  testthat::expect_true(max(GOA2018SS$nages) > min(GOA2018SS$nages))

  pars <- suppressWarnings(suppressMessages(
    Rceattle::build_params(Rceattle::switch_check(GOA2018SS))
  ))

  testthat::expect_true(all(is.finite(pars$log_M1)))
  testthat::expect_equal(sum(!is.finite(pars$log_M1)), 0L)
})


testthat::test_that("a padding cell takes its own species' M1", {
  # Padding holds the species' own value rather than the array's initial 1.
  # A parameter sharing a map level with a padding cell starts at the mean over
  # that level (TMB:::updateMap is tapply(..., mean)), so a 1 there pulls an
  # estimated M1 toward 1.0 per year; this keeps the mean exact.
  data("GOA2018SS", package = "Rceattle", envir = environment())
  d <- Rceattle::switch_check(GOA2018SS)
  pars <- suppressWarnings(suppressMessages(Rceattle::build_params(d)))
  row1 <- as.numeric(d$M1_base[1, seq_len(d$nages[1]) + 2])

  # Real cells come from M1_base, on the log scale.
  testthat::expect_equal(unname(pars$log_M1[1, 1, seq_len(d$nages[1])]),
                         log(row1))

  # Species 1 is one-sex, so sex index 2 is padding: it mirrors sex 1.
  testthat::expect_equal(unname(pars$log_M1[1, 2, ]),
                         unname(pars$log_M1[1, 1, ]))

  # Species 1 has 10 ages against max_age 21: ages 11+ take its oldest age.
  testthat::expect_equal(unname(pars$log_M1[1, 1, 11]),
                         log(row1[d$nages[1]]))
  # Not the array's initial 1, which is what diluted the start.
  testthat::expect_false(isTRUE(all.equal(unname(pars$log_M1[1, 1, 11]), 0)))

  # A two-sex species keeps its sexes distinct.
  testthat::expect_false(isTRUE(all.equal(unname(pars$log_M1[2, 1, 1]),
                                          unname(pars$log_M1[2, 2, 1]))))
})


testthat::test_that("the two M1 fill paths agree cell for cell", {
  # fit_mod(updateM1 = TRUE) refills log_M1 from M1_base. It must produce
  # exactly what build_params() does, or the same workbook gives two different
  # models depending on which path ran.
  testthat::skip_on_cran()
  data("GOA2018SS", package = "Rceattle", envir = environment())

  built <- suppressWarnings(suppressMessages(
    Rceattle::build_params(Rceattle::switch_check(GOA2018SS))
  ))

  fit <- suppressWarnings(suppressMessages(
    Rceattle::fit_mod(data_list = GOA2018SS,
                      inits = built,
                      M1Fun = Rceattle::build_M1(updateM1 = TRUE),
                      estimateMode = "DebugBuild",
                      msmMode = 0,
                      fit_control = Rceattle::fit_control(verbose = 0))
  ))

  testthat::expect_true(all(is.finite(fit$initial_params$log_M1)))
  testthat::expect_equal(fit$initial_params$log_M1, built$log_M1)
})


testthat::test_that("a missing M1_base row is refused, not fitted as M = 1", {
  # A (species, sex) the model has but M1_base does not name kept the array's
  # initial 1, i.e. a residual M of 1.0 per year. That is inside
  # build_bounds()'s [log(0.001), log(2)], so nothing downstream caught it.
  data("BS2017SS", package = "Rceattle", envir = environment())
  d <- BS2017SS
  d$M1_base <- d$M1_base[as.numeric(as.character(d$M1_base$Species)) != 2, ]

  testthat::expect_error(
    suppressWarnings(suppressMessages(
      Rceattle::build_params(Rceattle::switch_check(d))
    )),
    "M1_base has no row for"
  )
})


testthat::test_that("an M1 linkage may supply the level over a blank M1_base", {
  # The finiteness check runs AFTER the linkage initial-value pass for this
  # reason: an intercept `init` writes the level over every real age, so the
  # workbook is allowed to leave M1_base blank there. Checking at the fill
  # refused this, which is what an earlier draft did.
  testthat::skip_on_cran()
  data("BS2017SS", package = "Rceattle", envir = environment())
  d <- BS2017SS
  ac <- grep("^Age", names(d$M1_base))
  d$M1_base[1, ac[1:12]] <- NA

  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = "DebugBuild",
    msmMode = 0, random_rec = FALSE,
    M1Fun = Rceattle::build_M1(M1_model = 1, linkages = list(
      M1 = Rceattle::linkage_spec(formula = ~ 1, by = ~ species, species = 1L,
                                  init = list("(Intercept)" = 0.3)))),
    fit_control = Rceattle::fit_control(verbose = 0))))

  testthat::expect_true(all(is.finite(fit$initial_params$log_M1)))
  # The linkage's init is the level, on the log scale.
  testthat::expect_equal(unname(exp(fit$initial_params$log_M1[1, 1, 1])), 0.3,
                         tolerance = 1e-8)
})


testthat::test_that("updateM1 with an estimated M1 gives a finite objective", {
  # The padding sex cells take a map index under M1_model >= 1, so the -Inf the
  # old updateM1 fill left there reached the AD tape: GOA2018SS returned NaN
  # for M1_model 1, 2 and 3. Pinned loosely -- the point is finiteness, not the
  # value, which the golden references do not cover for this path.
  testthat::skip_on_cran()
  data("GOA2018SS", package = "Rceattle", envir = environment())

  for (mm in 1:3) {
    fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
      data_list = GOA2018SS, inits = NULL, file = NULL,
      estimateMode = "DebugBuild", msmMode = 0, niter = 3, random_rec = FALSE,
      M1Fun = Rceattle::build_M1(M1_model = mm, updateM1 = TRUE),
      fit_control = Rceattle::fit_control(verbose = 0))))
    testthat::expect_true(all(is.finite(fit$initial_params$log_M1)))
    testthat::expect_true(is.finite(fit$obj$fn(fit$obj$par)))
  }
})
