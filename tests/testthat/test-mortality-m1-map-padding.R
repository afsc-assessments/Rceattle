# log_M1 is dimensioned to the WIDEST species -- [nspp, max(nsex), max(nages)]
# -- so a one-sex species in a model whose max(nsex) is 2 leaves padding
# sex cells.
# build_map_m1() used to index them: `map_list$log_M1[sp, , 1:nages_sp]` with a
# bare comma under M1_model 1 and 4, and `[sp, 2, ] <- [sp, 1, ]` under
# M1_model 3, put the padding in the SAME map level as the real cells.
#
# That matters because TMB does not take a shared parameter's start from the
# first cell of its level -- `TMB:::updateMap()` is
#
#     tapply(parameter.entry, map.entry, mean)
#
# so the start is the MEAN over the level, and a padding cell sitting at
# log(1) = 0 drags it toward a residual M of 1.0 per year.
#
# Measured on GOA2018SS (nsex c(1, 2, 1), nages c(10, 21, 12)) at M1_model = 1:
# pollock started at 0.637454 against an M1_base geometric mean of 0.406348
# (+57%, and exactly sqrt(0.406348)), cod at 0.710571 against 0.504911 (+41%).
# Species 2 was correct throughout -- its two REAL sexes legitimately share one
# level under a sex-invariant model, so 0.264575 = sqrt(0.20 * 0.35) is right.
#
# It is a STARTING value -- the template reads `sex < nsex(sp)`, so the shared
# parameter's MLE is identified by the real sex -- but the optimizer does not
# always reach that MLE from a start displaced this far. Fitting GOA2018SS at
# M1_model = 1 cold with phase = FALSE and newtonsteps = 0, the two starts
# converge to different points: objective 12849.030601 vs 12838.508102, M1
# species 1 0.38983395 vs 0.28725155, SSB at endyr 2018 362,066 vs 330,402 mt
# (-8.75%). BOTH are WARN by convergence_diagnostics() -- max|gradient| 2.2e-03
# and 1.7e-03 against an OK tier of 1e-03 -- so neither is a converged fit, and
# the surface carries at least four local minima over this range. Under
# phase = TRUE the two become bit-identical on this dataset.
#
# Phasing is NOT a general escape, and the live case is the reason to say so:
# the GOA multispecies assessment runs three fits, and its headline
# multispecies one (`R/02_fit_models.R`, `ms_mod`) sets this M1_model with
# `phase = FALSE`. A warm start does not help either -- the dilution happens
# when MakeADFun averages the map level, so it applies to whatever log_M1 the
# inits hold.
#
# The reason to fix it is that the start is simply wrong, and which optimum an
# unphased fit reaches is then arbitrary.
#
# The four golden references all run M1_model = 0, where the whole array is
# mapped out, so none of them covers any of this.

m1_geometric_mean <- function(d, sp) {
  row <- which(as.numeric(as.character(d$M1_base$Species)) == sp)[1]
  exp(mean(log(as.numeric(d$M1_base[row, seq_len(d$nages[sp]) + 2]))))
}

m1_map_levels <- function(fit, nspp, max_sex, max_age) {
  mp <- fit$obj$env$map$log_M1
  dim(mp) <- c(nspp, max_sex, max_age)
  mp
}


testthat::test_that("a padding sex cell takes no log_M1 map index", {
  testthat::skip_on_cran()
  data("GOA2018SS", package = "Rceattle", envir = environment())
  d <- suppressMessages(Rceattle::switch_check(GOA2018SS))

  # All five fixed-effect models, not just the three sex/age ones: the
  # environmentally-driven 4 and 5 write log_M1 through the same indices and
  # were edited with them.
  expect_levels <- c(3L, 4L, 64L, 3L, 4L)
  for (mm in 1:5) {
    fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
      data_list = GOA2018SS, inits = NULL, file = NULL,
      estimateMode = "DebugBuild", msmMode = 0, niter = 3, random_rec = FALSE,
      M1Fun = Rceattle::build_M1(M1_model = mm),
      fit_control = Rceattle::fit_control(verbose = 0))))

    mp <- m1_map_levels(fit, d$nspp, max(d$nsex), max(d$nages))
    lv <- unique(as.integer(mp)[!is.na(as.integer(mp))])
    testthat::expect_equal(length(lv), expect_levels[mm],
      info = paste("M1_model", mm, "free log_M1 levels"))
    for (sp in seq_len(d$nspp)) {
      padding <- setdiff(seq_len(max(d$nsex)), seq_len(d$nsex[sp]))
      for (sx in padding) {
        testthat::expect_true(all(is.na(mp[sp, sx, ])),
          info = paste("M1_model", mm, "species", sp, "padding sex", sx))
      }
    }
  }
})


testthat::test_that("the shared M1 starts at its own M1_base, undiluted", {
  testthat::skip_on_cran()
  data("GOA2018SS", package = "Rceattle", envir = environment())
  d <- suppressMessages(Rceattle::switch_check(GOA2018SS))

  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = GOA2018SS, inits = NULL, file = NULL,
    estimateMode = "DebugBuild", msmMode = 0, niter = 3, random_rec = FALSE,
    M1Fun = Rceattle::build_M1(M1_model = 1),
    fit_control = Rceattle::fit_control(verbose = 0))))

  starts <- exp(fit$obj$par[names(fit$obj$par) == "log_M1"])
  testthat::expect_length(starts, 3L)

  # Species 1 and 3 are one-sex: the start is that species' own geometric mean.
  testthat::expect_equal(unname(starts[1]), m1_geometric_mean(d, 1),
                         tolerance = 1e-8)
  testthat::expect_equal(unname(starts[3]), m1_geometric_mean(d, 3),
                         tolerance = 1e-8)
  # Not the diluted values, which were exactly the square roots.
  testthat::expect_false(isTRUE(all.equal(unname(starts[1]),
                                          sqrt(m1_geometric_mean(d, 1)))))

  # Species 2 is two-sex, so one sex-invariant M1 legitimately averages its two
  # REAL sexes: sqrt(0.20 * 0.35).
  testthat::expect_equal(unname(starts[2]), sqrt(0.20 * 0.35), tolerance = 1e-8)
})


testthat::test_that("M1_model 2 starts each species at its own M1_base", {
  # NOTE the level COUNT here is 4 on both sides of the fix, because fit_mod()
  # downgrades M1_model 2 to 1 for a single-sex species before the map is built
  # (R/6-fit_mod.R, "sex-specific -> sex-invariant for 1-sex model"), so the
  # raw M1_model = 2 branch is never reached through a fit. What discriminates
  # is the start values below, which come through the downgraded branch. The
  # count is asserted only to pin that the downgrade keeps happening.
  testthat::skip_on_cran()
  data("GOA2018SS", package = "Rceattle", envir = environment())
  d <- suppressMessages(Rceattle::switch_check(GOA2018SS))

  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = GOA2018SS, inits = NULL, file = NULL,
    estimateMode = "DebugBuild", msmMode = 0, niter = 3, random_rec = FALSE,
    M1Fun = Rceattle::build_M1(M1_model = 2),
    fit_control = Rceattle::fit_control(verbose = 0))))

  # One level for each one-sex species, two for the two-sex one.
  testthat::expect_length(fit$obj$par[names(fit$obj$par) == "log_M1"], 4L)

  starts <- exp(fit$obj$par[names(fit$obj$par) == "log_M1"])
  testthat::expect_equal(unname(starts[1]), m1_geometric_mean(d, 1),
                         tolerance = 1e-8)
  testthat::expect_equal(sort(unname(starts[2:3])), c(0.20, 0.35),
                         tolerance = 1e-8)
  testthat::expect_equal(unname(starts[4]), m1_geometric_mean(d, 3),
                         tolerance = 1e-8)
})
