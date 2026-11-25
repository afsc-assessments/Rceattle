# =============================================================================
# init_dev padding in build_params().
#
# The cpp indexes init_dev as init_dev(sp, age - 1) over age = 1:(nages - 1) in
# every initMode -- the free-parameter branch, the equilibrium / non-equilibrium
# branches, and the init_dev penalty alike -- so columns nages:max_age are
# unused padding for ALL modes and carry the -999 sentinel.
#
# This was previously split in two, with the nages column gated on
# `data_list$initMode > 0`. By that point switch_check() has resolved initMode
# to its canonical *string*, so the test compared "FreeParams" > "0" -- TRUE by
# lexicographic collation -- and the gate was always open. These tests pin the
# invariant the gate was silently providing, so a future alias beginning with a
# digit (or any reintroduced string comparison) cannot quietly change it.
# =============================================================================

.idev <- function(mode, dat) {
  d <- dat
  d$initMode <- mode
  suppressMessages(suppressWarnings(Rceattle::build_params(d)))$init_dev
}

testthat::test_that("init_dev padding is identical for every initMode", {
  data(BS2017SS, package = "Rceattle")
  ref <- .idev(names(Rceattle:::initMode_map)[1], BS2017SS) == -999

  # Every string alias pads the same columns as the first ...
  for (m in names(Rceattle:::initMode_map)) {
    testthat::expect_identical(.idev(m, BS2017SS) == -999, ref, info = m)
  }
  # ... and so does the equivalent integer code, which switch_check() resolves.
  for (i in unname(Rceattle:::initMode_map)) {
    testthat::expect_identical(.idev(i, BS2017SS) == -999, ref, info = paste("int", i))
  }
})

testthat::test_that("init_dev starts at 0 except under FreeParams, and -999 from nages on", {
  data(BS2017SS, package = "Rceattle")
  modes <- names(Rceattle:::initMode_map)
  for (m in modes) {
    idev <- .idev(m, BS2017SS)
    for (sp in seq_len(BS2017SS$nspp)) {
      nages_sp <- BS2017SS$nages[sp]
      read <- idev[sp, seq_len(nages_sp - 1)]          # columns the cpp reads
      if (m == "FreeParams") {
        testthat::expect_true(all(is.finite(read) & read > 0), info = m)
      } else {
        testthat::expect_true(all(read == 0), info = m)
      }
      testthat::expect_true(all(idev[sp, nages_sp:ncol(idev)] == -999), info = m)
    }
  }
})

testthat::test_that("FreeParams starts at the unfished equilibrium the template computes", {
  # Under FreeParams init_dev is log numbers at age, so its start should be the
  # population initMode = "Equilibrium" builds from the same R0 and M1.
  data(BS2017SS, package = "Rceattle")
  q <- function(mode) suppressMessages(suppressWarnings(
    Rceattle::fit_mod(BS2017SS, estimateMode = 3, initMode = mode, msmMode = 0,
                      fit_control = fit_control(phase = FALSE, verbose = 0))))
  eq <- q("Equilibrium"); fp <- q("FreeParams")
  for (sp in seq_len(BS2017SS$nspp)) {
    ages <- 2:BS2017SS$nages[sp]                          # age indices the start sets
    testthat::expect_equal(fp$quantities$N_at_age[sp, 1, ages, 1],
                           eq$quantities$N_at_age[sp, 1, ages, 1], tolerance = 1e-12)
  }
})

testthat::test_that("FreeParams keeps a 0 start under a single-species stock-recruit curve", {
  # There R_init comes from alpha and beta, not the R0 the equilibrium start reads.
  data(BS2017SS, package = "Rceattle")
  fp <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    BS2017SS, estimateMode = 3, initMode = "FreeParams", msmMode = 0,
    recFun = build_srr(srr_fun = 2), fit_control = fit_control(phase = FALSE, verbose = 0))))
  idev <- fp$estimated_params$init_dev
  for (sp in seq_len(BS2017SS$nspp)) {
    testthat::expect_true(all(idev[sp, seq_len(BS2017SS$nages[sp] - 1)] == 0))
  }
})

testthat::test_that("FreeParams start reads M1 after an M1 linkage init sets it", {
  # M1_base may be blank when an M1 linkage intercept `init` supplies the level.
  data(BS2017SS, package = "Rceattle")
  d <- BS2017SS
  d$M1_base[1, grep("^Age", names(d$M1_base))] <- NA
  m1f <- build_M1(M1_model = 1, linkages = list(M1 = linkage_spec(
    ~ 1, by = ~ species, species = 1L, init = list(`(Intercept)` = 0.3))))
  q <- function(mode) suppressMessages(suppressWarnings(
    Rceattle::fit_mod(d, estimateMode = 3, initMode = mode, msmMode = 0, M1Fun = m1f,
                      fit_control = fit_control(phase = FALSE, verbose = 0))))
  eq <- q("Equilibrium"); fp <- q("FreeParams")
  testthat::expect_true(is.finite(fp$quantities$jnll))
  ages <- 2:d$nages[1]
  testthat::expect_equal(fp$quantities$N_at_age[1, 1, ages, 1],
                         eq$quantities$N_at_age[1, 1, ages, 1], tolerance = 1e-12)
})

testthat::test_that("padding respects ragged nages in a multispecies model", {
  data(BS2017MS, package = "Rceattle")
  idev <- .idev("FreeParams", BS2017MS)
  max_age <- ncol(idev)
  # BS2017MS carries species of differing nages, so the -999 boundary must move
  # per species rather than sitting at a single column.
  for (sp in seq_len(BS2017MS$nspp)) {
    nages_sp <- BS2017MS$nages[sp]
    testthat::expect_true(all(idev[sp, seq_len(nages_sp - 1)] > 0))
    testthat::expect_true(all(idev[sp, nages_sp:max_age] == -999))
  }
})
