# Every M1 random-effect family in build_map_m1() nests its deviation writes
# inside an M1_model test. Where no arm fires, log_M1_dev stays entirely NA --
# but the standard deviation and the AR1 correlations were freed anyway,
# outside that guard, so the density evaluated against a vector of zeros.
#
# M1_dev_log_sd has no entry in build_bounds(), so it keeps the generic
# [-Inf, Inf]. The contribution is n * (log(2*pi)/2 + log sigma) with n the
# element count, so the objective is LINEAR AND UNBOUNDED BELOW in it: the
# gradient is exactly n, and driving log sigma down reduces the objective
# without limit. Measured on BS2017SS at the starting sigma = 1, jnll row 16
# ("M random effects"), with the element count per family:
#
#   M1_re 1/4 (age)        n = nages      11.03 / 11.03 / 19.30
#   M1_re 2/5 (year)       n = nyrs_hind  35.84 each
#   M1_re 3/6 (age x year) n = nages*nyrs 430.06 / 430.06 / 752.61
#
# So the cost is family-specific, not one number: `Rceattle-models/EBS
# pollock/2024/06-time-varying-M.R` records 56.06 = 61 * log(2*pi)/2, which is
# the year family on a 61-year hindcast. That script's configuration is
# build_M1(M1_model = "fixed", M1_re = "iid_year") -- and "fixed" is
# M1_model = 0, build_M1()'s own default.
#
# Two changes follow. M1_model = 0 holds log_M1 at the input schedule, so there
# is no free fixed effect for the deviations to be absorbed into and they are
# the most identifiable of any level -- that arm is now open rather than
# refused. Where no arm can fire at all (M1_model 3, 4, 5) the combination is
# refused: asking for time-varying M and getting constant M is a different
# model. Through fit_mod() that is 18 of the 42 (M1_model, M1_re) pairs --
# M1_model has six levels, 0 to 5 -- measured identically on GOA2018SS
# (nsex c(1, 2, 1)) and BS2017MS (all one-sex). A direct build_map() call on a
# single-sex species also refuses M1_model = 2 for all six families, which
# fit_mod() never reaches because it downgrades 2 to 1 there first.
#
# The refusal reads whether an arm fired rather than restating the supported
# pairs, so opening an arm legalises its combination with no second registry.
# It skips a species with estDynamics > 0, because build_map_fixed_natage()
# maps that species' whole M1 random-effect block out afterwards, so the free
# sd cannot arise and refusing would reject a correct model.

prep <- function(mm, re, est = 0) {
  data("GOA2018SS", package = "Rceattle", envir = environment())
  d <- GOA2018SS
  d$M1_model <- if (length(mm) == 1) rep(mm, 3) else mm
  d$M1_re    <- if (length(re) == 1) rep(re, 3) else re
  d <- suppressMessages(Rceattle::switch_check(d))
  d$estDynamics <- if (length(est) == 1) rep(est, d$nspp) else est
  # build_map() reads these per species; fit_mod() extends them first.
  d$suitMode <- rep(0, d$nspp)
  d$growth_model <- rep(0, d$nspp)
  d
}

m1_free <- function(d) {
  p <- suppressWarnings(suppressMessages(Rceattle::build_params(d)))
  m <- suppressWarnings(suppressMessages(Rceattle::build_map(d, p)))
  list(dev = sum(!is.na(m$mapList$log_M1_dev)),
       # Cells, not parameters: M1_dev_log_sd[sp,] writes one level into both
       # sex cells, so 6 cells on GOA2018SS are 3 free parameters.
       sd_cells = sum(!is.na(m$mapList$M1_dev_log_sd)),
       sd_pars = length(unique(stats::na.omit(
         as.integer(m$mapList$M1_dev_log_sd)))),
       rho_cells = sum(!is.na(m$mapList$M1_rho)))
}


testthat::test_that("a random effect that no arm can serve is refused", {
  testthat::skip_on_cran()
  # M1_model 3 (sex- and age-specific), 4 and 5 (environmentally driven) have
  # no deviation arm in any family.
  for (mm in 3:5) {
    for (re in 1:6) {
      testthat::expect_error(m1_free(prep(mm, re)),
        "is not implemented for M1_model",
        info = paste("M1_model", mm, "M1_re", re))
    }
  }
})


testthat::test_that("M1_model = 0 keeps its deviations", {
  # The case the EBS pollock script wanted, and build_M1()'s default level.
  # M1_base is fixed input, so the deviations are the only free description of
  # M's shape -- nothing can absorb them, which makes this the most
  # identifiable level rather than one to refuse. An earlier version of this
  # change refused all six.
  testthat::skip_on_cran()
  want_rho <- c(`1` = 0L, `2` = 0L, `3` = 0L, `4` = 6L, `5` = 6L, `6` = 12L)
  for (re in 1:6) {
    got <- m1_free(prep(0, re))
    lab <- paste("M1_model 0 M1_re", re)
    testthat::expect_gt(got$dev, 0L)
    testthat::expect_equal(got$sd_pars, 3L, info = lab)
    testthat::expect_equal(got$rho_cells, want_rho[[as.character(re)]],
                           info = lab)
  }
  # Identical wiring to the estimated sex-invariant level: the deviations do
  # not care whether the level above them is fixed or free.
  for (re in 1:6) {
    testthat::expect_equal(m1_free(prep(0, re))$dev,
                           m1_free(prep(1, re))$dev,
                           info = paste("M1_re", re))
  }
})


testthat::test_that("a canonical M1_re string frees the deviations it names", {
  # switch_check() now canonicalises M1_re to its integer code. Without that a
  # workbook string reached build_map_m1() as a character, matched no arm --
  # `"iid_year" %in% c(2, 5)` is FALSE -- and every deviation was mapped out
  # with no message, which is the same silent constant-M this file is about.
  testthat::skip_on_cran()
  testthat::expect_equal(m1_free(prep(1, "iid_year"))$dev,
                         m1_free(prep(1, 2))$dev)
  testthat::expect_equal(m1_free(prep(1, "ar1_age_year"))$rho_cells, 12L)
  # "none" is 0, and must not be read as "greater than zero" by a string
  # comparison -- "none" > 0 is TRUE lexicographically.
  none <- m1_free(prep(1, "none"))
  testthat::expect_equal(none$dev, 0L)
  testthat::expect_equal(none$sd_cells, 0L)
  testthat::expect_equal(none$rho_cells, 0L)
})


testthat::test_that("a fixed-numbers species is not refused", {
  # build_map_fixed_natage() maps log_M1_dev, M1_dev_log_sd and M1_rho out for
  # an estDynamics > 0 species, after build_map_m1() runs -- so the free sd
  # cannot arise and the refusal must not fire. A fixed-numbers predator is
  # the common multispecies setup.
  testthat::skip_on_cran()
  got <- m1_free(prep(c(1, 1, 3), 2, est = c(0, 0, 1)))
  testthat::expect_gt(got$dev, 0L)
  # Species 3's sd is mapped out by the later helper, so 4 cells not 6.
  testthat::expect_equal(got$sd_cells, 4L)
})


testthat::test_that("M1_re = 0 frees no deviation, sd or correlation", {
  # Passes either side of the change -- with M1_re = 0 no family block
  # executes, so the sd was never freed there. Here to pin that the refusal
  # does not fire on the default and that the no-random-effects case stays
  # exactly empty.
  testthat::skip_on_cran()
  for (mm in 0:5) {
    got <- m1_free(prep(mm, 0))
    lab <- paste("M1_model", mm, "M1_re 0")
    testthat::expect_equal(got$dev, 0L, info = lab)
    testthat::expect_equal(got$sd_cells, 0L, info = lab)
    testthat::expect_equal(got$rho_cells, 0L, info = lab)
  }
})


testthat::test_that("the refusal does not read a caller's log_M1_dev", {
  # The guard is `&&`, not `&`: the right-hand side subscripts log_M1_dev
  # [sp,,,], and build_map() is exported, so a caller passing a 3-D block with
  # no random effect requested must not be touched.
  testthat::skip_on_cran()
  d <- prep(1, 0)
  p <- suppressWarnings(suppressMessages(Rceattle::build_params(d)))
  p$log_M1_dev <- array(0, dim = c(d$nspp, max(d$nages), 2))
  testthat::expect_no_error(
    suppressWarnings(suppressMessages(Rceattle::build_map(d, p))))
})


testthat::test_that("M1_model = 0 deviations are informed by the data", {
  # The reason to open this arm rather than refuse it. ceattle.cpp's M1
  # assembly adds log_M1_dev unconditionally --
  #
  #   M1_at_age = exp(log_M1 + log_M1_dev + ... )
  #
  # with no M1_model gate -- so under the input-schedule level the template
  # already reads the deviations; only the map was withholding them. And
  # log_M1 is mapped out there, so nothing competes to absorb them: a freed
  # parameter the data cannot move would be worse than the refusal it
  # replaced.
  testthat::skip_on_cran()
  data("GOA2018SS", package = "Rceattle", envir = environment())
  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = GOA2018SS, inits = NULL, file = NULL,
    estimateMode = "DebugBuild", msmMode = 0, niter = 3, random_rec = FALSE,
    M1Fun = Rceattle::build_M1(M1_model = 0, M1_re = 2),
    fit_control = Rceattle::fit_control(verbose = 0))))

  p <- fit$obj$env$last.par.best
  dev <- which(names(p) == "log_M1_dev")
  testthat::expect_equal(length(dev), 126L)
  # The input schedule stays fixed, so the deviations carry M's shape alone.
  testthat::expect_equal(sum(names(p) == "log_M1"), 0L)

  # They reach the likelihood, and the data can move them.
  testthat::expect_true(any(fit$obj$gr(p)[dev] != 0))
  q <- p
  q[dev[1]] <- 0.25
  testthat::expect_false(isTRUE(all.equal(fit$obj$fn(p), fit$obj$fn(q))))
})


testthat::test_that("the refused set through fit_mod() is M1_model 3, 4, 5", {
  # fit_mod() is the path a user takes, and it extends and downgrades the
  # switches first -- so the reachable refused set is not the same as the one a
  # direct build_map() call sees. Measured on both a ragged-nsex and an
  # all-one-sex dataset: identical, and no other error anywhere in the grid.
  testthat::skip_on_cran()
  for (ds in c("GOA2018SS", "BS2017MS")) {
    e <- new.env()
    utils::data(list = ds, package = "Rceattle", envir = e)
    dl <- get(ds, envir = e)
    refused <- character()
    other <- character()
    for (mm in 0:5) {
      for (re in 0:6) {
        msg <- tryCatch({
          suppressWarnings(suppressMessages(Rceattle::fit_mod(
            data_list = dl, inits = NULL, file = NULL,
            estimateMode = "DebugBuild", msmMode = 0, niter = 3,
            random_rec = FALSE,
            M1Fun = Rceattle::build_M1(M1_model = mm, M1_re = re),
            fit_control = Rceattle::fit_control(verbose = 0))))
          NULL
        }, error = function(e) conditionMessage(e))
        key <- paste0(mm, "x", re)
        if (!is.null(msg)) {
          if (grepl("is not implemented for M1_model", msg)) {
            refused <- c(refused, key)
          } else {
            other <- c(other, paste0(key, ": ", msg))
          }
        }
      }
    }
    want <- as.vector(outer(3:5, 1:6, function(a, b) paste0(a, "x", b)))
    testthat::expect_setequal(refused, want)
    testthat::expect_equal(length(refused), 18L, info = ds)
    # Nothing in the grid may fail for any OTHER reason.
    testthat::expect_equal(other, character(0), info = ds)
  }
})


testthat::test_that("a factor or numeric-string switch still resolves", {
  # switch_check() canonicalises M1_model and M1_re through .map_switch(), the
  # resolver every other per-species switch here uses -- not build_M1()'s
  # stricter .coerce_M1_arg(), which refuses a factor, a numeric-looking "1"
  # and an NA. read.csv(stringsAsFactors = TRUE) produces factors, which is
  # why test-switches-map-switch-factor.R exists.
  testthat::skip_on_cran()
  data("GOA2018SS", package = "Rceattle", envir = environment())
  resolve <- function(v) {
    d <- GOA2018SS
    d$M1_model <- v
    d$M1_re <- rep(0, 3)
    suppressMessages(Rceattle::switch_check(d))$M1_model
  }
  testthat::expect_equal(as.numeric(resolve(c(1L, 1L, 1L))), c(1, 1, 1))
  testthat::expect_equal(as.numeric(resolve(c("1", "1", "1"))), c(1, 1, 1))
  testthat::expect_equal(as.numeric(resolve(factor(c(1, 1, 1)))), c(1, 1, 1))
  testthat::expect_equal(as.numeric(resolve(rep("sex_age_invariant", 3))),
                         c(1, 1, 1))
  # An NA passes through rather than erroring, as for every other switch.
  testthat::expect_true(is.na(as.numeric(resolve(c(1L, NA, 1L)))[2]))
})


testthat::test_that("no family writes a deviation into a padding sex", {
  # M1_re 1/4's M1_model == 2 arm lacked the `nsex_sp == 2 &` guard that 2/5
  # and 3/6 carry, so a single-sex species took deviations in its padding sex
  # slot -- free and Laplace-integrated, but scored by no density and read by
  # no cell. All four families now agree, so the whole configuration refuses
  # on a single-sex species rather than half of it freeing dead parameters.
  testthat::skip_on_cran()
  data("GOA2018SS", package = "Rceattle", envir = environment())
  for (re in c(1, 4)) {
    d <- GOA2018SS
    d$M1_model <- rep(2, 3)
    d$M1_re <- rep(re, 3)
    d <- suppressMessages(Rceattle::switch_check(d))
    d$suitMode <- rep(0, d$nspp)
    d$growth_model <- rep(0, d$nspp)
    testthat::expect_error(
      suppressWarnings(suppressMessages({
        p <- Rceattle::build_params(d)
        Rceattle::build_map(d, p)
      })),
      "is not implemented for M1_model", info = paste("M1_re", re))
  }
  # On an all-two-sex model the same configuration is fine, and nothing lands
  # in a padding slot because there is none.
  d2 <- GOA2018SS
  d2$M1_model <- rep(2, 3)
  d2$M1_re <- rep(1, 3)
  d2 <- suppressMessages(Rceattle::switch_check(d2))
  d2$nsex <- rep(2, 3)
  d2$suitMode <- rep(0, d2$nspp)
  d2$growth_model <- rep(0, d2$nspp)
  p2 <- suppressWarnings(suppressMessages(Rceattle::build_params(d2)))
  mp <- suppressWarnings(suppressMessages(
    Rceattle::build_map(d2, p2)))$mapList$log_M1_dev
  testthat::expect_gt(sum(!is.na(mp)), 0L)
})
