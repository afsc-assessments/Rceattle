# build_map_growth() assigned log_growth_pars a FIXED STRIDE of 4 per species,
# but a two-sex species needs 8 slots: females took `(sp-1)*4 + 1:3` and males
# `(sp-1)*4 + 5:7`. So species sp's MALE block was species sp+1's FEMALE block
# -- `5,6,7` for both sp1-males and sp2-females -- the same TMB parameter,
# estimated once for both, with 3 (von Bertalanffy) or 4 (Richards) free
# parameters silently gone. Richards had it too, at `1:4` / `5:8`.
#
# Measured collisions under the old arithmetic, growth on for every species:
#
#   nsex c(1,2,1)   9 distinct of 12   3 parameters lost
#   nsex c(2,2)     9 of 12            3 lost
#   nsex c(2,1)     6 of 9             3 lost
#   nsex c(2,2,2)  12 of 18            6 lost
#   nsex c(1,2)     9 of 9             none: a TRAILING two-sex species is safe
#
# Not reachable by anything shipping: no bundled dataset sets `growth_model`
# (the schema default is 0), all four golden references run 0, and every
# multispecies sibling config uses `build_growth(fun = "empirical")`.
#
# The suite DOES fit multispecies parametric growth -- `test-dynamics-brps.R`
# and `test-composition-dm-diet-weight.R` both do, on a `make_msm_test_data()`
# fixture. They are immune because that helper sets `nsex <- rep(1, nspp)`
# (`helpers-make-msm-data.R:332`): the collision needs a species with two
# sexes, and a trailing one is harmless. So what hid this was single SEX, not
# single species. It arms the moment parametric growth goes on a multispecies
# model with a non-trailing two-sex species, which is the shape of the GOA
# multispecies assessment (arrowtooth is two-sex at species 2 of 3).
#
# The fix is a running counter, as build_map_m1() already uses. These tests
# drive build_map() directly rather than through a fit, because reaching the
# branch through fit_mod() needs parametric growth on a ragged multispecies
# model and the point is the index arithmetic.

growth_levels <- function(mp, nsex, npar) {
  out <- integer(0)
  for (sp in seq_along(nsex)) {
    for (sx in seq_len(nsex[sp])) {
      v <- mp[sp, sx, seq_len(npar)]
      out <- c(out, v[!is.na(v)])
    }
  }
  out
}


testthat::test_that("no growth parameter is shared across a species boundary", {
  testthat::skip_on_cran()
  data("BS2017MS", package = "Rceattle", envir = environment())

  for (gm in c(1, 2)) {
    npar <- if (gm == 1) 3L else 4L
    for (nsex_case in list(c(1, 2, 1), c(2, 2, 2), c(2, 1, 1))) {
      d <- suppressMessages(Rceattle::switch_check(BS2017MS))
      d$nsex <- nsex_case
      d$growth_model <- rep(gm, length(nsex_case))
      # build_map() indexes these per species; fit_mod() extends them to nspp
      # before calling it, so a direct call has to do the same.
      d$suitMode <- rep(0, length(nsex_case))

      pars <- suppressWarnings(suppressMessages(Rceattle::build_params(d)))
      mp <- suppressWarnings(suppressMessages(
        Rceattle::build_map(d, pars)))$mapList$log_growth_pars

      lab <- paste("growth_model", gm, "nsex",
                   paste(nsex_case, collapse = ","))
      lv <- growth_levels(mp, nsex_case, npar)
      # One distinct level per real (species, sex) block: nothing shared.
      testthat::expect_equal(length(unique(lv)), length(lv),
        info = lab)
      testthat::expect_equal(length(lv), npar * sum(nsex_case), info = lab)
    }
  }
})


testthat::test_that("a padding sex takes no growth parameter", {
  # Pins an invariant that already held -- the old code also wrote `[sp, 2, ]`
  # only under `nsex_sp == 2` -- so this block passes either side of the stride
  # fix. It is here because the running counter is the kind of change that
  # could start writing a padding cell by accident.
  testthat::skip_on_cran()
  data("BS2017MS", package = "Rceattle", envir = environment())
  d <- suppressMessages(Rceattle::switch_check(BS2017MS))
  d$nsex <- c(1, 2, 1)
  d$growth_model <- rep(1, 3)
  d$suitMode <- rep(0, 3)

  pars <- suppressWarnings(suppressMessages(Rceattle::build_params(d)))
  mp <- suppressWarnings(suppressMessages(
    Rceattle::build_map(d, pars)))$mapList$log_growth_pars

  # Species 1 and 3 are one-sex, so sex index 2 is padding.
  testthat::expect_true(all(is.na(mp[1, 2, ])))
  testthat::expect_true(all(is.na(mp[3, 2, ])))
  # Slots past the model's own parameter count are never estimated.
  testthat::expect_true(all(is.na(mp[, , 4])))
})


testthat::test_that("CAAL endpoints survive ragged bin counts per species", {
  # log_growth_pars' L1 and L-infinity start at the smallest and largest length
  # with CAAL data. Taking them by pivoting on the per-species bin ORDINAL
  # named the columns Bin1..Bin<n>, so a 3-bin and a 5-bin species produced
  # `Bin1, Bin3, Bin5` -- three value columns for a two-column target, with NA
  # for whichever species lacked that ordinal. build_params() then died with
  # "number of items to replace is not a multiple of replacement length".
  testthat::skip_on_cran()
  data("BS2017MS", package = "Rceattle", envir = environment())
  d <- suppressMessages(Rceattle::switch_check(BS2017MS))
  d$growth_model <- rep(1, d$nspp)

  # Ragged on purpose: 3 lengths for species 1, 5 for species 2.
  d$caal_data <- data.frame(
    Species = c(rep(1, 3), rep(2, 5)),
    Length  = c(10, 20, 30, 12, 18, 24, 30, 36),
    Year = 1979L, Sex = 0L, Fleet_code = 1L, Age = 1L, Sample_size = 1,
    stringsAsFactors = FALSE)

  pars <- suppressWarnings(suppressMessages(Rceattle::build_params(d)))

  # Each species' endpoints are its OWN min and max, not another's.
  gp <- pars$log_growth_pars
  testthat::expect_equal(unname(exp(gp[1, 1, 2:3])), c(10, 30))
  testthat::expect_equal(unname(exp(gp[2, 1, 2:3])), c(12, 36))
  testthat::expect_true(all(is.finite(gp[1:2, 1, 2:3])))
})


testthat::test_that("a two-sex multispecies fit keeps all its growth pars", {
  # One level up from build_map(): a configuration a user can construct, fitted
  # through fit_mod(). The old stride fused species 1's males with species 2's
  # females, so obj$par was 3 short -- 132 against 135 -- and the objective at
  # the STARTING values was identical either way, because this fixture gives
  # every species the same K/L1/Linf and the mean-collapse is then a no-op.
  # That identical start is the mechanism of the silence: the fused parameter
  # only diverges once the optimizer moves, so a test comparing obj$fn() at the
  # start would see nothing. Pin the parameter count instead.
  testthat::skip_on_cran()
  sim <- make_msm_test_data(nspp = 2, years = 1:20)
  dd <- sim$data_list
  dd$nsex <- c(2, 1)

  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = dd, estimateMode = "DebugBuild", msmMode = 0,
    random_rec = FALSE,
    growthFun = Rceattle::build_growth(fun = "vonBertalanffy"),
    fit_control = Rceattle::fit_control(verbose = 0))))

  mp <- fit$obj$env$map$log_growth_pars
  lv <- unique(as.integer(mp)[!is.na(as.integer(mp))])
  # 3 parameters for each of sp1-female, sp1-male, sp2-female.
  testthat::expect_equal(length(lv), 9L)
  testthat::expect_equal(length(fit$obj$par), 135L)
  testthat::expect_true(is.finite(fit$obj$fn(fit$obj$par)))
})


testthat::test_that("an unusable CAAL length range is refused", {
  # L1 and L-infinity start at the smallest and largest length with CAAL data,
  # so a species with one distinct length starts them equal. length_sd_at_age()
  # interpolates the length-at-age SD as
  # sd0 + (sd1 - sd0) / (linf - l1) * (len - l1), and build_params() starts
  # both growth_log_sd at 0, so that is 0/0 = NaN for every age above age_L1
  # and the NaN reaches the age-length key and the likelihood. Refused rather
  # than warned: a single length bin cannot inform a growth curve.
  testthat::skip_on_cran()
  data("whamGrowthData", package = "Rceattle", envir = environment())
  d <- suppressMessages(Rceattle::switch_check(whamGrowthData))

  # As shipped: 65 distinct lengths, so it builds.
  testthat::expect_no_error(
    suppressWarnings(suppressMessages(Rceattle::build_params(d))))

  flat <- d
  flat$caal_data$Length <- 50
  testthat::expect_error(
    suppressWarnings(suppressMessages(Rceattle::build_params(flat))),
    "one distinct CAAL length")

  # A blank Length is tolerated while the species keeps a usable range...
  one_na <- d
  one_na$caal_data$Length[1] <- NA
  testthat::expect_no_error(
    suppressWarnings(suppressMessages(Rceattle::build_params(one_na))))

  # ...but a species with no finite length at all is refused by name, rather
  # than starting L1 at NA.
  all_na <- d
  all_na$caal_data$Length <- NA_real_
  testthat::expect_error(
    suppressWarnings(suppressMessages(Rceattle::build_params(all_na))),
    "no usable length range for species 1")
})
