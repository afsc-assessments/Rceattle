# The per-sex apical selectivity offset (5.38.0): a selectivity linkage on
# `apical` scales one sex's whole curve by exp(log_sel_apical), after the form
# and before the shared normalizer. log_F is shared by the sexes, so only the
# male:female ratio is identified; the checks below pin what the model refuses.
#
# GOAatf is the only bundled two-sex fit that converges (inst/dev/TRAPS.md);
# estimateMode = 3 builds and evaluates without fitting.

apical_data <- function(flt = 3L, form = "Logistic", scope = "AcrossSexes",
                        norm_bin = "Max") {
  d <- Rceattle::GOAatf
  d$fleet_control$Selectivity[flt]       <- form
  d$fleet_control$Sel_norm_scope[flt]    <- scope
  d$fleet_control$Sel_norm_bin[flt]      <- norm_bin
  # GOAatf ships Sel_curve_pen2 = 200 on this fleet. On the AR1 forms that slot
  # is a correlation, not a weight, and data_check() stops on a magnitude above
  # 10 -- with "AR1" in the message, which is the string the AR1 refusal below
  # matches. Zeroing it makes that test fail for its own reason.
  d$fleet_control$Sel_curve_pen2[flt]    <- 0
  d
}

apical_build <- function(d, spec = NULL, inits = NULL, map = NULL) {
  suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = d, inits = inits, map = map, msmMode = 0, estimateMode = 3,
    selFun = spec,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))
}

male_offset <- function(flt = 3L, ...) {
  Rceattle::build_selectivity(linkages = list(
    apical = Rceattle::linkage_spec(~ 1, by = ~ fleet + sex, fleet = flt,
                                    sex = "male", ...)))
}

testthat::test_that("the encoder, the accumulator and the prior re-target agree on apical", {
  testthat::expect_identical(Rceattle:::LINKAGE_PARAM_CODES$sel[["apical"]], 5L)
  testthat::expect_true("apical" %in% Rceattle:::SEL_LINKAGE_PARAMS)
  src <- testthat::test_path("..", "..", "src", "TMB")
  testthat::skip_if(!dir.exists(src), "src/TMB not available")
  lk  <- paste(readLines(file.path(src, "linkage.hpp"), warn = FALSE), collapse = "\n")
  cpp <- paste(readLines(file.path(src, "ceattle.cpp"), warn = FALSE), collapse = "\n")
  testthat::expect_match(lk,  "param == 5", fixed = TRUE)
  testthat::expect_match(cpp, "linkage_param(i) == 5", fixed = TRUE)
  testthat::expect_match(cpp, "else if (param == 5)", fixed = TRUE)
})

testthat::test_that("an apical linkage frees the named sex's cell and changes nothing at its start", {
  testthat::skip_on_cran()
  d  <- apical_data()
  m0 <- apical_build(d)
  m  <- apical_build(d, male_offset())

  testthat::expect_true(all(is.na(m0$map$mapList$log_sel_apical)))
  testthat::expect_identical(sum(names(m$obj$par) == "log_sel_apical"), 1L)
  testthat::expect_true(is.na(m$map$mapList$log_sel_apical[3, 1]))
  testthat::expect_false(is.na(m$map$mapList$log_sel_apical[3, 2]))
  # exp(0) = 1: the multiplier is on the tape but the fit is the same number.
  testthat::expect_equal(m$obj$fn(), m0$obj$fn())
  testthat::expect_equal(m$quantities$sel_at_age, m0$quantities$sel_at_age)
})

testthat::test_that("the offset scales the named sex and survives across-sex and no normalization", {
  testthat::skip_on_cran()
  for (norm_bin in c("Max", "Off")) {
    d <- apical_data(norm_bin = norm_bin)
    m <- apical_build(d, male_offset())
    inits <- m$initial_params
    inits$log_sel_apical[3, 2] <- log(0.5)
    m2 <- apical_build(d, male_offset(), inits = inits)
    s0 <- m$quantities$sel_at_age[3, , , 1]
    s  <- m2$quantities$sel_at_age[3, , , 1]
    # Same logistic shape at the start, so the male curve is half the female one.
    testthat::expect_equal(s[1, ], s0[1, ])
    testthat::expect_equal(s[2, ], 0.5 * s0[2, ])
    testthat::expect_equal(max(s[2, ]) / max(s[1, ]), 0.5)
  }
})

testthat::test_that("an apical linkage the model cannot identify is refused", {
  testthat::skip_on_cran()
  d <- apical_data()

  # No sex named: both sexes would carry it.
  no_sex <- Rceattle::build_selectivity(linkages = list(
    apical = Rceattle::linkage_spec(~ 1, by = ~ fleet, fleet = 3L)))
  testthat::expect_error(apical_build(d, no_sex), "names no sex")

  # Both sexes named through the strata.
  both <- Rceattle::build_selectivity(linkages = list(
    apical = Rceattle::linkage_spec(~ 1, by = ~ fleet + sex, fleet = 3L)))
  testthat::expect_error(apical_build(d, both), "both sexes")

  # Both sexes named across two rows whose fleet keys differ: the same ridge,
  # reached through the union of named sexes per fleet.
  two_rows <- Rceattle::build_selectivity(linkages = list(
    apical = Rceattle::linkage_spec(~ 1, by = ~ fleet + sex, fleet = 3L, sex = "male"),
    apical = Rceattle::linkage_spec(~ 1, by = ~ fleet + sex, fleet = 3L, sex = "female")))
  testthat::expect_error(apical_build(d, two_rows), "both sexes")

  # No fleet named: one cell per fleet would be freed under a single prior.
  no_fleet <- Rceattle::build_selectivity(linkages = list(
    apical = Rceattle::linkage_spec(~ 1, by = ~ sex, sex = "male")))
  testthat::expect_error(apical_build(d, no_fleet), "names no fleet")

  # An identity-link offset below -1 would make selectivity negative.
  ident <- Rceattle::build_selectivity(linkages = list(
    apical = Rceattle::linkage_spec(~ 1, by = ~ fleet + sex, fleet = 3L, sex = "male",
                                    link = "identity")))
  testthat::expect_error(apical_build(d, ident), "identity")

  # Within-sex normalization divides the offset out.
  testthat::expect_error(apical_build(apical_data(scope = "WithinSex"), male_offset()),
                         "WithinSex")

  # The AR1 forms estimate a per-sex level in sel_coff already. Through
  # fit_mod(), so the refusal is proven on the path a user takes rather than
  # against the helper alone. Matched on the apical clause, not on "AR1":
  # data_check() also says "AR1" when it rejects Sel_curve_pen2 as an
  # out-of-range correlation on these forms, and that message would satisfy a
  # looser pattern whether or not this refusal existed.
  testthat::expect_error(apical_build(apical_data(form = "2DAR1"), male_offset()),
                         "apical selectivity linkage on fleet\\(s\\) .* with an AR1")

  # A one-sex species: the offset is the common level log_F already carries.
  d1 <- Rceattle::GOApollock
  d1$fleet_control$Time_varying_sel[8]  <- "Off"
  one <- Rceattle::build_selectivity(linkages = list(
    apical = Rceattle::linkage_spec(~ 1, by = ~ fleet + sex, fleet = 8L, sex = 1L)))
  testthat::expect_error(apical_build(d1, one), "one-sex")

  # A Fixed (input) curve has no height to offset; checked at the table level
  # because data_check() asks for emp_sel first.
  fc <- d$fleet_control
  fc$Selectivity[3] <- "Fixed"
  ap <- data.frame(process = "sel", param = "apical", fleet = 3L, sex = 2L,
                   stringsAsFactors = FALSE)
  testthat::expect_error(Rceattle:::.check_sel_apical_rows(ap, fc, d$nsex), "Fixed")
})

testthat::test_that("mirrored fleets share one apical offset", {
  testthat::skip_on_cran()
  d <- Rceattle::GOAatf
  # GOAatf ships the length-comp survey (fleet 2) as a mirror of the bottom
  # trawl (fleet 1, Selectivity_index 1), both Logistic, and no Sel_norm_scope
  # column (switch_check() fills the default).
  d$fleet_control$Sel_norm_scope <- NA
  for (f in 1:2) {
    d$fleet_control$Sel_norm_scope[f] <- "AcrossSexes"
    d$fleet_control$Sel_norm_bin[f]   <- "Max"
  }
  m <- apical_build(d, male_offset(flt = 1L))
  testthat::expect_identical(m$map$mapList$log_sel_apical[2, ],
                             m$map$mapList$log_sel_apical[1, ])
  testthat::expect_identical(sum(names(m$obj$par) == "log_sel_apical"), 1L)
  # On the mirror itself the block is the lead's, so the linkage is refused.
  testthat::expect_error(apical_build(d, male_offset(flt = 2L)), "lead fleet")
})

# Selectivity_index is a group key, not a fleet code. The lead is the group's
# first fleet that is not "Off", and a group of one shares nothing -- so the
# refusal must follow .shared_block_lead(), not `Selectivity_index != Fleet_code`.
# Reading it as a fleet code refused GOAatf fleet 3 (a group of one keyed 2) and,
# worse, passed the follower of a group keyed on its second member while
# adjust_map_shared_params() then mapped the offset off, pinning it at 1 unscored.
testthat::test_that("the apical refusal follows the group lead, not the index value", {
  testthat::skip_on_cran()

  # A group of one whose key is not its own Fleet_code owns its block.
  d <- Rceattle::GOAatf
  testthat::expect_identical(d$fleet_control$Selectivity_index[3], 2)
  testthat::expect_silent(
    Rceattle:::.check_sel_apical_rows(
      data.frame(fleet = 3L, sex = 2L, link = "log", est_phase = NA_integer_),
      suppressMessages(Rceattle::switch_check(d))$fleet_control, nsex = d$nsex))

  # A group of two keyed on its SECOND member: fleet 1 leads, fleet 2 follows.
  d2 <- Rceattle::GOAatf
  d2$fleet_control$Selectivity_index[1:2] <- 2
  fc <- suppressMessages(Rceattle::switch_check(d2))$fleet_control
  row <- function(f) data.frame(fleet = f, sex = 2L, link = "log", est_phase = NA_integer_)
  testthat::expect_silent(Rceattle:::.check_sel_apical_rows(row(1L), fc, nsex = d2$nsex))
  testthat::expect_error(Rceattle:::.check_sel_apical_rows(row(2L), fc, nsex = d2$nsex),
                         "lead fleet")
})

testthat::test_that("an apical offset on an Off fleet is refused", {
  testthat::skip_on_cran()
  # Nothing is fit to an "Off" fleet, so the offset is a flat direction: the
  # Hessian is singular and getsd fails with nothing naming the cause.
  d <- Rceattle::GOAatf
  d$fleet_control$Fleet_type[3] <- "Off"
  fc <- suppressMessages(Rceattle::switch_check(d))$fleet_control
  testthat::expect_error(
    Rceattle:::.check_sel_apical_rows(
      data.frame(fleet = 3L, sex = 2L, link = "log", est_phase = NA_integer_),
      fc, nsex = d$nsex),
    "Fleet_type")
})

testthat::test_that("inits and a stored map that predate log_sel_apical still refit", {
  testthat::skip_on_cran()
  d  <- apical_data()
  m0 <- apical_build(d)
  inits <- m0$estimated_params; inits$log_sel_apical <- NULL
  map   <- m0$map
  map$mapList$log_sel_apical <- NULL; map$mapFactor$log_sel_apical <- NULL
  testthat::expect_message(
    m1 <- suppressWarnings(Rceattle::fit_mod(
      data_list = d, inits = inits, map = map, msmMode = 0, estimateMode = 3,
      fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))),
    "predate")
  testthat::expect_equal(m1$obj$fn(), m0$obj$fn())
  testthat::expect_true(all(is.na(m1$map$mapList$log_sel_apical)))
})

testthat::test_that("`by` defaults to fleet + sex, and an uninformed offset warns", {
  testthat::skip_on_cran()
  d <- apical_data()
  # `by` omitted: the apical default is ~ fleet + sex, so naming the fleet and
  # the sex is the whole specification.
  auto <- Rceattle::build_selectivity(linkages = list(
    apical = Rceattle::linkage_spec(~ 1, fleet = 3L, sex = "male")))
  m <- apical_build(d, auto)
  testthat::expect_identical(sum(names(m$obj$par) == "log_sel_apical"), 1L)
  testthat::expect_false(is.na(m$map$mapList$log_sel_apical[3, 2]))

  # No joint-sex compositions on that fleet and no prior: nothing informs the
  # sexes' ratio, and the fit would still report a number.
  d2 <- d
  d2$comp_data$Sex[d2$comp_data$Fleet_code == 3] <- 1L
  bare <- function(dat, spec) suppressMessages(Rceattle::fit_mod(
    data_list = dat, msmMode = 0, estimateMode = 3, selFun = spec,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0)))
  testthat::expect_warning(bare(d2, male_offset()), "no joint-sex composition")
  # A prior is the documented answer, so it is silent.
  testthat::expect_no_warning(
    bare(d2, male_offset(priors = list(intercept = lognormal(0, 0.5)))),
    message = "no joint-sex composition")
  # With joint comps it is silent too.
  testthat::expect_no_warning(bare(d, male_offset()), message = "no joint-sex composition")
  # And a refit is silent: .refit_like() sets quiet_data_check, and the first
  # fit already raised it. The refusals stay on either way.
  testthat::expect_no_warning(
    suppressMessages(Rceattle::fit_mod(
      data_list = d2, msmMode = 0, estimateMode = 3, selFun = male_offset(),
      quiet_data_check = TRUE,
      fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))),
    message = "no joint-sex composition")
})


testthat::test_that("a refit through .refit_like() keeps the apical linkage", {
  testthat::skip_on_cran()
  # retrospective(), jitter(), self_test(), profile() and run_mse() all refit
  # through .refit_like(), which rebuilds selFun from data_list$sel_linkages.
  d <- apical_data()
  m <- apical_build(d, male_offset())
  rf <- suppressMessages(suppressWarnings(Rceattle:::.refit_like(
    data_list = m$data_list, inits = m$estimated_params, estimateMode = 3,
    map = m$map)))
  testthat::expect_identical(sum(rf$data_list$linkage_table$param == "apical"), 1L)
  testthat::expect_identical(sum(names(rf$obj$par) == "log_sel_apical"), 1L)
})


testthat::test_that("a fitted apical offset moves the male:female ratio off 1", {
  testthat::skip_on_cran()
  d <- apical_data()
  spec <- male_offset(priors = list(intercept = lognormal(0, 0.5)))
  fit <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = d, msmMode = 0, estimateMode = "Hindcast", selFun = spec,
    fit_control = Rceattle::fit_control(phase = TRUE, getsd = FALSE, verbose = 0))))
  s <- fit$quantities$sel_at_age[3, , , 1]
  ratio <- max(s[2, ]) / max(s[1, ])
  testthat::expect_true(is.finite(fit$opt$objective))
  testthat::expect_lt(max(abs(fit$obj$gr())), 0.1)
  testthat::expect_false(isTRUE(all.equal(ratio, 1)))
  testthat::expect_false(isTRUE(all.equal(fit$estimated_params$log_sel_apical[3, 2], 0)))
  # AcrossSexes: one pooled reference, so the more-selected sex peaks at 1 and
  # the other sits below it. (5.38.0 on GOAatf: ratio 2.15, log offset 0.77,
  # SE 0.21, objective 356.56 against 364.01 without the offset.)
  testthat::expect_equal(max(s), 1)
})
