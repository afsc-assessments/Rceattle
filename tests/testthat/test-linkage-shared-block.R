# Fleets sharing a `Selectivity_index` or a `Catchability_index` estimate ONE
# parameter block, and `adjust_map_shared_params()` reconciles 15 by-fleet map
# slices onto the group's donor so the whole group uses one value. The linkage
# coefficients were never among them: `beta_linkage` is a flat
# `PARAMETER_VECTOR` with no fleet dimension, and the linkage table keys its
# rows on `fleet`, so one design column on a shared group became one FREE
# coefficient per member fleet. Estimated separately the copies diverge, and
# fleets declared to mirror each other end up with different realised
# selectivity, with no warning.
#
# The GOA cod bridge session measured the consequence on a real model: 47
# selectivity block columns became 59 linkage rows, and the duplicated copies
# diverged by up to 20.89 on the log scale, with a non-invertible Hessian and a
# fit ~33 nats below SS3's optimum on FEWER declared parameters.
#
# Reproduced here on bundled data instead. `GOA2018SS` fleets 9 and 10 are
# `ATF_bottom_trawl` and `ATF_bottom_trawl_length_comp` -- the same gear
# observed two ways -- both `Logistic`, sharing BOTH `Selectivity_index` 8 and
# `Catchability_index` 8.
#
# Measured on that pair with a real optimized fit (`estimateMode = "Hindcast"`,
# `phase = FALSE`), a `~ cut(Year, 2)` block column on `inf_asc`:
#
#   before   2 coefficients, -0.232615 and -0.352109; selectivity differs by
#            0.123 at its worst age; objective 12879.8748
#   after    1 coefficient, -0.240043, between the two; selectivity identical
#            (max abs diff 0); objective 12879.9970
#
# The objective RISES by 0.1222 nats, which is the only direction it can go:
# one fewer free parameter. That is the price of the invariant, and it is the
# same shape as the GOA cod case at smaller scale.
#
# THE TIE IS NOT ENOUGH ON ITS OWN. One coefficient makes the group's offsets
# equal only where every member HAS one. The base parameter is shared, but the
# offset the coefficient scales accumulates into a per-fleet tensor --
# `inf_offset(param, flt, sex, yr)` in `linkage.hpp` -- so a row naming one
# member moves that member alone. Fitted on the same pair:
#
#   spec names          coefficients   max |sel(9) - sel(10)|
#   both fleets                  1     0          <- the invariant
#   the lead only (fleet 9)      1     0.265
#   the follower only (fleet 10) 1     0.248
#
# So "put it on the lead fleet", which is right for a prior or the apical
# offset because those write the shared base, is WRONG for a design column.
# Partial coverage is refused instead, and `.stop_if_mirrored_block_linkage()`
# carries the three refusals that make the tie safe: partial coverage, a
# random-effect linkage (which cannot be tied at all -- the members land in
# separate sigma groups, so one map level would score the same deviations twice
# and fit two SDs to identical data; measured at 84 slots and 2 sigmas on this
# pair, where a mirrored group owes 42 and 1), and two specs that disagree on
# `link`, `bounds` or `init`, where the surviving value would otherwise be
# whichever row sits first in the table rather than the donor's.
#
# `est_phase` is deliberately NOT a conflict: a held donor holding the whole
# group is the documented rule and what an SS3 bridge needs where the reference
# model fixes some blocks. Block 4 pins that.
#
# Two things this fixture cannot show, stated rather than implied:
#
#   * `estimateMode = "DebugBuild"` cannot witness the divergence. Nothing is
#     optimized, so every coefficient sits at its start and the two fleets agree
#     trivially whether or not they share a parameter. An "identical output"
#     assertion there passes on BROKEN code -- it did, before this was caught --
#     so the output invariant is asserted only after a real fit. The map-level
#     count and the refusals are safe at DebugBuild and are kept there because
#     it is fast.
#   * The intercept column is NOT a witness either. `build_map_linkages()` pins
#     intercept rows to NA, so an intercept-only spec yields two rows that
#     cannot diverge because neither is estimated. Hence the time-block column.
#   * There is no bundled fixture for the CATCHABILITY half's fitted output.
#     The only mirrored pair available has `Catchability = "Fixed"`, so a q
#     linkage on it is correctly refused before it can reach the map. The fix
#     and all three refusals cover `process %in% c("sel", "q")` by the same
#     code path; block 6 pins the q refusal, and only the selectivity half's
#     fitted equality is measured.
#
# `by` defaults to `~ fleet` for both `sel` and `q` (`.default_stratum()`), so
# this was the DEFAULT path rather than a rare spelling. A spec with no
# `fleet =` does not avoid it -- `fleet =` is a filter and `NULL` means every
# fleet -- the shared-block spelling is an explicit `by = NULL`.

.sbl_grp <- c(9L, 10L)

.sbl_fit <- function(spec_fleet = .sbl_grp) {
  data("GOA2018SS", package = "Rceattle", envir = environment())
  suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = GOA2018SS, inits = NULL, estimateMode = "DebugBuild",
    msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(verbose = 0),
    selFun = Rceattle::build_selectivity(linkages = list(
      inf_asc = Rceattle::linkage_spec(~ cut(Year, 2), fleet = spec_fleet))))))
}

# The block design column, as model.matrix() names it on this dataset.
.sbl_blk <- "cut(Year, 2)(2e+03,2.02e+03]"

.sbl_build <- function(specs) {
  data("GOA2018SS", package = "Rceattle", envir = environment())
  suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = GOA2018SS, inits = NULL, estimateMode = "DebugBuild",
    msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(verbose = 0),
    selFun = Rceattle::build_selectivity(linkages = list(inf_asc = specs)))))
}


testthat::test_that("the fixture really is a mirrored pair of live fleets", {
  testthat::skip_on_cran()
  data("GOA2018SS", package = "Rceattle", envir = environment())
  fc <- suppressMessages(Rceattle::switch_check(GOA2018SS))$fleet_control
  # A group containing an `Off` fleet cannot exhibit the defect: there is only
  # one estimated member, so there are no two coefficients to diverge. That is
  # why GOApollock's shared group (fleet 1 + the Off fleet 7) is not used here,
  # and why the guards skip a group with fewer than two live members.
  testthat::expect_equal(unname(fc$Selectivity_index[.sbl_grp]), c(8, 8))
  testthat::expect_equal(unname(fc$Catchability_index[.sbl_grp]), c(8, 8))
  testthat::expect_true(all(as.character(fc$Fleet_type[.sbl_grp]) != "Off"))
  testthat::expect_equal(unname(as.character(fc$Selectivity[.sbl_grp])),
                         c("Logistic", "Logistic"))
})


testthat::test_that("a shared block gets ONE linkage coefficient", {
  testthat::skip_on_cran()
  m   <- .sbl_fit()
  tbl <- as.data.frame(unclass(m$data_list$linkage_table))
  mp  <- as.integer(m$obj$env$map$beta_linkage)
  gr  <- which(tbl$fleet %in% .sbl_grp)

  # Two rows per design column, one per member fleet: the table still keys on
  # fleet, which is what makes the map the place to reconcile.
  testthat::expect_equal(length(gr), 4L)
  est <- gr[!is.na(mp[gr])]
  # The intercept rows are pinned, so only the block column is estimated.
  testthat::expect_equal(length(est), 2L)
  # The defect was 2 distinct levels here.
  testthat::expect_equal(length(unique(mp[est])), 1L)
  testthat::expect_equal(sum(names(m$obj$par) == "beta_linkage"), 1L)
})


testthat::test_that("the mirrored pair's fitted selectivity is identical", {
  # The invariant a user cares about: one parameter is the mechanism, identical
  # output is the point. This needs a REAL fit -- at DebugBuild the coefficients
  # sit at their starts and the two fleets agree trivially, which is how an
  # earlier version of this block passed on broken code.
  testthat::skip_on_cran()
  data("GOA2018SS", package = "Rceattle", envir = environment())
  m <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = GOA2018SS, inits = NULL, estimateMode = "Hindcast",
    msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0,
                                        getsd = FALSE),
    selFun = Rceattle::build_selectivity(linkages = list(
      inf_asc = Rceattle::linkage_spec(~ cut(Year, 2), fleet = .sbl_grp))))))

  # One estimated coefficient, not two. Without the reconcile these were
  # -0.232615 and -0.352109.
  b <- m$opt$par[names(m$opt$par) == "beta_linkage"]
  testthat::expect_length(b, 1L)
  # The value, so a change in which row leads is visible rather than silent.
  testthat::expect_equal(unname(b), -0.240043, tolerance = 1e-4)
  testthat::expect_equal(m$opt$objective, 12879.9970, tolerance = 1e-6)

  s <- m$quantities$sel_at_age
  # Shape before value: an out-of-range slice is empty, and all.equal() on two
  # empty arrays is TRUE.
  testthat::expect_gte(dim(s)[1], max(.sbl_grp))
  a1 <- s[.sbl_grp[1], 1, , ]
  a2 <- s[.sbl_grp[2], 1, , ]
  testthat::expect_gt(length(a1), 0L)
  # Was 0.123 at its worst age.
  testthat::expect_equal(a1, a2, tolerance = 0)
})


testthat::test_that("a HELD donor holds the whole group", {
  # Reported by the GOA cod bridge session. Tying a follower only when the
  # donor's map entry was non-NA left the group untied wherever the donor row is
  # fixed -- `est_phase = 0`, or a pinned intercept -- so the follower kept a
  # FREE coefficient against a held donor. That is the same divergence on a
  # subset of rows, and it contradicts the rule `.shared_block_lead()` states:
  # a value set on the donor is what the whole group uses.
  #
  # On the SS3-bridged GOA Pacific cod model, where SS3 fixes two of `Srv`'s
  # block replacements at phase -5, 2 of `Srv_ae1`'s 18 selectivity linkage rows
  # stayed divergent. Two specs for one parameter at different phases reach it
  # on bundled data -- and that configuration must keep BUILDING, which is why
  # `est_phase` is not one of the fields the conflict check refuses.
  testthat::skip_on_cran()
  m <- .sbl_build(list(
    Rceattle::linkage_spec(~ cut(Year, 2), fleet = .sbl_grp[1], est_phase = 0),
    Rceattle::linkage_spec(~ cut(Year, 2), fleet = .sbl_grp[2], est_phase = 1)))
  tbl <- as.data.frame(unclass(m$data_list$linkage_table))
  mp  <- as.integer(m$obj$env$map$beta_linkage)
  gr  <- which(tbl$fleet %in% .sbl_grp)
  # The donor's rows really are held, so the fixture exercises the branch.
  testthat::expect_true(any(tbl$est_phase[gr] == 0))
  # Every row on the group is held: the follower kept a free level before.
  testthat::expect_true(all(is.na(mp[gr])))
  testthat::expect_equal(sum(names(m$obj$par) == "beta_linkage"), 0L)
  s <- m$quantities$sel_at_age
  testthat::expect_equal(s[.sbl_grp[1], 1, , ], s[.sbl_grp[2], 1, , ],
                         tolerance = 0)
})


testthat::test_that("partial coverage of a mirrored block is refused", {
  # The tie makes the group's offsets EQUAL; it does not give a member an offset
  # it has no row for. Naming one fleet of the pair fitted to one coefficient
  # and selectivity differing by 0.265 (lead named) / 0.248 (follower named) --
  # the mirroring broken in the opposite direction from the original defect.
  #
  # Asserted as a refusal rather than as a number because the number is what the
  # refusal exists to prevent. Both spellings must be named in the message, so a
  # user who wrote one of them is told which fleets are missing.
  testthat::skip_on_cran()
  data("GOA2018SS", package = "Rceattle", envir = environment())
  fc <- suppressMessages(Rceattle::switch_check(GOA2018SS))$fleet_control
  for (i in 1:2) {
    testthat::expect_error(.sbl_fit(spec_fleet = .sbl_grp[i]),
                           "share Selectivity_index 8", info = i)
    # The fleet NOT named is the one the message must report as missing.
    testthat::expect_error(.sbl_fit(spec_fleet = .sbl_grp[i]),
                           fc$Fleet_name[.sbl_grp[3L - i]], fixed = TRUE,
                           info = i)
  }
  # Naming the whole group is the fix the message recommends, and it builds.
  testthat::expect_equal(
    sum(names(.sbl_fit(spec_fleet = .sbl_grp)$obj$par) == "beta_linkage"), 1L)
  # Two separate specs covering the group between them also build: coverage is
  # per design column across the table, not per spec.
  testthat::expect_equal(sum(names(.sbl_build(list(
    Rceattle::linkage_spec(~ cut(Year, 2), fleet = .sbl_grp[1]),
    Rceattle::linkage_spec(~ cut(Year, 2), fleet = .sbl_grp[2])
  ))$obj$par) == "beta_linkage"), 1L)
})


testthat::test_that("a random-effect linkage on a mirrored block is refused", {
  # An RE row's `beta_linkage` is pinned at 0 and its deviation lives in
  # `beta_linkage_re`, indexed by `re_index` -- which `encode_linkage_for_tmb()`
  # asserts is a bijection over ROWS ("every slot is referenced by exactly one
  # row, so each has one data path and one density term"). So two mirrored
  # fleets get two independent deviation series: measured at 84 RE rows in two
  # sigma groups of 42, with 2 free `log_sigma_linkage` levels, where the lead
  # alone gives 42 rows in one group with 1.
  #
  # Not tied the way the fixed coefficients are, because the members land in
  # SEPARATE sigma groups: one map level would leave each group's density
  # scoring the same 42 deviations, counting the prior twice and fitting two SDs
  # to identical data. Collapsing the block into one RE group is an encoder
  # change, so this is refused until that exists.
  testthat::skip_on_cran()
  for (f in list(~ (1 | Year), ~ rw(1 | Year), ~ ar1(1 | Year))) {
    testthat::expect_error(
      .sbl_build(Rceattle::linkage_spec(f, fleet = .sbl_grp)),
      "random-effect selectivity linkage", info = deparse1(f))
    # Refused for the LEAD alone too: one named fleet still gets a series the
    # rest of the group does not.
    testthat::expect_error(
      .sbl_build(Rceattle::linkage_spec(f, fleet = .sbl_grp[1])),
      "random-effect selectivity linkage", info = deparse1(f))
  }
  # The same guard covers the catchability half, which has no fitted fixture:
  # fleet 9/10 share Catchability_index 8, and a q linkage there is refused
  # before the `Catchability = "Fixed"` check can speak.
  data("GOA2018SS", package = "Rceattle", envir = environment())
  testthat::expect_error(
    Rceattle:::.stop_if_mirrored_block_linkage(
      Rceattle:::linkage_row(process = "q", param = "q", X_col = 2L,
                             fleet = .sbl_grp[2], re_struct = "us",
                             re_group = "Year"),
      suppressMessages(Rceattle::switch_check(GOA2018SS))$fleet_control, "q"),
    "random-effect catchability linkage")
})


testthat::test_that("specs disagreeing on the shared coefficient are refused", {
  # The group's rows collapse to ONE coefficient, so a field that gives that
  # coefficient its meaning or its constraint cannot differ between them. Left
  # to merge: the surviving bound is whichever row sits first in the table (the
  # `match()` first-occurrence rule in fit_mod's bounds reduction), not the
  # donor's -- so writing the follower's spec first silently wins -- and two
  # `init` values start at their MEAN, which is neither.
  #
  # Measured before the refusal: bounds (-5, 5) and (-1, 1) merged to one
  # parameter keeping one bound; init 2.0 and -2.0 started at 0.
  testthat::skip_on_cran()
  LS  <- Rceattle::linkage_spec
  one <- function(x) { z <- list(x); names(z) <- .sbl_blk; z }
  # One spec per fleet, differing only in the field under test.
  spec <- function(i, extra) {
    do.call(LS, c(list(~ cut(Year, 2), fleet = .sbl_grp[i]), extra))
  }
  pair <- function(a, b) .sbl_build(list(spec(1L, a), spec(2L, b)))

  testthat::expect_error(
    pair(list(bounds = one(c(-5, 5))), list(bounds = one(c(-1, 1)))),
    "disagree on `bounds`")
  # The PAIR is reported, not just the lower value: `bounds = c(-5, 5)` is one
  # argument, and naming only -5 reads as a different disagreement.
  testthat::expect_error(
    pair(list(bounds = one(c(-5, 5))), list(bounds = one(c(-1, 1)))),
    "[-5, 5] vs [-1, 1]", fixed = TRUE)
  testthat::expect_error(pair(list(init = one(2)), list(init = one(-2))),
                         "disagree on `init`")
  # A default on one side and a value on the other is the same ambiguity.
  testthat::expect_error(pair(list(), list(bounds = one(c(-1, 1)))),
                         "disagree on `bounds`")
  testthat::expect_error(pair(list(link = "log"), list(link = "identity")),
                         "disagree on `link`")

  # Agreeing specs still build, so the check refuses the conflict and not the
  # configuration.
  testthat::expect_equal(sum(names(.sbl_build(list(
    LS(~ cut(Year, 2), fleet = .sbl_grp[1], bounds = one(c(-5, 5))),
    LS(~ cut(Year, 2), fleet = .sbl_grp[2], bounds = one(c(-5, 5)))
  ))$obj$par) == "beta_linkage"), 1L)
})


testthat::test_that("a prior stays on one row of the group", {
  # The fixed-beta prior loop runs over every ROW of the linkage table with no
  # lead gate, so with the coefficients tied a prior named on both fleets would
  # be counted once per member. `build_selectivity()` already refuses that --
  # "the shared block would be penalized once per sharing fleet" -- and this
  # pins that the refusal still fires with the tie in place, and that the
  # correct spelling (prior on the lead, a plain row on the follower) builds to
  # one coefficient.
  testthat::skip_on_cran()
  LS  <- Rceattle::linkage_spec
  pri <- list(Rceattle::prior_normal(0, 1)); names(pri) <- .sbl_blk
  testthat::expect_error(
    .sbl_build(LS(~ cut(Year, 2), fleet = .sbl_grp, priors = pri)),
    "penalized once per sharing fleet")
  testthat::expect_equal(sum(names(.sbl_build(list(
    LS(~ cut(Year, 2), fleet = .sbl_grp[1], priors = pri),
    LS(~ cut(Year, 2), fleet = .sbl_grp[2])
  ))$obj$par) == "beta_linkage"), 1L)
})
