# Fleets sharing a `Selectivity_index` / `Catchability_index` estimate ONE
# parameter block. `beta_linkage` has no fleet dimension, so
# `adjust_map_shared_params()` could not reach it and each member got its own
# free coefficient: fleets declared to mirror each other fitted different
# curves, silently. `build_map_linkages()` now ties them, and
# `.stop_if_mirrored_block_linkage()` refuses the four configurations the tie
# cannot honour. `inst/dev/TRAPS.md`, "Shared parameter blocks", carries the
# mechanism and every measured number; this header keeps only what a reader
# needs to maintain the file.
#
# Fixture: `GOA2018SS` fleets 9 and 10, the same gear observed two ways, both
# `Logistic`, sharing both indices. Fitted with a `~ cut(Year, 2)` block column
# on `inf_asc`, the tie takes 2 coefficients to 1 and selectivity from 0.123
# apart to identical, at a cost of 0.1222 nats -- the only direction one fewer
# free parameter can go.
#
# Three things this fixture cannot show, stated so they are not assumed:
#
#   * `estimateMode = "DebugBuild"` cannot witness divergence -- nothing is
#     optimized, so the two fleets agree trivially whether or not they share a
#     parameter, and an "identical output" assertion there passes on BROKEN
#     code, as an earlier version of this file did. The output invariant is
#     asserted only after a real fit; the map and refusal checks stay at
#     DebugBuild because they are fast.
#   * An intercept-only spec cannot witness it either: those rows are pinned to
#     NA, so neither is estimated. Hence the time-block column.
#   * No bundled pair has an estimated `q`, so the catchability half's refusals
#     are pinned but its fitted equality is not.
#
# `by` defaults to `~ fleet`, so this was the DEFAULT path; the shared-block
# spelling is an explicit `by = NULL`.

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
  # A group containing an `Off` fleet cannot exhibit the OUTPUT divergence:
  # there is only one estimated member, so there are no two curves to differ.
  # That is why GOApollock's shared group (fleet 1 + the Off fleet 7) is not
  # used here, and why COVERAGE is owed only to a member whose block is
  # estimated. The agreement rules still run on such a group, because the map
  # ties a follower's coefficient whatever its Fleet_type -- see the Off-member
  # block below.
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
  # One parameter is the mechanism; identical output is the point. Needs a REAL
  # fit -- at DebugBuild the coefficients sit at their starts and agree
  # trivially, which is how an earlier version passed on broken code.
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
  # Tying a follower only where the donor's entry was non-NA left the group
  # untied wherever the donor row is fixed, so the follower kept a free
  # coefficient against a held donor. Two specs at different phases reach it,
  # and that configuration must keep BUILDING -- which is why `est_phase` is not
  # one of the fields the conflict check refuses.
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
  # The tie makes the group's offsets equal; it does not give a member an offset
  # it has no row for. Naming one fleet fitted to one coefficient and 0.265 of
  # divergence -- the mirroring broken in the other direction. Asserted as a
  # refusal because the number is what the refusal prevents.
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
  # An RE row's deviation lives in `beta_linkage_re`, indexed by `re_index`,
  # which the encoder asserts is a bijection over ROWS -- so each named fleet
  # gets its own series and its own SD (84 slots in two sigma groups of 42,
  # against 42 and 1 for the lead alone). Not tied, because separate sigma
  # groups would score the same deviations twice.
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
  # One coefficient, so a field giving it its meaning or constraint cannot
  # differ. Merged, the surviving bound was whichever row sat first in the table
  # -- so writing the follower's spec first silently won -- and two `init`s
  # started at their mean.
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


testthat::test_that("two rows for one design column on one fleet are refused", {
  # The tie equalises COEFFICIENTS, not row counts: the offset is added once per
  # row, so a member with two rows moves twice as far on the same coefficient --
  # one reported parameter and two curves. The mirror shape broke the tie
  # instead, because the old guard required the donor to own exactly one match.
  testthat::skip_on_cran()
  LS  <- Rceattle::linkage_spec
  dup <- "rows for design column"
  testthat::expect_error(
    .sbl_build(list(LS(~ cut(Year, 2), fleet = .sbl_grp),
                    LS(~ cut(Year, 2), fleet = .sbl_grp[2]))),
    dup)
  testthat::expect_error(
    .sbl_build(list(LS(~ cut(Year, 2), fleet = .sbl_grp),
                    LS(~ cut(Year, 2), fleet = .sbl_grp[1]))),
    dup)
  # The SS3-bridge shape: the donor named twice at different phases. This is the
  # case the earlier guard swallowed, and the one whose two free levels left the
  # mirrored pair 0.295 apart.
  testthat::expect_error(
    .sbl_build(list(
      LS(~ cut(Year, 2), fleet = .sbl_grp[1], est_phase = 0),
      LS(~ cut(Year, 2), fleet = .sbl_grp[1], est_phase = 1),
      LS(~ cut(Year, 2), fleet = .sbl_grp[2], est_phase = 1))),
    dup)
  # A single spec that merely repeats a fleet in its filter is NOT a duplicate:
  # the strata expand to one row per fleet.
  testthat::expect_equal(sum(names(.sbl_fit(
    spec_fleet = c(.sbl_grp, .sbl_grp[2]))$obj$par) == "beta_linkage"), 1L)
})


testthat::test_that("an intercept row is not held to coverage", {
  # An intercept re-targets the shared base and is pinned to NA, so it carries no
  # per-fleet offset. Holding it to coverage made an intercept `init` on a
  # mirrored block unreachable by every spelling.
  testthat::skip_on_cran()
  data("GOA2018SS", package = "Rceattle", envir = environment())
  m <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = GOA2018SS, inits = NULL, estimateMode = "DebugBuild",
    msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(verbose = 0),
    selFun = Rceattle::build_selectivity(linkages = list(
      inf_asc = Rceattle::linkage_spec(~ 1, fleet = .sbl_grp[1]))))))
  tbl <- as.data.frame(unclass(m$data_list$linkage_table))
  gr  <- which(tbl$fleet %in% .sbl_grp)
  testthat::expect_equal(length(gr), 1L)
  testthat::expect_true(all(Rceattle:::.is_pinned_intercept(tbl[gr, ])))
  # Pinned, so it frees nothing -- which is why coverage is meaningless for it.
  testthat::expect_true(all(is.na(as.integer(m$obj$env$map$beta_linkage)[gr])))
})


testthat::test_that("a member not estimating the block is not owed one", {
  # Coverage is owed only to a member whose block is estimated. Demanding it from
  # every live member deadlocked a q linkage whose follower holds q fixed:
  # naming the lead failed coverage, naming both failed the `Catchability` check.
  testthat::skip_on_cran()
  data("GOA2018SS", package = "Rceattle", envir = environment())
  d <- GOA2018SS
  lead <- 12L; nonq <- 14L
  fc <- suppressMessages(Rceattle::switch_check(d))$fleet_control
  # Fixture validity: the follower really is a fleet whose q is not estimated.
  testthat::expect_equal(as.character(fc$Catchability[lead]), "Estimated")
  testthat::expect_true(is.na(fc$Catchability[nonq]))
  d$fleet_control$Catchability_index[nonq] <-
    d$fleet_control$Catchability_index[lead]
  m <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = d, inits = NULL, estimateMode = "DebugBuild", msmMode = 0,
    random_rec = FALSE, fit_control = Rceattle::fit_control(verbose = 0),
    qFun = Rceattle::build_catchability(linkages = list(
      q = Rceattle::linkage_spec(~ cut(Year, 2), fleet = lead))))))
  testthat::expect_equal(sum(names(m$obj$par) == "beta_linkage"), 1L)
})


testthat::test_that("a group with an Off member is still reconciled", {
  # The map ties a follower whatever its `Fleet_type`, so the agreement rules
  # must run on the same set. Skipping a group with one live member let the `Off`
  # fleet's bound win by table position.
  testthat::skip_on_cran()
  data("GOA2018SS", package = "Rceattle", envir = environment())
  fc <- suppressMessages(Rceattle::switch_check(d <- GOA2018SS))$fleet_control
  off_grp <- c(7L, 1L)   # Off member first, which is how the bound was taken
  testthat::expect_equal(unname(fc$Selectivity_index[off_grp]), c(1, 1))
  testthat::expect_equal(as.character(fc$Fleet_type[off_grp[1]]), "Off")
  LS <- Rceattle::linkage_spec
  one <- function(x) { z <- list(x); names(z) <- .sbl_blk; z }
  testthat::expect_error(
    .sbl_build(list(LS(~ cut(Year, 2), fleet = off_grp[1],
                       bounds = one(c(-1, 1))),
                    LS(~ cut(Year, 2), fleet = off_grp[2],
                       bounds = one(c(-5, 5))))),
    "disagree on `bounds`")
  # Coverage, though, is NOT owed to the Off member: its curve is not estimated.
  testthat::expect_equal(sum(names(.sbl_build(
    LS(~ cut(Year, 2), fleet = off_grp[2]))$obj$par) == "beta_linkage"), 1L)
})


testthat::test_that("a q prior on two members of a block is refused", {
  # The fixed-beta prior loop runs over every ROW with no lead gate, so with the
  # coefficients tied a prior named on two members is counted twice on the one
  # they share -- dividing its stated variance by two. `build_selectivity()` has
  # refused the selectivity equivalent since 5.42.0; the q half had none, and the
  # coverage rule is what pushes a caller toward naming both.
  testthat::skip_on_cran()
  data("GOA2018SS", package = "Rceattle", envir = environment())
  d <- GOA2018SS
  grp <- c(12L, 13L)
  d$fleet_control$Catchability_index[grp[2]] <-
    d$fleet_control$Catchability_index[grp[1]]
  qbuild <- function(specs) suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = d, inits = NULL, estimateMode = "DebugBuild", msmMode = 0,
    random_rec = FALSE, fit_control = Rceattle::fit_control(verbose = 0),
    qFun = Rceattle::build_catchability(linkages = list(q = specs)))))
  LS  <- Rceattle::linkage_spec
  # Discover the block column name rather than hard-coding a second copy of it.
  qtbl <- as.data.frame(unclass(
    qbuild(LS(~ cut(Year, 2), fleet = grp))$data_list$linkage_table))
  blk <- setdiff(unique(qtbl$design_col), "(Intercept)")
  testthat::expect_length(blk, 1L)
  pri <- list(Rceattle::prior_normal(0.4, 0.1)); names(pri) <- blk
  testthat::expect_error(
    qbuild(LS(~ cut(Year, 2), fleet = grp, priors = pri)),
    "once per sharing fleet")
  # The accepted spelling: prior on one member, plain row for the other.
  testthat::expect_equal(sum(names(qbuild(list(
    LS(~ cut(Year, 2), fleet = grp[1], priors = pri),
    LS(~ cut(Year, 2), fleet = grp[2])))$obj$par) == "beta_linkage"), 1L)
})


testthat::test_that("a prior stays on one row of the group", {
  # With the coefficients tied, a prior named on both fleets would be counted
  # twice on one coefficient. `build_selectivity()` refuses that already; this
  # pins that it still fires, and that the correct spelling builds.
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
