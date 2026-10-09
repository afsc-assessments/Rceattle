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
# Three things this fixture cannot show, stated rather than implied:
#
#   * `estimateMode = "DebugBuild"` cannot witness the divergence. Nothing is
#     optimized, so every coefficient sits at its start and the two fleets agree
#     trivially whether or not they share a parameter. An "identical output"
#     assertion there passes on BROKEN code -- it did, before this was caught --
#     so the output invariant is asserted only after a real fit. The map-level
#     count is safe at DebugBuild and is kept there because it is fast.
#   * The intercept column is NOT a witness either. `build_map_linkages()` pins
#     intercept rows to NA, so an intercept-only spec yields two rows that
#     cannot diverge because neither is estimated. Hence the time-block column.
#   * There is no bundled fixture for the CATCHABILITY half. The only mirrored
#     pair available has `Catchability = "Fixed"`, so a q linkage on it is
#     correctly refused before it can reach the map. The q hole is real by
#     inspection -- same table, same flat vector, same `.shared_block_lead()` --
#     and the fix covers `process %in% c("sel", "q")`, but only the selectivity
#     half is pinned by a test.
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


testthat::test_that("the fixture really is a mirrored pair of live fleets", {
  testthat::skip_on_cran()
  data("GOA2018SS", package = "Rceattle", envir = environment())
  fc <- suppressMessages(Rceattle::switch_check(GOA2018SS))$fleet_control
  # A group containing an `Off` fleet cannot exhibit the defect: there is only
  # one estimated member, so there are no two coefficients to diverge. That is
  # why GOApollock's shared group (fleet 1 + the Off fleet 7) is not used here.
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


testthat::test_that("the group shares one coefficient however the spec names it", {
  # `.stop_if_shared_block()` is NOT a general refusal: it is called only from
  # the intercept init push (`2-build_params.R`) and the bounds push
  # (`4-build_parameter_bounds.R`), so a plain formula naming a single follower
  # is accepted. An earlier version of this block asserted it threw, which was
  # wrong about the code rather than the other way round.
  #
  # What must hold is the invariant, in all three spellings a user might write:
  # whether the spec names the donor, the follower, or both, the group ends up
  # with ONE coefficient.
  testthat::skip_on_cran()
  for (flts in list(.sbl_grp, .sbl_grp[1], .sbl_grp[2])) {
    m <- .sbl_fit(spec_fleet = flts)
    testthat::expect_equal(sum(names(m$obj$par) == "beta_linkage"), 1L,
                           info = paste(flts, collapse = ","))
  }
})
