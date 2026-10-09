# `Time_varying_sel = "RandomWalkAscending"` walks the ASCENDING limb only, and
# `build_map()` assigns deviate indices for it on `DoubleLogistic` alone -- which
# is the switch's own definition ("random walk on ascending portion of double
# logistic only", R/3-build_map.R).
#
# On every other form the combination fell through both arms of that form's
# branch, so no deviate index was assigned, every deviate stayed NA, and a
# STATIC curve fitted where the workbook asked for a walk. No condition was
# raised. Measured on the GOA pollock fishery (fleet 8), switching only
# `Selectivity` and `Time_varying_sel`:
#
#   form                Off        RandomWalk   RandomWalkAscending
#   DoubleLogistic      220 / 0    412 / 192    316 / 96   <- implemented
#   DoubleNormal        220 / 0    412 / 192    220 / 0    <- silent static fit
#   Logistic            218 / 0    314 / 96     218 / 0    <- silent static fit
#   DescendingLogistic  218 / 0    314 / 96     218 / 0    <- silent static fit
#
# (npar / free selectivity deviates. Fleets 1-7 of that dataset are all
# `Time_varying_sel = "Off"`, so every selectivity deviate in `obj$par` belongs
# to fleet 8.)
#
# It was worse under `random_sel = TRUE`: the deviation sd is freed for any form
# that is not `Fixed` whenever `Time_varying_sel` is one of IID / AR1 /
# RandomWalk / RandomWalkAscending, so `sel_dev_log_sd` became a free parameter
# scaling deviations that were all mapped out -- the same stray-hyperparameter
# shape as `M1_re` under an unimplemented `M1_model` (5.53.0).
#
# Refused rather than implemented, because the same model is already reachable:
# a random-effect linkage on the ascending parameters -- `peak` and `sigma_asc`
# for DoubleNormal, `inf_asc` and `slp_asc` for the logistic family, which alias
# the same slots the template scores -- gives a time-varying ascending limb with
# the descending one fixed. `Time_varying_sel` on a parametric form is
# soft-deprecated in favour of those linkages, so implementing mode 5 for a
# second form would add a route to a mechanism already on its way out.
#
# `2DAR1` / `3DAR1` are exempt: they estimate their field for every bin and year
# regardless of `Time_varying_sel`, and `build_map()` says so with a warning.
#
# Nothing shipping asked for it. Of the 15 bundled datasets, 11 carry a
# `fleet_control` directly and only `GOA2018SS` and `GOApollock` set
# `RandomWalkAscending`, both on the `DoubleLogistic` GOA pollock fishery; no
# script in the sibling assessment repos sets it at all.
#
# `AR1` is refused package-wide ("Time_varying_sel = 'AR1' is removed"), so the
# `c("IID", "AR1", "RandomWalk")` sets in `build_map()` name a value that can no
# longer arrive.

# The GOA pollock fishery, which ships DoubleLogistic + RandomWalkAscending.
.rwa_fleet <- 8L

.rwa_build <- function(sel, tv, ...) {
  data("GOApollock", package = "Rceattle", envir = environment())
  d <- GOApollock
  d$fleet_control$Selectivity[.rwa_fleet] <- sel
  d$fleet_control$Time_varying_sel[.rwa_fleet] <- tv
  suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = d, inits = NULL, estimateMode = "DebugBuild", msmMode = 0,
    random_rec = FALSE, fit_control = Rceattle::fit_control(verbose = 0), ...)))
}

.rwa_devs <- function(fit) {
  nm <- names(fit$obj$par)
  c(npar = length(nm), devs = sum(nm %in% c("sel_inf_dev", "log_sel_slp_dev")))
}


testthat::test_that("RandomWalkAscending is refused on the forms that drop it", {
  testthat::skip_on_cran()   # builds a real fit_mod() object
  # Named in the message so a user can act on it, and the form named too.
  for (sel in c("DoubleNormal", "Logistic", "DescendingLogistic")) {
    testthat::expect_error(.rwa_build(sel, "RandomWalkAscending"),
                           "RandomWalkAscending", info = sel)
    # Anchored: a bare "Logistic" also matches "'DoubleLogistic' selectivity
    # only" in the message, so it would pass without naming the fleet's form.
    testthat::expect_error(.rwa_build(sel, "RandomWalkAscending"),
                           paste0("On '", sel, "'"), info = sel)
  }
  # And still builds on the one form that implements it, with its deviates.
  got <- .rwa_devs(.rwa_build("DoubleLogistic", "RandomWalkAscending"))
  off <- .rwa_devs(.rwa_build("DoubleLogistic", "Off"))
  testthat::expect_equal(unname(got[["devs"]]), 96L)
  # The DIFFERENCE, not the absolute count: a parameter added anywhere else
  # would otherwise fail this file with a selectivity-sounding message.
  testthat::expect_equal(unname(got[["npar"]] - off[["npar"]]), 96L)
})


testthat::test_that("the refused forms still take RandomWalk, which they do map", {
  testthat::skip_on_cran()   # builds a real fit_mod() object
  # The refusal must not take away the structure these forms DO support, which
  # is the whole walk rather than its ascending limb.
  for (sel in c("DoubleNormal", "Logistic", "DescendingLogistic")) {
    got <- .rwa_devs(.rwa_build(sel, "RandomWalk"))
    testthat::expect_gt(got[["devs"]], 0L)
    off <- .rwa_devs(.rwa_build(sel, "Off"))
    testthat::expect_equal(unname(off[["devs"]]), 0L)
    # A walk must free strictly more than Off, which is what the dropped
    # combination failed to do.
    testthat::expect_gt(got[["npar"]], off[["npar"]])
  }
})


testthat::test_that("no selectivity form can silently drop RandomWalkAscending", {
  # Derived from the source, so a NEW form added without a RandomWalkAscending
  # arm is caught here rather than by a user. For each `sel_type` branch in
  # build_map(), does its tv_sel membership test name RandomWalkAscending? If
  # not, data_check() must refuse the combination -- either by this rule or by
  # that form's own per-form restriction.
  f <- c("R/3-build_map.R", testthat::test_path("..", "..", "R", "3-build_map.R"))
  f <- f[file.exists(f)]
  testthat::skip_if(length(f) == 0, "R/3-build_map.R not available")
  src <- readLines(f[1], warn = FALSE)

  guards <- grep('sel_type\\s*(==|%in%)', src)
  tests  <- grep('tv_sel\\s*%in%', src)
  # Pinned, not floored: an `expect_gt(., 5L)` would not notice a branch
  # rewritten as `switch()` or its mode set hoisted into a named vector, either
  # of which drops it out of the denominator and leaves this block green while
  # measuring less. Update deliberately.
  testthat::expect_equal(length(guards), 12L)
  testthat::expect_equal(length(tests), 9L)

  implements <- character()
  for (ln in tests) {
    owner <- guards[guards < ln]
    if (!length(owner)) next
    form <- regmatches(src[max(owner)],
                       gregexpr('"([A-Za-z0-9]+)"', src[max(owner)]))[[1]]
    form <- gsub('"', '', form)
    if (grepl("RandomWalkAscending", src[ln])) implements <- c(implements, form)
  }
  # Exactly one form assigns RandomWalkAscending deviate indices.
  testthat::expect_equal(sort(unique(implements)), "DoubleLogistic")
  # Every form the refusal names must be a real one, so a renamed form cannot
  # leave the rule pointing at nothing.
  testthat::expect_true(all(c("DoubleLogistic", "DoubleNormal", "Logistic",
                              "DescendingLogistic", "2DAR1", "3DAR1") %in%
                            names(sel_map)))
})
