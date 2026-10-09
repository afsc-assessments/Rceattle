# =============================================================================
# An intercept `init` sets the BASE parameter's starting value, for every
# process the linkage grammar covers.
#
# An intercept-bearing formula emits an "(Intercept)" row whose coefficient is
# mapped out at 0; the base parameter carries the level. So `init` on that
# column is a request about the base parameter, and it means the same thing
# whichever process the linkage is attached to. Values are given on the
# parameter's natural scale, and land logged wherever the base parameter is
# stored logged -- everything except sel_inf, whose inflections are ages or
# length midpoints held unlogged.
#
# These are pinned per process because each writes into a different tensor, and
# a process omitted from the push accepts the init, records it in the linkage
# table, and silently does nothing with it -- which is what q, selectivity and
# composition all did before.
# =============================================================================

.ibp_data <- function() {
  d <- Rceattle::BS2017SS
  yrs <- d$styr:d$projyr
  d$env_data <- data.frame(Year = yrs,
                           temp = as.numeric(scale(sin(seq_along(yrs) / 4))))
  d
}

.ibp_fit <- function(d, ...) {
  suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, estimateMode = 3, msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(getsd = FALSE, phase = FALSE,
                                        verbose = 0), ...)))
}

.ibp_spec <- function(init, ...) {
  Rceattle::linkage_spec(~ temp, init = list(intercept = init), ...)
}

# Intercept only, for the two shared-block blocks below. `~ temp` also emits a
# `temp` SLOPE column, and a slope reaches one fleet's offset tensor rather than
# the group's shared base -- so naming one member of a mirrored block gives that
# member a time-varying offset the others never get, which is refused since
# 5.57.0 (`test-linkage-shared-block.R` measures the 0.265 of divergence it
# caused). These two blocks are about where an intercept `init` lands, so they
# drop the covariate rather than carry a configuration they do not test.
.ibp_spec_int <- function(init, ...) {
  Rceattle::linkage_spec(~ 1, init = list(intercept = init), ...)
}


testthat::test_that("natural mortality takes an intercept init", {
  testthat::skip_on_cran(); testthat::skip_if_not_installed("TMB")
  fit <- .ibp_fit(.ibp_data(), M1Fun = Rceattle::build_M1(
    M1_model = 1, linkages = list(
      M1 = .ibp_spec(0.25, by = ~ species, species = 1L))))
  testthat::expect_equal(as.numeric(fit$estimated_params$log_M1[1, 1, 1]),
                         log(0.25))
})


testthat::test_that("catchability takes an intercept init", {
  testthat::skip_on_cran(); testthat::skip_if_not_installed("TMB")
  # Fleet 7 is the only Estimated-q survey in BS2017SS.
  fit <- .ibp_fit(.ibp_data(), qFun = Rceattle::build_catchability(
    linkages = list(q = .ibp_spec(0.25, by = ~ fleet, fleet = 7L))))
  testthat::expect_equal(as.numeric(fit$estimated_params$index_log_q[7]),
                         log(0.25))
})


testthat::test_that("selectivity takes an intercept init on both scales", {
  testthat::skip_on_cran(); testthat::skip_if_not_installed("TMB")
  d <- .ibp_data()
  fc <- Rceattle::switch_check(Rceattle::clean_data(d))$fleet_control
  flt <- fc$Fleet_code[which(fc$Selectivity == "Logistic")[1]]
  testthat::skip_if(is.na(flt))

  # Slopes live on log_sel_slp and are stored logged ...
  fit <- .ibp_fit(d, selFun = Rceattle::build_selectivity(linkages = list(
    slp_asc = .ibp_spec(0.25, by = ~ fleet, fleet = flt))))
  testthat::expect_equal(
    as.numeric(fit$estimated_params$log_sel_slp[1, flt, 1]), log(0.25))

  # ... inflections live on sel_inf and are not.
  fit <- .ibp_fit(d, selFun = Rceattle::build_selectivity(linkages = list(
    inf_asc = .ibp_spec(4, by = ~ fleet, fleet = flt))))
  testthat::expect_equal(
    as.numeric(fit$estimated_params$sel_inf[1, flt, 1]), 4)
})


testthat::test_that("a composition weight takes an intercept init", {
  testthat::skip_on_cran(); testthat::skip_if_not_installed("TMB")
  d <- .ibp_data()
  comp_flts <- sort(unique(d$comp_data$Fleet_code[d$comp_data$Fleet_code > 0]))
  testthat::skip_if(length(comp_flts) == 0)
  d$fleet_control$Comp_distribution <- "MultinomialAFSC"
  d$fleet_control$Comp_distribution[comp_flts] <- "DirichletMultinomial"

  fit <- .ibp_fit(d, compFun = Rceattle::build_composition(linkages = list(
    theta_comp = Rceattle::linkage_spec(~ 1, by = ~ fleet, fleet = comp_flts[1],
                                        init = list(intercept = 0.25)))))
  # The DM weight is exp(parameter), so a natural-scale 0.25 lands as log(0.25).
  testthat::expect_equal(
    as.numeric(fit$estimated_params$comp_weights[comp_flts[1]]), log(0.25))
  # Fleets the linkage did not name keep the column value.
  others <- setdiff(comp_flts, comp_flts[1])
  if (length(others)) {
    testthat::expect_equal(
      unname(fit$estimated_params$comp_weights[others]),
      unname(d$fleet_control$Comp_weights[others]))
  }
})


testthat::test_that("growth takes an intercept init", {
  testthat::skip_on_cran(); testthat::skip_if_not_installed("TMB")
  # Growth estimation needs conditional-age-at-length data, which the shared
  # fixture supplies.
  e <- new.env(parent = asNamespace("Rceattle"))
  for (f in list.files(testthat::test_path(), "^helper", full.names = TRUE)) {
    sys.source(f, e)
  }
  d <- e$make_test_data(growth = "vonBertalanffy")
  d$env_data <- data.frame(Year = d$styr:d$projyr, temp = 0)

  fit <- .ibp_fit(d, growthFun = Rceattle::build_growth(
    fun = "vonBertalanffy",
    linkages = list(K = .ibp_spec(0.25, by = ~ species, species = 1L))))
  testthat::expect_equal(
    as.numeric(fit$estimated_params$log_growth_pars[1, 1, 1]), log(0.25))
})


testthat::test_that("intercept bounds reach the base parameter too", {
  testthat::skip_on_cran(); testthat::skip_if_not_installed("TMB")
  d <- .ibp_data()
  # The range has to contain the fleet's default catchability, or build_bounds()
  # rejects the starting values -- which is itself evidence the push landed.
  fit <- .ibp_fit(d, qFun = Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ temp, by = ~ fleet, fleet = 7L,
                               bounds = list(intercept = c(1e-6, 1))))))
  testthat::expect_equal(as.numeric(fit$bounds$lower$index_log_q[7]), log(1e-6))
  testthat::expect_equal(as.numeric(fit$bounds$upper$index_log_q[7]), log(1))
})


# ---------------------------------------------------------------------------
# Refusals. Each of these would otherwise write a value on the wrong scale, or
# onto a parameter shared with another fleet -- silently, and into quantities
# that set catch advice.
# ---------------------------------------------------------------------------

testthat::test_that("a sel_inf slot that is not natural-scale is refused", {
  # Slot 2 is an inflection only for the logistic family. DoubleNormal stores
  # logit(right_floor) there, LogisticPM a log age-1 selectivity: a natural
  # 0.2 would silently become plogis(0.2) = 0.55.
  nat <- Rceattle:::.sel_inf_is_natural
  testthat::expect_true(nat(1L, "DoubleNormal"))
  testthat::expect_true(nat(2L, "Logistic"))
  testthat::expect_true(nat(2L, "DescendingLogistic"))
  testthat::expect_false(nat(2L, "DoubleNormal"))
  testthat::expect_false(nat(2L, "LogisticPM"))

  testthat::skip_on_cran(); testthat::skip_if_not_installed("TMB")
  d <- .ibp_data()
  fc <- Rceattle::switch_check(Rceattle::clean_data(d))$fleet_control
  flt <- fc$Fleet_code[which(fc$Selectivity == "Logistic")[1]]
  testthat::skip_if(is.na(flt))
  d$fleet_control$Selectivity[flt] <- "DoubleNormal"
  testthat::expect_error(
    .ibp_fit(d, selFun = Rceattle::build_selectivity(linkages = list(
      right_floor = .ibp_spec(0.2, by = ~ fleet, fleet = flt)))),
    "natural scale")
})


# A shared block, built here rather than borrowed from a bundled model, so the
# donor and the follower are known: both survey fleets are Logistic, and giving
# fleet 5 fleet 4's Selectivity_index makes them estimate one block led by 4.
.ibp_shared <- function() {
  d <- .ibp_data()
  d$fleet_control$Selectivity_index[5] <- d$fleet_control$Selectivity_index[4]
  d
}


testthat::test_that("an init on the follower of a shared block is refused", {
  testthat::skip_on_cran(); testthat::skip_if_not_installed("TMB")
  # Fleet 5's map slice is overwritten with fleet 4's, so a value set on 5 is
  # dropped without a word -- and it is the catch advice that moves.
  testthat::expect_error(
    .ibp_fit(.ibp_shared(), selFun = Rceattle::build_selectivity(linkages = list(
      slp_asc = .ibp_spec_int(0.5, by = ~ fleet, fleet = 5L)))),
    "mirrors fleet 4")
})


testthat::test_that("an init on the donor of a shared block lands", {
  testthat::skip_on_cran(); testthat::skip_if_not_installed("TMB")
  # Setting the donor is the correct way to set a shared block, so it must not
  # be refused -- the guard has to tell donor from follower, not merely detect
  # that a fleet belongs to a group.
  fit <- .ibp_fit(.ibp_shared(), selFun = Rceattle::build_selectivity(
    linkages = list(slp_asc = .ibp_spec_int(0.5, by = ~ fleet, fleet = 4L))))
  testthat::expect_equal(
    as.numeric(fit$estimated_params$log_sel_slp[1, 4, 1]), log(0.5))

  # It lands in the donor's ARRAY CELL, which is what this block is named for.
  # It does NOT land on the parameter the optimizer starts from, and that is
  # worth pinning so the name is not read as more than it says: the push writes
  # the donor's cell only, fleet 5 keeps the build default 0.5, and the shared
  # map level starts at `TMB:::updateMap()`'s mean of the two -- -0.096574,
  # where log(0.5) is -0.693147. The documented shared-block trap, with the
  # intercept push injecting per FLEET rather than per block
  # (`inst/dev/CLEANUP_BACKLOG.md`).
  testthat::expect_equal(
    as.numeric(fit$estimated_params$log_sel_slp[1, 5, 1]), 0.5)
  pv <- fit$obj$env$parList(fit$obj$par)$log_sel_slp
  testthat::expect_equal(as.numeric(pv[1, 4, 1]), mean(c(log(0.5), 0.5)),
                         tolerance = 1e-6)
  testthat::expect_equal(as.numeric(pv[1, 4, 1]), as.numeric(pv[1, 5, 1]))
})


testthat::test_that("a group of one is not treated as shared", {
  # Catchability_index counts survey indices, so it rarely equals the fleet
  # code: BS2017SS fleet 7 carries index 4 and shares with nobody. Reading the
  # key as a fleet code refused every survey but the first.
  d <- Rceattle::switch_check(Rceattle::clean_data(.ibp_data()))
  testthat::expect_true(is.na(
    Rceattle:::.shared_block_lead(d, 7L, "q")))
  testthat::expect_true(is.na(
    Rceattle:::.shared_block_lead(d, 4L, "sel")))

  # And the bundled multi-fleet example, whose two shared selectivity groups
  # are what this has to get right.
  data(GOA2018SS, package = "Rceattle", envir = environment())
  g <- Rceattle::switch_check(Rceattle::clean_data(GOA2018SS))
  fc <- g$fleet_control
  for (key in unique(fc$Selectivity_index[!is.na(fc$Selectivity_index)])) {
    rows <- which(fc$Selectivity_index == key)
    leads <- vapply(rows, function(f) Rceattle:::.shared_block_lead(g, f, "sel"),
                    integer(1))
    # Exactly one member of any group is the donor, and it is not a follower of
    # anything; a singleton group has no follower at all.
    testthat::expect_equal(sum(is.na(leads)), 1L)
    if (length(rows) > 1) {
      testthat::expect_true(all(stats::na.omit(leads) == rows[is.na(leads)]))
    }
  }

  # An unstratified linkage names no fleet, so it covers the whole group at
  # once and has nothing to disambiguate.
  testthat::expect_silent(Rceattle:::.stop_if_shared_block(g, NA, "sel", "slp_asc"))
})


testthat::test_that("an init of zero on a logged parameter is refused", {
  testthat::skip_on_cran(); testthat::skip_if_not_installed("TMB")
  d <- .ibp_data()
  fc <- Rceattle::switch_check(Rceattle::clean_data(d))$fleet_control
  flt <- fc$Fleet_code[which(fc$Selectivity == "Logistic")[1]]
  # log(0) is -Inf, which reaches TMB as a NaN objective with nothing naming
  # the cause. The bounds half of the contract already rejected this.
  testthat::expect_error(
    .ibp_fit(d, selFun = Rceattle::build_selectivity(linkages = list(
      slp_asc = .ibp_spec(0, by = ~ fleet, fleet = flt)))),
    "not > 0")
})


# A pinned (Intercept) coefficient is 0 by construction and mapped out, so a
# warm start must not leave a fitted value on it. `inits` from a slope-only
# spec has the SAME beta_linkage length as an intercept-only one, so fit_mod()'s
# length guard passes and nothing else would notice: measured at M1 6.69x and
# held there for the whole fit, with no message.
testthat::test_that("a warm start cannot leave a fitted value on a pinned intercept", {
  testthat::skip_on_cran(); testthat::skip_if_not_installed("TMB")
  d <- .ibp_data()

  # Fit 1: slope-only, so there is no (Intercept) row and the coefficient the
  # optimiser moves is the temp slope.
  slope <- .ibp_fit(d, M1Fun = Rceattle::build_M1(linkages = list(
    M1 = Rceattle::linkage_spec(~ 0 + temp))))
  stale <- slope$estimated_params
  stale$beta_linkage[] <- 1.9         # fitted slopes, standing in for real ones

  # Fit 2: intercept-only on the same parameter. Both specs are `by = ~ species`
  # by default, so beta_linkage has the SAME length either way and the length
  # guard cannot catch the mismatch.
  int_spec <- Rceattle::build_M1(linkages = list(
    M1 = Rceattle::linkage_spec(~ 1)))
  cold <- .ibp_fit(d, M1Fun = int_spec)
  warm <- .ibp_fit(d, M1Fun = int_spec, inits = stale)

  n_rows <- nrow(warm$data_list$linkage_table)
  testthat::expect_equal(n_rows, length(stale$beta_linkage))
  testthat::expect_equal(unname(warm$estimated_params$beta_linkage),
                         rep(0, n_rows), tolerance = 1e-12)
  testthat::expect_equal(warm$quantities$M1_at_age, cold$quantities$M1_at_age,
                         tolerance = 1e-10)
})
