# A retrospective peel must not shrink a random effect's process SD.
#
# Provenance: 5.15.0 stopped pinning the peeled tail of a random-effect block,
# because a pinned deviation is still scored by its density and "the deviation
# was exactly zero" is the strongest possible evidence for a small SD -- worth
# -6.6% on sigma at 5 peels, and monotone in peel depth, so it becomes a trend in
# the quantity Mohn's rho measures. Selectivity was exempted wholesale at the
# time because SOME selectivity forms are penalized with a bare sum of squares
# rather than a normalized density, and freeing one of those biases the SD the
# other way (-k*log(sigma) in the objective, unbounded above).
#
# The exemption was too broad: the exempt forms are a property of the FLEET, not
# of the parameter block, and 2DAR1/3DAR1/Hake fleets are scored by proper
# densities. These tests pin the per-fleet split so neither bias returns.
#
# The mapping this asserts lives in ceattle.cpp slots 4-5; .rce_sel_unnormalized_rows()
# carries the term-by-term table.

testthat::test_that("normalized selectivity forms are freed, ADMB-lineage ones pinned", {
  fc <- function(forms) {
    data.frame(Fleet_code = seq_along(forms),
               Selectivity = forms,
               Selectivity_index = seq_along(forms),
               stringsAsFactors = FALSE)
  }

  # Proper densities: dnorm (Hake), SCALE(SEPARABLE(AR1, AR1)) (2DAR1), GMRF (3DAR1).
  out <- Rceattle:::.rce_sel_unnormalized_rows(fc(c("Hake", "2DAR1", "3DAR1")))
  testthat::expect_length(out$sel_coff_dev, 0L)
  testthat::expect_length(out$limb, 0L)

  # Bare sums of squares: the ADMB/AMAK lineage.
  out <- Rceattle:::.rce_sel_unnormalized_rows(
    fc(c("NonParametric", "NonParametricPM", "LogisticPM")))
  testthat::expect_equal(out$sel_coff_dev, 1:3)
  # Only LogisticPM penalizes its limb deviates without a normalizing constant.
  testthat::expect_equal(out$limb, 3L)

  # Mixed model: the 3DAR1 fleet must be freed even though a NonParametricPM
  # fleet sits beside it. This is the case the block-wide exemption got wrong.
  out <- Rceattle:::.rce_sel_unnormalized_rows(
    fc(c("3DAR1", "NonParametricPM", "2DAR1")))
  testthat::expect_equal(out$sel_coff_dev, 2L)
  testthat::expect_length(out$limb, 0L)
})

testthat::test_that("integer switch codes read the same as names", {
  sel_map <- Rceattle:::.rce_allowed_map("Selectivity")
  forms <- c("3DAR1", "NonParametricPM", "Hake")
  by_name <- Rceattle:::.rce_sel_unnormalized_rows(
    data.frame(Selectivity = forms, Selectivity_index = seq_along(forms),
               stringsAsFactors = FALSE))
  by_code <- Rceattle:::.rce_sel_unnormalized_rows(
    data.frame(Selectivity = unname(sel_map[forms]),
               Selectivity_index = seq_along(forms)))
  testthat::expect_equal(by_name, by_code)
})

testthat::test_that("fleets sharing a Selectivity_index are pinned or freed together", {
  # Fleets sharing an index share ONE parameter block through shared map levels.
  # Pinning one member and freeing another would leave the pinned fleet's cells
  # at their starting value while the shared level moved, so two fleets that must
  # have identical selectivity would stop sharing it in the peeled years.
  out <- Rceattle:::.rce_sel_unnormalized_rows(
    data.frame(Selectivity = c("3DAR1", "NonParametricPM", "3DAR1"),
               Selectivity_index = c(1, 2, 1),
               stringsAsFactors = FALSE))
  testthat::expect_equal(out$sel_coff_dev, 2L)

  # ... and a group with one unnormalized member pins the whole group.
  out <- Rceattle:::.rce_sel_unnormalized_rows(
    data.frame(Selectivity = c("3DAR1", "NonParametricPM"),
               Selectivity_index = c(7, 7),
               stringsAsFactors = FALSE))
  testthat::expect_equal(out$sel_coff_dev, 1:2)
})

testthat::test_that("a retired or renamed Selectivity form errors rather than freeing a pinned tail", {
  # An empty code set reads as "nothing is unnormalized", which would free a tail
  # that must stay pinned -- silently, and only under the forms that carry the
  # ADMB penalties. Fail loudly instead.
  local_mocked <- function(expr) {
    testthat::with_mocked_bindings(
      .rce_allowed_map = function(col) c(Fixed = 0L, Logistic = 1L),
      code = expr, .package = "Rceattle")
  }
  local_mocked(testthat::expect_error(
    Rceattle:::.rce_sel_unnormalized_rows(
      data.frame(Selectivity = "Logistic", Selectivity_index = 1)),
    "no longer names"))
})


# The map a peel actually runs under. .rce_peel_map() is the production code
# retrospective() calls; driving it directly is the only way to assert what a
# peel pins without first finding a model that converges, and the selectivity
# forms this distinguishes are exactly the ones hardest to fit.

# n_flt fleets, 2 sexes, 3 bins, nyrs hindcast years + nproj projection years.
# Every cell estimated, so anything NA afterwards was pinned by the function.
fake_map <- function(n_flt = 3L, nyrs = 10L, nproj = 2L, nsex = 2L, nbin = 3L) {
  seqarr <- function(...) array(seq_len(prod(...)), dim = c(...))
  ml <- list(
    rec_dev         = seqarr(2L, nyrs + nproj),
    log_M1_dev      = seqarr(2L, nsex, 4L, nyrs + nproj),
    index_q_dev     = seqarr(n_flt, nyrs),
    log_sel_slp_dev = seqarr(2L, n_flt, nsex, nyrs),
    sel_inf_dev     = seqarr(2L, n_flt, nsex, nyrs),
    sel_coff_dev    = seqarr(n_flt, nsex, nbin, nyrs)
  )
  list(mapList = ml, mapFactor = lapply(ml, factor))
}

fc_of <- function(forms) {
  data.frame(Fleet_code = seq_along(forms), Selectivity = forms,
             Selectivity_index = seq_along(forms), stringsAsFactors = FALSE)
}

# Free cells in the peeled years, per fleet.
peeled_free <- function(map, nm, fleet_dim, peel_yrs) {
  a <- map$mapList[[nm]]
  vapply(seq_len(dim(a)[fleet_dim]), function(f) {
    sl <- if (fleet_dim == 1L) a[f, , , peel_yrs, drop = FALSE]
          else a[, f, , peel_yrs, drop = FALSE]
    sum(!is.na(sl))
  }, integer(1))
}

testthat::test_that("a random-effect selectivity block with a proper density is NOT pinned", {
  NYRS <- 10L; KEEP <- 6L; PEEL <- (KEEP + 1L):NYRS
  m <- Rceattle:::.rce_peel_map(
    fake_map(nyrs = NYRS),
    random_vars   = c("sel_coff_dev", "log_sel_slp_dev", "sel_inf_dev"),
    fleet_control = fc_of(c("3DAR1", "2DAR1", "Hake")),
    nyrs_peel = KEEP, nyrs = NYRS, nyrs_proj = NYRS + 2L)

  # This is the regression: every one of these fleets is scored by a normalized
  # density, so the peel must leave the tail free for the Laplace approximation
  # to integrate out. Pinning it shrinks sel_dev_log_sd with peel depth.
  testthat::expect_true(all(peeled_free(m, "sel_coff_dev", 1L, PEEL) > 0))
  testthat::expect_true(all(peeled_free(m, "log_sel_slp_dev", 2L, PEEL) > 0))
})

testthat::test_that("an ADMB-lineage selectivity block IS pinned even as a random effect", {
  NYRS <- 10L; KEEP <- 6L; PEEL <- (KEEP + 1L):NYRS
  m <- Rceattle:::.rce_peel_map(
    fake_map(nyrs = NYRS),
    random_vars   = c("sel_coff_dev", "log_sel_slp_dev", "sel_inf_dev"),
    fleet_control = fc_of(c("NonParametricPM", "NonParametric", "LogisticPM")),
    nyrs_peel = KEEP, nyrs = NYRS, nyrs_proj = NYRS + 2L)

  # Freeing a bare-SSQ tail adds -k*log(sigma) to the objective, which drives the
  # SD up without bound -- the same bias as pinning, in the other direction.
  testthat::expect_true(all(peeled_free(m, "sel_coff_dev", 1L, PEEL) == 0))
  # Only LogisticPM (fleet 3) penalizes its limb deviates without a constant.
  testthat::expect_equal(peeled_free(m, "sel_inf_dev", 2L, PEEL)[3], 0L)
  testthat::expect_true(all(peeled_free(m, "sel_inf_dev", 2L, PEEL)[1:2] > 0))
})

testthat::test_that("one model can free one fleet and pin another", {
  NYRS <- 10L; KEEP <- 6L; PEEL <- (KEEP + 1L):NYRS
  m <- Rceattle:::.rce_peel_map(
    fake_map(nyrs = NYRS),
    random_vars   = "sel_coff_dev",
    fleet_control = fc_of(c("3DAR1", "NonParametricPM", "Hake")),
    nyrs_peel = KEEP, nyrs = NYRS, nyrs_proj = NYRS + 2L)
  free <- peeled_free(m, "sel_coff_dev", 1L, PEEL)
  testthat::expect_true(free[1] > 0)   # 3DAR1 -- GMRF
  testthat::expect_equal(free[2], 0L)  # NonParametricPM -- bare SSQ
  testthat::expect_true(free[3] > 0)   # Hake -- dnorm
})

testthat::test_that("a fixed-effect block is pinned whatever its selectivity form", {
  NYRS <- 10L; KEEP <- 6L; PEEL <- (KEEP + 1L):NYRS
  # random_sel = FALSE: nothing is a random effect, so every tail pins as it
  # always did. This is the path every model without random_sel takes.
  m <- Rceattle:::.rce_peel_map(
    fake_map(nyrs = NYRS), random_vars = character(0),
    fleet_control = fc_of(c("3DAR1", "2DAR1", "Hake")),
    nyrs_peel = KEEP, nyrs = NYRS, nyrs_proj = NYRS + 2L)
  testthat::expect_true(all(peeled_free(m, "sel_coff_dev", 1L, PEEL) == 0))
  testthat::expect_true(all(peeled_free(m, "sel_inf_dev", 2L, PEEL) == 0))
  testthat::expect_true(all(is.na(m$mapList$rec_dev[, PEEL])))
})

testthat::test_that("rec_dev and log_M1_dev follow random_vars, and pin through the projection", {
  NYRS <- 10L; KEEP <- 6L; NPROJ <- 2L
  PEEL_ALL <- (KEEP + 1L):(NYRS + NPROJ)
  free <- Rceattle:::.rce_peel_map(
    fake_map(n_flt = 1L, nyrs = NYRS, nproj = NPROJ),
    random_vars   = c("rec_dev", "log_M1_dev"),
    fleet_control = fc_of("Logistic"),
    nyrs_peel = KEEP, nyrs = NYRS, nyrs_proj = NYRS + NPROJ)
  testthat::expect_true(all(!is.na(free$mapList$rec_dev[, PEEL_ALL])))

  pinned <- Rceattle:::.rce_peel_map(
    fake_map(n_flt = 1L, nyrs = NYRS, nproj = NPROJ),
    random_vars   = character(0),
    fleet_control = fc_of("Logistic"),
    nyrs_peel = KEEP, nyrs = NYRS, nyrs_proj = NYRS + NPROJ)
  # Pinned through the PROJECTION, not just the peeled hindcast years.
  testthat::expect_true(all(is.na(pinned$mapList$rec_dev[, PEEL_ALL])))
  testthat::expect_true(all(is.na(pinned$mapList$log_M1_dev[, , , PEEL_ALL])))
})

testthat::test_that("retained years are never touched", {
  NYRS <- 10L; KEEP <- 6L
  m <- Rceattle:::.rce_peel_map(
    fake_map(nyrs = NYRS), random_vars = character(0),
    fleet_control = fc_of(c("3DAR1", "NonParametricPM", "Hake")),
    nyrs_peel = KEEP, nyrs = NYRS, nyrs_proj = NYRS + 2L)
  testthat::expect_true(all(!is.na(m$mapList$sel_coff_dev[, , , seq_len(KEEP)])))
  testthat::expect_true(all(!is.na(m$mapList$rec_dev[, seq_len(KEEP)])))
})

# NULL (unrecorded) and character(0) (recorded, and there were none) pin the
# same deviations -- the map has no third option -- but they are NOT the same
# situation, and only the first is worth telling anyone about. The warning lives
# in retrospective(), raised once before the peels are dispatched, because a
# warning raised in here is discarded by the parallel worker that runs the peel;
# see test-functions-retrospective.R. This helper stays silent for both.
testthat::test_that("an unrecorded or empty random_vars pins everything, silently", {
  NYRS <- 10L; KEEP <- 6L; PEEL <- (KEEP + 1L):NYRS
  call_with <- function(rv) Rceattle:::.rce_peel_map(
    fake_map(n_flt = 1L, nyrs = NYRS), random_vars = rv,
    fleet_control = fc_of("Hake"),
    nyrs_peel = KEEP, nyrs = NYRS, nyrs_proj = NYRS + 2L)

  m_null <- testthat::expect_no_warning(call_with(NULL))
  m_none <- testthat::expect_no_warning(call_with(character(0)))

  testthat::expect_true(all(is.na(m_null$mapList$sel_coff_dev[, , , PEEL])))
  testthat::expect_true(all(is.na(m_null$mapList$rec_dev[, PEEL])))
  # The two agree block for block, so nothing downstream has to distinguish them.
  testthat::expect_equal(m_null$mapList, m_none$mapList)
  testthat::expect_equal(m_null$mapFactor, m_none$mapFactor)
})

testthat::test_that("a selectivity block whose fleet dimension disagrees is an error", {
  NYRS <- 10L
  testthat::expect_error(
    Rceattle:::.rce_peel_map(
      fake_map(n_flt = 3L, nyrs = NYRS), random_vars = "sel_coff_dev",
      fleet_control = fc_of(c("Hake", "Hake")),   # 2 fleets, map has 3
      nyrs_peel = 6L, nyrs = NYRS, nyrs_proj = NYRS + 2L),
    "cannot tell which fleet's deviations")
})
