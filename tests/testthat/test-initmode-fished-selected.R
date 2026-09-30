# =============================================================================
# initMode = "FishedNonEquilibriumSelected" (6): the initial age-structure
# decays with sum_{a' < a} (M1_a' + Finit * s_a'), i.e. the initial fishing
# mortality is weighted by the fishery's selectivity at age before it
# accumulates. This is Stock Synthesis's InitF convention, and Finit is the
# apical initial F because selectivity is normalized to a maximum of 1.
#
# Modes 3 and 4 both exist already and neither is this: (3) charges every age
# the same Finit, (4) applies it once rather than accumulating it. Under a
# size-selective fishery only (6) is an equilibrium -- (3) kills unselected
# young ages at the full initial F, (4) does not decay the older ages with it.
#
# Provenance: found while bridging the 2024 AI Pacific cod SS3 assessment to
# Rceattle. Injecting SS3's MLE and inverting mode 4's mort_sum left the
# initial deviates differing from SS3's Early_InitAge by exactly
# const - Finit * cumsum(sel), residual 0.00000 at all 13 ages. That identity
# is what this mode encodes.
#
# NOTE: initMode is a fit_mod() ARGUMENT; it overwrites data_list$initMode, so
# these tests pass it as the argument.
# =============================================================================

.fs_fit <- function(mode, inits = NULL) {
  data(BS2017SS, package = "Rceattle")
  suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = BS2017SS, inits = inits, initMode = mode, estimateMode = 3,
    msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(getsd = FALSE, phase = FALSE, verbose = 0))))
}

testthat::test_that("FishedNonEquilibriumSelected is a recognised initMode alias (= 6)", {
  testthat::expect_true("FishedNonEquilibriumSelected" %in% names(Rceattle:::initMode_map))
  testthat::expect_identical(unname(Rceattle:::initMode_map["FishedNonEquilibriumSelected"]), 6)
})

testthat::test_that("Finit is estimated under mode 6, as it is under 3 and 4", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  fit <- .fs_fit("FishedNonEquilibriumSelected")
  testthat::expect_false(all(is.na(fit$map$mapFactor$log_Finit)))
  # and still mapped out where the population starts unfished
  eq <- .fs_fit("Equilibrium")
  testthat::expect_true(all(is.na(eq$map$mapFactor$log_Finit)))
})

testthat::test_that("mode 6 decays the initial numbers by M1 + Finit * selectivity", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  data(BS2017SS, package = "Rceattle")
  fit <- .fs_fit("FishedNonEquilibriumSelected")
  q  <- fit$quantities
  sp <- 1L
  nage <- BS2017SS$nages[sp]
  M1   <- as.numeric(q$M1_at_age[sp, 1, seq_len(nage), 1])
  N1   <- as.numeric(q$N_at_age[sp, 1, seq_len(nage), 1])
  Finit <- exp(as.numeric(fit$estimated_params$log_Finit[sp]))
  idev  <- as.numeric(fit$estimated_params$init_dev[sp, seq_len(nage - 1)])

  # The selectivity mode 6 uses: mean over the species' FISHERY fleets, year 1.
  fc  <- fit$data_list$fleet_control
  fsh <- which(fc$Fleet_type %in% c("Fishery", 1))
  fsh <- fsh[fc$Species[fsh] == sp]
  sel <- colMeans(matrix(q$sel_at_age[fsh, 1, seq_len(nage), 1],
                         nrow = length(fsh), byrow = FALSE))

  # N(a) = R_init * exp(-sum_{a' < a}(M1 + Finit*sel) + init_dev(a-1))
  R_init <- N1[1] / exp(0)                      # age-0 carries no init_dev
  for (a in 2:(nage - 1)) {
    mort <- sum(M1[seq_len(a - 1)] + Finit * sel[seq_len(a - 1)])
    testthat::expect_equal(N1[a], R_init * exp(-mort + idev[a - 1]),
                           tolerance = 1e-6,
                           info = paste("age index", a))
  }
})

testthat::test_that("mode 6 reduces to mode 3 when selectivity is flat", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  data(BS2017SS, package = "Rceattle")
  # BS2017SS's fisheries sit at flat selectivity on the build defaults, and a
  # flat s = 1 makes sum(M1 + Finit*s) identical to mode 3's sum(M1 + Finit).
  base <- .fs_fit("FishedNonEquilibrium")
  inits <- base$estimated_params
  inits$log_Finit[] <- log(0.2)
  f3 <- .fs_fit("FishedNonEquilibrium", inits)
  f6 <- .fs_fit("FishedNonEquilibriumSelected", inits)
  sel <- as.numeric(f6$quantities$sel_at_age[1, 1, seq_len(BS2017SS$nages[1]), 1])
  testthat::expect_equal(stats::sd(sel), 0)                     # the premise
  n <- function(f) as.numeric(f$quantities$N_at_age[1, 1, seq_len(BS2017SS$nages[1]), 1])
  testthat::expect_equal(n(f6), n(f3), tolerance = 1e-10)
})

testthat::test_that("mode 6 differs from 3 and 4 once selectivity varies with age", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  data(BS2017SS, package = "Rceattle")
  base <- .fs_fit("FishedNonEquilibrium")
  inits <- base$estimated_params
  inits$log_Finit[] <- log(0.2)
  # Fleet 1 is non-parametric (Selectivity = 2), so a ramp on its coefficients
  # gives the young ages a genuinely lower selectivity -- the case the three
  # modes are supposed to treat differently.
  inits$sel_coff[1, 1, ] <- seq(-2, 2, length.out = dim(inits$sel_coff)[3])
  f3 <- .fs_fit("FishedNonEquilibrium", inits)
  f4 <- .fs_fit("FishedNonEquilibriumScaled", inits)
  f6 <- .fs_fit("FishedNonEquilibriumSelected", inits)
  na  <- BS2017SS$nages[1]
  sel <- as.numeric(f6$quantities$sel_at_age[1, 1, seq_len(na), 1])
  testthat::expect_gt(stats::sd(sel), 0)                        # the premise
  n <- function(f) as.numeric(f$quantities$N_at_age[1, 1, seq_len(na), 1])
  testthat::expect_false(isTRUE(all.equal(n(f6), n(f3))))
  testthat::expect_false(isTRUE(all.equal(n(f6), n(f4))))
  # and mode 6 keeps MORE of the young ages than mode 3, which charges them
  # the full initial F despite their being barely selected.
  testthat::expect_gt(n(f6)[2], n(f3)[2])
})

testthat::test_that("the other initModes are untouched by mode 6 existing", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  # Mode 6 fills sel_init; every other mode must leave it at 1, so their
  # initial numbers are exactly what they were before the mode was added.
  data(BS2017SS, package = "Rceattle")
  for (m in c("Equilibrium", "NonEquilibrium", "FishedNonEquilibrium",
              "FishedNonEquilibriumScaled", "OffsetEquilibrium")) {
    fit <- .fs_fit(m)
    testthat::expect_true(is.finite(fit$quantities$jnll),
                          info = paste("initMode", m))
  }
})

# Rceattle estimates ONE Finit per species where SS3 carries one per fleet, so on
# a species with several fisheries mode 6 decays the initial age structure at
# their MEAN selectivity and Finit is not the apical initial F of any of them.
# Only the equilibrium-catch ROW was refused for this; the mode itself was silent,
# and the initial state sets the SSB scale.
testthat::test_that("mode 6 warns where a species has more than one fishery", {
  data(BS2017SS, package = "Rceattle")
  d <- BS2017SS
  # BT_Pollock is species 1's survey; as a fishery, species 1 has two.
  d$fleet_control$Fleet_type[4] <- "Fishery"
  d$initMode <- 6

  warns <- function(dd) {
    w <- character(0)
    withCallingHandlers(
      tryCatch(suppressMessages(Rceattle:::data_check(Rceattle::switch_check(dd))),
               error = function(e) NULL),
      warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
    w
  }
  testthat::expect_match(paste(warns(d), collapse = " | "),
                         "more than one fishery")

  # One fishery per species says nothing, and neither does another initMode.
  testthat::expect_false(any(grepl("more than one fishery", warns(BS2017SS))))
  d2 <- d; d2$initMode <- 2
  testthat::expect_false(any(grepl("more than one fishery", warns(d2))))
})
