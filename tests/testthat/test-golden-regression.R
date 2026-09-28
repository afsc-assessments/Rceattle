# =============================================================================
# Golden-reference regression: the four reference fits (Bering Sea + Gulf of
# Alaska, single- and multi-species) must reproduce their pinned objective
# functions. This is the committed backstop for "numeric changes need golden
# equivalence" -- any edit that moves a reference fit fails here.
#
# The multispecies fits are warm-started from their single-species MLEs (the
# predation likelihood is non-convex, so the start point selects the local
# optimum); the GOA multispecies fit uses fixed M. See the /golden-check
# command for the recipe these constants were generated from.
# =============================================================================

# The GOA reference objectives changed in 5.33.0 (goa_ss -0.0150, goa_ms -0.0030):
# GOA2018SS fleet 2 has a q prior, which is now mean-centred under bias_adjust_proc.
#
# WHAT THIS PINS, AND WHY IT NO LONGER RE-OPTIMIZES. The question golden exists
# to answer is "is this still the same likelihood". Re-fitting answered a second
# question as well -- "does a cold phased optimization find the lower of two
# minima on this machine" -- and that one has a platform-dependent answer:
# `goa_ss` has a second local minimum 52.9 units up, a one-ULP change in a single
# log_F gradient element is enough for nlminb to reach it, and `goa_ms` inherits
# the basin through its warm start. The two GOA references therefore moved
# together on the CI runners while reproducing exactly on local macOS, and the
# job that runs this file was red on `main` for weeks over a summation order.
#
# The likelihood is now evaluated AT the reference parameters, which removes the
# optimizer from the gate entirely: `fit_mod(estimateMode = 3)` builds the object
# and evaluates without optimizing, so any change to the model moves the
# objective while a change to the optimizer's path cannot. Measured at 5.45.0:
# all four reproduce their pinned objective to 1e-12 absolute (1e-16 relative)
# and sit at max|gradient| between 1.1e-11 and 3.3e-06, against the 1e-4
# threshold. The multispecies pair was the open question -- predation iterates
# inside fit_mod() -- and the iteration reaches the same state from the
# reference parameters, so they pin as cleanly as the single-species pair.
#
# The cold-start path is still exercised below, on convergence only. Regenerate
# the fixture with tools/verify/regenerate-golden-reference.R, and only when a
# model change is intended.
testthat::test_that("the likelihood at the reference parameters is unchanged", {
  testthat::skip_on_cran()
  # covr instruments the TMB model at -O0. The evaluation here no longer depends
  # on an optimizer path, so the historical reason for this skip is gone, but a
  # -O0 build's arithmetic is still not what the reference was taken under.
  testthat::skip_on_covr()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  ref <- readRDS(testthat::test_path("fixtures", "golden-reference.rds"))

  # Each reference's own configuration, matching the recipe the fixture script
  # fits: the multispecies pair warm-starts from its single-species MLEs and the
  # GOA multispecies fit holds M fixed.
  spec <- list(
    ss     = list(data = Rceattle::BS2017SS,  msmMode = 0, niter = 3),
    ms     = list(data = Rceattle::BS2017MS,  msmMode = 1, niter = 5),
    goa_ss = list(data = Rceattle::GOA2018SS, msmMode = 0, niter = 3),
    goa_ms = list(data = Rceattle::GOA2018SS, msmMode = 1, niter = 3))

  for (m in names(spec)) {
    s <- spec[[m]]
    fit <- Rceattle::fit_mod(
      data_list = s$data, inits = ref[[m]]$params, file = NULL,
      estimateMode = 3, niter = s$niter, random_rec = FALSE,
      msmMode = s$msmMode, suitMode = 0,
      fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0))

    par <- fit$obj$par
    testthat::expect_equal(as.numeric(fit$obj$fn(par)), ref[[m]]$objective,
                           tolerance = 1e-8, info = m)
    # And the reference is still a stationary point of it, so the pin cannot
    # drift to a number that merely reproduces at a non-stationary point.
    testthat::expect_lt(max(abs(fit$obj$gr(par))), 1e-4)
  }
})

# Convergence from a cold start, with NO objective pinned. This is the half that
# used to make the gate platform-dependent: which of `goa_ss`'s two minima is
# reached varies with summation order, and both are stationary points, so the
# gradient holds on either. A failure here means the optimizer or the phasing
# stopped working, not that the model changed.
testthat::test_that("a cold start still converges", {
  testthat::skip_on_cran()
  testthat::skip_on_covr()
  testthat::skip_if_not_installed("TMB")

  fc <- function(...) Rceattle::fit_control(getsd = FALSE, verbose = 0,
                                            newtonsteps = 3, ...)

  ss <- Rceattle::fit_mod(data_list = Rceattle::BS2017SS, file = NULL,
    inits = NULL, estimateMode = 0, random_rec = FALSE, msmMode = 0,
    fit_control = fc(phase = TRUE))
  goa_ss <- Rceattle::fit_mod(data_list = Rceattle::GOA2018SS, file = NULL,
    inits = NULL, estimateMode = 0, random_rec = FALSE, msmMode = 0,
    fit_control = fc(phase = TRUE))
  ms <- Rceattle::fit_mod(data_list = Rceattle::BS2017MS,
    inits = ss$estimated_params, file = NULL, estimateMode = 0, niter = 5,
    random_rec = FALSE, msmMode = 1, suitMode = 0, fit_control = fc())
  goa_ms <- Rceattle::fit_mod(data_list = Rceattle::GOA2018SS,
    inits = goa_ss$estimated_params, file = NULL, estimateMode = 0, niter = 3,
    random_rec = FALSE, msmMode = 1, suitMode = 0, fit_control = fc(phase = TRUE))

  for (f in list(ss = ss, ms = ms, goa_ss = goa_ss, goa_ms = goa_ms))
    testthat::expect_lt(max(abs(f$obj$gr(f$opt$par))), 1e-4)
})

# =============================================================================
# Golden reference by CONFIGURATION.
#
# The four fits above share one configuration: plain multinomial compositions, a
# lognormal index, no tail accumulation, fixed recruitment. So the likelihood
# branches the ecosystem models actually run had no golden cover, and a change
# confined to one of them could not fail here. Plain Multinomial is not even the
# schema default for Comp_distribution -- MultinomialAFSC is -- so the branch a
# model gets when the column is left out was unpinned as well. Each block below
# switches on a single branch, on the same BS2017SS data, so a failure points at
# the branch rather than at the model.
#
# Still unpinned after this, and worth adding if they ever break: CAAL, length
# compositions, nonparametric selectivity, and any multispecies combination of
# the branches below -- both multispecies references run the default setup.
#
# Each pinned value was measured under exactly the fit_control settings its
# block uses, and checked for basin robustness first. A pin whose basin is not
# unique fails on another platform's optimizer path for reasons unrelated to any
# code change, which is what the multi-OS check already sees on the phased
# integration fits. For the three cheap fits that check was six jittered starts
# (sd = 0.1), which spread the objective by ~1e-11 relative, matching the four
# references above. The random-recruitment fit is too slow to jitter six times;
# its evidence is instead that the phased and unphased optimizer paths agree to
# 1e-13 relative -- path robustness rather than start-point robustness, which is
# weaker, so treat a failure there as worth investigating before re-pinning.
#
# Dirichlet-multinomial is pinned WITH a prior on theta, and cannot be pinned
# without one. BS2017SS's compositions do not identify the DM theta on their
# own: as theta grows the Dirichlet-multinomial tends to the multinomial, so
# with uninformative comps there is no interior optimum. Bare DM shows exactly
# that -- phased it stops at max|grad| ~ 2 on 'comp_weights' with the
# convergence diagnostic reporting FAIL, unphased it reaches only 2.9e-4, and
# phased and unphased land 3 objective units apart (3.2e-4 relative, ~300x the
# tolerance here), so a pinned value would track the phasing path rather than
# the model. Its jitter spread is 5.7e-8: passable today, but ~2000x looser
# than every other reference here and the first margin to go on another
# platform. The gamma prior below identifies theta instead of masking it --
# theta lands finite and per-fleet (0.54 to 17.0) rather than running to the
# multinomial limit, the gradient falls to 3.9e-11, and the jitter spread
# tightens to 1.7e-11, in line with the rest of this file. The prior path is
# itself untested elsewhere at the fit level, so this pin covers the linkage
# grammar's intercept-prior machinery as well as the DM likelihood.
# =============================================================================

golden_fit_control <- function(...) {
  Rceattle::fit_control(getsd = FALSE, verbose = 0, ...)
}

testthat::test_that("composition tail accumulation reproduces its pinned objective", {
  testthat::skip_on_cran()
  testthat::skip_on_covr()
  testthat::skip_if_not_installed("TMB")

  # Fold ages 1-2 into age 3 and the oldest two ages into age nages-2 on every
  # fleet, so both tails of the composition likelihood are exercised. The bin
  # columns are indices on the fleet's own dimension, so the old-tail bin is
  # taken per species rather than as one number.
  d <- Rceattle::BS2017SS
  spp <- d$fleet_control$Species
  d$fleet_control$Comp_accum_young <- 3L
  d$fleet_control$Comp_accum_old   <- as.integer(d$nages[spp] - 2L)

  fit <- Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 0,
    random_rec = FALSE, msmMode = 0,
    fit_control = golden_fit_control(phase = TRUE, newtonsteps = 3))

  testthat::expect_equal(fit$opt$objective, 6667.9360362185, tolerance = 1e-6)
  testthat::expect_lt(max(abs(fit$obj$gr(fit$opt$par))), 1e-4)
})

testthat::test_that("the MVN and Normal index families reproduce their pinned objective", {
  testthat::skip_on_cran()
  testthat::skip_on_covr()
  testthat::skip_if_not_installed("TMB")

  d <- Rceattle::BS2017SS
  d$fleet_control$Index_distribution    <- "Lognormal"
  d$fleet_control$Index_distribution[4] <- "MVN"      # BT_Pollock, 36 years
  d$fleet_control$Index_distribution[7] <- "Normal"   # EIT_Pollock
  # AnalyticalArith is the q form the covariance likelihood is paired with (see
  # section 8.2b of ceattle.cpp), and the pairing EBS pollock's BTS survey runs.
  # Leaving q freely estimated instead is not a configuration anyone uses, and
  # it converges markedly worse -- see the note above the gradient check below.
  d$fleet_control$Catchability[4]       <- "AnalyticalArith"

  # Compound-symmetric Sigma on the OBSERVATION scale -- the covariance index
  # likelihood works on the natural scale, not the log scale. Built from the
  # fleet's own observations and log sds so it is deterministic and scaled like
  # the data it describes; rho = 0.3 makes the off-diagonals load-bearing, so a
  # change that ignored them would move the objective.
  i <- d$index_data[d$index_data$Fleet_code == 4, ]
  sds <- i$Observation * i$Log_sd
  R <- matrix(0.3, length(sds), length(sds))
  diag(R) <- 1
  Sigma <- diag(sds) %*% R %*% diag(sds)
  dimnames(Sigma) <- list(as.character(i$Year), as.character(i$Year))
  d$index_cov <- list(BT_Pollock = Sigma)

  # Unphased with heavy polishing: phasing leaves this configuration twice as
  # loose (4.8e-5 against 2.5e-5), and the extra Newton steps are what reach a
  # stationary point at all.
  fit <- Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 0,
    random_rec = FALSE, msmMode = 0,
    fit_control = golden_fit_control(phase = FALSE, newtonsteps = 15,
                                     loopnum = 3))

  testthat::expect_equal(fit$opt$objective, 17567.4350624190, tolerance = 1e-6)
  # This fit converges less tightly than the others here (2.5e-5 against
  # ~1e-11), for a reason worth recording rather than tuning away: index_log_q
  # and rec_pars carry EQUAL gradient at the optimum, which is the signature of
  # the q-versus-population-scale ridge -- raising catchability while lowering
  # abundance leaves the predicted index unchanged. Compound symmetry in Sigma
  # then softens precisely that direction, by 1 + (n-1)*rho = 11.5x here, which
  # is why it shows up under the covariance likelihood and not under the
  # lognormal reference. It still clears the same 1e-4 gate as the rest.
  testthat::expect_lt(max(abs(fit$obj$gr(fit$opt$par))), 1e-4)
})

testthat::test_that("Dirichlet-multinomial comps under a theta prior reproduce their pinned objective", {
  testthat::skip_on_cran()
  testthat::skip_on_covr()
  testthat::skip_if_not_installed("TMB")

  d <- Rceattle::BS2017SS
  d$fleet_control$Comp_distribution <- "DirichletMultinomial"

  # Prior-only process: the DM weights have no covariate or random-effect form,
  # so the intercept re-targets the log DM scalar and a gamma prior on the
  # natural theta falls out of the shared intercept-prior loop. gamma(2, 0.5)
  # matches the shape used in test-linkage-composition.R.
  cfun <- Rceattle::build_composition(linkages = list(
    theta_comp = Rceattle::linkage_spec(
      ~ 1, by = ~ fleet,
      priors = list(`(Intercept)` = gamma(2, 0.5)))))

  fit <- Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 0,
    random_rec = FALSE, msmMode = 0, compFun = cfun,
    fit_control = golden_fit_control(phase = FALSE, newtonsteps = 3))

  testthat::expect_equal(fit$opt$objective, 9293.4046676646, tolerance = 1e-6)
  testthat::expect_lt(max(abs(fit$obj$gr(fit$opt$par))), 1e-4)

  # The prior has to identify theta, not merely regularise the objective: an
  # unidentified DM runs theta toward the multinomial limit. Finite, per-fleet
  # values are the evidence that it is pinned by data plus prior.
  theta <- exp(fit$opt$par[names(fit$opt$par) == "comp_weights"])
  testthat::expect_true(all(is.finite(theta)))
  testthat::expect_true(all(theta > 0.1 & theta < 100))
})

testthat::test_that("the default MultinomialAFSC composition likelihood reproduces its pinned objective", {
  testthat::skip_on_cran()
  testthat::skip_on_covr()
  testthat::skip_if_not_installed("TMB")

  # The branch a model gets when Comp_distribution is left out of fleet_control,
  # and the one EBS pollock runs. Phased and unphased agree to ten decimals here,
  # so the pin does not depend on the optimizer path.
  d <- Rceattle::BS2017SS
  d$fleet_control$Comp_distribution <- "MultinomialAFSC"

  fit <- Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 0,
    random_rec = FALSE, msmMode = 0,
    fit_control = golden_fit_control(phase = TRUE, newtonsteps = 3))

  testthat::expect_equal(fit$opt$objective, 4272.2320896464, tolerance = 1e-6)
  testthat::expect_lt(max(abs(fit$obj$gr(fit$opt$par))), 1e-4)
})

testthat::test_that("random recruitment reproduces its pinned objective", {
  testthat::skip_on_cran()
  testthat::skip_on_covr()
  testthat::skip_if_not_installed("TMB")

  # The only golden cover of the Laplace path, which every random_rec = TRUE
  # assessment runs. It is also the slowest fit in the suite at ~4.5 minutes;
  # unphasing does not help (315 s, same optimum), so it stays phased.
  fit <- Rceattle::fit_mod(
    data_list = Rceattle::BS2017SS, inits = NULL, file = NULL,
    estimateMode = 0, random_rec = TRUE, msmMode = 0,
    fit_control = golden_fit_control(phase = TRUE, newtonsteps = 3))

  # as.numeric(): under the Laplace approximation the objective comes back
  # carrying the log-determinant's `logarithm` attribute, which expect_equal
  # compares and would fail on even though the value matches.
  testthat::expect_equal(as.numeric(fit$opt$objective), 10401.6243383277,
                         tolerance = 1e-6)
  # Marginal gradient, the random effects having been integrated out.
  testthat::expect_lt(max(abs(fit$obj$gr(fit$opt$par))), 1e-4)
})
