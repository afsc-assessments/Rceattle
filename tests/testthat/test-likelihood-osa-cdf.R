# One-step-ahead residuals via the conditional CDF (method = "cdf").
#
# Under this method TMB does not treat the observation as a free variable: it
# asks the template for log F(x) and log(1 - F(x)) of each conditional, gated by
# keep.cdf_lower / keep.cdf_upper, and returns qnorm(F(x)). No conditional mean
# is formed, so nothing can leave the support of a bounded observation.
#
# The conditional for a composition bin is a binomial, written through the
# regularized incomplete beta so it is defined at the fractional counts
# composition data carry:  P(X <= x) = 1 - I_p(x + 1, n - x).  The aggregate
# series are normal (lognormal catch and index on the log scale, "Normal" and
# "MVN" on their own), so theirs is pnorm.
#
# The distributional acceptance test -- residualizing simulated data at the
# parameters that generated it, where the answer is exactly N(0, 1) -- lives in
# tools/verify/verify-osa-cdf.R, which is far too slow for the suite. What is
# here is everything that pins the implementation: the identity the CDF is built
# on, the bin alignment through the keep reorder, the invariance of the
# objective, and the two families that cannot use this method.

# The continuous binomial CDF the template uses, written in R. Used below to
# check the template's own conditional decomposition bin by bin.
.cdf_binom <- function(x, n, p) 1 - stats::pbeta(p, x + 1, n - x)

# Conditional-binomial CDF of every bin of one composition, in the same
# sequential order dmultinom_osa() walks: bin i against everything still
# unplaced, at probability p_i / (1 - p_used).
.cdf_comp_ref <- function(x, p) {
  p <- p / sum(p)
  n <- sum(x)
  used <- 0
  out <- numeric(length(x))
  for (i in seq_along(x)) {
    if (i == length(x)) { out[i] <- 1; break }        # fixed by sum-to-N
    out[i] <- .cdf_binom(x[i], n, p[i] / (1 - used))
    n    <- n - x[i]
    used <- used + p[i]
  }
  out
}


testthat::test_that("the continuous binomial CDF agrees with pbinom at whole counts", {
  # The identity the composition CDF rests on, checked against R's own binomial
  # CDF where both are defined. Nothing model-specific: if this fails the
  # formula is wrong, not the plumbing.
  # Absolute, not relative. `1 - pbeta(p, x + 1, n - x)` cancels to exactly 0
  # once the lower tail falls below about 1e-16, where pbinom() still computes
  # it; that is the mechanism behind the lower-tail censoring the method carries
  # (the template floors the CDF at 4e-16 anyway, so a bin down there is
  # reported at the ceiling either way). Everything above the floor must agree.
  for (n in c(5, 40, 500)) {
    for (p in c(0.001, 0.05, 0.4, 0.9, 0.999)) {
      x <- unique(round(seq(0, n - 1, length.out = 12)))
      testthat::expect_lt(max(abs(.cdf_binom(x, n, p) - stats::pbinom(x, n, p))),
                          1e-12, label = paste("n =", n, "p =", p))
    }
  }
  # And that it interpolates monotonically between them, which is what makes it
  # usable at a fractional count.
  f <- .cdf_binom(seq(0, 9, by = 0.05), 10, 0.3)
  testthat::expect_true(all(diff(f) >= 0))
})


testthat::test_that("osa_mode 2 leaves the objective and every likelihood slot alone", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # cdf_lower and cdf_upper are zero except inside a oneStepPredict(method =
  # "cdf") call, so adding the CDF terms must not move anything. Checked slot by
  # slot rather than on the total, so two slots cannot cancel.
  dat <- make_test_data(nyrs = 10, nages = 6, seed = 11)
  fit <- Rceattle::fit_mod(dat, file = NULL, estimateMode = 3, msmMode = 0,
                           fit_control = fit_control(phase = FALSE, verbose = 0))
  osa_dat <- Rceattle:::build_osa_data(fit$obj$env$data, build_osa = TRUE)

  o1 <- Rceattle:::.osa_build_obj(fit, osa_dat, osa_mode = 1)
  o2 <- Rceattle:::.osa_build_obj(fit, osa_dat, osa_mode = 2)
  testthat::expect_equal(o2$env$data$osa_mode, 2)

  testthat::expect_identical(o2$fn(o2$par), o1$fn(o1$par))
  testthat::expect_identical(o2$report(o2$par)$jnll_comp,
                             o1$report(o1$par)$jnll_comp)
})


testthat::test_that("the cdf gates survive the keep reorder", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # `dmultinom_osa()` reorders its bins so kept ones come first (osa_order), and
  # keep.cdf_lower / keep.cdf_upper are separate members of the data_indicator
  # that have to be reordered with them. Getting that wrong reads the CDF off
  # the wrong bin and still returns a plausible number.
  #
  # It cannot be reached through osa_residuals(). osa_order() branches on keep,
  # keep is a PARAMETER under oneStepPredict, so the branch is taped once at the
  # initial values -- 1 on every bin that will ever be kept, 0 on the rest -- and
  # within a composition that is 1 everywhere but the last bin. The order is
  # always the identity and `conditional =` does not change it, because a
  # conditioned bin is also kept at tape time.
  #
  # So supply the keep matrix directly, as constant data, with a pattern the
  # package never produces: keep only bin 3. osa_order() then puts bin 3 first
  # and the CDF gate on bin 1 has to travel to the second slot.
  dat <- make_test_data(nyrs = 10, nages = 6, seed = 18)
  fit <- suppressWarnings(Rceattle::fit_mod(
    dat, file = NULL, estimateMode = 1, msmMode = 0,
    fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE)))
  osa_dat <- Rceattle:::build_osa_data(fit$obj$env$data, build_osa = TRUE)
  ctl <- osa_dat$obs_ctl
  testthat::skip_if(!any(ctl$source == "comp"), "no composition observations")

  grp <- ctl[ctl$source == "comp", ]
  grp <- grp[grp$group_id == grp$group_id[1], ]
  grp <- grp[order(grp$bin_index), ]
  n <- nrow(grp)
  testthat::skip_if(n < 4, "composition too short to permute")

  x <- as.numeric(osa_dat$obsvec[grp$obs_pos + 1L])
  p <- as.numeric(fit$quantities$comp_hat[grp$data_row[1], seq_len(n)]) +
    osa_dat$comp_offset

  # Rebuild in OSA-CDF mode with `keep` as a fully-mapped (constant) parameter,
  # so the taped branch sees exactly the pattern given here.
  obj0 <- Rceattle:::.osa_build_obj(fit, osa_dat, osa_mode = 2)
  nobs <- length(osa_dat$obsvec)
  at   <- grp$obs_pos + 1L                    # this composition's obsvec slots
  objective <- function(k, cl, cu) {
    kp <- cbind(k, cl, cu)
    o <- TMB::MakeADFun(
      data       = obj0$env$data,
      parameters = c(obj0$env$parList(par = obj0$env$last.par.best),
                     list(keep = kp)),
      map        = c(obj0$env$map, list(keep = factor(rep(NA, length(kp))))),
      random     = if (length(obj0$env$random))
        unique(names(obj0$env$par)[obj0$env$random]) else NULL,
      DLL = fit$TMBfilename, silent = TRUE)
    o$fn(o$par)
  }
  z <- rep(0, nobs)
  keep3 <- z; keep3[at[3]] <- 1                # keep ONLY bin 3
  lo <- z; lo[at[1]] <- 1                      # lower-CDF gate on bin 1
  hi <- z; hi[at[1]] <- 1
  base <- objective(keep3, z, z)
  Fx <- 1 / (1 + exp((objective(keep3, lo, z) - base) -
                       (objective(keep3, z, hi) - base)))

  # Bin 3 is sorted to the front, so bin 1 is read second: against the count and
  # the probability still unspent after bin 3.
  o   <- c(3L, setdiff(seq_len(n), 3L))
  ref <- .cdf_comp_ref(x[o], p[o])[which(o == 1L)]
  testthat::expect_equal(Fx, ref, tolerance = 1e-7)

  # ... and that is a different number from the unpermuted answer, so a gate
  # left behind by the reorder would fail here.
  testthat::expect_false(isTRUE(all.equal(ref, .cdf_comp_ref(x, p)[1],
                                          tolerance = 1e-4)))
})


testthat::test_that("splitting a call into groups conditions on what came before it", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # Compositions need `discrete = TRUE` and the aggregate series do not, and
  # oneStepPredict takes one setting per call, so they go in separate calls.
  # Every call is handed the rows earlier in the sequence as `conditional`;
  # without that, oneStepPredict marks them unconditional and zeroes their data
  # terms, which on a random-effects model changes the conditional distribution
  # of the latent states and moves the residuals.
  dat <- make_test_data(nyrs = 14, nages = 8, seed = 31)
  fit <- suppressWarnings(Rceattle::fit_mod(
    dat, file = NULL, estimateMode = 1, msmMode = 0, random_rec = TRUE,
    fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE)))
  testthat::skip_if(length(fit$obj$env$random) == 0, "no random effects")

  osa_dat <- Rceattle:::build_osa_data(fit$obj$env$data, build_osa = TRUE)
  sel <- osa_dat$obs_ctl[osa_dat$obs_ctl$source %in% c("index", "catch", "comp") &
                           !osa_dat$obs_ctl$is_last_bin, ]
  sel <- sel[order(match(sel$source, c("index", "catch", "comp")), sel$year,
                   sel$fleet_code, sel$bin_index), ]
  pos <- sel$obs_pos + 1L
  agg <- which(sel$source != "comp")
  cmp <- which(sel$source == "comp")
  testthat::skip_if(length(cmp) < 3 || length(agg) < 3, "not enough of each type")

  obj <- Rceattle:::.osa_build_obj(fit, osa_dat, osa_mode = 2)
  osp <- function(...) TMB::oneStepPredict(obj, "obsvec", "keep", method = "cdf",
                                           discrete = FALSE, parallel = FALSE,
                                           trace = FALSE, ...)$residual
  one   <- osp(subset = pos)[cmp]
  split <- osp(subset = pos[cmp], conditional = pos[agg])
  naive <- osp(subset = pos[cmp])

  # The split is invisible ...
  testthat::expect_equal(split, one, tolerance = 1e-9)
  # ... and would not have been without `conditional`, so this is not vacuous.
  testthat::expect_gt(max(abs(naive - one)), 1e-4)
})


testthat::test_that("the template's conditional CDF is the right one, bin by bin", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # The trap this guards. keep.cdf_lower / keep.cdf_upper are separate members
  # of the data_indicator, and dmultinom_osa() reorders its bins (osa_order) --
  # so they have to be reordered with it. A gate left on the wrong bin still
  # returns plausible residuals; it just returns the wrong bin's. Comparing
  # TMB's Fx against the conditional decomposition computed independently in R
  # catches that, and catches a wrong binomial size or conditional probability
  # with it.
  dat <- make_test_data(nyrs = 10, nages = 6, seed = 12)
  fit <- suppressWarnings(Rceattle::fit_mod(
    dat, file = NULL, estimateMode = 1, msmMode = 0,
    fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE)))

  osa_dat <- Rceattle:::build_osa_data(fit$obj$env$data, build_osa = TRUE)
  ctl <- osa_dat$obs_ctl
  testthat::skip_if(!any(ctl$source == "comp"), "no composition observations")

  # One composition (all its bins, including the last), in bin order.
  grp <- ctl[ctl$source == "comp", ]
  grp <- grp[grp$group_id == grp$group_id[1], ]
  grp <- grp[order(grp$bin_index), ]
  n   <- nrow(grp)

  x <- as.numeric(osa_dat$obsvec[grp$obs_pos + 1L])
  p <- as.numeric(fit$quantities$comp_hat[grp$data_row[1], seq_len(n)]) +
    osa_dat$comp_offset
  ref <- .cdf_comp_ref(x, p)

  obj <- Rceattle:::.osa_build_obj(fit, osa_dat, osa_mode = 2)
  res <- TMB::oneStepPredict(obj, "obsvec", "keep", method = "cdf",
                             subset = grp$obs_pos[-n] + 1L, discrete = FALSE,
                             parallel = FALSE, trace = FALSE)

  testthat::expect_equal(res$Fx, ref[-n], tolerance = 1e-7)
  # px is the conditional point mass, so it must reproduce the density too.
  testthat::expect_true(all(res$px > 0 & res$px <= 1))
  # ... and the residual is the plain probability integral transform of Fx.
  testthat::expect_equal(res$residual, stats::qnorm(res$Fx), tolerance = 1e-10)
})


testthat::test_that("cdf reproduces the Gaussian residual on the aggregate series", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # Catch and index are genuinely normal on the scale obsvec stores them in, so
  # the Gaussian method's conditional mode IS the mean and both methods must
  # return the same number. Disagreement here is an implementation bug, not a
  # better method. The tolerance is set by oneStepGaussianOffMode's own
  # numerics, not by ours.
  dat <- make_test_data(nyrs = 12, nages = 6, seed = 13)
  fit <- suppressWarnings(Rceattle::fit_mod(
    dat, file = NULL, estimateMode = 1, msmMode = 0,
    fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE)))

  g <- suppressWarnings(suppressMessages(Rceattle::osa_residuals(
    fit, source = c("index", "catch"), parallel = FALSE)))
  c_ <- suppressWarnings(suppressMessages(Rceattle::osa_residuals(
    fit, source = c("index", "catch"), method = "cdf", parallel = FALSE)))

  testthat::expect_equal(nrow(g), nrow(c_))
  testthat::expect_equal(g$residual, c_$residual, tolerance = 1e-4)
  testthat::expect_true(all(is.finite(c_$residual)))
})


testthat::test_that("cdf reports no conditional mean, and says so on the object", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  dat <- make_test_data(nyrs = 10, nages = 6, seed = 14)
  fit <- suppressWarnings(Rceattle::fit_mod(
    dat, file = NULL, estimateMode = 1, msmMode = 0,
    fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE)))

  o <- suppressWarnings(suppressMessages(Rceattle::osa_residuals(
    fit, source = c("index", "comp"), method = "cdf", parallel = FALSE)))

  # No Gaussian mean is formed, so predicted and sd are NA rather than a value
  # borrowed from somewhere it does not mean the same thing.
  testthat::expect_true(all(is.na(o$predicted)))
  testthat::expect_true(all(is.na(o$sd)))
  # ... which is exactly why an expected count can no longer go negative.
  testthat::expect_true(all(!is.na(o$observed)))
  testthat::expect_true(all(is.finite(o$residual)))
  # A CDF read in double precision cannot resolve past the last value below 1.
  testthat::expect_lte(max(abs(o$residual)), 8.042)

  # `discrete` resolves per method: the randomized quantile residual under cdf,
  # and the previous behaviour under everything else.
  testthat::expect_true(attr(o, "discrete"))
  g <- suppressWarnings(suppressMessages(Rceattle::osa_residuals(
    fit, source = "comp", parallel = FALSE)))
  testthat::expect_false(attr(g, "discrete"))
  testthat::expect_identical(attr(o, "method"), "cdf")

  # Randomized, so the seed moves the composition residuals and nothing else.
  o2 <- suppressWarnings(suppressMessages(Rceattle::osa_residuals(
    fit, source = c("index", "comp"), method = "cdf", parallel = FALSE,
    seed = 4321)))
  testthat::expect_false(isTRUE(all.equal(o$residual[o$source == "comp"],
                                          o2$residual[o2$source == "comp"])))
  testthat::expect_equal(o$residual[o$source == "index"],
                         o2$residual[o2$source == "index"])
})


testthat::test_that("a Dirichlet-multinomial fleet is moved off cdf, loudly", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # The beta-binomial has no closed-form CDF, so the template supplies no CDF
  # term for a D-M composition. That does not fail loudly on its own -- both
  # tails come back equal, giving Fx = 0.5 and a residual of exactly 0 for every
  # bin -- so osa_residuals() has to route those fleets elsewhere. This is the
  # test that the routing happens and that the silent-zero case cannot ship.
  dat <- make_test_data(nyrs = 10, nages = 6, seed = 15)
  dat$fleet_control$Comp_distribution <- "DirichletMultinomial"
  fit <- suppressWarnings(Rceattle::fit_mod(
    dat, file = NULL, estimateMode = 1, msmMode = 0,
    fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE)))

  testthat::expect_message(
    o <- suppressWarnings(Rceattle::osa_residuals(
      fit, source = "comp", method = "cdf", parallel = FALSE)),
    "Dirichlet-multinomial")

  m <- attr(o, "method")
  testthat::expect_equal(unname(m[["default"]]), "cdf")
  testthat::expect_equal(unname(m[["DirichletMultinomial"]]),
                         "oneStepGaussianOffMode")
  # The fallback ran a real method: these are exactly the residuals the package
  # default gives, not the identical-zero a missing CDF term would produce. (An
  # sd threshold would not do -- this fixture's D-M weight runs away, so the fit
  # is near-exact and the residuals are legitimately tiny.)
  g <- suppressWarnings(suppressMessages(Rceattle::osa_residuals(
    fit, source = "comp", parallel = FALSE)))
  testthat::expect_equal(o$residual, g$residual)
  testthat::expect_false(all(abs(o$residual) < 1e-8))

  # The fallback runs a Gaussian method, which DOES return a conditional mode --
  # and it is one that can go negative. It is blanked, so `predicted` means the
  # same thing (nothing) on every row of a "cdf" object.
  testthat::expect_true(all(is.na(o$predicted)))
  testthat::expect_true(all(is.na(o$sd)))
  # Every fleet here is Dirichlet-multinomial, so every composition went to the
  # Gaussian fallback and NOTHING was randomized -- the flag has to say so. It
  # reports what was actually done, not what `discrete` resolved to.
  testthat::expect_false(attr(o, "discrete"))
  # `all = FALSE` here would assert only that NOT EVERY line matches, which is
  # true of any multi-line output -- the assertion has to be that no line does.
  testthat::expect_false(any(grepl("randomized", utils::capture.output(print(o)))))
})


testthat::test_that("cdf residualizes a predator diet composition", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # Diet goes through the same conditional-binomial decomposition as an age
  # composition, from its own call site in the C++ (the stomach likelihood), and
  # is the one source `source = "all"` reaches that nothing else here covers.
  # Its "counts" are prey-weight proportion * number of stomachs, so this is also
  # where discrete = TRUE has the least to do with counting anything -- the model
  # simulates them as multinomial draws, which is what the transform reads.
  nyrs <- 10L; nspp <- 2L
  Fmort  <- c(seq(0.02, 0.3, length.out = nyrs / 2),
              seq(0.3, 0.05, length.out = nyrs / 2))
  Fmort2 <- seq(0.02, 0.3, length.out = nyrs)
  sim <- make_msm_test_data(
    years   = seq_len(nyrs),
    Fmort   = matrix(c(Fmort, Fmort2), nspp, nyrs, byrow = TRUE),
    gam_a   = c(1, 0.1), gam_b = rep(0.15, nspp),
    log_phi = matrix(c(-5, 0.5, -10, -2), nspp, nspp, byrow = TRUE))

  msm <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    sim$data_list, file = NULL, estimateMode = 1, msmMode = 1, suitMode = 4,
    niter = 3, initMode = "NonEquilibrium",
    fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE))))
  ctl <- Rceattle:::build_osa_data(msm$obj$env$data, build_osa = TRUE)$obs_ctl
  testthat::skip_if(!any(ctl$source == "diet"), "no diet observations")

  d <- suppressWarnings(suppressMessages(Rceattle::osa_residuals(
    msm, source = "diet", method = "cdf", parallel = FALSE)))
  testthat::expect_gt(nrow(d), 0)
  testthat::expect_true(all(is.finite(d$residual)))
  testthat::expect_true(all(is.na(d$predicted)))
  # Not the degenerate 0 a missing CDF term would give for every bin.
  testthat::expect_gt(stats::sd(d$residual), 0.05)
  testthat::expect_lte(max(abs(d$residual)), 8.042)
})


testthat::test_that("cdf residualizes a multinomial composition where the Gaussian mean goes negative", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # Issue #108: at a small composition sample size the Gaussian conditional mean
  # lands below zero on the sparse bins, which is impossible for an expected
  # count and biases those rows' residuals positive. Under cdf no mean is formed
  # at all, so the whole failure mode is out of reach -- and the residuals are
  # still finite and non-degenerate, which "no negative values" alone would not
  # tell you.
  # A composition sample size of 1 puts every bin count near zero on this
  # fixture, which is what drives the conditional mode below the support: all 10
  # of its expected counts come back negative under the Gaussian method.
  dat <- make_test_data(nyrs = 10, nages = 6, seed = 16)
  dat$comp_data$Sample_size <- 1
  fit <- suppressWarnings(Rceattle::fit_mod(
    dat, file = NULL, estimateMode = 1, msmMode = 0,
    fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE)))

  g <- suppressWarnings(suppressMessages(Rceattle::osa_residuals(
    fit, source = "comp", parallel = FALSE)))
  testthat::expect_true(any(g$predicted < 0, na.rm = TRUE))

  c_ <- suppressWarnings(suppressMessages(Rceattle::osa_residuals(
    fit, source = "comp", method = "cdf", parallel = FALSE)))
  testthat::expect_true(all(is.na(c_$predicted)))
  testthat::expect_true(all(is.finite(c_$residual)))
  testthat::expect_gt(stats::sd(c_$residual), 0.1)
})


testthat::test_that("cdf residualizes conditional age-at-length", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # CAAL is a second composition source with its own call site and its own
  # likelihood-family vector (`caal_ll_type`), and nothing else here reaches it.
  dat <- make_test_data(nyrs = 10, nages = 6, seed = 51, growth = "vonBertalanffy")
  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    dat, file = NULL, estimateMode = 1, msmMode = 0,
    growthFun = build_growth(fun = "vonBertalanffy"),
    fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE))))

  o <- suppressWarnings(suppressMessages(Rceattle::osa_residuals(
    fit, source = "caal", method = "cdf", parallel = FALSE)))
  testthat::expect_gt(nrow(o), 0)
  testthat::expect_setequal(unique(o$source), "caal")
  testthat::expect_true(all(is.finite(o$residual)))
  testthat::expect_true(all(is.na(o$predicted)))
  testthat::expect_gt(stats::sd(o$residual), 0.05)   # not the degenerate 0
  testthat::expect_lte(max(abs(o$residual)), 8.042)
})


testthat::test_that("cdf residualizes a state-space covariate observation", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # The `ecov` source: a QAR1 catchability linkage whose AR1 latent is measured
  # by an observed covariate. Its CDF sits in the linkage block, not with the
  # aggregate series, and it is the only data term inside a jnll_comp row whose
  # name says "random effects".
  dat <- make_test_data(nyrs = 12, nages = 6, seed = 52)
  dat$fleet_control$Catchability[1] <- "Estimated"
  dat$env_data <- data.frame(Year = dat$styr:dat$endyr,
                             qcov = as.numeric(scale(seq_len(12))))
  qf <- build_catchability(linkages = list(
    q = linkage_spec(~ ar1(1 | Year), by = ~ fleet, fleet = 1,
                     observe = "qcov", obs_sd = 0.3)))
  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    dat, file = NULL, qFun = qf, estimateMode = 1, msmMode = 0,
    fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE))))

  o <- suppressWarnings(suppressMessages(Rceattle::osa_residuals(
    fit, source = "ecov", method = "cdf", parallel = FALSE)))
  testthat::expect_equal(nrow(o), 12L)
  testthat::expect_setequal(unique(o$source), "ecov")
  testthat::expect_true(all(is.finite(o$residual)))
  testthat::expect_true(all(is.na(o$predicted)))
  testthat::expect_gt(stats::sd(o$residual), 0.05)

  # It does NOT match a Gaussian method here, and that is expected rather than a
  # defect: the latent state is a random effect, so oneStepPredict integrates the
  # CDF over it by Laplace, and that integrand is a Gaussian times a sigmoid
  # rather than a density. For a Gaussian observation the Gaussian methods
  # integrate a density and are exact -- see ?osa_residuals. Pinned loosely, as
  # a record that the two differ by an amount of this order.
  g <- suppressWarnings(suppressMessages(Rceattle::osa_residuals(
    fit, source = "ecov", method = "oneStepGaussian", parallel = FALSE)))
  d <- max(abs(o$residual - g$residual))
  testthat::expect_gt(d, 1e-3)
  testthat::expect_lt(d, 1)
})


testthat::test_that("a failed Laplace solve does not cost every observation after it", {
  # What this pins is the recovery MECHANISM -- splice, progress guard, bounds,
  # and which message it emits -- not that the mechanism fixes the failure it was
  # written for. It does not: on BS2017SS with random recruitment, 1879 of 4538
  # composition residuals are non-finite under cdf, and redoing the tail leaves
  # 1879. See .osa_retry_tail() for why (depth of conditioning, not a poisoned
  # warm start) and ?osa_residuals for the limitation.
  #
  # Unit-tested on stubs rather than on a model: reproducing the failure needs a
  # few thousand observations over a hundred random effects, which is hours.
  mk <- function(resid) data.frame(.row = seq_along(resid), observed = 1,
                                   predicted = NA_real_, sd = NA_real_,
                                   residual = resid)

  # A clean tail is recovered in one pass, and the recovered values land in the
  # right rows.
  res <- mk(c(1:5, rep(NaN, 5)))
  testthat::expect_message(
    got <- Rceattle:::.osa_retry_tail(res, function(from) mk(rep(99, 11 - from))),
    "recovering 5 of 5")
  testthat::expect_equal(got$residual, c(1:5, rep(99, 5)))

  # Partial progress is kept and retried: the first pass fixes rows 6-8, the
  # second the rest.
  n_call <- 0L
  res <- mk(c(1:5, rep(NaN, 5)))
  got <- suppressMessages(Rceattle:::.osa_retry_tail(res, function(from) {
    n_call <<- n_call + 1L
    out <- rep(if (n_call == 1L) NaN else 7, 11 - from)
    if (n_call == 1L) out[1:3] <- 50
    mk(out)
  }))
  testthat::expect_true(all(is.finite(got$residual)))
  testthat::expect_equal(got$residual, c(1:5, 50, 50, 50, 7, 7))

  # A retry that makes no progress stops instead of looping to the bound, and the
  # rows it could not fix stay non-finite rather than being silently dropped.
  n_call <- 0L
  res <- mk(c(1, NaN, NaN))
  testthat::expect_message(
    got <- Rceattle:::.osa_retry_tail(res, function(from) {
      n_call <<- n_call + 1L
      mk(rep(NaN, 4 - from))
    }),
    "recovered none")          # and it must NOT claim to have recovered anything
  testthat::expect_equal(n_call, 1L)
  testthat::expect_equal(sum(!is.finite(got$residual)), 2L)

  # Nothing to do is a no-op, and says nothing.
  clean <- mk(1:4)
  testthat::expect_no_message(
    testthat::expect_identical(Rceattle:::.osa_retry_tail(clean, function(from)
      stop("must not be called")), clean))
})


testthat::test_that("an unknown method is rejected before any model is built", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  dat <- make_test_data(nyrs = 8, nages = 5, seed = 17)
  fit <- Rceattle::fit_mod(dat, file = NULL, estimateMode = 3, msmMode = 0,
                           fit_control = fit_control(phase = FALSE, verbose = 0))
  fit$data_list$estimateMode <- 1
  testthat::expect_error(
    Rceattle::osa_residuals(fit, source = "index", method = "cfd"),
    "should be one of")
})


testthat::test_that("`discrete` is validated before any model is built", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # NULL means "choose per method", so the check cannot simply be is.logical():
  # a non-logical has to be refused rather than resolving to one by coercion.
  dat <- make_test_data(nyrs = 8, nages = 5, seed = 17)
  fit <- Rceattle::fit_mod(dat, file = NULL, estimateMode = 3, msmMode = 0,
                           fit_control = fit_control(phase = FALSE, verbose = 0))
  fit$data_list$estimateMode <- 1
  for (bad in list("yes", 1L, c(TRUE, FALSE), NA)) {
    testthat::expect_error(
      Rceattle::osa_residuals(fit, source = "index", method = "cdf",
                              discrete = bad),
      "TRUE, FALSE, or NULL")
  }
})


testthat::test_that("turning `discrete` off under cdf says what it costs", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # discrete = FALSE leaves the step in qnorm(F(x)) for a count observation,
  # which is the worst-calibrated of the three composition options measured
  # (mean +0.610, sd 1.262 against a known standard normal). It is a legitimate
  # thing to ask for -- the verification harness does -- so it runs, but it is
  # announced, like every other resolution in the function that moves numbers.
  dat <- make_test_data(nyrs = 10, nages = 6, seed = 18)
  fit <- suppressWarnings(Rceattle::fit_mod(
    dat, file = NULL, estimateMode = 1, msmMode = 0,
    fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE)))

  testthat::expect_message(
    suppressWarnings(Rceattle::osa_residuals(
      fit, source = "comp", method = "cdf", discrete = FALSE, parallel = FALSE)),
    "discrete = FALSE")
  # Not on the default, which already is TRUE there ...
  testthat::expect_no_message(suppressWarnings(Rceattle::osa_residuals(
    fit, source = "comp", method = "cdf", parallel = FALSE)))
  # ... nor on a Gaussian method, whose compositions were always continuous.
  # Scoped to the discrete notice: leaving `method` unset also raises the
  # default-method notice below, which is a different message.
  testthat::expect_no_message(suppressWarnings(Rceattle::osa_residuals(
    fit, source = "comp", discrete = FALSE, parallel = FALSE)),
    message = "discrete = FALSE")
})

# The package default is biased on composition data by its own scoring table, so
# a caller who never chose a method is told which one they got. Naming any
# method, including the default, is taken as a choice and stays silent.
testthat::test_that("an unset method says so on composition sources", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  dat <- make_test_data(nyrs = 10, nages = 6, seed = 18)
  fit <- suppressWarnings(Rceattle::fit_mod(
    dat, file = NULL, estimateMode = 1, msmMode = 0,
    fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE)))
  testthat::expect_message(
    suppressWarnings(Rceattle::osa_residuals(fit, source = "comp", parallel = FALSE)),
    "biased on composition data")
  testthat::expect_no_message(
    suppressWarnings(Rceattle::osa_residuals(
      fit, source = "comp", method = "oneStepGaussianOffMode", parallel = FALSE)),
    message = "biased on composition data")
  # Aggregate series are Gaussian, so the default is not the biased choice there.
  testthat::expect_no_message(
    suppressWarnings(Rceattle::osa_residuals(fit, source = "index", parallel = FALSE)),
    message = "biased on composition data")
})


testthat::test_that("one Dirichlet-multinomial fleet does not take the others off cdf", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # The fallback is per OBSERVATION, not per call. A whole-model decision would
  # be wrong either way here: routing everything to the Gaussian method throws
  # away the residuals this method exists to give, and routing everything to cdf
  # returns exactly 0 for the D-M bins, since a missing CDF term leaves the two
  # tails equal. Both fleets are in one obsvec, so this is the case that says
  # which it does.
  dat <- make_test_data(nyrs = 10, nages = 6, seed = 15)
  dat$fleet_control$Comp_distribution[dat$fleet_control$Fleet_code == 1] <-
    "DirichletMultinomial"
  fit <- suppressWarnings(Rceattle::fit_mod(
    dat, file = NULL, estimateMode = 1, msmMode = 0,
    fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE)))

  testthat::expect_message(
    o <- suppressWarnings(Rceattle::osa_residuals(
      fit, source = "comp", method = "cdf", parallel = FALSE)),
    "Dirichlet-multinomial")
  g <- suppressWarnings(suppressMessages(Rceattle::osa_residuals(
    fit, source = "comp", parallel = FALSE)))

  dm  <- o$fleet == 1
  mn  <- o$fleet == 2
  testthat::expect_true(any(dm) && any(mn))
  # The D-M rows are the package default's, to the bit.
  testthat::expect_equal(o$residual[dm], g$residual[dm])
  # The multinomial rows are NOT: they went through the conditional CDF, and
  # they are not the identical zero a missing CDF term would have given.
  testthat::expect_false(isTRUE(all.equal(o$residual[mn], g$residual[mn])))
  testthat::expect_false(all(abs(o$residual[mn]) < 1e-8))

  # The attributes report the split rather than one global answer.
  testthat::expect_equal(unname(attr(o, "method")[["default"]]), "cdf")
  testthat::expect_equal(unname(attr(o, "method")[["DirichletMultinomial"]]),
                         "oneStepGaussianOffMode")
  d <- attr(o, "discrete")
  testthat::expect_true(unname(d[["default"]]))
  testthat::expect_false(unname(d[["DirichletMultinomial"]]))
})
