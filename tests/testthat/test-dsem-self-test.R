# self_test() on a DSEM.
#
# WHAT A SELF TEST MEANS HERE, because it is easy to over-read. With the default
# process = FALSE, sim_mod() draws new OBSERVATIONS and holds the latent states
# at their fitted values, so the test asks whether the model recovers itself
# under fresh observation error -- it says nothing about how well the SEM
# structure is identified across recruitment realizations. process = TRUE
# re-realizes the field from the precision that scored the fit, which is the
# form that does ask that question.
#
# A covariate is the exception to "the process is held": under any family but
# "fixed" it is a latent state observed with error, so its OBSERVATION is
# redrawn whenever the observations are, process error or not.
#
# What these tests pin: the simulated data really is different, the refit
# re-estimates the DSEM rather than echoing the parent back, and the replicates
# converge.

testthat::test_that("self_test() refits a DSEM on simulated data", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("dsem")

  d <- Rceattle::BS2017SS
  yrs <- d$styr:d$endyr
  set.seed(1)
  d$env_data <- data.frame(
    Year = yrs, temp = as.numeric(scale(cumsum(stats::rnorm(length(yrs))))))
  sem <- "
    recdevs1 -> recdevs1, 1, rho_R,     0.3
    temp     -> recdevs1, 0, temp_to_R, 0.2
    recdevs1 <-> recdevs1, 0, sigmaR1,  0.6
    recdevs2 <-> recdevs2, 0, sigmaR2,  0.6
    recdevs3 <-> recdevs3, 0, sigmaR3,  0.6
  "
  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 1,
    random_rec = TRUE, msmMode = 0,
    dsem = Rceattle::build_DSEM(sem = sem, family = "fixed"),
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0))))

  # The simulated data must actually differ, or the "self test" refits the very
  # data set it started from and passes for free.
  sim <- suppressWarnings(suppressMessages(
    Rceattle::sim_mod(fit, simulate = TRUE)))
  testthat::expect_false(isTRUE(all.equal(sim$index_data$Observation,
                                          d$index_data$Observation)))

  # self_test() returns the model list DIRECTLY, named Sim_1..Sim_k -- not
  # wrapped in $Rceattle_list the way retrospective() and jitter() are.
  # The SEM maps out rec_dev, which build_map() would free: that is not a hold the
  # author made, so the redrawn-deviation warning must not name it.
  testthat::expect_true(all(is.na(fit$map$mapFactor$rec_dev)))
  testthat::expect_false(any(c("rec_dev", "R_log_sd") %in%
                               names(Rceattle:::.map_holds(fit))))
  warns <- character()
  st <- withCallingHandlers(
    suppressMessages(
      Rceattle::self_test(fit, nsim = 2, seed = 11, cores = 1, getsd = FALSE)),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  testthat::expect_false(any(grepl("redrew deviations", warns)))
  testthat::expect_gt(length(st), 0L)

  bp <- fit$estimated_params$dsem_beta_z
  for (nm in names(st)) {
    b <- st[[nm]]$estimated_params$dsem_beta_z
    testthat::expect_equal(length(b), length(bp))
    # Re-estimated, not carried over: a refit that silently reused the parent's
    # DSEM would return these unchanged.
    testthat::expect_gt(max(abs(b - bp)), 1e-6)
    testthat::expect_true(is.finite(st[[nm]]$opt$objective))
  }
})

# `process` accepts FALSE/TRUE, "none"/"all"/"dynamics"/"observation", and any
# subset of the five process names. An earlier guard tested isTRUE(process),
# which is FALSE for every character spelling -- so process = "recruitment" and
# process = "all" walked past it and returned a data set whose recruitment was
# never redrawn, silently. Under "all" it was worse than silent: the warning
# enumerated the OTHER un-drawn processes and so read as confirmation that
# recruitment HAD been drawn.

testthat::test_that("redrawing recruitment on a DSEM works by every spelling", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("dsem")

  d <- Rceattle::BS2017SS
  d$env_data <- data.frame(Year = d$styr:d$endyr, BT = 0)
  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 3,
    random_rec = TRUE, msmMode = 0, dsem = Rceattle::build_DSEM(),
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0))))

  # Every spelling that resolves to "redraw recruitment" must reach the GMRF
  # draw. `process` accepts TRUE, "all", "dynamics", "recruitment" and a vector
  # containing it, and isTRUE() is FALSE for all but the first -- so a guard
  # written against the raw argument rather than the resolved state vector lets
  # four of these five through silently.
  drawn <- list()
  for (p in list(TRUE, "all", "dynamics", "recruitment",
                 c("recruitment", "M"))) {
    lab <- paste(deparse(p), collapse = " ")
    set.seed(4)
    sim <- suppressWarnings(suppressMessages(
      Rceattle::sim_mod(fit, simulate = TRUE, process = p)))
    testthat::expect_s3_class(sim, "Rceattle_data")
    drawn[[lab]] <- sim$index_data$Observation
  }

  # ...and the draw is a DRAW that REACHES THE DYNAMICS. Comparing two seeds
  # proves nothing here: the observation draw alone makes them differ, so that
  # assertion passed while the substituted field was being discarded wholesale
  # by a second draw inside the template. Substitute a field that cannot be
  # confused with a draw and check rec_dev moves with it, at one fixed seed.
  obj <- Rceattle:::.sim_obj(fit)
  par <- obj$env$last.par.best
  shifted <- par; shifted[names(par) == "dsem_x_tj"] <-
    par[names(par) == "dsem_x_tj"] + 100
  st <- Rceattle:::.sim_state_codes(TRUE)
  set.seed(4); base <- Rceattle:::.sim_draw(obj, state = st, par = par)
  set.seed(4); moved <- Rceattle:::.sim_draw(obj, state = st, par = shifted)
  testthat::expect_false(isTRUE(all.equal(base$rec_dev_sim, moved$rec_dev_sim)))
  testthat::expect_gt(max(abs(base$rec_dev_sim - moved$rec_dev_sim)), 50)

  # ...and the deviations that generated the data come back, so a recovery
  # statistic has something to score against. Absent, a caller reads it as "no
  # process error was simulated" -- which is what the template's zeroed mask
  # made it look like.
  set.seed(4)
  sim2 <- suppressWarnings(suppressMessages(
    Rceattle::sim_mod(fit, simulate = TRUE, process = TRUE)))
  truth <- attr(sim2, "process_sim")
  testthat::expect_false(is.null(truth))
  testthat::expect_true(all(c("rec_dev", "rec_dev_drawn",
                              "dsem_x_tj", "dsem_x_tj_drawn") %in% names(truth)))
  testthat::expect_gt(sum(truth$rec_dev_drawn), 0L)
  # The mask marks the cells the draw touched; everything else is a fitted
  # value, and scoring recovery over it would report recovery of something that
  # was never simulated.
  testthat::expect_equal(dim(truth$rec_dev_drawn), dim(truth$rec_dev))
  # Against the deviations the FIT used, not estimated_params$rec_dev -- under a
  # DSEM that block is mapped out and sits at zero, so it would only assert the
  # draw is non-zero.
  fitted_dev <- Rceattle:::.sim_draw(
    obj, state = Rceattle:::.sim_state_codes(FALSE))$rec_dev_sim
  testthat::expect_gt(
    max(abs(truth$rec_dev[truth$rec_dev_drawn] -
              fitted_dev[truth$rec_dev_drawn])), 1e-6)

  # ...and the spellings that do NOT ask for recruitment must still work.
  for (p in list(FALSE, "none", "observation", "M")) {
    testthat::expect_error(
      suppressWarnings(suppressMessages(
        Rceattle::sim_mod(fit, simulate = TRUE, process = p))),
      NA, info = paste(deparse(p), collapse = " "))
  }
})

testthat::test_that("a naive DSEM matches a non-DSEM fit at the DEFAULT bias correction", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("dsem")

  # This used to assert a WARNING that R0 was not comparable, because the
  # lognormal bias correction was not applied to the DSEM deviations. It is now,
  # so the warning is gone and the property it warned about is the thing to
  # test. R0 is the number that mattered: SSB absorbs the offset and looked
  # fine, while R0 came out 20-51% low -- and dynamic B0 and the Tier-3 B40%
  # proxy are keyed to R0.
  #
  # Deliberately at the DEFAULT bias_adjust_proc. The other equivalence test
  # sets it FALSE so the two models share a parameterization, which is right for
  # that comparison but leaves the SHIPPED default untested -- and the default
  # is where the defect lived.
  d <- Rceattle::BS2017SS
  d$env_data <- data.frame(Year = d$styr:d$endyr, BT = 0)
  fc <- Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0)
  fd <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 1,
    random_rec = TRUE, msmMode = 0, dsem = Rceattle::build_DSEM(),
    fit_control = fc)))
  fp <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 1,
    random_rec = TRUE, msmMode = 0, fit_control = fc)))

  testthat::expect_equal(fd$opt$objective, fp$opt$objective, tolerance = 1e-6)
  testthat::expect_equal(as.numeric(fd$quantities$R_sd),
                         as.numeric(fp$quantities$R_sd), tolerance = 1e-4)
  # The one that was 20-51% out.
  testthat::expect_equal(as.numeric(fd$quantities$R0)[1:d$nspp],
                         as.numeric(fp$quantities$R0)[1:d$nspp],
                         tolerance = 1e-3)

  # The correction is built on the marginal variance, which must be readable
  # from R or it cannot be checked. For an IID sem it equals R_sd^2 exactly.
  mv <- fd$quantities$dsem_margvar_tj
  testthat::expect_false(is.null(mv))
  testthat::expect_equal(as.numeric(mv[round(nrow(mv) / 2), 1]),
                         as.numeric(fd$quantities$R_sd[1])^2, tolerance = 1e-6)
})

testthat::test_that("remove_F(), reweight_comps() and osa_residuals() keep the DSEM", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("dsem")

  d <- Rceattle::BS2017SS
  d$env_data <- data.frame(Year = d$styr:d$endyr, BT = 0)
  fc <- Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0,
                              bias_adjust_proc = FALSE)
  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 0,
    random_rec = TRUE, msmMode = 0, dsem = Rceattle::build_DSEM(),
    fit_control = fc)))
  bp <- as.numeric(fit$estimated_params$dsem_beta_z)
  # The value the bug produced, so the assertions below name what they exclude.
  START <- 0.707107

  rf <- suppressWarnings(suppressMessages(Rceattle::remove_F(fit)))
  # remove_F() rebuilds at estimateMode = 3 and re-estimates nothing, so the
  # DSEM must come back EXACTLY, not merely close.
  testthat::expect_equal(as.numeric(rf$estimated_params$dsem_beta_z), bp,
                         tolerance = 1e-10)

  # reweight_comps() RE-OPTIMIZES (estimateMode = 0), so a start value is not
  # observable in its answer -- asserting "beta_z is not 0.707107" would pass
  # even with the warm-start bug present, because any start reaches the same
  # optimum. Verified by reproducing the bug: resetting every dsem_* block to the
  # template and refitting returned 0.9888/1.3321/0.7976 anyway. So assert the
  # property that IS at risk: the loop's answer must be the fit implied by the
  # weights it converged on, not merely a converged fit.
  rw <- suppressWarnings(suppressMessages(
    Rceattle::reweight_comps(fit, n_iter = 2, verbose = FALSE)))
  br <- as.numeric(rw$estimated_params$dsem_beta_z)
  testthat::expect_equal(length(br), length(bp))

  d2 <- d
  d2$fleet_control$Comp_weights <- rw$data_list$fleet_control$Comp_weights
  scratch <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d2, inits = NULL, file = NULL, estimateMode = 0,
    random_rec = TRUE, msmMode = 0, dsem = Rceattle::build_DSEM(),
    fit_control = fc)))
  testthat::expect_equal(br, as.numeric(scratch$estimated_params$dsem_beta_z),
                         tolerance = 1e-2)

  # osa_residuals(): one-step-ahead residuals for the DATA. The random effects
  # are integrated out whatever their structure, so a DSEM must give the same
  # index residuals as the non-DSEM fit of the same model.
  fp <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 0,
    random_rec = TRUE, msmMode = 0, fit_control = fc)))
  od <- suppressWarnings(suppressMessages(
    Rceattle::osa_residuals(fit, source = "index", parallel = FALSE)))
  op <- suppressWarnings(suppressMessages(
    Rceattle::osa_residuals(fp, source = "index", parallel = FALSE)))
  testthat::expect_true(all(is.finite(od$residual)))
  testthat::expect_equal(od$residual, op$residual, tolerance = 1e-3)
})
