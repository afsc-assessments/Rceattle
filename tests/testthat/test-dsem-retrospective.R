# retrospective() on a DSEM.
#
# A peel does not shorten the model: retrospective() sets endyr_peel and turns
# off DATA after it, so the latent states still span every year. They stay free
# and in `random`, and the Laplace approximation integrates the peeled-year
# states out against the GMRF prior. That is the peeled marginal likelihood.
#
# Mirroring what rec_dev does -- zero the tail, map it out -- would be wrong.
# Pinning is inert for an INDEPENDENT deviate, but DSEM states are coupled
# through the RAM, so pinned zeros stay in the quadratic form and shrink the
# terminal retained state (1 + rho^2 inflated precision). Terminal recruitment
# drives terminal SSB, so Mohn's rho would measure the peel's own artifact.
#
# The sem below is LAGGED and CARRIES A COVARIATE on purpose. A default
# build_DSEM() sem is IID and has no covariate column, and can detect neither
# failure: pinning is inert without coupling, and a covariate bug needs one.

dsem_retro_data <- function() {
  d <- Rceattle::BS2017SS
  yrs <- d$styr:d$endyr
  set.seed(1)
  d$env_data <- data.frame(
    Year = yrs, temp = as.numeric(scale(cumsum(stats::rnorm(length(yrs))))))
  d
}

# No `temp <-> temp`: the builder adds V[temp] itself, and giving an exogenous
# variable two variance terms makes the precision singular and the objective NaN.
DSEM_RETRO_SEM <- "
  recdevs1 -> recdevs1, 1, rho_R,     0.3
  temp     -> recdevs1, 0, temp_to_R, 0.2
  recdevs1 <-> recdevs1, 0, sigmaR1,  0.6
  recdevs2 <-> recdevs2, 0, sigmaR2,  0.6
  recdevs3 <-> recdevs3, 0, sigmaR3,  0.6
"

testthat::test_that("a DSEM peel reproduces a model that genuinely ends at endyr_peel", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("dsem")

  d  <- dsem_retro_data()
  fc <- Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0)
  mk <- function() Rceattle::build_DSEM(sem = DSEM_RETRO_SEM, family = "fixed")

  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 1,
    random_rec = TRUE, msmMode = 0, dsem = mk(), fit_control = fc)))

  retro <- suppressWarnings(suppressMessages(
    Rceattle::retrospective(fit, peels = 1, cores = 1, getsd = FALSE)))
  ml <- retro$Rceattle_list
  # rev(c(list(Rceattle), peels)): the LAST element is the input object, so an
  # "unpeeled peel matches the parent" assertion would compare a fit to itself.
  testthat::expect_equal(length(ml), 2L)
  pk <- ml[[1]]
  endyr_peel <- pk$data_list$endyr_peel
  testthat::expect_equal(endyr_peel, d$endyr - 1L)

  # A model that genuinely ends at endyr_peel: the peeled year does not exist.
  # Integrating a state out must equal never having had it, for a Gaussian
  # process -- so the retained years must match. A peel that PINNED its
  # peeled-year state would inject information this model does not have.
  dt <- d
  dt$index_data <- dt$index_data[dt$index_data$Year <= endyr_peel, ]
  dt$comp_data  <- dt$comp_data[dt$comp_data$Year <= endyr_peel, ]
  dt$catch_data$Catch[dt$catch_data$Year > endyr_peel] <- 0
  dt$endyr <- endyr_peel
  dt$projyr <- endyr_peel
  dt$env_data <- dt$env_data[dt$env_data$Year <= endyr_peel, ]

  ends_at_peel <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = dt, inits = NULL, file = NULL, estimateMode = 1,
    random_rec = TRUE, msmMode = 0, dsem = mk(), fit_control = fc)))

  keep <- seq_len(endyr_peel - d$styr + 1L)
  sp <- 1L   # the species carrying the lag and the covariate
  R_peel <- pk$quantities$R[sp, keep]
  R_ends <- ends_at_peel$quantities$R[sp, keep]
  R_par  <- fit$quantities$R[sp, keep]

  # Marginalized: matches the shorter model to optimizer tolerance.
  testthat::expect_equal(as.numeric(R_peel), as.numeric(R_ends), tolerance = 1e-3)
  testthat::expect_equal(pk$quantities$ssb[sp, keep],
                         ends_at_peel$quantities$ssb[sp, keep], tolerance = 1e-3)

  # ...and the comparison is not vacuous: the peel must DIFFER from the parent,
  # which is the retrospective signal itself.
  testthat::expect_gt(max(abs(R_peel - R_par)) / max(abs(R_par)), 0.01)
})

testthat::test_that("a DSEM peel keeps its covariate and reports a finite Mohn's rho", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("dsem")

  d  <- dsem_retro_data()
  fc <- Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0)
  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 1,
    random_rec = TRUE, msmMode = 0,
    dsem = Rceattle::build_DSEM(sem = DSEM_RETRO_SEM, family = "fixed"),
    fit_control = fc)))

  retro <- suppressWarnings(suppressMessages(
    Rceattle::retrospective(fit, peels = 1, cores = 1, getsd = FALSE)))
  pk <- retro$Rceattle_list[[1]]

  # Under family = "fixed" the covariate column of x_tj IS the environmental
  # data, held fixed by the map. Zeroing x_tj to "peel" it would delete the
  # covariate outright rather than withhold it.
  xt <- as.matrix(pk$estimated_params$dsem_x_tj)
  testthat::expect_equal(ncol(xt), 4L)              # 3 recdev columns + temp
  testthat::expect_gt(max(abs(xt[, 4])), 1e-6)

  rho <- retro$mohns
  ssb_rho <- rho[rho$Object == "ssb", -(1:3), drop = FALSE]
  testthat::expect_true(any(is.finite(as.matrix(ssb_rho))))
})


# Comparing projection methods on ONE DSEM fit.
#
# This is the comparison hindcast_skill() exists for: does the SEM's correlation
# structure forecast the peeled years better than mean recruitment? Both answers
# come from the same fit, because `forecast_rec` only changes the FORECAST refit
# -- the peeled hindcast is identical either way -- so the projection rule is
# isolated from everything else.
#
# What makes this worth a test is that the comparison is silently vacuous
# whenever proj_mean_rec = TRUE: it takes precedence over the model's own
# process, build_srr() defaults it to TRUE, and then both settings return the
# same forecast. Nothing in the suite asserted the two settings differ at all, so
# a regression that quietly made "model" behave like "mean" would have left every
# DSEM projection comparison reading "no difference" and passed.
#
# The sem is LAGGED and carries a covariate (see DSEM_RETRO_SEM), so there is a
# correlation structure for the forecast to propagate. An IID sem would show much
# less and could not detect the failure.
testthat::test_that("forecast_rec tells the DSEM's projection from mean recruitment", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("dsem")

  d  <- dsem_retro_data()
  fc <- Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0)

  # proj_mean_rec = FALSE is the whole point: with TRUE the two settings agree
  # by construction and the assertions below would be measuring nothing.
  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 1,
    random_rec = TRUE, msmMode = 0,
    dsem = Rceattle::build_DSEM(sem = DSEM_RETRO_SEM, family = "fixed"),
    recFun = Rceattle::build_srr(srr_fun = 0, proj_mean_rec = FALSE),
    fit_control = fc)))
  testthat::expect_false(isTRUE(as.logical(fit$data_list$proj_mean_rec)))

  endyr <- fit$data_list$endyr
  styr  <- fit$data_list$styr
  depth <- 3L
  nm    <- paste0("Year_", endyr - depth)

  # A repeated depth asks for that one peel and nothing shallower, which keeps
  # this to two peel fits per setting rather than six.
  retro <- function(frec) suppressWarnings(suppressMessages(
    Rceattle::retrospective(fit, peels = c(depth, depth), cores = 1,
                            getsd = FALSE, forecast_rec = frec)))
  r_mean  <- retro("mean")
  r_model <- retro("model")

  testthat::skip_if(is.null(r_mean$Rceattle_list[[nm]]) ||
                    is.null(r_model$Rceattle_list[[nm]]),
                    "the 3-year DSEM peel did not converge")
  p_mean  <- r_mean$Rceattle_list[[nm]]
  p_model <- r_model$Rceattle_list[[nm]]

  hind <- seq_len(endyr - depth - styr + 1L)       # styr:endyr_peel
  fore <- (endyr - depth - styr + 2L):(endyr - styr + 1L)
  spp1 <- fit$data_list$spnames[1]

  # The peeled hindcast is the same fit either way: `forecast_rec` is read only
  # when the forecast refit's starting recruitment is written. If this ever
  # fails, the setting is reaching the hindcast and the comparison is no longer
  # isolating the projection rule.
  testthat::expect_equal(p_mean$quantities$R[, hind],
                         p_model$quantities$R[, hind])

  # And the forecast years genuinely differ: the SEM's lagged and covariate
  # paths carry the terminal state forward, where "mean" flattens it to the
  # hindcast average. Species 1 is the one the sem gives a lag and a covariate.
  R_mean  <- p_mean$quantities$R[1, fore]
  R_model <- p_model$quantities$R[1, fore]
  testthat::expect_gt(max(abs(R_model - R_mean)) / max(abs(R_mean)), 0.01)

  # "mean" really is flat across the forecast years, which is what the DSEM is
  # being compared against.
  testthat::expect_lt(stats::sd(R_mean) / mean(R_mean), 1e-6)

  # Both are scoreable, so the comparison can actually be read off a MASE.
  #
  # Scored on RECRUITMENT, not SSB. Over a 3-year horizon the recruits the two
  # rules disagree about are not mature yet, so SSB barely moves: measured on
  # this fixture, a 46.8% difference in forecast recruitment moved forecast SSB
  # 0.03%, and the SSB MASE agreed to 4 significant figures for two of the three
  # species. An SSB-scored assertion would therefore pass or fail on rounding
  # rather than on the projection rule.
  skill_mean  <- Rceattle::hindcast_skill(fit, retro = r_mean,  quantity = "R")
  skill_model <- Rceattle::hindcast_skill(fit, retro = r_model, quantity = "R")
  testthat::expect_true(all(is.finite(skill_mean$mase$mae_forecast)))
  testthat::expect_true(all(is.finite(skill_model$mase$mae_forecast)))

  # Species 1 is the one the sem gives a lag and a covariate, so it is the row
  # the two rules must disagree on.
  # Asserted at the DEEPEST horizon, against the separation actually measured
  # there -- not at horizon 1, and not as a max over horizons.
  #
  # A lag-1 SEM's one-step forecast is rho times the terminal deviation, so it
  # leaves from near the mean and separates as it decays. Measured on this
  # fixture, on the statistic ACTUALLY asserted here -- the relative difference
  # in mae_forecast for species 1 -- that is 0.34% at h = 1, 46.2% at h = 2 and
  # 96.7% at h = 3. (The often-quoted 0.2% / 46.8% pair is the difference in
  # forecast RECRUITMENT, which is not what this line tests.) Horizon 1 asks the
  # question where the answer is smallest, so assert at the deepest.
  #
  # The threshold is 10%, an order of magnitude under the 96.7% measured. Species
  # 1 is the SEM-linked one; arrowtooth measures 8.5% at h = 3, UNDER the
  # threshold, so this depends on spnames[1] remaining the linked species.
  pick <- function(sk) {
    z <- sk$mase[sk$mase$species == spp1, ]
    stats::setNames(z$mae_forecast, z$years_ahead)
  }
  mf_mean  <- pick(skill_mean)
  mf_model <- pick(skill_model)
  h <- intersect(names(mf_mean), names(mf_model))
  # Guarded: with no shared horizon, max(integer(0)) is -Inf and the lookup
  # below throws "subscript out of bounds" instead of failing as an assertion.
  testthat::expect_gt(length(h), 0L)
  testthat::skip_if(length(h) == 0L, "no shared horizon to compare at")
  deepest <- as.character(max(as.integer(h)))
  testthat::expect_gt(
    abs(mf_model[[deepest]] - mf_mean[[deepest]]) / mf_mean[[deepest]], 0.10)
})


# The counterpart: with proj_mean_rec = TRUE the comparison IS vacuous, and
# retrospective() has to say so. The warning used to be raised inside the
# per-peel closure under `i == 1L`, which meant it was discarded by a PSOCK
# worker under the default `cores` and skipped entirely by `peels = 2:10`. It is
# now settled once before dispatch, so a subset request still gets it.
testthat::test_that("retrospective warns that forecast_rec is inert, whatever peels are asked for", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("dsem")

  d  <- dsem_retro_data()
  fc <- Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0)
  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 1,
    random_rec = TRUE, msmMode = 0,
    dsem = Rceattle::build_DSEM(sem = DSEM_RETRO_SEM, family = "fixed"),
    recFun = Rceattle::build_srr(srr_fun = 0, proj_mean_rec = TRUE),
    fit_control = fc)))
  testthat::expect_true(isTRUE(as.logical(fit$data_list$proj_mean_rec)))

  # Raised before any peel is fitted, so a depth set that never includes 1 still
  # gets it. Only the warning is under test, so the peels need not converge.
  testthat::expect_warning(
    suppressMessages(try(Rceattle::retrospective(
      fit, peels = 2:3, cores = 1, getsd = FALSE, forecast_rec = "model"),
      silent = TRUE)),
    "inert on this fit")

  # Not raised when the caller did not ask for the model's own process. Collect
  # the warnings and look for this one rather than asserting none at all: a peel
  # drop or a warm-start note is unrelated and must not fail this.
  warned <- character(0)
  withCallingHandlers(
    suppressMessages(try(Rceattle::retrospective(
      fit, peels = c(3, 3), cores = 1, getsd = FALSE, forecast_rec = "mean"),
      silent = TRUE)),
    warning = function(cnd) {
      warned <<- c(warned, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    })
  testthat::expect_false(any(grepl("inert on this fit", warned)))
})
