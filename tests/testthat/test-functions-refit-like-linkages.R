# =============================================================================
# A refit reproduces the source model's linkage structure.
#
# fit_mod() takes qFun / selFun / compFun as the source of truth for
# data_list$q_linkages / sel_linkages / comp_linkages, so .refit_like() rebuilds
# all three from the source model. That is what makes a refit the same model:
# the linkage table sizes the beta_linkage parameter block, which has to stay
# matched to the parameters carried over in `inits`.
#
# Every diagnostic that re-fits -- retrospective(), jitter(), self_test(),
# profile(), run_mse(), remove_F(), sample_rec() -- routes through
# .refit_like(), so this covers all of them.
# =============================================================================

.refit_linkage_env_data <- function(d) {
  yrs <- d$styr:d$projyr
  data.frame(Year = yrs, temp = as.numeric(scale(sin(seq_along(yrs) / 4))))
}

# A model carrying one of each linkage type routed through qFun / selFun /
# compFun: a catchability covariate, a selectivity covariate, and a prior-only
# composition (Dirichlet-multinomial) weight.
.refit_linkage_fit <- function() {
  d <- Rceattle::BS2017SS
  d$env_data <- .refit_linkage_env_data(d)

  # Fleet 7 (EIT_Pollock) is the only Estimated-q survey.
  qfun <- Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ temp, by = ~ fleet, fleet = 7L)))

  fc <- Rceattle::switch_check(Rceattle::clean_data(d))$fleet_control
  logistic_flt <- fc$Fleet_code[which(fc$Selectivity == "Logistic")[1]]
  selfun <- Rceattle::build_selectivity(linkages = list(
    inf_asc = Rceattle::linkage_spec(~ temp, by = ~ fleet,
                                     fleet = logistic_flt)))

  comp_flts <- sort(unique(d$comp_data$Fleet_code[d$comp_data$Fleet_code > 0]))
  d$fleet_control$Comp_distribution[comp_flts] <- "DirichletMultinomial"
  compfun <- Rceattle::build_composition(linkages = list(
    theta_comp = Rceattle::linkage_spec(
      ~ 1, by = ~ fleet, fleet = comp_flts,
      priors = list(`(Intercept)` = Rceattle::prior_lognormal(0, 2)))))

  suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, qFun = qfun, selFun = selfun, compFun = compfun,
    estimateMode = 3, msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0))))
}


testthat::test_that(".refit_like() preserves q / sel / comp linkages", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  fit <- .refit_linkage_fit()
  # The source model really does carry all three (guards the fixture itself).
  testthat::expect_false(is.null(fit$data_list$q_linkages))
  testthat::expect_false(is.null(fit$data_list$sel_linkages))
  testthat::expect_false(is.null(fit$data_list$comp_linkages))

  refit <- suppressWarnings(suppressMessages(Rceattle:::.refit_like(
    data_list = fit$data_list, inits = fit$estimated_params,
    estimateMode = 3)))

  for (nm in c("q_linkages", "sel_linkages", "comp_linkages")) {
    testthat::expect_equal(refit$data_list[[nm]], fit$data_list[[nm]],
                           info = nm)
  }
})


testthat::test_that("a linkage-carrying model refits to the same objective", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  fit <- .refit_linkage_fit()
  refit <- suppressWarnings(suppressMessages(Rceattle:::.refit_like(
    data_list = fit$data_list, inits = fit$estimated_params,
    estimateMode = 3)))

  # The linkage table sizes beta_linkage, so the parameter shapes are the
  # sharpest signal that the structure survived the refit.
  testthat::expect_identical(lengths(refit$estimated_params),
                             lengths(fit$estimated_params))
  testthat::expect_gt(length(refit$estimated_params$beta_linkage), 0)
  testthat::expect_equal(nrow(refit$data_list$linkage_table),
                         nrow(fit$data_list$linkage_table))

  # Same structure, same starting parameters => same objective.
  testthat::expect_equal(refit$obj$fn(refit$obj$par),
                         fit$obj$fn(fit$obj$par),
                         tolerance = 1e-10)
})


testthat::test_that("sample_rec(update_model = TRUE) keeps the model's linkages", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # sample_rec() rebuilds the model to propagate the new projection
  # recruitment, and routes that rebuild through .refit_like().
  fit <- .refit_linkage_fit()
  # sample_rec() rebuilds a fit that ran no projection at estimateMode = 3, then
  # restores the source value; label the source differently so that restore is
  # observable rather than trivially equal to the build-only mode.
  fit$data_list$estimateMode <- 1L

  resampled <- suppressWarnings(suppressMessages(
    Rceattle::sample_rec(fit, sample_rec = FALSE, update_model = TRUE)))

  for (nm in c("q_linkages", "sel_linkages", "comp_linkages")) {
    # Guard the fixture too: without this, a model that silently lost its
    # linkages would compare NULL to NULL and pass while covering nothing.
    testthat::expect_false(is.null(fit$data_list[[nm]]), info = nm)
    testthat::expect_equal(resampled$data_list[[nm]], fit$data_list[[nm]],
                           info = nm)
  }
  testthat::expect_identical(lengths(resampled$estimated_params),
                             lengths(fit$estimated_params))
  testthat::expect_gt(length(resampled$estimated_params$beta_linkage), 0)
  testthat::expect_equal(resampled$data_list$estimateMode, 1L)
})


testthat::test_that("a refit preserves a penalized (integrate = FALSE) linkage", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # .refit_like() rebuilds qFun/selFun/compFun from the stored spec objects, so
  # `integrate` rides along as a spec field. If it did not, a refit would quietly
  # switch the deviations from a penalized fixed effect to an integrated random
  # effect -- a different model -- and every diagnostic that re-fits
  # (retrospective, jitter, self_test, profile, run_mse, remove_F) would be
  # comparing against something the user never asked for.
  d <- Rceattle::BS2017SS
  qfun <- Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ rw(1 | Year), by = ~ fleet, fleet = 7L,
                               link = "log", init = list(sigma = 0.05),
                               integrate = FALSE)))
  fit <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = d, file = NULL, inits = NULL, estimateMode = 1,
    random_rec = FALSE, msmMode = 0, qFun = qfun,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0))))

  refit <- suppressWarnings(suppressMessages(Rceattle:::.refit_like(
    data_list = fit$data_list, inits = fit$estimated_params,
    estimateMode = 1)))

  # Still penalized: the flag survived, the deviations are still a fixed effect,
  # and nothing migrated into the Laplace approximation.
  testthat::expect_false(refit$data_list$q_linkages$q$integrate)
  testthat::expect_length(refit$obj$env$random, 0)
  testthat::expect_true("beta_linkage_re_pen" %in% names(refit$obj$par))
  testthat::expect_equal(
    length(refit$estimated_params$beta_linkage_re_pen),
    length(fit$estimated_params$beta_linkage_re_pen))

  # And it is the same model numerically. Compare the OPTIMIZED objectives, not
  # obj$fn(obj$par): obj$par holds the start values, and the refit starts from
  # the source model's MLEs while the source started from the defaults.
  testthat::expect_equal(refit$opt$objective, fit$opt$objective,
                         tolerance = 1e-8)
})
