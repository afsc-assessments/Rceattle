# M1_model = 3 (sex- and age-specific M1) with M1_re 2 or 5 (deviations that
# vary by YEAR and are constant over ages) is now available. 5.53.0 refused it
# along with the rest of M1_model 3/4/5, because no arm wrote its deviations.
#
# Why this one and not the others, which stay refused:
#
#   * M1_re 1/4 vary by AGE. An age-varying deviation and an age-specific level
#     enter the likelihood only as their sum, so the mode puts the deviations
#     at zero and absorbs everything into log_M1 -- not identifiable. The
#     identifiable version of that model already exists as M1_model = 1 with
#     M1_re = 4.
#   * M1_re 3/6 vary by age AND year, so they carry the same age-dimension
#     collinearity.
#   * M1_model 4/5 are in .M1_DEPRECATED_MODELS, and the linkage grammar says
#     the same thing as linkage_spec(~ temp + (1 | Year)).
#
# A by-year deviation that is constant over ages is orthogonal to an
# age-specific time-constant level apart from the mean, which the zero-mean
# random effect assigns to log_M1. That is WHAM's age-specific-mean M with
# ar1_y deviations.
#
# Shared across sexes, like the other levels in that family: the template
# scores num_re_sexes fields, which is 1 unless M1_model == 2
# (src/TMB/ceattle.cpp), so one field is what M1_model 3 owes.
#
# Measured on GOA2018SS: log_M1 keeps its 64 age-specific parameters, 126
# deviations are freed (42 hindcast years x 3 species), 3 sds, and the
# deviations are a Laplace random block of 126 rather than penalised fixed
# effects. Identifiability here is the analytic argument above plus non-zero
# gradients; a simulation-recovery study would be stronger and is not done.

fit_m1 <- function(mm, re) {
  data("GOA2018SS", package = "Rceattle", envir = environment())
  suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = GOA2018SS, inits = NULL, file = NULL,
    estimateMode = "DebugBuild", msmMode = 0, niter = 3, random_rec = FALSE,
    M1Fun = Rceattle::build_M1(M1_model = mm, M1_re = re),
    fit_control = Rceattle::fit_control(verbose = 0))))
}


testthat::test_that("M1_model 3 estimates by-year M deviations", {
  testthat::skip_on_cran()
  want_rho <- c(`2` = 0L, `5` = 6L)   # 0 for IID, 6 for the year AR1
  for (re in c(2, 5)) {
    fit <- fit_m1(3, re)
    mp <- fit$obj$env$map
    lab <- paste("M1_re", re)

    # The age-specific level survives alongside the deviations.
    testthat::expect_equal(sum(!is.na(mp$log_M1)), 64L, info = lab)
    # 42 hindcast years x 3 species, shared across sexes.
    testthat::expect_equal(
      length(unique(stats::na.omit(as.integer(mp$log_M1_dev)))), 126L,
      info = lab)
    testthat::expect_equal(
      length(unique(stats::na.omit(as.integer(mp$M1_dev_log_sd)))), 3L,
      info = lab)
    testthat::expect_equal(sum(!is.na(mp$M1_rho)), want_rho[[as.character(re)]],
                           info = lab)
    # Integrated, not penalised: the deviations are the Laplace block.
    testthat::expect_equal(length(fit$obj$env$random), 126L, info = lab)
  }
})


testthat::test_that("the by-year deviations are informed by the data", {
  # Freeing a parameter the data cannot move would be worse than refusing it:
  # a zero-gradient row gives a singular Hessian, and newtonsteps solves with
  # it. Every deviation carries gradient here, and moving one moves the
  # objective.
  testthat::skip_on_cran()
  fit <- fit_m1(3, 2)
  p <- fit$obj$env$last.par.best
  g <- fit$obj$gr(p)
  dev <- which(names(p) == "log_M1_dev")

  testthat::expect_true(all(is.finite(g)))
  testthat::expect_false(any(g[dev] == 0))
  q <- p
  q[dev[1]] <- 0.2
  testthat::expect_false(isTRUE(all.equal(fit$obj$fn(p), fit$obj$fn(q))))
})


testthat::test_that("the by-age M1_model 3 families stay refused", {
  # The non-identifiable half of the decision. If a future change opens one of
  # these, it should have to delete this block and say why.
  testthat::skip_on_cran()
  for (re in c(1, 3, 4, 6)) {
    testthat::expect_error(fit_m1(3, re),
      "is not implemented for M1_model",
      info = paste("M1_re", re))
  }
  # And M1_model 4/5 stay refused for every family.
  for (mm in c(4, 5)) {
    for (re in c(2, 5)) {
      testthat::expect_error(fit_m1(mm, re),
        "is not implemented for M1_model",
        info = paste("M1_model", mm, "M1_re", re))
    }
  }
})
