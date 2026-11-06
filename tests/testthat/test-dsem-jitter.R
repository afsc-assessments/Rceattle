# jitter() on a DSEM.
#
# jitter() needs nothing special for a DSEM, unlike retrospective()'s peel: it
# refits under the stored map, which carries the DSEM blocks fit_mod() recorded,
# and its perturbation only touches entries whose map is not NA. Under a DSEM that lands exactly where it should -- the path coefficients
# and the latent recruitment-deviation columns are jittered, the covariate
# columns of dsem_x_tj are not, because with family = "fixed" those ARE the
# environmental data and perturbing them would jitter the DATA.
#
# The trap these tests are written against: "it ran and the fits converged" is
# NOT evidence the starting values moved. A jitter that silently skipped every
# DSEM parameter would return the same converged fits and the same objective.
# So assert on the STARTS, not only the ends.

testthat::test_that("jitter() perturbs the DSEM starts but never the covariate data", {
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

  # Which x_tj columns are the fixed covariate, and which the latent recdevs?
  mm <- matrix(as.numeric(fit$map$mapList$dsem_x_tj),
               nrow = nrow(fit$estimated_params$dsem_x_tj))
  fixed_cols <- which(apply(mm, 2, function(z) all(is.na(z))))
  testthat::expect_equal(length(fixed_cols), 1L)     # temp

  jit <- suppressWarnings(suppressMessages(
    Rceattle::jitter(fit, njitter = 2, sd = 0.1, seed = 7, cores = 1,
                     getsd = FALSE)))
  ml <- jit$Rceattle_list
  testthat::expect_gt(length(ml), 0L)

  parent_start <- fit$initial_params
  for (m in ml) {
    # THE assertion: the DSEM starting values actually moved. Without this the
    # test passes just as happily on a jitter that ignored every DSEM block.
    testthat::expect_gt(
      max(abs(as.numeric(m$initial_params$dsem_beta_z) -
                as.numeric(parent_start$dsem_beta_z))), 1e-8)

    xs <- as.matrix(m$initial_params$dsem_x_tj)
    xp <- as.matrix(parent_start$dsem_x_tj)
    # Latent recdev columns move...
    testthat::expect_gt(max(abs(xs[, -fixed_cols] - xp[, -fixed_cols])), 1e-8)
    # ...and the covariate does not: that is data, not a starting value.
    testthat::expect_equal(xs[, fixed_cols], xp[, fixed_cols], tolerance = 1e-12)
  }

  # A jitter is a convergence diagnostic, so the perturbed starts should return
  # to the same optimum rather than scatter.
  testthat::expect_true(all(is.finite(jit$nll)))
  testthat::expect_lt(max(abs(jit$nll - fit$opt$objective)) /
                        abs(fit$opt$objective), 1e-3)
})
