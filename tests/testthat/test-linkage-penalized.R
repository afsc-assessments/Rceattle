# =============================================================================
# Penalized fixed-effect linkage deviations: linkage_spec(integrate = FALSE).
#
# A random-effect term's deviations are normally integrated out by the Laplace
# approximation. `integrate = FALSE` instead estimates them as a plain penalized
# fixed effect, which is what the ADMB/AMAK reference models (and the legacy
# Time_varying_sel / Time_varying_q switches) do. It is permitted only when the
# deviation SD is fixed: a variance is consistently estimated only by
# integrating, so estimating the deviations AND their SD jointly as fixed
# effects is degenerate -- the likelihood improves without bound as both go to
# zero. The guards below are the constructor-level half of that contract; the
# pooled-group half lives in .re_group_table().
# =============================================================================

testthat::test_that("integrate must be a single TRUE/FALSE", {
  for (bad in list("yes", c(TRUE, FALSE), NA, 1)) {
    testthat::expect_error(
      Rceattle::linkage_spec(~ rw(1 | Year), init = list(sigma = 0.05),
                             integrate = bad),
      "single TRUE/FALSE")
  }
})

testthat::test_that("integrate = FALSE requires a fixed deviation SD", {
  # No sigma supplied: the SD would be estimated alongside the deviations.
  testthat::expect_error(
    Rceattle::linkage_spec(~ rw(1 | Year), integrate = FALSE),
    "requires a fixed deviation SD")

  # A sigma PRIOR still estimates the SD -- shrunk, but not fixed -- so it does
  # not satisfy the contract even though `init` is present.
  testthat::expect_error(
    Rceattle::linkage_spec(~ rw(1 | Year), init = list(sigma = 0.05),
                           priors = list(sigma = lognormal(log(0.05), 0.3)),
                           integrate = FALSE),
    "requires a fixed deviation SD")
})

testthat::test_that("integrate = FALSE on ar1 also requires a fixed rho", {
  testthat::expect_error(
    Rceattle::linkage_spec(~ ar1(1 | Year), init = list(sigma = 0.05),
                           integrate = FALSE),
    "requires a fixed correlation")
  testthat::expect_error(
    Rceattle::linkage_spec(~ ar1(1 | Year), init = list(sigma = 0.05, rho = 0.5),
                           priors = list(rho = normal(0, 0.3)),
                           integrate = FALSE),
    "requires a fixed correlation")
  # Both pinned: accepted.
  spec <- Rceattle::linkage_spec(~ ar1(1 | Year),
                                 init = list(sigma = 0.05, rho = 0.5),
                                 integrate = FALSE)
  testthat::expect_false(spec$integrate)
})

testthat::test_that("est_phase cannot silently fail to fix a random-effect term", {
  # est_phase reaches only beta_linkage; a random effect's deviations live in
  # beta_linkage_re / beta_linkage_re_pen, so est_phase = 0 used to leave every
  # deviation estimated -- the opposite of what the caller asked for.
  for (spec_fn in list(
    function() Rceattle::linkage_spec(~ rw(1 | Year), init = list(sigma = 0.05),
                                      integrate = FALSE, est_phase = 0L),
    function() Rceattle::linkage_spec(~ (1 | Year), est_phase = 0L),
    function() Rceattle::linkage_spec(~ temp + rw(1 | Year), est_phase = 0L))) {
    testthat::expect_error(spec_fn(), "cannot fix a random-effect term")
  }

  # A purely fixed-effect linkage still honours est_phase = 0, and a phase of 1
  # (or above) is accepted on a random-effect term.
  testthat::expect_equal(
    Rceattle::linkage_spec(~ temp, est_phase = 0L)$est_phase, 0L)
  testthat::expect_equal(
    Rceattle::linkage_spec(~ rw(1 | Year), est_phase = 2L)$est_phase, 2L)
})

testthat::test_that("integrate = FALSE is rejected without a random-effect term", {
  testthat::expect_error(
    Rceattle::linkage_spec(~ temp, init = list(sigma = 0.05), integrate = FALSE),
    "random-effect term")
})

testthat::test_that("integrate = FALSE is rejected with observe", {
  # An observed latent state is a state-space model: the state must be
  # integrated for the observation to identify it.
  testthat::expect_error(
    Rceattle::linkage_spec(~ ar1(1 | Year), observe = "qcov", obs_sd = 0.2,
                           init = list(sigma = 0.05, rho = 0.5),
                           integrate = FALSE),
    "cannot be combined with `observe`")
})

testthat::test_that("integrate defaults to TRUE and round-trips onto the spec", {
  testthat::expect_true(Rceattle::linkage_spec(~ rw(1 | Year))$integrate)
  testthat::expect_false(
    Rceattle::linkage_spec(~ rw(1 | Year), init = list(sigma = 0.05),
                           integrate = FALSE)$integrate)
})

testthat::test_that("the linkage table carries re_integrate on RE rows only", {
  env <- data.frame(Year = 2000:2009)
  pool <- function(spec) {
    Rceattle:::pool_linkages(list(q = list(q = spec)), env_data = env,
                             strata = list(fleet = 1L))$table
  }

  pen <- pool(Rceattle::linkage_spec(~ rw(1 | Year), by = ~ fleet, fleet = 1L,
                                     init = list(sigma = 0.05),
                                     integrate = FALSE))
  is_re <- !is.na(pen$re_struct)
  testthat::expect_true(any(is_re))
  testthat::expect_true(all(pen$re_integrate[is_re] == FALSE))
  # Fixed rows (the intercept that carries the walk's level) stay NA.
  testthat::expect_true(all(is.na(pen$re_integrate[!is_re])))

  int <- pool(Rceattle::linkage_spec(~ rw(1 | Year), by = ~ fleet, fleet = 1L))
  testthat::expect_true(all(int$re_integrate[!is.na(int$re_struct)] == TRUE))
})

testthat::test_that("the group registry is unchanged for integrated models", {
  # Back-compat: adding the flag must not renumber a single existing sigma group.
  env <- data.frame(Year = 2000:2004)
  tbl <- Rceattle:::pool_linkages(
    list(q = list(q = Rceattle::linkage_spec(~ rw(1 | Year), by = ~ fleet,
                                             fleet = 1L))),
    env_data = env, strata = list(fleet = 1L))$table
  re <- !is.na(tbl$re_struct)
  testthat::expect_equal(unique(tbl$sigma_index[re]), 0L)
  testthat::expect_setequal(tbl$re_index[re], seq_len(sum(re)) - 1L)

  gt <- Rceattle:::.re_group_table(tbl)
  testthat::expect_true(gt$integrate)
  testthat::expect_false(gt$sigma_fixed)   # estimated SD is fine when integrated
})

testthat::test_that("a penalized group carries integrate = FALSE and a fixed SD", {
  env <- data.frame(Year = 2000:2004)
  tbl <- Rceattle:::pool_linkages(
    list(q = list(q = Rceattle::linkage_spec(~ rw(1 | Year), by = ~ fleet,
                                             fleet = 1L,
                                             init = list(sigma = 0.05),
                                             integrate = FALSE))),
    env_data = env, strata = list(fleet = 1L))$table
  gt <- Rceattle:::.re_group_table(tbl)
  testthat::expect_false(gt$integrate)
  testthat::expect_true(gt$sigma_fixed)
  testthat::expect_equal(gt$sigma_start, 0.05)
})

testthat::test_that("a penalized group with an estimated SD is refused at pooling", {
  # The constructor cannot catch this: each spec is individually legal, and only
  # the pooled group reveals that the SD is estimated for a penalized effect.
  # Build the table by hand to reach .re_group_table() with that combination.
  env <- data.frame(Year = 2000:2004)
  tbl <- Rceattle:::pool_linkages(
    list(q = list(q = Rceattle::linkage_spec(~ rw(1 | Year), by = ~ fleet,
                                             fleet = 1L,
                                             init = list(sigma = 0.05),
                                             integrate = FALSE))),
    env_data = env, strata = list(fleet = 1L))$table
  tbl$re_sigma_init <- NA_real_          # drop the fix -> SD becomes estimated
  testthat::expect_error(Rceattle:::.re_group_table(tbl),
                         "require a fixed deviation SD")
})

testthat::test_that("a group cannot be both integrated and penalized", {
  env <- data.frame(Year = 2000:2004)
  tbl <- Rceattle:::pool_linkages(
    list(q = list(q = Rceattle::linkage_spec(~ rw(1 | Year), by = ~ fleet,
                                             fleet = 1L,
                                             init = list(sigma = 0.05),
                                             integrate = FALSE))),
    env_data = env, strata = list(fleet = 1L))$table
  re <- which(!is.na(tbl$re_struct))
  tbl$re_integrate[re[1]] <- TRUE        # straddle one row across both vectors
  testthat::expect_error(Rceattle:::.assign_re_registry(tbl),
                         "both integrated and penalized")
})

# ---- end-to-end: the penalized path reproduces the legacy penalized fit ------

# Both models are built at estimateMode = 3 and evaluated at the SAME deviations,
# so the comparison is exact and free of optimizer noise.
.pen_build <- function(dat, qfun = NULL) {
  suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = dat, file = NULL, inits = NULL, estimateMode = 3,
    random_rec = FALSE, msmMode = 0, qFun = qfun,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0))))
}

testthat::test_that("a penalized rw() reproduces legacy random-walk catchability", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # Legacy Time_varying_q = "RandomWalk" IS a penalized fixed effect with a fixed
  # SD: index_q_dev_log_sd is estimable only for Catchability == "AR1", the first
  # deviate is pinned, the density is a first difference, and the deviate enters
  # additively on the log scale exactly where q_linkage_offset does. So the two
  # formulations are the same model and must agree exactly -- this is the check
  # that the penalized density is neither double-counted nor dropped.
  FLT   <- 7L        # EIT_Pollock, the only Estimated-q survey in BS2017SS
  SIGMA <- 0.05
  d0    <- Rceattle::BS2017SS
  nyr   <- length(d0$styr:d0$endyr)

  d_leg <- d0
  d_leg$fleet_control$Time_varying_q[FLT]    <- "RandomWalk"
  d_leg$fleet_control$Time_varying_q_sd[FLT] <- SIGMA
  f_leg <- .pen_build(d_leg)

  d_gr <- d0
  d_gr$fleet_control$Time_varying_q[FLT] <- "Off"
  f_gr <- .pen_build(d_gr, Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ rw(1 | Year), by = ~ fleet, fleet = FLT,
                               link = "log", init = list(sigma = SIGMA),
                               integrate = FALSE))))

  # The penalized deviations are ordinary fixed effects, not random effects.
  testthat::expect_true("beta_linkage_re_pen" %in% names(f_gr$obj$par))
  testthat::expect_false("beta_linkage_re" %in% names(f_gr$obj$par))
  testthat::expect_length(f_gr$obj$env$random, 0)
  testthat::expect_equal(length(f_gr$obj$par), length(f_leg$obj$par))

  set.seed(7)
  dev <- c(0, stats::rnorm(nyr - 1, 0, SIGMA))   # both pin the first deviate

  p_leg <- f_leg$obj$par
  p_gr  <- f_gr$obj$par
  i_leg <- names(p_leg) == "index_q_dev"
  i_gr  <- names(p_gr)  == "beta_linkage_re_pen"
  testthat::expect_equal(sum(i_leg), nyr - 1L)
  testthat::expect_equal(sum(i_gr),  nyr - 1L)
  p_leg[i_leg] <- dev[-1]
  p_gr[i_gr]   <- dev[-1]

  testthat::expect_equal(f_gr$obj$fn(p_gr), f_leg$obj$fn(p_leg), tolerance = 1e-12)

  # ... and the density lands in the linkage row at the value written by hand.
  expected <- -sum(stats::dnorm(diff(dev), 0, SIGMA, log = TRUE))
  r_gr <- f_gr$obj$report(p_gr)
  testthat::expect_equal(sum(r_gr$jnll_comp[linkage_re_row(), ]), expected,
                         tolerance = 1e-10)
})

testthat::test_that("integrating or penalizing changes the objective, not the density", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # The single clearest statement of what the feature is: both treatments carry
  # the SAME row-20 density; they differ only in whether TMB integrates it out.
  FLT   <- 7L
  SIGMA <- 0.05
  d     <- Rceattle::BS2017SS
  spec  <- function(integrate) Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ rw(1 | Year), by = ~ fleet, fleet = FLT,
                               link = "log", init = list(sigma = SIGMA),
                               integrate = integrate)))
  f_int <- .pen_build(d, spec(TRUE))
  f_pen <- .pen_build(d, spec(FALSE))

  testthat::expect_true("beta_linkage_re" %in% names(f_int$obj$env$par[f_int$obj$env$random]))
  testthat::expect_length(f_pen$obj$env$random, 0)

  nyr <- length(d$styr:d$endyr)
  set.seed(11)
  dev <- c(0, stats::rnorm(nyr - 1, 0, SIGMA))

  # last.par omits the map-pinned first deviate, which stays at its init of 0 --
  # so supply the remaining nyr - 1 and let the pin carry dev[1] = 0.
  fill <- function(obj, nm) {
    p <- obj$env$last.par
    testthat::expect_equal(sum(names(p) == nm), length(dev) - 1L)
    p[names(p) == nm] <- dev[-1]
    p
  }
  r_int <- f_int$obj$report(fill(f_int$obj, "beta_linkage_re"))
  r_pen <- f_pen$obj$report(fill(f_pen$obj, "beta_linkage_re_pen"))

  row <- linkage_re_row()
  testthat::expect_equal(sum(r_int$jnll_comp[row, ]),
                         sum(r_pen$jnll_comp[row, ]), tolerance = 1e-12)
  testthat::expect_equal(sum(r_pen$jnll_comp[row, ]),
                         -sum(stats::dnorm(diff(dev), 0, SIGMA, log = TRUE)),
                         tolerance = 1e-10)
})

testthat::test_that("one model can mix an integrated and a penalized group", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # The motivating requirement: TMB selects random effects by parameter NAME, so
  # mixing the two treatments in one fit is exactly what the second parameter
  # vector exists for.
  d <- Rceattle::BS2017SS
  SIG_Q <- 0.05
  SIG_S <- 0.10
  qfun <- Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ rw(1 | Year), by = ~ fleet, fleet = 7L,
                               link = "log", init = list(sigma = SIG_Q),
                               integrate = FALSE)))
  selfun <- Rceattle::build_selectivity(linkages = list(
    inf_asc = Rceattle::linkage_spec(~ rw(1 | Year), by = ~ fleet, fleet = 4L,
                                     init = list(sigma = SIG_S))))
  f <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = d, file = NULL, inits = NULL, estimateMode = 3,
    random_rec = FALSE, msmMode = 0, qFun = qfun, selFun = selfun,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0))))

  n_int <- length(f$estimated_params$beta_linkage_re)
  n_pen <- length(f$estimated_params$beta_linkage_re_pen)
  testthat::expect_gt(n_int, 0)
  testthat::expect_gt(n_pen, 0)
  # Every slot is stored exactly once, across the two vectors.
  testthat::expect_equal(n_int + n_pen,
                         sum(!is.na(f$data_list$linkage_table$re_index)))

  # Only the integrated vector is in the random set; only the penalized one is a
  # free fixed effect.
  testthat::expect_true("beta_linkage_re" %in%
                        names(f$obj$env$par[f$obj$env$random]))
  testthat::expect_false("beta_linkage_re_pen" %in%
                         names(f$obj$env$par[f$obj$env$random]))
  testthat::expect_true("beta_linkage_re_pen" %in% names(f$obj$par))

  # Row 20 carries both groups' densities: fill both vectors and check the total
  # against the two hand-written random-walk densities.
  set.seed(3)
  p <- f$obj$env$last.par
  d_int <- c(0, stats::rnorm(n_int - 1, 0, SIG_S))
  d_pen <- c(0, stats::rnorm(n_pen - 1, 0, SIG_Q))
  # Both walks pin their first deviate, so last.par holds n - 1 of each.
  p[names(p) == "beta_linkage_re"]     <- d_int[-1]
  p[names(p) == "beta_linkage_re_pen"] <- d_pen[-1]
  rep <- f$obj$report(p)
  testthat::expect_equal(
    sum(rep$jnll_comp[linkage_re_row(), ]),
    -sum(stats::dnorm(diff(d_int), 0, SIG_S, log = TRUE)) -
      sum(stats::dnorm(diff(d_pen), 0, SIG_Q, log = TRUE)),
    tolerance = 1e-8)
})

testthat::test_that("integrate survives a save_config round trip", {
  spec <- Rceattle::linkage_spec(~ rw(1 | Year), by = ~ fleet, fleet = 1L,
                                 init = list(sigma = 0.05), integrate = FALSE)
  lst <- Rceattle:::.rce_spec_to_list(spec)
  testthat::expect_false(lst$integrate)
  testthat::expect_false(Rceattle:::.rce_spec_from_list(lst)$integrate)
})


testthat::test_that("re_pos locates a deviation in the vector that holds it", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  # re_index is the GLOBAL slot across both parameter vectors; re_pos is the
  # position within the one that actually holds the deviation. In a model with
  # only one treatment they coincide, which is exactly why indexing by re_index
  # looks correct until a mixed model appears -- and then silently writes to the
  # wrong deviation. Pin the divergence so the two cannot be conflated again.
  d <- Rceattle::BS2017SS
  f <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = d, file = NULL, inits = NULL, estimateMode = 3,
    random_rec = FALSE, msmMode = 0,
    qFun = Rceattle::build_catchability(linkages = list(
      q = Rceattle::linkage_spec(~ rw(1 | Year), by = ~ fleet, fleet = 7L,
                                 link = "log"))),                     # integrated
    selFun = Rceattle::build_selectivity(linkages = list(
      inf_asc = Rceattle::linkage_spec(~ rw(1 | Year), by = ~ fleet, fleet = 4L,
                                       init = list(sigma = 0.1),
                                       integrate = FALSE))),          # penalized
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE,
                                        verbose = 0))))

  tb <- as.data.frame(f$data_list$linkage_table)
  tb <- tb[!is.na(tb$re_index), ]
  int <- tb[tb$re_integrate, ]
  pen <- tb[!tb$re_integrate, ]
  testthat::expect_gt(nrow(int), 0)
  testthat::expect_gt(nrow(pen), 0)

  # Each vector is indexed densely from 0 by re_pos ...
  testthat::expect_setequal(int$re_pos, seq_len(nrow(int)) - 1L)
  testthat::expect_setequal(pen$re_pos, seq_len(nrow(pen)) - 1L)
  testthat::expect_equal(max(int$re_pos) + 1L,
                         length(f$estimated_params$beta_linkage_re))
  testthat::expect_equal(max(pen$re_pos) + 1L,
                         length(f$estimated_params$beta_linkage_re_pen))

  # ... and for the penalized group that is NOT the global slot, which is the
  # whole point: re_index would run past the end of beta_linkage_re_pen.
  testthat::expect_false(identical(pen$re_pos, pen$re_index))
  testthat::expect_gt(max(pen$re_index) + 1L,
                      length(f$estimated_params$beta_linkage_re_pen))
})


testthat::test_that("printing a spec shows link and the penalized treatment", {
  # Both change the model: `link` decides whether an offset is added to the
  # target or multiplies it, and `integrate` whether the deviations are
  # integrated out or penalized. Neither should have to be read off the call.
  pen <- Rceattle::linkage_spec(~ rw(1 | Year), init = list(sigma = 0.05),
                                integrate = FALSE)
  testthat::expect_output(print(pen), "link:")
  testthat::expect_output(print(pen), "penalized fixed effect")

  int <- Rceattle::linkage_spec(~ rw(1 | Year))
  testthat::expect_output(print(int), "link:")
  out <- utils::capture.output(print(int))
  testthat::expect_false(any(grepl("penalized", out)))
})
