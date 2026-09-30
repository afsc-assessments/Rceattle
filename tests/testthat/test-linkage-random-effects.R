# IID random-effect linkage density: `~ (1 | group)` expands to one deviation
# per level, damped by a N(0, sigma) density in jnll_comp row 20 (0-based; the
# "Linkage random effects" row). These cover the registry/encoding wiring
# (fast) and the end-to-end fit + density (skip_on_cran).

# ---- registry + encoding (fast, no fit) -------------------------------------

testthat::test_that("an IID linkage builds a bijective RE registry", {
  env <- data.frame(Year = 2000:2009)
  spec <- Rceattle::linkage_spec(~ (1 | Year), param = "q", by = ~ fleet,
                                 fleet = 1L, link = "log")
  pool <- Rceattle:::pool_linkages(list(q = list(q = spec)), env_data = env,
                                   strata = list(fleet = 1L))
  tbl <- pool$table
  re  <- tbl[!is.na(tbl$re_index), ]

  testthat::expect_equal(nrow(re), 10L)                 # one slot per year
  testthat::expect_setequal(re$re_index, 0:9)           # bijection
  testthat::expect_true(all(re$sigma_index == 0L))      # single group
  testthat::expect_false(anyNA(re$re_time))             # numeric time captured
  # fixed rows carry no registry entry
  testthat::expect_true(all(is.na(tbl$re_index[is.na(tbl$re_struct)])))
})

testthat::test_that("distinct fleets get distinct sigma groups", {
  env <- data.frame(Year = 2000:2004)
  s1 <- Rceattle::linkage_spec(~ (1 | Year), param = "q", by = ~ fleet, fleet = 1L)
  s2 <- Rceattle::linkage_spec(~ (1 | Year), param = "q", by = ~ fleet, fleet = 2L)
  pool <- Rceattle:::pool_linkages(list(q = list(q = list(s1, s2))),
                                   env_data = env, strata = list(fleet = 1:2))
  enc <- Rceattle:::encode_linkage_for_tmb(pool$table, pool$X)
  re  <- pool$table[!is.na(pool$table$re_index), ]

  testthat::expect_equal(max(re$sigma_index) + 1L, 2L)  # two groups
  # per-slot sigma group matches the owning row
  expect_sigma <- integer(nrow(re))
  expect_sigma[re$re_index + 1L] <- re$sigma_index
  testthat::expect_identical(enc$linkage_re_sigma, as.integer(expect_sigma))
})

testthat::test_that("a fixed-effect prior on an RE column is rejected", {
  env <- data.frame(Year = 2000:2004)
  spec <- Rceattle::linkage_spec(~ (1 | Year), param = "q", by = ~ fleet, fleet = 1L)
  pool <- Rceattle:::pool_linkages(list(q = list(q = spec)), env_data = env,
                                   strata = list(fleet = 1L))
  bad <- pool$table
  i <- which(!is.na(bad$re_struct))[1]
  bad$prior_family[i] <- "normal"; bad$prior_p1[i] <- 0; bad$prior_p2[i] <- 1
  testthat::expect_error(
    Rceattle:::encode_linkage_for_tmb(bad, pool$X), "fixed-effect prior")
})

testthat::test_that("env_data must align to the model years for linkages", {
  # Linkages apply env_data positionally (row r -> model year styr + r - 1), so
  # a Year column must start at styr and be contiguous.
  ok <- data.frame(Year = 1979:2000, temp = 0)
  testthat::expect_silent(Rceattle:::.check_env_data_years(ok, 1979L))
  # starts after styr
  testthat::expect_error(
    Rceattle:::.check_env_data_years(data.frame(Year = 1985:2000), 1979L),
    "must start at the model start year")
  # interior gap (1990 missing)
  testthat::expect_error(
    Rceattle:::.check_env_data_years(
      data.frame(Year = c(1979:1989, 1991:2000)), 1979L),
    "gaps")
  # a shorter-than-projection but aligned env_data is fine
  testthat::expect_silent(Rceattle:::.check_env_data_years(data.frame(Year = 1979:1985), 1979L))
  # no Year column -> positional contract assumed, no error
  testthat::expect_silent(Rceattle:::.check_env_data_years(data.frame(temp = 1:5), 1979L))
})

# ---- end-to-end density (skip_on_cran) --------------------------------------

.re_row <- function(fit) {
  which(rownames(fit$quantities$jnll_comp) == "Linkage random effects")
}

testthat::test_that("the RE density equals -sum(dnorm(dev, 0, sigma)) exactly", {
  # The strongest check of the density math, and it needs no optimisation:
  # build the object (estimateMode = 3), evaluate it once at a KNOWN deviate
  # pattern + SD, and confirm the reported density row is the mathematical
  # definition to machine precision. (The legacy M1_re = 2 density is the same
  # SCALE(AR1(rho = 0), sigma) form, i.e. the identical dnorm sum.)
  testthat::skip_on_cran()
  d <- Rceattle::BS2017SS
  qfun <- Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ (1 | Year), by = ~ fleet, fleet = 7L)))
  obj <- Rceattle::fit_mod(data_list = d, estimateMode = 3, msmMode = 0,
    qFun = qfun,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))$obj

  p <- obj$env$last.par
  n_re <- sum(names(p) == "beta_linkage_re")
  set.seed(1)
  dev  <- stats::rnorm(n_re, 0, 0.4)
  lsig <- log(0.3)
  p[names(p) == "beta_linkage_re"]   <- dev
  p[names(p) == "log_sigma_linkage"] <- lsig
  rep <- obj$report(p)
  # raw REPORT has no rownames; the density is the last jnll_comp row
  # (cpp row 20, 0-based).
  reported <- sum(rep$jnll_comp[linkage_re_row(), ])
  testthat::expect_equal(reported,
    -sum(stats::dnorm(dev, 0, exp(lsig), log = TRUE)), tolerance = 1e-8)
})

testthat::test_that("a model without a random linkage keeps density row 20 at 0", {
  testthat::skip_on_cran()
  fit <- Rceattle::fit_mod(
    data_list = Rceattle::BS2017SS, estimateMode = 1, msmMode = 0,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))
  testthat::expect_equal(sum(fit$quantities$jnll_comp[.re_row(fit), ]), 0)
})

testthat::test_that("an IID (1|Year) q linkage builds finite and estimates a density", {
  testthat::skip_on_cran()
  d <- Rceattle::BS2017SS                      # env_data already has Year
  survey_flt <- 7L  # EIT_Pollock, the only Estimated-q survey (Fixed q is rejected)
  qfun <- Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ (1 | Year), by = ~ fleet, fleet = survey_flt)))
  ctl <- Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0)

  # build-only: obj$fn / obj$gr finite (no NaN in the Laplace tape)
  b3 <- Rceattle::fit_mod(data_list = d, estimateMode = 3, msmMode = 0,
                          qFun = qfun, fit_control = ctl)
  testthat::expect_true(is.finite(b3$obj$fn()))
  testthat::expect_true(all(is.finite(b3$obj$gr())))

  # fit: one deviation per hindcast year, an active finite density, finite sigma
  fit <- Rceattle::fit_mod(data_list = d, estimateMode = 1, msmMode = 0,
                           qFun = qfun, fit_control = ctl)
  n_hind <- d$endyr - d$styr + 1L
  testthat::expect_equal(length(fit$estimated_params$beta_linkage_re), n_hind)
  testthat::expect_true(sum(fit$quantities$jnll_comp[.re_row(fit), ]) != 0)
  sigma <- exp(fit$estimated_params$log_sigma_linkage)
  testthat::expect_true(is.finite(sigma) && sigma > 0)
})

testthat::test_that("a phased (1|Year) q linkage holds its RE-SD and matches the direct fit", {
  testthat::skip_on_cran()
  # Regression for the phaser RE-SD hold-list: during phasing the deviates are
  # estimated as penalised fixed effects, so log_sigma_linkage must be HELD, not
  # estimated jointly with them. Without the hold, the near-saturated field and the
  # free SD reach a joint optimum that collapses sigma -> ~0; the final Laplace fit
  # stays trapped there and silently returns the RE-off objective (the linkage
  # variation vanishes) -- a wrong-quota misfit, not a crash. So the discriminating
  # check is that the phased sigma and objective MATCH the direct (unphased) fit.
  d <- Rceattle::BS2017SS
  qfun <- Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ (1 | Year), by = ~ fleet, fleet = 7L)))

  fit_ph <- Rceattle::fit_mod(data_list = d, estimateMode = 1, msmMode = 0,
    qFun = qfun, fit_control = Rceattle::fit_control(phase = TRUE, getsd = FALSE, verbose = 0))
  fit_dir <- Rceattle::fit_mod(data_list = d, estimateMode = 1, msmMode = 0,
    qFun = qfun, fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))

  sigma_ph  <- exp(fit_ph$estimated_params$log_sigma_linkage)
  sigma_dir <- exp(fit_dir$estimated_params$log_sigma_linkage)
  # sigma must NOT collapse (the bug drives it to ~1e-9) -- it recovers the direct MLE.
  testthat::expect_equal(sigma_ph, sigma_dir, tolerance = 1e-2)
  testthat::expect_gt(sigma_ph, 0.05)
  # phasing is a warm-start, not a different model: same optimum.
  testthat::expect_equal(fit_ph$opt$objective, fit_dir$opt$objective, tolerance = 1e-3)
})

testthat::test_that("rw(1|Year) fits, pins the first deviate, and uses a difference density", {
  testthat::skip_on_cran()
  d <- Rceattle::BS2017SS
  qfun <- Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ rw(1 | Year), by = ~ fleet, fleet = 7L)))
  fit <- Rceattle::fit_mod(data_list = d, estimateMode = 1, msmMode = 0,
    qFun = qfun, fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))
  dev   <- fit$estimated_params$beta_linkage_re
  sigma <- exp(fit$estimated_params$log_sigma_linkage)
  # first deviate pinned at 0 (the walk's level is carried by the base q)
  testthat::expect_equal(dev[1], 0, tolerance = 1e-10)
  # the density is N(0, sigma) on successive first differences (deviates are in
  # elapsed-time order)
  testthat::expect_equal(sum(fit$quantities$jnll_comp[.re_row(fit), ]),
    -sum(stats::dnorm(diff(dev), 0, sigma, log = TRUE)), tolerance = 1e-7)
})

testthat::test_that("ar1(1|Year) fits with a stationary AR1 density and estimable rho", {
  testthat::skip_on_cran()
  d <- Rceattle::BS2017SS
  qfun <- Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ ar1(1 | Year), by = ~ fleet, fleet = 7L)))
  fit <- Rceattle::fit_mod(data_list = d, estimateMode = 1, msmMode = 0,
    qFun = qfun, fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))
  re    <- fit$estimated_params$beta_linkage_re
  sigma <- exp(fit$estimated_params$log_sigma_linkage)
  rho   <- tanh(fit$estimated_params$trans_rho_linkage)   # rho_trans(x) == tanh(x)
  testthat::expect_true(is.finite(rho) && abs(rho) < 1)
  # exact stationary AR1 negative log density, MARGINAL SD = sigma:
  n <- length(re)
  nll <- -stats::dnorm(re[1], 0, sigma, log = TRUE) -
          sum(stats::dnorm(re[2:n], rho * re[1:(n - 1)], sigma * sqrt(1 - rho^2), log = TRUE))
  testthat::expect_equal(sum(fit$quantities$jnll_comp[.re_row(fit), ]), nll, tolerance = 1e-6)
})

testthat::test_that("init = list(rho = v) fixes the ar1 correlation at v", {
  testthat::skip_on_cran()
  d <- Rceattle::BS2017SS
  qfun <- Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ ar1(1 | Year), by = ~ fleet, fleet = 7L,
                               init = list(rho = 0.5))))
  fit <- Rceattle::fit_mod(data_list = d, estimateMode = 1, msmMode = 0,
    qFun = qfun, fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))
  testthat::expect_equal(tanh(fit$estimated_params$trans_rho_linkage), 0.5, tolerance = 1e-6)
})

testthat::test_that("an IID q linkage recovers a known deviation SD (simulate-refit self-consistency)", {
  testthat::skip_on_cran()
  # Beyond pinning the density VALUE (tests above): simulate a survey index from
  # KNOWN year deviations on log-q, refit, and confirm the estimated SD *and* the
  # deviations are recovered. This catches an estimation bias the value-pins cannot
  # -- e.g. a sigma_hat that ignores the truth (a real failure mode during dev).
  d    <- Rceattle::BS2017SS
  qfun <- Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ (1 | Year), by = ~ fleet, fleet = 7L)))
  ctl  <- Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0)
  obj  <- Rceattle::fit_mod(data_list = d, estimateMode = 1, msmMode = 0,
                            qFun = qfun, fit_control = ctl)$obj
  i7   <- which(d$index_data$Fleet_code == 7L)
  n_re <- sum(names(obj$env$last.par) == "beta_linkage_re")

  recover <- function(sig_true) {
    set.seed(7); z <- stats::rnorm(n_re, 0, sig_true)
    par <- obj$env$last.par.best
    par[names(par) == "beta_linkage_re"] <- z       # inject known deviations
    set.seed(7); sim <- obj$simulate(par)           # index generated from z (+ obs error)
    dsim <- d; dsim$index_data$Observation[i7] <- as.numeric(sim$index_hat)[i7]
    fit <- suppressWarnings(Rceattle::fit_mod(data_list = dsim, estimateMode = 1,
             msmMode = 0, qFun = qfun, fit_control = ctl))
    list(sig_hat = exp(fit$estimated_params$log_sigma_linkage), emp = stats::sd(z),
         cor = stats::cor(fit$estimated_params$beta_linkage_re, z))
  }
  lo <- recover(0.4); hi <- recover(0.7)

  # SD recovered within a generous band of the empirical SD of the injected devs
  testthat::expect_gt(lo$sig_hat / lo$emp, 0.55); testthat::expect_lt(lo$sig_hat / lo$emp, 1.5)
  testthat::expect_gt(hi$sig_hat / hi$emp, 0.55); testthat::expect_lt(hi$sig_hat / hi$emp, 1.5)
  testthat::expect_gt(hi$sig_hat, lo$sig_hat)       # sigma_hat tracks the truth
  testthat::expect_gt(lo$cor, 0.5); testthat::expect_gt(hi$cor, 0.5)  # devs recovered, not just spread
})

testthat::test_that("a (1|Year) q linkage with sigma pinned near 0 reproduces the no-linkage fit", {
  testthat::skip_on_cran()
  d   <- Rceattle::BS2017SS
  ctl <- Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0)
  base <- Rceattle::fit_mod(data_list = d, estimateMode = 1, msmMode = 0,
                            fit_control = ctl)$opt$objective
  pin  <- Rceattle::fit_mod(data_list = d, estimateMode = 1, msmMode = 0,
            qFun = Rceattle::build_catchability(linkages = list(
              q = Rceattle::linkage_spec(~ (1 | Year), by = ~ fleet, fleet = 7L,
                                         init = list(sigma = 1e-4)))),
            fit_control = ctl)
  testthat::expect_equal(pin$opt$objective, base, tolerance = 1e-3)
  testthat::expect_lt(max(abs(pin$estimated_params$beta_linkage_re)), 1e-3)
})

testthat::test_that("Rogers QAR1: observed ar1 adds an observation term and a beta effect", {
  # linkage_spec(observe=, obs_sd=) makes the ar1 latent a state-space covariate
  # (Rogers et al. 2024): row 20 carries the AR1 process density PLUS a Gaussian
  # observation of the latent, and the latent enters q scaled by an estimated
  # beta. Verified build-only (no optimisation) at a known parameter vector.
  testthat::skip_on_cran()
  d <- Rceattle::BS2017SS
  qfun <- Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ ar1(1 | Year), by = ~ fleet, fleet = 7L,
                               observe = "BTempC", obs_sd = 0.5)))
  o <- Rceattle::fit_mod(data_list = d, estimateMode = 3, msmMode = 0, qFun = qfun,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))$obj
  p <- o$env$last.par
  n_re <- sum(names(p) == "beta_linkage_re")
  set.seed(5)
  re <- stats::rnorm(n_re, 0, 0.4); lsig <- log(0.3); tr <- atanh(0.6); bobs <- 0.8
  p[names(p) == "beta_linkage_re"]   <- re
  p[names(p) == "log_sigma_linkage"] <- lsig
  p[names(p) == "trans_rho_linkage"] <- tr
  p[names(p) == "beta_linkage_obs"]  <- bobs
  rep <- o$report(p)
  sigma <- exp(lsig); rho <- tanh(tr); obs_series <- d$env_data$BTempC
  ar1_nll <- -stats::dnorm(re[1], 0, sigma, log = TRUE) -
    sum(stats::dnorm(re[2:n_re], rho * re[1:(n_re - 1)], sigma * sqrt(1 - rho^2), log = TRUE))
  obs_nll <- -sum(stats::dnorm(obs_series, re, 0.5, log = TRUE))
  testthat::expect_equal(sum(rep$jnll_comp[linkage_re_row(), ]),
                         ar1_nll + obs_nll, tolerance = 1e-6)
  # the latent enters q scaled by beta
  testthat::expect_equal(as.numeric(rep$q_linkage_offset[7, seq_len(n_re)]),
                         bobs * re, tolerance = 1e-8)
})

testthat::test_that("Rogers QAR1: obs_sd_est = TRUE estimates the observation SD", {
  # The build-only test above pins the density at a fixed obs_sd; this exercises
  # the ESTIMATION path (log_obs_sd_linkage free). The obs SD is weakly identified
  # and collapses toward 0 on a smooth covariate (documented degeneracy), so assert
  # the path runs and yields a finite, genuinely-estimated value -- not a target.
  testthat::skip_on_cran()
  d <- Rceattle::BS2017SS
  qfun <- Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ ar1(1 | Year), by = ~ fleet, fleet = 7L,
          observe = "BTempC", obs_sd = 0.5, obs_sd_est = TRUE)))
  fit <- suppressWarnings(Rceattle::fit_mod(data_list = d, estimateMode = 1,
           msmMode = 0, qFun = qfun,
           fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0)))
  lo <- fit$estimated_params$log_obs_sd_linkage
  testthat::expect_equal(length(lo), 1L)                       # one per observed group
  testthat::expect_true(is.finite(fit$opt$objective))
  testthat::expect_true(is.finite(exp(lo)) && exp(lo) > 0)
  testthat::expect_true(any(abs(lo - log(0.5)) > 1e-6))        # genuinely estimated, not held
})

testthat::test_that("observe requires an ar1 term and a valid column", {
  env <- data.frame(Year = 2000:2004, temp = stats::rnorm(5), fleet = 1L)
  testthat::expect_error(
    Rceattle:::materialize_linkage(
      Rceattle::linkage_spec(~ (1 | Year), param = "q", by = ~ fleet, fleet = 1L,
                             observe = "temp", obs_sd = 0.5),
      "q", env, list(fleet = 1L)),
    "requires an")
  testthat::expect_error(
    Rceattle::linkage_spec(~ ar1(1 | Year), param = "q", observe = "temp"),
    "obs_sd")
})

testthat::test_that("rw()/ar1() require a numeric grouping variable", {
  env <- data.frame(Year = letters[1:5])
  testthat::expect_error(
    Rceattle:::pool_linkages(list(q = list(q = Rceattle::linkage_spec(
      ~ rw(1 | Year), param = "q", by = ~ fleet, fleet = 1L))),
      env_data = env, strata = list(fleet = 1L)),
    "numeric grouping")
})

testthat::test_that("rw()/ar1() require equally-spaced time (reject gaps)", {
  # A gap (missing 2002) would be silently treated as a unit lag -- wrong
  # variance (rw) / correlation (ar1). It must error, not mis-specify.
  gap <- data.frame(Year = c(2000, 2001, 2003, 2004), fleet = 1L)
  for (form in c(~ rw(1 | Year), ~ ar1(1 | Year))) {
    testthat::expect_error(
      Rceattle:::pool_linkages(list(q = list(q = Rceattle::linkage_spec(
        form, param = "q", by = ~ fleet, fleet = 1L))),
        env_data = gap, strata = list(fleet = 1L)),
      "equally-spaced")
  }
  # a regularly-spaced (non-unit) grid is fine
  reg <- data.frame(Year = c(2000, 2002, 2004, 2006), fleet = 1L)
  testthat::expect_s3_class(
    Rceattle:::pool_linkages(list(q = list(q = Rceattle::linkage_spec(
      ~ ar1(1 | Year), param = "q", by = ~ fleet, fleet = 1L))),
      env_data = reg, strata = list(fleet = 1L))$table,
    "Rceattle_linkage_table")
})

testthat::test_that("init = list(sigma = v) fixes the deviation SD at v", {
  testthat::skip_on_cran()
  d <- Rceattle::BS2017SS
  survey_flt <- 7L  # EIT_Pollock, the only Estimated-q survey (Fixed q is rejected)
  qfun <- Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ (1 | Year), by = ~ fleet, fleet = survey_flt,
                               init = list(sigma = 0.1))))
  fit <- Rceattle::fit_mod(data_list = d, estimateMode = 1, msmMode = 0,
    qFun = qfun, fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))
  # the SD is held exactly at the input value, and the density is still active
  testthat::expect_equal(exp(fit$estimated_params$log_sigma_linkage), 0.1, tolerance = 1e-9)
  testthat::expect_true(sum(fit$quantities$jnll_comp[.re_row(fit), ]) != 0)
})

testthat::test_that("priors = list(sigma = ...) puts a prior on the deviation SD", {
  testthat::skip_on_cran()
  d <- Rceattle::BS2017SS
  survey_flt <- 7L  # EIT_Pollock, the only Estimated-q survey (Fixed q is rejected)
  qfun <- Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ (1 | Year), by = ~ fleet, fleet = survey_flt,
                               priors = list(sigma = lognormal(log(0.05), 0.3)))))
  fit <- Rceattle::fit_mod(data_list = d, estimateMode = 1, msmMode = 0,
    qFun = qfun, fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))
  pr_row <- which(rownames(fit$quantities$jnll_comp) == "Linkage-table priors")
  testthat::expect_true(sum(fit$quantities$jnll_comp[pr_row, ]) != 0)   # prior contributes
  testthat::expect_true(is.finite(exp(fit$estimated_params$log_sigma_linkage)))
})


# fit_mod()'s discontinuous-likelihood check compares nlminb's objective against
# a fresh evaluation. Both must be the MARGINAL objective: obj$report()$jnll is
# the joint nll at the random-effect mode, a Laplace correction away, which made
# the warning unconditional on every model carrying a random effect.
testthat::test_that("a smooth random-effect model is not called discontinuous", {
  testthat::skip_on_cran()
  d <- Rceattle::BS2017SS
  survey_flt <- 7L  # EIT_Pollock, the only Estimated-q survey (Fixed q is rejected)
  qfun <- Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ (1 | Year), by = ~ fleet, fleet = survey_flt)))

  # getsd = TRUE: the check only runs where an sdreport exists.
  msgs <- testthat::capture_messages(
    fit <- Rceattle::fit_mod(data_list = d, estimateMode = 1, msmMode = 0,
      qFun = qfun,
      fit_control = Rceattle::fit_control(phase = FALSE, getsd = TRUE, verbose = 0)))

  testthat::expect_true(length(fit$obj$env$random) > 0)   # the check is reachable
  testthat::expect_false(any(grepl("discontinuous", msgs)))

  # The joint nll is more than rel_tol away, so the old comparison flagged this
  # model every time. The marginal gap is under it.
  jnll <- fit$obj$report(fit$obj$env$last.par.best)$jnll
  testthat::expect_gt(abs(fit$opt$objective - jnll), 1)
  testthat::expect_lt(Rceattle:::.rce_marginal_gap(fit$obj, fit$opt), 1)
})


# The gap has to be measured against the object the optimizer ran on.
# fit_mod() rebuilds obj for the projection, and again under
# projection_uncertainty, and build_hcr_map() maps every hindcast parameter off,
# so those rebuilds disagree with `opt` on which parameters are random. Measured
# on this model: opt$objective 10207.488519 against 10228.940923 from the
# all-parameters-on rebuild, a 21.45 gap on a fit with nothing wrong with it.
testthat::test_that("projection_uncertainty does not fake a discontinuity", {
  testthat::skip_on_cran()
  d <- Rceattle::BS2017SS
  d$fleet_control$Proj_F_proportion[d$fleet_control$Fleet_type == 1] <- 1
  qfun <- Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ (1 | Year), by = ~ fleet, fleet = 7L)))

  msgs <- testthat::capture_messages(
    fit <- Rceattle::fit_mod(data_list = d, estimateMode = 0, msmMode = 0,
      qFun = qfun,
      HCR = Rceattle::build_hcr(HCR = "NPFMC", Ptarget = 0.4, Plimit = 0.2),
      fit_control = Rceattle::fit_control(phase = FALSE, getsd = TRUE,
                                          projection_uncertainty = TRUE,
                                          verbose = 0)))

  # The rebuild carries the random effects the projection fit did not.
  testthat::expect_true(length(fit$obj$env$random) > 0)
  testthat::expect_gt(Rceattle:::.rce_marginal_gap(fit$obj, fit$opt), 1)
  testthat::expect_false(any(grepl("discontinuous", msgs)))
})
