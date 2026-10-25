# =============================================================================
# Intercept-row priors must be re-targeted to the BASE parameter and contribute
# to the likelihood (jnll slot "Linkage-table priors" / row 20).
#
# The existing prior tests (test-growth-linkage-species.R,
# test-linkage-priors.R) only exercise priors on SLOPE terms (temp, PDO),
# where the prior is evaluated directly on beta_linkage. For an (Intercept) row
# the linkage coefficient is mapped out (beta_linkage(i) == 0) and the level is
# carried by the base parameter (log_growth_pars / log_M1). The cpp (section 13.13 /
# row 20) re-targets the prior onto that base parameter and back-transforms
# through the log link:
#   normal/gamma/beta : evaluated on exp(base_par)   (natural scale)
#   lognormal         : evaluated on base_par        (already on the log scale)
# These tests pin that behaviour so a prior on K / L1 / Linf / M actually
# constrains the realized parameter (and isn't silently evaluated at exp(0)=1).
#
# Fixtures use nspp = 2 to match make_msm_test_data's defaults.
# =============================================================================

# Integration test (fits a CEATTLE model): skipped on CRAN to keep R CMD check
# fast; the coverage workflow runs it in full (NOT_CRAN=true). See README.md.
testthat::skip_on_cran()

# jnll row for the linkage-table priors, robust to row-name vs index.
linkage_prior_nll <- function(fit) {
  jc <- fit$quantities$jnll_comp
  rn <- rownames(jc)
  if (!is.null(rn) && "Linkage-table priors" %in% rn) {
    return(sum(jc["Linkage-table priors", ]))
  }
  sum(jc[20, ])
}

make_two_spp <- function(seed = 101) {
  set.seed(seed)
  nyrs <- 30
  Fmort <- matrix(c(seq(0.02, 0.3, length.out = nyrs / 2),
                    seq(0.3, 0.05, length.out = nyrs / 2),
                    seq(0.02, 0.3, length.out = nyrs)),
                  nrow = 2, ncol = nyrs, byrow = TRUE)
  make_msm_test_data(years = 1:nyrs, Fmort = Fmort,
                     log_phi = matrix(-Inf, 2, 2))$data_list
}


testthat::test_that("growth (Intercept) priors are evaluated on log_growth_pars (natural scale)", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  sim_data <- make_two_spp(101)

  # Per-species base levels for the estimable VB params, with normal priors on
  # the NATURAL scale centred away from the init so the NLL is non-trivial.
  Linf1 <- 90; Linf2 <- 50; Linf_mu <- 95; Linf_sd <- 3
  growth_spec <- Rceattle::build_growth(
    fun = "vonBertalanffy",
    linkages = list(
      Linf = list(
        Rceattle::linkage_spec(
          formula = ~ 1, by = ~ species, species = 1L,
          init   = list("(Intercept)" = Linf1),
          priors = list("(Intercept)" = normal(Linf_mu, Linf_sd))),
        Rceattle::linkage_spec(
          formula = ~ 1, by = ~ species, species = 2L,
          init   = list("(Intercept)" = Linf2),
          priors = list("(Intercept)" = normal(Linf_mu, Linf_sd)))
      )
    )
  )

  fit <- suppressMessages(Rceattle::fit_mod(
    data_list = sim_data, growthFun = growth_spec, estimateMode = 3,
    msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0)
  ))

  # 1. The base parameter carries the (Intercept) level on the log scale.
  lgp <- fit$estimated_params$log_growth_pars
  testthat::expect_equal(exp(as.numeric(lgp[1, 1, 3])), Linf1, tolerance = 1e-8)
  testthat::expect_equal(exp(as.numeric(lgp[2, 1, 3])), Linf2, tolerance = 1e-8)

  # 2. (Intercept) coefficients in beta_linkage are mapped out at 0.
  testthat::expect_true(all(fit$estimated_params$beta_linkage == 0))

  # 3. Prior NLL evaluated on the NATURAL-SCALE base parameter, summed over spp.
  expected <- -(stats::dnorm(Linf1, Linf_mu, Linf_sd, log = TRUE) +
                stats::dnorm(Linf2, Linf_mu, Linf_sd, log = TRUE))
  testthat::expect_equal(linkage_prior_nll(fit), expected, tolerance = 1e-6)

  # 4. Negative control: NOT evaluated on the mapped-out beta_linkage = 0
  #    (which would give exp(0) = 1).
  broken <- -2 * stats::dnorm(1, Linf_mu, Linf_sd, log = TRUE)
  testthat::expect_false(isTRUE(all.equal(linkage_prior_nll(fit), broken)))
})


testthat::test_that("lognormal (Intercept) prior is evaluated on the log-scale base parameter", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  sim_data <- make_two_spp(102)
  K1 <- 0.25; K2 <- 0.30; K_meanlog <- log(0.20); K_sdlog <- 0.30

  growth_spec <- Rceattle::build_growth(
    fun = "vonBertalanffy",
    linkages = list(
      K = list(
        Rceattle::linkage_spec(
          formula = ~ 1, by = ~ species, species = 1L,
          init = list("(Intercept)" = K1),
          priors = list("(Intercept)" = lognormal(K_meanlog, K_sdlog))),
        Rceattle::linkage_spec(
          formula = ~ 1, by = ~ species, species = 2L,
          init = list("(Intercept)" = K2),
          priors = list("(Intercept)" = lognormal(K_meanlog, K_sdlog)))
      )
    )
  )

  fit <- suppressMessages(Rceattle::fit_mod(
    data_list = sim_data, growthFun = growth_spec, estimateMode = 3,
    msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0)
  ))

  # lognormal -> dnorm(log(b_nat), meanlog - sdlog^2/2, sdlog) on log_growth_pars,
  # mean-centred under the default bias_adjust_proc = TRUE.
  mu <- K_meanlog - K_sdlog^2 / 2
  expected <- -(stats::dnorm(log(K1), mu, K_sdlog, log = TRUE) +
                stats::dnorm(log(K2), mu, K_sdlog, log = TRUE))
  testthat::expect_equal(linkage_prior_nll(fit), expected, tolerance = 1e-6)
})


testthat::test_that("M (Intercept) prior is evaluated on log_M1 (natural scale)", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  sim_data <- make_two_spp(103)
  M1v <- 0.30; M2v <- 0.20; M_mu <- 0.25; M_sd <- 0.05

  M1Fun <- Rceattle::build_M1(
    M1_model = 1,
    linkages = list(M1 = list(
      Rceattle::linkage_spec(
        formula = ~ 1, by = ~ species, species = 1L,
        init = list("(Intercept)" = M1v),
        priors = list("(Intercept)" = normal(M_mu, M_sd))),
      Rceattle::linkage_spec(
        formula = ~ 1, by = ~ species, species = 2L,
        init = list("(Intercept)" = M2v),
        priors = list("(Intercept)" = normal(M_mu, M_sd)))
    ))
  )

  fit <- suppressMessages(Rceattle::fit_mod(
    data_list = sim_data, M1Fun = M1Fun, estimateMode = 3,
    msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0)
  ))

  # Base M carries the (Intercept) level on the log scale (per species).
  testthat::expect_equal(exp(as.numeric(fit$estimated_params$log_M1[1, 1, 1])),
                         M1v, tolerance = 1e-8)
  testthat::expect_equal(exp(as.numeric(fit$estimated_params$log_M1[2, 1, 1])),
                         M2v, tolerance = 1e-8)
  expected <- -(stats::dnorm(M1v, M_mu, M_sd, log = TRUE) +
                stats::dnorm(M2v, M_mu, M_sd, log = TRUE))
  testthat::expect_equal(linkage_prior_nll(fit), expected, tolerance = 1e-6)
})


testthat::test_that("moving the base parameter changes the (Intercept) prior NLL", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  sim_data <- make_two_spp(104)
  Linf_mu <- 95; Linf_sd <- 3
  growth_spec <- Rceattle::build_growth(
    fun = "vonBertalanffy",
    linkages = list(
      Linf = list(
        Rceattle::linkage_spec(
          formula = ~ 1, by = ~ species, species = 1L,
          init = list("(Intercept)" = 90),
          priors = list("(Intercept)" = normal(Linf_mu, Linf_sd))),
        Rceattle::linkage_spec(
          formula = ~ 1, by = ~ species, species = 2L,
          init = list("(Intercept)" = 90),
          priors = list("(Intercept)" = normal(Linf_mu, Linf_sd)))
      )
    )
  )

  fit0 <- suppressMessages(Rceattle::fit_mod(
    data_list = sim_data, growthFun = growth_spec, estimateMode = 3,
    msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0)
  ))

  # Move both base Linf to exactly the prior mean -> prior NLL hits its minimum.
  inits <- fit0$estimated_params
  inits$log_growth_pars[1, 1, 3] <- log(Linf_mu)
  inits$log_growth_pars[2, 1, 3] <- log(Linf_mu)
  fit1 <- suppressMessages(Rceattle::fit_mod(
    data_list = sim_data, inits = inits, growthFun = growth_spec,
    estimateMode = 3, msmMode = 0, random_rec = FALSE,
    fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0)
  ))

  testthat::expect_lt(linkage_prior_nll(fit1), linkage_prior_nll(fit0))
  testthat::expect_equal(linkage_prior_nll(fit1),
                         -2 * stats::dnorm(Linf_mu, Linf_mu, Linf_sd, log = TRUE),
                         tolerance = 1e-6)
})
