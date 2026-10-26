# Integration test (fits a CEATTLE model): skipped on CRAN to keep R CMD check
# fast; the coverage workflow runs it in full (NOT_CRAN=true). See README.md.
testthat::skip_on_cran()

testthat::test_that("Test SB0 under mean recruitment", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Use helper factory for deterministic minimal test data
  #source(file.path("tests", "testthat", "helpers.R"))
  nages = 5
  R0 = 10
  nyrs = 8
  dat <- make_test_data(nyrs = nyrs, nprojyrs = 20, nages = nages, seed = 42)

  # Set params
  ss_run <- suppressMessages(
    Rceattle::fit_mod(data_list = dat,
                      estimateMode = 3, # Don't estimate
                      msmMode = 0, # Single species mode
                      fit_control = fit_control(getsd = FALSE, 
                        verbose = 1))
  )
  inits <- ss_run$estimated_params
  inits$rec_pars[,1] <- R0
  inits$R_log_sd <- 0
  inits$log_F[] <- -999 # No fishing
  inits$index_log_q[] <- 0 # Set q to 1

  # Set logistic params
  inits$log_sel_slp[] <- -Inf
  inits$sel_inf[] <- 0     # Females
  inits$rec_dev[] <- 0.02 # Adding .1 to R0

  # Run
  ss_run <- Rceattle::fit_mod(data_list = dat,
                              inits = inits, # Initial parameters from starting
                              file = NULL, # Don't save
                              recFun = build_srr(
                                proj_mean_rec = TRUE, # Project using mean rec over hindcast
                              ),
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              fit_control = fit_control(getsd = FALSE, 
                                verbose = 0))

  # Calculate SB0
  M <- dat$M1_base$Age1
  wt <- dat$weight$Age1
  sex_ratio <- dat$sex_ratio$Age1
  pmature <- dat$maturity$Age1

  natage = rep(1, nages)
  for(age in 2:nages){
    natage[age] <- natage[age-1] * exp(-M)
  }

  natage[nages] <-  natage[nages-1] * as.numeric(exp(-M)) / as.numeric(1-exp(-M))
  SB0 <- sum(natage) * exp(R0 + 0.02) * sex_ratio

  # Check SB0
  testthat::expect_equal(as.numeric(ss_run$quantities$SB0[1,nyrs + 20]), SB0, tolerance = 0.0001)
})

testthat::test_that("Test SB0 under R0", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Use helper factory for deterministic minimal test data
  #source(file.path("tests", "testthat", "helpers.R"))
  nages = 5
  R0 = 10
  nyrs = 8
  dat <- make_test_data(nyrs = nyrs, nprojyrs = 20, nages = nages, seed = 42)

  # Set params
  ss_run <- suppressMessages(
    Rceattle::fit_mod(data_list = dat,
                      estimateMode = 3, # Don't estimate
                      msmMode = 0, # Single species mode
                      fit_control = fit_control(getsd = FALSE, 
                        verbose = 1))
  )
  inits <- ss_run$estimated_params
  inits$rec_pars[,1] <- R0
  inits$R_log_sd <- 0
  inits$log_F[] <- -999 # No fishing
  inits$index_log_q[] <- 0 # Set q to 1

  # Set logistic params
  inits$log_sel_slp[] <- -Inf
  inits$sel_inf[] <- 0     # Females
  inits$rec_dev[] <- 0.02 # Adding .1 to R0

  # Run
  ss_run <- Rceattle::fit_mod(data_list = dat,
                              inits = inits, # Initial parameters from starting
                              file = NULL, # Don't save
                              recFun = build_srr(
                                proj_mean_rec = FALSE, # Project using R0 or SRR over hindcast
                              ),
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              fit_control = fit_control(getsd = FALSE, 
                                verbose = 0))

  # Calculate SB0
  M <- dat$M1_base$Age1
  wt <- dat$weight$Age1
  sex_ratio <- dat$sex_ratio$Age1
  pmature <- dat$maturity$Age1

  natage = rep(1, nages)
  for(age in 2:nages){
    natage[age] <- natage[age-1] * exp(-M)
  }

  natage[nages] <-  natage[nages-1] * as.numeric(exp(-M)) / as.numeric(1-exp(-M))
  SB0 <- sum(natage) * exp(R0) * sex_ratio

  # Check SB0
  testthat::expect_equal(as.numeric(ss_run$quantities$SB0[1,nyrs + 20]), SB0, tolerance = 0.0001)
})


testthat::test_that("Test SB40 under mean recruitment", {
  data("BS2017SS") # ?BS2017SS for more information on the data
  BS2017SS$fleet_control$Proj_F_proportion <- rep(1, 7)

  # -- NPFMC Tier 3
  ss_run_Tier3 <- Rceattle::fit_mod(data_list = BS2017SS,
                                    inits = NULL, # Initial parameters from ss_run
                                    estimateMode = 0, # Fit the hindcast and the HCR projection
                                    HCR = build_hcr(HCR = 5, # Tier3 HCR
                                                    Ftarget = 0.4, # F40%
                                                    Flimit = 0.35, # F35%
                                                    Plimit = 0.2, # No fishing when SB<SB20
                                                    Alpha = 0.05),
                                    msmMode = 0, # Single species mode
                                    fit_control = fit_control(getsd = FALSE, 
                                      phase = TRUE,
                                      verbose = 1))

  testthat::expect_equal(as.numeric((ss_run_Tier3$quantities$SBF/ss_run_Tier3$quantities$SB0)[,72]), rep(0.4, 3), tolerance = 0.0001)
  testthat::expect_equal(as.numeric((ss_run_Tier3$quantities$ssb/ss_run_Tier3$quantities$SB0)[,72]), rep(0.4, 3), tolerance = 0.0001) # Default to mean R
})

testthat::test_that("Test SB40 under R0", {
  data("BS2017SS") # ?BS2017SS for more information on the data
  BS2017SS$fleet_control$Proj_F_proportion <- rep(1, 7)

  # -- NPFMC Tier 3
  ss_run_Tier3 <- Rceattle::fit_mod(data_list = BS2017SS,
                                    inits = NULL, # Initial parameters from ss_run
                                    estimateMode = 0, # Fit the hindcast and the HCR projection
                                    HCR = build_hcr(HCR = 5, # Tier3 HCR
                                                    Ftarget = 0.4, # F40%
                                                    Flimit = 0.35, # F35%
                                                    Plimit = 0.2, # No fishing when SB<SB20
                                                    Alpha = 0.05),
                                    recFun = build_srr(proj_mean_rec = FALSE),
                                    msmMode = 0, # Single species mode
                                    fit_control = fit_control(getsd = FALSE, 
                                      phase = TRUE,
                                      verbose = 1))

  testthat::expect_equal(as.numeric((ss_run_Tier3$quantities$SBF/ss_run_Tier3$quantities$SB0)[,72]), rep(0.4, 3), tolerance = 0.0001)
  testthat::expect_equal(as.numeric((ss_run_Tier3$quantities$ssb/ss_run_Tier3$quantities$SB0)[,72]), rep(0.4, 3), tolerance = 0.0001)
})


testthat::test_that("Test SPR0 calculation", {
  library(Rceattle)
  data("GOA2018SS") # Single-species data. ?BS2017SS for more information on the data

  # Change WT and M for SSB
  GOA2018SS$weight[,6:35] <- 1
  GOA2018SS$M1_base[,3:32] <- 0.2
  GOA2018SS$maturity[,-1] <- 1
  GOA2018SS$sex_ratio[,-1] <- 0.5
  GOA2018SS$spawn_month <- rep(0,3)

  # Specify recruitment
  yrs <- GOA2018SS$styr:GOA2018SS$endyr
  nyrs <- length(yrs)
  R0 = 10:12
  Rdev <- stats::rnorm(nyrs)

  # Set params
  ss_run <- suppressMessages(
    Rceattle::fit_mod(data_list = GOA2018SS,
                      estimateMode = 3, # Don't estimate
                      msmMode = 0, # Single species mode
                      fit_control = fit_control(getsd = FALSE, 
                        verbose = 1))
  )
  inits <- ss_run$estimated_params
  inits$rec_pars[,1] <- R0

  # Run
  ss_run <- Rceattle::fit_mod(data_list = GOA2018SS,
                              inits = inits, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              initMode = "NonEquilibrium",
                              fit_control = fit_control(getsd = FALSE, 
                                verbose = 1))

  # Calculate SPR
  M <- 0.2
  wt <- 1
  sex_ratio <- 0.5
  pmature <- 1

  natage = rep(1, 10)
  for(age in 2:9){
    natage[age] <- natage[age-1] * exp(-M)
  }

  natage[10] <-  natage[9] * as.numeric(exp(-M)) / as.numeric(1-exp(-M))
  ssb_at_age <- natage * wt * pmature * sex_ratio
  sum(ssb_at_age)

  # Check SPRO
  testthat::expect_equal(as.numeric(ss_run$quantities$SPR0),  rep(sum(ssb_at_age), 3), tolerance = 0.0001)
})


# Female SSB rebuilt from the reported numbers- and weight-at-age, for all three
# species and all 74 years, hindcast AND projection. Two terms no other test
# could see, each with the error it produces if dropped:
#
#   * The spawn-month discount `exp(-Z * spawn_month/12)`. GOA pollock spawns in
#     mid-March (spawn_month 2.52); removing the factor moves its SSB by up to
#     24.8%. Species 2 and 3 spawn at month 0, where the factor is identically 1,
#     so they are the controls that show the term is read PER SPECIES.
#   * `spawn_output`'s sex branch (`ceattle.cpp` ~1166 and ~940): `mature_females`
#     is maturity times sex_ratio on a ONE-sex species and maturity alone on a
#     two-sex one. Arrowtooth is GOA2018SS's two-sex species; applying sex_ratio
#     there too moves its SSB by 182%.
#
# Projection years take `spawn_output(sp, age, nyrs_hind - 1)` rather than year
# `yr` (`ceattle.cpp` ~2442 against ~2088). That is NOT what this test pins:
# `weight_hat` already carries the last hindcast year forward, so its projection
# columns equal column `nh` exactly here, and the two readings agree. Using year
# `y` throughout is therefore correct and the distinction is invisible on this
# fixture.
testthat::test_that("SSB is taken at the spawn month in every year and species", {
  testthat::skip_if_not_installed("TMB")
  g <- Rceattle::GOA2018SS
  fit <- suppressMessages(suppressWarnings(
    Rceattle::fit_mod(data_list = g, estimateMode = 3, msmMode = 0,
                      fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0))))
  q  <- fit$quantities
  nh <- g$endyr - g$styr + 1
  ny <- ncol(q$ssb)
  testthat::expect_gt(ny, nh)            # there ARE projection years to cover
  # By name, not `-1`: the age columns are the convention here, but `M1_base` is
  # `Species, Sex, Age1...`, so a leading-column count is not portable in this
  # data_list.
  agecols <- grep("^Age", names(g$maturity))
  testthat::expect_gt(length(agecols), 0L)

  for (sp in seq_len(g$nspp)) {
    na  <- g$nages[sp]
    sm  <- g$spawn_month[sp]
    mat <- as.numeric(g$maturity[g$maturity[["Species"]] == sp, agecols])[1:na]
    sr  <- as.numeric(g$sex_ratio[g$sex_ratio[["Species"]] == sp, agecols])[1:na]
    spawners <- if (g$nsex[sp] == 1) mat * sr else mat
    # weight_hat row 2*sp is this species' spawning weight: wt_idx_ssb is
    # `2 * sp + 1` on a 0-based species index.
    want <- vapply(seq_len(ny), function(y) sum(
      q$N_at_age[sp, 1, 1:na, y] * exp(-q$Z_at_age[sp, 1, 1:na, y] * sm / 12) *
        q$weight_hat[2 * sp, 1, 1:na, y] * spawners), numeric(1))
    testthat::expect_equal(unname(q$ssb[sp, ]), want, tolerance = 1e-10,
                           info = paste("species", sp))
  }

  # The premise for the species that carries the timing term: without it the
  # reconstruction is wrong, so the assertion above is not passing on a fixture
  # where the discount is inert.
  na <- g$nages[1]
  mat <- as.numeric(g$maturity[g$maturity[["Species"]] == 1, agecols])[1:na]
  sr  <- as.numeric(g$sex_ratio[g$sex_ratio[["Species"]] == 1, agecols])[1:na]
  testthat::expect_gt(g$spawn_month[1], 0)
  undiscounted <- vapply(seq_len(ny), function(y) sum(
    q$N_at_age[1, 1, 1:na, y] * q$weight_hat[2, 1, 1:na, y] * mat * sr), numeric(1))
  testthat::expect_gt(max(abs(unname(q$ssb[1, ]) / undiscounted - 1)), 0.1)
  # And the two-sex species really is two-sex, so its branch was exercised.
  testthat::expect_equal(unname(g$nsex), c(1, 2, 1))
})


testthat::test_that("Test hindcast the same across different HCRs/BRPs", {
  data("BS2017SS") # ?BS2017SS for more information on the data
  BS2017SS$fleet_control$Proj_F_proportion <- rep(1, 7)

  # -- Base
  ss_run <- Rceattle::fit_mod(data_list = BS2017SS,
                              estimateMode = 0,
                              msmMode = 0, # Single species mode
                              fit_control = fit_control(getsd = FALSE, 
                                phase = TRUE,
                                verbose = 0))

  # -- Constant F as a percentage of SB0
  ss_run_fb0 <- Rceattle::fit_mod(data_list = BS2017SS,
                                  inits = ss_run$estimated_params, # Initial parameters from ss_run
                                  estimateMode = 0, # Fit the hindcast and the HCR projection
                                  HCR = build_hcr(HCR = 3, # F achieving a target % of SSB0
                                                  DynamicHCR = FALSE, # Static reference points
                                                  Ftarget = 0.4), # F that achieves 40% SB0
                                  msmMode = 0, # Single species mode
                                  fit_control = fit_control(getsd = FALSE, 
                                    verbose = 0))


  ss_run_dynamicfb0 <- Rceattle::fit_mod(data_list = BS2017SS,
                                         inits = ss_run$estimated_params, # Initial parameters from ss_run
                                         estimateMode = 0, # Fit the hindcast and the HCR projection
                                         HCR = build_hcr(HCR = 3, # F achieving a target % of SSB0
                                                         DynamicHCR = TRUE, # Use dynamic reference points
                                                         Ftarget = 0.4), # F that achieves 40% SB0
                                         msmMode = 0, # Single species mode
                                         fit_control = fit_control(getsd = FALSE, 
                                           verbose = 0))


  # -- Constant Fspr
  ss_run_Fspr <- Rceattle::fit_mod(data_list = BS2017SS,
                                   inits = ss_run$estimated_params, # Initial parameters from ss_run
                                   estimateMode = 0, # Fit the hindcast and the HCR projection
                                   HCR = build_hcr(HCR = 4, # Constant target Fspr
                                                   Ftarget = 0.4 # F40%
                                   ),
                                   msmMode = 0, # Single species mode
                                   fit_control = fit_control(getsd = FALSE, 
                                     verbose = 0))


  ss_run_dynamicFspr <- Rceattle::fit_mod(data_list = BS2017SS,
                                          inits = ss_run$estimated_params, # Initial parameters from ss_run
                                          estimateMode = 0, # Fit the hindcast and the HCR projection
                                          HCR = build_hcr(HCR = 4, # Constant target Fspr
                                                          DynamicHCR = TRUE, # Use dynamic reference points
                                                          Ftarget = 0.4 # F40%
                                          ),
                                          msmMode = 0, # Single species mode
                                          fit_control = fit_control(getsd = FALSE, 
                                            verbose = 0))


  # -- NPFMC Tier 3
  ss_run_Tier3 <- Rceattle::fit_mod(data_list = BS2017SS,
                                    inits = ss_run$estimated_params, # Initial parameters from ss_run
                                    estimateMode = 0, # Fit the hindcast and the HCR projection
                                    HCR = build_hcr(HCR = 5, # Tier3 HCR
                                                    Ftarget = 0.4, # F40%
                                                    Flimit = 0.35, # F35%
                                                    Plimit = 0.2, # No fishing when SB<SB20
                                                    Alpha = 0.05),
                                    msmMode = 0, # Single species mode
                                    fit_control = fit_control(getsd = FALSE, 
                                      verbose = 0))


  ss_run_dynamicTier3 <- Rceattle::fit_mod(data_list = BS2017SS,
                                           inits = ss_run$estimated_params, # Initial parameters from ss_run
                                           estimateMode = 0, # Fit the hindcast and the HCR projection
                                           HCR = build_hcr(HCR = 5, # Tier3 HCR
                                                           DynamicHCR = TRUE, # Use dynamic reference points
                                                           Ftarget = 0.4, # F40%
                                                           Flimit = 0.35, # F35%
                                                           Plimit = 0.2, # No fishing when SB<SB20
                                                           Alpha = 0.05),
                                           msmMode = 0, # Single species mode
                                           fit_control = fit_control(getsd = FALSE, 
                                             verbose = 0))

  # -- PFMC Category 1
  ss_run_Cat1 <- Rceattle::fit_mod(data_list = BS2017SS,
                                   inits = ss_run$estimated_params, # Initial parameters from ss_run
                                   estimateMode = 0, # Fit the hindcast and the HCR projection
                                   HCR = build_hcr(HCR = 6, # Cat 1 HCR
                                                   Flimit = 0.45, # F45%
                                                   Ptarget = 0.4, # Target is 40% B0
                                                   Plimit = 0.1, # No fishing when SB<SB10
                                                   Pstar = 0.45,
                                                   Sigma = 0.5),
                                   msmMode = 0, # Single species mode
                                   fit_control = fit_control(getsd = FALSE, 
                                     verbose = 0))

  ss_run_dynamicCat1 <- Rceattle::fit_mod(data_list = BS2017SS,
                                          inits = ss_run$estimated_params, # Initial parameters from ss_run
                                          estimateMode = 0, # Fit the hindcast and the HCR projection
                                          HCR = build_hcr(HCR = 6, # Cat 1 HCR
                                                          DynamicHCR = TRUE, # Use dynamic reference points
                                                          Flimit = 0.45, # F45%
                                                          Ptarget = 0.4, # Target is 40% SB0
                                                          Plimit = 0.1, # No fishing when SB<SB10
                                                          Pstar = 0.45,
                                                          Sigma = 0.5),
                                          msmMode = 0, # Single species mode
                                          fit_control = fit_control(getsd = FALSE, 
                                            verbose = 0))

  # -- SESSF Tier 1
  ss_run_Tier1 <- Rceattle::fit_mod(data_list = BS2017SS,
                                    inits = ss_run$estimated_params, # Initial parameters from ss_run
                                    estimateMode = 0, # Fit the hindcast and the HCR projection
                                    HCR = build_hcr(HCR = 7, # Tier 1 HCR
                                                    Ftarget = 0.48, # F40%
                                                    Flimit = 0.20, # F20%
                                                    Ptarget = 0.35, # Target is 35% SSB0
                                                    Plimit = 0.20, # No fishing when B<B20
                                    ),
                                    msmMode = 0, # Single species mode
                                    fit_control = fit_control(getsd = FALSE, 
                                      verbose = 0))


  ss_run_dynamicTier1 <- Rceattle::fit_mod(data_list = BS2017SS,
                                           inits = ss_run$estimated_params, # Initial parameters from ss_run
                                           estimateMode = 0, # Fit the hindcast and the HCR projection
                                           HCR = build_hcr(HCR = 7, # Tier 1 HCR
                                                           DynamicHCR = TRUE,
                                                           Ftarget = 0.48, # F40%
                                                           Flimit = 0.20, # F20%
                                                           Ptarget = 0.35, # Target is 35% SSB0
                                                           Plimit = 0.20, # No fishing when B<B20
                                           ),
                                           msmMode = 0, # Single species mode
                                           fit_control = fit_control(getsd = FALSE, 
                                             verbose = 0))

  # -- Plot
  mod_list <- list(ss_run, ss_run_fb0, ss_run_Fspr, ss_run_Tier3, ss_run_Cat1, ss_run_Tier1, ss_run_dynamicfb0, ss_run_dynamicFspr, ss_run_dynamicTier3, ss_run_dynamicCat1, ss_run_dynamicTier1)
  nyrs <- length(BS2017SS$styr:BS2017SS$endyr)
  nyrs_proj <- length(BS2017SS$styr:BS2017SS$projyr)
  proj_yrs <- (nyrs+2):nyrs_proj

  # Test hindcast is the same
  for(i in 2:length(mod_list)){
    testthat::expect_equal(ss_run$quantities$ssb[,1:nyrs], mod_list[[i]]$quantities$ssb[,1:nyrs])
  }

  # Test forecast is different
  for(i in 2:length(mod_list)){
    expect_all_true(c(ss_run$quantities$ssb[,proj_yrs] != mod_list[[i]]$quantities$ssb[,proj_yrs]))
  }
})


testthat::test_that("Test mean recruitment calculation", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Prepare small deterministic dataset using helper
  #source(file.path("tests", "testthat", "helpers.R"))

  # 1) Set up simulation
  nyrs = 30
  nspp = 2
  Fmort <- c(seq(0.02, 0.3, length.out = nyrs/2), seq(0.3, 0.05, length.out = nyrs/2))
  Fmort2 <- seq(0.02, 0.3, length.out = nyrs)
  log_phi = matrix(-Inf, nspp, nspp, byrow = TRUE)

  # First, simulate some data for the model
  set.seed(123)
  sim <- make_msm_test_data(
    years = 1:nyrs,
    Fmort = matrix(c(Fmort, Fmort2), nspp, nyrs, byrow = TRUE),

    # Multispecies bits
    log_phi = log_phi # set to -Inf so no predation
  )

  # sum(sim$model_quantities$B_eaten_as_prey)

  # Set up Rceattle data
  simData <- sim$data_list


  # Fit multi-species
  # make_msm_test_data() produces a parametric-growth fixture (CAAL
  # populated, growth_model = 1), so fit_mod() must be told to use
  # vonBertalanffy growth -- otherwise it defaults to empirical and
  # data_check rejects the empirical-growth + CAAL combination.
  ss_run <- suppressMessages(
    Rceattle::fit_mod(data_list = simData,
                      estimateMode = 3, # Don't estimate
                      growthFun = build_growth(fun = "vonBertalanffy"),
                      msmMode = 0, # Single species mode
                      fit_control = fit_control(getsd = FALSE, 
                        verbose = 1))
  )
  inits <- ss_run$estimated_params
  inits$sel_inf[1,,1] <- c(3,6,2.5,4)
  inits$log_sel_slp[1,,1] <- log(c(2,2.5,2,2.5))
  inits$log_F[2,] <- log(Fmort)
  inits$log_F[4,] <- log(Fmort2)
  inits$rec_pars[,1] <- log(c(1e2, 1e3))
  inits$index_log_q[] <- log(1)
  inits$R_log_sd[] <- log(1)
  inits$rec_dev[,1:30] <- sim$model_quantities$rec_devs
  inits$init_dev[,1:14] <- sim$model_quantities$init_devs

  ss_run <- Rceattle::fit_mod(data_list = simData,
                              inits = inits, # Initial parameters from inits
                              file = NULL, # Don't save
                              estimateMode = 3, # Don't estimate
                              growthFun = build_growth(fun = "vonBertalanffy"),
                              recFun = build_srr(
                                proj_mean_rec = TRUE, # Project using mean rec over hindcast
                              ),
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0,
                              initMode = "NonEquilibrium",
                              fit_control = fit_control(getsd = FALSE, 
                                phase = FALSE,
                                verbose = 0))

  # Recruitment
  nyrs <- length(simData$styr:simData$endyr)
  testthat::expect_equal(as.numeric(sim$model_quantities$NAA[,1,]), as.numeric(ss_run$quantities$R[,1:nyrs]))

  # Mean R
  testthat::expect_equal(as.numeric(ss_run$quantities$avg_R), as.numeric(rowMeans(ss_run$quantities$R[,1:nyrs])), tolerance = 1e-6)
})


testthat::test_that("Invalid HCRs/BRPs", {
  data("BS2017SS") # ?BS2017SS for more information on the data
  BS2017SS$fleet_control$Proj_F_proportion <- rep(1, 7)

  testthat::expect_error(Rceattle::fit_mod(data_list = BS2017SS,
                                           estimateMode = 0, # Fit the hindcast and the HCR projection
                                           HCR = build_hcr(HCR = 10,
                                                           DynamicHCR = FALSE, # Static reference points
                                                           Ftarget = 0.4), # F that achieves 40% SB0
                                           msmMode = 0, # Single species mode
                                           fit_control = fit_control(getsd = FALSE, 
                                             verbose = 0))
  )


  testthat::expect_error(Rceattle::fit_mod(data_list = BS2017SS,
                                           estimateMode = 0, # Fit the hindcast and the HCR projection
                                           HCR = build_hcr(HCR = "BadHCR",
                                                           DynamicHCR = FALSE, # Static reference points
                                                           Ftarget = 0.4), # F that achieves 40% SB0
                                           msmMode = 0, # Single species mode
                                           fit_control = fit_control(getsd = FALSE, 
                                             verbose = 0))
  )
})


testthat::test_that("No F prop", {
  data("BS2017SS") # ?BS2017SS for more information on the data
  BS2017SS$fleet_control$Proj_F_proportion <- rep(0, 7)

  testthat::expect_error(Rceattle::fit_mod(data_list = BS2017SS,
                                           estimateMode = 0, # Fit the hindcast and the HCR projection
                                           HCR = build_hcr(HCR = 5, # Tier3 HCR
                                                           Ftarget = 0.4, # F40%
                                                           Flimit = 0.35, # F35%
                                                           Plimit = 0.2, # No fishing when SB<SB20
                                                           Alpha = 0.05),
                                           msmMode = 0, # Single species mode
                                           fit_control = fit_control(getsd = FALSE, 
                                             verbose = 0))
  )
}
)
