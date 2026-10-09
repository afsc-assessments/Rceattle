# Integration test (fits a CEATTLE model): skipped on CRAN to keep R CMD check
# fast; the coverage workflow runs it in full (NOT_CRAN=true). See README.md.
testthat::skip_on_cran()

testthat::test_that("Sex-specific logistic selectivity divided by max sel (across sexes)", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")


  data("GOA2018SS") # Single-species data. ?BS2017SS for more information on the data
  rows_use <- which(GOA2018SS$fleet_control$Selectivity != 0 & GOA2018SS$fleet_control$Fleet_type != 0)
  GOA2018SS$fleet_control$Selectivity <-0
  GOA2018SS$fleet_control$Selectivity[rows_use] <- "Logistic" # Age-based logistic
  # Fleet 8 ships RandomWalkAscending, which freed no deviates on a
  # Logistic curve, so this block has always fitted a static one. These
  # tests are about normalization, not time variation.
  GOA2018SS$fleet_control$Time_varying_sel[rows_use] <- "Off"
  GOA2018SS$fleet_control$Selectivity_index <- 1:nrow(GOA2018SS$fleet_control)
  GOA2018SS$fleet_control$Bin_first_selected <- 1
  GOA2018SS$fleet_control$Sel_norm_bin <- -999


  # Specify logistic selectivity
  inf = 10; alpha = 0.5
  ages <- 1:21
  sel <- 1/(1+exp(-alpha*(ages-inf)))
  sel2 <- 1/(1+exp(-alpha*(ages-inf-1)))
  # curve(1/(1+exp(-alpha*(x-inf))), from = 0, to = 21)

  # Set params to logistic
  mod0 <- suppressMessages( fit_mod(data_list = GOA2018SS, inits = NULL, estimateMode = 3, random_rec = FALSE, msmMode = 0, fit_control = fit_control(verbose = 0)) )
  inits <- mod0$estimated_params
  inits$log_sel_slp[1,,] <- log(alpha)
  inits$sel_inf[1,,] <- inf     # Females
  inits$sel_inf[1,9:11,2] <- inf + 1 # Males

  # Run
  ss_run <- suppressMessages(
    Rceattle::fit_mod(data_list = GOA2018SS,
                      inits = inits, # Initial parameters = 0
                      file = NULL, # Don't save
                      estimateMode = 3, # Don't estimate
                      random_rec = FALSE, # No random recruitment
                      msmMode = 0, # Single species mode
                      fit_control = fit_control(
                        verbose = 1))
  )

  # Check selectivity
  # - Pollock
  apply(ss_run$quantities$sel_at_age[1,1,1:10,], 2, function(x) testthat::expect_equal(as.numeric(x), sel[1:10]/max(sel[1:10]), tolerance = 0.0001))


  # - ATF
  apply(ss_run$quantities$sel_at_age[9,1,1:21,], 2, function(x) testthat::expect_equal(as.numeric(x), sel[1:21]/max(c(sel, sel2)), tolerance = 0.0001)) # Females
  apply(ss_run$quantities$sel_at_age[9,2,1:21,], 2, function(x) testthat::expect_equal(as.numeric(x), sel2[1:21]/max(c(sel, sel2)), tolerance = 0.0001)) # Males
})


testthat::test_that("Sex-specific logistic selectivity not normalized", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("GOA2018SS") # Single-species data. ?BS2017SS for more information on the data
  rows_use <- which(GOA2018SS$fleet_control$Selectivity != 0 & GOA2018SS$fleet_control$Fleet_type != 0)
  GOA2018SS$fleet_control$Selectivity <-0
  GOA2018SS$fleet_control$Selectivity[rows_use] <- "Logistic" # Age-based logistic
  # Fleet 8 ships RandomWalkAscending, which freed no deviates on a
  # Logistic curve, so this block has always fitted a static one. These
  # tests are about normalization, not time variation.
  GOA2018SS$fleet_control$Time_varying_sel[rows_use] <- "Off"
  GOA2018SS$fleet_control$Selectivity_index <- 1:nrow(GOA2018SS$fleet_control)
  GOA2018SS$fleet_control$Bin_first_selected <- 1
  GOA2018SS$fleet_control$Sel_norm_bin <- NA

  # Specify logistic selectivity
  mod0 <- suppressMessages( fit_mod(data_list = GOA2018SS, inits = NULL, estimateMode = 3, random_rec = FALSE, msmMode = 0, fit_control = fit_control(verbose = 0)) )
  inits <- mod0$estimated_params
  inf = 13; alpha = 0.2
  ages <- 1:21
  sel <- 1/(1+exp(-alpha*(ages-inf)))
  sel2 <- 1/(1+exp(-alpha*(ages-inf-1)))
  # curve(1/(1+exp(-alpha*(x-inf))), from = 0, to = 21)

  # Set params to logistic
  inits$log_sel_slp[1,,] <- log(alpha)
  inits$sel_inf[1,,] <- inf     # Females
  inits$sel_inf[1,9:11,2] <- inf + 1 # Males

  # Run
  ss_run <- suppressMessages(
    Rceattle::fit_mod(data_list = GOA2018SS,
                      inits = inits, # Initial parameters = 0
                      file = NULL, # Don't save
                      estimateMode = 3, # Don't estimate
                      random_rec = FALSE, # No random recruitment
                      msmMode = 0, # Single species mode
                      fit_control = fit_control(
                        verbose = 1))
  )

  # - Devs
  testthat::expect_equal(c(ss_run$map$mapList$log_sel_slp_dev), as.numeric(rep(NA, 2688)))
  testthat::expect_equal(c(ss_run$map$mapList$sel_inf_dev), as.numeric(rep(NA, 2688)))
  testthat::expect_equal(as.numeric(ss_run$map$mapList$sel_dev_log_sd), as.numeric(rep(NA, 16))) # Dev sigma


  # Check selectivity
  # - Pollock
  apply(ss_run$quantities$sel_at_age[1,1,1:10,], 2, function(x) testthat::expect_equal(as.numeric(x), sel[1:10], tolerance = 0.0001))


  # - ATF
  apply(ss_run$quantities$sel_at_age[9,1,1:21,], 2, function(x) testthat::expect_equal(as.numeric(x), sel[1:21], tolerance = 0.0001)) # Females
  apply(ss_run$quantities$sel_at_age[9,2,1:21,], 2, function(x) testthat::expect_equal(as.numeric(x), sel2[1:21], tolerance = 0.0001)) # Males
})


testthat::test_that("Sex-invariant logistic selectivity divided by sel-at-age", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("GOA2018SS") # Single-species data. ?BS2017SS for more information on the data
  rows_use <- which(GOA2018SS$fleet_control$Selectivity != 0 & GOA2018SS$fleet_control$Fleet_type != 0)
  GOA2018SS$fleet_control$Selectivity <-0
  GOA2018SS$fleet_control$Selectivity[rows_use] <- "Logistic" # Age-based logistic
  # Fleet 8 ships RandomWalkAscending, which freed no deviates on a
  # Logistic curve, so this block has always fitted a static one. These
  # tests are about normalization, not time variation.
  GOA2018SS$fleet_control$Time_varying_sel[rows_use] <- "Off"
  GOA2018SS$fleet_control$Bin_first_selected <- 1
  GOA2018SS$fleet_control$Sel_norm_bin <- 7

  # Specify logistic selectivity
  inf = 10; alpha = 0.5
  ages <- 1:21
  sel <- 1/(1+exp(-alpha*(ages-inf)))
  # curve(1/(1+exp(-alpha*(x-inf))), from = 0, to = 21)

  # Set params to logistic
  mod0 <- suppressMessages( fit_mod(data_list = GOA2018SS, inits = NULL, estimateMode = 3, random_rec = FALSE, msmMode = 0, fit_control = fit_control(verbose = 0)) )
  inits <- mod0$estimated_params
  inits$log_sel_slp[1,,] <- log(alpha)
  inits$sel_inf[1,,] <- inf

  # Run
  ss_run <- suppressMessages(
    Rceattle::fit_mod(data_list = GOA2018SS,
                      inits = inits, # Initial parameters = 0
                      file = NULL, # Don't save
                      estimateMode = 3, # Don't estimate
                      random_rec = FALSE, # No random recruitment
                      msmMode = 0, # Single species mode
                      fit_control = fit_control(
                        verbose = 1))
  )

  # Check selectivity
  # - Pollock
  apply(ss_run$quantities$sel_at_age[1,1,1:10,], 2, function(x) testthat::expect_equal(as.numeric(x), sel[1:10]/sel[7], tolerance = 0.0001))


  # - ATF
  apply(ss_run$quantities$sel_at_age[9,1,1:21,], 2, function(x) testthat::expect_equal(as.numeric(x), sel[1:21]/sel[7], tolerance = 0.0001))
  apply(ss_run$quantities$sel_at_age[9,2,1:21,], 2, function(x) testthat::expect_equal(as.numeric(x), sel[1:21]/sel[7], tolerance = 0.0001))
})



testthat::test_that("Sex-invariant logistic selectivity divided by sel-at-age-RANGE", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("GOA2018SS") # Single-species data. ?BS2017SS for more information on the data
  rows_use <- which(GOA2018SS$fleet_control$Selectivity != 0 & GOA2018SS$fleet_control$Fleet_type != 0)
  GOA2018SS$fleet_control$Selectivity <-0
  GOA2018SS$fleet_control$Selectivity[rows_use] <- "Logistic" # Age-based logistic
  # Fleet 8 ships RandomWalkAscending, which freed no deviates on a
  # Logistic curve, so this block has always fitted a static one. These
  # tests are about normalization, not time variation.
  GOA2018SS$fleet_control$Time_varying_sel[rows_use] <- "Off"
  GOA2018SS$fleet_control$Bin_first_selected <- 1
  GOA2018SS$fleet_control$Sel_norm_bin <- 7
  GOA2018SS$fleet_control$Sel_norm_bin_upper <- 9

  # Specify logistic selectivity
  inf = 10; alpha = 0.5
  ages <- 1:21
  sel <- 1/(1+exp(-alpha*(ages-inf)))
  # curve(1/(1+exp(-alpha*(x-inf))), from = 0, to = 21)

  # Set params to logistic
  mod0 <- suppressMessages( fit_mod(data_list = GOA2018SS, inits = NULL, estimateMode = 3, random_rec = FALSE, msmMode = 0, fit_control = fit_control(verbose = 0)) )
  inits <- mod0$estimated_params
  inits$log_sel_slp[1,,] <- log(alpha)
  inits$sel_inf[1,,] <- inf

  # Run
  ss_run <- suppressMessages(
    Rceattle::fit_mod(data_list = GOA2018SS,
                      inits = inits, # Initial parameters = 0
                      file = NULL, # Don't save
                      estimateMode = 3, # Don't estimate
                      random_rec = FALSE, # No random recruitment
                      msmMode = 0, # Single species mode
                      fit_control = fit_control(
                        verbose = 1))
  )

  # Check selectivity
  # - Pollock
  apply(ss_run$quantities$sel_at_age[1,1,1:10,], 2, function(x) testthat::expect_equal(as.numeric(x), sel[1:10]/mean(sel[7:9]), tolerance = 0.0001))


  # - ATF
  apply(ss_run$quantities$sel_at_age[9,1,1:21,], 2, function(x) testthat::expect_equal(as.numeric(x), sel[1:21]/mean(sel[7:9]), tolerance = 0.0001))
  apply(ss_run$quantities$sel_at_age[9,2,1:21,], 2, function(x) testthat::expect_equal(as.numeric(x), sel[1:21]/mean(sel[7:9]), tolerance = 0.0001))
})


testthat::test_that("Sex-specific logistic selectivity divided by sel-at-age-RANGE", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("GOA2018SS") # Single-species data. ?BS2017SS for more information on the data
  rows_use <- which(GOA2018SS$fleet_control$Selectivity != 0 & GOA2018SS$fleet_control$Fleet_type != 0)
  GOA2018SS$fleet_control$Selectivity <-0
  GOA2018SS$fleet_control$Selectivity[rows_use] <- "Logistic" # Age-based logistic
  # Fleet 8 ships RandomWalkAscending, which freed no deviates on a
  # Logistic curve, so this block has always fitted a static one. These
  # tests are about normalization, not time variation.
  GOA2018SS$fleet_control$Time_varying_sel[rows_use] <- "Off"
  GOA2018SS$fleet_control$Bin_first_selected <- 1
  # Each sex divided by its OWN mean over the range. Before v5.8.0 a named bin
  # always meant this; it is now Sel_norm_scope = "WithinSex", and the pooled
  # default is the companion test below. Species 9 (arrowtooth) is the two-sex
  # one, so it is the only place the two differ.
  GOA2018SS$fleet_control$Sel_norm_scope <- "WithinSex"
  GOA2018SS$fleet_control$Sel_norm_bin <- 7
  GOA2018SS$fleet_control$Sel_norm_bin_upper <- 9

  # Specify logistic selectivity
  inf = 10; alpha = 0.5
  ages <- 1:21
  sel <- 1/(1+exp(-alpha*(ages-inf)))
  sel2 <- 1/(1+exp(-(alpha+1)*(ages-inf)))
  # curve(1/(1+exp(-alpha*(x-inf))), from = 0, to = 21)

  # Set params to logistic
  mod0 <- suppressMessages( fit_mod(data_list = GOA2018SS, inits = NULL, estimateMode = 3, random_rec = FALSE, msmMode = 0, fit_control = fit_control(verbose = 0)) )
  inits <- mod0$estimated_params
  inits$log_sel_slp[1,,1] <- log(alpha)     # Females
  inits$log_sel_slp[1,,2] <- log(alpha + 1) # Males
  inits$sel_inf[1,,] <- inf

  # Run
  ss_run <- suppressMessages(
    Rceattle::fit_mod(data_list = GOA2018SS,
                      inits = inits, # Initial parameters = 0
                      file = NULL, # Don't save
                      estimateMode = 3, # Don't estimate
                      random_rec = FALSE, # No random recruitment
                      msmMode = 0, # Single species mode
                      fit_control = fit_control(
                        verbose = 1))
  )

  # Check selectivity
  # - Pollock
  apply(ss_run$quantities$sel_at_age[1,1,1:10,], 2, function(x) testthat::expect_equal(as.numeric(x), sel[1:10]/mean(sel[7:9]), tolerance = 0.0001))


  # - ATF
  apply(ss_run$quantities$sel_at_age[9,1,1:21,], 2, function(x) testthat::expect_equal(as.numeric(x), sel[1:21]/mean(sel[7:9]), tolerance = 0.0001))
  apply(ss_run$quantities$sel_at_age[9,2,1:21,], 2, function(x) testthat::expect_equal(as.numeric(x), sel2[1:21]/mean(sel2[7:9]), tolerance = 0.0001))
})


testthat::test_that("AcrossSexes pools the sel-at-age-RANGE reference over both sexes", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Same setup as the WithinSex test above, under the v5.8.0 default: one
  # reference pooled over the sexes, so the less-selected sex keeps its relative
  # level instead of being lifted to 1. Males are the steeper curve here, so
  # their mean over ages 7-9 is the smaller of the two and the pooled reference
  # is the females'.
  data("GOA2018SS")
  rows_use <- which(GOA2018SS$fleet_control$Selectivity != 0 & GOA2018SS$fleet_control$Fleet_type != 0)
  GOA2018SS$fleet_control$Selectivity <- 0
  GOA2018SS$fleet_control$Selectivity[rows_use] <- "Logistic"
  # Fleet 8 ships RandomWalkAscending, which freed no deviates on a
  # Logistic curve, so this block has always fitted a static one. These
  # tests are about normalization, not time variation.
  GOA2018SS$fleet_control$Time_varying_sel[rows_use] <- "Off"
  GOA2018SS$fleet_control$Bin_first_selected <- 1
  GOA2018SS$fleet_control$Sel_norm_scope <- "AcrossSexes"
  GOA2018SS$fleet_control$Sel_norm_bin <- 7
  GOA2018SS$fleet_control$Sel_norm_bin_upper <- 9

  inf <- 10; alpha <- 0.5
  ages <- 1:21
  sel  <- 1/(1+exp(-alpha*(ages-inf)))        # females
  sel2 <- 1/(1+exp(-(alpha+1)*(ages-inf)))    # males, steeper

  mod0 <- suppressMessages(fit_mod(data_list = GOA2018SS, inits = NULL,
    estimateMode = 3, random_rec = FALSE, msmMode = 0,
    fit_control = fit_control(verbose = 0)))
  inits <- mod0$estimated_params
  inits$log_sel_slp[1,,1] <- log(alpha)
  inits$log_sel_slp[1,,2] <- log(alpha + 1)
  inits$sel_inf[1,,] <- inf

  ss_run <- suppressMessages(Rceattle::fit_mod(data_list = GOA2018SS,
    inits = inits, file = NULL, estimateMode = 3, random_rec = FALSE,
    msmMode = 0, fit_control = fit_control(verbose = 0)))

  pooled <- max(mean(sel[7:9]), mean(sel2[7:9]))
  expect_gt(pooled, mean(sel2[7:9]))          # the pooling actually bites

  apply(ss_run$quantities$sel_at_age[9,1,1:21,], 2, function(x)
    testthat::expect_equal(as.numeric(x), sel[1:21]/pooled, tolerance = 0.0001))
  apply(ss_run$quantities$sel_at_age[9,2,1:21,], 2, function(x)
    testthat::expect_equal(as.numeric(x), sel2[1:21]/pooled, tolerance = 0.0001))
})



testthat::test_that("Sex-invariant time-varying logistic selectivity divided by max sel", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("GOA2018SS") # Single-species data. ?BS2017SS for more information on the data
  rows_use <- which(GOA2018SS$fleet_control$Selectivity != 0 & GOA2018SS$fleet_control$Fleet_type != 0)
  GOA2018SS$fleet_control$Selectivity <-0
  GOA2018SS$fleet_control$Selectivity[rows_use] <- "Logistic" # Age-based logistic
  GOA2018SS$fleet_control$Time_varying_sel <- 1
  GOA2018SS$fleet_control$Time_varying_sel_sd <- 1
  GOA2018SS$fleet_control$Bin_first_selected <- 1
  GOA2018SS$fleet_control$Sel_norm_bin <- -999

  # Specify logistic selectivity
  nyrs <- GOA2018SS$styr:GOA2018SS$endyr
  inf = 10
  inf_dev <- rnorm(nyrs)
  log_slp_dev <- rnorm(nyrs)

  alpha = 0.5
  ages <- 1:21
  sel <- apply(cbind(log_slp_dev, inf_dev), 1, function(x) 1/(1+exp(-alpha*exp(x[1]) * (ages - inf - x[2]))))

  # Set params to logistic
  mod0 <- suppressMessages( fit_mod(data_list = GOA2018SS, inits = NULL, estimateMode = 3, random_rec = FALSE, msmMode = 0, fit_control = fit_control(verbose = 0)) )
  inits <- mod0$estimated_params
  inits$log_sel_slp[1,,] <- log(alpha)
  inits$sel_inf[1,,] <- inf
  for(i in 1:dim(inits$log_sel_slp_dev[1,,,])[1]){
    for(j in 1:dim(inits$log_sel_slp_dev[1,,,])[2]){
      inits$log_sel_slp_dev[1,i,j,] <- log_slp_dev
      inits$sel_inf_dev[1,i,j,] <- inf_dev
    }
  }

  # Run
  ss_run <- suppressMessages(
    Rceattle::fit_mod(data_list = GOA2018SS,
                      inits = inits, # Initial parameters = 0
                      file = NULL, # Don't save
                      estimateMode = 3, # Don't estimate
                      random_rec = FALSE, # No random recruitment
                      msmMode = 0, # Single species mode
                      fit_control = fit_control(
                        verbose = 0))
  )

  # Check selectivity
  # - Pollock
  testthat::expect_equal(as.numeric(ss_run$quantities$sel_at_age[1,1,1:10,1:42]), as.numeric(apply(sel[1:10,], 2, function(x) x/max(x))), tolerance = 0.0001)


  # - ATF
  testthat::expect_equal(as.numeric(ss_run$quantities$sel_at_age[9,1,1:21,1:42]), as.numeric(apply(sel[1:21,], 2, function(x) x/max(x))), tolerance = 0.0001) # Females
  testthat::expect_equal(as.numeric(ss_run$quantities$sel_at_age[9,2,1:21,1:42]), as.numeric(apply(sel[1:21,], 2, function(x) x/max(x))), tolerance = 0.0001) # Males
})


testthat::test_that("Normalize by max for each fishery and year across bins, and sexes", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("GOA2018SS") # Single-species data. ?BS2017SS for more information on the data
  rows_use <- which(GOA2018SS$fleet_control$Selectivity != 0 & GOA2018SS$fleet_control$Fleet_type != 0)
  GOA2018SS$fleet_control$Selectivity <-0
  GOA2018SS$fleet_control$Selectivity[rows_use] <- "Logistic" # Age-based logistic
  GOA2018SS$fleet_control$Time_varying_sel <- 1
  GOA2018SS$fleet_control$Time_varying_sel_sd <- 1
  GOA2018SS$fleet_control$Bin_first_selected <- 1
  GOA2018SS$fleet_control$Sel_norm_bin <- -999

  # Specify logistic selectivity
  nyrs <- length(GOA2018SS$styr:GOA2018SS$endyr)
  inf = 10
  inf_dev <- rnorm(nyrs)
  log_slp_dev <- rnorm(nyrs)

  alpha = 0.5
  ages <- 1:21
  sel <- apply(cbind(log_slp_dev, inf_dev), 1, function(x) 1/(1+exp(-alpha*exp(x[1]) * (ages - inf - x[2])))) # Females
  sel2 <- apply(cbind(log_slp_dev, inf_dev), 1, function(x) 1/(1+exp(-(alpha+1)*exp(x[1]) * (ages - inf - x[2])))) # Males

  # Set params to logistic
  mod0 <- suppressMessages( fit_mod(data_list = GOA2018SS, inits = NULL, estimateMode = 3, random_rec = FALSE, msmMode = 0, fit_control = fit_control(verbose = 0)) )
  inits <- mod0$estimated_params
  inits$log_sel_slp[1,,1] <- log(alpha) # Females
  inits$log_sel_slp[1,,2] <- log(alpha+1) # Males
  inits$sel_inf[1,,] <- inf
  for(i in 1:dim(inits$log_sel_slp_dev[1,,,])[1]){
    for(j in 1:dim(inits$log_sel_slp_dev[1,,,])[2]){
      inits$log_sel_slp_dev[1,i,j,] <- log_slp_dev
      inits$sel_inf_dev[1,i,j,] <- inf_dev
    }
  }

  # Run
  ss_run <- suppressMessages(
    Rceattle::fit_mod(data_list = GOA2018SS,
                      inits = inits, # Initial parameters = 0
                      file = NULL, # Don't save
                      estimateMode = 3, # Don't estimate
                      random_rec = FALSE, # No random recruitment
                      msmMode = 0, # Single species mode
                      fit_control = fit_control(
                        verbose = 0))
  )


  # Normalize by max for each year across bins, and sexes
  sel2_norm <- sel2
  sel_norm <- sel
  for(i in 1:nyrs){
    for(a in ages){
      sel_norm[a,i] <- sel[a,i] / max(c(sel[,i], sel2[,i]))
      sel2_norm[a,i] <- sel2[a,i] / max(c(sel[,i], sel2[,i]))
    }
  }

  # Check selectivity
  # - Pollock
  testthat::expect_equal(as.numeric(ss_run$quantities$sel_at_age[1,1,1:10,1:42]), as.numeric(apply(sel[1:10,], 2, function(x) x/max(x))), tolerance = 0.0001)


  # - ATF
  testthat::expect_equal(as.numeric(ss_run$quantities$sel_at_age[9,1,1:21,1:42]), as.numeric(sel_norm[1:21,]), tolerance = 0.0001) # Females
  testthat::expect_equal(as.numeric(ss_run$quantities$sel_at_age[9,2,1:21,1:42]), as.numeric(sel2_norm[1:21,]), tolerance = 0.0001) # Males
})

