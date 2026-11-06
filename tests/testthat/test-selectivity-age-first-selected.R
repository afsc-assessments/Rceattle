# Integration test (fits a CEATTLE model): skipped on CRAN to keep R CMD check
# fast; the coverage workflow runs it in full (NOT_CRAN=true). See README.md.
testthat::skip_on_cran()

testthat::test_that("Test age-based non-parametric selectivity bin-first selected", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  data("GOA2018SS") # Single-species data. ?BS2017SS for more information on the data

  # Adjust data
  n_sel_bins <- 8
  rows_use <- which(GOA2018SS$fleet_control$Selectivity != 0 & GOA2018SS$fleet_control$Fleet_type != 0)
  GOA2018SS$fleet_control$Fleet_type[-rows_use] <- 0
  GOA2018SS$fleet_control$Selectivity <- "Logistic" # Age-based logistic
  GOA2018SS$fleet_control$N_sel_bins <- n_sel_bins
  GOA2018SS$fleet_control$Sel_curve_pen1 <- 5
  GOA2018SS$fleet_control$Sel_curve_pen2 <- 10
  GOA2018SS$fleet_control$Sel_norm_bin <- NA
  GOA2018SS$fleet_control$Age_first_selected <- 3
  GOA2018SS$fleet_control$Time_varying_sel <- 0
  GOA2018SS$fleet_control$Time_varying_sel_sd <- 1


  # Set params
  mod0 <- suppressMessages( fit_mod(data_list = GOA2018SS, inits = NULL, estimateMode = 3, random_rec = FALSE, msmMode = 0, fit_control = fit_control(verbose = 0)) )
  inits <- mod0$estimated_params
  log_selcoffs <- rnorm(n_sel_bins)
  log_selcoffs2 <- rnorm(n_sel_bins)
  inits$sel_coff[,1,1:8] <- rep(log_selcoffs, each = dim(inits$sel_coff)[1])
  inits$sel_coff[9:11,2,1:8] <- rep(log_selcoffs2, each = 3) # Males

  # Run
  ss_run <- Rceattle::fit_mod(data_list = GOA2018SS,
                              inits = inits, # Initial parameters = 0
                              file = NULL, # Don't save
                              estimateMode = 3, # Don't estimate
                              random_rec = FALSE, # No random recruitment
                              msmMode = 0, # Single species mode
                              fit_control = fit_control(
                                verbose = 0))

  # Check selectivity is 0 for bins 1-2
  output <- c(ss_run$quantities$sel_at_age[,1,1:2,])
  testthat::expect_equal(output, rep(0, length(output)))

  # Map for sel coff is NA?
  outputmap <- c(ss_run$map$mapList$sel_coff[,,1:2])
  testthat::expect_equal(outputmap, as.numeric(rep(NA, length(outputmap))))
})


testthat::test_that("Time-varying double logistic selectivity bin first selected", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Data
  data("GOA2018SS")
  rows_use <- which(GOA2018SS$fleet_control$Selectivity != 0 & GOA2018SS$fleet_control$Fleet_type != 0)
  GOA2018SS$fleet_control$Fleet_type[-rows_use] <-0
  GOA2018SS$fleet_control$Selectivity <- "DoubleLogistic" # Age-based double logistic
  GOA2018SS$fleet_control$Selectivity_index <- 1:nrow(GOA2018SS$fleet_control)
  GOA2018SS$fleet_control$Time_varying_sel <- 1
  GOA2018SS$fleet_control$Time_varying_sel_sd <- 1 # Note that inf does 4*sd prior, should adapt to scale-invariant sd
  #FIXME: move to scale invariant setup
  GOA2018SS$fleet_control$Age_first_selected <- 4
  GOA2018SS$fleet_control$Sel_norm_bin <- NA # Do not normalize
  GOA2018SS$catch_data$Catch <- 1e6 # If catch is zero, sel devs are turned off

  # Specify double logistic selectivity
  nyrs <- length(GOA2018SS$styr:GOA2018SS$endyr)
  inf = 7
  inf2 = 15
  inf_dev <- rnorm(nyrs)
  inf_dev_desc <- rnorm(nyrs)
  log_slp_dev <- rnorm(nyrs)
  log_slp_dev_desc <- rnorm(nyrs)

  alpha = 0.5
  ages <- 1:21
  sel <- apply(cbind(log_slp_dev, inf_dev, log_slp_dev_desc, inf_dev_desc), 1,
               function(x) 1/(1+exp(-(alpha)*exp(x[1]) * (ages - inf - x[2]))) * (1/(1+exp((alpha)*exp(x[3]) * (ages - inf2 - x[4])))))
  sel2 <- apply(cbind(log_slp_dev, inf_dev, log_slp_dev_desc, inf_dev_desc), 1,
                function(x) 1/(1+exp(-(alpha+1)*exp(x[1]) * (ages - inf - x[2]))) * (1/(1+exp((alpha+1)*exp(x[3]) * (ages - inf2 - x[4])))))


  # Set params to double logistic
  mod0 <- suppressMessages( fit_mod(data_list = GOA2018SS, inits = NULL, estimateMode = 3, random_rec = FALSE, random_sel = TRUE, msmMode = 0, fit_control = fit_control(verbose = 0)) )
  inits <- mod0$estimated_params
  inits$log_sel_slp[,,1] <- log(alpha) # Females
  inits$log_sel_slp[,,2] <- log(alpha+1) # Males
  inits$sel_inf[1,,] <- inf
  inits$sel_inf[2,,] <- inf2
  for(i in 1:dim(inits$log_sel_slp_dev[1,,,])[1]){
    for(j in 1:dim(inits$log_sel_slp_dev[1,,,])[2]){
      inits$log_sel_slp_dev[1,i,j,] <- log_slp_dev
      inits$sel_inf_dev[1,i,j,] <- inf_dev
      inits$log_sel_slp_dev[2,i,j,] <- log_slp_dev_desc
      inits$sel_inf_dev[2,i,j,] <- inf_dev_desc
    }
  }

  # Run
  ss_run <- suppressMessages(
    Rceattle::fit_mod(data_list = GOA2018SS,
                      inits = inits, # Initial parameters = 0
                      file = NULL, # Don't save
                      estimateMode = 3, # Don't estimate
                      random_rec = FALSE, # No random recruitment
                      random_sel = TRUE, # Turn on laplace for sel devs
                      msmMode = 0, # Single species mode
                      fit_control = fit_control(
                        verbose = 0))
  )

  # Check selectivity is 0 for bins 1-2
  output <- c(ss_run$quantities$sel_at_age[,1,1:3,])
  testthat::expect_equal(output, rep(0, length(output)))
})


# Age_first_selected is an absolute age, and on an age-based fleet
# Bin_first_selected is its deprecated spelling, read as age minage + bin - 1.
# minage = 3 here, because at minage = 1 the two readings coincide (rule 11).
.afs_bins <- function(d) {
  suppressMessages(rearrange_data(switch_check(d)))$bin_first_selected
}

testthat::test_that("an age-based fleet reads Bin_first_selected as the age it names", {
  d <- make_test_data(nyrs = 6, nprojyrs = 2, nages = 6, minage = 3)
  d$fleet_control$Age_first_selected <- c(5, NA)
  testthat::expect_identical(.afs_bins(d), c(2L, 0L))       # age 5 is bin 3

  old <- d
  old$fleet_control$Age_first_selected <- NULL
  old$fleet_control$Bin_first_selected <- c(3, NA)
  testthat::expect_identical(.afs_bins(old), c(2L, 0L))

  both <- old
  both$fleet_control$Age_first_selected <- c(5, NA)          # the same age, twice
  testthat::expect_identical(.afs_bins(both), c(2L, 0L))
  both$fleet_control$Bin_first_selected <- c(5, NA)          # age 7, not age 5
  testthat::expect_error(.afs_bins(both), "deprecated spelling of 'Age_first_selected'",
                         fixed = TRUE)

  bad <- d
  bad$fleet_control$Age_first_selected <- c(2, NA)           # below minage
  testthat::expect_error(suppressMessages(data_check(switch_check(bad))),
                         "Age_first_selected (2) must be an age in 3:8", fixed = TRUE)
  bad$fleet_control$Age_first_selected <- c(4.5, NA)
  testthat::expect_error(suppressMessages(data_check(switch_check(bad))),
                         "Age_first_selected (4.5) must be a whole age", fixed = TRUE)
})


testthat::test_that("the old spelling is reported once per fit", {
  d <- make_test_data(nyrs = 6, nprojyrs = 2, nages = 5)
  d$fleet_control$Age_first_selected <- NULL
  d$fleet_control$Bin_first_selected <- 1
  msgs <- character(0)
  withCallingHandlers(
    suppressWarnings(fit_mod(data_list = d, estimateMode = 3, msmMode = 0,
                             fit_control = fit_control(verbose = 0))),
    message = function(cnd) {
      msgs <<- c(msgs, conditionMessage(cnd))
      invokeRestart("muffleMessage")
    })
  testthat::expect_equal(sum(grepl("'Bin_first_selected' on an age-based fleet is deprecated",
                                   msgs, fixed = TRUE)), 1)
})


testthat::test_that("Age_first_selected round-trips through write_data() and read_data()", {
  d <- GOA2018SS
  f <- tempfile(fileext = ".xlsx")
  on.exit(unlink(f))
  suppressMessages(write_data(d, f))
  rt <- suppressMessages(read_data(f))
  testthat::expect_identical(as.numeric(rt$fleet_control$Age_first_selected),
                             as.numeric(d$fleet_control$Age_first_selected))
  testthat::expect_true(any(!is.na(d$fleet_control$Age_first_selected)))
})


# A length-based fleet selects no fish younger than Age_first_selected at any
# length (a Stock Synthesis age-range mask, age pattern 11): those ages leave the
# fleet's own age-length key, and nothing else moves.
testthat::test_that("on a length-based fleet, younger ages leave the fleet's age-length key", {
  set.seed(11)
  d   <- make_msm_test_data()$data_list
  flt <- which(d$fleet_control$Species == 1 & d$fleet_control$Fleet_type == 1)[1]
  d$fleet_control$Selectivity_dimension[flt] <- "Length"
  d$fleet_control$Selectivity[flt]           <- "Logistic"
  d$fleet_control$Time_varying_sel[flt]      <- "Off"
  sp  <- d$fleet_control$Species[flt]
  na  <- d$nages[sp]
  fit <- function(d) suppressWarnings(suppressMessages(fit_mod(
    data_list = d, estimateMode = 3, msmMode = 0, random_rec = FALSE,
    growthFun = build_growth(fun = "vonBertalanffy"),
    fit_control = fit_control(phase = FALSE, verbose = 0, getsd = FALSE))))

  base <- fit(d)
  # A floor at minage is the model without one.
  d0 <- d
  d0$fleet_control$Age_first_selected <- d$minage[d$fleet_control$Species]
  testthat::expect_identical(fit(d0)$quantities$jnll_comp, base$quantities$jnll_comp)

  d3 <- d
  d3$fleet_control$Age_first_selected[flt] <- d$minage[sp] + 2
  m3 <- fit(d3)
  testthat::expect_identical(m3$data_list$fleet_control$Bin_first_selected,
                             d$fleet_control$Bin_first_selected)

  wt <- 2 * d$nspp + flt   # growth_matrix: two slots per species, then one per fleet
  gm <- m3$quantities$growth_matrix
  g0 <- base$quantities$growth_matrix
  testthat::expect_true(all(gm[wt, , 1:2, , ] == 0))
  testthat::expect_gt(sum(g0[wt, , 1:2, , ]), 0)           # they were in the key
  testthat::expect_equal(gm[wt, , 3:na, , ], g0[wt, , 3:na, , ])
  testthat::expect_equal(gm[-wt, , , , ], g0[-wt, , , , ])

  s3 <- m3$quantities$sel_at_age
  s0 <- base$quantities$sel_at_age
  testthat::expect_true(all(s3[flt, , 1:2, ] == 0))
  testthat::expect_gt(min(s0[flt, 1, 1:2, 1]), 0)
  testthat::expect_equal(s3[flt, , 3:na, ], s0[flt, , 3:na, ])
  testthat::expect_equal(s3[-flt, , , ], s0[-flt, , , ])
})


# Before 4.0.1 'Bin_first_selected' was called 'Age_first_selected', and 167
# sibling-repository workbooks still carry only that header. Read as the new
# column, a script's later `fc$Bin_first_selected[2] <- 2` made a fresh column
# that R recycled down the rows: SSB 435,170 against 439,165 on GOA2018SS.
testthat::test_that("a workbook from before the rename reads its column as Bin_first_selected", {
  testthat::skip_if_not_installed("writexl")
  d <- make_test_data(nyrs = 6, nprojyrs = 2, nages = 5, minage = 2)
  d$fleet_control$Age_first_selected <- NULL
  d$fleet_control$Bin_first_selected <- c(2, NA)
  f <- withr::local_tempfile(fileext = ".xlsx")
  suppressMessages(write_data(d, f))
  sheets <- readxl::excel_sheets(f)
  wb <- stats::setNames(
    lapply(sheets, function(s) as.data.frame(readxl::read_xlsx(f, sheet = s))), sheets)
  fc <- wb[["fleet_control"]]
  fc[["Age_first_selected"]] <- NULL
  names(fc)[names(fc) == "Bin_first_selected"] <- "Age_first_selected"
  wb[["fleet_control"]] <- fc
  f_old <- withr::local_tempfile(fileext = ".xlsx")
  writexl::write_xlsx(wb, f_old)

  testthat::expect_message(old <- read_data(f_old), "pre-4.0.1 name", fixed = TRUE)
  testthat::expect_false("Age_first_selected" %in% names(old$fleet_control))
  testthat::expect_identical(.afs_bins(old), c(1L, 0L))      # bin 2, age 3
})
