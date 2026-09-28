# Fleets sharing a Selectivity_index share ONE selectivity deviation sd
# (sel_dev_log_sd). TMB collapses a shared parameter to the mean of its members'
# starting values (updateMap(): tapply(par, map, mean)), and this one is held on
# the log scale, so the group starts at the GEOMETRIC MEAN of their
# Time_varying_sel_sd -- no fleet keeps the value in its own row.
#
# That only happens when the sd is estimated, which is random_sel = TRUE with a
# random-effect Time_varying_sel. With random_sel = FALSE the map is NA and each
# fleet keeps its own value, so warning there would be a false alarm. Both
# directions are pinned here, as is the value the warning reports.

make_shared_sel_data <- function(sd1 = 0.3, sd2 = 0.7, tv_sel = "IID") {
  d <- make_test_data()
  fc <- d$fleet_control
  # Put both fleets in one selectivity group with the same curve, so the only
  # thing they disagree on is the deviation sd.
  fc$Selectivity_index   <- c(1, 1)
  fc$Selectivity         <- rep(fc$Selectivity[1], 2)
  fc$Time_varying_sel    <- rep(tv_sel, 2)
  fc$Time_varying_sel_sd <- c(sd1, sd2)
  d$fleet_control <- fc
  d
}

build_shared_sel_map <- function(d, random_sel) {
  suppressMessages(
    Rceattle::fit_mod(d, estimateMode = 3, msmMode = 0, random_sel = random_sel)
  )
}


test_that("a shared Selectivity_index warns when an estimated sd differs", {
  d <- make_shared_sel_data(sd1 = 0.3, sd2 = 0.7)

  expect_warning(
    build_shared_sel_map(d, random_sel = TRUE),
    "different Time_varying_sel_sd"
  )

  # The message names the group and the value the group actually starts at, so
  # the assessor can see that it is neither fleet's. Collect every warning:
  # build_map raises others (the Time_varying_sel soft-deprecation) ahead of it.
  hit <- grep("Time_varying_sel_sd", collect_warnings(
    build_shared_sel_map(d, random_sel = TRUE)), value = TRUE)
  expect_length(hit, 1)
  expect_match(hit, "Selectivity_index 1")
  expect_match(hit, "(0.3, 0.7)", fixed = TRUE)
  expect_match(hit, "0.4583", fixed = TRUE)   # sqrt(0.3 * 0.7)
})


test_that("the reported start is the one TMB uses", {
  # The warning would be worse than silence if it named a value the model does
  # not start from. Read it back off obj$par, which is what TMB optimizes.
  d <- make_shared_sel_data(sd1 = 0.3, sd2 = 0.7)
  m <- suppressWarnings(build_shared_sel_map(d, random_sel = TRUE))

  par <- m$obj$env$par
  start <- exp(unname(par[names(par) == "sel_dev_log_sd"]))
  expect_length(start, 1)                       # one shared parameter
  expect_equal(start, sqrt(0.3 * 0.7), tolerance = 1e-10)
})


test_that("one warning per group, however many fleets share the index", {
  # Warned per group rather than per additional member: fit_mod() de-duplicates
  # build_map() warnings by message text, so a message that grew with each
  # member would slip past it once per retrospective peel and MSE simulation.
  #
  # C carries no data of its own, so clean_data() turns it Off and its sd is
  # mapped out -- it keeps 0.9 and contributes nothing, which is why the message
  # names only the two fleets that share the estimated parameter.
  d  <- make_shared_sel_data(sd1 = 0.3, sd2 = 0.7)
  fc <- d$fleet_control[c(1, 2, 2), ]
  fc$Fleet_code          <- seq_len(nrow(fc))
  fc$Fleet_name          <- c("A", "B", "C")
  fc$Fleet_type          <- c("Fishery", "Survey", "Survey")
  fc$Time_varying_sel_sd <- c(0.3, 0.7, 0.9)
  d$fleet_control <- fc

  hit <- grep("Time_varying_sel_sd", collect_warnings(
    build_shared_sel_map(d, random_sel = TRUE)), value = TRUE)
  expect_length(hit, 1)
  expect_match(hit, "(A, B)", fixed = TRUE)
  expect_no_match(hit, "0.9", fixed = TRUE)
})


test_that("an Off fleet in the group is not counted", {
  # A Fleet_type "Off" fleet has sel_dev_log_sd mapped out (in
  # build_map_f_and_data_weights(), AFTER the sharing pass), so it keeps its own
  # value, contributes nothing to the mean, and is nothing to warn about.
  d <- make_shared_sel_data(sd1 = 0.3, sd2 = 0.7)
  d$fleet_control$Fleet_type[2] <- "Off"

  hit <- grep("Time_varying_sel_sd", collect_warnings(
    build_shared_sel_map(d, random_sel = TRUE)), value = TRUE)
  expect_length(hit, 0)
})


test_that("agreeing fleets are silent", {
  d <- make_shared_sel_data(sd1 = 0.5, sd2 = 0.5)
  expect_no_warning(
    withCallingHandlers(
      build_shared_sel_map(d, random_sel = TRUE),
      warning = function(w) {
        if (!grepl("Time_varying_sel_sd", conditionMessage(w))) {
          invokeRestart("muffleWarning")
        }
      }
    )
  )
})


test_that("random_q estimates the catchability deviation sd, and random_sel the selectivity one", {
  # The two are symmetric: integrating the deviates out is what makes their sd a
  # marginal variance rather than one half of a degenerate joint mode. Without
  # the flag the sd stays at its fleet_control value.
  d <- make_test_data()
  d$fleet_control$Catchability[1]      <- "Estimated"
  d$fleet_control$Time_varying_q[1]    <- "IID"
  d$fleet_control$Time_varying_q_sd[1] <- 0.2

  est <- function(rq) {
    m <- suppressWarnings(suppressMessages(
      Rceattle::fit_mod(d, estimateMode = 3, msmMode = 0, random_q = rq)))
    m$map$mapFactor
  }

  expect_true(all(is.na(est(FALSE)$index_q_dev_log_sd)))
  expect_true(any(!is.na(est(TRUE)$index_q_dev_log_sd)))

  # index_q_log_sd is a prior sd the assessor sets; random_q must not free it.
  expect_true(all(is.na(est(TRUE)$index_q_log_sd)))
})


test_that("a fixed deviation sd is honoured per fleet, so it does not warn", {
  # random_sel = FALSE leaves sel_dev_log_sd mapped out. Nothing is shared and
  # nothing is discarded, so differing values are legitimate.
  d <- make_shared_sel_data(sd1 = 0.3, sd2 = 0.7)
  m <- build_shared_sel_map(d, random_sel = FALSE)
  expect_true(all(is.na(m$map$mapFactor$sel_dev_log_sd)))

  expect_no_warning(
    withCallingHandlers(
      build_shared_sel_map(d, random_sel = FALSE),
      warning = function(w) {
        if (!grepl("Time_varying_sel_sd", conditionMessage(w))) {
          invokeRestart("muffleWarning")
        }
      }
    )
  )
})


test_that("a shared Catchability_index reports the same way", {
  # The q side goes through the same helper, so the number it reports must also
  # be the one TMB starts from: index_q_dev_log_sd is on the log scale too.
  d <- make_test_data()
  fc <- d$fleet_control
  fc$Catchability_index <- c(1, 1)
  fc$Catchability       <- rep("Estimated", 2)
  fc$Time_varying_q     <- rep("IID", 2)
  fc$Time_varying_q_sd  <- c(0.2, 0.8)
  d$fleet_control <- fc

  m <- suppressWarnings(suppressMessages(
    Rceattle::fit_mod(d, estimateMode = 3, msmMode = 0, random_q = TRUE)))
  par <- m$obj$env$par
  start <- exp(unname(par[names(par) == "index_q_dev_log_sd"]))
  expect_length(start, 1)
  expect_equal(start, sqrt(0.2 * 0.8), tolerance = 1e-10)

  hit <- grep("Time_varying_q_sd", collect_warnings(suppressMessages(
    Rceattle::fit_mod(d, estimateMode = 3, msmMode = 0, random_q = TRUE))),
    value = TRUE)
  expect_length(hit, 1)
  expect_match(hit, "Catchability_index 1")
  expect_match(hit, "0.4", fixed = TRUE)
})


test_that("build_map() takes random_q directly, and still reads the data_list", {
  # random_rec and random_sel were arguments while random_q was only ever read
  # off data_list, so a direct build_map() call could not ask for it. It is an
  # argument now, defaulting to what fit_mod() stores, so both routes work.
  d <- make_test_data()
  d$fleet_control$Catchability[1]      <- "Estimated"
  d$fleet_control$Time_varying_q[1]    <- "IID"
  d$fleet_control$Time_varying_q_sd[1] <- 0.2

  # A build_map()-ready data_list: fit_mod() prepares one on the way to the map.
  fit <- suppressWarnings(suppressMessages(
    Rceattle::fit_mod(d, estimateMode = 3, msmMode = 0, random_q = FALSE)))
  dl <- fit$data_list
  dl$random_q <- NULL                          # as if fit_mod() never ran
  p  <- suppressWarnings(Rceattle::build_params(dl))

  est <- function(...) {
    m <- suppressWarnings(suppressMessages(Rceattle::build_map(dl, p, ...)))
    any(!is.na(m$mapList$index_q_dev_log_sd))
  }

  expect_false(est())                          # no flag anywhere
  expect_true(est(random_q = TRUE))            # passed directly
  expect_false(est(random_q = FALSE))          # explicitly off

  dl$random_q <- 1                             # what fit_mod() stores
  expect_true(est())                           # picked up by default
  expect_false(est(random_q = FALSE))          # an explicit argument still wins
})
