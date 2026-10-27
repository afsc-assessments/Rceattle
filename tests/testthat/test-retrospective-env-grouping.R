# A retrospective peel blanks post-peel env_data so it cannot read the future,
# but a random-effect grouping column is a level label, not an observation.
# Blanked, model.matrix() drops those rows from the random-effect design and the
# peel refit stopped with "number of rows of matrices must match".

testthat::test_that("a peel keeps a random-effect grouping column", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")

  d <- make_test_data(nyrs = 20, seed = 1)
  d$fleet_control$Catchability[1] <- "Estimated"
  yrs <- d$styr:d$endyr
  d$env_data <- data.frame(Year = yrs, regime = ifelse(yrs < d$styr + 12, 1, 2),
                           unused = seq_along(yrs))
  qf <- Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ (1 | regime), by = ~ fleet, fleet = 1L)))
  fit <- suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, inits = NULL, file = NULL, estimateMode = 0,
    random_rec = FALSE, msmMode = 0, qFun = qf,
    fit_control = Rceattle::fit_control(phase = FALSE, getsd = FALSE, verbose = 0))))

  ws <- character(0)
  r <- withCallingHandlers(suppressMessages(
    Rceattle::retrospective(fit, peels = 2, cores = 1, getsd = FALSE)),
    warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") })
  testthat::expect_length(r$Rceattle_list, 3L)   # both peels plus the parent
  # Kept post-peel values are named, as a fixed-effect covariate's are.
  testthat::expect_true(any(grepl("grouping column(s) regime", ws, fixed = TRUE)))

  # The deepest peel: its post-peel rows keep `regime` but still lose `unused`.
  env <- r$Rceattle_list[[1]]$data_list$env_data
  post <- env$Year > d$endyr - 2
  testthat::expect_true(any(post))
  testthat::expect_identical(env$regime, d$env_data$regime)
  testthat::expect_true(all(is.na(env$unused[post])))
})
